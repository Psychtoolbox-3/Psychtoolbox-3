/*
    PsychToolbox3/Source/OSX/PsychHID/PsychHIDStandardInterfaces.c

    PROJECTS: PsychHID

    PLATFORMS:  OSX

    AUTHORS:

    rwoods@ucla.edu                 rpw
    mario.kleiner.de@gmail.com      mk

    HISTORY:

    19.08.2007      rpw     Created the original implementation used before Psychtoolbox 3.0.12.
    2008 - 2014     mk      Various improvements and bug fixes to rpw's original implementation.
    04.10.2014      mk      Refactored and almost completely rewritten for PTB 3.0.12.

    TO DO:

*/

#include "PsychHIDStandardInterfaces.h"
#include <errno.h>

#define NUMDEVICEUSAGES 8

static  UInt32 modifierKeyState[PSYCH_HID_MAX_DEVICES];
static  double* psychHIDKbQueueFirstPress[PSYCH_HID_MAX_DEVICES];
static  double* psychHIDKbQueueFirstRelease[PSYCH_HID_MAX_DEVICES];
static  double* psychHIDKbQueueLastPress[PSYCH_HID_MAX_DEVICES];
static  double* psychHIDKbQueueLastRelease[PSYCH_HID_MAX_DEVICES];
static  int*    psychHIDKbQueueScanKeys[PSYCH_HID_MAX_DEVICES];
static  psych_bool psychHIDKbQueueActive[PSYCH_HID_MAX_DEVICES];
static  psych_bool queueIsAKeyboard[PSYCH_HID_MAX_DEVICES];
static  IOHIDQueueRef queue[PSYCH_HID_MAX_DEVICES];
static  CFRunLoopRef psychHIDKbQueueCFRunLoopRef[PSYCH_HID_MAX_DEVICES];
static  psych_thread KbQueueThread[PSYCH_HID_MAX_DEVICES];

static  psych_mutex KbQueueMutex;
static  psych_condition KbQueueCondition;
static  int ndevices = 0;
static  int deviceIndices[PSYCH_HID_MAX_DEVICES];
static  pRecDevice deviceRecords[PSYCH_HID_MAX_DEVICES];
static  int defaultKeyboardIndex = 0;
static  const UCKeyboardLayout *keyboardLayout = NULL;
static  UInt8 kbdType = 0;

// Return real deviceIndex (aka KbQueue index) and HID device record pRecDevice
// for a given PsychHID('Devices') style input deviceIndex. *deviceIndex is an
// in/out parameter:
int PsychHIDOSGetKbQueueDevice(int HIDdeviceIndex, pRecDevice *deviceRecord)
{
    psych_bool foundUserSpecifiedDevice = FALSE;
    int i;

    // A negative device number causes the default device to be used:
    if (HIDdeviceIndex >= 0) {
        // Make sure that the device number provided by the user is really a suitable input device:
        for (i = 0; i < ndevices; i++) {
            if ((foundUserSpecifiedDevice = (deviceIndices[i] == HIDdeviceIndex))) break;
        }

        if (!foundUserSpecifiedDevice)
            PsychErrorExitMsg(PsychError_user, "Specified device number is not a suitable keyboard type input device.");
    }
    else {
        // set the keyboard or keypad device to be the first keyboard device or, if no keyboard, the first keypad:
        i = defaultKeyboardIndex;
        if (ndevices == 0) PsychErrorExitMsg(PsychError_user, "Could not detect any keyboard type input devices!");
    }

    // Assign return deviceIndex and deviceRecord for selected i'th entry:
    if (deviceRecord) *deviceRecord = deviceRecords[i];
    return(deviceIndices[i]);
}

void PsychHIDInitializeHIDStandardInterfaces(void)
{
    long KbDeviceUsagePages[NUMDEVICEUSAGES] = {kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_Consumer};
    long KbDeviceUsages[NUMDEVICEUSAGES] = {kHIDUsage_GD_Keyboard, kHIDUsage_GD_Keypad, kHIDUsage_GD_Mouse, kHIDUsage_GD_Pointer, kHIDUsage_GD_Joystick, kHIDUsage_GD_GamePad, kHIDUsage_GD_MultiAxisController, 0x01};
    int  numDeviceUsages = NUMDEVICEUSAGES;
    int  rc, i;

    for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
        queue[i] = NULL;
        psychHIDKbQueueCFRunLoopRef[i] = NULL;
        modifierKeyState[i] = 0;
        KbQueueThread[i] = NULL;
        queueIsAKeyboard[i] = FALSE;
    }

    // Init keyboard queue arrays:
    memset(&psychHIDKbQueueFirstPress[0], 0, sizeof(psychHIDKbQueueFirstPress));
    memset(&psychHIDKbQueueFirstRelease[0], 0, sizeof(psychHIDKbQueueFirstRelease));
    memset(&psychHIDKbQueueLastPress[0], 0, sizeof(psychHIDKbQueueLastPress));
    memset(&psychHIDKbQueueLastRelease[0], 0, sizeof(psychHIDKbQueueLastRelease));
    memset(&psychHIDKbQueueActive[0], 0, sizeof(psychHIDKbQueueActive));
    memset(&psychHIDKbQueueScanKeys[0], 0, sizeof(psychHIDKbQueueScanKeys));

    // Enumerate all possible candidates for default keyboard device:
    PsychHIDGetDeviceListByUsages(numDeviceUsages, KbDeviceUsagePages, KbDeviceUsages, &ndevices, deviceIndices, deviceRecords);

    // Nothing?
    if (ndevices == 0) printf("PTB-WARNING: No keyboard, keypad, mouse, touchpad, joystick etc. button input devices detected! KbQueues and KbCheck won't work.\n");

    // Find default keyboard device. We start with the fail-safe assumption
    // that the first detected keyboard is really a keyboard, then iterate
    // over all candidates, first keyboards, then keypads, mice, pointing
    // devices, joysticks, gamepads. If none is valid, we will fall back to
    // device 0 as best "worst choice".
    //
    // In the past we simply used device 0, the first keyboard, but the iToys
    // company decided to make life difficult again when introducing the
    // MacBook"Pro" 2016 with a new iToy gadget, the oh so magic touchbar. This
    // touchbar enumerates as the first "keyboard" in the system, with no
    // transport, version, manufacturer, product name or serial number, but
    // instead fricking 1046 keys. Luckily this thing also has an invalid
    // locationID of zero, so we can filter it out with this code by only
    // accepting keyboards with a non-zero locationID. The 2nd keyboard is
    // then the actual keyboard of the MacBook"Pro" 2016
    defaultKeyboardIndex = 0;
    for (i = 0; i < ndevices; i++) {
        // Filter out "Logitech Gaming Mouse G502" which is falsely enumerated as keyboard by at least OSX 10.12:
        if (IOHIDDevice_GetVendorID(deviceRecords[i]) == 1133 && IOHIDDevice_GetProductID(deviceRecords[i]) == 49277)
            continue;

        // Filter out / ignore devices with an invalid locationID of zero,
        // like the brain-damaged TouchBar of the MacBook"Pro" 2016, which
        // enumerates as 1st keyboard with 1046 keys:
        if (IOHIDDevice_GetLocationID(deviceRecords[i]) != 0) {
            defaultKeyboardIndex = i;
            break;
        }
    }

    // Create keyboard queue mutex for later use:
    PsychInitMutex(&KbQueueMutex);
    PsychInitCondition(&KbQueueCondition, NULL);

    return;
}

void PsychHIDShutdownHIDStandardInterfaces(void)
{
    int i;

    // Release all keyboard queues:
    for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
        if (psychHIDKbQueueFirstPress[i]) {
            PsychHIDOSKbQueueRelease(i);
        }
    }

    // Release keyboard queue mutex:
    PsychDestroyMutex(&KbQueueMutex);
    PsychDestroyCondition(&KbQueueCondition);

    return;
}

PsychError PsychHIDOSKbCheck(int deviceIndex, double* scanList)
{
    pRecDevice deviceRecord;
    pRecElement currentElement, lastElement = NULL;
    int    debuglevel = 0;
    double *timeValueOutput, *isKeyDownOutput, *keyArrayOutput;
    double dummyKeyDown;
    double dummykeyArrayOutput[256];
    uint32_t usage, usagePage;
    int nout, value;

    // Choose the device index and its record:
    // Want default keyboard?
    if (deviceIndex != INT_MAX) {
        // No. Wants a specific one. Remap negative numbers to positive number
        // and low-level debug output enabled for this invocation:
        if (deviceIndex < 0) {
            debuglevel = 1;
            deviceIndex = -deviceIndex;
        }
    } else {
        // Yes. Wants default keyboard device:
        deviceIndex = -1;
    }

    // Do we finally have a valid keyboard or other suitable input device?
    // PsychHIDOSGetKbQueueDevice() will error out if no suitable device
    // for deviceIndex can be found. Otherwise it will return the HID
    // device record and remapped deviceIndex for use with KbCheck:
    deviceIndex = PsychHIDOSGetKbQueueDevice(deviceIndex, &deviceRecord);

    // Allocate and init out return arguments.
    //
    // Either alloc out the arguments, or redirect to
    // internal dummy variables. This to avoid mxMalloc() call overhead
    // inside the PsychAllocOutXXX() routines:
    nout = PsychGetNumOutputArgs();

    // keyDown flag:
    if (nout >= 1) {
        PsychAllocOutDoubleArg(1, FALSE, &isKeyDownOutput);
    } else {
        isKeyDownOutput = &dummyKeyDown;
    }
    *isKeyDownOutput= (double) FALSE;

    // key state vector:
    if (nout >= 3) {
        PsychAllocOutDoubleMatArg(3, FALSE, 1, 256, 1, &keyArrayOutput);
    }
    else {
        keyArrayOutput = &dummykeyArrayOutput[0];
    }
    memset((void*) keyArrayOutput, 0, sizeof(double) * 256);

    // Query timestamp:
    if (nout >= 2) {
        PsychAllocOutDoubleArg(2, FALSE, &timeValueOutput);

        // Get query timestamp:
        PsychGetPrecisionTimerSeconds(timeValueOutput);
    }

    // Make sure our keyboard query mechanism is not blocked for security reasons, e.g.,
    // secure password entry field active in another process, i.e., EnableSecureEventInput() active.
    if (PsychHIDWarnInputDisabled("PsychHID('KbCheck')")) return(PsychError_none);

    //step through the elements of the device.  Set flags in the return array for down keys.
    for (currentElement = HIDGetFirstDeviceElement(deviceRecord, kHIDElementTypeInput | kHIDElementTypeCollection);
         (currentElement != NULL) && (currentElement != lastElement);
         currentElement = HIDGetNextDeviceElement(currentElement, kHIDElementTypeInput | kHIDElementTypeCollection)) {
        // Keep track of last queried element:
        lastElement = currentElement;

        usage = IOHIDElementGetUsage(currentElement);
        usagePage = IOHIDElementGetUsagePage(currentElement);
        // printf("PTB-DEBUG: [KbCheck]: ce %p page %d usage: %d isArray: %d\n", currentElement, usagePage, usage, IOHIDElementIsArray(currentElement));

        if (IOHIDElementGetType(currentElement) == kIOHIDElementTypeCollection) {
            CFArrayRef children = IOHIDElementGetChildren(currentElement);
            if (!children) continue;

            CFIndex idx, cnt = CFArrayGetCount(children);
            for (idx = 0; idx < cnt; idx++) {
                IOHIDElementRef tIOHIDElementRef = (IOHIDElementRef) CFArrayGetValueAtIndex(children, idx);
                if (tIOHIDElementRef && ((IOHIDElementGetType(tIOHIDElementRef) == kIOHIDElementTypeInput_Button) ||
                                         (IOHIDElementGetType(tIOHIDElementRef) == kIOHIDElementTypeInput_ScanCodes))) {
                    usage = IOHIDElementGetUsage(tIOHIDElementRef);
                    if ((usage <= 256) && (usage >= 1) && ( (scanList == NULL) || (scanList[usage - 1] > 0) )) {
                        value = (int) IOHIDElement_GetValue(tIOHIDElementRef, kIOHIDValueScaleTypePhysical);
                        if (debuglevel > 0) printf("PTB-DEBUG: [KbCheck]: usage: %x value: %d \n", usage, value);
                        keyArrayOutput[usage - 1] = (value || (int) keyArrayOutput[usage - 1]);
                        *isKeyDownOutput = keyArrayOutput[usage - 1] || *isKeyDownOutput;
                    }
                }
            }

            // Done with this currentElement, which was a collection of buttons/keys.
            // Iterate to next currentElement:
            continue;
        }

        // Classic path, or 64-Bit path for non-collection elements:
        if (((usagePage == kHIDPage_KeyboardOrKeypad) || (usagePage == kHIDPage_Button)) && (usage <= 256) && (usage >= 1) &&
            ( (scanList == NULL) || (scanList[usage - 1] > 0) ) ) {
            value = (int) IOHIDElement_GetValue(currentElement, kIOHIDValueScaleTypePhysical);

            if (debuglevel > 0) printf("PTB-DEBUG: [KbCheck]: usage: %x value: %d \n", usage, value);
            keyArrayOutput[usage - 1]=(value || (int) keyArrayOutput[usage - 1]);
            *isKeyDownOutput= keyArrayOutput[usage - 1] || *isKeyDownOutput;
        }
    }

    return(PsychError_none);
}

/*  Table mapping HID usages in the KeyboardOrKeypad page to virtual key codes.
    Names given for keys refer to US layout.
    May be copied freely.
*/

// Taken from this MIT licensed software: <https://github.com/Ahruman/KeyNaming> according
// to above permission note.
/*  KeyNaming.cp
    Keynaming 2.2 implementation
    © 2001-2008 Jens Ayton <jens@ayton.se>, except where otherwise noted.

    Copyright © 2001–2008 Jens Ayton
*/

#define kVKC_Unknown    0xFFFF

enum
{
    /*  Virtual key codes handled specially (from: Inside Macintosh: Text, Appendix C)

    Function keys. These are handled by loading the relevant string IFF the
    character code 10 (kCC_FKey) is generated. This allows the function keys
    to be remapped.

    The Menu/Application key on Windows keyboards also generates the character
    code 10 by default. It isn't immediately clear whether the VKC is used for
    an F key in other circumstances, but apps generally seem to interpret it as
    a normal key. For now, I'm defining it as Menu Key. It's possible Unix-
    oriented keyboards use it with a different key cap, e.g. Meta.
    */
    kVKC_F1        = 122,
    kVKC_F2        = 120,
    kVKC_F3        = 99,
    kVKC_F4        = 118,
    kVKC_F5        = 96,
    kVKC_F6        = 97,
    kVKC_F7        = 98,
    kVKC_F8        = 100,
    kVKC_F9        = 101,
    kVKC_F10      = 109,
    kVKC_F11      = 103,
    kVKC_F12      = 111,
    kVKC_F13      = 105,
    kVKC_F14      = 107,
    kVKC_F15      = 113,
    kVKC_F16      = 106,
    kVKC_Menu      = 110,

    /*  Escape and clear are like the F keys, using the character code kCC_EscClr. */
    kVKC_Esc      = 53,
    kVKC_Clear      = 71,

    /*  The following are handled directly by recognising the virtual code. */
    kVKC_Space      = 49,
    kVKC_CapsLock    = 57,
    kVKC_Shift      = 56,
    kVKC_Option      = 58,
    kVKC_Control    = 59,
    kVKC_rShift      = 60,    /*  Right-hand modifiers; not implemented */
    kVKC_rOption    = 61,
    kVKC_rControl    = 62,
    kVKC_Command    = 55,
    kVKC_Return      = 36,
    kVKC_Backspace    = 51,    /*  Left delete */
    kVKC_Delete      = 117,    /*  right delete */
    kVKC_Help      = 114,
    kVKC_Home      = 115,
    kVKC_PgUp      = 116,
    kVKC_PgDn      = 121,
    kVKC_End      = 119,
    kVKC_LArrow      = 123,
    kVKC_RArrow      = 124,
    kVKC_UArrow      = 126,
    kVKC_DArrow      = 125,
    kVKC_KpdEnter    = 76,    /*  "normal" enter key */
    kVKC_KbdEnter    = 52,    /*  Powerbooks (and some early Macs) */
    kVKC_Fn        = 63,

    /*  Keypad keys. These are named by loading the string "Keypad %@" and
    replacing the %@ with the key's character name. The constant names
    do not correspond to the key caps in any KCHR that I'm aware of;
    they're just used to recognise the set of keys. Note that Enter and
    Clear aren't handled this way. */
    kVKC_Kpd_0      = 81,
    kVKC_Kpd_1      = 75,
    kVKC_Kpd_2      = 67,
    kVKC_Kpd_3      = 89,
    kVKC_Kpd_4      = 91,
    kVKC_Kpd_5      = 92,
    kVKC_Kpd_6      = 78,
    kVKC_Kpd_7      = 86,
    kVKC_Kpd_8      = 87,
    kVKC_Kpd_9      = 88,
    kVKC_Kpd_A      = 69,
    kVKC_Kpd_B      = 83,
    kVKC_Kpd_C      = 84,
    kVKC_Kpd_D      = 85,
    kVKC_Kpd_E      = 82,
    kVKC_Kpd_F      = 65,

    /* 2.1b5: values from new list in Event.h in OS X 10.5 */
    kVKC_VolumeUp    = 72,
    kVKC_VolumeDown    = 73,
    kVKC_Mute      = 74,

    kVKC_F17      = 64,
    kVKC_F18      = 79,
    kVKC_F19      = 80,
    kVKC_F20      = 90,

#if KEYNAMING_ENABLE_HID

    /*  Fake VKCs. These are used for HID usages that, to my knowledge, have
    no corresponding VKC. I use the range 0x6000 up; GetKeys() can't return
    in this range, and other VKC usages probably won't. */
    kFKC_base_      = 0x6000,
    kFKC_rCommand,
    /*  kFKC_Mute,      */
    /*  kFKC_VolumeDown,  */
    /*  kFKC_VolumeUp,    */
    kFKC_Power      = 0x6005,
    kFKC_Eject,
    /*  kFKC_F17,      */
    /*  kFKC_F18,      */
    /*  kFKC_F19,      */
    /*  kFKC_F20,      */
    kFKC_F21      = 0x600B,
    kFKC_F22,
    kFKC_F23,
    kFKC_F24,

#endif

    kEnumAtTheEndToSatisfyIrritableCompilers
};

static const uint16_t  kHID2VKC[] =
{
    kVKC_Unknown,    /* Reserved (no event indicated) */
    kVKC_Unknown,    /* ErrorRollOver */
    kVKC_Unknown,    /* POSTFail */
    kVKC_Unknown,    /* ErrorUndefined */
    0x00,        /* a and A */
    0x0B,        /* b and B */
    0x08,        /* ... */
    0x02,
    0x0E,
    0x03,
    0x05,
    0x04,
    0x22,
    0x26,
    0x28,
    0x25,
    0x2E,
    0x2D,
    0x1F,
    0x23,
    0x0C,
    0x0F,
    0x01,
    0x11,
    0x20,
    0x09,
    0x0D,
    0x07,
    0x10,
    0x06,        /* z and Z */
    0x12,        /* 1 */
    0x13,        /* 2 */
    0x14,        /* ... */
    0x15,
    0x17,
    0x16,
    0x1A,
    0x1C,
    0x19,        /* 9 */
    0x1D,        /* 0 */
    kVKC_Return,    /* Keyboard Return (ENTER) */
    kVKC_Esc,      /* Escape */
    kVKC_Backspace,    /* Delete (Backspace) */
    0x30,        /* Tab */
    kVKC_Space,      /* Space bar */
    0x1B,        /* - and _ */
    0x18,        /* = and + */
    0x21,        /* [ and { */
    0x1E,        /* ] and } */
    0x2A,        /* \ and | */
    kVKC_Unknown,    /* "Non-US # and ~" ?? */
    0x29,        /* ; and : */
    0x27,        /* ' and " */
    0x32,        /* ` and ~ */
    0x2B,        /* , and < */
    0x2F,        /* . and > */
    0x2C,        /* / and ? */
    kVKC_CapsLock,    /* Caps Lock */
    kVKC_F1,      /* F1 */
    kVKC_F2,      /* ... */
    kVKC_F3,
    kVKC_F4,
    kVKC_F5,
    kVKC_F6,
    kVKC_F7,
    kVKC_F8,
    kVKC_F9,
    kVKC_F10,
    kVKC_F11,
    kVKC_F12,      /* F12 */
    kVKC_Unknown,    /* Print Screen */
    kVKC_Unknown,    /* Scroll Lock */
    kVKC_Unknown,    /* Pause */
    kVKC_Unknown,    /* Insert */
    kVKC_Home,      /* Home */
    kVKC_PgUp,      /* Page Up */
    kVKC_Delete,    /* Delete Forward */
    kVKC_End,      /* End */
    kVKC_PgDn,      /* Page Down */
    kVKC_RArrow,    /* Right Arrow */
    kVKC_LArrow,    /* Left Arrow */
    kVKC_DArrow,    /* Down Arrow */
    kVKC_UArrow,    /* Up Arrow */
    kVKC_Clear,      /* Keypad Num Lock and Clear */
    0x4B,        /* Keypad / */
    0x43,        /* Keypad * */
    0x4E,        /* Keypad - */
    0x45,        /* Keypad + */
    kVKC_KpdEnter,    /* Keypad ENTER */
    0x53,        /* Keypad 1 */
    0x54,        /* Keypad 2 */
    0x55,        /* Keypad 3 */
    0x56,        /* Keypad 4 */
    0x57,        /* Keypad 5 */
    0x58,        /* Keypad 6 */
    0x59,        /* Keypad 7 */
    0x5B,        /* Keypad 8 */
    0x5C,        /* Keypad 9 */
    0x52,        /* Keypad 0 */
    0x41,        /* Keypad . */
    0x32,        /* "Keyboard Non-US \ and |" */
    kVKC_Unknown,    /* "Keyboard Application" (Windows key for Windows 95, and "Compose".) */
    kVKC_Unknown,    /* Keyboard Power (status, not key... but Apple doesn't seem to have read the spec properly) */
    0x51,        /* Keypad = */
    kVKC_F13,      /* F13 */
    kVKC_F14,      /* ... */
    kVKC_F15,
    kVKC_F16,
    kVKC_F17,
    kVKC_F18,
    kVKC_F19,
    kVKC_F20,
    kVKC_Unknown,
    kVKC_Unknown,
    kVKC_Unknown,
    kVKC_Unknown,      /* F24 */
    kVKC_Unknown,    /* Keyboard Execute */
    kVKC_Help,      /* Keyboard Help */
    kVKC_Unknown,    /* Keyboard Menu */
    kVKC_Unknown,    /* Keyboard Select */
    kVKC_Unknown,    /* Keyboard Stop */
    kVKC_Unknown,    /* Keyboard Again */
    kVKC_Unknown,    /* Keyboard Undo */
    kVKC_Unknown,    /* Keyboard Cut */
    kVKC_Unknown,    /* Keyboard Copy */
    kVKC_Unknown,    /* Keyboard Paste */
    kVKC_Unknown,    /* Keyboard Find */
    kVKC_Mute,      /* Keyboard Mute */
    kVKC_VolumeUp,    /* Keyboard Volume Up */
    kVKC_VolumeDown,  /* Keyboard Volume Down */
    kVKC_CapsLock,    /* Keyboard Locking Caps Lock */
    kVKC_Unknown,    /* Keyboard Locking Num Lock */
    kVKC_Unknown,    /* Keyboard Locking Scroll Lock */
    0x41,        /*  Keypad Comma ("Keypad Comma is the appropriate usage for the Brazilian
                        keypad period (.) key. This represents the closest possible  match, and
                        system software should do the correct mapping based on the current locale
                        setting." If strange stuff happens on a (physical) Brazilian keyboard,
                        I'd like to know about it. */
    0x51,        /* Keypad Equal Sign ("Used on AS/400 Keyboards.") */
    kVKC_Unknown,    /* Keyboard International1 (Brazilian / and ? key? Kanji?) */
    kVKC_Unknown,    /* Keyboard International2 (Kanji?) */
    kVKC_Unknown,    /* Keyboard International3 (Kanji?) */
    kVKC_Unknown,    /* Keyboard International4 (Kanji?) */
    kVKC_Unknown,    /* Keyboard International5 (Kanji?) */
    kVKC_Unknown,    /* Keyboard International6 (Kanji?) */
    kVKC_Unknown,    /* Keyboard International7 (Kanji?) */
    kVKC_Unknown,    /* Keyboard International8 (Kanji?) */
    kVKC_Unknown,    /* Keyboard International9 (Kanji?) */
    kVKC_Unknown,    /* Keyboard LANG1 (Hangul/English toggle) */
    kVKC_Unknown,    /* Keyboard LANG2 (Hanja conversion key) */
    kVKC_Unknown,    /* Keyboard LANG3 (Katakana key) */    // kVKC_Kana?
    kVKC_Unknown,    /* Keyboard LANG4 (Hirigana key) */
    kVKC_Unknown,    /* Keyboard LANG5 (Zenkaku/Hankaku key) */
    kVKC_Unknown,    /* Keyboard LANG6 */
    kVKC_Unknown,    /* Keyboard LANG7 */
    kVKC_Unknown,    /* Keyboard LANG8 */
    kVKC_Unknown,    /* Keyboard LANG9 */
    kVKC_Unknown,    /* Keyboard Alternate Erase ("Example, Erase-Eaze™ key.") */
    kVKC_Unknown,    /* Keyboard SysReq/Attention */
    kVKC_Unknown,    /* Keyboard Cancel */
    kVKC_Unknown,    /* Keyboard Clear */
    kVKC_Unknown,    /* Keyboard Prior */
    kVKC_Unknown,    /* Keyboard Return */
    kVKC_Unknown,    /* Keyboard Separator */
    kVKC_Unknown,    /* Keyboard Out */
    kVKC_Unknown,    /* Keyboard Oper */
    kVKC_Unknown,    /* Keyboard Clear/Again */
    kVKC_Unknown,    /* Keyboard CrSel/Props */
    kVKC_Unknown,    /* Keyboard ExSel */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Keypad 00 */
    kVKC_Unknown,    /* Keypad 000 */
    kVKC_Unknown,    /* Thousands Separator */
    kVKC_Unknown,    /* Decimal Separator */
    kVKC_Unknown,    /* Currency Unit */
    kVKC_Unknown,    /* Currency Sub-unit */
    kVKC_Unknown,    /* Keypad ( */
    kVKC_Unknown,    /* Keypad ) */
    kVKC_Unknown,    /* Keypad { */
    kVKC_Unknown,    /* Keypad } */
    kVKC_Unknown,    /* Keypad Tab */
    kVKC_Unknown,    /* Keypad Backspace */
    kVKC_Unknown,    /* Keypad A */
    kVKC_Unknown,    /* Keypad B */
    kVKC_Unknown,    /* Keypad C */
    kVKC_Unknown,    /* Keypad D */
    kVKC_Unknown,    /* Keypad E */
    kVKC_Unknown,    /* Keypad F */
    kVKC_Unknown,    /* Keypad XOR */
    kVKC_Unknown,    /* Keypad ^ */
    kVKC_Unknown,    /* Keypad % */
    kVKC_Unknown,    /* Keypad < */
    kVKC_Unknown,    /* Keypad > */
    kVKC_Unknown,    /* Keypad & */
    kVKC_Unknown,    /* Keypad && */
    kVKC_Unknown,    /* Keypad | */
    kVKC_Unknown,    /* Keypad || */
    kVKC_Unknown,    /* Keypad : */
    kVKC_Unknown,    /* Keypad # */
    kVKC_Unknown,    /* Keypad Space */
    kVKC_Unknown,    /* Keypad @ */
    kVKC_Unknown,    /* Keypad ! */
    kVKC_Unknown,    /* Keypad Memory Store */
    kVKC_Unknown,    /* Keypad Memory Recall */
    kVKC_Unknown,    /* Keypad Memory Clear */
    kVKC_Unknown,    /* Keypad Memory Add */
    kVKC_Unknown,    /* Keypad Memory Subtract */
    kVKC_Unknown,    /* Keypad Memory Multiply */
    kVKC_Unknown,    /* Keypad Memory Divide */
    kVKC_Unknown,    /* Keypad +/- */
    kVKC_Unknown,    /* Keypad Clear */
    kVKC_Unknown,    /* Keypad Clear Entry */
    kVKC_Unknown,    /* Keypad Binary */
    kVKC_Unknown,    /* Keypad Octal */
    kVKC_Unknown,    /* Keypad Decimal */
    kVKC_Unknown,    /* Keypad Hexadecimal */
    kVKC_Unknown,    /* Reserved */
    kVKC_Unknown,    /* Reserved */
    kVKC_Control,    /* Keyboard LeftControl */
    kVKC_Shift,      /* Keyboard LeftShift */
    kVKC_Option,    /* Keyboard LeftAlt */
    kVKC_Command,    /* Keyboard LeftGUI */
    kVKC_rControl,    /* Keyboard RightControl */
    kVKC_rShift,    /* Keyboard RightShift */
    kVKC_rOption,    /* Keyboard RightAlt */
    kVKC_Unknown    /* Keyboard RightGUI */
};
enum { kHID2VKCSize = sizeof kHID2VKC / sizeof kHID2VKC[0] };

// Convert uint64 time in nanoseconds to double value in GetSecs() seconds:
static double convertTime(uint64_t timeUInt64) {
    double timeDouble = (double) timeUInt64;
    return timeDouble / 1000000000;
}

static void PsychHIDKbQueueCallbackFunction(void *target, IOReturn result, void *sender)
{
    // This routine is executed each time the queue transitions from empty to non-empty
    // The CFRunLoop of the thread in KbQueueWorkerThreadMain() is the one that executes here:
    IOHIDQueueRef queue = (IOHIDQueueRef) sender;
    IOHIDValueRef valueRef = NULL;
    int deviceIndex = (int) target;
    double timestamp;
    int eventValue;
    long keysUsage = -1;
    PsychHIDEventRecord evt;

    result=kIOReturnError;
    if (!queue) return; // Nothing we can do because we can't access queue, (shouldn't happen)

    while (1) {
        // This function only gets called when queue transitions from empty to non-empty
        // Therefore, we must process all available events in this while loop before
        // it will be possible for this function to be notified again.
        if (valueRef) {
            CFRelease(valueRef);
            valueRef = NULL;
        }

        // Dequeue next event from queue in a polling non-blocking fashion:
        valueRef = IOHIDQueueCopyNextValueWithTimeout(queue, 0.0);

        // Done? Exit, if so:
        if (!valueRef) break;

        // Get event value, e.g., the key state of a key or button 1 = pressed, 0 = released:
        eventValue = IOHIDValueGetIntegerValue(valueRef);

        // Get usage value, ie., the identity of the key:
        IOHIDElementRef element = IOHIDValueGetElement(valueRef);
        keysUsage = IOHIDElementGetUsage(element);

        // Get double GetSecs timestamp, computed from returned uint64 nanoseconds timestamp:
        timestamp = convertTime(IOHIDValueGetTimeStamp(valueRef));

        // Don't bother with keysUsage of 0 (meaningless) or 1 (ErrorRollOver) for keyboards:
        if ((queueIsAKeyboard[deviceIndex]) && (keysUsage <= 1)) continue;

        // Clear ringbuffer event:
        memset(&evt, 0 , sizeof(evt));

        // Cooked key code defaults to "unhandled", and stays that way for anything but keyboards:
        evt.cookedEventCode = -1;

        // For real keyboards we can compute cooked key codes: Requires OSX 10.5 or later.
        if (queueIsAKeyboard[deviceIndex]) {
            // Keyboard(ish) device. We can handle this under some conditions.
            // Init to a default of handled, but unmappable/ignored keycode:
            evt.cookedEventCode = 0;

            // Keypress event code available in mapping table?
            if (keysUsage < kHID2VKCSize) {
                // Yes: We try to map this to a character code:

                // Step 1: Map HID usage value to virtual keycode via LUT:
                uint16_t vcKey = kHID2VKC[keysUsage];

                // Keep track of SHIFT keys as modifier keys: Bits 0 == Command, 1 == Shift, 2 == CapsLock, 3 == Alt/Option, 4 == CTRL
                if ((vcKey == kVKC_Shift || vcKey == kVKC_rShift) && (eventValue != 0)) modifierKeyState[deviceIndex] |=  (1 << 1);
                if ((vcKey == kVKC_Shift || vcKey == kVKC_rShift) && (eventValue == 0)) modifierKeyState[deviceIndex] &= ~(1 << 1);

                // Keep track of ALT keys as modifier keys:
                if ((vcKey == kVKC_Option || vcKey == kVKC_rOption) && (eventValue != 0)) modifierKeyState[deviceIndex] |=  (1 << 3);
                if ((vcKey == kVKC_Option || vcKey == kVKC_rOption) && (eventValue == 0)) modifierKeyState[deviceIndex] &= ~(1 << 3);

                // Keep track of CTRL keys as modifier keys:
                if ((vcKey == kVKC_Control || vcKey == kVKC_rControl) && (eventValue != 0)) modifierKeyState[deviceIndex] |=  (1 << 4);
                if ((vcKey == kVKC_Control || vcKey == kVKC_rControl) && (eventValue == 0)) modifierKeyState[deviceIndex] &= ~(1 << 4);

                // Was this a CTRL + C interrupt request?
                if ((eventValue != 0) && (vcKey == 0x08) && (modifierKeyState[deviceIndex] & (1 << 4))) {
                    // Yes: Tell the console input helper about it, so it can send interrupt
                    // signals to the runtime and reenable keyboard input if appropriate:
                    // Note: Not sure if the mutex exclusion is needed here, but better safe than sorry.
                    PsychLockMutex(&KbQueueMutex);
                    ConsoleInputHelper(-1);
                    PsychUnlockMutex(&KbQueueMutex);
                }

                // Key press?
                if (eventValue != 0) {
                    // Step 2: Translate virtual key code into unicode char:
                    // Ok, this is the usual horrifying complexity of Apple's system. We use code
                    // snippets found on StackOverflow, modified to suit our needs, e.g., we track
                    // modifier keys manually, at least left and right ALT and SHIFT keys. We don't
                    // care about other modifiers.
                    if (keyboardLayout) {
                        UInt32 deadKeyState = 0;
                        UniCharCount maxStringLength = 255;
                        UniCharCount actualStringLength = 0;
                        UniChar unicodeString[maxStringLength];

                        OSStatus status = UCKeyTranslate(keyboardLayout,
                                                        vcKey, kUCKeyActionDown, modifierKeyState[deviceIndex],
                                                        kbdType, 0,
                                                        &deadKeyState,
                                                        maxStringLength,
                                                        &actualStringLength, unicodeString);

                        if ((actualStringLength == 0) && deadKeyState) {
                            status = UCKeyTranslate(keyboardLayout,
                                                    kVK_Space, kUCKeyActionDown, 0,
                                                    kbdType, 0,
                                                    &deadKeyState,
                                                    maxStringLength,
                                                    &actualStringLength, unicodeString);
                        }

                        if((actualStringLength > 0) && (status == noErr)) {
                            // Assign final cooked / mapped keycode:
                            evt.cookedEventCode = (int) unicodeString[0];

                            // Send same keystroke character to console input helper.
                            // In kbqueue-based ListenChar(1) mode, the helper will
                            // inject/forward the character into the runtime:
                            // Note: ConsoleInputHelper() should be safe to call without
                            // mutex protection for >= 0 event codes.
                            ConsoleInputHelper(evt.cookedEventCode);
                        }
                    }
                }
            }
        }

        PsychLockMutex(&KbQueueMutex);

        // Update records of first and latest key presses and releases
        if (eventValue != 0) {
            if (psychHIDKbQueueFirstPress[deviceIndex]) {
                // First key press timestamp:
                if (psychHIDKbQueueFirstPress[deviceIndex][keysUsage-1] == 0) {
                    psychHIDKbQueueFirstPress[deviceIndex][keysUsage-1] = timestamp;
                }
            }

            if (psychHIDKbQueueLastPress[deviceIndex]) {
                // Last key press timestamp:
                psychHIDKbQueueLastPress[deviceIndex][keysUsage-1] = timestamp;
            }
            evt.status |= (1 << 0);
        }
        else {
            if (psychHIDKbQueueFirstRelease[deviceIndex]) {
                // First key release timestamp:
                if (psychHIDKbQueueFirstRelease[deviceIndex][keysUsage-1] == 0) psychHIDKbQueueFirstRelease[deviceIndex][keysUsage-1] = timestamp;
            }

            if (psychHIDKbQueueLastRelease[deviceIndex]) {
                // Last key release timestamp:
                psychHIDKbQueueLastRelease[deviceIndex][keysUsage-1] = timestamp;
            }
            evt.status &= ~(1 << 0);
        }

        // Update event buffer:
        evt.timestamp = timestamp;
        evt.rawEventCode = keysUsage;
        PsychHIDAddEventToEventBuffer(deviceIndex, &evt);

        // Tell waiting userspace (under KbQueueMutxex protection for better scheduling) something interesting has changed:
        PsychSignalCondition(&KbQueueCondition);

        PsychUnlockMutex(&KbQueueMutex);

        // Next while loop iteration to dequeue potentially more events:
    }

    // Done for this queue transition. Return to runloop.
}

// Async processing thread for keyboard events:
static void *KbQueueWorkerThreadMain(void *inarg) {
    int deviceIndex = (int) inarg;
    int rc;

    // Assign a name to ourselves, for debugging:
    PsychSetThreadName("PsychHIDKbQueue");

    // Switch ourselves (NULL) to RT scheduling: We promise to use / require at most (0+1) == 1 msec every
    // 10 msecs and allow for wakeup delay/jitter of up to 2 msecs -- perfectly reasonable, given that we
    // only do minimal << 1 msec processing, only at the timescale of human reaction times, and driven by
    // input devices with at least 4+/-4 msecs jitter at 8 msec USB polling frequency.
    if ((rc = PsychSetThreadPriority(NULL, 2, 0)) > 0) {
        printf("PsychHID: KbQueueCreate: Failed to switch to realtime priority [%s].\n", strerror(rc));
    }

    // Keep a global reference to the runloop, as we need it in KbQueueRelease to get this thread to exit:
    psychHIDKbQueueCFRunLoopRef[deviceIndex] = (CFRunLoopRef) CFRunLoopGetCurrent();
    CFRetain(psychHIDKbQueueCFRunLoopRef[deviceIndex]);

    // Add HID queue to current runloop:
    IOHIDQueueScheduleWithRunLoop(queue[deviceIndex], psychHIDKbQueueCFRunLoopRef[deviceIndex], kCFRunLoopDefaultMode);

    // Poor man's memory barrier:
    PsychLockMutex(&KbQueueMutex);
    PsychUnlockMutex(&KbQueueMutex);

    // Start the run loop, code execution will block here until run loop is stopped again by PsychHIDKbQueueRelease
    // The run loop will be responsible for executing the code in PsychHIDKbQueueCallbackFunction() whenever input
    // is available:
    CFRunLoopRun();

    // Remove HID queue from current runloop:
    IOHIDQueueUnscheduleFromRunLoop(queue[deviceIndex], psychHIDKbQueueCFRunLoopRef[deviceIndex], kCFRunLoopDefaultMode);

    // Done. Die peacefully:
    return(NULL);
}

int PsychHIDGetDefaultKbQueueDevice(void)
{
    // Return the default keyboard or keypad device as the first keyboard device or, if no keyboard, the first keypad:
    if (ndevices == 0) PsychErrorExitMsg(PsychError_user, "No KbQueue capable input devices detected! Game over.");
    return(deviceIndices[defaultKeyboardIndex]);
}

// Put element into the dictionary and into the queue:
PsychError PsychHIDOSKbElementAdd(IOHIDElementRef element, IOHIDQueueRef queue, int deviceIndex)
{
    // If at least one keyboard style device is detected, mark this queue as keyboard queue:
    if (IOHIDElementGetUsagePage(element) == kHIDPage_KeyboardOrKeypad) queueIsAKeyboard[deviceIndex] = TRUE;

    // Avoid redundant assignment to same keycode:
    if (IOHIDQueueContainsElement(queue, element)) {
        if (getenv("PSYCHHID_TELLME")) printf("--> Key %i Already assigned --> Skipping.\n", IOHIDElementGetUsage(element) - 1);
        return(PsychError_none);
    }

    if (getenv("PSYCHHID_TELLME")) {
        printf("--> Accepting key %i as new KbQueue element%s.\n", IOHIDElementGetUsage(element) - 1, (queueIsAKeyboard[deviceIndex]) ? " for a keyboard" : "");
    }

    // Put the element cookie into the queue:
    IOHIDQueueAddElement(queue, element);

    return(PsychError_none);
}

PsychError PsychHIDOSKbQueueCreate(int deviceIndex, int numScankeys, int* scanKeys, int numValuators, int numSlots, unsigned int flags, unsigned int windowHandle)
{
    pRecDevice deviceRecord;
    psych_bool verbose = getenv("PSYCHHID_TELLME") != NULL;
    // Initialize keyboardLayout and kbdType, used by the keyboard queue thread for mapping keycodes to cooked characters:
    TISInputSourceRef currentKeyboard = TISCopyCurrentKeyboardInputSource();
    CFDataRef uchr = (CFDataRef) ((currentKeyboard) ? TISGetInputSourceProperty(currentKeyboard, kTISPropertyUnicodeKeyLayoutData) : NULL);
    keyboardLayout = (const UCKeyboardLayout*) ((uchr) ? CFDataGetBytePtr(uchr) : NULL);
    kbdType = LMGetKbdType();

    (void) windowHandle;

    // Valid number of keys?
    if (scanKeys && (numScankeys != 256)) {
        PsychErrorExitMsg(PsychError_user, "Second argument to KbQueueCreate must be a vector with 256 elements.");
    }

    if (numValuators > 0)
        PsychErrorExitMsg(PsychError_unimplemented, "Valuators are not supported on macOS.");

    // Do we finally have a valid keyboard or other suitable input device?
    // PsychHIDOSGetKbQueueDevice() will error out if no suitable device
    // for deviceIndex can be found. Otherwise it will return the HID
    // device record and remapped deviceIndex for use with our KbQueues:
    deviceIndex = PsychHIDOSGetKbQueueDevice(deviceIndex, &deviceRecord);

    // Keyboard queue for this deviceIndex already created?
    if (psychHIDKbQueueFirstPress[deviceIndex]) {
        // Yep. Release it, so we can start from scratch:
        PsychHIDOSKbQueueRelease(deviceIndex);
    }

    // Allocate and zero-init memory for tracking key presses and key releases:
    psychHIDKbQueueFirstPress[deviceIndex]   = calloc(256, sizeof(double));
    psychHIDKbQueueFirstRelease[deviceIndex] = calloc(256, sizeof(double));
    psychHIDKbQueueLastPress[deviceIndex]    = calloc(256, sizeof(double));
    psychHIDKbQueueLastRelease[deviceIndex]  = calloc(256, sizeof(double));
    psychHIDKbQueueScanKeys[deviceIndex]     = calloc(256, sizeof(int));

    // Assign scanKeys vector, if any:
    if (scanKeys) {
        // Copy it:
        memcpy(psychHIDKbQueueScanKeys[deviceIndex], scanKeys, 256 * sizeof(int));
    } else {
        // None provided. Enable all keys by default:
        memset(psychHIDKbQueueScanKeys[deviceIndex], 1, 256 * sizeof(int));
    }

    // Create HIDQueue for device:
    queue[deviceIndex] = IOHIDQueueCreate(kCFAllocatorDefault, deviceRecord, 30, 0);
    if (NULL == queue[deviceIndex]) {
        PsychHIDOSKbQueueRelease(deviceIndex);
        PsychErrorExitMsg(PsychError_system, "Failed to create event queue for detecting key press.");
    }

    // Mark as a non-keyboard device, to start with:
    queueIsAKeyboard[deviceIndex] = FALSE;

    // Parse HID device to add all detected and selected keys:
    {
        // Add deviceRecord's elements to our queue, filtering unwanted keys via 'scanList'.
        // This code is almost identical to the enumeration code in PsychHIDKbCheck, to make sure we get
        // matching performance and behaviour and hopefully that it works on the latest problematic Apple
        // hardware, e.g., late 2013 MacBookAir and OSX 10.9:
        {
            uint32_t usage, usagePage;
            pRecElement currentElement, lastElement = NULL;

            // Step through the elements of the device and add matching ones:
            for (currentElement = HIDGetFirstDeviceElement(deviceRecord, kHIDElementTypeInput | kHIDElementTypeCollection);
                (currentElement != NULL) && (currentElement != lastElement);
                currentElement = HIDGetNextDeviceElement(currentElement, kHIDElementTypeInput | kHIDElementTypeCollection))
            {
                // Keep track of last queried element:
                lastElement = currentElement;

                usage     = IOHIDElementGetUsage(currentElement);
                usagePage = IOHIDElementGetUsagePage(currentElement);
                if (verbose) {
                    printf("PTB-DEBUG: [KbQueueCreate]: ce %p page %d usage: %d isArray: %d\n", currentElement, usagePage, usage, IOHIDElementIsArray(currentElement));
                }

                if (IOHIDElementGetType(currentElement) == kIOHIDElementTypeCollection) {
                    CFArrayRef children = IOHIDElementGetChildren(currentElement);
                    if (!children) continue;

                    CFIndex idx, cnt = CFArrayGetCount(children);
                    if (verbose) {
                        printf("PTB-DEBUG: [KbQueueCreate]: ce %p Collection with %d children:\n", currentElement, cnt);
                    }

                    for (idx = 0; idx < cnt; idx++) {
                        IOHIDElementRef tIOHIDElementRef = (IOHIDElementRef) CFArrayGetValueAtIndex(children, idx);
                        if (tIOHIDElementRef && ((IOHIDElementGetType(tIOHIDElementRef) == kIOHIDElementTypeInput_Button) ||
                                                (IOHIDElementGetType(tIOHIDElementRef) == kIOHIDElementTypeInput_ScanCodes))) {
                            usage = IOHIDElementGetUsage(tIOHIDElementRef);
                            usagePage = IOHIDElementGetUsagePage(tIOHIDElementRef);
                            if (verbose) {
                                printf("PTB-DEBUG: [KbQueueCreate]: ce %p child %d page %d usage: %d isArray: %d\n", currentElement, idx, usagePage, usage, IOHIDElementIsArray(tIOHIDElementRef));
                            }
                            if ((usage <= 256) && (usage >= 1) && ( (scanKeys == NULL) || (scanKeys[usage - 1] > 0) )) {
                                // Add it for use in keyboard queue:
                                PsychHIDOSKbElementAdd(tIOHIDElementRef, queue[deviceIndex], deviceIndex);
                            }
                        }
                    }

                    // Done with this currentElement, which was a collection of buttons/keys.
                    // Iterate to next currentElement:
                    continue;
                }

                // Classic path for non-collection elements:
                if(((usagePage == kHIDPage_KeyboardOrKeypad) || (usagePage == kHIDPage_Button)) && (usage <= 256) && (usage >= 1) &&
                ( (scanKeys == NULL) || (scanKeys[usage - 1] > 0) ) ) {
                    // Add it for use in keyboard queue:
                    PsychHIDOSKbElementAdd(currentElement, queue[deviceIndex], deviceIndex);
                }
            }
        }
    }

    // Register "queue empty -> non-empty transition" callback:
    IOHIDQueueRegisterValueAvailableCallback(queue[deviceIndex], (IOHIDCallback) PsychHIDKbQueueCallbackFunction, (void*) (long) deviceIndex);

    // Create event buffer:
    if (!PsychHIDCreateEventBuffer(deviceIndex, numValuators, numSlots)) {
        PsychHIDOSKbQueueRelease(deviceIndex);
        PsychErrorExitMsg(PsychError_system, "Failed to create keyboard queue for detecting key press.");
    }

    // Start the processing thread for this queue:
    if (PsychCreateThread(&KbQueueThread[deviceIndex], NULL, KbQueueWorkerThreadMain, (void*) (long) deviceIndex)) {
        KbQueueThread[deviceIndex] = NULL;
        PsychHIDOSKbQueueRelease(deviceIndex);

        printf("PsychHID-ERROR: Creation of keyboard queue thread for deviceIndex %i failed!\n", deviceIndex);
        PsychErrorExitMsg(PsychError_system, "Creation of keyboard queue background processing thread failed!");
    }

    // Thread successfully created. Wait for it to init psychHIDKbQueueCFRunLoopRef[] and then to actually
    // enter the runloop. This to avoid races if KbQueueCreate -> KbQueueRelease is called back-to-back.
    // This is not exactly textbook style optimal way to sync the threads, but efficiency doesn't matter
    // here, whereas not adding extra complexity for the "proper way" of interacting with Apples nightmare does:
    while ((!psychHIDKbQueueCFRunLoopRef[deviceIndex]) || (!CFRunLoopIsWaiting(psychHIDKbQueueCFRunLoopRef[deviceIndex]))) {
        // Poor man's memory barrier:
        PsychLockMutex(&KbQueueMutex);
        PsychUnlockMutex(&KbQueueMutex);

        // Wait a msec, thread should be there any moment:
        PsychYieldIntervalSeconds(0.001);
    }

    // Ready to use this keybord queue.
    return(PsychError_none);
}

void PsychHIDOSKbQueueRelease(int deviceIndex)
{
    // Get true keyboardqueue index assigned to deviceIndex from original user provided deviceIndex:
    deviceIndex = PsychHIDOSGetKbQueueDevice(deviceIndex, NULL);

    // Keyboard queue for this deviceIndex already exists?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Nothing to do then.
        return;
    }

    // Ok, we have a keyboard queue. Stop any operation on it first:
    PsychHIDOSKbQueueStop(deviceIndex);

    if (KbQueueThread[deviceIndex]) {
        // Not quite needed if things go to plan, but be extra cautious and make absolutely sure there actually is
        // a running CFRunLoop() which can respond to the CFRunLoopStop call below:
        while ((!psychHIDKbQueueCFRunLoopRef[deviceIndex]) || (!CFRunLoopIsWaiting(psychHIDKbQueueCFRunLoopRef[deviceIndex]))) {
            // Poor man's memory barrier:
            PsychLockMutex(&KbQueueMutex);
            PsychUnlockMutex(&KbQueueMutex);

            // Wait a msec, thread should be there any moment:
            PsychYieldIntervalSeconds(0.001);
        }

        // Stop the CFRunLoop, which will allow its associated thread to exit:
        CFRunLoopStop(psychHIDKbQueueCFRunLoopRef[deviceIndex]);

        // Shutdown the thread, wait for its termination:
        PsychDeleteThread(&KbQueueThread[deviceIndex]);
        KbQueueThread[deviceIndex] = NULL;

        // Release the CFRunLoop for this queue:
        CFRelease(psychHIDKbQueueCFRunLoopRef[deviceIndex]);
        psychHIDKbQueueCFRunLoopRef[deviceIndex] = NULL;
    }

    // Release queue object:
    if (queue[deviceIndex]) CFRelease(queue[deviceIndex]);
    queue[deviceIndex] = NULL;

    // Release its data structures:
    free(psychHIDKbQueueFirstPress[deviceIndex]); psychHIDKbQueueFirstPress[deviceIndex] = NULL;
    free(psychHIDKbQueueFirstRelease[deviceIndex]); psychHIDKbQueueFirstRelease[deviceIndex] = NULL;
    free(psychHIDKbQueueLastPress[deviceIndex]); psychHIDKbQueueLastPress[deviceIndex] = NULL;
    free(psychHIDKbQueueLastRelease[deviceIndex]); psychHIDKbQueueLastRelease[deviceIndex] = NULL;
    free(psychHIDKbQueueScanKeys[deviceIndex]); psychHIDKbQueueScanKeys[deviceIndex] = NULL;

    // Release kbqueue event buffer:
    PsychHIDDeleteEventBuffer(deviceIndex);

    // Done.
    return;
}

void PsychHIDOSKbQueueStop(int deviceIndex)
{
    // Get true keyboardqueue index assigned to deviceIndex from original user provided deviceIndex:
    deviceIndex = PsychHIDOSGetKbQueueDevice(deviceIndex, NULL);

    // Keyboard queue for this deviceIndex already exists?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Nothing to do then.
        return;
    }

    // Keyboard queue already stopped?
    if (!psychHIDKbQueueActive[deviceIndex]) return;

    // Stop event collection in the queue:
    IOHIDQueueStop(queue[deviceIndex]);

    // Mark queue logically stopped:
    psychHIDKbQueueActive[deviceIndex] = FALSE;

    return;
}

void PsychHIDOSKbQueueStart(int deviceIndex)
{
    // Get true keyboardqueue index assigned to deviceIndex from original user provided deviceIndex:
    deviceIndex = PsychHIDOSGetKbQueueDevice(deviceIndex, NULL);

    // Does Keyboard queue for this deviceIndex already exist?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Bad bad...
        printf("PsychHID-ERROR: Tried to start processing on non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No queue for that device yet!");
    }

    // Keyboard queue already stopped? Then we ain't nothing to do:
    if (psychHIDKbQueueActive[deviceIndex]) return;

    // Clear old content:
    PsychHIDOSKbQueueFlush(deviceIndex);

    // Start event collection in the queue:
    IOHIDQueueStart(queue[deviceIndex]);

    // Mark this queue as logically started:
    psychHIDKbQueueActive[deviceIndex] = TRUE;

    return;
}

void PsychHIDOSKbQueueFlush(int deviceIndex)
{
    // Get true keyboardqueue index assigned to deviceIndex from original user provided deviceIndex:
    deviceIndex = PsychHIDOSGetKbQueueDevice(deviceIndex, NULL);

    // Does Keyboard queue for this deviceIndex already exist?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Bad bad...
        printf("PsychHID-ERROR: Tried to flush non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No queue for that device yet!");
    }

    // Clear out current state for this queue:
    PsychLockMutex(&KbQueueMutex);
    memset(psychHIDKbQueueFirstPress[deviceIndex]   , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueFirstRelease[deviceIndex] , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueLastPress[deviceIndex]    , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueLastRelease[deviceIndex]  , 0, (256 * sizeof(double)));
    modifierKeyState[deviceIndex] = 0;
    PsychUnlockMutex(&KbQueueMutex);

    return;
}

void PsychHIDOSKbQueueCheck(int deviceIndex)
{
    double *hasKeyBeenDownOutput, *firstPressTimeOutput, *firstReleaseTimeOutput, *lastPressTimeOutput, *lastReleaseTimeOutput;
    psych_bool isFirstPressSpecified, isFirstReleaseSpecified, isLastPressSpecified, isLastReleaseSpecified;
    int i;

    // Get true keyboardqueue index assigned to deviceIndex from original user provided deviceIndex:
    deviceIndex = PsychHIDOSGetKbQueueDevice(deviceIndex, NULL);

    // Does Keyboard queue for this deviceIndex already exist?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Bad bad...
        printf("PsychHID-ERROR: Tried to check non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No queue for that device yet!");
    }

    // Allocate output
    PsychAllocOutDoubleArg(1, FALSE, &hasKeyBeenDownOutput);
    isFirstPressSpecified = PsychAllocOutDoubleMatArg(2, FALSE, 1, 256, 1, &firstPressTimeOutput);
    isFirstReleaseSpecified = PsychAllocOutDoubleMatArg(3, FALSE, 1, 256, 1, &firstReleaseTimeOutput);
    isLastPressSpecified = PsychAllocOutDoubleMatArg(4, FALSE, 1, 256, 1, &lastPressTimeOutput);
    isLastReleaseSpecified = PsychAllocOutDoubleMatArg(5, FALSE, 1, 256, 1, &lastReleaseTimeOutput);

    // Initialize output
    if(isFirstPressSpecified) memset((void*) firstPressTimeOutput, 0, sizeof(double) * 256);
    if(isFirstReleaseSpecified) memset((void*) firstReleaseTimeOutput, 0, sizeof(double) * 256);
    if(isLastPressSpecified) memset((void*) lastPressTimeOutput, 0, sizeof(double) * 256);
    if(isLastReleaseSpecified) memset((void*) lastReleaseTimeOutput, 0, sizeof(double) * 256);

    *hasKeyBeenDownOutput=0;

    // Compute and assign output:
    PsychLockMutex(&KbQueueMutex);

    for (i = 0; i < 256; i++) {
        double lastRelease  = psychHIDKbQueueLastRelease[deviceIndex][i];
        double lastPress    = psychHIDKbQueueLastPress[deviceIndex][i];
        double firstRelease = psychHIDKbQueueFirstRelease[deviceIndex][i];
        double firstPress   = psychHIDKbQueueFirstPress[deviceIndex][i];

        if (firstPress) {
            *hasKeyBeenDownOutput=1;
            if(isFirstPressSpecified) firstPressTimeOutput[i] = firstPress;
            psychHIDKbQueueFirstPress[deviceIndex][i] = 0;
        }

        if (firstRelease) {
            if(isFirstReleaseSpecified) firstReleaseTimeOutput[i] = firstRelease;
            psychHIDKbQueueFirstRelease[deviceIndex][i] = 0;
        }

        if (lastPress) {
            if(isLastPressSpecified) lastPressTimeOutput[i] = lastPress;
            psychHIDKbQueueLastPress[deviceIndex][i] = 0;
        }

        if (lastRelease) {
            if(isLastReleaseSpecified) lastReleaseTimeOutput[i] = lastRelease;
            psychHIDKbQueueLastRelease[deviceIndex][i] = 0;
        }
    }

    PsychUnlockMutex(&KbQueueMutex);

    return;
}

void PsychHIDOSKbTriggerWait(int deviceIndex, int numScankeys, int* scanKeys)
{
    PsychErrorExitMsg(PsychError_unimplemented, "PsychHID('KbTriggerWait') is no longer implemented. Use KbTriggerWait() instead!");
}
