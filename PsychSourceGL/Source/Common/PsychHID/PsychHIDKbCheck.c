/*
  PsychToolbox3/Source/Common/PsychHID/PsychHIDKbCheck.c		
  
  PROJECTS: PsychHID only.
  
  PLATFORMS:  All.
  
  AUTHORS:
  Allen.Ingling@nyu.edu             awi 
  mario.kleiner@tuebingen.mpg.de    mk

  HISTORY:
  5/12/03  awi		Created.
  12/17/09 rpw		Added keypad support
  07/28/11 mk           Refactored for multi-os support.
  
  TO DO:

*/

#include "PsychHID.h"

static char useString[]= "[keyIsDown,secs,keyCode]=PsychHID('KbCheck' [, deviceNumber][, scanList])";
static char synopsisString[] = 
        "Scan a keyboard, keypad, or other HID device with buttons, and return a vector of logical values indicating the "
        "state of each key.\n"
        "By default the first keyboard device (the one with the lowest device number) is "
        "scanned. If no keyboard is found, the first keypad device is scanned, followed by other "
        "devices, e.g., mice.  Optionally, the 'deviceNumber' of any keyboard or HID device may be specified.\n"
        "As checking all potentially 256 keys on a HID device is a time consuming process, "
        "which can easily take up to 1 msec on modern hardware, you can restrict the scan to a subset of "
        "the 256 keys by providing the optional 'scanList' parameter: 'scanList' must be a vector of 256 "
        "doubles, where the i'th element corresponds to the i'th key and a zero value means: Ignore this "
        "key during scan, whereas a positive non-zero value means: Scan this key.\n"
        "The PsychHID('KbCheck') implements the KbCheck command as provided by the  OS 9 Psychtoolbox. "
        "KbCheck is defined in the OS X Psychtoolbox and invokes PsychHID('KbCheck'). "
        "Always use KbCheck instead of directly calling PsychHID('KbCheck'), unless you have very good "
        "reasons to do otherwise and really know what you're doing!";
        
static char seeAlsoString[] = "";
 
PsychError PSYCHHIDKbCheck(void) 
{
    int					deviceIndex;
    int					m, n, p;
    double				*scanList = NULL;
    psych_bool                          isDeviceSpecified;

    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(2));
	
    // Choose the device index and its record
    isDeviceSpecified=PsychCopyInIntegerArg(1, FALSE, &deviceIndex);
    if (!isDeviceSpecified) {
        // set the keyboard or keypad device to be the first keyboard device or, if no keyboard, the first keypad
        deviceIndex = INT_MAX;
    }

    // Optional 2nd argument 'scanlist' provided?
    if (PsychAllocInDoubleMatArg(2, FALSE, &m, &n, &p, &scanList)) {
        // Yep. Matching size?
        if (p!=1 || m * n != 256) PsychErrorExitMsg(PsychError_user, "Provided 'scanList' parameter is not a vector of 256 doubles, as required!");
    }

    return(PsychHIDOSKbCheck(deviceIndex, scanList));
}

#if PSYCH_SYSTEM == PSYCH_OSX

#define NUMDEVICEUSAGES 7

PsychError PsychHIDOSKbCheck(int deviceIndex, double* scanList)
{
    pRecDevice			deviceRecord;
    pRecElement			currentElement, lastElement = NULL;
    int					i, debuglevel = 0;
    static int			numDeviceIndices = -1;
    int					numDeviceUsages=NUMDEVICEUSAGES;
	long				KbDeviceUsagePages[NUMDEVICEUSAGES]= {kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop};
	long				KbDeviceUsages[NUMDEVICEUSAGES]={kHIDUsage_GD_Keyboard, kHIDUsage_GD_Keypad, kHIDUsage_GD_Mouse, kHIDUsage_GD_Pointer, kHIDUsage_GD_Joystick, kHIDUsage_GD_GamePad, kHIDUsage_GD_MultiAxisController};
    static int			deviceIndices[PSYCH_HID_MAX_KEYBOARD_DEVICES]; 
    static pRecDevice	deviceRecords[PSYCH_HID_MAX_KEYBOARD_DEVICES];
    psych_bool			isDeviceSpecified, foundUserSpecifiedDevice;
    double				*timeValueOutput, *isKeyDownOutput, *keyArrayOutput;
    int					m, n, p, nout;
    double				dummyKeyDown;
    double				dummykeyArrayOutput[256];
    uint32_t            usage, usagePage;
    int                 value;

    // We query keyboard and keypad devices only on first invocation, then cache and recycle the data:
    if (numDeviceIndices == -1) {
        PsychHIDVerifyInit();
        PsychHIDGetDeviceListByUsages(numDeviceUsages, KbDeviceUsagePages, KbDeviceUsages, &numDeviceIndices, deviceIndices, deviceRecords);
    }
	
    // Choose the device index and its record
    isDeviceSpecified = (deviceIndex != INT_MAX);
    if(isDeviceSpecified){  //make sure that the device number provided by the user is really a keyboard or keypad.
        if (deviceIndex < 0) {
            debuglevel = 1;
            deviceIndex = -deviceIndex;
        }

        for(i=0;i<numDeviceIndices;i++){
            if(foundUserSpecifiedDevice=(deviceIndices[i]==deviceIndex))
                break;
        }
        if(!foundUserSpecifiedDevice)
            PsychErrorExitMsg(PsychError_user, "Specified device number is not a keyboard or keypad device.");
    } else { // set the keyboard or keypad device to be the first keyboard device or, if no keyboard, the first keypad
        i=0;
        if(numDeviceIndices==0)
            PsychErrorExitMsg(PsychError_user, "No keyboard or keypad devices detected.");
        else{
            deviceIndex=deviceIndices[i];
        }
    }
    deviceRecord=deviceRecords[i]; 

    // Allocate and init out return arguments.

    // Either alloc out the arguments, or redirect to
    // internal dummy variables. This to avoid mxMalloc() call overhead
    // inside the PsychAllocOutXXX() routines:
    nout = PsychGetNumNamedOutputArgs();

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
    for(currentElement=HIDGetFirstDeviceElement(deviceRecord, kHIDElementTypeInput | kHIDElementTypeCollection); 
        (currentElement != NULL) && (currentElement != lastElement);
        currentElement=HIDGetNextDeviceElement(currentElement, kHIDElementTypeInput | kHIDElementTypeCollection))
    {
        // Keep track of last queried element:
        lastElement = currentElement;
        
        #ifndef __LP64__
        usage = currentElement->usage;
        usagePage = currentElement->usagePage;
        #else
        usage = IOHIDElementGetUsage(currentElement);
        usagePage = IOHIDElementGetUsagePage(currentElement);
        // printf("PTB-DEBUG: [KbCheck]: ce %p page %d usage: %d isArray: %d\n", currentElement, usagePage, usage, IOHIDElementIsArray(currentElement));

        if (IOHIDElementGetType(currentElement) == kIOHIDElementTypeCollection) {
            CFArrayRef children = IOHIDElementGetChildren(currentElement);
            if (!children) continue;
            
            CFIndex idx, cnt = CFArrayGetCount(children);
            for ( idx = 0; idx < cnt; idx++ ) {
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
        #endif

        // Classic path, or 64-Bit path for non-collection elements:
        if(((usagePage == kHIDPage_KeyboardOrKeypad) || (usagePage == kHIDPage_Button)) && (usage <= 256) && (usage >= 1) &&
			( (scanList == NULL) || (scanList[usage - 1] > 0) ) ) {
            #ifndef __LP64__
            value = (int) HIDGetElementValue(deviceRecord, currentElement);
            #else
            value = (int) IOHIDElement_GetValue(currentElement, kIOHIDValueScaleTypePhysical);
            #endif

            if (debuglevel > 0) printf("PTB-DEBUG: [KbCheck]: usage: %x value: %d \n", usage, value);
            keyArrayOutput[usage - 1]=(value || (int) keyArrayOutput[usage - 1]);
            *isKeyDownOutput= keyArrayOutput[usage - 1] || *isKeyDownOutput; 
        }
    }

    return(PsychError_none);	
}

#endif
