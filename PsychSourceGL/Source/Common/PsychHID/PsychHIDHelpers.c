/*
  PsychToolbox/Source/Common/PsychHID/PsychHIDHelpers.c

  PROJECTS: PsychHID

  PLATFORMS:  All

  AUTHORS:
  Allen.Ingling@nyu.edu             awi
  mario.kleiner@tuebingen.mpg.de    mk

  HISTORY:
  5/05/03  awi        Created.
  4/19/05  dgp      cosmetic.
  8/23/07  rpw      added PsychHIDKbQueueRelease() to PsychHIDCleanup()
  4/04/09  mk        added support routines for generic USB devices and usbDeviceRecordBank.

  TO DO:

*/

#include "PsychHID.h"

#if PSYCH_SYSTEM != PSYCH_OSX

// Global variable which holds a reference to the last
// accessed hid_device* on non-OS/X. This is for error
// handling:
hid_device* last_hid_device = NULL;

// List with all known low level USB-HID devices from enumeration:
struct hid_device_info* hidlib_devices = NULL;
pRecDevice hid_devices = NULL;
hid_device* source[MAXDEVICEINDEXS];

#endif

#if PSYCH_SYSTEM == PSYCH_OSX
IOHIDDeviceInterface122** deviceInterfaces[MAXDEVICEINDEXS];
#endif

// Tracker used to maintain references to open generic USB devices.
// PsychUSBDeviceRecord is currently defined in PsychHID.h.
PsychUSBDeviceRecord usbDeviceRecordBank[PSYCH_HID_MAX_GENERIC_USB_DEVICES];

PsychHIDEventRecord* hidEventBuffer[PSYCH_HID_MAX_DEVICES];
unsigned int    hidEventBufferCapacity[PSYCH_HID_MAX_DEVICES];
unsigned int    hidEventBufferReadPos[PSYCH_HID_MAX_DEVICES];
unsigned int    hidEventBufferWritePos[PSYCH_HID_MAX_DEVICES];
psych_mutex     hidEventBufferMutex[PSYCH_HID_MAX_DEVICES];
psych_condition hidEventBufferCondition[PSYCH_HID_MAX_DEVICES];

/* PsychInitializePsychHID()
 *
 * Master init routine - Called at module load time / first time init.
 *
 */
void PsychInitializePsychHID(void)
{
    int i;

    // Initialize the generic USB tracker to "all off" state:
    for (i = 0; i < PSYCH_HID_MAX_GENERIC_USB_DEVICES; i++) {
        usbDeviceRecordBank[i].valid = 0;
    }

    // Setup event ringbuffers:
    for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
        hidEventBuffer[i] = NULL;
        hidEventBufferCapacity[i] = 10000; // Initial capacity of event buffer.
        hidEventBufferReadPos[i] = 0;
        hidEventBufferWritePos[i] = 0;
    }

#if PSYCH_SYSTEM == PSYCH_OSX
    for (i = 0; i < MAXDEVICEINDEXS; i++) deviceInterfaces[i] = NULL;

    // Try to load all bundles from Psychtoolbox/PsychHardware/
    // This loads the HID_Utilities.framework bundle if it is present. The whole point of it is
    // to allow our statically compiled-in version of the library to find the location of
    // the XML file with the database of (vendorId, productId) -> (VendorName, ProductName) and
    // (usagePage, usage) -> (usageName) mappings.
    //
    // In practice, the XML file only serves as a fallback, and one that doesn't contain much
    // useful info for mainstream products, only for a few niche products. Given its limited
    // value, i think we can refrain from shipping the framework as part of Psychtoolbox and
    // just provide the option to use it (== its XML file) if users decide to install it themselves.
    char tmpString[1024];

    sprintf(tmpString, "%sPsychHardware/", PsychRuntimeGetPsychtoolboxRoot(FALSE));
    CFStringRef urlString = CFStringCreateWithCString(kCFAllocatorDefault, tmpString, kCFStringEncodingASCII);
    CFURLRef directoryURL = CFURLCreateWithFileSystemPath(kCFAllocatorDefault, urlString, kCFURLPOSIXPathStyle, false);
    CFRelease(urlString);
    CFArrayRef bundleArray = CFBundleCreateBundlesFromDirectory(kCFAllocatorDefault, directoryURL, NULL);
    CFRelease(directoryURL);
    CFRelease(bundleArray);

#endif

    // Initialize OS specific interfaces and routines:
    PsychHIDInitializeHIDStandardInterfaces();

    // This sets up data structures for HID report reception inside PsychHIDReceiveReports.c:
    PsychHIDReleaseAllReportMemory();

    return;
}

/* PsychHIDGetFreeUSBDeviceSlot();
 *
 * Return a device record pointer to a free generic USB device
 * slot, as well as the associated numeric usbHandle.
 *
 * Abort with error if no more slots are free.
 */
PsychUSBDeviceRecord* PsychHIDGetFreeUSBDeviceSlot(int* usbHandle)
{
    int i;

    // Find the next available USB slot:
    for (i = 0; i < PSYCH_HID_MAX_GENERIC_USB_DEVICES; i++) {
        if (usbDeviceRecordBank[i].valid == 0) {
            *usbHandle = i;
            return(&(usbDeviceRecordBank[i]));
        }
    }

    // If we reach this point, then all slots are occupied: Fail!
    PsychErrorExitMsg(PsychError_user, "Unable to open another generic USB device! Too many devices open. Please close one and retry.");
    return(NULL);
}

/* PsychHIDGetFreeUSBDeviceSlot();
 *
 * Return a device record pointer to a free generic USB device
 * slot, as well as the associated numeric usbHandle.
 *
 * Abort with error if no more slots are free.
 */
PsychUSBDeviceRecord* PsychHIDGetUSBDevice(int usbHandle)
{
    // Child protection:
    if (usbHandle < 0 || usbHandle >= PSYCH_HID_MAX_GENERIC_USB_DEVICES) PsychErrorExitMsg(PsychError_user, "Invalid generic USB device handle provided! Handle outside valid range.");
    if (usbDeviceRecordBank[usbHandle].valid == 0) PsychErrorExitMsg(PsychError_user, "Invalid generic USB device handle provided! The handle doesn't correspond to an open device.");

    // Valid handle for slot corresponding to an open device. Return PsychUSBDeviceRecord* to it:
    return(&(usbDeviceRecordBank[usbHandle]));
}

void PsychHIDCloseAllUSBDevices(void)
{
    int i;
    for (i = 0; i < PSYCH_HID_MAX_GENERIC_USB_DEVICES; i++) {
        if (usbDeviceRecordBank[i].valid) {
            PsychHIDOSCloseUSBDevice(PsychHIDGetUSBDevice(i));
        }
    }
}

/*
    PsychHIDCleanup()

    Cleanup before flushing the mex file.
*/
PsychError PsychHIDCleanup(void)
{
    pRecDevice curdev = NULL;

    // Disable online help system:
    PsychClearGiveHelp();

    // Disable any kind of low-level stdin<->tty magic for character reception
    // or suppression in console mode (for octave and matlab -nojvm):
    ConsoleInputHelper(-10);

    // Shutdown USB-HID report low-level functions, e.g., for DAQ toolbox on OS/X:
    PsychHIDReceiveReportsCleanup(); // PsychHIDReceiveReport.c

    // Shutdown os specific interfaces and routines:
    PsychHIDShutdownHIDStandardInterfaces();

    // Release all other HID device data structures:
#if PSYCH_SYSTEM == PSYCH_OSX
    // Via Apple HIDUtils:
    int i;
    for (i = 0; i < MAXDEVICEINDEXS; i++) {
        if (deviceInterfaces[i]) {
            IOHIDDeviceInterface** interface = (IOHIDDeviceInterface**)deviceInterfaces[i];
            (*interface)->close(interface);
            (*interface)->Release(interface);
            deviceInterfaces[i] = NULL;
        }
    }

    if (HIDHaveDeviceList()) HIDReleaseDeviceList();
#else
    // Then our high-level list:
    while (hid_devices) {
        // Get current device record to release:
        curdev = hid_devices;

        // Advance to the next one for next loop iteration:
        hid_devices = hid_devices->pNext;

        // Interface attached aka device opened? If so we need to close the device handle:
        if (curdev->interface) hid_close((hid_device*)curdev->interface);

        // Release:
        free(curdev);
    }

    // Reset last hid device for error handling:
    last_hid_device = NULL;

    // Release the HIDLIB low-level device list:
    if (hidlib_devices) hid_free_enumeration(hidlib_devices);
    hidlib_devices = NULL;

    // Shutdown HIDAPI:
    hid_exit();

#endif

    // Close and release all open generic USB devices:
    PsychHIDCloseAllUSBDevices();

    return(PsychError_none);
}

/*
    PsychHIDGetDeviceRecordPtrFromIndex()

    The inverse of PsychHIDGetIndexFromRecord()

    Accept the index from the list of device records and return a pointer to the indicated record.  Externally the list is one-indexed.
*/
pRecDevice PsychHIDGetDeviceRecordPtrFromIndex(int deviceIndex)
{
    int        i;
    pRecDevice currentDevice = NULL;

    PsychHIDVerifyInit();
    i = 1;
    for (currentDevice = HIDGetFirstDevice(); currentDevice != NULL; currentDevice = HIDGetNextDevice(currentDevice)) {
        if (i==deviceIndex) {
#if PSYCH_SYSTEM != PSYCH_OSX
            if (!currentDevice->interface) {
                currentDevice->interface = (void*)hid_open_path(currentDevice->transport);
                if (!currentDevice->interface) PsychErrorExitMsg(PsychError_system, "HIDLIB Failed to open USB device!");

                // Set read ops on device to non-blocking:
                hid_set_nonblocking((hid_device*)currentDevice->interface, 1);
            }
#endif
            return currentDevice;
        }
        ++i;
    }

    PsychErrorExitMsg(PsychError_user, "Invalid device index specified. Has a device been unplugged? Try rebuilding the device list");
    return(NULL);  //make the compiler happy.
}

psych_bool PsychHIDCreateEventBuffer(int deviceIndex, int numValuators, int numSlots)
{
    unsigned int bufferSize;

    if (deviceIndex < 0) deviceIndex = PsychHIDGetDefaultKbQueueDevice();

    if (numSlots < 0) {
        printf("PTB-ERROR: PsychHIDCreateEventBuffer(): numSlots %i invalid. Must be at least 0.\n", numSlots);
        return(FALSE);
    }

    // Non-zero numSlots ==> Set new capacity, otherwise leave at default/last capacity:
    if (numSlots > 0)
        hidEventBufferCapacity[deviceIndex] = numSlots;

    bufferSize = hidEventBufferCapacity[deviceIndex];

    // Already created? If so, nothing to do:
    if (hidEventBuffer[deviceIndex] || (bufferSize < 1)) return(FALSE);

    if (numValuators > PSYCH_HID_MAX_VALUATORS) {
        printf("PTB-ERROR: PsychHIDCreateEventBuffer(): numValuators %i > current compile time maximum of %i!\n",
               numValuators, PSYCH_HID_MAX_VALUATORS);
        return(FALSE);
    }

    hidEventBuffer[deviceIndex] = (PsychHIDEventRecord*) calloc(sizeof(PsychHIDEventRecord), bufferSize);
    if (NULL==hidEventBuffer[deviceIndex]) {
        printf("PTB-ERROR: PsychHIDCreateEventBuffer(): Insufficient memory to create KbQueue event buffer!");
        return(FALSE);
    }

    // Prepare mutex for buffer:
    PsychInitMutex(&hidEventBufferMutex[deviceIndex]);
    PsychInitCondition(&hidEventBufferCondition[deviceIndex], NULL);

    // Init & Flush it:
    hidEventBufferWritePos[deviceIndex] = 0;
    PsychHIDFlushEventBuffer(deviceIndex);

    return(TRUE);
}

psych_bool PsychHIDDeleteEventBuffer(int deviceIndex)
{
    if (deviceIndex < 0) deviceIndex = PsychHIDGetDefaultKbQueueDevice();

    if (hidEventBuffer[deviceIndex]) {
        // Empty the buffer, reset read/writepointers:
        PsychHIDFlushEventBuffer(deviceIndex);

        // Release it:
        free(hidEventBuffer[deviceIndex]);
        hidEventBuffer[deviceIndex] = NULL;
        PsychDestroyMutex(&hidEventBufferMutex[deviceIndex]);
        PsychDestroyCondition(&hidEventBufferCondition[deviceIndex]);
    }

    return TRUE;
}

psych_bool PsychHIDFlushEventBuffer(int deviceIndex)
{
    if (deviceIndex < 0) deviceIndex = PsychHIDGetDefaultKbQueueDevice();

    if (!hidEventBuffer[deviceIndex]) return FALSE;

    PsychLockMutex(&hidEventBufferMutex[deviceIndex]);
    hidEventBufferReadPos[deviceIndex] = hidEventBufferWritePos[deviceIndex];
    PsychUnlockMutex(&hidEventBufferMutex[deviceIndex]);

    return TRUE;
}

/* Return number of events in buffer for 'deviceIndex':
 * flags == 0 -> All events.
 * flags &  1 -> Only keypress events with valid mapped ASCII CookedKey keycode.
 */
unsigned int PsychHIDAvailEventBuffer(int deviceIndex, unsigned int flags)
{
    unsigned int navail, i, j;

    if (deviceIndex < 0) deviceIndex = PsychHIDGetDefaultKbQueueDevice();

    if (!hidEventBuffer[deviceIndex]) return(0);

    PsychLockMutex(&hidEventBufferMutex[deviceIndex]);

    // Compute total number of available events by default:
    navail = hidEventBufferWritePos[deviceIndex] - hidEventBufferReadPos[deviceIndex];

    // Only count of valid "CookedKey" mapped keypress events, e.g., for use by CharAvail(), requested?
    if (flags & 1) {
        // Yes: Iterate over all available events and only count number of keypress events
        // with meaningful 'CookedKey' field:
        navail = 0;
        for (i = hidEventBufferReadPos[deviceIndex]; i < hidEventBufferWritePos[deviceIndex]; i++) {
            j = i % hidEventBufferCapacity[deviceIndex];
            if ((hidEventBuffer[deviceIndex][j].status & (1 << 0)) && (hidEventBuffer[deviceIndex][j].cookedEventCode > 0)) navail++;
        }
    }

    PsychUnlockMutex(&hidEventBufferMutex[deviceIndex]);

    return(navail);
}

int PsychHIDReturnEventFromEventBuffer(int deviceIndex, int outArgIndex, double maxWaitTimeSecs)
{
    unsigned int navail, j;
    PsychHIDEventRecord evt;
    PsychGenericScriptType *retevent;
    double* foo = NULL;
    PsychGenericScriptType *outMat;
    double *v;
    const char *FieldNames[] = { "Type", "Time", "Pressed", "Keycode", "CookedKey", "ButtonStates", "Motion", "X", "Y", "NormX", "NormY", "Valuators" };

    if (deviceIndex < 0) deviceIndex = PsychHIDGetDefaultKbQueueDevice();
    if (!hidEventBuffer[deviceIndex]) return(0);

    PsychLockMutex(&hidEventBufferMutex[deviceIndex]);
    navail = hidEventBufferWritePos[deviceIndex] - hidEventBufferReadPos[deviceIndex];

    // If nothing available and we're asked to wait for something, then wait:
    if ((navail == 0) && (maxWaitTimeSecs > 0)) {
        // Wait for something:
        PsychTimedWaitCondition(&hidEventBufferCondition[deviceIndex], &hidEventBufferMutex[deviceIndex], maxWaitTimeSecs);

        // Recompute number of available events:
        navail = hidEventBufferWritePos[deviceIndex] - hidEventBufferReadPos[deviceIndex];
    }

    // Check if anything available, copy it if so:
    if (navail) {
        memcpy(&evt, &(hidEventBuffer[deviceIndex][hidEventBufferReadPos[deviceIndex] % hidEventBufferCapacity[deviceIndex]]), sizeof(PsychHIDEventRecord));
        hidEventBufferReadPos[deviceIndex]++;
    }

    PsychUnlockMutex(&hidEventBufferMutex[deviceIndex]);

    if (navail) {
        // Return event struct:
        switch (evt.type) {
            case 0: // Press/Release
            case 1: // Motion/Valuator change
                PsychAllocOutStructArray(outArgIndex, kPsychArgOptional, -1, 12, FieldNames, &retevent);
                break;

            case 2: // Touch begin
            case 3: // Touch update/move
            case 4: // Touch end
            case 5: // Touch sequence compromised marker. If this one shows up - with magic touch point
                    // id 0xffffffff btw., then the user script knows the sequence was cut short / aborted
                    // by some higher priority consumer, e.g., some global gesture recognizer.
                PsychAllocOutStructArray(outArgIndex, kPsychArgOptional, -1, 12, FieldNames, &retevent);
                break;

            default:
                PsychErrorExitMsg(PsychError_internal, "Unhandled keyboard queue event type!");
        }

        PsychSetStructArrayDoubleElement("Type",         0, (double) evt.type,                        retevent);
        PsychSetStructArrayDoubleElement("Time",         0, evt.timestamp,                            retevent);
        PsychSetStructArrayDoubleElement("Pressed",      0, (double) (evt.status & (1 << 0)) ? 1 : 0, retevent);
        PsychSetStructArrayDoubleElement("Keycode",      0, (double) evt.rawEventCode,                retevent);
        PsychSetStructArrayDoubleElement("CookedKey",    0, (double) evt.cookedEventCode,             retevent);
        PsychSetStructArrayDoubleElement("ButtonStates", 0, (double) evt.buttonStates,                retevent);
        PsychSetStructArrayDoubleElement("Motion",       0, (double) (evt.status & (1 << 1)) ? 1 : 0, retevent);
        PsychSetStructArrayDoubleElement("X",            0, (double) evt.X,                           retevent);
        PsychSetStructArrayDoubleElement("Y",            0, (double) evt.Y,                           retevent);
        PsychSetStructArrayDoubleElement("NormX",        0, (double) evt.normX,                       retevent);
        PsychSetStructArrayDoubleElement("NormY",        0, (double) evt.normY,                       retevent);

        // Copy out all valuators (including redundant (X,Y) again:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, evt.numValuators, 1, &v, &outMat);
        for (j = 0; j < evt.numValuators; j++)
            *(v++) = (double) evt.valuators[j];
        PsychSetStructArrayNativeElement("Valuators", 0, outMat, retevent);

        return(navail - 1);
    }
    else {
        // Return empty matrix:
        PsychCopyOutDoubleMatArg(outArgIndex, kPsychArgOptional, 0, 0, 0, foo);
        return(0);
    }
}

PsychHIDEventRecord* PsychHIDLastTouchEventFromEventBuffer(int deviceIndex, int touchID)
{
    int nend, current;
    PsychHIDEventRecord *evt;

    if (!hidEventBuffer[deviceIndex]) return(0);

    PsychLockMutex(&hidEventBufferMutex[deviceIndex]);
    nend = (hidEventBufferWritePos[deviceIndex] - 1) % hidEventBufferCapacity[deviceIndex];
    current = nend;

    // Go backwards through all touch events until the most recent one with touchID is found:
    do {
        if ((hidEventBuffer[deviceIndex][current].type >= 2) &&
            (hidEventBuffer[deviceIndex][current].type <= 4) &&
            (hidEventBuffer[deviceIndex][current].rawEventCode == touchID))
            break;

        current = (current - 1) % hidEventBufferCapacity[deviceIndex];
    } while ((current != nend) && (current >= 0));

    if (hidEventBuffer[deviceIndex][current].rawEventCode == touchID)
        evt = &(hidEventBuffer[deviceIndex][current]);
    else
        evt = NULL;

    PsychUnlockMutex(&hidEventBufferMutex[deviceIndex]);

    return (evt);
}

int PsychHIDAddEventToEventBuffer(int deviceIndex, PsychHIDEventRecord* evt)
{
    unsigned int navail;

    if (deviceIndex < 0) deviceIndex = PsychHIDGetDefaultKbQueueDevice();

    if (!hidEventBuffer[deviceIndex]) return 0;

    PsychLockMutex(&hidEventBufferMutex[deviceIndex]);

    navail = hidEventBufferWritePos[deviceIndex] - hidEventBufferReadPos[deviceIndex];
    if (navail < hidEventBufferCapacity[deviceIndex]) {
        memcpy(&(hidEventBuffer[deviceIndex][hidEventBufferWritePos[deviceIndex] % hidEventBufferCapacity[deviceIndex]]), evt, sizeof(PsychHIDEventRecord));
        hidEventBufferWritePos[deviceIndex]++;

        // Announce new event to potential waiters:
        PsychSignalCondition(&hidEventBufferCondition[deviceIndex]);
    }
    else {
        printf("PsychHID: WARNING: KbQueue event buffer is full! Maximum capacity of %i elements reached, will discard future events.\n", hidEventBufferCapacity[deviceIndex]);
    }

    PsychUnlockMutex(&hidEventBufferMutex[deviceIndex]);

    return navail - 1;
}

// Platform specific code starts here:
// ===================================

#if PSYCH_SYSTEM == PSYCH_OSX

/*
    PSYCHHIDCheckInit()

    Check to see if we need to create the USB-HID device list. If it has not been created then create it.
*/
void PsychHIDVerifyInit(void)
{
    // GenericDesktop for mice, keyboards, gamepads etc. vendor defined for things like MCC USB-DAQ boxes...
    UInt32 inUsagePages[2] = {kHIDPage_GenericDesktop, kHIDPage_VendorDefinedStart};
    UInt32 inUsages[2] = {0, 0};
    int inNumDeviceTypes = 2;
    psych_bool success = TRUE;

    // Build HID device list if it doesn't already exist:
    if (!HIDHaveDeviceList()) success = (psych_bool)HIDBuildDeviceList(0, 0);

    // This check can only be made against the 64-Bit HID Utilities, as the older 32-Bit
    // version is even more crappy and can't report meaningful error status:
    if (!success) {
        printf("PsychHID-ERROR: Could not enumerate and attach to all HID devices (HIDBuildDeviceList(0,0) failed)!\n");
        printf("PsychHID-ERROR: One reason could be that some HID devices are already exclusively claimed by some 3rd party device drivers\n");
        printf("PsychHID-ERROR: or applications. I will now retry to only claim control of a hopefully safe subset of devices like standard\n");
        printf("PsychHID-ERROR: keyboards, mice, gamepads and supported USB-DAQ devices and other vendor defined devices and hope this goes better...\n");

        // User override for special devices provided?
        if (getenv("PSYCHHID_HIDPAGEOVERRIDE")) {
            // Override via environment variable provided: Use it in addition to the generic desktop input devices:
            inUsagePages[1] = atoi(getenv("PSYCHHID_HIDPAGEOVERRIDE"));
            inUsages[1] = atoi(getenv("PSYCHHID_HIDUSAGEOVERRIDE"));
            printf("PsychHID-INFO: User override: Generic desktop page devices + usagePage 0x%x and usage value 0x%x ...\n", inUsagePages[1], inUsages[1]);
        }

        success = (psych_bool)HIDBuildMultiDeviceList(inUsagePages, inUsages, inNumDeviceTypes);
        if (!success) {
            printf("PsychHID-ERROR: Nope. Excluding everything but a few known USB-DAQ devices and standard mice, keyboards etc. Retrying ...\n");
            inUsagePages[1] = kHIDPage_VendorDefinedStart;
            inUsages[1] = 1; // Known combo is usagepage 0xff00 + usage 0x1 for MCC USB 1208-FS USB HID-DAQ device.
            success = (psych_bool)HIDBuildMultiDeviceList(inUsagePages, inUsages, inNumDeviceTypes);
        }

        if (!success) {
            printf("PsychHID-ERROR: Nope. Excluding everything but standard mice, keyboards etc. Retrying ...\n");
            success = (psych_bool)HIDBuildMultiDeviceList(inUsagePages, inUsages, 1);
            if (!success) {
                printf("PsychHID-ERROR: Could not enumerate any HID devices (HIDBuildDeviceList() still failed even for standard generic desktop input devices)!\n");
                printf("PsychHID-ERROR: Reasons can be ranging from bugs in Apples HID software to a buggy HID device driver for some connected device,\n");
                printf("PsychHID-ERROR: to general operating system malfunction. A reboot or device driver update for 3rd party HID devices\n");
                printf("PsychHID-ERROR: maybe could help. Check the OSX system log for possible HID related error messages or hints. Aborting...\n");
                PsychErrorExitMsg(PsychError_system, "HID device enumeration failed due to malfunction in the OSX 64 Bit Apple HID Utilities framework.");
            }
        }
        else {
            printf("PsychHID-INFO: That worked. A subset of regular mouse, keyboard etc. input devices and maybe some vendor defined devices will be available at least.\n");
        }
    }

    // Double-Check to protect against pathetic Apple software:
    if (!HIDHaveDeviceList()) {
        printf("PsychHID-ERROR: Could not enumerate HID devices (HIDBuildDeviceList() success, but HIDHaveDeviceList() still failed)!\n");
        printf("PsychHID-ERROR: Reasons can be ranging from bugs in Apples HID software to a buggy HID device driver for some connected device,\n");
        printf("PsychHID-ERROR: to general operating system malfunction. A reboot or device driver update for 3rd party HID devices\n");
        printf("PsychHID-ERROR: maybe could help. Check the OSX system log for possible HID related error messages or hints. Aborting...\n");
        PsychErrorExitMsg(PsychError_system, "HID device enumeration failed due to malfunction in the OSX Apple HID Utilities framework (II).");
    }

    // Verify no security sensitive application is blocking our low-level access to HID devices:
    PsychHIDWarnInputDisabled(NULL);
}

/*
    PsychHIDWarnInputDisabled()

    Check if HID event input is disabled by external processes, e.g., due to
    secure password entry protection. Return TRUE and output a warning message
    to user if HID input won't work due to active security measures. Return FALSE
    and stay silent if PsychHID can work as expected.

*/
psych_bool PsychHIDWarnInputDisabled(const char* callerName)
{
    if (IsSecureEventInputEnabled()) {
        printf("PTB-WARNING: During %s: Some other running application is preventing me from accessing the keyboard/keypad/mouse/...!\n", (callerName) ? callerName : "PsychHID invocation");
        printf("PTB-WARNING: This is likely a security measure, e.g., to protect some active password entry field.\n");
        printf("PTB-WARNING: Please identify and quit the offending application. E.g., some versions of Firefox are known to cause such problems...\n");
        return TRUE;
    }

    return FALSE;
}

// Device interface setup & query: 10.5 Leopard and later:
//
// ---------------------------------
// This routines is largely transplanted - with modifications - from HID Utilities v1.0:
// Create and open an interface to device, required prior to extracting values or building queues
// Note: application now owns the device and must close and release it prior to exiting

//unsigned long HIDCreateOpenDeviceInterface (UInt32 hidDevice, pRecDevice pDevice)
static IOHIDDeviceInterface122** HIDCreateOpenDeviceInterface(pRecDevice pDevice)
{
    IOReturn result = kIOReturnSuccess;
    HRESULT plugInResult = S_OK;
    SInt32 score = 0;
    IOCFPlugInInterface** ppPlugInInterface = NULL;
    IOHIDDeviceInterface** interface = NULL;
    io_service_t hidDevice = 0L;
    mach_port_t port;

    // Get low-level device for given HIDDeviceRef from HID Utilities v2.0:
    hidDevice = AllocateHIDObjectFromIOHIDDeviceRef((IOHIDDeviceRef)pDevice);
    if (0L == hidDevice) {
        printf("PTB-ERROR: PsychHID: Failed to allocate low-level HID-Object for high-level HID device.\n");
        return(NULL);
    }

    // This code is borrowed from HID Utilities v1.0:
    // Create and open interface for IORegistry device:
    result = IOCreatePlugInInterfaceForService(hidDevice, kIOHIDDeviceUserClientTypeID,
                                               kIOCFPlugInInterfaceID, &ppPlugInInterface, &score);
    if ((kIOReturnSuccess == result) && ppPlugInInterface)
    {
        // Call a method of the intermediate plug-in to create the device interface
        plugInResult = (*ppPlugInInterface)->QueryInterface(ppPlugInInterface,
                                                            CFUUIDGetUUIDBytes(kIOHIDDeviceInterfaceID), (void *) &interface);
        if ((S_OK != plugInResult) || (NULL == interface))
            printf("PTB-ERROR: PsychHID: Couldn't query HID class device interface from plugInInterface: %x.\n", plugInResult);

        // Release ppPlugInInterface in any case - no longer needed:
        IODestroyPlugInInterface(ppPlugInInterface);
    }
    else
        printf("PTB-ERROR: PsychHID: Failed to create **plugInInterface via IOCreatePlugInInterfaceForService: %x.\n", result);

    // Successfully retrieved an interface?
    if (NULL != interface)
    {
        // Yes: Open it.
        result = (*interface)->open(interface, 0);
        if (kIOReturnSuccess != result) {
            printf("PTB-ERROR: PsychHID: Failed to open HID device low-level interface via open: %x.\n", result);
        }
        else {
            // Interface open. Create async mach port for it to deliver interrupt in endpoint data:
            result = (*interface)->createAsyncPort(interface, &port);
            if (kIOReturnSuccess != result) {
                printf("PTB-ERROR: PsychHID: Failed to create async port for HID device low-level interface: %x.\n", result);
            }
        }
    }

    // Release reference to high-level device via HID Utilities v2.0::
    if (hidDevice) FreeHIDObject(hidDevice);

    return((IOHIDDeviceInterface122**)interface);
}

IOHIDDeviceInterface122** PsychHIDGetDeviceInterfacePtrFromIndex(int deviceIndex)
{
    IOHIDDeviceInterface122 **interface = NULL;
    pRecDevice dev = PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
    if (dev) {
        // Interface already open?
        if (NULL == deviceInterfaces[deviceIndex]) {
            // No. Need to create it:
            deviceInterfaces[deviceIndex] = HIDCreateOpenDeviceInterface(dev);
        }

        // Assign it:
        interface = deviceInterfaces[deviceIndex];
    }

    return interface;
}

/*
    PsychHIDGetDeviceListByUsage()
*/
void PsychHIDGetDeviceListByUsage(long usagePage, long usage, int *numDeviceIndices, int *deviceIndices, pRecDevice *deviceRecords)
{
    pRecDevice currentDevice;
    int        currentDeviceIndex;

    PsychHIDVerifyInit();
    currentDeviceIndex = 0;
    *numDeviceIndices = 0;
    for (currentDevice = HIDGetFirstDevice(); currentDevice != NULL; currentDevice = HIDGetNextDevice(currentDevice)) {
        ++currentDeviceIndex;
        if (IOHIDDevice_GetUsagePage(currentDevice) == usagePage && IOHIDDevice_GetUsage(currentDevice) == usage) {
            deviceRecords[*numDeviceIndices] = currentDevice;
            deviceIndices[*numDeviceIndices] = currentDeviceIndex;  //the array is 0-indexed, devices are 1-indexed.
            ++(*numDeviceIndices);
        }
    }
}

/*
    PsychHIDGetDeviceListByUsages()
 */
void PsychHIDGetDeviceListByUsages(int numUsages, long *usagePages, long *usages, int *numDeviceIndices, int *deviceIndices, pRecDevice *deviceRecords)
{
    pRecDevice currentDevice;
    int        currentDeviceIndex;
    int        currentUsage;
    long       *usagePage;
    long       *usage;

    PsychHIDVerifyInit();
    *numDeviceIndices = 0;
    for (usagePage = usagePages, usage = usages, currentUsage = 0; currentUsage < numUsages; usagePage++, usage++, currentUsage++) {
        currentDeviceIndex = 0;
        for (currentDevice = HIDGetFirstDevice(); currentDevice != NULL; currentDevice = HIDGetNextDevice(currentDevice)) {
            ++currentDeviceIndex;
            if (IOHIDDevice_GetPrimaryUsagePage(currentDevice) == *usagePage && IOHIDDevice_GetPrimaryUsage(currentDevice) == *usage) {
                deviceRecords[*numDeviceIndices] = currentDevice;
                deviceIndices[*numDeviceIndices] = currentDeviceIndex;  //the array is 0-indexed, devices are 1-indexed.
                ++(*numDeviceIndices);
            }
        }
    }
}


/*
    PsychHIDGetIndexFromRecord()

    The inverse of PsychHIDGetDeviceRecordPtrFromIndex.

    This O(n) where n is the number of device elements.   We could make it O(1) if we modified
    the element structure in the HID Utilities library to include a field specifying the index of the element or
    collection.

    Note that if PsychHIDGetIndexFromRecord() is O(n) then its caller, PsychHIDGetCollections, is O(n^2) for each
    device, whereas if PsychHIDGetIndexFromRecord() is O(1) then psychHIDGetCollections becomes O(n) for each
    device.
*/
int PsychHIDGetIndexFromRecord(pRecDevice deviceRecord, pRecElement elementRecord, HIDElementTypeMask typeMask)
{
    int         elementIndex;
    pRecElement currentElement;

    if (elementRecord == NULL)
        return(0);
    elementIndex = 1;
    for (currentElement = HIDGetFirstDeviceElement(deviceRecord, typeMask);
        currentElement != elementRecord && currentElement != NULL;
        currentElement = HIDGetNextDeviceElement(currentElement, typeMask))
        ++elementIndex;
    if (currentElement == elementRecord)
        return elementIndex;
    else {
        PsychErrorExitMsg(PsychError_internal, "Element record not found within device record");
        return 0; //make the compiler happy
    }
}


pRecElement PsychHIDGetElementRecordFromDeviceRecordAndElementIndex(pRecDevice deviceRecord, int elementIndex)
{
    int         i;
    pRecElement currentElement;

    PsychHIDVerifyInit();
    i = 1;
    for (currentElement = HIDGetFirstDeviceElement(deviceRecord, kHIDElementTypeIO);
        currentElement != NULL;
        currentElement = HIDGetNextDeviceElement(currentElement, kHIDElementTypeIO))
    {
        if (i == elementIndex)
            return currentElement;
        ++i;
    }
    PsychErrorExitMsg(PsychError_internal, "Invalid device index specified.  Has a device has been unplugged? Try rebuilding the device list");
    return NULL;  //make the compiler happy.

}


pRecElement PsychHIDGetCollectionRecordFromDeviceRecordAndCollectionIndex(pRecDevice deviceRecord, int elementIndex)
{
    int         i;
    pRecElement currentElement;

    PsychHIDVerifyInit();
    i = 1;
    for (currentElement = HIDGetFirstDeviceElement(deviceRecord, kHIDElementTypeCollection);
        currentElement != NULL;
        currentElement = HIDGetNextDeviceElement(currentElement, kHIDElementTypeCollection))
    {
        if (i == elementIndex)
            return currentElement;
        ++i;
    }
    PsychErrorExitMsg(PsychError_internal, "Invalid collection index specified.  Has a device has been unplugged? Try rebuilding the device list");
    return NULL;  //make the compiler happy.

}


/*
        PsychHIDQueryOpenDeviceInterfaceFromDeviceIndex()

        Check the interface field of the libHIDUtilities device structure for NULLness.  libHIDUtilities.h seems to indicate that it is neccessary for application
        to invoke HIDCreateOpenDeviceInterface() before accessing a device.  However,
        1) libHIDUtilities provides no way to obtain a value for the required first argument to HIDCreateOpenDeviceInterface().
        2) Apple's example HID Explorer application does not call HIDCreateOpenDeviceInterface().
        3) Internally, libHIDUtilities invokes HIDCreateOpenDeviceInterface() itself when HIDBuildDeviceList() is called.

        Because the call lies within mysterious conditionals there is some uncertainty about whether HIDCreateOpenDeviceInterface() will always
        invoke HIDBuildDeviceList().  Therefore, PsychHID verifies that the device interface has been opened before accessing the elements of a device.

*/
psych_bool PsychHIDQueryOpenDeviceInterfaceFromDeviceIndex(int deviceIndex)
{
    PsychHIDVerifyInit();
    return PsychHIDGetDeviceInterfacePtrFromIndex(deviceIndex) != NULL;
}

void PsychHIDVerifyOpenDeviceInterfaceFromDeviceIndex(int deviceIndex)
{
    if (!PsychHIDQueryOpenDeviceInterfaceFromDeviceIndex(deviceIndex))
        PsychErrorExitMsg(PsychError_internal, "Device interface field is NULL.  libHIDUtilities failed to open the device interface ?");
}

/*
    PsychHIDGetTypeMaskStringFromTypeMask()

    Apple's HID Utilities library uses two different specificationos of the device type:
        - enumerated type HIDElementTypeMask, used to specify which element types to return.
        - unsigned long type identifier constants.

    The mask values will not mask the unsigned long type identifiers.
*/
void PsychHIDGetTypeMaskStringFromTypeMask(HIDElementTypeMask maskValue, char **pStr)
{
    char *maskNames[] = {"input", "output", "feature", "collection", "io", "all",};

    switch (maskValue)
    {
    case kHIDElementTypeInput:
        *pStr = maskNames[0];
        break;
    case kHIDElementTypeOutput:
        *pStr = maskNames[1];
        break;
    case kHIDElementTypeFeature:
        *pStr = maskNames[2];
        break;
    case kHIDElementTypeCollection:
        *pStr = maskNames[3];
        break;
    case kHIDElementTypeIO:
        *pStr = maskNames[4];
        break;
    case kHIDElementTypeAll:
        *pStr = maskNames[5];
        break;
    }
}


/*
    PsychHIDCountCollectionElements()

    Non-recursively count all elements of a collection which are of the specified type.

    HID element records hold three pointers to other element records: pPrevious, pChild and pSibling.  PsychHIDCountCollectionElements()
    operates on the theory that the members of a collection are its child and all of that child's siblings.


*/
int PsychHIDCountCollectionElements(pRecElement collectionRecord, HIDElementTypeMask elementTypeMask)
{
    pRecElement         currentElement;
    int                 numElements = 0;
    CFIndex             i, nmax;
    HIDElementTypeMask  currentElementMaskValue;

    CFArrayRef children = IOHIDElementGetChildren(collectionRecord);
    nmax = CFArrayGetCount(children);
    for (i = 0; i < nmax; i++) {
        currentElement = (pRecElement)CFArrayGetValueAtIndex(children, i);
        currentElementMaskValue = HIDConvertElementTypeToMask(IOHIDElementGetType(currentElement));
        if (currentElementMaskValue & elementTypeMask) ++numElements;
    }
    return numElements;
}


/*
    FindCollectionElements()

    Non-recursively return of a list of a collection's memember elements.

    HID element records hold three pointers to other element records: pPrevious, pChild and pSibling.  FindCollectionElements()
    operates on the theory that the members of a collection are its child and all of that child's siblings.

*/
int PsychHIDFindCollectionElements(pRecElement collectionRecord, HIDElementTypeMask elementTypeMask, pRecElement *collectionMembers, int maxListElements)
{
    pRecElement         currentElement;
    int                 numElements = 0;
    CFIndex             i, nmax;
    HIDElementTypeMask  currentElementMaskValue;

    CFArrayRef children = IOHIDElementGetChildren(collectionRecord);
    nmax = CFArrayGetCount(children);
    for (i = 0; i < nmax; i++) {
        currentElement = (pRecElement)CFArrayGetValueAtIndex(children, i);
        currentElementMaskValue = HIDConvertElementTypeToMask(IOHIDElementGetType(currentElement));
        if (currentElementMaskValue & elementTypeMask) {
            if (numElements == maxListElements) PsychErrorExitMsg(PsychError_internal, "Number of collection elements exceeds allocated storage space.");
            collectionMembers[numElements] = currentElement;
            ++numElements;
        }
    }

    return numElements;
}

#else

/* Linux and MS-Windows support via HIDLIB: */
/* ======================================== */

/*
    PSYCHHIDCheckInit()

    Check to see if we need to create the USB-HID device list. If it has not been created then create it.
*/
void PsychHIDVerifyInit(void)
{
    int busId, devId;
    pRecDevice currentDevice = NULL;
    struct hid_device_info* hid_dev = NULL;

    // If hid_devices list of all HID devices not yet initialized,
    // perform device enumeration:
    if (!hidlib_devices) {
        // Low-Level enumeration by HIDLIB:
        hidlib_devices = hid_enumerate(0x0, 0x0);

        // Build our own higher-level device list filled with info
        // from the low-level list:
        for (hid_dev = hidlib_devices; hid_dev != NULL; hid_dev = hid_dev->next) {
            // Allocate and zero-init high level struct currentDevice:
            currentDevice = calloc(1, sizeof(recDevice));

            // Copy low-level props to corresponding high-level props:
            currentDevice->usagePage = hid_dev->usage_page;
            currentDevice->usage = hid_dev->usage;
            // Abuse the "transport" string for the device path. On OS/X this just contains "USB":
            sprintf(&currentDevice->transport[0], "%s", hid_dev->path);
            currentDevice->vendorID = hid_dev->vendor_id;
            currentDevice->productID = hid_dev->product_id;
            currentDevice->version = hid_dev->release_number;
            if (hid_dev->manufacturer_string) wcstombs(&currentDevice->manufacturer[0], hid_dev->manufacturer_string, 256);
            if (hid_dev->product_string) wcstombs(&currentDevice->product[0], hid_dev->product_string, 256);
            if (hid_dev->serial_number) wcstombs(&currentDevice->serial[0], hid_dev->serial_number, 256);

            // Convert unique device path into unique numeric location id:
            if (PSYCH_SYSTEM == PSYCH_LINUX) {
                // Use USB bus-id and device-id as unique location identifier:
                sscanf(hid_dev->path, "%x:%x", &busId, &devId);
                currentDevice->locID = (double)((busId << 16) + (devId << 0));
            }

            if (PSYCH_SYSTEM == PSYCH_WINDOWS) {
                // Use device container id as unique location identifier.
                // This may only work on Windows-7+ and is a bit of a hack here,
                // the id is a GUID, nothing related to busId or devId. We init
                // devId with the hid_dev pointer value, to get a devId and thereby
                // location id in case proper parsing of a container id doesn't work:
                busId = 0;
                devId = (int)(size_t)hid_dev;
                if (strstr(hid_dev->path, "{")) sscanf(strstr(hid_dev->path, "{"), "{%x-%x", &busId, &devId);
                currentDevice->locID = (double)(((psych_uint64)busId << 32) + devId);
            }

            // Interface number is great for identifying DAQ devices, but not available
            // on OS/X, so this will be a Linux/Windows only thing.
            currentDevice->interfaceId = hid_dev->interface_number;

            // Enqueue record into linked list:
            currentDevice->pNext = hid_devices;
            hid_devices = currentDevice;
        }
    }

    return;
}

pRecDevice HIDGetFirstDevice(void)
{
    return hid_devices;
}

pRecDevice HIDGetNextDevice(pRecDevice pDevice)
{
    return pDevice->pNext;
}

/* HIDCountDevices(): Return count of all enumerated HID devices: */
psych_uint32 HIDCountDevices(void)
{
    pRecDevice cur_dev = hid_devices;
    psych_uint32 count = 0;

    while (cur_dev) {
        count++;
        cur_dev = cur_dev->pNext;
    }

    return count;
}

#endif
