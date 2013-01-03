/*
 PsychToolbox3/Source/Windows/PsychHID/PsychHIDStandardInterfaces.c
 
 PROJECTS: PsychHID only.
 
 PLATFORMS:  Windows.
 
 AUTHORS:
 mario.kleiner@tuebingen.mpg.de    mk
 
 HISTORY:
 9.08.2011     mk     Created.
 
 TO DO:
 * Classic KbCheck for multiple keyboards.
 * Gamepad support.
 * Code cleanup.
 * Maybe use 32-Bit timestamps delivered from the os instead of rolling our own?

*/

#include "PsychHIDStandardInterfaces.h"

#ifdef __cplusplus
extern "C" {
#endif

static int ndevices = 0;

// Pointer to DirectInput-8 interface:
LPDIRECTINPUT8 dinput = NULL;

// GUID identifiers of all enumerated keyboards:
typedef struct dinfo {
    GUID guidInstance;
    GUID guidProduct;
    psych_uint16 usagePage;
    psych_uint16 usageValue;
    psych_uint32 dwDevType;
    TCHAR tszInstanceName[MAX_PATH];
    TCHAR tszProductName[MAX_PATH];
} dinfo;

static dinfo info[PSYCH_HID_MAX_DEVICES];
LPDIRECTINPUTDEVICE8 x_dev[PSYCH_HID_MAX_DEVICES];

static	double* psychHIDKbQueueFirstPress[PSYCH_HID_MAX_DEVICES];
static	double* psychHIDKbQueueFirstRelease[PSYCH_HID_MAX_DEVICES];
static	double* psychHIDKbQueueLastPress[PSYCH_HID_MAX_DEVICES];
static	double* psychHIDKbQueueLastRelease[PSYCH_HID_MAX_DEVICES];
static  int*    psychHIDKbQueueScanKeys[PSYCH_HID_MAX_DEVICES];
static  psych_bool psychHIDKbQueueActive[PSYCH_HID_MAX_DEVICES];
static  psych_mutex KbQueueMutex;
static  psych_condition KbQueueCondition;
static	psych_bool  KbQueueThreadTerminate;
static  psych_thread KbQueueThread;
static  HANDLE hEvent;

// Enumeration callback from DirectInput-8 for keyboard enumeration.
// Called once for each detected keyboard. lpddi points to the following
// info structure:
//
//typedef struct DIDEVICEINSTANCE {
//    DWORD dwSize;
//    GUID guidInstance;
//    GUID guidProduct;
//    DWORD dwDevType;
//    TCHAR tszInstanceName[MAX_PATH];
//    TCHAR tszProductName[MAX_PATH];
//    GUID guidFFDriver;
//    WORD wUsagePage;
//    WORD wUsage;
//} DIDEVICEINSTANCE, *LPDIDEVICEINSTANCE;
BOOL keyboardEnumCallback(LPCDIDEVICEINSTANCE lpddi, LPVOID pvRef)
{
    // Useless assignment to make compiler happy:
    pvRef = NULL;
    
    // Copy relevant info into our own info array:
    info[ndevices].guidInstance = lpddi->guidInstance;
    info[ndevices].guidProduct = lpddi->guidProduct;
    info[ndevices].usagePage = lpddi->wUsagePage;
    info[ndevices].usageValue = lpddi->wUsage;
    info[ndevices].dwDevType = lpddi->dwDevType;
    memcpy(&info[ndevices].tszInstanceName[0], &(lpddi->tszInstanceName[0]), MAX_PATH * sizeof(TCHAR));
    memcpy(&info[ndevices].tszProductName[0], &(lpddi->tszProductName[0]), MAX_PATH * sizeof(TCHAR));
    
    ndevices++;
    
    // Done. Continue with enumeration, unless the capacity of our internal
    // array is exhausted:
    if (ndevices == PSYCH_HID_MAX_DEVICES) printf("PsychHID-WARNING: Number of detected HID input devices %i now equal to our maximum capacity. May miss some devices!\n", ndevices);    
    
    return((ndevices < PSYCH_HID_MAX_DEVICES) ? DIENUM_CONTINUE : DIENUM_STOP);
}

void PsychHIDInitializeHIDStandardInterfaces(void)
{
	int i;
    HRESULT rc;
    HINSTANCE modulehandle = NULL;
    dinput = NULL;
    ndevices = 0;

	// Init x_dev array:
	for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) x_dev[i] = NULL;
    
	// Init keyboard queue arrays:
	memset(&psychHIDKbQueueFirstPress[0], 0, sizeof(psychHIDKbQueueFirstPress));
	memset(&psychHIDKbQueueFirstRelease[0], 0, sizeof(psychHIDKbQueueFirstRelease));
	memset(&psychHIDKbQueueLastPress[0], 0, sizeof(psychHIDKbQueueLastPress));
	memset(&psychHIDKbQueueLastRelease[0], 0, sizeof(psychHIDKbQueueLastRelease));
	memset(&psychHIDKbQueueActive[0], 0, sizeof(psychHIDKbQueueActive));
	memset(&psychHIDKbQueueScanKeys[0], 0, sizeof(psychHIDKbQueueScanKeys));

    // We need the module instance handle of ourselves, ie., the PsychHID mex file to
    // open a DirectInput-8 interface, so the OS can apply backwards compatibility fixes
    // specific to the way our mex file DLL was built. For this we need the name of the
    // mex file, which is dependent on Octave vs. Matlab and 32-Bit vs. 64-Bit:
    #ifndef PTBOCTAVE3MEX
        // Matlab: 64-Bit or 32-Bit mex file?
        #if defined(__LP64__) || defined(_M_IA64) || defined(_WIN64)
            // 64-Bit:
            modulehandle = GetModuleHandle("PsychHID.mexw64");
        #else
            // 32-Bit:
            modulehandle = GetModuleHandle("PsychHID.mexw32");
        #endif
    #else
        // Octave: Same mex file file-extension for 32/64-Bit:
        modulehandle = GetModuleHandle("PsychHID.mex");
    #endif
    
    // If this doesn't work, try with application module handle as fallback. This works usually on
    // Windows XP/Vista/7, but may fail catastrophically on Windows-8 and later:
    if (NULL == modulehandle) {
        printf("PsychHID-WARNING: Could not get module handle to PsychHID mex file. Did you rename it? Please don't do that!\n");
        printf("PsychHID-WARNING: Will try application module handle as fallback. This may end badly, e.g., with a crash. Cross your fingers!\n");
        modulehandle = GetModuleHandle(NULL);
    }

    if (NULL == modulehandle) PsychErrorExitMsg(PsychError_system, "PsychHID: FATAL ERROR: Couldn't get module handle to create interface to Microsoft DirectInput-8! Game over!");
    
    // Open a DirectInput-8 interface:
    rc = DirectInput8Create(modulehandle, DIRECTINPUT_VERSION, IID_IDirectInput8, (LPVOID*)&dinput, NULL);
    
	if (DI_OK != rc) {
        printf("PsychHID-ERROR: Error return from DirectInput8Create: %x\n", (int) rc);
        if (rc == DIERR_OLDDIRECTINPUTVERSION) printf("PsychHID-ERROR: You need to install a more recent version of DirectX -- at least DirectX-8.\n");
		PsychErrorExitMsg(PsychError_system, "PsychHID: FATAL ERROR: Couldn't create interface to Microsoft DirectInput-8! Game over!");
	}
    
	// Enumerate all DirectInput keyboard(-like) devices:
    rc = dinput->EnumDevices(DI8DEVCLASS_KEYBOARD, (LPDIENUMDEVICESCALLBACK) keyboardEnumCallback, NULL, DIEDFL_ATTACHEDONLY | DIEDFL_INCLUDEHIDDEN);
    if (DI_OK != rc) {
        printf("PsychHID-ERROR: Error return from DirectInput8 EnumDevices(): %x! Game over!\n", (int) rc);
        goto out;
	}

	// Enumerate all DirectInput mouse(-like) devices:
    rc = dinput->EnumDevices(DI8DEVCLASS_POINTER, (LPDIENUMDEVICESCALLBACK) keyboardEnumCallback, NULL, DIEDFL_ATTACHEDONLY | DIEDFL_INCLUDEHIDDEN);
    if (DI_OK != rc) {
        printf("PsychHID-ERROR: Error return from DirectInput8 EnumDevices(): %x! Game over!\n", (int) rc);
        goto out;
	}

	// Enumerate all DirectInput joystick/gamepad(-like) devices:
    rc = dinput->EnumDevices(DI8DEVCLASS_GAMECTRL, (LPDIENUMDEVICESCALLBACK) keyboardEnumCallback, NULL, DIEDFL_ATTACHEDONLY | DIEDFL_INCLUDEHIDDEN);
    if (DI_OK != rc) {
        printf("PsychHID-ERROR: Error return from DirectInput8 EnumDevices(): %x! Game over!\n", (int) rc);
        goto out;
	}
    
	// Create keyboard queue mutex for later use:
	KbQueueThreadTerminate = FALSE;
	PsychInitMutex(&KbQueueMutex);
	PsychInitCondition(&KbQueueCondition, NULL);
    
    // Create event object for signalling device state changes:
    hEvent = CreateEvent(   NULL,	// default security attributes
                            FALSE,	// auto-reset event: This would need to be set TRUE for PsychBroadcastCondition() to work on Windows!
                            FALSE,	// initial state is nonsignaled
                            NULL	// no object name
                            ); 
    
    // Ready.
	return;
    
out:
        ndevices = 0;
    
	// Close our dedicated x-display connection and we are done:
    if (dinput) dinput->Release();
    dinput = NULL;
    
	PsychErrorExitMsg(PsychError_system, "PsychHID: FATAL ERROR: X Input extension version 2.0 or later not available! Game over!");	
}

static LPDIRECTINPUTDEVICE8 GetXDevice(int deviceIndex)
{
	if (deviceIndex < 0 || deviceIndex >= PSYCH_HID_MAX_DEVICES || deviceIndex >= ndevices) PsychErrorExitMsg(PsychError_user, "Invalid deviceIndex specified. No such device!");
	if (x_dev[deviceIndex] == NULL) {
        if (DI_OK != dinput->CreateDevice(info[deviceIndex].guidInstance, &(x_dev[deviceIndex]), NULL)) {
            PsychErrorExitMsg(PsychError_user, "Could not open connection to device. CreateDevice() failed!");
        }
    }
	return(x_dev[deviceIndex]);
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
    
	// Close all devices registered in x_dev array:
	for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
		if (x_dev[i]) x_dev[i]->Release();
		x_dev[i] = NULL;
	}
    
	// Release keyboard queue mutex:
	PsychDestroyMutex(&KbQueueMutex);
	PsychDestroyCondition(&KbQueueCondition);
	KbQueueThreadTerminate = FALSE;
    
    if (!CloseHandle(hEvent)) {
		printf("PsychHID-WARNING: Closing keyboard event handle failed!\n");
    }
    
    ndevices = 0;
    
	// Close our dedicated x-display connection and we are done:
    if (dinput) dinput->Release();
    dinput = NULL;
    
	return;
}

PsychError PsychHIDEnumerateHIDInputDevices(int deviceClass)
{
    const char *deviceFieldNames[]={"usagePageValue", "usageValue", "usageName", "index", "transport", "vendorID", "productID", "version", 
        "manufacturer", "product", "serialNumber", "locationID", "interfaceID", "totalElements", "features", "inputs", 
        "outputs", "collections", "axes", "buttons", "hats", "sliders", "dials", "wheels"};
    int numDeviceStructElements, numDeviceStructFieldNames=24, deviceIndex;
    PsychGenericScriptType	*deviceStruct;
    dinfo *dev;
    int i;
    int numKeys, numAxis;
    char *type = "";
    
    // Preparse: Count matching devices for deviceClass
    numDeviceStructElements = ndevices;
    
    // Alloc struct array of sufficient size:
    PsychAllocOutStructArray(1, kPsychArgOptional, numDeviceStructElements, numDeviceStructFieldNames, deviceFieldNames, &deviceStruct);
    deviceIndex = 0;
    
    // Return info:
    for(i = 0; i < ndevices; i++) {
        // Check i'th device:
        dev = &info[i];
        
        switch(dev->dwDevType & 0xff) {
            case DI8DEVTYPE_MOUSE:
			case DI8DEVTYPE_SCREENPOINTER:
				type = "slave pointer";
				if (dev->usagePage == 0) dev->usagePage = 1;
				if (dev->usageValue == 0) dev->usageValue = 2;
			break;

			case DI8DEVTYPE_KEYBOARD:
				type = "slave keyboard";
				if (dev->usagePage == 0) dev->usagePage = 1;
				if (dev->usageValue == 0) dev->usageValue = 6;
			break;
			
			case DI8DEVTYPE_JOYSTICK:
				type = "slave joystick";
				if (dev->usagePage == 0) dev->usagePage = 1;
				if (dev->usageValue == 0) dev->usageValue = 4;
			break;			
        }
        
        PsychSetStructArrayDoubleElement("usagePageValue",	deviceIndex, 	(double) dev->usagePage, deviceStruct);        
        PsychSetStructArrayDoubleElement("usageValue",	deviceIndex,        (double) dev->usageValue, deviceStruct);
        
        PsychSetStructArrayStringElement("usageName",		deviceIndex, 	type, deviceStruct);
        PsychSetStructArrayDoubleElement("index",           deviceIndex, 	(double) i, deviceStruct);
        PsychSetStructArrayStringElement("transport",		deviceIndex, 	dev->tszInstanceName, deviceStruct);
        PsychSetStructArrayStringElement("product",         deviceIndex, 	dev->tszProductName, deviceStruct);
        PsychSetStructArrayDoubleElement("locationID",		deviceIndex, 	(double) -1, deviceStruct);
        PsychSetStructArrayDoubleElement("interfaceID",		deviceIndex, 	(double) -1, deviceStruct);
        PsychSetStructArrayDoubleElement("productID",		deviceIndex, 	(double) -1, deviceStruct);
        
        //PsychSetStructArrayDoubleElement("vendorID",		deviceIndex, 	(double)currentDevice->vendorID, 	deviceStruct);
        //PsychSetStructArrayDoubleElement("version",		deviceIndex, 	(double)currentDevice->version, 	deviceStruct);
        //PsychSetStructArrayStringElement("manufacturer",	deviceIndex, 	currentDevice->manufacturer, 		deviceStruct);
        //PsychSetStructArrayStringElement("serialNumber",	deviceIndex, 	currentDevice->serial, 			deviceStruct);
        
        numKeys = numAxis = 0;
        
        PsychSetStructArrayDoubleElement("totalElements",	deviceIndex, 	(double) numKeys + numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("features",		deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("inputs",          deviceIndex, 	(double) numKeys + numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("outputs",         deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("collections",     deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("axes",		deviceIndex, 	(double) numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("buttons",		deviceIndex, 	(double) numKeys, deviceStruct);
        PsychSetStructArrayDoubleElement("hats",		deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("sliders",		deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("dials",		deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("wheels",		deviceIndex, 	(double) 0, deviceStruct);
        deviceIndex++;
    }
    
    return(PsychError_none);
}

int PsychHIDGetDefaultKbQueueDevice(void)
{
    // Return first enumerated keyboard (index == 0) if any available:
    if (ndevices > 0) return(0);
    
    // Nothing found? If so, abort:
    PsychErrorExitMsg(PsychError_user, "Could not find any useable keyboard device or other default input device!");

	// Utterly bogus return to make crappy Microsoft compiler shut up:
	return(0);
}

static unsigned int PsychHIDOSMapKey(unsigned int inkeycode)
{
	unsigned int keycode;

	// Must also shift by one count to account for difference 1-based vs. zero-based indexing:
	#ifndef MAPVK_VSC_TO_VK_EX
	#define MAPVK_VSC_TO_VK_EX 3
	#endif

	keycode = MapVirtualKeyEx(inkeycode, MAPVK_VSC_TO_VK_EX, GetKeyboardLayout(0));
	if (keycode == 0) {
		// Untranslated keycode: Use table.
		switch(inkeycode) {
			case DIK_RCONTROL:
				keycode = 163;
			break;
			case DIK_RALT:
				keycode = 165;
			break;
			case DIK_LWIN:
				keycode = 91;
			break;
			case DIK_RWIN:
				keycode = 92;
			break;
			case DIK_LEFT:
				keycode = 37;
			break;
			case DIK_RIGHT:
				keycode = 39;
			break;
			case DIK_UP:
				keycode = 38;
			break;
			case DIK_DOWN:
				keycode = 40;
			break;
			case DIK_PRIOR:
				keycode = 33;
			break;
			case DIK_NEXT:
				keycode = 34;
			break;
			case DIK_HOME:
				keycode = 36;
			break;
			case DIK_END:
				keycode = 35;
			break;
			case DIK_INSERT:
				keycode = 45;
			break;
			case DIK_DELETE:
				keycode = 46;
			break;
			case DIK_DIVIDE:
				keycode = 111;
			break;
			case DIK_NUMPADENTER:
				keycode = 13;
			break;
			case DIK_PAUSE:
				keycode = 19;
			break;
			case DIK_SYSRQ:
				keycode = 44;
			break;
			case DIK_APPS:
				keycode = 93;
			break;

			default:
				keycode = inkeycode & 0xff;
		}
	}
    
    // This is a table of overrides for the numeric keypad. It makes sure that
    // we always get keys from the keypad as if the "NumLock" key is active, aka
    // numbers instead of other stuff:
    switch(inkeycode) {
        case DIK_NUMPAD0:
            keycode = 96;
        break;
        case DIK_NUMPAD1:
            keycode = 97;
        break;
        case DIK_NUMPAD2:
            keycode = 98;
        break;
        case DIK_NUMPAD3:
            keycode = 99;
        break;
        case DIK_NUMPAD4:
            keycode = 100;
        break;
        case DIK_NUMPAD5:
            keycode = 101;
        break;
        case DIK_NUMPAD6:
            keycode = 102;
        break;
        case DIK_NUMPAD7:
            keycode = 103;
        break;
        case DIK_NUMPAD8:
            keycode = 104;
        break;
        case DIK_NUMPAD9:
            keycode = 105;
        break;
        case DIK_DECIMAL:
            keycode = 110;
        break;
    }

	// Adapt for 1 offset:
	if (keycode > 0) keycode--;

	return(keycode);
}

PsychError PsychHIDOSKbCheck(int deviceIndex, double* scanList)
{    
	psych_uint8 keys[1024];
    LPDIRECTINPUTDEVICE8 kb;
	unsigned int i, j;
	double* buttonStates;
	int keysdown;
	double timestamp;
    DWORD cbSize;

	if (deviceIndex == INT_MAX) {
		deviceIndex = PsychHIDGetDefaultKbQueueDevice();
		// Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
	}
    
	if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
		// Out of range index:
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
	}

	// Get DirectInput keyboard device:
	kb = GetXDevice(deviceIndex);    

	// Keyboard queue for this deviceIndex already exists?
	if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
		// No. Create one which accepts all keys:
		PsychHIDOSKbQueueCreate(deviceIndex, 0, NULL);
	}

	// Keyboard queue for this device active? If not, we need
	// to start it:
	if (!psychHIDKbQueueActive[deviceIndex]) {
		// Keyboard idle: Need to start it:
		PsychHIDOSKbQueueStart(deviceIndex);

		// Startup to first key delivery takes time. Wait for
		// 50 msecs to be on the safe side:
		PsychYieldIntervalSeconds(0.050);
	}

	// Size of state structure is device dependent:
	switch (info[deviceIndex].dwDevType & 0xff) {
		case DI8DEVTYPE_KEYBOARD:
			cbSize = 256;
		break;

		case DI8DEVTYPE_MOUSE:
		case DI8DEVTYPE_SCREENPOINTER:
			cbSize = sizeof(DIMOUSESTATE2);
		break;

		case DI8DEVTYPE_JOYSTICK:
			cbSize = sizeof(DIJOYSTATE2);
		break;

		default: // Unkown device. Fail.
			cbSize = 0;
	}


	// Query current state snapshot of keyboard:
	memset(keys, 0, sizeof(keys));
	if (DI_OK != kb->GetDeviceState(cbSize, (LPVOID) &keys[0])) {
		printf("PsychHID-ERROR: KbCheck for deviceIndex %i failed, because query of device failed!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "KbCheck failed!");
	}

	// Request current time of query:
	PsychGetAdjustedPrecisionTimerSeconds(&timestamp);

	// Reset overall key state to "none pressed":
	keysdown = 0;

	// Copy out timestamp:
	PsychCopyOutDoubleArg(2, kPsychArgOptional, timestamp);

	// Copy keyboard state:
	PsychAllocOutDoubleMatArg(3, kPsychArgOptional, 1, 256, 1, &buttonStates);
	for (i = 0; i < 256; i++) buttonStates[i] = 0;

	// Keyboard?
	if (cbSize == 256) {
		// Copy button state to output vector, apply scanlist mask, compute
		// resulting overall keysdown state. We ignore keyboard scancode zero and
		// start with 1 instead. We also ignore code 255. These are borderline codes
		// which may do weird things...
		for (i = 1; i < 255; i++) {
			// Compute target key slot for this scancode i:
			j = PsychHIDOSMapKey(i);

			// This key down?
			buttonStates[j] += (keys[i] > 0) ? 1 : 0;
			// Apply scanList mask, if any provided:
			if (scanList && (scanList[j] <= 0)) buttonStates[j] = 0;
			keysdown += (unsigned int) buttonStates[j];
		}
	}
	
	// Joystick?
	if (cbSize == sizeof(DIJOYSTATE2)) {
		// Copy button state to output vector, apply scanlist mask, compute
		// resulting overall keysdown state. There are 128 buttons at an offset:
		for (i = (8 * sizeof(LONG) + 4 * sizeof(DWORD)); i < (8 * sizeof(LONG) + 4 * sizeof(DWORD)) + 128; i++) {
			// Compute target key slot for this scancode i:
			j = i - (8 * sizeof(LONG) + 4 * sizeof(DWORD));

			// This key down?
			buttonStates[j] += (keys[i] > 0) ? 1 : 0;
			// Apply scanList mask, if any provided:
			if (scanList && (scanList[j] <= 0)) buttonStates[j] = 0;
			keysdown += (unsigned int) buttonStates[j];
		}
	}
	
	// Mouse?
	if (cbSize == sizeof(DIMOUSESTATE2)) {
		// Copy button state to output vector, apply scanlist mask, compute
		// resulting overall keysdown state. There are 8 buttons at an offset:
		for (i = (3 * sizeof(LONG)); i < (3 * sizeof(LONG)) + 8; i++) {
			// Compute target key slot for this scancode i:
			j = i - (3 * sizeof(LONG));

			// This key down?
			buttonStates[j] += (keys[i] > 0) ? 1 : 0;
			// Apply scanList mask, if any provided:
			if (scanList && (scanList[j] <= 0)) buttonStates[j] = 0;
			keysdown += (unsigned int) buttonStates[j];
		}
	}
	
	// Copy out overall keystate:
	PsychCopyOutDoubleArg(1, kPsychArgOptional, (keysdown > 0) ? 1 : 0);

	return(PsychError_none);
}

PsychError PsychHIDOSGamePadAxisQuery(int deviceIndex, int axisId, double* min, double* max, double* val, char* axisLabel)
{
    PsychErrorExitMsg(PsychError_unimplemented, "This function is not yet implemented for MS-Windows, sorry!");
	return(PsychError_none);
}

// This is the event dequeue & process function which updates
// Keyboard queue state. It can be called with 'blockingSinglepass'
// set to TRUE to process exactly one event, if called from the
// background keyboard queue processing thread. Alternatively it
// can be called synchronously from KbQueueCheck with a setting of FALSE
// to iterate over all available events and process them instantaneously:
void KbQueueProcessEvents(psych_bool blockingSinglepass)
{
    LPDIRECTINPUTDEVICE8 kb;
    DIDEVICEOBJECTDATA event;
    HRESULT rc;
    DWORD dwItems; 
	double tnow;
	unsigned int i, keycode, keystate;
 	PsychHIDEventRecord evt;
	WORD asciiValue[2];
	UCHAR keyboardState[256];

	while (1) {        
		// Single pass or multi-pass?
		if (blockingSinglepass) {
			// Wait until at least one event available and dequeue it:
            // We use a timeout of 100 msecs.
            WaitForSingleObject(hEvent, 100);
		} else {
			// Check if event available, dequeue it, if so. Abort
			// processing if no new event available, aka queue empty:
			// TODO if (!XCheckTypedEvent(thread_dpy, GenericEvent, &KbQueue_xevent)) break;
		}
        
		// Take timestamp:
		PsychGetAdjustedPrecisionTimerSeconds(&tnow);

        // Need the lock from here on:
        PsychLockMutex(&KbQueueMutex);
        
        // Do a sweep over all keyboard devices whose queues are active:
        for (i = 0; i < (unsigned int) ndevices; i++) {
            // Skip this one if inactive:
            if (!psychHIDKbQueueActive[i]) continue;
            
            // Check this device:
            kb = GetXDevice(i);
            
            // Fetch one item from the buffer:
            // event.dwTimeStamp = Timestamp in msecs of timeGetTime() timebase.
            // event.dwSequence = Sequence number.
            
            // Fetch from this device, item-by-item, until nothing more to fetch:
            while (TRUE) {
                // Try one fetch from this device:
                dwItems = 1;
                rc = kb->GetDeviceData(sizeof(DIDEVICEOBJECTDATA), &event, &dwItems, 0);
                
                // If failed or nothing more to fetch, break out of fetch loop:
                if (!SUCCEEDED(rc) || (0 == dwItems)) break;

				// Clear ringbuffer event:
				memset(&evt, 0 , sizeof(evt));

				// Init character code to "unmapped": It will stay like that for anything but real keyboards:
				evt.cookedEventCode = -1;

                // Map to key code and key state:
                keycode = event.dwOfs & 0xff;
                keystate = event.dwData & 0x80;
                
				// Remap keycode into target slot in our arrays, depending on input device:
                switch (info[i].dwDevType & 0xff) {
					case DI8DEVTYPE_KEYBOARD:
						// Try to map scancode to ascii character:
                        memset(keyboardState, 0, sizeof(keyboardState));
                        if (GetAsyncKeyState(VK_SHIFT)) keyboardState[VK_SHIFT] = 0xff;
                        
                        if ((1 == ToAsciiEx(MapVirtualKeyEx(keycode, 1, GetKeyboardLayout(0)), keycode, keyboardState, (LPWORD) &(asciiValue[0]), 0, GetKeyboardLayout(0)))) {
							// Mapped to single char: Return it as cooked keycode:
							evt.cookedEventCode = (int) (asciiValue[0] & 0xff);
						}
						else {
							// Could not map key to valid ascii character: Mark as "not mapped" aka zero:
							evt.cookedEventCode = 0;
						}

                        // Map scancode 'keycode' to virtual key code 'keycode':
						keycode = PsychHIDOSMapKey(keycode);
					break;

					case DI8DEVTYPE_MOUSE:
					case DI8DEVTYPE_SCREENPOINTER:
						// Button event? Otherwise skip it.
						if (keycode < 3 * sizeof(LONG)) continue;
						// Correct for buttons offset in data structure DIMOUSESTATE2:
						keycode -= 3 * sizeof(LONG);
					break;

					case DI8DEVTYPE_JOYSTICK:
						// Button event? Otherwise skip it.
						if (keycode < (8 * sizeof(LONG) + 4 * sizeof(DWORD))) continue;
						// Correct for buttons offset in data structure DIJOYSTATE2:
						keycode -= (8 * sizeof(LONG) + 4 * sizeof(DWORD));
						// Also skip if beyond button array:
						if (keycode >= 128) continue;
					break;
				
					default: // Unkown device -- Skip it.
						continue;
                }

                // This keyboard queue interested in this keycode?
                if (psychHIDKbQueueScanKeys[i][keycode] != 0) {
                    // Yes: The queue wants to receive info about this key event.
                    
                    // Press or release?
                    if (keystate) {
                        // Enqueue key press. Always in the "last press" array, because any
                        // press at this time is the best candidate for the last press.
                        // Only enqeue in "first press" if there wasn't any registered before,
                        // ie., the slot is so far empty:
                        if (psychHIDKbQueueFirstPress[i][keycode] == 0) psychHIDKbQueueFirstPress[i][keycode] = tnow;
                        psychHIDKbQueueLastPress[i][keycode] = tnow;
						evt.status |= (1 << 0);
                    } else {
                        // Enqueue key release. See logic above:
                        if (psychHIDKbQueueFirstRelease[i][keycode] == 0) psychHIDKbQueueFirstRelease[i][keycode] = tnow;
                        psychHIDKbQueueLastRelease[i][keycode] = tnow;
						evt.status &= ~(1 << 0);
						// Clear cooked keycode - We don't record key releases this way:
						if (evt.cookedEventCode > 0) evt.cookedEventCode = 0;
                    }

					// Update event buffer:
					evt.timestamp = tnow;
					evt.rawEventCode = keycode + 1;
					PsychHIDAddEventToEventBuffer(i, &evt);
                    
                    // Tell waiting userspace (under KbQueueMutex protection for better scheduling) something interesting has changed:
                    PsychSignalCondition(&KbQueueCondition);
                }
                // Next fetch iteration for this device...
            }
            // Check next device...
        }
        
        // Done with shared data access:
        PsychUnlockMutex(&KbQueueMutex);
        
		// Done if we were only supposed to handle one sweep, which we did:
		if (blockingSinglepass) break;
	}
    
	return;
}

// Async processing thread for keyboard events:
void* KbQueueWorkerThreadMain(void* dummy)
{
	int rc;
    
	// Try to raise our priority: We ask to switch ourselves (NULL) to priority class 1 aka
	// realtime scheduling, with a tweakPriority of +1, ie., raise the relative
	// priority level by +1 wrt. to the current level:
	if ((rc = PsychSetThreadPriority(NULL, 1, 1)) > 0) {
		printf("PsychHID: KbQueueStart: Failed to switch to realtime priority [%s].\n", strerror(rc));
	}
    
	while (1) {
		PsychLockMutex(&KbQueueMutex);
        
		// Check if we should terminate:
		if (KbQueueThreadTerminate) break;
        
		PsychUnlockMutex(&KbQueueMutex);
        
		// Perform event processing until no more events are pending:
		KbQueueProcessEvents(TRUE);
	}
    
	// Done. Unlock the mutex:
	PsychUnlockMutex(&KbQueueMutex);
    
	// printf("DEBUG: THREAD TERMINATING...\n"); fflush(NULL);
    
	// Return and terminate:
	return(NULL);
}

PsychError PsychHIDOSKbQueueCreate(int deviceIndex, int numScankeys, int* scanKeys)
{
	dinfo* dev = NULL;
    
	// Valid number of keys?
	if (scanKeys && (numScankeys != 256)) {
		PsychErrorExitMsg(PsychError_user, "Second argument to KbQueueCreate must be a vector with 256 elements.");
	}
    
	if (deviceIndex < 0) {
		deviceIndex = PsychHIDGetDefaultKbQueueDevice();
		// Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
	} else if (deviceIndex >= ndevices) {
		// Out of range index:
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
	} 
    
	// Do we finally have a valid keyboard?
	dev = &info[deviceIndex];
	if ((dev->dwDevType & 0xff) != DI8DEVTYPE_KEYBOARD) {
		//		PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. Not a keyboard device!");
	}
    
	// Keyboard queue for this deviceIndex already created?
	if (psychHIDKbQueueFirstPress[deviceIndex]) {
		// Yep. Release it, so we can start from scratch:
		PsychHIDOSKbQueueRelease(deviceIndex);
	}
    
	// Allocate and zero-init memory for tracking key presses and key releases:
	psychHIDKbQueueFirstPress[deviceIndex]   = (double*) calloc(256, sizeof(double));
	psychHIDKbQueueFirstRelease[deviceIndex] = (double*) calloc(256, sizeof(double));
	psychHIDKbQueueLastPress[deviceIndex]    = (double*) calloc(256, sizeof(double));
	psychHIDKbQueueLastRelease[deviceIndex]  = (double*) calloc(256, sizeof(double));
	psychHIDKbQueueScanKeys[deviceIndex]     = (int*) calloc(256, sizeof(int));
    
	// Assign scanKeys vector, if any:
	if (scanKeys) {
		// Copy it:
		memcpy(psychHIDKbQueueScanKeys[deviceIndex], scanKeys, 256 * sizeof(int));
	} else {
		// None provided. Enable all keys by default:
		memset(psychHIDKbQueueScanKeys[deviceIndex], 1, 256 * sizeof(int));        
	}

	// Create event buffer:
	PsychHIDCreateEventBuffer(deviceIndex);
    
	// Ready to use this keybord queue.
	return(PsychError_none);
}

void PsychHIDOSKbQueueRelease(int deviceIndex)
{
	if (deviceIndex < 0) {
		deviceIndex = PsychHIDGetDefaultKbQueueDevice();
		// Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
	}
    
	if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
		// Out of range index:
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
	}
    
	// Keyboard queue for this deviceIndex already exists?
	if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
		// No. Nothing to do then.
		return;
	}
    
	// Ok, we have a keyboard queue. Stop any operation on it first:
	PsychHIDOSKbQueueStop(deviceIndex);
    
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
    LPDIRECTINPUTDEVICE8 kb;
	psych_bool queueActive;
	int i;
    
	if (deviceIndex < 0) {
		deviceIndex = PsychHIDGetDefaultKbQueueDevice();
		// Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
	}
    
	if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
		// Out of range index:
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
	}
    
	// Keyboard queue for this deviceIndex already exists?
	if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
		// No. Nothing to do then.
		return;
	}
    
	// Keyboard queue already stopped?
	if (!psychHIDKbQueueActive[deviceIndex]) return;
    
    // Get device:
    kb = GetXDevice(deviceIndex);    
    
	// Queue is active. Stop it:
	PsychLockMutex(&KbQueueMutex);

	// Release the device:
    if (DI_OK != kb->Unacquire()) {
        PsychUnlockMutex(&KbQueueMutex);
		printf("PsychHID-ERROR: Tried to stop processing on keyboard queue for deviceIndex %i, but releasing device failed!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "Stopping keyboard queue failed!");
    }

    // Disable state-change event notifications:
    if (DI_OK != kb->SetEventNotification(NULL)) {
        PsychUnlockMutex(&KbQueueMutex);
		printf("PsychHID-ERROR: Tried to stop processing on keyboard queue for deviceIndex %i, but disabling device state notifications failed!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "Stopping keyboard queue failed!");
    }

	// Mark queue logically stopped:
	psychHIDKbQueueActive[deviceIndex] = FALSE;
    
	PsychUnlockMutex(&KbQueueMutex);
    
	// Was this the last active queue?
	queueActive = FALSE;
	for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
		queueActive |= psychHIDKbQueueActive[i];
	}
    
	// If more queues are active then we're done:
	if (queueActive) return;
    
	// No more active queues. Shutdown the common processing thread:
	PsychLockMutex(&KbQueueMutex);
    
	KbQueueThreadTerminate = TRUE;
    
	// Done.
	PsychUnlockMutex(&KbQueueMutex);
    
	// Shutdown the thread, wait for its termination:
	PsychDeleteThread(&KbQueueThread);
	KbQueueThreadTerminate = FALSE;
    
	// printf("DEBUG: THREAD JOINED.\n"); fflush(NULL);
    
	return;
}

void PsychHIDOSKbQueueStart(int deviceIndex)
{
    LPDIRECTINPUTDEVICE8 kb;
    DIPROPDWORD  dipdw;
	psych_bool queueActive;
	int i;
    
	if (deviceIndex < 0) {
		deviceIndex = PsychHIDGetDefaultKbQueueDevice();
		// Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
	}
    
	if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
		// Out of range index:
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
	}
    
	// Does Keyboard queue for this deviceIndex already exist?
	if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
		// No. Bad bad...
		printf("PsychHID-ERROR: Tried to start processing on non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No queue for that device yet!");
	}
    
	// Keyboard queue already stopped? Then we ain't nothing to do:
	if (psychHIDKbQueueActive[deviceIndex]) return;
    
	// Queue is inactive. Start it:
    
	// Will this be the first active queue, ie., aren't there any queues running so far?
	queueActive = FALSE;
	for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
		queueActive |= psychHIDKbQueueActive[i];
	}
    
	PsychLockMutex(&KbQueueMutex);
    
	// Clear out current state for this queue:
	memset(psychHIDKbQueueFirstPress[deviceIndex]   , 0, (256 * sizeof(double)));
	memset(psychHIDKbQueueFirstRelease[deviceIndex] , 0, (256 * sizeof(double)));
	memset(psychHIDKbQueueLastPress[deviceIndex]    , 0, (256 * sizeof(double)));
	memset(psychHIDKbQueueLastRelease[deviceIndex]  , 0, (256 * sizeof(double)));
    
	// Setup event mask, so events from our associated xinput device
	// get enqueued in our event queue:
    kb = GetXDevice(deviceIndex);
	
	// Device specific data format setup:
	switch (info[deviceIndex].dwDevType & 0xff) {
		case DI8DEVTYPE_KEYBOARD:
			if (DI_OK != kb->SetDataFormat(&c_dfDIKeyboard)) {
				PsychUnlockMutex(&KbQueueMutex);
				printf("PsychHID-ERROR: Tried to start processing on keyboard queue for deviceIndex %i, but setting dataformat failed!\n", deviceIndex);
				PsychErrorExitMsg(PsychError_user, "Starting keyboard queue failed!");
			}			
		break;
			
		case DI8DEVTYPE_MOUSE:
		case DI8DEVTYPE_SCREENPOINTER:
			if (DI_OK != kb->SetDataFormat(&c_dfDIMouse2)) {
				PsychUnlockMutex(&KbQueueMutex);
				printf("PsychHID-ERROR: Tried to start processing on keyboard queue for deviceIndex %i, but setting dataformat failed!\n", deviceIndex);
				PsychErrorExitMsg(PsychError_user, "Starting keyboard queue failed!");
			}			
		break;
			
		case DI8DEVTYPE_JOYSTICK:
			if (DI_OK != kb->SetDataFormat(&c_dfDIJoystick2)) {
				PsychUnlockMutex(&KbQueueMutex);
				printf("PsychHID-ERROR: Tried to start processing on keyboard queue for deviceIndex %i, but setting dataformat failed!\n", deviceIndex);
				PsychErrorExitMsg(PsychError_user, "Starting keyboard queue failed!");
			}			
		break;				
	}

    // Set device event buffer size to 256 elements:
    dipdw.diph.dwSize = sizeof(DIPROPDWORD);
    dipdw.diph.dwHeaderSize = sizeof(DIPROPHEADER);
    dipdw.diph.dwObj = 0;
    dipdw.diph.dwHow = DIPH_DEVICE;
    dipdw.dwData = 256;
    
    if (DI_OK != kb->SetProperty(DIPROP_BUFFERSIZE, &dipdw.diph)) {
        PsychUnlockMutex(&KbQueueMutex);
		printf("PsychHID-ERROR: Tried to start processing on keyboard queue for deviceIndex %i, but setting buffersize on device failed!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "Starting keyboard queue failed!");
    }
    
    // Enable state-change event notifications:
    if (DI_OK != kb->SetEventNotification(hEvent)) {
        PsychUnlockMutex(&KbQueueMutex);
		printf("PsychHID-ERROR: Tried to start processing on keyboard queue for deviceIndex %i, but setting device state notifications failed!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "Starting keyboard queue failed!");
    }
    
    if (DI_OK != kb->Acquire()) {
        PsychUnlockMutex(&KbQueueMutex);
		printf("PsychHID-ERROR: Tried to start processing on keyboard queue for deviceIndex %i, but acquiring device failed!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "Starting keyboard queue failed!");
    }
    
	// Mark this queue as logically started:
	psychHIDKbQueueActive[deviceIndex] = TRUE;
    
	// Queue started.
	PsychUnlockMutex(&KbQueueMutex);
    
	// If other queues are already active then we're done:
	if (queueActive) return;
    
	// No other active queues. We are the first one.
    
	// Start the common processing thread for all queues:
	PsychLockMutex(&KbQueueMutex);
	KbQueueThreadTerminate = FALSE;
    
	if (PsychCreateThread(&KbQueueThread, NULL, KbQueueWorkerThreadMain, NULL)) {
		// We are soo screwed:
        
		// Cleanup the mess:
		psychHIDKbQueueActive[deviceIndex] = FALSE;
		PsychUnlockMutex(&KbQueueMutex);
        
		// Whine a little bit:
		printf("PsychHID-ERROR: Start of keyboard queue processing failed!\n");
		PsychErrorExitMsg(PsychError_system, "Creation of keyboard queue background processing thread failed!");
	}
    
	// Up and running, we're done!
	PsychUnlockMutex(&KbQueueMutex);
    
	return;
}

void PsychHIDOSKbQueueFlush(int deviceIndex)
{
    LPDIRECTINPUTDEVICE8 kb;
	HRESULT rc;
    DWORD dwItems = INFINITE;
    
	if (deviceIndex < 0) {
		deviceIndex = PsychHIDGetDefaultKbQueueDevice();
		// Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
	}
    
	if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
		// Out of range index:
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
	}
    
	// Does Keyboard queue for this deviceIndex already exist?
	if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
		// No. Bad bad...
		printf("PsychHID-ERROR: Tried to flush non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No queue for that device yet!");
	}
    
    kb = GetXDevice(deviceIndex);    
    
	// Clear out current state for this queue:
	PsychLockMutex(&KbQueueMutex);
    
    // Flush device buffer:
    rc = kb->GetDeviceData(sizeof(DIDEVICEOBJECTDATA), NULL, &dwItems, 0);
    
    // Clear our buffer:
	memset(psychHIDKbQueueFirstPress[deviceIndex]   , 0, (256 * sizeof(double)));
	memset(psychHIDKbQueueFirstRelease[deviceIndex] , 0, (256 * sizeof(double)));
	memset(psychHIDKbQueueLastPress[deviceIndex]    , 0, (256 * sizeof(double)));
	memset(psychHIDKbQueueLastRelease[deviceIndex]  , 0, (256 * sizeof(double)));
    
	PsychUnlockMutex(&KbQueueMutex);
    
    return;
}

void PsychHIDOSKbQueueCheck(int deviceIndex)
{
	double *hasKeyBeenDownOutput, *firstPressTimeOutput, *firstReleaseTimeOutput, *lastPressTimeOutput, *lastReleaseTimeOutput;
	psych_bool isFirstPressSpecified, isFirstReleaseSpecified, isLastPressSpecified, isLastReleaseSpecified;
	int i;
    
	if (deviceIndex < 0) {
		deviceIndex = PsychHIDGetDefaultKbQueueDevice();
		// Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
	}
    
	if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
		// Out of range index:
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
	}
    
	// Does Keyboard queue for this deviceIndex already exist?
	if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
		// No. Bad bad...
		printf("PsychHID-ERROR: Tried to check non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No queue for that device yet!");
	}
    
	// Allocate output
	PsychAllocOutDoubleArg(1, kPsychArgOptional, &hasKeyBeenDownOutput);
	isFirstPressSpecified = PsychAllocOutDoubleMatArg(2, kPsychArgOptional, 1, 256, 1, &firstPressTimeOutput);
	isFirstReleaseSpecified = PsychAllocOutDoubleMatArg(3, kPsychArgOptional, 1, 256, 1, &firstReleaseTimeOutput);
	isLastPressSpecified = PsychAllocOutDoubleMatArg(4, kPsychArgOptional, 1, 256, 1, &lastPressTimeOutput);
	isLastReleaseSpecified = PsychAllocOutDoubleMatArg(5, kPsychArgOptional, 1, 256, 1, &lastReleaseTimeOutput);
    
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
    int keyMask[256];
    int i;
    double t, tc;
    
    if (deviceIndex < 0) {
        deviceIndex = PsychHIDGetDefaultKbQueueDevice();
        // Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
    }
    
    if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
        // Out of range index:
        PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
    }
    
    if(psychHIDKbQueueFirstPress[deviceIndex]) PsychErrorExitMsg(PsychError_user, "A queue for this device is already running, you must call KbQueueRelease() before invoking KbTriggerWait.");
    
    // Create a keyboard queue for this deviceIndex:
    memset(&keyMask[0], 0, sizeof(keyMask));
    for (i = 0; i < numScankeys; i++) {
        if (scanKeys[i] < 1 || scanKeys[i] > 256) PsychErrorExitMsg(PsychError_user, "Invalid entry for triggerKey specified. Not in valid range 1 - 256!");
        keyMask[scanKeys[i] - 1] = 1;
    }
    
    // Create keyboard queue with proper mask:
    PsychHIDOSKbQueueCreate(deviceIndex, 256, &keyMask[0]);
    PsychHIDOSKbQueueStart(deviceIndex);
    
    PsychLockMutex(&KbQueueMutex);
    
    // Scan for trigger key:
    while (1) {
        // Wait until something changes in a keyboard queue:
        PsychWaitCondition(&KbQueueCondition, &KbQueueMutex);
        
        // Check if our queue had one of the dedicated trigger keys pressed:
        for (i = 0; i < numScankeys; i++) {
            // Break out of scan loop if key pressed:
            if (psychHIDKbQueueFirstPress[deviceIndex][scanKeys[i] - 1] != 0) break;
        }
        
        // Triggerkey pressed?
        if ((i < numScankeys) && (psychHIDKbQueueFirstPress[deviceIndex][scanKeys[i] - 1] != 0)) break;
        
        // No change for our trigger keys. Repeat scan loop.
    }

    // If we reach this point, we know some triggerkey has been pressed. As we aborted
    // the scan on detection of the first pressed key, we can't be certain we caught the
    // key with the earliest key press, maybe one of the untested keys was pressed even
    // earlier. Therefore do another pass over all keys to find the pressed one with the
    // earliest (minimum) pressed time:
    t = DBL_MAX;
    for (i = 0; i < numScankeys; i++) {
        tc = psychHIDKbQueueFirstPress[deviceIndex][scanKeys[i] - 1];
        if ((tc != 0) && (tc <= t)) t = tc;
    }

    // Done. Release the lock:
    PsychUnlockMutex(&KbQueueMutex);
    
    // Stop and release the queue:
    PsychHIDOSKbQueueStop(deviceIndex);
    PsychHIDOSKbQueueRelease(deviceIndex);
    
    // Return timestamp:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, t);
    
    return;
}

#ifdef __cplusplus
}
#endif
