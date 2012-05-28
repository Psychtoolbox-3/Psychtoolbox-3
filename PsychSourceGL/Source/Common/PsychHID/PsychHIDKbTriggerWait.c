/*
	PsychtoolboxGL/Source/Common/PsychHID/PsychHIDKbTriggerWait.c		
  
	PROJECTS: 
	
		PsychHID only.
  
	PLATFORMS:  
	
		All.
  
	AUTHORS:
	
		rwoods@ucla.edu                     rpw
        mario.kleiner@tuebingen.mpg.de      mk
      
	HISTORY:
		8/9/07  rpw		Created.
		12/17/09 rpw	Added keypad support
  
	NOTES:
	
		This incorporates elements from PsychHIDKbCheck and PsychHIDKbWait, but is distinct in that:
			1) It waits for a specific trigger key rather than any keypress
			2) It uses a queue to check for the trigger rather than simply checking to see if the key is down
			
		The use of a queue assures that brief trigger events will not be missed by virtue of being too brief
		
		Requires Mac OS X 10.3 or later. The matlab function KbTriggerWait.m screens away Mac OS X 10.2
		and earlier, and this code does nothing to verify the Mac OS X version
		
		This routine will only listen to a single device. It would be possible to design an analogous routine
		that would listen for a particular trigger key (or set of trigger keys) on multiple devices, but 
		this would involve setting up a queue for each device and then polling each device in turn in the 
		loop that watches for the trigger (or setting callback functions).
		
		By nature, the oldest events are discarded from the queue once it reaches its full capacity.
		Consequently, the time returned by this routine is not guaranteed to be the time of the first
		event if multiple events arrive in between pollings. However, this is probably still more accurate
		than polling for the key being down. If it is desirable to guarantee that the very first event is
		the one used to extract the timestamp, it would be best to set an event callout on the queue.
		Short of this, setting the queue depth to a value larger than the maximum expected number of events
		(if this number is known) when creating the queue would suffice.
		
		If one wanted to call a function like this repeatedly with the same trigger while minimizing
		the computation time associated with setting up the queue, this funcion could be rewritten
		to declare queue static, adding an input argument to govern creating and setting up the queue, 
		starting the queue, flushing the queue, polling the queue and destroying the queue, but this would put
		the burden for remembering to destroy the queue back on the Matlab user.

		---
    

*/

#include "PsychHID.h"

static char useString[]= "secs=PsychHID('KbTriggerWait', KeysUsage, [deviceNumber])";
static char synopsisString[] = 
        "Scan a keyboard, keypad, or other HID device with buttons, and wait for a trigger key press.\n"
        "NOTE: This function should not be called directly! Use KbTriggerWait() instead.\n"
        "By default the first keyboard device (the one with the lowest device number) is "
        "scanned. If no keyboard is found, the first keypad device is scanned, followed by other "
        "devices, e.g., mice.  Optionally, the deviceNumber of any keyboard or HID device may be specified.\n"
        "The 'KeysUsage' parameter must specify the keycode of a single key to wait for on "
        "OS/X. On Linux and Windows, 'KeysUsage' can be a vector of trigger key codes and the "
        "wait will finish as soon as at least one of the keys in the vector is pressed.\n"
        "On MS-Windows XP and later, it is currently not possible to enumerate different keyboards and mice "
        "separately. Therefore the 'deviceNumber' argument is mostly useless for keyboards and mice, usually you can "
        "only check the system keyboard or mouse.\n"
        "OS/X on 64-Bit runtimes: This function is unsupported. You must use KbTriggerWait() instead.\n";

static char seeAlsoString[] = "";
 
PsychError PSYCHHIDKbTriggerWait(void) 
{
    int	deviceIndex;
    int numScankeys;
    int* scanKeys;
    
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));  	//Specify trigger key code and the deviceNumber of the keyboard or keypad to scan.  

	// Identify the mandatory trigger array:
    if (!PsychAllocInIntegerListArg(1, TRUE, &numScankeys, &scanKeys)){
        PsychErrorExitMsg(PsychError_user, "Keycode is required.");
    }
	
    // Get optional deviceIndex:
    if (!PsychCopyInIntegerArg(2, FALSE, &deviceIndex)) deviceIndex = -1;

    // Execute:
    PsychHIDOSKbTriggerWait(deviceIndex, numScankeys, scanKeys);

    return(PsychError_none);
}

#if PSYCH_SYSTEM == PSYCH_OSX

#include "PsychHIDKbQueue.h"
#define NUMDEVICEUSAGES 7
extern HIDDataRef hidDataRef;

void PsychHIDOSKbTriggerWait(int deviceIndex, int numScankeys, int* scanKeys)
{
    pRecDevice	deviceRecord;
    int			i, numDeviceIndices;
    long		KeysUsagePage;			// This is the usage page of the target element: A key on a keyboard/keypad or a button.
	long		KeysUsage;				// This will contain the key code of the trigger key
	long		KbDeviceUsagePages[NUMDEVICEUSAGES]= {kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop};
	long		KbDeviceUsages[NUMDEVICEUSAGES]={kHIDUsage_GD_Keyboard, kHIDUsage_GD_Keypad, kHIDUsage_GD_Mouse, kHIDUsage_GD_Pointer, kHIDUsage_GD_Joystick, kHIDUsage_GD_GamePad, kHIDUsage_GD_MultiAxisController};
    int 		numDeviceUsages=NUMDEVICEUSAGES;
    int			deviceIndices[PSYCH_HID_MAX_KEYBOARD_DEVICES]; 
    pRecDevice	deviceRecords[PSYCH_HID_MAX_KEYBOARD_DEVICES];
    psych_bool	isDeviceSpecified, foundUserSpecifiedDevice;
    double		*timeValueOutput;
	
	IOHIDQueueInterface **queue;
	HRESULT result;
	IOHIDDeviceInterface122** interface=NULL;	// This requires Mac OS X 10.3 or higher
	IOReturn success;	
	IOHIDElementCookie triggerCookie;
    uint32_t usage;

    // Assign single supported trigger key:
    if (numScankeys != 1) PsychErrorExitMsg(PsychError_user, "Sorry, the OS/X version of KbTriggerWait only supports one trigger key.");
    KeysUsage = (long) scanKeys[0];

    // Optional deviceIndex specified?
    isDeviceSpecified = (deviceIndex >= 0) ? TRUE : FALSE;
    
    PsychHIDVerifyInit();	
	if(hidDataRef!=NULL) PsychErrorExitMsg(PsychError_user, "A queue is already running, you must call KbQueueRelease() before invoking KbTriggerWait.");

    //Choose the device index and its record
    PsychHIDGetDeviceListByUsages(numDeviceUsages, KbDeviceUsagePages, KbDeviceUsages, &numDeviceIndices, deviceIndices, deviceRecords);  
    
    if(isDeviceSpecified){  //make sure that the device number provided by the user is really a keyboard or keypad.
        for(i=0;i<numDeviceIndices;i++){
            if(foundUserSpecifiedDevice=(deviceIndices[i]==deviceIndex))
                break;
        }
        if(!foundUserSpecifiedDevice)
            PsychErrorExitMsg(PsychError_user, "Specified device number is not a keyboard or keypad device.");
    }else{ // set the keyboard or keypad device to be the first keyboard device or, if no keyboard, the first keypad
        i=0;
        if(numDeviceIndices==0)
            PsychErrorExitMsg(PsychError_user, "No keyboard or keypad devices detected.");
        else{
            deviceIndex=deviceIndices[i];
        }
    }
    deviceRecord=deviceRecords[i]; 
#ifndef __LP64__
    usage = deviceRecord->usage;
#else
    usage = IOHIDDevice_GetUsage(deviceRecord);
#endif

	KeysUsagePage = ((usage == kHIDUsage_GD_Keyboard) || (usage == kHIDUsage_GD_Keypad)) ? kHIDPage_KeyboardOrKeypad : kHIDPage_Button;

    //Allocate and init out return arguments.
    PsychAllocOutDoubleArg(1, FALSE, &timeValueOutput);
	if(!timeValueOutput)
		PsychErrorExitMsg(PsychError_system, "Failed to allocate memory for output.");

	interface = PsychHIDGetDeviceInterfacePtrFromIndex(deviceIndex);
	if(!interface)
		PsychErrorExitMsg(PsychError_system, "Could not get interface to device.");
	
	// The following bracketed clause will get a cookie corresponding to the
	// trigger. If multiple keys were of interest, the code could be modified
	// trivially to iterate over an array of KeysUsage to generate an array of 
	// corresponding cookies
	{
		CFArrayRef elements=NULL;
		psych_bool usedDictionary=FALSE;
		{
			CFDictionaryRef dict=NULL;
		
			// The next few lines define a dictionary describing the key of interest
			// If they are omitted, the code will still work, but all elements will match
			// initially rather than just the one key of interest, so the code will be
			// slower since it will iterate through hundreds of keys
			CFStringRef keys[2] = {CFSTR("UsagePage"), CFSTR("Usage")};
			CFNumberRef values[2];
			values[0] = CFNumberCreate(kCFAllocatorDefault, kCFNumberIntType, &KeysUsagePage);
			values[1] = CFNumberCreate(kCFAllocatorDefault, kCFNumberIntType, &KeysUsage);
			if(values[0]!=NULL && values[1]!=NULL){
				// Even if they are NULL, it will be ok since dict can be NULL at the expense of some loss of efficiency
				dict = CFDictionaryCreate(kCFAllocatorDefault, (const void**)keys, (const void**)values, 2, &kCFCopyStringDictionaryKeyCallBacks, &kCFTypeDictionaryValueCallBacks);
			}
	
			// copyMatchinfElements requires IOHIDDeviceInterface122, thus Mac OS X 10.3 or higher
			// elements would have to be obtained directly from IORegistry for 10.2 or earlier
			// If dict==NULL, all elements will match, leading to some inefficiency
			success = (*interface)->copyMatchingElements(interface, dict, &elements);
		
			if(dict){
				usedDictionary=TRUE;
				CFRelease(dict);
			}

			if(values[0]) CFRelease(values[0]);
			if(values[1]) CFRelease(values[1]);
			
			if(!elements){
				PsychErrorExitMsg(PsychError_user, "Specified key code not found on device (I).");
			}
		}
		{
			// elements will only contain one element in this implementation, but has the
			// advantage of generalizing to future derived implementations that listen
			// for multiple keys
			CFIndex i;
			for (i=0; i<CFArrayGetCount(elements); i++)
			{
				long number;
				CFDictionaryRef element= CFArrayGetValueAtIndex(elements, i);
				CFTypeRef object;
				
				if(!element) continue;
				
				if(!usedDictionary){
				
					// Verify tht we are dealing with a keypad or keyboard
					object = CFDictionaryGetValue(element, CFSTR(kIOHIDElementUsageKey));
					if (object == 0 || CFGetTypeID(object) != CFNumberGetTypeID()) continue;
					if (!CFNumberGetValue((CFNumberRef) object, kCFNumberLongType,&number)) continue;
					if(number!=KeysUsage) continue;
				
					// See if element corresponds to the desired key
					object = CFDictionaryGetValue(element, CFSTR(kIOHIDElementUsagePageKey));
					if (object == 0 || CFGetTypeID(object) != CFNumberGetTypeID()) continue;
					if (!CFNumberGetValue((CFNumberRef) object, kCFNumberLongType, &number)) continue;
					if(number!=KeysUsagePage) continue;
				}
				
				// Get the cookie for this element
				object= (CFDictionaryGetValue(element, CFSTR(kIOHIDElementCookieKey)));
				if (object == 0 || CFGetTypeID(object) != CFNumberGetTypeID()) continue;
				if(!CFNumberGetValue((CFNumberRef) object, kCFNumberLongType, &number)) continue;
				triggerCookie = (IOHIDElementCookie) number;
				
				break;
			}
			if(CFArrayGetCount(elements)==i){
				CFRelease(elements);
				PsychErrorExitMsg(PsychError_user, "Specified key code not found on device (II).");
			}
			CFRelease(elements);
		}
	}

	// Allocate for the queue
	queue=(*interface)->allocQueue(interface);
	if(!queue){
		PsychErrorExitMsg(PsychError_system, "Failed to allocate event queue for detecting key press.");
	}
	
	// Create the queue
	result = (*queue)->create(queue, 0, 8);		// 8 events can be stored before the earliest event is lost
	if (kIOReturnSuccess != result){
		(*queue)->Release(queue);
		(*queue)=NULL;
		PsychErrorExitMsg(PsychError_system, "Failed to create event queue for detecting key press.");
	}
	
	// Add the trigger to the queue
	// If multiple keys were of interest, their cookies could be added in turn
	result = (*queue)->addElement(queue, triggerCookie, 0);
	if (kIOReturnSuccess != result){
		result = (*queue)->dispose(queue);
		(*queue)->Release(queue);
		(*queue)=NULL;
		PsychErrorExitMsg(PsychError_system, "Failed to add trigger key to event queues.");
	}
	
	// Start the queue
	result = (*queue)->start(queue);
	if (kIOReturnSuccess != result){
		result = (*queue)->dispose(queue);
		(*queue)->Release(queue);
		(*queue)=NULL;
		PsychErrorExitMsg(PsychError_system, "Failed to start event queues.");
	}
	 
	// Watch for the trigger
	{
		IOHIDEventStruct event;
		while(1){
			
			AbsoluteTime zeroTime = {0,0};
			result = (*queue)->getNextEvent(queue, &event, zeroTime, 0);
			if(kIOReturnSuccess==result) break;
			PsychWaitIntervalSeconds((double)0.004);  //surrender some time to other processes	
		
			// If it were of interest to trigger selectively on key press or key release,
			// this could be evaluated by checking event.value (zero versus non-zero)
			// but this would put more code inside the loop
			
			// If multiple keys are registered via addElement (not the case here), the
			// cookie for the key responsible for the event can be obtained from 
			// event.elementCookie
		}
		
		// If event.longValue is not NULL, the documentation indicates that it is up
		// to the caller to deallocate it. The circumstances under which a non-NULL
		// value would be generated are not specified. My guess is that some devices 
		// might return a 64-bit value (e.g., a tracking device coordinate).
		// Keys, having only two states, shouldn't need this, but check and free to
		// be safe		
		if ((event.longValueSize != 0) && (event.longValue != NULL)) free(event.longValue);
		
		// Set the time, using the same strategy as PsychTimeGlue's PsychGetPrecisionTimerSeconds
		// For code maintainability, it would be better if this conversion were performed
		// by a function in PsychTimeGlue
		{
			Nanoseconds timeNanoseconds=AbsoluteToNanoseconds(event.timestamp);
			UInt64 timeUInt64=UnsignedWideToUInt64(timeNanoseconds);
			double timeDouble=(double)timeUInt64;
			*timeValueOutput=timeDouble / 1000000000;
		}
	}
	
	// Clean up
	result = (*queue)->stop(queue);
	// Code from Apple sometimes removes elements from queue before disposing and sometimes does not
	// I can't see any reason to do so for a queue that's one line of code away from destruction
	// and therefore haven't
	result = (*queue)->dispose(queue);
	(*queue)->Release(queue);
	(*queue)=NULL;				// Just in case queue is redefined as static in the future
	
    return;
}

#endif
