/*
	PsychtoolboxGL/Source/Common/PsychHID/PsychHIDWait.c		
  
	PROJECTS: 
	
		PsychHID only.
  
	PLATFORMS:  
	
		Only OS X for now.  
  
	AUTHORS:
	
		Allen.Ingling@nyu.edu		awi 
		rwoods@ucla.edu				rpw
      
	HISTORY:
		6/5/03  awi		Created.
		4/7/05	awi		Added a wait to to the polling loop. 
		12/17/09 rpw	Added support for keypads
  
	NOTES:
	
		We add a compile flag to the project MachTimebase settings: "-Wno-long-double" turns off a warning caused by using type long double.  
		The warning is:
			use of `long double' type; its size may change in a future release (Long double usage is reported only once for each file.
			To disable this warning, use -Wno-long-double.)
			
		Turning off the warning is probably ok because we do not depend on long doubles being any particular size or precision; if does change,
		that's not going to break anything.  
			
		The -Wno-long-double flag is appended to the project setting "Other Warning Flags".

		---
    

*/

#include "PsychHID.h"
#if PSYCH_SYSTEM == PSYCH_OSX

#define NUMDEVICEUSAGES 7

static char useString[]= "secs=PsychHID('KbWait', [deviceNumber])";
static char synopsisString[] = 
		"THIS FUNCTION IS DEPRECATED AND SHOULD NOT BE USED! USE KbWait INSTEAD!\n\n"
        "Scan a keyboard or keypad device and return wait for a keypress "
        "By default the first keyboard device (the one with the lowest device number) is "
        "scanned. If no keyboard is found, the first keypad device is "
        "scanned.  Optionally, the device number of any keyboard or keypad may be specified "
        "The PsychHID('KbWait') implements the KbCheck command a provided by the  OS 9 Psychtoolbox."
        "KbWait is defined in the OS X Psychtoolbox and invokes PsychHID('KbWait').  For backwards "
        "compatability with earlier Psychtoolbox we recommend using KbWait instead of PsychHID('KbWait'). ";
        
static char seeAlsoString[] = "";

 
PsychError PSYCHHIDKbWait(void) 
{
    pRecDevice          deviceRecord;
    pRecElement			currentElement;
    int					i, deviceIndex, numDeviceIndices;
    int 				numDeviceUsages=NUMDEVICEUSAGES;
	long				KbDeviceUsagePages[NUMDEVICEUSAGES]= {kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop, kHIDPage_GenericDesktop};
	long				KbDeviceUsages[NUMDEVICEUSAGES]={kHIDUsage_GD_Keyboard, kHIDUsage_GD_Keypad, kHIDUsage_GD_Mouse, kHIDUsage_GD_Pointer, kHIDUsage_GD_Joystick, kHIDUsage_GD_GamePad, kHIDUsage_GD_MultiAxisController};
    int					deviceIndices[PSYCH_HID_MAX_KEYBOARD_DEVICES]; 
    pRecDevice			deviceRecords[PSYCH_HID_MAX_KEYBOARD_DEVICES];
    psych_bool 			isDeviceSpecified, foundUserSpecifiedDevice;
    double				*timeValueOutput;
    psych_bool 			isKeyDown;
    uint32_t            usage, usagePage;

    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));  	//Specifies the number of the keyboard or keypad to scan.  

    // We don't support this deprecated function on 64-Bit OSX anymore:
#ifndef __LP64__
    
    PsychHIDVerifyInit();
    
    //Choose the device index and its record
    PsychHIDGetDeviceListByUsages(numDeviceUsages, KbDeviceUsagePages, KbDeviceUsages, &numDeviceIndices, deviceIndices, deviceRecords);  
    isDeviceSpecified=PsychCopyInIntegerArg(1, FALSE, &deviceIndex);
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
    
    //Allocate and init out return arguments.
    PsychAllocOutDoubleArg(1, FALSE, &timeValueOutput);
          
    //step through the elements of the device.  Set flags in the return array for down keys.
    isKeyDown=FALSE;
    while(!isKeyDown){
		PsychWaitIntervalSeconds((double)0.004);  //surrender some time to other processes	
        for(currentElement=HIDGetFirstDeviceElement(deviceRecord, kHIDElementTypeInput); 
            currentElement != NULL; 
            currentElement=HIDGetNextDeviceElement(currentElement, kHIDElementTypeInput))
        {
            usage = currentElement->usage;
            usagePage = currentElement->usagePage;

            if(((usagePage == kHIDPage_KeyboardOrKeypad) || (usagePage == kHIDPage_Button)) && (usage <= 256) && (usage >= 1)){
                if (HIDGetElementValue(deviceRecord, currentElement)) {
                    isKeyDown=TRUE;
                    break;  //break out of inner for loop.  
                }
            }
        }
    }
    PsychGetPrecisionTimerSeconds(timeValueOutput);

#else
    PsychErrorExitMsg(PsychError_unimplemented, "This function is not supported by the 64-Bit version of PsychHID.");
#endif
    return(PsychError_none);	
}

#endif
