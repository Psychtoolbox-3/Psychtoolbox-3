/*
  PsychToolbox3/Source/Common/PsychHID/PsychHIDGetCalibratedState.c		
  
  PROJECTS: PsychHID only.
  
  PLATFORMS:  Only OS X so far.  
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
      
  HISTORY:
  5/12/03  awi		Created.
  
  TO DO:

*/

#include "PsychHID.h"

#if PSYCH_SYSTEM == PSYCH_OSX

static char useString[]= "elementState=PsychHID('CalibratedState', deviceNumber, elementNumber)";
static char synopsisString[] = 
        "Return the calibrated immediate state of the specified element on the specified device. "
        "Elements states have both a nominal and an actual range.  For joysticks and other rotary elements "
        "the actual range of reported state values is smaller than than the nominal range.  The HID driver stores the "
        "most extreme values reported by the element in the "
        "calMin and calMax fields of the element structure returned by PsychHID('Elements'); "
        "The actual state as returned by 'RawState' is scalled using calMin and calMax and returned "
        "by 'CalibratedState'. Insure that calMin and calMax hold true max and min values by swinging "
        "the device element to its poles before calling 'CalibratedState'."; 
        
static char seeAlsoString[] = "";
 
PsychError PSYCHHIDGetCalibratedState(void) 
{
    long 		elementValue;
    SInt32		calibratedState;
    pRecDevice  deviceRecord;
    pRecElement	elementRecord;
    int			deviceIndex, elementIndex;

    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    	 
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));
    
    PsychCopyInIntegerArg(1, TRUE, &deviceIndex);
    PsychCopyInIntegerArg(2, TRUE, &elementIndex);
    PsychHIDVerifyInit();
    deviceRecord= PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
    elementRecord= PsychHIDGetElementRecordFromDeviceRecordAndElementIndex(deviceRecord, elementIndex);
#ifndef __LP64__
    elementValue= HIDGetElementValue(deviceRecord, elementRecord);
    calibratedState= HIDCalibrateValue((SInt32)elementValue, elementRecord);
    PsychCopyOutDoubleArg(1, FALSE, (double)calibratedState);
#else
    PsychCopyOutDoubleArg(1, FALSE, IOHIDElement_GetValue(elementRecord, kIOHIDValueScaleTypeCalibrated));
#endif
    return(PsychError_none);	
}

#endif
