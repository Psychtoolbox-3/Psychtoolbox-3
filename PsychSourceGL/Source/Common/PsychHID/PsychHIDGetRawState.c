/*
  PsychToolbox3/Source/Common/PsychHID/PsychHIDGetRawState.c		

  PROJECTS: PsychHID only.

  PLATFORMS:  All.
  
  AUTHORS:

  Allen.Ingling@nyu.edu		  awi
  mario.kleiner@tuebingen.mpg.de  mk

  HISTORY:
  5/12/03  awi		Created.
  7/29/11  mk           Refactored for multi-platform.  
  TO DO:

*/

#include "PsychHID.h"

static char useString[]= "elementState=PsychHID('RawState', deviceNumber, elementNumber)";
static char synopsisString[] =
        "Return the immediate state of the specified element on the specified device.";
static char seeAlsoString[] = "";

PsychError PSYCHHIDGetRawState(void) 
{
    int			deviceIndex, elementIndex;
    double              value;

    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));

    PsychCopyInIntegerArg(1, TRUE, &deviceIndex);
    PsychCopyInIntegerArg(2, TRUE, &elementIndex);

    PsychHIDOSGamePadAxisQuery(deviceIndex, elementIndex, NULL, NULL, &value, NULL);

    if (PSYCH_SYSTEM != PSYCH_LINUX) PsychCopyOutDoubleArg(1, FALSE, value);

    return(PsychError_none);
}

#if PSYCH_SYSTEM == PSYCH_OSX

PsychError PsychHIDOSGamePadAxisQuery(int deviceIndex, int axisId, double* min, double* max, double* val, char* axisLabel)
{
    double          elementValue;
    pRecDevice      deviceRecord;
    pRecElement		elementRecord;
    int             elementIndex = axisId;

    PsychHIDVerifyInit();
    deviceRecord= PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
    elementRecord= PsychHIDGetElementRecordFromDeviceRecordAndElementIndex(deviceRecord, elementIndex);
#ifndef __LP64__
    elementValue = (double) HIDGetElementValue(deviceRecord, elementRecord);
#else
    elementValue = IOHIDElement_GetValue(elementRecord, kIOHIDValueScaleTypePhysical);
#endif
    if (val) *val = elementValue;
    return(PsychError_none);
}

#endif
