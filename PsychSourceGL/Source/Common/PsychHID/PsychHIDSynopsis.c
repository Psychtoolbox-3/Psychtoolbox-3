/*
  PsychToolbox2/Source/Common/PsychHID/PsychHIDSynopsis.c	
  
  AUTHORS:
  Allen.Ingling@nyu.edu             awi 
  mario.kleiner@tuebingen.mpg.de    mk
  
  PLATFORMS: All
  
  PROJECTS: PsychHID
   

  HISTORY:
  06/20/03  awi		Wrote it.  Based on ScreenSynopsis  
  
*/

#include "PsychHID.h"

void InitializeSynopsis(void); // I added this prototype to make the compiler happy. dgp.

#define MAX_SYNOPSIS_STRINGS 500  

//declare variables local to this file.  
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

void InitializeSynopsis(void)
{
	int i=0;
	const char **synopsis = synopsisSYNOPSIS;  // abbreviate the long name

	synopsis[i++] = "Usage:";

	synopsis[i++] = "\n\nGeneral commands:\n\n";
    synopsis[i++] = "rc = PsychHID('KeyboardHelper', commandCode)";
    
	synopsis[i++] = "\n\nSupport for generic USB-HID devices:\n\n";
	synopsis[i++] = "numberOfDevices=PsychHID('NumDevices')";
	synopsis[i++] = "numberOfElements=PsychHID('NumElements',deviceNumber)";
	synopsis[i++] = "numberOfCollections=PsychHID('NumCollections',deviceNumber)";
	synopsis[i++] = "devices=PsychHID('Devices' [, deviceClass])";
	synopsis[i++] = "elements=PsychHID('Elements',deviceNumber)";
	synopsis[i++] = "collections=PsychHID('Collections',deviceNumber)";
	synopsis[i++] = "elementState=PsychHID('RawState',deviceNumber,elementNumber)";
	synopsis[i++] = "elementState=PsychHID('CalibratedState',deviceNumber,elementNumber)";
	synopsis[i++] = "[keyIsDown,secs,keyCode]=PsychHID('KbCheck' [, deviceNumber][, scanList])";
	synopsis[i++] = "[report,err]=PsychHID('GetReport',deviceNumber,reportType,reportID,reportBytes)";
	synopsis[i++] = "err=PsychHID('SetReport',deviceNumber,reportType,reportID,report)";
	synopsis[i++] = "[reports,err]=PsychHID('GiveMeReports',deviceNumber,[reportBytes])";
	synopsis[i++] = "err=PsychHID('ReceiveReports',deviceNumber[,options])";
	synopsis[i++] = "err=PsychHID('ReceiveReportsStop',deviceNumber)";
    
	synopsis[i++] = "\n\nQueue based keyboard queries: See 'help KbQueueCreate' for explanations:\n\n";
	synopsis[i++] = "PsychHID('KbQueueCreate', [deviceNumber], [keyFlags])";
	synopsis[i++] = "PsychHID('KbQueueRelease' [, deviceIndex])"; 
	synopsis[i++] = "[navail] = PsychHID('KbQueueFlush' [, deviceIndex][, flushType=1])"; 
	synopsis[i++] = "PsychHID('KbQueueStart' [, deviceIndex])"; 
	synopsis[i++] = "PsychHID('KbQueueStop' [, deviceIndex])"; 
	synopsis[i++] = "[keyIsDown, firstKeyPressTimes, firstKeyReleaseTimes, lastKeyPressTimes, lastKeyReleaseTimes]=PsychHID('KbQueueCheck' [, deviceIndex])"; 
	synopsis[i++] = "secs=PsychHID('KbTriggerWait', KeysUsage, [deviceNumber])";
	synopsis[i++] = "[event, navail] = PsychHID('KbQueueGetEvent' [, deviceIndex][, maxWaitTimeSecs=0])";

	synopsis[i++] = "\n\nSupport for access to generic USB devices: See 'help ColorCal2' for one usage example:\n\n";
	synopsis[i++] = "usbHandle = PsychHID('OpenUSBDevice', vendorID, deviceID [, configurationId=0])";
	synopsis[i++] = "PsychHID('CloseUSBDevice' [, usbHandle])";
	synopsis[i++] = "outData = PsychHID('USBControlTransfer', usbHandle, bmRequestType, bRequest, wValue, wIndex, wLength, inData)";

	synopsis[i++] = NULL;  // this tells PsychDisplayPsychHIDSynopsis where to stop
	if (i > MAX_SYNOPSIS_STRINGS) {
		PrintfExit("%s: increase dimension of synopsis[] from %ld to at least %ld and recompile.",__FILE__,(long)MAX_SYNOPSIS_STRINGS,(long)i);
	}
}

PsychError PsychDisplayPsychHIDSynopsis(void)
{
	int i;	
	for (i = 0; synopsisSYNOPSIS[i] != NULL; i++) printf("%s\n",synopsisSYNOPSIS[i]);

	return(PsychError_none);
}
