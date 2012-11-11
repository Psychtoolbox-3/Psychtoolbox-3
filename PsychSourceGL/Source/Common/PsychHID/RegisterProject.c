/*
  PsychToolbox3/Source/Common/PsychHID/RegisterProject.c		
  
  PROJECTS: PsychHID
  
  PLATFORMS:  All
  
  AUTHORS:
  Allen.Ingling@nyu.edu             awi 
  mario.kleiner@tuebingen.mpg.de    mk
  
  HISTORY:
  4/16/03  awi		Created. 
  4/15/05  dgp		Added Get/SetReport. 
  8/23/07  rpw      Added PsychHIDKbQueue suite and PsychHIDKbTriggerWait

*/

#include "Psych.h"
#include "PsychHID.h"

PsychError PsychModuleInit(void)
{
	InitializeSynopsis();  //first initialize the list of synopsis strings.

	// Register the project function which is called when the module
	// is invoked with no arguments: In this case just print command overview.
	PsychRegister(NULL,  &PsychDisplayPsychHIDSynopsis);

	//register the module name
	PsychErrorExit(PsychRegister("PsychHID", NULL));

	//register the module exit function
	PsychRegisterExit(&PsychHIDCleanup);

	//register module subfunctions
	PsychErrorExit(PsychRegister("Version",  &MODULEVersion));
	PsychErrorExit(PsychRegister("NumDevices",  &PSYCHHIDGetNumDevices));
	PsychErrorExit(PsychRegister("Devices",  &PSYCHHIDGetDevices));

#if PSYCH_SYSTEM == PSYCH_OSX
    // OS/X only, so far:
	PsychErrorExit(PsychRegister("NumElements",  &PSYCHHIDGetNumElements));
	PsychErrorExit(PsychRegister("Elements",  &PSYCHHIDGetElements));
	PsychErrorExit(PsychRegister("NumCollections",  &PSYCHHIDGetNumCollections));
	PsychErrorExit(PsychRegister("Collections",  &PSYCHHIDGetCollections));
	PsychErrorExit(PsychRegister("CalibratedState",  &PSYCHHIDGetCalibratedState));

	PsychErrorExit(PsychRegister("KbWait",  &PSYCHHIDKbWait));
#endif

	PsychErrorExit(PsychRegister("KbTriggerWait", &PSYCHHIDKbTriggerWait));
	PsychErrorExit(PsychRegister("KbQueueCreate", &PSYCHHIDKbQueueCreate));
	PsychErrorExit(PsychRegister("KbQueueStart", &PSYCHHIDKbQueueStart));
	PsychErrorExit(PsychRegister("KbQueueStop", &PSYCHHIDKbQueueStop));
	PsychErrorExit(PsychRegister("KbQueueCheck", &PSYCHHIDKbQueueCheck));
	PsychErrorExit(PsychRegister("KbQueueFlush", &PSYCHHIDKbQueueFlush));
	PsychErrorExit(PsychRegister("KbQueueRelease", &PSYCHHIDKbQueueRelease));
	PsychErrorExit(PsychRegister("KbQueueGetEvent", &PSYCHHIDKbQueueGetEvent));

	PsychErrorExit(PsychRegister("RawState",  &PSYCHHIDGetRawState));
	PsychErrorExit(PsychRegister("KbCheck",  &PSYCHHIDKbCheck));

	PsychErrorExit(PsychRegister("GetReport",  &PSYCHHIDGetReport));
	PsychErrorExit(PsychRegister("ReceiveReports",  &PSYCHHIDReceiveReports));
	PsychErrorExit(PsychRegister("ReceiveReportsStop",  &PSYCHHIDReceiveReportsStop));
	PsychErrorExit(PsychRegister("GiveMeReports",  &PSYCHHIDGiveMeReports));
	PsychErrorExit(PsychRegister("SetReport",  &PSYCHHIDSetReport));
	PsychErrorExit(PsychRegister("OpenUSBDevice", &PSYCHHIDOpenUSBDevice));
	PsychErrorExit(PsychRegister("CloseUSBDevice", &PSYCHHIDCloseUSBDevice));
	PsychErrorExit(PsychRegister("USBControlTransfer", &PSYCHHIDUSBControlTransfer));
    PsychErrorExit(PsychRegister("KeyboardHelper", &PSYCHHIDKeyboardHelper));

	PsychSetModuleAuthorByInitials("awi");
	PsychSetModuleAuthorByInitials("dgp");
	PsychSetModuleAuthorByInitials("mk");
	PsychSetModuleAuthorByInitials("rpw");
	PsychSetModuleAuthorByInitials("cgb");

	// Initialize PsychHID core functions:
	PsychInitializePsychHID();

	return(PsychError_none);
}
