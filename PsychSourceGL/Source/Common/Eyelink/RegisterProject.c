/*

	/osxptb/trunk/PsychSourceGL/Source/Common/Eyelink/RegisterProject.c

	PROJECTS: Eyelink 
  
	AUTHORS:
		cburns@berkeley.edu				cdb
		E.Peters@ai.rug.nl				emp
		f.w.cornelissen@med.rug.nl		fwc
		e_flister@yahoo.com				edf
  
	PLATFORMS:	All.
    
	HISTORY:
		2003	  emp		created alpha version
		11/21/05  cdb		Adapted from alpha version.
		15/06/06  fwc		adapted and added functions
		19/10/06  fwc		added raw data function
		21/01/07  fwc		added new timing functions
		19/02/09  edf		added EyelinkGetFloatDataRaw
		22/03/09  edf		added EyelinkGetQueuedData

	TARGET LOCATION:

		Eyelink.mexmac resides in:
			EyelinkToolbox
*/

#include "RegisterProject.h"
#include "Psych.h"

PsychError PsychModuleInit(void)
{
	// Initialize Eyelink project variables
	giSystemInitialized = 0;
	
	// Register the project exit function:
	PsychErrorExit(PsychRegisterExit(&PsychEyelinkShutdown)); 
	
	// Register the project function which is called when the module
	// is invoked with no arguments:
	PsychErrorExitMsg(PsychRegister(NULL,  &PsychDisplayEyelinkSynopsis), "Failed to register the Eyelink synopsis function.");        
	
	// Register the module name:
	PsychErrorExitMsg(PsychRegister("Eyelink", NULL), "Failed to register Eyelink Module name!");

    // Report module version:
	PsychErrorExit(PsychRegister("Version",  &MODULEVersion));

	// Register sub-commands
	PsychErrorExit(PsychRegister("ButtonStates",		&EyelinkButtonStates));
	PsychErrorExit(PsychRegister("Command",				&EyelinkCommand));
	PsychErrorExit(PsychRegister("CurrentMode",			&EyelinkCurrentMode));
	PsychErrorExit(PsychRegister("ImageModeDisplay",	&EyelinkImageModeDisplay));
	PsychErrorExit(PsychRegister("InitializeDummy",		&EyelinkInitializeDummy));

	PsychErrorExit(PsychRegister("Initialize",			&EyelinkInitialize));
	PsychErrorExit(PsychRegister("IsConnected",			&EyelinkIsConnected));
	PsychErrorExit(PsychRegister("SetAddress",			&EyelinkSetAddress));
	PsychErrorExit(PsychRegister("NewestFloatSample",	&EyelinkNewestFloatSample));
	PsychErrorExit(PsychRegister("SendKeyButton",		&EyelinkSendKeyButton));
	PsychErrorExit(PsychRegister("Shutdown",			&EyelinkShutdown));
	PsychErrorExit(PsychRegister("StartRecording",		&EyelinkStartRecording));
	PsychErrorExit(PsychRegister("StartSetup",			&EyelinkStartSetup));
	PsychErrorExit(PsychRegister("StopRecording",		&EyelinkStopRecording));
	//added as of 15/06/06
	PsychErrorExit(PsychRegister("TargetCheck",				&EyelinkTargetCheck));
	PsychErrorExit(PsychRegister("AcceptTrigger",			&EyelinkAcceptTrigger));
	PsychErrorExit(PsychRegister("WaitForModeReady",		&EyelinkWaitForModeReady));
	PsychErrorExit(PsychRegister("CalResult",				&EyelinkCalResult));
	PsychErrorExit(PsychRegister("Message",					&EyelinkMessage));
	PsychErrorExit(PsychRegister("OpenFile",				&EyelinkOpenFile));
	PsychErrorExit(PsychRegister("CloseFile",				&EyelinkCloseFile));
	PsychErrorExit(PsychRegister("ReceiveFile",				&EyelinkReceiveFile));
	PsychErrorExit(PsychRegister("CheckRecording",			&EyelinkCheckRecording));
	PsychErrorExit(PsychRegister("NewFloatSampleAvailable",	&EyelinkNewFloatSampleAvailable));
	PsychErrorExit(PsychRegister("SetOfflineMode",			&EyelinkSetOfflineMode));
	PsychErrorExit(PsychRegister("DriftCorrStart",			&EyelinkDriftCorrStart));
	PsychErrorExit(PsychRegister("ApplyDriftCorr",			&EyelinkApplyDriftCorr));
	PsychErrorExit(PsychRegister("EyeAvailable",			&EyelinkEyeAvailable));
	PsychErrorExit(PsychRegister("ReadTime",				&EyelinkReadTime));
	PsychErrorExit(PsychRegister("RequestTime",				&EyelinkRequestTime));
	PsychErrorExit(PsychRegister("GetTrackerVersion",		&EyelinkGetTrackerVersion));

	//added as of 19/10/06
	PsychErrorExit(PsychRegister("NewestFloatSampleRaw",	&EyelinkNewestFloatSampleRaw));

	//added as of 30/10/06
	PsychErrorExit(PsychRegister("GetNextDataType",	&EyelinkGetNextDataType));
	PsychErrorExit(PsychRegister("GetNextData",	&EyelinkGetNextDataType));  // for compatibility reasons only
	PsychErrorExit(PsychRegister("GetFloatData",	&EyelinkGetFloatData));

	// added as of 21/01/07
	PsychErrorExit(PsychRegister("TrackerTime",				&EyelinkTrackerTime));
	PsychErrorExit(PsychRegister("TimeOffset",				&EyelinkTimeOffset));

	// added as of 19/02/09
	PsychErrorExit(PsychRegister("GetFloatDataRaw",	&EyelinkGetFloatDataRaw));    
    
	// added as of 22/03/09
	PsychErrorExit(PsychRegister("GetQueuedData", &EyelinkGetQueuedData));
	PsychErrorExit(PsychRegister("Verbosity", &EyelinkVerbosity));
	PsychErrorExit(PsychRegister("TestSuite", &EyelinkTestSuite));
	
	// NJ: Added as of 24/05/11
	PsychErrorExit(PsychRegister("ImageTransfer", &EyelinkImageTransfer));
	PsychErrorExit(PsychRegister("TrackerMode",	&EyelinkTrackerMode));
	PsychErrorExit(PsychRegister("CalMessage", &EyelinkCalMessage));
	PsychErrorExit(PsychRegister("ReadFromTracker", &EyelinkReadFromTracker));

	//register synopsis and named subfunctions.
	InitializeSynopsis();   //Scripting glue won't require this if the function takes no arguments.
	PsychSetModuleAuthorByInitials("emp");
	PsychSetModuleAuthorByInitials("fwc");
	PsychSetModuleAuthorByInitials("cdb");
	PsychSetModuleAuthorByInitials("mk");
	PsychSetModuleAuthorByInitials("edf");
	PsychSetModuleAuthorByInitials("nj");

	return(PsychError_none);
}
