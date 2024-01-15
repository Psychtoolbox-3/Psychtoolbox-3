/*

	/osxptb/trunk/PsychSourceGL/Source/OSX/Eyelink/EyelinkShutdown.c
  
	PROJECTS: Eyelink 
  
	AUTHORS:
		cburns@berkeley.edu				cdb
		E.Peters@ai.rug.nl				emp
		f.w.cornelissen@med.rug.nl		fwc
  
	PLATFORMS:	Currently only OS X  
    
	HISTORY:

		11/23/05  cdb		Created.

	TARGET LOCATION:

		Eyelink.mexmac resides in:
			EyelinkToolbox
*/

#include "PsychEyelink.h"
extern char eyelinkDisplayCallbackFunc[];

PsychError PsychEyelinkShutdown(void)
{
	int		iStatus		= -1;
	char	strMsg[256];

	if (giSystemInitialized) {
		// Zero-out return string:
		memset(strMsg, 0, sizeof(strMsg));
		
		// Disconnect if connected
		if (eyelink_is_connected()) {
			set_offline_mode();
			iStatus = eyelink_close(1);
			if (iStatus != 0) {
				sprintf(strMsg, "Eyelink shutdown: eyelink_close error %d\n", iStatus);
				PsychErrorExitMsg(PsychError_internal, strMsg);
			}
		}
		
		// Detach all callback hook functions:
		PsychEyelink_uninit_core_graphics();
		
		// If using default callback handler PsychEyelinkDispatchCallback, clear function
		// from persistence in memory
		if (!strcmp(eyelinkDisplayCallbackFunc, "PsychEyelinkDispatchCallback")) {
			PsychRuntimeEvaluateString("clear PsychEyelinkDispatchCallback");
		} else {
			mexPrintf("DEBUG: eyelinkDisplayCallbackFunc not match PsychEyelinkDispatchCallback\n");
		}

		// Close down eyelink and reset global flag
		close_eyelink_system();
		msec_delay(100);
		giSystemInitialized = 0;
	}
	
	return(PsychError_none);
}

PsychError EyelinkShutdown(void)
{
	static char useString[] = "Eyelink('Shutdown')";
	static char synopsisString[] = "Shuts down the eyelink";
	static char seeAlsoString[] = "";

	// Add help strings
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	
	// Output help if asked
	if(PsychIsGiveHelp()) {
		PsychGiveHelp();
		return(PsychError_none);
	}

	// Check arguments
	PsychErrorExit(PsychCapNumInputArgs(0));
	PsychErrorExit(PsychRequireNumInputArgs(0));
	PsychErrorExit(PsychCapNumOutputArgs(0));
	
	// Call actual shutdown routine:
	PsychEyelinkShutdown();

	return(PsychError_none);
}
