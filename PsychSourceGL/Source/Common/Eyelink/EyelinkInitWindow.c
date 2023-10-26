/*

	/osxptb/trunk/PsychSourceGL/Source/OSX/Eyelink/EyelinkInitWindow.c
  
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
//#include "WindowBank.h"

static char useString[] = "Eyelink('InitWindow')";
	
static char synopsisString[] = 
   "initializes eyelink global window, port, device & pixel depth\n" 
   "and rect variables";

static char seeAlsoString[] = "";

/*
ROUTINE: EYELINKinitwindow
PURPOSE:
   initializes eyelink global window, port, device & pixeldepth
   and rect variables */ 

PsychError EyelinkInitWindow(void)
{
	// Add help strings
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	
	// Output help if asked
	if(PsychIsGiveHelp()) {
		PsychGiveHelp();
		return(PsychError_none);
	}

	// Abort early with simple warning message, until it's actually implemented.
	mexPrintf("WARNING: Eyelink initwindow not implemented yet!  Function call ignored!\n");
	return(PsychError_none);

	
	/*
	// Check arguments
	PsychErrorExit(PsychCapNumInputArgs(0));
	PsychErrorExit(PsychRequireNumInputArgs(0));
	PsychErrorExit(PsychCapNumOutputArgs(0));    

	// Verify eyelink is up and running
	EyelinkSystemIsConnected();
	EyelinkSystemIsInitialized();
	
	return(PsychError_none);
	*/
}

