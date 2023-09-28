/*
 
	Psychtoolbox-3/PsychSourceGL/Source/Common/Eyelink/EyelinkEDKVersion.c
 
	PROJECTS: Eyelink 
 
	AUTHORS:
    brian@sr-research.com   br
 
	HISTORY:
 
    03/06/21	fbr		created
 
	TARGET LOCATION:

    Eyelink.mex* resides in:
    PsychHardware/EyelinkToolbox/EyelinkBasic
 
 */

#include "PsychEyelink.h"

static char useString[] = "[edkVersionString]  = Eyelink('EDKVersion')";

static char synopsisString[] =
"Returns the eyelink_core library version number.\n";

static char seeAlsoString[] = "";

/*
 ROUTINE: EyelinkDllVersion
 PURPOSE:
	Returns the eyelink_core library version number.
	returns: 0 if no reply, else time
 */

PsychError EyelinkEDKVersion(void)
{
	char strEdkVersion[40]="0.0.0.0";
	
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(0));
	PsychErrorExit(PsychRequireNumInputArgs(0));
	PsychErrorExit(PsychCapNumOutputArgs(1));
	
//	// Verify eyelink is up and running
//	EyelinkSystemIsConnected();
//	EyelinkSystemIsInitialized();
	
	eyelink_dll_version(strEdkVersion);
	
	PsychCopyOutCharArg(1, FALSE, strEdkVersion);
	
	return(PsychError_none);	
}
