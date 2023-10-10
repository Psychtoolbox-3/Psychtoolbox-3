/*
    Psychtoolbox-3/PsychSourceGL/Source/Common/Eyelink/EyelinkEDKVersion.c

    PROJECTS: Eyelink 

    AUTHORS:
    brian@sr-research.com   br

*/

#include "PsychEyelink.h"

static char useString[] = "[edkVersionString] = Eyelink('EDKVersion')";

static char synopsisString[] = "Returns a string with the eyelink_core library version number.\n";

static char seeAlsoString[] = "";

/*
ROUTINE: EyelinkDllVersion
PURPOSE:
    Returns the eyelink_core library version number.
    returns: 0 if no reply, else time
*/

PsychError EyelinkEDKVersion(void)
{
    char strEdkVersion[40];

    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumInputArgs(0));
    PsychErrorExit(PsychRequireNumInputArgs(0));
    PsychErrorExit(PsychCapNumOutputArgs(1));

    eyelink_dll_version(strEdkVersion);

    PsychCopyOutCharArg(1, FALSE, strEdkVersion);

    return(PsychError_none);	
}
