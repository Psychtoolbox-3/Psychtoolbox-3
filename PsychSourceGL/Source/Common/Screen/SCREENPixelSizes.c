/*
    Psychtoolbox3/Source/Common/SCREENPixelSizes.c

    AUTHORS:

    Allen.Ingling@nyu.edu		awi

    PLATFORMS:

    This file should build on any platform.

    HISTORY:

    12/03/02  awi		Created.

    DESCRIPTION:

    Return a list of the valid pixel sizes for the specified screen.

    TO DO:
 */

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "pixelSizes=Screen('PixelSizes', windowPtrOrScreenNumber)";
static char synopsisString[] =
"List all possible pixelSizes for this screen available at the current height and width settings. "
"pixelSizes refer to the color depth of single pixels in bits - also known as bpp or bits per pixel. "
"If you change height and with settings the avaialable depths might change.";
static char seeAlsoString[] = "Resolutions Resolution";

PsychError SCREENPixelSizes(void)
{
    int screenNumber;
    PsychDepthType depths;
    
    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));
    
    //get specified screen number.  This also sanity checks against the number of displays.
    PsychCopyInScreenNumberArg(kPsychUseDefaultArgPosition, TRUE, &screenNumber);
    
    //get the screen depths
    PsychInitDepthStruct(&depths);
    PsychGetScreenDepths(screenNumber, &depths);
    
    //Allocate a return matrix and load it with the depth values.
    PsychCopyOutDepthArg(1, FALSE, &depths);
    
    return(PsychError_none);
}
