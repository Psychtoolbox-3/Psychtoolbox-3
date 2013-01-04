/*
  Psychtoolbox3/Source/Common/SCREENPixelSize.c		
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
  
  PLATFORMS:	
  This file should build on any platform. 

  HISTORY:
  1/14/03  awi		Created. 
 
  DESCRIPTION:
  
  Return the current pixel size of a screen or window.

  TO DO:

*/

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "pixelSize=Screen('PixelSize', windowPtrOrScreenNumber)";
static char synopsisString[] =
    "Return the current pixel size of a screen or window. "
    "pixel size refers to the color depth of a pixel in bits - also known as bpp or bits per pixel, "
    "not to spatial size of the pixel. ";
static char seeAlsoString[] = "";

PsychError SCREENPixelSize(void)
{
    int 			screenNumber;
    double 			depth;
    PsychWindowRecordType	*windowRecord;

    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));

    //get specified screen number.
    if(PsychIsScreenNumberArg(1)){
        PsychCopyInScreenNumberArg(1, TRUE, &screenNumber);
        depth=(double)PsychGetScreenDepthValue(screenNumber);
    }else if(PsychIsWindowIndexArg(1)){
        PsychAllocInWindowRecordArg(1,TRUE,&windowRecord);
        depth=(double)windowRecord->depth;
    }else
        PsychErrorExit(PsychError_invalidNumdex);

    //Allocate a return matrix and load it with the depth values.
    PsychCopyOutDoubleArg(1, FALSE, depth);

    return(PsychError_none);
}
