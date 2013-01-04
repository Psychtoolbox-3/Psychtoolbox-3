/*
  SCREENRect.c		
  
  AUTHORS:
  
		Allen.Ingling@nyu.edu		awi 
  		mario.kleiner@tuebingen.mpg.de	mk
  
  PLATFORMS:
  
		All.

  HISTORY:
  
		1/14/03		awi		Created.   
		11/14/08	mk		Cleaned up.
  
  TO DO:

*/

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "rect=Screen('Rect', windowPointerOrScreenNumber [, realFBSize=0]);";
//                                             1                              2
static char synopsisString[] = 
	"Get local rect of window or screen. This has its top-left corner always at (0,0) "
	"and encodes the useable size of the window or screen. E.g., in certain stereo "
	"display modes or other special output modes, the actual useable window area for "
	"stimulus drawing may be much smaller than the real area occupied by the window. "
	"Example: In interleaved stereo modes, the effective useable height of a window "
	"is only half the real height of the window. Use this function to get the actual "
	"useable drawing area for a window or screen.\n"
    "If the optional 'realFBSize' flag is set to 1, then the function returns the "
    "real size of the windows framebuffer. This is mostly for Psychtoolbox internal "
    "use, not for regular user-code.\n";
static char seeAlsoString[] = "";

PsychError SCREENRect(void)  
{
	
	PsychWindowRecordType *windowRecord;
	int screenNumber;
	PsychRectType rect; 
    int realFBSize = 0;
    
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(2));		//The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1));	//Insist that the argument be present.   
	PsychErrorExit(PsychCapNumOutputArgs(1));		//The maximum number of outputs

    // Get optional 'realFBSize' flag: Defaults to zero.
    PsychCopyInIntegerArg(2, FALSE, &realFBSize);
    
	if(PsychIsScreenNumberArg(1)){
		PsychCopyInScreenNumberArg(1, TRUE, &screenNumber);
		PsychGetScreenRect(screenNumber, rect);
		PsychCopyOutRectArg(1, FALSE, rect);
	}else if(PsychIsWindowIndexArg(1)){
        PsychAllocInWindowRecordArg(1, TRUE, &windowRecord);
        PsychOSProcessEvents(windowRecord, 0);		
        PsychCopyOutRectArg(1, FALSE, (realFBSize) ? windowRecord->rect : windowRecord->clientrect);
	}else
		PsychErrorExitMsg(PsychError_user, "Argument was recognized as neither a window index nor a screen pointer");
    
	return(PsychError_none);
}
