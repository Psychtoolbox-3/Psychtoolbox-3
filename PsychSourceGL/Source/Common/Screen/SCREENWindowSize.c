/*
  SCREENWindowSize.c		
  
  AUTHORS:
  
		Allen.Ingling@nyu.edu           awi
        mario.kleiner@tuebingen.mpg.de  mk

  PLATFORMS:
  
		All

  HISTORY:
  
		2/26/05  awi		Created.  Inspired by Frans Cornelissen's script "WindowSize".   

*/

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "[width, height]=Screen('WindowSize', windowPointerOrScreenNumber [, realFBSize=0]);";
//                                                              1                              2
static char synopsisString[] = 
	"Return the width and height of a window or screen in units of pixels.\n"
    "By default, the useable size in pixels is returned, ie., the size of the "
    "area usercode can draw to. If the optional 'realFBSize' flag is set to 1, "
    "the real size of the windows framebuffer is returned. Those sizes can differ, "
    "e.g., because certain stereo display modes or high color precision display "
    "modes require adjustments. The 'realFBSize' 1 setting is mostly for Psychtoolbox "
    "internal use, not for regular user-code.\n";

static char seeAlsoString[] = "Screen('Rect')";	

PsychError SCREENWindowSize(void)  
{
	PsychWindowRecordType *windowRecord;
	int screenNumber;
	PsychRectType rect;
	double	rectWidth, rectHeight;
    int realFBSize = 0;
    
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(2));		//The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1));	//Insist that the argument be present.   
	PsychErrorExit(PsychCapNumOutputArgs(2));		//The maximum number of outputs

    // Get optional 'realFBSize' flag: Defaults to zero.
    PsychCopyInIntegerArg(2, FALSE, &realFBSize);

	if(PsychIsScreenNumberArg(1)){
		PsychCopyInScreenNumberArg(1, TRUE, &screenNumber);
		PsychGetScreenRect(screenNumber, rect);
		rectWidth=PsychGetWidthFromRect(rect);
		rectHeight=PsychGetHeightFromRect(rect);
		PsychCopyOutDoubleArg(1, kPsychArgOptional, rectWidth);
		PsychCopyOutDoubleArg(2, kPsychArgOptional, rectHeight);
	}else if(PsychIsWindowIndexArg(1)){
		PsychAllocInWindowRecordArg(1, TRUE, &windowRecord);
		PsychOSProcessEvents(windowRecord, 0);

		rectWidth=PsychGetWidthFromRect((realFBSize) ? windowRecord->rect : windowRecord->clientrect);
		rectHeight=PsychGetHeightFromRect((realFBSize) ? windowRecord->rect : windowRecord->clientrect);

		PsychCopyOutDoubleArg(1, kPsychArgOptional, rectWidth);
		PsychCopyOutDoubleArg(2, kPsychArgOptional, rectHeight);
	}else
		PsychErrorExitMsg(PsychError_user, "Argument was recognized as neither a window index nor a screen pointer");

	return(PsychError_none);
}

PsychError SCREENDisplaySize(void)  
{
    static char useString[] = "[width, height]=Screen('DisplaySize', ScreenNumber);";
    static char synopsisString[] = 
	"Return the physical width and height of the output display device associated with 'ScreenNumber'. "
        "The size is returned in units of millimeters as reported by the display device itself to the OS. "
        "On MacOS-X, if Extended Display Identification Data (EDID) for the display device is not "
        "available, the size is estimated by OS-X based on the device width and height in pixels, with an "
        "assumed resolution of 2.835 pixels/mm or 72 DPI, a reasonable guess for displays predating "
        "EDID support. On M$-Windows and GNU/Linux, the behaviour in case of missing EDID data is "
        "unknown. This function returns a width and height of zero if physical display size can't be "
        "queried from the operating system. Please handle the returned information with great caution. "
        "It is unclear how accurate the EDID data for typical monitors or video beamers really is. This "
        "information could be pretty unreliable and therefore misleading!";
    static char seeAlsoString[] = "";	    
    int screenNumber, Width, Height;
    
    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    //check for superfluous arguments
    PsychErrorExit(PsychCapNumInputArgs(1));		//The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1));	//Insist that the argument be present.   
    PsychErrorExit(PsychCapNumOutputArgs(2));		//The maximum number of outputs
    
    // Retrieve screen number:
    PsychCopyInScreenNumberArg(1, TRUE, &screenNumber);
    // Query physical size of this screen in millimeters:
    PsychGetDisplaySize(screenNumber, &Width, &Height);
    // Return it:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) Width);
    PsychCopyOutDoubleArg(2, kPsychArgOptional, (double) Height);
    
    return(PsychError_none);
}
