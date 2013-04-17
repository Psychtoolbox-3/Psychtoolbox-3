/*

	PsychSourceGL/Source/Common/Screen/SCREENClose.c	

	AUTHORS:

		Allen.Ingling@nyu.edu				awi 
		mario.kleiner at tuebingen.mpg.de	mk
	PLATFORMS:	

		All.  

	HISTORY:

		07/26/02  awi		Created from OpenWindow
		01/30/05  mk        Closes/Deletes all textures, if no windowOrTextureIndex provided. "CloseAllTextures" - Functionality.
		01/18/08  mk		Closes/Deletes all textures or offscreen windows if a vector of handles is passed.

	DESCRIPTION:
	
		Close named onscreen windows, offscreen windows, proxy windows or textures.
		Batch close all offscreen windows, proxies and textures with a single convenience call on request.

	TO DO:

*/

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('Close', [windowOrTextureIndex or list of textureIndices/offscreenWindowIndices]);";
static char synopsisString[] = 
	"Close an onscreen or offscreen window or a texture. If the optional windowOrTextureIndex isn't "
	"provided, then all textures and offscreen windows are closed/deleted while regular onscreen "
	"windows are left open. If you want to close a subset of your offscreen windows or textures, but "
	"not all of them you can also pass in a vector of texture/offscreen window handles and all handles "
	"in the vector will be closed. ";  
static char seeAlsoString[] = "OpenWindow, OpenOffscreenWindow";	 

PsychError SCREENClose(void)
{
	PsychWindowRecordType			*windowRecord;
	int								screenNumber;
	PsychWindowRecordType			**windowRecordArray;
	int								i, numWindows; 
	int								*winHandles;

	windowRecord=NULL;
	
	//all subfunctions should have these two lines.  
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

	PsychErrorExit(PsychCapNumInputArgs(1));      //The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(0));  //The minimum required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(0));     //The maximum number of outputs

	// First try to alloc in a whole list of handles:
	winHandles = NULL;
	numWindows = 0;
	PsychAllocInIntegerListArg(1, FALSE, &numWindows, &winHandles);

	// None, One or many handles?
	if (winHandles && (numWindows > 1)) {
		// Multiple window handles provided: Iterate over them and close them all:
		for(i=0; i < numWindows; i++) {
			// Iterate over all handles, ignore all but texture/offscreen window handles:
			if (IsWindowIndex(winHandles[i]) && (PsychError_none == FindWindowRecord(winHandles[i], &windowRecord)) &&
				(windowRecord->windowType == kPsychTexture)) {
				// Ok its a texture or offscreen window - We close it:
				PsychCloseWindow(windowRecord);
			}
		}

		return(PsychError_none);		
	}
	
	// Either none or exactly one window handle provided...
	
	// Try to get the window record or exit with an error if the windex was bogus.
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, kPsychArgOptional, &windowRecord);

	// Window handle of a specific window provided?
	if (windowRecord==NULL) {
		// No window handle provided: In this case, we close/destroy all textures:
		PsychCreateVolatileWindowRecordPointerList(&numWindows, &windowRecordArray);
		for(i=0;i<numWindows;i++) {
			if ((windowRecordArray[i]->windowType == kPsychTexture) && !(windowRecordArray[i]->specialflags & kPsychDontDeleteOnClose)) PsychCloseWindow(windowRecordArray[i]);			
		}

		PsychDestroyVolatileWindowRecordPointerList(windowRecordArray);
		return(PsychError_none);
	}

	// Window handle of a specific window or texture provided: Close it...
	if(PsychIsLastOnscreenWindow(windowRecord)){
		// Check for stale texture ressources and movies. Report to user if any:
		PsychRessourceCheckAndReminder(TRUE);	

		screenNumber=windowRecord->screenNumber;
		PsychCloseWindow(windowRecord);
		PsychReleaseScreen(screenNumber);
	}
	else {	
		PsychCloseWindow(windowRecord);
	}

	return(PsychError_none);
}
