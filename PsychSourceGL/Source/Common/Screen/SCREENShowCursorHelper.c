/*
	SCREENShowCursorHelper.c	
  
	AUTHORS:
	
		Allen.Ingling@nyu.edu		awi 
		mario.kleiner@tuebingen.mpg.de  mk

	PLATFORMS:
	
		All.
    
	HISTORY:
	
		12/14/02	awi		Created.
		10/12/04	awi		Changed "SCREEN" to "Screen" in useString.  Added ShowCursor to see also.  Fixed inappropriate "HideCursor" in this file.
		11/16/04	awi		Added "Helper" suffix. 
 
    
 
	DESCRIPTION:
  
		Hides the mouse pointer

  
	TO DO:
  
		The showcursor.m script should be modified to invoke hidecursor through screen only on OS X.  Otherwise it is a 
		separate mex file.  
  
		If we store a pointer to the window and screen lists in shared memory (accessible between mex files) then 
		we could implement this a a seperate mex file.  Seems more natural to do it from SCREEN though.  Let's wait and
		see what Windows and Linux are like before changeing it to a separate mex file.  
  
*/


#include "Screen.h"

#if PSYCH_SYSTEM == PSYCH_LINUX
#include <X11/X.h>
#include <X11/cursorfont.h>
#endif

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('ShowCursorHelper', screenIndex [, cursorshapeid][, mouseIndex]);";
//													  1			    2
static char synopsisString[] = 
	"This is a helper function called by ShowCursor.  Do not call Screen(\'ShowCursorHelper\'), use "
	"ShowCursor instead.\n"
	"Show the mouse pointer. If optional 'cursorshapeid' is given, select a specific cursor shape as well.\n"
	"If optional 'mouseIndex' is given, setup cursor for given master mouseIndex device (Linux only).\n";
static char seeAlsoString[] = "HideCursorHelper";
	 
PsychError SCREENShowCursorHelper(void) 
{
	int	screenNumber, cursorid, mouseIdx;
#if PSYCH_SYSTEM == PSYCH_LINUX
	Cursor  mycursor;
#endif
	//all subfunctions should have these two lines.  
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	PsychErrorExit(PsychCapNumInputArgs(3));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs
        
	PsychCopyInScreenNumberArg(1, TRUE, &screenNumber);

	mouseIdx = -1;
	PsychCopyInIntegerArg(3, FALSE, &mouseIdx);

	PsychShowCursor(screenNumber, mouseIdx);
	
	// Copy in optional cursor shape id argument: The default of -1 means to
	// not change cursor appearance. Any other positive value changes to an
	// OS dependent shape (== the mapping of numbers to shapes is OS dependent).
	cursorid = -1;
	PsychCopyInIntegerArg(2, FALSE, &cursorid);
	
	// Cursor change request?
	if (cursorid!=-1) {
		// Yes.
#if PSYCH_SYSTEM == PSYCH_OSX
		// OS/X:
		PsychCocoaSetThemeCursor(cursorid);
#endif

#if PSYCH_SYSTEM == PSYCH_LINUX
		// GNU/Linux with X11 windowing system:
		
		// Map screenNumber to X11 display handle and screenid:
		CGDirectDisplayID dpy;
		PsychGetCGDisplayIDFromScreenNumber(&dpy, screenNumber);
		// Create cursor spec from passed cursorid:
		mycursor = XCreateFontCursor(dpy, (unsigned int) cursorid);
		if (mouseIdx < 0) {
			// Set cursor for our window:
			PsychOSDefineX11Cursor(screenNumber, -1, mycursor);
		} else {
			// XInput cursor: Master pointers only.
			int nDevices;
			XIDeviceInfo* indevs = PsychGetInputDevicesForScreen(screenNumber, &nDevices);

			// Sanity check:
			if (NULL == indevs) PsychErrorExitMsg(PsychError_user, "Sorry, your system does not support individual mouse pointers.");
			if (mouseIdx >= nDevices) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");
			if (indevs[mouseIdx].use != XIMasterPointer) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such master cursor pointer.");

			PsychOSDefineX11Cursor(screenNumber, indevs[mouseIdx].deviceid, mycursor);
		}

		XFlush(dpy);

		// Done (hopefully).
#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS
		// Microsoft Windows:
		LPCTSTR lpCursorName;
		
		#ifndef IDC_HAND
		#define IDC_HAND MAKEINTRESOURCE(32649)
		#endif

		// Map provided cursor id to a Windows system id for such a cursor:
		switch(cursorid) {
			case 0:
				// Standard arrow cursor:
				lpCursorName = IDC_ARROW;
				break;

			case 1:
				// haircross cursor:
				lpCursorName = IDC_CROSS;
				break;

			case 2:
				// hand cursor:
				lpCursorName = IDC_HAND;
				break;

			case 3:
				// Arrows in all 4 directions cursor:
				lpCursorName = IDC_SIZEALL;
				break;

			case 4:
				// north-south cursor:
				lpCursorName = IDC_SIZENS;
				break;

			case 5:
				// east-west cursor:
				lpCursorName = IDC_SIZEWE;
				break;

			case 6:
				// hourglass cursor:
				lpCursorName = IDC_WAIT;
				break;

			case 7:
				// No cursor:
				lpCursorName = IDC_NO;
				break;

			default:
				// Default for unknown id is the standard arrow cursor:
				lpCursorName = IDC_ARROW;
		}
		
		// Load and set a cursor, based on the selected lpCursorName cursor id string:
		SetCursor(LoadCursor(NULL, lpCursorName));
#endif
		// End of cursor shape setup.
	}
	
	return(PsychError_none);	
}
