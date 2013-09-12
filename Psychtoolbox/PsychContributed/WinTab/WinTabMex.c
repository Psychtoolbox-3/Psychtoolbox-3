/*------------------------------------------------------------------------------
WinTabMex.c -- A simple Matlab MEX file for polling touch tablets that use the
WinTab API on Microsoft Windows.

WinTabMex.c is derived from the "Rule" example code from the wtkit126.zip
Wintab SDK from Wacom. You will need to download and install the freely available
WinTab SDK from Wacom in order to compile this MEX file.

Compile with:

mex -v WinTabMex.c -IPATHTOWTK\INCLUDE 
where PATHTOWTK must be the filesystem path to the installation folder of
the WintabSDK, aka where the file wintab.h is

Example: Assume you copied the WinTabMex.c file into the main folder of
the Wintab SDK, then you can simply cd into that directory in Matlab,
then type:

mex -v WinTabMex.c -IINCLUDE

If you want to compile on Octave, do this:

mex -v WinTabMex.c -DPTBOCTAVE3MEX -IINCLUDE


------------------------------------------------------------------------------

    WinTabMex.c is Copyright (C) 2008 - 2013 Mario Kleiner
    with contributions       (C) 2013 by Jason Friedman

    It is licensed to you under the MIT free and open source software license:

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to permit
    persons to whom the Software is furnished to do so, subject to the
    following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
    NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
    DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
    USE OR OTHER DEALINGS IN THE SOFTWARE.

------------------------------------------------------------------------------*/

/* Windows includes: */
#include <windows.h>
#include <stdlib.h>
#include <string.h>

/* WinTab includes: */
#include <wintab.h>

#ifdef USE_X_LIB
#include <wintabx.h>
#endif

/* Data format of returned packets: (x,y,z,buttons,serialno,timestamp,pressure): */
#define PACKETDATA	(PK_X | PK_Y | PK_Z | PK_BUTTONS | PK_SERIAL_NUMBER | PK_TIME | PK_STATUS | PK_CHANGED | PK_NORMAL_PRESSURE)
#define PACKETMODE	0
#include <pktdef.h>

/* Matlab includes: */
#include "mex.h"
#ifndef PTBOCTAVE3MEX
#include "matrix.h"
#endif

/* -------------------------------------------------------------------------- */
#define Inch2Cm	CASTFIX32(2.54)
#define Cm2Inch	CASTFIX32(1.0/2.54)
/* -------------------------------------------------------------------------- */

/* Local variables: */
static unsigned int initialized = 0;
static HCTX hTab = NULL;
static FIX32 scale[2];

// Define function pointer for the wintab functions we will use
typedef UINT ( API * WTINFOA ) ( UINT, UINT, LPVOID );
typedef UINT ( API * WTENABLE )  ( HCTX, BOOL );
typedef HCTX ( API * WTOPENA )( HWND, LPLOGCONTEXTA, BOOL );
typedef BOOL ( API * WTCLOSE ) ( HCTX );
typedef int  ( API * WTPACKETSGET ) (HCTX, int, LPVOID);
typedef BOOL ( API * WTQUEUESIZESET ) ( HCTX, int );
typedef int ( API * WTQUEUESIZEGET ) ( HCTX);

WTENABLE gpWTEnable = NULL;
WTOPENA gpWTOpenA = NULL;
WTINFOA gpWTInfoA = NULL;
WTCLOSE gpWTClose = NULL;
WTPACKETSGET gpWTPacketsGet = NULL;
WTQUEUESIZESET gpWTQueueSizeSet = NULL;
WTQUEUESIZEGET gpWTQueueSizeGet = NULL;

HINSTANCE ghWintab = NULL;

/* -------------------------------------------------------------------------- */
HCTX static NEAR TabletInit(HWND hWnd, FIX32 scale[])
{
	LOGCONTEXT lcMine;

	/* get default region */
	gpWTInfoA(WTI_DEFCONTEXT, 0, &lcMine);

	/* modify the digitizing region */
	strcpy(lcMine.lcName, "Psychtoolbox-3 WinTab MEX driver");
	lcMine.lcPktData = PACKETDATA;
	lcMine.lcPktMode = PACKETMODE;
	lcMine.lcMoveMask = 0;
	lcMine.lcBtnUpMask = lcMine.lcBtnDnMask;

	/* output in 1000ths of cm */
	lcMine.lcOutOrgX = lcMine.lcOutOrgY = 0;
	lcMine.lcOutExtX = INT(scale[0] * lcMine.lcInExtX);
	lcMine.lcOutExtY = INT(scale[1] * lcMine.lcInExtY);

	/* open the region */
	return gpWTOpenA(hWnd, &lcMine, FALSE);
}

/* -------------------------------------------------------------------------- */
/* return scaling factors in thousandths of cm per axis unit */
static void TabletScaling(FIX32 scale[])
{
	AXIS aXY[2];
	int i;
	UINT wDevice;

	/* get the data */
	gpWTInfoA(WTI_DEFCONTEXT, CTX_DEVICE, &wDevice);
	gpWTInfoA(WTI_DEVICES+wDevice, DVC_X, &aXY[0]);
	gpWTInfoA(WTI_DEVICES+wDevice, DVC_Y, &aXY[1]);

	/* calculate the scaling factors */
	for (i = 0; i < 2; i++) {
		FIX_DIV(scale[i], CASTFIX32(1000), aXY[i].axResolution);
		if (aXY[i].axUnits == TU_INCHES) {
			FIX_MUL(scale[i], scale[i], Inch2Cm);
		}
	}
}

/* -------------------------------------------------------------------------- */

/* Shutdown Wintab context, release all ressources: */
void shutDownWinTab(void)
{
	if (!initialized) return;
	
	/* Release WinTAB context: */
	gpWTClose(hTab);
	hTab = NULL;
	
	/* Release library: */
	#ifdef USE_X_LIB
	_UnlinkWintab();
	#endif

	if (ghWintab)
	  FreeLibrary( ghWintab );

	ghWintab = NULL;

	mexPrintf("WinTabMex shutdown complete.\n");
	
	/* Reset state to offline: */
	initialized = 0;
	
	return;
}

/* -------------------------------------------------------------------------- */

/* Load the wintab library and the functions used */
void loadWinTab(void) {
        ghWintab = LoadLibrary( "Wintab32.dll" );
	if ( !ghWintab ) {
	  mexErrMsgTxt("Unable to load Wintab32.dll library. Is it installed on your system?");
	}

	gpWTEnable = (WTENABLE) GetProcAddress( ghWintab, "WTEnable" );
	if ( !gpWTEnable ) { mexErrMsgTxt("Error load Wintab library."); }
	gpWTOpenA = (WTOPENA) GetProcAddress( ghWintab, "WTOpenA" );
	if ( !gpWTOpenA ) { mexErrMsgTxt("Error load Wintab library."); }
	gpWTInfoA = (WTINFOA) GetProcAddress( ghWintab, "WTInfoA" );
	if ( !gpWTInfoA ) { mexErrMsgTxt("Error load Wintab library."); }
	gpWTClose = (WTCLOSE) GetProcAddress( ghWintab, "WTClose" );
	if ( !gpWTClose ) { mexErrMsgTxt("Error load Wintab library."); }
	gpWTPacketsGet = (WTPACKETSGET) GetProcAddress( ghWintab, "WTPacketsGet" );
	if ( !gpWTPacketsGet ) { mexErrMsgTxt("Error load Wintab library."); }
	gpWTQueueSizeSet = (WTQUEUESIZESET) GetProcAddress( ghWintab, "WTQueueSizeSet" );
	if ( !gpWTQueueSizeSet ) { mexErrMsgTxt("Error load Wintab library."); }
	gpWTQueueSizeGet = (WTQUEUESIZEGET) GetProcAddress( ghWintab, "WTQueueSizeGet" );
	if ( !gpWTQueueSizeGet ) { mexErrMsgTxt("Error load Wintab library."); }
}

/* This is the main entry point from Matlab: */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	double* out;
	int cmd;
	
	/* Window handle to attach the WinTAB context to: */
	HWND hWnd;

	/* pkt buffer for single WinTab event packet: */
	PACKET pkt;

	/* Name of search string for parent window to attach to: */
	char winName[1000];
	int winid = 10;
	int show = 1;

	const char* me = mexFunctionName();

	if(nrhs<1) {
		mexPrintf("WinTabMex: A simple MEX file for driving touch tablets which are\n");
		mexPrintf("compatible with the WinTAB API on Microsoft Windows.\n\n");
		mexPrintf("(C) 2008-2013 by Mario Kleiner and Jason Friedman -- Licensed to you under MIT license.\n");
		mexPrintf("This file is part of Psychtoolbox-3 but should also work independently.\n");
		mexPrintf("\n");
		mexPrintf("Usage:\n\n");
		mexPrintf("%s(0, windowHandleOrNameString [, show=1]);\n", me);
		mexPrintf("- Initialize connection to tablet. You must provide either a Psychtoolbox onscreen window handle,\n");
		mexPrintf("or the name of the titlebar of another window, e.g., the Matlab window, in 'windowHandleOrNameString'.\n");
		mexPrintf("The optional 'show' flag if provided, can be set to 1 (default) to show, or zero to hide the window.\n");
		mexPrintf("\nAfter a successfull init, you can call the enable function below to start data collection.\n\n");
		mexPrintf("%s(1); - Shutdown the driver, close connection to tablet.\n\n", me);
		mexPrintf("%s(2); - Clear all pending tablet events and start data collection / recording into the event queue.\n\n", me);
		mexPrintf("%s(3); - Stop data collection / recording into the event queue.\n\n", me);
		mexPrintf("[old, new] = %s(4 [, reqSize]); - Report, and optionally try to resize the internal event queue to 'reqSize' entries.\n", me);
		mexPrintf("Report queue size before/after resize in 'old' and 'new'. You must check 'old' and 'new' for a valid size, as\n");
		mexPrintf("this method can fail if you asked for too much. If it failed, retry with smaller 'reqSize' until success.\n\n");
		mexPrintf("pkt = %s(5); - Try to fetch the oldest tablet event packet from the internal event queue.\n", me);
		mexPrintf("This polls for available data, returning at most one packet from the queue:\n");
		mexPrintf("If no new data is available, it will return an empty matrix, aka isempty(pkt) will be true.\n");
		mexPrintf("If data is available, it will return a 8-by-1 double vector with the following elements in the different rows:\n");
		mexPrintf("1 = x-position, 2 = y-position, 3 = z-position (in units of 1/1000th of centimeters).\n");
		mexPrintf("4 = Button state: If you do button = uint32(pkt(4)), you'll get an integer: Each bit in the integer describes\n");
		mexPrintf("state of one button. Mapping is device specific: E.g., bitget(button, 3) would tell you the state of the 3rd button.\n");
		mexPrintf("5 = Serial number. A monotonically increasing serial number. Should increase without gap, otherwise you lost events!\n");
		mexPrintf("6 = Timestamp. Time of when this tablet event was detected, in units of milliseconds. Baseline and accuracy are unknown.\n");
		mexPrintf("7 = Status: An encoded status code: status = uint32(pkt(7)). Then bitget(status, 1) will be non-zero if the tool has left\n");
		mexPrintf("    the active area of the tablet. bitget(status, 2) will be non-zero if packets were lost due to queue-overflow.\n");
		mexPrintf("    bitget(status, 3) will be non-zero if the pen is at the margin of the active area.\n");
		mexPrintf("8 = Changed: Info about what changed in this packet wrt. to previous packet: changed=uint32(pkt(8));\n");
		mexPrintf("    bitget(changed, 2) means that the 'Status' field has changed. bitget(changed, 7) == button changed,\n");
		mexPrintf("    bitget(changed, 8:10) means x, y or z position has changed.\n");
		mexPrintf("9 = Normal pressure: Info about the pressure along the tablet surface normal, or pressure on the pen tip, in arbitrary units.\n");
		mexPrintf("\n");		
		mexPrintf("\n");
		mexPrintf("This driver will only work if the WinTab libraries and driver are properly installed on your system. Good luck!\n\n\n");
        return;
	}
	
	// Load the library
	if (!ghWintab) loadWinTab();

	/* First argument must be the command code: */
	cmd = (int) mxGetScalar(prhs[0]);

	if (cmd !=0 && !initialized) mexErrMsgTxt("Tried to call WinTab function, but it isn't initialized yet! Call code 0 (Initialize) first!");
	
	switch(cmd) {
		case 0: // Init.
			if (initialized) mexErrMsgTxt("Tried to initialize WinTab, but it is already initialized! Call code 1 (Shutdown) first!");
			
			if (nrhs<2) mexErrMsgTxt("No window handle or window name provided as 2nd argument! This is required!");
			
			/* Register exit handler: */
			mexAtExit(shutDownWinTab);
			
			/* Try to get window handle of parent window: */
			if (mxIsDouble(prhs[1])) {
				/* PTB window handle: Translate into window name: */
				winid = (int) mxGetScalar(prhs[1]);
				sprintf(winName, "PTB Onscreen window [%i]:", winid);

				/* Find window with matching title, return window handle: */
				hWnd = FindWindow("PTB-OpenGL", winName);
			}
			else {
				/* Name string: Just copy it. */
				winName[0] = 0;
				mxGetString(prhs[1], winName, sizeof(winName)-1);

				/* Find window with matching title, return window handle: */
				hWnd = FindWindow(NULL, winName);
			}
			
			if (hWnd!=NULL) {
				mexPrintf("Will try to attach WinTab to Window with title %s .\n", winName);
				if(nrhs > 2) show = (int) mxGetScalar(prhs[2]);
				ShowWindow(hWnd, (show) ? SW_SHOW : SW_HIDE);
			}
			else {
				mexPrintf("While trying to attach WinTab to Window with title %s .\n", winName);
				mexErrMsgTxt("Failed to find a window with this title! WinTab init failed.\n");
			}
			
			/* Ok, we have a window handle. Initialize the Wintab context and connection to the tablet with this handle: */
			TabletScaling(scale);
			if ((hTab = TabletInit(hWnd, scale)) == NULL) mexErrMsgTxt("Failed to initialize WinTab context for window! WinTab init failed.\n");

			initialized = 1;
			
			// We're ready to go,.
			mexPrintf("WinTabMex tablet driver ready.\n");
		break;
		
		case 1: // Shutdown.
			shutDownWinTab();
		break;
		
		case 2: // Clear event queue and enable scanning.
			// Empty queue of pending events, if any:
			gpWTPacketsGet(hTab, gpWTQueueSizeGet(hTab), NULL);
			
			// Enable scanning:
			gpWTEnable(hTab, (BOOL) 1);
		break;
		
		case 3: // Disable scanning.
		        gpWTEnable(hTab, (BOOL) 0);
		break;
		
		case 4: // Resize queue:
			plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
			*mxGetPr(plhs[0]) = (double) gpWTQueueSizeGet(hTab);
			if (nrhs >= 2) gpWTQueueSizeSet(hTab, (int) mxGetScalar(prhs[1]));
			plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
			*mxGetPr(plhs[1]) = (double) gpWTQueueSizeGet(hTab);
		break; 
		
		case 5: // Query oldest event packet:
			// Any news from the tablet? Poll it:
			if (gpWTPacketsGet(hTab, 1, &pkt)) {
				// Yep. Return packet data to Matlab:
				plhs[0]=mxCreateDoubleMatrix(9, 1, mxREAL);
				out = mxGetPr(plhs[0]);
				
				out[0] = (double) pkt.pkX;
				out[1] = (double) pkt.pkY;
				out[2] = (double) pkt.pkZ;
				out[3] = (double) pkt.pkButtons;
				out[4] = (double) pkt.pkSerialNumber;
				out[5] = (double) pkt.pkTime;
				out[6] = (double) pkt.pkStatus;
				out[7] = (double) pkt.pkChanged;
				out[8] = (double) pkt.pkNormalPressure;
			}
			else {
				// Nope. Return an empty matrix to signal no news from WinTab:
				plhs[0]=mxCreateDoubleMatrix(0, 0, mxREAL);
			}
		break;
		
		default:
			mexErrMsgTxt("Invalid / Unknown command code provided as first argument!");
	}
	
	// Done.
	return;
}
