/*
	PsychToolbox3/Source/windows/Screen/PsychWindowGlue.c
	
	PLATFORMS:	
	
		This is the Windows  version only.  
				
	AUTHORS:
	
                Allen Ingling           awi             Allen.Ingling@nyu.edu
                Mario Kleiner           mk              mario.kleiner at tuebingen.mpg.de

	HISTORY:
	
		12/20/02		awi		Wrote it mostly by modifying SDL-specific refugees (from an experimental SDL-based Psychtoolbox).
		11/16/04		awi		Added description.
                 4/22/05                mk              Added support for OpenGL stereo windows and enhanced Flip-behaviour:
                                                        Flip'ing at specified deadline, retaining instead of clear'ing backbuffer during flip,
                                                        return of stimulus onset related timestamps, properly syncing to VBL.
                 4/29/05                mk              Additional safety checks for VBL sync in PsychOpenOnscreenWindow().
                 5/14/05                mk              Additional safety checks for insufficient gfx-hardware and multi-display setups,
                                                        failing beam-position queries. New refresh interval estimation code, reworked Flip.
                 5/19/05                mk              Extra check for 'flipwhen' values over 1000 secs in future: Abort.
                 5/30/05                mk              Check for Screen('Preference', 'SkipSyncTests', 1) -> Shortened tests, if set.
                 6/09/05                mk              Experimental support for busy-waiting for VBL and for multi-flip of stereo displays.
                 9/30/05                mk              Added PsychRealtimePriority for improving timing tests in PsychOpenWindow()
                 9/30/05                mk              Added check for Screen('Preference', 'VisualDebugLevel', level) -> Amount of vis. feedback.
                 10/10/05               mk              Important Bugfix for PsychRealtimePriority() - didn't switch back to non-RT priority!!
		 10/19/05		awi		Cast NULL to CGLPixelFormatAttribute type to make the compiler happy.
 
	DESCRIPTION:
	
		Functions in this file comprise an abstraction layer for probing and controlling window state, except for window content.  
		
		Each C function which implements a particular Screen subcommand should be platform neutral.  For example, the source to SCREENPixelSizes() 
		should be platform-neutral, despite that the calls in OS X and Windows to detect available pixel sizes are different.  The platform 
		specificity is abstracted out in C files which end it "Glue", for example PsychScreenGlue, PsychWindowGlue, PsychWindowTextClue.

	NOTES:
	
		Technical information about the Windows-Vista display driver model (WDDM 1.0) and later: http://msdn.microsoft.com/en-us/library/ms796755.aspx
		Contains details about driver interfaces, memory management, GPU scheduling, swap control etc.
		Official name: "Windows Vista Display Driver Model Design Guide"

		Useful information about the Vista et al. Direct-X Graphics kernel subsystem: http://msdn.microsoft.com/en-us/library/ee415671(VS.85).aspx
		OpenGL ICD's on Vista and later communicate with the GPU (or more accurately, the Display Miniport Driver of the GPU) through the
		Direct-X Graphics kernel subsystem. For that reason it is helpful to understand its internal working, as some of its properties and
		of the DirectX/Direct3D driver settings may also affect OpenGL rendering and stimulus presentation.

	TO DO: 

*/

#include "Screen.h"

#ifndef WS_EX_LAYERED
/* Define prototype for this function: */
BOOL SetLayeredWindowAttributes(      
    HWND hwnd,
    COLORREF crKey,
    BYTE bAlpha,
    DWORD dwFlags);
#endif

// Add missing defines from Vista SDK for the DWM, if not already defined:

// Need to enable 1-Byte alignment of structures:
#pragma pack(push,1)

typedef struct UNSIGNED_RATIO {
    UINT32 uiNumerator;
    UINT32 uiDenominator;
} UNSIGNED_RATIO;

typedef ULONGLONG DWM_FRAME_COUNT;
typedef ULONGLONG QPC_TIME;

typedef  struct DWM_TIMING_INFO {
    UINT32          cbSize;

    // Data on DWM composition overall
    
    // Monitor refresh rate
    UNSIGNED_RATIO  rateRefresh;

    // Actual period
    QPC_TIME        qpcRefreshPeriod;

    // composition rate     
    UNSIGNED_RATIO  rateCompose;

    // QPC time at a VSync interupt
    QPC_TIME        qpcVBlank;

    // DWM refresh count of the last vsync
    // DWM refresh count is a 64bit number where zero is
    // the first refresh the DWM woke up to process
    DWM_FRAME_COUNT cRefresh;

    // DX refresh count at the last Vsync Interupt
    // DX refresh count is a 32bit number with zero 
    // being the first refresh after the card was initialized
    // DX increments a counter when ever a VSync ISR is processed
    // It is possible for DX to miss VSyncs
    //
    // There is not a fixed mapping between DX and DWM refresh counts
    // because the DX will rollover and may miss VSync interupts
    UINT cDXRefresh;

    // QPC time at a compose time.  
    QPC_TIME        qpcCompose;

    // Frame number that was composed at qpcCompose
    DWM_FRAME_COUNT cFrame;

    // The present number DX uses to identify renderer frames
    UINT            cDXPresent;

    // Refresh count of the frame that was composed at qpcCompose
    DWM_FRAME_COUNT cRefreshFrame;


    // DWM frame number that was last submitted
    DWM_FRAME_COUNT cFrameSubmitted;

    // DX Present number that was last submitted
    UINT cDXPresentSubmitted;

    // DWM frame number that was last confirmed presented
    DWM_FRAME_COUNT cFrameConfirmed;

    // DX Present number that was last confirmed presented
    UINT cDXPresentConfirmed;

    // The target refresh count of the last
    // frame confirmed completed by the GPU
    DWM_FRAME_COUNT cRefreshConfirmed;

    // DX refresh count when the frame was confirmed presented
    UINT cDXRefreshConfirmed;

    // Number of frames the DWM presented late
    // AKA Glitches
    DWM_FRAME_COUNT          cFramesLate;
    
    // the number of composition frames that 
    // have been issued but not confirmed completed
    UINT          cFramesOutstanding;


    // Following fields are only relavent when an HWND is specified
    // Display frame


    // Last frame displayed
    DWM_FRAME_COUNT cFrameDisplayed;

    // QPC time of the composition pass when the frame was displayed
    QPC_TIME        qpcFrameDisplayed; 

    // Count of the VSync when the frame should have become visible
    DWM_FRAME_COUNT cRefreshFrameDisplayed;

    // Complete frames: DX has notified the DWM that the frame is done rendering

    // ID of the the last frame marked complete (starts at 0)
    DWM_FRAME_COUNT cFrameComplete;

    // QPC time when the last frame was marked complete
    QPC_TIME        qpcFrameComplete;

    // Pending frames:
    // The application has been submitted to DX but not completed by the GPU
 
    // ID of the the last frame marked pending (starts at 0)
    DWM_FRAME_COUNT cFramePending;

    // QPC time when the last frame was marked pending
    QPC_TIME        qpcFramePending;

    // number of unique frames displayed
    DWM_FRAME_COUNT cFramesDisplayed;

    // number of new completed frames that have been received
    DWM_FRAME_COUNT cFramesComplete;

     // number of new frames submitted to DX but not yet complete
    DWM_FRAME_COUNT cFramesPending;

    // number of frames available but not displayed, used or dropped
    DWM_FRAME_COUNT cFramesAvailable;

    // number of rendered frames that were never
    // displayed because composition occured too late
    DWM_FRAME_COUNT cFramesDropped;
    
    // number of times an old frame was composed 
    // when a new frame should have been used
    // but was not available
    DWM_FRAME_COUNT cFramesMissed;
    
    // the refresh at which the next frame is
    // scheduled to be displayed
    DWM_FRAME_COUNT cRefreshNextDisplayed;

    // the refresh at which the next DX present is 
    // scheduled to be displayed
    DWM_FRAME_COUNT cRefreshNextPresented;

    // The total number of refreshes worth of content
    // for this HWND that have been displayed by the DWM
    // since DwmSetPresentParameters was called
    DWM_FRAME_COUNT cRefreshesDisplayed;
	
    // The total number of refreshes worth of content
    // that have been presented by the application
    // since DwmSetPresentParameters was called
    DWM_FRAME_COUNT cRefreshesPresented;


    // The actual refresh # when content for this
    // window started to be displayed
    // it may be different than that requested
    // DwmSetPresentParameters
    DWM_FRAME_COUNT cRefreshStarted;

    // Total number of pixels DX redirected
    // to the DWM.
    // If Queueing is used the full buffer
    // is transfered on each present.
    // If not queuing it is possible only 
    // a dirty region is updated
    ULONGLONG  cPixelsReceived;

    // Total number of pixels drawn.
    // Does not take into account if
    // if the window is only partial drawn
    // do to clipping or dirty rect management 
    ULONGLONG  cPixelsDrawn;

    // The number of buffers in the flipchain
    // that are empty.   An application can 
    // present that number of times and guarantee 
    // it won't be blocked waiting for a buffer to 
    // become empty to present to
    DWM_FRAME_COUNT      cBuffersEmpty;

} DWM_TIMING_INFO;

typedef enum _DWM_SOURCE_FRAME_SAMPLING {
    DWM_SOURCE_FRAME_SAMPLING_POINT = 1,
    DWM_SOURCE_FRAME_SAMPLING_COVERAGE,
    DWM_SOURCE_FRAME_SAMPLING_LAST
} DWM_SOURCE_FRAME_SAMPLING;

typedef struct _DWM_PRESENT_PARAMETERS {
    UINT32 cbSize;
    BOOL fQueue;
    DWM_FRAME_COUNT cRefreshStart;
    UINT cBuffer;
    BOOL fUseSourceRate;
    UNSIGNED_RATIO rateSource;
    UINT cRefreshesPerFrame;
    DWM_SOURCE_FRAME_SAMPLING eSampling;
} DWM_PRESENT_PARAMETERS;

#pragma pack(pop)

// Application instance handle:
static HINSTANCE hInstance = 0;

// Number of currently open onscreen windows:
static int win32_windowcount = 0;

// Mouse button states:
static psych_bool mousebutton_l=FALSE;
static psych_bool mousebutton_m=FALSE;
static psych_bool mousebutton_r=FALSE;

// Module handle for the DWM library 'dwmapi.dll': Or 0 if unsupported.
HMODULE dwmlibrary = 0;

// dwmSupported is true if DWM is supported and library is linked:
static psych_bool dwmSupported = FALSE;

// DWM function definitions and procpointers:
typedef HRESULT (APIENTRY *DwmIsCompositionEnabledPROC)(BOOL *pfEnabled);
typedef HRESULT (APIENTRY *DwmEnableCompositionPROC)(UINT enable);
typedef HRESULT (APIENTRY *DwmEnableMMCSSPROC)(BOOL fEnableMMCSS);
typedef HRESULT (APIENTRY *DwmGetCompositionTimingInfoPROC)(HWND hwnd, DWM_TIMING_INFO* pTimingInfo);
typedef HRESULT (APIENTRY *DwmSetPresentParametersPROC)(HWND hwnd, DWM_PRESENT_PARAMETERS *pPresentParams);
typedef HRESULT (APIENTRY *DwmFlushPROC)(void);

DwmIsCompositionEnabledPROC			PsychDwmIsCompositionEnabled = NULL;
DwmEnableCompositionPROC			PsychDwmEnableComposition = NULL;
DwmEnableMMCSSPROC					PsychDwmEnableMMCSS = NULL;
DwmGetCompositionTimingInfoPROC		PsychDwmGetCompositionTimingInfo = NULL;
DwmSetPresentParametersPROC			PsychDwmSetPresentParameters = NULL;
DwmFlushPROC						PsychDwmFlush = NULL;

char hostwinName[512];
HWND hostwinHandle = 0;	

/** PsychHostWindowEnumFunc(): Called by EnumWindows() in PsychOpenOnscreenWindow()
 *  as callback function for window enumeration. This is meant to filter out the
 *  main window of Matlab or Octave from the enumerated candidate windows.
 *
 *  If it finds the target window, it assigns its handle to 'hostwinHandle' for
 *  later use and stops enumeration.
 */
BOOL CALLBACK PsychHostWindowEnumFunc(HWND hwnd, LPARAM passId)
{
	// Get text in title bar of Window hwnd as string for pattern matching:
    GetWindowText(hwnd, hostwinName, sizeof(hostwinName) - 1);
	if (PsychPrefStateGet_Verbosity() > 15) printf("PTB-DEBUG: Window enumeration: Pass %i, HWND = %p Current: %s\n", passId, hwnd, hostwinName);

	// Which runtime?
	#ifndef PTBOCTAVE3MEX
		// Running on Matlab: Use Matlab name matching:
		// Pass 1: Scan for Matlab in GUI mode:
		if ((passId == 1) && (strstr(hostwinName, "MATLAB  ") || strstr(hostwinName, "MATLAB R20"))) {
			// Found something that looks ok:
			hostwinHandle = hwnd;
			return(FALSE);
		}
		
		// Pass 2: Scan for Matlab in Console mode:
		if ((passId == 2) && strstr(hostwinName, "MATLAB Command Window")) {
			// Found something that looks ok:
			hostwinHandle = hwnd;
			return(FALSE);
		}
	#else
		// Running on Octave: Use Octave name matching:
		// Pass 1: Scan for Octave in Console mode or for QtOctave in GUI mode:
		if ((passId == 1) && strstr(hostwinName, "Octave")) {
			// Found something that looks ok:
			hostwinHandle = hwnd;
			return(FALSE);
		}
	#endif
	
	// Nothing yet? Continue enumeration:
    return(TRUE);
}


/** PsychRealtimePriority: Temporarily boost priority to highest available priority in M$-Windows.
    PsychRealtimePriority(true) enables realtime-scheduling (like Priority(2) would do in Matlab).
    PsychRealtimePriority(false) restores scheduling to the state before last invocation of PsychRealtimePriority(true),
    it undos whatever the previous switch did.

    We switch to RT scheduling during PsychGetMonitorRefreshInterval() and a few other timing tests in
    PsychOpenWindow() to reduce measurement jitter caused by possible interference of other tasks.
*/
psych_bool PsychRealtimePriority(psych_bool enable_realtime)
{
    HANDLE	   currentProcess;
    static psych_bool old_enable_realtime = FALSE;
    static DWORD   oldPriority = NORMAL_PRIORITY_CLASS;
    const  DWORD   realtime_class = REALTIME_PRIORITY_CLASS;

    if (old_enable_realtime == enable_realtime) {
        // No transition with respect to previous state -> Nothing to do.
        return(true);
    }
    
    // Transition requested:
    old_enable_realtime = enable_realtime;
    
    // Determine our processID:
    currentProcess = GetCurrentProcess();
    
    if (enable_realtime) {
      // Transition to realtime requested:
      
      // Get current scheduling policy and back it up for later restore:
      oldPriority = GetPriorityClass(currentProcess);
        
      // Check if realtime scheduling isn't already active.
      // If we are already in RT mode (e.g., Priority(2) call in Matlab), we skip the switch...
      if (oldPriority != realtime_class) {
			// RT scheduling not yet active -> Switch to it.
			SetPriorityClass(currentProcess, realtime_class);
			
			// Check if transition to realtime class worked. If not then retry
			// with high priority class -- better than nothing...
			if (GetPriorityClass(currentProcess)!=realtime_class) {
				if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: Could not enable REALTIME_PRIORITY_CLASS scheduling! Retrying with HIGH_PRIORITY_CLASS...\n");
				SetPriorityClass(currentProcess, HIGH_PRIORITY_CLASS);
				if ((GetPriorityClass(currentProcess)!=HIGH_PRIORITY_CLASS) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: Could not enable HIGH_PRIORITY_CLASS scheduling. OS malfunction? Timing may be noisy...\n");
			}
			
			// Additionally try to schedule us MMCSS: This will lift us roughly into the
			// same scheduling range as REALTIME_PRIORITY_CLASS, even if we are non-admin users
			// on Vista and Windows-7 and later, however with a scheduler safety net applied.
			PsychSetThreadPriority((psych_thread*) 0x1, 10, 0);			
      }
    }
    else {
      // Transition from RT to whatever-it-was-before scheduling requested: We just reestablish the backed-up old
      // policy: If the old policy wasn't Non-RT, then we don't switch back...
      SetPriorityClass(currentProcess, oldPriority);
	  
	  // Disable any MMCSS scheduling for us if new priority class is non-RT, bog-standard normal scheduling:
	  if (oldPriority == NORMAL_PRIORITY_CLASS) PsychSetThreadPriority((psych_thread*) 0x1, 0, 0);
    }

    // Success.
    return(TRUE);
}

// Perform OS specific processing of Window events:
void PsychOSProcessEvents(PsychWindowRecordType *windowRecord, int flags)
{
	POINT	lPoint;
	RECT	lRect;

	// Trigger event queue dispatch processing for GUI windows:
	if (windowRecord == NULL || windowRecord->specialflags & kPsychGUIWindow) PsychGetMouseButtonState(NULL);

	if (windowRecord == NULL) {
		// Done, so far...
		return;
	}
		
	// GUI windows need to behave GUIyee:
	if ((windowRecord->specialflags & kPsychGUIWindow) && PsychIsOnscreenWindow(windowRecord)) {
		// Update windows rect and globalrect, based on current size and location:
		lPoint.x = lPoint.y = 0;
		ClientToScreen(windowRecord->targetSpecific.windowHandle, &lPoint);
		GetClientRect(windowRecord->targetSpecific.windowHandle, &lRect);
		PsychMakeRect(windowRecord->globalrect, lPoint.x, lPoint.y, lPoint.x + lRect.right - 1, lPoint.y + lRect.bottom - 1);
		PsychNormalizeRect(windowRecord->globalrect, windowRecord->rect);
		PsychSetupClientRect(windowRecord);
		PsychSetupView(windowRecord, FALSE);
	}		
}

// Callback handler for Window manager: Handles some events
LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	static PAINTSTRUCT ps;
	PsychWindowRecordType	**windowRecordArray;
	int i, numWindows;
	LRESULT res;
	int verbosity = PsychPrefStateGet_Verbosity();
	
	if (verbosity > 6) printf("PTB-DEBUG: WndProc(): Called!\n");
	
	// What event happened?
	switch(uMsg) {
		case WM_MOUSEACTIVATE:
			// Mouseclick into our inactive window (non-fullscreen) received. Eat it:
			if (verbosity > 6) printf("PTB-DEBUG: WndProc(): MOUSE ACTIVATION!\n");
            
			// Default to ignore activation event:
			res = MA_NOACTIVATEANDEAT;
            
			// Scan the list of windows to find onscreen window with handle hWnd:
			PsychCreateVolatileWindowRecordPointerList(&numWindows, &windowRecordArray);
			for(i = 0; i < numWindows; i++) {
				if (PsychIsOnscreenWindow(windowRecordArray[i]) &&
					windowRecordArray[i]->targetSpecific.windowHandle == hWnd) {
                    // This is our window, and if it is a GUI window then perform activation:
                    if (windowRecordArray[i]->specialflags & kPsychGUIWindow) res = MA_ACTIVATEANDEAT;
                    break;
				}
			}
			PsychDestroyVolatileWindowRecordPointerList(windowRecordArray);

			return(res);
			break;
			
		case WM_SYSCOMMAND:
			// System command received: We intercept system commands that would start
			// the screensaver or put the display into powersaving sleep-mode:
			if (verbosity > 6) printf("PTB-DEBUG: WndProc(): WM_SYSCOMMAND!\n");
			switch(wParam) {
				case SC_SCREENSAVE:
				case SC_MONITORPOWER:
					return(0);
			}
			break;
			
		case WM_LBUTTONDOWN:
			// Left mouse button depressed:
			mousebutton_l = TRUE;
			break;
			
		case WM_LBUTTONUP:
			// Left mouse button released:
			mousebutton_l = FALSE;
			break;
			
		case WM_MBUTTONDOWN:
			// Middle mouse button depressed:
			mousebutton_m = TRUE;
			break;
			
		case WM_MBUTTONUP:
			// Middle mouse button released:
			mousebutton_m = FALSE;
			break;
			
		case WM_RBUTTONDOWN:
			// Right mouse button depressed:
			mousebutton_r = TRUE;
			break;
			
		case WM_RBUTTONUP:
			// Right mouse button released:
			mousebutton_r = FALSE;
			break;
			
		case WM_PAINT:
			// Repaint event: This happens if a previously covered non-fullscreen window
			// got uncovered, so part of it needs to be redrawn. PTB's rendering model
			// doesn't have a concept of redrawing a stimulus. We do nothing here.
			if (verbosity > 6) printf("PTB-DEBUG: WndProc(): WM_PAINT!\n");
			break;
			
		case WM_SIZE:
			// Window resize event: Only happens in debug-mode (non-fullscreen).
			// We resize the viewport accordingly and then trigger a repaint-op.
			if (verbosity > 6) printf("PTB-DEBUG: WndProc(): WM_SIZE!\n");			
			break;
			
		case WM_CLOSE:
			// WM_CLOSE falls through to WM_CHAR and emulates an Abort-key press.
			// -> Manually closing an onscreen window does the same as pressing the Abort-key.
			if (verbosity > 6) printf("PTB-DEBUG: WndProc(): WM_PAINT!\n");
			wParam='@';
		case WM_CHAR:
			// Character received. We only care about one key, the '@' key.
			// Pressing '@' will immediately close all onscreen windows, show
			// the cursor and such. It is the emergency stop key.
			if (wParam=='@') {
				// Emergency shutdown:
				printf("\nPTB-INFO: Master-Abort key '@' pressed by user.\n");
				printf("PTB-INFO: Enforcing script abortion and restoring desktop by executing Screen('CloseAll') now!\n");
				printf("PTB-INFO: Please ignore the false error message (INTERNAL PSYCHTOOLBOX ERROR) caused by this...\n");
				ScreenCloseAllWindows();
				return(0);
			}
			else {
				// The next failed good idea :-( Transmitting keystrokes (characters) that we received in the
				// queue of our onscreen windows (due to us having keyboard focus) to the runtime (Matlab or Octave)
				// by simply posting/injecting the character events into the host apps window event queue doesn't work.
				// Well, the actual posting/injection of events does work, but the runtimes don't like this "alien"
				// events out of nowhere and don't handle them as expected: Octave does receive them almost properly,
				// but loses repeated characters, e.g., sequence mmm turns into a single m. Matlab in nojvm mode doesn't
				// work at all this way. Matlab in jmv/gui mode does receive them only after huge delays, so it totally
				// breaks tons of existing GetChar() code, including our GetCharTest correctness test.
				//
				// For this reason, we disable this new method by setting the following if(FALSE):
				if (FALSE) {
						// Other character received. Can we redispatch it to our host applications window?
						if (hostwinHandle) {
							// We have its window handle! Post a message to its message queue. This way,
							// CharAvail() and GetChar() may be able to indirectly receive characters that
							// were typed into our onscreen window's due to us stealing the keyboard input focus.
							// Especially important on MS-Vista and later due to our need to steal keyboard focus:
							if (PsychPrefStateGet_Verbosity() > 6) printf("PTB-DEBUG: WndProc(): WM_CHAR '%c' (lParam = %i) redispatched to hostapp-window.\n", (char) wParam, (int) lParam);
							PostMessage(hostwinHandle, WM_CHAR, wParam, lParam);
							//SendMessage(hostwinHandle, WM_CHAR, wParam, lParam);
						}
				}
				if (PsychPrefStateGet_Verbosity() > 16) printf("PTB-DEBUG: WndProc(): WM_CHAR '%c' (lParam = %i) received.\n", (char) wParam, (int) lParam);
			}
			break;
	}
	
    return DefWindowProc(hWnd, uMsg, wParam, lParam);
}

/* PsychGetMouseButtonState: Returns current mouse button up-/down state. Called by SCREENGetMouseHelper. */
void PsychGetMouseButtonState(double* buttonArray)
{
    MSG msg;

	// Run our message dispatch loop until we've processed all events in the event queue for our onscreen windows:
	// We need to do this periodically so WindowsXP and later are convinced that our window/application is still
	// alive and responsive. Otherwise it might conclude trouble and replace our onscreen window by a ghost window
	// which is effectively dead and useless for OpenGL rendering and display. See "About Messages and Message Queues"
	// at MSDN: http://msdn.microsoft.com/en-us/library/ms644927(VS.85).aspx for explanation of this braindead
	// mechanism...
	while(PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {		
		// Translate event in case it is a keyboard event that needs translation:
		TranslateMessage(&msg);
		
		// Dispatch to the WndProc() function above for final event handling:
		DispatchMessage(&msg);
	}

	// The routine expects a preallocated 3-element double matrix for the three mouse-buttons. If
	// NULL is passed, it doesn't execute the actual mouse query, but only performs the event processing
	// above:
	if (NULL != buttonArray) {
		// GetAsyncKeyState() directly returns the current state of the physical mouse
		// buttons, independent of window system event processing, keyboard or mouse
		// focus and such. Much more robust than old approach:
		buttonArray[0]=(double) ((GetAsyncKeyState(VK_LBUTTON) & -32768) ? 1 : 0);
		buttonArray[1]=(double) ((GetAsyncKeyState(VK_MBUTTON) & -32768) ? 1 : 0);
		buttonArray[2]=(double) ((GetAsyncKeyState(VK_RBUTTON) & -32768) ? 1 : 0);
	}
	
	return;
}

psych_bool ChangeScreenResolution (int screenNumber, int width, int height, int bitsPerPixel, int fps)	// Change The Screen Resolution
{
	DEVMODE dmScreenSettings; // Device mode structure
	int rc;
	
	// Clear structure:
	memset (&dmScreenSettings, 0, sizeof (DEVMODE));
	dmScreenSettings.dmSize		= sizeof (DEVMODE);
	dmScreenSettings.dmDriverExtra	= 0;
	
	// Query current display settings and init struct with them:
	EnumDisplaySettings(PsychGetDisplayDeviceName(screenNumber), ENUM_CURRENT_SETTINGS, &dmScreenSettings);
	
	// Override current settings with the requested settings, if any:
	if (width>0)  dmScreenSettings.dmPelsWidth  = width;  // Select Screen Width
	if (height>0) dmScreenSettings.dmPelsHeight = height; // Select Screen Height
	if (bitsPerPixel>0) dmScreenSettings.dmBitsPerPel = bitsPerPixel; // Select Bits Per Pixel
	if (fps>0) dmScreenSettings.dmDisplayFrequency = fps; // Select display refresh rate in Hz
	
	// All provided values should be honored: We need to spec the refresh explicitely,
	// because otherwise the system will select the lowest fps for a given display mode.
	dmScreenSettings.dmFields = DM_BITSPERPEL | DM_PELSWIDTH | DM_PELSHEIGHT | DM_DISPLAYFREQUENCY;
	
	if (PsychPrefStateGet_Verbosity() > 4) {
		printf("PTB-DEBUG: Switching display for screen %i to fullscreen mode with resolution %i x %i with bpc = %i @ %i Hz.\n", screenNumber, (int) dmScreenSettings.dmPelsWidth, (int) dmScreenSettings.dmPelsHeight, (int) dmScreenSettings.dmBitsPerPel, (int) dmScreenSettings.dmDisplayFrequency);
	}
	
	#ifndef DISP_CHANGE_BADDUALVIEW
	#define DISP_CHANGE_BADDUALVIEW -6 
	#endif
	
	// Perform the change:
	if ((rc = (int) ChangeDisplaySettingsEx(PsychGetDisplayDeviceName(screenNumber), &dmScreenSettings, NULL, CDS_FULLSCREEN, NULL)) != DISP_CHANGE_SUCCESSFUL) {
		if (PsychPrefStateGet_Verbosity() > 4) {
			printf("PTB-DEBUG: Fullscreen switch of display for screen %i FAILED! Reason given by OS: ", screenNumber);
			switch(rc) {
				case DISP_CHANGE_BADDUALVIEW:
					printf("DISP_CHANGE_BADDUALVIEW	 The settings change was unsuccessful because system is DualView capable.\n");
				break;
				
				case DISP_CHANGE_BADFLAGS:
					printf("DISP_CHANGE_BADFLAGS	An invalid set of flags was passed in.\n");
				break;
				
				case DISP_CHANGE_BADMODE:
					printf("DISP_CHANGE_BADMODE	The graphics mode is not supported.\n");
				break;
				
				case DISP_CHANGE_BADPARAM:
					printf("DISP_CHANGE_BADPARAM	An invalid parameter was passed in. This can include an invalid flag or combination of flags.\n");
				break; 
				
				case DISP_CHANGE_FAILED:
					printf("DISP_CHANGE_FAILED	The display driver failed the specified graphics mode.\n");
				break;
				
				case DISP_CHANGE_RESTART:
					printf("DISP_CHANGE_RESTART	The computer must be restarted in order for the graphics mode to work.\n");
				break;

				case DISP_CHANGE_NOTUPDATED:
					printf("DISP_CHANGE_NOTUPDATED	 Unable to write settings to the registry.\n");
				break;

				default:
					printf("Unknown error condition.\n");
			}
		}
		
		return(FALSE);	// Display Change Failed, Return False
	}

	// Display Change Was Successful, Return True
	return(TRUE);
}

/* Internal method, called by PsychOSOpenOnscreenWindow(): Optimize present parameters
 * for Windows DWM compositor, to minimize its interference with our presentation timing
 * in the cases where we can't avoid DWM compositing completely.
 *
 * All we do is disable queued presentation mode, as this is deprecated and didn't ever
 * work well at all in the first place. Also we set the number of buffers in the buffer
 * queue to the absolute allowable minimum of 2 buffers, trying to force the DWM to not
 * queue up anything.
 *
 * This optimization is inspired by the one used by Google's ANGLE project for
 * OpenGL rendering over DirectX-10. ANGLE is used for the WebGL implementation of
 * the Chromium and Firefox web browsers on Windows, so we assume they know what they're
 * doing. See the setup code for the Surface::initialize() method in ANGLE's libEGL
 * implementation: <code.google.com/p/angleproject/source/browse/trunk/src/libEGL/Surface.cpp>
 *
 * In any case, this is just band-aid, as the only decent way to get good performance
 * and timing is to not have the DWM active on the stimulus presentation display. We
 * achieve this by actively disabling the DWM on Windows Vista and Windows-7. On Windows-8
 * we can't do that, but apparently the DWM auto-disables at least for fullscreen onscreen
 * stimulus displays, similar to what Linux compositors do when unredirect_fullscreen_windows
 * mode is requested. This optimization only applies to half-transparent "debug windows" or
 * non fullscreen windows, or when running on a Windows-8 + setup with broken DWM switching.
 */
static void OptimizeDWMParameters(PsychWindowRecordType *windowRecord)
{
    DWM_PRESENT_PARAMETERS dwmPresentParams;

    if (!PsychOSIsDWMEnabled(0)) return;
    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-DEBUG: Optimizing windows DWM present parameters. Using minimum queue length of 2 buffers.\n");

    memset(&dwmPresentParams, 0, sizeof(dwmPresentParams));
    dwmPresentParams.cbSize = sizeof(dwmPresentParams);
    dwmPresentParams.cBuffer = 2;

    // Call function if DWM is supported and enabled: We don't check its return code, as this
    // function is only really implemented on Windows-Vista, Windows-7 and Windows-8. It is
    // deprecated as of Windows-8 and support for it will be removed in Windows-9, ie., the
    // function will turn into a no-op and always return error code E_NOTSUP -- "Not supported".
    // See: http://msdn.microsoft.com/en-us/library/windows/desktop/hh994465%28v=vs.85%29.aspx
    // "Queued present model is being deprecated" on MSDN.
    PsychDwmSetPresentParameters(windowRecord->targetSpecific.windowHandle, &dwmPresentParams);

    return;
}

/*
    PsychOSOpenOnscreenWindow()
    
    Creates the pixel format and the context objects and then instantiates the context onto the screen.
    
    -The pixel format and the context are stored in the target specific field of the window recored.  Close
    should clean up by destroying both the pixel format and the context.
    
    -We mantain the context because it must be be made the current context by drawing functions to draw into 
    the specified window.
    
    -We maintain the pixel format object because there seems to be now way to retrieve that from the context.
    
    -To tell the caller to clean up PsychOSOpenOnscreenWindow returns FALSE if we fail to open the window. It 
    would be better to just issue an PsychErrorExit() and have that clean up everything allocated outside of
    PsychOpenOnscreenWindow().
*/
psych_bool PsychOSOpenOnscreenWindow(PsychScreenSettingsType *screenSettings, PsychWindowRecordType *windowRecord, int numBuffers, int stereomode, int conserveVRAM)
{
  char winname[100];
  RECT winRec;
  PsychRectType             screenrect;
  CGDirectDisplayID         cgDisplayID;
  int         pf;
  unsigned int nNumFormats;
  HDC         hDC;
  HWND        hWnd;
  WNDCLASS    wc;
  PIXELFORMATDESCRIPTOR pfd;
  int         attribs[58];
  int         attribcount;
  float       fattribs[2]={0,0};
  int x, y, width, height, i, bpc;
  int		  windowLevel;
  GLenum      glerr;
  DWORD		  flags;
  BOOL        compositorEnabled, compositorPostEnabled;
  
  psych_bool fullscreen = FALSE;
  DWORD windowStyle = WS_CLIPSIBLINGS | WS_CLIPCHILDREN;

  // The WS_EX_NOACTIVATE flag prevents the window from grabbing keyboard focus. That way,
  // the new Java-GetChar can do its job. Can't be used with GUI windows though:
  DWORD windowExtendedStyle = WS_EX_APPWINDOW | ((windowRecord->specialflags & kPsychGUIWindow) ? 0 : WS_EX_NOACTIVATE);

  if (PsychPrefStateGet_Verbosity()>6) {
    printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Entering Win32 specific window setup...\n");
    fflush(NULL);
  }

  // Retrieve windowLevel, an indicator of where non-fullscreen windows should
  // be located wrt. to other windows. 0 = Behind everything else, occluded by
  // everything else. 1 - 999 = At layer windowLevel -> Occludes stuff on layers "below" it.
  // 1000 - 1999 = At highest level, but partially translucent / alpha channel allows to make
  // regions transparent. 2000 or higher: Above everything, fully opaque, occludes everything.
  // 2000 is the default.
  windowLevel = PsychPrefStateGet_WindowShieldingLevel();
  
  // Shielding levels below 1500 will let mouse event through to underlying windows, i.e.,
  // the window is non-existent for the mouse:
  if ((windowLevel >= 1000) && (windowLevel < 1500)) windowExtendedStyle = windowExtendedStyle | WS_EX_TRANSPARENT;

    // Init to safe default:
    windowRecord->targetSpecific.glusercontextObject = NULL;
    windowRecord->targetSpecific.glswapcontextObject = NULL;
    
    // Map the logical screen number to a device handle for the corresponding
    // physical display device: CGDirectDisplayID is currently typedef'd to a
    // HDC windows hardware device context handle.
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenSettings->screenNumber);

    // The compositing desktop window manager (DWM) of MS-Vista and later needs special
    // treatment. Under normal operating conditions we want to disable it, because it
    // can severely interfere with display timing. In certain conditions we want to
    // keep it enabled. In any case we need access to the functions for query and setup
    // of DWM state. As we can't know beforehand if we're running on Vista, we need to
    // dynamically load and link the relevant functions:
    
    // First time invocation? hInstance is zero at this point, if so. We only execute
    // the whole dynamic library detection code once for performance reasons:
    if (!hInstance) {
        // First time. Check if we can load the dwmapi.dll:
		if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: Checking for DWM desktop compositor support... "); 
		
        dwmSupported = FALSE;
        dwmlibrary = LoadLibrary("dwmapi.dll");
        if (dwmlibrary) {
            // Load success. Dynamically bind the relevant functions:
			if (PsychPrefStateGet_Verbosity() > 5) printf(" ... DWM available on this Vista (or later) system. Binding controls ..."); 
            PsychDwmIsCompositionEnabled = (DwmIsCompositionEnabledPROC) GetProcAddress(dwmlibrary, "DwmIsCompositionEnabled");
            PsychDwmEnableComposition    = (DwmEnableCompositionPROC) GetProcAddress(dwmlibrary, "DwmEnableComposition");
            PsychDwmEnableMMCSS          = (DwmEnableMMCSSPROC) GetProcAddress(dwmlibrary, "DwmEnableMMCSS");
            PsychDwmGetCompositionTimingInfo = (DwmGetCompositionTimingInfoPROC) GetProcAddress(dwmlibrary, "DwmGetCompositionTimingInfo");
            PsychDwmSetPresentParameters = (DwmSetPresentParametersPROC) GetProcAddress(dwmlibrary, "DwmSetPresentParameters");
            PsychDwmFlush                = (DwmFlushPROC) GetProcAddress(dwmlibrary, "DwmFlush");
			
            if (PsychDwmIsCompositionEnabled && PsychDwmEnableComposition && PsychDwmEnableMMCSS && PsychDwmGetCompositionTimingInfo && PsychDwmSetPresentParameters && PsychDwmFlush) {
                // Mark dwm as supported:
                dwmSupported = TRUE;
				if (PsychPrefStateGet_Verbosity() > 5) printf(" ...done\n"); 
            }
			else {
				FreeLibrary(dwmlibrary);
				dwmlibrary = 0;
				if (PsychPrefStateGet_Verbosity() > 5) printf(" ...FAILED!\n"); 
				if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not attach to DWM window manager control library dwmapi.dll - Trouble ahead!!\n"); 
			}
        }
		else {
			if (PsychPrefStateGet_Verbosity() > 5) printf(" ... DWM desktop compositing window manager unsupported. Running on a Pre-Vista system.\n"); 
		}
    }
    
    // Check if this should be a fullscreen window:
    PsychGetGlobalScreenRect(screenSettings->screenNumber, screenrect);
    if (PsychMatchRect(screenrect, windowRecord->rect) && (windowLevel >= 2000)) {
      // This is supposed to be a fullscreen window with the dimensions of
      // the current display/desktop:
      fullscreen = TRUE;
    }
    else {
      // Window size different from current screen size:
      // A regular desktop window with borders and control icons is requested, e.g., for debugging:
      fullscreen = FALSE;
    }

    // DWM supported?
    if (dwmSupported) {
        // This is Vista, Windows-7, or a later system with DWM compositing window manager.
		
		// Check current enable state:
        if (PsychDwmIsCompositionEnabled(&compositorEnabled)) {
            // Failed to query state: Assume the worst, i.e., compositor on:
            compositorEnabled = TRUE;
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-WARNING: PsychOSOpenOnscreenWindow: Failed to query state of DWM Windows desktop compositor! Assuming it is ON!\n");
            }
        }
        else {
            if (PsychPrefStateGet_Verbosity() > 3) {
                printf("PTB-INFO: DWM desktop compositor is currently %s.\n", (compositorEnabled) ? "enabled" : "disabled");
            }
        }
		
		// Determine if we want the DWM to be enabled or disabled.
		// These are our rules in priority order:
		
		// Disable unconditionally if it is forced off via kPsychDisableAeroWDM conserveVRAM setting:
        if (conserveVRAM & kPsychDisableAeroDWM) {
            // Compositor shall be off:
            compositorEnabled = FALSE;            
			
            if (PsychPrefStateGet_Verbosity() > 2) {
				printf("PTB-INFO: Will disable DWM because the kPsychDisableAeroDWM flag was set via a call to Screen('Preference', 'ConserveVRAM').\n");
			}
        }
        else {
			// Enable unconditionally if it is forced on via kPsychUseAGLCompositorForFullscreenWindows conserveVRAM setting:
			if (conserveVRAM & kPsychUseAGLCompositorForFullscreenWindows) {
				// Compositor shall be on:
				compositorEnabled = TRUE;            
				
				if (PsychPrefStateGet_Verbosity() > 2) {
					printf("PTB-INFO: Will enable DWM because the kPsychUseCompositorForFullscreenWindows flag was set via a call to Screen('Preference', 'ConserveVRAM').\n");
				}
			}
			else {
				// Disable if a fullscreen window at windowShieldingLevel >= 2000 (i.e., standard, non-transparent, on top of everything else)
				// is requested:
				if (fullscreen && (windowLevel >= 2000)) {
					compositorEnabled = FALSE;            
					
					if (PsychPrefStateGet_Verbosity() > 2) {
						printf("PTB-INFO: Will disable DWM because a regular fullscreen onscreen window is opened -> We want best timing and performance.\n");
					}
				}
				else {
					// Enable if a window with windowShieldingLevel < 2000, i.e. a transparent debug window, is opened:
					if (windowLevel < 2000) {
						compositorEnabled = TRUE;
						
						if (PsychPrefStateGet_Verbosity() > 2) {
							printf("PTB-INFO: Will enable DWM because a special transparent onscreen window in debug mode is opened -> Reduced timing precision and performance!\n");
						}
					}
					else {
						// Regular "windowed" non-fullscreen GUI window and no special conditions or overrides.
						// We don't do anything wrt. DWM, just leave it as it is:
						goto dwmdontcare;
					}
				}
			}
        }

        // Set new compositor state:
        if (PsychDwmEnableComposition((compositorEnabled) ? 1 : 0)) {
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-WARNING: PsychOSOpenOnscreenWindow: Failed to change state of Windows desktop compositor! Expect timing or performance problems!!\n");
            }
        }

		// Retest state:
        if (PsychDwmIsCompositionEnabled(&compositorPostEnabled)) {
            // Failed to query state: Assume the best, i.e., compositor on:
            compositorPostEnabled = TRUE;
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-WARNING: PsychOSOpenOnscreenWindow: Failed to query state of Windows desktop compositor! Assuming it is ON!\n");
            }
        }
        else {
            if (PsychPrefStateGet_Verbosity() > 3) {
                printf("PTB-INFO: DWM desktop compositor is now %s.\n", (compositorPostEnabled) ? "enabled" : "disabled");
            }
        }
		
        // On Windows-8 and later, disabling the DWM is no longer possible, so we need to treat such
        // a failure to switch state as a benign problem and leave it to later code to decide if the
        // user should be notified of potential serious problems or not. Therefore we use relatively
        // modest language here for warnings:
		if (compositorPostEnabled != compositorEnabled) {
            if (PsychPrefStateGet_Verbosity() > 3) {
                printf("PTB-INFO: PsychOSOpenOnscreenWindow: Windows desktop compositor is not %s as requested!\n", (compositorEnabled) ? "enabled" : "disabled");
                if (compositorPostEnabled) printf("PTB-WARNING: PsychOSOpenOnscreenWindow: This may cause reduced or wrong stimulus onset timing and timestamping precision\n");
            }
		}
        
        // Switch compositor to MMCSS scheduling for good timing, if compositor shall be enabled:
        if (compositorPostEnabled && (PsychDwmEnableMMCSS(compositorPostEnabled))) {
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-WARNING: PsychOSOpenOnscreenWindow: Failed to switch Windows desktop compositor to realtime scheduling! Expect timing and performance problems!!\n");
            }
        }		
    }

// goto target jump label for skipping DWM state changes in the "Don't care" case:
dwmdontcare:

    // Wanna fullscreen?
    if (fullscreen) {
		// Switch system to fullscreen-mode without changing any settings:
		
		// Special case "dualdisplay window"? If number of screens is greater than 2, then we've
		// got a multi-display (at least dual-display) setup and screenNumber zero means to open
		// a window that covers both displays 1 and 2 in fullscreen mode:
		if (PsychGetNumDisplays()>2 && screenSettings->screenNumber == 0) {
			// Dual-Display fullscreen mode: Switch both display heads 1 and 2 into fullscreen mode:
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Dual-display operation on logical screen zero --> Switching both screens 1 and 2 into fullscreen mode.\n");
			fullscreen = ChangeScreenResolution(1, 0, 0, 0, 0) && ChangeScreenResolution(2, 0, 0, 0, 0);
		}
		else {
			// Single display case: Just change resolution on that screen:
			fullscreen = ChangeScreenResolution(screenSettings->screenNumber, 0, 0, 0, 0);
		}
    }

    // Special case for explicit multi-display setup under Windows when opening a window on
    // screen zero. We enforce the fullscreen - flag, aka a borderless top level window. This way,
    // if anything of our automatic full-desktop window emulation code goes wrong on exotic setups,
    // the user can still enforce a suitably positioned and sized borderless toplevel window.
    if (PsychGetNumDisplays()>2 && screenSettings->screenNumber == 0) fullscreen = TRUE;

    if (fullscreen) {
        windowStyle |= WS_POPUP;					// Set The WindowStyle To WS_POPUP (Popup Window without borders)
        windowExtendedStyle |= WS_EX_TOPMOST;		// Set The Extended Window Style To WS_EX_TOPMOST

        // Copy absolute screen location and area of window to 'globalrect',
        // so functions like Screen('GlobalRect') can still query the real
        // bounding gox of a window onscreen:
        PsychGetGlobalScreenRect(screenSettings->screenNumber, windowRecord->globalrect);

        // Mark this window as fullscreen window:
        windowRecord->specialflags |= kPsychIsFullscreenWindow;
    }
    else {
        // Only GUI windows have decorations. Non-GUI windows are border/decorationless:
        if (!windowRecord->specialflags & kPsychGUIWindow) {
            // Decorationless, borderless window:
            windowStyle |= WS_POPUP;
        }
        else {
            // GUI window: Needs title-bar, borders, resize handles, the whole bling:
            windowStyle |= WS_OVERLAPPEDWINDOW;
        }

        // Set The Extended Window Style To WS_EX_TOPMOST, ie., this window is in front of all other
        // windows all the time, unless windowLevel is smaller than 1000:
        if (windowLevel >= 1000) windowExtendedStyle |= WS_EX_TOPMOST;

        // If windowLevel is that of a transparent window, then try to enable support for transparent
        // windows:
        // Could also define _WIN32_WINNT >= 0x0500
        #ifndef WS_EX_LAYERED
        #define WS_EX_LAYERED           0x00080000
        #endif
        #ifndef LWA_ALPHA
        #define LWA_ALPHA               2
        #endif
        if ((windowLevel >= 1000) && (windowLevel <  2000)) windowExtendedStyle |= WS_EX_LAYERED;

        // Copy absolute screen location and area of window to 'globalrect',
        // so functions like Screen('GlobalRect') can still query the real
        // bounding gox of a window onscreen:
        PsychCopyRect(windowRecord->globalrect, windowRecord->rect);	  
    }

    // Define final position and size of window:
    x = (int) windowRecord->rect[kPsychLeft];
    y = (int) windowRecord->rect[kPsychTop];
    width  = (int) PsychGetWidthFromRect(windowRecord->rect);
    height = (int) PsychGetHeightFromRect(windowRecord->rect);

    // Assemble windows caption name from window index:
    sprintf(winname, "PTB Onscreen window [%i]:", (int) windowRecord->windowIndex);

	 if (PsychPrefStateGet_Verbosity()>6) {
		 printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Window parameters computed, display switch to fullscreen done (if needed). Registering window class...\n");
		 fflush(NULL);
	 }

	 // Register our own window class for Psychtoolbox onscreen windows:
	 // Only register the window class once - use hInstance as a flag.
	 if (!hInstance) {
		 // Invalidate name and handle of host applications main window:
		 hostwinName[0] = 0;
		 hostwinHandle  = NULL;

		 hInstance = GetModuleHandle(NULL);
		 wc.style         = ((windowLevel >= 1000) && (windowLevel <  2000)) ? 0 : CS_OWNDC;
		 wc.lpfnWndProc   = WndProc;
		 wc.cbClsExtra    = 0;
		 wc.cbWndExtra    = 0;
		 wc.hInstance     = hInstance;
		 wc.hIcon         = LoadIcon(hInstance, IDI_WINLOGO);
		 wc.hCursor       = NULL; // Set window class cursor to NULL. LoadCursor(NULL, IDC_ARROW);
		 wc.hbrBackground = NULL;
		 wc.lpszMenuName  = NULL;
		 wc.lpszClassName = "PTB-OpenGL";
		 
		 if (!RegisterClass(&wc)) {
			 hInstance = 0;
			 printf("\nPTB-ERROR[Register Windowclass failed]: Unknown error, Win32 specific.\n\n");
			 return(FALSE);
		 }
		 
		 // Try to find the window handle of the main window of our runtime environment, i.e.,
		 // the main Matlab GUI window, Matlab console window, or Octave console window. We
		 // use an enumeration and pattern matching procedure to find this. The purpose is to
		 // be able to post WM_CHAR messages received in our own WndProc() handler into the event queue
		 // of the host applications window. This way we can forward characters that were received by
		 // one of our onscreen window (e.g., due to it having the keyboard input focus) to the host
		 // app. Then the CharAvail() and GetChar() functions should be able to dequeue keystrokes from
		 // there and keep working properly if we steal keyboard focus:
		 if (PsychPrefStateGet_Verbosity() > 15) printf("PTB-DEBUG: Window enumeration running...\n");
		 for (i = 1; (hostwinHandle == NULL) && (i <= 2); i++) EnumWindows(PsychHostWindowEnumFunc, i);
		 if (hostwinHandle) {
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Window enumeration done. Our hostwindow is HWND=%p, Name: '%s'\n\n", hostwinHandle, hostwinName);
            
            // Find thread that created GUI windows aka hostwinHandle, then attach our own
            // input queues etc. to it, so they share their event processing. This may help
            // with GetChar() et al. functionality on Win-Vista and later, and may help making
            // HideCursor et al. working again. Or at least, it is my last hope before i'm out
            // of ideas...
            if (!AttachThreadInput(GetCurrentThreadId(), GetWindowThreadProcessId(hostwinHandle, NULL), TRUE) && (PsychPrefStateGet_Verbosity() > 1)) {
                printf("PTB-WARNING: Attaching our thread input to Matlab thread input failed. Expect GetChar() and Show/HideCursor() trouble.\n");
            }
		}
		else {
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-WARNING: Host application window enumeration failed! This may cause trouble with CharAvail() and GetChar() :-(\n\n");
		}
	 }

    // Adjust window bounds to account for the window borders if we are in non-fullscreen mode:
    if (!fullscreen) {
      winRec.left=x; winRec.top=y; winRec.right=x+width; winRec.bottom=y+height;
      AdjustWindowRectEx(&winRec, windowStyle, 0, windowExtendedStyle);
      x=winRec.left; y=winRec.top; width=winRec.right - winRec.left; height=winRec.bottom - winRec.top;
    }

	 if (PsychPrefStateGet_Verbosity()>6) {
		 printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Window class registered - Creating GDI window...\n");
		 fflush(NULL);
	 }

    // Window class registered: Create a window of this class with some specific properties:
    hWnd = CreateWindowEx(windowExtendedStyle,
			  "PTB-OpenGL",
			  winname,
			  windowStyle,
			  x, y, width, height, NULL, NULL, hInstance, NULL);

    if (hWnd == NULL) {
        printf("\nPTB-ERROR[CreateWindow() failed]: Unknown error, Win32 specific.\n\n");
        return(FALSE);
    }

	if (PsychPrefStateGet_Verbosity() > 4) {
		printf("PTB-DEBUG: Created onscreen window has position %i x %i and a size of %i x %i.\n", x, y, width, height);
	}
	
	// Setup transparency level for eligible non-fullscreen windows:
	if (!fullscreen && (windowLevel >= 1000) && (windowLevel < 2000)) {
		// For windowLevels between 1000 and 1999, make the window background transparent, so standard GUI
		// would be visible, wherever nothing is drawn, i.e., where alpha channel is zero.
		// Levels 1000 - 1499 and 1500 to 1999 map to a master opacity level of 0.0 - 1.0:
		SetLayeredWindowAttributes(hWnd, 0, (BYTE) ((((float) (windowLevel % 500)) / 499.0) * 255 + 0.5), LWA_ALPHA);
	}

    // Retrieve device context for the window:
    hDC = GetDC(hWnd);

	 if (PsychPrefStateGet_Verbosity()>6) {
		 printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Window created - Pixelformat selection...\n");
		 fflush(NULL);
	 }

    // Setup optional flags for pixelformat:
    flags = 0;
    // Init pfd to zero:
    memset(&pfd, 0, sizeof(pfd));
    attribcount = 0;

    attribs[attribcount++]=0x2001; // WGL_DRAW_TO_WINDOW_ARB
    attribs[attribcount++]=GL_TRUE;
    attribs[attribcount++]=0x2010; // WGL_SUPPORT_OPENGL_ARB
    attribs[attribcount++]=GL_TRUE;
    attribs[attribcount++]=0x2007; // WGL_SWAP_METHOD_ARB
    attribs[attribcount++]=0x2028; // WGL_SWAP_EXCHANGE_ARB
    attribs[attribcount++]=0x2013; // WGL_PIXEL_TYPE_ARB
    
    // Select either floating point or fixed point framebuffer:
    if (windowRecord->depth == 64 || windowRecord->depth == 128) {
      // Request a floating point drawable instead of a fixed-point one:
      attribs[attribcount++]=WGL_TYPE_RGBA_FLOAT_ARB;
    }
    else {
      // Request standard fixed point drawable:
      attribs[attribcount++]=0x202B; // WGL_TYPE_RGBA_ARB
	}
    
    // Select requested depth per color component 'bpc' for each channel:
    bpc = 8; // We default to 8 bpc == RGBA8
    if (windowRecord->depth == 30)  { bpc = 10; printf("PTB-INFO: Trying to enable at least 10 bpc fixed point framebuffer.\n"); }
    if (windowRecord->depth == 64)  { bpc = 16; printf("PTB-INFO: Trying to enable 16 bpc fixed point framebuffer.\n"); }
    if (windowRecord->depth == 128) { bpc = 32; printf("PTB-INFO: Trying to enable 32 bpc fixed point framebuffer.\n"); }
    
    // Set up color depth for each channel:
    attribs[attribcount++]=WGL_RED_BITS_ARB;
    attribs[attribcount++]=bpc;
    attribs[attribcount++]=WGL_GREEN_BITS_ARB;
    attribs[attribcount++]=bpc;
    attribs[attribcount++]=WGL_BLUE_BITS_ARB;
    attribs[attribcount++]=bpc;
    attribs[attribcount++]=WGL_ALPHA_BITS_ARB;
    // Alpha channel has only 2 bpc in the fixed point bpc=10 case, i.e. RGBA=1010102.
    attribs[attribcount++]=(bpc == 10) ? 2 : bpc;
    
    
    // Stereo display support: If stereo display output is requested with OpenGL native stereo,
    // we request a stereo-enabled rendering context.
    if(stereomode==kPsychOpenGLStereo) {
      flags = flags | PFD_STEREO;
      attribs[attribcount++]=0x2012; // WGL_STEREO_ARB
      attribs[attribcount++]=GL_TRUE;
    }
    
    // Double buffering requested?
    if(numBuffers>=2) {
      // Enable double-buffering:
      flags = flags | PFD_DOUBLEBUFFER;
      attribs[attribcount++]=0x2011; // WGL_DOUBLE_BUFFER_ARB
      attribs[attribcount++]=GL_TRUE;
      
      // AUX buffers for Flip-Operations needed?
      if ((conserveVRAM & kPsychDisableAUXBuffers) == 0) {
	// Allocate one or two (mono vs. stereo) AUX buffers for new "don't clear" mode of Screen('Flip'):
	// Not clearing the framebuffer after "Flip" is implemented by storing a backup-copy of
	// the backbuffer to AUXs before flip and restoring the content from AUXs after flip.
	pfd.cAuxBuffers=(stereomode==kPsychOpenGLStereo || stereomode==kPsychCompressedTLBRStereo || stereomode==kPsychCompressedTRBLStereo) ? 2 : 1;
	attribs[attribcount++]=0x2024; // WGL_AUX_BUFFERS_ARB
	attribs[attribcount++]=pfd.cAuxBuffers;
      }
    }
    
    //if (PSYCH_DEBUG == PSYCH_ON) printf("Device context is %p\n", hDC);
    
    // Build pixelformat descriptor:
    pfd.nSize        = sizeof(pfd);
    pfd.nVersion     = 1;
    pfd.dwFlags      = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_SWAP_EXCHANGE |flags;  // Want OpenGL capable window with bufferswap via page-flipping...
    pfd.iPixelType   = PFD_TYPE_RGBA; // Want a RGBA pixel format.
    pfd.cColorBits   = 32;            // 32 bpp at least...
    pfd.cAlphaBits   = (bpc == 10) ? 2 : 8;	// Want a 8 bit alpha-buffer, unless R10G10B10A2 pixelformat requested for native 10 bpc support.

    // Support for OpenGL 3D rendering requested?
    if (PsychPrefStateGet_3DGfx()) {
      // Yes. Allocate and attach a 24bit depth buffer and 8 bit stencil buffer:
      pfd.cDepthBits = 24;
      pfd.cStencilBits = 8;
      attribs[attribcount++]=0x2022; // WGL_DEPTH_BITS_ARB
      attribs[attribcount++]=24;
      attribs[attribcount++]=0x2023; // WGL_STENCIL_BITS_ARB
      attribs[attribcount++]=8;
	  // Alloc an accumulation buffer as well?
	  if (PsychPrefStateGet_3DGfx() & 2) {
		  // Yes: Alloc accum buffer, request 64 bpp, aka 16 bits integer per color component if possible:
		  attribs[attribcount++] = WGL_ACCUM_BITS_EXT;
		  attribs[attribcount++] = 64;
		  attribs[attribcount++] = WGL_ACCUM_RED_BITS_EXT;
		  attribs[attribcount++] = 16;
		  attribs[attribcount++] = WGL_ACCUM_GREEN_BITS_EXT;
		  attribs[attribcount++] = 16;
		  attribs[attribcount++] = WGL_ACCUM_BLUE_BITS_EXT;
		  attribs[attribcount++] = 16;
		  attribs[attribcount++] = WGL_ACCUM_ALPHA_BITS_EXT;
		  attribs[attribcount++] = 16;
		  pfd.cAccumBits      = 64;
		  pfd.cAccumRedBits   = 16;
		  pfd.cAccumGreenBits = 16;
		  pfd.cAccumBlueBits  = 16;
		  pfd.cAccumAlphaBits = 16;
	  }	  
    }
    
    // Multisampled Anti-Aliasing requested?
    if (windowRecord->multiSample > 0) {
      // Request a multisample buffer:
      attribs[attribcount++]= 0x2041; // WGL_SAMPLE_BUFFERS_ARB
      attribs[attribcount++]= 1;
      // Request at least multiSample samples per pixel:
      attribs[attribcount++]= 0x2042; // WGL_SAMPLES_ARB
      attribs[attribcount++]= windowRecord->multiSample;
    }
    
    // Finalize attribs-array:
    attribs[attribcount++]= 0;
    
    //if (PSYCH_DEBUG == PSYCH_ON) printf("Choosing pixelformat\n");
    
    // Create pixelformat:
    // This is typical Microsoft brain-damage: We first need to create a window the
    // conventional old way just to be able to get a handle to the new wglChoosePixelFormat
    // method, which will us - after destroying and recreating the new window - allow to
    // select the pixelformat we actually want!
    
    // Step 1: Choose pixelformat old-style:
    pf = ChoosePixelFormat(hDC, &pfd);
    
    // Do we have a valid pixelformat?
    if (pf == 0) {
      // Nope. We give up!
      ReleaseDC(hWnd, hDC);
      DestroyWindow(hWnd);      
      printf("\nPTB-ERROR[ChoosePixelFormat() failed]: Unknown error, Win32 specific.\n\n");
      return(FALSE);
    }
    
    // Yes. Set it:
    if (SetPixelFormat(hDC, pf, &pfd) == FALSE) {
      ReleaseDC(hWnd, hDC);
      DestroyWindow(hWnd);      
      
      printf("\nPTB-ERROR[SetPixelFormat() failed]: Unknown error, Win32 specific.\n\n");
      return(FALSE);
    }
    
	 if (PsychPrefStateGet_Verbosity()>6) {
		 printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: ChoosePixelFormat(), SetPixelFormat() done... Creating OpenGL context...\n");
		 fflush(NULL);
	 }

    // Ok, create and attach the rendering context.
    windowRecord->targetSpecific.contextObject = wglCreateContext(hDC);
    if (windowRecord->targetSpecific.contextObject == NULL) {
      ReleaseDC(hWnd, hDC);
      DestroyWindow(hWnd);
      
      printf("\nPTB-ERROR:[Context creation failed] Unknown, Win32 specific.\n\n");
      return(FALSE);
    }
    
    // Store the handles...
    windowRecord->targetSpecific.windowHandle = hWnd;
    windowRecord->targetSpecific.deviceContext = hDC;
    
	 if (PsychPrefStateGet_Verbosity()>6) {
		 printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Context created - Activating and binding context...\n");
		 fflush(NULL);
	 }

    // Activate the rendering context:
    PsychOSSetGLContext(windowRecord);
    
	 if (PsychPrefStateGet_Verbosity()>6) {
		 printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Online - glewInit()...\n");
		 fflush(NULL);
	 }

    // Ok, the OpenGL rendering context is up and running. Auto-detect and bind all
    // available OpenGL extensions via GLEW:
    glerr = glewInit();
    if (GLEW_OK != glerr)
      {
	/* Problem: glewInit failed, something is seriously wrong. */
	printf("\nPTB-ERROR[GLEW init failed: %s]: Please report this to the forum. Will try to continue, but may crash soon!\n\n", glewGetErrorString(glerr));
	fflush(NULL);
      }
    else {
      if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: Using GLEW version %s for automatic detection of OpenGL extensions...\n", glewGetString(GLEW_VERSION));
    }
    
    DescribePixelFormat(hDC, pf, sizeof(PIXELFORMATDESCRIPTOR), &pfd);
    
	// Check for pageflipping for bufferswaps and output warning if we don't get it:
    if (!(pfd.dwFlags & PFD_SWAP_EXCHANGE) && (pfd.dwFlags & PFD_SWAP_COPY) && (PsychPrefStateGet_Verbosity() > 1)) {
        printf("PTB-WARNING: Created onscreen window on screenid %i will probably not be able to use GPU pageflipping for\n", screenSettings->screenNumber);
        printf("PTB-WARNING: Screen('Flip')! May cause tearing artifacts and unreliable or wrong visual stimulus onset timestamping!\n");
    }
    
	// Special debug override for faulty drivers with non-working extension:
	if (conserveVRAM & kPsychOverrideWglChoosePixelformat) wglChoosePixelFormatARB = NULL;

	// Step 2: Ok, we have an OpenGL rendering context with all known extensions bound:
	// Do we have (or want) support for wglChoosePixelFormatARB?
   // We skip use of it if we can do without it, i.e., when we don't need unusual framebuffer
   // configs (like non 8bpc fixed) and we don't need multisampling. This is a work-around
   // for hardware that has trouble (=driver bugs) with wglChoosePixelformat() in some modes, e.g., the NVidia
	// Quadro gfx card, which fails to enable quad-buffered stereo when using wglChoosePixelformat(),
	// but does so perfectly well when using the old-style ChoosePixelFormat(). 
	if ((wglChoosePixelFormatARB == NULL) || ((bpc==8) && (windowRecord->multiSample <= 0))) {
	  // Failed (==NULL) or disabled via override.
	  if ((wglChoosePixelFormatARB == NULL) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: Could not bind wglChoosePixelFormat - Extension. Some features will be unavailable, e.g., Anti-Aliasing and high precision framebuffers.\n");
	}
	else {
	  // Supported. We destroy the rendering context and window, then recreate it with
	  // the wglChoosePixelFormat - method...

	 if (PsychPrefStateGet_Verbosity()>6) {
		 printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Need to use wglChoosePixelFormat() path due to multisampling or bpc > 8 enabled. Destroy and reinit sequence...\n");
		 fflush(NULL);
	 }

	  wglMakeCurrent(NULL, NULL);

	  // Delete rendering context:
	  wglDeleteContext(windowRecord->targetSpecific.contextObject);
	  windowRecord->targetSpecific.contextObject=NULL;

      // Release device context:
      ReleaseDC(windowRecord->targetSpecific.windowHandle, windowRecord->targetSpecific.deviceContext);
      windowRecord->targetSpecific.deviceContext=NULL;

      // Close & Destroy the window:
      DestroyWindow(windowRecord->targetSpecific.windowHandle);
      windowRecord->targetSpecific.windowHandle=NULL;

		// Ok, old window and stuff is dead. Create new window:
    	hWnd = CreateWindowEx(windowExtendedStyle, "PTB-OpenGL", winname, windowStyle,
			  						 x, y, width, height, NULL, NULL, hInstance, NULL);
    	if (hWnd == NULL) {
        printf("\nPTB-ERROR[CreateWindow() - II failed]: Unknown error, Win32 specific.\n\n");
        return(FALSE);
    	}

		// Setup transparency level for eligible non-fullscreen windows:
		if (!fullscreen && (windowLevel >= 1000) && (windowLevel < 2000)) {
			// For windowLevels between 1000 and 1999, make the window background transparent, so standard GUI
			// would be visible, wherever nothing is drawn, i.e., where alpha channel is zero.
			// Levels 1000 - 1499 and 1500 to 1999 map to a master opacity level of 0.0 - 1.0:
			SetLayeredWindowAttributes(hWnd, 0, (BYTE) ((((float) (windowLevel % 500)) / 499.0) * 255 + 0.5), LWA_ALPHA);
		}
		
	   // Retrieve device context for the window:
    	hDC = GetDC(hWnd);

		if (PsychPrefStateGet_Verbosity()>6) {
			printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Context, pixelformat, window destroyed, window recreated - Now wglChoosePixelformat()...\n");
			fflush(NULL);
		}
		
		pf = 0;
		nNumFormats=0;
		wglChoosePixelFormatARB(hDC, &attribs[0], NULL, 1, &pf, &nNumFormats);
      if (nNumFormats==0 && windowRecord->multiSample > 0) {
		 	// Failed. Probably due to too demanding multisample requirements: Lets lower them...
			for (i=0; i<attribcount && attribs[i]!=0x2042; i++);
			// Reduce requested number of samples/pixel and retry until we get one:
			while (nNumFormats==0 && windowRecord->multiSample > 0) {
				// printf("Failed for multisampling level %i nNum=%i\n", windowRecord->multiSample, nNumFormats);
				attribs[i+1]--;
	  			windowRecord->multiSample--;
				wglChoosePixelFormatARB(hDC, &attribs[0], NULL, 1, &pf, &nNumFormats);
			}
			// If we still do not get one with 0 samples per pixel, we can try to disable
			// multisampling completely:
			if (windowRecord->multiSample == 0 && nNumFormats==0) {
				// printf("Failed for multisampling level %i nNum=%i --> Disabling multisampling...\n", windowRecord->multiSample, nNumFormats);
				// We 0-Out all entries related to multi-sampling, including the
				// GL_SAMPLES_ARB and GL_SAMPLE_BUFFERS_ARB enums themselves:
				for (i=0; i<attribcount && attribs[i]!=0x2042; i++);
	  			attribs[i]=0;
	  			for (i=0; i<attribcount && attribs[i]!=0x2041; i++);
	  			attribs[i]=0;
	  			attribs[i+1]=0;
				wglChoosePixelFormatARB(hDC, &attribs[0], NULL, 1, &pf, &nNumFormats);
			}
      }

		if (nNumFormats==0 && numBuffers>=2) {
			// We still don't have a valid pixelformat, but double-buffering is enabled.
			// Let's try if we get one if we do not request any AUX-Buffers:
			for (i=0; i<attribcount && attribs[i]!=0x2024; i++);
			attribs[i+1] = 0; // Zero AUX-Buffers requested.
			wglChoosePixelFormatARB(hDC, &attribs[0], NULL, 1, &pf, &nNumFormats);
		}

		if (nNumFormats==0 && numBuffers>=2) {
			// We still don't have a valid pixelformat, but double-buffering is enabled.
			// Let's try if we get one if we do not request SWAP_EXCHANGED double buffering anymore:
			for (i=0; i<attribcount && attribs[i]!=0x2028; i++);
			attribs[i] = 0x202A; // WGL_SWAP_UNDEFINED_ARB
			wglChoosePixelFormatARB(hDC, &attribs[0], NULL, 1, &pf, &nNumFormats);
		}

		// Either we have a multisampled or non-multisampled format, or none. If we failed,
		// then we can not do anything anymore about it.

	   // Do we have a valid pixelformat?
      if (nNumFormats == 0) {
		   // Nope. We give up!
         ReleaseDC(hWnd, hDC);
         DestroyWindow(hWnd);      
         printf("\nPTB-ERROR[wglChoosePixelFormat() failed]: Unknown error, Win32 specific. Code: %i.\n\n", GetLastError());
         return(FALSE);
      }

      // Yes. Set it:
      if (SetPixelFormat(hDC, pf, &pfd) == FALSE) {
         ReleaseDC(hWnd, hDC);
         DestroyWindow(hWnd);      

         printf("\nPTB-ERROR[SetPixelFormat() - II failed]: Unknown error, Win32 specific.\n\n");
         return(FALSE);
      }

		if (PsychPrefStateGet_Verbosity()>6) {
			printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: New pixelformat chosen and set - About to recreate master opengl context...\n");
			fflush(NULL);
		}

      // Ok, create and attach the rendering context.
      windowRecord->targetSpecific.contextObject = wglCreateContext(hDC);
      if (windowRecord->targetSpecific.contextObject == NULL) {
         ReleaseDC(hWnd, hDC);
         DestroyWindow(hWnd);

         printf("\nPTB-ERROR:[Context creation II failed] Unknown, Win32 specific.\n\n");
         return(FALSE);
      }
    
      // Store the handles...
      windowRecord->targetSpecific.windowHandle = hWnd;
      windowRecord->targetSpecific.deviceContext = hDC;

		if (PsychPrefStateGet_Verbosity()>6) {
			printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Context recreated - Activating...\n");
			fflush(NULL);
		}

      // Activate the rendering context:
      PsychOSSetGLContext(windowRecord);

      DescribePixelFormat(hDC, pf, sizeof(PIXELFORMATDESCRIPTOR), &pfd);

		// Done with final window and OpenGL context setup. We've got our final context enabled.
	 }

	if (PsychPrefStateGet_Verbosity()>6) {
		printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Mastercontext created, activated and bound - Enabling multisampling if needed...\n");
		fflush(NULL);
	}
	
	 // Enable multisampling if this was requested:
    if (windowRecord->multiSample > 0) glEnable(0x809D); // 0x809D == GL_MULTISAMPLE_ARB
	 // Throw away any error-state this could have created on old hardware...

	 while((glerr = glGetError())!=GL_NO_ERROR) {
		 if (PsychPrefStateGet_Verbosity()>6) {
			 printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: After multisample enable - if any - glGetError reports error: %s..\n", gluErrorString(glerr));
			 fflush(NULL);
		 }		 
	 }

	 if (PsychPrefStateGet_Verbosity()>4) {
		 printf("PTB-DEBUG: Window and master OpenGL context creation finished.\n");
		 fflush(NULL);
	 }
     
     while((glerr = glGetError())!=GL_NO_ERROR) {
         if (PsychPrefStateGet_Verbosity()>6) {
             printf("PTB-DEBUG: Before slaveWindow context sharing: glGetError reports error: %s..\n", gluErrorString(glerr));
             fflush(NULL);
         }       
     }
     
     if (PsychPrefStateGet_Verbosity()>6) {
         printf("PTB-DEBUG: Before slaveWindow context sharing: glGetString reports %p pointer...\n", glGetString(GL_EXTENSIONS));
         fflush(NULL);
     }
     
     // Do we have a slaveWindow with which to share all object ressources like display lists, textures, FBO's and shaders?
     if (windowRecord->slaveWindow) {
         // Enable ressource sharing with slaveWindow context for this window:
         if (!wglShareLists(windowRecord->slaveWindow->targetSpecific.contextObject, windowRecord->targetSpecific.contextObject)) {
             // This is ugly, but not fatal...
             if (PsychPrefStateGet_Verbosity()>0) {
                 printf("\nPTB-WARNING[wglShareLists for slaveWindow context failed]: Ressource sharing with OpenGL context for slave window failed for unknown reasons. Dual-Window stereo may not work.\n\n");
             }      
         }
     }
     
     glFinish();
     while((glerr = glGetError())!=GL_NO_ERROR) {
         if (PsychPrefStateGet_Verbosity()>6) {
             printf("PTB-DEBUG: After slaveWindow context sharing: glGetError reports error: %s..\n", gluErrorString(glerr));
             fflush(NULL);
         }       
     }
     
     if (PsychPrefStateGet_Verbosity()>6) {
         printf("PTB-DEBUG: After slaveWindow context sharing: glGetString reports %p pointer...\n", glGetString(GL_EXTENSIONS));
         fflush(NULL);
     }
	 
    // External 3D graphics support enabled? Or OpenGL quad-buffered stereo enabled?
	 // For the former, we need this code for OpenGL state isolation. For the latter we
    // need this code as workaround for Windows brain-damage. For some reason it helps
    // to properly shutdown stereo contexts on Windows...
	// One can disable this mechanism via the flag kPsychDisableContextIsolation to account
	// for even more MS brain-damage.
	 if ((!(conserveVRAM & kPsychDisableContextIsolation)) && (PsychPrefStateGet_3DGfx() || stereomode == kPsychOpenGLStereo)) {
		// Yes. We need to create an extra OpenGL rendering context for the external
		// OpenGL code to provide optimal state-isolation. The context shares all
		// heavyweight ressources likes textures, FBOs, VBOs, PBOs, display lists and
		// starts off as an identical copy of PTB's context as of here.
		
		if (PsychPrefStateGet_Verbosity()>4) {
			printf("PTB-DEBUG: Setup of userspace rendering context...\n");
			fflush(NULL);
		}

		windowRecord->targetSpecific.glusercontextObject = wglCreateContext(hDC);
		if (windowRecord->targetSpecific.glusercontextObject == NULL) {
            ReleaseDC(hWnd, hDC);
            DestroyWindow(hWnd);
            windowRecord->targetSpecific.windowHandle = NULL;
            windowRecord->targetSpecific.deviceContext = NULL;

            printf("\nPTB-ERROR[UserContextCreation failed]: Creating a private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n");
            return(FALSE);
		}

		if (PsychPrefStateGet_Verbosity()>4) {
			printf("PTB-DEBUG: wglCreateContext() done. Now copying context state...\n");
			fflush(NULL);
		}

		if (PsychPrefStateGet_Verbosity()>6) {
			printf("PTB-DEBUG: Userspace context wglCreateContext() done: glGetString reports %p pointer...\n", glGetString(GL_EXTENSIONS));
			fflush(NULL);
		}
		
		glFinish();
		while((glerr = glGetError())!=GL_NO_ERROR) {
			if (PsychPrefStateGet_Verbosity()>6) {
				printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Pre-Switch to usercontext glGetError reports error: %s..\n", gluErrorString(glerr));
				fflush(NULL);
			}		 
		}

		// Special debug override for faulty drivers with context sharing setup bugs (NVidia GF8xxx with drivers from Dec. 2007):
		// Ok, we unconditionally disable this code-path, regardless what kPsychUseWindowsContextSharingWorkaround1 is.
		// The wglCopyContext() call isn't needed for proper operation, but causes crashes on recent drivers. By avoiding
		// it completely we don't lose functionality, but avoid noise on the forum about non-working stuff...
/*
		if (!(conserveVRAM & kPsychUseWindowsContextSharingWorkaround1)) {
		
			wglMakeCurrent(windowRecord->targetSpecific.deviceContext, windowRecord->targetSpecific.glusercontextObject);
			
			while((glerr = glGetError())!=GL_NO_ERROR) {
				if (PsychPrefStateGet_Verbosity()>6) {
					printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Post-Switch to usercontext glGetError reports error: %s..\n", gluErrorString(glerr));
					fflush(NULL);
				}		 
			}
			
			if (PsychPrefStateGet_Verbosity()>6) {
				printf("PTB-DEBUG: In Userspace context: glGetString reports %p pointer...\n", glGetString(GL_EXTENSIONS));
				fflush(NULL);
			}
			
			glFinish();
			wglMakeCurrent(windowRecord->targetSpecific.deviceContext, windowRecord->targetSpecific.contextObject);
			
			glFinish();
			while((glerr = glGetError())!=GL_NO_ERROR) {
				if (PsychPrefStateGet_Verbosity()>6) {
					printf("PTB-DEBUG: PsychOSOpenOnscreenWindow: Post-Switch to mastercontext glGetError reports error: %s..\n", gluErrorString(glerr));
					fflush(NULL);
				}		 
			}
			
			if (PsychPrefStateGet_Verbosity()>6) {
				printf("PTB-DEBUG: Back in mastercontext: glGetString reports %p pointer...\n", glGetString(GL_EXTENSIONS));
				fflush(NULL);
			}
			
			// Copy full state from our main context:
			if(!wglCopyContext(windowRecord->targetSpecific.contextObject, windowRecord->targetSpecific.glusercontextObject, GL_ALL_ATTRIB_BITS)) {
				// This is ugly, but not fatal...
				if (PsychPrefStateGet_Verbosity()>1) {
					printf("\nPTB-WARNING[wglCopyContext for user context failed]: Copying state to private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n");
				}
			}
			
			if (PsychPrefStateGet_Verbosity()>6) {
				printf("PTB-DEBUG: After wglCopyContext(): glGetString reports %p pointer...\n", glGetString(GL_EXTENSIONS));
				fflush(NULL);
			}
			
			if (PsychPrefStateGet_Verbosity()>4) {
				printf("PTB-DEBUG: wglCopyContext() done. Now enabling ressource sharing..\n");
				fflush(NULL);
			}
		}
		else {
			if (PsychPrefStateGet_Verbosity()>4) {
				printf("PTB-DEBUG: wglMakeCurrent(usercontext); wglMakeCurrent(mastercontext); wglCopyContext() skipped due to kPsychUseWindowsContextSharingWorkaround1 flag. Now enabling ressource sharing..\n");
				fflush(NULL);
			}
		}
*/
		
	   // Enable ressource sharing with master context for this window:
		if (!wglShareLists(windowRecord->targetSpecific.contextObject, windowRecord->targetSpecific.glusercontextObject)) {
			// This is ugly, but not fatal...
			if (PsychPrefStateGet_Verbosity()>1) {
				printf("\nPTB-WARNING[wglShareLists for user context failed]: Ressource sharing with private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n");
			}		
		}

		if (PsychPrefStateGet_Verbosity()>6) {
			printf("PTB-DEBUG: After wglShareLists(): glGetString reports %p pointer...\n", glGetString(GL_EXTENSIONS));
			fflush(NULL);
		}

		if (PsychPrefStateGet_Verbosity()>4) {
			printf("PTB-DEBUG: Userspace context setup done..\n");
			fflush(NULL);
		}
	 }

	 glFinish();

     // Setup dedicated swap context for async flips:
     windowRecord->targetSpecific.glswapcontextObject = wglCreateContext(hDC);
     if (windowRecord->targetSpecific.glswapcontextObject == NULL) {
         ReleaseDC(hWnd, hDC);
         DestroyWindow(hWnd);
         windowRecord->targetSpecific.windowHandle = NULL;
         windowRecord->targetSpecific.deviceContext = NULL;
         
         printf("\nPTB-ERROR[SwapContextCreation failed]: Creating a private OpenGL context for async flips failed for unknown reasons.\n\n");
         return(FALSE);
     }
     
     // Enable ressource sharing with master context for this context:
     if (!wglShareLists(windowRecord->targetSpecific.contextObject, windowRecord->targetSpecific.glswapcontextObject)) {
         // This is ugly, but not fatal...
         if (PsychPrefStateGet_Verbosity()>1) {
             printf("\nPTB-WARNING[wglShareLists for swap context failed]: Ressource sharing with private OpenGL context for async flips failed for unknown reasons.\n\n");
         }		
     }

	 if (PsychPrefStateGet_Verbosity()>4) {
		printf("PTB-DEBUG: Final low-level window setup: ShowWindow(), SetCapture(), diagnostics...\n");
		fflush(NULL);
	 }

	 if (PsychPrefStateGet_Verbosity()>4) {
		printf("PTB-DEBUG: OpenGL initialization of all master-/slave-/shared-/userspace contexts finished...\n");
		printf("PTB-DEBUG: Final low-level window setup: ShowWindow(), SetCapture(), diagnostics...\n");
		fflush(NULL);
	 }

    // Finally, show our new window:
    if (windowLevel != -1) ShowWindow(hWnd, SW_SHOW);

	// Following settings enforce the onscreen window being and staying the
	// foreground window with keyboard focus. However they interfere with our
	// Java based GetChar implementation, so use of GetChar and of this feature
	// is mutually exclusive. This mode is needed on some degenerated Windows-2000 and
	// Windows-XP systems under some weird conditions. It is also crucially needed on
	// MS-Vista, Windows-7 and later to achieve tear-free double-buffered stimulus onset via
	// page-flipping that is hardware-sychronized to VSYNC. Without this flag, the window is
	// not considered a "real", unoccluded fullscreen window, with the consequence that if
	// the DWM is enabled, it will not auto-suspend but stay in the loop and perform compositing
	// (--> bad performance), and if the DWM is disabled in "Classic mode", the DirectX
	// graphics kernel subsystem will prevent use of VSYNC'ed page-flipping, but revert to
	// VSYNC-IRQ triggered BitBlt's --> Tearing and timing trouble on systems with bad ISR timing!
	//
	// The documentation for the DirectX function "IDXGISwapChain::SetFullscreenState", to be found
	// at http://msdn.microsoft.com/en-us/library/ee421963(VS.85).aspx contains the conditions for
	// "true" fullscreen mode. Swapchains are the underlying DirectX mechanism for OpenGL double-buffers
	// and the bufferswap operations. A SwapBuffers() OpenGL call eventually results in a DirectX
	// Present() call to enqueue the backbuffer in the swapchain which will perform the actual swap.
	// Only swapchains in fullscreen mode will use page-flipping with hardware VSYNC triggers, this is
	// why it is important that our windows satisfy above documented fullscreen conditions.
	//
	//
	// Soo, these are the rules for enabling this mode:
	// 1. If the kPsychPreventForegroundWindow flag is set, we won't enable --> This is a master-off switch.
	// 2. If the kPsychEnforceForegroundWindow flag is set, we will enable unconditionally.
	// 3. If neither of the above is true, we will not enable on WinXP, Win2000 and earlier, but
	//    we will enable on Vista, Windows-7 and later iff a regular, non-transparent (windowLevel >= 2000)
	//    fullscreen window is opened.
	//
	// A special 'Preference'-'ConserveVRAM' flag allows
	// to enable this feature: kPsychPreventForegroundWindow
	if (!(conserveVRAM & kPsychPreventForegroundWindow) && ((conserveVRAM & kPsychEnforceForegroundWindow) || (fullscreen && (windowLevel >= 2000) && PsychIsMSVista()))) {
		// Give it higher priority than other applications windows:
		SetForegroundWindow(hWnd);

		// Set the focus on it:
		SetFocus(hWnd);
		
		// Make sure it is the topmost window:
		SetWindowPos(hWnd, HWND_TOPMOST, 0, 0, 0, 0, (SWP_NOMOVE | SWP_NOSIZE));
		
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Executed SetForegroundWindow() and SetFocus() on window to optimize pageflipping and timing.\n");
	}
	else {
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: No call to SetForegroundWindow() and SetFocus() on window, because trigger-conditions not satisfied.\n");
	}

    // Capture the window if it is a fullscreen one: This window will receive all
    // mouse move and mouse button press events. Important for GetMouse() to work
    // properly...
    if (fullscreen) SetCapture(hWnd);

	// Recheck for pageflipping for bufferswaps and output warning if we don't get it:
	pf = GetPixelFormat(windowRecord->targetSpecific.deviceContext);
    DescribePixelFormat(windowRecord->targetSpecific.deviceContext, pf, sizeof(PIXELFORMATDESCRIPTOR), &pfd);

    if (!(pfd.dwFlags & PFD_SWAP_EXCHANGE) && (pfd.dwFlags & PFD_SWAP_COPY) && (PsychPrefStateGet_Verbosity() > 1)) {
        printf("PTB-WARNING: Created onscreen window on screenid %i will probably not be able to use GPU pageflipping for\n", screenSettings->screenNumber);
        printf("PTB-WARNING: Screen('Flip')! May cause tearing artifacts and unreliable or wrong visual stimulus onset timestamping!\n");
    }
    
    // Increase our own open window counter:
    win32_windowcount++;

    // Some info for the user regarding non-fullscreen windows:
    if (!fullscreen) {
      printf("PTB-INFO: Most graphics cards will not support proper syncing to vertical retrace when\n");
      printf("PTB-INFO: running in windowed mode (non-fullscreen). If PTB aborts with 'Synchronization failure'\n");
      printf("PTB-INFO: you can disable the sync test via call to Screen('Preference', 'SkipSyncTests', 2); .\n");
      printf("PTB-INFO: You won't get proper stimulus onset timestamps though, so windowed mode may be of limited use.\n");
    }

    // Check for the VSYNC extension:
	if (wglSwapIntervalEXT == NULL) {
		wglSwapIntervalEXT = NULL;
		if (PsychPrefStateGet_Verbosity() > 1) {
			printf("PTB-WARNING: Your graphics driver doesn't allow me to control if bufferswaps should be synchronized to the vertical retrace!\n");
			printf("PTB-WARNING: This can cause massive stimulus timing problems, failure of the sync tests and calibrations and severe visual tearing artifacts!\n");
			printf("PTB-WARNING: Please update your display graphics driver as soon as possible to fix this and make sure this functionality is not disabled in\n");
			printf("PTB-WARNING: the display settings control panel of your graphics card.\n");
			printf("PTB-WARNING: If everything else fails, you can usually manually enable synchronization to vertical retrace somewhere in the display settings\n");
			printf("PTB-WARNING: control panel of your machine.\n");
		}
    }

    // If the DWM is enabled, try to optimize its presentation parameters for our purpose:
    OptimizeDWMParameters(windowRecord);
    
	// Enforce a one-shot GUI event queue dispatch via this dummy call to PsychGetMouseButtonState() to
	// make windows GUI event processing happy:
	PsychGetMouseButtonState(NULL);

    // Ok, we should be ready for OS independent setup...
	 if (PsychPrefStateGet_Verbosity()>4) {
		 printf("PTB-DEBUG: Final low-level window setup finished. Continuing with OS-independent setup.\n");
		 fflush(NULL);
	 }

    // Well Done!
    return(TRUE);
}


/*
    PsychOSOpenOffscreenWindow()
    
    Accept specifications for the offscreen window in the platform-neutral structures, convert to native CoreGraphics structures,
    create the surface, allocate a window record and record the window specifications and memory location there.
	
	TO DO:  We need to walk down the screen number and fill in the correct value for the benefit of TexturizeOffscreenWindow
*/
psych_bool PsychOSOpenOffscreenWindow(double *rect, int depth, PsychWindowRecordType **windowRecord)
{
  /*
    //PsychTargetSpecificWindowRecordType 	cgStuff;
    CGLPixelFormatAttribute 			attribs[5];
    //CGLPixelFormatObj					pixelFormatObj;
    long								numVirtualScreens;
    CGLError							error;
    int									windowWidth, windowHeight;
    int									depthBytes;

    //First allocate the window recored to store stuff into.  If we exit with an error PsychErrorExit() should
    //call PsychPurgeInvalidWindows which will clean up the window record. 
    PsychCreateWindowRecord(windowRecord);  		//this also fills the window index field.
    
    attribs[0]=kCGLPFAOffScreen;
    attribs[1]=kCGLPFAColorSize;
    attribs[2]=(CGLPixelFormatAttribute)depth;
    attribs[3]=(CGLPixelFormatAttribute)NULL;
    
    error=CGLChoosePixelFormat(attribs, &((*windowRecord)->targetSpecific.pixelFormatObject), &numVirtualScreens);
    error=CGLCreateContext((*windowRecord)->targetSpecific.pixelFormatObject, NULL, &((*windowRecord)->targetSpecific.contextObject));
	CGLSetCurrentContext((*windowRecord)->targetSpecific.contextObject);
	
    windowWidth=(int)PsychGetWidthFromRect(rect);
    windowHeight=(int) PsychGetHeightFromRect(rect);
	//This section looks wrong because it does not allocate enough memory to insure alignment on word bounaries, which presumably is
	//dicated by the pixel format.  
    depthBytes=depth / 8;
    (*windowRecord)->surfaceSizeBytes= windowWidth * windowHeight * depthBytes;
    (*windowRecord)->surface=malloc((*windowRecord)->surfaceSizeBytes);
    CGLSetOffScreen((*windowRecord)->targetSpecific.contextObject, windowWidth, windowHeight, windowWidth * depthBytes, (*windowRecord)->surface); 
    gluOrtho2D(rect[kPsychLeft], rect[kPsychRight], rect[kPsychBottom], rect[kPsychTop]);
          
    //Fill in the window record.
    (*windowRecord)->windowType=kPsychSystemMemoryOffscreen;
    (*windowRecord)->screenNumber=kPsychUnaffiliatedWindow;
    PsychCopyRect((*windowRecord)->rect, rect);
    (*windowRecord)->depth=depth;
	

    //mark the contents of the window record as valid.  Between the time it is created (always with PsychCreateWindowRecord) and when it is marked valid 
    //(with PsychSetWindowRecordValid) it is a potential victim of PsychPurgeInvalidWindows.  
    PsychSetWindowRecordValid(*windowRecord);
    return(TRUE);
  */

  // FIXME: Not yet implemented.
  return(FALSE);

}


void PsychOSCloseWindow(PsychWindowRecordType *windowRecord)
{
  // Detach rendering context again - just to be safe!
  wglMakeCurrent(NULL, NULL);

  // Delete rendering context:
  wglDeleteContext(windowRecord->targetSpecific.contextObject);
  windowRecord->targetSpecific.contextObject=NULL;

  // Delete swap context:
  wglDeleteContext(windowRecord->targetSpecific.glswapcontextObject);
  windowRecord->targetSpecific.glswapcontextObject=NULL;

  // Delete userspace context:
  if (windowRecord->targetSpecific.glusercontextObject) {
		wglDeleteContext(windowRecord->targetSpecific.glusercontextObject);
  		windowRecord->targetSpecific.glusercontextObject = NULL;
  }

  // Release device context:
  ReleaseDC(windowRecord->targetSpecific.windowHandle, windowRecord->targetSpecific.deviceContext);
  windowRecord->targetSpecific.deviceContext=NULL;

  // Release the capture, whatever that means...
  ReleaseCapture();

  // Close & Destroy the window:
  DestroyWindow(windowRecord->targetSpecific.windowHandle);
  windowRecord->targetSpecific.windowHandle=NULL;

  // Restore video settings from the defaults in the Windows registry:
  ChangeDisplaySettings(NULL, 0);

  // Was this the last window?
  win32_windowcount--;

  if (win32_windowcount<=0) {
    win32_windowcount=0;

    // Unregister our window class if it is still registered:
    if (hInstance) {
      UnregisterClass("PTB-OpenGL", hInstance);
      hInstance=NULL;
      
      // Detach from and release dwmapi.dll if loaded and attached:
      if (dwmSupported && dwmlibrary) {
		  // Reenable DWM if it was disabled (by us or others):
		  if (!PsychOSIsDWMEnabled(0)) {
	          // Enable compositor:
			  if (PsychDwmEnableComposition(1)) {
				  if (PsychPrefStateGet_Verbosity() > 1) {
					  printf("PTB-WARNING: PsychOSCloseWindow: Failed to reenable DWM Aero Windows desktop compositor!\n");
				  }
			  }
		  }
		  
		  // Detach and free:
          FreeLibrary(dwmlibrary);
          dwmlibrary = 0;
          dwmSupported = FALSE;
      }
    }
  }

  // Done.
  return;
}

/*
    PsychOSGetVBLTimeAndCount()

    Returns absolute system time of last VBL and current total count of VBL interrupts since
    startup of gfx-system for the given screen. Returns a time of -1 and a count of 0 if this
    feature is unavailable on the given OS/Hardware configuration.
*/
double  PsychOSGetVBLTimeAndCount(PsychWindowRecordType *windowRecord, psych_uint64* vblCount)
{
	DWM_TIMING_INFO	dwmtiming;	
	psych_uint64 ust, msc, sbc;
	CGDirectDisplayID displayID;
	HRESULT rc = 0xdeadbeef;
	unsigned int screenid = windowRecord->screenNumber;

	// Retrieve displayID, aka HDC for this screenid:
	PsychGetCGDisplayIDFromScreenNumber(&displayID, screenid);

	// Windows Vista DWM available, supported and enabled?
	dwmtiming.cbSize = sizeof(dwmtiming);
	
	if (PsychOSIsDWMEnabled(0) && (NULL != PsychDwmGetCompositionTimingInfo) &&
		((rc = PsychDwmGetCompositionTimingInfo(NULL, &dwmtiming)) == 0)) {
		// Yes. Supported, enabled, and we got timing info from it. Extract:
		
		// VBLCount of last VBL:
		*vblCount = (psych_uint64) dwmtiming.cDXRefresh;
		
		// VBLTime of last VBL in QPC, ie., as query performance counter 64-bit psych_uint64 value:
		ust = (psych_uint64) dwmtiming.qpcVBlank;

		if (PsychPrefStateGet_Verbosity() > 15) {
			 printf("PTB-DEBUG: VBLCount = %i : VBLTime = %f secs. ClockHz = %f. rc=%x\n", (int) *vblCount, PsychMapPrecisionTimerTicksToSeconds(ust), PsychGetKernelTimebaseFrequencyHz(), rc);
		}

		// Convert ust into regular GetSecs timestamp:
		return(PsychMapPrecisionTimerTicksToSeconds(ust));
	}
	else {
		 if (PsychOSIsDWMEnabled(0) && PsychPrefStateGet_Verbosity()>6) {
			 printf("PTB-DEBUG: Call to PsychDwmGetCompositionTimingInfo(%i) failed with rc = %x, GetLastError() = %i\n", dwmtiming.cbSize, rc, GetLastError());
		}
	}

	// DWM unsupported, unavailable, disabled or failed to query if we reach this point.
	// Let's try if we have more luck with OpenML support...
	
	// Ok, this will return VBL count and last VBL time via the OML GetSyncValuesOML call
	// if that extension is supported on this setup. As of mid 2009 i'm not aware of any
	// affordable graphics card that would support this extension, but who knows??
	if ((NULL != wglGetSyncValuesOML) && (wglGetSyncValuesOML((HDC) displayID, (INT64*) &ust, (INT64*) &msc, (INT64*) &sbc))) {
		*vblCount = msc;
		
		// Convert ust into regular GetSecs timestamp:
		return(PsychMapPrecisionTimerTicksToSeconds(ust));
	}
	else {
		// Unsupported on Windows so far :(
		*vblCount = 0;
		return(-1);
	}
}

/*
	PsychOSGetPresentationTimingInfo()

	Retrieve low-level counts and timestamps related to stimulus onset
	for the last presented frame via PsychOSFlipWindowBuffers().
	
	Returns true on success, false on error.

*/
psych_bool PsychOSGetPresentationTimingInfo(PsychWindowRecordType *windowRecord, psych_bool postSwap, unsigned int flags, psych_uint64* onsetVBLCount, double* onsetVBLTime, psych_uint64* frameId, double* compositionRate, int fullStateStructReturnArgPos)
{
	DWM_TIMING_INFO	dwmtiming;	
	psych_uint64 qpcFrameDisplayed;
	static double qpcfreq = -1;
	HRESULT rc1 = 0;

	// First time call?
	if (qpcfreq == -1) {
		// Query and assign QPC frequency: This just to make sure that our QPC timebase is working correctly:
		qpcfreq = PsychGetKernelTimebaseFrequencyHz();

		if (PsychPrefStateGet_Verbosity() > 15) {
			 printf("PTB-DEBUG: Initial ClockHz = %f : Health = %i\n", qpcfreq, PsychGetTimeBaseHealthiness());
		}
		
		// Reliable?
		if ((PsychGetTimeBaseHealthiness() & 1) || (qpcfreq < 10000)) {
			// QPC doesn't work correctly. We disable ourselves:
			qpcfreq = 0;
		}

		if (PsychPrefStateGet_Verbosity() > 15) {
			 printf("PTB-DEBUG: Final ClockHz = %f\n", qpcfreq);
		}
	}
	
	// Fail if qpcfreq == 0 ie. timing trouble:
	if (qpcfreq == 0) return(FALSE);

	// Windows Vista DWM available, supported and enabled?
	dwmtiming.cbSize = sizeof(dwmtiming);
	if ( PsychOSIsDWMEnabled(0) && ((rc1 = PsychDwmGetCompositionTimingInfo(NULL, &dwmtiming)) == 0) ) {
		// Yes. Supported, enabled, and we got valid timing info from it. Extract:

		// Only qpcRefreshPeriod requested?
		if (flags == 1) {
			// Yes. Only return this in compositionRate:
			*compositionRate = PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcRefreshPeriod);
			return(TRUE);
		}

		// qpcFrameDisplayed seems to correspond almost perfectly to our concept of stimulus onset time, so use this:
		qpcFrameDisplayed = (psych_uint64) dwmtiming.qpcFrameDisplayed;
		
		// Convert to GetSecs() time:
		*onsetVBLTime = PsychMapPrecisionTimerTicksToSeconds(qpcFrameDisplayed);
		
		// VBL count of stimulus onset:
		*onsetVBLCount = (psych_uint64) dwmtiming.cRefreshFrameDisplayed;
		
		// Assumed frame id: Basically a flip count - Very first bufferswap is id 0, and so on -- A unique
		// serial number for each Screen('Flip') request: This is the id of the most recently completed flip:
		*frameId = (psych_uint64) dwmtiming.cFramesDisplayed;

		// Current composition rate of the DWM: Ideally at least our video refresh rate.
		*compositionRate = (double) dwmtiming.rateCompose.uiNumerator / (double) dwmtiming.rateCompose.uiDenominator;

		if (PsychPrefStateGet_Verbosity() > 15) {
			printf("PTB-DEBUG: === PsychDwmGetCompositionTimingInfo returned data follows: ===\n\n");
			printf("qpcFrameDisplayed       : %15.6f \n", PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcFrameDisplayed));
			printf("qpcRefreshPeriod        : %15.6f \n", PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcRefreshPeriod));
			printf("qpcFrameComplete        : %15.6f \n", PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcFrameComplete));
			printf("cFrameComplete          : %i \n", (psych_uint64) dwmtiming.cFrameComplete);
			printf("cFramePending           : %i \n", (psych_uint64) dwmtiming.cFramePending);
			printf("cFramesDisplayed        : %i \n", (psych_uint64) dwmtiming.cFramesDisplayed);
			printf("qpcFramePending         : %15.6f \n", PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcFramePending));
			printf("cRefreshFrameDisplayed  : %i \n", (psych_uint64) dwmtiming.cRefreshFrameDisplayed);
			printf("cRefreshStarted         : %i \n", (psych_uint64) dwmtiming.cRefreshStarted);
			printf("cFramesLate             : %i \n", (psych_uint64) dwmtiming.cFramesLate);
			printf("cFramesDropped          : %i \n", (psych_uint64) dwmtiming.cFramesDropped);
			printf("cFramesMissed           : %i \n", (psych_uint64) dwmtiming.cFramesMissed);
			printf("cBuffersEmpty           : %i \n", (psych_uint64) dwmtiming.cBuffersEmpty);
			
			printf("PTB-DEBUG: === End of PsychDwmGetCompositionTimingInfo returned data.  ===\n\n");
		}
		
		// Return of complete info structure with all retrieved data requested by caller?
		if (fullStateStructReturnArgPos > 0) {
			// Yes. Create a structure with all info:
			const char *DWMGraphicsFieldNames[]={ "rateRefresh", "qpcRefreshPeriod", "rateCompose", "qpcVBlank", "cRefresh", "cDXRefresh",
				"qpcCompose", "cFrame", "cDXPresent", "cRefreshFrame", "cFrameSubmitted", "cDXPresentSubmitted",
				"cFrameConfirmed", "cDXPresentConfirmed", "cRefreshConfirmed", "cDXRefreshConfirmed", "cFramesLate",
				"cFramesOutstanding", "cFrameDisplayed", "qpcFrameDisplayed", "cRefreshFrameDisplayed", "cFrameComplete", 
				"qpcFrameComplete", "cFramePending", "qpcFramePending", "cFramesDisplayed", "cFramesComplete",
				"cFramesPending", "cFramesAvailable", "cFramesDropped", "cFramesMissed", "cRefreshNextDisplayed",
				"cRefreshNextPresented", "cRefreshesDisplayed", "cRefreshesPresented", "cRefreshStarted",
				"cPixelsReceived", "cPixelsDrawn", "cBuffersEmpty" };
			const int DWMGraphicsFieldCount = 39;
			PsychGenericScriptType	*s;
			
			// Alloc struct and return it as return argument at position 'fullStateStructReturnArgPos':
			PsychAllocOutStructArray(fullStateStructReturnArgPos, FALSE, 1, DWMGraphicsFieldCount, DWMGraphicsFieldNames, &s);
			PsychSetStructArrayDoubleElement("rateRefresh", 0, (double) dwmtiming.rateRefresh.uiNumerator / (double) dwmtiming.rateRefresh.uiDenominator, s);
			PsychSetStructArrayDoubleElement("qpcRefreshPeriod", 0, PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcRefreshPeriod), s);
			PsychSetStructArrayDoubleElement("rateCompose", 0, (double) dwmtiming.rateCompose.uiNumerator / (double) dwmtiming.rateCompose.uiDenominator, s);
			PsychSetStructArrayDoubleElement("qpcVBlank", 0, PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcVBlank), s);
			PsychSetStructArrayDoubleElement("cRefresh", 0, (double) (psych_int64) dwmtiming.cRefresh, s);
			PsychSetStructArrayDoubleElement("cDXRefresh", 0, (double) (psych_int64) dwmtiming.cDXRefresh, s);
			PsychSetStructArrayDoubleElement("qpcCompose", 0, PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcCompose), s);
			PsychSetStructArrayDoubleElement("cFrame", 0, (double) (psych_int64) dwmtiming.cFrame, s);
			PsychSetStructArrayDoubleElement("cDXPresent", 0, (double) (psych_int64) dwmtiming.cDXPresent, s);
			PsychSetStructArrayDoubleElement("cRefreshFrame", 0, (double) (psych_int64) dwmtiming.cRefreshFrame, s);
			PsychSetStructArrayDoubleElement("cFrameSubmitted", 0, (double) (psych_int64) dwmtiming.cFrameSubmitted, s);
			PsychSetStructArrayDoubleElement("cDXPresentSubmitted", 0, (double) (psych_int64) dwmtiming.cDXPresentSubmitted, s);
			PsychSetStructArrayDoubleElement("cFrameConfirmed", 0, (double) (psych_int64) dwmtiming.cFrameConfirmed, s);
			PsychSetStructArrayDoubleElement("cDXPresentConfirmed", 0, (double) (psych_int64) dwmtiming.cDXPresentConfirmed, s);
			PsychSetStructArrayDoubleElement("cRefreshConfirmed", 0, (double) (psych_int64) dwmtiming.cRefreshConfirmed, s);
			PsychSetStructArrayDoubleElement("cDXRefreshConfirmed", 0, (double) (psych_int64) dwmtiming.cDXRefreshConfirmed, s);
			PsychSetStructArrayDoubleElement("cFramesLate", 0, (double) (psych_int64) dwmtiming.cFramesLate, s);
			PsychSetStructArrayDoubleElement("cFramesOutstanding", 0, (double) (psych_int64) dwmtiming.cFramesOutstanding, s);
			PsychSetStructArrayDoubleElement("cFrameDisplayed", 0, (double) (psych_int64) dwmtiming.cFrameDisplayed, s);
			PsychSetStructArrayDoubleElement("qpcFrameDisplayed", 0, PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcFrameDisplayed), s);
			PsychSetStructArrayDoubleElement("cRefreshFrameDisplayed", 0, (double) (psych_int64) dwmtiming.cRefreshFrameDisplayed, s);
			PsychSetStructArrayDoubleElement("cFrameComplete", 0, (double) (psych_int64) dwmtiming.cFrameComplete, s);
			PsychSetStructArrayDoubleElement("qpcFrameComplete", 0, PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcFrameComplete), s);
			PsychSetStructArrayDoubleElement("cFramePending", 0, (double) (psych_int64) dwmtiming.cFramePending, s);
			PsychSetStructArrayDoubleElement("qpcFramePending", 0, PsychMapPrecisionTimerTicksToSeconds(dwmtiming.qpcFramePending), s);
			PsychSetStructArrayDoubleElement("cFramesDisplayed", 0, (double) (psych_int64) dwmtiming.cFramesDisplayed, s);
			PsychSetStructArrayDoubleElement("cFramesComplete", 0, (double) (psych_int64) dwmtiming.cFramesComplete, s);
			PsychSetStructArrayDoubleElement("cFramesPending", 0, (double) (psych_int64) dwmtiming.cFramesPending, s);
			PsychSetStructArrayDoubleElement("cFramesAvailable", 0, (double) (psych_int64) dwmtiming.cFramesAvailable, s);
			PsychSetStructArrayDoubleElement("cFramesDropped", 0, (double) (psych_int64) dwmtiming.cFramesDropped, s);
			PsychSetStructArrayDoubleElement("cFramesMissed", 0, (double) (psych_int64) dwmtiming.cFramesMissed, s);
			PsychSetStructArrayDoubleElement("cRefreshNextDisplayed", 0, (double) (psych_int64) dwmtiming.cRefreshNextDisplayed, s);
			PsychSetStructArrayDoubleElement("cRefreshNextPresented", 0, (double) (psych_int64) dwmtiming.cRefreshNextPresented, s);
			PsychSetStructArrayDoubleElement("cRefreshesDisplayed", 0, (double) (psych_int64) dwmtiming.cRefreshesDisplayed, s);
			PsychSetStructArrayDoubleElement("cRefreshesPresented", 0, (double) (psych_int64) dwmtiming.cRefreshesPresented, s);
			PsychSetStructArrayDoubleElement("cRefreshStarted", 0, (double) (psych_int64) dwmtiming.cRefreshStarted, s);
			PsychSetStructArrayDoubleElement("cPixelsReceived", 0, (double) (psych_int64) dwmtiming.cPixelsReceived, s);
			PsychSetStructArrayDoubleElement("cPixelsDrawn", 0, (double) (psych_int64) dwmtiming.cPixelsDrawn, s);
			PsychSetStructArrayDoubleElement("cBuffersEmpty", 0, (double) (psych_int64) dwmtiming.cBuffersEmpty, s);
		}
		
		// Return success:
		return(TRUE);
	}
	else {
		if (PsychOSIsDWMEnabled(0) && PsychPrefStateGet_Verbosity() > 6) {
			printf("PTB-DEBUG: Call to PsychDwmGetCompositionTimingInfo() failed with rc1 = %x, GetLastError() = %i\n", rc1, GetLastError());			
		}
	}

	// DWM Unsupported, DWM disabled or query failed:
	return(FALSE);
}

/* PsychOSIsDWMEnabled()
 *
 * Return current Desktop Window Manager (DWM) status. Zero for disabled, Non-Zero for enabled.
 */
int	PsychOSIsDWMEnabled(int screenNumber)
{
    DWM_TIMING_INFO	dwmtiming;
    BOOL compositorEnabled;
    psych_bool IsDWMEnabled;

    // screenNumber unused on MS-Windows:
    (void) screenNumber;
    
    // Need to init dwmtiming for our dummy-call to get composition timing info:
    dwmtiming.cbSize = sizeof(dwmtiming);

    // If DWM is nominally supported, query and assign its effective enable state. If the OS reports it is enabled,
    // check if it is really active, ie., compositing instead of just on standby. The PsychDwmIsCompositionEnabled()
    // reports nominal enable state, but an enabled compositor still (should) get out of our way and therefore be
    // effectively disabled if a fullscreen window is displayed. A query to PsychDwmGetCompositionTimingInfo() will
    // succeed if the DWM is really active, but fail if the DWM is on standby.
    IsDWMEnabled = (dwmSupported && (0 == PsychDwmIsCompositionEnabled(&compositorEnabled)) && compositorEnabled && (0 == PsychDwmGetCompositionTimingInfo(NULL, &dwmtiming)));
    return(IsDWMEnabled);
}

/* PsychOSSetPresentParameters()
 *
 * Set presentation timing parameters for DWM that affect PsychOSFlipWindowBuffers()
 * presentation timing.
 *
 * 'targetVBL'	VSYNC count at which the framebuffer (backbuffer) content of the next PsychOSFlipWindowBuffers()
 *				call after this call should be displayed - a target "deadline" for stimulus onset, expressed in
 *				video refreshes. As far as i understand, this value affects the first bufferswap after the call.
 *				After the first bufferswap, successive bufferswaps will be timed wrt. each other and the DWM
 *				cycle, according to the following present parameters.
 *
 * 'rateDuration'	If set to a positive value (>0), defines the target refresh rate for swapbuffers ops, ie.
 *					a value of 100 would want the DWM to schedule swaps at a rate of 100 Hz. If set to a value
 *					<=0, the abs(rateDuration) encodes the delta between successive bufferswaps in refresh cycles,
 *					e.g., -10 == Execute next bufferswap 10 refresh cycles after the previous one.
 *					The most useful settings for us would be either rateDuration == nominal video refresh rate of
 *					display (Hz), or -1 for "swap each refresh" or maybe 0 for "swap immediately" (if zero is a valid setting?).
 *
 * 'queueLength"	Number of available buffers in flipqueue, ie., how many frames can be queued up. Valid range is
 *					between 2 and 8 for Vista and Windows-7.
 *
 * Returns TRUE on success, FALSE on error, or if the DWM is disabled or unsupported on pre-Vista systems.
 *
 */
psych_bool PsychOSSetPresentParameters(PsychWindowRecordType *windowRecord, psych_uint64 targetVBL, unsigned int queueLength, double rateDuration)
{
	int rc;
	DWM_PRESENT_PARAMETERS dwmPresentParams;
	UNSIGNED_RATIO rateSource;
	
	// Map positive rateDuration to some refresh rate:
	if (rateDuration > 0) {
		rateSource.uiNumerator = (int) rateDuration;
	}
	else {
		rateSource.uiNumerator = 0;
	}
	rateSource.uiDenominator = 1;
	
	dwmPresentParams.cbSize = sizeof(dwmPresentParams);
	dwmPresentParams.fQueue = (targetVBL > 0) ? TRUE : FALSE;
	dwmPresentParams.cRefreshStart = targetVBL;
	dwmPresentParams.cBuffer = queueLength;
	dwmPresentParams.fUseSourceRate = (rateDuration > 0) ? TRUE : FALSE;
	dwmPresentParams.rateSource = rateSource;
	dwmPresentParams.cRefreshesPerFrame = (UINT) ((rateDuration <= 0) ? -rateDuration : 0);
	dwmPresentParams.eSampling = DWM_SOURCE_FRAME_SAMPLING_POINT;
	
	// Call function if DWM is supported and enabled:
	if (PsychOSIsDWMEnabled(0) && ((rc = PsychDwmSetPresentParameters(windowRecord->targetSpecific.windowHandle, &dwmPresentParams)) == 0)) return(TRUE);
	
	// DWM unsupported, disabled, or call failed:
	return(FALSE);
}

/* PsychOSGetSwapCompletionTimestamp()
 *
 * Retrieve a very precise timestamp of doublebuffer swap completion by means
 * of OS specific facilities. This function is optional. If the underlying
 * OS/driver/GPU combo doesn't support a high-precision, high-reliability method
 * to query such timestamps, the function should return -1 as a signal that it
 * is unsupported or (temporarily) unavailable. Higher level timestamping code
 * should use/prefer timestamps returned by this function over other timestamps
 * provided by other mechanisms if possible. Calling code must be prepared to
 * use alternate timestamping methods if this method fails or returns a -1
 * unsupported error. Calling code must expect this function to block until
 * swap completion.
 *
 * Input argument targetSBC: Swapbuffers count for which to wait for. A value
 * of zero means to block until all pending bufferswaps for windowRecord have
 * completed, then return the timestamp of the most recently completed swap.
 *
 * A value of zero is recommended.
 *
 * Returns: Highly precise and reliable swap completion timestamp in seconds of
 * system time in variable referenced by tSwap, and msc value of completed swap,
 * or a negative value on error (-1 == unsupported, -2 == Query failed).
 *
 */
psych_int64 PsychOSGetSwapCompletionTimestamp(PsychWindowRecordType *windowRecord, psych_int64 targetSBC, double* tSwap)
{
    DWM_TIMING_INFO	dwmtiming;
    HRESULT rc;
    psych_uint64 msc, ust;
    double tSwapMapped;

    // We only support (or need) Windows-OS specific swap timestamping if the DWM
    // is active for our onscreen windows display. Fallback to default timestamping
    // if DWM is inactive:
    if (!PsychOSIsDWMEnabled(0)) return(-1);

    // DWM active on at least one display. Check if this is a single-display setup,
    // in which case we can try to wait for and timestamp swap completion. On a
    // multi-display setup we can't do that 
    if (PsychGetNumDisplays() > 1) return(-1);

    // Flush the DWM queue for our applications outstanding DirectX rendering requests:
    PsychDwmFlush();

    // DWM active on a single-display setup. This means that our onscreen window(s)
    // will be affected by compositing. Let's be optimistic and assume that the
    // timing info provided by the DWM is relevant to our onscreen window and not
    // polluted / confounded by display updates of potential other unrelated windows
    // on the composited desktop. Worst case, the confounding Present() requests
    // should cause us to return a too late timestamp, ie., false detection of deadline
    // miss, which is the more acceptable error case.
    while (TRUE) {
        // Yield for one msec:
        PsychYieldIntervalSeconds(0.001);

        // Get current timing state: Trigger standard fallback on failure.
        dwmtiming.cbSize = sizeof(dwmtiming);
        if ((rc = PsychDwmGetCompositionTimingInfo(NULL, &dwmtiming)) != 0) {
        if (PsychPrefStateGet_Verbosity() > 6) printf("PTB-DEBUG: PsychOSGetSwapCompletionTimestamp failed: Call to PsychDwmGetCompositionTimingInfo() failed with rc = %x, GetLastError() = %i\n", rc, GetLastError());
            return(-1);
        }

        // Default targetSBC is Present request number of most recently submitted
        // DirectX Present() request. This assumes the latest Present() request was
        // triggered by the DWM in response to the composition pass which was triggered
        // by our OpenGL SwapBuffers() call for our onscreen window:
        if (targetSBC == 0) targetSBC = dwmtiming.cDXPresentSubmitted;

        // Do have all outstanding DirectX Present() requests up to our targetSBC
        // have completed, ie., their corresponding buffers were flipped onscreen?
        // This would signal effective swap completion for the target composition pass
        // which involves our onscreen windows bufferswap. If not, we repeat the polling
        // loop for another wait-iteration:
        if (targetSBC == dwmtiming.cDXPresentConfirmed) break;
    }

    // Present confirmed completed. To the best of our limited knowledge, the
    // DirectX refresh count of present confirmation should be the vblank count
    // of the vblank in which the swap completed, ie., the msc. Calculate the
    // difference between current vblank count and count of swap completion,
    // translate it into a time delta from the most recent vblank and subtract that
    // delte from the QPC timestamp of the most recent vblank. This should yield
    // our best estimate of swap completion time:
    msc = dwmtiming.cDXRefreshConfirmed;
    ust = dwmtiming.qpcVBlank - ((dwmtiming.cDXRefresh - dwmtiming.cDXRefreshConfirmed) * dwmtiming.qpcRefreshPeriod);

    // Translate to GetSecs time:
    tSwapMapped = PsychMapPrecisionTimerTicksToSeconds(ust);

    // A last consistency check. Known to trigger at least sometimes when running Windows-7 in a VM:
    if (msc < 0 || tSwapMapped < 0) return(-1);

    // Assign swap completion timestamp:
    if (tSwap) *tSwap = tSwapMapped;

    // Return msc of swap completion:
    return((psych_int64) msc);
}

/*
    PsychOSScheduleFlipWindowBuffers()
    
    Schedules a double buffer swap operation for given window at a given
	specific target time or target refresh count in a specified way.
	
	This uses OS specific API's and algorithms to schedule the asynchronous
	swap. This function is optional, target platforms are free to not implement
	it but simply return a "not supported" status code.
	
	Arguments:
	
	windowRecord - The window to be swapped.
	tWhen        - Requested target system time for swap. Swap shall happen at first
				   VSync >= tWhen.
	targetMSC	 - If non-zero, specifies target msc count for swap. Overrides tWhen.
	divisor, remainder - If set to non-zero, msc at swap must satisfy (msc % divisor) == remainder.
	specialFlags - Additional options. Unused so far.
	
	Return value:
	 
	Value greater than or equal to zero on success: The target msc for which swap is scheduled.
	Negative value: Error. Function failed. -1 == Function unsupported on current system configuration.
	-2 ... -x == Error condition.
	
*/
psych_int64 PsychOSScheduleFlipWindowBuffers(PsychWindowRecordType *windowRecord, double tWhen, psych_int64 targetMSC, psych_int64 divisor, psych_int64 remainder, unsigned int specialFlags)
{
	// On Windows this function is unsupported:
	return(-1);
}

/*
    PsychOSFlipWindowBuffers()
    
    Performs OS specific double buffer swap call.
*/
void PsychOSFlipWindowBuffers(PsychWindowRecordType *windowRecord)
{
	// Execute OS neutral bufferswap code first:
	PsychExecuteBufferSwapPrefix(windowRecord);

	// Trigger the "Front <-> Back buffer swap (flip) (on next vertical retrace)":
	SwapBuffers(windowRecord->targetSpecific.deviceContext);
}

/* Enable/disable syncing of buffer-swaps to vertical retrace. */
void PsychOSSetVBLSyncLevel(PsychWindowRecordType *windowRecord, int swapInterval)
{
  static unsigned int failcount = 0;
  
  // Enable rendering context of window:
  PsychSetGLContext(windowRecord);
  
  // Store new setting also in internal helper variable, e.g., to allow workarounds to work:
  windowRecord->vSynced = (swapInterval > 0) ? TRUE : FALSE;

  // Try to set requested swapInterval if swap-control extension is supported on
  // this windows machine. Otherwise this will be a no-op...
  if (wglSwapIntervalEXT) {
	if(!wglSwapIntervalEXT(swapInterval)) {
		failcount++;
		if (failcount <= 10) printf("PTB-ERROR: Setting wglSwapInterval(%i) failed! Expect severe display timing and display tearing problems!!! See 'help SyncTrouble' for more info.\n", swapInterval);
		if (failcount == 10) printf("PTB-ERROR: This error message won't repeat on subsequent failure...\n");
	}
	
	// Double check new setting:
	if ((NULL == wglGetSwapIntervalEXT) || (wglGetSwapIntervalEXT() != swapInterval)) {
		failcount++;
		if (failcount <= 10) {
			if (NULL == wglGetSwapIntervalEXT) {
				printf("PTB-ERROR: wglGetSwapIntervalEXT() unsupported. This is a bug in your graphics driver or system setup!! See 'help SyncTrouble' for more info.\n");
			}
			else {
				printf("PTB-ERROR: Mismatch between requested display swap interval %i and actual swap interval %i! Synchronization of Screen('Flip') to display refresh will likely malfunction!\n", swapInterval, wglGetSwapIntervalEXT());
				printf("PTB-ERROR: This is either a bug in your graphics driver or system setup, or some misconfiguration in the display setting control panel of your system. See 'help SyncTrouble' for more info.\n");
			}
		}
		if (failcount == 10) printf("PTB-ERROR: These error messages won't repeat on subsequent failure...\n");
	}
  }
  
  return;
}

/*
    PsychOSSetGLContext()
    
    Set the window to which GL drawing commands are sent.  
*/
void PsychOSSetGLContext(PsychWindowRecordType *windowRecord)
{
	if (wglGetCurrentContext() != windowRecord->targetSpecific.contextObject) {
		if (wglGetCurrentContext() != NULL) {
			// We need to glFlush the context before switching, otherwise race-conditions may occur:
			glFlush();
			
			// Need to unbind any FBO's in old context before switch, otherwise bad things can happen...
			if (glBindFramebufferEXT) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		}

		// Switch to new context:
		wglMakeCurrent(windowRecord->targetSpecific.deviceContext, windowRecord->targetSpecific.contextObject);		
	}
}

/*
    PsychOSUnsetGLContext()
    
    Clear the drawing context.  
*/
void PsychOSUnsetGLContext(PsychWindowRecordType* windowRecord)
{
    if (wglGetCurrentContext() != NULL) {
        // We need to glFlush the context before switching, otherwise race-conditions may occur:
        glFlush();
        
        // Need to unbind any FBO's in old context before switch, otherwise bad things can happen...
        if (glBindFramebufferEXT) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    }
    
    wglMakeCurrent(windowRecord->targetSpecific.deviceContext, NULL);
}

/* Same as PsychOSSetGLContext() but for selecting userspace rendering context,
 * optionally copying state from PTBs context.
 */
void PsychOSSetUserGLContext(PsychWindowRecordType *windowRecord, psych_bool copyfromPTBContext)
{
	// Child protection:
	if (windowRecord->targetSpecific.glusercontextObject == NULL) PsychErrorExitMsg(PsychError_user, "GL Userspace context unavailable! Call InitializeMatlabOpenGL *before* Screen('OpenWindow')!");
	
	if (copyfromPTBContext) {
  		wglMakeCurrent(windowRecord->targetSpecific.deviceContext, NULL);		
		wglCopyContext(windowRecord->targetSpecific.contextObject, windowRecord->targetSpecific.glusercontextObject, GL_ALL_ATTRIB_BITS);
	}
	
    // Setup new context if it isn't already setup. -> Avoid redundant context switch.
  	 if (wglGetCurrentContext() != windowRecord->targetSpecific.glusercontextObject) {
		 wglMakeCurrent(windowRecord->targetSpecific.deviceContext, windowRecord->targetSpecific.glusercontextObject);
    }
}

/* PsychOSSetupFrameLock - Check if framelock / swaplock support is available on
 * the given graphics system implementation and try to enable it for the given
 * pair of onscreen windows.
 *
 * If possible, will try to add slaveWindow to the swap group and/or swap barrier
 * of which masterWindow is already a member, putting slaveWindow into a swap-lock
 * with the masterWindow. If masterWindow isn't yet part of a swap group, create a
 * new swap group and attach masterWindow to it, before joining slaveWindow into the
 * new group. If masterWindow is part of a swap group and slaveWindow is NULL, then
 * remove masterWindow from the swap group.
 *
 * The swap lock mechanism used is operating system and GPU dependent. Many systems
 * will not support framelock/swaplock at all.
 *
 * Returns TRUE on success, FALSE on failure.
 */
psych_bool PsychOSSetupFrameLock(PsychWindowRecordType *masterWindow, PsychWindowRecordType *slaveWindow)
{
	GLuint maxGroups, maxBarriers, targetGroup;
	psych_bool rc = FALSE;

	// MS-Windows: Only NV_swap_group support. Try it.
	
	// NVidia swap group extension supported?
	if((wglewIsSupported("WGL_NV_swap_group") || glewIsSupported("WGL_NV_swap_group")) && (NULL != wglQueryMaxSwapGroupsNV)) {
		// Yes. Check if given GPU really supports it:
		if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: NV_swap_group supported. Querying available groups...\n");

		if (wglQueryMaxSwapGroupsNV(masterWindow->targetSpecific.deviceContext, &maxGroups, &maxBarriers) && (maxGroups > 0)) {
			// Yes. What to do?
			if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: NV_swap_group supported. Implementation supports up to %i swap groups. Trying to join or unjoin group.\n", maxGroups);

			if (NULL == slaveWindow) {
				// Asked to remove master from swap group:
				wglJoinSwapGroupNV(masterWindow->targetSpecific.deviceContext, 0);
				masterWindow->swapGroup = 0;
				return(TRUE);
			}
			else {
				// Non-NULL slaveWindow: Shall attach to swap group.
				// Master already part of a swap group?
				if (0 == masterWindow->swapGroup) {
					// Nope. Try to attach it to first available one:
					targetGroup = (GLuint) PsychFindFreeSwapGroupId(maxGroups);
					
					if ((targetGroup == 0) || !wglJoinSwapGroupNV(masterWindow->targetSpecific.deviceContext, targetGroup)) {
						// Failed!
						if (PsychPrefStateGet_Verbosity() > 1) {
							printf("PTB-WARNING: Tried to enable framelock support for master-slave window pair, but masterWindow failed to join swapgroup %i! Skipped.\n", targetGroup);
						}
						
						return(FALSE);
					}
					
					// Sucess for master!
					masterWindow->swapGroup = targetGroup;
				}
				
				// Now try to join the masters swapgroup with the slave:
				if (!wglJoinSwapGroupNV(slaveWindow->targetSpecific.deviceContext, masterWindow->swapGroup)) {
					// Failed!
					if (PsychPrefStateGet_Verbosity() > 1) {
						printf("PTB-WARNING: Tried to enable framelock support for master-slave window pair, but slaveWindow failed to join swapgroup %i of master! Skipped.\n", masterWindow->swapGroup);
					}
					
					return(FALSE);
				}
				
				// Success! Now both windows are in a common swapgroup and framelock should work!
				slaveWindow->swapGroup = masterWindow->swapGroup;
				
				if (PsychPrefStateGet_Verbosity() > 1) {
					printf("PTB-INFO: Framelock support for master-slave window pair via NV_swap_group extension enabled! Joined swap group %i.\n", masterWindow->swapGroup);
				}
				
				return(TRUE);
			}
		}
	}

	if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: NV_swap_group unsupported or join operation failed.\n");

	return(rc);
}
