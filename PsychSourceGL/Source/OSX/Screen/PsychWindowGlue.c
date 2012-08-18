/*
	PsychToolbox3/Source/OSX/Screen/PsychWindowGlue.c
	
	PLATFORMS:	
	
		This is the OS X Core Graphics version only.  
				
	AUTHORS:
	
		Allen Ingling		awi		Allen.Ingling@nyu.edu
        Mario Kleiner       mk      mario.kleiner at tuebingen.mpg.de

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
                 01/02/05               mk              Modified to only contain the OS-X specific code. All OS independent code has been moved to
                                                        Common/Screen/PsychWindowSupport.c
 
	DESCRIPTION:
	
		Functions in this file comprise an abstraction layer for probing and controlling window state, except for window content.  
		
		Each C function which implements a particular Screen subcommand should be platform neutral.  For example, the source to SCREENPixelSizes() 
		should be platform-neutral, despite that the calls in OS X and Windows to detect available pixel sizes are different.  The platform 
		specificity is abstracted out in C files which end it "Glue", for example PsychScreenGlue, PsychWindowGlue, PsychWindowTextClue.

	NOTES:
	
	TO DO: 
	
		¥ The "glue" files should should be suffixed with a platform name.  The original (bad) plan was to distingish platform-specific files with the same 
		name by their placement in a directory tree.
 
*/

#include "Screen.h"

// Need AGL headers for windowed and multi-screen rendering:
#include <AGL/agl.h>

// Includes for low-level access to IOKit Framebuffer device:
#include <CoreFoundation/CoreFoundation.h>
#include <CoreVideo/CoreVideo.h>
#include <ApplicationServices/ApplicationServices.h>
#include <IOKit/graphics/IOGraphicsLib.h>
#include <IOKit/graphics/IOFramebufferShared.h>
#include <libkern/OSAtomic.h>

// Disable warnings about deprecated API calls on OSX 10.7
// of which we are aware and that we can't remove as long as
// we need to stay compatible to 10.4 - 10.6
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

// Surrogate routines for missing implementations on MacOS/X 64-Bit:
#ifdef __LP64__

#include "PsychCocoaGlue.h"

OSErr DMGetGDeviceByDisplayID(
                        int   displayID,
                        GDHandle *      displayDevice,
                        Boolean         failToMain)
{
    // Return failure:
    return 0x1;
}

CGrafPtr GetWindowPort(WindowRef window) { return 0; }

#endif

static struct {
    io_connect_t        connect;
    StdFBShmem_t*       shmem;
    mach_vm_size_t      shmemSize;    
} fbsharedmem[kPsychMaxPossibleDisplays];   

static struct {
    psych_mutex         mutex;
    double              vblTimestamp;
    psych_uint64        vblCount;
} cvDisplayLinkData[kPsychMaxPossibleDisplays];

static CVDisplayLinkRef cvDisplayLink[kPsychMaxPossibleDisplays] = { NULL };
static int screenRefCount[kPsychMaxPossibleDisplays] = { 0 };

static SInt32 osMajor, osMinor;
static psych_bool useCoreVideoTimestamping;

// Display link callback: Needed so we can actually start the display link:
// Gets apparently called from a separate high-priority thread, close to vblank
// time. "inNow" is the timestamp of last vblank.
static CVReturn PsychCVDisplayLinkOutputCallback(CVDisplayLinkRef displayLink, const CVTimeStamp *inNow, const CVTimeStamp *inOutputTime,
                                                 CVOptionFlags flagsIn, CVOptionFlags *flagsOut, void *displayLinkContext)
{
    double tVBlank;
    CVTimeStamp tVbl;
    double tHost;
    
    // Retrieve screenId of associated display screen:
    int screenId = (int) (long int) displayLinkContext;
    
    // Translate CoreVideo inNow timestamp with time of last vbl from gpu time
    // to host system time, aka our GetSecs() timebase:
    memset(&tVbl, 0, sizeof(tVbl));
    tVbl.version = 0;
    tVbl.flags = kCVTimeStampHostTimeValid;
    CVDisplayLinkTranslateTime(displayLink, inNow, &tVbl);
    tVBlank = (double) tVbl.hostTime / (double) 1000000000;

    // Store timestamp in our shared data structure, also increment virtual vblank counter:
    PsychLockMutex(&(cvDisplayLinkData[screenId].mutex));
    cvDisplayLinkData[screenId].vblCount++;
    cvDisplayLinkData[screenId].vblTimestamp = tVBlank;
    PsychUnlockMutex(&(cvDisplayLinkData[screenId].mutex));
    
    // Low-level timestamp debugging requested?
    if (PsychPrefStateGet_Verbosity() > 20) {
        // Compare CV timestamps against host time for correctness check. We wait 4 msecs,
        // then take tHost and hopefully tHost will be at least 4 msecs later than the
        // computed vblank timestamp tVBlank:
        PsychWaitIntervalSeconds(0.004);
        PsychGetAdjustedPrecisionTimerSeconds(&tHost);
        
        // Caution: Don't run from Matlab GUI! This printf will crash Matlab otherwise.
        printf("CVCallback: %i : tHost = %lf secs, tVBlank = %lf secs. tHost - tVBlank = %lf secs.\n", screenId, tHost, tVBlank, tHost - tVBlank);
    }

    return(kCVReturnSuccess);
}

/** PsychRealtimePriority: Temporarily boost priority to THREAD_TIME_CONSTRAINT_POLICY.
    PsychRealtimePriority(true) enables realtime-scheduling (like Priority(9) would do in Matlab).
    PsychRealtimePriority(false) restores scheduling to the state before last invocation of PsychRealtimePriority(true),
    it undos whatever the previous switch did.

    We switch to RT scheduling during PsychGetMonitorRefreshInterval() and a few other timing tests in
    PsychOpenWindow() to reduce measurement jitter caused by possible interference of other tasks.
*/
psych_bool PsychRealtimePriority(psych_bool enable_realtime)
{
    psych_bool				isError;
    thread_policy_flavor_t	flavorConstant;
    int						kernError;
    task_t					threadID;
    thread_policy_t			threadPolicy;
    static thread_policy_t	old_threadPolicy;
    mach_msg_type_number_t	policyCount, policyCountFilled;
    static mach_msg_type_number_t	old_policyCountFilled;
    boolean_t				isDefault;
    
    static psych_bool old_enable_realtime = FALSE;
    static psych_bool oldModeWasStandard = FALSE;
    
    if (old_enable_realtime == enable_realtime) {
        // No transition with respect to previous state -> Nothing to do.
        return(TRUE);
    }
    
    // Transition requested:
    old_enable_realtime = enable_realtime;
    
    // Determine our threadID:
    threadID = mach_thread_self();
    
    if (enable_realtime) {
        // Transition to realtime requested:

        // Get current scheduling policy and its settings and back it up for later restore:
        old_threadPolicy = (thread_policy_t) malloc(sizeof(thread_time_constraint_policy_data_t));
        policyCount = THREAD_TIME_CONSTRAINT_POLICY_COUNT;
        old_policyCountFilled = policyCount;
        isDefault = FALSE;
        // We check if STANDARD_POLICY is active and query its settings, if so...
        kernError = thread_policy_get(threadID, THREAD_STANDARD_POLICY, old_threadPolicy, &old_policyCountFilled, &isDefault);
        if (kernError) {
            // Failed!
            old_enable_realtime = FALSE;
            free(old_threadPolicy);
			printf("PsychRealtimePriority: ERROR! COULDN'T QUERY CURRENT SCHEDULING SETTINGS!!!\n");
            return(FALSE);
        }
        
        // oldModeWasStandard == TRUE --> We need to revert to STANDARD POLICY later...
        oldModeWasStandard = !isDefault;

        // printf("PRE-RT: CURRENTLY IN %s mode\n", oldModeWasStandard ? "STANDARD" : "REALTIME");

		if (!oldModeWasStandard) {
			// We are already RT scheduled. Backup settings for later switch-back:
			policyCount = THREAD_TIME_CONSTRAINT_POLICY_COUNT;
			old_policyCountFilled = policyCount;
			isDefault = FALSE;
			// We check if STANDARD_POLICY is active and query its settings, if so...
			kernError = thread_policy_get(threadID, THREAD_TIME_CONSTRAINT_POLICY, old_threadPolicy, &old_policyCountFilled, &isDefault);
			if (kernError) {
				// Failed!
				old_enable_realtime = FALSE;
				free(old_threadPolicy);
				printf("PsychRealtimePriority: ERROR! COULDN'T QUERY CURRENT RT SCHEDULING SETTINGS!!!\n");
				return(FALSE);
			}
		}

		// Switch to our ultra-high priority realtime mode: Guaranteed up to 3 msecs of uninterrupted
		// runtime as soon as we want to run: Perfect for swap completion timestamping in refresh rate
		// calibration - our only use-case:
		PsychSetThreadPriority(NULL, 10, 2);
    }
    else {
        // Transition from RT to Non-RT scheduling requested: We just reestablish the backed-up old
        // policy:
		kernError = thread_policy_set(threadID, (oldModeWasStandard) ? THREAD_STANDARD_POLICY : THREAD_TIME_CONSTRAINT_POLICY, old_threadPolicy, old_policyCountFilled);
		if (kernError) {
			// Failed!
			old_enable_realtime = TRUE;
			free((void*) old_threadPolicy);
			
			printf("PsychRealtimePriority: ERROR! COULDN'T SWITCH BACK TO NON-RT SCHEDULING!!!\n");
			fflush(NULL);
			return(FALSE);
		}
        
        // Successfully switchted to RT-Scheduling:
        free((void*) old_threadPolicy);
    }

    // Success.
    return(TRUE);
}


/*
    PsychOSOpenOnscreenWindow()
    
    Creates the CGL pixel format and the CGL context objects and then instantiates the context onto the screen.
    
    -The pixel format and the context are stored in the target specific field of the window recored.  Close
    should clean up by destroying both the pixel format and the context.
    
    -We mantain the context because it must be be made the current context by drawing functions to draw into 
    the specified window.
    
    -We maintain the pixel format object because there seems to be now way to retrieve that from the context.
    
    -To tell the caller to clean up PsychOpenOnscreenWindow returns FALSE if we fail to open the window. It 
    would be better to just issue an PsychErrorExit() and have that clean up everything allocated outside of
    PsychOpenOnscreenWindow().
    
    MK: The new option 'stereomode' allows selection of stereo display instead of mono display:
    0 (default) == Old behaviour -> Monoscopic rendering context.
    >0          == Stereo display, where the number defines the type of stereo algorithm to use.
    =1          == Use OpenGL built-in stereo by creating a context/window with left- and right backbuffer.
    =2          == Use compressed frame stereo: Put both views into one framebuffer, one in top half, other in lower half.

*/
psych_bool PsychOSOpenOnscreenWindow(PsychScreenSettingsType *screenSettings, PsychWindowRecordType *windowRecord, int numBuffers, int stereomode, int conserveVRAM)
{
    CGLRendererInfoObj				rendererInfo;
    CGOpenGLDisplayMask 			displayMask;
    CGLError						error;
	OSStatus						err;
    CGDirectDisplayID				cgDisplayID;
    CGLPixelFormatAttribute			attribs[40];
    GLint							numVirtualScreens;
    GLboolean						isDoubleBuffer, isFloatBuffer;
    GLint bpc;
	GLenum glerr;
	PsychRectType					screenrect;
    int attribcount=0;
    int i;
	int								windowLevel;
	psych_bool						useAGL, AGLForFullscreen;
	WindowRef						carbonWindow = NULL;
	int								aglbufferid;
	AGLPixelFormat					pf = NULL;
 	AGLContext						glcontext = NULL;
	GDHandle						gdhDisplay;
	GDHandle						*gdhDisplayPtr = NULL; 
	int								gdhcount = 0;
    psych_bool                      useCocoa = FALSE;
	    
    // Query OS/X version:
	Gestalt(gestaltSystemVersionMajor, &osMajor);
	Gestalt(gestaltSystemVersionMinor, &osMinor);

	// NULL-out Carbon/Cocoa window handle, so this is well-defined in case of error:
	windowRecord->targetSpecific.windowHandle = NULL;

	// Retrieve windowLevel, an indicator of where non-fullscreen AGL windows should
	// be located wrt. to other windows. 0 = Behind everything else, occluded by
	// everything else. 1 - 999 = At layer windowLevel -> Occludes stuff on layers "below" it.
	// 1000 - 1999 = At highest level, but partially translucent / alpha channel allows to make
	// regions transparent. 2000 or higher: Above everything, fully opaque, occludes everything.
	// 2000 is the default.
	windowLevel = PsychPrefStateGet_WindowShieldingLevel();

	// Use AGL for fullscreen windows?
	AGLForFullscreen = (PsychPrefStateGet_ConserveVRAM() & kPsychUseAGLForFullscreenWindows) ? TRUE : FALSE;

	// Window rect provided which has a different size than screen?
	useAGL = TRUE;
	
	// We do not use windowed mode if the provided window rectangle either
	// matches the target screens rectangle (and therefore its exact size)
	// or its screens global rectangle.
	PsychGetScreenRect(screenSettings->screenNumber, screenrect);
	if (PsychMatchRect(screenrect, windowRecord->rect)) useAGL=FALSE;
	PsychGetGlobalScreenRect(screenSettings->screenNumber, screenrect);
	if (PsychMatchRect(screenrect, windowRecord->rect)) useAGL=FALSE;

	// GUI windows always use AGL:
	if (windowRecord->specialflags & kPsychGUIWindow) useAGL = TRUE;

	// FALSE means: This is a fullscreen window:
	if (!useAGL) {
	  // Mark this window as fullscreen window:
	  windowRecord->specialflags |= kPsychIsFullscreenWindow;
	}

	// Override for use on f$%#$Fd OS/X 10.5.3 - 10.5.8 with NVidia GF 8800 GPU's:
	// If requested, always use AGL API and regular (non-fullscreen) windows via the
	// Quartz compositor. These are just borderless fullscreen windows.
	// This works with the NVidia GF 8800 driver bug on dual-display setups,
	// but as Quartz is in full control of buffer swaps, the stimulus onset timing
	// is horrible, animations are jerky and all our high precision timestamping is
	// completely broken and pointless: This mode is only useful for slowly updating
	// mostly static stimuli with no timing requirements.
	if (!useAGL && (PsychPrefStateGet_ConserveVRAM() & kPsychUseAGLCompositorForFullscreenWindows)) {
		// Force use of AGL in windowed mode, ie., via Quartz compositor:
		useAGL = TRUE;
		
		// Force a window rectangle that matches the global screen rectangle for that windows screen:
		PsychCopyRect(windowRecord->rect, screenrect);
		
		// Disable all timestamping:
		PsychPrefStateSet_VBLTimestampingMode(-1);
		
		// Warn user about what's going on:
		if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Using desktop compositor for composited onscreen window creation: High precision timestamping disabled,\nany kind of visual stimulus onset timing wil be very unreliable!!\n");		
	}

	// Do we need to use windowed mode with AGL?
	if (useAGL) {
		// Yes. Need to create Carbon window and attach OpenGL to it via AGL:		
		if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Using desktop compositor for onscreen window creation...\n");
		
		// Create onscreen Carbon window of requested position and size:
		Rect winRect;
		winRect.top = (short) windowRecord->rect[kPsychTop];
		winRect.left = (short) windowRecord->rect[kPsychLeft];
		winRect.bottom = (short) windowRecord->rect[kPsychBottom];
		winRect.right = (short) windowRecord->rect[kPsychRight];

		WindowClass wclass;
		if (windowLevel >= 1000) {
			wclass = kOverlayWindowClass;
		}
		else {
			wclass = kSimpleWindowClass;
		}

		if (windowRecord->specialflags & kPsychGUIWindow) wclass = kDocumentWindowClass;

		// Additional attribs to set:
		WindowAttributes addAttribs;
		if (windowRecord->specialflags & kPsychGUIWindow) {
			if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Onscreen window is configured as regular GUI window.\n");
			addAttribs = kWindowNoUpdatesAttribute | kWindowNoActivatesAttribute | kWindowStandardHandlerAttribute;
		}
		else {
			addAttribs = kWindowNoUpdatesAttribute | kWindowNoActivatesAttribute;
		}
        
		// For levels 1000 to 1499, where the window is a partially transparent
		// overlay window with global alpha 0.0 - 1.0, we disable reception of mouse
		// events. --> Can move and click to windows behind the window!
		// A range 1500 to 1999 would also allow transparency, but block mouse events:
		if (windowLevel >= 1000 && windowLevel < 1500) addAttribs += kWindowIgnoreClicksAttribute;
		
        #ifndef __LP64__
            // 32-Bit Carbon path:
            if (noErr != CreateNewWindow(wclass, addAttribs, &winRect, &carbonWindow)) {
                printf("\nPTB-ERROR[CreateNewWindow failed]: Failed to open Carbon onscreen window\n\n");
                return(FALSE);
            }
        #else
            // Request use of Cocoa specific setup path on 10.6 and later:
            if ((osMajor > 10) || ((osMajor == 10) && (osMinor > 5))) useCocoa = TRUE;

            // Usercode wants transparency?
            psych_bool enableTransparentGL = ((windowLevel >= 1000) && (windowLevel < 2000)) ? TRUE : FALSE;

            // 64-Bit Cocoa path:
            if (PsychCocoaCreateWindow(windowRecord, screenrect, &winRect, wclass, addAttribs, &carbonWindow, enableTransparentGL)) {
                printf("\nPTB-ERROR[CreateNewWindow failed]: Failed to open Carbon onscreen window\n\n");
                return(FALSE);
            }
        
            if (enableTransparentGL) {
                // Setup of global window alpha value for transparency. This is premultiplied to
                // the individual per-pixel alpha values if transparency is enabled by Cocoa code.
                //
                // Levels 1000 - 1499 and 1500 to 1999 map to a master opacity level of 0.0 - 1.0:
                SetWindowAlpha(carbonWindow, ((float) (windowLevel % 500)) / 499.0);
            }        
        #endif

		// Show it! Unless a windowLevel of -1 requests hiding the window:
		if (windowLevel != -1) ShowWindow(carbonWindow);

		// Level zero means: Place behind all other windows:
		if (windowLevel == 0) SendBehind(carbonWindow, NULL);

		// Levels 1 to 998 define window levels for the group of the window. A level
		// of 999 would leave this to the system:
		if (windowLevel > 0 && windowLevel < 999) SetWindowGroupLevel(GetWindowGroup(carbonWindow), (SInt32) windowLevel);

		// Store window handle in windowRecord:
		windowRecord->targetSpecific.windowHandle = carbonWindow;
		
		// Copy absolute screen location and area of window to 'globalrect',
		// so functions like Screen('GlobalRect') can still query the real
		// bounding gox of a window onscreen:
		PsychCopyRect(windowRecord->globalrect, windowRecord->rect);
		
		// Disable the fullscreen flag, in case it was set:
		AGLForFullscreen = FALSE;
	}
	else {
		// No. Standard CGL setup for fullscreen single display windows:
		if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Using %s for fullscreen onscreen window creation...\n", (AGLForFullscreen) ? "AGL" : "CGL");
		
		// Copy absolute screen location and area of window to 'globalrect',
		// so functions like Screen('GlobalRect') can still query the real
		// bounding gox of a window onscreen:
		PsychGetGlobalScreenRect(screenSettings->screenNumber, windowRecord->globalrect);		
	}

    // Map screen number to physical display handle cgDisplayID:
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenSettings->screenNumber);
    displayMask=CGDisplayIDToOpenGLDisplayMask(cgDisplayID);

	if (useCocoa || (!useAGL && !AGLForFullscreen)) {
		if (!useCocoa) attribs[attribcount++]=kCGLPFAFullScreen;
		attribs[attribcount++]=kCGLPFADisplayMask;
		attribs[attribcount++]=displayMask;

		// 10 bit per component framebuffer requested (10-10-10-2)?
		/*	Disabled: We use our special kernel helper driver and PsychEnableNative10BitFramebuffer(); instead on Radeon Gfx et al.
		if (windowRecord->depth == 30) {
			// Request a 10 bit per color component framebuffer with 2 bit alpha channel:
			printf("PTB-INFO: Trying to enable 10 bpc framebuffer...\n");
			attribs[attribcount++]=kCGLPFANoRecovery;
			attribs[attribcount++]=kCGLPFAColorSize;
			attribs[attribcount++]=16*3;
			attribs[attribcount++]=kCGLPFAAlphaSize;
			attribs[attribcount++]=16;
		}
		*/
		
		// 16 bit per component, 64 bit framebuffer requested (16-16-16-16)?
		if (windowRecord->depth == 64) {
			// Request a floating point framebuffer in 16-bit half-float format, i.e., RGBA = 16 bits per component.
			printf("PTB-INFO: Trying to enable 16 bpc float framebuffer...\n");
			attribs[attribcount++]=kCGLPFAColorFloat;
			attribs[attribcount++]=kCGLPFAColorSize;
			attribs[attribcount++]=16*3;
			attribs[attribcount++]=kCGLPFAAlphaSize;
			attribs[attribcount++]=16;
		}
		
		// 32 bit per component, 128 bit framebuffer requested (32-32-32-32)?
		if (windowRecord->depth == 128) {
			// Request a floating point framebuffer in 32-bit float format, i.e., RGBA = 32 bits per component.
			printf("PTB-INFO: Trying to enable 32 bpc float framebuffer...\n");
			attribs[attribcount++]=kCGLPFAColorFloat;
			attribs[attribcount++]=kCGLPFAColorSize;
			attribs[attribcount++]=32*3;
			attribs[attribcount++]=kCGLPFAAlphaSize;
			attribs[attribcount++]=32;
		}
	}
	else {
		// AGL specific pixelformat setup:

		if (AGLForFullscreen) attribs[attribcount++] = AGL_FULLSCREEN;

		// Possible to request use of the Apple floating point software renderer:
		if (conserveVRAM & kPsychUseSoftwareRenderer) {
			#ifndef kCGLRendererGenericFloatID
			#define kCGLRendererGenericFloatID    0x00020400
			#endif
			attribs[attribcount++]=AGL_RENDERER_ID;
			attribs[attribcount++]=kCGLRendererGenericFloatID;
		}

		// 32 bit per component, 128 bit framebuffer requested (32-32-32-32)?
		if (windowRecord->depth == 128) {
			// Request a floating point framebuffer in 32-bit float format, i.e., RGBA = 32 bits per component.
			printf("PTB-INFO: Trying to enable 32 bpc float framebuffer...\n");
			attribs[attribcount++]= AGL_RGBA;
			attribs[attribcount++]= AGL_COLOR_FLOAT;
		}
		else {
			attribs[attribcount++]= AGL_RGBA;
			attribs[attribcount++]=	AGL_RED_SIZE;
			attribs[attribcount++]= 8;
			attribs[attribcount++]=	AGL_GREEN_SIZE;
			attribs[attribcount++]= 8;
			attribs[attribcount++]=	AGL_BLUE_SIZE;
			attribs[attribcount++]= 8;
			attribs[attribcount++]=	AGL_ALPHA_SIZE;
			attribs[attribcount++]= 8;
			attribs[attribcount++]= AGL_PIXEL_SIZE;
			attribs[attribcount++]= 32;
		}
	}
	
    // Support for 3D rendering requested?
    if (PsychPrefStateGet_3DGfx()) {
        // Yes. Allocate a 24-Bit depth and 8-Bit stencilbuffer for this purpose:
        attribs[attribcount++]=kCGLPFADepthSize;
        attribs[attribcount++]=24;
        attribs[attribcount++]=kCGLPFAStencilSize;
        attribs[attribcount++]=8;
		// Alloc an accumulation buffer as well?
		if (PsychPrefStateGet_3DGfx() & 2) {
			// Yes: Alloc accum buffer, request 64 bpp, aka 16 bits integer per color component if possible:
			if (useCocoa || (!useAGL && !AGLForFullscreen)) {
				attribs[attribcount++]=kCGLPFAAccumSize;
				attribs[attribcount++]=64;
			}
			else {
				attribs[attribcount++]=AGL_ACCUM_RED_SIZE;
				attribs[attribcount++]=16;
				attribs[attribcount++]=AGL_ACCUM_GREEN_SIZE;
				attribs[attribcount++]=16;
				attribs[attribcount++]=AGL_ACCUM_BLUE_SIZE;
				attribs[attribcount++]=16;
				attribs[attribcount++]=AGL_ACCUM_ALPHA_SIZE;
				attribs[attribcount++]=16;
			}
		}
    }
    if(numBuffers>=2){
        // Enable double-buffering:
        attribs[attribcount++]=kCGLPFADoubleBuffer;
        if ((conserveVRAM & kPsychDisableAUXBuffers) == 0) {
            // MK: Allocate one or two (mono vs. stereo) AUX buffers for new "don't clear" mode of Screen('Flip'):
            // Not clearing the framebuffer after "Flip" is implemented by storing a backup-copy of
            // the backbuffer to AUXs before flip and restoring the content from AUXs after flip.
            attribs[attribcount++]=kCGLPFAAuxBuffers;
            attribs[attribcount++]=(stereomode==kPsychOpenGLStereo || stereomode==kPsychCompressedTLBRStereo || stereomode==kPsychCompressedTRBLStereo) ? 2 : 1;
        }
    }

    // MK: Stereo display support: If stereo display output is requested with OpenGL native stereo,
    // we request a stereo-enabled rendering context.
    if(stereomode==kPsychOpenGLStereo) {
        attribs[attribcount++]=kCGLPFAStereo;
    }

    // Multisampled Anti-Aliasing requested?
    if (windowRecord->multiSample > 0) {
      // Request a multisample buffer:
      attribs[attribcount++]= kCGLPFASampleBuffers;
      attribs[attribcount++]= 1;
      // Request at least multiSample samples per pixel:
      attribs[attribcount++]= kCGLPFASamples;
      attribs[attribcount++]= windowRecord->multiSample;
    }

    // Finalize attribute array with NULL.
    attribs[attribcount]=(CGLPixelFormatAttribute)NULL;

    // Init to zero:
    windowRecord->targetSpecific.pixelFormatObject = NULL;
	windowRecord->targetSpecific.glusercontextObject = NULL;
	windowRecord->targetSpecific.glswapcontextObject = NULL;

	if (useCocoa || (!useAGL && !AGLForFullscreen)) {
		// Context setup for CGL (non-windowed, non-multiscreen mode):
		
		// First try in choosing a matching format for multisample mode:
		if (windowRecord->multiSample > 0) {
			error=CGLChoosePixelFormat(attribs, &(windowRecord->targetSpecific.pixelFormatObject), &numVirtualScreens);
			if (windowRecord->targetSpecific.pixelFormatObject==NULL && windowRecord->multiSample > 0) {
				// Failed. Probably due to too demanding multisample requirements: Lets lower them...
				for (i=0; i<attribcount && attribs[i]!=kCGLPFASamples; i++);
				while (windowRecord->targetSpecific.pixelFormatObject==NULL && windowRecord->multiSample > 0) {
					attribs[i+1]--;
					windowRecord->multiSample--;
					error=CGLChoosePixelFormat(attribs, &(windowRecord->targetSpecific.pixelFormatObject), &numVirtualScreens);
				}
				if (windowRecord->multiSample == 0 && windowRecord->targetSpecific.pixelFormatObject==NULL) {
					for (i=0; i<attribcount && attribs[i]!=kCGLPFASampleBuffers; i++);
					attribs[i+1]=0;
				}
			}
		}
		
		// Choose a matching display configuration and create the window and rendering context:
		// If one of these two fails, then the installed gfx hardware is not good enough to satisfy our
		// requirements, or we have massive ressource shortage in the system. -> Screwed up anyway, so we abort.
		if (windowRecord->targetSpecific.pixelFormatObject==NULL) error=CGLChoosePixelFormat(attribs, &(windowRecord->targetSpecific.pixelFormatObject), &numVirtualScreens);
		if (error) {
			printf("\nPTB-ERROR[ChoosePixelFormat failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", CGLErrorString(error));
			return(FALSE);
		}
		
		// Create an OpenGL rendering context with the selected pixelformat: Share its ressources with 'slaveWindow's context, if slaveWindow is non-NULL:
		error=CGLCreateContext(windowRecord->targetSpecific.pixelFormatObject, ((windowRecord->slaveWindow) ? windowRecord->slaveWindow->targetSpecific.contextObject : NULL), &(windowRecord->targetSpecific.contextObject));
		if (error) {
			printf("\nPTB-ERROR[ContextCreation failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", CGLErrorString(error));
			return(FALSE);
		}
		
		// Enable the OpenGL rendering context associated with our window:
		error=CGLSetCurrentContext(windowRecord->targetSpecific.contextObject);
		if (error) {
			printf("\nPTB-ERROR[SetCurrentContext failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", CGLErrorString(error));
			return(FALSE);
		}
		
        if (!useCocoa) {
            // Switch to fullscreen display: We don't support windowed display on OS-X
            error=CGLSetFullScreen(windowRecord->targetSpecific.contextObject);
            if (error) {
                printf("\nPTB-ERROR[CGLSetFullScreen failed: %s]:The specified display may not support the current color depth -\nPlease switch to 'Millions of Colors' in Display Settings.\n\n", CGLErrorString(error));
                CGLSetCurrentContext(NULL);
                return(FALSE);
            }
        }
        
		// NULL-out the AGL context field, just for safety...
		windowRecord->targetSpecific.deviceContext = NULL;
	}
	else {
		// Context setup for AGL (windowed or multiscreen mode or fullscreen mode):
        // Disabled on 64-Bit.
        #ifndef __LP64__
		// Fullscreen mode?
		if (AGLForFullscreen) {
			// Need to map to target display device:
			err = DMGetGDeviceByDisplayID(cgDisplayID, &gdhDisplay, FALSE);
			if (noErr != err) {
				printf("\nPTB-ERROR[DMGetGDeviceByDisplayID failed: %s]: Can't map screenId to valid AGL display device id.\n\n", aglErrorString(aglGetError()));
				return(FALSE);
			}
			
			// And set mapping:
			gdhDisplayPtr = &gdhDisplay;
			gdhcount = 1;
		}

		// First try in choosing a matching format for multisample mode:
		if (windowRecord->multiSample > 0) {
			pf = aglChoosePixelFormat(gdhDisplayPtr, gdhcount, (GLint*) attribs);
			if (pf==NULL && windowRecord->multiSample > 0) {
				// Failed. Probably due to too demanding multisample requirements: Lets lower them...
				for (i=0; i<attribcount && attribs[i]!=kCGLPFASamples; i++);
				while (pf==NULL && windowRecord->multiSample > 0) {
					attribs[i+1]--;
					windowRecord->multiSample--;
					pf = aglChoosePixelFormat(gdhDisplayPtr, gdhcount, (GLint*) attribs);
				}
				
				if (windowRecord->multiSample == 0 && pf==NULL) {
					for (i=0; i<attribcount && attribs[i]!=kCGLPFASampleBuffers; i++);
					attribs[i+1]=0;
				}
			}
		}

		// Choose matching pixelformat:
		pf = aglChoosePixelFormat(gdhDisplayPtr, gdhcount, (GLint*) attribs);
		if (pf == NULL) {
			printf("\nPTB-ERROR[aglChoosePixelFormat failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
			DisposeWindow(carbonWindow); 
			return(FALSE);
		}
		
		// Create OpenGL rendering context for format: Share its ressources with 'slaveWindow' if that's non-NULL:
		glcontext = aglCreateContext(pf, (windowRecord->slaveWindow) ? (AGLContext)(windowRecord->slaveWindow->targetSpecific.deviceContext) : NULL);
		if (glcontext == NULL) {
			printf("\nPTB-ERROR[aglCreateContext failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
			aglDestroyPixelFormat(pf);
			DisposeWindow(carbonWindow);
			return(FALSE);
		}
		
		// Make it the current rendering context:
		if (!aglSetCurrentContext(glcontext)) {
			printf("\nPTB-ERROR[aglSetCurrentContext failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
			aglDestroyContext(glcontext);
			aglDestroyPixelFormat(pf);
			DisposeWindow(carbonWindow);
			return(FALSE);
		}
		
		aglbufferid = (int) windowRecord->windowIndex;
		if (!aglSetInteger(glcontext, AGL_BUFFER_NAME, &aglbufferid)) {
			printf("\nPTB-ERROR[aglSetInteger failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
			aglSetCurrentContext(NULL);
			aglDestroyContext(glcontext);
			aglDestroyPixelFormat(pf);
			DisposeWindow(carbonWindow);
			return(FALSE);
		}

		if ((windowLevel >= 1000) && (windowLevel < 2000)) {
			// For windowLevels between 1000 and 1999, make the window background transparent, so standard GUI
			// would be visible, wherever nothing is drawn, i.e., where alpha channel is zero:
			i = 0;
			aglSetInteger(glcontext, AGL_SURFACE_OPACITY, &i);
			
			// Levels 1000 - 1499 and 1500 to 1999 map to a master opacity level of 0.0 - 1.0:
			SetWindowAlpha(carbonWindow, ((float) (windowLevel % 500)) / 499.0);
		}
		
		if (AGLForFullscreen) {
			// Attach context to our fullscreen device:
			if (!aglSetFullScreen(glcontext, 0, 0, 0, 0)) {
				printf("\nPTB-ERROR[aglSetFullScreen failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
				aglSetCurrentContext(NULL);
				aglDestroyContext(glcontext);
				aglDestroyPixelFormat(pf);
				return(FALSE);
			}
		}
		else {
			// Attach context to our Carbon windows drawable area:
			if (!aglSetDrawable(glcontext, GetWindowPort(carbonWindow))) {
				printf("\nPTB-ERROR[aglSetDrawable failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
				aglSetCurrentContext(NULL);
				aglDestroyContext(glcontext);
				aglDestroyPixelFormat(pf);
				DisposeWindow(carbonWindow);
				return(FALSE);
			}
		}
		
		// Ok, theoretically we should have a fully functional OpenGL rendering context, attached to a fully
		// functional and visible onscreen window at this point.
		
		// Query CGL context and pixelformat for the AGL context and pixelformat, so all further code can use them:
		if (!aglGetCGLPixelFormat(pf, (void**) &(windowRecord->targetSpecific.pixelFormatObject))) {
			printf("\nPTB-ERROR[aglGetCGLPixelFormat failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
			aglSetCurrentContext(NULL);
			aglDestroyContext(glcontext);
			aglDestroyPixelFormat(pf);
			DisposeWindow(carbonWindow);
			return(FALSE);
		}
		
		if (!aglGetCGLContext(glcontext, (void**) &(windowRecord->targetSpecific.contextObject))) {
			printf("\nPTB-ERROR[aglGetCGLContext failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
			aglSetCurrentContext(NULL);
			aglDestroyContext(glcontext);
			aglDestroyPixelFormat(pf);
			DisposeWindow(carbonWindow);
			return(FALSE);
		}
				
		// Store AGLContext for reference...
		windowRecord->targetSpecific.deviceContext = (void*) glcontext;
		
		if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: AGL context setup finished..\n");
        
        #endif // 32-Bit AGL setup path.
        
        // 64-Bit setup path, using Cocoa + NSOpenGLContext via Objective-C.
        // This is the legacy path for OSX 10.5 "Leopard", which only provides
        // basic functionality:
        #ifdef __LP64__
        if (PsychCocoaSetupAndAssignLegacyOpenGLContext(carbonWindow, windowRecord)) {
            printf("\nPTB-ERROR[Cocoa Legacy-OpenGL setup failed]: Setup failed for unknown reasons.\n\n");
            PsychCocoaDisposeWindow(windowRecord);
            return(FALSE);
        }
        #endif
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
		if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Using GLEW version %s for automatic detection of OpenGL extensions...\n", glewGetString(GLEW_VERSION));
	}
	
    // Enable multisampling if it was requested:
    if (windowRecord->multiSample > 0) glEnable(GL_MULTISAMPLE);
    
    // Double-check double buffer support:
    if(numBuffers>=2){
        //This doesn't work.  GL thinks that there are double buffers when we fail to get that because
        //their was insufficient video memory to open a back buffer.  
        isDoubleBuffer=false;
        glGetBooleanv(GL_DOUBLEBUFFER, &isDoubleBuffer);
        if(!isDoubleBuffer){
			if (useCocoa || !useAGL) {
				CGLDestroyPixelFormat(windowRecord->targetSpecific.pixelFormatObject);
				CGLSetCurrentContext(NULL);
				CGLClearDrawable(windowRecord->targetSpecific.contextObject ) ;
			}
			else {
				aglSetCurrentContext(NULL);
				aglDestroyContext(glcontext);
				aglDestroyPixelFormat(pf);
				DisposeWindow(carbonWindow);
			}
            return(FALSE);
        }
    }
    
    // External 3D graphics support enabled?
	if (PsychPrefStateGet_3DGfx()) {
		// Yes. We need to create an extra OpenGL rendering context for the external
		// OpenGL code to provide optimal state-isolation. The context shares all
		// heavyweight ressources likes textures, FBOs, VBOs, PBOs, display lists and
		// starts off as an identical copy of PTB's context as of here.

		if (useCocoa || !useAGL) {
			error=CGLCreateContext(windowRecord->targetSpecific.pixelFormatObject, windowRecord->targetSpecific.contextObject, &(windowRecord->targetSpecific.glusercontextObject));
			if (error) {
				printf("\nPTB-ERROR[UserContextCreation failed: %s]: Creating a private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n", CGLErrorString(error));
				return(FALSE);
			}
            
            if (!useCocoa) {
                // Attach it to our onscreen drawable:
                error=CGLSetFullScreen(windowRecord->targetSpecific.glusercontextObject);
                if (error) {
                    printf("\nPTB-ERROR[CGLSetFullScreen for user context failed: %s]: Attaching private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n", CGLErrorString(error));
                    CGLSetCurrentContext(NULL);
                    return(FALSE);
                }
            }
            
			// Copy full state from our main context:
			error = CGLCopyContext(windowRecord->targetSpecific.contextObject, windowRecord->targetSpecific.glusercontextObject, GL_ALL_ATTRIB_BITS);
			if (error) {
				printf("\nPTB-ERROR[CGLCopyContext for user context failed: %s]: Copying state to private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n", CGLErrorString(error));
				CGLSetCurrentContext(NULL);
				return(FALSE);
			}
		}
		else {
            // Disabled on 64-Bit
            #ifndef __LP64__
			AGLContext usercontext = NULL;
			usercontext = aglCreateContext(pf, glcontext);
			if (usercontext == NULL) {
				printf("\nPTB-ERROR[AGL-UserContextCreation failed: %s]: Creating a private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n", aglErrorString(aglGetError()));
				// Ok, this is dirty, but better than nothing...
				DisposeWindow(carbonWindow);
				return(FALSE);
			}

			if (!aglSetInteger(usercontext, AGL_BUFFER_NAME, &aglbufferid)) {
				printf("\nPTB-ERROR[aglSetInteger for user context failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
				DisposeWindow(carbonWindow);
				return(FALSE);
			}
			
			if (AGLForFullscreen) {
				// Attach context to our fullscreen device:
				if (!aglSetFullScreen(usercontext, 0, 0, 0, 0)) {
					printf("\nPTB-ERROR[aglSetFullScreen for user context failed: %s]: Attaching private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n",  aglErrorString(aglGetError()));
					// Ok, this is dirty, but better than nothing...
					return(FALSE);
				}
			}
			else {
				// Attach it to our onscreen drawable:
				if (!aglSetDrawable(usercontext, GetWindowPort(carbonWindow))) {
					printf("\nPTB-ERROR[aglSetDrawable for user context failed: %s]: Attaching private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n",  aglErrorString(aglGetError()));
					// Ok, this is dirty, but better than nothing...
					DisposeWindow(carbonWindow);
					return(FALSE);
				}
			}
			
			// Copy full state from our main context:
			if (!aglCopyContext(glcontext, usercontext, GL_ALL_ATTRIB_BITS)) {
				printf("\nPTB-ERROR[aglCopyContext for user context failed: %s]: Copying state to private OpenGL context for Matlab OpenGL failed for unknown reasons.\n\n",  aglErrorString(aglGetError()));
				// Ok, this is dirty, but better than nothing...
				DisposeWindow(carbonWindow);
				return(FALSE);
			}

			// Retrieve CGL context for userspace context:
			if (!aglGetCGLContext(usercontext, (void**) &(windowRecord->targetSpecific.glusercontextObject))) {
				printf("\nPTB-ERROR[aglGetCGLContext failed: %s]: Getting CGL userspace context for Matlab OpenGL failed for unknown reasons.\n\n", aglErrorString(aglGetError()));
				DisposeWindow(carbonWindow);
				return(FALSE);
			}
            #endif
		}
	}

    // Create glswapcontextObject - An OpenGL context for exclusive use by parallel
    // background threads, e.g., our thread for async flip operations:
    if (TRUE)
    {
		if (useCocoa || !useAGL) {
			error=CGLCreateContext(windowRecord->targetSpecific.pixelFormatObject, windowRecord->targetSpecific.contextObject, &(windowRecord->targetSpecific.glswapcontextObject));
			if (error) {
				printf("\nPTB-ERROR[SwapContextCreation failed: %s]: Creating a private OpenGL context for async-bufferswaps failed.\n\n", CGLErrorString(error));
				return(FALSE);
			}
            
            if (!useCocoa) {
                // Attach it to our onscreen drawable:
                error=CGLSetFullScreen(windowRecord->targetSpecific.glswapcontextObject);
                if (error) {
                    printf("\nPTB-ERROR[CGLSetFullScreen for swapcontext failed: %s]: Attaching private OpenGL context async-bufferswaps failed.\n\n", CGLErrorString(error));
                    CGLSetCurrentContext(NULL);
                    return(FALSE);
                }
            }
            
			// Copy full state from our main context:
			error = CGLCopyContext(windowRecord->targetSpecific.contextObject, windowRecord->targetSpecific.glswapcontextObject, GL_ALL_ATTRIB_BITS);
			if (error) {
				printf("\nPTB-ERROR[CGLCopyContext for swapcontext failed: %s]: Copying state to private OpenGL context for async-bufferswaps failed.\n\n", CGLErrorString(error));
				CGLSetCurrentContext(NULL);
				return(FALSE);
			}
		}
		else {
            // Disabled on 64-Bit
            #ifndef __LP64__
			AGLContext usercontext = NULL;
			usercontext = aglCreateContext(pf, glcontext);
			if (usercontext == NULL) {
				printf("\nPTB-ERROR[AGL-SwapContextCreation failed: %s]: Creating a private OpenGL context for async-bufferswaps failed for unknown reasons.\n\n", aglErrorString(aglGetError()));
				// Ok, this is dirty, but better than nothing...
				DisposeWindow(carbonWindow);
				return(FALSE);
			}

			if (!aglSetInteger(usercontext, AGL_BUFFER_NAME, &aglbufferid)) {
				printf("\nPTB-ERROR[aglSetInteger for swapcontext failed: %s]:The specified display may not support double buffering and/or stereo output. There could be insufficient video memory\n\n", aglErrorString(aglGetError()));
				DisposeWindow(carbonWindow);
				return(FALSE);
			}
			
			if (AGLForFullscreen) {
				// Attach context to our fullscreen device:
				if (!aglSetFullScreen(usercontext, 0, 0, 0, 0)) {
					printf("\nPTB-ERROR[aglSetFullScreen for swapcontext failed: %s]: Attaching private OpenGL context for async-bufferswaps failed for unknown reasons.\n\n",  aglErrorString(aglGetError()));
					// Ok, this is dirty, but better than nothing...
					return(FALSE);
				}
			}
			else {
				// Attach it to our onscreen drawable:
				if (!aglSetDrawable(usercontext, GetWindowPort(carbonWindow))) {
					printf("\nPTB-ERROR[aglSetDrawable for swapcontext failed: %s]: Attaching private OpenGL context for async-bufferswaps failed for unknown reasons.\n\n",  aglErrorString(aglGetError()));
					// Ok, this is dirty, but better than nothing...
					DisposeWindow(carbonWindow);
					return(FALSE);
				}
			}
			
			// Copy full state from our main context:
			if (!aglCopyContext(glcontext, usercontext, GL_ALL_ATTRIB_BITS)) {
				printf("\nPTB-ERROR[aglCopyContext for swapcontext failed: %s]: Copying state to private OpenGL context for async-bufferswaps failed for unknown reasons.\n\n",  aglErrorString(aglGetError()));
				// Ok, this is dirty, but better than nothing...
				DisposeWindow(carbonWindow);
				return(FALSE);
			}

			// Retrieve CGL context for userspace context:
			if (!aglGetCGLContext(usercontext, (void**) &(windowRecord->targetSpecific.glswapcontextObject))) {
				printf("\nPTB-ERROR[aglGetCGLContext failed: %s]: Getting CGL userspace context for async-bufferswaps failed for unknown reasons.\n\n", aglErrorString(aglGetError()));
				DisposeWindow(carbonWindow);
				return(FALSE);
			}
            #endif
		}
    }
    
	// Ok, if we reached this point and AGL is used, we should store its onscreen Carbon window handle:
	if (useCocoa || (useAGL && !AGLForFullscreen)) {
		windowRecord->targetSpecific.windowHandle = carbonWindow;
	}
	else {
		windowRecord->targetSpecific.windowHandle = NULL;
	}

    // 64-Bit setup path, using Cocoa + NSOpenGLContext wrapped around already
    // setup CGLContext via Objective-C.
    // This is the modern path for OSX 10.6 "Snow Leopard" and later:
    #ifdef __LP64__
    if (useCocoa) {
        if (PsychCocoaSetupAndAssignOpenGLContextsFromCGLContexts(carbonWindow, windowRecord)) {
            printf("\nPTB-ERROR[Cocoa OpenGL setup failed]: Setup failed for unknown reasons.\n\n");
            PsychCocoaDisposeWindow(windowRecord);
            return(FALSE);
        }
    }
    #endif

    // Initialize a low-level mapping of Framebuffer device data structures into
    // our address space: Needed for additional timing checks:
    if ((PsychPrefStateGet_VBLTimestampingMode() > 0) && (screenRefCount[screenSettings->screenNumber] == 0)) {
        // Initialize to safe default:
        fbsharedmem[screenSettings->screenNumber].shmem = NULL;

        // Get access to Mach service port for the physical display device associated
        // with this onscreen window and open our own connection to the port:
        if ((kIOReturnSuccess == IOServiceOpen(CGDisplayIOServicePort(cgDisplayID), mach_task_self(), kIOFBSharedConnectType, &(fbsharedmem[screenSettings->screenNumber].connect))) ||
            (kIOReturnSuccess == IOServiceOpen(CGDisplayIOServicePort(CGMainDisplayID()), mach_task_self(), kIOFBSharedConnectType, &(fbsharedmem[screenSettings->screenNumber].connect)))) {
            // Connection established.

			// Map the slice of device memory into our VM space:
			if (kIOReturnSuccess != IOConnectMapMemory(fbsharedmem[screenSettings->screenNumber].connect, kIOFBCursorMemory, mach_task_self(),
													   (mach_vm_address_t *) &(fbsharedmem[screenSettings->screenNumber].shmem),
													   &(fbsharedmem[screenSettings->screenNumber].shmemSize), kIOMapAnywhere)) {
				// Mapping failed!
				fbsharedmem[screenSettings->screenNumber].shmem = NULL;
				if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Failed to gain access to kernel-level vbl handler [IOConnectMapMemory()] - Fallback path for time stamping won't be available.\n");
			}
			else {
				if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Connection to kernel-level vbl handler established (shmem = %p).\n",  fbsharedmem[screenSettings->screenNumber].shmem);
			}
        }
        else {
            if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Failed to gain access to kernel-level vbl handler [IOServiceOpen()] - Fallback path for time stamping won't be available.\n");
        }
        
        // If the mapping worked, we have a pointer to the driver memory in .shmem, otherwise we have NULL:
    }
    
    // Use of CoreVideo is needed on 10.7 and later due to brokeness of the old method (thanks Apple!):
    if ((osMajor > 10) || (osMinor >= 7)) {
        useCoreVideoTimestamping = TRUE;
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Broken Apple OS/X 10.7 or later detected: Using CoreVideo timestamping instead of precise vbl-irq timestamping.\n");
    } else {
        useCoreVideoTimestamping = FALSE;
    }

    if (useCoreVideoTimestamping && (fbsharedmem[screenSettings->screenNumber].shmem) && (NULL == cvDisplayLink[screenSettings->screenNumber])) {
        // Create and start a CVDisplayLink for this screen.
        if (kCVReturnSuccess != CVDisplayLinkCreateWithCGDisplay(cgDisplayID, &cvDisplayLink[screenSettings->screenNumber])) {
            if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Failed to create CVDisplayLink for screenId %i. This may impair VBL timestamping.\n", screenSettings->screenNumber);
        } else {
            // Assign dummy output callback, as this is mandatory to get the link up and running:
            CVDisplayLinkSetOutputCallback(cvDisplayLink[screenSettings->screenNumber], &PsychCVDisplayLinkOutputCallback, (void*) (long int) screenSettings->screenNumber);
            
            // Setup shared data structure and mutex:
            memset(&cvDisplayLinkData[screenSettings->screenNumber], 0, sizeof(cvDisplayLinkData[screenSettings->screenNumber]));
            PsychInitMutex(&(cvDisplayLinkData[screenSettings->screenNumber].mutex));

            // Start the link:
            if (kCVReturnSuccess != CVDisplayLinkStart(cvDisplayLink[screenSettings->screenNumber])) {
                // Failed to start: Release it again and report error:
                CVDisplayLinkRelease(cvDisplayLink[screenSettings->screenNumber]);
                cvDisplayLink[screenSettings->screenNumber] = NULL;

                // Teardown shared data structure and mutex:
                PsychDestroyMutex(&(cvDisplayLinkData[screenSettings->screenNumber].mutex));
                useCoreVideoTimestamping = FALSE;
                
                if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Failed to start CVDisplayLink for screenId %i. This may impair VBL timestamping.\n", screenSettings->screenNumber);
            }
            else {
                // Display link started: Report some stuff for the fun of it...
                if (PsychPrefStateGet_Verbosity() > 3) {
                    // Wait for 50 msecs before query of video refresh from display link to give it a chance to start up:
                    PsychWaitIntervalSeconds(0.050);
                    
                    printf("PTB-INFO: CVDisplayLink for screen %i created to work around the brokenness of Apple Mac OS/X 10.7 and later:\n", screenSettings->screenNumber);
                    printf("PTB-INFO: Video refresh interval as measured by CoreVideo display link: %f msecs.\n", (float) CVDisplayLinkGetActualOutputVideoRefreshPeriod(cvDisplayLink[screenSettings->screenNumber]) * 1000.0);
                    CVTime outLatency = CVDisplayLinkGetOutputVideoLatency(cvDisplayLink[screenSettings->screenNumber]);
                    printf("PTB-INFO: Video display output delay as reported by CoreVideo display link: %f msecs.\n", screenSettings->screenNumber, (float) (((double) outLatency.timeValue / (double) outLatency.timeScale) * 1000.0));
                }
            }
        }
    }

    // Retain reference of this window to its screen:
    screenRefCount[screenSettings->screenNumber]++;

    // Done.
    return(TRUE);
}

/*
    PsychOSGetVBLTimeAndCount()

    Returns absolute system time of last VBL and current total count of VBL interrupts since
    startup of gfx-system for the given screen. Returns a time of -1 and a count of 0 if this
    feature is unavailable on the given OS/Hardware configuration.
*/
double PsychOSGetVBLTimeAndCount(PsychWindowRecordType *windowRecord, psych_uint64* vblCount)
{
    unsigned int screenid = windowRecord->screenNumber;
    psych_uint64 refvblcount;
    double tnow, t1, t2, cvRefresh, cvTime = 0;

    // Should we use CoreVideo display link timestamping?
    if (useCoreVideoTimestamping && cvDisplayLink[screenid]) {
        // Yes: Retrieve data from our shared data structure:
        PsychGetAdjustedPrecisionTimerSeconds(&tnow);
        
        PsychLockMutex(&(cvDisplayLinkData[screenid].mutex));
        *vblCount = cvDisplayLinkData[screenid].vblCount;
        cvTime = cvDisplayLinkData[screenid].vblTimestamp;
        PsychUnlockMutex(&(cvDisplayLinkData[screenid].mutex));

        // Current time more than 0.9 video refresh durations after queried
        // vblank timestamps cvTime?
        cvRefresh = CVDisplayLinkGetActualOutputVideoRefreshPeriod(cvDisplayLink[screenid]);
        if ((tnow - cvTime) > (0.9 * cvRefresh)) {
            // Yes: It is more likely that we fetched a stale cvTime timestamp from
            // the previous vblank instead of the current one, because we got called
            // before the display link callback could execute for the current vblank.
            // Assume this is the case and correct the returned timestamp and count.
            // Add one video refresh duration and count to make it so:
            cvTime += cvRefresh;
            *vblCount = *vblCount + 1;
        }

        // If timestamp debugging is off, we're done:
        if (PsychPrefStateGet_Verbosity() <= 19) return(cvTime);
    }

    // Do we have a valid shared mapping?
    if (fbsharedmem[screenid].shmem) {
		// We query each value twice and repeat this double-query until both readings of
		// both variables show the same values. This because our read access to this kernel
		// data structure is not protected by a lock, so the kernel might modify the values
		// while we are reading them, causing inconsistent readings between the two values and
		// within a value. We use memory barriers to prevent compiler optimizations or cache
		// coherency issues:
		do {
			// Retrieve absolute count of vbls since startup:
			*vblCount = (psych_uint64) fbsharedmem[screenid].shmem->vblCount;
			OSMemoryBarrier();
			t1 = ((double) UnsignedWideToUInt64(AbsoluteToNanoseconds(fbsharedmem[screenid].shmem->vblTime))) / 1000000000.0;
			OSMemoryBarrier();
			PsychWaitIntervalSeconds(0.000250);
			refvblcount = (psych_uint64) fbsharedmem[screenid].shmem->vblCount;
			OSMemoryBarrier();
			t2 = ((double) UnsignedWideToUInt64(AbsoluteToNanoseconds(fbsharedmem[screenid].shmem->vblTime))) / 1000000000.0;
			OSMemoryBarrier();
		} while ((*vblCount != refvblcount) || (t1 != t2));

        // Diagnostic check of CV timestamps against our host timestamps:
        // See comments in PsychCVDisplayLinkOutputCallback() on what to expect here...
        if (useCoreVideoTimestamping && (PsychPrefStateGet_Verbosity() > 19)) {
                t2 = cvTime - t1;
                printf("PTB-DEBUG: Difference CoreVideoTimestamp - vblTimestamp = %lf msecs.\n", 1000.0 * t2);
        }

        // Retrieve absolute system time of last retrace, convert into PTB standard time system and return it:
		return(t1);
    }
    else {
        // Unsupported :(
        *vblCount = 0;
        return(-1);
    }
}

/*
    PsychOSOpenOffscreenWindow()
    
    Accept specifications for the offscreen window in the platform-neutral structures, convert to native CoreGraphics structures,
    create the surface, allocate a window record and record the window specifications and memory location there.
	
	TO DO:  We need to walk down the screen number and fill in the correct value for the benefit of TexturizeOffscreenWindow
*/
psych_bool PsychOSOpenOffscreenWindow(double *rect, int depth, PsychWindowRecordType **windowRecord)
{

    //PsychTargetSpecificWindowRecordType 	cgStuff;
    CGLPixelFormatAttribute 			attribs[5];
    //CGLPixelFormatObj					pixelFormatObj;
    GLint								numVirtualScreens;
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
}


void PsychOSCloseWindow(PsychWindowRecordType *windowRecord)
{    
    CGDirectDisplayID				cgDisplayID;

    // Disable rendering context:
    CGLSetCurrentContext(NULL);
 
	if (windowRecord->targetSpecific.windowHandle == NULL) {
		// Shutdown sequence for CGL, i.e., no AGL+Carbon mode:
		// MK: Hack, needed to work around a "screen corruption on shutdown" bug.
		// When closing stereo display windows, it sometimes leads to a completely
		// messed up and unusable display.
		if (PsychIsOnscreenWindow(windowRecord)) {
			PsychReleaseScreen(windowRecord->screenNumber);
			// Destroy onscreen window, detach context:
			CGLClearDrawable(windowRecord->targetSpecific.contextObject);
			if (windowRecord->targetSpecific.glusercontextObject) CGLClearDrawable(windowRecord->targetSpecific.glusercontextObject);
			if (windowRecord->targetSpecific.glswapcontextObject) CGLClearDrawable(windowRecord->targetSpecific.glswapcontextObject);
			PsychCaptureScreen(windowRecord->screenNumber);
		}
	}

	// Destroy pixelformat object:
    CGLDestroyPixelFormat(windowRecord->targetSpecific.pixelFormatObject);
    
	// Destroy rendering context:
    #ifdef __LP64__
        // 64-Bit Path: We're on OSX 10.5 or later. Use new-style CGLReleaseContext(), which is needed for our
        // Cocoa windowed mode with NSOpenGLContext wrapped around CGLContext:

        // printf("PRERELEASE: Refcounts: mccgl=%i sccgl=%i\n", CGLGetContextRetainCount(windowRecord->targetSpecific.contextObject), CGLGetContextRetainCount(windowRecord->targetSpecific.glswapcontextObject));
        CGLReleaseContext(windowRecord->targetSpecific.contextObject);
        if (windowRecord->targetSpecific.glusercontextObject) CGLReleaseContext(windowRecord->targetSpecific.glusercontextObject);
        if (windowRecord->targetSpecific.glswapcontextObject) CGLReleaseContext(windowRecord->targetSpecific.glswapcontextObject);
    #else
        // 32-Bit Path: Use good ol' CGLDestroyContext():
        CGLDestroyContext(windowRecord->targetSpecific.contextObject);
        if (windowRecord->targetSpecific.glusercontextObject) CGLDestroyContext(windowRecord->targetSpecific.glusercontextObject);
        if (windowRecord->targetSpecific.glswapcontextObject) CGLDestroyContext(windowRecord->targetSpecific.glswapcontextObject);
    #endif
    
    // Disable low-level mapping of framebuffer cursor memory:
    if (fbsharedmem[windowRecord->screenNumber].shmem && (screenRefCount[windowRecord->screenNumber] == 1)) {
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Releasing shared memory mapping for screen %i.\n", windowRecord->screenNumber);

        IOConnectUnmapMemory(fbsharedmem[windowRecord->screenNumber].connect, kIOFBCursorMemory, mach_task_self(), (vm_address_t) fbsharedmem[windowRecord->screenNumber].shmem);
        fbsharedmem[windowRecord->screenNumber].shmem = NULL;
        
        // Close the service port:
        IOServiceClose(fbsharedmem[windowRecord->screenNumber].connect);
        
        // Shutdown and release CVDisplayLink for this windows screen, if any:
        if (cvDisplayLink[windowRecord->screenNumber]) {
            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Releasing CVDisplayLink for screen %i.\n", windowRecord->screenNumber);
            
            if (CVDisplayLinkIsRunning(cvDisplayLink[windowRecord->screenNumber])) CVDisplayLinkStop(cvDisplayLink[windowRecord->screenNumber]);
            PsychYieldIntervalSeconds(0.020);
            CVDisplayLinkRelease(cvDisplayLink[windowRecord->screenNumber]);
            cvDisplayLink[windowRecord->screenNumber] = NULL;

            // Teardown shared data structure and mutex:
            PsychDestroyMutex(&(cvDisplayLinkData[windowRecord->screenNumber].mutex));
        }
    }

    // Release reference of this window to its screen:
    screenRefCount[windowRecord->screenNumber]--;
    
    #ifdef __LP64__
        // Destroy Cocoa onscreen window, if any:
        if (windowRecord->targetSpecific.windowHandle) PsychCocoaDisposeWindow(windowRecord);
    #else
        // Destroy Carbon onscreen window, if any:
        if (windowRecord->targetSpecific.windowHandle) DisposeWindow(windowRecord->targetSpecific.windowHandle);
    #endif
    
    // printf("POSTWINDISPOSE: Refcounts: mccgl=%i sccgl=%i\n", CGLGetContextRetainCount(windowRecord->targetSpecific.contextObject), CGLGetContextRetainCount(windowRecord->targetSpecific.glswapcontextObject));

    return;
}

/* PsychOSGetSwapCompletionTimestamp()
 *
 * Retrieve a very precise timestamp of doublebuffer swap completion by means
 * of OS specific facilities. This function is optional. If the underlying
 * OS/drier/GPU combo doesn't support a high-precision, high-reliability method
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
	// Unsupported on OS/X:
	return(-1);
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
	// On OS/X this function is unsupported:
	return(-1);
}

/*
 * PsychOSFlipWindowBuffers() -- OS-X swapbuffers call.
 */
void PsychOSFlipWindowBuffers(PsychWindowRecordType *windowRecord)
{	
	CGLError			cglerr;
    psych_bool oldStyle = (PsychPrefStateGet_ConserveVRAM() & kPsychUseOldStyleAsyncFlips) ? TRUE : FALSE;
    
	// Execute OS neutral bufferswap code first:
	PsychExecuteBufferSwapPrefix(windowRecord);
	
    // Trigger the "Front <-> Back buffer swap (flip) (on next vertical retrace)":
    if ((cglerr = CGLFlushDrawable((oldStyle || PsychIsMasterThread()) ? windowRecord->targetSpecific.contextObject : windowRecord->targetSpecific.glswapcontextObject))) {
		// Failed! This is an internal OpenGL/CGL error. We can't do anything about it, just report it:
		printf("PTB-ERROR: Doublebuffer-Swap failed (probably during 'Flip')! Internal OpenGL subsystem/driver error: %s. System misconfigured or driver/operating system bug?!?\n", CGLErrorString(cglerr));
	}
}

/*
    PsychOSSetGLContext()
    
    Set the window to which GL drawing commands are sent.  
*/
void PsychOSSetGLContext(PsychWindowRecordType *windowRecord)
{
    // Setup new context if it isn't already setup. -> Avoid redundant context switch.
    if (CGLGetCurrentContext() != windowRecord->targetSpecific.contextObject) {
		if (CGLGetCurrentContext() != NULL) {
			// We need to glFlush the old context before switching, otherwise race-conditions may occur:
			glFlush();
			
			// Need to unbind any FBO's in old context before switch, otherwise bad things can happen...
			if (glBindFramebufferEXT) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		}
		
		// Switch to new context:
		CGLSetCurrentContext(windowRecord->targetSpecific.contextObject);		
    }
}

/* Same as PsychOSSetGLContext() but for selecting userspace rendering context,
 * optionally copying state from PTBs context.
 */
void PsychOSSetUserGLContext(PsychWindowRecordType *windowRecord, psych_bool copyfromPTBContext)
{
	// Child protection:
	if (windowRecord->targetSpecific.glusercontextObject == NULL) PsychErrorExitMsg(PsychError_user, "GL Userspace context unavailable! Call InitializeMatlabOpenGL *before* Screen('OpenWindow')!");
	
	if (copyfromPTBContext) {
		// Syncing of external contexts state with PTBs internal state requested. Do it:
		CGLSetCurrentContext(NULL);
		CGLCopyContext(windowRecord->targetSpecific.contextObject, windowRecord->targetSpecific.glusercontextObject, GL_ALL_ATTRIB_BITS);
	}
	
    // Setup new context if it isn't already setup. -> Avoid redundant context switch.
    if (CGLGetCurrentContext() != windowRecord->targetSpecific.glusercontextObject) {
        CGLSetCurrentContext(windowRecord->targetSpecific.glusercontextObject);
    }
}

/*
    PsychOSUnsetGLContext()
    
    Clear the drawing context.  
*/
void PsychOSUnsetGLContext(PsychWindowRecordType *windowRecord)
{
	if (CGLGetCurrentContext() != NULL) {
		// We need to glFlush the old context before switching, otherwise race-conditions may occur:
		glFlush();
		
		// Need to unbind any FBO's in old context before switch, otherwise bad things can happen...
		if (glBindFramebufferEXT) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	}

	// Detach totally:
    CGLSetCurrentContext(NULL);
}

/* PsychOSSetVBLSyncLevel - Enable or disable synchronization of bufferswaps to
 * onset of vertical retrace. This is the OS-X version...
 */
void PsychOSSetVBLSyncLevel(PsychWindowRecordType *windowRecord, int swapInterval)
{
    CGLError	error;
    GLint myinterval = (GLint) swapInterval;
    psych_bool oldStyle = (PsychPrefStateGet_ConserveVRAM() & kPsychUseOldStyleAsyncFlips) ? TRUE : FALSE;
    	
	// Store new setting also in internal helper variable, e.g., to allow workarounds to work:
	windowRecord->vSynced = (swapInterval > 0) ? TRUE : FALSE;
	
    error=CGLSetParameter((oldStyle || PsychIsMasterThread()) ? windowRecord->targetSpecific.contextObject : windowRecord->targetSpecific.glswapcontextObject, kCGLCPSwapInterval, &myinterval);
    if (error) {
        if (PsychPrefStateGet_Verbosity()>1) printf("\nPTB-WARNING: FAILED to %s synchronization to vertical retrace!\n\n", (swapInterval>0) ? "enable" : "disable");
    }

    error=CGLGetParameter((oldStyle || PsychIsMasterThread()) ? windowRecord->targetSpecific.contextObject : windowRecord->targetSpecific.glswapcontextObject, kCGLCPSwapInterval, &myinterval);
    if (error || (myinterval != (long) swapInterval)) {
        if (PsychPrefStateGet_Verbosity()>1) printf("\nPTB-WARNING: FAILED to %s synchronization to vertical retrace (System ignored setting)!\n\n", (swapInterval>0) ? "enable" : "disable");
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
	// On OS/X the situation is simple. This OS doesn't support framelock/swaplock at
	// all on any GPU:
	return(FALSE);
}

// Perform OS specific processing of Window events:
void PsychOSProcessEvents(PsychWindowRecordType *windowRecord, int flags)
{
	Rect globalBounds;

	// Trigger event queue dispatch processing for GUI windows:
	if (windowRecord == NULL) {
		// No op, so far...
		return;
	}
	
	// GUI windows need to behave GUIyee:
	if ((windowRecord->specialflags & kPsychGUIWindow) && PsychIsOnscreenWindow(windowRecord)) {
		// Update windows rect and globalrect, based on current size and location:
		GetWindowBounds(windowRecord->targetSpecific.windowHandle, kWindowContentRgn, &globalBounds);
		PsychMakeRect(windowRecord->globalrect, globalBounds.left, globalBounds.top, globalBounds.right, globalBounds.bottom);
		PsychNormalizeRect(windowRecord->globalrect, windowRecord->rect);
		PsychSetupClientRect(windowRecord);
		PsychSetupView(windowRecord, FALSE);
	}
}
