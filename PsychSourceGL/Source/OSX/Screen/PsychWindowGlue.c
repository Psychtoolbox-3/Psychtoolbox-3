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
                 10/19/05               awi             Cast NULL to CGLPixelFormatAttribute type to make the compiler happy.
                 01/02/05               mk              Modified to only contain the OS-X specific code. All OS independent code has been moved to
                                                        Common/Screen/PsychWindowSupport.c
 
	DESCRIPTION:
	
		Functions in this file comprise an abstraction layer for probing and controlling window state, except for window content.  
		
		Each C function which implements a particular Screen subcommand should be platform neutral.  For example, the source to SCREENPixelSizes() 
		should be platform-neutral, despite that the calls in OS X and Windows to detect available pixel sizes are different.  The platform 
		specificity is abstracted out in C files which end it "Glue", for example PsychScreenGlue, PsychWindowGlue, PsychWindowTextClue.

	NOTES:
	
	TO DO: 

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

#include "PsychCocoaGlue.h"

static struct {
    psych_mutex         mutex;
    double              vblTimestamp;
    psych_uint64        vblCount;
} cvDisplayLinkData[kPsychMaxPossibleDisplays];

static CVDisplayLinkRef cvDisplayLink[kPsychMaxPossibleDisplays] = { NULL };
static int screenRefCount[kPsychMaxPossibleDisplays] = { 0 };

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
    
    // Extra guard against shutdown races:
    if (NULL == cvDisplayLink[screenId]) return(kCVReturnSuccess);
    
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
    
    -To tell the caller to clean up PsychOSOpenOnscreenWindow returns FALSE if we fail to open the window. It 
    would be better to just issue an PsychErrorExit() and have that clean up everything allocated outside of
    PsychOpenOnscreenWindow().
    
    MK: The new option 'stereomode' allows selection of stereo display instead of mono display:
    0 (default) ==  Old behaviour -> Monoscopic rendering context.
    >0          ==  Stereo display, where the number defines the type of stereo algorithm to use.
    =1          ==  Use OpenGL built-in stereo by creating a context/window with left- and right backbuffer. This is
                    the only mode of interest here, as it requires use of a stereo capable OpenGL pixelformat. All other
                    stereo modes are implemented by PTB itself in a platform independent manner on top of a standard mono
                    context.
*/
psych_bool PsychOSOpenOnscreenWindow(PsychScreenSettingsType *screenSettings, PsychWindowRecordType *windowRecord, int numBuffers, int stereomode, int conserveVRAM)
{
    CGOpenGLDisplayMask             displayMask;
    CGLError                        error;
    CGDirectDisplayID               cgDisplayID;
    CGLPixelFormatAttribute         attribs[40];
    int                             attribcount;
    GLint                           numVirtualScreens;
    GLenum                          glerr;
    PsychRectType                   screenrect;
    int                             i;
    int                             windowLevel;
    long                            scw, sch;
    void*                           cocoaWindow = NULL;

    // Map screen number to physical display handle cgDisplayID:
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenSettings->screenNumber);
    displayMask = CGDisplayIDToOpenGLDisplayMask(cgDisplayID);

    // NULL-out Cocoa window handle, so this is well-defined in case of error:
    windowRecord->targetSpecific.windowHandle = NULL;

    // Retrieve windowLevel, an indicator of where non-CGL / non-fullscreen windows should
    // be located wrt. to other windows. -2 = Allow regular window manager control of stacking
    // order, visibility etc., -1 = Invisible/hidden, 0 = Behind everything else, occluded by
    // everything else. 1 - 999 = At layer 'windowLevel' -> Occludes stuff on layers "below" it.
    // 1000 - 1999 = At highest level, but partially translucent / alpha channel allows to make
    // regions transparent. Range 1000 - 1499 = transparent for mouse and keyboard, alpha 0-99.99%,
    // 1500-1599 = opaque for mouse and keyboard, alpha 0-99.99%, 2000 or higher: Above everything,
    // fully opaque, occludes everything, a typical fullscreen onscreen window. 2000 is the default.
    windowLevel = PsychPrefStateGet_WindowShieldingLevel();

    // Window rect provided which has same size as screen?
    // We do not use windowed mode if the provided window rectangle either
    // matches the target screens rectangle (and therefore its exact size)
    // or its screens global rectangle. In such cases we use CGL for better
    // low level control and to exclude the desktop compositor from interfering:
    PsychGetScreenRect(screenSettings->screenNumber, screenrect);
    if (PsychMatchRect(screenrect, windowRecord->rect)) windowRecord->specialflags |= kPsychIsFullscreenWindow;

    PsychGetGlobalScreenRect(screenSettings->screenNumber, screenrect);
    if (PsychMatchRect(screenrect, windowRecord->rect)) windowRecord->specialflags |= kPsychIsFullscreenWindow;

    if ((windowRecord->specialflags & kPsychIsFullscreenWindow) && (PsychPrefStateGet_Verbosity() > 3)) {
        printf("PTB-INFO: Always using Cocoa for fullscreen windows to work around graphics driver bugs in OSX.\n");
        printf("PTB-INFO: Presentation timing precision is not yet known for this configuration on most machines. Check your results.\n");
    }

    // Display for fullscreen window not captured? Timing precision is unclear in this mode. In theory the compositor should disable
    // itself for fullscreen windows on modern OSX versions. If it really does that, who knows?
    if ((windowRecord->specialflags & kPsychIsFullscreenWindow) && (PsychPrefStateGet_ConserveVRAM() & kPsychUseAGLCompositorForFullscreenWindows)) {
        // Force a window rectangle that matches the global screen rectangle for that windows screen:
        PsychCopyRect(windowRecord->rect, screenrect);

        // Warn user about what's going on:
        if (PsychPrefStateGet_Verbosity()>1) printf("PTB-INFO: No display capture / compositor lockout for fullscreen window. Timing precision is unknown.\n");
    }

    if ((windowRecord->specialflags & kPsychGUIWindow) && (PsychPrefStateGet_Verbosity() > 3)) {
        printf("PTB-INFO: Onscreen window is configured as regular GUI window.\n");
    }

    if ((windowRecord->specialflags & kPsychIsFullscreenWindow) && (PsychPrefStateGet_Verbosity() > 1) && CGDisplayIsInMirrorSet(cgDisplayID)) {
        printf("PTB-WARNING: Onscreen window is displayed on a display with display mirroring enabled. This will likely ruin visual timing.\n");
    }

    // Create onscreen Cocoa window of requested position and size:
    if (PsychCocoaCreateWindow(windowRecord, windowLevel, &cocoaWindow)) {
        printf("\nPTB-ERROR[CreateNewWindow failed]: Failed to open Cocoa onscreen window\n\n");
        return(FALSE);
    }

    // Transparent window requested?
    if ((windowLevel >= 1000) && (windowLevel < 2000)) {
        // Setup of global window alpha value for transparency. This is premultiplied to
        // the individual per-pixel alpha values if transparency is enabled by Cocoa code.
        //
        // Levels 1000 - 1499 and 1500 to 1999 map to a master opacity level of 0.0 - 1.0:
        PsychCocoaSetWindowAlpha(cocoaWindow, ((float) (windowLevel % 500)) / 499.0);
    }

    // Show it! Unless a windowLevel of -1 requests hiding the window:
    if (windowLevel != -1) PsychCocoaShowWindow(cocoaWindow);

    // If usercode wants a black startup screen then we add a pause of 0.5 seconds here
    // before proceeding. This will avoid a white flash at window open time, which might
    // be something the user wanted to avoid. Why does this help or is needed at all?
    // Nobody knows, but this is Apples ridiculous toy OS, so why even ask such questions?
    if (PsychPrefStateGet_VisualDebugLevel() < 4)
        PsychYieldIntervalSeconds(0.5);

    // Level zero means: Place behind all other windows:
    if (windowLevel == 0) PsychCocoaSendBehind(cocoaWindow);

    // Levels 1 to 999 define window levels for the group of the window.
    // A level of -2 would leave this to the system:
    if (windowLevel > 0 && windowLevel < 1000) PsychCocoaSetWindowLevel(cocoaWindow, windowLevel);

    // Is the target display captured for a fullscreen window?
    if (PsychIsScreenCaptured(screenSettings->screenNumber)) {
        // Yes. Make sure our window is above the shielding window level:
        PsychCocoaSetWindowLevel(cocoaWindow, (int) CGShieldingWindowLevel());
    }

    // Store window handle in windowRecord:
    windowRecord->targetSpecific.windowHandle = cocoaWindow;

    // Store vblank startline aka true height of physical display screen in pixels:
    PsychGetScreenPixelSize(screenSettings->screenNumber, &scw, &sch);
    windowRecord->VBL_Startline = (int) sch;

    // Define pixelformat attributes for OpenGL contexts:

    // No pixelformat attribs to start with:
    attribcount = 0;

    attribs[attribcount++]=kCGLPFADisplayMask;
    attribs[attribcount++]=displayMask;

    // 10 bit per component integer framebuffer requested (10-10-10-2)?
    if (windowRecord->depth == 30) {
        // Request a 10 bit per color component framebuffer with 2 bit alpha channel:
        printf("PTB-INFO: Trying to enable 10 bpc, 30 bit integer framebuffer...\n");
        attribs[attribcount++]=kCGLPFANoRecovery;
        attribs[attribcount++]=kCGLPFAMinimumPolicy;
        attribs[attribcount++]=kCGLPFAColorSize;
        attribs[attribcount++]=10*3;
        attribs[attribcount++]=kCGLPFAAlphaSize;
        attribs[attribcount++]=2;
    }

    // 11 bit per component integer framebuffer requested (11-11-10-0)?
    if (windowRecord->depth == 33) {
        // Request a ~ 11 bit per color component framebuffer without alpha channel:
        printf("PTB-INFO: Trying to enable 11 bpc, 32 bit integer framebuffer...\n");
        attribs[attribcount++]=kCGLPFANoRecovery;
        attribs[attribcount++]=kCGLPFAMinimumPolicy;
        attribs[attribcount++]=kCGLPFAColorSize;
        attribs[attribcount++]=32;
        attribs[attribcount++]=kCGLPFAAlphaSize;
        attribs[attribcount++]=0;
    }

    // 16 bit per component integer framebuffer requested (16-16-16-16)?
    if (windowRecord->depth == 48) {
        // Request a 16 bit per color component framebuffer:
        printf("PTB-INFO: Trying to enable 16 bpc, 64 bit integer framebuffer...\n");
        attribs[attribcount++]=kCGLPFANoRecovery;
        attribs[attribcount++]=kCGLPFAMinimumPolicy;
        attribs[attribcount++]=kCGLPFAColorSize;
        attribs[attribcount++]=16*3;
        attribs[attribcount++]=kCGLPFAAlphaSize;
        attribs[attribcount++]=16;
    }

    // 16 bit per component, 64 bit framebuffer requested (16-16-16-16)?
    if (windowRecord->depth == 64) {
        // Request a floating point framebuffer in 16-bit half-float format, i.e., RGBA = 16 bits per component.
        printf("PTB-INFO: Trying to enable 16 bpc float framebuffer...\n");
        attribs[attribcount++]=kCGLPFAColorFloat;
        attribs[attribcount++]=kCGLPFAMinimumPolicy;
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
        attribs[attribcount++]=kCGLPFAMinimumPolicy;
        attribs[attribcount++]=kCGLPFAColorSize;
        attribs[attribcount++]=32*3;
        attribs[attribcount++]=kCGLPFAAlphaSize;
        attribs[attribcount++]=32;
    }

    // Possible to request use of the Apple floating point software renderer:
    if (conserveVRAM & kPsychUseSoftwareRenderer) {
        #ifndef kCGLRendererGenericFloatID
        #define kCGLRendererGenericFloatID    0x00020400
        #endif

        attribs[attribcount++]=AGL_RENDERER_ID;
        attribs[attribcount++]=kCGLRendererGenericFloatID;
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
            attribs[attribcount++]=kCGLPFAAccumSize;
            attribs[attribcount++]=64;
        }
    }

    if(numBuffers>=2){
        // Enable double-buffering:
        attribs[attribcount++]=kCGLPFADoubleBuffer;
        if ((conserveVRAM & kPsychDisableAUXBuffers) == 0) {
            // Allocate one or two (for mono vs. stereo display) AUX buffers for "don't clear" mode of Screen('Flip'):
            // Not clearing the framebuffer after "Flip" is implemented by storing a backup-copy of
            // the backbuffer to AUXs before flip and restoring the content from AUXs after flip.
            // Unless the imaging pipeline is active, which doesn't need AUX buffers due to internal
            // storage of fb content in its drawbufferFBO's:
            attribs[attribcount++]=kCGLPFAAuxBuffers;
            attribs[attribcount++]=(stereomode==kPsychOpenGLStereo || stereomode==kPsychCompressedTLBRStereo || stereomode==kPsychCompressedTRBLStereo) ? 2 : 1;
        }
    }

    // If stereo display output is requested with OpenGL native stereo, request a stereo-enabled rendering context.
    // This is deprecated since 10.11 El Capitan, and in fact does no longer work - OpenGL quad buffered stereo is
    // dead on 10.11 on all tested GPU's from Intel, NVidia, AMD.
    if(stereomode==kPsychOpenGLStereo) {
        attribs[attribcount++] = kCGLPFAStereo;
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
    attribs[attribcount++]=(CGLPixelFormatAttribute)NULL;

    // Init to zero:
    windowRecord->targetSpecific.pixelFormatObject = NULL;
    windowRecord->targetSpecific.glusercontextObject = NULL;
    windowRecord->targetSpecific.glswapcontextObject = NULL;

    // Try to find matching pixelformat:
    error = CGLChoosePixelFormat(attribs, &(windowRecord->targetSpecific.pixelFormatObject), &numVirtualScreens);

    // No valid pixelformat found and stereo format requested?
    if ((error || (windowRecord->targetSpecific.pixelFormatObject == NULL)) && (stereomode == kPsychOpenGLStereo)) {
        // Yep: Stereo may be the culprit. Remove the stereo attribute by overwriting it with something
        // that is essentially a no-op, specifically kCGLPFAAccelerated which is supported by all real
        // renderers that might end up in this code-path:
        for (i = 0; i < attribcount && attribs[i] != kCGLPFAStereo; i++);
        attribs[i] = kCGLPFAAccelerated;
        
        // Retry query of pixelformat without request for native OpenGL quad-buffered stereo. If we succeed, we're
        // sort of ok, as the higher-level code will fallback to stereomode kPsychFrameSequentialStereo - our own
        // homegrown frame-sequential stereo support, which may be good enough.
        error = CGLChoosePixelFormat(attribs, &(windowRecord->targetSpecific.pixelFormatObject), &numVirtualScreens);
        if (error || (windowRecord->targetSpecific.pixelFormatObject == NULL)) {
            windowRecord->targetSpecific.pixelFormatObject = NULL;
            printf("\nPTB-ERROR[ChoosePixelFormat failed: %s]: Disabling OpenGL native quad-buffered stereo did not help. Moving on...\n\n", CGLErrorString(error));
        }
    }

    // Now try if choosing a matching format for a lower multisample mode helps to get unstuck:
    if (windowRecord->multiSample > 0) {
        if (windowRecord->targetSpecific.pixelFormatObject==NULL && windowRecord->multiSample > 0) {
            // Failed. Probably due to too demanding multisample requirements: Lets lower them...
            for (i = 0; i < attribcount && attribs[i] != kCGLPFASamples; i++);
            while (windowRecord->targetSpecific.pixelFormatObject == NULL && windowRecord->multiSample > 0) {
                attribs[i+1]--;
                windowRecord->multiSample--;
                error = CGLChoosePixelFormat(attribs, &(windowRecord->targetSpecific.pixelFormatObject), &numVirtualScreens);
            }

            if (windowRecord->multiSample == 0 && windowRecord->targetSpecific.pixelFormatObject == NULL) {
                // Ok, multisampling is now at zero and we still don't succeed. Disable multisampling completely:
                for (i=0; i<attribcount && attribs[i]!=kCGLPFASampleBuffers; i++);
                attribs[i+1] = 0;
                printf("\nPTB-ERROR[ChoosePixelFormat failed: %s]: Disabling multisample anti-aliasing did not help. Moving on...\n\n", CGLErrorString(error));
            }
        }
    }

    // Try choosing a matching display configuration again and create the window and rendering context:
    // If one of these two fails, then the installed gfx hardware is not good enough to satisfy our
    // requirements, or we have massive ressource shortage in the system. -> Screwed up anyway, so we abort.
    if (windowRecord->targetSpecific.pixelFormatObject == NULL) error = CGLChoosePixelFormat(attribs, &(windowRecord->targetSpecific.pixelFormatObject), &numVirtualScreens);
    if (error) {
        printf("\nPTB-ERROR[ChoosePixelFormat failed: %s]: Reason unknown. There could be insufficient video memory or a driver malfunction. Giving up.\n\n", CGLErrorString(error));
        return(FALSE);
    }

    // Create an OpenGL rendering context with the selected pixelformat: Share its ressources with 'slaveWindow's context, if slaveWindow is non-NULL.
    // If slaveWindow is non-NULL here, then slaveWindow is typically another onscreen window. Therefore this establishes OpenGL resource sharing across
    // different onscreen windows in a session, e.g., for multi-display operation:
    error=CGLCreateContext(windowRecord->targetSpecific.pixelFormatObject, ((windowRecord->slaveWindow) ? windowRecord->slaveWindow->targetSpecific.contextObject : NULL),
                           &(windowRecord->targetSpecific.contextObject));
    if (error) {
        printf("\nPTB-ERROR[ContextCreation failed: %s]: Could not create master OpenGL context for new onscreen window. Insufficient video memory?\n\n", CGLErrorString(error));
        return(FALSE);
    }

    // Enable the OpenGL rendering context associated with our window:
    error=CGLSetCurrentContext(windowRecord->targetSpecific.contextObject);
    if (error) {
        printf("\nPTB-ERROR[SetCurrentContext failed: %s]: Insufficient video memory\n\n", CGLErrorString(error));
        return(FALSE);
    }

    // NULL-out the AGL context field, just for safety...
    windowRecord->targetSpecific.deviceContext = NULL;

    // Ok, the master OpenGL rendering context for this new onscreen window is up and running.
    // Auto-detect and bind all available OpenGL extensions via GLEW:
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

    // External 3D graphics support enabled?
    if (PsychPrefStateGet_3DGfx()) {
        // Yes. We need to create an extra OpenGL rendering context for the external
        // OpenGL code to provide optimal state-isolation. The context shares all
        // heavyweight ressources likes textures, FBOs, VBOs, PBOs, shader, display lists and
        // starts off as an identical copy of PTB's context as of here.
        error=CGLCreateContext(windowRecord->targetSpecific.pixelFormatObject, windowRecord->targetSpecific.contextObject, &(windowRecord->targetSpecific.glusercontextObject));
        if (error) {
            printf("\nPTB-ERROR[UserContextCreation failed: %s]: Creating a private OpenGL context for userspace OpenGL failed.\n\n", CGLErrorString(error));
            return(FALSE);
        }
    }

    // Create glswapcontextObject - An OpenGL context for exclusive use by parallel background threads,
    // e.g., our thread for async flip operations and self-made frame-sequential stereo:
    error=CGLCreateContext(windowRecord->targetSpecific.pixelFormatObject, windowRecord->targetSpecific.contextObject, &(windowRecord->targetSpecific.glswapcontextObject));
    if (error) {
        printf("\nPTB-ERROR[SwapContextCreation failed: %s]: Creating a private OpenGL context for async-bufferswaps failed.\n\n", CGLErrorString(error));
        CGLSetCurrentContext(NULL);
        return(FALSE);
    }

    // Store Cocoa onscreen window handle:
    windowRecord->targetSpecific.windowHandle = cocoaWindow;

    // Objective-C setup path, using Cocoa + NSOpenGLContext wrapped around already
    // existing and setup CGLContext:
    if (PsychCocoaSetupAndAssignOpenGLContextsFromCGLContexts(cocoaWindow, windowRecord)) {
        printf("\nPTB-ERROR[Cocoa OpenGL setup failed]: Setup failed for unknown reasons.\n\n");
        PsychCocoaDisposeWindow(windowRecord);
        return(FALSE);
    }

    // Check for output display rotation enabled. Will impair timing/timestamping because
    // it uses the desktop compositor for a rotated copy blit, instead of via rotated crtc
    // scanout, as most crtc's don't support this in hardware:
    if ((((int) CGDisplayRotation(cgDisplayID)) != 0) && (PsychPrefStateGet_Verbosity() > 1)) {
        printf("PTB-WARNING: Your onscreen windows output display has rotation enabled. It is not displaying in upright orientation.\n");
        printf("PTB-WARNING: On most graphics cards this will cause unreliable stimulus presentation timing and timestamping.\n");
        printf("PTB-WARNING: If you want non-upright stimulus presentation, look at 'help PsychImaging' on how to achieve this in\n");
        printf("PTB-WARNING: a way that doesn't impair timing. The subfunctions 'FlipHorizontal' and 'FlipVertical' are what you probably need.\n");
    }

    // First reference to this screen by a window?
    if (screenRefCount[screenSettings->screenNumber] == 0) {
        // High precision timestamping enabled? If so, we need to setup the fallback
        // timestamping methods in case beamposition timestamping doesn't work:
        if (PsychPrefStateGet_VBLTimestampingMode() > 0) {
            // Use of CoreVideo is needed on 10.7 and later due to brokeness of the old method (thanks Apple!):
            if (PsychPrefStateGet_Verbosity() > 2) {
                printf("PTB-INFO: Will use fragile CoreVideo timestamping as fallback if beamposition timestamping doesn't work.\n");
                // Recommend use of kernel driver if it isn't installed already for all but Intel GPU's:
                if (!PsychOSIsKernelDriverAvailable(screenSettings->screenNumber) && !strstr((char*) glGetString(GL_VENDOR), "Intel")) {
                    printf("PTB-INFO: Installation of the PsychtoolboxKernelDriver is strongly recommended if you care about precise visual\n");
                    printf("PTB-INFO: onset timestamping or timing. See 'help PsychtoolboxKernelDriver' for installation instructions.\n");
                }
            }

            if (NULL == cvDisplayLink[screenSettings->screenNumber]) {
                // CoreVideo timestamping:
                //
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
        }
        else {
            // VBLtimestampingmode 0 or -1 -- No CoreVideo fallback for timestamping if beamposition timestamping is unavailable:
            // This is the new default as of Psychtoolbox 3.0.12 to avoid the buggy, crashy, unreliably CoreVideo fallback.

            // Recommend use of kernel driver if it isn't installed already for all but Intel GPU's:
            if (!PsychOSIsKernelDriverAvailable(screenSettings->screenNumber) && !strstr((char*) glGetString(GL_VENDOR), "Intel")) {
                printf("PTB-INFO: Installation of the PsychtoolboxKernelDriver is strongly recommended if you care about precise visual\n");
                printf("PTB-INFO: onset timestamping or timing. See 'help PsychtoolboxKernelDriver' for installation instructions.\n");
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
    double cvTime;

    // Should we use CoreVideo display link timestamping?
    if (cvDisplayLink[screenid]) {
        // Yes: Retrieve data from our shared data structure:

        PsychLockMutex(&(cvDisplayLinkData[screenid].mutex));
        *vblCount = cvDisplayLinkData[screenid].vblCount;
        cvTime = cvDisplayLinkData[screenid].vblTimestamp;
        PsychUnlockMutex(&(cvDisplayLinkData[screenid].mutex));

        return(cvTime);
    }
    else {
        // Unsupported:
        *vblCount = 0;
        return(-1);
    }
}

void PsychOSCloseWindow(PsychWindowRecordType *windowRecord)
{
    CGDirectDisplayID   cgDisplayID;

    // Disable rendering context:
    CGLSetCurrentContext(NULL);

    // Destroy pixelformat object:
    CGLDestroyPixelFormat(windowRecord->targetSpecific.pixelFormatObject);

    // Destroy rendering context:
    CGLReleaseContext(windowRecord->targetSpecific.contextObject);
    if (windowRecord->targetSpecific.glusercontextObject) CGLReleaseContext(windowRecord->targetSpecific.glusercontextObject);
    if (windowRecord->targetSpecific.glswapcontextObject) CGLReleaseContext(windowRecord->targetSpecific.glswapcontextObject);

    // Last reference to this screen? In that case we have to shutdown the fallback
    // vbl timestamping and vblank counting facilities for this screen:
    if (screenRefCount[windowRecord->screenNumber] == 1) {
        // Last one on this screen will be gone in a second.
        // Shutdown and release CVDisplayLink for this windows screen, if any:
        if (cvDisplayLink[windowRecord->screenNumber]) {
            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Releasing CVDisplayLink for screen %i.\n", windowRecord->screenNumber);

            if (CVDisplayLinkIsRunning(cvDisplayLink[windowRecord->screenNumber])) CVDisplayLinkStop(cvDisplayLink[windowRecord->screenNumber]);
            PsychYieldIntervalSeconds(0.1);
            CVDisplayLinkRelease(cvDisplayLink[windowRecord->screenNumber]);
            cvDisplayLink[windowRecord->screenNumber] = NULL;
            PsychYieldIntervalSeconds(0.1);

            // Teardown shared data structure and mutex:
            PsychDestroyMutex(&(cvDisplayLinkData[windowRecord->screenNumber].mutex));
        }
    }

    // Release reference of this window to its screen:
    screenRefCount[windowRecord->screenNumber]--;

    // Destroy Cocoa onscreen window, if any:
    if (windowRecord->targetSpecific.windowHandle) PsychCocoaDisposeWindow(windowRecord);
    windowRecord->targetSpecific.windowHandle = NULL;

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
	
	if (copyfromPTBContext && (PsychPrefStateGet_Verbosity() > 1)) {
        // This was deprecated by the iPhone company as of OSX 10.8 - who are we to question their wisdom? Luckily seldomly needed in practice:
        printf("PTB-WARNING: Screen('BeginOpenGL', windowPtr, 2) called to synchronize userspace context state with Screen state. This is unsupported on OSX. Code may misbehave!\n");
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
        
        // Detach totally:
        CGLSetCurrentContext(NULL);        
	}
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
    if (error || (myinterval != (GLint) swapInterval)) {
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
		PsychCocoaGetWindowBounds(windowRecord->targetSpecific.windowHandle, windowRecord->globalrect, windowRecord->rect);
		PsychSetupClientRect(windowRecord);
		PsychSetupView(windowRecord, FALSE);
	}
}
