/*
	SCREENOpenWindow.c		

	AUTHORS:

	Allen.Ingling@nyu.edu		awi 
	mario.kleiner@tuebingen.mpg.de  mk
  

	PLATFORMS:	All

	HISTORY:

		12/18/01	awi		Created.  Copied the Synopsis string from old version of psychtoolbox. 
		10/18/02	awi		Added defaults to allow for optional arguments.
		12/05/02	awi		Started over again for OS X without SDL.
		10/12/04	awi		In useString: changed "SCREEN" to "Screen", and moved commas to inside [].
		2/15/05		awi		Commented out glEnable(GL_BLEND) and mode settings.  
		04/03/05	mk		Added support for selecting binocular stereo output via native OpenGL.
 		11/14/06	mk		New onscreen windows blank to their background color after successfull init.
							Support for specification of pixelSize's for 10-10-10-2, 16-16-16-16 and
							32-32-32-32 framebuffers on supported hardware.
	TO DO:

*/

#include "Screen.h"

// Pointer to master onscreen window during setup phase of stereomode 10 (Dual-window stereo):
static PsychWindowRecordType* sharedContextWindow = NULL;

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] =  "[windowPtr,rect]=Screen('OpenWindow',windowPtrOrScreenNumber [,color] [,rect][,pixelSize][,numberOfBuffers][,stereomode][,multisample][,imagingmode][,specialFlags][,clientRect]);";
//                                                               1                         2        3      4           5                 6            7             8             9              10
static char synopsisString[] =
    "Open an onscreen window. Specify a screen by a windowPtr or a screenNumber (0 is "
    "the main screen, with menu bar). \"color\" is the clut index (scalar or [r g b] "
    "triplet) that you want to poke into each pixel; default is white. If supplied, "
    "\"rect\" must contain at least one pixel. If a windowPtr is supplied then \"rect\" "
    "is in the window's coordinates (origin at upper left), and defaults to the whole "
    "window. If a screenNumber is supplied then \"rect\" is in screen coordinates "
    "(origin at upper left), and defaults to the whole screen. (In all cases, "
    "subsequent references to this new window will use its coordinates: origin at its "
    "upper left.). Please note that while providing a \"rect\" parameter to open a normal "
    "window instead of a fullscreen window is convenient for debugging, but drawing performance, "
    "stimulus onset timing and onset timestamping may be impaired, so be careful.\n"
    "\"pixelSize\" sets the depth (in bits) of each pixel; default is to leave depth unchanged. "
    "\"numberOfBuffers\" is the number of buffers to use. Setting anything else than 2 will be "
    "useful for development/debugging of PTB itself but will mess up any real experiment. "
    "\"stereomode\" Type of stereo display algorithm to use: 0 (default) means: Monoscopic viewing. "
    "1 means: Stereo output via OpenGL on any stereo hardware that is supported by MacOS-X, e.g., the "
    "shutter glasses from CrystalView. 2 means: Left view compressed into top half, right view into bottom half. "
    "3 means left view compressed into bottom half, right view compressed into top half. 4 and 5 allow split "
    "screen display where left view is shown in left half, right view is shown in right half or the display. "
    "A value of 5 does the opposite (cross-fusion). Values of 6,7,8 and 9 enable Anaglyph stereo rendering "
    "of types left=Red, right=Green, vice versa and left=Red, right=Blue and vice versa. A value of 10 "
    "enables multi-window stereo: Open one window for left eye view, one for right eye view, treat both "
    "of them as one single stereo window. A value of 11 enables our own frame-sequential stereo mode. "
    "See StereoDemo.m for examples of usage of the different stereo "
    "modes. See ImagingStereoDemo.m for more advanced usage on modern hardware.\n"
    "\"multisample\" This parameter, if provided and set to a value greater than zero, enables automatic "
    "hardware anti-aliasing of the display: For each pixel, 'multisample' color samples are computed and "
    "combined into a single output pixel color. Higher numbers provide better quality but consume more "
    "video memory and lead to a reduction in framerate due to the higher computational demand. The maximum "
    "number of samples is hardware dependent. Psychtoolbox will silently clamp the number to the maximum "
    "supported by your hardware if you ask for too much. On very old hardware, the value will be ignored. "
    "Read 'help AntiAliasing' for more in-depth information about multi-sampling. "
    "\"imagingmode\" This optional parameter enables PTB's internal image processing pipeline. The pipeline is "
    "off by default. Read 'help PsychImaging' for information about typical use and benefits of this feature.\n"
    "\"specialFlags\" This optional parameter enables some special window behaviours if the sum of certain "
    "flags is passed. A currently supported flag is the symbolic constant kPsychGUIWindow. It enables windows "
    "to behave more like regular GUI windows on your system. See 'help kPsychGUIWindow' for more info.\n"
    "\"clientRect\" This optional parameter allows to define a size of the onscreen windows drawing area "
    "that is different from the actual size of the windows framebuffer. If set, then the imaging pipeline "
    "is started and a virtual framebuffer of the size of \"clientRect\" is created. Your code will draw "
    "into that framebuffer. At display time, the content of this virtual framebuffer will get scaled to "
    "the size of the true onscreen window, a process known as panel-scaling or panel-fitting. This allows "
    "to decouple the size of a stimulus as drawn by your code from the actual resolution of the display "
    "device. The feature is mostly useful if you need to run the same presentation code on different setups "
    "with different native resolutions. See the 'help PsychImaging' section about 'UsePanelFitter' for more info.\n"
    "\n"
    "Opening or closing a window takes about one to three seconds, depending on type of connected display. "
    "COMPATIBILITY TO OS-9 PTB: If you absolutely need to run old code for the old MacOS-9 or Windows "
    "Psychtoolbox, you can switch into a compatibility mode by adding the command "
    "Screen('Preference', 'EmulateOldPTB', 1) at the very top of your script. This will restore "
    "Offscreen windows and WaitBlanking functionality, but at the same time disable most of the new "
    "features of the OpenGL Psychtoolbox. Please do not write new experiment code in the old style! "
    "Emulation mode is pretty new and may contain significant bugs, so use with great caution!";  

static char seeAlsoString[] = "OpenOffscreenWindow, SelectStereoDrawBuffer, PanelFitter, Close, CloseAll";

PsychError SCREENOpenWindow(void)
{
    int						screenNumber, numWindowBuffers, stereomode, multiSample, imagingmode, specialflags;
    PsychRectType 			rect, screenrect, clientRect;
    PsychColorType			color;
    psych_bool				isArgThere, didWindowOpen, useAGL;
    PsychScreenSettingsType	screenSettings;
    PsychWindowRecordType	*windowRecord;
    PsychDepthType			specifiedDepth, possibleDepths, currentDepth, useDepth;
	int                     dummy1;
	double                  dummy2, dummy3, dummy4;

	psych_bool EmulateOldPTB = PsychPrefStateGet_EmulateOldPTB();

    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //cap the number of inputs
    PsychErrorExit(PsychCapNumInputArgs(10));   //The maximum number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(2));  //The maximum number of outputs

    //get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.  
    PsychCopyInScreenNumberArg(kPsychUseDefaultArgPosition, TRUE, &screenNumber);
    if(screenNumber==-1)
        PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen."); 

    /*
      The depth checking is ugly because of this stupid depth structure stuff.  
      Instead get a descriptor of the current video settings, change the depth field,
      and pass it to a validate function wich searches a list of valid video modes for the display.
      There seems to be no point in checking the depths alone because the legality of a particular
      depth depends on the other settings specified below.  Its probably best to wait until we have
      digested all settings and then test the full mode, declarin an invalid
      mode and not an invalid pixel size.  We could notice when the depth alone is specified 
      and in that case issue an invalid depth value.
     */  

    //find the PixelSize first because the color specifier depends on the screen depth.  
    PsychInitDepthStruct(&currentDepth);  //get the current depth
    PsychGetScreenDepth(screenNumber, &currentDepth);
    PsychInitDepthStruct(&possibleDepths); //get the possible depths
    PsychGetScreenDepths(screenNumber, &possibleDepths);

    #if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_WINDOWS
       // MK Experimental Hack: Add the special depth values 64 and 128 to the depth struct. This should 
       // allows for 16 bpc, 32 bpc floating point color buffers on the latest ATI and NVidia hardware.
	   // "Should" means: It doesn't really work with any current driver, but we leave the testcode in
	   // in the hope for future OS and driver releases ;-)
       // Unfortunately at this point of the init sequence, we are not able
       // to check if these formats are supported by the hardware. Ugly ugly ugly...
       PsychAddValueToDepthStruct(64, &possibleDepths);
       PsychAddValueToDepthStruct(128, &possibleDepths);
    #endif

//    #if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
		// On MacOS/X and Linux with ATI Radeon X1000/HD2000/HD3000 hardware and the special
		// kernel support driver installed, we should be able to configure the hardwares
		// framebuffers into ABGR2101010 mode, ie. 2 bits alpha, 10 bpc for red, green, blue.
		// This needs support from the imaging pipeline, or manually converted stimuli, as
		// the GPU doesn't format pixel data properly, only the CRTC scans out in that format.
		// Anyway, allow this setting on OS/X and Linux:
		
		// Update: Some FireGL cards (2008 and later) claim to support this on MS-Windows. Enable
		// this option on Windows as well, so it is at least testable:
		PsychAddValueToDepthStruct(30, &possibleDepths);
//    #endif


    PsychInitDepthStruct(&specifiedDepth); //get the requested depth and validate it.  
    isArgThere = PsychCopyInSingleDepthArg(4, FALSE, &specifiedDepth);

    PsychInitDepthStruct(&useDepth);
    if(isArgThere){ //if the argument is there check that the screen supports it...
        if(!PsychIsMemberDepthStruct(&specifiedDepth, &possibleDepths))
            PsychErrorExit(PsychError_invalidDepthArg);
        else
            PsychCopyDepthStruct(&useDepth, &specifiedDepth);
    }else //otherwise use the default
        PsychCopyDepthStruct(&useDepth, &currentDepth);

    // Initialize the rect argument to the screen rectangle:
    PsychGetGlobalScreenRect(screenNumber, rect); 	//get the rect describing the screen bounds.  This is the default Rect.  

    // Override it with a user supplied rect, if one was supplied:
    isArgThere=PsychCopyInRectArg(kPsychUseDefaultArgPosition, FALSE, rect );
    if (IsPsychRectEmpty(rect)) PsychErrorExitMsg(PsychError_user, "OpenWindow called with invalid (empty) rect argument.");

	if (PSYCH_SYSTEM == PSYCH_OSX) {
		// OS/X system: Need to decide if we use a Carbon window + AGL
		// or a fullscreen context with CGL:
		
		// Default to AGL, switch to CGL if below constraints are met:
		useAGL = TRUE;
	
		// Window rect provided which has a different size than screen?

		// We do not use windowed mode if the provided window rectangle either
		// matches the target screens rectangle (and therefore its exact size)
		// or its screens global rectangle.
		PsychGetScreenRect(screenNumber, screenrect);
		if (PsychMatchRect(screenrect, rect)) useAGL=FALSE;
		PsychGetGlobalScreenRect(screenNumber, screenrect);
		if (PsychMatchRect(screenrect, rect)) useAGL=FALSE;

		// Override for use with Quartz compositor: Must not capture/release screen, therefore
        // set useAGL = true to prevent screen capture/release:
		if ((PsychPrefStateGet_ConserveVRAM() & kPsychUseAGLCompositorForFullscreenWindows) ||
            (PsychPrefStateGet_WindowShieldingLevel() < 2000) ||
            (PsychPrefStateGet_ConserveVRAM() & kPsychUseAGLForFullscreenWindows)) {
            useAGL = TRUE;
        }
	}
	else {
		// Non OS/X system: Do not use AGL ;-)
		useAGL = FALSE;
	}
	
    //find the number of specified buffers. 
    numWindowBuffers=2;	
    PsychCopyInIntegerArg(5,FALSE,&numWindowBuffers);
    if(numWindowBuffers < 1 || numWindowBuffers > kPsychMaxNumberWindowBuffers) PsychErrorExit(PsychError_invalidNumberBuffersArg);

    stereomode=0;
    PsychCopyInIntegerArg(6,FALSE,&stereomode);
    if(stereomode < 0 || stereomode > 11) PsychErrorExitMsg(PsychError_user, "Invalid stereomode provided (Valid between 0 and 11).");
	if (stereomode!=0 && EmulateOldPTB) PsychErrorExitMsg(PsychError_user, "Sorry, stereo display functions are not supported in OS-9 PTB emulation mode.");

    multiSample=0;
    PsychCopyInIntegerArg(7,FALSE,&multiSample);
    if(multiSample < 0) PsychErrorExitMsg(PsychError_user, "Invalid multisample value provided (Valid are positive numbers >= 0).");
	if (multiSample!=0 && EmulateOldPTB) PsychErrorExitMsg(PsychError_user, "Sorry, anti-aliasing functions are not supported in OS-9 PTB emulation mode.");

	imagingmode=0;
    PsychCopyInIntegerArg(8,FALSE,&imagingmode);
    if(imagingmode < 0) PsychErrorExitMsg(PsychError_user, "Invalid imaging mode provided (See 'help PsychImagingMode' for usage info).");
	if (imagingmode!=0 && EmulateOldPTB) PsychErrorExitMsg(PsychError_user, "Sorry, imaging pipeline functions are not supported in OS-9 PTB emulation mode.");
	
	specialflags=0;
    PsychCopyInIntegerArg(9,FALSE,&specialflags);
    if (specialflags < 0 || (specialflags > 0 && specialflags!=kPsychGUIWindow)) PsychErrorExitMsg(PsychError_user, "Invalid 'specialflags' provided.");

    // Optional clientRect defined? If so, we need to enable our internal panel scaler and
    // the imaging pipeline to actually use the scaler:
    if (PsychCopyInRectArg(10, FALSE, clientRect)) {
        // clientRect given. The panelscaler integrated into the imaging pipeline will
        // scale all content from the size of the drawBufferFBO (our virtual framebuffer),
        // which is the size of the clientRect, to the true size of the onscreen windows
        // system framebuffer - appropriately tweaked for special display modes of course.
        
        // Validate clientRect:
        if (IsPsychRectEmpty(clientRect)) PsychErrorExitMsg(PsychError_user, "OpenWindow called with invalid (empty) 'clientRect' argument.");
        
        // Set special imagingmode flags to signal need for full imaging pipeline
        // and for the panel scaler. Used in PsychInitializeImagingPipeline() and
        // to make sure PsychOpenOnscreenWindow() gets called with a multisample value
        // of zero, so the system backbuffer isn't multisampled -- crucial for us!
        // This will also turn PsychSetupClientRect() into a no-op:
        imagingmode |= kPsychNeedFastBackingStore;
        imagingmode |= kPsychNeedGPUPanelFitter;        
    }
    
    
	// We require use of the imaging pipeline if stereomode for dualwindow display is requested.
	// This makes heavy use of FBO's and blit operations, so imaging pipeline is needed.
	if ((stereomode==kPsychDualWindowStereo) || (imagingmode & kPsychNeedDualWindowOutput)) {
		// Dual window stereo requested, but imaging pipeline not enabled. Enable it:
		imagingmode |= kPsychNeedFastBackingStore;
		if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Trying to enable imaging pipeline for dual-window stereo display mode or dual-window output mode...\n");
	}

    // We also require imaging pipeline if homegrown frameseq. stereo is requested. Need to do this here,
    // so the call below to PsychOpenOnscreenWindow() knows already about use of imaging pipe and can
    // do the right thing wrt. to multisampling. Most of the setup code for kPsychFrameSequentialStereo
    // follows after opening the window. Rationale: multisampling must be off on the system framebuffer,
    // otherwise we will get into invalid operating conditions for multisample resolve ops from within
    // imaging pipeline.
    if (stereomode == kPsychFrameSequentialStereo) imagingmode |= kPsychNeedFastBackingStore;
	
    //set the video mode to change the pixel size.  TO DO: Set the rect and the default color  
    PsychGetScreenSettings(screenNumber, &screenSettings);    
    PsychInitDepthStruct(&(screenSettings.depth));
    PsychCopyDepthStruct(&(screenSettings.depth), &useDepth);

    // Here is where all the work goes on:

    // If the screen is not already captured then to that:
    if(!PsychIsScreenCaptured(screenNumber) && !useAGL) {
        PsychCaptureScreen(screenNumber);
    }

#if PSYCH_SYSTEM == PSYCH_WINDOWS
    // On M$-Windows we currently only support - and therefore require >= 30 bpp color depth.
    if (PsychGetScreenDepthValue(screenNumber) < 30) {
		// Display running at less than 30 bpp. OpenWindow will fail on M$-Windows anyway, so let's abort
		// now.

		// Output warning text:
        printf("PTB-ERROR: Your display screen %i is not running at the required color depth of at least 30 bit.\n", screenNumber);
        printf("PTB-ERROR: The current setting is %i bit color depth..\n", PsychGetScreenDepthValue(screenNumber));
        printf("PTB-ERROR: This will not work on Microsoft Windows operating systems.\n");
        printf("PTB-ERROR: Please use the 'Display settings' control panel of Windows to change the color depth to\n");
        printf("PTB-ERROR: 32 bits per pixel ('True color' or 'Highest' setting) and then retry. It may be neccessary\n");
        printf("PTB-ERROR: to restart Matlab after applying the change...\n");
        fflush(NULL);

		// Release the captured screen:
		PsychRestoreScreenSettings(screenNumber);
		PsychReleaseScreen(screenNumber);

        // Reset master assignment to prepare possible further dual-window config operations:
		sharedContextWindow = NULL;

		// Abort with Matlab error:
		PsychErrorExitMsg(PsychError_user, "Insufficient color depth setting for display device (smaller than 30 bpp).");
    }

#endif

    //if (PSYCH_DEBUG == PSYCH_ON) printf("Entering PsychOpenOnscreenWindow\n");
    PsychCopyDepthStruct(&(screenSettings.depth), &useDepth);
	
	// Create the onscreen window and perform initialization of everything except
	// imaging pipeline and a few other special quirks. If sharedContextWindow is non-NULL,
	// the new window will share its OpenGL context ressources with sharedContextWindow.
	// This is typically used for dual-window stereo mode. Btw. If imaging pipeline is really
	// active, we force multiSample to zero: This way the system backbuffer / pixelformat
	// is enabled without multisampling support, as we do all the multisampling stuff ourselves
	// within the imaging pipeline with multisampled drawbuffer FBO's...
    didWindowOpen=PsychOpenOnscreenWindow(&screenSettings, &windowRecord, numWindowBuffers, stereomode, rect, ((imagingmode==0 || imagingmode==kPsychNeedFastOffscreenWindows) ? multiSample : 0), sharedContextWindow, specialflags);
    if (!didWindowOpen) {
        if (!useAGL) {
			PsychRestoreScreenSettings(screenNumber);
			PsychReleaseScreen(screenNumber);
		}

		// Reset master assignment to prepare possible further dual-window config operations:
		sharedContextWindow = NULL;

        // We use this dirty hack to exit with an error, but without printing
        // an error message. The specific error message has been printed in
        // PsychOpenOnscreenWindow() already..
        PsychErrMsgTxt("");
    }

    // Sufficient display depth for full alpha-blending and such?
    if (PsychGetScreenDepthValue(screenNumber) < 24) {
        // Nope. Output a little warning.
        printf("PTB-WARNING: Your display screen %i is not running at 24 bit color depth or higher.\n", screenNumber);
        printf("PTB-WARNING: The current setting is %i bit color depth..\n", PsychGetScreenDepthValue(screenNumber));
        printf("PTB-WARNING: This could cause failure to work correctly or visual artifacts in stimuli\n");
        printf("PTB-WARNING: that involve Alpha-Blending. It can also cause drastically reduced color resolution\n");
        printf("PTB-WARNING: for your stimuli! Please try to switch your display to 'True Color' (Windows)\n");
        printf("PTB-WARNING: our 'Millions of Colors' (MacOS-X) to get rid of this warning and the visual artifacts.\n");
        fflush(NULL);
    }
    
	// Define clear color: This depends on the color range of our onscreen window...
    isArgThere=PsychCopyInColorArg(kPsychUseDefaultArgPosition, FALSE, &color); //get from user
    if(!isArgThere) PsychLoadColorStruct(&color, kPsychIndexColor, PsychGetWhiteValueFromWindow(windowRecord)); //or use the default
    PsychCoerceColorMode(&color);

    // The imaging pipeline and graphics drivers had over 5 years of time to mature. As of 2012, imaging pipeline based
    // support for fast offscreen windows and for stereoscopic display modes is far superior in performance,
    // robustness, flexibility and convenience to the legacy method which was used in ptb by default so far.
    // Now it is 2012+ and we switch the defaults: If the GPU+driver combo supports it, and usercode doesn't
    // actively opt-out of it, we auto-enable use of FBO backed fast offscreen windows. We don't auto-enable
    // the full pipeline for stereoscopic display modes, but we print some recommendations to the user to
    // consider enabling the full pipeline for stereo display:
    if ((windowRecord->gfxcaps & kPsychGfxCapFBO) && !(PsychPrefStateGet_ConserveVRAM() & kPsychDontAutoEnableImagingPipeline)) {
        // Support for basic use of the PTB imaging pipeline and/or for fast offscreen windows
        // is available - a GPU + driver combo with support for OpenGL framebuffer objects with
        // at least RGBA8 format and rectangle rendertargets.
        // Usercode doesn't disallow automatic use of imaging pipeline or fast offscreen windows,
        // ie. it didn't set the kPsychDontAutoEnableImagingPipeline conserveVRAM flag.
        // Good!
        
        // We will therefore auto-enable use of fast offscreen windows:
        imagingmode |= kPsychNeedFastOffscreenWindows;

        // Is a stereomode requested which would benefit from enabling the full imaging pipeline?
        if (stereomode > 0) {
            if (((stereomode == kPsychOpenGLStereo) && !(windowRecord->gfxcaps & kPsychGfxCapNativeStereo)) || (stereomode == kPsychFrameSequentialStereo)) {
                // Native OpenGL quad-buffered frame-sequential stereo requested, but unsupported by gpu & driver.
                // Or use of our own method requested. We have FBO and framebuffer blit support, so we can roll our
                // own framesequential stereo by use of the imaging pipeline.
                
                // Sanity check: If multisampling is enabled and imaging pipeline isn't yet enabled,
                // enabling it now will cause trouble: With pipeline enabled, we need a system framebuffer
                // without multisampling, but the system framebuffer is multisampled, as the setup code
                // has already executed without knowing about this constraint. We can't go on with multisampling
                // at this point or malfunctions will happen. So if this check triggers, take the lesser of two
                // evils and disable multisampling and tell user how to resolve the problem properly:
                if ((imagingmode == kPsychNeedFastOffscreenWindows) && (multiSample > 0)) {
                    // Troublesome. Disable our own multisampling, as it clashes with the fact that
                    // the onscreen windows system framebuffer already is multisampled.
                    multiSample = 0;
                    windowRecord->multiSample = 0;
                    if (PsychPrefStateGet_Verbosity() > 1) {
                        printf("\nPTB-WARNING: You are trying to use frame-sequential stereo with multisample anti-aliasing, but you don't use\n");
                        printf("PTB-WARNING: PsychImaging('OpenWindow', ...) to do this. This mode is unsupported by your system,\n");
                        printf("PTB-WARNING: so i'm trying now to enable some workaround, which however is incompatible with multisample\n");
                        printf("PTB-WARNING: anti-aliasing. Will disable anti-aliasing now. If you don't like this, please change your code\n");
                        printf("PTB-WARNING: to use PsychImaging('OpenWindow',...) instead of Screen('OpenWindow',...); and the problem\n");
                        printf("PTB-WARNING: will be automatically resolved, ie., you can have frame-sequential stereo and anti-aliasing!\n\n");
                    }
                }
                
                // Enable basic imaging pipeline:
                imagingmode |= kPsychNeedFastBackingStore;
                
                // Override stereomode to our own homegrown implementation:
                stereomode = kPsychFrameSequentialStereo;
                windowRecord->stereomode = stereomode;
                
                if (PsychPrefStateGet_Verbosity() > 2) {
                    printf("\n");
                    printf("PTB-INFO: Your script requests use of frame-sequential stereo, but your graphics card\n");
                    printf("PTB-INFO: and driver doesn't support this. I will now fully enable the imaging pipeline\n");
                    printf("PTB-INFO: and use my own home-grown frame-sequential stereo implementation. Note that this\n");
                    printf("PTB-INFO: may not be as robust and high-performance as using a graphics card with native\n");
                    printf("PTB-INFO: frame-sequential stereo support. But let's see what i can do for you...\n\n");
                }
            }
            else {
                // Yes: Provide the user with recommendations to enable the pipeline.
                if (!(imagingmode & kPsychNeedFastBackingStore) && (PsychPrefStateGet_Verbosity() > 2)) {
                    printf("\n");
                    printf("PTB-INFO: Your script requests use of a stereoscopic display mode (stereomode = %i).\n", stereomode);
                    printf("PTB-INFO: Stereoscopic stimulus display is usually more flexible, convenient and robust if\n");
                    printf("PTB-INFO: the Psychtoolbox imaging pipeline is enabled. Your graphics card is capable\n");
                    printf("PTB-INFO: of using the pipeline but your script doesn't enable use of the pipeline.\n");
                    printf("PTB-INFO: I recommend you enable use of the pipeline for enhanced stereo stimulus display.\n");
                    printf("PTB-INFO: Have a look at the demoscript ImagingStereoDemo.m on how to do this.\n\n");
                }
            }
        }
    }

	// Query if OpenGL stereo is natively supported or if our own emulation mode will work:
	if ((((stereomode == kPsychOpenGLStereo) && !(windowRecord->gfxcaps & kPsychGfxCapNativeStereo)) || (stereomode == kPsychFrameSequentialStereo)) &&
	    (!(imagingmode & kPsychNeedFastBackingStore) || (windowRecord->stereomode != kPsychFrameSequentialStereo) || !(windowRecord->gfxcaps & kPsychGfxCapFBO))) {
		// OpenGL native stereo was requested, but is obviously not supported and we can't roll our own implementation either :-(
		printf("\nPTB-ERROR: Asked for OpenGL native stereo (frame-sequential mode) but this doesn't seem to be supported by your graphics hardware or driver.\n");
		printf("PTB-ERROR: Unfortunately using my own implementation via imaging pipeline did not work either, due to lack of hardware support, or because\n");
		printf("PTB-ERROR: did not allow me to auto-enable the pipeline and use this method. This means game over!\n");
		if (PSYCH_SYSTEM == PSYCH_OSX) {
			printf("PTB-ERROR: Frame-sequential stereo should be supported on all recent ATI/AMD and NVidia cards on OS/X, except for the Intel onboard chips,\n");
			printf("PTB-ERROR: at least in fullscreen mode with OS/X 10.5, and also mostly on OS/X 10.4. If it doesn't work, check for OS updates etc.\n\n");
		}
		else {
			printf("PTB-ERROR: Frame-sequential native stereo on Windows or Linux is usually only supported with the professional line of graphics cards\n");
			printf("PTB-ERROR: from NVidia and ATI/AMD, e.g., NVidia Quadro series or ATI Fire series. If you happen to have such a card, check\n");
			printf("PTB-ERROR: your driver settings and/or update your graphics driver.\n\n");
		}
		PsychErrMsgTxt("Frame-Sequential stereo display mode requested, but unsupported. Emulation unsupported as well. Game over!");
	}

	// Special setup code for dual window stereomode or output mode:
	if (stereomode == kPsychDualWindowStereo || (imagingmode & kPsychNeedDualWindowOutput)) {
		if (sharedContextWindow) {
			// This is creation & setup of the slave onscreen window, ie. the one
			// representing the right-eye or channel 1 view. This window doesn't do much. It
			// is not used or referenced in the users experiment script. It receives
			// its final image content during Screen('Flip') operation of the master
			// onscreen window, then gets flipped in sync with the master window.
			
			// Ok, we already have the slave window open and it shares its OpenGL context
			// with the master window. Reset its internal reference to the master:
			windowRecord->slaveWindow = NULL;
			
			// Reset imagingmode for this window prior to imaging pipeline setup. This
			// window is totally passive so it doesn't need the imaging pipeline.
			imagingmode = 0;

			// Assign this window to the master window as a slave:
			sharedContextWindow->slaveWindow = windowRecord;
			
			// Try to optionally enable framelock / swaplock extensions for the window-pair
			// if this is supported by the given system configuration. If supported, this
			// should guarantee perfect synchronization of bufferswaps across the window-pair:
			PsychOSSetupFrameLock(sharedContextWindow, windowRecord);

			// Reset master assignment to prepare possible further dual-window config operations:
			sharedContextWindow = NULL;

			// Activate the IdentitiyBlitChain for the slave window and add a single identity blit
			// operation to it: This is needed in PsychPreFlipOperations() for final copy of stimulus
			// image into this slave window:
			PsychPipelineAddBuiltinFunctionToHook(windowRecord, "IdentityBlitChain", "Builtin:IdentityBlit", INT_MAX, "");
			PsychPipelineEnableHook(windowRecord, "IdentityBlitChain");

			if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Created master-slave window relationship for dual-window stereo/output display mode...\n");

			// Special config finished. The master-slave combo should work from now on...			
		}
		else {
			// This is initial setup & creation of the master onscreen window, ie. the one
			// representing the left-eye or channel 0 view and doing all the heavy work, acting as a
			// proxy for both windows.
			
			// Not much to do here. Just store its windowRecord as a reference for creation
			// of the slave window. We'll need it for that purpose...
			sharedContextWindow = windowRecord;
		}
	}

	// Set special half-width flag for window if we are either in a dual-display/dual-view stereo mode or if
	// if is requested as part of the imagingMode flag. This will cause PTB 2D drawing routines and window size
	// query routines etc. to return an effective window width or window rect only half the real width.
	if (windowRecord->stereomode==kPsychFreeFusionStereo || windowRecord->stereomode==kPsychFreeCrossFusionStereo || (imagingmode & kPsychHalfWidthWindow)) {
		windowRecord->specialflags = windowRecord->specialflags | kPsychHalfWidthWindow;
		imagingmode = imagingmode & (~kPsychHalfWidthWindow);
	}

    // Similar handling for twice-width windows: Used for certain packed-pixels (2 stimulus pixels in one fb pixel) formats:
	if (imagingmode & kPsychTwiceWidthWindow) {
		windowRecord->specialflags = windowRecord->specialflags | kPsychTwiceWidthWindow;
		imagingmode = imagingmode & (~kPsychTwiceWidthWindow);
	}

	// Similar handling for windows of half the real height, except that none of our built-in stereo modes requires these,
	// so this is only done on request from external code via the imagingmode flag kPsychHalfHeightWindow.
	// One use of this is when using interleaved line stereo mode (PsychImaging(...'InterleavedLineStereo')) where windows
	// only have a useable net height of half their physical height:
	if (imagingmode & kPsychHalfHeightWindow) {
		windowRecord->specialflags = windowRecord->specialflags | kPsychHalfHeightWindow;
		imagingmode = imagingmode & (~kPsychHalfHeightWindow);
	}

    // Optional clientRect defined? If so, we need to enable our internal panel scaler and
    // the imaging pipeline to actually use the scaler:
    // This is part II, after part I happened above, before opening the window. This
    // weirdness / redundancy is needed to resolve our chicken & egg problem with
    // multisampling...
    if (imagingmode & kPsychNeedGPUPanelFitter) {
        // clientRect given. The panelscaler integrated into the imaging pipeline will
        // scale all content from the size of the drawBufferFBO (our virtual framebuffer),
        // which is the size of the clientRect, to the true size of the onscreen windows
        // system framebuffer - appropriately tweaked for special display modes of course.
        
        // Set it as "official" window client rectangle, whose size is reported
        // by default by functions like Screen('Rect'), Screen('WindowSize') or the
        // returned winRect of Screen('OpenWindow'):
        PsychNormalizeRect(clientRect, windowRecord->clientrect);
        PsychCopyRect(clientRect, windowRecord->clientrect);

        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Trying to enable my builtin panel-fitter on user request.\n");
    }
    else {
        // No specific clientRect given - the default case.
        
        // Define windows clientrect. It is a copy of windows rect, but stretched or compressed
        // to twice or half the width or height of the windows rect, depending on the special size
        // flags. clientrect is used as reference for all size query functions Screen('Rect'), Screen('WindowSize')
        // and for all Screen 2D drawing functions:
        PsychSetupClientRect(windowRecord);
    }
    
	// Initialize internal image processing pipeline if requested:
	if (numWindowBuffers > 1) PsychInitializeImagingPipeline(windowRecord, imagingmode, multiSample);

	if (imagingmode & kPsychNeedGPUPanelFitter) {
        // Setup default panelfitter parameters: This is a scaled blit, which does not
        // preserve the aspect-ratio of the virtual framebuffer, unless by pure chance
        // the aspect ratios of source and target are already identical:
        windowRecord->panelFitterParams[0] = 0; // srcX0
        windowRecord->panelFitterParams[1] = 0; // srcY0
        windowRecord->panelFitterParams[2] = (int) PsychGetWidthFromRect(clientRect);  // srcX1
        windowRecord->panelFitterParams[3] = (int) PsychGetHeightFromRect(clientRect); // srcY1
        
        windowRecord->panelFitterParams[4] = 0; // dstX0
        windowRecord->panelFitterParams[5] = 0; // dstY0
        windowRecord->panelFitterParams[6] = (int) windowRecord->fboTable[windowRecord->inputBufferFBO[0]]->width;  // dstX1
        windowRecord->panelFitterParams[7] = (int) windowRecord->fboTable[windowRecord->inputBufferFBO[0]]->height; // dstY1
    }

	// On OS-X, if we are in quad-buffered frame sequential stereo mode, we automatically generate
	// blue-line-sync style sync lines for use with stereo shutter glasses. We don't do this
	// by default on Windows or Linux: These systems either don't have stereo capable hardware,
	// or they have some and its drivers already take care of sync signal generation.
	if (((PSYCH_SYSTEM == PSYCH_OSX) && (windowRecord->stereomode == kPsychOpenGLStereo)) || (windowRecord->stereomode == kPsychFrameSequentialStereo)) {
		if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Enabling internal blue line sync renderer for quad-buffered stereo...\n");
		PsychPipelineAddBuiltinFunctionToHook(windowRecord, "LeftFinalizerBlitChain", "Builtin:RenderStereoSyncLine", INT_MAX, "");
		PsychPipelineEnableHook(windowRecord, "LeftFinalizerBlitChain");		
		PsychPipelineAddBuiltinFunctionToHook(windowRecord, "RightFinalizerBlitChain", "Builtin:RenderStereoSyncLine", INT_MAX, "");
		PsychPipelineEnableHook(windowRecord, "RightFinalizerBlitChain");		
	}

	// Activate new onscreen window for userspace drawing: If imaging pipeline is active, this
	// will bind the correct rendertargets for the first time. We soft-reset first to get
	// into a defined state:
	PsychSetDrawingTarget((PsychWindowRecordType*) 0x1);
	PsychSetDrawingTarget(windowRecord);

    // Set the clear color and perform a backbuffer-clear:
    PsychConvertColorToDoubleVector(&color, windowRecord, windowRecord->clearColor);
    PsychGLClear(windowRecord);

    // Mark end of drawing op. This is needed for single buffered drawing:
    PsychFlushGL(windowRecord);

    // Make sure no OpenGL errors happened up to this point:
    PsychTestForGLErrors();

    // If we are in logo-startup mode (former blue-screen mode) and double-buffering
    // is enabled, then do an initial bufferswap & clear, so the display starts in
    // the user selected background color instead of staying at the blue screen or
    // logo display until the Matlab script first calls 'Flip'.
    if (((PsychPrefStateGet_VisualDebugLevel()>=4) || (windowRecord->stereomode > 0)) && numWindowBuffers>=2) {
      // Do three immediate bufferswaps by an internal call to Screen('Flip'). This will also
      // take care of clearing the backbuffer in preparation of first userspace drawing
      // commands and such. We need up-to 3 calls to clear triple-buffered setups from framebuffer junk.
      PsychFlipWindowBuffers(windowRecord, 0, 0, 0, 0, &dummy1, &dummy2, &dummy3, &dummy4);
      PsychFlipWindowBuffers(windowRecord, 0, 0, 0, 0, &dummy1, &dummy2, &dummy3, &dummy4);
      PsychFlipWindowBuffers(windowRecord, 0, 0, 0, 0, &dummy1, &dummy2, &dummy3, &dummy4);
      // Display now shows background color, so user knows that PTB's 'OpenWindow'
      // procedure is successfully finished.
    }

    PsychTestForGLErrors();

    // Homegrown frame-sequential strereo mode on a EGL backed window active?
    if ((windowRecord->stereomode == kPsychFrameSequentialStereo) && (windowRecord->specialflags & kPsychIsEGLWindow)) {
        // Detach the OpenGL context from window surface. The following PsychSetDrawingTarget()
        // command will rebind the context as a first step, but it will not attach it to the
        // windowing system framebuffer surface (== the associated EGLSurface) anymore due to
        // the selected frame-sequential stereo mode on EGL. This will allow the background
        // frame-sequential stereo swapper-thread to bind its context to the surface and all
        // will be good. Rationale: With a EGL windowing system backend, only at most one context
        // is allowed to attach to a EGL surface (window framebuffer) at a given time. Because the
        // stereo thread needs to bind its context permanently to the surface to do its job, we
        // must make sure going forward that we'll never ever bind our contexts from the master-thread
        // to 'windowRecord's surface again, or bad things will happen:
        PsychSetDrawingTarget((PsychWindowRecordType*) 0x1);
        PsychOSUnsetGLContext(windowRecord);
        windowRecord->specialflags |= kPsychSurfacelessContexts;
        PsychSetDrawingTarget(windowRecord);
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Preventing master-thread contexts from future binding to this windows (%i) backing surface.\n", windowRecord->windowIndex);

        // TODO: Ideally we should check if KHR_surfaceless_context extension is supported,
        // because otherwise this won't work and we should reject use of multi-threaded ops
        // like frame-sequential stereo or async flips as unsupported on a given system.
        // Doing this check is a bit difficult at the moment without querying the EGL
        // extension string - we don't want to introduce a hard dependency on libEGL at
        // this point of development...
    }

    // Reset flipcounter and missed flip deadline counter to zero:
    windowRecord->flipCount = 0;
    windowRecord->nr_missed_deadlines = 0;
	
    //Return the window index and the rect argument.
    PsychCopyOutDoubleArg(1, FALSE, windowRecord->windowIndex);

    // Optionally return the windows clientrect:
    PsychCopyOutRectArg(2, FALSE, windowRecord->clientrect);

    return(PsychError_none);   
}

PsychError SCREENPanelFitter(void)
{
    static char useString1[] = "oldParams = Screen('PanelFitter', windowPtr [, newParams]);";
    static char synopsisString1[] =
        "Change operating parameters of builtin panel fitter.\n\n"
        "The size of the source framebuffer is given by the 'clientRect' parameter in Screen('OpenWindow'), "
        "the size of the destination framebuffer is given by the 'rect' parameter in that function. "
        "Default panel fitter behaviour is to rescale the source content to completely fit into the "
        "destination buffer, something that may not preserve aspect-ratio unless care is taken by the "
        "user to make sure source and destination framebuffer have already the same aspect ratio.\n"
        "This function allows to define new src and dst rectangles, thereby implicitely defining scaling "
        "and filtering properties. It optionally takes new settings in 'newParams' and returns old settings "
        "in 'oldParams'. The parameters are 8-element vectors of format\n"
        "params = [srcX0, srcY0, srcX1, srcY1, dstX0, dstY0, dstX1, dstY1];\n"
        "These tuples define top-left and bottom-right (x,y) corners of the source and destination "
        "rectangles for the (scaled)blit.\n"
        "You usually won't call this function directly, but leave the job to a higher-level setup "
        "routine, e.g., PsychImaging() and its 'UsePanelFitter' setup code.\n\n";
    static char seeAlsoString1[] = "OpenWindow";

    PsychWindowRecordType   *windowRecord;
    double* outParams;
    int*    newParams;
    int     count, i;
    
    // All sub functions should have these two lines
    PsychPushHelp(useString1, synopsisString1, seeAlsoString1);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    //cap the number of inputs
    PsychErrorExit(PsychCapNumInputArgs(2));   //The maximum number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));  //The maximum number of outputs

    // Get window record:
    PsychAllocInWindowRecordArg(1, TRUE, &windowRecord);
    
    // Return optional fitter settings:
    PsychAllocOutDoubleMatArg(1, FALSE, 1, 8, 1, &outParams);
    for (i = 0; i < 8; i++) outParams[i] = (double) windowRecord->panelFitterParams[i];
    
    // Get optional new panelFitter settings:
    if (PsychAllocInIntegerListArg(2, FALSE, &count, &newParams)) {
        if (count != 8) PsychErrorExitMsg(PsychError_user, "'newParams' must be a vector with 8 integer elements.");
        for (i = 0; i < count; i++) windowRecord->panelFitterParams[i] = newParams[i];

        // Fallback path and problematic new config setting?
        if (!(windowRecord->gfxcaps & kPsychGfxCapFBOBlit) && (PsychPrefStateGet_Verbosity() > 1) &&
            (windowRecord->panelFitterParams[0] != 0 || windowRecord->panelFitterParams[1] != 0 ||
             windowRecord->panelFitterParams[2] != (int) PsychGetWidthFromRect(windowRecord->clientrect) ||
             windowRecord->panelFitterParams[3] != (int) PsychGetHeightFromRect(windowRecord->clientrect))) {
            // Fallback path for panelFitter in use and sourceRegion is not == full clientRect. This is an
            // unsupported setting with the fallback, which will cause wrong results. Warn user:
            printf("PTB-WARNING: Selected a non-default srcRegion in Screen('PanelFitter'). This is not supported when\n");
            printf("PTB-WARNING: the fallback path for the panel fitter is in use. Expect distorted visual stimuli!\n");
        }
    }

    return(PsychError_none);
}
