/*
  Psychtoolbox3/Source/Common/Screen/SCREENResolutions.c		
  
  AUTHORS:

  mario.kleiner at tuebingen.mpg.de  mk
  
  PLATFORMS:	

  This file should build on any platform. 

  HISTORY:
  9/23/07  mk		Created. 
 
  DESCRIPTION:
  
  Contains the subfunctions 'ConfigureDisplay', 'Resolutions' and 'Resolution'.
  Enumerates and returns all supported display modes - combinations of resolution, refresh rate and color depth.
  Allows to query and set display mode.
  
  TO DO:

*/

#include "Screen.h"

#if PSYCH_SYSTEM == PSYCH_OSX
// Need these includes to make setting display brightness work:
#include <IOKit/graphics/IOGraphicsLib.h>
#endif

const char *FieldNames[]={"width", "height", "pixelSize", "hz"};

PsychError SCREENConfigureDisplay(void)
{
	static char useString[] = "oldSettings = Screen('ConfigureDisplay', setting, screenNumber, outputId [, newwidth][, newheight][, newHz][, newX][, newY]);";
	static char synopsisString[] =	"Query or change 'setting' for display output 'outputId' of screen 'screenNumber'.\n"
					"This function allows you to configure the different attached display outputs of a virtual screen.\n"
					"Optionally sets new settings for that output.\n"
					"Possible values for subfunction parameter 'setting':\n"
					"'Brightness': Return or set brightness of an attached display device. Many displays and systems don't support this function. "
                    " brightness values are in the range 0.0 to 1.0 from darkest to brightest. Returns old brightness setting.\n"
					"'NumberOutputs': Return number of active separate display outputs for given screen 'screenNumber'.\n"
					"'Capture': Capture output 'outputId' of a screen 'screenNumber' for exclusive use by Psychtoolbox.\n"
					"'Release': Release output 'outputId' of a screen 'screenNumber' from exclusive use by Psychtoolbox.\n"
					"Please note that 'Capture' and 'Release' are automatically applied to a display output as appropriate "
					"whenever a fullscreen onscreen window is opened or closed on a display. You usually don't need to call "
					"these functions yourself, in fact, you even shouldn't call them yourself usually, as wrong use of the "
					"functions can cause graphics malfunctions or even a freeze of the display or GUI, forcing you to kill "
					"Matlab or Octave as a last resort to regain control over your displays. The main purpose of these "
					"functions is to temporarily block all graphics operations or display updates on displays that are not "
					"used for visual stimulation, ie., not covered by fullscreen onscreen windows. This prevents other running "
					"applications from interfering with your experiment script by preventing them to take up valuable resources "
					"on the graphics card. Display capture is not supported on all operating systems in all modes of operation. "
					"If the functions get called on an unsupported configuration, they silently return and simply do nothing."
					"\n\n"
					"'Scanout': Retrieve or set scanout parameters for a given output 'outputId' of screen 'screenNumber'. "
					"Returns a struct 'oldSettings' with the current settings for that output. Only supported on Linux.\n"
					"It returns and accepts the following optional parameters:\n"
					"* Display resolution \"newwidth\" x \"newheight\", and nominal refresh rate \"newHz\". "
					"* Panning ('newX','newY') - The location of the top-left corner of the display in the framebuffer.\n"
					"Providing invalid or incompatible settings will raise an error.\n";

	static char seeAlsoString[] = "Screen('Resolutions'), Screen('Resolution');";

	const char *OutputFieldNames[]={"width", "height", "pixelSize", "hz", "xStart", "yStart"};
	char *settingName = NULL;
	int screenNumber, outputId;

    #if PSYCH_SYSTEM == PSYCH_LINUX
	PsychGenericScriptType *oldResStructArray;
	int newWidth, newHeight, newHz, newBpp, newX, newY;
	psych_bool rc;
    #endif

	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) { PsychGiveHelp(); return(PsychError_none); };

	// Check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumOutputArgs(1));
	PsychErrorExit(PsychCapNumInputArgs(8));
    
	// Get name of parameter class:
	PsychAllocInCharArg(1, kPsychArgRequired, &settingName);

	// Usercode wants to change display brightness:
	if (PsychMatch(settingName, "Brightness")) {
		// OS/X specific section:
		#if PSYCH_SYSTEM == PSYCH_OSX
		outputId = -1;
        
		// Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.  
		if (PsychCopyInScreenNumberArg(2, FALSE, &screenNumber)) {
			if(screenNumber==-1) PsychErrorExitMsg(PsychError_user, "Invalid screen number."); 
			outputId = screenNumber;
		}
        
		// Get outputId:
		PsychCopyInIntegerArg(3, FALSE, &outputId);
		if (outputId < 0 || outputId >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid video output specified!");
        
		CGDirectDisplayID displayID;
		CGDisplayErr err;
		float brightness;
		double nbrightness;
        
		PsychGetCGDisplayIDFromScreenNumber(&displayID, outputId);
		io_service_t service = CGDisplayIOServicePort(displayID);
        
		// Return current brightness value:
		err = IODisplayGetFloatParameter(service, kNilOptions, CFSTR(kIODisplayBrightnessKey), &brightness);
		if (err != kIOReturnSuccess) PsychErrorExitMsg(PsychError_user, "Failed to query current display brightness from system. Unsupported on this system.");
		PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) brightness);
        
		// Optionally set new brightness value:
		if (PsychCopyInDoubleArg(4, FALSE, &nbrightness)) {
			// Clamp to valid range:
			if (nbrightness < 0.0) nbrightness = 0.0;
			if (nbrightness > 1.0) nbrightness = 1.0;
            
			// Set it:
			err = IODisplaySetFloatParameter(service, kNilOptions, CFSTR(kIODisplayBrightnessKey), (float) nbrightness);
			if (err != kIOReturnSuccess) PsychErrorExitMsg(PsychError_user, "Failed to set new display brightness. Unsupported on this system?");
		}        
		#endif

		#if PSYCH_SYSTEM == PSYCH_LINUX

		static Atom backlight, backlight_new, backlight_legacy;
        
		// Get the screen number from the windowPtrOrScreenNumber. This also checks to make sure that the specified screen exists.
		PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
		if (screenNumber==-1) PsychErrorExitMsg(PsychError_user, "Invalid screen number."); 

		// Get outputId:
		outputId = 0;
		PsychCopyInIntegerArg(3, FALSE, &outputId);
        
		Display *dpy;
		int screen;
		float brightness;
		double nbrightness;

		// Map screenNumber and outputIdx to dpy, rootwindow and RandR output:
		PsychGetCGDisplayIDFromScreenNumber(&dpy, screenNumber);

        PsychLockDisplay();
		backlight_new    = XInternAtom(dpy, "Backlight", True);
		backlight_legacy = XInternAtom(dpy, "BACKLIGHT", True);
        PsychUnlockDisplay();

		if (backlight_new == None && backlight_legacy == None) {
			PsychErrorExitMsg(PsychError_user, "Failed to query current display brightness from system. System does not support brightness query and setting.");
		}

		screen = PsychGetXScreenIdForScreen(screenNumber);

        PsychLockDisplay();
		Window root = RootWindow(dpy, screen);
		XRRScreenResources *resources = XRRGetScreenResources(dpy, root);
        PsychUnlockDisplay();

		if (resources == NULL) PsychErrorExitMsg(PsychError_user, "Failed to query current display brightness from system. Feature not supported.");
		if (outputId < 0 || outputId >= resources->noutput) PsychErrorExitMsg(PsychError_user, "Invalid video output specified!");
		RROutput output = resources->outputs[outputId];

		// Query current brightness of output, if possible. Bail otherwise:
		unsigned long nitems;
		unsigned long bytes_after;
		unsigned char *prop;
		Atom actual_type;
		int actual_format;
		long value;
		XRRPropertyInfo *info;
		double cur, new;
		double min, max;
    
        PsychLockDisplay();

		backlight = backlight_new;
		if (!backlight || XRRGetOutputProperty(dpy, output, backlight, 0, 4, False, False, None, &actual_type, &actual_format, &nitems, &bytes_after, &prop) != Success) {
			backlight = backlight_legacy;
			if (!backlight || XRRGetOutputProperty(dpy, output, backlight, 0, 4, False, False, None, &actual_type, &actual_format, &nitems, &bytes_after, &prop) != Success) {
                PsychUnlockDisplay();
				PsychErrorExitMsg(PsychError_user, "Failed to query current display brightness from system. Unsupported feature?");
            }
		}

		if (actual_type != XA_INTEGER || nitems != 1 || actual_format != 32) {
			XFree(prop);
            PsychUnlockDisplay();
			PsychErrorExitMsg(PsychError_user, "Failed to query current display brightness from system. Unsupported feature?");
		} else {
			value = *((long *) prop);
			XFree(prop);
		}

		info = XRRQueryOutputProperty(dpy, output, backlight);
		if (info && (info->range) && (info->num_values == 2)) {
			min = info->values[0];
			max = info->values[1];
			cur = (double) value;

			// Return current value:
			PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) ((cur - min) / (max - min)));

			// Optionally set new brightness value:
			if (PsychCopyInDoubleArg(4, FALSE, &nbrightness)) {
				// Clamp to valid range:
				if (nbrightness < 0.0) nbrightness = 0.0;
				if (nbrightness > 1.0) nbrightness = 1.0;            

				// Map to raw value:
				new = (nbrightness * (max - min)) + min;

				// Be paranoid, clamp again:
				if (new > max) new = max;
				if (new < min) new = min;

				// Set new brightness value:
				value = (long) new;
				XRRChangeOutputProperty(dpy, output, backlight, XA_INTEGER, 32, PropModeReplace, (unsigned char *) &value, 1);

				// Kick the server to do its job quickly:
				XFlush(dpy);
			}
			XFree(info);
		} else {
			if (info) XFree(info);
            PsychUnlockDisplay();
			PsychErrorExitMsg(PsychError_user, "Failed to query current display brightness from system. Unsupported feature?");
		}

		XRRFreeScreenResources(resources);

        PsychUnlockDisplay();

		#endif

		#if PSYCH_SYSTEM == PSYCH_WINDOWS
		PsychErrorExitMsg(PsychError_unimplemented, "Sorry, this function is not implemented on MS-Windows.");
		#endif

		return(PsychError_none);
	}
    
	// Usercode wants to know number of outputs for a screen?
	if (PsychMatch(settingName, "NumberOutputs")) {
		// Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.  
		PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
		if(screenNumber==-1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number."); 

		// Count and return assigned outputs for screen:
		for (outputId = 0; (outputId < kPsychMaxPossibleCrtcs) && (PsychScreenToHead(screenNumber, outputId) >= 0); outputId++);
		PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) outputId);

		return(PsychError_none);
	}

	// Usercode wants to capture a display screen?
	if (PsychMatch(settingName, "Capture")) {
		// Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.  
		PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
		if(screenNumber==-1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number."); 

		// This functions is only supported on OSX for the time being, a no-op on other platforms:
		if ((PSYCH_SYSTEM == PSYCH_OSX) && !PsychIsScreenCaptured(screenNumber)) PsychCaptureScreen(screenNumber);

		return(PsychError_none);
	}

	// Usercode wants to release a display screen?
	if (PsychMatch(settingName, "Release")) {
		// Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.  
		PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
		if(screenNumber==-1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number."); 

		// This functions is only supported on OSX for the time being, a no-op on other platforms:
		if ((PSYCH_SYSTEM == PSYCH_OSX) && PsychIsScreenCaptured(screenNumber)) PsychReleaseScreen(screenNumber);

		return(PsychError_none);
	}

#if PSYCH_SYSTEM != PSYCH_LINUX
	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, this function is only supported on Linux.");
#else

	if(!PsychMatch(settingName, "Scanout")) PsychErrorExitMsg(PsychError_user, "Unknown 'setting' name provided. Typo?");

	// Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.  
	PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
	if(screenNumber==-1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number."); 

	PsychCopyInIntegerArg(3, TRUE, &outputId);
	if(outputId < 0 || outputId >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid display output index provided."); 

	// Create a structure and populate it.
	PsychAllocOutStructArray(1, FALSE, 1, 6, OutputFieldNames, &oldResStructArray);

	// Query current video mode of this output:
	XRRCrtcInfo *crtc_info = NULL;

    PsychLockDisplay();
	XRRModeInfo *mode = PsychOSGetModeLine(screenNumber, outputId, &crtc_info);
    PsychUnlockDisplay();

	if (NULL == mode) PsychErrorExitMsg(PsychError_user, "Could not query video mode for this output. Invalid outputId or unsupported function on this system?");

	// Get (x,y) top-left corner of crtc's viewport -- panning info:
	PsychSetStructArrayDoubleElement("xStart", 0, (double) crtc_info->x, oldResStructArray);
	PsychSetStructArrayDoubleElement("yStart", 0, (double) crtc_info->y, oldResStructArray);
	XRRFreeCrtcInfo(crtc_info);

	// Query and return resolution:
	newWidth = (int) mode->width;
	newHeight = (int) mode->height;
	PsychSetStructArrayDoubleElement("width", 0, newWidth, oldResStructArray);
	PsychSetStructArrayDoubleElement("height", 0, newHeight, oldResStructArray);

	// Query and return refresh rate:

    // Modeline with plausible values returned by RandR?
    if (mode && (mode->hTotal > mode->width) && (mode->vTotal > mode->height)) {
        // Yes: use RandR results:
        newHz = (int) (PsychOSVRefreshFromMode(mode) + 0.5);
    } else {
        // No: Fallback to old style - Refresh rate of primary output for this screen, following
        //     the "whatever works" method:
        newHz = (int) (PsychGetNominalFramerate(screenNumber) + 0.5);
    }

	PsychSetStructArrayDoubleElement("hz", 0, newHz, oldResStructArray);
	
	// Query and return current display depth:
	newBpp = PsychGetScreenDepthValue(screenNumber);
	PsychSetStructArrayDoubleElement("pixelSize", 0, newBpp, oldResStructArray);

	// Any new settings provided? Otherwise we skip this:
	if(PsychGetNumInputArgs() > 3) {
		// Make sure we're not called while onscreen windows are open:
		if ((PsychCountOpenWindows(kPsychSingleBufferOnscreen) + PsychCountOpenWindows(kPsychDoubleBufferOnscreen)) > 0) PsychErrorExitMsg(PsychError_user, "Tried to change video display settings via Screen('ConfigureDisplay'); while onscreen windows were open! Not allowed.");

		// Copy in optional new settings:
		PsychCopyInIntegerArg(4, FALSE, &newWidth);
		PsychCopyInIntegerArg(5, FALSE, &newHeight);
		PsychCopyInIntegerArg(6, FALSE, &newHz);
		newX = newY = -1;
		PsychCopyInIntegerArg(7, FALSE, &newX);
		PsychCopyInIntegerArg(8, FALSE, &newY);

		rc = PsychOSSetOutputConfig(screenNumber, outputId, newWidth, newHeight, newHz, newX, newY);
		// Check if successfull:
		if (!rc) PsychErrorExitMsg(PsychError_user, "Invalid or mutually incompatible video settings requested!\nOne or more of the values are invalid or unsupported by your display device.");
	}

#endif
	return(PsychError_none);	
}

PsychError SCREENResolution(void) 
{
	static char useString[] = "oldResolution=Screen('Resolution', screenNumber [, newwidth][, newheight][, newHz][, newPixelSize][, specialMode]);";
	static char synopsisString[] =	"Query or change display settings for screen \"screenNumber\".\n"
					"Returns a struct \"oldResolutions\" with the current settings for screen "
					"resolution, refresh rate and pixel depth. Optionally sets new settings "
					"for screen resolution \"newwidth\" x \"newheight\", refresh rate \"newHz\" "
					"and framebuffer pixel depth \"newPixelSize\". Providing invalid or incompatible "
					"settings will raise an error. Especially the color depth \"newPixelSize\" should "
					"usually not be set to anything else than its default of 32 bpp or 24 bpp. "
					"Other settings can impair alpha-blending on some systems, a setting of 16 bpp "
					"will disable alpha-blending and create drastically reduced color resolution of "
					" 5 bits per color channel. "
					"A setting of 8 bpp is not supported at all on MacOS/X and will create artifacts "
					"on all other systems. Use a size of 32 bpp even for clut animation. This function "
					"may not work on all MS-Windows setups, your mileage may vary...\n"
					"On Linux the function only switches display settings in the conventional sense on "
					"a single-display setup. On a multi-display setup, this function only changes the "
					"total size of the framebuffer, ie., 'newwidth' and 'newheight', the other "
					"parameters are silently ignored. On Linux, the video settings of each individual display, "
					"e.g., resolution, video refresh rate, panning, are queried and changed via the "
					"Screen('ConfigureDisplay') function instead. This allows for much more flexibility, "
					"e.g., you can have a framebuffer bigger than the combined resolution of all displays "
					"and only show a fraction of it. You can change the relative position of all physical "
					"displays, configure \"mirror modes\", \"side by side\", or \"on top of each other\" "
					"display configurations.\n"
					"Psychtoolbox will automatically restore the systems display resolution to the "
					"system settings made via the display control panel as soon as either your script "
					"finishes by closing all its windows, or by some error. Terminating Matlab due to "
					"quit command will also restore the system preference settings. On a multi-display "
					"Linux setup, display settings are never automatically restored.\n"
					"If you call this command without ever opening onscreen windows and closing them "
					"at some point, Psychtoolbox will not restore display settings automatically.\n"
					"You can query a list of all supported combinations of display settings via the "
					"Screen('Resolutions') command. \"specialMode\" is a flag you must not touch, "
					"unless you really know what you're doing, that's why we don't tell you its purpose.";

	static char seeAlsoString[] = "Screen('Resolutions')";
	
	PsychGenericScriptType *oldResStructArray;
	PsychScreenSettingsType screenSettings;
	PsychDepthType useDepth;
	int screenNumber, specialMode;
	int newWidth, newHeight, newHz, newBpp;
	long pnewWidth, pnewHeight;
	psych_bool rc;
	
	// Purpose of 'specialMode': If bit zero is set, then its possible to switch display settings while
	// onscreen windows - possibly fullscreen windows and display capture - are open/active. This is mostly
	// useful for temporarilly changing display framerate, e.g., from a wanted value of x Hz to a intermittent
	// value of y Hz, then back to x Hz. Such an approach may be useful for "drift-syncing" multiple displays
	// whose video refresh cycles are not perfectly in sync (not yet clear if this will really work), but that's
	// about the only useful purpose. Changing other settings - or changing any settings permanently - while
	// windows are already open will likely subvert display calibration and lead to severe timing problems
	// and possible worse things - like complete hangs of the graphics subsystem and the need for a hard
	// machine reset!
	// Other bits of specialMode will trigger other sync related actions - implementation and semantics may
	// change without notice!
	
	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) { PsychGiveHelp(); return(PsychError_none); };

	// Check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumOutputArgs(1));
	PsychErrorExit(PsychCapNumInputArgs(6));
    
	// Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.  
	PsychCopyInScreenNumberArg(kPsychUseDefaultArgPosition, TRUE, &screenNumber);
	if(screenNumber==-1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number."); 

	// Create a structure and populate it.
	PsychAllocOutStructArray(1, FALSE, 1, 4, FieldNames, &oldResStructArray);

	// Query and return resolution:
	PsychGetScreenSize(screenNumber, &pnewWidth, &pnewHeight);
	newWidth = (int) pnewWidth;
	newHeight = (int) pnewHeight;
	PsychSetStructArrayDoubleElement("width", 0, newWidth, oldResStructArray);
	PsychSetStructArrayDoubleElement("height", 0, newHeight, oldResStructArray);

	// Query and return refresh rate:
	newHz = (int) (PsychGetNominalFramerate(screenNumber) + 0.5);
	PsychSetStructArrayDoubleElement("hz", 0, newHz, oldResStructArray);
	
	// Query and return current display depth:
	newBpp = PsychGetScreenDepthValue(screenNumber);
	PsychSetStructArrayDoubleElement("pixelSize", 0, newBpp, oldResStructArray);

	// Any new settings provided? Otherwise we skip this:
	if(PsychGetNumInputArgs() > 1) {
		// Get optional specialMode flag, default to zero:
		specialMode = 0;
		PsychCopyInIntegerArg(6, FALSE, &specialMode);

		// Make sure we're not called while onscreen windows are open. Only allow calling with open onscreen windows if Bit zero of specialMode is set:
		if (((PsychCountOpenWindows(kPsychSingleBufferOnscreen) + PsychCountOpenWindows(kPsychDoubleBufferOnscreen)) > 0) && !(specialMode & 1)) PsychErrorExitMsg(PsychError_user, "Tried to change video display settings via Screen('Resolutions'); while onscreen windows were open! Not allowed.");

		// Copy in optional new settings:
		PsychCopyInIntegerArg(2, FALSE, &newWidth);
		PsychCopyInIntegerArg(3, FALSE, &newHeight);
		PsychCopyInIntegerArg(4, FALSE, &newHz);
		PsychCopyInDepthValueArg(5, FALSE, &newBpp);
		PsychInitDepthStruct(&useDepth);
		PsychAddValueToDepthStruct(newBpp, &useDepth);
		
		// Switch to new display mode, according to these specs:
		PsychGetScreenSettings(screenNumber, &screenSettings);
		screenSettings.rect[kPsychLeft] = 0;
		screenSettings.rect[kPsychTop] = 0;
		screenSettings.rect[kPsychRight] = newWidth;
		screenSettings.rect[kPsychBottom] = newHeight;
		screenSettings.nominalFrameRate = newHz;
		PsychInitDepthStruct(&(screenSettings.depth));
		PsychCopyDepthStruct(&(screenSettings.depth), &useDepth);

		// Perform sanity check:
		if ((PSYCH_SYSTEM != PSYCH_LINUX) && !PsychCheckVideoSettings(&screenSettings)) PsychErrorExitMsg(PsychError_user, "Invalid or mutually incompatible video settings requested!\nOne or more of the values are invalid or unsupported by your display device.");
		
		// Perform actual switch:
		if (!(specialMode & 1)) PsychCaptureScreen(screenNumber);
		rc = PsychSetScreenSettings( TRUE, &screenSettings);
		if (!(specialMode & 1)) PsychReleaseScreen(screenNumber);

		// Check if successfull:
		if (!rc) PsychErrorExitMsg(PsychError_user, "Invalid or mutually incompatible video settings requested!\nOne or more of the values are invalid or unsupported by your display device.");
	}

	return(PsychError_none);	
}

PsychError SCREENResolutions(void) 
{
	static char useString[] = "resolutions = Screen('Resolutions', screenNumber [, outputId]);";
	static char synopsisString[] =	"Query a list of all supported and valid display settings for screen "
					"\"screenNumber\" and display output \"outputId\". If the optional "
					"'outputId' is omitted, the unified settings of the screen are returned. "
					"Currently 'outputId' is only honored on Linux, ignored on other systems.\n"
					"You can set your display to one of the supported "
					"combinations of settings via the Screen('Resolution') command.\n"
					"Returns an array of structs \"resolutions\", where "
					"each element in the array is a struct that describes one valid "
					"combination of resolution, color depth and refresh rate. Fields are "
					"self explanatory.\n"
					"Please note that unless you have good reason to do so, especially the "
					"color depth value \"newPixelSize\" should usually not be changed. "
					"Usually it is 32 bpp or 24 bpp. "
					"A setting of 16 bpp will disable alpha-blending and create "
					"drastically reduced color resolution of 5 bits per color channel. A "
					"setting of 8 bpp is not supported at all on MacOS/X and will create artifacts "
					"on all other systems. Use a size of 32 bpp even for clut animation. This function "
					"may not work on all MS-Windows setups, your mileage may vary.\n"
					"Please note that there are a couple of helper functions in the PsychOneLiners "
					"directory of Psychtoolbox which can simplify the task of switching resolutions. "
					"That functions are probably more convenient to use than the low-level Screen "
					"functions for display settings... ";
	static char seeAlsoString[] = "Screen('Resolution')";
	
	PsychGenericScriptType *ResStructArray;
	int screenNumber, outputid, nrmodes, i;
	long *widths, *heights, *hz, *bpp;
	
	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) { PsychGiveHelp(); return(PsychError_none); };

	// Check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumOutputArgs(1));
	PsychErrorExit(PsychCapNumInputArgs(2));
    
	// Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.  
	PsychCopyInScreenNumberArg(kPsychUseDefaultArgPosition, TRUE, &screenNumber);
	if(screenNumber==-1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number."); 

	outputid = -1;
	PsychCopyInIntegerArg(2, kPsychArgOptional, &outputid);

	// Query all supported settings:
	nrmodes = PsychGetAllSupportedScreenSettings(screenNumber, outputid, &widths, &heights, &hz, &bpp);
	
	// Create a structure array and populate it.
	PsychAllocOutStructArray(1, FALSE, nrmodes, 4, FieldNames, &ResStructArray);

	for (i=0; i<nrmodes; i++) {
		PsychSetStructArrayDoubleElement("width", i, (double) widths[i], ResStructArray);
		PsychSetStructArrayDoubleElement("height", i, (double) heights[i], ResStructArray);
		PsychSetStructArrayDoubleElement("hz", i, (double) hz[i], ResStructArray);
		PsychSetStructArrayDoubleElement("pixelSize", i, (double) bpp[i], ResStructArray);
	}
	
	// Ready.
	return(PsychError_none);	
}
