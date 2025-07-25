/*
    Psychtoolbox3/Source/Common/Screen/SCREENResolutions.c

    AUTHORS:

    mario.kleiner.de@gmail.com      mk

    PLATFORMS:

    This file should build on any platform.

    HISTORY:

    9/23/07  mk Created.

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
# include <ApplicationServices/ApplicationServices.h>

/*
    The code for using System private interfaces (SPI's) for getting/setting display
    brightness on Apple's trainwreck operating system was adapted from the code of
    Nicholas Riley at https://github.com/nriley/brightness which is licensed under
    the following BSD-2 clause license (Thank you!):

    Copyright (c) 2014-2019, Nicholas Riley
    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification,
    are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice, this
    list of conditions and the following disclaimer in the documentation and/or
    other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
     LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* As of macOS 10.12.4, brightness set by public IOKit API is
 overridden by CoreDisplay's brightness (to support Night Shift). In
 addition, using CoreDisplay to get brightness supports additional
 display types, e.g. the 2015 iMac's internal display.
 The below functions in CoreDisplay seem to work to adjust the
 "user" brightness (like dragging the slider in System Preferences
 or using function keys).  The symbols below are listed in the .tbd
 file in CoreDisplay.framework so it is at least more "public" than
 a symbol in a private framework, though there are no public headers
 distributed for this framework. */
extern double CoreDisplay_Display_GetUserBrightness(CGDirectDisplayID id) __attribute__((weak_import));
extern void CoreDisplay_Display_SetUserBrightness(CGDirectDisplayID id, double brightness) __attribute__((weak_import));

/*  Some issues with the above CoreDisplay functions include:
     - There's no way to tell if setting the brightness was successful
     - There's no way to tell if a brightness of 1 means that the
     brightness is actually 1, or if there's no adjustable brightness
     - Brightness changes aren't reflected in System Preferences
     immediately
     Fixing these means using the private DisplayServices.framework.  Be
     even more careful about these.
 */
extern _Bool DisplayServicesCanChangeBrightness(CGDirectDisplayID id) __attribute__((weak_import));
extern void DisplayServicesBrightnessChanged(CGDirectDisplayID id, double brightness) __attribute__((weak_import));

extern int DisplayServicesGetBrightness(CGDirectDisplayID id, float *brightness) __attribute__((weak_import));
extern int DisplayServicesSetBrightness(CGDirectDisplayID id, float brightness) __attribute__((weak_import));
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
                    " To find out if your system supports brightness queries and setting, call this subfunction without specifying a brightness. "
                    " The function will either return a >= brightness value, or -1 if the system doesn't support brightness control. "
                    " If you try to set brightness on a system that doesn't support it, the function will abort with an error. "
                    " Brightness values are in the range 0.0 to 1.0 from darkest to brightest. Returns old brightness setting.\n"
                    "'AutoBrightnessMode': Set the AutoBrightness mode on supported systems. Doesn't work anywhere at the moment. "
                    "If you call the function without specifying a new mode, it will return -1 if setting AutoBrightness is "
                    "unsupported on the system, otherwise it will return the value 2 to confirm setting AutoBrightness should be possible.\n"
                    "If you try to set AutoBrightness on a system that doesn't support it, the function will abort with an error. "
                    "Otherwise it will return 2 to confirm it tried to set AutoBrightness mode and then try to set the mode you "
                    "specified: A mode of 1 tries to enable AutoBrightness, a mode of 0 tries to disable AutoBrightness on the given output. "
                    "Please note that there isn't any feedback if switching AutoBrightness worked, or what the current or previous enable "
                    "state of AutoBrightness is or was.\n"
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
                    "'Dithering' Control digital display dithering on supported GPUs. This is currently only supported on AMD "
                    "graphics cards under Linux and OSX. 'screenNumber' will apply the setting on all video outputs connected "
                    "to screen 'screenNumber'. The next setting is the 'ditherEnable' flag: A value of zero will disable dithering. "
                    "A non-zero value will reenable dithering if it was previous disabled by Psychtoolbox. If dithering wasn't "
                    "disabled, then instead the given value will be written into the dither control register to modify dithering "
                    "settings. The meaning of specific values is GPU and vendor specific and only useful for experts.\n"
                    "Example call: Screen('ConfigureDisplay', 'Dithering', screenNumber, ditherEnable); "
                    "\n\n"
                    "'RequestMinimumOutputPrecision' Try to request a specific minimum output depth 'minBpc' for video signals on a "
                    "X-Screen under Linux/X11. This is currently only supported on Linux/X11, and it does not guarantee that "
                    "you will get the requested minimum precision, as various hardware constraints in the graphics card, cabling, "
                    "display device, and available memory bandwidth can lead to a lower precision than requested. In the end it "
                    "is only a polite request to the graphics driver. Example call:\n"
                    "Screen('ConfigureDisplay', 'RequestMinimumOutputPrecision', screenNumber, minBpc); "
                    "\n\n"
                    "'FineGrainedSwitchRefreshRate': Switch the video refresh rate of a given output 'outputId' of screen 'screenNumber'.\n"
                    "On Linux/X11 this function tries to quickly switch the video refresh rate of an output in a seamless way "
                    "with no or only minimum visual disruption, in a tiny fraction of a second, e.g., at most one video refresh cycle. "
                    "Allow for very fine-grained choice of refresh rate, e.g., fractions of 1 Hz. This generally only works well on some "
                    "graphics cards and Linux kernels while driving Variable Refresh Rate capable displays (FreeSync, GSync, DP adaptive sync).\n"
                    "One limitation at the moment is that in a multi-display setup, your visual stimulation window must be "
                    "fully or at least predominantly covering the display area of the primary video output monitor, as the "
                    "refresh rate of only that monitor will be assigned for Screen('Flip', ...) scheduling of stimulus onset time. "
                    "Otherwise only immediate flips at each video refresh cycle will work with reliable timing. "
                    "See 'help VRRSupport' for possibly more background and setup info.\n"
                    "As of the year 2025 this works with FreeSync capable AMD graphics cards and display devices connected via DisplayPort. "
                    "Other graphics cards or display devices (even HDMI VRR!) may not allow fast seamless switching. Example call:\n"
                    "actualHz = Screen('ConfigureDisplay', 'FineGrainedSwitchRefreshRate', screenNumber, outputId, requestedHz); "
                    "\n\n"
                    "'Scanout': Retrieve or set scanout parameters for a given output 'outputId' of screen 'screenNumber'. "
                    "Returns a struct 'oldSettings' with the current settings for that output. Only supported on Linux.\n"
                    "It returns and accepts the following optional parameters:\n\n"
                    "* Display resolution \"newwidth\" x \"newheight\", and nominal refresh rate \"newHz\".\n"
                    "* Panning ('newX','newY') - The location of the top-left corner of the display in the framebuffer.\n"
                    "Readonly properties:\n\n"
                    "* Display size \"displayWidthMM\" x \"displayHeightMM\" in millimeters, as reported by attached display.\n"
                    "* Video output name \"name\".\n"
                    "* Video output handle \"outputHandle\". On Linux/X11 this is the numeric XID of the associated RandR video output.\n"
                    "\nProviding invalid or incompatible settings will raise an error.\n";

    static char seeAlsoString[] = "Screen('Resolutions'), Screen('Resolution');";

    const char *OutputFieldNames[] = { "width", "height", "pixelSize", "hz", "xStart", "yStart", "name", "outputHandle", "displayWidthMM", "displayHeightMM" };
    char *settingName = NULL;
    int screenNumber, outputId, ditherEnable;

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

    if (PsychMatch(settingName, "AutoBrightnessMode")) {
        // Unsupported on all operating systems, including on macOS, where it never ever worked, to my knowledge.
        // Return -1 for "not supported" if this was just a query, but fail if code tries to set a new value.
        // This way a failure to change autobrightness won't go unnoticed, but code can query first if setting
        // autobrightness is supported before it trusts results:
        if (PsychIsArgPresent(PsychArgIn, 4))
            PsychErrorExitMsg(PsychError_user, "Failed to change current display auto-brightness mode. Unsupported on this system.");
        else
            PsychCopyOutDoubleArg(1, kPsychArgOptional, -1);

        return(PsychError_none);
    }

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

        // Query possible and meaningful?
        if ((DisplayServicesCanChangeBrightness != NULL) && !DisplayServicesCanChangeBrightness(displayID)) {
            // Nope. This means "unsupported". Return -1 for "not supported" if this was just a query,
            // but fail if code tried to set a new value. This way a failure to set a brightness level
            // won't go unnoticed, but code can query first if setting/getting brightness is supported
            // before it trusts results:
            goto brightness_unsupported;
        }

        // Return current brightness value: Prefer SPI (System private interface) for this,
        // as it isn't as fucked up as Apple's public implementation:
        if ((DisplayServicesGetBrightness != NULL) && ((err = DisplayServicesGetBrightness(displayID, &brightness)) == 0)) {
            if (PsychPrefStateGet_Verbosity() > 4)
                printf("PTB-INFO: Current display brightness %f queried via system SPI DisplayServicesGetBrightness().\n", brightness);
        }
        else if (CoreDisplay_Display_GetUserBrightness != NULL) {
            if (PsychPrefStateGet_Verbosity() > 4)
                printf("PTB-INFO: Display brightness system SPI DisplayServicesGetBrightness() unsupported [%p] or failed [%i]. Trying IntelMac SPI fallback.\n",
                       DisplayServicesGetBrightness, (int) err);

            // Supported. Query it:
            brightness = CoreDisplay_Display_GetUserBrightness(displayID);
        }
        else {
            goto brightness_unsupported;
        }

        // Return old / current brightness value:
        PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) brightness);

        // Optionally set new brightness value:
        if (PsychCopyInDoubleArg(4, FALSE, &nbrightness)) {
            // Clamp to valid range:
            if (nbrightness < 0.0) nbrightness = 0.0;
            if (nbrightness > 1.0) nbrightness = 1.0;

            // Set it: Prefer SPI (System private interface) for this, as it isn't as fucked up as Apple's public implementation:
            if ((DisplayServicesSetBrightness != NULL) && ((err = DisplayServicesSetBrightness(displayID, (float) nbrightness)) == 0)) {
                if (PsychPrefStateGet_Verbosity() > 4)
                    printf("PTB-INFO: New display brightness %f set via system SPI DisplayServicesSetBrightness().\n", nbrightness);
            }
            else if (CoreDisplay_Display_SetUserBrightness != NULL) {
                if (PsychPrefStateGet_Verbosity() > 4)
                    printf("PTB-INFO: Display brightness system SPI DisplayServicesSetBrightness() unsupported [%p] or failed [%i]. Trying IntelMac SPI fallback.\n",
                           DisplayServicesSetBrightness, (int) err);

                // Fallback, only works on IntelMacs, no way to figure out if it worked:
                CoreDisplay_Display_SetUserBrightness(displayID, nbrightness);

                // Notify UI, this function only exists on Intel Macs, not on Apple Silicon:
                #ifdef __amd64__
                if (DisplayServicesBrightnessChanged != NULL)
                    DisplayServicesBrightnessChanged(displayID, nbrightness);
                #endif
            }
            else {
                PsychErrorExitMsg(PsychError_user, "Failed to set new display brightness due to unsupported system functionality.");
            }
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

        CGDirectDisplayID dpy;
        int screen;
        double nbrightness;
        XRRScreenResources *resources = NULL;

        // Map screenNumber and outputIdx to dpy, rootwindow and RandR output:
        PsychGetCGDisplayIDFromScreenNumber(&dpy, screenNumber);

        PsychLockDisplay();
        backlight_new    = XInternAtom(dpy, "Backlight", True);
        backlight_legacy = XInternAtom(dpy, "BACKLIGHT", True);
        PsychUnlockDisplay();

        if (backlight_new == None && backlight_legacy == None)
            goto brightness_unsupported;

        screen = PsychGetXScreenIdForScreen(screenNumber);

        PsychLockDisplay();
        Window root = RootWindow(dpy, screen);
        resources = XRRGetScreenResources(dpy, root);
        PsychUnlockDisplay();

        if (resources == NULL)
            goto brightness_unsupported;

        if (outputId < 0 || outputId >= resources->noutput)
            PsychErrorExitMsg(PsychError_user, "Invalid video output specified!");

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
                goto brightness_unsupported;
            }
        }

        if (actual_type != XA_INTEGER || nitems != 1 || actual_format != 32) {
            XFree(prop);
            PsychUnlockDisplay();
            goto brightness_unsupported;
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
            goto brightness_unsupported;
        }

        XRRFreeScreenResources(resources);

        PsychUnlockDisplay();

        #endif

        #if PSYCH_SYSTEM == PSYCH_WINDOWS
        goto brightness_unsupported;
        #endif

        return(PsychError_none);

        brightness_unsupported:
            #if PSYCH_SYSTEM == PSYCH_LINUX
            if (resources)
                XRRFreeScreenResources(resources);
            #endif

            // Failed. This means "unsupported". Return -1 for "not supported" if this was just a query,
            // but fail if code tried to set a new value. This way a failure to set a brightness level
            // won't go unnoticed, but code can query first if setting/getting brightness is supported
            // before it trusts results:
            if (PsychIsArgPresent(PsychArgIn, 4))
                PsychErrorExitMsg(PsychError_user, "Failed to query and/or set current display brightness from system. Unsupported on this system.");
            else
                PsychCopyOutDoubleArg(1, kPsychArgOptional, -1);

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

    // Usercode wants to manually control dithering on an output?
    if (PsychMatch(settingName, "Dithering")) {
        // Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.
        PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
        if (screenNumber==-1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number.");

        // Get ditherEnable flag: 0 = Disable dithering. Other values may either trigger reenable, or
        // dither register setup with given value:
        PsychCopyInIntegerArg(3, TRUE, &ditherEnable);

        // Enable or disable dithering on all video outputs associated with screen 'screenNumber':
        PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) PsychSetOutputDithering(NULL, screenNumber, ditherEnable));

        return(PsychError_none);
    }

#if PSYCH_SYSTEM != PSYCH_LINUX
    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, unless you made a typo in the function name, this function is likely only supported on Linux.");
#else

    // Usercode wants to manually request a higher video output bit depth on an output?
    if (PsychMatch(settingName, "RequestMinimumOutputPrecision")) {
        int minBpc;

        // Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.
        PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
        if (screenNumber == -1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number.");

        // Get minimum bpc for video output:
        PsychCopyInIntegerArg(3, TRUE, &minBpc);
        if (minBpc < 8) PsychErrorExitMsg(PsychError_user, "Invalid minimum output bit depth requested. Must be at least 8.");

        // Running under Linux/X11 native?
        #if !defined(PTB_USE_WAFFLE)
            // Request at least minBpc video signal precision on all video outputs associated with screen 'screenNumber':
            // Note: Actual precision can be lower, subject to constraints of the video output, cable, display sink and hardware,
            // also due to video bandwidth constraints.
            PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) PsychOSEnsureMinimumOutputPrecision(screenNumber, minBpc));
        #else
            // Nope: Report failure back:
            PsychCopyOutDoubleArg(1, kPsychArgOptional, 0);
        #endif

        return(PsychError_none);
    }

    if (PsychMatch(settingName, "FineGrainedSwitchRefreshRate")) {
        double requestedHz;

        // Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.
        PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
        if (screenNumber == -1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number.");

        PsychCopyInIntegerArg(3, TRUE, &outputId);
        if(outputId < 0 || outputId >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid display output index provided.");

        // Get minimum bpc for video output:
        PsychCopyInDoubleArg(4, TRUE, &requestedHz);
        if (requestedHz < 1) PsychErrorExitMsg(PsychError_user, "Invalid video refresh rate requested. Must be at least 1 Hz.");

        // Running under Linux/X11 native?
        #if !defined(PTB_USE_WAFFLE)
        PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychOSRandRCreateAndSetFRRVRRMode(screenNumber, outputId, requestedHz));
        #else
        // Nope: Report failure back:
        PsychCopyOutDoubleArg(1, kPsychArgOptional, 0);
        #endif

        return(PsychError_none);
    }

    if(!PsychMatch(settingName, "Scanout")) PsychErrorExitMsg(PsychError_user, "Unknown 'setting' name provided. Typo?");

    // Get the screen number from the windowPtrOrScreenNumber.  This also checks to make sure that the specified screen exists.
    PsychCopyInScreenNumberArg(2, TRUE, &screenNumber);
    if(screenNumber==-1) PsychErrorExitMsg(PsychError_user, "The specified onscreen window has no ancestral screen or invalid screen number.");

    PsychCopyInIntegerArg(3, TRUE, &outputId);
    if(outputId < 0 || outputId >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid display output index provided.");

    // Create a structure and populate it.
    PsychAllocOutStructArray(1, FALSE, -1, 10, OutputFieldNames, &oldResStructArray);

    // Query and return output name, and physical size of display in mm:
    unsigned long mm_width, mm_height;
    RROutput output_xid;

    PsychSetStructArrayStringElement("name", 0, (char*) PsychOSGetOutputProps(screenNumber, outputId, FALSE, &mm_width, &mm_height, &output_xid), oldResStructArray);
    PsychSetStructArrayDoubleElement("outputHandle", 0, (double) output_xid, oldResStructArray);
    PsychSetStructArrayDoubleElement("displayWidthMM", 0, mm_width, oldResStructArray);
    PsychSetStructArrayDoubleElement("displayHeightMM", 0, mm_height, oldResStructArray);

    // Query current video mode of this output:
    XRRCrtcInfo *crtc_info = NULL;

    PsychLockDisplay();
    XRRModeInfo *mode = PsychOSGetModeLine(screenNumber, outputId, &crtc_info);
    PsychUnlockDisplay();

    // Get (x,y) top-left corner of crtc's viewport -- panning info:
    PsychSetStructArrayDoubleElement("xStart", 0, (double) ((crtc_info) ? crtc_info->x : 0), oldResStructArray);
    PsychSetStructArrayDoubleElement("yStart", 0, (double) ((crtc_info) ? crtc_info->y : 0), oldResStructArray);
    if (crtc_info)
        XRRFreeCrtcInfo(crtc_info);

    crtc_info = NULL;

    // Note: mode would be NULL on a disabled output.

    // Query and return resolution:
    newWidth = (int) ((mode) ? mode->width : 0);
    newHeight = (int) ((mode) ? mode->height : 0);
    PsychSetStructArrayDoubleElement("width", 0, newWidth, oldResStructArray);
    PsychSetStructArrayDoubleElement("height", 0, newHeight, oldResStructArray);

    // Query and return refresh rate:

    // Modeline with plausible values returned by RandR?
    if (mode && (mode->hTotal >= mode->width) && (mode->vTotal >= mode->height)) {
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
    static char synopsisString[] = "Query or change display settings for screen \"screenNumber\".\n"
                    "Returns a struct \"oldResolutions\" with the current settings for screen "
                    "resolution, refresh rate and pixel depth. Optionally sets new settings "
                    "for screen resolution \"newwidth\" x \"newheight\", refresh rate \"newHz\" "
                    "and framebuffer pixel depth \"newPixelSize\". Providing invalid or incompatible "
                    "settings will raise an error. Especially the color depth \"newPixelSize\" should "
                    "usually not be set to anything else than its default of 32 bpp or 24 bpp, iow. it "
                    "should be best left alone. "
                    "Other settings can impair alpha-blending on some systems, a setting of 16 bpp "
                    "will disable alpha-blending and create drastically reduced color resolution of "
                    " 5 bits per color channel. "
                    "A setting of 8 bpp is not supported at all on MacOS/X and will create artifacts "
                    "on all other systems. Use a size of 32 bpp even for clut animation. This function "
                    "may not work on all MS-Windows setups, your mileage may vary...\n"
                    "On Linux/X11 with the classic X-Server, the function only switches display settings "
                    "in the conventional sense on a single-display setup, or on a multi-X-screen setup with "
                    "one display per X-Screen. On a multi-display-per-X-Screen setup, this function only "
                    "changes the total size of the framebuffer, ie., 'newwidth' and 'newheight', the other "
                    "parameters are silently ignored. You can create X-Screens and fullscreen onscreen "
                    "windows that are bigger than what is actually displayable on your attached displays. "
                    "This allows for some interesting special tricks if you know what you are doing, but can "
                    "also cause quite a bit of confusion. We recommend using this function on setups with "
                    "a single display or a single display for each X-Screen (as listed by Screen('Screens'). "
                    "On Linux/X11 setups with multiple attached displays per X-Screen, the video settings "
                    "of each individual display, e.g., resolution, video refresh rate, panning, can be queried "
                    "and changed via the Screen('ConfigureDisplay') function instead. This allows for more "
                    "flexibility, e.g., you can have a framebuffer bigger than the combined resolution of all "
                    "displays and only show a fraction of it. You can change the relative position of all "
                    "physical displays, configure \"mirror modes\", \"side by side\", or \"on top of each other\" "
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
    int oldhead;

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
    PsychAllocOutStructArray(1, FALSE, -1, 4, FieldNames, &oldResStructArray);

    // Query and return resolution:
    PsychGetScreenPixelSize(screenNumber, &pnewWidth, &pnewHeight);
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

        // So specialMode flag & 2 means on Linux to not change
        // the actual output resolution of a screen, but only to
        // resize the screen, ie., its framebuffer, without changing
        // scanout resolution / viewport etc. This would happen anyway
        // on multi-display setups, but to make sure it can also be
        // enforced on single-display setup, we use specialMode & 2:
        if ((PSYCH_SYSTEM == PSYCH_LINUX) && (specialMode & 2)) {
            // We fool the PsychSetScreenSettings() routine into thinking
            // this is at least a dual-display setup by temporarilly
            // faking a 2nd head, then doing the screen resize, then
            // restoring the true multi-display configuration:
            oldhead = PsychScreenToHead(screenNumber, 1);
            PsychSetScreenToHead(screenNumber, 0, 1);
            rc = PsychSetScreenSettings(TRUE, &screenSettings);
            PsychSetScreenToHead(screenNumber, oldhead, 1);
        }
        else {
            rc = PsychSetScreenSettings(TRUE, &screenSettings);
        }
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
