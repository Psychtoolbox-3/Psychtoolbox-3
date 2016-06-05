/*
    PsychToolbox3/Source/OSX/Screen/PsychScreenGlue.c

    PLATFORMS:

        This is the OS X Core Graphics version only.

    AUTHORS:

        Allen Ingling       awi     Allen.Ingling@nyu.edu
        Mario Kleiner       mk      mario.kleiner.de@gmail.com

    HISTORY:

        10/20/02        awi         Wrote it mostly by adding in SDL-specific refugeess (from an experimental SDL Psychtoolbox) from header and source files.
        11/16/04        awi         added  PsychGetGlobalScreenRect.  Enhanced DESCRIPTION text.
        04/29/05        mk          Return id of primary display for displays in hardware mirroring sets.
        12/29/06        mk          Implement query code for DAC output resolution on OS-X, finally...

    DESCRIPTION:

        Functions in this file comprise an abstraction layer for probing and controlling screen state.

        Each C function which implements a particular Screen subcommand should be platform neutral.  For example, the source to SCREENPixelSizes()
        should be platform-neutral, despite that the calls in OS X and Windows to detect available pixel sizes are
        different.  The platform specificity is abstracted out in C files which end it "Glue", for example PsychScreenGlue, PsychWindowGlue,
        PsychWindowTextClue.

        In addition to glue functions for windows and screen there are functions which implement shared functionality between between Screen commands,
        such as ScreenTypes.c and WindowBank.c.

    NOTES:

*/

#include "Screen.h"

// Include IOKIT support for connecting to the kernel level support driver:
#include <IOKit/IOKitLib.h>
// Include shared data structures and definitions between kernel driver and us:
#include "PsychUserKernelShared.h"
// Include specifications of the GPU registers:
#include "PsychGraphicsCardRegisterSpecs.h"
// Include for mouse cursor control via Cocoa:
#include "PsychCocoaGlue.h"

// Defined in PsychGraphicsHardwareHALSupport.c, but accessed and initialized here:
extern unsigned int crtcoff[kPsychMaxPossibleCrtcs];

// file local variables:
unsigned int  activeGPU = 0;
unsigned int  fDeviceType[kPsychMaxPossibleDisplays];
unsigned int  fCardType[kPsychMaxPossibleDisplays];
unsigned int  fPCIDeviceId[kPsychMaxPossibleDisplays];
unsigned int  fNumDisplayHeads[kPsychMaxPossibleDisplays];

// Maybe use NULLs in the settings arrays to mark entries invalid instead of using psych_bool flags in a different array.
static psych_bool           displayLockSettingsFlags[kPsychMaxPossibleDisplays];
static CGDisplayModeRef     displayOriginalCGSettings[kPsychMaxPossibleDisplays];       //these track the original video state before the Psychtoolbox changed it.
static psych_bool           displayOriginalCGSettingsValid[kPsychMaxPossibleDisplays];
static CGDisplayModeRef     displayOverlayedCGSettings[kPsychMaxPossibleDisplays];      //these track settings overlayed with 'Resolutions'.
static psych_bool           displayOverlayedCGSettingsValid[kPsychMaxPossibleDisplays];
static CGDisplayCount       numDisplays, numPhysicalDisplays;
static CGDirectDisplayID    displayCGIDs[kPsychMaxPossibleDisplays];
static CGDirectDisplayID    displayOnlineCGIDs[kPsychMaxPossibleDisplays];

// List of service connect handles for connecting to the kernel support driver (if any):
static int                  numKernelDrivers = 0;
static io_connect_t         displayConnectHandles[kPsychMaxPossibleDisplays];

//file local functions
void InitCGDisplayIDList(void);
void PsychLockScreenSettings(int screenNumber);
void PsychUnlockScreenSettings(int screenNumber);
psych_bool PsychCheckScreenSettingsLock(int screenNumber);
psych_bool PsychGetCGModeFromVideoSetting(CGDisplayModeRef *cgMode, PsychScreenSettingsType *setting);
void InitPsychtoolboxKernelDriverInterface(void);
kern_return_t PsychOSKDDispatchCommand(io_connect_t connect, const PsychKDCommandStruct* inStruct, PsychKDCommandStruct* outStruct, unsigned int* status);
io_connect_t PsychOSCheckKDAvailable(int screenId, unsigned int * status);
int PsychOSKDGetBeamposition(int screenId);
void PsychDisplayReconfigurationCallBack (CGDirectDisplayID display, CGDisplayChangeSummaryFlags flags, void *userInfo);
static void PsychOSKDGetGPUInfo(io_connect_t connect, int slot);
unsigned int PsychOSKDGetRevision(io_connect_t connect);

kern_return_t CallKDSimpleMethod(io_connect_t connect, uint32_t index)
{
    kern_return_t result;
    uint32_t outputCnt = 0;

    // IOConnectCallScalarMethod replaces IOConnectMethodScalarIScalarO
    // in OS/X 10.5 and later, and it has an incompatible interface in general,
    // but luckily we only use this method in its most simple form with 0 inputs
    // and outputs. For this special case we have a very simple 1-to-1 mapping:
    result = IOConnectCallScalarMethod((mach_port_t) connect, index, NULL, 0, NULL, &outputCnt);

    return(result);
}

//Initialization functions
void InitializePsychDisplayGlue(void)
{
    int i;

    //init the display settings flags.
    for(i=0;i<kPsychMaxPossibleDisplays;i++){
        displayLockSettingsFlags[i]=FALSE;
        displayOriginalCGSettingsValid[i]=FALSE;
        displayOverlayedCGSettingsValid[i]=FALSE;
        displayConnectHandles[i]=0;
    }

    #ifdef PTBOCTAVE3MEX
        // Restrict the latest idiotic hack to Octave on OSX 10.11+
        int major, minor, patchlevel;
        PsychCocoaGetOSXVersion(&major, &minor, &patchlevel);
        if ((major > 10) || (minor >= 11)) {
            // Another tribute to the most idiotic OS in existence: Redirect the stderr
            // stream, so OSX 10.11.0 El Capitans broken logger can't flood us with
            // pointless warning messages anymore. These unsolicited and unexpected
            // spam messages caused Octave's GUI to lock up completely after the first
            // run of a PTB script. The softer method of installing an asl filter didn't
            // work, so we do it the brute force way and just pray this will not cause
            // larger problems and side effects somewhere else (haha, hope against hope,
            // this would be the first time a dirty hack wouldn't bite us when dealing with
            // Apples crappy products).
            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: Redirecting stderr to work around broken OSX 10.11. This may have unpleasant side-effects.\n");
            freopen(((getenv("PSYCH_REDIRECT_STDERR_PATH")) ? getenv("PSYCH_REDIRECT_STDERR_PATH") : "/dev/null"), "a", stderr);
        }
    #endif

    // Init the list of Core Graphics display IDs.
    InitCGDisplayIDList();

    // Setup screenId -> display head mappings:
    PsychInitScreenToHeadMappings(PsychGetNumDisplays());

    // Register a display reconfiguration callback:
    // CGDisplayRegisterReconfigurationCallback(PsychDisplayReconfigurationCallBack, NULL);

    // Attach to kernel-level Psychtoolbox graphics card interface driver if possible
    // *and* allowed by settings, setup all relevant mappings:
    InitPsychtoolboxKernelDriverInterface();

    // Prevent OSX 10.9+ "AppNap" power saving and timer coalescing etc.:
    PsychCocoaPreventAppNap(TRUE);
}

void PsychCleanupDisplayGlue(void)
{
    int i;

    // Shutdown connection to kernel level driver, if any exists:
    PsychOSShutdownPsychtoolboxKernelDriverInterface();

    // Unregister our display reconfiguration callback:
    // CGDisplayRemoveReconfigurationCallback(PsychDisplayReconfigurationCallBack, NULL);

    // Release retained display mode objects:
    for(i = 0; i < kPsychMaxPossibleDisplays; i++){
        if (displayOriginalCGSettingsValid[i]) CGDisplayModeRelease(displayOriginalCGSettings[i]);
        if (displayOverlayedCGSettingsValid[i]) CGDisplayModeRelease(displayOverlayedCGSettings[i]);
        displayOriginalCGSettingsValid[i] = FALSE;
        displayOverlayedCGSettingsValid[i] = FALSE;
    }

    // Release font database:
    PsychFreeFontList();

    // Allow OSX 10.9+ "AppNap" power saving and timer coalescing etc.:
    PsychCocoaPreventAppNap(FALSE);
}

// PsychDisplayReconfigurationCallBack() respond to display reconfiguration events. This function is
// currently unused and not hooked up - Just left here for reference. It might have been involved in
// problems on OSX 10.9 and 10.10 in that releasing a display left the system hanging with a black
// screen for long periods of time. This possibly due to some broken OSX event handling interacting
// with this callback...
void PsychDisplayReconfigurationCallBack(CGDirectDisplayID display, CGDisplayChangeSummaryFlags flags, void *userInfo)
{
    (void) userInfo;

    // Provide feedback at verbosity level 4 or higher:
    if (PsychPrefStateGet_Verbosity() > 3) {
        if (flags & kCGDisplayBeginConfigurationFlag) printf("PTB-INFO: Display reconfiguration for display %p in progress...\n", display);
        if (flags & ~kCGDisplayBeginConfigurationFlag) printf("PTB-INFO: Reconfiguration for display %p finished: Flags = %i. Reenumerating all displays.\n", display, flags);
    }

    // Display reconfiguration finished?
    if (flags & ~kCGDisplayBeginConfigurationFlag) {
        // Yes: Reenumerate all our displays.
        InitCGDisplayIDList();
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Display reenumeration done.\n");
    }

    return;
}

void InitCGDisplayIDList(void)
{
    CGDisplayErr error;

    // Safe-Init to zero available displays:
    numDisplays = numPhysicalDisplays = 0;

    error = CGGetActiveDisplayList(kPsychMaxPossibleDisplays, displayCGIDs, &numDisplays);
    if (error) {
        printf("PTB-CRITICAL: CGGetActiveDisplayList failed to enumerate displays! Screen() will be mostly dysfunctional!\n");
        numDisplays = numPhysicalDisplays = 0;
        return;
    }

    // Also enumerate physical displays:
    error = CGGetOnlineDisplayList(kPsychMaxPossibleDisplays, displayOnlineCGIDs, &numPhysicalDisplays);
    if (error) {
        printf("PTB-CRITICAL: CGGetOnlineDisplayList failed to enumerate displays! Screen() will be mostly dysfunctional!\n");
        numDisplays = numPhysicalDisplays = 0;
        return;
    }
}

void PsychGetCGDisplayIDFromScreenNumber(CGDirectDisplayID *displayID, int screenNumber)
{
    if(screenNumber>= (int) numDisplays) PsychErrorExit(PsychError_invalidScumber);

    if (screenNumber < 0) {
        // Special case: Physical displays handle: Put back into positive range and
        // correct for 1-based external indexing:
        screenNumber = (-1 * screenNumber) - 1;
        if (screenNumber >= (int) numPhysicalDisplays) PsychErrorExitMsg(PsychError_user, "Invalid physical screenNumber provided! Higher than number of connected physical displays!");

        // Valid range: Map it:
        *displayID=displayOnlineCGIDs[screenNumber];

        return;
    }

    // Standard case: Logical displays:

    // MK: We return the id of the primary display of the hardware-mirror set to which
    // the display for 'screenNumber' belongs to. This will be the same display on
    // single display setups. On dual-display setups, it will return the ID of the
    // display we are really syncing in Screen('Flip'). This is important for querying
    // the rasterbeam position of the correct display in such setups.
    //
    // I'm not sure if this is the best place for performing this lookup, but
    // at least it should be safe to do here...
    *displayID=CGDisplayPrimaryDisplay(displayCGIDs[screenNumber]);

    // Protection against Apples stupidity... - our our's if we made false assumptions ;-)
    if (CGDisplayUnitNumber(*displayID)!=CGDisplayUnitNumber(displayCGIDs[screenNumber])) {
        mexPrintf("PTB-DEBUG : ACTIVE DISPLAY <-> PRIMARY DISPLAY MISMATCH FOR SCREEN %i!!!!\n", screenNumber);
    }

    return;
}

/*  About locking display settings:

    SCREENOpenWindow and SCREENOpenOffscreenWindow  set the lock when opening  windows and
    SCREENCloseWindow unsets upon the close of the last of a screen's windows. PsychSetVideoSettings checks for a lock
    before changing the settings.  Anything (SCREENOpenWindow or SCREENResolutions) which attemps to change
    the display settings should report that attempts to change a dipslay's settings are not allowed when its windows are open.

    PsychSetVideoSettings() may be called by either SCREENOpenWindow or by Resolutions().  If called by Resolutions it both sets the video settings
    and caches the video settings so that subsequent calls to OpenWindow can use the cached mode regardless of whether interceding calls to OpenWindow
    have changed the display settings or reverted to the virgin display settings by closing.  SCREENOpenWindow should thus invoke the cached mode
    settings if they have been specified and not current actual display settings.

*/
void PsychLockScreenSettings(int screenNumber)
{
    displayLockSettingsFlags[screenNumber]=TRUE;
}

void PsychUnlockScreenSettings(int screenNumber)
{
    displayLockSettingsFlags[screenNumber]=FALSE;
}

psych_bool PsychCheckScreenSettingsLock(int screenNumber)
{
    return(displayLockSettingsFlags[screenNumber]);
}

/* Because capture and lock will always be used in conjuction, capture calls lock, and SCREENOpenWindow must only call Capture and Release */
void PsychCaptureScreen(int screenNumber)
{
    CGDisplayErr  error;

    if (screenNumber >= numDisplays) PsychErrorExit(PsychError_invalidScumber);

    error = CGDisplayCapture(displayCGIDs[screenNumber]);
    if (error) PsychErrorExitMsg(PsychError_internal, "Unable to capture display");

    PsychLockScreenSettings(screenNumber);

    return;
}

/*
 *  PsychReleaseScreen()
 */
void PsychReleaseScreen(int screenNumber)
{
    CGDisplayErr  error;

    if (screenNumber >= numDisplays) PsychErrorExit(PsychError_invalidScumber);

    error = CGDisplayRelease(displayCGIDs[screenNumber]);
    if (error) PsychErrorExitMsg(PsychError_internal, "Unable to release display");

    PsychUnlockScreenSettings(screenNumber);

    // Try to restore keyboard input focus to whatever window had focus before
    // the CGDisplayCapture()/CGDisplayRelease(). Turns out to be a bit unreliable,
    // and of limited use when it works, but anyway...
    PsychCocoaSetUserFocusWindow(NULL);

    return;
}

psych_bool PsychIsScreenCaptured(screenNumber)
{
    return(PsychCheckScreenSettingsLock(screenNumber));
}

//Read display parameters.
/*
    PsychGetNumDisplays()
    Get the number of video displays connected to the system.
*/
int PsychGetNumDisplays(void)
{
    return((int)numDisplays);
}

/* This is only defined on OS/X for now: */
int PsychGetNumPhysicalDisplays(void)
{
    return((int) numPhysicalDisplays);
}

/* Returns actual net color bits per pixel (excluding alpha or padding): */
static int getDisplayBitsPerPixel(CGDisplayModeRef mode)
{
    int bpp = 0;
    CFStringRef n = CGDisplayModeCopyPixelEncoding(mode);

    if (CFStringCompare(n, CFSTR(kIO30BitDirectPixels), kCFCompareCaseInsensitive) == kCFCompareEqualTo) {
        bpp = 30;
    }
    else if (CFStringCompare(n, CFSTR(IO32BitDirectPixels), kCFCompareCaseInsensitive) == kCFCompareEqualTo) {
        bpp = 24;
    }
    else if (CFStringCompare(n, CFSTR(IO16BitDirectPixels), kCFCompareCaseInsensitive) == kCFCompareEqualTo) {
        bpp = 15;
    }
    else if (CFStringCompare(n, CFSTR(IO8BitIndexedPixels), kCFCompareCaseInsensitive) == kCFCompareEqualTo) {
        bpp = 8;
    }

    CFRelease(n);

    return(bpp);
}

void PsychGetScreenDepths(int screenNumber, PsychDepthType *depths)
{
    CGDisplayModeRef currentMode, tempMode;
    CFArrayRef modeList;
    int i, numPossibleModes;
    long currentWidth, currentHeight, tempWidth, tempHeight, tempDepth;
    double currentFrequency, tempFrequency;

    if (screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);

    // Get the current display mode.  We will want to match against width and hz when looking for available depths.
    currentMode = CGDisplayCopyDisplayMode(displayCGIDs[screenNumber]);
    currentWidth = (long) CGDisplayModeGetWidth(currentMode);
    currentHeight = (long) CGDisplayModeGetHeight(currentMode);
    currentFrequency = CGDisplayModeGetRefreshRate(currentMode);
    CGDisplayModeRelease(currentMode);

    // Build options dictionary to make sure we also get HiDPI / Retina scaled modes:
    i = 1;
    CFNumberRef number = CFNumberCreate(kCFAllocatorDefault, kCFNumberIntType, &i);
    CFStringRef key = kCGDisplayShowDuplicateLowResolutionModes;
    CFDictionaryRef options = CFDictionaryCreate(kCFAllocatorDefault, (const void**) &key, (const void**) &number, 1, NULL, NULL);
    CFRelease(number);

    // Get a list of available modes for the specified display, including HiDPI/Retina
    // modes, as requested by the options dictionary:
    modeList = CGDisplayCopyAllDisplayModes(displayCGIDs[screenNumber], options);
    CFRelease(options);

    numPossibleModes = CFArrayGetCount(modeList);
    for (i = 0; i < numPossibleModes; i++) {
        tempMode = (CGDisplayModeRef) CFArrayGetValueAtIndex(modeList,i);
        tempWidth = (long) CGDisplayModeGetWidth(tempMode);
        tempHeight = (long) CGDisplayModeGetHeight(tempMode);
        tempFrequency = CGDisplayModeGetRefreshRate(tempMode);
        tempDepth = getDisplayBitsPerPixel(tempMode);

        if (currentWidth == tempWidth && currentHeight == tempHeight && currentFrequency == tempFrequency) {
            PsychAddValueToDepthStruct((int) tempDepth, depths);
        }

        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: PsychGetScreenDepths(): mode %i : w x h = %i x %i, fps = %f, depth = %i\n", i, tempWidth, tempHeight, tempFrequency, tempDepth);
    }

    CFRelease(modeList);

    // At least one match?
    if (PsychGetNumDepthsFromStruct(depths) < 1) {
        // Yes, this should not ever happen on a sane operating system, but this is
        // OSX, so it does. Observed on a 2010 MacBookPro with OSX 10.7.5 on a external
        // panel. Output a warning and fake entries for the most common pixel sizes:
        PsychAddValueToDepthStruct(16, depths);
        PsychAddValueToDepthStruct(32, depths);
        if (PsychPrefStateGet_Verbosity() > 1) {
            printf("PTB-WARNING: Broken MacOS/X detected. It misreports (== omits some) available video modes and thereby returns empty display depths due to matching failure.\n");
            printf("PTB-WARNING: Will try to workaround this by creating a fake list of available display depths of 16 bpp and 32 bpp. Expect potential trouble further on...\n");
        }
    }
}

/*   PsychGetAllSupportedScreenSettings()
 *
 *     Queries the display system for a list of all supported display modes, ie. all valid combinations
 *     of resolution, pixeldepth and refresh rate. Allocates temporary arrays for storage of this list
 *     and returns it to the calling routine. This function is basically only used by Screen('Resolutions').
 */
int PsychGetAllSupportedScreenSettings(int screenNumber, int outputId, long** widths, long** heights, long** hz, long** bpp)
{
    CGDisplayModeRef tempMode;
    CFArrayRef modeList;
    int i, numPossibleModes;
    long tempWidth, tempHeight, tempDepth;
    double tempFrequency;

    if (screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);

    // Build options dictionary to make sure we also get HiDPI / Retina scaled modes:
    i = 1;
    CFNumberRef number = CFNumberCreate(kCFAllocatorDefault, kCFNumberIntType, &i);
    CFStringRef key = kCGDisplayShowDuplicateLowResolutionModes;
    CFDictionaryRef options = CFDictionaryCreate(kCFAllocatorDefault, (const void**) &key, (const void**) &number, 1, NULL, NULL);
    CFRelease(number);

    // Get a list of available modes for the specified display, including HiDPI/Retina
    // modes, as requested by the options dictionary:
    modeList = CGDisplayCopyAllDisplayModes(displayCGIDs[screenNumber], options);
    CFRelease(options);
    numPossibleModes = CFArrayGetCount(modeList);

    // Allocate output arrays: These will get auto-released at exit from Screen():
    *widths = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *heights = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *hz = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *bpp = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));

    // Fetch modes and store into arrays:
    for(i = 0; i < numPossibleModes; i++) {
        tempWidth = tempHeight = tempFrequency = tempDepth = 0;

        tempMode = (CGDisplayModeRef) CFArrayGetValueAtIndex(modeList,i);
        tempWidth = (long) CGDisplayModeGetWidth(tempMode);
        tempHeight = (long) CGDisplayModeGetHeight(tempMode);
        tempFrequency = CGDisplayModeGetRefreshRate(tempMode);
        tempDepth = getDisplayBitsPerPixel(tempMode);

        (*widths)[i] = tempWidth;
        (*heights)[i] = tempHeight;
        (*hz)[i] = (long) (tempFrequency + 0.5);
        (*bpp)[i] = tempDepth;
    }

    CFRelease(modeList);

    return(numPossibleModes);
}

/* Only returned a cgMode which needs to be CGDisplayModeRelease()'ed on success, ie., if return true */
psych_bool PsychGetCGModeFromVideoSetting(CGDisplayModeRef *cgMode, PsychScreenSettingsType *setting)
{
    CGDisplayModeRef tempMode;
    CFArrayRef modeList;
    int i, numPossibleModes;
    long width, height, depth, tempWidth, tempHeight, tempDepth;
    double tempFrequency, frameRate;

    if(setting->screenNumber>=numDisplays) {
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetCGModeFromVideoSetting() is out of range.");
    }

    width = (long) PsychGetWidthFromRect(setting->rect);
    height = (long) PsychGetHeightFromRect(setting->rect);
    depth = (long) PsychGetValueFromDepthStruct(0,&(setting->depth));
    frameRate = (double) setting->nominalFrameRate;

    // Build options dictionary to make sure we also get HiDPI / Retina scaled modes:
    i = 1;
    CFNumberRef number = CFNumberCreate(kCFAllocatorDefault, kCFNumberIntType, &i);
    CFStringRef key = kCGDisplayShowDuplicateLowResolutionModes;
    CFDictionaryRef options = CFDictionaryCreate(kCFAllocatorDefault, (const void**) &key, (const void**) &number, 1, NULL, NULL);
    CFRelease(number);

    // Get a list of available modes for the specified display, including HiDPI/Retina
    // modes, as requested by the options dictionary:
    modeList = CGDisplayCopyAllDisplayModes(displayCGIDs[setting->screenNumber], options);
    CFRelease(options);
    numPossibleModes = CFArrayGetCount(modeList);

    // Fetch modes and store into arrays:
    for(i = 0; i < numPossibleModes; i++) {
        *cgMode = (CGDisplayModeRef) CFArrayGetValueAtIndex(modeList,i);

        tempWidth = (long) CGDisplayModeGetPixelWidth(*cgMode);
        tempHeight = (long) CGDisplayModeGetPixelHeight(*cgMode);
        tempFrequency = CGDisplayModeGetRefreshRate(*cgMode);
        tempDepth = getDisplayBitsPerPixel(*cgMode);

        // Match? Be lenient for frameRate, as OSX operates with fractional doubles, but our current api only with integral values.
        if ((width == tempWidth) && (height == tempHeight) && (fabs(frameRate - tempFrequency) <= 1) && (depth == tempDepth)) {
            CGDisplayModeRetain(*cgMode);
            CFRelease(modeList);
            return(TRUE);
        }
    }

    // Failed.
    CFRelease(modeList);
    *cgMode = NULL;

    return(FALSE);
}

/*
    PsychCheckVideoSettings()

    Check all available video display modes for the specified screen number and return true if the
    settings are valid and false otherwise.
*/
psych_bool PsychCheckVideoSettings(PsychScreenSettingsType *setting)
{
    CGDisplayModeRef cgMode;
    if (PsychGetCGModeFromVideoSetting(&cgMode, setting)) {
        CGDisplayModeRelease(cgMode);
        return(TRUE);
    }

    return(FALSE);
}

/*
 *  PsychGetScreenDepth()
 *  The caller must allocate and initialize the depth struct.
 */
void PsychGetScreenDepth(int screenNumber, PsychDepthType *depth)
{
    if (screenNumber >= numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");
    CGDisplayModeRef currentMode = CGDisplayCopyDisplayMode(displayCGIDs[screenNumber]);
    PsychAddValueToDepthStruct((int) getDisplayBitsPerPixel(currentMode), depth);
    CGDisplayModeRelease(currentMode);
}

int PsychGetScreenDepthValue(int screenNumber)
{
    PsychDepthType    depthStruct;

    PsychInitDepthStruct(&depthStruct);
    PsychGetScreenDepth(screenNumber, &depthStruct);
    return(PsychGetValueFromDepthStruct(0,&depthStruct));
}

float PsychGetNominalFramerate(int screenNumber)
{
    double currentFrequency;

    if (PsychPrefStateGet_ConserveVRAM() & kPsychIgnoreNominalFramerate) return(0);

    if (screenNumber >= numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");

    CGDisplayModeRef currentMode = CGDisplayCopyDisplayMode(displayCGIDs[screenNumber]);
    currentFrequency = CGDisplayModeGetRefreshRate(currentMode);
    CGDisplayModeRelease(currentMode);
    return((float) currentFrequency);
}

// Report video mode size in pixels, as opposed to logical units (points) like the
// cross-platform PsychGetScreenSize(). Need this to handle Apples Retina / HiDPI:
void PsychGetScreenPixelSize(int screenNumber, long *width, long *height)
{
    if (screenNumber >= numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");

    CGDisplayModeRef currentMode = CGDisplayCopyDisplayMode(displayCGIDs[screenNumber]);
    *width = (long) CGDisplayModeGetPixelWidth(currentMode);
    *height = (long) CGDisplayModeGetPixelHeight(currentMode);
    CGDisplayModeRelease(currentMode);
}

void PsychGetScreenSize(int screenNumber, long *width, long *height)
{
    if (screenNumber >= numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");

    CGDisplayModeRef currentMode = CGDisplayCopyDisplayMode(displayCGIDs[screenNumber]);
    *width = (long) CGDisplayModeGetWidth(currentMode);
    *height = (long) CGDisplayModeGetHeight(currentMode);
    CGDisplayModeRelease(currentMode);
}

/* Returns the physical display size as reported by OS-X: */
void PsychGetDisplaySize(int screenNumber, int *width, int *height)
{
    CGSize physSize;
    if(screenNumber>=numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetDisplaySize() is out of range");
    physSize = CGDisplayScreenSize(displayCGIDs[screenNumber]);
    *width = (int) physSize.width;
    *height = (int) physSize.height;
}

void PsychGetGlobalScreenRect(int screenNumber, double *rect)
{
    CGDirectDisplayID    displayID;
    CGRect                cgRect;
    double                rLeft, rRight, rTop, rBottom;

    if(screenNumber>=numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetGlobalScreenRect() is out of range");

    displayID=displayCGIDs[screenNumber];
    cgRect=CGDisplayBounds(displayID);
    rLeft=cgRect.origin.x;
    rTop=cgRect.origin.y;
    rRight=cgRect.origin.x + cgRect.size.width;
    rBottom=cgRect.origin.y + cgRect.size.height;
    PsychMakeRect(rect, rLeft, rTop, rRight, rBottom);
}

void PsychGetScreenRect(int screenNumber, double *rect)
{
    long width, height;

    PsychGetScreenSize(screenNumber, &width, &height);
    rect[kPsychLeft]=0;
    rect[kPsychTop]=0;
    rect[kPsychRight]=(int)width;
    rect[kPsychBottom]=(int)height;
}

PsychColorModeType PsychGetScreenMode(int screenNumber)
{
    PsychDepthType depth;

    PsychInitDepthStruct(&depth);
    PsychGetScreenDepth(screenNumber, &depth);
    return(PsychGetColorModeFromDepthStruct(&depth));
}

/*
    PsychGetDacBitsFromDisplay()

    Return output resolution of video DAC in bits per color component.
    We return a safe default of 8 bpc if we can't query the real value.
*/
int PsychGetDacBitsFromDisplay(int screenNumber)
{
    CGDirectDisplayID    displayID;
    CFMutableDictionaryRef properties;
    CFNumberRef cfGammaWidth;
    SInt32 dacbits;
    io_service_t displayService;
    kern_return_t kr;

    // Retrieve display handle for screen:
    PsychGetCGDisplayIDFromScreenNumber(&displayID, screenNumber);

    // Retrieve low-level IOKit service port for this display:
    displayService = CGDisplayIOServicePort(displayID);

    // Obtain the properties from that service
    kr = IORegistryEntryCreateCFProperties(displayService, &properties, NULL, 0);
    if((kr == kIOReturnSuccess) && ((cfGammaWidth = (CFNumberRef) CFDictionaryGetValue(properties, CFSTR(kIOFBGammaWidthKey)))!=NULL)) {
        CFNumberGetValue(cfGammaWidth, kCFNumberSInt32Type, &dacbits);
        CFRelease(properties);
        return((int) dacbits);
    }
    else {
        // Failed! Return safe 8 bits...
        CFRelease(properties);
        if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Failed to query resolution of video DAC for screen %i! Will return safe default of 8 bits.\n", screenNumber);
        return(8);
    }
}

/*
    PsychGetVideoSettings()

    Fills a structure describing the screen settings such as x, y, depth, frequency, etc.
    Consider inverting the calling sequence so that this function is at the bottom of call hierarchy.
*/
void PsychGetScreenSettings(int screenNumber, PsychScreenSettingsType *settings)
{
    settings->screenNumber=screenNumber;
    PsychGetScreenRect(screenNumber, settings->rect);
    PsychInitDepthStruct(&(settings->depth));
    PsychGetScreenDepth(screenNumber, &(settings->depth));
    settings->mode=PsychGetColorModeFromDepthStruct(&(settings->depth));
    settings->nominalFrameRate=PsychGetNominalFramerate(screenNumber);
    // settings->dacbits=PsychGetDacBitsFromDisplay(screenNumber);
}

//Set display parameters

/*
    PsychSetScreenSettings()

    Accept a PsychScreenSettingsType structure holding a video mode and set the display mode accordingly.

    If we can not change the display settings because of a lock (set by open window or close window) then return false.

    SCREENOpenWindow should capture the display before it sets the video mode.  If it doesn't, then PsychSetVideoSettings will
    detect that and exit with an error.  SCREENClose should uncapture the display.

    The duties of SCREENOpenWindow are:
    -Lock the screen which serves the purpose of preventing changes in video setting with open Windows.
    -Capture the display which gives the application synchronous control of display parameters.

    TO DO: for 8-bit palletized mode there is probably more work to do.

*/
psych_bool PsychSetScreenSettings(psych_bool cacheSettings, PsychScreenSettingsType *settings)
{
    CGDisplayModeRef    cgMode;
    psych_bool          isValid, isCaptured;
    CGDisplayErr        error;
    float               *redlut, *greenlut, *bluelut;
    int                 numlutslots;

    if(settings->screenNumber >= numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychSetScreenSettings() is out of range");

    //store the original display mode if this is the first time we have called this function.  The psychtoolbox will disregard changes in
    //the screen state made through the control panel after the Psychtoolbox was launched. That is, OpenWindow will by default continue to
    //open windows with finder settings which were in place at the first call of OpenWindow.  That's not intuitive, but not much of a problem
    //either.
    if(!displayOriginalCGSettingsValid[settings->screenNumber]) {
        displayOriginalCGSettings[settings->screenNumber] = CGDisplayCopyDisplayMode(displayCGIDs[settings->screenNumber]);
        displayOriginalCGSettingsValid[settings->screenNumber] = TRUE;
    }

    //Find core graphics video settings which correspond to settings as specified withing by an abstracted psychsettings structure.
    isValid = PsychGetCGModeFromVideoSetting(&cgMode, settings);
    if (!isValid) PsychErrorExitMsg(PsychError_internal, "Attempt to set invalid video settings");

    //If the caller passed cache settings (then it is SCREENResolutions) and we should cache the current video mode settings for this display.  These
    //are cached in the form of CoreGraphics settings and not Psychtoolbox video settings.  The only caller which should pass a set cache flag is
    //SCREENResolutions
    if (cacheSettings) {
        displayOverlayedCGSettings[settings->screenNumber] = cgMode;
        displayOverlayedCGSettingsValid[settings->screenNumber] = TRUE;
    }
    else CGDisplayModeRelease(cgMode);

    //Check to make sure that this display is captured, which OpenWindow should have done.  If it has not been done, then exit with an error.
    isCaptured = PsychIsScreenCaptured(settings->screenNumber);
    if(!isCaptured) PsychErrorExitMsg(PsychError_internal, "Attempt to change video settings without capturing the display");

    // Readback a backup copy of the current gamma table:
    PsychReadNormalizedGammaTable(settings->screenNumber, -1, &numlutslots, &redlut, &greenlut, &bluelut);

    //Change the display mode.
    CGDisplayConfigRef configRef;
    error = CGBeginDisplayConfiguration(&configRef);
    if (kCGErrorSuccess == error) error = CGConfigureDisplayWithDisplayMode(configRef, displayCGIDs[settings->screenNumber], cgMode, NULL);
    if (kCGErrorSuccess == error) {
        error = CGCompleteDisplayConfiguration(configRef, kCGConfigureForAppOnly);
    }
    else {
        CGCancelDisplayConfiguration(configRef);
    }

    // And now we rest for 3 seconds after video mode switch to honor all the mentally troubled
    // people working at Apple. At least OSX 10.9 will restore the graphics cards gamma tables to
    // the user session default (or whatever) within 2 seconds after a video mode switch,
    // thereby silently undoing the effect of any of our gamma table setup operations, should they
    // happen within a 2 seconds interval to a modeset operation. To protect against this, we will
    // pause execution for 3 seconds.
    PsychYieldIntervalSeconds(3.0);

    // Restore pre-modeswitch gamma table to undo the brain-damage of the OSX modeswitch implementation:
    PsychLoadNormalizedGammaTable(settings->screenNumber, -1, numlutslots, redlut, greenlut, bluelut);

    return((kCGErrorSuccess == error) ? TRUE : FALSE);
}

/*
    PsychRestoreScreenSettings()

    Restores video settings to the state set by the finder.  Returns TRUE if the settings can be restored or false if they
    can not be restored because a lock is in effect, which would mean that there are still open windows.

*/
psych_bool PsychRestoreScreenSettings(int screenNumber)
{
    psych_bool          isCaptured;
    CGDisplayErr        error;

    if (screenNumber>=numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychRestoreScreenSettings() is out of range");

    //Check to make sure that the original graphics settings were cached. If not, it means that the settings were never changed, so we can just
    //return true.
    if (!displayOriginalCGSettingsValid[screenNumber]) return(true);

    //Check to make sure that this display is captured, which OpenWindow should have done.  If it has not been done, then exit with an error.
    isCaptured = PsychIsScreenCaptured(screenNumber);
    if (!isCaptured) printf("PTB-ERROR: Attempt to change video settings without capturing the display on screen %i.\n", screenNumber);

    // Change the display mode.
    CGDisplayConfigRef configRef;
    error = CGBeginDisplayConfiguration(&configRef);
    if (kCGErrorSuccess == error) error = CGConfigureDisplayWithDisplayMode(configRef, displayCGIDs[screenNumber], displayOriginalCGSettings[screenNumber], NULL);
    if (kCGErrorSuccess == error) {
        error = CGCompleteDisplayConfiguration(configRef, kCGConfigureForAppOnly);
    }
    else {
        CGCancelDisplayConfiguration(configRef);
    }

    if (error) printf("PTB-ERROR: Unable to restore original video mode setting on screen %i.\n", screenNumber);

    return(true);
}

void PsychHideCursor(int screenNumber, int deviceIdx)
{
    CGDisplayErr        error;
    CGDirectDisplayID   cgDisplayID;

    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);
    error=CGDisplayHideCursor(cgDisplayID);
    if (error)
        PsychErrorExit(PsychError_internal);
}

void PsychShowCursor(int screenNumber, int deviceIdx)
{
    CGDisplayErr        error;
    CGDirectDisplayID   cgDisplayID;

    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);
    error=CGDisplayShowCursor(cgDisplayID);
    if(error)
        PsychErrorExit(PsychError_internal);
}

void PsychPositionCursor(int screenNumber, int x, int y, int deviceIdx)
{
    CGDirectDisplayID   cgDisplayID;
    CGPoint             point;

    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);
    point.x=(float)x;
    point.y=(float)y;

    // Move cursor. These functions return non-zero error status on failure.
    //
    // Note: This "utterly senseless" call to CGAssociateMouseAndMouseCursorPosition(true) is here to
    // avoid suppression of delivery of mouse input events from connected mice to the system
    // for a suppression interval of 250 msecs after CGDisplayMoveCursorToPoint(), aka posting
    // some mouse positioning Quartz event. Of course this is in direct contradiction to Apple docs
    // which state that "CGDisplayMoveCursorToPoint() works "without generating or posting an event"",
    // but hey, Apple logic!
    //
    // The solution implemented here was found in this StackExchange posting:
    // http://stackoverflow.com/questions/10196603/using-cgeventsourcesetlocaleventssuppressioninterval-instead-of-the-deprecated
    //
    // A fine read to see the sad state of Apple docs - the stupid, it hurts!
    //
    // Without this hack, the mouse cursor would not respond to actual mouse movement by the user for
    // a time interval of 250 msecs after executing SetMouse()! Fixes PTB forum bug from message #19365.
    if (CGDisplayMoveCursorToPoint(cgDisplayID, point) || CGAssociateMouseAndMouseCursorPosition(true))
        PsychErrorExit(PsychError_internal);
}

/*
    PsychReadNormalizedGammaTable()

    TO DO: This should probably be changed so that the caller allocates the memory.
    TO DO: Adopt a naming convention which distinguishes between functions which allocate memory for return variables from those which do not.
            For example, PsychReadNormalizedGammaTable() vs. PsychGetNormalizedGammaTable().

*/
void PsychReadNormalizedGammaTable(int screenNumber, int outputId, int *numEntries, float **redTable, float **greenTable, float **blueTable)
{
    CGDirectDisplayID   cgDisplayID;
    static float        localRed[4096], localGreen[4096], localBlue[4096];
    CGDisplayErr        error;
    uint32_t            sampleCount, realLutCapacity;

    *redTable=localRed; *greenTable=localGreen; *blueTable=localBlue;
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);
    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: ReadNormalizedGammatable: screenid %i mapped to CGDisplayId %p.\n", screenNumber, cgDisplayID);

    // Another hack for the brain-damaged OSX. Tell the OS we can't take more than what it claims
    // the real LUT size is, to avoid it doing stupid interpolation tricks which seem to be buggy
    // in at least OSX 10.11.
    realLutCapacity = CGDisplayGammaTableCapacity(cgDisplayID);
    if (realLutCapacity > 4096) {
        if (PsychPrefStateGet_Verbosity() > 1)
            printf("PTB-WARNING: ReadNormalizedGammatable: OSX reported LUT size %i exceeds our current limit of 4096 slots, clamping! Please tell the PTB developers about this.\n",
                   realLutCapacity);
        realLutCapacity = 4096;
    }

    error=CGGetDisplayTransferByTable(cgDisplayID, realLutCapacity, *redTable, *greenTable, *blueTable, &sampleCount);
    *numEntries=(int)sampleCount;
    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: ReadNormalizedGammatable: numEntries = %i.\n", *numEntries);
}

unsigned int PsychLoadNormalizedGammaTable(int screenNumber, int outputId, int numEntries, float *redTable, float *greenTable, float *blueTable)
{
    CGDisplayErr error = 0;
    CGDirectDisplayID cgDisplayID;
    uint8_t byteLUT[256];
    int i;

    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);
    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: LoadNormalizedGammatable: screenid %i mapped to CGDisplayId %p.\n", screenNumber, cgDisplayID);

    // More than one row in table?
    if (numEntries > 1) {
        // Yes: This is the regular case. We upload a 0.0 - 1.0 encoded table with numEntries slots. The OS will
        // interpolate inbetween our slots if the number of required slots in the GPU doesn't match the numEntries
        // we provided:
        error = CGSetDisplayTransferByTable(cgDisplayID, numEntries, redTable, greenTable, blueTable);
        if (error && (PsychPrefStateGet_Verbosity() > 0)) printf("PTB-ERROR: Failed to update the gamma tables for screen %i in call to CGSetDisplayTransferByTable()!\n", screenNumber);
    }
    else {
        if (numEntries <= 0) {
            // No: Special case 0-Slot table. We shall upload an identity CLUT:
            for (i = 0; i < 256; i++) byteLUT[i] = (uint8_t) i;
            error = CGSetDisplayTransferByByteTable(cgDisplayID, 256, &byteLUT[0], &byteLUT[0], &byteLUT[0]);
            if (error && (PsychPrefStateGet_Verbosity() > 0)) printf("PTB-ERROR: Failed to upload identity gamma tables for screen %i in call to CGSetDisplayTransferByByteTable()!\n", screenNumber);
        }
        else {
            // No: Special case 1-slot table: Interpret red, green and blue value of single slot as min, max and gamma value for CGSetDisplayTransferByFormula()
            error = CGSetDisplayTransferByFormula(cgDisplayID, redTable[0], greenTable[0], blueTable[0], redTable[0], greenTable[0], blueTable[0], redTable[0], greenTable[0], blueTable[0]);
            if (error && (PsychPrefStateGet_Verbosity() > 0)) printf("PTB-ERROR: Failed to upload computed gamma tables for screen %i in call to CGSetDisplayTransferByFormula()!\n", screenNumber);
        }
    }

    // Return success status:
    return((error == 0) ? 1 : 0);
}

// Return true (non-zero) if a desktop compositor is likely active on screen 'screenNumber':
int PsychOSIsDWMEnabled(int screenNumber)
{
    // Only way to disable compositing on OSX is to capture the screen
    // for exclusive fullscreen use, so composition state is the negation
    // of capture state:
    return(!PsychIsScreenCaptured(screenNumber));
}

// PsychGetDisplayBeamPosition() contains the implementation of display beamposition queries.
int PsychGetDisplayBeamPosition(CGDirectDisplayID cgDisplayId, int screenNumber)
{
    // As of PTB 3.0.12, always use and return beamposition as measured by our own
    // PsychtoolboxKernelDriver implementation. Since at least OSX 10.7, our own
    // implementation is much more reliable, robust, precise and low-overhead.
    // Beamposition queries weren't supported by OSX natively for Intel and AMD gpu's
    // on IntelMacs ever, they were severely buggy on OSX 10.8+ for the remaining
    // NVidia cards where they sort of "worked", and the functionality has been
    // effectively removed starting with OSX 10.9, turning into a no-op on all gpus.
    return(PsychOSKDGetBeamposition(screenNumber));
}

// Try to attach to kernel level ptb support driver and setup everything, if it works:
void InitPsychtoolboxKernelDriverInterface(void)
{
    kern_return_t       kernResult;
    io_service_t        service;
    io_connect_t        connect;
    io_iterator_t       iterator;
    CFDictionaryRef     classToMatch;
    int                 i, revision;

    // Reset to zero open drivers to start with:
    numKernelDrivers = 0;

    // Select first instance (index 0) as active GPU/KernelDriver by default:
    activeGPU = 0;

    // Setup matching criterion to find our driver in the IORegistry device tree:
    classToMatch = IOServiceMatching(kMyDriversIOKitClassName);
    if (classToMatch == NULL) {
        printf("PTB-DEBUG: IOServiceMatching() for Psychtoolbox kernel support driver returned a NULL dictionary. Kernel driver support disabled.\n");
        return;
    }

    // This creates an io_iterator_t of all instances of our driver that exist in the I/O Registry. Each installed graphics card
    // will get its own instance of a driver. The iterator allows to iterate over all instances:
    kernResult = IOServiceGetMatchingServices(kIOMasterPortDefault, classToMatch, &iterator);
    if (kernResult != KERN_SUCCESS) {
        printf("PTB-DEBUG: IOServiceGetMatchingServices for Psychtoolbox kernel support driver returned 0x%08x. Kernel driver support disabled.\n", kernResult);
        return;
    }

    // In a polished final version we would want to handle the case where more than one gfx-card is attached.
    // The iterator would return multiple instances of our driver and we need to decide which one to connect to.
    // For now, we do not handle this case but instead just get the first item from the iterator.
    while ((service = IOIteratorNext(iterator)) != IO_OBJECT_NULL) {
        // Instantiate a connection to the user client.

        // This call will cause the user client to be instantiated. It returns an io_connect_t handle
        // that is used for all subsequent calls to the user client.
        connect = IO_OBJECT_NULL;
        kernResult = IOServiceOpen(service, mach_task_self(), 0, &connect);
        if ((kernResult != KERN_SUCCESS) || (connect == IO_OBJECT_NULL)) {
            printf("PTB-DEBUG: IOServiceOpen for driver instance %i returned 0x%08x. Not using this instance...\n", service, kernResult);
        }
        else {
            // This is an example of calling our user client's openUserClient method.
            kernResult = CallKDSimpleMethod(connect, kMyUserClientOpen);
            if (kernResult != KERN_SUCCESS) {
                // Release connection:
                IOServiceClose(connect);
                connect = IO_OBJECT_NULL;
                printf("PTB-DEBUG: CallKDSimpleMethod for driver instance %i returned 0x%08x. Kernel driver support disabled.\n", service, kernResult);
                if (kernResult == kIOReturnExclusiveAccess) printf("PTB-DEBUG: Please check if other applications (e.g., other open Matlab or Octave instances) use the driver already.\n");
            }
        }

        // Release the io_service_t now that we're done with it.
        IOObjectRelease(service);

        if (connect != IO_OBJECT_NULL) {
            // Final success!

            // Query driver revision: We disconnect and don't use the driver if it
            // doesn't provide the required minimum revision number for its API:
            revision = (int) PsychOSKDGetRevision(connect);
            if (revision < 0) {
                printf("PTB-ERROR: The currently loaded PsychtoolboxKernelDriver.kext is outdated!\n");
                printf("PTB-ERROR: Its revision number is %i, but we require a minimum revision of 0, better 1.\n", revision);
                printf("PTB-ERROR: Please uninstall the current driver and reinstall the latest one delivered\n");
                printf("PTB-ERROR: with your Psychtoolbox (see 'help PsychtoolboxKernelDriver').\n");
                printf("PTB-ERROR: Driver support disabled for now, special functions not available.\n");

                // Call shutdown method:
                kern_return_t kernResult = CallKDSimpleMethod(connect, kMyUserClientClose);
                if (kernResult == KERN_SUCCESS) {
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: CallKDSimpleMethod Closing was successfull.\n");
                }
                else {
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: CallKDSimpleMethod Closing failed with kernel return code 0x%08x.\n\n", kernResult);
                }

                // Close IOService:
                kernResult = IOServiceClose(connect);
                if (kernResult == KERN_SUCCESS) {
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: IOServiceClose() was successfull.\n");
                }
                else {
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: IOServiceClose returned 0x%08x\n\n", kernResult);
                }

                goto error_abort;
            }

            if (PsychPrefStateGet_Verbosity() > 2) {
                printf("PTB-INFO: Connection to Psychtoolbox kernel support driver instance #%i (Revision %i) established.\n", numKernelDrivers, revision);
            }

            // Store the connect handle for this instance:
            displayConnectHandles[numKernelDrivers] = connect;

            // Query and assign GPU info:
            PsychOSKDGetGPUInfo(connect, numKernelDrivers);

            // Skip Intel gpu's, unless the PSYCH_ALLOW_DANGEROUS env variable is set:
            // Intel IGP's have a design defect which can cause machine hard lockup if multiple
            // regs are accessed simultaneously! As we can't serialize our MMIO reads with the
            // kms-driver, using our MMIO code on Intel is unsafe. Horrible crashes are reported
            // against Haswell on the freedesktop bug tracker for this issue.
            if ((fDeviceType[numKernelDrivers] == kPsychIntelIGP) && !getenv("PSYCH_ALLOW_DANGEROUS")) {
                if (PsychPrefStateGet_Verbosity() > 2) {
                    printf("PTB-INFO: Disconnecting from kernel driver instance #%i for detected Intel GPU for safety reasons. setenv('PSYCH_ALLOW_DANGEROUS', '1') to override.\n", numKernelDrivers);
                }

                displayConnectHandles[numKernelDrivers] = IO_OBJECT_NULL;

                // Call shutdown method:
                kern_return_t kernResult = CallKDSimpleMethod(connect, kMyUserClientClose);
                if (kernResult == KERN_SUCCESS) {
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: CallKDSimpleMethod Closing was successfull.\n");
                }
                else {
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: CallKDSimpleMethod Closing failed with kernel return code 0x%08x.\n\n", kernResult);
                }

                // Close IOService:
                kernResult = IOServiceClose(connect);
                if (kernResult == KERN_SUCCESS) {
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: IOServiceClose() was successfull.\n");
                }
                else {
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: IOServiceClose returned 0x%08x\n\n", kernResult);
                }

                connect = IO_OBJECT_NULL;

                // Closed connection to Intel-instance of the driver. Skip to beginning of
                // enumeration loop to see if we get some alternative, e.g., discrete GPU:
                continue;
            }

            // A word of warning is due for users of outdated Rev. 0 kernel drivers on AMD/ATI GPU's with pre DCE-4 display engine: They shall upgrade or suffer.
            if ((fDeviceType[numKernelDrivers] == kPsychRadeon) && (fCardType[numKernelDrivers] < 40) && (revision < 1) && (PsychPrefStateGet_Verbosity() > 1)) {
                printf("PTB-INFO: You use an outdated Psychtoolbox kernel driver of revision %i. Please upgrade to the latest driver of at least revision 1.\n", revision);
                printf("PTB-INFO: With the outdated driver, robustness of > 8 bits per color displays (10 bit framebuffer, Bits+, Datapixx etc.) will be limited \n");
                printf("PTB-INFO: on your AMD graphics card due to limitations of this driver. You will need to install the latest 64-Bit driver, which will\n");
                printf("PTB-INFO: require you to run a 64-Bit kernel on your machine. The 32-Bit driver is no longer updated, so it will stay outdated.\n");
            }

            // Need to setup crtc offsets separate for use in PsychGraphicsHardwareHALSupport.c on AMD hw:
            if (fDeviceType[numKernelDrivers] == kPsychRadeon) {
                // Setup for DCE-4/5/6/8:
                if ((fCardType[numKernelDrivers] == 40) || (fCardType[numKernelDrivers] == 50) || (fCardType[numKernelDrivers] == 60) || (fCardType[numKernelDrivers] == 80)) {
                    // Offset of crtc blocks of evergreen gpu's for each of the six possible crtc's:
                    crtcoff[0] = EVERGREEN_CRTC0_REGISTER_OFFSET;
                    crtcoff[1] = EVERGREEN_CRTC1_REGISTER_OFFSET;
                    crtcoff[2] = EVERGREEN_CRTC2_REGISTER_OFFSET;
                    crtcoff[3] = EVERGREEN_CRTC3_REGISTER_OFFSET;
                    crtcoff[4] = EVERGREEN_CRTC4_REGISTER_OFFSET;
                    crtcoff[5] = EVERGREEN_CRTC5_REGISTER_OFFSET;
                }

                // Setup for DCE-10/11:
                if ((fCardType[numKernelDrivers] == 100) || (fCardType[numKernelDrivers] == 110)) {
                    // DCE-10/11 of the "Volcanic Islands" gpu family uses (mostly) the same register specs,
                    // but the offsets for the different CRTC blocks are different wrt. to pre DCE-10. Therefore
                    // need to initialize the offsets differently. Also, some of these parts seem to support up
                    // to 7 display engines instead of the old limit of 6 engines:

                    // Offset of crtc blocks of Volcanic Islands DCE-10/11 gpu's for each of the possible crtc's:
                    crtcoff[0] = DCE10_CRTC0_REGISTER_OFFSET;
                    crtcoff[1] = DCE10_CRTC1_REGISTER_OFFSET;
                    crtcoff[2] = DCE10_CRTC2_REGISTER_OFFSET;
                    crtcoff[3] = DCE10_CRTC3_REGISTER_OFFSET;
                    crtcoff[4] = DCE10_CRTC4_REGISTER_OFFSET;
                    crtcoff[5] = DCE10_CRTC5_REGISTER_OFFSET;
                    crtcoff[6] = DCE10_CRTC6_REGISTER_OFFSET;
                }
            }

            // Increment instance count by one:
            numKernelDrivers++;
        }
        // Iterate to next GPU / Driver:
    }

    // Is this a hybrid graphics system with two GPUs, ie. an integrated Intel IGP and
    // a discrete NVidia or AMD GPU? If so, is our currently selected default activeGPU
    // (instance 0) the Intel IGP? If so, then that's probably not what we want and the
    // mapping of kernel driver instances to GPUs got mangled/switched at driver load.
    // OSX default behaviour is to power up and switch to the high-performance discrete
    // non-Intel GPU as soon as PTB creates its first OpenGL rendering context while opening
    // its first onscreen window, then sticking to that GPU for the remainder of the Matlab/
    // Octave session. For this reason we should also default to choosing the discrete GPU
    // for all low-level operations. Check if we do so and make it so:
    if ((numKernelDrivers == 2) && (fDeviceType[activeGPU] == kPsychIntelIGP)) {
        activeGPU = 1 - activeGPU;
        OSMemoryBarrier();
        if (PsychPrefStateGet_Verbosity() > 2) {
            printf("PTB-INFO: Switching to kernel driver instance #%i in hybrid graphics system, assuming i am attached to discrete non-Intel GPU.", activeGPU);
        }

        // PsychOSKDGetBeamposition() has a way of recovering from a wrong choice here.
        // If the Intel IGP should be actually used as the GPU of choice, e.g., because
        // the system is low on battery power, some other condition prevents use of the
        // discrete GPU, or the user has forcefully enabled the IGP via some 3rd party
        // tool, e.g., via "gfxCardStatus", then the discrete GPU will be powered down at
        // time of beamposition query --> query results bogus results on a non-Intel GPU in
        // a hybrid graphics system --> the function will switch to the alternative GPU and
        // try to recover.
        // This should hopefully catch the special case of Intel IGP, and as a side-effect
        // also recover from wrong mapping if we are on a older generation hybrid gfx system
        // with two NVidia or AMD GPUs, one high perf, one low perf.
        //
        // This whole logic does not cover multi-gpu systems without graphics switching in
        // any way, ie., MacPro's with multi-gpus are not automatically treated correctly,
        // but MacBookPro's with hybrid graphics have a decent chance of working.
    }

    // A bit of a hack for now: Allow usercode to select which gpu in a multi-gpu
    // system should be used for low-level mmio based features. If the environment
    // variable PSYCH_USE_GPUIDX is set to a number, it will try to use that GPU:
    // TODO: Replace this by true multi-gpu support and - far in the future? -
    // automatic mapping of screens to gpu's:
    if (getenv("PSYCH_USE_GPUIDX")) {
        activeGPU = atoi(getenv("PSYCH_USE_GPUIDX"));
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Will try to use GPU number %i for low-level access during this session, as requested by usercode.\n", activeGPU);
    }

    // Perform auto-detection of screen to head mappings:
    PsychAutoDetectScreenToHeadMappings(fNumDisplayHeads[activeGPU]);

error_abort:

    // Release the io_iterator_t now that we're done with it.
    IOObjectRelease(iterator);

    // Done.
    return;
}

/*
 * Return identifying information about GPU for a given screen screenNumber:
 *
 * Returns TRUE on success, and the actual info in int variables, FALSE if info
 * not available:
 * Input: screenNumber of the screen for which to query GPU.
 *
 * Output: All optional - NULL == Don't return info.
 *
 * gpuMaintype = Basically what vendor.
 * gpuMinortype = Vendor specific id meaningful to us to define a certain class or generation of hardware.
 * pciDeviceId = The PCI device id.
 * numDisplayHeads = Maximum number of crtc's.
 *
 */
psych_bool PsychGetGPUSpecs(int screenNumber, int* gpuMaintype, int* gpuMinortype, int* pciDeviceId, int* numDisplayHeads)
{
    if (!PsychOSIsKernelDriverAvailable(screenNumber)) return(FALSE);

    if (gpuMaintype) *gpuMaintype = fDeviceType[activeGPU];
    if (gpuMinortype) *gpuMinortype = fCardType[activeGPU];
    if (pciDeviceId) *pciDeviceId = fPCIDeviceId[activeGPU];

    if (numDisplayHeads) *numDisplayHeads = fNumDisplayHeads[activeGPU];

    return(TRUE);
}

// Try to detach to kernel level ptb support driver and tear down everything:
void PsychOSShutdownPsychtoolboxKernelDriverInterface(void)
{
    io_connect_t connect;
    kern_return_t kernResult;
    int i;

    for (i = 0; i < numKernelDrivers; i++) {
        connect = displayConnectHandles[i];

        // Close IOService:
        kernResult = IOServiceClose(connect);
        if (kernResult == KERN_SUCCESS) {
            if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: IOServiceClose() for driver instance %i was successfull.\n", i);
        }
        else {
            if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: IOServiceClose returned 0x%08x for driver instance %i.\n\n", kernResult, i);
        }
    }

    // Ok, whatever happened, we're detached (for good or bad):
    numKernelDrivers = 0;

    return;
}

psych_bool PsychOSIsKernelDriverAvailable(int screenId)
{
    return((numKernelDrivers > 0) ? TRUE : FALSE);
}

io_connect_t PsychOSCheckKDAvailable(int screenId, unsigned int * status)
{
    io_connect_t connect = displayConnectHandles[activeGPU];

    if (numKernelDrivers <= 0) {
        if (status) *status = kIOReturnNotFound;
        return(0);
    }

    if (!connect) {
        if (status) *status = kIOReturnNoDevice;
        if (PsychPrefStateGet_Verbosity() > 6) {
            printf("PTB-DEBUGINFO: Could not access kernel driver connection %i for screenId %i - No such connection.\n", activeGPU, screenId);
        }
        return(0);
    }

    if (status) *status = kIOReturnSuccess;
    return(connect);
}

kern_return_t PsychOSKDDispatchCommand(io_connect_t connect, const PsychKDCommandStruct* inStruct, PsychKDCommandStruct* outStruct, unsigned int* status)
{
    size_t outputStructCnt = sizeof(PsychKDCommandStruct);

    kern_return_t kernResult = IOConnectCallStructMethod((mach_port_t) connect, kPsychKDDispatchCommand, (const void*) inStruct, sizeof(PsychKDCommandStruct), outStruct, &outputStructCnt);

    if (status) *status = kernResult;
    if (kernResult != kIOReturnSuccess) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Kernel driver command dispatch failure for code %lx (Kernel error code: %lx).\n", inStruct->command, kernResult);
    }

    return kernResult;
}

unsigned int PsychOSKDReadRegister(int screenId, unsigned int offset, unsigned int* status)
{
    // Have syncCommand locally defined, ie. on threads local stack: Important for thread-safety, e.g., for async-flip etc.:
    PsychKDCommandStruct syncCommand;

    // Check availability of connection:
    io_connect_t connect;
    if (!(connect = PsychOSCheckKDAvailable(screenId, status))) return(0xffffffff);

    // Set command code for register read:
    syncCommand.command = kPsychKDReadRegister;

    // Register offset is arg 0:
    syncCommand.inOutArgs[0] = offset;

    // Issue request to driver:
    kern_return_t kernResult = PsychOSKDDispatchCommand(connect, &syncCommand, &syncCommand, status);
    if (kernResult != KERN_SUCCESS) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Kernel driver register read failed for register %lx (Kernel error code: %lx).\n", offset, kernResult);
        // A value of 0xffffffff signals failure:
        return(0xffffffff);
    }

    // Return readback register value:
    return((int) syncCommand.inOutArgs[0]);
}

unsigned int PsychOSKDWriteRegister(int screenId, unsigned int offset, unsigned int value, unsigned int* status)
{
    // Have syncCommand locally defined, ie. on threads local stack: Important for thread-safety, e.g., for async-flip etc.:
    PsychKDCommandStruct syncCommand;

    // Check availability of connection:
    io_connect_t connect;
    if (!(connect = PsychOSCheckKDAvailable(screenId, status))) return(0xffffffff);

    // Set command code for display sync:
    syncCommand.command = KPsychKDWriteRegister;
    syncCommand.inOutArgs[0] = offset;
    syncCommand.inOutArgs[1] = value;

    // Issue request to driver:
    kern_return_t kernResult = PsychOSKDDispatchCommand(connect, &syncCommand, &syncCommand, status);
    if (kernResult != KERN_SUCCESS) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Kernel driver register write failed for register %lx, value %lx (Kernel error code: %lx).\n", offset, value, kernResult);
        // A value of 0xffffffff signals failure:
        return(0xffffffff);
    }

    // Return success:
    return(0);
}

// Synchronize display screens video refresh cycle. See PsychSynchronizeDisplayScreens() for help and details...
PsychError PsychOSSynchronizeDisplayScreens(int *numScreens, int* screenIds, int* residuals, unsigned int syncMethod, double syncTimeOut, int allowedResidual)
{
    // Have syncCommand locally defined, ie. on threads local stack: Important for thread-safety, e.g., for async-flip etc.:
    PsychKDCommandStruct syncCommand;

    int screenId = 0;
    double abortTimeOut, now;
    int residual;

    // Check availability of connection:
    io_connect_t connect;
    unsigned int status;
    kern_return_t kernResult;

    // No support for other methods than fast hard sync:
    if (syncMethod > 1) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not execute display resync operation with requested non hard sync method. Not supported for this setup and settings.\n");
        return(PsychError_unimplemented);
    }

    // The current implementation only supports syncing all heads of a single card
    if (*numScreens <= 0) {
        // Resync all displays requested: Choose screenID zero for connect handle:
        screenId = 0;
    }
    else {
        // Resync of specific display requested: We only support resync of all heads of a single multi-head card,
        // therefore just choose the screenId of the passed master-screen for resync handle:
        screenId = screenIds[0];
    }

    if (!(connect = PsychOSCheckKDAvailable(screenId, &status))) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not execute display resync operation for master screenId %i. Not supported for this setup and settings.\n", screenId);
        return(PsychError_unimplemented);
    }

    // Setup deadline for abortion or repeated retries:
    PsychGetAdjustedPrecisionTimerSeconds(&abortTimeOut);
    abortTimeOut+=syncTimeOut;
    residual = INT_MAX;

    // Repeat until timeout or good enough result:
    do {
        // If this isn't the first try, wait 0.5 secs before retry:
        if (residual != INT_MAX) PsychWaitIntervalSeconds(0.5);

        residual = INT_MAX;

        // Set command code for display sync:
        syncCommand.command = kPsychKDFastSyncAllHeads;

        // Issue request to driver:
        kernResult = PsychOSKDDispatchCommand(connect, &syncCommand, &syncCommand, &status);
        if (kernResult == KERN_SUCCESS) {
            residual = (int) syncCommand.inOutArgs[0];
            if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Graphics display heads resynchronized. Residual vertical beamposition error is %ld scanlines.\n", residual);
        }
        else {
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Graphics display head synchronization failed (Kernel error code: %lx).\n", kernResult);
            break;
        }

        // Timestamp:
        PsychGetAdjustedPrecisionTimerSeconds(&now);
    } while ((now < abortTimeOut) && (abs(residual) > allowedResidual) && (kernResult == KERN_SUCCESS));

    // Return residual value if wanted:
    if (residuals) {
        residuals[0] = residual;
    }

    if (abs(residual) > allowedResidual) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Failed to synchronize heads down to the allowable residual of +/- %i scanlines. Final residual %i lines.\n", allowedResidual, residual);
    }

    // TODO: Error handling not really worked out...
    if (residual == INT_MAX || kernResult != KERN_SUCCESS) return(PsychError_system);

    // Done.
    return(PsychError_none);
}

int PsychOSKDGetBeamposition(int screenId)
{
    int beampos, vblbias, vbltotal;

    // Have syncCommand locally defined, ie. on threads local stack: Important for thread-safety, e.g., for async-flip etc.:
    PsychKDCommandStruct syncCommand;

    // Check availability of connection:
    io_connect_t connect;
    if (!(connect = PsychOSCheckKDAvailable(screenId, NULL))) {
        // Beampos queries unavailable:
        if (PsychPrefStateGet_Verbosity() > 21) printf("PTB-DEBUG: Kernel driver based beamposition queries unavailable for screenId %i.\n", screenId);
        return(-1);
    }

    // Set command code for beamposition query:
    syncCommand.command = kPsychKDGetBeamposition;

    // Assign headid for this screen:
    syncCommand.inOutArgs[0] = PsychScreenToCrtcId(screenId, 0);

    // Issue request:
    kern_return_t kernResult = PsychOSKDDispatchCommand(connect,  &syncCommand, &syncCommand, NULL);
    if (kernResult != KERN_SUCCESS) {
        if (PsychPrefStateGet_Verbosity() > 6) printf("PTB-ERROR: Kernel driver beamposition query failed (Kernel error code: %lx).\n", kernResult);
        // A value of -1 signals beamposition queries unsupported or invalid:
        return(-1);
    }

    beampos = (int) syncCommand.inOutArgs[0];

    // Reasonable result? If beampos is very large (ie., much larger than the display)
    // and this is a hybrid-graphics machine with more than one GPU (numKernelDrivers > 1)
    // then the most likely reason for the result is that we are querying the wrong kernel
    // driver instance, and therefore the wrong GPU. E.g., NVidia discrete GPUs report a value
    // of 0xffff aka 65535 here if they are powered down and disconnected from the bus due to
    // the integrated Intel card being active. We try if this is the case, if the activeGPU is
    // not the Intel integrated IGP, ie., it is a discrete GPU from NVidia or AMD.
    if ((beampos > 16384) && (numKernelDrivers == 2) && (fDeviceType[activeGPU] != kPsychIntelIGP)) {
        // Probably wrong, powered-down, discrete GPU. Switch the GPU/Driver instance to other GPU
        // for all future operations in the hope that the other GPU is actually operational:
        activeGPU = 1 - activeGPU;

        // Make sure our SMP sibling cores get notified:
        OSMemoryBarrier();

        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Switching kernel driver interface to alternative GPU %i in system until Screen() reset.\n", activeGPU);

        // Call ourselves recursively in the hope we get a better result this time.
        // If this also fails, it is game over:
        return(PsychOSKDGetBeamposition(screenId));
    }

    // Catch still invalid values and map to "unsupported":
    if (beampos > 16384) return(-1);

    // Apply corrective offsets if any (i.e., if non-zero):
    PsychGetBeamposCorrection(screenId, &vblbias, &vbltotal);
    beampos = beampos - vblbias;
    if (beampos < 0) beampos = vbltotal + beampos;

    // Return queried position:
    return(beampos);
}

// Try to change hardware dither mode on GPU:
void PsychOSKDSetDitherMode(int screenId, unsigned int ditherOn)
{
    // Have syncCommand locally defined, ie. on threads local stack: Important for thread-safety, e.g., for async-flip etc.:
    PsychKDCommandStruct syncCommand;

    // Check availability of connection:
    io_connect_t connect;
    if (!(connect = PsychOSCheckKDAvailable(((screenId >= 0) ? screenId : 0), NULL))) {
        // Dither control unavailable:
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Kernel driver based dither control unavailable for screenId %i.\n", screenId);
        return;
    }

    // Set command code for dither control:
    syncCommand.command = kPsychKDSetDitherMode;

    // Assign headid for this screen:
    syncCommand.inOutArgs[0] = (unsigned int) ((screenId >= 0) ? PsychScreenToCrtcId(screenId, 0) : -1 * screenId);

    // Assign dither setting:
    syncCommand.inOutArgs[1] = ditherOn;

    // Issue request:
    kern_return_t kernResult = PsychOSKDDispatchCommand(connect,  &syncCommand, &syncCommand, NULL);
    if (kernResult != KERN_SUCCESS) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Kernel driver dither control call failed (Kernel error code: %lx).\n", kernResult);
        return;
    }

    return;
}

unsigned int PsychOSKDGetRevision(io_connect_t connect)
{
    PsychKDCommandStruct syncCommand;
    IOByteCount structSize1 = sizeof(PsychKDCommandStruct);

    // Set command code for revision query:
    syncCommand.command = kPsychKDGetRevision;

    // Issue request:
    kern_return_t kernResult = PsychOSKDDispatchCommand(connect, &syncCommand, &syncCommand, NULL);
    if (kernResult != KERN_SUCCESS) {
        printf("PTB-ERROR: Kernel driver revision read failed (Kernel error code: %lx).\n", kernResult);
        // A value of 0xffffffff signals failure:
        return(0xffffffff);
    }

    // Return revision:
    return((unsigned int) syncCommand.inOutArgs[0]);
}

static void PsychOSKDGetGPUInfo(io_connect_t connect, int slot)
{
    PsychKDCommandStruct syncCommand;
    IOByteCount structSize1 = sizeof(PsychKDCommandStruct);

    // Set command code for gpu info query:
    syncCommand.command = kPsychKDGetGPUInfo;

    // Issue request:
    kern_return_t kernResult = PsychOSKDDispatchCommand(connect, &syncCommand, &syncCommand, NULL);
    if (kernResult != KERN_SUCCESS) {
        printf("PTB-ERROR: Kernel driver gpu info read failed (Kernel error code: %lx).\n", kernResult);
        // A value of 0xffffffff signals failure:
        return;
    }

    // Assign:
    fDeviceType[slot] = syncCommand.inOutArgs[0];
    fPCIDeviceId[slot] = syncCommand.inOutArgs[1];
    fCardType[slot] = syncCommand.inOutArgs[2];
    fNumDisplayHeads[slot] = syncCommand.inOutArgs[3];

    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: GPU %i - GPU-Vendor %i, PCIId %x, GPU-Type %i [x%x], numDisplayHeads %i.\n", slot, fDeviceType[slot], fPCIDeviceId[slot], fCardType[slot], fCardType[slot], fNumDisplayHeads[slot]);
    return;
}

unsigned int PsychOSKDGetLUTState(int screenId, unsigned int head, unsigned int debug)
{
    PsychKDCommandStruct syncCommand;
    IOByteCount structSize1 = sizeof(PsychKDCommandStruct);

    // Check availability of connection:
    io_connect_t connect;
    if (!(connect = PsychOSCheckKDAvailable(screenId, NULL))) return(0xffffffff);

    // Set command code for LUT check:
    syncCommand.command = kPsychKDGetLUTState;
    syncCommand.inOutArgs[0] = head;
    syncCommand.inOutArgs[1] = debug;

    // Issue request:
    kern_return_t kernResult = PsychOSKDDispatchCommand(connect, &syncCommand, &syncCommand, NULL);
    if (kernResult != KERN_SUCCESS) {
        printf("PTB-ERROR: Kernel driver lut state read failed (Kernel error code: %lx).\n", kernResult);
        // A value of 0xffffffff signals failure:
        return(0xffffffff);
    }

    // Return lut state:
    return((unsigned int) syncCommand.inOutArgs[0]);
}

unsigned int PsychOSKDLoadIdentityLUT(int screenId, unsigned int head)
{
    PsychKDCommandStruct syncCommand;
    IOByteCount structSize1 = sizeof(PsychKDCommandStruct);

    // Check availability of connection:
    io_connect_t connect;
    if (!(connect = PsychOSCheckKDAvailable(screenId, NULL))) return(0);

    // Set command code for identity LUT load:
    syncCommand.command = kPsychKDSetIdentityLUT;
    syncCommand.inOutArgs[0] = head;

    // Issue request:
    kern_return_t kernResult = PsychOSKDDispatchCommand(connect, &syncCommand, &syncCommand, NULL);
    if (kernResult != KERN_SUCCESS) {
        printf("PTB-ERROR: Kernel driver identity lut setup failed (Kernel error code: %lx).\n", kernResult);
        // A value of 0 signals failure:
        return(0);
    }

    // Return lut setup rc:
    return((unsigned int) syncCommand.inOutArgs[0]);
}
