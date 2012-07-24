/*
	PsychToolbox3/Source/windows/Screen/PsychScreenGlue.c
	
	PLATFORMS:	
	
		This is the Windows version only.  
				
	AUTHORS:
	
	Mario Kleiner		mk		mario.kleiner at tuebingen.mpg.de

	HISTORY:
	
	12/27/05                mk              Wrote it. Derived from Allens OSX version.
        							
	DESCRIPTION:
	
		Functions in this file comprise an abstraction layer for probing and controlling screen state.  
		
		Each C function which implements a particular Screen subcommand should be platform neutral.  For example, the source to SCREENPixelSizes() 
		should be platform-neutral, despite that the calls in OS X and Windows to detect available pixel sizes are
		different.  The platform specificity is abstracted out in C files which end it "Glue", for example PsychScreenGlue, PsychWindowGlue, 
		PsychWindowTextClue.
	
		In addition to glue functions for windows and screen there are functions which implement shared functionality between between Screen commands,
		such as ScreenTypes.c and WindowBank.c. 
			
	NOTES:
	
	TO DO: 
	
	* DirectDraw objects are not yet destroyed at Screen flush time. This will leak a bit of memory...
*/


#include "Screen.h"
#include <Windows.h>

#ifdef PTBOCTAVE3MEX
// For some weird reason, this isn't defined on Octave-3 builds, although
// all other relevant information that would usually go along with this in
// MultiMon.h (MSVC) or multimon.h (DirectX-SDK) is present.
#define MONITOR_DEFAULTTONULL       0x00000000
#define MONITOR_DEFAULTTOPRIMARY    0x00000001
#define MONITOR_DEFAULTTONEAREST    0x00000002

#define MONITORINFOF_PRIMARY        0x00000001

#endif

// We need to define this and include Multimon.h to allow for enumeration of multiple display screens:
#define COMPILE_MULTIMON_STUBS
#include <Multimon.h>

// Include DirectDraw header for access to the GetScanLine() function:
#include <ddraw.h>

// Module handle for the DirectDraw library 'ddraw.dll': Or 0 if unsupported.
HMODULE ddrawlibrary = 0;

typedef HRESULT (WINAPI * LPDIRECTDRAWCREATE)(GUID FAR *lpGUID, LPDIRECTDRAW FAR *lplpDD, IUnknown FAR *pUnkOuter);
LPDIRECTDRAWCREATE      PsychDirectDrawCreate = NULL;
LPDIRECTDRAWENUMERATEEX PsychDirectDrawEnumerateEx = NULL;

// file local variables

// Maybe use NULLs in the settings arrays to mark entries invalid instead of using psych_bool flags in a different array.   
static psych_bool			displayLockSettingsFlags[kPsychMaxPossibleDisplays];
static CFDictionaryRef		displayOriginalCGSettings[kPsychMaxPossibleDisplays];        	//these track the original video state before the Psychtoolbox changed it.  
static psych_bool			displayOriginalCGSettingsValid[kPsychMaxPossibleDisplays];
static CFDictionaryRef		displayOverlayedCGSettings[kPsychMaxPossibleDisplays];        	//these track settings overlayed with 'Resolutions'.  
static psych_bool			displayOverlayedCGSettingsValid[kPsychMaxPossibleDisplays];
static CGDisplayCount 		numDisplays;
static CGDirectDisplayID 	displayCGIDs[kPsychMaxPossibleDisplays];				// DC device contexts for displays.
static char*                displayDeviceName[kPsychMaxPossibleDisplays];			// Windows internal monitor device name. Default display has NULL
static int 	                displayDeviceStartX[kPsychMaxPossibleDisplays];			// Top-Left corner of display on virtual screen. Default display has (0,0).
static int 	                displayDeviceStartY[kPsychMaxPossibleDisplays];
static LPDIRECTDRAW         displayDeviceDDrawObject[kPsychMaxPossibleDisplays];	// Pointer to associated DirectDraw object, if any. NULL otherwise.
static HMONITOR				displayDevicehMonitor[kPsychMaxPossibleDisplays];		// HMonitor handle of associated monitor, NULL otherwise.
static GUID					displayDeviceGUID[kPsychMaxPossibleDisplays];			// GUID for DDRAW display device, unless displayDevicelpGUIDValid < 2.
static int					displayDeviceGUIDValid[kPsychMaxPossibleDisplays];		// GUID for DDRAW display device valid? 0 = No, 1 = Must be NULL, 2 = Valid.
static int					ddrawnumDisplays;

static psych_bool enableVBLBeamposWorkaround = FALSE;	// Is the special workaround for beamposition queries needed?

//file local functions
void InitCGDisplayIDList(void);
void PsychLockScreenSettings(int screenNumber);
void PsychUnlockScreenSettings(int screenNumber);
psych_bool PsychCheckScreenSettingsLock(int screenNumber);
psych_bool PsychGetCGModeFromVideoSetting(CFDictionaryRef *cgMode, PsychScreenSettingsType *setting);
void PsychTestDDrawBeampositionQueries(int screenNumber);

// This is actually a function in PsychWindowGlue.c, we redefine the prototype here to make compiler happy:
extern psych_bool ChangeScreenResolution (int screenNumber, int width, int height, int bitsPerPixel, int fps);

//Initialization functions
void InitializePsychDisplayGlue(void)
{
    int i;
    
    //init the display settings flags.
    for(i=0;i<kPsychMaxPossibleDisplays;i++){
        displayLockSettingsFlags[i] = FALSE;
        displayOriginalCGSettingsValid[i] = FALSE;
        displayOverlayedCGSettingsValid[i] = FALSE;
		displayDeviceDDrawObject[i] = NULL;
		displayDevicehMonitor[i] = NULL;
		displayDeviceGUIDValid[i] = 0;
    }
 
	// Disable beampos workaround by default:
	enableVBLBeamposWorkaround = FALSE;
	
    //init the list of Core Graphics display IDs.
    InitCGDisplayIDList();
}

// This callback function is called by Windows EnumDisplayMonitors() function for each
// detected display device: We can happily ignore all provided parameters, except for the
// hMonitor struct which contains the Windows internal name for the detected display. We
// need to pass this name string to a variety of Windows-Functions to refer to the monitor
// of interest.
BOOL CALLBACK MonitorEnumProc(HMONITOR hMonitor, HDC hdcMonitor, LPRECT lprcMonitor, LPARAM dwData);
BOOL WINAPI PsychDirectDrawEnumProc(GUID FAR* lpGUID, LPSTR lpDriverDescription, LPSTR lpDriverName, LPVOID displayIdx, HMONITOR hMonitor);

BOOL CALLBACK MonitorEnumProc(HMONITOR hMonitor, HDC hdcMonitor, LPRECT lprcMonitor, LPARAM dwData)
{
	MONITORINFOEX moninfo;

	// hMonitor is the handle to the monitor info. Resolve it to a moninfo information struct:
	moninfo.cbSize = sizeof(MONITORINFOEX);
	GetMonitorInfo(hMonitor, (LPMONITORINFO) &moninfo);

	// Query and copy the display device name into our own screenNumber->Name mapping array:
	displayDeviceName[numDisplays] = (char*) malloc(256);
	strncpy(displayDeviceName[numDisplays], moninfo.szDevice, 256);

	// Query and copy the top-left corner of the monitor:
	displayDeviceStartX[numDisplays] = moninfo.rcMonitor.left;
	displayDeviceStartY[numDisplays] = moninfo.rcMonitor.top;

	// Create a device context for this display and store it in our displayCGIDs array:
	displayCGIDs[numDisplays] = CreateDC(displayDeviceName[numDisplays], displayDeviceName[numDisplays], NULL, NULL);

	// Store HMONITOR handle for DirectDraw device enumeration:
	displayDevicehMonitor[numDisplays] = hMonitor;

	// Increase global counter of available separate displays:
	numDisplays++;

	// Return and ask system to continue display enumeration:
	return(TRUE);
}

// Callback for DirectDraw display enumeration which complements monitor enumeration above:
// Called by PsychDirectDrawEnumerateEx() below...
BOOL WINAPI PsychDirectDrawEnumProc(GUID FAR* lpGUID, LPSTR lpDriverDescription, LPSTR lpDriverName, LPVOID displayIdx, HMONITOR hMonitor)
{
	// Retrieve index into our display array for the display to enumerate / map in this callback:
	int idx = (int) displayIdx;
	
	// Is the passed in hMonitor identical to the hMonitor of target display "idx"?
	if (hMonitor == displayDevicehMonitor[idx]) {
		// Bingo! Store the GUID globally unique display device handle for this display. We'll need this
		// to create the matching DDRAW device for that display:
		if (lpGUID) {
			// Matched, with unique GUID: Assign it.
			displayDeviceGUID[idx] = *lpGUID;
			displayDeviceGUIDValid[idx] = 2;
		}
		else {
			// Matched, but lpGUID is NULL:
			displayDeviceGUIDValid[idx] = 1;
		}
		// printf("PTB-DEBUG: In DDRAW enumeration: Screenid %i: '%s' and '%s' -> lpGUID %p.\n", idx, (char*) lpDriverDescription, (char*) lpDriverName, lpGUID);
		
		// Finalize enumeration pass:
		return(FALSE);
	}
	else {
		// printf("PTB-DEBUG: In enum call for screenid %i: FAIL for '%s' , '%s'\n", idx, (char*) lpDriverDescription, (char*) lpDriverName);
	}

	// Nope. Continue enumeration:
	return(TRUE);
}

void InitCGDisplayIDList(void)
{
  int i, w1, w2, h1, h2;
  psych_uint32 beampos = 100000;
  HRESULT rc;

  // We always provide the full (virtual) desktop as screen number zero. This way,
  // queries to screen 0 will always provide the global settings and dimensions of
  // the full desktop (either single display, or extended desktop on multi-display systems).
  // Opening an onscreen window on screen 0 will always yield a window covering the full
  // desktop, possibly spanning multiple physical display devices. Very useful for standard
  // dual display stereo applications in stereoMode 4.
  ddrawnumDisplays = 0;
  numDisplays=1;
  displayCGIDs[0]=GetDC(GetDesktopWindow());
  displayDeviceName[0] = NULL;
  displayDeviceStartX[0] = 0;
  displayDeviceStartY[0] = 0;

  // Since Windows Vista/7 we need to manually load the DirectDraw runtime DLL and bind
  // the functions we need manually. Thanks Microsoft for creating pointless extra work!
  ddrawlibrary = LoadLibrary("ddraw.dll");
  if (ddrawlibrary) {
    // Load success. Dynamically bind the relevant functions:
    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: DirectDraw DLL available. Binding controls ...");

    PsychDirectDrawCreate      = (LPDIRECTDRAWCREATE) GetProcAddress(ddrawlibrary, "DirectDrawCreate");
    PsychDirectDrawEnumerateEx = (LPDIRECTDRAWENUMERATEEX) GetProcAddress(ddrawlibrary, "DirectDrawEnumerateExA");

    if (PsychDirectDrawCreate && PsychDirectDrawEnumerateEx) {
      // Mark DirectDraw as supported:
      if (PsychPrefStateGet_Verbosity() > 5) printf(" ...done\n");
    }
    else {
      FreeLibrary(ddrawlibrary);
      ddrawlibrary = 0;
        if (PsychPrefStateGet_Verbosity() > 0) {
            printf("PTB-WARNING: Could not attach to DirectDraw library ddraw.dll - Trouble ahead!\n");
            printf("PTB-WARNING: DirectDrawCreate() = %p : DirectDrawEnumerateEx() = %p\n", PsychDirectDrawCreate, PsychDirectDrawEnumerateEx);
        }
    }
  }
  else {
    if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: DirectDraw runtime DLL unsupported! Trouble ahead!\n"); 
  }

  // Create a DirectDraw object for the primary display device, i.e. the single display on a
  // single display setup or the display device corresponding to the desktop on a multi-display setup:
  if ((ddrawlibrary == 0) || (PsychDirectDrawCreate((GUID FAR *) DDCREATE_HARDWAREONLY, &(displayDeviceDDrawObject[0]), NULL) != DD_OK)) {
    // Failed to create Direct Draw object:
    displayDeviceDDrawObject[0]=NULL;
    printf("PTB-WARNING: Failed to create DirectDraw interface for primary display. Won't be able to generate high-precision 'Flip' timestamps.\n");
  }
  else {
    rc=IDirectDraw_GetScanLine(displayDeviceDDrawObject[0], (LPDWORD) &beampos);
    if (rc!=DD_OK && rc!=DDERR_VERTICALBLANKINPROGRESS) {
		// Beamposition query failed :(
		switch(rc) {
			case DDERR_UNSUPPORTED:
				printf("PTB-INFO: Beamposition query unsupported on this system configuration.\n");
			break;

			default:
				printf("PTB-INFO: Beamposition query test failed: Direct Draw Error.\n");
		}
    }
	else {
		// Mark at least basic beamposition queries for primary display as ok:
		ddrawnumDisplays = 1;
	}
  }
  // Now call M$-Windows monitor enumeration routine. It will call our callback-function
  // MonitorEnumProc() for each detected display device...
  EnumDisplayMonitors(NULL, NULL, MonitorEnumProc, 0);

  // Only one additional display found?
  if (numDisplays <=2) {
    // Single display mode: Enumeration only found the single monitor that constitutes
    // the desktop on a single display setup, so screen 1 is identical to screen 0.
    // In that case, we release screen 1 again as it would be redundant...
    if (numDisplays == 2) DeleteDC(displayCGIDs[1]);
    numDisplays = 1;
  }
  else {
    // At least two different displays enumerated: This is a multi-display setup.
    // Screen 0 is the full desktop. Screens i=1,2,...,n are display monitors 1 to n.
    // Output some info to the command window:
    printf("PTB-INFO: Multi-display setup in explicit multi-display mode detected. Using the following mapping:\n");
    printf("PTB-INFO: Screen 0 corresponds to the full Windows desktop area. Useful for stereo presentations in stereomode=4 ...\n");
    for (i=1; i<numDisplays; i++) {
      printf("PTB-INFO: Screen %i corresponds to the display area of the monitor with the Windows-internal name %s ...\n", i, displayDeviceName[i]);
    }
    
    // Check for sane display dimensions: Our emulation trick for creating a display spanning screen 0 will only work,
    // if the first two physical displays 1 and 2 are of the same resolution/size and if they are arranged so that both
    // together are suitable for a horizontal desktop spanning window, aka they touch each other at their horizontal borders.
    // If the user has a more exotic display setup, e.g., triple-display or monitors with different resolution, (s)he can still
    // use the [rect] parameter when opening a window on screen 0 to enforce arbitrary onscreen window position and size. We'll
    // always create a borderless window on screen 0 when in multidisplay mode...
    PsychGetScreenSize(1, &w1, &h1);
    PsychGetScreenSize(2, &w2, &h2);
    if (w1!=w2 || h1!=h2) {
		printf("PTB-WARNING: Screens 1 and 2 do have non-equal width and height. This will probably cause wrong behaviour\n");
		printf("PTB-WARNING: when trying to open a fullscreen window on Screen 0 that is supposed to fully cover displays 1 and 2.\n");
		printf("PTB-WARNING: In that case, either change your screen settings to matching size and refresh rate in the display control\n");
		printf("PTB-WARNING: panel and retry after a 'clear all', or manually specify a [rect] parameter for the window in the\n");
		printf("PTB-WARNING: Screen('OpenWindow', 0, [color], [rect]); command to manually enforce creation of a proper onscreen window.\n");
	 }
    
    printf("\n"); fflush(NULL);

	// Try to enumerate all DirectDraw devices for all our enumerated display monitors:
	// We'd like to create one DDRAW object per display, so we have per-display beamposition queries available.
	for (i = 1; i < numDisplays; i++) {
		// Enumerate for display screen 'i':
		rc = 0;

		if (PsychDirectDrawEnumerateEx && ((rc = PsychDirectDrawEnumerateEx(PsychDirectDrawEnumProc, (LPVOID) i, DDENUM_ATTACHEDSECONDARYDEVICES)) == DD_OK) && (displayDeviceGUIDValid[i] > 1)) {
			// Success: Create corresponding DDRAW device interface for this screen:
			if ((rc = PsychDirectDrawCreate(&displayDeviceGUID[i], &(displayDeviceDDrawObject[i]), NULL)) != DD_OK) {
				// Failed to create Direct Draw object for this screen:
				// We use the primary display DDRAW object [0] as a fallback:
				displayDeviceDDrawObject[i] = displayDeviceDDrawObject[0];
				printf("PTB-WARNING: Failed to create DirectDraw interface specific to display screen %i [rc = %x]. Will use primary display DirectDraw object as backup and hope for the best.\n", i, rc);
				printf("PTB-WARNING: This means that beamposition queries for high-precision timestamping may not work correctly on your multi-display setup. We'll see...\n");
			}
			else {
				// Mark this displays beamposition query mechanism as ready:
				ddrawnumDisplays++;
			}
		}
		else {
			// Failed to enumerate/detect Direct Draw object for this screen:
			// We use the primary display DDRAW object [0] as a fallback:
			displayDeviceDDrawObject[i] = displayDeviceDDrawObject[0];
			printf("PTB-WARNING: Failed to detect specific DirectDraw interface GUID for display screen %i [rc = %x]. Will use primary display DirectDraw object as backup and hope for the best.\n", i, rc);
			printf("PTB-WARNING: This means that beamposition queries for high-precision timestamping may not work correctly on your multi-display setup. We'll see...\n");
		}
		// Done with this one. Init next display...
	}

	// Do we have per-display beamposition queries available via corresponding DDRAW objects?
	// If ddrawnumDisplays < numDisplays, then we missed some out. In that case use of the beamposition
	// mechanism may be unsafe and we have to use our heuristic to decide if we cope or fail...
	if (ddrawnumDisplays < numDisplays) {
		// Individual per-display queries unavailable! We basically have one query for all displays, which is
		// only safe to use if all displays are configured in exactly the same way on Windows Vista and later:
		if ((numDisplays == 3) && (w1 == w2) && (h1 == h2) && (PsychGetNominalFramerate(1) == PsychGetNominalFramerate(2)) && PsychIsMSVista()) {
			// Vista baby! Keep common beampos mechanism enabled:
			printf("PTB-INFO: This setup is running on a MS-Vista system or later with two displays of same resolution and refresh rate.\n");
			printf("PTB-INFO: I believe it is safe under this conditions to keep high-precision timestamping enabled, but you should run the\n");
			printf("PTB-INFO: 'PerceptualVBLSyncTest' script at least once after any big change in display configuration to verify this.\n\n");
		}
		else {
			// On a multi-display setup in explicit multi-display mode, we disable beamposition queries by default for now.
			// Users are free to override the default though...
			PsychPrefStateSet_VBLTimestampingMode(-1);
			printf("PTB-INFO: This setup is either a pre-Vista (Win2k or WinXP) setup, or it has more than two displays connected, or\n");
			printf("PTB-INFO: it is a dual-display Vista/Windows-7 setup, but display settings for the displays are not identical.\n\n");
			printf("PTB-INFO: Please note that beamposition queries for accurate Screen('Flip') timestamping are disabled by\n");
			printf("PTB-INFO: default on such MS-Windows multi-display setups. If you want to use them, first run the 'PerceptualVBLSyncTest'\n");
			printf("PTB-INFO: script to verify they're working correctly on your setup. Then you can add the command:\n");
			printf("PTB-INFO: Screen('Preference', 'VBLTimestampingMode', 1); at the top of your script to manually enable them.\n");
			printf("PTB-INFO: Usually beamposition queries work correctly if both of your displays are set to the same resolution,\n");
			printf("PTB-INFO: color depths and video refresh rate, but you *must verify this*.\n");
			printf("PTB-INFO: Make also sure that the 'primary monitor' in the display settings is set to your main stimulus \n");
			printf("PTB-INFO: presentation display for more reliable timing and generally less problems.\n");
		}
	}
  }

  // Setup screenId -> display head mappings:
  PsychInitScreenToHeadMappings(numDisplays);

  // Ready.
  return;
}

void PsychGetCGDisplayIDFromScreenNumber(CGDirectDisplayID *displayID, int screenNumber)
{
    if(screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);
    *displayID=displayCGIDs[screenNumber];
}

char* PsychGetDisplayDeviceName(int screenNumber)
{
    if(screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);
	 return(displayDeviceName[screenNumber]);
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
	// Hack hack hack! This special screenNumber == -1 setting with kPsychUseBeampositionQueryWorkaround
	// enables the workaround as a shortcut, bypassing the regular route via PsychTestDDrawBeampositionQueries().
	// This is called solely from within PsychFlipWindowBuffers() as a last resort auto-workaround:
	if ((screenNumber == -1) && (PsychPrefStateGet_ConserveVRAM() & kPsychUseBeampositionQueryWorkaround)) {
		enableVBLBeamposWorkaround = TRUE;
		return;
	}

    if(screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);
    // For now, this is a "No operation" on Windows. Don't know if facilities for
    // capturing the screen exist or are necessary at all. There's a SetCapture() call
    // for HWND window handles, but that applies to windows, not screens and i don't
    // know, what it actually does...
    PsychLockScreenSettings(screenNumber);

	// However we do call the testing/workaround routine for working beamposition queries here:
	PsychTestDDrawBeampositionQueries(screenNumber);
}

/*
    PsychReleaseScreen()    
*/
void PsychReleaseScreen(int screenNumber)
{	
    CGDisplayErr  error=0;
    
    if(screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);
    PsychUnlockScreenSettings(screenNumber);
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

void PsychGetScreenDepths(int screenNumber, PsychDepthType *depths)
{
    int i, rc;
    DEVMODE result;
    
    if(screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);

    // Query all available modes for this display device and retrieve their
    // depth values:
    i=0;
    do {
        // Query next setting (i) for current (NULL) display:
        result.dmSize = sizeof(DEVMODE);
        result.dmDriverExtra = 0;
        rc = EnumDisplaySettings(PsychGetDisplayDeviceName(screenNumber), i, &result);
        i++;

        // Valid setting returned?
        if (rc) {
            // Yes. Add its depth-value to our depth struct:
            PsychAddValueToDepthStruct((int) result.dmBitsPerPel, depths);
        }
    } while (rc!=0);
    return;    
}

/*   PsychGetAllSupportedScreenSettings()
 *
 *	 Queries the display system for a list of all supported display modes, ie. all valid combinations
 *	 of resolution, pixeldepth and refresh rate. Allocates temporary arrays for storage of this list
 *	 and returns it to the calling routine. This function is basically only used by Screen('Resolutions').
 */
int PsychGetAllSupportedScreenSettings(int screenNumber, int outputId, long** widths, long** heights, long** hz, long** bpp)
{
    int i, rc, numPossibleModes;
    DEVMODE result;

    if(screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);

	// First pass: How many modes are supported?
	i=-1;
	do {
        result.dmSize = sizeof(DEVMODE);
        result.dmDriverExtra = 0;
        rc = EnumDisplaySettings(PsychGetDisplayDeviceName(screenNumber), i, &result);
        i++;
	} while (rc!=0);

    // Get a list of avialable modes for the specified display:
    numPossibleModes= i;
	
	// Allocate output arrays: These will get auto-released at exit
	// from Screen():
	*widths = (long*) PsychMallocTemp(numPossibleModes * sizeof(int));
	*heights = (long*) PsychMallocTemp(numPossibleModes * sizeof(int));
	*hz = (long*) PsychMallocTemp(numPossibleModes * sizeof(int));
	*bpp = (long*) PsychMallocTemp(numPossibleModes * sizeof(int));
	
	// Fetch modes and store into arrays:
    for(i=0; i < numPossibleModes; i++) {
        result.dmSize = sizeof(DEVMODE);
        result.dmDriverExtra = 0;
        rc = EnumDisplaySettings(PsychGetDisplayDeviceName(screenNumber), i, &result);
		(*widths)[i] = (long)  result.dmPelsWidth;
		(*heights)[i] = (long)  result.dmPelsHeight;
		(*hz)[i] = (long) result.dmDisplayFrequency;
		(*bpp)[i] = (long) result.dmBitsPerPel;
    }

	return(numPossibleModes);
}

/*
    static PsychGetCGModeFromVideoSettings()
   
*/
psych_bool  PsychGetCGModeFromVideoSetting(CFDictionaryRef *cgMode, PsychScreenSettingsType *setting)
{
  /*
    FIXME - We just return a 1.

    CFArrayRef modeList;
    CFNumberRef n;
    int i, numPossibleModes;
    long width, height, depth, frameRate, tempWidth, tempHeight, tempDepth,  tempFrameRate;
    
    if(setting->screenNumber>=numDisplays)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetScreenDepths() is out of range"); //also checked within SCREENPixelSizes
        
    //adjust parameter formats
    width=(long)PsychGetWidthFromRect(setting->rect);
    height=(long)PsychGetHeightFromRect(setting->rect);
    depth=(long)PsychGetValueFromDepthStruct(0,&(setting->depth));
    frameRate=(long)setting->nominalFrameRate;

    //get a list of avialable modes for the specified display and iterate over the list looking for our mode.
    modeList = CGDisplayAvailableModes(displayCGIDs[setting->screenNumber]);
    numPossibleModes= CFArrayGetCount(modeList);
    for(i=0;i<numPossibleModes;i++){
        *cgMode = CFArrayGetValueAtIndex(modeList,i);			
        n=CFDictionaryGetValue(*cgMode, kCGDisplayWidth);		//width
        CFNumberGetValue(n,kCFNumberLongType, &tempWidth);
        n=CFDictionaryGetValue(*cgMode, kCGDisplayHeight);		//height
        CFNumberGetValue(n,kCFNumberLongType, &tempHeight);
        n=CFDictionaryGetValue(*cgMode, kCGDisplayRefreshRate);	//frequency
        CFNumberGetValue(n, kCFNumberLongType, &tempFrameRate) ;
        n=CFDictionaryGetValue(*cgMode, kCGDisplayBitsPerPixel);	//depth
        CFNumberGetValue(n, kCFNumberLongType, &tempDepth) ;
        if(width==tempWidth && height==tempHeight && frameRate==tempFrameRate && depth==tempDepth)
            return(TRUE);
    }
    return(FALSE);    

  */

  // Dummy assignment:
  *cgMode = 1;
  return(TRUE);
}


/*
    PsychCheckVideoSettings()
    
    Check all available video display modes for the specified screen number and return true if the 
    settings are valid and false otherwise.
*/
psych_bool PsychCheckVideoSettings(PsychScreenSettingsType *setting)
{
        CFDictionaryRef cgMode;       
        return(PsychGetCGModeFromVideoSetting(&cgMode, setting));
}

/*
    PsychGetScreenDepth()
    
    The caller must allocate and initialize the depth struct. 
*/
void PsychGetScreenDepth(int screenNumber, PsychDepthType *depth)
{    
    if(screenNumber>=numDisplays)
        PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range"); //also checked within SCREENPixelSizes
    PsychAddValueToDepthStruct((int) GetDeviceCaps(displayCGIDs[screenNumber], BITSPIXEL), depth);
}

int PsychGetScreenDepthValue(int screenNumber)
{
    PsychDepthType	depthStruct;
    
    PsychInitDepthStruct(&depthStruct);
    PsychGetScreenDepth(screenNumber, &depthStruct);
    return(PsychGetValueFromDepthStruct(0,&depthStruct));
}


int PsychGetNominalFramerate(int screenNumber)
{
  if (PsychPrefStateGet_ConserveVRAM() & kPsychIgnoreNominalFramerate) return(0);

  if(screenNumber>=numDisplays)
    PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetScreenDepths() is out of range"); 
  return(GetDeviceCaps(displayCGIDs[screenNumber], VREFRESH));
}

/* Returns the physical display size as reported by Windows: */
void PsychGetDisplaySize(int screenNumber, int *width, int *height)
{
    if(screenNumber>=numDisplays)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetDisplaySize() is out of range");
    *width = (int) GetDeviceCaps(displayCGIDs[screenNumber], HORZSIZE);
    *height = (int) GetDeviceCaps(displayCGIDs[screenNumber], VERTSIZE); 
}

void PsychGetScreenSize(int screenNumber, long *width, long *height)
{
  if(screenNumber>=numDisplays) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetScreenDepths() is out of range"); 
  *width = GetDeviceCaps(displayCGIDs[screenNumber], HORZRES);
  *height = GetDeviceCaps(displayCGIDs[screenNumber], VERTRES);
  // If we are in multi-display mode and size of screen 0 (our virtual full desktop) is
  // requested, then we return a size which is twice the desktop width, thereby providing
  // the proper dimensions for a full desktop spanning onscreen window on screen 0.
  if (numDisplays>2 && screenNumber == 0) *width = *width * 2;
}


void PsychGetGlobalScreenRect(int screenNumber, double *rect)
{
    long width, height; 
	
    PsychGetScreenSize(screenNumber, &width, &height);
    rect[kPsychLeft]= (int) displayDeviceStartX[screenNumber];
    rect[kPsychTop]= (int) displayDeviceStartY[screenNumber];
    rect[kPsychRight]=rect[kPsychLeft] + (int) width;
    rect[kPsychBottom]=rect[kPsychTop] + (int) height;	
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
    Its probably better to read this directly from the CG renderer info than to infer it from the pixel size
*/	
int PsychGetNumScreenPlanes(int screenNumber)
{    
    return((PsychGetScreenDepthValue(screenNumber)>24) ? 4 : 3);
}



/*
	This is a place holder for a function which uncovers the number of dacbits.  To be filled in at a later date.
	If you know that your card supports >8 then you can fill that in the PsychtPreferences and the psychtoolbox
	will act accordingly.
	
	There seems to be no way to uncover the dacbits programatically.  According to apple CoreGraphics
	sends a 16-bit word and the driver throws out whatever it chooses not to use.
		
	For now we just use 8 to avoid false precision.  
	
	If we can uncover the video card model then  we can implement a table lookup of video card model to number of dacbits.  
*/
int PsychGetDacBitsFromDisplay(int screenNumber)
{
	return(8);
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
    //settings->dacbits=PsychGetDacBits(screenNumber);
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
    CFDictionaryRef     cgMode;
    psych_bool          isValid, isCaptured;

    if(settings->screenNumber>=numDisplays)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychSetScreenSettings() is out of range");

    //Check for a lock which means onscreen or offscreen windows tied to this screen are currently open.
    // MK: Disabled : if(PsychCheckScreenSettingsLock(settings->screenNumber)) return(false);  //calling function should issue an error for attempt to change display settings while windows were open.
    
    //store the original display mode if this is the first time we have called this function.  The psychtoolbox will disregard changes in 
    //the screen state made through the control panel after the Psychtoolbox was launched. That is, OpenWindow will by default continue to 
    //open windows with finder settings which were in place at the first call of OpenWindow.  That's not intuitive, but not much of a problem
    //either. 
    if(!displayOriginalCGSettingsValid[settings->screenNumber]){
      displayOriginalCGSettings[settings->screenNumber]= 1; // FIXME!!! CGDisplayCurrentMode(displayCGIDs[settings->screenNumber]);
      displayOriginalCGSettingsValid[settings->screenNumber]=TRUE;
    }
    
    //Find core graphics video settings which correspond to settings as specified withing by an abstracted psychsettings structure.  
    isValid=PsychGetCGModeFromVideoSetting(&cgMode, settings);
    if(!isValid){
        PsychErrorExitMsg(PsychError_internal, "Attempt to set invalid video settings"); 
        //this is an internal error because the caller is expected to check first. 
    }
    
    //If the caller passed cache settings (then it is SCREENResolutions) and we should cache the current video mode settings for this display.  These
    //are cached in the form of CoreGraphics settings and not Psychtoolbox video settings.  The only caller which should pass a set cache flag is 
    //SCREENResolutions
    if(cacheSettings){
        displayOverlayedCGSettings[settings->screenNumber]=cgMode;
        displayOverlayedCGSettingsValid[settings->screenNumber]=TRUE;
    }
    
    //Check to make sure that this display is captured, which OpenWindow should have done.  If it has not been done, then exit with an error.  
    isCaptured=PsychIsScreenCaptured(settings->screenNumber);
    if(!isCaptured) PsychErrorExitMsg(PsychError_internal, "Attempt to change video settings without capturing the display");
        
    // Change the display mode.   
	// We do this in PsychOSOpenWindow() if necessary for fullscreen-mode, but without changing any settings except
	// switch to fullscreen mode. Here we call PsychOSOpenWindow's helper routine to really change video settings:
	return(ChangeScreenResolution(settings->screenNumber, (int) PsychGetWidthFromRect(settings->rect), (int) PsychGetHeightFromRect(settings->rect), PsychGetValueFromDepthStruct(0,&(settings->depth)), settings->nominalFrameRate));
}

/*
    PsychRestoreVideoSettings()
    
    Restores video settings to the state set by the finder.  Returns TRUE if the settings can be restored or false if they 
    can not be restored because a lock is in effect, which would mean that there are still open windows.    
    
*/
psych_bool PsychRestoreScreenSettings(int screenNumber)
{
    psych_bool 			isCaptured;
    CGDisplayErr 		error=0;


    if(screenNumber>=numDisplays)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetScreenDepths() is out of range"); //also checked within SCREENPixelSizes

    //Check for a lock which means onscreen or offscreen windows tied to this screen are currently open.
    // Disabled: if(PsychCheckScreenSettingsLock(screenNumber)) return(false);  //calling function will issue error for attempt to change display settings while windows were open.
    
    //Check to make sure that the original graphics settings were cached.  If not, it means that the settings were never changed, so we can just
    //return true. 
    if(!displayOriginalCGSettingsValid[screenNumber])
        return(true);
    
    //Check to make sure that this display is captured, which OpenWindow should have done.  If it has not been done, then exit with an error.  
    isCaptured=PsychIsScreenCaptured(screenNumber);
    if(!isCaptured) PsychErrorExitMsg(PsychError_internal, "Attempt to change video settings without capturing the display");
    
    // Restore video settings from the defaults in the Windows registry:
	if (PsychGetNumDisplays() > 2 && screenNumber == 0) {
		// Special case: screenNumber 0 on a dual-display setup or multi-display setup. This is a "virtual" screen,
		// which means the displays for screen 1 and 2. Therefore we need to restore both heads 1 and 2 settings:
		ChangeDisplaySettingsEx(PsychGetDisplayDeviceName(1), NULL, NULL, 0, NULL);
		ChangeDisplaySettingsEx(PsychGetDisplayDeviceName(2), NULL, NULL, 0, NULL);
	}
	else {
		// Standard case: Just reset given screen:
		ChangeDisplaySettingsEx(PsychGetDisplayDeviceName(screenNumber), NULL, NULL, 0, NULL);
	}

    return(true);
}


void PsychHideCursor(int screenNumber, int deviceIdx)
{
  // Hide the mouse cursor: We ignore the screenNumber as Windows
  // doesn't allow to set the cursor per screen anyway.
  while(ShowCursor(FALSE)>=0);
  return;
}

void PsychShowCursor(int screenNumber, int deviceIdx)
{
  // Show the mouse cursor: We ignore the screenNumber as Windows
  // doesn't allow to set the cursor per screen anyway.
  while(ShowCursor(TRUE)<0);
}

void PsychPositionCursor(int screenNumber, int x, int y, int deviceIdx)
{
  // Reposition the mouse cursor: We ignore the screenNumber as Windows
  // doesn't allow to set the cursor per screen anyway.
  if (SetCursorPos(x,y) == 0) PsychErrorExitMsg(PsychError_internal, "Couldn't position the mouse cursor! (SetCursorPos() failed).");
}

/*
    PsychReadNormalizedGammaTable()
    
    TO DO: This should probably be changed so that the caller allocates the memory.
    TO DO: Adopt a naming convention which distinguishes between functions which allocate memory for return variables from those which do not.
            For example, PsychReadNormalizedGammaTable() vs. PsychGetNormalizedGammaTable().
    
*/
void PsychReadNormalizedGammaTable(int screenNumber, int outputId, int *numEntries, float **redTable, float **greenTable, float **blueTable)
{
    CGDirectDisplayID	cgDisplayID;
    static  float localRed[256], localGreen[256], localBlue[256];

	 // Windows hardware LUT has 3 tables for R,G,B, 256 slots each, concatenated to one table.
    // Each entry is a 16-bit word with the n most significant bits used for an n-bit DAC.
	 psych_uint16	gammaTable[256 * 3]; 
    psych_bool    ok;
	 int     i;        

    // Query OS for gamma table:
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);
    ok=GetDeviceGammaRamp(cgDisplayID, &gammaTable);
	 if (!ok) PsychErrorExitMsg(PsychError_internal, "Failed to query the hardware gamma table from graphics adapter!");

	 // Convert concatenated table into three separate tables, map 16-bit values into
    // 0-1 normalized floats
    *redTable=localRed; *greenTable=localGreen; *blueTable=localBlue;
    for (i=0; i<256; i++) localRed[i]   = ((float) gammaTable[i]) / 65535.0f;
    for (i=0; i<256; i++) localGreen[i] = ((float) gammaTable[i+256]) / 65535.0f;
    for (i=0; i<256; i++) localBlue[i]  = ((float) gammaTable[i+512]) / 65535.0f;

	 // The LUT's always have 256 slots for the 8-bit framebuffer:
    *numEntries= 256;
}

unsigned int PsychLoadNormalizedGammaTable(int screenNumber, int outputId, int numEntries, float *redTable, float *greenTable, float *blueTable)
{
    psych_bool 	ok; 
    CGDirectDisplayID	cgDisplayID;
	 int     i;        
	 // Windows hardware LUT has 3 tables for R,G,B, 256 slots each, concatenated to one table.
     // Each entry is a 16-bit word with the n most significant bits used for an n-bit DAC.
	 psych_uint16	gammaTable[256 * 3]; 

     // Special case empty 0-slot table provided? That means to load an identity
     // gamma table and setup the GPU for identity pass-through from framebuffer to
     // encoders. This is unsupported on Windows, so return the 0xffffffff "unsupported" code.
     if (numEntries == 0) return(0xffffffff);
     
	 // Table must have 256 slots!
	 if (numEntries!=256) PsychErrorExitMsg(PsychError_user, "Loadable hardware gamma tables must have 256 slots!");    

	 // Convert input table to Windows specific gammaTable:
    for (i=0; i<256; i++) gammaTable[i]     = (int)(redTable[i]   * 65535.0f + 0.5f);
    for (i=0; i<256; i++) gammaTable[i+256] = (int)(greenTable[i] * 65535.0f + 0.5f);
    for (i=0; i<256; i++) gammaTable[i+512] = (int)(blueTable[i]  * 65535.0f + 0.5f);

	 // Set new gammaTable: On M$-Windows, we retry up to 10 times before giving up, because some
	 // buggy Windows graphics drivers seem to fail on first invocation of SetDeviceGammaRamp(), just
	 // to succeed on a 2nd invocation!
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);
	 ok=FALSE;
    for (i=0; i<10 && !ok; i++) ok=SetDeviceGammaRamp(cgDisplayID, &gammaTable);
	 if (!ok) PsychErrorExitMsg(PsychError_user, "Failed to upload the hardware gamma table into graphics adapter! Read the help for explanation...");

    // Return "success":
     return(1);
}

// Beamposition queries on Windows are implemented via the DirectDraw-7 interface. It provides
// the GetScanLine() method for beamposition queries, and also functions like WaitForVerticalBlank()
// and VerticalBlankStatus(). This is supposed to work on Windows 2000 and later.
// See http://msdn.microsoft.com/archive/default.asp?url=/archive/en-us/ddraw7/directdraw7/ddref_2n5j.asp
//
// This seems to work on all recent NVidia, ATI and Intel gfx-chips.
//
// On MS-Windows, the display specifier (cgDisplayId, screenNumber) is ignored! We always query
// the "default" pipeline due to current lack of support for mapping display specs to pipes.
// The choice of "default" pipeline is OS/Driver dependent and depends on the display settings
// and choices made by the user in the Windows display settings panel, usually the "primary display",
// which is reasonable (always correct on single display setup, user selectable on dual-display setup
// by assignment of primary monitor, irrelevant in dual-display desktop spanning mode for binocular
// stimulation, as both display pipes are synchronized in that case). 
//
// PsychGetDisplayBeampPosition() contains the implementation of display beamposition queries.
// It requires both, a cgDisplayID handle, and a logical screenNumber and uses one of both for
// deciding which display pipe to query, whatever of both is more efficient or suitable for the
// host platform -- This is ugly, but neccessary, because the mapping with only one of these
// specifiers would be either ambigous (wrong results!) or usage would be inefficient and slow
// (bad for such a time critical low level call!). On some systems it may even ignore the arguments,
// because it's not capable of querying different pipes - ie., it will always query a hard-coded pipe.
//
int PsychGetDisplayBeamPosition(CGDirectDisplayID cgDisplayId, int screenNumber)
{
	HRESULT rc;
	int vblbias, vbltotal;
	psych_uint32 ubeampos = 0;
	int	beampos = -1;
	
	// Apply remapping of screenId's to heads, if any: Usually identity mapping.
	screenNumber = PsychScreenToHead(screenNumber, 0);
	
	if(displayDeviceDDrawObject[screenNumber]) {
		// We have a Direct draw object: Try to use GetScanLine():
		if (enableVBLBeamposWorkaround) {
			// Beamposition queries don't work within vertical blank interval. We need to
			// spin-wait in a busy waiting loop as long as the query status reports
			// DDERR_VERTICALBLANKINPROGRESS aka "display in retrace state":
			while((rc=IDirectDraw_GetScanLine(displayDeviceDDrawObject[screenNumber], (LPDWORD) &ubeampos)) == DDERR_VERTICALBLANKINPROGRESS);
		}
		else {
			// No known problems with query in VBL. Do a one-time query:
			rc=IDirectDraw_GetScanLine(displayDeviceDDrawObject[screenNumber], (LPDWORD) &ubeampos);
		}
		
		// Mask returned ubeampos with 0xffff: This will not allow any
		// beamposition greater than 65535 scanlines to be returned.
		// A reasonable limit for the foreseeable future:
		ubeampos &= 0xffff;
		
		beampos = (int) ubeampos;
		
		// Apply corrective offsets if any (i.e., if non-zero):
		PsychGetBeamposCorrection(screenNumber, &vblbias, &vbltotal);
		beampos = beampos - vblbias;
		if (beampos < 0) beampos = vbltotal + beampos;
		
		// Valid result? If so, return it. Otherwise fall-through to error return...
		if (rc==DD_OK || rc==DDERR_VERTICALBLANKINPROGRESS) return(beampos);
	}
	
	// Direct Draw unavailable or function unsupported, or hardware
	// doesn't support query under given configuration or error condition:
	// We return -1 as an indicator to high-level routines that we don't
	// know the rasterbeam position.
	return(-1);
}

// Do not have kernel driver support for MS-Windows:
psych_bool PsychOSIsKernelDriverAvailable(int screenId)
{
	return(FALSE);
}

// Internal helper routine, called from PsychCaptureScreen(), which in turn is called from
// SCREENOpenWindow(). Give the DDRAW beamposition query routines a serious workout to find
// some possible bugs. Enable proper workarounds if bugs should be encountered:
void PsychTestDDrawBeampositionQueries(int screenNumber)
{
	int w1, w2, h1, h2, vbldetectcount, bogusvaluecount, bogusvalueinvblcount, totalcount;
	psych_uint32 maxvpos;
	psych_uint32 beampos = 100000;
	HRESULT rc;
	double tdeadline, now;
	int verbosity = PsychPrefStateGet_Verbosity();
	
	// Check how beamposition query behaves inside the vertical blanking interval:
	if((displayDeviceDDrawObject[PsychScreenToHead(screenNumber, 0)]) && (PsychPrefStateGet_VBLTimestampingMode()>=0)) {
		// We have a Direct draw object and beampos queries are enabled: Try to test GetScanLine():
		
		// First find reference height values for display, aka start of vertical blank.
		// Valid values inside the vertical blank area should be between h1 and 1.25*h1
		// or between h2 and 1.25*h2, i.e., we assume the VBL doesn't extend more than
		// 25% beyond the active areas height, so we should not see any values outside
		// that area if DDERR_VERTICALBLANKINPROGRESS is reported. If we see higher values,
		// or lower values, then that likely means that the gfx-driver doesn't report meaningful
		// beampos values during VBL, but does overwrite our default 0xdeadbeef canary with some
		// random trash.
		
		// Do all screens have their own individual beampos query DDRAW object assigned?
		if ((ddrawnumDisplays == numDisplays) || (numDisplays == 1)) {
			// Yes. Our reference maxvpos is the screen height of the target screen 'screenNumber':
			PsychGetScreenSize(screenNumber, &w1, &h1);
			maxvpos = (psych_uint32) h1;
		}
		else {
			// No. A common beampos query object for all displays. We choose the bigger vertical
			// height of the first two displays:
			PsychGetScreenSize(1, &w1, &h1);
			PsychGetScreenSize(2, &w2, &h2);
			maxvpos = (psych_uint32) ((h1 > h2) ? h1 : h2);
		}
		
		// maxvpos is maximum valid scanline:
		maxvpos = (psych_uint32) (1.25 * ((float) maxvpos));
		
		// We measure for 500 msecs -- That is at least 30 video refresh cycles:
		vbldetectcount = 0; 
		bogusvaluecount = 0;
		bogusvalueinvblcount = 0;
		totalcount = 0;
		
		PsychGetAdjustedPrecisionTimerSeconds(&tdeadline);
		now = tdeadline;
		tdeadline+=0.500;
		
		while (now < tdeadline) {
			// Update timestamp:
			PsychGetAdjustedPrecisionTimerSeconds(&now);
			
			// Query beam position:
			beampos = 0xdeadbeef;
			rc=IDirectDraw_GetScanLine(displayDeviceDDrawObject[PsychScreenToHead(screenNumber, 0)], (LPDWORD) &beampos);
			if (rc==DD_OK || rc==DDERR_VERTICALBLANKINPROGRESS) {
				// Some sample returned...
				totalcount++;
				
				// Bogus value, ie., beampos not updated at all or out of range?
				if ((beampos != 0xdeadbeef) && ((beampos < 0) || (beampos > maxvpos))) beampos = 0xdead;
				
				// Update total count of bogus values:
				if (beampos == 0xdeadbeef || beampos == 0xdead) bogusvaluecount++;
				
				if (rc == DDERR_VERTICALBLANKINPROGRESS) {
					// Update count of VBLs:
					vbldetectcount++;
					
					// Update total count of bogus values while VBL active:
					if (beampos == 0xdeadbeef || beampos == 0xdead) bogusvalueinvblcount++;
				}
			}
		}
		
		// Ok, run the numbers:
		if (bogusvaluecount > 0) {
			// Some failed queries:
			if (verbosity > 1) {
				printf("PTB-WARNING: Startup test of beamposition queries for high precision timestamping detected problems on your graphics card + driver combo.\n");
				if (verbosity > 3) {
					printf("PTB-WARNING: During initial test of a total of %i queries, a total of %i beamposition queries reported bogus values!\n", totalcount, bogusvaluecount);
					printf("PTB-WARNING: During a total of %i queries within blanking interval, %i reported wrong values.\n", vbldetectcount, bogusvalueinvblcount);
				}
			}
			
			// All queries failed in VBL?
			if (bogusvalueinvblcount == bogusvaluecount) {
				// Yes. This allows for an ugly but workable workaround:
				if (verbosity > 2) {
					printf("PTB-INFO: Seems beamposition queries are unsupported whenever the display is in a vertical retrace state.\n");
					printf("PTB-INFO: I will enable a work-around which should satisfy the timing needs of most applications, ie.\n");
					printf("PTB-INFO: robust and jitter-free deterministic timestamps. However, all reported timestamps may have a\n");
					printf("PTB-INFO: constant offset of up to the duration of the vertical blanking interval (max 1 msec worst case)\n");
					printf("PTB-INFO: of your display. Read 'help BeampositionQueries' on how you can remove that bias by some manual\n");
					printf("PTB-INFO: intervention, should you need that last bit of accuracy. A second drawback is increased cpu load,\n");
					printf("PTB-INFO: which may cause louder air-fan noise due to heat production and lower battery runtime of laptops.\n");
				}
				
				// Enable beampos workaround by default:
				enableVBLBeamposWorkaround = TRUE;
			}
			else {
				// No. Some bogus values also happen outside VBL. In this case, our workaround won't work
				// and this problem is unfixable for now. The only safe and sane solution is to totally
				// disable beamposition based timestamping:
				
				// Disable beampos workaround:
				enableVBLBeamposWorkaround = FALSE;

				// Disable beampos queries:
				PsychPrefStateSet_VBLTimestampingMode(-1);
				
				if (verbosity > 1) {
					printf("PTB-WARNING: Some of the queries fail even outside the vertical blank interval, so no effective workaround for this driver bug exists.\n");
					printf("PTB-WARNING: This renders beamposition queries pretty useless -- Disabling high precision timestamping for now.\n");
					printf("PTB-WARNING: Please report this message with a description of your graphics card, operating system and video driver to\n");
					printf("PTB-WARNING: the Psychtoolbox forum. Maybe the gathered information allows for some work-around in\n");
					printf("PTB-WARNING: future PTB releases to get high precision timestamping back on your setup.\n\n");
				}
			}
		}
		else {
			// No apparent problems detected. Should the workaround be forced on?
			if (PsychPrefStateGet_ConserveVRAM() & kPsychUseBeampositionQueryWorkaround) {
				// Force it on:
				if (verbosity > 2) {
					printf("PTB-INFO: Automatic startup test of beamposition queries couldn't detect any beamposition problems on your system,\n");
					printf("PTB-INFO: but by setting the special flag Screen('Preference', 'ConserveVRAM', 4096); somewhere in your script,\n");
					printf("PTB-INFO: you requested that the workaround for beamposition problems is always unconditionally enabled.\n");
					printf("PTB-INFO: Therefore i will enable the work-around: This should satisfy the timing needs of most applications, ie.\n");
					printf("PTB-INFO: robust and jitter-free deterministic timestamps. However, all reported timestamps may have a\n");
					printf("PTB-INFO: constant offset of up to the duration of the vertical blanking interval (max 1 msec worst case)\n");
					printf("PTB-INFO: of your display. Read 'help BeampositionQueries' on how you can remove that bias by some manual\n");
					printf("PTB-INFO: intervention, should you need that last bit of accuracy. A second drawback is increased cpu load,\n");
					printf("PTB-INFO: which may cause louder air-fan noise due to heat production and lower battery runtime of laptops.\n");
				}
				
				// Enable beampos workaround by default:
				enableVBLBeamposWorkaround = TRUE;				
			}
		}
	} // Beamposition queries supported and requested.
	
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
  // A no-op on Windows:
  return(FALSE);

  /*
  if (!PsychOSIsKernelDriverAvailable(screenNumber)) return(FALSE);

  if (gpuMaintype) *gpuMaintype = fDeviceType;
  if (gpuMinortype) *gpuMinortype = fCardType;
  if (pciDeviceId) *pciDeviceId = fPCIDeviceId;
  if (numDisplayHeads) *numDisplayHeads = fNumDisplayHeads;

  return(TRUE);
  */
}
