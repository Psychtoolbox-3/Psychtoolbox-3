/*
	ScreenPreferenceState.c		
  
	AUTHORS:

    Allen.Ingling@nyu.edu           awi
    mario.kleiner@tuebingen.mpg.de  mk
  
	PLATFORMS:	
	
		All.
    

	HISTORY:

		2/28/04  awi		Created.   
		5/30/05  mk           New preference setting screenSkipSyncTests: If set to TRUE, self-test and calibration will be as
							short and sloppy as possible and some problematic conditions will not cause an abort.
							One can trade-off correctness and accuracy of stimulus presentation against level of annoyance
							during development/debugging of new experiments.
		9/30/05  mk           new setting VisualDebugLevel: Defines how much visual feedback PTB should give about errors and
							state: 0=none, 1=only errors, 2=also warnings, 3=also infos, 4=also blue bootup screen, 5=also visual test sheets.
		3/7/06   awi        Added state for new preference flag SuppressAllWarnings. 
  
	DESCRIPTION:
  
		This file holds state for the Screen Preference command.
  
   
  
*/


#include "Screen.h"

// Default textrenderer already defined? E.g., Compile time options?
#ifndef PTB_DEFAULT_TEXTRENDERER
// Nope: Default to the high-quality OS specific renderer.
// This setting of 1 maps to GDI renderer on Windows, ATSU on OS/X and the FTGL plugin on Linux, if available.
// A setting of zero would map to ATSU on OS/X and display list renderers on Windows and Linux, as they are faster.
// A setting of 2 would map to the FTGL plugin on Linux, OS/X and Windows, if available.
#define PTB_DEFAULT_TEXTRENDERER 1
#endif

// Default video capture engine is GStreamer (3):
#define PTB_DEFAULTVIDCAPENGINE 3

#if PSYCH_SYSTEM == PSYCH_LINUX
#define INITIAL_DEFAULT_FONT_NAME		"Times"
#define INITIAL_DEFAULT_FONT_SIZE		24
#define INITIAL_DEFAULT_FONT_STYLE		0
#endif

#if PSYCH_SYSTEM == PSYCH_OSX
#define INITIAL_DEFAULT_FONT_NAME		"Geneva"
#define INITIAL_DEFAULT_FONT_SIZE		12
#define INITIAL_DEFAULT_FONT_STYLE		0
#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS
#define INITIAL_DEFAULT_FONT_NAME		"Courier New"
#define INITIAL_DEFAULT_FONT_SIZE		18
#define INITIAL_DEFAULT_FONT_STYLE		0
#endif

//PsychTable preference state
static int								psychTableVersion;				//there is no psych table yet, this is provided for the future. 
static char								PsychTableCreator[]="Screen";   //there is no psych table yet, this is provided for the future.
//Text preference state
static int								defaultTextYPositionIsBaseline; // Use new style of text positioning by default: y-pos is top of text.
#define MAX_DEFAULT_FONT_NAME_LENGTH    256
static char								defaultFontName[MAX_DEFAULT_FONT_NAME_LENGTH];
static int								defaultTextSize;
static int								defaultTextStyle;				// 0=normal,1=bold,2=italic,4=underline,8=outline,32=condense,64=extend
static psych_bool                       textAlphaBlending;
static int								textAntiAliasing;				// -1=System defined (don't care), 0=Always off, 1=Always on.
static int								textRenderer;					// 0=Default OS specific (fast one), 1=OS specific High quality.
static int                              screenSkipSyncTests;			// 0=Do full synctests, abort on failure, 1=Reduced tests, continue with warning, 2=Skip'em
//Debug preference state
static psych_bool						TimeMakeTextureFlag;
static int								screenVisualDebugLevel;
static int                              screenConserveVRAM;
// If EmulateOldPTB is set to true, then try to behave like the old OS-9 PTB:
static psych_bool                       EmulateOldPTB;
// Support for real 3D rendering enabled? Any non-zero value enables 3D rendering, a setting of 1 with defaults, values > 1 enable additional features. Disabled by default.
static int								Enable_3d_gfx;
static int                              screenVBLTimestampingMode;
static int								screenVBLEndlineOverride;		// Manual override for VTOTAL - Endline of VBL. -1 means "auto-detect" this is the default.
static int								videoCaptureEngineId;			// Default video capture engine: 0 = Quicktime, 1 = LibDC1394 Firewire, 2 = ARVideo.
static int								windowShieldingLevel;			// Level of priority of windowed onscreen window wrt. other windows:
																		// From 0 for "behind everything" to 2000 for "in front of everything. Exact meaning of
																		// number is OS specific. This value is used at window open time for each window.
static double							frameRectLadderCorrection;		// Tweak factor to apply in SCREENFrameRect.c for different GPU's.
static psych_bool						suppressAllWarnings;

// General level of verbosity:
// 0 = Shut up.
// 1 = Errors only.
// 2 = Errors and Warnings.
// 3 = Above + Basic Info, e.g., window properties and copyright. This is the default.
// 4 = Lot's more debug output and infos.
// 5 = Extreme debug output.
static int								Verbosity;

// Thresholds used during VBL sync tests and refresh rate measurement and calibration
// in PsychOpenOnscreenWindow()
static double							sync_maxStddev;					// Maximum standard deviation from mean duration in seconds.
static double							sync_maxDeviation;				// Maximum deviation (in percent) between measured and OS reported reference frame duration.
static double							sync_maxDuration;				// Maximum duration of a calibration run in seconds.
static int								sync_minSamples;				// Minimum number of valid measurement samples needed.

static int                              useGStreamer;                   // Use GStreamer for multi-media processing? 1==yes.

//All state checking goes through accessors located in this file.

// Called by Screen init code first: Sets up all default values after a
// (re-)load of the MEX file, ie., at first invocation, or after a clear all, clear mex etc. 
void PrepareScreenPreferences(void)
{
	//set the fake psych table version
	psychTableVersion=20;
	sprintf(PsychTableCreator, "Screen");
	defaultTextYPositionIsBaseline=0;
	defaultTextSize=INITIAL_DEFAULT_FONT_SIZE;
	defaultTextStyle=INITIAL_DEFAULT_FONT_STYLE;
	textAlphaBlending=FALSE;
	textAntiAliasing=-1;
	textRenderer=PTB_DEFAULT_TEXTRENDERER;
	screenSkipSyncTests=0;
	TimeMakeTextureFlag=FALSE;
	screenVisualDebugLevel=4;
	screenConserveVRAM=0;
	EmulateOldPTB=FALSE;
	Enable_3d_gfx=0;
	screenVBLTimestampingMode = ((PSYCH_SYSTEM == PSYCH_LINUX) || (PSYCH_SYSTEM == PSYCH_WINDOWS)) ? 4 : 1;
	screenVBLEndlineOverride=-1;
	videoCaptureEngineId=PTB_DEFAULTVIDCAPENGINE;
	windowShieldingLevel=2000;
	frameRectLadderCorrection=-1.0;
	suppressAllWarnings=FALSE;

	// Default level of verbosity is 3:
    Verbosity=3;

    // Early override via environment variable, if defined:
    if (getenv("PSYCH_SCREEN_VERBOSITY")) {
        Verbosity = atoi(getenv("PSYCH_SCREEN_VERBOSITY"));
    }

	// Default synctest settings: 200 usec allowable max standard deviation from measured
	// mean flip duration, at least 50 valid sync samples, at most 10% deviation between
	// measured duration and reference duration (os reported or other), at most 5 seconds
	// worst-case duration per calibration run:
	PsychPrefStateSet_SynctestThresholds(0.000200, 50, 0.1, 5);
	
	// Initialize our locale setting for multibyte/singlebyte to unicode character conversion
	// for Screen('DrawText') et al. to be the current default system locale, as defined by
	// system settings and environment variables at startup of Matlab/Octave:
	// N.B. This function is special as the affected state setting and routines are not defined
	// here, but inside the Screen text handling routines, currently in SCREENDrawText.c
	PsychSetUnicodeTextConversionLocale("");
	PsychPrefStateSet_DefaultFontName(INITIAL_DEFAULT_FONT_NAME);

	// Don't use GStreamer by default on 32-bit builds for OS/X and Windows:
	PsychPrefStateSet_UseGStreamer(0);

	return;
}

//Accessors


//****************************************************************************************************************
//PsychTable preferences

/*
preference: PsychTableVersion
*/
int PsychPrefStateGet_PsychTableVersion(void)
{
	return(psychTableVersion);
}

/*
preference: PsychTableCreator
*/
const char *PsychPrefStateGet_PsychTableCreator(void)
{
	return(PsychTableCreator);
}

//****************************************************************************************************************
//Text and Font preferences

// If set to true, then the y-position specified in Screen('DrawText') defines
// the baseline of the text, not the top of the text. Default is off -> top of text.
// The default changes to 'on' -> Baseline if emulation of old PTB enabled on OS-X.
// This setting can be overriden case-by-case with the optional 7th argument to 'DrawText':
int PsychPrefStateGet_TextYPositionIsBaseline(void)
{
	return(defaultTextYPositionIsBaseline);
}

void PsychPrefStateSet_TextYPositionIsBaseline(int textPosIsBaseline)
{
	defaultTextYPositionIsBaseline = (textPosIsBaseline > 0) ? 1 : 0;
}

int PsychPrefStateGet_TextAntiAliasing(void)
{
	return(textAntiAliasing);
}

void PsychPrefStateSet_TextAntiAliasing(int mode)
{
	textAntiAliasing = mode;
}

int PsychPrefStateGet_TextRenderer(void)
{
	return(textRenderer);
}

void PsychPrefStateSet_TextRenderer(int mode)
{
	textRenderer = mode;
}

/*
preference: DefaultFontName
*/
void PsychPrefStateGet_DefaultFontName(const char **fontName )
{
	*fontName=defaultFontName;
}


void PsychPrefStateSet_DefaultFontName(const char *newName)
{
	if(strlen(newName)+1 > MAX_DEFAULT_FONT_NAME_LENGTH)
		PsychErrorExitMsg(PsychError_user, "Attempt to set a default font name using a string >255 characters");
	strcpy(defaultFontName, newName);
}

/*
preference: DefaultTextSize
*/
int PsychPrefStateGet_DefaultTextSize(void)
{
	return(defaultTextSize);
}

void PsychPrefStateSet_DefaultTextSize(int textSize)
{
	defaultTextSize=textSize;
}

/*
preference: DefaultTextStyle
*/
int PsychPrefStateGet_DefaultTextStyle(void)
{
	return(defaultTextStyle);
}

void PsychPrefStateSet_DefaultTextStyle(int textStyle)
{
	defaultTextStyle=textStyle;
}

/*
preference: TextAlphaBlending
*/
psych_bool PsychPrefStateGet_TextAlphaBlending(void)
{
	return(textAlphaBlending);
}

void PsychPrefStateSet_TextAlphaBlending(psych_bool enableFlag)
{
	textAlphaBlending=enableFlag;
}

// Screen self-test and calibration preferences for syncing to VBL and such...  
int PsychPrefStateGet_SkipSyncTests(void)
{
    return(screenSkipSyncTests);
}

void PsychPrefStateSet_SkipSyncTests(int level)
{
    screenSkipSyncTests = level;
}

// Screen debug level for output of visual feedback signals...  
int PsychPrefStateGet_VisualDebugLevel(void)
{
    return(screenVisualDebugLevel);
}

void PsychPrefStateSet_VisualDebugLevel(int level)
{
    screenVisualDebugLevel=level;
}

// Get and set mode of operation for Screen('Flip')'s VBL
// and bufferswap timestamping. This parameter is better left
// at default setting by most users. Useful for debugging/testing
// and in special cases. Currently meaningless on Linux
// and therefore silently ignored.
//
// -1= Always use uncorrected timestamps.
// 0 = low-level, kernel-based timestamping always off: Use beampos method if available, otherwise use uncorrected timestamps.
// 1 = Automatic, on demand: Use beampos method if available, use the kernel-level method if beampos method unsupported, e.g., IntelMacs.
// 2 = Always use kernel-level method, but only as a consistency check for beampos method -- For the super-paranoid and for testing.
// 3 = Always use kernel-level method, this method overrides everything else. This is the opposite of 1.
// 4 = Use operating system specific timestamping method: On Linux this is OpenML OML_sync_control, on Windows this is DwmGetCompositionTimingInfo().
//
int PsychPrefStateGet_VBLTimestampingMode(void)
{
    return(screenVBLTimestampingMode);
}

void PsychPrefStateSet_VBLTimestampingMode(int level)
{
    screenVBLTimestampingMode = level;
}

// Override for display endline aka VTOTAL:
int PsychPrefStateGet_VBLEndlineOverride(void)
{
    return(screenVBLEndlineOverride);
}

void PsychPrefStateSet_VBLEndlineOverride(int level)
{
    screenVBLEndlineOverride = level;
}

// Settings for conserving VRAM usage by disabling certain features.
// Current constants (defined in ScreenTypes.h) that can be or'ed together:
// kPsychDisableAUXBuffers = Don't use AUX buffers, disable associated features.
// kPsychDontCacheTextures = Use Client storage for RAM->DMA texturing instead of VRAM texturing.
// Defaults to zero = Make full use of VRAM.
int PsychPrefStateGet_ConserveVRAM(void)
{
    return(screenConserveVRAM);
}

void PsychPrefStateSet_ConserveVRAM(int level)
{
    screenConserveVRAM = level;
}

// If EmulateOldPTB is set to true (default is false) at startup, then we try to
// behave like the old OS-9 and Windows Psychtoolboxes. We don't use double-buffering
// and SCREEN('Flip') for stimulus onset- offset and timing, but good'ol
// Screen('WaitBlanking') in conjunction with drawing commands that execute
// immediately on invocation (glFinish()) and a single-buffered onscreen window.
// Offscreen windows are also available in a compatible way.
//
// Note: Internally we still use a double-buffered context, but the front buffer is
// the drawing/reading target for all commands and the backbuffer is used as a
// scratchpad buffer for Offscreen window handling.
psych_bool PsychPrefStateGet_EmulateOldPTB(void)
{
    return(EmulateOldPTB);
}

void PsychPrefStateSet_EmulateOldPTB(psych_bool level)
{
    EmulateOldPTB = level;
	// When emulation for old PTB gets enabled, we change the default for
	// text baseline to 'on' -- The behaviour of old PTB.
	#if PSYCH_SYSTEM == PSYCH_OSX
		if (EmulateOldPTB>0) defaultTextYPositionIsBaseline=1;
	#endif
}

// Enable switch for 3D graphics support. If set to true, PTB will allocate stencil-
// and depth-buffers additionally to the AUX and Colorbuffers and perform additional
// bookkeeping to make sure we can do real 3D rendering and interface to external
// OpenGL mexfiles like, e.g., moglcore...
int PsychPrefStateGet_3DGfx(void)
{
    return(Enable_3d_gfx);
}

void PsychPrefStateSet_3DGfx(int level)
{
	if (level < 0) PsychErrorExitMsg(PsychError_invalidIntegerArg, "3D graphics preference setting must be a non-negative integer value! You passed a negative one.");
    Enable_3d_gfx = level;
}

// Default for selection of video capture engine:
int PsychPrefStateGet_VideoCaptureEngine(void)
{
	return(videoCaptureEngineId);
}

void PsychPrefStateSet_VideoCaptureEngine(int mode)
{
	videoCaptureEngineId = mode;
}

int PsychPrefStateGet_WindowShieldingLevel(void)
{
    return(windowShieldingLevel);
}

void PsychPrefStateSet_WindowShieldingLevel(int level)
{
    windowShieldingLevel = level;
}

// Correction tweak offset for proper Screen('FrameRect') behaviour:
void PsychPrefStateSet_FrameRectCorrection(double level)
{
	frameRectLadderCorrection = level;
}

double PsychPrefStateGet_FrameRectCorrection(void)
{
	return(frameRectLadderCorrection);
}

// Tweakable parameters for VBL sync tests and refresh rate calibration:
void PsychPrefStateSet_SynctestThresholds(double maxStddev, int minSamples, double maxDeviation, double maxDuration)
{
	if (maxStddev <= 0.0) PsychErrorExitMsg(PsychError_user, "Invalid (<= 0) maximum standard deviation provided!");
	if (maxDeviation <= 0.0 || maxDeviation > 1.0) PsychErrorExitMsg(PsychError_user, "Invalid (not in range 0.0 - 1.0) maximum deviation from reference framerate provided!");
	if (maxDuration <= 0.0) PsychErrorExitMsg(PsychError_user, "Invalid (<= 0) maximum duration of calibration procedure provided!");
	if (minSamples < 0) PsychErrorExitMsg(PsychError_user, "Invalid (negative) minimum number of samples provided!");
	sync_maxStddev = maxStddev;
	sync_maxDeviation = maxDeviation;
	sync_maxDuration = maxDuration;
	sync_minSamples = minSamples;
}

void PsychPrefStateGet_SynctestThresholds(double* maxStddev, int* minSamples, double* maxDeviation, double* maxDuration)
{
	*maxStddev    = sync_maxStddev;
	*maxDeviation = sync_maxDeviation; 
	*maxDuration  = sync_maxDuration; 
	*minSamples   = sync_minSamples; 
}

//****************************************************************************************************************
//Debug preferences

psych_bool PsychPrefStateGet_DebugMakeTexture(void)
{
	return(TimeMakeTextureFlag);
}

void PsychPrefStateSet_DebugMakeTexture(psych_bool setFlag)
{
	TimeMakeTextureFlag=setFlag;
}

psych_bool PsychPrefStateGet_SuppressAllWarnings(void)
{
	return(suppressAllWarnings);
}

void PsychPrefStateSet_SuppressAllWarnings(psych_bool setFlag)
{
	suppressAllWarnings=setFlag;
}

// Level of verbosity:
int PsychPrefStateGet_Verbosity(void)
{
	// Clamp level of verbosity to max. 1 if suppression of warnings is requested:
	if (suppressAllWarnings && Verbosity>1) return(1); 
	return(Verbosity);
}

void PsychPrefStateSet_Verbosity(int level)
{
	Verbosity = level;
}

void PsychPrefStateSet_UseGStreamer(int value)
{
	useGStreamer = value;
}

int PsychPrefStateGet_UseGStreamer(void)
{
	return(useGStreamer);
}

// Screen -> Head mappings: These are special, because the default
// mapping gets initialized during display initialization, and
// the actual mappings are stored in PsychGraphicsHardwareHALSupport.c,
// so these are just wrappers around the true accessor functions:
void PsychPrefStateSet_ScreenToHead(int screenId, int headId, int crtcId, int rankId)
{
	PsychSetScreenToHead(screenId, headId, rankId);
    PsychSetScreenToCrtcId(screenId, crtcId, rankId);
}

int PsychPrefStateGet_ScreenToHead(int screenId, int rankId)
{
	return(PsychScreenToHead(screenId, rankId));
}

int PsychPrefStateGet_ScreenToCrtcId(int screenId, int rankId)
{
	return(PsychScreenToCrtcId(screenId, rankId));
}
