/*
 *        PsychToolbox3/Source/Common/Screen/ScreenTypes.h
 *
 *        PLATFORMS: All
 *
 *        AUTHORS:
 *              Allen Ingling   awi     Allen.Ingling@nyu.edu
 *              Mario Kleiner   mk      mario.kleiner.de@gmail.com
 *
 *        HISTORY:
 *        09/06/02        awi        wrote it.
 *        11/26/02        awi        turned PsychDepthType into an int from enum.
 *        01/06/07        mk         Change data type for Color type from int to double to accomodate all the HDR stuff and
 *                                   get rid of most color conversion routines, which were useless and a bad idea from day 1.
 *
 *        DESCRIPTION:
 *
 *        Platform-neutral Typedefs for the Screen Psychtoolbox project.  This file should appear ahead of all
 *        other files on the Screen.h include list.
 *
 *        This file should contain only platform-neutral typedefs.  Platform specfic stuff all goes into
 *        files within PsychToolbox3/Source/PLATFORMNAME/Screen/.  Currently those files are PsychSurfaceGlue.c
 *        and PsychSurfaceGlue.h
 *
 *        TO DO:
 *
 *        -The depth stuff is stupid.  Fix it: break out the number of values and an array from that dumb structure and abolish the structure.
 *
 *        -PsychColorModeType seems mostly useless, it contains no more information than the depth alone.  If it is always the case that:
 *        8-bits = indexed mode
 *        16-bits = 5 bits of red , 5 bits of green 5 bits of blue.
 *        24-bits = 8 bits of red, 8 bits of green, 8 bits of blue.
 *        32-bits = 8 bits of red, 8 bits of green, 8 bits of blue.
 *        Then we don't need the mode as a seperate and purely redundant representation.  However, I think somewhere a DirectX structure defined an 8-bit direct mode,
 *        so we might want to consider verifying that before abolishing PsychColorModeType.  Still, it might be better to have the depth and like PsychIsColorModeIndexed
 *        instead of PsychColorModeType and depth, because the former is more orthagonal.
 */

//begin include once
#ifndef PSYCH_IS_INCLUDED_ScreenTypes
#define PSYCH_IS_INCLUDED_ScreenTypes

//constants
#define kPsychMaxPossiblePixelDepths    10
#define kPsychMaxPossibleDisplays       10
#define kPsychMaxPossibleCrtcs          10      // Maximum number of physical display scanout engines (crtc's) per screen.
#define kPsychMaxNumberWindowBuffers    2       // The maximum upper bound. Really it depends on display hardware and should be probed.

//includes
#include "Screen.h"

//project typedefs
typedef int PsychWindowIndexType;
typedef int PsychNumdexType;                                    //"numdex" is screen NUMber or window inDEX.
typedef int PsychDepthArrayType[kPsychMaxPossiblePixelDepths];

typedef struct {
    int                        numDepths;
    PsychDepthArrayType        depths;
} PsychDepthType;

typedef enum {
    kPsychFrontBuffer,
    kPsychBackBuffer,
} PsychWindowBufferType;

typedef enum {
    kPsychNoWindow =                    0,
    kPsychAnyWindow =                   1,            //  specify if we want to match to any window type.
    kPsychSingleBufferOnscreen =        2,
    kPsychDoubleBufferOnscreen =        3,
    kPsychVideoMemoryOffscreen =        4,
    kPsychSystemMemoryOffscreen =       5,
    kPsychTexture =                     6,
    kPsychProxyWindow =                 7
} PsychWindowType; //naming conventions dictate "PsychWindowTypeType" but thats silly.

typedef enum  {
    kPsychMonoscopic   = 0,             // Monoscopic display aka non-stereo. This is the default setting.
    kPsychOpenGLStereo = 1,             // Stereo display via built-in native OpenGL stereo facilities.
    kPsychCompressedTLBRStereo = 2,     // Vertical compression stereo with assingment top=left-eye, bottom=right-eye.
    kPsychCompressedTRBLStereo = 3,     // Vertical compression stereo with assingment top=right-eye, bottom=left-eye.
    kPsychFreeFusionStereo = 4,         // Stereo via free fusion of left and right view on a horizontally split screen.
    kPsychFreeCrossFusionStereo = 5,    // Stereo via free cross-fusion of left and right view on a horizontally split screen.
    kPsychAnaglyphRGStereo = 6,         // Stereo via color filter glasses (Anaglyph principle): Red-Green
    kPsychAnaglyphGRStereo = 7,         // Stereo via color filter glasses (Anaglyph principle): Green-Red
    kPsychAnaglyphRBStereo = 8,         // Stereo via color filter glasses (Anaglyph principle): Red-Blue
    kPsychAnaglyphBRStereo = 9,         // Stereo via color filter glasses (Anaglyph principle): Blue-Red
    kPsychDualWindowStereo = 10,        // Stereo for dual-display setups, where each view is rendered into its own onscreen window. (Imaging pipe only!)
    kPsychFrameSequentialStereo = 11,   // Stereo for frame-sequential display, but using PTB's own implementation instead of native OpenGL mode.
    kPsychDualStreamStereo = 12         // Stereo mode which needs two separate streams but is not used for frame-sequential stereo, native OpenGL or our own.
} PsychStereoDisplayType;

typedef enum {
    kPsychUnknownColor = 0,
    kPsychRGBColor,             // means PsychRGBColorType
    kPsychRGBAColor,            // means RGBColorType with alpha channel
    kPsychIndexColor            // means PsychIndexColorType
} PsychColorModeType;

typedef struct {
    double r;
    double g;
    double b;
    double a;
} PsychRGBAColorType;

typedef struct {
    double r;
    double g;
    double b;
} PsychRGBColorType;

typedef struct {
    double i;    //color index
} PsychIndexColorType;

typedef struct {
    PsychColorModeType mode;
    union {
        PsychRGBAColorType      rgba;
        PsychRGBColorType       rgb;
        PsychIndexColorType     index;
    } value;
} PsychColorType;

typedef struct {
    int                 screenNumber;
    PsychRectType       rect;
    PsychDepthType      depth;
    PsychColorModeType  mode;
    int                 nominalFrameRate;
} PsychScreenSettingsType;

typedef struct {
    int foo;
    //put more screen info here.
    //this is
} PsychScreenRecordType;

//indices into 3D or 4D matrix array returned by GetImage
#define kPsychRedPlane          0
#define kPsychGreenPlane        1
#define kPsychBluePlane         2
#define kPsychAlphaPlane        3
#define kPsychIndexPlane        0

// Constants that can be ored together for the "ConserveVRAM"
// Screen Preference setting:
// Do not use AUX buffers, disable associated features.
#define kPsychDisableAUXBuffers 1
// Use Client storage for RAM->DMA texturing instead of VRAM texturing.
#define kPsychDontCacheTextures 2
// Do not use wglChoosePixelformat - extension on M$-Windows:
#define kPsychOverrideWglChoosePixelformat 4
// Do not use separate GL contexts for PTB vs. MOGL on f%%%% M$-Windows and its broken drivers:
#define kPsychDisableContextIsolation 8
// Do not attach stencil buffer attachments to framebuffer objects in order to work around broken drivers:
#define kPsychDontAttachStencilToFBO 16
// Do not share heavyweight ressources like shaders, VBO's, FBO's, display lists, textures etc. between
// OpenGL contexts:
#define kPsychDontShareContextRessources 32
// Request explicit use of a software implemented renderer instead of the GPU. This may or may not
// be supported on a specific platform and is silently ignored if unsupported. On MacOS/X, this
// selects the Apple software floating point renderer -- Allows float rendering on non-capable hardware,
// e.g., for testing etc.
#define kPsychUseSoftwareRenderer 64

// Request that Windows PTB calls SetForegroundWindow() and SetFocus() on new onscreen windows.
// This makes them stay more reliably in foreground on especially f%%#$@ Windows systems, but
// interferes with the Jave GetChar implementation, so its only enabled optionally:
#define kPsychEnforceForegroundWindow 128

// Possible workaround for NVidia driver bug: Skip a few commands inside creation code for
// usermode OpenGL contexts -- The skipped commands are non-essential...
#define kPsychUseWindowsContextSharingWorkaround1 256

// Hint to PTB that any calls should be avoided that can cause a CPU-GPU synchronization, e.g.,
// error checking and handling calls. This is a hint: May be overriden at high debug levels etc.
#define kPsychAvoidCPUGPUSync 512

// Tell PTB to use the opposite texture format of what its auto-detection thinks is optimal:
#define kPsychTextureUploadFormatOverride 1024

// Tell PTB to avoid EXT_framebuffer_blit if possible, e.g., in SCREENCopyWindow():
#define kPsychAvoidFramebufferBlitIfPossible 2048

// On MS-Windows, always unconditionally enable the work-around for broken beamposition
// queries, regardless if the automatic test detects problems or not:
#define kPsychUseBeampositionQueryWorkaround 4096

// On OS/X, enforce use of Cocoa/NSOpenGL API for GL context setup insted of CGL API, even
// if the requested window is a fullscreen window:
#define kPsychUseAGLForFullscreenWindows 8192

// On OS/X enforce use of regular Quartz composited windows for fullscreen
// displays via Cocoa, instead of Cocoa or CGL fullscreen contexts. On Windows Vista
// and later it will also enforce use of the DWM. A future version may enforce
// use of Compiz on Linux as well:
#define kPsychUseAGLCompositorForFullscreenWindows 16384

// Wait for onset of VBL before emitting a PsychOSFlipWindowBuffers() request, by
// use of beamposition queries and busy-spin-waiting. This to work-around broken
// support of sync-to-vbl on some operating systems (e.g., Leopard + NVidia) under
// certain conditions.
#define kPsychBusyWaitForVBLBeforeBufferSwapRequest 32768

// Don't use native beamposition query mechanism but either our own
// homegrown solution, or no beampos query at all:
#define kPsychDontUseNativeBeamposQuery 65536

// Disable the Aero DWM desktop composition manager on Windows Vista and later:
#define kPsychDisableAeroDWM 131072

// Disable SetForegroundWindow() and SetFocus() calls on Windows for onscreen windows:
#define kPsychPreventForegroundWindow 262144

// Override default setting for OpenML swap scheduling: Disable OpenML scheduling.
#define kPsychDisableOpenMLScheduling 524288

// Apply our "bypass gamma LUT's in 10 bit framebuffer mode" hack also on FireGL/FirePro: 2^20
#define kPsychBypassLUTFor10BitFramebuffer 1048576

// Apply 10 bit framebuffer mode hack, regardless what automatic detection says: 2^21
#define kPsychEnforce10BitFramebufferHack 2097152

// Do not return OS provided nominal framerate, but zero instead, to override defective OS queries: 2^22
#define kPsychIgnoreNominalFramerate 4194304

// Do not use new-style asynchronous flips with dedicated OpenGL swapbuffer context per flipper thread:
#define kPsychUseOldStyleAsyncFlips (1 << 23)

// Do not auto-enable fast offscreen window support or full imaging pipeline for stereo modes
// by default on GPU's which support it:
#define kPsychDontAutoEnableImagingPipeline (1 << 24)

// Use old-style setup of the override_redirect flag for Linux X11 windows, as it was done until end of 2012?
#define kPsychOldStyleOverrideRedirect (1 << 25)

// Force use of native beamposition query mechanism instead
// of our own homegrown solution, or no beampos query at all:
#define kPsychForceUseNativeBeamposQuery (1 << 26)

// Force use of OpenML timestamp workaround for overriding kms-pageflip events:
// Not needed on any shipping OS release atm., but left here in case it is ever
// needed.
#define kPsychForceOpenMLTSWorkaround (1 << 27)

// Force setting of the kPsychAssumeVCColorsGood flag for graphics capabilities, ie.,
// assume the gpu can process/interpolate vertex colors at full 32bpc float precision.
#define kPsychAssumeGfxCapVCGood (1 << 28)

// Skip wait until scanout out-of-vblank before issuing swaprequest:
#define kPsychSkipOutOfVblankWait (1 << 29)

//function protoptypes

//Accessors for PsychDepthType
void            PsychInitDepthStruct(PsychDepthType *depth);
int             PsychGetNumDepthsFromStruct(PsychDepthType *depth);
int             PsychGetValueFromDepthStruct(int index, PsychDepthType *depth);
void            PsychAddValueToDepthStruct(int value, PsychDepthType *depth);
psych_bool      PsychIsMemberDepthStruct(PsychDepthType *depth, PsychDepthType *depthSet);
void            PsychCopyDepthStruct(PsychDepthType *toDepth, PsychDepthType *fromDepth);

// Useless function, always returns RGBA mode, the only one we support.
PsychColorModeType      PsychGetColorModeFromDepthStruct(PsychDepthType *depth);

// Accessors for PsychColorType:
// Convert some color into our PsychColorType: This only used internally to set default
// colors, usually white or black: Careful! It only accepts double values and screws with
// anything else, due to some macro preprocessor magic inside:
void PsychLoadColorStruct(PsychColorType *color, PsychColorModeType mode,  ...);

// Convert from whatever color representation "color" is to RGBA color spec, the only
// spec we use for internal color handling:
void PsychCoerceColorMode(PsychColorType *color);

//end include once
#endif
