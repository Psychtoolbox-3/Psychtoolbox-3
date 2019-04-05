/*
    PsychToolbox2/Source/Common/WindowBank.h

    AUTHORS:

        Allen.Ingling@nyu.edu               awi
        mario.kleiner.de@gmail.com          mk

    PLATFORMS:

        All. Only the PsychTargetSpecificWindowRecordType struct
        is different for the different operating systems.

    HISTORY:

        07/18/02  awi           Wrote it.
        04/22/05  mk            Added new flags to PsychWindowRecordType for stereo mode, skipped frames, IFI estimate and timestamping
                                - These functions are mostly used by Screen('Flip') and Screen('DrawingFinished').
        07/22/05  mk            Removed constants for max number of windows. We resize dynamically now (see WindowBank.c)
        10/11/05  mk            Support for special Quicktime movie textures in .targetSpecific part of PsychWindowRecord added.
        12/27/05  mk            Added new targetSpecific - struct for the Win32 port of PTB.

    DESCRIPTION:

        WindowBank contains functions for storing records of open  onscreen and offscreen windows.

        We don't use real window pointers when we return a window pointer to the scripting
        environment.  Instead, we return an index which the PsychToolbox maps to a window
        pointer. Since the mapping is always from the index passed from the scripting environment,
        we can implement this with array lookup.

        The window pointer or screen number is frequently passed as a single number, and by using
        our own indices we make things easier for ourselves when disambiguating the meaning of a
        windowPointerOrScreenNumber value, because the value itself is a direct and unambiguouis
        reference.

*/

//begin include once
#ifndef PSYCH_IS_INCLUDED_WindowBank
#define PSYCH_IS_INCLUDED_WindowBank

/*
    includes
*/

#include "Screen.h"

//constants
#define PSYCH_MAX_SCREENS               10        //the total possible number of screens
#define PSYCH_LAST_SCREEN                9        //the highest possible screen index
#define PSYCH_FIRST_SCREEN               0        //the lowest possible screen index
#define PSYCH_ALLOC_SCREEN_RECORDS      10        //length of static array allocated to hold screen pointers.
#define PSYCH_FIRST_WINDOW              10        //the lowest possible windox index


//constants
#define PSYCH_INVALID_WINDEX                 -1
#define PSYCH_INVALID_SCUMBER                -1

// Maximum number of different hook chains:
#define MAX_SCREEN_HOOKS 17

// Maximum number of slots in windowRecords fboTable:
#define MAX_FBOTABLE_SLOTS 2+2+3+4+2

// Type of hook function attached to a specific hook chain slot:
#define kPsychShaderFunc    0
#define kPsychCFunc         1
#define kPsychMFunc         2
#define kPsychBuiltinFunc   3

// Detected capabilities of the gfx-hardware, as interrogated by PsychDetectAndAssignGfxCapabilities()
// at onscreen window creation time and stored in windowRecord->gfxcaps as part of a bitfield:
#define kPsychGfxCapFBO                     1            // Hw supports OpenGL FBOs as rendertargets.
#define kPsychGfxCapFPTex16                 2            // Hw supports 16bpc float textures.
#define kPsychGfxCapFPTex32                 4            // Hw supports 32bpc float textures.
#define kPsychGfxCapFPFBO16                 8            // Hw supports 16bpc float FBOs.
#define kPsychGfxCapFPFBO32                 16           // Hw supports 32bpc float FBOs.
#define kPsychGfxCapFPFilter16              32           // Hw supports bilinear filtering of 16bpc float textures.
#define kPsychGfxCapFPFilter32              64           // Hw supports bilinear filtering of 32bpc float textures.
#define kPsychGfxCapFPBlend16               128          // Hw supports alpha blending on 16 bpc float FBOs.
#define kPsychGfxCapFPBlend32               256          // Hw supports alpha blending on 32 bpc float FBOs.
#define kPsychGfxCapFP32Shading             512          // Hw supports IEEE 32 bit float precision throughout shaders.
#define kPsychGfxCapVCGood                  1024         // Hw supports unclamped vertex colors of high (at least 16 bit effective) precision.
#define kPsychGfxCapFBOMultisample          2048         // Hw supports multisampled rendering into FBO's, aka EXT_framebuffer_multisample.
#define kPsychGfxCapFBOBlit                 4096         // Hw supports blitting between FBO's, aka EXT_framebuffer_blit.
#define kPsychGfxCapNeedsUnsignedByteRGBATextureUpload 8192        // Hw requires use of GL_UNSIGNED_BYTE instead of GL_UNSIGNED_INT_8_8_8_8_REV for optimal RGBA8 texture upload.
#define kPsychGfxCapSupportsOpenML          16384        // System supports OML_sync_control extension of OpenML for precisely timed bufferswaps and stimulus onset timestamping.
#define kPsychGfxCapSNTex16                 32768        // Hw supports 16 bit signed normalized integer textures.
#define kPsychGfxCapUYVYTexture             65536        // Hw supports UYVY encoded textures. Used for GStreamer video capture/playback engine optimizations.
#define kPsychGfxCapNativeStereo            (1 << 17)    // Hw supports native OpenGL quad-buffered stereo (frame-sequential etc.).
#define kPsychGfxCapNPOTTex                 (1 << 18)    // Hw supports non-power-of-twp GL_TEXTURE_2D textures.
#define kPsychGfxCapSupportsBufferAge       (1 << 19)    // Hw supports EXT_buffer_age extension.
#define kPsychGfxCapFBOScaledResolveBlit    (1 << 20)    // Hw supports simultaneous multisample resolve and rescaling in one framebuffer blit.
#define kPsychGfxCapSmoothPrimitives        (1 << 21)    // Hw supports GL_POINT_SMOOTH, LINE smooth etc.
#define kPsychGfxCapFloatToIntRound         (1 << 22)    // Hw, when converting float colors to fixed point colors, rounds to nearest instead of truncating.

// Definition of flags for imagingMode of Image processing pipeline.
// These are used internally, but need to be exposed to Matlab as well.
#define kPsychNeedFastBackingStore 1        // Any FBO's needed at all?
#define kPsychNeedImageProcessing  2        // Any image processing needed at all?
#define kPsychNeedOutputConversion 4        // Output conversion blit needed?
#define kPsychNeedSeparateStreams  8        // Quad-buffered stereo needed?
#define kPsychNeedStereoMergeOp    16       // Merging stereo -> single fb needed?
#define kPsychNeed32BPCFloat       32       // float32 FBO's needed?
#define kPsychNeed16BPCFloat       64       // float16 FBO's needed?
#define kPsychNeed16BPCFixed       128      // fixed point 16 bpc FBO's needed?
#define kPsychNeedDualPass         256      // At least support for dual-pass processing needed?
#define kPsychNeedMultiPass        512      // Support for multi-pass processing needed?
#define kPsychNeedFastOffscreenWindows 1024 // Only support for fast FBO-backed Offscreen windows, everything else off.
#define kPsychHalfWidthWindow          2048 // This flag is also used as 'specialflag' for onscreen windows. Ask for windows with half-width, e.g., for dualview stereo...
#define kPsychUse32BPCFloatAsap        4096 // This flag asks to use 32 bpc float FBOs for all stages of the pipeline after the initial drawBufferFBO. The initial drawBufferFBO
                                            // should also be a 32 bpc float FBO if the hardware reliably supports framebuffer blending in 32bpc. Otherwise it should be a 16 bpc
                                            // float FBO to allow for hardware accelerated framebuffer blending.
#define kPsychNeedDualWindowOutput    16384 // Apply imaging pipelines output formatting to two visual output streams and distribute result to two onscreen
                                            // windows. This allows e.g., to drive two physical displays from one input stimulus image, but applying different
                                            // post-processing to the two outputs. An example is display mirroring or special HDR displays.
#define kPsychNeedGPUPanelFitter      32768     // This flag is used for imagingMode to signal need/use of the GPU based imaging pipeline panel fitter.
#define kPsychNeedOtherStreamInput (1 << 17)    // This flag signals that the input image for the other view channel of a stereo pipeline should be bound to texture unit 1 as 2nd input.
#define kPsychNeedRetinaResolution (1 << 18)    // Do not auto-enable panelfitter on Retina display, ie., user-framebuffer shall be full native Retina resolution.
#define kPsychNeedClientRectNoFitter (1 << 19)  // Do not use the panelfitter, but still use provided clientRect to define window geometry as seen by client code.
#define kPsychNeedFinalizedFBOSinks   (1 << 20) // finalizedFBO's must be backed by real OpenGL FBO's as image sinks, not by the system backbuffer, as would be the default.
#define kPsychUseExternalSinkTextures (1 << 21) // finalizedFBO's have backing OpenGL FBO's which use externally created and injected OpenGL textures as color attachment (aka rendertargets).
#define kPsychSinkIsMSAACapable       (1 << 22) // External sinks for finalizedFBO color buffer textures can process GL_TEXTURE_2D_MULTISAMPLE textures as input, not only GL_TEXTURE_2D textures.
#define kPsychEnableSRGBRendering     (1 << 23) // Use OpenGL sRGB rendering and blending if the target color framebuffer is in sRGB format.

// 'specialflags' fields, partially shared with imagingmode flags:
#define kPsychUseTextureMatrixForRotation   1       // Setting for 'specialflags' field of windowRecords that describe textures. If set, drawtexture routine should implement
                                                    // rotated drawing of textures via texture matrix, not via modelview matrix. To be set as flag in 'DrawTexture(s)'
#define kPsychDontDoRotation                2       // Setting for 'specialflags' field of windowRecords that describe textures. If set, drawtexture routine should implement
                                                    // rotated drawing of textures via shader, not via matrices, ie., just pass rotation angle to shader. To be set as flag in 'DrawTexture(s)'
#define kPsychHalfHeightWindow              8192    // This flag is also used as 'specialflag' for onscreen windows. Ask for windows with half-height, e.g., for interleaved stereo...
#define kPsychNative10bpcFBActive           1024    // Setting for 'specialflags' field of windowRecords: Means that this windowRecord is attached to a native 10bpc system framebuffer
                                                    // and needs some special handling it init, shutdown and during operation.
// Value 32768 is defined in ScreenTypes.h as kPsychBusyWaitForVBLBeforeBufferSwapRequest and is also used as a 'specialflags' setting to define
// this behaviour on a per-window basis.
#define kPsychTwiceWidthWindow              65536   // This flag is also used as 'specialflag' for onscreen windows. Ask for windows with twice-width, e.g., for packed pixel modes.

#define kPsychIsFullscreenWindow            4       // 'specialflags' setting 4 means: This is a fullscreen window.
#define kPsychNeedOpenMLWorkaround1         8       // 'specialflags' setting 8 means: This needs the special workarounds for slightly broken OpenML sync control ext.
#define kPsychOpenMLDefective               16      // 'specialflags' setting 16 means: Don't use OpenML sync control extension for timestamping anywhere.

#define kPsychGUIWindow                     32      // 'specialflags' setting 32 means: This window should behave like a regular GUI window, e.g, allow moving it.
#define kPsychPlanarTexture                 64      // 'specialflags' setting 64: This texture uses a planar storage format instead of pixel-interleaved.
#define kPsychDontAutoGenMipMaps            128     // 'specialflags' setting 128: This texture shall not auto-generate its mip-map chain on demand.
#define kPsychIsX11Window                   256     // 'specialflags' setting 256: This window is living on a Linux X11/XServer backend.
#define kPsychIsGLXWindow                   512     // 'specialflags' setting 512: This window is living on a Linux X11/GLX backend.
#define kPsychIsEGLWindow                   4096    // 'specialflags' setting 4096: This window is living on a EGL backend (X11/Wayland/GBM/Android/...)
#define kPsychSurfacelessContexts           16384   // 'specialflags' setting 16384: This windows main context and userspace contexts must not attach to windowing system framebuffer surfaces.
#define kPsychDontDeleteOnClose             (1 << 17) // 'specialflags' setting 2^17: Do not close this texture/offscreen window on a call to Screen('Close'), only if explicitely closed by handle.
#define kPsychNeedPostSwapLockedFlush       (1 << 18) // 'specialflags' setting 2^18: Window needs display lock protected pixelwrite+flush on framebuffer immediately after bufferswap.
#define kPsychGUIWindowWMPositioned         (1 << 19) // 'specialflags' setting 2^19: GUI window is positioned onscreen by window manager, users position spec is ignored.
#define kPsychNeedOpenMLTSWorkaround        (1 << 20) // 'specialflags' setting 2^20: KMS pageflip completion events are faulty on a FOSS driver: Use glXGetSyncValuesOML workaround.
#define kPsychClockPrecisionOneTimeWarningDone (1 << 21) // 'specialflags' setting 2^21: Signals that the one-time warning wrt. imprecise visual onset timestamp was issued.
#define kPsychIsDRI3Window                  (1 << 22) // 'specialflags' setting 2^22: This X11 window uses DRI3/Present for visual stimulus presentation.
#define kPsychBufferAgeWarningDone          (1 << 23) // 'specialflags' setting 2^23: One time warning for non-double-buffering due to ext_buffer_age queries already done.
#define kPsychSafeForDRI3                   (1 << 24) // 'specialflags' setting 2^24: This window is considered safe for use with DRI3/Present, given X-Server and Mesa version in use.
#define kPsychTripleWidthWindow             (1 << 25) // This flag is also used as 'imagingmode' for onscreen windows. Ask for windows with triple-width, e.g., for packed pixel modes.
#define kPsychFbOverrideSizeActive          (1 << 26) // This window has a windowRecord->rect set via the 'OpenWindow' fbOverrideRect, and nothing should ever change the rect during runtime.
#define kPsychSkipVsyncForFlipOnce          (1 << 27) // 'specialflags': Perform next flip on this window without VSYNC.
#define kPsychSkipTimestampingForFlipOnce   (1 << 28) // 'specialflags': Perform next flip on this window without waiting for swap completion + timestamping + timing correctness checks.
#define kPsychSkipSwapForFlipOnce           (1 << 29) // 'specialflags': Perform next flip on this window without actually performing the OpenGL bufferswap, iow. don't present to the onscreen window.
#define kPsychSkipWaitForFlipOnce           (1 << 30) // 'specialflags': Perform next flip on this window without waiting until the 'when' target time for the flip.
#define kPsychUseFineGrainedOnset           (1ULL << 31) // 'specialflags': Schedule flips not on a fixed refresh interval, but use some scheduling with finer time granularity if possible.
#define kPsychNeedVBODouble12Workaround     (1ULL << 32) // 'specialflags': Gfx driver bug makes < 2 component vertex attribute buffers problematic if GL_DOUBLE is used for submission.

// The following numbers are allocated to imagingMode flag above: A (S) means, shared with specialFlags:
// 1,2,4,8,16,32,64,128,256,512,1024,S-2048,4096,S-8192,16384,32768,S-65536,2^17,2^18,2^19,2^20,2^21,2^22,2^23,S-2^25. --> Flags of 2^24 as well as 2^26 and higher are available...

// The following numbers are allocated to specialFlags flag above: A (S) means, shared with imagingMode:
// 1,2,4,8,16,32,64,128,256,512,1024,S-2048,4096,S-8192, 16384, 32768, S-65536,2^17,2^18,2^19,2^20,2^21,2^22,2^23,2^24,S-2^25,2^26,2^27,2^28,2^29,2^30,2^31,2^32. --> Flags of 2^33 and higher are available...

// Definition of a single hook function spec:
typedef struct PsychHookFunction*   PtrPsychHookFunction;
typedef struct PsychHookFunction {
    PtrPsychHookFunction    next;
    char*                   idString;
    int                     hookfunctype;
    char*                   pString1;
    void*                   cprocfunc;
    unsigned int            shaderid;
    unsigned int            luttexid1;
} PsychHookFunction;

// Definition of an OpenGL Framebuffer object (FBO) for internal use.
typedef struct PsychFBO {
    GLuint                  fboid;          // Handle to FBO.
    GLuint                  coltexid;       // Texture handle for color buffer texture (color attachment zero). A multisampled Renderbuffer handle, if multisample > 0.
    GLuint                  ztexid;         // Texture handle for z-Buffer texture, if any. Zero otherwise.
    GLuint                  stexid;         // Texture handle for stencil-Buffer texture, if any. Zero otherwise.
    GLenum                  format;         // fboInternalFormat of the colorbuffer0 attachment.
    int                     width;          // Width of FBO.
    int                     height;         // Height of FBO.
    int                     multisample;    // Multisampling level of FBO: 0 == No multisampling. > 0 means Multisampled.
    GLenum                  textarget;      // Type of texture target for texture coltexid (GL_TEXTURE_RECTANGLE_EXT or GL_TEXTURE_2D etc.)
} PsychFBO;

// Typedefs for WindowRecord in WindowBank.h

// This support structure for async flips is supported on all non-Windows platforms, aka all Unix platforms:
// It gets attached to the asyncFlipInfo* of a windowRecord whenever async flips are used.
typedef struct PsychFlipInfoStruct {
    unsigned char           asyncstate;         // Current execution state of flip.
    unsigned char           flipperState;
    // Parameters that specify the flip request: See prototype of PsychFlipWindowBuffers() for meaning of arguments:
    int                     opmode;
    int                     multiflip;
    int                     vbl_synclevel;
    int                     dont_clear;
    double                  flipwhen;
    // Return arguments with results of flip: See prototype of PsychFlipWindowBuffers() for meaning of arguments:
    int                     beamPosAtFlip;
    double                  miss_estimate;
    double                  time_at_flipend;
    double                  time_at_onset;
    double                  vbl_timestamp;

    psych_thread            flipperThread;      // Thread handle for background flipping thread.
    psych_mutex             performFlipLock;    // Primary lock.
    psych_condition         flipperGoGoGo;      // Signalling condition variable to trigger execution of a flip request by the flipper thread.
} PsychFlipInfoStruct;


#if PSYCH_SYSTEM == PSYCH_OSX
// Definition of OS-X core graphics and Core OpenGL handles:
typedef struct{
    CGLContextObj       contextObject;
    CGLPixelFormatObj   pixelFormatObject;
    CGLContextObj       glusercontextObject;    // OpenGL context for userspace rendering code, e.g., moglcore...
    CGLContextObj       glswapcontextObject;    // OpenGL context for performing doublebuffer swaps in PsychFlipWindowBuffers().
    void*               deviceContext;          // Pointer to an AGLContext object, or a NULL-pointer.
    // NSWindow* type stored in void* to avoid "Cocoa/Objective-C pollution" in this header file.
    void*               windowHandle;           // Handle for Cocoa window when using windowed mode. (NULL in non-windowed mode).

    // NSOpenGLContext* type stored in void* to avoid "Cocoa/Objective-C pollution" in this header file.
    void*               nsmasterContext;        // Cocoa OpenGL master context.
    void*               nsswapContext;          // Cocoa OpenGL async swap context.
    void*               nsuserContext;          // Cocoa OpenGL userspace rendering context.
} PsychTargetSpecificWindowRecordType;
#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS
// Definition of Win32 Window handles, device handles and OpenGL contexts
typedef struct{
    HGLRC                     contextObject;          // OpenGL rendering context.
    HDC                       deviceContext;          // Device context of the window.
    HWND                      windowHandle;           // The window handle.
    PIXELFORMATDESCRIPTOR     pixelFormatObject;      // The context's pixel format object.
    HGLRC                     glusercontextObject;    // OpenGL context for userspace rendering code, e.g., moglcore...
    HGLRC                     glswapcontextObject;    // OpenGL context for performing doublebuffer swaps in PsychFlipWindowBuffers().
} PsychTargetSpecificWindowRecordType;
#endif

#if PSYCH_SYSTEM == PSYCH_LINUX

#if defined(PTB_USE_WAYLAND_PRESENT) || defined(PTB_USE_WAYLAND)
#include <wayland-client.h>
#endif

#ifdef PTB_USE_WAFFLE
// Definition of Linux/Waffle specific information:

#ifdef PTB_USE_WAYLAND
// For the Linux Waffle + Wayland backend:
typedef struct {
    struct waffle_context*    contextObject;                  // Primary OpenGL rendering context.
    int                       pixelFormatObject;              // Just here for compatibility. Its a dummy entry without meaning.
    struct waffle_display*    deviceContext;                  // Pointer to the Waffle display connection.
    void*                     privDpy;                        // Pointer to the private X11 display connection for non-OpenGL ops.
    struct waffle_window*     windowHandle;                   // Handle to the Waffle onscreen window handle.
    struct wl_surface*        xwindowHandle;                  // Associated Wayland "window", if any.
    struct waffle_context*    glusercontextObject;            // OpenGL context for userspace rendering code, e.g., moglcore...
    struct waffle_context*    glswapcontextObject;            // OpenGL context for performing doublebuffer swaps in PsychFlipWindowBuffers().
    struct wl_list            presentation_feedback_list;     // Used for Wayland backend presentation_feedback extension to queue feedback events.
} PsychTargetSpecificWindowRecordType;
#else
// For the Linux Waffle generic backend:
typedef struct {
    struct waffle_context*    contextObject;                  // GLX OpenGL rendering context.
    int                       pixelFormatObject;              // Just here for compatibility. Its a dummy entry without meaning.
    struct waffle_display*    deviceContext;                  // Pointer to the X11 display connection.
    Display*                  privDpy;                        // Pointer to the private X11 display connection for non-OpenGL ops.
    struct waffle_window*     windowHandle;                   // Handle to the onscreen window.
    Window                    xwindowHandle;                  // Associated X-Window if any.
    struct waffle_context*    glusercontextObject;            // OpenGL context for userspace rendering code, e.g., moglcore...
    struct waffle_context*    glswapcontextObject;            // OpenGL context for performing doublebuffer swaps in PsychFlipWindowBuffers().
#ifdef PTB_USE_WAYLAND_PRESENT
    struct wl_list            presentation_feedback_list;     // Used for Wayland backend presentation_feedback extension to queue feedback events.
#endif
} PsychTargetSpecificWindowRecordType;
#endif

#else
// Definition of Linux/X11 specific information:
typedef struct{
    GLXContext        contextObject;       // GLX OpenGL rendering context.
    int               pixelFormatObject;   // Just here for compatibility. Its a dummy entry without meaning.
    Display*          deviceContext;       // Pointer to the X11 display connection.
    Display*          privDpy;             // Pointer to the private X11 display connection for non-OpenGL ops.
    GLXWindow         windowHandle;        // Handle to the onscreen window.
    Window            xwindowHandle;       // Associated X-Window if any.
    GLXContext        glusercontextObject; // OpenGL context for userspace rendering code, e.g., moglcore...
    GLXContext        glswapcontextObject; // OpenGL context for performing doublebuffer swaps in PsychFlipWindowBuffers().
} PsychTargetSpecificWindowRecordType;
#endif

// End of Linux targetSpecific struct.
#endif

#define kPsychUnaffiliatedWindow    -1  // valid value for screenNumber field of a window record meaning that that pixel format
                                        // and alignment of the window are not set to match those of any display surface.

typedef struct _PsychWindowRecordType_ *PsychWindowRecordPntrType;

//typedefs for the window bank.  We use the same structure for both windows and textures.
typedef struct _PsychWindowRecordType_{

    //need to be divided up according to use for textures, windows, or both.
    PsychWindowType             windowType;
    int                         winsysType;         // Windowing/Display system backend type: 0 = Classic (=X11/GLX on Linux), or a Waffle backend type if waffle is used.
    int                         glApiType;          // Type of OpenGL rendering API in use: 0 = Classic desktop OpenGL-1/2/3/4, 10 = GL-ES1.0, 20 = GL-ES2.0, 30 = GL-ES3.0 ...
    int                         screenNumber;       // kPsychUnaffiliated is -1 and means the offscreen window is unaffiliated.
    PsychWindowIndexType        windowIndex;
    void                        *surface;
    size_t                      surfaceSizeBytes;   // Estimate of used system memory in bytes. Only used for accounting and debugging output.
    PsychRectType               rect;               // Bounding rectangle of true window framebuffer -- Normalized to always have top-left corner in (0,0)!
    PsychRectType               globalrect;         // Same as rect, but not normalized -- Describes real window location in global desktop reference frame.
    PsychRectType               clientrect;         // Bounding rectangle of *useable* framebuffer size for client stimulus drawing code. This should be returned
                                                    // to calling client code as 'Rect', 'WindowSize' etc. and used by Screen 2D drawing functions as reference
                                                    // window size. Normally clientrect == rect, but for special display modes, e.g., stereo modes 2-5, Bits+
                                                    // Color++ mode etc., the drawable area for stimuli may be of different size than the physical window backbuffer,
                                                    // e.g., half or twice as wide or high than the true scanout buffer. This is usually due to special pixel
                                                    // encodings for video scanout (two display pixels packed into one framebuffer pixel, one display pixel
                                                    // encoded by two adjacent framebuffer pixels (color++), special ultra-wide or stretched/compressed stereo displays...
    psych_bool                  isValid;            // between when we allocate the record and when we fill in values.
    int                         depth;              // Number of bits per pixel -- Attention: This is often misleading and not really reliable!
    int                         nrchannels;
    int                         bpc;                // Number of bits per color component: Typically 8 for LDR, but can be 16 or 32 for float drawables.
    int                         multiSample;        // Number of samples to use for anti-aliasing of each drawn pixel:
                                                    // 0=Standard GL, n>0 : Use Multisampling, try to get at least n samples/pixel.
    //used only when this structure holds a window
    psych_bool                  actualEnableBlending;
    GLenum                      actualSourceAlphaBlendingFactor;
    GLenum                      actualDestinationAlphaBlendingFactor;
    GLenum                      nextSourceAlphaBlendingFactor;
    GLenum                      nextDestinationAlphaBlendingFactor;

    //Used only when this structure holds a texture:
    GLuint                      *textureMemory;
    size_t                      textureMemorySizeBytes;
    GLuint                      textureNumber;
    int                         textureOrientation;     // Orientation of texture data in internal storage. Defines texcoord assingment.
    int                         texturecache_slot;      // Reference of cache structure for this texture, if any...
    GLenum                      texturetarget;          // Explicit target type of texture (GL_TEXTURE_2D, ...)
    // The following three are only used for injecting special textures into PTB, e.g., High Dynamic range textures in floating point format.
    // They default to zero, which means: Derive texture representation from depth.
    GLint                       textureinternalformat;  // Explicit definition of glinternalformat for texture creation.
    GLenum                      textureexternalformat;  // Explicit definitin of external format for texture creation.
    GLenum                      textureexternaltype;    // Explicit definition of data type for texture creation.
    GLint                       textureFilterShader;    // Optional GLSL program handle for a shader to apply during PsychBlitTextureToDisplay().
    GLint                       textureLookupShader;    // Optional GLSL handle for nearest neighbour texture drawing shader.
    GLint                       textureByteAligned;     // 0 = No knowledge about byte alignment of texture data. > 1, texture rows are x byte aligned.
    GLint                       texturePlanarShader[4]; // Optional GLSL program handles for shaders to apply to planar storage textures - 4 handles for 4 possible channel counts.
    GLint                       textureI420PlanarShader; // Optional GLSL program handle for shader to convert a YUV-I420 planar texture into a standard RGBA8 texture.
    GLint                       textureI800PlanarShader; // Optional GLSL program handle for shader to convert a Y8-I800 planar texture into a standard RGBA8 texture.
    GLint                       multiSampleFetchShader;  // Optional GLSL program handler for shader to fetch from multisample texture.

    psych_bool                  needsViewportSetup;     // Set on userspace OpenGL contexts of onscreen windows to signal need for glViewport setup and other one-time
                                                        // stuff on first Screen('BeginOpenGL'). Also (ab)used for textures and offscreen windows to track "dirty" state.

    //line stipple attributes, for windows not textures.
    psych_bool                  stippleEnabled;
    GLushort                    stipplePattern;
    GLint                       stippleFactor;

    GLboolean                   colorMask[4];           // psych_bool 4 element array which encodes the glColorMask() for this window.
    unsigned int                gfxcaps;                // Bitfield of gfx-cards capabilities and limitations: See constants kPsychGfxCapXXXX above.
    psych_uint64                specialflags;           // Container for all kind of special flags...
    char                        gpuCoreId[8];           // ASCII string with name of GPU core, as detected by PsychDetectAndAssignGfxCapabilities().
    int                         maxTextureSize;         // Maximum size of a texture supported by GPU.
    int                         stereomode;             // MK: Is this a stereo window? 0=non-stereo, >0 == specific type of stero.
    int                         stereodrawbuffer;       // MK: Which drawbuffer is active in stereo? 0=left, 1=right, 2=none
    PsychWindowRecordPntrType   slaveWindow;            // MK: In stereomode 10 (dual-window stereo) Either NULL or windowrecord of right view window.
    PsychWindowRecordPntrType   parentWindow;           // MK: Ptr. to windowRecord of the parent window, or NULL if this window doesn't have a parent.
    int                         targetFlipFieldType;    // MK: Usually == -1 (=Don't care). Can select that bufferswap should always happen in even frames (=0) or odd frames (=1). Useful for frame sequential stereo.
    psych_bool                  auxbuffer_dirty[2];     // MK: State of auxbuffers 0 and 1: Dirty or not? (For stereo algs.)
    int                         nrIFISamples;           // MK: nrIFISamples and IFIRunningSum are used to calculate an
    double                      IFIRunningSum;          // MK: accurate estimate of the real interframe interval (IFI) in Flip.
    double                      time_at_last_vbl;       // MK: Timestamp (system-time) at last VBL detected by Flip. This is the same as a returned vbltimestamp of stimulus onset from 'Flip'
    double                      VideoRefreshInterval;   // MK: Estimated video refresh interval of display. Can be different to IFI.
    double                      ifi_beamestimate;       // MK: Yet another video refresh estimate, based on beamposition method (or 0 if invalid).
    int                         VBL_Endline;            // MK: Estimated scanline which marks end of VBL area.
    int                         VBL_Startline;          // MK: Start scanline of VBL area.
    psych_bool                  PipelineFlushDone;      // MK: Will be set by SCREENDrawingFinished to signal pipeline flush.
    psych_bool                  backBufferBackupDone;   // MK: Will be set by SCREENDrawingFinished to signal backbuffer backup.
    psych_bool                  vSynced;                // MK: Flag that stores VSYNC enable state: TRUE = Sync to VBL, FALSE = Don't.
    int                         hybridGraphics;         // "true" on a hybrid graphics laptop that needs special treatment: 1=DRI3/Present renderoffload, 2=Optimus NVidia prop.
    int                         nr_missed_deadlines;    // MK: Counter, incremented by Flip if it detects a missed/skipped frame.
    int                         flipCount;              // Counter of total number of finished flip operations - A swapcounter.
    double                      rawtime_at_swapcompletion; // Raw timestamp of swapcompletion (result without high-precision timestamping).
    double                      time_at_swaprequest;    // Timestamp taken immediately before call to PsychOSFlipWindowBuffers(); - Before swaprequest submission.
    double                      time_post_swaprequest;  // Timestamp taken immediately after call PsychOSFlipWindowBuffers();
    double                      postflip_vbltimestamp;  // Optional timestamp taken after flip completion via PsychGetVBLTimeAndCount();
    double                      osbuiltin_swaptime;     // Optional timestamp of swap completion computed via PsychOSGetSwapCompletionTimestamp();
    double                      gpuRenderTime;          // GPU time spent on rendering. Only returned if a query object is successfully generated.
    GLuint                      gpuRenderTimeQuery;     // Handle to the GPU time query object. 0 if none assigned.
    psych_int64                 reference_ust;          // UST reference timestamp of vblank with count reference_msc from OpenML. (Optional)
    psych_int64                 reference_msc;          // MSC reference vblank count from OpenML. (Optional)
    psych_int64                 reference_sbc;          // SBC reference swapbuffers count from OpenML. (Optional)
    psych_int64                 target_sbc;             // Target SBC value for next glXWaitForSbcOML() call from OpenML. (Optional)
    psych_int64                 lastSwaptarget_msc;     // Target MSC value for which most recent swap was scheduled by DRM/DRI2 from OpenML. (Optional)
    psych_int64                 submitted_sbc;          // Running count of submitted swaprequests - e.g., used by Wayland. (Optional)
    int                         swapevents_enabled;     // State of swap events: 0 = Disabled/Unsupported, 1 = Enabled for use by usercode, 2 = Enabled for use by us.
    int                         swapcompletiontype;     // Type of last completed swap: 0 = Unknown, 1 = Pageflip, 2 = Exchange, 3 = Copy.

    // Pointers to temporary arrays with gamma tables to upload to the gfx-card at next Screen('Flip'):
    // They default to NULL and get possibly set in Screen('LoadNormalizedGammaTable'):
    float*                      inRedTable;
    float*                      inGreenTable;
    float*                      inBlueTable;
    int                         inTableSize;            // Number of slots in the LUT tables.
    int                         loadGammaTableOnNextFlip;   // Type of upload operation: 0 = None, 1 = Load on next Flip via OS gamma table routines, then reset flag.

    // Settings for the image processing and hook callback pipeline: See PsychImagingPipelineSupport.hc for definition and implementation:
    int                         applyColorRangeToDoubleInputMakeTexture;    // Should colorRange also affect uint8 textures created from double input in Screen('MakeTexture')?
    double                      colorRange;                                 // Maximum allowable color component value. See SCREENColorRange.c for explanation.
    GLuint                      unclampedDrawShader;                        // Handle of GLSL shader object for drawing of non-texture stims without vertex color clamping. Zero by default.
    GLuint                      defaultDrawShader;                          // Default GLSL shader object for drawing of non-texture stims. Zero by default.
    GLuint                      smoothPointShader;                          // GLSL shader to implement point smoothing via point sprites.
    double                      currentColor[4];                            // Current unclamped but colorrange remapped RGBA drawcolor for whatever drawop, as spec'd by PsychSetGLColor().
    double                      clearColor[4];                              // Window clear color (as GL double vector) to use in PsychGLClear();
    int                         imagingMode;                                // Master mode switch for imaging and callback hook pipeline.
    PtrPsychHookFunction        HookChain[MAX_SCREEN_HOOKS];                // Array of pointers to the hook-chains for different hooks.
    psych_bool                  HookChainEnabled[MAX_SCREEN_HOOKS];         // Array of Booleans to en-/disable single chains temporarily.

    // Indices into our FBO table: The special value -1 means: Don't use.
    int                         drawBufferFBO[2];                   // Storage for drawing FBOs: These are the targets of all drawing operations before
                                                                    // Screen('DrawingFinished') or Screen('Flip') is called. They are read-only wrt.
                                                                    // to all following processing stages: 0=Left eye (or mono) channel, 1=Right eye channel.
    int                         inputBufferFBO[2];                  // Storage for input FBOs: These are the input to all further pipeline stages. They are
                                                                    // initialized with the content of the drawBufferFBO's, either directly, by setting both
                                                                    // equal, so that drawBufferFBO->inputBufferFBO is a no-op, or - if some special purpose
                                                                    // preprocessing is needed - the drawBufferFBO's are copied to these. The main application
                                                                    // is multisample anti-aliasing: drawBufferFBO's must be multisampled, but rest of pipeline
                                                                    // must not be multisampled, so we perform the multisample-resolve operation during copy
                                                                    // from drawBuffers to inputBuffers.
    int                         processedDrawBufferFBO[3];          // These contain the final output of all per-view channel postprocessing operations:
                                                                    // 0=Left eye (or mono) channel, 1=Right eye channel, 2=Temporary bounce buffer for iterative
                                                                    // multi-pass processing. These provide the input for the stereo merger in stereo modes that
                                                                    // require merging of the two views, e.g., anaglyph stereo.
    int                         preConversionFBO[4];                // preConversion FBO's: FBO zero/one are the targets for any stereo merge operations. FBO two and three are
                                                                    // (optionally) temporary bounce buffers for multipass post processing.
    int                         finalizedFBO[2];                    // This is the final framebuffer: Usually the system backbuffer, but could be something special.

    PsychFBO*                   fboTable[MAX_FBOTABLE_SLOTS];       // This array contains pointers to the FBO structs which are referenced by the indices above.
    int                         fboCount;                           // This contains the number of FBO's in fboTable.

    // Cached handles for display lists -- used for recycling in compute intense drawing functions:
    GLuint                      fillOvalDisplayList;
    GLuint                      frameOvalDisplayList;

    // Pointer to double-array of auxiliary parameters for bound shaders - or NULL by default.
    double*                     auxShaderParams;
    int                         auxShaderParamsCount;

    //Used only when this structure holds a window:
    //platform specific stuff goes within the targetSpecific structure.  Defined in PsychVideoGlue and accessors are in PsychWindowGlue.c
    //Only use abstracted accessors on this structure, otherwise you will break platform portability.
    PsychTargetSpecificWindowRecordType     targetSpecific;  //defined within

    PsychFlipInfoStruct*        flipInfo;                           // This is either a NULL-Ptr, or it points to structure with all the info
                                                                    // needed for implementing asynchronous flip operations. Its always NULL on
                                                                    // MS-Windows, non-NULL on Linux/OSX as soon as async flips are used at least once.
                                                                    // See SCREENFlip.c and flipping routines in PsychWindowSupport.c for more details...

    psych_uint64                gpu_preflip_Surfaces[2];            // Framebuffer addresses of the primary-/secondary surfaces on GPU before flip.

    // Support for framelock / swaplock / output lock / genlock via swap groups / swap barriers extensions:
    GLuint                      swapGroup;                          // Swap group handle of swap group for this window, zero if none assigned.
    GLuint                      swapBarrier;                        // Swap barrier handle of swap barrier for this window, zero if none assigned.

    GLint                       panelFitterParams[11];              // Parameters used for panel fitting.

    double                      text2DMatrix[2][3];                 // 2D affine transform matrix for text.

    double*                     proj;                               // Override projection matrices (4x4x[1or2]) or NULL for default projection.

    // Used only when this structure holds a window:
    // CAUTION FIXME TODO: Due to some pretty ugly circular include dependencies in the #include chain of
    // PTB, this field can not be used in files that #define PSYCH_DONT_INCLUDE_TEXTATTRIBUTES_IN_WINDOWRECORD,
    // e.g., PsychGraphicsHardwareHALSupport.c/h. Reason: The order of #include's in Screen.h prevents the
    // definition of this field of being available when PsychGraphicsHardwareHALSupport.h is included,
    // and we can't move that file down the include chain due to other nasty dependency issues!
    // So in such files, which need definition of PsychWindowRecordType, but are not related to
    // any text handling, we redef this to a void*. Important: This *must* be the last definition
    // of the struct PsychWindowRecordType, otherwise memory-alignment of the datastructure would
    // change in a corruptive way!!! So use the PSYCH_DONT_INCLUDE_TEXTATTRIBUTES_IN_WINDOWRECORD with
    // great caution! A proper fix would be to redesign the include dependency chain - either the whole
    // design, or at least the parts related to PsychTextAttributes, but that's unlikely to be a pain-free
    // procedure - it will cause severe breakage and the need for lot's of fixups and auditing, so defer
    // this to a later point in time :-(
    #ifndef PSYCH_DONT_INCLUDE_TEXTATTRIBUTES_IN_WINDOWRECORD
    PsychTextAttributes         textAttributes;
    #else
    void*                       textAttributes;
    #endif
} PsychWindowRecordType;


//various window and screen related
void                    PsychFindScreenWindowFromScreenNumber(int screenNumber, PsychWindowRecordType **winRec);
PsychError              InitWindowBank(void);
PsychError              CloseWindowBank(void);
void                    PsychCreateWindowRecord(PsychWindowRecordType **winRec);
PsychError              FreeWindowRecordFromIndex(PsychWindowIndexType windowIndex);
PsychError              FreeWindowRecordFromPntr(PsychWindowRecordType *winRec);
int                     PsychCountOpenWindows(PsychWindowType winType);
psych_bool              IsValidScreenNumber(PsychNumdexType numdex);
psych_bool              IsValidScreenNumberOrUnaffiliated(PsychNumdexType numdex);
psych_bool              IsWindowIndex(PsychNumdexType numdex);
psych_bool              IsScreenNumberUnaffiliated(PsychNumdexType numdex);
PsychError              FindWindowRecord(PsychWindowIndexType windowIndex, PsychWindowRecordType **windowRecord);
PsychError              FindScreenRecord(int screenNumber, PsychScreenRecordType **screenRecord);
psych_bool              PsychIsLastOnscreenWindow(PsychWindowRecordType *windowRecord);
psych_bool              PsychIsLastOnscreenWindowOnScreen(PsychWindowRecordType *windowRecord);
void                    PsychCreateVolatileWindowRecordPointerList(int *numWindows, PsychWindowRecordType ***pointerList);
void                    PsychDestroyVolatileWindowRecordPointerList(PsychWindowRecordType **pointerList);
void                    PsychAssignParentWindow(PsychWindowRecordType *childWin, PsychWindowRecordType *parentWin);
PsychWindowRecordType*  PsychGetParentWindow(PsychWindowRecordType *windowRecord);

//end include once
#endif
