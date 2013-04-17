/*
  PsychToolbox3/Source/Linux/Screen/PsychWindowGlueWaffle.c
    
  PLATFORMS:    
    
  This is the Linux/Waffle version only.

  AUTHORS:

  Mario Kleiner      mk    mario.kleiner at tuebingen.mpg.de

  HISTORY:
    
  2/12/13            mk    Created - Derived from Linux/X11-GLX version.

  DESCRIPTION:
    
  Functions in this file comprise an abstraction layer for probing and controlling window state, except for window content.  
        
  Each C function which implements a particular Screen subcommand should be platform neutral.  For example, the source to SCREENPixelSizes() 
  should be platform-neutral, despite that the calls in OS X and Linux to detect available pixel sizes are different.  The platform 
  specificity is abstracted out in C files which end it "Glue", for example PsychScreenGlue, PsychWindowGlue, PsychWindowTextClue.

  NOTES:
    
  Preformatted via: indent -linux -l240 -i4 PsychWindowGlueWaffle.c

  TO DO: 
     
*/

#ifdef PTB_USE_WAFFLE

#include "Screen.h"

#include <waffle.h>
#include <waffle_glx.h>
#include <waffle_x11_egl.h>
#include <waffle_gbm.h>
#include <waffle_wayland.h>

// First time init?
static psych_bool firstTime = TRUE;

// Use X11 / GLX backend?
static psych_bool useX11 = FALSE;
static psych_bool useGLX = FALSE;

// Use GLX version 1.3 setup code? Enabled INTEL_SWAP_EVENTS and other goodies...
static psych_bool useGLX13 = FALSE;

// Event base for GLX extension:
static int glx_error_base, glx_event_base;

// Number of currently open onscreen windows:
static int open_windowcount = 0;

// Shared waffle display connection handle for whole session:
static struct waffle_display *wdpy = NULL;

// Also share native underlying EGL display for EGL based backends (NULL otherwise):
static EGLDisplay egl_display = NULL;

// Tracking of currently bound OpenGL rendering context for master-thread:
static struct waffle_context *currentContext = NULL;

// Forward define prototype for glewContextInit(), which is normally not a public function:
GLenum glewContextInit(void);

/*
  PsychOSOpenOnscreenWindow()
    
  Creates the pixel format and the context objects and then instantiates the context onto the screen.
    
  -The pixel format and the context are stored in the target specific field of the window recored.  Close
  should clean up by destroying both the pixel format and the context.
    
  -We mantain the context because it must be be made the current context by drawing functions to draw into 
  the specified window.
    
  -We maintain the pixel format object because there seems to be now way to retrieve that from the context.
    
  -To tell the caller to clean up PsychOSOpenOnscreenWindow returns FALSE if we fail to open the window. It 
  would be better to just issue an PsychErrorExit() and have that clean up everything allocated outside of
  PsychOpenOnscreenWindow().
*/
psych_bool PsychOSOpenOnscreenWindow(PsychScreenSettingsType * screenSettings, PsychWindowRecordType * windowRecord, int numBuffers, int stereomode, int conserveVRAM)
{
    PsychRectType screenrect;
    CGDirectDisplayID dpy;
    int scrnum;
    XSetWindowAttributes attr;
    unsigned long mask;
    Window win = 0;
    GLXWindow glxwindow;
    XVisualInfo *visinfo = NULL;
    int i, x, y, width, height, nrconfigs, buffdepth;
    GLenum glerr;
    int32_t attrib[41];
    int attribcount = 0;
    int stereoenableattrib = 0;
    int depth, bpc;
    int windowLevel;
    int major, minor;
    int xfixes_event_base1, xfixes_event_base2;
    psych_bool xfixes_available = FALSE;
    psych_bool newstyle_setup = FALSE;
    int32_t opengl_api;
    char backendname[16];
    char backendname2[16];
    int32_t oldBackend;
    struct waffle_config *config;
    struct waffle_window *window;
    struct waffle_context *ctx;

    // Define default rendering backend:
    #ifdef PTB_USE_GLES1
    opengl_api = WAFFLE_CONTEXT_OPENGL_ES1;
    #else
    opengl_api = WAFFLE_CONTEXT_OPENGL;
    #endif

    // Define waffle window system backend to use by default: It is our good'ol X11/GLX,
    // unless EGL is explicitely requested at compile time:
    static int32_t init_attrs[3] = {
        WAFFLE_PLATFORM,
        #ifdef PTB_USE_EGL
        WAFFLE_PLATFORM_X11_EGL,
        #else
        WAFFLE_PLATFORM_GLX,
        #endif
        0,
    };

    // Map the logical screen number to the corresponding X11 display connection handle
    // for the corresponding X-Server connection.
    PsychGetCGDisplayIDFromScreenNumber(&dpy, screenSettings->screenNumber);

    // TODO FIXME: We currently don't have any way of selecting the target X-Screen 'scrnum' for
    // our window, as Waffle does not yet support selection of X-Screen. It always opens
    // windows on the display's default screen. Therefore this is mostly a dead placeholder for now:
    scrnum = PsychGetXScreenIdForScreen(screenSettings->screenNumber);

    // First time invocation?
    if (firstTime) {
        // Initialize waffle for selected display system backend:
        if (PsychPrefStateGet_Verbosity() > 2) {
            printf("PTB-INFO: Using FOSS Waffle display backend library, written by Chad Versace, Copyright 2012 Intel.\n");
        }

        // Override default windowing system backend selection with requested type, if any requested:
        if (getenv("PSYCH_USE_DISPLAY_BACKEND")) {
            if (!strcmp(getenv("PSYCH_USE_DISPLAY_BACKEND"), "glx")) windowRecord->winsysType = (int) WAFFLE_PLATFORM_GLX; 
            if (!strcmp(getenv("PSYCH_USE_DISPLAY_BACKEND"), "x11egl")) windowRecord->winsysType = (int) WAFFLE_PLATFORM_X11_EGL; 
            if (!strcmp(getenv("PSYCH_USE_DISPLAY_BACKEND"), "wayland")) windowRecord->winsysType = (int) WAFFLE_PLATFORM_WAYLAND; 
            if (!strcmp(getenv("PSYCH_USE_DISPLAY_BACKEND"), "gbm")) windowRecord->winsysType = (int) WAFFLE_PLATFORM_GBM; 
            if (!strcmp(getenv("PSYCH_USE_DISPLAY_BACKEND"), "android")) windowRecord->winsysType = (int) WAFFLE_PLATFORM_ANDROID; 
        }

        // If any backend chosen from calling code or by env-variable, assign it as new requested choice,
        // but backup current setting for possible later restore:
        oldBackend = init_attrs[1];
        if (windowRecord->winsysType > 0) init_attrs[1] = (int32_t) windowRecord->winsysType;

        if (!waffle_init(init_attrs) && (waffle_error_get_code() != WAFFLE_ERROR_ALREADY_INITIALIZED)) {
            // Failed to init with requested backend. Try different other backends:
            if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not initialize requested Waffle display backend, trying fallbacks: %s.\n", waffle_error_to_string(waffle_error_get_code()));

            // X11/GLX (classic):
            init_attrs[1] = WAFFLE_PLATFORM_GLX;
            if (!waffle_init(init_attrs)) {
                // X11/EGL:
                init_attrs[1] = WAFFLE_PLATFORM_X11_EGL;
                if (!waffle_init(init_attrs)) {
                    // Wayland:
                    init_attrs[1] = WAFFLE_PLATFORM_WAYLAND;
                    if (!waffle_init(init_attrs)) {
                        // GBM/DRM/KMS raw framebuffer:
                        init_attrs[1] = WAFFLE_PLATFORM_GBM;
                        if (!waffle_init(init_attrs)) {
                            // Android:
                            init_attrs[1] = WAFFLE_PLATFORM_ANDROID;
                            if (!waffle_init(init_attrs)) {
                                // Final fail:
                                if (PsychPrefStateGet_Verbosity() > 0) {
                                    printf("PTB-ERROR: Could not initialize any Waffle display backend - Error: %s.\n", waffle_error_to_string(waffle_error_get_code()));
                                    printf("PTB-ERROR: Try to fix the reason for the error, then restart Octave/Matlab, then retry.\n");
                                }
                                return(FALSE);
                            }
                        }
                    }
                }
            }
        }
        else if (waffle_error_get_code() == WAFFLE_ERROR_ALREADY_INITIALIZED) {
            // waffle_init() skipped/no-opped because waffle was already initialized with a
            // backend and this choice is permanent until the host runtime (Matlab/Octave) is
            // restarted. Was a different backend chosen? If so we warn the user that this is
            // not possible without a restart of the runtime environment:
            if ((oldBackend != init_attrs[1]) && (PsychPrefStateGet_Verbosity() > 1)) {
                printf("PTB-WARNING: Tried to choose a different display backend, but this is only possible on first invocation in a session.\n");
                printf("PTB-WARNING: Sticking to previous initial choice. You must restart Octave/Matlab before you can choose a different backend.\n");
            }

            // Restore actual backend choice from backup:
            init_attrs[1] = oldBackend;
        }
        else {
            // Success. init_attrs[1] contains correct choice, set oldBackend to it for completeness:
            oldBackend = init_attrs[1];
        }

        // First-Time init done:
        firstTime = FALSE;
    }

    // Set windowing system backend type to truly selected type:
    windowRecord->winsysType = (int) init_attrs[1];

    // Translate spec to human readable name and spec string:
    switch (init_attrs[1]) {
    case WAFFLE_PLATFORM_GLX:
        sprintf(backendname, "X11/GLX");
        sprintf(backendname2, "glx");
        break;

    case WAFFLE_PLATFORM_X11_EGL:
        sprintf(backendname, "X11/EGL");
        sprintf(backendname2, "x11egl");
        break;

    case WAFFLE_PLATFORM_WAYLAND:
        sprintf(backendname, "Wayland/EGL");
        sprintf(backendname2, "wayland");
        break;

    case WAFFLE_PLATFORM_GBM:
        sprintf(backendname, "GBM/EGL");
        sprintf(backendname2, "gbm");
        break;

    case WAFFLE_PLATFORM_ANDROID:
        sprintf(backendname, "Android/EGL");
        sprintf(backendname2, "android");
        break;
    }

    // Announce actual choice of backend to runtime environment. This is a marker
    // to, e.g., moglcore, so it can adapt its context/gl setup:
    setenv("PSYCH_USE_DISPLAY_BACKEND", backendname2, 1);

    if (PsychPrefStateGet_Verbosity() > 2) {
        printf("PTB-INFO: Waffle display backend '%s' initialized [%s].\n", backendname, backendname2);
    }

    // Waffle display connection not yet initialized? Create one if not yet done. Currently
    // we are limited to one connection per session, until we move waffle init and connection
    // init to the screen glue and manage those as we do with classic X11 display connections.
    // Display connection gets closed on last onscreen window close, (re)opened on first window open:
    if (NULL == wdpy) {
        // Connect it to the chosen or default display device:
        if ((windowRecord->winsysType == WAFFLE_PLATFORM_GLX) || (windowRecord->winsysType == WAFFLE_PLATFORM_X11_EGL)) {
            // X11 backend: Use X11 display as specified by our Psychtoolbox screenId:
            wdpy = waffle_display_connect(DisplayString(dpy));
            if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Trying to connect Waffle to display '%s'.\n", DisplayString(dpy));
        }
        else {
            // Other backend: Environment variable or NULL aka auto-selected default, if variable is undefined:
            wdpy = waffle_display_connect(getenv("PSYCH_WAFFLE_DISPLAY"));
            if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Trying to connect Waffle to display '%s'.\n", getenv("PSYCH_WAFFLE_DISPLAY"));
        }
    }

    if (!wdpy) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Could not connect Waffle to display: %s.\n",
                                                      waffle_error_to_string(waffle_error_get_code()));
        return(FALSE);
    }

    // Extract EGL_Display for backends which use EGL:
    union waffle_native_display* wafflenatdis = waffle_display_get_native(wdpy);
    if (wafflenatdis) {
        switch (windowRecord->winsysType) {
            case WAFFLE_PLATFORM_GLX:
                egl_display = NULL;
                break;

            case WAFFLE_PLATFORM_X11_EGL:
                egl_display = wafflenatdis->x11_egl->egl_display;
                break;

            case WAFFLE_PLATFORM_WAYLAND:
                egl_display = wafflenatdis->wayland->egl_display;
                break;

            case WAFFLE_PLATFORM_GBM:
                egl_display = wafflenatdis->gbm->egl_display;
                break;

            case WAFFLE_PLATFORM_ANDROID:
                egl_display = NULL;
                break;

            default:
                egl_display = NULL;
        }

        free(wafflenatdis);
        wafflenatdis = NULL;
    }
    else {
        egl_display = NULL;
    }

    // Allow usercode to override our default backend choice via environment variable:
    if (getenv("PSYCH_USE_GFX_BACKEND")) {
        if (!strcmp(getenv("PSYCH_USE_GFX_BACKEND"), "gl")) opengl_api = WAFFLE_CONTEXT_OPENGL;
        if (!strcmp(getenv("PSYCH_USE_GFX_BACKEND"), "gles1")) opengl_api = WAFFLE_CONTEXT_OPENGL_ES1;
        if (!strcmp(getenv("PSYCH_USE_GFX_BACKEND"), "gles2")) opengl_api = WAFFLE_CONTEXT_OPENGL_ES2;
        if (!strcmp(getenv("PSYCH_USE_GFX_BACKEND"), "gles3")) opengl_api = WAFFLE_CONTEXT_OPENGL_ES3;
    }

    switch (opengl_api) {
        case WAFFLE_CONTEXT_OPENGL:
            sprintf(backendname, "OpenGL");
            break;

        case WAFFLE_CONTEXT_OPENGL_ES1:
            sprintf(backendname, "OpenGL-ES 1");
            break;

        case WAFFLE_CONTEXT_OPENGL_ES2:
            sprintf(backendname, "OpenGL-ES 2");
            break;

        case WAFFLE_CONTEXT_OPENGL_ES3:
            sprintf(backendname, "OpenGL-ES 3");
            break;

        default:
            sprintf(backendname, "UNKNOWN");
    }

    // Try the requested backend, try alternate backends on failure:
    if (!waffle_display_supports_context_api(wdpy, opengl_api)) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Selected Waffle display backend does not support requested OpenGL rendering API '%s': %s. Trying fallbacks...\n",
                                                      backendname, waffle_error_to_string(waffle_error_get_code()));
        // Try fallbacks: OpenGL > OpenGL-ES1 > OpenGL-ES2 > OpenGL-ES3
        if (waffle_display_supports_context_api(wdpy, WAFFLE_CONTEXT_OPENGL)) {
            opengl_api = WAFFLE_CONTEXT_OPENGL;
        }
        else if (waffle_display_supports_context_api(wdpy, WAFFLE_CONTEXT_OPENGL_ES1)) {
            opengl_api = WAFFLE_CONTEXT_OPENGL_ES1;
        }
        else if (waffle_display_supports_context_api(wdpy, WAFFLE_CONTEXT_OPENGL_ES2)) {
            opengl_api = WAFFLE_CONTEXT_OPENGL_ES2;
        }
        else if (waffle_display_supports_context_api(wdpy, WAFFLE_CONTEXT_OPENGL_ES3)) {
            opengl_api = WAFFLE_CONTEXT_OPENGL_ES3;
        }
        else {
            // Game over:
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Waffle display backend does not support any OpenGL rendering API: %s. Game over!\n",
                                                          waffle_error_to_string(waffle_error_get_code()));
            return(FALSE);
        }
    }

    // Assign human readable name to selected OpenGL rendering api. Also assign
    // type flag for higher-level code:
    switch (opengl_api) {
        case WAFFLE_CONTEXT_OPENGL:
            sprintf(backendname, "OpenGL");
            sprintf(backendname2, "gl");
            windowRecord->glApiType = 0;
            break;

        case WAFFLE_CONTEXT_OPENGL_ES1:
            sprintf(backendname, "OpenGL-ES 1");
            sprintf(backendname2, "gles1");
            windowRecord->glApiType = 10;
            break;

        case WAFFLE_CONTEXT_OPENGL_ES2:
            sprintf(backendname, "OpenGL-ES 2");
            sprintf(backendname2, "gles2");
            windowRecord->glApiType = 20;
            break;

        case WAFFLE_CONTEXT_OPENGL_ES3:
            sprintf(backendname, "OpenGL-ES 3");
            sprintf(backendname2, "gles3");
            windowRecord->glApiType = 30;
            break;

        default:
            sprintf(backendname, "UNKNOWN");
            sprintf(backendname2, "X");
            windowRecord->glApiType = 0;
    }

    // Tell environment about our final rendering backend choice:
    setenv("PSYCH_USE_GFX_BACKEND", backendname2, 1);

    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Waffle display backend connected to OpenGL rendering API '%s' [%s].\n", backendname, backendname2);

    // OpenGL embedded subset api in use?
    if (windowRecord->glApiType > 0) {
        // Yes: Try to disable all OpenGL error checking for now during initial OpenGL-ES support bringup.
        // We know there will be many errors due to incompatibilities but just try to keep going and rely
        // on external trace and debug tools to find and fix errors at this point:
        PsychPrefStateSet_ConserveVRAM(PsychPrefStateGet_ConserveVRAM() | kPsychAvoidCPUGPUSync);
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Will try to disable/suppress all internal OpenGL error reporting/handling for OpenGL-ES operation.\n");
    }

    // Enable X11 specific setup code if running on X11 backend:
    useX11 = ((windowRecord->winsysType == WAFFLE_PLATFORM_GLX) || (windowRecord->winsysType == WAFFLE_PLATFORM_X11_EGL)) ? TRUE : FALSE;
    if (useX11) windowRecord->specialflags |= kPsychIsX11Window;

    // Ditto for GLX:
    useGLX = (windowRecord->winsysType == WAFFLE_PLATFORM_GLX) ? TRUE : FALSE;

    // Mark window as either GLX or EGL controlled:
    windowRecord->specialflags |= ((useGLX) ? kPsychIsGLXWindow : kPsychIsEGLWindow);

    // Retrieve windowLevel, an indicator of where non-fullscreen windows should
    // be located wrt. to other windows. 0 = Behind everything else, occluded by
    // everything else. 1 - 999 = At layer windowLevel -> Occludes stuff on layers "below" it.
    // 1000 - 1999 = At highest level, but partially translucent / alpha channel allows to make
    // regions transparent. 2000 or higher: Above everything, fully opaque, occludes everything.
    // 2000 is the default.
    windowLevel = PsychPrefStateGet_WindowShieldingLevel();

    // Init userspace GL context to safe default:
    windowRecord->targetSpecific.glusercontextObject = NULL;
    windowRecord->targetSpecific.glswapcontextObject = NULL;

    // Default to use of one shared x-display connection "dpy" for all onscreen windows
    // on a given x-display and x-screen:
    windowRecord->targetSpecific.privDpy = dpy;

    // XFixes extension version 2.0 or later available and initialized?
    if (useX11 && XFixesQueryExtension(dpy, &xfixes_event_base1, &xfixes_event_base2) && XFixesQueryVersion(dpy, &major, &minor) && (major >= 2)) xfixes_available = TRUE;
    major = minor = 0;

    if (useX11 && useGLX) {
        // Init GLX extension, get its version, determine if at least V1.3 supported:
        useGLX13 = (glXQueryExtension(dpy, &glx_error_base, &glx_event_base) && glXQueryVersion(dpy, &major, &minor) && ((major > 1) || ((major == 1) && (minor >= 3))));

        // Initialze GLX-1.3 protocol support. Use if possible:
        glXChooseFBConfig = (PFNGLXCHOOSEFBCONFIGPROC) glXGetProcAddressARB("glXChooseFBConfig");
        glXGetFBConfigAttrib = (PFNGLXGETFBCONFIGATTRIBPROC) glXGetProcAddressARB("glXGetFBConfigAttrib");
        glXGetVisualFromFBConfig = (PFNGLXGETVISUALFROMFBCONFIGPROC) glXGetProcAddressARB("glXGetVisualFromFBConfig");
        glXCreateWindow = (PFNGLXCREATEWINDOWPROC) glXGetProcAddressARB("glXCreateWindow");
        glXCreateNewContext = (PFNGLXCREATENEWCONTEXTPROC) glXGetProcAddressARB("glXCreateNewContext");
        glXDestroyWindow = (PFNGLXDESTROYWINDOWPROC) glXGetProcAddressARB("glXDestroyWindow");
        glXSelectEvent = (PFNGLXSELECTEVENTPROC) glXGetProcAddressARB("glXSelectEvent");
        glXGetSelectedEvent = (PFNGLXGETSELECTEDEVENTPROC) glXGetProcAddressARB("glXGetSelectedEvent");

        // Check if everything we need from GLX-1.3 is supported:
        if (!useGLX13 || !glXChooseFBConfig || !glXGetVisualFromFBConfig || !glXCreateWindow || !glXCreateNewContext || !glXDestroyWindow || !glXSelectEvent || !glXGetSelectedEvent || !glXGetFBConfigAttrib) {
            useGLX13 = FALSE;
            if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Not using GLX-1.3 extension. Unsupported? Some features may be disabled.\n");
        } else {
            useGLX13 = TRUE;
        }
    }
    else {
        // GLX backend not used, therefore also GLX-1.3 disabled:
        useGLX13 = FALSE;
    }

    // Check if this should be a fullscreen window, and if not, what its dimensions
    // should be:
    PsychGetScreenRect(screenSettings->screenNumber, screenrect);
    if (PsychMatchRect(screenrect, windowRecord->rect)) {
        // This is supposed to be a fullscreen window with the dimensions of
        // the current display/desktop:
        x = 0;
        y = 0;
        width = PsychGetWidthFromRect(screenrect);
        height = PsychGetHeightFromRect(screenrect);

        // Mark this window as fullscreen window:
        windowRecord->specialflags |= kPsychIsFullscreenWindow;

        // Copy absolute screen location and area of window to 'globalrect',
        // so functions like Screen('GlobalRect') can still query the real
        // bounding gox of a window onscreen:
        PsychGetGlobalScreenRect(screenSettings->screenNumber, windowRecord->globalrect);
    } else {
        // Window size different from current screen size:
        // A regular desktop window with borders and control icons is requested, e.g., for debugging:
        // Extract settings:
        x = windowRecord->rect[kPsychLeft];
        y = windowRecord->rect[kPsychTop];
        width = PsychGetWidthFromRect(windowRecord->rect);
        height = PsychGetHeightFromRect(windowRecord->rect);

        // Copy absolute screen location and area of window to 'globalrect',
        // so functions like Screen('GlobalRect') can still query the real
        // bounding gox of a window onscreen:
        PsychCopyRect(windowRecord->globalrect, windowRecord->rect);
    }

    // Select OpenGL context API for use with this window:
    attrib[attribcount++] = WAFFLE_CONTEXT_API;
    attrib[attribcount++] = opengl_api;

    // Which display depth is supported?
    depth = PsychGetValueFromDepthStruct(0, &(screenSettings->depth));

    // Select requested depth per color component 'bpc' for each channel:
    bpc = 8; // We default to 8 bpc == RGBA8
    if (windowRecord->depth == 30) {
        bpc = 10;
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Trying to enable at least 10 bpc fixed point framebuffer.\n");
    }
    if (windowRecord->depth == 64) {
        bpc = 16;
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Trying to enable 16 bpc fixed point framebuffer.\n");
    }
    if (windowRecord->depth == 128) {
        bpc = 32;
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Trying to enable 32 bpc fixed point framebuffer.\n");
    }

    attrib[attribcount++] = WAFFLE_RED_SIZE;   // Setup requested minimum depth of each color channel:
    attrib[attribcount++] = (depth > 16) ? bpc : 1;
    attrib[attribcount++] = WAFFLE_GREEN_SIZE;
    attrib[attribcount++] = (depth > 16) ? bpc : 1;
    attrib[attribcount++] = WAFFLE_BLUE_SIZE;
    attrib[attribcount++] = (depth > 16) ? bpc : 1;
    attrib[attribcount++] = WAFFLE_ALPHA_SIZE;
    // Alpha channel needs special treatment:
    if (bpc != 10) {
        // Non 10 bpc drawable: Request a 'bpc' alpha channel if the underlying framebuffer
        // is in true-color mode ( >= 24 bpp format). If framebuffer is in 16 bpp mode, we
        // don't have/request an alpha channel at all:
        attrib[attribcount++] = (depth > 16) ? bpc : 0; // In 16 bit mode, we don't request an alpha-channel.
    } else {
        // 10 bpc drawable: We have a 32 bpp pixel format with R10G10B10 10 bpc per color channel.
        // There are at most 2 bits left for the alpha channel, so we request an alpha channel with
        // minimum size 1 bit --> Will likely translate into a 2 bit alpha channel:
        attrib[attribcount++] = 1;
    }

    // Stereo display support: If stereo display output is requested with OpenGL native stereo,
    // we request a stereo-enabled rendering context. Or at least we would, if Waffle would support this.
    if (stereomode == kPsychOpenGLStereo) {
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: OpenGL native quad-buffered stereomode requested, but this is not supported by the Waffle backend.\n");
    }

    // Multisampling support:
    if (windowRecord->multiSample > 0) {
        // Request a multisample buffer:
        attrib[attribcount++] = WAFFLE_SAMPLE_BUFFERS;
        attrib[attribcount++] = 1;

        // Request at least multiSample samples per pixel:
        attrib[attribcount++] = WAFFLE_SAMPLES;
        attrib[attribcount++] = windowRecord->multiSample;
    }

    // Support for OpenGL 3D rendering requested?
    if (PsychPrefStateGet_3DGfx()) {
        // Yes. Allocate and attach a 24 bit depth buffer and 8 bit stencil buffer:
        attrib[attribcount++] = WAFFLE_DEPTH_SIZE;
        attrib[attribcount++] = 24;
        attrib[attribcount++] = WAFFLE_STENCIL_SIZE;
        attrib[attribcount++] = 8;

        // Alloc an accumulation buffer as well?
        if (PsychPrefStateGet_3DGfx() & 2) {
            // Yes: Allocate an accumulation buffer if possible:
            attrib[attribcount++] = WAFFLE_ACCUM_BUFFER;
            attrib[attribcount++] = 1;
        }
    }

    // Double buffering requested?
    attrib[attribcount++] = WAFFLE_DOUBLE_BUFFERED;
    attrib[attribcount++] = (numBuffers >= 2) ? 1 : 0;

    // Finalize attrib array:
    attrib[attribcount++] = 0;

    // Choose waffle configuration for attrib's - the equivalent of
    // a pixelformat or framebuffer config in GLX speak:
    config = waffle_config_choose(wdpy, attrib);
   
    if (!config) {
        // Failed to find matching visual: Could it be related to request for unsupported native 10 bpc framebuffer?
        if ((windowRecord->depth == 30) && (bpc == 10)) {
            // 10 bpc framebuffer requested: Let's see if we can get a visual by lowering our demand to 8 bpc:
            for (i = 0; i < attribcount && attrib[i] != WAFFLE_RED_SIZE; i++);
            attrib[i + 1] = 8;
            for (i = 0; i < attribcount && attrib[i] != WAFFLE_GREEN_SIZE; i++);
            attrib[i + 1] = 8;
            for (i = 0; i < attribcount && attrib[i] != WAFFLE_BLUE_SIZE; i++);
            attrib[i + 1] = 8;
            for (i = 0; i < attribcount && attrib[i] != WAFFLE_ALPHA_SIZE; i++);
            attrib[i + 1] = 1;

            // Retry:
            config = waffle_config_choose(wdpy, attrib);
        }
    }

    if (!config) {
        // Failed to find matching visual: Could it be related to multisampling?
        if (windowRecord->multiSample > 0) {
            // Multisampling requested: Let's see if we can get a visual by
            // lowering our demand:
            for (i = 0; i < attribcount && attrib[i] != WAFFLE_SAMPLES; i++);
            while (!config && windowRecord->multiSample > 0) {
                attrib[i + 1]--;
                windowRecord->multiSample--;
                config = waffle_config_choose(wdpy, attrib);
            }

            // Either we have a valid visual at this point or we still fail despite
            // requesting zero samples.
            if (!config) {
                // We still fail. Disable multisampling by requesting zero multisample buffers:
                for (i = 0; i < attribcount && attrib[i] != WAFFLE_SAMPLE_BUFFERS; i++);
                windowRecord->multiSample = 0;
                attrib[i + 1] = 0;
                config = waffle_config_choose(wdpy, attrib);
            }
        }

        // Worked? If not, see if we can lower our 3D settings if in 3D mode:
        if (!config && PsychPrefStateGet_3DGfx()) {
            // Ok, retry with a 16 bit depth buffer...
            for (i = 0; i < attribcount && attrib[i] != WAFFLE_DEPTH_SIZE; i++);
            if (attrib[i] == WAFFLE_DEPTH_SIZE && i < attribcount) attrib[i + 1] = 16;
            if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Have to use 16 bit depth buffer instead of 24 bit buffer due to limitations of your gfx-hardware or driver. Accuracy of 3D-Gfx may be limited...\n");
            config = waffle_config_choose(wdpy, attrib);

            if (!config) {
                // Failed again. Retry with disabled stencil buffer:
                if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Have to disable stencil buffer due to limitations of your gfx-hardware or driver. Some 3D Gfx algorithms may fail...\n");
                for (i = 0; i < attribcount && attrib[i] != WAFFLE_STENCIL_SIZE; i++);
                if (attrib[i] == WAFFLE_STENCIL_SIZE && i < attribcount) attrib[i + 1] = 0;
                config = waffle_config_choose(wdpy, attrib);
            }
        }
    }

    if (!config) {
        // Ok, nothing worked. Require the absolute minimum configuration:
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not get any display configuration matching your minimum requirements. Trying to get the bare minimum config of a RGB framebuffer of arbitrarily low resolution now.\n");
        attribcount = 0;
        attrib[attribcount++] = WAFFLE_CONTEXT_API;
        attrib[attribcount++] = opengl_api;
        attrib[attribcount++] = WAFFLE_RED_SIZE;
        attrib[attribcount++] = 1;
        attrib[attribcount++] = WAFFLE_GREEN_SIZE;
        attrib[attribcount++] = 1;
        attrib[attribcount++] = WAFFLE_BLUE_SIZE;
        attrib[attribcount++] = 1;
        attrib[attribcount++] = WAFFLE_ALPHA_SIZE;
        attrib[attribcount++] = WAFFLE_DONT_CARE;
        config = waffle_config_choose(wdpy, attrib);
    }

    // Finally?
    if (!config) {
        // Game over:
        if (PsychPrefStateGet_Verbosity() > 0) printf("\nPTB-ERROR[waffle_config_choose() failed]: Couldn't get any suitable visual from display backend.\n\n");
        return (FALSE);
    }

    // Old style of override_redirect handling requested? This was used until beginning 2013
    // and worked well for us, but it prevents the windowmanager from seeing properties on
    // our windows which allow us to control desktop composition, e.g., on KDE/KWIN and GNOME-3/Mutter,
    // as well as on other wm's compliant with latest ICCCM spec:
    // Ok, for now we only use the new-style path if we are running under KDE/KWin and user
    // doesn't explicitely override/forbid that choice. Otherwise we use the old path, as
    // that seems to perform better, at least on tested Unity/compiz, GNOME3-Shell and LXDE/OpenBox.
    if ((PsychPrefStateGet_ConserveVRAM() & kPsychOldStyleOverrideRedirect) || !getenv("KDE_FULL_SESSION")) {
        // Old style: Always override_redirect to lock out window manager, except when a real "GUI-Window"
        // is requested, which needs to behave and be treated like any other desktop app window:
        attr.override_redirect = (windowRecord->specialflags & kPsychGUIWindow) ? 0 : 1;
    } else {
        // New style: override_redirect by default:
        newstyle_setup = TRUE;

        attr.override_redirect = 1;

        // Don't override if it is a GUI window, for some reasons as in classic path:
        if (windowRecord->specialflags & kPsychGUIWindow) attr.override_redirect = 0;

        // Don't override if it is a fullscreen window. The NETM_FULLSCREEN state should
        // take care of fullscreen windows nicely without need to override. Although we
        // could override for transparent fullscreen windows if we were extra paranoid,
        // we *must not* override for opaque (standard) fullscreen windows, because that
        // would prevent us from disabling desktop composition on our window by use of
        // special properties -- the window manager wouldn't notice those properties/requests
        // on a override_redirected window.
        // Now here's a catch: While the WM picks up our window properties this way, and at
        // least KWIN honors them by disabling composition, the WM also takes this opportunity
        // to "misplace" our onscreen window on multi-display setups. It follows various rules on
        // where a fullscreen window should be placed on a multi-display setup, but while these
        // rules are pretty sensible for regular desktop users, they are absolutely not what we
        // want. How to solve? We let the WM have its way during window creation, so it picks up
        // our window props. Later on (see below) after the Window is mapped and the WM satisified,
        // we set the override_redirect flag to lock out the WM, then move the window to its proper
        // location, no that the WM can't interfere anymore.
        if (windowRecord->specialflags & kPsychIsFullscreenWindow) attr.override_redirect = 0;
    }
    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Using %s-style override-redirect setup path for onscreen window creation.\n", (newstyle_setup) ? "new" : "old");

    // Create our onscreen window:
    window = waffle_window_create(config, width, height);

    if (useX11) XSync(dpy, False);

    // Setup X11 window properties:
    if (useX11) {
        // Retrieve underlying native X11 window:
        union waffle_native_window *wafflewin = waffle_window_get_native(window);
        if (windowRecord->winsysType == WAFFLE_PLATFORM_GLX) {
            struct waffle_glx_window* glx = wafflewin->glx;
            win = (Window) glx->xlib_window;
        }
        else {
            win = (Window) ((struct waffle_x11_egl_window*) (wafflewin->x11_egl))->xlib_window;
        }

        free(wafflewin);
        wafflewin = NULL;

        // Set hints and properties:
        {
            XSizeHints sizehints;
            sizehints.x = x;
            sizehints.y = y;
            sizehints.width = width;
            sizehints.height = height;
            sizehints.flags = USSize | USPosition;
            XSetNormalHints(dpy, win, &sizehints);
            XSetStandardProperties(dpy, win, "PTB Onscreen window", "PTB Onscreen window", None, (char **)NULL, 0, &sizehints);
        }

        // Setup window transparency for user input (keyboard and mouse events):
        if (xfixes_available && (windowLevel < 1500)) {
            // Define region as an empty input region:
            XserverRegion region = XFixesCreateRegion(dpy, NULL, 0);

            // Assign as region in which window receives input events, thereby
            // setting the input region to empty, so the window is transparent
            // to any input events like key presses or mouse clicks:
            XFixesSetWindowShapeRegion(dpy, win, ShapeInput, 0, 0, region);

            // Destroy region after assignment:
            XFixesDestroyRegion(dpy, region);
        }
    }

    if (useX11) XSync(dpy, False);

    // Make sure a potential slaveWindow of us resides on the same X-Screen == has same screenNumber as us,
    // otherwise trying to perform OpenGL context resource sharing would end badly:
    if ((windowRecord->slaveWindow) && (windowRecord->slaveWindow->screenNumber != screenSettings->screenNumber)) {
        // Ohoh! Let's abort with some more helpful error message than a simple hard application crash:
        if (PsychPrefStateGet_Verbosity() > 0) printf("\nPTB-ERROR:[waffle_context_create() resource sharing] Our peer window resides on a different screen, which is forbidden. Aborting.\n\n");
        return (FALSE);
    }

    // Create associated OpenGL rendering context: We use ressource
    // sharing of textures, display lists, FBO's and shaders if 'slaveWindow'
    // is assigned for that purpose as master-window.
    ctx = waffle_context_create(config, ((windowRecord->slaveWindow) ? windowRecord->slaveWindow->targetSpecific.contextObject : NULL));
    if (!ctx) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("\nPTB-ERROR:[waffle_context_create() failed] OpenGL context creation failed: %s\n\n", waffle_error_to_string(waffle_error_get_code()));
        return (FALSE);
    }

    // Store the handles:

    // windowHandle is a waffle_window*:
    windowRecord->targetSpecific.windowHandle = window;

    // xwindowHandle stores the underlying X-Window:
    windowRecord->targetSpecific.xwindowHandle = win;

    // waffle_display*
    windowRecord->targetSpecific.deviceContext = wdpy;

    // waffle_context*
    windowRecord->targetSpecific.contextObject = ctx;

    // Create rendering context for async flips with identical config as main context, share all heavyweight ressources with it:
    windowRecord->targetSpecific.glswapcontextObject = waffle_context_create(config, windowRecord->targetSpecific.contextObject);
    if (windowRecord->targetSpecific.glswapcontextObject == NULL) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("\nPTB-ERROR[SwapContextCreation failed]: Creating a private OpenGL context for async flips failed for unknown reasons.\n\n");
        return (FALSE);
    }

    // External 3D graphics support enabled?
    if (PsychPrefStateGet_3DGfx()) {
        // Yes. We need to create an extra OpenGL rendering context for the external
        // OpenGL code to provide optimal state-isolation. The context shares all
        // heavyweight ressources likes textures, FBOs, VBOs, PBOs, display lists and
        // starts off as an identical copy of PTB's context as of here.

        // Create rendering context with identical config as main context, share all heavyweight ressources with it:
        windowRecord->targetSpecific.glusercontextObject = waffle_context_create(config, windowRecord->targetSpecific.contextObject);
        if (windowRecord->targetSpecific.glusercontextObject == NULL) {
            if (PsychPrefStateGet_Verbosity() > 0) printf("\nPTB-ERROR[UserContextCreation failed]: Creating a private OpenGL context for MOGL failed for unknown reasons.\n\n");
            return (FALSE);
        }
    }

    // Release config info:
    waffle_config_destroy(config);

    if (useX11) XSync(dpy, False);

    // Setup window transparency:
    if ((windowLevel >= 1000) && (windowLevel < 2000)) {
        // For windowLevels between 1000 and 1999, make the window background transparent, so standard GUI
        // would be visible, wherever nothing is drawn, i.e., where alpha channel is zero:

        // Levels 1000 - 1499 and 1500 to 1999 map to a master opacity level of 0.0 - 1.0:        
        unsigned int opacity = (unsigned int)(0xffffffff * (((float)(windowLevel % 500)) / 499.0));

        if (useX11) {
            // Get handle on opacity property of X11:
            Atom atom_window_opacity = XInternAtom(dpy, "_NET_WM_WINDOW_OPACITY", False);

            // Assign new value for property:
            XChangeProperty(dpy, win, atom_window_opacity, XA_CARDINAL, 32, PropModeReplace, (unsigned char *)&opacity, 1);
        }
    }

    // Is this a non-transparent (fully opaque), non-GUI, fullscreen onscreen window?
    if (!((windowLevel >= 1000) && (windowLevel < 2000)) && !(windowRecord->specialflags & kPsychGUIWindow) && (windowRecord->specialflags & kPsychIsFullscreenWindow)) {
        // Yes. This is a standard stimulus presentation window which should get best
        // timing precision and performance for stimulus presentation. We don't want
        // any desktop composition to interfere with it, so it is eligible for direct
        // page-flipping (unredirected). If we are running under KDE's KWin desktop
        // manager, then we can explicitely ask KWin to disable compositing while
        // our onscreen window is open, by setting a special NETWM property on the window.
        // This approach has just become a NETWM standard that should work with other
        // compositors in the future, e.g., Mutter/GNOME-3 as of 18th December 2012.
        //
        // On other compositors, e.g., compiz / unity et al. this problem is solved by
        // asking them to unredirect_fullscreen_windows, as done by PsychGPUControl.m during
        // installation of PTB.
        //
        // Btw. for other properties that KDE supports/understands see function create_netwm_atoms()
        // in file netwm.cpp, e.g., at http://code.woboq.org/kde/kdelibs/kdeui/windowmanagement/netwm.cpp.html
        //
        if (useX11) {
            // Set KDE-4 specific property: This is supported since around KWin 4.6, since July 2011:
            unsigned int dontcomposite = 1;
            Atom atom_window_dontcomposite = XInternAtom(dpy, "_KDE_NET_WM_BLOCK_COMPOSITING", False);

            // Assign new value for property:
            XChangeProperty(dpy, win, atom_window_dontcomposite, XA_CARDINAL, 32, PropModeReplace, (unsigned char *)&dontcomposite, 1);

            // Set the standardized NETWM property. This is supported in Mutter (== GNOME-3) since
            // 18. December 2012 (see last comment/patch in https://bugzilla.gnome.org/show_bug.cgi?id=683020 ),
            // and will supposedly get supported by other compositing window managers in the future as well,
            // e.g., future KWin/KDE releases or possibly Unity/Compiz:
            dontcomposite = 1;
            Atom atom_window_dontcomposite2 = XInternAtom(dpy, "_NET_WM_BYPASS_COMPOSITOR", False);

            // Assign new value for property:
            XChangeProperty(dpy, win, atom_window_dontcomposite2, XA_CARDINAL, 32, PropModeReplace, (unsigned char *)&dontcomposite, 1);
        }
    }

    // Is this a non-GUI fullscreen window? If so, set the fullscreen NETWM property:
    if (useX11 && !(windowRecord->specialflags & kPsychGUIWindow) && (windowRecord->specialflags & kPsychIsFullscreenWindow)) {
        // Yes. Set the fullscreen state hint. Any well behaved window manager should understand this as
        // a request to turn the window into a completely decorationless fullscreen window, something
        // similar to what you'd get with override_redirect = 1, just in a less dirty way, that still
        // allows us to communicate our wishes, e.g., wrt. desktop composition, to the window manager.
        Atom stateFullscreen = XInternAtom(dpy, "_NET_WM_STATE_FULLSCREEN", False);
        XChangeProperty(dpy, win, XInternAtom(dpy, "_NET_WM_STATE", False), XA_ATOM, 32, PropModeReplace, (unsigned char *)&stateFullscreen, 1);
    }

    if (useX11) {
        // For some reason we need to use unsigned long and long here instead of
        // int32_t etc., despite the fact that on a 64-Bit build, a long is 64-Bit
        // and on a 32-Bit build, a long is 32-Bit, whereas the XChangeProperty()
        // request says a single unit is 32-Bits? Anyway, it works correctly on a
        // 64-Bit build, so this seems to be magically ok.
        struct MwmHints {
            unsigned long flags;
            unsigned long functions;
            unsigned long decorations;
            long input_mode;
            unsigned long status;
        };

        enum {
            MWM_HINTS_FUNCTIONS = (1L << 0),
            MWM_HINTS_DECORATIONS = (1L << 1),

            MWM_FUNC_ALL = (1L << 0),
            MWM_FUNC_RESIZE = (1L << 1),
            MWM_FUNC_MOVE = (1L << 2),
            MWM_FUNC_MINIMIZE = (1L << 3),
            MWM_FUNC_MAXIMIZE = (1L << 4),
            MWM_FUNC_CLOSE = (1L << 5)
        };

        Atom mwmHintsProperty = XInternAtom(dpy, "_MOTIF_WM_HINTS", False);

        struct MwmHints hints;
        memset(&hints, 0, sizeof(hints));
        hints.flags = MWM_HINTS_DECORATIONS | MWM_HINTS_FUNCTIONS;

        // If this window is a GUI window then enable all window decorations and
        // manipulations, except for the window close button, which would wreak havoc:
        if (windowRecord->specialflags & kPsychGUIWindow) {
            hints.decorations = MWM_FUNC_ALL;
            hints.functions = MWM_FUNC_RESIZE | MWM_FUNC_MOVE | MWM_FUNC_MINIMIZE | MWM_FUNC_MAXIMIZE;
        }

        XChangeProperty(dpy, win, mwmHintsProperty, mwmHintsProperty, 32, PropModeReplace, (unsigned char *)&hints, sizeof(hints) / sizeof(long));

        // For windowLevels of at least 500, tell window manager to try to keep
        // our window above most other windows, by setting the state to WM_STATE_ABOVE:
        if (windowLevel >= 500) {
            Atom stateAbove = XInternAtom(dpy, "_NET_WM_STATE_ABOVE", False);
            XChangeProperty(dpy, win, XInternAtom(dpy, "_NET_WM_STATE", False), XA_ATOM, 32, PropModeAppend, (unsigned char *)&stateAbove, 1);
        }
    }
            
    if (useX11) XSync(dpy, False);

    // Show our new window:
    waffle_window_show(window);

    if (useX11) {
        // Spin-Wait for it to be really mapped:
        while (0) {
            XEvent ev;
            XNextEvent(dpy, &ev);
            if (ev.type == MapNotify) break;
            PsychYieldIntervalSeconds(0.001);
        }
    
        XSync(dpy, False);

        // If windowLevel is zero, lower it to the bottom of the stack of windows:
        if (windowLevel <= 0) XLowerWindow(dpy, win);

        // Setup window transparency for user input (keyboard and mouse events):
        if (windowLevel < 1500) {
            // Need to try to be transparent for keyboard events and mouse clicks:
            XSetInputFocus(dpy, PointerRoot, RevertToPointerRoot, CurrentTime);
        }

        // Is this a non-GUI window?
        if (!(windowRecord->specialflags & kPsychGUIWindow)) {
            // Yes. As Waffle didn't override_redirect it during creation and mapping, so the
            // WM could pick up our special window properties, the WM very likely misplaced
            // it on the screen at a very ergonomic location for desktop users which doesn't
            // suit our needs. Let's fix this:

            // Wait for pending ops to complete...
            XSync(dpy, False);

            // First we override_redirect it to lock out the WM from further manipulations:
            attr.override_redirect = 1;
            mask = CWOverrideRedirect;
            XChangeWindowAttributes(dpy, win, mask, &attr);

            // Wait for override to complete...
            XSync(dpy, False);
        }

        // Then we move it to its proper location. This step is needed with Waffle even on
        // GUI windows, as Waffle has no way for us to specify window (x,y) start location:
        XMoveWindow(dpy, win, x, y);
        
        // Make sure it reaches its target position:
        XSync(dpy, False);    
    }

    // Ok, the onscreen window is ready on the screen. Time for OpenGL setup...

    // Activate the associated rendering context:
    waffle_make_current(wdpy, window, ctx);

    // Ok, the OpenGL rendering context is up and running. Auto-detect and bind all
    // available OpenGL extensions via GLEW. Must be careful to only call GLX independent
    // init code if we are not using the X11/GLX backend:
    glerr = (useGLX) ? glewInit() : glewContextInit();
    if (GLEW_OK != glerr) {
        /* Problem: glewInit failed, something is seriously wrong. */
        if (PsychPrefStateGet_Verbosity() > 0) printf("\nPTB-ERROR[GLEW init failed: %s]: Please report this to the forum. Will try to continue, but may crash soon!\n\n", glewGetErrorString(glerr));
        fflush(NULL);
    } else {
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Using GLEW version %s for automatic detection of OpenGL extensions...\n", glewGetString(GLEW_VERSION));
    }

    fflush(NULL);

    // Increase our own open window counter:
    open_windowcount++;

    // Disable X-Windows screensavers:
    if (useX11 && (open_windowcount == 1)) {
        XSync(dpy, False);

        // First window. Disable future use of screensaver:
        XSetScreenSaver(dpy, 0, 0, DefaultBlanking, DefaultExposures);

        // If the screensaver is currently running, forcefully shut it down:
        XForceScreenSaver(dpy, ScreenSaverReset);

        XSync(dpy, False);
    }

    // Some info for the user regarding non-fullscreen mode and sync problems:
    if (!(windowRecord->specialflags & kPsychIsFullscreenWindow) && (PsychPrefStateGet_Verbosity() > 2)) {
        printf("PTB-INFO: Many graphics cards do not support proper timing and timestamping of visual stimulus onset\n");
        printf("PTB-INFO: when running in windowed mode (non-fullscreen). If PTB aborts with 'Synchronization failure'\n");
        printf("PTB-INFO: you can disable the sync test via call to Screen('Preference', 'SkipSyncTests', 2); .\n");
        printf("PTB-INFO: You won't get proper stimulus onset timestamps though, so windowed mode may be of limited use.\n");
    }
    fflush(NULL);

    // Check for availability of VSYNC extension:
    if (useX11 && useGLX) {
        // First we try if the MESA variant of the swap control extensions is available. It has two advantages:
        // First, it also provides a function to query the current swap interval. Second it allows to set a
        // zero swap interval to dynamically disable sync to retrace, just as on OS/X and Windows:
        if (strstr(glXQueryExtensionsString(dpy, scrnum), "GLX_MESA_swap_control")) {
            // Bingo! Bind Mesa variant of setup call to sgi setup call, just to simplify the code
            // that actually uses the setup call -- no special cases or extra code needed there :-)
            // This special glXSwapIntervalSGI() call will simply accept an input value of zero for
            // disabling vsync'ed bufferswaps as a valid input parameter:
            glXSwapIntervalSGI = (PFNGLXSWAPINTERVALSGIPROC) glXGetProcAddressARB("glXSwapIntervalMESA");

            // Additionally bind the Mesa query call:
            glXGetSwapIntervalMESA = (PFNGLXGETSWAPINTERVALMESAPROC) glXGetProcAddressARB("glXGetSwapIntervalMESA");
            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Using GLX_MESA_swap_control extension for control of vsync.\n");
        } else {
            // Unsupported. Disable the get call:
            glXGetSwapIntervalMESA = NULL;
        }

        // Special case: Buggy ATI driver: Supports the VSync extension and glXSwapIntervalSGI, but provides the
        // wrong extension namestring "WGL_EXT_swap_control" (from MS-Windows!), so GLEW doesn't auto-detect and
        // bind the extension. If this special case is present, we do it here manually ourselves:
        if ((glXSwapIntervalSGI == NULL) && (strstr(glGetString(GL_EXTENSIONS), "WGL_EXT_swap_control") != NULL)) {
            // Looks so: Bind manually...
            glXSwapIntervalSGI = (PFNGLXSWAPINTERVALSGIPROC) glXGetProcAddressARB("glXSwapIntervalSGI");
        }

        // Extension finally supported?
        if (glXSwapIntervalSGI == NULL || (strstr(glXQueryExtensionsString(dpy, scrnum), "GLX_SGI_swap_control") == NULL &&
                                           strstr(glGetString(GL_EXTENSIONS), "WGL_EXT_swap_control") == NULL && strstr(glXQueryExtensionsString(dpy, scrnum), "GLX_MESA_swap_control") == NULL)) {
            // No, total failure to bind extension:
            glXSwapIntervalSGI = NULL;
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-WARNING: Your graphics driver doesn't allow me to control syncing wrt. vertical retrace!\n");
                printf("PTB-WARNING: Please update your display graphics driver as soon as possible to fix this.\n");
                printf("PTB-WARNING: Until then, you can manually enable syncing to VBL somehow in a manner that is\n");
                printf("PTB-WARNING: dependent on the type of gfx-card and driver. Google is your friend...\n");
            }
        }
        fflush(NULL);
    }
    else if (!egl_display && !useGLX && (PsychPrefStateGet_Verbosity() > 1)) {
        printf("PTB-WARNING: The current display backend does not allow me to control synchronization of bufferswap to vertical retrace.\n"); 
    }

    // First opened onscreen window? If so, we try to map GPU MMIO registers
    // to enable beamposition based timestamping and other special goodies:
    if (open_windowcount == 1) PsychScreenMapRadeonCntlMemory();

    // Ok, we should be ready for OS independent setup...
    fflush(NULL);

    // Wait for X-Server to settle...
    if (useX11) XSync(dpy, False);

    // Wait 250 msecs extra to give desktop compositor a chance to settle:
    PsychYieldIntervalSeconds(0.25);

    // Retrieve modeline of current video mode on primary crtc for the screen to which
    // this onscreen window is assigned. Could also query useful info about crtc, but let's not
    // overdo it in the first iteration...
    XRRCrtcInfo *crtc_info = NULL;
    XRRModeInfo *mode = PsychOSGetModeLine(screenSettings->screenNumber, 0, &crtc_info);
    if (mode) {
        // Assign modes display height aka vactive or vdisplay as startline of vblank interval:
        windowRecord->VBL_Startline = mode->height;

        // Assign vbl endline as vtotal - 1:
        windowRecord->VBL_Endline = mode->vTotal - 1;

        // Check for output display rotation enabled. Will likely impair timing/timestamping
        // because it uses copy-swaps via an intermediate shadow framebuffer to do rotation
        // during copy-swap blit, instead of via rotated crtc scanout, as most crtc's don't
        // support this in hardware:
        if ((crtc_info->rotation != RR_Rotate_0) && (PsychPrefStateGet_Verbosity() > 1)) {
            printf("PTB-WARNING: Your primary output display has hardware rotation enabled. It is not displaying in upright orientation.\n");
            printf("PTB-WARNING: On many graphics cards, this will cause unreliable stimulus presentation timing and timestamping.\n");
            printf("PTB-WARNING: If you want non-upright stimulus presentation, look at 'help PsychImaging' on how to achieve this in\n");
            printf("PTB-WARNING: a way that doesn't impair timing. The subfunctions 'FlipHorizontal' and 'FlipVertical' are what you probably need.\n");
        }

        XRRFreeCrtcInfo(crtc_info);
    }

    if (useX11) XSync(dpy, False);

    // Try to enable swap event delivery to us:
    if (PsychOSSwapCompletionLogging(windowRecord, 2, 0) && (PsychPrefStateGet_Verbosity() > 3)) {
        printf("PTB-INFO: INTEL_swap_event support for additional swap completion correctness checks enabled.\n");
    }

    if (useX11) XSync(dpy, False);

    // OpenGL-ES bindings needed?
    if (PsychIsGLES(windowRecord)) {
        // Yes: Allow userspace to force rebinding:
        psych_bool forceRebind = (getenv("PSYCH_FORCE_REBINDGLES")) ? TRUE : FALSE;

        // Which version of OpenGL-ES? Set proper selector:
        int32_t apidl = 0;

        if (windowRecord->glApiType < 20) {
            apidl = WAFFLE_DL_OPENGL_ES1;
        }
        else if (windowRecord->glApiType < 30) {
            apidl = WAFFLE_DL_OPENGL_ES2;
        }
        else if (windowRecord->glApiType < 40) {
            apidl = WAFFLE_DL_OPENGL_ES3;
        }
        else printf("PTB-WARNING: Unsupported OpenGL-ES API version >= 4.0! Update the code, this may end badly!!\n");

        // Some diagnostics for forceRebind mode:
        if (forceRebind) printf("PTB-DEBUG: Forced rebinding of OpenGL-ES extensions requested.\n");
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Pre-GLES-Rebind: glBindFramebufferEXT = %p , glBindFramebuffer = %p \n",
                                                      (void*) glBindFramebufferEXT, (void*) waffle_dl_sym(WAFFLE_DL_OPENGL_ES2, "glBindFramebuffer"));

        // Try to rebind OpenGL-ES specific entry points to standard desktop OpenGL entry
        // points of compatible syntax and semantics, so we don't need to rewrite lots of
        // desktop GL support code just to account for syntactic sugar.

        // Needed in imaging pipeline: ES implementations don't expose the ARB suffixed version anymore,
        // as glActiveTexture() is part of core-spec:
        if (NULL == glActiveTextureARB || forceRebind) glActiveTextureARB = waffle_dl_sym(apidl, "glActiveTexture");

        // Fun with NVidia ES implementation, which defines glOrthof() but not glOrthofOES() on their
        // desktop drivers for GeForce and Quadro and their Tegra drivers for embedded.
        // It is crucial to force rebind here, especially on embedded, as the entry point might
        // be already bound to a dead zombie desktop GL implementation lib:
        glOrthofOES = waffle_dl_sym(WAFFLE_DL_OPENGL_ES1, "glOrthofOES");
        if (NULL == glOrthofOES) glOrthofOES = waffle_dl_sym(WAFFLE_DL_OPENGL_ES1, "glOrthof");
        if (NULL == glOrthofOES) {
            printf("PTB-ERROR: NO glOrthofOES() or glOrthof() available under OpenGL-ES api! This will not work with OpenGL-ES! Aborting.\n");
            return(FALSE);
        }

        // Bind glBlitFrameBuffer() if provided by extension or by core OpenGL-ES3:
        if (strstr(glGetString(GL_EXTENSIONS), "GL_NV_framebuffer_blit") || (windowRecord->glApiType >= 30)) {
            glBlitFramebufferEXT = waffle_dl_sym(apidl, "glBlitFramebufferNV");
            if (NULL == glBlitFramebufferEXT) glBlitFramebufferEXT = waffle_dl_sym(WAFFLE_DL_OPENGL_ES3, "glBlitFramebuffer");
            if (NULL == glBlitFramebufferEXT) printf("PTB-ERROR: Could not bind glBlitFramebuffer() function as expected! This may end in a crash!\n");
        }

        // Framebuffer objects supported? As extension on ES-1, core on ES-2+
        if (strstr(glGetString(GL_EXTENSIONS), "framebuffer_object") || (windowRecord->glApiType >= 20)) {
            // FBO extension or core code binding:
            // Some of this is a mess, as apparently each ES implementation vendor does its own non-standard thing. In theory it should be
            // enough to query for the "OES" suffixed entry points on the api as defined by apidl. In practice, three different
            // implementations (MESA, NVidia binary desktop GL drivers, NVidia binary Tegra drivers) required three different types
            // of hacks to get it working. We try the following strategy:

            // First try to bind ES-1 standard compliant: OES suffix and apidl as selected. Then retry without OES suffix with higher
            // apidl versions:
            if (apidl == WAFFLE_DL_OPENGL_ES1) {
                // OpenGL-ES 1: Framebuffer objects are only available as extension, so entry points should be OES-suffixed. Try with OES suffix:
                if (NULL == glBindFramebufferEXT || forceRebind) glBindFramebufferEXT = waffle_dl_sym(apidl, "glBindFramebufferOES");
                if (NULL == glDeleteFramebuffersEXT || forceRebind) glDeleteFramebuffersEXT = waffle_dl_sym(apidl, "glDeleteFramebuffersOES");
                if (NULL == glGenFramebuffersEXT || forceRebind) glGenFramebuffersEXT = waffle_dl_sym(apidl, "glGenFramebuffersOES");
                if (NULL == glIsFramebufferEXT || forceRebind) glIsFramebufferEXT = waffle_dl_sym(apidl, "glIsFramebufferOES");
                if (NULL == glCheckFramebufferStatusEXT || forceRebind) glCheckFramebufferStatusEXT = waffle_dl_sym(apidl, "glCheckFramebufferStatusOES");

                if (NULL == glFramebufferTexture2DEXT || forceRebind) glFramebufferTexture2DEXT = waffle_dl_sym(apidl, "glFramebufferTexture2DOES");
                if (NULL == glFramebufferRenderbufferEXT || forceRebind) glFramebufferRenderbufferEXT = waffle_dl_sym(apidl, "glFramebufferRenderbufferOES");
                if (NULL == glGetFramebufferAttachmentParameterivEXT || forceRebind) glGetFramebufferAttachmentParameterivEXT = waffle_dl_sym(apidl, "glGetFramebufferAttachmentParameterivOES");
                if (NULL == glGenerateMipmapEXT || forceRebind) glGenerateMipmapEXT = waffle_dl_sym(apidl, "glGenerateMipmapOES");

                if (NULL == glIsRenderbufferEXT || forceRebind) glIsRenderbufferEXT = waffle_dl_sym(apidl, "glIsRenderbufferOES");
                if (NULL == glBindRenderbufferEXT || forceRebind) glBindRenderbufferEXT = waffle_dl_sym(apidl, "glBindRenderbufferOES");
                if (NULL == glDeleteRenderbuffersEXT || forceRebind) glDeleteRenderbuffersEXT = waffle_dl_sym(apidl, "glDeleteRenderbuffersOES");
                if (NULL == glGenRenderbuffersEXT || forceRebind) glGenRenderbuffersEXT = waffle_dl_sym(apidl, "glGenRenderbuffersOES");
                if (NULL == glRenderbufferStorageEXT || forceRebind) glRenderbufferStorageEXT = waffle_dl_sym(apidl, "glRenderbufferStorageOES");
                if (NULL == glGetRenderbufferParameterivEXT || forceRebind) glGetRenderbufferParameterivEXT = waffle_dl_sym(apidl, "glGetRenderbufferParameterivOES");
            }

            // Try binding without OES suffix: We do this for OpenGL-ES2 and later, as FBO's are part of the core,
            // and on OpenGL-ES1 if the "OES" binding failed, e.g., on the non standard compliant NVidia Tegra binary drivers:
            if ((apidl != WAFFLE_DL_OPENGL_ES1) || (NULL == glBindFramebufferEXT)) {
                // NVidia weirdness?
                if (apidl == WAFFLE_DL_OPENGL_ES1) {
                    // Yes. Their Tegra OpenGL-ES1 library internally loads their OpenGL-ES2 library behind our back
                    // and then uses their ES-2 functionality to emulate ES-1 behaviour, ie., their ES-1 lib is a wrapper
                    // around ES-2 libs, using self-defined shaders etc. to replicate fixed-pipeline behaviour. Unfortunately
                    // their emulation is not perfectly transparent, as their ES-1 libs claim FBO support, but get FBO entry
                    // points from ES-2 libs, thereby their entry points don't have the ES-1 compliant OES suffixes and our
                    // binding procedure for ES-1 above fails miserably. What to do? We query their ES-2 implementation/libs
                    // directly and bind ES2 entry points and hope this works and doesn't die a miserable death:
                    apidl = WAFFLE_DL_OPENGL_ES2;

                    // Of course they may pull the same trick with ES-2 vs. ES-3 some time, so if their ES-2 libs don't
                    // expose the most basic FBO entry point, we try with ES3:
                    if (NULL == waffle_dl_sym(apidl, "glBindFramebuffer")) apidl = WAFFLE_DL_OPENGL_ES3;
                }

                if (NULL == glBindFramebufferEXT || forceRebind) glBindFramebufferEXT = waffle_dl_sym(apidl, "glBindFramebuffer");
                if (NULL == glDeleteFramebuffersEXT || forceRebind) glDeleteFramebuffersEXT = waffle_dl_sym(apidl, "glDeleteFramebuffers");
                if (NULL == glGenFramebuffersEXT || forceRebind) glGenFramebuffersEXT = waffle_dl_sym(apidl, "glGenFramebuffers");
                if (NULL == glIsFramebufferEXT || forceRebind) glIsFramebufferEXT = waffle_dl_sym(apidl, "glIsFramebuffer");
                if (NULL == glCheckFramebufferStatusEXT || forceRebind) glCheckFramebufferStatusEXT = waffle_dl_sym(apidl, "glCheckFramebufferStatus");

                if (NULL == glFramebufferTexture2DEXT || forceRebind) glFramebufferTexture2DEXT = waffle_dl_sym(apidl, "glFramebufferTexture2D");
                if (NULL == glFramebufferRenderbufferEXT || forceRebind) glFramebufferRenderbufferEXT = waffle_dl_sym(apidl, "glFramebufferRenderbuffer");
                if (NULL == glGetFramebufferAttachmentParameterivEXT || forceRebind) glGetFramebufferAttachmentParameterivEXT = waffle_dl_sym(apidl, "glGetFramebufferAttachmentParameteriv");
                if (NULL == glGenerateMipmapEXT || forceRebind) glGenerateMipmapEXT = waffle_dl_sym(apidl, "glGenerateMipmap");

                if (NULL == glIsRenderbufferEXT || forceRebind) glIsRenderbufferEXT = waffle_dl_sym(apidl, "glIsRenderbuffer");
                if (NULL == glBindRenderbufferEXT || forceRebind) glBindRenderbufferEXT = waffle_dl_sym(apidl, "glBindRenderbuffer");
                if (NULL == glDeleteRenderbuffersEXT || forceRebind) glDeleteRenderbuffersEXT = waffle_dl_sym(apidl, "glDeleteRenderbuffers");
                if (NULL == glGenRenderbuffersEXT || forceRebind) glGenRenderbuffersEXT = waffle_dl_sym(apidl, "glGenRenderbuffers");
                if (NULL == glRenderbufferStorageEXT || forceRebind) glRenderbufferStorageEXT = waffle_dl_sym(apidl, "glRenderbufferStorage");
                if (NULL == glGetRenderbufferParameterivEXT || forceRebind) glGetRenderbufferParameterivEXT = waffle_dl_sym(apidl, "glGetRenderbufferParameteriv");
            }
        } // End of FBO binding.
    } // End of OpenGL-ES binding.

    return(TRUE);
}

void PsychOSCloseWindow(PsychWindowRecordType * windowRecord)
{
    Display *dpy = windowRecord->targetSpecific.privDpy;

    // Check if we are trying to close the window after it had an "odd" (== non-even)
    // number of bufferswaps. If so, we execute one last bufferswap to make the count
    // even. This means that if this window was swapped via page-flipping, the system
    // should end with the same backbuffer-frontbuffer assignment as the one prior
    // to opening the window. This may help sidestep certain bugs in compositing desktop
    // managers (e.g., Compiz).
    /*    if (PsychOSGetPostSwapSBC(windowRecord) % 2) {
        // Uneven count. Submit a swapbuffers request and wait for it to truly finish:

        // We have to rebind the OpenGL context for this swapbuffers call to work around some
        // mesa bug for intel drivers which would cause a crash without context:
        glXMakeCurrent(dpy, windowRecord->targetSpecific.windowHandle, windowRecord->targetSpecific.contextObject);

        PsychOSFlipWindowBuffers(windowRecord);
        PsychOSGetPostSwapSBC(windowRecord);
    }

    if (PsychPrefStateGet_Verbosity() > 5) {
        printf("PTB-DEBUG:PsychOSCloseWindow: Closing with a final swapbuffers count of %i.\n", (int)PsychOSGetPostSwapSBC(windowRecord));
    }
    */

    if (useX11) XSync(dpy, False);

    // Detach OpenGL rendering context again - just to be safe!
    waffle_make_current(windowRecord->targetSpecific.deviceContext, NULL, NULL);
    currentContext = NULL;

    // Delete rendering context:
    waffle_context_destroy(windowRecord->targetSpecific.contextObject);
    windowRecord->targetSpecific.contextObject = NULL;

    // Delete swap context:
    waffle_context_destroy(windowRecord->targetSpecific.glswapcontextObject);
    windowRecord->targetSpecific.glswapcontextObject = NULL;

    // Delete userspace context, if any:
    if (windowRecord->targetSpecific.glusercontextObject) {
        waffle_context_destroy(windowRecord->targetSpecific.glusercontextObject);
        windowRecord->targetSpecific.glusercontextObject = NULL;
    }

    // Close & Destroy the window:
    waffle_window_destroy(windowRecord->targetSpecific.windowHandle);
    windowRecord->targetSpecific.windowHandle = 0;
    windowRecord->targetSpecific.xwindowHandle = 0;

    // Decrement global count of open onscreen windows:
    open_windowcount--;

    // Was this the last window?
    if (open_windowcount <= 0) {
        open_windowcount = 0;

        // Release our shared waffle display connection:
        waffle_display_disconnect(windowRecord->targetSpecific.deviceContext);
        wdpy = NULL;

        if (useX11) XSync(dpy, False);

        // (Re-)enable X-Windows screensavers if they were enabled before opening windows:
        // Set screensaver to previous settings, potentially enabling it:
        if (useX11) XSetScreenSaver(dpy, -1, 0, DefaultBlanking, DefaultExposures);

        // Unmap/release possibly mapped device memory: Defined in PsychScreenGlue.c
        PsychScreenUnmapDeviceMemory();
    }

    windowRecord->targetSpecific.deviceContext = NULL;

    // Done.
    return;
}

/*
    PsychOSGetVBLTimeAndCount()

    Returns absolute system time of last VBL and current total count of VBL interrupts since
    startup of gfx-system for the given screen. Returns a time of -1 and a count of 0 if this
    feature is unavailable on the given OS/Hardware configuration.
*/
double PsychOSGetVBLTimeAndCount(PsychWindowRecordType *windowRecord, psych_uint64* vblCount)
{
    // Unsupported:
    *vblCount = 0;
    return(-1);
}


/* PsychOSGetSwapCompletionTimestamp()
 *
 * Retrieve a very precise timestamp of doublebuffer swap completion by means
 * of OS specific facilities. This function is optional. If the underlying
 * OS/driver/GPU combo doesn't support a high-precision, high-reliability method
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
 * Returns: Precise and reliable swap completion timestamp in seconds of
 * system time in variable referenced by 'tSwap', and msc value of completed swap,
 * or a negative value on error (-1 == unsupported, -2/-3 == Query failed).
 *
 */
psych_int64 PsychOSGetSwapCompletionTimestamp(PsychWindowRecordType *windowRecord, psych_int64 targetSBC, double* tSwap)
{
	// Return -1 "unsupported":
	return(-1);
}

/* PsychOSInitializeOpenML() - Linux specific function.
 *
 * Performs basic initialization of the OpenML OML_sync_control extension.
 * Performs basic and extended correctness tests and disables extension if it
 * is unreliable, or enables workarounds for partially broken extensions.
 *
 * Called from PsychDetectAndAssignGfxCapabilities() as part of the PsychOpenOffscreenWindow()
 * procedure for a window with OpenML support.
 *
 */
void PsychOSInitializeOpenML(PsychWindowRecordType *windowRecord)
{
    // Unsupported:

    // Disable OpenML swap scheduling:
    windowRecord->gfxcaps &= ~kPsychGfxCapSupportsOpenML;

    // OpenML timestamping in PsychOSGetSwapCompletionTimestamp() and PsychOSGetVBLTimeAndCount() disabled:
    windowRecord->specialflags |= kPsychOpenMLDefective;

    return;
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
	specialFlags - Additional options, a bit field consisting of single bits that can be or'ed together:
				   1 = Constrain swaps to even msc values, 2 = Constrain swaps to odd msc values. (Used for frame-seq. stereo field selection)
	
	Return value:
	 
	Value greater than or equal to zero on success: The target msc for which swap is scheduled.
	Negative value: Error. Function failed. -1 == Function unsupported on current system configuration.
	-2 ... -x == Error condition.
	
*/
psych_int64 PsychOSScheduleFlipWindowBuffers(PsychWindowRecordType *windowRecord, double tWhen, psych_int64 targetMSC, psych_int64 divisor, psych_int64 remainder, unsigned int specialFlags)
{
    // Unsupported:
    return(-1);
}

/*
  PsychOSFlipWindowBuffers()
    
  Performs OS specific double buffer swap call.
*/
void PsychOSFlipWindowBuffers(PsychWindowRecordType *windowRecord)
{
    // Execute OS neutral bufferswap code first:
    PsychExecuteBufferSwapPrefix(windowRecord);

    // Trigger the "Front <-> Back buffer swap (flip) (on next vertical retrace)":
    waffle_window_swap_buffers(windowRecord->targetSpecific.windowHandle);
    windowRecord->target_sbc = 0;

    return;
}

/* Enable/disable syncing of buffer-swaps to vertical retrace. */
void PsychOSSetVBLSyncLevel(PsychWindowRecordType *windowRecord, int swapInterval)
{
    int error, myinterval;

    // Enable rendering context of window (no-ops internally when not called from master-thread):
    PsychSetGLContext(windowRecord);

    // Store new setting also in internal helper variable, e.g., to allow workarounds to work:
    windowRecord->vSynced = (swapInterval > 0) ? TRUE : FALSE;

    // Classix X11/GLX?
    if (useGLX) {
        XSync(windowRecord->targetSpecific.privDpy, False);

        // Try to set requested swapInterval if swap-control extension is supported on
        // this Linux machine. Otherwise this will be a no-op...
        // Note: On Mesa, glXSwapIntervalSGI() is actually a redirected call to glXSwapIntervalMESA()!
        if (glXSwapIntervalSGI) {
            error = glXSwapIntervalSGI(swapInterval);
            if (error) {
                if (PsychPrefStateGet_Verbosity() > 1) printf("\nPTB-WARNING: FAILED to %s synchronization to vertical retrace!\n\n", (swapInterval > 0) ? "enable" : "disable");
            }
        }

        // If Mesa query is supported, double-check if the system accepted our settings:
        if (glXGetSwapIntervalMESA) {
            myinterval = glXGetSwapIntervalMESA();
            if (myinterval != swapInterval) {
                if (PsychPrefStateGet_Verbosity() > 1) printf("\nPTB-WARNING: FAILED to %s synchronization to vertical retrace (System ignored setting [Req %i != Actual %i])!\n\n", (swapInterval > 0) ? "enable" : "disable", swapInterval, myinterval);
            }
        }
        XSync(windowRecord->targetSpecific.privDpy, False);
    }
    else {
        // EGL-based backend in use:
        if (egl_display) {
        #ifdef PTB_USE_EGL
            if (!eglSwapInterval(egl_display, (EGLint) swapInterval)) {
                if (PsychPrefStateGet_Verbosity() > 1) printf("\nPTB-WARNING: FAILED to %s synchronization to vertical retrace!\n\n", (swapInterval > 0) ? "enable" : "disable");
            }
        #endif
        }
    }

    return;
}

/*
  PsychOSSetGLContext()
    
  Set the window to which GL drawing commands are sent.  
*/
void PsychOSSetGLContext(PsychWindowRecordType * windowRecord)
{
    // Only change context if not already proper context bound:
    // This is very important not only for performance, but also to avoid harmful
    // glBindFrameBufferEXT(0) calls without switching away from context -- something
    // that would completely mess up imaging pipeline state!
    if (currentContext != windowRecord->targetSpecific.contextObject) {
        if (currentContext) {
            // We need to glFlush the context before switching, otherwise race-conditions may occur:
            glFlush();

            // Need to unbind any FBO's in old context before switch, otherwise bad things can happen...
            if (glBindFramebufferEXT) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
        }

        if (useX11) XSync(windowRecord->targetSpecific.privDpy, False);

        // Switch to new context: Skip surface bind if forbidden by higher-level code, e.g., for multi-threaded EGL use.
        if (!waffle_make_current(windowRecord->targetSpecific.deviceContext,
                                 (PsychGetParentWindow(windowRecord)->specialflags & kPsychSurfacelessContexts) ? EGL_NO_SURFACE : windowRecord->targetSpecific.windowHandle,
                                 windowRecord->targetSpecific.contextObject) && (PsychPrefStateGet_Verbosity() > 0)) {
            printf("PTB-ERROR: PsychOSSetGLContext(): Failed to bind master context of window %i: [%s] - Expect trouble!\n", windowRecord->windowIndex, waffle_error_to_string(waffle_error_get_code()));
        }
        else {
            // Update context tracking:
            currentContext = windowRecord->targetSpecific.contextObject;
        }
    }

    return;
}

/*
  PsychOSUnsetGLContext()
    
  Clear the drawing context.  
*/
void PsychOSUnsetGLContext(PsychWindowRecordType * windowRecord)
{
    if (PsychIsMasterThread()) {
        if (currentContext) {
            // We need to glFlush the context before switching, otherwise race-conditions may occur:
            glFlush();

            // Need to unbind any FBO's in old context before switch, otherwise bad things can happen...
            if (glBindFramebufferEXT) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
        }

        // Reset tracking state:
        currentContext = NULL;
    }

    if (useX11) XSync(windowRecord->targetSpecific.privDpy, False);

    // Detach context:
    waffle_make_current(windowRecord->targetSpecific.deviceContext, NULL, NULL);

    return;
}

/* Same as PsychOSSetGLContext() but for selecting userspace rendering context,
 * optionally copying state from PTBs context.
 */
void PsychOSSetUserGLContext(PsychWindowRecordType * windowRecord, psych_bool copyfromPTBContext)
{
    // Child protection:
    if (windowRecord->targetSpecific.glusercontextObject == NULL)
        PsychErrorExitMsg(PsychError_user, "GL Userspace context unavailable! Call InitializeMatlabOpenGL *before* Screen('OpenWindow')!");

    // Only change context if not already proper context bound:
    // This is very important not only for performance, but also to avoid harmful
    // glBindFrameBufferEXT(0) calls without switching away from context -- something
    // that would completely mess up imaging pipeline state!
    if (currentContext != windowRecord->targetSpecific.glusercontextObject) {

        if (useX11) XSync(windowRecord->targetSpecific.privDpy, False);

        // Copy from context unsupported on Waffle backends:
        if (copyfromPTBContext && (PsychPrefStateGet_Verbosity() > 1)) {
            printf("PTB-WARNING: Tried to set the 'sharecontext' flag to 2 in Screen('BeginOpenGL'), but this is not supported with Waffle display backends. Ignored!\n");
        }

        // Bind it: Skip surface bind if forbidden by higher-level code, e.g., for multi-threaded EGL use.
        if (!waffle_make_current(windowRecord->targetSpecific.deviceContext,
                                 (PsychGetParentWindow(windowRecord)->specialflags & kPsychSurfacelessContexts) ? EGL_NO_SURFACE : windowRecord->targetSpecific.windowHandle,
                                 windowRecord->targetSpecific.glusercontextObject) && (PsychPrefStateGet_Verbosity() > 0)) {
            printf("PTB-ERROR: PsychOSSetUserGLContext(): Failed to bind userspace gl-context of window %i: [%s] - Expect trouble!\n", windowRecord->windowIndex, waffle_error_to_string(waffle_error_get_code()));
        }
        else {
            // Update context tracking:
            currentContext = windowRecord->targetSpecific.glusercontextObject;
        }
    }

    return;
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
	if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: NV_swap_group and GLX_SGIX_swap_group unsupported.\n");
	return(FALSE);
}

psych_bool PsychOSSwapCompletionLogging(PsychWindowRecordType *windowRecord, int cmd, int aux1)
{
    // Failed to enable swap events, because they're unsupported:
    windowRecord->swapevents_enabled = 0;
    return(FALSE);
}

/* End of PTB_USE_WAFFLE */
#endif

