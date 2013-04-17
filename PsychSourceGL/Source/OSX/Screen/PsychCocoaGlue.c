/*
 
    PsychToolbox3/Source/OSX/Screen/PsychCocoaGlue.c

    PLATFORMS:	

    This is the OS X 64-Bit Cocoa version only.
 
    TARGETS:
 
    Screen

    AUTHORS:

    Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de
 
    DESCRIPTION:
 
    Glue code for window management for the 64-Bit Psychtoolbox,
    using Objective-C wrappers to use Cocoa.
 
    These functions are called by PsychWindowGlue.c. Some are
    specific to Cocoa window management. Some are drop-in replacements
    for missing Carbon 32-Bit functions.
 
    NOTES:
 
    The setup code for multiple OpenGL contexts per Cocoa window makes
    use of functions only supported on OSX 10.6 "Snow Leopard" and later.
    Therefore, full windowed mode functionality is only available with
    10.6 and later. 10.5 only provides restricted basic 2D/3D drawing
    contexts with hard-coded pixelformat parameters and no resource sharing
    across onscreen windows.

*/

/* Only use Cocoa glue on 64-Bit OSX for now */
#ifdef __LP64__

#include "Screen.h"
#include "PsychCocoaGlue.h"

#include <ApplicationServices/ApplicationServices.h>
#include <Cocoa/Cocoa.h>

PsychError PsychCocoaCreateWindow(PsychWindowRecordType *windowRecord,
                           PsychRectType      screenRect,
                           const Rect *       contentBounds,
                           WindowClass        wclass,
                           WindowAttributes   addAttribs,
                           WindowRef *        outWindow,
                           psych_bool         enableTransparentGL)
{
    char windowTitle[100];
    NSWindow *cocoaWindow;

    // Zero-Init NSOpenGLContext-Pointers for our private Cocoa OpenGL contexts:
    windowRecord->targetSpecific.nsmasterContext = NULL;
    windowRecord->targetSpecific.nsswapContext = NULL;
    windowRecord->targetSpecific.nsuserContext = NULL;
    
    // Include onscreen window index in title:
    sprintf(windowTitle, "PTB Onscreen Window [%i]:", windowRecord->windowIndex);
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    // Initialize the Cocoa application object, connect to CoreGraphics-Server:
    // Can be called many times, as redundant calls are ignored.
    NSApplicationLoad();

    // Define size of client area - the actual stimulus display area:
    // The window itself will resize and reposition itself so that the size of
    // the content area is preserved/honored, adjusting for the thickness of window decorations.
    NSRect windowRect = NSMakeRect(0, 0, (contentBounds->right - contentBounds->left), (contentBounds->bottom - contentBounds->top));

    NSUInteger windowStyle = 0;
    if (windowRecord->specialflags & kPsychGUIWindow) {
        // GUI window:
        windowStyle = NSTitledWindowMask|NSClosableWindowMask|NSResizableWindowMask|NSMiniaturizableWindowMask|NSTexturedBackgroundWindowMask;
    }
    else {
        // Pure non-GUI visual stimulus window:
        windowStyle = NSBorderlessWindowMask;
    }

    cocoaWindow = [[NSWindow alloc] initWithContentRect:windowRect styleMask:windowStyle backing:NSBackingStoreBuffered defer:false];
    if (cocoaWindow == nil) {
        printf("PTB-ERROR:PsychCocoaCreateWindow(): Could not create Cocoa-Window!\n");
        // Return failure:
        return(PsychError_system);
    }

    [cocoaWindow setTitle:[NSString stringWithUTF8String:windowTitle]];
    
    if (enableTransparentGL) {
        // Set window as non-opaque, with a transparent window background color:
        // This together with the OpenGL context setup for transparency allows the OpenGL
        // colorbuffer alpha-channel to determine window opacity at a per-pixel level, so
        // experimental code has full control over transpareny if it wishes so. By default,
        // Screen() clears the backbuffer to an opaque white with alpha 1.0, so by default
        // windows do appear solid, unless usercode explicitely does something else:
        [cocoaWindow setOpaque:false];
        [cocoaWindow setBackgroundColor:[NSColor colorWithDeviceWhite:0.0 alpha:0.0]];
    }
    
    // Make window "transparent" for mouse events like clicks and drags, if requested:
    if (addAttribs & kWindowIgnoreClicksAttribute) [cocoaWindow setIgnoresMouseEvents:true];

    // In non-GUI mode we want the window to be above all other regular windows, so the
    // stimulus doesn't get occluded. If we make ourselves transparent to mouse clicks, we
    // must be above all other windows, as otherwise any mouse-click that "goes through"
    // to an underlying window will raise that window above ours and we get occluded, ie.,
    // any actual passed-through mouse-click would defeat the purpose of pass-through mode:
    if (!(windowRecord->specialflags & kPsychGUIWindow) || (addAttribs & kWindowIgnoreClicksAttribute)) {
        // Set level of window to be in front of every regular window:
        [cocoaWindow setLevel:NSScreenSaverWindowLevel];
    }

    // Disable auto-flushing of drawed content to frontbuffer:
    [cocoaWindow disableFlushWindow];

    // Position the window. Origin is bottom-left of screen, as opposed to Carbon / PTB origin
    // of top-left. Therefore need to invert the vertical position. Cocoa only takes our request
    // as a hint. It tries to position as requested, but places the window differently if required
    // to make sure the full windowRect content area is displayed. It doesn't allow the window to
    // overlap the menu bar or dock area by default.
    NSPoint winPosition = NSMakePoint(contentBounds->left, screenRect[kPsychBottom] - contentBounds->top);
    [cocoaWindow setFrameTopLeftPoint:winPosition];
    
    // Query and translate content rect of final window to a PTB rect:
    NSRect clientRect = [cocoaWindow contentRectForFrameRect:[cocoaWindow frame]];
    PsychMakeRect(windowRecord->rect, clientRect.origin.x, screenRect[kPsychBottom] - (clientRect.origin.y + clientRect.size.height), clientRect.origin.x + clientRect.size.width, screenRect[kPsychBottom] - clientRect.origin.y);

    // Drain the pool:
    [pool drain];
    
    // Return window pointer, packed into an old-school Carbon window ref:
    *outWindow = (WindowRef) cocoaWindow;

    // Return success:
    return(PsychError_none);
}

OSStatus
GetWindowBounds(
                WindowRef          window,
                WindowRegionCode   regionCode,
                Rect *             globalBounds)
{
    NSWindow *cocoaWindow = (NSWindow*) window;
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    // Query and translate content rect of final window to a PTB rect:
    NSRect clientRect = [cocoaWindow contentRectForFrameRect:[cocoaWindow frame]];
    
    globalBounds->left = clientRect.origin.x;
    globalBounds->right = clientRect.origin.x + clientRect.size.width;
    NSRect screenRect = [[cocoaWindow screen] frame];
    globalBounds->top = screenRect.size.height - (clientRect.origin.y + clientRect.size.height);
    globalBounds->bottom = globalBounds->top + clientRect.size.height;

    // Drain the pool:
    [pool drain];
    
    return(0);
}

void ShowWindow(WindowRef window)
{
    PsychCocoaShowWindow(window);
}

WindowGroupRef GetWindowGroup(WindowRef inWindow)
{
    return((WindowGroupRef) inWindow);
}

pid_t GetHostingWindowsPID(void)
{
    pid_t pid = (pid_t) 0;
    CFIndex i;
    CFNumberRef numRef;
    char winName[256];
    psych_bool found = FALSE;
    psych_bool verbose = (PsychPrefStateGet_Verbosity() > 5) ? TRUE : FALSE;

    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    
    CFArrayRef windowList = CGWindowListCopyWindowInfo(kCGWindowListOptionOnScreenOnly | kCGWindowListExcludeDesktopElements, kCGNullWindowID);
    if (!windowList) goto hwinpidout;	
    
    const CFIndex kSize = CFArrayGetCount(windowList);
    
    for (i = 0; i < kSize; ++i) {
        CFDictionaryRef d = (CFDictionaryRef) CFArrayGetValueAtIndex(windowList, i);

        // Get process id pid of window owner:
        numRef = (CFNumberRef) CFDictionaryGetValue(d, kCGWindowOwnerPID);
        if (numRef) {
            int val;
            CFNumberGetValue(numRef, kCFNumberIntType, &val);
            pid = (pid_t) val;
            if (verbose) printf("OwnerPID: %i\n", val);
        }
        
        numRef = (CFNumberRef) CFDictionaryGetValue(d, kCGWindowLayer);
        if (numRef) {
            int val;
            CFNumberGetValue(numRef, kCFNumberIntType, &val);
            if (verbose) printf("WindowLevel: %i  (ShieldingWindow %i)\n", val, CGShieldingWindowLevel());
        }

        // Get window name of specific window. Rarely set by apps:
        winName[0] = 0;
        CFStringRef nameRef = (CFStringRef) CFDictionaryGetValue(d, kCGWindowName);
        if (nameRef) {
            const char* name = CFStringGetCStringPtr(nameRef, kCFStringEncodingMacRoman);
            if (name && verbose) printf("WindowName: %s\n", name);
            if (name) snprintf(winName, sizeof(winName), "%s", name);
        }

        // Get name of owner process/app:
        CFStringRef nameOwnerRef = (CFStringRef) CFDictionaryGetValue(d, kCGWindowOwnerName);
        if (nameOwnerRef) {
            const char* name = CFStringGetCStringPtr(nameOwnerRef, kCFStringEncodingMacRoman);
            if (name && verbose) printf("WindowOwnerName: %s\n", name);
            if (name && ((strstr(name, "X11") && strstr(winName, "xterm")) || strstr(name, "Terminal") || strstr(name, "MATLAB"))) {
                // Matched either X11 xterm, or a OSX native Terminal or MATLAB GUI. These are candidates for the
                // hosting windows of our matlab, matlab -nojvm or octave console session. As windows are returned
                // in front-to-back order, the first match here is a candidate window that is on top of
                // the visible window stack. This is our best candidate for the command window, assuming
                // it is frontmost as the user just interacted with it. Therefore, aborting the search
                // on the first match is the most robust heuristic i can think of, given that the name
                // strings do not contain any info if a specific window hosts our session.
                found = TRUE;
                
                // pid contains the pid of the owning process.
                break;
            }
        }        
    }
    
    CFRelease(windowList);
    
hwinpidout:
    
    // Drain the pool:
    [pool drain];

    if (found) {
        if (verbose) printf("TARGETWINDOWNAME: '%s' with pid %i.\n", winName, pid);
    }
    else pid = 0;
    
    return(pid);
}

// SetUserFocusWindow() is a drop-in replacement for Carbon's function.
OSStatus SetUserFocusWindow(WindowRef inWindow)
{
    NSWindow* focusWindow = (NSWindow*) inWindow;
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    // Special flag: Try to restore main apps focus:
    if (inWindow == kUserFocusAuto) {
        focusWindow = [[NSApplication sharedApplication] mainWindow];
    }

    // Direct keyboard input focus to window 'inWindow':
    [focusWindow makeKeyAndOrderFront: nil];

    // Special handle NULL provided? Try to regain keyboard focus rambo-style for
    // our hosting window for octave / matlab -nojvm in terminal window:
    if (inWindow == NULL) {
        // This works to give keyboard focus to a process other than our (Matlab/Octave) runtime, if
        // the process id (pid_t) of the process is known and valid for a GUI app. E.g., passing in
        // the pid of the XServer process X11.app or the Konsole.app will restore the xterm'inal windows
        // or Terminal windows keyboard focus after a CGDisplayRelease() call, and thereby to the
        // octave / matlab -nojvm process which is hosted by those windows.
        //
        // Problem: Finding the pid requires iterating and filtering over all windows and name matching for
        // all possible candidates, and a shielding window from CGDisplayCapture() will still prevent keyboard
        // input, even if the window has input focus...
        pid_t pid = GetHostingWindowsPID();
        
        // Also, the required NSRunningApplication class is unsupported on 64-Bit OSX 10.5, so we need to
        // dynamically bind it and no-op if it is unsupported:
        Class nsRunningAppClass = NSClassFromString(@"NSRunningApplication");
        
        if (pid && (nsRunningAppClass != NULL)) {
            NSRunningApplication* motherapp = [nsRunningAppClass runningApplicationWithProcessIdentifier: pid];    
            [motherapp activateWithOptions: NSApplicationActivateIgnoringOtherApps];
        }
    }

    // Drain the pool:
    [pool drain];
    
    return((OSStatus) 0);
}

// GetUserFocusWindow() is a drop-in replacement for Carbon's function.
// TODO FIXME: This totally doesn't work for some reason. It always
// returns a focusWindow of NULL.
WindowRef GetUserFocusWindow(void)
{
    NSWindow* focusWindow = NULL;
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    // Retrieve pointer to current keyWindow - the window that receives
    // key events aka the window with keyboard input focus. Or 'nil' if
    // no such window exists:
    focusWindow = [[NSApplication sharedApplication] keyWindow];
    
    // Drain the pool:
    [pool drain];

    return((WindowRef) focusWindow);
}

void PsychCocoaDisposeWindow(PsychWindowRecordType *windowRecord)
{
    NSWindow *cocoaWindow = (NSWindow*) windowRecord->targetSpecific.windowHandle;
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    // Release NSOpenGLContext's - this will also release the wrapped
    // CGLContext's and finally really destroy them:
    if (windowRecord->targetSpecific.nsmasterContext) [((NSOpenGLContext*) windowRecord->targetSpecific.nsmasterContext) release];
    if (windowRecord->targetSpecific.nsswapContext) [((NSOpenGLContext*) windowRecord->targetSpecific.nsswapContext) release];
    if (windowRecord->targetSpecific.nsuserContext) [((NSOpenGLContext*) windowRecord->targetSpecific.nsuserContext) release];

    // Zero-Out the contexts after release:
    windowRecord->targetSpecific.nsmasterContext = NULL;
    windowRecord->targetSpecific.nsswapContext = NULL;
    windowRecord->targetSpecific.nsuserContext = NULL;

    // Close window. This will also release the associated contentView:
    [cocoaWindow close];

    // Drain the pool:
    [pool drain];
    
    return;
}

void PsychCocoaShowWindow(WindowRef window)
{
    NSWindow *cocoaWindow = (NSWindow*) window;
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    // Bring to front:
    [cocoaWindow orderFrontRegardless];    

    // Show window:
    [cocoaWindow display];
    
    // Drain the pool:
    [pool drain];
    
    return;   
}

psych_bool PsychCocoaSetupAndAssignOpenGLContextsFromCGLContexts(WindowRef window, PsychWindowRecordType *windowRecord)
{
    GLint opaque = 0;
    NSOpenGLContext *masterContext = NULL;
    NSOpenGLContext *gluserContext = NULL;
    NSOpenGLContext *glswapContext = NULL;
    
    NSWindow *cocoaWindow = (NSWindow*) window;
    
    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: MacOSX 10.6+ primary window display path enabled. Should be fully functional.\n");

    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    // Enable opacity for OpenGL contexts if underlying window is opaque:
    if ([cocoaWindow isOpaque] == true) opaque = 1;
    
    // Build NSOpenGLContexts as wrappers around existing CGLContexts already
    // created in calling routine:
    masterContext = [[NSOpenGLContext alloc] initWithCGLContextObj: windowRecord->targetSpecific.contextObject];
    [masterContext setValues:&opaque forParameter:NSOpenGLCPSurfaceOpacity];
    [masterContext setView:[cocoaWindow contentView]];
    [masterContext makeCurrentContext];
    // printf("MasterContext created & Made current!\n");
    
    // Ditto for potential gl userspace rendering context:
    if (windowRecord->targetSpecific.glusercontextObject) {
        gluserContext = [[NSOpenGLContext alloc] initWithCGLContextObj: windowRecord->targetSpecific.glusercontextObject];
        [gluserContext setValues:&opaque forParameter:NSOpenGLCPSurfaceOpacity];
        [gluserContext setView:[cocoaWindow contentView]];
        // printf("GLUserContext created!\n");
    }
    
    // Ditto for potential glswapcontext for async flips and frame sequential stereo:
    if (windowRecord->targetSpecific.glswapcontextObject) {
        glswapContext = [[NSOpenGLContext alloc] initWithCGLContextObj: windowRecord->targetSpecific.glswapcontextObject];
        [glswapContext setValues:&opaque forParameter:NSOpenGLCPSurfaceOpacity];
        [glswapContext setView:[cocoaWindow contentView]];
        // printf("GLSwapContext created!\n");
    }

    // printf("Refcounts: window=%i , view=%i , mc=%i mccgl=%i sc=%i sccgl=%i\n", [cocoaWindow retainCount], [[cocoaWindow contentView] retainCount], [masterContext retainCount], CGLGetContextRetainCount(windowRecord->targetSpecific.contextObject), [glswapContext retainCount], CGLGetContextRetainCount(windowRecord->targetSpecific.glswapcontextObject));
    
    // Assign contexts for use in window close sequence later on:
    windowRecord->targetSpecific.nsmasterContext = (void*) masterContext;
    windowRecord->targetSpecific.nsswapContext = (void*) glswapContext;
    windowRecord->targetSpecific.nsuserContext = (void*) gluserContext;

    // Drain the pool:
    [pool drain];
    
    // Return success:
    return(false);
}


psych_bool PsychCocoaSetupAndAssignLegacyOpenGLContext(WindowRef window, PsychWindowRecordType *windowRecord)
{
    GLint opaque = 0;
    NSOpenGLContext *masterContext = NULL;
    NSOpenGLContext *gluserContext = NULL;
    NSOpenGLContext *glswapContext = NULL;
    
    NSWindow *cocoaWindow = (NSWindow*) window;
    
    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: MacOSX 10.5 legacy window display mode for backward compatibility enabled. Functionality will be limited.\n");
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    // Enable opacity for OpenGL contexts if underlying window is opaque:
    if ([cocoaWindow isOpaque] == true) opaque = 1;

    // Define a pixelformat for the context. We use a hard-coded "one size fits all"
    // format which is used in legacy mode for OSX 10.5 and covers the most common use cases:
    NSOpenGLPixelFormatAttribute attrs[] = {
        NSOpenGLPFADoubleBuffer,
        NSOpenGLPFADepthSize, 24,
        NSOpenGLPFAStencilSize, 8,
        0
    };
    
    NSOpenGLPixelFormat* cocoaPixelFormat = [[NSOpenGLPixelFormat alloc] initWithAttributes:attrs];
    if (cocoaPixelFormat == nil) {
        printf("PTB-ERROR: Could not create NSOpenGLPixelFormat!\n");
        return(true);
    }
    
    // Build NSOpenGLContexts as wrappers around existing CGLContexts already
    // created in calling routine:
    masterContext = [[NSOpenGLContext alloc] initWithFormat: cocoaPixelFormat shareContext:nil];
    [masterContext setValues:&opaque forParameter:NSOpenGLCPSurfaceOpacity];
    [masterContext setView:[cocoaWindow contentView]];
    [masterContext makeCurrentContext];
    windowRecord->targetSpecific.contextObject = [masterContext CGLContextObj];
    CGLRetainContext(windowRecord->targetSpecific.contextObject);
    // printf("Legacy:MasterContext created & Made current!\n");
    
    // Ditto for potential glswapcontext for async flips and frame sequential stereo:
    glswapContext = [[NSOpenGLContext alloc] initWithFormat: cocoaPixelFormat shareContext:masterContext];
    [glswapContext setValues:&opaque forParameter:NSOpenGLCPSurfaceOpacity];
    [glswapContext setView:[cocoaWindow contentView]];
    windowRecord->targetSpecific.glswapcontextObject = [glswapContext CGLContextObj];
    CGLRetainContext(windowRecord->targetSpecific.glswapcontextObject);
    // printf("Legacy:GLSwapContext created!\n");
    
    // Ditto for potential gl userspace rendering context:
    gluserContext = [[NSOpenGLContext alloc] initWithFormat: cocoaPixelFormat shareContext:masterContext];
    [gluserContext setValues:&opaque forParameter:NSOpenGLCPSurfaceOpacity];
    [gluserContext setView:[cocoaWindow contentView]];
    windowRecord->targetSpecific.glusercontextObject = [gluserContext CGLContextObj];
    CGLRetainContext(windowRecord->targetSpecific.glusercontextObject);
    // printf("Legacy:GLUserContext created!\n");
    
    // Don't need the pixelformat object anymore:
    [cocoaPixelFormat release];
    
    // printf("Legacy:Refcounts: window=%i , view=%i , mc=%i mccgl=%i sc=%i sccgl=%i\n", [cocoaWindow retainCount], [[cocoaWindow contentView] retainCount], [masterContext retainCount],CGLGetContextRetainCount(windowRecord->targetSpecific.contextObject), [glswapContext retainCount], CGLGetContextRetainCount(windowRecord->targetSpecific.glswapcontextObject));
           
    // Assign contexts for use in window close sequence later on:
    windowRecord->targetSpecific.nsmasterContext = (void*) masterContext;
    windowRecord->targetSpecific.nsswapContext = (void*) glswapContext;
    windowRecord->targetSpecific.nsuserContext = (void*) gluserContext;

    // Drain the pool:
    [pool drain];

    return(false);
}

void SendBehind(WindowRef   window,
                WindowRef   behindWindow)
{
    NSWindow *cocoaWindow = (NSWindow*) window;
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    
    // Move window behind all others:
    [cocoaWindow orderBack:nil];
    
    // Drain the pool:
    [pool drain];
    
    return;    
}

OSStatus SetWindowGroupLevel(
                             WindowGroupRef   inGroup,
                             SInt32           inLevel)
{
    NSWindow *cocoaWindow = (NSWindow*) inGroup;
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    
    // Set level of window:
    [cocoaWindow setLevel:inLevel];
    
    // Drain the pool:
    [pool drain];
    
    return 0;
}

OSStatus SetWindowAlpha(
                        WindowRef   inWindow,
                        CGFloat     inAlpha)
{
    NSWindow *cocoaWindow = (NSWindow*) inWindow;
    
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    
    // Set Window transparency:
    [cocoaWindow setAlphaValue: (float) inAlpha];
    
    // Drain the pool:
    [pool drain];
    
    return(0);
}

OSStatus SetThemeCursor(ThemeCursor inCursor)
{
    // Allocate auto release pool:
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];

    switch(inCursor) {
        case kThemeArrowCursor:
            [[NSCursor arrowCursor] set];
        break;
            
        case kThemeIBeamCursor:
            [[NSCursor IBeamCursor] set];
        break;

        case kThemeCrossCursor:
            [[NSCursor crosshairCursor] set];
        break;

        case kThemePointingHandCursor:
            [[NSCursor pointingHandCursor] set];
        break;

        default:
            // Failed - Unknown cursor type:
            return(1);
    }

    // Drain the pool:
    [pool drain];
    
    return(0);
}

#endif
