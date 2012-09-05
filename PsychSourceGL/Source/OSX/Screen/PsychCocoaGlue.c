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
