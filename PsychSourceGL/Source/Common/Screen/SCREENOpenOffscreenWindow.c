/*

	SCREENOpenOffscreenWindow.c		

	AUTHORS:

		Allen.Ingling@nyu.edu               awi 
		mario.kleiner at tuebingen.mpg.de   mk
  
	PLATFORMS:	

        All.
		
	HISTORY:

		1/18/02		awi		Created.  Copied the Synopsis string from old version of psychtoolbox.
		2/18/03		awi		Set the screen number for all offscreen windows to -1.   
							What use does the screen number serve for offscren windows ?
							¥Because we use glDrawPixels and glReadPixels to move pixels to and from surfaces, and those handle 
							alignment on the surface, knowing the parent window is not necessary for the purpose of matching 
							alignments in copy operations. 
							¥The bit depths do need to match for CopyWindow, and we could use the parent window as an indicator of 
							matching depth.  However, if we want to match depths, it is better to compare depths than to use 
							the parent window number as a proxy for depth; It is more transparent and more general because it handles copying 			
							between two onscreen windows of equal depth but which do not share the same parent window.

		10/12/04	awi		In useString: moved commas to inside [].
		2/15/05		awi		Commented out glEnable(GL_BLEND) 
		1/16/06     mk      Rewritten to initialize offscreen windows as if they were textures with a constant colored background instead of content.
		1/07/07		mk		Setup path for FBO backed Offscreen windows added: Used if imaging pipeline is active.
		3/20/11		mk		Make 64-bit clean.

	TO DO:

*/

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] =  "[windowPtr,rect]=Screen('OpenOffscreenWindow',windowPtrOrScreenNumber [,color] [,rect] [,pixelSize] [,specialFlags] [,multiSample]);";
//                                                                        1                         2        3       4          5				6

static char synopsisString[] =
    "Open an offscreen window. This is simply an OpenGL texture that is treated "
    "as a window, so you can draw to it. Offscreen windows should be used "
    "to keep old code from OS-9 Psychtoolbox working, or if you need an offscreen "
	"drawing canvas for fast drawing, that you can later use as a texture. "
    "For quickly displaying premade Matlab image matrices, use the 'MakeTexture' command instead. "
    "It allow for significantly higher drawing speeds.\nYou can specify a screen "
    "(any windowPtr or a screenNumber>=0) or no screen(screenNumber=-1), but any real screen "
    "must already have an open Screen window when you call OpenOffscreenWindow.\n"
    "\"color\" is the clut index (scalar or [r g b] triplet or [r g b a] quadruple) that you want to poke into "
    "each pixel as initial background color; default is white.\n"
	"\"rect\" specifies the size of the offscreen window "
    "If supplied, \"rect\" must contain at least one pixel. If a windowPtr is supplied, "
    "then \"rect\" defaults to the whole window. "
    "If a screenNumber is supplied then \"rect\" defaults to the whole screen. If a screenNumber of "
    "-1 is supplied, then \"rect\" defaults to the size of the main screen. (In all "
    "cases, subsequent references to this new offscreen window will use its "
    "coordinates: origin at its upper left.)\n"
	"\"pixelSize\" sets the depth (in bits) of "
    "each pixel. If you specify no screen (screenNumber=-1) then the "
    "default pixelSize is 32, but you can specify any legal depth: 8, 16, 24, 32. "
    "A pixelSize of 0 or [] is replaced by the default of 32 bits per pixel. If you specify "
    "a screen number of windowPtr, then the default depth is that of the screen or window. "
	"If you run your script with the imaging pipeline enabled (imagingmode flag > 0 in "
	"Screen('OpenWindow'), then all Offscreen windows always have 4 color channels RGBA and "
	"the selectable depths are additionally 64 bits or 128 bits, corresponding to 16 bits or "
	"32 bits floating point precision per color component. If 64 bits are selected but the "
	"hardware does not support this in float precision, a 15 bit precision per color channel "
	"signed integer format will be tried instead. On OpenGL-ES hardware, only the 32 bpc float "
    "type or 8 bit integer type is supported, therefore a pixelSize of more than 32 will always "
    "get silently upgraded to 128 bits per pixel, if possible.\n"
	"'specialFlags' optional parameter to set special properties, defaults to zero. "
    "If you set it to 1 then the offscreen window is created in GL_TEXTURE_2D format if possible. Use of "
    "GL_TEXTURE_2D format is currently not automatically compatible with use of specialFlags setting 2.\n"
    "If you set 'specialFlags' to 2 then the offscreen window will be drawn with especially high precision, see "
	"specialFlags setting of 2 in help for Screen('DrawTexture') for more explanation.\n"
    "A 'specialFlags' == 8 will prevent automatic mipmap-generation for GL_TEXTURE_2D textures.\n"
    "A 'specialFlags' == 32 will prevent automatic closing of the offscreen window by a call to Screen('Close');\n"
	"'multiSample' optional number of samples to use for anti-aliased drawing: This defaults "
	"to zero if omitted, ie., no anti-aliasing is performed when drawing into this offscreen "
	"window. If you set a positive non-zero number of samples and your system supports "
	"anti-aliased drawing to offscreen windows and the imaging pipeline is active, then "
	"Screen will try to allocate the offscreen window with at least the requested number "
	"of samples per pixel for anti-aliasing, but gracefully fall back to lower numbers if "
	"the hardware isn't capable of handling the requested number. Please note that anti-aliased "
	"offscreen windows can't be directly used for transformations in Screen('TransformTexture') "
	"or for drawing via Screen('DrawTexture') -- this is a hardware limitation. To display or "
	"use the content of an anti-aliased offscreen window you must first create a normal offscreen "
	"window of the same size and color format (same 'rect' and 'pixelSize' parameter), then use "
	"Screen('CopyWindow', antialiasedWindowhandle, normalWindowhandle); to copy the content to the "
	"normal offscreen window. This will perform the actual conversion into anti-aliased and "
	"displayable content.\n"
	"If the imaging pipeline is disabled, the 'multiSample' parameter will be silently ignored.\n\n"
    "NOTE: Screen's windows are known only to Screen and must be closed by it, e.g., "
    "Screen('Close', w). Matlab knows nothing about Screen's windows, so the Matlab "
    "CLOSE command won't work on Screen's windows. ";

static char seeAlsoString[] = "OpenWindow";

PsychError SCREENOpenOffscreenWindow(void) 
{
    int						screenNumber, depth, targetScreenNumber;
    PsychRectType			rect;
    PsychColorType			color;
    PsychWindowRecordType	*exampleWindowRecord, *windowRecord, *targetWindow;
    psych_bool				wasColorSupplied;
    char*					texturePointer;
    size_t					xSize, ySize, nbytes;
    psych_bool				bigendian;
	GLubyte					*rpb;
    int						ix;
	GLenum					fboInternalFormat;
	psych_bool				needzbuffer;
	psych_bool				overridedepth = FALSE;
	int						usefloatformat = 0;
	int						specialFlags = 0;
	int						multiSample = 0;
	
    // Detect endianity (byte-order) of machine:
    ix=255;
    rpb=(GLubyte*) &ix;
    bigendian = ( *rpb == 255 ) ? FALSE : TRUE;
    ix = 0; rpb = NULL;

    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //cap the number of inputs
    PsychErrorExit(PsychCapNumInputArgs(6));   //The maximum number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(2));  //The maximum number of outputs

    //1-User supplies a window ptr 2-User supplies a screen number 3-User supplies rect and pixel size
    if(PsychIsWindowIndexArg(1)){
        PsychAllocInWindowRecordArg(1, TRUE, &exampleWindowRecord);
		// Assign normalized copy of example windows rect -- Top-Left corner is always (0,0)
		PsychNormalizeRect(exampleWindowRecord->clientrect, rect);

        // We copy depth only from exampleWindow if it is a offscreen window (=texture). Copying from
        // onscreen windows doesn't make sense, e.g. depth=16 for onscreen means RGBA8 window, but it
        // would map onto a LUMINANCE+ALPHA texture for the offscreen window! We always use 32 bit RGBA8
        // in such a case.
        depth=(PsychIsOffscreenWindow(exampleWindowRecord)) ? exampleWindowRecord->depth : 32;
		// unless it is a FBO backed onscreen window in imaging mode: Then we can use the depth from it.
		if (exampleWindowRecord->imagingMode & kPsychNeedFastBackingStore || exampleWindowRecord->imagingMode & kPsychNeedFastOffscreenWindows) depth = exampleWindowRecord->depth;
        targetScreenNumber=exampleWindowRecord->screenNumber;
        targetWindow=exampleWindowRecord;
    } else if(PsychIsScreenNumberArg(1)){
        PsychCopyInScreenNumberArg(1, TRUE, &screenNumber);
        PsychGetScreenRect(screenNumber, rect);
        depth=32; // Always use RGBA8 in this case! See above...
        targetScreenNumber=screenNumber;
        targetWindow=NULL;
    } else if(PsychIsUnaffiliatedScreenNumberArg(1)){  //that means -1 or maybe also NaN if we add that option.  
        // Default to a depth of 32 bpp:
        depth = 32;
        targetWindow = NULL;
        // Get first open onscreen window as target window:
        PsychFindScreenWindowFromScreenNumber(kPsychUnaffiliatedWindow, &targetWindow);
		if (targetWindow == NULL) PsychErrorExitMsg(PsychError_user, "Could not find any open onscreen window to act as parent for this offscreen window. Open an onscreen window first!");
		targetScreenNumber = targetWindow->screenNumber;
        PsychGetScreenRect(targetScreenNumber, rect);
    } else {
        PsychErrorExit(PsychError_invalidNumdex);
    }

    if (targetWindow==NULL) {
        // Get target window of screen:
        PsychFindScreenWindowFromScreenNumber(targetScreenNumber, &targetWindow);
		if (targetWindow == NULL) PsychErrorExitMsg(PsychError_user, "Could not find any open onscreen window to act as parent for this offscreen window. Open an onscreen window first!");
		targetScreenNumber = targetWindow->screenNumber;
    }
    
    //Depth and rect argument supplied as arguments override those inherited from reference screen or window.
    //Note that PsychCopyIn* prefix means that value will not be overwritten if the arguments are not present.
    PsychCopyInRectArg(3,FALSE, rect);
    if (IsPsychRectEmpty(rect)) PsychErrorExitMsg(PsychError_user, "Invalid rect value provided: Empty rects are not allowed.");

	// Copy in optional depth: This gets overriden in many ways if imaging pipeline is on:
    if (PsychCopyInIntegerArg(4,FALSE, &depth)) overridedepth = TRUE;

    // If any of the no longer supported values 0, 1, 2 or 4 is provided, we
    // silently switch to 32 bits per pixel, which is the safest and fastest setting:
    if (depth==0 || depth==1 || depth==2 || depth==4) depth=32;

    // Final sanity check:
	if (!(targetWindow->imagingMode & kPsychNeedFastOffscreenWindows) && !(targetWindow->imagingMode & kPsychNeedFastBackingStore) && (depth==64 || depth==128)) {
      PsychErrorExitMsg(PsychError_user, "Invalid depth value provided. Must be 8 bpp, 16 bpp, 24 bpp or 32 bpp, unless you enable the imaging pipeline, which provides you with more options!");
	}
	
    if (depth!=8 && depth!=16 && depth!=24 && depth!=32 && depth!=64 && depth!=128) {
      PsychErrorExitMsg(PsychError_user, "Invalid depth value provided. Must be 8 bpp, 16 bpp, 24 bpp, 32 bpp, or if imagingmode is enabled also 64 bpp or 128 bpp!");
    }

	// If the imaging pipeline is enabled for the associated onscreen window and fast backing store, aka FBO's
	// is requested, then we only accept depths of at least 32 bit, i.e. RGBA windows. We override any lower
	// precision spec. This is because some common hardware only supports rendering to RGBA textures, not to
	// RGB, LA or Luminance textures.
	if ((targetWindow->imagingMode & kPsychNeedFastBackingStore || targetWindow->imagingMode & kPsychNeedFastOffscreenWindows) && (depth < 32)) depth = 32;

    // Find the color for the window background.  
    wasColorSupplied=PsychCopyInColorArg(kPsychUseDefaultArgPosition, FALSE, &color);
	
	// If none provided, use a proper white-value for this window:
    if(!wasColorSupplied) PsychLoadColorStruct(&color, kPsychIndexColor, PsychGetWhiteValueFromWindow(targetWindow));  

    // Get the optional specialmode flag:
    PsychCopyInIntegerArg(5, FALSE, &specialFlags);

    // OpenGL-ES only supports GL_TEXTURE_2D targets, so enforce these via flags setting 1:
    if (PsychIsGLES(targetWindow)) specialFlags |= 1;

	// This command converts whatever color we got into RGBA format:
    PsychCoerceColorMode(&color);

    // printf("R=%i G=%i B=%i A=%i I=%i", color.value.rgba.r, color.value.rgba.g,color.value.rgba.b,color.value.rgba.a,color.value.index); 
    // First allocate the offscreen window record to store stuff into. If we exit with an error PsychErrorExit() should
    // call PsychPurgeInvalidWindows which will clean up the window record. 
    PsychCreateWindowRecord(&windowRecord);  // This also fills the window index field.

    // This offscreen window is implemented as a Psychtoolbox texture:
    windowRecord->windowType=kPsychTexture;

    // We need to assign the screen number of the onscreen-window, so PsychCreateTexture()
    // can query the size of the screen/onscreen-window...
    windowRecord->screenNumber = targetScreenNumber;

    // Assign the computed depth:
    windowRecord->depth=depth;

	// Default number of channels:
	windowRecord->nrchannels=depth / 8;

    // Assign the computed rect, but normalize it to start with top-left at (0,0):
    PsychNormalizeRect(rect, windowRecord->rect);

    // Client rect of an offscreen window is always == rect of it:
    PsychCopyRect(windowRecord->clientrect, windowRecord->rect);
    
	// Until here no OpenGL commands executed. Now we need a valid context: Set targetWindow
	// as drawing target. This will perform neccessary context-switch and all backbuffer
	// backup/restore/whatever operations to make sure we can do what we want without
	// possibly screwing any offscreen windows and bindings:
	if (PsychIsOnscreenWindow(targetWindow) || PsychIsOffscreenWindow(targetWindow)) {
		// This is a possible on-/offscreen drawingtarget:
		PsychSetDrawingTarget(targetWindow);
	}
	else {
		// This must be a proxy-window object: Can't transition to it!
		
		// But we can safe-reset the current drawingtarget...
		PsychSetDrawingTarget((PsychWindowRecordType*) 0x1);
		
		// ...and then switch to the OpenGL context of the 'targetWindow' proxy object:
		PsychSetGLContext(targetWindow);

		// Ok, framebuffer and bindings are safe and disabled, context is set. We
		// should be safe to continue with the proxy...
	}
	
	// From here on we have a defined context and state. We can detach the drawing target whenever
	// we want, as everything is backed up somewhere for later reinit.
	
	// Create offscreen window either new style as FBO, or old style as texture:
	if ((targetWindow->imagingMode & kPsychNeedFastBackingStore) || (targetWindow->imagingMode & kPsychNeedFastOffscreenWindows)) {
		// Imaging mode for this window enabled: Use new way of creating the offscreen window:
		
		// We safely unbind any FBO bindings and drawingtargets:
		PsychSetDrawingTarget((PsychWindowRecordType*) 0x1);
		
		// Overriden for imagingmode: There we always have 4 channels...
		windowRecord->nrchannels=4;

		// Start off with standard 8 bpc fixed point:
		fboInternalFormat = GL_RGBA8; windowRecord->depth=32; usefloatformat = 0;
		
		// Need 16 bpc fixed point precision?
		if (targetWindow->imagingMode & kPsychNeed16BPCFixed) {
			fboInternalFormat = (targetWindow->gfxcaps & kPsychGfxCapSNTex16) ? GL_RGBA16_SNORM : GL_RGBA16;
			windowRecord->depth=64;
			usefloatformat = 0;
		}
		
		// Need 16 bpc floating point precision?
		if (targetWindow->imagingMode & kPsychNeed16BPCFloat) { fboInternalFormat = GL_RGBA_FLOAT16_APPLE; windowRecord->depth=64; usefloatformat = 1; }
		
		// Need 32 bpc floating point precision?
		if (targetWindow->imagingMode & kPsychNeed32BPCFloat) { fboInternalFormat = GL_RGBA_FLOAT32_APPLE; windowRecord->depth=128; usefloatformat = 2; }
		
		// Override depth value provided?
		if (overridedepth) {
			// Manual depth specified: Override with that depth:
			switch(depth) {
				case 32:
					fboInternalFormat = GL_RGBA8; windowRecord->depth=32; usefloatformat = 0;
				break;

				case 64:
					fboInternalFormat = GL_RGBA_FLOAT16_APPLE; windowRecord->depth=64; usefloatformat = 1;
					// Need fallback for lack of float 16 support?
					if (!(targetWindow->gfxcaps & kPsychGfxCapFPTex16) && !PsychIsGLES(targetWindow)) {
						// Yes. Try 16 bit signed normalized texture instead:
						if (PsychPrefStateGet_Verbosity() > 4)
							printf("PTB-INFO:OpenOffscreenWindow: Code requested 16 bpc float precision, but this is unsupported. Trying to use 15 bit snorm precision instead.\n");
						fboInternalFormat = GL_RGBA16_SNORM; windowRecord->depth=64; usefloatformat = 0;
						if (!(targetWindow->gfxcaps & kPsychGfxCapSNTex16)) {
							printf("PTB-ERROR:OpenOffscreenWindow: Code requested 16 bpc float precision, but this is unsupported by this graphics card.\n");
							printf("PTB-ERROR:OpenOffscreenWindow: Tried to use 16 bit snorm format instead, but failed as this is unsupported as well.\n");
						}
					}
				break;

				case 128:
					fboInternalFormat = GL_RGBA_FLOAT32_APPLE; windowRecord->depth=128; usefloatformat = 2;
				break;
				
				default:
					fboInternalFormat = GL_RGBA8; windowRecord->depth=32; usefloatformat = 0;
			}			
		}
		
        // Floating point framebuffer on OpenGL-ES requested?
        if (PsychIsGLES(targetWindow) && (usefloatformat > 0)) {
            // Yes. We only support 32 bpc float framebuffers with alpha-blending. On less supportive hardware we fail:
            if (!(targetWindow->gfxcaps & kPsychGfxCapFPTex32) || !(targetWindow->gfxcaps & kPsychGfxCapFPFBO32)) {
                PsychErrorExitMsg(PsychError_user, "Sorry, the requested offscreen window color resolution of 32 bpc floating point is not supported by your graphics card. Game over.");
            }

            // Supported. Upgrade requested format to 32 bpc float, whatever it was before:
            fboInternalFormat = GL_RGBA_FLOAT32_APPLE; windowRecord->depth=128; usefloatformat = 2;
        }

		// Do we need additional depth buffer attachments?
		needzbuffer = (PsychPrefStateGet_3DGfx()>0) ? TRUE : FALSE;
		
		// Copy in optional multiSample argument: It defaults to zero, aka multisampling disabled.
		PsychCopyInIntegerArg(6, FALSE, &multiSample);
		if (multiSample < 0) PsychErrorExitMsg(PsychError_user, "Invalid negative multiSample level provided!");

		// Multisampled anti-aliasing requested?
		if (multiSample > 0) {
			// Yep. Supported by GPU?
			if (!(targetWindow->gfxcaps & kPsychGfxCapFBOMultisample)) {
				// No. We fall back to non-multisampled mode:
				multiSample = 0;
				
				// Tell user if warnings enabled:
				if (PsychPrefStateGet_Verbosity() > 1) {
					printf("PTB-WARNING: You requested stimulus anti-aliasing via multisampling by setting the multiSample parameter of Screen('OpenOffscreenWindow', ...) to a non-zero value.\n");
					printf("PTB-WARNING: You also requested use of the imaging pipeline. Unfortunately, your combination of operating system, graphics hardware and driver does not\n");
					printf("PTB-WARNING: support simultaneous use of the imaging pipeline and multisampled anti-aliasing.\n");
					printf("PTB-WARNING: Will therefore continue without anti-aliasing...\n\n");
					printf("PTB-WARNING: A driver upgrade may resolve this issue. Users of MacOS-X need at least OS/X 10.5.2 Leopard for support on recent ATI hardware.\n\n");
				}
			}
		}

		// Allocate framebuffer object for this Offscreen window:
		if (!PsychCreateFBO(&(windowRecord->fboTable[0]), fboInternalFormat, needzbuffer, (int) PsychGetWidthFromRect(rect), (int) PsychGetHeightFromRect(rect), multiSample, specialFlags)) {
			// Failed!
			PsychErrorExitMsg(PsychError_user, "Creation of Offscreen window in imagingmode failed for some reason :(");
		}

		// Assign this FBO as drawBuffer for mono channel of our Offscreen window:
		windowRecord->drawBufferFBO[0] = 0;
		windowRecord->fboCount = 1;
		
		// Assign it as texture as well:
		windowRecord->textureNumber = windowRecord->fboTable[0]->coltexid;
		windowRecord->textureMemorySizeBytes = 0;
		windowRecord->textureMemory = NULL;
		windowRecord->texturetarget = (specialFlags & 0x1) ? GL_TEXTURE_2D : GL_TEXTURE_RECTANGLE_EXT;
		windowRecord->surfaceSizeBytes = (size_t) (PsychGetWidthFromRect(rect) * PsychGetHeightFromRect(rect) * (windowRecord->depth / 8));

		// Set bpc for FBO backed offscreen window:
		windowRecord->bpc = (int) (windowRecord->depth / 4);

		// Initial setup done, continues below after some shared code...
	}
	else {
		// Traditional texture creation code:
		
		// Special case for alpha-channel: DBL_MAX signals maximum alpha
		// value requested. In our own code we need to manually map this to
		// the maximum uint8 alpha value of 255:
		if (color.value.rgba.a == DBL_MAX) color.value.rgba.a = 255;
		
		// Allocate the texture memory:
		// We only allocate the amount really needed for given format, aka numMatrixPlanes - Bytes per pixel.
		xSize = (size_t) PsychGetWidthFromRect(rect);
		ySize = (size_t) PsychGetHeightFromRect(rect);
		windowRecord->textureMemorySizeBytes = ((size_t) (depth/8)) * xSize * ySize;
		windowRecord->textureMemory = malloc(windowRecord->textureMemorySizeBytes);
		texturePointer=(char*) windowRecord->textureMemory;
		// printf("depth=%i xsize=%i ysize=%i mem=%i ptr=%p", depth, xSize, ySize, windowRecord->textureMemorySizeBytes, texturePointer);
		// Fill with requested background color:
		nbytes=0;
		switch (depth) {
			case 8: // Pure LUMINANCE texture:
				memset((void*) texturePointer, (int) color.value.rgba.r, windowRecord->textureMemorySizeBytes);
				break;
				
			case 16: // LUMINANCE + ALPHA
				while (nbytes < windowRecord->textureMemorySizeBytes) {
					*(texturePointer++) = (psych_uint8) color.value.rgba.r;
					*(texturePointer++) = (psych_uint8) color.value.rgba.a;
					nbytes+=2;
				}
				break;
				
			case 24: // RGB:
				while (nbytes < windowRecord->textureMemorySizeBytes) {
					*(texturePointer++) = (psych_uint8) color.value.rgba.r;
					*(texturePointer++) = (psych_uint8) color.value.rgba.g;
					*(texturePointer++) = (psych_uint8) color.value.rgba.b;
					nbytes+=3;
				}
				break;        
				
			case 32: // RGBA
				if (bigendian) {
					// Code for big-endian machines, e.g., PowerPC:
					while (nbytes < windowRecord->textureMemorySizeBytes) {
						*(texturePointer++) = (psych_uint8) color.value.rgba.a;
						*(texturePointer++) = (psych_uint8) color.value.rgba.r;
						*(texturePointer++) = (psych_uint8) color.value.rgba.g;
						*(texturePointer++) = (psych_uint8) color.value.rgba.b;
						nbytes+=4;
					}
				}
				else {
					// Code for little-endian machines, e.g., IntelPC, IntelMAC, aka Pentium.
					while (nbytes < windowRecord->textureMemorySizeBytes) {
						*(texturePointer++) = (psych_uint8) color.value.rgba.b;
						*(texturePointer++) = (psych_uint8) color.value.rgba.g;
						*(texturePointer++) = (psych_uint8) color.value.rgba.r;
						*(texturePointer++) = (psych_uint8) color.value.rgba.a;
						nbytes+=4;
					}
				}
				break;
		}
	}
	
	// Shared setup code for FBO vs. non-FBO Offscreen windows:
	
	// Assign parent window and copy its inheritable properties:
	PsychAssignParentWindow(windowRecord, targetWindow);
	
    // Texture orientation is type 2 aka upright, non-transposed aka Offscreen window:
    windowRecord->textureOrientation = 2;
    
	if ((windowRecord->imagingMode & kPsychNeedFastBackingStore) || (windowRecord->imagingMode & kPsychNeedFastOffscreenWindows)) {
		// Last step for FBO backed Offscreen window: Clear it to its background color:
		PsychSetDrawingTarget(windowRecord);

		// Set default draw shader:
		PsychSetShader(windowRecord, -1);
	
		// Set background fill color:
		PsychSetGLColor(&color, windowRecord);

		// Setup alpha-blending:
		PsychUpdateAlphaBlendingFactorLazily(windowRecord);

		// Fullscreen fill of a non-onscreen window:
		PsychGLRect(windowRecord->rect);

		// Multisampling requested? If so, we need to enable it:
		if (multiSample > 0) {
			glEnable(GL_MULTISAMPLE);
			while (glGetError() != GL_NO_ERROR);
		}
		
		// Ready. Unbind it.
		PsychSetDrawingTarget(NULL);		
	}
	else {
		// Old-style setup for non-FBO Offscreen windows:
        
        // Special texture format?
		if (specialFlags & 0x1) windowRecord->texturetarget = GL_TEXTURE_2D;
        
		// Let's create and bind a new texture object and fill it with our new texture data.
		PsychCreateTexture(windowRecord);
    }

	// Assign GLSL filter-/lookup-shaders if needed:
	PsychAssignHighPrecisionTextureShaders(windowRecord, targetWindow, usefloatformat, (specialFlags & 2) ? 1 : 0);
	
    // specialFlags setting 8? Disable auto-mipmap generation:
    if (specialFlags & 0x8) windowRecord->specialflags |= kPsychDontAutoGenMipMaps;    

    // A specialFlags setting of 32? Protect texture against deletion via Screen('Close') without providing a explicit handle:
    if (specialFlags & 32) windowRecord->specialflags |= kPsychDontDeleteOnClose;    

    // Window ready. Mark it valid and return handle to userspace:
    PsychSetWindowRecordValid(windowRecord);
    
    //Return the window index and the rect argument.
    PsychCopyOutDoubleArg(1, FALSE, windowRecord->windowIndex);
    PsychCopyOutRectArg(2, FALSE, rect);

    // Ready.
    return(PsychError_none);
}
