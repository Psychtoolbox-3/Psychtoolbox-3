/*
	SCREENFillRect.c		
  
	AUTHORS:

		Allen.Ingling@nyu.edu           awi 
        mario.kleiner@tuebingen.mpg.de  mk
 
	PLATFORMS:	
	
		All.
    

	HISTORY:

		09/12/02	awi		Created.  
		10/12/04	awi		In useString: changed "SCREEN" to "Screen", and moved commas to inside [].
		2/25/05		awi		Relocated PsychSetGLContext() to outside condtional, it only executed for small rects.
							glClearColor() now sets variable alpha, not static at 1.0 (255). 
							Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
 
 
	TO DO:
  
		FillRect should detect when it is called without a rect argument and use glClearColor to clear the 
		window to the specified color because  gl implementations may optinally accellerate that. 
  
		Break out the actual rendering commands into a seperate function because OpenWindow and OpenOffscreenWindow
		need to call FillRect to initialize the window to the correct color.  
  
*/


#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('FillRect', windowPtr [,color] [,rect] )";                                          
//                                            1           2        3      
static char synopsisString[] = 
	"Fill \"rect\". \"color\" is the clut index (scalar or [r g b] triplet or [r g b a] quadruple) "
	"that you want to poke into each pixel;  default produces white with the standard CLUT for this "
	"window's pixelSize. Default \"rect\" is entire window, so you can use this function to clear "
	"the window. Please note that clearing the entire window will set the background color of the "
	"window to the clear color, ie., future Screen('Flip') commands will clear to the new background "
	"clear color specified in Screen('FillRect').\n"
	"Instead of filling one rectangle, you can also specify a list of multiple rectangles to be "
	"filled - this is much faster when you need to draw many rectangles per frame. To fill n "
	"rectangles, provide \"rect\" as a 4 rows by n columns matrix, each column specifying one "
	"rectangle, e.g., rect(1,5)=left border of 5th rectange, rect(2,5)=top border of 5th rectangle, "
	"rect(3,5)=right border of 5th rectangle, rect(4,5)=bottom border of 5th rectangle. If the "
	"rectangles should have different colors, then provide \"color\" as a 3 or 4 row by n column "
	"matrix, the i'th column specifiying the color of the i'th rectangle. ";
static char seeAlsoString[] = "FrameRect";	

PsychError SCREENFillRect(void)  
{
	
	PsychColorType					color;
	PsychRectType					rect;
	PsychWindowRecordType			*windowRecord;
	double							whiteValue;
	psych_bool						isArgThere, isScreenRect;
    double							*xy, *colors;
	unsigned char					*bytecolors;
	int								numRects, i, nc, mc, nrsize;

	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(4));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs

	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(1, kPsychArgRequired, &windowRecord);
	
	// Query, allocate and copy in all vectors...
	numRects = 4;
	nrsize = 0;
	colors = NULL;
	bytecolors = NULL;

	// The negative position -3 means: xy coords are expected at position 3, but they are optional.
	// NULL means - don't want a size's vector.
	PsychPrepareRenderBatch(windowRecord, -3, &numRects, &xy, 2, &nc, &mc, &colors, &bytecolors, 0, &nrsize, NULL, FALSE);
	isScreenRect=FALSE;
	
	// Only up to one rect provided?
	if (numRects <= 1) {
		// Get the rect and draw it
		isArgThere=PsychCopyInRectArg(kPsychUseDefaultArgPosition, FALSE, rect);	
		isScreenRect= !isArgThere || isArgThere && PsychMatchRect(rect, windowRecord->clientrect);
		if (isArgThere && IsPsychRectEmpty(rect)) return(PsychError_none);
	}

	// There are tons of cases where we can't use glClear() to do a fullscreen fill. E.g., it
	// doesn't work on anything else than onscreen windows, doesn't work in many stereo modes
	// unless imaging pipe is fully active, and doesn't work when color clamping is disabled
	// by use of our own GLSL shader based workaround. If in doubt, we use a fullscreen rectangle
	// blit to do the clear -- Slower and disables tons of gfx-driver optimizations like fast clears
	// and framebuffer compression / Hyper-z stuff, but at least safe and well-defined:
	if(isScreenRect && PsychIsOnscreenWindow(windowRecord) && !((windowRecord->defaultDrawShader != 0) && (windowRecord->defaultDrawShader == windowRecord->unclampedDrawShader)) &&
	   (windowRecord->stereomode < kPsychAnaglyphRGStereo || windowRecord->stereomode > kPsychAnaglyphBRStereo || (windowRecord->imagingMode > 0 && windowRecord->imagingMode != kPsychNeedFastOffscreenWindows))){
		// Fullscreen rect fill which in GL is a special case which may be accelerated.
		// We only use this fast-path on real onscreen windows, not on textures or
		// offscreen windows.
		
		//Get the color argument or use the default, then coerce to the form determened by the window depth.
		isArgThere=PsychCopyInColorArg(2, FALSE, &color);
		if(!isArgThere){
			whiteValue=PsychGetWhiteValueFromWindow(windowRecord);
			PsychLoadColorStruct(&color, kPsychIndexColor, whiteValue ); //index mode will coerce to any other.
		}
		PsychCoerceColorMode( &color);

		// Disable any active shaders:
		PsychSetShader(windowRecord, 0);

		PsychConvertColorToDoubleVector(&color, windowRecord, windowRecord->clearColor);
		PsychGLClear(windowRecord);

		// Fixup possible low-level framebuffer layout changes caused by commands above this point. Needed from native 10bpc FB support to work reliably.
		PsychFixupNative10BitFramebufferEnableAfterEndOfSceneMarker(windowRecord);

		// Done with fast clear path.
	}else{
	  // Subregion fill or fullscreen fill into offscreen window or texture: Draw a colored rect.

	  // At this point if we have a single color value, it is already set by PsychPrepareRenderBatch,
	  // so we can skip this: PsychSetGLColor(&color, windowRecord);
	  // that function has also set the proper shader, if any...
	  
	  // Fullscreen or partial fill?
	  if (isScreenRect) {
	    // Fullscreen fill of a (non-)onscreen window:
		
		// Draw a rect in the clear color:
	    PsychGLRect(windowRecord->clientrect);
		
		// If this was an onscreen window which couldn't be glClear()'ed for some other reason,
		// we need to set the glClearColor() in the windowRecord:
		if (PsychIsOnscreenWindow(windowRecord)) {
			//Get the color argument or use the default, then coerce to the form determened by the window depth.  
			isArgThere=PsychCopyInColorArg(2, FALSE, &color);
			if(!isArgThere){
				whiteValue=PsychGetWhiteValueFromWindow(windowRecord);
				PsychLoadColorStruct(&color, kPsychIndexColor, whiteValue ); //index mode will coerce to any other.
			}
			PsychCoerceColorMode( &color);
			PsychConvertColorToDoubleVector(&color, windowRecord, windowRecord->clearColor);
		}
	  } else {
	    // Partial fill: Draw provided rects:
		if (numRects>1) {
			// Multiple rects provided: Draw the whole batch:
			for (i=0; i<numRects; i++) {
				// Per rect color provided?
				if (nc>1) {
					// Yes. Set color for this specific rect:
					PsychSetArrayColor(windowRecord, i, mc, colors, bytecolors);
				}
				
				// Submit rect for drawing:
				if (!IsPsychRectEmpty(rect)) PsychGLRect(&(xy[i*4]));
			}
		}
		else {
			// Single partial screen rect provided: Draw it.
			PsychGLRect(rect);
		}
	  }
	}
	
	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);
	
 	//All psychfunctions require this.
	return(PsychError_none);
}
