/*
	SCREENFrameRect.c		

	AUTHORS:

	Allen.Ingling@nyu.edu			awi 
	mario.kleiner@tuebingen.mpg.de	mk

	PLATFORMS:	
	
	All.

	HISTORY:

		07/23/04	awi		Created.
		10/12/04	awi		In useString: moved commas to inside [].
		2/25/05		awi		Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
		  

	TO DO:

*/


#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('FrameRect', windowPtr [,color] [,rect] [,penWidth]);";
//                                             1           2        3       4          
static char synopsisString[] = 
	"Draw the outline of a rectangle \"rect\"."
	" \"color\" is the clut index (scalar or [r g b] triplet or [r g b a] quadruple) "
	"that you want to poke into each pixel;  default produces white with the standard CLUT for this "
	"window's pixelSize. Default \"rect\" is entire window.\n"
	"Instead of framing one rectangle, you can also specify a list of multiple rectangles to be "
	"framed - this is much faster when you need to draw many rectangles per frame. To frame n "
	"rectangles, provide \"rect\" as a 4 rows by n columns matrix, each column specifying one "
	"rectangle, e.g., rect(1,5)=left border of 5th rectange, rect(2,5)=top border of 5th rectangle, "
	"rect(3,5)=right border of 5th rectangle, rect(4,5)=bottom border of 5th rectangle. If the "
	"rectangles should have different colors, then provide \"color\" as a 3 or 4 row by n column "
	"matrix, the i'th column specifiying the color of the i'th rectangle. You can also specify a "
    "single pen size or a vector pen sizes (one for each rectangle) to control the thickness of "
	"the lines used to frame the rect. Default \"penWidth\" is 1 pixel unit. ";
	
static char seeAlsoString[] = "FillRect";	

PsychError SCREENFrameRect(void)  
{	
	PsychRectType					rect;
	PsychWindowRecordType			*windowRecord;
	psych_bool						isArgThere;
	double							penSize, lf, fudge;
    double							*xy, *colors, *penSizes;
	unsigned char					*bytecolors;
	int								numRects, i, j, nc, mc, nrsize;

	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(4));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs

	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(1, kPsychArgRequired, &windowRecord);

	// Get tweakable correction factor for framerect: Not an option on non-classic GL anymore.
	lf = (PsychIsGLClassic(windowRecord)) ? PsychPrefStateGet_FrameRectCorrection() : -1;
	
	// Query, allocate and copy in all vectors...
	numRects = 4;
	nrsize = 0;
	colors = NULL;
	bytecolors = NULL;
	penSizes = NULL;
	
	// The negative position -3 means: xy coords are expected at position 3, but they are optional.
	// NULL means - don't want a size's vector.
	PsychPrepareRenderBatch(windowRecord, -3, &numRects, &xy, 2, &nc, &mc, &colors, &bytecolors, 4, &nrsize, &penSizes, FALSE);

	// Default rect is fullscreen:
	PsychCopyRect(rect, windowRecord->clientrect);

	// Only up to one rect provided?
	if (numRects <= 1) {
		// Get the rect and draw it
		isArgThere=PsychCopyInRectArg(kPsychUseDefaultArgPosition, FALSE, rect);	
		if (isArgThere && IsPsychRectEmpty(rect)) return(PsychError_none);
		numRects = 1;
	}

	// Pen size starts as "undefined", just to make sure it gets initially set:
	penSize = -DBL_MAX;
	
	// Framed rect drawing loop:
	for (i=0; i<numRects; i++) {
		// Multiple rects to draw or single iteration to draw provided rect?
		if (numRects > 1) {
			// Multi-Rect drawing: Assign next rect from array and setup corresponding
			// color and penSize...

			// Assign rect:
			rect[kPsychLeft] = xy[i*4 + 0];
			rect[kPsychTop] = xy[i*4 + 1];
			rect[kPsychRight] = xy[i*4 + 2];
			rect[kPsychBottom] = xy[i*4 + 3];

			// Per rect color provided?
			if (nc>1) {
				// Yes. Set color for this specific rect:
				PsychSetArrayColor(windowRecord, i, mc, colors, bytecolors);
			}
		}
		else {
			// Only one single rect to draw in this single loop iteration.
			// The rect is already set up in 'rect', and the drawing color has
			// been set as well by PsychPrepareRenderBatch(). penSize has been
			// set by that routine as well in penSizes[0], so we don't have
			// anything to do here...
			// NO OP.
		}
		
		j = (nrsize > 1) ? i : 0; 
		
        if (penSizes[j] != penSize) {
			penSize = penSizes[j];
			if (lf != -1) glLineWidth((GLfloat) penSize);
		}
		
		if (IsPsychRectEmpty(rect)) continue;

		if (lf == -1) {
			// New style rendering: More robust against variations in GPU implementations:
			fudge = penSize;
            GLRECTd(rect[kPsychLeft], rect[kPsychTop], rect[kPsychRight], rect[kPsychTop] + fudge);
			GLRECTd(rect[kPsychLeft], rect[kPsychBottom], rect[kPsychRight], rect[kPsychBottom] - fudge);
			GLRECTd(rect[kPsychLeft], rect[kPsychTop]+fudge, rect[kPsychLeft]+fudge, rect[kPsychBottom]-fudge);
			GLRECTd(rect[kPsychRight]-fudge, rect[kPsychTop]+fudge, rect[kPsychRight], rect[kPsychBottom]-fudge);
		}
		else {
			// Old style: Has a couple of problems in corner cases. Left for now as reference...
			if (penSize > 1) {
				// Width > 1
				
				fudge = (penSize > 1) ? lf * penSize/2 : 0.0;
				
				glBegin(GL_LINES);
				// Draw 4 separate segments, extend the left and right
				// vertical segments by half a penWidth.
				glVertex2d(rect[kPsychLeft], rect[kPsychTop] - fudge);
				glVertex2d(rect[kPsychLeft], rect[kPsychBottom] + fudge);
				glVertex2d(rect[kPsychRight], rect[kPsychTop]);
				glVertex2d(rect[kPsychLeft], rect[kPsychTop]);
				glVertex2d(rect[kPsychRight], rect[kPsychBottom] + fudge);
				glVertex2d(rect[kPsychRight], rect[kPsychTop] - fudge);
				glVertex2d(rect[kPsychRight], rect[kPsychBottom]);
				glVertex2d(rect[kPsychLeft], rect[kPsychBottom]);
				glEnd();			
			}
			else {
				// Width <= 1: Simple case...
				glBegin(GL_LINE_LOOP);
				glVertex2d(rect[kPsychLeft], rect[kPsychBottom]);
				glVertex2d(rect[kPsychLeft], rect[kPsychTop]);
				glVertex2d(rect[kPsychRight], rect[kPsychTop]);
				glVertex2d(rect[kPsychRight], rect[kPsychBottom]);
				glEnd();
			}
		}
		// Next rect...
	}
	
	// Need to reset line width?
	if (penSize!=1 && lf!=-1) glLineWidth(1);

	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);

	return(PsychError_none);
}
