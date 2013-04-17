/*
    SCREENFrameOval.c		
  
    AUTHORS:
    
		Allen.Ingling@nyu.edu				awi
		mario.kleiner@tuebingen.mpg.de		mk
  
    PLATFORMS:
	
		All.

    HISTORY:
	
		mm/dd/yy	
		
		10/10/03	awi		Created.  Based on FillOval.
		10/12/04	awi		In useString: changed "SCREEN" to "Screen", and moved commas to inside [].
		1/15/05		awi		Removed GL_BLEND setting a MK's suggestion. 
		1/25/05		awi		Really removed GL_BLEND.  Correction provide by mk.
		2/25/05		awi		Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
		6/14/09      mk		Add batch-drawing support, just as with FillOval et al.

    TO DO:
    
    BUGS:  
    
		The pen width is not uniform along the oval if it is elongated.  This is because we scale the line thickness with the oval.
		Here are some ideas about how to fix this:
			-The OS9 Psychtoolbox used drawing commands which mapped directly onto QuickDraw commands.  The same drawing commands are a bad fit for OpenGL.
			After version 1.0 we should migrate the Psychtoolbox to a new set of drawing commands which map closely onto GL commands.  Meanwhile, for backwards
			compatability with older scripts, we might consider implementing the old drawing commands with actual QuickDraw or perhaps Quartz calls.  See 
			Apple's QA1011 for how to combine QuickDraw with CGDirectDisplay:
			http://developer.apple.com/qa/qa2001/qa1011.html

			-We could composite an oval by drawing a smaller background oval within a larger forground oval and setting alpha on the smaller oval to transparent and
			on the larger oval to opaque. The problem with this is that we actally have to allocate a compositing window since, unless there is some clever way to 
			do this in the destination window without disturbing he existing image, by careful selectin of copy modes.
			
			-We could implement (or borrow) an oval drawing routine.  This seems like a bad solution.
			
		MK: All these proposals are not workable solutions, because they would either be awfully slow,
			or have significant side-effects on things like HDR drawing or alpha blending.
			
			On modern GPU's a shader based solution would be perfect and fast.

*/

#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('FrameOval', windowPtr [,color] [,rect] [,penWidth] [,penHeight] [,penMode]);";
//                                             1           2        3      4            5            6
static char synopsisString[] = 
            "Draw the outline of an oval inscribed in \"rect\". \"color\" is the clut index (scalar or [r g b] "
            "triplet) that you want to poke into each pixel; default produces white with the "
            "standard CLUT for this window's pixelSize. Default \"rect\" is entire window.\n"
            "Default pen size is (1,1). The penMode argument is ignored. The penWidth must "
            "equal the penHeight.  If non-equal arguments are given, FrameOval will choose the maximum "
            "value of both. The pen width will be non-uniform for non-circular ovals, this is a known "
			"limitation.\n"
			"Instead of drawing one oval, you can also specify a list of multiple ovals to be "
			"drawn - this is faster when you need to draw many ovals per frame. To draw n "
			"ovals, provide \"rect\" as a 4 rows by n columns matrix, each column specifying one "
			"oval, e.g., rect(1,5)=left border of 5th oval, rect(2,5)=top border of 5th oval, "
			"rect(3,5)=right border of 5th oval, rect(4,5)=bottom border of 5th oval. If the "
			"ovals should have different colors, then provide \"color\" as a 3 or 4 row by n column "
			"matrix, the i'th column specifiying the color of the i'th oval.\n"
			"If drawing multiple ovals at once, both the penHeight and penMode arguments are completely "
			"ignored! The penWidth argument defines the size of each drawn oval. You can either specify "
			"one common penWidth for all ovals, or provide a per-oval penWidth.\n";

static char seeAlsoString[] = "FillOval";	
            
PsychError SCREENFrameOval(void)  
{
	PsychRectType			rect;
	double					numSlices, outerRadius, xScale, yScale, xTranslate, yTranslate, rectY, rectX, penWidth, penHeight, penSize, innerRadius;
	PsychWindowRecordType	*windowRecord;
	psych_bool				isArgThere, isclassic;
    double					*xy, *colors;
	unsigned char			*bytecolors;
	double*					penSizes;
	int						numRects, i, nc, mc, nrsize;
	GLUquadricObj			*diskQuadric;

	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);}
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(6));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs

	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);

	// Query, allocate and copy in all vectors...
	numRects = 4;
	nrsize = 0;
	colors = NULL;
	bytecolors = NULL;
	mc = nc = 0;
	
	// The negative position -3 means: xy coords are expected at position 3, but they are optional.
	// NULL means - don't want a size's vector.
	PsychPrepareRenderBatch(windowRecord, -3, &numRects, &xy, 2, &nc, &mc, &colors, &bytecolors, 4, &nrsize, &penSizes, FALSE);
    isclassic = PsychIsGLClassic(windowRecord);

	// Only up to one rect provided?
	if (numRects <= 1) {
		// Get the oval and draw it:
		PsychCopyRect(rect, windowRecord->clientrect);
		isArgThere=PsychCopyInRectArg(kPsychUseDefaultArgPosition, FALSE, rect);	
		if (isArgThere && IsPsychRectEmpty(rect)) return(PsychError_none);
		numRects = 1;

		// Get the pen width and height arguments
		penWidth=1;
		penHeight=1;
		PsychCopyInDoubleArg(4, FALSE, &penWidth);
		PsychCopyInDoubleArg(5, FALSE, &penHeight);
		penSize = (penWidth > penHeight) ? penWidth : penHeight;
	}
	else {
		// Multiple ovals provided. Set up the first one:
        PsychCopyRect(rect, &xy[0]);
		penSize = penSizes[0];
	}

	// Create quadric object:
	if (isclassic) diskQuadric = gluNewQuadric();

	// Draw all ovals (one or multiple):
	for (i=0; i < numRects;) {
		// Per oval color provided? If so then set it up. If only one common color
		// was provided then PsychPrepareRenderBatch() has already set it up.
		if (nc>1) {
			// Yes. Set color for this specific item:
			PsychSetArrayColor(windowRecord, i, mc, colors, bytecolors);
		}

		// Per oval penSize provided? If so, set it up. Otherwise keep at default size
		// common for all ovals, set by code outside loop:
		if (nrsize > 1) penSize = penSizes[i];

		// Compute drawing parameters for ellipse:
		if (!IsPsychRectEmpty(rect)) {
			//The glu disk object location and size with a  center point and a radius,   
			//whereas FrameOval accepts a bounding rect. Converting from one set of parameters
			//to the other we should careful what we do for rects size of even number of pixels in length.
			PsychGetCenterFromRectAbsolute(rect, &xTranslate, &yTranslate);
			rectY=PsychGetHeightFromRect(rect);
			rectX=PsychGetWidthFromRect(rect);
			if(rectX == rectY){
				xScale=1; 
				yScale=1;
				outerRadius=rectX/2;
			}else if(rectX > rectY){ 
				xScale=1;
				yScale=rectY/rectX;
				outerRadius=rectX/2;
			}else if(rectY > rectX){
				yScale=1;
				xScale=rectX/rectY;
				outerRadius=rectY/2;
			}
			
			numSlices   =   3.14159265358979323846  * 2 * outerRadius;
			innerRadius = outerRadius - penSize;
			innerRadius = (innerRadius < 0) ? 0 : innerRadius;         

            if (isclassic) {
                // Draw: Set up position, scale and size via matrix transform:
                glPushMatrix();
                glTranslated(xTranslate, yTranslate, 0);
                glScaled(xScale, yScale, 1);

                // Compute disk quadric for given params: This is awfully slow and would
                // benefit a lot from shader magic on modern GPUs:
                gluDisk(diskQuadric, innerRadius, outerRadius, (int) numSlices, 1);
                glPopMatrix();
            }
            else {
                PsychDrawDisc(windowRecord, (float) xTranslate, (float) yTranslate, (float) innerRadius, (float) outerRadius, (int) numSlices, (float) xScale, (float) yScale, 0, 360);
            }
		}
		
		// Done with this one. Set up the next one, if any...
		i++;
		if (i < numRects) {
            PsychCopyRect(rect, &xy[i*4]);
        }

		// Next oval.
	}

	// Release quadric object:
	if (isclassic) gluDeleteQuadric(diskQuadric);

	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);
	
 	// All Psychfunctions require this.
	return(PsychError_none);
}
