/*
	SCREENFillOval.c		
  
	AUTHORS:

		Allen.Ingling@nyu.edu				awi
		mario.kleiner@tuebingen.mpg.de		mk
  
    PLATFORMS:
	
		All.

	HISTORY:
	
		mm/dd/yy
	
		10/10/03	awi		Created.  Based on SCREENFillRect and SCREENgluDisk.
		10/12/03	awi		Changed help string for new argument order.
		10/12/04	awi		In useString: changed "SCREEN" to "Screen", and moved commas to inside [].
		1/15/05		awi		Removed GL_BLEND setting a MK's suggestion.  
		2/25/05		awi		Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').

	TO DO:

*/

#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('FillOval', windowPtr [,color] [,rect] [,perfectUpToMaxDiameter]);";
static char synopsisString[] = 
        "Fills an ellipse with the given color, inscribed within \"rect\".\"color\" is the "
        "clut index (scalar or [r g b] triplet) that you want to poke into each pixel; "
        "default produces white with the standard CLUT for this window's pixelSize. "
        "Default rect is whole window.\n"
		"Instead of filling one oval, you can also specify a list of multiple ovals to be "
		"filled - this is much faster when you need to draw many ovals per frame. To fill n "
		"ovals, provide \"rect\" as a 4 rows by n columns matrix, each column specifying one "
		"oval, e.g., rect(1,5)=left border of 5th oval, rect(2,5)=top border of 5th oval, "
		"rect(3,5)=right border of 5th oval, rect(4,5)=bottom border of 5th oval. If the "
		"ovals should have different colors, then provide \"color\" as a 3 or 4 row by n column "
		"matrix, the i'th column specifiying the color of the i'th oval.\n"
		"The optional parameter 'perfectUpToMaxDiameter' allows to specify the maximum diameter "
		"in pixels for which a filled oval should look perfect. By default, the maximum diameter "
		"is chosen to be the full display size, so all ovals will look perfect, at a possible "
		"speed penalty. If you know your ovals will never be bigger than a certain diameter, "
		"you can provide that diameter as a hint via 'perfectUpToMaxDiameter' to allow for "
		"some potential speedup when drawing filled ovals. ";

static char seeAlsoString[] = "FrameOval";	

PsychError SCREENFillOval(void)  
{	
	PsychRectType			rect;
	double					numSlices, radius, xScale, yScale, xTranslate, yTranslate, rectY, rectX;
	PsychWindowRecordType	*windowRecord;
	psych_bool				isArgThere, isclassic;
    double					*xy, *colors;
	unsigned char			*bytecolors;
	int						numRects, i, nc, mc, nrsize;
	GLUquadricObj			*diskQuadric;
	double					perfectUpToMaxDiameter;
	static double			perfectUpToMaxDiameterOld = 0;

	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);}
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(4));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs

	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
    isclassic = PsychIsGLClassic(windowRecord);

	perfectUpToMaxDiameter = PsychGetWidthFromRect(windowRecord->clientrect);
	if (PsychGetHeightFromRect(windowRecord->clientrect) < perfectUpToMaxDiameter) perfectUpToMaxDiameter = PsychGetHeightFromRect(windowRecord->clientrect);
	PsychCopyInDoubleArg(4, kPsychArgOptional, &perfectUpToMaxDiameter);
	
    if ((perfectUpToMaxDiameter != perfectUpToMaxDiameterOld) || (windowRecord->fillOvalDisplayList == 0)) {
        perfectUpToMaxDiameterOld = perfectUpToMaxDiameter;

        // Compute number of subdivisions (slices) to provide a perfect oval, i.e., one subdivision for each
        // distance unit on the circumference of the oval.
        numSlices = 3.14159265358979323846 * perfectUpToMaxDiameter;

        // Destroy old display list so it gets rebuilt with the new numSlices setting:
        if (isclassic && (windowRecord->fillOvalDisplayList != 0)) {
            glDeleteLists(windowRecord->fillOvalDisplayList, 1);
            windowRecord->fillOvalDisplayList = 0;
        }
    }

    // Already cached display list for filled ovals for this windowRecord available?
    if (isclassic && (windowRecord->fillOvalDisplayList == 0)) {
        // Nope. Create our prototypical filled oval:
        // Generate a filled disk of that radius and subdivision and store it in a display list:
        diskQuadric=gluNewQuadric();
        windowRecord->fillOvalDisplayList = glGenLists(1);
        glNewList(windowRecord->fillOvalDisplayList, GL_COMPILE);
        gluDisk(diskQuadric, 0, 1, (int) numSlices, 1);
        glEndList();
        gluDeleteQuadric(diskQuadric);
        // Display list ready for use in this and all future drawing calls for this windowRecord.
    }

	// Query, allocate and copy in all vectors...
	numRects = 4;
	nrsize = 0;
	colors = NULL;
	bytecolors = NULL;
	mc = nc = 0;
	
	// The negative position -3 means: xy coords are expected at position 3, but they are optional.
	// NULL means - don't want a size's vector.
	PsychPrepareRenderBatch(windowRecord, -3, &numRects, &xy, 2, &nc, &mc, &colors, &bytecolors, 0, &nrsize, NULL, FALSE);

	// Only up to one rect provided?
	if (numRects <= 1) {
		// Get the oval and draw it:
		PsychCopyRect(rect, windowRecord->clientrect);
		isArgThere=PsychCopyInRectArg(kPsychUseDefaultArgPosition, FALSE, rect);
		if (isArgThere && IsPsychRectEmpty(rect)) return(PsychError_none);
		numRects = 1;
	}
	else {
		// Multiple ovals provided. Set up the first one:
		PsychCopyRect(rect, &xy[0]);
	}

	// Draw all ovals (one or multiple):
	for (i = 0; i < numRects;) {
		// Per oval color provided? If so then set it up. If only one common color
		// was provided then PsychPrepareRenderBatch() has already set it up.
		if (nc>1) {
			// Yes. Set color for this specific item:
			PsychSetArrayColor(windowRecord, i, mc, colors, bytecolors);
		}

		// Compute drawing parameters for ellipse:
		if (!IsPsychRectEmpty(rect)) {
			//The glu disk object location and size with a  center point and a radius,   
			//whereas FillOval accepts a bounding rect.   Converting from one set of parameters
			//to the other we should careful what we do for rects size of even number of pixels in length.
			PsychGetCenterFromRectAbsolute(rect, &xTranslate, &yTranslate);
			rectY=PsychGetHeightFromRect(rect);
			rectX=PsychGetWidthFromRect(rect);
			if(rectX == rectY){
				xScale=1; 
				yScale=1;
				radius=rectX/2;
			}else if(rectX > rectY){ 
				xScale=1;
				yScale=rectY/rectX;
				radius=rectX/2;
			}else if(rectY > rectX){
				yScale=1;
				xScale=rectX/rectY;
				radius=rectY/2;
			}

            if (isclassic) {
                // Draw: Set up position, scale and size via matrix transform:
                glPushMatrix();
                glTranslatef((float) xTranslate, (float) yTranslate, (float) 0);
                glScalef((float) (xScale * radius), (float) (yScale * radius), (float) 1);
                // Draw cached disk object (stored in display list):
                glCallList(windowRecord->fillOvalDisplayList);
                glPopMatrix();
            }
            else {
                PsychDrawDisc(windowRecord, (float) xTranslate, (float) yTranslate, (float) 0, (float) radius, (int) numSlices, (float) xScale, (float) yScale, 0, 360);
            }
		}
		
		// Done with this one. Set up the next one, if any...
		i++;
		if (i < numRects) PsychCopyRect(rect, &xy[i*4]);

		// Next oval.
	}
	
	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);

 	//All psychfunctions require this.
	return(PsychError_none);
}
