/*
SCREENDrawLines.c	
 
	AUTHORS:
 
		kas@princeton.edu				kas		Keith Schneider
		fcalabro@bu.edu					fjc		Finnegan Calabro
		mario.kleiner@tuebingen.mpg.de                  mk		Mario Kleiner
		dtaylor@ski.org					dgt		Douglas Taylor
 
	PLATFORMS:	
	
		All.  
 
	HISTORY:
	
		mm/dd/yy   
 
		12/13/04	kas		Created  
		12/17/04	fjc		added antialiased points
		12/17/04	kas		added dot_type flag
		1/11/04		awi		Merged into pychtoolbox.org distribution.
		1/01/05		mk		Performance optimizations and small bug-fixes
								- Fixes some minor bugs (or possible bugs, they would only
								trigger in some usage cases): GL-Calls were missing to undo
								parameter changes like glPointSize,
								glEnable(GL_POINT_SMOOTH). This can affect later
								invocations of DrawDots, e.g., Call DrawDots with point
								smoothing on --> Smoothing gets enabled, after that it
								would not be possible to disable point smoothing again -->
								any following DrawDots call would always render smooth
								points, regardless of setting.
								
								- Uses an optimized way of drawing the dots: Instead of
								looping over the 2-column vector of dots in a for-loop and
								submitting each Point via glVertex2d(), it now *requires* a
								2-row vector of dots as input and then renders all dots
								with one function call to glDrawArrays(). This is the most
								efficient way of passing a huge number of primitives
								(Points, rectangles, ...) to the OpenGL-System, avoiding
								for-loop and function call overhead and allowing possible
								optimizations by the graphics driver.								
		1/12/05		awi		Merged mk's improvements into psychtoolbox.org distribution.
		3/8/05		dgt		Created it from SCREENDrawDots.c
								- Uses the same vertex array method of drawing a lot of lines.
								Line drawing is not as optimised as point drawing but this
								method has the advantage of passing an rgb or rgba color
								array to use as the lines are drawn.
		3/23/05		dgt		Merged mk's improvements into SCREENDrawLines.c
		4/22/05     mk      Small bug fix (size = PsychMallocTemp.....)
		12/4/06		mk		Rewrite to make it functional again and to implement a similar
							syntax to Screen('DrawDots').

 */

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('DrawLines', windowPtr, xy [,width] [,colors] [,center] [,smooth]);";
//                                             1          2    3        4         5         6
static char synopsisString[] = 
"Quickly draw an array of lines into the specified window \"windowPtr\".\n"
"\"xy\" is a two-row vector containing the x and y coordinates of the line segments: Pairs of consecutive "
"columns define (x,y) positions of the starts and ends of line segments. All positions are relative "
"to \"center\" (default center is [0 0]). \"width\" is either a scalar with the global width for "
"all lines in pixels (default is 1), or a vector with one separate width value for each separate line. "
"\"colors\" is either a single global color argument for all lines, or an array of rgb or rgba "
"color values for each line, where each column corresponds to the color of the corresponding line start or "
"endpoint in the xy position argument. If you specify different colors for the start- and endpoint of a "
"line segment, PTB will generate a smooth transition of colors along the line via linear interpolation. "
"The default color is white if colors is omitted. \"smooth\" is a flag that determines whether lines "
"should be smoothed: 0 (default) no smoothing, 1 smoothing (with anti-aliasing), 2 = high quality smoothing. "
"If you use smoothing, you'll also need to set a proper blending mode with Screen('BlendFunction').";

static char seeAlsoString[] = "BlendFunction";	 

PsychError SCREENDrawLines(void)  
{
	PsychWindowRecordType		*windowRecord;
	int							m,n,p, smooth;
	int							nrsize, nrvertices, mc, nc, i;
	psych_bool                  isArgThere, usecolorvector, isdoublecolors, isuint8colors;
	double						*xy, *size, *center, *dot_type, *colors;
	unsigned char               *bytecolors;
	float						linesizerange[2];
    float						*sizef;
    psych_bool                  lenient = FALSE;

	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(6));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs
	
	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(1, kPsychArgRequired, &windowRecord);
	
	// Query, allocate and copy in all vectors...
	nrvertices = 2;
	nrsize = 1;

	colors = NULL;
	bytecolors = NULL;

	PsychPrepareRenderBatch(windowRecord, 2, &nrvertices, &xy, 4, &nc, &mc, &colors, &bytecolors, 3, &nrsize, &size, (GL_FLOAT == PsychGLFloatType(windowRecord)));
	isdoublecolors = (colors) ? TRUE:FALSE;
	isuint8colors  = (bytecolors) ? TRUE:FALSE;
	usecolorvector = (nc>1) ? TRUE:FALSE;

    // Assign sizef as float-type array of sizes, if float mode active, NULL otherwise:
    sizef = (GL_FLOAT == PsychGLFloatType(windowRecord)) ? (float*) size : NULL;

	// Get center argument
	isArgThere = PsychIsArgPresent(PsychArgIn, 5);
	if(!isArgThere){
		center = (double *) PsychMallocTemp(2 * sizeof(double));
		center[0] = 0;
		center[1] = 0;
	} else {
		PsychAllocInDoubleMatArg(5, TRUE, &m, &n, &p, &center);
		if(p!=1 || n!=2 || m!=1) PsychErrorExitMsg(PsychError_user, "center must be a 1-by-2 vector");
	}
	
	// Get smooth argument
	isArgThere = PsychIsArgPresent(PsychArgIn, 6);
	if(!isArgThere){
		smooth = 0;
	} else {
		PsychAllocInDoubleMatArg(6, TRUE, &m, &n, &p, &dot_type);
		smooth = (int) dot_type[0];
		if(p!=1 || n!=1 || m!=1 || (smooth < 0 || smooth > 2)) PsychErrorExitMsg(PsychError_user, "smooth must be 0, 1 or 2.");
	}

	// Child-protection: Alpha blending needs to be enabled for smoothing to work:
	if (smooth > 0 && windowRecord->actualEnableBlending!=TRUE) {
		PsychErrorExitMsg(PsychError_user, "Line smoothing doesn't work with alpha-blending disabled! See Screen('BlendFunction') on how to enable it.");
	}

	// turn on antialiasing to draw anti-aliased lines:
	if(smooth) {
        glGetFloatv(GL_SMOOTH_LINE_WIDTH_RANGE, (GLfloat*) &linesizerange);
        glEnable(GL_LINE_SMOOTH);
		// A smooth level of 2 requests highest quality line smoothing:
		glHint(GL_LINE_SMOOTH_HINT, (smooth > 1) ? GL_NICEST : GL_DONT_CARE);
    }
    else {
        glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, (GLfloat*) &linesizerange);
    }

    // Does ES-GPU only support a fixed line width of 1 pixel?
    if ((linesizerange[0] == linesizerange[1]) && (linesizerange[0] <= 1) && PsychIsGLES(windowRecord)) {
        // Yes. Not much point bailing on this, as it should be easily visible
        // during testing of a studies code on a OpenGL-ES device.
        lenient = TRUE;
    }

	if (!lenient && ((sizef && (sizef[0] > linesizerange[1] || sizef[0] < linesizerange[0])) ||
                     (!sizef && (size[0] > linesizerange[1] || size[0] < linesizerange[0])))) {
		printf("PTB-ERROR: You requested a line width of %f units, which is not in the range (%f to %f) supported by your graphics hardware.\n",
			   (sizef) ? sizef[0] : size[0], linesizerange[0], linesizerange[1]);
		PsychErrorExitMsg(PsychError_user, "Unsupported line width requested.");
	}

	// Set global width of lines:
	if (!lenient) glLineWidth((sizef) ? sizef[0] : (float) size[0]);

	// Setup modelview matrix to perform translation by 'center':
	glMatrixMode(GL_MODELVIEW);	
	
	// Make a backup copy of the matrix:
	glPushMatrix();
	
	// Apply a global translation of (center(x,y)) pixels to all following lines:
	glTranslatef((float) center[0], (float) center[1], (float) 0);
	
	// Render the array of 2D-Lines - Efficient version:
	// This command sequence allows fast processing of whole arrays
	// of vertices (or lines, in this case). It saves the call overhead
	// associated with the original implementation below and is potentially
	// optimized in specific OpenGL implementations.
	
	// Pass a pointer to the start of the arrays:
	glVertexPointer(2, PSYCHGLFLOAT, 0, &xy[0]);
	
	if (usecolorvector) {
		PsychSetupVertexColorArrays(windowRecord, TRUE, mc, colors, bytecolors);
	}

	// Enable fast rendering of arrays:
	glEnableClientState(GL_VERTEX_ARRAY);

	if (nrsize==1) {
		// Common line-width for all lines: Render all lines, starting at line 0:
		glDrawArrays(GL_LINES, 0, nrvertices);
	}
	else {
		// Different line-width per line: Need to manually loop through this mess:
		for (i=0; i < nrvertices/2; i++) {
            if (!lenient && ((sizef && (sizef[i] > linesizerange[1] || sizef[i] < linesizerange[0])) ||
                             (!sizef && (size[i] > linesizerange[1] || size[i] < linesizerange[0])))) {
                printf("PTB-ERROR: You requested a line width of %f units, which is not in the range (%f to %f) supported by your graphics hardware.\n",
                       (sizef) ? sizef[i] : size[i], linesizerange[0], linesizerange[1]);
                PsychErrorExitMsg(PsychError_user, "Unsupported line width requested.");
            }
            
            if (!lenient) glLineWidth((sizef) ? sizef[i] : (float) size[i]);

            // Render line:
            glDrawArrays(GL_LINES, i * 2, 2);
		}
	}
	
	// Disable fast rendering of arrays:
	glDisableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2, PSYCHGLFLOAT, 0, NULL);

	if (usecolorvector) PsychSetupVertexColorArrays(windowRecord, FALSE, 0, NULL, NULL);
	
	// Restore old matrix from backup copy, undoing the global translation:
	glPopMatrix();
	
	// Turn off anti-aliasing:
	if(smooth) glDisable(GL_LINE_SMOOTH);
	
	// Reset line width to 1.0:
	glLineWidth(1);
	
	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);
	
 	//All psychfunctions require this.
	return(PsychError_none);
}
