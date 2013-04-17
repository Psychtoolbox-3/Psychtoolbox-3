/*
	SCREENDrawDots.c	
 
	AUTHORS:
 
		kas@princeton.edu				kas		Keith Schneider
		fcalabro@bu.edu					fjc		Finnegan Calabro
		mario.kleiner@tuebingen.mpg.de  mk		Mario Kleiner
 
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
		1/30/05     mk      Bug-fixes: Adding PsychMallocTemp()'s and a few checks to avoid memory corruption
							crashes.
		2/25/05     awi     Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
		3/22/05     mk      Added possibility to spec vectors with individual color and size spec per dot.
		4/29/05     mk      Bugfix for color vectors: They should also take values in range 0-255 instead of 0.0-1.0.
		11/14/06    mk      We now also accept color vectors in uint8 format and pass them directly for higher efficiency.
		
	TO DO:
 
 */


#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('DrawDots', windowPtr, xy [,size] [,color] [,center] [,dot_type]);";
//                                            1          2    3       4        5         6
static char synopsisString[] = 
"Quickly draw an array of dots.  "
"\"xy\" is a two-row vector containing the x and y coordinates of the dot centers, "
"relative to \"center\" (default center is [0 0]).  "
"\"size\" is the width of each dot in pixels (default is 1). "
"Instead of a common size for all dots you can also provide a "
"vector which defines a different dot size for each dot.  "
"\"color\" is the the clut index (scalar or [r g b a] vector) "
"that you want to poke into each dot pixel (default is black).  "
"Instead of a single \"color\" you can also provide a 3 or 4 row vector,"
"which specifies an individual RGB or RGBA color for each corresponding point."
"\"dot_type\" is a flag that determines what type of dot is drawn: "
"0 (default) squares, 1 circles (with anti-aliasing), 2 circles (with high-quality "
"anti-aliasing, if supported by your hardware). "
"If you use dot_type = 1 you'll also need to set a proper blending mode with the "
"Screen('BlendFunction') command!";  
static char seeAlsoString[] = "BlendFunction";	 

PsychError SCREENDrawDots(void)  
{
	PsychWindowRecordType                   *windowRecord;
	int                                     m,n,p,mc,nc,idot_type;
	int                                     i, nrpoints, nrsize;
	psych_bool                              isArgThere, usecolorvector, isdoublecolors, isuint8colors;
	double									*xy, *size, *center, *dot_type, *colors;
    float                                   *sizef;
	unsigned char                           *bytecolors;
	GLfloat									pointsizerange[2];
    psych_bool                              lenient = FALSE;

	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	// Check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(6));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs
	
	// Get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(1, kPsychArgRequired, &windowRecord);
	
	// Query, allocate and copy in all vectors...
	nrpoints = 2;
	nrsize = 0;
	colors = NULL;
	bytecolors = NULL;

	PsychPrepareRenderBatch(windowRecord, 2, &nrpoints, &xy, 4, &nc, &mc, &colors, &bytecolors, 3, &nrsize, &size, (GL_FLOAT == PsychGLFloatType(windowRecord)));
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
	
	// Get dot_type argument
	isArgThere = PsychIsArgPresent(PsychArgIn, 6);
	if(!isArgThere){
		idot_type = 0;
	} else {
		PsychAllocInDoubleMatArg(6, TRUE, &m, &n, &p, &dot_type);
		if(p!=1 || n!=1 || m!=1 || (dot_type[0]<0 || dot_type[0]>2))
			PsychErrorExitMsg(PsychError_user, "dot_type must be 0, 1 or 2");
		idot_type = (int) dot_type[0];
	}
	
	// Turn on antialiasing to draw circles	
	if(idot_type) {
		glEnable(GL_POINT_SMOOTH);
		glGetFloatv(GL_POINT_SIZE_RANGE, (GLfloat*) &pointsizerange);
		// A dot type of 2 requests for highest quality point smoothing:
		glHint(GL_POINT_SMOOTH_HINT, (idot_type>1) ? GL_NICEST : GL_DONT_CARE);
	}
	else {
		glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, (GLfloat*) &pointsizerange);
	}
	
    // Does ES-GPU only support a fixed point diameter of 1 pixel?
    if ((pointsizerange[1] <= 1) && PsychIsGLES(windowRecord)) {
        // Yes. Not much point bailing on this, as it should be easily visible
        // during testing of a studies code on a OpenGL-ES device.
        lenient = TRUE;
    }

	// Set size of a single dot:
	if (!lenient && ((sizef && (sizef[0] > pointsizerange[1] || sizef[0] < pointsizerange[0])) ||
                     (!sizef && (size[0] > pointsizerange[1] || size[0] < pointsizerange[0])))) {
		printf("PTB-ERROR: You requested a point size of %f units, which is not in the range (%f to %f) supported by your graphics hardware.\n",
			   (sizef) ? sizef[0] : size[0], pointsizerange[0], pointsizerange[1]);
		PsychErrorExitMsg(PsychError_user, "Unsupported point size requested in Screen('DrawDots').");
	}
	
	// Setup initial common point size for all points:
    if (!lenient) glPointSize((sizef) ? sizef[0] : (float) size[0]);
	
	// Setup modelview matrix to perform translation by 'center':
	glMatrixMode(GL_MODELVIEW);
	
	// Make a backup copy of the matrix:
	glPushMatrix();
	
	// Apply a global translation of (center(x,y)) pixels to all following points:
	glTranslatef((float) center[0], (float) center[1], 0);
	
	// Render the array of 2D-Points - Efficient version:
	// This command sequence allows fast processing of whole arrays
	// of vertices (or points, in this case). It saves the call overhead
	// associated with the original implementation below and is potentially
	// optimized in specific OpenGL implementations.
	
	// Pass a pointer to the start of the point-coordinate array:
	glVertexPointer(2, PSYCHGLFLOAT, 0, &xy[0]);
	
	// Enable fast rendering of arrays:
	glEnableClientState(GL_VERTEX_ARRAY);
	
	if (usecolorvector) {
		PsychSetupVertexColorArrays(windowRecord, TRUE, mc, colors, bytecolors);
	}
	
	// Render all n points, starting at point 0, render them as POINTS:
	if (nrsize==1) {
		// One common point size for all dots provided. Good! This is very efficiently
		// done with one single render-call:
		glDrawArrays(GL_POINTS, 0, nrpoints);
	}
	else {
		// Different size for each dot provided: We have to do One GL - call per dot.
		// This is *pretty inefficient* and should be reimplemented in the future via
		// Point-Sprite extensions, cleverly used display lists or via vertex-shaders...
		// For now we do it the stupid way:
		for (i=0; i<nrpoints; i++) {
            if (!lenient && ((sizef && (sizef[i] > pointsizerange[1] || sizef[i] < pointsizerange[0])) ||
                             (!sizef && (size[i] > pointsizerange[1] || size[i] < pointsizerange[0])))) {
				printf("PTB-ERROR: You requested a point size of %f units, which is not in the range (%f to %f) supported by your graphics hardware.\n",
					   (sizef) ? sizef[i] : size[i], pointsizerange[0], pointsizerange[1]);
				PsychErrorExitMsg(PsychError_user, "Unsupported point size requested in Screen('DrawDots').");
			}
			
			// Setup point size for this point:
			if (!lenient) glPointSize((sizef) ? sizef[i] : (float) size[i]);
			
			// Render point:
			glDrawArrays(GL_POINTS, i, 1);
		}
	}
	
	// Disable fast rendering of arrays:
	glDisableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2, PSYCHGLFLOAT, 0, NULL);
	
	if (usecolorvector) PsychSetupVertexColorArrays(windowRecord, FALSE, 0, NULL, NULL);
	
	// Restore old matrix from backup copy, undoing the global translation:
	glPopMatrix();
	
	// turn off antialiasing again
	if(idot_type) glDisable(GL_POINT_SMOOTH);
	
	// Reset pointsize to 1.0
	glPointSize(1);
	
	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);
	
 	//All psychfunctions require this.
	return(PsychError_none);
}
