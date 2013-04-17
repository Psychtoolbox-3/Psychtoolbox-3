/*
	SCREENFillPoly.c		
 
	AUTHORS:
 
	Allen.Ingling@nyu.edu			awi 
	mario.kleiner@tuebingen.mpg.de	mk

	PLATFORMS:	
	
	All.
 
	HISTORY:
 
	07/23/04	awi		Created.
	10/12/04	awi		In useString: moved commas to inside [].
	01/12/05     mk     Added a slow-path that draws concave and self-intersecting polygons correctly.
	02/25/05	awi		Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
	11/01/08	 mk		Improved speed of slow-path. Still pretty slow -> Most time spent inside gluTesselator(), nothing we could do.
 
	TO DO:
 
	NOTES:
 
	 It would be faster to use one call to glVertex2Dv instead of looping over glVertex2D except that we would have to reshape the matrix. 
	 We have to reshape the matrix because glVertex assumes that values are stored in memory pairwise but MATLAB packs elements columnwise  
	 and FillPolly assumes that a vertex is a row.  We could migrate to the new order except that 2x2 would be ambiguous.  If speed ever 
	 becomes an issue then we could supply a preference to specify vertex columns instead of vertex rows and use glVertex2D.
 */

#include "Screen.h"

#if PSYCH_SYSTEM == PSYCH_LINUX
#define GLUTESSCBCASTER (_GLUfuncptr)
#else
#define GLUTESSCBCASTER
#endif

// Our pointer to the used GLUtesselator object:
static GLUtesselator		*tess = NULL;

// Cache datastructures to reduce or avoid many memory allocation/deallocations:
static GLdouble				*combinerCache = NULL;
static unsigned int			combinerCacheSize = 0;
static unsigned int			combinerCacheSlot = 0;

static double*				tempv = NULL;
static int					tempvsize = 0;

// Callback-Routines for the GLU-Tesselator functions used on the FillPoly - Slow - path:
void APIENTRY PsychtcbBegin(GLenum prim)
{
    PsychGLBegin(NULL, prim);
}

void APIENTRY PsychtcbVertex(void *data)
{
    GLdouble *v = (GLdouble*) data;
    PsychGLVertex4f(NULL, (float) v[0], (float) v[1], (float) v[2], (float) 1);
}

void APIENTRY PsychtcbEnd(void)
{
    PsychGLEnd(NULL);
}

void APIENTRY PsychtcbCombine(GLdouble c[3], void *d[4], GLfloat w[4], void **out)
{
    GLdouble *nv;

	// Free slots available?
	if (combinerCacheSlot >= combinerCacheSize) {
		// Nope. Need to alloc another cache for up to another 1000 elements:
		combinerCacheSize = 1000;
		combinerCacheSlot = 0;
		combinerCache = (GLdouble *) PsychMallocTemp(sizeof(GLdouble) * 3 * combinerCacheSize);
		if (NULL == combinerCache) PsychErrorExitMsg(PsychError_outofMemory, "Out of memory condition in Screen('FillPoly')! Not enough space.");
	}

    nv = (GLdouble *) &(combinerCache[combinerCacheSlot * 3]);
    nv[0] = c[0];
    nv[1] = c[1];
    nv[2] = c[2];
    *out = nv;

	combinerCacheSlot++;
}

// Cleanup routine for our tesselators and other data structures. Called from
// ScreenExit.c at Screen shutdown. May be called without OpenGL active! Don't
// use any GL calls here, just plain C-level operations!!
void PsychCleanupSCREENFillPoly(void)
{
	// Release tesselator object and associated data structures, if any:
	if (tess) {
		gluDeleteTess(tess);
		tess=NULL;
	}

	if (tempv) {
		free(tempv);
		tempv = NULL;
		tempvsize = 0;
	}
	
	return;
}

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('FillPoly', windowPtr [,color], pointList [, isConvex]);";
//                                            1           2       3			4
static char synopsisString[] = 
"Fill polygon. \"color\" is the clut index (scalar or [r g b] or [r g b a] vector) that you "
"want to poke into each pixel; default produces white. \"pointList\" is a matrix: each row specifies the (x,y) "
"coordinates of a vertex.\n"
"The optional flag 'isConvex' allows you to tell the routine if the polygon is convex, (value of 1) "
"or if it is concave (value of 0), allowing the routine to skip the test for convexity, thereby "
"saving a bit of computation time in timing sensitive scripts.\n"
"Drawing filled polygons is a rather compute intense and slow operation. In general, drawing "
"convex polygons is pretty fast, drawing concave, but not self-intersecting polygons is much "
"slower (> 30x slower), and drawing concave self-intersecting polygons is extremely slow, "
"e.g. > 2000 times slower. If you have to draw very irregular polygons with may points, "
"it might be a good idea to preprocess them in some way and maybe break them up into "
"a sequence of more convex/regular polygons before submitting them to 'FillPoly'. Or "
"you may want to use some custom written drawing function for your purpose which is "
"optimized for drawing your type of polygons. ";

static char seeAlsoString[] = "FramePoly";	

PsychError SCREENFillPoly(void)  
{	
	PsychColorType				color;
	PsychWindowRecordType		*windowRecord;
	double						whiteValue;
	int							i, mSize, nSize, pSize;
	psych_bool					isArgThere;
	double						*pointList;
	double						isConvex;
	int							j,k;
	int							flag;
	double						z;
	
	combinerCacheSlot = 0;
	combinerCacheSize = 0;
	combinerCache = NULL;
	
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(4));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs
	
	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(1, kPsychArgRequired, &windowRecord);
	
	//Get the color argument or use the default, then coerce to the form determened by the window depth.  
	isArgThere=PsychCopyInColorArg(2, FALSE, &color);
	if(!isArgThere){
		whiteValue=PsychGetWhiteValueFromWindow(windowRecord);
		PsychLoadColorStruct(&color, kPsychIndexColor, whiteValue ); //index mode will coerce to any other.
	}
 	PsychCoerceColorMode( &color);
	
	//get the list of pairs and validate.  
	PsychAllocInDoubleMatArg(3, kPsychArgRequired, &mSize, &nSize, &pSize, &pointList);
	if(nSize!=2) PsychErrorExitMsg(PsychError_user, "Width of pointList must be 2");
	if(mSize<3)  PsychErrorExitMsg(PsychError_user, "Polygons must consist of at least 3 points; M dimension of pointList was < 3!");
	if(pSize>1)  PsychErrorExitMsg(PsychError_user, "pointList must be a 2D matrix, not a 3D matrix!");
	
	isConvex = -1;
	PsychCopyInDoubleArg(4, kPsychArgOptional, &isConvex);
	
    // On non-OpenGL1/2 we always force isConvex to zero, so the GLU tesselator is
    // always used. This because the tesselator only emits GL_TRIANGLES and GL_TRIANGLE_STRIP
    // and GL_TRIANGLE_FANS primitives which are supported on all current OpenGL API's, whereas
    // or "classic" fast-path needs GL_POLYGONS, which are only supported on classic OpenGL1/2:
    if (!PsychIsGLClassic(windowRecord)) isConvex = 0;

	// Enable this windowRecords framebuffer as current drawingtarget:
	PsychSetDrawingTarget(windowRecord);
	
	// Set default drawshader:
	PsychSetShader(windowRecord, -1);
	
	PsychUpdateAlphaBlendingFactorLazily(windowRecord);		 
	PsychSetGLColor(&color, windowRecord);
	
	///////// Test for convexity ////////
	// This algorithm checks, if the polygon is definitely convex, or not.
	// We take the slow-path, if polygon is non-convex or if we can't prove
	// that it is convex.
	//
	// Algorithm adapted from: http://astronomy.swin.edu.au/~pbourke/geometry/clockwise/
	// Which was written by Paul Bourke, 1998.
	//
	// -> This webpage explains the mathematical principle behind the test and provides
	// a C-Source file which has been adapted for use here.
	//	
	if (isConvex == -1) {
		flag = 0;
		for (i=0; i < mSize; i++) {
			j = (i + 1) % mSize;
			k = (i + 2) % mSize;
			z  = (pointList[j] - pointList[i]) * (pointList[k+mSize] - pointList[j+mSize]);
			z -= (pointList[j+mSize] - pointList[i+mSize]) * (pointList[k] - pointList[j]);
			
			if (z < 0) {
				flag |= 1;
			}
			else if (z > 0) {
				flag |= 2;
			}
			
			if (flag == 3) {
				// This is definitely a CONCAVE polygon --> not Convex --> Take slow but safe path.
				break;
			}
		}
		
		if (flag!=0 && flag!=3) {
			// This is a convex polygon --> Take fast path.
			isConvex = 1;
		}
		else {
			// This is a complex polygon --> can't determine if it is convex or not --> Take slow but safe path.
			isConvex = 0;
		}
	}
			
	////// Switch between fast path and slow path, depending on convexity of polygon:
	if (isConvex > 0) {
		// Convex, non-self-intersecting polygon - Take the fast-path:
		glBegin(GL_POLYGON);
		for(i=0;i<mSize;i++) glVertex2d((GLdouble)pointList[i], (GLdouble)pointList[i+mSize]);
		glEnd();
	}
	else {
		// Possibly concave and/or self-intersecting polygon - At least we couldn't prove it is convex.
		// Take the slow, but safe, path using GLU-Tesselators to break it up into a couple of convex, simple
		// polygons:
		
		// Create and initialize a new GLU-Tesselator object, if needed:
		if (NULL == tess) {
			// Create tesselator:
			tess = gluNewTess();
			if (NULL == tess) PsychErrorExitMsg(PsychError_outofMemory, "Out of memory condition in Screen('FillPoly')! Not enough space.");

			// Assign our callback-functions:
			gluTessCallback(tess, GLU_TESS_BEGIN, GLUTESSCBCASTER PsychtcbBegin);
			gluTessCallback(tess, GLU_TESS_VERTEX, GLUTESSCBCASTER PsychtcbVertex);
			gluTessCallback(tess, GLU_TESS_END, GLUTESSCBCASTER PsychtcbEnd);
			gluTessCallback(tess, GLU_TESS_COMBINE, GLUTESSCBCASTER PsychtcbCombine);

			// Define all tesselated polygons to lie in the x-y plane:
			gluTessNormal(tess, 0, 0, 1);
		}

		// We need to hold the values in a temporary array:
		if (tempvsize < mSize) {
			tempvsize = ((mSize / 1000) + 1) * 1000;
			tempv = (double*) realloc((void*) tempv, sizeof(double) * 3 * tempvsize);
			if (NULL == tempv) PsychErrorExitMsg(PsychError_outofMemory, "Out of memory condition in Screen('FillPoly')! Not enough space.");
		}

		// Now submit our Polygon for tesselation:
		gluTessBeginPolygon(tess, NULL);
		gluTessBeginContour(tess);

		for(i=0; i < mSize; i++) {
			tempv[i*3]=(GLdouble) pointList[i];
			tempv[i*3+1]=(GLdouble) pointList[i+mSize];
			tempv[i*3+2]=0;
			gluTessVertex(tess, (GLdouble*) &(tempv[i*3]), (void*) &(tempv[i*3]));
		}
		
		// Process, finalize and render it by calling our callback-functions:
		gluTessEndContour(tess);
		gluTessEndPolygon (tess);
		
		// Done with drawing the filled polygon. (Slow-Path)
	}
	
	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);
	
	// printf("CombinerCalls %i out of %i allocated.\n", combinerCacheSlot, combinerCacheSize);

	return(PsychError_none);
}
