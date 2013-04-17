/*
	SCREENFramePoly.c		

	AUTHORS:

	Allen.Ingling@nyu.edu			awi 
	mario.kleiner@tuebingen.mpg.de	mk

	PLATFORMS:	
	
	All.

	HISTORY:

		07/24/04	awi		Created.
		10/12/04	awi		In useString: moved commas to inside [].
		2/25/05		awi		Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
		
	TO DO:

	NOTES:

		It would be faster to use one call to glVertex2Dv instead of looping over glVertex2D except that we would have to reshape the matrix. 
		We have to reshape the matrix because glVertex assumes that values are stored in memory pairwise but MATLAB packs elements columnwise  
		and FillPolly assumes that a vertex is a row.  We could migrate to the new order except that 2x2 would be ambiguous.  If speed ever 
		becomes an issue then we could supply a preference to specify vertex columns instead of vertex rows and use glVertex2D.

*/


#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('FramePoly', windowPtr [,color], pointList [,penWidth]);";
//                                             1           2       3           4
static char synopsisString[] = 
	"Draw a polygon frame. \"color\" is the clut index (scalar or [r g b a] vector) that you "
	"want to poke into each pixel; default produces white. \"pointList\" is a matrix: each row specifies the (x,y) "
	"coordinates of a vertex. ";
	
static char seeAlsoString[] = "FillPoly";	

PsychError SCREENFramePoly(void)  
{
	
	PsychColorType					color;
	PsychWindowRecordType			*windowRecord;
	double							whiteValue;
	int								i, mSize, nSize, pSize;
	psych_bool						isArgThere;
	double							penSize, *pointList;
    
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
	if(nSize!=2)
		PsychErrorExitMsg(PsychError_user, "Width of point list must be 2");
	if(mSize<3)
		PsychErrorExitMsg(PsychError_user, "Polygons must consist of at leat 3 vertices; M dimension of pointList was <3");
	if(pSize>1)
		PsychErrorExitMsg(PsychError_user, "Polygons must consist of at leat 3 vertices; M dimension of pointList was <3");
	
	//get the pen size
	penSize=1;
	PsychCopyInDoubleArg(4, kPsychArgOptional, &penSize);
	
	// Enable this windowRecords framebuffer as current drawingtarget:
	PsychSetDrawingTarget(windowRecord);

	// Set default draw shader:
	PsychSetShader(windowRecord, -1);
	
	glLineWidth((GLfloat)penSize);

	PsychUpdateAlphaBlendingFactorLazily(windowRecord);
	PsychSetGLColor(&color, windowRecord);
	GLBEGIN(GL_LINE_LOOP);
		for(i = 0; i < mSize; i++)
			GLVERTEX2d((GLdouble)pointList[i], (GLdouble)pointList[i+mSize]);
	GLEND();

	glLineWidth((GLfloat) 1);

	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);

	return(PsychError_none);
}
