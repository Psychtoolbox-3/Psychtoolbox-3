/*
	SCREENDrawLine.c		

	AUTHORS:

	Allen.Ingling@nyu.edu           awi
    mario.kleiner@tuebingen.mpg.de  mk

	PLATFORMS:	

	All.

	HISTORY:

		07/23/04	awi		Created.   
		10/12/04	awi		In useString: moved commas to inside [].
		11/16/04    awi		Fixed a bug in DrawLine where dX should have been dY.
		2/25/05		awi		Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').

	TO DO:

*/

#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('DrawLine', windowPtr [,color], fromH, fromV, toH, toV [,penWidth]);";
//                                            1           2       3      4      5    6     7          
static char synopsisString[] = 
    "Draw a line. \"color\" is the clut index (scalar or [r g b a] vector) that "
    "you want to poke into each pixel; default produces black. "
    "\"fromH\" and \"fromV\" are the starting x and y positions, respectively. "
    "\"toH\" and \"toV\" are the ending x and y positions, respectively. "
    "Default \"penWidth\" is 1. ";
	
static char seeAlsoString[] = "DrawLines";	

PsychError SCREENDrawLine(void)  
{
	
	PsychColorType					color;
	PsychWindowRecordType			*windowRecord;
	double							whiteValue;
	psych_bool						isArgThere;
	double							sX, sY, dX, dY, penSize;
	float                           linesizerange[2];

	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(7));   //The maximum number of inputs
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
        
	//get source and destination X and Y values
	PsychCopyInDoubleArg(3, kPsychArgRequired, &sX);
	PsychCopyInDoubleArg(4, kPsychArgRequired, &sY);
	PsychCopyInDoubleArg(5, kPsychArgRequired, &dX);
	PsychCopyInDoubleArg(6, kPsychArgRequired, &dY);
	
	//get and set the pen size
	penSize=1;
	PsychCopyInDoubleArg(7, kPsychArgOptional, &penSize);
	
	// Enable this windowRecords framebuffer as current drawingtarget:
	PsychSetDrawingTarget(windowRecord);

	// Set default draw shader:
	PsychSetShader(windowRecord, -1);

    glGetFloatv(GL_LINE_WIDTH_RANGE, (GLfloat*) &linesizerange);
    if (penSize < linesizerange[0] || penSize > linesizerange[1]) {
        // Does ES-GPU only support a fixed line width of 1 pixel?
        if ((linesizerange[0] == linesizerange[1]) && (linesizerange[0] <= 1) && PsychIsGLES(windowRecord)) {
            // Yes. Not much point bailing on this, as it should be easily visible
            // during testing of a studies code on a OpenGL-ES device.
            penSize = 1;
        }
        else {
            printf("PTB-ERROR: You requested a line width of %f units, which is not in the range (%f to %f) supported by your graphics hardware.\n",
                   penSize, linesizerange[0], linesizerange[1]);
            PsychErrorExitMsg(PsychError_user, "Unsupported line width requested.");
        }
	}

	glLineWidth((GLfloat) penSize);

	PsychUpdateAlphaBlendingFactorLazily(windowRecord);
	PsychSetGLColor(&color, windowRecord);
	GLBEGIN(GL_LINES);
		GLVERTEX2f((GLfloat)sX, (GLfloat)sY);
		GLVERTEX2f((GLfloat)dX, (GLfloat)dY);
	GLEND();
	
	glLineWidth((GLfloat) 1);

	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);

	return(PsychError_none);
}
