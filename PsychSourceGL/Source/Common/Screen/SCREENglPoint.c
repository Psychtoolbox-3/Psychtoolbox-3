/*
	SCREENglPoint.c		
  
	AUTHORS:

	Allen.Ingling@nyu.edu			awi 
	mario.kleiner@tuebingen.mpg.de	mk

	PLATFORMS:	
	
	All.

	HISTORY:
  
		10/10/03	awi		Created.  It's experimental.
		1/14/05		awi		Added calls to undo glpoint smooth and point size changes at MK's suggestion.  
		2/25/05		awi		Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
   
	TO DO:
  

*/

#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "Screen('glPoint', windowPtr, color, x, y [,size]);";
static char synopsisString[] = 
	"Draw a point at the specified location 'x', 'y' with the requested optional diameter 'size' "
	"in the requested 'color'. For drawing of many points, e.g., random dot fields, use "
	"Screen('DrawDots') instead. ";

static char seeAlsoString[] = "FrameRect";	

PsychError SCREENglPoint(void)  
{
	
	PsychColorType			color;
	double					*xPosition, *yPosition, dotSize;
	PsychWindowRecordType	*windowRecord;
	double					whiteValue;
	psych_bool				isArgThere;
    
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(5));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs

	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
        
	//Get the color argument or use the default, then coerce to the form determened by the window depth.  
	isArgThere=PsychCopyInColorArg(kPsychUseDefaultArgPosition, FALSE, &color);
		if(!isArgThere){
			whiteValue=PsychGetWhiteValueFromWindow(windowRecord);
			PsychLoadColorStruct(&color, kPsychIndexColor, whiteValue ); //index mode will coerce to any other.
		}

 	PsychCoerceColorMode( &color);
        
	//get the x and y position values. 
	PsychAllocInDoubleArg(3, TRUE,  &xPosition);
	PsychAllocInDoubleArg(4, TRUE,  &yPosition);
	dotSize=1;	//set the default
	PsychCopyInDoubleArg(5, FALSE, &dotSize);

	// Enable this windowRecords framebuffer as current drawingtarget:
	PsychSetDrawingTarget(windowRecord);

	// Set default draw shader:
	PsychSetShader(windowRecord, -1);

	PsychUpdateAlphaBlendingFactorLazily(windowRecord);
	PsychSetGLColor(&color, windowRecord);
	glEnable(GL_POINT_SMOOTH);
	glPointSize((float) dotSize);
	GLBEGIN(GL_POINTS);
		GLVERTEX2d(*xPosition, *yPosition);
	GLEND();
	glDisable(GL_POINT_SMOOTH);
	glPointSize(1);

	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);

 	//All psychfunctions require this.
	return(PsychError_none);
}
