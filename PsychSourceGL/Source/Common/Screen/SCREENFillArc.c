/*
	SCREENFillArc.cpp		
 
	Draw a filled arc, like a part of a pie-chart.
 
	AUTHORS:
 
	Allen.Ingling@nyu.edu                   awi
	mario.kleiner at tuebingen.mpg.de       mk
 
	PLATFORMS:
	
	Should compile and work for all platforms.
 
	NOTES:
 
	HISTORY:
 
	12/11/05        mk		Created. Derived from Allen Inglings gluDisc.   
 
	TO DO:
 
 */


#include "Screen.h"

// The routine that does the real work - see at bottom of file:
void PsychRenderArc(unsigned int mode);

PsychError SCREENDrawArc(void)
{
    // If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
    static char useString[] = "Screen('DrawArc',windowPtr,[color],[rect],startAngle,arcAngle)";
    //                                          1         2       3      4          5
    static char synopsisString[] = 
        "Draw an arc inscribed within the rect. 'color' is the clut index (scalar "
        "or [r g b] triplet) that you want to poke into each pixel; default produces "
        "black with the standard CLUT for this window's pixelSize. Default 'rect' is "
        "entire window. Angles are measured clockwise from vertical.";
    static char seeAlsoString[] = "FrameArc FillArc";	
	
    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    //check for superfluous arguments
    PsychErrorExit(PsychCapNumInputArgs(5));        //The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(3));    //The minimum number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));       //The maximum number of outputs
	
    // Render arc of type 1 - Just the outline.
    PsychRenderArc(1);
    
    return(PsychError_none);
}

PsychError SCREENFrameArc(void)
{
    // If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
    static char useString[] = "Screen('FrameArc',windowPtr,[color],[rect],startAngle,arcAngle, [penWidth], [penHeight], [penMode])";
    //                                          1         2       3      4          5           6           7           8
    static char synopsisString[] = 
        "Draw an arc inscribed within the rect. 'color' is the clut index (scalar "
        "or [r g b] triplet) that you want to poke into each pixel; default produces "
        "black with the standard CLUT for this window's pixelSize. Default 'rect' is "
        "entire window. Angles are measured clockwise from vertical. 'penWidth' and "
        "'penHeight' are the width and height of the pen to use. penWidth must equal "
        "penHeight and the 'penMode' argument is currently ignored.";
    static char seeAlsoString[] = "DrawArc FillArc";	
    
    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    //check for superfluous arguments
    PsychErrorExit(PsychCapNumInputArgs(8));        //The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(3));    //The minimum number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));       //The maximum number of outputs
    
    // Render arc of type 2 - Just the outline with penWidth.
    PsychRenderArc(2);
    
    return(PsychError_none);
}


PsychError SCREENFillArc(void)  
{	
	// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
	static char useString[] = "Screen('FillArc',windowPtr,[color],[rect],startAngle,arcAngle)";
	//                                          1         2       3      4          5
	static char synopsisString[] = 
		"Draw a filled arc inscribed within the rect. 'color' is the clut index (scalar "
		"or [r g b a] triplet) that you want to poke into each pixel; default produces "
        "black with the standard CLUT for this window's pixelSize. Default 'rect' is "
		"entire window. Angles are measured clockwise from vertical.";
	static char seeAlsoString[] = "DrawArc FrameArc";	
	
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(5));        //The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(3));    //The minimum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));       //The maximum number of outputs
	
	// Render arc of type 3 - Filled arc.
	PsychRenderArc(3);
	
	return(PsychError_none);
}

void PsychRenderArc(unsigned int mode)
{
	PsychColorType			color;
	PsychRectType           rect;
	double					*startAngle, *arcAngle, *penWidth, *penHeight;
	PsychWindowRecordType	*windowRecord;
	double					whiteValue;
	double                  dotSize;
	psych_bool				isArgThere;
	GLUquadric              *diskQuadric = NULL;
	double cx, cy, w, h;
	
	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
	
	//Get the color argument or use the default, then coerce to the form determened by the window depth.  
	isArgThere=PsychCopyInColorArg(kPsychUseDefaultArgPosition, FALSE, &color);
	if(!isArgThere){
		whiteValue=PsychGetWhiteValueFromWindow(windowRecord);
		PsychLoadColorStruct(&color, kPsychIndexColor, whiteValue ); //index mode will coerce to any other.
	}
 	PsychCoerceColorMode( &color);
	
	// Get the rect to which the object should be inscribed: Default is "full screen"
	PsychMakeRect(rect, 0, 0, PsychGetWidthFromRect(windowRecord->clientrect), PsychGetHeightFromRect(windowRecord->clientrect));
	PsychCopyInRectArg(3, FALSE, rect);
	if (IsPsychRectEmpty(rect)) return;
	
	w=PsychGetWidthFromRect(rect);
	h=PsychGetHeightFromRect(rect);
	
	PsychGetCenterFromRectAbsolute(rect, &cx, &cy);
	if (w==0 || h==0) PsychErrorExitMsg(PsychError_user, "Invalid rect (width or height equals zero) provided!");
	
	// Get start angle: 
	PsychAllocInDoubleArg(4, TRUE,  &startAngle);
	PsychAllocInDoubleArg(5, TRUE,  &arcAngle);
	
	if (mode==2) {
		// Get pen width and height:
		penWidth=NULL;
		penHeight=NULL;
		PsychAllocInDoubleArg(6, FALSE,  &penWidth);
		PsychAllocInDoubleArg(7, FALSE,  &penHeight);
		// Check if penWidth and penHeight spec'd. If so, they
		// need to be equal:
		if (penWidth && penHeight && (*penWidth!=*penHeight)) {
			PsychErrorExitMsg(PsychError_user, "penWidth and penHeight must be equal on OS-X if both are specified!");
		}
		dotSize=1;
		if (penWidth) dotSize = *penWidth;
		if (penHeight) dotSize = *penHeight;
	}
	
	// Enable this windowRecords framebuffer as current drawingtarget:
	PsychSetDrawingTarget(windowRecord);
	
	// Set default drawshader:
	PsychSetShader(windowRecord, -1);
	
	PsychUpdateAlphaBlendingFactorLazily(windowRecord);
	PsychSetGLColor(&color,  windowRecord);
	
    if (PsychIsGLClassic(windowRecord)) {
        // Backup our modelview matrix:
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
	
        // Position disk at center of rect:
        glTranslated(cx, cy, 0);
	
        // Scale in order to fit to rect in case w!=h:
        glScaled(1.0, -h/w, 1.0);
	
        // Draw filled partial disk:
        diskQuadric=gluNewQuadric();
	
        switch (mode) {
		case 1: // One pixel thin arc: InnerRadius = OuterRadius - 1
			gluPartialDisk(diskQuadric, (w/2) - 1.0, w/2, (int) w, 2, *startAngle, *arcAngle);
			break;
		case 2: // dotSize thick arc:  InnerRadius = OuterRadius - dotsize
			gluPartialDisk(diskQuadric, (dotSize < (w/2)) ? (w/2) - dotSize : 0, w/2, (int) w, 2, *startAngle, *arcAngle);
			break;
		case 3: // Filled arc:
			gluPartialDisk(diskQuadric, 0, w/2, (int) w, 1, *startAngle, *arcAngle);
			break;
        }
	
        gluDeleteQuadric(diskQuadric);
	
        // Restore old matrix:
        glPopMatrix();
	}
    else {
        switch (mode) {
		case 1: // One pixel thin arc: InnerRadius = OuterRadius - 1
            PsychDrawDisc(windowRecord, (float) cx, (float) cy, (w/2) - 1.0, w/2, (int) w, (float) 1, (float) -h/w, (float) *startAngle, (float) *arcAngle);
			break;
		case 2: // dotSize thick arc:  InnerRadius = OuterRadius - dotsize
            PsychDrawDisc(windowRecord, (float) cx, (float) cy, (dotSize < (w/2)) ? (w/2) - dotSize : 0, w/2, (int) w, (float) 1, (float) -h/w, (float) *startAngle, (float) *arcAngle);
			break;
		case 3: // Filled arc:
            PsychDrawDisc(windowRecord, (float) cx, (float) cy, 0, w/2, (int) w, (float) 1, (float) -h/w, (float) *startAngle, (float) *arcAngle);
			break;
        }
    }

	// Mark end of drawing op. This is needed for single buffered drawing:
	PsychFlushGL(windowRecord);
	
	return;
}
