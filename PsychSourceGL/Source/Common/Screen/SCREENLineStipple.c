/*
	SCREENLineStipple.c		
  
	AUTHORS:

		Allen.Ingling@nyu.edu		awi 
  
	PLATFORMS:	
	
		Only OS X for now.
    

	HISTORY:

		03/18/05	awi		Created.  

	TO DO:
  
*/


#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "[stippleEnabled stippleFactor stipleVector]=Screen('LineStipple', windowPtr, [enable], [factor], [stipplePattern] )";                                          
//                                                                                           1           2         3         4      
static char synopsisString[] = 
	"Set 'DrawLine' to draw dashed lines.  Argument \"stippleArray\" is a 16-element logical array specifying the stipple pattern. "
	"Argument \"factor\"  scales the pattern by an integer mutiple.  By default, factor is 1 and stipplePattern is "
	"[0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1]"; 
	
static char seeAlsoString[] = "FrameRect";	

static void ConvertShortToStippleArray(PsychNativeBooleanType *patternArray, GLushort shortPattern);
static void ConvertStippleArrayToShort(psych_bool *patternArray, GLushort *shortPattern);

PsychError SCREENLineStipple(void)  
{
	
	PsychWindowRecordType			*winRec;
	static GLushort					stipplePatternTemp;
	psych_bool							isFactorThere, isPatternThere, isFlagThere, didChange;
	double							*newFactor;
	psych_bool							*newPatternArray;
	int								numInputVectorElements;
	PsychNativeBooleanType			*oldPatternArray;
	psych_bool							*newEnableFlag;
								
    
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(4));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(3));  //The maximum number of outputs

	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(1, kPsychArgRequired, &winRec);
	
	//return existing values
	PsychCopyOutFlagArg(1, kPsychArgOptional, winRec->stippleEnabled);
	PsychCopyOutDoubleArg(2,  kPsychArgOptional, (double)winRec->stippleFactor); 
	PsychAllocOutBooleanMatArg(3, kPsychArgOptional, 1, 16, 0, &oldPatternArray);
	ConvertShortToStippleArray(oldPatternArray, winRec->stipplePattern);
	
	//read in new values
	didChange=FALSE;
	isFlagThere=PsychAllocInFlagArg(2, kPsychArgOptional, &newEnableFlag);
	if(isFlagThere && *newEnableFlag != winRec->stippleEnabled){
		didChange=TRUE;
		winRec->stippleEnabled=*newEnableFlag;
	}
	
	isFactorThere=PsychAllocInDoubleArg(3, kPsychArgOptional, &newFactor);
	if(isFactorThere && (GLint)(*newFactor) != winRec->stippleFactor){
		didChange=TRUE;
		winRec->stippleFactor=(GLint)(*newFactor);
	}
	
	//NOTE: fix PsychAllocInFlagArgVector so that it limits to numInputVectorElements.  
	isPatternThere=PsychAllocInFlagArgVector(4, kPsychArgOptional, &numInputVectorElements, &newPatternArray);
	if(isPatternThere){
		if(numInputVectorElements != 16)
			PsychErrorExitMsg(PsychError_inputMatrixIllegalDimensionSize, "Argument \"stipplePattern\" should be a vector of 16 elements in size");
		ConvertStippleArrayToShort(newPatternArray, &stipplePatternTemp);
		if(stipplePatternTemp != winRec->stipplePattern){
			didChange=TRUE;
			winRec->stipplePattern=stipplePatternTemp;
		}
	}
	
	//Update GL context according to new settings. 
	if(didChange){   //avoids unnecessary context switches
		PsychSetGLContext(winRec);
		glLineStipple(winRec->stippleFactor, winRec->stipplePattern);
		if(winRec->stippleEnabled)
			glEnable(GL_LINE_STIPPLE);
		else
			glDisable(GL_LINE_STIPPLE);
	}
   
 	//All psychfunctions require this.
	return(PsychError_none);
}


static void ConvertShortToStippleArray(PsychNativeBooleanType *patternArray, GLushort shortPattern)
{
	int i;
	
	for(i=0;i<16;i++){
		if(shortPattern & 1)
			patternArray[i]=1;
		else
			patternArray[i]=0;
		shortPattern=shortPattern>>1;
	}
}


static void ConvertStippleArrayToShort(psych_bool *patternArray, GLushort *shortPattern)
{
	int i;
	
	*shortPattern=0;
	for(i=0;i<16;i++){
		*shortPattern=*shortPattern<<1;
		if(patternArray[15-i])
			*shortPattern=*shortPattern | 1;
	}
}
