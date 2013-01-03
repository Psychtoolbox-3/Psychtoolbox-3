/*
PsychAlphaBlending.c

PLATFORMS:	

	Only OS X for now.  
			
AUTHORS:
	
	Allen Ingling		awi		Allen.Ingling@nyu.edu

HISTORY:
	
	1/7/04		awi		Wrote it 
						
DESCRIPTION:

	Functions for getting and setting the alpha blending rule in windows.  To assure consistent treatement of state ONLY functions in this file should
	be used as a gateway to modify Screen state.
	
*/
	

#include "Screen.h"


#define DEFAULT_SRC_ALPHA_FACTOR		GL_SRC_ALPHA
#define DEFAULT_DST_ALPHA_FACTOR		GL_ONE

#define NUM_BLENDING_MODE_CONSTANTS		11

char *blendingModeStrings[]={
	"GL_ZERO",
	"GL_ONE",
	"GL_SRC_COLOR",
	"GL_ONE_MINUS_SRC_COLOR",
	"GL_DST_COLOR",
	"GL_ONE_MINUS_DST_COLOR",
	"GL_SRC_ALPHA",
	"GL_ONE_MINUS_SRC_ALPHA",
	"GL_DST_ALPHA",
	"GL_ONE_MINUS_DST_ALPHA",
	"GL_SRC_ALPHA_SATURATE"
};

GLenum blendingModeConstants[]={
	GL_ZERO,
	GL_ONE,
	GL_SRC_COLOR,
	GL_ONE_MINUS_SRC_COLOR,
	GL_DST_COLOR,
	GL_ONE_MINUS_DST_COLOR,
	GL_SRC_ALPHA,
	GL_ONE_MINUS_SRC_ALPHA,
	GL_DST_ALPHA,
	GL_ONE_MINUS_DST_ALPHA,
	GL_SRC_ALPHA_SATURATE 
};


/* 
	PsychValidateBlendingConstantForSource()
	
	Constants are specified for both the source and destination.  Not all constants are available for both.  This function checks
	a constnat to see if it may be specified for the source.
	
	Return TRUE if the choice if valid and false otherwise.
	
*/
psych_bool PsychValidateBlendingConstantForSource(GLenum sourceFactor)
{
	psych_bool isValid;

	switch(sourceFactor)
	{
		case GL_ZERO:					isValid=TRUE;  break;
		case GL_ONE:					isValid=TRUE;  break;	
		case GL_SRC_COLOR:				isValid=FALSE;  break;
		case GL_ONE_MINUS_SRC_COLOR:	isValid=FALSE;  break;
		case GL_DST_COLOR:				isValid=TRUE;  break;
		case GL_ONE_MINUS_DST_COLOR:	isValid=TRUE;  break;
		case GL_SRC_ALPHA:				isValid=TRUE;  break;
		case GL_ONE_MINUS_SRC_ALPHA:	isValid=TRUE;  break;
		case GL_DST_ALPHA:				isValid=TRUE;  break;
		case GL_ONE_MINUS_DST_ALPHA:	isValid=TRUE;  break;
		case GL_SRC_ALPHA_SATURATE:		isValid=TRUE;  break; 
		default: PsychErrorExitMsg(PsychError_internal, "Failed to find alpha blending factor constant when validating for source specification");
	}
	return(isValid);
}


/* 
	PsychValidateBlendingConstantForSource()
	
	Constants are specified for both the source and destination.  Not all constants are available for both.  This function checks
	a constant to see if it may be specified for the destination.
	
	Return TRUE if the choice if valid and false otherwise.
  
*/
psych_bool PsychValidateBlendingConstantForDestination(GLenum destinationFactor)
{
	psych_bool isValid;
	
	switch(destinationFactor)
	{
		case GL_ZERO:					isValid=TRUE;  break;
		case GL_ONE:					isValid=TRUE;  break;	
		case GL_SRC_COLOR:				isValid=TRUE;  break;
		case GL_ONE_MINUS_SRC_COLOR:	isValid=TRUE;  break;
		case GL_DST_COLOR:				isValid=FALSE;  break;
		case GL_ONE_MINUS_DST_COLOR:	isValid=FALSE;  break;
		case GL_SRC_ALPHA:				isValid=TRUE;  break;
		case GL_ONE_MINUS_SRC_ALPHA:	isValid=TRUE;  break;
		case GL_DST_ALPHA:				isValid=TRUE;  break;
		case GL_ONE_MINUS_DST_ALPHA:	isValid=TRUE;  break;
		case GL_SRC_ALPHA_SATURATE:		isValid=FALSE;  break; 
		default: PsychErrorExitMsg(PsychError_internal, "Failed to find alpha blending factor constant when validating for destination specification");
	}
	return(isValid);
}




/*  
	PsychGetGLBlendConstantFromString()
	
	Lookup the constant from the string.  Return true if we find it.
*/
int PsychGetAlphaBlendingFactorConstantFromString(char *blendString, GLenum *blendConstant)
{
	int		i;
	
	for(i=0;i<NUM_BLENDING_MODE_CONSTANTS;i++){
		if(PsychMatch(blendingModeStrings[i], blendString)){
			*blendConstant=blendingModeConstants[i];
			return(TRUE);
		}
	}
	return(FALSE);
}


/*
	stringSize	PsychGetGLBlendStringFromConstant(GLenum blendConstant, NULL);
	foundIt		PsychGetGLBlendSTringFromConstnat(GLenum blendConstant, char *blendString)

	Lookup the string from the constant.  The caller allocates the return string.
	to find the length of the return string, first call lookup with NULL for the 
	return.  Call a second time and pass the pointer to your allocated space.  
	
	PsychGetBlendStringFromConstant will always return 0 if it fails to find the 
	constant.
*/
int PsychGetAlphaBlendingFactorStringFromConstant(GLenum blendConstant, char *blendString)
{
	int		i;
	
	for(i=0;i<NUM_BLENDING_MODE_CONSTANTS;i++){
		if(blendConstant==blendingModeConstants[i]){
			if(blendString!=NULL)
				strcpy(blendString, blendingModeStrings[i]);
			return((int) strlen(blendingModeStrings[i]) + 1);
		}
	}
	return(0);
}



/*
	PsychInitWindowRecordAlphaBlendingFields()
	
	Called when we create a new window.
*/
void PsychInitWindowRecordAlphaBlendingFactors(PsychWindowRecordType *winRec)
{
	winRec->actualEnableBlending=FALSE;
	winRec->actualSourceAlphaBlendingFactor=GL_ONE;				
	winRec->actualDestinationAlphaBlendingFactor=GL_ZERO;
	winRec->nextSourceAlphaBlendingFactor=GL_ONE;
	winRec->nextDestinationAlphaBlendingFactor=GL_ZERO;
	
}

/*
	PsychGetBlendConstantsFromWindow()
	
	Get the blending constants which will be used if we draw into the window.
*/
void PsychGetAlphaBlendingFactorsFromWindow(PsychWindowRecordType *winRec, GLenum *oldSource, GLenum *oldDestination)
{
	*oldSource=winRec->nextSourceAlphaBlendingFactor;
	*oldDestination=winRec->nextDestinationAlphaBlendingFactor;
}


/* 
	PsychStoreAlphaBlendingFactorsForWindow()
	
	To avoid unnecessary context switching, only store the desired blending factors in the window record.  Defer calls to glBlendFunc()
	until a drawing function, which must change context, invokes PsychUpdateAlphaBlendingFactorLazily() which in turn calls glBlendFunc()
	the blending mode has changed.  
*/
void PsychStoreAlphaBlendingFactorsForWindow(PsychWindowRecordType *winRec, GLenum sourceBlendConstant, GLenum destinationBlendConstant)
{
	winRec->nextSourceAlphaBlendingFactor=sourceBlendConstant;
	winRec->nextDestinationAlphaBlendingFactor=destinationBlendConstant;
}


/*
	PsychUpdateAlphaBlendingFactorLazily()
	
	Sets both, the OpenGL blending function and mode, and the OpenGL colorbuffer writemask for given
	onscreen- or offscreen window 'winRec', based on cached settings made via Screen('Blendfunction').
	
	Although the function is called xxxLazily() it doesn't do this lazily, but always. Deferred updates
	were the cause of many subtile nasty bugs in the past, so better safe than sorry.
*/
void PsychUpdateAlphaBlendingFactorLazily(PsychWindowRecordType *winRec)
{
	// Enable alpha-blending whenever some blendfunction other than GL_ONE, GL_ZERO is selected, disable otherwise:
	if (winRec->nextSourceAlphaBlendingFactor != GL_ONE || winRec->nextDestinationAlphaBlendingFactor != GL_ZERO) {
		// Alpha blending requested. Possible on this hardware for current drawable?
		if ((winRec->bpc < 16) || ((winRec->bpc == 16) && (winRec->imagingMode & kPsychNeed16BPCFixed)) ||
		    ((winRec->bpc == 16) && (winRec->gfxcaps & kPsychGfxCapFPBlend16)) || ((winRec->bpc == 32) && (winRec->gfxcaps & kPsychGfxCapFPBlend32)) ) {
			// Alpha blending should work - Set it up:
			winRec->actualEnableBlending=TRUE;
			glEnable(GL_BLEND);
			winRec->actualSourceAlphaBlendingFactor=winRec->nextSourceAlphaBlendingFactor;
			winRec->actualDestinationAlphaBlendingFactor = winRec->nextDestinationAlphaBlendingFactor;
			glBlendFunc(winRec->actualSourceAlphaBlendingFactor, winRec->actualDestinationAlphaBlendingFactor);
		}
		else {
			// Alpha blending unsupported for this drawable. Disable it:
			winRec->actualEnableBlending=FALSE;
			glDisable(GL_BLEND);
			glBlendFunc(GL_ONE, GL_ZERO);
		}
	}
	else {
		winRec->actualEnableBlending=FALSE;
		glDisable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ZERO);
	}
	
	// Set GL color writemask:
	glColorMask(winRec->colorMask[0], winRec->colorMask[1], winRec->colorMask[2], winRec->colorMask[3]);
}
