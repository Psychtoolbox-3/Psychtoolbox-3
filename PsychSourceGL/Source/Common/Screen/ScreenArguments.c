/*
	PsychToolbox2/Source/Common/Screen/ScreenArguments.c
	
	PLATFORMS:	Windows
				MacOS9
			
	
	AUTHORS:
	Allen Ingling		awi		Allen.Ingling@nyu.edu

	HISTORY:
	09/09/02			awi		wrote it.  
        10/20/03			awi		fixed a bug in PsychCopyInColorArg where we failed to save alph value correctly in the color struct.
	11/14/06                        mk              Colors are now also accepted in uint8 format. Empty rects are no longer rejected as invalid.
	DESCRIPTION:
	
	Functions which get and put arguments used by Screen.  	
        
        TO DO
        
	Return the Window index.  Before returning it check to be sure that an argument was supplied
	and that it refers to a currently opened window.  If either is not true then exit with an 
	error.   
*/

#include "Screen.h"


psych_bool PsychCopyInWindowIndexArg(int position, psych_bool required, PsychWindowIndexType *windowIndex)
{

	double *arg;
        psych_bool isThere;
	
	if(position==kPsychUseDefaultArgPosition)
		position = kPsychDefaultNumdexArgPosition; 
        isThere=PsychAllocInDoubleArg(position,required,&arg);
        if(!isThere)
            return(FALSE);
	*windowIndex = (PsychWindowIndexType)*arg;
	if(IsWindowIndex(*windowIndex))
		return(TRUE);
	else{
		PsychErrorExitMsg(PsychError_invalidWindex,NULL);
		return(FALSE);  //only to satisfy the compiler with a return statement.
	}
}


/* 
        PsychAllocInWindowRecordArg()
        
        Lookup the window record using the window index argument. 
        
        Because this looks up the window record using its somewhat different
        from other "PsychAllocIn*Arg" functions, but we use the same naming scheme nonetheless.
        From the caller's point of view looking up a prexisting record and allocating a new
        one should be equivalent.  
*/
psych_bool PsychAllocInWindowRecordArg(int position, psych_bool required, PsychWindowRecordType **winRec)
{
	PsychWindowIndexType windowIndex;
	double arg;
	psych_bool isThere;
        
	if(position==kPsychUseDefaultArgPosition)
            position = kPsychDefaultNumdexArgPosition; 	
	isThere=PsychCopyInDoubleArg(position,required,&arg);
        if(!isThere)
            return(FALSE);
	windowIndex = (PsychWindowIndexType)arg;
	PsychErrorExit(FindWindowRecord(windowIndex, winRec));
        return(TRUE);
}


/*
    PsychCopyInScreenNumberArg()
    
    Automaticially derive the screen number from a window index if provided.
    Otherwise return the screen number if that is what is provided.
        
*/
psych_bool PsychCopyInScreenNumberArg(int position, psych_bool required, int *screenNumber)
{
	PsychNumdexType numdex;
	PsychWindowRecordType *winRec;
	double arg;
        psych_bool isThere;

	if(position==kPsychUseDefaultArgPosition)
		position = kPsychDefaultNumdexArgPosition; 	
	isThere=PsychCopyInDoubleArg(position,required,&arg);
        if(!isThere)
            return(FALSE);
	numdex = (PsychNumdexType)arg;
	if(IsWindowIndex(numdex)){ 
		//it's a window index, so get the window record and from that get the screen number.  
		FindWindowRecord((PsychWindowIndexType)numdex, &winRec);
		*screenNumber=winRec->screenNumber;
                return(TRUE);
	}else if(IsValidScreenNumber(numdex)){
		//it's a screen number, so just return it.
		*screenNumber=(int)numdex;
                return(TRUE);
	}else{
		//we were passed something that is neither a window index nor a screen number so issue an error.
		PsychErrorExitMsg(PsychError_invalidNumdex,NULL);
		return(FALSE);
	}
}


/*
    PsychCopyInScreenNumberOrUnaffiliatedArg()
    
    Like PsychCopyInScreenNumberArg but we will also accept -1 to indicate no particular screen.
        
*/
psych_bool PsychCopyInScreenNumberOrUnaffiliatedArg(int position, psych_bool required, int *screenNumber)
{
	PsychNumdexType numdex;
	PsychWindowRecordType *winRec;
	double arg;
        psych_bool isThere;

	if(position==kPsychUseDefaultArgPosition)
		position = kPsychDefaultNumdexArgPosition; 	
	isThere=PsychCopyInDoubleArg(position,required,&arg);
        if(!isThere)
            return(FALSE);
	numdex = (PsychNumdexType)arg;
	if(IsWindowIndex(numdex)){ 
		//it's a window index, so get the window record and from that get the screen number.  
		FindWindowRecord((PsychWindowIndexType)numdex, &winRec);
		*screenNumber=winRec->screenNumber;
                return(TRUE);
	}else if(IsValidScreenNumberOrUnaffiliated(numdex)){
		//it's a screen number, so just return it.
		*screenNumber=(int)numdex;
                return(TRUE);
	}else{
		//we were passed something that is neither a window index nor the token meaning unaffiliated window nor a screen number so issue an error.
		PsychErrorExitMsg(PsychError_invalidNumdex,NULL);
		return(FALSE);
	}
}


/*
    PsychAllocInScreenNumberArg()
    
    This automaticially derives the screen number from a window index if that is what is passed.
    
    Write PsychCopyInRawScreenNumberArg for cases where we might not want to do that. 
*/
psych_bool PsychAllocInScreenRecordArg(int position, psych_bool required, PsychScreenRecordType **screenRecord)
{
	int screenNumber;
        psych_bool isThere;
	
	if(position==kPsychUseDefaultArgPosition)
		position = kPsychDefaultNumdexArgPosition; 
        isThere=PsychCopyInScreenNumberArg(position, required, &screenNumber);
        if(!isThere)
            return(FALSE);
	PsychErrorExitMsg(FindScreenRecord(screenNumber, screenRecord),NULL);
	return(TRUE);	
}

		
psych_bool PsychIsWindowIndexArg(int position)
{
	PsychNumdexType numdex;
	double arg;

	if(position==kPsychUseDefaultArgPosition)
		position =kPsychDefaultNumdexArgPosition; 	
	if(!PsychCopyInDoubleArg(position,FALSE,&arg))
            return(FALSE);
	numdex = (PsychNumdexType)arg;
	return(IsWindowIndex(numdex));
	
}

/*
	PsychIsScreenNumberOrUnaffiliatedArg()
	
	Returns true iff the argument in the specified position is either a window argumen or one of the key values
	which means not to associate a window with any particular display.  For now that value is -1 but  within
	the PTB environment we could also alow NaN because that makes more sense than -1.
*/
psych_bool PsychIsScreenNumberOrUnaffiliatedArg(int position)
{
	PsychNumdexType numdex;
	double arg;

	if(position==kPsychUseDefaultArgPosition)
		position = kPsychDefaultNumdexArgPosition; 	
	if(!PsychCopyInDoubleArg(position,FALSE,&arg))
            return(FALSE);
	numdex = (PsychNumdexType)arg;
	return(IsValidScreenNumber(numdex) || numdex==kPsychUnaffiliatedWindow);
	
}

/*
	PsychIsUnaffiliatedScreenNumberArg()
	
	Returns true iff the argument in the specified position indicatesnot to associate a window with any particular
	display.  For now that value is -1 but  within the PTB environment we could also alow NaN because that makes 
	more sense than -1.
*/
psych_bool PsychIsUnaffiliatedScreenNumberArg(int position)
{
	PsychNumdexType numdex;
	double arg;

	if(position==kPsychUseDefaultArgPosition)
		position = kPsychDefaultNumdexArgPosition; 	
	if(!PsychCopyInDoubleArg(position,FALSE,&arg))
            return(FALSE);
	numdex = (PsychNumdexType)arg;  //cast is to int.  We should get rid of this obfuscatory numdex type. 
	return(numdex==kPsychUnaffiliatedWindow);
	
}



psych_bool PsychIsScreenNumberArg(int position)
{
	PsychNumdexType numdex;
	double arg;

	if(position==kPsychUseDefaultArgPosition)
		position = kPsychDefaultNumdexArgPosition; 	
	if(!PsychCopyInDoubleArg(position,kPsychArgAnything,&arg))
            return(FALSE);
	numdex = (PsychNumdexType)arg;
	return(IsValidScreenNumber(numdex));
	
}



/* put and get color and rect specifiers directly from arguments

-Naming scheme:
	-Outputing return arguments:
		- "PsychAllocOut*Arg" : allocate and set a pointer to volatile memory to be filled with returned information by the caller.
		- "PsychCopyOut*Arg : accept a pointer to ouput values and fill in the return matrix memory with the values.
	-Reading input arguments:
		- "PsychAllocIn*Arg" : set a pointer to volatile memory allocated by "PsychAllocIn*Arg" and holding the input value.
		- "PsychCopyIn*Arg" : accept a pointer to memory and fill in that memory with the input argument values.     
		
*/



/*
	Accept a pointer to a Psychtoolbox color specifier and fill it with the 
	color information supplied at the specified argument position in module call .
	
	The behavior depends on the value of required:
	
	1. Required = TRUE
		A. If argument is not present exit with an error PsychError_invalidColorArg.
		B. If argument is present and valid then load it into *color and return true.
	2. Required = FALSE 
		A. If argument is not present then don't touch *color and return false.
		B. If argument is present and valid then load it into *color and return true.
*/
psych_bool PsychCopyInColorArg(int position, psych_bool required, PsychColorType *color)
{
	int i,m,n,p,argSize;
	psych_bool isArg;
	double dummyColor[4];
	double *colorArgMat=NULL;
	unsigned char *colorArgMatBytes=NULL;
	
	if(position == kPsychUseDefaultArgPosition)
		position = kPsychDefaultColorArgPosition;
	isArg = PsychIsArgPresent(PsychArgIn, position);
	if(!isArg){
		if(required)
			PsychErrorExitMsg(PsychError_user, "No color argument supplied"); //1A
		else
			return(FALSE);	//2A
	}

	// Try to retrieve double-matrix:
	if (!PsychAllocInDoubleMatArg(position, kPsychArgAnything, &m, &n, &p, &colorArgMat)) {
	  // No double matrix: Try to retrieve uint8 matrix:
	  if (!PsychAllocInUnsignedByteMatArg(position, TRUE, &m, &n, &p, &colorArgMatBytes)) {
	    PsychErrorExitMsg(PsychError_user, "No color argument or invalid color argument supplied");
	  }

	  // Color as uint8 received: Convert to double.
	  if(p!=1) PsychErrorExit(PsychError_invalidColorArg);
	  argSize = m*n;
	  for(i=0; i<argSize; i++) dummyColor[i] = (double) colorArgMatBytes[i];
	  colorArgMat = (double*) (&dummyColor[0]);
	}

	if(p!=1) PsychErrorExit(PsychError_invalidColorArg);
	argSize = m*n;

	if(argSize==4){
		color->mode = kPsychRGBAColor;
		color->value.rgba.r = colorArgMat[0];
		color->value.rgba.g = colorArgMat[1];
		color->value.rgba.b = colorArgMat[2];
		color->value.rgba.a = colorArgMat[3];
                return(TRUE);	//1B, 2B	
	}if(argSize==3){
		color->mode = kPsychRGBColor;
		color->value.rgb.r = colorArgMat[0];
		color->value.rgb.g = colorArgMat[1];
		color->value.rgb.b = colorArgMat[2];
		return(TRUE);	//1B, 2B	
	}else if(argSize==1){
		color->mode = kPsychIndexColor;
		color->value.index.i = colorArgMat[0];
		return(TRUE); //1B, 2B	
	}else{ 
		PsychErrorExit(PsychError_invalidColorArg);
		return(FALSE);
	}
}


		


/*
	Accept a pointer to a Psychtoolbox color specifier and return it from the module
	call in the specified argument position. 
	
	The return value indicates whether the return argument was present.   
	
*/


psych_bool PsychCopyOutColorArg(int position, psych_bool required, PsychColorType *color)
{
	double *colorArgMat;
	
	if(color->mode == kPsychRGBAColor){
		if(!PsychAllocOutDoubleMatArg(position, required, 1, 4, 0, &colorArgMat))
			return(FALSE);
		colorArgMat[0] = (double)color->value.rgba.r;
		colorArgMat[1] = (double)color->value.rgba.g;
		colorArgMat[2] = (double)color->value.rgba.b;
		colorArgMat[3] = (double)color->value.rgba.a;
		return(TRUE);
	}else if(color->mode == kPsychRGBColor){
		if(!PsychAllocOutDoubleMatArg(position, required, 1, 3, 0, &colorArgMat))
			return(FALSE);
		colorArgMat[0] = (double)color->value.rgb.r;
		colorArgMat[1] = (double)color->value.rgb.g;
		colorArgMat[2] = (double)color->value.rgb.b;
		return(TRUE);
        }else if(color->mode == kPsychIndexColor){
		if(!PsychAllocOutDoubleMatArg(position, required, 1, 1, 0, &colorArgMat))
			return(FALSE);
		colorArgMat[0] = (double)color->value.index.i;
		return(TRUE);
        }
        PsychErrorExitMsg(PsychError_internal, "Unrecognized color mode");
        return(FALSE);  //make the compiler happy
}


/*
	Accept a pointer to a Psychtoolbox rect specifier and fill it with the rect
	information supplied at the specified argument position in module call.
	
	The behavior depends on the value of required:
	
	1. Required = TRUE
		A. If argument is not present exit with an error PsychError_invalidColorArg.
		B. If argument is present and valid then load it into *rect and return true.
	2. Required = FALSE 
		A. If argument is not present then don't touch *color and return false.
		B. If argument is present and valid then load it into *rect and return true.

*/
psych_bool PsychCopyInRectArg(int position, psych_bool required, PsychRectType rect)
{
	int m,n,p,argSize;
	psych_bool isArg;
	double *rectArgMat=NULL;
	
	if(position == kPsychUseDefaultArgPosition)
		position = kPsychDefaultRectArgPosition;
	isArg = PsychIsArgPresent(PsychArgIn, position);
	if(!isArg){
		if(required)
			PsychErrorExitMsg(PsychError_user, "Required rect argument missing."); //1A
		else
			return(FALSE);	//2A
	}
	PsychAllocInDoubleMatArg(position, TRUE, &m, &n, &p, &rectArgMat);
	if(p!=1)
		PsychErrorExitMsg(PsychError_invalidRectArg, ">2 D array passed as rect");
	argSize = m*n;
	if(argSize!=4)
		PsychErrorExitMsg(PsychError_invalidRectArg, "rect argument not 4 elements in size");
	memcpy(rect,rectArgMat,sizeof(PsychRectType));
	if(!ValidatePsychRect(rect)){
		PsychErrorExitMsg(PsychError_invalidRectArg, "impossible rect values given");
		return(FALSE);
	}		
	return(TRUE); //1B, 2B
}

/*
	Accept a pointer to a Psychtoolbox rect specifier and return it from the module
	call in the specified argument position. 
	
	The return value indicates whether the return argument was present. 
	
	We could add a PsychAllocRectArg statement also.    
	
*/
psych_bool PsychCopyOutRectArg(int position, psych_bool required, PsychRectType rect)
{
	double *rectArgMat;
	
	if(!PsychAllocOutDoubleMatArg(position, required, 1, 4, 0, &rectArgMat))
		return(FALSE); //optional argument was omitted
	memcpy(rectArgMat,rect,sizeof(PsychRectType));
	return(TRUE);
}

psych_bool PsychAllocOutRectArg(int position, psych_bool required, const double **rect)
{
	double *rectArg;
	
	if(!PsychAllocOutDoubleMatArg(position, required, 1, 4, 0, &rectArg))
		return(FALSE); //optional argument was omitted
	*rect=rectArg;
	return(TRUE);
}

/*
    PsychCopyInSingleDepthArg()
    
    The depth argument must already be allocated and initialized. For reading a 1x1 matrix holding a depth.  
    Use PsychCopyInMultiDepthArg() for reading in lists of depths.
*/

psych_bool PsychCopyInSingleDepthArg(int position, psych_bool required, PsychDepthType *depth)
{
    double *depthArray;
    int m,n,p;
    
    
    position= kPsychUseDefaultArgPosition ?  kPsychDefaultDepthArgPosition : position;
    if(!PsychAllocInDoubleMatArg(position, required, &m, &n, &p, &depthArray))  //shouldn't this check the dimensions of the matrix ? 
        return(FALSE);
    if(m!=1 || n!=1 || p!=1)
        PsychErrorExitMsg(PsychError_invalidDepthArg, "1x1 depth argument expected");
    PsychAddValueToDepthStruct((int) depthArray[0], depth); 
    return(TRUE);
}


psych_bool	PsychCopyInDepthValueArg(int position, psych_bool required, int *depth)
{
    psych_bool	isThere;
    double	value;
    
    if(isThere=PsychCopyInDoubleArg(position,required, &value)){
        *depth=(int)value; 
        if(!(*depth==8 || *depth==16 || *depth==24 || *depth == 32))
            PsychErrorExitMsg(PsychError_invalidDepthArg, "Illegal depth value");
    }
    return(isThere);

}



psych_bool PsychCopyOutDepthArg(int position, psych_bool required, PsychDepthType *depths)
{
    double *depthsArray;
    int i;
    
    if(!PsychAllocOutDoubleMatArg(position, required, 1, PsychGetNumDepthsFromStruct(depths), 0, &depthsArray))
        return(FALSE); //optional argument was omitted 
    for(i=0;i<PsychGetNumDepthsFromStruct(depths);i++)
        depthsArray[i]=(double)PsychGetValueFromDepthStruct(i,depths);
    return(TRUE);
}
