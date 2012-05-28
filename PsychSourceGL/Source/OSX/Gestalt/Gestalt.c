/*
	Gestalt.c		

	PROJECTS:
	
		Gestalt only

	PLATFORMS:  
	
		Only OS X

	AUTHORS:
	
		Allen.Ingling@nyu.edu           awi
        mario.kleiner@tuebingen.mpg.de  mk

	DESCRIPTION:
	
		Gestalt is a work-alike replacement for the Gestalt function build into MATLAB 5.  Mathworks omitted it from MATLAB 6.  For backwards  compatability
		with some OS 9 MATLAB scripts, we supply it in the OS X Psychtoolbox.  The MATLAB 5 help file for Gestalt:
		
			 GESTALT Macintosh Gestalt function.
				GESTALTBITS = GESTALT('SELECTOR') passes the four-character string
				SELECTOR to the Macintosh Operating System function Gestalt.  For
				details about Gestalt, refer to Chapter 1 of Inside Macintosh:
				Operating System Utilities. 
			 
				The result, a 32-bit integer, is stored bitwise in GESTALTBITS.  Thus,
				the least significant bit of the result is GESTALTBITS(32), while the
				most significant bit is GESTALTBITS(1).
			 
				GESTALT is only available on the Macintosh.
			 
				Example:
				   After executing the following line:
			 
					  gestaltbits = gestalt('sysa')
			 
				   gestaltbits(32) will be 1 if run from a 680x0-based Macintosh, while
				   gestaltbits(31) will be 1 if run from a PowerPC-based Macintosh. 
		
	HISTORY:
		
		12/6/04		awi		Created.
		12/7/04		awi		Reteurn the error code instead of exiting with an error. 
		
  
  
*/


#include "Gestalt.h"


//
PsychError GESTALTGestalt(void) 
{
	char					*selectorCString;
	OSType					selectorConstant;
	CFStringRef				selectorCFString;
	OSErr					callError;
	SInt32					responseLong, i;
	PsychNativeBooleanType	*responseArray;
	psych_bool				returnResponse;

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));
	
    PsychAllocInCharArg(1, kPsychArgRequired, &selectorCString);
	if( strlen(selectorCString) != 4)
		PsychErrorExitMsg(PsychError_user, "The selector code must be a four-character string");
	selectorCFString=CFStringCreateWithCString (NULL, selectorCString, kCFStringEncodingUTF8);
	selectorConstant= UTGetOSTypeFromString(selectorCFString); 
	callError= Gestalt(selectorConstant, &responseLong);
	CFRelease(selectorCFString);
	if(callError)
		PsychCopyOutDoubleArg(1, kPsychArgOptional, (double)callError);
	else{
		returnResponse=PsychAllocOutBooleanMatArg(1, kPsychArgOptional, 1, 32, 1, &responseArray);
		if(returnResponse){
			for(i=0;i<32;i++){
				if(1<<i & responseLong) 
					responseArray[31-i]=TRUE;
				else
					responseArray[31-i]=FALSE;
			}
		}
	}
    return(PsychError_none);	
}




