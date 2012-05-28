/*
	PsychToolbox2/Source/OSX/FontInfo/FONTSFonts.c		

	PROJECTS: 
	
		FontInfo only.  

	AUTHORS:
	
        Allen Ingling		awi		Allen.Ingling@nyu.edu
        Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de

	PLATFORMS:	
	
		Mac OS X Only

	HISTORY:

		1/07/04		awi		Created. Derived from my earlier ATS incarnation.
		3/7/06		awi		Changed references from "Font" to "FontInfo".  The function had been previously renamed, but not all references updated. 		

	DESCRIPTION: 

	NOTES:

		It's not clear with what to fill the "Name" field of the returned struct because a font has so many names.  Ultimately we want the names
		to match on all platforms, so perhaps we should wait to see which names are available in Windows and then choose one wich is available on both
		platforms.  For now we use the Font Manager name.
  

*/ 

#include "FontInfo.h"


static char useString[] = "fontList=FontInfo('Fonts')";
static char synopsisString[] = 
	"Returns a list of structs, one struct per available font, holding "
    "information about each font.";
static char seeAlsoString[] = "NumFonts";

PsychError FONTSFonts(void) 
{
    PsychGenericScriptType 	*nativeStructArray;
    int				arrayIndex, numFonts;
    PsychFontStructType		**fontPointerList, *fontElement;
    
    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check for  required and superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(0));
        
    //build a list of pointers to all the the font records in the font list, then hand the list of pointers to PsychCopyFontRecordsToNativeStructArray()
    //to get a native struct array of font records and return it to the scripting environment.
    arrayIndex=0;
    numFonts=PsychGetFontListLength();
    fontPointerList=(PsychFontStructType**)malloc(numFonts * sizeof(PsychFontStructType*));
    for(fontElement=PsychGetFontListHead();fontElement;fontElement=fontElement->next)
        fontPointerList[arrayIndex++]=fontElement;
    PsychCopyFontRecordsToNativeStructArray(numFonts, fontPointerList, &nativeStructArray);  
    free((void*)fontPointerList);
    PsychAssignOutStructArray(1, FALSE, nativeStructArray);
    
    return(PsychError_none);	
}
