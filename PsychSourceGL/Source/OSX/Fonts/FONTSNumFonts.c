/*
	PsychToolbox2/Source/OSX/FontInfo/FONTSNumFonts.c		

	PROJECTS: 
	
		FontInfo only.  

	AUTHORS:
	
		Allen.Ingling@nyu.edu		awi 

	PLATFORMS:
	
		Mac OS X Only

	HISTORY:
	
		10/23/03  	awi		Created.
		3/7/06		awi		Changed references from "Font" to "FontInfo".  The function had been previously renamed, but not all references updated.		 

	DESCRIPTION: 

	  
*/ 

#include "FontInfo.h"


static char useString[] = "numFonts=FontInfo('NumFonts')";
static char synopsisString[] = 
	"Returns the number of available fonts on this system.";
static char seeAlsoString[] = "FontFamilies";
	 

PsychError FONTSNumFonts(void) 
{
    double		*numFonts;
	
    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(0));
    
    PsychAllocOutDoubleArg(1, FALSE,  &numFonts);
    *numFonts=(double)PsychGetFontListLength();
        
    return(PsychError_none);	
}
