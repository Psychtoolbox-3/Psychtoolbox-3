/*
	PsychToolbox2/Source/OSX/FontInfo/ScreenFontGlue.h
	
	PLATFORMS:	
	
		OS X  
				
	AUTHORS:
	
    Allen Ingling		awi		Allen.Ingling@nyu.edu
    Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
		7/8/04		awi		Split it of from PsychFontGlue because PsychWindowRecordType in PsychSetATSUStyleAttributes
							was causing trouble for the FontInfo project
		3/7/06		awi		Changed references from "Font" to "FontInfo".  The function had been previously renamed, but not all references updated.
							
	DESCRIPTION:
        

*/

//include once
#ifndef PSYCH_IS_INCLUDED_ScreenFontGlue
#define PSYCH_IS_INCLUDED_ScreenFontGlue

#include "Psych.h"

// Only include this function when building Screen(), but not FontInfo():
#ifdef PTBMODULE_Screen
void PsychSetATSUStyleAttributesFromPsychWindowRecord(ATSUStyle atsuStyle,  PsychWindowRecordType *winRec);
#endif

#endif
