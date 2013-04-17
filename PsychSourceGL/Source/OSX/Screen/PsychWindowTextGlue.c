/*
	PsychToolbox3/Source/OSX/Screen/PsychWindowTextGlue.c
	
	PLATFORMS:	
	
		This is the OS X Core Graphics version.  
				
	AUTHORS:

		Allen Ingling		awi		Allen.Ingling@nyu.edu
        Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de

	HISTORY:
		
		11/18/03		awi		Wrote it 
							
	DESCRIPTION:
        
	NOTES: 
        
	TO DO:
	
		Some of this stuff looks platform neutral so we should consider moving that part of it outside of the the glue layer which is
        only for functions which abstract up platform-specific calls into platform neutral calls invoked by Screen subcommands.
        

*/

#include "Screen.h"


const char *PsychTextDrawingModeNames[]={"TextFill", "TextStroke", "TextFillStroke", "TextInvisible", 
                                "FillClip", "StrokeClip", "TextFillStrokeClip", "TextClip" };
                                

const PsychTextDrawingModeType PsychTextDrawingModes[]= {kPsychTextFill, kPsychTextStroke, kPsychTextFillStroke, kPsychTextInvisible,
                                                                    kPsychTextFillClip, kPsychTextStrokeClip, kPsychTextFillStrokeClip, kPsychTextClip};


void PsychInitTextRecordSettings(PsychTextAttributes *settings)
{
	const char*	tryFontName;
	psych_bool	foundFont;
	PsychFontStructType	*initFontRecord;
	PsychPrefStateGet_DefaultFontName(&tryFontName);

	settings->DisplayList=0;
    settings->textMode=kPsychTextFill;
    settings->textPositionX=0;
    settings->textPositionY=0;
	settings->textSize= PsychPrefStateGet_DefaultTextSize();
	settings->textStyle= PsychPrefStateGet_DefaultTextStyle();	// 0=normal,1=bold,2=italic,4=underline,8=outline,32=condense,64=extend

	// Set default font number and name:
    strcpy((char*) settings->textFontName, "");
	settings->textFontNumber = 0;
    
	/* to initialize the font record to coherent settings, we choose a default font and lookup the matching number */
	foundFont=PsychGetFontRecordFromFontFamilyNameAndFontStyle((char*) tryFontName, settings->textStyle, &initFontRecord);
	if(foundFont) {
        // Use detected font name and number:
        strcpy((char*) settings->textFontName, tryFontName);
        settings->textFontNumber=initFontRecord->fontNumber;
    }
    else {
		PsychErrorExitMsg(PsychError_internal,"Failed to initialize the window record because the default font for DrawText, Geneva, was not found.");
    }
	
	PsychLoadColorStruct(&(settings->textColor), kPsychIndexColor,  0.0);  //index type which may be coerced into anything.
	PsychLoadColorStruct(&(settings->textBackgroundColor), kPsychRGBAColor, 0.0, 0.0, 0.0, 0.0); // Assign black with zero alpha -- transparent.
	PsychCoerceColorMode(&(settings->textColor));
	PsychCoerceColorMode(&(settings->textBackgroundColor));
}


/* 
    PsychGetTextDrawingModeNameFromTextDrawingModeConstant()
*/
void PsychGetTextDrawingModeNameFromTextDrawingModeConstant(char *modeNameStr, int modeNameStrSize, PsychTextDrawingModeType mode)
{
    int i;
    for(i=0; i<kPsychNumTextDrawingModes; i++){
        if(mode==PsychTextDrawingModes[i]){
            strncpy(modeNameStr, PsychTextDrawingModeNames[i], 255); 
            modeNameStr[255]='\0';
            return;
        }
    }
    PsychErrorExitMsg(PsychError_internal, "Unrecognized text drawing mode constant");
}

/*
    PsychGetTextDrawingModeConstantFromTextDrawingModeName()
    
    Convert a string specifying a drawing mode into a constant specifing a drawing mode.  If we fail to find
    the string in the list of allowable strings then return true to indicate error.  Otherwise return FALSE to 
    indicate no error. 
*/
psych_bool PsychGetTextDrawingModeConstantFromTextDrawingModeName(PsychTextDrawingModeType *mode, char *modeNameStr)
{
    int i;
    
    for(i=0; i<kPsychNumTextDrawingModes; i++){
        if(!strcmp(modeNameStr, PsychTextDrawingModeNames[i])) {
            *mode=PsychTextDrawingModes[i];
            return(FALSE);
        } //if
    } //for
    return(TRUE);
}
