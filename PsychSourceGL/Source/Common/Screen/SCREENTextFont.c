/*
	SCREENTextFont.c	

	AUTHORS:
	
		Allen.Ingling@nyu.edu			awi
		mario.kleiner@tuebingen.mpg.de	mk
		
	PLATFORMS:	
	
		All.

	HISTORY:

		11/19/03	awi		Wrote it.

	DESCRIPTION:

		Sets the text font for the specified window record.

	NOTES:

*/

#include "Screen.h"
// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] ="[oldFontName,oldFontNumber]=Screen('TextFont', windowPtr [,fontNameOrNumber]);";
//                                                           0          1           2
static char synopsisString[] = 
    "Get/Set the font for future text draws in this window.\n"
	" You can ask what the current font is, or specify the desired font by number "
	"or by name (e.g. 'Helvetica'). Font numbers are not consistent from Mac to Mac, "
	"and they aren't supported but silently ignored on MS-Windows and Linux, so use "
	"font names for reliability and portability. Font numbers are only available for "
	"backward compatibility to old OS-9 Psychtoolbox versions.\n"
	"The font name can be a string of at most 255 characters length, e.g. 'Helvetica', "
	"or a list containing one string of at most 255 characters, e.g. {'Helvetica'}. "
	"The default font depends on the operating system and is selected for good readability. "
	"You can query and change it via a call to Screen('Preference', 'DefaultFontName').\n"
	"It's ok to request a non-existent font on OS/X; this will have no effect. If you care, call "
    "TextFont again to find out whether you got the font you requested. See FontDemo.\n"
    "On Linux you can either provide a font name - PTB will select the closest matching available "
    "font for that name, size and style requirements - or you can start the fontName with a dash '-' "
    "followed by a full FontConfig font specifier string which encodes all kinds of properties. The "
    "'fc-list' command under Linux allows you to query all available fonts.\n"
    "Depending on the selected text renderer, Linux can be picky about the supplied fonts - "
    "if you request a non-existent font with the old text renderer, the DrawText command will fail "
    "with an error message. However, with the default FTGL text renderer on Linux, Linux will be "
    "lenient in its font selection. If you call Screen('TextFont') after you've drawn some text, the "
    "command will return the font family name of the true selected font. ";

static char seeAlsoString[] = "";

PsychError SCREENTextFont(void) 
{
    psych_bool			doSetByName, doSetByNumber, foundFont;
    PsychWindowRecordType	*windowRecord;
#if PSYCH_SYSTEM == PSYCH_OSX
    PsychFontStructType		*fontRecord;
#endif
    int				oldTextFontNumber, inputTextFontNumber;
    char			*oldTextFontName, *inputTextFontName;
    
    //all subfunctions should have these two lines.  
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    //check for valid number of arguments
    PsychErrorExit(PsychRequireNumInputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));   	
    PsychErrorExit(PsychCapNumOutputArgs(2)); 
    
    //Get the window record
    PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
    
    //Save the old text font value and return it.
    oldTextFontNumber=windowRecord->textAttributes.textFontNumber;
    PsychCopyOutDoubleArg(2, FALSE, (double)oldTextFontNumber);
    oldTextFontName=(char*) windowRecord->textAttributes.textFontName;
    PsychCopyOutCharArg(1, FALSE, oldTextFontName); 
    
	
    //Fetch and set the new font if specified by name or number
    PsychCheckInputArgType(2, kPsychArgOptional, PsychArgType_double | PsychArgType_char);  //if the argument is there check that it is the right type.
    doSetByNumber= PsychCopyInIntegerArg(2, kPsychArgAnything, &inputTextFontNumber);
    doSetByName= PsychAllocInCharArg(2, kPsychArgAnything, &inputTextFontName);
    foundFont=0;
#if PSYCH_SYSTEM == PSYCH_OSX
    if(doSetByNumber)
        foundFont=PsychGetFontRecordFromFontNumber(inputTextFontNumber, &fontRecord);
    if(doSetByName)
        foundFont=PsychGetFontRecordFromFontFamilyNameAndFontStyle(inputTextFontName, windowRecord->textAttributes.textStyle, &fontRecord);
    if(foundFont) {
        strncpy((char*) windowRecord->textAttributes.textFontName, (char*) fontRecord->fontFMFamilyName, 255);
        windowRecord->textAttributes.textFontNumber= fontRecord->fontNumber;
    }
	else {
		// Font not found. Is this textrenderer 2 with a font given by name?
		if (doSetByName && (PsychPrefStateGet_TextRenderer() > 1)) {
			// Yes: Must be a special font specifier string for the renderer plugin. Just assign it directly:
			strncpy((char*) windowRecord->textAttributes.textFontName, inputTextFontName, 255);

			// Don't have a valid fontNumber: Just assign a zero...
			windowRecord->textAttributes.textFontNumber = 0;
		}
	}
    
    return(PsychError_none);
#else
    // Special case for MS-Windows and Linux:
    if(doSetByNumber) printf("PTB-WARNING: Sorry, selecting font by number in Screen('TextFont') is not yet supported on Windows or Linux. Command ignored.\n");
    if(doSetByName && (strncmp(windowRecord->textAttributes.textFontName, inputTextFontName, 255 )!=0)) {
      strncpy(windowRecord->textAttributes.textFontName, inputTextFontName, 255);
      windowRecord->textAttributes.textFontNumber= 0;
      // Set the rebuild flag:
      windowRecord->textAttributes.needsRebuild=TRUE;
    }

    return(PsychError_none);
#endif
}
