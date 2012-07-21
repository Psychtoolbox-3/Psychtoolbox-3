/*
	SCREENTextStyle.c	
  
	AUTHORS:

		Allen.Ingling@nyu.edu			awi
		mario.kleiner@tuebingen.mpg.de	mk
  
	PLATFORMS:	
	
		All.
    
	HISTORY:
	
		11/19/03	awi		Wrote it.
		10/12/04	awi		In useString: changed "SCREEN" to "Screen", and moved commas to inside [].

	DESCRIPTION:
  
		Sets the text style for the specified window record.
  
	NOTES:

*/

#include "Screen.h"

static char useString[] ="oldStyle=Screen('TextStyle', windowPtr [,style]);";
//                                         0           1           2
static char synopsisString[] = 
    "Get/set the font style for future text draws in this window. Useful values for "
    "style follow; they may be OR'd. See Inside Mac TextFace() for more. On M$-Windows "
    "and GNU/Linux, only a subset of these settings is honored - All settings are accepted, "
    "but some of them are silently ignored on Windows and Linux. On Windows the default "
	"text style is bold, because this increases readability on that system.\n"
    "0=normal,1=bold,2=italic,4=underline,8=outline,32=condense,64=extend.\n"
	"You can assign a default font style for new windows via a call to "
	"Screen('Preference', 'DefaultFontStyle'). The initial default font "
	"style is operating system dependent.\n"
	"Not all fonts support all style settings. Unsupported settings for the currently "
	"selected font will be silently ignored.";

static char seeAlsoString[] = "";

PsychError SCREENTextStyle(void) 
{
    psych_bool					doSetStyle;
    PsychWindowRecordType		*windowRecord;
    int							oldTextStyle, newTextStyle;

#if PSYCH_SYSTEM == PSYCH_OSX
    psych_bool                  foundFont;
    PsychFontStructType			*fontRecord;
#endif

    //all subfunctions should have these two lines.  
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    //check for valid number of arguments
    PsychErrorExit(PsychRequireNumInputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));   	
    PsychErrorExit(PsychCapNumOutputArgs(1)); 
    
    //Get the window record
    PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
    
    //Save the old text size value and return it.
    oldTextStyle=windowRecord->textAttributes.textStyle;
    PsychCopyOutDoubleArg(1, FALSE, (double)oldTextStyle);
    
    //Fetch and set the new size if it is specified. 
    doSetStyle= PsychCopyInIntegerArg(2, FALSE, &newTextStyle);
    if (doSetStyle) {
      windowRecord->textAttributes.needsRebuild|=(windowRecord->textAttributes.textStyle != newTextStyle) ? TRUE : FALSE;
      windowRecord->textAttributes.textStyle=newTextStyle;
	  
	  #if PSYCH_SYSTEM == PSYCH_OSX
	  // Need to update font name and number from changed style on OS/X:
	  foundFont = PsychGetFontRecordFromFontFamilyNameAndFontStyle((char*) windowRecord->textAttributes.textFontName, windowRecord->textAttributes.textStyle, &fontRecord);
	  if (foundFont) {
		strncpy((char*) windowRecord->textAttributes.textFontName, (const char*) fontRecord->fontFMFamilyName, 255);
		windowRecord->textAttributes.textFontNumber= fontRecord->fontNumber;
	  }
	  else {
		// Failed! Revert to old setting:
		windowRecord->textAttributes.textStyle = oldTextStyle;
	  }
	  
	  #endif
    }
	
    return(PsychError_none);
}
