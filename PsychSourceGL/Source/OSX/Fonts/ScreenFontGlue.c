/*
	PsychToolbox2/Source/OSX/FontInfo/ScreenFontGlue.c
	
	PLATFORMS:	This is the OS X  version.  
				
	AUTHORS:
	
		Allen Ingling		awi		Allen.Ingling@nyu.edu
        Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
		7/8/04		awi		Split it off from PsychFontGlue because PsychWindowRecordType in PsychSetATSUStyleAttributes
							was causing trouble for the FontInfo project
		3/7/06		awi		Changed references from "Font" to "FontInfo".  The function had been previously renamed, but not all references updated.
							
	DESCRIPTION:
        
*/

#include "ScreenFontGlue.h"

// Only include this function when building Screen(), but not FontInfo():
#ifdef PTBMODULE_Screen

// Disable warnings about deprecated API calls on OSX 10.7
// of which we are aware and that we can't remove as long as
// we need to stay compatible to 10.4 - 10.6
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

/*
	PsychSetATSUStyleAttributesFromPsychWindowRecord()
	
	Use the Font Manager attributes for the given window, which are stored in the PsychFontAtributeStruct within the window's record, to set the 
	attributes of the ATSU style.
	
	For now we only set the font ID, the color, and the size. Boldness, underline, condensed and extended we assume are brought along
	with the font ID , though perhaps not and we will have to set those also that.
	
	PsychSetATSUStyleAttributes() requires that the font number within the specified window record indicate a
	valid font.  
*/

void PsychSetATSUStyleAttributesFromPsychWindowRecord(ATSUStyle atsuStyle,  PsychWindowRecordType *winRec)
{
	PsychFontStructType		*psychFontRecord;
	int						psychColorSize;
	//for ATSU attributes
	ATSUFontID				atsuFontID;
	Fixed					atsuFontSize;
	ATSURGBAlphaColor		atsuFontColor;
	ATSStyleRenderingOptions atsuRenderOptions;
	GLdouble				colorVector[4];
	OSStatus				callError;
	ATSUAttributeTag		aaTags[] =  {kATSUFontTag, kATSUSizeTag, kATSURGBAlphaColorTag, kATSUStyleRenderingOptionsTag };
	ByteCount				aaSizes[] = {sizeof(ATSUFontID), sizeof(Fixed), sizeof(ATSURGBAlphaColor), sizeof(ATSStyleRenderingOptions) };
	ATSUAttributeValuePtr   aaValue[] = {&atsuFontID, &atsuFontSize, &atsuFontColor, &atsuRenderOptions};
	
	//set the font index
	PsychGetFontRecordFromFontNumber(winRec->textAttributes.textFontNumber, &psychFontRecord);
	if(psychFontRecord==NULL)
		PsychErrorExitMsg(PsychError_internal, "Failed to lookup the font from the font number");
	atsuFontID=psychFontRecord->fontFMRef;

	//set the font size
	atsuFontSize=Long2Fix((long)(winRec->textAttributes.textSize));

	//set the color
	PsychCoerceColorMode(&(winRec->textAttributes.textColor));
	PsychConvertColorToDoubleVector(&(winRec->textAttributes.textColor), winRec, colorVector);
	atsuFontColor.red=(float)colorVector[0];
	atsuFontColor.green=(float)colorVector[1];
	atsuFontColor.blue=(float)colorVector[2];
	atsuFontColor.alpha=(float)colorVector[3];
	
	// Set anti-aliasing mode (on/off): Default is to leave it up to the system to decide when
	// to anti-alias and when not. But this flag allows to always force anti-aliasin on or off:
	atsuRenderOptions = kATSStyleNoOptions;
	if(PsychPrefStateGet_TextAntiAliasing()==0) atsuRenderOptions = kATSStyleNoAntiAliasing; 
	if(PsychPrefStateGet_TextAntiAliasing()> 0) atsuRenderOptions = kATSStyleApplyAntiAliasing; 
	
	//assign attributes to the style object
	callError=ATSUSetAttributes(atsuStyle, 4, aaTags, aaSizes, aaValue);
}

#endif
