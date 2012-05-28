/*
	PsychToolbox2/Source/OSX/FontInfo/FontInfo.h		

	PROJECTS: 
	
		FontInfo only.  

	AUTHORS:
	
        Allen Ingling		awi		Allen.Ingling@nyu.edu
        Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de

	PLATFORMS:
	
		OS X


	HISTORY:
	
		10/23/03  	awi		Created.
		3/7/06		awi		Changed references from "Font" to "FontInfo".  The function had been previously renamed, but not all references updated. 		 

	NOTES:

		What the OS 9 Psychtoolbox refers to as a "font"  Apple's Font Manager (FM) and Apple Type Services (ATS) APIs term a "font family".  
		In Font Mangager and ATS, a 'font' is a unique combination of a font family and a font style.

		Both ATS and Font Manager will enumerate both fonts and font families.  However, there appears to be no way in ATS to either 1) find the
		constituent fonts of a font family or 2) find the style of a font.  ATS alone is therefore useless for uncovring font information in 
		a form which matches the OS 9 Psychtoolbox functions.  

		Font Mangager on the other hand will both uncover the consituent fonts of a font family and reports the style of fonts.  To get the 
		font family and style of a particular font use FMGetFontFamilyInstanceFromFont():
			OSStatus	FMGetFontFamilyInstanceFromFont (FMFont iFont,
															FMFontFamily *oFontFamily,
															FMFontStyle *oStyle);
															
		The Psychtoolbox, given a font family and style, must uncover a font reference.  One way to  do that is to enumerate fonts until
		we find the font of the specified family and style.  A variation of this is to store values found by the font enumeration within
		a linked list of records.

		An alternative to enumerating fonts is to enumerate fonts families to find the family of the searched-for font and and then to enumerate that 
		family's constituent fonts to find the searched-for font with matching style.    

		Also, it does not matter whether we use ATS or FontManager to enumerate fonts becasue we can convert one type of reference into 
		the other.  In particular if we enumerate with ATS then we can obtain a FontManager referernce using FMGetATSFontRefFromFont();
			 ATSFontRef	FMGetATSFontRefFromFont(FMFont iFont);

		Note that ultimately we need to pass quartz an ATS font reference but we can get that by converting the Font Manager
		reference to an ATS reference using FMGetATSFontRefFromFont()

		If we enumerate ATS fonts then additionals information which we can uncover by converting them to Font Manager fonts is:
			-The associated font family.
			-The associcated font family info. 
			-The font stye.
				   

		About "activating" fonts.  This does not seem to be anything useful.  It seems that if we choose the correct context and scope then 
		we get the default font set and we don't need to worry about activating additional fonts.  This seems to be the domain of font utilites
		and does not concern us.  

		Apple says:
			You can control which fonts are available to your users by activating and
			deactivating fonts. Fonts are activated and deactivated in groups defined by their
			representation in the file system in the file formats supported by ATS for fonts.
		
	TO DO:

*/ 

//begin include once 
#ifndef PSYCH_IS_INCLUDED_FontInfo
#define PSYCH_IS_INCLUDED_FontInfo


//project includes 
#include "PsychFontGlue.h"

//prototypes for functions defined in this moduel

//Psych typical module support stuff
void InitializeSynopsis();
PsychError FontsExitFunction(void);
PsychError PsychDisplayATSSynopsis(void);
PsychError MODULEVersion(void); 

//module specific functions
PsychError FONTSNumFonts(void);
PsychError FONTSFonts(void);

//end include once
#endif
