/*
	PsychToolbox2/Source/OSX/FontInfo/PsychFontGlue.h
	
	PLATFORMS:	
	
		OS X  
				
	AUTHORS:
	
         Allen Ingling		awi		Allen.Ingling@nyu.edu
         Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
		11/24/03	awi		Wrote it.
		3/7/06		awi		Changed references from "Font" to "FontInfo".  The function had been previously renamed, but not all references updated. 
							
	DESCRIPTION:
        

*/

//include once
#ifndef PSYCH_IS_INCLUDED_PsychFontGlue
#define PSYCH_IS_INCLUDED_PsychFontGlue

#include "Psych.h"


typedef struct _PsychFontLocaleStructType_{
        Str255			language;
        Str255			languageVariant;
        Str255			region;
        Str255			regionVariant;
        Str255			fullName;
} PsychFontLocaleStructType;


typedef struct _PsychFontMetricsStructType_{
        double			ascent;
        double			descent;
        double		        leading;
        double			avgAdvanceWidth;
        double			minLeftSideBearing;
        double			minRightSideBearing;
        double      		stemWidth;
        double     		stemHeight;
        double      		capHeight;
        double        		xHeight;
        double    		italicAngle;
        double			underlinePosition;
        double			underlineThickness;
} PsychFontMetricsStructType;    

#define 	kPsychMaxFontFileNameChars		1024
typedef struct _PsychFontStructType_	*PsychFontStructPtrType;
typedef struct _PsychFontStructType_{

        //list management
        int								fontNumber;
        PsychFontStructPtrType			next;
 
        //Names associated with the font.  We could also add the ATS
        Str255					fontFMName;
        Str255					fontFMFamilyName;
        Str255					fontFamilyQuickDrawName;
        Str255					fontPostScriptName;
        char					fontFile[kPsychMaxFontFileNameChars];  	
			
        //depricated font style stuff for compatability with OS 9 scripts
        FMFontStyle				fontFMStyle;
        int						fontFMNumStyles;
        
        //Retain the the Font Manager (FM) and Apple Type Services (ATS) references to the font.
        //We could use these to tie font families to font names, for example 
        // FontInfo('GetFontsFromFontFamilyName') or FontInfo('GetFontFamilyFromFont');
        ATSFontRef					fontATSRef;
        FMFont						fontFMRef;
        ATSFontFamilyRef			fontFamilyATSRef;
        FMFontFamily				fontFamilyFMRef;
        
        //font's language and country
        PsychFontLocaleStructType		locale;			
        
        //font metrics
        PsychFontMetricsStructType		horizontalMetrics;
        PsychFontMetricsStructType		verticalMetrics;
} PsychFontStructType;

// function prototypes//functions for ATSU

//functions for handling psych font lists
int 			PsychFreeFontList(void);
PsychFontStructPtrType	PsychGetFontListHead(void);
int 			PsychGetFontListLength(void);
psych_bool			PsychGetFontRecordFromFontNumber(int fontIndex, PsychFontStructType **fontStruct);
psych_bool			PsychGetFontRecordFromFontFamilyNameAndFontStyle(char *fontName, FMFontStyle fontStyle, PsychFontStructType **fontStruct);
void 			PsychCopyFontRecordsToNativeStructArray(int numFonts, PsychFontStructType **fontStructs, PsychGenericScriptType **nativeStructArray);  

//functions for dealing with Font Manager styles
int				PsychFindNumFMFontStylesFromStyle(FMFontStyle fmStyleFlag);
void 			PsychGetFMFontStyleNameFromIndex(int styleIndex, FMFontStyle fontStyle, char *styleName, int styleNameLength);

//end include once
#endif
