/*
	PsychToolbox2/Source/OSX/FontInfo/PsychFontGlue.c
	
	PLATFORMS:	
	
		OS X
				
	AUTHORS:
	
		Allen Ingling		awi		Allen.Ingling@nyu.edu
        Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
		11/24/03		awi		Wrote it 
		3/7/06			awi		Changed references from "Font" to "FontInfo".  The function had been previously renamed, but not all references updated.		
							
	DESCRIPTION:
        

*/
#include "PsychFontGlue.h"

#define 	PsychFMIterationScope			kFMLocalIterationScope
#define 	kPsychMaxFontFileNameChars		1024

//fixes discrepancy between Apple documentation & Apple ATS.h header file.
#define kATSOptionFlagsUnrestrictedScope		kATSOptionFlagsUnRestrictedScope	

/*
Note: When you iterate using a local context and an unrestricted scope, you
enumerate the default font families. This includes all globally activated font
families and those activated locally to your application.  In other words, its
what's available to your app.  
*/
#define PSYCH_ATS_ITERATOR_CONTEXT			kATSFontContextLocal
#define PSYCH_ATS_ITERATOR_SCOPE			kATSOptionFlagsUnrestrictedScope				

//declare file local functions
static void						PsychInitFontList(void);
static PsychFontStructPtrType	PsychFontListHeadKeeper(psych_bool set, PsychFontStructPtrType value);

#ifndef PTBMODULE_Screen
// Provide surrogate for FontInfo() build:
int PsychPrefStateGet_Verbosity(void)
{
    // Shut up all time:
    return(0);
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions for maintaining and accessing psych font lists.

/*
    PsychFontListHeadKeeper()
    
    A file-local function which enforces limited access to the font list head by encapsulating it within a function.
    
    Only PsychGetFontListHead(), PsychInitFontList() and PsychFreeFontList should call PsychFontListHeadKeeper().
    
    Other functions needing access to the font list should get the list head only through 
    PsychGetFontListHead().
    
*/
static PsychFontStructPtrType	PsychFontListHeadKeeper(psych_bool set, PsychFontStructPtrType value)
{
    static PsychFontStructPtrType		fontListHead=NULL;

    if(set){
        if(value && fontListHead)
            PsychErrorExitMsg(PsychError_internal, "Attempt to set new font list head when one is already set."); 
        fontListHead=value;
        return(fontListHead);
    }else
        return(fontListHead);
}



/*
    GetFontListHead(void)
    
    Get the list head of the font list which we built.  To insure that the font list has been initialized before
    it is accessed  We go through GetFontListHead() instead of allowing direct access to the font list variable.   
*/
PsychFontStructPtrType	PsychGetFontListHead(void)
{
    PsychFontStructPtrType	fontListHead;
    
    fontListHead=PsychFontListHeadKeeper(FALSE, NULL);
    if(!fontListHead){
        PsychInitFontList();
        fontListHead=PsychFontListHeadKeeper(FALSE, NULL);
    }
    return(fontListHead);
}





/*
    PsychFreeFontList()
    
    Free memory used to hold the Psychtoolbox's font list.  Set the list head to NULL when we
    are done to flag an empty list. 
*/
int PsychFreeFontList(void)
{
    PsychFontStructPtrType	current, previous;
    int				numFreed;
    
    numFreed=0;
    previous=NULL;
    current=PsychFontListHeadKeeper(FALSE, NULL); //get the font list head.
    while(current){
        if(previous){		
            free(previous);
            ++numFreed;
        }
        previous=current;
        current=current->next;
    }
    if(previous){		
        free(previous);			
        ++numFreed;
    }
    PsychFontListHeadKeeper(TRUE, NULL); //set the font list head to null to mark an empty list.
        
    return(numFreed);
}


/*
    PsychInitFontList()
    
    Build a list of records describing installed fonts.
    
    
*/
void PsychInitFontList(void)
{
    ATSFontRef			tempATSFontRef;
#ifdef __LP64__
    CTFontRef           tempCTFontRef;
#endif
    //for font structures
    PsychFontStructPtrType	fontListHead, fontRecord, previousFontRecord;
    //for ATI font iteration
    ATSFontIterator		fontIterator;
    OSStatus			halt;
    //for font field names
    CFStringRef 		cfFontName;
    int                 i, j;
    psych_bool			resultOK;
    //for font file
    OSStatus			osStatus;			
    FSSpec              fontFileSpec;
    FSRef               fontFileRef;
    //for the font metrics
    ATSFontMetrics		horizontalMetrics;
    ATSFontMetrics		verticalMetrics;
    //for info from Font Manager 
    FMFontStyle			fmStyle;
    OSStatus			fmStatus;
    Str255				fmFontFamilyNamePString;
    //whatever
    Str255              fontFamilyQuickDrawNamePString;
    TextEncoding		textEncoding;
    OSStatus			scriptInfoOK;
    ScriptCode			scriptCode;
    LangCode			languageCode;
    OSStatus			localOK;
    LocaleRef			locale;
	psych_bool			trouble = FALSE;
	psych_bool			reportTrouble = TRUE;

    fontListHead=PsychFontListHeadKeeper(FALSE, NULL); //get the font list head.
    if(fontListHead) PsychErrorExitMsg(PsychError_internal, "Attempt to set new font list head when one is already set.");
        
    fontRecord = NULL;
    fontIterator = NULL;
    halt = ATSFontIteratorCreate(PSYCH_ATS_ITERATOR_CONTEXT, NULL, NULL, PSYCH_ATS_ITERATOR_SCOPE, &fontIterator);
    i = 0;
	
    while (halt==noErr) {
        // Give repair hints early. Experience shows we might crash during enumeration of a
        // corrupt OSX font database, so make sure we get out the helpful message as early as possible. Doing
        // this (just) at the end of enumeration might be too late - we might never get there...
        // However, allow user to suppress the hint, as this happens quite regularly on 32-Bit OSX:
        if (reportTrouble && trouble && PsychPrefStateGet_Verbosity() > 2) {
            reportTrouble = FALSE; // Only show this hint once, not for every invalid font.
            printf("\nPTB-HINT: =============================================================================================================================\n");
            printf("PTB-HINT: At least one font on this system has issues and can not be accessed by Psychtoolbox. If you want to know which font(s) make\n");
            printf("PTB-HINT: trouble, do a 'clear all' and rerun your script with Screen()'s verbosity level set to at least 4 for more diagnostic output.\n");
#ifndef __LP64__
            printf("PTB-HINT: On 32-Bit Matlab under OSX it is quite normal that a few fonts can't be enumerated and used properly and there is no known\n");
            printf("PTB-HINT: way to fix this, so the following tips will not help you to resolve this issues, but upgrading to 64-Bit Psychtoolbox helps.\n");
#endif
            printf("PTB-HINT: The following tips may help you to resolve font issues:\n");
            printf("PTB-HINT: Go to the Application folder and open the 'Font Book' application. It allows you to check and repair your font database.\n");
            printf("PTB-HINT: Run its 'Validate' function on all installed fonts. Another thing you could try is downloading and running the free\n");
            printf("PTB-HINT: FontNuke application (Google will find it for you) to regenerate corrupt OSX font caches. Good luck!\n");
            printf("PTB-HINT: You can suppress this hint by choosing a verbosity level for Screen() of 2 or lower.\n");
            printf("PTB-HINT: =============================================================================================================================\n\n");
        }

        halt=ATSFontIteratorNext(fontIterator, &tempATSFontRef);
        if(halt==noErr){
            //create a new  font  font structure.  Set the next field  to NULL as  soon as we allocate the font so that if 
            //we break with an error then we can find the end when we  walk down the linked list. 
            fontRecord=(PsychFontStructPtrType) calloc(1, sizeof(PsychFontStructType));
            fontRecord->next=NULL;

            //Get  FM and ATS font and font family references from the ATS font reference, which we get from iteration.
            fontRecord->fontATSRef=tempATSFontRef; 
            fontRecord->fontFMRef=FMGetFontFromATSFontRef(fontRecord->fontATSRef);

#ifndef __LP64__
            fmStatus=FMGetFontFamilyInstanceFromFont(fontRecord->fontFMRef, &(fontRecord->fontFamilyFMRef), &fmStyle);
            fontRecord->fontFamilyATSRef = FMGetATSFontFamilyRefFromFontFamily(fontRecord->fontFamilyFMRef);
#else            
            // Create CTFont from given ATSFontRef. Available since OSX 10.5
            tempCTFontRef = CTFontCreateWithPlatformFont(fontRecord->fontATSRef, 0.0, NULL, NULL);
            if (tempCTFontRef) {
                // Get font family name from CTFont:
                CFStringRef cfFamilyName = CTFontCopyFamilyName(tempCTFontRef);
                
                // Retrieve symbolic traits of font -- the closest equivalent of the fmStyle from the
                // good'ol fontManager:
                CTFontSymbolicTraits ctTraits = CTFontGetSymbolicTraits(tempCTFontRef);
                
                // Remap new trait constants to old constants for later Screen('TextStyle') matching.
                fmStyle = 0;
                if (ctTraits & kCTFontBoldTrait) fmStyle |= 1;
                if (ctTraits & kCTFontItalicTrait) fmStyle |= 2;
                if (ctTraits & kCTFontCondensedTrait) fmStyle |= 32;
                if (ctTraits & kCTFontExpandedTrait) fmStyle |= 64;
                
                // CTFont no longer needed:
                CFRelease(tempCTFontRef);
                
                // Convert to C-String and assign:
                resultOK = cfFamilyName && CFStringGetCString(cfFamilyName, (char*) fontRecord->fontFMFamilyName, 255, kCFStringEncodingASCII);
                if(!resultOK){
                    if (cfFamilyName) CFRelease(cfFamilyName);
                    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to retrieve font family name for font... Defective font?!? Skipped this entry...\n");
                    trouble = TRUE;
                    continue;
                }
                
                // Get ATSRef for font family:
                fontRecord->fontFamilyATSRef = ATSFontFamilyFindFromName(cfFamilyName, kATSOptionFlagsDefault);
                CFRelease(cfFamilyName);
            }
            else {
                if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to retrieve CTFontRef for font... Defective font?!? Skipped this entry...\n");
                trouble = TRUE;
                continue;
            }
#endif

            //get the font name and set the the corresponding field of the struct
            if (ATSFontGetName(fontRecord->fontATSRef, kATSOptionFlagsDefault, &cfFontName)!=noErr) {
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to query font name in ATSFontGetName! OS-X font handling screwed up?!? Skipped this entry...\n");
				trouble = TRUE;
				continue;
            }
            
            resultOK = cfFontName && CFStringGetCString(cfFontName, (char*) fontRecord->fontFMName, 255, kCFStringEncodingASCII);
            if(!resultOK){
                if (cfFontName) CFRelease(cfFontName);
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to convert fontFMName CF string to char string. Defective font?!? Skipped this entry...\n");
				trouble = TRUE;
				continue;
            }
            CFRelease(cfFontName);

            //get the font postscript name and set the corresponding field of the struct
            if (ATSFontGetPostScriptName(fontRecord->fontATSRef, kATSOptionFlagsDefault, &cfFontName)!=noErr) {
                if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: The following font makes trouble: %s. Please REMOVE the offending font file from your font folders and restart Matlab. Skipped entry for now...\n", fontRecord->fontFMName);
				trouble = TRUE;
				continue;
            }

            resultOK = cfFontName && CFStringGetCString(cfFontName, (char*) fontRecord->fontPostScriptName, 255, kCFStringEncodingASCII); //kCFStringEncodingASCII matches MATLAB for 0-127
            if(!resultOK){
                if (cfFontName) CFRelease(cfFontName);
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to convert fontPostScriptName CF string to char string for font %s. Defective font?!? Skipped this entry...\n", fontRecord->fontFMName);
				trouble = TRUE;
				continue;
            }
            CFRelease(cfFontName);

            // Get the QuickDraw name of the font:
            fontRecord->fontFamilyQuickDrawName[0] = 0;
            fontFamilyQuickDrawNamePString[0] = 0;
            ATSFontFamilyGetQuickDrawName(fontRecord->fontFamilyATSRef, fontFamilyQuickDrawNamePString);
            #ifndef __LP64__
            CopyPascalStringToC(fontFamilyQuickDrawNamePString, (char*) fontRecord->fontFamilyQuickDrawName);
            #else
            for (j = 0; j < fontFamilyQuickDrawNamePString[0]; j++) fontRecord->fontFamilyQuickDrawName[j] = fontFamilyQuickDrawNamePString[j+1];
            fontRecord->fontFamilyQuickDrawName[j] = 0;
            #endif

            #ifndef __LP64__
            // Get the font file used for this font
            osStatus = ATSFontGetFileSpecification(fontRecord->fontATSRef, &fontFileSpec);
            if(osStatus != noErr) {
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to get the font file specifier for font %s. Defective font?!? Skipped this entry...\n", fontRecord->fontFMName);
				trouble = TRUE;
				continue;
			}

            FSpMakeFSRef(&fontFileSpec, &fontFileRef);
            #else
            // 64-Bit version, available since OSX 10.5:
            osStatus = ATSFontGetFileReference(fontRecord->fontATSRef, &fontFileRef);
            if(osStatus != noErr) {
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to get the font file specifier for font %s. Defective font?!? Skipped this entry...\n", fontRecord->fontFMName);
				trouble = TRUE;
				continue;
			}
            #endif

            osStatus= FSRefMakePath(&fontFileRef, (UInt8*) fontRecord->fontFile, (UInt32)(kPsychMaxFontFileNameChars - 1));
            if(osStatus!=noErr){
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to get the font file path for font %s. Defective font?!? Skipped this entry...\n", fontRecord->fontFMName);
				trouble = TRUE;
				continue;
            }

            //get the font metrics of this font.
            //Explicit copy between fields to make it clear what is going one, will likely have to mix & match between native and Psych structures for different
            //platforms.
            ATSFontGetHorizontalMetrics(fontRecord->fontATSRef, kATSOptionFlagsDefault, &horizontalMetrics);
            ATSFontGetVerticalMetrics(fontRecord->fontATSRef, kATSOptionFlagsDefault, &verticalMetrics);
            //horizontal
            fontRecord->horizontalMetrics.ascent=			horizontalMetrics.ascent;
            fontRecord->horizontalMetrics.descent=			horizontalMetrics.descent;
            fontRecord->horizontalMetrics.leading=			horizontalMetrics.leading;
            fontRecord->horizontalMetrics.avgAdvanceWidth=		horizontalMetrics.avgAdvanceWidth;
            fontRecord->horizontalMetrics.minLeftSideBearing=		horizontalMetrics.minLeftSideBearing;
            fontRecord->horizontalMetrics.minRightSideBearing=		horizontalMetrics.minRightSideBearing;
            fontRecord->horizontalMetrics.stemWidth=			horizontalMetrics.stemWidth;
            fontRecord->horizontalMetrics.stemHeight=			horizontalMetrics.stemHeight;
            fontRecord->horizontalMetrics.capHeight=			horizontalMetrics.capHeight;
            fontRecord->horizontalMetrics.xHeight=			horizontalMetrics.xHeight;
            fontRecord->horizontalMetrics.italicAngle=			horizontalMetrics.italicAngle;
            fontRecord->horizontalMetrics.underlinePosition=		horizontalMetrics.underlinePosition;
            fontRecord->horizontalMetrics.underlineThickness=		horizontalMetrics.underlineThickness;
            fontRecord->horizontalMetrics.underlineThickness=		horizontalMetrics.underlineThickness;
            //vertical
            fontRecord->verticalMetrics.ascent=				verticalMetrics.ascent;
            fontRecord->verticalMetrics.descent=			verticalMetrics.descent;
            fontRecord->verticalMetrics.leading=			verticalMetrics.leading;
            fontRecord->verticalMetrics.avgAdvanceWidth=		verticalMetrics.avgAdvanceWidth;
            fontRecord->verticalMetrics.minLeftSideBearing=		verticalMetrics.minLeftSideBearing;
            fontRecord->verticalMetrics.minRightSideBearing=		verticalMetrics.minRightSideBearing;
            fontRecord->verticalMetrics.stemWidth=			verticalMetrics.stemWidth;
            fontRecord->verticalMetrics.stemHeight=			verticalMetrics.stemHeight;
            fontRecord->verticalMetrics.capHeight=			verticalMetrics.capHeight;
            fontRecord->verticalMetrics.xHeight=			verticalMetrics.xHeight;
            fontRecord->verticalMetrics.italicAngle=			verticalMetrics.italicAngle;
            fontRecord->verticalMetrics.underlinePosition=		verticalMetrics.underlinePosition;
            fontRecord->verticalMetrics.underlineThickness=		verticalMetrics.underlineThickness;
            fontRecord->verticalMetrics.underlineThickness=		verticalMetrics.underlineThickness;

            // On 32-Bit use Font Manager to get the FM  font family name font style.
            // This has been already done above for 64-Bit:
            #ifndef __LP64__
            fmStatus=FMGetFontFamilyName(fontRecord->fontFamilyFMRef, fmFontFamilyNamePString);
            if(fmStatus!=noErr){
                // Sadly, this failure is now expected behaviour for a few fonts under 32-Bit runtime when used on modern OSX versions, e.g., 10.5 and later,
                // see forum message #16109 for the symptoms:
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: In font initialization: Failed to get the fontFMFamilyName for font %s. Defective font?!? Skipped this entry...\n", fontRecord->fontFMName);
				trouble = TRUE;
				continue;
            }
            CopyPascalStringToC(fmFontFamilyNamePString, (char*) fontRecord->fontFMFamilyName);
            #endif
            
            fontRecord->fontFMStyle=fmStyle;
            fontRecord->fontFMNumStyles=PsychFindNumFMFontStylesFromStyle(fmStyle);
            fontRecord->fontFMNumStyles= fontRecord->fontFMNumStyles ? fontRecord->fontFMNumStyles : 1; //because the name is "normal" even if there are no styles.  
            // Get the locale info which is a property of the font family:
            // No error checking is done here, because many (most?) fonts miss the information,
            // so we would error-out all the time and this is non-critical for us:
            textEncoding=ATSFontFamilyGetEncoding(fontRecord->fontFamilyATSRef);
            scriptInfoOK=RevertTextEncodingToScriptInfo(textEncoding, &scriptCode, &languageCode, NULL);
            localOK=LocaleRefFromLangOrRegionCode(languageCode, kTextRegionDontCare, &locale); 
            localOK |= LocaleRefGetPartString(locale, kLocaleLanguageMask, 255, (char*) fontRecord->locale.language);			fontRecord->locale.language[255]='\0'; 
            localOK |= LocaleRefGetPartString(locale, kLocaleLanguageVariantMask, 255, (char*) fontRecord->locale.languageVariant);	fontRecord->locale.languageVariant[255]='\0';
            localOK |= LocaleRefGetPartString(locale, kLocaleRegionMask, 255, (char*) fontRecord->locale.region);			fontRecord->locale.region[255]='\0';
            localOK |= LocaleRefGetPartString(locale, kLocaleRegionVariantMask, 255, (char*) fontRecord->locale.regionVariant);	fontRecord->locale.regionVariant[255]='\0';
            localOK |= LocaleRefGetPartString(locale, kLocaleAllPartsMask, 255, (char*) fontRecord->locale.fullName);		fontRecord->locale.fullName[255]='\0';
            
			// Init for fontRecord (nearly) finished.
			
			// Set this fontRecord as head of font-list, or enqueue it in existing list:
            if(i==0) {
                PsychFontListHeadKeeper(TRUE, fontRecord); 
            }
			else {
                previousFontRecord->next=fontRecord;
            }
			
			// Set the font number field of the struct
            fontRecord->fontNumber=i+1;

            // Increment the font index and update the next font pointer
            ++i;
            previousFontRecord=fontRecord;
        }else if(halt == kATSIterationScopeModified){
            //exit because the font database changed during this loop.
            PsychFreeFontList();
            if (fontIterator) ATSFontIteratorRelease(&fontIterator);
            PsychErrorExitMsg(PsychError_system, "The system font database was modified during font list setup. Please 'clear all' and restart your script.");
        }
		// Next parse iteration in system font database...
    }
	
    if (fontIterator) ATSFontIteratorRelease(&fontIterator);
    
    if (halt != kATSIterationCompleted) {
        PsychFreeFontList();
        trouble = TRUE;
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Font iteration during enumeration terminated prematurely. OS-X Font database corrupted?!?");
    }
	
	// Did we get a hand on at least one font?
	if (i==0) {
        PsychFreeFontList();
        trouble = TRUE;
		if (PsychPrefStateGet_Verbosity() > 0) {
			printf("PTB-ERROR: In font initialization: Could not even retrieve one valid font from the system! The OS-X font database must be corrupt.\n");
			printf("PTB-ERROR: Will try to continue but will likely abort if your code tries to call any of the font handling or text drawing functions.\n");
		}
	}
	
    // If there was some trouble and it wasn't reported yet, then report it now at the end if either
    // verbosity level is at least 3 (our default level), or if the error was especially serious and
    // lead to a premature abortion of font iteration or complete failure to find even a single valid
    // font. In the latter cases, we must report the trouble, regardless of verbosity level. In the
    // former case, probably only a few fonts had trouble, so we allow the user to suppress such messages
    // by lowering the verbosity to warning level or lower:
	if (reportTrouble && trouble && ((PsychPrefStateGet_Verbosity() > 2) || (halt != kATSIterationCompleted) || (i == 0))) {
		printf("PTB-HINT: =============================================================================================================================\n");
        printf("PTB-HINT: At least one font on this system has issues and can not be accessed by Psychtoolbox. If you want to know which font(s) make\n");
        printf("PTB-HINT: trouble, do a 'clear all' and rerun your script with Screen()'s verbosity level set to at least 4 for more diagnostic output.\n");
        printf("PTB-HINT: The following tips may help you to resolve font issues:\n");
		printf("PTB-HINT: Go to the Application folder and open the 'Font Book' application. It allows you to check and repair your font database.\n");
        printf("PTB-HINT: Run its 'Validate' function on all installed fonts. Another thing you could try is downloading and running the free\n");
        printf("PTB-HINT: FontNuke application (Google will find it for you) to regenerate corrupt OSX font caches. Good luck!\n");
        printf("PTB-HINT: In case of non-fatal errors, you can suppress this hint by choosing a verbosity level for Screen() of 2 or lower.\n");
		printf("PTB-HINT: =============================================================================================================================\n");
	} 

    // Font database ready for use.
	return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions which use font lists.


/*
    PsychGetFontListLength()
    
    Walk down the list and count the number of elements in the list
*/
int PsychGetFontListLength(void)
{
    PsychFontStructPtrType		current;
    int 				i;
    
    i=0;
    for(current=PsychGetFontListHead();  current; current=current->next)
        ++i;
    return(i);
}




/*
    PsychGetFontRecordFromFontNumber()
    
    Accept a font number.  Lookup the specified font record in the list of font records.  If the font is found in the 
    list then set the fontStruct pointer to point to the font record.  Otherwise set fontStruct to NULL and return FALSE.   
        
    Font numbers are something the Psychtoolbox assigns when building its font list.  They are not genine properties of the 
    font.
*/
psych_bool	PsychGetFontRecordFromFontNumber(int fontIndex, PsychFontStructType **fontStruct)
{
    PsychFontStructPtrType	current;
    psych_bool			found;

    found=0;
    for(current=PsychGetFontListHead();  current; current=current->next){
        if(current->fontNumber==fontIndex){
            found=1;
            *fontStruct=current;
            break;
        }
    }
    if(!found)
        *fontStruct=NULL;

    return(found);
}


/*
    PsychGetFontRecordFromFontNameAndFontStyle()
    
    Accept a font family name and a font style.  Lookup the specified font in the list of font records.  If the font is found in the 
    list then set the fontStruct pointer to point to the font record.  Otherwise set fontStruct to NULL and return FALSE.    
*/
psych_bool	PsychGetFontRecordFromFontFamilyNameAndFontStyle(char *fontFamilyName, FMFontStyle fontStyle, PsychFontStructType **fontStruct)
{
    PsychFontStructPtrType	current;
    psych_bool			nameMatch, styleMatch, found;

    found=0;
    for(current=PsychGetFontListHead();  current; current=current->next){
        nameMatch= !(strcmp((char*) current->fontFMFamilyName, fontFamilyName));
        styleMatch= current->fontFMStyle==fontStyle;
        found= nameMatch && styleMatch;
        if(found){
            *fontStruct=current;
            break;
        }
    }
    if(!found)
        *fontStruct=NULL;
    return(found);
}


/*
    PsychGetNumFontRecordsFromFontFamilyName()
    
    Accept a font family name string and a pointer to strorage.  Set the contents of the storage to hold 
    an array of pointers to font records for the fonts which are members of the  specified font family.
    
    Return the number of elements.
    
    Passing NULL in for pFontStructArray will cause PsychGetFontRecordsFromFontFamilyName() to return the tally of 
    the family's member fonts without storing the pointers. To determine the size of the pFontStructArray array 
    storage, first call PsychGetFontRecordsFromFontFamilyName, passing NULL in for pFontStructArray.
    
*/
int	PsychMemberFontsFromFontFamilyName(char *fontFamilyName, PsychFontStructPtrType *pFontStructArray)
{
    int				i;
    psych_bool			nameMatch;
    PsychFontStructPtrType	current;
    
    i=0;
    for(current=PsychGetFontListHead();  current; current=current->next){
        nameMatch= !(strcmp((char*) current->fontFMFamilyName, fontFamilyName));
        if(nameMatch){
            if(pFontStructArray)
                pFontStructArray[i]=current;
            i++;
        }
    }
    return(i);
}


/*
   PsychCopyFontRecordToNativeStruct()
   
   Accept an array of pointers to psych font structures and the lenght of the array.  
   Allocate and fill a native font struct array for the outside scripting environment   
*/    
void PsychCopyFontRecordsToNativeStructArray(int numFonts, PsychFontStructType **fontStructs, PsychGenericScriptType **fontStructArray)  
{
    int				i, arrayIndex;
    PsychFontStructPtrType	fontElement;
    PsychGenericScriptType	*styleNameList;
    Str255			styleName;
    //for the outer font struct
    int				numFontFieldNames=10;
    const char 			*fontFieldNames[]={"number", "name", "QuickDrawName", "PostScriptName", "familyName", "styleCode", 
                                                            "styleNames", "file",  "metrics", "locale" };
    //for the metrics meta struct
    int				numMetaMetricsFieldNames=2;
    const char 			*metaMetricsFieldNames[]={"verticalMetrics", "horizontalMetrics"};
    PsychGenericScriptType	*metaMetrics;
    //for the metrics structs
    int				numMetricsFieldNames=13;
    const char 			*metricsFieldNames[]={"ascent", "descent", "leading", "avgAdvanceWidth", "minLeftSideBearing", "minRightSideBearing",
                                                    "stemWidth","stemHeight","capHeight","xHeight","italicAngle","underlinePosition","underlineThickness" };
    PsychGenericScriptType	*horizontalMetrics, *verticalMetrics;
    //for locale
    int				numLocaleFieldNames=5;
    const char 			*localeFieldNames[]={"language", "languageVariant", "region", "regionVariant", "fullName"};
    PsychGenericScriptType	*locale;

    
    //convert the linked list of font structures to a struct array for the scripting environment.
    PsychAllocOutStructArray(kPsychNoArgReturn, FALSE, numFonts, numFontFieldNames, fontFieldNames, fontStructArray);
    for(arrayIndex=0; arrayIndex<numFonts; ++arrayIndex){
        fontElement=fontStructs[arrayIndex];
        PsychSetStructArrayDoubleElement("number", arrayIndex, fontElement->fontNumber, *fontStructArray);
        PsychSetStructArrayStringElement("name", arrayIndex, (char*) fontElement->fontFMName, *fontStructArray);
        PsychSetStructArrayStringElement("QuickDrawName", arrayIndex, (char*) fontElement->fontFamilyQuickDrawName, *fontStructArray);
        PsychSetStructArrayStringElement("PostScriptName", arrayIndex, (char*) fontElement->fontPostScriptName, *fontStructArray);
        PsychSetStructArrayStringElement("familyName", arrayIndex, (char*) fontElement->fontFMFamilyName, *fontStructArray);
        //style
        PsychSetStructArrayDoubleElement("styleCode", arrayIndex, (double)fontElement->fontFMStyle, *fontStructArray);
        PsychAllocOutCellVector(kPsychNoArgReturn, FALSE, fontElement->fontFMNumStyles,  &styleNameList);
        for(i=0;i<fontElement->fontFMNumStyles;i++){
            PsychGetFMFontStyleNameFromIndex(i, fontElement->fontFMStyle, (char*) styleName, 255); 
            PsychSetCellVectorStringElement(i, (const char*) styleName, styleNameList);
        }
        PsychSetStructArrayNativeElement("styleNames", arrayIndex, styleNameList, *fontStructArray);
        PsychSetStructArrayStringElement("file", arrayIndex, fontElement->fontFile, *fontStructArray);
        //metrics
        PsychAllocOutStructArray(kPsychNoArgReturn, FALSE, 1, numMetaMetricsFieldNames, metaMetricsFieldNames,  &metaMetrics);
        PsychAllocOutStructArray(kPsychNoArgReturn, FALSE, 1, numMetricsFieldNames, metricsFieldNames, &horizontalMetrics);
        PsychAllocOutStructArray(kPsychNoArgReturn, FALSE, 1, numMetricsFieldNames, metricsFieldNames,  &verticalMetrics);
            //horizontal
            PsychSetStructArrayDoubleElement("ascent", 0, fontElement->horizontalMetrics.ascent, horizontalMetrics);
            PsychSetStructArrayDoubleElement("descent", 0, fontElement->horizontalMetrics.descent, horizontalMetrics);
            PsychSetStructArrayDoubleElement("leading", 0, fontElement->horizontalMetrics.leading, horizontalMetrics);
            PsychSetStructArrayDoubleElement("avgAdvanceWidth", 0, fontElement->horizontalMetrics.avgAdvanceWidth, horizontalMetrics);
            PsychSetStructArrayDoubleElement("minLeftSideBearing", 0, fontElement->horizontalMetrics.minLeftSideBearing, horizontalMetrics);
            PsychSetStructArrayDoubleElement("minRightSideBearing", 0, fontElement->horizontalMetrics.minRightSideBearing, horizontalMetrics);
            PsychSetStructArrayDoubleElement("stemWidth", 0, fontElement->horizontalMetrics.stemWidth, horizontalMetrics);
            PsychSetStructArrayDoubleElement("stemHeight", 0, fontElement->horizontalMetrics.stemHeight, horizontalMetrics);
            PsychSetStructArrayDoubleElement("capHeight", 0, fontElement->horizontalMetrics.capHeight, horizontalMetrics);
            PsychSetStructArrayDoubleElement("xHeight", 0, fontElement->horizontalMetrics.xHeight, horizontalMetrics);
            PsychSetStructArrayDoubleElement("italicAngle", 0, fontElement->horizontalMetrics.italicAngle, horizontalMetrics);
            PsychSetStructArrayDoubleElement("underlinePosition", 0, fontElement->horizontalMetrics.underlinePosition, horizontalMetrics);
            PsychSetStructArrayDoubleElement("underlineThickness", 0, fontElement->horizontalMetrics.underlineThickness, horizontalMetrics);
            //vertical
            PsychSetStructArrayDoubleElement("ascent", 0, fontElement->verticalMetrics.ascent, verticalMetrics);
            PsychSetStructArrayDoubleElement("descent", 0, fontElement->verticalMetrics.descent, verticalMetrics);
            PsychSetStructArrayDoubleElement("leading", 0, fontElement->verticalMetrics.leading, verticalMetrics);
            PsychSetStructArrayDoubleElement("avgAdvanceWidth", 0, fontElement->verticalMetrics.avgAdvanceWidth, verticalMetrics);
            PsychSetStructArrayDoubleElement("minLeftSideBearing", 0, fontElement->verticalMetrics.minLeftSideBearing, verticalMetrics);
            PsychSetStructArrayDoubleElement("minRightSideBearing", 0, fontElement->verticalMetrics.minRightSideBearing, verticalMetrics);
            PsychSetStructArrayDoubleElement("stemWidth", 0, fontElement->verticalMetrics.stemWidth, verticalMetrics);
            PsychSetStructArrayDoubleElement("stemHeight", 0, fontElement->verticalMetrics.stemHeight, verticalMetrics);
            PsychSetStructArrayDoubleElement("capHeight", 0, fontElement->verticalMetrics.capHeight, verticalMetrics);
            PsychSetStructArrayDoubleElement("xHeight", 0, fontElement->verticalMetrics.xHeight, verticalMetrics);
            PsychSetStructArrayDoubleElement("italicAngle", 0, fontElement->verticalMetrics.italicAngle, verticalMetrics);
            PsychSetStructArrayDoubleElement("underlinePosition", 0, fontElement->verticalMetrics.underlinePosition, verticalMetrics);
            PsychSetStructArrayDoubleElement("underlineThickness", 0, fontElement->verticalMetrics.underlineThickness, verticalMetrics);
        PsychSetStructArrayNativeElement("horizontalMetrics", 0, horizontalMetrics, metaMetrics);
        PsychSetStructArrayNativeElement("verticalMetrics", 0, verticalMetrics, metaMetrics);
        PsychSetStructArrayNativeElement("metrics", arrayIndex, metaMetrics, *fontStructArray);
        //locale
        PsychAllocOutStructArray(kPsychNoArgReturn, FALSE, 1, numLocaleFieldNames, localeFieldNames,  &locale);
            PsychSetStructArrayStringElement("language", 0, (char*) fontElement->locale.language, locale);
            PsychSetStructArrayStringElement("languageVariant", 0, (char*) fontElement->locale.languageVariant, locale);
            PsychSetStructArrayStringElement("region", 0, (char*) fontElement->locale.region, locale);
            PsychSetStructArrayStringElement("regionVariant", 0, (char*) fontElement->locale.regionVariant, locale);
            PsychSetStructArrayStringElement("fullName", 0, (char*) fontElement->locale.fullName, locale);
        PsychSetStructArrayNativeElement("locale", arrayIndex, locale, *fontStructArray);
    }
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//For dealing with Font Manager font styles.

#define		kPsychFMMaxStyles	9		//the max number of non-normal styles, the maximum index of FontStyleNameTable[].  
static char *FontStyleNameTable[] = {"normal",	"bold",	"italic", "underline", "outline", "unknownStyle3", "condense", "extend", "unknownStyle1", "unknownStyle2"};

/*  
    PsychFindNumFMFontStylesFromStyle()

    The minimum is 0 which means "normal" style, for which there is no flag.  (It is implied by the absense of other flags.)
    
*/
int	PsychFindNumFMFontStylesFromStyle(FMFontStyle fmStyleFlag)
{
    int			i, numStyles;
    FMFontStyle		testBit;
    
    numStyles=0;
    testBit=1;
    for(i=0;i<8;i++){
        testBit=1<<i;
        if(fmStyleFlag & testBit)
            ++numStyles;
    }
    return(numStyles);
}


 
                    
/*
    PsychGetFontStyleNameFromIndex()
    
    Get a string naming the nth font style which is set within an FMFontStyle value, where n is a number between 0  and the value
    returned by PsychFindNumFMFontStylesFromStyle() -1 when called on the FMFontStyle value.
    
    styleNameLength is the length of the string in characters, not the length in bytes.  Length in bytes = length in chars + 1.
    
*/
void PsychGetFMFontStyleNameFromIndex(int styleIndex, FMFontStyle fontStyle, char *styleName, int styleNameLength)
{
    int		i, totalStyles, styleCounter;
    
    
    //count up through the set bits until we reach the styleIndex then return the name of that bit position.
    //if we count to the end and find no bits then return "normal".  
    totalStyles=PsychFindNumFMFontStylesFromStyle(fontStyle);
    if(styleIndex > totalStyles)
        PsychErrorExitMsg(PsychError_internal, "Index value exceeds maximum");
    styleCounter=-1;
    for(i=0;i<=kPsychFMMaxStyles;i++){
        if( 1<<i & fontStyle)
            ++styleCounter;
        if(styleCounter == styleIndex){
            strncpy(styleName, FontStyleNameTable[i+1], styleNameLength);
            styleName[styleNameLength]='\0';
            return;
        }
    }
    if(styleIndex != 0)	
        //if the style is normal then there can be only one style and its index must be 0.    
        PsychErrorExitMsg(PsychError_internal, "indexed style not present in style value");
    strncpy(styleName, FontStyleNameTable[0], styleNameLength);
    styleName[styleNameLength]='\0';
    return;
} 
