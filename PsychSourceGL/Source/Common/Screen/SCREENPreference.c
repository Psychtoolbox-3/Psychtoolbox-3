/*
	SCREENPreference.c		
  
	AUTHORS:

		Allen.Ingling@nyu.edu		awi 
  
	PLATFORMS:
	
		Mac OS X MATLAB only.
    
	HISTORY:
  
		2/28/04		awi		Created.
		10/29/04	awi		Fixed two warnings, one where a const char should have been been char, the other the other way around.
		1/20/05		awi		Added a preference for debugging MakeTexture 
		5/30/05		mk		New preference setting screenSkipSyncTests -- allow skipping/shortening of some internal tests.
		5/30/05		mk		New preference setting screenVisualDebugLevel.
		3/07/05		awi		New preference SuppressAllWarnings.
                11/15/06        mk              New preference vbl & flip timestamping mode.
 
	DESCRIPTION:
  
		Read and set preferences.
  
		Note that we do don't do sanity checking on arguments provided in conjunction for unsuported preferences.  That might be conveivably be a mildly
		useful feature for catching cross-platform errors if we detect argument type errors before detecting the error of using a platform unsupported preference. 
  
*/


#include "Screen.h"

#if PSYCH_SYSTEM != PSYCH_WINDOWS

// for getpid()
#include <sys/types.h>
#include <unistd.h>

#else

// Windows doesn't have getpid() - We create a pseudo-getpid:
psych_uint64 getpid(void)
{
  HANDLE p;
  p=GetCurrentProcess();
  return((psych_uint64) p);
}
#endif

/*
OS 9 Psychtoolbox preferences:

	Supported
oldBool=Screen('Preference', 'IgnoreCase' [,psych_bool])
tick0Secs=Screen('Preference', 'Tick0Secs', tick0Secs)
psychTableVersion=Screen('Preference', 'PsychTableVersion')
mexFunctionName=Screen('Preference', 'PsychTableCreator')
proc=Screen('Preference', 'Process', signature)   //returns PID and does not accept the signature 


	Read supported but will not accept set:
oldBool=Screen('Preference','Backgrounding',[psych_bool])
oldSecondsMultiplier=Screen('Preference', 'SecondsMultiplier' [,secondsMultiplier])
	
	Unsupported and we should
oldBool=Screen('Preference' ',AllowMirroring' [,psych_bool])
psych_bool=Screen('Preference' ',MirroringIsOn')
oldBits=Screen(windowPtrOrScreenNumber, 'Preference', 'DacBits', bits) 

	Unsupported
available=Screen('Preference','Available')
nextProc=Screen('Preference','NextProcess',[proc])
oldBool=Screen('Preference','UseNewHideMenuBar',[psych_bool])

	Unsupported OS 9 Mac video driver specific 
psych_bool=Screen(windowPtrOrScreenNumber,'Preference','FixedClut')
oldBool=Screen(windowPtrOrScreenNumber,'Preference','SetClutDriverWaitsForBlanking',[psych_bool])
oldBool=Screen(windowPtrOrScreenNumber,'Preference','AskSetClutDriverToWaitForBlanking',[psych_bool])
oldValue=Screen(windowPtrOrScreenNumber,'Preference','SetClutSAI',[value])
oldBool=Screen(windowPtrOrScreenNumber,'Preference','SetClutCallsWaitBlanking',[psych_bool])
oldBool=Screen(windowPtrOrScreenNumber,'Preference','SetClutPunchesBlankingClock',[psych_bool])
oldBool=Screen(windowPtrOrScreenNumber,'Preference','InterruptPunchesBlankingClock',[psych_bool])
oldPriority=Screen(windowPtrOrScreenNumber,'Preference','MaxPriorityForBlankingInterrupt',[priority])
oldBool=Screen(windowPtrOrScreenNumber,'Preference','WaitBlankingAlwaysCallsSetClut',[psych_bool])
oldSecs=Screen(windowPtrOrScreenNumber,'Preference','BlankingDuration',[secs])
oldN=Screen(windowPtrOrScreenNumber,'Preference','MinimumEntriesForSetClutToWaitForBlanking',[n])
oldPriority=Screen(windowPtrOrScreenNumber,'Preference','MinimumSetClutPriority',[priority])
oldPriority=Screen(windowPtrOrScreenNumber,'Preference','MaximumSetClutPriority',[priority])
oldBool=Screen(windowPtrOrScreenNumber,'Preference','DipPriorityAfterSetClut',[psych_bool])
oldBool=Screen(windowPtrOrScreenNumber,'Preference','UsesHighGammaBits',[psych_bool])

*/

static char useString[] = "oldPreferenceValue = Screen('Preference', preferenceName, [newPreferenceValue])";
//                                                         0             1                2
static char synopsisString[] = 
	"Get or set a Psychtoolbox preference."
	"Preference settings are global - they affect all operations of a module until changed."
	""
	"\noldBool = Screen('Preference', 'IgnoreCase', [psych_bool]);"
	"\ntick0Secs = Screen('Preference', 'Tick0Secs', tick0Secs);"
	"\npsychTableVersion = Screen('Preference', 'PsychTableVersion');"
	"\nmexFunctionName = Screen('Preference', 'PsychTableCreator');"
	"\nproc = Screen('Preference', 'Process', signature);"
	"\nproc = Screen('Preference', 'DebugMakeTexture', enableDebugging);"
	"\noldEnableFlag = Screen('Preference', 'TextAlphaBlending', [enableFlag]);"
	"\noldSize = Screen('Preference', 'DefaultFontSize', [fontSize]);"
	"\noldStyleFlag = Screen('Preference', 'DefaultFontStyle', [styleFlag]);"
	"\noldfontName = Screen('Preference', 'DefaultFontName', [fontName]);"
	"\noldEnableFlag = Screen('Preference', 'DefaultTextYPositionIsBaseline', [enableFlag]);"
	"\noldEnableFlag = Screen('Preference', 'TextAntiAliasing', [enableFlag=-1 (System setting), 0 = Disable, 1 = Enable, 2 = EnableHighQuality]);"
	"\noldEnableFlag = Screen('Preference', 'TextRenderer', [enableFlag=0 (Default OS-specific [fast]), 1 = HighQ OS-specific]);"
	"\noldLocaleNameString = Screen('Preference', 'TextEncodingLocale', [newLocalenNameString]);"
	"\noldEnableFlag = Screen('Preference', 'SkipSyncTests', [enableFlag]);"
	"\n[maxStddev, minSamples, maxDeviation, maxDuration] = Screen('Preference', 'SyncTestSettings' [, maxStddev=0.001 secs][, minSamples=50][, maxDeviation=0.1][, maxDuration=5 secs]);"
	"\noldEnableFlag = Screen('Preference', 'FrameRectCorrection', [enableFlag=1]);"
	"\noldLevel = Screen('Preference', 'VisualDebugLevel', level);"
	"\n\nWorkaround flags to work around all kind of deficient drivers and hardware:\n"
	"See 'help ConserveVRAMSettings' for settings and their effect.\n"
	"\noldMode = Screen('Preference', 'ConserveVRAM', mode);"
	"\n\nActivate compatibility mode: Try to behave like the old MacOS-9 Psychtoolbox:"
	"\noldEnableFlag = Screen('Preference', 'EmulateOldPTB', [enableFlag]);\n"
	"\noldEnableFlags = Screen('Preference', 'Enable3DGraphics', [enableFlags]);"
	"\noldEnableFlag = Screen('Preference', 'SuppressAllWarnings', [enableFlag]);"
	"\noldMode = Screen('Preference', 'VBLTimestampingMode', [newmode]);"
	"\noldVTOTAL = Screen('Preference', 'VBLEndlineOverride' [, newVTOTAL]);"
	"\noldMode = Screen('Preference', 'DefaultVideocaptureEngine', [newmode (0=Quicktime-SequenceGrabbers, 1=LibDC1394-Firewire, 2=LibARVideo, 3=GStreamer)]);"
	"\noldMode = Screen('Preference', 'OverrideMultimediaEngine', [newmode (0=Legacy-Quicktime, 1=GStreamer)]);"
	"\noldLevel = Screen('Preference', 'WindowShieldingLevel', [newLevel (0 = Behind all other windows - 2000 = In front of all other windows, the default)]);"
	"\nresiduals = Screen('Preference', 'SynchronizeDisplays', syncMethod [, screenId]);"
	"\noldMappings = Screen('Preference', 'ScreenToHead', screenId [, newHeadId, newCrtcId][, rank=0]);"

	"\noldLevel = Screen('Preference', 'Verbosity' [,level]);";

static char seeAlsoString[] = "";	


#define kPsychNumUnsupportedMacVideoPreferences			14
static char *unsupportedMacVideoPreferenceNames[]={
	"FixedClut" 
	"SetClutDriverWaitsForBlanking",
	"AskSetClutDriverToWaitForBlanking",
	"SetClutSAI SetClutCallsWaitBlanking",
	"SetClutPunchesBlankingClock",
	"InterruptPunchesBlankingClock",
	"MaxPriorityForBlankingInterrupt",
	"WaitBlankingAlwaysCallsSetClut",
	"BlankingDuration",
	"MinimumEntriesForSetClutToWaitForBlanking",
	"MinimumSetClutPriority",
	"MaximumSetClutPriority",
	"DipPriorityAfterSetClut",
	"UsesHighGammaBits"
};

#define kPsychNumUnsupportedMacNonVideoPreferences			6
static char *unsupportedMacNonVideoPreferenceNames[]={
	"AllowMirroring",
	"MirroringIsOn",
	"DacBits",
	"Available",
	"NextProcess",
	"UseNewHideMenuBar"
};

PsychError SCREENPreference(void)  
{

	PsychArgFormatType		arg1Type;
	char					*preferenceName, *newFontName;
	const char				*tableCreator, *oldDefaultFontName;
	psych_bool				preferenceNameArgumentValid, booleanInput, ignoreCase, tempFlag, textAlphaBlendingFlag, suppressAllWarningsFlag;
	int						numInputArgs, i, newFontStyleNumber, newFontSize, tempInt, tempInt2, tempInt3, tempInt4;
	double					returnDoubleValue, inputDoubleValue;
	double					maxStddev, maxDeviation, maxDuration;
	int						minSamples;
    double                  *dheads = NULL;

	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check for superfluous or missing arguments
	PsychErrorExit(PsychCapNumInputArgs(5));			
	PsychErrorExit(PsychRequireNumInputArgs(1));		   
	PsychErrorExit(PsychCapNumOutputArgs(4));
	
	numInputArgs=PsychGetNumInputArgs();
	arg1Type=PsychGetArgType(1);
	preferenceNameArgumentValid=FALSE;
	
	//Cases which require both a window pointer or screen number and preference name.  Argument 1 is the wposn and argument 2 is the preference name.
	if( numInputArgs >= 2 && (PsychIsScreenNumberArg(1) || PsychIsScreenNumberArg(1)) && PsychGetArgType(2)==PsychArgType_char ){
		PsychAllocInCharArg(2, kPsychArgRequired, &preferenceName);
		//preferences which require window pointer or screen number argument which we DO NOT support  
		for(i=0;i<kPsychNumUnsupportedMacVideoPreferences;i++){
			if(PsychMatch(preferenceName, unsupportedMacVideoPreferenceNames[i]))
				PsychErrorExit(PsychError_unsupportedOS9Preference);
		}
		//insert here conditionals  to act on prefernces which accept a window pointer or screen number argument which we DO support.
		PsychErrorExit(PsychError_unrecognizedPreferenceName);
	}

	//Cases which do not require a wposn.  Argument 1 is the preference name.  if present Argument 2 is the new value
	if(arg1Type==PsychArgType_char){
		PsychAllocInCharArg(1, kPsychArgRequired, &preferenceName);
		//Preferernces which we do not support and which do not require a wposn
		for(i=0;i<kPsychNumUnsupportedMacNonVideoPreferences;i++){
			if(PsychMatch(preferenceName, unsupportedMacNonVideoPreferenceNames[i]))
				PsychErrorExit(PsychError_unsupportedOS9Preference);
		}
		//Preferences which we do support
		if(PsychMatch(preferenceName, "IgnoreCase")){
			ignoreCase=!PsychIsPsychMatchCaseSensitive();
			PsychCopyOutFlagArg(1, kPsychArgOptional, ignoreCase);
			if(numInputArgs==2){
				PsychCopyInFlagArg(2, kPsychArgRequired, &booleanInput);
				PsychSetPsychMatchCaseSenstive(!booleanInput);			
			}
			preferenceNameArgumentValid=TRUE;
		}else 
		if(PsychMatch(preferenceName, "Tick0Secs")){
			if(PsychCopyInDoubleArg(2, kPsychArgOptional, &inputDoubleValue) && inputDoubleValue==PsychGetNanValue())
				PsychEstimateGetSecsValueAtTickCountZero();
			returnDoubleValue=PsychGetEstimatedSecsValueAtTickCountZero();
			PsychCopyOutDoubleArg(1, kPsychArgOptional, returnDoubleValue);
			preferenceNameArgumentValid=TRUE;
		}else 
		if(PsychMatch(preferenceName, "PsychTableVersion")){
			if(numInputArgs==2)
				PsychErrorExit(PsychError_extraInputArg);
			PsychCopyOutDoubleArg(1, kPsychArgOptional, (double)PsychPrefStateGet_PsychTableVersion());
			preferenceNameArgumentValid=TRUE;
		}else 
		if(PsychMatch(preferenceName, "PsychTableCreator")){
			if(numInputArgs==2)
				PsychErrorExit(PsychError_extraInputArg);
			tableCreator=PsychPrefStateGet_PsychTableCreator();
			PsychCopyOutCharArg(1, kPsychArgOptional, tableCreator);
			preferenceNameArgumentValid=TRUE;
		}else 
		if(PsychMatch(preferenceName, "Process")){
			if(numInputArgs==2)
				PsychErrorExit(PsychError_extraInputArg);
			PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) (psych_int64) getpid());
			preferenceNameArgumentValid=TRUE;


		}else 
		if(PsychMatch(preferenceName, "DefaultFontName")){
			PsychPrefStateGet_DefaultFontName(&oldDefaultFontName);
			PsychCopyOutCharArg(1, kPsychArgOptional, oldDefaultFontName);
			if(numInputArgs==2){
				PsychAllocInCharArg(2, kPsychArgRequired, &newFontName);
				PsychPrefStateSet_DefaultFontName(newFontName);
			}
			preferenceNameArgumentValid=TRUE;			
		}else 
		if(PsychMatch(preferenceName, "TextEncodingLocale")){
			PsychCopyOutCharArg(1, kPsychArgOptional, PsychGetUnicodeTextConversionLocale());
			if(numInputArgs==2){
				PsychAllocInCharArg(2, kPsychArgRequired, &newFontName);
				if (!PsychSetUnicodeTextConversionLocale(newFontName)) PsychErrorExitMsg(PsychError_user, "Setting the 'TextEncodingLocale' failed, most likely because you provided an invalid/unknown locale setting string.");
			}
			preferenceNameArgumentValid=TRUE;			
		}else 
		if(PsychMatch(preferenceName, "DefaultFontStyle")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_DefaultTextStyle());
			if(numInputArgs==2){
				PsychCopyInIntegerArg(2, kPsychArgRequired, &newFontStyleNumber);
				PsychPrefStateSet_DefaultTextStyle(newFontStyleNumber);
			}
			preferenceNameArgumentValid=TRUE;
		}else
		if(PsychMatch(preferenceName, "OverrideMultimediaEngine")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_UseGStreamer());
			if(numInputArgs==2){
				PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
				PsychPrefStateSet_UseGStreamer(tempInt);
			}
			preferenceNameArgumentValid=TRUE;
		}else
		if(PsychMatch(preferenceName, "DefaultTextYPositionIsBaseline")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_TextYPositionIsBaseline());
			if(numInputArgs==2){
				PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
				PsychPrefStateSet_TextYPositionIsBaseline(tempInt);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
		if(PsychMatch(preferenceName, "TextAntiAliasing")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_TextAntiAliasing());
			if(numInputArgs==2){
				PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
				PsychPrefStateSet_TextAntiAliasing(tempInt);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
		if(PsychMatch(preferenceName, "TextRenderer")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_TextRenderer());
			if(numInputArgs==2){
				PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
				PsychPrefStateSet_TextRenderer(tempInt);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
		if(PsychMatch(preferenceName, "DefaultFontSize")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_DefaultTextSize());
			if(numInputArgs==2){
				PsychCopyInIntegerArg(2, kPsychArgRequired, &newFontSize);
				PsychPrefStateSet_DefaultTextSize(newFontSize);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
		if(PsychMatch(preferenceName, "DebugMakeTexture")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_DebugMakeTexture());
			if(numInputArgs==2){
				PsychCopyInFlagArg(2, kPsychArgRequired, &tempFlag);
				PsychPrefStateSet_DebugMakeTexture(tempFlag);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "SkipSyncTests")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_SkipSyncTests());
			if(numInputArgs==2){
				PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
				PsychPrefStateSet_SkipSyncTests(tempInt);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "VisualDebugLevel")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_VisualDebugLevel());
			if(numInputArgs==2){
				PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
				PsychPrefStateSet_VisualDebugLevel(tempInt);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "VBLTimestampingMode")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_VBLTimestampingMode());
			if(numInputArgs>=2){
                            PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
                            PsychPrefStateSet_VBLTimestampingMode(tempInt);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "SyncTestSettings")){
			PsychPrefStateGet_SynctestThresholds(&maxStddev, &minSamples, &maxDeviation, &maxDuration);
			PsychCopyOutDoubleArg(1, kPsychArgOptional, maxStddev);
			PsychCopyOutDoubleArg(2, kPsychArgOptional, minSamples);
			PsychCopyOutDoubleArg(3, kPsychArgOptional, maxDeviation);
			PsychCopyOutDoubleArg(4, kPsychArgOptional, maxDuration);
			if(numInputArgs>=2){
							PsychCopyInDoubleArg( 2, kPsychArgOptional, &maxStddev);
                            PsychCopyInIntegerArg(3, kPsychArgOptional, &minSamples);
							PsychCopyInDoubleArg( 4, kPsychArgOptional, &maxDeviation);
							PsychCopyInDoubleArg( 5, kPsychArgOptional, &maxDuration);
							PsychPrefStateSet_SynctestThresholds(maxStddev, minSamples, maxDeviation, maxDuration);
			}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "VBLEndlineOverride")){
			PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_VBLEndlineOverride());
			if(numInputArgs>=2){
                            PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
                            PsychPrefStateSet_VBLEndlineOverride(tempInt);
			}
			preferenceNameArgumentValid=TRUE;
		}else
			if(PsychMatch(preferenceName, "DefaultVideocaptureEngine")){
				PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_VideoCaptureEngine());
				if(numInputArgs==2){
					PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
					PsychPrefStateSet_VideoCaptureEngine(tempInt);
				}
			preferenceNameArgumentValid=TRUE;
		}else
			if(PsychMatch(preferenceName, "WindowShieldingLevel")){
				PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_WindowShieldingLevel());
				if(numInputArgs==2){
					PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
					PsychPrefStateSet_WindowShieldingLevel(tempInt);
				}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "ConserveVRAM") || PsychMatch(preferenceName, "Workarounds1")){
					PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_ConserveVRAM());
					if(numInputArgs==2){
						PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
						PsychPrefStateSet_ConserveVRAM(tempInt);
					}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "Verbosity")){
					PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_Verbosity());
					if(numInputArgs==2){
						PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
						PsychPrefStateSet_Verbosity(tempInt);
					}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "FrameRectCorrection")){
					PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_FrameRectCorrection());
					if(numInputArgs==2){
						PsychCopyInDoubleArg(2, kPsychArgRequired, &inputDoubleValue);
						PsychPrefStateSet_FrameRectCorrection(inputDoubleValue);
					}
			preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "EmulateOldPTB")){
				PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_EmulateOldPTB());
				if(numInputArgs==2){
					PsychCopyInFlagArg(2, kPsychArgRequired, &tempFlag);
					PsychPrefStateSet_EmulateOldPTB(tempFlag);
				}
				preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "Enable3DGraphics")){
				PsychCopyOutDoubleArg(1, kPsychArgOptional, PsychPrefStateGet_3DGfx());
				if(numInputArgs==2){
					PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
					PsychPrefStateSet_3DGfx(tempInt);
				}
				preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "TextAlphaBlending")){
				textAlphaBlendingFlag=PsychPrefStateGet_TextAlphaBlending();
				PsychCopyOutFlagArg(1, kPsychArgOptional, textAlphaBlendingFlag);
				if(numInputArgs==2){
					PsychCopyInFlagArg(2, kPsychArgRequired, &booleanInput);
					PsychPrefStateSet_TextAlphaBlending(booleanInput);
				}
				preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "SuppressAllWarnings")){
				suppressAllWarningsFlag=PsychPrefStateGet_SuppressAllWarnings();
				PsychCopyOutFlagArg(1, kPsychArgOptional, suppressAllWarningsFlag);
				if(numInputArgs==2){
					PsychCopyInFlagArg(2, kPsychArgRequired, &booleanInput);
					PsychPrefStateSet_SuppressAllWarnings(booleanInput);
				}
				preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "SynchronizeDisplays")){
				if(numInputArgs >= 2) {
					// This is a special call: It currently doesn't set a preference setting,
					// but instead triggers an instantaneous synchronization of all available
					// display heads, if possible. We may have a more clever and "standard" interface
					// interface for this later on, but for first tests this will do.
					// Syncmethod is hard-coded to 0 -> Use whatever's available to sync.
					// timeout for retries is 5.0 seconds.
					// Acceptable residual offset is +/- 2 scanlines.
					// Returns the real residual offset after sync.
					PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
					if (!PsychCopyInIntegerArg(3, kPsychArgOptional, &tempInt3)) {
                        // No screenId specified: Resync default screen or whatever...
                        tempInt2 = 0;
                        if (PsychSynchronizeDisplayScreens(&tempInt2, NULL, &tempInt, tempInt, 5.0, 2)!=PsychError_none) PsychErrorExitMsg(PsychError_user, "Sync failed for reasons mentioned above.");
                    } else {
                        // Specific screenId provided: Resync crtc's associated with this screenId if possible:
                        tempInt2 = 1;
                        if (PsychSynchronizeDisplayScreens(&tempInt2, &tempInt3, &tempInt, tempInt, 5.0, 2)!=PsychError_none) PsychErrorExitMsg(PsychError_user, "Sync failed for reasons mentioned above.");
                    }
                    
					PsychCopyOutDoubleArg(1, kPsychArgOptional, tempInt);
				}
				preferenceNameArgumentValid=TRUE;
		}else 
			if(PsychMatch(preferenceName, "ScreenToHead")){
				// screenId is required:
				PsychCopyInIntegerArg(2, kPsychArgRequired, &tempInt);
				if (tempInt < 0 || tempInt >= PsychGetNumDisplays() || tempInt >= kPsychMaxPossibleDisplays) PsychErrorExitMsg(PsychError_user, "Invalid screenId provided. Out of valid range!");

				// Return old mappings for this screenId:
                for (tempInt2 = 0; (tempInt2 < kPsychMaxPossibleCrtcs) && (PsychPrefStateGet_ScreenToHead(tempInt, tempInt2) >= 0); tempInt2++);
                PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 2, tempInt2, 1, &dheads);
                
                tempInt4 = 0;
                for (tempInt3 = 0; tempInt3 < tempInt2; tempInt3++) {
                    dheads[tempInt4++] = (double) PsychPrefStateGet_ScreenToHead(tempInt, tempInt3);
                    dheads[tempInt4++] = (double) PsychPrefStateGet_ScreenToCrtcId(tempInt, tempInt3);
                }
                
                // Optionally retrieve and set new mappings for this screenId:
				if(numInputArgs>=3) {
					// Set new headId for screenId:
					PsychCopyInIntegerArg(3, kPsychArgRequired, &tempInt2);
					if (tempInt2 < 0) PsychErrorExitMsg(PsychError_user, "Invalid negative headId provided!");

					// Set new crtcId for screenId:
					PsychCopyInIntegerArg(4, kPsychArgRequired, &tempInt3);
					if (tempInt3 < 0) PsychErrorExitMsg(PsychError_user, "Invalid negative crtcId provided!");

                    // Assign primary head by default (index 0), but allow optionally others as well:
                    tempInt4 = 0;
					PsychCopyInIntegerArg(5, kPsychArgOptional, &tempInt4);
					if (tempInt4 < 0 || tempInt4 >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid rankId provided! Too many heads for one screen!");

					PsychPrefStateSet_ScreenToHead(tempInt, tempInt2, tempInt3, tempInt4);
				}
				preferenceNameArgumentValid=TRUE;
		}else 
			PsychErrorExit(PsychError_unrecognizedPreferenceName);
	}
	
	if(!preferenceNameArgumentValid)
		PsychErrorExitMsg(PsychError_user, "Invalid arguments to preferences command");
		
	return(PsychError_none);
}
