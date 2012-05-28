/*
	Psychtoolbox3/Source/Common/SCREENReadNormalizedGammaTable.c		

	AUTHORS:

	Allen.Ingling@nyu.edu		awi 
	mario.kleiner@tuebingen.mpg.de  mk

	PLATFORMS:

	All.

	HISTORY:

	1/27/03		awi		Created. 

	DESCRIPTION:

	TO DO:
*/

#include "Screen.h"

#if PSYCH_SYSTEM == PSYCH_OSX
// Disable warnings about deprecated API calls on OSX 10.7
// of which we are aware and that we can't remove as long as
// we need to stay compatible to 10.4 - 10.6
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "[gammatable, dacbits, reallutsize] = Screen('ReadNormalizedGammaTable', windowPtrOrScreenNumber [, physicalDisplay]);";
static char synopsisString[] = 
		"Reads and returns the gamma table 'gammatable' of the specified screen or window 'windowPtrOrScreenNumber'.\n"
		"Returns the output resolution "
		"of the video DAC as optional second argument 'dacbits'. Will return dacbits=8 as a safe default if it "
		"is unable to query the real resolution of the DAC. Currently only OS-X reports the real DAC size. "
		"Will return the real number of slots in the hardware lookup table in optional return argument 'reallutsize'. "
		"Currently only OS-X and Linux report the real LUT size.\n"
		"On MacOS-X, the optional 'physicalDisplay' flag can be set to 1, zero is the default. In this case, "
		"the 'windowPtrOrScreenNumber' argument (which then must be a real screen number, not a window index) "
		"selects among physically present display devices, instead of logical devices. "
		"This is important if you want to assign different gamma-tables to multiple displays in a 'clone' or "
		"'mirror mode' configuration, as there is only one logical display, but multiple physical displays, mirroring "
		"each other. Please note that screen numbering is different for physical vs. logical displays. For a list of "
		"physical display indices, call Screen('Screens', 1);\n"
		"On GNU/Linux, the optional 'physicalDisplay' parameter selects the video output from which the gamma "
		"table should be read in multi-display mode. On Linux a screen can output to multiple video displays, "
		"therefore this parameter allows to query the individual gamma tables for each display. The default setting "
		"is -1, which means to read the gamma table of the primary output of the given screen.\n"
		"See help for Screen('LoadNormalizedGammaTable'); for infos about the format of the returned table "
		"and for further explanations regarding gamma tables.";

static char seeAlsoString[] = "";

PsychError SCREENReadNormalizedGammaTable(void)
{
    int		i, screenNumber, numEntries, reallutsize, physicalDisplay, outputId;
    float 	*redTable, *greenTable, *blueTable;
    double	*gammaTable;	

    //all subfunctions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(2));

    // Get optional physicalDisplay argument - It defaults to zero:
    physicalDisplay = -1;
    PsychCopyInIntegerArg(2, FALSE, &physicalDisplay);

    // Read in the screen number:
    // On OS/X we also accept screen indices for physical displays (as opposed to active dispays).
    // This only makes a difference in mirror-mode, where there is only 1 active display, but that
    // corresponds to two physical displays which can have different gamma setting requirements:
    if ((PSYCH_SYSTEM == PSYCH_OSX) && (physicalDisplay > 0)) {
        PsychCopyInIntegerArg(1, TRUE, &screenNumber);
	if (screenNumber < 1) PsychErrorExitMsg(PsychError_user, "A 'screenNumber' that is smaller than one provided, although 'physicalDisplay' flag set. This is not allowed!");

	// Invert screenNumber as a sign its a physical display, not an active display:
	screenNumber = -1 * screenNumber;
    }
    else {
        PsychCopyInScreenNumberArg(1, TRUE, &screenNumber);
    }

    if ((PSYCH_SYSTEM == PSYCH_LINUX) && (physicalDisplay > -1)) {
	// Affect one specific display output for given screen:
	outputId = physicalDisplay;
    }
    else {
	// Other OS'es, and Linux with default setting: Affect all outputs
	// for a screen.
	outputId = -1;
    }

    // Retrieve gamma table:
    PsychReadNormalizedGammaTable(screenNumber, outputId, &numEntries, &redTable, &greenTable, &blueTable);
	
    // Copy it out to runtime:
    PsychAllocOutDoubleMatArg(1, FALSE, numEntries, 3, 0, &gammaTable);

    for(i=0;i<numEntries;i++){
        gammaTable[PsychIndexElementFrom3DArray(numEntries, 3, 0, i, 0, 0)]=(double)redTable[i];
        gammaTable[PsychIndexElementFrom3DArray(numEntries, 3, 0, i, 1, 0)]=(double)greenTable[i];
        gammaTable[PsychIndexElementFrom3DArray(numEntries, 3, 0, i, 2, 0)]=(double)blueTable[i];
    }

    // Copy out optional DAC resolution value:
    PsychCopyOutDoubleArg(2, FALSE, (double) PsychGetDacBitsFromDisplay(screenNumber));
	
    // We default to the assumption that the real size of the hardware LUT is identical to
    // the size of the returned LUT:
    reallutsize = numEntries;
	
    #if PSYCH_SYSTEM == PSYCH_OSX
		// On OS-X we query the real LUT size from the OS and return that value:
		CGDirectDisplayID	displayID;
		CFMutableDictionaryRef properties;
		CFNumberRef cfGammaLength;
		SInt32 lutslotcount;
		io_service_t displayService;
		kern_return_t kr;
		CFMutableArrayRef framebufferTimings0 = 0;
		CFDataRef framebufferTimings1 = 0;
		IODetailedTimingInformationV2 *framebufferTiming = NULL;
		
		// Retrieve display handle for screen:
		PsychGetCGDisplayIDFromScreenNumber(&displayID, screenNumber);
		
		if (PsychPrefStateGet_Verbosity()>5) printf("PTB-DEBUG: Screen %i has framebuffer address %p.\n", screenNumber, CGDisplayBaseAddress(displayID));

		// Retrieve low-level IOKit service port for this display:
		displayService = CGDisplayIOServicePort(displayID);
				
		// Obtain the properties from that service
		kr = IORegistryEntryCreateCFProperties(displayService, &properties, NULL, 0);
		if((kr == kIOReturnSuccess) && ((cfGammaLength = (CFNumberRef) CFDictionaryGetValue(properties, CFSTR(kIOFBGammaCountKey)))!=NULL))
		{
			CFNumberGetValue(cfGammaLength, kCFNumberSInt32Type, &lutslotcount);
			CFRelease(properties);
			reallutsize = (int) lutslotcount;
		}
		else {
			// Failed!
			if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Failed to query real size of video LUT for screen %i! Will return safe default of %i slots.\n", screenNumber, reallutsize);
		}	

		if (PsychPrefStateGet_Verbosity()>9) {			
			CFDictionaryRef currentMode;
			CFNumberRef n;
			int modeId;
			currentMode = CGDisplayCurrentMode(displayID);
			n=CFDictionaryGetValue(currentMode, kCGDisplayMode);
			CFNumberGetValue(n, kCFNumberIntType, &modeId);
			printf("Current mode has id %i\n\n", modeId);
			kr = IORegistryEntryCreateCFProperties(displayService, &properties, NULL, 0);
			if((kr == kIOReturnSuccess) && ((framebufferTimings0 = (CFMutableArrayRef) CFDictionaryGetValue(properties, CFSTR(kIOFBDetailedTimingsKey) ) )!=NULL))
			{
				for (i=0; i<CFArrayGetCount(framebufferTimings0); i++) {
					if ((framebufferTimings1 = CFArrayGetValueAtIndex(framebufferTimings0, i)) != NULL) {
						if ((framebufferTiming = (IODetailedTimingInformationV2*) CFDataGetBytePtr(framebufferTimings1)) != NULL) {
							printf("[%i] : VActive =  %li, VBL = %li, VSYNC = %li, VSYNCWIDTH = %li , VBORDERBOT = %li, VTOTAL = %li \n", i, framebufferTiming->verticalActive, framebufferTiming->verticalBlanking, framebufferTiming->verticalSyncOffset, framebufferTiming->verticalSyncPulseWidth, framebufferTiming->verticalBorderBottom, framebufferTiming->verticalActive + framebufferTiming->verticalBlanking);
						}
					}
				}

				CFRelease(properties);
			}
			else {
				// Failed!
				if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Failed to query STUFF for screen %i --> %p!\n", screenNumber, properties);
			}	
		}
		


    #endif
	
    // Copy out optional real LUT size (number of slots):
    PsychCopyOutDoubleArg(3, FALSE, (double) reallutsize);

    return(PsychError_none);
}
