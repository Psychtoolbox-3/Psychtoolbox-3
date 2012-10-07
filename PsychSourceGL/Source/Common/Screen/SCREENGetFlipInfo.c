/*
  SCREENGetFlipInfo.c		
  
  AUTHORS:

  mario.kleiner at tuebingen.mpg.de  mk
  
  PLATFORMS:	All

  HISTORY:

  5.09.2011	mk		Created.
 
  DESCRIPTION:
  
  Returns information about execution of flips of a specific window.
  
  NOTES:

  Be careful with length of struct field names! Only names up to 31 characters are
  supported by Matlab 5.x (and maybe 6.x -- untested). Larger names cause matching
  failure!

*/

#include "Screen.h"

static char useString[] = "info = Screen('GetFlipInfo', windowPtr [, infoType=0] [, auxArg1]);";
static char synopsisString[] = 
	"Returns a struct with miscellaneous info about finished flips on the specified onscreen window.\n"
	"\n"
	"This function is currently only supported on Linux with the free graphics drivers.\n\n"
	"The function allows you to enable logging of timestamps and other status information "
	"about all completed bufferswaps, as triggered via Screen('Flip'), Screen('AsyncFlipBegin') etc. "
	"Whenever a flip completes, a little info struct is stored in a internal queue. This function "
	"allows you to fetch these info structs from the queue. It allows you to enable and disable "
	"logging of these info structs. Logging is disabled by default.\n\n"
	"\"windowPtr\" is the handle of the onscreen window for which info should be returned.\n"
	"\"infoType\" If left out or set to zero, the flip handle for the last invocation of flip is returned. "
	"E.g., if a value of 53 is returned, then the info struct with a info.SwapbuffersCount field of 53 "
	"will contain the info for that last invocation of flip.\n"
	"This allows to associate specific flips with the returned logged timestamps and other info.\n"
	"If set to 1, logging of flip completion info is enabled.\n"
	"If set to 2, logging of flip completion info is disabled.\n"
	"If set to 3, the oldest stored flip completion info is returned in a struct 'info'.\n\n"
	"The info struct contains the following fields:\n"
	"----------------------------------------------\n\n"
	"OnsetTime: Visual stimulus onset time of the completed Screen('Flip') operation.\n"
	"OnsetVBLCount: Video refresh cycle count when the flip completed.\n"
	"SwapbuffersCount: Serial number of this info struct. Corresponds to the handle returned for 'infoType' zero.\n"
	"SwapType: How was the flip executed? Low level info about strategy chosen by GPU.\n"
	"Note: Currently only PAGEFLIP flips are considered to have reliable timing and trustworthy timestamps!\n\n";

static char seeAlsoString[] = "OpenWindow, Flip, NominalFrameRate";

PsychError SCREENGetFlipInfo(void) 
{
#if PSYCH_SYSTEM == PSYCH_LINUX
	PsychWindowRecordType *windowRecord;
	int infoType = 0, retIntArg;
	double auxArg1, auxArg2, auxArg3;
#endif

	// All subfunctions should have these two lines.  
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

#if PSYCH_SYSTEM == PSYCH_LINUX

	PsychErrorExit(PsychCapNumInputArgs(2));     //The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1)); //The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(1));    //The maximum number of outputs

	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
	if (!PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_user, "Invalid 'windowPtr' specified. Not an onscreen window!");

	// Query infoType flag: Defaults to zero.
	PsychCopyInIntegerArg(2, FALSE, &infoType);
	if (infoType < 0 || infoType > 3) PsychErrorExitMsg(PsychError_user, "Invalid 'infoType' argument specified! Valid are 0, 1, 2, 3.");

	// Type 0: Return SBC handle of last scheduled flip:
	if (infoType == 0) {
		PsychCopyOutDoubleArg(1, FALSE, (double) windowRecord->target_sbc);
		return(PsychError_none);
	}

	// Type 1: Enable logging of swap completion events:
	if (infoType == 1) {
		PsychOSSwapCompletionLogging(windowRecord, 1, 0);
		return(PsychError_none);
	}

	// Type 2: Disable logging of swap completion events:
	if (infoType == 2) {
		PsychOSSwapCompletionLogging(windowRecord, 0, 0);
		return(PsychError_none);
	}

	// Type 3: Fetch logged swap completion events:
	if (infoType == 3) {
		PsychOSSwapCompletionLogging(windowRecord, 3, 1);
		return(PsychError_none);
	}
#else
	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, this function is only supported on Linux.");
#endif

	// Done.
	return(PsychError_none);
}
