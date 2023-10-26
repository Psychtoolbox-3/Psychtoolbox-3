/*

	/osxptb/trunk/PsychSourceGL/Source/OSX/Eyelink/EyelinkImageModeDisplay.c
  
	PROJECTS: Eyelink 
  
	AUTHORS:
		cburns@berkeley.edu				cdb
		E.Peters@ai.rug.nl				emp
		f.w.cornelissen@med.rug.nl		fwc
  
	PLATFORMS:	Currently only OS X  
    
	HISTORY:

		11/23/05  cdb		Created.
		27/03/09  edf       added waiting loop and check that we made it into the mode

	TARGET LOCATION:

		Eyelink.mexmac resides in:
			EyelinkToolbox
*/

#include "PsychEyelink.h"

#define WAIT_MS 10
#define MAX_LOOPS 1000
#define ERR_BUFF_LEN 1000

static char useString[] = "[result =] Eyelink('ImageModeDisplay')";

static char synopsisString[] =
	"This handles display of the EyeLink camera images. "
	"While in imaging mode, it continuously requests "
	"and displays the current camera image. "
	"It also displays the camera name and threshold setting. "
	"Keys on the subject PC keyboard are sent to the tracker "
	"so the experimenter can use it during setup. "
	"It will exit when the tracker leaves "
	"imaging mode or disconnects. "
    "RETURNS: 0 if OK, 32767 (Ox7FFF or TERMINATE_KEY) if pressed, -1 if disconnect";

static char seeAlsoString[] = "";

/*
ROUTINE: EYELINKimagemodedisplay
PURPOSE:
   uses INT16 image_mode_display(void);

	This handles display of the EyeLink camera images
	While in imaging mode, it contiuously requests
	and displays the current camera image
	It also displays the camera name and threshold setting
	Keys on the subject PC keyboard are sent to the tracker
	so the experimenter can use it during setup.
	It will exit when the tracker leaves
	imaging mode or discannects

    RETURNS: 0 if OK, TERMINATE_KEY if pressed, -1 if disconnect*/

PsychError EyelinkImageModeDisplay(void)
{
	int		iResult		= 0;
	
	// Add help strings
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	
	// Output help if asked
	if(PsychIsGiveHelp()) {
		PsychGiveHelp();
		return(PsychError_none);
	}
		
	// Check arguments
	PsychErrorExit(PsychCapNumInputArgs(0));
	PsychErrorExit(PsychRequireNumInputArgs(0));
	PsychErrorExit(PsychCapNumOutputArgs(1));
	
	// Verify eyelink is up and running
	EyelinkSystemIsConnected();
	EyelinkSystemIsInitialized();
	
	// Optionally dump the whole hookfunctions struct:
	if (Verbosity() > 5) { printf("Eyelink-Debug: ImageModeDisplay: PreOp: \n"); PsychEyelink_dumpHookfunctions(); }
	
	iResult = image_mode_display();
	
	// Optionally dump the whole hookfunctions struct:
	if (Verbosity() > 5) { printf("Eyelink-Debug: ImageModeDisplay: PostOp: \n"); PsychEyelink_dumpHookfunctions(); }

	// Assign output arg if available
	PsychCopyOutDoubleArg(1, FALSE, iResult);
	
	return(PsychError_none);
}
