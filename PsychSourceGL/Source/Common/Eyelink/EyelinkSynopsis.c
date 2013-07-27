/*

	PsychSourceGL/Source/OSX/Eyelink/EyelinkSynopsis.c

	PROJECTS: Eyelink 
  
	AUTHORS:
		cburns@berkeley.edu				cdb
		E.Peters@ai.rug.nl				emp
		f.w.cornelissen@rug.nl          fwc

	PLATFORMS:	All.

	HISTORY:

		11/21/05  cdb		Created based on old EyeMex.c and ScreenSynopsis.c.
		15/06/06  fwc		Added few functions and small changes to synopsis.

	TARGET LOCATION:

		Eyelink.mexmac resides in:
			EyelinkToolbox
*/

#include "PsychEyelink.h"

#define MAX_SYNOPSIS_STRINGS 500  

//declare variables local to this file.  
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

void InitializeSynopsis()
{
	int i=0;
	const char **synopsis = synopsisSYNOPSIS;  //abbreviate the long name

	synopsis[i++] = "\n% This is the main function of the Eyelink toolbox";
	synopsis[i++] = "Usage:";
    synopsis[i++] = "\n% For general advice, try:";
    synopsis[i++] = "help Eyelink";
    synopsis[i++] = "\n% For a more detailed explanation of any Eyelink function, just add a question mark \"?\".";
    synopsis[i++] = "% E.g. for 'Initialize', try either of these equivalent forms:";
    synopsis[i++] = "Eyelink('Initialize?')";
    synopsis[i++] = "Eyelink Initialize?";
    synopsis[i++] = "\n% If you think you've found a bug, please report it";
    synopsis[i++] = "on the forum, see: http://psychtoolbox.org/\n";
    
	// Init or close eyelink
	synopsis[i++] = "\n% Initialize or shutdown Eyelink connection:";
	synopsis[i++] = "[status =] Eyelink('Initialize' [, displayCallbackFunction])";
	synopsis[i++] = "[status =] Eyelink('InitializeDummy' [, displayCallbackFunction])";
	synopsis[i++] = "[status =] Eyelink('IsConnected')";
	synopsis[i++] = "[status =] Eyelink('SetAddress', ipaddress);";
	synopsis[i++] = "Eyelink('Shutdown')";
	synopsis[i++] = "oldlevel = Eyelink('Verbosity' [,level]);";
	synopsis[i++] = "Eyelink('TestSuite')";
	synopsis[i++] = "[status =] Eyelink('OpenFile','filename')";
	synopsis[i++] = "[status =] Eyelink('CloseFile')";
	synopsis[i++] = "[status =] Eyelink('ReceiveFile',['filename'], ['dest'], ['dest_is_path'])";
	
	// Calibration
	synopsis[i++] = "\n% Calibration:";
	synopsis[i++] = "[result =] Eyelink('StartSetup' [, stype=0])";
	synopsis[i++] = "[status = ] Eyelink('DriftCorrStart', x, y [,dtype=0][, dodraw=1][, allow_setup=0])";
	synopsis[i++] = "[result = ] Eyelink('ApplyDriftCorr')";
	synopsis[i++] = "[result, tx, ty] = Eyelink('TargetCheck')";
	synopsis[i++] = "[result = ] Eyelink('AcceptTrigger')";
	synopsis[i++] = "[result, messageString =] Eyelink('CalMessage')";
	
	// Start or stop recording, Data acquisition
	synopsis[i++] = "\n% Start or stop recording and acquiring data:";
	synopsis[i++] = "[startrecording_error =] Eyelink('StartRecording' [,file_samples, file_events, link_samples, link_events] )";
	synopsis[i++] = "Eyelink('Stoprecording')";
	synopsis[i++] = "error = Eyelink('CheckRecording')";
	synopsis[i++] =  "eyeused = Eyelink('EyeAvailable')";
	synopsis[i++] = "NewOrOld = Eyelink('NewFloatSampleAvailable')";
	synopsis[i++] = "sample = Eyelink('NewestFloatSample')";
	synopsis[i++] = "[sample, raw] = Eyelink('NewestFloatSampleRaw' [, eye])";
	synopsis[i++] = "type = Eyelink('GetNextDataType')";
	synopsis[i++]  = "item = Eyelink('GetFloatData', type)";
	synopsis[i++]  = "[item, raw] = Eyelink('GetFloatDataRaw', type [, eye])";
	synopsis[i++]  = "[samples, events, drained] = Eyelink('GetQueuedData'[, eye])";
    
	// Misc eyelink communication:
	synopsis[i++] = "\n% Miscellaneous functions to communicate with Eyelink:";
	synopsis[i++] = "result = Eyelink('ButtonStates')";
	synopsis[i++] = "[status =] Eyelink('Command', 'formatstring', [...])";
	synopsis[i++] = "[status =] Eyelink('Message', 'formatstring', [...])";
	synopsis[i++] = "[result =] Eyelink('SendKeyButton', code, mods, state)";
	synopsis[i++] = "[time =] Eyelink('TrackerTime')";
	synopsis[i++] = "[offset =] Eyelink('TimeOffset')";
	synopsis[i++] = "[status =] Eyelink('RequestTime')";
	synopsis[i++] = "[time =] Eyelink('ReadTime')";
	synopsis[i++] = "[mode =] Eyelink('TrackerMode')";
	synopsis[i++] = "[result, reply =] Eyelink('ReadFromTracker', VariableName)";

	synopsis[i++] = "\n% Miscellaneous Eyelink functions:";
	synopsis[i++] = "[result =] Eyelink('WaitForModeReady', maxwait)";
	synopsis[i++] = "[result =] Eyelink('ImageModeDisplay')";
	synopsis[i++] = "mode = Eyelink('CurrentMode')";
	synopsis[i++] = "result = Eyelink('CalResult')";
	synopsis[i++] = "Eyelink('SetOfflineMode')";
	synopsis[i++] = "[version, versionString]  = Eyelink('GetTrackerVersion')";
	synopsis[i++] = "[time =] Eyelink('TrackerTime')";
	synopsis[i++] = "[offset =] Eyelink('TimeOffset')";
	synopsis[i++] = "[status] = Eyelink('ImageTransfer', imagePath [, xPosition=0][, yPosition=0][, width=0][, height=0][, trackerXPosition=0][, trackerYPosition=0][, xferoptions=0])";

	// Place Holder
	synopsis[i++] = "\n\n\n\n% EyelinkToolbox version for the OpenGL PsychToolbox";
	synopsis[i++] = "% The EyelinkToolbox was developed by:\n";
	synopsis[i++] = "\tFrans Cornelissen";
	synopsis[i++] = "\tEnno Peters";
	synopsis[i++] = "\tJohn Palmer";
	synopsis[i++] = "\tChris Burns";
	synopsis[i++] = "\tMario Kleiner";
	synopsis[i++] = "\tErik Flister";
	synopsis[i++] = "\tNuha Jabakhanji";
	
	synopsis[i++] = NULL;  //this tells PsychDisplayScreenSynopsis where to stop

	if (i > MAX_SYNOPSIS_STRINGS) {
		PrintfExit("%s: increase dimension of synopsis[] from %ld to at least %ld and recompile.", __FILE__, (long) MAX_SYNOPSIS_STRINGS, (long) i);
	}
}

PsychError PsychDisplayEyelinkSynopsis(void)
{
	int i;
	
	for (i = 0; synopsisSYNOPSIS[i] != NULL; i++) {
		printf("%s\n", synopsisSYNOPSIS[i]);
	}
		
	return(PsychError_none);
}
