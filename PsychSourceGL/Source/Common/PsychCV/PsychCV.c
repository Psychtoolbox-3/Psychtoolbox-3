/*
	PsychToolbox3/Source/Common/PsychCV/PsychCV.c
	
	PLATFORMS:	All

	AUTHORS:
	
	Mario Kleiner   mk      mario.kleiner at tuebingen.mpg.de
	
	HISTORY:
	
	5.01.08			mk		wrote it.  
	
	DESCRIPTION:
	
	A container for all kind of miscellaneous C or high-perf OpenGL routines somehow
	loosely related to computer vision and image processing tasks.
	
	Usually user code won't call this module directly, but it will call high level
	M-File code which in turn "outsources" parts of computationally very expensive
	tasks to this module.
	
	The initial implementation contains parts of the OpenEyes C-based computer vision
	eyetracker: The original pure C based, pure Linux implementation was split up in
	parts that can be handled by Psychtoolbox Screen() (image acquisition and visualization),
	high-level M file code for general flow control, and compute-intense C routines and
	routines that utilize OpenCV - these parts are inside PsychCV.
	
	The initial implementation will also contain other interesting bits of computer vision
	code taken from OpenCV.

	IMPORTANT:
	
	PsychCV will be built with CVEyeTracker / OpenEyes support and OpenCV library support
	enabled if the #define PSYCHCV_USE_OPENCV 1 is defined somewhere in the build scripts!
	
	By default (if this is omitted), it will build without OpenCV / OpenEyes support.
	The reason is that support of OpenCV incurs a large number of library dependencies which
	make installation quite a hazzle.
	
	As long as OpenEyes support is not in a enduser useable state anyway, we can save
	ourselves and the users from this dependency hell.

*/


#include "PsychCV.h"

// Declare variables local to this file.  

#ifndef PI
#define PI 3.141592654
#endif

// Help/Usage strings:
#define MAX_SYNOPSIS_STRINGS 50  
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

// Level of verbosity:
unsigned int  verbosity = 4;

// Status: Initialized or not?
static psych_bool psychCVInitialized = FALSE;

// Extern shutdown function for AR toolkit, defined in PsychCVARToolkit.hc
extern void PsychCVARExit(void);


/* Logger callback function to output OpenCV et al. debug messages at 'verbosity' > 5. */
void PsychCVLogger(const char* msg)
{
	if (verbosity > 5) printf("PTB-DEBUG: OpenCV says: %s", msg);
	return;
}

// Setup help/usage strings:
void InitializeSynopsis()
{
	int i=0;
	const char **synopsis = synopsisSYNOPSIS;  //abbreviate the long name
	
	synopsis[i++] = "PsychCV - Helper module for miscellaneous stuff related to OpenCV and/or computer vision:\n\n";
	synopsis[i++] = "This module contains a large amount of 3rd party software, so here are the credits for those parts:\n";
	#ifdef PSYCHCV_USE_OPENCV
	synopsis[i++] = "OpenEyes eye tracking is based mostly on a modified version of the OpenEyes toolkit,";
	synopsis[i++] = "http://hcvl.hci.iastate.edu/openEyes which was written by Dongheng Li, Derrick Parkhurst,";
	synopsis[i++] = "Jason Babcock and David Winfield.";
	synopsis[i++] = "OpenEyes toolkit is licensed to you under GPL v2.\n";
	synopsis[i++] = "PsychCV is also dependent on the OpenCV Open computer vision library: http://sourceforge.net/projects/opencvlibrary/ \n";
	synopsis[i++] = "which is licensed under the BSD license.\n";
	#endif
	synopsis[i++] = "The 'ARxxx' subfunctions are based on the ARToolkit: http://www.hitl.washington.edu/artoolkit/";
	synopsis[i++] = "ARToolkit is licensed to you under GPL v2.\n";
	synopsis[i++] = "\nGeneral information and settings:\n";
	synopsis[i++] = "version = PsychCV('Version');";
	synopsis[i++] = "oldlevel = PsychCV('Verbosity' [,level]);";
	synopsis[i++] = "\nHelper functions for memory buffer copies:\n";
	synopsis[i++] = "PsychCV('CopyMatrixToMemBuffer', matrix, memBufferPtr);";
	#ifdef PSYCHCV_USE_OPENCV
	synopsis[i++] = "\nSupport for the OpenEyes computer vision based eye tracker:\n";
	synopsis[i++] = "[EyeImageMemBuffer, EyeColorImageMemBuffer, SceneImageMemBuffer, ThresholdImageMemBuffer, EllipseImageMemBuffer] = PsychCV('OpenEyesInitialize', handle [, eyeChannels] [, eyeWidth][, eyeHeight][, sceneWidth][, sceneHeight][, logfilename]);";
	synopsis[i++] = "PsychCV('OpenEyesShutdown', handle);";
 	synopsis[i++] = "[oldSettings, ...] = PsychCV('OpenEyesParameters', handle [, pupilEdgeThreshold][, starburstRays][, minFeatureCandidates][, corneaWindowSize][, edgeThreshold][, gaussWidth][, maxPupilEccentricity] [, initialAngleSpread] [, fanoutAngle1] [, fanoutAngle2] [, featuresPerRay] [, specialFlags]);";
	synopsis[i++] = "EyeResult = PsychCV('OpenEyesTrackEyePosition', handle [, mode] [, px], [, py]);";
	#endif
	synopsis[i++] = "\nSupport for the ARToolkit computer vision based 3D marker tracking library:\n";
	synopsis[i++] = "[SceneImageMemBuffer, glProjectionMatrix, DebugImageMemBuffer] = PsychCV('ARInitialize', cameraCalibFilename, imgWidth, imgHeight, imgChannels [, imgFormat]);";
	synopsis[i++] = "PsychCV('ARShutdown');";
	synopsis[i++] = "[markerId] = PsychCV('ARLoadMarker', markerFilename [, isMultiMarker][, patt_width][, patt_center_x][, patt_center_y]);";
	synopsis[i++] = "[templateMatchingInColor, imageProcessingFullSized, imageProcessingIdeal, trackingWithPCA] = PsychCV('ARTrackerSettings' [, templateMatchingInColor][, imageProcessingFullSized][, imageProcessingIdeal][, trackingWithPCA]);";
	synopsis[i++] = "[detectedMarkers] = PsychCV('ARDetectMarkers'[, markerSubset][, threshold] [, infoType]);";
	synopsis[i++] = "[scale, minDist, maxDist] = PsychCV('ARRenderSettings' [, scale][, minDist][, maxDist]);";
	synopsis[i++] = "PsychCV('ARRenderImage');";
//	synopsis[i++] 

	synopsis[i++] = NULL;  //this tells PsychDisplayScreenSynopsis where to stop
	if (i > MAX_SYNOPSIS_STRINGS) {
		PrintfExit("%s: increase dimension of synopsis[] from %ld to at least %ld and recompile.",__FILE__,(long)MAX_SYNOPSIS_STRINGS,(long)i);
	}
}

// Show help/usage strings:
PsychError PSYCHCVDisplaySynopsis(void)
{
	int i;
	
	for (i = 0; synopsisSYNOPSIS[i] != NULL; i++)
		printf("%s\n",synopsisSYNOPSIS[i]);
	
	return(PsychError_none);
}

// Module exit function: Shutdown what's running, clean up after yourself...
PsychError PsychCVExit(void)
{
	// Shutdown only if we are online:
	if (psychCVInitialized) {
		// Detach our callback function for low-level debug output:
		// TODO...
		
		// Perform AR toolkit shutdown, if needed:
		PsychCVARExit();
		
		// Mark us dead:
		psychCVInitialized = FALSE;
	}
	
	return(PsychError_none);
}

// Startup function: Called at module load time / First time invocation:
void PsychCVInitialize(void)
{
	// PortAudio already initialized?
	if (!psychCVInitialized) {
		// Setup callback function for low-level debug output:
		// TODO...

		// We're alive :-)
		psychCVInitialized = TRUE;
	}
}

/* PsychCV('Verbosity') - Set level of verbosity.
 */
PsychError PSYCHCVVerbosity(void) 
{
 	static char useString[] = "oldlevel = PsychCV('Verbosity' [,level]);";
	static char synopsisString[] = 
		"Set level of verbosity for error/warning/status messages. 'level' optional, new level "
		"of verbosity. 'oldlevel' is the old level of verbosity. The following levels are "
		"supported: 0 = Shut up. 1 = Print errors, 2 = Print also warnings, 3 = Print also some info, "
		"4 = Print more useful info (default), >5 = Be very verbose (mostly for debugging the driver itself). ";		
	static char seeAlsoString[] = "";	 
	
	int level= -1;

	// Setup online help: 
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
	
	PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(1));	 // The maximum number of outputs

	PsychCopyInIntegerArg(1, kPsychArgOptional, &level);
	if (level < -1) PsychErrorExitMsg(PsychError_user, "Invalid level of verbosity provided. Valid are levels of zero and greater.");
	
	// Return current/old level:
	PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) verbosity);

	// Set new level, if one was provided:
	if (level > -1) verbosity = level;

	return(PsychError_none);
}

/* Runtime parameters for the OpenEyes computer vision based eyetracker: See source of cvEyeTracker.cc for
 * meaning of the parameters. They are just cached here in PsychCV.c for simplicity of internal implementation.
 *
 * They get initialized to their defaults (according to original OpenEyes implementation) in PsychCVOpenEyesInitialize().
 * After that they're cached, returned during queries from the cached local copy, and send to cvEyeTracker.cc when they're
 * changed in PSYCHCVOpenEyesParameters().
 */

static int pupilEdgeThreshold, starburstRays, minFeatureCandidates, corneaWindowSize, edgeThreshold, gaussWidth;
static double maxPupilEccentricity, initialAngleSpread;
static double fanoutAngle1, fanoutAngle2;
static int featuresPerRay, specialFlags;

#ifdef PSYCHCV_USE_OPENCV

/* PsychCV('OpenEyesInitialize') - Initialize a new tracking session with OpenEyes:
 */
PsychError PSYCHCVOpenEyesInitialize(void)
{

 	static char useString[] = "[EyeImageMemBuffer, EyeColorImageMemBuffer, SceneImageMemBuffer, ThresholdImageMemBuffer, EllipseImageMemBuffer] = PsychCV('OpenEyesInitialize', handle [, eyeChannels] [, eyeWidth][, eyeHeight][, sceneWidth][, sceneHeight][, logfilename]);";
	static char synopsisString[] = 
		"Set level of verbosity for error/warning/status messages. 'level' optional, new level "
		"of verbosity. 'oldlevel' is the old level of verbosity. The following levels are "
		"supported: 0 = Shut up. 1 = Print errors, 2 = Print also warnings, 3 = Print also some info, "
		"4 = Print more useful info (default), >5 = Be very verbose (mostly for debugging the driver itself). ";		
	static char seeAlsoString[] = "";	 
	
	int				handle = -1;
	int				eyewidth, eyeheight, eyechannels, scenewidth, sceneheight;
	char*			logfilename = NULL;
	void*			eyeImage = NULL;
	void*			eyeColorImage = NULL;
	void*			sceneImage = NULL;
	void*			thresholdImage = NULL;
	void*			ellipseImage = NULL;
	
	// Setup online help: 
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
	
	PsychErrorExit(PsychCapNumInputArgs(7));     // The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(5));	 // The maximum number of outputs

	// Get tracker handle: This is not used yet, just here for future extensions:
	PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
	if (handle < 0) PsychErrorExitMsg(PsychError_user, "Invalid trackerhandle 'handle' provided. Valid handles are of value zero or greater.");
	
	eyechannels = 1;
	PsychCopyInIntegerArg(2, kPsychArgOptional, &eyechannels);
	if (eyechannels < 1 || eyechannels > 4) PsychErrorExitMsg(PsychError_invalidRectArg, "Invalid eye image channels provided. Must be between 1 and 4!");

	// Get optional eye ROI: We default to 640x480 pixels ROI:
	eyewidth = 640;
	PsychCopyInIntegerArg(3, kPsychArgOptional, &eyewidth);
	if (eyewidth < 1) PsychErrorExitMsg(PsychError_invalidRectArg, "Invalid eye image width provided. Must be at least 1 pixel!");

	eyeheight = 480;
	PsychCopyInIntegerArg(4, kPsychArgOptional, &eyeheight);
	if (eyeheight < 1) PsychErrorExitMsg(PsychError_invalidRectArg, "Invalid eye image height provided. Must be at least 1 pixel!");

	scenewidth = 640;
	PsychCopyInIntegerArg(5, kPsychArgOptional, &scenewidth);
	if (scenewidth < 1) PsychErrorExitMsg(PsychError_invalidRectArg, "Invalid scene image width provided. Must be at least 1 pixel!");

	sceneheight = 480;
	PsychCopyInIntegerArg(6, kPsychArgOptional, &sceneheight);
	if (sceneheight < 1) PsychErrorExitMsg(PsychError_invalidRectArg, "Invalid scene image height provided. Must be at least 1 pixel!");
	
	// Get optional name of tracker logfile. If it is empty, we disable logging:
	PsychAllocInCharArg(7, kPsychArgOptional, &logfilename);
	
	// Call tracker init call:
	if (!cvEyeTrackerInitialize(logfilename, eyewidth, eyeheight, eyechannels, &eyeImage, &eyeColorImage, scenewidth, sceneheight, &sceneImage, &ellipseImage, &thresholdImage)) {
		// Failed!
		PsychErrorExitMsg(PsychError_system, "OpenEyes tracker setup failed.");
	}
	
	// Return double-encoded void* memory pointer to eye image input/output buffer:
	PsychCopyOutDoubleArg(1, kPsychArgRequired, PsychPtrToDouble(eyeImage));

	// Return double-encoded void* memory pointer to eye color image input/output buffer:
	PsychCopyOutDoubleArg(2, kPsychArgRequired, PsychPtrToDouble(eyeColorImage));

	// Return double-encoded void* memory pointer to eye image input/output buffer:
	PsychCopyOutDoubleArg(3, kPsychArgOptional, PsychPtrToDouble(sceneImage));

	// Return double-encoded void* memory pointer to eye image input/output buffer:
	PsychCopyOutDoubleArg(4, kPsychArgOptional, PsychPtrToDouble(thresholdImage));

	// Return double-encoded void* memory pointer to eye image input/output buffer:
	PsychCopyOutDoubleArg(5, kPsychArgOptional, PsychPtrToDouble(ellipseImage));

	// Setup default parameters for OpenEyes:
	pupilEdgeThreshold = 20;
	starburstRays = 18;
	minFeatureCandidates = 10;
	corneaWindowSize = 301;
	edgeThreshold = 20;
	gaussWidth = 5;
	maxPupilEccentricity = 1.2;	// Original is 2.0 (!)
	initialAngleSpread = 360.0;
	fanoutAngle1 = 0.0;
	fanoutAngle2 = 180.0;
	featuresPerRay = 1;
	specialFlags = 0x0;
	
	// Commit default parameters to tracker:
	cvEyeTrackerSetParameters(pupilEdgeThreshold, starburstRays, minFeatureCandidates, corneaWindowSize, edgeThreshold, gaussWidth, maxPupilEccentricity, initialAngleSpread * PI/180,
							  fanoutAngle1 * PI/180, fanoutAngle2 * PI/180, featuresPerRay, specialFlags);
	
	return(PsychError_none);
}
 
/* PsychCV('OpenEyesShutdown') - Shutdown a tracking session with OpenEyes:
 */
PsychError PSYCHCVOpenEyesShutdown(void)
{
 	static char useString[] = "PsychCV('OpenEyesShutdown', handle);";
	static char synopsisString[] = 
		"Set level of verbosity for error/warning/status messages. 'level' optional, new level "
		"of verbosity. 'oldlevel' is the old level of verbosity. The following levels are "
		"supported: 0 = Shut up. 1 = Print errors, 2 = Print also warnings, 3 = Print also some info, "
		"4 = Print more useful info (default), >5 = Be very verbose (mostly for debugging the driver itself). ";		
	static char seeAlsoString[] = "";	 
	
	int				handle = -1;
	
	// Setup online help: 
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
	
	PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(0));	 // The maximum number of outputs

	// Get tracker handle: This is not used yet, just here for future extensions:
	PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
	if (handle < 0) PsychErrorExitMsg(PsychError_user, "Invalid trackerhandle 'handle' provided. Valid handles are of value zero or greater.");
	
	// Call tracker shutdown call:
	if (!cvEyeTrackerShutdown()) {
		// Failed!
		PsychErrorExitMsg(PsychError_system, "OpenEyes tracker shutdown failed.");
	}

	return(PsychError_none);
}

/* PsychCV('OpenEyesTrackEyePosition') - Perform a single eye tracking cycle with OpenEyes and return result:
 */
PsychError PSYCHCVOpenEyesTrackEyePosition(void)
{
 	static char useString[] = "EyeResult = PsychCV('OpenEyesTrackEyePosition', handle [, mode][, px][, py][, c1][, c2]);";
	static char synopsisString[] = 
		"Set level of verbosity for error/warning/status messages. 'level' optional, new level "
		"of verbosity. 'oldlevel' is the old level of verbosity. The following levels are "
		"supported: 0 = Shut up. 1 = Print errors, 2 = Print also warnings, 3 = Print also some info, "
		"4 = Print more useful info (default), >5 = Be very verbose (mostly for debugging the driver itself). ";		
	static char seeAlsoString[] = "";	 

	PsychGenericScriptType 	*eyeStruct;
	const char *FieldNames[]={"GazeX", "GazeY", "PupilX", "PupilY", "CorneaX", "CorneaY", "Valid", "Count", "Timestamp"};
	const int FieldCount = 9;
	
	int					handle = -1;
	int					mode;
	double				px, py, minArea, maxArea;
	psych_bool				useGUI = FALSE;
	PsychCVEyeResult	eyeResult;

	// Setup online help: 
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
	
	PsychErrorExit(PsychCapNumInputArgs(6));     // The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(1));	 // The maximum number of outputs

	// Get tracker handle: This is not used yet, just here for future extensions:
	PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
	if (handle < 0) PsychErrorExitMsg(PsychError_user, "Invalid trackerhandle 'handle' provided. Valid handles are of value zero or greater.");
	
	mode = 0;
	PsychCopyInIntegerArg(2, kPsychArgOptional, &mode);
	if (mode < 0 || mode > 6) PsychErrorExitMsg(PsychError_user, "Invalid mode provided. Valid modes between 0 and 4.");

	if (mode == 0 || mode == 3) {
		// Call standard tracker processing cycle:
		if (!cvEyeTrackerExecuteTrackingCycle(&eyeResult, useGUI)) {
			// Failed!
			PsychErrorExitMsg(PsychError_system, "OpenEyes tracking cycle failed.");
		}
	}
	
	// Any calibration mode requested?
	if (mode > 0) {
		if (mode == 1) cvEyeTrackerRecalibrate(FALSE);					// Activate scene calibration.
		if (mode == 2) cvEyeTrackerRecalibrate(TRUE);					// Reset scene calibration to "None".
		
		if (mode > 2) {
			// Need (px, py) location either for eyepos init or for setting of calibration points:
			PsychCopyInDoubleArg(3, kPsychArgRequired, &px);
			PsychCopyInDoubleArg(4, kPsychArgRequired, &py);
		
			if (mode == 3) cvEyeTrackerSetPupilLocation((int) px, (int) py);		// Set Starburst start location (px, py) for approximate pupil center.
			if (mode == 4) cvEyeTrackerAddCalibrationPoint((int) px, (int) py);		// Add (px, py) as a new eye -> scene calibration point.

			if (mode == 5) {
				minArea = 0.0;
				PsychCopyInDoubleArg(5, kPsychArgOptional, &minArea);
				maxArea = 1000.0 * 1000.0;
				PsychCopyInDoubleArg(6, kPsychArgOptional, &maxArea);
				cvEyeTrackerSetRansacConstraints(px, py, minArea, maxArea);	// Add minDist=px, maxDist=py as distance constraint for features.
			}
			
			if (mode == 6) {
				// Set manual override reference point (cornea x/y):
				cvEyeTrackerSetOverrideReferencePoint((int) px, (int) py);
			}
		}
	}

	// Create result return struct:
	PsychAllocOutStructArray(1, kPsychArgOptional, -1, FieldCount, FieldNames, &eyeStruct);
	PsychSetStructArrayDoubleElement("GazeX", 0, eyeResult.gaze_x, eyeStruct);
	PsychSetStructArrayDoubleElement("GazeY", 0, eyeResult.gaze_y, eyeStruct);
	PsychSetStructArrayDoubleElement("PupilX", 0, eyeResult.pupil_x, eyeStruct);
	PsychSetStructArrayDoubleElement("PupilY", 0, eyeResult.pupil_y, eyeStruct);
	PsychSetStructArrayDoubleElement("CorneaX", 0, eyeResult.cornea_x, eyeStruct);
	PsychSetStructArrayDoubleElement("CorneaY", 0, eyeResult.cornea_y, eyeStruct);
	PsychSetStructArrayBooleanElement("Valid", 0, eyeResult.validresult, eyeStruct);
	PsychSetStructArrayDoubleElement("Count", 0, eyeResult.trackcount, eyeStruct);
	PsychSetStructArrayDoubleElement("Timestamp", 0, eyeResult.timestamp, eyeStruct);

	return(PsychError_none);
}

/* PsychCV('OpenEyesParameters') - Change parameters of a tracking session with OpenEyes:
 */
PsychError PSYCHCVOpenEyesParameters(void)
{

 	static char useString[] = "[oldSettings, ...] = PsychCV('OpenEyesParameters', handle [, pupilEdgeThreshold][, starburstRays][, minFeatureCandidates][, corneaWindowSize][, edgeThreshold][, gaussWidth][, maxPupilEccentricity] [, initialAngleSpread] [, fanoutAngle1] [, fanoutAngle2] [, featuresPerRay] [, specialFlags]);";
	static char synopsisString[] = 
		"Set level of verbosity for error/warning/status messages. 'level' optional, new level "
		"of verbosity. 'oldlevel' is the old level of verbosity. The following levels are "
		"supported: 0 = Shut up. 1 = Print errors, 2 = Print also warnings, 3 = Print also some info, "
		"4 = Print more useful info (default), >5 = Be very verbose (mostly for debugging the driver itself). ";		
	static char seeAlsoString[] = "";	 
	
	int	handle = -1;
	
	// Setup online help: 
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
	
	PsychErrorExit(PsychCapNumInputArgs(13));     // The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(12));	 // The maximum number of outputs

	// Get tracker handle: This is not used yet, just here for future extensions:
	PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
	if (handle < 0) PsychErrorExitMsg(PsychError_user, "Invalid trackerhandle 'handle' provided. Valid handles are of value zero or greater.");
	
	// Firt copy out optional return values: The old / current settings of the tracker parameters:
	PsychCopyOutDoubleArg(1, kPsychArgOptional, pupilEdgeThreshold);
	PsychCopyOutDoubleArg(2, kPsychArgOptional, starburstRays);
	PsychCopyOutDoubleArg(3, kPsychArgOptional, minFeatureCandidates);
	PsychCopyOutDoubleArg(4, kPsychArgOptional, corneaWindowSize);
	PsychCopyOutDoubleArg(5, kPsychArgOptional, edgeThreshold);
	PsychCopyOutDoubleArg(6, kPsychArgOptional, gaussWidth);
	PsychCopyOutDoubleArg(7, kPsychArgOptional, maxPupilEccentricity);
	PsychCopyOutDoubleArg(8, kPsychArgOptional, initialAngleSpread);
	PsychCopyOutDoubleArg(9, kPsychArgOptional, fanoutAngle1);
	PsychCopyOutDoubleArg(10, kPsychArgOptional, fanoutAngle2);
	PsychCopyOutDoubleArg(11, kPsychArgOptional, featuresPerRay);
	PsychCopyOutDoubleArg(12, kPsychArgOptional, specialFlags);
	
	// Get optional parameters. The defaults are the original settings of OpenEyes, set in the initialization call...
	PsychCopyInIntegerArg(2, kPsychArgOptional, &pupilEdgeThreshold);
	if (pupilEdgeThreshold < 1) PsychErrorExitMsg(PsychError_user, "Invalid pupilEdgeThreshold provided. Must be at least 1 unit!");

	PsychCopyInIntegerArg(3, kPsychArgOptional, &starburstRays);
	if (starburstRays < 2) PsychErrorExitMsg(PsychError_user, "Invalid number of starburstRays provided. Must be at least 2 rays!");

	PsychCopyInIntegerArg(4, kPsychArgOptional, &minFeatureCandidates);
	if (minFeatureCandidates < 1) PsychErrorExitMsg(PsychError_user, "Invalid number of minFeatureCandidates provided. Must be at least 1 candidate!");

	PsychCopyInIntegerArg(5, kPsychArgOptional, &corneaWindowSize);
	if (corneaWindowSize < 3) PsychErrorExitMsg(PsychError_user, "Invalid corneaWindowSize provided. Must be at least 3 pixels!");
	
	PsychCopyInIntegerArg(6, kPsychArgOptional, &edgeThreshold);
	if (edgeThreshold < 5) PsychErrorExitMsg(PsychError_user, "Invalid edgeThreshold provided. Must be at least 5 units!");
	
	PsychCopyInIntegerArg(7, kPsychArgOptional, &gaussWidth);
	if (gaussWidth < 0) PsychErrorExitMsg(PsychError_user, "Invalid gaussWidth provided. Must be at least 0 units!");
	if (gaussWidth % 2 == 0) gaussWidth += 1; // Make it odd if it's even!
	
	PsychCopyInIntegerArg(7, kPsychArgOptional, &gaussWidth);
	if (gaussWidth < 0) PsychErrorExitMsg(PsychError_user, "Invalid gaussWidth provided. Must be at least 0 units!");

	PsychCopyInDoubleArg(8, kPsychArgOptional, &maxPupilEccentricity);
	if (maxPupilEccentricity < 1.0) PsychErrorExitMsg(PsychError_user, "Invalid maxPupilEccentricity provided. Must be at least 1.0 units!");

	PsychCopyInDoubleArg(9, kPsychArgOptional, &initialAngleSpread);
	if (initialAngleSpread > 360 || initialAngleSpread < 0) PsychErrorExitMsg(PsychError_user, "Invalid initialAngleSpread provided. Must be between 0 and 360 deg.!");

	// Fanout midline angle in degrees for beam fan 1: This parameter is unlimited.
	PsychCopyInDoubleArg(10, kPsychArgOptional, &fanoutAngle1);

	// Fanout midline angle in degrees for beam fan 2: This parameter is unlimited.
	PsychCopyInDoubleArg(11, kPsychArgOptional, &fanoutAngle2);

	PsychCopyInIntegerArg(12, kPsychArgOptional, &featuresPerRay);
	if (featuresPerRay < 1) PsychErrorExitMsg(PsychError_user, "Invalid featuresPerRay provided. Must be at least 1!");

	PsychCopyInIntegerArg(13, kPsychArgOptional, &specialFlags);
	if (specialFlags < 0) PsychErrorExitMsg(PsychError_user, "Invalid specialFlags provided. Must be at least 0!");

	// Commit new parameters to tracker:
	cvEyeTrackerSetParameters(pupilEdgeThreshold, starburstRays, minFeatureCandidates, corneaWindowSize, edgeThreshold, gaussWidth, maxPupilEccentricity, initialAngleSpread * PI/180,
								fanoutAngle1 * PI/180, fanoutAngle2 * PI/180, featuresPerRay, specialFlags);

	return(PsychError_none);
}

#endif

PsychError PSYCHCVCopyMatrixToMemBuffer(void)
{
 	static char useString[] = "PsychCV('CopyMatrixToMemBuffer', matrix, memBufferPtr);";
	static char synopsisString[] = 
		"Copies a Matlab/Octave uint8 or double matrix 'matrix' into a C memory buffer, "
		"whose memory (void*) pointer is encoded as a double scalar value 'memBufferPtr'. "
		"The target buffer must be of sufficient size, no checking is performed! This is "
		"essentially a C memcpy() operation, so use with caution, or Matlab/Octave will crash!";
	static char seeAlsoString[] = "";	 
	
	int				m, n, p;
	unsigned char*	inMatrixUint8;
	void*			srcptr;
	double			doubleptr;
	void*			dstptr;
	int				elsize;
	
	// Setup online help: 
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
	
	PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(2)); // The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(0));	 // The maximum number of outputs

	// Get the uint8 or double matrix from Matlab/Octave:
	m = n = p = 0;
	if (PsychAllocInUnsignedByteMatArg(1, kPsychArgAnything, &m, &n, &p, (unsigned char**) &srcptr)) {
		elsize = sizeof(unsigned char);
	}
	else {
		if (PsychAllocInDoubleMatArg(1, kPsychArgAnything, &m, &n, &p, (double**) &srcptr)) {
			elsize = sizeof(double);
		}
		else {
			PsychErrorExitMsg(PsychError_user, "Invalid input matrix specified. Must be a uint8 or double matrix!");
		}
	}

	// Check dimensions:
	if (p < 1) p = 1;
	if (m < 1 || n < 1) PsychErrorExitMsg(PsychError_user, "Invalid input matrix specified. Must have non-zero row- and column count!");

	// Get the double-encoded void* to the destination buffer:
	PsychCopyInDoubleArg(2, kPsychArgRequired, &doubleptr);
	dstptr = PsychDoubleToPtr(doubleptr);
	
	// Perform memcpy to target buffer:
	memcpy(dstptr, (const void*) srcptr, elsize * m * n* p);
	
	return(PsychError_none);
}
