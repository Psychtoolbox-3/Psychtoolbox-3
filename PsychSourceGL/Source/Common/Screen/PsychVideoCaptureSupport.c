/*
	Common/Screen/PsychVideoCaptureSupport.c
	
	PLATFORMS:	
	
		This is the OS independent version (for now: Should work on OS-X and Windows)  
		A GNU/Linux specific version is stored in the /Linux/ folder. It has the
		same API - and therefore the same header file, but a pretty different
		implementation.

	AUTHORS:
	
		Mario Kleiner           mk              mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
        DESCRIPTION:
	
		Psychtoolbox functions for dealing with video capture devices.
	
	NOTES:

*/

#include "Screen.h"
#include <float.h>

// Forward declaration of internal helper function:
void PsychDeleteAllCaptureDevices(void);

// Record which defines capture engine independent state for a capture device:
typedef struct {
	int engineId;		// Type of capture engine: -1 == Free slot, 0 == Quicktime, 1 == LibDC, 2 == ARVideo, 3 == GStreamer
} PsychMasterVidcapRecordType;

static PsychMasterVidcapRecordType mastervidcapRecordBANK[PSYCH_MAX_CAPTUREDEVICES];
static int numCaptureRecords = 0;

static psych_bool firsttime = TRUE;

/*
 *     PsychVideoCaptureInit() -- Initialize video capture subsystem.
 *     This routine is called by Screen's RegisterProject.c PsychModuleInit()
 *     routine at Screen load-time. It clears out the vidcapRecordBANK to
 *     bring the subsystem into a clean initial state.
 */
void PsychVideoCaptureInit(void)
{
	// Initialize mastervidcapRecordBANK with NULL-entries:
	int i;
	for (i=0; i < PSYCH_MAX_CAPTUREDEVICES; i++) {
		mastervidcapRecordBANK[i].engineId = -1;
	}
	numCaptureRecords = 0;

	// Initialize the different capture engines:
	#ifdef PTBVIDEOCAPTURE_LIBDC
	PsychDCVideoCaptureInit();
	#endif

	#ifdef PTB_USE_GSTREAMER
	PsychGSVideoCaptureInit();
	#endif

	return;
}

void PsychEnumerateVideoSources(int engineId, int outPos)
{
	psych_bool dispatched = FALSE;

	#ifdef PTBVIDEOCAPTURE_LIBDC
	if (engineId == 1) {
		PsychErrorExitMsg(PsychError_user, "The DC1394 IIDC firewire videocapture engine does not support enumeration of video devices yet, sorry.");
	}
	#endif

	#ifdef PTB_USE_GSTREAMER
	if (engineId == 3) {
		// GStreamer device enumeration:
		PsychGSEnumerateVideoSources(outPos, -1);
		dispatched = TRUE;
	}
	#endif

	// Unsupported engine requested?
	if (!dispatched) PsychErrorExitMsg(PsychError_user, "The requested video capture engine is not supported on your system, either not at all, or has been disabled at compile time.");

	return;
}

/*
 *      PsychOpenVideoCaptureDevice() -- Create a video capture object.
 *
 *      This function tries to open a video capture device and return the associated captureHandle for it.
 *
 *		engineId = Type of video capture engine to use: 0 == Quicktime sequence grabbers, 1 == LibDC1394-V2 IIDC Firewire.
 *      win = Pointer to window record of associated onscreen window.
 *      deviceIndex = Index of the grabber device. (Currently ignored)
 *      capturehandle = handle to the new capture object.
 *      capturerectangle = If non-NULL a ptr to a PsychRectangle which contains the ROI for capture.
 *      The following arguments are currently ignored on Windows and OS-X:
 *      reqdepth = Number of layers for captured output textures. (0=Don't care, 1=LUMINANCE8, 2=LUMINANCE8_ALPHA8, 3=RGB8, 4=RGBA8)
 *      num_dmabuffers = Number of buffers in the ringbuffer queue (e.g., DMA buffers) - This is OS specific. Zero = Don't care.
 *      allow_lowperf_fallback = If set to 1 then PTB can use a slower, low-performance fallback path to get nasty devices working.
 *      targetmoviefilename = NULL == Only live capture, non-NULL == Pointer to char-string with name of target QT file for video recording.
 *      recordingflags = Only used for recording: Request audio recording, ram recording vs. disk recording and such...
 *	   	// Query optional movie recording flags:
 *		// 0 = Record video, stream to disk immediately (slower, but unlimited recording duration).
 *		// 1 = Record video, stream to memory, then at end of recording to disk (limited duration by RAM size, but faster).
 *		// 2 = Record audio as well.
 *
 */
psych_bool PsychOpenVideoCaptureDevice(int engineId, PsychWindowRecordType *win, int deviceIndex, int* capturehandle, double* capturerectangle,
				 int reqdepth, int num_dmabuffers, int allow_lowperf_fallback, char* targetmoviefilename, unsigned int recordingflags)
{
	int i, slotid;
	psych_bool dispatched = FALSE;
	*capturehandle = -1;

	// Sanity checking:
	if (!PsychIsOnscreenWindow(win)) {
		PsychErrorExitMsg(PsychError_user, "Provided windowPtr is not an onscreen window.");
	}

	if (numCaptureRecords >= PSYCH_MAX_CAPTUREDEVICES) {
		PsychErrorExitMsg(PsychError_user, "Allowed maximum number of simultaneously open capture devices exceeded!");
	}

	// Search first free slot in mastervidcapRecordBANK:
	for (i=0; (i < PSYCH_MAX_CAPTUREDEVICES) && (mastervidcapRecordBANK[i].engineId != -1); i++);
        if (i>=PSYCH_MAX_CAPTUREDEVICES) {
		PsychErrorExitMsg(PsychError_user, "Allowed maximum number of simultaneously open capture devices exceeded!");
	}

	// Slot slotid will contain the record for our new capture object:
	slotid=i;
	
	// Decide which engine to use and dispatch into proper open function:	
	#ifdef PTBVIDEOCAPTURE_LIBDC
	if (engineId == 1) {
		// LibDC1394 video capture:
		if (!PsychDCOpenVideoCaptureDevice(slotid, win, deviceIndex, capturehandle, capturerectangle, reqdepth,
						   num_dmabuffers, allow_lowperf_fallback, targetmoviefilename, recordingflags)) {
			// Probably won't ever reach this point due to error handling triggered in subfunction... anyway...
			return(FALSE);
		}
		dispatched = TRUE;
	}
	#endif
	
	#ifdef PTB_USE_GSTREAMER
	if (engineId == 3) {
		// GStreamer video capture:
		if (!PsychGSOpenVideoCaptureDevice(slotid, win, deviceIndex, capturehandle, capturerectangle, reqdepth,
						   num_dmabuffers, allow_lowperf_fallback, targetmoviefilename, recordingflags)) {
			// Probably won't ever reach this point due to error handling triggered in subfunction... anyway...
			return(FALSE);
		}
		dispatched = TRUE;
	}
	#endif

	// Unsupported engine requested?
	if (!dispatched) PsychErrorExitMsg(PsychError_user, "The requested video capture engine is not supported on your system, either not at all, or has been disabled at compile time.");
	
	// Ok, new capture device for requested capture engine created...
	
	// Assign new record:
	mastervidcapRecordBANK[slotid].engineId = engineId;    

	// Assign final handle:
	*capturehandle = slotid;

	// Increase counter:
	numCaptureRecords++;

	// Ready.
	return(TRUE);
}

/*
 *  PsychCloseVideoCaptureDevice() -- Close a capture device and release all associated ressources.
 */
void PsychCloseVideoCaptureDevice(int capturehandle)
{
	if (capturehandle < 0 || capturehandle >= PSYCH_MAX_CAPTUREDEVICES) {
		PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided!");
	}
	
	if (mastervidcapRecordBANK[capturehandle].engineId == -1) {
		PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided. No capture device associated with this handle !!!");
	}
        
	// Call engine specific method:
	#ifdef PTBVIDEOCAPTURE_LIBDC
	if (mastervidcapRecordBANK[capturehandle].engineId == 1) PsychDCCloseVideoCaptureDevice(capturehandle);
	#endif
	
	#ifdef PTB_USE_GSTREAMER
	if (mastervidcapRecordBANK[capturehandle].engineId == 3) PsychGSCloseVideoCaptureDevice(capturehandle);	
	#endif

	// Release record:
	mastervidcapRecordBANK[capturehandle].engineId = -1;
	
	// Decrease counter:
	if (numCaptureRecords>0) numCaptureRecords--;
        
	return;
}

/*
 *  PsychDeleteAllCaptureDevices() -- Delete all capture objects and release all associated ressources.
 */
void PsychDeleteAllCaptureDevices(void)
{
    int i;
    for (i=0; i<PSYCH_MAX_CAPTUREDEVICES; i++) {
        if (mastervidcapRecordBANK[i].engineId !=-1) PsychCloseVideoCaptureDevice(i);
    }
    return;
}


/*
 *  PsychGetTextureFromCapture() -- Create an OpenGL texturemap from a specific videoframe from given capture object.
 *
 *  win = Window pointer of onscreen window for which a OpenGL texture should be created.
 *  capturehandle = Handle to the capture object.
 *  checkForImage = >0 == Just check if new image available, 0 == really retrieve the image, blocking if necessary.
 *                   2 == Check for new image, block inside this function (if possible) if no image available.
 *  timeindex = This parameter is currently ignored and reserved for future use.
 *  out_texture = Pointer to the Psychtoolbox texture-record where the new texture should be stored.
 *  presentation_timestamp = A ptr to a double variable, where the presentation timestamp of the returned frame should be stored.
 *  summed_intensity = An optional ptr to a double variable. If non-NULL, then sum of intensities over all channels is calculated and returned.
 *  outrawbuffer = An optional ptr to a memory buffer of sufficient size. If non-NULL, the buffer will be filled with the captured raw image data, e.g., for use inside Matlab or whatever...
 *  Returns Number of pending or dropped frames after fetch on success (>=0), -1 if no new image available yet, -2 if no new image available and there won't be any in future.
 */
int PsychGetTextureFromCapture(PsychWindowRecordType *win, int capturehandle, int checkForImage, double timeindex, PsychWindowRecordType *out_texture, double *presentation_timestamp, double* summed_intensity, rawcapimgdata* outrawbuffer)
{        
	// Sanity checks:
	if (capturehandle < 0 || capturehandle >= PSYCH_MAX_CAPTUREDEVICES || mastervidcapRecordBANK[capturehandle].engineId == -1) {
		PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided.");
	}
    
	// Call engine specific method:
	#ifdef PTBVIDEOCAPTURE_LIBDC
	if (mastervidcapRecordBANK[capturehandle].engineId == 1) return(PsychDCGetTextureFromCapture(win, capturehandle, checkForImage, timeindex, out_texture, presentation_timestamp, summed_intensity, outrawbuffer));
	#endif
	
	#ifdef PTB_USE_GSTREAMER
	if (mastervidcapRecordBANK[capturehandle].engineId == 3) return(PsychGSGetTextureFromCapture(win, capturehandle, checkForImage, timeindex, out_texture, presentation_timestamp, summed_intensity, outrawbuffer));	
	#endif
	
    return(-2);
}

/*
 *  PsychVideoCaptureRate() - Start- and stop video capture.
 *
 *  capturehandle = Grabber to start-/stop.
 *  playbackrate = zero == Stop capture, non-zero == Capture
 *  dropframes = At 'start': Decide if low latency capture shall be used. At 'stop' If zero, don't
 *               discard pending buffers in internal capture queue.
 *  startattime = Deadline (in system time) to wait for before real start of capture.
 *  Returns Number of dropped frames during capture.
 */
int PsychVideoCaptureRate(int capturehandle, double capturerate, int dropframes, double* startattime)
{
	if (capturehandle < 0 || capturehandle >= PSYCH_MAX_CAPTUREDEVICES) {
		PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided!");
	}
        
	// Call engine specific method:
	#ifdef PTBVIDEOCAPTURE_LIBDC
	if (mastervidcapRecordBANK[capturehandle].engineId == 1) return(PsychDCVideoCaptureRate(capturehandle, capturerate, dropframes, startattime));
	#endif

	#ifdef PTB_USE_GSTREAMER
	if (mastervidcapRecordBANK[capturehandle].engineId == 3) return(PsychGSVideoCaptureRate(capturehandle, capturerate, dropframes, startattime));	
	#endif

	return(0);
}

/* Set capture device specific parameters:
 * On OS-X and Windows (and therefore in this implementation) this is currently a no-op, until
 * we find out how to do this with the Sequence-Grabber API.
 */
double PsychVideoCaptureSetParameter(int capturehandle, const char* pname, double value)
{
	// Valid handle provided?
	if (capturehandle < 0 || capturehandle >= PSYCH_MAX_CAPTUREDEVICES) {
		PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided!");
	}
	
	// Call engine specific method:
	#ifdef PTBVIDEOCAPTURE_LIBDC
	if (mastervidcapRecordBANK[capturehandle].engineId == 1) return(PsychDCVideoCaptureSetParameter(capturehandle, pname, value)); 
	#endif

	#ifdef PTB_USE_GSTREAMER
	if (mastervidcapRecordBANK[capturehandle].engineId == 3) return(PsychGSVideoCaptureSetParameter(capturehandle, pname, value));	
	#endif

	return(0);
}

/*
 *  void PsychExitVideoCapture() - Shutdown handler.
 *
 *  This routine is called by Screen('CloseAll') and on clear Screen time to
 *  do final cleanup. It deletes all capture objects
 *
 */
void PsychExitVideoCapture(void)
{
	// Release all capture devices:
	PsychDeleteAllCaptureDevices();
    
	// Call engine specific method:
	#ifdef PTBVIDEOCAPTURE_LIBDC
	PsychDCExitVideoCapture();
	#endif

	#ifdef PTB_USE_GSTREAMER
	PsychGSExitVideoCapture();	
	#endif

	firsttime = TRUE;
	return;
}
