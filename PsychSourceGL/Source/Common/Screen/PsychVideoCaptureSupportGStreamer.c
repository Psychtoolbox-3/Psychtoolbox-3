/*
	PsychSourceGL/Source/Common/Screen/PsychVideoCaptureSupportGStreamer.c
 
	PLATFORMS:	
	
	GNU/Linux, Apple MacOS/X and MS-Windows.
 
	AUTHORS:
	
	Mario Kleiner           mk              mario.kleiner@tuebingen.mpg.de
 
	HISTORY:
	
	9.01.2011               Created initial version.
	8.04.2011               Make video/audio recording work ok.
	5.06.2011               Make video/audio recording good enough
                            for initial release on Linux.

	DESCRIPTION:
	
	This is the videocapture engine based on the free software (LGPL'ed)
	GStreamer multimedia framework. It supports video capture, sound capture and
	recording of video and/or sound.
	 
	NOTES:

	The following works well so far (tested on Linux with Webcams):
 
	* Video device enumeration and selection (tested with iSight, USB webcams, DV Firewire cams, ...)

        * Simultaneous capture from 2 video cams.

        * Videocapture -> Live texture display. All formats L, LA, RGB, RGBA, YUYV.

        * Video-only recording, also optionally combined with live texture display.

	* Video + audio recording (also with optional live display).

        * Timestamping, GetSecs or pipeline running time timebase.

        * Low latency (dropframes) capture, and non-dropframes capture.

        * Selection of resolution and framerate.

        * Definition of ROI's and application to video feed and/or video recording.

        * Video codec selection.

        * Provides interface to change camera settings like exposure time, gain, contrast, etc.

	=> Most functionality for typical everyday tasks works perfect or reasonably well.
    => Some issues for special case apps persist, as written below.

	TODO:

	The following problems/limitations exist, which need to be fixed at some point:

        * Some codecs (e.g., huffyuv and h263) don't work yet. Some others show low quality or
          performance. Need to optimize parameters.

 */

#ifdef PTB_USE_GSTREAMER

#include "Screen.h"
#include <float.h>
#include <locale.h>
#include "PsychVideoCaptureSupport.h"

// These are the includes for GStreamer:
#include <glib.h>
#include <gst/gst.h>
#include <gst/app/gstappsink.h>
#include <gst/interfaces/propertyprobe.h>
#include <gst/interfaces/colorbalance.h>

// Compile-Time enable GStreamer encoding profile support by default on
// MacOSX and on 64-Bit Windows. These are the platforms where we know
// for sure that the installed GStreamer runtimes support this feature,
// because 64-Bit Windows == GStreamer-SDK and OSX == GStreamer-SDK or
// Homebrews GStreamer. Both Homebrew installs and the SDK do support this.
//
// 32-Bit Windows requires the OSSBuilds GStreamer which is too old to support this, but
// there we use camerabin1 + old setup path anyway, so no loss.
//
// Linux is mixed: Recent distros support it, old ones don't. I'm simply too lazy atm.
// to implement dynamic detection of support and dynamic linking on Linux, and use of
// camerabin2 + these profiles doesn't have any advantage over use of camerabin1 + the old
// codec setup, so as the old-school method is supported everywhere on Linux, we simply
// default to camerabin1 capture+recording and old style setup on Linux and safe ourselves
// a bit of hassle. Once we support a GStreamer V1.0 multi-media backend, we can implement
// this unconditionally without any pain or compatibility issues.
#if (PSYCH_SYSTEM == PSYCH_OSX) || ((PSYCH_SYSTEM == PSYCH_WINDOWS) && (defined(__LP64__) || defined(_M_IA64) || defined(_WIN64)))
#define PTB_USE_GSTENCODINGPROFILES 1
#endif

// For camerabin2 based video encoding:
#ifdef PTB_USE_GSTENCODINGPROFILES
#include <gst/pbutils/encoding-profile.h>
#endif

// Name of the presets for camerabin2 video recording with encoding profiles:
const char video_preset_name[] = "ptb3_videorecording_presets_video";
const char audio_preset_name[] = "ptb3_videorecording_presets_audio";
const char muxer_preset_name[] = "ptb3_videorecording_presets_muxer";

// 0 = No camerabin, 1 = camerabin, 2 = camerabin2:
static unsigned int usecamerabin = 1;

#define PSYCH_MAX_VIDSRC    256
PsychVideosourceRecordType *devices = NULL;
int ntotal = 0;

// Record which defines all state for a capture device:
typedef struct {
	int valid;                        // Is this a valid device record? zero == Invalid.
	psych_mutex mutex;
	psych_condition condition;
	int frameAvail;                   // Number of frames in videosink.
	int preRollAvail;
	GstElement *camera;               // Ptr to a GStreamer camera object that holds the internal state for such cams.
	GMainLoop *VideoContext;          // Message bus context for delivery of status/error/warning messages by GStreamer.
	GstElement *videosink;            // Our appsink for retrieval of live video feeds --> PTB texture conversion.
	GstElement *videosource;          // The videosource, encapsulating the actual video capture device.
	GstElement *videowrappersrc;      // The wrappercambinsrc used for encapsulating videosource for camerabin2 capture.
	GstElement *videoenc;             // Video encoder for video recording.
	GstElement *videorate_filter;     // Video framerate converter to guarantee constant framerate and a/v sync for recording.
	GstClockTime lastSavedBaseTime;   // Last time the pipeline was put into PLAYING state. Saved at StopCapture time before stop.
	int nrAudioTracks;
	int nrVideoTracks;
	int dropframes;                   // 1 == Always deliver most recent frame in FIFO, even if dropping of frames is neccessary.
	unsigned char* scratchbuffer;     // Scratch buffer for YUV->RGB conversion.
	int reqpixeldepth;                // Requested depth of single pixel in output texture.
	int pixeldepth;                   // Depth of single pixel from grabber in bits.
	int num_dmabuffers;               // Number of DMA ringbuffers to use in DMA capture.
	int nrframes;                     // Total count of decompressed images.
	double fps;                       // Acquisition framerate of capture device.
	int width;                        // Width x height of captured images from video source (prior to cropping etc.)
	int height;                       // This is effectively the configured resolution (width x height pixels) of the video source.
	int frame_width;                  // Width x height of frames returned from videosink for OpenGL texture conversion or
	int frame_height;                 // image return in raw image buffer.
	double last_pts;                  // Capture timestamp of previous frame.
	double current_pts;               // Capture timestamp of current frame.
	int current_dropped;              // Dropped count for this fetch cycle...
	int nr_droppedframes;             // Counter for dropped frames.
	int grabber_active;               // Grabber running?
	int recording_active;             // Movie file recording requested?
	unsigned int recordingflags;      // recordingflags, as passed to 'OpenCaptureDevice'.
	PsychRectType roirect;            // Region of interest rectangle - denotes subarea of full video capture area.
	double avg_decompresstime;        // Average time spent in Quicktime/Sequence Grabber decompressor.
	double avg_gfxtime;               // Average time spent in GWorld --> OpenGL texture conversion and statistics.
	int nrgfxframes;                  // Count of fetched textures.
	char* targetmoviefilename;        // Filename of a movie file to record.
	char* cameraFriendlyName;         // Camera friendly device name.
	char videosourcename[100];        // Plugin name of the videosource plugin.
} PsychVidcapRecordType;

static PsychVidcapRecordType vidcapRecordBANK[PSYCH_MAX_CAPTUREDEVICES];
static int numCaptureRecords = 0;
static psych_bool gs_firsttime = TRUE;
double gs_startupTime = 0.0;

// Forward declaration of internal helper function:
void PsychGSDeleteAllCaptureDevices(void);
int PsychGSDrainBufferQueue(PsychVidcapRecordType* capdev, int numFramesToDrain, unsigned int flags);


/*    PsychGetGSVidcapRecord() -- Given a handle, return ptr to video capture record.
 *    --> Internal helper function of PsychVideoCaptureSupport.
 */
PsychVidcapRecordType* PsychGetGSVidcapRecord(int deviceIndex)
{
	// Sanity checks:
	if (deviceIndex < 0) {
		PsychErrorExitMsg(PsychError_user, "Invalid (negative) deviceIndex for video capture device passed!");
	}
	
	if (numCaptureRecords >= PSYCH_MAX_CAPTUREDEVICES) {
		PsychErrorExitMsg(PsychError_user, "Invalid deviceIndex for video capture device passed. Index exceeds number of registered devices!");
	}
	
	if (!vidcapRecordBANK[deviceIndex].valid) {
		PsychErrorExitMsg(PsychError_user, "Invalid deviceIndex for video capture device passed. No such device open!");
	}
	
	// Ok, we have a valid device record, return a ptr to it:
	return(&vidcapRecordBANK[deviceIndex]);
}

#if PSYCH_SYSTEM == PSYCH_OSX
extern gboolean	gst_init_check(int *argc, char **argv[], GError **err) __attribute__((weak_import));
#endif

/* Internal: Check if GStreamer is already initialized. Initialize it,
 * if neccessary.
 */
void PsychGSCheckInit(const char* engineName)
{
	GError      *error = NULL;
	GstClock    *system_clock = NULL;

	if (gs_firsttime) {
		// First time invocation:
        
        // Check if GStreamer is properly installed and (can be) dynamically loaded and linked:
        #if PSYCH_SYSTEM != PSYCH_LINUX
            // On Windows and OSX, we need to delay-load the GStreamer libraries. This loading
            // and linking will automatically happen downstream. However, if delay loading
            // would fail, we would end up with a crash! For that reason, on MS-Windows, we
            // try to load the DLL, just to probe if the real load/link/bind op later on will
            // likely succeed. If the following LoadLibrary() call fails and returns NULL,
            // then we know we would end up crashing. We check for two versions of the dll's, as
            // different GStreamer runtime distributions use different filenames for the dll's.
            //
            // On OSX we check if the gst_init_check() function is defined, aka non-NULL. The
            // OSX linker sets the symbol to NULL if dynamic weak linking during runtime failed.
            // On failure we'll output some helpful error-message instead:
            #if PSYCH_SYSTEM == PSYCH_WINDOWS
                if (((NULL == LoadLibrary("libgstreamer-0.10.dll")) || (NULL == LoadLibrary("libgstapp-0.10.dll"))) &&
                    ((NULL == LoadLibrary("libgstreamer-0.10-0.dll")) || (NULL == LoadLibrary("libgstapp-0.10-0.dll")))) {
            #endif
            #if PSYCH_SYSTEM == PSYCH_OSX
                if (NULL == gst_init_check) {
            #endif
                // Failed:
                printf("\n\n\n");
                printf("PTB-ERROR: Tried to startup GStreamer multi-media framework. This didn't work, because one\n");
                printf("PTB-ERROR: of the required GStreamer runtime libraries failed to load, probably because it\n");
                printf("PTB-ERROR: could not be found, could not be accessed (e.g., due to permission problems),\n");
                printf("PTB-ERROR: or most likely because GStreamer isn't installed on this machine at all.\n");
                printf("PTB-ERROR: Another reason could be that you have GStreamer version 1.0 instead of the required\n");
                printf("PTB-ERROR: version 0.10 installed. The version 1 series is not yet supported.\n\n");
                #if PSYCH_SYSTEM == PSYCH_WINDOWS
                    printf("PTB-ERROR: The system returned error code %d.\n", GetLastError());
                #endif
                printf("PTB-ERROR: Please read the help by typing 'help GStreamer' for installation and troubleshooting\n");
                printf("PTB-ERROR: instructions.\n\n");
                printf("PTB-ERROR: Due to failed GStreamer initialization, the %s engine is disabled for this session.\n\n", engineName);

                // Quicktime supported on this setup?
                #if (PSYCH_SYSTEM != PSYCH_LINUX) && defined(PSYCHQTAVAIL)
                    // Yes. Give user a hint about this alternative, at least for movie playback or video capture,
                    // but not for movie writing:
                    if (NULL == strstr(engineName, "movie writing")) {
                        printf("PTB-TIP: As a stop-gap measure until you've installed or fixed GStreamer on your system,\n");
                        printf("PTB-TIP: you could try to use the legacy Quicktime based %s engine instead via use of the\n", engineName);
                        printf("PTB-TIP: override Screen('Preference', ...); switches 'DefaultVideoCaptureEngine' and\n");
                        printf("PTB-TIP: 'OverrideMultimediaEngine'.\n\n");
                    }
                #endif
                PsychErrorExitMsg(PsychError_user, "GStreamer initialization failed due to library loading problems. Aborted.");
            }
        #endif

		// Initialize GStreamer:
		if(!gst_init_check(NULL, NULL, &error)) {
			if (error) {
				printf("PTB-ERROR: GStreamer initialization failed with error: %s\n", (char*) error->message);
				g_error_free(error);
			}
			else {
				printf("PTB-ERROR: GStreamer initialization failed for unknown reason\n");
			}

			printf("PTB-ERROR: Due to failed GStreamer initialization, the %s engine is out of order.\n", engineName);
			PsychErrorExitMsg(PsychError_system, "GStreamer initialization failed! Aborted.");
		}

        // gs_startupTime is added to all timestamps from GStreamer to compensate for 
        // clock zero offset wrt. to our GetSecs() time:
        if (PSYCH_SYSTEM == PSYCH_WINDOWS) {
            // Windows: Zero-Point is time of GStreamer startup aka sometime
            // during execution of gst_init_check() above. Current system time
            // is our best approximation so far:
            PsychGetAdjustedPrecisionTimerSeconds(&gs_startupTime);
        }
        else {
            // Other OS: Zero-Point of GStreamer clock is identical to Zero-Point
            // of our GetSecs() clock, so apply zero-correction:
            // TODO FIXME: Check if this assumption holds on OSX for the Video capture engine.
            // We can only test this once we have 
            gs_startupTime = 0.0;
        }

        // Select opmode of GStreamers master clock:
        // We use monotonic clock on Windows and OS/X, as these correspond to the
        // clocks we use for GetSecs(), but realtime clock on Linux:
        system_clock = gst_system_clock_obtain();
        if (system_clock) {
            g_object_set(G_OBJECT(system_clock), "clock-type", ((PSYCH_SYSTEM == PSYCH_LINUX) ? GST_CLOCK_TYPE_REALTIME : GST_CLOCK_TYPE_MONOTONIC), NULL);
        }
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Using GStreamer version '%s'.\n", (char*) gst_version_string());

        // Reset firsttime flag:
        gs_firsttime = FALSE;
	}	
}

/*
 *     PsychGSVideoCaptureInit() -- Initialize video capture subsystem.
 *     This routine is called by Screen's RegisterProject.c PsychModuleInit()
 *     routine at Screen load-time. It clears out the vidcapRecordBANK to
 *     bring the subsystem into a clean initial state.
 */
void PsychGSVideoCaptureInit(void)
{
	// Initialize vidcapRecordBANK with NULL-entries:
	int i;
	for (i=0; i < PSYCH_MAX_CAPTUREDEVICES; i++) {
		vidcapRecordBANK[i].valid = 0;
	}    
	numCaptureRecords = 0;
	
	return;
}

/*
 *  void PsychGSExitVideoCapture() - Shutdown handler.
 *
 *  This routine is called by Screen('CloseAll') and on clear Screen time to
 *  do final cleanup. It deletes all capture objects
 *
 */
void PsychGSExitVideoCapture(void)
{
	// Release all capture devices
	PsychGSDeleteAllCaptureDevices();

	// Shutdown GStreamer framework, release all resources.
	// This is usually not needed/used. For memory leak
	// checking only!
	if (FALSE) {
		gs_firsttime = TRUE;
		gst_deinit();
	}

	return;
}

/*
 *  PsychGSDeleteAllCaptureDevices() -- Delete all capture objects and release all associated ressources.
 */
void PsychGSDeleteAllCaptureDevices(void)
{
	int i;
	for (i=0; i<PSYCH_MAX_CAPTUREDEVICES; i++) {
		if (vidcapRecordBANK[i].valid) PsychGSCloseVideoCaptureDevice(i);
	}
	return;
}


/* Perform one context loop iteration (for bus message handling) if doWait == false,
 * or two seconds worth of iterations if doWait == true. This drives the message-bus
 * callback, so needs to be performed to get any error reporting etc.
 */
int PsychGSProcessVideoContext(GMainLoop *loop, psych_bool doWait)
{
	double tdeadline, tnow;
	PsychGetAdjustedPrecisionTimerSeconds(&tdeadline);
	tnow = tdeadline;
	tdeadline+=2.0;

	if (NULL == loop) return(0);

	while (doWait && (tnow < tdeadline)) {
		// Perform non-blocking work iteration:
		if (!g_main_context_iteration(g_main_loop_get_context(loop), false)) PsychYieldIntervalSeconds(0.010);

		// Update time:
		PsychGetAdjustedPrecisionTimerSeconds(&tnow);
	}

	// Perform one more work iteration of the event context, but don't block:
	return(g_main_context_iteration(g_main_loop_get_context(loop), false));
}

/* Initiate pipeline state changes: Startup, Preroll, Playback, Pause, Standby, Shutdown. */
static psych_bool PsychVideoPipelineSetState(GstElement* camera, GstState state, double timeoutSecs)
{
    GstState			state_pending;
    GstStateChangeReturn	rcstate;

    gst_element_set_state(camera, state);

    // Non-Blocking, async?
    if (timeoutSecs < 0) return(TRUE);
 
    // Wait for up to timeoutSecs for state change to complete or fail:
    rcstate = gst_element_get_state(camera, &state, &state_pending, (GstClockTime) (timeoutSecs * 1e9));
    switch(rcstate) {
	case GST_STATE_CHANGE_SUCCESS:
		//printf("PTB-DEBUG: Statechange completed with GST_STATE_CHANGE_SUCCESS.\n");
	break;

	case GST_STATE_CHANGE_ASYNC:
		printf("PTB-INFO: Statechange in progress with GST_STATE_CHANGE_ASYNC.\n");
	break;

	case GST_STATE_CHANGE_NO_PREROLL:
		//printf("PTB-INFO: Statechange completed with GST_STATE_CHANGE_NO_PREROLL.\n");
	break;

	case GST_STATE_CHANGE_FAILURE:
		printf("PTB-ERROR: Statechange failed with GST_STATE_CHANGE_FAILURE!\n");
		return(FALSE);
	break;

	default:
		printf("PTB-ERROR: Unknown state-change result in preroll.\n");
		return(FALSE);
    }

    return(TRUE);
}

/* Receive messages from the playback pipeline message bus and handle them: */
static gboolean PsychVideoBusCallback(GstBus *bus, GstMessage *msg, gpointer dataptr)
{
  PsychVidcapRecordType* capdev = (PsychVidcapRecordType*) dataptr;

  switch (GST_MESSAGE_TYPE (msg)) {
    case GST_MESSAGE_EOS:
	if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Videobus: Message EOS received.\n"); fflush(NULL);
    break;

    case GST_MESSAGE_WARNING: {
      gchar  *debug;
      GError *error;

      gst_message_parse_warning(msg, &error, &debug);
      if (PsychPrefStateGet_Verbosity() > 3) { 
	      printf("PTB-WARNING: GStreamer videocapture engine reports this warning:\n"
		     "             Warning from element %s: %s\n", GST_OBJECT_NAME(msg->src), error->message);
	      printf("             Additional debug info: %s.\n", (debug) ? debug : "None");
      }

      g_free(debug);
      g_error_free(error);
      break;
    }

    case GST_MESSAGE_ERROR: {
      gchar  *debug;
      GError *error;

      gst_message_parse_error(msg, &error, &debug);
      if (PsychPrefStateGet_Verbosity() > 0) { 
	      printf("PTB-ERROR: GStreamer videocapture engine reports this error:\n"
		     "           Error from element %s: %s\n", GST_OBJECT_NAME(msg->src), error->message);
	      printf("           Additional debug info: %s.\n\n", (debug) ? debug : "None");

	      // Special tips for the challenged:
	      if (strstr(error->message, "property") || (debug && strstr(debug, "property"))) {
		      // Bailed due to unsupported x264enc parameter "speed-preset". Can be solved by upgrading
		      // GStreamer or the OS or the VideoCodec= override:
		      printf("PTB-TIP: The reason this failed is because your GStreamer codec installation is too outdated.\n");
		      printf("PTB-TIP: Either upgrade your GStreamer (plugin) installation to a more recent version,\n");
		      printf("PTB-TIP: or upgrade your operating system (e.g., Ubuntu 10.10 'Maverick Meercat' and later are fine).\n");
		      printf("PTB-TIP: A recent GStreamer installation is required to use all features and get optimal performance.\n");
		      printf("PTB-TIP: As a workaround, you can manually specify all codec settings, leaving out the unsupported\n");
		      printf("PTB-TIP: option. See 'help VideoRecording' on how to do that.\n\n");
	      }

	      if ((error->domain == GST_RESOURCE_ERROR) && (error->code != GST_RESOURCE_ERROR_NOT_FOUND)) {
		      printf("           This means that there was some problem with opening the video device (permissions etc.).\n\n");
	      }

	      if ((error->domain == GST_RESOURCE_ERROR) && (error->code == GST_RESOURCE_ERROR_NOT_FOUND)) {
		      printf("           This means that no such video device with the given name could be found.\n\n");
	      }
      }

      g_free(debug);
      g_error_free(error);
      break;
    }

    default:
      break;
  }

  return TRUE;
}


/* Called at each end-of-stream event at end of playback: */
static void PsychEOSCallback(GstAppSink *sink, gpointer user_data)
{
	PsychVidcapRecordType* capdev = (PsychVidcapRecordType*) user_data;

	PsychLockMutex(&capdev->mutex);
	printf("PTB-DEBUG: Videosink reached EOS.\n");
    // Signal EOS to trigger an abort of device open sequence:
    capdev->grabber_active = true;
	PsychUnlockMutex(&capdev->mutex);

	return;
}

static void PsychProbeBufferProps(GstBuffer *videoBuffer, int *w, int *h, double *fps)
{
	GstCaps                *caps;
	GstStructure	       *str;
	gint		       rate1, rate2;
	rate1 = rate2 = 0;
	str = NULL;
    
	if (videoBuffer) {
		caps = gst_buffer_get_caps(videoBuffer);
		if (caps) str = gst_caps_get_structure(caps, 0);
		if (fps) {
			if (str) gst_structure_get_fraction(str, "framerate", &rate1, &rate2);
			if (rate2 > 0) {
				*fps = ((double) rate1) / ((double) rate2);
			} else {
				*fps = 0.0;
			}
		}
		if (w && str) gst_structure_get_int(str,"width", w);
		if (h && str) gst_structure_get_int(str,"height", h);
		if (caps) gst_caps_unref(caps);
	}

	return;
}

/* Called whenever pipeline goes into pause mode.
 * Signals/handles arrival of preroll buffers. Used to detect/signal when
 * new videobuffers are available in non-playback mode, when device becomes ready.
 */
static GstFlowReturn PsychNewPrerollCallback(GstAppSink *sink, gpointer user_data)
{
	GstBuffer              *videoBuffer;

	PsychVidcapRecordType* capdev = (PsychVidcapRecordType*) user_data;

	PsychLockMutex(&capdev->mutex);
	videoBuffer = gst_app_sink_pull_preroll(GST_APP_SINK(capdev->videosink));
	if (videoBuffer) {
		PsychProbeBufferProps(videoBuffer, NULL, NULL, &capdev->fps);
	}

	if (PsychPrefStateGet_Verbosity() > 5) {
		printf("PTB-DEBUG: New PrerollBuffer received. fps = %f\n", capdev->fps);
		fflush(NULL);
	}

	capdev->preRollAvail++;

	PsychUnlockMutex(&capdev->mutex);

	return(GST_FLOW_OK);
}

/* Called whenever pipeline is in active playback and a new video frame arrives.
 * Used to detect/signal when new videobuffers are available in playback mode.
 */
static GstFlowReturn PsychNewBufferCallback(GstAppSink *sink, gpointer user_data)
{
	PsychVidcapRecordType* capdev = (PsychVidcapRecordType*) user_data;

	PsychLockMutex(&capdev->mutex);
	capdev->frameAvail++;
	if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: New Buffer received. %i\n", capdev->frameAvail);
	PsychSignalCondition(&capdev->condition);
	PsychUnlockMutex(&capdev->mutex);

	return(GST_FLOW_OK);
}

/* Not used by us, but needs to be defined as no-op anyway: */
static GstFlowReturn PsychNewBufferListCallback(GstAppSink *sink, gpointer user_data)
{
	PsychVidcapRecordType* capdev = (PsychVidcapRecordType*) user_data;

	PsychLockMutex(&capdev->mutex);
	//printf("PTB-DEBUG: New Bufferlist received.\n");
	PsychUnlockMutex(&capdev->mutex);

	return(GST_FLOW_OK);
}

/* Not used by us, but needs to be defined as no-op anyway: */
static void PsychDestroyNotifyCallback(gpointer user_data)
{
	return;
}

static GstAppSinkCallbacks videosinkCallbacks = {
    PsychEOSCallback,
    PsychNewPrerollCallback,
    PsychNewBufferCallback,
    PsychNewBufferListCallback,
    {NULL}
};

/*
 *  PsychGSCloseVideoCaptureDevice() -- Close a capture device and release all associated ressources.
 */
void PsychGSCloseVideoCaptureDevice(int capturehandle)
{
	// Retrieve device record for handle:
	PsychVidcapRecordType* capdev = PsychGetGSVidcapRecord(capturehandle);

	// Make sure GStreamer is ready:
	PsychGSCheckInit("videocapture");

	if (capdev->camera) {
		// Stop capture immediately if it is still running:
		PsychGSVideoCaptureRate(capturehandle, 0, 0, NULL);
	
		// Close & Shutdown camera, release ressources:
		// Stop video capture immediately:
		PsychVideoPipelineSetState(capdev->camera, GST_STATE_NULL, 20.0);

		// Delete camera for this handle:
		gst_object_unref(GST_OBJECT(capdev->camera));
		capdev->camera=NULL;
	}

	// Delete video context:
	if (capdev->VideoContext) g_main_loop_unref(capdev->VideoContext);
	capdev->VideoContext = NULL;
	
	PsychDestroyMutex(&capdev->mutex);
	PsychDestroyCondition(&capdev->condition);
	
	// Was our videosink (aka appsink) detached during device operation,
	// because videorecording was active but live feedback disabled, i.e.,
	// a fakesink was attached to camerabin? If so then camerabin didn't
	// own our videosink at destruction time and didn't auto-delete it.
	// We need to manually delete our orphaned videosink:
	if ((capdev->recording_active) && (capdev->recordingflags & 4) && (capdev->videosink)) {
		gst_object_unref(GST_OBJECT(capdev->videosink));
	}
	capdev->videosink = NULL;

	// This has been auto-destructed (hopefully) by camerabin:
	capdev->videosource = NULL;
    
	if (capdev->targetmoviefilename) free(capdev->targetmoviefilename);
	capdev->targetmoviefilename = NULL;

	if (capdev->cameraFriendlyName) free(capdev->cameraFriendlyName);
	capdev->cameraFriendlyName = NULL;

	// Invalidate device record to free up this slot in the array:
	capdev->valid = 0;
    
	// Decrease counter of open capture devices:
	if (numCaptureRecords>0) numCaptureRecords--;
	
	// Done.
	return;
}

/* Helper routine for PsychGSEnumerateVideoSources()
 *
 * Probes video source plugin 'srcname' [with class index 'classIndex' and
 * name 'className', using the property 'devHandlePropName' and probe strategy
 * as given by 'flags'. Adds all detected video input devices for that plugin
 * to the global 'devices' array and increases ntotal count accordingly.
 */
void PsychGSEnumerateVideoSourceType(const char* srcname, int classIndex, const char* className, const char* devHandlePropName, unsigned int flags)
{
	int					i, n, nmaxp, dopoke;
	char				port_str[64];
	char				class_str[64];
	int					inputIndex;
	GstElement			*videosource = NULL;
	GstPropertyProbe	*probe = NULL;
	GValueArray			*viddevs = NULL;
	GValue				*dev = NULL;
	char				*device_name = NULL;
	gchar               *pstring = NULL;
	psych_uint64		deviceURI = 0;
	
	// Info to the user:
	if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to probe %s as video source...\n", srcname);
	
	// Create video source plugin, define class name:
	videosource = gst_element_factory_make(srcname, "ptb_videosource");
	sprintf(class_str, "%s", className);

	// Nothing to do if no such video plugin available:
	if (!videosource) return;

	// No property probe interface for dc1394src or qtkitvideosrc, but "enumeration by trying" requested?
	// This is what we need to do for IIDC IEEE-1394 video sources via the dc1394src,
	// as it doesn't support property probe interface:
	if (!strcmp(srcname, "dc1394src") || !strcmp(srcname, "qtkitvideosrc")) {
		// Try a reasonable range of cameras, e.g., up to 100 cameras:
		n = 0;
        nmaxp  = (!strcmp(srcname, "dc1394src")) ? 100 : 5;
        dopoke = (!strcmp(srcname, "dc1394src")) ? 1   : 0;
		for (i = 0; i < nmaxp; i++) {
            // Only really probe if dopoke:
            if (dopoke) {
                // Set suspected camera id (select i'th camera on bus):
                g_object_set(G_OBJECT(videosource), devHandlePropName, i, NULL);
                // Try to set it to "paused" state, which should fail if no such
                // camera is connected:
                if (gst_element_set_state(videosource, GST_STATE_PAUSED) == GST_STATE_CHANGE_FAILURE) {
                    // No such camera connected. Game over, no need to probe further non-existent cams:
                    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: No camera %i connected to %s video input bus. Probe finished.\n", i, srcname);
                    break;
                }
            }

			// i'th camera exists. Probe and assign:
			inputIndex = i;
			devices[ntotal].deviceIndex = classIndex * 10000 + inputIndex;
			devices[ntotal].classIndex = classIndex;
			devices[ntotal].inputIndex = inputIndex;
			sprintf(devices[ntotal].deviceClassName, "%s", className);
			sprintf(devices[ntotal].deviceHandle, "%i", inputIndex);
			sprintf(devices[ntotal].deviceSelectorProperty, "%s", devHandlePropName);
			sprintf(devices[ntotal].deviceVideoPlugin, "%s", srcname);
			sprintf(devices[ntotal].deviceName, "%i", inputIndex);
			sprintf(devices[ntotal].devicePath, "%i", inputIndex);
			sprintf(devices[ntotal].device, "%i", inputIndex);
			devices[ntotal].deviceURI = inputIndex;

			ntotal++;
			n++;

			// Reset this cam:
			if (dopoke) gst_element_set_state(videosource, GST_STATE_READY);

			if (dopoke && (PsychPrefStateGet_Verbosity() > 4)) printf("PTB-INFO: %i'th %s camera enumerated.\n", inputIndex, srcname);
		}

		// Release videosource:
		gst_element_set_state(videosource, GST_STATE_NULL);
		gst_object_unref(GST_OBJECT(videosource));

		// Any success?
		if (n == 0) {
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: No video devices to enumerate for plugin '%s'.\n", class_str);
		}

		return;
	}

	// Generate property probe for videosource:
	// Make sure videosource implements property probe interface to avoid useless warning clutter,
	// unless on Linux, where this query actually causes an assertion on some GStreamer versions...
	if ((PSYCH_SYSTEM == PSYCH_LINUX) || gst_element_implements_interface(videosource, GST_TYPE_PROPERTY_PROBE)) {
		probe = GST_PROPERTY_PROBE(videosource);
	}
	else {
		probe = NULL;
	}
	
	if (probe) {
		if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Probing '%s' property...\n", devHandlePropName);
		viddevs = gst_property_probe_probe_and_get_values_name(probe, devHandlePropName);

		if (viddevs) {
			// Assign count of available devices:
			n = viddevs->n_values;
			
			// Iterate all available devices:
			for(i = 0; i < n; i++) {
				// Get device handle:
				dev = g_value_array_get_nth(viddevs, i);
				if (G_VALUE_HOLDS_STRING(dev)) {
					sprintf(port_str, "%s", (const char*) g_value_get_string(dev));

					// Select it for further probing:
					g_object_set(G_OBJECT(videosource), devHandlePropName, port_str, NULL);
				}

				if (G_VALUE_HOLDS_UINT64(dev)) {
					deviceURI = g_value_get_uint64(dev);

					// Select it for further probing:
					g_object_set(G_OBJECT(videosource), devHandlePropName, deviceURI, NULL);
					sprintf(port_str, "%llu", deviceURI);
				}

				inputIndex = i;
				devices[ntotal].deviceIndex = classIndex * 10000 + inputIndex;
				devices[ntotal].classIndex = classIndex;
				devices[ntotal].inputIndex = inputIndex;
				sprintf(devices[ntotal].deviceClassName, "%s", className);
				sprintf(devices[ntotal].deviceHandle, "%s", port_str);
				sprintf(devices[ntotal].deviceSelectorProperty, "%s", devHandlePropName);
				sprintf(devices[ntotal].deviceVideoPlugin, "%s", srcname);

				// Query and assign device specific parameters:
				pstring = NULL;
				if (g_object_class_find_property(G_OBJECT_GET_CLASS(videosource), "device")) {
					g_object_get(G_OBJECT(videosource), "device", &pstring, NULL);
				}
				
				if (pstring) {
					sprintf(devices[ntotal].device, "%s", pstring);
					g_free(pstring);
				}
				else {
					sprintf(devices[ntotal].device, "%s", port_str);
				}
				
				pstring = NULL;
				if (g_object_class_find_property(G_OBJECT_GET_CLASS(videosource), "device-path")) {
					g_object_get(G_OBJECT(videosource), "device-path", &pstring, NULL);
				}
				
				if (pstring) {
					sprintf(devices[ntotal].devicePath, "%s", pstring);
					g_free(pstring);
				}
				else {
					sprintf(devices[ntotal].devicePath, "%s", port_str);
				}
				
				pstring = NULL;
				if (g_object_class_find_property(G_OBJECT_GET_CLASS(videosource), "device-name")) {
					g_object_get(G_OBJECT(videosource), "device-name", &pstring, NULL);
				}
				
				if (pstring) {
					sprintf(devices[ntotal].deviceName, "%s", pstring);
					g_free(pstring);
				}
				else {
					sprintf(devices[ntotal].deviceName, "%s", port_str);
				}

				pstring = NULL;
				if (g_object_class_find_property(G_OBJECT_GET_CLASS(videosource), "guid")) {
					g_object_get(G_OBJECT(videosource), "guid", &deviceURI, NULL);
					devices[ntotal].deviceURI = deviceURI;
				} else {
					devices[ntotal].deviceURI = 0;
				}

				pstring = NULL;
				
				// Probe next device...
				ntotal++;
				
				if (ntotal >= PSYCH_MAX_VIDSRC) {
					printf("PTB-WARNING: Maximum number of allowable video sources during enumeration %i exceeded! Aborting enumeration.\n", PSYCH_MAX_VIDSRC);
					break;
				}
			}

			g_value_array_free(viddevs);
		}
		else {
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Video plugin '%s' doesn't provide any video devices to enumerate.\n", class_str);
		}
	}
	else {
		if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Video plugin '%s' doesn't support probing. No video devices to enumerate for this plugin.\n", class_str);
	}

	// Release videosource:
	gst_element_set_state(videosource, GST_STATE_NULL);
	gst_object_unref(GST_OBJECT(videosource));
	
	return;
}

/* PsychGSEnumerateVideoSources(int outPos, int deviceIndex);  -- Internal.
 *
 * Enumerates all connected and supported video sources into an internal
 * array "devices".
 *
 * If deviceIndex >= 0 : Returns pointer to PsychVideosourceRecordType struct
 *                       with info about the detected device with index 'deviceIndex'
 *                       or NULL if no such device exists. The pointer is valid until
 *                       the Screen module returns control to the runtime - then it
 *                       will get deallocated and must not be accessed anymore!
 *
 * If deviceIndex < 0 : Returns NULL to caller, returns a struct array to runtime
 *                      environment return argument position 'outPos' with all info
 *                      about the detected sources.
 */
PsychVideosourceRecordType* PsychGSEnumerateVideoSources(int outPos, int deviceIndex)
{
	PsychGenericScriptType 	*devs;
	const char *FieldNames[]={"DeviceIndex", "ClassIndex", "InputIndex", "ClassName", "InputHandle", "Device", "DevicePath", "DeviceName", "GUID", "DevicePlugin", "DeviceSelectorProperty" };

	int					i;
	GstElement			*videosource = NULL;
	GstPropertyProbe	*probe = NULL;
	GValueArray			*viddevs = NULL;
	GValue				*dev = NULL;
	char				*device_name = NULL;
	gchar               *pstring = NULL;
	PsychVideosourceRecordType *mydevice = NULL;
	
	// Make sure GStreamer is ready:
	PsychGSCheckInit("videocapture");

	// Allocate temporary space for enumerated devices:
	devices = (PsychVideosourceRecordType*) PsychCallocTemp(PSYCH_MAX_VIDSRC, sizeof(PsychVideosourceRecordType));
	ntotal  = 0;

	// Linux specific setup path:
	if (PSYCH_SYSTEM == PSYCH_LINUX) {
		// Try Video4Linux-II camera source:
		PsychGSEnumerateVideoSourceType("v4l2camsrc", 1, "Video4Linux2-CameraSource", "device", 0);

		// Try standard Video4Linux-II source:
		PsychGSEnumerateVideoSourceType("v4l2src", 2, "Video4Linux2", "device", 0);
	}

	if (PSYCH_SYSTEM == PSYCH_WINDOWS) {
		// Try Windows kernel streaming source:
		PsychGSEnumerateVideoSourceType("ksvideosrc", 1, "Windows WDM kernel streaming", "device-name", 0);

		// Use DirectShow to probe:
		PsychGSEnumerateVideoSourceType("dshowvideosrc", 2, "DirectShow", "device-name", 0);
	}
	
	if (PSYCH_SYSTEM == PSYCH_OSX) {
		// Try OSX Quicktime sequence grabber video source for 32-Bit systems with Quicktime-7:
		PsychGSEnumerateVideoSourceType("osxvideosrc", 1, "OSXQuicktimeSequenceGrabber", "device", 0);

		// Try OSX MIO video source: Unless we're under Octave, where some weird bug/interaction
		// would cause a crash in the miovideosrc plugin if we tried, so we don't try on Octave.
		#ifndef PTBOCTAVE3MEX
		PsychGSEnumerateVideoSourceType("miovideosrc", 3, "OSXMIOVideoSource", "device-name", 0);
		#endif

		// Try OSX QTKit video source for 64-Bit systems with Quicktime-X aka QTKit:
		PsychGSEnumerateVideoSourceType("qtkitvideosrc", 2, "OSXQuicktimeKitVideoSource", "device-index", 1);        
	}
	
	// Try DV-Cameras:
	PsychGSEnumerateVideoSourceType("dv1394src", 5, "DV1394", "guid", 0);
	
	// Try HDV-Cameras:
	PsychGSEnumerateVideoSourceType("hdv1394src", 6, "HDV1394", "guid", 0);
	
	// Try IIDC-1394 Cameras:
	// Does not work, no property probe interface:
	PsychGSEnumerateVideoSourceType("dc1394src", 7, "1394-IIDC", "camera-number", 1);

	// Try GeniCam-Cameras via aravis plugin:
	PsychGSEnumerateVideoSourceType("aravissrc", 8, "GeniCam-Aravis", "camera-name", 0);
    
	if (ntotal <= 0) {
		if (PsychPrefStateGet_Verbosity() > 4) {
            printf("PTB-INFO: Could not detect any supported video devices on this system.\n");
            printf("PTB-INFO: Trying to fake an auto-detected default device and a test video source...\n");
        }
        
        // Create a fake entry for the autovideosrc:
        devices[0].deviceIndex = 0;
        sprintf(devices[0].deviceVideoPlugin, "%s", "autovideosrc");
        sprintf(devices[0].deviceSelectorProperty, "%s", "");
        sprintf(devices[0].deviceHandle, "%s", "");
        
        // Create a fake entry for the videotestsrc:
        devices[1].deviceIndex = 1;
        sprintf(devices[1].deviceVideoPlugin, "%s", "videotestsrc");
        sprintf(devices[1].deviceSelectorProperty, "%s", "");
        sprintf(devices[1].deviceHandle, "%s", "");

        ntotal= 2;
	}

    // Add fake entry for deviceIndex zero, as a copy of the first real entry:
    memcpy(&devices[ntotal], &devices[0], sizeof(PsychVideosourceRecordType));
    devices[ntotal].deviceIndex = 0;
    ntotal++;
    
	// Have enumerated devices:
	if (deviceIndex >= 0) {
		// Yes: Return device name for that index:
		for (i=0; i < ntotal; i++) {
			if (devices[i].deviceIndex == deviceIndex) {
				mydevice = &devices[i];
				break;
			}
		}
	}
	else {
		// No: Code wants us to return struct array with all enumerated
		// devices to userspace:
		
		// Create output struct array with n output slots:
		PsychAllocOutStructArray(outPos, TRUE, ntotal, 11, FieldNames, &devs);
		
		// Iterate all available devices:
		for(i = 0; i < ntotal; i++) {
			PsychSetStructArrayDoubleElement("DeviceIndex", i, devices[i].deviceIndex, devs);
			PsychSetStructArrayDoubleElement("ClassIndex", i, devices[i].classIndex, devs);
			PsychSetStructArrayDoubleElement("InputIndex", i, devices[i].inputIndex, devs);
			PsychSetStructArrayStringElement("ClassName", i, devices[i].deviceClassName , devs);
			PsychSetStructArrayStringElement("InputHandle", i, devices[i].deviceHandle, devs);
			PsychSetStructArrayStringElement("Device", i, devices[i].device, devs);
			PsychSetStructArrayStringElement("DevicePath", i, devices[i].devicePath, devs);
			PsychSetStructArrayStringElement("DeviceName", i, devices[i].deviceName, devs);
			PsychSetStructArrayStringElement("GUID", i,  (devices[i].deviceURI > 0) ? devices[i].deviceHandle : "0", devs);
			PsychSetStructArrayStringElement("DevicePlugin", i, devices[i].deviceVideoPlugin, devs);
			PsychSetStructArrayStringElement("DeviceSelectorProperty", i, devices[i].deviceSelectorProperty, devs);
		}
	}

	// Done. Return device struct if assigned:
	return(mydevice);
}

psych_bool PsychGSGetResolutionAndFPSForSpec(PsychVidcapRecordType *capdev, int* width, int* height, double* fps, int reqdepth)
{
	GstCaps         *caps = NULL;
	GstStructure	*str;
	gint			qwidth, qheight;
	gint            qbpp;
	gint			rate1 = 0, rate2 = 1;
	gint			twidth = -1, theight = -1;
	gint            maxpixelarea = -1;
	double          tfps = 0.0;
	int             i;

	if (!usecamerabin) {
		// No camerabin, no way to query this stuff. Just fail
		printf("PTB-BUG: PsychGSGetResolutionAndFPSForSpec() called while !usecamerabin !?! Returning failure.\n");
		return(FALSE);
	}

	// Camerabin in use, we can actually enumerate and match.

	// Query caps of videosource and extract supported video capture modes:
    if (usecamerabin == 1) {
        // Camerabin1:
        g_object_get(G_OBJECT(capdev->camera), "video-source-caps", &caps, NULL);
    }
    else {
        // Camerabin2:
        g_object_get(G_OBJECT(capdev->camera), "viewfinder-supported-caps", &caps, NULL);        
    }

	if (caps) {
		if (PsychPrefStateGet_Verbosity() > 4)
			printf("PTB-DEBUG: Videosource caps are: %" GST_PTR_FORMAT "\n\n", caps);

		// Iterate through all supported video capture modes:
		for (i = 0; i < (int) gst_caps_get_size(caps); i++) {
			str = gst_caps_get_structure(caps, i);

            // Set a default of 1 pixel, in case query doesn't return anything:
            qwidth = 1;
			gst_structure_get_int(str, "width", &qwidth);
            qheight = 1;
			gst_structure_get_int(str, "height", &qheight);
			gst_structure_get_int(str, "bpp", &qbpp);
			gst_structure_get_fraction(str, "framerate", &rate1, &rate2);
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Videosource cap %i: w = %i h = %i\n", i, qwidth, qheight);

			// Is this detection of default resolution, or validation of a
			// given resolution?
			if ((*width == -1) && (*height == -1)) {
				// Auto-Detection of optimal default resolution.
				// Need to find the one with highest resolution (== pixel area):
				if (qwidth * qheight > maxpixelarea) {
					// A new favorite with max pixel area:
					maxpixelarea = qwidth * qheight;
					twidth = qwidth;
					theight = qheight;
					tfps = (double) rate1 / (double) rate2;
				}
			}
			else {
				// Validation: Reject/Skip modes which don't support requested resolution.
				if (((*width != (int) qwidth) && (qwidth > 1)) || ((*height != (int) qheight) && (qheight > 1))) continue;
                if ((qwidth == 1) && (qheight == 1)) continue;

				// Acceptable mode for requested resolution and framerate. Set it:
				maxpixelarea = qwidth * qheight;
                twidth = qwidth;
				theight = qheight;
				tfps = (double) rate1 / (double) rate2;
			}
		}

		gst_caps_unref(caps);

		// Any matching mode found?
		if (twidth == -1) {
			// No! The requested resolution + fps + pixelformat combo
			// is not supported:
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Could not validate video source resolution %i x %i. Returning failure.\n", *width, *height);
			return(FALSE);
		}

        // Special case DV video source, which has defined horizontal resolution of 720 pixels,
        // but no defined vertical resolution? Two possible values: 576 (PAL) or 480 (NTSC).
        if ((twidth == 720) && (theight == 1)) {
            // If auto-detection requested, then assume PAL is the correct choice,
            // ie., 576 pixels height. If just validation is requested, just pass-through
            // whatever was passed in if it was 480 or 576 pixels:
            if (*height == -1) {
                // Auto-Detect: Assume 576 pixels PAL:
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Special case auto-detect DV video source resolution. Defaulting to 720 x 576 PAL.\n");
                theight = 576;
            }
            else if (*height == 576 || *height == 480) {
                // Potentially valid height value: Accept.
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Special case DV video source resolution 720 x %i. Assuming ok.\n", *height);
                theight = *height;
            }
            else {
                // Impossible value: Reject.
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Could not validate DV video source resolution %i x %i. Returning failure.\n", *width, *height);
                return(FALSE);                
            }
        }
        
		// Yes. Return settings:
		*fps = tfps;
		*width = twidth;
		*height = theight;

		if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Will use auto-detected or validated video source resolution %i x %i.\n", *width, *height);

		return(TRUE);
	} else {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychGSGetResolutionAndFPSForSpec(): Capability query to video source failed! Returning failure.\n");
		return(FALSE);
	}
}

// Parse codecSpec string, check for certain element specs. If found, parse them
// and create corresponding GstElement*, otherwise return NULL.
GstElement* CreateGStreamerElementFromString(const char* codecSpec, const char* typeSpec, char* codecParams)
{
	char *codecPipelineSpec, *codecPipelineEnd;
	GstElement* element = NULL;

	if (strstr(codecSpec, typeSpec)) {
		// Find start of codec spec string:
		codecPipelineSpec = strstr(codecSpec, typeSpec);

		// Cut off the prefix:
		codecPipelineSpec += strlen(typeSpec);

		// Copy the remaining string, so we can mess with it:
		codecPipelineSpec = strdup(codecPipelineSpec);

		// Search for end-of-codecspec marker:
		codecPipelineEnd = strstr(codecPipelineSpec, ":::");

		// If any, null-terminate at start of marker:
		if (codecPipelineEnd) *codecPipelineEnd = 0;

		// codecPipelineSpec is now a null-terminated string which only
		// contains the GStreamer syntax description of the video codec and
		// its parameters.

		// Parse it and create a corresponding bin for use as encoder element:
		// Set GError* to NULL: Real men don't do error handling.
		element = gst_parse_bin_from_description((const gchar *) codecPipelineSpec, TRUE, NULL);
		if (element == NULL) {
			// Oopsie: Failed!
			if (PsychPrefStateGet_Verbosity() > 1) {
				printf("PTB-WARNING: Failed to create an encoder element of type '%s' from the following passed parameter string:\n", typeSpec);
				printf("PTB-WARNING: %s\n", codecPipelineSpec);
				printf("PTB-WARNING: Will revert to default settings for this element. This will likely fail soon...\n");
				printf("PTB-WARNING: Full parameter string was:\n");
				printf("PTB-WARNING: %s\n", codecSpec);
			}
		}
		else {
			// Success!
			if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Element '%s' created from spec '%s'.\n", typeSpec, codecPipelineSpec);
		}
		sprintf(codecParams, "%s", codecPipelineSpec);
		free(codecPipelineSpec);
	}

	return(element);
}

// If encoding profiles are supported, use this function to create one from given specs:
#ifdef PTB_USE_GSTENCODINGPROFILES
GstEncodingProfile* PsychCreateEncodingProfile(const char* muxer, const char* videocodec, const char* audiocodec)
{
    GstEncodingContainerProfile *prof = NULL;
    GstCaps *caps;
    
    // Create media container from 'muxer' spec:
    caps = gst_caps_from_string(muxer);
    prof = gst_encoding_container_profile_new("Psychtoolbox audio/video",
                                              "Psychtoolbox recorded audio and video",
                                              caps, muxer_preset_name);
    gst_caps_unref(caps);

    // Create video encoding profile from 'videocodec' spec:
    caps = gst_caps_from_string(videocodec);
    gst_encoding_container_profile_add_profile(prof, (GstEncodingProfile*) gst_encoding_video_profile_new(caps, video_preset_name, NULL, 0));
    gst_caps_unref(caps);

    // Create audio encoding profile from 'audiocodec' spec, if any. Otherwise omit audio encoding:
    if (audiocodec) {
        caps = gst_caps_from_string(audiocodec);
        gst_encoding_container_profile_add_profile(prof, (GstEncodingProfile*) gst_encoding_audio_profile_new(caps, audio_preset_name, NULL, 0));
        gst_caps_unref(caps);
    }
    
    return((GstEncodingProfile*) prof);
}
#endif

// Save all preset values of all plugins (encoders, muxers) contained in 'element' bin under
// the preset name 'preset_name':
void PsychSaveEncodingPresetsForElement(GstElement* element, const char* preset_name)
{
    GstElement* plugin;
    psych_bool done = FALSE;
    GstIterator* it = NULL;

    // Bin or single plugin?
    if (!GST_IS_BIN(element)) {
        // Handling of single elements:
        plugin = element;
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: In pipeline: %s presets for encoder plugin: %s\n", GST_IS_PRESET(plugin) ? "Saving" : "No", (const char*) gst_object_get_name(GST_OBJECT(plugin)));
        if (GST_IS_PRESET(plugin)) gst_preset_save_preset(GST_PRESET(plugin), preset_name);
        return;
    }
    
    // Handling of elements which are bins, ie., containers for other elements:
    it = gst_bin_iterate_all_by_interface(element, GST_TYPE_PRESET); // Retrieve all plugins in element with GstPreset interface:

    // Iterate over eligible plugins:
    while (it && !done) {
        plugin = NULL;
        switch (gst_iterator_next(it, (void**) &plugin)) {
            case GST_ITERATOR_OK:
                if (plugin) {
                    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: In pipeline: %s presets for encoder plugin: %s\n", GST_IS_PRESET(plugin) ? "Saving" : "No", (const char*) gst_object_get_name(GST_OBJECT(plugin)));
                    
                    // Double check the 'plugin' supports the GstPreset interface. Save its presets if so:
                    if (GST_IS_PRESET(plugin)) gst_preset_save_preset(GST_PRESET(plugin), preset_name);

                    gst_object_unref(plugin);
                }
                break;
                
            case GST_ITERATOR_RESYNC:
                gst_iterator_resync(it);
                break;
                
            case GST_ITERATOR_DONE:
                done = TRUE;
                break;
                
            default:
                if (plugin) gst_object_unref(plugin);
        }
    }
    
    if (it) gst_iterator_free(it);
    it = NULL;

    return;
}
    
psych_bool PsychSetupRecordingPipeFromString(PsychVidcapRecordType* capdev, char* codecSpec, char* outCodecName, psych_bool launchline, psych_bool forCamerabin, psych_bool soundForVideoRecording)
{
	GstElement *camera = NULL;
	GstElement *some_element = NULL;
	GstElement *audio_enc = NULL;
	GstElement *audio_src = NULL;
	GstElement *muxer_elt = NULL;

	char *poption = NULL;
	char *codecName = NULL;
	char *codecSep  = NULL;
	char muxerProfile[1000];
	char audioProfile[1000];
	char muxer[1000];
	char audiocodec[1000];
	char videocodec[1000];
	char codecoption[1000];
	char audiosrc[1000];
    char aactype[20];

	int nrAudioChannels;
	int audioFreq;

	int interlaced = -1;
	int keyFrameInterval = -1;
	float videoQuality = -1;
	float audioQuality = -1;
	int audioBitrate = -1;
	int videoBitrate = -1;
	int bigFiles = -1;
	int fastStart = -1;
	int indexItemsSec = -1;
	int profile = -1;

	// Use GStreamer encoding profiles instead of low-level do-it-yourself setup?
	psych_bool use_profiles = (forCamerabin && (usecamerabin == 2)) ? TRUE : FALSE;

	// Require setup of an audio-codec for audio recording or audio track writing?
	psych_bool use_audio = (soundForVideoRecording || strstr(codecSpec, "AddAudioTrack")) ? TRUE : FALSE;

    // Use of audio encoder?
    if (use_audio) {
        // Yes. Try if faac AAC encoder is available, choose it as default choice for
        // AAC encoding. If unavailable, choose ffenc_aac instead as fallback:
        audio_enc = gst_parse_bin_from_description("faac", TRUE, NULL);
        if (audio_enc) {
            gst_object_unref(G_OBJECT(audio_enc));
            audio_enc = NULL;
            sprintf(aactype, "faac");
        }
        else {
            sprintf(aactype, "ffenc_aac");
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: faac audio encoder not available. Using ffenc_aac instead.\n");
        }
    }
    
    memset(muxerProfile, 0, sizeof(muxerProfile));
	memset(muxer, 0, sizeof(muxer));
	memset(audiosrc, 0, sizeof(audiosrc));
	memset(audiocodec, 0, sizeof(audiocodec));
    memset(audioProfile, 0, sizeof(audioProfile));
	memset(videocodec, 0, sizeof(videocodec));

	// Get camera object - the camerabin1/2:
	camera = capdev->camera;

	// Parse and assign high-level properties if any:
    // Must switch numeric input format to C-Language style, so decimal point
    // gets correctly recognized!
    setlocale(LC_NUMERIC, "C");

	codecSep = strstr(codecSpec, "Profile=");
	if (codecSep) {
		sscanf(codecSep, "Profile=%i", &profile);
	}

	codecSep = strstr(codecSpec, "Interlaced=");
	if (codecSep) {
		sscanf(codecSep, "Interlaced=%i", &interlaced);
	}

	codecSep = strstr(codecSpec, "Keyframe=");
	if (codecSep) {
		sscanf(codecSep, "Keyframe=%i", &keyFrameInterval);
	}

	// Audio bit rate in kBit/Sec
	codecSep = strstr(codecSpec, "Audiobitrate=");
	if (codecSep) {
		sscanf(codecSep, "Audiobitrate=%i", &audioBitrate);
	}

	// Video bit rate in kBit/Sec
	codecSep = strstr(codecSpec, "Videobitrate=");
	if (codecSep) {
		sscanf(codecSep, "Videobitrate=%i", &videoBitrate);
	}

	codecSep = strstr(codecSpec, "Bigfiles=");
	if (codecSep) {
		sscanf(codecSep, "Bigfiles=%i", &bigFiles);
	}

	codecSep = strstr(codecSpec, "Faststart=");
	if (codecSep) {
		sscanf(codecSep, "Faststart=%i", &fastStart);
	}

	codecSep = strstr(codecSpec, "Timeresolution=");
	if (codecSep) {
		sscanf(codecSep, "Timeresolution=%i", &indexItemsSec);
	}

	codecSep = strstr(codecSpec, "Videoquality=");
	if (codecSep) {
		sscanf(codecSep, "Videoquality=%f", &videoQuality);
	}

	codecSep = strstr(codecSpec, "Audioquality=");
	if (codecSep) {
		sscanf(codecSep, "Audioquality=%f", &audioQuality);
	}

	// Create matching video encoder for codecSpec:
	// Codecs are sorted by suitability for high quality realtime video recording, so if usercode
	// doesn't specify a codec, and thereby requests use of the recommended default codec 'DEFAULTenc',
	// we will try to choose the best codec, then on failure fallback to the 2nd best, 3rd best, etc.

	// Video encoder from parameter string?
    some_element = NULL;
	if ((use_profiles && strstr(codecSpec, "VideoCodec=")) || ((some_element = CreateGStreamerElementFromString(codecSpec, "VideoCodec=", videocodec)) != NULL)) {
		// Yes. Assign it as our encoder:
		capdev->videoenc = some_element;
        if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            

		// Need to extract the actual name of the codec from the codecSpec string:
		codecName = strstr(codecSpec, "VideoCodec=");
		codecName+= strlen("VideoCodec=");
		codecName = strdup(codecName);
		codecSep = codecName;
		while ((*codecSep > 0) && (*codecSep != ' ')) codecSep++;
		*codecSep = 0;

		sprintf(outCodecName, "%s", codecName);
		free(codecName);
		codecName = NULL;
	}

	// Start with H264 encoder as default:
	if (strstr(codecSpec, "x264enc") || strstr(codecSpec, "1635148593") || (strstr(codecSpec, "DEFAULTenc") && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/mpeg,mpegversion=4");
        sprintf(muxerProfile, "video/x-msvideo");
        sprintf(outCodecName, "video/x-h264");
        sprintf(audiocodec, "AudioCodec=%s ", aactype); // Need to use faac MPEG-4 audio encoder.
        sprintf(muxer, "qtmux");                        // Need to use Quicktime-Multiplexer.
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=x264enc ");

			// Interlace flag specified?
			if (interlaced >= 0) {
				// Enable or disable handling of interlaced input video stream:
				sprintf(codecoption, "interlaced=%i ", (interlaced) ? 1 : 0);
				strcat(videocodec, codecoption);
			}

			// Keyframe interval specified?
			if (keyFrameInterval >= 0) {
				// Assign maximum distance between key frames:
				sprintf(codecoption, "key-int-max=%i ", keyFrameInterval);
				strcat(videocodec, codecoption);
			} else {
				// Default to a keyframe at least every 30 frames:
				strcat(videocodec, "key-int-max=30 ");
			}

			// Encoding profile specified?
			if (profile >= 0) {
				// Assign profile:
				sprintf(codecoption, "profile=%i ", profile);
				strcat(videocodec, codecoption);
			} else {
				// Default to "High" profile: 640 x 480 @ 30 fps possible:
				strcat(videocodec, "profile=3 ");
			}

			// Quality vs. Speed tradeoff specified?
			if (videoQuality >= 0) {
				// Yes: Map quality vs. speed scalar to 0-10 number for speed preset:
				if (videoQuality > 1) videoQuality = 1;

				// Speed-Quality profile: 0 = "None"
				// 1 = "Ultra fast", 2 = "Super fast", 3 = "Very fast",  4 = "Faster", 5 = "Fast",
				// 6 = "Medium" - the default,7 = "Slow", 8 = "Slower", 9 = "Very slow", 10 = "Placebo"
				sprintf(codecoption, "speed-preset=%i ", (int) (videoQuality * 10.0 + 0.5));
				strcat(videocodec, codecoption);
			}
			else {
				// No: Use the fastest speed at lowest quality:
				strcat(videocodec, "speed-preset=1 ");
			}
            
			// Bitrate specified?
			if (videoBitrate >= 0) {
				sprintf(codecoption, "bitrate=%i ", (int) videoBitrate);
				strcat(videocodec, codecoption);
			}

			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);

            // If specific codec couldn't be created but encoding profiles are in use then we "fake" a
            // 0x1 codec to suppress warning messages and error/fallback handling. The encoding profile
            // will take care of this, we just can't set specific encoding options:
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;
		}

		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'x264enc' H.264 video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "x264enc");
		}
	}

	// Try then Apple OSX specific H264 encoder:
	if (strstr(codecSpec, "vtenc_h264") || (strstr(codecSpec, "DEFAULTenc") && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/mpeg,mpegversion=4");
        sprintf(muxerProfile, "video/x-msvideo");
        sprintf(outCodecName, "video/x-h264");
        sprintf(audiocodec, "AudioCodec=%s ", aactype); // Need to use faac MPEG-4 audio encoder.
        sprintf(muxer, "qtmux");                 // Need to use Quicktime-Multiplexer.
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=vtenc_h264 ");
            
			// Bitrate specified?
			if (videoBitrate >= 0) {
				sprintf(codecoption, "bitrate=%i ", (int) videoBitrate);
				strcat(videocodec, codecoption);
			}
            
			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            
		}
        
		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'vtenc_h264' H.264 MacOSX specific video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "vtenc_h264");
		}
	}

	// Then xvidenc - MPEG4 in a AVI container: High quality, handles 640 x 480 @ 30 fps on
	// a 4 year old MacBookPro Core2Duo 2.2 Ghz with about 70% - 100% cpu load, depending on settings.
	if (strstr(codecSpec, "xvidenc") || strstr(codecSpec, "1836070006") || (strstr(codecSpec, "DEFAULTenc") && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/mpeg,mpegversion=4");
        sprintf(muxerProfile, "video/x-msvideo");
        sprintf(outCodecName, "video/x-xvid");
        sprintf(audiocodec, "AudioCodec=%s ", aactype); // Need to use faac MPEG-4 audio encoder.
        sprintf(muxer, "avimux");                // Need to use avi-multiplexer.
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=xvidenc ");

			// Interlace flag specified?
			if (interlaced >= 0) {
				// Enable or disable handling of interlaced input video stream:
				sprintf(codecoption, "interlaced=%i ", (interlaced) ? 1 : 0);
				strcat(videocodec, codecoption);
			}

			// Keyframe interval specified?
			if (keyFrameInterval >= 0) {
				// Assign maximum distance between key frames:
				sprintf(codecoption, "max-key-interval=%i ", keyFrameInterval);
				strcat(videocodec, codecoption);
			}

			// Quality vs. Speed tradeoff specified?
			if (videoBitrate >= 0) {
				sprintf(codecoption, "bitrate=%i ", (int) videoBitrate * 1024);
				strcat(videocodec, codecoption);
			}

			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            
		}

		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'xvidenc' xvid/mpeg-4 video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "xvidenc");
		}
	}

	// ffenc_mpeg4 also creates MPEG4, but at lower quality - much more blocky etc.
	if (strstr(codecSpec, "ffenc_mpeg4") || ((strstr(codecSpec, "DEFAULTenc") || strstr(codecSpec, "1836070006")) && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/mpeg,mpegversion=4");
        sprintf(muxerProfile, "video/x-msvideo");
        sprintf(outCodecName, "video/mpeg,mpegversion=4");
        sprintf(audiocodec, "AudioCodec=%s ", aactype); // Need to use faac MPEG-4 audio encoder.
        sprintf(muxer, "avimux");                // Need to use avi-multiplexer.
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=ffenc_mpeg4 ");

			// Interlace flag specified?
			if (interlaced >= 0) {
				// Enable or disable handling of interlaced input video stream:
				sprintf(codecoption, "interlaced=%i ", (interlaced) ? 1 : 0);
				strcat(videocodec, codecoption);
			}

			// Keyframe interval specified?
			if (keyFrameInterval >= 0) {
				// Assign maximum distance between key frames:
				sprintf(codecoption, "max-key-interval=%i ", keyFrameInterval);
				strcat(videocodec, codecoption);
			}

			// Quality vs. Speed tradeoff specified?
			if (videoBitrate >= 0) {
				sprintf(codecoption, "bitrate=%i ", (int) videoBitrate * 1024);
				strcat(videocodec, codecoption);
			}

			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            
		}

		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'ffenc_mpeg4' mpeg-4 video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "ffenc_mpeg4");
		}
	}

	// Theora + Ogg vorbis audio in .ogv container:
	if (strstr(codecSpec, "theoraenc") || (strstr(codecSpec, "DEFAULTenc") && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/x-vorbis");
        sprintf(muxerProfile, "application/ogg");
        sprintf(outCodecName, "video/x-theora");
        sprintf(audiocodec, "AudioCodec=vorbisenc "); // Need to use Ogg Vorbis audio encoder.
        sprintf(muxer, "oggmux");                // Need to use Ogg-multiplexer.
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=theoraenc drop-frames=0 speed-level=2 ");

			// Keyframe interval specified?
			if (keyFrameInterval >= 0) {
				// Assign maximum distance between key frames:
				sprintf(codecoption, "keyframe-auto=0 keyframe-force=%i ", keyFrameInterval);
				strcat(videocodec, codecoption);
			}

			// Quality vs. Speed tradeoff specified?
			if (videoQuality >= 0) {
				// Yes: Map quality vs. speed scalar to 0-63 number:
				if (videoQuality > 1) videoQuality = 1;

				sprintf(codecoption, "quality=%i ", (int) (videoQuality * 63.0 + 0.5));
				strcat(videocodec, codecoption);
			}

			// Quality vs. Speed tradeoff specified?
			if (videoBitrate >= 0) {
				sprintf(codecoption, "bitrate=%i ", (int) videoBitrate);
				strcat(videocodec, codecoption);
			}

			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            
		}

		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'theoraenc' video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "theoraenc");
		}
	}

	// VP-8 in WebM or Matroska container:
	if (strstr(codecSpec, "vp8enc") || (strstr(codecSpec, "DEFAULTenc") && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/x-vorbis");
        sprintf(muxerProfile, "video/webm");
        sprintf(outCodecName, "video/x-vp8");
        sprintf(audiocodec, "AudioCodec=vorbisenc "); // Need to use Ogg Vorbis audio encoder.
        sprintf(muxer, "webmmux");                    // Default to WebM multiplexer.
        
		// Need to use matroska/webm-multiplexer:
		if (strstr(codecSpec, "_matroska")) {
            sprintf(muxer, "matroskamux");
            sprintf(muxerProfile, "video/x-matroska");
        }
        
		if (strstr(codecSpec, "_webm")) {
            sprintf(muxer, "webmmux");
            sprintf(muxerProfile, "video/webm");
        }
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=vp8enc ");

			// Keyframe interval specified?
			if (keyFrameInterval >= 0) {
				// Assign maximum distance between key frames:
				sprintf(codecoption, "max-keyframe-distance=%i ", keyFrameInterval);
				strcat(videocodec, codecoption);
			}

			// Quality vs. Speed tradeoff specified?
			if (videoQuality >= 0) {
				// Yes: Map quality vs. speed scalar to 0-10 number:
				if (videoQuality > 1) videoQuality = 1;

				sprintf(codecoption, "quality=%i ", (int) (videoQuality * 10.0 + 0.5));
				strcat(videocodec, codecoption);
			}

			// Quality vs. Speed tradeoff specified?
			if (videoBitrate >= 0) {
				sprintf(codecoption, "bitrate=%i ", (int) videoBitrate * 1024);
				strcat(videocodec, codecoption);
			}

			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            
		}

		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'vp8enc' VP-8 video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "vp8enc");
		}
	}

	// Raw YUV: Works with camerabin1, but only for fixed hard-coded size 640 x 480 by now. Always works with camerabin2.
	if (strstr(codecSpec, "yuvraw") || (strstr(codecSpec, "DEFAULTenc") && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/mpeg,mpegversion=4");
        sprintf(muxerProfile, "video/x-msvideo");
        sprintf(outCodecName, "video/x-raw-yuv");
        sprintf(audiocodec, "AudioCodec=%s ", aactype); // Need to use faac MPEG-4 audio encoder.
        sprintf(muxer, "avimux");                // Need to use AVI-Multiplexer.
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=capsfilter caps=\"video/x-raw-yuv, format=(fourcc)I420, width=(int)640, height=(int)480\" ");

			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            
		}

		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'capsfilter' YUV pass-through video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "yuvraw");
		}
	}

	// H263 -- Does not work well yet.
	if (strstr(codecSpec, "ffenc_h263p") || strstr(codecSpec, "1748121139") || (strstr(codecSpec, "DEFAULTenc") && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/mpeg,mpegversion=4");
        sprintf(muxerProfile, "video/x-msvideo");
        sprintf(outCodecName, "video/x-h263");
        sprintf(audiocodec, "AudioCodec=%s ", aactype); // Need to use faac MPEG-4 audio encoder.
        sprintf(muxer, "qtmux");                 // Need to use Quicktime-Multiplexer.
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=ffenc_h263p ");

			// Interlace flag specified?
			if (interlaced >= 0) {
				// Enable or disable handling of interlaced input video stream:
				sprintf(codecoption, "interlaced=%i ", (interlaced) ? 1 : 0);
				strcat(videocodec, codecoption);
			}

			// Keyframe interval specified?
			if (keyFrameInterval >= 0) {
				// Assign maximum distance between key frames:
				sprintf(codecoption, "max-key-interval=%i ", keyFrameInterval);
				strcat(videocodec, codecoption);
			}

			// Quality vs. Speed tradeoff specified?
			if (videoBitrate >= 0) {
				sprintf(codecoption, "bitrate=%i ", (int) videoBitrate * 1024);
				strcat(videocodec, codecoption);
			}

			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            
		}

		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'ffenc_h263p' H.263 video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "ffenc_h263p");
		}
	}

	// Raw Huffman encoded YUV: Works with camerabin2.
	if (strstr(codecSpec, "huffyuv") || (strstr(codecSpec, "DEFAULTenc") && !capdev->videoenc)) {
		// Define recommended (compatible) audioencoder/muxer and their default options:
        sprintf(audioProfile, "audio/mpeg,mpegversion=4");
        sprintf(muxerProfile, "video/x-msvideo");
        sprintf(outCodecName, "video/x-huffyuv");
        sprintf(audiocodec, "AudioCodec=%s ", aactype); // Need to use faac MPEG-4 audio encoder.
        sprintf(muxer, "avimux");           // Need to use AVI-Multiplexer.
        
		// Videoencoder not yet created? If so, we have to do it now:
		if (!capdev->videoenc) {
			// Not yet created. Create full codec & option string from high level properties,
			// if any, then create based on string:
			sprintf(videocodec, "VideoCodec=ffenc_huffyuv ");

			// Keyframe interval specified?
			if (keyFrameInterval >= 0) {
				// Assign maximum distance between key frames:
				sprintf(codecoption, "gop-size=%i ", keyFrameInterval);
				strcat(videocodec, codecoption);
			}

			// Quality vs. Speed tradeoff specified?
			if (videoBitrate >= 0) {
				sprintf(codecoption, "bitrate=%i ", (int) videoBitrate * 1024);
				strcat(videocodec, codecoption);
			}

			// Create videocodec from options string:
			capdev->videoenc = CreateGStreamerElementFromString(videocodec, "VideoCodec=", videocodec);
            if (!capdev->videoenc && use_profiles) capdev->videoenc = (GstElement*) 0x1;            
		}

		if (!capdev->videoenc) {
			printf("PTB-WARNING: Failed to create 'ffenc_huffyuv' compressed YUV lossless video encoder! Does not seem to be installed on your system?\n");
		}
		else {
			if (!use_profiles) sprintf(outCodecName, "ffenc_huffyuv");
		}
	}

    // Use of audio encoding for any means requested? If so, perform codec and audio source setup:
    if (use_audio) {
        // Audio encoder from parameter string? Or new-style encoding profiles in use?

        // Encoding profiles in use? Then we need to parse optional AudioProfile= string ourselves for the MIME spec,
        // unless none specified. In that case audioProfile will already contain the auto-setup MIME string from
        // video codec setup above.
        if (use_profiles && strstr(codecSpec, "AudioProfile=")) {
            poption = strstr(codecSpec, "AudioProfile=");
            poption+= strlen("AudioProfile=");
            
            // Store audiocodec name in 'audiocodec':
            sprintf(audioProfile, "%s", poption);
            
            // Terminate audiocodec string if a ::: marker is encountered:
            poption = strstr(audioProfile, ":::");
            if (poption) *poption = 0;
            
            // audioProfile now contains audio profile MIME spec.
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Selecting audio encoder according to: %s\n", audioProfile);
        }
        
        if ((audio_enc = CreateGStreamerElementFromString(codecSpec, "AudioCodec=", audiocodec)) != NULL) {
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Created audio encoder according to: %s\n", audiocodec);
        } else {
            // No: Build from preset as defined by video codec and some high level settings:
            // This path is hit if encoding-profiles aren't used, so we need to do our own
            // codec setup (ie., camerabin1 or classic movie writing) and usercode didn't
            // specify a codec and its settings manually via "AudioCodec=" parameter.
            
            // Audio quality flag specified?
            if (audioQuality >= 0) {
                if (audioQuality > 1) audioQuality = 1;
                
                if (strstr(audiocodec, "faac") || strstr(audiocodec, "ffenc_aac")) {
                    // Map quality range 0.0 - 1.0 to bitrate range 0 - 320000 bits/sec:
                    sprintf(codecoption, "bitrate=%i ", (int) (audioQuality * 320000));
                    strcat(audiocodec, codecoption);
                }
                
                if (strstr(audiocodec, "vorbisenc")) {
                    // Assign quality range 0.0 - 1.0 directly:
                    sprintf(codecoption, "quality=%f ", audioQuality);
                    strcat(audiocodec, codecoption);
                }
            }
            
            // Audio bitrate specified?
            if (audioBitrate >= 0) {
                if (strstr(audiocodec, "faac") || strstr(audiocodec, "ffenc_aac")) {
                    sprintf(codecoption, "bitrate=%i ", audioBitrate * 1000);
                    strcat(audiocodec, codecoption);
                }
                
                if (strstr(audiocodec, "vorbisenc")) {
                    sprintf(codecoption, "managed=1 bitrate=%i ", audioBitrate * 1000);
                    strcat(audiocodec, codecoption);
                }
            }
            
            // Create audio encoder:
            if ((audio_enc = CreateGStreamerElementFromString(audiocodec, "AudioCodec=", audiocodec)) != NULL) {
                // Yes: Assign it.
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Created audio encoder according to: %s\n", audiocodec);
            } else {
                if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Failed to create requested audio encoder [%s]! Falling back to default encoder%s.\n", audiocodec, (use_profiles) ? " for profile" : "");
            }
        }
        
        // Audio source from parameter string? Most of the time usercode would not want
        // to change type/assignment/source parameters of a source. Also, the available
        // settings are common enough across sources that it doesn't make sense for us to
        // introduce our own high-level parameters. Therefore we only provide the low-level
        // interface here.
        // Interesting sources would be: alsasrc, osssrc, oss4src, pulsesrc, jackaudiosrc,
        // autoaudiosrc, gconfaudiosrc.
        //
        // Interesting settings would be:
        // device = The audio device name.
        // server = Audio server name for Jack and PulseAudio.
        // slave-method = Type of syncing to master clock.
        if (soundForVideoRecording) audio_src = CreateGStreamerElementFromString(codecSpec, "AudioSource=", audiosrc);
    }
    else {
        // No audio encoding/writing: Disable all audio related stuff.
        audio_src = NULL;
        audio_enc = NULL;
        audiocodec[0] = 0;
        audioProfile[0] = 0;
    }
    
	// Multiplexer from parameter string?
	if (strstr(codecSpec, "Muxer=")) {
		// Must create it without help of CreateGStreamerElementFromString() as the muxer has
		// 2 sink pads, but our method can only handle 1 sink pad. Therefore just create a
		// muxer by name, setting it to its default settings:
		poption = strstr(codecSpec, "Muxer=");
		poption+= strlen("Muxer=");

		// Store muxer name in 'muxer':
		sprintf(muxer, "%s", poption);

		// Terminate muxer string if a ::: marker is encountered:
		poption = strstr(muxer, ":::");
		if (poption) *poption = 0;
        
        // muxer now contains muxer spec.
	}

    // Muxer MIME profile specified?
	if (strstr(codecSpec, "MuxerProfile=")) {
		poption = strstr(codecSpec, "MuxerProfile=");
		poption+= strlen("MuxerProfile=");
        
		// Store muxer name in 'muxer':
		sprintf(muxerProfile, "%s", poption);
        
		// Terminate muxer string if a ::: marker is encountered:
		poption = strstr(muxerProfile, ":::");
		if (poption) *poption = 0;
        
        // muxerProfile now contains container spec.
	}
    
    // Build multiplexer from 'muxer' spec:
    muxer_elt = gst_element_factory_make(muxer, "ptbvideomuxer0");

	// Special muxer-specific setup of mux parameters:
	if (strstr(muxer, "avimux")) {
		// Big file writing support (> 2GB) is on by default. This allows to change it:
		if (bigFiles >= 0) {
			g_object_set(muxer_elt, "bigfile", (bigFiles > 0) ? 1 : 0, NULL);
			sprintf(codecoption, " bigfile=%i", (bigFiles > 0) ? 1 : 0);
			strcat(muxer, codecoption);
		}
	}
	
	if (strstr(muxer, "qtmux")) {
		// Big file writing support (> 2GB) is on by default. This allows to change it:
		if (bigFiles >= 0) {
			g_object_set(muxer_elt, "large-file", (bigFiles > 0) ? 1 : 0, NULL);
			sprintf(codecoption, " large-file=%i", (bigFiles > 0) ? 1 : 0);
			strcat(muxer, codecoption);
		}

		if (fastStart >= 0) {
			// Enable/Disable fast start support for low-latency start of movie playback:
			g_object_set(muxer_elt, "faststart", (fastStart > 0) ? 1 : 0, NULL);
			sprintf(codecoption, " faststart=%i", (fastStart > 0) ? 1 : 0);
			strcat(muxer, codecoption);
		}
		else {
			// Enable fast start support for low-latency start of movie playback:
            // FIXME: Unless we're on Windows where this would cause a crash in qtmux!
            // TODO: Investigate & Fix root cause. For now we just disable faststart
            // on Windows.
            fastStart = (PSYCH_SYSTEM == PSYCH_WINDOWS) ? 0 : 1;
			g_object_set(muxer_elt, "faststart", fastStart, NULL);
			sprintf(codecoption, " faststart=%i", fastStart);
			strcat(muxer, codecoption);
		}

		if (indexItemsSec >= 0) {
			// Write seek index entry for each 1/indexItemsSec of movie time:
			g_object_set(muxer_elt, "movie-timescale", indexItemsSec, NULL);
			sprintf(codecoption, " movie-timescale=%i", indexItemsSec);
			strcat(muxer, codecoption);
		} else {
			// Write seek index entry for each millisecond of movie time:
			g_object_set(muxer_elt, "movie-timescale", 1000, NULL);
			sprintf(codecoption, " movie-timescale=%i", 1000);
			strcat(muxer, codecoption);
		}
	}

	if (strstr(muxer, "matroskamux") || strstr(muxer, "webmmux")) {
		if (indexItemsSec >= 0) {
			// Write seek index entry for each 1/indexItemsSec of movie time:
			g_object_set(muxer_elt, "min-index-interval", (int)  (1e9 / indexItemsSec), NULL);
			sprintf(codecoption, " min-index-interval=%i", (int) (1e9 / indexItemsSec));
			strcat(muxer, codecoption);
		} else {
			// Write seek index entry for each millisecond (1e6 nanoseconds) of movie time:
			g_object_set(muxer_elt, "min-index-interval", (int)  1e6, NULL);
			sprintf(codecoption, " min-index-interval=%i", (int) 1e6);
			strcat(muxer, codecoption);
		}
	}

	// Still no video codec? Then this is game over: FIXME Small memory leak here on error exit.
	if (capdev->videoenc == NULL) PsychErrorExitMsg(PsychError_user, "Could not find or setup requested video codec or any fallback codec for video recording. Aborted.");
	if (use_audio && (!audio_enc && !use_profiles)) PsychErrorExitMsg(PsychError_user, "Could not find or setup requested audio codec for recording. Aborted.");
	if (!muxer_elt && !use_profiles) PsychErrorExitMsg(PsychError_user, "Could not find or setup requested audio-video multiplexer for recording. Aborted.");

    // Invalidate our no longer needed fake videoenc entry if profiles are in use:
    if (use_profiles && (capdev->videoenc == (GstElement*) 0x1)) capdev->videoenc = NULL;
    
	if (PsychPrefStateGet_Verbosity() > 3) {
		if (use_audio) printf("PTB-INFO: Audiosource: %s\n", audiosrc);
		if (use_audio) printf("PTB-INFO: Audiocodec : %s\n", audiocodec);
        if (use_audio && use_profiles) printf("PTB-INFO: Audio profile: %s\n", audioProfile);
		printf("PTB-INFO: Videocodec : %s\n", videocodec);
        if (use_profiles) printf("PTB-INFO: Video profile: %s\n", outCodecName);
		printf("PTB-INFO: Multiplexer: %s\n", muxer);
        if (use_profiles) printf("PTB-INFO: Container profile: %s\n", muxerProfile);
	}

	// Setup for video recording via generated codecs or encoding profiles for use with camerabin-1/2?
	if (!launchline || forCamerabin) {
        if (usecamerabin == 1) {
            // Attach our created objects to camerabin1:
            g_object_set(camera, "video-encoder", capdev->videoenc, NULL);
            if (audio_enc) g_object_set(camera, "audio-encoder", audio_enc, NULL);
            g_object_set(camera, "video-muxer", muxer_elt, NULL);
            if (audio_src) g_object_set(camera, "audio-source", audio_src, NULL);
        }
        else {
            // Camerabin2 in use:
            // Release our objects, if any. But first store their parameter settings
            // as presets for use by our video/audio/mux profiles.
            if (muxer_elt) {
                PsychSaveEncodingPresetsForElement(muxer_elt, muxer_preset_name);
                gst_object_unref(G_OBJECT(muxer_elt));
                muxer_elt = NULL;
            }
            
            if (audio_enc) {
                PsychSaveEncodingPresetsForElement(audio_enc, audio_preset_name);
                gst_object_unref(G_OBJECT(audio_enc));
                audio_enc = NULL;
            }
            
            if (capdev->videoenc) {
                // Need to black-list certain encoders, for which preset saving is buggy and would crash:
                if (!strstr(videocodec, "x264enc")) {
                    // Not black-listed - Save presets:
                    PsychSaveEncodingPresetsForElement(capdev->videoenc, video_preset_name);
                }
                else {
                    // Black listed:
                    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Saving of encoder presets for black-listed video encoder '%s' skipped.\n", videocodec);
                }
                gst_object_unref(G_OBJECT(capdev->videoenc));
                capdev->videoenc = NULL;
            }
            
            // Assign special audio source, if any:
            if (audio_src) g_object_set(camera, "audio-source", audio_src, NULL);
            
            // Assign video (and audio) encoding and muxing profile:
            #ifdef PTB_USE_GSTENCODINGPROFILES
            g_object_set(camera, "video-profile", PsychCreateEncodingProfile(muxerProfile, outCodecName, (soundForVideoRecording) ? audioProfile : NULL), NULL);
            #else
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-WARNING: Using default (low-performance) Ogg Theora/Vorbis video/audio encoding for camerabin2 video\n");
                printf("PTB-WARNING: recording, because your GStreamer installation lacks support for selectable encoding profiles.\n");
            }
            #endif
        }
	} else {
		// Setup for launch-line based movie writing : Release our objects, we have our launch line:
		if (capdev->videoenc) { gst_object_unref(G_OBJECT(capdev->videoenc)); capdev->videoenc = NULL; }
		if (audio_enc) gst_object_unref(G_OBJECT(audio_enc));
		if (audio_src) gst_object_unref(G_OBJECT(audio_src));
		if (muxer_elt) gst_object_unref(G_OBJECT(muxer_elt));

		// Build gst-launch style GStreamer pipeline spec string:

		// Special keyword for audio tracks in written movies?
		if ((poption=strstr(codecSpec, "AddAudioTrack"))) {
			// Audio and Video:
			poption+=strlen("AddAudioTrack");
			if (sscanf(poption, "=%i@%i", &nrAudioChannels, &audioFreq) != 2) {
				// Assign default 1 mono-channel count and 48 kHz frequency:
				nrAudioChannels = 1;
				audioFreq = 48000;
			}

			if (nrAudioChannels < 1 || nrAudioChannels > 256 || audioFreq < 1 || audioFreq > 200000) {
				printf("PTB-ERROR: Invalid parameters in 'AddAudioTrack='! Either %i audio channels outside valid range 1-256, or %i audioFreq outside allowable range 1 - 200000 Hz.\n",
				       nrAudioChannels, audioFreq);
				PsychErrorExitMsg(PsychError_user, "Invalid audio recording parameters provided");
			}

			// If no "AudioSource=" was provided then configure for appsrc feeding
			// from Screen('AddAudioBufferToMovie') into movie:
			if (strlen(audiosrc) == 0) sprintf(audiosrc, "appsrc name=ptbaudioappsrc do-timestamp=0 stream-type=0 max-bytes=0 block=1 is-live=0 emit-signals=0 caps=\"audio/x-raw-float, endianness=(int)1234, width=(int)32, channels=(int)%i, rate=(int)%i\" ! audioconvert ! queue", nrAudioChannels, audioFreq);

			// We add bits to feed from 'audiosrc' into 'audiocodec' into the common muxer of video and audio stream:
			sprintf(outCodecName, " %s ! ptbvideomuxer0. %s ! %s ! ptbvideomuxer0. %s name=ptbvideomuxer0 ", videocodec, audiosrc, audiocodec, muxer);            
		} else {
			// Video only:
            
			// We feed the output of 'videocodec' directly into the muxer:
			sprintf(outCodecName, " %s ! %s ", videocodec, muxer);
		}
        
		// Example launch line for audio+video recording:
		// gst-launch-0.10 -e v4l2src device=/dev/video0 ! videorate ! ffmpegcolorspace ! "video/x-raw-yuv, format=(fourcc)I420, width=(int)640, height=(int)480", framerate=15/1 ! queue ! x264enc speed-preset=1 ! muxout. autoaudiosrc ! faac ! muxout. avimux name=muxout ! filesink location=~/Desktop/record.avi
		// Example without audio, pure video:
		// gst-launch-0.10 -e v4l2src device=/dev/video0 ! videorate ! ffmpegcolorspace ! "video/x-raw-yuv, format=(fourcc)I420, width=(int)640, height=(int)480", framerate=15/1 ! queue ! x264enc speed-preset=1 ! muxout. avimux name=muxout ! filesink location=~/Desktop/record.avi
		// Better:
		// gst-launch-0.10 -e v4l2src device=/dev/video0 ! videorate ! ffmpegcolorspace ! "video/x-raw-yuv, format=(fourcc)I420, width=(int)640, height=(int)480", framerate=15/1 ! queue ! x264enc speed-preset=1 ! avimux ! filesink location=~/Desktop/record.avi
		//
		// Raw data:
		// st-launch-0.10 -e v4l2src device=/dev/video0 ! videorate ! ffmpegcolorspace ! "video/x-raw-yuv, format=(fourcc)I420, width=(int)640, height=(int)480", framerate=15/1 ! queue ! identity ! avimux ! filesink location=~/Desktop/record.avi
	}

	return(TRUE);
}

/* PsychGetCodecLaunchLineFromString() - Helper function for GStreamer based movie writing.
 *
 * Take a 'codecSpec' string with user provided codec settings string, get
 * a proper gst-launch style segment for the codecs and muxers back.
 * launchString must be preallocated and of sufficient size (approx. 1Kb).
 * Returns TRUE on success, FALSE on failure.
 */
psych_bool PsychGetCodecLaunchLineFromString(char* codecSpec, char* launchString)
{
	// Need a fake capdev to pass to PsychSetupRecordingPipeFromString():
	PsychVidcapRecordType dummydev;
	memset(&dummydev, 0, sizeof(dummydev));

	// Pass off the job and just return its results. This is shared with encoder/muxer
	// creation, setup, configuration and testing for the video recording engine:
	return(PsychSetupRecordingPipeFromString(&dummydev, codecSpec, launchString, TRUE, FALSE, FALSE));
}

/*  
*      PsychGSOpenVideoCaptureDevice() -- Create a video capture object.
*
*      This function tries to open and initialize a connection to a camera
*      and returns the associated captureHandle for it.
*
*      slotid = Number of slot in vidcapRecordBANK[] array to use for this camera.
*      win = Pointer to window record of associated onscreen window.
*      deviceIndex = Index of the grabber device.
*      capturehandle = handle to the new capture object.
*      capturerectangle = If non-NULL a ptr to a PsychRectangle which contains the ROI for capture. Special roi [0 0 w h] selects resolution w x h on device.
*      reqdepth = Number of layers for captured output textures. (0=Don't care, 1=LUMINANCE8, 2=LUMINANCE8_ALPHA8, 3=RGB8, 4=RGBA8, 5=YCBCR)
*      num_dmabuffers = Number of buffers to queue internally before dropping buffers. Zero = Don't drop buffers, fill up whole memory if neccessary.
*      allow_lowperf_fallback = If set to 1 then PTB can use less capable fallback path on setups which don't support the 'camerabin' plugin.
*      targetmoviefilename = Filename of movie file to record (if any) if deviceIndex >=0,
*                            special spec string for a video capture device if deviceIndex < 0.
*      recordingflags = Special flags to control internal operation (harddisc recording yes/no, live video during recording?, audio recording? workarounds?)
*/
psych_bool PsychGSOpenVideoCaptureDevice(int slotid, PsychWindowRecordType *win, int deviceIndex, int* capturehandle, double* capturerectangle,
								   int reqdepth, int num_dmabuffers, int allow_lowperf_fallback, char* targetmoviefilename, unsigned int recordingflags)
{
    GstCaps         *colorcaps;
    GstCaps         *vfcaps, *reccaps;
    GstElement      *camera = NULL;
    GMainLoop       *VideoContext = NULL;
    GstBus          *bus = NULL;
    GstElement      *videosink = NULL;
    GstElement      *videosource = NULL;
    GstElement      *videosource_filter = NULL;
    GstElement      *videocrop_filter = NULL;
    GstElement      *videowrappersrc = NULL;
    GstPad          *pad, *peerpad;
    GstPad          *dvpad;
    GstCaps         *caps;
    GstStructure    *str;
    gint            width, height;
    gint            rate1, rate2;
    gint            twidth, theight;
    char            *codecSpec;
    char            codecName[10000];
    PsychVidcapRecordType   *capdev = NULL;
    char            config[1000];
    char            tmpstr[1000];
    char            device_name[1000];
    char            plugin_name[1000];
    char            prop_name[1000];
    gchar           *pstring = NULL; 
    PsychVideosourceRecordType  *theDevice = NULL;
    psych_bool      overrideFrameSize = FALSE;

    config[0] = 0;
    tmpstr[0] = 0;
    device_name[0] = 0;
    plugin_name[0] = 0;
    prop_name[0] = 0;

	// Init capturehandle to none:
	*capturehandle = -1;
	
	// Make sure GStreamer is ready:
	PsychGSCheckInit("videocapture");
	
	// Map deviceIndex of requested video source to device name:
	if (deviceIndex >= 0) {
		// Get device name for given deviceIndex from video device
		// enumeration (or NULL if no such device):
		theDevice = PsychGSEnumerateVideoSources(-1, deviceIndex);
		if (NULL == theDevice) {
			printf("PTB-ERROR: There isn't any video capture device available for provided deviceIndex %i.\n", deviceIndex);
			PsychErrorExitMsg(PsychError_user, "Invalid deviceIndex provided. No such video source. Aborted.");
		}
		
		// Assign name:
		sprintf(device_name, "%s", theDevice->deviceHandle);
		sprintf(plugin_name, "%s", theDevice->deviceVideoPlugin);
		sprintf(prop_name, "%s", theDevice->deviceSelectorProperty);
		if (theDevice->deviceURI > 0) {
			sprintf(device_name, "%llu", theDevice->deviceURI);
		}

		if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Trying to open video capture device with deviceIndex %i [%s].\n", deviceIndex, device_name);
	}

    // Slot 'slotid' will contain the record for our new capture object:

    // Zero-out new record in video bank:
    memset(&vidcapRecordBANK[slotid], 0, sizeof(PsychVidcapRecordType));
    
    // Initialize new record:
    vidcapRecordBANK[slotid].valid = 1;
    
    // Retrieve device record for slotid:
    capdev = PsychGetGSVidcapRecord(slotid);
	
    capdev->camera = NULL;
    capdev->grabber_active = 0;
    capdev->scratchbuffer = NULL;        
    
    // Selection of pixel depths:
    if (reqdepth == 4 || reqdepth == 0) {
	    // Default is RGB 32 bit, aka RGBA8:
	    reqdepth = 4;
    }
    else {
        // Only other supported format is RGB24 bit:
        switch (reqdepth) {
            case 2:
                // A no-go: Instead we use 1 channel luminance8:
                if (PsychPrefStateGet_Verbosity()>1)
                    printf("PTB-WARNING: Video capture engine doesn't support requested Luminance+Alpha format. Will revert to pure luminance instead...\n");
                reqdepth = 1;
            break;

            case 1: // Accept as is: L8   aka Luminance 8 bit.
            case 3: // Accept as is: RGB8 aka RGB 24 bit.
            break;

            case 5: // Accept as YVYU.
                if (!(win->gfxcaps & kPsychGfxCapUYVYTexture)) {
                    // Usercode requested type 5 - UYVY textures, but GPU does not support them.
                    // Fallback to type 4 - RGBA8 textures:
                    reqdepth = 4;
                    if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Requested YUYV-422 texture format for video capture unsupported by GPU. Falling back to RGBA8 format.\n");
                }
            break;

            case 6: // Accept as YUV-I420
                if (!(win->gfxcaps & kPsychGfxCapFBO) || !PsychAssignPlanarI420TextureShader(NULL, win)) {
                    // Usercode requested type 6 - I420 textures, but GPU does not support them.
                    // Fallback to type 4 - RGBA8 textures:
                    reqdepth = 4;
                    if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Requested YUV-I420 texture format for video capture unsupported by GPU. Falling back to RGBA8 format.\n");
                }
            break;

            default:
            // Unknown format:
            PsychErrorExitMsg(PsychError_user, "You requested an invalid image depths (not one of 0, 1, 2, 3, 4, 5 or 6). Aborted.");
        }
    }

    // Requested output texture pixel depth in layers:
    capdev->reqpixeldepth = reqdepth;
    capdev->pixeldepth = reqdepth * 8;

    // Assign number of dma buffers to use:
    capdev->num_dmabuffers = num_dmabuffers;

    PsychInitMutex(&vidcapRecordBANK[slotid].mutex);
    PsychInitCondition(&vidcapRecordBANK[slotid].condition, NULL);
    
    // Try to open and initialize camera according to given settings:
    // Create video capture pipeline with camerabin plugin, unless this
    // is forcefully disabled by setting the env variable "PSYCH_FORCE_CAMERABIN2":
    usecamerabin = 1;
    if (!getenv("PSYCH_FORCE_CAMERABIN2")) camera = gst_element_factory_make ("camerabin", "ptbvideocapturepipeline");

    if (NULL == camera) {
        // Failed to create camerabin plugin. Retry with new camerabin2:
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Failed to create camerabin plugin. Retrying with camerabin2 ...\n");
        camera = gst_element_factory_make ("camerabin2", "ptbvideocapturepipeline");
        usecamerabin = 2;

        if (camera) {
            // Set basic opmode to 2 for video capture/recording, instead of default 1
            // for still image capture:
            g_object_set(G_OBJECT(camera), "mode", 2, NULL);
            
            // Workaround the continued brokeness of GStreamer-SDK, as of 2013.6. We create a dummy image/jpeg MIME
            // encoding profile and assign it for use with still image capture, despite the fact that we never use
            // still image capture. However this is neccessary because by default, camerabin2 has a image/jpeg container
            // profile assigned which doesn't work due to apparently missing jpeg muxer plugins (jifmux?) in the SDK. Our
            // MIME profile is valid but only defines image encoding, no container, so no muxer plugin needed -> problem solved.
            // I don't know if this would work for actual still image capture, but our workarounds purpose is only to prevent
            // error handling/error abort during camerabin2 init from kicking in, so we don't care:
            #ifdef PTB_USE_GSTENCODINGPROFILES
            caps = gst_caps_new_simple("image/jpeg", NULL);
            // We probably leak the created GstEncodingProfile here, but will be a small leak...
            g_object_set(G_OBJECT(camera), "image-profile", gst_encoding_video_profile_new(caps, NULL, NULL, 0), NULL);
            gst_caps_unref(caps); caps = NULL;
            #endif
        }
    }
    
    sprintf(config, "%s", device_name);

    // Camerabin disabled or creation failed?
    if (NULL == camera) {
		// Failed or disabled: Use fallback playbin2 implementation.
		usecamerabin = 0;
		
		if (!allow_lowperf_fallback)
			PsychErrorExitMsg(PsychError_user, "Failed to create high-performance video capture pipeline and script doesn't allow fallback! Aborted.");
		
		camera = gst_element_factory_make ("playbin2", "ptbvideocapturepipeline");
		
		// Assign a specific input video source URI name if possible:
		if (theDevice && (theDevice->classIndex  < 5)) sprintf(config, "v4l2://%s", device_name);
		// DV1394 and HDV1394 only encode the port property in their URI, not the GUID. Therefore
		// we can not select a camera by GUID if we have to use playbin2 instead of camerabin.
		// Do the best we can do: Select the generic "default camera" URI for now:
		if (theDevice && (theDevice->classIndex == 5)) sprintf(config, "dv://");
		if (theDevice && (theDevice->classIndex == 6)) sprintf(config, "hdv://");
		
		if (PsychPrefStateGet_Verbosity() > 1) {
			printf("PTB-WARNING: Could not use GStreamer 'camerabin' or 'camerabin2' plugins for videocapture. Will use 'playbin2' as fallback. Most\n");
			printf("PTB-WARNING: features, e.g., video recording, ROI and video resolution selection, are not supported in fallback mode.\n");
		}
    }
	
    // Pipeline creation failed?
    if (NULL == camera) PsychErrorExitMsg(PsychError_user, "Failed to create video capture pipeline! Fallback capture engine is also unsupported. Aborted.");

    // Assign new record in videobank:
    vidcapRecordBANK[slotid].camera = camera;
    vidcapRecordBANK[slotid].frameAvail = 0;

    // Enforce use of the system clock for this pipeline instead of leaving it to the pipeline
    // to choose a proper clock automatically:
    // TODO FIXME: To use or not to use? Preliminare results on Linux indicate it
    // doesn't make a difference, so skip it for now: gst_pipeline_use_clock(camera, gst_system_clock_obtain());

    // Need to select a video source for this capture pipeline:
    if (!usecamerabin) {
		// Fallback path with playbin2: All video source parameters are encoded in "URI"
		// property as a string:
		
		// Assign name and configuration parameters of video capture device to open:
		g_object_set(G_OBJECT(camera), "uri", config, NULL);
    }
    else {
		// High performance path with camerabin. Build an appropriate video source and set up its input:
		videosource = NULL;
		
		if (deviceIndex >= 0) {
			// Create proper videosource plugin if possible:
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach '%s' as video source...\n", plugin_name);
            
            if (strstr(plugin_name, "dv1394src")) {
                // dv1394src and hdv1394src need special treatment:
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device with guid '%llu' as video input [Property %s].\n", theDevice->deviceURI, prop_name);                
                sprintf(tmpstr, "%s typefind=true do-timestamp=true guid=%llu ! dvdemux ! dvdec name=ptbdvsource", plugin_name, theDevice->deviceURI);
                videosource = gst_parse_bin_from_description((const gchar *) tmpstr, FALSE, NULL);                
            }
            else {
                // Standard path:
                videosource = gst_element_factory_make(plugin_name, "ptb_videosource");
            }
            
			if (!videosource) {
				if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed! We are out of options and will probably fail soon.\n");
				PsychErrorExitMsg(PsychError_system, "GStreamer failed to find a suitable video source! Game over.");
			}

			// Attach correct video input device to it:
            if (strstr(plugin_name, "dv1394src")) {
                dvpad = gst_element_get_static_pad(gst_bin_get_by_name(GST_BIN(videosource), "ptbdvsource"), "src");
                gst_element_add_pad(videosource, gst_ghost_pad_new("src", dvpad));
                gst_object_unref(GST_OBJECT(dvpad));
                dvpad = NULL;                
            }
            else if ((!strcmp(plugin_name, "dc1394src") || !strcmp(plugin_name, "qtkitvideosrc")) && (prop_name[0] != 0)) {
				// DC1394 source or QTKITVideosource:
				if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device with guid '%llu' as video input [Property %s].\n", theDevice->deviceURI, prop_name);
				g_object_set(G_OBJECT(videosource), prop_name, (int) theDevice->deviceURI, NULL);
			} else {
				if ((config[0] != 0) && (prop_name[0] != 0) && (theDevice->deviceURI == 0)) {
					if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device '%s' as video input [Property %s].\n", config, prop_name);
					g_object_set(G_OBJECT(videosource), prop_name, config, NULL);
				}

				if ((theDevice->deviceURI > 0) && (prop_name[0] != 0)) {
					if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device with guid '%llu' as video input [Property %s].\n", theDevice->deviceURI, prop_name);
					g_object_set(G_OBJECT(videosource), prop_name, theDevice->deviceURI, NULL);
				}
			}           
		}
		
		// MS-Windows specific setup path:
		if ((PSYCH_SYSTEM == PSYCH_WINDOWS) && (deviceIndex > -8)) {
			if (deviceIndex < 0) {
				// Non-Firewire video source selected:
				
				//if ((deviceIndex >= 10000) || (deviceIndex == -1) || (deviceIndex == -2)) {
				if ((deviceIndex == -1) || (deviceIndex == -2)) {
					// First try Kernel-Streaming based video source for low-latency capture:
					// if (deviceIndex >= 10000) deviceIndex = deviceIndex - 10000;
					if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach ksvideosrc as video source...\n");
					videosource = gst_element_factory_make("ksvideosrc", "ptb_videosource");
				}

				if (videosource) {
					// Kernel streaming video source:
					
					// Specific deviceIndex requested, instead of auto-select?
					if (deviceIndex < 0) {
						// Fetch optional targetmoviename parameter as name spec string:
						if (targetmoviefilename == NULL) PsychErrorExitMsg(PsychError_user, "You set 'deviceIndex' to a negative value, but didn't provide the required device name string in the 'moviename' argument! Aborted.");

						// Assign:
						strcat(config, targetmoviefilename);
						
						switch(deviceIndex) {
							case -1:
								// Human friendly device name provided:
								if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device '%s' as video input.\n", config);
								g_object_set(G_OBJECT(videosource), "device-name", config, NULL);
								break;
								
							case -2:
								// DirectShow device path provided:
								if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device at path '%s' as video input.\n", config);
								g_object_set(G_OBJECT(videosource), "device-path", config, NULL);
								break;
								
							default:
								PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' provided! Not a value of -1 or -2. Aborted.");
						}
					}
				}

				// No kernel streaming video source available?
				if (!videosource) {
					// No. Try a Directshow video source instead:
					if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach dshowvideosrc as video source...\n");
					videosource = gst_element_factory_make("dshowvideosrc", "ptb_videosource");

					if (videosource) {
						// Specific deviceIndex requested, instead of auto-select?
						if (deviceIndex < 0) {
							if (deviceIndex != -5) {
								// Fetch optional targetmoviename parameter as name spec string:
								if (targetmoviefilename == NULL) PsychErrorExitMsg(PsychError_user, "You set 'deviceIndex' to a negative value, but didn't provide the required device name string in the 'moviename' argument! Aborted.");
								
								// Assign:
								strcat(config, targetmoviefilename);
							}
							
							switch(deviceIndex) {
								case -1:
								case -3:
									// Human friendly device name provided:
									if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device '%s' as video input.\n", config);
									g_object_set(G_OBJECT(videosource), "device-name", config, NULL);
									break;
									
								case -2:
								case -4:
									// DirectShow device path provided:
									if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device at path '%s' as video input.\n", config);
									g_object_set(G_OBJECT(videosource), "device", config, NULL);
									break;									
								case -5:
									if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach default DirectShow video device as video input.\n");
									break;
								default:
									PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' provided! Aborted.");
							}
						}
						else {
							// Human friendly device name provided:
							if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device '%s' as video input.\n", config);
							g_object_set(G_OBJECT(videosource), "device-name", config, NULL);							
						}
					}
				}
				// Have a standard windows video source ready.
			}
			
			// Still no video source available?
			if (!videosource) {
				if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed! We are out of options and will probably fail soon.\n");
				PsychErrorExitMsg(PsychError_system, "GStreamer failed to find a suitable video source! Game over.");
			}
		} // End of MS-Windows Video source creation.
		
		// MacOS/X specific setup path:
		if ((PSYCH_SYSTEM == PSYCH_OSX) && (deviceIndex > -8)) {
			if (deviceIndex < 0) {
				if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach osxvideosrc as video source...\n");
				videosource = gst_element_factory_make("osxvideosrc", "ptb_videosource");
				
				if (!videosource) {
					if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed! We are out of options and will probably fail soon.\n");
					PsychErrorExitMsg(PsychError_system, "GStreamer failed to find a suitable video source! Game over.");
				}
				
				// Fetch optional targetmoviename parameter as name spec string:
				if (deviceIndex != -5) {
					// Fetch optional targetmoviename parameter as name spec string:
					if (targetmoviefilename == NULL) PsychErrorExitMsg(PsychError_user, "You set 'deviceIndex' to a negative value, but didn't provide the required device name string in the 'moviename' argument! Aborted.");
					
					// Assign:
					strcat(config, targetmoviefilename);
				}
				
				switch(deviceIndex) {
					case -1:
					case -3:
						// Human friendly device name provided:
						if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device '%s' as video input.\n", config);
						g_object_set(G_OBJECT(videosource), "device-name", config, NULL);
						break;
						
					case -2:
					case -4:
						// System device path provided:
						if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach video device '%s' as video input.\n", config);
						g_object_set(G_OBJECT(videosource), "device", config, NULL);
						break;

					case -5:
					if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach default Quicktime OS/X video device as video input.\n");
						break;
						
					default:
						PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' provided! Not a value of -1 or -2. Aborted.");
				}
			}
		} // End of OS/X Video source creation.

        // OS independent special sources with no support for enumeration:
        
        // aravissrc for GeniCam cams: Doesn't support enumeration, just opening of cams by name:
        if (deviceIndex == -8) {

            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach aravissrc as GeniCam video source...\n");

            sprintf(plugin_name, "aravissrc");
            videosource = gst_element_factory_make("aravissrc", "ptb_videosource");
            
            if (!videosource) {
                if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed to create aravissrc GeniCam video source! We are out of options and will probably fail soon.\n");
                PsychErrorExitMsg(PsychError_system, "GStreamer failed to find a suitable video source! Game over.");
            }
            

            // Fetch optional targetmoviename parameter as name spec string:
            if (targetmoviefilename) {
                // Assign:
                strcat(config, targetmoviefilename);

                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach GeniCam video device '%s' as video input.\n", config);
                g_object_set(G_OBJECT(videosource), "camera-name", config, NULL);
            }
            else {
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to attach default Aravis GeniCam video device as video input.\n");
            }
        }

		// Some plugins need typefind'ing dhowvideosrc for sure, but we also set it for aravissrc to be safe:
		if (strstr(plugin_name, "dshowvideosrc") || strstr(plugin_name, "aravissrc")) g_object_set(G_OBJECT(videosource), "typefind", 1, NULL);
		
		// Enable timestamping by videosource, unless its been done already for a dv1394src:
		if (!strstr(plugin_name, "dv1394src")) g_object_set(G_OBJECT(videosource), "do-timestamp", 1, NULL);

		// Assign video source to pipeline:
        if (usecamerabin == 1) {
            // Attach directly to camerabin aka camerabin1:
            g_object_set(camera, "video-source", videosource, NULL);
        }
        else {
            // Attach indirectly to camerabin2 via a camerawrappersrc:
            videowrappersrc = gst_element_factory_make ("wrappercamerabinsrc", "ptbwrappervideosrc0");
            g_object_set(videowrappersrc, "video-source", videosource, NULL);
            g_object_set(camera, "camera-source", videowrappersrc, NULL);
        }
	}

    // Assign message context, message bus and message callback for
    // the pipeline to report events and state changes, errors etc.:    
    VideoContext = g_main_loop_new (NULL, FALSE);
    vidcapRecordBANK[slotid].VideoContext = VideoContext;
    bus = gst_pipeline_get_bus(GST_PIPELINE(camera));
    gst_bus_add_watch(bus, PsychVideoBusCallback, &(vidcapRecordBANK[slotid]));
    gst_object_unref(bus);

    // Name of target movie file for video and audio recording specified?
    if ((deviceIndex >= 0) && targetmoviefilename) {
	    // Assign it to camerabin to perform video recording:
	    if (!usecamerabin)
		    PsychErrorExitMsg(PsychError_user, "You requested video recording, but current fallback video engine doesn't support this. Aborted.");

	    // Codec settings or type specified?
	    if ((codecSpec = strstr(targetmoviefilename, ":CodecSettings="))) {
		    // Replace ':' with a zero in targetmoviefilename, so it gets null-terminated
		    // and only points to the actual movie filename:
		    *codecSpec = 0;

		    // Move after null-terminator:
		    codecSpec++;

		    // Replace the ':CodecSettings=' with the special keyword 'DEFAULTenc', so
		    // so the default video codec is chosen, but the given settings override its
		    // default parameters.
		    strncpy(codecSpec, "DEFAULTenc    ", strlen("DEFAULTenc    "));

		    if (strlen(codecSpec) == 0) PsychErrorExitMsg(PsychError_user, "Invalid (empty) :CodecSettings= parameter specified. Aborted.");
	    } else if ((codecSpec = strstr(targetmoviefilename, ":CodecType="))) {
		    // Replace ':' with a zero in targetmoviefilename, so it gets null-terminated
		    // and only points to the actual movie filename:
		    *codecSpec = 0;
		    // Advance codecSpec to point to the actual codec spec string:
		    codecSpec+= 11;

		    if (strlen(codecSpec) == 0) PsychErrorExitMsg(PsychError_user, "Invalid (empty) :CodecType= parameter specified. Aborted.");
	    } else {
		    // No codec specified: Use our default encoder, the one that's been shown to
		    // produce good results:
		    codecSpec = strdup("DEFAULTenc");
	    }

	    // Create matching video encoder for codecSpec:
	    if (PsychSetupRecordingPipeFromString(capdev, codecSpec, codecName, FALSE, TRUE, (recordingflags & 2) ? TRUE : FALSE)) {
		    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Video%s recording into file [%s] enabled for device %i. Codec is [%s].\n",
								  ((recordingflags & 2) ? " and audio" : ""), targetmoviefilename, deviceIndex, codecName);
		    if (strcmp(codecSpec, "DEFAULTenc") == 0) free(codecSpec);
		    codecSpec = NULL;
	    } else {
		    if (strcmp(codecSpec, "DEFAULTenc") == 0) free(codecSpec);
		    PsychErrorExitMsg(PsychError_system, "Setup of video recording failed. Reason hopefully given above.");
	    }

	    capdev->targetmoviefilename = strdup(targetmoviefilename); 
	    capdev->recording_active = TRUE;
    } else {
	    capdev->recording_active = FALSE;
    }

    // Create a appsink named "ptbsink0" for attachment as the destination video-sink for
    // all video content. This allows us to get hold of the videoframe buffers for
    // converting them into PTB OpenGL textures:
    videosink = gst_element_factory_make ("appsink", "ptbsink0");
    if (!videosink) {
        printf("PTB-ERROR: Failed to create video-sink appsink ptbsink!\n");
        PsychGSProcessVideoContext(vidcapRecordBANK[slotid].VideoContext, TRUE);
        PsychErrorExitMsg(PsychError_system, "Opening the videocapture device failed. Reason hopefully given above.");
    };

    vidcapRecordBANK[slotid].videosink = videosink;

    // Our OpenGL texture creation routine needs GL_BGRA8 data in GL_UNSIGNED_8_8_8_8_REV
    // format, but the pipeline usually delivers YUV data in planar format. Therefore
    // need to perform colorspace/colorformat conversion. These colorcaps, as assigned to
    // the videosink will try to nudge the video source to deliver data in our requested
    // format. If this ain't possible, they will enforce creation of a colorspace converter
    // inbetween the video source and our videosink:
    switch (reqdepth) {
        case 4:
            colorcaps = gst_caps_new_simple("video/x-raw-rgb",
                                            "bpp", G_TYPE_INT, capdev->pixeldepth,
                                            "depth", G_TYPE_INT, capdev->pixeldepth,
                                            "alpha_mask", G_TYPE_INT, 0x000000FF,
                                            "red_mask", G_TYPE_INT,   0x0000FF00,
                                            "green_mask", G_TYPE_INT, 0x00FF0000,
                                            "blue_mask", G_TYPE_INT,  0xFF000000,
                                            NULL);
        break;

        case 3:
            colorcaps = gst_caps_new_simple("video/x-raw-rgb",
                                            "bpp", G_TYPE_INT, capdev->pixeldepth,
                                            "depth", G_TYPE_INT, capdev->pixeldepth,
                                            "red_mask", G_TYPE_INT,   0x00FF0000,
                                            "green_mask", G_TYPE_INT, 0x0000FF00,
                                            "blue_mask", G_TYPE_INT,  0x000000FF,
                                            NULL);
        break;

        case 2:
        case 1:
            colorcaps = gst_caps_new_simple("video/x-raw-gray",
                                            "bpp", G_TYPE_INT, capdev->pixeldepth,
                                            NULL);
        break;

        case 5: // YUYV-422 packed pixel encoding:
            colorcaps = gst_caps_new_simple("video/x-raw-yuv",
                                            "format", GST_TYPE_FOURCC, GST_MAKE_FOURCC('U', 'Y', 'V', 'Y'),
                                            NULL);
            reqdepth = 2;
            capdev->reqpixeldepth = 2;
            capdev->pixeldepth = 16;
        break;

        case 6: // YUV-I420 planar pixel encoding:
            colorcaps = gst_caps_new_simple("video/x-raw-yuv",
                                            "format", GST_TYPE_FOURCC, GST_MAKE_FOURCC('I', '4', '2', '0'),
                                            NULL);
            reqdepth = 2;
            capdev->reqpixeldepth = 2;
            capdev->pixeldepth = 12;
        break;

        default:
        PsychErrorExitMsg(PsychError_internal, "Unknown reqdepth parameter received!");            
    }

    // Assign 'colorcaps' as caps to our videosink. This marks the videosink so
    // that it can only receive video image data in the format defined by colorcaps,
    // i.e., a format that is easy to consume for OpenGL's texture creation on std.
    // gpu's. It is the job of the video pipeline's autoplugger to plug in proper
    // color & format conversion plugins to satisfy videosink's needs.
    gst_app_sink_set_caps(GST_APP_SINK(videosink), colorcaps);

    // ROI rectangle specified and use of it supported by camerabin?
    if (usecamerabin && capturerectangle) {
        if ((capturerectangle[kPsychLeft] == 0) && (capturerectangle[kPsychTop] == 0)) {
            // roi = [0 0 w h] --> Specs a target capture resolution.
            // Extract wanted width and height and use it as target capture resolution:
            twidth  = (int) PsychGetWidthFromRect(capturerectangle);
            theight = (int) PsychGetHeightFromRect(capturerectangle);
        } else {
            // roi = [l t r b] --> Specs a ROI to crop out of full res capture.

            // Don't change videocapture resolution -- Leave it at auto-detected settings:
            twidth  = -1;
            theight = -1;

            // Create videocrop filter to crop away everything outside the defined ROI:
            videocrop_filter = gst_element_factory_make ("videocrop", "ptbvideocropfilter");
            if (!videocrop_filter) {
                // Disable capturerectangle, so we revert to full device default resolution:
                capturerectangle = NULL;
                if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Selection of specified video ROI not supported by setup. Using full resolution.\n");
            }
        }
    }
    else {
        if (capturerectangle) {
            // ROI spec'd but camerabin not supported. This is a no-go.
            capturerectangle = NULL;
            if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Selection of video ROI's or capture resolution not supported by your setup. Using full resolution.\n");
        }

        // No user specified resolution. Rely on auto-detection:
        twidth  = -1;
        theight = -1;
    }

    // Assign our special appsink 'videosink' as video-sink of the pipeline:
    if (!usecamerabin) {
        if (capdev->recording_active) PsychErrorExitMsg(PsychError_user, "Video recording requested, but this isn't supported on this setup, sorry.");
        g_object_set(G_OBJECT(camera), "video-sink", videosink, NULL);

        // Invalidate capdev->width,height/fps etc. so it gets auto-detected
        // by preroll op:
        capdev->width = capdev->height = 0;
        capdev->fps = 0;
    } else {
	    // Attach our appsink as videosink:
	    g_object_set(G_OBJECT(camera), "viewfinder-sink", videosink, NULL);

	    // Setup pipeline for video recording or pure capture:
	    if (!capdev->recording_active) {
		    // Pure video capture, no recording: Optimize pipeline for this case:
		    g_object_set(G_OBJECT(camera), "flags", (usecamerabin == 1) ? 1+2+4 : 0, NULL);
            
            // Setup colorcaps for camerabin1, do it later for camerabin2:
            if (usecamerabin == 1) g_object_set(G_OBJECT(camera), "filter-caps", colorcaps, NULL);
            
	    } else {
		    // Video recording (with optional capture). Setup pipeline:
            if (usecamerabin == 1) {
                g_object_set(G_OBJECT(camera),
                                // Only enable sound encoding if "audio recording" flag 2 is set in
                                // recordingflags. Otherwise add flags 0x20 to disable audio encoding:
                                "flags", 0x04 + 0x02 + ((recordingflags & 2) ? 0x10 : 0x20),
                                NULL);
            }
            else {
                // camerabin2:
                g_object_set(G_OBJECT(camera),
                             "flags", 0,
                             NULL);
                
                // Audio recording requested?
                if (recordingflags & 2) {
                    // Yes:
                    g_object_set(G_OBJECT(camera),
                                 "mute", FALSE,
                                 NULL);                    
                }
                else {
                    // No:
                    g_object_set(G_OBJECT(camera),
                                 "mute", TRUE,
                                 NULL);
                }
            }
	    }

	    // Create and use a videorate converter always when video recording is active (because without it
	    // we can get choppy or truncated video recordings and/or audio-video sync problems), and
	    // optionally in pure video feedback mode if explicitely requested by usercode via
	    // recordingflags & 128, as this is usually not what one wants for live video processing or
	    // feedback.
	    if (capdev->recording_active || (recordingflags & 128)) {
		    // Attach a videorate converter as source filter. It should take care of adjusting
		    // for fluctuating framerates from the videosource device by dropping or duplicating
		    // videoframes and adjusting timestamps as needed to generate a video stream with the
		    // selected target output framerate:
		    capdev->videorate_filter = gst_element_factory_make ("videorate", "ptbvideoratefilter");
		    if (capdev->videorate_filter == NULL) {
			    if (PsychPrefStateGet_Verbosity() > 1) {
				    printf("PTB-WARNING: Could not create video rate conversion filter. Therefore can't compensate\n");
				    printf("PTB-WARNING: for possible fluctuations in videocapture framerate. This can cause\n");
				    printf("PTB-WARNING: audio-video sync problems or wrong timing/speed in recorded videos!\n");
			    }
		    }
		    else {
                if (usecamerabin == 2) {
                    // Camerabin 2 path:
                    //
                    // By default we apply the videorate converter upstream at the videosource,
                    // so it affects both video recording and live video feedback. Usercode can
                    // restrict conversion to the recorded video stream only by setting
                    // recordingflags & 256:
                    if (recordingflags & 256) {
                        // Apply to video recording only -- attach to video-filter:
                        g_object_set(G_OBJECT(camera), "video-filter", capdev->videorate_filter, NULL);
                        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Restricting video framerate conversion to recorded video.\n");                        
                    }
                    else {
                        // Attach as video-source-filter -- upstream right at the source:
                        videosource_filter = capdev->videorate_filter;
                        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Video framerate conversion applies to recorded video and live video.\n");
                    }
                }
                else {
                    // Camerabin 1 path:
                    //
                    // By default we apply the videorate converter upstream at the videosource,
                    // so it affects both video recording and live video feedback. Usercode can
                    // restrict conversion to the recorded video stream only by setting
                    // recordingflags & 256:
                    if (recordingflags & 256) {
                        // Attach to video-post-processing -- immediately before video encoder:
                        g_object_set(G_OBJECT(camera), "video-post-processing", capdev->videorate_filter, NULL);
                        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Restricting video framerate conversion to recorded video.\n"); 
                    }
                    else {
                        // Attach as video-source-filter -- upstream right at the source:
                        videosource_filter = capdev->videorate_filter;
                        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Video framerate conversion applies to recorded video and live video.\n"); 
                    }
                }
		    }
	    }

        // If we are using camerabin2, we need to READY the pipeline, otherwise the queries
        // and validations below will fail due to lack of available video capture caps:
        if (usecamerabin == 2) {
            // Ready the pipeline:
            if (!PsychVideoPipelineSetState(camera, GST_STATE_READY, 30.0)) {
                PsychGSProcessVideoContext(vidcapRecordBANK[slotid].VideoContext, TRUE);
                PsychErrorExitMsg(PsychError_user, "In OpenVideoCapture: Opening the video capture device failed during camerabin2 pipeline zero -> ready. Reason given above.");
            }
        }

	    // Specific capture resolution requested?
	    if ((twidth != -1) && (theight != -1)) {
            // Yes. Validate and request it.
            
            // The usual crap. Enumeration of supported resolutions doesn't work with various video sources, so
            // we skip validation and trust blindly that the usercode is right if this is one of the non-enumerating
            // video sources:
            if (!strstr(plugin_name, "dshowvideosrc") && !strstr(plugin_name, "autovideosrc") && !strstr(plugin_name, "videotestsrc") &&
                !strstr(plugin_name, "aravissrc")) {
                // Query camera if it supports the requested resolution:
                capdev->fps = -1;
                if (!PsychGSGetResolutionAndFPSForSpec(capdev, &twidth, &theight, &capdev->fps, reqdepth)) {
                    // Unsupported resolution. Game over!
                    printf("PTB-ERROR: Video capture device %i doesn't support requested video capture resolution of %i x %i pixels! Aborted.\n", slotid, twidth, theight);
                    PsychErrorExitMsg(PsychError_user, "Failed to open video device at requested video resolution.");
                }
            }
            
            // Resolution supported. Request it:
            if (usecamerabin == 1) {
                g_object_set(G_OBJECT(camera),
                             "video-capture-width", twidth,
                             "video-capture-height", theight,
                             NULL);
            }
            else {
                // Setup colorcaps for camerabin2:
                gst_caps_set_simple(colorcaps, "width", G_TYPE_INT, twidth, "height", G_TYPE_INT, theight, NULL);
                g_object_set(G_OBJECT(camera), "viewfinder-caps", colorcaps, NULL);
                if (capdev->recording_active) {
                    g_object_set(G_OBJECT(camera), "video-capture-caps", colorcaps, NULL);
                }
            }
            
            // Assign requested and validated resolution as capture resolution of video source:
            capdev->width = twidth;
            capdev->height = theight;
            capdev->fps = 0;
	    }
	    else {
            // No, usercode wants auto-detected default resolution. Query
            // highest possible resolution (maximum pixel area) and choose that:
            capdev->width = -1;
            capdev->height = -1;
            capdev->fps = -1;

            // Auto-Detection doesn't work with various video source plugins. Skip it and hope that later probing code will do the job:
            if (!strstr(plugin_name, "dshowvideosrc") && !strstr(plugin_name, "autovideosrc") && !strstr(plugin_name, "videotestsrc") &&
                !strstr(plugin_name, "aravissrc")) {
                // Ask camera to provide auto-detected parameters:
                if (!PsychGSGetResolutionAndFPSForSpec(capdev, &capdev->width, &capdev->height, &capdev->fps, reqdepth)) {
                    // Unsupported resolution. Game over!
                    printf("PTB-ERROR: Auto-Detection of optimal video resolution on video capture device %i failed! Aborted.\n", slotid);
                    PsychErrorExitMsg(PsychError_user, "Failed to open video device with auto-detected video resolution.");
                }
                
                // Resolution supported. Request it:
                if (usecamerabin == 1) {
                    g_object_set(G_OBJECT(camera),
                                 "video-capture-width", capdev->width,
                                 "video-capture-height", capdev->height,
                                 NULL);
                }
                else {
                    // Setup colorcaps for camerabin2:
                    gst_caps_set_simple(colorcaps, "width", G_TYPE_INT, capdev->width, "height", G_TYPE_INT, capdev->height, NULL);
                    g_object_set(G_OBJECT(camera), "viewfinder-caps", colorcaps, NULL);
                    if (capdev->recording_active) {
                        g_object_set(G_OBJECT(camera), "video-capture-caps", colorcaps, NULL);
                    }
                }
            }
            else {
                // Source without enumeration capability. Set "don't know" values and hope the fallback code below
                // does a better job at guessing the true source resolution:
                capdev->width  = 0;
                capdev->height = 0;

                // ROI defined? We can't handle this, because we don't know the true video capture resolution
                // which would be needed at this point in the setup path due to the broken enumeration.
                if (capturerectangle) {
                    capturerectangle = NULL;
                    overrideFrameSize = TRUE;

                    // Delete useless videocrop element if any:
                    if (videocrop_filter) gst_object_unref(G_OBJECT(videocrop_filter));
                    videocrop_filter = NULL;

                    if (PsychPrefStateGet_Verbosity() > 1) {
                        printf("PTB-WARNING: Usercode specified a 'roirectangle' with a ROI for video cropping, but this system setup\n");
                        printf("PTB-WARNING: doesn't support this. Ignoring 'roirectangle' and reverting to full video capture resolution.\n");
                    }
                }
            }
            
            // Reset capdev->fps to neutral zero: capdev->width and capdev->height are already auto-assigned.
            capdev->fps = 0;
	    }
        
        // We can only assign videosource_filters (for videorate plugin) or videocrop_filter
        // to camerabin2 while the pipeline is in NULL state. Therefore make it so:
        if ((usecamerabin == 2) && (videosource_filter || videocrop_filter)) {
            if (!PsychVideoPipelineSetState(camera, GST_STATE_NULL, 30.0)) {
                PsychGSProcessVideoContext(vidcapRecordBANK[slotid].VideoContext, TRUE);
                PsychErrorExitMsg(PsychError_user, "In OpenVideoCapture: Opening the video capture device failed during intermediate camerabin2 pipeline ready -> zero transition. Reason given above.");
            }
        }

	    // videocrop_filter for ROI processing available?
	    if (videocrop_filter) {
		    // Setup final cropping region:
		    g_object_set(G_OBJECT(videocrop_filter),
				 "left",   (int) capturerectangle[kPsychLeft],
				 "top",    (int) capturerectangle[kPsychTop],
				 "right",  capdev->width - (int) capturerectangle[kPsychRight],
				 "bottom", capdev->height - (int) capturerectangle[kPsychBottom],
				 NULL);

            if (usecamerabin == 1) {
                // Camerabin1 path:

                // Is video recording active and ROI cropping for recording explicitely enabled
                // via recordingflags setting 512?
                if (capdev->recording_active && (recordingflags & 512)) {
                    // Yes. Attach filter to video-post-processing. This way it affects
                    // recorded video, but unfortunately also *always* the viewfinder video feed.
                    g_object_set(G_OBJECT(camera), "video-post-processing", videocrop_filter, NULL);
                    overrideFrameSize = FALSE;

                    if (PsychPrefStateGet_Verbosity() > 1) {
                        printf("PTB-WARNING: Application of ROI's to recorded video is somewhat limited:\n");
                        printf("PTB-WARNING: The ROI is correctly applied to the recorded video, but it\n");
                        printf("PTB-WARNING: also is applied to the live video feed, regardless if you\n");
                        printf("PTB-WARNING: want this or not. This is a limitation of camerabin1 capture.\n");
                    }
                }
                else {
                    // No. Either no videorecording active -- pure live feedback, or it is
                    // active, but cropping of recorded video is not enabled.
                    // This means only the viewfinder, aka appsink aka our live video textures
                    // shall be cropped to ROI, if at all. This works as expected and is
                    // thankfully the common case for most applications of ROI's.

                    // Check if live feed is enabled and ROI cropping for it is not disabled:
                    if (!(recordingflags & 4) && !(recordingflags & 1024)) {
                        // Yes. Viewfinder feed shall be cropped: Attach to viewfinder-filter:
                        g_object_set(G_OBJECT(camera), "viewfinder-filter", videocrop_filter, NULL);

                        // No stupid hack needed in this case:
                        overrideFrameSize = FALSE;
                    }
                    else {
                        // Disable application of ROI to video buffers:
                        overrideFrameSize = TRUE;
                    }
                }
            }

            if (usecamerabin == 2) {
                // Camerabin2 path:
                
                // Check if live feed is enabled and ROI cropping for it is not disabled:
                if (!(recordingflags & 1024) && !(recordingflags & 4)) {
                    g_object_set(G_OBJECT(camera), "viewfinder-filter", videocrop_filter, NULL);

                    // We need to adapt the caps of the viewfinder to expect image data of the
                    // expected size after cropping takes place:
                    vfcaps = gst_caps_copy(colorcaps);
                    gst_caps_set_simple(vfcaps, "width", G_TYPE_INT, (int) PsychGetWidthFromRect(capturerectangle),
                                                "height", G_TYPE_INT, (int) PsychGetHeightFromRect(capturerectangle), NULL);
                    g_object_set(G_OBJECT(camera), "viewfinder-caps", vfcaps, NULL);
                    gst_caps_unref(vfcaps);

                    // No stupid hack needed in this case:
                    overrideFrameSize = FALSE;
                }
                else {
                    // Disable application of ROI to video buffers:
                    overrideFrameSize = TRUE;
                }

                // Is videorecording active and cropping for it requested?
                if ((recordingflags & 512) && capdev->recording_active) {
                    // Need to attach to "video-filter", which affects video encoding and recording.
                    // We need to setup a 2nd cropping filter with identical settings for this, unless
                    // the 1st one isn't used as the viewfinder-filter.
                    if ((recordingflags & 1024) || (recordingflags & 4)) {
                        // Live feed disabled or no cropping for live feed wanted. The videocrop_filter
                        // is unused and we can use it here for cropping the videorecording by
                        // attaching to video-filter:
                        g_object_set(G_OBJECT(camera), "video-filter", videocrop_filter, NULL);
                    }
                    else {
                        // videocrop_filter already used for cropping the live feed aka viewfinder.
                        // Generate a new videocrop element, set it up identically, attach it to
                        // video-filter:
                        videocrop_filter = gst_element_factory_make ("videocrop", "ptbvideoreccropfilter");
                        if (!videocrop_filter) {
                            if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Failed to apply video ROI to videorecording! Using full resolution.\n");
                        }
                        else {
                            g_object_set(G_OBJECT(videocrop_filter),
                                 "left",   (int) capturerectangle[kPsychLeft],
                                 "top",    (int) capturerectangle[kPsychTop],
                                 "right",  capdev->width - (int) capturerectangle[kPsychRight],
                                 "bottom", capdev->height - (int) capturerectangle[kPsychBottom],
                                 NULL);
                            g_object_set(G_OBJECT(camera), "video-filter", videocrop_filter, NULL);
                        }
                    }
                    
                    // Cropping applies to video frames which go to the encoder. Adapt video encoder bin caps
                    // to expect image data of the size after cropping takes place:
                    reccaps = gst_caps_copy(colorcaps);
                    gst_caps_set_simple(reccaps, "width", G_TYPE_INT, (int) PsychGetWidthFromRect(capturerectangle),
                                                 "height", G_TYPE_INT, (int) PsychGetHeightFromRect(capturerectangle), NULL);
                    g_object_set(G_OBJECT(camera), "video-capture-caps", reccaps, NULL);
                    gst_caps_unref(reccaps);          
                }
            }

            // End of video cropping setup for ROI implementation.
	    }

	    // Attach videosource filter upstream at video source, if any. This way it affects
        // both video recording and video live feed:
	    if (videosource_filter) {
            if (usecamerabin == 1) {
                g_object_set(G_OBJECT(camera), "video-source-filter", videosource_filter, NULL);
            }
            else {
                // Camerabin2: Can't attach to camerabin2 directly, but to wrappercamerasrc:
                g_object_set(G_OBJECT(videowrappersrc), "video-source-filter", videosource_filter, NULL);                
            }
	    }
    }

    gst_caps_unref(colorcaps);

    // Disable asynchronous transition to GST_PAUSED state for video sink when simultaneous
    // video recording is active: TODO FIXME: This seems to have no effect whatsoever....
    if (capdev->recording_active) g_object_set(G_OBJECT(videosink), "async", FALSE, NULL);

    // Disable internal queuing of the last buffer, as we don't ever use the "last-buffer" property:
    g_object_set(G_OBJECT(videosink), "enable-last-buffer", FALSE, NULL);

    // Get the pad from the final sink for probing width x height of video frames and nominal framerate of video source:	
    pad = gst_element_get_pad(videosink, "sink");

    // Install callbacks used by the videosink (appsink) to announce various events:
    gst_app_sink_set_callbacks(GST_APP_SINK(videosink), &videosinkCallbacks, &(vidcapRecordBANK[slotid]), PsychDestroyNotifyCallback);

    // Specific maximum number of buffers to use?
    if (capdev->num_dmabuffers) {
	    // Yes: Drop frames if callback can't pull buffers fast enough:
	    gst_app_sink_set_drop(GST_APP_SINK(videosink), TRUE);
	    
	    // Allownum_dmabuffers to be queued before dropping buffers:
	    gst_app_sink_set_max_buffers(GST_APP_SINK(videosink), capdev->num_dmabuffers);
    }
    else {
	    // No: Leave queue length at default and disable buffer dropping:
	    gst_app_sink_set_drop(GST_APP_SINK(videosink), FALSE);
    }

    PsychGSProcessVideoContext(vidcapRecordBANK[slotid].VideoContext, FALSE);

    if ((capdev->recording_active) && (recordingflags & 4)) {
	    // No live feedback (flags 4) implies no preroll (flags 8):
	    // This because a preroll would bring the pipeline from null state to
	    // at least ready state and thereby prevent us from detaching our
	    // appsink and replacing it with a fakesink. This would cause our
	    // appsink to collect video buffers which wouldn't get pulled/drained,
	    // effectively filling up all memory with "leaked" buffers until
	    // out of memory abort happens after a minute or so.
	    //
	    // Not prerolling -> pipeline in null state at fakesink attach time
	    // -> problem solved.
	    recordingflags |= 8;

	    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Pipeline preroll skipped because only harddisc recording (recordingflags & 4).\n");
    }

    if (strstr(plugin_name, "dv1394src")) {
	    recordingflags |= 8;        
	    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Pipeline preroll skipped because DV capture via dv1394src or hdv1394src is active.\n");
    }

    // Only preroll if prerolling not disabled by recordingflag 8,
    // or if we use the fallback path, which utterly needs this:
    if (!(recordingflags & 8) || !usecamerabin) {
	    // Ready the pipeline:
	    if (!PsychVideoPipelineSetState(camera, GST_STATE_READY, 30.0)) {
		    PsychGSProcessVideoContext(vidcapRecordBANK[slotid].VideoContext, TRUE);
		    PsychErrorExitMsg(PsychError_user, "In OpenVideoCapture: Opening the video capture device failed during pipeline zero -> ready. Reason given above.");
	    }

	    // Preroll the pipeline:
	    if (!PsychVideoPipelineSetState(camera, GST_STATE_PLAYING, 30.0)) {
		    PsychGSProcessVideoContext(vidcapRecordBANK[slotid].VideoContext, TRUE);
		    PsychErrorExitMsg(PsychError_user, "In OpenVideoCapture: Opening the video capture device failed during pipeline preroll ready->playing. Reason given above.");
	    }
    }

    if (!usecamerabin) {
	    g_object_get (G_OBJECT(camera),
			  "source", &videosource,
			  NULL);
    } else {
	    g_object_get (G_OBJECT(camera),
              (usecamerabin == 1) ? "video-source" : "camera-source", &videosource,
			  NULL);
    }

    // Query number of available video and audio channels on capture device:
    if (!usecamerabin) {
	    g_object_get(G_OBJECT(camera),
		         "n-video", &vidcapRecordBANK[slotid].nrVideoTracks,
		         "n-audio", &vidcapRecordBANK[slotid].nrAudioTracks,
			 NULL);
    } else {
        // There's always at least a video channel:
        vidcapRecordBANK[slotid].nrVideoTracks = 1;
    }

    // We need a valid onscreen window handle for real video playback:
    if ((NULL == win) && (vidcapRecordBANK[slotid].nrVideoTracks > 0) && !(recordingflags & 4)) {
	    PsychErrorExitMsg(PsychError_user, "No windowPtr to an onscreen window provided. Must do so for sources with video channels!");
    }
 
    PsychGSProcessVideoContext(vidcapRecordBANK[slotid].VideoContext, FALSE);

    // Assign harmless initial settings for fps and frame size:
    rate1 = 0;
    rate2 = 1;
    width = height = 0;

    // Videotrack available and return of video data enabled?
    if (vidcapRecordBANK[slotid].nrVideoTracks > 0) {
        // Yes: Query video frame size and framerate of device:
        peerpad = gst_pad_get_peer(pad);
        caps = NULL;
        if (GST_IS_PAD(peerpad)) caps = gst_pad_get_negotiated_caps(peerpad);
        if (caps) {
            str=gst_caps_get_structure(caps,0);

            /* Get some data about the frame */
            gst_structure_get_int(str,"width",&width);
            gst_structure_get_int(str,"height",&height);
            gst_structure_get_fraction(str, "framerate", &rate1, &rate2);
            gst_caps_unref(caps);
            if (PsychPrefStateGet_Verbosity() > 4) printf("Negotiated videosink res: %i x %i\n", width, height);
        } else {
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: No frame info available after preroll.\n");
        }
    }

    // Release the pad:
    gst_object_unref(pad);

    // Our camera should be ready: Assign final handle.
    *capturehandle = slotid;

    // Increase counter of open capture devices:
    numCaptureRecords++;

    // Assign our first guess on video device resolution, only already
    // setup by detection code above:
    if ((capdev->width == 0) && (capdev->height == 0)) {
        capdev->width = width;
        capdev->height = height;
    }
    
    // Sanity check: grabber_active should be false, because we initialized it so.
    // It will be true if the videosink went EOS during the open and preroll sequence
    // due to something going wrong (see PsychEOSCallback()). We try to catch this here
    // and abort the open sequence instead of later hanging the whole app...
    PsychLockMutex(&capdev->mutex);
    if (capdev->grabber_active) {
        capdev->grabber_active = false;
        PsychUnlockMutex(&capdev->mutex);
        PsychErrorExitMsg(PsychError_user, "In OpenVideoCapture: Opening the video capture device failed during preroll due to premature EOS.");
    }
    PsychUnlockMutex(&capdev->mutex);

    // Reset framecounter:
    capdev->nrframes = 0;
    capdev->grabber_active = 0;

    rate1 = 0;
    rate2 = 1;
    width = height = 0;

    // Query true properties of attached video source:
    videosource = NULL;
	
    if (!usecamerabin) {
	    g_object_get (G_OBJECT(camera),
			  "source", &videosource,
			  NULL);
	    capdev->videoenc = NULL;
    } else {
        if (usecamerabin == 1) {
            g_object_get (G_OBJECT(camera), "video-source", &videosource, NULL);
        }
        else {
            g_object_get (G_OBJECT(videowrappersrc), "video-source", &videosource, NULL);
        }
        
        if (usecamerabin == 1) {
            g_object_get (G_OBJECT(camera),
                  "video-encoder", &capdev->videoenc,
                  NULL);
        }
        else {
            // camerabin2: Doesn't use explicit video encoder, but video encoding profiles...
            capdev->videoenc = NULL;
        }
    }

    if (videosource) {
        pstring = NULL; 
        if (g_object_class_find_property(G_OBJECT_GET_CLASS(videosource), "device")) {
            g_object_get(G_OBJECT(videosource), "device", &pstring, NULL);
            if (pstring) {
                if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Camera device name is '%s'.\n", pstring);
                g_free(pstring); pstring = NULL;
            }
        }

        if (g_object_class_find_property(G_OBJECT_GET_CLASS(videosource), "device-path")) {
            g_object_get(G_OBJECT(videosource), "device-path", &pstring, NULL);
            if (pstring) {
                if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Camera device-path is '%s'.\n", pstring);
                g_free(pstring); pstring = NULL;
            }
        }

        if (g_object_class_find_property(G_OBJECT_GET_CLASS(videosource), "device-name")) {
            g_object_get(G_OBJECT(videosource), "device-name", &pstring, NULL);
            if (pstring) {
                if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Camera friendly device-name is '%s'.\n", pstring);
                capdev->cameraFriendlyName = strdup(pstring);
                g_free(pstring); pstring = NULL;
            }
            else {
                capdev->cameraFriendlyName = strdup("Unknown");
            }
        }

        if (g_object_class_find_property(G_OBJECT_GET_CLASS(videosource), "camera-name")) {
            g_object_get(G_OBJECT(videosource), "camera-name", &pstring, NULL);
            if (pstring) {
                if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Camera friendly device-name is '%s'.\n", pstring);
                capdev->cameraFriendlyName = strdup(pstring);
                g_free(pstring); pstring = NULL;
            }
            else {
                capdev->cameraFriendlyName = strdup("Unknown");
            }
        }
            
        // Get the pad from the src pad of the source for probing width x height
        // of video frames and nominal framerate of video source:	
        pad = gst_element_get_pad(videosource, "src");

        // Videotrack available?
        if (vidcapRecordBANK[slotid].nrVideoTracks > 0) {
            // Yes: Query video frame size and framerate of device:
            peerpad = gst_pad_get_peer(pad);
            caps = NULL;
            if (GST_IS_PAD(peerpad)) caps = gst_pad_get_negotiated_caps(peerpad);
            if (caps) {
                str=gst_caps_get_structure(caps,0);

                /* Get some data about the frame */
                gst_structure_get_int(str,"width",&width);
                gst_structure_get_int(str,"height",&height);
                gst_structure_get_fraction(str, "framerate", &rate1, &rate2);
                if (PsychPrefStateGet_Verbosity() > 4) printf("Negotiated videosource w = %i h = %i fps = %f\n", width, height, rate1/rate2);
                gst_caps_unref(caps);
            } else {
                 if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: No frame info for video source available after preroll.\n");	
            }
        }

        // Release the pad:
        gst_object_unref(pad);		
    }

    // If correct width x height not yet assigned, retry assignment from
    // results of videosource query:
    if (capdev->width == 0) capdev->width = width;
    if (capdev->height == 0) capdev->height = height;

    // If we totally failed to query width and height, set what usercode
    // requested - if it requested anything:
    if (capdev->width == 0) capdev->width = (twidth > 0) ? twidth : 0;
    if (capdev->height == 0) capdev->height = (theight > 0) ? theight : 0;

    // Create final user-visible effective ROI:
    if (capturerectangle && ((twidth == -1) && (theight == -1))) {
	    // Usercode specified a ROI, videosource is left at its full default
	    // resolution, ROI is applied via videocrop_filter. Specified ROI is
	    // effective ROI, so just copy it:
	    PsychCopyRect(capdev->roirect, capturerectangle);
    } else {
	    // Usercode specified no ROI, but a target device resolution or no
	    // resolution at all. Videosource is configured for a certain resolution,
	    // as stored in capdev->width x height. Our effective ROI is simply
	    // the full videosource resolution:
	    PsychMakeRect(capdev->roirect, 0, 0, capdev->width, capdev->height);
    }

    // The size of returned video frames for texture creation or image buffer return
    // etc., is equal to the width x height of the ROI:
    capdev->frame_width  = (int) PsychGetWidthFromRect(capdev->roirect);
    capdev->frame_height = (int) PsychGetHeightFromRect(capdev->roirect);

    // Do we need to force the frame_width x frame_height to full capture
    // resolution?
    if (overrideFrameSize) {
	    // Yes.
	    capdev->frame_width  = capdev->width;
	    capdev->frame_height = capdev->height;
    }

    // If we prerolled before, we need to undo its effects:
    if (!(recordingflags & 8) || !usecamerabin) {
	    // Pause the pipeline again:
	    if (!PsychVideoPipelineSetState(camera, (usecamerabin) ? GST_STATE_READY : GST_STATE_PAUSED, 30.0)) {
		    PsychGSProcessVideoContext(vidcapRecordBANK[slotid].VideoContext, TRUE);
		    PsychErrorExitMsg(PsychError_user, "In OpenVideoCapture: Opening the video capture device failed during preroll playing -> pause. Reason given above.");
	    }
    }

    if ((capdev->recording_active) && (recordingflags & 4)) {
	    // Only disc recording, no data return: Attach a fakesink to effectively
	    // disable most of viewfinder processing:
	    videosink = gst_element_factory_make ("fakesink", "ptbfakesink0");
	    g_object_set(G_OBJECT(camera), "viewfinder-sink", videosink, NULL);
    }

    // Obsolete on Linux, hopefully also on Windows, but needed on MacOS/X:
    // Make sure our capture timestamps are in system time (GetSecs() reference time)
    // unless recordingflags 64 are set, in which case timestamps are in pipeline
    // running time, or whatever GStreamer finds convenient.
    if ((PSYCH_SYSTEM == PSYCH_OSX) && !(recordingflags & 64)) {
	    // Set defined start time for pipeline - Used for cts timestamps:
	    gst_element_set_base_time(camera, GST_CLOCK_TIME_NONE);
	    gst_element_set_start_time(camera, GST_CLOCK_TIME_NONE);
    }

    // Assign final recordingflags:
    vidcapRecordBANK[slotid].recordingflags = recordingflags;

    // Assign plugin name:
    sprintf(&capdev->videosourcename[0], "%s", plugin_name);
    
    // Store a pointer to the videosource plugin:
    capdev->videosource = videosource;
    capdev->videowrappersrc = videowrappersrc;
    
    if (PsychPrefStateGet_Verbosity() > 2) {
	    printf("PTB-INFO: Camera %i opened [Source resolution width x height = %i x %i, video image size %i x %i]\n",
		   slotid, capdev->width, capdev->height, capdev->frame_width, capdev->frame_height);
    }

    // Sanity check: grabber_active should be false, because we initialized it so.
    // It will be true if the videosink went EOS during the open and preroll sequence
    // due to something going wrong (see PsychEOSCallback()). We try to catch this here
    // and abort the open sequence instead of later hanging the whole app...
    PsychLockMutex(&capdev->mutex);
    if (capdev->grabber_active) {
        capdev->grabber_active = false;
        PsychUnlockMutex(&capdev->mutex);
        PsychErrorExitMsg(PsychError_user, "In OpenVideoCapture: Opening the video capture device failed during preroll due to premature EOS.");
    }
    PsychUnlockMutex(&capdev->mutex);

    return(TRUE);
}

/* Internal helper:
 * Drain up to 'numFramesToDrain' videobuffers from the videosink of 'capdev'.
 * 'flags' modify drain behaviour. Unused so far.
 * Return number of drained buffers.
 */
int PsychGSDrainBufferQueue(PsychVidcapRecordType* capdev, int numFramesToDrain, unsigned int flags)
{
	GstBuffer *videoBuffer = NULL;
	int drainedCount = 0;

	// Drain while anything available, but at most numFramesToDrain frames.
	while (GST_IS_APP_SINK(capdev->videosink) && !gst_app_sink_is_eos(GST_APP_SINK(capdev->videosink))
		&& (capdev->frameAvail > 0) && (numFramesToDrain > drainedCount)) {
		capdev->frameAvail--;
		videoBuffer = gst_app_sink_pull_buffer(GST_APP_SINK(capdev->videosink));
		gst_buffer_unref(videoBuffer);
		videoBuffer = NULL;
		drainedCount++;
	}

	return(drainedCount);
}

/* 
*  PsychGSVideoCaptureRate() - Start- and stop video capture.
*
*  capturehandle = Grabber to start-/stop.
*  capturerate = zero == Stop capture, non-zero == Capture at given rate if possible.
*  dropframes = 0 - Always deliver oldest frame in DMA ringbuffer. 1 - Always deliver newest frame.
*               --> 1 == drop frames in ringbuffer if behind -- low-latency capture.
*
*               If dropframes is 1 at 'stop' time, then discard all pending buffers in queue,
*               otherwise retain them so they can be fetched after stop.
*
*  startattime = Deadline (in system time) for which to wait before real start of capture.
*  Returns Number of dropped frames during capture.
*/
int PsychGSVideoCaptureRate(int capturehandle, double capturerate, int dropframes, double* startattime)
{
	GstElement              *camera = NULL;
	GstBuffer               *videoBuffer = NULL;
	guint64                 nrInFrames, nrOutFrames, nrDroppedFrames, nrDuplicatedFrames;
    GstCaps                 *caps = NULL;
    GstCaps                 *capsi = NULL;
    GstCaps                 *capss = NULL;
    GstCaps                 *capsr = NULL;
    psych_bool              fps_matched = FALSE;
    int                     idx, idx2, fps_n, fps_d;
    double                  fpsmin, fpsmax;
	int                     dropped = 0;
	int                     drainedCount;
	float                   framerate = 0;

	// Retrieve device record for handle:
	PsychVidcapRecordType* capdev = PsychGetGSVidcapRecord(capturehandle);
	camera = capdev->camera;

	// Make sure GStreamer is ready:
	PsychGSCheckInit("videocapture");

	PsychGSProcessVideoContext(capdev->VideoContext, FALSE);

	// Start- or stop capture?
	if (capturerate > 0) {
		// Start capture:
		if (capdev->grabber_active) PsychErrorExitMsg(PsychError_user, "You tried to start video capture, but capture is already started!");
		
		// Reset statistics:
		capdev->last_pts = -1.0;
		capdev->nr_droppedframes = 0;
		capdev->lastSavedBaseTime = 0;

		// Framedropping in the sense we define it is not supported by libGStreamer, so we implement it ourselves.
		// Store the 'dropframes' flag in our capdev struct, so the PsychGSGetTextureFromCapture()
		// knows how to handle this:
		capdev->dropframes = (dropframes > 0) ? 1 : 0;

		if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: StartVideoCapture: Draining any queued stale video buffers before start of capture.\n");
		drainedCount = PsychGSDrainBufferQueue(capdev, INT_MAX, 0);
		if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: StartVideoCapture: Drain of %i buffers completed.\n", drainedCount);

		// Only do state transitions for camerabin:
		if (usecamerabin) {
			// Start the video capture for this camera.
			if (!PsychVideoPipelineSetState(camera, GST_STATE_READY, 10.0)) {
				// Failed!
				PsychGSProcessVideoContext(capdev->VideoContext, FALSE);
				PsychErrorExitMsg(PsychError_user, "Failure in pipeline transition null -> ready - Start of video capture failed!");
			}

            // Camerabin2 in use and requested capturerate is not DBL_MAX aka "auto-selected default"?
            if ((usecamerabin == 2) && (capturerate < DBL_MAX)) {
                // Camerabin 2 with specific capturerate requested.
                // Modify caps for camerabin2 - Add a "framerate" property matching our requirements, but first
                // validate the requested capturerate against supported capture rates of this device at its
                // current settings for resolution, color format and color depths.

                // Get current capture caps (so far without a framerate) - Store a writable copy of them in caps:
                g_object_get(G_OBJECT(camera), "viewfinder-caps", &capsi, NULL);
                caps = gst_caps_copy(capsi);
                capsr = gst_caps_copy(capsi);
                gst_caps_unref(capsi);
                
                // Get list of supported capture caps for this device in capss:
                g_object_get(G_OBJECT(camera), "viewfinder-supported-caps", &capss, NULL);
                
                // Intersect with current caps to find the subset of possible caps, given the
                // already set/fixed resolution and color format/depth. This essentially leaves
                // us with caps representing the available framerates in capsi:
                gst_caps_set_simple(capsr, "width", G_TYPE_INT, capdev->width, "height", G_TYPE_INT, capdev->height, NULL);
                capsi = gst_caps_intersect(capsr, capss);
                gst_caps_unref(capsr);
                gst_caps_unref(capss);
                
                // Print 'em:
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Videosource intersected caps are: %" GST_PTR_FORMAT "\n\n", capsi);
                
                // Match requested framerate against supported framerates:
                for (idx = 0; idx < (int) gst_caps_get_size(capsi); idx++) {
                    // For idx'th cap, get framerate value with all framerates of that cap:
                    GstStructure *capsstruct = gst_caps_get_structure (capsi, idx);
                    const GValue* framerates = gst_structure_get_value(capsstruct, "framerate");
                    
                    // framerates can be in the format of a single fraction, a list of fractions, or
                    // an allowable range of fractions:
                    if (GST_VALUE_HOLDS_FRACTION(framerates)) {
                        if (gst_structure_get_fraction (capsstruct, "framerate", &fps_n, &fps_d)) {
                            if (fabs((int)(capturerate + 0.5) - ((double) fps_n / (double) fps_d)) < 1.0) fps_matched = TRUE;
                            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Caps %i : FPS %f Hz.\n", idx, (float) fps_n / (float) fps_d);
                        }
                    }
                    else if (GST_VALUE_HOLDS_LIST(framerates)) {
                        for (idx2 = 0; idx2 < (int) gst_value_list_get_size(framerates); idx2++) {
                            const GValue* value = gst_value_list_get_value (framerates, idx2);
                            fps_n = gst_value_get_fraction_numerator(value);
                            fps_d = gst_value_get_fraction_denominator(value);
                            if (fabs((int)(capturerate + 0.5) - ((double) fps_n / (double) fps_d)) < 1.0) fps_matched = TRUE;
                            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: %i-%i : FPS %f Hz.\n", idx, idx2, (float) fps_n / (float) fps_d);
                        }
                    }
                    else if (GST_VALUE_HOLDS_FRACTION_RANGE(framerates)) {
                        const GValue* frmin = gst_value_get_fraction_range_min(framerates);
                        const GValue* frmax = gst_value_get_fraction_range_max(framerates);
                        fps_n = gst_value_get_fraction_numerator(frmin);
                        fps_d = gst_value_get_fraction_denominator(frmin);
                        fpsmin = (double) fps_n / (double) fps_d;
                        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: %i: FPS min %f - ", idx, fpsmin);

                        fps_n = gst_value_get_fraction_numerator(frmax);
                        fps_d = gst_value_get_fraction_denominator(frmax);
                        fpsmax = (double) fps_n / (double) fps_d;
                        if (PsychPrefStateGet_Verbosity() > 4) printf("max %f\n", fpsmax);

                        if (((int)(capturerate + 0.5) >= fpsmin) && ((int)(capturerate + 0.5) <= fpsmax)) fps_matched = TRUE;
                    }
                }
                gst_caps_unref(capsi);
                
                // Requested capturerate supported?
                if (fps_matched) {
                    // Yes! Add it to our caps object:
                    gst_caps_set_simple(caps, "framerate", GST_TYPE_FRACTION, (int)(capturerate + 0.5), 1, NULL);
                    
                    // Set caps object and thereby capture/recording framerate:
                    g_object_set(G_OBJECT(camera), "viewfinder-caps", caps, NULL);

                    // Free "used up" caps object with writable viewfinder caps:
                    gst_caps_unref(caps);
                    
                    // If recording is active we want to add the "framerate" property to the recording
                    // caps as well:
                    if (capdev->recording_active) {
                        // Need to create a fresh copy of the current recording caps, so it is mutable
                        // and we can add the new "framerate" constraint/property.                        
                        g_object_get(G_OBJECT(camera), "video-capture-caps", &capsi, NULL);
                        caps = gst_caps_copy(capsi);
                        gst_caps_unref(capsi);

                        // Assign our validated recording framerate to recording caps:
                        gst_caps_set_simple(caps, "framerate", GST_TYPE_FRACTION, (int)(capturerate + 0.5), 1, NULL);

                        // Assign new caps as video recording caps:
                        g_object_set(G_OBJECT(camera), "video-capture-caps", caps, NULL);

                        // Free "used up" caps object with writable video recording caps:
                        gst_caps_unref(caps);
                    }
                }
                else {
                    // No. Play it safe and just refrain from setting a capture framerate. This will
                    // run the capture and recording at whatever the default framerate of the device is.
                    // However, warn user about unsupported setting and our fallback solution:
                    if (PsychPrefStateGet_Verbosity() > 1) {
                        printf("PTB-WARNING: Video device %i does not support requested framerate %i fps at current settings for video resolution\n", capturehandle, (int)(capturerate + 0.5));
                        printf("PTB-WARNING: and color depth. Will workaround this by trying to capture at the default framerate of the device for current settings.\n");
                    }

                    // Free unneeded caps object with writable viewfinder caps:
                    gst_caps_unref(caps);
                }
            }

			if (!PsychVideoPipelineSetState(camera, GST_STATE_PAUSED, 10.0)) {
				// Failed!
				PsychGSProcessVideoContext(capdev->VideoContext, FALSE);
				PsychErrorExitMsg(PsychError_user, "Failure in pipeline transition ready -> paused - Start of video capture failed!");
			}
		}

		// Setup of capture framerate & resolution:
		if (usecamerabin) {
			// Set requested capture/recording video resolution and framerate:
            if (usecamerabin == 1) {
                // Camerabin 1:
                // Map special capturerate value DBL_MAX to a fps nominator of zero. This
                // asks the engine to capture at the maximum supported framerate for given
                // format and resolution:
                g_signal_emit_by_name (G_OBJECT(camera),
                                       "set-video-resolution-fps", capdev->width, capdev->height,
                                       ((capturerate < DBL_MAX) ? (int)(capturerate + 0.5) : 0), 1);
            }
		}
		else {
			// Set playback rate in non camerabin configuration:
			// We only set it if the special value DBL_MAX isn't provided. DBL_MAX means
			// to run at maximum supported fps. As we can't pass this request to playbin2,
			// we do the 2nd best and don't ask for a specific fps at all, hoping that playbin2
			// will then configure to a save high speed default value -- or at least to something
			// supported and safe:
			if (capturerate < DBL_MAX) {
				gst_element_seek(camera, capturerate, GST_FORMAT_TIME, GST_SEEK_FLAG_FLUSH | GST_SEEK_FLAG_ACCURATE, GST_SEEK_TYPE_NONE, 0, GST_SEEK_TYPE_NONE, 0);
			}
		}

		// Init parameters:
		capdev->last_pts = -1.0;
		capdev->nr_droppedframes = 0;
		capdev->frameAvail = 0;
		capdev->preRollAvail = 0;

		// Wait until start deadline reached:
		if (*startattime != 0) PsychWaitUntilSeconds(*startattime);
		
		// Start actual video capture and/or recording:
		if (PsychPrefStateGet_Verbosity()>5) printf("PTB-DEBUG: Starting capture...\n");

		if (!PsychVideoPipelineSetState(camera, GST_STATE_PLAYING, 10.0)) {
			// Failed!
			PsychGSProcessVideoContext(capdev->VideoContext, FALSE);
			PsychErrorExitMsg(PsychError_user, "Failure in pipeline transition paused -> playing - Start of video capture failed!");
		}
		
		// Start video recording if requested:
		if (usecamerabin && capdev->recording_active) {
			if (PsychPrefStateGet_Verbosity()>5) printf("PTB-DEBUG: Starting recording...\n");
			g_object_set(G_OBJECT(camera), "mode", (usecamerabin == 1) ? 1 : 2, NULL);
			if (PsychPrefStateGet_Verbosity()>5) printf("PTB-DEBUG: Recording 1 started...\n");
			g_object_set(G_OBJECT(camera), (usecamerabin == 1) ? "filename" : "location", capdev->targetmoviefilename, NULL);
			if (PsychPrefStateGet_Verbosity()>5) printf("PTB-DEBUG: Recording 2 started...\n");
			g_signal_emit_by_name (camera, (usecamerabin == 1) ? "capture-start" : "start-capture", 0, 0);
			if (PsychPrefStateGet_Verbosity()>5) printf("PTB-DEBUG: Recording started...\n");
		}

        // Should we dump the whole encoding pipeline graph to a file for visualization
        // with GraphViz? This can be controlled via PsychTweak('GStreamerDumpFilterGraph' dirname);
        if (getenv("GST_DEBUG_DUMP_DOT_DIR")) {
            // Dump complete capture/recording filter graph to a .dot file for later visualization with GraphViz:
            printf("PTB-DEBUG: Dumping video capture/recording graph to directory %s.\n", getenv("GST_DEBUG_DUMP_DOT_DIR"));
            GST_DEBUG_BIN_TO_DOT_FILE_WITH_TS(GST_BIN(camera), GST_DEBUG_GRAPH_SHOW_ALL, "PsychVideoCaptureGraph");
        }

		// Wait for real start of capture, i.e., arrival of 1st captured
		// video buffer:
		PsychLockMutex(&capdev->mutex);

		// This wait is just here to clear potentially pending conditions on
		// MS-Windows. Therefore only wait for 1.1 msecs at most:
		PsychTimedWaitCondition(&capdev->condition, &capdev->mutex, 0.0011);

		// Wait for first frame to become available, but only if return of videodata is
		// actually enabled:
		if (!(capdev->recordingflags & 4)) {
			while (capdev->frameAvail == 0) {
				if (PsychPrefStateGet_Verbosity() > 5) {
					printf("PTB-DEBUG: Waiting for real start: fA = %i pA = %i fps=%f\n", capdev->frameAvail, capdev->preRollAvail, capdev->fps);
					fflush(NULL);
				}
				PsychTimedWaitCondition(&capdev->condition, &capdev->mutex, 10.0);
			}
		}

		// Assign a plausible framerate if none assigned properly:
		if (capdev->fps <= 0) capdev->fps = capturerate;

		// Ok, capture is now started:
		capdev->grabber_active = 1;

		PsychUnlockMutex(&capdev->mutex);

		// Record real start time:
		PsychGetAdjustedPrecisionTimerSeconds(startattime);

		// Some processing time for message dispatch:
		PsychGSProcessVideoContext(capdev->VideoContext, FALSE);

		if (PsychPrefStateGet_Verbosity()>5) printf("PTB-DEBUG: Capture engine fully running...\n");
		
		/* // Allocate conversion buffer if needed for YUV->RGB conversions.
		if (capdev->pixeldepth == -1) {
			// Not used at the moment!!
			// Software conversion of YUV -> RGB needed. Allocate a proper scratch-buffer:
			capdev->scratchbuffer = malloc(capdev->frame_width * capdev->frame_height * 3);
		}
		*/

		if(PsychPrefStateGet_Verbosity() > 3) {
			printf("PTB-INFO: Capture started on device %i - Input video resolution %i x %i - Framerate: %f fps.\n",
			       capturehandle, capdev->width, capdev->height, capdev->fps);
		}
	}
	else {
		// Stop capture:
		if (capdev->grabber_active) {
			// Store a backup copy of pipeline basetime for use in offline frame fetch:
			capdev->lastSavedBaseTime = gst_element_get_base_time(camera);

			// Video framerate converter attached?
			if (capdev->videorate_filter) {
				// Query its stats in dropped and duplicated frames:
				g_object_get(G_OBJECT(capdev->videorate_filter),
					     "in", &nrInFrames,
					     "out", &nrOutFrames,
					     "drop", &nrDroppedFrames,
					     "duplicate", &nrDuplicatedFrames,
					     NULL);
			}

			// Stop video recording if requested:
			if (usecamerabin && capdev->recording_active) {
				PsychGSProcessVideoContext(capdev->VideoContext, FALSE);

				if(PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: StopVideoCapture: Stopping video recording...\n");
				g_signal_emit_by_name (camera, (usecamerabin == 1) ? "capture-stop" : "stop-capture", 0);

				if(PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: StopVideoCapture: Stopping pipeline [playing -> paused]\n");
				if (!PsychVideoPipelineSetState(camera, GST_STATE_PAUSED, 10.0)) {
					if(PsychPrefStateGet_Verbosity() > 0) {
						PsychGSProcessVideoContext(capdev->VideoContext, FALSE);
						printf("PTB-ERROR: StopVideoCapture: Unable to pause recording pipeline! Prepare for trouble!\n");
					}
				}

				if(PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: StopVideoCapture: Videorecording stopped.\n");
			}

			// Stop pipeline, bring it back into READY state:
			if(PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: StopVideoCapture: Fully stopping and shutting down capture pipeline.\n");
			if (!PsychVideoPipelineSetState(camera, (usecamerabin) ? GST_STATE_READY : GST_STATE_PAUSED, 10.0)) {
				if(PsychPrefStateGet_Verbosity() > 0) {
					PsychGSProcessVideoContext(capdev->VideoContext, FALSE);
					printf("PTB-ERROR: StopVideoCapture: Unable to stop capture pipeline [Transition 1 to ready state failed]! Prepare for trouble!\n");
				}
			}
			else {
				// Capture stopped. More cleanup work:
				PsychGSProcessVideoContext(capdev->VideoContext, FALSE);
				if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: StopVideoCapture: Stopped.\n");

				// Drain any queued buffers, if requested:
				if (dropframes) {
					drainedCount = PsychGSDrainBufferQueue(capdev, INT_MAX, 0);
					if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: StopVideoCapture: Drain of %i buffers completed. GStreamer idle, cleaning up...\n", drainedCount);
				}
			}

			// Ok, capture is now stopped.
			capdev->grabber_active = 0;
			
			if (capdev->scratchbuffer) {
				// Release scratch-buffer:
				free(capdev->scratchbuffer);
				capdev->scratchbuffer = NULL;
			}

			if (PsychPrefStateGet_Verbosity() > 3) {
				// Output count of dropped frames:
				if ((dropped=capdev->nr_droppedframes) > 0) {
					printf("PTB-INFO: Video live capture dropped at least %i frames on device %i to keep capture running in sync with realtime.\n", dropped, capturehandle); 
				}

				if (capdev->videorate_filter) {
					printf("PTB-INFO: Framerate stabilization received %i frames, delivered %i frames, had to drop %i frames and duplicate %i frames\n",
					       (int) nrInFrames, (int) nrOutFrames, (int) nrDroppedFrames, (int) nrDuplicatedFrames);
					printf("PTB-INFO: to compensate for framerate fluctuations or mismatch between expected vs. real framerate on device %i.\n", capturehandle);
				}

				if (capdev->nrframes>0) capdev->avg_decompresstime/= (double) capdev->nrframes;
				printf("PTB-INFO: Average time spent in video decompressor (waiting/polling for new frames) was %f milliseconds.\n", (float) capdev->avg_decompresstime * 1000.0f);
				if (capdev->nrgfxframes>0) capdev->avg_gfxtime/= (double) capdev->nrgfxframes;
				printf("PTB-INFO: Average time spent in GetCapturedImage (intensity calculation Video->OpenGL texture conversion) was %f milliseconds.\n",  (float) capdev->avg_gfxtime * 1000.0f);
			}
		}
	}
	
	fflush(NULL);

	// Reset framecounters and statistics:
	capdev->nrframes = 0;
	capdev->avg_decompresstime = 0;
	capdev->nrgfxframes = 0;
	capdev->avg_gfxtime = 0;
	
	// Return either the real capture framerate (at start of capture) or count of dropped frames - at end of capture.
	return((capturerate!=0) ? (int) (capdev->fps + 0.5) : dropped);
}


/*  
*  PsychGSGetTextureFromCapture() -- Create an OpenGL texturemap from a specific videoframe from given capture object.
*
*  win = Window pointer of onscreen window for which a OpenGL texture should be created.
*  capturehandle = Handle to the capture object.
*  checkForImage = >0 == Just check if new image available, 0 == really retrieve the image, blocking if necessary.
*                   2 == Check for new image, block inside this function (if possible) if no image available.
*
*  timeindex = This parameter is currently ignored and reserved for future use.
*  out_texture = Pointer to the Psychtoolbox texture-record where the new texture should be stored.
*  presentation_timestamp = A ptr to a double variable, where the presentation timestamp of the returned frame should be stored.
*  summed_intensity = An optional ptr to a double variable. If non-NULL, then sum of intensities over all channels is calculated and returned.
*  outrawbuffer = An optional ptr to a memory buffer of sufficient size. If non-NULL, the buffer will be filled with the captured raw image data, e.g., for use inside Matlab or whatever...
*  Returns Number of pending or dropped frames after fetch on success (>=0), -1 if no new image available yet, -2 if no new image available and there won't be any in future.
*/
int PsychGSGetTextureFromCapture(PsychWindowRecordType *win, int capturehandle, int checkForImage, double timeindex,
								 PsychWindowRecordType *out_texture, double *presentation_timestamp, double* summed_intensity, rawcapimgdata* outrawbuffer)
{
    PsychVidcapRecordType *capdev;
    GstBuffer *videoBuffer = NULL;
    gint64 bufferIndex;
    double deltaT = 0;
    GstClockTime baseTime;

    int waitforframe;
    int w, h;
    unsigned int intensity = 0;
    unsigned int count, i, bpp;
    unsigned char* pixptr;
    psych_bool newframe = FALSE;
    double tstart, tend;
    int nrdropped = 0;
    unsigned char* input_image = NULL;


    // Make sure GStreamer is ready:
    PsychGSCheckInit("videocapture");

    // Take start timestamp for timing stats:
    PsychGetAdjustedPrecisionTimerSeconds(&tstart);
	
    // Retrieve device record for handle:
    capdev = PsychGetGSVidcapRecord(capturehandle);

    // Allow context task to do its internal bookkeeping and cleanup work:
    PsychGSProcessVideoContext(capdev->VideoContext, FALSE);

    // If this is a pure audio capture with no video channels, we always return failed,
    // as those certainly don't have movie frames associated.
    if ((capdev->recordingflags & 4) || (capdev->nrVideoTracks == 0)) return((checkForImage) ? -1 : FALSE);

    if (!PsychIsOnscreenWindow(win)) {
        PsychErrorExitMsg(PsychError_user, "You need to pass a handle to an onscreen window, not to something else!");
    }

    // Compute width and height for later creation of textures etc. Need to do this here,
    // so we can return the values for raw data retrieval:
    w = capdev->frame_width;
    h = capdev->frame_height;

    // Size of a single pixel in bytes:
    bpp = capdev->reqpixeldepth;
	
    // If a outrawbuffer struct is provided, we fill it with info needed to allocate a
    // sufficient memory buffer for returned raw image data later on:
    if (outrawbuffer) {
	    outrawbuffer->w = w;
	    outrawbuffer->h = h;
	    outrawbuffer->depth = bpp;
    }
	
    waitforframe = (checkForImage > 1) ? 1:0; // Blocking wait for new image requested?
	
    // A checkForImage 4 means "no op" with the GStreamer capture engine: This is meant to drive
    // a single threaded movie recording engine, ie., grant processing time to it.
    // Our GStreamer engine doesn't need this as it is highly multi-threaded and does all
    // relevant work in the background.
    if (checkForImage == 4) return(0);
	
    // Should we just check for new image?
    if (checkForImage) {
	    // Reset current dropped count to zero:
	    capdev->current_dropped = 0;
	    
	    PsychLockMutex(&capdev->mutex);
	    if (!capdev->frameAvail) {
		    // No new frame available yet:

		    if (capdev->grabber_active == 0) {
			    // Grabber stopped. We'll never get a new image:
			    PsychUnlockMutex(&capdev->mutex);
			    return(-2);
		    }

		    // Blocking wait requested?
		    if (!waitforframe) {
			    // No blocking wait, just return "no new frame yet":
			    PsychUnlockMutex(&capdev->mutex);
			    //printf("PTB-DEBUG: NO NEW FRAME\n");
			    return(-1);
		    }

		    // No new frame available. Perform a blocking wait:
		    PsychTimedWaitCondition(&capdev->condition, &capdev->mutex, 10.0);

		    // Recheck:
		    if (!capdev->frameAvail) {
			    // Game over! Wait timed out after 10 secs.
			    PsychUnlockMutex(&capdev->mutex);
			    if (PsychPrefStateGet_Verbosity()>4) printf("PTB-ERROR: In blocking wait: No new video frame received after timeout of 10 seconds! Aborting fetch.\n");
			    return(-1);
		    }

		    // Have at least one new frame after wait :-)
	    }
	    
	    // New frame available:

	    // Store count of currently queued frames (in addition to the one just fetched).
	    // This is an indication of how well the users script is keeping up with the video stream,
	    // technically the number of frames that would need to be dropped to keep in sync with the
	    // stream.
	    capdev->current_dropped = capdev->frameAvail - 1;

	    //printf("PTB-DEBUG: NEW FRAME %d\n", capdev->frameAvail);
	    
	    // Ok, at least one new frame ready. If more than one frame has queued up and
	    // we are in 'dropframes' mode, ie. we should always deliver the most recent available
	    // frame, then we quickly fetch & discard all queued frames except the last one.
	    while((capdev->dropframes) && ((int) capdev->frameAvail > 1)) {
		    // On frame to dequeue:
		    capdev->frameAvail--;
		    PsychUnlockMutex(&capdev->mutex);

		    // Pull next buffer which is to be dropped:
		    videoBuffer = gst_app_sink_pull_buffer(GST_APP_SINK(capdev->videosink));

		    // Release the capture buffer. Return it to the DMA ringbuffer pool:
		    gst_buffer_unref(videoBuffer);
		    videoBuffer = NULL;

		    // Retry if there are more pending frames to drop:
		    PsychLockMutex(&capdev->mutex);
	    }

	    // Release lock:
	    PsychUnlockMutex(&capdev->mutex);

	    // Update stats for decompression:
	    PsychGetAdjustedPrecisionTimerSeconds(&tend);
	    
	    // Increase counter of decompressed frames:
	    capdev->nrframes++;
	    
	    // Update avg. decompress time:
	    capdev->avg_decompresstime+= (tend - tstart);
	    
	    // Return availability status: At least one new frame available :-)
	    return(0);
    }
    
    // This point is only reached if checkForImage == FALSE, which only happens
    // if a new frame is available in our buffer:    
    PsychLockMutex(&capdev->mutex);

    //printf("PTB-DEBUG: Blocking fetch start %d\n", capdev->frameAvail);
    if (!capdev->frameAvail) {
		// No new frame available but grabber active. Perform a blocking wait:
		if (capdev->grabber_active) PsychTimedWaitCondition(&capdev->condition, &capdev->mutex, 10.0);
		
		// Recheck:
		if (!capdev->frameAvail) {
			// Game over! Wait timed out after 10 secs.
			PsychUnlockMutex(&capdev->mutex);
			if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: No new video frame received after timeout of 10 seconds! Something's wrong. Aborting fetch.\n");
			return(-1);
		}
		
		// At this point we should have at least one frame available.
    }
	
    // We're here with at least one frame available and the mutex lock held.
    if (PsychPrefStateGet_Verbosity()>6) printf("PTB-DEBUG: Pulling from videosink, %d buffers avail...\n", capdev->frameAvail);

    // One less frame available after our fetch:
    capdev->frameAvail--;
    
    // Clamp frameAvail to queue lengths, unless queue length is set to zero, which means "unlimited":
    if ((int) gst_app_sink_get_max_buffers(GST_APP_SINK(capdev->videosink)) < capdev->frameAvail) {
	    if (gst_app_sink_get_max_buffers(GST_APP_SINK(capdev->videosink)) > 0) {
		    capdev->frameAvail = gst_app_sink_get_max_buffers(GST_APP_SINK(capdev->videosink));
	    }
    }
    if (PsychPrefStateGet_Verbosity()>6) printf("PTB-DEBUG: Post-Pulling from videosink, %d buffers avail...\n", capdev->frameAvail);
    
    // This will pull the oldest video buffer from the videosink. It would block if none were available,
    // but that won't happen as we wouldn't reach this statement if none were available. It would return
    // NULL if the stream would be EOS or the pipeline off, but that shouldn't ever happen:
    videoBuffer = gst_app_sink_pull_buffer(GST_APP_SINK(capdev->videosink));

    // We can unlock early, thanks to videosink's internal buffering:
    PsychUnlockMutex(&capdev->mutex);
    
    if (videoBuffer) {
	    // Assign pointer to videoBuffer's data directly: Avoids one full data copy compared to oldstyle method.
	    // input_image points to the image buffer in our cam:
	    input_image = (unsigned char*) (GLuint*) GST_BUFFER_DATA(videoBuffer);
	    
	    // Assign pts presentation timestamp in pipeline stream time and convert to seconds:
	    if (capdev->recordingflags & 64) {
		    // Retrieve raw buffer timestamp - pipeline running time.
		    capdev->current_pts = (double) GST_BUFFER_TIMESTAMP(videoBuffer) / (double) 1e9;
	    } else {
		    // Add base time to convert running time buffer timestamp into absolute time:
		    baseTime = gst_element_get_base_time(capdev->camera);
		    if (baseTime == 0) baseTime = capdev->lastSavedBaseTime;

		    capdev->current_pts = (double) (GST_BUFFER_TIMESTAMP(videoBuffer) + baseTime) / (double) 1e9;

            // Apply corrective offset for GStreamer clock base zero point:
            capdev->current_pts+= gs_startupTime;     
        }
        
	    deltaT = 0.0;
	    if (GST_CLOCK_TIME_IS_VALID(GST_BUFFER_DURATION(videoBuffer)))
		    deltaT = (double) GST_BUFFER_DURATION(videoBuffer) / (double) 1e9;

	    if (PsychPrefStateGet_Verbosity() > 6) {
		    PsychProbeBufferProps(videoBuffer, NULL, NULL, &capdev->fps);
		    printf("Bufferprobe: newfps = %f altfps = %f\n", capdev->fps, (float) ((deltaT > 0) ? 1.0 / deltaT : 0.0));
	    }

	    bufferIndex = GST_BUFFER_OFFSET(videoBuffer);
    } else {
	    if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: No new video frame received in gst_app_sink_pull_buffer! Something's wrong. Aborting fetch.\n");
	    return(-1);
    }
    if (PsychPrefStateGet_Verbosity()>6) printf("PTB-DEBUG: ...done.\n");

    // Presentation timestamp requested?
    if (presentation_timestamp) {
	    // Return it:
	    *presentation_timestamp = capdev->current_pts;
    }
	
    // Synchronous texture fetch: Copy content of capture buffer into a texture:
    // =========================================================================
	
    // Do we want to do something with the image data and have a
    // scratch buffer for color conversion alloc'ed?
    if ((capdev->scratchbuffer) && ((out_texture) || (summed_intensity) || (outrawbuffer))) {
		// Yes. Perform color-conversion YUV->RGB from cameras DMA buffer
		// into the scratch buffer and set scratch buffer as source for
		// all further operations:
		memcpy(capdev->scratchbuffer, input_image, w * h * bpp);
		
		// Ok, at this point we should have a RGB8 texture image ready in scratch_buffer.
		// Set scratch buffer as our new image source for all further processing:
		input_image = (unsigned char*) capdev->scratchbuffer;
    }
	
    // Only setup if really a texture is requested (non-benchmarking mode):
    if (out_texture) {
	    // Activate OpenGL context of target window:
	    PsychSetGLContext(win);
	    
        #if PSYCH_SYSTEM == PSYCH_OSX
            // Explicitely disable Apple's Client storage extensions. For now they are not really useful to us.
            glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, GL_FALSE);
        #endif

        // Set default rect of texture:
	    PsychMakeRect(out_texture->rect, 0, 0, w, h);    
        PsychCopyRect(out_texture->clientrect, out_texture->rect);

	    // Set NULL - special texture object as part of the PTB texture record:
	    out_texture->targetSpecific.QuickTimeGLTexture = NULL;
	    
	    // Set texture orientation as if it were an inverted Offscreen window: Upside-down.
	    out_texture->textureOrientation = 3;
	    
	    // Setup a pointer to our buffer as texture data pointer: Setting memsize to zero
	    // prevents unwanted free() operation in PsychDeleteTexture...
	    out_texture->textureMemorySizeBytes = 0;
	    
	    // Set default texture depth: Could be 8, 16, 24 or 32 bpp.
	    out_texture->depth = capdev->reqpixeldepth * 8;
	    
	    // 4-channel textures are aligned on 4 Byte boundaries because texels are RGBA8. If they
        // have an even number of pixels (width) per row, they are even 8 Byte aligned. Otherwise
        // we play safe and assume no alignment, ie., 1 Byte alignment:
	    out_texture->textureByteAligned = (capdev->reqpixeldepth == 4) ? ((w % 2) ? 4 : 8) : 1;

	    // This will retrieve an OpenGL compatible pointer to the pixel data and assign it to our texmemptr:
	    out_texture->textureMemory = (GLuint*) input_image;
	    
	    // Special case depths == 2, aka YCBCR texture?
	    if ((capdev->reqpixeldepth == 2) && (capdev->pixeldepth == 16) && (win->gfxcaps & kPsychGfxCapUYVYTexture)) {
            // GPU supports UYVY textures and we get data in that YCbCr format. Tell
            // texture creation routine to use this optimized format:
            if (!glewIsSupported("GL_APPLE_ycbcr_422")) {
                // No support for more powerful Apple extension. Use Linux MESA extension:
                out_texture->textureinternalformat = GL_YCBCR_MESA;
                out_texture->textureexternalformat = GL_YCBCR_MESA;
            } else {
                // Apple extension supported:
                out_texture->textureinternalformat = GL_RGB8;
                out_texture->textureexternalformat = GL_YCBCR_422_APPLE;
            }
            // Same enumerant for Apple and Mesa:
            out_texture->textureexternaltype   = GL_UNSIGNED_SHORT_8_8_MESA;

            // Number of effective channels is 3 for RGB8:
            out_texture->nrchannels = 3;
            
            // And 24 bpp depth:
            out_texture->depth = 24;

            // Byte alignment: For even number of pixels, assume at least 4 Byte alignment due to packing of 2 effective
            // pixels into one 32-Bit packet, maybe even 8 Byte alignment if divideable by 4. For other width's, assume
            // no alignment ie., 1 Byte:
            out_texture->textureByteAligned = (w % 2) ? 1 : ((w % 4) ? 4 : 8);	    
        }

        // YUV I420 planar pixel upload requested?
	    if ((capdev->reqpixeldepth == 2) && (capdev->pixeldepth == 12)) {
            // We encode I420 planar data inside a 8 bit per pixel luminance texture of
            // 1.5x times the height of the video frame. First the "Y" luminance plane
            // is stored at full 1 sample per pixel resolution with 8 bits. Then a 0.25x
            // height slice with "U" Cr chrominance data at half the horizontal and vertical
            // resolution aka 1 sample per 2x2 pixel quad. Then a 0.25x height slice with "V"
            // Cb chrominance data at 1 sample per 2x2 pixel quad resolution. As such the texture
            // appears to OpenGL as a normal LUMINANCE8 texture. Conversion of the planar format
            // into useable RGBA8 pixel fragments will happen during rendering via a suitable fragment
            // shader. The net gain of this is that we effectively only need 1.5 Bytes per pixel instead
            // of 3 Bytes for RGB8 or 4 Bytes for RGBA8:
			out_texture->textureexternaltype   = GL_UNSIGNED_BYTE;
			out_texture->textureexternalformat = GL_LUMINANCE;
			out_texture->textureinternalformat = GL_LUMINANCE8;
            
			// Define a rect of 1.5 times the video frame height, so PsychCreateTexture() will source
            // the whole input data buffer:
			PsychMakeRect(out_texture->rect, 0, 0, w, h * 1.5);
            
            // Check if 1.5x height texture fits within hardware limits of this GPU:
            if (h * 1.5 > win->maxTextureSize) PsychErrorExitMsg(PsychError_user, "Videoframe size too big for this graphics card and pixelFormat! Please retry with a pixeldepth of 4 in 'OpenVideoCapture'.");
            
            // Create planar "I420 inside L8" texture:
            PsychCreateTexture(out_texture);
            
            // Restore rect and clientrect of texture to effective size of video frame:
            PsychMakeRect(out_texture->rect, 0, 0, w, h);
            PsychCopyRect(out_texture->clientrect, out_texture->rect);
            
            // Mark texture as planar encoded, so proper conversion shader gets applied during
            // call to PsychNormalizeTextureOrientation(), prior to any render-to-texture operation, e.g.,
            // if used as an offscreen window, or as a participant of a Screen('TransformTexture') call:
            out_texture->specialflags |= kPsychPlanarTexture;
            
            // Assign special filter shader for sampling and color-space conversion of the
            // planar texture during drawing or PsychNormalizeTextureOrientation():
            if (!PsychAssignPlanarI420TextureShader(out_texture, win)) PsychErrorExitMsg(PsychError_user, "Assignment of I420 video decoding shader failed during video texture creation!");
            
            // Number of effective channels is 3 for RGB8:
            out_texture->nrchannels = 3;
            
            // And 24 bpp depth:
            out_texture->depth = 24;
            
            // Byte alignment: Assume no alignment for now:
            out_texture->textureByteAligned = 1;
        }
        else {
            // Let PsychCreateTexture() do the rest of the job of creating, setting up and
            // filling an OpenGL texture with content:
            PsychCreateTexture(out_texture);
        }

        // Immediate conversion of texture into normalized orientation and format requested
        // by usercode?
        if (capdev->recordingflags & 2048) {
            // Transform out_texture video texture into a normalized, upright texture if it isn't already in
            // that format. We require this standard orientation for simplified shader design.
            PsychSetShader(win, 0);
            PsychNormalizeTextureOrientation(out_texture);            
        }

	    // Ready to use the texture...
    }
    
    // Sum of pixel intensities requested?
    if(summed_intensity) {
	    pixptr = (unsigned char*) input_image;
	    count  = w * h * bpp;
	    for (i=0; i<count; i++) intensity+=(unsigned int) pixptr[i];
	    *summed_intensity = ((double) intensity) / w / h / bpp;
    }
    
    // Raw data requested?
    if (outrawbuffer) {
	    // Copy it out:
	    outrawbuffer->w = w;
	    outrawbuffer->h = h;
	    outrawbuffer->depth = bpp;
	    count = (w * h * outrawbuffer->depth);
	    memcpy(outrawbuffer->data, (const void *) input_image, count);
    }
	
    // Release the capture buffer. Return it to the DMA ringbuffer pool:
    gst_buffer_unref(videoBuffer);
    videoBuffer = NULL;

    // Update total count of dropped (or pending) frames:
    capdev->nr_droppedframes += capdev->current_dropped;
    nrdropped = capdev->current_dropped;
    capdev->current_dropped = 0;
	
    // Timestamping:
    PsychGetAdjustedPrecisionTimerSeconds(&tend);
	
    // Increase counter of retrieved textures:
    capdev->nrgfxframes++;
	
    // Update average time spent in texture conversion:
    capdev->avg_gfxtime+= (tend - tstart);
    
    // We're successfully done! Return number of dropped (or pending in DMA ringbuffer) frames:
    return(nrdropped);
}

/* Set capture device specific parameters:
* Currently, the named parameters are a subset of the parameters supported by the
* IIDC specification, mapped to more convenient names.
*
* Input: pname = Name string to specify the parameter.
*        value = Either DBL_MAX to not set but only query the parameter, or some other
*                value, that we try to set in the Firewire camera.
*
* Returns: Old value of the setting
*/
double PsychGSVideoCaptureSetParameter(int capturehandle, const char* pname, double value)
{
	int intval, oldintval;
	float oldfvalue = FLT_MAX;
	double oldvalue = DBL_MAX; // Initialize return value to the "unknown/unsupported" default.
	psych_bool assigned = false;
	psych_bool present  = false;
	GstColorBalance* cb = NULL;
	GList* cl = NULL;
	GList* iter = NULL;
	GstColorBalanceChannel* cc = NULL;

	// Retrieve device record for handle:
	PsychVidcapRecordType* capdev = PsychGetGSVidcapRecord(capturehandle);

	// Make sure GStreamer is ready:
	PsychGSCheckInit("videocapture");

	oldintval = 0xFFFFFFFF;
	
	// Round value to integer:
	intval = (int) (value + 0.5);
	
	// Set a new target movie name for video recordings:
	if (strstr(pname, "SetNewMoviename=")) {
		// Find start of movie namestring and assign to pname:
		pname = strstr(pname, "=");
		pname++;

		// Child protection:
		if (!capdev->recording_active) {
			if (PsychPrefStateGet_Verbosity() > 1) {
				printf("PTB-WARNING: Tried to change name of target movie file on device %i, but recording not enabled on that device! Ignored.\n", capturehandle);
			}
			return(-2);
		}

		// Can't reassign new codec without reopening the device:
		if (strstr(pname, ":CodecType")) {
			*(strstr(pname, ":CodecType")) = 0;
			if (PsychPrefStateGet_Verbosity() > 1) {
				printf("PTB-WARNING: Tried to change recording codec on device %i, but this isn't possible without reopening the device. Ignored.\n", capturehandle);
			}
		}

		// Release old movie name:
		if (capdev->targetmoviefilename) free(capdev->targetmoviefilename);
		capdev->targetmoviefilename = NULL;

		// Assign new movie name:
		capdev->targetmoviefilename = strdup(pname);

		if (PsychPrefStateGet_Verbosity() > 2) {
			printf("PTB-INFO: Changed name of movie file for recording on device %i to '%s'.\n", capturehandle, pname);
		}

		return(0);
	}

	// Check parameter name pname and call the appropriate subroutine:
	if (strcmp(pname, "TriggerCount")==0 || strcmp(pname, "WaitTriggerCount")==0) {
		// Query of cameras internal trigger counter or waiting for a specific
		// value in the counter requested. Trigger counters are special features,
		// (so called "Smart Features" or "Advanced Features" in the IIDC spec)
		// which are only available on selected cameras. They are not currently
		// available on GStreamer.
		return(-2);
	}
	
    // Check if GstColorBalanceInterface is supported and assign it for use downstream. Probe
    // different providers: camerabin1 (should support it), camerabin2 (doesn't at this point in time),
    // the wrappercamerabinsrc of camerabin2 (doesn't at this point in time), the video source attached
    // to the wrappercamerabinsrc (currently the video4linux2 source does support it):
    if (usecamerabin) {
        cb = NULL;
        // Probe camerabin1 / camerabin2:
        if (!cb && gst_element_implements_interface(capdev->camera, GST_TYPE_COLOR_BALANCE)) cb = GST_COLOR_BALANCE(capdev->camera);
        // If fail, probe wrappercamerabinsrc:
        if (!cb && gst_element_implements_interface(capdev->videowrappersrc, GST_TYPE_COLOR_BALANCE)) cb = GST_COLOR_BALANCE(capdev->videowrappersrc);
        // If fail, probe videosource itself:
        if (!cb && gst_element_implements_interface(capdev->videosource, GST_TYPE_COLOR_BALANCE)) cb = GST_COLOR_BALANCE(capdev->videosource);
        // If fail, game over:
        if (!cb && (PsychPrefStateGet_Verbosity() > 3)) printf("PTB-WARNING: Video device %i does not suppport GstColorBalance interface as expected.\n", capturehandle);
	}

	if (strcmp(pname, "PrintParameters")==0) {
		// Special command: List and print all features...
		printf("PTB-INFO: The video source provides the following controllable parameters:\n");
		printf("PTB-INFO: ----------------------------------------------------------------\n\n");
		printf("PTB-INFO: Optional parameters - may or may not be supported:\n");
		printf("PTB-INFO: Shutter, Aperture, EVCompensation, Flickermode, Whitebalancemode,\n");
		printf("PTB-INFO: Flashmode, Scenemode, Focusmode\n\n");
		printf("PTB-INFO: These are definitely supported by the connected camera:\n");


		if (cb) {
			// Enumerate all color balance channels:
			cl = (GList*) gst_color_balance_list_channels(cb);
			for (iter = g_list_first(cl); iter != NULL ; iter = g_list_next(iter)) {
				cc = (GstColorBalanceChannel*) iter->data;
				printf("PTB-INFO: '%s'\t\t min=%i\t : max=%i\t : current=%i\n", (char*) cc->label,
				(int) cc->min_value, (int) cc->max_value, (int) gst_color_balance_get_value(cb, cc));
			}
		}
		printf("PTB-INFO: ----------------------------------------------------------------\n\n");
		return(0);
	}

	// Return current framerate:
	if (strcmp(pname, "GetFramerate")==0) {
		PsychCopyOutDoubleArg(1, FALSE, capdev->fps);
		return(0);
	}
	
	// Return current ROI of camera, as requested (and potentially modified during
	// PsychOpenCaptureDevice(). This is a read-only parameter, as the ROI can
	// only be set during Screen('OpenVideoCapture').
	if (strcmp(pname, "GetROI")==0) {
		PsychCopyOutRectArg(1, FALSE, capdev->roirect);
		return(0);
	}
	
	// Return vendor name string:
	if (strcmp(pname, "GetVendorname")==0) {
		PsychCopyOutCharArg(1, FALSE, "Unknown Vendor");
		return(0);
	}
	
	// Return model name string:
	if (strcmp(pname, "GetModelname")==0) {
		PsychCopyOutCharArg(1, FALSE, capdev->cameraFriendlyName);
		return(0);
	}

	// All code below this check is for camerabin only:
	if (!usecamerabin) {
		// No camerabin, no way to query this stuff. Just fail
		return(DBL_MAX);
	}

	if (strstr(pname, "Gain") && strstr(&capdev->videosourcename[0], "aravissrc")) {
        // Special case aravissrc plugin?
        
		// Query old "gain" setting, which is integer and in units of dB Dezibel:
		g_object_get(capdev->videosource, "gain", &oldintval, NULL);
		oldvalue = (double) oldintval;
        
		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) {
            value = DBL_MAX;
            g_object_set(capdev->videosource, "gain-auto", TRUE, NULL);
        }
        
		// Optionally set new setting and switch to manual mode:
		if (value != DBL_MAX) {
            g_object_set(capdev->videosource, "gain-auto", FALSE, NULL);
            g_object_set(capdev->videosource, "gain", (int) (value + 0.5), NULL);
        }
		return(oldvalue);
	}

	if (strstr(&capdev->videosourcename[0], "aravissrc")) {
        // aravissrc plugin supports the following additional parameters:
        if (strstr(pname, "offset-x") || strstr(pname, "offset-y") || strstr(pname, "h-binning") || strstr(pname, "v-binning")) {
            // Query old setting, which is an integer:
            g_object_get(capdev->videosource, pname, &oldintval, NULL);
            oldvalue = (double) oldintval;

            // Optionally set new setting:
            if (value != DBL_MAX) {
                g_object_set(capdev->videosource, pname, (int) (value + 0.5), NULL);
            }
            return(oldvalue);
        }
	}
    
	if (strstr(pname, "Shutter")!=0) {
        // Special case aravissrc plugin?
        if (strstr(&capdev->videosourcename[0], "aravissrc")) {
            // exposure is a double value, expressing exposure time in microseconds (1e6):
            g_object_get(capdev->videosource, "exposure", &oldvalue, NULL);
            oldvalue = oldvalue / 1e6;
            
            // Reset to auto-mode, if requested:
            if (strstr(pname, "Auto")) {
                value = DBL_MAX;
                g_object_set(capdev->videosource, "exposure-auto", TRUE, NULL);
            }
            
            // Optionally set new setting:
            if (value != DBL_MAX) {
                g_object_set(capdev->videosource, "exposure-auto", FALSE, NULL);
                g_object_set(capdev->videosource, "exposure", (double) (value * 1e6), NULL);
            }
            return(oldvalue);
        }
        
        // Standard case: Exposure is an integer in nanoseconds (1e9):
        
		// Query old "exposure" setting, which is duration of shutter open:
		g_object_get(capdev->camera, "exposure", &oldintval, NULL);
		oldvalue = (double) oldintval / 1e9;

		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) value = 0;

		// Optionally set new setting:
		if (value != DBL_MAX) g_object_set(capdev->camera, "exposure", (int) (value * 1e9), NULL);
		return(oldvalue);
	}

	if (strstr(pname, "Aperture")) {
		// Query old "aperture" setting, which is the amount of lens opening:
		g_object_get(capdev->camera, "aperture", &oldintval, NULL);
		oldvalue = (double) oldintval;

		// Optionally set new setting:
		if (intval < 0) intval = 0;
		if (intval > 255) intval = 255;

		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) { intval = 0; value = 0; }

		if (value != DBL_MAX) g_object_set(capdev->camera, "aperture", intval, NULL);
		return(oldvalue);
	}

	if (strstr(pname, "EVCompensation")!=0) {
		// Query old "ev-compensation" setting:
		g_object_get(capdev->camera, "ev-compensation", &oldfvalue, NULL);
		oldvalue = (double) oldfvalue;

		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) value = 0;

		// Optionally set new setting:
		if (value != DBL_MAX) g_object_set(capdev->camera, "ev-compensation", (float) value, NULL);
		return(oldvalue);
	}

	if (strstr(pname, "Flickermode")) {
		// Query old "flicker-mode" setting:
		g_object_get(capdev->camera, "flicker-mode", &oldintval, NULL);
		oldvalue = (double) oldintval;

		// Optionally set new setting:
		if (intval < 0) intval = 0;
		if (intval > 3) intval = 3;

		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) { intval = 3; value = 0; }

		if (value != DBL_MAX) g_object_set(capdev->camera, "flicker-mode", intval, NULL);
		return(oldvalue);
	}

	if (strstr(pname, "Whitebalancemode")) {
		// Query old "white-balance-mode" setting:
		g_object_get(capdev->camera, "white-balance-mode", &oldintval, NULL);
		oldvalue = (double) oldintval;

		// Optionally set new setting:
		if (intval < 0) intval = 0;
		if (intval > 5) intval = 5;

		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) { intval = 0; value = 0; }

		if (value != DBL_MAX) g_object_set(capdev->camera, "white-balance-mode", intval, NULL);
		return(oldvalue);
	}

	if (strstr(pname, "Focusmode")) {
		// Query old "focus-mode" setting:
		g_object_get(capdev->camera, "focus-mode", &oldintval, NULL);
		oldvalue = (double) oldintval;

		// Optionally set new setting:
		if (intval < 0) intval = 0;
		if (intval > 7) intval = 7;

		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) { intval = 0; value = 0; }

		if (value != DBL_MAX) g_object_set(capdev->camera, "focus-mode", intval, NULL);
		return(oldvalue);
	}

	if (strstr(pname, "Flashmode")) {
		// Query old "flash-mode" setting:
		g_object_get(capdev->camera, "flash-mode", &oldintval, NULL);
		oldvalue = (double) oldintval;


		// Optionally set new setting:
		if (intval < 0) intval = 0;
		if (intval > 4) intval = 4;

		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) { intval = 0; value = 0; }

		if (value != DBL_MAX) g_object_set(capdev->camera, "flash-mode", intval, NULL);
		return(oldvalue);
	}

	if (strstr(pname, "Scenemode")) {
		// Query old "scene-mode" setting:
		g_object_get(capdev->camera, "scene-mode", &oldintval, NULL);
		oldvalue = (double) oldintval;

		// Optionally set new setting:
		if (intval < 0) intval = 0;
		if (intval > 6) intval = 6;

		// Reset to auto-mode, if requested:
		if (strstr(pname, "Auto")) { intval = 6; value = 0; }

		if (value != DBL_MAX) g_object_set(capdev->camera, "scene-mode", intval, NULL);
		return(oldvalue);
	}

	// Not yet matched? Try if it matches one of the color channel properties
	// from the color balance interface.
	if (cb) {
		// Search all color balance channels:
		cl = (GList*) gst_color_balance_list_channels(cb);
		for (iter = g_list_first(cl); iter != NULL ; iter = g_list_next(iter)) {
			cc = (GstColorBalanceChannel*) iter->data;

			// Match?
			if (strcmp((const char*) cc->label, (const char*) pname) == 0) {
				assigned = TRUE;

				// Query and return old setting:
				oldvalue = (double) gst_color_balance_get_value(cb, cc);

				// Optionally assign new setting:
				if (intval < cc->min_value) intval = cc->min_value;
				if (intval > cc->max_value) intval = cc->max_value;
				if (value != DBL_MAX) gst_color_balance_set_value(cb, cc, intval);
			}
		}

		if (!assigned) {
			if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Screen('SetVideoCaptureParameter', ...) called with unknown parameter %s. Ignored...\n", pname);
		}

		return(oldvalue);
	}

	return(DBL_MAX);
}

#endif
