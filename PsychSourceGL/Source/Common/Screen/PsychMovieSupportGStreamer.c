/*
	PsychSourceGL/Source/Common/Screen/PsychMovieSupportGStreamer.c
	
	PLATFORMS:	All with PTB_USE_GSTREAMER defined.

	AUTHORS:

	mario.kleiner@tuebingen.mpg.de		mk	Mario Kleiner

	HISTORY:

        28.11.2010    mk      Wrote it.

	DESCRIPTION:
	
	Psychtoolbox functions for dealing with movies. This is the operating system independent
	version which uses the GStreamer media framework.

	These PsychGSxxx functions are called from the dispatcher in
	Common/Screen/PsychMovieSupport.[hc].

	TODO:

		- Fix frame-based seeking: Time base seeking works well, frame based not so much.
		- Check if the 'drop' property + max_queue property of the appsink could be used
		  in a creative way to synchronize 1st played frame/sound with 1st texture fetch.
		- dto. other uses settings for queue length.
		- Make auto-rewind and reverse playback more robust.
		- Preload into RAM - implement if possible, although preroll seems to be sufficient.
		- Avoid spin-wait polling from calling high-level code when waiting for new frames.
		- Check check check exact timing, precision, robustness, performance...

*/

#ifdef PTB_USE_GSTREAMER

#include "Screen.h"
#include <glib.h>
#include "PsychMovieSupportGStreamer.h"
#include <gst/gst.h>
#include <gst/app/gstappsink.h>

// Need to define this for playbin2 as it is not defined
// in any header file: (Expected behaviour - not a bug)
typedef enum {
  GST_PLAY_FLAG_VIDEO         = (1 << 0),
  GST_PLAY_FLAG_AUDIO         = (1 << 1),
  GST_PLAY_FLAG_TEXT          = (1 << 2),
  GST_PLAY_FLAG_VIS           = (1 << 3),
  GST_PLAY_FLAG_SOFT_VOLUME   = (1 << 4),
  GST_PLAY_FLAG_NATIVE_AUDIO  = (1 << 5),
  GST_PLAY_FLAG_NATIVE_VIDEO  = (1 << 6),
  GST_PLAY_FLAG_DOWNLOAD      = (1 << 7),
  GST_PLAY_FLAG_BUFFERING     = (1 << 8),
  GST_PLAY_FLAG_DEINTERLACE   = (1 << 9)
} GstPlayFlags;

static const psych_bool oldstyle = FALSE;
static psych_bool useYUVDecode = FALSE;

#define PSYCH_MAX_MOVIES 100
    
typedef struct {
    psych_mutex		mutex;
    psych_condition     condition;
    double		pts;
    GstElement		*theMovie;
    GMainLoop		*MovieContext;
    GstElement          *videosink;
    unsigned char	*imageBuffer;
    int			frameAvail;
    int                 preRollAvail;
    double		rate;
    int                 startPending;
    int                 endOfFetch;
    int			loopflag;
    double		movieduration;
    int			nrframes;
    double		fps;
    int			width;
    int			height;
    double              aspectRatio;
    double		last_pts;
    int			nr_droppedframes;
    int                 nrAudioTracks;
    int                 nrVideoTracks;
    char                movieLocation[FILENAME_MAX];
    char                movieName[FILENAME_MAX];
    GLuint		cached_texture;
} PsychMovieRecordType;

static PsychMovieRecordType movieRecordBANK[PSYCH_MAX_MOVIES];
static int numMovieRecords = 0;
static psych_bool firsttime = TRUE;

#if PSYCH_SYSTEM == PSYCH_OSX
extern void g_thread_init(GThreadFunctions *vtable) __attribute__((weak_import));
#endif

/*
 *     PsychGSMovieInit() -- Initialize movie subsystem.
 *     This routine is called by Screen's RegisterProject.c PsychModuleInit()
 *     routine at Screen load-time. It clears out the movieRecordBANK to
 *     bring the subsystem into a clean initial state.
 */
void PsychGSMovieInit(void)
{
    // Initialize movieRecordBANK with NULL-entries:
    int i;
    for (i=0; i < PSYCH_MAX_MOVIES; i++) {
        memset(&movieRecordBANK[i], 0, sizeof(PsychMovieRecordType));
    }    
    numMovieRecords = 0;

    #if PSYCH_SYSTEM == PSYCH_WINDOWS
    // On Windows, we need to delay-load the GLib DLL's. This loading
    // and linking will automatically happen downstream. However, if delay loading
    // would fail, we would end up with a crash! For that reason, we try here to
    // load the DLL, just to probe if the real load/link/bind op later on will
    // likely succeed. If the following LoadLibrary() call fails and returns NULL,
    // then we know we would end up crashing. Therefore we'll output some helpful
    // error-message instead:
    if ((NULL == LoadLibrary("libgstreamer-0.10.dll")) || (NULL == LoadLibrary("libgstapp-0.10.dll"))) {
        // Failed: GLib and its threading support isn't installed. This means that
        // GStreamer won't work as the relevant .dll's are missing on the system.
        // We silently return, skpipping the GLib init, as it is completely valid
        // for a Windows installation to not have GStreamer installed at all.
        return;
    }
    #endif

    #if PSYCH_SYSTEM == PSYCH_OSX
    // On OSX we also delay-load GLib and other GStreamer related stuff, a process
    // called "weak-linking" in Apple terminology. Use of the -weak_library XXX option
    // in our Makefile marks GLib, GStreamer et al. as weak-linked. The dynamic linker
    // will link the libraries at first use and resolve all symbols. If a symbol can't
    // get resolved, because libraries can't get loaded (missing on system, wrong permissions
    // et.c) or are too old and don't contain the symbol, the linker will define the symbol
    // as a NULL function pointer.
    //
    // Test here if GLib's g_thread_init() symbol is defined and take it as indicator
    // that GLib could get linked properly.
    if (g_thread_init == NULL) {
        // Failed: GLib and its threading support isn't installed. This means that
        // GStreamer won't work as the relevant libraries are missing on the system.
        // We silently return, skpipping the GLib init, as it is completely valid
        // for a OSX installation to not have GStreamer installed at all. On first active
        // use of a GStreamer dependent function, Screen() will detect lack of GStreamer
        // and error-out with a helpful error message to the user.
        return;
    }
    #endif
    
    // Initialize GLib's threading system early:
    g_thread_init(NULL);

    return;
}

int PsychGSGetMovieCount(void) {
	return(numMovieRecords);
}

/* Perform one context loop iteration (for bus message handling) if doWait == false,
 * or two seconds worth of iterations if doWait == true. This drives the message-bus
 * callback, so needs to be performed to get any error reporting etc.
 */
int PsychGSProcessMovieContext(GMainLoop *loop, psych_bool doWait)
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
static psych_bool PsychMoviePipelineSetState(GstElement* theMovie, GstState state, double timeoutSecs)
{
    GstState			state_pending;
    GstStateChangeReturn	rcstate;

    gst_element_set_state(theMovie, state);

    // Non-Blocking, async?
    if (timeoutSecs < 0) return(TRUE);
 
    // Wait for up to timeoutSecs for state change to complete or fail:
    rcstate = gst_element_get_state(theMovie, &state, &state_pending, (GstClockTime) (timeoutSecs * 1e9));
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
gboolean PsychMovieBusCallback(GstBus *bus, GstMessage *msg, gpointer dataptr)
{
  PsychMovieRecordType* movie = (PsychMovieRecordType*) dataptr;

  switch (GST_MESSAGE_TYPE (msg)) {
    case GST_MESSAGE_SEGMENT_DONE:
      // We usually receive segment done message instead of eos if looped playback is active and
      // the end of the stream is approaching, so we fallthrough to message eos for rewinding...
      if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: Message SEGMENT_DONE received.\n");

    case GST_MESSAGE_EOS: {
      // Rewind at end of movie if looped playback enabled:
      if ((GST_MESSAGE_TYPE (msg) == GST_MESSAGE_EOS) && (PsychPrefStateGet_Verbosity() > 5)) printf("PTB-DEBUG: Message EOS received.\n");

      if ((movie->loopflag & 0x1) && (movie->rate != 0)) {
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: End of iteration in active looped playback reached: Rewinding...\n");

        // Seek: We don't GST_SEEK_FLAG_FLUSH here, so the rewinding is smooth because we don't throw away buffers queued in the pipeline.
        if (movie->rate > 0) {
          if (!gst_element_seek(movie->theMovie, movie->rate, GST_FORMAT_TIME, GST_SEEK_FLAG_ACCURATE | GST_SEEK_FLAG_SEGMENT, GST_SEEK_TYPE_SET, 0, GST_SEEK_TYPE_NONE, GST_CLOCK_TIME_NONE)) {
            printf("PTB-DEBUG: Rewinding video in forward playback failed!\n");
          }
        }
        else {
          if (!gst_element_seek(movie->theMovie, movie->rate, GST_FORMAT_TIME, GST_SEEK_FLAG_ACCURATE | GST_SEEK_FLAG_SEGMENT, GST_SEEK_TYPE_NONE, GST_CLOCK_TIME_NONE, GST_SEEK_TYPE_END, 0)) {
            printf("PTB-DEBUG: Rewinding video in reverse playback failed!\n");
          }
        }

        // Block until seek completed, failed, or timeout of 10 seconds reached:
        gst_element_get_state(movie->theMovie, NULL, NULL, (GstClockTime) (10 * 1e9));
      }

      break;
    }

    case GST_MESSAGE_BUFFERING: {
      // Pipeline is buffering data, e.g., during network streaming playback.
      // Print some optional status info:
      gint percent = 0;
      gst_message_parse_buffering(msg, &percent);
      if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Movie '%s', buffering video data: %i percent done ...\n", movie->movieName, (int) percent);
      break;
    }

    case GST_MESSAGE_WARNING: {
      gchar  *debug;
      GError *error;

      gst_message_parse_warning(msg, &error, &debug);

      if (PsychPrefStateGet_Verbosity() > 3) {
	      printf("PTB-WARNING: GStreamer movie playback engine reports this warning:\n"
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
	      // Most common case, "File not found" error? If so, we provide a pretty-printed error message:
	      if ((error->domain == GST_RESOURCE_ERROR) && (error->code == GST_RESOURCE_ERROR_NOT_FOUND)) {
		      printf("PTB-ERROR: Could not open movie file [%s] for playback! No such moviefile with the given path and filename.\n",
			     movie->movieName);
		      printf("PTB-ERROR: The specific file URI of the missing movie was: %s.\n", movie->movieLocation);
	      }
	      else {
		      // Nope, something more special. Provide detailed GStreamer error output:
		      printf("PTB-ERROR: GStreamer movie playback engine reports this error:\n"
			     "           Error from element %s: %s\n", GST_OBJECT_NAME(msg->src), error->message);
		      printf("           Additional debug info: %s.\n\n", (debug) ? debug : "None");

		      // And some interpretation for our technically challenged users ;-):
		      if ((error->domain == GST_RESOURCE_ERROR) && (error->code != GST_RESOURCE_ERROR_NOT_FOUND)) {
			      printf("           This means that there was some problem with reading the movie file (permissions etc.).\n\n");
		      }
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

/* Video data arrived callback: Purely for documentation, because only used if oldstyle == true, that is *never*. */
static gboolean PsychHaveVideoDataCallback(GstPad *pad, GstBuffer *buffer, gpointer dataptr)
{
	unsigned int alloc_size;
	PsychMovieRecordType* movie = (PsychMovieRecordType*) dataptr;
	
	PsychLockMutex(&movie->mutex);

	if (movie->rate == 0) {
		PsychUnlockMutex(&movie->mutex);
		return(TRUE);
	}

	/* Perform onetime-init for the buffer */
	if (NULL == movie->imageBuffer) {
		// Allocate the buffer:
		alloc_size = buffer->size;
		if ((int) buffer->size < movie->width * movie->height * 4) {
			alloc_size = movie->width * movie->height * 4;
			printf("PTB-DEBUG: Overriding unsafe buffer size of %d bytes with %d bytes.\n", buffer->size, alloc_size);
		} 
		// printf("PTB-DEBUG: Allocating image buffer of %d bytes.\n", alloc_size);
		movie->imageBuffer = calloc(1, alloc_size);
	}

	// Copy new image data to our buffer:
	memcpy(movie->imageBuffer, buffer->data, buffer->size);
	movie->frameAvail++;
        // printf("PTB-DEBUG: New frame %d [size %d] %lf.\n", movie->frameAvail, buffer->size, (double) buffer->timestamp / (double) 1e9);
	
	// Fetch presentation timestamp and convert to seconds:
	movie->pts = (double) buffer->timestamp / (double) 1e9;

	PsychSignalCondition(&movie->condition);
	PsychUnlockMutex(&movie->mutex);

	return(TRUE);
}

/* Called at each end-of-stream event at end of playback: */
static void PsychEOSCallback(GstAppSink *sink, gpointer user_data)
{
	PsychMovieRecordType* movie = (PsychMovieRecordType*) user_data;

	//PsychLockMutex(&movie->mutex);
	//printf("PTB-DEBUG: Videosink reached EOS.\n");
	//PsychSignalCondition(&movie->condition);
	//PsychUnlockMutex(&movie->mutex);
	return;
}

/* Called whenever an active seek has completed or pipeline goes into pause.
 * Signals/handles arrival of preroll buffers. Used to detect/signal when
 * new videobuffers are available in non-playback mode.
 */
static GstFlowReturn PsychNewPrerollCallback(GstAppSink *sink, gpointer user_data)
{
	PsychMovieRecordType* movie = (PsychMovieRecordType*) user_data;

	PsychLockMutex(&movie->mutex);
	//printf("PTB-DEBUG: New PrerollBuffer received.\n");
	movie->preRollAvail++;
	PsychSignalCondition(&movie->condition);
	PsychUnlockMutex(&movie->mutex);

	return(GST_FLOW_OK);
}

/* Called whenever pipeline is in active playback and a new video frame arrives.
 * Used to detect/signal when new videobuffers are available in playback mode.
 */
static GstFlowReturn PsychNewBufferCallback(GstAppSink *sink, gpointer user_data)
{
	PsychMovieRecordType* movie = (PsychMovieRecordType*) user_data;

	PsychLockMutex(&movie->mutex);
	//printf("PTB-DEBUG: New Buffer received.\n");
	movie->frameAvail++;
	PsychSignalCondition(&movie->condition);
	PsychUnlockMutex(&movie->mutex);

	return(GST_FLOW_OK);
}

/* Not used by us, but needs to be defined as no-op anyway: */
static GstFlowReturn PsychNewBufferListCallback(GstAppSink *sink, gpointer user_data)
{
	PsychMovieRecordType* movie = (PsychMovieRecordType*) user_data;

	//PsychLockMutex(&movie->mutex);
	//printf("PTB-DEBUG: New Bufferlist received.\n");
	//PsychSignalCondition(&movie->condition);
	//PsychUnlockMutex(&movie->mutex);

	return(GST_FLOW_OK);
}

/* Not used by us, but needs to be defined as no-op anyway: */
static void PsychDestroyNotifyCallback(gpointer user_data)
{
	return;
}

/* This callback is called when the pipeline is about to finish playback
 * of the current movie stream. If looped playback is enabled, this needs
 * to trigger a repetition by rescheduling the movie URI for playback.
 *
 * Allows gapless playback, but doesn't work reliable on all media types.
 *
 */
static void PsychMovieAboutToFinishCB(GstElement *theMovie, gpointer user_data)
{
	PsychMovieRecordType* movie = (PsychMovieRecordType*) user_data;
	if ((movie->loopflag & 0x2) && (movie->rate != 0)) {
		g_object_set(G_OBJECT(theMovie), "uri", movie->movieLocation, NULL);
		if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: About-to-finish received: Rewinding...\n");
	}

	return;
}

/* Not used, didn't work, but left here in case we find a use for it in the future. */
static void PsychMessageErrorCB(GstBus *bus, GstMessage *msg)
{
      gchar  *debug;
      GError *error;

      gst_message_parse_error (msg, &error, &debug);
      g_free (debug);

      printf("PTB-BUSERROR: %s\n", error->message);
      g_error_free (error);
	return;
}

static GstAppSinkCallbacks videosinkCallbacks = {
    PsychEOSCallback,
    PsychNewPrerollCallback,
    PsychNewBufferCallback,
    PsychNewBufferListCallback
};

/*
 *      PsychGSCreateMovie() -- Create a movie object.
 *
 *      This function tries to open a moviefile (with or without audio/video tracks)
 *      and create an associated movie object for it.
 *
 *      win = Pointer to window record of associated onscreen window.
 *      moviename = char* with the name of the moviefile.
 *      preloadSecs = How many seconds of the movie should be preloaded/prefetched into RAM at movie open time?
 *      moviehandle = handle to the new movie.
 *      asyncFlag = As passed to 'OpenMovie'
 *      specialFlags1 = As passed to 'OpenMovie'
 */
void PsychGSCreateMovie(PsychWindowRecordType *win, const char* moviename, double preloadSecs, int* moviehandle, int asyncFlag, int specialFlags1)
{
    GstCaps                     *colorcaps;
    GstElement			*theMovie = NULL;
    GstElement			*videocodec = NULL;
    GMainLoop			*MovieContext = NULL;
    GstBus			*bus = NULL;
    GstFormat			fmt;
    GstElement			*videosink = NULL;
    gint64			length_format;
    GstPad			*pad, *peerpad;
    const GstCaps		*caps;
    GstStructure		*str;
    gint			width,height;
    gint			rate1, rate2;
    int				i, slotid;
    int				max_video_threads;
    GError			*error = NULL;
    char			movieLocation[FILENAME_MAX];
    psych_bool			trueValue = TRUE;
    char			msgerr[10000];
    char			errdesc[1000];
    psych_bool			printErrors;
    GstIterator			*it;
    psych_bool			done;
    GstPlayFlags		playflags = 0;

    // Suppress output of error-messages if moviehandle == 1000. That means we
    // run in our own Posix-Thread, not in the Matlab-Thread. Printing via Matlabs
    // printing facilities would likely cause a terrible crash.
    printErrors = (*moviehandle == -1000) ? FALSE : TRUE;

    // Gapless playback requested? Normally *moviehandle is == -1, so a positive
    // handle requests this mode and defines the actual handle of the movie to use:
    if (*moviehandle >= 0) {
	// Queueing a new moviename of a movie to play next: This only works
	// for already opened/created movies whose pipeline is at least in
	// READY state, better PAUSED or PLAYING. Validate preconditions:

	// Valid handle for existing movie?
	if (*moviehandle < 0 || *moviehandle >= PSYCH_MAX_MOVIES) {
	    PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
	}
        
	// Fetch references to objects we need:
	theMovie = movieRecordBANK[*moviehandle].theMovie;    
	if (theMovie == NULL) {
	    PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
	}

	// Ok, this means we have a handle to an existing, fully operational
	// playback pipeline. Convert moviename to a valid URL and queue it:

	// Create name-string for moviename: If an URI qualifier is at the beginning,
	// we're fine and just pass the URI as-is. Otherwise we add the file:// URI prefix.
	if (strstr(moviename, "://") || ((strstr(moviename, "v4l") == moviename) && strstr(moviename, "//"))) {
	    snprintf(movieLocation, sizeof(movieLocation)-1, "%s", moviename);
	} else {
	    snprintf(movieLocation, sizeof(movieLocation)-1, "file:///%s", moviename);
	}

	strncpy(movieRecordBANK[*moviehandle].movieLocation, movieLocation, FILENAME_MAX);
	strncpy(movieRecordBANK[*moviehandle].movieName, moviename, FILENAME_MAX);

	// Assign name of movie to play to pipeline. If the pipeline is not in playing
	// state, this will switch to the specified movieLocation immediately. If it
	// is playing, it will switch to it at the end of the current playback iteration:
	g_object_set(G_OBJECT(theMovie), "uri", movieLocation, NULL);

	// Ready.
	return;
    }
    
    // Set movie handle to "failed" initially:
    *moviehandle = -1;

    // We start GStreamer only on first invocation.
    if (firsttime) {        
        // Initialize GStreamer: The routine is defined in PsychVideoCaptureSupportGStreamer.c
        PsychGSCheckInit("movie playback");

        // Enable use of YUV textures for movie playback on supported GPUs if environment variable
        // is defined. YUV mode defaults to "off", because as of 1st December 2011, at least H264
        // to YUV decoding was much slower than bog standard decoding to RGBA8 -- much to my surprise.
        if (getenv("PSYCHTOOLBOX_USE_YUV_MOVIEDECODING") || (specialFlags1 & 1)) useYUVDecode = TRUE;

        firsttime = FALSE;
    }

    if (win && !PsychIsOnscreenWindow(win)) {
        if (printErrors) PsychErrorExitMsg(PsychError_user, "Provided windowPtr is not an onscreen window."); else return;
    }

    if (NULL == moviename) {
        if (printErrors) PsychErrorExitMsg(PsychError_internal, "NULL-Ptr instead of moviename passed!"); else return;
    }

    if (numMovieRecords >= PSYCH_MAX_MOVIES) {
        *moviehandle = -2;
        if (printErrors) PsychErrorExitMsg(PsychError_user, "Allowed maximum number of simultaneously open movies exceeded!"); else return;
    }

    // Search first free slot in movieRecordBANK:
    for (i=0; (i < PSYCH_MAX_MOVIES) && (movieRecordBANK[i].theMovie); i++);
    if (i>=PSYCH_MAX_MOVIES) {
        *moviehandle = -2;
        if (printErrors) PsychErrorExitMsg(PsychError_user, "Allowed maximum number of simultaneously open movies exceeded!"); else return;
    }

    // Slot slotid will contain the movie record for our new movie object:
    slotid=i;

    // Zero-out new record in moviebank:
    memset(&movieRecordBANK[slotid], 0, sizeof(PsychMovieRecordType));
    
    // Create name-string for moviename: If an URI qualifier is at the beginning,
    // we're fine and just pass the URI as-is. Otherwise we add the file:// URI prefix.
    if (strstr(moviename, "://") || ((strstr(moviename, "v4l") == moviename) && strstr(moviename, "//"))) {
	snprintf(movieLocation, sizeof(movieLocation)-1, "%s", moviename);
    } else {
	snprintf(movieLocation, sizeof(movieLocation)-1, "file:///%s", moviename);
    }
    strncpy(movieRecordBANK[slotid].movieLocation, movieLocation, FILENAME_MAX);
    strncpy(movieRecordBANK[slotid].movieName, moviename, FILENAME_MAX);

    // Create movie playback pipeline:
    if (TRUE) {
	// Use playbin2:
	theMovie = gst_element_factory_make ("playbin2", "ptbmovieplaybackpipeline");

	// Assign name of movie to play:
	g_object_set(G_OBJECT(theMovie), "uri", movieLocation, NULL);

	// Default flags for playbin: Decode video and deinterlace it if needed:
	playflags = GST_PLAY_FLAG_VIDEO  | GST_PLAY_FLAG_DEINTERLACE;

	// Decode and play audio by default, with software audio volume control, unless specialFlags setting 2 enabled:
	if (!(specialFlags1 & 2)) playflags |= GST_PLAY_FLAG_AUDIO | GST_PLAY_FLAG_SOFT_VOLUME;

	// Enable network buffering for network videos of at least 10 seconds, or preloadSecs seconds,
	// whatever is bigger.

	// Setup without any buffering and caching (aka preloadSecs == 0) requested?
	// Note: For now treat the default preloadSecs value 1 as a zero -> No buffering and caching.
	// Why? Because the usefulness of this extra setup is not yet proven and the specific choice
	// of buffering parameters may need a bit of tuning. We don't want to cause regressions in
	// performance of existing scripts, so we stick to the GStreamer default buffering behaviour
	// until more time has been spent tuning & testing this setup code.
	if ((preloadSecs != 0) && (preloadSecs != 1)) {
	    // No: Use internal buffering/caching [BUFFERING] of demultiplexed/parsed data, e.g., for fast
	    // recycling during looped video playback, random access out-of-order frame fetching, fast
	    // seeking and reverse playback:
	    playflags |= GST_PLAY_FLAG_BUFFERING;

	    // Ok, this is ugly: Some movie formats, when streamed from the internet, need progressive
	    // download buffering to work without problems, whereas other formats will cause problems
	    // with progressive download buffering. So far we know that some .mov Quicktime movies, e.g.,
	    // Apple's commercials need it, whereas some .webm movies choke on it. Let's be optimistic
	    // and assume it works with everything except .webm. Also provide secret cheat code == -2
	    // to override the blacklisting of .webm to allow for further experiments:
	    if ((preloadSecs == -2) || (!strstr(moviename, ".webm"))) {
		// Want some local progressive download buffering [DOWNLOAD] for network video streams,
		// as temporary file on local filesystem:
		playflags |= GST_PLAY_FLAG_DOWNLOAD;
	    }

	    // Undo our cheat-code if used: Map to 10 seconds preload time:
	    if (preloadSecs == -2) preloadSecs = 10;

	    // Setting maximum size of internal RAM ringbuffer supported? (since v0.10.31)
	    if (g_object_class_find_property(G_OBJECT_GET_CLASS(theMovie), "ring-buffer-max-size")) {
		// Supported. The ringbuffer is disabled by default, we enable it with a certain maximum
		// size in bytes. For preloadSecs == -1, aka "unlimited buffering", we set it to its
		// allowable maximum of G_MAXUINT == 4 GB. For a given finite preloadSecs we have
		// to set something reasonable. Set it to preloadSecs buffer duration (in seconds) multiplied
		// by some assumed maximum datarate in bytes/second. We use 4e6 bytes, which is roughly
		// 4 MB/sec. Why? This is a generously padded value, assuming a max. fps of 60 Hz, max.
		// resolution 1920x1080p HD video + HD audio. Numbers are based on the bit rates of
		// a HD movie trailer (Warner Brothers "I am Legend" public HD movie trailer), which has
		// 7887 kbits/s for 1920x816 H264/AVC progessive scan video at 24 fps and 258 kbits/s for
		// MPEG-4 AAC audio in Surround 5.1 format with 48 kHz sampling rate. This upscaled to
		// research use and padded should give a good value for our purpose. Also at a default
		// preloadSecs value of 1 second, this wastes at most 4 MB for buffering - a safe default:
		g_object_set(G_OBJECT(theMovie), "ring-buffer-max-size", ((preloadSecs == -1) ? G_MAXUINT : (guint64) (preloadSecs * (double) 4e6)), NULL);
		if (PsychPrefStateGet_Verbosity() > 4) {
		    printf("PTB-INFO: Playback for movie %i will use adapted RAM ring-buffer-max-size of %f MB.\n", slotid,
			   (float) (((double) ((preloadSecs == -1) ? G_MAXUINT : preloadSecs * (double) 4e6)) / 1024.0 / 1024.0));
		}
	    }

	    // Setting of maximum buffer duration for network video stream playback:
	    if (preloadSecs == -1) {
		// "Unlimited" - Set maximum buffering size to G_MAXINT == 2 GB.
		g_object_set(G_OBJECT(theMovie), "buffer-size", (gint) G_MAXINT, NULL);
	    }
	    else {
		// Limited - Set maximum buffer-duration to preloadSecs, the playbin will derive
		// a proper maximum buffering size from duration and streaming bitrate:
		g_object_set(G_OBJECT(theMovie), "buffer-duration", (gint64) (preloadSecs * (double) 1e9), NULL);
	    }

	    if (PsychPrefStateGet_Verbosity() > 4) {
		printf("PTB-INFO: Playback for movie %i will use RAM buffering. Additional prebuffering for network streams is\n", slotid);
		printf("PTB-INFO: limited to %f %s.\n", (preloadSecs == -1) ? 2 : preloadSecs, (preloadSecs == -1) ? "GB" : "seconds");
		if (playflags & GST_PLAY_FLAG_DOWNLOAD) printf("PTB-INFO: Network video streams will be additionally cached to the filesystem.\n");
	    }

	    // All in all, we can end up with up to 6*x GB RAM and 6 GB disc consumption for the "unlimited" setting,
	    // about 4*x MB RAM and 4 MB disc consumption for the default setting of 1, and preloadSecs multiples of
	    // that for a given value. x is an unkown factor, depending on which internal plugins maintain ringbuffers,
	    // but assume x somewhere between 1 and maybe 4.
	}
	else {
	    if (PsychPrefStateGet_Verbosity() > 4) {
		printf("PTB-INFO: Playback for movie %i will not use additional buffering or caching due to 'preloadSecs' setting %i.\n", slotid, (int) preloadSecs);
	    }
	}

	// Setup final playback control flags:
	g_object_set(G_OBJECT(theMovie), "flags", playflags , NULL);

	// Connect callback to about-to-finish signal: Signal is emitted as soon as
	// end of current playback iteration is approaching. The callback checks if
	// looped playback is requested. If so, it schedules a new playback iteration.
	g_signal_connect(G_OBJECT(theMovie), "about-to-finish", G_CALLBACK(PsychMovieAboutToFinishCB), &(movieRecordBANK[slotid]));
    }
    else {
	// Self-Assembled pipeline: Does not work for some not yet investigated reason,
	// but is not needed anyway, so we disable it and just leave it for documentation,
	// in case it will be needed in the future:
	sprintf(movieLocation, "filesrc location='%s' ! qtdemux ! queue ! ffdec_h264 ! ffmpegcolorspace ! appsink name=ptbsink0", moviename);
	theMovie = gst_parse_launch((const gchar*) movieLocation, NULL);
	videosink = gst_bin_get_by_name(GST_BIN(theMovie), "ptbsink0");
	printf("LAUNCHLINE[%p]: %s\n", videosink, movieLocation);
    }

    // Assign message context, message bus and message callback for
    // the pipeline to report events and state changes, errors etc.:    
    MovieContext = g_main_loop_new (NULL, FALSE);
    movieRecordBANK[slotid].MovieContext = MovieContext;
    bus = gst_pipeline_get_bus(GST_PIPELINE(theMovie));
    // Didn't work: g_signal_connect (G_OBJECT(bus), "message::error", G_CALLBACK(PsychMessageErrorCB), NULL);
    //              g_signal_connect (G_OBJECT(bus), "message::warning", G_CALLBACK(PsychMessageErrorCB), NULL);
    gst_bus_add_watch(bus, PsychMovieBusCallback, &(movieRecordBANK[slotid]));
    gst_object_unref(bus);

    // Assign a fakesink named "ptbsink0" as destination video-sink for
    // all video content. This allows us to get hold of the video frame buffers for
    // converting them into PTB OpenGL textures:
    if (!videosink) videosink = gst_element_factory_make ("appsink", "ptbsink0");
    if (!videosink) {
	printf("PTB-ERROR: Failed to create video-sink appsink ptbsink!\n");
	PsychGSProcessMovieContext(movieRecordBANK[slotid].MovieContext, TRUE);
	PsychErrorExitMsg(PsychError_system, "Opening the movie failed. Reason hopefully given above.");
    }

    movieRecordBANK[slotid].videosink = videosink;

    // Our OpenGL texture creation routine usually needs GL_BGRA8 data in G_UNSIGNED_8_8_8_8_REV
    // format, but the pipeline usually delivers YUV data in planar format. Therefore
    // need to perform colorspace/colorformat conversion. We build a little videobin
    // which consists of a ffmpegcolorspace converter plugin connected to our appsink
    // plugin which will deliver video data to us for conversion into textures.
    // The "sink" pad of the converter plugin is connected as the "sink" pad of our
    // videobin, and the videobin is connected to the video-sink output of the pipeline,
    // thereby receiving decoded video data. We place a videocaps filter inbetween the
    // converter and the appsink to enforce a color format conversion to the "colorcaps"
    // we need. colorcaps define the needed data format for efficient conversion into
    // a RGBA8 texture. Some GPU + driver combos do support direct handling of UYVU YCrCb
    // data as textures. If we are on such a GPU we request yuv UYVU data and upload it
    // directly in this format to the GPU. This more efficient both for GStreamers decode
    // pipeline, and the later Videobuffer -> OpenGL texture conversion:
    if (win && (win->gfxcaps & kPsychGfxCapUYVYTexture) && useYUVDecode) {
	// GPU supports handling and decoding of UYVY type yuv textures: We use these,
	// as they are more efficient to decode and handle by typical video codecs:
	colorcaps = gst_caps_new_simple ( "video/x-raw-yuv",
					  "format", GST_TYPE_FOURCC, GST_MAKE_FOURCC('U', 'Y', 'V', 'Y'),
					  NULL);
	if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Movie playback for movie %i will use UYVY YCrCb textures for optimized decode and rendering.\n", slotid);
    } else {
	// GPU does not support yuv textures. Need to go brute-force and convert
	// video into RGBA8 format:
	colorcaps = gst_caps_new_simple ( "video/x-raw-rgb",
					  "bpp", G_TYPE_INT, 32,
					  "depth", G_TYPE_INT, 32,
					  "alpha_mask", G_TYPE_INT, 0x000000FF,
					  "red_mask", G_TYPE_INT,   0x0000FF00,
					  "green_mask", G_TYPE_INT, 0x00FF0000,
					  "blue_mask", G_TYPE_INT,  0xFF000000,
					  NULL);
	if ((PsychPrefStateGet_Verbosity() > 3) && useYUVDecode) printf("PTB-INFO: Movie playback for movie %i will use RGBA8 textures due to lack of YUV texture support on GPU.\n", slotid);
	if ((PsychPrefStateGet_Verbosity() > 3) && !useYUVDecode) printf("PTB-INFO: Movie playback for movie %i will use RGBA8 textures.\n", slotid);
    }

    /*
    // Old style method: Only left here for documentation to show how one can create
    // video sub-pipelines via bin's and connect them to each other via ghostpads: 

    GstElement *videobin = gst_bin_new ("video_output_bin");
    GstElement *videocon = gst_element_factory_make ("ffmpegcolorspace", "color_converter");
    gst_bin_add_many(GST_BIN(videobin), videocon, videosink, NULL);

    GstPad *ghostpad = gst_ghost_pad_new("Video_Ghostsink", gst_element_get_pad(videocon, "sink"));
    gst_element_add_pad(videobin, ghostpad);

    gst_element_link_filtered(videocon, videosink, colorcaps);

    // Assign our special videobin as video-sink of the pipeline:
    g_object_set(G_OBJECT(theMovie), "video-sink", videobin, NULL);
    */

    // New style method: Leaves the freedom of choice of color converter (if any)
    // to the auto-plugger.

    // Assign 'colorcaps' as caps to our videosink. This marks the videosink so
    // that it can only receive video image data in the format defined by colorcaps,
    // i.e., a format that is easy to consume for OpenGL's texture creation on std.
    // gpu's. It is the job of the video pipeline's autoplugger to plug in proper
    // color & format conversion plugins to satisfy videosink's needs.
    gst_app_sink_set_caps(GST_APP_SINK(videosink), colorcaps);

    // Assign our special appsink 'videosink' as video-sink of the pipeline:
    g_object_set(G_OBJECT(theMovie), "video-sink", videosink, NULL);
    gst_caps_unref(colorcaps);

    // Get the pad from the final sink for probing width x height of movie frames and nominal framerate of movie:
    pad = gst_element_get_pad(videosink, "sink");

    PsychGSProcessMovieContext(movieRecordBANK[slotid].MovieContext, FALSE);

    // Should we preroll / preload?
    // MK: Actually, *always* preroll via transition to PAUSED. The other path
    // is seldomly used, therefore not well tested, bound to cause regressions on
    // future code changes. Also just readying the pipeline causes problems with
    // movie duration queries and other property queries, so probably not worth
    // the trouble...
    // The preloadSecs parameter is used instead above to control buffering behaviour.
    if (TRUE || (preloadSecs > 0) || (preloadSecs == -1)) {
	// Preload / Preroll the pipeline:
	if (!PsychMoviePipelineSetState(theMovie, GST_STATE_PAUSED, 30.0)) {
		PsychGSProcessMovieContext(movieRecordBANK[slotid].MovieContext, TRUE);
		PsychErrorExitMsg(PsychError_user, "In OpenMovie: Opening the movie failed I. Reason given above.");
	}
    } else {
	// Ready the pipeline:
	if (!PsychMoviePipelineSetState(theMovie, GST_STATE_READY, 30.0)) {
		PsychGSProcessMovieContext(movieRecordBANK[slotid].MovieContext, TRUE);
		PsychErrorExitMsg(PsychError_user, "In OpenMovie: Opening the movie failed II. Reason given above.");
	}
    }

    // Check if a multi-threaded decoder from is used: If so, set number of processing
    // threads to use: By default many codecs would only use one single thread on any system,
    // even if they are multi-threading capable.
    it = gst_bin_iterate_recurse(GST_BIN(theMovie));
    done = FALSE;
    videocodec = NULL;

    while (!done) {
	switch (gst_iterator_next(it, (void**) &videocodec)) {
	    case GST_ITERATOR_OK:
		if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In pipeline: Child element name: %s\n", (const char*) gst_object_get_name(GST_OBJECT(videocodec)));
		//if (strstr((const char*) gst_object_get_name(GST_OBJECT(videocodec)), "h264")) {
		if (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "max-threads")) {
		    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Found video decoder element %s.\n", (const char*) gst_object_get_name(GST_OBJECT(videocodec)));
		    done = TRUE;
		} else {
		    gst_object_unref(videocodec);
		    videocodec = NULL;
		}
	    break;

	    case GST_ITERATOR_RESYNC:
	        gst_iterator_resync(it);
	    break;

	    case GST_ITERATOR_DONE:
		done = TRUE;
	    break;
       }
    }

    gst_iterator_free(it);
    it = NULL;

    if (videocodec && (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "max-threads"))) {
        max_video_threads = 1;
	g_object_get(G_OBJECT(videocodec), "max-threads", &max_video_threads, NULL);
	if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Movie playback for movie %i uses video decoder with a default maximum number of %i processing threads.\n", slotid, max_video_threads);

	// Set max_threads to 0: This means to auto-detect the optimal number of threads.
	if (getenv("PSYCHTOOLBOX_MAX_VIDEODECODER_THREADS")) {
	    max_video_threads = atoi(getenv("PSYCHTOOLBOX_MAX_VIDEODECODER_THREADS"));
	    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Setting video decoder to use a maximum of %i processing threads.\n", max_video_threads);
	} else {
	    max_video_threads = 0;
	    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Setting video decoder to use auto-selected optimal number of processing threads.\n");
	}

	// Ready the video codec, so a new max thread count can be set:
	if (!PsychMoviePipelineSetState(videocodec, GST_STATE_READY, 30.0)) {
		PsychGSProcessMovieContext(movieRecordBANK[slotid].MovieContext, TRUE);
		PsychErrorExitMsg(PsychError_user, "In OpenMovie: Opening the movie failed III. Reason given above.");
	}    

	g_object_set(G_OBJECT(videocodec), "max-threads", max_video_threads, NULL);

	// Pause the video codec, so the new max thread count is accepted:
	if (!PsychMoviePipelineSetState(videocodec, GST_STATE_PAUSED, 30.0)) {
		PsychGSProcessMovieContext(movieRecordBANK[slotid].MovieContext, TRUE);
		PsychErrorExitMsg(PsychError_user, "In OpenMovie: Opening the movie failed IV. Reason given above.");
	}

	g_object_get(G_OBJECT(videocodec), "max-threads", &max_video_threads, NULL);
	if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Movie playback for movie %i uses video decoder with a current maximum number of %i processing threads.\n", slotid, max_video_threads);
    }

    // Release reference to videocodec:
    if (videocodec) gst_object_unref(videocodec);
    videocodec = NULL;

    // Query number of available video and audio tracks in movie:
    g_object_get (G_OBJECT(theMovie),
               "n-video", &movieRecordBANK[slotid].nrVideoTracks,
               "n-audio", &movieRecordBANK[slotid].nrAudioTracks,
                NULL);

    // We need a valid onscreen window handle for real video playback:
    if ((NULL == win) && (movieRecordBANK[slotid].nrVideoTracks > 0)) {
        if (printErrors) PsychErrorExitMsg(PsychError_user, "No windowPtr to an onscreen window provided. Must do so for movies with videotrack!"); else return;
    }
 
    PsychGSProcessMovieContext(movieRecordBANK[slotid].MovieContext, FALSE);

    PsychInitMutex(&movieRecordBANK[slotid].mutex);
    PsychInitCondition(&movieRecordBANK[slotid].condition, NULL);

    if (oldstyle) {
	// Install the probe callback for reception of video frames from engine at the sink-pad itself:
	gst_pad_add_buffer_probe(pad, G_CALLBACK(PsychHaveVideoDataCallback), &(movieRecordBANK[slotid]));
    } else {
	// Install callbacks used by the videosink (appsink) to announce various events:
	gst_app_sink_set_callbacks(GST_APP_SINK(videosink), &videosinkCallbacks, &(movieRecordBANK[slotid]), PsychDestroyNotifyCallback);
    }

    // Assign harmless initial settings for fps and frame size:
    rate1 = 0;
    rate2 = 1;
    width = height = 0;

    // Videotrack available?
    if (movieRecordBANK[slotid].nrVideoTracks > 0) {
	// Yes: Query size and framerate of movie:
	peerpad = gst_pad_get_peer(pad);
	caps=gst_pad_get_negotiated_caps(peerpad);
	if (caps) {
		str=gst_caps_get_structure(caps,0);

		/* Get some data about the frame */
		rate1 = 1; rate2 = 1;
		gst_structure_get_fraction(str, "pixel-aspect-ratio", &rate1, &rate2);
		movieRecordBANK[slotid].aspectRatio = (double) rate1 / (double) rate2;
		gst_structure_get_int(str,"width",&width);
		gst_structure_get_int(str,"height",&height);
		rate1 = 0; rate2 = 1;
		gst_structure_get_fraction(str, "framerate", &rate1, &rate2);

	 } else {
		printf("PTB-DEBUG: No frame info available after preroll.\n");	
	 }
    }

    if (strstr(moviename, "v4l2:")) {
	// Special case: The "movie" is actually a video4linux2 live source.
	// Need to make parameters up for now, so it to work as "movie":
	rate1 = 30; width = 640; height = 480;
	movieRecordBANK[slotid].nrVideoTracks = 1;

	// Uglyness at its best ;-)
	if (strstr(moviename, "320")) { width = 320; height = 240; };
    }

    // Release the pad:
    gst_object_unref(pad);

    // Assign new record in moviebank:
    movieRecordBANK[slotid].theMovie = theMovie;
    movieRecordBANK[slotid].loopflag = 0;
    movieRecordBANK[slotid].frameAvail = 0;
    movieRecordBANK[slotid].imageBuffer = NULL;
    movieRecordBANK[slotid].startPending = 0;
    movieRecordBANK[slotid].endOfFetch = 0;

    *moviehandle = slotid;

    // Increase counter:
    numMovieRecords++;

    // Compute basic movie properties - Duration and fps as well as image size:
    
    // Retrieve duration in seconds:
    fmt = GST_FORMAT_TIME;
    if (gst_element_query_duration(theMovie, &fmt, &length_format)) {
	// This returns nsecs, so convert to seconds:
    	movieRecordBANK[slotid].movieduration = (double) length_format / (double) 1e9;
	//printf("PTB-DEBUG: Duration of movie %i [%s] is %lf seconds.\n", slotid, moviename, movieRecordBANK[slotid].movieduration);
    } else {
	movieRecordBANK[slotid].movieduration = DBL_MAX;
	printf("PTB-WARNING: Could not query duration of movie %i [%s] in seconds. Returning infinity.\n", slotid, moviename);
    }

    // Assign expected framerate, assuming a linear spacing between frames:
    movieRecordBANK[slotid].fps = (double) rate1 / (double) rate2;
    //printf("PTB-DEBUG: Framerate fps of movie %i [%s] is %lf fps.\n", slotid, moviename, movieRecordBANK[slotid].fps);

    // Drop frames if callback can't pull buffers fast enough, unless ascynFlags & 4 is set:
    // This together with a max queue lengths of 1 allows to
    // maintain audio-video sync by framedropping if needed.
    gst_app_sink_set_drop(GST_APP_SINK(videosink), (asyncFlag & 4) ? FALSE : TRUE);

	// Buffering of decoded video frames requested?
	if (asyncFlag & 4) {
		// Yes: If a specific preloadSecs and a valid fps playback framerate is available, we
		// set the maximum buffer capacity to the number of frames corresponding to the given 'preloadSecs'.
		// Otherwise we set it to zero, which means "unlimited capacity", ie., until RAM full:
		gst_app_sink_set_max_buffers(GST_APP_SINK(videosink), ((movieRecordBANK[slotid].fps > 0) && (preloadSecs >= 0)) ? ((int) (movieRecordBANK[slotid].fps * preloadSecs) + 1) : 0);
	}
	else {
		// No: Only allow one queued buffer before dropping, to avoid optimal audio-video sync:
		gst_app_sink_set_max_buffers(GST_APP_SINK(videosink), 1);
	}
	
    // Compute framecount from fps and duration:
    movieRecordBANK[slotid].nrframes = (int)(movieRecordBANK[slotid].fps * movieRecordBANK[slotid].movieduration + 0.5);
    //printf("PTB-DEBUG: Number of frames in movie %i [%s] is %i.\n", slotid, moviename, movieRecordBANK[slotid].nrframes);

    // Define size of images in movie:
    movieRecordBANK[slotid].width = width;
    movieRecordBANK[slotid].height = height;

    // Ready to rock!
    return;
}

/*
 *  PsychGSGetMovieInfos() - Return basic information about a movie.
 *
 *  framecount = Total number of video frames in the movie, determined by counting.
 *  durationsecs = Total playback duration of the movie, in seconds.
 *  framerate = Estimated video playback framerate in frames per second (fps).
 *  width = Width of movie images in pixels.
 *  height = Height of movie images in pixels.
 *  nrdroppedframes = Total count of videoframes that had to be dropped during last movie playback,
 *                    in order to keep the movie synced with the realtime clock.
 */
void PsychGSGetMovieInfos(int moviehandle, int* width, int* height, int* framecount, double* durationsecs, double* framerate, int* nrdroppedframes, double* aspectRatio)
{
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
    
    if (movieRecordBANK[moviehandle].theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }

    if (framecount) *framecount = movieRecordBANK[moviehandle].nrframes;
    if (durationsecs) *durationsecs = movieRecordBANK[moviehandle].movieduration;
    if (framerate) *framerate = movieRecordBANK[moviehandle].fps;
    if (nrdroppedframes) *nrdroppedframes = movieRecordBANK[moviehandle].nr_droppedframes;
    if (width) *width = movieRecordBANK[moviehandle].width; 
    if (height) *height = movieRecordBANK[moviehandle].height; 
    if (aspectRatio) *aspectRatio = movieRecordBANK[moviehandle].aspectRatio;

    return;
}

/*
 *  PsychGSDeleteMovie() -- Delete a movie object and release all associated ressources.
 */
void PsychGSDeleteMovie(int moviehandle)
{
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
    
    if (movieRecordBANK[moviehandle].theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }
        
    // Stop movie playback immediately:
    PsychMoviePipelineSetState(movieRecordBANK[moviehandle].theMovie, GST_STATE_NULL, 20.0);

    // Delete movieobject for this handle:
    gst_object_unref(GST_OBJECT(movieRecordBANK[moviehandle].theMovie));
    movieRecordBANK[moviehandle].theMovie=NULL;

    // Delete visual context for this movie:
    movieRecordBANK[moviehandle].MovieContext = NULL;

    PsychDestroyMutex(&movieRecordBANK[moviehandle].mutex);
    PsychDestroyCondition(&movieRecordBANK[moviehandle].condition);

    free(movieRecordBANK[moviehandle].imageBuffer);
    movieRecordBANK[moviehandle].imageBuffer = NULL;
    movieRecordBANK[moviehandle].videosink = NULL;

	// Recycled texture in texture cache?
    if (movieRecordBANK[moviehandle].cached_texture > 0) {
		// Yes. Release it.
		glDeleteTextures(1, &(movieRecordBANK[moviehandle].cached_texture));
		movieRecordBANK[moviehandle].cached_texture = 0;
	}

    // Decrease counter:
    if (numMovieRecords>0) numMovieRecords--;
        
    return;
}

/*
 *  PsychGSDeleteAllMovies() -- Delete all movie objects and release all associated ressources.
 */
void PsychGSDeleteAllMovies(void)
{
    int i;
    for (i=0; i<PSYCH_MAX_MOVIES; i++) {
        if (movieRecordBANK[i].theMovie) PsychGSDeleteMovie(i);
    }
    return;
}

/*
 *  PsychGSGetTextureFromMovie() -- Create an OpenGL texture map from a specific videoframe from given movie object.
 *
 *  win = Window pointer of onscreen window for which a OpenGL texture should be created.
 *  moviehandle = Handle to the movie object.
 *  checkForImage = 0 == Retrieve the image, blocking until error timeout if necessary.
 *                  1 == Check for new image in polling fashion.
 *                  2 == Check for new image in blocking fashion. Wait up to 5 seconds blocking for a new frame.
 *  timeindex = When not in playback mode, this allows specification of a requested frame by presentation time.
 *              If set to -1, or if in realtime playback mode, this parameter is ignored and the next video frame is returned.
 *  out_texture = Pointer to the Psychtoolbox texture-record where the new texture should be stored.
 *  presentation_timestamp = A ptr to a double variable, where the presentation timestamp of the returned frame should be stored.
 *
 *  Returns true (1) on success, false (0) if no new image available, -1 if no new image available and there won't be any in future.
 */
int PsychGSGetTextureFromMovie(PsychWindowRecordType *win, int moviehandle, int checkForImage, double timeindex,
			     PsychWindowRecordType *out_texture, double *presentation_timestamp)
{
    GstElement			*theMovie;
    unsigned int		failcount=0;
    double			rate;
    double			targetdelta, realdelta, frames;
    GstBuffer                   *videoBuffer = NULL;
    gint64		        bufferIndex;
    double                      deltaT = 0;
    GstEvent                    *event;
    static double               tStart = 0;
    double                      tNow;
    double                      preT, postT;

    if (!PsychIsOnscreenWindow(win)) {
        PsychErrorExitMsg(PsychError_user, "Need onscreen window ptr!!!");
    }
    
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided.");
    }
    
    if ((timeindex!=-1) && (timeindex < 0 || timeindex >= 100000.0)) {
        PsychErrorExitMsg(PsychError_user, "Invalid timeindex provided.");
    }
    
    if (NULL == out_texture && !checkForImage) {
        PsychErrorExitMsg(PsychError_internal, "NULL-Ptr instead of out_texture ptr passed!!!");
    }
    
    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle.");
    }

    // Deferred start of movie playback requested? This so if movie is supposed to be
    // actively playing (rate != 0) and the startPending flag marks a pending deferred start:
    if ((movieRecordBANK[moviehandle].rate != 0) && movieRecordBANK[moviehandle].startPending) {
        // Deferred start: Reset flag, start pipeline with a max timeout of 1 second:
        movieRecordBANK[moviehandle].startPending = 0;
        PsychMoviePipelineSetState(theMovie, GST_STATE_PLAYING, 1);
        // This point is reached after either the pipeline is fully started, or the
        // timeout has elapsed. In the latter case, a GST_STATE_CHANGE_ASYNC message
        // is printed and start of pipeline continues asynchronously. No big deal for
        // us, as we'll simply block in the rest of the texture fetch (checkForImage) path
        // until the first frame is ready and audio playback has started. The main purpose
        // of setting a reasonable timeout above is to avoid cluttering the console with
        // status messages (timeout big enough for common case) but allow user to interrupt
        // ops that take too long (timeout small enough to avoid long user-perceived exec-hangs).
        // 1 Second is used to cater to the common case of playing files from disc, but coping
        // with multi-second delays for network streaming (buffering delays in preroll).
    }

    // Allow context task to do its internal bookkeeping and cleanup work:
    PsychGSProcessMovieContext(movieRecordBANK[moviehandle].MovieContext, FALSE);

    // If this is a pure audio "movie" with no video tracks, we always return failed,
    // as those certainly don't have movie frames associated.
    if (movieRecordBANK[moviehandle].nrVideoTracks == 0) return((checkForImage) ? -1 : FALSE);

    // Get current playback rate:
    rate = movieRecordBANK[moviehandle].rate;

    // Is movie actively playing (automatic async playback, possibly with synced sound)?
    // If so, then we ignore the 'timeindex' parameter, because the automatic playback
    // process determines which frames should be delivered to PTB when. This function will
    // simply wait or poll for arrival/presence of a new frame that hasn't been fetched
    // in previous calls.
    if (0 == rate) {
        // Movie playback inactive. We are in "manual" mode: No automatic async playback,
        // no synced audio output. The user just wants to manually fetch movie frames into
        // textures for manual playback in a standard Matlab-loop.
		
		// First pass - checking for new image?
		if (checkForImage) {
			// Image for specific point in time requested?
			if (timeindex >= 0) {
				// Yes. We try to retrieve the next possible image for requested timeindex.
				// Seek to target timeindex:
				PsychGSSetMovieTimeIndex(moviehandle, timeindex, FALSE);
			}
			// Check for frame availability happens down there in the shared check code...
		}
    }
	
    // Should we just check for new image? If so, just return availability status:
    if (checkForImage) {
		// Take reference timestamps of fetch start:
		if (tStart == 0) PsychGetAdjustedPrecisionTimerSeconds(&tStart);
		PsychLockMutex(&movieRecordBANK[moviehandle].mutex);
		
		if ((((0 != rate) && movieRecordBANK[moviehandle].frameAvail) || ((0 == rate) && movieRecordBANK[moviehandle].preRollAvail)) &&
			!gst_app_sink_is_eos(GST_APP_SINK(movieRecordBANK[moviehandle].videosink))) {
			// New frame available. Unlock and report success:
			//printf("PTB-DEBUG: NEW FRAME %d\n", movieRecordBANK[moviehandle].frameAvail);
			PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
			return(TRUE);
		}
		
		// Is this the special case of a movie without video, but only sound? In that case
		// we always return a 'false' because there ain't no image to return. We check this
		// indirectly - If the imageBuffer is NULL then the video callback hasn't been called.
		if (oldstyle && (NULL == movieRecordBANK[moviehandle].imageBuffer)) {
			PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
			return(FALSE);
		}
		
		// None available. Any chance there will be one in the future?
        if (((rate != 0) && gst_app_sink_is_eos(GST_APP_SINK(movieRecordBANK[moviehandle].videosink)) && (movieRecordBANK[moviehandle].loopflag == 0)) ||
            ((rate == 0) && (movieRecordBANK[moviehandle].endOfFetch))) {
			// No new frame available and there won't be any in the future, because this is a non-looping
			// movie that has reached its end.
			movieRecordBANK[moviehandle].endOfFetch = 0;
			PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
			return(-1);
        }
        else {
			// No new frame available yet:
			PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
			//printf("PTB-DEBUG: NO NEW FRAME\n");
			
			// In the polling check, we return with statue "no new frame yet" aka false:
			if (checkForImage < 2) return(FALSE);
			
			// Otherwise (blocking/waiting check) we fall-through the wait code below...
        }
    }
	
    // If we reach this point, then an image fetch is requested. If no new data
    // is available we shall block:

    PsychLockMutex(&movieRecordBANK[moviehandle].mutex);
    // printf("PTB-DEBUG: Blocking fetch start %d\n", movieRecordBANK[moviehandle].frameAvail);

    if (((0 != rate) && !movieRecordBANK[moviehandle].frameAvail) ||
		((0 == rate) && !movieRecordBANK[moviehandle].preRollAvail)) {
		// No new frame available. Perform a blocking wait with timeout of 0.5 seconds:
		PsychTimedWaitCondition(&movieRecordBANK[moviehandle].condition, &movieRecordBANK[moviehandle].mutex, 0.5);
		
		// Recheck:
		if (((0 != rate) && !movieRecordBANK[moviehandle].frameAvail) ||
			((0 == rate) && !movieRecordBANK[moviehandle].preRollAvail)) {
			// Wait timed out after 0.5 secs.
			PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
			if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: No frame received after timed blocking wait of of 0.5 seconds.\n");

			// This is the end of a "up to 0.5 seconds blocking wait" style checkForImage of type 2.
			// Return "no new frame available yet". The calling code will retry the wait until its own
			// higher master timeout value is reached:
			return(FALSE);
		}
		
		// At this point we should have at least one frame available.
        // printf("PTB-DEBUG: After blocking fetch start %d\n", movieRecordBANK[moviehandle].frameAvail);
    }
	
    // We're here with at least one frame available and the mutex lock held.
	// Was this a pure "blocking check for new image"?
	if (checkForImage) {
		// Yes. Unlock mutex and signal success to caller - A new frame is ready.
		PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
		return(TRUE);
	}

	// If we reach this point, then at least 1 frame should be available and we are
	// asked to fetch it now and return it as a new OpenGL texture. The mutex is locked:

    // Preroll case is simple:
    movieRecordBANK[moviehandle].preRollAvail = 0;

    // Perform texture fetch & creation:
    if (oldstyle) {
		// Reset frame available flag:
		movieRecordBANK[moviehandle].frameAvail = 0;
		
		// This will retrieve an OpenGL compatible pointer to the pixel data and assign it to our texmemptr:
		out_texture->textureMemory = (GLuint*) movieRecordBANK[moviehandle].imageBuffer;
    } else {
	// Active playback mode?
	if (0 != rate) {
		// Active playback mode:
		if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: Pulling buffer from videosink, %d buffers decoded and queued.\n", movieRecordBANK[moviehandle].frameAvail);

		// Clamp frameAvail to maximum queue capacity, unless queue capacity is zero == "unlimited" capacity:
		if (((int) gst_app_sink_get_max_buffers(GST_APP_SINK(movieRecordBANK[moviehandle].videosink)) < movieRecordBANK[moviehandle].frameAvail) &&
			(gst_app_sink_get_max_buffers(GST_APP_SINK(movieRecordBANK[moviehandle].videosink)) > 0)) {
			movieRecordBANK[moviehandle].frameAvail = (int) gst_app_sink_get_max_buffers(GST_APP_SINK(movieRecordBANK[moviehandle].videosink));
		}

		// One less frame available after our fetch:
		movieRecordBANK[moviehandle].frameAvail--;

		// We can unlock early, thanks to videosink's internal buffering: XXX FIXME: Perfectly race-free to do this before the pull?
		PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);

		// This will pull the oldest video buffer from the videosink. It would block if none were available,
		// but that won't happen as we wouldn't reach this statement if none were available. It would return
		// NULL if the stream would be EOS or the pipeline off, but that shouldn't ever happen:
		videoBuffer = gst_app_sink_pull_buffer(GST_APP_SINK(movieRecordBANK[moviehandle].videosink));
	} else {
		// Passive fetch mode: Use prerolled buffers after seek:
		// These are available even after eos...

		// We can unlock early, thanks to videosink's internal buffering: XXX FIXME: Perfectly race-free to do this before the pull?
		PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);

		videoBuffer = gst_app_sink_pull_preroll(GST_APP_SINK(movieRecordBANK[moviehandle].videosink));
	}

	if (videoBuffer) {
		// Assign pointer to videoBuffer's data directly: Avoids one full data copy compared to oldstyle method.
		out_texture->textureMemory = (GLuint*) GST_BUFFER_DATA(videoBuffer);

		// Assign pts presentation timestamp in pipeline stream time and convert to seconds:
		movieRecordBANK[moviehandle].pts = (double) GST_BUFFER_TIMESTAMP(videoBuffer) / (double) 1e9;
		if (GST_CLOCK_TIME_IS_VALID(GST_BUFFER_DURATION(videoBuffer)))
			deltaT = (double) GST_BUFFER_DURATION(videoBuffer) / (double) 1e9;
		bufferIndex = GST_BUFFER_OFFSET(videoBuffer);
		if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: pts %f secs, dT %f secs, bufferId %i.\n", movieRecordBANK[moviehandle].pts, deltaT, (int) bufferIndex);
	} else {
		printf("PTB-ERROR: No new video frame received in gst_app_sink_pull_buffer! Something's wrong. Aborting fetch.\n");
		return(FALSE);
	}
	if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: ...done.\n");
    }

    PsychGetAdjustedPrecisionTimerSeconds(&tNow);
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Start of frame query to decode completion: %f msecs.\n", (tNow - tStart) * 1000.0);
    tStart = tNow;

    // Assign presentation_timestamp:
    if (presentation_timestamp) *presentation_timestamp = movieRecordBANK[moviehandle].pts;

    // Activate OpenGL context of target window:
    PsychSetGLContext(win);

    #if PSYCH_SYSTEM == PSYCH_OSX
    // Explicitely disable Apple's Client storage extensions. For now they are not really useful to us.
    glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, GL_FALSE);
    #endif

    // Build a standard PTB texture record:
    PsychMakeRect(out_texture->rect, 0, 0, movieRecordBANK[moviehandle].width, movieRecordBANK[moviehandle].height);    

    // Set NULL - special texture object as part of the PTB texture record:
    out_texture->targetSpecific.QuickTimeGLTexture = NULL;

    // Set texture orientation as if it were an inverted Offscreen window: Upside-down.
    out_texture->textureOrientation = 3;

    // We use zero client storage memory bytes:
    out_texture->textureMemorySizeBytes = 0;

    // Textures are aligned on 4 Byte boundaries because texels are RGBA8:
    out_texture->textureByteAligned = 4;

    // Assign texturehandle of our cached texture, if any, so it gets recycled now:
    out_texture->textureNumber = movieRecordBANK[moviehandle].cached_texture;

    if ((win->gfxcaps & kPsychGfxCapUYVYTexture) && useYUVDecode) {
	// GPU supports UYVY textures and we get data in that YCbCr format. Tell
	// texture creation routine to use this optimized format:
	if (!glewIsSupported("GL_APPLE_ycbcr_422")) {
	    // No support for more powerful Apple extension. Use Linux MESA extension:
	    out_texture->textureinternalformat = GL_YCBCR_MESA;
	    out_texture->textureexternalformat = GL_YCBCR_MESA;
	} else {
	    // Apple extension supported:
	    out_texture->textureinternalformat = GL_RGB;
	    out_texture->textureexternalformat = GL_YCBCR_422_APPLE;
	}
	// Same enumerant for Apple and Mesa:
	out_texture->textureexternaltype   = GL_UNSIGNED_SHORT_8_8_MESA;
    }

    // Let PsychCreateTexture() do the rest of the job of creating, setting up and
    // filling an OpenGL texture with content:
    PsychCreateTexture(out_texture);

    // After PsychCreateTexture() the cached texture object from our cache is used
    // and no longer available for recycling. We mark the cache as empty:
    // It will be filled with a new textureid for recycling if a texture gets
    // deleted in PsychMovieDeleteTexture()....
    movieRecordBANK[moviehandle].cached_texture = 0;

    PsychGetAdjustedPrecisionTimerSeconds(&tNow);
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Decode completion to texture created: %f msecs.\n", (tNow - tStart) * 1000.0);
    tStart = tNow;

    // Detection of dropped frames: This is a heuristic. We'll see how well it works out...
    // TODO: GstBuffer videoBuffer provides special flags that should allow to do a more
    // robust job, although nothing's wrong with the current approach per se...
    if (rate && presentation_timestamp) {
        // Try to check for dropped frames in playback mode:

        // Expected delta between successive presentation timestamps:
	// This is not dependent on playback rate, as it measures time in the
	// GStreamer movies timeline == Assuming 1x playback rate.
        targetdelta = 1.0f / movieRecordBANK[moviehandle].fps;

        // Compute real delta, given rate and playback direction:
        if (rate > 0) {
            realdelta = *presentation_timestamp - movieRecordBANK[moviehandle].last_pts;
            if (realdelta < 0) realdelta = 0;
        }
        else {
            realdelta = -1.0 * (*presentation_timestamp - movieRecordBANK[moviehandle].last_pts);
            if (realdelta < 0) realdelta = 0;
        }
        
        frames = realdelta / targetdelta;
        // Dropped frames?
        if (frames > 1 && movieRecordBANK[moviehandle].last_pts >= 0) {
            movieRecordBANK[moviehandle].nr_droppedframes += (int) (frames - 1 + 0.5);
        }

        movieRecordBANK[moviehandle].last_pts = *presentation_timestamp;
    }

    // Unlock.
    if (oldstyle) {
	PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
    } else {
	gst_buffer_unref(videoBuffer);
	videoBuffer = NULL;
    }
    
    // Manually advance movie time, if in fetch mode:
    if (0 == rate) {
	// We are in manual fetch mode: Need to manually advance movie to next
	// media sample:
	movieRecordBANK[moviehandle].endOfFetch = 0;
	preT = PsychGSGetMovieTimeIndex(moviehandle);
	event = gst_event_new_step(GST_FORMAT_BUFFERS, 1, 1.0, TRUE, FALSE);
	if (!gst_element_send_event(theMovie, event)) printf("PTB-DEBUG: In single-step seek I - Failed.\n");

	// Block until seek completed, failed, or timeout of 30 seconds reached:
	if (GST_STATE_CHANGE_FAILURE == gst_element_get_state(theMovie, NULL, NULL, (GstClockTime) (30 * 1e9))) printf("PTB-DEBUG: In single-step seek I - Failed.\n");
	postT = PsychGSGetMovieTimeIndex(moviehandle);

	if (PsychPrefStateGet_Verbosity() > 6) printf("PTB-DEBUG: Movie fetch advance: preT %f   postT %f  DELTA %lf %s\n", preT, postT, postT - preT, (postT - preT < 0.001) ? "SAME" : "DIFF");

	// Signal end-of-fetch if time no longer progresses signficiantly:
	if (postT - preT < 0.001) movieRecordBANK[moviehandle].endOfFetch = 1;
    }

    PsychGetAdjustedPrecisionTimerSeconds(&tNow);
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Texture created to fetch completion: %f msecs.\n", (tNow - tStart) * 1000.0);

    // Reset tStart for next fetch cycle:
    tStart = 0;

    return(TRUE);
}

/*
 *  PsychGSFreeMovieTexture() - Release texture memory for a Quicktime texture.
 *
 *  This routine is called by PsychDeleteTexture() in PsychTextureSupport.c
 *  It performs the special cleanup necessary for Quicktime created textures.
 *  As this ain't Quicktime but GStreamer there ain't nothing to do for us.
 */
void PsychGSFreeMovieTexture(PsychWindowRecordType *win)
{
	// Is this a GStreamer movietexture? If not, just skip this routine.
	if (win->windowType!=kPsychTexture || win->textureOrientation != 3 || win->texturecache_slot < 0) return;

	// Movie texture: Check if we can move it into our recycler cache
	// for later reuse...
	if (movieRecordBANK[win->texturecache_slot].cached_texture == 0) {
		// Cache free. Put this texture object into it for later reuse:
		movieRecordBANK[win->texturecache_slot].cached_texture = win->textureNumber;

		// 0-out the textureNumber so our standard cleanup routine (glDeleteTextures) gets
   	 	// skipped - if we wouldn't do this, our caching scheme would screw up.
		win->textureNumber = 0;
	}
	else {
		// Cache already occupied. We don't do anything but leave the cleanup work for
		// this texture to the standard PsychDeleteTexture() routine...
	}
	
    return;
}

/*
 *  PsychGSPlaybackRate() - Start- and stop movieplayback, set playback parameters.
 *
 *  moviehandle = Movie to start-/stop.
 *  playbackrate = zero == Stop playback, non-zero == Play movie with spec. rate,
 *                 e.g., 1 = forward, 2 = double speed forward, -1 = backward, ...
 *  loop = 0 = Play once. 1 = Loop, aka rewind at end of movie and restart.
 *  soundvolume = 0 == Mute sound playback, between 0.0 and 1.0 == Set volume to 0 - 100 %.
 *  Returns Number of dropped frames to keep playback in sync.
 */
int PsychGSPlaybackRate(int moviehandle, double playbackrate, int loop, double soundvolume)
{
    int			dropped = 0;
    GstElement		*theMovie = NULL;
    double timeindex;
    
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
        
    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;    
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }
    
    if (playbackrate != 0) {
        // Start playback of movie:

	// Set volume and mute state for audio:
	g_object_set(G_OBJECT(theMovie), "mute", (soundvolume <= 0) ? TRUE : FALSE, NULL);
	g_object_set(G_OBJECT(theMovie), "volume", soundvolume, NULL);

	// Set playback rate: An explicit seek to the position we are already (supposed to be)
	// is needed to avoid jumps in movies with bad encoding or keyframe placement:
	timeindex = PsychGSGetMovieTimeIndex(moviehandle);

	if (playbackrate > 0) {
		gst_element_seek(theMovie, playbackrate, GST_FORMAT_TIME, GST_SEEK_FLAG_FLUSH | GST_SEEK_FLAG_ACCURATE | ((loop & 0x1) ? GST_SEEK_FLAG_SEGMENT : 0), GST_SEEK_TYPE_SET,
				 (gint64) (timeindex * (double) 1e9), GST_SEEK_TYPE_NONE, GST_CLOCK_TIME_NONE);
	}
	else {
		gst_element_seek(theMovie, playbackrate, GST_FORMAT_TIME, GST_SEEK_FLAG_FLUSH | GST_SEEK_FLAG_ACCURATE | ((loop & 0x1) ? GST_SEEK_FLAG_SEGMENT : 0), GST_SEEK_TYPE_SET,
				 0, GST_SEEK_TYPE_SET, (gint64) (timeindex * (double) 1e9));
	}

        movieRecordBANK[moviehandle].loopflag = loop;
        movieRecordBANK[moviehandle].last_pts = -1.0;
        movieRecordBANK[moviehandle].nr_droppedframes = 0;
	movieRecordBANK[moviehandle].rate = playbackrate;
	movieRecordBANK[moviehandle].frameAvail = 0;
	movieRecordBANK[moviehandle].preRollAvail = 0;

	// Is this a movie with actual videotracks and frame-dropping on videosink full enabled?
	if ((movieRecordBANK[moviehandle].nrVideoTracks > 0) && gst_app_sink_get_drop(GST_APP_SINK(movieRecordBANK[moviehandle].videosink))) {
	    // Yes: We only schedule deferred start of playback at first Screen('GetMovieImage')
	    // frame fetch. This to avoid dropped frames due to random delays between
	    // call to Screen('PlayMovie') and Screen('GetMovieImage'):
	    movieRecordBANK[moviehandle].startPending = 1;
	}
	else {
	    // Only soundtrack or framedropping disabled with videotracks - Start it immediately:
	    movieRecordBANK[moviehandle].startPending = 0;
	    PsychMoviePipelineSetState(theMovie, GST_STATE_PLAYING, 10.0);
	    PsychGSProcessMovieContext(movieRecordBANK[moviehandle].MovieContext, FALSE);
	}
    }
    else {
	// Stop playback of movie:
	movieRecordBANK[moviehandle].rate = 0;
	movieRecordBANK[moviehandle].startPending = 0;
	movieRecordBANK[moviehandle].loopflag = 0;
	movieRecordBANK[moviehandle].endOfFetch = 0;

	PsychMoviePipelineSetState(theMovie, GST_STATE_PAUSED, 10.0);
	PsychGSProcessMovieContext(movieRecordBANK[moviehandle].MovieContext, FALSE);

        // Output count of dropped frames:
        if ((dropped=movieRecordBANK[moviehandle].nr_droppedframes) > 0) {
            if (PsychPrefStateGet_Verbosity()>2) {
		printf("PTB-INFO: Movie playback had to drop %i frames of movie %i to keep playback in sync.\n", movieRecordBANK[moviehandle].nr_droppedframes, moviehandle);
	    }
        }
    }

    return(dropped);
}

/*
 *  void PsychGSExitMovies() - Shutdown handler.
 *
 *  This routine is called by Screen('CloseAll') and on clear Screen time to
 *  do final cleanup. It deletes all textures and releases all movie objects.
 *
 */
void PsychGSExitMovies(void)
{
    PsychWindowRecordType	**windowRecordArray;
    int				i, numWindows; 
    
    // Release all Quicktime related OpenGL textures:
    PsychCreateVolatileWindowRecordPointerList(&numWindows, &windowRecordArray);
    for(i=0; i<numWindows; i++) {
        // Delete all Quicktime textures:
        if ((windowRecordArray[i]->windowType == kPsychTexture) && (windowRecordArray[i]->targetSpecific.QuickTimeGLTexture !=NULL)) { 
            PsychCloseWindow(windowRecordArray[i]);
        }
    }
    PsychDestroyVolatileWindowRecordPointerList(windowRecordArray);
    
    // Release all movies:
    PsychGSDeleteAllMovies();

    firsttime = TRUE;
    
    return;
}

/*
 *  PsychGSGetMovieTimeIndex()  -- Return current playback time of movie.
 */
double PsychGSGetMovieTimeIndex(int moviehandle)
{
    GstElement		*theMovie = NULL;
    GstFormat		fmt;
    gint64		pos_nsecs;

    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
    
    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;    
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }

    fmt = GST_FORMAT_TIME;
    if (!gst_element_query_position(theMovie, &fmt, &pos_nsecs)) {
	if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not query position in movie %i in seconds. Returning zero.\n", moviehandle);
	pos_nsecs = 0;
    }

    // Retrieve timeindex:
    return((double) pos_nsecs / (double) 1e9);
}

/*
 *  PsychGSSetMovieTimeIndex()  -- Set current playback time of movie, perform active seek if needed.
 */
double PsychGSSetMovieTimeIndex(int moviehandle, double timeindex, psych_bool indexIsFrames)
{
    GstElement		*theMovie;
    double		oldtime;
    gint64		targetIndex;
    GstEvent            *event;
    GstSeekFlags        flags;

    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
    
    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }
    
    // Retrieve current timeindex:
    oldtime = PsychGSGetMovieTimeIndex(moviehandle);

    // NOTE: We could use GST_SEEK_FLAG_SKIP to allow framedropping on fast forward/reverse playback...
    flags = GST_SEEK_FLAG_FLUSH | GST_SEEK_FLAG_ACCURATE;

    // Need segment seek flag for seek during active looped playback:
    flags|= ((movieRecordBANK[moviehandle].rate != 0) && (movieRecordBANK[moviehandle].loopflag & 0x1)) ? GST_SEEK_FLAG_SEGMENT : 0;

    // Index based or target time based seeking?
    if (indexIsFrames) {
	// Index based seeking:		
	targetIndex = (gint64) (timeindex + 0.5);

	// Simple seek, videobuffer (index) oriented, with pipeline flush and accurate seek,
	// i.e., not locked to keyframes, but frame-accurate:
	if (!gst_element_seek_simple(theMovie, GST_FORMAT_DEFAULT, flags, targetIndex)) {
		// Failed: This can happen on various movie formats as not all codecs and formats support frame-based seeks.
		// Fallback to time-based seek by faking a target time for given targetIndex:
		timeindex = (double) targetIndex / (double) movieRecordBANK[moviehandle].fps;

		if (PsychPrefStateGet_Verbosity() > 1) {
			printf("PTB-WARNING: Could not seek to frame index %i via frame-based seeking in movie %i.\n", (int) targetIndex, moviehandle);
			printf("PTB-WARNING: Will do a time-based seek to approximately equivalent time %f seconds instead.\n", timeindex);
			printf("PTB-WARNING: Not all movie formats support frame-based seeking. Please change your movie format for better precision.\n");
		}

		if (!gst_element_seek_simple(theMovie, GST_FORMAT_TIME, flags, (gint64) (timeindex * (double) 1e9)) &&
		    (PsychPrefStateGet_Verbosity() > 1)) {
			printf("PTB-WARNING: Time-based seek failed as well! Something is wrong with this movie!\n");
		}
	}
    }
    else {
	// Time based seeking:
	// Set new timeindex as time in seconds:

	// Simple seek, time-oriented, with pipeline flush and accurate seek,
	// i.e., not locked to keyframes, but frame-accurate:
	if (!gst_element_seek_simple(theMovie, GST_FORMAT_TIME, flags, (gint64) (timeindex * (double) 1e9)) &&
	    (PsychPrefStateGet_Verbosity() > 1)) {
		printf("PTB-WARNING: Time-based seek to %f seconds in movie %i failed. Something is wrong with this movie or the target time.\n", timeindex, moviehandle);
	}
    }

    // Block until seek completed, failed or timeout of 30 seconds reached:
    if (GST_STATE_CHANGE_FAILURE == gst_element_get_state(theMovie, NULL, NULL, (GstClockTime) (30 * 1e9)) &&
        (PsychPrefStateGet_Verbosity() > 1)) {
            printf("PTB-WARNING: SetTimeIndex on movie %i failed. Something is wrong with this movie or the target position. [Statechange-Failure in seek]\n", moviehandle);
            printf("PTB-WARNING: Requested target position was %f %s. This could happen if the movie is not efficiently seekable and a timeout was hit.\n",
                   timeindex, (indexIsFrames) ? "frames" : "seconds");
    }

    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Seeked to position %f secs in movie %i.\n", PsychGSGetMovieTimeIndex(moviehandle), moviehandle);

    // Reset fetch flag:
    movieRecordBANK[moviehandle].endOfFetch = 0;

    // Return old time value of previous position:
    return(oldtime);
}

// #ifdef PTB_USE_GSTREAMER
#endif
