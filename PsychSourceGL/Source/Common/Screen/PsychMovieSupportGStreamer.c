/*
    PsychSourceGL/Source/Common/Screen/PsychMovieSupportGStreamer.c

    PLATFORMS:  All

    AUTHORS:

    mario.kleiner.de@gmail.com      mk  Mario Kleiner

    HISTORY:

        28.11.2010    mk      Wrote it.
        20.08.2014    mk      Ported to GStreamer-1.4.x and later.

    DESCRIPTION:

    Psychtoolbox functions for dealing with movies. This is the operating system independent
    version which uses the GStreamer media framework, version 1.4 or later.

    These PsychGSxxx functions are called from the dispatcher in
    Common/Screen/PsychMovieSupport.[hc].

*/

#ifdef PTB_USE_GSTREAMER

#include "Screen.h"
#include <glib.h>
#include "PsychMovieSupportGStreamer.h"
#include <gst/gst.h>

#if GST_CHECK_VERSION(1,0,0)
#include <gst/app/gstappsink.h>
#include <gst/video/video.h>

// Need to define this for playbin as it is not defined
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

#define PSYCH_MAX_MOVIES 100

typedef struct {
    psych_mutex         mutex;
    psych_condition     condition;
    double              pts;
    GstElement          *theMovie;
    GMainLoop           *MovieContext;
    GstElement          *videosink;
    PsychWindowRecordType* parentRecord;
    unsigned char       *imageBuffer;
    int                 frameAvail;
    int                 preRollAvail;
    double              rate;
    int                 startPending;
    int                 endOfFetch;
    int                 specialFlags1;
    int                 pixelFormat;
    int                 loopflag;
    double              movieduration;
    int                 nrframes;
    double              fps;
    int                 width;
    int                 height;
    double              aspectRatio;
    int                 bitdepth;
    double              last_pts;
    int                 nr_droppedframes;
    int                 nrAudioTracks;
    int                 nrVideoTracks;
    char                movieLocation[FILENAME_MAX];
    char                movieName[FILENAME_MAX];
    GLuint              cached_texture;
} PsychMovieRecordType;

static PsychMovieRecordType movieRecordBANK[PSYCH_MAX_MOVIES];
static int numMovieRecords = 0;
static psych_bool firsttime = TRUE;

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

    // Note: This is deprecated and not needed anymore on GLib 2.31.0 and later, as
    // GLib's threading system auto-initializes on first use since that version. We
    // keep it for now to stay compatible to older systems, e.g., Ubuntu 10.04 LTS,
    // conditionally on the GLib version we build against:
    #if !GLIB_CHECK_VERSION (2, 31, 0)
        // Initialize GLib's threading system early:
        if (!g_thread_supported()) g_thread_init(NULL);
    #endif

    return;
}

int PsychGSGetMovieCount(void) {
    return(numMovieRecords);
}

// Forward declaration:
static gboolean PsychMovieBusCallback(GstBus *bus, GstMessage *msg, gpointer dataptr);

/* Perform context loop iterations (for bus message handling) if doWait == false,
 * as long as there is work to do, or at least two seconds worth of iterations
 * if doWait == true. This drives the message-bus callback, so needs to be
 * performed to get any error reporting etc.
 */
int PsychGSProcessMovieContext(PsychMovieRecordType* movie, psych_bool doWait)
{
    GstBus* bus;
    GstMessage *msg;
    psych_bool workdone = FALSE;
    double tdeadline, tnow;
    PsychGetAdjustedPrecisionTimerSeconds(&tdeadline);
    tnow = tdeadline;
    tdeadline+=2.0;

    // New style:
    bus = gst_pipeline_get_bus(GST_PIPELINE(movie->theMovie));
    msg = NULL;

    // If doWait, try to perform iterations until 2 seconds elapsed or at least one event handled:
    while (doWait && (tnow < tdeadline) && !gst_bus_have_pending(bus)) {
        // Update time:
        PsychYieldIntervalSeconds(0.010);
        PsychGetAdjustedPrecisionTimerSeconds(&tnow);
    }

    msg = gst_bus_pop(bus);
    while (msg) {
        workdone = TRUE;
        PsychMovieBusCallback(bus, msg, movie);
        gst_message_unref(msg);
        msg = gst_bus_pop(bus);
    }

    gst_object_unref(bus);

    return(workdone);
}

/* Initiate pipeline state changes: Startup, Preroll, Playback, Pause, Standby, Shutdown. */
static psych_bool PsychMoviePipelineSetState(GstElement* theMovie, GstState state, double timeoutSecs)
{
    GstState state_pending;
    GstStateChangeReturn rcstate;

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

psych_bool PsychIsMovieSeekable(PsychMovieRecordType* movie)
{
    GstQuery *query;
    gint64 start, end;
    gboolean seekable = FALSE;

    query = gst_query_new_seeking(GST_FORMAT_TIME);
    if (gst_element_query(movie->theMovie, query)) {
        gst_query_parse_seeking(query, NULL, &seekable, &start, &end);
        if (seekable) {
            if (PsychPrefStateGet_Verbosity() > 4) {
                printf("PTB-DEBUG: Seeking is enabled from %" GST_TIME_FORMAT " to %" GST_TIME_FORMAT "\n",
                        GST_TIME_ARGS (start), GST_TIME_ARGS (end));
            }
        }
        else {
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Seeking is disabled for this movie stream.\n");
        }
    }
    else {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Seeking query failed!\n");
    }
    gst_query_unref(query);

    return((psych_bool) seekable);
}

/* Receive messages from the playback pipeline message bus and handle them: */
static gboolean PsychMovieBusCallback(GstBus *bus, GstMessage *msg, gpointer dataptr)
{
    GstSeekFlags rewindFlags = 0;
    PsychMovieRecordType* movie = (PsychMovieRecordType*) dataptr;

    (void) bus;

    if (PsychPrefStateGet_Verbosity() > 11) printf("PTB-DEBUG: PsychMovieBusCallback: Msg source name and type: %s : %s\n", GST_MESSAGE_SRC_NAME(msg), GST_MESSAGE_TYPE_NAME(msg));

    switch (GST_MESSAGE_TYPE (msg)) {
        case GST_MESSAGE_SEGMENT_DONE:
            // We usually receive segment done message instead of eos if looped playback is active and
            // the end of the stream is approaching, so we fallthrough to message eos for rewinding...
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: PsychMovieBusCallback: Message SEGMENT_DONE received.\n");
            // Fall through.
        case GST_MESSAGE_EOS: {
            // Rewind at end of movie if looped playback enabled:
            if ((GST_MESSAGE_TYPE (msg) == GST_MESSAGE_EOS) && (PsychPrefStateGet_Verbosity() > 4)) printf("PTB-DEBUG: PsychMovieBusCallback: Message EOS received.\n");

            // Looping via seek requested (method 0x1) and playback active?
            if ((movie->loopflag & 0x1) && (movie->rate != 0)) {
                // Perform loop via rewind via seek:
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: PsychMovieBusCallback: End of iteration in active looped playback reached: Rewinding...\n");

                // Seek: We normally don't GST_SEEK_FLAG_FLUSH here, so the rewinding is smooth because we don't throw away buffers queued in the pipeline,
                // unless we are at the end of the stream (EOS), so there ain't anything queued in the pipeline, or code requests an explicit pipeline flush via flag 0x8.
                // This seems to make no sense (why flush an already EOS - empty pipeline?) but is neccessary for some movies with sound on some systems:
                if ((GST_MESSAGE_TYPE(msg) == GST_MESSAGE_EOS) || (movie->loopflag & 0x8)) rewindFlags |= GST_SEEK_FLAG_FLUSH;

                // On some movies and configurations, we need a segment seek as indicated by flag 0x4:
                if (movie->loopflag & 0x4) rewindFlags |= GST_SEEK_FLAG_SEGMENT;

                // Seek method depends on playback direction:
                if (movie->rate > 0) {
                    if (!gst_element_seek(movie->theMovie, movie->rate, GST_FORMAT_TIME, GST_SEEK_FLAG_ACCURATE | rewindFlags, GST_SEEK_TYPE_SET, 0, GST_SEEK_TYPE_NONE, GST_CLOCK_TIME_NONE)) {
                        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Rewinding video in forward playback failed!\n");
                    }
                }
                else {
                    if (!gst_element_seek(movie->theMovie, movie->rate, GST_FORMAT_TIME, GST_SEEK_FLAG_ACCURATE | rewindFlags, GST_SEEK_TYPE_NONE, GST_CLOCK_TIME_NONE, GST_SEEK_TYPE_END, 0)) {
                        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Rewinding video in reverse playback failed!\n");
                    }
                }
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
                    printf("PTB-ERROR: Could not open movie file [%s] for playback! No such moviefile with the given absolute path and filename.\n",
                           movie->movieName);
                    printf("PTB-ERROR: Please note that you *must* provide an absolute path and filename for your movie file, filename alone won't work.\n");
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

/* Called at each end-of-stream event at end of playback: */
static void PsychEOSCallback(GstAppSink *sink, gpointer user_data)
{
    (void) sink, (void) user_data;

    //PsychMovieRecordType* movie = (PsychMovieRecordType*) user_data;
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
    (void) sink;

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
    (void) sink;

    PsychLockMutex(&movie->mutex);
    //printf("PTB-DEBUG: New Buffer received.\n");
    movie->frameAvail++;
    PsychSignalCondition(&movie->condition);
    PsychUnlockMutex(&movie->mutex);

    return(GST_FLOW_OK);
}

/* Not used by us, but needs to be defined as no-op anyway: */
/* // There are only 3 function pointers in GstAppSinkCallbacks now
static GstFlowReturn PsychNewBufferListCallback(GstAppSink *sink, gpointer user_data)
{
    (void) sink, (void) user_data;

    //PsychMovieRecordType* movie = (PsychMovieRecordType*) user_data;
    //PsychLockMutex(&movie->mutex);
    //printf("PTB-DEBUG: New Bufferlist received.\n");
    //PsychSignalCondition(&movie->condition);
    //PsychUnlockMutex(&movie->mutex);

    return(GST_FLOW_OK);
}
*/

/* Not used by us, but needs to be defined as no-op anyway: */
static void PsychDestroyNotifyCallback(gpointer user_data)
{
    (void) user_data;

    return;
}

/* This callback is called when the pipeline is about to finish playback
 * of the current movie stream. If looped playback via method 0x2 is enabled,
 * this needs to trigger a repetition by rescheduling the movie URI for playback.
 *
 * Allows gapless playback, but doesn't work reliable on all media types.
 *
 */
static void PsychMovieAboutToFinishCB(GstElement *theMovie, gpointer user_data)
{
    PsychMovieRecordType* movie = (PsychMovieRecordType*) user_data;

    // Loop method 0x2 active?
    if ((movie->loopflag & 0x2) && (movie->rate != 0)) {
        g_object_set(G_OBJECT(theMovie), "uri", movie->movieLocation, NULL);
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: About-to-finish received: Rewinding via uri method.\n");
    }

    return;
}

/* Not used, didn't work, but left here in case we find a use for it in the future. */
/*
static void PsychMessageErrorCB(GstBus *bus, GstMessage *msg)
{
    gchar  *debug;
    GError *error;

    (void) bus;

    gst_message_parse_error (msg, &error, &debug);
    g_free (debug);

    printf("PTB-BUSERROR: %s\n", error->message);
    g_error_free (error);
    return;
}
*/

static GstAppSinkCallbacks videosinkCallbacks = {
    PsychEOSCallback,
    PsychNewPrerollCallback,
    PsychNewBufferCallback,
    {0}
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
 *      pixelFormat = As passed to 'OpenMovie'
 *      maxNumberThreads = Maximum number of decode threads to use (0 = auto, 1 = One, ...);
 *      movieOptions = Options string with additional options for movie playback.
 */
void PsychGSCreateMovie(PsychWindowRecordType *win, const char* moviename, double preloadSecs, int* moviehandle, int asyncFlag, int specialFlags1, int pixelFormat, int maxNumberThreads, char* movieOptions)
{
    GValue          item = G_VALUE_INIT;
    GstCaps         *colorcaps = NULL;
    GstElement      *theMovie = NULL;
    GstElement      *videocodec = NULL;
    GstElement      *videosink = NULL;
    GstElement      *audiosink;
    gchar*          pstring;
    gint64          length_format;
    GstPad          *pad, *peerpad;
    const GstCaps   *caps;
    GstStructure    *str;
    gint            width,height;
    gint            rate1, rate2;
    int             i, slotid;
    int             max_video_threads;
    char            movieLocation[FILENAME_MAX];
    psych_bool      printErrors;
    GstIterator     *it;
    psych_bool      done;
    GstPlayFlags    playflags = 0;
    psych_bool      needCodecSetup = FALSE;

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

        firsttime = FALSE;
    }

    if (win && !PsychIsOnscreenWindow(win)) {
        if (printErrors)
            PsychErrorExitMsg(PsychError_user, "Provided windowPtr is not an onscreen window.");
        else
            return;
    }

    // As a side effect of some PsychGSCheckInit() some broken GStreamer runtimes can change
    // the OpenGL context binding behind our back to some GStreamer internal context.
    // Make sure our own context is bound after return from PsychGSCheckInit() to protect
    // against the state bleeding this would cause:
    if (win) PsychSetGLContext(win);

    if (NULL == moviename) {
        if (printErrors)
            PsychErrorExitMsg(PsychError_internal, "NULL-Ptr instead of moviename passed!");
        else
            return;
    }

    if (numMovieRecords >= PSYCH_MAX_MOVIES) {
        *moviehandle = -2;
        if (printErrors)
            PsychErrorExitMsg(PsychError_user, "Allowed maximum number of simultaneously open movies exceeded!");
        else
            return;
    }

    // Search first free slot in movieRecordBANK:
    for (i = 0; (i < PSYCH_MAX_MOVIES) && (movieRecordBANK[i].theMovie); i++);
    if (i >= PSYCH_MAX_MOVIES) {
        *moviehandle = -2;
        if (printErrors)
            PsychErrorExitMsg(PsychError_user, "Allowed maximum number of simultaneously open movies exceeded!");
        else
            return;
    }

    // Slot slotid will contain the movie record for our new movie object:
    slotid=i;

    // Zero-out new record in moviebank:
    memset(&movieRecordBANK[slotid], 0, sizeof(PsychMovieRecordType));

    // Store specialFlags1 from open call:
    movieRecordBANK[slotid].specialFlags1 = specialFlags1;

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
        // Use playbin:
        theMovie = gst_element_factory_make("playbin", "ptbmovieplaybackpipeline");
        movieRecordBANK[slotid].theMovie = theMovie;
        if (theMovie == NULL) {
            printf("PTB-ERROR: Failed to create GStreamer playbin element! Your GStreamer installation is\n");
            printf("PTB-ERROR: incomplete or damaged and misses at least the gst-plugins-base set of plugins!\n");
            if (printErrors)
                PsychErrorExitMsg(PsychError_system, "Opening the movie failed. GStreamer configuration problem.");
            else
                return;
        }

        // Assign name of movie to play:
        g_object_set(G_OBJECT(theMovie), "uri", movieLocation, NULL);

        // Default flags for playbin: Decode video ...
        playflags = GST_PLAY_FLAG_VIDEO;

        // ... and deinterlace it if needed, unless prevented by specialFlags setting 256:
        if (!(specialFlags1 & 256)) playflags|= GST_PLAY_FLAG_DEINTERLACE;

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
        movieRecordBANK[slotid].theMovie = theMovie;
        videosink = gst_bin_get_by_name(GST_BIN(theMovie), "ptbsink0");
        printf("LAUNCHLINE[%p]: %s\n", videosink, movieLocation);
    }

    // Assign a fakesink named "ptbsink0" as destination video-sink for
    // all video content. This allows us to get hold of the video frame buffers for
    // converting them into PTB OpenGL textures:
    if (!videosink)
        videosink = gst_element_factory_make ("appsink", "ptbsink0");

    if (!videosink) {
        printf("PTB-ERROR: Failed to create video-sink appsink ptbsink! Your GStreamer installation is\n");
        printf("PTB-ERROR: incomplete or damaged and misses at least the gst-plugins-base set of plugins!\n");
        PsychGSProcessMovieContext(&(movieRecordBANK[slotid]), TRUE);
        if (printErrors)
            PsychErrorExitMsg(PsychError_system, "Opening the movie failed. Reason hopefully given above.");
        else
            return;
    }

    movieRecordBANK[slotid].videosink = videosink;

    // Setting flag 1 in specialFlags1 is equivalent to setting pixelFormat == 5, to retain
    // backwards compatibility to previous ptb releases:
    if (specialFlags1 & 0x1) pixelFormat = 5;

    // Assign initial pixelFormat to use, as requested by usercode:
    movieRecordBANK[slotid].pixelFormat = pixelFormat;

    // Default to 8 bpc bitdepth as a starter:
    movieRecordBANK[slotid].bitdepth = 8;

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
    // a RGBA8 texture. Some GPU + driver combos do support direct handling of UYVU YCrCb 4:2:2
    // packed pixel data as textures. If we are on such a GPU we request UYVY data and upload it
    // directly in this format to the GPU. This more efficient both for the GStreamers decode
    // pipeline, and the later Videobuffer -> OpenGL texture conversion:
    if (win && (win->gfxcaps & kPsychGfxCapUYVYTexture) && (movieRecordBANK[slotid].pixelFormat == 5)) {
        // GPU supports handling and decoding of UYVY type yuv textures: We use these,
        // as they are more efficient to decode and handle by typical video codecs:
        colorcaps = gst_caps_new_simple ("video/x-raw",
                                         "format", G_TYPE_STRING, "UYVY",
                                         NULL);

        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Movie playback for movie %i will use UYVY YCrCb 4:2:2 textures for optimized decode and rendering.\n", slotid);
    }
    else if ((movieRecordBANK[slotid].pixelFormat == 6) && win && (win->gfxcaps & kPsychGfxCapFBO) && PsychAssignPlanarI420TextureShader(NULL, win)) {
        // Usercode wants YUV I420 planar encoded format and GPU suppports needed fragment shaders and FBO's.
        // Ask for I420 decoded video. This is the native output format of HuffYUV and H264 codecs, so using it
        // allows to skip colorspace conversion in GStreamer. The format is also highly efficient for texture
        // creation and upload to the GPU, but requires a fragment shader for colorspace conversion during drawing:
        colorcaps = gst_caps_new_simple ("video/x-raw",
                                         "format", G_TYPE_STRING, "I420",
                                         NULL);
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Movie playback for movie %i will use YUV-I420 planar textures for optimized decode and rendering.\n", slotid);
    }
    else if (((movieRecordBANK[slotid].pixelFormat == 7) || (movieRecordBANK[slotid].pixelFormat == 8)) && win && (win->gfxcaps & kPsychGfxCapFBO) && PsychAssignPlanarI800TextureShader(NULL, win)) {
        // Usercode wants Y8/Y800 planar encoded format and GPU suppports needed fragment shaders and FBO's.
        // Ask for I420 or Y800 decoded video. I420 is the native output format of HuffYUV and H264 codecs, so using it
        // allows to skip colorspace conversion in GStreamer. The format is also highly efficient for texture
        // creation and upload to the GPU, but requires a fragment shader for colorspace conversion during drawing:
        // Note: The FOURCC 'Y800' is equivalent to 'Y8  ' and 'GREY' as far as we know. As of June 2012, using the
        // Y800 decoding doesn't have any performance benefits compared to I420 decoding, actually performance is a
        // little bit lower. Apparently the video codecs don't take the Y800 format into account, ie., they don't skip
        // decoding of the chroma components, instead they probably decode as usual and the colorspace conversion then
        // throws away the unwanted chroma planes, causing possible extra overhead for discarding this data. We leave
        // the option in anyway, because there may be codecs (possibly in future GStreamer versions) that can take
        // advantage of Y800 format for higher performance.
        colorcaps = gst_caps_new_simple ("video/x-raw",
                                         "format", G_TYPE_STRING,
                                         (movieRecordBANK[slotid].pixelFormat == 8) ? "Y800" : "I420",
                                         NULL);
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Movie playback for movie %i will use Y8-I800 planar textures for optimized decode and rendering.\n", slotid);
    }
    else {
        // GPU does not support yuv textures or shader based decoding. Need to go brute-force and convert
        // video into RGBA8 format:

        // Force unsupportable formats to RGBA8 aka format 4, except for formats 7/8, which map
        // to 1 == L8, and map 2 == LA8 to 1 == L8:
        switch (movieRecordBANK[slotid].pixelFormat) {
            case 1: // 1 is fine.
                break;
            case 3: // 3 is handled after switch-case.
                break;
            case 4: // 4 is fine.
                break;
            case 2:
            case 7:
            case 8:
                movieRecordBANK[slotid].pixelFormat = 1;
                break;
            case 5:
            case 6:
                movieRecordBANK[slotid].pixelFormat = 4;
                break;
            case 9:     // 9 and 10 are fine:
            case 10:
                break;

            default:
                if (printErrors)
                    PsychErrorExitMsg(PsychError_user, "Invalid 'pixelFormat' parameter specified!");
                else
                    // Revert to something safe, RGBA8, as we can not error abort here:
                    movieRecordBANK[slotid].pixelFormat = 4;
                break;
        }

        // Map 3 == RGB8 to 4 == RGBA8, unless it is our special proprietary 16 bpc encoding:
        if ((movieRecordBANK[slotid].pixelFormat == 3) && !(specialFlags1 & 512)) movieRecordBANK[slotid].pixelFormat = 4;

        // At this point we can have any of these pixelFormats: 1, 3, 4, 9, 10. Handle them:

        if (movieRecordBANK[slotid].pixelFormat == 4) {
            // Use RGBA8 format:
            colorcaps = gst_caps_new_simple("video/x-raw",
                                            "format", G_TYPE_STRING, "BGRA",
                                            NULL);
            if ((PsychPrefStateGet_Verbosity() > 3) && (pixelFormat == 5)) printf("PTB-INFO: Movie playback for movie %i will use RGBA8 textures due to lack of YUV-422 texture support on GPU.\n", slotid);
            if ((PsychPrefStateGet_Verbosity() > 3) && (pixelFormat == 6)) printf("PTB-INFO: Movie playback for movie %i will use RGBA8 textures due to lack of YUV-I420 support on GPU.\n", slotid);
            if ((PsychPrefStateGet_Verbosity() > 3) && ((pixelFormat == 7) || (pixelFormat == 8))) printf("PTB-INFO: Movie playback for movie %i will use L8 textures due to lack of Y8-I800 support on GPU.\n", slotid);

            if ((PsychPrefStateGet_Verbosity() > 3) && !(pixelFormat < 5)) printf("PTB-INFO: Movie playback for movie %i will use RGBA8 textures.\n", slotid);
        }

        if ((movieRecordBANK[slotid].pixelFormat == 1) && !(specialFlags1 & 512)) {
            // Use LUMINANCE8 format:
            colorcaps = gst_caps_new_simple("video/x-raw",
                                            "format", G_TYPE_STRING, "GRAY8",
                                            NULL);
            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Movie playback for movie %i will use L8 luminance textures.\n", slotid);
        }

        // Psychtoolbox proprietary 16 bpc pixelformat for 1 or 3 channel data?
        if ((pixelFormat == 1 || pixelFormat == 3) && (specialFlags1 & 512)) {
            // Yes. Need to always decode as RGB8 24 bpp: Texture creation will then handle this further.
            colorcaps = gst_caps_new_simple("video/x-raw",
                                            "format", G_TYPE_STRING, "RGB",
                                            NULL);
            movieRecordBANK[slotid].pixelFormat = pixelFormat;
        }

        if (movieRecordBANK[slotid].pixelFormat == 9) {
            // Use GRAY 16 bpc format for 16 bpc decoding:
            colorcaps = gst_caps_new_simple("video/x-raw",
                                            "format", G_TYPE_STRING, "GRAY16_LE",
                                            NULL);

            // Switch to 16 bpc bitdepth and single channel pixelFormat:
            movieRecordBANK[slotid].bitdepth = 16;
            movieRecordBANK[slotid].pixelFormat = 1;

            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Movie playback for movie %i will use 16 bpc LUMINANCE 16 float textures.\n", slotid);
        }

        if (movieRecordBANK[slotid].pixelFormat == 10) {
            // Use ARGB64 16 bpc per color channel format for 16 bpc decoding:
            colorcaps = gst_caps_new_simple("video/x-raw",
                                            "format", G_TYPE_STRING, "ARGB64",
                                            NULL);

            // Switch to 16 bpc bitdepth and 4-channel RGBA pixelFormat:
            movieRecordBANK[slotid].bitdepth = 16;
            movieRecordBANK[slotid].pixelFormat = 4;

            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Movie playback for movie %i will use 16 bpc RGBA 16 bpc float textures.\n", slotid);
        }
    }

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
    pad = gst_element_get_static_pad(videosink, "sink");

    PsychGSProcessMovieContext(&(movieRecordBANK[slotid]), FALSE);

    // Attach custom made audio sink?
    if ((pstring = strstr(movieOptions, "AudioSink="))) {
        pstring += strlen("AudioSink=");
        pstring = strdup(pstring);
        if (strstr(pstring, ":::") != NULL) *(strstr(pstring, ":::")) = 0;
        audiosink = gst_parse_bin_from_description((const gchar *) pstring, TRUE, NULL);
        if (audiosink) {
            g_object_set(G_OBJECT(theMovie), "audio-sink", audiosink, NULL);
            audiosink = NULL;
            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Attached custom audio sink for playback of movie. Spec: '%s'\n", pstring);
            free(pstring); pstring = NULL;
        }
        else {
            printf("PTB-ERROR: Could not create requested audio sink for playback of movie! Failing sink spec was: '%s'\n", pstring);
            free(pstring); pstring = NULL;
            if (printErrors)
                PsychErrorExitMsg(PsychError_user, "Failed to create custom audio sink for movie playback.");
            else
                return;
        }
    }

    // Preload / Preroll the pipeline:
    if (!PsychMoviePipelineSetState(theMovie, GST_STATE_PAUSED, 30.0)) {
        PsychGSProcessMovieContext(&(movieRecordBANK[slotid]), TRUE);
        if (printErrors)
            PsychErrorExitMsg(PsychError_user, "In OpenMovie: Opening the movie failed I. Reason given above.");
        else
            return;
    }

    // Set video decoder parameters:

    // First we try to find the video decoder plugin in the media graph. We iterate over all
    // elements and try to find one whose properties match known properties of video codecs:
    // Note: The typeid/name or classid/name proved useless for finding the video codec, as
    // the video codecs are not derived from some "video codec class", but each codec defines
    // its own type and class.
    it = gst_bin_iterate_recurse(GST_BIN(theMovie));
    done = FALSE;
    videocodec = NULL;

    while (!done) {
        switch (gst_iterator_next(it, &item)) {
            case GST_ITERATOR_OK:
                videocodec = g_value_peek_pointer(&item);
                if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In pipeline: Child element name: %s\n", (const char*) gst_object_get_name(GST_OBJECT(videocodec)));
                // Our current match critera for video codecs: Having at least one of the properties "max-threads" or "lowres":
                if ((g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "max-threads")) ||
                    (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "lowres")) ||
                    (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "debug-mv")) ||
                    (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "skip-frame"))) {
                    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Found video decoder element %s.\n", (const char*) gst_object_get_name(GST_OBJECT(videocodec)));
                    done = TRUE;
                } else {
                    videocodec = NULL;
                    g_value_reset(&item);
                }
            break;

            case GST_ITERATOR_RESYNC:
                gst_iterator_resync(it);
            break;

            case GST_ITERATOR_DONE:
                done = TRUE;
            break;

            default:
                videocodec = NULL;
        }
    }

    g_value_unset(&item);
    gst_iterator_free(it);
    it = NULL;

    // Check if some codec properties need to be changed.
    // This happens if usercode provides some supported non-default override parameter for the codec,
    // or if the codec is multi-threaded and usercode wants us to configure its multi-threading behaviour:
    needCodecSetup = FALSE;
    if (videocodec &&  ((g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "max-threads") && (maxNumberThreads > -1)) ||
                       (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "lowres") && (specialFlags1 & (0))) || /* MK: 'lowres' disabled for now. */
                       (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "debug-mv") && (specialFlags1 & 4)) ||
                       (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "skip-frame") && (specialFlags1 & 8))
                       )) {
        needCodecSetup = TRUE;
    }

    // Set videocodec state to "ready" if parameter change is needed, as the codec only
    // accepts the new settings in that state:
    if (needCodecSetup) {
        // Ready the video codec, so a new max thread count or other parameters can be set:
        if (!PsychMoviePipelineSetState(videocodec, GST_STATE_READY, 30.0)) {
            PsychGSProcessMovieContext(&(movieRecordBANK[slotid]), TRUE);
            if (printErrors)
                PsychErrorExitMsg(PsychError_user, "In OpenMovie: Opening the movie failed III. Reason given above.");
            else
                return;
        }
    }

    // Drawing of motion vectors requested by usercode specialflags1 flag 4? If so, enable it:
    if (needCodecSetup && (specialFlags1 & 4) && (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "debug-mv"))) {
        g_object_set(G_OBJECT(videocodec), "debug-mv", 1, NULL);
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Playback for movie %i will display motion vectors for visualizaton of optic flow.\n", slotid);
    }

    // Skipping of B-Frames video decoding requested by usercode specialflags1 flag 8? If so, enable skipping:
    if (needCodecSetup && (specialFlags1 & 8) && (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "skip-frame"))) {
        g_object_set(G_OBJECT(videocodec), "skip-frame", 1, NULL);
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Playback for movie %i will skip B-Frames during video decoding for higher speed.\n", slotid);
    }

    // Multi-threaded codec and usercode requests setup? If so, set its multi-threading behaviour:
    // By default many codecs would only use one single thread on any system, even if they are multi-threading capable.
    if (needCodecSetup && (g_object_class_find_property(G_OBJECT_GET_CLASS(videocodec), "max-threads")) && (maxNumberThreads > -1)) {
        max_video_threads = 1;
        g_object_get(G_OBJECT(videocodec), "max-threads", &max_video_threads, NULL);
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Movie playback for movie %i uses video decoder with a default maximum number of %i processing threads.\n", slotid, max_video_threads);

        // Specific number of threads requested, or zero for auto-select?
        if (maxNumberThreads > 0) {
            // Specific value provided: Use it.
            max_video_threads = maxNumberThreads;
            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Setting video decoder to use a maximum of %i processing threads.\n", max_video_threads);
        } else {
            // Default behaviour: A settig of zero asks GStreamer to auto-detect the optimal
            // setting for the given host computer -- typically number of threads == number of processor cores:
            max_video_threads = 0;
            if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Setting video decoder to use auto-selected optimal number of processing threads.\n");
        }

        // Assign new setting:
        g_object_set(G_OBJECT(videocodec), "max-threads", max_video_threads, NULL);

        // Requery:
        g_object_get(G_OBJECT(videocodec), "max-threads", &max_video_threads, NULL);
        if (PsychPrefStateGet_Verbosity() > 4) {
            if (max_video_threads != 0) {
                printf("PTB-INFO: Movie playback for movie %i uses video decoder with a current maximum number of %i processing threads.\n", slotid, max_video_threads);
            } else {
                printf("PTB-INFO: Movie playback for movie %i uses video decoder with a current auto-detected optimal number of processing threads.\n", slotid);
            }
        }
    }

    // Bring codec back to paused state, so it is ready to do its job with the
    // new codec parameters set:
    if (needCodecSetup) {
        // Pause the video codec, so the new max thread count is accepted:
        if (!PsychMoviePipelineSetState(videocodec, GST_STATE_PAUSED, 30.0)) {
            PsychGSProcessMovieContext(&(movieRecordBANK[slotid]), TRUE);
            if (printErrors)
                PsychErrorExitMsg(PsychError_user, "In OpenMovie: Opening the movie failed IV. Reason given above.");
            else
                return;
        }
    }

    // NULL out videocodec, we must not unref it, as we didn't aquire or own private ref:
    videocodec = NULL;

    // Query number of available video and audio tracks in movie:
    g_object_get (G_OBJECT(theMovie),
               "n-video", &movieRecordBANK[slotid].nrVideoTracks,
               "n-audio", &movieRecordBANK[slotid].nrAudioTracks,
                NULL);

    // We need a valid onscreen window handle for real video playback:
    if ((NULL == win) && (movieRecordBANK[slotid].nrVideoTracks > 0)) {
        if (printErrors)
            PsychErrorExitMsg(PsychError_user, "No windowPtr to an onscreen window provided. Must do so for movies with videotrack!");
        else
            return;
    }

    PsychGSProcessMovieContext(&(movieRecordBANK[slotid]), FALSE);

    PsychInitMutex(&movieRecordBANK[slotid].mutex);
    PsychInitCondition(&movieRecordBANK[slotid].condition, NULL);

    // Install callbacks used by the videosink (appsink) to announce various events:
    gst_app_sink_set_callbacks(GST_APP_SINK(videosink), &videosinkCallbacks, &(movieRecordBANK[slotid]), PsychDestroyNotifyCallback);

    // Assign harmless initial settings for fps and frame size:
    rate1 = 0;
    rate2 = 1;
    width = height = 0;

    // Videotrack available?
    if (movieRecordBANK[slotid].nrVideoTracks > 0) {
        // Yes: Query size and framerate of movie:
        peerpad = gst_pad_get_peer(pad);
        caps=gst_pad_get_current_caps(peerpad);
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
    if (gst_element_query_duration(theMovie, GST_FORMAT_TIME, &length_format)) {
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
        gst_app_sink_set_max_buffers(GST_APP_SINK(videosink),
                                     ((movieRecordBANK[slotid].fps > 0) && (preloadSecs >= 0)) ? ((int) (movieRecordBANK[slotid].fps * preloadSecs) + 1) : 0);
    }
    else {
        // No: Only allow one queued buffer before dropping, to avoid optimal audio-video sync:
        gst_app_sink_set_max_buffers(GST_APP_SINK(videosink), 1);
    }

    // Compute framecount from fps and duration:
    movieRecordBANK[slotid].nrframes = (int)(movieRecordBANK[slotid].fps * movieRecordBANK[slotid].movieduration + 0.5);
    //printf("PTB-DEBUG: Number of frames in movie %i [%s] is %i.\n", slotid, moviename, movieRecordBANK[slotid].nrframes);

    // Is this movie supposed to be encoded in Psychtoolbox special proprietary "16 bpc stuffed into 8 bpc" format?
    if (specialFlags1 & 512) {
        // Yes. Invert the hacks applied during encoding/writing of movie:

        // Only 1 layer and 3 layer are supported:
        if (pixelFormat != 1 && pixelFormat !=3) {
            if (printErrors)
                PsychErrorExitMsg(PsychError_user, "You specified 'specialFlags1' setting 512 for Psychtoolbox proprietary 16 bpc decoding, but pixelFormat is not 1 or 3 as required for this decoding method!");
            else
                return;
        }

        // Set bitdepth of movie to 16 bpc for later texture creation from decoded video frames:
        movieRecordBANK[slotid].bitdepth = 16;

        // 1-layer: 1 bpc 16 gray pixel stored in 2/3 RGB8 pixel. This was achieved by multiplying
        // height by 2/3 in encoding, so invert by multiplying with 3/2:
        if (pixelFormat == 1) height = height * 3 / 2;

        // 3-layer: 1 RGB16 pixel stored in two adjacent RGB8 pixels. This was achieved by doubling
        // width, so undo by dividing width by 2:
        if (pixelFormat == 3) width = width / 2;

        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Playing back movie in Psychtoolbox proprietary 16 bpc %i channel encoding.\n", pixelFormat);
    }

    // Make sure the input format for raw Bayer sensor data is actually 1 layer grayscale, and that PTB for this OS supports debayering:
    if (specialFlags1 & 1024) {
        if (movieRecordBANK[slotid].pixelFormat != 1) {
            if (printErrors)
                PsychErrorExitMsg(PsychError_user, "specialFlags1 & 1024 set to indicate this movie consists of raw Bayer sensor data, but pixelFormat is not == 1, as required!");
            else
                return;
        }

        // Abort early if libdc1394 support is not available on this configuration:
        #ifndef PTBVIDEOCAPTURE_LIBDC
        if (printErrors)
            PsychErrorExitMsg(PsychError_user, "Sorry, Bayer filtering of video frames, as requested by specialFlags1 setting & 1024, is not supported on this operating system.");
        else
            return;
        #endif
    }

    // Define size of images in movie:
    movieRecordBANK[slotid].width = width;
    movieRecordBANK[slotid].height = height;

    // Assign parent window record, for use in movie deletion code:
    movieRecordBANK[slotid].parentRecord = win;

    // Should we dump the whole decoding pipeline graph to a file for visualization
    // with GraphViz? This can be controlled via PsychTweak('GStreamerDumpFilterGraph' dirname);
    if (getenv("GST_DEBUG_DUMP_DOT_DIR")) {
        // Dump complete decoding filter graph to a .dot file for later visualization with GraphViz:
        printf("PTB-DEBUG: Dumping movie decoder graph for movie %s to directory %s.\n", moviename, getenv("GST_DEBUG_DUMP_DOT_DIR"));
        GST_DEBUG_BIN_TO_DOT_FILE_WITH_TS(GST_BIN(movieRecordBANK[slotid].theMovie), GST_DEBUG_GRAPH_SHOW_ALL, "PsychMoviePlaybackGraph");
    }

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
    PsychWindowRecordType **windowRecordArray;
    int i, numWindows;

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
    if ((movieRecordBANK[moviehandle].parentRecord) && (movieRecordBANK[moviehandle].cached_texture > 0)) {
        // Yes. Release it.
        PsychSetGLContext(movieRecordBANK[moviehandle].parentRecord);
        glDeleteTextures(1, &(movieRecordBANK[moviehandle].cached_texture));
        movieRecordBANK[moviehandle].cached_texture = 0;
    }

    // Delete all references to us in textures originally originating from us:
    PsychCreateVolatileWindowRecordPointerList(&numWindows, &windowRecordArray);
    for(i = 0; i < numWindows; i++) {
        if ((windowRecordArray[i]->windowType == kPsychTexture) && (windowRecordArray[i]->texturecache_slot == moviehandle)) {
            // This one is referencing us. Reset its reference to "undefined" to detach it from us:
            windowRecordArray[i]->texturecache_slot = -1;
        }
    }
    PsychDestroyVolatileWindowRecordPointerList(windowRecordArray);

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
    GstElement      *theMovie;
    double          rate;
    double          targetdelta, realdelta, frames;
    GstBuffer       *videoBuffer = NULL;
    GstSample       *videoSample = NULL;
    gint64          bufferIndex;
    double          deltaT = 0;
    GstEvent        *event;
    static double   tStart = 0;
    double          tNow;
    double          preT, postT;
    unsigned char*  releaseMemPtr = NULL;
#if PSYCH_SYSTEM == PSYCH_WINDOWS
    #pragma warning( disable : 4068 )
#endif
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wmissing-field-initializers"
    GstMapInfo      mapinfo = GST_MAP_INFO_INIT;
    #pragma GCC diagnostic pop

    if (!PsychIsOnscreenWindow(win)) {
        PsychErrorExitMsg(PsychError_user, "Need onscreen window ptr!!!");
    }

    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided.");
    }

    if ((timeindex!=-1) && (timeindex < 0 || timeindex >= 100000.0)) {
        PsychErrorExitMsg(PsychError_user, "Invalid timeindex provided.");
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
    PsychGSProcessMovieContext(&(movieRecordBANK[moviehandle]), FALSE);

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

        // Allow context task to do its internal bookkeeping and cleanup work:
        PsychGSProcessMovieContext(&(movieRecordBANK[moviehandle]), FALSE);

        // Recheck:
        if (((0 != rate) && !movieRecordBANK[moviehandle].frameAvail) ||
            ((0 == rate) && !movieRecordBANK[moviehandle].preRollAvail)) {
            // Wait timed out after 0.5 secs.
            PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
            if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: No frame received after timed blocking wait of 0.5 seconds.\n");

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
        videoSample = gst_app_sink_pull_sample(GST_APP_SINK(movieRecordBANK[moviehandle].videosink));
    } else {
        // Passive fetch mode: Use prerolled buffers after seek:
        // These are available even after eos...

        // We can unlock early, thanks to videosink's internal buffering: XXX FIXME: Perfectly race-free to do this before the pull?
        PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);
        videoSample = gst_app_sink_pull_preroll(GST_APP_SINK(movieRecordBANK[moviehandle].videosink));
    }

    // Sample received?
    if (videoSample) {
        // Get pointer to buffer - no ownership transfer, no unref needed:
        videoBuffer = gst_sample_get_buffer(videoSample);

        // Assign pts presentation timestamp in pipeline stream time and convert to seconds:
        movieRecordBANK[moviehandle].pts = (double) GST_BUFFER_PTS(videoBuffer) / (double) 1e9;

        // Iff forward playback is active and a target timeindex was specified and this buffer is not at least of
        // that timeindex and at least one more buffer is queued, then skip this buffer, pull the next one and check
        // if that one meets the required pts:
        while ((rate > 0) && (timeindex >= 0) && (movieRecordBANK[moviehandle].pts < timeindex) && (movieRecordBANK[moviehandle].frameAvail > 0)) {
            // Tell user about reason for rejecting this buffer:
            if (PsychPrefStateGet_Verbosity() > 5) {
                printf("PTB-DEBUG: Fast-Skipped buffer id %i with pts %f secs < targetpts %f secs.\n", (int) GST_BUFFER_OFFSET(videoBuffer), movieRecordBANK[moviehandle].pts, timeindex);
            }

            // Decrement available frame counter:
            PsychLockMutex(&movieRecordBANK[moviehandle].mutex);
            movieRecordBANK[moviehandle].frameAvail--;
            PsychUnlockMutex(&movieRecordBANK[moviehandle].mutex);

            // Return the unused sample to queue:
            gst_sample_unref(videoSample);

            // Pull the next one. As frameAvail was > 0 at check-time, we know there is at least one pending,
            // so there shouldn't be a danger of hanging here:
            videoSample = gst_app_sink_pull_sample(GST_APP_SINK(movieRecordBANK[moviehandle].videosink));
            if (NULL == videoSample) {
                // This should never happen!
                printf("PTB-ERROR: No new video frame received in gst_app_sink_pull_sample skipper loop! Something's wrong. Aborting fetch.\n");
                return(FALSE);
            }

            // Get pointer to buffer - no ownership transfer, no unref needed:
            videoBuffer = gst_sample_get_buffer(videoSample);

            // Assign updated pts presentation timestamp of new candidate in pipeline stream time and convert to seconds:
            movieRecordBANK[moviehandle].pts = (double) GST_BUFFER_PTS(videoBuffer) / (double) 1e9;

            // Recheck if this is a better candidate...
        }

        // Compute timedelta and bufferindex:
        if (GST_CLOCK_TIME_IS_VALID(GST_BUFFER_DURATION(videoBuffer)))
            deltaT = (double) GST_BUFFER_DURATION(videoBuffer) / (double) 1e9;
        bufferIndex = GST_BUFFER_OFFSET(videoBuffer);

        if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: pts %f secs, dT %f secs, bufferId %i.\n", movieRecordBANK[moviehandle].pts, deltaT, (int) bufferIndex);

        // Assign pointer to videoBuffer's data directly:
        if (out_texture) {
            // Map the buffers memory for reading:
            if (!gst_buffer_map(videoBuffer, &mapinfo, GST_MAP_READ)) {
                printf("PTB-ERROR: Failed to map video data of movie frame! Something's wrong. Aborting fetch.\n");
                gst_sample_unref(videoSample);
                videoBuffer = NULL;
                return(FALSE);
            }
            out_texture->textureMemory = (GLuint*) mapinfo.data;
        }
    } else {
        printf("PTB-ERROR: No new video frame received in gst_app_sink_pull_sample! Something's wrong. Aborting fetch.\n");
        return(FALSE);
    }
    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: ...done.\n");

    PsychGetAdjustedPrecisionTimerSeconds(&tNow);
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Start of frame query to decode completion: %f msecs.\n", (tNow - tStart) * 1000.0);
    tStart = tNow;

    // Assign presentation_timestamp:
    if (presentation_timestamp) *presentation_timestamp = movieRecordBANK[moviehandle].pts;

    // Only create actual OpenGL texture if out_texture is non-NULL. Otherwise we're
    // just skipping this. Useful for benchmarks, fast forward seeking, etc.
    if (out_texture) {
        // Activate OpenGL context of target window:
        PsychSetGLContext(win);

        #if PSYCH_SYSTEM == PSYCH_OSX
        // Explicitely disable Apple's Client storage extensions. For now they are not really useful to us.
        glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, GL_FALSE);
        #endif

        // Build a standard PTB texture record:
        PsychMakeRect(out_texture->rect, 0, 0, movieRecordBANK[moviehandle].width, movieRecordBANK[moviehandle].height);

        // Set texture orientation as if it were an inverted Offscreen window: Upside-down.
        out_texture->textureOrientation = 3;

        // We use zero client storage memory bytes:
        out_texture->textureMemorySizeBytes = 0;

        // Assign default number of effective color channels:
        out_texture->nrchannels = movieRecordBANK[moviehandle].pixelFormat;

        // Is this grayscale movie actually a Bayer-encoded RGB movie? specialFlags1 & 1024 would indicate that:
        if (movieRecordBANK[moviehandle].specialFlags1 & 1024) {
            // This is Bayer raw sensor data which needs to get decoded into full RGB images.
            #ifdef PTBVIDEOCAPTURE_LIBDC
                // Ok, need to convert this grayscale image which actually contains raw Bayer sensor data into
                // a RGB image. Need to perform software Bayer filtering via libdc1394 Debayering routines.
                out_texture->textureMemory = (GLuint*) PsychDCDebayerFrame((unsigned char*) (out_texture->textureMemory), movieRecordBANK[moviehandle].width, movieRecordBANK[moviehandle].height, movieRecordBANK[moviehandle].bitdepth);

                // Return failure if Debayering did not work:
                if (out_texture->textureMemory == NULL) {
                    gst_buffer_unmap(videoBuffer, &mapinfo);
                    gst_sample_unref(videoSample);
                    videoBuffer = NULL;
                    return(FALSE);
                }

                releaseMemPtr = (unsigned char*) out_texture->textureMemory;
                out_texture->nrchannels = 3; // Always 3 for RGB.
            #else
                // Won't ever reach this, as already Screen('OpenMovie') would have bailed out
                // if libdc1394 is not supported.
                return(FALSE);
            #endif
        }

        // Assign default depth according to number of channels:
        out_texture->depth = out_texture->nrchannels * movieRecordBANK[moviehandle].bitdepth;

        if (out_texture->nrchannels < 4) {
            // For 1-3 channel textures, play safe, don't assume alignment:
            out_texture->textureByteAligned = 1;
        }
        else {
            // 4 channel format:
            // Textures are aligned on at least 4 Byte boundaries because texels are RGBA8. For
            // frames of even-numbered pixel width, we can even get 8 Byte alignment:
            out_texture->textureByteAligned = (movieRecordBANK[moviehandle].width % 2) ? 4 : 8;
        }

        // Assign texturehandle of our cached texture, if any, so it gets recycled now:
        out_texture->textureNumber = movieRecordBANK[moviehandle].cached_texture;

        // Mark this texture as originating from us, ie., our moviehandle, so texture recycling
        // actually gets used:
        out_texture->texturecache_slot = moviehandle;

        // YUV 422 packed pixel upload requested?
        if ((win->gfxcaps & kPsychGfxCapUYVYTexture) && (movieRecordBANK[moviehandle].pixelFormat == 5)) {
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
            out_texture->textureexternaltype = GL_UNSIGNED_SHORT_8_8_MESA;

            // Number of effective channels is 3 for RGB8:
            out_texture->nrchannels = 3;

            // And 24 bpp depth:
            out_texture->depth = 24;

            // Byte alignment: For even number of pixels, assume at least 4 Byte alignment due to packing of 2 effective
            // pixels into one 32-Bit packet, maybe even 8 Byte alignment if divideable by 4. For other width's, assume
            // no alignment ie., 1 Byte:
            out_texture->textureByteAligned = (movieRecordBANK[moviehandle].width % 2) ? 1 : ((movieRecordBANK[moviehandle].width % 4) ? 4 : 8);
        }

        // Upload of a "pseudo YUV" planar texture with only 8 bits Y component requested?
        if ((movieRecordBANK[moviehandle].pixelFormat == 7) || (movieRecordBANK[moviehandle].pixelFormat == 8)) {
            // We encode Y luminance data inside a 8 bit per pixel luminance texture. The
            // "Y" luminance plane is stored at full 1 sample per pixel resolution with 8 bits.
            // As such the texture appears to OpenGL as a normal LUMINANCE8 texture. Conversion of the Y
            // luminance data into useable RGBA8 pixel fragments will happen during rendering via a suitable fragment
            // shader. The net gain of this is that we can skip any kind of cpu based colorspace conversion
            // for video formats/codecs which provide YUV data, by offloading the conversion to the GPU:
            out_texture->textureinternalformat = 0;

            // Mark texture as planar encoded, so proper conversion shader gets applied during
            // call to PsychNormalizeTextureOrientation(), prior to any render-to-texture operation, e.g.,
            // if used as an offscreen window, or as a participant of a Screen('TransformTexture') call:
            out_texture->specialflags |= kPsychPlanarTexture;

            // Assign special filter shader for Y8 -> RGBA8 color-space conversion of the
            // planar texture during drawing or PsychNormalizeTextureOrientation():
            if (!PsychAssignPlanarI800TextureShader(out_texture, win)) PsychErrorExitMsg(PsychError_user, "Assignment of Y8-Y800 video decoding shader failed during movie texture creation!");

            // Number of effective channels is 1 for L8:
            out_texture->nrchannels = 1;

            // And 8 bpp depth: This will trigger bog-standard LUMINANCE8 texture creation in PsychCreateTexture():
            out_texture->depth = 8;

            // Byte alignment - Only depends on width of an image row, given the 1 Byte per pixel data:
            out_texture->textureByteAligned = 1;
            if (movieRecordBANK[moviehandle].width % 2 == 0) out_texture->textureByteAligned = 2;
            if (movieRecordBANK[moviehandle].width % 4 == 0) out_texture->textureByteAligned = 4;
            if (movieRecordBANK[moviehandle].width % 8 == 0) out_texture->textureByteAligned = 8;
        }

        // YUV I420 planar pixel upload requested?
        if (movieRecordBANK[moviehandle].pixelFormat == 6) {
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
            PsychMakeRect(out_texture->rect, 0, 0, movieRecordBANK[moviehandle].width, movieRecordBANK[moviehandle].height * 1.5);

            // Check if 1.5x height texture fits within hardware limits of this GPU:
            if (movieRecordBANK[moviehandle].height * 1.5 > win->maxTextureSize) PsychErrorExitMsg(PsychError_user, "Videoframe size too big for this graphics card and pixelFormat! Please retry with a pixelFormat of 4 in 'OpenMovie'.");

            // Byte alignment: Assume no alignment for now:
            out_texture->textureByteAligned = 1;

            // Create planar "I420 inside L8" texture:
            PsychCreateTexture(out_texture);

            // Restore rect and clientrect of texture to effective size of video frame:
            PsychMakeRect(out_texture->rect, 0, 0, movieRecordBANK[moviehandle].width, movieRecordBANK[moviehandle].height);
            PsychCopyRect(out_texture->clientrect, out_texture->rect);

            // Mark texture as planar encoded, so proper conversion shader gets applied during
            // call to PsychNormalizeTextureOrientation(), prior to any render-to-texture operation, e.g.,
            // if used as an offscreen window, or as a participant of a Screen('TransformTexture') call:
            out_texture->specialflags |= kPsychPlanarTexture;

            // Assign special filter shader for sampling and color-space conversion of the
            // planar texture during drawing or PsychNormalizeTextureOrientation():
            if (!PsychAssignPlanarI420TextureShader(out_texture, win)) PsychErrorExitMsg(PsychError_user, "Assignment of I420 video decoding shader failed during movie texture creation!");

            // Number of effective channels is 3 for RGB8:
            out_texture->nrchannels = 3;

            // And 24 bpp depth:
            out_texture->depth = 24;
        }
        else if (movieRecordBANK[moviehandle].bitdepth > 8) {
            // Is this a > 8 bpc image format? If not, we ain't nothing more to prepare.
            // If yes, we need to use a high precision floating point texture to represent
            // the > 8 bpc image payload without loss of image information:

            // highbitthreshold: If the net bpc value is greater than this, then use 32bpc floats
            // instead of 16 bpc half-floats, because 16 bpc would not be sufficient to represent
            // more than highbitthreshold bits faithfully:
            const int highbitthreshold = 11;
            unsigned int w = movieRecordBANK[moviehandle].width;

            // 9 - 16 bpc color/luminance resolution:
            out_texture->depth = out_texture->nrchannels * ((movieRecordBANK[moviehandle].bitdepth > highbitthreshold) ? 32 : 16);

            // Byte alignment: Assume at least 2 Byte alignment due to 16 bit per component aka 2 Byte input:
            out_texture->textureByteAligned = 2;

            if (out_texture->nrchannels == 1) {
                // 1 layer Luminance:
                out_texture->textureinternalformat = (movieRecordBANK[moviehandle].bitdepth > highbitthreshold) ? GL_LUMINANCE_FLOAT32_APPLE : GL_LUMINANCE_FLOAT16_APPLE;
                out_texture->textureexternalformat = GL_LUMINANCE;
                // Override for missing floating point texture support: Try to use 16 bit fixed point signed normalized textures [-1.0 ; 1.0] resolved at 15 bits:
                if (!(win->gfxcaps & kPsychGfxCapFPTex16)) out_texture->textureinternalformat = GL_LUMINANCE16_SNORM;
                out_texture->textureByteAligned = (w % 2) ? 2 : ((w % 4) ? 4 : 8);
            }
            else if (out_texture->nrchannels == 3) {
                // 3 layer RGB:
                out_texture->textureinternalformat = (movieRecordBANK[moviehandle].bitdepth > highbitthreshold) ? GL_RGB_FLOAT32_APPLE : GL_RGB_FLOAT16_APPLE;
                out_texture->textureexternalformat = GL_RGB;
                // Override for missing floating point texture support: Try to use 16 bit fixed point signed normalized textures [-1.0 ; 1.0] resolved at 15 bits:
                if (!(win->gfxcaps & kPsychGfxCapFPTex16)) out_texture->textureinternalformat = GL_RGB16_SNORM;
                out_texture->textureByteAligned = (w % 2) ? 2 : ((w % 4) ? 4 : 8);
            }
            else {
                // 4 layer RGBA:
                out_texture->textureinternalformat = (movieRecordBANK[moviehandle].bitdepth > highbitthreshold) ? GL_RGBA_FLOAT32_APPLE : GL_RGBA_FLOAT16_APPLE;
                out_texture->textureexternalformat = GL_RGBA;
                // Override for missing floating point texture support: Try to use 16 bit fixed point signed normalized textures [-1.0 ; 1.0] resolved at 15 bits:
                if (!(win->gfxcaps & kPsychGfxCapFPTex16)) out_texture->textureinternalformat = GL_RGBA16_SNORM;
                // Always 8 Byte aligned:
                out_texture->textureByteAligned = 8;
            }

            // External datatype is 16 bit unsigned integer, each color component encoded in a 16 bit value:
            out_texture->textureexternaltype = GL_UNSIGNED_SHORT;

            // Scale input data, so highest significant bit of payload is in bit 16:
            glPixelTransferi(GL_RED_SCALE,   1 << (16 - movieRecordBANK[moviehandle].bitdepth));
            glPixelTransferi(GL_GREEN_SCALE, 1 << (16 - movieRecordBANK[moviehandle].bitdepth));
            glPixelTransferi(GL_BLUE_SCALE,  1 << (16 - movieRecordBANK[moviehandle].bitdepth));

            // Let PsychCreateTexture() do the rest of the job of creating, setting up and
            // filling an OpenGL texture with content:
            PsychCreateTexture(out_texture);

            // Undo scaling:
            glPixelTransferi(GL_RED_SCALE, 1);
            glPixelTransferi(GL_GREEN_SCALE, 1);
            glPixelTransferi(GL_BLUE_SCALE, 1);
        }
        else {
            // Let PsychCreateTexture() do the rest of the job of creating, setting up and
            // filling an OpenGL texture with content:
            PsychCreateTexture(out_texture);
        }

        // Release buffer for target RGB debayered image, if any:
        if ((movieRecordBANK[moviehandle].specialFlags1 & 1024) && releaseMemPtr) free(releaseMemPtr);

        // NULL-out the texture memory pointer after PsychCreateTexture(). This is not strictly
        // needed, as PsychCreateTexture() did it already, but we add it here as an annotation
        // to make it obvious during code correctness review that we won't touch or free() the
        // video memory buffer anymore, which is owned and only memory-managed by GStreamer:
        out_texture->textureMemory = NULL;

        // After PsychCreateTexture() the cached texture object from our cache is used
        // and no longer available for recycling. We mark the cache as empty:
        // It will be filled with a new textureid for recycling if a texture gets
        // deleted in PsychMovieDeleteTexture()....
        movieRecordBANK[moviehandle].cached_texture = 0;

        // Does usercode want immediate conversion of texture into standard RGBA8 packed pixel
        // upright format for use as a render-target? If so, do it:
        if (movieRecordBANK[moviehandle].specialFlags1 & 16) {
            // Transform out_texture video texture into a normalized, upright texture if it isn't already in
            // that format. We require this standard orientation for simplified shader design.
            PsychSetShader(win, 0);
            PsychNormalizeTextureOrientation(out_texture);
        }

        PsychGetAdjustedPrecisionTimerSeconds(&tNow);
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Decode completion to texture created: %f msecs.\n", (tNow - tStart) * 1000.0);
        tStart = tNow;

        // End of texture creation code.
    }

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
    gst_buffer_unmap(videoBuffer, &mapinfo);
    gst_sample_unref(videoSample);
    videoBuffer = NULL;

    // Manually advance movie time, if in fetch mode:
    if (0 == rate) {
        // We are in manual fetch mode: Need to manually advance movie to next
        // media sample:
        movieRecordBANK[moviehandle].endOfFetch = 0;
        preT = PsychGSGetMovieTimeIndex(moviehandle);
        event = gst_event_new_step(GST_FORMAT_BUFFERS, 1, 1.0, TRUE, FALSE);
        // Send the seek event *only* to the videosink. This follows recommendations from GStreamer SDK tutorial 13 (Playback speed) to
        // not send to high level playbin itself, as that would propagate to all sinks and trigger multiple seeks. While this was not
        // ever a problem in the past on Linux or with upstream GStreamer, it caused deadlocks, timeouts and seek failures when done
        // with the GStreamer SDK on some movie files that have audio tracks, e.g., our standard demo movie! Sending only to videosink
        // fixes this problem:
        if (!gst_element_send_event(movieRecordBANK[moviehandle].videosink, event)) printf("PTB-DEBUG: In single-step seek I - Failed.\n");

        // Block until seek completed, failed, or timeout of 10 seconds reached:
        if (GST_STATE_CHANGE_SUCCESS != gst_element_get_state(theMovie, NULL, NULL, (GstClockTime) (10 * 1e9))) printf("PTB-DEBUG: In single-step seek II - Failed.\n");
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
 *  PsychGSFreeMovieTexture() - Release texture memory for a texture.
 *
 *  This routine is called by PsychDeleteTexture() in PsychTextureSupport.c
 *  It performs the special cleanup necessary for cached movie textures.
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
    GstElement *audiosink, *actual_audiosink;
    gchar* pstring;
    int    dropped = 0;
    GstElement *theMovie = NULL;
    double timeindex;
    GstSeekFlags seekFlags = 0;

    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }

    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }

    // Try to set movie playback rate to value identical to current value?
    if (playbackrate == movieRecordBANK[moviehandle].rate) {
        // Yes: This would be a no-op, except we allow to change the sound output volume
        // dynamically and on-the-fly with low overhead this way:

        // Set volume and mute state for audio:
        g_object_set(G_OBJECT(theMovie), "mute", (soundvolume <= 0) ? TRUE : FALSE, NULL);
        g_object_set(G_OBJECT(theMovie), "volume", soundvolume, NULL);

        // Done. Return success status code:
        return(0);
    }

    if (playbackrate != 0) {
        // Start playback of movie:

        // Set volume and mute state for audio:
        g_object_set(G_OBJECT(theMovie), "mute", (soundvolume <= 0) ? TRUE : FALSE, NULL);
        g_object_set(G_OBJECT(theMovie), "volume", soundvolume, NULL);

        // Set playback rate: An explicit seek to the position we are already (supposed to be)
        // is needed to avoid jumps in movies with bad encoding or keyframe placement:
        timeindex = PsychGSGetMovieTimeIndex(moviehandle);

        // Which loop setting?
        if (loop <= 0) {
            // Looped playback disabled. Set to well defined off value zero:
            loop = 0;
        }
        else {
            // Looped playback requested. With default settings (==1)?
            // Otherwise we'll just pass on any special != 1 setting as a
            // user-override:
            if (loop == 1) {
                // Playback with defaults. Apply default setup + specialFlags1 quirks:

                // specialFlags & 32? Use 'uri' injection method for looped playback, instead of seek method:
                if (movieRecordBANK[moviehandle].specialFlags1 &  32) loop = 2;

                // specialFlags & 64? Use segment seeks.
                if (movieRecordBANK[moviehandle].specialFlags1 &  64) loop |= 4;

                // specialFlags & 128? Use pipeline flushing seeks
                if (movieRecordBANK[moviehandle].specialFlags1 & 128) loop |= 8;
            }
        }

        // On some movies and configurations, we need a segment seek as indicated by flag 0x4:
        if (loop & 0x4) seekFlags |= GST_SEEK_FLAG_SEGMENT;

        if (playbackrate > 0) {
            gst_element_seek(theMovie, playbackrate, GST_FORMAT_TIME, GST_SEEK_FLAG_FLUSH | GST_SEEK_FLAG_ACCURATE | ((loop & 0x1) ? seekFlags : 0), GST_SEEK_TYPE_SET,
                             (gint64) (timeindex * (double) 1e9), GST_SEEK_TYPE_NONE, GST_CLOCK_TIME_NONE);
        }
        else {
            gst_element_seek(theMovie, playbackrate, GST_FORMAT_TIME, GST_SEEK_FLAG_FLUSH | GST_SEEK_FLAG_ACCURATE | ((loop & 0x1) ? seekFlags : 0), GST_SEEK_TYPE_SET,
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
            PsychGSProcessMovieContext(&(movieRecordBANK[moviehandle]), FALSE);
        }
    }
    else {
        // Stop playback of movie:
        movieRecordBANK[moviehandle].rate = 0;
        movieRecordBANK[moviehandle].startPending = 0;
        movieRecordBANK[moviehandle].loopflag = 0;
        movieRecordBANK[moviehandle].endOfFetch = 0;

        // Print name of audio sink - the output device which was actually playing the sound, if requested:
        // This is a Linux only feature, as GStreamer for MS-Windows doesn't support such queries at all,
        // and GStreamer for OSX doesn't expose the information in a way that would be in any way meaningful for us.
        if ((PSYCH_SYSTEM == PSYCH_LINUX) && (PsychPrefStateGet_Verbosity() > 3)) {
            audiosink = NULL;
            if (g_object_class_find_property(G_OBJECT_GET_CLASS(theMovie), "audio-sink")) {
                g_object_get(G_OBJECT(theMovie), "audio-sink", &audiosink, NULL);
            }

            if (audiosink) {
                actual_audiosink = NULL;
                actual_audiosink = (GST_IS_CHILD_PROXY(audiosink)) ? ((GstElement*) gst_child_proxy_get_child_by_index(GST_CHILD_PROXY(audiosink), 0)) : audiosink;
                if (actual_audiosink) {
                    if (g_object_class_find_property(G_OBJECT_GET_CLASS(actual_audiosink), "device")) {
                        pstring = NULL;
                        g_object_get(G_OBJECT(actual_audiosink), "device", &pstring, NULL);
                        if (pstring) {
                            printf("PTB-INFO: Audio output device name for movie playback was '%s'", pstring);
                            g_free(pstring); pstring = NULL;
                        }
                    }

                    if (g_object_class_find_property(G_OBJECT_GET_CLASS(actual_audiosink), "device-name")) {
                        pstring = NULL;
                        g_object_get(G_OBJECT(actual_audiosink), "device-name", &pstring, NULL);
                        if (pstring) {
                            printf(" [%s].", pstring);
                            g_free(pstring); pstring = NULL;
                        }
                    }

                    printf("\n");
                    if (actual_audiosink != audiosink) gst_object_unref(actual_audiosink);
                }
                gst_object_unref(audiosink);
            }
        }

        PsychMoviePipelineSetState(theMovie, GST_STATE_PAUSED, 10.0);
        PsychGSProcessMovieContext(&(movieRecordBANK[moviehandle]), FALSE);

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
 *  do final cleanup. It releases all movie objects.
 *
 */
void PsychGSExitMovies(void)
{
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
    GstElement  *theMovie = NULL;
    gint64      pos_nsecs;

    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }

    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }

    if (!gst_element_query_position(theMovie, GST_FORMAT_TIME, &pos_nsecs)) {
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
    GstElement      *theMovie;
    double          oldtime;
    gint64          targetIndex;
    GstSeekFlags    flags;

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

    // Need segment seek flag for seek during active looped playback if also flag 0x4 is set:
    if ((movieRecordBANK[moviehandle].rate != 0) && (movieRecordBANK[moviehandle].loopflag & 0x1) && (movieRecordBANK[moviehandle].loopflag & 0x4)) {
        flags |= GST_SEEK_FLAG_SEGMENT;
    }

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
    if (GST_STATE_CHANGE_FAILURE == gst_element_get_state(theMovie, NULL, NULL, (GstClockTime) (30 * 1e9)) && (PsychPrefStateGet_Verbosity() > 1)) {
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

// #if GST_CHECK_VERSION(1,0,0)
#endif
// #ifdef PTB_USE_GSTREAMER
#endif
