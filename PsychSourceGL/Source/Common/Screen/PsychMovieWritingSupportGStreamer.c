/*
	Common/Screen/PsychMovieWritingSupportGStreamer.c

	PLATFORMS:

		GNU/Linux and MS-Windows with Octave or Matlab R2007a and later for now.
		Will be supported on OS/X as well in the future.

		This is the movie editing and writing/creation engine based on the
		GStreamer multimedia framework.

	AUTHORS:

		Mario Kleiner           mk              mario.kleiner@tuebingen.mpg.de

	HISTORY:

		06/06/11		mk		Wrote it.

	DESCRIPTION:

		Psychtoolbox functions for dealing with GStreamer movie editing.

	NOTES:

*/

#include "Screen.h"

// GStreamer support enabled? If so, we use GStreamer based movie writing.
#ifdef PTB_USE_GSTREAMER

// GStreamer includes:
#include <gst/gst.h>

// PsychGetCodecLaunchLineFromString() - Helper function for GStreamer based movie writing.
// Defined in PsychVideoCaptureSupport.h: psych_bool PsychGetCodecLaunchLineFromString(char* codecSpec, char* launchString);

// GStreamer implementation of movie writing support:

// Record which defines all state for a capture device:
typedef struct {
	volatile psych_bool                             eos;
	GMainLoop*                                      Context;
	GstElement*                                     Movie;
	GstElement*                                     ptbvideoappsrc;
	GstElement*                                     ptbaudioappsrc;
	GstBus*                                         bus;
	GstBuffer*                                      PixMap;
	guint32                                         CodecType;
	char                                            File[FILENAME_MAX];
	int                                             height;
	int                                             width;
} PsychMovieWriterRecordType;

static PsychMovieWriterRecordType moviewriterRecordBANK[PSYCH_MAX_MOVIEWRITERDEVICES];
static int moviewritercount = 0;
static psych_bool firsttime = TRUE;

/* Perform one context loop iteration (for bus message handling) if doWait == false,
 * or two seconds worth of iterations if doWait == true. This drives the message-bus
 * callback, so needs to be performed to get any error reporting etc.
 */
static int PsychGSProcessMovieContext(GMainLoop *loop, psych_bool doWait)
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

void PsychMovieWritingInit(void)
{
	int i;
	
	for (i = 0; i < PSYCH_MAX_MOVIEWRITERDEVICES; i++) {
		memset(&(moviewriterRecordBANK[i]), 0, sizeof(PsychMovieWriterRecordType));
	}
	
	moviewritercount = 0;
	firsttime = TRUE;
	return;
}

void PsychDeleteAllMovieWriters(void)
{
	int i;
		
	for (i = 0; i < PSYCH_MAX_MOVIEWRITERDEVICES; i++) {
		if (moviewriterRecordBANK[i].Movie) PsychFinalizeNewMovieFile(i);
	}
}

void PsychExitMovieWriting(void)
{
	if (firsttime) return;

	PsychDeleteAllMovieWriters();
	firsttime = TRUE;
	return;
}

PsychMovieWriterRecordType* PsychGetMovieWriter(int moviehandle, psych_bool unsafe)
{
	if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIEWRITERDEVICES) PsychErrorExitMsg(PsychError_user, "Invalid handle for moviewriter provided!");
	if (!unsafe && (NULL == moviewriterRecordBANK[moviehandle].Movie)) PsychErrorExitMsg(PsychError_user, "Invalid handle for moviewriter provided! No such writer open.");
	return(&(moviewriterRecordBANK[moviehandle]));
}

unsigned char*	PsychGetVideoFrameForMoviePtr(int moviehandle, unsigned int* twidth, unsigned int* theight)
{
	PsychMovieWriterRecordType* pwriterRec = PsychGetMovieWriter(moviehandle, FALSE);

	// Buffer already created?
	if (NULL == pwriterRec->PixMap) {
		// No. Let's create a suitable one:
		pwriterRec->PixMap = gst_buffer_try_new_and_alloc(pwriterRec->width * pwriterRec->height * 4 * 1);

		// Out of memory condition!
		if (NULL == pwriterRec->PixMap) return(NULL);
	}

	// Double-check:
	if (NULL == GST_BUFFER_DATA(pwriterRec->PixMap)) return(NULL);

	*twidth  = pwriterRec->width;
	*theight = pwriterRec->height;
	return((unsigned char*) GST_BUFFER_DATA(pwriterRec->PixMap));
}

int PsychAddVideoFrameToMovie(int moviehandle, int frameDurationUnits, psych_bool isUpsideDown)
{
	PsychMovieWriterRecordType* pwriterRec = PsychGetMovieWriter(moviehandle, FALSE);
    GstBuffer*          refBuffer = NULL;
    GstBuffer*          curBuffer = NULL;
	GstFlowReturn       ret;
	int                 x, y, w, h;
	unsigned char*		pixptr;
	unsigned int*		wordptr;
	unsigned int		*wordptr2, *wordptr1;
	unsigned int		dummy;
    int                 bframeDurationUnits = frameDurationUnits;
    
	if (NULL == pwriterRec->ptbvideoappsrc) return(0);
	if (NULL == pwriterRec->PixMap) return(0);

	if ((frameDurationUnits < 1) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In AddFrameToMovie: Negative or zero 'frameduration' %i units for moviehandle %i provided! Sounds like trouble ahead.\n", frameDurationUnits, moviehandle);

    pixptr   = (unsigned char*) GST_BUFFER_DATA(pwriterRec->PixMap);
    wordptr  = (unsigned int*)  GST_BUFFER_DATA(pwriterRec->PixMap);

	// Draw testpattern: Disabled at compile-time by default:
	if (FALSE) {
		for (y = 0; y < pwriterRec->height; y++) {
			for (x = 0; x < pwriterRec->width; x++) {
				*(pixptr++) = (unsigned char) 255; // alpha
				*(pixptr++) = (unsigned char) y; // Red
				*(pixptr++) = (unsigned char) x; // Green
				*(pixptr++) = (unsigned char) 0; // Blue
			}
		}
	}
	
	// Imagebuffer is upside-down: Need to flip it vertically:
	if (isUpsideDown) {
		h = pwriterRec->height;
		w = pwriterRec->width;
		wordptr1 = wordptr;
		for (y = 0; y < h/2; y++) {
			wordptr2 = wordptr;
			wordptr2 += ((h - 1 - y) * w);
			for (x = 0; x < w; x++) {
				dummy = *wordptr1;
				*(wordptr1++) = *wordptr2;
				*(wordptr2++) = dummy;				 
			}
		}
	}

    // Make backup copy of buffer for replication if needed:
    if (frameDurationUnits > 1) {
        refBuffer = gst_buffer_copy(pwriterRec->PixMap);
    }

	// Add encoded buffer to movie:
	g_signal_emit_by_name(pwriterRec->ptbvideoappsrc, "push-buffer", pwriterRec->PixMap, &ret);
    
	// Unref it - it is now owned and memory managed by the pipeline:
	gst_buffer_unref(pwriterRec->PixMap);

	// Drop our handle to it, so we can allocate a new one on demand:
	pwriterRec->PixMap = NULL;

    // A dumb way to let this buffer last frameDurationUnits > 1 instead of
    // the default typical frameDurationUnits == 1. We simply create and
    // push frameDurationUnits-1 extra identicaly copies of the buffer. This
    // will not win the "computational efficiency award 2011", but should be robust.
    if (frameDurationUnits > 1) {
        // One already done:
        frameDurationUnits--;
        
        // Repeat for remaining ones:
        while ((frameDurationUnits > 0) && (ret == GST_FLOW_OK)) {
            // Create a new copy:
            curBuffer = gst_buffer_copy(refBuffer);
            
            // Add copied buffer to movie:
            g_signal_emit_by_name(pwriterRec->ptbvideoappsrc, "push-buffer", curBuffer, &ret);
    
            // Unref it - it is now owned and memory managed by the pipeline:
            gst_buffer_unref(curBuffer);
            curBuffer = NULL;

            // One less...
            frameDurationUnits--;
        }
    }
    
    if (refBuffer) gst_buffer_unref(refBuffer);
    refBuffer = NULL;

	if (ret != GST_FLOW_OK) {
		// Oopsie! Error encountered - Abort.
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In AddFrameToMovie: Adding current frame to moviehandle %i failed [push-buffer returned error code %i]!\n", moviehandle, (int) ret);
		return((int) ret);
	}

	PsychGSProcessMovieContext(pwriterRec->Context, FALSE);

	if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG:In AddFrameToMovie: Added new videoframe with %i units duration and upsidedown = %i to moviehandle %i.\n", bframeDurationUnits, (int) isUpsideDown, moviehandle);

	// Return success:
	return((int) ret);
}

psych_bool PsychAddAudioBufferToMovie(int moviehandle, unsigned int nrChannels, unsigned int nrSamples, double* buffer)
{
	PsychMovieWriterRecordType* pwriterRec = PsychGetMovieWriter(moviehandle, FALSE);

	GstFlowReturn       ret;
	float*              fwordptr;
	float               v;
	unsigned int        n, i;
	GstBuffer*          pushBuffer;
    
	// Child protection: Audio writing enabled for this movie?
	if (NULL == pwriterRec->ptbaudioappsrc) {
		PsychErrorExitMsg(PsychError_user, "Tried to add audio data to a movie which was created without an audio track.");
	}
    
	// nrChannels and nrSamples are already validated by high level code.
	// Just calculate total sample count and required buffer size:
	n = nrChannels * nrSamples;
    
	// Create GstBuffer for audio data:
	pushBuffer = gst_buffer_try_new_and_alloc(n * sizeof(float));

	// Out of memory condition!
	if (NULL == pushBuffer) {
		PsychErrorExitMsg(PsychError_outofMemory, "Out of memory when trying to add audio data to movie! (Part I)");
		return(FALSE);
	}

	// Double-check:
	if (NULL == GST_BUFFER_DATA(pushBuffer)) {
		PsychErrorExitMsg(PsychError_outofMemory, "Out of memory when trying to add audio data to movie! (Part II)");
		return(FALSE);
	}
    
	// Convert and copy sample data:
	fwordptr = (float*) GST_BUFFER_DATA(pushBuffer);
	for (i = 0; i < n; i++) {
		// Fetch and convert from double to float:
		v = (float) *(buffer++);;

		// Clip:
		if (v < -1.0) v = -1.0;
		if (v > +1.0) v = +1.0;

		// Push to float buffer:
		*(fwordptr++) = v; 
	}
    
	// Add encoded buffer to movie:
	g_signal_emit_by_name(pwriterRec->ptbaudioappsrc, "push-buffer", pushBuffer, &ret);

	// Unref it - it is now owned and memory managed by the pipeline:
	gst_buffer_unref(pushBuffer);
	pushBuffer = NULL;

	if (ret != GST_FLOW_OK) {
		// Oopsie! Error encountered - Abort.
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In AddAudioBufferToMovie: Adding current audio buffer to moviehandle %i failed [push-buffer returned error code %i]!\n", moviehandle, (int) ret);
		return(FALSE);
	}

	// Do a bit of event processing for handling of potential GStreamer messages:
	PsychGSProcessMovieContext(pwriterRec->Context, FALSE);

	if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG:In AddAudioBufferToMovie: Added new audio buffer to moviehandle %i.\n", moviehandle);

	// Return success:
	return(TRUE);
}


/* Initiate pipeline state changes: Startup, Preroll, Playback, Pause, Standby, Shutdown. */
static psych_bool PsychMoviePipelineSetState(GstElement* camera, GstState state, double timeoutSecs)
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
		printf("PTB-INFO: Statechange completed with GST_STATE_CHANGE_NO_PREROLL.\n");
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

/* Receive messages from the pipeline message bus and handle them: */
static gboolean PsychMovieBusCallback(GstBus *bus, GstMessage *msg, gpointer dataptr)
{
  PsychMovieWriterRecordType* dev = (PsychMovieWriterRecordType*) dataptr;

  switch (GST_MESSAGE_TYPE (msg)) {
    case GST_MESSAGE_EOS:
	if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Videobus: Message EOS received.\n");
	dev->eos = TRUE;
    break;

    case GST_MESSAGE_WARNING: {
      gchar  *debug;
      GError *error;

      gst_message_parse_warning(msg, &error, &debug);
      if (PsychPrefStateGet_Verbosity() > 3) { 
	      printf("PTB-WARNING: GStreamer movie writing engine reports this warning:\n"
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
	      printf("PTB-ERROR: GStreamer movie writing engine reports this error:\n"
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

int PsychCreateNewMovieFile(char* moviefile, int width, int height, double framerate, char* movieoptions)
{
	PsychMovieWriterRecordType*             pwriterRec = NULL;
	int                                     moviehandle = 0;
	GError                                  *myErr = NULL;
	char*                                   poption;
	char                                    codecString[1000];
	char                                    launchString[10000];
	int                                     dummyInt;
	float                                   dummyFloat;
	char                                    myfourcc[5];
	psych_bool                              doAudio = FALSE;

	// Still capacity left?
	if (moviewritercount >= PSYCH_MAX_MOVIEWRITERDEVICES) PsychErrorExitMsg(PsychError_user, "Maximum number of movie writers exceeded. Please close some first!");

	// Find first free (i.e., NULL) slot and assign moviehandle:
	while ((pwriterRec = PsychGetMovieWriter(moviehandle, TRUE)) && pwriterRec->Movie) moviehandle++;

	if (firsttime) {
		// Make sure GStreamer is ready:
		PsychGSCheckInit("movie writing");
		firsttime = FALSE;
	}

	// Store movie filename:
	strcpy(pwriterRec->File, moviefile);

	// Store width, height:
	pwriterRec->height  = height;
	pwriterRec->width   = width;
	pwriterRec->eos     = FALSE;

	// If no movieoptions specified, create default string for default
	// codec selection and configuration:
	if (strlen(movieoptions) == 0) {
		// No options provided. Select default encoder with default settings:
		movieoptions = strdup("DEFAULTenc");
	} else if ((poption = strstr(movieoptions, ":CodecSettings="))) {
		// Replace ':' with a zero in movieoptions, so it gets null-terminated:
		movieoptions = poption;
		*movieoptions = 0;

		// Move after null-terminator:
		movieoptions++;

		// Replace the ':CodecSettings=' with the special keyword 'DEFAULTenc', so
		// so the default video codec is chosen, but the given settings override its
		// default parameters.
		strncpy(movieoptions, "DEFAULTenc    ", strlen("DEFAULTenc    "));

		if (strlen(movieoptions) == 0) PsychErrorExitMsg(PsychError_user, "Invalid (empty) :CodecSettings= parameter specified. Aborted.");
	} else if ((poption = strstr(movieoptions, ":CodecType="))) {
		// Replace ':' with a zero in movieoptions, so it gets null-terminated
		// and only points to the actual movie filename:
		movieoptions = poption;
		*movieoptions = 0;

		// Advance movieoptions to point to the actual codec spec string:
		movieoptions+= 11;

		if (strlen(movieoptions) == 0) PsychErrorExitMsg(PsychError_user, "Invalid (empty) :CodecType= parameter specified. Aborted.");
	}

	// Assign numeric 32-bit FOURCC equivalent code to select codec:
	// This is optional. We default to kH264CodecType:
	if ((poption = strstr(movieoptions, "CodecFOURCCId="))) {
		if (sscanf(poption, "CodecFOURCCId=%i", &dummyInt) == 1) {
			pwriterRec->CodecType = dummyInt;
			if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Codec with FOURCC numeric id %i [%" GST_FOURCC_FORMAT "] requested for encoding of movie %i [%s].\n", dummyInt, GST_FOURCC_ARGS(dummyInt), moviehandle, moviefile);
			if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Codec selection by FOURCC not yet supported. FOURCC code ignored!\n");            
		}
		else PsychErrorExitMsg(PsychError_user, "Invalid CodecFOURCCId= parameter provided in movieoptions parameter. Parse error!");
	}

	// Assign 4 character string FOURCC code to select codec:
	if ((poption = strstr(movieoptions, "CodecFOURCC="))) {
		if (sscanf(poption, "CodecFOURCC=%c%c%c%c", &myfourcc[0], &myfourcc[1], &myfourcc[2], &myfourcc[3]) == 4) {
			myfourcc[4] = 0;
			dummyInt = (int) GST_STR_FOURCC (myfourcc);
			pwriterRec->CodecType = dummyInt;
			if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Codec with FOURCC numeric id %i [%" GST_FOURCC_FORMAT "] requested for encoding of movie %i [%s].\n", dummyInt, GST_FOURCC_ARGS(dummyInt), moviehandle, moviefile);
			if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Codec selection by FOURCC not yet supported. FOURCC code ignored!\n");            
		}
		else PsychErrorExitMsg(PsychError_user, "Invalid CodecFOURCC= parameter provided in movieoptions parameter. Must be exactly 4 characters! Parse error!");
	}

	// Assign numeric encoding quality level:
	// This is optional. We default to "normal quality":
	if ((poption = strstr(movieoptions, "EncodingQuality="))) {
		if ((sscanf(poption, "EncodingQuality=%f", &dummyFloat) == 1) && (dummyFloat >= 0) && (dummyFloat <= 1)) {
			// Map floating point quality level between 0.0 and 1.0 to 10 discrete levels:
			if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Encoding quality level %f selected for encoding of movie %i [%s].\n", dummyFloat, moviehandle, moviefile);

			// Rewrite "EncodingQuality=" string into "VideoQuality=" string, with proper
			// padding:      "EncodingQuality="
			// This way EncodingQuality in Quicktime lingo corresponds to
			// VideoQuality in GStreamer lingo:
			strncpy(poption, "   Videoquality=", strlen("   Videoquality="));
		}
		else PsychErrorExitMsg(PsychError_user, "Invalid EncodingQuality= parameter provided in movieoptions parameter. Parse error or out of valid 0 - 1 range!");
	}

	// Check for valid parameters. Also warn if some parameters are borderline for certain codecs:
	if ((framerate < 1) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In CreateMovie: Negative or zero 'framerate' %f units for moviehandle %i provided! Sounds like trouble ahead.\n", (float) framerate, moviehandle);
	if (width < 1) PsychErrorExitMsg(PsychError_user, "In CreateMovie: Invalid zero or negative 'width' for video frame size provided!");
	if ((width < 4) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In CreateMovie: 'width' of %i pixels for moviehandle %i provided! Some video codecs may malfunction with such a small width.\n", width, moviehandle);
	if ((width % 4 != 0) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In CreateMovie: 'width' of %i pixels for moviehandle %i provided! Some video codecs may malfunction with a width which is not a multiple of 4 or 16.\n", width, moviehandle);
	if (height < 1) PsychErrorExitMsg(PsychError_user, "In CreateMovie: Invalid zero or negative 'height' for video frame size provided!");
	if ((height < 4) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In CreateMovie: 'height' of %i pixels for moviehandle %i provided! Some video codecs may malfunction with such a small height.\n", height, moviehandle);

	// Full GStreamer launch line a la gst-launch command provided?
	if (strstr(movieoptions, "gst-launch")) {
		// Yes: We use movieoptions directly as launch line:
		movieoptions = strstr(movieoptions, "gst-launch");
        
		// Move string pointer behind the "gst-launch" word (plus a blank):
		movieoptions+= strlen("gst-launch ");
        
		// Can directly use this:
		sprintf(launchString, "%s", movieoptions);

		// With audio track?
		if (strstr(movieoptions, "name=ptbaudioappsrc")) doAudio = TRUE;
	}
	else {
		// No: Do our own parsing and setup:

		// Find the gst-launch style string for codecs and muxers:
		if (!PsychGetCodecLaunchLineFromString(movieoptions, &(codecString[0]))) {
			// No config for this format possible:
			if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating movie file with handle %i [%s] failed: Could not find matching codec setup.\n", moviehandle, moviefile);
			goto bail;
		}
        
		// With audio track?
		if (strstr(movieoptions, "AddAudioTrack")) doAudio = TRUE;
        
		// Build final launch string:
		if (doAudio) {
			// Video and audio:
			sprintf(launchString, "appsrc name=ptbvideoappsrc do-timestamp=0 stream-type=0 max-bytes=0 block=1 is-live=0 emit-signals=0 ! capsfilter caps=\"video/x-raw-rgb, bpp=(int)32, depth=(int)32, endianess=(int)4321, red_mask=(int)16711680, green_mask=(int)65280, blue_mask=(int)255, width=(int)%i, height=(int)%i, framerate=%i/1 \" ! videorate ! ffmpegcolorspace ! %s ! filesink name=ptbfilesink async=0 location=%s ", width, height, ((int) (framerate + 0.5)), codecString, moviefile);
		} else {
			// Video only:
			sprintf(launchString, "appsrc name=ptbvideoappsrc do-timestamp=0 stream-type=0 max-bytes=0 block=1 is-live=0 emit-signals=0 ! capsfilter caps=\"video/x-raw-rgb, bpp=(int)32, depth=(int)32, endianess=(int)4321, red_mask=(int)16711680, green_mask=(int)65280, blue_mask=(int)255, width=(int)%i, height=(int)%i, framerate=%i/1 \" ! videorate ! ffmpegcolorspace ! %s ! filesink name=ptbfilesink async=0 location=%s ", width, height, ((int) (framerate + 0.5)), codecString, moviefile);
		}
	}
        
	// Create a movie file for the destination movie:
	if (PsychPrefStateGet_Verbosity() > 3) {
		printf("PTB-INFO: Movie writing pipeline gst-launch line (without the -e option required on the command line!) is:\n");
		printf("gst-launch %s\n", launchString);
	}

	// Build pipeline from launch string:
	pwriterRec->Movie = gst_parse_launch((const gchar*) launchString, &myErr);
	if ((NULL == pwriterRec->Movie) || myErr) {
		if (PsychPrefStateGet_Verbosity() > 0) {
			printf("PTB-ERROR: In CreateMovie: Creating movie file with handle %i [%s] failed: Could not build pipeline.\n", moviehandle, moviefile);
			printf("PTB-ERROR: Parameters were: %s\n", movieoptions);
			printf("PTB-ERROR: Launch string was: %s\n", launchString);
			printf("PTB-ERROR: GStreamer error message was: %s\n", (char*) myErr->message);

		      // Special tips for the challenged:
		      if (strstr(myErr->message, "property")) {
			      // Bailed due to unsupported x264enc parameter "speed-preset" or "profile". Can be solved by upgrading
			      // GStreamer or the OS or the VideoCodec= override:
			      printf("PTB-TIP: The reason this failed is because your GStreamer codec installation is too outdated.\n");
			      printf("PTB-TIP: Either upgrade your GStreamer (plugin) installation to a more recent version,\n");
			      printf("PTB-TIP: or upgrade your operating system (e.g., Ubuntu 10.10 'Maverick Meercat' and later are fine).\n");
			      printf("PTB-TIP: A recent GStreamer installation is required to use all features and get optimal performance.\n");
			      printf("PTB-TIP: As a workaround, you can manually specify all codec settings, leaving out the unsupported\n");
			      printf("PTB-TIP: option. See 'help VideoRecording' on how to do that.\n\n");
		      }
		}

		goto bail;
	}

	// Get handle to ptbvideoappsrc:
	pwriterRec->ptbvideoappsrc = gst_bin_get_by_name(GST_BIN(pwriterRec->Movie), (const gchar *) "ptbvideoappsrc");
	if (NULL == pwriterRec->ptbvideoappsrc) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: In CreateMovie: Creating movie file with handle %i [%s] failed: Could not find ptbvideoappsrc pipeline element.\n", moviehandle, moviefile);
		goto bail;
	}

	// Get handle to ptbaudioappsrc:
	pwriterRec->ptbaudioappsrc = gst_bin_get_by_name(GST_BIN(pwriterRec->Movie), (const gchar *) "ptbaudioappsrc");
	if (doAudio && (NULL == pwriterRec->ptbaudioappsrc)) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: In CreateMovie: Creating movie file with handle %i [%s] failed: Could not find ptbaudioappsrc pipeline element.\n", moviehandle, moviefile);
		goto bail;
	}

	pwriterRec->Context = g_main_loop_new (NULL, FALSE);
	pwriterRec->bus = gst_pipeline_get_bus (GST_PIPELINE(pwriterRec->Movie));
	gst_bus_add_watch(pwriterRec->bus, (GstBusFunc) PsychMovieBusCallback, pwriterRec);
	gst_object_unref(pwriterRec->bus);

	// Start the pipeline:
	if (!PsychMoviePipelineSetState(pwriterRec->Movie, GST_STATE_PLAYING, 10)) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: In CreateMovie: Creating movie file with handle %i [%s] failed:  Failed to start movie encoding pipeline!\n", moviehandle, moviefile);
		goto bail;
	}

	PsychGSProcessMovieContext(pwriterRec->Context, FALSE);

	// Increment count of open movie writers:
	moviewritercount++;
	
	if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Moviehandle %i successfully opened for movie writing into file '%s'.\n", moviehandle, moviefile);

    // Should we dump the whole encoding pipeline graph to a file for visualization
    // with GraphViz? This can be controlled via PsychTweak('GStreamerDumpFilterGraph' dirname);
    if (getenv("GST_DEBUG_DUMP_DOT_DIR")) {
        // Dump complete encoding filter graph to a .dot file for later visualization with GraphViz:
        printf("PTB-DEBUG: Dumping movie encoder graph for movie %s to directory %s.\n", moviefile, getenv("GST_DEBUG_DUMP_DOT_DIR"));
        GST_DEBUG_BIN_TO_DOT_FILE_WITH_TS(GST_BIN(pwriterRec->Movie), GST_DEBUG_GRAPH_SHOW_ALL, "PsychMovieWritingGraph");
    }

	// Return new handle:
	return(moviehandle);

bail:
	if (pwriterRec->ptbvideoappsrc) gst_object_unref(GST_OBJECT(pwriterRec->ptbvideoappsrc));
	pwriterRec->ptbvideoappsrc = NULL;

	if (pwriterRec->ptbaudioappsrc) gst_object_unref(GST_OBJECT(pwriterRec->ptbaudioappsrc));
	pwriterRec->ptbaudioappsrc = NULL;

	if (pwriterRec->Movie) gst_object_unref(GST_OBJECT(pwriterRec->Movie));
	pwriterRec->Movie = NULL;

	if (pwriterRec->Context) g_main_loop_unref(pwriterRec->Context);
	pwriterRec->Context = NULL;

	// Return failure:
	return(-1);
}

int PsychFinalizeNewMovieFile(int movieHandle)
{
	int myErr = 0;
	GError *ret = NULL;

	PsychMovieWriterRecordType* pwriterRec = PsychGetMovieWriter(movieHandle, FALSE);

	if (NULL == pwriterRec->ptbvideoappsrc) return(0);

	// Release any pending buffers:
	if (pwriterRec->PixMap) gst_buffer_unref(pwriterRec->PixMap);
	pwriterRec->PixMap = NULL;

	PsychGSProcessMovieContext(pwriterRec->Context, FALSE);

	// Send EOS signal downstream:
	g_signal_emit_by_name(pwriterRec->ptbvideoappsrc, "end-of-stream", &ret);
	if (pwriterRec->ptbaudioappsrc) g_signal_emit_by_name(pwriterRec->ptbaudioappsrc, "end-of-stream", &ret);

	// Wait for eos flag to turn TRUE due to bus callback receiving the
	// downstream EOS event that we just sent out:
	while (!pwriterRec->eos) {
		PsychGSProcessMovieContext(pwriterRec->Context, FALSE);
		PsychYieldIntervalSeconds(0.010);
	}

	// Yield another 10 msecs after EOS signalled, just to be safe:
	PsychYieldIntervalSeconds(0.010);

	// Pause the encoding pipeline:
	if (!PsychMoviePipelineSetState(pwriterRec->Movie, GST_STATE_PAUSED, 10)) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed to pause movie encoding pipeline at close time!!\n");
	}

	// Stop the encoding pipeline:
	if (!PsychMoviePipelineSetState(pwriterRec->Movie, GST_STATE_READY, 10)) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed to stop movie encoding pipeline at close time!!\n");
	}

	// Shutdown and release encoding pipeline:
	if (!PsychMoviePipelineSetState(pwriterRec->Movie, GST_STATE_NULL, 10)) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed to shutdown movie encoding pipeline at close time!!\n");
	}

	gst_object_unref(GST_OBJECT(pwriterRec->Movie));
	pwriterRec->Movie = NULL;

	if (pwriterRec->ptbvideoappsrc) gst_object_unref(GST_OBJECT(pwriterRec->ptbvideoappsrc));
	pwriterRec->ptbvideoappsrc = NULL;

	if (pwriterRec->ptbaudioappsrc) gst_object_unref(GST_OBJECT(pwriterRec->ptbaudioappsrc));
	pwriterRec->ptbaudioappsrc = NULL;

	// Delete video context:
	if (pwriterRec->Context) g_main_loop_unref(pwriterRec->Context);
	pwriterRec->Context = NULL;

	// gst_object_unref(GST_OBJECT(pwriterRec->filesink));
	// pwriterRec->filesink = NULL;

	// Decrement count of active writers:
	moviewritercount--;

	// Return success/fail status:
	if ((myErr == 0) && (PsychPrefStateGet_Verbosity() > 3)) printf("PTB-INFO: Moviehandle %i successfully closed and movie written to filesystem.\n", movieHandle);

	return(myErr == 0);
}

// End of GStreamer routines.
#endif
