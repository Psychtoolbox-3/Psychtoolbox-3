/*
	PsychSourceGL/Source/Common/Screen/PsychMovieSupport.c
	
	PLATFORMS:	
	
		All.
				
	AUTHORS:
	
		Mario Kleiner           mk              mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
        DESCRIPTION:
	
		Psychtoolbox functions for dealing with movies.

		This is the master dispatcher. Function herein get called by
		other parts of Psychtoolbox, most notably the SCREENxxx functions
		for movie playback. It dispatches into either the QuickTime based
		or GStreamer based implementations of the movie playback functions,
		depending on OS/Availability/Chosen configuration.

		On OS/X and Windows with 32-bit Matlab/Octave, usually Quicktime is
		used as movie playback engine. If GStreamer is supported on that
		platforms and enabled via Screen('Preference','OverrideMultimediaEngine', 1);
		then GStreamer is used.

		On all other platforms it is either GStreamer or nothing, as Quicktime
		isn't supported. These platforms are Window and OS/X with 64-bit Matlab/Octave,
		and Linux. As of 30.11.2010, only Linux is supported, 64-bit OS/X and
		Windows are planned.

	NOTES:

*/

#include "Screen.h"

#ifdef PSYCHQTAVAIL
#include "PsychMovieSupportQuickTime.h"
#endif

#if defined(__LP64__) || defined(_M_IA64)
#define PSYCHOTHER64BIT 1
#else
#define PSYCHOTHER64BIT 0
#endif

#ifdef PTB_USE_GSTREAMER
#include "PsychMovieSupportGStreamer.h"
#define USE_GSTREAMER 1
#else
#define USE_GSTREAMER 0
#endif

static psych_bool firstTime = TRUE;
static psych_bool doUsegs   = FALSE;

static psych_bool usegs(void) {
	// First invocation since module load time?
	if (firstTime) {
		// Yes. Need to probe which engine to use:
		firstTime = FALSE;

		// Default to Quicktime instead of GStreamer,
		// override in detection code below if appropriate:
		doUsegs = FALSE;

		// We always use GStreamer if we are running on
		// Linux, or 64-Bit builds on OS/X or Windows, as
		// these systems only support GStreamer, but they
		// support it consistently:
		if ((PSYCH_SYSTEM == PSYCH_LINUX) || PSYCHOTHER64BIT) {
			// Yep: Unconditionally use GStreamer:
			doUsegs = TRUE;
		} else {
			// This is a 32-bit build on Windows or OS/X.
			// We use GStreamer if it is supported and usercode
			// wants to use it, according to preference setting:
			if (USE_GSTREAMER && (PsychPrefStateGet_UseGStreamer()==1)) {
				doUsegs = TRUE;
			}
		}
        
        // Signal use of GStreamer to userspace via preference setting:
        if (doUsegs) PsychPrefStateSet_UseGStreamer(1);
	}

	// Return cached engine use flag:
	return(doUsegs);
}

/*
 *     PsychMovieInit() -- Initialize movie subsystem.
 *     This routine is called by Screen's RegisterProject.c PsychModuleInit()
 *     routine at Screen load-time. It clears out the movieRecordBANK to
 *     bring the subsystem into a clean initial state.
 */
void PsychMovieInit(void)
{
	#ifdef PTB_USE_GSTREAMER
	// This will also first-time initialize Glib's threading system:
	PsychGSMovieInit();
	#endif

	#ifdef PSYCHQTAVAIL
	PsychQTMovieInit();
	#endif

	// Reset firstTime flag:
	firstTime = TRUE;
	return;
}

int PsychGetMovieCount(void)
{
	if (usegs()) {
	#ifdef PTB_USE_GSTREAMER
	return(PsychGSGetMovieCount());
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	return(PsychQTGetMovieCount());
	#endif
	}

	return(0);
}

/*
 *      PsychAsyncCreateMovie() -- Open a movie file in the background.
 *
 *      This function is called by SCREENOpenMovie as main-function of
 *      a new Posix-Thread for background movie loading. It simply calls
 *      PsychCreateMovie(), waiting for it to return a new moviehandle.
 *      Then it returns those info and terminates.
 *      -> By calling PsychCreateMovie from the run-function of a dedicated
 *      Posix-Thread which runs independent of the main Matlab/PTB Thread with
 *      non-realtime priority, we can do the work of opening a Quicktime movie
 *      in the background, hopefully not affecting the timing of the main PTB
 *      thread too much.
 */
void* PsychAsyncCreateMovie(void* inmovieinfo)
{
    int rc;
	PsychAsyncMovieInfo* movieinfo = (PsychAsyncMovieInfo*) inmovieinfo;
	
    // The special value -1000 tells PsychCreateMovie to not output any error-
    // messages as this could easily crash Matlab/Octave.
    int mymoviehandle=-1000;

    // Reduce our scheduling priority to the minimum value of zero, so
    // we do not interfere too much with the PTB main thread:
    if ((rc=PsychSetThreadPriority(NULL, 0, 0))!=0) {
        printf("PTB-WARNING: In PsychAsyncCreateMovie(): Failed to lower my priority to non-realtime [System errorcode %i]. Expect timing problems for movie playback!", rc);
    }
    
    // Execute our normal OpenMovie function: This does the hard work:
    PsychCreateMovie(&(movieinfo->windowRecord), movieinfo->moviename, movieinfo->preloadSecs, &mymoviehandle, movieinfo->asyncFlag, movieinfo->specialFlags1);
	
    // Ok, either we have a moviehandle to a valid movie, or we failed, which would
    // be signalled to the calling function via some negative moviehandle:
    movieinfo->moviehandle = mymoviehandle; // Return moviehandle.
    movieinfo->asyncstate = 2; // Set state to "Completed"
    
    // Exit from the routine. This will automatically terminate our Thread.
    return(NULL);
}

/*
 *      PsychCreateMovie() -- Create a movie object.
 *
 *      This function tries to open a Quicktime-Moviefile and create an
 *      associated movie object for it.
 *
 *      win = Pointer to window record of associated onscreen window.
 *      moviename = char* with the name of the moviefile.
 *      moviehandle = handle to the new movie.
 */
void PsychCreateMovie(PsychWindowRecordType *win, const char* moviename, double preloadSecs, int* moviehandle, int asyncFlag, int specialFlags1)
{
	if (usegs()) {
	#ifdef PTB_USE_GSTREAMER
	PsychGSCreateMovie(win, moviename, preloadSecs, moviehandle, asyncFlag, specialFlags1);
	return;
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	PsychQTCreateMovie(win, moviename, preloadSecs, moviehandle, specialFlags1);
	return;
	#endif
	}

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
}

/*
 *  PsychGetMovieInfo() - Return basic information about a movie.
 *
 *  framecount = Total number of video frames in the movie, determined by counting.
 *  durationsecs = Total playback duration of the movie, in seconds.
 *  framerate = Estimated video playback framerate in frames per second (fps).
 *  width = Width of movie images in pixels.
 *  height = Height of movie images in pixels.
 *  nrdroppedframes = Total count of videoframes that had to be dropped during last movie playback,
 *                    in order to keep the movie synced with the realtime clock.
 *  aspectRatio = Pixel aspect ratio of encoded pixels in movie frame. (GStreamer only).
 */
void PsychGetMovieInfos(int moviehandle, int* width, int* height, int* framecount, double* durationsecs, double* framerate, int* nrdroppedframes, double* aspectRatio)
{
	if (usegs()) {
        #ifdef PTB_USE_GSTREAMER
	PsychGSGetMovieInfos(moviehandle, width, height, framecount, durationsecs, framerate, nrdroppedframes, aspectRatio);
	return;
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	PsychQTGetMovieInfos(moviehandle, width, height, framecount, durationsecs, framerate, nrdroppedframes);
	if (aspectRatio) *aspectRatio = 1.0;
	return;
	#endif
	}

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
}

/*
 *  PsychDeleteMovie() -- Delete a movie object and release all associated ressources.
 */
void PsychDeleteMovie(int moviehandle)
{
	if (usegs()) {
        #ifdef PTB_USE_GSTREAMER
	PsychGSDeleteMovie(moviehandle);
	return;
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	PsychQTDeleteMovie(moviehandle);
	return;
	#endif
	}

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
}

/*
 *  PsychDeleteAllMovies() -- Delete all movie objects and release all associated ressources.
 */
void PsychDeleteAllMovies(void)
{
	if (usegs()) {
        #ifdef PTB_USE_GSTREAMER
	PsychGSDeleteAllMovies();
	return;
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	PsychQTDeleteAllMovies();
	return;
	#endif
	}

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
}

/*
 *  PsychGetTextureFromMovie() -- Create an OpenGL texture map from a specific videoframe from given movie object.
 *
 *  win = Window pointer of onscreen window for which a OpenGL texture should be created.
 *  moviehandle = Handle to the movie object.
 *  checkForImage = true == Just check if new image available, false == really retrieve the image, blocking if necessary.
 *  timeindex = When not in playback mode, this allows specification of a requested frame by presentation time.
 *              If set to -1, or if in realtime playback mode, this parameter is ignored and the next video frame is returned.
 *  out_texture = Pointer to the Psychtoolbox texture-record where the new texture should be stored.
 *  presentation_timestamp = A ptr to a double variable, where the presentation timestamp of the returned frame should be stored.
 *
 *  Returns true (1) on success, false (0) if no new image available, -1 if no new image available and there won't be any in future.
 */
int PsychGetTextureFromMovie(PsychWindowRecordType *win, int moviehandle, int checkForImage, double timeindex, PsychWindowRecordType *out_texture, double *presentation_timestamp)
{
	if (usegs()) {
	#ifdef PTB_USE_GSTREAMER
	return(PsychGSGetTextureFromMovie(win, moviehandle, checkForImage, timeindex, out_texture, presentation_timestamp));
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	return(PsychQTGetTextureFromMovie(win, moviehandle, checkForImage, timeindex, out_texture, presentation_timestamp));
	#endif
	}

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
	return(-1);
}

/*
 *  PsychFreeMovieTexture() - Release texture memory for a Quicktime texture.
 *
 *  This routine is called by PsychDeleteTexture() in PsychTextureSupport.c
 *  It performs the special cleanup necessary for Quicktime created textures.
 *
 */
void PsychFreeMovieTexture(PsychWindowRecordType *win)
{
	if (usegs()) {
	#ifdef PTB_USE_GSTREAMER
	PsychGSFreeMovieTexture(win);
	return;
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	PsychQTFreeMovieTexture(win);
	return;
	#endif
	}
}

/*
 *  PsychPlaybackRate() - Start- and stop movieplayback, set playback parameters.
 *
 *  moviehandle = Movie to start-/stop.
 *  playbackrate = zero == Stop playback, non-zero == Play movie with spec. rate,
 *                 e.g., 1 = forward, 2 = double speed forward, -1 = backward, ...
 *  loop = 0 = Play once. 1 = Loop, aka rewind at end of movie and restart.
 *  soundvolume = 0 == Mute sound playback, between 0.0 and 1.0 == Set volume to 0 - 100 %.
 *  Returns Number of dropped frames to keep playback in sync.
 */
int PsychPlaybackRate(int moviehandle, double playbackrate, int loop, double soundvolume)
{
	if (usegs()) {
        #ifdef PTB_USE_GSTREAMER
	return(PsychGSPlaybackRate(moviehandle, playbackrate, loop, soundvolume));
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	return(PsychQTPlaybackRate(moviehandle, playbackrate, loop, soundvolume));
	#endif
	}

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
	return(0);
}

/*
 *  void PsychExitMovies() - Shutdown handler.
 *
 *  This routine is called by Screen('CloseAll') and on clear Screen time to
 *  do final cleanup. It deletes all Quicktime textures and releases all Quicktime
 *  movie objects. Then it shuts down the Quicktime subsystem.
 *
 */
void PsychExitMovies(void)
{
	if (usegs()) {
        #ifdef PTB_USE_GSTREAMER
	PsychGSExitMovies();
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	PsychQTExitMovies();
	#endif
	}

	// Reset firstTime flag:
	firstTime = TRUE;
	return;
}

/*
 *  PsychGetMovieTimeIndex()  -- Return current playback time of movie.
 */
double PsychGetMovieTimeIndex(int moviehandle)
{

	if (usegs()) {
        #ifdef PTB_USE_GSTREAMER
	return(PsychGSGetMovieTimeIndex(moviehandle));
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	return(PsychQTGetMovieTimeIndex(moviehandle));
	#endif
	}

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
	return(0.0);
}

/*
 *  PsychSetMovieTimeIndex()  -- Set current playback time of movie.
 */
double PsychSetMovieTimeIndex(int moviehandle, double timeindex, psych_bool indexIsFrames)
{
	if (usegs()) {
        #ifdef PTB_USE_GSTREAMER
	return(PsychGSSetMovieTimeIndex(moviehandle, timeindex, indexIsFrames));
	#endif
	} else {
	#ifdef PSYCHQTAVAIL
	return(PsychQTSetMovieTimeIndex(moviehandle, timeindex, indexIsFrames));
	#endif
	}

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
	return(0.0);
}
