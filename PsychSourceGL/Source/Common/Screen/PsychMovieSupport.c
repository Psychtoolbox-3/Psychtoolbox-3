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
		for movie playback.
 
		On all platforms it is currently GStreamer or nothing.

	NOTES:

*/

#include "Screen.h"

#ifdef PTB_USE_GSTREAMER
#include "PsychMovieSupportGStreamer.h"
#endif

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

	return;
}

int PsychGetMovieCount(void)
{
	#ifdef PTB_USE_GSTREAMER
	return(PsychGSGetMovieCount());
	#endif

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
 *      non-realtime priority, we can do the work of opening a movie
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
    PsychCreateMovie(&(movieinfo->windowRecord), movieinfo->moviename, movieinfo->preloadSecs, &mymoviehandle, movieinfo->asyncFlag, movieinfo->specialFlags1, movieinfo->pixelFormat, movieinfo->maxNumberThreads);
	
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
 *      This function tries to open a moviefile and create an
 *      associated movie object for it.
 *
 *      win = Pointer to window record of associated onscreen window.
 *      moviename = char* with the name of the moviefile.
 *      moviehandle = handle to the new movie.
 */
void PsychCreateMovie(PsychWindowRecordType *win, const char* moviename, double preloadSecs, int* moviehandle, int asyncFlag, int specialFlags1, int pixelFormat, int maxNumberThreads)
{
    #ifdef PTB_USE_GSTREAMER
    PsychGSCreateMovie(win, moviename, preloadSecs, moviehandle, asyncFlag, specialFlags1, pixelFormat, maxNumberThreads);
    return;
    #endif

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
    #ifdef PTB_USE_GSTREAMER
    PsychGSGetMovieInfos(moviehandle, width, height, framecount, durationsecs, framerate, nrdroppedframes, aspectRatio);
    return;
    #endif

    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
}

/*
 *  PsychDeleteMovie() -- Delete a movie object and release all associated ressources.
 */
void PsychDeleteMovie(int moviehandle)
{
    #ifdef PTB_USE_GSTREAMER
    PsychGSDeleteMovie(moviehandle);
    return;
    #endif

	PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
}

/*
 *  PsychDeleteAllMovies() -- Delete all movie objects and release all associated ressources.
 */
void PsychDeleteAllMovies(void)
{
    #ifdef PTB_USE_GSTREAMER
    PsychGSDeleteAllMovies();
    return;
    #endif

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
    #ifdef PTB_USE_GSTREAMER
    return(PsychGSGetTextureFromMovie(win, moviehandle, checkForImage, timeindex, out_texture, presentation_timestamp));
    #endif

    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
    return(-1);
}

/*
 *  PsychFreeMovieTexture() - Release texture memory for a movie texture.
 *
 *  This routine is called by PsychDeleteTexture() in PsychTextureSupport.c
 *  It performs the special cleanup necessary for movie cached textures.
 *
 */
void PsychFreeMovieTexture(PsychWindowRecordType *win)
{
    #ifdef PTB_USE_GSTREAMER
    PsychGSFreeMovieTexture(win);
    return;
    #endif
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
    #ifdef PTB_USE_GSTREAMER
    return(PsychGSPlaybackRate(moviehandle, playbackrate, loop, soundvolume));
    #endif

    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
    return(0);
}

/*
 *  void PsychExitMovies() - Shutdown handler.
 *
 *  This routine is called by Screen('CloseAll') and on clear Screen time to
 *  do final cleanup. It deletes all textures and releases all movie objects.
 *  Then it shuts down the movie subsystem.
 *
 */
void PsychExitMovies(void)
{
    #ifdef PTB_USE_GSTREAMER
    PsychGSExitMovies();
    #endif

    return;
}

/*
 *  PsychGetMovieTimeIndex()  -- Return current playback time of movie.
 */
double PsychGetMovieTimeIndex(int moviehandle)
{
    #ifdef PTB_USE_GSTREAMER
    return(PsychGSGetMovieTimeIndex(moviehandle));
    #endif

    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
    return(0.0);
}

/*
 *  PsychSetMovieTimeIndex()  -- Set current playback time of movie.
 */
double PsychSetMovieTimeIndex(int moviehandle, double timeindex, psych_bool indexIsFrames)
{
    #ifdef PTB_USE_GSTREAMER
    return(PsychGSSetMovieTimeIndex(moviehandle, timeindex, indexIsFrames));
    #endif

    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Movie playback support not supported on your configuration.");
    return(0.0);
}
