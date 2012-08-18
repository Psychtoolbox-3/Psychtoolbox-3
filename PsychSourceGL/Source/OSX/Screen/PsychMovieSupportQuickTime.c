/*
	PsychMovieSupportQuickTime.c
	
	PLATFORMS:	
	
		This is the OS X Core Graphics version.  
				
	AUTHORS:
	
		Mario Kleiner           mk              mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
        DESCRIPTION:
	
		Psychtoolbox functions for dealing with movies.
	
	NOTES:

*/

#include "PsychMovieSupportQuickTime.h"

#ifdef PSYCHQTAVAIL

// MK: Experimental switch: If set to 1, then we use old-style
// GWorlds for Quicktime rendering and manually convert them
// into OpenGL textures. This should default to zero for releases
// This is a test to find out if we can enable QT support on M$-Windows
// despite the braindamage of Apples Windows QT-7 SDK...
#if PSYCH_SYSTEM == PSYCH_WINDOWS
// We typedef all missing functions on Windows away...
typedef void* QTVisualContextRef;
Boolean QTVisualContextIsNewImageAvailable(void* a, void* b) { return(false); }
OSErr QTVisualContextCopyImageForTime(void* a, void* b, void* c) { return(noErr) };
void QTVisualContextRelease(void* a) { return; }
void QTVisualContextTask(void *a) { return; }
void CVOpenGLTextureRelease(void* a) { return; }
GLuint CVOpenGLTextureGetName(void* a) { return(0); }
void CVOpenGLTextureGetCleanTexCoords(void* a, float* b, float* c, float* d, float* e) { return; }

// ... and enable old style GWorld rendering, so this doesn't do any harm...
#define PSYCH_USE_QT_GWORLDS 1

#else

// OS-X: We don't use GWorlds but the new and shiny QTVisualcontexts and such...
#define PSYCH_USE_QT_GWORLDS 0
#include <pthread.h>
#include <sched.h>

#endif 

#define PSYCH_MAX_MOVIES 100

typedef struct {
    Movie   theMovie;
    QTVisualContextRef  QTMovieContext;
    QTAudioContextRef   QTAudioContext;
    GWorldPtr           QTMovieGWorld;
    PsychWindowRecordType* parentRecord;
    int loopflag;
    double movieduration;
    int nrframes;
    double fps;
    int width;
    int height;
    double last_pts;
    int nr_droppedframes;
    int specialFlags1;
} PsychMovieRecordType;

static PsychMovieRecordType movieRecordBANK[PSYCH_MAX_MOVIES];
static int numMovieRecords = 0;
static psych_bool firsttime = TRUE;

/*
 *     PsychQTMovieInit() -- Initialize movie subsystem.
 *     This routine is called by Screen's RegisterProject.c PsychModuleInit()
 *     routine at Screen load-time. It clears out the movieRecordBANK to
 *     bring the subsystem into a clean initial state.
 */
void PsychQTMovieInit(void)
{
    // Initialize movieRecordBANK with NULL-entries:
    int i;
    for (i=0; i < PSYCH_MAX_MOVIES; i++) {
        movieRecordBANK[i].theMovie = NULL;
        movieRecordBANK[i].QTMovieContext = NULL;
        movieRecordBANK[i].QTAudioContext = NULL;
        movieRecordBANK[i].QTMovieGWorld = NULL;
        movieRecordBANK[i].loopflag = 0;
    }    
    numMovieRecords = 0;
    
    return;
}

/** Internal helper function: Returns fps rate of movie and optionally
 *  the total number of video frames in the movie. Framecount is determined
 *  by stepping through the whole movie and counting frames. This can take
 *  significant time on big movie files.
 *
 *  Always returns fps as a double. Only counts and returns full framecount,
 *  if *nrframes is non-NULL.
 */
double PsychDetermineMovieFramecountAndFps(Movie theMovie, int* nrframes)
{
    // Count total number of videoframes: This code is derived from Apple
    // example code.
    long		myCount = -1;
    short		myFlags;
    TimeValue           myTime = 0;
    TimeValue           myDuration = 0;
    OSType		myTypes[1];
    // We want video samples.
    myTypes[0] = VisualMediaCharacteristic;
    // We want to begin with the first frame in the movie:
    myFlags = nextTimeStep + nextTimeEdgeOK;
    
    // We count either the first 3 frames if nrframes==NULL aka only
    // fps requested, or if framecount is requested, we count all frames.
    while (myTime >= 0 && (myCount<2 || nrframes!=NULL)) {
        myCount++;        
        // look for the next frame in the track; when there are no more frames,
        // myTime is set to -1, so we'll exit the while loop
        GetMovieNextInterestingTime(theMovie, myFlags, 1, myTypes, myTime, FloatToFixed(1), &myTime, &myDuration);        
        // after the first interesting time, don't include the time we're currently at
        myFlags = nextTimeStep;
    }    
    
    // Return optional count of frames:
    if (nrframes) *nrframes = (int) myCount;
    
    GoToBeginningOfMovie(theMovie);
    MoviesTask(theMovie, 0);
    
    // Compute and return frame rate in fps as (Ticks per second / Duration of single frame in ticks): 
    return((double) GetMovieTimeScale(theMovie) / (double) myDuration);    
}

int PsychQTGetMovieCount(void) {
	return(numMovieRecords);
}

/*
 *      PsychQTCreateMovie() -- Create a movie object.
 *
 *      This function tries to open a Quicktime-Moviefile and create an
 *      associated movie object for it.
 *
 *      win = Pointer to window record of associated onscreen window.
 *      moviename = char* with the name of the moviefile.
 *      preloadSecs = How many seconds of the movie should be preloaded/prefetched into RAM at movie open time?
 *      moviehandle = handle to the new movie.
 */
void PsychQTCreateMovie(PsychWindowRecordType *win, const char* moviename, double preloadSecs, int* moviehandle, int specialFlags1)
{
    Movie theMovie = NULL;
    QTVisualContextRef QTMovieContext = NULL;
    QTAudioContextRef  QTAudioContext = NULL;
    int i, slotid;
    OSErr error;
    CFStringRef movieLocation;
	CFURLRef movieURLLocation;
    CFStringRef coreAudioDeviceUID;
    psych_bool trueValue = TRUE;
    QTNewMoviePropertyElement newMovieProperties[4] = {0};
    int propcount = 0;
    char msgerr[10000];
    char errdesc[1000];
    Rect movierect;
    psych_bool printErrors;

    // Suppress output of error-messages if moviehandle == 1000. That means we
    // run in our own Posix-Thread, not in the Matlab-Thread. Printing via Matlabs
    // printing facilities would likely cause a terrible crash.
    printErrors = (*moviehandle == -1000) ? FALSE : TRUE;
    
    // Set movie handle to "failed" initially:
    *moviehandle = -1;

    // We startup the Quicktime subsystem only on first invocation.
    if (firsttime) {
#if PSYCH_SYSTEM == PSYCH_WINDOWS
        // Initialize Quicktime for Windows compatibility layer: This will fail if
        // QT isn't installed on the Windows machine...
        error = InitializeQTML(0);
        if (error!=noErr) {
            if (printErrors) {
                PsychErrorExitMsg(PsychError_internal, "Quicktime Media Layer initialization failed: Quicktime not properly installed?!?");
            } else return;

        }
#endif

        // Initialize Quicktime-Subsystem:
        error = EnterMovies();
        if (error!=noErr) {
            if (printErrors) PsychErrorExitMsg(PsychError_internal, "Quicktime EnterMovies() failed!!!"); else return;
        }
        firsttime = FALSE;
    }
    
    if (!PsychIsOnscreenWindow(win)) {
        if (printErrors) PsychErrorExitMsg(PsychError_user, "Provided windowPtr is not an onscreen window."); else return;
    }

    if (NULL==moviename) {
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
    
    // Create name-string for moviename:
    movieLocation = CFStringCreateWithCString (kCFAllocatorDefault, moviename, kCFStringEncodingASCII);

    // Zero-out new record in moviebank:
    movieRecordBANK[slotid].theMovie=NULL;    
    movieRecordBANK[slotid].QTMovieContext=NULL;    
    movieRecordBANK[slotid].QTAudioContext=NULL;
    movieRecordBANK[slotid].QTMovieGWorld=NULL;
    movieRecordBANK[slotid].parentRecord=NULL;

    // Assign flags:
    movieRecordBANK[slotid].specialFlags1 = specialFlags1;
    
    if (!PSYCH_USE_QT_GWORLDS) {
        // Create QTGLTextureContext:
#if PSYCH_SYSTEM != PSYCH_WINDOWS
        error = QTOpenGLTextureContextCreate (kCFAllocatorDefault,
                                                    win->targetSpecific.contextObject,
                                                    win->targetSpecific.pixelFormatObject,
                                                    NULL,
                                                    &QTMovieContext);
#endif
        if (error!=noErr) {
            if (printErrors) PsychErrorExitMsg(PsychError_internal, "OpenGL Quicktime visual context creation failed!!!"); else return;
        }        
    }

    // The Movie location 
    newMovieProperties[propcount].propClass = kQTPropertyClass_DataLocation;
	if (strstr(moviename, "http:") || strstr(moviename, "ftp:")) {
		// Open movie from URL, e.g., http- or ftp- server:
		movieURLLocation = CFURLCreateWithString(kCFAllocatorDefault, movieLocation, NULL);
		newMovieProperties[propcount].propID = kQTDataLocationPropertyID_CFURL;
		newMovieProperties[propcount].propValueSize = sizeof(movieURLLocation);
		newMovieProperties[propcount++].propValueAddress = (void*) &movieURLLocation;
	}
	else {
		// Open movie file from filesystem:
		newMovieProperties[propcount].propID = kQTDataLocationPropertyID_CFStringPosixPath;
		newMovieProperties[propcount].propValueSize = sizeof(CFStringRef);
		newMovieProperties[propcount++].propValueAddress = &movieLocation;
    }
	
    if (!PSYCH_USE_QT_GWORLDS) {
        // The Movie visual context
        newMovieProperties[propcount].propClass = kQTPropertyClass_Context;
        newMovieProperties[propcount].propID = kQTContextPropertyID_VisualContext;
        newMovieProperties[propcount].propValueSize = sizeof(QTVisualContextRef);
        newMovieProperties[propcount++].propValueAddress = &QTMovieContext;
    }
    
	// Also create and attach an audio decoding and playback context, unless forbidden
	// by specialFlags1 setting 0x2:
	// MK: Ok, this does not work, as the brain-damaged Apple QT implementation simply assigns
	// a default audio context if we skip this setup step :-( - It also isn't possible to
	// remove the context later via SetMovieAudioContext(NULL), and the wisdom of Apple's
	// Quicktime forum tells us it is really a no-go. Retarded tech from a retarded company...
    if (!(specialFlags1 & 0x2)) {
        // Create QTAudioContext for default CoreAudio device:
        coreAudioDeviceUID = NULL; // Use default audio-output device.
        error =QTAudioContextCreateForAudioDevice (kCFAllocatorDefault,
                                                   coreAudioDeviceUID,
                                                   NULL,
                                                   &QTAudioContext);
        if (error!=noErr) {
            if (printErrors) PsychErrorExitMsg(PsychError_internal, "Quicktime audio context creation failed!!!"); else return;
        }
        
        // The Movie audio context
        newMovieProperties[propcount].propClass = kQTPropertyClass_Context;
        newMovieProperties[propcount].propID = kQTContextPropertyID_AudioContext;
        newMovieProperties[propcount].propValueSize = sizeof(QTAudioContextRef);
        newMovieProperties[propcount++].propValueAddress = &QTAudioContext;
    }
    
    // The Movie active
    newMovieProperties[propcount].propClass = kQTPropertyClass_NewMovieProperty;
    newMovieProperties[propcount].propID = kQTNewMoviePropertyID_Active;
    newMovieProperties[propcount].propValueSize = sizeof(trueValue);
    newMovieProperties[propcount++].propValueAddress = &trueValue;
    
    // Instantiate the Movie
    error = NewMovieFromProperties(propcount, newMovieProperties, 0, NULL, &theMovie);
    if (error!=noErr) {
        QTVisualContextRelease(QTMovieContext);
        if (QTAudioContext) QTAudioContextRelease(QTAudioContext);
        switch(error) {
            case -2000:
            case -50:
            case -43:
                sprintf(errdesc, "File not found.");
            break;
            
            case -2048:
                sprintf(errdesc, "This is not a file that Quicktime understands.");
            break;
            
            case -2003:
                sprintf(errdesc, "Can't find media handler (codec) for this movie.");
            break;
            
            default:
                sprintf(errdesc, "Unknown: Check http://developer.apple.com/documentation/QuickTime/APIREF/ErrorCodes.htm#//apple_ref/doc/constant_group/Error_Codes");
        }
        
        sprintf(msgerr, "Couldn't load movie %s! Quicktime error code %i [%s]", moviename, (int) error, errdesc);
        *moviehandle = (int) error;
        if (printErrors) PsychErrorExitMsg(PsychError_user, msgerr); else return;
    }
    
    CFRelease(movieLocation);

    if (PSYCH_USE_QT_GWORLDS) {
        // Determine size of images in movie:
        GetMovieBox(theMovie, &movierect);
        
		// Only create a GWorld if movie frames contain at least 1 pixel. This way we skip GWorld
		// setup on "movies" which only consist of sound tracks.
		if ((movierect.right - movierect.left != 0) && (movierect.bottom - movierect.top != 0)) {
			// Create GWorld for this movie object:
			// error = QTNewGWorld(&movieRecordBANK[slotid].QTMovieGWorld, k32ABGRPixelFormat, &movierect,  NULL, NULL, 0);
			error = QTNewGWorld(&movieRecordBANK[slotid].QTMovieGWorld, 0, &movierect,  NULL, NULL, 0);
			if (error!=noErr) {
				if (QTAudioContext) QTAudioContextRelease(QTAudioContext);
				DisposeMovie(movieRecordBANK[slotid].theMovie);
				movieRecordBANK[slotid].theMovie=NULL;    
				if (printErrors) PsychErrorExitMsg(PsychError_internal, "Quicktime GWorld creation failed!!!"); else return;
			}
			
			// Attach this GWorld as rendering target for Quicktime:
			SetMovieGWorld(theMovie, movieRecordBANK[slotid].QTMovieGWorld, NULL);
		}
    }
    
    // Preload preloadSecs seconds of movie into system RAM for faster playback:
	if (preloadSecs > 0) LoadMovieIntoRam(theMovie, 0, ((long) preloadSecs + 0.5) * GetMovieTimeScale(theMovie),  keepInRam);
	// Special setting - 1 means: Load whole movie into RAM:
	if (preloadSecs == -1) LoadMovieIntoRam(theMovie, 0, GetMovieDuration(theMovie),  keepInRam);

    // We don't preroll: Didn't help for async playback, but leads to failure in
    // manual playback mode: PrerollMovie(theMovie, 0, FloatToFixed(1));

    // MoviesTask() it to make sure start of plaback will be as stutter-free as possible:
    MoviesTask(theMovie, 10000);
    
    // Assign new record in moviebank:
    movieRecordBANK[slotid].theMovie=theMovie;    
    movieRecordBANK[slotid].QTMovieContext=QTMovieContext;    
    movieRecordBANK[slotid].QTAudioContext=QTAudioContext;
    movieRecordBANK[slotid].parentRecord=win;
    movieRecordBANK[slotid].loopflag = 0;
    *moviehandle = slotid;

    // Increase counter:
    numMovieRecords++;

    // Compute basic movie properties - Duration and fps as well as image size:
    
    // Compute duration in seconds:
    movieRecordBANK[slotid].movieduration = (double) GetMovieDuration(theMovie) / (double) GetMovieTimeScale(theMovie);

    // Compute expected framerate, assuming a linear spacing between frames: It is derived as
    // reciprocal of the duration of the first video frame in the movie:
    movieRecordBANK[slotid].fps = PsychDetermineMovieFramecountAndFps(theMovie, NULL);

    // Determine size of images in movie:
    GetMovieBox(theMovie, &movierect);
    movieRecordBANK[slotid].width = movierect.right - movierect.left;
    movieRecordBANK[slotid].height = movierect.bottom - movierect.top;
    
    // We set nrframes == -1 to indicate that this value is not yet available.
    // Will do counting on first query for this parameter as it is very time-consuming:
    movieRecordBANK[slotid].nrframes = -1;
    
    return;
}

/*
 *  PsychQTGetMovieInfos() - Return basic information about a movie.
 *
 *  framecount = Total number of video frames in the movie, determined by counting.
 *  durationsecs = Total playback duration of the movie, in seconds.
 *  framerate = Estimated video playback framerate in frames per second (fps).
 *  width = Width of movie images in pixels.
 *  height = Height of movie images in pixels.
 *  nrdroppedframes = Total count of videoframes that had to be dropped during last movie playback,
 *                    in order to keep the movie synced with the realtime clock.
 */
void PsychQTGetMovieInfos(int moviehandle, int* width, int* height, int* framecount, double* durationsecs, double* framerate, int* nrdroppedframes)
{
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
    
    if (movieRecordBANK[moviehandle].theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }

    // Framecount requested?
    if (framecount) {
        // Try to fetch from internal data:
        *framecount = movieRecordBANK[moviehandle].nrframes;
        // Internal entry invalid?
        if (*framecount<0) {
            // Framecount unavailable: We need to trigger active counting:
            PsychDetermineMovieFramecountAndFps(movieRecordBANK[moviehandle].theMovie, &(movieRecordBANK[moviehandle].nrframes));
            *framecount = movieRecordBANK[moviehandle].nrframes;
        }
    }

    if (durationsecs) *durationsecs = movieRecordBANK[moviehandle].movieduration;
    if (framerate) *framerate = movieRecordBANK[moviehandle].fps;
    if (nrdroppedframes) *nrdroppedframes = movieRecordBANK[moviehandle].nr_droppedframes;
    if (width) *width = movieRecordBANK[moviehandle].width; 
    if (height) *height = movieRecordBANK[moviehandle].height; 
    return;
}

/*
 *  PsychQTDeleteMovie() -- Delete a movie object and release all associated ressources.
 */
void PsychQTDeleteMovie(int moviehandle)
{
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
    
    if (movieRecordBANK[moviehandle].theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }
    
    // Make sure OpenGL context of our parent window record is active:
    PsychSetGLContext(movieRecordBANK[moviehandle].parentRecord);
    
    // Stop movie playback immediately:
    MoviesTask(movieRecordBANK[moviehandle].theMovie, 0);
    StopMovie(movieRecordBANK[moviehandle].theMovie);
    glFinish();
    QTVisualContextTask(movieRecordBANK[moviehandle].QTMovieContext);
    MoviesTask(movieRecordBANK[moviehandle].theMovie, 0);
    QTVisualContextTask(movieRecordBANK[moviehandle].QTMovieContext);
    glFinish();
        
    // Delete movieobject for this handle:
    DisposeMovie(movieRecordBANK[moviehandle].theMovie);
    movieRecordBANK[moviehandle].theMovie=NULL;    

    // Delete GWorld if any:
    if (movieRecordBANK[moviehandle].QTMovieGWorld) DisposeGWorld(movieRecordBANK[moviehandle].QTMovieGWorld);
    movieRecordBANK[moviehandle].QTMovieGWorld = NULL;

    // Delete visual context for this movie:
    if (movieRecordBANK[moviehandle].QTMovieContext) QTVisualContextRelease(movieRecordBANK[moviehandle].QTMovieContext);
    movieRecordBANK[moviehandle].QTMovieContext = NULL;

    // Delete audio context for this movie:
    if (movieRecordBANK[moviehandle].QTAudioContext) QTAudioContextRelease(movieRecordBANK[moviehandle].QTAudioContext);
    movieRecordBANK[moviehandle].QTAudioContext = NULL;
    
    // Decrease counter:
    if (numMovieRecords>0) numMovieRecords--;
        
    return;
}

/*
 *  PsychQTDeleteAllMovies() -- Delete all movie objects and release all associated ressources.
 */
void PsychQTDeleteAllMovies(void)
{
    int i;
    for (i=0; i<PSYCH_MAX_MOVIES; i++) {
        if (movieRecordBANK[i].theMovie) PsychQTDeleteMovie(i);
    }
    return;
}


/*
 *  PsychQTGetTextureFromMovie() -- Create an OpenGL texture map from a specific videoframe from given movie object.
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
int PsychQTGetTextureFromMovie(PsychWindowRecordType *win, int moviehandle, int checkForImage, double timeindex, PsychWindowRecordType *out_texture, double *presentation_timestamp)
{
	static TimeValue myNextTimeCached = -2;
	static TimeValue nextFramesTimeCached = -2;
    TimeValue		myCurrTime;
    TimeValue		myNextTime;
    TimeValue       nextFramesTime=0;
    short		myFlags;
    OSType		myTypes[1];
    OSErr		error = noErr;
    Movie               theMovie;
    CVOpenGLTextureRef newImage = NULL;
    QTVisualContextRef  theMoviecontext;
    unsigned int failcount=0;
    float lowerLeft[2];
    float lowerRight[2];    
    float upperRight[2];    
    float upperLeft[2];
    GLuint texid;
    Rect rect;
    float rate;
    double targetdelta, realdelta, frames;
	PsychRectType outRect;

    if (!PsychIsOnscreenWindow(win)) {
        PsychErrorExitMsg(PsychError_user, "Need onscreen window ptr!!!");
    }
    
    // Activate OpenGL context of target window:
    PsychSetGLContext(win);

    // Explicitely disable Apple's Client storage extensions. For now they are not really useful to us.
    glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, GL_FALSE);
    
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided.");
    }
    
    if ((timeindex!=-1) && (timeindex < 0 || timeindex >= 10000.0)) {
        PsychErrorExitMsg(PsychError_user, "Invalid timeindex provided.");
    }
    
    if (NULL == out_texture && !checkForImage) {
        PsychErrorExitMsg(PsychError_user, "Sorry, use of a 'specialFlags2' setting of 2 for skipping texture creation is not supported by the Quicktime engine.");
    }
    
    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;
    theMoviecontext = movieRecordBANK[moviehandle].QTMovieContext;

    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle.");
    }

    // Check if end of movie is reached. Rewind, if so...
    if (IsMovieDone(theMovie) && movieRecordBANK[moviehandle].loopflag > 0) {
        if (GetMovieRate(theMovie)>0) {
            GoToBeginningOfMovie(theMovie);
        } else {
            GoToEndOfMovie(theMovie);
        }
    }
    
    // Is movie actively playing (automatic async playback, possibly with synced sound)?
    // If so, then we ignore the 'timeindex' parameter, because the automatic playback
    // process determines which frames should be delivered to PTB when. This function will
    // simply wait or poll for arrival/presence of a new frame that hasn't been fetched
    // in previous calls.
    if (0 == GetMovieRate(theMovie)) {
        // Movie playback inactive. We are in "manual" mode: No automatic async playback,
        // no synced audio output. The user just wants to manually fetch movie frames into
        // textures for manual playback in a standard Matlab-loop.

		// First pass - checking for new image?
		if (checkForImage) {
			// Image for specific point in time requested?
			if (timeindex >= 0) {
				// Yes. We try to retrieve the next possible image for requested timeindex.
				myCurrTime = (TimeValue) ((timeindex * (double) GetMovieTimeScale(theMovie)) + 0.5f);
			}
			else {
				// No. We just retrieve the next frame, given the current movie time.
				myCurrTime = GetMovieTime(theMovie, NULL);
			}
            
			// Retrieve timeindex of the closest image sample after myCurrTime:
			myFlags = nextTimeStep + nextTimeEdgeOK;	// We want the next frame in the movie's media.
			myTypes[0] = VisualMediaCharacteristic;		// We want video samples.
			GetMovieNextInterestingTime(theMovie, myFlags, 1, myTypes, myCurrTime, FloatToFixed(1), &myNextTime, &nextFramesTime);
			error = GetMoviesError();
			if (error != noErr) {
				PsychErrorExitMsg(PsychError_internal, "Failed to fetch texture from movie for given timeindex!");
			}
			
			// Found useful event?
			if (myNextTime == -1) {
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-WARNING: Bogus timevalue in movie track for movie %i. Trying to keep going.\n", moviehandle);
				
				// No. Just push timestamp to current time plus a little bit in the hope
				// this will get us unstuck:
				myNextTime = myCurrTime + (TimeValue) 1;
				nextFramesTime = (TimeValue) 0;
			}
			
			if (myNextTime != myNextTimeCached) {
				// Set movies current time to myNextTime, so the next frame will be fetched from there:
				SetMovieTimeValue(theMovie, myNextTime);
				
				// nextFramesTime is the timeindex to which we need to advance for retrieval of next frame: (see code below)
				nextFramesTime=myNextTime + nextFramesTime;
				
				if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: Current timevalue in movie track for movie %i is %lf secs.\n", moviehandle, (double) myNextTime / (double) GetMovieTimeScale(theMovie));
				if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: Next timevalue in movie track for movie %i is %lf secs.\n", moviehandle, (double) nextFramesTime / (double) GetMovieTimeScale(theMovie));
				
				// Cache values for 2nd pass:
				myNextTimeCached = myNextTime;
				nextFramesTimeCached = nextFramesTime;
			}
			else {
				// Somehow got stuck? Do nothing...
				if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: Seem to be a bit stuck at timevalue [for movie %i] of %lf secs. Nudging a bit forward...\n", moviehandle, (double) myNextTime / (double) GetMovieTimeScale(theMovie));
				// Nudge the timeindex a bit forware in the hope that this helps:
				SetMovieTimeValue(theMovie, GetMovieTime(theMovie, NULL) + 1);
			}
		}
		else {
			// This is the 2nd pass: Image fetching. Use cached values from first pass:
			// Caching in a static works because we're always called immediately for 2nd
			// pass after successfull return from 1st pass, and we're not multi-threaded,
			// i.e., don't need to be reentrant or thread-safe here:
			myNextTime = myNextTimeCached;
			nextFramesTime = nextFramesTimeCached;
			myNextTimeCached = -2;
		}
	}
    else {
        // myNextTime unavailable if in autoplayback-mode:
        myNextTime=-1;
    }
    
    // Presentation timestamp requested?
    if (presentation_timestamp) {
        // Already available?
        if (myNextTime==-1) {
            // Retrieve the exact presentation timestamp of the retrieved frame (in movietime):
            myFlags = nextTimeStep + nextTimeEdgeOK;            // We want the next frame in the movie's media.
            myTypes[0] = VisualMediaCharacteristic;		// We want video samples.
                                                                // We search backward for the closest available image for the current time. Either we get the current time
                                                                // if we happen to fetch a frame exactly when it becomes ready, or we get a bit earlier timestamp, which is
                                                                // the optimal presentation timestamp for this frame:
            GetMovieNextInterestingTime(theMovie, myFlags, 1, myTypes, GetMovieTime(theMovie, NULL), FloatToFixed(-1), &myNextTime, NULL);
        }
        // Convert pts (in Quicktime ticks) to pts in seconds since start of movie and return it:
        *presentation_timestamp = (double) myNextTime / (double) GetMovieTimeScale(theMovie);
    }

    // Allow quicktime visual context task to do its internal bookkeeping and cleanup work:
    if (theMoviecontext) QTVisualContextTask(theMoviecontext);

    // Perform decompress-operation:
    if (checkForImage) MoviesTask(theMovie, 0);
    
    // Should we just check for new image? If so, just return availability status:
    if (checkForImage) {
        if (PSYCH_USE_QT_GWORLDS) {
            // We use GWorlds. In this case we either suceed immediately due to the
            // synchronous nature of GWorld rendering, or we fail completely at end
            // of non-looping movie:
            if (IsMovieDone(theMovie) && movieRecordBANK[moviehandle].loopflag == 0) {
                // No new frame available and there won't be any in the future, because this is a non-looping
                // movie that has reached its end.
                return(-1);
            }
            
            // Is this the special case of a movie without video, but only sound? In that case,
			// we always return a 'false' because there ain't no image to return.
			if (movieRecordBANK[moviehandle].QTMovieGWorld == NULL) return(false);
			
			// Success!
            return(true);
        }
        
        // Code which uses QTVisualContextTasks...
        if (QTVisualContextIsNewImageAvailable(theMoviecontext, NULL)) {
            // New frame ready!
            return(true);
        }
        else if (IsMovieDone(theMovie) && movieRecordBANK[moviehandle].loopflag == 0) {
            // No new frame available and there won't be any in the future, because this is a non-looping
            // movie that has reached its end.
            return(-1);
        }
        else {
            // No new frame available yet:
            return(false);
        }
    }
    
    if (!PSYCH_USE_QT_GWORLDS) {
        // Blocking wait-code for non-GWorld mode:
        // Try up to 1000 iterations for arrival of requested image data in wait-mode:
        failcount=0;
        while ((failcount < 1000) && !QTVisualContextIsNewImageAvailable(theMoviecontext, NULL)) {
            PsychWaitIntervalSeconds(0.005);
            MoviesTask(theMovie, 0);
            failcount++;
        }
        
        // No new frame available and there won't be any in the future, because this is a non-looping
        // movie that has reached its end.
        if ((failcount>=1000) && IsMovieDone(theMovie) && (movieRecordBANK[moviehandle].loopflag == 0)) {
            return(-1);
        }
        
        // Fetch new OpenGL texture with the new movie image frame:
        error = QTVisualContextCopyImageForTime(theMoviecontext, kCFAllocatorDefault, NULL, &newImage);
        if ((error!=noErr) || newImage == NULL) {
            // We fail gracefully, but not totally silently:
            if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: In GetMovieImage: OpenGL<->Quicktime texture fetch failed! Premature end of movie or file corrupt?!?\n");
            return(-1);
        }
    
        // Disable client storage, if it was enabled:
        glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, GL_FALSE);
        
        // Build a standard PTB texture record:    
        CVOpenGLTextureGetCleanTexCoords (newImage, lowerLeft, lowerRight, upperRight, upperLeft);
        texid = CVOpenGLTextureGetName(newImage);
        
        // Assign texture rectangle:
        PsychMakeRect(outRect, upperLeft[0], upperLeft[1], lowerRight[0], lowerRight[1]);    
        
        // Set texture orientation as if it were an inverted Offscreen window: Upside-down.
        out_texture->textureOrientation = (CVOpenGLTextureIsFlipped(newImage)) ? 3 : 4;

        // Assign OpenGL texture id:
        out_texture->textureNumber = texid;
        
        // Store special texture object as part of the PTB texture record:
        out_texture->targetSpecific.QuickTimeGLTexture = newImage;
    }
    else {
        // Synchronous texture fetch code for GWorld rendering mode:
        // At this point, the GWorld should contain the source image for creating a
        // standard OpenGL texture:
        
        // Disable client storage, if it was enabled:
        glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, GL_FALSE);
        
        // Build a standard PTB texture record:    

        // Assign texture rectangle:
        GetMovieBox(theMovie, &rect);

        // Hack: Need to extend rect by 4 pixels, because GWorlds are 4 pixels-aligned via
        // image row padding:
        rect.right = rect.right + 4;
        PsychMakeRect(out_texture->rect, rect.left, rect.top, rect.right, rect.bottom);    
        
        // Set NULL - special texture object as part of the PTB texture record:
        out_texture->targetSpecific.QuickTimeGLTexture = NULL;

        // Set texture orientation as if it were an inverted Offscreen window: Upside-down.
        out_texture->textureOrientation = 3;
        
        // Setup a pointer to our GWorld as texture data pointer:
        out_texture->textureMemorySizeBytes = 0;

		// Quicktime textures are aligned on 4 Byte boundaries:
		out_texture->textureByteAligned = 4;

        // Lock GWorld:
        if(!LockPixels(GetGWorldPixMap(movieRecordBANK[moviehandle].QTMovieGWorld))) {
            // Locking surface failed! We abort.
            PsychErrorExitMsg(PsychError_internal, "PsychQTGetTextureFromMovie(): Locking GWorld pixmap surface failed!!!");
        }
        
        // This will retrieve an OpenGL compatible pointer to the GWorlds pixel data and assign it to our texmemptr:
        out_texture->textureMemory = (GLuint*) GetPixBaseAddr(GetGWorldPixMap(movieRecordBANK[moviehandle].QTMovieGWorld));
            
        // Let PsychCreateTexture() do the rest of the job of creating, setting up and
        // filling an OpenGL texture with GWorlds content:
        PsychCreateTexture(out_texture);
        
        // Undo hack from above after texture creation: Now we need the real width of the
        // texture for proper texture coordinate assignments in drawing code et al.
        rect.right = rect.right - 4;
        PsychMakeRect(outRect, rect.left, rect.top, rect.right, rect.bottom);    

        // Unlock GWorld surface. We do a glFinish() before, for safety reasons...
        //glFinish();
        UnlockPixels(GetGWorldPixMap(movieRecordBANK[moviehandle].QTMovieGWorld));

        // Ready to use the texture... We're done.
    }
    
	// Normalize texture rectangle and assign it:
	PsychNormalizeRect(outRect, out_texture->rect);
    PsychCopyRect(out_texture->clientrect, out_texture->rect);
	
    rate = FixedToFloat(GetMovieRate(theMovie));
    
    // Detection of dropped frames: This is a heuristic. We'll see how well it works out...
    if (rate && presentation_timestamp) {
        // Try to check for dropped frames in playback mode:

        // Expected delta between successive presentation timestamps:
        targetdelta = 1.0f / (movieRecordBANK[moviehandle].fps * rate);

        // Compute real delta, given rate and playback direction:
        if (rate>0) {
            realdelta = *presentation_timestamp - movieRecordBANK[moviehandle].last_pts;
            if (realdelta<0) realdelta = 0;
        }
        else {
            realdelta = -1.0 * (*presentation_timestamp - movieRecordBANK[moviehandle].last_pts);
            if (realdelta<0) realdelta = 0;
        }
        
        frames = realdelta / targetdelta;
        // Dropped frames?
        if (frames > 1 && movieRecordBANK[moviehandle].last_pts>=0) {
            movieRecordBANK[moviehandle].nr_droppedframes += (int) (frames - 1 + 0.5);
        }

        movieRecordBANK[moviehandle].last_pts = *presentation_timestamp;
    }
    
    // Manually advance movie time, if in fetch mode:
    if (0 == GetMovieRate(theMovie)) {
        // We are in manual fetch mode: Need to manually advance movie time to next
        // media sample:
		if (nextFramesTime == myNextTime) {
			// Invalid value? Try to hack something that gets us unstuck:
			myNextTime = GetMovieTime(theMovie, NULL);
			nextFramesTime = myNextTime + (TimeValue) 1;
		}

        SetMovieTimeValue(theMovie, nextFramesTime);        
    }
    
    // Check if end of movie is reached. Rewind, if so...
    if (IsMovieDone(theMovie) && movieRecordBANK[moviehandle].loopflag > 0) {
        if (GetMovieRate(theMovie)>0) {
            GoToBeginningOfMovie(theMovie);
        } else {
            GoToEndOfMovie(theMovie);
        }
    }

    return(TRUE);
}

/*
 *  PsychQTFreeMovieTexture() - Release texture memory for a Quicktime texture.
 *
 *  This routine is called by PsychDeleteTexture() in PsychTextureSupport.c
 *  It performs the special cleanup necessary for Quicktime created textures.
 *
 */
void PsychQTFreeMovieTexture(PsychWindowRecordType *win)
{
    // Fetch special Quicktime texture handle...
    CVOpenGLTextureRef theTexture = win->targetSpecific.QuickTimeGLTexture;

    // ...is this a quicktime movietexture? If not, just skip this routine.
    if (win->windowType!=kPsychTexture || NULL == theTexture) return;
    
    // Quicktime movie texture: Release it via special CoreVideo release function:
    CVOpenGLTextureRelease(theTexture);
    win->targetSpecific.QuickTimeGLTexture = NULL;
    
    // 0-out the textureNumber so our standard cleanup routine (glDeleteTextures) gets
    // skipped - if we don't do this, Quicktime & CoreVideo will leak memory like hell!!!
    win->textureNumber = 0;

    return;
}

/*
 *  PsychQTPlaybackRate() - Start- and stop movieplayback, set playback parameters.
 *
 *  moviehandle = Movie to start-/stop.
 *  playbackrate = zero == Stop playback, non-zero == Play movie with spec. rate,
 *                 e.g., 1 = forward, 2 = double speed forward, -1 = backward, ...
 *  loop = 0 = Play once. 1 = Loop, aka rewind at end of movie and restart.
 *  soundvolume = 0 == Mute sound playback, between 0.0 and 1.0 == Set volume to 0 - 100 %.
 *  Returns Number of dropped frames to keep playback in sync.
 */
int PsychQTPlaybackRate(int moviehandle, double playbackrate, int loop, double soundvolume)
{
    int dropped = 0;
    Movie   theMovie;
    
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
        
    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;    
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }

    // Fake the "no-audio decoding" request by muting the sound output,
    // the best we can do on this retarded QT implementation:
    if (movieRecordBANK[moviehandle].specialFlags1 & 0x2) soundvolume = 0;
    
    if (playbackrate != 0) {
        // Start playback of movie:
        SetMovieAudioMute(theMovie, (soundvolume==0) ? TRUE : FALSE, 0);
        SetMovieVolume(theMovie, (short) (soundvolume * 255.0));
        movieRecordBANK[moviehandle].loopflag = loop;
        movieRecordBANK[moviehandle].last_pts = -1.0;
        movieRecordBANK[moviehandle].nr_droppedframes = 0;
        SetMoviePreferredRate(theMovie, FloatToFixed(playbackrate));
        StartMovie(theMovie);
        MoviesTask(theMovie, 10000);
    }
    else {
        // Stop playback of movie:
        StopMovie(theMovie);
        QTVisualContextTask(movieRecordBANK[moviehandle].QTMovieContext);

        // Output count of dropped frames:
        if ((dropped=movieRecordBANK[moviehandle].nr_droppedframes) > 0) {
            if (PsychPrefStateGet_Verbosity()>2) printf("PTB-INFO: Movie playback had to drop %i frames of movie %i to keep playback in sync.\n", movieRecordBANK[moviehandle].nr_droppedframes, moviehandle); 
        }
    }
    
    return(dropped);
}

/*
 *  void PsychQTExitMovies() - Shutdown handler.
 *
 *  This routine is called by Screen('CloseAll') and on clear Screen time to
 *  do final cleanup. It deletes all Quicktime textures and releases all Quicktime
 *  movie objects.
 *
 */
void PsychQTExitMovies(void)
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
    PsychQTDeleteAllMovies();
    
    // Shutdown Quicktime toolbox: We skip this, because according to Apple its not necessary,
    // and for some reason it reliably hangs Matlab, so one has to force-quit it :-(
    // Don't do this: ExitMovies();

#if PSYCH_SYSTEM == PSYCH_WINDOWS
    // Shutdown Quicktime core system:
    ExitMovies();
    
    // Shutdown Quicktime for Windows compatibility layer:
    TerminateQTML();
#endif
    
    return;
}

/*
 *  PsychQTGetMovieTimeIndex()  -- Return current playback time of movie.
 */
double PsychQTGetMovieTimeIndex(int moviehandle)
{
    Movie   theMovie;
    
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
    
    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;    
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }

    // Retrieve timeindex:
    return((double) GetMovieTime(theMovie, NULL) / (double) GetMovieTimeScale(theMovie));
}

/*
 *  PsychQTSetMovieTimeIndex()  -- Set current playback time of movie.
 */
double PsychQTSetMovieTimeIndex(int moviehandle, double timeindex, psych_bool indexIsFrames)
{
    Movie		theMovie;
    double		oldtime;
	long		targetIndex, myIndex;
    short		myFlags;
    TimeValue	myTime;
    OSType		myTypes[1];
    
    if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIES) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided!");
    }
    
    // Fetch references to objects we need:
    theMovie = movieRecordBANK[moviehandle].theMovie;    
    if (theMovie == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid moviehandle provided. No movie associated with this handle !!!");
    }
    
    // Retrieve current timeindex:
    oldtime = (double) GetMovieTime(theMovie, NULL) / (double) GetMovieTimeScale(theMovie);
    
	// Index based or target time based seeking?
	if (indexIsFrames) {
		// Index based seeking:
		
		// Seek to given targetIndex:
		targetIndex = (long) (timeindex + 0.5);

		// We want video samples.
		myTypes[0] = VisualMediaCharacteristic;
		
		// We want to begin with the first frame in the movie:
		myFlags = nextTimeStep + nextTimeEdgeOK;
		
		// Start with iteration at beginning:
		myTime = 0;
		myIndex = -1;
		
		// We iterate until end of movie (myTime < 0) or targetIndex reached:
		while ((myTime >= 0) && (myIndex < targetIndex)) {
			// Increment our index position:
			myIndex++;
			
			// Look for the next frame in the track; when there are no more frames,
			// myTime is set to -1, so we'll exit the while loop
			GetMovieNextInterestingTime(theMovie, myFlags, 1, myTypes, myTime, FloatToFixed(1), &myTime, NULL);

			// after the first interesting time, don't include the time we're currently at
			myFlags = nextTimeStep;
		}    
		
		// Valid time for existing target frame?
		if (myTime >= 0) {
			// Yes. Seek to it:
			SetMovieTimeValue(theMovie, myTime);
		}

		// Done with seek.
	}
	else {
		// Time based seeking:

		// Set new timeindex as time in seconds:
		SetMovieTimeValue(theMovie, (TimeValue) (((timeindex * (double) GetMovieTimeScale(theMovie))) + 0.5f));

		// Done with seek.
	}

    // Check if end of movie is reached. Rewind, if so...
    if (IsMovieDone(theMovie) && movieRecordBANK[moviehandle].loopflag > 0) {
        if (GetMovieRate(theMovie) > 0) {
            GoToBeginningOfMovie(theMovie);
        } else {
            GoToEndOfMovie(theMovie);
        }
    }

	// Yield some processing time to Quicktime to update properly:
    MoviesTask(theMovie, 0);
    
    // Return old time value of previous position:
    return(oldtime);
}

#endif
