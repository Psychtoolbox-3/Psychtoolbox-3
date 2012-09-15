/*
  Psychtoolbox3/Source/Common/SCREENOpenMovie.c		
  
  AUTHORS:

    mario.kleiner at tuebingen.mpg.de   mk
  
  PLATFORMS:	

  This file should build on any platform. 

  HISTORY:

  10/23/05  mk		Created. 
 
  DESCRIPTION:
  
  Open a named movie file from the filesystem, create and initialize a corresponding movie object
  and return a handle to it to MATLAB space.
 
  On OS-X and Windows, all movie/multimedia handling functions are implemented via the Apple Quicktime API,
  version 7 or later. On Linux we use GStreamer: We need at least the gstreamer core library and base plugin
  set installed for it to work, but more plugins never hurt.
  
  TO DO:

*/

#include "Screen.h"

static char useString[] = "[ moviePtr [duration] [fps] [width] [height] [count] [aspectRatio]]=Screen('OpenMovie', windowPtr, moviefile [, async=0] [, preloadSecs=1] [, specialFlags1=0][, pixelFormat=4][, maxNumberThreads=-1]);";
static char synopsisString[] = 
		"Try to open the multimediafile 'moviefile' for playback in onscreen window 'windowPtr' and "
        "return a handle 'moviePtr' on success.\nOn OS-X and Windows, media files are handled by use of "
        "Apples Quicktime-7 API. On Linux, and optionally Windows, the GStreamer multimedia framework is used. "
        "The following movie properties are optionally returned: 'duration' Total duration of movie in seconds. "
        "'fps' Video playback framerate, assuming a linear spacing of videoframes in time. There may "
        "exist exotic movie formats which don't have this linear spacing. In that case, 'fps' would "
        "return bogus values and the check for skipped frames would report bogus values as well. "
        "'width' Width of the images contained in the movie. 'height' Height of the images.\n"
        "'count' Total number of videoframes in the movie. Determined by counting, so querying 'count' "
        "can significantly increase the execution time of this command.\n"
		"'aspectRatio' Pixel aspect ratio of pixels in the video frames. Typically 1.0 for square pixels. "
		"Aspect ratio is only detected with the GStreamer playback engine. Quicktime always returns 1.0.\n"
        "If you want to play multiple movies in succession with lowest possible delay inbetween the movies "
        "then you can ask PTB to load a movie in the background while another movie is still playing: "
        "Call this function with the 'async' flag set to 1. This will initiate the background load operation. "
        "After some sufficient time has passed, you can call the 'OpenMovie' function again, this time with "
        "the 'async' flag set to zero. Now the function will return a valid movie handle for playback. Background "
        "loading of movies currently does only work well with movies that don't have sound, unless you use Linux.\n"
        "If you use the GStreamer playback engine and all your movies have exactly the same format and only differ "
        "in duration and content, but not in image size, color depth or format, or fps, then you can also use "
        "an aync setting of 2 and provide the 'moviePtr' handle of an already opened movie in the 'preloadSecs' "
        "parameter. This will queue the movie 'moviefile' as a successor to the currently playing moviefile in "
        "'moviePtr'. Queuing movies this way is more efficient than async flag setting 1, although also more "
        "restricted.\n"
		"If the 'async' flag also contains the number 4 or is equal to 4, then movie playback will not automatically "
		"drop video frames to preserve audio-video sync in case fetching and display of video frames by your script "
		"is delayed or too slow. This has the disadvantage that you'll need to take care of audio-video sync and framerate "
		"control yourself by proper comparison of movie presentation timestamps and GetSecs or Screen('Flip') timestamps. "
		"The advantage is, that after start of playback the playback engine can internally predecode and buffer up to "
		"'preloadSecs' seconds worth of video and audio data. This may allow complex movies to play more smoothly or "
		"at higher framerates. This option is only supported with the GStreamer playback engine.\n"
		"'preloadSecs' This optional parameter allows to ask Screen() to load at least the first 'preloadSecs' "
		"seconds of the movie into system RAM. By default, with Quicktime playback, the first second of the "
		"movie file is loaded into RAM. This potentially allows for more stutter-free playback, but your mileage "
		"may vary, depending on movie format, storage medium and lots of other factors. In most cases, the default "
		"setting is perfectly sufficient. The special setting -1 means: Load whole movie into RAM. Caution: Long "
		"movies may cause your system to run low on memory and have disastrous effects on playback performance! "
		"Also, the exact type of buffering applied depends a lot on the movie playback engine and movie format, ie., "
		"what the parameter exactly does can be different for Quicktime vs. GStreamer, but it always affects the "
		"buffering behaviour and capacity of buffering in some meaningful way.\n"
		"'specialFlags1' Optional flags, numbers to be added together: 1 = Use YUV video decoding instead of RGBA, if "
		"supported by movie codec and GPU - May be more efficient. 2 = Don't decode and use sound - May be more efficient. "
		"On Linux you may need to specify a setting of 2 if you try to use movie playback at the same time as "
		"PsychPortAudio sound output, otherwise movie playback may hang. A flag of 4 will draw motion vectors on top "
        "of decoded video frames, for debugging or entertainment. A flag of 8 will ask the video decoder to skip all "
        "B-Frames during decoding to reduce processor load on very slow machines. Not all codecs may support flags 4 or "
        "8, in which case these flags are silently ignored. A flag of 16 asks Screen to convert all video textures "
        "immediately into a format which makes them useable as offscreen windows, and for the Screen('TransformTexture') "
        "function as well as for drawing them with your own custom GLSL shaders. Normally this conversion would be "
        "deferred until needed, ie. it would get skipped if you would just draw the texture regularly. If you know "
        "already that you want to use the texture with one of the given functions, manually triggering the conversion "
        "via this flag may be a bit more efficient - or convenient if you want to use your own GLSL shaders.\n"
        "'pixelFormat' optional argument specifying the pixel format of decoded video frames. Not all possible valid "
        "values are supported by all video codecs, graphics cards and operating systems. If an unsupported format is "
        "requested, Screen() will try to choose the closest matching format that meets or exceeds the specified format, "
        "at a performance or efficiency penalty. If no sufficiently close match is possible without severely degraded "
        "performance or other restrictions, the function will abort with an error. The following formats are supported "
        "on some setups: 1 = Luminance/Greyscale image, 2 = Luminance+Alpha, 3 = RGB 8 bit per channel, 4 = RGBA8, "
        "5 = YUV 4:2:2 packed pixel format on some graphics hardware, 6 = YUV-I420 planar format, using GLSL shaders "
        "for color space conversion on suitable graphics cards. 7 or 8 = Y8-Y800 planar format, using GLSL shaders. "
        "The always supported default is '4' == RGBA8 format. "
        "A setting of 6 (for color) or 7/8 (for grayscale) for selection of YUV-I420/Y8-Y800 format, as supported by at "
        "least the H264 and HuffYUV video codecs on any GPU with "
        "shader support, can be especially efficient for fast playback of high resolution video. As this format uses "
        "shaders for post-processing, it should be fast for texture drawing, but can incur significant overhead if you "
        "try to draw into a texture of this format, or try to post-process it via Screen('TransformTexture'). If you try "
        "to attach your own shaders to such a texture during Screen('DrawTexture'), you will need to implement color "
        "conversion yourself in your shaders, as your shaders would override Screen's builtin color conversion shader.\n"
        "'maxNumberThreads' Optional parameter which allows to set the maximum number of parallel processing threads "
        "that should be used by multi-threaded video codecs to decode the movie. The parameter has no effect on single "
        "threaded codecs and default behaviour is to let the codec do whatever it wants. A setting of zero tells the "
        "codec to use multi-threaded decoding with a number of threads that is auto-selected to be optimal for your given "
        "computer. A number n greater zero asks the codec to use at most n threads for decoding. The most safe choice is "
        "to not specify this parameter - this should work even with problematic movie formats. If you need higher playback "
        "performance, e.g., for high resolution video or high framerate playback, you should set the parameter to zero to "
        "leave the optimal choice to the video codec. This should work flawlessly with well encoded high quality movie files "
        "and can provide a significant performance boost on multi-core computers. Specify a discrete non-zero number of threads "
        "if you want to benefit from multi-core decoding but want to prevent movie playback from using up all available computation "
        "power, e.g., because you want to run some other timing-sensitive tasks in parallel and want to make sure to leave some processor "
        "cores dedicated to them.\n"
        "CAUTION: On 32-Bit OS/X, some movie files, e.g., MPEG-1 movies sometimes cause Matlab to hang. This seems to be "
        "a bad interaction between parts of Apples Quicktime toolkit and Matlabs Java Virtual Machine (JVM). "
        "If you experience stability problems, please start Matlab with JVM and desktop disabled, e.g., "
        "with the command: 'matlab -nojvm'. An example command sequence in a terminal window could be: "
        "/Applications/MATLAB701/bin/matlab -nojvm ";

static char seeAlsoString[] = "CloseMovie PlayMovie GetMovieImage GetMovieTimeIndex SetMovieTimeIndex";

PsychAsyncMovieInfo asyncmovieinfo;

PsychError SCREENOpenMovie(void) 
{
        PsychWindowRecordType                   *windowRecord;
        char                                    *moviefile;
        int                                     moviehandle = -1;
        int                                     framecount;
        double                                  durationsecs;
        double                                  framerate;
        double                                  aspectRatio;
        int                                     width;
        int                                     height;
        int                                     asyncFlag = 0;
		int                                     specialFlags1 = 0;
        static psych_bool                       firstTime = TRUE;
        double                                  preloadSecs = 1;
        int                                     rc;
        int                                     pixelFormat = 4;
        int                                     maxNumberThreads = -1;

        if (firstTime) {
            // Setup asyncopeninfo on first invocation:
            firstTime = FALSE;
            asyncmovieinfo.asyncstate = 0; // State = No async open in progress.
        }
        
        // All sub functions should have these two lines
        PsychPushHelp(useString, synopsisString, seeAlsoString);
        if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

        PsychErrorExit(PsychCapNumInputArgs(7));            // Max. 7 input args.
        PsychErrorExit(PsychRequireNumInputArgs(1));        // Min. 1 input args required.
        PsychErrorExit(PsychCapNumOutputArgs(7));           // Max. 7 output args.
        
        // Get the window record from the window record argument and get info from the window record
        windowRecord = NULL;
        PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, FALSE, &windowRecord);
        // Only onscreen windows allowed:
        if(windowRecord && !PsychIsOnscreenWindow(windowRecord)) {
            PsychErrorExitMsg(PsychError_user, "OpenMovie called on something else than an onscreen window.");
        }
        
        // Get the movie name string:
        moviefile = NULL;
        PsychAllocInCharArg(2, kPsychArgRequired, &moviefile);

        // Get the (optional) asyncFlag:
        PsychCopyInIntegerArg(3, FALSE, &asyncFlag);

        PsychCopyInDoubleArg(4, FALSE, &preloadSecs);
        if (preloadSecs < 0 && preloadSecs!= -1 && preloadSecs!= -2) PsychErrorExitMsg(PsychError_user, "OpenMovie called with invalid (negative, but not equal -1) 'preloadSecs' argument!");

        // Get the (optional) specialFlags1:
        PsychCopyInIntegerArg(5, FALSE, &specialFlags1);
		if (specialFlags1 < 0) PsychErrorExitMsg(PsychError_user, "OpenMovie called with invalid 'specialFlags1' setting! Only positive values allowed.");

        // Get the (optional) pixelFormat:
        PsychCopyInIntegerArg(6, FALSE, &pixelFormat);
        if (pixelFormat < 1 || pixelFormat > 8) PsychErrorExitMsg(PsychError_user, "OpenMovie called with invalid 'pixelFormat' setting! Only values 1 to 8 are allowed.");

        // Get the (optional) maxNumberThreads:
        PsychCopyInIntegerArg(7, FALSE, &maxNumberThreads);
        if (maxNumberThreads < -1) PsychErrorExitMsg(PsychError_user, "OpenMovie called with invalid 'maxNumberThreads' setting! Only values of -1 or greater are allowed.");
    
        // Queueing of a new movie for seamless playback requested?
        if (asyncFlag & 2) {
            // Yes. Do a special call, just passing the moviename of the next
            // movie to play. Pass the relevant moviehandle as retrieved from
            // preloadSecs:
            moviehandle = (int) preloadSecs;
            preloadSecs = 0;
            PsychCreateMovie(windowRecord, moviefile, preloadSecs, &moviehandle, asyncFlag, specialFlags1, pixelFormat, maxNumberThreads);
            if (moviehandle == -1) PsychErrorExitMsg(PsychError_user, "Could not queue new moviefile for gapless playback.");
            return(PsychError_none);
        }

        // Asynchronous Open operation in progress or requested?
        if ((asyncmovieinfo.asyncstate == 0) && !(asyncFlag & 1)) {
            // No. We should just synchronously open the movie:

            // Try to open the named 'moviefile' and create & initialize a corresponding movie object.
            // A MATLAB handle to the movie object is returned upon successfull operation.
            PsychCreateMovie(windowRecord, moviefile, preloadSecs, &moviehandle, asyncFlag, specialFlags1, pixelFormat, maxNumberThreads);
        }
        else {
            // Asynchronous open operation requested or running:
            switch(asyncmovieinfo.asyncstate) {
                case 0: // No async open running, but async open requested
                    // Fill all information needed for opening the movie into the info struct:
                    asyncmovieinfo.asyncstate = 1; // Mark state as "Operation in progress"
                    asyncmovieinfo.moviename = strdup(moviefile);
                    asyncmovieinfo.preloadSecs = preloadSecs;
                    asyncmovieinfo.asyncFlag = asyncFlag;
                    asyncmovieinfo.specialFlags1 = specialFlags1;
                    asyncmovieinfo.pixelFormat = pixelFormat;
					asyncmovieinfo.maxNumberThreads = maxNumberThreads;
                    
                    if (windowRecord) {
                        memcpy(&asyncmovieinfo.windowRecord, windowRecord, sizeof(PsychWindowRecordType));
                    } else {
                        memcpy(&asyncmovieinfo.windowRecord, 0, sizeof(PsychWindowRecordType));
                    }

                    asyncmovieinfo.moviehandle = -1;

                    // Increase our scheduling priority to basic RT priority: This way we should get
                    // more cpu time for our PTB main thread than the async. background prefetch-thread:
					// On Windows we must not go higher than basePriority 1 (HIGH PRIORITY) or bad interference can happen.
					// On OS/X we use basePriority 2 for robust realtime, using up to (4+1) == 5 msecs of time in every 10 msecs slice, allowing for up to 1 msec jitter/latency for ops.
					// On Linux we just use standard basePriority 2 RT-FIFO scheduling and trust the os to do the right thing.
                    if ((rc=PsychSetThreadPriority(NULL, ((PSYCH_SYSTEM == PSYCH_WINDOWS) ? 1 : 2), ((PSYCH_SYSTEM == PSYCH_OSX) ? 4 : 0)))!=0) {
                        printf("PTB-WARNING: In OpenMovie(): Failed to raise priority of main thread [System error %i]. Expect movie timing problems.\n", rc);
                    }

                    // Start our own movie loader Posix-Thread:
                    PsychCreateThread(&asyncmovieinfo.pid, NULL, PsychAsyncCreateMovie, &asyncmovieinfo);
                    
                    // Async movie open initiated. We return control to host environment:
                    return(PsychError_none);
                break;
                    
                case 1: // Async open operation in progress, but not yet finished.
                    // Should we wait for completion or just return?
                    if (asyncFlag & 1) {
                        // Async poll requested. We just return -1 to signal that open isn't finished yet:
                        PsychCopyOutDoubleArg(1, TRUE, -1);
                        return(PsychError_none);
                    }
                    // We fall through to case 2 - Wait for "Load operation successfully finished."
                
                case 2: // Async open operation successfully finished. Parse asyncinfo struct and return it to host environment:
                    // We need to join our terminated worker thread to release its ressources. If the worker-thread
                    // isn't done yet (fallthrough from case 1 for sync. wait), this join will block us until worker
                    // completes:
                    PsychDeleteThread(&asyncmovieinfo.pid);

                    asyncmovieinfo.asyncstate = 0; // Reset state to idle:
                    moviehandle = asyncmovieinfo.moviehandle;
                    
                    // Movie successfully opened?
                    if (moviehandle < 0) {
                        // Movie loading failed for some reason.
                        printf("PTB-ERROR: When trying to asynchronously load movie %s, the operation failed: ", asyncmovieinfo.moviename);
                        #if PSYCH_SYSTEM == PSYCH_OSX
                        switch(moviehandle) {
                            case -2000:
                            case -50:
                            case -43:
                                printf("File not found.");
                                break;
                                
                            case -2048:
                                printf("This is not a file that Quicktime understands.");
                                break;
                                
                            case -2003:
                                printf("Can't find media handler (codec) for this movie.");
                                break;
                                
                            case -2:
                                printf("Maximum allowed number of simultaneously open movie files exceeded!");
                                break;
                                
                            case -1:
                                printf("Internal error: Failure in PTB's movie playback engine!");
                                break;
                                
                            default:
                                printf("Unknown error (Quicktime error %i): Check http://developer.apple.com/documentation/QuickTime/APIREF/ErrorCodes.htm#//apple_ref/doc/constant_group/Error_Codes", moviehandle);
                        }
                        printf("\n\n");
                        #endif

                        PsychErrorExitMsg(PsychError_user, "Asynchronous loading of the Quicktime movie failed.");
                    }
                    
                    // We can fall out of the switch statement and continue with the standard synchronous load code as if
                    // the movie had been loaded synchronously.
                break;
                default:
                    PsychErrorExitMsg(PsychError_internal, "Unhandled async movie state condition encountered! BUG!!");
            }
        }

        // Upon sucessfull completion, we'll have a valid handle in 'moviehandle'.
        // Return it to Matlab-world:
        PsychCopyOutDoubleArg(1, TRUE, (double) moviehandle);

        // Retrieve infos about new movie:
        
        // Is the "count" output argument (total number of frames) requested by user?
        if (PsychGetNumOutputArgs() > 5) {
            // Yes. Query the framecount (expensive!) and return it:
            PsychGetMovieInfos(moviehandle, &width, &height, &framecount, &durationsecs, &framerate, NULL, &aspectRatio);
            PsychCopyOutDoubleArg(6, TRUE, (double) framecount);
        }
        else {
            // No. Don't compute and return it.
            PsychGetMovieInfos(moviehandle, &width, &height, NULL, &durationsecs, &framerate, NULL, &aspectRatio);
        }

        PsychCopyOutDoubleArg(2, FALSE, (double) durationsecs);
        PsychCopyOutDoubleArg(3, FALSE, (double) framerate);
        PsychCopyOutDoubleArg(4, FALSE, (double) width);
        PsychCopyOutDoubleArg(5, FALSE, (double) height);
        PsychCopyOutDoubleArg(7, FALSE, (double) aspectRatio);

	// Ready!
	return(PsychError_none);
}

// Functions for movie creation/editing/writing:

PsychError SCREENFinalizeMovie(void)
{
	static char useString[] = "Screen('FinalizeMovie', moviePtr);";
	static char synopsisString[] = "Finish creating a new movie file with handle 'moviePtr' and store it to filesystem.\n";
	static char seeAlsoString[] = "CreateMovie AddFrameToMovie CloseMovie PlayMovie GetMovieImage GetMovieTimeIndex SetMovieTimeIndex";

	int			moviehandle = -1;

	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};
	
	PsychErrorExit(PsychCapNumInputArgs(1));            // Max. 3 input args.
	PsychErrorExit(PsychRequireNumInputArgs(1));        // Min. 2 input args required.
	PsychErrorExit(PsychCapNumOutputArgs(0));           // Max. 1 output args.

	// Get the moviehandle:
	PsychCopyInIntegerArg(1, kPsychArgRequired, &moviehandle);
	
	// Finalize the movie:
	if (!PsychFinalizeNewMovieFile(moviehandle)) {
		printf("See http://developer.apple.com/documentation/QuickTime/APIREF/ErrorCodes.htm#//apple_ref/doc/constant_group/Error_Codes.\n\n");
		PsychErrorExitMsg(PsychError_user, "FinalizeMovie failed for reason mentioned above.");
	}

	return(PsychError_none);
}

PsychError SCREENCreateMovie(void)
{
	static char useString[] = "moviePtr = Screen('CreateMovie', windowPtr, movieFile [, width][, height][, frameRate=30][, movieOptions]);";
	static char synopsisString[] = 
		"Create a new movie file with filename 'movieFile' and according to given 'movieOptions'.\n"
		"The function returns a handle 'moviePtr' to the file.\n"
        "On OS/X and on MS-Windows with Matlab versions prior to R2007a, Apple Quicktime is used for movie writing. "
        "On GNU/Linux and on MS-Windows with recent versions of Matlab or GNU/Octave, the GStreamer multimedia "
        "framework is used for movie writing. GStreamer is generally more advanced and offers more functionality. "
		"Currently only single-track video encoding is supported. Audio encoding is not supported on OS/X.\n"
        "See 'Screen AddAudioBufferToMovie?' on how to add audio tracks to movies via GStreamer.\n"
        "\n"
		"Movie creation is a 3 step procedure:\n"
		"1. Create a movie and define encoding options via 'CreateMovie'.\n"
		"2. Add video and audio data to the movie via calls to 'AddFrameToMovie' et al.\n"
		"3. Finalize and close the movie via a call to 'FinalizeMovie'.\n\n"
		"All following parameters are optional and have reasonable defaults:\n\n"
		"'width' Width of movie video frames in pixels. Defaults to width of window 'windowPtr'.\n"
		"'height' Height of movie video frames in pixels. Defaults to height of window 'windowPtr'.\n"
		"'frameRate' Playback framerate of movie. Defaults to 30 fps. Technically this is not the "
		"playback framerate but the granularity in 1/frameRate seconds with which the duration of "
		"a single movie frame can be specified. When you call 'AddFrameToMovie', there's an optional "
		"parameter 'frameDuration' which defaults to one. The parameter defines the display duration "
		"of that frame as the fraction 'frameDuration' / 'frameRate' seconds, so 'frameRate' defines "
		"the denominator of that term. However, for a default 'frameDuration' of one, this is equivalent "
		"to the 'frameRate' of the movie, at least if you leave everything at defaults.\n\n"
		"'movieoptions' a textstring which allows to define additional parameters via keyword=parm pairs. "
        "For GStreamer movie writing, you can provide the same options as for GStreamer video recording. "
        "See 'help VideoRecording' for supported options and tips.\n"
		"Keywords unknown to a certain implementation or codec will be silently ignored:\n"
		"EncodingQuality=x Set encoding quality to value x, in the range 0.0 for lowest movie quality to "
		"1.0 for highest quality. Default is 0.5 = normal quality. 1.0 usually provides lossless encoding.\n"
        "On systems with Quicktime movie writing, codecs can be selected by FOURCC codes. This is not supported "
        "on GStreamer setups for now:\n"
		"CodecFOURCCId=id FOURCC id. The FOURCC of a desired video codec as a number. Defaults to H.264 codec.\n"
		"Choice of codec and quality defines a tradeoff between filesize, quality, processing demand and speed, "
		"as well as on which target devices you'll be able to play your movie.\n"
		"CodecFOURCC=xxxx FOURCC as a four character text string instead of a number.\n"
		"\n";

	static char seeAlsoString[] = "FinalizeMovie AddFrameToMovie CloseMovie PlayMovie GetMovieImage GetMovieTimeIndex SetMovieTimeIndex";
	
	PsychWindowRecordType					*windowRecord;
	char                                    *moviefile;
	char									*movieOptions;
	int                                     moviehandle = -1;
	double                                  framerate = 30.0;
	int                                     width;
	int                                     height;
	char									defaultOptions[2] = "";
	
	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};
	
	PsychErrorExit(PsychCapNumInputArgs(6));            // Max. 6 input args.
	PsychErrorExit(PsychRequireNumInputArgs(2));        // Min. 2 input args required.
	PsychErrorExit(PsychCapNumOutputArgs(1));           // Max. 1 output args.
	
	// Get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
	// Only onscreen windows allowed:
	if(!PsychIsOnscreenWindow(windowRecord)) {
		PsychErrorExitMsg(PsychError_user, "CreateMovie called on something else than an onscreen window.");
	}
	
	// Get the movie name string:
	moviefile = NULL;
	PsychAllocInCharArg(2, kPsychArgRequired, &moviefile);
	
	// Get the optional size:
	// Default Width and Height of movie frames is derived from size of window:
	width  = (int) PsychGetWidthFromRect(windowRecord->clientrect);
	height = (int) PsychGetHeightFromRect(windowRecord->clientrect);
	PsychCopyInIntegerArg(3, kPsychArgOptional, &width);
	PsychCopyInIntegerArg(4, kPsychArgOptional, &height);
	
	// Get the optional framerate:
	PsychCopyInDoubleArg(5, kPsychArgOptional, &framerate);
	
	// Get the optional options string:
	movieOptions = defaultOptions;
	PsychAllocInCharArg(6, kPsychArgOptional, &movieOptions);

	// Create movie of given size and framerate with given options:
	moviehandle = PsychCreateNewMovieFile(moviefile, width, height, framerate, movieOptions);
	if (0 > moviehandle) {
		#ifndef PTB_USE_GSTREAMER
			printf("See http://developer.apple.com/documentation/QuickTime/APIREF/ErrorCodes.htm#//apple_ref/doc/constant_group/Error_Codes.\n\n");
		#endif
		PsychErrorExitMsg(PsychError_user, "CreateMovie failed for reason mentioned above.");
	}
	
	// Return handle to it:
	PsychCopyOutDoubleArg(1, FALSE, (double) moviehandle);
	
	return(PsychError_none);
}

PsychError SCREENAddAudioBufferToMovie(void)
{
	static char useString[] = "Screen('AddAudioBufferToMovie', moviePtr, audioBuffer);";
	static char synopsisString[] = 
		"Add a buffer filled with audio data samples to movie 'moviePtr'.\n"
        "This function is only supported with the GStreamer based movie writing functions. "
        "It doesn't work on MS-Windows with Matlab versions before R2007a and it doesn't work "
        "on Apple OS/X yet.\n"
        "The movie must have been created in 'CreateMovie' with an options string that "
        "enables writing of an audio track into the movie, otherwise this function will fail.\n"
        "You enable writing of audio tracks by adding the keyword 'AddAudioTrack' to the options string.\n"
        "Alternatively, if your options string is a gst-launch style pipeline description, it must contain "
        "one pipeline element with a name option of 'name=ptbaudioappsrc'.\n"
        "'audioBuffer' must be 'numChannels' rows by 'numSamples' columns double matrix of audio data. "
        "Each row encodes one audio channel, each column element in a row encodes a sample. "
        "E.g., a 2-by-48000 matrix would encode 48000 samples for a two channel stereo sound track.\n"
        "Sample values must lie in the range between -1.0 and +1.0.\n"
        "The audio buffer is converted into a movie specific sound format and then appended to "
        "the audio samples already stored in the audio track.\n"
		"\n";

	static char seeAlsoString[] = "FinalizeMovie AddFrameToMovie CloseMovie PlayMovie GetMovieImage GetMovieTimeIndex SetMovieTimeIndex";
	
	int     moviehandle = -1;
    int     m, n, p;
    double* buffer;
	
	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};
	
	PsychErrorExit(PsychCapNumInputArgs(2));            // Max. 2 input args.
	PsychErrorExit(PsychRequireNumInputArgs(2));        // Min. 2 input args required.
	PsychErrorExit(PsychCapNumOutputArgs(0));           // Max. 0 output args.
	
    // Get movie handle:
	PsychCopyInIntegerArg(1, kPsychArgRequired, &moviehandle);
    
    // And audio date buffer:
	PsychAllocInDoubleMatArg(2, kPsychArgRequired, &m, &n, &p, &buffer);
    if (p!=1 || m < 1 || n < 1) PsychErrorExitMsg(PsychError_user, "Invalid audioBuffer provided. Must be a 2D matrix with at least one row and at least one column!");

    // Pass audio data to movie writing engine:
    PsychAddAudioBufferToMovie(moviehandle, m, n, buffer);
    
	return(PsychError_none);
}
