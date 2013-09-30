/*
Psychtoolbox3/Source/Common/SCREENGetMovieImage.c		
 
 AUTHORS:
 mario.kleiner at tuebingen.mpg.de   mk
 
 PLATFORMS:	
 This file should build on any platform. 
 
 HISTORY:
 10/23/05  mk		Created. 
 
 DESCRIPTION:
 
 Fetch an image from the specified movie object and create an OpenGL texture out of it.
 Return a handle to the texture.

 TO DO:
 
 */

#include "Screen.h"

static char useString[] = "[ texturePtr [timeindex]]=Screen('GetMovieImage', windowPtr, moviePtr, [waitForImage=1], [fortimeindex], [specialFlags = 0] [, specialFlags2 = 0]);";
static char synopsisString[] = 
"Try to fetch a new texture image from movie object 'moviePtr' for visual playback/display in onscreen window 'windowPtr' and "
"return a texture-handle 'texturePtr' on successfull completion. 'waitForImage' If set to 1 (default), the function will wait "
"until the image becomes available. If set to zero, the function will just poll for a new image. If none is ready, it will return "
"a texturePtr of zero, or -1 if none will become ready because movie has reached its end and is not in loop mode.\n"
"'fortimeindex' Don't request the next image, but the image closest to time 'fortimeindex' in seconds. 'fortimeindex' is only used "
"if either playback is stopped in any case, or if forward playback is active while using the GStreamer playback engine and "
"fetching video frames asynchronously due to opening the movie with async flag set and sufficient buffering enabled on suitable "
"movie file formats. In all other cases the 'fortimeindex' is silently ignored. Therefore you should always check if the timestamp "
"of the returned movie frame actually satisfies the given 'fortimeindex' value if you choose to use 'fortimeindex'.\n"
"The (optional) return value 'timeindex' contains the exact time when the returned image should be displayed wrt. to the "
"start of the movie - a presentation timestamp. \n"
"'specialFlags' (optional) encodes special requirements for the returned texture. A setting of 1 will try to create the texture "
"as a GL_TEXTURE_2D texture. However this is not well supported with a setting of 2 (see below) or with use of a special 'pixelFormat' > 4 in "
"Screen('OpenMovie'). A setting of 2 will request that the texture "
"is prepared for drawing it with highest possible spatial precision. See explanation of 'specialFlags' == 2 in Screen('MakeTexture') "
"for details. A 'specialFlags' == 8 will prevent automatic mipmap-generation for GL_TEXTURE_2D textures. A 'specialflags' == 32 "
"will prevent closing the texture by a call to Screen('Close');\n"
"'specialFlags2' (optional) More special flags: A setting of 1 tells the function to not return presentation timestamps in 'timeindex'. "
"This means that 'timeindex' will be always returned as zero and that the built-in detector for skipped frames is disabled as well. This "
"may (or may not) save a little bit of computation time during playback of very demanding movies on lower end systems, your mileage will "
"vary, depending on many factors.\n"
"A setting of 2 will skip creation and return of an actual video texture, instead a texture handle of 1 will be returned and no texture "
"gets created. This is useful for fast fine-grained forward seeking in the movie by skipping single frames, and for benchmarking.\n"
;

static char seeAlsoString[] = "CloseMovie PlayMovie GetMovieImage GetMovieTimeIndex SetMovieTimeIndex";

PsychError SCREENGetMovieImage(void) 
{
    PsychWindowRecordType	*windowRecord;
    PsychWindowRecordType	*textureRecord;
    PsychRectType		rect;
    double			deadline, tnow;
    int                         moviehandle = -1;
    int                         waitForImage = TRUE;
    double                      requestedTimeIndex = -1;
    double                      presentation_timestamp = 0;
    int				rc=0;
    int				specialFlags = 0;
    int				specialFlags2 = 0;
	
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};
    
    PsychErrorExit(PsychCapNumInputArgs(6));            // Max. 6 input args.
    PsychErrorExit(PsychRequireNumInputArgs(2));        // Min. 2 input args required.
    PsychErrorExit(PsychCapNumOutputArgs(2));           // Max. 2 output args.
    
    // Get the window record from the window record argument and get info from the window record
    PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
    // Only onscreen windows allowed:
    if(!PsychIsOnscreenWindow(windowRecord)) {
        PsychErrorExitMsg(PsychError_user, "GetMovieImage called on something else than an onscreen window.");
    }
    
    // Get the movie handle:
    PsychCopyInIntegerArg(2, TRUE, &moviehandle);
    if (moviehandle==-1) {
        PsychErrorExitMsg(PsychError_user, "GetMovieImage called without valid handle to a movie object.");
    }

    // Get the 'waitForImage' flag: If waitForImage == true == 1, we'll do a blocking wait for
    // arrival of a new image for playback. Otherwise we will return with a 0-Handle if there
    // isn't any new image available.
    PsychCopyInIntegerArg(3, FALSE, &waitForImage);
    
    // Get the requested timeindex for the frame. The default is -1, which means: Get the next image,
    // according to current movie playback time.
    PsychCopyInDoubleArg(4, FALSE, &requestedTimeIndex);
    
    // Get the optional specialFlags flag:
    PsychCopyInIntegerArg(5, FALSE, &specialFlags);

    // Get the optional specialFlags2 flag:
    PsychCopyInIntegerArg(6, FALSE, &specialFlags2);

    PsychGetAdjustedPrecisionTimerSeconds(&deadline);
    deadline += 5;

    while (rc==0) {
		// With some backends (GStreamer), we can use a checkForImage of type 2 if a blocking wait is
		// requested. Type 2 actually blocks inside the backend, to minimize cpu load and achieve lowest signalling delay
		// once a new frame becomes available. Type 1 only polls.
        rc = PsychGetTextureFromMovie(windowRecord, moviehandle, (waitForImage != 0) ? 2 : 1, requestedTimeIndex, NULL, NULL);
		PsychGetAdjustedPrecisionTimerSeconds(&tnow);
        if (rc<0 || ((tnow > deadline) && (waitForImage != 0))) {
            // No image available and there won't be any in the future, because the movie has reached
            // its end and we are not in looped playback mode:
            if (tnow > deadline) printf("PTB-ERROR: In Screen('GetMovieImage') for movie %i: Timed out while waiting for new frame after 5 seconds!\n", moviehandle);

            // No new texture available: Return a negative handle:
            PsychCopyOutDoubleArg(1, TRUE, -1);
            // ...and an invalid timestamp:
            PsychCopyOutDoubleArg(2, FALSE, -1);
            // Ready!
            return(PsychError_none);
        }
        else if (rc==0 && waitForImage == 0) {
            // We should just poll once and no new texture available: Return a null-handle:
            PsychCopyOutDoubleArg(1, TRUE, 0);
            // ...and an invalid timestamp:
            PsychCopyOutDoubleArg(2, FALSE, -1);
            // Ready!
            return(PsychError_none);
        }
        else if (rc==0 && waitForImage != 0) {
            // No new texture available yet. Just sleep a bit and then retry...
			// This is thankfully only needed / used for the deprecated QT backend.
            PsychYieldIntervalSeconds(0.001);
        }
    }

    // New image available: Go ahead...
    
    // Skipping of actual texture creation requested via specialFlags2 setting 2?
    if (specialFlags2 & 2) {
        // Yes. We skip OpenGL texture creation and just return some dummy value:

        // Still need to do a dummy-fetch from the movie object to retrieve a potential
        // presentation timestamp and to make sure the videobuffers get properly dequeued, so
        // the playback engine stays happy. We pass a textureRecord pointer of NULL to tell the
        // engine to skip most work:
        PsychGetTextureFromMovie(windowRecord, moviehandle, FALSE, requestedTimeIndex, NULL, ((specialFlags2 & 1) ? NULL : &presentation_timestamp));
        
        // Return positive pseudo-texture handle:
        PsychCopyOutDoubleArg(1, TRUE, 1);
        
        // Return presentation timestamp for this image:
        PsychCopyOutDoubleArg(2, FALSE, presentation_timestamp);
        
        // Ready!
        return(PsychError_none);        
    }
    
    // Create a texture record.  Really just a window record adapted for textures.  
    PsychCreateWindowRecord(&textureRecord);	// This also fills the window index field.
    // Set mode to 'Texture':
    textureRecord->windowType=kPsychTexture;
    // We need to assign the screen number of the onscreen-window.
    textureRecord->screenNumber=windowRecord->screenNumber;
    // It is always a 32 bit texture for movie textures:
    textureRecord->depth=32;
    textureRecord->nrchannels = 4;

    // Create default rectangle which describes the dimensions of the image. Will be overwritten
    // later on.
    PsychMakeRect(rect, 0, 0, 10, 10);
    PsychCopyRect(textureRecord->rect, rect);
    
    // Other setup stuff:
    textureRecord->textureMemorySizeBytes= 0;
    textureRecord->textureMemory=NULL;

    // Assign parent window and copy its inheritable properties:
    PsychAssignParentWindow(textureRecord, windowRecord);

    // Special texture format requested?
    if (specialFlags & 0x1) textureRecord->texturetarget = GL_TEXTURE_2D;
    
    // specialFlags setting 8? Disable auto-mipmap generation:
    if (specialFlags & 0x8) textureRecord->specialflags |= kPsychDontAutoGenMipMaps;    

    // A specialFlags setting of 32? Protect texture against deletion via Screen('Close') without providing a explicit handle:
    if (specialFlags & 32) textureRecord->specialflags |= kPsychDontDeleteOnClose;    

    // Try to fetch an image from the movie object and return it as texture:
    PsychGetTextureFromMovie(windowRecord, moviehandle, FALSE, requestedTimeIndex, textureRecord, ((specialFlags2 & 1) ? NULL : &presentation_timestamp));

    // Assign GLSL filter-/lookup-shaders if needed: usefloatformat is always == 0 as
    // our current movie engine implementations only return 8 bpc fixed textures.
    // The 'userRequest' flag is set if specialmode flag is set to 8.
    PsychAssignHighPrecisionTextureShaders(textureRecord, windowRecord, 0, (specialFlags & 2) ? 1 : 0);

    // Texture ready for consumption. Mark it valid and return handle to userspace:
    PsychSetWindowRecordValid(textureRecord);
    PsychCopyOutDoubleArg(1, TRUE, textureRecord->windowIndex);

    // Return presentation timestamp for this image:
    PsychCopyOutDoubleArg(2, FALSE, presentation_timestamp);
    
    // Ready!
    return(PsychError_none);
}
