/*
Psychtoolbox3/Source/Common/SCREENGetCapturedImage.c		
 
 AUTHORS:
 mario.kleiner at tuebingen.mpg.de   mk
 
 PLATFORMS:	
 This file should build on any platform. 
 
 HISTORY:
 2/7/06  mk		Created. 
 
 DESCRIPTION:
 
 Fetch an image from the specified capture object and create an OpenGL texture out of it.
 Return a handle to the texture.
 
 TO DO:

 */

#include "Screen.h"

static char useString[] = "[ texturePtr [capturetimestamp] [droppedcount] [summed_intensityOrRawImageMatrix]]=Screen('GetCapturedImage', windowPtr, capturePtr [, waitForImage=1] [,oldTexture] [,specialmode] [,targetmemptr]);";
static char synopsisString[] = 
"Try to fetch a new image from video capture device 'capturePtr' for visual playback/display in window 'windowPtr' and "
"return a texture-handle 'texturePtr' on successfull completion. 'waitForImage' If set to 1 (default), the function will wait "
"until the image becomes available. If set to zero, the function will just poll for a new image. If none is ready, it will return "
"a texturePtr of zero, or -1 if none will become ready because capture has been stopped. Set 'waitForImage' = 2 if you just want to get the "
"summed intensity or timestamp of the image. 'waitForImage' = 3 will behave as a setting of 2, but it will only poll, not block. "
"Note: Settings 2 and 3 will not return textures, therefore the returned 'texturePtr' will always be zero, irrespective of availabilty "
"of new frames. You'll have to check the other return arguments for indication that a frame has arrived, e.g., non-zero summed_intensity.\n"
"The setting "
"'waitForImage' = 4 is special: It will not block, nor will it return any information, but will make sure the capture engine keeps running. "
"This is useful if you want PTB to perform harddisk recording of video footage in the background without any need to access information about "
"the capture process or the actual video data from within Matlab. This is the least ressource consuming way of doing this.\n"
"'oldTexture' (if provided) allows you to pass the texture handle of an already existing texture. "
"Psychtoolbox will reuse that texture by overwriting its previous image content with the new image, instead of creating a completely new "
"texture for the new image. Use of this ''recycling facility'' may allow for higher capture framerates and lower latencies on low-end "
"graphics hardware in some cases. Providing a value of 'oldTexture'=0 is the same as leaving it out. The optional argument 'specialmode' "
"allows to request special treatment of textures. Currently, specialmode = 1 will ask PTB to use GL_TEXTURE_2D textures instead of other "
"formats. This is sometimes less efficient, unless you want to do realtime blurring of images. If you set specialmode = 2, then "
"the optional return value 'summed_intensityOrRawImageMatrix' will not return the summed pixel intensity, but a Matlab uint8 Matrix with "
"the captured raw image data for direct use within Matlab, e.g., via the image processing toolbox. If you set 'specialmode = 4 and "
"provide a double-encoded memory pointer in 'targetmemptr', then PTB will copy the raw image data into that buffer. The buffer is "
"expected to be of sufficient size, otherwise a crash will occur (Experts only!).\n"
"A 'specialmode' == 8 will require high-precision drawing, see the specialFlag == 2 setting in Screen('MakeTexture') for a "
"description of its meaning. A 'specialMode' == 16 will prevent automatic mipmap-generation for GL_TEXTURE_2D textures. "
"A 'specialMode' == 32 will prevent closing of the texture by a call to Screen('Close');\n"
"'capturetimestamp' contains the system time when the returned image was captured. This timestamp has been verified to "
"be very precise on Linux with suitable professional IIDC 1394 firewire cameras when the dc1394 capture engine is used. "
"The same may be true for OS/X, although this hasn't been extensively tested. If other operating systems, capture engines "
"or cameras are used, the timestamp may become just a rough approximation of unknown precision.\n"
"The (optional) return value 'droppedcount' contains the "
"number of captured frames that had to be dropped to keep in sync with realtime or due to internal shortage of buffer memory. "
"If you didn't specify the 'dropFrames' flag in Screen('StartVideoCapture') then it will report the number of pending "
"buffers which can be fetched, i.e., how many more buffers are queued up for delivery.\n"
"The (optional) return value 'summed_intensityOrRawImageMatrix' contains the sum of all pixel intensity values of all channels of the image - some measure of overall brightness. "
"Only query this value if you really need it, its computation is time consuming.";

static char seeAlsoString[] = "CloseVideoCapture StartVideoCapture StopVideoCapture GetCapturedImage";

PsychError SCREENGetCapturedImage(void) 
{
    PsychWindowRecordType		*windowRecord;
    PsychWindowRecordType		*textureRecord;
    PsychRectType				rect;
    double                      summed_intensity;
    int                         capturehandle = -1;
    int                         waitForImage = TRUE;
    int                         specialmode = 0;
	double						timeout, tnow;
    double                      presentation_timestamp = 0;
    int							rc=-1;
    double						targetmemptr = 0;
	double*						tsummed = NULL;
	psych_uint8					*targetmatrixptr = NULL;
	static rawcapimgdata		rawCaptureBuffer = {0, 0, 0, NULL};

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};
    
    PsychErrorExit(PsychCapNumInputArgs(6));            // Max. 6 input args.
    PsychErrorExit(PsychRequireNumInputArgs(2));        // Min. 2 input args required.
    PsychErrorExit(PsychCapNumOutputArgs(4));           // Max. 4 output args.
    
    // Get the window record from the window record argument and get info from the window record
    PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
    // Only onscreen windows allowed:
    if(!PsychIsOnscreenWindow(windowRecord) && !PsychIsOffscreenWindow(windowRecord)) {
        PsychErrorExitMsg(PsychError_user, "GetCapturedImage called on something else than an onscreen window or offscreen window.");
    }
    
    // Get the handle:
    PsychCopyInIntegerArg(2, TRUE, &capturehandle);
    if (capturehandle==-1) {
        PsychErrorExitMsg(PsychError_user, "GetCapturedImage called without valid handle to a capture object.");
    }

    // Get the 'waitForImage' flag: If waitForImage == true == 1, we'll do a blocking wait for
    // arrival of a new image. Otherwise we will return with a 0-Handle if there
    // isn't any new image available.
    PsychCopyInIntegerArg(3, FALSE, &waitForImage);

	// Special case waitForImage == 4? This would ask to call into the capture driver, but
	// not wait for any image to arrive and not return any information. This is only useful
	// on OS/X and Windows when using the capture engine for video recording to harddisk. In
	// that case we are not interested at all in the captured live video, we just want it to
	// get written to harddisk in the background. To keep the video encoder going, we need to
	// call its SGIdle() routine and waitForImage==4 does just that, call SGIdle().
	if (waitForImage == 4) {
		// Perform the null-call to the capture engine, ie a SGIdle() on OS/X and Windows:
		PsychGetTextureFromCapture(windowRecord, capturehandle, 4, 0.0, NULL, NULL, NULL, NULL);
		// Done. Nothing to return...
		return(PsychError_none);
	}
	
    // Get the optional textureRecord for the optional texture handle. If the calling script
    // provides the texture handle of an existing Psychtoolbox texture that has a matching
    // format, then that texture is recycled by overwriting its previous content with the
    // image data from the new captured image. This can save some overhead for texture destruction
    // and recreation. While this is probably not noticeable on mid- to high-end gfx cards with
    // rectangle texture support, it can provide a significant speedup on low-end gfx cards with
    // only power-of-two texture support.
    textureRecord = NULL;
    if ((PsychGetNumInputArgs()>=4) && PsychIsWindowIndexArg(4)) PsychAllocInWindowRecordArg(4, FALSE, &textureRecord);
    
    // Get the optional specialmode flag:
    PsychCopyInIntegerArg(5, FALSE, &specialmode);

	// Set a 10 second maximum timeout for waiting for new frames:
	PsychGetAdjustedPrecisionTimerSeconds(&timeout);
	timeout+=10;

    while (rc==-1) {		
      // We pass a checkForImage value of 2 if waitForImage>0. This way we can signal if we are in polling or blocking mode.
      // With the libdc1394 engine this allows to do a real blocking wait in the driver -- much more efficient than the spin-waiting approach!
      rc = PsychGetTextureFromCapture(windowRecord, capturehandle, ((waitForImage>0 && waitForImage<3) ? 2 : 1), 0.0, NULL, &presentation_timestamp, NULL, &rawCaptureBuffer);
		PsychGetAdjustedPrecisionTimerSeconds(&tnow);
        if (rc==-2 || (tnow > timeout)) {
            // No image available and there won't be any in the future, because capture has been stopped or there is a timeout:
			if (tnow > timeout) printf("PTB-WARNING: In Screen('GetCapturedImage') timed out waiting for a new frame. No video data in over 10 seconds!\n");

            // No new texture available: Return a negative handle:
            PsychCopyOutDoubleArg(1, TRUE, -1);
            // ...and an invalid timestamp:
            PsychCopyOutDoubleArg(2, FALSE, -1);
            PsychCopyOutDoubleArg(3, FALSE, 0);
            PsychCopyOutDoubleArg(4, FALSE, 0);

            // Ready!
            return(PsychError_none);
        }
        else if (rc==-1 && (waitForImage == 0 || waitForImage == 3)) {
            // We should just poll once and no new texture available: Return a null-handle:
            PsychCopyOutDoubleArg(1, TRUE, 0);
            // ...and the current timestamp:
            PsychCopyOutDoubleArg(2, FALSE, presentation_timestamp);
            PsychCopyOutDoubleArg(3, FALSE, 0);
            PsychCopyOutDoubleArg(4, FALSE, 0);

            // Ready!
            return(PsychError_none);
        }
        else if (rc==-1 && waitForImage != 0) {
            // No new texture available yet. Just sleep a bit and then retry...
            PsychYieldIntervalSeconds(0.002);
        }
    }

    // rc == 0 --> New image available: Go ahead...
    if (waitForImage!=2 && waitForImage!=3) {
      // Ok, we need a texture for the image. Did script provide an old one for recycling?
      if (textureRecord) {
	// Old texture provided for reuse? Some basic sanity check: Everything else is
	// up to the lower level PsychGetTextureFromCapture() routine.
        if(!PsychIsOffscreenWindow(textureRecord)) {
	  PsychErrorExitMsg(PsychError_user, "GetCapturedImage provided with something else than a texture as fourth call parameter.");
        }	
      }
      else {
        // No old texture provided: Create a new texture record:
        PsychCreateWindowRecord(&textureRecord);

        // Set mode to 'Texture':
        textureRecord->windowType=kPsychTexture;

        // We need to assign the screen number of the onscreen-window.
        textureRecord->screenNumber=windowRecord->screenNumber;

        // It defaults to a 32 bit texture for captured images. On Linux, this will be overriden,
		// if optimized formats exist for our purpose:
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

        // Set textureNumber to zero, which means "Not cached, do not recycle"
        // Todo: Texture recycling like in PsychMovieSupport for higher efficiency!
        textureRecord->textureNumber = 0;
      }

      // Power-of-two texture requested?
      if (specialmode & 0x01) {
		// Yes. Spec it:
		textureRecord->texturetarget = GL_TEXTURE_2D;
      }
    }
    else {
        // Just want to return summed_intensity and timestamp, not real texture...
        textureRecord = NULL;
    }

	// Default to no calculation of summed image intensity:
	tsummed = NULL;
    if ((PsychGetNumOutputArgs() > 3) && !(specialmode & 0x2)) {
        // Return sum of pixel intensities for all channels of this image: Need to
		// assign the output pointer for this to happen:
		tsummed = &summed_intensity;
	}

    // Try to fetch an image from the capture object and return it as texture:
	targetmatrixptr = NULL;
	
	// Shall we return a Matlab matrix?
	if ((PsychGetNumOutputArgs() > 3) && (specialmode & 0x2)) {
		// We shall return a matrix with raw image data. Allocate a uint8 matrix
		// of sufficient size:
		PsychAllocOutUnsignedByteMatArg(4, TRUE, rawCaptureBuffer.depth, rawCaptureBuffer.w, rawCaptureBuffer.h, &targetmatrixptr);
		tsummed = NULL;
	}
	
	// Shall we return data into preallocated memory buffer?
	if (specialmode & 0x4) {
		// Copy in memory address (which itself is encoded in a double value):
		PsychCopyInDoubleArg(6, TRUE, &targetmemptr);
		targetmatrixptr = (psych_uint8*) PsychDoubleToPtr(targetmemptr);
	}
	
	if (targetmatrixptr == NULL) {
		// Standard fetch of a texture and its timestamp:
		rc = PsychGetTextureFromCapture(windowRecord, capturehandle, 0, 0.0, textureRecord, &presentation_timestamp, tsummed, NULL);
	}
	else {
		// Fetch of a memory raw image buffer + timestamp + possibly a texture:
		rawCaptureBuffer.data = (void*) targetmatrixptr;
		rc = PsychGetTextureFromCapture(windowRecord, capturehandle, 0, 0.0, textureRecord, &presentation_timestamp, tsummed, &rawCaptureBuffer);			
	}
	
    if (tsummed) {
        // Return sum of pixel intensities for all channels of this image:
        PsychCopyOutDoubleArg(4, FALSE, summed_intensity);
    }

    // Real texture requested?
    if (textureRecord) {
        // Texture ready for consumption.

		// Assign GLSL filter-/lookup-shaders if needed: usefloatformat is always == 0 as
		// our current capture engine implementations only return 8 bpc fixed textures.
		// The 'userRequest' flag is set if specialmode flag is set to 8.
		PsychAssignHighPrecisionTextureShaders(textureRecord, windowRecord, 0, (specialmode & 8) ? 1 : 0);

        // specialmode setting 16? Disable auto-mipmap generation:
        if (specialmode & 16) textureRecord->specialflags |= kPsychDontAutoGenMipMaps;    

        // A specialFlags setting of 32? Protect texture against deletion via Screen('Close') without providing a explicit handle:
        if (specialmode & 32) textureRecord->specialflags |= kPsychDontDeleteOnClose;    

		// Mark it valid and return handle to userspace:
        PsychSetWindowRecordValid(textureRecord);
        PsychCopyOutDoubleArg(1, TRUE, textureRecord->windowIndex);
    }
    else {
        PsychCopyOutDoubleArg(1, TRUE, 0);
    }
    
    // Return presentation timestamp for this image:
    PsychCopyOutDoubleArg(2, FALSE, presentation_timestamp);

    // Return count of pending frames in buffers or of dropped frames:
    PsychCopyOutDoubleArg(3, FALSE, (double) rc);

    // Ready!
    return(PsychError_none);
}

