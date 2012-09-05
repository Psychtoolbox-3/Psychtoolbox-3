/*
	SCREENGetImage.c		
  
	AUTHORS:
		Allen.Ingling@nyu.edu			awi 
		mario.kleiner@tuebingen.mpg.de	mk

	PLATFORMS:	
	
		All.
    
	HISTORY:
	
		01/08/03  	awi		Created.
		10/12/04	awi		In useString: moved commas to inside [].
		03/20/11	mk		Made 64-bit clean.

	TO DO:
    
		We could probably speed this up by using a single call to glReadPixels instead of three.  It doesn't matter what we pack the result into,
		since we have to iterate over and copy out from the memory which glReadPixels fills so to rearrange correclty into the return matrix. None
		of the glReadPixels modes will fill the return matrix in the same format that MATLAB wants.  
*/



#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] =  "imageArray=Screen('GetImage', windowPtr [,rect] [,bufferName] [,floatprecision=0] [,nrchannels=3])";
//                                                        1           2       3				4				   5

static char synopsisString[] =
"Slowly copy an image from a window or texture to Matlab/Octave, by default returning a uint8 array.\n\n"
"Calling this function on an onscreen window while an asynchronous flip is pending on "
"the window due to Screen('AsyncFlipBegin') is not allowed! Finalize such flips first. "
"Readback of other onscreen or offscreen windows or textures is possible during async "
"flip, but discouraged because it will have a significant impact on performance.\n\n"
"The returned imageArray by default has three layers, i.e. it's an RGB image. "
"\"windowPtr\" is the handle of the onscreen window, offscreen window or texture whose image "
"should be returned. "
"\"rect\" is the rectangular subregion to copy, and its default is the whole window. "
"Matlab/Octave will complain if you try to do math on a uint8 array, so you may need "
"to use DOUBLE to convert it, e.g. imageArray/255 will produce an error, but "
"double(imageArray)/255 is ok. Also see Screen 'PutImage' and 'CopyWindow'.\n\n"
"\"bufferName\" is a string specifying the buffer from which to copy the image: "
"The 'bufferName' argument is meaningless for offscreen windows and textures and "
"will be silently ignored. For onscreen windows, it defaults to 'frontBuffer', i.e., "
"it returns the image that your subject would see at that moment. If frame-sequential "
"stereo mode is enabled, 'frontLeftBuffer' returns what your subject would see in its "
"left eye, 'frontRightBuffer' returns the subjects right-eye view. If double-buffering "
"is enabled, you can also return the 'backBuffer', i.e. what your subject will see after "
"the next Screen('Flip') command, and for frame-sequential stereo also 'backLeftBuffer' "
"and 'backRightBuffer' respectively. 'aux0Buffer' - 'aux3Buffer' returns the content of "
"OpenGL AUX buffers 0 to 3. Only query the AUX buffers if you know what you are doing, "
"otherwise your script will crash. This is mostly meant for internal debugging of PTB. "
"If the imaging pipeline is enabled, you can also return the content of the unprocessed "
"backbuffer, ie. before processing by the pipeline, by requesting 'drawBuffer'.\n"
"If the imaging pipeline is enabled, querying the 'backBuffer' will only give you up to date "
"results after you called Screen('Flip') or Screen('AsyncFlipBegin') or Screen('DrawingFinished'), "
"otherwise you may encounter stale results. If our own homegrown frame-sequential stereo mode "
"is in use, querying the 'backLeftBuffer' and 'backRightBuffer' is well defined, after a Screen('Flip') "
"or Screen('DrawingFinished'), but querying the 'frontLeftBuffer' or 'frontRightBuffer' or "
"'frontBuffer' will result in an assignment of left- and right- stereo images that is mostly "
"based on chance, e.g., you may get accidentally swapped left- and rightBuffer images!\n\n"
"\"floatprecision\" If you set this optional flag to 1, the image data will be returned "
"as a double precision matrix instead of a uint8 matrix. Please note that normal image "
"data will be returned in the normalized range 0.0 to 1.0 instead of 0 - 255. Floating "
"point readback is only beneficial when reading back floating point precision textures, "
"offscreen windows or the framebuffer when the imaging pipeline is active and HDR mode "
"is selected (ie. more than 8bpc framebuffer).\n"
"\"nrchannels\" Number of color channels to return. By default, 3 channels (RGB) are "
"returned. Specify 1 for Red/Luminance only, 2 for Red+Green or Luminance+Alpha, 3 for "
"RGB and 4 for RGBA.\n\n";

static char useString2[] = "Screen('AddFrameToMovie', windowPtr [,rect] [,bufferName] [,moviePtr=0] [,frameduration=1])";
//                                                    1           2       3				4			  5

static char synopsisString2[] =
"Get an image from a window or texture and add it as a new video frame to a movie.\n\n"
"Calling this function on an onscreen window while an asynchronous flip is pending on "
"the window due to Screen('AsyncFlipBegin') is not allowed! Finalize such flips first. "
"Readback of other onscreen or offscreen windows or textures is possible during async "
"flip, but discouraged because it will have a significant impact on performance.\n\n"
"\"windowPtr\" is the handle of the onscreen window, offscreen window or texture whose image "
"should be added.\n\n"
"\"rect\" is the rectangular subregion to get, and its default is the whole window. "
"The function will only honor the top-left corner of the 'rect' argument, but ignore "
"the bottom-right corner and thereby the width and height as defined by the 'rect'. "
"The size of added video frames is fully defined as the fixed size of movie frames as "
"specified at movie creation time via the Screen('CreateMovie') call.\n\n"
"\"bufferName\" is a string specifying the buffer from which to copy the image: "
"The 'bufferName' argument is meaningless for offscreen windows and textures and "
"will be silently ignored. For onscreen windows, it defaults to 'frontBuffer', i.e., "
"it returns the image that your subject would see at that moment. If frame-sequential "
"stereo mode is enabled, 'frontLeftBuffer' returns what your subject would see in its "
"left eye, 'frontRightBuffer' returns the subjects right-eye view. If double-buffering "
"is enabled, you can also return the 'backBuffer', i.e. what your subject will see after "
"the next Screen('Flip') command, and for frame-sequential stereo also 'backLeftBuffer' "
"and 'backRightBuffer' respectively. 'aux0Buffer' - 'aux3Buffer' returns the content of "
"OpenGL AUX buffers 0 to 3. Only query the AUX buffers if you know what you are doing, "
"otherwise your script will crash. This is mostly meant for internal debugging of PTB. "
"If the imaging pipeline is enabled, you can also return the content of the unprocessed "
"backbuffer, ie. before processing by the pipeline, by requesting 'drawBuffer'. \n\n"
"If the imaging pipeline is enabled, querying the 'backBuffer' will only give you up to date "
"results after you called Screen('Flip') or Screen('AsyncFlipBegin') or Screen('DrawingFinished'), "
"otherwise you may encounter stale results. If our own homegrown frame-sequential stereo mode "
"is in use, querying the 'backLeftBuffer' and 'backRightBuffer' is well defined, after a Screen('Flip') "
"or Screen('DrawingFinished'), but querying the 'frontLeftBuffer' or 'frontRightBuffer' or "
"'frontBuffer' will result in an assignment of left- and right- stereo images that is mostly "
"based on chance, e.g., you may get accidentally swapped left- and rightBuffer images!\n\n"
"\"moviePtr\" is the optional handle to the movie to which the video frame should be "
"added. You can get this handle from the Screen('CreateMovie') function when creating "
"the movie. By default frames are added to the first movie created with Screen('CreateMovie').\n\n"
"\"frameduration\" optionally defines the display duration of the added video frame in "
"units of movie frame intervals. See the help for 'CreateMovie' for further explanation of "
"\"frameduration\".\n\n"
"Movie images are always stored as uint8 images with 8 bits resolution per pixel color component. "
"Images are always stored as four channel RGBA frames.\n\n"
"See Screen('CreateMovie?') for help on movie creation.\n";

static char seeAlsoString[] = "PutImage CopyWindow CreateMovie FinalizeMovie";
	
// This also works as 'AddFrameToMovie', as almost all code is shared with 'GetImage'.
// Only difference is where the fetched pixeldata is sent: To the movie encoder or to
// a matlab/octave matrix.
PsychError SCREENGetImage(void) 
{
	PsychRectType   windowRect, sampleRect;
	int 			nrchannels, invertedY;
	size_t			ix, iy, sampleRectWidth, sampleRectHeight, redReturnIndex, greenReturnIndex, blueReturnIndex, alphaReturnIndex, planeSize;
	int				viewid;
	psych_uint8 	*returnArrayBase, *redPlane;
	float 			*dredPlane;
	double 			*returnArrayBaseDouble;
	PsychWindowRecordType	*windowRecord;
	GLboolean		isDoubleBuffer, isStereo;
	char*           buffername = NULL;
	psych_bool		floatprecision = FALSE;
	GLenum			whichBuffer = 0; 
	int				frameduration = 1;
	int				moviehandle = 0;
	unsigned int	twidth, theight;
	unsigned char*	framepixels;

	// Called as 2nd personality "AddFrameToMovie" ?
	psych_bool isAddMovieFrame = PsychMatch(PsychGetFunctionName(), "AddFrameToMovie");

	// All sub functions should have these two lines
	if (isAddMovieFrame) {
		PsychPushHelp(useString2, synopsisString2, seeAlsoString);
	}
	else {
		PsychPushHelp(useString, synopsisString, seeAlsoString);
	}
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//cap the numbers of inputs and outputs
	PsychErrorExit(PsychCapNumInputArgs(5));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(1));  //The maximum number of outputs
	
	// Get windowRecord for this window:
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);

	// Make sure we don't execute on an onscreen window with pending async flip, as this would interfere
	// by touching the system backbuffer -> Impaired timing of the flip thread and undefined readback
	// of image data due to racing with the ops of the flipperthread on the same drawable.
	//
	// Note: It would be possible to allow drawBuffer readback if the drawBuffer is not multi-sampled
	// or if we can safely multisample-resolve without touching the backbuffer, but checking for all
	// special cases adds ugly complexity and is not really worth the effort, so we don't allow this.
	//
	// If this passes then PsychSetDrawingTarget() below will trigger additional validations to check
	// if execution of 'GetImage' is allowed under the current conditions for offscreen windows and
	// textures:
	if (PsychIsOnscreenWindow(windowRecord) && (windowRecord->flipInfo->asyncstate > 0)) {
		PsychErrorExitMsg(PsychError_user, "Calling this function on an onscreen window with a pending asynchronous flip is not allowed!");
	}

	// Set window as drawingtarget: Even important if this binding is changed later on!
	// We need to make sure all needed transitions are done - esp. in non-imaging mode,
	// so backbuffer is in a useable state:
	PsychSetDrawingTarget(windowRecord);
	
	// Disable shaders:
	PsychSetShader(windowRecord, 0);

	// Soft-Reset drawingtarget. This is important to make sure no FBO's are bound,
	// otherwise the following glGets for GL_DOUBLEBUFFER and GL_STEREO will retrieve
	// wrong results, leading to totally wrong read buffer assignments down the road!!
	PsychSetDrawingTarget((PsychWindowRecordType*) 0x1);

	glGetBooleanv(GL_DOUBLEBUFFER, &isDoubleBuffer);
	glGetBooleanv(GL_STEREO, &isStereo);

    // Force "quad-buffered" stereo mode if our own homegrown implementation is active:
    if (windowRecord->stereomode == kPsychFrameSequentialStereo) isStereo = TRUE;
    
	// Assign read buffer:
	if(PsychIsOnscreenWindow(windowRecord)) {
		// Onscreen window: We read from the front- or front-left buffer by default.
		// This works on single-buffered and double buffered contexts in a consistent fashion:
		
		// Copy in optional override buffer name:
		PsychAllocInCharArg(3, FALSE, &buffername);
		
		// Override buffer name provided?
		if (buffername) {
			// Which one is it?
			
			// "frontBuffer" is always a valid choice:
			if (PsychMatch(buffername, "frontBuffer")) whichBuffer = GL_FRONT;
			// Allow selection of left- or right front stereo buffer in stereo mode:
			if (PsychMatch(buffername, "frontLeftBuffer") && isStereo) whichBuffer = GL_FRONT_LEFT;
			if (PsychMatch(buffername, "frontRightBuffer") && isStereo) whichBuffer = GL_FRONT_RIGHT;
			// Allow selection of backbuffer in double-buffered mode:
			if (PsychMatch(buffername, "backBuffer") && isDoubleBuffer) whichBuffer = GL_BACK;
			// Allow selection of left- or right back stereo buffer in stereo mode:
			if (PsychMatch(buffername, "backLeftBuffer") && isStereo && isDoubleBuffer) whichBuffer = GL_BACK_LEFT;
			if (PsychMatch(buffername, "backRightBuffer") && isStereo && isDoubleBuffer) whichBuffer = GL_BACK_RIGHT;
			// Allow AUX buffer access for debug purposes:
			if (PsychMatch(buffername, "aux0Buffer")) whichBuffer = GL_AUX0;
			if (PsychMatch(buffername, "aux1Buffer")) whichBuffer = GL_AUX1;
			if (PsychMatch(buffername, "aux2Buffer")) whichBuffer = GL_AUX2;
			if (PsychMatch(buffername, "aux3Buffer")) whichBuffer = GL_AUX3;			
		}
		else {
			// Default is frontbuffer:
			whichBuffer = GL_FRONT;
		}
	}
	else {
		// Offscreen window or texture: They only have one buffer, which is the
		// backbuffer in double-buffered mode and the frontbuffer in single buffered mode:
		whichBuffer=(isDoubleBuffer) ? GL_BACK : GL_FRONT;
	}
	
	// Enable this windowRecords framebuffer as current drawingtarget. This should
	// also allow us to "GetImage" from Offscreen windows:
	if ((windowRecord->imagingMode & kPsychNeedFastBackingStore) || (windowRecord->imagingMode & kPsychNeedFastOffscreenWindows)) {
		// Special case: Imaging pipeline active - We need to activate system framebuffer
		// so we really read the content of the framebuffer and not of some FBO:
		if (PsychIsOnscreenWindow(windowRecord)) {
			// It's an onscreen window:
            
            // Homegrown frame-sequential stereo active? Need to remap some stuff:
            if (windowRecord->stereomode == kPsychFrameSequentialStereo) {
                // Back/Front buffers map to backleft/frontleft buffers:
                if (whichBuffer == GL_BACK) whichBuffer = GL_BACK_LEFT;
                if (whichBuffer == GL_FRONT) whichBuffer = GL_FRONT_LEFT;
                
                // Special case: Want to read from stereo front buffer?
                if ((whichBuffer == GL_FRONT_LEFT) || (whichBuffer == GL_FRONT_RIGHT)) {
                    // These don't really exist in our homegrown implementation. Their equivalents are the
                    // regular system front/backbuffers. Due to the bufferswaps happening every video
                    // refresh cycle and the complex logic on when and how to blit finalizedFBOs into
                    // the system buffers and the asynchronous execution of the parallel flipper thread,
                    // we don't know which buffer (GL_BACK or GL_FRONT) corresponds to the leftFront or
                    // rightFront buffer. Let's be stupid and just return the current front buffer for
                    // FRONT_LEFT and the current back buffer for FRONT_RIGHT, but warn user about the
                    // ambiguity:
                    whichBuffer = (whichBuffer == GL_FRONT_LEFT) ? GL_FRONT : GL_BACK;
                    
                    if (PsychPrefStateGet_Verbosity() > 2) {
                        printf("PTB-WARNING: In Screen('GetImage'): You selected retrieval of one of the stereo front buffers, while our homegrown frame-sequential\n");
                        printf("PTB-WARNING: In Screen('GetImage'): stereo display mode is active. This will impair presentation timing and may cause flicker. The\n");
                        printf("PTB-WARNING: In Screen('GetImage'): mapping of 'frontLeftBuffer' and 'frontRightBuffer' to actual stimulus content is very ambiguous\n");
                        printf("PTB-WARNING: In Screen('GetImage'): in this mode. You may therefore end up with the content of the wrong buffer returned! Check results\n");
                        printf("PTB-WARNING: In Screen('GetImage'): carefully! Better read from 'backLeftBuffer' or 'backRightBuffer' for well defined results.\n\n");
                    }
                }
            }
            
            // Homegrown frame-sequential stereo active and backleft or backright buffer requested?
            if (((whichBuffer == GL_BACK_LEFT) || (whichBuffer == GL_BACK_RIGHT)) && (windowRecord->stereomode == kPsychFrameSequentialStereo)) {
                // We can get the equivalent of the backLeft/RightBuffer from the finalizedFBO's in this mode. Get their content:                
				viewid = (whichBuffer == GL_BACK_RIGHT) ? 1 : 0;
				whichBuffer = GL_COLOR_ATTACHMENT0_EXT;
                
                // Bind finalizedFBO as framebuffer to read from:
                glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, windowRecord->fboTable[windowRecord->finalizedFBO[viewid]]->fboid);
                
                // Make sure binding gets released at end of routine:
                viewid = -1;
                
            } // No frame-sequential stereo: Full imaging pipeline active and one of the drawBuffer's requested?
            else if (buffername && (PsychMatch(buffername, "drawBuffer")) && (windowRecord->imagingMode & kPsychNeedFastBackingStore)) {
				// Activate drawBufferFBO:
				PsychSetDrawingTarget(windowRecord);
				whichBuffer = GL_COLOR_ATTACHMENT0_EXT;
				
				// Is the drawBufferFBO multisampled?
				viewid = (((windowRecord->stereomode > 0) && (windowRecord->stereodrawbuffer == 1)) ? 1 : 0);
				if (windowRecord->fboTable[windowRecord->drawBufferFBO[viewid]]->multisample > 0) {
					// It is! We can't read from a multisampled FBO. Need to perform a multisample resolve operation and read
					// from the resolved unisample buffer instead. This is only safe if the unisample buffer is either a dedicated
					// FBO, or - in case its the final system backbuffer etc. - if preflip operations haven't been performed yet.
					// If non dedicated buffer (aka finalizedFBO) and preflip ops have already happened, then the backbuffer contains
					// final content for an upcoming Screen('Flip') and we can't use (and therefore taint) that buffer.
					if ((windowRecord->inputBufferFBO[viewid] == windowRecord->finalizedFBO[viewid]) && (windowRecord->backBufferBackupDone)) {
						// Target for resolve is finalized FBO (probably system backbuffer) and preflip ops have run already. We
						// can't do the resolve op, as this would screw up the backbuffer with the final stimulus:
						printf("PTB-ERROR: Tried to 'GetImage' from a multisampled 'drawBuffer', but can't perform anti-aliasing pass due to\n");
						printf("PTB-ERROR: lack of a dedicated resolve buffer.\n");
						printf("PTB-ERROR: You can get what you wanted by either one of two options:\n");
						printf("PTB-ERROR: Either enable a processing stage in the imaging pipeline, even if you don't need it, e.g., by setting\n");
						printf("PTB-ERROR: the imagingmode argument in the 'OpenWindow' call to kPsychNeedImageProcessing. This will create a\n");
						printf("PTB-ERROR: suitable resolve buffer. Or place the 'GetImage' call before any Screen('DrawingFinished') call, then\n");
						printf("PTB-ERROR: i can (ab-)use the system backbuffer as a temporary resolve buffer.\n\n");
						PsychErrorExitMsg(PsychError_user, "Tried to 'GetImage' from a multi-sampled 'drawBuffer'. Unsupported operation under given conditions.");						
					}
					else {
						// Ok, the inputBufferFBO is a suitable temporary resolve buffer. Perform a multisample resolve blit to it:
						// A simple glBlitFramebufferEXT() call will do the copy & downsample operation:
						glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, windowRecord->fboTable[windowRecord->drawBufferFBO[viewid]]->fboid);
						glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, windowRecord->fboTable[windowRecord->inputBufferFBO[viewid]]->fboid);
						glBlitFramebufferEXT(0, 0, windowRecord->fboTable[windowRecord->inputBufferFBO[viewid]]->width, windowRecord->fboTable[windowRecord->inputBufferFBO[viewid]]->height,
											 0, 0, windowRecord->fboTable[windowRecord->inputBufferFBO[viewid]]->width, windowRecord->fboTable[windowRecord->inputBufferFBO[viewid]]->height,
											 GL_COLOR_BUFFER_BIT, GL_NEAREST);

						// Bind inputBuffer as framebuffer:
						glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, windowRecord->fboTable[windowRecord->inputBufferFBO[viewid]]->fboid);
						viewid = -1;
					}
				}
			}
			else {
				// No: Activate system framebuffer:
				PsychSetDrawingTarget(NULL);
			}
		}
		else {
			// Offscreen window or texture: Select drawing target as usual,
			// but set color attachment as read buffer:
			PsychSetDrawingTarget(windowRecord);
			whichBuffer = GL_COLOR_ATTACHMENT0_EXT;

			// We do not support multisampled readout:
			if (windowRecord->fboTable[windowRecord->drawBufferFBO[0]]->multisample > 0) {
				printf("PTB-ERROR: You tried to Screen('GetImage', ...); from an offscreen window or texture which has multisample anti-aliasing enabled.\n");
				printf("PTB-ERROR: This operation is not supported. You must first use Screen('CopyWindow') to create a non-multisampled copy of the\n");
				printf("PTB-ERROR: texture or offscreen window, then use 'GetImage' on that copy. The copy will be anti-aliased, so you'll get what you\n");
				printf("PTB-ERROR: wanted with a bit more effort. Sorry for the inconvenience, but this is mostly a hardware limitation.\n\n");
				
				PsychErrorExitMsg(PsychError_user, "Tried to 'GetImage' from a multi-sampled texture or offscreen window. Unsupported operation.");
			}
		}
	}
	else {
		// Normal case: No FBO based imaging - Select drawing target as usual:
		PsychSetDrawingTarget(windowRecord);
	}
	
	// Select requested read buffer, after some double-check:
	if (whichBuffer == 0) PsychErrorExitMsg(PsychError_user, "Invalid or unknown 'bufferName' argument provided.");
	glReadBuffer(whichBuffer);

	if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In Screen('GetImage'): GL-Readbuffer whichBuffer = %i\n", whichBuffer);

    if (whichBuffer == GL_COLOR_ATTACHMENT0_EXT) {
        // FBO of texture / offscreen window / onscreen drawBuffer/inputBuffer
        // has size of clientrect -- potentially larger or smaller than backbuffer:
        PsychCopyRect(windowRect, windowRecord->clientrect);
    }
    else {
        // Non-FBO backed texture / offscreen window / onscreen window has size
        // of raw rect (==clientrect for non-onscreen, == backbuffer size for onscreen):
        PsychCopyRect(windowRect, windowRecord->rect);
    }

	// Retrieve optional read rectangle:    
	if(!PsychCopyInRectArg(2, FALSE, sampleRect)) PsychCopyRect(sampleRect, windowRect);
    
	if (IsPsychRectEmpty(sampleRect)) return(PsychError_none);

	// Compute sampling rectangle:
	if ((PsychGetWidthFromRect(sampleRect) >= INT_MAX) || (PsychGetHeightFromRect(sampleRect) >= INT_MAX)) {
		PsychErrorExitMsg(PsychError_user, "Too big 'rect' argument provided. Both width and height of the rect must not exceed 2^31 pixels!");
	}

	sampleRectWidth = (size_t) PsychGetWidthFromRect(sampleRect);
	sampleRectHeight= (size_t) PsychGetHeightFromRect(sampleRect);

	// Regular image fetch to runtime, or adding to a movie?
	if (!isAddMovieFrame) {
		// Regular fetch:

		// Get optional floatprecision flag: We return data with float-precision if
		// this flag is set. By default we return uint8 data:
		PsychCopyInFlagArg(4, FALSE, &floatprecision);
		
		// Get the optional number of channels flag: By default we return 3 channels,
		// the Red, Green, and blue color channel:
		nrchannels = 3;
		PsychCopyInIntegerArg(5, FALSE, &nrchannels);
		if (nrchannels < 1 || nrchannels > 4) PsychErrorExitMsg(PsychError_user, "Number of requested channels 'nrchannels' must be between 1 and 4!");
		
		if (!floatprecision) {
			// Readback of standard 8bpc uint8 pixels:  
			PsychAllocOutUnsignedByteMatArg(1, TRUE, (int) sampleRectHeight, (int) sampleRectWidth, (int) nrchannels, &returnArrayBase);
			redPlane  = (psych_uint8*) PsychMallocTemp((size_t) nrchannels * sampleRectWidth * sampleRectHeight);
			planeSize = sampleRectWidth * sampleRectHeight;

			glPixelStorei(GL_PACK_ALIGNMENT,1);
			invertedY = (int) (windowRect[kPsychBottom] - sampleRect[kPsychBottom]);

			if (nrchannels==1) glReadPixels((int) sampleRect[kPsychLeft], invertedY, (int) sampleRectWidth, (int) sampleRectHeight, GL_RED, GL_UNSIGNED_BYTE, redPlane); 
			if (nrchannels==2) glReadPixels((int) sampleRect[kPsychLeft], invertedY, (int) sampleRectWidth, (int) sampleRectHeight, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, redPlane);
			if (nrchannels==3) glReadPixels((int) sampleRect[kPsychLeft], invertedY, (int) sampleRectWidth, (int) sampleRectHeight, GL_RGB, GL_UNSIGNED_BYTE, redPlane);
			if (nrchannels==4) glReadPixels((int) sampleRect[kPsychLeft], invertedY, (int) sampleRectWidth, (int) sampleRectHeight, GL_RGBA, GL_UNSIGNED_BYTE, redPlane);
			
			//in one pass transpose and flip what we read with glReadPixels before returning.  
			//-glReadPixels insists on filling up memory in sequence by reading the screen row-wise whearas Matlab reads up memory into columns.
			//-the Psychtoolbox screen as setup by gluOrtho puts 0,0 at the top left of the window but glReadPixels always believes that it's at the bottom left.     
			for(ix=0; ix < sampleRectWidth; ix++){
				for(iy=0; iy < sampleRectHeight; iy++){
					// Compute write-indices for returned data:
					redReturnIndex=PsychIndexElementFrom3DArray(sampleRectHeight, sampleRectWidth, nrchannels, iy, ix, 0);
					greenReturnIndex=PsychIndexElementFrom3DArray(sampleRectHeight, sampleRectWidth,  nrchannels, iy, ix, 1);
					blueReturnIndex=PsychIndexElementFrom3DArray(sampleRectHeight, sampleRectWidth,  nrchannels, iy, ix, 2);
					alphaReturnIndex=PsychIndexElementFrom3DArray(sampleRectHeight, sampleRectWidth,  nrchannels, iy, ix, 3);
					
					// Always return RED/LUMINANCE channel:
					returnArrayBase[redReturnIndex] = redPlane[(ix + ((sampleRectHeight-1) - iy ) * sampleRectWidth) * (size_t) nrchannels + 0];  
					// Other channels on demand:
					if (nrchannels>1) returnArrayBase[greenReturnIndex] = redPlane[(ix + ((sampleRectHeight-1) - iy ) * sampleRectWidth) * (size_t) nrchannels + 1];
					if (nrchannels>2) returnArrayBase[blueReturnIndex]  = redPlane[(ix + ((sampleRectHeight-1) - iy ) * sampleRectWidth) * (size_t) nrchannels + 2];
					if (nrchannels>3) returnArrayBase[alphaReturnIndex] = redPlane[(ix + ((sampleRectHeight-1) - iy ) * sampleRectWidth) * (size_t) nrchannels + 3];
				}
			}		
		}
		else {
			// Readback of standard 32bpc float pixels into a double matrix:  
			PsychAllocOutDoubleMatArg(1, TRUE, (int) sampleRectHeight, (int) sampleRectWidth, (int) nrchannels, &returnArrayBaseDouble);
			dredPlane = (float*) PsychMallocTemp((size_t) nrchannels * sizeof(float) * sampleRectWidth * sampleRectHeight);
			planeSize = sampleRectWidth * sampleRectHeight * sizeof(float);

			glPixelStorei(GL_PACK_ALIGNMENT, 1);
			invertedY = (int) (windowRect[kPsychBottom]-sampleRect[kPsychBottom]);

			if (nrchannels==1) glReadPixels((int) sampleRect[kPsychLeft], invertedY, (int) sampleRectWidth, (int) sampleRectHeight, GL_RED, GL_FLOAT, dredPlane); 
			if (nrchannels==2) glReadPixels((int) sampleRect[kPsychLeft], invertedY, (int) sampleRectWidth, (int) sampleRectHeight, GL_LUMINANCE_ALPHA, GL_FLOAT, dredPlane);
			if (nrchannels==3) glReadPixels((int) sampleRect[kPsychLeft], invertedY, (int) sampleRectWidth, (int) sampleRectHeight, GL_RGB, GL_FLOAT, dredPlane);
			if (nrchannels==4) glReadPixels((int) sampleRect[kPsychLeft], invertedY, (int) sampleRectWidth, (int) sampleRectHeight, GL_RGBA, GL_FLOAT, dredPlane);
			
			//in one pass transpose and flip what we read with glReadPixels before returning.  
			//-glReadPixels insists on filling up memory in sequence by reading the screen row-wise whearas Matlab reads up memory into columns.
			//-the Psychtoolbox screen as setup by gluOrtho puts 0,0 at the top left of the window but glReadPixels always believes that it's at the bottom left.     
			for(ix=0; ix < sampleRectWidth; ix++){
				for(iy=0; iy < sampleRectHeight; iy++){
					// Compute write-indices for returned data:
					redReturnIndex=PsychIndexElementFrom3DArray(sampleRectHeight, sampleRectWidth, nrchannels, iy, ix, 0);
					greenReturnIndex=PsychIndexElementFrom3DArray(sampleRectHeight, sampleRectWidth,  nrchannels, iy, ix, 1);
					blueReturnIndex=PsychIndexElementFrom3DArray(sampleRectHeight, sampleRectWidth,  nrchannels, iy, ix, 2);
					alphaReturnIndex=PsychIndexElementFrom3DArray(sampleRectHeight, sampleRectWidth,  nrchannels, iy, ix, 3);
					
					// Always return RED/LUMINANCE channel:
					returnArrayBaseDouble[redReturnIndex] = dredPlane[(ix + ((sampleRectHeight-1) - iy ) * sampleRectWidth) * (size_t) nrchannels + 0];  
					// Other channels on demand:
					if (nrchannels>1) returnArrayBaseDouble[greenReturnIndex] = dredPlane[(ix + ((sampleRectHeight-1) - iy ) * sampleRectWidth) * (size_t) nrchannels + 1];
					if (nrchannels>2) returnArrayBaseDouble[blueReturnIndex]  = dredPlane[(ix + ((sampleRectHeight-1) - iy ) * sampleRectWidth) * (size_t) nrchannels + 2];
					if (nrchannels>3) returnArrayBaseDouble[alphaReturnIndex] = dredPlane[(ix + ((sampleRectHeight-1) - iy ) * sampleRectWidth) * (size_t) nrchannels + 3];
				}
			}		
		}
	}
	
	if (isAddMovieFrame) {
		// Adding of image to a movie requested:
		
		// Get optional moviehandle:
		moviehandle = 0;
		PsychCopyInIntegerArg(4, FALSE, &moviehandle);
		if (moviehandle < 0) PsychErrorExitMsg(PsychError_user, "Provided 'moviehandle' is negative. Must be greater or equal to zero!");
		
		// Get optional frameduration:
		frameduration = 1;
		PsychCopyInIntegerArg(5, FALSE, &frameduration);
		if (frameduration < 1) PsychErrorExitMsg(PsychError_user, "Number of requested framedurations 'frameduration' is negative. Must be greater than zero!");
		
		framepixels = PsychGetVideoFrameForMoviePtr(moviehandle, &twidth, &theight);
		if (framepixels) {
			glPixelStorei(GL_PACK_ALIGNMENT,1);
			invertedY = (int) (windowRect[kPsychBottom] - sampleRect[kPsychBottom]);
			
			glReadPixels((int) sampleRect[kPsychLeft], invertedY, twidth, theight, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8, framepixels);
			if (PsychAddVideoFrameToMovie(moviehandle, frameduration, TRUE) != 0) {
				printf("See http://developer.apple.com/documentation/QuickTime/APIREF/ErrorCodes.htm#//apple_ref/doc/constant_group/Error_Codes.\n\n");
				PsychErrorExitMsg(PsychError_user, "AddFrameToMovie failed with error above!");
			}
		}
		else {
			PsychErrorExitMsg(PsychError_user, "Invalid 'moviePtr' provided. Doesn't correspond to a movie open for recording!");
		}
	}
	
	if (viewid == -1) {
		// Need to reset framebuffer binding to get rid of the inputBufferFBO which is bound due to
		// multisample resolve ops, or of other special FBO bindings --> Activate system framebuffer:
		PsychSetDrawingTarget(NULL);		
	}

	return(PsychError_none);
}
