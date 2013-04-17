/*
	SCREENSelectStereoDrawBuffer.c		
  
	AUTHORS:

		mario.kleiner at tuebingen.mpg.de 		mk 
  
	PLATFORMS:	
	
		All.
    

	HISTORY:

		04/03/05	mk		Created.
		10/30/05        mk              Synopsis string changed.

	DESCRIPTION:
  
		Selects the target buffer for drawing commands on a stereoscopic display:
                All drawing commands after this command will apply to the selected buffer.
  
	TO DO:  

*/


#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "currentbuffer = Screen('SelectStereoDrawBuffer', windowPtr [, bufferid] [, param1]);";
static char synopsisString[] = 
		"Select the target buffer for drawing commands in stereo display mode and or return "
		"the current/old buffer selection in 'currentbuffer'.\n"
		"This function only applies to stereo mode, it does nothing in mono mode. "
		"\"windowPtr\" is the pointer to the onscreen stereo window. "
		"\"bufferid\" (optional) is either == 0 for selecting the left-eye buffer or == 1 for "
		"selecting the right-eye buffer. You need to call this command after each "
		"Screen('Flip') command or after drawing to an offscreen window again in order "
		"to reestablish your selection of draw buffer, otherwise the results of drawing "
		"operations will be undefined and most probably not what you want.\n"
		"\"param1\" (optional) a parameter whose meaning depends on the active stereo mode:\n"
		"In stereoModes 1 and 11 (frame sequential stereo) it allows to select if stimulus onset should "
		"happen in an even video refresh interval (value 0) or in an odd interval (value 1). "
		"Even and odd intervals correspond to either left- or right-eye view, so this allows "
		"to choose if you want to have stimulus onset on left eye or right eye.\n"
		"If you want to use a stereo display mode, we recommend enabling the imaging "
		"pipeline as well. The imaging pipeline "
		"will only work with recent hardware, but it allows for more reliable stereo "
		"operation especially when using low-level OpenGL for real 3D drawing. It also "
		"allows to parameterize aspects of stereo presentation, e.g., gains and other "
		"settings of anaglyph stereo.";

static char seeAlsoString[] = "OpenWindow Flip";	 

PsychError SCREENSelectStereoDrawBuffer(void) 
{
	PsychWindowRecordType *windowRecord;
	int bufferid, param1;
	
	//all subfunctions should have these two lines.  
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	PsychErrorExit(PsychCapNumInputArgs(3));     //The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1)); //The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(1));    //The maximum number of outputs
        
	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
        
	if(!PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_user, "Tried to select stereo draw buffer on something else than a onscreen window.");
            	
	// Return (optionally) current buffer id:
	PsychCopyOutDoubleArg(1, FALSE, (double) windowRecord->stereodrawbuffer);
	
	// Get the optional buffer id (0==left, 1==right). If none provided then we're done and can return immediately:
	if (!PsychCopyInIntegerArg(2, FALSE, &bufferid)) return(PsychError_none);
	
	// Bufferid provided...
	if (bufferid<0 || bufferid>1) PsychErrorExitMsg(PsychError_user, "Invalid bufferid provided: Must be 0 for left-eye or 1 for right-eye buffer.");
	
	// Trying to select other than left buffer on mono-window?
	if(windowRecord->windowType!=kPsychDoubleBufferOnscreen || windowRecord->stereomode == kPsychMonoscopic) {
		// Yes. Reset to left buffer, which is used for mono mode:
		bufferid = 0;
	}
	
	// Store assignment in windowRecord:
	windowRecord->stereodrawbuffer = bufferid;
	
	// Copy in optional param1 argument:
	param1 = -10000;
	PsychCopyInIntegerArg(3, FALSE, &param1);
	if (param1!=-10000) {
		if ((windowRecord->stereomode == kPsychOpenGLStereo) || (windowRecord->stereomode == kPsychFrameSequentialStereo)) {
			if (param1<-1 || param1>1) PsychErrorExitMsg(PsychError_user, "Invalid fieldid provided: Must be -1 for don't care, 0 for even field or 1 for odd field.");
			// Valid: Assign it.
			windowRecord->targetFlipFieldType = param1;
		}
	}

	// If the imaging pipeline is active, then we're done.
	if (windowRecord->imagingMode & kPsychNeedFastBackingStore) {
		// Enable this windowRecords framebuffer as current drawingtarget. In imagingmode this will also
		// select the proper backing FBO:
		if (PsychGetDrawingTarget() == windowRecord) PsychSetDrawingTarget(NULL);
		PsychSetDrawingTarget(windowRecord);
		
		// Done for imaging mode:
		return(PsychError_none);
	}

	// The following code is only used for non-imaging mode operations:	
	PsychSetDrawingTarget(windowRecord);
	PsychSwitchFixedFunctionStereoDrawbuffer(windowRecord);
	return(PsychError_none);
}

/*	This routine performs the setup work for switching between left- and right-eye
	drawbuffers in stereo mode when using the fixed function pipeline, i.e., imagingpipeline
	is disabled.
	
	It is called by SCREENSelectStereoDrawBuffer() directly for setup of PTB's internal 2D
	drawing context. It's also called by SCREENBeginOpenGL() when performing a transition to
	the userspace 3D OpenGL drawing context, so both 2D and 3D contexts stay in sync for
	stereo drawing.
	
	Setup of drawbuffers in imaging mode is much simpler, automatically performed in
	PsychSetDrawingTarget() whenever needed by a simple switch of the drawBuffer FBO's.
	
 */
void PsychSwitchFixedFunctionStereoDrawbuffer(PsychWindowRecordType *windowRecord)
{
	int screenwidth, screenheight;
	int bufferid=windowRecord->stereodrawbuffer;

	if(windowRecord->windowType!=kPsychDoubleBufferOnscreen || windowRecord->stereomode == kPsychMonoscopic) {
		// Trying to select the draw target buffer on a non-stereo window: We just reset it to monoscopic default.
		glDrawBuffer(GL_BACK);
		return;
	}
	
	// OpenGL native stereo?
	if (windowRecord->stereomode==kPsychOpenGLStereo) {
		// OpenGL native stereo via separate back-buffers: Select target draw buffer:
		switch(bufferid) {
			case 0:
				glDrawBuffer(GL_BACK_LEFT);
				break;
			case 1:
				glDrawBuffer(GL_BACK_RIGHT);
				break;
		}		
		return;
	}
	
	// Vertical compression stereo?
	if (windowRecord->stereomode==kPsychCompressedTLBRStereo || windowRecord->stereomode==kPsychCompressedTRBLStereo) {
		PsychSwitchCompressedStereoDrawBuffer(windowRecord, bufferid);
		return;
	}
	
	// "Free fusion" stereo? Simply place views side-by-side, downscaled by a factor of 2 in horizontal dimension...
	if (windowRecord->stereomode==kPsychFreeFusionStereo || windowRecord->stereomode==kPsychFreeCrossFusionStereo) {
	    // Switch between drawing into left- and right-half of the single framebuffer:
	    screenwidth=(int) PsychGetWidthFromRect(windowRecord->rect);
		screenheight=(int) PsychGetHeightFromRect(windowRecord->rect);
		
		// Cross fusion instead of fusion requested? Switch left-right if so:
		if (windowRecord->stereomode==kPsychFreeCrossFusionStereo) bufferid = 1 - bufferid;
		
		// Change of projection matrix is only applied if we are not executed in a userspace
		// GL context. Us changing the matrix here in conjunction with MOGL would defeat
		// 3D perpective projection code and thereby interfere with real 3D perspective rendering:
		if (!PsychIsUserspaceRendering()) {
			// Setup projection matrix for ortho-projection of full window height, but only
			// half window width:
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
            if (!PsychIsGLES(windowRecord)) {
                gluOrtho2D(0, screenwidth/2, windowRecord->rect[kPsychBottom], windowRecord->rect[kPsychTop]);
            }
            else {
                glOrthofOES((GLfloat) 0, (GLfloat) screenwidth/2, (GLfloat) windowRecord->rect[kPsychBottom], (GLfloat) windowRecord->rect[kPsychTop], (GLfloat) -1, (GLfloat) 1);
            }

			// Switch back to modelview matrix, but leave it unaltered:
			glMatrixMode(GL_MODELVIEW);
		}
	    // When using this stereo modes, we are most likely running on a
	    // dual display setup with desktop set to "horizontal spanning mode". In this mode,
	    // we get a virtual onscreen window that is at least twice as wide as its height and
	    // the left half of the window is displayed on the left monitor, the right half of the
	    // window is displayed on the right monitor. To make good use of the space, we only
	    // scale the viewport horizontally to half the window width, but keep it at full height:
	    // All Screen subfunctions that report the size of the onscreen window to Matlab/Octave/...
	    // will report it to be only half of its real width, so experiment code can adapt to it. 
	    switch(bufferid) {
			case 0:
				// Place viewport in the left half of screen:
				glViewport(0, 0, screenwidth/2, screenheight);
				glScissor(0, 0, screenwidth/2, screenheight);
				break;
				
			case 1:
				// Place viewport in the right half of screen:
				glViewport(screenwidth/2, 0, screenwidth/2, screenheight);
				glScissor(screenwidth/2, 0, screenwidth/2, screenheight);
				break;
		}
		return;
	}
	
	// And now for the Anaglyph stereo modes, were the left-eye vs. right-eye images are encoded in
	// the separate color-channels of the same display. We do this via the OpenGL writemask, which
	// allows to selectively enable/disable write operations to the different color channels:
	// The alpha channel is always enabled, the red,gree,blue channels are depending on mode and
	// bufferid conditionally enabled/disabled:	
	switch (windowRecord->stereomode) {
		case kPsychAnaglyphRGStereo:
			glColorMask(bufferid==0, bufferid==1, FALSE, TRUE);
			windowRecord->colorMask[0] = (bufferid==0) ? GL_TRUE : GL_FALSE;
			windowRecord->colorMask[1] = (bufferid==1) ? GL_TRUE : GL_FALSE;
			windowRecord->colorMask[2] = GL_FALSE;
			windowRecord->colorMask[3] = GL_TRUE;					
            break;
			
		case kPsychAnaglyphGRStereo:
			glColorMask(bufferid==1, bufferid==0, FALSE, TRUE);
			windowRecord->colorMask[0] = (bufferid==1) ? GL_TRUE : GL_FALSE;
			windowRecord->colorMask[1] = (bufferid==0) ? GL_TRUE : GL_FALSE;
			windowRecord->colorMask[2] = GL_FALSE;
			windowRecord->colorMask[3] = GL_TRUE;
			
            break;
			
		case kPsychAnaglyphRBStereo:
			glColorMask(bufferid==0, FALSE, bufferid==1, TRUE);
			windowRecord->colorMask[0] = (bufferid==0) ? GL_TRUE : GL_FALSE;
			windowRecord->colorMask[2] = (bufferid==1) ? GL_TRUE : GL_FALSE;
			windowRecord->colorMask[1] = GL_FALSE;			
			windowRecord->colorMask[3] = GL_TRUE;
            break;
			
		case kPsychAnaglyphBRStereo:
			glColorMask(bufferid==1, FALSE, bufferid==0, TRUE);
			windowRecord->colorMask[0] = (bufferid==1) ? GL_TRUE : GL_FALSE;
			windowRecord->colorMask[2] = (bufferid==0) ? GL_TRUE : GL_FALSE;
			windowRecord->colorMask[1] = GL_FALSE;			
			windowRecord->colorMask[3] = GL_TRUE;
            break;
	}
	return;
}

/* 
 * PsychSwitchCompressedStereoDrawBuffer  -- Handle buffer transitions
 *
 * If in stereomode for vertically compressed stereo display, switch between the
 * two drawing buffers via some AUX buffer copy magic...
 *
 * Vertical compressed stereo works as follows:
 *
 * The left-eye view is drawn at full resolution into the backbuffer as if it would
 * be a regular monoscopic stimulus on a full-res mono-display. When user wants to
 * switch to right-eye drawing, this routine makes a full scale backup copy of the
 * framebuffer into AUX buffer 0 --> AUX0 contains the full res left-eye image.
 * Now the right eye image is drawn into the backbuffer. When switch back to left-eye
 * view occurs, the backbuffer is copied into AUX1 buffer.
 *
 * All possible transitions between the two eye buffers are handled by copying data
 * from backbuffer to AUX0/1 and back. The net result is that at Screen('DrawingFinished')
 * time or Screen('Flip') time, AUX0 contains the left-eye image, AUX1 contains the right-eye
 * image. Additional logic in PsychPreFlipOperations() makes sure this is always the case,
 * regardless what the user does. The content of both AUX buffers is now copied back into
 * the backbuffer, AUX0 into top half, AUX1 into bottom half of screen (or vice versa, depending
 * on mode). Both buffers are vertically scaled to half resolution, so both fit onto screen
 * simultaneously. This is done by PsychComposeCompressedStereoBuffer().
 *
 */
int PsychSwitchCompressedStereoDrawBuffer(PsychWindowRecordType *windowRecord, int newbuffer)
{
    // Query screen dimension:
    int screenwidth=(int) PsychGetWidthFromRect(windowRecord->rect);
    int screenheight=(int) PsychGetHeightFromRect(windowRecord->rect);

    // Query currently active buffer:
    int oldbuffer = windowRecord->stereodrawbuffer;

    // Transition necessary?
    if (oldbuffer != newbuffer) {
        // Real transition requested...

		// Disable shader:
		PsychSetShader(windowRecord, 0);
	
        glDisable(GL_BLEND);
        
        // Set transform matrix to well-defined state:
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        //glLoadIdentity();
		glPixelZoom(1,1);

        // Do we have a drawbuffer active at the moment?
        if (oldbuffer!=2) {
            // Stereo drawbuffer is active: Content of our backbuffer represents
            // the content of this drawbuffer. We need to save a backup copy in
            // corresponding AUX-Buffer.
            glReadBuffer(GL_BACK);
            glDrawBuffer((oldbuffer==0) ? GL_AUX0 : GL_AUX1);
            glRasterPos2i(0, screenheight);
            glCopyPixels(0, 0, screenwidth, screenheight, GL_COLOR);            
            // Mark this AUX buffer as "dirty" it contains real image content
            // instead of just background color.
            windowRecord->auxbuffer_dirty[oldbuffer] = TRUE;
        }
        
        glDrawBuffer(GL_BACK);
        
        // Ok, old content backed-up if necessary. Now we switch to new buffer...
        if (newbuffer!=2) {
            // Switch to real buffer requested. Check if corresponding AUX buffer
            // is dirty, aka contains real content instead of just background color.
            if (windowRecord->auxbuffer_dirty[newbuffer]) {
                // AUX buffer contains real content. Copy its content back into
                // backbuffer:
                glReadBuffer((newbuffer==0) ? GL_AUX0 : GL_AUX1);
                glRasterPos2i(0, screenheight);
                glCopyPixels(0, 0, screenwidth, screenheight, GL_COLOR);
                glReadBuffer(GL_BACK);
            }
            else {
                // AUX buffer is clean, aka just contains background color. We do
                // a clear screen as its faster than an AUX->BACK copy.
                PsychGLClear(windowRecord);
            }
        }
        
        glEnable(GL_BLEND);

        // Restore transform matrix:
        glPopMatrix();
        
        // Store new state:
        windowRecord->stereodrawbuffer = newbuffer;
        
        // Transition finished. Backbuffers content represents new draw buffers content.
    }
    
    // Done. Return id of previous buffer (0=left, 1=right, 2=none).
    return(oldbuffer);
}

/*
 * PsychComposeCompressedStereoBuffer - Final compositing for compressed stereo.
 *
 * This routine copies both AUX buffers (0 and 1) back into the backbuffer, each of them
 * vertically scaled/compressed by a factor of 2. Its called by PsychPreFlipOperations().
 */
void PsychComposeCompressedStereoBuffer(PsychWindowRecordType *windowRecord)
{
  /*
    if (FALSE) {
        // Upload, setup and enable Anaglyph stereo fragment shader:
        const float redgain=1.0, greengain=0.7, bluegain=0.0;
        
        // This is the shader source code:
        const char anaglyphshader[] =
        "!!ARBfp1.0 "
        "PARAM ColorToGrayWeights = { 0.3, 0.59, 0.11, 1.0 }; "
        "TEMP luminance; "
        "TEMP incolor;"
        //"MOV incolor, fragment.color;"
        "TEX incolor, fragment.texcoord[0], texture[0], RECT;"
        "DP3 luminance, incolor, ColorToGrayWeights; "
        "MUL result.color.rgb, luminance, program.env[0]; "
        "MOV result.color.a, ColorToGrayWeights.a; "
        //"MOV result.color.a, fragment.color.a; "
        "END";
        
        // Upload and compile shader:
        PsychTestForGLErrors();
        glProgramStringARB(GL_FRAGMENT_PROGRAM_ARB, GL_PROGRAM_FORMAT_ASCII_ARB, strlen(anaglyphshader), anaglyphshader);
        PsychTestForGLErrors();
        // Setup the rgb-gains as global parameters for the shader:
        glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB, 0, redgain, greengain, bluegain, 0.0);
        // Enable the shader:
        glEnable(GL_FRAGMENT_PROGRAM_ARB);
    }
  */
    // Query screen dimension:
    int screenwidth=(int) PsychGetWidthFromRect(windowRecord->rect);
    int screenheight=(int) PsychGetHeightFromRect(windowRecord->rect);

    // When entering this routine, the modelview matrix is already set to identity and
    // the proper OpenGL context is active.

	// Disable shader:
	PsychSetShader(windowRecord, 0);

    // Set up zoom for vertical compression:
    glPixelZoom(1, 0.5f);
    glDrawBuffer(GL_BACK);
    glDisable(GL_BLEND);
    
    // Draw left view aka AUX0:
    glReadBuffer(GL_AUX0);
    glRasterPos2i(0, (windowRecord->stereomode==kPsychCompressedTLBRStereo) ? screenheight/2 : screenheight);
    glCopyPixels(0, 0, screenwidth, screenheight, GL_COLOR);

    // Draw right view aka AUX1:
    glReadBuffer(GL_AUX1);
    glRasterPos2i(0, (windowRecord->stereomode==kPsychCompressedTLBRStereo) ? screenheight : screenheight/2);
    glCopyPixels(0, 0, screenwidth, screenheight, GL_COLOR);
    
    // Restore settings:
    glReadBuffer(GL_BACK);
    glPixelZoom(1,1);
    glEnable(GL_BLEND);

    // Unconditionally disable fragment shaders:
    // glDisable(GL_FRAGMENT_PROGRAM_ARB);

    // Done.
    return;
}
