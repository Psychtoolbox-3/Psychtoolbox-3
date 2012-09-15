/*
	PsychToolbox3/Source/Common/Screen/SCREENCopyWindow.c		

	AUTHORS:
	
	Allen.Ingling@nyu.edu				awi 
	Mario.Kleiner@tuebingen.mpg.de		mk

	PLATFORMS:	
		
	All.  
    
	HISTORY:
	
	02/19/03	awi		Created.
	03/11/04	awi		Modified for textures
	1/11/05		awi		Cosmetic
	1/14/05		awi		added glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE) at mk's suggestion;
	1/25/05		awi		Relocated glTexEnvf below glBindTexture.  Fix provided by mk.  
	1/22/05		 mk		Completely rewritten for the new OffscreenWindow implementation.
 
	TO DO:

*/

#include "Screen.h"

static char useString[] =  "Screen('CopyWindow',srcWindowPtr,dstWindowPtr,[srcRect],[dstRect],[copyMode])";
static char synopsisString[] =  "Copy images, quickly, between two windows (on- or off- screen). "
                                "srcRect and dstRect are set to the size of windows srcWindowPtr and dstWindowPtr "
                                "by default. [copyMode] is accepted as input but currently ignored. "
                                "CopyWindow is mostly here for compatibility to PTB-2. If you want to "
                                "copy images really quickly, use the 'MakeTexture' and 'DrawTexture' commands."
                                " They also allow for rotated drawing and advanced blending operations. "
								"The current CopyWindow implementation has a couple of restrictions: "
								"One can't copy from an offscreen window into the -same- offscreen window. "
								"One can't copy from an onscreen window into a -different- onscreen window. "
								"Sizes of sourceRect and targetRect need to match for Onscreen->Offscreen copy. ";

static char seeAlsoString[] = "PutImage, GetImage, OpenOffscreenWindow, MakeTexture, DrawTexture";

PsychError SCREENCopyWindow(void) 
{
	PsychRectType			sourceRect, targetRect;
	PsychWindowRecordType	*sourceWin, *targetWin;
	GLuint					srcFBO, dstFBO;
	GLenum					glerr;
	
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

	//cap the number of inputs
	PsychErrorExit(PsychCapNumInputArgs(5));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs
        
	//get parameters for the source window:
	PsychAllocInWindowRecordArg(1, TRUE, &sourceWin);
	PsychCopyRect(sourceRect, sourceWin->clientrect);
	
	PsychCopyInRectArg(3, FALSE, sourceRect);
	if (IsPsychRectEmpty(sourceRect)) return(PsychError_none);

	//get paramters for the target window:
	PsychAllocInWindowRecordArg(2, TRUE, &targetWin);

	// By default, the targetRect is equal to the sourceRect, but centered in
	// the target window.
	PsychCopyRect(targetRect, targetWin->clientrect);

	PsychCopyInRectArg(4, FALSE, targetRect);
	if (IsPsychRectEmpty(targetRect)) return(PsychError_none);

	if (0) {
		printf("SourceRect: %f %f %f %f  ---> TargetRect: %f %f %f %f\n", sourceRect[0], sourceRect[1],
             sourceRect[2], sourceRect[3], targetRect[0], targetRect[1],targetRect[2],targetRect[3]);
	}

	// Validate rectangles:
	if (!ValidatePsychRect(sourceRect) || sourceRect[kPsychLeft]<sourceWin->rect[kPsychLeft] ||
		sourceRect[kPsychTop]<sourceWin->rect[kPsychTop] || sourceRect[kPsychRight]>sourceWin->rect[kPsychRight] ||
		sourceRect[kPsychBottom]>sourceWin->rect[kPsychBottom]) {
		PsychErrorExitMsg(PsychError_user, "Invalid source rectangle specified - (Partially) outside of source window.");
	}
	
	if (!ValidatePsychRect(targetRect) || targetRect[kPsychLeft]<targetWin->rect[kPsychLeft] ||
		targetRect[kPsychTop]<targetWin->rect[kPsychTop] || targetRect[kPsychRight]>targetWin->rect[kPsychRight] ||
		targetRect[kPsychBottom]>targetWin->rect[kPsychBottom]) {
		PsychErrorExitMsg(PsychError_user, "Invalid target rectangle specified - (Partially) outside of target window.");
	}
	
	if (!(PsychPrefStateGet_ConserveVRAM() & kPsychAvoidCPUGPUSync)) PsychTestForGLErrors();

	// Does this GL implementation support the EXT_framebuffer_blit extension for fast blitting between
	// framebuffers? And is the imaging pipeline active? And is the kPsychAvoidFramebufferBlitIfPossible not set?
	if ((sourceWin->gfxcaps & kPsychGfxCapFBOBlit) && (targetWin->gfxcaps & kPsychGfxCapFBOBlit) &&
		(sourceWin->imagingMode > 0) && (targetWin->imagingMode > 0) && !(PsychPrefStateGet_ConserveVRAM() & kPsychAvoidFramebufferBlitIfPossible)) {
		// Yes :-) -- This simplifies the CopyWindow implementation to a simple framebuffer blit op,
		// regardless what the source- or destination is:
		
		// Set each windows framebuffer as a drawingtarget once: This is just to make sure all of them
		// have proper FBO's attached:
		PsychSetDrawingTarget(sourceWin);
		PsychSetDrawingTarget(targetWin);
		
		// Soft-Reset drawing target - Detach anything bound or setup:
		PsychSetDrawingTarget((PsychWindowRecordType*) 0x1);
		
		// Find source framebuffer:
		if (sourceWin->fboCount == 0) {
			// No FBO's attached to sourceWin: Must be a system framebuffer, e.g., if imagingMode == kPsychNeedFastOffscreenWindows and
			// this is the onscreen window system framebuffer. Assign system framebuffer handle:
			srcFBO = 0;
		}
		else {
			// FBO's attached: Want drawBufferFBO 0 or 1 - 1 for right eye view in stereomode, 0 for
			// everything else: left eye view, monoscopic buffer, offscreen windows / textures:
			if ((sourceWin->stereomode > 0) && (sourceWin->stereodrawbuffer == 1)) {
				// We are in stereo mode and want to access the right-eye channel. Bind FBO-1
				srcFBO = sourceWin->fboTable[sourceWin->drawBufferFBO[1]]->fboid;
			}
			else {
				srcFBO = sourceWin->fboTable[sourceWin->drawBufferFBO[0]]->fboid;
			}
		}

		// Find target framebuffer:
		if (targetWin->fboCount == 0) {
			// No FBO's attached to targetWin: Must be a system framebuffer, e.g., if imagingMode == kPsychNeedFastOffscreenWindows and
			// this is the onscreen window system framebuffer. Assign system framebuffer handle:
			dstFBO = 0;
		}
		else {
			// FBO's attached: Want drawBufferFBO 0 or 1 - 1 for right eye view in stereomode, 0 for
			// everything else: left eye view, monoscopic buffer, offscreen windows / textures:
			if ((targetWin->stereomode > 0) && (targetWin->stereodrawbuffer == 1)) {
				// We are in stereo mode and want to access the right-eye channel. Bind FBO-1
				dstFBO = targetWin->fboTable[targetWin->drawBufferFBO[1]]->fboid;
			}
			else {
				dstFBO = targetWin->fboTable[targetWin->drawBufferFBO[0]]->fboid;
			}
		}

		// Bind read- / write- framebuffers:
		glBindFramebufferEXT(GL_READ_FRAMEBUFFER_EXT, srcFBO);
		glBindFramebufferEXT(GL_DRAW_FRAMEBUFFER_EXT, dstFBO);

		// Perform blit-operation: If blitting from a multisampled FBO into a non-multisampled one, this will also perform the
		// multisample resolve operation:
		glBlitFramebufferEXT((GLint) sourceRect[kPsychLeft], (GLint) (PsychGetHeightFromRect(sourceWin->rect) - sourceRect[kPsychBottom]), (GLint) sourceRect[kPsychRight], (GLint) (PsychGetHeightFromRect(sourceWin->rect) - sourceRect[kPsychTop]),
							 (GLint) targetRect[kPsychLeft], (GLint) (PsychGetHeightFromRect(targetWin->rect) - targetRect[kPsychBottom]), (GLint) targetRect[kPsychRight], (GLint) (PsychGetHeightFromRect(targetWin->rect) - targetRect[kPsychTop]),
							 GL_COLOR_BUFFER_BIT, GL_NEAREST);

		if (PsychPrefStateGet_Verbosity() > 5) {
			printf("FBB-SRC: X0 = %f Y0 = %f X1 = %f Y1 = %f \n", sourceRect[kPsychLeft], PsychGetHeightFromRect(sourceWin->rect) - sourceRect[kPsychBottom], sourceRect[kPsychRight], PsychGetHeightFromRect(sourceWin->rect) - sourceRect[kPsychTop]);
			printf("FBB-DST: X0 = %f Y0 = %f X1 = %f Y1 = %f \n", targetRect[kPsychLeft], PsychGetHeightFromRect(targetWin->rect) - targetRect[kPsychBottom], targetRect[kPsychRight], PsychGetHeightFromRect(targetWin->rect) - targetRect[kPsychTop]);
		}
		
		if (!(PsychPrefStateGet_ConserveVRAM() & kPsychAvoidCPUGPUSync)) {
			if ((glerr = glGetError())!=GL_NO_ERROR) {			

				// Reset framebuffer binding to something safe - The system framebuffer:
				glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
			
				if((glerr == GL_INVALID_OPERATION) && (PsychGetWidthFromRect(sourceRect) != PsychGetWidthFromRect(targetRect) ||
													   PsychGetHeightFromRect(sourceRect) != PsychGetHeightFromRect(targetRect))) {
					// Non-matching sizes. Make sure we do not operate on multisampled stuff
					PsychErrorExitMsg(PsychError_user, "CopyWindow failed: Most likely cause: You tried to copy a multi-sampled window into a non-multisampled window, but there is a size mismatch of sourceRect and targetRect. Matching size is required for such copies.");
				}
				else {
					if (glerr == GL_INVALID_OPERATION) {
						PsychErrorExitMsg(PsychError_user, "CopyWindow failed: Most likely cause: You tried to copy from a multi-sampled window into another multisampled window, but there is a mismatch between the multiSample levels of both. Identical multiSample values are required for such copies.");
					}
					else {
						printf("CopyWindow failed for unresolved reason: OpenGL says after call to glBlitFramebufferEXT(): %s\n", gluErrorString(glerr));
						PsychErrorExitMsg(PsychError_user, "CopyWindow failed for unresolved reason.");
					}
				}
			}
		}
		
		// Reset framebuffer binding to something safe - The system framebuffer:
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		
		// Just to make sure we catch invalid values:
		if (!(PsychPrefStateGet_ConserveVRAM() & kPsychAvoidCPUGPUSync)) PsychTestForGLErrors();
		
		// Done.
		return(PsychError_none);
	}

	// We have four possible combinations for copy ops:
	// Onscreen -> Onscreen
	// Onscreen -> Texture
	// Texture  -> Texture
	// Texture  -> Onscreen
        
	// Texture -> something copy? (Offscreen to Offscreen or Offscreen to Onscreen)
	// This should work for both, copies from a texture/offscreen window to a different texture/offscreen window/onscreen window,
	// and for copies of a subregion of a texture/offscreen window into a non-overlapping subregion of the texture/offscreen window
	// itself:
	if (sourceWin->windowType == kPsychTexture) {
		// Bind targetWin (texture or onscreen windows framebuffer) as
		// drawing target and just blit texture into it. Binding is done implicitely
		
		if ((sourceWin == targetWin) && (targetWin->imagingMode > 0)) {
			// Copy of a subregion of an offscreen window into itself while imaging pipe active, ie. FBO storage: This is actually the same
			// as on onscreen -> onscreen copy, just with the targetWin FBO bound.
			
			// Set target windows framebuffer as drawing target:
			PsychSetDrawingTarget(targetWin);

			// Disable alpha-blending:
			glDisable(GL_BLEND);
			
			// Disable any shading during copy-op:
			PsychSetShader(targetWin, 0);

			// Start position for pixel write is:
			glRasterPos2f((float) targetRect[kPsychLeft], (float) targetRect[kPsychBottom]);
			
			// Zoom factor if rectangle sizes don't match:
			glPixelZoom((float) (PsychGetWidthFromRect(targetRect) / PsychGetWidthFromRect(sourceRect)), (float) (PsychGetHeightFromRect(targetRect) / PsychGetHeightFromRect(sourceRect)));
			
			// Perform pixel copy operation:
			glCopyPixels((int) sourceRect[kPsychLeft], (int) (PsychGetHeightFromRect(sourceWin->rect) - sourceRect[kPsychBottom]), (int) PsychGetWidthFromRect(sourceRect), (int) PsychGetHeightFromRect(sourceRect), GL_COLOR);
			
			// That's it.
			glPixelZoom(1,1);
			
			// Flush drawing commands and wait for render-completion in single-buffer mode:
			PsychFlushGL(targetWin);				
		}
		else {
				// Sourcewin != Targetwin and/or imaging pipe (FBO storage) not used. We blit the
				// backing texture into itself, aka into its representation inside the system
				// backbuffer. The blit routine will setup proper bindings:

				// Disable alpha-blending:
				glDisable(GL_BLEND);

				// We use filterMode == 1 aka Bilinear filtering, so we get nice texture copies
				// if size of sourceRect and targetRect don't match and some scaling is needed.
				// We maybe could map the copyMode argument into some filterMode settings, but
				// i don't know the spec of copyMode, so ...
				PsychBlitTextureToDisplay(sourceWin, targetWin, sourceRect, targetRect, 0, 1, 1);
				
				// That's it.
				
				// Flush drawing commands and wait for render-completion in single-buffer mode:
				PsychFlushGL(targetWin);
		}
	}

	// Onscreen to texture copy?
	if (PsychIsOnscreenWindow(sourceWin) && PsychIsOffscreenWindow(targetWin)) {
		// With the current implemenation we can't zoom if sizes of sourceRect and targetRect don't
		// match: Only one-to-one copy possible...
		if(PsychGetWidthFromRect(sourceRect) != PsychGetWidthFromRect(targetRect) ||
		   PsychGetHeightFromRect(sourceRect) != PsychGetHeightFromRect(targetRect)) {
				// Non-matching sizes. We can't perform requested scaled copy :(
				PsychErrorExitMsg(PsychError_user, "Size mismatch of sourceRect and targetRect. Matching size is required for Onscreen to Offscreen copies. Sorry.");
		}
		
		// Update selected textures content:
		// Looks weird but we need the framebuffer of sourceWin:
		PsychSetDrawingTarget(sourceWin);

		// Disable alpha-blending:
		glDisable(GL_BLEND);
		
		// Disable any shading during copy-op:
		PsychSetShader(sourceWin, 0);

		// Texture objects are shared across contexts, so doesn't matter if targetWin's texture actually
		// belongs to the bound context of sourceWin:
		glBindTexture(PsychGetTextureTarget(targetWin), targetWin->textureNumber);
		
		// Copy into texture:
		glCopyTexSubImage2D(PsychGetTextureTarget(targetWin), 0, (int) targetRect[kPsychLeft], (int) (PsychGetHeightFromRect(targetWin->rect) - targetRect[kPsychBottom]), (int) sourceRect[kPsychLeft], (int) (PsychGetHeightFromRect(sourceWin->rect) - sourceRect[kPsychBottom]),
							(int) PsychGetWidthFromRect(sourceRect), (int) PsychGetHeightFromRect(sourceRect));
		
		// Unbind texture object:
		glBindTexture(PsychGetTextureTarget(targetWin), 0);
		
		// That's it.
		glPixelZoom(1,1);
	}

	// Onscreen to Onscreen copy?
	if (PsychIsOnscreenWindow(sourceWin) && PsychIsOnscreenWindow(targetWin)) {
		// Currently only works for copies of subregion -> subregion inside same onscreen window,
		// not across different onscreen windows! TODO: Only possible with EXT_framebuffer_blit
		if (sourceWin != targetWin) PsychErrorExitMsg(PsychError_user, "Sorry, the current implementation only supports copies within the same onscreen window, not accross onscreen windows.");

		// Set target windows framebuffer as drawing target:
		PsychSetDrawingTarget(targetWin);
		
		// Disable alpha-blending:
		glDisable(GL_BLEND);
		
		// Disable any shading during copy-op:
		PsychSetShader(targetWin, 0);

		// Start position for pixel write is:
		glRasterPos2f((float) targetRect[kPsychLeft], (float) targetRect[kPsychBottom]);
		
		// Zoom factor if rectangle sizes don't match:
		glPixelZoom((float) (PsychGetWidthFromRect(targetRect) / PsychGetWidthFromRect(sourceRect)), (float) (PsychGetHeightFromRect(targetRect) / PsychGetHeightFromRect(sourceRect)));
		
		// Perform pixel copy operation:
		glCopyPixels((int) sourceRect[kPsychLeft], (int) (PsychGetHeightFromRect(sourceWin->rect) - sourceRect[kPsychBottom]), (int) PsychGetWidthFromRect(sourceRect), (int) PsychGetHeightFromRect(sourceRect), GL_COLOR);
		
		// That's it.
		glPixelZoom(1,1);
		
		// Flush drawing commands and wait for render-completion in single-buffer mode:
		PsychFlushGL(targetWin);
	}

	// Just to make sure we catch invalid values:
	if (!(PsychPrefStateGet_ConserveVRAM() & kPsychAvoidCPUGPUSync)) PsychTestForGLErrors();

	// Done.
	return(PsychError_none);
}
