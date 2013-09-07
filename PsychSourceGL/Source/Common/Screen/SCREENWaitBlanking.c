/*
	SCREENWaitBlanking.c		
  
	AUTHORS:

		mario.kleiner at tuebingen.mpg.de 		mk 
  
	PLATFORMS:	
	
		All.
    

	HISTORY:

		01/23/06	mk		Created.
        05/31/11    mk      Add 3rd method of waiting for vblank, based on vblank counter
                            queries on supported systems (OS/X and Linux) to work around
                            bugs in bufferswap wait method on OS/X "Snow Leopard".
                            Code will try beampos waits, fallback to vblank counter waits,
                            then fallback to swapbuffers waits.

	DESCRIPTION:
  
                Waits for beginning of vertical retrace. Blocks Matlabs or PTB's execution until then.
                Screen('WaitBlanking') is here mostly for backwards-compatibility with old OS-9 and WinPTB.
                It is not the recommended way of doing things for new code!
 
	TO DO:  

*/


#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "framesSinceLastWait = Screen('WaitBlanking', windowPtr [, waitFrames]);";
static char synopsisString[] = 
        "Wait for specified number of monitor refresh intervals, stopping PTB's "
        "execution until then. Select waitFrames=1 (or omit it, since that's the default) "
        "to wait for the beginning of the next frame. "
        "\"windowPtr\" is the pointer to the onscreen window for which we should wait for. "
        "framesSinceLastWait contains the number of video refresh intervals from the last "
        "time Screen('WaitBlanking') or Screen('Flip') returned until return from "
        "this call to WaitBlanking. Please note that this function is only provided to keep old code "
        "from OS-9 PTB running. Use the Screen('Flip') command for all new code as this "
        "allows for much higher accuracy and reliability of stimulus timing and enables "
        "a huge number of new and very useful features! "
        "COMPATIBILITY TO OS-9 PTB: If you absolutely need to run old code for the old MacOS-9 or Windows "
        "Psychtoolbox, you can switch into a compatibility mode by adding the command "
        "Screen('Preference', 'EmulateOldPTB', 1) at the very top of your script. This will restore "
        "Offscreen windows and WaitBlanking functionality, but at the same time disable most of the new "
        "features of the OpenGL Psychtoolbox. Please do not write new experiment code in the old style! "
        "Emulation mode is pretty new and may contain significant bugs, so use with great caution!";  

static char seeAlsoString[] = "OpenWindow Flip Screen('Preference', 'EmulateOldPTB', 1)";

PsychError SCREENWaitBlanking(void) 
{
    PsychWindowRecordType *windowRecord;
    int waitFrames, framesWaited;
    double tvbl, ifi;
    long screenwidth, screenheight;
    int vbl_startline, beampos, lastline;
    psych_uint64 vblCount, vblRefCount;
    CGDirectDisplayID	cgDisplayID;
    GLint read_buffer, draw_buffer;
    
    // All subfunctions should have these two lines.  
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    PsychErrorExit(PsychCapNumInputArgs(2));     //The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); //The required number of inputs	
    PsychErrorExit(PsychCapNumOutputArgs(1));    //The maximum number of outputs
    
    // Get the window record from the window record argument and get info from the window record
    PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
    
    if(!PsychIsOnscreenWindow(windowRecord))
        PsychErrorExitMsg(PsychError_user, "Tried to call 'WaitBlanking' on something else than an onscreen window!");
    
    // Get the number of frames to wait:
    waitFrames = 0;
    PsychCopyInIntegerArg(2, FALSE, &waitFrames);

	// We default to wait at least one interval if no argument supplied:
    waitFrames = (waitFrames < 1) ? 1 : waitFrames;
    
    // Enable this windowRecords framebuffer as current drawingtarget:
    // This is needed to make sure that Offscreen windows work propely.
    PsychSetDrawingTarget(windowRecord);
    
    // Retrieve display handle for beamposition queries:
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, windowRecord->screenNumber);
    
    // Retrieve final vbl_startline, aka physical height of the display in pixels:
    PsychGetScreenSize(windowRecord->screenNumber, &screenwidth, &screenheight);
    vbl_startline = (int) screenheight;
    
    // Query duration of a monitor refresh interval: We try to use the measured interval,
	// but fallback of the nominal value reported by the operating system if necessary:
    if ((ifi = windowRecord->VideoRefreshInterval)<=0) {
        if (PsychGetNominalFramerate(windowRecord->screenNumber) > 0) {
            // Valid nominal framerate returned by OS: Calculate nominal IFI from it.
            ifi = 1.0 / ((double) PsychGetNominalFramerate(windowRecord->screenNumber));        
        }
        else {
            // No reasonable value available! We fallback to an assumed 60 Hz refresh...
            ifi = 1.0 / 60.0;
        }
    }
    
    // Query vblcount to test if this method works correctly:
    PsychOSGetVBLTimeAndCount(windowRecord, &vblRefCount);
    
    // Check if beamposition queries are supported by this OS and working properly:
    if (-1 != PsychGetDisplayBeamPosition(cgDisplayID, windowRecord->screenNumber) && windowRecord->VBL_Endline >= 0) {
        // Beamposition queries supported and fine. We can wait for VBL without bufferswap-tricks:
        // We query the rasterbeamposition and compare it
        // to the known values for the VBL area. If we enter VBL, we take a timestamp and return -
        // or wait for the next VBL if waitFrames>0
        
        // Query current beamposition when entering WaitBlanking:
        beampos = PsychGetDisplayBeamPosition(cgDisplayID, windowRecord->screenNumber);
        // Are we in VBL when entering WaitBlanking? If so, we should wait for one additional frame,
        // because by definition, WaitBlanking should always wait for at least one monitor refresh
        // interval...
        if ((beampos<=windowRecord->VBL_Endline) && (beampos>=vbl_startline)) waitFrames++;
        
        while(waitFrames > 0) {
            // Enough time for a sleep? If the beam is far away from VBL area, we try to sleep to
            // yield some CPU time to other processes in the system -- we are nice citizens ;)
            beampos = PsychGetDisplayBeamPosition(cgDisplayID, windowRecord->screenNumber);
            while (( ((float)(vbl_startline - beampos)) / (float) windowRecord->VBL_Endline * ifi) > 0.003) {
                // At least 3 milliseconds left until retrace. We sleep for 1 millisecond.
                PsychWaitIntervalSeconds(0.001);
                beampos = PsychGetDisplayBeamPosition(cgDisplayID, windowRecord->screenNumber);
            }
            
            // Less than 3 ms away from retrace. Busy-Wait for retrace...
			lastline = PsychGetDisplayBeamPosition(cgDisplayID, windowRecord->screenNumber);
			beampos = lastline;
			while ((beampos < vbl_startline) && (beampos >= lastline)) {
				lastline = beampos;
				beampos = (long) PsychGetDisplayBeamPosition(cgDisplayID, windowRecord->screenNumber);
			} 
			
            // Retrace! Take system timestamp of VBL onset:
            PsychGetAdjustedPrecisionTimerSeconds(&tvbl);
            
            // If this wasn't the last frame to wait, we need to wait for end of retrace before
            // repeating the loop, because otherwise we would detect the same VBL and skip frames.
            // If it is the last frame, we skip it and return as quickly as possible to save the
            // Matlab script some extra Millisecond for drawing...
            if (waitFrames>1) { 
                beampos = vbl_startline;
                while ((beampos<=windowRecord->VBL_Endline) && (beampos>=vbl_startline)) { beampos = PsychGetDisplayBeamPosition(cgDisplayID, windowRecord->screenNumber); };
            }
            
            // Done with this refresh interval...
            // Decrement remaining number of frames to wait:
            waitFrames--;
        }
    }
    else if (vblRefCount > 0) {
        // Display beamposition queries unsupported, but vblank count queries seem to work. Try those.
        // Should work on Linux and OS/X:
        while(waitFrames > 0) {
            vblCount = vblRefCount;

            // Wait for next vblank counter increment - aka start of next frame (its vblank):
            while (vblCount == vblRefCount) {
                // Requery:
                tvbl = PsychOSGetVBLTimeAndCount(windowRecord, &vblCount);
                // Yield at least 100 usecs. This is accurate as this code-path
                // only executes on OS/X and Linux, never on Windows (as of 01/06/2011):
                PsychYieldIntervalSeconds(0.000100);
            }

            vblRefCount = vblCount;
            
            // Done with this refresh interval...
            // Decrement remaining number of frames to wait:
            waitFrames--;
        }
    }
    else {            
        // Other methods unsupported. We use the doublebuffer swap method of waiting for retrace.
        //
        // Working principle: On each frame, we first copy the content of the (user visible) frontbuffer into the backbuffer.
        // Then we ask the OS to perform a front-backbuffer swap on next vertical retrace and go to sleep via glFinish() et al.
        // until the OS signals swap-completion. This way PTB's/Matlabs execution will stall until VBL, when swap happens and
        // we get woken up. We repeat this procedure 'waitFrames' times, then we take a high precision timestamp and exit the
        // Waitblanking loop. As backbuffer and frontbuffer are identical (due to the copy) at swap time, the visual display
        // won't change at all for the subject.
        // This method should work reliably, but it has one drawback: A random wakeup delay (scheduling jitter) is added after
        // retrace has been entered, so Waitblanking returns only after the beam has left retrace state on older hardware.
        // This means a bit less time (1 ms?) for stimulus drawing on Windows than on OS-X where Waitblanking returns faster. 
        
        // Child protection:
        if (windowRecord->windowType != kPsychDoubleBufferOnscreen) {
            PsychErrorExitMsg(PsychError_internal, "WaitBlanking tried to perform swap-waiting on a single buffered window!");
        }
        
        // Setup buffers for front->back copy op:
        
        // Backup old read- writebuffer assignments:
        glGetIntegerv(GL_READ_BUFFER, &read_buffer);
        glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
        
        // Set read- and writebuffer properly:
        glReadBuffer(GL_FRONT);
        glDrawBuffer(GL_BACK);
        
        // Reset viewport to full-screen default:
        glViewport(0, 0, screenwidth, screenheight);                
        glScissor(0, 0, screenwidth, screenheight);                
        
        // Reset color buffer writemask to "All enabled":
        glColorMask(TRUE, TRUE, TRUE, TRUE);
        glDisable(GL_BLEND);
        glPixelZoom(1,1);

		// Disable draw shader:
		PsychSetShader(windowRecord, 0);

        // Swap-Waiting loop for 'waitFrames' frames:
        while(waitFrames > 0) {
            // Copy current content of front buffer into backbuffer:
            glRasterPos2i(0, screenheight);
            glCopyPixels(0, 0, screenwidth, screenheight, GL_COLOR);            
            
            // Ok, front- and backbuffer are now identical, so a bufferswap
            // will be a visual no-op.
            
            // Enable beamsyncing of bufferswaps to VBL:
            PsychOSSetVBLSyncLevel(windowRecord, 1);                
            
            // Trigger bufferswap in sync with retrace:
            PsychOSFlipWindowBuffers(windowRecord);
            
            // Wait for swap-completion, aka beginning of VBL:
            PsychWaitPixelSyncToken(windowRecord);
            
            // VBL happened - Take system timestamp:
            PsychGetAdjustedPrecisionTimerSeconds(&tvbl);
            
            // This code-chunk is an alternate way of syncing, only used for debugging:
            if (false) {
                // Disable beamsyncing of bufferswaps to VBL:
                PsychOSSetVBLSyncLevel(windowRecord, 0);                
                // Swap buffers immediately without vsync:
                PsychOSFlipWindowBuffers(windowRecord);
            }
            
            // Decrement remaining number of frames to wait:
            waitFrames--;
        } // Do it again...
        
        // Enable beamsyncing of bufferswaps to VBL:
        PsychOSSetVBLSyncLevel(windowRecord, 1);                
        
        // Restore assignment of read- writebuffers and such:
        glEnable(GL_BLEND);
        glReadBuffer(read_buffer);
        glDrawBuffer(draw_buffer);        
        // Done with Windows waitblanking...
    }
    
    // Compute number of frames waited: It is timestamp of return of this waitblanking minus
    // timestamp of return of last waitblanking, divided by duration of a monitor refresh
    // interval, mathematically rounded to an integral number:
    framesWaited = (int) (((tvbl - windowRecord->time_at_last_vbl) / ifi) + 0.5f);
    
    // Update timestamp for next invocation of Waitblanking:
    windowRecord->time_at_last_vbl = tvbl;
    
    // Return optional number of frames waited:
    PsychCopyOutDoubleArg(1, FALSE, (double) framesWaited);
    
    // Done.
    return(PsychError_none);
}
