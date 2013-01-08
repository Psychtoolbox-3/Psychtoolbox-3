/*
	PsychToolbox3/Source/Common/Screen/PsychGraphicsHardwareHALSupport.c		

	AUTHORS:
	
		mario.kleiner@tuebingen.mpg.de		mk

	PLATFORMS:	
	
		All. However with dependencies on OS specific glue-layers which are mostly Linux/OSX for now...

	HISTORY:
	
	01/12/2008	mk		Created.
	
	DESCRIPTION:

	This file is a container for miscellaneous routines that take advantage of specific low level
	features of graphics/related hardware and the target operating system to achieve special tasks.
	
	Most of the routines here are more tied to specific displays (screens) than to windows and usually
	only a subset of these routines is available for a specific system configuration with a specific
	model of graphics card. Other layers of PTB should not rely on these routines being supported on
	a given system config and should be prepared to have fallback-implementations.
	
	Many of the features are experimental in nature!

	TO DO:

*/

#include "PsychGraphicsHardwareHALSupport.h"

// Include specifications of the GPU registers:
#include "PsychGraphicsCardRegisterSpecs.h"

// Maps screenid's to Graphics hardware pipelines: Used to choose pipeline for beampos-queries and similar
// GPU crtc specific stuff. Each screen can have up to kPsychMaxPossibleCrtcs assigned. Slot 0 contains the
// primary crtc, used for beamposition timestamping, framerate queries etc. A -1 value in a slot terminates
// the sequence of assigned crtc's.
static int	displayScreensToCrtcIds[kPsychMaxPossibleDisplays][kPsychMaxPossibleCrtcs];
static int	displayScreensToPipes[kPsychMaxPossibleDisplays][kPsychMaxPossibleCrtcs];
static int  numScreenMappings = 0;
static psych_bool displayScreensToPipesUserOverride = FALSE;
static psych_bool displayScreensToPipesAutoDetected = FALSE;

// Corrective values for beamposition queries to correct for any constant and systematic offsets in
// the scanline positions returned by low-level code:
static int	screenBeampositionBias[kPsychMaxPossibleDisplays];
static int	screenBeampositionVTotal[kPsychMaxPossibleDisplays];

/* PsychSynchronizeDisplayScreens() -- (Try to) synchronize display refresh cycles of multiple displays
 *
 * This tries whatever method is available/appropriate/or requested to synchronize the video refresh
 * cycles of multiple graphics cards physical display heads -- corresponding to PTB logical Screens.
 *
 * The method may or may not be supported on a specific OS/gfx-card combo. It will return a PsychError_unimplemented
 * if it can't do what core wants.
 *
 * numScreens	=	Ptr to the number of display screens to sync. If numScreens>0, all screens with the screenIds stored
 *					in the integer array 'screenIds' will be synched. If numScreens == 0, PTB will try to sync all
 *					available screens in the system. On return, the location will contain the count of synced screens.
 *
 * screenIds	=	Either a list with 'numScreens' screenIds for the screens to sync, or NULL if numScreens == 0.
 *
 * residuals	=	List with 'numScreens' (on return) values indicating the residual sync error wrt. to the first
 *					screen (the reference). Ideally all items should contain zero for perfect sync on return.
 *
 * syncMethod	=	Numeric Id for the sync method to use: 0 = Don't care, whatever is appropriate. 1 = Only hard
 *					sync, which is fast and reliable if supported. 2 = Soft sync via drift-syncing. More to come...
 *
 * syncTimeOut	=	If some non-immediate method is requested/chosen, it should give up after syncTimeOut seconds if
 *					it doesn't manage to bring the displays in sync in that timeframe.
 *
 * allowedResidual = How many scanlines offset after sync are acceptable? Will retry until syncTimeOut if criterion not met.
 *
 */
PsychError PsychSynchronizeDisplayScreens(int *numScreens, int* screenIds, int* residuals, unsigned int syncMethod, double syncTimeOut, int allowedResidual)
{
	// Currently, we only support a hard, immediate sync of all display heads of a single dual-head gfx-card,
	// so we ignore most of our arguments. Well, still check them for validity, but then ignore them after
	// successfull validation ;-)
	
	if (numScreens == NULL) PsychErrorExitMsg(PsychError_internal, "NULL-Ptr passed as numScreens argument!");
	if (*numScreens < 0 || *numScreens > PsychGetNumDisplays()) PsychErrorExitMsg(PsychError_internal, "Invalid number passed as numScreens argument! (Negative or more than available screens)");
	if (syncMethod < 0 || syncMethod > 2) PsychErrorExitMsg(PsychError_internal, "Invalid syncMethod argument passed!");
	if (syncTimeOut < 0) PsychErrorExitMsg(PsychError_internal, "Invalid (negative) syncTimeOut argument passed!");
	if (allowedResidual < 0) PsychErrorExitMsg(PsychError_internal, "Invalid (negative) allowedResidual argument passed!");
	
	// System support:
	#if PSYCH_SYSTEM == PSYCH_WINDOWS
		if(PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Synchronization of graphics display heads requested, but this is not supported on MS-Windows.\n");
		return(PsychError_unimplemented);
	#endif
	
	#if PSYCH_SYSTEM == PSYCH_LINUX
		// Dispatch to routine in PsychScreenGlue.c Linux:
		return(PsychOSSynchronizeDisplayScreens(numScreens, screenIds, residuals, syncMethod, syncTimeOut, allowedResidual));
	#endif
	
	#if PSYCH_SYSTEM == PSYCH_OSX
		// Dispatch to routine in PsychScreenGlue.c OSX:
		return(PsychOSSynchronizeDisplayScreens(numScreens, screenIds, residuals, syncMethod, syncTimeOut, allowedResidual));
	#endif
	
	// Often not reached...
	return(PsychError_none);
}

/* PsychSetOutputDithering() - Control bit depth control and dithering on digital display output encoder:
 * 
 * This function enables or disables bit depths truncation or dithering of digital display output ports of supported
 * graphics hardware. Currently the ATI Radeon X1000/HD2000/HD3000/HD4000/HD5000 and later cards should allow this.
 *
 * This needs support from the Psychtoolbox kernel level support driver for low-level register reads
 * and writes to the GPU registers.
 *
 *
 * 'windowRecord'	Is used to find the Id of the screen for which mode should be changed. If set to NULL then...
 * 'screenId'       ... is used to determine the screenId for the screen. Otherwise 'screenId' is ignored.
 * 'ditherEnable'   Zero = Disable any dithering. Non-Zero Reenable dithering after it has been disabled by us,
 *                  or if it wasn't disabled beforehand, enable it with a control mode as specified by the numeric
 *                  value of 'ditherEnable'. The value is GPU specific.
 *
 */
psych_bool  PsychSetOutputDithering(PsychWindowRecordType* windowRecord, int screenId, unsigned int ditherEnable)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX

	// Child protection:
	if (windowRecord && !PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_internal, "Invalid non-onscreen windowRecord provided!!!");
	
	// Either screenid from windowRecord or as passed in:
	if (windowRecord) screenId = windowRecord->screenNumber;
    
    // Do the call:
    PsychOSKDSetDitherMode(screenId, ditherEnable);

    return(TRUE);
#else
	// This cool stuff not supported on the uncool Windows OS:
    if(PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: GPU dithering control requested, but this is not supported on MS-Windows.\n");
	return(FALSE);
#endif
}

/* PsychSetGPUIdentityPassthrough() - Control identity passthrough of framebuffer 8 bpc pixel values to encoders/connectors:
 * 
 * This function enables or disables bit depths truncation or dithering of digital display output ports of supported
 * graphics hardware, and optionally loads a identity LUT into the hardware and configures other parts of the GPU's
 * color management for untampered passthrough of framebuffer pixels.
 * Currently the ATI Radeon X1000/HD2000/HD3000/HD4000/HD5000/HD6000 and later cards should allow this.
 *
 * This needs support from the Psychtoolbox kernel level support driver for low-level register reads
 * and writes to the GPU registers.
 *
 *
 * 'windowRecord'	Is used to find the Id of the screen for which mode should be changed. If set to NULL then...
 * 'screenId'       ... is used to determine the screenId for the screen. Otherwise 'screenId' is ignored.
 * 'passthroughEnable' Zero = Disable passthrough: Currently only reenables dithering, otherwise a no-op. 
 *                     1 = Enable passthrough, if possible.
 *
 * Returns:
 *
 * 0xffffffff if feature unsupported by given OS/Driver/GPU combo.
 * 0 = On failure to establish passthrough.
 * 1 = On partial success: Dithering disabled and identityt LUT loaded, but other GPU color transformation
 *                         features may not be configured optimally for passthrough.
 * 2 = On full success, as far as can be determined by software.
 *
 */
unsigned int PsychSetGPUIdentityPassthrough(PsychWindowRecordType* windowRecord, int screenId, psych_bool passthroughEnable)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
    unsigned int rc, rcret;
    int head, iter;
    
    // Init return code to "unsupported":
    rcret = 0xffffffff;
    
	// Child protection:
	if (windowRecord && !PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_internal, "Invalid non-onscreen windowRecord provided!!!");
	
	// Either screenid from windowRecord or as passed in:
	if (windowRecord) screenId = windowRecord->screenNumber;
    
    // Check if kernel driver is enabled, otherwise this won't work:
    if (!PsychOSIsKernelDriverAvailable(screenId)) {
        if(PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: GPU framebuffer passthrough setup requested, but this is not supported without kernel driver.\n");
        return(0xffffffff);
    }
    
    // Try to enable or disable dithering on display:
    PsychSetOutputDithering(windowRecord, screenId, (passthroughEnable) ? 0 : 1);
    
    // We're done if this an actual passthrough disable, as a full disable isn't yet implemented:
    if (!passthroughEnable) return(0);
    
    // Start with head undefined:
    head = -1;
    
    for (iter = 0; iter < kPsychMaxPossibleCrtcs; iter++) {
        if (screenId >= 0) {
            // Positive screenId: Apply to all crtc's for this screenId:
            
            // Is there an iter'th crtc assigned to this screen?
            head = PsychScreenToCrtcId(screenId, iter);
            
            // If end of list of associated crtc's for this screenId reached, then we're done:
            if (head < 0) break;
        }
        else {
            // Negative screenId -> Only affect one head defined by screenId:
            if (head < 0) {
                // Setup single target head in this iteration:
                head = -screenId;
            }
            else {
                // Single target head already set up: We're done:
                break;
            }
        }
        
        // Check if remaining GPU is already configured for untampered identity passthrough:
        rc = PsychOSKDGetLUTState(screenId, head, (PsychPrefStateGet_Verbosity() > 4) ? 1 : 0);
        if(PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: [screen %i, head %i] 1st LUT query rc = %i.\n", screenId, head, rc);
        if (rc == 0xffffffff) {
            // Unsupported for this GPU. We're done:
            if(PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: GPU framebuffer passthrough setup on screenid %i, head %i requested, but this is not supported on this GPU.\n", screenId, head);
            rcret = 0xffffffff;
            continue;
        }
        
        // Perfect identity passthrough already configured?
        if (rc == 2) {
            // Yes. We're successfully done!
            if(PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: GPU framebuffer passthrough setup I completed on screenid %i, head %i. Perfect passthrough should work now.\n", screenId, head);
            rcret = 2;
            continue;
        }
        
        // No. Try to setup GPU for passthrough:
        if (!PsychOSKDLoadIdentityLUT(screenId, head)) {
            // Failed.
            if(PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: GPU framebuffer passthrough setup requested on screenid %i, head %i, but setup failed.\n", screenId, head);
            rcret = 0xffffffff;
            continue;
        }
        
        // Make sure, GPU's gamma table can settle by waiting 250 msecs:
        PsychYieldIntervalSeconds(0.250);
        
        // Setup supposedly successfully finished. Re-Query state:
        rc = PsychOSKDGetLUTState(screenId, head, (PsychPrefStateGet_Verbosity() > 4) ? 1 : 0);
        if(PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: [screen %i, head %i] 2nd LUT query rc = %i.\n", screenId, head, rc);
        
        // Perfect identity passthrough now configured?
        if (rc == 2) {
            // Yes. We're successfully done!
            if(PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: GPU framebuffer passthrough setup II completed on screenid %i, head %i. Perfect passthrough should work now.\n", screenId, head);
            rcret = 2;
            continue;
        }
        
        if (rc == 3) {
            // Not quite. We've done what we could. A perfect identity LUT is setup, but the rest of the hw
            // isn't in that a great shape. This may or may not be good enough...
            if(PsychPrefStateGet_Verbosity() > 3) {
                printf("PTB-INFO: GPU framebuffer passthrough setup II completed on screenid %i, head %i. Sort of ok passthrough achieved. May or may not work.\n", screenId, head);
                printf("PTB-INFO: A perfect identity gamma table is loaded, but the other GPU color transformation settings are still suboptimal.\n");
            }
            rcret = 1;
            continue;
        }
        
        // Ok, we failed.
        if(PsychPrefStateGet_Verbosity() > 3) {
            printf("PTB-INFO: GPU framebuffer passthrough setup II completed on screenid %i, head %i. Failed to establish identity passthrough!\n", screenId, head);
            printf("PTB-INFO: Could not upload a perfect identity LUT. May still work due to hopefully disabled dithering, who knows?\n");
        }

        rcret = 0;
    }
    
    return(rcret);
#else
	// This cool stuff not supported on the uncool Windows OS:
    if(PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: GPU framebuffer passthrough setup requested, but this is not supported on MS-Windows.\n");
	return(0xffffffff);
#endif
}

/* PsychEnableNative10BitFramebuffer()  - Enable/Disable native 10 bpc RGB framebuffer modes.
 *
 * This function enables or disables the native ARGB2101010 framebuffer readout mode of supported
 * graphics hardware. Currently the ATI Radeon X1000/HD2000/HD3000 and later cards should allow this.
 *
 * This needs support from the Psychtoolbox kernel level support driver for low-level register reads
 * and writes to the GPU registers.
 *
 * 'windowRecord'	Is used to find the Id of the screen for which mode should be changed, as well as enable
 *					flags to see if a change is required at all, and the OpenGL context for some specific
 *					fixups. A value of NULL will try to apply the operation to all heads, but may only work
 *					for *disabling* 10 bpc mode, not for enabling it -- Mostly useful for a master reset to
 *					system default, e.g., as part of error handling or Screen shutdown handling.
 * 'enable'   True = Enable ABGR2101010 support, False = Disable ARGB2101010 support, reenable ARGB8888 support. 
 *
 */
psych_bool	PsychEnableNative10BitFramebuffer(PsychWindowRecordType* windowRecord, psych_bool enable)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
	int i, si, ei, headid, headiter, screenId;
	unsigned int lutreg, ctlreg, value, status;
	int gpuMaintype, gpuMinortype;

	// Child protection:
	if (windowRecord && !PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_internal, "Invalid non-onscreen windowRecord provided!!!");
	
	// Either screenid from windowRecord or our special -1 "all Screens" Id:
	screenId = (windowRecord) ? windowRecord->screenNumber : -1;

	// We only support Radeon GPU's with AVIVO display engine, aka DCE-1, nothing more recent:
	if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, NULL) ||
	    (gpuMaintype != kPsychRadeon) || (gpuMinortype > 0x10)) {
	  return(FALSE);
	}
	
	// Define range of screens: Either a single specific one, or all:
	si = (screenId!=-1) ? screenId   : 0;
	ei = (screenId!=-1) ? screenId+1 : PsychGetNumDisplays();

	// Loop over all target screens:
	for (i = si; i < ei; i++) {
		// Iterate over range of all assigned heads for this screenId 'i' and reconfigure them:
		for (headiter = 0; headiter < kPsychMaxPossibleCrtcs; headiter++) {
            // Map screenid to headid for headiter'th head:
            headid = PsychScreenToCrtcId(i, headiter);

            // We're done as soon as we encounter invalid negative headid.
            if (headid < 0) break;

			// Select Radeon HW registers for corresponding heads:
			lutreg = (headid == 0) ? RADEON_D1GRPH_LUT_SEL : RADEON_D2GRPH_LUT_SEL;
			ctlreg = (headid == 0) ? RADEON_D1GRPH_CONTROL : RADEON_D2GRPH_CONTROL;

			// Enable or Disable?
			if (enable) {
				// Enable:
			
				// Switch hardware LUT's to bypass mode:
				// We set bit 8 to enable "bypass LUT in 2101010 mode":
				value = PsychOSKDReadRegister(screenId, lutreg, &status);
				if (status) {
					printf("PTB-ERROR: Failed to set 10 bit framebuffer mode (LUTReg read failed).\n");
					return(false);
				}

				// Set the bypass bit:
				value = value | 0x1 << 8;

				PsychOSKDWriteRegister(screenId, lutreg, value, &status);
				if (status) {
					printf("PTB-ERROR: Failed to set 10 bit framebuffer mode (LUTReg write failed).\n");
					return(false);
				}
			
				// Only reconfigure framebuffer scanout if this is really our true Native10bpc hack:
				// This is usually skipped on FireGL/FirePro GPU's as their drivers do it already...
				if (windowRecord->specialflags & kPsychNative10bpcFBActive) {
					// Switch CRTC to ABGR2101010 readout mode:
					// We set bit 8 to enable that mode:
					value = PsychOSKDReadRegister(screenId, ctlreg, &status);
					if (status) {
						printf("PTB-ERROR: Failed to set 10 bit framebuffer mode (CTLReg read failed).\n");
						return(false);
					}
                
					// Set 2101010 mode bit:
					value = value | 0x1 << 8;
                
					PsychOSKDWriteRegister(screenId, ctlreg, value, &status);
					if (status) {
						printf("PTB-ERROR: Failed to set 10 bit framebuffer mode (CTLReg write failed).\n");
						return(false);
					}
				}
            
				// Pipe should be in 10 bpc mode now...
				if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: System framebuffer switched to ARGB2101010 mode for screen %i [head %i].\n", i, headid);
			} else {
				// Disable:

				// Only reconfigure framebuffer scanout if this is really our true Native10bpc hack:
				// This is usually skipped on FireGL/FirePro GPU's as their drivers do it already...
				if (windowRecord->specialflags & kPsychNative10bpcFBActive) {
					// Switch CRTC to ABGR8888 readout mode:
					// We clear bit 8 to enable that mode:
					value = PsychOSKDReadRegister(screenId, ctlreg, &status);
					if (status) {
						// This codepath gets always called in PsychCloseWindow(), so we should skip it
						// silently if register read fails, as this is expected on MS-Windows and on all
						// non-Radeon hardware and if kernel driver isn't loaded:
						if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-ERROR: Failed to set 8 bit framebuffer mode (CTLReg read failed).\n");
						return(false);
					}
					else if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In disable 10bpc: Readreg. ctlreg yields %lx\n", value);
                
					// Clear 2101010 mode bit:
					value = value & ~(0x1 << 8);
                
					PsychOSKDWriteRegister(screenId, ctlreg, value, &status);
					if (status) {
						printf("PTB-ERROR: Failed to set 8 bit framebuffer mode (CTLReg write failed).\n");
						return(false);
					}
					else if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In disable 10bpc: ctlreg reset\n");
                
					// Wait 500 msecs for GPU to settle:
					PsychWaitIntervalSeconds(0.5);
				}
            
				// Switch hardware LUT's to standard mapping mode:
				// We clear bit 8 to disable "bypass LUT in 2101010 mode":
				value = PsychOSKDReadRegister(screenId, lutreg, &status);
				if (status) {
					printf("PTB-ERROR: Failed to set 8 bit framebuffer mode (LUTReg read failed).\n");
					return(false);
				}
				else if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In disable 10bpc: Readreg. lutreg yields %lx\n", value);

				// Clear LUT bypass bit:
				value = value & ~(0x1 << 8);

				PsychOSKDWriteRegister(screenId, lutreg, value, &status);
				if (status) {
					printf("PTB-ERROR: Failed to set 8 bit framebuffer mode (LUTReg write failed).\n");
					return(false);
				}
				else if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In disable 10bpc: lutreg reset\n");

				// Pipe should be in 8 bpc mode now...
				if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: System framebuffer switched to standard ARGB8888 mode for screen %i [head %i].\n", i, headid);
			}
		} // Next display head...
	} // Next screenId.

	// Done.
	return(TRUE);
	
#else
	// This cool stuff not supported on the uncool Windows OS:
	return(FALSE);
#endif
}

/* PsychFixupNative10BitFramebufferEnableAfterEndOfSceneMarker()
 *
 * Undo changes made by the graphics driver to the framebuffer pixel format control register
 * as part of an OpenGL/Graphics op that marks "End of Scene", e.g., a glClear() command, that
 * would revert the framebuffers opmode to standard 8bpc mode and thereby kill our 10 bpc mode
 * setting.
 *
 * This routine *must* be called after each such problematic "End of scene" marker command like
 * glClear(). The routine does nothing if 10bpc mode is not enabled/requested for the corresponding
 * display head associated with the given onscreen window. It rewrites the control register on
 * 10 bpc configured windows to basically undo the unwanted change of the gfx-driver *before*
 * a vertical retrace cycle starts, ie., before that changes take effect (The register is double-
 * buffered and latched to update only at VSYNC time, so we just have to be quick enough).
 *
 *
 * Expected Sequence of operations is:
 * 1. Some EOS command like glClear() issued.
 * 2. EOS command schedules ctrl register update to "bad" value at next VSYNC.
 * 3. This routine gets called, detects need for fixup, glGetError() waits for "2." to finish.
 * 4. This routine undos the "bad" value change request by overwriting the latch with our
 *    "good" value --> Scheduled for next VSYNC. Then it returns...
 * 5. At next VSYNC or old "good" value is overwritten/latched with our new old "good" value,
 *    --> "good value" persists, framebuffer stays in 2101010 configuration --> All good.
 *
 * So far the theory, let's see if this really works in real world...
 *
 * This is not needed in Carbon+AGL windowed mode, as the driver doesnt' mess with the control
 * register there, but that mode has its own share of drawback, e.g., generally reduced performance
 * and less robust stimulus onset timing and timestamping... Life's full of tradeoffs...
 */
void PsychFixupNative10BitFramebufferEnableAfterEndOfSceneMarker(PsychWindowRecordType* windowRecord)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX

	int headiter, headid, screenId;
	unsigned int ctlreg;
	int gpuMaintype, gpuMinortype;

	// Fixup needed? Only if 10bpc mode is supposed to be active! Early exit if not:
	if (!(windowRecord->specialflags & kPsychNative10bpcFBActive)) return;

	// Map windows screen to gfx-headid aka register subset. TODO : We'll need something better,
	// more generic, abstracted out for the future, but as a starter this will do:
	screenId = windowRecord->screenNumber;

	// We only support Radeon GPU's with AVIVO display engine, aka DCE-1, nothing more recent:
	if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, NULL) ||
	    (gpuMaintype != kPsychRadeon) || (gpuMinortype > 0x10)) {
	  return;
	}

	// This command must be called with the OpenGL context of the given windowRecord active, so
	// we can rely on glGetError() waiting for the correct pipeline to settle! Wait via glGetError()
	// for the end-of-scene marker to finish completely, so our register write happens after
	// the "wrong" register write of that command. glFinish() doesn't work here for unknown
	// reasons - probably it waits too long or whatever. Pretty shaky this stuff...
	glGetError();
	
	// Ok, now rewrite the double-buffered (latched) register with our "good" value for keeping
	// the 10 bpc framebuffer online:
	
	// Iterate over range of all assigned heads for this screenId 'i' and reconfigure them:
	for (headiter = 0; headiter < kPsychMaxPossibleCrtcs; headiter++) {
	  // Map screenid to headid for headiter'th head:
	  headid = PsychScreenToCrtcId(screenId, headiter);
        
	  // We're done as soon as we encounter invalid negative headid.
	  if (headid < 0) break;
        
	  ctlreg = (headid == 0) ? RADEON_D1GRPH_CONTROL : RADEON_D2GRPH_CONTROL;
        
	  // One-liner read-modify-write op, which simply sets bit 8 of the register - the "Enable 2101010 mode" bit:
	  PsychOSKDWriteRegister(screenId, ctlreg, (0x1 << 8) | PsychOSKDReadRegister(screenId, ctlreg, NULL), NULL);
        
	  // Debug output, if wanted:
	  if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychFixupNative10BitFramebufferEnableAfterEndOfSceneMarker(): ARGB2101010 bit set on screen %i, head %i.\n", screenId, headid);
	}
#endif

	// Done.
	return;
}

/* Stores content of GPU's surface address registers of the surfaces that
 * correspond to the windowRecords frontBuffers. Only called inside
 * PsychFlipWindowBuffers() immediately before triggering a double-buffer
 * swap. The values are used as reference values: If another readout of
 * these registers shows values different from the ones stored preflip,
 * then that is a certain indicator that bufferswap has happened.
 */
void PsychStoreGPUSurfaceAddresses(PsychWindowRecordType* windowRecord)
{

#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
	int gpuMaintype, gpuMinortype;

	// If we are called, we know that 'windowRecord' is an onscreen window.
	int screenId = windowRecord->screenNumber;

	// Just need to check if GPU low-level access is supported:
	if (!PsychOSIsKernelDriverAvailable(screenId)) return;

	// We only support Radeon GPU's with AVIVO display engine, aka DCE-1, nothing more recent:
	if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, NULL) ||
	    (gpuMaintype != kPsychRadeon) || (gpuMinortype > 0x10)) {
	  return;
	}
	
	// Driver is online: Read the registers, but only for primary crtc in a multi-crtc config:
	windowRecord->gpu_preflip_Surfaces[0] = PsychOSKDReadRegister(screenId, (PsychScreenToCrtcId(screenId, 0) < 1) ? RADEON_D1GRPH_PRIMARY_SURFACE_ADDRESS : RADEON_D2GRPH_PRIMARY_SURFACE_ADDRESS, NULL);
	windowRecord->gpu_preflip_Surfaces[1] = PsychOSKDReadRegister(screenId, (PsychScreenToCrtcId(screenId, 0) < 1) ? RADEON_D1GRPH_SECONDARY_SURFACE_ADDRESS : RADEON_D2GRPH_SECONDARY_SURFACE_ADDRESS, NULL);

#endif

	return;
}

/*  PsychWaitForBufferswapPendingOrFinished()
 *  Waits until a bufferswap for window windowRecord has either already happened or
 *  bufferswap is certain.
 *  Input values:
 *  windowRecord struct of onscreen window to monitor.
 *  timestamp    = Deadline for abortion of flip detection at input.
 *
 *  Return values:
 *  timestamp    = System time at polling loop exit.
 *  beamposition = Beamposition at polling loop exit.
 *
 *  Return value: FALSE if swap happened already, TRUE if swap is imminent.
 */
psych_bool PsychWaitForBufferswapPendingOrFinished(PsychWindowRecordType* windowRecord, double* timestamp, int *beamposition)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
  int gpuMaintype, gpuMinortype;
  
  CGDirectDisplayID displayID;
  unsigned int primarySurface, secondarySurface;
  unsigned int updateStatus;
  double deadline = *timestamp;

  // If we are called, we know that 'windowRecord' is an onscreen window.
  int screenId = windowRecord->screenNumber;

  // Retrieve display id and screen size spec that is needed later...
  PsychGetCGDisplayIDFromScreenNumber(&displayID, screenId);

#define RADEON_D1GRPH_UPDATE	0x6144
#define RADEON_D2GRPH_UPDATE	0x6944
#define RADEON_SURFACE_UPDATE_PENDING 4
#define RADEON_SURFACE_UPDATE_TAKEN   8

	// Just need to check if GPU low-level access is supported:
	if (!PsychOSIsKernelDriverAvailable(screenId)) return(FALSE);

	// We only support Radeon GPU's with AVIVO display engine, aka DCE-1, nothing more recent:
	if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, NULL) ||
	    (gpuMaintype != kPsychRadeon) || (gpuMinortype > 0x10)) {
	  return(FALSE);
	}
	
	// Driver is online. Enter polling loop:
	while (TRUE) {
		// Read surface address registers:
		primarySurface   = PsychOSKDReadRegister(screenId, (PsychScreenToCrtcId(screenId, 0) < 1) ? RADEON_D1GRPH_PRIMARY_SURFACE_ADDRESS : RADEON_D2GRPH_PRIMARY_SURFACE_ADDRESS, NULL);
		secondarySurface = PsychOSKDReadRegister(screenId, (PsychScreenToCrtcId(screenId, 0) < 1) ? RADEON_D1GRPH_SECONDARY_SURFACE_ADDRESS : RADEON_D2GRPH_SECONDARY_SURFACE_ADDRESS, NULL);

		// Read update status registers:
		updateStatus     = PsychOSKDReadRegister(screenId, (PsychScreenToCrtcId(screenId, 0) < 1) ? RADEON_D1GRPH_UPDATE : RADEON_D2GRPH_UPDATE, NULL);

		PsychGetAdjustedPrecisionTimerSeconds(timestamp);

		if (primarySurface!=windowRecord->gpu_preflip_Surfaces[0] || secondarySurface!=windowRecord->gpu_preflip_Surfaces[1] || (updateStatus & (RADEON_SURFACE_UPDATE_PENDING | RADEON_SURFACE_UPDATE_TAKEN)) || (*timestamp > deadline)) {
			// Abort condition: Exit loop.
			break;
		}
		
		if (PsychPrefStateGet_Verbosity() > 9) {
			printf("PTB-DEBUG: Head %i: primarySurface=%p : secondarySurface=%p : updateStatus=%i\n", PsychScreenToCrtcId(screenId, 0), primarySurface, secondarySurface, updateStatus);
		}

		// Sleep slacky at least 200 microseconds, then retry:
		PsychYieldIntervalSeconds(0.0002);
	};
	
	// Take timestamp and beamposition:
	*beamposition = PsychGetDisplayBeamPosition(displayID, screenId);
	PsychGetAdjustedPrecisionTimerSeconds(timestamp);

	// Exit due to timeout?
	if (*timestamp > deadline) {
		// Mark timestamp as invalid due to timeout:
		*timestamp = -1;
	}
	
	// Return FALSE if bufferswap happened already, TRUE if swap is still pending:
	return((updateStatus & RADEON_SURFACE_UPDATE_PENDING) ? TRUE : FALSE);
#else
	// On Windows, we always return "swap happened":
	return(FALSE);
#endif
}

/* PsychGetNVidiaGPUType()
 *
 * Decodes hw register of NVidia GPU into GPU core id / chip family:
 * Returns 0 for unknown card, otherwise xx for NV_xx:
 */
unsigned int PsychGetNVidiaGPUType(PsychWindowRecordType* windowRecord)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
	psych_uint32 chipset, card_type;

	// Get hardware id code from gpu register:
	psych_uint32 reg0 = PsychOSKDReadRegister((windowRecord) ? windowRecord->screenNumber : 0, NV03_PMC_BOOT_0, NULL);
	
	/* We're dealing with >=NV10 */
	if ((reg0 & 0x0f000000) > 0) {
		/* Bit 27-20 contain the architecture in hex */
		chipset = (reg0 & 0xff00000) >> 20;
		/* NV04 or NV05 */
	} else if ((reg0 & 0xff00fff0) == 0x20004000) {
		if (reg0 & 0x00f00000)
			chipset = 0x05;
		else
			chipset = 0x04;
	} else
		chipset = 0xff;
	
	switch (chipset & 0xf0) {
		case 0x00:
			// NV_04/05: RivaTNT , RivaTNT2
			card_type = 0x04;
			break;
		case 0x10:
		case 0x20:
		case 0x30:
			// NV30 or earlier: GeForce-5 / GeForceFX and earlier:
			card_type = chipset & 0xf0;
			break;
		case 0x40:
		case 0x60:
			// NV40: GeForce6/7 series:
			card_type = 0x40;
			break;
		case 0x50:
		case 0x80:
		case 0x90:
		case 0xa0:
			// NV50: GeForce8/9/G100-G300.
			card_type = 0x50;
			break;
		case 0xc0:
			// Fermi: GeForce G400/500 series:
			card_type = 0xc0;
			break;
		case 0xe0:
			// Kepler: GeForce G600+ series:
			card_type = 0xe0;
			break;
		default:
			printf("PTB-DEBUG: Unknown NVidia chipset 0x%08x \n", reg0);
			card_type = 0x00;
	}
	
	return(card_type);
#else
	return(0);
#endif
}

/* PsychScreenToHead() - Map PTB screenId to GPU output headId (aka pipeId):
 *
 * See explanations for PsychScreenToCrtcId() to understand what this is good for!
 *
 * screenId = PTB screen index.
 * rankId = Select which head in a multi-head config. rankId 0 == Primary output.
 * A return value of -1 for a given rankId means that no such output is assigned,
 * it terminates the array.
 */
int	PsychScreenToHead(int screenId, int rankId)
{
	return(displayScreensToPipes[screenId][rankId]);
}

/* PsychSetScreenToHead() - Change mapping of a PTB screenId to GPU headId: */
void PsychSetScreenToHead(int screenId, int headId, int rankId)
{
    // Assign new mapping:
	displayScreensToPipes[screenId][rankId] = headId;
    
    // Mark mappings as user-defined instead of auto-detected/default-setup:
    displayScreensToPipesUserOverride = TRUE;
}

/* PsychScreenToCrtcId()
 *
 * Map PTB screenId and output head id to the index of the associated low-level
 * crtc scanout engine of the GPU: rankId selects which output head (0 = primary).
 *
 * PsychScreenToHead() returns the os-specific identifier of a specific
 * display output head, e.g., a display connector. On Windows and OS/X this is currently
 * simply a running number: 0 for the first display output, 1 for the second etc. On
 * Linux/X11 this is the X11 RandR extension protocol XID of the crtc associated
 * with a given display output, which allows to use the RandR extension to address
 * specific crtc's and do things like query and set video mode of a crtc (resolution,
 * video refresh rate), viewport of a crtc, rotation, mirroring state and other
 * geometric transforms, backlight and dithering settings etc. A XID of zero, which means
 * "invalid/not assigned" gets mapped to -1 for compatibility reasons in PTB.
 *
 * PsychScreenToCrtcId() returns the operating system independent, but gpu-specific index
 * of the low-level crtc display scanout engine associated with a display output. The
 * naming convention here is purely Psychtoolbox specific, as this index is used for
 * low-level direct access to GPU MMIO control registers via PTB's own magic. Values
 * are -1 for "not assigned/invalid" and then 0, 1, 2, ... for scanout engine zero, one,
 * two, ... These numbers are mapped in a gpu specific way to the addresses and offsets
 * of low-level control registers of the GPU hardware.
 *
 * Unfortunately, operating systems don't provide any well defined means to find out the
 * mapping between PsychScreenToHead() "high-level" output id's and PsychScreenToCrtcId()
 * low-level crtc id's, so the mapping gets determined at Screen() startup time via some more
 * or less clever heuristics which should do the right thing(tm) for common display and gpu
 * setups, but may fail on exotic configs. To cope with those, manual overrides are provided to
 * usercode, so the user can hopefully figure out correct mappings via trial and error.
 */
int	PsychScreenToCrtcId(int screenId, int rankId)
{
	return(displayScreensToCrtcIds[screenId][rankId]);
}

void PsychSetScreenToCrtcId(int screenId, int crtcId, int rankId)
{
    // Assign new mapping:
	displayScreensToCrtcIds[screenId][rankId] = crtcId;
}

/* PsychInitScreenToHeadMappings() - Setup initial mapping for 'numDisplays' displays:
 *
 * Called from end of InitCGDisplayIDList() during os-specific display initialization.
 *
 * 1. Starts with an identity mapping screen 0 -> (head 0 / crtcid 0), screen 1 -> (head 1 / crtcid 1) ...
 *
 * 2. Allows override of low-level crtc id mapping of the first output of a screen via
 *    environment variable "PSYCHTOOLBOX_PIPEMAPPINGS".
 *
 *    Format is: One character (a number between "0" and "9") for each screenid,
 *    e.g., "021" would map screenid 0 to crtcid 0, screenid 1 to crtcid 2 and screenid 2 to crtcid 1.
 *
 * 3. This mapping can be overriden via Screen('Preference', 'ScreenToHead') setting.
 *
 */
void PsychInitScreenToHeadMappings(int numDisplays)
{
    int i, j;
    char* ptbpipelines = NULL;
    
    displayScreensToPipesAutoDetected = FALSE;
    
    // Setup default identity one-to-one mapping:
    for(i = 0; i < kPsychMaxPossibleDisplays; i++){
	displayScreensToPipes[i][0] = i;
        displayScreensToCrtcIds[i][0] = i;
   
	for (j = 1; j < kPsychMaxPossibleCrtcs; j++) {
            displayScreensToPipes[i][j] = -1;
            displayScreensToCrtcIds[i][j] = -1;   
        }

	// We also setup beamposition bias values to "neutral defaults":
	screenBeampositionBias[i] = 0;
	screenBeampositionVTotal[i] = 0;
    }
	
    // Did user provide an override for the screenid --> pipeline mapping?
    ptbpipelines = getenv("PSYCHTOOLBOX_PIPEMAPPINGS");
    if (ptbpipelines) {
        // The default is "012...", ie screen 0 = pipe 0, 1 = pipe 1, 2 =pipe 2, n = pipe n
        for (i = 0; (i < (int) strlen(ptbpipelines)) && (i < kPsychMaxPossibleDisplays); i++) {
            PsychSetScreenToCrtcId(i, (((ptbpipelines[i] - 0x30) >=0) && ((ptbpipelines[i] - 0x30) < 10)) ? (ptbpipelines[i] - 0x30) : -1, 0);
        }
    }

    // Store number of mapping entries internally:
    numScreenMappings = numDisplays;
}

// Try to auto-detect screen to head mappings if possible and not yet overriden by usercode:
void PsychAutoDetectScreenToHeadMappings(int maxHeads)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX

    float nullTable[256];
    int screenId, headId, numEntries;
    float *redTable, *greenTable, *blueTable;

    // MK FIXME TODO: DISABLED FOR NOW!
    // As far as i understand, the all-zero gamma tables that are loaded into the crtc's by this routine get
    // "sticky". Our low-level lut readback code reads "all zeros" gamma tables back long after the gamma tables
    // have been restored to their original settings by the OS high level code -- the gpu is "lying" to us
    // about its hardware state :-(
    // This needs more careful examination. Until then we disable auto-detection. With the hard-coded default
    // settings - which are often correct and user-tweakable - the identity passthrough setup code for devices
    // like Bits+ and Datapixx seems to work ok.
    return;
    
    // If user / usercode has provided manual mapping, i.e., overriden the
    // default identity mapping, then we don't do anything, but accept the
    // users choice instead. Also skip this if it has been successfully executed
    // already:
    if (displayScreensToPipesUserOverride || displayScreensToPipesAutoDetected) return;
    
    // nullTable is our all-zero gamma table:
    memset(&nullTable[0], 0, sizeof(nullTable));

    // Ok, iterate over all logical screens and try to update
    // their mapping:
    for (screenId = 0; screenId < numScreenMappings; screenId++) {
        // Kernel driver for this screenId enabled? Otherwise we skip it:
        if (!PsychOSIsKernelDriverAvailable(screenId)) continue;
        
        // Yes. Perform detection sequence:
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Trying to detect screenId to display head mapping for screenid %i ...", screenId);
        
        // Retrieve current gamma table. Need to back it up internally:
        PsychReadNormalizedGammaTable(screenId, -1, &numEntries, &redTable, &greenTable, &blueTable);
        
        // Now load an all-zero gamma table for that screen:
        PsychLoadNormalizedGammaTable(screenId, -1, 256, nullTable, nullTable, nullTable);
        
        // Wait for 100 msecs, so the gamma table has actually settled (e.g., if its update was
        // delayed until next vblank on a >= 20 Hz display):
        PsychYieldIntervalSeconds(0.100);
        
        // Check all display heads to find the null table:
        for (headId = 0; headId < maxHeads; headId++) {
            if (PsychOSKDGetLUTState(screenId, headId, 0) == 1) {
                // Got it. Store mapping:
                displayScreensToPipes[screenId][0] = headId;
                
                // Done with searching:
                if (PsychPrefStateGet_Verbosity() > 2) printf(" found headId %i.", headId);
                break;
            }
        } 
        
        // Now restore original gamma table for that screen:
        PsychLoadNormalizedGammaTable(screenId, -1, numEntries, redTable, greenTable, blueTable);
        
        // Wait for 100 msecs, so the gamma table has actually settled (e.g., if its update was
        // delayed until next vblank on a >= 20 Hz display):
        PsychYieldIntervalSeconds(0.100);        

        if (PsychPrefStateGet_Verbosity() > 2) printf(" Done.\n");
    }
    
    // Done.
    displayScreensToPipesAutoDetected = TRUE;

#endif

    return;
}

/* PsychGetBeamposCorrection() -- Get corrective beamposition values.
 * Some GPU's and drivers don't return the true vertical scanout position on
 * query, but a value that is offset by a constant value (for a give display mode).
 * This function returns corrective values to apply to the GPU returned values
 * to get the "true scanoutposition" for timestamping etc.
 *
 * Proper values are setup via PsychSetBeamposCorrection() from high-level startup code
 * if needed. Otherwise they are set to (0,0), so the correction is an effective no-op.
 *
 * truebeampos = measuredbeampos - *vblbias;
 * if (truebeampos < 0) truebeampos = *vbltotal + truebeampos;
 *
 */
void PsychGetBeamposCorrection(int screenId, int *vblbias, int *vbltotal)
{
	*vblbias  = screenBeampositionBias[screenId];
	*vbltotal = screenBeampositionVTotal[screenId];
}

/* PsychSetBeamposCorrection() -- Set corrective beamposition values.
 * Called from high-level setup/calibration code at onscreen window open time.
 */
void PsychSetBeamposCorrection(int screenId, int vblbias, int vbltotal)
{
	// Need head id of primary crtc of this screen for auto-detection:
	int crtcid = PsychScreenToCrtcId(screenId, 0);
	
	// Auto-Detection of correct values requested? A valid OpenGL context must
	// be bound for this to work or we will crash horribly:
	if (((unsigned int) vblbias == 0xffffffff) && ((unsigned int) vbltotal == 0xffffffff)) {
		// First set'em to neutral safe values in case we fail our auto-detect:
		vblbias  = 0;
		vbltotal = 0;
		
		// Can do this on NVidia GPU's >= NV-50 if low-level access (PTB kernel driver or equivalent) is enabled:
		if ((strstr((char*) glGetString(GL_VENDOR), "NVIDIA") || strstr((char*) glGetString(GL_VENDOR), "nouveau") ||
			strstr((char*) glGetString(GL_RENDERER), "NVIDIA") || strstr((char*) glGetString(GL_RENDERER), "nouveau")) &&
			PsychOSIsKernelDriverAvailable(screenId)) {

			// Need to read different regs for NV-50 and later:
			if ((PsychGetNVidiaGPUType(NULL) >= 0x50) || (PsychGetNVidiaGPUType(NULL) == 0x0)) {
				// Auto-Detection. Read values directly from NV-50 class and later hardware:
				//
				// SYNC_START_TO_BLANK_END 16 bit high-word in CRTC_VAL block of NV50_PDISPLAY on NV-50 encodes
				// length of interval from vsync start line to vblank end line. This is the corrective offset we
				// need to subtract from read out scanline position to get true scanline position.
				// Hardware registers "scanline position" measures positive distance from vsync start line (== "scanline 0").
				// The low-word likely encodes hsyncstart to hblank end length in pixels, but we're not interested in that,
				// so we shift and mask it out:
				#if PSYCH_SYSTEM != PSYCH_WINDOWS
				vblbias = (int) ((PsychOSKDReadRegister(crtcid, 0x610000 + 0xa00 + 0xe8 + (crtcid * 0x540), NULL) >> 16) & 0xFFFF);

				// DISPLAY_TOTAL: Encodes VTOTAL in high-word, HTOTAL in low-word. Get the VTOTAL in high word:
				vbltotal = (int) ((PsychOSKDReadRegister(crtcid, 0x610000 + 0xa00 + 0xf8 + (crtcid * 0x540), NULL) >> 16) & 0xFFFF);
				#endif
			} else {
				// Auto-Detection. Read values directly from pre-NV-50 class hardware:
				// We only get VTOTAL and assume a bias value of zero, which seems to be
				// the case according to measurements on NV-40 and NV-30 gpu's:
				#if PSYCH_SYSTEM != PSYCH_WINDOWS
				vblbias = 0;

				// FP_TOTAL 0x804 relative to PRAMDAC base 0x680000 with stride 0x2000: Encodes VTOTAL in low-word:
				vbltotal = (int) ((PsychOSKDReadRegister(crtcid, 0x680000 + 0x804 + ((crtcid > 0) ? 0x2000 : 0), NULL)) & 0xFFFF) + 1;
				#endif
			}
		}

		if ((strstr((char*) glGetString(GL_VENDOR), "INTEL") || strstr((char*) glGetString(GL_VENDOR), "Intel") ||
			strstr((char*) glGetString(GL_RENDERER), "INTEL") || strstr((char*) glGetString(GL_RENDERER), "Intel")) &&
			PsychOSIsKernelDriverAvailable(screenId)) {
			#if PSYCH_SYSTEM != PSYCH_WINDOWS
			vblbias = 0;

			// VTOTAL at 0x6000C with stride 0x1000: Encodes VTOTAL in upper 16 bit word masked with 0x1fff :
			vbltotal = (int) 1 + ((PsychOSKDReadRegister(crtcid, 0x6000c + (crtcid * 0x1000), NULL) >> 16) & 0x1FFF);

			// Decode VBL_START and VBL_END for debug purposes:
			if (PsychPrefStateGet_Verbosity() > 5) {
				unsigned int vbl_start, vbl_end, vbl;
				vbl = PsychOSKDReadRegister(crtcid, 0x60010 + (crtcid * 0x1000), NULL);
				vbl_start = vbl & 0x1fff;
				vbl_end   = (vbl >> 16) & 0x1FFF;
				printf("PTB-DEBUG: Screen %i [head %i]: vbl_start = %i  vbl_end = %i.\n", screenId, crtcid, (int) vbl_start, (int) vbl_end);
			}
			#endif
		}
	}

	// Feedback is good:
	if (((vblbias != 0) || (vbltotal != 0)) && (PsychPrefStateGet_Verbosity() > 3)) {
		printf("PTB-INFO: Screen %i [head %i]: Applying beamposition corrective offsets: vblbias = %i, vbltotal = %i.\n", screenId, crtcid, vblbias, vbltotal);
	}

	// Assign:
	screenBeampositionBias[screenId] = vblbias;
	screenBeampositionVTotal[screenId] = vbltotal;
}
