/*
    PsychToolbox3/Source/Common/Screen/PsychGraphicsHardwareHALSupport.c

    AUTHORS:

        mario.kleiner.de@gmail.com  mk

    PLATFORMS:

        All. However with dependencies on OS specific glue-layers which are mostly Linux/OSX for now...

    HISTORY:

    01/12/2008  mk      Created.

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

// Array with register offsets of the CRTCs used by AMD/ATI gpus.
// Initialized by OS specific screen glue, accessed from different locations:
unsigned int crtcoff[kPsychMaxPossibleCrtcs];

// Maps screenid's to Graphics hardware pipelines: Used to choose pipeline for beampos-queries and similar
// GPU crtc specific stuff. Each screen can have up to kPsychMaxPossibleCrtcs assigned. Slot 0 contains the
// primary crtc, used for beamposition timestamping, framerate queries etc. A -1 value in a slot terminates
// the sequence of assigned crtc's.
static int  displayScreensToCrtcIds[kPsychMaxPossibleDisplays][kPsychMaxPossibleCrtcs];
static int  displayScreensToPipes[kPsychMaxPossibleDisplays][kPsychMaxPossibleCrtcs];
static psych_bool displayScreensToCrtcIdsUserOverride = FALSE;
static psych_bool displayScreensToPipesAutoDetected = FALSE;

// Corrective values for beamposition queries to correct for any constant and systematic offsets in
// the scanline positions returned by low-level code:
static int  screenBeampositionBias[kPsychMaxPossibleDisplays];
static int  screenBeampositionVTotal[kPsychMaxPossibleDisplays];

/* PsychSynchronizeDisplayScreens() -- (Try to) synchronize display refresh cycles of multiple displays
 *
 * This tries whatever method is available/appropriate/or requested to synchronize the video refresh
 * cycles of multiple graphics cards physical display heads -- corresponding to PTB logical Screens.
 *
 * The method may or may not be supported on a specific OS/gfx-card combo. It will return a PsychError_unimplemented
 * if it can't do what core wants.
 *
 * numScreens   =   Ptr to the number of display screens to sync. If numScreens>0, all screens with the screenIds stored
 *                  in the integer array 'screenIds' will be synched. If numScreens == 0, PTB will try to sync all
 *                  available screens in the system. On return, the location will contain the count of synced screens.
 *
 * screenIds    =   Either a list with 'numScreens' screenIds for the screens to sync, or NULL if numScreens == 0.
 *
 * residuals    =   List with 'numScreens' (on return) values indicating the residual sync error wrt. to the first
 *                  screen (the reference). Ideally all items should contain zero for perfect sync on return.
 *
 * syncMethod   =   Numeric Id for the sync method to use: 0 = Don't care, whatever is appropriate. 1 = Only hard
 *                  sync, which is fast and reliable if supported. 2 = Soft sync via drift-syncing. More to come...
 *
 * syncTimeOut  =   If some non-immediate method is requested/chosen, it should give up after syncTimeOut seconds if
 *                  it doesn't manage to bring the displays in sync in that timeframe.
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
    if (syncMethod > 2) PsychErrorExitMsg(PsychError_internal, "Invalid syncMethod argument passed!");
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
psych_bool PsychSetOutputDithering(PsychWindowRecordType* windowRecord, int screenId, unsigned int ditherEnable)
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
 * 'changeDithering' FALSE = Don't touch dither control, TRUE = Control dithering enable/disable if possible.
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
unsigned int PsychSetGPUIdentityPassthrough(PsychWindowRecordType* windowRecord, int screenId, psych_bool passthroughEnable, psych_bool changeDithering)
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

    #ifdef PTB_USE_WAYLAND
    // Wayland/colord specific identity setup requested?
    if (getenv("PSYCH_USE_COLORD_IDENTITYLUT")) {
        // Try to use experimental code for profiling inhibition on Wayland + colord setups.
        // Call conventional code below on failure:
        if (PsychWaylandProfilingInhibit(screenId, passthroughEnable)) {
            return(2);
        }
    }
    #endif

    // Check if kernel driver is enabled, otherwise this won't work:
    if (!PsychOSIsKernelDriverAvailable(screenId)) {
        if(PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: GPU framebuffer passthrough setup requested, but this is not supported without kernel driver.\n");
        return(0xffffffff);
    }

    // Try to enable or disable dithering on display:
    if (changeDithering) PsychSetOutputDithering(windowRecord, screenId, (passthroughEnable) ? 0 : 1);

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

/* PsychEnableNative10BitFramebuffer()  - Enable/Disable native >= 10 bpc RGB framebuffer modes.
 *
 * This function enables or disables the native high bit depth framebuffer readout modes of supported
 * graphics hardware. Currently the ATI Radeon X1000, HD2000 and later cards do allow this.
 *
 * This needs support from the Psychtoolbox kernel level support driver for low-level register reads
 * and writes to the GPU registers.
 *
 * 'windowRecord'   Is used to find the Id of the screen for which mode should be changed, as well as enable
 *                  flags to see if a change is required at all, and the OpenGL context for some specific
 *                  fixups. A value of NULL will try to apply the operation to all heads, but may only work
 *                  for *disabling* 10 bpc mode, not for enabling it -- Mostly useful for a master reset to
 *                  system default, e.g., as part of error handling or Screen shutdown handling.
 * 'enable'   True = Enable high bit depth support, False = Disable high bit depth support, reenable ARGB8888 support.
 *
 */
psych_bool PsychEnableNative10BitFramebuffer(PsychWindowRecordType* windowRecord, psych_bool enable)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
    int i, si, ei, headid, headiter, screenId;
    unsigned int lutreg, ctlreg, value, status;
    int gpuMaintype, gpuMinortype, fNumDisplayHeads;

    // Child protection:
    if (windowRecord && !PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_internal, "Invalid non-onscreen windowRecord provided!!!");

    // Either screenid from windowRecord or our special -1 "all Screens" Id:
    screenId = (windowRecord) ? windowRecord->screenNumber : -1;

    // We only support Radeon GPU's, nothing else:
    if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, &fNumDisplayHeads) ||
        (gpuMaintype != kPsychRadeon) || (gpuMinortype >= 0xffff)) {
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
            if (gpuMinortype < 40) {
                // DCE-3 and earlier:
                lutreg = (headid == 0) ? RADEON_D1GRPH_LUT_SEL : RADEON_D2GRPH_LUT_SEL;
                ctlreg = (headid == 0) ? RADEON_D1GRPH_CONTROL : RADEON_D2GRPH_CONTROL;
            }
            else {
                // DCE-4 and later:
                if (headid > fNumDisplayHeads - 1) {
                    printf("PTB-ERROR: Invalid headId %i (greater than max %i) provided for DCE-4+ display engine!\n", headid, fNumDisplayHeads - 1);
                    return(false);
                }

                lutreg = EVERGREEN_DC_LUT_10BIT_BYPASS + crtcoff[headid];
                ctlreg = EVERGREEN_GRPH_CONTROL + crtcoff[headid];
            }

            // Enable or Disable?
            if (enable) {
                // Enable:

                // Switch hardware LUT's to bypass mode:
                // We set bit 8 to enable "bypass LUT in 2101010 mode", aka bypass in all high bit depth modes:
                value = PsychOSKDReadRegister(screenId, lutreg, &status);
                if (status) {
                    printf("PTB-ERROR: Failed to enable lut bypass framebuffer mode (LUTReg read failed).\n");
                    return(false);
                }

                // Set the bypass bit:
                value = value | 0x1 << 8;

                PsychOSKDWriteRegister(screenId, lutreg, value, &status);
                if (status) {
                    printf("PTB-ERROR: Failed to enable lut bypass framebuffer mode (LUTReg write failed).\n");
                    return(false);
                }

                // Only reconfigure framebuffer scanout if this is really our true Native10bpc hack:
                // This is usually skipped on FireGL/FirePro GPU's as their drivers do it already...
                if (windowRecord->specialflags & kPsychNative10bpcFBActive) {
                    // Switch CRTC to ABGR2101010 or BGR101111 or RGBA16161616 mode:
                    value = PsychOSKDReadRegister(screenId, ctlreg, &status);
                    if (status) {
                        printf("PTB-ERROR: Failed to set high bit depth framebuffer mode (CTLReg read failed).\n");
                        return(false);
                    }

                    if (windowRecord->depth == 48) {
                        // Set bit 0 to switch from 32 bpp to 64 bpp framebuffer:
                        value = value | (1 << 0);
                    }
                    else {
                        // Set 2101010 or 101111 mode bits:
                        value = value | (((windowRecord->depth == 33) ? 7 : 1) << 8);
                    }

                    PsychOSKDWriteRegister(screenId, ctlreg, value, &status);
                    if (status) {
                        printf("PTB-ERROR: Failed to set high bit depth framebuffer mode (CTLReg write failed).\n");
                        return(false);
                    }

                    // Pipe should be in 10 bpc mode now...
                    if (PsychPrefStateGet_Verbosity() > 2) {
                        printf("PTB-INFO: System framebuffer switched to %s mode for screen %i [head %i].\n", (windowRecord->depth == 48) ? "RGBA16161616" : (windowRecord->depth == 33) ? "BGR101111" : "ARGB2101010", i, headid);
                    }
                }
            } else {
                // Disable:

                // Only reconfigure framebuffer scanout if this is really our true Native10bpc hack:
                // This is usually skipped on FireGL/FirePro GPU's as their drivers do it already...
                if (windowRecord->specialflags & kPsychNative10bpcFBActive) {
                    // Switch CRTC to ABGR8888 readout mode:
                    // We clear bits 0 and 8:10 to enable that mode:
                    value = PsychOSKDReadRegister(screenId, ctlreg, &status);
                    if (status) {
                        // This codepath gets always called in PsychCloseWindow(), so we should skip it
                        // silently if register read fails, as this is expected on MS-Windows and on all
                        // non-Radeon hardware and if kernel driver isn't loaded:
                        if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-ERROR: Failed to set 8 bit framebuffer mode (CTLReg read failed).\n");
                        return(false);
                    }
                    else if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In disable high bit depth: Readreg. ctlreg yields %lx\n", value);

                    // Clear bits 0, 8:10 :
                    value = value & ~((0x7 << 8) | 0x1);

                    PsychOSKDWriteRegister(screenId, ctlreg, value, &status);
                    if (status) {
                        printf("PTB-ERROR: Failed to set 8 bit framebuffer mode (CTLReg write failed).\n");
                        return(false);
                    }
                    else if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In disable high bit depth: ctlreg reset\n");

                    // Wait 500 msecs for GPU to settle:
                    PsychWaitIntervalSeconds(0.5);

                    // Pipe should be in 8 bpc mode now...
                    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: System framebuffer switched to standard ARGB8888 mode for screen %i [head %i].\n", i, headid);
                }

                // Switch hardware LUT's to standard mapping mode:
                // We clear bit 8 to disable "bypass LUT in 2101010 mode":
                value = PsychOSKDReadRegister(screenId, lutreg, &status);
                if (status) {
                    printf("PTB-ERROR: Failed to disable lut bypass framebuffer mode (LUTReg read failed).\n");
                    return(false);
                }
                else if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In disable lut bypass: Readreg. lutreg yields %lx\n", value);

                // Clear LUT bypass bit:
                value = value & ~(0x1 << 8);

                PsychOSKDWriteRegister(screenId, lutreg, value, &status);
                if (status) {
                    printf("PTB-ERROR: Failed to disable lut bypass framebuffer mode (LUTReg write failed).\n");
                    return(false);
                }
                else if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In disable lut bypass: lutreg reset\n");
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
 * would revert the framebuffers opmode to standard 8bpc mode and thereby kill our >= 10 bpc mode
 * setting.
 *
 * This routine *must* be called after each such problematic "End of scene" marker command like
 * glClear(). The routine does nothing if 10bpc mode is not enabled/requested for the corresponding
 * display head associated with the given onscreen window. It rewrites the control register on
 * >= 10 bpc configured windows to basically undo the unwanted change of the gfx-driver *before*
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
 *    --> "good value" persists, framebuffer stays in high bpc configuration --> All good.
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
    unsigned int ctlreg, lutreg, val1, val2;
    int gpuMaintype, gpuMinortype, fNumDisplayHeads;

    // Fixup needed? Only if 10bpc mode is supposed to be active! Early exit if not:
    if (!(windowRecord->specialflags & kPsychNative10bpcFBActive)) return;

    // Map windows screen to gfx-headid aka register subset. TODO : We'll need something better,
    // more generic, abstracted out for the future, but as a starter this will do:
    screenId = windowRecord->screenNumber;

    // We only support Radeon GPU's, nothing else:
    if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, &fNumDisplayHeads) ||
        (gpuMaintype != kPsychRadeon) || (gpuMinortype >= 0xffff)) {
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

        // Select Radeon HW registers for corresponding heads:
        if (gpuMinortype < 40) {
            // DCE-3 and earlier:
            lutreg = (headid == 0) ? RADEON_D1GRPH_LUT_SEL : RADEON_D2GRPH_LUT_SEL;
            ctlreg = (headid == 0) ? RADEON_D1GRPH_CONTROL : RADEON_D2GRPH_CONTROL;
        }
        else {
            // DCE-4 and later:
            if (headid > fNumDisplayHeads - 1) {
                printf("PTB-ERROR: Invalid headId %i (greater than max %i) provided for DCE-4+ display engine!\n", headid, fNumDisplayHeads - 1);
                return;
            }

            lutreg = EVERGREEN_DC_LUT_10BIT_BYPASS + crtcoff[headid];
            ctlreg = EVERGREEN_GRPH_CONTROL + crtcoff[headid];
        }

        // Get current state of registers at high debug levels:
        if (PsychPrefStateGet_Verbosity() > 9) {
            val1 = PsychOSKDReadRegister(screenId, lutreg, NULL);
            val2 = PsychOSKDReadRegister(screenId, ctlreg, NULL);
            printf("PTB-DEBUG: PsychFixupNative10BitFramebufferEnableAfterEndOfSceneMarker(): Screen %i, head %i: LUT = %i [%i], GRPHCONT = %i [%i]\n", screenId, headid, val1 & (0x1 << 8), val1, (val2 & (0x7 << 8)) >> 8, val2);
        }

        // One-liner read-modify-write op, which simply sets bit 8 of the register - the "10 bit LUT bypass" bit:
        PsychOSKDWriteRegister(screenId, lutreg, (0x1 << 8) | PsychOSKDReadRegister(screenId, lutreg, NULL), NULL);

        if (windowRecord->depth == 48) {
            // One-liner read-modify-write op, which simply sets bit 0 of the register - the "Enable 64 bpp mode" bit:
            PsychOSKDWriteRegister(screenId, ctlreg, (0x1) | PsychOSKDReadRegister(screenId, ctlreg, NULL), NULL);
        }
        else {
            // One-liner read-modify-write op, which simply sets bit 8 or bit 8:10 of the register - the "Enable 2101010 mode or Enable 101111 mode " bit:
            PsychOSKDWriteRegister(screenId, ctlreg, (((windowRecord->depth == 33) ? 7 : 1) << 8) | PsychOSKDReadRegister(screenId, ctlreg, NULL), NULL);
        }

        // Debug output, if wanted:
        if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychFixupNative10BitFramebufferEnableAfterEndOfSceneMarker(): Control bits set on screen %i, head %i.\n", screenId, headid);
    }
#endif

    // Done.
    return;
}

/* PsychGetCurrentGPUSurfaceAddresses() - Get current scanout surface addresses
 *
 * Tries to get current addresses of primary and secondary scanout buffers and
 * the pending status of pending pageflips if any.
 *
 * primarySurface - Pointer to 64-Bit target for frontBuffers address.
 * secondarySurface - Pointer to a potential secondary buffer, e.g., for frame-sequential stereo.
 * updatePending - TRUE if a flip has been queued and is still pending. FALSE otherwise or "don't know"
 *
 * Returns TRUE on success, FALSE if given GPU isn't supported for such queries.
 *
 */
psych_bool PsychGetCurrentGPUSurfaceAddresses(PsychWindowRecordType* windowRecord, psych_uint64* primarySurface, psych_uint64* secondarySurface, psych_bool* updatePending)
{
    #if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
        int gpuMaintype, gpuMinortype, fNumDisplayHeads;
        unsigned int updateStatus;

        // If we are called, we know that 'windowRecord' is an onscreen window.
        int screenId = windowRecord->screenNumber;
        int headid = PsychScreenToCrtcId(screenId, 0);

        // Just need to check if GPU low-level access is supported:
        if (!PsychOSIsKernelDriverAvailable(screenId)) return(FALSE);

        // We only support AMD Radeon/Fire GPU's, nothing else:
        if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, &fNumDisplayHeads) || (gpuMaintype != kPsychRadeon)) {
            return(FALSE);
        }

        // Driver is online: Read the registers, but only for primary crtc in a multi-crtc config:
        if  (gpuMinortype < 40) {
            // Pre DCE-4: AVIVO class display hardware:
            *primarySurface = (psych_uint64) PsychOSKDReadRegister(screenId, (PsychScreenToCrtcId(screenId, 0) < 1) ? RADEON_D1GRPH_PRIMARY_SURFACE_ADDRESS : RADEON_D2GRPH_PRIMARY_SURFACE_ADDRESS, NULL);
            *secondarySurface = (psych_uint64) PsychOSKDReadRegister(screenId, (PsychScreenToCrtcId(screenId, 0) < 1) ? RADEON_D1GRPH_SECONDARY_SURFACE_ADDRESS : RADEON_D2GRPH_SECONDARY_SURFACE_ADDRESS, NULL);
            updateStatus = PsychOSKDReadRegister(screenId, (PsychScreenToCrtcId(screenId, 0) < 1) ? RADEON_D1GRPH_UPDATE : RADEON_D2GRPH_UPDATE, NULL);

            *updatePending = (updateStatus & RADEON_SURFACE_UPDATE_PENDING) ? TRUE : FALSE;
        }

        if  (gpuMinortype >= 40) {
            // DCE-4 or later display hardware:
            if (headid < 0 || headid > fNumDisplayHeads - 1) {
                if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: PsychGetCurrentGPUSurfaceAddresses(): Invalid headId %i (greater than max %i) provided for DCE-4+ display engine!\n", headid, fNumDisplayHeads - 1);
                return(FALSE);
            }

            // DCE 4 and later has 64-Bit addresses split in high/low 32-bit regs:
            *primarySurface =   ((psych_uint64) PsychOSKDReadRegister(screenId, EVERGREEN_GRPH_PRIMARY_SURFACE_ADDRESS_HIGH + crtcoff[headid], NULL) << 32) +
                                ((psych_uint64) PsychOSKDReadRegister(screenId, EVERGREEN_GRPH_PRIMARY_SURFACE_ADDRESS + crtcoff[headid], NULL));
            *secondarySurface = ((psych_uint64) PsychOSKDReadRegister(screenId, EVERGREEN_GRPH_SECONDARY_SURFACE_ADDRESS_HIGH + crtcoff[headid], NULL) << 32) +
                                ((psych_uint64) PsychOSKDReadRegister(screenId, EVERGREEN_GRPH_SECONDARY_SURFACE_ADDRESS + crtcoff[headid], NULL));
            updateStatus =      PsychOSKDReadRegister(screenId, EVERGREEN_GRPH_UPDATE + crtcoff[headid], NULL);

            *updatePending = (updateStatus & EVERGREEN_GRPH_SURFACE_UPDATE_PENDING) ? TRUE : FALSE;
        }

        if (PsychPrefStateGet_Verbosity() > 14) {
            printf("PTB-DEBUG: Screen %i: Head %i: primarySurface=%p : secondarySurface=%p : updateStatus=%i\n", screenId, PsychScreenToCrtcId(screenId, 0), *primarySurface, *secondarySurface, updateStatus);
        }

        // Success:
        return(TRUE);
    #else
        // Not supported:
        return(FALSE);
    #endif
}

/* Stores content of GPU's surface address registers of the surfaces that
 * correspond to the windowRecords frontBuffers. Only called inside
 * PsychExecuteBufferSwapPrefix() immediately before triggering a double-buffer
 * swap. The values are used as reference values: If another readout of
 * these registers shows values different from the ones stored preflip,
 * then that is a certain indicator that bufferswap via pageflip has happened
 * or will happen.
 */
void PsychStoreGPUSurfaceAddresses(PsychWindowRecordType* windowRecord)
{
    psych_bool updatePending;
    PsychGetCurrentGPUSurfaceAddresses(windowRecord, &(windowRecord->gpu_preflip_Surfaces[0]), &(windowRecord->gpu_preflip_Surfaces[1]), &updatePending);
    return;
}

/* PsychIsGPUPageflipUsed() - Is a pageflip used on the GPU for buffer swaps at this moment?
 *
 * This routine compares preflip scanout addresses, as gathered via a previous PsychStoreGPUSurfaceAddresses()
 * call prior to scheduling a swap, with the current addresses and update status. It should only be called
 * after we detected bufferswap completion to check if the swap happened via pageflip and therefore our
 * completion detection and timestamping is trustworthy, or if the swap happened by some other means like
 * compositor or copyswap and therefore our results wrt. completion or timestamping are not trustworthy -
 * at least not for conventional timestamping as used on OSX, or Linux without special OS support.
 *
 * The interesting scenario is if - after detection of swap completion by our conventional standard method
 * for use with proprietary graphics drivers - the surface scanout addresses have changed and the flip is
 * confirmed finished. In this case we can be somewhat certain that we triggered the pageflip and it completed,
 * ie. the results for timestamping are trustworthy. This is indicated by return value 2. If a flip is used
 * but still pending (value 1) although our code assumes swap has completed then a pageflip was likely queued
 * by the desktop compositor and is still pending -> Timestamping not trustworthy. A value of 0 could indicate
 * copyswap or a compositor to which we sent our updated composition source surface and posted damage, but which
 * hasn't yet picked up on it or at least hasn't performed the full composition pass + queueing a pageflip.
 *
 * Ergo: For checking the trustworthiness of swap completion timestamping, the only "good" result is a
 * return value of 2, a value of 0 or 1 is considered bad for timing, a value of -1 is non-diagnostic.
 *
 * As of beginning of March 2015, this routine can only be used with some confidence on Linux for conventional
 * timestamping with the proprietary drivers and X11, as we know how X11 compositors work on Linux and what
 * to expect. Use with FOSS graphics stack or on Wayland is not needed as we have much better facilities there.
 * Additionally the PsychGetCurrentGPUSurfaceAddresses() support code is limited to AMD gpu's, so the only
 * interesting/valid use cases are Linux/X11 + proprietary AMD Catalyst driver for some clever handling, and
 * OSX + PsychtoolboxKernelDriver + AMD gpu for purely diagnostic use for manual diagnostic, not automatic
 * problem solving!
 *
 * Return values:
 * -1  == Unknown / Query unsupported.
 *  0  == No.
 *  1  == Yes, and the flip has been queued but is still pending, not finished.
 *  2  == Yes, and the flip is finished.
 *
 */
int PsychIsGPUPageflipUsed(PsychWindowRecordType* windowRecord)
{
    psych_uint64 primarySurface, secondarySurface;
    psych_bool updatePending;

    if (!PsychGetCurrentGPUSurfaceAddresses(windowRecord, &primarySurface, &secondarySurface, &updatePending)) {
        // Query not possible/supported: Return "I don't know" value -1:
        return(-1);
    }

    // Scanout addresses changed since last PsychStoreGPUSurfaceAddresses() call?
    // That would mean a pageflip was either queued or already executed, in any
    // case, pageflip is used for bufferswap:
    if (primarySurface != windowRecord->gpu_preflip_Surfaces[0] || secondarySurface != windowRecord->gpu_preflip_Surfaces[1]) {
        // Pageflip in use. Still pending (queued but not completed) or already completed?
        if (updatePending) return(1);
        return(2);
    }

    // Nope, scanout hasn't changed: Assume copyswap/blit etc.
    return(0);
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
    double deadline = *timestamp;

    // If we are called, we know that 'windowRecord' is an onscreen window.
    int screenId = windowRecord->screenNumber;

    // Retrieve display id and screen size spec that is needed later...
    PsychGetCGDisplayIDFromScreenNumber(&displayID, screenId);

    // Just need to check if GPU low-level access is supported:
    if (!PsychOSIsKernelDriverAvailable(screenId)) return(FALSE);

    // We only support AMD GPU's:
    if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, NULL) || (gpuMaintype != kPsychRadeon)) {
        return(FALSE);
    }

    // Driver is online. Enter polling loop:
    while (TRUE) {
        PsychGetAdjustedPrecisionTimerSeconds(timestamp);

        if ((PsychIsGPUPageflipUsed(windowRecord) > 0) || (*timestamp > deadline)) {
            // Abort condition: Exit loop.
            break;
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
    return((PsychIsGPUPageflipUsed(windowRecord) == 1) ? TRUE : FALSE);
#else
    // On Windows, we always return "swap happened":
    return(FALSE);
#endif
}

/* PsychGetNVidiaGPUType()
 *
 * Decodes hw register of NVidia GPU into GPU core id / chip family:
 * Returns 0 for unknown card, otherwise xx for NV_xx:
 *
 * Reference Linux nouveau-kms driver implementation in:
 * nouveau/core/engine/device/base.c: nouveau_devobj_ctor()
 */
unsigned int PsychGetNVidiaGPUType(PsychWindowRecordType* windowRecord)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
    psych_uint32 chipset, card_type;

    // Get hardware id code from gpu register:
    psych_uint32 reg0 = PsychOSKDReadRegister((windowRecord) ? windowRecord->screenNumber : 0, NV03_PMC_BOOT_0, NULL);

    // Special case: Register read delivers "hardware powered down" value, because this
    // GPU is offline, e.g., in a hybrid-graphics setup. Return a special code to signal
    // to caller that it doesn't need to bother anymore with this GPU:
    if (reg0 == 0xffffffff) return(0xffffffff);

    /* We're dealing with >=NV10 */
    if ((reg0 & 0x1f000000) > 0) {
        /* Bit 28-20 contain the architecture in hex */
        chipset = (reg0 & 0x1ff00000) >> 20;
        /* NV04 or NV05 */
    } else if ((reg0 & 0xff00fff0) == 0x20004000) {
        if (reg0 & 0x00f00000)
            chipset = 0x05;
        else
            chipset = 0x04;
    } else
        chipset = 0xff;

    switch (chipset & 0x1f0) {
        case 0x000:
            // NV_04/05: RivaTNT , RivaTNT2
            card_type = 0x04;
            break;
        case 0x010:
        case 0x020:
        case 0x030:
            // NV30 or earlier: GeForce-5 / GeForceFX and earlier:
            card_type = chipset & 0xf0;
            break;
        case 0x040:
        case 0x060:
            // NV40: GeForce6/7 series:
            card_type = 0x040;
            break;
        case 0x050:
        case 0x080:
        case 0x090:
        case 0x0a0:
            // NV50: GeForce8/9/G100-G300.
            card_type = 0x050;
            break;
        case 0x0c0:
            // Fermi: GeForce G400/500 series:
            card_type = 0x0c0;
            break;
        case 0x0d0:
            // Fermi: But with 3rd gen scanout engine, but still only 2 CRTC's:
            card_type = 0x0d0;
            break;
        case 0x0e0:
        case 0x0f0:
        case 0x100:
            // Kepler: GeForce G600+ series: 3rd gen scanout engine, but now up to 4 CRTC's.
            card_type = 0x0e0;
            break;
        case 0x110:
        case 0x120:
            // Maxwell: GeForce 750+ series: 3rd gen scanout engine, up to 4 CRTC's.
            card_type = 0x110;
            break;
        case 0x130:
            // Pascal: GeForce 1000+ series: 3rd gen scanout engine, up to 4 CRTC's.
            card_type = 0x130;
            break;

        default:
            printf("PTB-DEBUG: Unknown NVidia chipset 0x%08x - Assuming latest generation.\n", reg0);
            card_type = 0x000;
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
int PsychScreenToHead(int screenId, int rankId)
{
    return(displayScreensToPipes[screenId][rankId]);
}

/* PsychSetScreenToHead() - Change mapping of a PTB screenId to GPU headId: */
void PsychSetScreenToHead(int screenId, int headId, int rankId)
{
    // Assign new mapping:
    displayScreensToPipes[screenId][rankId] = headId;
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
int PsychScreenToCrtcId(int screenId, int rankId)
{
    return(displayScreensToCrtcIds[screenId][rankId]);
}

void PsychSetScreenToCrtcId(int screenId, int crtcId, int rankId)
{
    // Assign new mapping:
    displayScreensToCrtcIds[screenId][rankId] = crtcId;

    // Mark mappings as user-defined instead of auto-detected/default-setup:
    displayScreensToCrtcIdsUserOverride = TRUE;
}

void PsychResetCrtcIdUserOverride(void)
{
    displayScreensToCrtcIdsUserOverride = FALSE;
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
    (void) numDisplays;

    displayScreensToPipesAutoDetected = FALSE;

    // Setup default identity one-to-one mapping:
    for(i = 0; i < kPsychMaxPossibleDisplays; i++) {
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
}

// Try to auto-detect screen to head mappings if possible and not yet overriden by usercode:
void PsychAutoDetectScreenToHeadMappings(int maxHeads)
{
#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
    int gpuMaintype, gpuMinortype;
    float nullTable[256];
    int screenId, outputId, headId, numEntries;
    float *redTable, *greenTable, *blueTable;

    // If user / usercode has provided manual mapping, i.e., overriden the
    // default identity mapping, then we don't do anything, but accept the
    // users choice instead. Also skip this if it has been successfully executed
    // already:
    if (displayScreensToCrtcIdsUserOverride || displayScreensToPipesAutoDetected) return;

    // nullTable is our all-zero gamma table:
    memset(&nullTable[0], 0, sizeof(nullTable));

    // Ok, iterate over all logical screens and try to update
    // their mapping:
    for (screenId = 0; screenId < PsychGetNumDisplays(); screenId++) {
        // Kernel driver for this screenId enabled? Otherwise we skip it:
        if (!PsychOSIsKernelDriverAvailable(screenId)) continue;

        // We only support AMD/ATI gpus at the moment, nothing else:
        if (!PsychGetGPUSpecs(screenId, &gpuMaintype, &gpuMinortype, NULL, NULL) ||
            (gpuMaintype != kPsychRadeon) || (gpuMinortype >= 0xffff)) {
            continue;
        }

        // Yes. Perform detection sequence:
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Trying to detect screenId to display head mapping for screenid %i ...\n", screenId);

        // Iterate over all outputs for this screen:
        for (outputId = 0; outputId < kPsychMaxPossibleCrtcs; outputId++) {
            // Abort iteration if last output reached:
            if (PsychScreenToHead(screenId, outputId) < 0)
                break;

            if (PsychPrefStateGet_Verbosity() > 3) printf(" ... probing for video output %i ...\n", outputId);

            // Retrieve current gamma table. Need to back it up internally:
            PsychReadNormalizedGammaTable(screenId, outputId, &numEntries, &redTable, &greenTable, &blueTable);

            // Now load an all-zero gamma table for that screen:
            PsychLoadNormalizedGammaTable(screenId, outputId, 256, nullTable, nullTable, nullTable);

            // Wait for 100 msecs, so the gamma table has actually settled (e.g., if its update was
            // delayed until next vblank on a >= 20 Hz display):
            PsychYieldIntervalSeconds(0.100);

            // Check all display heads to find the null table:
            for (headId = 0; headId < maxHeads; headId++) {
                if (PsychOSKDGetLUTState(screenId, headId, 0) == 1) {
                    // Got it. Store mapping:
                    displayScreensToCrtcIds[screenId][outputId] = headId;

                    // Done with searching:
                    if (PsychPrefStateGet_Verbosity() > 3) printf(" found GPU hardware headId %i for output %i.\n", headId, outputId);
                    break;
                }
            }

            if ((headId == maxHeads) && (PsychPrefStateGet_Verbosity() > 3))
                printf(" detection failed for output %i. Sticking with heuristics result.\n", outputId);

            // Now restore original gamma table for that screen:
            PsychLoadNormalizedGammaTable(screenId, outputId, numEntries, redTable, greenTable, blueTable);

            // Wait for 100 msecs, so the gamma table has actually settled (e.g., if its update was
            // delayed until next vblank on a >= 20 Hz display):
            PsychYieldIntervalSeconds(0.100);

            if (PsychPrefStateGet_Verbosity() > 3) printf(" Done.\n");

            // Next output for this screen.
        }
        // Next screen.
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
    int gpuMaintype;
    int crtcid = PsychScreenToCrtcId(screenId, 0);

    // Auto-Detection of correct values requested?
    if (((unsigned int) vblbias == 0xffffffff) && ((unsigned int) vbltotal == 0xffffffff)) {
        // First set'em to neutral safe values in case we fail our auto-detect:
        vblbias  = 0;
        vbltotal = 0;

        // Get model of display gpu, which provides beamposition:
        gpuMaintype = kPsychUnknown;
        PsychGetGPUSpecs(screenId, &gpuMaintype, NULL, NULL, NULL);

        // Can do this on NVidia GPU's >= NV-50 if low-level access (PTB kernel driver or equivalent) is enabled:
        if ((gpuMaintype == kPsychGeForce) && PsychOSIsKernelDriverAvailable(screenId)) {
            // Need to read different regs, depending on GPU generation:
            if ((PsychGetNVidiaGPUType(NULL) >= 0x0d0) || (PsychGetNVidiaGPUType(NULL) == 0x0)) {
                // Auto-Detection. Read values directly from NV-D0 / E0-"Kepler" class and later hardware:
                //
                #if PSYCH_SYSTEM != PSYCH_WINDOWS
                // VBLANKE end line of vertical blank - smaller than VBLANKS. Subtract VBLANKE + 1 to normalize to "scanline zero is start of active scanout":
                vblbias = (int) ((PsychOSKDReadRegister(crtcid, 0x64041c + (crtcid * 0x300), NULL) >> 16) & 0xFFFF) + 1;

                // DISPLAY_TOTAL: Encodes VTOTAL in high-word, HTOTAL in low-word. Get the VTOTAL in high word:
                vbltotal = (int) ((PsychOSKDReadRegister(crtcid, 0x640414 + (crtcid * 0x300), NULL) >> 16) & 0xFFFF);

                // Decode VBL_START and VBL_END and VACTIVE for debug purposes:
                if (PsychPrefStateGet_Verbosity() > 5) {
                    unsigned int vbl_start, vbl_end, vactive;
                    vbl_start = (int) ((PsychOSKDReadRegister(crtcid, 0x640420 + (crtcid * 0x300), NULL) >> 16) & 0xFFFF);
                    vbl_end   = (int) ((PsychOSKDReadRegister(crtcid, 0x64041c + (crtcid * 0x300), NULL) >> 16) & 0xFFFF);
                    vactive   = (int) ((PsychOSKDReadRegister(crtcid, 0x640414 + (crtcid * 0x300), NULL) >> 16) & 0xFFFF);
                    printf("PTB-DEBUG: Screen %i [head %i]: vbl_start = %i  vbl_end = %i  vactive = %i.\n", screenId, crtcid, (int) vbl_start, (int) vbl_end, (int) vactive);
                }
                #endif
            }
            else if (PsychGetNVidiaGPUType(NULL) >= 0x50) {
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

                // Decode VBL_START and VBL_END and VACTIVE for debug purposes:
                if (PsychPrefStateGet_Verbosity() > 5) {
                    unsigned int vbl_start, vbl_end, vactive;
                    vbl_start = (int) ((PsychOSKDReadRegister(crtcid, 0x610af4 + (crtcid * 0x540), NULL) >> 16) & 0xFFFF);
                    vbl_end   = (int) ((PsychOSKDReadRegister(crtcid, 0x610aec + (crtcid * 0x540), NULL) >> 16) & 0xFFFF);
                    vactive   = (int) ((PsychOSKDReadRegister(crtcid, 0x610afc + (crtcid * 0x540), NULL) >> 16) & 0xFFFF);
                    printf("PTB-DEBUG: Screen %i [head %i]: vbl_start = %i  vbl_end = %i  vactive = %i.\n", screenId, crtcid, (int) vbl_start, (int) vbl_end, (int) vactive);
                }
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

        if ((gpuMaintype == kPsychIntelIGP) && PsychOSIsKernelDriverAvailable(screenId)) {
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
