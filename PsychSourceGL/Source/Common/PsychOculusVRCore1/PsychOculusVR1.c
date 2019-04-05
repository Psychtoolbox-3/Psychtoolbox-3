/*
 * PsychToolbox/Source/Common/PsychOculusVRCore1/PsychOculusVR1.c
 *
 * PROJECTS: PsychOculusVRCore1 only.
 *
 * AUTHORS:
 *
 * mario.kleiner.de@gmail.com   mk
 *
 * PLATFORMS:   All which are supported by the Oculus SDK/Runtime V 1.16+
 *
 * DESCRIPTION:
 *
 * A Psychtoolbox driver for the Oculus VR virtual reality
 * head sets, using the OculusVR 1.16 SDK and runtime and later.
 *
 * Copyright (c) 2018 Mario Kleiner. Licensed under the MIT license:
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "PsychOculusVR1.h"
#include "PsychTimeGlue.h"

// Need this _USE_MATH_DEFINES so MSVC 2010 knows M_PI
#define _USE_MATH_DEFINES
#include <math.h>

// Includes from Oculus SDK 1:
#include "OVR_CAPI.h"
#include "OVR_CAPI_GL.h"
#include "Extras/OVR_CAPI_Util.h"

// Number of maximum simultaneously open VR devices:
#define MAX_PSYCH_OCULUS_DEVS 10
#define MAX_SYNOPSIS_STRINGS 50

//declare variables local to this file.
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

// Our device record:
typedef struct PsychOculusDevice {
    psych_bool          closing;
    psych_bool          multiThreaded;
    psych_thread        presenterThread;
    psych_mutex         presenterLock;
    psych_condition     presentedSignal;
    ovrSession          hmd;
    ovrHmdDesc          hmdDesc;
    ovrTextureSwapChain textureSwapChain[2];
    int                 textureSwapChainLength;
    ovrMirrorTexture    mirrorTexture;
    psych_bool          isStereo;
    int                 textureWidth;
    int                 textureHeight;
    psych_bool          isTracking;
    ovrSizei            texSize[2];
    ovrFovPort          ofov[2];
    ovrEyeRenderDesc    eyeRenderDesc[2];
    uint32_t            frameIndex;
    int                 commitFrameIndex;
    int                 needSubmit;
    ovrPosef            outEyePoses[2];
    double              frameDuration;
    double              sensorSampleTime;
    double              lastPresentExecTime;
    double              scheduledPresentExecTime;
    double              VRtimeoutSecs;
    ovrPerfStats        perfStats;
} PsychOculusDevice;

PsychOculusDevice oculusdevices[MAX_PSYCH_OCULUS_DEVS];
static int available_devices = 0;
static unsigned int devicecount = 0;
static int verbosity = 3;
static psych_bool initialized = FALSE;
static ovrErrorInfo errorInfo;

static void* PresenterThreadMain(void* psychOculusDeviceToCast);
static double PresentExecute(PsychOculusDevice *oculus, psych_bool commitTextures, psych_bool onlyCommit);

void InitializeSynopsis(void)
{
    int i = 0;
    const char **synopsis = synopsisSYNOPSIS;

    synopsis[i++] = "PsychOculusVRCore1 - A Psychtoolbox driver for Oculus VR hardware.\n";
    synopsis[i++] = "This driver allows to control devices supported by the Oculus runtime V1.16 and higher.\n";
    synopsis[i++] = "Copyright (c) 2018, 2019 Mario Kleiner.\n";
    synopsis[i++] = "The PsychOculusVRCore1 driver is licensed to you under the terms of the MIT license, with the following restriction:\n";
    synopsis[i++] = "Uses the Oculus SDK which is Copyright © Facebook Technologies, LLC and its affiliates. All rights reserved.\n";
    synopsis[i++] = "See 'help License.txt' in the Psychtoolbox root folder for more details.\n";
    synopsis[i++] = "\n";
    synopsis[i++] = "Usage:";
    synopsis[i++] = "\n";
    synopsis[i++] = "Functions used by regular user scripts (mostly indirectly via PsychVRHMD() or PsychOculusVR1()):\n";
    synopsis[i++] = "\n";
    synopsis[i++] = "oldVerbosity = PsychOculusVRCore1('Verbosity' [, verbosity]);";
    synopsis[i++] = "numHMDs = PsychOculusVRCore1('GetCount');";
    synopsis[i++] = "[oculusPtr, modelName, resolutionX, resolutionY, refreshHz, controllerTypes] = PsychOculusVRCore1('Open' [, deviceIndex=0][, multiThreaded=0]);";
    synopsis[i++] = "PsychOculusVRCore1('Close' [, oculusPtr]);";
    synopsis[i++] = "PsychOculusVRCore1('SetHUDState', oculusPtr , mode);";
    synopsis[i++] = "[isVisible, playboundsxyz, outerboundsxyz] = PsychOculusVRCore1('VRAreaBoundary', oculusPtr [, requestVisible]);";
    synopsis[i++] = "[isTriggering, closestDistance, closestPointxyz, surfaceNormal] = PsychOculusVRCore1('TestVRBoundary', oculusPtr, trackedDeviceType, boundaryType);";
    synopsis[i++] = "[isTriggering, closestDistance, closestPointxyz, surfaceNormal] = PsychOculusVRCore1('TestVRBoundaryPoint', oculusPtr, pointxyz, boundaryType);";
    synopsis[i++] = "success = PsychOculusVRCore1('RecenterTrackingOrigin', oculusPtr);";
    synopsis[i++] = "oldType = PsychOculusVRCore1('TrackingOriginType', oculusPtr [, newType]);";
    synopsis[i++] = "PsychOculusVRCore1('Start', oculusPtr);";
    synopsis[i++] = "PsychOculusVRCore1('Stop', oculusPtr);";
    synopsis[i++] = "[state, touch] = PsychOculusVRCore1('GetTrackingState', oculusPtr [, predictionTime=nextFrame]);";
    synopsis[i++] = "input = PsychOculusVRCore1('GetInputState', oculusPtr, controllerType);";
    synopsis[i++] = "pulseEndTime = PsychOculusVRCore1('HapticPulse', oculusPtr, controllerType [, duration=2.5][, freq=1.0][, amplitude=1.0]);";
    synopsis[i++] = "[projL, projR] = PsychOculusVRCore1('GetStaticRenderParameters', oculusPtr [, clipNear=0.01][, clipFar=10000.0]);";
    synopsis[i++] = "[eyePoseL, eyePoseR, tracked, frameTiming] = PsychOculusVRCore1('StartRender', oculusPtr);";
    synopsis[i++] = "[eyePose, eyeIndex] = PsychOculusVRCore1('GetEyePose', oculusPtr, renderPass);";
    synopsis[i++] = "[adaptiveGpuPerformanceScale, frameStats, anyFrameStatsDropped, aswIsAvailable] = PsychOculusVRCore1('GetPerformanceStats', oculusPtr);";
    synopsis[i++] = "oldValues = PsychOculusVRCore1('FloatsProperty', oculusPtr, propertyName [, newValues]);";
    synopsis[i++] = "oldValues = PsychOculusVRCore1('FloatProperty', oculusPtr, propertyName [, newValue]);";
    synopsis[i++] = "oldString = PsychOculusVRCore1('StringProperty', oculusPtr, propertyName [, defaultString][, newString]);";
    synopsis[i++] = "\n";
    synopsis[i++] = "Functions usually only used internally by Psychtoolbox:\n";
    synopsis[i++] = "\n";
    synopsis[i++] = "[width, height, fovPort] = PsychOculusVRCore1('GetFovTextureSize', oculusPtr, eye [, fov=[HMDRecommended]][, pixelsPerDisplay=1]);";
    synopsis[i++] = "[width, height, numTextures] = PsychOculusVRCore1('CreateRenderTextureChain', oculusPtr, eye, width, height, floatFormat);";
    synopsis[i++] = "texObjectHandle = PsychOculusVRCore1('GetNextTextureHandle', oculusPtr, eye);";
    synopsis[i++] = "texObjectHandle = PsychOculusVRCore1('CreateMirrorTexture', oculusPtr, width, height);";
    synopsis[i++] = "[width, height, viewPx, viewPy, viewPw, viewPh, pptax, pptay, hmdShiftx, hmdShifty, hmdShiftz] = PsychOculusVRCore1('GetUndistortionParameters', oculusPtr, eye [, inputWidth][, inputHeight][, fov]);";
    synopsis[i++] = "trackers = PsychOculusVRCore1('GetTrackersState', oculusPtr);";
    synopsis[i++] = "PsychOculusVRCore1('EndFrameRender', oculusPtr, targetPresentTime);";
    synopsis[i++] = "[frameTiming, tPredictedOnset, referenceFrameIndex] = PsychOculusVRCore1('PresentFrame', oculusPtr [, doTimestamp=0][, when=0]);";
    synopsis[i++] = NULL;  //this tells PsychOculusVRDisplaySynopsis where to stop

    if (i > MAX_SYNOPSIS_STRINGS) {
        PrintfExit("%s: Increase dimension of synopsis[] from %ld to at least %ld and recompile.", __FILE__, (long) MAX_SYNOPSIS_STRINGS, (long) i);
    }
}

PsychError PsychOculusVR1DisplaySynopsis(void)
{
    int i;

    for (i = 0; synopsisSYNOPSIS[i] != NULL; i++)
        printf("%s\n",synopsisSYNOPSIS[i]);

    return(PsychError_none);
}

static double deg2rad(double deg)
{
    return deg / 180.0 * M_PI;
}

static double rad2deg(double rad)
{
    return rad / M_PI * 180.0;
}

PsychOculusDevice* PsychGetOculus(int handle, psych_bool dontfail)
{
    if (handle < 1 || handle > MAX_PSYCH_OCULUS_DEVS || oculusdevices[handle-1].hmd == NULL) {
        if (!dontfail) {
            printf("PTB-ERROR: Invalid Oculus device handle %i passed. No such device open.\n", handle);
            PsychErrorExitMsg(PsychError_user, "Invalid Oculus handle.");
        }

        return(NULL);
    }

    return(&(oculusdevices[handle-1]));
}

static void OVR_CDECL PsychOculusLogCB(uintptr_t userData, int level, const char* message)
{
    (void) userData;

    // Can't print from thread other than main thread at least in
    // recent Matlab versions like R2018a, so suppress debug logging
    // when called from separate VR runtime thread:
    if (!PsychIsMasterThread())
        return;

    if ((level == ovrLogLevel_Error && verbosity > 0) ||
        (level == ovrLogLevel_Info && verbosity > 2)  ||
        (level == ovrLogLevel_Debug && verbosity > 4)) {
            printf("PsychOculusVRCore1-Runtime: %s\n", message);
    }
}

void PsychOculusVRCheckInit(psych_bool dontfail)
{
    char clientName[40 * 10];
    ovrResult result;
    ovrInitParams iparms;
    memset(&iparms, 0, sizeof(iparms));
    iparms.Flags = ovrInit_RequestVersion | ovrInit_Debug; // Use debug libraries. TODO: Remove for final release!
    iparms.RequestedMinorVersion = OVR_MINOR_VERSION;
    iparms.LogCallback = PsychOculusLogCB;
    iparms.UserData = 0;  // Userdata pointer, currently NULL.
    iparms.ConnectionTimeoutMS = 0; // Default timeout.

    // Already initialized? No op then.
    if (initialized)
        return;

    // Initialize Oculus VR runtime with default parameters:
    result = ovr_Initialize(&iparms);
    if (OVR_SUCCESS(result)) {
        ovrDetectResult dresult;

        if (verbosity >= 3) printf("PsychOculusVRCore1-INFO: Oculus VR runtime (version '%s') initialized.\n", ovr_GetVersionString());

        // Setup a unique identifying name for us:
        sprintf(clientName, "EngineName: PsychOculusVR1\nEngineVersion: %s\nEnginePluginName: Date %s\n", PsychGetVersionString(), PsychGetBuildDate());
        ovr_IdentifyClient(clientName);

        // Poll for existence of a HMD (detection timeout = 0 msecs == poll only).
        // The 1.16 SDK still seems to be unable to enumerate more than 1 HMD,
        // so we still can't count how many there are:
        dresult = ovr_Detect(0);
        available_devices = (dresult.IsOculusHMDConnected) ? 1 : 0;
        if (!dresult.IsOculusServiceRunning) {
            available_devices = -1;
            if (verbosity >= 2) printf("PsychOculusVRCore1-WARNING: Could not connect to Oculus VR server process yet. Did you forget to start it?\n");
        }

        if (verbosity >= 3) printf("PsychOculusVRCore1-INFO: At startup there are %i Oculus HMDs available.\n", available_devices);
        initialized = TRUE;
    }
    else {
        if (!dontfail) {
            printf("PsychOculusVRCore1-ERROR: ovr_Initialize failed: %i\n", result);
            PsychErrorExitMsg(PsychError_system, "PsychOculusVRCore1-ERROR: Initialization of VR runtime failed. Driver disabled!");
        }
    }
}

void PsychOculusStop(int handle)
{
    int rc;
    PsychOculusDevice* oculus;
    oculus = PsychGetOculus(handle, TRUE);
    if (NULL == oculus || !oculus->isTracking) return;

    if (verbosity >= 4) printf("PsychOculusVRCore1-INFO: Tracking stopped on device with handle %i.\n", handle);

    PsychLockMutex(&(oculus->presenterLock));

    // Mark tracking as stopped:
    oculus->isTracking = FALSE;

    // Initialize eye poses to an identity mapping, where both eyes
    // just stare straight ahead from the origin. All values need to
    // be zero - provided by memset zero init below - except
    // for the w component of the Orientation quaternion.
    // This provides proper head-locking of frames submitted to the VR
    // compositor if rendering without head tracking is in use, e.g., for
    // standard 2D stereo rendering, movie playback etc., or for (ab)use
    // of the HMD as a "strapped onto the head" mono- or stereo-monitor:
    memset(oculus->outEyePoses, 0, 2 * sizeof(oculus->outEyePoses[0]));
    oculus->outEyePoses[0].Orientation.w = 1;
    oculus->outEyePoses[1].Orientation.w = 1;

    // Need to start presenterThread if this is the first Present operation:
    if ((oculus->multiThreaded) && (oculus->presenterThread == (psych_thread) NULL)) {
        // Create and startup thread:
        if ((rc = PsychCreateThread(&(oculus->presenterThread), NULL, PresenterThreadMain, (void*) oculus))) {
            PsychUnlockMutex(&(oculus->presenterLock));
            oculus->presenterThread = (psych_thread) NULL;
            printf("PsychOculusVRCore1-ERROR: Could not create internal presenterThread  [%s].\n", strerror(rc));
            PsychErrorExitMsg(PsychError_system, "Insufficient system resources for thread creation as part of VR compositor init!");
        }
    }

    PsychUnlockMutex(&(oculus->presenterLock));

    return;
}

void PsychOculusClose(int handle)
{
    int rc;
    PsychOculusDevice* oculus;
    oculus = PsychGetOculus(handle, TRUE);
    if (NULL == oculus) return;

    // presenterThread shutdown: Ask thread to terminate, wait for thread termination, cleanup and release the thread:
    PsychLockMutex(&(oculus->presenterLock));
    oculus->closing = TRUE;
    PsychUnlockMutex(&(oculus->presenterLock));

    if (oculus->presenterThread != (psych_thread) NULL) {
        if (verbosity > 5)
            printf("PTB-DEBUG: Waiting (join()ing) for helper thread of HMD %p to finish up. If this doesn't happen quickly, you'll have to kill Octave...\n", oculus);

        PsychDeleteThread(&(oculus->presenterThread));
    }

    // Ok, thread is dead. Mark it as such:
    oculus->presenterThread = (psych_thread) NULL;
    oculus->closing = FALSE;

    // PresentExecute() a last time on the main thread:
    PresentExecute(oculus, FALSE, FALSE);

    // Destroy/Release texture swap chain to compositor:
    if (oculus->textureSwapChain[0]) {
        ovr_DestroyTextureSwapChain(oculus->hmd, oculus->textureSwapChain[0]);
        oculus->textureSwapChain[0] = NULL;
    }

    if (oculus->isStereo && oculus->textureSwapChain[1]) {
        ovr_DestroyTextureSwapChain(oculus->hmd, oculus->textureSwapChain[1]);
        oculus->textureSwapChain[1] = NULL;
    }

    // Destroy mirror texture if any:
    if (oculus->mirrorTexture) {
        ovr_DestroyMirrorTexture(oculus->hmd, oculus->mirrorTexture);
        oculus->mirrorTexture = NULL;
    }

    // Close the HMD:
    ovr_Destroy(oculus->hmd);
    oculus->hmd = NULL;

    if ((rc = PsychDestroyMutex(&(oculus->presenterLock)))) {
        printf("PsychOculusVRCore1-WARNING: In PsychOculusClose(): Could not destroy presenterLock mutex lock [%s].\n", strerror(rc));
        printf("PsychOculusVRCore1-WARNING: This will cause resource leakage. Maybe you should better exit and restart Octave?");
    }

    if ((rc=PsychDestroyCondition(&(oculus->presentedSignal)))) {
        printf("PsychOculusVRCore1-WARNING: In PsychOculusClose(): Could not destroy presentedSignal condition variable [%s].\n", strerror(rc));
        printf("PsychOculusVRCore1-WARNING: This will cause resource leakage. Maybe you should better exit and restart Octave?");
    }

    if (verbosity >= 4) printf("PsychOculusVRCore1-INFO: Closed Oculus HMD with handle %i.\n", handle);

    // Done with this device:
    devicecount--;
}

void PsychOculusVR1Init(void) {
    int handle;

    for (handle = 0 ; handle < MAX_PSYCH_OCULUS_DEVS; handle++)
        oculusdevices[handle].hmd = NULL;

    available_devices = 0;
    devicecount = 0;
    initialized = FALSE;
}

PsychError PsychOculusVR1ShutDown(void) {
    int handle;

    if (initialized) {
        for (handle = 0 ; handle < MAX_PSYCH_OCULUS_DEVS; handle++)
            PsychOculusClose(handle);

        // Shutdown runtime:
        ovr_Shutdown();

        if (verbosity >= 4) printf("PsychOculusVRCore1-INFO: Oculus VR runtime shutdown complete.\n");
    }
    initialized = FALSE;

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1Verbosity(void)
{
    static char useString[] = "oldVerbosity = PsychOculusVRCore1('Verbosity' [, verbosity]);";
    static char synopsisString[] = "Returns and optionally sets level of 'verbosity' for driver debug output.\n"
                                   "'verbosity' = New level of verbosity: 0 = Silent, 1 = Errors only, 2 = Warnings, 3 = Info, 4 = Debug.\n";
    static char seeAlsoString[] = "";

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));

    // Return optional old verbosity setting:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, verbosity);

    // Get optional new verbosity setting:
    PsychCopyInIntegerArg(1, kPsychArgOptional, &verbosity);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetCount(void)
{
    static char useString[] = "numHMDs = PsychOculusVRCore1('GetCount');";
    static char synopsisString[] =  "Returns count of currently connected HMDs.\n"
                                    "Returns -1 if the runtime or server couldn't get initialized.\n";
    static char seeAlsoString[] = "Open";
    ovrDetectResult dresult;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if( PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(0));

    dresult = ovr_Detect(0);
    available_devices = (dresult.IsOculusHMDConnected) ? 1 : 0;
    if (!dresult.IsOculusServiceRunning) {
        available_devices = -1;
        if (verbosity >= 2) printf("PsychOculusVRCore1-WARNING: Could not connect to Oculus VR server process yet. Did you forget to start it?\n");
    }

    PsychCopyOutDoubleArg(1, kPsychArgOptional, available_devices);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1Open(void)
{
    static char useString[] = "[oculusPtr, modelName, resolutionX, resolutionY, refreshHz, controllerTypes]] = PsychOculusVRCore1('Open' [, deviceIndex=0][, multiThreaded=0]);";
    //                          1          2          3            4            5          6                                                1                2
    static char synopsisString[] =
        "Open connection to Oculus VR HMD, return a 'oculusPtr' handle to it.\n\n"
        "The call tries to open the HMD with index 'deviceIndex', or the "
        "first detected HMD, if 'deviceIndex' is omitted. Please note that currently "
        "only one single HMD is supported, so this "
        "'deviceIndex' is redundant at the moment, given that zero is the only valid value.\n"
        "'multiThreaded' if provided as non-zero value, will use an asynchronous presenter thread "
        "to improve stimulus scheduling. Highly experimental: Does not work in many cases!\n"
        "The returned handle can be passed to the other subfunctions to operate the device.\n"
        "A second return argument 'modelName' returns the model name string of the Oculus device.\n"
        "'resolutionX' and 'resolutionY' return the HMD display panels horizontal and vertical resolution.\n"
        "'refreshHz' is the nominal refresh rate of the display in Hz.\n"
        "'controllerTypes' A bit mask of OVR.ControllerType_XXX flags describing the currently "
        "connected input controllers.\n";

    static char seeAlsoString[] = "GetCount Close";

    ovrDetectResult result;
    PsychOculusDevice* oculus;
    int deviceIndex = 0;
    int handle = 0;
    int multiThreaded = 0;
    int rc;
    unsigned int oldCaps;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(6));
    PsychErrorExit(PsychCapNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Find a free device slot:
    for (handle = 0; (handle < MAX_PSYCH_OCULUS_DEVS) && oculusdevices[handle].hmd; handle++);
    if ((handle >= MAX_PSYCH_OCULUS_DEVS) || oculusdevices[handle].hmd) PsychErrorExitMsg(PsychError_internal, "Maximum number of simultaneously open Oculus VR devices reached.");

    // Get optional Oculus device index:
    PsychCopyInIntegerArg(1, kPsychArgOptional, &deviceIndex);

    // Don't support anything than a single "default" OculusVR Rift yet - A limitation of the current SDK:
    if (deviceIndex < 0)
        PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' provided. Must be greater or equal to zero!");

    // Get optional multiThreaded arg:
    PsychCopyInIntegerArg(2, kPsychArgOptional, &multiThreaded);
    if ((multiThreaded != 0) && (multiThreaded != 1))
        PsychErrorExitMsg(PsychError_user, "Invalid 'multiThreaded' flag provided. Must be 0 or 1!");

    result = ovr_Detect(0);
    available_devices = (result.IsOculusHMDConnected) ? 1 : 0;
    if (!result.IsOculusServiceRunning) {
        available_devices = -1;
        if (verbosity >= 2) printf("PsychOculusVRCore1-WARNING: Could not connect to Oculus VR server process yet. Did you forget to start it?\n");
    }

    if (deviceIndex >= available_devices) {
        if (verbosity >= 0) printf("PsychOculusVRCore1-ERROR: Invalid deviceIndex %i >= number of available HMDs %i.\n", deviceIndex, available_devices);
        PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' provided. Not enough HMDs available!");
    }

    // Zero init device structure:
    memset(&oculusdevices[handle], 0, sizeof(PsychOculusDevice));

    oculus = &oculusdevices[handle];

    // Try to open real or emulated HMD with deviceIndex:
    if (deviceIndex >= 0) {
        // The real thing:
        ovrGraphicsLuid Luid;
        ovrResult rc = ovr_Create(&(oculusdevices[handle].hmd), &Luid);
        if (OVR_FAILURE(rc) || (NULL == oculusdevices[handle].hmd)) {
            if (verbosity >= 0) {
                printf("PsychOculusVRCore1-ERROR: Failed to connect to Oculus Rift with deviceIndex %i. This could mean that the device\n", deviceIndex);
                printf("PsychOculusVRCore1-ERROR: is already in use by another application or driver.\n");
            }
            PsychErrorExitMsg(PsychError_user, "Could not connect to Rift device with given 'deviceIndex'! [ovr_Create() failed]");
        }
        else if (verbosity >= 3) {
            printf("PsychOculusVRCore1-INFO: Opened Oculus Rift with deviceIndex %i as handle %i.\n", deviceIndex, handle + 1);
        }
    }

    // Query current enabled caps:
    oculus->hmdDesc = ovr_GetHmdDesc(oculus->hmd);
    oldCaps = oculus->hmdDesc.AvailableHmdCaps;

    // Stats for nerds:
    if (verbosity >= 3) {
        printf("PsychOculusVRCore1-INFO: Product: %s - Manufacturer: %s - SerialNo: %s [VID: 0x%x PID: 0x%x]\n",
               oculus->hmdDesc.ProductName, oculus->hmdDesc.Manufacturer, oculus->hmdDesc.SerialNumber, (int) oculus->hmdDesc.VendorId, (int) oculus->hmdDesc.ProductId);
        printf("PsychOculusVRCore1-INFO: Firmware version: %i.%i\n", (int) oculus->hmdDesc.FirmwareMajor, (int) oculus->hmdDesc.FirmwareMinor);
        printf("PsychOculusVRCore1-INFO: Panel size in pixels w x h = %i x %i - Refresh rate = %f fps\n", oculus->hmdDesc.Resolution.w, oculus->hmdDesc.Resolution.h, oculus->hmdDesc.DisplayRefreshRate);
        printf("PsychOculusVRCore1-INFO: Caps: Debug device=%i\n", (oldCaps & ovrHmdCap_DebugDevice) ? 1 : 0);
        printf("PsychOculusVRCore1-INFO: Fields of view (Fov): [left, right, up, down] in degrees follow:\n");
        printf("PsychOculusVRCore1-INFO: Max Left eye Fov:     [%f, %f, %f, %f]\n", rad2deg(oculus->hmdDesc.MaxEyeFov[0].LeftTan), rad2deg(oculus->hmdDesc.MaxEyeFov[0].RightTan), rad2deg(oculus->hmdDesc.MaxEyeFov[0].UpTan), rad2deg(oculus->hmdDesc.MaxEyeFov[0].DownTan));
        printf("PsychOculusVRCore1-INFO: Max Right eye Fov:    [%f, %f, %f, %f]\n", rad2deg(oculus->hmdDesc.MaxEyeFov[1].LeftTan), rad2deg(oculus->hmdDesc.MaxEyeFov[1].RightTan), rad2deg(oculus->hmdDesc.MaxEyeFov[1].UpTan), rad2deg(oculus->hmdDesc.MaxEyeFov[1].DownTan));
        printf("PsychOculusVRCore1-INFO: Def Left eye Fov:     [%f, %f, %f, %f]\n", rad2deg(oculus->hmdDesc.DefaultEyeFov[0].LeftTan), rad2deg(oculus->hmdDesc.DefaultEyeFov[0].RightTan), rad2deg(oculus->hmdDesc.DefaultEyeFov[0].UpTan), rad2deg(oculus->hmdDesc.DefaultEyeFov[0].DownTan));
        printf("PsychOculusVRCore1-INFO: Def Right eye Fov:    [%f, %f, %f, %f]\n", rad2deg(oculus->hmdDesc.DefaultEyeFov[1].LeftTan), rad2deg(oculus->hmdDesc.DefaultEyeFov[1].RightTan), rad2deg(oculus->hmdDesc.DefaultEyeFov[1].UpTan), rad2deg(oculus->hmdDesc.DefaultEyeFov[1].DownTan));
        printf("PsychOculusVRCore1-INFO: ----------------------------------------------------------------------------------\n");
    }

    // Initialize eye poses to an identity mapping, where both eyes
    // just stare straight ahead from the origin. All values need to
    // be zero - provided by zero init of the oculus struct - except
    // for the w component of the Orientation quaternion.
    // This provides proper head-locking of frames submitted to the VR
    // compositor if rendering without head tracking is in use, e.g., for
    // standard 2D stereo rendering, movie playback etc., or for (ab)use
    // of the HMD as a "strapped onto the head" mono- or stereo-monitor:
    oculus->outEyePoses[0].Orientation.w = 1;
    oculus->outEyePoses[1].Orientation.w = 1;

    // Store video frame duration on the HMD:
    oculus->frameDuration = 1.0 / (double) oculus->hmdDesc.DisplayRefreshRate;

    // Assume the timeout for the compositor thinking we are unresponsive is 2 HMD frame durations:
    oculus->VRtimeoutSecs = 2 * oculus->frameDuration;

    // Assign multi-threading mode:
    oculus->multiThreaded = (psych_bool) multiThreaded;

    // Initialize the mutex lock:
    if ((rc = PsychInitMutex(&(oculus->presenterLock)))) {
        printf("PsychOculusVRCore1-ERROR: Could not create internal presenterLock mutex lock [%s].\n", strerror(rc));
        PsychErrorExitMsg(PsychError_system, "Insufficient system resources for mutex creation as part of HMD open operation!");
    }

    if ((rc = PsychInitCondition(&(oculus->presentedSignal), NULL))) {
        printf("PsychOculusVRCore1-ERROR: Could not create internal presentedSignal condition variable [%s].\n", strerror(rc));
        PsychErrorExitMsg(PsychError_system, "Insufficient system ressources for condition variable creation as part of HMD open operation!");
    }

    // Reset performance statistics for this session:
    ovr_ResetPerfStats(oculus->hmd);

    // Increment count of open devices:
    devicecount++;

    // Return device handle: We use 1-based handle indexing to make life easier for Octave/Matlab:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, handle + 1);

    // Return product name:
    PsychCopyOutCharArg(2, kPsychArgOptional, (const char*) oculus->hmdDesc.ProductName);

    // Horizontal and vertical panel resolution:
    PsychCopyOutDoubleArg(3, kPsychArgOptional, oculus->hmdDesc.Resolution.w);
    PsychCopyOutDoubleArg(4, kPsychArgOptional, oculus->hmdDesc.Resolution.h);

    // Panel refresh rate:
    PsychCopyOutDoubleArg(5, kPsychArgOptional, oculus->hmdDesc.DisplayRefreshRate);

    // Connected controllers:
    PsychCopyOutDoubleArg(6, kPsychArgOptional, (double) ovr_GetConnectedControllerTypes(oculus->hmd));

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1Close(void)
{
    static char useString[] = "PsychOculusVRCore1('Close' [, oculusPtr]);";
    //                                                      1
    static char synopsisString[] =
        "Close connection to Oculus Rift device 'oculusPtr'. Do nothing if no such device is open.\n"
        "If the optional 'oculusPtr' is omitted, then close all open devices and shutdown the driver, "
        "ie. perform the same cleanup as if 'clear PsychOculusVRCore1' would be executed.\n";
    static char seeAlsoString[] = "Open";

    int handle = -1;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(0));
    PsychErrorExit(PsychCapNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get optional device handle:
    PsychCopyInIntegerArg(1, kPsychArgOptional, &handle);

    if (handle >= 1) {
        // Close device:
        PsychOculusClose(handle);
    }
    else {
        // No handle provided: Close all devices, shutdown driver.
        PsychOculusVR1ShutDown();
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1RecenterTrackingOrigin(void)
{
    static char useString[] = "success = PsychOculusVRCore1('RecenterTrackingOrigin', oculusPtr);";
    //                         1                                                      1
    static char synopsisString[] =
        "Recenter the tracking origin for Oculus device 'oculusPtr'.\n"
        "This uses the current pose of the HMD to redefine the tracking origin.\n"
        "Returns 1 on success, 0 on failure.\n"
        "A reason for failure can be that the HMD is not roughly level, but "
        "instead is facing upward or downward, which is not allowed.\n";
    static char seeAlsoString[] = "GetTrackersState TrackingOriginType";

    int handle;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    if (OVR_FAILURE(ovr_RecenterTrackingOrigin(oculus->hmd))) {
        if (verbosity > 0) {
            ovr_GetLastErrorInfo(&errorInfo);
            printf("PsychOculusVRCore1-ERROR: Recentering the tracking origin failed! %s\n", errorInfo.ErrorString);
        }
        PsychCopyOutDoubleArg(1, kPsychArgOptional, 0);
    }
    else {
        if (verbosity > 2) {
            printf("PsychOculusVRCore1-INFO: Recentered the tracking origin for HMD.\n");
        }
        PsychCopyOutDoubleArg(1, kPsychArgOptional, 1);

        ovr_ClearShouldRecenterFlag(oculus->hmd);
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1TrackingOriginType(void)
{
    static char useString[] = "oldType = PsychOculusVRCore1('TrackingOriginType', oculusPtr [, newType]);";
    //                         1                                                  1            2
    static char synopsisString[] =
        "Specify the type of tracking origin for Oculus device 'oculusPtr'.\n\n"
        "This returns the current type of tracking origin in 'oldType'.\n\n"
        "Optionally you can specify a new tracking origin type as 'newType'. "
        "Type must be either:\n\n"
        "0 = Origin is at eye height (HMD height).\n"
        "1 = Origin is at floor height.\n\n"
        "The eye height or floor height gets defined by the system during "
        "calls to 'RecenterTrackingOrigin' and during sensor calibration in "
        "the Oculus GUI application.\n";
    static char seeAlsoString[] = "RecenterTrackingOrigin GetTrackersState";

    int handle;
    int originType;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Query and return old setting:
    originType = (int) ovr_GetTrackingOriginType(oculus->hmd);
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) originType);

    // New origin type provided?
    if (PsychCopyInIntegerArg(2, kPsychArgOptional, &originType)) {
        if (originType < 0 || originType > 1)
            PsychErrorExitMsg(PsychError_user, "Invalid 'newType' for tracking origin type specified. Must be 0 or 1.");

        if (OVR_FAILURE(ovr_SetTrackingOriginType(oculus->hmd, (ovrTrackingOrigin) originType))) {
            if (verbosity > 0) {
                ovr_GetLastErrorInfo(&errorInfo);
                printf("PsychOculusVRCore1-ERROR: Setting new tracking origin type failed! %s\n", errorInfo.ErrorString);
            }
            PsychErrorExitMsg(PsychError_user, "Setting new tracking origin type failed.");
        }
        else {
            if (verbosity > 3)
                printf("PsychOculusVRCore1-INFO: Set new tracking origin type for HMD to %i.\n", originType);
        }
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1VRAreaBoundary(void)
{
    static char useString[] = "[isVisible, playboundsxyz, outerboundsxyz] = PsychOculusVRCore1('VRAreaBoundary', oculusPtr [, requestVisible]);";
    //                          1          2              3                                                      1            2
    static char synopsisString[] =
        "Request visualization of the VR play area boundary for Oculus device 'oculusPtr' and returns its current extents.\n\n"
        "'requestVisible' 1 = Request showing the boundary area markers, 0 = Don't request showing the markers.\n"
        "The driver can not prevent the boundaries to be visualized if some external setting asks for "
        "their visibility. It can cancel its own request for visibility though via 'requestVisible' setting 0.\n\n"
        "Returns in 'isVisible' the current visibility status of the VR area boundaries.\n\n"
        "'playboundsxyz' is a 3-by-n matrix defining the play area boundaries. Each column represents "
        "the [x;y;z] coordinates of one 3D definition point. Connecting successive points by line segments "
        "defines the boundary, as projected onto the floor. Points are listed in clock-wise direction. "
        "An empty return argument means that the play area is so far undefined.\n\n"
        "'outerboundsxyz' defines the outer area boundaries in the same way as 'playboundsxyz'.\n";
    static char seeAlsoString[] = "TestVRBoundary TestVRBoundaryPoint";

    int handle, requestVisible;
    ovrBool isVisible;
    PsychOculusDevice *oculus;
    int outFloorPointsCount, i;
    ovrVector3f *outFloorPoints;
    double *xyz;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get optional 'requestVisible' flag:
    if (PsychCopyInIntegerArg(2, kPsychArgOptional, &requestVisible)) {
        if (OVR_FAILURE(ovr_RequestBoundaryVisible(oculus->hmd, (requestVisible) ? TRUE : FALSE)))
            PsychErrorExitMsg(PsychError_system, "Failed to set request for visibility of VR boundary area.");
    }

    // Copy out current visibility status:
    if (OVR_FAILURE(ovr_GetBoundaryVisible(oculus->hmd, &isVisible)) && (verbosity > 0))
        printf("PsychOculusVRCore1-ERROR:VRAreaBoundary: Could not query VR area boundary visibility status! Results is undefined.\n");
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) ((isVisible) ? 1 : 0));

    // Copy out play area boundaries: A matrix of 2D points which define the boundary
    // in clock-wise order at floor height:
    if (OVR_FAILURE(ovr_GetBoundaryGeometry(oculus->hmd, ovrBoundary_PlayArea, (ovrVector3f*) NULL, &outFloorPointsCount))) {
        PsychErrorExitMsg(PsychError_system, "Failed to get VR boundary area I.");
    }

    outFloorPoints = (ovrVector3f*) PsychMallocTemp(outFloorPointsCount * sizeof(ovrVector3f));
    if (OVR_FAILURE(ovr_GetBoundaryGeometry(oculus->hmd, ovrBoundary_PlayArea, outFloorPoints, &outFloorPointsCount))) {
        PsychErrorExitMsg(PsychError_system, "Failed to get VR boundary area II.");
    }

    PsychAllocOutDoubleMatArg(2, kPsychArgOptional, 3, outFloorPointsCount, 1, &xyz);
    for (i = 0 ; i < outFloorPointsCount; i++) {
        *(xyz++) = (double) outFloorPoints[i].x;
        *(xyz++) = (double) outFloorPoints[i].y;
        *(xyz++) = (double) outFloorPoints[i].z;
    }

    // Copy out outer area boundaries: A matrix of 3D points which define the boundary
    // in clock-wise order at floor height:
    if (OVR_FAILURE(ovr_GetBoundaryGeometry(oculus->hmd, ovrBoundary_Outer, (ovrVector3f*) NULL, &outFloorPointsCount))) {
        PsychErrorExitMsg(PsychError_system, "Failed to get VR boundary area III.");
    }

    outFloorPoints = (ovrVector3f*) PsychMallocTemp(outFloorPointsCount * sizeof(ovrVector3f));
    if (OVR_FAILURE(ovr_GetBoundaryGeometry(oculus->hmd, ovrBoundary_Outer , outFloorPoints, &outFloorPointsCount))) {
        PsychErrorExitMsg(PsychError_system, "Failed to get VR boundary area IV.");
    }

    PsychAllocOutDoubleMatArg(3, kPsychArgOptional, 3, outFloorPointsCount, 1, &xyz);
    for (i = 0 ; i < outFloorPointsCount; i++) {
        *(xyz++) = (double) outFloorPoints[i].x;
        *(xyz++) = (double) outFloorPoints[i].y;
        *(xyz++) = (double) outFloorPoints[i].z;
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1TestVRBoundary(void)
{
    static char useString[] = "[isTriggering, closestDistance, closestPointxyz, surfaceNormal] = PsychOculusVRCore1('TestVRBoundary', oculusPtr, trackedDeviceType, boundaryType);";
    //                          1             2                3                4                                                     1          2                  3
    static char synopsisString[] =
        "Return if a tracked device associated with Oculus device 'oculusPtr' is colliding with VR area boundaries.\n\n"
        "'trackedDeviceType' Bitmask of tracked devices to check for collision. See OVR.TrackedDeviceType_XXX constants.\n\n"
        "'boundaryType' Type of VR area boundary to test against: 0 = Play area, 1 = Outer area.\n\n"
        "Return values, assuming the specified device and areas are active and/or defined:\n"
        "'isTriggering' The boundary is triggered, ie. the guardian system shows boundaries to the subject, "
        "because one of the tracked devices is close enough to a VR area boundary.\n"
        "'closestDistance' Closest distance to specified VR area boundary.\n"
        "'closestPointxyz' [x;y;z] coordinates of the closest point on the VR area boundary.\n"
        "'surfaceNormal' [nx;ny;nz] surface normal of the VR area boundary surface at 'closestPointxyz'.\n\n"
        "Note that the guardian system must be set up properly for this function to return meaningful results.\n";
    static char seeAlsoString[] = "VRAreaBoundary TestVRBoundaryPoint";

    int handle, trackedDeviceType, boundaryType;
    PsychOculusDevice *oculus;
    ovrBoundaryTestResult outTestResult;
    double *xyz;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(4));
    PsychErrorExit(PsychCapNumInputArgs(3));
    PsychErrorExit(PsychRequireNumInputArgs(3));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Type of device:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &trackedDeviceType);

    // Type of boundary:
    PsychCopyInIntegerArg(3, kPsychArgRequired, &boundaryType);
    if (boundaryType < 0 || boundaryType > 1)
        PsychErrorExitMsg(PsychError_user, "Invalid 'boundaryType' provided. Must be 0 or 1.");

    // Check for collision:
    if (OVR_FAILURE(ovr_TestBoundary(oculus->hmd, (ovrTrackedDeviceType) trackedDeviceType,
                    (boundaryType == 1) ? ovrBoundary_Outer : ovrBoundary_PlayArea,
                    &outTestResult))) {
        PsychErrorExitMsg(PsychError_system, "ovr_TestBoundary() failed.");
    }

    // Collision detected?
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) ((outTestResult.IsTriggering) ? 1 : 0));

    // Distance to closest point - the one that triggers:
    PsychCopyOutDoubleArg(2, kPsychArgOptional, (double) outTestResult.ClosestDistance);

    // 3D closest point:
    PsychAllocOutDoubleMatArg(3, kPsychArgOptional, 3, 1, 1, &xyz);
    *(xyz++) = (double) outTestResult.ClosestPoint.x;
    *(xyz++) = (double) outTestResult.ClosestPoint.y;
    *(xyz++) = (double) outTestResult.ClosestPoint.z;

    // 3D surface normal (nx,ny,nz) of closest point:
    PsychAllocOutDoubleMatArg(4, kPsychArgOptional, 3, 1, 1, &xyz);
    *(xyz++) = (double) outTestResult.ClosestPointNormal.x;
    *(xyz++) = (double) outTestResult.ClosestPointNormal.y;
    *(xyz++) = (double) outTestResult.ClosestPointNormal.z;

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1TestVRBoundaryPoint(void)
{
    static char useString[] = "[isTriggering, closestDistance, closestPointxyz, surfaceNormal] = PsychOculusVRCore1('TestVRBoundaryPoint', oculusPtr, pointxyz, boundaryType);";
    //                          1             2                3                4                                                          1          2         3
    static char synopsisString[] =
        "Return if a 3D point is colliding with VR area boundaries associated with Oculus device 'oculusPtr'.\n\n"
        "'pointxyz' [x,y,z] vector defining the 3D point position to test against the boundaries.\n\n"
        "'boundaryType' Type of VR area boundary to test against: 0 = Play area, 1 = Outer area.\n\n"
        "Return values, assuming the specified area is defined:\n"
        "'isTriggering' The boundary is triggered, ie. the guardian system would show boundaries for given 'pointxyz'.\n"
        "'closestDistance' Closest distance of point to specified VR area boundary.\n"
        "'closestPointxyz' [x;y;z] coordinates of the closest point on the VR area boundary.\n"
        "'surfaceNormal' [nx;ny;nz] surface normal of the VR area boundary surface at 'closestPointxyz'.\n\n"
        "Note that the guardian system must be set up properly for this function to return meaningful results.\n";
    static char seeAlsoString[] = "VRAreaBoundary TestVRBoundary";

    int handle, boundaryType;
    PsychOculusDevice *oculus;
    ovrBoundaryTestResult outTestResult;
    ovrVector3f point;
    int m, n, p;
    double *xyz;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(4));
    PsychErrorExit(PsychCapNumInputArgs(3));
    PsychErrorExit(PsychRequireNumInputArgs(3));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Test point:
    PsychAllocInDoubleMatArg(2, kPsychArgRequired, &m, &n, &p, &xyz);
    if (m * n * p != 3)
        PsychErrorExitMsg(PsychError_user, "Invalid 'pointxyz' provided. Must be a 3D vector [x,y,z].");

    point.x = (float) *(xyz++);
    point.y = (float) *(xyz++);
    point.z = (float) *(xyz++);

    // Type of boundary:
    PsychCopyInIntegerArg(3, kPsychArgRequired, &boundaryType);
    if (boundaryType < 0 || boundaryType > 1)
        PsychErrorExitMsg(PsychError_user, "Invalid 'boundaryType' provided. Must be 0 or 1.");

    // Check for collision:
    if (OVR_FAILURE(ovr_TestBoundaryPoint(oculus->hmd, (const ovrVector3f*) &point,
                    (boundaryType == 1) ? ovrBoundary_Outer : ovrBoundary_PlayArea,
                    &outTestResult))) {
        PsychErrorExitMsg(PsychError_system, "ovr_TestBoundaryPoint() failed.");
    }

    // Collision detected?
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) ((outTestResult.IsTriggering) ? 1 : 0));

    // Distance to closest point - the one that triggers:
    PsychCopyOutDoubleArg(2, kPsychArgOptional, (double) outTestResult.ClosestDistance);

    // 3D closest point:
    PsychAllocOutDoubleMatArg(3, kPsychArgOptional, 3, 1, 1, &xyz);
    *(xyz++) = (double) outTestResult.ClosestPoint.x;
    *(xyz++) = (double) outTestResult.ClosestPoint.y;
    *(xyz++) = (double) outTestResult.ClosestPoint.z;

    // 3D surface normal (nx,ny,nz) of closest point:
    PsychAllocOutDoubleMatArg(4, kPsychArgOptional, 3, 1, 1, &xyz);
    *(xyz++) = (double) outTestResult.ClosestPointNormal.x;
    *(xyz++) = (double) outTestResult.ClosestPointNormal.y;
    *(xyz++) = (double) outTestResult.ClosestPointNormal.z;

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1Start(void)
{
    static char useString[] = "PsychOculusVRCore1('Start', oculusPtr);";
    //                                                     1
    static char synopsisString[] =
        "Start head orientation and position tracking operation on Oculus device 'oculusPtr'.\n\n";
    static char seeAlsoString[] = "Stop";

    int handle;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(0));
    PsychErrorExit(PsychCapNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    if (oculus->isTracking) {
        if (verbosity >= 0)
            printf("PsychOculusVRCore1-ERROR: Tried to start tracking on device %i, but tracking is already started.\n", handle);
        PsychErrorExitMsg(PsychError_user, "Tried to start tracking on HMD, but tracking already active.");
    }

    if (verbosity >= 4)
        printf("PsychOculusVRCore1-INFO: Tracking started on device with handle %i.\n", handle);

    // presenterThread shutdown: Ask thread to terminate, wait for thread termination, cleanup and release the thread:
    PsychLockMutex(&(oculus->presenterLock));
    oculus->closing = TRUE;
    PsychUnlockMutex(&(oculus->presenterLock));

    if (oculus->presenterThread != (psych_thread) NULL) {
        if (verbosity > 5)
            printf("PTB-DEBUG: Waiting (join()ing) for helper thread of HMD %p to finish up. If this doesn't happen quickly, you'll have to kill Octave...\n", oculus);

        PsychDeleteThread(&(oculus->presenterThread));
    }

    // Ok, thread is dead. Mark it as such:
    oculus->presenterThread = (psych_thread) NULL;
    oculus->closing = FALSE;
    oculus->isTracking = TRUE;

    // Tracking is running.
    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1Stop(void)
{
    static char useString[] = "PsychOculusVRCore1('Stop', oculusPtr);";
    static char synopsisString[] =
        "Stop head tracking operation on Oculus device 'oculusPtr'.\n\n";
    static char seeAlsoString[] = "Start";

    int handle;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(0));
    PsychErrorExit(PsychCapNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);

    // Stop device:
    PsychOculusStop(handle);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetTrackingState(void)
{
    static char useString[] = "[state, touch] = PsychOculusVRCore1('GetTrackingState', oculusPtr [, predictionTime=nextFrame]);";
    //                          1      2                                               1            2
    static char synopsisString[] =
        "Return current state of head position and orientation tracking for Oculus device 'oculusPtr'.\n"
        "Head position and orientation is predicted for target time 'predictionTime' in seconds if provided, "
        "based on the latest measurements from the tracking hardware. If 'predictionTime' is set to zero, "
        "then no prediction is performed and the current state based on latest measurements is returned.\n"
        "If 'predictionTime' is omitted, then the prediction is performed for the mid-point of the next "
        "possible video frame of the HMD, ie. the most likely presentation time for immediately rendered images.\n\n"
        "'state' is a struct with fields reporting the following values:\n"
        "'Time' = Time in seconds of returned tracking state.\n"
        "'Status' = Tracking status flags:\n"
        "+1 = Head orientation tracked,\n"
        "+2 = Head position tracked,\n"
        "+4 = At least one active tracking camera has a valid pose for absolute tracking,\n"
        "+32 = At least one tracking camera is connected and online,\n"
        "+128 = HMD display is connected, available and actually on users head, displaying our content.\n\n"
        "'SessionState' = VR session status flags, added together:\n"
        "+1  = Our rendering goes to the HMD, ie. we have control over it. If the Oculus store app would "
        "be in control, or if the Health and Safety warning would display, this flag would be missing.\n"
        "+2  = HMD is present and active.\n"
        "+4  = HMD is strapped onto users head. A Rift CV1 would switch off/blank if not on the head.\n"
        "+8  = DisplayLost condition! Some hardware/software malfunction, need to completely quit to recover.\n"
        "+16 = ShouldQuit The user interface asks us to voluntarily terminate this session.\n"
        "+32 = ShouldRecenter = The user interface asks us to recenter/recalibrate our tracking origin.\n\n"
        "'HeadPose' = Head position [x, y, z] in meters and rotation as quaternion [rx, ry, rz, rw], all as a vector [x,y,z,rx,ry,rz,rw].\n"
        "'HeadLinearSpeed' = Head linear velocity [vx,vy,vz] in meters/sec.\n"
        "'HeadAngularSpeed' = Head angular velocity [rx,ry,rz] in radians/sec.\n"
        "'HeadLinearAcceleration' = Head linear acceleration [ax,ay,az] in meters/sec^2.\n"
        "'HeadAngularAcceleration' = Head angular acceleration [rax,ray,raz] in radians/sec^2.\n"
        "'CalibratedOrigin' = The pose of the world coordinate system origin during last calibration. "
        "Units and format are like 'HeadPose' ie. a vector [x,y,z,rx,ry,rz,rw].\n"
        "\n"
        "Touch controller position and orientation:\n\n"
        "The return argument 'touch' is a struct array with 2 structs. touch(1) contains info about "
        "the tracking state and tracked pose of the left hand (= left touch controller) of the user, "
        "touch(2) contains info about the right hand (= right touch controller) of the user.\n"
        "The structs have very similar structure to the head (= HMD) tracking data returned by 'state':\n\n"
        "'Time' = Time in seconds of returned hand/controller tracking state.\n"
        "'Status' = Tracking status flags:\n"
        "0  = No tracking info for hand/controller, ie. no Oculus touch sensor connected.\n"
        "+1 = Hand orientation tracked,\n"
        "+2 = Hand position tracked,\n"
        "'HandPose' = Position and orientation of the hand, in usual [x,y,z,rx,ry,rz,rw] vector form as with 'HeadPose'.\n"
        "'HandLinearSpeed' = Hand linear velocity [vx,vy,vz] in meters/sec.\n"
        "'HandAngularSpeed' = Hand angular velocity [rx,ry,rz] in radians/sec.\n"
        "'HandLinearAcceleration' = Hand linear acceleration [ax,ay,az] in meters/sec^2.\n"
        "'HandAngularAcceleration' = Hand angular acceleration [rax,ray,raz] in radians/sec^2.\n"
        "\n";

    static char seeAlsoString[] = "Start Stop GetTrackersState GetInputState";

    PsychGenericScriptType *status;
    const char *FieldNames1[] = {"Time", "Status", "SessionState", "HeadPose", "HeadLinearSpeed", "HeadAngularSpeed",
                                 "HeadLinearAcceleration", "HeadAngularAcceleration", "CalibratedOrigin"};
    const int FieldCount1 = 9;

    const char *FieldNames2[] = {"Time", "Status", "HandPose", "HandLinearSpeed", "HandAngularSpeed",
                                 "HandLinearAcceleration", "HandAngularAcceleration"};
    const int FieldCount2 = 7;

    PsychGenericScriptType *outMat;
    double *v;
    int handle, trackerCount, i;
    double predictionTime;
    int StatusFlags = 0;
    PsychOculusDevice *oculus;
    ovrTrackingState state;
    ovrVector3f HmdToEyeOffset[2];
    ovrTrackerPose trackerPose;
    ovrSessionStatus sessionStatus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(2));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    if (OVR_FAILURE(ovr_GetSessionStatus(oculus->hmd, &sessionStatus))) {
        ovr_GetLastErrorInfo(&errorInfo);
        if (verbosity > 0) printf("PsychOculusVRCore1-ERROR: ovr_GetSessionStatus failed: %s\n", errorInfo.ErrorString);
        PsychErrorExitMsg(PsychError_system, "Failed to get current session status from VR compositor.");
    }

    PsychLockMutex(&(oculus->presenterLock));

    // Get optional target time for predicted tracking state. Default to the
    // predicted state for the predicted mid-point of the next video frame:
    if (!PsychCopyInDoubleArg(2, kPsychArgOptional, &predictionTime))
        predictionTime = ovr_GetPredictedDisplayTime(oculus->hmd, oculus->frameIndex);

    // Get current tracking status info at time predictionTime. Mark this point
    // as time from which motion to photon latency is measured (latencymarker = TRUE):
    state = ovr_GetTrackingState(oculus->hmd, predictionTime, TRUE);
    oculus->sensorSampleTime = ovr_GetTimeInSeconds();

    // Translate to per eye position and orientation:
    HmdToEyeOffset[0] = oculus->eyeRenderDesc[0].HmdToEyeOffset;
    HmdToEyeOffset[1] = oculus->eyeRenderDesc[1].HmdToEyeOffset;
    ovr_CalcEyePoses(state.HeadPose.ThePose, HmdToEyeOffset, oculus->outEyePoses);
    PsychUnlockMutex(&(oculus->presenterLock));

    // Print out tracking status:
    if (verbosity >= 4) {
        printf("PsychOculusVRCore1-INFO: Tracking state predicted for device %i at time %f.\n", handle, predictionTime);
        printf("PsychOculusVRCore1-INFO: Time %f : Status %i\n", state.HeadPose.TimeInSeconds, state.StatusFlags);
        printf("PsychOculusVRCore1-INFO: HeadPose: Position    [x,y,z]   = [%f, %f, %f]\n", state.HeadPose.ThePose.Position.x, state.HeadPose.ThePose.Position.y, state.HeadPose.ThePose.Position.z);
        printf("PsychOculusVRCore1-INFO: HeadPose: Orientation [x,y,z,w] = [%f, %f, %f, %f]\n", state.HeadPose.ThePose.Orientation.x, state.HeadPose.ThePose.Orientation.y, state.HeadPose.ThePose.Orientation.z, state.HeadPose.ThePose.Orientation.w);
    }

    PsychAllocOutStructArray(1, kPsychArgOptional, 1, FieldCount1, FieldNames1, &status);

    PsychSetStructArrayDoubleElement("Time", 0, state.HeadPose.TimeInSeconds, status);

    StatusFlags = state.StatusFlags & (ovrStatus_OrientationTracked | ovrStatus_PositionTracked);

    // HMD present, connected and online and on users head, displaying us?
    if (sessionStatus.HmdPresent && sessionStatus.IsVisible && sessionStatus.HmdMounted)
        StatusFlags |= 128;

    // At least one tracker connected?
    trackerCount = ovr_GetTrackerCount(oculus->hmd);
    for (i = 0; i < trackerCount; i++) {
        trackerPose = ovr_GetTrackerPose(oculus->hmd, i);
        StatusFlags |= trackerPose.TrackerFlags & (ovrTracker_Connected | ovrTracker_PoseTracked);
    }

    // Return head and general tracking status flags:
    PsychSetStructArrayDoubleElement("Status", 0, StatusFlags, status);

    // Return sesstion status flags:
    StatusFlags = 0;
    if (sessionStatus.IsVisible) StatusFlags |= 1;
    if (sessionStatus.HmdPresent) StatusFlags |= 2;
    if (sessionStatus.HmdMounted) StatusFlags |= 4;
    if (sessionStatus.DisplayLost) StatusFlags |= 8;
    if (sessionStatus.ShouldQuit) StatusFlags |= 16;
    if (sessionStatus.ShouldRecenter) StatusFlags |= 32;

    PsychSetStructArrayDoubleElement("SessionState", 0, StatusFlags, status);

    // Head pose:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 7, 1, &v, &outMat);
    v[0] = state.HeadPose.ThePose.Position.x;
    v[1] = state.HeadPose.ThePose.Position.y;
    v[2] = state.HeadPose.ThePose.Position.z;

    v[3] = state.HeadPose.ThePose.Orientation.x;
    v[4] = state.HeadPose.ThePose.Orientation.y;
    v[5] = state.HeadPose.ThePose.Orientation.z;
    v[6] = state.HeadPose.ThePose.Orientation.w;
    PsychSetStructArrayNativeElement("HeadPose", 0, outMat, status);

    // Linear velocity:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
    v[0] = state.HeadPose.LinearVelocity.x;
    v[1] = state.HeadPose.LinearVelocity.y;
    v[2] = state.HeadPose.LinearVelocity.z;
    PsychSetStructArrayNativeElement("HeadLinearSpeed", 0, outMat, status);

    // Angular velocity:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
    v[0] = state.HeadPose.AngularVelocity.x;
    v[1] = state.HeadPose.AngularVelocity.y;
    v[2] = state.HeadPose.AngularVelocity.z;
    PsychSetStructArrayNativeElement("HeadAngularSpeed", 0, outMat, status);

    // Linear acceleration:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
    v[0] = state.HeadPose.LinearAcceleration.x;
    v[1] = state.HeadPose.LinearAcceleration.y;
    v[2] = state.HeadPose.LinearAcceleration.z;
    PsychSetStructArrayNativeElement("HeadLinearAcceleration", 0, outMat, status);

    // Angular acceleration:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
    v[0] = state.HeadPose.AngularAcceleration.x;
    v[1] = state.HeadPose.AngularAcceleration.y;
    v[2] = state.HeadPose.AngularAcceleration.z;
    PsychSetStructArrayNativeElement("HeadAngularAcceleration", 0, outMat, status);

    // CalibratedOrigin:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 7, 1, &v, &outMat);
    v[0] = state.CalibratedOrigin.Position.x;
    v[1] = state.CalibratedOrigin.Position.y;
    v[2] = state.CalibratedOrigin.Position.z;

    v[3] = state.CalibratedOrigin.Orientation.x;
    v[4] = state.CalibratedOrigin.Orientation.y;
    v[5] = state.CalibratedOrigin.Orientation.z;
    v[6] = state.CalibratedOrigin.Orientation.w;
    PsychSetStructArrayNativeElement("CalibratedOrigin", 0, outMat, status);

    // Now the tracking info from the Oculus touch controllers 0 and 1 for left
    // and right hand, in a separate struct array:
    PsychAllocOutStructArray(2, kPsychArgOptional, 2, FieldCount2, FieldNames2, &status);

    for (i = 0; i < 2; i++) {
        // Timestamp for when this tracking info is valid:
        PsychSetStructArrayDoubleElement("Time", i, state.HandPoses[i].TimeInSeconds, status);

        // Hand / touch controller tracking state:
        StatusFlags = state.HandStatusFlags[i] & (ovrStatus_OrientationTracked | ovrStatus_PositionTracked);
        PsychSetStructArrayDoubleElement("Status", i, StatusFlags, status);

        // Hand pose:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 7, 1, &v, &outMat);
        v[0] = state.HandPoses[i].ThePose.Position.x;
        v[1] = state.HandPoses[i].ThePose.Position.y;
        v[2] = state.HandPoses[i].ThePose.Position.z;

        v[3] = state.HandPoses[i].ThePose.Orientation.x;
        v[4] = state.HandPoses[i].ThePose.Orientation.y;
        v[5] = state.HandPoses[i].ThePose.Orientation.z;
        v[6] = state.HandPoses[i].ThePose.Orientation.w;
        PsychSetStructArrayNativeElement("HandPose", i, outMat, status);

        // Linear velocity:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
        v[0] = state.HandPoses[i].LinearVelocity.x;
        v[1] = state.HandPoses[i].LinearVelocity.y;
        v[2] = state.HandPoses[i].LinearVelocity.z;
        PsychSetStructArrayNativeElement("HandLinearSpeed", i, outMat, status);

        // Angular velocity:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
        v[0] = state.HandPoses[i].AngularVelocity.x;
        v[1] = state.HandPoses[i].AngularVelocity.y;
        v[2] = state.HandPoses[i].AngularVelocity.z;
        PsychSetStructArrayNativeElement("HandAngularSpeed", i, outMat, status);

        // Linear acceleration:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
        v[0] = state.HandPoses[i].LinearAcceleration.x;
        v[1] = state.HandPoses[i].LinearAcceleration.y;
        v[2] = state.HandPoses[i].LinearAcceleration.z;
        PsychSetStructArrayNativeElement("HandLinearAcceleration", i, outMat, status);

        // Angular acceleration:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
        v[0] = state.HandPoses[i].AngularAcceleration.x;
        v[1] = state.HandPoses[i].AngularAcceleration.y;
        v[2] = state.HandPoses[i].AngularAcceleration.z;
        PsychSetStructArrayNativeElement("HandAngularAcceleration", i, outMat, status);
    }
    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetInputState(void)
{
    static char useString[] = "input = PsychOculusVRCore1('GetInputState', oculusPtr, controllerType);";
    //                         1                                           1          2
    static char synopsisString[] =
        "Return current state of input device 'controllerType' associated with Oculus device 'oculusPtr'.\n\n"
        "'controllerType' can be one of the follwing values:\n"
        "OVR.ControllerType_LTouch = Left touch controller (Left tracked hand).\n"
        "OVR.ControllerType_RTouch = Right touch controller (Right tracked hand).\n"
        "OVR.ControllerType_Remote = Oculus remote control.\n"
        "OVR.ControllerType_XBox = XBox controller.\n"
        "OVR.ControllerType_Active = Whatever controller is connected and active.\n"
        "\n"
        "'input' is a struct with fields reporting the following status values of the controller:\n"
        "'Valid' = 1 if 'input' contains valid results, 0 if input status is invalid/unavailable.\n"
        "'Time' = Time in seconds when controller state was last updated.\n"
        "'Buttons' = Vector with each positions value corresponding to a specifc button being pressed (1) "
        "or released (0). The OVR.Button_XXX constants map button names to vector indices (like KbName() "
        "does for KbCheck()).\n"
        "'Touches' = Vector with touch values as described by the OVR.Touch_XXX constants. Works like 'Buttons'.\n"
        "'Trigger'(1/2) = Left (1) and Right (2) trigger: Value range 0.0 - 1.0, filtered and with dead-zone.\n"
        "'TriggerNoDeadzone'(1/2) = Left (1) and Right (2) trigger: Value range 0.0 - 1.0, filtered.\n"
        "'TriggerRaw'(1/2) = Left (1) and Right (2) trigger: Value range 0.0 - 1.0, raw values unfiltered.\n"
        "'Grip'(1/2) = Left (1) and Right (2) grip button: Value range 0.0 - 1.0, filtered and with dead-zone.\n"
        "'GripNoDeadzone'(1/2) = Left (1) and Right (2) grip button: Value range 0.0 - 1.0, filtered.\n"
        "'GripRaw'(1/2) = Left (1) and Right (2) grip button: Value range 0.0 - 1.0, raw values unfiltered.\n"
        "'Thumbstick' = 2x2 matrix: Column 1 contains left thumbsticks [x;y] axis values, column 2 contains "
        "right sticks [x;y] axis values. Values are in range -1 to +1, filtered and with deadzone applied.\n"
        "'ThumbstickNoDeadzone' = Like 'Thumbstick', filtered, but without a deadzone applied.\n"
        "'ThumbstickRaw' = 'Thumbstick' raw date without deadzone or filtering applied.\n"
        "\n";

    static char seeAlsoString[] = "Start Stop GetTrackedState GetTrackersState";

    PsychGenericScriptType *status;
    const char *FieldNames[] = { "Valid", "Time", "Buttons", "Touches", "Trigger", "Grip", "TriggerNoDeadzone",
                                 "GripNoDeadzone", "TriggerRaw", "GripRaw", "Thumbstick",
                                 "ThumbstickNoDeadzone", "ThumbstickRaw" };
    const int FieldCount = 13;

    PsychGenericScriptType *outMat;
    double *v;
    int handle, i;
    unsigned long controllerType;
    double controllerTypeD;
    PsychOculusDevice *oculus;
    ovrInputState state;
    ovrSessionStatus sessionStatus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Controller type:
    PsychCopyInDoubleArg(2, kPsychArgRequired, &controllerTypeD);
    controllerType = (unsigned long) controllerTypeD;

    PsychAllocOutStructArray(1, kPsychArgOptional, 1, FieldCount, FieldNames, &status);

    // Check session status if we have VR input focus:
    if (OVR_FAILURE(ovr_GetSessionStatus(oculus->hmd, &sessionStatus))) {
        ovr_GetLastErrorInfo(&errorInfo);
        if (verbosity > 0) printf("PsychOculusVRCore1-ERROR: ovr_GetSessionStatus failed: %s\n", errorInfo.ErrorString);
        PsychErrorExitMsg(PsychError_system, "Failed to get current session status from VR compositor.");
    }

    if (sessionStatus.IsVisible) {
        // Have VR input focus, try to get input state:
        if (OVR_FAILURE(ovr_GetInputState(oculus->hmd, (ovrControllerType) controllerType, &state))) {
            ovr_GetLastErrorInfo(&errorInfo);
            if (verbosity > 0)
                printf("PsychOculusVRCore1-ERROR: ovr_GetInputState() for controller %f 0x%lx failed: %s\n",
                       controllerTypeD, controllerType, errorInfo.ErrorString);
            PsychErrorExitMsg(PsychError_system, "Failed to get some controller input status.");
        }

        // Mark as valid:
        PsychSetStructArrayDoubleElement("Valid", 0, 1, status);
    }
    else {
        // No VR input focus -> No valid input -> Mark as invalid:
        PsychSetStructArrayDoubleElement("Valid", 0, 0, status);

        // Fail gracefully, by returning an all-zero input struct:
        memset(&state, 0, sizeof(state));
    }

    // Controller update time:
    PsychSetStructArrayDoubleElement("Time", 0, state.TimeInSeconds, status);

    // Button states:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 32, 1, &v, &outMat);
    for (i = 0; i < 32; i++)
        v[i] = (state.Buttons & (1 << i)) ? 1 : 0;
    PsychSetStructArrayNativeElement("Buttons", 0, outMat, status);

    // Touch states:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 32, 1, &v, &outMat);
    for (i = 0; i < 32; i++)
        v[i] = (state.Touches & (1 << i)) ? 1 : 0;
    PsychSetStructArrayNativeElement("Touches", 0, outMat, status);

    // Trigger left/right:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
    v[0] = state.IndexTrigger[0];
    v[1] = state.IndexTrigger[1];
    PsychSetStructArrayNativeElement("Trigger", 0, outMat, status);

    // Grip left/right:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
    v[0] = state.HandTrigger[0];
    v[1] = state.HandTrigger[1];
    PsychSetStructArrayNativeElement("Grip", 0, outMat, status);

    // TriggerNoDeadzone left/right:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
    v[0] = state.IndexTriggerNoDeadzone[0];
    v[1] = state.IndexTriggerNoDeadzone[1];
    PsychSetStructArrayNativeElement("TriggerNoDeadzone", 0, outMat, status);

    // GripNoDeadzone left/right:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
    v[0] = state.HandTriggerNoDeadzone[0];
    v[1] = state.HandTriggerNoDeadzone[1];
    PsychSetStructArrayNativeElement("GripNoDeadzone", 0, outMat, status);

    // TriggerRaw left/right:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
    v[0] = state.IndexTriggerRaw[0];
    v[1] = state.IndexTriggerRaw[1];
    PsychSetStructArrayNativeElement("TriggerRaw", 0, outMat, status);

    // GripRaw left/right:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
    v[0] = state.HandTriggerRaw[0];
    v[1] = state.HandTriggerRaw[1];
    PsychSetStructArrayNativeElement("GripRaw", 0, outMat, status);

    // Thumbstick: column 1 = left hand, column 2 = right hand. row 1 = x, row 2= y:
    v = NULL;
    PsychAllocateNativeDoubleMat(2, 2, 1, &v, &outMat);
    v[0] = state.Thumbstick[0].x;
    v[1] = state.Thumbstick[0].y;
    v[2] = state.Thumbstick[1].x;
    v[3] = state.Thumbstick[1].y;
    PsychSetStructArrayNativeElement("Thumbstick", 0, outMat, status);

    // ThumbstickNoDeadzone:
    v = NULL;
    PsychAllocateNativeDoubleMat(2, 2, 1, &v, &outMat);
    v[0] = state.ThumbstickNoDeadzone[0].x;
    v[1] = state.ThumbstickNoDeadzone[0].y;
    v[2] = state.ThumbstickNoDeadzone[1].x;
    v[3] = state.ThumbstickNoDeadzone[1].y;
    PsychSetStructArrayNativeElement("ThumbstickNoDeadzone", 0, outMat, status);

    // ThumbstickRaw:
    v = NULL;
    PsychAllocateNativeDoubleMat(2, 2, 1, &v, &outMat);
    v[0] = state.ThumbstickRaw[0].x;
    v[1] = state.ThumbstickRaw[0].y;
    v[2] = state.ThumbstickRaw[1].x;
    v[3] = state.ThumbstickRaw[1].y;
    PsychSetStructArrayNativeElement("ThumbstickRaw", 0, outMat, status);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetFovTextureSize(void)
{
    static char useString[] = "[width, height, fovPort] = PsychOculusVRCore1('GetFovTextureSize', oculusPtr, eye [, fov=[HMDRecommended]][, pixelsPerDisplay=1]);";
    //                          1      2       3                                                  1          2      3                       4
    static char synopsisString[] =
    "Return recommended size of client renderbuffers for Oculus device 'oculusPtr'.\n"
    "'eye' which eye to provide the size for: 0 = Left, 1 = Right.\n"
    "'fov' Optional field of view in degrees, from line of sight: [leftdeg, rightdeg, updeg, downdeg]. "
    "If 'fov' is omitted, the Oculus runtime will be asked for a good default field of view and that "
    "will be used. The field of view may be dependent on the settings in the Oculus user profile of the "
    "currently selected user.\n"
    "'pixelsPerDisplay' Ratio of the number of render target pixels to display pixels at the center "
    "of distortion. Defaults to 1.0 if omitted. Lower values can improve performance, at lower quality.\n"
    "\n"
    "Return values are 'width' for minimum recommended width of framebuffer in pixels and "
    "'height' for minimum recommended height of framebuffer in pixels. 'fovPort' is the field of view "
    "in degrees finally used for calculation of 'width' x 'height'.\n";
    static char seeAlsoString[] = "GetUndistortionParameters";

    int handle, eyeIndex;
    int n, m, p;
    PsychOculusDevice *oculus;
    double *fov;
    double pixelsPerDisplay;
    double *outFov;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(4));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get eye index - left = 0, right = 1:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &eyeIndex);
    if (eyeIndex < 0 || eyeIndex > 1) PsychErrorExitMsg(PsychError_user, "Invalid 'eye' specified. Must be 0 or 1 for left- or right eye.");

    // Get optional field of view in degrees in left,right,up,down direction from line of sight:
    if (PsychAllocInDoubleMatArg(3, kPsychArgOptional, &n, &m, &p, &fov)) {
        // Validate and assign:
        if (n * m * p != 4) PsychErrorExitMsg(PsychError_user, "Invalid 'fov' specified. Must be a 4-component vector of form [leftdeg, rightdeg, updeg, downdeg].");
        oculus->ofov[eyeIndex].LeftTan  = (float) tan(deg2rad(fov[0]));
        oculus->ofov[eyeIndex].RightTan = (float) tan(deg2rad(fov[1]));
        oculus->ofov[eyeIndex].UpTan    = (float) tan(deg2rad(fov[2]));
        oculus->ofov[eyeIndex].DownTan  = (float) tan(deg2rad(fov[3]));
    }
    else {
        // None specified: Ask the runtime for good defaults.
        oculus->ofov[eyeIndex] = oculus->hmdDesc.DefaultEyeFov[eyeIndex];
    }

    // Get optional pixelsPerDisplay parameter:
    pixelsPerDisplay = 1.0;
    PsychCopyInDoubleArg(4, kPsychArgOptional, &pixelsPerDisplay);
    if (pixelsPerDisplay <= 0.0) PsychErrorExitMsg(PsychError_user, "Invalid 'pixelsPerDisplay' specified. Must be greater than zero.");

    // Ask the api for optimal texture size, aka the size of the client draw buffer:
    oculus->texSize[eyeIndex] = ovr_GetFovTextureSize(oculus->hmd, (ovrEyeType) eyeIndex, oculus->ofov[eyeIndex], (float) pixelsPerDisplay);

    // Return recommended width and height of drawBuffer:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, oculus->texSize[eyeIndex].w);
    PsychCopyOutDoubleArg(2, kPsychArgOptional, oculus->texSize[eyeIndex].h);

    PsychAllocOutDoubleMatArg(3, kPsychArgOptional, 4, 1, 1, &outFov);
    outFov[0] = rad2deg(atan(oculus->ofov[eyeIndex].LeftTan));
    outFov[1] = rad2deg(atan(oculus->ofov[eyeIndex].RightTan));
    outFov[2] = rad2deg(atan(oculus->ofov[eyeIndex].UpTan));
    outFov[3] = rad2deg(atan(oculus->ofov[eyeIndex].DownTan));

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1CreateRenderTextureChain(void)
{
    static char useString[] = "[width, height, numTextures] = PsychOculusVRCore1('CreateRenderTextureChain', oculusPtr, eye, width, height, floatFormat);";
    //                          1      2       3                                                             1          2    3      4       5
    static char synopsisString[] =
    "Create texture present chains for Oculus device 'oculusPtr'.\n\n"
    "'eye' Eye for which chain should get created: 0 = Left/Mono, 1 = Right.\n"
    "If only a chain for eye = 0 is created then the driver operates in monoscopic "
    "presentation mode for use with Screen() stereomode 0, showing the same mono image to both "
    "eyes. If a 2nd chain for eye = 1 is created then the driver switches to stereoscopic "
    "presentation mode for use with Screen() stereomode 12, presenting separate images to the left "
    "and right eye.\n\n"
    "'width' and 'height' are the width x height of the texture into which Psychtoolbox "
    "Screen() image processing pipeline will render the output image of an eye for submission "
    "to the VR compositor. Left and right eye must use identical 'width' and 'height'.\n\n"
    "'floatFormat' Texture format: 0 = RGBA8 sRGB format for sRGB rendering and output. 1 = 16 bpc "
    "half-float RGBA16F in linear format.\n\n"
    "'numTextures' returns the total number of compositor textures in the swap chain.\n"
    "\n"
    "Return values are 'width' for selected width of output texture in pixels and "
    "'height' for height of output texture in pixels.\n";
    static char seeAlsoString[] = "GetNextTextureHandle";

    int handle, eyeIndex;
    int width, height, out_Length, floatFormat;
    PsychOculusDevice *oculus;
    ovrTextureSwapChainDesc chainDesc;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(5));
    PsychErrorExit(PsychRequireNumInputArgs(5));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get eye:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &eyeIndex);
    if (eyeIndex < 0 || eyeIndex > 1)
        PsychErrorExitMsg(PsychError_user, "Invalid 'eye' specified. Must be 0 or 1 for left- or right eye.");

    // Get texture dimensions:
    PsychCopyInIntegerArg(3, kPsychArgRequired, &width);
    if (width < 1)
        PsychErrorExitMsg(PsychError_user, "Invalid width, smaller than 1 texel!");

    PsychCopyInIntegerArg(4, kPsychArgRequired, &height);
    if (height < 1)
        PsychErrorExitMsg(PsychError_user, "Invalid height, smaller than 1 texel!");

    // Get texture format: 0 = sRGB RGBA8, 1 = linear RGBA16F half-float:
    PsychCopyInIntegerArg(5, kPsychArgRequired, &floatFormat);
    if (floatFormat < 0 || floatFormat > 1)
        PsychErrorExitMsg(PsychError_user, "Invalid floatFormat flag provided. Must be 0 or 1.");

    if (oculus->textureSwapChain[eyeIndex])
        PsychErrorExitMsg(PsychError_user, "Tried to create already created texture swap chain for given eye.");

    if (eyeIndex > 0 && (width != oculus->textureWidth || height != oculus->textureHeight))
        PsychErrorExitMsg(PsychError_user, "Given width x height for 2nd eye does not match width x height of 1st eye, as required.");

    // Build OpenGL texture chain descriptor:
    memset(&chainDesc, 0, sizeof(chainDesc));
    chainDesc.Type = ovrTexture_2D;
    chainDesc.Format = (floatFormat == 1) ? OVR_FORMAT_R16G16B16A16_FLOAT : OVR_FORMAT_R8G8B8A8_UNORM_SRGB;
    chainDesc.ArraySize = 1;
    chainDesc.MipLevels = 1;
    chainDesc.SampleCount = 1;
    chainDesc.StaticImage = ovrFalse;
    chainDesc.MiscFlags = ovrTextureMisc_None;
    chainDesc.BindFlags = ovrTextureBind_None;
    chainDesc.Width = width;
    chainDesc.Height = height;

    // Create texture swap chain:
    PsychLockMutex(&(oculus->presenterLock));

    if (OVR_FAILURE(ovr_CreateTextureSwapChainGL(oculus->hmd, &chainDesc, &oculus->textureSwapChain[eyeIndex]))) {
        ovr_GetLastErrorInfo(&errorInfo);
        if (verbosity > 0) printf("PsychOculusVRCore1-ERROR: ovr_CreateTextureSwapChainGL failed: %s\n", errorInfo.ErrorString);
        PsychUnlockMutex(&(oculus->presenterLock));
        PsychErrorExitMsg(PsychError_system, "Failed to create texture swap chain for VR compositor.");
    }

    // Mark driver as in stereo mode if a swap chain for the right eye was created:
    if (eyeIndex > 0) {
        if (verbosity > 2) printf("PsychOculusVRCore1-INFO: Right eye swap chain created. Switching to stereo mode.\n");
        oculus->isStereo = TRUE;
    }

    ovr_GetTextureSwapChainLength(oculus->hmd, oculus->textureSwapChain[eyeIndex], &out_Length);
    oculus->textureSwapChainLength = out_Length;

    if (verbosity > 3)
        printf("PsychOculusVRCore1-INFO: Allocated texture swap chain has %i buffers.\n", out_Length);

    // Assign total texture buffer width/height for the frame submission later on:
    oculus->textureWidth = chainDesc.Width;
    oculus->textureHeight = chainDesc.Height;

    PsychUnlockMutex(&(oculus->presenterLock));

    // Return recommended width and height of drawBuffer:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, chainDesc.Width);
    PsychCopyOutDoubleArg(2, kPsychArgOptional, chainDesc.Height);

    // Return number of textures in swap chain:
    PsychCopyOutDoubleArg(3, kPsychArgOptional, out_Length);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetNextTextureHandle(void)
{
    static char useString[] = "texObjectHandle = PsychOculusVRCore1('GetNextTextureHandle', oculusPtr, eye);";
    //                         1                                                            1          2
    static char synopsisString[] =
    "Retrieve OpenGL texture object handle for next target texture for Oculus device 'oculusPtr'.\n"
    "'eye' Eye for which handle of next texture should be returned: 0 = Left/Mono, 1 = Right.\n"
    "Returns a GL_TEXTURE2D texture object name/handle in 'texObjectHandle' for the texture "
    "to which the next VR frame should be rendered. Returns -1 if busy.\n";
    static char seeAlsoString[] = "CreateRenderTextureChain";

    int handle, eyeIndex;
    unsigned int texObjectHandle;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get eye:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &eyeIndex);
    if (eyeIndex < 0 || eyeIndex > 1)
        PsychErrorExitMsg(PsychError_user, "Invalid 'eye' specified. Must be 0 or 1 for left- or right eye.");

    if (eyeIndex > 0 && !(oculus->isStereo))
        PsychErrorExitMsg(PsychError_user, "Invalid 'eye' specified. Must be 0, as mono display mode is selected.");

    // Get next free/non-busy buffer for this eyes texture swap chain:
    PsychLockMutex(&(oculus->presenterLock));
    /*
    if (OVR_FAILURE(ovr_GetTextureSwapChainCurrentIndex(oculus->hmd, oculus->textureSwapChain[eyeIndex], &out_Index))) {
        ovr_GetLastErrorInfo(&errorInfo);
        if (verbosity > 0) printf("PsychOculusVRCore1-ERROR: eye %i ovr_GetTextureSwapChainCurrentIndex failed: %s\n", eyeIndex, errorInfo.ErrorString);
        PsychUnlockMutex(&(oculus->presenterLock));
        PsychErrorExitMsg(PsychError_system, "Failed to retrieve next OpenGL texture from swap chain.");
    }
    else if (verbosity > 2)
        printf("PsychOculusVRCore1-DEBUG: eye %i - swapchain index %i\n", eyeIndex, out_Index);
    */

    if (OVR_FAILURE(ovr_GetTextureSwapChainBufferGL(oculus->hmd, oculus->textureSwapChain[eyeIndex], -1, &texObjectHandle))) {
        ovr_GetLastErrorInfo(&errorInfo);
        if (verbosity > 0) printf("PsychOculusVRCore1-ERROR: eye %i ovr_GetTextureSwapChainBufferGL failed: %s\n", eyeIndex, errorInfo.ErrorString);
        PsychUnlockMutex(&(oculus->presenterLock));
        PsychErrorExitMsg(PsychError_system, "Failed to retrieve next OpenGL texture from swap chain.");
    }

    PsychUnlockMutex(&(oculus->presenterLock));

    // Return texture object handle:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) texObjectHandle);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1CreateMirrorTexture(void)
{
    static char useString[] = "texObjectHandle = PsychOculusVRCore1('CreateMirrorTexture', oculusPtr, width, height);";
    //                         1                                                           1          2      3
    static char synopsisString[] =
    "Create mirror texture for Oculus device 'oculusPtr'.\n"
    "A mirror texture is a auto-updating texture which receives the same image content "
    "which is sent to the HMD for display. It can be used for displaying a debug image "
    "on the regular windowing system, ie., inside a regular Psychtoolbox onscreen window.\n"
    "'width' and 'height' are the width x height of the texture into which "
    "the VR compositor will render the output image which is also sent to the HMD.\n"
    "'texObjectHandle' returns the OpenGL texture handle for accessing the mirror texture.\n";
    static char seeAlsoString[] = "";

    int handle;
    int width, height;
    unsigned int texObjectHandle;
    PsychOculusDevice *oculus;
    ovrMirrorTextureDesc mirrorDesc;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(3));
    PsychErrorExit(PsychRequireNumInputArgs(3));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get texture dimensions:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &width);
    if (width < 1)
        PsychErrorExitMsg(PsychError_user, "Invalid width, smaller than 1 texel!");

    PsychCopyInIntegerArg(3, kPsychArgRequired, &height);
    if (height < 1)
        PsychErrorExitMsg(PsychError_user, "Invalid height, smaller than 1 texel!");

    if (oculus->mirrorTexture)
        PsychErrorExitMsg(PsychError_user, "Tried to create already created mirror texture.");

    // Build OpenGL texture chain descriptor:
    memset(&mirrorDesc, 0, sizeof(mirrorDesc));
    mirrorDesc.Format = OVR_FORMAT_R8G8B8A8_UNORM_SRGB;
    mirrorDesc.MiscFlags = ovrTextureMisc_None;
    mirrorDesc.Width = width;
    mirrorDesc.Height = height;

    // Create mirror texture:
    PsychLockMutex(&(oculus->presenterLock));

    if (OVR_FAILURE(ovr_CreateMirrorTextureGL(oculus->hmd, &mirrorDesc, &oculus->mirrorTexture))) {
        ovr_GetLastErrorInfo(&errorInfo);
        if (verbosity > 0) printf("PsychOculusVRCore1-ERROR: ovr_CreateMirrorTextureGL failed: %s\n", errorInfo.ErrorString);
        PsychUnlockMutex(&(oculus->presenterLock));
        PsychErrorExitMsg(PsychError_system, "Failed to create mirror texture for VR compositor.");
    }

    if (OVR_FAILURE(ovr_GetMirrorTextureBufferGL(oculus->hmd, oculus->mirrorTexture, &texObjectHandle))) {
        ovr_GetLastErrorInfo(&errorInfo);
        if (verbosity > 0) printf("PsychOculusVRCore1-ERROR: ovr_GetMirrorTextureBufferGL failed: %s\n", errorInfo.ErrorString);
        PsychUnlockMutex(&(oculus->presenterLock));
        PsychErrorExitMsg(PsychError_system, "Failed to get OpenGL texture handle of mirror texture for VR compositor.");
    }

    PsychUnlockMutex(&(oculus->presenterLock));

    // Return texture object handle:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) texObjectHandle);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetUndistortionParameters(void)
{
    static char useString[] = "[width, height, viewPx, viewPy, viewPw, viewPh, pptax, pptay, hmdShiftx, hmdShifty, hmdShiftz] = PsychOculusVRCore1('GetUndistortionParameters', oculusPtr, eye [, inputWidth][, inputHeight][, fov]);";
    //                          1      2       3       4       5       6       7      8      9          10         11                                                           1          2      3             4              5
    static char synopsisString[] =
    "Return parameters needed for rendering for output on Oculus device 'oculusPtr'.\n"
    "'eye' which eye to provide the data: 0 = Left, 1 = Right.\n"
    "'inputWidth' = Width of the rendered input image buffer in pixels.\n"
    "'inputHeight' = Height of the rendered input image buffer in pixels.\n"
    "'fov' Optional field of view in degrees, from line of sight: [leftdeg, rightdeg, updeg, downdeg]. "
    "You can pass in the 'fovPort' value returned from PsychOculusVR1('GetFovTextureSize'); "
    "Defaults to whatever has been set for the given eye in the last call to PsychOculusVR1('GetFovTextureSize').\n"
    "\n"
    "Return values:\n"
    "[width, height] = Width and height of client renderbuffers in pixels. Same as the provided 'inputWidth' and 'inputHeight'.\n"
    "[viewPx, viewPy, viewPw, viewPh] Render viewport [x,y,w,h] start (x,y) position, width and height. Mostly useless atm., as "
    "viewPx and viewPy are always zero, viewPw and viewPh are identical to width and height in the current driver design.\n"
    "[pptax pptay] = Pixels per tangens angle at display center in x direction and y direction.\n"
    "[hmdShiftx, hmdShifty, hmdShiftz] = HmdToEyeOffset 3D translation vector. Defines the location of the optical center of the eye "
    "relative to the origin of the local head reference frame, ie. the tracked head position.\n";
    static char seeAlsoString[] = "GetFovTextureSize";

    int handle, eyeIndex;
    PsychOculusDevice *oculus;
    int m, n, p;
    double *fov;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(11));
    PsychErrorExit(PsychCapNumInputArgs(5));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get eye index - left = 0, right = 1:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &eyeIndex);
    if (eyeIndex < 0 || eyeIndex > 1) PsychErrorExitMsg(PsychError_user, "Invalid 'eye' specified. Must be 0 or 1 for left- or right eye.");

    // Get input texture width:
    if (PsychCopyInIntegerArg(3, kPsychArgOptional, &(oculus->texSize[eyeIndex].w))) {
        if (oculus->texSize[eyeIndex].w < 1) PsychErrorExitMsg(PsychError_user, "Invalid 'inputWidth' specified. Must be greater than zero.");
    }

    // Get input texture height:
    if (PsychCopyInIntegerArg(4, kPsychArgOptional, &(oculus->texSize[eyeIndex].h))) {
        if (oculus->texSize[eyeIndex].h < 1) PsychErrorExitMsg(PsychError_user, "Invalid 'inputHeight' specified. Must be greater than zero.");
    }

    // Get optional field of view in degrees in left,right,up,down direction from line of sight:
    if (PsychAllocInDoubleMatArg(5, kPsychArgOptional, &n, &m, &p, &fov)) {
        // Validate and assign:
        if (n * m * p != 4) PsychErrorExitMsg(PsychError_user, "Invalid 'fov' specified. Must be a 4-component vector of form [leftdeg, rightdeg, updeg, downdeg].");
        oculus->ofov[eyeIndex].LeftTan  = (float) tan(deg2rad(fov[0]));
        oculus->ofov[eyeIndex].RightTan = (float) tan(deg2rad(fov[1]));
        oculus->ofov[eyeIndex].UpTan    = (float) tan(deg2rad(fov[2]));
        oculus->ofov[eyeIndex].DownTan  = (float) tan(deg2rad(fov[3]));
    }

    // Return width and height of input texture - Just mirror out what we got:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, oculus->texSize[eyeIndex].w);
    PsychCopyOutDoubleArg(2, kPsychArgOptional, oculus->texSize[eyeIndex].h);

    // Get eye render description for this eye:
    oculus->eyeRenderDesc[eyeIndex] = ovr_GetRenderDesc(oculus->hmd, (ovrEyeType) eyeIndex, oculus->ofov[eyeIndex]);

    if (verbosity > 3) {
        printf("PsychOculusVRCore1-INFO: For HMD %i, eye %i - RenderDescription:\n", handle, eyeIndex);
        printf("PsychOculusVRCore1-INFO: FoV: %f %f %f %f - %f %f %f %f\n", oculus->ofov[eyeIndex].LeftTan, oculus->ofov[eyeIndex].RightTan, oculus->ofov[eyeIndex].UpTan, oculus->ofov[eyeIndex].DownTan, oculus->eyeRenderDesc[eyeIndex].Fov.LeftTan, oculus->eyeRenderDesc[eyeIndex].Fov.RightTan, oculus->eyeRenderDesc[eyeIndex].Fov.UpTan, oculus->eyeRenderDesc[eyeIndex].Fov.DownTan);
        printf("PsychOculusVRCore1-INFO: DistortedViewport: [x,y,w,h] = [%i, %i, %i, %i]\n", oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Pos.x, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Pos.y, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Size.w, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Size.h);
        printf("PsychOculusVRCore1-INFO: PixelsPerTanAngleAtCenter: %f x %f\n", oculus->eyeRenderDesc[eyeIndex].PixelsPerTanAngleAtCenter.x, oculus->eyeRenderDesc[eyeIndex].PixelsPerTanAngleAtCenter.y);
        printf("PsychOculusVRCore1-INFO: HmdToEyeViewOffset: [x,y,z] = [%f, %f, %f]\n", oculus->eyeRenderDesc[eyeIndex].HmdToEyeOffset.x, oculus->eyeRenderDesc[eyeIndex].HmdToEyeOffset.y, oculus->eyeRenderDesc[eyeIndex].HmdToEyeOffset.z);
    }

    // DistortedViewport [x,y,w,h]:
    PsychCopyOutDoubleArg(3, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Pos.x);
    PsychCopyOutDoubleArg(4, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Pos.y);
    PsychCopyOutDoubleArg(5, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Size.w);
    PsychCopyOutDoubleArg(6, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Size.h);

    // PixelsPerTanAngleAtCenter:
    PsychCopyOutDoubleArg(7, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].PixelsPerTanAngleAtCenter.x);
    PsychCopyOutDoubleArg(8, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].PixelsPerTanAngleAtCenter.y);

    // HmdToEyeViewOffset: [x,y,z]:
    PsychCopyOutDoubleArg(9,  kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].HmdToEyeOffset.x);
    PsychCopyOutDoubleArg(10, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].HmdToEyeOffset.y);
    PsychCopyOutDoubleArg(11, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].HmdToEyeOffset.z);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetTrackersState(void)
{
    static char useString[] = "trackers = PsychOculusVRCore1('GetTrackersState', oculusPtr);";
    //                         1                                                 1
    static char synopsisString[] =
    "Return info about all connected trackers for Oculus device 'oculusPtr'.\n\n"
    "'trackers' is an array of structs, with one struct for each connected tracker camera. "
    "Each struct contains the following fields, describing static and dynamic properties "
    "of the corresponding tracking camera:\n\n"
    "'Status': Sum of +4 = The tracking camera has a valid pose, +32 = The camera is present and online.\n"
    "'CameraPose' as vector with position [tx,ty,tz] in meters and orientation quaternion [rx,ry,rz,rw] in radians, "
    "concatenated together as [tx,ty,tz,rx,ry,rz,rw].\n"
    "'LeveledCameraPose' Like 'CameraPose', but aligned to the gravity vector of the world.\n"
    "'CameraFrustumHVFov' Horizontal and vertical field of view of the tracking camera in radians.\n"
    "'CameraFrustumNearFarZInMeters' Near and far limit of the camera view frustum in meters.\n"
    "\n";
    static char seeAlsoString[] = "GetTrackingState";

    PsychGenericScriptType *status;
    const char *FieldNames[] = {"Status", "CameraPose", "LeveledCameraPose", "CameraFrustumHVFov", "CameraFrustumNearFarZInMeters"};
    const int FieldCount = 5;
    PsychGenericScriptType *outMat;
    double *v;

    int handle, trackerCount, i;
    PsychOculusDevice *oculus;
    int StatusFlags = 0;
    ovrTrackerPose trackerPose;
    ovrTrackerDesc trackerDesc;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Number of available tracker cameras:
    trackerCount = ovr_GetTrackerCount(oculus->hmd);

    PsychAllocOutStructArray(1, kPsychArgOptional, trackerCount, FieldCount, FieldNames, &status);

    for (i = 0; i < trackerCount; i++) {
        trackerPose = ovr_GetTrackerPose(oculus->hmd, (unsigned int) i);
        StatusFlags = trackerPose.TrackerFlags & (ovrTracker_Connected | ovrTracker_PoseTracked);

        // Return head and general tracking status flags:
        PsychSetStructArrayDoubleElement("Status", i, StatusFlags, status);

        // Camera pose:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 7, 1, &v, &outMat);
        v[0] = trackerPose.Pose.Position.x;
        v[1] = trackerPose.Pose.Position.y;
        v[2] = trackerPose.Pose.Position.z;

        v[3] = trackerPose.Pose.Orientation.x;
        v[4] = trackerPose.Pose.Orientation.y;
        v[5] = trackerPose.Pose.Orientation.z;
        v[6] = trackerPose.Pose.Orientation.w;
        PsychSetStructArrayNativeElement("CameraPose", i, outMat, status);

        // Camera leveled pose:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 7, 1, &v, &outMat);
        v[0] = trackerPose.LeveledPose.Position.x;
        v[1] = trackerPose.LeveledPose.Position.y;
        v[2] = trackerPose.LeveledPose.Position.z;

        v[3] = trackerPose.LeveledPose.Orientation.x;
        v[4] = trackerPose.LeveledPose.Orientation.y;
        v[5] = trackerPose.LeveledPose.Orientation.z;
        v[6] = trackerPose.LeveledPose.Orientation.w;
        PsychSetStructArrayNativeElement("LeveledCameraPose", i, outMat, status);

        trackerDesc = ovr_GetTrackerDesc(oculus->hmd, (unsigned int) i);

        // Camera frustum HFov and VFov in radians:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
        v[0] = trackerDesc.FrustumHFovInRadians;
        v[1] = trackerDesc.FrustumVFovInRadians;
        PsychSetStructArrayNativeElement("CameraFrustumHVFov", i, outMat, status);

        // Camera frustum near and far clip plane in meters:
        v = NULL;
        PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
        v[0] = trackerDesc.FrustumNearZInMeters;
        v[1] = trackerDesc.FrustumFarZInMeters;
        PsychSetStructArrayNativeElement("CameraFrustumNearFarZInMeters", i, outMat, status);
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetStaticRenderParameters(void)
{
    static char useString[] = "[projL, projR] = PsychOculusVRCore1('GetStaticRenderParameters', oculusPtr [, clipNear=0.01][, clipFar=10000.0]);";
    //                          1      2                                                        1            2                3
    static char synopsisString[] =
    "Retrieve static rendering parameters for Oculus device 'oculusPtr' at current settings.\n"
    "'clipNear' Optional near clipping plane for OpenGL. Defaults to 0.01.\n"
    "'clipFar' Optional far clipping plane for OpenGL. Defaults to 10000.0.\n"
    "\nReturn arguments:\n\n"
    "'projL' is the 4x4 OpenGL projection matrix for the left eye rendering.\n"
    "'projR' is the 4x4 OpenGL projection matrix for the right eye rendering.\n"
    "Please note that projL and projR are usually identical for typical rendering scenarios.\n";
    static char seeAlsoString[] = "";

    int handle;
    PsychOculusDevice *oculus;
    ovrMatrix4f M;
    int i, j;
    double clip_near, clip_far;
    double *outM;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(2));
    PsychErrorExit(PsychCapNumInputArgs(3));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get optional near clipping plane:
    clip_near = 0.01;
    PsychCopyInDoubleArg(2, kPsychArgOptional, &clip_near);

    // Get optional far clipping plane:
    clip_far = 10000.0;
    PsychCopyInDoubleArg(3, kPsychArgOptional, &clip_far);

    // Return left projection matrix as return argument 1:
    M = ovrMatrix4f_Projection(oculus->eyeRenderDesc[0].Fov, (float) clip_near, (float) clip_far, ovrProjection_ClipRangeOpenGL);
    PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 4, 4, 1, &outM);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            *(outM++) = (double) M.M[j][i];

    // Return right projection matrix as return argument 2:
    M = ovrMatrix4f_Projection(oculus->eyeRenderDesc[1].Fov, (float) clip_near, (float) clip_far, ovrProjection_ClipRangeOpenGL);
    PsychAllocOutDoubleMatArg(2, kPsychArgOptional, 4, 4, 1, &outM);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            *(outM++) = (double) M.M[j][i];

    return(PsychError_none);
}


PsychError PSYCHOCULUSVR1StartRender(void)
{
    static char useString[] = "[eyePoseL, eyePoseR] = PsychOculusVRCore1('StartRender', oculusPtr);";
    //                          1         2                                             1
    static char synopsisString[] =
    "Mark start of a new 3D head tracked render cycle for Oculus device 'oculusPtr'.\n"
    "Return values are the vectors which define the two eye cameras positions and orientations "
    "for the left eye and right eye 'eyePoseL' and 'eyePoseR'. The vectors are of form "
    "[tx, ty, tz, rx, ry, rz, rw] - A 3 component 3D position followed by a 4 component rotation "
    "quaternion.\n"
    "\n";
    static char seeAlsoString[] = "GetEyePose GetTrackingState EndFrameRender";

    int handle;
    PsychOculusDevice *oculus;
    ovrVector3f HmdToEyeOffset[2];
    double *outM;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(2));
    PsychErrorExit(PsychCapNumInputArgs(1));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get current eye poses for both eyes:
    HmdToEyeOffset[0] = oculus->eyeRenderDesc[0].HmdToEyeOffset;
    HmdToEyeOffset[1] = oculus->eyeRenderDesc[1].HmdToEyeOffset;

    PsychLockMutex(&(oculus->presenterLock));
    ovr_GetEyePoses(oculus->hmd, 0, FALSE, HmdToEyeOffset, oculus->outEyePoses, &oculus->sensorSampleTime);
    PsychUnlockMutex(&(oculus->presenterLock));

    // Left eye pose as raw data:
    PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 1, 7, 1, &outM);

    // Position (x,y,z):
    outM[0] = oculus->outEyePoses[0].Position.x;
    outM[1] = oculus->outEyePoses[0].Position.y;
    outM[2] = oculus->outEyePoses[0].Position.z;

    // Orientation as a quaternion (x,y,z,w):
    outM[3] = oculus->outEyePoses[0].Orientation.x;
    outM[4] = oculus->outEyePoses[0].Orientation.y;
    outM[5] = oculus->outEyePoses[0].Orientation.z;
    outM[6] = oculus->outEyePoses[0].Orientation.w;

    // Right eye pose as raw data:
    PsychAllocOutDoubleMatArg(2, kPsychArgOptional, 1, 7, 1, &outM);

    // Position (x,y,z):
    outM[0] = oculus->outEyePoses[1].Position.x;
    outM[1] = oculus->outEyePoses[1].Position.y;
    outM[2] = oculus->outEyePoses[1].Position.z;

    // Orientation as a quaternion (x,y,z,w):
    outM[3] = oculus->outEyePoses[1].Orientation.x;
    outM[4] = oculus->outEyePoses[1].Orientation.y;
    outM[5] = oculus->outEyePoses[1].Orientation.z;
    outM[6] = oculus->outEyePoses[1].Orientation.w;

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1EndFrameRender(void)
{
    static char useString[] = "PsychOculusVRCore1('EndFrameRender', oculusPtr, targetPresentTime);";
    //                                                              1          2
    static char synopsisString[] =
    "Mark end of a render cycle for Oculus HMD device 'oculusPtr'.\n\n"
    "'targetPresentTime' the target time for presentation of this frame on the HMD.\n"
    "You usually won't call this function yourself, but Screen('Flip') will "
    "call it automatically for you at the appropriate moment.\n";
    static char seeAlsoString[] = "StartRender PresentFrame";

    int handle;
    double targetPresentTime;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(0));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get expected presentation time for the frame that just finished rendering:
    PsychCopyInDoubleArg(2, kPsychArgRequired, &targetPresentTime);

    // Update the expected present timestamp to the predicted one:
    PsychLockMutex(&(oculus->presenterLock));
    oculus->scheduledPresentExecTime = targetPresentTime;
    PsychUnlockMutex(&(oculus->presenterLock));

    return(PsychError_none);
}

// Execute VR present operation, possibly committing new content to the swapchains:
// Must be called with the presenterLock locked!
// Called by idle presenterThread to keep VR compositor timeout handling from kicking in,
// and directly from PSYCHOCULUSVR1PresentFrame when userspace wants to present new content.
static double PresentExecute(PsychOculusDevice *oculus, psych_bool commitTextures, psych_bool inInit)
{
    int eyeIndex, rc;
    psych_bool success = TRUE;
    double tPredictedOnset = 0;
    ovrResult result;
    ovrLayerEyeFov layer0;
    const ovrLayerHeader *layers[1] = { &layer0.Header };

    // printf("PresentExecute-1 %i = %f\n", commitTextures, ovr_GetTimeInSeconds());

    // Is this a present with updated visual content - called from "userspace" runtime?
    if (commitTextures) {
        // Free capacity in swapchains?
        while ((oculus->needSubmit >= oculus->textureSwapChainLength - 1) &&
               (!oculus->isTracking) && (oculus->presenterThread != (psych_thread) NULL)) {
            // No, and the presenterThread is running and should make free space soon. Wait for it to do its job:
            if ((rc = PsychWaitCondition(&(oculus->presentedSignal), &(oculus->presenterLock)))) {
                // Failed: Log it in a hopefully not too unsafe way:
                printf("PsychOculusVRCore1-ERROR: Waitcondition on presentedSignal trigger failed  [%s].\n", strerror(rc));
            }
        }

        // Commit current target textures to chain:
        for (eyeIndex = 0; eyeIndex < ((oculus->isStereo) ? 2 : 1); eyeIndex++) {
            if (OVR_FAILURE(ovr_CommitTextureSwapChain(oculus->hmd, oculus->textureSwapChain[eyeIndex]))) {
                success = FALSE;
                ovr_GetLastErrorInfo(&errorInfo);
                if (verbosity > 0)
                    printf("PsychOculusVRCore1-ERROR: eye %i ovr_CommitTextureSwapChain() failed: %s\n", eyeIndex, errorInfo.ErrorString);
            }
        }

        if (!success)
            goto present_fail;

        // Mark the need for a ovr_SubmitFrame() call before we can commit textures again:
        oculus->needSubmit++;
    }

    // printf("PresentExecute-2 %i - %f\n", commitTextures, ovr_GetTimeInSeconds());

    // Always submit updated textures when called from presenter thread. If called
    // from usercode 'PresentFrame' only submit if the usercode is in full charge,
    // ie. HMD head tracking is active and driving render timing in a fast closed loop:
    if (inInit || !commitTextures || oculus->isTracking || (oculus->presenterThread == (psych_thread) NULL)) {
        // If HMD tracking is disabled then set sensorSampleTime to "now" - reasonable:
        if (!oculus->isTracking)
            oculus->sensorSampleTime = ovr_GetTimeInSeconds();

        // Setup layer headers:
        layer0.Header.Type = ovrLayerType_EyeFov;
        layer0.Header.Flags = ovrLayerFlag_HighQuality | ovrLayerFlag_TextureOriginAtBottomLeft;

        // Use head locked layer if head tracking is disabled, so our frames stay in a fixed
        // position wrt. the HMD - the HMD simply acts as a strapped on mono- or stereo monitor:
        layer0.Header.Flags |= ((!oculus->isTracking) ? ovrLayerFlag_HeadLocked : 0);

        layer0.ColorTexture[0] = oculus->textureSwapChain[0];
        layer0.ColorTexture[1] = (oculus->isStereo) ? oculus->textureSwapChain[1] : NULL;
        layer0.Viewport[0].Pos.x = 0;
        layer0.Viewport[0].Pos.y = 0;
        layer0.Viewport[0].Size.w = oculus->textureWidth;
        layer0.Viewport[0].Size.h = oculus->textureHeight;
        layer0.Viewport[1].Pos.x = 0;
        layer0.Viewport[1].Pos.y = 0;
        layer0.Viewport[1].Size.w = oculus->textureWidth;
        layer0.Viewport[1].Size.h = oculus->textureHeight;
        layer0.Fov[0] = oculus->ofov[0];
        layer0.Fov[1] = oculus->ofov[1];
        layer0.RenderPose[0] = oculus->outEyePoses[0];
        layer0.RenderPose[1] = oculus->outEyePoses[1];
        layer0.SensorSampleTime = oculus->sensorSampleTime;

        tPredictedOnset = ovr_GetPredictedDisplayTime(oculus->hmd, oculus->frameIndex);
        if (verbosity > 3)
            printf("PRESENT[%i]: Last %f, Next %f, dT = %f msecs\n", commitTextures, oculus->lastPresentExecTime, tPredictedOnset, 1000 * (tPredictedOnset - oculus->lastPresentExecTime));

        // Submit frame to compositor for display at earliest possible time:
        result = ovr_SubmitFrame(oculus->hmd, oculus->frameIndex, NULL, layers, 1);
        if (OVR_FAILURE(result)) {
            success = FALSE;
            ovr_GetLastErrorInfo(&errorInfo);
            if (verbosity > 0)
                printf("PsychOculusVRCore1-ERROR: ovr_SubmitFrame() failed: %s\n", errorInfo.ErrorString);
        }
        else {
            // Update onset time again after successful submit:
            tPredictedOnset = ovr_GetPredictedDisplayTime(oculus->hmd, oculus->frameIndex);

            // Record frameIndex at which new visual content will be committed:
            if (oculus->perfStats.FrameStatsCount > 0) {
                oculus->commitFrameIndex = oculus->perfStats.FrameStats[0].AppFrameIndex + 1;
            }

            oculus->frameIndex++;
            if (oculus->needSubmit > 0)
                oculus->needSubmit--;
        }

        if (result == ovrSuccess_NotVisible) {
            // Submitted frame does not display - a no-op submit:
            if (verbosity > 3)
                printf("PsychOculusVRCore1:WARNING: Frame %i not presented to HMD. HMD removed from subjects head? [Time %f]\n",
                       oculus->frameIndex, tPredictedOnset);
        }

        // Update the lastPresentExecTime timestamp:
        oculus->lastPresentExecTime = tPredictedOnset;

        // printf("PresentExecute-3 %i - %f\n", commitTextures, ovr_GetTimeInSeconds());
        if (success) {
            double tDeadline;

            // Retrieve performance stats, mostly for our timestamping:
            tDeadline = ovr_GetTimeInSeconds() + 0.005;
            do {
                result = ovr_GetPerfStats(oculus->hmd, &oculus->perfStats);
                if (OVR_FAILURE(result)) {
                    ovr_GetLastErrorInfo(&errorInfo);
                    if (verbosity > 0)
                        printf("PsychOculusVRCore1-ERROR: ovr_GetPerfStats() failed: %s\n", errorInfo.ErrorString);
                }
                else if ((verbosity > 3) && (oculus->perfStats.FrameStatsCount > 0)) {
                    printf("PsychOculusVRCore1: DEBUG: Oldest AppFrameIndex %i vs. commitFrameIndex %i\n",
                           oculus->perfStats.FrameStats[oculus->perfStats.FrameStatsCount-1].AppFrameIndex,
                           oculus->commitFrameIndex);
                }

                if (oculus->perfStats.FrameStatsCount < 1)
                    PsychYieldIntervalSeconds(0.001);
            } while ((oculus->frameIndex > 0) && (ovr_GetTimeInSeconds() < tDeadline) &&
                     ((oculus->perfStats.FrameStatsCount < 1) || (oculus->perfStats.FrameStats[oculus->perfStats.FrameStatsCount-1].AppFrameIndex < oculus->commitFrameIndex)));

            if (verbosity > 3) {
                printf("\nPsychOculusVRCore1: [%f msecs headroom] DEBUG: FrameStatsCount=%i, AnyFrameStatsDropped=%i\n",
                       1000 * (tDeadline - ovr_GetTimeInSeconds()), oculus->perfStats.FrameStatsCount,
                       (int) oculus->perfStats.AnyFrameStatsDropped);
            }

            // Are we running on the presenterThread?
            if (!commitTextures && oculus->perfStats.FrameStatsCount) {
                // Yes. Signal end of a new submitFrame -> timestamping cycle to main-thread:
                if ((rc = PsychSignalCondition(&(oculus->presentedSignal)))) {
                    printf("PsychOculusVRCore1-ERROR: PsychSignalCondition() trigger operation failed  [%s].\n", strerror(rc));
                }
            }
        }
    }

    // printf("PresentExecute-4 %i - %f\n", commitTextures, ovr_GetTimeInSeconds());

present_fail:

    return(success ? tPredictedOnset : -1);
}

static void* PresenterThreadMain(void* psychOculusDeviceToCast)
{
    int rc;
    PsychOculusDevice* oculus = (PsychOculusDevice*) psychOculusDeviceToCast;

    // Assign a name to ourselves, for debugging:
    PsychSetThreadName("PsychOculusVR1PresenterThread");

    // VR compositor timeout prevention loop: Repeats infinitely, well, not infinitely,
    // but until we receive a shutdown request and terminate ourselves...
    while (TRUE) {
        // Try to lock, block until available if not available:
        if ((rc = PsychLockMutex(&(oculus->presenterLock)))) {
            // This could potentially kill the runtime, as we're printing from outside the main interpreter thread.
            // Use fprintf() instead of the overloaded printf() (aka mexPrintf()) in the hope that we don't
            // wreak havoc -- maybe it goes to the system log, which should be safer...
            fprintf(stderr, "PsychOculusVRCore1-ERROR: In PresenterThreadMain(): First mutex_lock in init failed  [%s].\n", strerror(rc));
            return(NULL);
        }

        // Shall we terminate? If so break out of our main loop:
        if (oculus->closing) {
            PsychUnlockMutex(&(oculus->presenterLock));
            break;
        }

        // Non-tracked presentation mode for pure monoscopic or stereoscopic presentation
        // without use of head tracking - The HMD is just used as a strapped on stereo or
        // mono monitor, not in a closed-loop setup. Subject will just sit on a chair and
        // likely not move the head much at all. Here we expect that usercode may pass
        // target presentation times 'when' into the Screen('Flip', window, when, ...);
        // command to ask for stimulus presentation in the (potentially far) future, ie.,
        // more than a HMD video refresh duration into the future. We want to prevent the
        // compositor timeout mechanisms and performance optimizations like Asynchronous space
        // warp (ASW) to trigger, as that would totally disrupt the intended experiment timing.
        //
        // Our purpose is to prevent timeouts by executing ovr_SubmitFrame() at least once per
        // HMD refresh interval / compositor duty cycle, with the old/current frame content in
        // the swapchains.
        if ((PresentExecute(oculus, FALSE, FALSE) < 0) && verbosity > 0) {
            fprintf(stderr, "PsychOculusVRCore1-ERROR: In PresenterThreadMain(): PresentExecute() for timeout prevention failed!\n");
        }

        // Unlock and sleep for a frameDuration:
        PsychUnlockMutex(&(oculus->presenterLock));
        PsychYieldIntervalSeconds(oculus->frameDuration);

        // Next dispatch loop iteration...
    }

    // Exit path from thread at thread termination. Go and die peacefully...
    return(NULL);
}

PsychError PSYCHOCULUSVR1PresentFrame(void)
{
    static char useString[] = "[frameTiming, tPredictedOnset, referenceFrameIndex] = PsychOculusVRCore1('PresentFrame', oculusPtr [, doTimestamp=0][, when=0]);";
    //                          1            2                3                                                         1            2                3
    static char synopsisString[] =
    "Present last rendered frame to Oculus HMD device 'oculusPtr'.\n\n"
    "This will commit the current set of 2D textures with new rendered content "
    "to the texture swapchains, for consumption by the VR runtime/compositor. "
    "The swapchains will advance, providing new unused image textures as new "
    "render targets for the next rendering cycle.\n\n"
    "You usually won't call this function yourself, but Screen('Flip') "
    "will call it automatically for you at the appropriate moment.\n\n"
    "Returns a 'frameTiming' struct with information about the presentation "
    "timing for this presentation cycle.\n\n"
    "'doTimestamp' If set to 1, perform timestamping of stimulus onset on the HMD, "
    "or at least try to estimate such onset time. If set to 0, do nothing timestamping-wise.\n\n"
    "'when' If provided, defines the target presentation time, as provided by Screen('Flip', win, when); "
    "a value of zero, or omission, means to present as soon as possible.\n\n"
    "'frameTiming' is an array of structs, where each struct contains presentation "
    "timing and performance info for one presented frame, with the most recently "
    "presented frame in frameTiming(1), older ones at higher indices if the script "
    "falls behind the compositor in Screen('Flip')ping frames in time for the HMD.\n"
    "If successful timestamping happened, 'referenceFrameIndex' will be the index into "
    "the 'frameTiming' array to the struct which corresponds to the just presented frame, "
    "otherwise zero is returned.\n\n"
    "'tPredictedOnset' is another predicted onset time for the just presented frame, based on "
    "very cheap and fast time-stamping which may be much less reliable though.\n\n"
    "The following fields are currently supported in each struct of 'frameTiming':\n\n"
    "HmdVsyncIndex: Vsync counter of the HMD, increments at each refresh.\n"
    "AppFrameIndex: Running count of 'Flip'ped frames. Increments at each Screen('Flip').\n"
    "HMDTime: Time in the HMDs time base in seconds.\n"
    "GetSecsTime: Time in Psychtoolbox standard GetSecs timebase\n"
    "AppMotionToPhotonLatency: Latency between frame submission (or last head tracker query) "
    "and presentation of the associated frame.\n"
    "RecentSensorSampleTime: Baseline for AppMotionToPhotonLatency, if app can keep up with compositor.\n"
    "StimulusOnsetTime: Estimated visual stimulus onset time for the frame - the mid-point of scanout!\n"
    "VBlankTime: Screen('Flip') vblTime equivalent.\n"
    "AppQueueAheadTime: How much ahead the application has to queue a frame for presentation. In some sense "
    "a measure of compositor lag. 0 msecs = No lag, presenting at next VBLANK boundary. 1 refresh duration = "
    "1 frame lag, etc.\n\n"
    "'tPredictedOnset' Predicted onset time for the just submitted frame, according to compositor.\n"
    "\n\n";
    static char seeAlsoString[] = "EndFrameRender";

    PsychGenericScriptType *frameT;
    const char *FieldNames[] = {"HmdVsyncIndex", "AppFrameIndex", "HMDTime", "GetSecsTime", "AppMotionToPhotonLatency",
                                "RecentSensorSampleTime", "StimulusOnsetTime", "VBlankTime", "AppQueueAheadTime" };
    const int FieldCount = 9;

    int handle, rc;
    int doTimestamp = 0;
    int referenceFrameIndex = -1;
    double tNow, tHMD, tStimOnset, tVBL, tPredictedOnset;
    double tWhen;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(3));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get optional timestamping flag:
    PsychCopyInIntegerArg(2, kPsychArgOptional, &doTimestamp);

    // Get optional presentation target time:
    tWhen = 0;
    PsychCopyInDoubleArg(3, kPsychArgOptional, &tWhen);

    // Single-threaded operation, all done on this thread?
    if (!oculus->multiThreaded) {
        // Yes. Go into a check -> submit -> check loop until target time is
        // exceeded, then present an updated frame:
        while ((tPredictedOnset = ovr_GetPredictedDisplayTime(oculus->hmd, oculus->frameIndex)) < tWhen + oculus->frameDuration / 2) {
            if (verbosity > 4)
                printf("PsychOculusVRCore1:PresentFrame-DEBUG: In busy-present-wait loop: pred %f < %f tWhen\n",
                       tPredictedOnset, tWhen);
            PresentExecute(oculus, FALSE, FALSE);
        }
    }

    // Execute the present operation with the presenterThread locked out.
    // Invalidate any scheduledPresentExecTime after such a Present:
    PsychLockMutex(&(oculus->presenterLock));
    oculus->scheduledPresentExecTime = tWhen;
    tPredictedOnset = PresentExecute(oculus, TRUE, FALSE);
    oculus->scheduledPresentExecTime = -DBL_MAX;
    PsychUnlockMutex(&(oculus->presenterLock));

    if ((tPredictedOnset < 0) && (verbosity > 0))
        printf("PsychOculusVRCore1-ERROR: Failed to present new frames to VR compositor.\n");

    PsychGetAdjustedPrecisionTimerSeconds(&tNow);
    tHMD = ovr_GetTimeInSeconds();

    if (doTimestamp) {
        int i;
        uint32_t timestampedFrameIndex = 0;

        PsychLockMutex(&(oculus->presenterLock));
        while (timestampedFrameIndex == 0) {
            for (i = oculus->perfStats.FrameStatsCount - 1; i >= 0; i--) {
                if (oculus->perfStats.FrameStats[i].AppFrameIndex >= oculus->commitFrameIndex) {
                    timestampedFrameIndex = oculus->perfStats.FrameStats[i].AppFrameIndex;
                    referenceFrameIndex = i;
                    break;
                }
            }

            if (timestampedFrameIndex == 0) {
                if (oculus->presenterThread != (psych_thread) NULL) {
                    // Multi-Threaded: Wait efficiently for signalling of new frame submitted by presenter thread:
                    if ((rc = PsychWaitCondition(&(oculus->presentedSignal), &(oculus->presenterLock)))) {
                        // Failed: Log it in a hopefully not too unsafe way:
                        printf("PsychOculusVRCore1-ERROR: Waitcondition II on presentedSignal trigger failed  [%s].\n", strerror(rc));
                    }
                }
                else {
                    // Not multi-threaded. If we ain't got no result yet, we won't
                    // get any in the future:
                    break;
                }
            }
        }

        if (oculus->perfStats.FrameStatsCount > 0) {
            PsychAllocOutStructArray(1, kPsychArgOptional, oculus->perfStats.FrameStatsCount, FieldCount, FieldNames, &frameT);
            for (i = 0; i < oculus->perfStats.FrameStatsCount; i++) {
                // HMD Vsync counter:
                PsychSetStructArrayDoubleElement("HmdVsyncIndex", i, oculus->perfStats.FrameStats[i].HmdVsyncIndex, frameT);

                // Our sbc counter:
                PsychSetStructArrayDoubleElement("AppFrameIndex", i, oculus->perfStats.FrameStats[i].AppFrameIndex, frameT);

                // Clock-Sync PTB timebase vs. Oculus timebase:
                PsychGetAdjustedPrecisionTimerSeconds(&tNow);
                tHMD = ovr_GetTimeInSeconds();
                PsychSetStructArrayDoubleElement("HMDTime", i, tHMD, frameT);
                PsychSetStructArrayDoubleElement("GetSecsTime", i, tNow, frameT);

                // Time between oculus->sensorSampleTime and visual onset (video frame midpoint of scanout):
                PsychSetStructArrayDoubleElement("AppMotionToPhotonLatency", i, oculus->perfStats.FrameStats[i].AppMotionToPhotonLatency, frameT);
                PsychSetStructArrayDoubleElement("RecentSensorSampleTime", i, oculus->sensorSampleTime, frameT);

                // Compute absolute stimulus onset (mid-point), remap to PTB GetSecs time:
                tStimOnset = (tNow - tHMD)  + oculus->sensorSampleTime + oculus->perfStats.FrameStats[i].AppMotionToPhotonLatency;
                PsychSetStructArrayDoubleElement("StimulusOnsetTime", i, tStimOnset, frameT);

                // Compute virtual start of VBLANK time as stimulus onset - half a HMD video refresh duration:
                tVBL = tStimOnset - 0.5 * oculus->frameDuration;
                PsychSetStructArrayDoubleElement("VBlankTime", i, tVBL, frameT);

                // Queue ahead for application. Citation from the v 1.16 SDK docs:
                // "Amount of queue-ahead in seconds provided to the app based on performance and overlap of CPU & GPU utilization
                // A value of 0.0 would mean the CPU & GPU workload is being completed in 1 frame's worth of time, while
                // 11 ms (on the CV1) of queue ahead would indicate that the app's CPU workload for the next frame is
                // overlapping the app's GPU workload for the current frame."
                PsychSetStructArrayDoubleElement("AppQueueAheadTime", i, oculus->perfStats.FrameStats[i].AppQueueAheadTime, frameT);
            }
        }

        PsychUnlockMutex(&(oculus->presenterLock));
    }
    else {
        // Return empty 'frameTiming' info struct array, so our calling scripting environment does not bomb out:
        PsychAllocOutStructArray(1, kPsychArgOptional, 0, FieldCount, FieldNames, &frameT);
    }

    // Copy out predicted onset time for the just emitted frame:
    tPredictedOnset = tPredictedOnset + (tNow - tHMD) - 0.5 * oculus->frameDuration;
    PsychCopyOutDoubleArg(2, kPsychArgOptional, tPredictedOnset);

    PsychCopyOutDoubleArg(3, kPsychArgOptional, (double) referenceFrameIndex + 1);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetEyePose(void)
{
    static char useString[] = "[eyePose, eyeIndex] = PsychOculusVRCore1('GetEyePose', oculusPtr, renderPass);";
    //                          1        2                                            1          2
    static char synopsisString[] =
    "Return current predicted pose vector for an eye for Oculus device 'oculusPtr'.\n"
    "'renderPass' is the view render pass for which to provide the data: 0 = First pass, 1 = Second pass.\n"
    "Return value is the vector 'eyePose' which defines the position and orientation for the eye corresponding "
    "to the requested renderPass ie. 'eyePose' = [posX, posY, posZ, rotX, rotY, rotZ, rotW].\n"
    "The second return value is the 'eyeIndex', the index of the eye whose view should be rendered. This would "
    "be 0 for left eye, and 1 for right eye, and could be used to select the target render view via, e.g.,\n"
    "Screen('SelectStereoDrawBuffer', window, eyeIndex);\n"
    "Which 'eyeIndex' corresponds to the first or second 'renderPass', ie., if the left eye should be rendered "
    "first, or if the right eye should be rendered first, depends on the visual scanning order of the HMD "
    "display panel - is it top to bottom, left to right, or right to left? This function provides that optimized "
    "mapping. Using this function to query the parameters for render setup of an eye can provide a bit more "
    "accuracy in rendering, at the expense of more complex usercode. If you want the more simple solution, just "
    "omit this function and use the eyePose's provided for both eyes by the 'StartRender' method.\n";

    static char seeAlsoString[] = "StartRender";

    int handle, renderPass;
    PsychOculusDevice *oculus;
    ovrVector3f HmdToEyeOffset[2];
    ovrEyeType eye;
    double *outM;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(2));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get renderPass:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &renderPass);
    if (renderPass < 0 || renderPass > 1) PsychErrorExitMsg(PsychError_user, "Invalid 'renderPass' specified. Must be 0 or 1 for first or second pass.");

    // Get eye pose:
    eye = renderPass;
    HmdToEyeOffset[0] = oculus->eyeRenderDesc[0].HmdToEyeOffset;
    HmdToEyeOffset[1] = oculus->eyeRenderDesc[1].HmdToEyeOffset;

    PsychLockMutex(&(oculus->presenterLock));
    ovr_GetEyePoses(oculus->hmd, 0, FALSE, HmdToEyeOffset, oculus->outEyePoses, &oculus->sensorSampleTime);
    PsychUnlockMutex(&(oculus->presenterLock));

    // Eye pose as raw data:
    PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 1, 7, 1, &outM);

    // Position (x,y,z):
    outM[0] = oculus->outEyePoses[eye].Position.x;
    outM[1] = oculus->outEyePoses[eye].Position.y;
    outM[2] = oculus->outEyePoses[eye].Position.z;

    // Orientation as a quaternion (x,y,z,w):
    outM[3] = oculus->outEyePoses[eye].Orientation.x;
    outM[4] = oculus->outEyePoses[eye].Orientation.y;
    outM[5] = oculus->outEyePoses[eye].Orientation.z;
    outM[6] = oculus->outEyePoses[eye].Orientation.w;

    // Copy out preferred eye render order for info:
    PsychCopyOutDoubleArg(2, kPsychArgOptional, (double) eye);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1SetHUDState(void)
{
    static char useString[] = "PsychOculusVRCore1('SetHUDState', oculusPtr , mode);";
    //                                                           1           2
    static char synopsisString[] =
    "Set mode of operation for the performance HUD for Oculus device 'oculusPtr'.\n"
    "'mode' can be one of the following settings:\n"
    "0 = Head up display HUD off ie. invisible.\n"
    "1 = HUD shows performance summary.\n"
    "2 = HUD shows latency timing.\n"
    "3 = HUD shows application render timing.\n"
    "4 = HUD shows VR compositor render timing.\n"
    "5 = HUD shows Version information of different components.\n"
    "6 = HUD shows Asynchronous time warp (ASW) stats.\n"
    "Higher numbers may do something useful in the future.\n";

    static char seeAlsoString[] = "";
    int handle, mode;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(0));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get the new performance HUD mode:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &mode);
    if (mode < 0)
        PsychErrorExitMsg(PsychError_user, "Invalid performance HUD mode passed. Must be >= 0.");

    // Apply new mode:
    ovr_SetInt(oculus->hmd, OVR_PERF_HUD_MODE, mode);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1GetPerformanceStats(void)
{
    static char useString[] = "[adaptiveGpuPerformanceScale, frameStats, anyFrameStatsDropped, aswIsAvailable] = PsychOculusVRCore1('GetPerformanceStats', oculusPtr);";
    //                          1                            2           3                     4                                                           1
    static char synopsisString[] =
    "Query performance statistics for HMD device 'oculusPtr'.\n\n"
    "Return values:\n\n"
    "'adaptiveGpuPerformanceScale' is an edge-filtered value that a caller can use to adjust "
    "the graphics quality of the application to keep the GPU utilization in check. The value "
    "is calculated as:\n\n"
    "   (desired_GPU_utilization / current_GPU_utilization)\n\n"
    "As such, when this value is 1.0, the GPU is doing the right amount of work for the app. "
    "Lower values mean the app needs to pull back on the GPU utilization.\n"
    "If the app is going to directly drive render-target resolution using this value, then "
    "be sure to take the square-root of the value before scaling the resolution with it. "
    "Changing render target resolutions however is one of the many things an app can do "
    "increase or decrease the amount of GPU utilization.\n"
    "Since AdaptiveGpuPerformanceScale is edge-filtered and does not change rapidly "
    "(i.e. reports non-1.0 values once every couple of seconds) the app can make the "
    "necessary adjustments and then keep watching the value to see if it has been satisfied.\n\n"
    "'aswIsAvailable' If 1 then Asynchronous Space Warp (ASW) is supported on your system, 0 if not.\n\n"
    "'anyFrameStatsDropped' If 1 then the 'frameStats' struct array does not contain complete per frame "
    "performance statistics, because you did not call this function frequently enough. In general, calling "
    "at the native video refresh rate of the HMD is required for best results.\n\n"
    "'frameStats' is a struct array with recorded performance statistics, going some frames backwards in time.\n"
    "The number of frames this goes back is implementation dependent. Array slot 1 contains the stats of the "
    "most recent frame, slot 2 the 2nd oldest etc. Each struct contains many fields with info about the frame "
    "presentation latency, cpu and gpu utilization, VR compositor performance, ASW state and performance etc. Consult Oculus "
    "SDK online documentation to learn about the meaning of the fields.\n";

    static char seeAlsoString[] = "";

    PsychGenericScriptType *frameT;
    const char *FieldNames[] = {"HmdVsyncIndex", "AppFrameIndex", "AppDroppedFrameCount", "AppMotionToPhotonLatency", "AppQueueAheadTime", "AppCpuElapsedTime",
                                "AppGpuElapsedTime", "CompositorFrameIndex", "CompositorDroppedFrameCount", "CompositorLatency",
                                "CompositorCpuElapsedTime", "CompositorGpuElapsedTime", "CompositorCpuStartToGpuEndElapsedTime",
                                "CompositorGpuEndToVsyncElapsedTime", "AswIsActive", "AswActivatedToggleCount",
                                "AswPresentedFrameCount", "AswFailedFrameCount" };
    const int FieldCount = 18;

    int handle, i;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(4));
    PsychErrorExit(PsychCapNumInputArgs(1));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Global parameters:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) oculus->perfStats.AdaptiveGpuPerformanceScale);
    PsychCopyOutDoubleArg(3, kPsychArgOptional, (double) oculus->perfStats.AnyFrameStatsDropped);
    PsychCopyOutDoubleArg(4, kPsychArgOptional, (double) oculus->perfStats.AswIsAvailable);

    // Struct array with per frame perfStats for the last couple of frames:
    PsychAllocOutStructArray(2, kPsychArgOptional, oculus->perfStats.FrameStatsCount, FieldCount, FieldNames, &frameT);
    for (i = 0; i < oculus->perfStats.FrameStatsCount; i++) {
        PsychSetStructArrayDoubleElement("HmdVsyncIndex", i, oculus->perfStats.FrameStats[i].HmdVsyncIndex, frameT);
        PsychSetStructArrayDoubleElement("AppFrameIndex", i, oculus->perfStats.FrameStats[i].AppFrameIndex, frameT);
        PsychSetStructArrayDoubleElement("AppDroppedFrameCount", i, oculus->perfStats.FrameStats[i].AppDroppedFrameCount, frameT);
        PsychSetStructArrayDoubleElement("AppMotionToPhotonLatency", i, oculus->perfStats.FrameStats[i].AppMotionToPhotonLatency, frameT);
        PsychSetStructArrayDoubleElement("AppQueueAheadTime", i, oculus->perfStats.FrameStats[i].AppQueueAheadTime, frameT);
        PsychSetStructArrayDoubleElement("AppCpuElapsedTime", i, oculus->perfStats.FrameStats[i].AppCpuElapsedTime, frameT);
        PsychSetStructArrayDoubleElement("AppGpuElapsedTime", i, oculus->perfStats.FrameStats[i].AppGpuElapsedTime, frameT);
        PsychSetStructArrayDoubleElement("CompositorFrameIndex", i, oculus->perfStats.FrameStats[i].CompositorFrameIndex, frameT);
        PsychSetStructArrayDoubleElement("CompositorDroppedFrameCount", i, oculus->perfStats.FrameStats[i].CompositorDroppedFrameCount, frameT);
        PsychSetStructArrayDoubleElement("CompositorLatency", i, oculus->perfStats.FrameStats[i].CompositorLatency, frameT);
        PsychSetStructArrayDoubleElement("CompositorCpuElapsedTime", i, oculus->perfStats.FrameStats[i].CompositorCpuElapsedTime, frameT);
        PsychSetStructArrayDoubleElement("CompositorGpuElapsedTime", i, oculus->perfStats.FrameStats[i].CompositorGpuElapsedTime, frameT);
        PsychSetStructArrayDoubleElement("CompositorCpuStartToGpuEndElapsedTime", i, oculus->perfStats.FrameStats[i].CompositorCpuStartToGpuEndElapsedTime, frameT);
        PsychSetStructArrayDoubleElement("CompositorGpuEndToVsyncElapsedTime", i, oculus->perfStats.FrameStats[i].CompositorGpuEndToVsyncElapsedTime, frameT);
        PsychSetStructArrayDoubleElement("AswIsActive", i, oculus->perfStats.FrameStats[i].AswIsActive, frameT);
        PsychSetStructArrayDoubleElement("AswActivatedToggleCount", i, oculus->perfStats.FrameStats[i].AswActivatedToggleCount, frameT);
        PsychSetStructArrayDoubleElement("AswPresentedFrameCount", i, oculus->perfStats.FrameStats[i].AswPresentedFrameCount, frameT);
        PsychSetStructArrayDoubleElement("AswFailedFrameCount", i, oculus->perfStats.FrameStats[i].AswFailedFrameCount, frameT);
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1HapticPulse(void)
{
    static char useString[] = "pulseEndTime = PsychOculusVRCore1('HapticPulse', oculusPtr, controllerType [, duration=2.5][, freq=1.0][, amplitude=1.0]);";
    //                         1                                                1          2                 3               4           5
    static char synopsisString[] =
    "Execute a haptic feedback pulse on controller 'controllerType' associated with Oculus device 'oculusPtr'.\n\n"
    "'duration' is the duration of the pulse in seconds. If omitted, the function returns immediately and "
    "the maximum duration of 2.5 seconds is used, unless you call the function again with 'freq' = 0, to cancel the "
    "currently active pulse. Otherwise, if a 'duration' other than 2.5 seconds is specified, the function executes "
    "for the specified duration, and may block execution of your script for that time span on some controller types, "
    "returning the absolute time when the pulse is expected to end.\n\n"
    "'freq' Frequency of the vibration in normalized 0.0 - 1.0 range. Currently only values 0, 0.5 and 1.0 "
    "will be used. 0 = Disable ongoing pulse immediately. Values other than 0, 0.5, or 1 will be clamped to the nearest value.\n\n"
    "'amplitude' Normalized amplitude in range 0.0 - 1.0\n\n"
    "The return argument 'pulseEndTime' contains the absolute time in seconds when the pulse is expected "
    "to end, as estimated at the time of calling the function. The precision and accuracy of pulse timing "
    "is not known.\n";

    static char seeAlsoString[] = "";
    int handle, controllerType;
    PsychOculusDevice *oculus;
    double duration, freq, amplitude, pulseEndTime;
    ovrResult result;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(5));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get the new performance HUD mode:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &controllerType);

    // Duration:
    duration = 2.5;
    PsychCopyInDoubleArg(3, kPsychArgOptional, &duration);
    if (duration < 0)
        PsychErrorExitMsg(PsychError_user, "Invalid negative 'duration' in seconds specified. Must be positive.");

    if ((duration > 2.5) && (verbosity > 1))
        printf("PsychOculusVRCore1-WARNING: 'HapticPulse' of %f seconds duration requested, but currently duration is limited to a maximum of 2.5 seconds. Clamping...\n",
                duration);

    if (duration > 2.5)
        duration = 2.5;

    freq = 1.0;
    PsychCopyInDoubleArg(4, kPsychArgOptional, &freq);
    if (freq < 0.0 || freq > 1.0)
        PsychErrorExitMsg(PsychError_user, "Invalid 'freq' frequency specified. Must be in range [0.0 ; 1.0].");

    amplitude = 1.0;
    PsychCopyInDoubleArg(5, kPsychArgOptional, &amplitude);
    if (amplitude < 0.0 || amplitude > 1.0)
        PsychErrorExitMsg(PsychError_user, "Invalid 'amplitude' specified. Must be in range [0.0 ; 1.0].");

    // Execute pulse:
    result = ovr_SetControllerVibration(oculus->hmd, (ovrControllerType) controllerType, (float) freq, (float) amplitude);
    if (OVR_FAILURE(result)) {
        ovr_GetLastErrorInfo(&errorInfo);
        if (verbosity > 0)
            printf("PsychOculusVRCore1-ERROR: ovr_SetControllerVibration() failed: %s\n", errorInfo.ErrorString);
        PsychErrorExitMsg(PsychError_system, "Failed to initiate haptic feedback pulse.");
    }

    if ((result == ovrSuccess_DeviceUnavailable) && (verbosity > 1))
        printf("PsychOculusVRCore1-WARNING: 'HapticPulse' will go nowhere, as suitable controller of type %i is not connected.\n", controllerType);

    if ((result == ovrSuccess) && (verbosity > 3))
        printf("PsychOculusVRCore1-INFO: 'HapticPulse' of duration %f secs, freq %f, amplitude %f for controller of type %i started.\n",
               duration, freq, amplitude, controllerType);

    // Predict "off" time:
    PsychGetAdjustedPrecisionTimerSeconds(&pulseEndTime);
    pulseEndTime += duration;
    PsychCopyOutDoubleArg(1, kPsychArgOptional, pulseEndTime);

    // Pulse of predefined duration requested?
    if ((freq != 0) && (duration < 2.5)) {
        // Yes. Wait until expected end time, then stop the pulse:
        PsychWaitUntilSeconds(pulseEndTime);
        result = ovr_SetControllerVibration(oculus->hmd, (ovrControllerType) controllerType, (float) 0, (float) 0);
        if (OVR_FAILURE(result)) {
            ovr_GetLastErrorInfo(&errorInfo);
            if (verbosity > 0)
                printf("PsychOculusVRCore1-ERROR: ovr_SetControllerVibration() failed to stop vibration: %s\n", errorInfo.ErrorString);
            PsychErrorExitMsg(PsychError_system, "Failed to stop haptic feedback pulse.");
        }
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1FloatsProperty(void)
{
    static char useString[] = "oldValues = PsychOculusVRCore1('FloatsProperty', oculusPtr, propertyName [, newValues]);";
    //                         1                                                1          2               3
    static char synopsisString[] =
        "Return or set a floating point value vector property associated with Oculus device 'oculusPtr'.\n\n"
        "'propertyName' String with name of the floats property to set/get. See OVR.KEY_XXX constants.\n\n"
        "'newValues' Optional vector with new floating point values to assign.\n\n"
        "Returns a vector 'oldValues' with the current (or old) floating point values associated with the property.\n"
        "Currently up to 100 elements can be retrieved or set for a property vector.\n";
    static char seeAlsoString[] = "";

    int handle, i, m, n, p, count;
    PsychOculusDevice *oculus;
    char *propertyName;
    float fvalues[100];
    double *values;
    psych_bool isscalar = FALSE;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(3));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Property name string:
    PsychAllocInCharArg(2, kPsychArgRequired, &propertyName);
    if (strlen(propertyName) == 0)
        PsychErrorExitMsg(PsychError_user, "Invalid, empty 'propertyName' provided.");

    // Get current / old values in fvalues, and their count in count, at most 100 elements:
    count = ovr_GetFloatArray(oculus->hmd, (const char*) propertyName, fvalues, 100);

    // Copy them out to userspace:
    PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 1, count, 1, &values);
    for (i = 0; i < count; i++)
        *(values++) = (double) fvalues[i];

    // Optional new values provided?
    if (PsychAllocInDoubleMatArg(3, kPsychArgOptional, &m, &n, &p, &values)) {
        n = m * n * p;
        if (n > 100)
            PsychErrorExitMsg(PsychError_user, "Too many 'newValues' provided for property. Can not take more than 100.");

        if ((count > 0) && (n != count) && (verbosity > 1))
            printf("PsychOculusVRCore1-WARNING: 'FloatsProperty' '%s' had %i values, but now trying to assign different number %i. Possible problem?\n",
                   propertyName, count, n);

        for (i = 0; i < n; i++)
            fvalues[i] = (float) values[i];

        if (!ovr_SetFloatArray(oculus->hmd, (const char*) propertyName, fvalues, n)) {
            if (verbosity > 0)
                printf("PsychOculusVRCore1-ERROR: 'FloatsProperty' '%s' can not accept new values, as it is read-only!\n", propertyName);
            PsychErrorExitMsg(PsychError_user, "Tried to set new values for a read-only property");
        }
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1FloatProperty(void)
{
    static char useString[] = "oldValues = PsychOculusVRCore1('FloatProperty', oculusPtr, propertyName [, newValue]);";
    //                         1                                               1          2               3
    static char synopsisString[] =
        "Return or set a floating point value property associated with Oculus device 'oculusPtr'.\n\n"
        "'propertyName' String with name of the float property to set/get. See OVR.KEY_XXX constants.\n\n"
        "'newValue' Optional new floating point value to assign.\n\n"
        "Returns a scalar 'oldValue' with the current (or old) floating point value associated with the property.\n";
    static char seeAlsoString[] = "";

    int handle;
    PsychOculusDevice *oculus;
    char *propertyName;
    float fvalue;
    double value;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(3));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Property name string:
    PsychAllocInCharArg(2, kPsychArgRequired, &propertyName);
    if (strlen(propertyName) == 0)
        PsychErrorExitMsg(PsychError_user, "Invalid, empty 'propertyName' provided.");

    // Get float property:
    fvalue = ovr_GetFloat(oculus->hmd, (const char*) propertyName, FLT_MAX);

    // Copy it out to userspace:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) fvalue);

    // Optional new values provided?
    if (PsychCopyInDoubleArg(3, kPsychArgOptional, &value)) {
        if (!ovr_SetFloat(oculus->hmd, (const char*) propertyName, (float) value)) {
            if (verbosity > 0)
                printf("PsychOculusVRCore1-ERROR: 'FloatProperty' '%s' can not accept new value, as it is read-only!\n", propertyName);
            PsychErrorExitMsg(PsychError_user, "Tried to set new value for a read-only property");
        }
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVR1StringProperty(void)
{
    static char useString[] = "oldString = PsychOculusVRCore1('StringProperty', oculusPtr, propertyName [, defaultString][, newString]);";
    //                         1                                                1          2               3                4
    static char synopsisString[] =
        "Return or set a string property associated with Oculus device 'oculusPtr'.\n\n"
        "'propertyName' String with name of the string property to set/get. See OVR.KEY_XXX constants.\n\n"
        "'defaultString' A string to return if string 'propertyName' is not defined yet.\n\n"
        "'newString' Optional new text string to assign.\n\n"
        "Returns a string 'oldString' with the current (or old, or 'defaultString') string associated with the property.\n";
    static char seeAlsoString[] = "";

    int handle;
    PsychOculusDevice *oculus;
    char *propertyName, *defaultVal;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(4));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Property name string:
    PsychAllocInCharArg(2, kPsychArgRequired, &propertyName);
    if (strlen(propertyName) == 0)
        PsychErrorExitMsg(PsychError_user, "Invalid, empty 'propertyName' provided.");

    // Optional default string:
    if (!PsychAllocInCharArg(3, kPsychArgOptional, &defaultVal))
        defaultVal = NULL;

    // Return current / old string:
    PsychCopyOutCharArg(1, kPsychArgOptional, ovr_GetString(oculus->hmd, (const char*) propertyName, (const char*) defaultVal));

    // Optional new values provided?
    if (PsychAllocInCharArg(4, kPsychArgOptional, &defaultVal)) {
        if (!ovr_SetString(oculus->hmd, (const char*) propertyName, (const char*) defaultVal)) {
            if (verbosity > 0)
                printf("PsychOculusVRCore1-ERROR: 'StringProperty' '%s' can not accept new value '%s', as it is read-only!\n", propertyName, defaultVal);
            PsychErrorExitMsg(PsychError_user, "Tried to set new values for a read-only string property");
        }
    }

    return(PsychError_none);
}
