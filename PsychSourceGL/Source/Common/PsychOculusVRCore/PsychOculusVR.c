/*
 * PsychToolbox/Source/Common/PsychOculusVRCore/PsychOculusVR.c
 *
 * PROJECTS: PsychOculusVRCore only.
 *
 * AUTHORS:
 *
 * mario.kleiner.de@gmail.com   mk
 *
 * PLATFORMS:   All.
 *
 * HISTORY:
 *
 * 1.09.2015   mk      Created.
 *
 * DESCRIPTION:
 *
 * A Psychtoolbox driver for the Oculus VR virtual reality
 * head sets. The initial version will support the Rift DK2,
 * and possibly the old Rift DK1, although that hasn't been
 * tested.
 *
 */

#include "PsychOculusVR.h"

// Need this _USE_MATH_DEFINES so MSVC 2010 knows M_PI
#define _USE_MATH_DEFINES
#include <math.h>

// Includes from Oculus SDK 0.5:
#include "OVR_CAPI.h"

#ifdef  __cplusplus
// Only needed for C++ math support routines in OVR:: namespace.
#include "Extras/OVR_Math.h"
#endif

// Number of maximum simultaneously open VR devices:
#define MAX_PSYCH_OCULUS_DEVS 10
#define MAX_SYNOPSIS_STRINGS 40

//declare variables local to this file.
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

// Our device record:
typedef struct PsychOculusDevice {
    ovrHmd hmd;
    psych_bool isTracking;
    ovrSizei texSize[2];
    ovrFovPort ofov[2];
    ovrEyeRenderDesc eyeRenderDesc[2];
    ovrDistortionMesh eyeDistortionMesh[2];
    ovrVector2f UVScaleOffset[2][2];
    ovrMatrix4f timeWarpMatrices[2];
    ovrPosef headPose[2];
    ovrFrameTiming frameTiming;
    uint32_t frameIndex;
    ovrPosef outEyePoses[2];
    unsigned char rgbColorOut[3];
    psych_bool latencyTestActive;
} PsychOculusDevice;

PsychOculusDevice oculusdevices[MAX_PSYCH_OCULUS_DEVS];
static int available_devices = 0;
static unsigned int devicecount = 0;
static int verbosity = 3;
static psych_bool initialized = FALSE;

void InitializeSynopsis(void)
{
    int i = 0;
    const char **synopsis = synopsisSYNOPSIS;

    synopsis[i++] = "PsychOculusVRCore - A Psychtoolbox driver for Oculus VR hardware.\n";
    synopsis[i++] = "This driver allows to control Oculus Rift DK1/DK2 and hopefully in the future also future Oculus devices.\n";
    synopsis[i++] = "The PsychOculusVRCore driver is licensed to you under the terms of the MIT license.";
    synopsis[i++] = "See 'help License.txt' in the Psychtoolbox root folder for more details.\n";
    synopsis[i++] = "\n";
    synopsis[i++] = "The driver currently requires the Oculus VR runtime version 0.5.0.1 to work.\n";
    synopsis[i++] = "\n";
    synopsis[i++] = "Usage:";
    synopsis[i++] = "\n";
    synopsis[i++] = "Functions used by regular user scripts:\n";
    synopsis[i++] = "\n";
    synopsis[i++] = "oldVerbosity = PsychOculusVRCore('Verbosity' [, verbosity]);";
    synopsis[i++] = "numHMDs = PsychOculusVRCore('GetCount');";
    synopsis[i++] = "[oculusPtr, modelName] = PsychOculusVRCore('Open' [, deviceIndex=0]);";
    synopsis[i++] = "PsychOculusVRCore('Close' [, oculusPtr]);";
    synopsis[i++] = "showHSW = PsychOculusVRCore('GetHSWState', oculusPtr [, dismiss]);";
    synopsis[i++] = "oldPersistence = PsychOculusVRCore('SetLowPersistence', oculusPtr [, lowPersistence]);";
    synopsis[i++] = "oldDynamicPrediction = PsychOculusVRCore('SetDynamicPrediction', oculusPtr [, dynamicPrediction]);";
    synopsis[i++] = "PsychOculusVRCore('Start', oculusPtr);";
    synopsis[i++] = "PsychOculusVRCore('Stop', oculusPtr);";
    synopsis[i++] = "state = PsychOculusVRCore('GetTrackingState', oculusPtr [, predictionTime=0]);";
    synopsis[i++] = "[projL, projR] = PsychOculusVRCore('GetStaticRenderParameters', oculusPtr [, clipNear=0.01][, clipFar=10000.0]);";
    synopsis[i++] = "[eyePoseL, eyePoseR, tracked, frameTiming] = PsychOculusVRCore('StartRender', oculusPtr);";
    synopsis[i++] = "[eyePose, eyeIndex] = PsychOculusVRCore('GetEyePose', oculusPtr, renderPass);";
    synopsis[i++] = "\n";
    synopsis[i++] = "Functions usually only used internally by Psychtoolbox:\n";
    synopsis[i++] = "\n";
    synopsis[i++] = "[width, height, fovPort] = PsychOculusVRCore('GetFovTextureSize', oculusPtr, eye [, fov=[HMDRecommended]][, pixelsPerDisplay=1]);";
    synopsis[i++] = "[width, height, viewPx, viewPy, viewPw, viewPh, pptax, pptay, hmdShiftx, hmdShifty, hmdShiftz, meshVertices, meshIndices, uvScaleX, uvScaleY, uvOffsetX, uvOffsetY] = PsychOculusVRCore('GetUndistortionParameters', oculusPtr, eye [, inputWidth][, inputHeight][, fov]);";
    synopsis[i++] = "[eyeRotStartMatrix, eyeRotEndMatrix] = PsychOculusVRCore('GetEyeTimewarpMatrices', oculusPtr, eye [, waitForTimewarpPoint=0]);";
    synopsis[i++] = "PsychOculusVRCore('EndFrameTiming', oculusPtr);";
    synopsis[i++] = "result = PsychOculusVRCore('LatencyTester', oculusPtr, cmd);";
    synopsis[i++] = NULL;  //this tells PsychOculusVRDisplaySynopsis where to stop

    if (i > MAX_SYNOPSIS_STRINGS) {
        PrintfExit("%s: Increase dimension of synopsis[] from %ld to at least %ld and recompile.", __FILE__, (long) MAX_SYNOPSIS_STRINGS, (long) i);
    }
}

PsychError PsychOculusVRDisplaySynopsis(void)
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

static void PsychOculusLogCB(int level, const char* message)
{
    if ((level == ovrLogLevel_Error && verbosity > 0) ||
        (level == ovrLogLevel_Info && verbosity > 2)  ||
        (level == ovrLogLevel_Debug && verbosity > 4)) {
            printf("PsychOculusVRCore-Runtime: %s\n", message);
    }
}

void PsychOculusVRCheckInit(psych_bool dontfail)
{
    ovrInitParams iparms;
    memset(&iparms, 0, sizeof(iparms));
    iparms.Flags = ovrInit_ForceNoDebug;
    iparms.LogCallback = PsychOculusLogCB;

    // Already initialized? No op then.
    if (initialized) return;

    // Initialize Oculus VR runtime with default parameters:
    if (ovr_Initialize(&iparms)) {
        if (verbosity >= 3) printf("PsychOculusVRCore-INFO: Oculus VR runtime version '%s' initialized.\n", ovr_GetVersionString());

        // Get count of available devices:
        available_devices = ovrHmd_Detect();
        if (available_devices < 0) {
            available_devices = 0;
            if (verbosity >= 2) printf("PsychOculusVRCore-WARNING: Could not connect to Oculus VR server process yet. Did you forget to start it?\n");
        }

        if (verbosity >= 3) printf("PsychOculusVRCore-INFO: At startup there are %i Oculus HMDs available.\n", available_devices);
        initialized = TRUE;
    }
    else {
        if (!dontfail)
            PsychErrorExitMsg(PsychError_system, "PsychOculusVRCore-ERROR: Initialization of VR runtime failed. Driver disabled!");
    }
}

void PsychOculusStop(int handle)
{
    PsychOculusDevice* oculus;
    oculus = PsychGetOculus(handle, TRUE);
    if (NULL == oculus || !oculus->isTracking) return;

    // Request stop of tracking:
    if (!ovrHmd_ConfigureTracking(oculus->hmd, 0, 0)) {
        if (verbosity >= 0) printf("PsychOculusVRCore-ERROR: Failed to stop tracking on device with handle %i [%s].\n", handle, ovrHmd_GetLastError(oculus->hmd));
        PsychErrorExitMsg(PsychError_system, "Stop of Oculus HMD tracking failed for reason given above.");
    }
    else if (verbosity >= 4) printf("PsychOculusVRCore-INFO: Tracking stopped on device with handle %i.\n", handle);

    oculus->isTracking = FALSE;

    return;
}

void PsychOculusClose(int handle)
{
    PsychOculusDevice* oculus;
    oculus = PsychGetOculus(handle, TRUE);
    if (NULL == oculus) return;

    // Stop device:
    PsychOculusStop(handle);

    // Release distortion meshes, if any:
    if (oculus->eyeDistortionMesh[0].pVertexData) {
        ovrHmd_DestroyDistortionMesh(&(oculus->eyeDistortionMesh[0]));
    }

    if (oculus->eyeDistortionMesh[1].pVertexData) {
        ovrHmd_DestroyDistortionMesh(&(oculus->eyeDistortionMesh[1]));
    }

    // Close the HMD:
    ovrHmd_Destroy(oculus->hmd);
    oculus->hmd = NULL;
    if (verbosity >= 4) printf("PsychOculusVRCore-INFO: Closed Oculus HMD with handle %i.\n", handle);

    // Done with this device:
    devicecount--;
}

void PsychOculusVRInit(void) {
    int handle;

    for (handle = 0 ; handle < MAX_PSYCH_OCULUS_DEVS; handle++)
        oculusdevices[handle].hmd = NULL;

    available_devices = 0;
    devicecount = 0;
    initialized = FALSE;
}

PsychError PsychOculusVRShutDown(void) {
    int handle;

    if (initialized) {
        for (handle = 0 ; handle < MAX_PSYCH_OCULUS_DEVS; handle++)
            PsychOculusClose(handle);

        // Shutdown runtime:
        ovr_Shutdown();

        if (verbosity >= 4) printf("PsychOculusVRCore-INFO: Oculus VR runtime shutdown complete.\n");
    }
    initialized = FALSE;

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRVerbosity(void)
{
    static char useString[] = "oldVerbosity = PsychOculusVRCore('Verbosity' [, verbosity]);";
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

PsychError PSYCHOCULUSVRGetCount(void)
{
    static char useString[] = "numHMDs = PsychOculusVRCore('GetCount');";
    static char synopsisString[] =  "Returns count of currently connected HMDs.\n"
                                    "Returns -1 if the runtime or server couldn't get initialized.\n";
    static char seeAlsoString[] = "Open";

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if( PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(0));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(TRUE);
    if (!initialized) {
        available_devices = -1;
    }
    else {
        available_devices = ovrHmd_Detect();
        if (available_devices < 0) {
            available_devices = -1;
            if (verbosity >= 2) printf("PsychOculusVRCore-WARNING: Could not connect to Oculus VR server process yet. Did you forget to start it?\n");
        }
    }

    PsychCopyOutDoubleArg(1, kPsychArgOptional, available_devices);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVROpen(void)
{
    static char useString[] = "[oculusPtr, modelName] = PsychOculusVRCore('Open' [, deviceIndex=0]);";
    //                          1          2                                        1
    static char synopsisString[] =
        "Open connection to Oculus VR HMD, return a 'oculusPtr' handle to it.\n\n"
        "The call tries to open the HMD with index 'deviceIndex', or the "
        "first detected HMD, if 'deviceIndex' is omitted. You can pass in a 'deviceIndex' "
        "of -1 to open an emulated HMD. It doesn't provide any sensor input, but allows "
        "some basic testing and debugging of VR software nonetheless.\n"
        "The returned handle can be passed to the other subfunctions to operate the device.\n"
        "A second return argument 'modelName' returns the model name string of the Oculus device.\n";

    static char seeAlsoString[] = "GetCount Close";

    PsychOculusDevice* oculus;
    int deviceIndex = 0;
    int handle = 0;
    unsigned int oldCaps;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(2));
    PsychErrorExit(PsychCapNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Find a free device slot:
    for (handle = 0; (handle < MAX_PSYCH_OCULUS_DEVS) && oculusdevices[handle].hmd; handle++);
    if ((handle >= MAX_PSYCH_OCULUS_DEVS) || oculusdevices[handle].hmd) PsychErrorExitMsg(PsychError_internal, "Maximum number of simultaneously open Oculus VR devices reached.");

    // Get optional Oculus device index:
    PsychCopyInIntegerArg(1, kPsychArgOptional, &deviceIndex);

    // Don't support anything than a single "default" OculusVR Rift yet - A limitation of the current SDK:
    if (deviceIndex < -1) PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' provided. Must be greater or equal to zero!");

    available_devices = ovrHmd_Detect();
    if (available_devices < 0) {
        available_devices = 0;
        if (verbosity >= 2) printf("PsychOculusVRCore-WARNING: Could not connect to Oculus VR server process yet. Did you forget to start it?\n");
    }

    if ((deviceIndex >= 0) && (deviceIndex >= available_devices)) {
        if (verbosity >= 0) printf("PsychOculusVRCore-ERROR: Invalid deviceIndex %i >= number of available HMDs %i.\n", deviceIndex, available_devices);
        PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' provided. Not enough HMDs available!");
    }

    // Zero init device structure:
    memset(&oculusdevices[handle], 0, sizeof(PsychOculusDevice));

    oculus = &oculusdevices[handle];

    // Try to open real or emulated HMD with deviceIndex:
    if (deviceIndex >= 0) {
        // The real thing:
        oculusdevices[handle].hmd = ovrHmd_Create(deviceIndex);
        if (NULL == oculusdevices[handle].hmd) {
            if (verbosity >= 0) {
                printf("PsychOculusVRCore-ERROR: Failed to connect to Oculus Rift with deviceIndex %i. This could mean that the device\n", deviceIndex);
                printf("PsychOculusVRCore-ERROR: is already in use by another application or driver.\n");
            }
            PsychErrorExitMsg(PsychError_user, "Could not connect to Rift device with given 'deviceIndex'! [ovrHmd_Create() failed]");
        }
        else if (verbosity >= 3) {
            printf("PsychOculusVRCore-INFO: Opened Oculus Rift with deviceIndex %i as handle %i.\n", deviceIndex, handle + 1);
        }
    }
    else {
        // Emulated: Simulate a Rift DK2.
        oculusdevices[handle].hmd = ovrHmd_CreateDebug(ovrHmd_DK2);
        if (verbosity >= 3) printf("PsychOculusVRCore-INFO: Opened an emulated Oculus Rift DK2 as handle %i.\n", handle + 1);
    }

    // Query current enabled caps:
    oldCaps = ovrHmd_GetEnabledCaps(oculus->hmd);

    // Stats for nerds:
    if (verbosity >= 3) {
        printf("PsychOculusVRCore-INFO: Product: %s - Manufacturer: %s - SerialNo: %s [VID: 0x%x PID: 0x%x]\n",
               oculus->hmd->ProductName, oculus->hmd->Manufacturer, oculus->hmd->SerialNumber, (int) oculus->hmd->VendorId, (int) oculus->hmd->ProductId);
        printf("PsychOculusVRCore-INFO: Firmware version: %i.%i\n", (int) oculus->hmd->FirmwareMajor, (int) oculus->hmd->FirmwareMinor);
        printf("PsychOculusVRCore-INFO: CameraFrustumHFovInRadians: %f - CameraFrustumVFovInRadians: %f\n", oculus->hmd->CameraFrustumHFovInRadians, oculus->hmd->CameraFrustumVFovInRadians);
        printf("PsychOculusVRCore-INFO: CameraFrustumNearZInMeters: %f - CameraFrustumFarZInMeters:  %f\n", oculus->hmd->CameraFrustumNearZInMeters, oculus->hmd->CameraFrustumFarZInMeters);
        printf("PsychOculusVRCore-INFO: Panel size in pixels w x h = %i x %i [WindowPos %i x %i]\n", oculus->hmd->Resolution.w, oculus->hmd->Resolution.h, oculus->hmd->WindowsPos.x, oculus->hmd->WindowsPos.y);
        printf("PsychOculusVRCore-INFO: DisplayDeviceName: %s\n", oculus->hmd->DisplayDeviceName);
        printf("PsychOculusVRCore-INFO: Caps: LowPersistence=%i : DynamicPrediction=%i\n", (oldCaps & ovrHmdCap_LowPersistence) ? 1 : 0, (oldCaps & ovrHmdCap_DynamicPrediction) ? 1 : 0);
        printf("PsychOculusVRCore-INFO: ----------------------------------------------------------------------------------\n");
    }

    // Increment count of open devices:
    devicecount++;

    // Return device handle: We use 1-based handle indexing to make life easier for Octave/Matlab:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, handle + 1);

    // Return product name:
    PsychCopyOutCharArg(2, kPsychArgOptional, (const char*) oculus->hmd->ProductName);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRClose(void)
{
    static char useString[] = "PsychOculusVRCore('Close' [, oculusPtr]);";
    //                                                      1
    static char synopsisString[] =
        "Close connection to Oculus Rift device 'oculusPtr'. Do nothing if no such device is open.\n"
        "If the optional 'oculusPtr' is omitted, then close all open devices and shutdown the driver, "
        "ie. perform the same cleanup as if 'clear PsychOculusVRCore' would be executed.\n";
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
        PsychOculusVRShutDown();
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRSetDynamicPrediction(void)
{
    static char useString[] = "oldDynamicPrediction = PsychOculusVRCore('SetDynamicPrediction', oculusPtr [, dynamicPrediction]);";
    //                                                                                          1            2
    static char synopsisString[] =
        "Enable or disable dynamic prediction mode on Oculus device 'oculusPtr'.\n"
        "'dynamicPrediction' 1 = Enable dynamic prediction, 0 = Disable dynamic prediction.\n"
        "Returns previous 'dynamic prediction ' setting.\n";
    static char seeAlsoString[] = "";

    int handle, dynamicPrediction;
    PsychOculusDevice *oculus;
    unsigned int oldCaps;

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

    // Query current enabled caps:
    oldCaps = ovrHmd_GetEnabledCaps(oculus->hmd);
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) (oldCaps & ovrHmdCap_DynamicPrediction) ? 1 : 0);

    // Set new enabled HMD caps:
    if (PsychCopyInIntegerArg(2, kPsychArgOptional, &dynamicPrediction)) {
        oldCaps &= ~ovrHmdCap_DynamicPrediction;
        ovrHmd_SetEnabledCaps(oculus->hmd, oldCaps | ((dynamicPrediction > 0) ? ovrHmdCap_DynamicPrediction : 0));
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRSetLowPersistence(void)
{
    static char useString[] = "oldPersistence = PsychOculusVRCore('SetLowPersistence', oculusPtr [, lowPersistence]);";
    //                                                                                 1            2
    static char synopsisString[] =
        "Enable or disable low persistence mode on display panel of Oculus device 'oculusPtr'.\n"
        "'lowPersistence' 1 = Enable low persistence, 0 = Disable low persistence.\n"
        "In low persistence mode, the pixels will only light up for a fraction of a video "
        "refresh duration, thereby reducing motion blur due to smooth pursuit eye movements "
        "or other shuttering effects. Brightness of the display will be reduced though, and "
        "flicker sensitive people may perceive more flicker.\n"
        "Returns previous 'oldPersistence' setting.\n";
    static char seeAlsoString[] = "";

    int handle, lowPersistence;
    PsychOculusDevice *oculus;
    unsigned int oldCaps;

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

    // Query current enabled caps:
    oldCaps = ovrHmd_GetEnabledCaps(oculus->hmd);
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) (oldCaps & ovrHmdCap_LowPersistence) ? 1 : 0);

    // Set new enabled HMD caps:
    if (PsychCopyInIntegerArg(2, kPsychArgOptional, &lowPersistence)) {
        oldCaps &= ~ovrHmdCap_LowPersistence;
        ovrHmd_SetEnabledCaps(oculus->hmd, oldCaps | ((lowPersistence > 0) ? ovrHmdCap_LowPersistence : 0));
    }

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRStart(void)
{
    static char useString[] = "PsychOculusVRCore('Start', oculusPtr);";
    //                                                1
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
        if (verbosity >= 0) printf("PsychOculusVRCore-ERROR: Tried to start tracking on device %i, but tracking is already started.\n", handle);
        PsychErrorExitMsg(PsychError_user, "Tried to start tracking on HMD, but tracking already active.");
    }

    // Request start of tracking for retrieval of head orientation and position, with drift correction, e.g., via magnetometer.
    // Do not fail if retrieval of any of this information isn't supported by the given hardware, ie., the required set of caps is empty == 0.
    // Rift DK1 only had orientation tracking, with magnetometer based drift correction. Rift DK2 also has vision based position tracking and
    // drift correction. The software emulated Rift has none of these and just returns a "static" head. This will start tracking:
    if (!ovrHmd_ConfigureTracking(oculus->hmd, ovrTrackingCap_Orientation | ovrTrackingCap_MagYawCorrection | ovrTrackingCap_Position, 0)) {
        if (verbosity >= 0) printf("PsychOculusVRCore-ERROR: Failed to start tracking on device with handle %i [%s].\n", handle, ovrHmd_GetLastError(oculus->hmd));
        PsychErrorExitMsg(PsychError_system, "Start of Oculus HMD tracking failed for reason given above.");
    }
    else if (verbosity >= 4) printf("PsychOculusVRCore-INFO: Tracking started on device with handle %i.\n", handle);

    oculus->frameIndex = 0;
    ovrHmd_ResetFrameTiming(oculus->hmd, oculus->frameIndex);

    oculus->isTracking = TRUE;

    // Tracking is running.
    return(PsychError_none);
}

PsychError PSYCHOCULUSVRStop(void)
{
    static char useString[] = "PsychOculusVRCore('Stop', oculusPtr);";
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

PsychError PSYCHOCULUSVRGetTrackingState(void)
{
    static char useString[] = "state = PsychOculusVRCore('GetTrackingState', oculusPtr [, predictionTime=0]);";
    //                         1                                             2            3
    static char synopsisString[] =
        "Return current state of head position and orientation tracking for Oculus device 'oculusPtr'.\n"
        "Head position and orientation is predicted for target time 'predictionTime' in seconds if provided, "
        "based on the latest measurements from the tracking hardware. If 'predictionTime' is omitted or set "
        "to zero, then no prediction is performed and the current state based on latest measurements is returned.\n\n"
        "'state' is a struct with fields reporting the following values:\n"
        "'Time' = Time in seconds of predicted tracking state.\n"
        "'Status' = Tracking status flags. +1 = Head orientation tracked, +2 = Head position tracked, +4 = Camera pose tracked "
        "+32 = Position tracking hardware connected, +128 = HMD display is connected and available.\n"
        "'HeadPose' = Head position [x, y, z] in meters and rotation as quaternion [rx, ry, rz, rw], all as a vector [x,y,z,rx,ry,rz,rw].\n"
        "'HeadLinearSpeed' = Linear velocity [vx,vy,vz] in meters/sec.\n"
        "'HeadAngularSpeed' = Angular velocity [rx,ry,rz] in radians/sec.\n"
        "'HeadLinearAcceleration' = Linear acceleration [ax,ay,az] in meters/sec^2.\n"
        "'HeadAngularAcceleration' = Angular acceleration [rax,ray,raz] in radians/sec^2.\n"
        "'CameraPose' as vector with position and orientation quaternion, like 'HeadPose'.\n"
        "'LeveledCameraPose' Like 'CameraPose' but aligned to the gravity vector of the world.\n"
        "'CameraFrustumHVFov' Horizontal and vertical field of view of the tracking camera in radians.\n"
        "'CameraFrustumNearFarZInMeters' Near and far limit of the camera view frustum in meters.\n"
        "'LastCameraFrameCounter' Last camera framecounter value of tracking camera.\n"
        "'RawSensorAcceleration' = Raw measured accelerometer reading in m/sec^2.\n"
        "'RawSensorGyroRate' = Raw gyrometer reading in rad/s.\n"
        "'RawMagnetometer' = Raw magnetic field in gauss.\n"
        "'SensorTemperature' = Sensor temperature in degrees Celsius.\n"
        "'IMUReadoutTime' = Readout time of the last IMU sample in seconds.\n"
        "\n";

    static char seeAlsoString[] = "Start Stop";

    PsychGenericScriptType *status;
    const char *FieldNames[] = {"Time", "Status", "HeadPose", "HeadLinearSpeed", "HeadAngularSpeed", "HeadLinearAcceleration", "HeadAngularAcceleration",
                                "CameraPose", "LeveledCameraPose", "LastCameraFrameCounter", "RawSensorAcceleration", "RawSensorGyroRate", "RawMagnetometer",
                                "SensorTemperature", "IMUReadoutTime", "CameraFrustumHVFov", "CameraFrustumNearFarZInMeters"};
    const int FieldCount = 17;
    PsychGenericScriptType *outMat;
    double *v;

    int handle;
    double predictionTime = 0.0;
    PsychOculusDevice *oculus;
    ovrTrackingState state;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    PsychCopyInDoubleArg(2, kPsychArgOptional, &predictionTime);

    // Get current tracking status info at time 0 ie., current measurements:
    state = ovrHmd_GetTrackingState(oculus->hmd, predictionTime);

    // Print out tracking status:
    if (verbosity >= 4) {
        printf("PsychOculusVRCore-INFO: Tracking state predicted for device %i at time %f.\n", handle, predictionTime);
        printf("PsychOculusVRCore-INFO: LastCameraFrameCounter = %i : Time %f : Status %i\n", state.LastCameraFrameCounter, state.HeadPose.TimeInSeconds, state.StatusFlags);
        printf("PsychOculusVRCore-INFO: HeadPose: Position    [x,y,z]   = [%f, %f, %f]\n", state.HeadPose.ThePose.Position.x, state.HeadPose.ThePose.Position.y, state.HeadPose.ThePose.Position.z);
        printf("PsychOculusVRCore-INFO: HeadPose: Orientation [x,y,z,w] = [%f, %f, %f, %f]\n", state.HeadPose.ThePose.Orientation.x, state.HeadPose.ThePose.Orientation.y, state.HeadPose.ThePose.Orientation.z, state.HeadPose.ThePose.Orientation.w);
    }

    PsychAllocOutStructArray(1, kPsychArgOptional, -1, FieldCount, FieldNames, &status);

    PsychSetStructArrayDoubleElement("Time", 0, state.HeadPose.TimeInSeconds, status);
    PsychSetStructArrayDoubleElement("Status", 0, state.StatusFlags, status);

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
    PsychSetStructArrayNativeElement("HeadPose", 0,	outMat, status);

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
    PsychSetStructArrayNativeElement("HeadAngularSpeed", 0,	outMat, status);

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

    // Camera pose:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 7, 1, &v, &outMat);
    v[0] = state.CameraPose.Position.x;
    v[1] = state.CameraPose.Position.y;
    v[2] = state.CameraPose.Position.z;

    v[3] = state.CameraPose.Orientation.x;
    v[4] = state.CameraPose.Orientation.y;
    v[5] = state.CameraPose.Orientation.z;
    v[6] = state.CameraPose.Orientation.w;
    PsychSetStructArrayNativeElement("CameraPose", 0, outMat, status);

    // Camera leveled pose:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 7, 1, &v, &outMat);
    v[0] = state.LeveledCameraPose.Position.x;
    v[1] = state.LeveledCameraPose.Position.y;
    v[2] = state.LeveledCameraPose.Position.z;

    v[3] = state.LeveledCameraPose.Orientation.x;
    v[4] = state.LeveledCameraPose.Orientation.y;
    v[5] = state.LeveledCameraPose.Orientation.z;
    v[6] = state.LeveledCameraPose.Orientation.w;
    PsychSetStructArrayNativeElement("LeveledCameraPose", 0, outMat, status);

    // LastCameraFrameCounter:
    PsychSetStructArrayDoubleElement("LastCameraFrameCounter", 0, (double) state.LastCameraFrameCounter, status);

    // RawSensorAcceleration:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
    v[0] = state.RawSensorData.Accelerometer.x;
    v[1] = state.RawSensorData.Accelerometer.y;
    v[2] = state.RawSensorData.Accelerometer.z;
    PsychSetStructArrayNativeElement("RawSensorAcceleration", 0, outMat, status);

    // RawSensorGyroRate:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
    v[0] = state.RawSensorData.Gyro.x;
    v[1] = state.RawSensorData.Gyro.y;
    v[2] = state.RawSensorData.Gyro.z;
    PsychSetStructArrayNativeElement("RawSensorGyroRate", 0, outMat, status);

    // RawMagnetometer:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 3, 1, &v, &outMat);
    v[0] = state.RawSensorData.Magnetometer.x;
    v[1] = state.RawSensorData.Magnetometer.y;
    v[2] = state.RawSensorData.Magnetometer.z;
    PsychSetStructArrayNativeElement("RawMagnetometer", 0, outMat, status);

    // SensorTemperature:
    PsychSetStructArrayDoubleElement("SensorTemperature", 0, (double) state.RawSensorData.Temperature, status);

    // IMU readout time:
    PsychSetStructArrayDoubleElement("IMUReadoutTime", 0, (double) state.RawSensorData.TimeInSeconds, status);

    // Camera frustum HFov and VFov in radians:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
    v[0] = oculus->hmd->CameraFrustumHFovInRadians;
    v[1] = oculus->hmd->CameraFrustumVFovInRadians;
    PsychSetStructArrayNativeElement("CameraFrustumHVFov", 0, outMat, status);

    // Camera frustum near and far clip plane in meters:
    v = NULL;
    PsychAllocateNativeDoubleMat(1, 2, 1, &v, &outMat);
    v[0] = oculus->hmd->CameraFrustumNearZInMeters;
    v[1] = oculus->hmd->CameraFrustumFarZInMeters;
    PsychSetStructArrayNativeElement("CameraFrustumNearFarZInMeters", 0, outMat, status);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRGetFovTextureSize(void)
{
    static char useString[] = "[width, height, fovPort] = PsychOculusVRCore('GetFovTextureSize', oculusPtr, eye [, fov=[HMDRecommended]][, pixelsPerDisplay=1]);";
    //                          1      2       3                                                 1          2      3                       4
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
        oculus->ofov[eyeIndex] = oculus->hmd->DefaultEyeFov[eyeIndex];
    }

    // Get optional pixelsPerDisplay parameter:
    pixelsPerDisplay = 1.0;
    PsychCopyInDoubleArg(4, kPsychArgOptional, &pixelsPerDisplay);
    if (pixelsPerDisplay <= 0.0) PsychErrorExitMsg(PsychError_user, "Invalid 'pixelsPerDisplay' specified. Must be greater than zero.");

    // Ask the api for optimal texture size, aka the size of the client draw buffer:
    oculus->texSize[eyeIndex] = ovrHmd_GetFovTextureSize(oculus->hmd, (ovrEyeType) eyeIndex, oculus->ofov[eyeIndex], (float) pixelsPerDisplay);

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

PsychError PSYCHOCULUSVRGetUndistortionParameters(void)
{
    static char useString[] = "[width, height, viewPx, viewPy, viewPw, viewPh, pptax, pptay, hmdShiftx, hmdShifty, hmdShiftz, meshVertices, meshIndices, uvScaleX, uvScaleY, uvOffsetX, uvOffsetY] = PsychOculusVRCore('GetUndistortionParameters', oculusPtr, eye [, inputWidth][, inputHeight][, fov]);";
    //                          1      2       3       4       5       6       7      8      9          10         11         12            13           14        15        16         17                                                              1          2    3           4              5
    static char synopsisString[] =
    "Return parameters needed for rendering and undistortion for Oculus device 'oculusPtr'.\n"
    "'eye' which eye to provide the data: 0 = Left, 1 = Right.\n"
    "'inputWidth' = Width of the rendered input image buffer in pixels.\n"
    "'inputHeight' = Height of the rendered input image buffer in pixels.\n"
    "'fov' Optional field of view in degrees, from line of sight: [leftdeg, rightdeg, updeg, downdeg]. "
    "You can pass in the 'fovPort' value returned from PsychOculusVR('GetFovTextureSize'); "
    "Defaults to whatever has been set for the given eye in the last call to PsychOculusVR('GetFovTextureSize'); if omitted.\n"
    "\n"
    "Return values:\n"
    "[width, height] = Width and height of client renderbuffers in pixels. Same as the provided 'inputWidth' and 'inputHeight'.\n"
    "[viewPx, viewPy, viewPw, viewPh] Render viewport [x,y,w,h] start (x,y) position, width and height. Mostly useless atm., as "
    "viewPx and viewPy are always zero, viewPw and viewPh are identical to width and height in the current driver design.\n"
    "[pptax pptay] = Pixels per tangens angle at display center in x direction and y direction.\n"
    "[hmdShiftx, hmdShifty, hmdShiftz] = HmdToEyeViewOffset 3D translation vector. Defines the location of the optical center of the eye "
    "relative to the origin of the local head reference frame, ie. the tracked head position.\n"
    "[meshVertices meshIndices] = Matrices with the definitions of the geometric undistortion triangle meshes. 'meshIndices' is "
    "a vector of 0-based indices into the vertex mesh 'meshVertices' for building a mesh out of triangles. 'meshVertices' is a "
    "matrix with as many columns as mesh vertices, each column encoding one vertex. Each row encodes for one vertex property. See "
    "the Oculus developers manual for explanation of the format.\n"
    "[uvScaleX, uvScaleY, uvOffsetX, uvOffsetY] Scaling factors and offsets for texture coordinates in normalized device coordinates.\n";
    static char seeAlsoString[] = "GetFovTextureSize";

    int handle, eyeIndex;
    PsychOculusDevice *oculus;
    int n, m, p, i;
    double *fov;
    unsigned int distortionCaps;
    ovrDistortionVertex* pVertexData;
    unsigned short* pIndexData;
    double *outVertexMesh, *outIndexMesh;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(17));
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
    oculus->eyeRenderDesc[eyeIndex] = ovrHmd_GetRenderDesc(oculus->hmd, (ovrEyeType) eyeIndex, oculus->ofov[eyeIndex]);

    // Override viewport size with the provided input texture size, as ovrHmd_GetRenderDesc() makes the
    // assumption that we render to a half-width viewport of the input renderbuffer, when we actually
    // use 2 separate input renderbuffers of appropriate size:
    oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Pos.x = 0;
    oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Pos.y = 0;
    oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Size.w = oculus->texSize[eyeIndex].w;
    oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Size.h = oculus->texSize[eyeIndex].h;

    if (verbosity > 3) {
        printf("PsychOculusVRCore-INFO: For HMD %i, eye %i - RenderDescription:\n", handle, eyeIndex);
        printf("PsychOculusVRCore-INFO: FoV: %f %f %f %f - %f %f %f %f\n", oculus->ofov[eyeIndex].LeftTan, oculus->ofov[eyeIndex].RightTan, oculus->ofov[eyeIndex].UpTan, oculus->ofov[eyeIndex].DownTan, oculus->eyeRenderDesc[eyeIndex].Fov.LeftTan, oculus->eyeRenderDesc[eyeIndex].Fov.RightTan, oculus->eyeRenderDesc[eyeIndex].Fov.UpTan, oculus->eyeRenderDesc[eyeIndex].Fov.DownTan);
        printf("PsychOculusVRCore-INFO: DistortedViewport: [x,y,w,h] = [%i, %i, %i, %i]\n", oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Pos.x, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Pos.y, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Size.w, oculus->eyeRenderDesc[eyeIndex].DistortedViewport.Size.h);
        printf("PsychOculusVRCore-INFO: PixelsPerTanAngleAtCenter: %f x %f\n", oculus->eyeRenderDesc[eyeIndex].PixelsPerTanAngleAtCenter.x, oculus->eyeRenderDesc[eyeIndex].PixelsPerTanAngleAtCenter.y);
        printf("PsychOculusVRCore-INFO: HmdToEyeViewOffset: [x,y,z] = [%f, %f, %f]\n", oculus->eyeRenderDesc[eyeIndex].HmdToEyeViewOffset.x, oculus->eyeRenderDesc[eyeIndex].HmdToEyeViewOffset.y, oculus->eyeRenderDesc[eyeIndex].HmdToEyeViewOffset.z);
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
    PsychCopyOutDoubleArg(9,  kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].HmdToEyeViewOffset.x);
    PsychCopyOutDoubleArg(10, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].HmdToEyeViewOffset.y);
    PsychCopyOutDoubleArg(11, kPsychArgOptional, oculus->eyeRenderDesc[eyeIndex].HmdToEyeViewOffset.z);

    // See enum ovrDistortionCaps: ovrDistortionCap_TimeWarp | ovrDistortionCap_Vignette | ovrDistortionCap_Overdrive
    // distortionCaps = ovrDistortionCap_LinuxDevFullscreen | ovrDistortionCap_HqDistortion | ovrDistortionCap_FlipInput;
    // TODO: distortionCaps are not used by SDK version 0.5 for ovrHmd_CreateDistortionMesh(), but maybe by later
    // versions? Need to recheck once we use higher SDK versions.
    distortionCaps = 0;
    if (!ovrHmd_CreateDistortionMesh(oculus->hmd, oculus->eyeRenderDesc[eyeIndex].Eye, oculus->eyeRenderDesc[eyeIndex].Fov, distortionCaps, &(oculus->eyeDistortionMesh[eyeIndex]))) {
        if (verbosity > 0) printf("PsychOculusVRCore-ERROR: Failed to compute distortion mesh for HMD %i, eye %i: [%s]\n", handle, eyeIndex, ovrHmd_GetLastError(oculus->hmd));
        PsychErrorExitMsg(PsychError_system, "Failed to compute distortion mesh for eye.");
    }

    if (verbosity > 3) {
        printf("PsychOculusVRCore-INFO: Distortion mesh has %i vertices, %i indices for triangles.\n",
               oculus->eyeDistortionMesh[eyeIndex].VertexCount,
               oculus->eyeDistortionMesh[eyeIndex].IndexCount);
    }

    // Return vertex data for the distortion mesh:

    // Each mesh has 10 parameters per vertex:
    m = 10;
    // For given number of vertices:
    n = oculus->eyeDistortionMesh[eyeIndex].VertexCount;
    // And one layer for a 2D matrix:
    p = 1;
    PsychAllocOutDoubleMatArg(12, kPsychArgOptional, m, n, p, &outVertexMesh);

    pVertexData = oculus->eyeDistortionMesh[eyeIndex].pVertexData;
    for (i = 0; i < n; i++) {
        // Store i'th column for i'th vertex:

        // output vertex 2D (x,y) position:
        *(outVertexMesh++) = pVertexData->ScreenPosNDC.x;
        *(outVertexMesh++) = pVertexData->ScreenPosNDC.y;

        // Timewarp lerp factor:
        *(outVertexMesh++) = pVertexData->TimeWarpFactor;

        // Vignette fade factor:
        *(outVertexMesh++) = pVertexData->VignetteFactor;

        // The tangents of the horizontal and vertical eye angles for the red channel.
        *(outVertexMesh++) = pVertexData->TanEyeAnglesR.x;
        *(outVertexMesh++) = pVertexData->TanEyeAnglesR.y * -1.0;

        // The tangents of the horizontal and vertical eye angles for the green channel.
        *(outVertexMesh++) = pVertexData->TanEyeAnglesG.x;
        *(outVertexMesh++) = pVertexData->TanEyeAnglesG.y * -1.0;

        // The tangents of the horizontal and vertical eye angles for the blue channel.
        *(outVertexMesh++) = pVertexData->TanEyeAnglesB.x;
        *(outVertexMesh++) = pVertexData->TanEyeAnglesB.y * -1.0;

        // Advance to next vertex in mesh:
        pVertexData++;
    }

    // Return index data for the distortion mesh: The mesh is composed of triangles.
    m = p = 1;
    n = oculus->eyeDistortionMesh[eyeIndex].IndexCount;
    PsychAllocOutDoubleMatArg(13, kPsychArgOptional, m, n, p, &outIndexMesh);

    pIndexData = oculus->eyeDistortionMesh[eyeIndex].pIndexData;
    for (i = 0; i < n; i++) {
        // Store i'th index - Convert uint16 to double:
        *(outIndexMesh++) = (double) *(pIndexData++);
    }

    // Get UV texture sampling scale and offset:
    ovrHmd_GetRenderScaleAndOffset(oculus->eyeRenderDesc[eyeIndex].Fov,
                                   oculus->texSize[eyeIndex], oculus->eyeRenderDesc[eyeIndex].DistortedViewport,
                                   (ovrVector2f*) &(oculus->UVScaleOffset[eyeIndex]));

    // EyeToSourceUVScale:
    PsychCopyOutDoubleArg(14, kPsychArgOptional, oculus->UVScaleOffset[eyeIndex][0].x);
    PsychCopyOutDoubleArg(15, kPsychArgOptional, oculus->UVScaleOffset[eyeIndex][0].y);

    // EyeToSourceUVOffset:
    PsychCopyOutDoubleArg(16, kPsychArgOptional, oculus->UVScaleOffset[eyeIndex][1].x);
    PsychCopyOutDoubleArg(17, kPsychArgOptional, oculus->UVScaleOffset[eyeIndex][1].y);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRGetEyeTimewarpMatrices(void)
{
    static char useString[] = "[eyeRotStartMatrix, eyeRotEndMatrix] = PsychOculusVRCore('GetEyeTimewarpMatrices', oculusPtr, eye [, waitForTimewarpPoint=0]);";
    //                          1                  2                                                              1          2      3
    static char synopsisString[] =
    "Return eye warp rotation matrices for timewarped undistortion for Oculus device 'oculusPtr'.\n"
    "'eye' which eye to provide the data: 0 = Left, 1 = Right.\n"
    "'waitForTimewarpPoint' If set to 1, stall execution of calling thread until next time warp point is reached. Defaults to zero.\n"
    "Return values are 4x4 'eyeRotStartMatrix' and 'eyeRotEndMatrix' for given eye.\n";
    static char seeAlsoString[] = "";

    int handle, eyeIndex, waitForTimewarpPoint;
    PsychOculusDevice *oculus;
    ovrEyeType eye;
    double tNow;
    int i, j;
    double *startMatrix, *endMatrix;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(2));
    PsychErrorExit(PsychCapNumInputArgs(3));
    PsychErrorExit(PsychRequireNumInputArgs(2));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Get eye index - left = 0, right = 1:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &eyeIndex);
    if (eyeIndex < 0 || eyeIndex > 1) PsychErrorExitMsg(PsychError_user, "Invalid 'eye' specified. Must be 0 or 1 for left- or right eye.");

    // Get waitForTimewarpPoint:
    waitForTimewarpPoint = 0;
    PsychCopyInIntegerArg(3, kPsychArgOptional, &waitForTimewarpPoint);
    if (waitForTimewarpPoint < 0 || waitForTimewarpPoint > 1) PsychErrorExitMsg(PsychError_user, "Invalid 'waitForTimewarpPoint' specified. Must be 0 or 1.");

    if (waitForTimewarpPoint) {
        // Wait till time-warp point to reduce latency.
        tNow = ovr_GetTimeInSeconds();
        if (verbosity > 3)
            printf("PsychOculusVRCore-INFO: Waiting for %f msecs until next TimewarpPointSeconds %f.\n", 1000.0 * (oculus->frameTiming.TimewarpPointSeconds - tNow), oculus->frameTiming.TimewarpPointSeconds);

        ovr_WaitTillTime(oculus->frameTiming.TimewarpPointSeconds);
    }

    eye = oculus->hmd->EyeRenderOrder[eyeIndex];
    oculus->headPose[eye] = ovrHmd_GetHmdPosePerEye(oculus->hmd, eye);

    ovrHmd_GetEyeTimewarpMatrices(oculus->hmd, (ovrEyeType) eyeIndex,
                                  oculus->headPose[eyeIndex],
                                  oculus->timeWarpMatrices);

    PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 4, 4, 1, &startMatrix);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            *(startMatrix++) = (double) oculus->timeWarpMatrices[0].M[j][i];

    PsychAllocOutDoubleMatArg(2, kPsychArgOptional, 4, 4, 1, &endMatrix);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            *(endMatrix++) = (double) oculus->timeWarpMatrices[1].M[j][i];

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRGetStaticRenderParameters(void)
{
    static char useString[] = "[projL, projR] = PsychOculusVRCore('GetStaticRenderParameters', oculusPtr [, clipNear=0.01][, clipFar=10000.0]);";
    //                          1      2                                                       1            2                3
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
    M = ovrMatrix4f_Projection(oculus->eyeRenderDesc[0].Fov, (float) clip_near, (float) clip_far, ovrProjection_RightHanded | ovrProjection_ClipRangeOpenGL);
    PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 4, 4, 1, &outM);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            *(outM++) = (double) M.M[j][i];

    // Return right projection matrix as return argument 2:
    M = ovrMatrix4f_Projection(oculus->eyeRenderDesc[1].Fov, (float) clip_near, (float) clip_far, ovrProjection_RightHanded | ovrProjection_ClipRangeOpenGL);
    PsychAllocOutDoubleMatArg(2, kPsychArgOptional, 4, 4, 1, &outM);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            *(outM++) = (double) M.M[j][i];

    return(PsychError_none);
}


PsychError PSYCHOCULUSVRStartRender(void)
{
    static char useString[] = "[eyePoseL, eyePoseR, tracked, frameTiming] = PsychOculusVRCore('StartRender', oculusPtr);";
    //                          1         2         3        4                                               1
    static char synopsisString[] =
    "Mark start of a new 3D head tracked render cycle for Oculus device 'oculusPtr'.\n"
    "Return values are the vectors which define the two eye cameras positions and orientations "
    "for the left eye and right eye 'eyePoseL' and 'eyePoseR'. The vectors are of form "
    "[tx, ty, tz, rx, ry, rz, rw] - A 3 component 3D position followed by a 4 component rotation "
    "quaternion.\n"
    "'tracked' Tracking status flags: 0 = Head not tracked at the moment. 1 = Head orientation tracked. "
    "2 = Head position tracked (DK2 and later). 3 = Head position and orientation tracked (DK2 and later). "
    "4 = Camera pose tracked, 7 = 1+2+4 = Camera pose and Head position and orientation tracked.\n"
    "\n"
    "'frameTiming' Vector with predicted timing information for this frame with following elements:\n"
    "[1] = DeltaSeconds since last frame.\n"
    "[2] = ThisFrameSeconds start of scanout of this frame.\n"
    "[3] = TimewarpPointSeconds\n"
    "[4] = NextFrameSeconds\n"
    "[5] = ScanoutMidpointSeconds\n"
    "[6] = EyeScanoutSeconds[0]\n"
    "[7] = EyeScanoutSeconds[1]\n"
    "\n";
    static char seeAlsoString[] = "GetEyePose EndFrameTiming";

    int handle;
    PsychOculusDevice *oculus;
    ovrTrackingState os;
    ovrVector3f hmdToEyeViewOffset[2];
    #ifdef  __cplusplus
    int i, j;
    #endif
    double *outM;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(6));
    PsychErrorExit(PsychCapNumInputArgs(1));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Mark beginning of frame rendering. This takes timstamps and stuff:
    oculus->frameTiming = ovrHmd_BeginFrameTiming(oculus->hmd, oculus->frameIndex);

    // Get current eye poses for both eyes:
    hmdToEyeViewOffset[0] = oculus->eyeRenderDesc[0].HmdToEyeViewOffset;
    hmdToEyeViewOffset[1] = oculus->eyeRenderDesc[1].HmdToEyeViewOffset;
    ovrHmd_GetEyePoses(oculus->hmd, oculus->frameIndex, hmdToEyeViewOffset,
                       oculus->outEyePoses, &os);

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

    // Report tracking state of head orientation, head position and camera pose:
    PsychCopyOutDoubleArg(3, kPsychArgOptional, (double)
                          (os.StatusFlags & (ovrStatus_OrientationTracked | ovrStatus_PositionTracked | ovrStatus_CameraPoseTracked)));

    // Report frame timing info:
    PsychAllocOutDoubleMatArg(4, kPsychArgOptional, 1, 7, 1, &outM);
    outM[0] = oculus->frameTiming.DeltaSeconds;
    outM[1] = oculus->frameTiming.ThisFrameSeconds;
    outM[2] = oculus->frameTiming.TimewarpPointSeconds;
    outM[3] = oculus->frameTiming.NextFrameSeconds;
    outM[4] = oculus->frameTiming.ScanoutMidpointSeconds;
    outM[5] = oculus->frameTiming.EyeScanoutSeconds[0];
    outM[6] = oculus->frameTiming.EyeScanoutSeconds[1];
    
    // Only compile this if C++ compiled:
    #ifdef  __cplusplus
    {
        // Invert camera translation vectors:
        oculus->outEyePoses[0].Position.x *= -1;
        oculus->outEyePoses[0].Position.y *= -1;
        oculus->outEyePoses[0].Position.z *= -1;

        oculus->outEyePoses[1].Position.x *= -1;
        oculus->outEyePoses[1].Position.y *= -1;
        oculus->outEyePoses[1].Position.z *= -1;

        // Return left modelview matrix as return argument 3:
        OVR::Quatf orientationL = OVR::Quatf(oculus->outEyePoses[0].Orientation);
        OVR::Matrix4f viewL = OVR::Matrix4f(orientationL.Inverted()) * OVR::Matrix4f::Translation(oculus->outEyePoses[0].Position);

        PsychAllocOutDoubleMatArg(5, kPsychArgOptional, 4, 4, 1, &outM);
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                *(outM++) = (double) viewL.M[j][i];

        // Return right modelview matrix as return argument 4:
        OVR::Quatf orientationR = OVR::Quatf(oculus->outEyePoses[1].Orientation);
        OVR::Matrix4f viewR = OVR::Matrix4f(orientationR.Inverted()) * OVR::Matrix4f::Translation(oculus->outEyePoses[1].Position);

        PsychAllocOutDoubleMatArg(6, kPsychArgOptional, 4, 4, 1, &outM);
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                *(outM++) = (double) viewR.M[j][i];
    }
    #endif

    return(PsychError_none);
}

PsychError PSYCHOCULUSVREndFrameTiming(void)
{
    static char useString[] = "PsychOculusVRCore('EndFrameTiming', oculusPtr);";
    //                                                             1
    static char synopsisString[] =
    "Mark end of a new 3D head tracked render and presentation cycle for Oculus device 'oculusPtr'.\n"
    "You usually won't call this function yourself, but Screen('Flip') will call it automatically for you "
    "at the appropriate moment.\n";
    static char seeAlsoString[] = "StartRender";

    int handle;
    PsychOculusDevice *oculus;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(0));
    PsychErrorExit(PsychCapNumInputArgs(1));
    PsychErrorExit(PsychRequireNumInputArgs(1));

    // Make sure driver is initialized:
    PsychOculusVRCheckInit(FALSE);

    // Get device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);
    oculus = PsychGetOculus(handle, FALSE);

    // Tell runtime that the last rendered frame was just presented onscreen,
    // aka OpenGL bufferswap completed:
    // Note: This also runs processing for the DK2 latency tester.
    ovrHmd_EndFrameTiming(oculus->hmd);
    oculus->frameIndex++;

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRGetEyePose(void)
{
    static char useString[] = "[eyePose, eyeIndex] = PsychOculusVRCore('GetEyePose', oculusPtr, renderPass);";
    //                          1        2                                           1          2
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
    eye = oculus->hmd->EyeRenderOrder[renderPass];
    oculus->headPose[eye] = ovrHmd_GetHmdPosePerEye(oculus->hmd, eye);

    // Eye pose as raw data:
    PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 1, 7, 1, &outM);

    // Position (x,y,z):
    outM[0] = oculus->headPose[eye].Position.x;
    outM[1] = oculus->headPose[eye].Position.y;
    outM[2] = oculus->headPose[eye].Position.z;

    // Orientation as a quaternion (x,y,z,w):
    outM[3] = oculus->headPose[eye].Orientation.x;
    outM[4] = oculus->headPose[eye].Orientation.y;
    outM[5] = oculus->headPose[eye].Orientation.z;
    outM[6] = oculus->headPose[eye].Orientation.w;

    // Copy out preferred eye render order for info:
    PsychCopyOutDoubleArg(2, kPsychArgOptional, (double) eye);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRGetHSWState(void)
{
    static char useString[] = "showHSW = PsychOculusVRCore('GetHSWState', oculusPtr [, dismiss]);";
    //                         1                                          1            2
    static char synopsisString[] =
    "Return if health and safety warnings popup should be displayed for Oculus device 'oculusPtr'.\n"
    "'showHSW' is 1 if the popup should be displayed at the moment, 0 otherwise.\n"
    "You can set the optional parameter 'dismiss' to 1 to tell the runtime that the user "
    "dismissed the health and safety warning popup.\n"
    "Display of these warnings is mandated by the Oculus VR terms of the license.\n";

    static char seeAlsoString[] = "";

    int handle, dismiss;
    PsychOculusDevice *oculus;
    ovrHSWDisplayState hasWarningState;

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

    // Get optional 'dismiss' flag:
    if (PsychCopyInIntegerArg(2, kPsychArgOptional, &dismiss) && dismiss) {
        // Tell runtime the user accepted and dismissed the HSD display:
        ovrHmd_DismissHSWDisplay(oculus->hmd);
    }

    ovrHmd_GetHSWDisplayState(oculus->hmd, &hasWarningState);
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) (hasWarningState.Displayed) ? 1 : 0);

    return(PsychError_none);
}

PsychError PSYCHOCULUSVRLatencyTester(void)
{
    static char useString[] = "result = PsychOculusVRCore('LatencyTester', oculusPtr, cmd);";
    //                          1                                          1          2
    static char synopsisString[] =
    "Drive or query latency tester for Oculus Rift DK2 device 'oculusPtr'.\n"
    "'cmd' is the command code:\n"
    "0 = Get next test draw color from 0 - 255 as [r,g,b] vector.\n"
    "1 = Get most recent latency test results in seconds as a 5 component vector 'latency':\n"
    "latency[1] = LatencyRender:     (seconds) Last time between render IMU sample and scanout.\n"
    "latency[2] = LatencyTimewarp:   (seconds) Last time between timewarp IMU sample and scanout.\n"
    "latency[3] = LatencyPostPresent (seconds) Average time between Vsync and scanout.\n"
    "latency[4] = ErrorRender        (seconds) Last error in render predicted scanout time.\n"
    "latency[5] = ErrorTimewarp      (seconds) Last error in timewarp predicted scanout time.\n"
    "\n"
    "If 'cmd' can't get executed, then an empty result [] is returned.\n"
    "\n";

    static char seeAlsoString[] = "";

    int handle, cmd;
    float latencies[5];
    PsychOculusDevice *oculus;
    double *outM;

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

    // Get cmd code:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &cmd);

    if (cmd == 0) {
        // Get next latency test marker pixel color, update the runtimes internal
        // target match color with it, return color to us for our client rendering:
        // Return required clear color for latency test quad if latencytest is enabled, empty otherwise:
        if (ovrHmd_GetLatencyTest2DrawColor(oculus->hmd, &(oculus->rgbColorOut[0]))) {
            PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 1, 3, 1, &outM);
            outM[0] = (double) oculus->rgbColorOut[0];
            outM[1] = (double) oculus->rgbColorOut[1];
            outM[2] = (double) oculus->rgbColorOut[2];
            oculus->latencyTestActive = TRUE;
        }
        else {
            PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 0, 0, 1, &outM);
            oculus->latencyTestActive = FALSE;
        }
    }

    if (cmd == 1) {
        // Shall we do processing for the latency tester and return results?
        if (oculus->latencyTestActive) {
            // Seems so:
            if (ovrHmd_GetFloatArray(oculus->hmd, "DK2Latency", latencies, 5) == 5) {
                if (verbosity > 3) {
                    printf("PsychOculusVRCore-EndFrameTiming: Latency test result is: %f : %f : %f : %f : %f\n",
                           latencies[0] * 1000.0, latencies[1] * 1000.0, latencies[2] * 1000.0, latencies[3] * 1000.0,
                           latencies[4] * 1000.0);
                }

                // Return results:
                PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 1, 5, 1, &outM);
                outM[0] = (double) latencies[0];
                outM[1] = (double) latencies[1];
                outM[2] = (double) latencies[2];
                outM[3] = (double) latencies[3];
                outM[4] = (double) latencies[4];
            }
        }
        else {
            PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 0, 0, 1, &outM);
        }
    }

    return(PsychError_none);
}
