/*
 * PsychToolbox/Source/Common/PsychOculusVRCore1/RegisterProject.c
 *
 * PROJECTS: PsychOculusVRCore1 only.
 *
 * AUTHORS:
 *
 * mario.kleiner.de@gmail.com   mk
 *
 * PLATFORMS:   All.
 *
 * DESCRIPTION:
 *
 * A Psychtoolbox driver for the Oculus VR virtual reality
 * head sets, using the OculusVR 1.0 SDK and runtime and later.
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

//begin include once
#include "Psych.h"
#include "PsychOculusVR1.h"

PsychError PsychModuleInit(void)
{
    // Register the project exit function
    PsychErrorExit(PsychRegisterExit(&PsychOculusVR1ShutDown));

    // Register the project function which is called when the module
    // is invoked with no named subfunction:
    PsychErrorExit(PsychRegister(NULL,  &PsychOculusVR1DisplaySynopsis));

    PsychErrorExit(PsychRegister("Open", &PSYCHOCULUSVR1Open));
    PsychErrorExit(PsychRegister("Close", &PSYCHOCULUSVR1Close));
    PsychErrorExit(PsychRegister("Start", &PSYCHOCULUSVR1Start));
    PsychErrorExit(PsychRegister("Stop", &PSYCHOCULUSVR1Stop));
    PsychErrorExit(PsychRegister("GetCount", &PSYCHOCULUSVR1GetCount));
    PsychErrorExit(PsychRegister("GetTrackingState", &PSYCHOCULUSVR1GetTrackingState));
    PsychErrorExit(PsychRegister("GetInputState", &PSYCHOCULUSVR1GetInputState));
    PsychErrorExit(PsychRegister("GetFovTextureSize", &PSYCHOCULUSVR1GetFovTextureSize));
    PsychErrorExit(PsychRegister("GetUndistortionParameters", &PSYCHOCULUSVR1GetUndistortionParameters));
    PsychErrorExit(PsychRegister("GetTrackersState", &PSYCHOCULUSVR1GetTrackersState));
    PsychErrorExit(PsychRegister("StartRender", &PSYCHOCULUSVR1StartRender));
    PsychErrorExit(PsychRegister("GetStaticRenderParameters", &PSYCHOCULUSVR1GetStaticRenderParameters));
    PsychErrorExit(PsychRegister("EndFrameRender", &PSYCHOCULUSVR1EndFrameRender));
    PsychErrorExit(PsychRegister("PresentFrame", &PSYCHOCULUSVR1PresentFrame));
    PsychErrorExit(PsychRegister("Verbosity", &PSYCHOCULUSVR1Verbosity));
    PsychErrorExit(PsychRegister("GetEyePose", &PSYCHOCULUSVR1GetEyePose));
    PsychErrorExit(PsychRegister("SetHUDState", &PSYCHOCULUSVR1SetHUDState));
    PsychErrorExit(PsychRegister("VRAreaBoundary", &PSYCHOCULUSVR1VRAreaBoundary));
    PsychErrorExit(PsychRegister("TestVRBoundary", &PSYCHOCULUSVR1TestVRBoundary));
    PsychErrorExit(PsychRegister("TestVRBoundaryPoint", &PSYCHOCULUSVR1TestVRBoundaryPoint));
    PsychErrorExit(PsychRegister("RecenterTrackingOrigin", &PSYCHOCULUSVR1RecenterTrackingOrigin));
    PsychErrorExit(PsychRegister("TrackingOriginType", &PSYCHOCULUSVR1TrackingOriginType));
    PsychErrorExit(PsychRegister("GetPerformanceStats", &PSYCHOCULUSVR1GetPerformanceStats));
    PsychErrorExit(PsychRegister("HapticPulse", &PSYCHOCULUSVR1HapticPulse));
    PsychErrorExit(PsychRegister("FloatsProperty", &PSYCHOCULUSVR1FloatsProperty));
    PsychErrorExit(PsychRegister("FloatProperty", &PSYCHOCULUSVR1FloatProperty));
    PsychErrorExit(PsychRegister("StringProperty", &PSYCHOCULUSVR1StringProperty));
    //PsychErrorExit(PsychRegister("", &PSYCHOCULUSVR1));
    //PsychErrorExit(PsychRegister("", &PSYCHOCULUSVR1));
    //PsychErrorExit(PsychRegister("", &PSYCHOCULUSVR1));
    //PsychErrorExit(PsychRegister("", &PSYCHOCULUSVR1));
    //PsychErrorExit(PsychRegister("", &PSYCHOCULUSVR1));
    //PsychErrorExit(PsychRegister("", &PSYCHOCULUSVR1));

    PsychErrorExit(PsychRegister("CreateRenderTextureChain", &PSYCHOCULUSVR1CreateRenderTextureChain));
    PsychErrorExit(PsychRegister("GetNextTextureHandle", &PSYCHOCULUSVR1GetNextTextureHandle));
    PsychErrorExit(PsychRegister("CreateMirrorTexture", &PSYCHOCULUSVR1CreateMirrorTexture));

    // Report the version
    PsychErrorExit(PsychRegister("Version", &MODULEVersion));

    // Register the module name & authors:
    PsychErrorExit(PsychRegister("PsychOculusVRCore1", NULL));
    PsychSetModuleAuthorByInitials("mk");

    // Register synopsis and named subfunctions.
    InitializeSynopsis();

    // Preinit everything:
    PsychOculusVR1Init();

    return(PsychError_none);
}
