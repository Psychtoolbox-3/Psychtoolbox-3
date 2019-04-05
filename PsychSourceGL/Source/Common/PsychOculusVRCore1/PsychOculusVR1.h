/*
 * PsychToolbox/Source/Common/PsychOculusVRCore1/PsychOculusVR1.h
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
 * head sets, using the OculusVR 1.0 SDK and runtime.
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

// begin include once
#ifndef PSYCH_IS_INCLUDED_PsychOculusVR1
#define PSYCH_IS_INCLUDED_PsychOculusVR1

//project includes
#include "Psych.h"

#ifdef  __cplusplus
extern "C" {
#endif

//function prototypes
PsychError MODULEVersion(void);
void InitializeSynopsis(void);
PsychError PsychOculusVR1DisplaySynopsis(void);

void PsychOculusVR1Init(void);
PsychError PsychOculusVR1ShutDown(void);

PsychError PSYCHOCULUSVR1Verbosity(void);
PsychError PSYCHOCULUSVR1Open(void);
PsychError PSYCHOCULUSVR1Close(void);
PsychError PSYCHOCULUSVR1SetHUDState(void);
PsychError PSYCHOCULUSVR1Start(void);
PsychError PSYCHOCULUSVR1Stop(void);
PsychError PSYCHOCULUSVR1GetCount(void);
PsychError PSYCHOCULUSVR1GetTrackingState(void);
PsychError PSYCHOCULUSVR1GetInputState(void);
PsychError PSYCHOCULUSVR1GetFovTextureSize(void);
PsychError PSYCHOCULUSVR1CreateRenderTextureChain(void);
PsychError PSYCHOCULUSVR1GetNextTextureHandle(void);
PsychError PSYCHOCULUSVR1CreateMirrorTexture(void);
PsychError PSYCHOCULUSVR1GetUndistortionParameters(void);
PsychError PSYCHOCULUSVR1GetStaticRenderParameters(void);
PsychError PSYCHOCULUSVR1GetTrackersState(void);
PsychError PSYCHOCULUSVR1StartRender(void);
PsychError PSYCHOCULUSVR1EndFrameRender(void);
PsychError PSYCHOCULUSVR1PresentFrame(void);
PsychError PSYCHOCULUSVR1GetEyePose(void);
PsychError PSYCHOCULUSVR1VRAreaBoundary(void);
PsychError PSYCHOCULUSVR1TestVRBoundary(void);
PsychError PSYCHOCULUSVR1TestVRBoundaryPoint(void);
PsychError PSYCHOCULUSVR1RecenterTrackingOrigin(void);
PsychError PSYCHOCULUSVR1TrackingOriginType(void);
PsychError PSYCHOCULUSVR1GetPerformanceStats(void);
PsychError PSYCHOCULUSVR1HapticPulse(void);
PsychError PSYCHOCULUSVR1FloatsProperty(void);
PsychError PSYCHOCULUSVR1FloatProperty(void);
PsychError PSYCHOCULUSVR1StringProperty(void);

#ifdef  __cplusplus
}
#endif

//end include once
#endif
