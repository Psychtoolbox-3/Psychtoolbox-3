/*
        PsychToolbox3/Source/Common/PsychPortAudio/PsychPortAudio.h

        PLATFORMS:    All

        AUTHORS:

        Mario Kleiner   mk      mario.kleiner.de@gmail.com

        HISTORY:

        21.03.07        mk      wrote it.

        NOTES:

        This needs the portaudio.h header file, which should be located in /usr/include/portaudio.h
        after a standard ./configure; make; make install installation on Linux & OS/X.
*/

//begin include once
#ifndef PSYCH_IS_INCLUDED_PsychPortAudio
#define PSYCH_IS_INCLUDED_PsychPortAudio

//project includes
#include "Psych.h"
#include "PsychTimeGlue.h"
#include "portaudio.h"

// Internal helper functions:

// Module exit function: Stop and close all audio streams, terminate PortAudio...
PsychError PsychPortAudioExit(void);

// Init function: Can be called anytime to make sure that PortAudio is online...
void PsychPortAudioInitialize(void);

// Show command overview:
PsychError PSYCHPORTAUDIODisplaySynopsis(void);
const char** InitializeSynopsis(void);

//function prototypes:

// Version:
PsychError MODULEVersion(void);
// Level of verbosity:
PsychError PSYCHPORTAUDIOVerbosity(void);

// Open audio device:
PsychError PSYCHPORTAUDIOOpen(void);
// Open virtual audio slave device:
PsychError PSYCHPORTAUDIOOpenSlave(void);
// Close audio device, shutdown PortAudio if last device is closed:
PsychError PSYCHPORTAUDIOClose(void);
// Fill audio output buffer with data:
PsychError PSYCHPORTAUDIOFillAudioBuffer(void);
// Refill audio output buffer with data:
PsychError PSYCHPORTAUDIORefillBuffer(void);
// Start or schedule start of device:
PsychError PSYCHPORTAUDIOStartAudioDevice(void);
// Reschedule start of audio device:
PsychError PSYCHPORTAUDIORescheduleStart(void);
// Stop device immediately:
PsychError PSYCHPORTAUDIOStopAudioDevice(void);
// Enumerate audio devices:
PsychError PSYCHPORTAUDIOGetDevices(void);
// Return status of device:
PsychError PSYCHPORTAUDIOGetStatus(void);
// Set a manual bias for the latencies we operate on.
PsychError PSYCHPORTAUDIOLatencyBias(void);
// Retrieve buffer with captured audio data:
PsychError PSYCHPORTAUDIOGetAudioData(void);
// Select general run mode for audio device:
PsychError PSYCHPORTAUDIORunMode(void);
// Select sample loop for audio device:
PsychError PSYCHPORTAUDIOSetLoop(void);
// Set tuneable low-level engine parameters:
PsychError PSYCHPORTAUDIOEngineTunables(void);
// Get count of open audio devices:
PsychError PSYCHPORTAUDIOGetOpenDeviceCount(void);
// Alloc/Dealloc Enable/Disable playback schedules:
PsychError PSYCHPORTAUDIOUseSchedule(void);
// Add sot to playback schedule:
PsychError PSYCHPORTAUDIOAddToSchedule(void);
// Create and fill dynamic audio buffer:
PsychError PSYCHPORTAUDIOCreateBuffer(void);
// Delete dynamic audio buffer:
PsychError PSYCHPORTAUDIODeleteBuffer(void);
// Change device opMode at runtime:
PsychError PSYCHPORTAUDIOSetOpMode(void);
// Enable/Disable or reconfigure direct input monitoring:
PsychError PSYCHPORTAUDIODirectInputMonitoring(void);
// Set per-device volume:
PsychError PSYCHPORTAUDIOVolume(void);
//end include once
#endif
