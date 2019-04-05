/*
    PsychToolbox3/Source/Common/PsychPortAudio/RegisterProject.c

    PROJECTS: PsychPortAudio only.

    AUTHORS:

    Mario Kleiner   mk      mario.kleiner.de@gmail.com

    PLATFORMS:      All.

    HISTORY:

    21.3.2007       mk      Created.

*/

#include "Psych.h"
#include "RegisterProject.h"
#include "PsychPortAudio.h"

PsychError PsychModuleInit(void)
{
    //register the project exit function
    PsychErrorExit(PsychRegisterExit(&PsychPortAudioExit));

    // Register the project function which is called when the module
    // is invoked with no arguments:
    PsychErrorExit(PsychRegister(NULL,  &PSYCHPORTAUDIODisplaySynopsis));

    // Report the version
    PsychErrorExit(PsychRegister("Version",  &MODULEVersion));

    // Register the module name
    PsychErrorExit(PsychRegister("PsychPortAudio", NULL));

    // Register synopsis and named subfunctions.
    PsychErrorExit(PsychRegister("Verbosity", &PSYCHPORTAUDIOVerbosity));
    PsychErrorExit(PsychRegister("Open", &PSYCHPORTAUDIOOpen));
    PsychErrorExit(PsychRegister("OpenSlave", &PSYCHPORTAUDIOOpenSlave));
    PsychErrorExit(PsychRegister("Close", &PSYCHPORTAUDIOClose));
    PsychErrorExit(PsychRegister("Start", &PSYCHPORTAUDIOStartAudioDevice));
    PsychErrorExit(PsychRegister("RescheduleStart", &PSYCHPORTAUDIORescheduleStart));
    PsychErrorExit(PsychRegister("Stop", &PSYCHPORTAUDIOStopAudioDevice));
    PsychErrorExit(PsychRegister("FillBuffer", &PSYCHPORTAUDIOFillAudioBuffer));
    PsychErrorExit(PsychRegister("RefillBuffer", &PSYCHPORTAUDIORefillBuffer));
    PsychErrorExit(PsychRegister("GetDevices", &PSYCHPORTAUDIOGetDevices));
    PsychErrorExit(PsychRegister("GetStatus", &PSYCHPORTAUDIOGetStatus));
    PsychErrorExit(PsychRegister("LatencyBias", &PSYCHPORTAUDIOLatencyBias));
    PsychErrorExit(PsychRegister("GetAudioData", &PSYCHPORTAUDIOGetAudioData));
    PsychErrorExit(PsychRegister("RunMode", &PSYCHPORTAUDIORunMode));
    PsychErrorExit(PsychRegister("SetLoop", &PSYCHPORTAUDIOSetLoop));
    PsychErrorExit(PsychRegister("EngineTunables", &PSYCHPORTAUDIOEngineTunables));
    PsychErrorExit(PsychRegister("GetOpenDeviceCount", &PSYCHPORTAUDIOGetOpenDeviceCount));
    PsychErrorExit(PsychRegister("UseSchedule", &PSYCHPORTAUDIOUseSchedule));
    PsychErrorExit(PsychRegister("AddToSchedule", &PSYCHPORTAUDIOAddToSchedule));
    PsychErrorExit(PsychRegister("CreateBuffer", &PSYCHPORTAUDIOCreateBuffer));
    PsychErrorExit(PsychRegister("DeleteBuffer", &PSYCHPORTAUDIODeleteBuffer));
    PsychErrorExit(PsychRegister("SetOpMode", &PSYCHPORTAUDIOSetOpMode));
    PsychErrorExit(PsychRegister("DirectInputMonitoring", &PSYCHPORTAUDIODirectInputMonitoring));
    PsychErrorExit(PsychRegister("Volume", &PSYCHPORTAUDIOVolume));

    // Setup synopsis help strings:
    InitializeSynopsis();   //Scripting glue won't require this if the function takes no arguments.

    // Setup module author:
    PsychSetModuleAuthorByInitials("mk");

    // Call wait-routine for 0.1 secs: This to initialize the time glue on MS-Windows,
    // so the first call to a timing function won't delay:
    PsychWaitIntervalSeconds(0.1);

    // Startup finished.
    return(PsychError_none);
}
