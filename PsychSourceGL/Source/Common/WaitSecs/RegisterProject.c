/*
 * PsychToolbox3/Source/Common/WaitSecs/RegisterProject.c
 *
 * PROJECTS: WaitSecs only.
 *
 * AUTHORS:
 *
 *      Allen.Ingling@nyu.edu           awi
 *      mario.kleiner.de@gmail.com      mk
 *
 * PLATFORMS: All
 *
 * HISTORY:
 * 8/20/02      awi Created.
 * 1/20/02      awi Created derived the GetSecs version from the Screen version.
 * 1/2/08       mk  Add subfunction for waiting until absolute time, and return of wakeup time.
 *
 */

//begin include once
#include "Psych.h"
#include "WaitSecs.h"

PsychError PsychModuleInit(void)
{
    // Register the project exit function
    PsychErrorExit(PsychRegisterExit(NULL));

    // Register the project function which is called when the module
    // is invoked with no named subfunction:
    PsychErrorExit(PsychRegister(NULL,  &WAITSECSWaitSecs));

    // Wait until specific deadline:
    PsychErrorExit(PsychRegister("UntilTime", &WAITSECSWaitUntilSecs));
    PsychErrorExit(PsychRegister("YieldSecs", &WAITSECSYieldSecs));

    // Report the version
    PsychErrorExit(PsychRegister("Version", &MODULEVersion));

    // Register the module name
    PsychErrorExit(PsychRegister("WaitSecs", NULL));

    PsychSetModuleAuthorByInitials("awi");
    PsychSetModuleAuthorByInitials("mk");

    // Register synopsis and named subfunctions.
    InitializeSynopsis();

    return(PsychError_none);
}
