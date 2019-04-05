/*
 *    PsychSourceGL/Source/Common/GetSecs/RegisterProject.c
 *
 *    PROJECTS:
 *
 *        GetSecs only
 *
 *    AUTHORS:
 *
 *        Allen.Ingling@nyu.edu           awi
 *        mario.kleiner.de@gmail.com      mk
 *
 *    PLATFORMS:
 *
 *        All.
 *
 *    HISTORY:
 *
 *        1/20/02         awi     Derived the GetSecs project from Screen .
 *        4/6/05          awi     Updated header comments.
 */

//begin include once

#include "Psych.h"
#include "GetSecs.h"

PsychError PsychModuleInit(void)
{
    //register the project exit function
    PsychErrorExit(PsychRegisterExit(NULL));

    //register the project function which is called when the module
    //is invoked with no arguments:
    PsychErrorExit(PsychRegister(NULL,  &GETSECSGetSecs));

    PsychErrorExit(PsychRegister("Version",  &MODULEVersion));
    PsychErrorExit(PsychRegister("AllClocks",  &GETSECSAllClocks));

    //register the module name
    PsychErrorExit(PsychRegister("GetSecs", NULL));

    PsychSetModuleAuthorByInitials("awi");
    PsychSetModuleAuthorByInitials("mk");

    // Register synopsis and named subfunctions.
    InitializeSynopsis();

    return(PsychError_none);
}
