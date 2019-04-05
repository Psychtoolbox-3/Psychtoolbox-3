/*
    PsychToolbox3/Source/Common/WaitSecs/WaitSecs.h

    PLATFORMS:  All

    AUTHORS:

    Allen Ingling   awi     Allen.Ingling@nyu.edu

    HISTORY:
    1/20/02         awi     wrote it.
*/

//begin include once
#ifndef PSYCH_IS_INCLUDED_WaitSecs
#define PSYCH_IS_INCLUDED_WaitSecs

//project includes
#include "Psych.h"
#include "PsychTimeGlue.h"

//function prototypes
const char** InitializeSynopsis(void);
PsychError MODULEVersion(void);         //MODULEVersion.c
PsychError WAITSECSWaitSecs(void);
PsychError WAITSECSWaitUntilSecs(void);
PsychError WAITSECSYieldSecs(void);

//end include once
#endif
