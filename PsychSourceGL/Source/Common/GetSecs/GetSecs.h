/*
 *    PsychSourceGL/Source/Common/GetSecs/GetSecs.h
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
#ifndef PSYCH_IS_INCLUDED_GetSecs
#define PSYCH_IS_INCLUDED_GetSecs

//project includes
#include "Psych.h"
#include "PsychTimeGlue.h"

//function prototypes
const char** InitializeSynopsis(void);
PsychError MODULEVersion(void);
PsychError GETSECSGetSecs(void);
PsychError GETSECSAllClocks(void);

//end include once
#endif
