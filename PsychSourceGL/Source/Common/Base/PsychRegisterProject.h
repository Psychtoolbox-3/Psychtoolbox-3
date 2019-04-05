/*
  PsychToolbox2/Source/Common/PsychRegisterProject.h

  PROJECTS: All.

  AUTHORS:
    Allen.Ingling@nyu.edu             awi
    mario.kleiner.de@gmail.com        mk

  PLATFORMS: All.

  HISTORY:
    8/23/02  awi        Created.

*/

//begin include once
#ifndef PSYCH_IS_INCLUDED_PsychRegisterProject
#define PSYCH_IS_INCLUDED_PsychRegisterProject

#include "Psych.h"

#define PSYCH_MAX_FUNCTION_NAME_LENGTH 64
#define PSYCH_MAX_FUNCTIONS 256

typedef struct
{
    char name[PSYCH_MAX_FUNCTION_NAME_LENGTH+1];  // +1 for term null
    PsychFunctionPtr function;
} PsychFunctionTableEntry;

PsychError PsychDescribeModuleFunctions(void);
PsychError PsychRegister(char *name,  PsychFunctionPtr func);
PsychError PsychRegisterExit(PsychFunctionPtr exitFunc);
void PsychResetRegistry(void);
PsychFunctionPtr PsychGetProjectFunction(char *command);
char *PsychGetFunctionName(void);
char *PsychGetModuleName(void);
char *PsychGetBuildDate(void);
char *PsychGetBuildTime(void);
int PsychGetBuildNumber(void);
char *PsychGetVersionString(void);
PsychFunctionPtr PsychGetProjectExitFunction(void);

//end include once
#endif
