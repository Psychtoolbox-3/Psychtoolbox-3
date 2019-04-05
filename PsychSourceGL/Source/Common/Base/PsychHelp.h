/*
  PsychToolbox3/Source/Common/PsychHelp.h

  AUTHORS:

    Allen.Ingling@nyu.edu           awi
    mario.kleiner.de@gmail.com      mk

  PLATFORMS: All

  PROJECTS:

    08/19/02    awi All

*/

//begin include once
#ifndef PSYCH_IS_INCLUDED_PsychHelp
#define PSYCH_IS_INCLUDED_PsychHelp

#include "Psych.h"

void PsychSetGiveHelp(void);
void PsychClearGiveHelp(void);
psych_bool PsychIsGiveHelp(void);
void PsychPushHelp(char *functionUse, char *functionSynopsis, char *functionSeeAlso);
void PsychGiveHelp(void);
void PsychGiveUsage(void);
void PsychGiveUsageExit(void);
char *PsychGetFunctionName(void);
void PsychOneShotReturnHelp(void);
const char* PsychBuildSynopsisString(const char* modulename);

//end include once
#endif
