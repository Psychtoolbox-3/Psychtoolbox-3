/*
  PsychToolbox3/Source/Common/PsychInit.c

  AUTHORS:

    Allen.Ingling@nyu.edu                 awi
    mario.kleiner.de@gmail.com            mk

  PLATFORMS: All

  PROJECTS:  All

  HISTORY:

  08/25/02  awi     wrote it.
  03/24/08  mk      Add call to PsychExitTimeGlue() and some cleanup.

  DESCRIPTION:

    Holds the master initialization and shutdown for the Psychtoolbox function.
    Sub-section hold their own inits which are called by the master
    init function contained herein.  After the Psychtoolbox library
    has initialzed itself in then invokes the project init which
    must be named PsychModuleInit().

    Each PsychToolbox module should register its subfuctions within
    its PsychModuleInit().

*/

#include "Psych.h"

PsychError PsychInit(void)
{
    // Assign distinctive name to master thread:
    PsychSetThreadName("PTB mainthread");

    // Record thread id of master thread:
    PsychInitMasterThreadId();

    // First init Psychtoolbox libraries:
    InitPsychError();
    InitPsychAuthorList();
    PsychInitTimeGlue();

    // Registration of the Psychtoolbox exit function is
    // done in PsychScriptingGlueXXX.c because how that is done is
    // specific to the scripting environment. Note that registration
    // of the project exit function is done within the project init.

    // Then call the project init.
    PsychModuleInit();

    return(PsychError_none);
}

/* PsychExit is the function invoked last before the module is
   purged. It is abstracted to be unspecific to the scripting
   language. The XXX language-specific versions are named
   PsychExitGlue() and they are found in PsychScriptingXXX.c
 */
PsychError PsychExit(void)
{
    PsychFunctionPtr projectExit;

    projectExit = PsychGetProjectExitFunction();
    if (projectExit != NULL) (*projectExit)();

    // Put whatever cleanup of the Psychtoolbox is required here.
    PsychExitTimeGlue();

    // Reset / Clear function and module name registry:
    PsychResetRegistry();

    return(PsychError_none);
}
