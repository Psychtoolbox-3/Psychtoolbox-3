/*
    PsychPlatform.h

    PLATFORMS:

        Only OS X

    AUTHORS:

        Allen Ingling        awi        Allen.Ingling@nyu.edu

    HISTORY:

        7/16/02        awi        Wrote it.  Split off from PsychPlatform.h
        11/3/04        awi        Corrected errors in comments and improved DESCRIPTION and TO DO comments.

    DESCRIPTION:

        PsychPlatform.h contains constant definitions asserting conditions
        specific to the the OS X version of Screen.  If you
        use a conditional macro wich branches according to one of the constants
        in this file, then that conditional macro does not belong here.

        Note that this should be the ONLY Psychtoolbox file which is conditionally
        included by platform.  Exceptions to this rule might be:
        -project.def files included in VC++ files
        -StdAfx pre-compiled header files included in VC++ files.

    TO DO:

        For now the project path specifies which platform version of this file to
        include, with each version of this file defining constants which identify
        the platform.  A smarter way to do this would be to use a single version
        of this file which detects contants defined within the compilation
        environment, for example as a flag passed to the compiler.
*/

#include "PsychPlatformConstants.h"

//these control build switches
#define PSYCH_SYSTEM        PSYCH_OSX

#ifndef PSYCH_LANGUAGE
#define PSYCH_LANGUAGE      PSYCH_MATLAB
#endif

#define PSYCH_DEBUG         PSYCH_ON
