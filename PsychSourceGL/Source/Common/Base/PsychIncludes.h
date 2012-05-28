/*
  Psychtoolbox2/Source/Common/PsychIncludes.h		
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
  
  PLATFORMS: All 
  
  PROJECTS:
  07/16/02	awi		Screen on MacOS9
   

  HISTORY:
  07/16/02  awi		Pulled out of PsychPlatform.h
  11/15/02  awi		Added includes for OSX.  
  
  DESCRIPTION:
  PsychIncludes.h includes all C, system, and language binding 
  header files which a Psychtoolbox library 
  would require.
  
  This file should ONLY be included by PsychConstants.h 

  T0 DO: 
	  

*/

#ifndef PSYCH_IS_INCLUDED_PsychIncludes
#define PSYCH_IS_INCLUDED_PsychIncludes

#include "PsychPlatform.h"

//Includes dependent on runtime environment:
#if PSYCH_LANGUAGE == PSYCH_MATLAB
    // Include File with Matlab interface API definitions:
    #include "mex.h"
#endif

// Platform independent include for glew: This is a catch-all
// for all OpenGL definitions and functions, currently up to
// OpenGL 2.1:
#if defined(PTBMODULE_Screen) || defined(PTBMODULE_FontInfo)
#include "../Screen/glew.h"
#endif

//platform dependent includes stage 1
#if PSYCH_SYSTEM == PSYCH_LINUX
	// We need __USE_UNIX98, so pthread.h defines/supports the mutex policy
	// attribute for mutex priority inheritance for our realtime threads.
	// For some reason this gets undefined in mex.h at least when building
	// on Octave 3.2.4. Scary scary...
	#ifndef __USE_UNIX98
	#define __USE_UNIX98
	// For testing only: #warning __USE_UNIX98 undefined. Redefining it.
	#endif

	#include <stdint.h>
	// This is the new glew include for GLX extension support:
	#include "../Screen/glxew.h"
	
        #include <X11/Xlib.h>
        #include <X11/keysym.h>
        #include <sys/time.h>
        #include <unistd.h>
        #include <stdarg.h>
	#include <pthread.h>

	#ifndef _POSIX_THREAD_PRIO_INHERIT
	#error This build system does not support pthread_mutexattr_setprotocol(&attr, PTHREAD_PRIO_INHERIT); Fix your build system!
	#endif
#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS
	// Need to define #define _WIN32_WINNT as 0x0400 so we can use TryEnterCriticalSection() call for PsychTryLockMutex() implementation.
    // We set WINVER and _WIN32_WINNT to 0x0500, which requires Windows 2000 or later as target system:
	#define _WIN32_WINNT 0x0500
    #define WINVER       0x0500
	#include <windows.h>
	#include "../Screen/wglew.h"

    #ifndef INT64_MAX
    #define INT64_MAX _I64_MAX
    #endif
    
    #ifndef INT32_MAX
    #define INT32_MAX _I32_MAX
    #endif

#elif PSYCH_SYSTEM == PSYCH_OS9
	#include <Types.h>
	#include <MacTypes.h>
	#include <Events.h>
	#include <stdarg.h>
	#include <string.h>
	#include <ctype.h>
	#include "gl.h"
	#include "glu.h"
#elif PSYCH_SYSTEM == PSYCH_OSX
	//includes for Apple OpenGL
    #include <stdint.h>
    #include <sys/time.h>
	#include <pthread.h>

    // If we are included from PsychScriptingGlue.cc, which
    // is indicated by PTBINSCRIPTINGGLUE, then we must only
    // include MacTypes.h, not the rest of Carbon, ApplicationServices
    // et al. -> Some of the Octave header files conflict with Apple
    // system header files --> do not include the Apple headers when
    // building Scriptingglue - they are not needed but would prevent
    // the build.
    #ifndef PTBINSCRIPTINGGLUE
    // File included during compile of some PTB file other than ScriptingGlue...
    #include <Carbon/Carbon.h>
    #include <CoreServices/CoreServices.h>
    #include <ApplicationServices/ApplicationServices.h>
    #include <CoreAudio/HostTime.h>
    #include <CoreVideo/CoreVideo.h>
    #else
    // File included from ScriptingGlue - only import minimal amount of headers...
    // MK TODO FIXME! OCTAVE!
    #include <CoreServices/CoreServices.h>
    // Does not work: #include <CarbonCore.framework/MacTypes.h>
    // Used to work on Tiger/32-Bit, but does not work now: #include <CoreServices/Frameworks/CarbonCore.framework/Headers/MacTypes.h>
    #endif

    #if defined(PTBMODULE_Screen) || defined(PTBMODULE_FontInfo)
    #include <OpenGL/OpenGL.h>
    #endif
#endif 

//C standard library headers
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <float.h>

#if PSYCH_LANGUAGE == PSYCH_OCTAVE
// File with pseudo MATLAB interface API definitions:
// Emulates missing Matlab functions and structures...
// Stored in the Source/Octave subfolder...
    #include <octavemex.h>
#endif

//SDL_getenv.h includes a macro which redefines getenv() to be
//SDL getenv().   This then conflicts with the statement "using 
//std::getenv" within stdlib.h.  To avoid the conflict, include 
//SDL.h only after the standard includes above. 
//platform dependent includes stage 2
#if PSYCH_SYSTEM == PSYCH_WINDOWS

#elif PSYCH_SYSTEM == PSYCH_OS9
	#include "SDL.h"
#endif

//end include once
#endif
