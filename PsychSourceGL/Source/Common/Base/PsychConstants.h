/*
        Psychtoolbox3/Source/Common/PsychConstants.h            
  
        AUTHORS:
 
        Allen.Ingling@nyu.edu                 awi 
        mario.kleiner@tuebingen.mpg.de        mk
 
        PLATFORMS:      All

        HISTORY:
        7/16/02         awi     Split off from PsychProject.h
  
        PROJECTS
        7/16/02         awi     MacOS9 Screen.  

        DESCRIPTION:
        PsychConstants.h defines constants which abstract out platform-specific types.  Note
        That all such abstract types are prefixed with "psych_" as in "psych_uint32".   

        This file also includes platform-specific macros.
*/

//begin include once 
#ifndef PSYCH_IS_INCLUDED_PsychConstants
#define PSYCH_IS_INCLUDED_PsychConstants

#define __STDC_LIMIT_MACROS 1

#if PSYCH_SYSTEM == PSYCH_WINDOWS
// Try to cut down compile time on Windows by only including important headers:
#define WINDOWS_LEAN_AND_MEAN
// Otherwise no math constants like M_PI defined in math.h et al. under MSVC:
#define _USE_MATH_DEFINES
#endif

//bring in the standard c and system headers 
#include "PsychIncludes.h"

#if PSYCH_LANGUAGE == PSYCH_MATLAB
    #undef printf
    #define printf mexPrintf
#endif

//platform dependent macro defines 
#if PSYCH_SYSTEM == PSYCH_WINDOWS
    // Define snprintf as _snprintf to take Microsoft brain-damage into account:
    #ifndef snprintf
    #define snprintf _snprintf
    #endif

    // MSVC wants _strdup instead of strdup:
    #ifdef _MSC_VER
    #define strdup _strdup
    #endif

#elif PSYCH_SYSTEM == PSYCH_OSX
#elif PSYCH_SYSTEM == PSYCH_LINUX
#endif 

#ifndef FALSE
#define FALSE   0
#define TRUE    1
#endif

#ifndef false
#define false FALSE
#endif

#ifndef true
#define true TRUE
#endif

// Define our own base psych_bool type psych_bool to be an unsigned char,
// i.e., 1 byte per value on each platform:
typedef unsigned char		psych_bool;

//abstract up simple data types. 
#if PSYCH_SYSTEM == PSYCH_LINUX
        typedef int64_t                         psych_int64;
        typedef unsigned long long              psych_uint64;
        typedef unsigned int                    psych_uint32;
        typedef unsigned char                   psych_uint8;
        typedef unsigned short                  psych_uint16;
        typedef char                            Str255[256];

        // We don't have these types for Linux, so we provide a little hack to
        // make the compiler happy:
        typedef int CGDisplayCount;

        #ifdef PTB_USE_WAYLAND
        // CGDirectDisplayID is typedef'd to a void* connection handle for now:
        typedef void* CGDirectDisplayID;
        #else
        // CGDirectDisplayID is typedef'd to a X11 display connection handle:
        typedef Display* CGDirectDisplayID;
        #endif
        typedef int CGDisplayErr;
        typedef unsigned int CFDictionaryRef;

        // Datatype for condition variables:
        typedef pthread_cond_t  psych_condition;
        // Datatype for Mutex Locks:
        typedef pthread_mutex_t psych_mutex;
        // Datatype for threads:
        typedef pthread_t   psych_thread;
        typedef pthread_t   psych_threadid;
#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS
        typedef LONGLONG                        psych_int64;
        typedef ULONGLONG                       psych_uint64;
        typedef DWORD                           psych_uint32;
        typedef BYTE                            psych_uint8;
        typedef WORD                            psych_uint16;

        // The Microsoft Visual C compiler doesn't know about the __func__ keyword, but it knows __FUNCTION__ instead:
        #ifdef _MSC_VER
        #define __func__ __FUNCTION__
        #endif

        // Hack to make compiler happy:
        typedef int CGDisplayCount;
        typedef HDC CGDirectDisplayID;
        typedef int CGDisplayErr;

        // Define missing types on Windows:
        typedef char         Str255[256];
        typedef unsigned int CFDictionaryRef;

        // Datatype for condition variables:
        typedef HANDLE  psych_condition;
        // Datatype for Mutex Locks:
        typedef CRITICAL_SECTION  psych_mutex;
        // Datatype for threads:
        typedef struct psych_threadstruct {
          HANDLE  handle;               // Handle to actual thread (NULL == Invalid).
          DWORD   threadId;             // Unique numeric id (0 = Invalid.)
          psych_condition terminateReq; // Condition/Event object to signal the request for termination.
          HANDLE  taskHandleMMCS;       // Handle to task for MMCSS scheduled thread, NULL otherwise.
        } psych_threadstruct;
        typedef struct psych_threadstruct*  psych_thread;

        typedef psych_uint32  psych_threadid;

#elif PSYCH_SYSTEM == PSYCH_OSX
        typedef UInt8         psych_uint8;
        typedef UInt16        psych_uint16;
        typedef UInt32        psych_uint32;
        typedef unsigned long long  psych_uint64;
        typedef long long     psych_int64;

        // Datatype for condition variables:
        typedef pthread_cond_t  psych_condition;
        // Datatype for Mutex Locks:
        typedef pthread_mutex_t psych_mutex;
        // Datatype for threads:
        typedef pthread_t     psych_thread;
        typedef pthread_t     psych_threadid;
#endif

#if PSYCH_LANGUAGE == PSYCH_MATLAB
        typedef const mxArray CONSTmxArray;
        #define PsychGenericScriptType mxArray
        typedef mxLogical PsychNativeBooleanType; 
#endif

#define EXP
#ifndef APIENTRY
#define APIENTRY
#endif

//end include once
#endif
