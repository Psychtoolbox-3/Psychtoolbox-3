#ifndef _MOGLTYPE_H
#define _MOGLTYPE_H

/*
 * mogltype.h -- common definitions for gl/glu and glm modules
 *
 * 09-Dec-2005 -- created (RFM)
 * 24-Mar-2011 -- Make 64-bit clean (MK).
 * 27-Mar-2011 -- Remove obsolete and totally bitrotten Octave-2 support (MK).
 * 01-Jul-2012 -- Kill Matlab R11 support (MK).
 * 17-Oct-2015 -- For Octave-4 we need the gluUnProject4 definition like on Matlab (MK).
 */

#define PSYCH_MATLAB 0
#define PSYCH_LANGUAGE PSYCH_MATLAB

// Need to include uchar.h before mex.h for definition of
// char16_t, which seems to be missing on recent compilers
// by default. Otherwise compile failure at least under
// Matlab R2012a with current Linux or OSX. On Linux this
// solves it. On OSX 10.9 we instead -Dchar16_t=UINT16_T in mexopts.sh
// and are fine again.
#ifdef LINUX
#include <uchar.h>
#endif

// Include mex.h with MEX - API definition for Matlab:
#include "mex.h"

#include <stdio.h>
#include <limits.h>

#define printf mexPrintf
#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifdef MACOSX
/* Shut up the OpenGL related deprecation warnings on Apples trainwreck: */
#define GL_SILENCE_DEPRECATION
#endif

/* glew.h is part of GLEW library for automatic detection and binding of
   OpenGL core functionality and extensions.
 */
#include "GL/glew.h"

/* Includes specific to MacOS-X version of mogl: */
#ifdef MACOSX

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#ifdef MATLAB_MEX_FILE 
#include <AGL/agl.h>
#endif
#define CALLCONV

#endif

/* Includes specific to GNU/Linux version of mogl: */
#ifdef LINUX

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "GL/glxew.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#define CALLCONV

#endif

/* Includes specific to M$-Windows version of mogl: */
#ifdef WINDOWS

#include <stdlib.h>
#include <string.h>
#include "GL/wglew.h"
#include <GL/glut.h>
#define CALLCONV __stdcall

/* Hacks to get Windows version running: */
double gluCheckExtension(const GLubyte* a, const GLubyte* b);
double gluBuild1DMipmapLevels(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, void* a9);
double gluBuild2DMipmapLevels(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9, void* a10);
double gluBuild3DMipmapLevels(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9, double a10, void* a11);
double gluBuild3DMipmaps(double a1, double a2, double a3, double a4, double a5, double a6, double a7, void* a8);
double gluUnProject4(double a1, double a2, double a3, double a4, double* a5, double* a6, int* a7, double a8, double a9, double* a10, double* a11, double* a12, double* a13);
#endif

// Mapping of scalar buffer offset value (in units of bytes) to an
// equivalent memory void*.
void* moglScalarToPtrOffset(const mxArray *m);

// Function prototype for exit with printf(), like mogl_glunsupported...
void mogl_printfexit(const char* str);

// Function prototype for error handler for unsupported gl-Functions.
void mogl_glunsupported(const char* fname);

// typedef for command map entries
typedef struct cmdhandler {
    char *cmdstr;
    void (*cmdfn)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
} cmdhandler;

// Definitions for GLU-Tesselators:
#define MOGLDEFMYTESS mogl_tess_struct* mytess = (mogl_tess_struct*) PsychDoubleToPtr(mxGetScalar(prhs[0]));
// Old style: 32-bit systems only: #define MOGLDEFMYTESS mogl_tess_struct* mytess = (mogl_tess_struct*) (unsigned int) mxGetScalar(prhs[0]);

#define MAX_TESSCBNAME 32

typedef struct mogl_tess_struct {
    GLUtesselator*  glutesselator;
    int             userData;
    double          polygondata;
    double*         destructBuffer;
    size_t          destructSize;
    size_t          maxdestructSize;
    size_t          destructCount;
    size_t          nrElements;
    char            nGLU_TESS_BEGIN[MAX_TESSCBNAME];
    char            nGLU_TESS_BEGIN_DATA[MAX_TESSCBNAME];
    char            nGLU_TESS_EDGE_FLAG[MAX_TESSCBNAME];
    char            nGLU_TESS_EDGE_FLAG_DATA[MAX_TESSCBNAME];
    char            nGLU_TESS_VERTEX[MAX_TESSCBNAME];
    char            nGLU_TESS_VERTEX_DATA[MAX_TESSCBNAME];
    char            nGLU_TESS_END[MAX_TESSCBNAME];
    char            nGLU_TESS_END_DATA[MAX_TESSCBNAME];
    char            nGLU_TESS_COMBINE[MAX_TESSCBNAME];
    char            nGLU_TESS_COMBINE_DATA[MAX_TESSCBNAME];
    char            nGLU_TESS_ERROR[MAX_TESSCBNAME];
    char            nGLU_TESS_ERROR_DATA[MAX_TESSCBNAME];    
} mogl_tess_struct;

// Our own little buffer memory manager: This code is adapted from
// Psychtoolbox's PsychMemory.h/.c routines:

// Definition of unsigned int 64 bit datatype for Windows vs. Unix
#ifndef WINDOWS
typedef unsigned long long int psych_uint64;
#else
typedef ULONGLONG psych_uint64;
#endif

// Definition of how a memory pointer is encoded in the runtime env.
// Currently we encode pointers as double's:
#define psych_RuntimePtrClass mxDOUBLE_CLASS
typedef double psych_RuntimePtrType;

// Convert a double value (which encodes a memory address) into a ptr:
void*  PsychDoubleToPtr(volatile double dptr);

// Convert a memory address pointer into a double value:
double PsychPtrToDouble(void* ptr);

// Return the size of the buffer pointed to by ptr in bytes.
// CAUTION: Only works with buffers allocated via PsychMallocTemp()
// or PsychCallocTemp(). Will segfault, crash & burn with other pointers!
size_t PsychGetBufferSizeForPtr(void* ptr);

// Allocate and zero-out n elements of memory, each size bytes big.
void *PsychCallocTemp(size_t n, size_t size, int mlist);

// Allocate n bytes of memory:
void *PsychMallocTemp(size_t n, int mlist);

// PsychFreeTemp frees memory allocated with PsychM(C)allocTemp().
// This is not strictly needed as our memory manager will free
// the memory anyway when this file is flushed:
void PsychFreeTemp(void* ptr, int mlist);

// Master cleanup routine: Frees all allocated memory.
void PsychFreeAllTempMemory(int mlist);

#endif
