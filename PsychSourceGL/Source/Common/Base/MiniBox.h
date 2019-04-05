/*
 * MiniBox.h
 *
 * AUTHORS:
 *
 * Allen Ingling        awi     Allen.W.Ingling@nyu.edu
 * Mario Kleiner        mk      mario.kleiner.de@gmail.com
 *
 * PLATFORMS: All
 *
 * HISTORY:
 *  11/25/01        awi     Created file.
 *  07/16/02        awi     Included PsychConstants.h and added include once.
 *   3/19/10        mk      Cosmetic and make 64-bit clean.
 *
 */

// begin include once
#ifndef PSYCH_IS_INCLUDED_MiniBox
#define PSYCH_IS_INCLUDED_MiniBox

#include "Psych.h"

// define mex-derived data types.
#if PSYCH_LANGUAGE == PSYCH_MATLAB
typedef void (*MexFunctionPtr)(int nlhs, mxArray *plhs[], int nrhs, CONSTmxArray *prhs[]);
#endif

extern int PrintfExit(const char *format,...);
char *BreakLines(char *string,long lineLength);
psych_bool PsychMatch(char *s1,char *s2);
char *int2str(psych_int64 num);
size_t PsychIndexElementFrom2DArray(size_t mDim/*|Y|*/, size_t nDim/*|X|*/, size_t m/*y*/, size_t n/*x*/);
size_t PsychIndexElementFrom3DArray(size_t mDim/*|Y|*/, size_t nDim/*|X|*/, size_t pDim/*|Z|*/, size_t m/*y*/, size_t n/*x*/, size_t p/*z*/);
size_t PsychIndexPlaneFrom3DArray(size_t mDim, size_t nDim, size_t pDim, size_t planeIndex);
psych_int64 maxInt(psych_int64 a, psych_int64 b);  //might be better to define max as a macro...
psych_bool PsychIsIntegerInDouble(double *value);
psych_bool PsychIsInteger64InDouble(double *value);

psych_bool PsychIsPsychMatchCaseSensitive(void);
void PsychSetPsychMatchCaseSenstive(psych_bool arg);

void PsychInitMasterThreadId(void);
psych_threadid PsychGetMasterThreadId(void);
psych_bool PsychIsMasterThread(void);

// end include once
#endif
