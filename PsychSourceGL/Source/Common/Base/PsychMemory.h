/*
  Psychtoolbox/Source/Common/PsychMemory.h

  AUTHORS:

  Allen.Ingling@nyu.edu             awi
  mario.kleiner.de@gmail.com        mk

  PLATFORMS: All

  PROJECTS: All

  HISTORY:

  09/04/02  awi     Wrote it.
  05/10/06  mk      Added our own allocator for Octave-Port.
  03/19/11  mk      Make 64-bit clean.  

*/

//begin include once
#ifndef PSYCH_IS_INCLUDED_PsychMemory
#define PSYCH_IS_INCLUDED_PsychMemory

// Convert a double value (which encodes a memory address) into a ptr:
void* PsychDoubleToPtr(volatile double dptr);

// Convert a memory address pointer into a double value:
double PsychPtrToDouble(void* ptr);

//allocate memory which is valid between until control returns to
//the calling module:
void *PsychCallocTemp(size_t n, size_t size);
void *PsychMallocTemp(size_t n);

//free memory
#if PSYCH_LANGUAGE == PSYCH_MATLAB
    #define PsychFreeTemp  mxFree
#else

// Our own allocator / memory manager is used.

// PsychFreeTemp frees memory allocated with PsychM(C)allocTemp().
// This is not strictly needed as our memory manager will free
// the memory anyway when returning control to Matlab/Octave et al.
void PsychFreeTemp(void* inptr);

// Master cleanup routine: Frees all allocated memory.
void PsychFreeAllTempMemory(void);

#endif

//allocate memory which is valid while the module is loaded

/*
    For Matlab this should 
    1- allocate memory with mxCalloc or mxMalloc
    2- call mexMakeMemoryPersistent on the pointer
    3- store the pointer in a linked list of pointers 
        to memory which should be freed with MexAtExit
        when the module is purged

    do something similar on cleanup for environments where we use
    genuine malloc and calloc and free instead of the matlab 
    versions.
*/

//end include once
#endif
