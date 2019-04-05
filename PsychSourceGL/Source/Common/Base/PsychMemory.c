/*
    Psychtoolbox3/Source/Common/PsychMemory.c

    AUTHORS:

    Allen.Ingling@nyu.edu           awi
    mario.kleiner.de@gmail.com      mk

    PLATFORMS: All

    PROJECTS: All

    HISTORY:

    09/04/02  awi   Wrote it.
    03/19/11  mk    Make 64-bit clean.

*/

#include "Psych.h"

// Pointer to the start of our list of allocated temporary memory buffers:
static size_t* PsychTempMemHead = NULL;

// Total count of allocated memory in Bytes:
static size_t totalTempMemAllocated = 0;

// Convert a double value (which encodes a memory address) into a ptr:
void*  PsychDoubleToPtr(volatile double dptr)
{
    volatile psych_uint64* iptr = (psych_uint64*) &dptr;
    volatile psych_uint64 ival = *iptr;
    return((void*) ((size_t) ival));
}

// Convert a memory address pointer into a double value:
double PsychPtrToDouble(void* ptr)
{
    volatile psych_uint64 ival = (psych_uint64) ((size_t) ptr);
    volatile double* dptr = (double*) &ival;
    volatile double outval = *dptr;
    return(outval);
}

#if PSYCH_LANGUAGE == PSYCH_MATLAB

// If running on Matlab, we use Matlab's memory manager...
void *PsychCallocTemp(size_t n, size_t size)
{
    void *ret;

    if(NULL==(ret=mxCalloc(n, size))){
        if(size * n != 0)
        PsychErrorExitMsg(PsychError_outofMemory, NULL);
    }

    return(ret);
}

void *PsychMallocTemp(size_t n)
{
    void *ret;

    if(NULL==(ret=mxMalloc(n))){
        if(n!=0)
        PsychErrorExitMsg(PsychError_outofMemory,NULL);
    }

    return(ret);
}

#else

// If not running on Matlab, we use our own allocator...

#define PTBTEMPMEMDEC(n) totalTempMemAllocated -=(n)

// Enqueues a new record into our linked list of temp. memory buffers.
// Returns the memory pointer to be passed to rest of Psychtoolbox.
static void* PsychEnqueueTempMemory(size_t* p, size_t n)
{
    // Add current buffer-head ptr as next-pointer to our new buffer:
    *p = (size_t) PsychTempMemHead;

    // Set our buffer as new head of list:
    PsychTempMemHead = p;

    // Add allocated buffer size as 2nd element:
    p++;
    *p = n;

    // Accounting:
    totalTempMemAllocated += n;
    // printf("TEMPMALLOCED %i Bytes,  new = %i.\n", n, totalTempMemAllocated); fflush(NULL);

    // Increment p again to get real start of user-visible buffer:
    p++;

    // Return ptr:
    return((void*) p);
}

void *PsychCallocTemp(size_t n, size_t size)
{
    size_t *ret;

    // MK: This could create an overflow if product n * size is
    // bigger than length of a unsigned long int --> Only
    // happens if more than 4 GB of RAM are allocated at once.
    // --> Improbable for PTB, unless someones trying a buffer
    // overflow attack -- PTB would lose there badly anyway...
    size_t realsize = n * size + sizeof(void*) + sizeof(size_t);

    // realsize has extra bytes allocated for our little header...  
    if (NULL == (ret = (size_t*) calloc((size_t) 1, realsize))) {
        PsychErrorExitMsg(PsychError_outofMemory, NULL);
    }

    // Need to enqueue memory buffer...
    return(PsychEnqueueTempMemory(ret, realsize));
}

void *PsychMallocTemp(size_t n)
{
    size_t *ret;

    // Allocate some extra bytes for our little header...
    n = n + sizeof(void*) + sizeof(size_t);

    if (NULL == (ret = (size_t*) malloc(n))) {
        PsychErrorExitMsg(PsychError_outofMemory,NULL);
    }

    // Need to enqueue memory buffer...
    return(PsychEnqueueTempMemory(ret, n));
}

// Free a single spec'd temp memory buffer.
// TODO Note: The current implementation of our allocator
// uses a single-linked list, which has O(1) cost for
// allocating memory (Optimal!) and O(n) cost for freeing
// all allocated memory (Optimal!), but it has up to
// O(n) cost for deleting a single memory buffer as well,
// be n the total number of allocated buffers. This is
// worst-case upper bound. If PsychFreeTemp() is used a
// lot on long buffer lists, this will incur significant
// overhead! A better implementation would use a double-
// linked list or even a binary tree or hash structure,
// but for now this has to be good enough(TM).
void PsychFreeTemp(void* inptr)
{
    size_t *ptr = (size_t*) inptr;
    size_t *ptrbackup = ptr;
    size_t *next = PsychTempMemHead;
    size_t *prevptr = NULL;

    if (ptr == NULL)
        return;

    // Convert ptb supplied pointer (in)ptr into real start
    // of our buffer, including our header:
    ptr-= 2;

    if (ptr == NULL)
        return;

    if (PsychTempMemHead == ptr) {
        // Special case: ptr is first buffer in queue. Dequeue:
        PsychTempMemHead = (size_t*) *PsychTempMemHead;

        // Some accounting:
        PTBTEMPMEMDEC(ptr[1]);

        // Release it:
        free(ptr);

        return;
    }

    // ptr valid and not first buffer in queue.
    // Walk the whole buffer list until we encounter our buffer:
    while (next != NULL && next != ptr) {
        prevptr = next;
        next = (size_t*) *next;
    }

    // Done with search loop. Did we find our buffer?
    if (next == ptr) {
        // Found it! Set next-ptr of previous buffer to next-ptr
        // of this buffer to dequeue from list:
        *prevptr = *next;

        // Some accounting:
        PTBTEMPMEMDEC(next[1]);

        // Release:
        free(ptr);

        // Done.
        return;
    }

    // Oops.: Did not find matching buffer to pointer --> Trouble!
    printf("PTB-BUG: In PsychFreeTemp: Tried to free non-existent temporary membuffer %p!!! Ignored.\n", ptrbackup);
    fflush(NULL);
    return;
}

// Master cleanup routine: Frees all allocated memory:
void PsychFreeAllTempMemory(void)
{
    size_t* p = NULL;
    size_t* next = PsychTempMemHead;

    // Walk our whole buffer list and release all buffers on it:
    while (next != NULL) {
        // next points to current buffer to release. Make a copy of
        // next:
        p = next;

        // Update next to point to the next buffer to release:
        next = (size_t*) *p;

        // Some accounting:
        PTBTEMPMEMDEC(p[1]);

        // Release buffer p:
        free(p);

        // We're done with this buffer, next points to next one to release
        // or is NULL if all released...
    }

    // Done. NULL-out the list start ptr:
    PsychTempMemHead = NULL;

    // Sanity check:
    if (totalTempMemAllocated != 0) {
        // Cannot use PsychErrorXXX Routines here, because this is outside
        // the jumpbuffer context for our error-routines. Could lead to
        // infinite recursion!!!
        printf("PTB-CRITICAL BUG: Inconsistency detected in temporary memory allocator!\n");
        printf("PTB-CRITICAL BUG: totalTempMemAllocated = %li after PsychFreeAllTempMemory()!!!!\n", totalTempMemAllocated);
        fflush(NULL);

        // Reset to defined state.
        totalTempMemAllocated = 0;
    }

    return;
}

#endif
