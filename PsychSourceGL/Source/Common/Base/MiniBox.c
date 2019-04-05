/*
 * MiniBox.cpp
 *
 * DESCRIPTION:
 *
 *  The minimal set of functions necessary to support SCREEN style error messages, help and
 *  argument list parsing. Originally copied from David Brainard and Denis Pelli's
 *  Psychtoolbox project, with some modification.
 *  It has since evolved into a miscellaneous category.
 *
 * HISTORY:
 *
 *    08/22/01    awi       Windows
 *    11/15/01    awi       Mac OS 9
 *    08/20/01    awi       Created file.
 *    08/20/01    awi       Replaced fancy memory allocation in PrintfExit with static array.
 *    11/26/01    awi       Replaced references to ActiveWireTable in PsychMatch with ProjectTable
 *                          This is so we can share the MiniBox file between ActiveWire and
 *                          Joystick projects.
 *    12/17/01    awi       Commented out include "joystick.h" and included stdafx.h for GetChar project.
 *    12/17/01    awi       All of the grief with ordering the header inclusion and porbability was
 *                          caused by PsychMatch using the project table to read the case-sensitive
 *                          character comparision flag from the psych table.  That's a pretty useless
 *                          feature anyway, so I just commented that out to make GetChar compile.
 *    1/20/04     awi       Cosmetic.
 *    3/19/10      mk       Cosmetic and make 64-bit clean.
 *
 */

#include "Psych.h"
#include <ctype.h>

// The handle of the masterthread - The Matlab/Octave/Python/PTB main interpreter thread:
static psych_threadid masterthread = (psych_threadid) 0;

#define MAX_PRINTFEXIT_LEN  2000

static psych_bool isPsychMatchCaseSensitive = FALSE;

/* PrintfExit used some fancy stuff to allocate the memory
 * which holds the the error message string passed to
 * mexErrMsgTxt.  I just used a static array here. -awi
 */
int PrintfExit(const char *format,...)
{
    char s[MAX_PRINTFEXIT_LEN];
    va_list args;
    int i;

    va_start(args,format);
    i = vsnprintf(s, MAX_PRINTFEXIT_LEN - 1, format, args);
    va_end(args);
    if(i + 1 > MAX_PRINTFEXIT_LEN)
        printf("PTB-ERROR: Buffer overrun in PrintfExit\n");
    PsychErrMsgTxt(s);

    return 0;
}

char *BreakLines(char *string,long lineLength)
{
    long i,leftMargin,rightMargin,length;
    int here;

    leftMargin=0;
    length = (long) strlen(string);
    while(1){
        rightMargin=leftMargin+lineLength;
        if(rightMargin>=length)return string;        /* successful completion */
            here=0;
        if(!here)for(i=leftMargin;i<rightMargin;i++)if(string[i]=='\n'){
            here=1;
            break;
        }
        if(!here)for(;i>=leftMargin;i--)if(string[i]==' ' || string[i]=='\n'){
            here=1;
            break;
        }
        if(!here)for(i=leftMargin;i<length;i++)if(string[i]==' ' || string[i]=='\n'){
            here=1;
            break;
        }
        if(!here)return string;
        string[i]='\n';
        leftMargin=i+1;
    }
}

psych_bool PsychIsPsychMatchCaseSensitive(void)
{
    return(isPsychMatchCaseSensitive);
}

void PsychSetPsychMatchCaseSenstive(psych_bool arg)
{
    isPsychMatchCaseSensitive=arg;
}

// Compare two strings for equality. Ignore case if Psychtoolbox preferences ignore case is set.
psych_bool PsychMatch(char *s1,char *s2)
{
    char a;

    if(!isPsychMatchCaseSensitive){
        do{
            a=*s1++;
            if(tolower(a)!=tolower(*s2++))return 0;
        }while(a!='\0');
        return 1;
    }else
        return strcmp(s1,s2)==0;
}

char *int2str(psych_int64 num)
{
    static char numStr[256];
    #if PSYCH_SYSTEM != PSYCH_WINDOWS
    sprintf(numStr, "%lld", (long long int) num);
    #else
    // TODO FIXME AUDIT 64BIT : Figure out a way to handle psych_int64 printing on Windows:
    sprintf(numStr, "%d", (int) num);
    #endif
    return(numStr);
}

size_t PsychIndexElementFrom2DArray(size_t mDim/*|Y|*/, size_t nDim/*|X|*/, size_t m/*y*/, size_t n/*x*/)
{
    (void) nDim;
    return(n*mDim + m);
}

size_t PsychIndexElementFrom3DArray(size_t mDim/*|Y|*/, size_t nDim/*|X|*/, size_t pDim/*|Z|*/, size_t m/*y*/, size_t n/*x*/, size_t p/*z*/)
{
    (void) pDim;
    return(p*mDim*nDim + n*mDim + m);  //planeindex * planesize + columnindex * columsize + rowindex
}

size_t PsychIndexPlaneFrom3DArray(size_t mDim, size_t nDim, size_t pDim, size_t planeIndex)
{
    (void) pDim;
    return(planeIndex*mDim*nDim);
}

psych_int64 maxInt(psych_int64 a, psych_int64 b)
{
    if(a>b)
        return(a);
    return(b);
}

/*
 *    PsychIsIntegerInDouble(double *value)
 *
 *    If the value stored in the specified double does not have a fractional part an the value is within
 *    the bounds of the signed/unsigned integer type then return TRUE.
 *    We allow to store 32-Bit unsigned int values inside 32-Bit signed int's for this validation,
 *    so uint32's can be passed, e.g., via PsychCopyInIntegerArg(). They would wrap to negative in
 *    the returned int32 if they exceed +INT_MAX, but allowing to cast forward and backward between
 *    uint32 and int32 has value for access to some hardware api's, e.g., forum message #17256.
 */
psych_bool PsychIsIntegerInDouble(double *value)
{
    return((*value >= INT_MIN) && (*value <= (double) 0xffffffff) && (floor(*value) == *value));
}

/* Check if it is a 64 bit integer (psych_int64) packed into a double:
 * This check will already fail for any integer greater than about 2^52
 * as double can't represent them accurately anymore.
 */
psych_bool PsychIsInteger64InDouble(double *value)
{
    return((*value >= -9.22337203685478e+18) && (*value <= 9.22337203685478e+18) && (floor(*value) == *value));
}

void PsychInitMasterThreadId(void)
{
    // Assign unique id of this thread (the Matlab/Octave/Python main interpreter thread)
    // as masterthread. This masterthread is the only one allowed to execute certain code
    // that may only be safe for execution on the main thread, e.g., stuff that interacts
    // with the interpreter / host runtime api.
    masterthread = PsychGetThreadId();
}

psych_threadid PsychGetMasterThreadId(void)
{
    return(masterthread);
}

psych_bool PsychIsMasterThread(void)
{
    return(PsychIsCurrentThreadEqualToId(masterthread));
}
