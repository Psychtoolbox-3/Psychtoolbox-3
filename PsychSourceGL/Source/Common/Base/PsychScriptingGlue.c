/*
 * PsychSourceGL/Source/Common/Base/PsychScriptingGlue.c
 *
 * AUTHORS:
 *
 * Allen.Ingling@nyu.edu        awi
 * mario.kleiner.de@gmail.com   mk
 *
 * PLATFORMS: All -- Glue layer to the runtime environment. Runtime environment independent parts.
 *
 * HISTORY:
 * 05/07/06   mk    Derived from Allen's PsychScriptingGlue
 * 11/08/06   mk    Fixes for Matlab beta on IntelMac -- Replace mxGetPr() by mxGetData()
 *                  or mxGetScalar() in places where this is appropriate. Using mxGetPr()
 *                  in the debug-build of the Matlab beta triggers an assertion when
 *                  passing a non-double array to mxGetPr().
 *
 * DESCRIPTION:
 *
 * ScriptingGlue defines abstracted functions to pass values
 * between the calling scripting environment and the PsychToolbox.
 *
 * NOTES:
 *
 * About default arguments:  In previous versions of the Psychtoolbox any matrix of size m*n=0
 * stood for the the "default" matrix.  When passed as an argument, it indicated that the
 * default value for that argument should be used.  This is useful when "omitting" intervening
 * arguments.
 *
 * Because each SCREEN subfunction interpreted arguments independently the ambiguities which
 * are discussed below did not have to be addressed but in the subfunctions which they arrose.
 * (which might be none).  The introduction of abstracted functions in ScriptingGlue mandates
 * a uniform policy for resloving ambiguities.
 *
 * Sometimes we want to pass an argument of size 0x0 and mean argument of size 0x0, not the
 * default matrix.  So ScriptingGlue functions which retrieve input arguments can not safetly
 * interpret an empty matrix as the default matrix.
 *
 * The problem is not as bad as it seems, because we can pass an empty
 * numerical matrix, "[]" when a string argument is expected, or pass an empty string "''" when
 * a numerical argument is expected.  Only in the case when an argument may be either a string or a number,
 * and 0 size arguments of both types are meaningful do we have a problem.  The case does not seem likely ever
 * to arise.
 *
 * For users, having two default arguments, '' and [],  and having to decide which to use depending on the
 * type of argument accepted, complicates the use of default arguments unpleasantly.  Furthermore, empty strings
 * are meaninful as strings, but empty numerical matrices are rarely meaninful as matrices. (why is that?)
 * Therefore, the best policy for ScriptingGlue functions would be: ScriptingGlue  functions which
 * retrieve string arguments will only interpret [] as the default matrix and will interpret '' as
 * the empty string.  ScriptingGlue functions which retrieve numerical arguments will accept either
 * [] or '' to be the empty string.
 *
 * So [] when passed for a number is always interpreted as the default matrix,
 * [] is the only value which stands for default when passed for a string,  Therefore, we can
 * reduce this further and accept only [] to stand for default, simplifing the users's decision of when to
 * use '' and when to use [], by ALWAYS using [].
 *
 * So in conclusion:
 * -[] and only [] means the default matrix.
 * -If you want a user to pass [] to mean a 0x0 matrix, too bad, you can't do that.
 * All ScriptingGlue functions will report that the argument was not present if the user
 * passes [].
 *
 */

// During inclusion of Psych.h, we define the special flag PTBINSCRIPTINGGLUE. This
// will cause some of the system headers in Psych.h not to be included during build
// of PsychScriptingGlue.c:
#define PTBINSCRIPTINGGLUE 1
#include "Psych.h"
#undef PTBINSCRIPTINGGLUE

////Static functions local to ScriptingGlue.c.
// _____________________________________________________________________________________

void InitializeSynopsis(char *synopsis[],int maxStrings);

#define MAX_SYNOPSIS 100
#define MAX_CMD_NAME_LENGTH 100

//Static variables local to ScriptingGlue.c.  The convention is to append a abbreviation in all
//caps of the C file name to the variable name.

static psych_bool subfunctionsEnabledGLUE = FALSE;

// Forward declaration for GNU/Linux compile:
#ifdef __cplusplus
    extern "C" void ScreenCloseAllWindows(void);
#else
    void ScreenCloseAllWindows(void);
#endif

/*
 *  functions for enabling and testing subfunction mode
 */
void PsychEnableSubfunctions(void)
{
    subfunctionsEnabledGLUE = TRUE;
}


psych_bool PsychAreSubfunctionsEnabled(void)
{
    return(subfunctionsEnabledGLUE);
}


PsychError PsychSetSpecifiedArgDescriptor(int                       position,
                                          PsychArgDirectionType     direction,
                                          PsychArgFormatType        type,
                                          PsychArgRequirementType   isRequired,
                                          psych_int64               mDimMin,        // minimum minimum is 1   |
                                          psych_int64               mDimMax,        // minimum maximum is 1, maximum maximum is -1 meaning infinity
                                          psych_int64               nDimMin,        // minimum minimum is 1   |
                                          psych_int64               nDimMax,        // minimum maximum is 1, maximum maximum is -1 meaning infinity
                                          psych_int64               pDimMin,        // minimum minimum is 0
                                          psych_int64               pDimMax)        // minimum maximum is 0, maximum maximum is -1 meaning infinity
{
    PsychArgDescriptorType d;

    // Check size of output dimensions if this is an output operation:
    if (direction == PsychArgOut) {
        // Do not exceed index size limits of hw/os/build architecture,
        // be it 32 bit or 64 bit:
        if (((size_t) mDimMin > SIZE_MAX) || ((size_t) mDimMax > SIZE_MAX) ||
            ((size_t) nDimMin > SIZE_MAX) || ((size_t) nDimMax > SIZE_MAX) ||
            ((size_t) pDimMin > SIZE_MAX) || ((size_t) pDimMax > SIZE_MAX)) {

            printf("PTB-ERROR: Tried to return a vector or matrix whose size along at least one dimension\n");
            printf("PTB-ERROR: exceeds the maximum supported number of elements.\n");

            if (sizeof(size_t) == 4) {
                printf("PTB-ERROR: This is a limitation of all 32 bit versions of Psychtoolbox.\n");
                printf("PTB-ERROR: You'd need to use a Psychtoolbox for 64-bit Matlab or 64-bit Octave\n");
                printf("PTB-ERROR: on a 64-bit operating system to get rid of this limit.\n");
            }

            PsychErrorExitMsg(PsychError_user, "One of the dimensions of a returned matrix or vector exceeds maximum number of elements. This is not supported on your setup!");
        }

        // Limits ok for given hw/os/build architecture. Check if they're ok for the
        // scripting environment as well:
        PsychCheckSizeLimits((size_t) mDimMin, (size_t) nDimMin, (size_t) pDimMin);
        PsychCheckSizeLimits((size_t) mDimMax, (size_t) nDimMax, (size_t) pDimMax);
    }

    d.position = position;
    d.direction = direction;
    d.type = type;
    // d.isThere                //field set only in the received are descriptor, not in the specified argument descriptor
    d.isRequired = isRequired;  //field set only in the specified arg descritor, not in the received argument descriptot.
    d.mDimMin = mDimMin;
    d.mDimMax = mDimMax;
    d.nDimMin = nDimMin;
    d.nDimMax = nDimMax;
    d.pDimMin = pDimMin;
    d.pDimMax = pDimMax;

    //NOTE that we are not setting the d.numDims field because that is inferred from pDimMin and pDimMax and the 3 dim cap.
    PsychStoreArgDescriptor(&d,NULL);
    return(PsychError_none);
}


/*
 *    PsychAcceptInputArgumentDecider()
 *
 *    This is a subroutine of Psychtoolbox functions such as PsychCopyInDoubleArg() which read in arguments to Psychtoolbox functino
 *    passed from the scripting environment.
 *
 *    Accept one constant specifying whether an argument is either required, optional, or anything will be allowed and another constant
 *    specifying how the provided argument agrees with the specified argument.  Based on the relationship between those constants either:
 *
 *        * Return TRUE indicating that the caller should read in the argument and itself return TRUE to indicate that the argument has been read.
 *        * Return FALSE indicating that the caller should ignore the argument and itself return FALSE to indicate that the argument was not read.
 *        * Exit to the calling environment with an error to indicate that the provided argument did not match the requested argument and that
 *          it was required to match.
 *
 *    The domain of supplied arguments is:
 *
 *    matchError:
 *        PsychError_internal                  -Internal Psychtoolbox error
 *        PsychError_invalidArg_absent         -There was no argument provided
 *        PsychError_invalidArg_type           -The argument was present but not the specified type
 *        PsychError_invalidArg_size           -The argument was presnet and the specified type but not the specified size
 *        PsychError_none                      -The argument matched the specified argument
 *
 *    isRequired:
 *        kPsychArgRequired                    - the argument must be present and must match the specified descriptor
 *        kPsychArgOptional                    - the argument must either be absent or must be present and match the specified descriptor
 *        kPsychArgAnything                    - the argument can be absent or anything
 *
 */
psych_bool PsychAcceptInputArgumentDecider(PsychArgRequirementType isRequired, PsychError matchError)
{
    if (isRequired==kPsychArgRequired) {
        if (matchError)
            PsychErrorExit(matchError);
        else
            return(TRUE);
    } else if (isRequired==kPsychArgOptional) {
        if (matchError==PsychError_invalidArg_absent)
            return(FALSE);
        else if (matchError)
            PsychErrorExit(matchError);
        else
            return(TRUE);
    } else if (isRequired==kPsychArgAnything) {
        if (!matchError)
            return(TRUE);
        else if (matchError==PsychError_invalidArg_absent)
            return(FALSE);
        else if (matchError==PsychError_invalidArg_type)
            return(FALSE);
        else if (matchError==PsychError_invalidArg_size)
            return(FALSE);
        else
            PsychErrorExit(matchError);
    }
    PsychErrorExitMsg(PsychError_internal, "Reached end of function unexpectedly");
    return(FALSE);            //make the compiler happy
}


/*
 *
 *    PsychAcceptOutputArgumentDecider()
 *
 *    This is a subroutine of Psychtoolbox functions such as PsychCopyCopyDoubleArg() which output arguments from Psychtoolbox functions
 *    back to the scripting environment.
 *
 */
psych_bool PsychAcceptOutputArgumentDecider(PsychArgRequirementType isRequired, PsychError matchError)
{
    if (isRequired==kPsychArgRequired) {
        if (matchError)
            PsychErrorExit(matchError);                     //the argument was required and absent so exit with an error. Or there was some other error.
        else
            return(TRUE);                                   //the argument was required and present so go read it.
    } else if (isRequired==kPsychArgOptional) {
        if (!matchError)
            return(TRUE);                                   //the argument was optional and present so go read it.
        else if (matchError==PsychError_invalidArg_absent)
            return(FALSE);                                  //the argument as optional and absent so dont' read  it.
        else if (matchError)
            PsychErrorExit(matchError);                     //there was some other error
    } else if (isRequired==kPsychArgAnything)
        PsychErrorExitMsg(PsychError_internal, "kPsychArgAnything argument passed to an output function.  Use kPsychArgOptional");
    else
        PsychErrorExit(PsychError_internal);

    PsychErrorExitMsg(PsychError_internal, "End of function reached unexpectedly");
    return(FALSE);        //make the compiler happy
}


/*
 *    PsychMatchDescriptors()
 *
 *    Compare descriptors for specified and received arguments. Return a mismatch error if they are
 *    incompatible, otherwise return a no error.
 *
 *    PsychMatchDescriptors compares:
 *        The argument type
 *        The argument size
 *        Argument presense
 *
 *    PsychMatchDescripts can return any of the following values describing the relationship between an
 *    argument provided from the scripting environment and argument requested by a Psychtoolbox module:
 *        PsychError_internal               -Internal Psychtoolbox error
 *        PsychError_invalidArg_absent      -There was no argument provided
 *        PsychError_invalidArg_type        -The argument was present but not the specified type
 *        PsychError_invalidArg_size        -The argument was presnet and the specified type but not the specified size
 *        PsychError_none                   -The argument matched the specified argument
 *
 *    This function should be enhnaced to report the nature of the disagrement
 */
PsychError PsychMatchDescriptors(void)
{
    PsychArgDescriptorType *specified, *received;

    PsychGetArgDescriptor(&specified, &received);

    //check for various bogus conditions resulting only from Psychtoolbox bugs and issue an internal error.
    if (specified->position != received->position)
        PsychErrorExit(PsychError_internal);
    if (specified->direction != received->direction)
        PsychErrorExit(PsychError_internal);

    if (specified->direction==PsychArgOut) {
        if (received->isThere==kPsychArgPresent || received->isThere==kPsychArgFixed)
            return(PsychError_none);
        else
            return(PsychError_invalidArg_absent);
    }
    if (specified->direction==PsychArgIn) {
        if (received->isThere==kPsychArgAbsent)
            return(PsychError_invalidArg_absent);
        //otherwise the argument is present and we proceed to the argument type and size checking block below
    }

    //if we get to here it means that an input argument was supplied.  Check if it agrees in type and size with the specified arg and return
    // an error type accordingly
    if (!(specified->type & received->type))
        return(PsychError_invalidArg_type);
    if (received->mDimMin != received->mDimMax || received->nDimMin != received->nDimMax ||  received->pDimMin != received->pDimMax)
        PsychErrorExit(PsychError_internal);    //unnecessary mandate
        if (received->mDimMin < specified->mDimMin)
            return(PsychError_invalidArg_size);
        if (received->nDimMin < specified->nDimMin)
            return(PsychError_invalidArg_size);
        if (specified->pDimMin != kPsychUnusedArrayDimension && received->pDimMin < specified->pDimMin)
            return(PsychError_invalidArg_size);
        if (specified->mDimMax != kPsychUnboundedArraySize && received->mDimMax > specified->mDimMax)
            return(PsychError_invalidArg_size);
        if (specified->nDimMax != kPsychUnboundedArraySize && received->nDimMax > specified->nDimMax)
            return(PsychError_invalidArg_size);
        if (specified->pDimMax != kPsychUnusedArrayDimension && specified->pDimMax != kPsychUnboundedArraySize && received->pDimMax > specified->pDimMax)
            return(PsychError_invalidArg_size);
        if (received->numDims > 3)  //we don't allow matrices with more than 3 dimensions.
            return(PsychError_invalidArg_size);

        //if we get to here it means that  the block above it means
        return(PsychError_none);
}


//functions for project access to module call arguments
//___________________________________________________________________________________________


//functions which query the number and nature of supplied arguments
PsychError PsychCapNumInputArgs(int maxInputs)
{
    if (PsychGetNumInputArgs() > maxInputs)
        return(PsychError_extraInputArg);
    else
        return(PsychError_none);
}


PsychError PsychRequireNumInputArgs(int minInputs)
{
    if (PsychGetNumInputArgs() < minInputs)
        return(PsychError_missingInputArg);
    else
        return(PsychError_none);
}


/*
 *    The argument is present if anything was supplied, including the default matrix
 */
psych_bool PsychIsArgReallyPresent(PsychArgDirectionType direction, int position)
{
    return(direction == PsychArgOut ? PsychGetNumOutputArgs() >= position : PsychGetNumInputArgs() >= position);
}


/*
 *    PsychCheckInputArgType()
 *
 *    Check that the input argument at the specifid position matches at least one of the types passed in the argType
 *    argument. If the argument violates the proscription exit with an error.  Otherwise return a psych_bool indicating
 *    whether the argument was present.
 */
psych_bool PsychCheckInputArgType(int position, PsychArgRequirementType isRequired, PsychArgFormatType argType)
{
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, argType, isRequired, 0,kPsychUnboundedArraySize,0,kPsychUnboundedArraySize,0,kPsychUnboundedArraySize);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    return(acceptArg);
}


/*functions with input arguments.
 * ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 * ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 */

/*
 * A)input argument mandatory:
 *
 *    1)input argument not present:         exit with error.
 *    2)input argument present:             set *array to the input matrix, *m, *n, and *p to its dimensions, return TRUE.
 * B)input argument optional:
 *
 *    1)input argument not present:         return FALSE
 *    2)input argument present:             set *array to the input matrix, *m, *n, and *p to its dimensions, return TRUE.
 *
 */

/*
 *    PsychAllocOutFlagListArg()
 *
 *    This seems silly.  Find out where its used and consider using an array of booleans instead.  Probably the best thing
 *    is just to transparently map arrays of booleans to logical arrays MATLAB.
 *
 *    In Matlab our psych_bool flags are actually doubles.  This will not be so in all scripting languages.  We disguise the
 *    implementation of psych_bool flags within the scripting envrironment by making the flag list opaque and
 *    providing accessor fucntions PsychLoadFlagListElement, PsychSetFlagListElement, and PsychClearFlagListElement.
 *
 *    TO DO: maybe this should return a logical array instead of a bunch of doubles.  Itwould be better for modern versions
 *    of MATLAB which store doubles as bytes internally.
 *
 *
 */
psych_bool PsychAllocOutFlagListArg(int position, PsychArgRequirementType isRequired, int numElements, PsychFlagListType *flagList)
{
    return(PsychAllocOutDoubleMatArg(position, isRequired, (int)1, numElements, (int)0, flagList));
}


void PsychLoadFlagListElement(int index, psych_bool value, PsychFlagListType flagList)
{
    flagList[index]=(double)value;
}


void PsychSetFlagListElement(int index, PsychFlagListType flagList)
{
    flagList[index]=(double)1;
}


void PsychClearFlagListElement(int index, PsychFlagListType flagList)
{
    flagList[index]=(double)0;
}
