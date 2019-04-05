/*
 *  PsychToolbox2/Source/Common/PsychCellGlue.c
 *
 *  Allen.Ingling@nyu.edu               awi
 *  mario.kleiner.de@gmail.com          mk
 *
 *  PLATFORMS: All
 *
 *  PROJECTS: All
 *
 *  HISTORY:
 *
 *    11/17/03  awi        wrote it.
 *    03/28/11   mk        Make 64-bit clean.
 *
 */

#include "Psych.h"

#if PSYCH_LANGUAGE == PSYCH_MATLAB

/*
 *    PsychAllocOutCellVector()
 *
 *    -If argument is optional we allocate the structure even if the argument is not present.  If this behavior bothers you,
 *    then check within your code for the presense of a return argument before creating the struct array.  We
 *    allocate space regardeless of whether the argument is present because this is consistant with other "PsychAllocOut*"
 *    functions which behave this way because in some situations subfunctions might derive returned results from values
 *    stored in an optional argument.
 *
 *    -If position is -1 then don't attempt to return the created structure to the calling environment.  Instead just
 *    allocate the structure and return it in pStruct.  This is how to create a structure which is embeded within another
 *    structure using PsychSetStructArrayStructArray().  Note that we use -1 as the flag and not NULL because NULL is 0 and
 *    0 is reserved for future use as a reference to the subfunction name, of if none then the function name.
 */
psych_bool PsychAllocOutCellVector(int position,
                                   PsychArgRequirementType isRequired,
                                   int numElements,
                                   PsychGenericScriptType **pCell)
{
    mxArray **mxArrayOut;
    mwSize cellArrayNumDims = 2;
    mwSize cellArrayDims[2];
    PsychError matchError;
    psych_bool putOut;

    cellArrayDims[0] = 1;
    cellArrayDims[1] = numElements;

    if (position != kPsychNoArgReturn) {  //Return the result to both the C caller and the scripting environment.
        PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
        PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_cellArray, isRequired, 1,1,numElements,numElements,0,0);
        *pCell = mxCreateCellArray(cellArrayNumDims, cellArrayDims);
        mxArrayOut = PsychGetOutArgMxPtr(position);
        matchError = PsychMatchDescriptors();
        putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
        if (putOut)
            *mxArrayOut=*pCell;
        return(putOut);
    } else { //Return the result only to the C caller, not to the calling environment.  Ignore "required".
        *pCell = mxCreateCellArray(cellArrayNumDims, cellArrayDims);
        return(TRUE);
    }
}


/*
 *    PsychSetCellVectorStringElement()
 *
 *    The variable "index", the index of the element within the struct array, is zero-indexed.
 */
void PsychSetCellVectorStringElement(int index,
                                     const char *text,
                                     PsychGenericScriptType *cellVector)
{
    size_t numElements;
    psych_bool isCell;
    mxArray *mxFieldValue;

    //check for bogus arguments
    numElements = mxGetM(cellVector) * mxGetN(cellVector);
    if ((size_t) index >= numElements)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a cell array field at an out-of-bounds index");

    isCell = mxIsCell(cellVector);
    if (!isCell)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a cell within a non-existent cell array.");

    //do stuff
    mxFieldValue = mxCreateString(text);
    mxSetCell(cellVector, (mwIndex) index, mxFieldValue);
    if (PSYCH_LANGUAGE == PSYCH_OCTAVE) mxDestroyArray(mxFieldValue);
}


/*
 *    PsychSetCellVectorDoubleElement() - UNUSED
 *
 *    Note: The variable "index" is zero-indexed.
void PsychSetCellVectorDoubleElement(int index,
                                     double value,
                                     PsychGenericScriptType *cellVector)
{
    size_t numElements;
    psych_bool isCell;
    mxArray *mxFieldValue;

    //check for bogus arguments
    numElements = mxGetM(cellVector) * mxGetN(cellVector);
    if ((size_t) index >= numElements)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a cell array field at an out-of-bounds index");

    isCell = mxIsCell(cellVector);
    if (!isCell)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent cell array.");

    //do stuff
    mxFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(mxFieldValue)[0]= value;
    mxSetCell(cellVector, (mwIndex) index, mxFieldValue);
    if (PSYCH_LANGUAGE == PSYCH_OCTAVE) mxDestroyArray(mxFieldValue);
}
 */


/*
 *    PsychSetCellVectorNativeElement() - UNUSED
void PsychSetCellVectorNativeElement(int index,
                                     PsychGenericScriptType *pNativeElement,
                                     PsychGenericScriptType *cellVector)
{
    size_t numElements;
    psych_bool isCell;

    //check for bogus arguments
    numElements = mxGetM(cellVector) * mxGetN(cellVector);
    if ((size_t) index >= numElements)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a cell array field at an out-of-bounds index");

    isCell = mxIsCell(cellVector);
    if (!isCell)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    //do stuff
    mxSetCell(cellVector, (mwIndex) index, pNativeElement);
}
*/


/*
 *    PsychAllocInNativeCellVector() - UNUSED
psych_bool PsychAllocInNativeCellVector(int position, PsychArgRequirementType isRequired, const PsychGenericScriptType **cellVector)
{
    PsychError matchError;
    psych_bool acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_cellArray, isRequired, 1, kPsychUnboundedArraySize,1, kPsychUnboundedArraySize,0, 1);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg)
        *cellVector = PsychGetInArgMxPtr(position);
    return(acceptArg);
}
*/


/*
 *    PsychAllocInNativeString() - UNUSED
psych_bool PsychAllocInNativeString(int position, PsychArgRequirementType isRequired, const PsychGenericScriptType **nativeString)
{
    PsychError matchError;
    psych_bool acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_char, isRequired, 1, kPsychUnboundedArraySize,1, kPsychUnboundedArraySize,0, 1);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg)
        *nativeString = PsychGetInArgMxPtr(position);
    return(acceptArg);
}
*/

#endif
