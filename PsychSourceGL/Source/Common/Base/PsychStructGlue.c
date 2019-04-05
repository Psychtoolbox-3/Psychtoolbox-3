/*
    PsychToolbox2/Source/Common/PsychStructGlue.c

    AUTHORS:

    Allen.Ingling@nyu.edu       awi
    mario.kleiner.de@gmail.com  mk

    PLATFORMS: All

    PROJECTS: All

    HISTORY:
    12/31/02  awi   wrote it.
    03/28/11   mk   Make 64-bit clean.

    DESCRIPTION:

    PsychStructGlue defines abstracted functions to create structs passed
    between the calling environment and the PsychToolbox.

    TO DO:

    -All "PsychAllocOut*" functions should be modified to accept -1 as the position argument and
    in that case allocate an mxArray which may be embedded within a struct.  We then have only one
    function for settings structure array fields for all types.  This requires the addition of another
    pointer argument which points to the native struct representation, we could accept that optionally
    only when we specify -1 as the argument number.

    -PsychSetStructArray* functions should check that the named field is present in the struct
    and exit gracefully with an error.

    -Consider changing PsychSetStructArrayStructElement() do de-allocate the inner struct which it was
    passed.  Pass the pointer in indirectly and set it to point to the field embedded within the struct.
*/

#include "Psych.h"

#if PSYCH_LANGUAGE == PSYCH_MATLAB

// functions for outputting structs
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
    PsychAllocOutStructArray()

    -If argument is optional we allocate the structure even if the argument is not present.  If this bothers you,
    then check within the subfunction for the presense of a return argument before creating the struct array.  We
    allocate space regardeless of whether the argument is present because this is consistant with other "PsychAllocOut*"
    functions which behave this way because in some situations subfunctions might derive returned results from values
    stored in an optional argument.

    -If position is -1 then don't attempt to return the created structure to the calling environment.  Instead just
    allocate the structure and return it in pStruct.  This is how to create a structure which is embeded within another
    structure using PsychSetStructArrayStructArray().  Note that we use -1 as the flag and not NULL because NULL is 0 and
    0 is reserved for future use as a reference to the subfunction name, of if none then the function name.

    The special value numElements == -1 means to return a single struct, instead of a struct array.

*/
psych_bool PsychAllocOutStructArray(int position, 
                                    PsychArgRequirementType isRequired, 
                                    int numElements,
                                    int numFields, 
                                    const char **fieldNames,  
                                    PsychGenericScriptType **pStruct)
{
    mxArray **mxArrayOut;
    mwSize structArrayNumDims = 2;
    mwSize structArrayDims[2];
    PsychError matchError;
    psych_bool putOut;

    // -1 is special and means "a single struct", which is the same as a 1-element
    // struct array in Octave/Matlab. This -1 special case is needed for other scripting
    // environments that don't have mex api semantics:
    if (numElements == -1)
        numElements = 1;

    structArrayDims[0] = 1;
    structArrayDims[1] = numElements;

    if (position != kPsychNoArgReturn) {  //Return the result to both the C caller and the scripting environment.
        PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
        PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_structArray, isRequired, 1,1,numElements,numElements,0,0);
        *pStruct = mxCreateStructArray(structArrayNumDims, structArrayDims, numFields, fieldNames);
        matchError = PsychMatchDescriptors();
        putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
        if (putOut) {
            mxArrayOut = PsychGetOutArgMxPtr(position);
            *mxArrayOut = *pStruct;
        }
        return(putOut);
    } else{ //Return the result only to the C caller.  Ignore "required".    
        *pStruct = mxCreateStructArray(structArrayNumDims, structArrayDims, numFields, fieldNames);
        return(TRUE);
    }
}


/*  FONTS on OSX only.
    PsychAssignOutStructArray()
    Accept a pointer to a struct array and Assign the struct array to be the designated return variable.
*/
psych_bool PsychAssignOutStructArray( int position,
                                      PsychArgRequirementType isRequired,
                                      PsychGenericScriptType *pStruct)
{
    mxArray **mxArrayOut;
    PsychError matchError;
    psych_bool putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_structArray, isRequired, 1,1,0,kPsychUnboundedArraySize,0,0);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxArrayOut = PsychGetOutArgMxPtr(position);
        *mxArrayOut = pStruct;
    }
    return(putOut);
}

// functions for filling in struct elements by type 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
    PsychSetStructArrayStringElement()
    The variable "index", the index of the element within the struct array, is zero-indexed.
*/
void PsychSetStructArrayStringElement(const char *fieldName,
                                      int index,
                                      char *text,
                                      PsychGenericScriptType *pStruct)
{
    int fieldNumber;
    size_t numElements;
    psych_bool isStruct;
    mxArray *mxFieldValue;
    char errmsg[256];

    //check for bogus arguments
    numElements = mxGetM(pStruct) * mxGetN(pStruct);
    if ((size_t) index >= numElements)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a structure field at an out-of-bounds index");

    fieldNumber = mxGetFieldNumber(pStruct, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    isStruct = mxIsStruct(pStruct);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    //do stuff
    mxFieldValue = mxCreateString(text);
    mxSetField(pStruct, (mwIndex) index, fieldName, mxFieldValue);
    if (PSYCH_LANGUAGE == PSYCH_OCTAVE) mxDestroyArray(mxFieldValue);
}


/*
    PsychSetStructArrayDoubleElement()
    Note: The variable "index" is zero-indexed.
*/
void PsychSetStructArrayDoubleElement(const char *fieldName,
                                      int index,
                                      double value,
                                      PsychGenericScriptType *pStruct)
{
    int fieldNumber;
    size_t numElements;
    psych_bool isStruct;
    mxArray *mxFieldValue;
    char errmsg[256];

    //check for bogus arguments
    numElements = mxGetM(pStruct) * mxGetN(pStruct);
    if ((size_t) index >= numElements)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a structure field at an out-of-bounds index");

    fieldNumber = mxGetFieldNumber(pStruct, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    isStruct = mxIsStruct(pStruct);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    //do stuff
    mxFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxGetPr(mxFieldValue)[0] = value;
    mxSetField(pStruct, (mwIndex) index, fieldName, mxFieldValue);
    if (PSYCH_LANGUAGE == PSYCH_OCTAVE) mxDestroyArray(mxFieldValue);
}


/*
    PsychSetStructArrayBooleanElement()

    Note: The variable "index" is zero-indexed.
*/
void PsychSetStructArrayBooleanElement( const char *fieldName,
                                        int index,
                                        psych_bool state,
                                        PsychGenericScriptType *pStruct)
{
    int fieldNumber;
    size_t numElements;
    psych_bool isStruct;
    mxArray *mxFieldValue;
    char errmsg[256];

    //check for bogus arguments
    numElements = mxGetM(pStruct) * mxGetN(pStruct);
    if ((size_t) index >= numElements)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a structure field at an out-of-bounds index");

    fieldNumber = mxGetFieldNumber(pStruct, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    isStruct = mxIsStruct(pStruct);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    //do stuff
    mxFieldValue = mxCreateLogicalMatrix(1, 1);
    mxGetLogicals(mxFieldValue)[0] = state;
    mxSetField(pStruct, (mwIndex) index, fieldName, mxFieldValue);
    if (PSYCH_LANGUAGE == PSYCH_OCTAVE) mxDestroyArray(mxFieldValue);
}


/*
    PsychSetStructArrayStructElement()
*/
void PsychSetStructArrayStructElement(const char *fieldName,
                                      int index,
                                      PsychGenericScriptType *pStructInner,
                                      PsychGenericScriptType *pStructOuter)
{
    int fieldNumber;
    size_t numElements;
    psych_bool isStruct;
    char errmsg[256];

    //check for bogus arguments
    numElements = mxGetM(pStructOuter) * mxGetN(pStructOuter);
    if ((size_t) index >= numElements)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a structure field at an out-of-bounds index");

    fieldNumber = mxGetFieldNumber(pStructOuter, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    isStruct = mxIsStruct(pStructInner);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a struct field to a non-existent structure.");

    isStruct = mxIsStruct(pStructOuter);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    //do stuff
    mxSetField(pStructOuter, (mwIndex) index, fieldName, pStructInner); 
    if (PSYCH_LANGUAGE == PSYCH_OCTAVE) mxDestroyArray(pStructInner);
}


/*
    PsychSetStructArrayNativeElement()
 */
void PsychSetStructArrayNativeElement(const char *fieldName,
                                      int index,
                                      PsychGenericScriptType *pNativeElement,
                                      PsychGenericScriptType *pStructArray)
{
    int fieldNumber;
    size_t numElements;
    psych_bool isStruct;
    char errmsg[256];

    //check for bogus arguments
    numElements = mxGetM(pStructArray) * mxGetN(pStructArray);
    if ((size_t) index >= numElements)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a structure field at an out-of-bounds index");

    fieldNumber = mxGetFieldNumber(pStructArray, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    isStruct = mxIsStruct(pStructArray);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    //do stuff
    mxSetField(pStructArray, (mwIndex) index, fieldName, pNativeElement);
}

#endif
