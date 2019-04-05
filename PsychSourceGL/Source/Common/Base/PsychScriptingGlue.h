/*
  PsychToolbox3/Source/Common/PsychScriptingGlue.h

  AUTHORS:

  Allen.Ingling@nyu.edu         awi
  mario.kleiner.de@gmail.com    mk

  PLATFORMS: All.

  PROJECTS: All.

  HISTORY:

  08/12/02      awi     wrote it.

  DESCRIPTION:

    ScriptingGlue defines abstracted functions to pass values
    between the calling environment and the PsychToolbox.

  TO DO:

*/

//begin include once 
#ifndef PSYCH_IS_INCLUDED_ScriptingGlue
#define PSYCH_IS_INCLUDED_ScriptingGlue

#include "Psych.h"

#define kPsychUseDefaultArgPosition -1   //see ScreenArguments.h for more argument position defaults.  We should get rid of this.
#define kPsychNoArgReturn           -1

//typedefs
typedef double  *PsychFlagListType;

//NaN
double PsychGetNanValue(void);

// Simple function evaluation by scripting environment via feval() style functions:
int PsychRuntimeEvaluateString(const char* cmdstring);

// Put variable in native PsychGenericScriptType into workspace of runtime:
int PsychRuntimePutVariable(const char* workspace, const char* variable, PsychGenericScriptType* pcontent);

// Get read-only variable in native PsychGenericScriptType from workspace of runtime:
psych_bool PsychRuntimeGetVariablePtr(const char* workspace, const char* variable, PsychGenericScriptType** pcontent);

// Try to retrieve filesystem path to Psychtoolbox root or per user configfolder (the result from PsychtoolboxRoot() in Matlab/Octave) from runtime:
const char* PsychRuntimeGetPsychtoolboxRoot(psych_bool getConfigDir);

// Tell scripting glue to use/assume a C programming language memory layout for exchanging
// multi-dimensional (== 2D, 3D, n-D) matrices with the scripting environment if that layout
// promises higher efficiency and performance in data exchange. This is an opt-in, requesting
// C-layout if 'tryEnableCMemoryLayout' = TRUE, otherwise standard Fortran layout is assumed.
// The default is Fortran layout if this function does not get called, and it resets to Fortran
// layout at each return of control to the calling scripting environment. Iow. it is a per-
// module subfunction-call opt-in.
// The function returns TRUE if C memory layout is engaged, otherwise FALSE is returned.
// The caller may have to adjust its own data processing according to the returned value,
// unless the function is called with tryEnableCMemoryLayout = FALSE or not called at all, in
// which case Fortran layout is the thing.
psych_bool PsychUseCMemoryLayoutIfOptimal(psych_bool tryEnableCMemoryLayout);

//for memory pointers (void*):
psych_bool PsychCopyInPointerArg(int position, PsychArgRequirementType isRequired, void **ptr);
psych_bool PsychCopyOutPointerArg(int position, PsychArgRequirementType isRequired, void* ptr);

//for integers
psych_bool PsychCopyInIntegerArg(int position, PsychArgRequirementType isRequired, int *value);
psych_bool PsychCopyInIntegerArg64(int position,  PsychArgRequirementType isRequired, psych_int64 *value);
psych_bool PsychAllocInIntegerListArg(int position, PsychArgRequirementType isRequired, int *numElements, int **array);

//for float's aka singles:
psych_bool PsychAllocInFloatMatArg64(int position, PsychArgRequirementType isRequired, psych_int64 *m, psych_int64 *n, psych_int64 *p, float **array);
psych_bool PsychAllocInFloatMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, float **array);
psych_bool PsychAllocOutFloatMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, float **array);

//for doubles
psych_bool PsychCopyInDoubleArg(int position, PsychArgRequirementType isRequired, double *value);
psych_bool PsychAllocInDoubleArg(int position, PsychArgRequirementType isRequired, double **value);
psych_bool PsychAllocInDoubleMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, double **array);
psych_bool PsychAllocInDoubleMatArg64(int position, PsychArgRequirementType isRequired, psych_int64 *m, psych_int64 *n, psych_int64 *p, double **array);
psych_bool PsychCopyOutDoubleArg(int position, PsychArgRequirementType isRequired, double value);
psych_bool PsychAllocOutDoubleArg(int position, PsychArgRequirementType isRequired, double **value);
psych_bool PsychAllocOutDoubleMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, double **array);
psych_bool PsychCopyOutDoubleMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, double *fromArray);
// PsychAllocateNativeXXXMat() is for use with cell arrays and structs. The right way to do this is to use the normal function for returning
// values, detect if the position is -1, and if so accept the optional "nativeElement" value.
void PsychAllocateNativeDoubleMat(psych_int64 m, psych_int64 n, psych_int64 p, double **cArray, PsychGenericScriptType **nativeElement);
void PsychAllocateNativeUnsignedByteMat(psych_int64 m, psych_int64 n, psych_int64 p, psych_uint8 **cArray, PsychGenericScriptType **nativeElement);

// for unsigned 16 bit integer:
psych_bool PsychCopyOutUnsignedInt16MatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, psych_uint16 *fromArray);
psych_bool PsychAllocOutUnsignedInt16MatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, psych_uint16 **array);

//for psych_bool.  These should be consolidated with the flags below.
psych_bool PsychAllocOutBooleanMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, PsychNativeBooleanType **array);
psych_bool PsychCopyOutBooleanArg(int position, PsychArgRequirementType isRequired, PsychNativeBooleanType value);

//for flags.
psych_bool PsychCopyInFlagArg(int position, PsychArgRequirementType isRequired, psych_bool *argVal);
psych_bool PsychAllocInFlagArg(int position, PsychArgRequirementType isRequired, psych_bool **argVal);
psych_bool PsychAllocInFlagArgVector(int position,  PsychArgRequirementType isRequired, int *numElements, psych_bool **argVal);
psych_bool PsychCopyOutFlagArg(int position, PsychArgRequirementType isRequired, psych_bool argVal);
psych_bool PsychAllocOutFlagListArg(int position, PsychArgRequirementType isRequired, int numElements, PsychFlagListType *flagList); //flag lists are opaque.
void PsychLoadFlagListElement(int index, psych_bool value, PsychFlagListType flagList);
void PsychSetFlagListElement(int index, PsychFlagListType flagList);
void PsychClearFlagListElement(int index, PsychFlagListType flagList);

//for bytes
psych_bool PsychAllocInUnsignedByteMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, unsigned char **array);
psych_bool PsychAllocOutUnsignedByteMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, psych_uint8 **array);

//for strings
psych_bool PsychAllocInCharArg(int position, PsychArgRequirementType isRequired, char **str);
psych_bool PsychCopyOutCharArg(int position, PsychArgRequirementType isRequired, const char *str);

//query and govern argumuments.  Use these sparingly, usually you can let the "PsychAlloc*" and "PsychCopy*" functions above will do the work for you.  
int PsychGetNumInputArgs(void);
int PsychGetNumOutputArgs(void);
PsychError PsychCapNumInputArgs(int maxInputs);
PsychError PsychRequireNumInputArgs(int minInputs);
PsychError PsychCapNumOutputArgs(int maxNamedOutputs);
int PsychGetNumNamedOutputArgs(void);
psych_bool PsychIsArgPresent(PsychArgDirectionType direction, int position);
psych_bool PsychIsArgReallyPresent(PsychArgDirectionType direction, int position);
PsychArgFormatType PsychGetArgType(int position); //this is for inputs because outputs are unspecified
size_t PsychGetArgM(int position);
size_t PsychGetArgN(int position);
size_t PsychGetArgP(int position);
void PsychErrMsgTxt(char *s);
void PsychProcessErrorInScripting(PsychError error, const char* message);
void PsychEnableSubfunctions(void);
psych_bool PsychAreSubfunctionsEnabled(void);
psych_bool PsychCheckInputArgType(int position, PsychArgRequirementType isRequired, PsychArgFormatType argType);


// For the benefit of PsychStructGlue and PsychCellGlue. Don't use these unless you are writing more glue libraries. 
// They should probably be moved to a separate header file.
PsychError PsychSetReceivedArgDescriptor(int argNum, psych_bool allow64BitSizes, PsychArgDirectionType direction);
PsychError PsychSetSpecifiedArgDescriptor(  int position,
                                            PsychArgDirectionType direction,
                                            PsychArgFormatType type,
                                            PsychArgRequirementType isRequired,
                                            psych_int64	mDimMin,    // minimum minimum is 1   |
                                            psych_int64	mDimMax,    // minimum maximum is 1, maximum maximum is -1 meaning infinity
                                            psych_int64	nDimMin,    // minimum minimum is 1   |
                                            psych_int64	nDimMax,    // minimum maximum is 1, maximum maximum is -1 meaning infinity
                                            psych_int64	pDimMin,    // minimum minimum is 0
                                            psych_int64	pDimMax);   // minimum maximum is 0, maximum maximum is -1 meaning infinity
psych_bool PsychAcceptInputArgumentDecider(PsychArgRequirementType isRequired, PsychError matchError);
psych_bool PsychAcceptOutputArgumentDecider(PsychArgRequirementType isRequired, PsychError matchError);
PsychError PsychMatchDescriptors(void);

void PsychCheckSizeLimits(psych_int64 m, psych_int64 n, psych_int64 p);

// Low-level accessors to input/output parameters from the runtime:
const PsychGenericScriptType *PsychGetInArgPtr(int position);

#if PSYCH_LANGUAGE == PSYCH_MATLAB
mxArray **PsychGetOutArgMxPtr(int position);
const mxArray *PsychGetInArgMxPtr(int position);
#endif

#if PSYCH_LANGUAGE == PSYCH_PYTHON
int mxGetString(PyObject* arrayPtr, char* outstring, int outstringsize);
double mxGetScalar(const PyObject* arrayPtr);
PyObject* mxGetField(const PyObject* structArray, int index, const char* fieldName);
PyObject** PsychGetOutArgPyPtr(int position);
const PyObject *PsychGetInArgPyPtr(int position);
PyObject* PsychScriptingGluePythonDispatch(PyObject* self, PyObject* args);
const char* PsychGetPyModuleFilename(void);
#endif

//end include once
#endif
