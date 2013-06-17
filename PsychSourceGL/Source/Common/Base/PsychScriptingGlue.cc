
/*
  PsychToolbox2/Source/Common/PsychOctaveScriptingGlue.cc		
  
  AUTHORS:
  Allen.Ingling@nyu.edu		     awi 
  mario kleiner at tuebingen.mpg.de  mk
 
  PLATFORMS: All -- Glue layer for GNU/Octave and Matlab runtime environment.
  
  HISTORY:
  05/07/06 	mk	Derived from Allen's PsychScriptingGlue
  11/08/06      mk      Fixes for Matlab beta on IntelMac -- Replace mxGetPr() by mxGetData()
                        or mxGetScalar() in places where this is appropriate. Using mxGetPr()
			in the debug-build of the Matlab beta triggers an assertion when
			passing a non-double array to mxGetPr().

  DESCRIPTION:
  
	ScriptingGlue defines abstracted functions to pass values 
	between the calling Matlab- or GNU/Octave environment and the PsychToolbox. 
  
  NOTES:
  
  About default arguments:  In previous versions of the Psychtoolbox any matrix of size m*n=0 
  stood for the the "default" matrix.  When passed as an argument, it indicated that the
  default value for that argument should be used.  This is useful when "omitting" intervening
  arguments. 
  
  Because each SCREEN subfunction interpreted arguments independently the ambiguities which 
  are discussed below did not have to be addressed but in the subfunctions which they arrose. 
  (which might be none).  The introduction of abstracted functions in ScriptingGlue mandates 
  a uniform policy for resloving ambiguities.   
  
  Sometimes we want to pass an argument of size 0x0 and mean argument of size 0x0, not the 
  default matrix.  So ScriptingGlue functions which retrieve input arguments can not safetly 
  interpret an empty matrix as the default matrix.
  
  The problem is not as bad as it seems, because we can pass an empty 
  numerical matrix, "[]" when a string argument is expected, or pass an empty string "''" when 
  a numerical argument is expected.  Only in the case when an argument may be either a string or a number,
  and 0 size arguments of both types are meaningful do we have a problem.  The case does not seem likely ever
  to arise. 
  
  For users, having two default arguments, '' and [],  and having to decide which to use depending on the 
  type of argument accepted, complicates the use of default arguments unpleasantly.  Furthermore, empty strings
  are meaninful as strings, but empty numerical matrices are rarely meaninful as matrices. (why is that?)
  Therefore, the best policy for ScriptingGlue functions would be: ScriptingGlue  functions which 
  retrieve string arguments will only interpret [] as the default matrix and will interpret '' as
  the empty string.  ScriptingGlue functions which retrieve numerical arguments will accept either
  [] or '' to be the empty string.  
  
  So [] when passed for a number is always interpreted as the default matrix,  
  [] is the only value which stands for default when passed for a string,  Therefore, we can 
  reduce this further and accept only [] to stand for default, simplifing the users's decision of when to
  use '' and when to use [], by ALWAYS using [].  
  
  So in conclusion:
   -[] and only [] means the default matrix. 
   -If you want a user to pass [] to mean a 0x0 matrix, too bad, you can't do that.  
   All ScriptingGlue functions will report that the argument was not present if the user
   passes [].    
  
       
  
  TO DO: 
    
    - baseFunctionInvoked and PsychSubfunctionEnabled are redundent, keep only baseFunctionInvoked
  	
        
    Less Important:
    
        -Expand for use with N dimensional arrays.  
  	The final required argument to these functions should be the number of dimensions and 
  	thereafter optional function arguments can give the size of each dimension.
  	
  	

*/

// We have to use #ifdef PTBOCTAVE here, because PSYCH_LANGUAGE is not yet defined.
// Definition will happen a few lines downwards in Psych.h, but we can't reorder code
// here :(
#ifdef PTBOCTAVE

// I dont know why, but it is *absolutely critical* that octave/oct.h is included
// before *any* other header file, esp. Psych.h, otherwise the C++ compiler f%%2!s up
// completely!
#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/ov-struct.h>
#include <octave/ov-cell.h>

#endif 

// During inclusion of Psych.h, we define the special flag PTBINSCRIPTINGGLUE. This
// will cause some of the system headers in Psych.h not to be included during build
// of PsychScriptingGlue.cc - They are not needed for this and they conflict with
// GNU/Octave header files!
#define PTBINSCRIPTINGGLUE 1
    #include "Psych.h"
#undef PTBINSCRIPTINGGLUE

// Define this to 1 if you want lots of debug-output for the Octave-Scripting glue.
#define DEBUG_PTBOCTAVEGLUE 0

// Special hacks to allow Psychtoolbox to build for GNU-OCTAVE:
#if PSYCH_LANGUAGE == PSYCH_OCTAVE

#include <string.h>
#include <setjmp.h>

// This jump-buffer stores CPU- and stackstate at the position
// where our octFunction() dispatcher actually starts executing
// the functions or subfunctions code. That is the point where
// a PsychErrorExit() or PsychErrorExitMsg() will return control...
jmp_buf jmpbuffer;

// Error exit handler: Replacement for Matlabs MEX-handler:
// Prints the error-string with Octaves error printing facilities,
// sets Octave error state and longjmp's to the cleanup routine
// at the end of our octFunction dispatcher...
void mexErrMsgTxt(const char* s) {
  if (s && strlen(s)>0) {
    error("%s: %s", mexFunctionName,s);
  }
  else error("");

  // Use the jump-buffer to unwind the stack...
  longjmp(jmpbuffer, 1);
}

// Interface to Octave's printf...
void mexPrintf(const char* fmt, ...)
{
  va_list args;
  va_start (args, fmt);
  octave_vformat(octave_diary, fmt, args);
  octave_vformat(octave_stdout, fmt, args);
  va_end(args);
}

void* mxMalloc(int size)
{
  return(PsychMallocTemp((unsigned long) size));
}

void* mxCalloc(int size, int numelements)
{
  return(PsychCallocTemp((unsigned long) size, (unsigned long) numelements));
}

void mxFree(void* p)
{
  PsychFreeTemp(p);
}

int mexCallMATLAB(const int nargout, mxArray* argout[], 
		  const int nargin, const mxArray* argin[],
		  const char* fname)
{
  PsychErrorExitMsg(PsychError_unimplemented, "FATAL Error: Internal call to mexCallMATLAB(), which is not yet implemented on GNU/Octave port!");
}

double mxGetNaN(void)
{
  return(lo_ieee_nan_value());
}

mxArray* mxCreateNumericArray(int numDims, int dimArray[], int arraytype, int realorcomplex)
{
  mxArray* retval;
  int rows, cols, layers;
  if (numDims>3) PsychErrorExitMsg(PsychError_internal, "FATAL Error: mxCreateNumericArray: Tried to create matrix with more than 3 dimensions!");
  rows = dimArray[0];
  cols = (numDims>1) ? dimArray[1] : 1;
  layers = (numDims>2) ? dimArray[2] : 1;
  dim_vector mydims((numDims>2) ? dim_vector(rows, cols, layers) : dim_vector(rows, cols));

  // Allocate our mxArray-Struct:
  retval = (mxArray*) PsychMallocTemp(sizeof(mxArray));

  // Create corresponding octave_value object for requested type and size of matrix.
  // Retrieve raw pointer to contained data and store it in our mxArray struct as well...

  if (arraytype==mxUINT8_CLASS) {
    if (DEBUG_PTBOCTAVEGLUE) printf("NEW UINT8 MATRIX: %i,%i,%i\n", rows, cols, layers); fflush(NULL);
    // Create empty uint8ND-Array of type mxREAL...
    uint8NDArray m(mydims);
    // Retrieve a pointer to internal representation. As m is new
    // this won't trigger a deep-copy.
    retval->d = (void*) m.data();
    if (DEBUG_PTBOCTAVEGLUE) printf("M-DATA %p\n", retval->d); fflush(NULL);
    // Build a new oct_value object from Matrix m: This is a
    // shallow-copy.
    octave_value* ovp = new octave_value();
    *ovp = m;
    retval->o = (void*) ovp;
    // At this point we can safely destroy Matrix m, as the new
    // octave_object holds a reference to its representation.    
  }
  else if (arraytype==mxDOUBLE_CLASS && rows*cols*layers > 1) {
    if (DEBUG_PTBOCTAVEGLUE) printf("NEW DOUBLE MATRIX: %i,%i,%i\n", rows, cols, layers); fflush(NULL);
    // Create empty ND-Array of type mxREAL...
    NDArray m(mydims);
    // Retrieve a pointer to internal representation. As m is new
    // this won't trigger a deep-copy.
    retval->d = (void*) m.data();
    if (DEBUG_PTBOCTAVEGLUE) printf("M-DATA %p\n", retval->d); fflush(NULL);
    // Build a new oct_value object from Matrix m: This is a
    // shallow-copy.
    octave_value* ovp = new octave_value();
    *ovp = m;
    retval->o = (void*) ovp;
    // At this point we can safely destroy Matrix m, as the new
    // octave_object holds a reference to its representation.
  }
  else if (arraytype==mxDOUBLE_CLASS && rows*cols*layers == 1) {
    if (DEBUG_PTBOCTAVEGLUE) printf("NEW SCALAR:\n"); fflush(NULL);
    // This is a scalar value:
    retval->o = (void*) new octave_value(0.0);
    double* dp = (double*) PsychMallocTemp(sizeof(double));
    retval->d = (void*) dp;
  }
  else if (arraytype==mxDOUBLE_CLASS && rows*cols*layers == 0) {
    // Special case: Empty matrix.
    if (DEBUG_PTBOCTAVEGLUE) printf("NEW EMPTY DOUBLE MATRIX:\n"); fflush(NULL);
    retval->o = (void*) new octave_value(Matrix(0,0));
    retval->d = NULL;
  }
  else if (arraytype==mxLOGICAL_CLASS) {
    if (DEBUG_PTBOCTAVEGLUE) printf("NEW BOOLMATRIX: %i, %i\n", rows, cols, layers); fflush(NULL);
    // Create empty double-matrix of type mxREAL...
    if (layers>1) PsychErrorExitMsg(PsychError_internal, "In mxCreateNumericArray: Tried to allocate a 3D psych_bool matrix!?! Unsupported.");

    boolMatrix m(rows, cols);
    // Retrieve a pointer to internal representation. As m is new
    // this won't trigger a deep-copy.
    retval->d = (void*) m.data();
    // Build a new oct_value object from Matrix m: This is a
    // shallow-copy.
    octave_value* ovp = new octave_value();
    *ovp = m;
    retval->o = (void*) ovp;
    // At this point we can safely destroy Matrix m, as the new
    // octave_object holds a reference to its representation.
  }
  else {
    PsychErrorExitMsg(PsychError_internal, "FATAL Error: mxCreateNumericArray: Unknown matrix type requested!");
  }

  return(retval);
}

mxArray* mxCreateDoubleMatrix(int rows, int cols, int realorcomplex)
{
  int dims[2];
  dims[0]=rows;
  dims[1]=cols;
  return(mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL));  
}

mxArray* mxCreateLogicalMatrix(int rows, int cols)
{
  int dims[2];
  dims[0]=rows;
  dims[1]=cols;
  return(mxCreateNumericArray(2, dims, mxLOGICAL_CLASS, mxREAL));
}

mxArray* mxCreateString(const char* instring)
{
  mxArray* retval;

  // Allocate our mxArray-Struct:
  retval = (mxArray*) PsychMallocTemp(sizeof(mxArray));
  
  retval->d = NULL;
  retval->o = (void*) new octave_value(instring);
  return(retval);
}

void* mxGetData(const mxArray* arrayPtr)
{
  return(arrayPtr->d);
}

double mxGetScalar(const mxArray* arrayPtr)
{
  return((double) (mxGetPr(arrayPtr)[0]));
}

double* mxGetPr(const mxArray* arrayPtr)
{
  return((double*) mxGetData(arrayPtr));
}

psych_bool* mxGetLogicals(const mxArray* arrayPtr)
{
  return((psych_bool*) mxGetData(arrayPtr));
}
#define GETOCTPTR(x) ((octave_value*) (x)->o)

int mxGetM(const mxArray* arrayPtr)
{
  return(GETOCTPTR(arrayPtr)->rows());
}

int mxGetN(const mxArray* arrayPtr)
{
  return(GETOCTPTR(arrayPtr)->columns());
}

int mxGetNumberOfDimensions(const mxArray* arrayPtr)
{
  return(GETOCTPTR(arrayPtr)->ndims());
}

int* mxGetDimensions(const mxArray* arrayPtr)
{
  dim_vector vdim = GETOCTPTR(arrayPtr)->dims();
  int* dims = (int*) PsychMallocTemp(3*sizeof(int));
  dims[0] = vdim(0);
  dims[1] = vdim(1);
  if (GETOCTPTR(arrayPtr)->ndims()>2) dims[2] = vdim(2); else dims[2]=1;
  return(dims);
}

int mxGetString(const mxArray* arrayPtr, char* outstring, int outstringsize)
{
  if (!mxIsChar(arrayPtr)) PsychErrorExitMsg(PsychError_internal, "FATAL Error: Tried to convert a non-string into a string!");
  return(((snprintf(outstring, outstringsize, "%s", GETOCTPTR(arrayPtr)->string_value().c_str()))>=0) ? 0 : 1);
}

void mxDestroyArray(mxArray *arrayPtr)
{
  // Destroy a mxArray:
  if (arrayPtr == NULL) return;

  // We only need to destroy the octave_value object referenced by arrayPtr,
  // because possible data buffers referenced by the ->d field and the
  // mxArray struct itself are allocted via PsychMallocTemp() anyway, so
  // they get automatically released when exiting our octFile...
  octave_value* ov = (octave_value*) arrayPtr->o;
  if (ov) delete(ov);
  arrayPtr->o = NULL;
  return;
}

mxArray* mxCreateStructArray(int numDims, int* ArrayDims, int numFields, const char** fieldNames)
{
  mxArray* retval;

  if (numDims>2 || numDims<1) PsychErrorExitMsg(PsychError_unimplemented, "FATAL Error: mxCreateStructArray: Anything else than 1D or 2D Struct-Arrays is not supported!");
  if (numFields<1) PsychErrorExitMsg(PsychError_internal, "FATAL Error: mxCreateStructArray: numFields < 1 ?!?");

  // Our dimension vector:
  dim_vector mydims((numDims>1) ? dim_vector(ArrayDims[0], ArrayDims[1]) : dim_vector(ArrayDims[0]));

  // Create cell array of requested dimensionality and size as template. This
  // will be cloned numFields - times to create one such array per field.
  class Cell myCell(mydims);

  // Create an Octave_map(): A Octave_map is an associative map that associates keys (namestrings)
  // with values (which are Cell-Arrays). For each named field in our struct array, we insert a
  // key->value pair into the Octave_map. The key is the field name, the value is a copy of our
  // cell array 'myCell', where myCell has the requested dimensions.
  // Lookup or assignment will work as follows: The fieldName is used as search-key to lookup the
  // cell array associated with that key. The fieldIndex is then used to index into that cell array
  // and get or set the octave_value object which encodes the value of that field at that index...
  Octave_map mymap(std::string(fieldNames[0]), myCell);  // Assign first key in constructor.

  // Assign all remaining key->Cell pairs...
  for (int i=1; i<numFields; i++) mymap.assign(std::string(fieldNames[i]), myCell);

  // Build our mxArray wrapper:
  retval = (mxArray*) PsychMallocTemp(sizeof(mxArray));

  // Fill it: Assign our map.
  octave_value* ovp = new octave_value(mymap);
  retval->o = (void*) ovp;
  retval->d = NULL;
  
  return(retval);
}

int mxGetFieldNumber(const mxArray* structArray, const char* fieldName)
{
  if(!mxIsStruct(structArray)) {
    PsychErrorExitMsg(PsychError_internal, "FATAL Error: mxGetFieldNumber: Tried to manipulate something other than a struct-Array!");
  }

  // Retrieve map:
  octave_value* ov = (octave_value*) structArray->o;
  Octave_map om = ov->map_value();

  // Find idx of key.
  std::string mykey(fieldName);
  for(int i=0; i<om.length(); i++) if (mykey == om.keys()(i)) return(i);

  // No such key :(
  return(-1);
}

void mxSetField(mxArray* pStructOuter, int index, const char* fieldName, mxArray* pStructInner)
{
  if(!mxIsStruct(pStructOuter)) {
    PsychErrorExitMsg(PsychError_internal, "FATAL Error: mxSetField: Tried to manipulate something other than a struct-Array!");
  }

  // Retrieve map:
  octave_value* ov = (octave_value*) pStructOuter->o;
  Octave_map om = ov->map_value();

  // Retrieve object:
  octave_value* iv = (octave_value*) pStructInner->o;
  if (iv->is_real_type() && iv->is_scalar_type()) {
    // Our special case. Do the extra work...
    *iv=octave_value(*((double*) pStructInner->d));
  }

  // Assign our object:
  om.contents(std::string(fieldName))(index) = *iv;

  // Reassign our map:
  *ov = octave_value(om);
}

mxArray* mxCreateCellArray(int numDims, int* ArrayDims)
{
  mxArray* retval;

  // Allocate our mxArray-Struct:
  retval = (mxArray*) PsychMallocTemp(sizeof(mxArray));

  if (numDims>2) PsychErrorExitMsg(PsychError_unimplemented, "FATAL Error: mxCreateCellArray: 3D Cell Arrays are not supported yet on GNU/Octave build!");

  // Create dimension vector:
  dim_vector mydims((numDims>1) ? dim_vector(ArrayDims[0], ArrayDims[1]) : dim_vector(ArrayDims[0]));
  
  // Create Cell object:
  Cell myCell(mydims);
  retval->o = (void*) new octave_value(myCell);
  retval->d = NULL;

  // Done.
  return(retval);
}

void mxSetCell(PsychGenericScriptType *cellVector, int index, mxArray* mxFieldValue)
{
  if (!mxIsCell((mxArray*) cellVector)) {
    PsychErrorExitMsg(PsychError_internal, "FATAL Error: mxSetCell: Tried to manipulate something other than a cell-vector!");
  }

  // Get a local (shallow) copy of the current real cellVector:
  octave_value* cv = (octave_value*) cellVector->o;
  Cell mycell = cv->cell_value();

  // Assign new mxFieldValue:
  octave_value* ov = (octave_value*) mxFieldValue->o;
  if (ov->is_real_type() && ov->is_scalar_type()) {
    // Our special case. Do the extra work...
    *ov=octave_value(*((double*) mxFieldValue->d));
  }

  mycell(index)=*ov;

  // Assign modified vector:
  *cv = mycell;

  return;
}

void mxSetLogical(mxArray* dummy)
{
  // This is a no-op on Octave build, because it is not needed anywhere...
  return;
}

int mxIsLogical(const mxArray* a)
{
  return(GETOCTPTR(a)->is_bool_type());
}

int mxIsCell(const mxArray* a)
{
  return(GETOCTPTR(a)->is_cell());
}

int mxIsStruct(const mxArray* a)
{
  return(GETOCTPTR(a)->is_map());
}

int mxIsNumeric(const mxArray* a)
{
  return(GETOCTPTR(a)->is_numeric_type());
}

int mxIsChar(const mxArray* a)
{
  return(GETOCTPTR(a)->is_char_matrix() || GETOCTPTR(a)->is_string());
}

int mxIsDouble(const mxArray* a)
{
  return(GETOCTPTR(a)->is_real_type());
}

int mxIsUint8(const mxArray* a)
{
  if (DEBUG_PTBOCTAVEGLUE) printf("BYTESIZE %i NUMEL %i --> mxIsUint8 = %s\n", GETOCTPTR(a)->byte_size(), GETOCTPTR(a)->numel(),
	 ((mxIsNumeric(a) && (GETOCTPTR(a)->byte_size() / GETOCTPTR(a)->numel() == 1))) ? "TRUE":"FALSE");

  return((mxIsNumeric(a) && (GETOCTPTR(a)->byte_size() / GETOCTPTR(a)->numel() == 1)) ? TRUE : FALSE);
}

int mxIsUint16(const mxArray* a)
{
  return((mxIsNumeric(a) && (GETOCTPTR(a)->byte_size() / GETOCTPTR(a)->numel() == 2)) ? TRUE : FALSE);
  //  int t=GETOCTPTR(a)->type_id();
  //return(t==octave_uint16_matrix::static_type_id());
}

int mxIsUint32(const mxArray* a)
{
  return((mxIsNumeric(a) && (GETOCTPTR(a)->byte_size() / GETOCTPTR(a)->numel() == 4)) ? TRUE : FALSE);

  //  int t=GETOCTPTR(a)->type_id();
  // return(t==octave_uint32_matrix::static_type_id());
}

int mxIsInt8(const mxArray* a)
{
  return(mxIsUint8(a));
  //  int t=GETOCTPTR(a)->type_id();
  //return(t==octave_int8_matrix::static_type_id());
}

int mxIsInt16(const mxArray* a)
{
  return(mxIsUint16(a));
  //  int t=GETOCTPTR(a)->type_id();
  //return(t==octave_int16_matrix::static_type_id());
}

int mxIsInt32(const mxArray* a)
{
  return(mxIsUint32(a));
  //int t=GETOCTPTR(a)->type_id();
  //return(t==octave_int32_matrix::static_type_id());
}

#endif

////Static functions local to ScriptingGlue.c.  
#if PSYCH_LANGUAGE == PSYCH_MATLAB
void InitializeSynopsis(char *synopsis[],int maxStrings);
#endif 



// _____________________________________________________________________________________
// for Matlab
#if PSYCH_LANGUAGE == PSYCH_MATLAB | PSYCH_LANGUAGE == PSYCH_OCTAVE

#define MAX_SYNOPSIS 100
#define MAX_CMD_NAME_LENGTH 100

//Static variables local to ScriptingGlue.c.  The convention is to append a abbreviation in all
//caps of the C file name to the variable name.   
int nlhsGLUE;  // Number of requested return arguments.
int nrhsGLUE;  // Number of provided call arguments.

#if PSYCH_LANGUAGE == PSYCH_MATLAB
static mxArray **plhsGLUE;       // A pointer to the plhs array passed to the MexFunction entry point  
static CONSTmxArray **prhsGLUE; // A pointer to the prhs array passed to the MexFunction entry point
#endif

#if PSYCH_LANGUAGE == PSYCH_OCTAVE
#define MAX_OUTPUT_ARGS 100
#define MAX_INPUT_ARGS 100
static mxArray* plhsGLUE[MAX_OUTPUT_ARGS]; // An array of pointers to the octave return arguments.
static mxArray* prhsGLUE[MAX_INPUT_ARGS];  // An array of pointers to the octave call arguments.
extern const char *mexFunctionName; // This gets initialized by Octave wrapper to contain our function name.
#endif

static psych_bool nameFirstGLUE;
//static PsychFunctionPtr psychExitFunctionGLUE=NULL; 
static psych_bool subfunctionsEnabledGLUE=FALSE;
static psych_bool baseFunctionInvoked=FALSE;
static void PsychExitGlue(void);

//local function declarations
static psych_bool PsychIsEmptyMat(CONSTmxArray *mat);
static psych_bool PsychIsDefaultMat(CONSTmxArray *mat);
static mwSize mxGetP(const mxArray *array_ptr);
static mwSize mxGetNOnly(const mxArray *arrayPtr);
static mxArray *mxCreateDoubleMatrix3D(psych_int64 m, psych_int64 n, psych_int64 p);

//declarations for functions exported from code module
#if PSYCH_LANGUAGE == PSYCH_MATLAB
EXP void mexFunction(int nlhs, mxArray *plhs[], int nrhs, CONSTmxArray *prhs[]);
#endif

// firstTime: This flag defines if this is the first invocation of the module
// since it was (re-)loaded:
static psych_bool firstTime = TRUE;

#if PSYCH_LANGUAGE == PSYCH_OCTAVE
PsychError PsychExitOctaveGlue(void);
static psych_bool jettisoned = FALSE;
#endif

// Forward declaration for GNU/Linux compile:
#ifdef __cplusplus
extern "C" void ScreenCloseAllWindows(void);
#else
void ScreenCloseAllWindows(void);
#endif

/*

	Main entry point for Matlab and Octave. Serves as a dispatch and handles
	first time initialization.
	
	EXP is a macro defined within Psychtoolbox source to be nothing
	except on win where it is the declaration which tells the linker to 
	make the function visible from outside the DLL. 
        
        The subfunction dispatcher can operate in either of two modes depending
        on whether the module has registed subfunctions, or only a single "base" 
        function.  
        
        subfunction mode:  
        The examines the  first and second 
        arguments for a string naming a module subfunction.  If it finds in either of those
        two arguments a string naming a module subfunctoin, then it looks up the approproate 
        function pointer and invokes that function.  Before invoking the function the dispatcher
        removes the function name argument form the list of argumnets which was passed to the 
        module.  
                
        base mode:  The dispatcher always invokes the same one subfunction and without
        alterinng the list of arguments.
        
        Modules should now register in subfunction mode to support the build-in 'version' command.
        
*/
#if PSYCH_LANGUAGE == PSYCH_MATLAB
EXP void mexFunction(int nlhs, mxArray *plhs[], int nrhs, CONSTmxArray *prhs[])
#endif
#if PSYCH_LANGUAGE == PSYCH_OCTAVE
EXP octave_value_list octFunction(const octave_value_list& prhs, const int nlhs)
#endif
{
	psych_bool errorcondition = FALSE;
	psych_bool isArgThere[2], isArgEmptyMat[2], isArgText[2], isArgFunction[2];
	PsychFunctionPtr fArg[2], baseFunction;
	char argString[2][MAX_CMD_NAME_LENGTH];
	int i; 
	const mxArray* tmparg = NULL; // mxArray is mxArray under MATLAB but #defined to octave_value on OCTAVE build.
	#if PSYCH_LANGUAGE == PSYCH_OCTAVE
	  // plhs is our octave_value_list of return values:
	  octave_value tmpval;      // Temporary, needed in parser below...
	  octave_value_list plhs;   // Our list of left-hand-side return values...
	  int nrhs = prhs.length();

	  // Child protection: Is someone trying to call us after we've shut down already?
	  if (jettisoned) {
	    // Yep! Stupido...
	    error("%s: Tried to call the module after it has been jettisoned!!! You need to do a 'clear %s;' now. Bug in Psychtoolbox?!?",
		  mexFunctionName, mexFunctionName);
	    return(plhs);
	  }
	#endif
 
	#if PSYCH_LANGUAGE == PSYCH_OCTAVE
	  // Save CPU-state and stack at this position in 'jmpbuffer'. If any further code
	  // calls an error-exit function like PsychErrorExit() or PsychErrorExitMsg() then
	  // the corresponding longjmp() call in our mexErrMsgTxt() implementation (see top of file)
	  // will unwind the stack and restore stack-state and CPU state to the saved values in
	  // jmpbuffer --> We will end up at this setjmp() call again, with a cleaned up stack and
	  // CPU state, but setjmp will return a non-zero error code, signaling the abnormal abortion.
	  if (setjmp(jmpbuffer)!=0) {
	    // PsychErrorExit() or friends called! The CPU and stack are restored to a sane state.
	    // Call our cleanup-routine to release memory that is PsychMallocTemp()'ed and to other
	    // error-handling...
	    errorcondition = TRUE;
	    goto octFunctionCleanup;
	  }
	#endif

	// Initialization
	if (firstTime) {
		
		//call the Psychtoolbox init function, which inits the Psychtoolbox and calls the project init. 
		PsychInit();
		
		//register the exit function, which calls PsychProjectExit() to clean up for the project then
		//calls whatever to clean up for all of Psych.h layer.

		#if PSYCH_LANGUAGE == PSYCH_MATLAB
		// Under Matlab we use the mexAtExit() Mex-API function to register our PsychExitGlue() routine.
		// Whenever Matlab wants to flush our module (to reload it, or in response to Matlab-Shutdown,
		// 'clear MODULENAME', 'clear mex' or 'clear all' command) it first calls our PsychExitGlue(),
		// then unloads the module from memory...
		mexAtExit(&PsychExitGlue);
		#endif

		#if PSYCH_LANGUAGE == PSYCH_OCTAVE
		// Octave (as of Version 2.1.73) does not seem to support a way to register such a
		// cleanup handler, so we use the following trick: We tell octave to lock our OCT file
		// into memory, so it can not be clear'ed out of memory by Octave with the standard clear
		// command. Then we register a new Module subfunction 'JettisonModuleHelper': If this
		// subcommand is called, it will call our PsychExitGlue() cleanup routine, then unlock
		// ourselves from memory, now that it is safe to flush us. We provide special scripts
		// clearall.m, clearoct.m, clearmex.m and clearMODULENAME.m that do what clear all,
		// clear mex and clear MODULENAME would do on Matlab, by simply calling the
		// MODULENAME('JettisonModuleHelper'); function, followed by a clear MODULENAME; command.
		// --> User has same functionality with nearly same syntax and should be safe on Octave
		// as well.
		PsychRegister("JettisonModuleHelper",  &PsychExitOctaveGlue);

		// Lock ourselves into Octaves runtime environment so we can't get clear'ed out easily:
		mlock(std::string(mexFunctionName));
		#endif

		// Register hidden helper function: This one dumps all registered subfunctions of
		// a module into a struct array of text strings. Needed by our automatic documentation
		// generator script to find out about subfunctions of a module:
		PsychRegister((char*) "DescribeModuleFunctionsHelper",  &PsychDescribeModuleFunctions);

		firstTime = FALSE;
	}
	
	// Store away call arguments for use by language-neutral accessor functions in ScriptingGlue.c
	#if PSYCH_LANGUAGE == PSYCH_MATLAB
	nlhsGLUE = nlhs;
	nrhsGLUE = nrhs;
	plhsGLUE = plhs;
	prhsGLUE = prhs;
	#endif

	#if PSYCH_LANGUAGE == PSYCH_OCTAVE

	// NULL-init our pointer array of call value pointers prhsGLUE:
	memset(&prhsGLUE[0], 0, sizeof(prhsGLUE));

	// Setup our prhsGLUE array of call argument pointers:
	// We make copies of prhs to simplify the rest of PsychScriptingGlue. This copy is not
	// as expensive as it might look, because Octave objects are all implemented via
	// "Copy-on-write" --> Only a pointer is copied as long as we don't modify the data.
	// MK: TODO FIXME -- Should we keep an extra array octave_value dummy[MAX_INPUT_ARGS];
	// around, assign to that dummy[i]=prhs(i); and set ptrs to it prhsGLUE[i]=&dummy[i];
	// This would require more memory, but prevent possible heap-fragmentation due to
	// lots of new()/delete() calls on each invocation of the OCT-Function --> possible
	// slow-down over time, could be confused with memory leaks???
	for(int i=0; i<nrhs && i<MAX_INPUT_ARGS; i++) {
	  // Create and assign our mxArray-Struct:
	  prhsGLUE[i] = (mxArray*) PsychMallocTemp(sizeof(mxArray));

	  // Extract data-pointer to each prhs(i) octave_value and store a type-casted version
	  // which is optimal for us.
	  if (prhs(i).is_string() || prhs(i).is_char_matrix()) {
	    // A string object:
	    if (DEBUG_PTBOCTAVEGLUE) printf("INPUT %i: STRING\n", i); fflush(NULL);

	    // Strings do not have a need for a data-ptr. Just copy the octave_value object...
	    prhsGLUE[i]->d = NULL;
	    prhsGLUE[i]->o = (void*) new octave_value(prhs(i));  // Refcont now >= 2
	    // Done.
	  } 
	  else if (prhs(i).is_real_type() && !prhs(i).is_scalar_type()) {
	    // A N-Dimensional Array:
	    if (DEBUG_PTBOCTAVEGLUE) printf("TYPE NAME %s\n", prhs(i).type_name().c_str()); fflush(NULL);

	    // Is it an uint8 or int8 NDArray?
	    if (strstr(prhs(i).type_name().c_str(), "int8")) {
	      // Seems to be an uint8 or int8 NDArray: Create an optimized uint8 object of it:
	      if (DEBUG_PTBOCTAVEGLUE) printf("INPUT %i: UINT8-MATRIX\n", i); fflush(NULL);

	      // Create intermediate representation m: This is a shallow-copy...
	      const uint8NDArray m(prhs(i).uint8_array_value()); // Refcount now >=2

	      // Get internal dataptr from it:        // This triggers a deep-copy :(
	      prhsGLUE[i]->d = (void*) m.data();      // Refcount now == 1
	      
	      // Create a shallow backup copy of corresponding octave_value...
	      octave_value* ovptr = new octave_value();
	      *ovptr = m;
	      prhsGLUE[i]->o = (void*) ovptr;  // Refcont now == 2
	      
	      // As soon as m gets destructed by leaving this if-branch,
	      // the refcount will drop to == 1...
	      
	      // Done.
	    }
	    else {
	      // Seems to be a non-uint8 NDArray, i.e. psych_bool type or double type.
	      if (DEBUG_PTBOCTAVEGLUE) printf("INPUT %i: DOUBLE-MATRIX\n", i); fflush(NULL);

	      // We create a generic double NDArray from it...

	      // Create intermediate representation m: This is a shallow-copy...
	      const NDArray m(prhs(i).array_value()); // Refcount now >=2

	      // Get internal dataptr from it:        // This triggers a deep-copy :(
	      prhsGLUE[i]->d = (void*) m.data();      // Refcount now == 1
	      
	      // Create a shallow backup copy of corresponding octave_value...
	      octave_value* ovptr = new octave_value();
	      *ovptr = m;
	      prhsGLUE[i]->o = (void*) ovptr;  // Refcont now == 2
	      
	      // As soon as m gets destructed by leaving this if-branch,
	      // the refcount will drop to == 1...
	      
	      // Done.
	    }
	  } else if (prhs(i).is_real_type() && prhs(i).is_scalar_type()) {

	    // A double or integer scalar value:
	    if (DEBUG_PTBOCTAVEGLUE) printf("INPUT %i: SCALAR\n", i); fflush(NULL);
	    prhsGLUE[i]->o = (void*) new octave_value(prhs(i));
	    // Special case: We allocate our own double value and store a
	    // copy of the value in it.
	    double* m = (double*) PsychMallocTemp(sizeof(double));
	    *m = prhs(i).double_value();
	    prhsGLUE[i]->d = (void*) m;
	  }
	  else {
	    // Unkown argument type that we can't handle :(
	    // We abort with a reasonable error message:
	    prhsGLUE[i]=NULL;
	    // We do, however, give an extra warning, as this could be Octave related...
	    printf("PTB-WARNING: One of the values in the argument list was not recognized.\n");
	    printf("PTB-WARNING: If your script runs well on Matlab then this may be a limitation or\n");
	    printf("PTB-WARNING: bug in the GNU/Octave version of Psychtoolbox :( ...\n");
	    PsychErrorExitMsg(PsychError_unimplemented, "Unrecognized argument in list of command parameters.");
	  }
	}

	// NULL-out our pointer array of return value pointers plhsGLUE:
	memset(&plhsGLUE[0], 0, sizeof(plhsGLUE));

	nrhsGLUE = nrhs;
	nlhsGLUE = nlhs;
	#endif

        baseFunctionInvoked=FALSE;

	//if no subfunctions have been registered by the project then just invoke the project base function
	//if one of those has been registered.
	if(!PsychAreSubfunctionsEnabled()){
		baseFunction = PsychGetProjectFunction(NULL);
		if(baseFunction != NULL){
                        baseFunctionInvoked=TRUE;
			(*baseFunction)();  //invoke the unnamed function
		}else
			PrintfExit("Project base function invoked but no base function registered");
	}else{ //subfunctions are enabled so pull out the function name string and invoke it.
		//assess the nature of first and second arguments for finding the name of the sub function.  
		for(i=0;i<2;i++)
		{
			isArgThere[i] = (nrhs>i) && (prhsGLUE[i]);
			#if PSYCH_LANGUAGE == PSYCH_MATLAB
			if (isArgThere[i]) tmparg = prhs[i]; else tmparg = NULL;
			#endif

			#if PSYCH_LANGUAGE == PSYCH_OCTAVE
			if (isArgThere[i]) { tmparg = prhsGLUE[i]; } else { tmparg = NULL; }
			#endif

			isArgEmptyMat[i] = isArgThere[i] ? mxGetM(tmparg)==0 || mxGetN(tmparg)==0 : FALSE;  
			isArgText[i] = isArgThere[i] ? mxIsChar(tmparg) : FALSE;
			if(isArgText[i]){
				mxGetString(tmparg,argString[i],sizeof(argString[i]));
				fArg[i]=PsychGetProjectFunction(argString[i]);
			}
			isArgFunction[i] = isArgText[i] ? fArg[i] != NULL : FALSE;
		}
		 


		//figure out which of the two arguments might be the function name and either invoke it or exit with error
		//if we can't find one.  

		if(!isArgThere[0] && !isArgThere[1]){ //no arguments passed so execute the base function 	
			baseFunction = PsychGetProjectFunction(NULL);
			if(baseFunction != NULL){
                                baseFunctionInvoked=TRUE;
				(*baseFunction)();
			}else
				PrintfExit("Project base function invoked but no base function registered");
		}
		// (!isArgThere[0] && isArgEmptyMat[1]) --disallowed
		// (!isArgThere[0] && isArgText[1])     --disallowed
		// (!isArgThere[0] && !isArgText[1]     --disallowed except in case of !isArgThere[0] caught above. 

		else if(isArgEmptyMat[0] && !isArgThere[1])
			PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state A)");
		else if(isArgEmptyMat[0] && isArgEmptyMat[1])
			PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state B)");
		else if(isArgEmptyMat[0] && isArgText[1]){
			if(isArgFunction[1]){
				nameFirstGLUE = FALSE;
				(*(fArg[1]))();
			}
			else
				PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state C)");
		}
		else if(isArgEmptyMat[0] && !isArgText[1])
			PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state D)");
			
		else if(isArgText[0] && !isArgThere[1]){
			if(isArgFunction[0]){
				nameFirstGLUE = TRUE;
				(*(fArg[0]))();
			}else{ //when we receive a first argument  wich is a string and it is  not recognized as a function name then call the default function 
			/*
                        else
				PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state E)");
                        */
                            baseFunction = PsychGetProjectFunction(NULL);
                            if(baseFunction != NULL){
                                baseFunctionInvoked=TRUE;
				(*baseFunction)();
                            }else
				PrintfExit("Project base function invoked but no base function registered");
                        }
                            
		}
		else if(isArgText[0] && isArgEmptyMat[1]){
			if(isArgFunction[0]){
				nameFirstGLUE = TRUE;
				(*(fArg[0]))();
			}
			else
				PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state F)");
		}
		else if(isArgText[0] && isArgText[1]){
			if(isArgFunction[0] && !isArgFunction[1]){ //the first argument is the function name
				nameFirstGLUE = TRUE;
				(*(fArg[0]))();
			}
			else if(!isArgFunction[0] && isArgFunction[1]){ //the second argument is the function name
				nameFirstGLUE = FALSE;
				(*(fArg[1]))();
			}
			else if(!isArgFunction[0] && !isArgFunction[1]){ //neither argument is a function name
                            //PrintfExit("Invalid command (error state G)");
                            baseFunction = PsychGetProjectFunction(NULL);
                            if(baseFunction != NULL){
                                baseFunctionInvoked=TRUE;
				(*baseFunction)();
                            }else
				PrintfExit("Project base function invoked but no base function registered");
                        }
			else if(isArgFunction[0] && isArgFunction[1]) //both arguments are function names
				PrintfExit("Passed two function names");
		}
		else if(isArgText[0] && !isArgText[1]){
			if(isArgFunction[0]){
				nameFirstGLUE = TRUE;
				(*(fArg[0]))();
			}
			else
				PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state H)");
		}

		else if(!isArgText[0] && !isArgThere[1]){  //this was modified for MODULEVersion with WaitSecs.
                    //PrintfExit("Invalid command (error state H)");
                    baseFunction = PsychGetProjectFunction(NULL);
                    if(baseFunction != NULL){
                        baseFunctionInvoked=TRUE;
                        (*baseFunction)();  //invoke the unnamed function
                    }else
                        PrintfExit("Project base function invoked but no base function registered");
                }
		else if(!isArgText[0] && isArgEmptyMat[1])
			PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state I)");
		else if(!isArgText[0] && isArgText[1])
		{
			if(isArgFunction[1]){
				nameFirstGLUE = FALSE;
				(*(fArg[1]))();
			}
			else
				PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state J)");
		}
		else if(!isArgText[0] && !isArgText[1]){  //this was modified for Priority.
                    //PrintfExit("Invalid command (error state K)");
                    baseFunction = PsychGetProjectFunction(NULL);
                    if(baseFunction != NULL){
                        baseFunctionInvoked=TRUE;
                        (*baseFunction)();  //invoke the unnamed function
                    }else
                        PrintfExit("Project base function invoked but no base function registered");
                }

	} //close else			

	// If we reach this point of execution under Matlab, then we're successfully done.

#if PSYCH_LANGUAGE == PSYCH_OCTAVE
	// If we reach this point of execution under Octave, then we're done, possibly due to
	// error abort. Let's first do the memory management cleanup work necessary on Octave.
	// This is either done due to successfull execution or via jump to octFunctionCleanup:
	// in case of error-abort:
octFunctionCleanup:

	// Release our own prhsGLUE array...
	// Release memory for scalar types is done by PsychFreeAllTempMemory(); 
	for(int i=0; i<nrhs && i<MAX_INPUT_ARGS; i++) if(prhsGLUE[i]) {
	  delete(((octave_value*)(prhsGLUE[i]->o)));
	  prhsGLUE[i]=NULL;	  
	}

	// "Copy" our octave-value's into the output array: If nlhs should be
	// zero (Octave-Script does not expect any return arguments), but our
	// subfunction has assigned a return argument in slot 0 anyway, then
	// we return that argument and release our own temp-memory. This
	// provides Matlab-semantic, where unsolicited return arguments are
	// printed anyway as content of the "ans" variable.
	for(i=0; (i==0 && plhsGLUE[0]!=NULL) || (i<nlhs && i<MAX_OUTPUT_ARGS); i++) {
	  if (plhsGLUE[i]) {
	    plhs(i) = *((octave_value*)(plhsGLUE[i]->o));
	    if (plhs(i).is_scalar_type()) {
	      // Special case: Scalar. Need to override with our double-ptrs value:
	      double* svalue = (double*) plhsGLUE[i]->d;
	      plhs(i) = octave_value((double) *svalue);
	    }

	    // Delete our own octave_value object. All relevant data has been
	    // copied via "copy-on-write" into plhs(i) already:
 	    delete(((octave_value*)(plhsGLUE[i]->o)));

	    // We don't need to free() the PsychMallocTemp()'ed object pointed to
	    // by the d-Ptr, nor do we need to free the mxArray-Struct. This is done
	    // below in PsychFreeAllTempMemory(). Just NULL-out the array slot:
	    plhsGLUE[i]=NULL;
	  }
	}

	// Release all memory allocated via PsychMallocTemp():
	PsychFreeAllTempMemory();

	// Is this a successfull return?
	if (errorcondition) {
	  // Nope - Error return, either due to some PTB detected error or due to
	  // the user pressing the CTRL+C key combo. Try to call PTB's
	  // Screen('CloseAll') command to close the display, at least if this is
	  // the Screen module.

		// Is this the Screen() module?
		if (strcmp(PsychGetModuleName(), "Screen")==0) {
			// Yes. We directly call our close and cleanup routine:
			#ifdef PTBMODULE_Screen
			ScreenCloseAllWindows();
			#endif
		} else {
			// Nope. This is a Psychtoolbox OCT file other than Screen.
			// We can't call directly, but we can call the 'sca' command
			// from Octave:
			PsychRuntimeEvaluateString("Screen('CloseAll');");
		}
	}

	// Return our octave_value_list of returned values in any case and yield control
	// back to Octave:
	return(plhs);
#endif
}

#if PSYCH_LANGUAGE == PSYCH_OCTAVE

/*      Call PsychExitGlue(), followed by unlocking the module:
	- Needed to safely remove modules on GNU/Octave.
*/
PsychError PsychExitOctaveGlue(void)
{
	// Debug output:
	if (DEBUG_PTBOCTAVEGLUE) printf("PTB-INFO: Jettisoning submodule %s ...\n", mexFunctionName); fflush(NULL);

	// Call our regular exit routines to clean up and release all ressources:
	PsychExitGlue();

	// Mark ourselves (via global variable "jettisoned") as shut-down. Any
	// further invocations of the module without previously clear'ing and
	// reloading it will be prevented.
	jettisoned = TRUE;

	// Unlock ourselves from Octaves runtime environment so we can get clear'ed out:
	munlock(std::string(mexFunctionName));

	// Done. Return control to Octave - It will now remove us from its process-space - RIP.
	return(PsychError_none);
}

#endif

/*
	Just call the abstracted PsychExit function.  This might seem dumb, but its necessary to 
	isolate the scripting language dependent stuff from the rest of the toolbox.  
	
*/
void PsychExitGlue(void)
{
    // Reset firstTime flag, so reloading the module leads to a full init cycle.
    // This is especially important on Octave, as it doesn't fully unload/reload
    // mex file modules at "clear" time, so the static firstTime variable won't
    // get automatically reset to its inital load default of TRUE:
    // firstTime = TRUE;
    
    // Perform platform independent shutdown:
	PsychErrorExitMsg(PsychExit(),NULL);
    
    // And we are dead! Now the runtime will flush us from process memory,
    // at least on Matlab. In any case no further invocation will happen
    // until reload.
}


/*
	Return the mxArray pointer to the specified position.  Note that we have some special rules for 
	for numbering the positions: 
	
	0 - This is always the command string or NULL if the project does not register a 
	    dispatch function and does accept subcommands.  If the function does accept sub
	    commands, in Matlab those may be passed in either the first or second position, but
	    PsychGetArgPtr() will always return the command as the 0th. 
	    
	1 - This is the first argument among the arguments which are not the subfunction name itself.
	    It can occur in either the first or second position of the argument list, depending on
	    in which of those two positions the function name itself appears.
	    
	2.. These positions are numbered correctly 
	
	TO DO:  
	
	
	2 - this function should be used by the one which gets the function name.    
		

	Arguments are numbered 0..n.  

		-The 0th argument is a pointer to the mxArray holding
		the subfunction name string if we are in subfunction mode.  
	
		-The 0th argument is undefined if not in subfunction mode.  
		
		-The 1st argument is the argument of the 1st and 2nd which is not
		 the subfunction name if in subfunction mode.
		 
		-The 1st argument is the first argument if not in subfunction mode.
		
		-The 2nd-nth arguments are always the 2nd-nth arguments. 		
*/
//we return NULL if a postion without an arg is specified.
const mxArray *PsychGetInArgMxPtr(int position)
{	


	if(PsychAreSubfunctionsEnabled() && !baseFunctionInvoked){ //when in subfunction mode
		if(position < nrhsGLUE){ //an argument was passed in the correct position.
			if(position == 0){ //caller wants the function name argument.
				if(nameFirstGLUE)
					return(prhsGLUE[0]);
				else
					return(prhsGLUE[1]);
			}else if(position == 1){ //they want the "first" argument.    
				if(nameFirstGLUE)
					return(prhsGLUE[1]);
				else
					return(prhsGLUE[0]);
			}else
				return(prhsGLUE[position]);
		}else
			return(NULL); 
	}else{ //when not in subfunction mode and the base function is not invoked.  
		if(position <= nrhsGLUE)
			return(prhsGLUE[position-1]);
		else
			return(NULL);
	}
}

mxArray **PsychGetOutArgMxPtr(int position)
{	

	if(position==1 || (position>0 && position<=nlhsGLUE)){ //an ouput argument was supplied at the specified location
		return(&(plhsGLUE[position-1]));
	}else
		return(NULL);
}

/*
	functions for enabling and testing subfunction mode
*/
void PsychEnableSubfunctions(void)
{
	subfunctionsEnabledGLUE = TRUE;
}


psych_bool PsychAreSubfunctionsEnabled(void)
{
	return(subfunctionsEnabledGLUE);
}

/* PsychCheckmWSizeLimits(size_t m, size_t n, size_t p)
 *
 * Makes sure matrix/vector dimensions stay within the limits imposed
 * by a specific Psychtoolbox build for a specific 32 bit or 64 bit
 * version of Matlab or Octave.
 *
 * 32 bit builds are always limited to 2^31-1 elements per dimension at
 * most. In reality the limit is much lower, this is an upper bound.
 *
 * 64 bit builds may or may not allow bigger limits, depending if
 * the runtime (Matlab or Octave) supports a 64 bit mwSize type. Even
 * then there are some limits imposed by available system memory and
 * processor architecture specific limits, e.g., many cpu's only
 * truly support 48 bit of memory, not the true 64 bit.
 *
 */
void PsychCheckmWSizeLimits(psych_int64 m, psych_int64 n, psych_int64 p)
{
	// No problem if mwSize has capacity to contain size_t:
	if (sizeof(mwSize) >= sizeof(size_t)) return;
	
	// Mismatch: mwSize too small to contain size_t. This means
	// mwSize is only defined as a 32-bit int on a system where
	// size_t is a 64 bit unsigned int. As long as the provided
	// size_t values are smaller than 2^31 we're good and don't
	// screw up when downcasting. Otherwis we'd overflow and
	// disaster would happen:
	if ((m < INT_MAX) && (n < INT_MAX) && (p < INT_MAX)) return;
	
	// Ok, this is a no-go :( Abort with some hopefully helpful
	// error message:
	printf("PTB-ERROR: Tried to return a vector or matrix whose size along at least one dimension\n");
	printf("PTB-ERROR: exceeds the maximum supported size of 2^31 - 1 elements.\n");
	if (sizeof(size_t) == 4) {
		printf("PTB-ERROR: This is a limitation of all 32 bit versions of Psychtoolbox.\n");
		printf("PTB-ERROR: You'd need to use a Psychtoolbox for 64-bit Matlab or 64-bit Octave\n");
		printf("PTB-ERROR: on a 64-bit operating system to get rid of this limit.\n");
	}
	else {
		printf("PTB-ERROR: This is a limitation of your version of Octave or Matlab.\n");
		printf("PTB-ERROR: You'd need to use a Psychtoolbox for 64-bit Matlab or 64-bit Octave\n");
		printf("PTB-ERROR: on a 64-bit operating system to get rid of this limit.\n");
	}
	PsychErrorExitMsg(PsychError_user, "One of the dimensions of a returned matrix or vector exceeds 2^31-1 elements. This is not supported on your setup!");
}

/*
	Get the third array dimension which we call "P".  mxGetP should act just like mxGetM and mxGetN.
        	
	The abstracted Psychtoolbox API supports matrices with up to 3 dimensions.     
*/
static mwSize mxGetP(const mxArray *arrayPtr)
{
	const mwSize *dimArray;
	
	if (mxGetNumberOfDimensions(arrayPtr) < 3) {
		return(1);
	}
	dimArray = (const mwSize*) mxGetDimensions(arrayPtr);
	return dimArray[2];
}


/*
	Get the 2nd array dimension.
        
	The Mex API's mxGetN is sometimes undersirable because it returns the product of all dimensions above 1.  Our mxGetNOnly only returns N, for when you need that.       
	
	The abstracted Psychtoolbox API supports matrices with up to 3 dimensions.     
*/
static mwSize mxGetNOnly(const mxArray *arrayPtr)
{
	const mwSize *dimArray;
	
	dimArray = (const mwSize*) mxGetDimensions(arrayPtr);
	return dimArray[1];
}


/*
    mxCreateDoubleMatrix3D()
    
    Create a 2D or 3D matrix of doubles. 
	
    Requirements are that m>0, n>0, p>=0.  
*/
mxArray *mxCreateDoubleMatrix3D(psych_int64 m, psych_int64 n, psych_int64 p)
{
	int numDims;
	mwSize dimArray[3];
	
	if(m==0 || n==0 ){
		dimArray[0]=0;dimArray[1]=0;dimArray[2]=0;	//this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices. 
	}else{
		PsychCheckmWSizeLimits(m,n,p);
		dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
	}
	numDims= (p==0 || p==1) ? 2 : 3;
	
	return mxCreateNumericArray(numDims, (mwSize*) dimArray, mxDOUBLE_CLASS, mxREAL);		
}

/*
    mxCreateFloatMatrix3D()
    
    Create a 2D or 3D matrix of floats. 
	
    Requirements are that m>0, n>0, p>=0.  
*/
mxArray *mxCreateFloatMatrix3D(size_t m, size_t n, size_t p)
{
	int numDims;
	mwSize dimArray[3];
	
	if(m==0 || n==0 ){
		dimArray[0]=0;dimArray[1]=0;dimArray[2]=0;	//this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices. 
	}else{
		PsychCheckmWSizeLimits(m,n,p);
		dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
	}
	numDims= (p==0 || p==1) ? 2 : 3;
	
	return mxCreateNumericArray(numDims, (mwSize*) dimArray, mxSINGLE_CLASS, mxREAL);		
}

/*
    mxCreateNativeBooleanMatrix3D()
    
    Create a 2D or 3D matrix of native psych_bool types. 
	
    Requirements are that m>0, n>0, p>=0.  
*/
mxArray *mxCreateNativeBooleanMatrix3D(size_t m, size_t n, size_t p)
{
	int			numDims;
	mwSize		dimArray[3];
	mxArray		*newArray;
	
	if(m==0 || n==0 ){
		dimArray[0]=0;dimArray[1]=0;dimArray[2]=0;	//this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices. 
	}else{
		PsychCheckmWSizeLimits(m,n,p);
		dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
	}
	numDims = (p==0 || p==1) ? 2 : 3;
	newArray = mxCreateNumericArray(numDims, (mwSize*) dimArray, mxLOGICAL_CLASS, mxREAL);
	
#if mxLOGICAL_CLASS == mxUINT8_CLASS
#if PSYCH_SYSTEM == PSYCH_LINUX
	// Manually set the flag to logical for Matlab versions < 6.5
	// MK DISABLED for Matlab > 6.5 mxSetLogical(newArray);
#endif
#endif
	
	return(newArray);		
}


/*
	Create a 2D or 3D matrix of ubytes.  
	
	Requirements are that m>0, n>0, p>=0.  
*/
mxArray *mxCreateByteMatrix3D(size_t m, size_t n, size_t p)
{
	int numDims;
	mwSize dimArray[3];
	
	if(m==0 || n==0 ){
		dimArray[0]=0;dimArray[1]=0;dimArray[2]=0; //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
	}else{
		PsychCheckmWSizeLimits(m,n,p);
		dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
	}
	numDims= (p==0 || p==1) ? 2 : 3;
	return mxCreateNumericArray(numDims, (mwSize*) dimArray, mxUINT8_CLASS, mxREAL);	
} 
 


/*
	Print string s and return return control to the calling environment.
*/
void PsychErrMsgTxt(char *s)
{
	// If this is the Matlab version then we try to close onscreen windows here.
	// The Octave version has a different implementation of this facility...
	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		// Is this the Screen() module?
		if (strcmp(PsychGetModuleName(), "Screen")==0) {
			// Yes. We directly call our close and cleanup routine:
			#ifdef PTBMODULE_Screen
			ScreenCloseAllWindows();
			#endif
		} else {
			// Nope. This is a Psychtoolbox MEX file other than Screen.
			// We can't call directly, but we can call the 'sca' command
			// from Matlab:
			PsychRuntimeEvaluateString("Screen('CloseAll');");
		}
	#endif
	
	// Call the Matlab- or Octave error printing and error handling facilities:
	mexErrMsgTxt((s && (strlen(s) > 0)) ? s : "See error message printed above.");
}


/*
	classify the mxArray element format using Pyschtoolbox argument type names
	
*/
static PsychArgFormatType PsychGetTypeFromMxPtr(const mxArray *mxPtr)
{
	PsychArgFormatType format;

	if(PsychIsDefaultMat(mxPtr))
		format = PsychArgType_default;
	else if(mxIsUint8(mxPtr))
		format = PsychArgType_uint8;
	else if(mxIsUint16(mxPtr))
		format = PsychArgType_uint16;
	else if(mxIsUint32(mxPtr))
		format = PsychArgType_uint32;
	else if(mxIsInt8(mxPtr))
		format = PsychArgType_int8;
	else if(mxIsInt16(mxPtr))
		format = PsychArgType_int16;
	else if(mxIsInt32(mxPtr))
		format = PsychArgType_int32;
	else if(mxIsDouble(mxPtr))
		format = PsychArgType_double;
	else if(mxIsSingle(mxPtr))
		format = PsychArgType_single;
	else if(mxIsChar(mxPtr))
		format = PsychArgType_char;
	else if(mxIsCell(mxPtr))
		format = PsychArgType_cellArray;
	else if(mxIsLogical(mxPtr))
		format = PsychArgType_boolean;  // This is tricky because MATLAB abstracts "logicals" conditionally on platform. Depending on OS, MATLAB implements booleans with either 8-bit or 64-bit values.  
	else 
		format = PsychArgType_unclassified;

	return format;	
			
}


/*
    PsychSetReceivedArgDescriptor()
    
    Accept an argument number and direction value (input or output).  Examine the specified argument and fill in an argument
    descriptor struture.  Ask a retainer function to store the descriptor. 
    
*/
PsychError PsychSetReceivedArgDescriptor(int argNum, psych_bool allow64BitSizes, PsychArgDirectionType direction)
{
	PsychArgDescriptorType d;
	int numNamedOutputs, numOutputs;
		
	const mxArray *mxPtr;

	d.position = argNum;
	d.direction = direction;	
	if(direction == PsychArgIn){
		mxPtr = PsychGetInArgMxPtr(argNum);
		d.isThere = (mxPtr && !PsychIsDefaultMat(mxPtr)) ? kPsychArgPresent : kPsychArgAbsent; 
		if(d.isThere == kPsychArgPresent){ //the argument is there so fill in the rest of the description
			d.numDims = (int) mxGetNumberOfDimensions(mxPtr);

			// If the calling function doesn't allow 64 bit sized input argument dimensions, then we check if
			// the input has a size within the positive signed integer range, i.e., at most INT_MAX elements
			// per dimension. Functions which can handle bigger inputs need to declare this explicitely by
			// setting allow64BitSizes == TRUE:
			if ((!allow64BitSizes) && ((mxGetM(mxPtr) >= INT_MAX) || (mxGetNOnly(mxPtr) >= INT_MAX) || (mxGetP(mxPtr) >= INT_MAX))) {
				printf("PTB-ERROR: %i. input argument exceeds allowable maximum size of 2^31 - 1 elements\n", argNum);
				printf("PTB-ERROR: in at least one dimension. Psychtoolbox can't handle such huge matrices or vectors.\n");
				PsychErrorExitMsg(PsychError_user, "Input argument exceeds maximum supported count of 2^31 - 1 elements!");
			}
			
			d.mDimMin = d.mDimMax = (psych_int64) mxGetM(mxPtr);
			d.nDimMin = d.nDimMax = (psych_int64) mxGetNOnly(mxPtr);
			d.pDimMin = d.pDimMax = (psych_int64) mxGetP(mxPtr);
			d.type = PsychGetTypeFromMxPtr(mxPtr);
		}
	}
	else{ //(direction == PsychArgOut)
		numNamedOutputs = PsychGetNumNamedOutputArgs();
		numOutputs = PsychGetNumOutputArgs();
		if(numNamedOutputs >=argNum)
			d.isThere = kPsychArgPresent;
		else if(numOutputs >=argNum)
			d.isThere = kPsychArgFixed;
		else
			d.isThere = kPsychArgAbsent;
	}
	PsychStoreArgDescriptor(NULL,&d);
	return(PsychError_none);	
							
}

PsychError PsychSetSpecifiedArgDescriptor(	int			position,
                                                        PsychArgDirectionType 	direction,
                                                        PsychArgFormatType 	type,
                                                        PsychArgRequirementType	isRequired,
                                                        psych_int64	mDimMin,		// minimum minimum is 1   |   
                                                        psych_int64	mDimMax, 		// minimum maximum is 1, maximum maximum is -1 meaning infinity
                                                        psych_int64	nDimMin,		// minimum minimum is 1   |   
                                                        psych_int64	nDimMax,		// minimum maximum is 1, maximum maximum is -1 meaning infinity
                                                        psych_int64	pDimMin,	    // minimum minimum is 0
                                                        psych_int64	pDimMax)		// minimum maximum is 0, maximum maximum is -1 meaning infinity
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
		PsychCheckmWSizeLimits((size_t) mDimMin, (size_t) nDimMin, (size_t) pDimMin);
		PsychCheckmWSizeLimits((size_t) mDimMax, (size_t) nDimMax, (size_t) pDimMax);
	}

	d.position = position;
	d.direction = direction;
	d.type = type;
	//d.isThere 			//field set only in the received are descriptor, not in the specified argument descriptor
	d.isRequired = isRequired;	//field set only in the specified arg descritor, not in the received argument descriptot.
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
	PsychAcceptInputArgumentDecider()
	
	This is a subroutine of Psychtoolbox functions such as PsychCopyInDoubleArg() which read in arguments to Psychtoolbox functino 
	passed from the scripting environment.  
	
	Accept one constant specifying whether an argument is either required, optional, or anything will be allowed and another constant
	specifying how the provided argument agrees with the specified argument.  Based on the relationship between those constants either:
	
		� Return TRUE indicating that the caller should read in the argument and itself return TRUE to indicate that the argument has been read.
		� Return FALSE indicating that the caller should ignore the argument and itself return FALSE to indicate that the argument was not read.
		� Exit to the calling environment with an error to indicate that the provided argument did not match the requested argument and that
		it was required to match.
		
	
	The domain of supplied arguments is: 
	
	matchError:
		PsychError_internal					-Internal Psychtoolbox error
		PsychError_invalidArg_absent		-There was no argument provided
		PsychError_invalidArg_type			-The argument was present but not the specified type
		PsychError_invalidArg_size			-The argument was presnet and the specified type but not the specified size
		PsychError_none						-The argument matched the specified argument
		
	isRequired:
		kPsychArgRequired					- the argument must be present and must match the specified descriptor
		kPsychArgOptional					- the argument must either be absent or must be present and match the specified descriptor
		kPsychArgAnything					- the argument can be absent or anything

*/
psych_bool PsychAcceptInputArgumentDecider(PsychArgRequirementType isRequired, PsychError matchError)
{
	if(isRequired==kPsychArgRequired){
		if(matchError)
			PsychErrorExit(matchError); 
		else
			return(TRUE);
	}else if(isRequired==kPsychArgOptional){
		if(matchError==PsychError_invalidArg_absent)
			return(FALSE);
		else if(matchError)
			PsychErrorExit(matchError);
		else 
			return(TRUE);
	}else if(isRequired==kPsychArgAnything){
		if(!matchError)
			return(TRUE);
		else if(matchError==PsychError_invalidArg_absent)
			return(FALSE);
		else if(matchError==PsychError_invalidArg_type)
			return(FALSE);
		else if(matchError==PsychError_invalidArg_size)
			return(FALSE);
		else
			PsychErrorExit(matchError);
	}
	PsychErrorExitMsg(PsychError_internal, "Reached end of function unexpectedly");
	return(FALSE);			//make the compiler happy
}



/*

	PsychAcceptOutputArgumentDecider()
	
	This is a subroutine of Psychtoolbox functions such as PsychCopyCopyDoubleArg() which output arguments from Psychtoolbox functions 
	back to the scripting environment.  

*/
psych_bool PsychAcceptOutputArgumentDecider(PsychArgRequirementType isRequired, PsychError matchError)
{

	if(isRequired==kPsychArgRequired){
		if(matchError)
			PsychErrorExit(matchError);							//the argument was required and absent so exit with an error. Or there was some other error.
		else 
			return(TRUE);										//the argument was required and present so go read it. 
	}else if(isRequired==kPsychArgOptional){
		if(!matchError)
			return(TRUE);										//the argument was optional and present so go read it.  
		else if(matchError==PsychError_invalidArg_absent)
			return(FALSE);										//the argument as optional and absent so dont' read  it. 
		else  if(matchError)
			PsychErrorExit(matchError);							//there was some other error
	}else if(isRequired==kPsychArgAnything) 
		PsychErrorExitMsg(PsychError_internal, "kPsychArgAnything argument passed to an output function.  Use kPsychArgOptional");
	else
		PsychErrorExit(PsychError_internal);
	
	PsychErrorExitMsg(PsychError_internal, "End of function reached unexpectedly");
	return(FALSE);		//make the compiler happy
}


/*
    PsychMatchDescriptors()
    
	Compare descriptors for specified and received arguments. Return a mismatch error if they are 
	incompatible, otherwise return a no error.
	
	PsychMatchDescriptors compares:
		The argument type
		The argument size
		Argument presense 
	
	PsychMatchDescripts can return any of the following values describing the relationship between an
	argument provided from the scripting environment and argument requested by a Psychtoolbox module:
		PsychError_internal					-Internal Psychtoolbox error
		PsychError_invalidArg_absent		-There was no argument provided
		PsychError_invalidArg_type			-The argument was present but not the specified type
		PsychError_invalidArg_size			-The argument was presnet and the specified type but not the specified size
		PsychError_none						-The argument matched the specified argument
		
    This function should be enhnaced to report the nature of the disagrement
*/


PsychError PsychMatchDescriptors(void)
{
	PsychArgDescriptorType *specified, *received;

	PsychGetArgDescriptor(&specified, &received);
	
	//check for various bogus conditions resulting only from Psychtoolbox bugs and issue an internal error. 
	if(specified->position != received->position)
		PsychErrorExit(PsychError_internal);
	if(specified->direction != received->direction)
		PsychErrorExit(PsychError_internal);
	
	if(specified->direction==PsychArgOut) {
		if(received->isThere==kPsychArgPresent || received->isThere==kPsychArgFixed)
			return(PsychError_none);
		else
			return(PsychError_invalidArg_absent);
	}
	if(specified->direction==PsychArgIn){
		if(received->isThere==kPsychArgAbsent)  
			return(PsychError_invalidArg_absent);
		//otherwise the argument is present and we proceed to the argument type and size checking block below 
	}

	//if we get to here it means that an input argument was supplied.  Check if it agrees in type and size with the specified arg and return 
	// an error type accordingly
	if(!(specified->type & received->type))
		return(PsychError_invalidArg_type);
	if(received->mDimMin != received->mDimMax || received->nDimMin != received->nDimMax ||  received->pDimMin != received->pDimMax)  
		PsychErrorExit(PsychError_internal);	//unnecessary mandate  
	if(received->mDimMin < specified->mDimMin)
		return(PsychError_invalidArg_size);
	if(received->nDimMin < specified->nDimMin)
		return(PsychError_invalidArg_size);
	if(specified->pDimMin != kPsychUnusedArrayDimension && received->pDimMin < specified->pDimMin)
		return(PsychError_invalidArg_size);
	if(specified->mDimMax != kPsychUnboundedArraySize && received->mDimMax > specified->mDimMax) 
		return(PsychError_invalidArg_size);
	if(specified->nDimMax != kPsychUnboundedArraySize && received->nDimMax > specified->nDimMax) 
		return(PsychError_invalidArg_size);
	if(specified->pDimMax != kPsychUnusedArrayDimension && specified->pDimMax != kPsychUnboundedArraySize && received->pDimMax > specified->pDimMax) 
		return(PsychError_invalidArg_size);
	if(received->numDims > 3)  //we don't allow matrices with more than 3 dimensions.
		return(PsychError_invalidArg_size);

	//if we get to here it means that  the block above it means 
	return(PsychError_none);	
}




PsychError PsychMatchDescriptorsOld(void)
{
	PsychArgDescriptorType *specified, *received;

	PsychGetArgDescriptor(&specified, &received);
	
	//check for various bogus conditions resulting only from Psychtoolbox bugs and issue an internal error
	if(specified->position != received->position)
		PsychErrorExit(PsychError_internal);
	if(specified->direction != received->direction)
		PsychErrorExit(PsychError_internal);
	
	switch(specified->direction) {
		case PsychArgOut:  
			switch(specified->isRequired){
				case kPsychArgRequired:
					switch(received->isThere){
						case kPsychArgPresent:		
							goto exitOk;					//both the argument is present and the return variable is named within the calling environment.
						case kPsychArgFixed:		    
							goto exitOk;					//the argument is present but a return variable is not named within the calling environment. In MATLAB this can only be the 1st return argument.
						case kPsychArgAbsent:				
							return(PsychError_invalidArg_absent);	//neither a return argument is present nor a return variable is named within the calling environment.
					}
				case kPsychArgOptional:	case kPsychArgAnything:			
					switch(received->isThere){
						case kPsychArgPresent:
							goto exitOk;
						case kPsychArgFixed:
							goto exitOk;
						case kPsychArgAbsent:
							goto exitOk;
					}
			}
			break;
		case PsychArgIn:
			switch(specified->isRequired){
				case kPsychArgRequired:
					switch(received->isThere){
						case kPsychArgPresent:
							break;						//we still need to comppare the actual type and size to specifications.
						case kPsychArgFixed:  
							PsychErrorExitMsg(PsychError_internal,"The input argument descriptor specifies a fixed argument, this property is unallowed for inputs.");
						case kPsychArgAbsent:
							return(PsychError_invalidArg_absent);
					}
					break;
				case kPsychArgOptional: case kPsychArgAnything: 
					switch(received->isThere){
						case kPsychArgPresent:
							break;						//we still need to comppare the actual type and size to specifications.
						case kPsychArgFixed:
							PsychErrorExitMsg(PsychError_internal,"The input argument descriptor secifies a fixed argument, this property is unallowed for inputs.");
			 			case kPsychArgAbsent:
							goto exitOk;				//we do NOT need to compare the actual type and size to specifications. 
					}
					break;

			}
			//if we get to here we are assured that both an input argument was specified and there is one there. In this block we compare they type and size of
			//specified and provide arguments.  For output arguments we skip over this block because those are not assigned types by the calling environment.  
			if(!(specified->type & received->type))
				return(PsychError_invalidArg_type);
			if(received->mDimMin != received->mDimMax || received->nDimMin != received->nDimMax ||  received->pDimMin != received->pDimMax)  
				PsychErrorExit(PsychError_internal);	//unnecessary mandate  
			if(received->mDimMin < specified->mDimMin)
				return(PsychError_invalidArg_size);
			if(received->nDimMin < specified->nDimMin)
				return(PsychError_invalidArg_size);
			if(received->pDimMin < specified->pDimMin)
				return(PsychError_invalidArg_size);
			if(specified->mDimMax != kPsychUnboundedArraySize && received->mDimMax > specified->mDimMax) 
				return(PsychError_invalidArg_size);
			if(specified->nDimMax != kPsychUnboundedArraySize && received->nDimMax > specified->nDimMax) 
				return(PsychError_invalidArg_size);
			if(specified->pDimMax != kPsychUnboundedArraySize && received->pDimMax > specified->pDimMax) 
				return(PsychError_invalidArg_size);
			if(received->numDims > 3)  //we don't allow matrices with more than 3 dimensions.
				return(PsychError_invalidArg_size);
			break; 
	}
		
	exitOk: 
	return(PsychError_none);	
}






//local function definitions for ScriptingGlue.cpp
//___________________________________________________________________________________________


psych_bool PsychIsDefaultMat(CONSTmxArray *mat)
{
	return (PsychIsEmptyMat(mat) && !mxIsChar(mat));
}

psych_bool PsychIsEmptyMat(CONSTmxArray *mat)
{
	return(mxGetM(mat)==0 || mxGetN(mat)==0);
}


//functions for project access to module call arguments (MATLAB)
//___________________________________________________________________________________________


//functions which query the number and nature of supplied arguments

/* 
	PsychGetNumInputArgs()

	-The count excludes the command argument and includes ALL arguments supplied, including 
	default arguments.  
	
	-For the time being, the only way to check if all required arguments are supplied in the 
	general case of mixed required and optional arguments is to check each individually. Probably
	the best way to to fix this is to employ a description of which are required and which optional
	and compare that against what was passed to the subfunction.
*/
int PsychGetNumInputArgs(void)
{
	if(PsychAreSubfunctionsEnabled() && !baseFunctionInvoked) //this should probably be just baseFunctionInvoked wo PsychSubfunctionEnabled.
		return(nrhsGLUE-1);
	else
		return(nrhsGLUE);
}

int PsychGetNumOutputArgs(void)
{
	return(nlhsGLUE==0 ? 1 : nlhsGLUE);
}

int PsychGetNumNamedOutputArgs(void)
{
	return(nlhsGLUE);
}

PsychError PsychCapNumInputArgs(int maxInputs)
{
    if(PsychGetNumInputArgs() > maxInputs)
            return(PsychError_extraInputArg);
    else
            return(PsychError_none);
}

PsychError PsychRequireNumInputArgs(int minInputs)
{
    if(PsychGetNumInputArgs() < minInputs)
            return(PsychError_missingInputArg);
    else
            return(PsychError_none);
    
}

PsychError PsychCapNumOutputArgs(int maxNamedOutputs)
{
	if(PsychGetNumNamedOutputArgs() > maxNamedOutputs)
		return(PsychError_extraOutputArg);
	else
		return(PsychError_none);
}



/*
	The argument is not present if a default m*n=0 matrix was supplied, '' or []	
*/
psych_bool PsychIsArgPresent(PsychArgDirectionType direction, int position)
{
	int numArgs;
	
	if(direction==PsychArgOut){
		return((psych_bool)(PsychGetNumOutputArgs()>=position));
	}else{
		if((numArgs=PsychGetNumInputArgs())>=position)
			return(!(PsychIsDefaultMat(PsychGetInArgMxPtr(position)))); //check if its default
		else
			return(FALSE);
	}
}

/*
	The argument is present if anything was supplied, including the default matrix
*/
psych_bool PsychIsArgReallyPresent(PsychArgDirectionType direction, int position)
{
	
	return(direction==PsychArgOut ? PsychGetNumOutputArgs()>=position : PsychGetNumInputArgs()>=position);
}



PsychArgFormatType PsychGetArgType(int position) //this is for inputs because outputs are unspecified by the calling environment.
{
	if(!(PsychIsArgReallyPresent(PsychArgIn, position)))
		return(PsychArgType_none);
	
	return(PsychGetTypeFromMxPtr(PsychGetInArgMxPtr(position)));	
}

size_t PsychGetArgM(int position)
{
	if(!(PsychIsArgPresent(PsychArgIn, position)))
		PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
	return( mxGetM(PsychGetInArgMxPtr(position)));
}

size_t PsychGetArgN(int position)
{
	if(!(PsychIsArgPresent(PsychArgIn, position)))
		PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
	return( mxGetNOnly(PsychGetInArgMxPtr(position)));
}


size_t PsychGetArgP(int position)
{
	if(!(PsychIsArgPresent(PsychArgIn, position)))
		PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
	return( mxGetP(PsychGetInArgMxPtr(position)));
}


/*
    PyschCheckInputArgType()
    
    Check that the input argument at the specifid position matches at least one of the types passed in the argType
    argument.  If the argument violates the proscription exit with an error.  Otherwise return a psych_bool indicating
    whether the argument was present.   
    
*/
psych_bool PsychCheckInputArgType(int position, PsychArgRequirementType isRequired, PsychArgFormatType argType)
{
	PsychError		matchError;
	psych_bool			acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, argType, isRequired, 0,kPsychUnboundedArraySize,0,kPsychUnboundedArraySize,0,kPsychUnboundedArraySize);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    return(acceptArg);
}



 
/*functions which output arguments.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
-Naming scheme:
	-Outputing return arguments:
		- "PsychAllocOut*Arg" : allocate and set a pointer to volatile memory to be filled with returned information by the caller.
		- "PsychCopyOut*Arg : accept a pointer to ouput values and fill in the return matrix memory with the values.
	-Reading input arguments:
		- "PsychAllocIn*Arg" : set a pointer to volatile memory allocated by "PsychAllocIn*Arg" and holding the input value.
		- "PsychCopyIn*Arg" : accept a pointer to memory and fill in that memory with the input argument values.     

-These all supply their own dynamic memory now, even functions which return arguments, and, in the case of
 Put functions,  even when those arguments are not present !  If you don't want the function to go allocating
 memory for an unsupplied return argument, detect the presense of that argument from within your script and
 conditionally invoke PsychPut*Arg.  This is a feature which allows you to ignore the presense of a return 
 argument in the case where memory which holds a return argument serves other purposes.   
 
-All dynamic memory provided by these functions is volatile, that is, it is lost when the mex module exits and
returns control to the Matlab environemnt.  To make it non volatile, call Psych??? on it.  

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/



psych_bool PsychCopyOutDoubleArg(int position, PsychArgRequirementType isRequired, double value)
{
	mxArray **mxpp;
	PsychError matchError;
	psych_bool putOut;
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double,  isRequired, 1,1,1,1,0,0);
	matchError=PsychMatchDescriptors();
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		(*mxpp) = mxCreateDoubleMatrix(1,1,mxREAL);
		mxGetPr((*mxpp))[0] = value;
	}
	return(putOut);
}



/*
    PsychAllocOutDoubleArg_2()
    
    usage:
    psych_bool PsychAllocOutDoubleArg_2(int position, PsychArgRequirementType isRequired, double **value)
    psych_bool PsychAllocOutDoubleArg_2(int position, PsychArgRequirementType isRequired, double **value, PsychGenericScriptType **nativeDouble)
    
    PsychAllocOutDoubleArg_2() is an experimental enhanced version of PsychAllocOutDoubleArg which will accept the kPsychNoArgReturn  
    constant in the position argument and then return via the optional 4th input a pointer to a native scripting type which holds the 
    double.
    
    Having a reference to the native type allows us to embed doubles withing cell arrays and structs and to pass doubles as arguments to functions
    called within MATLAB from a mex file.
    
    PsychAllocOutDoubleArg_2() should be backwards compatible with PsychAllocOutDoubleArg and could supplant that function.      
    
*/
/*
psych_bool PsychAllocOutDoubleArg_2(int position, PsychArgRequirementType isRequired, double **value, ...)
{
	mxArray **mxpp;
        va_list ap;
        
        if(position != kPsychNoArgReturn){
            PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
            PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, 1,1,1,1,0,0);
            PsychErrorExit(PsychMatchDescriptors());    
            mxpp = PsychGetOutArgMxPtr(position);
            if(mxpp == NULL){  //Here we allocated memory even if the return argument is not present.  Controversial.  
                    *value= (double *)mxMalloc(sizeof(double));
                    return(FALSE); 
            }
            else{
                    *mxpp = mxCreateDoubleMatrix3D(1,1,0);
                    *value = mxGetPrPtr(*mxpp);
                    return(TRUE);   
            }
        }else{
            va_start(ap, value);
            *(mxArray**)ap=mxCreateDoubleMatrix3D(1,1,0);
            *value = mxGetPrPtr(*(mxArray**)ap);
            va_end(ap);
            return(TRUE);
        }
}  
*/


  
psych_bool PsychAllocOutDoubleArg(int position, PsychArgRequirementType isRequired, double **value)
{
	mxArray			**mxpp;
	PsychError		matchError;
	psych_bool			putOut;
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, 1,1,1,1,0,0);
	matchError=PsychMatchDescriptors();
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		*mxpp = mxCreateDoubleMatrix3D(1,1,0);
		*value = mxGetPr(*mxpp);
	}else{
		mxpp = PsychGetOutArgMxPtr(position);
		*value= (double *)mxMalloc(sizeof(double));
	}
	return(putOut);
}


/* 
PsychAllocOutDoubleMatArg()

A)return argument mandatory:
	1)return argument not present: 		exit with an error.
	2)return argument present: 		allocate an output matrix and set return arg pointer. Set *array to the array within the new matrix. Return TRUE.  
B)return argument optional:
	1)return argument not present:  	return FALSE to indicate absent return argument.  Create an array.   Set *array to the new array. 
	2)return argument present:	 	allocate an output matrix and set return arg. pointer. Set *array to the array within the new matrix.  Return TRUE.   
*/
psych_bool PsychAllocOutDoubleMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, double **array)
{
	mxArray			**mxpp;
	PsychError		matchError;
	psych_bool		putOut;
	
	PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, m,m,n,n,p,p);
	matchError=PsychMatchDescriptors();
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		*mxpp = mxCreateDoubleMatrix3D(m,n,p);
		*array = mxGetPr(*mxpp);
	}else
		*array= (double *) mxMalloc(sizeof(double) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
	return(putOut);
}



/* 
PsychAllocOutFloatMatArg()

This function allocates out a matrix of single precision floating point type,
that is C data type 32-bit float or Matlab/Octave data type single().

A)return argument mandatory:
	1)return argument not present: 		exit with an error.
	2)return argument present: 		allocate an output matrix and set return arg pointer. Set *array to the array within the new matrix. Return TRUE.  
B)return argument optional:
	1)return argument not present:  	return FALSE to indicate absent return argument.  Create an array.   Set *array to the new array. 
	2)return argument present:	 	allocate an output matrix and set return arg. pointer. Set *array to the array within the new matrix.  Return TRUE.   
*/
psych_bool PsychAllocOutFloatMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, float **array)
{
	mxArray			**mxpp;
	PsychError		matchError;
	psych_bool		putOut;
	
	PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_single, isRequired, m,m,n,n,p,p);
	matchError=PsychMatchDescriptors();
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		*mxpp = mxCreateFloatMatrix3D((size_t) m, (size_t) n, (size_t) p);
		*array = (float*) mxGetData(*mxpp);
	}else
		*array = (float*) mxMalloc(sizeof(float) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
	return(putOut);
}



/*
    PsychCopyOutBooleanArg()
*/
psych_bool PsychCopyOutBooleanArg(int position, PsychArgRequirementType isRequired, PsychNativeBooleanType value)
{
	mxArray			**mxpp;
	PsychError		matchError;
	psych_bool			putOut;
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_boolean, isRequired, 1,1,1,1,0,0);
	matchError=PsychMatchDescriptors();
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		(*mxpp) = mxCreateLogicalMatrix(1,1);
		mxGetLogicals((*mxpp))[0] = value;
	}
	return(putOut);
}



/*
    PsychAllocOutBooleanArg()
*/
psych_bool PsychAllocOutBooleanArg(int position, PsychArgRequirementType isRequired, PsychNativeBooleanType **value)
{
	mxArray **mxpp;
	PsychError		matchError;
	psych_bool			putOut;
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_boolean, isRequired, 1,1,1,1,0,0);
	matchError=PsychMatchDescriptors(); 
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		(*mxpp) = mxCreateLogicalMatrix(1,1);
		*value = mxGetLogicals((*mxpp));
	}else{
		mxpp = PsychGetOutArgMxPtr(position);
		*value= (PsychNativeBooleanType *)mxMalloc(sizeof(PsychNativeBooleanType));
	}
	return(putOut);
}    


/* 
    PsychAllocOutBooleanMatArg()

    A)return argument mandatory:
	1)return argument not present: 		exit with an error.
	2)return argument present: 		allocate an output matrix and set return arg pointer. Set *array to the array within the new matrix. Return TRUE.  
    B)return argument optional:
	1)return argument not present:  	return FALSE to indicate absent return argument.  Create an array.   Set *array to the new array. 
	2)return argument present:	 	allocate an output matrix and set return arg. pointer. Set *array to the array within the new matrix.  Return TRUE.   
*/
psych_bool PsychAllocOutBooleanMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, PsychNativeBooleanType **array)
{
	mxArray			**mxpp;
	PsychError		matchError;
	psych_bool		putOut;
	
	PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_boolean, isRequired, m,m,n,n,p,p);
	matchError=PsychMatchDescriptors(); 
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		*mxpp = mxCreateNativeBooleanMatrix3D((size_t) m, (size_t) n, (size_t) p);
		*array = (PsychNativeBooleanType *)mxGetLogicals(*mxpp);
	}else{
		*array= (PsychNativeBooleanType *) mxMalloc(sizeof(PsychNativeBooleanType) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
	}
	return(putOut);
}




/* 
    PsychAllocOutUnsignedByteMatArg()
    
    Like PsychAllocOutDoubleMatArg() execept for unsigned bytes instead of doubles.  
*/
psych_bool PsychAllocOutUnsignedByteMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, psych_uint8 **array)
{
	mxArray			**mxpp;
	PsychError		matchError;
	psych_bool		putOut;
	
	PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_uint8, isRequired, m,m,n,n,p,p);
	matchError=PsychMatchDescriptors(); 
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		*mxpp = mxCreateByteMatrix3D((size_t) m, (size_t) n, (size_t) p);
		*array = (psych_uint8 *)mxGetData(*mxpp);
	}else{
		*array= (psych_uint8 *) mxMalloc(sizeof(psych_uint8) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
	}
	return(putOut);
}



psych_bool PsychCopyOutDoubleMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, double *fromArray)
{
	mxArray **mxpp;
	double *toArray;
	PsychError		matchError;
	psych_bool		putOut;
	
	PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, m,m,n,n,p,p);
	matchError=PsychMatchDescriptors(); 
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		*mxpp = mxCreateDoubleMatrix3D(m,n,p);
		toArray = mxGetPr(*mxpp);
		//copy the input array to the output array now
		memcpy(toArray, fromArray, sizeof(double) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
	}
	return(putOut);
}

psych_bool PsychCopyOutUnsignedInt16MatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, psych_uint16 *fromArray)
{
	mxArray **mxpp;
	psych_uint16 *toArray;
	PsychError		matchError;
	psych_bool		putOut;
	mwSize dimArray[3];
	int numDims;

	// Compute output array dimensions:
	if (m<=0 || n<=0) {
		dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;	//this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices. 
	} else {
		PsychCheckmWSizeLimits(m,n,p);
		dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
	}
	numDims = (p == 0 || p == 1) ? 2 : 3;
	
	PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_uint16, isRequired, m,m,n,n,p,p);
	matchError=PsychMatchDescriptors();
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		*mxpp = mxCreateNumericArray(numDims, (mwSize*) dimArray, mxUINT16_CLASS, mxREAL);
		toArray = (psych_uint16*) mxGetData(*mxpp);
        
		//copy the input array to the output array now
		memcpy(toArray, fromArray, sizeof(psych_uint16) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
	}
	return(putOut);
}

/*
	PsychCopyOutCharArg()

	Accept a null terminated string and return it in the specified position.  
	  
*/
psych_bool PsychCopyOutCharArg(int position, PsychArgRequirementType isRequired, const char *str)
{
	mxArray **mxpp;
	PsychError		matchError;
	psych_bool			putOut;	

	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_char, isRequired, 0, strlen(str),0,strlen(str),0,0);
	matchError=PsychMatchDescriptors(); 
	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		*mxpp = mxCreateString(str);
	}
	return(putOut);
}




/*functions which input arguments.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/



/*
A)input argument mandatory:
 
	1)input argument not present: 		exit with error.
	2)input argument present: 			set *array to the input matrix, *m, *n, and *p to its dimensions, return TRUE.    
B)input argument optional:

	1)input argument not present: 		return FALSE
	2)input argument present: 			set *array to the input matrix, *m, *n, and *p to its dimensions, return TRUE.    

*/
// TO DO: Needs to be updated for kPsychArgAnything
psych_bool PsychAllocInDoubleMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, double **array)
{
    const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool			acceptArg;
    
    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,-1,1,-1,0,-1);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		*m = (int) mxGetM(mxPtr);
		*n = (int) mxGetNOnly(mxPtr);
		*p = (int) mxGetP(mxPtr);
		*array=mxGetPr(mxPtr);
	}
	return(acceptArg);
}

/* Alloc-in double matrix, but allow for 64-bit dimension specs. */
psych_bool PsychAllocInDoubleMatArg64(int position, PsychArgRequirementType isRequired, psych_int64 *m, psych_int64 *n, psych_int64 *p, double **array)
{
    const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool		acceptArg;
    
    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,-1,1,-1,0,-1);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		*m = (psych_int64) mxGetM(mxPtr);
		*n = (psych_int64) mxGetNOnly(mxPtr);
		*p = (psych_int64) mxGetP(mxPtr);
		*array=mxGetPr(mxPtr);
	}
	return(acceptArg);
}

/* Like PsychAllocInFloatMatArg64, but with 32-Bit int type size return-arguments. */
psych_bool PsychAllocInFloatMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, float **array)
{
    psych_int64 mb, nb, pb;
    psych_bool rc = PsychAllocInFloatMatArg64(position, isRequired, &mb, &nb, &pb, array);
    *m = (int) mb;
    *n = (int) nb;
    *p = (int) pb;
    return(rc);
}

/*

Alloc-In a single precision floating point matrix, i.e. a matrix of
C data type 32 bit float, aka Matlab/Octave data type single().
This function allows to alloc in matrices with more than 2^32 elements
per matrix dimension on 64 bit systems. Therefore the returned size
descriptors must be psych_int64 variables, not int variables or bad things
will happen.

If the function receives a double() precision input matrix instead of the
expected single() precision matrix, it will automatically create a temporary
copy, with all values copied/converted from double to single data type, aka
double -> float cast. This is transparent to the caller, so it can always
operate on a returned float matrix -- at a performance penalty for the extra
copy of course.

A)input argument mandatory:
 
	1)input argument not present: 		exit with error.
	2)input argument present: 			set *array to the input matrix, *m, *n, and *p to its dimensions, return TRUE.    
B)input argument optional:

	1)input argument not present: 		return FALSE
	2)input argument present: 			set *array to the input matrix, *m, *n, and *p to its dimensions, return TRUE.    

*/
psych_bool PsychAllocInFloatMatArg64(int position, PsychArgRequirementType isRequired, psych_int64 *m, psych_int64 *n, psych_int64 *p, float **array)
{
    const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool		acceptArg;
    double*         arrayD;
    float*          arrayF;
    psych_int64     i;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_single, isRequired, 1,-1,1,-1,0,-1);
	matchError=PsychMatchDescriptors();

    // Argument provided, but not of required float type?
    if (matchError == PsychError_invalidArg_type) {
        // If the input type is double precision floating point, then we convert
        // it here into single precision floating point via a temporary buffer.
        // This is used for functions which absolutely need float input, e.g.,
        // OpenGL-ES rendering code, but should accept double input from usercode
        // so usercode doesn not need to be specifically ported for OpenGL-ES platforms.
        // Performance may suffer somwehat though...
        PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,-1,1,-1,0,-1);
        matchError=PsychMatchDescriptors();
        acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
        if(acceptArg){
            mxPtr = PsychGetInArgMxPtr(position);
            *m = (psych_int64) mxGetM(mxPtr);
            *n = (psych_int64) mxGetNOnly(mxPtr);
            *p = (psych_int64) mxGetP(mxPtr);

            // Get a double pointer to the double input data matrix:
            arrayD = (double*) mxGetData(mxPtr);

            // Allocate temporary float input matrix. It will get deallocated
            // automatically at return to runtime:
            *array = (float*) PsychMallocTemp(sizeof(float) * (*m) * (*n) * (*p));
            arrayF = *array;

            // Copy/Convert loop:
            for (i = (*m) * (*n) * (*p); i > 0; i--) *(arrayF++) = (float) *(arrayD++);            
        }

        return(acceptArg);
    }

    // Regular path: Matching float (aka single()) matrix/vector provided:
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		*m = (psych_int64) mxGetM(mxPtr);
		*n = (psych_int64) mxGetNOnly(mxPtr);
		*p = (psych_int64) mxGetP(mxPtr);
		*array = (float*) mxGetData(mxPtr);
	}
	return(acceptArg);
}



/*
	PsychAllocInIntegerListArg()
	
	In a scriptiong language such as MATLAB where numbers are almost always stored as doubles, this function is useful to check
	that the value input is an integer value stored within a double type.
	
	Otherwise it just here to imitate the version written for other scripting languages.
*/
psych_bool PsychAllocInIntegerListArg(int position, PsychArgRequirementType isRequired, int *numElements, int **array)
{
    int m, n, p, i; 
    double *doubleMatrix;
    psych_bool isThere; 

    isThere=PsychAllocInDoubleMatArg(position, isRequired, &m, &n, &p, &doubleMatrix);
    if(!isThere)
        return(FALSE);
    p= (p==0) ? 1 : p;
	
	if ((psych_uint64) m * (psych_uint64) n * (psych_uint64) p >= INT_MAX) {
		printf("PTB-ERROR: %i th input argument has more than 2^31 - 1 elements! This is not supported.\n", position);
		return(FALSE);
	}
	
    *numElements = m * n * p;    				
    *array = (int*) mxMalloc((size_t) *numElements * sizeof(int));
    for(i=0; i < *numElements; i++){
        if(!PsychIsIntegerInDouble(doubleMatrix+i))
            PsychErrorExit(PsychError_invalidIntegerArg);
        (*array)[i] = (int) doubleMatrix[i];
    }
    return(TRUE);
}



/*
    PsychAllocInByteMatArg()
    
    Like PsychAllocInDoubleMatArg() except it returns an array of unsigned bytes.  
*/
psych_bool PsychAllocInUnsignedByteMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, unsigned char **array)
{
	const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool			acceptArg;

	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_uint8, isRequired, 1,-1,1,-1,0,-1);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		*m = (int) mxGetM(mxPtr);
		*n = (int) mxGetNOnly(mxPtr);
		*p = (int) mxGetP(mxPtr);
		*array=(unsigned char *)mxGetData(mxPtr);
	}
	return(acceptArg);
}

			 



/* 
	PsychCopyInDoubleArg()
	
	For 1x1 double.
 
	Return in *value a double passed in the specified position, or signal an error if there is no 
	double there and the argument is required, or don't change "value" if the argument is optional
	and none is supplied.  
	
    Note that if the argument is optional and ommitted PsychGetDoubleArg won't overwrite *value, allowing 
    for specification of default values within project files without checking for their
    presense and conditinally filing in values.  
*/
// TO DO: Needs to be updated for kPsychArgAnything
psych_bool PsychCopyInDoubleArg(int position, PsychArgRequirementType isRequired, double *value)
{
	const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool			acceptArg;
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,1,1,1,1,1);
	matchError=PsychMatchDescriptors();

	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		*value=mxGetPr(mxPtr)[0]; 
	}

	return(acceptArg); 
}


/*  
    Like PsychCopyInDoubleArg() with the additional restriction that the passed value not have a fractoinal componenet
    and that the it fit within thebounds of a C integer
    
    We could also accept matlab native integer types by specifying a conjunction of those as the third argument 
    in the PsychSetSpecifiedArgDescriptor() call, but why bother ?    
*/
psych_bool PsychCopyInIntegerArg(int position,  PsychArgRequirementType isRequired, int *value)
{
	const mxArray 	*mxPtr;
	double			tempDouble;
	PsychError		matchError;
	psych_bool			acceptArg;

	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,1,1,1,1,1);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		tempDouble=mxGetPr(mxPtr)[0];
		if(!PsychIsIntegerInDouble(&tempDouble))
			PsychErrorExit(PsychError_invalidIntegerArg);
		*value=(int)tempDouble;
	}
	return(acceptArg);
}

psych_bool PsychCopyInIntegerArg64(int position,  PsychArgRequirementType isRequired, psych_int64 *value)
{
	const mxArray 	*mxPtr;
	double			tempDouble;
	PsychError		matchError;
	psych_bool		acceptArg;

	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,1,1,1,1,1);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		tempDouble=mxGetPr(mxPtr)[0];
		if(!PsychIsInteger64InDouble(&tempDouble))
			PsychErrorExit(PsychError_invalidIntegerArg);
		*value= (psych_int64) tempDouble;
	}
	return(acceptArg);
}


/*
    PsychAllocInDoubleArg()
     
*/
psych_bool PsychAllocInDoubleArg(int position, PsychArgRequirementType isRequired, double **value)
{
	const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool			acceptArg;
	
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,1,1,1,1,1);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		*value=mxGetPr(mxPtr);
	}
	return(acceptArg);
}



/*
	PsychAllocInCharArg()
	
	Reads in a string and sets *str to point to the string.
	
	This function violates the rule for AllocIn fuctions that if the argument is optional and absent we allocate 
	space.  That turns out to be an unuseful feature anyway, so we should probably get ride of it.

	The second argument has been modified to passively accept, without error, an argument in the specified position of non-character type.  
          
        0	kPsychArgOptional  Permit either an argument of the specified type or no argument.  An argument of any a different type is an error.
        1	kPsychArgRequired  Permit only an argument of the specifed type.  If no argument is present, exit with error.
        2	kPsychArgAnything  Permit any argument type without error, but only read the specified type. 
		
*/
psych_bool PsychAllocInCharArg(int position, PsychArgRequirementType isRequired, char **str)
{
	const mxArray 	*mxPtr;
	int				status;
	psych_uint64	strLen;	
	PsychError		matchError;
	psych_bool		acceptArg;

	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_char, isRequired, 0, kPsychUnboundedArraySize ,0, kPsychUnboundedArraySize, 0 , 1);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr  = PsychGetInArgMxPtr(position);
		strLen = ((psych_uint64) mxGetM(mxPtr) * (psych_uint64) mxGetNOnly(mxPtr) * (psych_uint64) sizeof(mxChar)) + 1;
		if (strLen >= INT_MAX) PsychErrorExitMsg(PsychError_user, "Tried to pass in a string with more than 2^31 - 1 characters. Unsupported!");
		*str   = (char *) PsychCallocTemp((size_t) strLen, sizeof(char));
		status = mxGetString(mxPtr, *str, (mwSize) strLen); 
		if(status!=0)
			PsychErrorExitMsg(PsychError_internal, "mxGetString failed to get the string");
	}
	return(acceptArg);
}



/*
	Get a psych_bool flag from the specified argument position.  The matlab type can be be psych_bool, uint8, or
	char.  If the numerical value is equal to zero or if its empty then the flag is FALSE, otherwise the
	flag is TRUE.
	
	PsychGetFlagArg returns TRUE if the argument was present and false otherwise:
	
	A- Argument required
		1- Argument is present: load *argVal and return TRUE 
		2- Argument is absent: exit with an error
	B- Argument is optional
		1- Argument is present: load *argVal and return TRUE 
		2- Argument is absent: leave *argVal alone and return FALSE

	Note: if we modify PsychGetDoubleArg to accept all types and coerce them, then we could simplify by 
	calling that instead of doing all of this stuff...
		
*/
psych_bool PsychAllocInFlagArg(int position,  PsychArgRequirementType isRequired, psych_bool **argVal)
{
	const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool			acceptArg;

	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType)(PsychArgType_double|PsychArgType_char|PsychArgType_uint8|PsychArgType_boolean), 
									isRequired, 1,1,1,1,kPsychUnusedArrayDimension,kPsychUnusedArrayDimension);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		//unlike other PsychAllocIn* functions, here we allocate new memory and copy the input to it rather than simply returning a pointer into the received array.
		//That's because we want the booleans returned to the caller by PsychAllocInFlagArg() to alwyas be 8-bit booleans, yet we accept as flags either 64-bit double, char, 
		//or logical type.  Restricting to logical type would be a nuisance in the MATLAB environment and does not solve the problem because on some platforms MATLAB
		//uses for logicals 64-bit doubles and on others 8-bit booleans (check your MATLAB mex/mx header files).     
	    *argVal = (psych_bool *)mxMalloc(sizeof(psych_bool));
		mxPtr = PsychGetInArgMxPtr(position);
		if(mxIsLogical(mxPtr)){
			if(mxGetLogicals(mxPtr)[0])
				**argVal=(psych_bool)1;
			else
				**argVal=(psych_bool)0;
		}else{	
			if(mxGetScalar(mxPtr))
				**argVal=(psych_bool)1;
			else
				**argVal=(psych_bool)0;
		}
	}
	return(acceptArg);    //the argument was not present (and optional).	
}


psych_bool PsychAllocInFlagArgVector(int position,  PsychArgRequirementType isRequired, int *numElements, psych_bool **argVal)
{
	const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool			acceptArg;
	int				i;

	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	// MK: Disabled. Doesn't work without conversion of mxGetData into many subcases...
	// PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType)(PsychArgType_double | PsychArgType_char | PsychArgType_uint8 | PsychArgType_boolean), 
	//		       isRequired, 1, kPsychUnboundedArraySize, 1, kPsychUnboundedArraySize, kPsychUnusedArrayDimension, kPsychUnusedArrayDimension);

	// Ok. Let's see if anybody ever complains about this...
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType)(PsychArgType_double | PsychArgType_boolean), 
				       isRequired, 1, kPsychUnboundedArraySize, 1, kPsychUnboundedArraySize, kPsychUnusedArrayDimension, kPsychUnusedArrayDimension);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		if ((psych_uint64) mxGetM(mxPtr) * (psych_uint64) mxGetN(mxPtr) >= INT_MAX) {
			printf("PTB-ERROR: %i th input argument has more than 2^31 - 1 elements! This is not supported.\n", position);
			*numElements = 0;
			return(FALSE);
		}
		
		*numElements = (int) (mxGetM(mxPtr) * mxGetN(mxPtr));
		//unlike other PsychAllocIn* functions, here we allocate new memory and copy the input to it rather than simply returning a pointer into the received array.
		//That's because we want the booleans returned to the caller by PsychAllocInFlagArgVector() to alwyas be 8-bit booleans, yet we accept as flags either 64-bit double, char, 
		//or logical type.  Restricting to logical type would be a nuisance in the MATLAB environment and does not solve the problem because on some platforms MATLAB
		//uses for logicals 64-bit doubles and on others 8-bit booleans (check your MATLAB mex/mx header files).     		
	    *argVal = (psych_bool *)mxMalloc(sizeof(psych_bool) * ((size_t) *numElements));
		for(i=0; i< *numElements;i++){
			if(mxIsLogical(mxPtr)){
				if(mxGetLogicals(mxPtr)[i])
					(*argVal)[i]=(psych_bool)1;
				else
					(*argVal)[i]=(psych_bool)0;
			}else{
				if(mxGetPr(mxPtr)[i])
					(*argVal)[i]=(psych_bool)1;
				else
					(*argVal)[i]=(psych_bool)0;
			}
		}
	}
	return(acceptArg);    //the argument was not present (and optional).	
}


/*
	PsychCopyInFlagArg()
*/
psych_bool PsychCopyInFlagArg(int position, PsychArgRequirementType isRequired, psych_bool *argVal)
{
	const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool			acceptArg;
	
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType)(PsychArgType_double|PsychArgType_char|PsychArgType_uint8|PsychArgType_boolean), 
									isRequired, 1,1,1,1,kPsychUnusedArrayDimension,kPsychUnusedArrayDimension);
	matchError=PsychMatchDescriptors();
	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		if(mxIsLogical(mxPtr)){
			if(mxGetLogicals(mxPtr)[0])
				*argVal=(psych_bool)1;
			else
				*argVal=(psych_bool)0;
		}else{	
			if(mxGetScalar(mxPtr))
				*argVal=(psych_bool)1;
			else
				*argVal=(psych_bool)0;
		}	
	}
	return(acceptArg);    //the argument was not present (and optional).	
}



psych_bool PsychCopyOutFlagArg(int position, PsychArgRequirementType isRequired, psych_bool argVal)
{
	return(PsychCopyOutDoubleArg(position, isRequired, (double)argVal));
}


/*
    PsychAllocOutFlagListArg()
	
	This seems silly.  Find out where its used and consider using an array of booleans instead.  Probably the best thing
	is just to transparently map arrays of booleans to logical arrays MATLAB.  
    
    In Matlab our psych_bool flags are actually doubles.  This will not be so in all scripting languages.  We disguise the 
    implementation of psych_bool flags within the scripting envrironment by making the flag list opaque and
    providing accessor fucntions PsychLoadFlagListElement, PsychSetFlagListElement, and PsychClearFlagListElement.
    
    TO DO: maybe this should return a logical array instead of a bunch of doubles.  Itwould be better for modern versions
	of MATLAB which store doubles as bytes internally.  
	

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
	


// functions which allocate native types without assigning them to return arguments.
// this is useful for embedding native structures within each other. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* 
    PsychAllocateNativeDoubleMat()
    
    Create an opaque native matrix.   Return both 
        - Its handle,  which is specified when nesting the native matrix nesting withing other native types.
        - A handle to the C array of doubles enclosed by the native type.
        
    If (*cArray != NULL) we copy m*n*p elements from cArray into the native matrix, otherwise not.
	In any case, *cArray will point to the C array of doubles enclosed by the native type in the end.

*/
void 	PsychAllocateNativeDoubleMat(psych_int64 m, psych_int64 n, psych_int64 p, double **cArray, PsychGenericScriptType **nativeElement)
{
    double *cArrayTemp;

	PsychCheckmWSizeLimits(m, n, p);
    *nativeElement = mxCreateDoubleMatrix3D(m,n,p);
    cArrayTemp = mxGetPr(*nativeElement);
    if(*cArray != NULL) memcpy(cArrayTemp, *cArray, sizeof(double) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    *cArray=cArrayTemp; 
}


double PsychGetNanValue(void)
{
	return(mxGetNaN());
}

/* PsychAllocInCharFromNativeArg()
 *
 * Given a pointer to a native PsychGenericScriptType datatype which represents character
 * strings in the runtime's native encoding, try to extract a standard char-string of it and
 * return it in the referenced char *str. Return TRUE on success, FALSE on failure, e.g., because
 * the nativeCharElement didn't contain a parseable string.
 *
 */
psych_bool PsychAllocInCharFromNativeArg(PsychGenericScriptType *nativeCharElement, char **str)
{
	mxArray	*mxPtr;
	int		status;
	psych_uint64 strLen;	

	*str = NULL;
	mxPtr  = (mxArray*) nativeCharElement;
	strLen = ((psych_uint64) mxGetM(mxPtr) * (psych_uint64) mxGetNOnly(mxPtr) * (psych_uint64) sizeof(mxChar)) + 1;
	if (strLen >= INT_MAX) PsychErrorExitMsg(PsychError_user, "Tried to pass in a string with more than 2^31 - 1 characters. Unsupported!");
	
	*str   = (char *) PsychCallocTemp((size_t) strLen, sizeof(char));
	status = mxGetString(mxPtr, *str, (mwSize) strLen); 
	if(status!=0) return(FALSE);
	return(TRUE);
}

/* PsychRuntimeGetPsychtoolboxRoot()
 *
 * Try to retrieve filesystem path to Psychtoolbox root folder (the result from PsychtoolboxRoot() in Matlab/Octave)
 * or users configuration folder (the result from PsychtoolboxConfigDir() in Matlab/Octave) from runtime. The result,
 * if any, will be cached for later fast lookup.
 *
 * getConfigDir = TRUE => Return PsychtoolboxConfigDir().
 * getConfigDir = FALSE => Return PsychtoolboxRoot(). 
 *
 * This function may fail to retrieve the path, in which case it returns an empty null-terminated string, i.e., strlen() == 0.
 * On successfull recovery of the path, returns a const char* to a readonly string which encodes the path.
 *
 */
const char* PsychRuntimeGetPsychtoolboxRoot(psych_bool getConfigDir)
{
	static psych_bool firstTime = TRUE;
	static char	psychtoolboxRootPath[FILENAME_MAX+1];
	static char psychtoolboxConfigPath[FILENAME_MAX+1];
	PsychGenericScriptType* myPathvar = NULL;
	char*	myPathvarChar = NULL;

	if (firstTime) {
		// Reset firstTime flag:
		firstTime = FALSE;
		
		// Init to null-terminated empty strings, so it is well-defined in case of error:
		psychtoolboxRootPath[0] = 0;
		psychtoolboxConfigPath[0] = 0;

		// Call into runtime to persuade it to create a global ptb_RootPath variable with
		// the path to the root folder: This will return zero on success.
		if (0 == PsychRuntimeEvaluateString("global ptb_RootPath; ptb_RootPath = ''; if exist('PsychtoolboxRoot', 'file'), ptb_RootPath = PsychtoolboxRoot; end;")) {
			// Success. Try to retrieve global ptb_RootPath from runtime:
			if (PsychRuntimeGetVariablePtr("global", "ptb_RootPath", &myPathvar)) {
				// Success, got a read-only pointer to variable. Try to convert into char* string:
				if (PsychAllocInCharFromNativeArg(myPathvar, &myPathvarChar)) {
					// Success. myPathvarChar is a temporary char string, so copy it to our persistent memory:
					strncpy(psychtoolboxRootPath, myPathvarChar, FILENAME_MAX);
				}
			}
		}
		
		// At this point we did our best and psychtoolboxRootPath is valid: Either a path string,
		// or an empty string signalling failure to get the path.

		// Same game again for PsychtoolboxConfigDir:
		if (0 == PsychRuntimeEvaluateString("global ptb_ConfigPath; ptb_ConfigPath = ''; if exist('PsychtoolboxConfigDir', 'file'), ptb_ConfigPath = PsychtoolboxConfigDir; end;")) {
			// Success. Try to retrieve global ptb_RootPath from runtime:
			if (PsychRuntimeGetVariablePtr("global", "ptb_ConfigPath", &myPathvar)) {
				// Success, got a read-only pointer to variable. Try to convert into char* string:
				if (PsychAllocInCharFromNativeArg(myPathvar, &myPathvarChar)) {
					// Success. myPathvarChar is a temporary char string, so copy it to our persistent memory:
					strncpy(psychtoolboxConfigPath, myPathvarChar, FILENAME_MAX);
				}
			}
		}

	}

	// Return whatever we've got:
	return((getConfigDir) ? &psychtoolboxConfigPath[0] : &psychtoolboxRootPath[0]);
}

/* PsychCopyInPointerArg() - Copy in a void* memory pointer which is
 * encoded as a 32 bit or 64 bit unsigned integer, depending if this
 * is a 32 bit or 64 bit build of Psychtoolbox.
 */
psych_bool PsychCopyInPointerArg(int position, PsychArgRequirementType isRequired, void **ptr)
{
	const mxArray 	*mxPtr;
	PsychError		matchError;
	psych_bool		acceptArg;
	psych_bool		is64Bit;
	
	// 64 bit build?
	is64Bit = sizeof(size_t) > 4;
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
	PsychSetSpecifiedArgDescriptor(position, PsychArgIn, ((is64Bit) ? PsychArgType_uint64 : PsychArgType_uint32), isRequired, 1,1,1,1,1,1);
	matchError=PsychMatchDescriptors();

	acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
	if(acceptArg){
		mxPtr = PsychGetInArgMxPtr(position);
		
		if (is64Bit) {
			*ptr = (void*) (size_t) (((psych_uint64*) mxGetData(mxPtr))[0]);
		} else {
			*ptr = (void*) (size_t) (((psych_uint32*) mxGetData(mxPtr))[0]);
		}
	}

	return(acceptArg); 
}


/* PsychCopyOutPointerArg() - Copy out a void* memory pointer which gets
 * encoded as a 32 bit or 64 bit unsigned integer, depending if this
 * is a 32 bit or 64 bit build of Psychtoolbox.
 */
psych_bool PsychCopyOutPointerArg(int position, PsychArgRequirementType isRequired, void* ptr)
{
	mxArray **mxpp;
	PsychError matchError;
	psych_bool putOut;
	psych_bool is64Bit;
	mwSize dimArray[2];
	int numDims = 2;
	dimArray[0] = dimArray[1] = 1;

	// 64 bit build?
	is64Bit = sizeof(size_t) > 4;
	
	PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
	PsychSetSpecifiedArgDescriptor(position, PsychArgOut, ((is64Bit) ? PsychArgType_uint64 : PsychArgType_uint32), isRequired, 1,1,1,1,0,0);
	matchError=PsychMatchDescriptors();

	putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
	if(putOut){
		mxpp = PsychGetOutArgMxPtr(position);
		
		if (is64Bit) {
			*mxpp = mxCreateNumericArray(numDims, (mwSize*) dimArray, mxUINT64_CLASS, mxREAL);
			((psych_uint64*) mxGetData(*mxpp))[0] = (psych_uint64) ((size_t) ptr);
		} else {
			*mxpp = mxCreateNumericArray(numDims, (mwSize*) dimArray, mxUINT32_CLASS, mxREAL);
			((psych_uint32*) mxGetData(*mxpp))[0] = (psych_uint32) ((size_t) ptr);
		}
	}

	return(putOut);
}

/* PsychRuntimePutVariable()
 *
 * CAUTION: Not supported on Matlab R2006b or earlier on the MS-Windows platform! Function will fail on such
 *			runtimes. PTB functions that rely on this function will not work on pre R2007a on Windows!
 *
 * Copy a given native variable of type PsychGenericScriptType, e.g., as created by PsychAllocateNativeDoubleMat()
 * in case of a double matrix, as a new variable into a specified workspace.
 *
 * 'workspace'	Namestring of workspace: "base" copy to base workspace. "caller" copy into calling functions workspace,
 *				'global' create new global variable with given name.
 *
 *				CAUTION:	Matlab and Octave show different behaviour when using the "caller" workspace! It is strongly
 *							recommended to avoid the "caller" workspace to avoid ugly compatibility bugs!!
 *
 * 'variable'	Name of the new variable.
 *
 * 'pcontent'	The actual content that should be copied into the variable.
 *
 *
 * Example: You want to create a double matrix with (m,n,p) rows/cols/layers as a variable 'myvar' in the base
 *          workspace and initialize it with content from the double array mycontent:
 *
 *          PsychGenericScriptType* newvar = NULL;
 *			double* newvarcontent = mycontent; // mycontent is double* onto existing data.
 *          PsychAllocateNativeDoubleMat(m, n, p, &newvarcontent, &newvar);
 *			At this point, newvar contains the content of 'mycontent' and 'newvarcontent' points to
 *			the copy. You could alter mycontent now without affecting the content of newvarcontent or newvar.
 *
 *			Create the corresponding variable in the base workspace:
 *			PsychRuntimePutVariable("base", "myvar", newvar);
 *
 *          The calling M-File etc. can access the content newvarcontent under the variable name 'myvar'.
 *
 *			As usual, the double matrix newvarcontent will be auto-destroyed when returning to the runtime,
 *			but the variable 'myvar' will remain valid until it goes out of scope.
 *
 * Returns zero on success, non-zero on failure.
 */
int PsychRuntimePutVariable(const char* workspace, const char* variable, PsychGenericScriptType* pcontent)
{
	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		return(mexPutVariable(workspace, variable, pcontent));
	#else
		PsychErrorExitMsg(PsychError_unimplemented, "Function PsychRuntimePutVariable() not yet supported for this runtime system!");
	#endif
} 

/* PsychRuntimeGetVariable()
 *
 * CAUTION: Not supported on Matlab R2006b or earlier on the MS-Windows platform! Function will fail on such
 *			runtimes. PTB functions that rely on this function will not work on pre R2007a on Windows!
 *
 * Copy a given native variable of type PsychGenericScriptType, as a new variable from a specified workspace.
 *
 * 'workspace'	Namestring of workspace: "base" copy from base workspace. "caller" copy from calling functions workspace,
 *				'global' get global variable with given name.
 *
 *				CAUTION:	Matlab and Octave show different behaviour when using the "caller" workspace! It is strongly
 *							recommended to avoid the "caller" workspace to avoid ugly compatibility bugs!!
 *
 * 'variable'	Name of the variable to get a copy of.
 *
 * 'pcontent'	Pointer to a PsychGenericScriptType* where the location of the new variables content should be stored.
 *				The pointed to pointer will be set to NULL on failure.
 *
 * Returns TRUE on success, FALSE on failure.
 */
psych_bool PsychRuntimeGetVariable(const char* workspace, const char* variable, PsychGenericScriptType** pcontent)
{
	// Init to empty default:
	*pcontent = NULL;

	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		*pcontent = mexGetVariable(workspace, variable);

		// Return true on success, false on failure:
		return((*pcontent) ? TRUE : FALSE);
	#else
		PsychErrorExitMsg(PsychError_unimplemented, "Function PsychRuntimeGetVariable() not yet supported for this runtime system!");
	#endif
} 

/* PsychRuntimeGetVariablePtr()
 *
 * CAUTION: Not supported on Matlab R2006b or earlier on the MS-Windows platform! Function will fail on such
 *			runtimes. PTB functions that rely on this function will not work on pre R2007a on Windows!
 *
 * Retrieve a *read-only* pointer to a given native variable of type PsychGenericScriptType in the specified workspace.
 * The variable is not copied, just referenced, so you *must not modify/write to the location* only perform read access!
 *
 * 'workspace'	Namestring of workspace: "base" get from base workspace. "caller" get from calling functions workspace,
 *				'global' get global variable with given name.
 *
 *				CAUTION:	Matlab and Octave show different behaviour when using the "caller" workspace! It is strongly
 *							recommended to avoid the "caller" workspace to avoid ugly compatibility bugs!!
 *
 * 'variable'	Name of the variable to get a reference.
 *
 * 'pcontent'	Pointer to a PsychGenericScriptType* where the location of the variables content should be stored.
 *				The pointed to pointer will be set to NULL on failure.
 *
 * Returns TRUE on success, FALSE on failure.
 */
psych_bool PsychRuntimeGetVariablePtr(const char* workspace, const char* variable, PsychGenericScriptType** pcontent)
{
	// Init to empty default:
	*pcontent = NULL;

	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		*pcontent = (PsychGenericScriptType*) mexGetVariablePtr(workspace, variable);

		// Return true on success, false on failure:
		return((*pcontent) ? TRUE : FALSE);
	#else
		PsychErrorExitMsg(PsychError_unimplemented, "Function PsychRuntimeGetVariablePtr() not yet supported for this runtime system!");
	#endif
} 

/* PsychRuntimeEvaluateString()
 *
 * Simple function evaluation by scripting environment via feval() style functions.
 * This asks the runtime environment to execute/evaluate the given string 'cmdstring',
 * passing no return arguments back, except an error code.
 * 
 */
int	PsychRuntimeEvaluateString(const char* cmdstring)
{
	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		return(mexEvalString(cmdstring));
	#else
		int parse_status;
		// Call Octaves evalstring-parser:
		eval_string(std::string(cmdstring), (bool) TRUE, parse_status);
		return(parse_status);
	#endif
}

//end of Matlab & Octave only stuff.
#endif
