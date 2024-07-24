/*
 * PsychSourceGL/Source/Common/Base/PythonGlue/PsychScriptingGluePython.c
 *
 * AUTHORS:
 *
 * mario.kleiner.de@gmail.com   mk
 *
 * PLATFORMS: All -- Glue layer for CPython runtime environment.
 *
 * HISTORY:
 *
 * 19-June-2018     mk  Derived from PsychScriptingGlueMatlab.c
 *
 * DESCRIPTION:
 *
 * PsychScriptingGluePython defines abstracted functions to pass values
 * between the calling Python environment and the PsychToolbox CPython modules.
 *
 * Copyright (c) 2018 Mario Kleiner.
 *
 * MIT license:
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

// During inclusion of Psych.h, we define the special flag PTBINSCRIPTINGGLUE. This
// will cause some of the system headers in Psych.h not to be included during build
// of PsychScriptingGluePython.c:
#define PTBINSCRIPTINGGLUE 1
#include "Psych.h"
#undef PTBINSCRIPTINGGLUE

#if PSYCH_LANGUAGE == PSYCH_PYTHON

// Import NumPy array handling functions: Require at least NumPy v 1.13.0, released
// in June 2017:
#define NPY_NO_DEPRECATED_API NPY_1_13_API_VERSION
#include <numpy/arrayobject.h>

#if defined(Py_LIMITED_API) && defined(__GNUC__) && (__GNUC__ < 10) && !defined(__clang__)
#error \
"This version of gcc has a bug and cannot compile this code \
with Py_LIMITED_API enabled. Either build without Py_LIMITED_API \
or upgrade gcc to a version >= 10.1 and try again. Ubuntu 20.04, \
e.g., ships a suitable compiler package named gcc-10. \
See https://github.com/numpy/numpy/issues/16970 for details."
#endif

// Define this to 1 if you want lots of debug-output for the Python-Scripting glue.
#define DEBUG_PTBPYTHONGLUE 0

// Special hacks to allow Psychtoolbox to build on Python, stubbing out Mex Api replacements.

#include <string.h>
#include <setjmp.h>

////Static functions local to ScriptingGluePython.c.
// _____________________________________________________________________________________

#define MAX_CMD_NAME_LENGTH 100

//Static variables local to ScriptingGluePython.c.  The convention is to append a abbreviation in all
//caps of the C file name to the variable name.

// nameFirstGLUE, baseFunctionInvoked, nlhsGLUE, nrhsGLUE, plhsGLUE, prhsGLUE
// are state which has to be maintained for each invocation of a mex module.
// If a modules calls itself recursively, this state has to be maintained for
// each recursive call level. We implement a little homemade stack for those
// variables. Maximum stack depth and therefore maximum recursion level for
// reentrant recursive calls is MAX_RECURSIONLEVEL. The variable recLevel
// keeps track of the current call recursion level and acts as a "stack pointer".
// It gets incremented by each entry to mexFunction() and decremented by each
// regular exit from mexFunction(). On error abort or modules reload it needs
// to get reset to initial -1 state:
#define MAX_RECURSIONLEVEL 5
#define MAX_INPUT_ARGS 100
#define MAX_OUTPUT_ARGS 100

static psych_bool nameFirstGLUE[MAX_RECURSIONLEVEL];
static psych_bool baseFunctionInvoked[MAX_RECURSIONLEVEL];
static psych_bool use_C_memory_layout[MAX_RECURSIONLEVEL];

static int nlhsGLUE[MAX_RECURSIONLEVEL];  // Number of requested return arguments.
static int nrhsGLUE[MAX_RECURSIONLEVEL];  // Number of provided call arguments.

static PyObject* plhsGLUE[MAX_RECURSIONLEVEL][MAX_OUTPUT_ARGS];             // An array of pointers to the Python return arguments.
static PyObject* prhsGLUE[MAX_RECURSIONLEVEL][MAX_INPUT_ARGS];              // An array of pointers to the Python call arguments.
static psych_bool prhsNeedsConversion[MAX_RECURSIONLEVEL][MAX_INPUT_ARGS];  // prhsGLUE needs one-time conversion to NumPy array?

static int recLevel = -1;
static psych_bool psych_recursion_debug = FALSE;
static int psych_refcount_debug = 0;

// Our own module object:
static PyObject *module = NULL;

// Full filesystem path/name to the library (DLL/dylib/libso) that defines this module:
static char modulefilename[FILENAME_MAX];

// MODULE INITIALIZATION FOR PYTHON:
// =================================

#define PPYINIT(...) _PPYINIT(__VA_ARGS__)
#define _PPYNAME(n) #n
#define PPYNAME(...) _PPYNAME(__VA_ARGS__)

static PyMethodDef GlobalPythonMethodsTable[] = {
    {PPYNAME(PTBMODULENAME), PsychScriptingGluePythonDispatch, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

// Python 3 init code:
#define _PPYINIT(n) PyMODINIT_FUNC PyInit_ ## n(void)

// Defined in PsychScriptingGluePython.c
PsychError PsychExitPythonGlue(void);

/* PythonModuleCleanup() - Call Python specific cleanup function.
 *
 * This cleanup function is only called on Python 3, and as far as i
 * understand only at interpreter shutdown time, ie. when calling quit()
 * or pressing CTRL+D. reload()ing of extension modules seems to be not
 * possible in Python 2 and 3, as described in PEP 0498:
 * https://www.python.org/dev/peps/pep-0489/#id29
 *
 * For this reason, this function is of limited value, but implemented
 * anyway for completeness and future reference.
 *
 */
void PythonModuleCleanup(void* userptr)
{
    (void) userptr;
    (void) PsychExitPythonGlue();
}

static struct PyModuleDef module_definition = {
    PyModuleDef_HEAD_INIT,                                                      // Base instance.
    PPYNAME(PTBMODULENAME),                                                     // Module name.
    "The " PPYNAME(PTBMODULENAME) " Psychtoolbox module for Python 3.\n"        // Help text.
    "Copyright (c) 2018-2023 Mario Kleiner.\n"
    "Copyright (c) 2022-2023 Alex Forrence.\n"
    "Licensed under the MIT license.",
    -1,                                                                         // -1 = No sub-interpreter support: https://docs.python.org/3/c-api/module.html#c.PyModuleDef
    GlobalPythonMethodsTable,                                                   // Function dispatch table, shared with Python 2.
    NULL,                                                                       // m_slots
    NULL,                                                                       // m_traverse
    NULL,                                                                       // m_clear
    PythonModuleCleanup                                                         // m_free = PythonModuleCleanup, cleanup at module destruction.
};

// This is the entry point - module init function, called at module import:
// PTBMODULENAME is -DPTBMODULENAME myname defined by the build script to the
// name of the module, e.g., GetSecs.
PPYINIT(PTBMODULENAME)
{
    modulefilename[0] = 0;

    // Add a help string with module synopsis to 1st function - our main dispatch function:
    GlobalPythonMethodsTable[0].ml_doc = PsychBuildSynopsisString(PPYNAME(PTBMODULENAME));

    // Initialize module:
    module = PyModule_Create(&module_definition);
    return(module);
}

// END OF MODULE INITIALIZATION FOR PYTHON:

// Return filename of the module definition file - the shared library:
const char* PsychGetPyModuleFilename(void)
{
    // Get full filesystem path/name of the module definition file, ie. the library:
    if (module && !modulefilename[0]) {
        PyObject *fname = PyModule_GetFilenameObject(module);

        if (fname)
            mxGetString(fname, modulefilename, sizeof(modulefilename) - 1);
        else
            sprintf(modulefilename, "%s", PyModule_GetFilename(module));
        Py_XDECREF(fname);
    }

    return(&modulefilename[0]);
}

// This jump-buffer stores CPU- and stackstate at the position
// where our octFunction() dispatcher actually starts executing
// the functions or subfunctions code. That is the point where
// a PsychErrorExit() or PsychErrorExitMsg() will return control...
jmp_buf jmpbuffer[MAX_RECURSIONLEVEL];

// Error exit handler:
// Prints the error-string with CPythons printing facilities, and then longjmp's
// to the cleanup routine at the end of our PsychScriptingGluePythonDispatch()
// dispatcher.
void mexErrMsgTxt(const char* s) {
    if (s && strlen(s) > 0)
        printf("%s:%s: %s\n", PsychGetModuleName(), PsychGetFunctionName(), s);
    else
        printf("%s:%s\n", PsychGetModuleName(), PsychGetFunctionName());

    // Use the jump-buffer to unwind the stack...
    longjmp(jmpbuffer[recLevel], 1);
}

// Interface to printf... TODO Used anywhere?
void mexPrintf(const char* fmt, ...)
{
    va_list args;
    va_start (args, fmt);
    vfprintf(stderr, fmt, args);
    vfprintf(stdout, fmt, args);
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

double mxGetNaN(void)
{
    /* TODO FIXME Stop-gap? Better solution in Python? */
    return(nan(""));
}

int mxIsLogical(const PyObject* a)
{
    return(PyBool_Check(a) || PyArray_ISBOOL((const PyArrayObject*) a));
}

int mxIsCell(const PyObject* a)
{
    return(PyTuple_Check(a));
}

int mxIsStruct(const PyObject* a)
{
    return((PyList_Check(a) && (PyList_Size((PyObject*) a) > 0) && PyDict_Check(PyList_GetItem((PyObject*) a, 0))) ||
            PyDict_Check(a));
}

int mxIsNumeric(const PyObject* a)
{
    return(PyFloat_Check(a) || PyLong_Check(a));
}

int mxIsChar(const PyObject* a)
{
    return(PyUnicode_Check(a) || PyBytes_Check(a));
}

int mxIsSingle(const PyObject* a)
{
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_FLOAT);
}

int mxIsDouble(const PyObject* a)
{
    //return(PyArray_ISFLOAT((const PyArrayObject*) a));
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_DOUBLE);
}

int mxIsUint8(const PyObject* a)
{
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_UINT8);
}

int mxIsUint16(const PyObject* a)
{
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_UINT16);
}

int mxIsUint32(const PyObject* a)
{
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_UINT32);
}

int mxIsUint64(const PyObject* a)
{
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_UINT64);
}

int mxIsInt8(const PyObject* a)
{
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_INT8);
}

int mxIsInt16(const PyObject* a)
{
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_INT16);
}

int mxIsInt32(const PyObject* a)
{
//    return(PyArray_ISINTEGER((PyArrayObject*) a));
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_INT32);
}

int mxIsInt64(const PyObject* a)
{
    return(PyArray_TYPE((PyArrayObject*) a) == NPY_INT64);
}

int PsychGetNumTypeFromArgType(PsychArgFormatType type)
{
    switch(type) {
        case PsychArgType_uint8:
            return(NPY_UINT8);

        case PsychArgType_uint16:
            return(NPY_UINT16);

        case PsychArgType_uint32:
            return(NPY_UINT32);

        case PsychArgType_uint64:
            return(NPY_UINT64);

        case PsychArgType_int8:
            return(NPY_INT8);

        case PsychArgType_int16:
            return(NPY_INT16);

        case PsychArgType_int32:
            return(NPY_INT32);

        case PsychArgType_int64:
            return(NPY_INT64);

        case PsychArgType_single:
            return(NPY_FLOAT);

        case PsychArgType_double:
            return(NPY_DOUBLE);

        case PsychArgType_boolean:
            return(NPY_BOOL); // 1-Byte

        case PsychArgType_char:
            return(NPY_STRING);

        case PsychArgType_cellArray:
            return(NPY_OBJECT);

        case PsychArgType_structArray:
            return(NPY_OBJECT);

        default:
            printf("PTB-CRITICAL: Can not map PsychArgFormatType %i to NumPy type!\n", type);
            PsychErrorExitMsg(PsychError_invalidArg_type, "Unknown PsychArgFormatType encountered. Don't know how to map it to NumPy.");
            return(NPY_NOTYPE);
    }
}

PyObject* mxCreateNumericArray(int numDims, ptbSize dimArray[], PsychArgFormatType arraytype)
{
    int typenum = PsychGetNumTypeFromArgType(arraytype);

    // Create empty/uninitialized array in (0) C contiguous style if use_C_memory_layout,
    // else create in (1) Fortran contiguous style:
    return(PyArray_EMPTY(numDims, (npy_intp*) dimArray, typenum, (use_C_memory_layout[recLevel]) ? 0 : 1));
}

PyObject* mxCreateDoubleMatrix(ptbSize rows, ptbSize cols)
{
    ptbSize dims[2];
    dims[0] = rows;
    dims[1] = cols;

    return(mxCreateNumericArray(2, dims, PsychArgType_double));
}

PyObject* mxCreateLogicalMatrix(ptbSize rows, ptbSize cols)
{
    ptbSize dims[2];
    dims[0] = rows;
    dims[1] = cols;
    return(mxCreateNumericArray(2, dims, PsychArgType_boolean));
}

PyObject* mxCreateString(const char* instring)
{
    PyObject* ret;

    if (!instring)
        return(PyUnicode_FromString("NULL"));

    // Try decoding from UTF-8:
    ret = PyUnicode_FromString(instring);
    PyErr_Clear();

    // On Windows, some low-level ANSI api's, e.g., input device enumeration in
    // PsychHID, return strings encoded in Windows codepage CP_ACP, the system
    // active code page. See the following Microsoft documentation for details:
    //
    // https://docs.microsoft.com/en-us/windows/desktop/Intl/conventions-for-function-prototypes
    //
    // PyUnicode_DecodeMBCS() decodes such multibyte character strings encoded in
    // the CP_ACP current system code page into Unicode. Therefore this is our
    // 1st fallback if UTF-8 decoding fails on Windows:
    // Thanks to Hiroyuki Sogo for finding this solution!
    #if PSYCH_SYSTEM == PSYCH_WINDOWS
    if (!ret) {
        ret = PyUnicode_DecodeMBCS(instring, strlen(instring), NULL);
        PyErr_Clear();
    }
    #endif

    if (!ret) {
        // Try decoding assuming current system locale setting:
        ret = PyUnicode_DecodeLocale(instring, "surrogateescape");
        PyErr_Clear();

        // If Py_LIMITED_API is enabled, it means we're compiling on
        // at least Python 3.7. If that's the case, we can skip the
        // next two checks, which deal with incompatibilties in earlier
        // Python versions
        #ifndef Py_LIMITED_API
        // Retry with strict error handler, because of backwards incompatible
        // change in Python 3.6 -> 3.7 (sigh):
        if (!ret) {
            ret = PyUnicode_DecodeLocale(instring, "strict");
            PyErr_Clear();
        }

        // Retry with Python startup locale. There were backwards incompatible
        // changes in Python 3.0 -> 3.2 -> 3.6, often only on Windows (sigh):
        if (!ret) {
            ret = PyUnicode_DecodeFSDefault(instring);
            PyErr_Clear();
        }
        #endif
    }

    if (!ret) {
        // Final attempt: Decode as Latin-1, corresponding to the first 256
        // Unicode codepoints. According to CPython source, this can essentially
        // only fail on system out-of-memory, and either returns a Latin-1
        // string if instring was Latin-1, or it will return whatever it is
        // crammed into Latin-1, ie. Byte values in range 0-255. Good enough for
        // debugging by calling script:
        ret = PyUnicode_DecodeLatin1(instring, strlen(instring), NULL);
        PyErr_Clear();
    }

    // Better than crashing:
    if (!ret)
        ret = PyUnicode_FromString("WARNING: INVALID UNDECODABLE STRING!");

    // Clear potential Unicode exceptions set by above:
    PyErr_Clear();

    return(ret);
}

void* mxGetData(const PyObject* arrayPtr)
{
    return(PyArray_DATA((PyArrayObject*)  arrayPtr));
}

double* mxGetPr(const PyObject* arrayPtr)
{
    return(PyArray_DATA((PyArrayObject*)  arrayPtr));
}

double mxGetScalar(const PyObject* arrayPtr)
{
    return(PyFloat_AsDouble((PyObject*) arrayPtr));
}

PsychNativeBooleanType* mxGetLogicals(const PyObject* arrayPtr)
{
    // PsychNativeBooleanType == psych_bool == unsigned char == 1 Byte on
    // Python with NumPy:
    return((PsychNativeBooleanType*) mxGetData(arrayPtr));
}

ptbSize mxGetNumberOfDimensions(const PyObject* arrayPtr)
{
    if (!PyArray_Check(arrayPtr))
        return(0);

    return((ptbSize) PyArray_NDIM((const PyArrayObject*) arrayPtr));
}

ptbSize mxGetM(const PyObject* arrayPtr)
{
    if (mxGetNumberOfDimensions(arrayPtr) < 1)
        return(1);

    return((ptbSize) PyArray_DIM((const PyArrayObject*)  arrayPtr, 0));
}

ptbSize mxGetN(const PyObject* arrayPtr)
{
    if (mxGetNumberOfDimensions(arrayPtr) < 2)
        return(1);

    return((ptbSize) PyArray_DIM((const PyArrayObject*)  arrayPtr, 1));
}

/*
 *    Get the 2nd array dimension.
 *
 *    The Mex API's mxGetN is sometimes undersirable because it returns the product of all dimensions above 1.  Our mxGetNOnly only returns N, for when you need that.
 *
 *    The abstracted Psychtoolbox API supports matrices with up to 3 dimensions.
 */
static ptbSize mxGetNOnly(const PyObject *arrayPtr)
{
    return(mxGetN(arrayPtr));
}

/*
 *    Get the third array dimension which we call "P".  mxGetP should act just like mxGetM and mxGetN.
 *
 *    The abstracted Psychtoolbox API supports matrices with up to 3 dimensions.
 */
static ptbSize mxGetP(const PyObject *arrayPtr)
{
    if (mxGetNumberOfDimensions(arrayPtr) < 3) {
        //printf("P %i\n", 1);
        return(1);
    }

    return((ptbSize) PyArray_DIM((const PyArrayObject*)  arrayPtr, 2));
}

int mxGetString(PyObject* arrayPtr, char* outstring, int outstringsize)
{
    if (!mxIsChar(arrayPtr))
        PsychErrorExitMsg(PsychError_internal, "FATAL Error: Tried to convert a non-string into a string!");

    // Python 3: No PyObject_Unicode(), distinguish unicode input vs. bytes 8-bit legacy string input:
    if (PyUnicode_Check(arrayPtr))
        // Provide it as Latin1 8-bit "bytes" string from unicode, giving a new reference:
        arrayPtr = PyUnicode_AsLatin1String(arrayPtr);
    else
        // Is already a 8-bit "bytes" string. Increment refcount, to counteract decref below:
        Py_INCREF(arrayPtr);

    // Got a 8-bit "bytes" string?
    if (arrayPtr) {
        // Extract as const char* C-style string - copy to return char array:
        int rc = ((snprintf(outstring, outstringsize, "%s", PyBytes_AsString(arrayPtr))) >= 0) ? 0 : 1;

        Py_DECREF(arrayPtr);
        return(rc);
    }
    else
        return(1);
}

void mxDestroyArray(PyObject *arrayPtr)
{
    // Destroy a PyObject:
    if (arrayPtr == NULL) return;

    fprintf(stderr, "WARN WARN UNIMPLEMENTED: mxDestroyArray()\n");

	return;
}

PyObject* mxCreateStructArray(int numDims, ptbSize* ArrayDims, int numFields, const char** fieldNames)
{
    int i, j, n;
    PyObject* retval = NULL;

    if (numDims != 1)
        PsychErrorExitMsg(PsychError_unimplemented, "Error: mxCreateStructArray: Anything else than 1D Struct-Array is not supported!");

    if (numFields < 1)
        PsychErrorExitMsg(PsychError_internal, "Error: mxCreateStructArray: numFields < 1 ?!?");

    n = (int) ArrayDims[0];

    if (n < -1)
        PsychErrorExitMsg(PsychError_internal, "Error: mxCreateStructArray: Negative number of array elements requested?!?");

    // Create a list of objects - the struct array - with each object being a
    // dictionary that contains fieldNames as keys, and PyObjects as the actual
    // values -- iow. a single "array of structs":

    // Create to-be-returned list that makes up the struct array, except if n == -1,
    // in which case we don't return a list (~ array), but just the single dict (~ struct):
    if (n != -1)
        retval = PyList_New((Py_ssize_t) n);

    // Create one dictionary for each slot:
    for (i = 0; i < abs(n); i++) {
        PyObject* slotdict = PyDict_New();

        // Create all fields for all fieldNames for this slots dictionary:
        for (j = 0; j < numFields; j++) {
            // Init value with Py_None:
            Py_INCREF(Py_None);
            if (PyDict_SetItemString(slotdict, fieldNames[j], Py_None))
                PsychErrorExitMsg(PsychError_internal, "Error: mxCreateStructArray: Failed to init struct-Array slot with item!");
        }

        // For n >=  0, assign to i'th slot of returned list retval.
        // For n == -1, directly return our one and only slotdict as retval:
        if (n > -1)
            PyList_SetItem(retval, i, slotdict);
        else
            retval = slotdict;
    }

    return(retval);
}

PyObject* mxGetField(const PyObject* structArray, int index, const char* fieldName)
{
    if (!mxIsStruct(structArray))
        PsychErrorExitMsg(PsychError_internal, "Error: mxGetField: Tried to manipulate something other than a struct-Array!");

    // Different code-path for single element structArray aka a dict, vs. multi-element
    // structArray aka list of dicts:
    if (!PyDict_Check(structArray)) {
        if (index >= PyList_Size((PyObject*) structArray))
            PsychErrorExitMsg(PsychError_internal, "Error: mxGetField: Index exceeds size of struct-Array!");

        return(PyDict_GetItemString(PyList_GetItem((PyObject*) structArray, index), fieldName));
    }
    else {
        if (index != 0)
            PsychErrorExitMsg(PsychError_internal, "Error: mxGetField: Index exceeds size of struct-Array!");

        return(PyDict_GetItemString((PyObject*) structArray, fieldName));
    }
}

int mxIsField(PyObject* structArray, const char* fieldName)
{
    if (!mxIsStruct(structArray))
        PsychErrorExitMsg(PsychError_internal, "Error: mxIsField: Tried to manipulate something other than a struct-Array!");

    // Ok, cheating. If the fieldName exists, we always returns a field number of 1,
    // otherwise we return -1, so this function can only check if a fieldName is valid.
    // But then, that's all that this function is used for inside our implementation,
    // so we should be fine.
    if (mxGetField(structArray, 0, fieldName))
        return(1);

    // No such field :(
    return(-1);
}

// The mxSetField() function unconditionally steals the reference to
// pStructInner from our caller. In general mxSetField() is only called by
// the PsychSetStructArrayXXXXXElement() functions, and these generate a new
// reference to whatever they assign, so now the reference in this slot+field is
// the only existing one, and pStructInner will get deleted, once pStructOuter
// has fulfilled its purpose inside the Python script and goes out of scope/gets
// deleted there:
void mxSetField(PyObject* pStructOuter, int index, const char* fieldName, PyObject* pStructInner)
{
    if (psych_refcount_debug && pStructInner)
        printf("PTB-DEBUG: In mxSetField: refcount of external object %p at enter is %li. %s\n",
               pStructInner, Py_REFCNT(pStructInner),
               (Py_REFCNT(pStructInner) > 1) ? "MIGHT leak if caller does not take care." : "");

    if (!mxIsStruct(pStructOuter)) {
        Py_XDECREF(pStructInner);
        PsychErrorExitMsg(PsychError_internal, "Error: mxSetField: Tried to manipulate something other than a struct-Array!");
    }

    PyObject *arraySlot;

    // Different code-path for single element structArray aka a dict, vs. multi-element
    // structArray aka list of dicts:
    if (PyList_Check(pStructOuter)) {
        if (index >= PyList_Size(pStructOuter)) {
            Py_XDECREF(pStructInner);
            PsychErrorExitMsg(PsychError_internal, "Error: mxSetField: Index exceeds size of struct-Array!");
        }

        // Get dictionary for slot 'index':
        arraySlot = PyList_GetItem(pStructOuter, index);
    }
    else {
        // Single-slot array: pStructOuter is already the single struct aka dict:
        arraySlot = pStructOuter;
    }

    // Assign pStructInner as new value of field fieldName in slot index:
    // This will drop the refcount of the previous value in that slot+field,
    // and increase the refcount of pStructInner by one, iow. it doesn't steal
    // a reference, but get our own one:
    if (PyDict_SetItemString(arraySlot, fieldName, pStructInner)) {
        Py_XDECREF(pStructInner);
        PsychErrorExitMsg(PsychError_internal, "Error: mxSetField: PyDict_SetItemString() failed!");
    }

    // The mxSetField() function unconditionally steals the reference to
    // pStructInner from our caller. In general mxSetField() is only called by
    // the PsychSetStructArrayXXXXXElement() functions, and these generate a new
    // reference to whatever they assign, so now the reference in this slot+field is
    // the only existing one, and pStructInner will get deleted, once pStructOuter
    // has fulfilled its purpose inside the Python script and goes out of scope/gets
    // deleted there:
    Py_XDECREF(pStructInner);
}

/* UNUSED, but here for reference */
PyObject* mxCreateCellArray(int numDims, ptbSize* ArrayDims)
{
    PyObject* retval;

    if (numDims > 1)
        PsychErrorExitMsg(PsychError_unimplemented, "Error: mxCreateCellArray: 2D Cell Arrays are not supported on Python build!");

    // Allocate our PyObject-Struct:
    retval = PyTuple_New((Py_ssize_t) ArrayDims[0]);

    // Done.
    return(retval);
}

/* UNUSED, but here for reference */
void mxSetCell(PsychGenericScriptType *cellVector, ptbIndex index, PyObject* mxFieldValue)
{
    if (!mxIsCell(cellVector))
        PsychErrorExitMsg(PsychError_internal, "Error: mxSetCell: Tried to manipulate something other than a cell-vector!");

    if (index >= (ptbIndex) PyTuple_Size(cellVector))
        PsychErrorExitMsg(PsychError_internal, "Error: mxSetCell: index tried to index beyond lenght of cell-vector!");

    PyTuple_SetItem(cellVector, index, mxFieldValue);

    return;
}


void mxSetLogical(PyObject* dummy)
{
    // This is a no-op, because it is not needed anywhere...
    return;
}

//local function declarations
static psych_bool PsychIsDefaultMat(const PyObject *mat);
static ptbSize mxGetP(const PyObject *array_ptr);
static ptbSize mxGetNOnly(const PyObject *arrayPtr);
static PyObject *mxCreateDoubleMatrix3D(psych_int64 m, psych_int64 n, psych_int64 p);

// firstTime: This flag defines if this is the first invocation of the module
// since it was (re-)loaded:
static psych_bool firstTime = TRUE;

PsychError PsychExitPythonGlue(void);
void ScreenCloseAllWindows(void);

void* init_numpy(void)
{
    import_array();
    return(NULL);
}

void PsychExitRecursion(void)
{
    if (recLevel < 0) {
        printf("PTB-CRITICAL: Recursion stack underflow in module %s! Brace for impact!\n", PsychGetModuleName());
        return;
    }

    if (psych_recursion_debug) printf("PTB-DEBUG: Module %s leaving recursive call level %i.\n", PsychGetModuleName(), recLevel);

    // Done with this call recursion level:
    recLevel--;
}


/*     PsychUseCMemoryLayoutIfOptimal() - Opt in to data exchange memory layout optimizations.
 *
 *     Tell scripting glue to use/assume a C programming language memory layout for exchanging
 *     multi-dimensional (== 2D, 3D, n-D) matrices with the scripting environment if that layout
 *     promises higher efficiency and performance in data exchange. This is an opt-in, requesting
 *     C-layout if 'tryEnableCMemoryLayout' = TRUE, otherwise standard Fortran layout is assumed.
 *     The default is Fortran layout if this function does not get called, and it resets to Fortran
 *     layout at each return of control to the calling scripting environment. Iow. it is a per-
 *     module subfunction-call opt-in.
 *     The function returns TRUE if C memory layout is engaged, otherwise FALSE is returned.
 *     The caller may have to adjust its own data processing according to the returned value,
 *     unless the function is called with tryEnableCMemoryLayout = FALSE or not called at all, in
 *     which case Fortran layout is the thing.
 *
 *     tryEnableCMemoryLayout = FALSE (default) Fortran classic style, TRUE = C-style.
 *
 *     Returns: TRUE if C-style is to be used, FALSE (default) if Fortran classic is to be used.
 *
 */
psych_bool PsychUseCMemoryLayoutIfOptimal(psych_bool tryEnableCMemoryLayout)
{
    // Python - more specifically NumPy - uses C programming language memory layout for
    // its n-D matrices and arrays. Using Fortran classic style requires memory layout
    // conversion when exchanging >= 2D data with Python/NumPy, which can increase
    // required memory bandwith, memory consumption, lead to cache trashing and therefore
    // reduced performance in input/output argument passing. Skipping that is beneficial
    // at least for large n-D arrays, so we always choose C-memory layout if caller does
    // opt-in to allow us to do it:
    use_C_memory_layout[recLevel] = tryEnableCMemoryLayout;

    if (DEBUG_PTBPYTHONGLUE)
        printf("PTB-DEBUG:%s:%s: %s C-Memory layout for NumPy based data exchange.\n",
               PsychGetModuleName(), PsychGetFunctionName(),
               (tryEnableCMemoryLayout) ? "Enabling" : "Disabling");

    return(tryEnableCMemoryLayout);
}


/*
 *
 *    Main entry point for Python runtime. Serves as a dispatch and handles
 *    first time initialization.
 *
 *        The subfunction dispatcher can operate in either of two modes depending
 *        on whether the module has registed subfunctions, or only a single "base"
 *        function.
 *
 *        subfunction mode:
 *        The examines the  first and second
 *        arguments for a string naming a module subfunction.  If it finds in either of those
 *        two arguments a string naming a module subfunctoin, then it looks up the approproate
 *        function pointer and invokes that function.  Before invoking the function the dispatcher
 *        removes the function name argument form the list of argumnets which was passed to the
 *        module.
 *
 *        base mode:  The dispatcher always invokes the same one subfunction and without
 *        alterinng the list of arguments.
 *
 *        Modules should now register in subfunction mode to support the build-in 'version' command.
 *
 */
PyObject* PsychScriptingGluePythonDispatch(PyObject* self, PyObject* args)
{
    psych_bool          isArgThere[2], isArgEmptyMat[2], isArgText[2], isArgFunction[2];
    PsychFunctionPtr    fArg[2], baseFunction;
    char                argString[2][MAX_CMD_NAME_LENGTH];
    PyObject*           tmparg = NULL;
    PyObject*           plhs = NULL;
    int                 i;
    int                 nrhs = (int) PyTuple_Size(args);

    if (!PyTuple_Check(args)) {
        printf("FAIL FAIL FAIL!\n");
        return(NULL);
    }

    // Initialization
    if (firstTime) {
        // Reset call recursion level to startup default:
        recLevel = -1;
        psych_recursion_debug = FALSE;
        psych_refcount_debug = 0;

        if (getenv("PSYCH_RECURSION_DEBUG")) psych_recursion_debug = TRUE;
        if (getenv("PSYCH_REFCOUNT_DEBUG")) psych_refcount_debug = atoi(getenv("PSYCH_REFCOUNT_DEBUG"));

        // Initialize NumPy array extension for use in *this compilation unit* only:
        (void) init_numpy();

        // Call the Psychtoolbox init function, which inits the Psychtoolbox and calls the project init.
        PsychInit();

        // Hard to believe, but apparently true: Python does not allow unloading extension modules!
        // As a workaround, we register a subfunction "Shutdown" that allows to trigger a manual module
        // shutdown, although not a true unload. At least this might allow module state resets if done
        // carefully:
        PsychRegister("Shutdown",  &PsychExitPythonGlue);

        // Register hidden helper function: This one dumps all registered subfunctions of
        // a module into a struct array of text strings. Needed by our automatic documentation
        // generator script to find out about subfunctions of a module:
        PsychRegister((char*) "DescribeModuleFunctionsHelper",  &PsychDescribeModuleFunctions);

        firstTime = FALSE;
    }

    // Increment call recursion level for this invocation of the module:
    recLevel++;
    if (recLevel >= MAX_RECURSIONLEVEL) {
        // Maximum level exceeded!
        printf("PTB-CRITICAL: Maximum recursion level %i for recursive calls into module '%s' exceeded!\n", recLevel, PsychGetModuleName());
        printf("PTB-CRITICAL: Aborting call sequence. Check code for recursion bugs!\n");
        recLevel--;
        return(NULL);
    }

    if (psych_recursion_debug) printf("PTB-DEBUG: Module %s entering recursive call level %i.\n", PsychGetModuleName(), recLevel);

    // Default to not using C memory layout, but classic (backwards compatible) Fortran layout:
    use_C_memory_layout[recLevel] = FALSE;

    // Save CPU-state and stack at this position in 'jmpbuffer'. If any further code
    // calls an error-exit function like PsychErrorExit() or PsychErrorExitMsg() then
    // the corresponding longjmp() call in our mexErrMsgTxt() implementation (see top of file)
    // will unwind the stack and restore stack-state and CPU state to the saved values in
    // jmpbuffer --> We will end up at this setjmp() call again, with a cleaned up stack and
    // CPU state, but setjmp will return a non-zero error code, signaling the abnormal abortion.
    if (setjmp(jmpbuffer[recLevel]) != 0) {
        // PsychErrorExit() or friends called! The CPU and stack are restored to a sane state.
        // Call our cleanup-routine to release memory that is PsychMallocTemp()'ed and to other
        // error-handling...
        goto PythonFunctionCleanup;
    }

    nrhsGLUE[recLevel] = nrhs;
    for (i = 0; i < nrhs; i++) {
        tmparg = PyTuple_GetItem(args, i);
        prhsGLUE[recLevel][i] = tmparg;

        // Empty args, strings and structs are special - handled directly the Python way.
        // Everything else goes through NumPy C-Interfaces:
        if ((tmparg == NULL) || (tmparg == Py_None) || mxIsChar(tmparg) || mxIsStruct(tmparg)) {
            // This object is already in our desired format:
            prhsNeedsConversion[recLevel][i] = FALSE;
        }
        else {
            // This object needs to get converted into a NumPy array of a suitable format for us.
            // We will do that on first access by client via PsychGetInArgPyPtr(). We can't convert
            // here, as we don't know yet if we need to convert into C-Memory layout or Fortran layout,
            // so just note down the need for lazy conversion on first "true" access:
            prhsNeedsConversion[recLevel][i] = TRUE;

            // Bump refcount if the input object is already a NumPy array, as our PythonFunctionCleanup
            // would decref it on normal exit or error exit. If it isn't a NumPy array yet, then no
            // need to bump, as a lazy conversion later on first access will turn it into a NumPy
            // array and then do the bumping. Or the argument may never get evaluated, so it would stay
            // a non NumPy array and don't need the bump/unbump treatment:
            if (PyArray_Check(prhsGLUE[recLevel][i]))
                Py_INCREF(prhsGLUE[recLevel][i]);
        }
    }

    // Set number of output arguments to "unknown" == -1, as we don't know yet:
    nlhsGLUE[recLevel] = -1;

    // NULL-out our pointer array of return value pointers plhsGLUE[recLevel]:
    memset(&plhsGLUE[recLevel][0], 0, sizeof(plhsGLUE[recLevel]));

    baseFunctionInvoked[recLevel] = FALSE;

    // If no subfunctions have been registered by the project then just invoke the project base function
    // If one of those has been registered.
    if (!PsychAreSubfunctionsEnabled()) {
        baseFunction = PsychGetProjectFunction(NULL);
        if (baseFunction != NULL) {
            baseFunctionInvoked[recLevel] = TRUE;
            (*baseFunction)();  //invoke the unnamed function
        } else
            PsychErrorExitMsg(PsychError_internal, "Project base function invoked but no base function registered");
    } else { // Subfunctions are enabled so pull out the function name string and invoke it.
        // Assess the nature of first and second arguments for finding the name of the sub function.
        for (i = 0; i < 2; i++) {
            isArgThere[i] = (nrhs > i) && (prhsGLUE[recLevel][i]);
            if (isArgThere[i]) {
                tmparg = prhsGLUE[recLevel][i];
            } else {
                tmparg = NULL;
            }

            isArgEmptyMat[i] = isArgThere[i] ? PsychIsDefaultMat(tmparg) : FALSE;
            isArgText[i] = isArgThere[i] ? mxIsChar(tmparg) : FALSE;
            if (isArgText[i]) {
                if (mxGetString(tmparg, argString[i], sizeof(argString[i])))
                    PsychErrorExitMsg(PsychError_user, "Invalid subcommand argument passed in. Not a command name!");

                // Empty subfunction command strings map to "empty matrix", for Mex/Matlab/Octave compatibility:
                if (strlen(argString[i]) == 0)
                    isArgEmptyMat[i] = TRUE;

                // Only consider 2nd arg as subfunction if 1st arg isn't already a subfunction:
                if ((i == 0) || (!isArgFunction[0])) {
                    fArg[i] = PsychGetProjectFunction(argString[i]);
                }
                else fArg[i] = NULL; // 1st arg is subfunction, so 2nd arg can't be as well.
            }
            isArgFunction[i] = isArgText[i] ? fArg[i] != NULL : FALSE;
        }

        // Figure out which of the two arguments might be the function name and either invoke it or exit with error
        // if we can't find one.
        if (!isArgThere[0] && !isArgThere[1]) { //no arguments passed so execute the base function
            baseFunction = PsychGetProjectFunction(NULL);
            if (baseFunction != NULL) {
                baseFunctionInvoked[recLevel]=TRUE;
                (*baseFunction)();
            } else
                PsychErrorExitMsg(PsychError_unimplemented, "Project base function invoked but no base function registered");
        }
        else if (isArgEmptyMat[0] && !isArgThere[1]) {
            PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state A)");
        }
        else if (isArgEmptyMat[0] && isArgEmptyMat[1]) {
            PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state B)");
        }
        else if (isArgEmptyMat[0] && isArgText[1]) {
            if (isArgFunction[1]) {
                nameFirstGLUE[recLevel] = FALSE;
                (*(fArg[1]))();
            }
            else
                PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state C)");
        }
        else if (isArgEmptyMat[0] && !isArgText[1]) {
            PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state D)");
        }
        else if (isArgText[0] && !isArgThere[1]) {
            if (isArgFunction[0]) {
                nameFirstGLUE[recLevel] = TRUE;
                (*(fArg[0]))();
            } else {
                // When we receive a first argument  which is a string and it is not recognized as a function name then call the default function
                // first to hopefully print a synopsis on a subfunctions-enabled module, then abort with "Unknown subfunction name".
                baseFunction = PsychGetProjectFunction(NULL);
                if (baseFunction != NULL) {
                    baseFunctionInvoked[recLevel]=TRUE;
                    (*baseFunction)();
                    PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state E)");
                } else
                    PsychErrorExitMsg(PsychError_unimplemented, "Project base function invoked but no base function registered");
            }
        }
        else if (isArgText[0] && isArgEmptyMat[1]) {
            if (isArgFunction[0]) {
                nameFirstGLUE[recLevel] = TRUE;
                (*(fArg[0]))();
            }
            else
                PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state F)");
        }
        else if (isArgText[0] && isArgText[1]) {
            if (isArgFunction[0] && !isArgFunction[1]) { //the first argument is the function name
                nameFirstGLUE[recLevel] = TRUE;
                (*(fArg[0]))();
            }
            else if (!isArgFunction[0] && isArgFunction[1]) { //the second argument is the function name
                nameFirstGLUE[recLevel] = FALSE;
                (*(fArg[1]))();
            }
            else if (!isArgFunction[0] && !isArgFunction[1]) { //neither argument is a function name
                //PrintfExit("Invalid command (error state G)");
                baseFunction = PsychGetProjectFunction(NULL);
                if (baseFunction != NULL) {
                    baseFunctionInvoked[recLevel]=TRUE;
                    (*baseFunction)();
                } else
                    PsychErrorExitMsg(PsychError_unimplemented, "Project base function invoked but no base function registered");
            }
            else if (isArgFunction[0] && isArgFunction[1]) //both arguments are function names
                PsychErrorExitMsg(PsychError_user, "Passed two function names");
        }
        else if (isArgText[0] && !isArgText[1]) {
            if (isArgFunction[0]) {
                nameFirstGLUE[recLevel] = TRUE;
                (*(fArg[0]))();
            }
            else
                PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state H)");
        }
        else if (!isArgText[0] && !isArgThere[1]) {  //this was modified for MODULEVersion with WaitSecs.
            //PrintfExit("Invalid command (error state H)");
            baseFunction = PsychGetProjectFunction(NULL);
            if (baseFunction != NULL) {
                baseFunctionInvoked[recLevel]=TRUE;
                (*baseFunction)();  //invoke the unnamed function
            } else
                PsychErrorExitMsg(PsychError_unimplemented, "Project base function invoked but no base function registered");
        }
        else if (!isArgText[0] && isArgEmptyMat[1]) {
            PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state I)");
        }
        else if (!isArgText[0] && isArgText[1]) {
            if (isArgFunction[1]) {
                nameFirstGLUE[recLevel] = FALSE;
                (*(fArg[1]))();
            } else
                PsychErrorExitMsg(PsychError_user, "Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state J)");
        }
        else if (!isArgText[0] && !isArgText[1]) {  //this was modified for Priority.
            //PrintfExit("Invalid command (error state K)");
            baseFunction = PsychGetProjectFunction(NULL);
            if (baseFunction != NULL) {
                baseFunctionInvoked[recLevel]=TRUE;
                (*baseFunction)();  //invoke the unnamed function
            } else
                PsychErrorExitMsg(PsychError_unimplemented, "Project base function invoked but no base function registered");
        }
    } //close else

    // If we reach this point of execution, then we're successfully done with function execution
    // and just need to return return arguments and clean up:
    if (psych_refcount_debug) {
        for (i = 0; i < MAX_OUTPUT_ARGS; i++) {
            if (plhsGLUE[recLevel][i] && (Py_REFCNT(plhsGLUE[recLevel][i]) >= psych_refcount_debug))
                printf("PTB-DEBUG: At non-error exit of PsychScriptingGluePythonDispatch: Refcount of plhsGLUE[recLevel %i][arg %i] = %li.\n",
                       recLevel, i, Py_REFCNT(plhsGLUE[recLevel][i]));
        }
    }

    // Find the true number of arguments to return in the return tuple:
    if (nlhsGLUE[recLevel] < 0) {
        for (i = 0; i < MAX_OUTPUT_ARGS; i++) {
            if (plhsGLUE[recLevel][i])
                nlhsGLUE[recLevel] = i + 1;
        }
    }

    // Multi-value return?
    if (nlhsGLUE[recLevel] > 1) {
        // Create an output tuple of suitable size:
        plhs = PyTuple_New((Py_ssize_t) nlhsGLUE[recLevel]);
        if (NULL == plhs)
            PsychErrorExitMsg(PsychError_internal, "PTB-CRITICAL: Failed to create output arg return tuple!!\n");

        // "Copy" our return values into the output tuple: If nlhs should be
        // zero (Python-Script does not expect any return arguments), but our
        // subfunction has assigned a return argument in slot 0 anyway, then
        // we return that argument and release our own temp-memory. This
        // provides "Matlab"-semantic, where a first unsolicited return argument
        // is printed anyway to the console for diagnostic purpose:
        for (i = 0; (i == 0 && plhsGLUE[recLevel][0] != NULL) || (i < nlhsGLUE[recLevel]); i++) {
            if (plhsGLUE[recLevel][i]) {
                // Assign return argument to proper slot of tuple:
                if (PyTuple_SetItem(plhs, (Py_ssize_t) i, plhsGLUE[recLevel][i])) {
                    printf("PTB-CRITICAL: Could not insert return argument for slot %i of output tuple!\n", i);
                    PsychErrorExitMsg(PsychError_internal, "PTB-CRITICAL: PyTuple_SetItem() failed.\n");
                }

                // NULL-out the array slot, only the output plhs tuple has a reference to
                // the output PyObject argument in slot i:
                plhsGLUE[recLevel][i] = NULL;
            }
            else {
                printf("PTB-DEBUG: Return argument for slot %i of output tuple not defined!\n", i);

                // Ref and assign empty return argument to proper slot of tuple:
                Py_INCREF(Py_None);
                if (PyTuple_SetItem(plhs, (Py_ssize_t) i, Py_None)) {
                    printf("PTB-CRITICAL: Could not insert return argument for slot %i of output tuple!\n", i);
                    PsychErrorExitMsg(PsychError_internal, "PTB-CRITICAL: PyTuple_SetItem() failed.\n");
                }
            }
        }
    }
    else if ((nlhsGLUE[recLevel] == 1) && (plhsGLUE[recLevel][0] != NULL)) {
        // Single return argument:
        plhs = plhsGLUE[recLevel][0];
        plhsGLUE[recLevel][0] = NULL;
    }
    else {
        // No return value at all:
        Py_INCREF(Py_None);
        plhs = Py_None;
    }

PythonFunctionCleanup:
    // The following code is executed both at end of normal execution, and also
    // during an error return. It has to do the common cleanup work:

    // Release references to NumPy PyArrays, as the PyObject -> PyArray code always
    // returns a new reference which we should get rid off, now that we don't need
    // it anymore:
    for (i = 0; i < nrhs; i++) {
        if (PyArray_Check(prhsGLUE[recLevel][i]))
            Py_XDECREF(prhsGLUE[recLevel][i]);

        prhsGLUE[recLevel][i] = NULL;
    }

    // Release "orphaned" output arguments that haven't been returned to the interpreter,
    // e.g., because of a PythonFunctionCleanup - error return:
    for (i = 0; i < MAX_OUTPUT_ARGS; i++) {
        if (psych_refcount_debug && plhsGLUE[recLevel][i] && (Py_REFCNT(plhsGLUE[recLevel][i]) >= psych_refcount_debug))
            printf("PTB-DEBUG: Orphaned output argument at cleanup: Refcount of plhsGLUE[recLevel %i][arg %i] = %li --> unref --> %li.\n",
                   recLevel, i, Py_REFCNT(plhsGLUE[recLevel][i]), Py_REFCNT(plhsGLUE[recLevel][i]) - 1);

        Py_XDECREF(plhsGLUE[recLevel][i]);
        plhsGLUE[recLevel][i] = NULL;
    }

    // Release all memory allocated via PsychMallocTemp():
    PsychFreeAllTempMemory();

    // Reset to not using C memory layout, but classic (backwards compatible) Fortran layout:
    use_C_memory_layout[recLevel] = FALSE;

    // Done with this call recursion level:
    PsychExitRecursion();

    // Return PyObject tuple with all return arguments:
    return(plhs);
}


/*  Call PsychExitGlue():
 *  Needed to safely reset modules in Python.
 */
PsychError PsychExitPythonGlue(void)
{
    // Debug output:
    if (DEBUG_PTBPYTHONGLUE)
        printf("PTB-INFO: Jettisoning submodule %s ...\n", PsychGetModuleName());

    // Mark ourselves as not yet initialized:
    firstTime = TRUE;

    // Call our regular exit routines to clean up and release all ressources:
    PsychErrorExitMsg(PsychExit(), NULL);

    // Done. Return control to Python:
    return(PsychError_none);
}


// PsychPyArgGet() helper for PsychGetInArgPyPtr() aka PsychGetInArgPtr():
// Does lazy, on-demand, one-time conversion of numeric data types to
// corresponding NumPy array data types of a suitable memory layout for
// our purpose. Adjusts reference counts of old/new instance of input arg
// accordingly and updates prhsGLUE cached copy. Does Fortran -> C conversion
// as needed.
static PyObject* PsychPyArgGet(int position)
{
    PyObject *ret = prhsGLUE[recLevel][position];

    if (psych_refcount_debug && (Py_REFCNT(prhsGLUE[recLevel][position]) >= psych_refcount_debug))
        printf("PTB-DEBUG:%s:PsychPyArgGet: Before convert: Refcount of prhsGLUE[recLevel %i][arg %i] = %li.\n",
               PsychGetFunctionName(), recLevel, position, Py_REFCNT(prhsGLUE[recLevel][position]));

    // Does this input argument need conversion to a NumPy array of suitable format?
    if (prhsNeedsConversion[recLevel][position]) {
        // Yes: Reset "needs conversion" flag and do the one-time conversion:
        prhsNeedsConversion[recLevel][position] = FALSE;

        // Convert it, either into C memory layout or Fortran memory layout. This
        // gives us a *new* reference to a NumPy array in *any* case, even if it was
        // an identity assignment because in-ret was already a suitable NumPy
        // array - in that case the refcount of original prhsGLUE got bumped by one.
        ret = PyArray_FROM_OF(ret, ((use_C_memory_layout[recLevel]) ? NPY_ARRAY_IN_ARRAY : NPY_ARRAY_IN_FARRAY));

        // If prhsGLUE was already a NumPy array, then its refcount got bumped in
        // the initial assignment code in PsychScriptingGluePythonDispatch(), so
        // need to undo the bump here, now that we have a new 'ret' reference,
        // before we assign a potentially new ret as the new prhsGLUE:
        if (PyArray_Check(prhsGLUE[recLevel][position]))
            Py_DECREF(prhsGLUE[recLevel][position]);

        if (DEBUG_PTBPYTHONGLUE || psych_refcount_debug)
            printf("PTB-DEBUG:%s:PsychPyArgGet: Arg %i: Conversion to NumPy array of %s triggered [refcount now %li]: %s\n",
                   PsychGetFunctionName(),
                   position, use_C_memory_layout[recLevel] ? "C layout" : "Fortran layout", Py_REFCNT(ret),
                   (ret == prhsGLUE[recLevel][position]) ? "No-Op passthrough." : "New object.");

        // Now that the new ret is a NumPy array in the wanted final format,
        // assign it, so we do not need to repeat the conversion on future access:
        prhsGLUE[recLevel][position] = ret;

        // At the end of this ballet, if this was a no-op conversion, then the
        // reference count of prhsGLUE should be unchanged. If it was a real
        // conversion, then only the Python interpreter should hold references
        // original prhsGLUE input argument aka the 'args' input tuple argument
        // to the PsychScriptingGluePythonDispatch(), and our new updated prhsGLUE,
        // aka ret is a new refcount == 1 NumPy array of suitable format and memory
        // layout, which we can now use and will dispose of in PythonFunctionCleanup
        // by dropping its refcount by one to zero.
    }

    if (psych_refcount_debug && (Py_REFCNT(prhsGLUE[recLevel][position]) >= psych_refcount_debug))
        printf("PTB-DEBUG:%s:PsychPyArgGet: After  convert: Refcount of prhsGLUE[recLevel %i][arg %i] = %li.\n",
                PsychGetFunctionName(), recLevel, position, Py_REFCNT(prhsGLUE[recLevel][position]));

    return(ret);
}


/*
 *    Return the PyObject pointer to the specified position. Note that we have some special rules for
 *    numbering the positions:
 *
 *    0 - This is always the command string or NULL if the project does not register a
 *        dispatch function and does accept subcommands.  If the function does accept sub
 *        commands, those may be passed in either the first or second position, but
 *        PsychGetArgPtr() will always return the command as the 0th.
 *
 *    1 - This is the first argument among the arguments which are not the subfunction name itself.
 *        It can occur in either the first or second position of the argument list, depending on
 *        in which of those two positions the function name itself appears.
 *
 *    2.. These positions are numbered correctly
 *
 *    TO DO:
 *
 *    2 - this function should be used by the one which gets the function name.
 *
 *
 *    Arguments are numbered 0..n.
 *
 *        - The 0th argument is a pointer to the PyObject holding
 *          the subfunction name string if we are in subfunction mode.
 *
 *        - The 0th argument is undefined if not in subfunction mode.
 *
 *        - The 1st argument is the argument of the 1st and 2nd which is not
 *          the subfunction name if in subfunction mode.
 *
 *        - The 1st argument is the first argument if not in subfunction mode.
 *
 *        - The 2nd-nth arguments are always the 2nd-nth arguments.
 */
const PyObject *PsychGetInArgPyPtr(int position)
{
    if (PsychAreSubfunctionsEnabled() && !baseFunctionInvoked[recLevel]) { //when in subfunction mode
        if (position < nrhsGLUE[recLevel]) { //an argument was passed in the correct position.
            if (position == 0) { //caller wants the function name argument.
                if (nameFirstGLUE[recLevel])
                    return(PsychPyArgGet(0));
                else
                    return(PsychPyArgGet(1));
            } else if (position == 1) { //they want the "first" argument.
                if (nameFirstGLUE[recLevel])
                    return(PsychPyArgGet(1));
                else
                    return(PsychPyArgGet(0));
            } else
                return(PsychPyArgGet(position));
        } else
            return(NULL);
    } else { //when not in subfunction mode and the base function is not invoked.
        if (position <= nrhsGLUE[recLevel])
            return(PsychPyArgGet(position-1));
        else
            return(NULL);
    }
}


PyObject **PsychGetOutArgPyPtr(int position)
{
    // Output argument PyObject
    if ((position == 1) ||
        ((position > 0) && ((position <= nlhsGLUE[recLevel]) || (nlhsGLUE[recLevel] == -1)))) {
        return(&(plhsGLUE[recLevel][position-1]));
    } else {
        printf("PTB-CRITICAL: PsychGetOutArgPyPtr() invalid position %i referenced [nlhs=%i], returning NULL!\n", position, nlhsGLUE[recLevel]);
        return(NULL);
    }
}


const PsychGenericScriptType *PsychGetInArgPtr(int position)
{
    return((const PsychGenericScriptType*) PsychGetInArgPyPtr(position));
}


/* PsychCheckSizeLimits(size_t m, size_t n, size_t p)
 *
 * Makes sure matrix/vector dimensions stay within the limits imposed
 * by a specific Psychtoolbox build for a specific 32 bit or 64 bit
 * version of Python.
 *
 * 32 bit builds are always limited to 2^31-1 elements per dimension at
 * most. In reality the limit is much lower, this is an upper bound.
 *
 * 64 bit builds may or may not allow bigger limits, depending if
 * the Python runtime supports a 64 bit ptbSize type. Even
 * then there are some limits imposed by available system memory and
 * processor architecture specific limits, e.g., many cpu's only
 * truly support 48 bit of memory, not the true 64 bit.
 *
 */
void PsychCheckSizeLimits(psych_int64 m, psych_int64 n, psych_int64 p)
{
    // No problem if ptbSize has capacity to contain size_t:
    if (sizeof(ptbSize) >= sizeof(size_t)) return;

    // Mismatch: ptbSize too small to contain size_t. This means
    // ptbSize is only defined as a 32-bit int on a system where
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
        printf("PTB-ERROR: You'd need to use a Psychtoolbox for 64-bit Python on a 64-Bit\n");
        printf("PTB-ERROR: operating system to get rid of this limit.\n");
    }
    else {
        printf("PTB-ERROR: This is a limitation of your version of Python.\n");
        printf("PTB-ERROR: You'd need to use a Psychtoolbox for 64-bit Python\n");
        printf("PTB-ERROR: on a 64-bit operating system to get rid of this limit.\n");
    }
    PsychErrorExitMsg(PsychError_user, "One of the dimensions of a returned matrix or vector exceeds 2^31-1 elements. This is not supported on your setup!");
}


/*
 *    mxCreateDoubleMatrix3D()
 *
 *    Create a 2D or 3D matrix of doubles.
 *
 *    Requirements are that m>0, n>0, p>=0.
 */
PyObject *mxCreateDoubleMatrix3D(psych_int64 m, psych_int64 n, psych_int64 p)
{
    int numDims;
    ptbSize dimArray[3];

    if (m == 0 || n == 0) {
        dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m, n, p);
        dimArray[0] = (ptbSize) m; dimArray[1] = (ptbSize) n; dimArray[2] = (ptbSize) p;
    }

    numDims = (p==0 || p==1) ? 2 : 3;

    return(mxCreateNumericArray(numDims, (ptbSize*) dimArray, PsychArgType_double));
}


/*
 *    mxCreateFloatMatrix3D()
 *
 *    Create a 2D or 3D matrix of floats.
 *
 *    Requirements are that m>0, n>0, p>=0.
 */
PyObject *mxCreateFloatMatrix3D(size_t m, size_t n, size_t p)
{
    int numDims;
    ptbSize dimArray[3];

    if (m == 0 || n == 0) {
        dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m, n, p);
        dimArray[0] = (ptbSize) m; dimArray[1] = (ptbSize) n; dimArray[2] = (ptbSize) p;
    }

    numDims = (p==0 || p==1) ? 2 : 3;

    return(mxCreateNumericArray(numDims, (ptbSize*) dimArray, PsychArgType_single));
}


/*
 *    mxCreateNativeBooleanMatrix3D()
 *
 *    Create a 2D or 3D matrix of native psych_bool types.
 *
 *    Requirements are that m>0, n>0, p>=0.
 */
PyObject *mxCreateNativeBooleanMatrix3D(size_t m, size_t n, size_t p)
{
    int numDims;
    ptbSize dimArray[3];

    if (m == 0 || n == 0) {
        dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m, n, p);
        dimArray[0] = (ptbSize) m; dimArray[1] = (ptbSize) n; dimArray[2] = (ptbSize) p;
    }

    numDims = (p==0 || p==1) ? 2 : 3;

    return(mxCreateNumericArray(numDims, (ptbSize*) dimArray, PsychArgType_boolean));
}


/*
 *    Create a 2D or 3D matrix of ubytes.
 *
 *    Requirements are that m>0, n>0, p>=0.
 */
PyObject *mxCreateByteMatrix3D(size_t m, size_t n, size_t p)
{
    int numDims;
    ptbSize dimArray[3];

    if (m == 0 || n == 0) {
        dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m, n, p);
        dimArray[0] = (ptbSize) m; dimArray[1] = (ptbSize) n; dimArray[2] = (ptbSize) p;
    }

    numDims = (p==0 || p==1) ? 2 : 3;

    return(mxCreateNumericArray(numDims, (ptbSize*) dimArray, PsychArgType_uint8));
}


static PyObject* PyExc[PsychError_last + 1] = { 0 };

/*
 * Python implementation of Python specific error handling.
 * Set proper Python exception state.
 */
void PsychProcessErrorInScripting(PsychError error, const char* message)
{
    PyObject *exception;

    if (PyExc[PsychError_invalidArg_absent] == NULL) {
        PyExc[PsychError_none] =                                 NULL;

        PyExc[PsychError_invalidArg_absent] =                    PyExc_SyntaxError;
        PyExc[PsychError_invalidArg_extra] =                     PyExc_SyntaxError;
        PyExc[PsychError_invalidArg_type] =                      PyExc_TypeError;
        PyExc[PsychError_invalidArg_size] =                      PyExc_ValueError;

        PyExc[PsychError_extraInputArg] =                        PyExc_SyntaxError;
        PyExc[PsychError_missingInputArg] =                      PyExc_SyntaxError;
        PyExc[PsychError_extraOutputArg] =                       PyExc_SyntaxError;
        PyExc[PsychError_missingOutputArg] =                     PyExc_SyntaxError;

        PyExc[PsychError_toomanyWin] =                           PyExc_MemoryError;
        PyExc[PsychError_outofMemory] =                          PyExc_MemoryError;
        PyExc[PsychError_scumberNotWindex] =                     PyExc_ValueError;
        PyExc[PsychError_windexNotScumber] =                     PyExc_ValueError;
        PyExc[PsychError_invalidWindex] =                        PyExc_IndexError;
        PyExc[PsychError_invalidIntegerArg] =                    PyExc_ValueError;
        PyExc[PsychError_invalidScumber] =                       PyExc_IndexError;
        PyExc[PsychError_invalidNumdex] =                        PyExc_IndexError;
        PyExc[PsychError_invalidColorArg] =                      PyExc_ValueError;
        PyExc[PsychError_invalidDepthArg] =                      PyExc_ValueError;
        PyExc[PsychError_invalidRectArg] =                       PyExc_ValueError;
        PyExc[PsychError_invalidNumberBuffersArg] =              PyExc_ValueError;
        PyExc[PsychError_nullWinRecPntr] =                       PyExc_RuntimeError;
        PyExc[PsychError_registerLimit] =                        PyExc_MemoryError;
        PyExc[PsychError_registered] =                           PyExc_RuntimeError;
        PyExc[PsychError_longString] =                           PyExc_ValueError;
        PyExc[PsychError_longStringPassed] =                     PyExc_ValueError;
        PyExc[PsychError_unimplemented] =                        PyExc_NotImplementedError;
        PyExc[PsychError_internal] =                             PyExc_RuntimeError;
        PyExc[PsychError_invalidArgRef] =                        PyExc_ValueError;
        PyExc[PsychError_OpenGL] =                               PyExc_EnvironmentError;
        PyExc[PsychError_system] =                               PyExc_EnvironmentError;
        PyExc[PsychError_InvalidWindowRecord] =                  PyExc_ValueError;
        PyExc[PsychError_unsupportedVideoMode] =                 PyExc_ValueError;
        PyExc[PsychError_user] =                                 PyExc_Exception;
        PyExc[PsychError_unrecognizedPreferenceName] =           PyExc_NameError;
        PyExc[PsychError_unsupportedOS9Preference] =             PyExc_NameError;
        PyExc[PsychError_inputMatrixIllegalDimensionSize] =      PyExc_ValueError;
        PyExc[PsychError_stringOverrun] =                        PyExc_BufferError;
        PyExc[PsychErorr_argumentValueOutOfRange] =              PyExc_ValueError;
    };

    // No error? Clear exception/error state:
    if (error == PsychError_none) {
        PyErr_Clear();
        return;
    }

    // Trust that the Python C-API or our own code has already set
    // a proper exception and error state if any error state is set:
    if (PyErr_Occurred())
        return;

    if (error == PsychError_system) {
        if (errno) {
            PyErr_SetFromErrno(PyExc_OSError);
            errno = 0;
            return;
        }
        else {
            #if PSYCH_SYSTEM == PSYCH_WINDOWS
                PyErr_SetFromWindowsErr(0);
                return;
            #endif
        }
        // Otherwise fall through to generic handling...
    }

    // Map to half-way suitable exception type and return with or without
    // additional error message:
    exception = PyExc[error];
    if (message != NULL)
        PyErr_SetString(exception, message);
    else
        PyErr_SetNone(exception);

    return;
}

/*
 *    Print string s and return return control to the calling environment.
 */
void PsychErrMsgTxt(char *s)
{
    PsychGenericScriptType *pcontent = NULL;

    // Is this the Screen() module?
    if (strcmp(PsychGetModuleName(), "Screen") == 0) {
        // Yes. We directly call our close and cleanup routine:
        #ifdef PTBMODULE_Screen
            ScreenCloseAllWindows();
        #endif
    } else {
        // Nope. This is a module other than Screen. Try to call Screen('Close')
        // via scripting environemnt:
        if (PsychRuntimeGetVariablePtr("global", "Screen", &pcontent)) {
            // Is it a function wrapper? Then call it Octave-style:
            if (!strcmp(PyEval_GetFuncName(pcontent), "Screen") && !strcmp(PyEval_GetFuncDesc(pcontent), "()"))
                PsychRuntimeEvaluateString("Screen('CloseAll');");
            else // Nope. Is it a module? Then call it module-style:
                if (!strcmp(PyEval_GetFuncName(pcontent), "module"))
                    PsychRuntimeEvaluateString("Screen.Screen('CloseAll');");
                else
                    printf("PsychErrMsgTxt: Failed to call Screen('CloseAll') - Weird signature, not the Screen module?!?\n");
        }
    }

    // Call the error printing and error handling facilities:
    mexErrMsgTxt((s && (strlen(s) > 0)) ? s : "See error message printed above.");
}


/*
 *    classify the PyObject element format using Pyschtoolbox argument type names
 *
 */
static PsychArgFormatType PsychGetTypeFromPyPtr(const PyObject *ppyPtr)
{
    PsychArgFormatType format;

    // First check for "empty" default argument:
    if (PsychIsDefaultMat(ppyPtr))
        format = PsychArgType_default;
    // then for string, as we use Python strings here, not NumPy array strings.
    else if (mxIsChar(ppyPtr))
        format = PsychArgType_char;
    // then for struct, as we use Python structs here, not NumPy array strings.
    else if (mxIsStruct(ppyPtr))
        format = PsychArgType_structArray;
    // then everything else, safely assuming it is a NumPy array object:
    else if (mxIsUint8(ppyPtr))
        format = PsychArgType_uint8;
    else if (mxIsUint16(ppyPtr))
        format = PsychArgType_uint16;
    else if (mxIsUint32(ppyPtr))
        format = PsychArgType_uint32;
    else if (mxIsUint64(ppyPtr))
        format = PsychArgType_uint64;
    else if (mxIsInt8(ppyPtr))
        format = PsychArgType_int8;
    else if (mxIsInt16(ppyPtr))
        format = PsychArgType_int16;
    else if (mxIsInt32(ppyPtr))
        format = PsychArgType_int32;
    else if (mxIsInt64(ppyPtr))
        format = PsychArgType_int64;
    else if (mxIsDouble(ppyPtr))
        format = PsychArgType_double;
    else if (mxIsSingle(ppyPtr))
        format = PsychArgType_single;
    else if (mxIsCell(ppyPtr))
        format = PsychArgType_cellArray;
    else if (mxIsLogical(ppyPtr))
        format = PsychArgType_boolean;
    else
        format = PsychArgType_unclassified;

    return format;
}


/*
 *    PsychSetReceivedArgDescriptor()
 *
 *    Accept an argument number and direction value (input or output).  Examine the specified argument and fill in an argument
 *    descriptor struture.  Ask a retainer function to store the descriptor.
 *
 */
PsychError PsychSetReceivedArgDescriptor(int argNum, psych_bool allow64BitSizes, PsychArgDirectionType direction)
{
    PsychArgDescriptorType  d;
    int                     numNamedOutputs, numOutputs;
    const PyObject          *ppyPtr;

    d.position = argNum;
    d.direction = direction;
    if (direction == PsychArgIn) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(argNum);
        d.isThere = (ppyPtr && !PsychIsDefaultMat(ppyPtr)) ? kPsychArgPresent : kPsychArgAbsent;
        if (d.isThere == kPsychArgPresent) { //the argument is there so fill in the rest of the description
            d.numDims = (int) mxGetNumberOfDimensions(ppyPtr);

            // If the calling function doesn't allow 64 bit sized input argument dimensions, then we check if
            // the input has a size within the positive signed integer range, i.e., at most INT_MAX elements
            // per dimension. Functions which can handle bigger inputs need to declare this explicitely by
            // setting allow64BitSizes == TRUE:
            if ((!allow64BitSizes) && ((mxGetM(ppyPtr) >= INT_MAX) || (mxGetNOnly(ppyPtr) >= INT_MAX) || (mxGetP(ppyPtr) >= INT_MAX))) {
                printf("PTB-ERROR: %i. Input argument exceeds allowable maximum size of 2^31 - 1 elements\n", argNum);
                printf("PTB-ERROR: in at least one dimension. Psychtoolbox can't handle such huge matrices or vectors.\n");
                PsychErrorExitMsg(PsychError_user, "Input argument exceeds maximum supported count of 2^31 - 1 elements!");
            }

            d.mDimMin = d.mDimMax = (psych_int64) mxGetM(ppyPtr);
            d.nDimMin = d.nDimMax = (psych_int64) mxGetNOnly(ppyPtr);
            d.pDimMin = d.pDimMax = (psych_int64) mxGetP(ppyPtr);
            d.type = PsychGetTypeFromPyPtr(ppyPtr);
        }
    }
    else { //(direction == PsychArgOut)
        numNamedOutputs = PsychGetNumNamedOutputArgs();
        numOutputs = PsychGetNumOutputArgs();
        if (numNamedOutputs >=argNum)
            d.isThere = kPsychArgPresent;
        else if (numOutputs >=argNum)
            d.isThere = kPsychArgFixed;
        else
            d.isThere = kPsychArgAbsent;
    }

    PsychStoreArgDescriptor(NULL,&d);
    return(PsychError_none);
}


//local function definitions for ScriptingGluePython.c
//___________________________________________________________________________________________


psych_bool PsychIsDefaultMat(const PyObject *mat)
{
    return ((mat == Py_None) ||
            (PyList_Check(mat) && (PyList_Size((PyObject *) mat) == 0)) ||
            (PyArray_Check(mat) && ((PyArray_SIZE((PyArrayObject*) mat) == 0) || (PyArray_IsZeroDim(mat) && !PyArray_CheckScalar(mat)))));
}


//functions for project access to module call arguments
//___________________________________________________________________________________________


//functions which query the number and nature of supplied arguments

/*
 *    PsychGetNumInputArgs()
 *
 *    -The count excludes the command argument and includes ALL arguments supplied, including
 *    default arguments.
 *
 *    -For the time being, the only way to check if all required arguments are supplied in the
 *    general case of mixed required and optional arguments is to check each individually. Probably
 *    the best way to to fix this is to employ a description of which are required and which optional
 *    and compare that against what was passed to the subfunction.
 */
int PsychGetNumInputArgs(void)
{
    if (PsychAreSubfunctionsEnabled() && !baseFunctionInvoked[recLevel]) //this should probably be just baseFunctionInvoked[recLevel] wo PsychSubfunctionEnabled.
        return(nrhsGLUE[recLevel]-1);
    else
        return(nrhsGLUE[recLevel]);
}


// Careful! Here's a catch: I don't know of a straightforward way to actually
// find out how many return arguments (output args) the calling Python runtime
// expects from us. Apparently that's already a major problem in the Python
// language itself for native Python scripting code.
int PsychGetNumOutputArgs(void)
{
    return((nlhsGLUE[recLevel] <= 0) ? ((nlhsGLUE[recLevel] == 0) ? 1 : MAX_OUTPUT_ARGS) : nlhsGLUE[recLevel]);
}


int PsychGetNumNamedOutputArgs(void)
{
    return(nlhsGLUE[recLevel]);
}


PsychError PsychCapNumOutputArgs(int maxNamedOutputs)
{
    // Number of output args already known?
    if ((PsychGetNumNamedOutputArgs() < 0) && (maxNamedOutputs < MAX_OUTPUT_ARGS)) {
        // No. We use the maxNamedOutputs passed in here as the clue to what the
        // correct value is, given that at least the Python 2 api may be incapable
        // of easily telling us this magic value:
        // TODO FIXME Don't use this as hint, always detect true number at function exit, always succeed: nlhsGLUE[recLevel] = maxNamedOutputs;
        return(PsychError_none);
    }

    // Maximum number of expected outputs known, so enforce them:
    if (PsychGetNumNamedOutputArgs() > maxNamedOutputs)
        return(PsychError_extraOutputArg);
    else
        return(PsychError_none);
}


/*
 *    The argument is not present if a default m*n=0 matrix was supplied, '' or []
 */
psych_bool PsychIsArgPresent(PsychArgDirectionType direction, int position)
{
    int numArgs;

    if (direction==PsychArgOut) {
        return((psych_bool)(PsychGetNumOutputArgs()>=position));
    } else {
        if ((numArgs=PsychGetNumInputArgs())>=position)
            return(!(PsychIsDefaultMat((PyObject*) PsychGetInArgPyPtr(position)))); //check if its default
        else
            return(FALSE);
    }
}


PsychArgFormatType PsychGetArgType(int position) //this is for inputs because outputs are unspecified by the calling environment.
{
    if (!(PsychIsArgReallyPresent(PsychArgIn, position)))
        return(PsychArgType_none);

    return(PsychGetTypeFromPyPtr((PyObject*) PsychGetInArgPyPtr(position)));
}


size_t PsychGetArgM(int position)
{
    if (!(PsychIsArgPresent(PsychArgIn, position)))
        PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
    return( mxGetM((PyObject*) PsychGetInArgPyPtr(position)));
}


size_t PsychGetArgN(int position)
{
    if (!(PsychIsArgPresent(PsychArgIn, position)))
        PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
    return( mxGetNOnly((PyObject*) PsychGetInArgPyPtr(position)));
}


size_t PsychGetArgP(int position)
{
    if (!(PsychIsArgPresent(PsychArgIn, position)))
        PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
    return( mxGetP((PyObject*) PsychGetInArgPyPtr(position)));
}


/*functions which output arguments.
 * ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 * ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 * -Naming scheme:
 *    -Outputing return arguments:
 *        - "PsychAllocOut*Arg" : allocate and set a pointer to volatile memory to be filled with returned information by the caller.
 *        - "PsychCopyOut*Arg : accept a pointer to ouput values and fill in the return matrix memory with the values.
 *    -Reading input arguments:
 *        - "PsychAllocIn*Arg" : set a pointer to volatile memory allocated by "PsychAllocIn*Arg" and holding the input value.
 *        - "PsychCopyIn*Arg" : accept a pointer to memory and fill in that memory with the input argument values.
 *
 * -These all supply their own dynamic memory now, even functions which return arguments, and, in the case of
 * Put functions,  even when those arguments are not present !  If you don't want the function to go allocating
 * memory for an unsupplied return argument, detect the presense of that argument from within your script and
 * conditionally invoke PsychPut*Arg.  This is a feature which allows you to ignore the presense of a return
 * argument in the case where memory which holds a return argument serves other purposes.
 *
 * -All dynamic memory provided by these functions is volatile, that is, it is lost when the mex module exits and
 * returns control to the scripting environemnt.  To make it non volatile, call Psych??? on it.
 *
 * ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 */

psych_bool PsychCopyOutDoubleArg(int position, PsychArgRequirementType isRequired, double value)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double,  isRequired, 1, 1, 1, 1, 0, 0);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = PyFloat_FromDouble(value);
    }
    return(putOut);
}

psych_bool PsychAllocOutDoubleArg(int position, PsychArgRequirementType isRequired, double **value)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, 1, 1, 1, 1, 0, 0);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateDoubleMatrix3D(1, 1, 0);
        *value = mxGetData(*mxpp);
    } else {
        *value = (double *) mxMalloc(sizeof(double));
    }

    return(putOut);
}


/*
 * PsychAllocOutDoubleMatArg()
 *
 * A)return argument mandatory:
 *    1)return argument not present:     exit with an error.
 *    2)return argument present:         allocate an output matrix and set return arg pointer. Set *array to the array within the new matrix. Return TRUE.
 * B)return argument optional:
 *    1)return argument not present:     return FALSE to indicate absent return argument.  Create an array.   Set *array to the new array.
 *    2)return argument present:         allocate an output matrix and set return arg. pointer. Set *array to the array within the new matrix.  Return TRUE.
 */
psych_bool PsychAllocOutDoubleMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, double **array)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, m, m, n, n, p, p);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateDoubleMatrix3D(m, n, p);
        *array = (double*) mxGetData(*mxpp);
    } else
        *array = (double*) mxMalloc(sizeof(double) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));

    return(putOut);
}


/*
 * PsychAllocOutFloatMatArg()
 *
 * This function allocates out a matrix of single precision floating point type,
 * that is C data type 32-bit float.
 *
 * A)return argument mandatory:
 *    1)return argument not present:         exit with an error.
 *    2)return argument present:         allocate an output matrix and set return arg pointer. Set *array to the array within the new matrix. Return TRUE.
 * B)return argument optional:
 *    1)return argument not present:      return FALSE to indicate absent return argument.  Create an array.   Set *array to the new array.
 *    2)return argument present:         allocate an output matrix and set return arg. pointer. Set *array to the array within the new matrix.  Return TRUE.
 */
psych_bool PsychAllocOutFloatMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, float **array)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_single, isRequired, m, m, n, n, p, p);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateFloatMatrix3D((size_t) m, (size_t) n, (size_t) p);
        *array = (float*) mxGetData(*mxpp);
    } else
        *array = (float*) mxMalloc(sizeof(float) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));

    return(putOut);
}


/*
 *    PsychCopyOutBooleanArg()
 */
psych_bool PsychCopyOutBooleanArg(int position, PsychArgRequirementType isRequired, PsychNativeBooleanType value)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_boolean, isRequired, 1, 1, 1, 1, 0, 0);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = PyBool_FromLong((long) value);
    }

    return(putOut);
}


/*
 *    PsychAllocOutBooleanMatArg()
 *
 *    A)return argument mandatory:
 *    1)return argument not present:         exit with an error.
 *    2)return argument present:         allocate an output matrix and set return arg pointer. Set *array to the array within the new matrix. Return TRUE.
 *    B)return argument optional:
 *    1)return argument not present:      return FALSE to indicate absent return argument.  Create an array.   Set *array to the new array.
 *    2)return argument present:         allocate an output matrix and set return arg. pointer. Set *array to the array within the new matrix.  Return TRUE.
 */
psych_bool PsychAllocOutBooleanMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, PsychNativeBooleanType **array)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_boolean, isRequired, m, m, n, n, p, p);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateNativeBooleanMatrix3D((size_t) m, (size_t) n, (size_t) p);
        *array = (PsychNativeBooleanType *) mxGetLogicals(*mxpp);
    } else {
        *array = (PsychNativeBooleanType *) mxMalloc(sizeof(PsychNativeBooleanType) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    }

    return(putOut);
}


/*
 *    PsychAllocOutUnsignedByteMatArg()
 *
 *    Like PsychAllocOutDoubleMatArg() execept for unsigned bytes instead of doubles.
 */
psych_bool PsychAllocOutUnsignedByteMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, psych_uint8 **array)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_uint8, isRequired, m, m, n, n, p, p);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateByteMatrix3D((size_t) m, (size_t) n, (size_t) p);
        *array = (psych_uint8 *) mxGetData(*mxpp);
    } else {
        *array = (psych_uint8 *) mxMalloc(sizeof(psych_uint8) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    }

    return(putOut);
}


/*
 *    PsychAllocOutUnsignedInt16MatArg()
 *
 *    Like PsychAllocOutDoubleMatArg() execept for unsigned shorts instead of doubles.
 */
psych_bool PsychAllocOutUnsignedInt16MatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, psych_uint16 **array)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;
    ptbSize         dimArray[3];
    int             numDims;

    // Compute output array dimensions:
    if (m <= 0 || n <= 0) {
        dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;  //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m, n, p);
        dimArray[0] = (ptbSize) m; dimArray[1] = (ptbSize) n; dimArray[2] = (ptbSize) p;
    }

    numDims = (p == 0 || p == 1) ? 2 : 3;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_uint16, isRequired, m, m, n, n, p, p);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateNumericArray(numDims, (ptbSize*) dimArray, PsychArgType_uint16);
        *array = (psych_uint16 *) mxGetData(*mxpp);
    } else {
        *array = (psych_uint16 *) mxMalloc(sizeof(psych_uint16) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    }

    return(putOut);
}


psych_bool PsychCopyOutDoubleMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, double *fromArray)
{
    PyObject    **mxpp;
    double      *toArray;
    PsychError  matchError;
    psych_bool  putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, m, m, n, n, p, p);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateDoubleMatrix3D(m, n, p);
        toArray = mxGetData(*mxpp);
        //copy the input array to the output array now
        memcpy(toArray, fromArray, sizeof(double) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    }

    return(putOut);
}


psych_bool PsychCopyOutUnsignedInt16MatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, psych_uint16 *fromArray)
{
    PyObject        **mxpp;
    psych_uint16    *toArray;
    PsychError      matchError;
    psych_bool      putOut;
    ptbSize         dimArray[3];
    int             numDims;

    // Compute output array dimensions:
    if (m <= 0 || n <= 0) {
        dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m, n, p);
        dimArray[0] = (ptbSize) m; dimArray[1] = (ptbSize) n; dimArray[2] = (ptbSize) p;
    }

    numDims = (p == 0 || p == 1) ? 2 : 3;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_uint16, isRequired, m, m, n, n, p, p);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateNumericArray(numDims, (ptbSize*) dimArray, PsychArgType_uint16);
        toArray = (psych_uint16*) mxGetData(*mxpp);

        //copy the input array to the output array now
        memcpy(toArray, fromArray, sizeof(psych_uint16) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    }

    return(putOut);
}


/*
 *    PsychCopyOutCharArg()
 *
 *    Accept a null terminated string and return it in the specified position.
 *
 */
psych_bool PsychCopyOutCharArg(int position, PsychArgRequirementType isRequired, const char *str)
{
    PyObject    **mxpp;
    PsychError  matchError;
    psych_bool  putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_char, isRequired, 0, strlen(str), 0, strlen(str), 0, 0);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = mxCreateString(str);
    }

    return(putOut);
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
psych_bool PsychAllocInDoubleMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, double **array)
{
    psych_int64 mb, nb, pb;
    psych_bool rc = PsychAllocInDoubleMatArg64(position, isRequired, &mb, &nb, &pb, array);
    *m = (int) mb;
    *n = (int) nb;
    *p = (int) pb;
    return(rc);
}


/* Alloc-in double matrix, but allow for 64-bit dimension specs. */
psych_bool PsychAllocInDoubleMatArg64(int position, PsychArgRequirementType isRequired, psych_int64 *m, psych_int64 *n, psych_int64 *p, double **array)
{
    const PyObject    *ppyPtr;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1, -1, 1, -1, 0, -1);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
        *m = (psych_int64) mxGetM(ppyPtr);
        *n = (psych_int64) mxGetNOnly(ppyPtr);
        *p = (psych_int64) mxGetP(ppyPtr);
        *array = mxGetData(ppyPtr);
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
 *
 * Alloc-In a single precision floating point matrix, i.e. a matrix of
 * C data type 32 bit float.
 * This function allows to alloc in matrices with more than 2^32 elements
 * per matrix dimension on 64 bit systems. Therefore the returned size
 * descriptors must be psych_int64 variables, not int variables or bad things
 * will happen.
 *
 * If the function receives a double() precision input matrix instead of the
 * expected single() precision matrix, it will automatically create a temporary
 * copy, with all values copied/converted from double to single data type, aka
 * double -> float cast. This is transparent to the caller, so it can always
 * operate on a returned float matrix -- at a performance penalty for the extra
 * copy of course.
 *
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
psych_bool PsychAllocInFloatMatArg64(int position, PsychArgRequirementType isRequired, psych_int64 *m, psych_int64 *n, psych_int64 *p, float **array)
{
    const PyObject  *ppyPtr;
    PsychError      matchError;
    psych_bool      acceptArg;
    double*         arrayD;
    float*          arrayF;
    psych_int64     i;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_single, isRequired, 1, -1, 1, -1, 0, -1);
    matchError = PsychMatchDescriptors();

    // Argument provided, but not of required float type?
    if (matchError == PsychError_invalidArg_type) {
        // If the input type is double precision floating point, then we convert
        // it here into single precision floating point via a temporary buffer.
        // This is used for functions which absolutely need float input, e.g.,
        // OpenGL-ES rendering code, but should accept double input from usercode
        // so usercode doesn not need to be specifically ported for OpenGL-ES platforms.
        // Performance may suffer somwehat though...
        PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1, -1, 1, -1, 0, -1);
        matchError = PsychMatchDescriptors();
        acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
        if (acceptArg) {
            ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
            *m = (psych_int64) mxGetM(ppyPtr);
            *n = (psych_int64) mxGetNOnly(ppyPtr);
            *p = (psych_int64) mxGetP(ppyPtr);

            // Get a double pointer to the double input data matrix:
            arrayD = (double*) mxGetData(ppyPtr);

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
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
        *m = (psych_int64) mxGetM(ppyPtr);
        *n = (psych_int64) mxGetNOnly(ppyPtr);
        *p = (psych_int64) mxGetP(ppyPtr);
        *array = (float*) mxGetData(ppyPtr);
    }
    return(acceptArg);
}


/*
 *    PsychAllocInIntegerListArg()
 *
 *    In a scriptiong language such as Python where numbers are almost always stored as doubles, this function is useful to check
 *    that the value input is an integer value stored within a double type.
 *
 *    Otherwise it just here to imitate the version written for other scripting languages.
 */
psych_bool PsychAllocInIntegerListArg(int position, PsychArgRequirementType isRequired, int *numElements, int **array)
{
    int                 m, n, p, i;
    double              *doubleMatrix;
    psych_int64         *int64Matrix;
    psych_bool          isThere;
    const PyObject      *ppyPtr;
    PsychError          matchError;
    psych_bool          acceptArg;

    // Try to get native integer matrix of int32:
    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_int32, isRequired, 1, -1, 1, -1, 0, -1);
    matchError = PsychMatchDescriptors();
    if (matchError == PsychError_none) {
        acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
        if (acceptArg) {
            ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
            m = (int) mxGetM(ppyPtr);
            n = (int) mxGetNOnly(ppyPtr);
            p = (int) mxGetP(ppyPtr);
            p = (p == 0) ? 1 : p;

            if ((psych_uint64) m * (psych_uint64) n * (psych_uint64) p >= INT_MAX) {
                printf("PTB-ERROR: %i th input argument has more than 2^31 - 1 elements! This is not supported.\n", position);
                return(FALSE);
            }

            *numElements = m * n * p;
            *array = (int*) mxGetData(ppyPtr);
        }
        return(acceptArg);
    }

    // Retry to get native integer matrix of int64:
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_int64, isRequired, 1, -1, 1, -1, 0, -1);
    matchError = PsychMatchDescriptors();
    if (matchError == PsychError_none) {
        acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
        if (acceptArg) {
            ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
            m = (int) mxGetM(ppyPtr);
            n = (int) mxGetNOnly(ppyPtr);
            p = (int) mxGetP(ppyPtr);
            p = (p == 0) ? 1 : p;

            if ((psych_uint64) m * (psych_uint64) n * (psych_uint64) p >= INT_MAX) {
                printf("PTB-ERROR: %i th input argument has more than 2^31 - 1 elements! This is not supported.\n", position);
                return(FALSE);
            }

            *numElements = m * n * p;
            int64Matrix = (psych_int64 *) mxGetData(ppyPtr);
            *array = (int*) mxMalloc((size_t) *numElements * sizeof(int));
            for (i = 0; i < *numElements; i++) {
                (*array)[i] = (int) int64Matrix[i];
            }
        }
        return(acceptArg);
    }

    // Nope. Try falling back to double matrix and convert to integers:
    isThere = PsychAllocInDoubleMatArg(position, isRequired, &m, &n, &p, &doubleMatrix);
    if (!isThere)
        return(FALSE);

    p = (p == 0) ? 1 : p;

    if ((psych_uint64) m * (psych_uint64) n * (psych_uint64) p >= INT_MAX) {
        printf("PTB-ERROR: %i th input argument has more than 2^31 - 1 elements! This is not supported.\n", position);
        return(FALSE);
    }

    *numElements = m * n * p;
    *array = (int*) mxMalloc((size_t) *numElements * sizeof(int));
    for (i = 0; i < *numElements; i++) {
        if (!PsychIsIntegerInDouble(doubleMatrix + i))
            PsychErrorExit(PsychError_invalidIntegerArg);

        (*array)[i] = (int) doubleMatrix[i];
    }

    return(TRUE);
}


/*
 *    PsychAllocInByteMatArg()
 *
 *    Like PsychAllocInDoubleMatArg() except it returns an array of unsigned bytes.
 */
psych_bool PsychAllocInUnsignedByteMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, unsigned char **array)
{
    const PyObject    *ppyPtr;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_uint8, isRequired, 1, -1, 1, -1, 0, -1);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
        *m = (int) mxGetM(ppyPtr);
        *n = (int) mxGetNOnly(ppyPtr);
        *p = (int) mxGetP(ppyPtr);
        *array = (unsigned char *) mxGetData(ppyPtr);
    }
    return(acceptArg);
}


/*
 *    PsychCopyInDoubleArg()
 *
 *    For 1x1 double.
 *
 *    Return in *value a double passed in the specified position, or signal an error if there is no
 *    double there and the argument is required, or don't change "value" if the argument is optional
 *    and none is supplied.
 *
 *    Note that if the argument is optional and ommitted PsychGetDoubleArg won't overwrite *value, allowing
 *    for specification of default values within project files without checking for their
 *    presense and conditinally filing in values.
 */
psych_bool PsychCopyInDoubleArg(int position, PsychArgRequirementType isRequired, double *value)
{
    PyObject        *ppyPtr;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double | PsychArgType_int32 | PsychArgType_int64, isRequired, 1,1,1,1,1,1);
    matchError = PsychMatchDescriptors();

    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
        *value = PyFloat_AsDouble(ppyPtr);

        if (PyErr_Occurred())
            PsychErrorExit(PsychError_invalidArg_type);
    }

    return(acceptArg);
}


/*
 *    Like PsychCopyInDoubleArg() with the additional restriction that the passed value must
 *    not have a fractional component and that it fits within the bounds of a C integer.
 *
 */
psych_bool PsychCopyInIntegerArg(int position, PsychArgRequirementType isRequired, int *value)
{
    PyObject        *ppyPtr;
    double          tempDouble;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double | PsychArgType_int32 | PsychArgType_int64, isRequired, 1,1,1,1,1,1);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);

        if (PyLong_Check(ppyPtr)) {
            *value = (int) PyLong_AsLong(ppyPtr);
            if (PyErr_Occurred())
                PsychErrorExit(PsychError_invalidIntegerArg);
        }
        else {
            tempDouble = PyFloat_AsDouble(ppyPtr);
            if (PyErr_Occurred() || !PsychIsIntegerInDouble(&tempDouble))
                PsychErrorExit(PsychError_invalidIntegerArg);
            *value = (int) tempDouble;
        }
    }

    return(acceptArg);
}


psych_bool PsychCopyInIntegerArg64(int position,  PsychArgRequirementType isRequired, psych_int64 *value)
{
    PyObject        *ppyPtr;
    double          tempDouble;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double | PsychArgType_int32 | PsychArgType_int64, isRequired, 1, 1, 1, 1, 1, 1);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);

        if (PyLong_Check(ppyPtr)) {
            *value = (psych_int64) PyLong_AsLongLong(ppyPtr);
            if (PyErr_Occurred())
                PsychErrorExit(PsychError_invalidIntegerArg);
        }
        else {
            tempDouble = PyFloat_AsDouble(ppyPtr);
            if (PyErr_Occurred() || !PsychIsInteger64InDouble(&tempDouble))
                PsychErrorExit(PsychError_invalidIntegerArg);
            *value = (psych_int64) tempDouble;
        }
    }

    return(acceptArg);
}


/*
 *    PsychAllocInDoubleArg()
 */
psych_bool PsychAllocInDoubleArg(int position, PsychArgRequirementType isRequired, double **value)
{
    int             m, n, p;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1, 1, 1, 1, 1, 1);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg)
        acceptArg = PsychAllocInDoubleMatArg(position, isRequired, &m, &n, &p, value);

    return(acceptArg);
}


/*
 *    PsychAllocInCharArg()
 *
 *    Reads in a string and sets *str to point to the string.
 *
 *    This function violates the rule for AllocIn fuctions that if the argument is optional and absent we allocate
 *    space. That turns out to be an unuseful feature anyway, so we should probably get ride of it.
 *
 *    The second argument has been modified to passively accept, without error, an argument in the specified position of non-character type.
 *
 *        0    kPsychArgOptional  Permit either an argument of the specified type or no argument.  An argument of any a different type is an error.
 *        1    kPsychArgRequired  Permit only an argument of the specifed type.  If no argument is present, exit with error.
 *        2    kPsychArgAnything  Permit any argument type without error, but only read the specified type.
 *
 */
psych_bool PsychAllocInCharArg(int position, PsychArgRequirementType isRequired, char **str)
{
    PyObject        *ppyPtr;
    int             status;
    psych_uint64    strLen;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_char, isRequired, 0, kPsychUnboundedArraySize ,0, kPsychUnboundedArraySize, 0 , 1);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
        if (PyUnicode_Check(ppyPtr))
            strLen = (psych_uint64) PyUnicode_GetLength(ppyPtr) + 1;
        else
            strLen = (psych_uint64) PyBytes_Size(ppyPtr) + 1;

        if (strLen >= INT_MAX)
            PsychErrorExitMsg(PsychError_user, "Tried to pass in a string with more than 2^31 - 1 characters. Unsupported!");

        *str = (char *) PsychCallocTemp((size_t) strLen, sizeof(char));
        status = mxGetString(ppyPtr, *str, (ptbSize) strLen);
        if (status != 0)
            PsychErrorExitMsg(PsychError_internal, "mxGetString failed to get the string");
    }
    return(acceptArg);
}


/*
 *    Get a psych_bool flag from the specified argument position. The type can be be psych_bool, uint8, or
 *    char.  If the numerical value is equal to zero or if its empty then the flag is FALSE, otherwise the
 *    flag is TRUE.
 *
 *    PsychGetFlagArg returns TRUE if the argument was present and false otherwise:
 *
 *    A- Argument required
 *        1- Argument is present: load *argVal and return TRUE
 *        2- Argument is absent: exit with an error
 *    B- Argument is optional
 *        1- Argument is present: load *argVal and return TRUE
 *        2- Argument is absent: leave *argVal alone and return FALSE
 *
 *    Note: if we modify PsychGetDoubleArg to accept all types and coerce them, then we could simplify by
 *    calling that instead of doing all of this stuff...
 *
 */
psych_bool PsychAllocInFlagArg(int position,  PsychArgRequirementType isRequired, psych_bool **argVal)
{
    psych_bool      acceptArg;
    psych_bool      value;

    acceptArg = PsychCopyInFlagArg(position, isRequired, &value);
    if (acceptArg) {
        *argVal = (psych_bool *) mxMalloc(sizeof(psych_bool));
        **argVal = value;
    }

    return(acceptArg);
}


/*
 *    PsychCopyInFlagArg()
 */
psych_bool PsychCopyInFlagArg(int position, PsychArgRequirementType isRequired, psych_bool *argVal)
{
    PyObject        *ppyPtr;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double | PsychArgType_int32 | PsychArgType_int64 | PsychArgType_boolean,
                                   isRequired, 1, 1, 1, 1, kPsychUnusedArrayDimension, kPsychUnusedArrayDimension);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
        if (PyBool_Check(ppyPtr)) {
            *argVal = (psych_bool) (ppyPtr == Py_True);
        }
        else if (PyArray_ISBOOL((const PyArrayObject*) ppyPtr)) {
            if (mxGetLogicals(ppyPtr)[0])
                *argVal = (psych_bool) 1;
            else
                *argVal = (psych_bool) 0;
        } else {
            if (mxGetScalar(ppyPtr))
                *argVal = (psych_bool) 1;
            else
                *argVal = (psych_bool) 0;
        }
    }

    return(acceptArg);
}


/* TODO FIXME */
psych_bool PsychAllocInFlagArgVector(int position,  PsychArgRequirementType isRequired, int *numElements, psych_bool **argVal)
{
    PyObject       *ppyPtr;
    PsychError     matchError;
    psych_bool     acceptArg;
    int            i;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType) (PsychArgType_double | PsychArgType_boolean),
                                   isRequired, 1, kPsychUnboundedArraySize, 1, kPsychUnboundedArraySize, kPsychUnusedArrayDimension, kPsychUnusedArrayDimension);
    matchError = PsychMatchDescriptors();
    acceptArg = PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
        if ((psych_uint64) mxGetM(ppyPtr) * (psych_uint64) mxGetN(ppyPtr) >= INT_MAX) {
            printf("PTB-ERROR: %i th input argument has more than 2^31 - 1 elements! This is not supported.\n", position);
            *numElements = 0;
            return(FALSE);
        }

        *numElements = (int) (mxGetM(ppyPtr) * mxGetN(ppyPtr));

        // Unlike other PsychAllocIn* functions, here we allocate new memory and copy the input to it rather than simply returning a pointer into the received array.
        // That's because we want the booleans returned to the caller by PsychAllocInFlagArgVector() to alwyas be psych_bool's, yet we accept as flags either 64-bit
        // doubles or native logical/boolean type.
        *argVal = (psych_bool*) mxMalloc(sizeof(psych_bool) * ((size_t) *numElements));
        for (i = 0; i < *numElements; i++) {
            if (mxIsLogical(ppyPtr)) {
                if (mxGetLogicals(ppyPtr)[i])
                    (*argVal)[i] = (psych_bool) 1;
                else
                    (*argVal)[i] = (psych_bool) 0;
            } else {
                if (((double*) mxGetData(ppyPtr))[i])
                    (*argVal)[i] = (psych_bool) 1;
                else
                    (*argVal)[i] = (psych_bool) 0;
            }
        }
    }
    return(acceptArg);    //the argument was not present (and optional).
}


psych_bool PsychCopyOutFlagArg(int position, PsychArgRequirementType isRequired, psych_bool argVal)
{
    return(PsychCopyOutDoubleArg(position, isRequired, (double) argVal));
}


// functions which allocate native types without assigning them to return arguments.
// this is useful for embedding native structures within each other.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 *    PsychAllocateNativeDoubleMat()
 *
 *    Create an opaque native double matrix.   Return both
 *        - Its handle,  which is specified when nesting the native matrix nesting withing other native types.
 *        - A handle to the C array of doubles enclosed by the native type.
 *
 *    In any case, *cArray will point to the C array of doubles enclosed by the native type in the end.
 *
 */
void PsychAllocateNativeDoubleMat(psych_int64 m, psych_int64 n, psych_int64 p, double **cArray, PsychGenericScriptType **nativeElement)
{
    double *cArrayTemp;

    PsychCheckSizeLimits(m, n, p);
    *nativeElement = mxCreateDoubleMatrix3D(m,n,p);
    cArrayTemp = mxGetData(*nativeElement);
    *cArray = cArrayTemp;
}


/*
 *    PsychAllocateNativeUnsignedByteMat()
 *
 *    Create an opaque native byte matrix.   Return both
 *        - Its handle,  which is specified when nesting the native matrix nesting withing other native types.
 *        - A handle to the C array of psych_uint8's enclosed by the native type.
 *
 *    In any case, *cArray will point to the C array of psych_uint8's enclosed by the native type in the end.
 *
 */
void PsychAllocateNativeUnsignedByteMat(psych_int64 m, psych_int64 n, psych_int64 p, psych_uint8 **cArray, PsychGenericScriptType **nativeElement)
{
    psych_uint8 *cArrayTemp;

    PsychCheckSizeLimits(m, n, p);
    *nativeElement = mxCreateByteMatrix3D(m,n,p);
    cArrayTemp = mxGetData(*nativeElement);
    *cArray = cArrayTemp;
}


double PsychGetNanValue(void)
{
    return(mxGetNaN());
}


/* PsychRuntimeGetPsychtoolboxRoot() TODO FIXME
 *
 * Try to retrieve filesystem path to Psychtoolbox root folder (the result from PsychtoolboxRoot() in Matlab/Octave)
 * or users configuration folder (the result from PsychtoolboxConfigDir() in Matlab/Octave) from runtime. The result,
 * if any, will be cached for later fast lookup.
 *
 * getConfigDir = TRUE => Return PsychtoolboxConfigDir().
 * getConfigDir = FALSE => Return PsychtoolboxRoot().
 *
 * This function may fail to retrieve the path, in which case it returns an empty null-terminated string, i.e., strlen() == 0.
 * On successful recovery of the path, returns a const char* to a readonly string which encodes the path.
 *
 */
const char* PsychRuntimeGetPsychtoolboxRoot(psych_bool getConfigDir)
{
    static char         psychtoolboxRootPath[FILENAME_MAX+1];
    static char         psychtoolboxConfigPath[FILENAME_MAX+1];
/*
    static psych_bool   firstTime = TRUE;
    char*               myPathvarChar = NULL;
    PyObject            *plhs[1]; // Capture the runtime result of PsychtoolboxRoot/ConfigDir

    if (firstTime) {
        // Reset firstTime flag:
        firstTime = FALSE;

        // Init to null-terminated empty strings, so it is well-defined in case of error:
        psychtoolboxRootPath[0] = 0;
        psychtoolboxConfigPath[0] = 0;

        // Call into runtime to get the path to the root folder: This will return 0 on success.
        // A non-zero return value probably means that the script wasn't in the path.
        if (0 == Psych_mexCallMATLAB(1, plhs, 0, NULL, "PsychtoolboxRoot")) {
            myPathvarChar = mxArrayToString(plhs[0]);
            if (myPathvarChar) {
                strncpy(psychtoolboxRootPath, myPathvarChar, FILENAME_MAX);
                mxFree(myPathvarChar);
            }
        }
        mxDestroyArray(plhs[0]);

        // At this point we did our best and psychtoolboxRootPath is valid: Either a path string,
        // or an empty string signalling failure to get the path.

        // Same game again for PsychtoolboxConfigDir:
        if (0 == Psych_mexCallMATLAB(1, plhs, 0, NULL, "PsychtoolboxConfigDir")) {
            myPathvarChar = mxArrayToString(plhs[0]);
            if (myPathvarChar) {
                strncpy(psychtoolboxConfigPath, myPathvarChar, FILENAME_MAX);
                mxFree(myPathvarChar);
            }
        }
        mxDestroyArray(plhs[0]);
    }
*/
    // Return whatever we've got:
    return((getConfigDir) ? &psychtoolboxConfigPath[0] : &psychtoolboxRootPath[0]);
}


/* PsychCopyInPointerArg() - Copy in a void* memory pointer. */
psych_bool PsychCopyInPointerArg(int position, PsychArgRequirementType isRequired, void **ptr)
{
    PyObject        *ppyPtr;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_unclassified, isRequired, 1,1,1,1,1,1);
    matchError = PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        ppyPtr = (PyObject*) PsychGetInArgPyPtr(position);
        *ptr = (void*) PyLong_AsVoidPtr(ppyPtr);
    }

    return(acceptArg);
}


/* PsychCopyOutPointerArg() - Copy out a void* memory pointer. */
psych_bool PsychCopyOutPointerArg(int position, PsychArgRequirementType isRequired, void* ptr)
{
    PyObject        **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_unclassified, isRequired, 1,1,1,1,0,0);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgPyPtr(position);
        *mxpp = PyLong_FromVoidPtr(ptr);
    }

    return(putOut);
}


/* PsychRuntimePutVariable()
 *
 * Copy a given native variable of type PsychGenericScriptType, e.g., as created by PsychAllocateNativeDoubleMat()
 * in case of a double matrix, as a new variable into a specified workspace.
 *
 * 'workspace'   Namestring of workspace: "base" assign to base workspace. "caller" assign to calling functions workspace,
 *               "global" assign to global variable with given name.
 *
 *               CAUTION: * Python only handles "global" and "caller" atm., "base" is treated like "caller".
 *
 *                        * The calling Python *code* can only directly access "global" variables created via PsychRuntimePutVariable().
 *                          It can access "caller" local variables only via the dict accessor locals()['variable'], and the locals() may
 *                          be behaviour not to be relied on, according to some docs, as it is an implementation detail subject to change.
 *
 *                          The reason for this weirdness is that Python code is compiled into bytecode at function load time, and the
 *                          set of available local variable names is fixed at compile time, so dynamically injected "caller" variables from
 *                          us would be invisible to the compiled bytecode.
 *
 *                        * Python code executed from *within* our module via PsychRuntimeEvaluateString() can access both
 *                          "caller" local, and "global" global variables though, due to the way we implement PsychRuntimeEvaluateString().
 *
 *                        * In fact, the safest choice may be to use "caller" for calling Python code from *within* our module.
 *
 * 'variable'    Name of the new variable.
 *
 * 'pcontent'    The actual content that should be copied into the variable.
 *
 *
 * Example: You want to create a double matrix with (m,n,p) rows/cols/layers as a variable 'myvar' in the global
 *          workspace and initialize it with content from the double array mycontent:
 *
 *          PsychGenericScriptType* newvar = NULL;
 *          double* newvarcontent = mycontent; // mycontent is double* onto existing data.
 *          PsychAllocateNativeDoubleMat(m, n, p, &newvarcontent, &newvar);
 *          At this point, newvar contains the content of 'mycontent' and 'newvarcontent' points to
 *          the copy. You could alter mycontent now without affecting the content of newvarcontent or newvar.
 *
 *          Create the corresponding variable 'myvar' in the global workspace:
 *          PsychRuntimePutVariable("global", "myvar", newvar);
 *
 *          The calling M-File etc. can access the content newvarcontent under the variable name 'myvar'.
 *
 *          As usual, the double matrix newvarcontent will be auto-destroyed when returning to the runtime,
 *          but the variable 'myvar' will remain valid until it goes out of scope.
 *
 *          This function steals a reference to pcontent, regardless of success or failure.
 *
 * Returns zero on success, non-zero on failure.
 */
int PsychRuntimePutVariable(const char* workspace, const char* variable, PsychGenericScriptType* pcontent)
{
    int rc = 1;

    // Get a borrowed reference to the dicts with global and local variables for the calling frame:
    PyObject *dict = strcmp(workspace, "global") ? PyEval_GetLocals() : PyEval_GetGlobals();
    if (dict && PyDict_Check(dict)) {
        // Try to add our new or updated variable, which will increment its refcount:
        rc = PyDict_SetItemString(dict, variable, pcontent);
    }

    if (psych_refcount_debug && pcontent)
        printf("PTB-DEBUG: In mxSetField: refcount of external object %p at enter is %li. %s\n",
               pcontent, Py_REFCNT(pcontent),
               (Py_REFCNT(pcontent) > 1) ? "MIGHT leak if caller does not take care." : "");

    // Drop one reference, so the function steals a reference to pcontent and the
    // calling client is no longer responsible for managing it, regardless of success
    // or failure:
    Py_XDECREF(pcontent);

    return(rc);
}


/* PsychRuntimeGetVariablePtr()
 *
 * Retrieve a *read-only* pointer to a given native variable of type PsychGenericScriptType in the specified workspace.
 * The variable is not copied, just referenced, so you *must not modify/write to the location* only perform read access!
 *
 * 'workspace'   Namestring of workspace: "base" get from base workspace. "caller" get from calling functions workspace,
 *               "global" get global variable with given name.
 *
 *               CAUTION: Python only handles "global" and "caller" atm., "base" is treated like "caller".
 *                        In fact, the safest choice may be to use "caller".
 *
 * 'variable'    Name of the variable to get a reference.
 *
 * 'pcontent'    Pointer to a PsychGenericScriptType* where the location of the variables content should be stored.
 *               The pointed-to pointer will be set to NULL on failure.
 *
 * Returns TRUE on success, FALSE on failure.
 */
psych_bool PsychRuntimeGetVariablePtr(const char* workspace, const char* variable, PsychGenericScriptType** pcontent)
{
    // Init to empty default:
    *pcontent = NULL;

    // Get a borrowed reference to the dicts with global and local variables for the calling frame:
    PyObject *dict = strcmp(workspace, "global") ? PyEval_GetLocals() : PyEval_GetGlobals();
    if (dict && PyDict_Check(dict)) {
        // Get borrowed reference to variable:
        *pcontent = (PsychGenericScriptType*) PyDict_GetItemString(dict, variable);
    }

    // Return true on success, false on failure:
    return((*pcontent) ? TRUE : FALSE);
}

/* PsychRuntimeEvaluateString()
 *
 * Simple function evaluation by the Python scripting environment.
 * This asks the runtime environment to execute/evaluate the given string 'cmdstring',
 * passing no return arguments back, except an error code.
 */
int PsychRuntimeEvaluateString(const char* cmdstring)
{
    PyObject* code = Py_CompileString(cmdstring, "PTB", Py_file_input);
    if (code) {
        PyObject* res = PyEval_EvalCode(code, PyEval_GetGlobals(), PyEval_GetLocals());
        Py_DECREF(code);
        if (res) {
            Py_DECREF(res);
            return(0);
        }
    }
    // Failed:
    return(-1);
}

// functions for outputting structs
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
    PsychAllocOutStructArray()

    -If argument is optional we allocate the structure even if the argument is not present.  If this bothers you,
    then check within the subfunction for the presense of a return argument before creating the struct array.  We
    allocate space regardless of whether the argument is present because this is consistant with other "PsychAllocOut*"
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
    PyObject        **mxArrayOut;
    ptbSize         structArrayNumDims = 1;
    ptbSize         structArrayDims[1];
    PsychError      matchError;
    psych_bool      putOut;

    structArrayDims[0] = numElements;

    // Handle special case numElements == -1 which is ~ == 1 for validation:
    numElements = abs(numElements);

    if (position != kPsychNoArgReturn) {
        // Return the result to both the C caller and the scripting environment.
        PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
        PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_structArray, isRequired, 1, 1, numElements, numElements, 0, 0);
        matchError = PsychMatchDescriptors();
        putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);

        *pStruct = mxCreateStructArray(structArrayNumDims, structArrayDims, numFields, fieldNames);

        if (putOut) {
            mxArrayOut = PsychGetOutArgPyPtr(position);
            *mxArrayOut = *pStruct;
        }

        return(putOut);
    } else{
        // Return the result only to the C caller. Ignore "required".
        *pStruct = mxCreateStructArray(structArrayNumDims, structArrayDims, numFields, fieldNames);
        return(TRUE);
    }
}


/*
    PsychAssignOutStructArray()

    Accept a pointer to a struct array and assign the struct array to be the
    designated return variable.
 */
psych_bool PsychAssignOutStructArray(int position,
                                     PsychArgRequirementType isRequired,
                                     PsychGenericScriptType *pStruct)
{
    PyObject        **mxArrayOut;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_structArray, isRequired, 1, 1, 0, kPsychUnboundedArraySize, 0, 0);
    matchError = PsychMatchDescriptors();
    putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);
    // putOut is always true in Python implementation, so no leakage, as interpreter will own
    // the only reference:
    if (putOut) {
        mxArrayOut = PsychGetOutArgPyPtr(position);
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
    int         fieldNumber;
    psych_bool  isStruct;
    PyObject    *mxFieldValue;
    char        errmsg[256];

    isStruct = mxIsStruct(pStruct);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    fieldNumber = mxIsField(pStruct, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    mxFieldValue = mxCreateString(text);

    // mxSetField steals the only reference to mxFieldValue, so we are done with it:
    mxSetField(pStruct, (ptbIndex) index, fieldName, mxFieldValue);
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
    int         fieldNumber;
    psych_bool  isStruct;
    PyObject    *mxFieldValue;
    char        errmsg[256];

    isStruct = mxIsStruct(pStruct);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    fieldNumber = mxIsField(pStruct, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    mxFieldValue = PyFloat_FromDouble(value);

    // mxSetField steals the only reference to mxFieldValue, so we are done with it:
    mxSetField(pStruct, (ptbIndex) index, fieldName, mxFieldValue);
}


/*
    PsychSetStructArrayUnsignedInt64Element()
    Note: The variable "index" is zero-indexed.
*/
void PsychSetStructArrayUnsignedInt64Element(const char *fieldName,
                                             int index,
                                             psych_uint64 value,
                                             PsychGenericScriptType *pStruct)
{
    int         fieldNumber;
    psych_bool  isStruct;
    PyObject    *mxFieldValue;
    char        errmsg[256];

    isStruct = mxIsStruct(pStruct);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    fieldNumber = mxIsField(pStruct, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    mxFieldValue = PyLong_FromUnsignedLongLong(value);

    // mxSetField steals the only reference to mxFieldValue, so we are done with it:
    mxSetField(pStruct, (ptbIndex) index, fieldName, mxFieldValue);
}


/*
    PsychSetStructArrayBooleanElement()

    Note: The variable "index" is zero-indexed.
*/
void PsychSetStructArrayBooleanElement(const char *fieldName,
                                       int index,
                                       psych_bool state,
                                       PsychGenericScriptType *pStruct)
{
    int         fieldNumber;
    psych_bool  isStruct;
    PyObject    *mxFieldValue;
    char        errmsg[256];

    isStruct = mxIsStruct(pStruct);
    if (!isStruct)
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");

    fieldNumber = mxIsField(pStruct, fieldName);
    if (fieldNumber == -1) {
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    mxFieldValue = PyBool_FromLong((long) state);

    // mxSetField steals the only reference to mxFieldValue, so we are done with it:
    mxSetField(pStruct, (ptbIndex) index, fieldName, mxFieldValue);
}


/*
 * PsychSetStructArrayStructElement()
 *
 * Assign a struct(-array) another struct(-array)'s field.
 *
 * Irrespective of success or failure, this steals the reference to pStructInner, so caller
 * does not have to take care of releasing pStructInner and should not touch it anymore after
 * calling this function, unless it aquired an own extra reference beforehand.
 *
 */
void PsychSetStructArrayStructElement(const char *fieldName,
                                      int index,
                                      PsychGenericScriptType *pStructInner,
                                      PsychGenericScriptType *pStructOuter)
{
    int         fieldNumber;
    psych_bool  isStruct;
    char        errmsg[256];

    isStruct = mxIsStruct(pStructOuter);
    if (!isStruct) {
        Py_XDECREF(pStructInner);
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");
    }

    isStruct = mxIsStruct(pStructInner);
    if (!isStruct) {
        Py_XDECREF(pStructInner);
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a struct field to a non-existent structure.");
    }

    fieldNumber = mxIsField(pStructOuter, fieldName);
    if (fieldNumber == -1) {
        Py_XDECREF(pStructInner);
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    // mxSetField steals the reference to pStructInner, so we are done with it:
    mxSetField(pStructOuter, (ptbIndex) index, fieldName, pStructInner);
}


/*
 * PsychSetStructArrayNativeElement()
 *
 * Assign a native element like a matrix (PsychAllocateNativeDoubleMat/PsychAllocateNativeUnsignedByteMat)
 * to a struct(-array) field.
 *
 * Irrespective of success or failure, this steals the reference to pNativeElement, so caller
 * does not have to take care of releasing pNativeElement and should not touch it anymore after
 * calling this function, unless it aquired an own extra reference beforehand.
 *
 */
void PsychSetStructArrayNativeElement(const char *fieldName,
                                      int index,
                                      PsychGenericScriptType *pNativeElement,
                                      PsychGenericScriptType *pStructArray)
{
    int         fieldNumber;
    psych_bool  isStruct;
    char        errmsg[256];

    isStruct = mxIsStruct(pStructArray);
    if (!isStruct) {
        Py_XDECREF(pNativeElement);
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a field within a non-existent structure.");
    }

    fieldNumber = mxIsField(pStructArray, fieldName);
    if (fieldNumber == -1) {
        Py_XDECREF(pNativeElement);
        sprintf(errmsg, "Attempt to set a non-existent structure name field: %s", fieldName);
        PsychErrorExitMsg(PsychError_internal, errmsg);
    }

    // mxSetField steals the reference to pNativeElement, so we are done with it:
    mxSetField(pStructArray, (ptbIndex) index, fieldName, pNativeElement);
}


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
    PyObject    **mxArrayOut;
    PsychError  matchError;
    psych_bool  putOut;

    if (position != kPsychNoArgReturn) {
        // Return the result to both the C caller and the scripting environment:
        PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
        PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_cellArray, isRequired, 1, 1, numElements, numElements, 0, 0);
        matchError = PsychMatchDescriptors();
        putOut = PsychAcceptOutputArgumentDecider(isRequired, matchError);

        *pCell = PyTuple_New((Py_ssize_t) numElements);

        // putOut is always TRUE under Python, so no leakage, as interpreter will own
        // the only reference:
        if (putOut) {
            mxArrayOut = PsychGetOutArgPyPtr(position);
            *mxArrayOut = *pCell;
        }

        return(putOut);
    } else {
        // Return the result only to the C caller, not to the calling environment. Ignore "required":
        // CAUTION: Could leak if caller doesn't assign it to an output-assignment, e.g., return arg
        // or struct(array):
        *pCell = PyTuple_New((Py_ssize_t) numElements);
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
    PyObject *mxFieldValue;

    // Check for bogus arguments
    if (!PyTuple_Check(cellVector))
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a cell within a non-existent cell array.");

    if ((ptbIndex) index >= (ptbIndex) PyTuple_Size(cellVector))
        PsychErrorExitMsg(PsychError_internal, "Attempt to set a cell array field at an out-of-bounds index");

    // mxFieldValue refcount is 1, PyTuple_SetItem below steals the reference, so
    // now the cellVector (=interpreter on eventual successfull out-assignment) owns the only
    // reference and we don't leak.
    mxFieldValue = mxCreateString(text);
    PyTuple_SetItem(cellVector, index, mxFieldValue);
}

// End of Python only stuff.
#endif
