/*
 * PsychSourceGL/Source/Common/Base/PsychScriptingGlueMatlab.c
 *
 * AUTHORS:
 *
 * Allen.Ingling@nyu.edu        awi
 * mario.kleiner.de@gmail.com   mk
 *
 * PLATFORMS: All -- Glue layer for GNU/Octave and Matlab runtime environment.
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
 * between the calling Matlab- or GNU/Octave environment and the PsychToolbox.
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
 *
 * TO DO:
 *
 * - baseFunctionInvoked[recLevel] and PsychSubfunctionEnabled are redundent, keep only baseFunctionInvoked[recLevel]
 *
 */

// During inclusion of Psych.h, we define the special flag PTBINSCRIPTINGGLUE. This
// will cause some of the system headers in Psych.h not to be included during build
// of PsychScriptingGlue.c:
#define PTBINSCRIPTINGGLUE 1
#include "Psych.h"
#undef PTBINSCRIPTINGGLUE

////Static functions local to PsychScriptingGlueMatlab.c.
// _____________________________________________________________________________________
// for Matlab/Octave
#if PSYCH_LANGUAGE == PSYCH_MATLAB

void InitializeSynopsis(char *synopsis[],int maxStrings);

#define MAX_SYNOPSIS 100
#define MAX_CMD_NAME_LENGTH 100

//Static variables local to ScriptingGlue.c.  The convention is to append a abbreviation in all
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
static int recLevel = -1;
static psych_bool psych_recursion_debug = FALSE;

static psych_bool nameFirstGLUE[MAX_RECURSIONLEVEL];
static psych_bool baseFunctionInvoked[MAX_RECURSIONLEVEL];

static int nlhsGLUE[MAX_RECURSIONLEVEL];  // Number of requested return arguments.
static int nrhsGLUE[MAX_RECURSIONLEVEL];  // Number of provided call arguments.
static mxArray **plhsGLUE[MAX_RECURSIONLEVEL];       // A pointer to the plhs array passed to the MexFunction entry point
static CONSTmxArray **prhsGLUE[MAX_RECURSIONLEVEL]; // A pointer to the prhs array passed to the MexFunction entry point

static void PsychExitGlue(void);

//local function declarations
static psych_bool PsychIsEmptyMat(CONSTmxArray *mat);
static psych_bool PsychIsDefaultMat(CONSTmxArray *mat);
static mwSize mxGetP(const mxArray *array_ptr);
static mwSize mxGetNOnly(const mxArray *arrayPtr);
static mxArray *mxCreateDoubleMatrix3D(psych_int64 m, psych_int64 n, psych_int64 p);

//declarations for functions exported from code module
EXP void mexFunction(int nlhs, mxArray *plhs[], int nrhs, CONSTmxArray *prhs[]);

// firstTime: This flag defines if this is the first invocation of the module
// since it was (re-)loaded:
static psych_bool firstTime = TRUE;

void ScreenCloseAllWindows(void);

// No-Op implementation of Mex specific error handling:
void PsychProcessErrorInScripting(PsychError error, const char* message)
{
    (void) error;
    (void) message;

    return;
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
    // Mex API based environemnts, ie., Octave and Matlab, use Fortran memory layout for
    // their n-D matrices and arrays, ergo classic style is most optimal. Therefore we
    // never opt-in to C-style, as that would be counterproductive:
    (void) tryEnableCMemoryLayout;

    return(FALSE);
}


/*
 *
 *    Main entry point for Matlab and Octave. Serves as a dispatch and handles
 *    first time initialization.
 *
 *    EXP is a macro defined within Psychtoolbox source to be nothing
 *    except on win where it is the declaration which tells the linker to
 *    make the function visible from outside the DLL.
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
EXP void mexFunction(int nlhs, mxArray *plhs[], int nrhs, CONSTmxArray *prhs[])
{
    psych_bool          isArgThere[2], isArgEmptyMat[2], isArgText[2], isArgFunction[2];
    PsychFunctionPtr    fArg[2], baseFunction;
    char                argString[2][MAX_CMD_NAME_LENGTH];
    int                 i;
    const mxArray*      tmparg = NULL; // mxArray is mxArray under MATLAB but #defined to octave_value on OCTAVE build.

    // Initialization
    if (firstTime) {
        // Reset call recursion level to startup default:
        recLevel = -1;
        psych_recursion_debug = FALSE;

        if (getenv("PSYCH_RECURSION_DEBUG")) psych_recursion_debug = TRUE;

        //call the Psychtoolbox init function, which inits the Psychtoolbox and calls the project init.
        PsychInit();

        //register the exit function, which calls PsychProjectExit() to clean up for the project then
        //calls whatever to clean up for all of Psych.h layer.

        // Under Matlab we use the mexAtExit() Mex-API function to register our PsychExitGlue() routine.
        // Whenever Matlab wants to flush our module (to reload it, or in response to Matlab-Shutdown,
        // 'clear MODULENAME', 'clear mex' or 'clear all' command) it first calls our PsychExitGlue(),
        // then unloads the module from memory...
        mexAtExit(&PsychExitGlue);

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
        PsychErrorExitMsg(PsychError_internal, "Module call recursion limit exceeded");
    }

    if (psych_recursion_debug) printf("PTB-DEBUG: Module %s entering recursive call level %i.\n", PsychGetModuleName(), recLevel);

    // Store away call arguments for use by language-neutral accessor functions in ScriptingGlue.c
    nlhsGLUE[recLevel] = nlhs;
    nrhsGLUE[recLevel] = nrhs;
    plhsGLUE[recLevel] = plhs;
    prhsGLUE[recLevel] = prhs;

    baseFunctionInvoked[recLevel]=FALSE;

    //if no subfunctions have been registered by the project then just invoke the project base function
    //if one of those has been registered.
    if (!PsychAreSubfunctionsEnabled()) {
        baseFunction = PsychGetProjectFunction(NULL);
        if (baseFunction != NULL) {
            baseFunctionInvoked[recLevel]=TRUE;
            (*baseFunction)();  //invoke the unnamed function
        } else
            PrintfExit("Project base function invoked but no base function registered");
    } else { //subfunctions are enabled so pull out the function name string and invoke it.
        //assess the nature of first and second arguments for finding the name of the sub function.
        for (i = 0; i < 2; i++) {
            isArgThere[i] = (nrhs>i) && (prhsGLUE[recLevel][i]);
            if (isArgThere[i]) tmparg = prhs[i]; else tmparg = NULL;

            isArgEmptyMat[i] = isArgThere[i] ? mxGetM(tmparg)==0 || mxGetN(tmparg)==0 : FALSE;
            isArgText[i] = isArgThere[i] ? mxIsChar(tmparg) : FALSE;
            if (isArgText[i]) {
                mxGetString(tmparg,argString[i],sizeof(argString[i]));
                // Only consider 2nd arg as subfunction if 1st arg isn't already a subfunction:
                if ((i == 0) || (!isArgFunction[0])) {
                    fArg[i]=PsychGetProjectFunction(argString[i]);
                }
                else fArg[i] = NULL; // 1st arg is subfunction, so 2nd arg can't be as well.
            }
            isArgFunction[i] = isArgText[i] ? fArg[i] != NULL : FALSE;
        }

        //figure out which of the two arguments might be the function name and either invoke it or exit with error
        //if we can't find one.
        if (!isArgThere[0] && !isArgThere[1]) { //no arguments passed so execute the base function
            baseFunction = PsychGetProjectFunction(NULL);
            if (baseFunction != NULL) {
                baseFunctionInvoked[recLevel]=TRUE;
                (*baseFunction)();
            } else
                PrintfExit("Project base function invoked but no base function registered");
        }
        // (!isArgThere[0] && isArgEmptyMat[1]) --disallowed
        // (!isArgThere[0] && isArgText[1])     --disallowed
        // (!isArgThere[0] && !isArgText[1]     --disallowed except in case of !isArgThere[0] caught above.

        else if (isArgEmptyMat[0] && !isArgThere[1])
            PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state A)");
        else if (isArgEmptyMat[0] && isArgEmptyMat[1])
            PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state B)");
        else if (isArgEmptyMat[0] && isArgText[1]) {
            if (isArgFunction[1]) {
                nameFirstGLUE[recLevel] = FALSE;
                (*(fArg[1]))();
            }
            else
                PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state C)");
        }
        else if (isArgEmptyMat[0] && !isArgText[1])
            PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state D)");
        else if (isArgText[0] && !isArgThere[1]) {
            if (isArgFunction[0]) {
                nameFirstGLUE[recLevel] = TRUE;
                (*(fArg[0]))();
            } else { //when we receive a first argument  wich is a string and it is  not recognized as a function name then call the default function
                baseFunction = PsychGetProjectFunction(NULL);
                if (baseFunction != NULL) {
                    baseFunctionInvoked[recLevel]=TRUE;
                    (*baseFunction)();
                } else
                    PrintfExit("Project base function invoked but no base function registered");
            }
        }
        else if (isArgText[0] && isArgEmptyMat[1]) {
            if (isArgFunction[0]) {
                nameFirstGLUE[recLevel] = TRUE;
                (*(fArg[0]))();
            }
            else
                PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state F)");
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
                    PrintfExit("Project base function invoked but no base function registered");
            }
            else if (isArgFunction[0] && isArgFunction[1]) //both arguments are function names
                PrintfExit("Passed two function names");
        }
        else if (isArgText[0] && !isArgText[1]) {
            if (isArgFunction[0]) {
                nameFirstGLUE[recLevel] = TRUE;
                (*(fArg[0]))();
            }
            else
                PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state H)");
        }
        else if (!isArgText[0] && !isArgThere[1]) {  //this was modified for MODULEVersion with WaitSecs.
            //PrintfExit("Invalid command (error state H)");
            baseFunction = PsychGetProjectFunction(NULL);
            if (baseFunction != NULL) {
                baseFunctionInvoked[recLevel]=TRUE;
                (*baseFunction)();  //invoke the unnamed function
            } else
                PrintfExit("Project base function invoked but no base function registered");
        }
        else if (!isArgText[0] && isArgEmptyMat[1])
            PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state I)");
        else if (!isArgText[0] && isArgText[1])
        {
            if (isArgFunction[1]) {
                nameFirstGLUE[recLevel] = FALSE;
                (*(fArg[1]))();
            } else
                PrintfExit("Unknown or invalid subfunction name - Typo? Check spelling of the function name.  (error state J)");
        }
        else if (!isArgText[0] && !isArgText[1]) {  //this was modified for Priority.
            //PrintfExit("Invalid command (error state K)");
            baseFunction = PsychGetProjectFunction(NULL);
            if (baseFunction != NULL) {
                baseFunctionInvoked[recLevel]=TRUE;
                (*baseFunction)();  //invoke the unnamed function
            } else
                PrintfExit("Project base function invoked but no base function registered");
        }
    } //close else

    PsychExitRecursion();
}


/*
 *    Just call the abstracted PsychExit function.  This might seem dumb, but its necessary to
 *    isolate the scripting language dependent stuff from the rest of the toolbox.
 *
 */
void PsychExitGlue(void)
{
    // Perform platform independent shutdown:
    PsychErrorExitMsg(PsychExit(),NULL);

    // And we are dead. Now the runtime will flush us from process memory,
    // at least on Matlab and Octave 3.7+. In any case no further invocation will happen
    // until reload.
}


/*
 *    Return the mxArray pointer to the specified position.  Note that we have some special rules for
 *    for numbering the positions:
 *
 *    0 - This is always the command string or NULL if the project does not register a
 *        dispatch function and does accept subcommands.  If the function does accept sub
 *        commands, in Matlab those may be passed in either the first or second position, but
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
 *
 *    2 - this function should be used by the one which gets the function name.
 *
 *
 *    Arguments are numbered 0..n.
 *
 *        -The 0th argument is a pointer to the mxArray holding
 *        the subfunction name string if we are in subfunction mode.
 *
 *        -The 0th argument is undefined if not in subfunction mode.
 *
 *        -The 1st argument is the argument of the 1st and 2nd which is not
 *        the subfunction name if in subfunction mode.
 *
 *        -The 1st argument is the first argument if not in subfunction mode.
 *
 *        -The 2nd-nth arguments are always the 2nd-nth arguments.
 */
//we return NULL if a postion without an arg is specified.
const mxArray *PsychGetInArgMxPtr(int position)
{
    if (PsychAreSubfunctionsEnabled() && !baseFunctionInvoked[recLevel]) { //when in subfunction mode
        if (position < nrhsGLUE[recLevel]) { //an argument was passed in the correct position.
            if (position == 0) { //caller wants the function name argument.
                if (nameFirstGLUE[recLevel])
                    return(prhsGLUE[recLevel][0]);
                else
                    return(prhsGLUE[recLevel][1]);
            } else if (position == 1) { //they want the "first" argument.
                if (nameFirstGLUE[recLevel])
                    return(prhsGLUE[recLevel][1]);
                else
                    return(prhsGLUE[recLevel][0]);
            } else
                return(prhsGLUE[recLevel][position]);
        } else
            return(NULL);
    } else { //when not in subfunction mode and the base function is not invoked.
        if (position <= nrhsGLUE[recLevel])
            return(prhsGLUE[recLevel][position-1]);
        else
            return(NULL);
    }
}


mxArray **PsychGetOutArgMxPtr(int position)
{
    if (position==1 || (position>0 && position<=nlhsGLUE[recLevel])) { //an ouput argument was supplied at the specified location
        return(&(plhsGLUE[recLevel][position-1]));
    } else
        return(NULL);
}


const PsychGenericScriptType *PsychGetInArgPtr(int position)
{
    return((const PsychGenericScriptType*) PsychGetInArgMxPtr(position));
}


/* PsychCheckSizeLimits(size_t m, size_t n, size_t p)
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
void PsychCheckSizeLimits(psych_int64 m, psych_int64 n, psych_int64 p)
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
 *    Get the third array dimension which we call "P".  mxGetP should act just like mxGetM and mxGetN.
 *
 *    The abstracted Psychtoolbox API supports matrices with up to 3 dimensions.
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
 *    Get the 2nd array dimension.
 *
 *    The Mex API's mxGetN is sometimes undersirable because it returns the product of all dimensions above 1.  Our mxGetNOnly only returns N, for when you need that.
 *
 *    The abstracted Psychtoolbox API supports matrices with up to 3 dimensions.
 */
static mwSize mxGetNOnly(const mxArray *arrayPtr)
{
    const mwSize *dimArray;

    dimArray = (const mwSize*) mxGetDimensions(arrayPtr);
    return dimArray[1];
}


/*
 *    mxCreateDoubleMatrix3D()
 *
 *    Create a 2D or 3D matrix of doubles.
 *
 *    Requirements are that m>0, n>0, p>=0.
 */
mxArray *mxCreateDoubleMatrix3D(psych_int64 m, psych_int64 n, psych_int64 p)
{
    int numDims;
    mwSize dimArray[3];

    if (m==0 || n==0 ) {
        dimArray[0]=0;dimArray[1]=0;dimArray[2]=0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m,n,p);
        dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
    }
    numDims= (p==0 || p==1) ? 2 : 3;

    return mxCreateNumericArray(numDims, (mwSize*) dimArray, mxDOUBLE_CLASS, mxREAL);
}


/*
 *    mxCreateFloatMatrix3D()
 *
 *    Create a 2D or 3D matrix of floats.
 *
 *    Requirements are that m>0, n>0, p>=0.
 */
mxArray *mxCreateFloatMatrix3D(size_t m, size_t n, size_t p)
{
    int numDims;
    mwSize dimArray[3];

    if (m==0 || n==0 ) {
        dimArray[0]=0;dimArray[1]=0;dimArray[2]=0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m,n,p);
        dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
    }
    numDims= (p==0 || p==1) ? 2 : 3;

    return mxCreateNumericArray(numDims, (mwSize*) dimArray, mxSINGLE_CLASS, mxREAL);
}


/*
 *    mxCreateNativeBooleanMatrix3D()
 *
 *    Create a 2D or 3D matrix of native psych_bool types.
 *
 *    Requirements are that m>0, n>0, p>=0.
 */
mxArray *mxCreateNativeBooleanMatrix3D(size_t m, size_t n, size_t p)
{
    int         numDims;
    mwSize      dimArray[3];
    mxArray     *newArray;

    if (m==0 || n==0 ) {
        dimArray[0]=0;dimArray[1]=0;dimArray[2]=0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m,n,p);
        dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
    }
    numDims = (p==0 || p==1) ? 2 : 3;
    newArray = mxCreateNumericArray(numDims, (mwSize*) dimArray, mxLOGICAL_CLASS, mxREAL);

    return(newArray);
}


/*
 *    Create a 2D or 3D matrix of ubytes.
 *
 *    Requirements are that m>0, n>0, p>=0.
 */
mxArray *mxCreateByteMatrix3D(size_t m, size_t n, size_t p)
{
    int numDims;
    mwSize dimArray[3];

    if (m==0 || n==0 ) {
        dimArray[0]=0;dimArray[1]=0;dimArray[2]=0; //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m,n,p);
        dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
    }
    numDims= (p==0 || p==1) ? 2 : 3;
    return mxCreateNumericArray(numDims, (mwSize*) dimArray, mxUINT8_CLASS, mxREAL);
}


/*
 *    Print string s and return return control to the calling environment.
 */
void PsychErrMsgTxt(char *s)
{
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

    PsychExitRecursion();

    // Call the Matlab- or Octave error printing and error handling facilities:
    mexErrMsgTxt((s && (strlen(s) > 0)) ? s : "See error message printed above.");
}


/*
 *    classify the mxArray element format using Pyschtoolbox argument type names
 *
 */
static PsychArgFormatType PsychGetTypeFromMxPtr(const mxArray *mxPtr)
{
    PsychArgFormatType format;

    if (PsychIsDefaultMat(mxPtr))
        format = PsychArgType_default;
    else if (mxIsUint8(mxPtr))
        format = PsychArgType_uint8;
    else if (mxIsUint16(mxPtr))
        format = PsychArgType_uint16;
    else if (mxIsUint32(mxPtr))
        format = PsychArgType_uint32;
    else if (mxIsInt8(mxPtr))
        format = PsychArgType_int8;
    else if (mxIsInt16(mxPtr))
        format = PsychArgType_int16;
    else if (mxIsInt32(mxPtr))
        format = PsychArgType_int32;
    else if (mxIsDouble(mxPtr))
        format = PsychArgType_double;
    else if (mxIsSingle(mxPtr))
        format = PsychArgType_single;
    else if (mxIsChar(mxPtr))
        format = PsychArgType_char;
    else if (mxIsCell(mxPtr))
        format = PsychArgType_cellArray;
    else if (mxIsLogical(mxPtr))
        // This is tricky because MATLAB abstracts "logicals" conditionally on platform.
        // Depending on OS, MATLAB implements booleans with either 8-bit or 64-bit values.
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
    const mxArray           *mxPtr;

    d.position = argNum;
    d.direction = direction;
    if (direction == PsychArgIn) {
        mxPtr = PsychGetInArgMxPtr(argNum);
        d.isThere = (mxPtr && !PsychIsDefaultMat(mxPtr)) ? kPsychArgPresent : kPsychArgAbsent;
        if (d.isThere == kPsychArgPresent) { //the argument is there so fill in the rest of the description
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


//local function definitions for ScriptingGlue.c
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


int PsychGetNumOutputArgs(void)
{
    return(nlhsGLUE[recLevel]==0 ? 1 : nlhsGLUE[recLevel]);
}


int PsychGetNumNamedOutputArgs(void)
{
    return(nlhsGLUE[recLevel]);
}


PsychError PsychCapNumOutputArgs(int maxNamedOutputs)
{
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
            return(!(PsychIsDefaultMat(PsychGetInArgMxPtr(position)))); //check if its default
        else
            return(FALSE);
    }
}


PsychArgFormatType PsychGetArgType(int position) //this is for inputs because outputs are unspecified by the calling environment.
{
    if (!(PsychIsArgReallyPresent(PsychArgIn, position)))
        return(PsychArgType_none);

    return(PsychGetTypeFromMxPtr(PsychGetInArgMxPtr(position)));
}


size_t PsychGetArgM(int position)
{
    if (!(PsychIsArgPresent(PsychArgIn, position)))
        PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
    return( mxGetM(PsychGetInArgMxPtr(position)));
}


size_t PsychGetArgN(int position)
{
    if (!(PsychIsArgPresent(PsychArgIn, position)))
        PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
    return( mxGetNOnly(PsychGetInArgMxPtr(position)));
}


size_t PsychGetArgP(int position)
{
    if (!(PsychIsArgPresent(PsychArgIn, position)))
        PsychErrorExitMsg(PsychError_invalidArgRef,NULL);
    return( mxGetP(PsychGetInArgMxPtr(position)));
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
 * returns control to the Matlab environemnt.  To make it non volatile, call Psych??? on it.
 *
 * ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 */

psych_bool PsychCopyOutDoubleArg(int position, PsychArgRequirementType isRequired, double value)
{
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double,  isRequired, 1,1,1,1,0,0);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        (*mxpp) = mxCreateDoubleMatrix(1,1,mxREAL);
        mxGetPr((*mxpp))[0] = value;
    }
    return(putOut);
}

psych_bool PsychAllocOutDoubleArg(int position, PsychArgRequirementType isRequired, double **value)
{
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, 1,1,1,1,0,0);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        *mxpp = mxCreateDoubleMatrix3D(1,1,0);
        *value = mxGetPr(*mxpp);
    } else {
        mxpp = PsychGetOutArgMxPtr(position);
        *value= (double *)mxMalloc(sizeof(double));
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
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, m,m,n,n,p,p);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        *mxpp = mxCreateDoubleMatrix3D(m,n,p);
        *array = mxGetPr(*mxpp);
    }else
        *array= (double *) mxMalloc(sizeof(double) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    return(putOut);
}


/*
 * PsychAllocOutFloatMatArg()
 *
 * This function allocates out a matrix of single precision floating point type,
 * that is C data type 32-bit float or Matlab/Octave data type single().
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
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_single, isRequired, m,m,n,n,p,p);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        *mxpp = mxCreateFloatMatrix3D((size_t) m, (size_t) n, (size_t) p);
        *array = (float*) mxGetData(*mxpp);
    }else
        *array = (float*) mxMalloc(sizeof(float) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    return(putOut);
}


/*
 *    PsychCopyOutBooleanArg()
 */
psych_bool PsychCopyOutBooleanArg(int position, PsychArgRequirementType isRequired, PsychNativeBooleanType value)
{
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_boolean, isRequired, 1,1,1,1,0,0);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        (*mxpp) = mxCreateLogicalMatrix(1,1);
        mxGetLogicals((*mxpp))[0] = value;
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
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_boolean, isRequired, m,m,n,n,p,p);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        *mxpp = mxCreateNativeBooleanMatrix3D((size_t) m, (size_t) n, (size_t) p);
        *array = (PsychNativeBooleanType *)mxGetLogicals(*mxpp);
    } else {
        *array= (PsychNativeBooleanType *) mxMalloc(sizeof(PsychNativeBooleanType) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
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
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_uint8, isRequired, m,m,n,n,p,p);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        *mxpp = mxCreateByteMatrix3D((size_t) m, (size_t) n, (size_t) p);
        *array = (psych_uint8 *)mxGetData(*mxpp);
    } else {
        *array= (psych_uint8 *) mxMalloc(sizeof(psych_uint8) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
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
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;
    mwSize          dimArray[3];
    int             numDims;

    // Compute output array dimensions:
    if (m<=0 || n<=0) {
        dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;  //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m,n,p);
        dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
    }
    numDims = (p == 0 || p == 1) ? 2 : 3;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_uint16, isRequired, m,m,n,n,p,p);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        *mxpp = mxCreateNumericArray(numDims, (mwSize*) dimArray, mxUINT16_CLASS, mxREAL);
        *array = (psych_uint16 *)mxGetData(*mxpp);
    } else {
        *array= (psych_uint16 *) mxMalloc(sizeof(psych_uint16) * (size_t) m * (size_t) n * (size_t) maxInt(1,p));
    }
    return(putOut);
}


psych_bool PsychCopyOutDoubleMatArg(int position, PsychArgRequirementType isRequired, psych_int64 m, psych_int64 n, psych_int64 p, double *fromArray)
{
    mxArray     **mxpp;
    double      *toArray;
    PsychError  matchError;
    psych_bool  putOut;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_double, isRequired, m,m,n,n,p,p);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
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
    mxArray         **mxpp;
    psych_uint16    *toArray;
    PsychError      matchError;
    psych_bool      putOut;
    mwSize          dimArray[3];
    int             numDims;

    // Compute output array dimensions:
    if (m<=0 || n<=0) {
        dimArray[0] = 0; dimArray[1] = 0; dimArray[2] = 0;    //this prevents a 0x1 or 1x0 empty matrix, we want 0x0 for empty matrices.
    } else {
        PsychCheckSizeLimits(m,n,p);
        dimArray[0] = (mwSize) m; dimArray[1] = (mwSize) n; dimArray[2] = (mwSize) p;
    }
    numDims = (p == 0 || p == 1) ? 2 : 3;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_uint16, isRequired, m,m,n,n,p,p);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
        *mxpp = mxCreateNumericArray(numDims, (mwSize*) dimArray, mxUINT16_CLASS, mxREAL);
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
    mxArray     **mxpp;
    PsychError  matchError;
    psych_bool  putOut;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, PsychArgType_char, isRequired, 0, strlen(str),0,strlen(str),0,0);
    matchError=PsychMatchDescriptors();
    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
        mxpp = PsychGetOutArgMxPtr(position);
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
// TO DO: Needs to be updated for kPsychArgAnything
psych_bool PsychAllocInDoubleMatArg(int position, PsychArgRequirementType isRequired, int *m, int *n, int *p, double **array)
{
    const mxArray   *mxPtr;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,-1,1,-1,0,-1);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
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
    const mxArray     *mxPtr;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, TRUE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,-1,1,-1,0,-1);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
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
 *
 * Alloc-In a single precision floating point matrix, i.e. a matrix of
 * C data type 32 bit float, aka Matlab/Octave data type single().
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
    const mxArray   *mxPtr;
    PsychError      matchError;
    psych_bool      acceptArg;
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
        if (acceptArg) {
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
    if (acceptArg) {
        mxPtr = PsychGetInArgMxPtr(position);
        *m = (psych_int64) mxGetM(mxPtr);
        *n = (psych_int64) mxGetNOnly(mxPtr);
        *p = (psych_int64) mxGetP(mxPtr);
        *array = (float*) mxGetData(mxPtr);
    }
    return(acceptArg);
}


/*
 *    PsychAllocInIntegerListArg()
 *
 *    In a scriptiong language such as MATLAB where numbers are almost always stored as doubles, this function is useful to check
 *    that the value input is an integer value stored within a double type.
 *
 *    Otherwise it just here to imitate the version written for other scripting languages.
 */
psych_bool PsychAllocInIntegerListArg(int position, PsychArgRequirementType isRequired, int *numElements, int **array)
{
    int         m, n, p, i;
    double      *doubleMatrix;
    psych_bool  isThere;

    isThere=PsychAllocInDoubleMatArg(position, isRequired, &m, &n, &p, &doubleMatrix);
    if (!isThere)
        return(FALSE);
    p= (p==0) ? 1 : p;

    if ((psych_uint64) m * (psych_uint64) n * (psych_uint64) p >= INT_MAX) {
        printf("PTB-ERROR: %i th input argument has more than 2^31 - 1 elements! This is not supported.\n", position);
        return(FALSE);
    }

    *numElements = m * n * p;
    *array = (int*) mxMalloc((size_t) *numElements * sizeof(int));
    for (i=0; i < *numElements; i++) {
        if (!PsychIsIntegerInDouble(doubleMatrix+i))
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
    const mxArray     *mxPtr;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_uint8, isRequired, 1,-1,1,-1,0,-1);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        mxPtr = PsychGetInArgMxPtr(position);
        *m = (int) mxGetM(mxPtr);
        *n = (int) mxGetNOnly(mxPtr);
        *p = (int) mxGetP(mxPtr);
        *array=(unsigned char *)mxGetData(mxPtr);
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
// TO DO: Needs to be updated for kPsychArgAnything
psych_bool PsychCopyInDoubleArg(int position, PsychArgRequirementType isRequired, double *value)
{
    const mxArray   *mxPtr;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,1,1,1,1,1);
    matchError=PsychMatchDescriptors();

    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        mxPtr = PsychGetInArgMxPtr(position);
        *value=mxGetPr(mxPtr)[0];
    }

    return(acceptArg);
}


/*
 *    Like PsychCopyInDoubleArg() with the additional restriction that the passed value not have a fractoinal componenet
 *    and that the it fit within thebounds of a C integer
 *
 *    We could also accept matlab native integer types by specifying a conjunction of those as the third argument
 *    in the PsychSetSpecifiedArgDescriptor() call, but why bother ?
 */
psych_bool PsychCopyInIntegerArg(int position,  PsychArgRequirementType isRequired, int *value)
{
    const mxArray     *mxPtr;
    double            tempDouble;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,1,1,1,1,1);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        mxPtr = PsychGetInArgMxPtr(position);
        tempDouble=mxGetPr(mxPtr)[0];
        if (!PsychIsIntegerInDouble(&tempDouble))
            PsychErrorExit(PsychError_invalidIntegerArg);
        *value=(int)tempDouble;
    }
    return(acceptArg);
}


psych_bool PsychCopyInIntegerArg64(int position,  PsychArgRequirementType isRequired, psych_int64 *value)
{
    const mxArray     *mxPtr;
    double            tempDouble;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,1,1,1,1,1);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        mxPtr = PsychGetInArgMxPtr(position);
        tempDouble=mxGetPr(mxPtr)[0];
        if (!PsychIsInteger64InDouble(&tempDouble))
            PsychErrorExit(PsychError_invalidIntegerArg);
        *value= (psych_int64) tempDouble;
    }
    return(acceptArg);
}


/*
 *    PsychAllocInDoubleArg()
 *
 */
psych_bool PsychAllocInDoubleArg(int position, PsychArgRequirementType isRequired, double **value)
{
    const mxArray     *mxPtr;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_double, isRequired, 1,1,1,1,1,1);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        mxPtr = PsychGetInArgMxPtr(position);
        *value=mxGetPr(mxPtr);
    }
    return(acceptArg);
}


/*
 *    PsychAllocInCharArg()
 *
 *    Reads in a string and sets *str to point to the string.
 *
 *    This function violates the rule for AllocIn fuctions that if the argument is optional and absent we allocate
 *    space.  That turns out to be an unuseful feature anyway, so we should probably get ride of it.
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
    const mxArray   *mxPtr;
    int             status;
    psych_uint64    strLen;
    PsychError      matchError;
    psych_bool      acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, PsychArgType_char, isRequired, 0, kPsychUnboundedArraySize ,0, kPsychUnboundedArraySize, 0 , 1);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        mxPtr  = PsychGetInArgMxPtr(position);
        strLen = ((psych_uint64) mxGetM(mxPtr) * (psych_uint64) mxGetNOnly(mxPtr) * (psych_uint64) sizeof(mxChar)) + 1;
        if (strLen >= INT_MAX) PsychErrorExitMsg(PsychError_user, "Tried to pass in a string with more than 2^31 - 1 characters. Unsupported!");
        *str   = (char *) PsychCallocTemp((size_t) strLen, sizeof(char));
        status = mxGetString(mxPtr, *str, (mwSize) strLen);
        if (status!=0)
            PsychErrorExitMsg(PsychError_internal, "mxGetString failed to get the string");
    }
    return(acceptArg);
}


/*
 *    Get a psych_bool flag from the specified argument position.  The matlab type can be be psych_bool, uint8, or
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
    const mxArray     *mxPtr;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType)(PsychArgType_double|PsychArgType_char|PsychArgType_uint8|PsychArgType_boolean),
                                   isRequired, 1,1,1,1,kPsychUnusedArrayDimension,kPsychUnusedArrayDimension);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        //unlike other PsychAllocIn* functions, here we allocate new memory and copy the input to it rather than simply returning a pointer into the received array.
        //That's because we want the booleans returned to the caller by PsychAllocInFlagArg() to alwyas be 8-bit booleans, yet we accept as flags either 64-bit double, char,
        //or logical type.  Restricting to logical type would be a nuisance in the MATLAB environment and does not solve the problem because on some platforms MATLAB
        //uses for logicals 64-bit doubles and on others 8-bit booleans (check your MATLAB mex/mx header files).
        *argVal = (psych_bool *)mxMalloc(sizeof(psych_bool));
        mxPtr = PsychGetInArgMxPtr(position);
        if (mxIsLogical(mxPtr)) {
            if (mxGetLogicals(mxPtr)[0])
                **argVal=(psych_bool)1;
            else
                **argVal=(psych_bool)0;
        } else {
            if (mxGetScalar(mxPtr))
                **argVal=(psych_bool)1;
            else
                **argVal=(psych_bool)0;
        }
    }
    return(acceptArg);    //the argument was not present (and optional).
}


psych_bool PsychAllocInFlagArgVector(int position,  PsychArgRequirementType isRequired, int *numElements, psych_bool **argVal)
{
    const mxArray     *mxPtr;
    PsychError        matchError;
    psych_bool        acceptArg;
    int               i;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    // MK: Disabled. Doesn't work without conversion of mxGetData into many subcases...
    // PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType)(PsychArgType_double | PsychArgType_char | PsychArgType_uint8 | PsychArgType_boolean),
    //               isRequired, 1, kPsychUnboundedArraySize, 1, kPsychUnboundedArraySize, kPsychUnusedArrayDimension, kPsychUnusedArrayDimension);

    // Ok. Let's see if anybody ever complains about this...
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType)(PsychArgType_double | PsychArgType_boolean),
                                   isRequired, 1, kPsychUnboundedArraySize, 1, kPsychUnboundedArraySize, kPsychUnusedArrayDimension, kPsychUnusedArrayDimension);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
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
        for (i=0; i< *numElements;i++) {
            if (mxIsLogical(mxPtr)) {
                if (mxGetLogicals(mxPtr)[i])
                    (*argVal)[i]=(psych_bool)1;
                else
                    (*argVal)[i]=(psych_bool)0;
            } else {
                if (mxGetPr(mxPtr)[i])
                    (*argVal)[i]=(psych_bool)1;
                else
                    (*argVal)[i]=(psych_bool)0;
            }
        }
    }
    return(acceptArg);    //the argument was not present (and optional).
}


/*
 *    PsychCopyInFlagArg()
 */
psych_bool PsychCopyInFlagArg(int position, PsychArgRequirementType isRequired, psych_bool *argVal)
{
    const mxArray     *mxPtr;
    PsychError        matchError;
    psych_bool        acceptArg;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, (PsychArgFormatType)(PsychArgType_double|PsychArgType_char|PsychArgType_uint8|PsychArgType_boolean),
                                   isRequired, 1,1,1,1,kPsychUnusedArrayDimension,kPsychUnusedArrayDimension);
    matchError=PsychMatchDescriptors();
    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
        mxPtr = PsychGetInArgMxPtr(position);
        if (mxIsLogical(mxPtr)) {
            if (mxGetLogicals(mxPtr)[0])
                *argVal=(psych_bool)1;
            else
                *argVal=(psych_bool)0;
        } else {
            if (mxGetScalar(mxPtr))
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


// functions which allocate native types without assigning them to return arguments.
// this is useful for embedding native structures within each other.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 *    PsychAllocateNativeDoubleMat()
 *
 *    Create an opaque native matrix.   Return both
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
    cArrayTemp = mxGetPr(*nativeElement);
    *cArray=cArrayTemp;
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
    static psych_bool   firstTime = TRUE;
    static char         psychtoolboxRootPath[FILENAME_MAX+1];
    static char         psychtoolboxConfigPath[FILENAME_MAX+1];
    char*               myPathvarChar = NULL;
    mxArray             *plhs[1]; // Capture the runtime result of PsychtoolboxRoot/ConfigDir

    if (firstTime) {
        // Reset firstTime flag:
        firstTime = FALSE;

        // Init to null-terminated empty strings, so it is well-defined in case of error:
        psychtoolboxRootPath[0] = 0;
        psychtoolboxConfigPath[0] = 0;

        mexSetTrapFlag(1); // Tell Octave/MATLAB that we'll handle exceptions
        // We could have used mexCallMATLABWithTrap below, but Octave doesn't support it.
        // Unfortunately, MATLAB plans to deprecate mexSetTrapFlag in lieu of *WithTrap, so
        // we'll need to patch Octave or make some sort of wrapper function in here.

        // Call into runtime to get the path to the root folder: This will return 0 on success.
        // A non-zero return value probably means that the script wasn't in the path. When that
        // happens, there will be an error in the command window, but control stays with the mex
        // file (thanks to mexSetTrapFlag(1) above) and it'll continue to run.
        if (0 == mexCallMATLAB(1, plhs, 0, NULL, "PsychtoolboxRoot")) {
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
        if (0 == mexCallMATLAB(1, plhs, 0, NULL, "PsychtoolboxConfigDir")) {
            myPathvarChar = mxArrayToString(plhs[0]);
            if (myPathvarChar) {
                strncpy(psychtoolboxConfigPath, myPathvarChar, FILENAME_MAX);
                mxFree(myPathvarChar);
            }
        }
        mxDestroyArray(plhs[0]);

        mexSetTrapFlag(0); // Tell the runtime that we'll no longer handle exceptions
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
    const mxArray     *mxPtr;
    PsychError        matchError;
    psych_bool        acceptArg;
    psych_bool        is64Bit;

    // 64 bit build?
    is64Bit = sizeof(size_t) > 4;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgIn);
    PsychSetSpecifiedArgDescriptor(position, PsychArgIn, ((is64Bit) ? PsychArgType_uint64 : PsychArgType_uint32), isRequired, 1,1,1,1,1,1);
    matchError=PsychMatchDescriptors();

    acceptArg=PsychAcceptInputArgumentDecider(isRequired, matchError);
    if (acceptArg) {
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
    mxArray         **mxpp;
    PsychError      matchError;
    psych_bool      putOut;
    psych_bool      is64Bit;
    mwSize          dimArray[2];
    int             numDims = 2;

    dimArray[0] = dimArray[1] = 1;

    // 64 bit build?
    is64Bit = sizeof(size_t) > 4;

    PsychSetReceivedArgDescriptor(position, FALSE, PsychArgOut);
    PsychSetSpecifiedArgDescriptor(position, PsychArgOut, ((is64Bit) ? PsychArgType_uint64 : PsychArgType_uint32), isRequired, 1,1,1,1,0,0);
    matchError=PsychMatchDescriptors();

    putOut=PsychAcceptOutputArgumentDecider(isRequired, matchError);
    if (putOut) {
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
 * Copy a given native variable of type PsychGenericScriptType, e.g., as created by PsychAllocateNativeDoubleMat()
 * in case of a double matrix, as a new variable into a specified workspace.
 *
 * 'workspace'  Namestring of workspace: "base" copy to base workspace. "caller" copy into calling functions workspace,
 *              "global" create new global variable with given name.
 *
 *              CAUTION: Some non-Octave/Matlab scripting environments may only handle "global" and "caller" atm.,
 *                       "base" is treated like "caller", so better avoid "base" for maximum portability.
 *                       In fact, the safest choice may be to use "caller".
 *
 * 'variable'   Name of the new variable.
 *
 * 'pcontent'   The actual content that should be copied into the variable.
 *
 *
 * Example: You want to create a double matrix with (m,n,p) rows/cols/layers as a variable 'myvar' in the base
 *          workspace and initialize it with content from the double array mycontent:
 *
 *          PsychGenericScriptType* newvar = NULL;
 *          double* newvarcontent = mycontent; // mycontent is double* onto existing data.
 *          PsychAllocateNativeDoubleMat(m, n, p, &newvarcontent, &newvar);
 *          At this point, newvar contains the content of 'mycontent' and 'newvarcontent' points to
 *          the copy. You could alter mycontent now without affecting the content of newvarcontent or newvar.
 *
 *          Create the corresponding variable in the base workspace:
 *          PsychRuntimePutVariable("base", "myvar", newvar);
 *
 *          The calling M-File etc. can access the content newvarcontent under the variable name 'myvar'.
 *
 *          As usual, the double matrix newvarcontent will be auto-destroyed when returning to the runtime,
 *          but the variable 'myvar' will remain valid until it goes out of scope.
 *
 * Returns zero on success, non-zero on failure.
 */
int PsychRuntimePutVariable(const char* workspace, const char* variable, PsychGenericScriptType* pcontent)
{
    return(mexPutVariable(workspace, variable, pcontent));
}


/* PsychRuntimeGetVariablePtr()
 *
 * Retrieve a *read-only* pointer to a given native variable of type PsychGenericScriptType in the specified workspace.
 * The variable is not copied, just referenced, so you *must not modify/write to the location* only perform read access!
 *
 * 'workspace'  Namestring of workspace: "base" get from base workspace. "caller" get from calling functions workspace,
 *              "global" get global variable with given name.
 *
 *              CAUTION: Some non-Octave/Matlab scripting environments may only handle "global" and "caller" atm.,
 *                       "base" is treated like "caller", so better avoid "base" for maximum portability.
 *                       In fact, the safest choice may be to use "caller".
 *
 * 'variable'   Name of the variable to get a reference.
 *
 * 'pcontent'   Pointer to a PsychGenericScriptType* where the location of the variables content should be stored.
 *              The pointed to pointer will be set to NULL on failure.
 *
 * Returns TRUE on success, FALSE on failure.
 */
psych_bool PsychRuntimeGetVariablePtr(const char* workspace, const char* variable, PsychGenericScriptType** pcontent)
{
    *pcontent = (PsychGenericScriptType*) mexGetVariablePtr(workspace, variable);

    // Return true on success, false on failure:
    return((*pcontent) ? TRUE : FALSE);
}

/* PsychRuntimeEvaluateString()
 *
 * Simple function evaluation by scripting environment via feval() style functions.
 * This asks the runtime environment to execute/evaluate the given string 'cmdstring',
 * passing no return arguments back, except an error code.
 *
 */
int PsychRuntimeEvaluateString(const char* cmdstring)
{
    return(mexEvalString(cmdstring));
}

//end of Matlab & Octave only stuff.
#endif
