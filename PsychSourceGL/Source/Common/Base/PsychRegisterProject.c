/*
  PsychToolbox3/Source/Common/PsychRegisterProject.cpp		
  
  PROJECTS: All.  
  
  AUTHORS:
	Allen.Ingling@nyu.edu              awi 
	mario.kleiner@tuebingen.mpg.de     mk

  PLATFORMS: All.    

  HISTORY:
  8/23/02  awi		Created. 
 
*/

#include "Psych.h"

//file static variable definitions
static PsychFunctionPtr exitFunctionREGISTER = NULL;
static PsychFunctionPtr baseFunctionREGISTER = NULL;
static PsychFunctionTableEntry functionTableREGISTER[PSYCH_MAX_FUNCTIONS];
static char ModuleNameREGISTER[PSYCH_MAX_FUNCTION_NAME_LENGTH+1]; //+1 for term null
static char *currentFunctionNameREGISTER;
static int numFunctionsREGISTER = 0;
static psych_bool nameRegistered = FALSE;


//file static function declarations
static PsychError PsychRegisterModuleName(char *name);
static PsychError PsychRegisterBase(PsychFunctionPtr baseFunc);

/* PsychResetRegistry()
 *
 * Called by exit glue. Resets all registered functions to empty init default.
 *
 */
void PsychResetRegistry(void)
{
    exitFunctionREGISTER = NULL;
    baseFunctionREGISTER = NULL;
    ModuleNameREGISTER[0] = 0;
    currentFunctionNameREGISTER = NULL;
    numFunctionsREGISTER = 0;
    nameRegistered = FALSE;
    memset(&functionTableREGISTER[0], 0, sizeof(functionTableREGISTER));
}

/*  This function is called by the special subfunction 'DescribeModuleFunctionsHelper'.
 *  It returns the full table of registered subfunction names as an array.
 *  Useful for automatic documentation generators...
 *
 */
PsychError PsychDescribeModuleFunctions(void)
{
	static char useString[] = "subfunctionNames = Modulename('DescribeModuleFunctionsHelper' [, mode] [, subfunctionName]);";
	static char synopsisString[] = "Return a cell array of strings naming all subfunctions supported by this module if "
								   "the optional 'subfunctionName' argument is omitted. If 'subfunctionName' is a "
								   "string with a valid subfunction name for the module, and 'mode' is 1, return a 3 element cell "
								   "array of strings which describe the detailed syntax, help and see also strings "
								   "for that function - the text you'd get for Modulename('subfunctionName?'); ";
	static char seeAlsoString[] = "";

	char*						subfname;
    int							i;
    PsychGenericScriptType		*cellVector;
	PsychFunctionPtr			fcn;
	
    //all subfunctions should have these two lines.  
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    // Check for valid number of arguments
    PsychErrorExit(PsychCapNumInputArgs(2));   	
    PsychErrorExit(PsychCapNumOutputArgs(1)); 
    
	PsychCopyInIntegerArg(1, FALSE, &i);
	
	subfname = NULL;
	PsychAllocInCharArg(2, FALSE, &subfname);
	
	if (subfname) {
		// Find named subfunction in functin table:
		fcn = PsychGetProjectFunction(subfname);
		if (fcn) {
			// This will trigger the subfunction 'fcn' to execute its PsychGiveHelp() method, then
			// return. PsychGiveHelp() will return the help text to Runtime as cell array of strings,
			// instead of printing to runtimes console:
			PsychSetGiveHelp();
			PsychOneShotReturnHelp();
			(*fcn)();
			PsychClearGiveHelp();
		}
	}
	else {
		// Generate a cell array to return
		PsychAllocOutCellVector(1, FALSE, numFunctionsREGISTER,  &cellVector);
		for(i=0; i < numFunctionsREGISTER; i++) PsychSetCellVectorStringElement(i, functionTableREGISTER[i].name, cellVector);
	}
	        
    return(PsychError_none);
}

/*
	This function is called by the project to register project functions.
	
	See the documentation in the example RegisterProject.c file.  
	
	It registers this elements of the project:
	
	-subfunctions
	-the project base function
	-the project name

*/

//for use by projects
PsychError PsychRegister(char *name,  PsychFunctionPtr func)
{
	int i;

	//check to see if name is null which means we register the module base function.  
	if(name==NULL){
		if(func==NULL)
			return(PsychError_internal);
		return(PsychRegisterBase(func));
	}
	
	//check to see if the function is null which means we register the project name.
	if(func==NULL)
		return(PsychRegisterModuleName(name)); 
	

	//check to see if we have space left in the registry
	if(numFunctionsREGISTER == PSYCH_MAX_FUNCTIONS)
		return(PsychError_registerLimit);
		
	//check to see if the name has already been registered
	for(i=0;i<PSYCH_MAX_FUNCTIONS;i++){
		if(strcmp(name, functionTableREGISTER[i].name)==0)
			return(PsychError_registered);
	}
	
	//register the function and enable the subfunction dispatcher
	functionTableREGISTER[numFunctionsREGISTER].function = func;
	if(strlen(name) > PSYCH_MAX_FUNCTION_NAME_LENGTH)
		return(PsychError_longString);
	strcpy(functionTableREGISTER[numFunctionsREGISTER].name, name);
	++numFunctionsREGISTER;
	PsychEnableSubfunctions();
	return(PsychError_none);
}

//for use by projects
PsychError PsychRegisterExit(PsychFunctionPtr exitFunc)
{
	if(exitFunctionREGISTER == NULL){
		exitFunctionREGISTER = exitFunc;
		return(PsychError_none);
	}else
		return(PsychError_registered);
}

/* If null is passed then we return the base function, otherwise 
	we return the named subfuction or NULL if we can't find it.
	
	Also we store a pointer to the current function name string or 
	set that pointer to NULL if we are returning the project base function. 
*/
PsychFunctionPtr PsychGetProjectFunction(char *command)
{
	int i; 

	//return the project base function
	if(command==NULL){
		currentFunctionNameREGISTER = NULL;
		return(baseFunctionREGISTER);
	}
	// See if help is being requested
	if (command[strlen(command)-1] == '?') {
		PsychSetGiveHelp();
		command[strlen(command)-1]=0;
	}else
		PsychClearGiveHelp();
	
	//lookup the function in the table
	for(i=0;i<numFunctionsREGISTER;i++){
		if(PsychMatch(functionTableREGISTER[i].name, command)){
			currentFunctionNameREGISTER = functionTableREGISTER[i].name;
			return(functionTableREGISTER[i].function);
		}
	}

	// Unknown command.
	return NULL;
}

//for use by projects
char *PsychGetFunctionName(void)
{	
	static char noFunction[] = ""; //the string used to name the current subfunction when 
								   //there is none:

	if(currentFunctionNameREGISTER == NULL)
		return(noFunction);
	else
		return(currentFunctionNameREGISTER);
}

//for use by projects
char *PsychGetModuleName(void)
{
	return(ModuleNameREGISTER);
}

PsychFunctionPtr PsychGetProjectExitFunction(void)
{
	return(exitFunctionREGISTER);
}


//file static function definitions


/* 
	Store away the project name when the project init registers it. 
	This function is not called directly by the project init, it is 
	called by PsychRegister if PsychRegister figures out that the 
	project init is registering the project name.
*/
static PsychError PsychRegisterModuleName(char *name)
{
	if(nameRegistered)
		return(PsychError_registered);
	else{
		if(strlen(name) > PSYCH_MAX_FUNCTION_NAME_LENGTH)
			return(PsychError_longString);
		strcpy(ModuleNameREGISTER, name); 
		nameRegistered=TRUE;
	}
	return(PsychError_none);
}

/* 
	Store a pointer to the project base function in a static variable.
	This function is not called directly by the project init, it is 
	called by PsychRegister if PsychRegister figures out that the 
	project init is registering the project base function.
*/
 
static PsychError PsychRegisterBase(PsychFunctionPtr baseFunc)
{
	if(baseFunctionREGISTER != NULL)
		return(PsychError_registered);
	else
		baseFunctionREGISTER = baseFunc;
		
	return(PsychError_none);
}
