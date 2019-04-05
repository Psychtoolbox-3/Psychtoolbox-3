/*
    PsychToolbox3/Source/Common/PsychError.c

    AUTHORS:

        Allen.Ingling@nyu.edu           awi
        mario.kleiner.de@gmail.com      mk

    PLATFORMS:

        All.

    HISTORY:

        8/19/02  awi    Wrote it.
        11/16/04 awi    Added PsychErorr_argumentValueOutOfRange.

*/

#include "Psych.h"
#define MAX_PSYCH_ERRORS 100

//file static variables
static char *errorStringsERROR[MAX_PSYCH_ERRORS];
static char *errorNameStringsERROR[MAX_PSYCH_ERRORS];
static char *ArgTypeStringsERROR[PsychArgType_NUMTYPES];
static psych_bool usageErrorFlagsERROR[MAX_PSYCH_ERRORS];
static PsychArgDescriptorType specifiedArgERROR, receivedArgERROR;


//file static functions (these don't need the "Psych" prefix because they are file local)
static void InitErrorNameStrings(void);
static void InitErrorStrings(void);
static void InitUsageErrorFlags(void);
static void InitArgFormatTypeStrings(void);


PsychError InitPsychError(void)
{
    InitErrorNameStrings();
    InitErrorStrings();
    InitUsageErrorFlags();
    InitArgFormatTypeStrings();
    return(PsychError_none);
}


static void InitErrorNameStrings(void)
{
    errorNameStringsERROR[PsychError_none]=                             "PsychError_none";

    errorNameStringsERROR[PsychError_invalidArg_absent]=                "PsychError_invalidArg_absent";
    errorNameStringsERROR[PsychError_invalidArg_extra]=                 "PsychError_invalidArg_extra";
    errorNameStringsERROR[PsychError_invalidArg_type]=                  "PsychError_invalidArg_type";
    errorNameStringsERROR[PsychError_invalidArg_size]=                  "PsychError_invalidArg_size";

    errorNameStringsERROR[PsychError_extraInputArg]=                    "PsychError_extraInputArg";
    errorNameStringsERROR[PsychError_missingInputArg]=                  "PsychError_missingInputArg";
    errorNameStringsERROR[PsychError_extraOutputArg]=                   "PsychError_extraOutputArg";
    errorNameStringsERROR[PsychError_missingOutputArg]=                 "PsychError_extraOutputArg";

    errorNameStringsERROR[PsychError_toomanyWin]=                       "PsychError_toomanyWin";
    errorNameStringsERROR[PsychError_outofMemory]=                      "PsychError_outofMemory";
    errorNameStringsERROR[PsychError_scumberNotWindex]=                 "PsychError_scumberNotWindex";
    errorNameStringsERROR[PsychError_windexNotScumber]=                 "PsychError_windexNotScumber";
    errorNameStringsERROR[PsychError_invalidIntegerArg]=                "PsychError_invalidIntegerArg";
    errorNameStringsERROR[PsychError_invalidWindex]=                    "PsychError_invalidWindex";
    errorNameStringsERROR[PsychError_invalidScumber]=                   "PsychError_invalidScumber";
    errorNameStringsERROR[PsychError_invalidNumdex]=                    "PsychError_invalidNumdex";
    errorNameStringsERROR[PsychError_invalidColorArg]=                  "PsychError_invalidColorArg";
    errorNameStringsERROR[PsychError_invalidDepthArg]=                  "PsychError_invalidDepthArg";
    errorNameStringsERROR[PsychError_invalidRectArg]=                   "PsychError_invalidRectArg";
    errorNameStringsERROR[PsychError_invalidNumberBuffersArg]=          "PsychError_invalidNumberBuffersArg";
    errorNameStringsERROR[PsychError_nullWinRecPntr]=                   "PsychError_nullWinRecPntr";
    errorNameStringsERROR[PsychError_registerLimit]=                    "PsychError_registerLimit";
    errorNameStringsERROR[PsychError_registered]=                       "PsychError_registered";
    errorNameStringsERROR[PsychError_longString]=                       "PsychError_longString";
    errorNameStringsERROR[PsychError_longStringPassed]=                 "PsychError_longStringPassed";
    errorNameStringsERROR[PsychError_unimplemented]=                    "PsychError_unimplemented";
    errorNameStringsERROR[PsychError_internal]=                         "PsychError_internal";
    errorNameStringsERROR[PsychError_system]=                           "PsychError_system";
    errorNameStringsERROR[PsychError_invalidArgRef]=                    "Psycherror_invalidArgRef";
    errorNameStringsERROR[PsychError_OpenGL]=                           "Psycherror_OpenGL";
    errorNameStringsERROR[PsychError_SDL]=                              "PsychError_SDL";
    errorNameStringsERROR[PsychError_SurfaceLockFailed]=                "PsychError_SurfaceLockFailed";
    errorNameStringsERROR[PsychError_SurfaceAlreadyLocked]=             "PsychError_SurfaceAlreadyLocked";
    errorNameStringsERROR[PsychError_InvalidWindowRecord]=              "PsychError_InvalidWindowRecord";
    errorNameStringsERROR[PsychError_unsupportedVideoMode]=             "PsychError_unsupportedVideoMode";
    errorNameStringsERROR[PsychError_user]=                             "PsychError_generalUser";
    errorNameStringsERROR[PsychError_unrecognizedPreferenceName]=       "PsychError_unrecognizedPreferenceName";
    errorNameStringsERROR[PsychError_unsupportedOS9Preference]=         "PsychError_unsupportedOS9Preference";
    errorNameStringsERROR[PsychError_inputMatrixIllegalDimensionSize]=  "PsychError_inputMatrixIllegalDimensionSize";
    errorNameStringsERROR[PsychError_stringOverrun]=                    "PsychError_stringOverrun";
    errorNameStringsERROR[PsychErorr_argumentValueOutOfRange]=          "PsychErorr_argumentValueOutOfRange";
}


static void InitErrorStrings(void)
{
    errorStringsERROR[PsychError_none]=                                 "No Error";

    errorStringsERROR[PsychError_invalidArg_absent]=                    "Missing argument";
    errorStringsERROR[PsychError_invalidArg_extra]=                     "Extra argument supplied";
    errorStringsERROR[PsychError_invalidArg_type]=                      "Incorrect argument type supplied";
    errorStringsERROR[PsychError_invalidArg_size]=                      "Wrong size argument supplied";

    errorStringsERROR[PsychError_extraInputArg]=                        "Extra input argument described";
    errorStringsERROR[PsychError_missingInputArg]=                      "Missing input argument";
    errorStringsERROR[PsychError_extraOutputArg]=                       "Extra output argument supplied";
    errorStringsERROR[PsychError_missingOutputArg]=                     "Missing output argument";

    errorStringsERROR[PsychError_toomanyWin]=                           "Exceeded window Limit";
    errorStringsERROR[PsychError_outofMemory]=                          "Out of memory";
    errorStringsERROR[PsychError_scumberNotWindex]=                     "Window index expected, screen number received";
    errorStringsERROR[PsychError_windexNotScumber]=                     "Screen number expected, window index received";
    errorStringsERROR[PsychError_invalidIntegerArg]=                    "Supplied value was either fractional or outside of allowed range.";
    errorStringsERROR[PsychError_invalidWindex]=                        "Invalid Window (or Texture) Index provided: It doesn't correspond to an open window or texture.\nDid you close it accidentally via Screen('Close') or Screen('CloseAll') ?";
    errorStringsERROR[PsychError_invalidScumber]=                       "Invalid Screen Number";
    errorStringsERROR[PsychError_invalidNumdex]=                        "Invalid or missing Screen Number or Window Index";
    errorStringsERROR[PsychError_invalidColorArg]=                      "Invalid Color Argument";
    errorStringsERROR[PsychError_invalidDepthArg]=                      "Invalid pixel size argument";
    errorStringsERROR[PsychError_invalidRectArg]=                       "Invalid Rect Argument";
    errorStringsERROR[PsychError_invalidNumberBuffersArg]=              "Invalid number of buffers specified";
    errorStringsERROR[PsychError_nullWinRecPntr]=                       "Window Record Pointer is NULL";
    errorStringsERROR[PsychError_registerLimit]=                        "Too many subfunctions registered";
    errorStringsERROR[PsychError_registered]=                           "Attempt to reregister a function";
    errorStringsERROR[PsychError_longString]=                           "String exceeds specified length limit";
    errorStringsERROR[PsychError_longStringPassed]=                     "String exceeds specified length limit";
    errorStringsERROR[PsychError_unimplemented]=                        "Feature not yet implemented";
    errorStringsERROR[PsychError_internal]=                             "Unspecified error, probably a Psychtoolbox bug";
    errorStringsERROR[PsychError_system]=                               "Error reported by a system call";
    errorStringsERROR[PsychError_invalidArgRef]=                        "Reference to nonexistent argument";
    errorStringsERROR[PsychError_OpenGL]=                               "OpenGL call returned an error";
    errorStringsERROR[PsychError_SDL]=                                  "SDL call failed.";
    errorStringsERROR[PsychError_SurfaceLockFailed]=                    "Could not lock screen surface.";
    errorStringsERROR[PsychError_SurfaceAlreadyLocked]=                 "Attempt to lock screen surface while it was already locked.";
    errorStringsERROR[PsychError_InvalidWindowRecord]=                  "An Invalid window record was referenced.";
    errorStringsERROR[PsychError_unsupportedVideoMode]=                 "The specified video mode is not supported by the specified display.";
    errorStringsERROR[PsychError_user]=                                 "Usage error";
    errorStringsERROR[PsychError_unrecognizedPreferenceName]=           "Unrecognized preference name";
    errorStringsERROR[PsychError_unsupportedOS9Preference]=             "OS 9 preference unsupported on OS X";
    errorStringsERROR[PsychError_inputMatrixIllegalDimensionSize]=      "Supplied matrix argument is the wrong size";
    errorStringsERROR[PsychError_stringOverrun]=                        "String longer than allocated storage space";
    errorStringsERROR[PsychErorr_argumentValueOutOfRange]=              "A numerical argument value was outside of the allowable range.";
}

/*
    All errors are classifed as either internal or usage errors.  Internal errors are caused
    presumably by a bug in the Psychtoolbox or project, usage errors are cuased by the user
    passing a bogus argument.  The classification determines the style in which PsychErrorExitMsg
    reports the bug.
*/
#define            kPsychUserError                    TRUE
#define            kPsychBug                        FALSE

static void InitUsageErrorFlags(void)
{

    //TRUE means that the user did something wrong.
    //FALSE means that it's an internal Psychtoolbox or OS error.
    usageErrorFlagsERROR[PsychError_none]=                                  kPsychBug;

    usageErrorFlagsERROR[PsychError_invalidArg_absent]=                     kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidArg_extra]=                      kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidArg_type]=                       kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidArg_size]=                       kPsychUserError;

    usageErrorFlagsERROR[PsychError_extraInputArg]=                         kPsychUserError;
    usageErrorFlagsERROR[PsychError_missingInputArg]=                       kPsychUserError;
    usageErrorFlagsERROR[PsychError_extraOutputArg]=                        kPsychUserError;
    usageErrorFlagsERROR[PsychError_missingOutputArg]=                      kPsychUserError;

    usageErrorFlagsERROR[PsychError_toomanyWin]=                            kPsychBug;
    usageErrorFlagsERROR[PsychError_outofMemory]=                           kPsychBug;
    usageErrorFlagsERROR[PsychError_scumberNotWindex]=                      kPsychUserError;
    usageErrorFlagsERROR[PsychError_windexNotScumber]=                      kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidWindex]=                         kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidIntegerArg]=                     kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidScumber]=                        kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidNumdex]=                         kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidColorArg]=                       kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidDepthArg]=                       kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidRectArg]=                        kPsychUserError;
    usageErrorFlagsERROR[PsychError_invalidNumberBuffersArg]=               kPsychUserError;
    usageErrorFlagsERROR[PsychError_nullWinRecPntr]=                        kPsychBug;
    usageErrorFlagsERROR[PsychError_registerLimit]=                         kPsychBug;
    usageErrorFlagsERROR[PsychError_registered]=                            kPsychBug;
    usageErrorFlagsERROR[PsychError_longString]=                            kPsychBug;
    usageErrorFlagsERROR[PsychError_longStringPassed]=                      kPsychUserError;
    usageErrorFlagsERROR[PsychError_unimplemented]=                         kPsychBug;
    usageErrorFlagsERROR[PsychError_internal]=                              kPsychBug;
    usageErrorFlagsERROR[PsychError_internal]=                              kPsychBug;
    usageErrorFlagsERROR[PsychError_invalidArgRef]=                         kPsychBug;
    usageErrorFlagsERROR[PsychError_OpenGL]=                                kPsychBug;
    usageErrorFlagsERROR[PsychError_SDL]=                                   kPsychBug;
    usageErrorFlagsERROR[PsychError_SurfaceLockFailed]=                     kPsychBug;
    usageErrorFlagsERROR[PsychError_SurfaceAlreadyLocked]=                  kPsychBug;
    usageErrorFlagsERROR[PsychError_InvalidWindowRecord]=                   kPsychBug;
    usageErrorFlagsERROR[PsychError_unsupportedVideoMode]=                  kPsychUserError;
    usageErrorFlagsERROR[PsychError_user]=                                  kPsychUserError;
    usageErrorFlagsERROR[PsychError_unrecognizedPreferenceName]=            kPsychUserError;
    usageErrorFlagsERROR[PsychError_unsupportedOS9Preference]=              kPsychUserError;
    usageErrorFlagsERROR[PsychError_inputMatrixIllegalDimensionSize]=       kPsychUserError;
    usageErrorFlagsERROR[PsychError_stringOverrun]=                         kPsychBug;
    usageErrorFlagsERROR[PsychErorr_argumentValueOutOfRange]=               kPsychUserError;
}


/*
    TO DO: These names depend on the scripting language and should be moved to ScriptingGlue
*/
static void InitArgFormatTypeStrings(void)
{
    ArgTypeStringsERROR[0] =    "unspecified";
    ArgTypeStringsERROR[1] =    "unclassified";
    ArgTypeStringsERROR[2] =    "char";
    ArgTypeStringsERROR[3] =    "uint8";
    ArgTypeStringsERROR[4] =    "uint16";
    ArgTypeStringsERROR[5] =    "uint32";
    ArgTypeStringsERROR[6] =    "int8";
    ArgTypeStringsERROR[7] =    "int16";
    ArgTypeStringsERROR[8] =    "int32";
    ArgTypeStringsERROR[9] =    "double";
    ArgTypeStringsERROR[10] =   "boolean";
    ArgTypeStringsERROR[11] =   "struct";
    ArgTypeStringsERROR[12] =   "cell";
    ArgTypeStringsERROR[13] =   "single";
    ArgTypeStringsERROR[14] =   "uint64";
    ArgTypeStringsERROR[15] =   "default";
}


/*
    format specifiers can consist of multiple type specifications when specifing
    allowed inputs.  This function returns the number of types which compose a
    type specification and  a list of the strings naming those types.
*/
int PsychDecomposeArgFormat(PsychArgFormatType argType,  const char **names)
{
    int i, numTypes=0;
    PsychArgFormatType strIndex;
      psych_bool isElement;

    //we are given a composite of types so find which to name in the string
    for(i=0;i<PsychArgType_NUMTYPES;i++){
        strIndex = (PsychArgFormatType)pow(2,i);
        isElement = (psych_bool)(((int)argType & (int)strIndex) > 0);
        if(isElement){
            if(names!=NULL) //in case we only want a count and pass NULL for names
                names[numTypes]=ArgTypeStringsERROR[i];
            ++numTypes;
        }
    }
    return numTypes;
}

/*
    If there is no error then just return().  If there is an error then
    print the function name which should have been pushed using PsychPushHelp.
    and also print the string describing the type of error and then return
    control to the scripting environment.

    extraErrorString can provide addional debugging info.

    TO DO:
    See "to do" list above.  This could be improved to print more specific error
    information, for example whith an invalid argument it could print which argument
    was invalid or which was missing.

    PsychErrorExitMsg should be renamed PsychErrorExitC and accept another argument which
    is the string naming the file or function where the error occurred. A macro named PsychErrorExitMsg
    should be provided which calls PsychErrorExitMsg and plugs in the ANSI predifined macros giving the
    file and function where the error occurred.  PsychErrorExitC only prints that if the error type
    is internal.

*/
void PsychErrorExitC(PsychError error,
                     const char *extraErrorString,
                     int lineNum,
                     const char *funcName,
                     const char *fileName)
{
    char *functionName;
    PsychArgDescriptorType *specified, *received;
    int i, numTypes;
    const char *typeStrings[PsychArgType_NUMTYPES];

    // This function must be implemented by the scripting glue and allows to
    // process the error in some frontend specific way, e.g., set exception
    // state, or simply ignore it:
    PsychProcessErrorInScripting(error, extraErrorString ? extraErrorString : errorStringsERROR[error]);

    //if the error is type none then just return
    if(error==PsychError_none)
        return;

    //if the error is an internal error then display copious info
    if(!usageErrorFlagsERROR[error]){
        printf("INTERNAL PSYCHTOOLBOX ERROR\n");
        printf("\terror:                %s\n",errorNameStringsERROR[error]);
        printf("\tgeneral description:  %s\n",errorStringsERROR[error]);
        if(extraErrorString != NULL)
            printf("\tspecific description: %s\n",extraErrorString);
        printf("\tmodule name:          %s\n",PsychGetModuleName());
        printf("\tsubfunction call:     %s\n",PsychGetFunctionName());
        printf("\tfile name:            %s\n",fileName);
        printf("\tfunction name:        %s\n",funcName);
        printf("\tline number:          %d\n",lineNum);
        PsychErrMsgTxt(NULL);  //exit the module
    }else{
        //if the error is usage error then
        functionName = PsychGetFunctionName();
        printf("Error in function %s: ",functionName);
        printf("\t%s\n",errorStringsERROR[error]);
        if(extraErrorString !=NULL)
            printf("%s\n",extraErrorString);
        //attempts to fetch or return arguments cache descriptions of the the desired
        //and returned values in static variables.  If attempts to set or retrieve an
        //argument resulted in an error value PsychError_invalidArg_*, then we can retrieve that info
        //about the specified and desired argument to generate a detailed error message.
        if(error == PsychError_invalidArg_absent || error == PsychError_invalidArg_extra || error == PsychError_invalidArg_type || error == PsychError_invalidArg_size){
            PsychGetArgDescriptor(&specified, &received);
            if(specified != NULL && received != NULL){  //why would these be null ?
                //for specified
                printf("\tDiscrepancy between a specified and supplied argument:\n");
                printf("\t\tSpecified Argument Description:\n");
                printf("\t\t\tdirection: %s\n", (specified->direction==PsychArgIn) ? "Input" : "Output");  //input  or output
                printf("\t\t\tposition: %d\n", specified->position);  //position
                printf("\t\t\tpresence:");
                switch(specified->isThere){
                    case kPsychArgAbsent: printf("forbidden");break;
                    case kPsychArgPresent: printf("required");break;
                    case kPsychArgFixed: printf("fixed");break;
                }
                printf("\n");

                printf("\t\t\tformats allowed:\n");
                numTypes=PsychDecomposeArgFormat(specified->type, typeStrings);
                for(i=0;i<numTypes;i++)
                    printf("\t\t\t\t%s\n",typeStrings[i]);

                printf("\t\t\tminimum M: %s\n", int2str(specified->mDimMin));
                printf("\t\t\tmaximum M: %s\n", (specified->mDimMax==kPsychUnboundedArraySize) ? "unbounded" : int2str(specified->mDimMax));
                printf("\t\t\tminimum N: %s\n", int2str(specified->nDimMin));
                printf("\t\t\tmaximum N: %s\n", (specified->nDimMax==kPsychUnboundedArraySize) ? "unbounded" : int2str(specified->nDimMax));
                printf("\t\t\tminimum P: %s\n", int2str(specified->pDimMin));
                printf("\t\t\tmaximum P: %s\n", (specified->pDimMax==kPsychUnboundedArraySize) ? "unbounded" : int2str(specified->pDimMax));

                //For received
                printf("\t\tProvided Argument Description:\n");
                printf("\t\t\tdirection: %s\n", (received->direction== PsychArgIn) ? "Input" : "Output");  //input  or output
                printf("\t\t\tposition: %d\n", received->position);  //position
                printf("\t\t\tpresence: ");
                switch(received->isThere){
                    case kPsychArgAbsent: printf("absent\n");break;
                    case kPsychArgPresent: printf("present\n");break;
                    case kPsychArgFixed: printf("fixed\n");break;
                }
                if(received->isThere!=kPsychArgAbsent && received->direction==PsychArgIn){
                    printf("\t\t\tformat:");

                    numTypes=PsychDecomposeArgFormat(received->type, typeStrings); //there should only be one, but check for bugs
                    for(i=0;i<numTypes;i++)
                        printf("%s\n",typeStrings[i]);
                    printf("\t\t\t number of dimensions: %d\n", received->numDims);
                    printf("\t\t\t\t M: %s\n", int2str(received->mDimMin));
                    printf("\t\t\t\t N: %s\n", int2str(received->nDimMin));
                    printf("\t\t\t\t P: %s\n", int2str(received->pDimMin));
                }
            }
        }
        PsychGiveUsageExit();
    }
}


/*
    copy out from structures specified by pointers and store in static variable
*/
PsychError PsychStoreArgDescriptor(PsychArgDescriptorType *specified,
                                   PsychArgDescriptorType *received)
{
    if(specified != NULL)
        memcpy(&specifiedArgERROR, specified, sizeof(PsychArgDescriptorType));

    if(received != NULL)
        memcpy(&receivedArgERROR, received, sizeof(PsychArgDescriptorType));

    return(PsychError_none);
}


/*
    set referenced pointers to the structures
*/
PsychError PsychGetArgDescriptor(PsychArgDescriptorType **specified,
                                 PsychArgDescriptorType **received)
{
    if(specified!=NULL)
        *specified = &specifiedArgERROR;
    if(received!=NULL)
        *received = &receivedArgERROR;

    return(PsychError_none);
}


/*
    PsychDumpArgDescriptors has not yet been tested.

    Note that argument descriptors are described in three dimensions:
    1- Specified or Received.
    2- Input or Output
    3- Position
*/
void PsychDumpArgDescriptors(void)
{
    printf("\n");
    printf("__________________________________________");
    /* display the specified argument descriptor. */
    printf("SPECIFIED ARGUMENT DESCRIPTOR:\n");
    printf("position:     %d\n", specifiedArgERROR.position);
    printf("direction:      ");
    switch(specifiedArgERROR.direction){
        case PsychArgIn: printf("PsychArgIn\n"); break;
        case PsychArgOut: printf("PsychArgOut\n"); break;
    }
    printf("isThere:      ");
    switch(specifiedArgERROR.isThere){
        case kPsychArgAbsent: printf("kPsychArgAbsent\n"); break;
        case kPsychArgPresent: printf("kPsychArgPresent\n"); break;
        case kPsychArgFixed: printf("kPsychArgFixed\n"); break;
    }

    printf("type:      %d\n", (int)specifiedArgERROR.type); //this should be fixed to display a type name
    printf("mDimMin:      %s\n", int2str(specifiedArgERROR.mDimMin));
    printf("mDimMax:      %s\n", int2str(specifiedArgERROR.mDimMax));
    printf("nDimMin:      %s\n", int2str(specifiedArgERROR.nDimMin));
    printf("nDimMax:      %s\n", int2str(specifiedArgERROR.nDimMax));
    printf("pDimMin:      %s\n", int2str(specifiedArgERROR.pDimMin));
    printf("pDimMax:      %s\n", int2str(specifiedArgERROR.pDimMax));

    /* display the received argument descriptor. */
    printf("RECEIVED ARGUMENT DESCRIPTOR:\n");
    printf("position:     %d\n", receivedArgERROR.position);
    printf("direction:      ");
    switch(receivedArgERROR.direction){
        case PsychArgIn: printf("PsychArgIn\n"); break;
        case PsychArgOut: printf("PsychArgOut\n"); break;
    }

    printf("isThere:      ");
    switch(receivedArgERROR.isThere){
        case kPsychArgAbsent: printf("kPsychArgAbsent\n"); break;
        case kPsychArgPresent: printf("kPsychArgPresent\n"); break;
        case kPsychArgFixed: printf("kPsychArgFixed\n"); break;
    }

    if((receivedArgERROR.isThere==kPsychArgPresent) && (receivedArgERROR.direction==PsychArgIn)){
        printf("type:      %d\n", (int)receivedArgERROR.type);
        printf("numDims:      %s\n", int2str(receivedArgERROR.numDims));
        printf("mDimMin:      %s\n", int2str(receivedArgERROR.mDimMin));
        printf("mDimMax:      %s\n", int2str(receivedArgERROR.mDimMax));
        printf("nDimMin:      %s\n", int2str(receivedArgERROR.nDimMin));
        printf("nDimMax:      %s\n", int2str(receivedArgERROR.nDimMax));
        printf("pDimMin:      %s\n", int2str(receivedArgERROR.pDimMin));
        printf("pDimMax:      %s\n", int2str(receivedArgERROR.pDimMax));
    }

    printf("__________________________________________");
    printf("\n");
}
