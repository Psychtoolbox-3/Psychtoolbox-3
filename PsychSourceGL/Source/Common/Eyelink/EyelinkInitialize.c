/*
 
	PsychSourceGL/Source/Common/Eyelink/EyelinkInitialize.c
 
	PROJECTS: Eyelink 
 
	AUTHORS:
	 cburns@berkeley.edu			cdb
	 E.Peters@ai.rug.nl				emp
	 f.w.cornelissen@med.rug.nl		fwc
 
	PLATFORMS:	All.
 
	HISTORY:
 
	 11/23/05  cdb		Created.
	 28/06/06	fwc		added EyelinkInitializeDummy function. Unlike in the OS9 version, 
						I decided to make it a seperate function rather than option in EyelinkInitialize.
						Eyelink manual actually seems to advice to use eyelink_open_connection()
						rather than the way we do it now.
	 29/06/06	fwc		fixed EyelinkInitializeDummy (didn't set giSystemInitialized)
	 30-10-06	fwc		no longer use PsychErrorExitMsg to report error, otherwise we cannot connect in dummy mode later on
	 4-4-09     mk+edf  added registration of callback
 
	TARGET LOCATION:
 
	 Eyelink.mexmac resides in:
	 PsychHardware/EyelinkToolbox/EyelinkBasic
 
	NOTES:
	 This is a scaled down version compared to the OS9 version that Enno and Frans wrote.
	 I hardcode in the buffersize and ignore other params.

*/

#include "PsychEyelink.h"

char useString[] = "[status =] Eyelink('Initialize' [, displayCallbackFunction])";
char useDummyString[] = "[status =] Eyelink('InitializeDummy' [, displayCallbackFunction])";
char useTestString[] = "Eyelink('TestSuite')";

static char synopsisString[] = 
"Initializes Eyelink and Ethernet system.\n"
"Opens tracker connection, reports any problems.\n"
"The optional argument 'displayCallbackFunction' registers the name "
"of an m-file on the path that will handle callbacks from eyelink "
"(typically 'PsychEyelinkDispatchCallback' to display the camera "
"image in PTB). No callbacks will be generated if that argument "
"is omitted.\n"
"Returns: 0 if OK, -1 if error";

static char synopsisDummyString[] = 
"Initializes Eyelink in dummy mode, useful for debugging\n"
"The optional argument 'displayCallbackFunction' registers the name "
"of an m-file on the path that will handle callbacks from eyelink "
"(typically 'PsychEyelinkDispatchCallback' to display the camera "
"image in PTB). No callbacks will be generated if that argument "
"is omitted.\n"
"Returns: 0 if OK, -1 if error";

static char synopsisTestString[] = 
"Perform test routines, e.g., display of synthetic test image, "
"exercising the callback functions.\n";

static char seeAlsoString[] = "";

#define ERR_BUFF_LEN 1000

// Check if optional callback passed, and if so, if it is valid.
// If passed and valid, enable callbacks:
void CheckAndAssignCallback(int argpos)
{
	char    errMsg[ERR_BUFF_LEN];
	char*   callbackString;
	PsychGenericScriptType	*input[1];
	PsychGenericScriptType	*output[1];
	
	// Initialize graphics callbacks for eye camera et al.:
	if (PsychAllocInCharArg(argpos, FALSE, &callbackString) && strlen(callbackString) > 0) {
		// Callback string passed. Check if it corresponds to a valid function on
		// Matlab's path:
		input[0] = mxCreateString(callbackString);
		output[0]=NULL;

		// mexCallMatlab should be safe here, as not called from within eyelink callbacks:
		if(mexCallMATLAB(1, output, 1, input, "exist")) {
			mxDestroyArray(input[0]);
			PsychErrorExitMsg(PsychError_system, "Fatal error calling runtime's 'exist' function!");
		}

		mxDestroyArray(input[0]);
		
		if(mxGetScalar(output[0]) <= 0) {
			mxDestroyArray(output[0]);
			snprintf(errMsg, ERR_BUFF_LEN, "Eyelink: The provided callback function '%s' is not a defined function, MEX file or M-File on your path as required!", callbackString);
			PsychErrorExitMsg(PsychError_user, errMsg);
		}

		mxDestroyArray(output[0]);
		
		// Everything good: Initialize callbacks and enable them:
		PsychEyelink_init_core_graphics(callbackString);
	}
	else {
		// No callback passed: Disable use of callbacks:
		PsychEyelink_uninit_core_graphics();
    }
	
	return;
}

PsychError EyelinkInitialize(void)
{
	int		iBufferSize = 20000; // Hardcode default
	int		iStatus		= -1;
	
	// Add help strings
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	
	// Output help if asked
	if(PsychIsGiveHelp()) {
		PsychGiveHelp();
		return(PsychError_none);
	}
	
	// Check arguments
	PsychErrorExit(PsychCapNumInputArgs(1));
	PsychErrorExit(PsychRequireNumInputArgs(0));
	PsychErrorExit(PsychCapNumOutputArgs(1));
	
	if(giSystemInitialized == 1) {
		// Is there a psych function for this?
		mexPrintf("Eyelink already initialized!\n");
		iStatus = 0;
	} else {
		
		iStatus = (int) open_eyelink_system((UINT16) iBufferSize, NULL);
		// If link is established, open connection
		if(iStatus != 0) {
			msec_delay(300);
			iStatus = eyelink_open();
		}
		
		// Check for errors and report
		if(iStatus != 0) {
			close_eyelink_system();
			if(iStatus == CONNECT_TIMEOUT_FAILED) {
			//	PsychErrorExitMsg(PsychError_user, "EyeLink: Initialize:  Cannot connect to EyeLink\n");
				mexPrintf("EyeLink: Initialize:  Cannot connect to EyeLink\n");
			} else if(iStatus == WRONG_LINK_VERSION) {
			//	PsychErrorExitMsg(PsychError_user, "EyeLink: Initialize:  Link version mismatch to EyeLink tracker\n");
				mexPrintf( "EyeLink: Initialize:  Link version mismatch to EyeLink tracker\n");
			} else if(iStatus == LINK_INITIALIZE_FAILED) {
			//	PsychErrorExitMsg(PsychError_user, "EyeLink: Initialize:  Cannot initialize link\n");
				mexPrintf( "EyeLink: Initialize:  Cannot initialize link\n");
			}
		} else {
			giSystemInitialized = 1;
		}
	}

	// Check for optional callbackString as argument 1, sanity check and enable callbacks, if provided:
	CheckAndAssignCallback(1);

	// Copy output arg
	PsychCopyOutDoubleArg(1, FALSE, iStatus);
	
	return(PsychError_none);
	
}


PsychError EyelinkInitializeDummy(void)
{
	int		iStatus		= -1;
	
	// Add help strings
	PsychPushHelp(useDummyString, synopsisDummyString, seeAlsoString);
	
	// Output help if asked
	if(PsychIsGiveHelp()) {
		PsychGiveHelp();
		return(PsychError_none);
	}
	
	// Check arguments
	PsychErrorExit(PsychCapNumInputArgs(1));
	PsychErrorExit(PsychRequireNumInputArgs(0));
	PsychErrorExit(PsychCapNumOutputArgs(1));
	
	if(giSystemInitialized == 1) {
		// Is there a psych function for this?
		mexPrintf("Eyelink already initialized!\n");
		iStatus = 0;
	} else {
		
		mexPrintf("Eyelink: Opening Eyelink in DUMMY mode\n");

		#if PSYCH_SYSTEM == PSYCH_WINDOWS
			// Need to open an eyelink connection on Windows, otherwise the beast will crash!
			iStatus = (int) open_eyelink_connection(-1);
			if (iStatus == 0) iStatus = (int) eyelink_dummy_open();
		#else
			iStatus = (int) eyelink_dummy_open();
		#endif

		// Check for errors and report
		if(iStatus != 0) {
			close_eyelink_system();
			if(iStatus == CONNECT_TIMEOUT_FAILED) {
				PsychErrorExitMsg(PsychError_user, "EyeLink: InitializeDummy:  Cannot connect to EyeLink\n");
			} else if(iStatus == WRONG_LINK_VERSION) {
				PsychErrorExitMsg(PsychError_user, "EyeLink: InitializeDummy:  Link version mismatch to EyeLink tracker\n");
			} else if(iStatus == LINK_INITIALIZE_FAILED) {
				PsychErrorExitMsg(PsychError_user, "EyeLink: InitializeDummy:  Cannot initialize link\n");
			}
		} else {
			giSystemInitialized = 1;
		}				
	}

	// Check for optional callbackString as argument 1, sanity check and enable callbacks, if provided:
	CheckAndAssignCallback(1);
	
	// Copy output arg
	PsychCopyOutDoubleArg(1, FALSE, iStatus);

	return(PsychError_none);	
}

PsychError EyelinkTestSuite(void)
{
	// Add help strings
	PsychPushHelp(useTestString, synopsisTestString, seeAlsoString);
	
	// Output help if asked
	if(PsychIsGiveHelp()) {
		PsychGiveHelp();
		return(PsychError_none);
	}
	
	// Check arguments
	PsychErrorExit(PsychCapNumInputArgs(0));
	PsychErrorExit(PsychRequireNumInputArgs(0));
	PsychErrorExit(PsychCapNumOutputArgs(0));
	
	// Execute internal tests:
	PsychEyelink_TestEyeImage();
	
	return(PsychError_none);	
}

PsychError EyelinkSetAddress(void)
{
    static char useAddressString[] = "[status =] Eyelink('SetAddress', ipaddress);";
    static char synopsisAddressString[] =
    "Set IP address of eyelink tracker computer to connect to.\n"
    "Call this before opening a tracker connection if you want to use "
    "a non-default IP address for the tracker computer.\n"
    "'ipaddress' is a character string containing the tracker address.\n"
    "Returns: 0 if OK, -1 on error\n";
    static char seeAlsoAddressString[] = "Initialize";
    
    char* ipAddress;
    int	iStatus = -1;
    
    // Add help strings:
    PsychPushHelp(useAddressString, synopsisAddressString, seeAlsoAddressString);
    
    // Output help if asked:
    if (PsychIsGiveHelp()) {
        PsychGiveHelp();
        return(PsychError_none);
    }
    
    // Check arguments:
    PsychErrorExit(PsychCapNumInputArgs(1));
    PsychErrorExit(PsychRequireNumInputArgs(1));
    PsychErrorExit(PsychCapNumOutputArgs(1));
    
    // Must call this function before opening a tracker connection:
    if (giSystemInitialized == 1) PsychErrorExitMsg(PsychError_user, "Tried to set target tracker host ip address, but connection already initialized!");
    
    // Get mandatory string with IP address spec:
    PsychAllocInCharArg(1, kPsychArgRequired, &ipAddress);
    
    // Assign new address:
    iStatus = set_eyelink_address(ipAddress);
    
    // Copy out status:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, iStatus);
    
    return(PsychError_none);
}
