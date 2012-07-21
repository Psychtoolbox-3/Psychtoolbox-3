/*
	SCREENHookFunction.c	
  
    AUTHORS:
    
		Mario.Kleiner@tuebingen.mpg.de  mk
  
    PLATFORMS:	
	
		All.
    
    HISTORY:
    
		12/05/06	mk		Wrote it.
	
    DESCRIPTION:

		Manage special callback functions (processing hooks). Screen allows to "hook" specific
		OpenGL GLSL shaders, C callback functions or Matlab M-File functions into its internal
		processing chain at well defined points in a (hopefully) well defined manner. This allows
		to extend or customize Screens behaviour. Useful applications: Automatic transparent
		on-the-fly image processing on stimuli just before stimulus onset, e.g., filtering.
		Customization of drawing commands. Advanced stereo stimulus output algorithms. Support for
		special output devices like BrightSide HDR, Bits++, Pelli & Zhang attenuator..., interfacing
		with 3rd party hardware like trigger systems...
		
		We provide well defined hook points, each one with a unique descriptive name. Each hook is
		associated with a queue or chain of callback functions which is executed from the beginning
		to the end whenever that hook is processed. A callback can be a GLSL OpenGL shader program
		for image processing, a Matlab M-Function for high-level processing or a C callback function
		for complex low-level or realtime operations like triggering some output device or similar.
		
		This Screen subfunction provides subfunctions to add callbacks to the chain, reset the chain,
		enable or disable the chain, query properties of single hooks or dump the whole chain state
		to the command window for debug purpose.
  
    NOTES:

    TO DO:

*/

#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "[ret1, ret2, ...] = Screen('HookFunction', windowPtr, 'Subcommand', 'HookName', arg1, arg2, arg3, arg4);";
//																		1           2             3          4     5     6     7
static char synopsisString[] = 
    "Manage Screen processing hook chains. Hook chains are a way to extend PTBs behaviour with plugins for generic processing or "
	"fast image processing. They should allow to interface PTB seamlessly with 3rd party special response collection devices or "
	"stimulus display devices. They also allow to perform high performance image processing, utilizing modern graphics hardwares "
	"huge computational bandwidth. Read 'help PsychGLImageProcessing' for more infos."
	"\n\n"
	"\nSubsubcommands and their syntax: \n\n"
	"Screen('HookFunction', [windowPtr], 'ListAll'); \n"
	"Print out a listing of all supported 'HookName' hook points, including a short synopsis on what they do, to the Matlab console "
	"in a human readable format - Useful as a reference for coding. "
	"The 'windowPtr' argument is optional. If left out (replaced by empty [] brackets) it prints all implemented chains."
	"\n\n"
	"[slot idstring blittercfg voidptr glslid luttexid insertString] = Screen('HookFunction', windowPtr, 'Query', hookname, slotnameOrIndex);\n"
	"Query information about a specific command slot in a specific hook processing chain: "
	"'hookname' is the name of the chain to query, e.g., 'StereoCompositingBlit' for the stereo processing chain. Use the subcommand "
	"'ListAll' for a printout of all available processing hooks and a short help on them. "
	"'slotnameOrIndex', either the symbolic name of the requested slot or the index in the hook-chain: Indices start with 0. Names can "
	"be assigned to slots when you add a processing slot with the 'Append' or 'Prepend' command, or they are system assigned names like "
	", e.g., 'StereoCompositingShader' for builtin shaders that are needed for stereo processing."
	"\n\n"
	"Return arguments, all optional: 'slot' is the index of the named slot in the chain to which the queried subfunction is assigned, or "
	"-1 if no such slot exists in the chain. 'glslid' numeric GLSL handle of the OpenGL GLSL shader object if the slot contains a "
	"GLSL shader for image processing on the GPU, 0 otherwise. 'luttexid' OpenGL texture handle of the first assigned lookup texture, "
	"0 if none assigned. 'voidptr' Memory pointer (encoded as double) to a C callback function, if one is assigned to this slot, 0 "
	"otherwise. 'blittercfg' either a parameter string with a meaning dependent of slot type, or the string 'NONE' if none assigned. "
	"'idstring' The symbolic name of this slot if any assigned, 'NONE' otherwise. 'insertString' the subcommand to use for reinserting "
	"this slot at the place it was after a deletion."
	"\n\n"
	"Screen('HookFunction', windowPtr, 'AppendShader', hookname, idString, glslid [, blittercfg] [luttexid1]); \n"
	"Append a new instruction slot to the end of hook processing chain 'hookname', assign the symbolic name 'idString' for later query "
	"by the 'Query' command. 'glslid' must be the name of a valid GLSL program object which defines the algorithm to apply on the GPU. "
	"'blittercfg' optional string with configuration commands for the shader. 'luttexid' Optional handle to an OpenGL texture which is "
	"used to encode lookup tables for the shader."
	"\n\n"
	"Screen('HookFunction', windowPtr, 'PrependShader', hookname idString, glslid [, blittercfg] [luttexid1]); \n"
	"Same as 'AppendShader' but add shader slot to beginning of the hook chain. It's recommended that you prepend slots instead of "
	"appending them, because PTB itself may add special slots at the end of a chain."
	"\n\n"
	"Screen('HookFunction', windowPtr, 'AppendCFunction', hookname, idString, voidfunctionptr); \n"
	"Screen('HookFunction', windowPtr, 'PrependCFunction', hookname, idString, voidfunctionptr); \n"
	"Attach a C callable function to the chain. voidfunctionptr is a double value which encodes a memory pointer to the function in "
	"memory. Encoding of the pointer and requirements to the function are non-trivial, so this is for expert developers only!"
	"\n\n"
	"Screen('HookFunction', windowPtr, 'AppendMFunction', hookname, idString, fevalstring); \n"
	"Screen('HookFunction', windowPtr, 'PrependMFunction', hookname, idString, fevalstring); \n"
	"Add a Matlab callable function to the chain. 'fevalstring' is the function call string: It will be passed to "
	"and evaluated by Matlabs or Octaves feval() function, so it has to work with that function. It is not allowed to "
	"return any return arguments. Some special names in the string will be replaced by PTB with internal settings, think "
	"of it as a macro replacement. Caution: The called function is not allowed to call any Screen() commands or you'll "
	"likely get undefined behaviour or a crash -- Screen is not reentrant!"
	"\n\n"
	"Screen('HookFunction', windowPtr, 'AppendBuiltin', hookname, idString, builtincmd); \n"
	"Screen('HookFunction', windowPtr, 'PrependBuiltin', hookname, idString, builtincmd); \n"
	"Add a call to a PTB built-in function 'builtincmd'."
	"\n\n"
	"You can also insert a hook function somewhere at a specific slot index via:\n"
	"Screen('HookFunction', windowPtr, 'InsertAtXXXYYY', ...);\n"
	"where XXX is a numeric slot id and YYY is the type specifier. Example: Insert a Builtin function\n"
	"at index 4 in hook 'UserDefinedBlit':\n"
	"Screen('HookFunction', windowPtr, 'InsertAt4Builtin', 'UserDefinedBlit', ...);\n"
	"\n\n"
	"Screen('Hookfunction', windowPtr, 'Remove', hookname, slotindex);\n"
	"Remove slot at index 'slotindex' in hookchain 'hookname'. The slot after this slot will move up by one.\n"
	"\n\n"
	"Screen('HookFunction', windowPtr, 'Enable', hookname); \n"
	"Screen('HookFunction', windowPtr, 'Disable', hookname); \n"
	"Enable or disable a specific hook chain. Chains are disabled until you enable them, with the exception of a few "
	"internal chains that get initialized and enabled by PTB itself, e.g., stereo algorithm chain. Disabled chains "
	"are not processed."
	"\n\n"
	"Screen('HookFunction', windowPtr, 'Reset', hookname); \n"
	"Reset a processing hook chain: All slots are deleted, resetting the chain to its startup state. Seldomly needed. "
	"\n\n"
	"Screen('HookFunction', windowPtr, 'Dump', hookname); \n"
	"Print out the full chain for hook 'hookname' to the Matlab console in a human readable format - Useful for debugging."
	"\n\n"
	"Screen('HookFunction', windowPtr, 'DumpAll'); \n"
	"Print out all chains for the given onscreen window 'windowPtr' to the Matlab console in a human readable format - Useful for debugging."
	"\n\n"
	"oldImagingMode = Screen('HookFunction', proxyPtr, 'ImagingMode' [, imagingMode]); \n"
	"Change or query imagingMode flags of provided proxy window 'proxyPtr' to 'imagingMode'. Proxy windows are used to define "
	"image processing operations, mostly for Screen('TransformTexture'). Returns old imaging mode."
	"\n\n"
	"General notes:\n\n"
	"* Hook chains are per onscreen window, so each window can have a different configuration and enable state.\n"
	"* Read all available documentation on the Psychtoolbox imaging pipeline in 'help PsychGLImageprocessing', the PsychDocumentation folder "
	"and on the Wiki before you make use of this function. Its way to complex to use it by guessing ;)\n";
	
static char seeAlsoString[] = "";

PsychError SCREENHookFunction(void) 
{
	PsychWindowRecordType	*windowRecord;
	char					numString[10];
	char					*cmdString, *hookString, *idString, *blitterString, *insertString;
	int						cmd, slotid, whereloc = 0;
	double					doubleptr;
	double					shaderid, luttexid1 = 0;

	blitterString = NULL;

    // All subfunctions should have these two lines.  
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
    
    PsychErrorExit(PsychCapNumInputArgs(7));   	
    PsychErrorExit(PsychRequireNumInputArgs(2)); 	
    PsychErrorExit(PsychCapNumOutputArgs(7));  

    // Get the subcommand string:
	PsychAllocInCharArg(2, kPsychArgRequired, &cmdString);
	
	// Subcommand dispatcher:
	cmd=0;
	if (strstr(cmdString, "Append"))  { cmd=1; whereloc = INT_MAX; }
	if (strstr(cmdString, "Prepend")) { cmd=2; whereloc = 0; }	
	if (strcmp(cmdString, "Reset")==0)   cmd=3;
	if (strcmp(cmdString, "Enable")==0)  cmd=4;
	if (strcmp(cmdString, "Disable")==0) cmd=5;
	if (strcmp(cmdString, "Query")==0)   cmd=6;
	if (strcmp(cmdString, "Dump")==0)    cmd=7;
	if (strcmp(cmdString, "DumpAll")==0) cmd=8;
	if (strcmp(cmdString, "ListAll")==0) cmd=9;
	if (strcmp(cmdString, "Edit")==0)   cmd=10;
	if (strcmp(cmdString, "ImagingMode")==0) cmd=11;
	if (strstr(cmdString, "InsertAt")) { cmd=12; whereloc = -1; sscanf(cmdString, "InsertAt%i", &whereloc); }
	if (strstr(cmdString, "Remove")) cmd=13;
	
	if(cmd==0) PsychErrorExitMsg(PsychError_user, "Unknown subcommand specified to 'HookFunction'.");
	if(whereloc < 0) PsychErrorExitMsg(PsychError_user, "Unknown/Invalid/Unparseable insert location specified to 'HookFunction' 'InsertAtXXX'.");
	
	// Need hook name?
	if(cmd!=9 && cmd!=8 && cmd!=11) {
		// Get it:
		PsychAllocInCharArg(3, kPsychArgRequired, &hookString);
	}
	
    // Get the window structure for the onscreen window.
	windowRecord = NULL;
    PsychAllocInWindowRecordArg(1, (cmd!=9) ? TRUE : FALSE, &windowRecord);
    
	switch(cmd) {
		case 1:  // Append:
		case 2:  // Prepend:
		case 12: // Insert at location 'whereloc':
			// Add a new hook function callback to chain, either at beginning or end.
			
			// What type of callback/handler is to be added?
			if(strstr(cmdString, "Shader")) {
				// GLSL shader program object:
				
				// The id string:
				PsychAllocInCharArg(4, kPsychArgRequired, &idString);
				
				// The shader object handle:
				PsychCopyInDoubleArg(5, TRUE, &shaderid);

				// The blitter config string:
				PsychAllocInCharArg(6, FALSE, &blitterString);

				// (Optionally) a texture handle for a texture to attach to 2nd unit:
				PsychCopyInDoubleArg(7, FALSE, &luttexid1);
				
				// Add shader: 
				PsychPipelineAddShaderToHook(windowRecord, hookString, idString, whereloc, (unsigned int) shaderid, blitterString, (unsigned int) luttexid1);
			}
			else if(strstr(cmdString, "CFunction")) {
				// C callback function:

				// First the id string:
				PsychAllocInCharArg(4, kPsychArgRequired, &idString);
				
				// Then the void* to the function, encoded as a double value:
				PsychCopyInDoubleArg(5, TRUE, &doubleptr);
				
				// Add the function void* to the chain:
				PsychPipelineAddCFunctionToHook(windowRecord, hookString, idString, whereloc, PsychDoubleToPtr(doubleptr));
				
			}
			else if(strstr(cmdString, "MFunction")) {
				// Matlab/Octave/Whatever runtime environment callback function:

				// First the id string:
				PsychAllocInCharArg(4, kPsychArgRequired, &idString);
				
				// Then the call string (for feval()) for the function:
				PsychAllocInCharArg(5, TRUE, &blitterString);
				
				// Add the function to the chain:
				PsychPipelineAddRuntimeFunctionToHook(windowRecord, hookString, idString, whereloc, blitterString);
				
			}
			else if(strstr(cmdString, "Builtin")) {
				// Built in special function:

				// First the id string:
				PsychAllocInCharArg(4, kPsychArgRequired, &idString);
				
				// Then the options-string:
				PsychAllocInCharArg(5, TRUE, &blitterString);
				
				// Add the function to the chain:
				PsychPipelineAddBuiltinFunctionToHook(windowRecord, hookString, idString, whereloc, blitterString);				
			}
			else {
				// Unknown?!?
				PsychErrorExitMsg(PsychError_user, "Unknown callback type specified to 'HookFunction'.");
			}
		break;
		
		case 3: // Reset hook-chain:
			PsychPipelineResetHook(windowRecord, hookString);
		break;

		case 4: // Enable hook-chain:
			PsychPipelineEnableHook(windowRecord, hookString);
		break;

		case 5: // Disable hook-chain:
			PsychPipelineDisableHook(windowRecord, hookString);
		break;
		
		case 6: // Query properties of a slot in a specific hook-chain:
			// Get the id string:
			if (PsychGetArgType(4)!=PsychArgType_char) {
				// No id string provided: Try to get numeric slot id and
				// turn it into an id string:
				PsychCopyInIntegerArg(4, TRUE, &slotid);
				sprintf(numString, "%i", slotid);
				idString = (char*) &numString[0];
				slotid = -1;
			}
			else {
				PsychAllocInCharArg(4, TRUE, &idString);
			}
			
			// Query everything that's there and copy it out:
			slotid = PsychPipelineQueryHookSlot(windowRecord, hookString, &insertString, &idString, &blitterString, &doubleptr, &shaderid, &luttexid1);
			
			// Copy out all infos:
			PsychCopyOutDoubleArg(1, FALSE, slotid);
			PsychCopyOutCharArg(2, FALSE, (idString) ? idString : "NONE");
			PsychCopyOutCharArg(3, FALSE, (blitterString) ? blitterString : "NONE");
			PsychCopyOutDoubleArg(4, FALSE, doubleptr);
			PsychCopyOutDoubleArg(5, FALSE, shaderid);
			PsychCopyOutDoubleArg(6, FALSE, luttexid1);
			PsychCopyOutCharArg(7, FALSE, (insertString) ? insertString : "NONE");
		break;

		case 7: // Dump specific hook-chain:
			PsychPipelineDumpHook(windowRecord, hookString);
		break;

		case 8: // Dump all hook-chains:
			PsychPipelineDumpAllHooks(windowRecord);
		break;
		
		case 9: // List all hook-chains:
			PsychPipelineListAllHooks(windowRecord);
		break;
		
		case 10: // Set properties of a slot in a specific hook-chain:
			// Get the id string:
			PsychAllocInCharArg(4, kPsychArgRequired, &idString);
			
			// Query everything that's there and copy it out:
			slotid = PsychPipelineQueryHookSlot(windowRecord, hookString, &insertString, &idString, &blitterString, &doubleptr, &shaderid, &luttexid1);
			
			if (slotid<0)  PsychErrorExitMsg(PsychError_user, "In 'Edit' No such hook slot in that hook chain for that object.");
			// Copy out all infos:
			PsychCopyOutDoubleArg(1, FALSE, slotid);
			PsychCopyOutCharArg(2, FALSE, (idString) ? idString : "NONE");
			PsychCopyOutCharArg(3, FALSE, (blitterString) ? blitterString : "NONE");
			PsychCopyOutDoubleArg(4, FALSE, doubleptr);
			PsychCopyOutDoubleArg(5, FALSE, shaderid);
			PsychCopyOutDoubleArg(6, FALSE, luttexid1);
		break;

		case 11: // Change imagingMode flags for a proxy window:
			if (windowRecord->windowType!=kPsychProxyWindow) PsychErrorExitMsg(PsychError_user, "In 'ImagingMode' Invalid windowPtr provided. Must be a proxy window!");
			PsychCopyOutDoubleArg(1, FALSE, (double) windowRecord->imagingMode);
			
			luttexid1 = -1;
			PsychCopyInDoubleArg(3, FALSE, &luttexid1);
			if (luttexid1 < 0 && luttexid1!=-1) PsychErrorExitMsg(PsychError_user, "In 'ImagingMode' Invalid imagingMode flags provided. Must be positive!");
			if (luttexid1!=-1) {
				windowRecord->imagingMode = (unsigned int) luttexid1;
			}
		break;
		
		// case 12 see at top.
		
		case 13: // Remove slot at given index.
			PsychCopyInIntegerArg(4, TRUE, &slotid);
			PsychPipelineDeleteHookSlot(windowRecord, hookString, slotid);
		break;
	}
	
    // Done.
    return(PsychError_none);
}
