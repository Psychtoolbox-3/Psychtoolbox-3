/*

	/osxptb/trunk/PsychSourceGL/Source/Common/Eyelink/PsychEyelink.c

	PROJECTS: Eyelink 
  
	AUTHORS:

		cburns@berkeley.edu				cdb
		E.Peters@ai.rug.nl				emp
		f.w.cornelissen@med.rug.nl		fwc
		mario.kleiner@tuebingen.mpg.de  mk
  
	PLATFORMS:	All.
    
	HISTORY:

		11/22/05  cdb		Created.
		29/06/06  fwc		fixed EyelinkSystemIsConnected to allow dummy mode connections
		15/03/09  mk		Added experimental support for eye camera image display.

	TARGET LOCATION:

		Eyelink.mexmac resides in:
			EyelinkToolbox
*/

#include "PsychEyelink.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////
// Global variables used throughout eyelink C files

int		giSystemInitialized = 0;
int		verbosity = 2;

// Callback string for eyelink display callback function:
char eyelinkDisplayCallbackFunc[1024];

// Memory pointer to malloc()'ed image pixel buffer that holds the
// image data for a RGBA8 texture with the most recent eye camera image:
static byte* eyeimage = NULL;

// Width x Height of eye camera image in pixels:
static int eyewidth  = 0;
static int eyeheight = 0;

// Color remapping palette table:
static unsigned int palmap32[256];
#define ERR_BUFF_LEN 1000

/* Declaration of callback functions defined later in this file: */
static INT16 ELCALLBACK  PsychEyelink_setup_image_display(INT16 width, INT16 height);
static void ELCALLBACK   PsychEyelink_exit_image_display(void);
static void ELCALLBACK   PsychEyelink_set_image_palette(INT16 ncolors, byte r[130], byte g[130], byte b[130]);
static void ELCALLBACK   PsychEyelink_draw_image_line(INT16 width, INT16 line, INT16 totlines, byte *pixels);

static INT16  ELCALLBACK PsychEyelink_setup_cal_display(void);
static void ELCALLBACK   PsychEyelink_exit_cal_display(void);
static void ELCALLBACK   PsychEyelink_clear_display(void);
static void ELCALLBACK   PsychEyelink_draw_cal_target(INT16 x, INT16 y);
static void ELCALLBACK   PsychEyelink_erase_cal_target(void);
static void ELCALLBACK   PsychEyelink_image_title(INT16 threshold, char *title);
static INT16 ELCALLBACK  PsychEyelink_get_input_key(InputEvent *keyinput);
static void ELCALLBACK   PsychEyelink_alert_printf_hook(const char *msg);
static void ELCALLBACK	 PsychEyelink_noop(void);

static void ELCALLBACK   PsychEyelink_cal_target_beep_hook(void);
static void ELCALLBACK	 PsychEyelink_cal_done_beep_hook(INT16 error);
static void ELCALLBACK	 PsychEyelink_dc_done_beep_hook(INT16 error);
static void ELCALLBACK   PsychEyelink_dc_target_beep_hook(void);

/////////////////////////////////////////////////////////////////////////
// Check if system is initialized
//
PsychError EyelinkSystemIsConnected(void)
{
	int iStatus=-9999;
	iStatus=eyelink_is_connected();
//	mexPrintf("EyelinkSystemIsConnected status %d ((iStatus==0)=%d)\n", iStatus, (iStatus==0) );
	if (iStatus==0) {
		PsychErrorExitMsg(PsychError_user, "Eyelink system is not connected!\n");
	}
/*	
	if (eyelink_is_connected()==0) {
		PsychErrorExitMsg(PsychError_user, "Eyelink system is not connected!\n");
	}
	*/
	return(PsychError_none);
}

/////////////////////////////////////////////////////////////////////////
// Check is system is initialized
//
PsychError EyelinkSystemIsInitialized(void)
{
	if (giSystemInitialized != 1) {
		PsychErrorExitMsg(PsychError_user, "Eyelink system is not initialized!\n");
	}
	return(PsychError_none);
}

/* Eyelink('Verbosity') - Set level of verbosity.
 */
PsychError EyelinkVerbosity(void) 
{
 	static char useString[] = "oldlevel = Eyelink('Verbosity' [,level]);";
	static char synopsisString[] = 
		"Set level of verbosity for error/warning/status messages. 'level' optional, new level "
		"of verbosity. 'oldlevel' is the old level of verbosity. The following levels are "
		"supported: 0 = Shut up. 1 = Print errors, 2 = Print also warnings, 3 = Print also some info, "
		"4 = Print more useful info (default), >5 = Be very verbose (mostly for debugging the driver itself). ";		
	static char seeAlsoString[] = " ";	 

	int level= -1;

	// Setup online help: 
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
	
	PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(1));	 // The maximum number of outputs

	PsychCopyInIntegerArg(1, kPsychArgOptional, &level);
	if (level < -1) PsychErrorExitMsg(PsychError_user, "Invalid level of verbosity provided. Valid are levels of zero and greater.");
	
	// Return current/old level:
	PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) verbosity);

	// Set new level, if one was provided:
	if (level > -1) verbosity = level;

	return(PsychError_none);
}

// Return level of verbosity:
int Verbosity(void) {
	return(verbosity);
}

// Parse printf() style format string and variable number of
// integer or string arguments into a printf() formatted
// string and return static pointer to the final string.
// Used, e.g., by Eyelink('Command') and Eyelink('Message'):
const char* PsychEyelinkParseToString(int startIdx)
{
	static char			strCommand[256];
	int				i			= 0;
	int				iNumInArgs		= 0;
	PsychArgFormatType	        psychArgType	        = PsychArgType_none;
	int                             iTempValue              = 0;
	char				*pstrTemp		= NULL;
	char				*pstrFormat		= NULL;
	char                            strFragment[256];
	char                            fSpec[256];
	int                             wIdx = 0;
	int                             argIdx;

	// Alloc and grab the input format string
	PsychAllocInCharArg(startIdx, TRUE, &pstrFormat);
	iNumInArgs = PsychGetNumInputArgs();   

	// Define start index of variable argument list:
	argIdx = startIdx + 1;

	// Clear strings
	memset(strCommand, 0, sizeof(strCommand));

	// Parse complete format string:
	while ((*pstrFormat != 0) && (wIdx < 255)) {
	  // Special character % detected?
	  if ((*pstrFormat != '%') || (*(pstrFormat+1) == '%')) {
	    // Easy: Regular char or escaped %. Just copy into target command string:

	    // Eat up the escape '%' character, if any:
	    if (pstrFormat == strstr(pstrFormat, "%%")) pstrFormat++;

	    // Copy escaped single % or regular character:
	    strCommand[wIdx++] = *(pstrFormat++);

	    // Next character...
	    continue;
	  }

	  // Special % char detected, which is not escaped, therefore
	  // a datatype format specifier follows immediately:

	  // Is there an argument available to match the format string spec?
	  if (iNumInArgs < argIdx) {
	    PsychErrorExitMsg(PsychError_user, "Number of supplied arguments does not match number of arguments required by format string!");
	  }

	  // Find end of actual parameter spec:
	  for (i = 0; (pstrFormat[i] > 0) && (pstrFormat[i] != ' '); i++);

	  // Copy format substring to fSpec:
	  memset(fSpec, 0, sizeof(fSpec));
	  strncpy(fSpec, pstrFormat, (i < 256) ? i : 255);

	  // Prepare output substring for writing:
	  memset(strFragment, 0, sizeof(strFragment));

	  // Check if input argument type matches parameter spec string
	  // and assign, if so, abort otherwise:
	  psychArgType = PsychGetArgType(argIdx);
	  switch(psychArgType) {
	    case PsychArgType_double:
	      if ((PsychGetArgM(argIdx) == 1) && (PsychGetArgN(argIdx) == 1)) {
		PsychCopyInIntegerArg(argIdx, TRUE, &iTempValue);
		
		// Got a int value. Was a int value expected?
		if (strstr(fSpec, "d") || strstr(fSpec, "i")) {
		  // Yes: Print into output string fragment:
		  snprintf(strFragment, 255, fSpec, iTempValue);
		} else {
		  // No: This is a mismatch - Game over:
		  PsychErrorExitMsg(PsychError_user, "Mismatch between provided scalar integer argument and expected argument!");
		}
	      } else {
		PsychGiveHelp();
		PsychErrorExitMsg(PsychError_user, "");
	      }
	      break;

	    case PsychArgType_char:
	      PsychAllocInCharArg(argIdx, TRUE, &pstrTemp);
	      // Got a string. Was a string expected?
	      if (strstr(fSpec, "s")) {
		// Yes: Print into output string fragment:
		snprintf(strFragment, 255, fSpec, pstrTemp);
	      } else {
		// No: This is a mismatch - Game over:
		PsychErrorExitMsg(PsychError_user, "Mismatch between provided character string and expected argument!");
	      }
	      break;

	    default:
	      PsychGiveHelp();
	      PsychErrorExitMsg(PsychError_user, "");
	      break;
	  }
	  
	  // If we made it here, then the strFragment is ready for
	  // joining:
	  if ((strlen(strCommand) + strlen(strFragment)) < 256) {
	    strcat(strCommand, strFragment);
	  } else {
	    // Break out of parser - Need to truncate:
	    break;
	  }

	  // Advance parse positions:
	  wIdx = strlen(strCommand);
	  pstrFormat += i;
	  argIdx++;

	  // Next parse iteration.
	}

	// Sanity check:
	if (*pstrFormat != 0) printf("Eyelink-Warning:Final overall command truncated to '%s'!\nMaximum of 255 characters allowed.\n", strCommand);

	// Return pointer to internally statically allocated final character string:
	return(strCommand);
}

// Initialize all callback hook functions for use by Eyelink runtime, e.g.,
// all the callbacks for eye camera image display:
void PsychEyelink_init_core_graphics(const char* callback)
{
	HOOKFCNS fcns;
	
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_init_core_graphics()\n");
	
	memset(&fcns, 0, sizeof(fcns));
	
	// Setup cam image callbacks:
	fcns.setup_image_display_hook = PsychEyelink_setup_image_display;
	fcns.exit_image_display_hook= PsychEyelink_exit_image_display;
	fcns.set_image_palette_hook = PsychEyelink_set_image_palette;
	fcns.image_title_hook       = PsychEyelink_image_title;
	fcns.draw_image_line_hook   = PsychEyelink_draw_image_line;

	// Setup calibration callbacks:
	fcns.setup_cal_display_hook = PsychEyelink_setup_cal_display;
	fcns.exit_cal_display_hook  = PsychEyelink_exit_cal_display;
	fcns.clear_cal_display_hook = PsychEyelink_clear_display;
	fcns.draw_cal_target_hook   = PsychEyelink_draw_cal_target;
	fcns.erase_cal_target_hook  = PsychEyelink_erase_cal_target;

	// Setup keyboard and I/O callbacks:
	fcns.get_input_key_hook     = PsychEyelink_get_input_key;
	fcns.alert_printf_hook      = PsychEyelink_alert_printf_hook;

	// Set auditory feedback callbacks:
	fcns.cal_target_beep_hook   = PsychEyelink_cal_target_beep_hook;
	fcns.cal_done_beep_hook		= PsychEyelink_cal_done_beep_hook;
	fcns.dc_target_beep_hook	= PsychEyelink_dc_target_beep_hook;
	fcns.dc_done_beep_hook		= PsychEyelink_dc_done_beep_hook;

	// Just set this to make eyelink-core happy:
	fcns.record_abort_hide_hook = PsychEyelink_noop;

	// Assign runtime environment display callback function:
	memset(eyelinkDisplayCallbackFunc, 0, sizeof(eyelinkDisplayCallbackFunc));

	snprintf(eyelinkDisplayCallbackFunc, sizeof(eyelinkDisplayCallbackFunc) - 1, "%s", callback);
//    #if PSYCH_SYSTEM != PSYCH_WINDOWS
//	snprintf(eyelinkDisplayCallbackFunc, sizeof(eyelinkDisplayCallbackFunc) - 1, "%s", callback);
//    #else
//	_snprintf(eyelinkDisplayCallbackFunc, sizeof(eyelinkDisplayCallbackFunc) - 1, "%s", callback);
//    #endif
    
	// Assign hooks to Eyelink runtime:
	setup_graphic_hook_functions(&fcns);
	
	// Optionally dump the whole hookfunctions struct:
	if (Verbosity() > 5) PsychEyelink_dumpHookfunctions();
	
	return;
}

// Disable all hook functions at shutdown time:
void PsychEyelink_uninit_core_graphics(void)
{
	HOOKFCNS fcns;
	
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_uninit_core_graphics()\n");
	
	memset(&fcns, 0, sizeof(fcns));
	setup_graphic_hook_functions(&fcns);

	// Optionally dump the whole hookfunctions struct:
	if (Verbosity() > 5) PsychEyelink_dumpHookfunctions();

	return;
}

void PsychEyelink_dumpHookfunctions(void)
{
	HOOKFCNS* pfcns = get_all_hook_functions();
	int i;
	
	printf("PsychEyelink: Dump of current Eyelink HOOKFCNS struct as byte array:\n\n");
	for (i=0; i < sizeof(HOOKFCNS); i++) printf(" %02x", (int)(((unsigned char*) pfcns)[i]));
	printf("\nPsychEyelink: Dump done\n\n");

	return;
}

void PsychEyelink_TestEyeImage(void)
{
	int i, x, y;
	byte r[256], g[256], b[256];
	byte scanline[640];
	InputEvent keyinput;
	
	// Pseudo-Eyelink camera image test:
	
	// Setup pseudo eye-display of 640 x 480 pixels via setup callback:
	PsychEyelink_setup_image_display(640, 480);
	
	// Build pseudo color LUT:
	for (i=0; i < 256; i++) {
		r[i]=i;
		g[i]=255 - i;
		b[i]=i * 2;
	}
	PsychEyelink_set_image_palette(256, r, g, b);

	// Set image title:
	PsychEyelink_image_title(1, "Foobar-O-Matic:");

	// Calibration beep:
	PsychEyelink_cal_target_beep_hook();
	
	// Run pseudo-display loop for 600 frames:
	for (i = 0; i < 600; i++) {
		// Draw calibration target:
		PsychEyelink_draw_cal_target(i, 200);

		// Fill buffer with image pattern:
		for (y=1; y <= 480; y++) {
			// Build y'th scanline:
			for (x=0; x < 640; x++) scanline[x] = (byte) ((x + y + i) % 256);
			
			// Submit y'th scanline:
			PsychEyelink_draw_image_line(640, y, 480, (byte*) &scanline);
		}
		
		// Check keyboard:
		keyinput.key.key = 0;
		if (PsychEyelink_get_input_key(&keyinput) > 0) {
			PsychEyelink_alert_printf_hook("Eyelink: Key detected.\n");
			// Break out of loop on keycode 41 or 27 == ESCAPE on OS/X or Windows.
			if (keyinput.key.key == 41 || keyinput.key.key == 27) break;
			if (keyinput.key.key == TERMINATE_KEY) {
				printf("Eyelink: TestSuite: WARNING: Abort code detected. Master abort.\n");
				break;
			}
		}
	}

	// Calibration end beep hook:
	PsychEyelink_cal_done_beep_hook(keyinput.key.key);
	
	// Tear down pseudo display:
	PsychEyelink_exit_image_display();

	// Test calibration display:
	PsychEyelink_setup_cal_display();
	PsychEyelink_dc_target_beep_hook();
	PsychEyelink_clear_display();
	PsychEyelink_dc_done_beep_hook(keyinput.key.key);
	
	PsychEyelink_exit_cal_display();

	return;
}

int PsychEyelinkCallRuntime(int cmd, int x, int y, char* msg)
{
	PsychGenericScriptType	*inputs[2];
	PsychGenericScriptType	*outputs[1];
	double* callargs;
	double rc;
	int retc;

	// Callbacks forcefully disabled by error-handling? Return with error code if so:
	if (0 == eyelinkDisplayCallbackFunc[0]) return(0xdeadbeef);

	// Create a Matlab double matrix with 4 elements: 1st is command code 
	// others are available for use specific to each command
	outputs[0]  = NULL;
	inputs[0]   = mxCreateDoubleMatrix(1, 4, mxREAL);
	callargs    = mxGetPr(inputs[0]);
	
	callargs[0] = (double) cmd; // Command code.
	callargs[1] = (double) x;
	callargs[2] = (double) y;

	if (msg != NULL) {
		inputs[1] = mxCreateString(msg);
	}
	else {
		inputs[1] = NULL;
	}
	
	// Call the runtime environment, on Matlab with the trap flag set, so control
	// returns to us on error, instead of returning to Matlab runtime system. Eyelink
	// runtime system doesn't like losing control and would crash otherwise!
	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		mexSetTrapFlag(1);
	#endif
	
	// Call the runtime environment:
	if ((retc = mexCallMATLAB((cmd == 2) ? 1 : 0, outputs, (inputs[1]) ? 2 : 1, inputs, eyelinkDisplayCallbackFunc)) > 0) {
		printf("EYELINK: WARNING! PsychEyelinkCallRuntime() Failed to call eyelink runtime callback function %s [rc = %i]!\n", eyelinkDisplayCallbackFunc, retc);
		printf("EYELINK: WARNING! Make sure that function is on your Matlab/Octave path and properly initialized.\n");
		printf("EYELINK: WARNING! May also be an error during execution of that function. Type ple at command prompt for error messages.\n");
		printf("EYELINK: WARNING! Auto-Disabling all callbacks to the runtime environment for safety reasons.\n");
		eyelinkDisplayCallbackFunc[0] = 0;
	}

	// Reset error handling to default on Matlab:
	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		mexSetTrapFlag(0);
	#endif

	// Release our matrix again:
	mxDestroyArray(inputs[0]);
	if (msg != NULL) mxDestroyArray(inputs[1]);
	
	if (outputs[0]) {
		rc = mxGetScalar(outputs[0]);
		mxDestroyArray(outputs[0]);
	}
	else {
		rc = 0;
	}
	
	return((int) rc);
}

// Callback functions, called by Eyelink runtime at various occassions, e.g,
// during tracker setup, drift correction/calibration etc.:
// =========================================================================

static void ELCALLBACK PsychEyelink_noop(void)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_noop()\n");
	
	// Done.
	return;
}

// PsychEyelink_setup_image_display() tells the width and height of the camera
// image in pixels.
static INT16 ELCALLBACK PsychEyelink_setup_image_display(INT16 width, INT16 height)
{
	
	
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_setup_image_display()\n");
	// Release any stale image buffer:
	if (eyeimage != NULL) free(eyeimage);
	
	// Reset everything to startup default:
	eyeimage  = NULL;
	eyewidth  = 0;
	eyeheight = 0;
	
	if (width < 1 || height < 1) {
		printf("EYELINK: WARNING! Invalid image dimensions (smaller than 1 pixel!) received from eyelink: Aborting image setup.\n");
		return(-1);
	}

	// Allocate an internal memory buffer of sufficient size to hold an image
	// of size width x height pixels:
	eyeimage = (byte*) malloc(sizeof(unsigned char) * 4 * width * height);
	if (eyeimage != NULL) {
		eyewidth  = width;
		eyeheight = height;
	}
	else {
		// Failed:
		return(-1);
	}
	
	// Tell callback about image dimensions fwiw:
	if (0xdeadbeef == PsychEyelinkCallRuntime(8, eyewidth, eyeheight, NULL)) {
		// Error condition. Return error to eyelink runtime:
		return(-1);
	}
	
	if (Verbosity() > 5) printf("Eyelink: Leaving PsychEyelink_setup_image_display()\n");

	// Done.
	return(0);
}

// PsychEyelink_exit_image_display() shuts down any camera image display:
static void ELCALLBACK PsychEyelink_exit_image_display(void)
{
	
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_exit_image_display()\n");

	// Release any allocated image buffer:
	if (eyeimage != NULL) free(eyeimage);
	
	// Reset everything to startup default:
	eyeimage  = NULL;
	eyewidth  = 0;
	eyeheight = 0;
	
	// Tell runtime to exit display: Command code 9.
	PsychEyelinkCallRuntime(9, 0, 0, NULL);
	
	// Done.
	return;
}
// added by NJ @ SR Research Sept 2010
#define UPSIDE 0
#define LEFTSIDE 1
#define RIGHTSIDE 2
#define DOWNSIDE 3
void drawSemiCircle(CrossHairInfo *chi, int left, int top, int dia, int side, int cindex)
{
	
	

	unsigned char r =0;
	unsigned char g =0;
	unsigned char b =0;
	int radius = dia/2;
	int x = left - 1;
	int y = top -1;
	unsigned int *v0;
	int x0,y0, ddF_x =1, ddF_y,f;
	

	if (eyeimage == NULL){
		return;
	}

	switch(cindex)
	{
		case CR_HAIR_COLOR:          r=g=b =255; break;//255,255,255
		case PUPIL_HAIR_COLOR:       r=g=b =255; break;//255,255,255
		case PUPIL_BOX_COLOR:			 g =255; break;//0,255,0
		case SEARCH_LIMIT_BOX_COLOR: 
		case MOUSE_CURSOR_COLOR:		r = 255; break;//255,0,0
	}
	
	v0 = (unsigned int*) (eyeimage);
	
	// implement manual clipping to mimic behavior on host pc
	
	if(side == UPSIDE)
	{
		
		x0 = left;
		y0 = top;
		radius = dia / 2;	
		y= radius;
		f = 1 - radius;
		ddF_y = -2 * radius;
		x =  0;
		y = radius;
		x0 = x0 + dia/2;
		y0 = y0 + dia/2;

		
		if (y0 < eyeheight && y0 > 0 && y0-radius < eyeheight && y0-radius > 0 && y0-radius < eyeheight &&  x0+radius < eyewidth && x0-radius >0 && x0+radius > 0 && x0-radius < eyewidth){
			v0[(eyewidth*eyeheight) - (y0-radius)*eyewidth+x0] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			v0[(eyewidth*eyeheight) - y0*eyewidth+(x0+radius)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			v0[(eyewidth*eyeheight) - y0*eyewidth+(x0-radius)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		}
		
		while( x < y ){	 
			
			if(f >= 0){
				y--;
				ddF_y += 2;
				f += ddF_y;
				
			}
			x++;
			ddF_x += 2;
			f += ddF_x;
			

			if (y0+y >0 && y0+y < y0 && y0+y < eyeheight) {
				
				if(x0+x < eyewidth)
				v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if (x0-x > 0 )
				v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			if (y0-y >0 && y0-y < y0 && y0-y < eyeheight){//if (y0-y < y0 && y0-y < eyeheight){
				
				if (x0+x < eyewidth)
				v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if (x0-x > 0 )
				v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			if (y0+x >0 && y0+x < y0 && y0+x < eyeheight){//if (y0+x < y0 && y0+x < eyeheight){
				
				if (x0+y < eyewidth)
				v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if((x0-y) >0)
				v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			if (y0-x > 0 && y0-x < y0 && y0-x < eyeheight){//if (y0-x < y0 && y0-x < eyeheight){
			
				if(x0+y < eyewidth)
				v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if(x0-y > 0)
				v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
		}
	}else if (side == DOWNSIDE){
		
		radius = dia / 2;
		y= radius;
		f = 1 - radius;
		ddF_y = -2 * radius;
		ddF_x =1;
		x =  0;
		y = radius;
		
		x0 = left;
		y0 = top;
		x0 = x0 + dia/2;
		
		if (y0+radius < eyeheight && y0+radius > 0 && x0 > 0 && x0 < eyewidth)
			v0[(eyewidth*eyeheight) - (y0+radius)*eyewidth+x0] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		if (y0 < eyeheight && y0 > 0 && x0+radius > 0 && x0+radius < eyewidth)
			v0[(eyewidth*eyeheight) - y0*eyewidth+(x0+radius)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		if (y0 < eyeheight && y0 > 0 && x0-radius > 0 && x0-radius < eyewidth)
			v0[(eyewidth*eyeheight) - y0*eyewidth+(x0-radius)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		
		
		while( x < y ){	 
			
			if(f >= 0){
				y--;
				ddF_y += 2;
				f += ddF_y;
				
			}
			x++;
			ddF_x += 2;
			f += ddF_x;
			
			if (y0+y >0 && y0+y > y0 && y0+y < eyeheight) {
				if ( x0+x < eyewidth)
					v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if ( x0-x > 0 )
					v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			if (y0-y > 0 && y0-y > y0 && y0-y < eyeheight){
				if ( x0+x < eyewidth)
					v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if ( x0-x > 0 )	
					v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			if (y0+x >0 && y0+x > y0 && y0+x < eyeheight){
				if ( x0+y < eyewidth )
					v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if ( x0-y > 0)
					v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			if (y0-x  >0 && y0-x > y0 && y0-x < eyeheight){
				if ( x0+y < eyewidth)
					v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if (x0-y > 0)
					v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
		}
		
		
	}else if (side == RIGHTSIDE){
		
		radius = dia / 2;
		y= radius;
		f = 1 - radius;
		ddF_y = -2 * radius;
		ddF_x =1;
		x =  0;
		y = radius;
		
		x0 = left;
		y0 = top;
		y0 = y0 + dia/2;
				
		
		if (x0 + radius < eyewidth && x0+radius > 0){
			v0[(eyewidth*eyeheight) - (y0+radius)*eyewidth+x0] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			v0[(eyewidth*eyeheight) - (y0+radius)*eyewidth+x0] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			v0[(eyewidth*eyeheight) - y0*eyewidth+(x0+radius)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		}
		
		while( x < y ){	 
			
			if(f >= 0){
				y--;
				ddF_y += 2;
				f += ddF_y;
				
			}
			x++;
			ddF_x += 2;
			f += ddF_x;
			
			if (x0+x  >0 && x0+x > x0 && x0+x < eyewidth) {
				if(y0+x < eyeheight)
					v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if(y0-x > 0)
					v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				
			}
			if (x0-x > 0 && x0-x > x0 && x0-x < eyewidth){
				if(y0+y < eyeheight)
					v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if(y0-y > 0)
					v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			if (x0+y > 0 && x0+y > x0 && x0+y < eyewidth){
				if(y0+x < eyeheight)
					v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if(y0-x > 0)
					v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);			
			}
			if (x0-y > 0 && x0-y > x0 && x0-y < eyewidth){
				if(y0+x < eyeheight)
					v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if(y0-x > 0)
					v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
		}
		
	}else if (side == LEFTSIDE){
		
		radius = dia / 2;
		y= radius;
		f = 1 - radius;
		ddF_y = -2 * radius;
		ddF_x =1;
		x =  0;
		y = radius;
		
		x0 = left;
		y0 = top;
		y0 = y0 + dia/2;
		x0 = x0 + dia/2; 

		if (x0 - radius > 0 && x0-radius < eyewidth){
		v0[(eyewidth*eyeheight) - (y0+radius)*eyewidth+x0] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		v0[(eyewidth*eyeheight) - (y0+radius)*eyewidth+x0] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		v0[(eyewidth*eyeheight) - y0*eyewidth+(x0-radius)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		}
		
		while( x < y ){	 
			
			if(f >= 0){
				y--;
				ddF_y += 2;
				f += ddF_y;
				
			}
			x++;
			ddF_x += 2;
			f += ddF_x;
			
			if (x0+x > 0 && x0+x < x0 && x0+x > 0) {
				v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);				
			}
			if (x0-x > 0 && x0-x < x0 && x0-x > 0){
				if(y0+x < eyeheight)
					v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if(y0-x > 0)	
					v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			if (x0+y > 0 && x0+y < x0 && x0+y > 0){
				if(y0+x < eyeheight)
					v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if(y0-x > 0)
					v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);			
			}
			if (x0-y > 0 && x0-y < x0 && x0-y > 0){
				if(y0+x < eyeheight)
					v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if(y0-x > 0)
					v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
		}
		
	}
	return;
}

// lack of graphics functions use primitive algorithm to draw circle by coloring in rgb in image memory 
void drawCircle(CrossHairInfo *chi, int x0, int y0, int width, int height, int cindex)
{
	
	

	unsigned char r =0;
	unsigned char g =0;
	unsigned char b =0;
	int x = 0, y, f;
	unsigned int *v0;
	int radius, ddF_x =1, ddF_y;
	
	if (eyeimage == NULL) return;

	switch(cindex){
		case CR_HAIR_COLOR:          r=g=b = 255; break;//255,255,255
		case PUPIL_HAIR_COLOR:       r=g=b = 255; break;//255,255,255
		case PUPIL_BOX_COLOR:			 g = 255; break;//0,255,0
		case SEARCH_LIMIT_BOX_COLOR: 
		case MOUSE_CURSOR_COLOR:		 r = 255; break;//255,0,0
	}

	
	v0 = (unsigned int*) (eyeimage);
	
	radius = width / 2;	
	y= radius;
	f = 1 - radius;
	
	ddF_y = -2 * radius;
	x =  0;
	y = radius;
	x0 = x0 + width/2;
	y0 = y0 + width/2;

	if (x0 - radius > 0 && x0+radius < eyewidth && y0+radius < eyeheight && y0-radius > 0){
		v0[(eyewidth*eyeheight) - (y0+radius)*eyewidth+x0] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		v0[(eyewidth*eyeheight) - (y0-radius)*eyewidth+x0] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		v0[(eyewidth*eyeheight) - y0*eyewidth+(x0+radius)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		v0[(eyewidth*eyeheight) - y0*eyewidth+(x0-radius)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
	}
	
	while( x < y ){	 
		
		if(f >= 0){
			y--;
			ddF_y += 2;
			f += ddF_y;			
		}
		x++;
		ddF_x += 2;
		f += ddF_x;		
		if (y0+y < eyeheight && y0+y > 0){
			if(x0+x < eyewidth)
				v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			if(x0-x  > 0)
				v0[(eyewidth*eyeheight) - (y0+y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		}
		if (y0-y < eyeheight && y0-y > 0){
			if(x0+x < eyewidth)	
				v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0+x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			if(x0-x  > 0)
				v0[(eyewidth*eyeheight) - (y0-y)*eyewidth+(x0-x)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		}
		if(y0+x < eyeheight && y0+x > 0){
			if(x0+y < eyewidth)
				v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			if(x0-y >0)
				v0[(eyewidth*eyeheight) - (y0+x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
		}
		if(y0-x < eyeheight && y0-x > 0){
			if(x0+y < eyewidth)
				v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0+y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			if(x0-y >0)
				v0[(eyewidth*eyeheight) - (y0-x)*eyewidth+(x0-y)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);		
		}
	}
	return;
}
// added by NJ @ SR Research Sept 2010
void drawLozenge(CrossHairInfo *chi, int x0, int y0, int width, int height, int cindex)
{
	
	unsigned char r =0;
	unsigned char g =0;
	unsigned char b =0;
	int x = 0, y;
	int y2, y1;
	unsigned int *v0;
	
	if (eyeimage == NULL) return;
	// clip to prevent memory issues and wrap around
		
	switch(cindex)	{
	
		case CR_HAIR_COLOR:          r=g=b = 255; break;//255,255,255
		case PUPIL_HAIR_COLOR:       r=g=b = 255; break;//255,255,255
		case PUPIL_BOX_COLOR:			 g = 255; break;//0,255,0
		case SEARCH_LIMIT_BOX_COLOR: 
		case MOUSE_CURSOR_COLOR:		 r = 255; break;//255,0,0
	}
	
if(eyeimage != NULL) {
	
	// Retrieve v0 as pointer to pixel row in output buffer:
	v0 = (unsigned int*) (eyeimage);
	
	// is it a circle?
	if(abs(width - height) < 4) 
	{
		drawCircle(chi, x0, y0, width, height, cindex);
			
	}else { // non. ligne
		
		int minwidth = width * (width < height) + height * (width >= height); //min(width,height);
		if (width == minwidth) // width was smaller
		{			
			y1 = y0+width/2;
			x = x0;
			y2 = y1+(height-width) + 1;
			// only draw vertical lines
			for (y=y1 ;y< y2;y++) {
				// be careful to clip or memory error occurs
				if (y < eyeheight && y > 0 && x > 0 && x < eyewidth)
					v0[(eyewidth*eyeheight) - y*eyewidth+x] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				if (y < eyeheight && y > 0 && (x+ width) > 0 && (x+width) < eyewidth)
					v0[(eyewidth*eyeheight) - y*eyewidth+(x+width)] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
				
			}	// Now draw the semi circles
			drawSemiCircle(chi,x0,y0,width,UPSIDE,cindex);
			drawSemiCircle(chi,x0,(y0+width/2+(height-width)),width,DOWNSIDE,cindex);
			
			
		}else{ // height smaller
			
			int x1 = x0+height/2;
			int x2 = x1+(width-height);
			y=y0;			 
			//horizontal			
			for (x=x1 ;x<= x2;x++) {
				if ( x < eyewidth && x > 0 && y > 0 && y < eyeheight)
					v0[(eyewidth*eyeheight) - y*eyewidth+x] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
			
			y = y0+height;
			for (x=x1 ;x<= x2;x++) {
				if ( x < eyewidth && x > 0 && y < eyeheight && y > 0)
					v0[(eyewidth*eyeheight) - y*eyewidth+x] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
						
			drawSemiCircle(chi,(x0+height/2+(width-height)),y0,height,RIGHTSIDE,cindex);
			drawSemiCircle(chi,x0,y0,height,LEFTSIDE,cindex);
			
		}	
	}
}
	return;
}

// added by NJ @ SR Research LTD
void drawLine(CrossHairInfo *chi, int x1, int y1, int x2, int y2, int cindex)
{
	
	unsigned char r =0;
	unsigned char g =0;
	unsigned char b =0;
	int dx, dy;
	int x, y, ch;
	INT16 xc[4],yc[4], enabled;
	unsigned int *v0;
	int xx1, xx2, yy1, yy2;
	
	if (eyeimage == NULL) return;
	// get camera channel. 2 = head.
	ch = get_image_xhair_data(xc, yc, &enabled);
	
	// clip if fail 
	if (ch == 2){
				
		if (x1<0) x1=0;
		if (x2<0) x2=0;
		if (y1<0) y1=0;
		if (y2<0) y2=0;
		
		if (x1>eyewidth-1) x1=eyewidth-1;
		if (x2>eyewidth-1) x2=eyewidth-1;
		if (y1>eyeheight-1) y1=eyeheight-1;
		if (y2>eyeheight-1) y2=eyeheight-1;
	}
	
	
	switch(cindex)
	{
		case CR_HAIR_COLOR:          r=g=b = 255; break;//255,255,255
		case PUPIL_HAIR_COLOR:       r=g=b = 255; break;//255,255,255
		case PUPIL_BOX_COLOR:			 g = 255; break;//0,255,0
		case SEARCH_LIMIT_BOX_COLOR: 
		case MOUSE_CURSOR_COLOR:		 r = 255; break;//255,0,0
	}
	// Memory pointer to malloc()'ed image pixel buffer that holds the
	// image data for a RGBA8 texture with the most recent eye camera image:
	// eyeimage is a global variable
	if(eyeimage != NULL) {
		
		// Retrieve v0 as pointer to pixel row in output buffer: image is upside down
		v0 = (unsigned int*) ( eyeimage);				
		
		dx = x2 - x1;
		dy = y2 - y1;			
		
		
		if(ch != 2){
			// never diagonal here y1 is always < y2
			for (y=y1 ;y< y2;y++) {
				x = x1 + (dx) * (y - y1)/(dy);	
				if (y < eyeheight && x < eyewidth && y >0 && x >0 && (eyewidth*eyeheight) - y*eyewidth+x > 0) 
					v0[(eyewidth*eyeheight) - y*eyewidth+x] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}		
		}
		
		//account for diagonal lines in binocular mode ( x1 and y1 may be > x2 and y2)
		
		// vertical 
		if (dx == 0 ){			
			if(y1>y2){
				yy1 = y2;
				yy2 = y1;
			}else {
				yy1 = y1;
				yy2 = y2;
			}
				
			for (y=yy1 ;y< yy2;y++) {
				x = x1 + (dx) * (y - y1)/(dy);			
				if (y < eyeheight && x < eyewidth && y >0 && x >0 && (eyewidth*eyeheight) - y*eyewidth+x > 0)
					v0[(eyewidth*eyeheight) - y*eyewidth+x] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}		
		}
		else  {	
			if(x1>x2){		
				xx1 = x2;
				xx2 = x1;
			}else {
				xx1 = x1;
				xx2 = x2;
			}
	
			for (x=xx1 ;x< xx2;x++) {
				y = y1 + (dy) * (x - x1)/(dx);
				if (y < eyeheight && x < eyewidth && y >0 && x >0 && (eyewidth*eyeheight) - y*eyewidth+x > 0)
					v0[(eyewidth*eyeheight) - y*eyewidth+x] = 0xFF000000 | ((unsigned int) b << 16) | ((unsigned int) g << 8) | ((unsigned int) r);
			}
		}
	}	
	return;
}


typedef void  (*GET_MOUSE_LOC)(CrossHairInfo *dt, int *x, int *y, int *state);
GET_MOUSE_LOC mouseLoc = NULL;

void  ELCALLTYPE  set_mouse_loc_callback(GET_MOUSE_LOC get_mouse_loc)
{
	mouseLoc = get_mouse_loc;
}


void getMouseState(CrossHairInfo *chi, int *rx, int *ry, int *rstate)
{
	float x =0;
	float y =0;
	PsychGenericScriptType			*inputs[1];
	PsychGenericScriptType			*outputs[1];
	double* callargs;
	double* outputargs;
	float ar[7];
	float w,h;
	int i;
	
	inputs[0]   = mxCreateDoubleMatrix(1, 4, mxREAL);
	callargs    = mxGetPr(inputs[0]);
	
	callargs[0] = 16; // 16 == Command code for mouse button event	

	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		mexSetTrapFlag(1);
	#endif
	
	mexCallMATLAB(1, outputs, 1, inputs, eyelinkDisplayCallbackFunc);

	
	// Reset error handling to default on Matlab:
	#if PSYCH_LANGUAGE == PSYCH_MATLAB
		mexSetTrapFlag(0);
	#endif
	
	outputargs = mxGetData(outputs[0]);
	for (i=0;i<7;i++){
		ar[i] = (int) outputargs[i];
	}
	
	// Release our matrix again:
	mxDestroyArray(inputs[0]);
	mxDestroyArray(outputs[0]);
	
	w = ar[0];
	h = ar[1];
	x = floor((ar[2] - ((w/2) - ar[4]/2)) * ((float)eyewidth/ar[4]));
	y = floor((ar[3] - ((h/2) - ar[5]/2)) * ((float)eyeheight/ar[5]));

	if(x>=0 && y >=0 && x <= eyewidth && y <= eyeheight)
	{
        *rx = (int)x;
		*ry = (int)y;
		*rstate =  (int)ar[6]; 
	}
	return;
}


// PsychEyelink_draw_image_line() retrieves exactly one scanline worth of eye camera
// image data. Once a full image has been received, it has to trigger the actual image
// display:
static void ELCALLBACK PsychEyelink_draw_image_line(INT16 width, INT16 line, INT16 totlines, byte *pixels)
{
	PsychGenericScriptType			*inputs[1];
	PsychGenericScriptType			*outputs[1];
	double* callargs;
	double teximage;
	static INT16 lastline = -1;
	static int wrapcount = 0;
	static double tlastwrap = 0.0;
	double tnow;
	int rc;
	byte* p;
	unsigned int *v0;
	short i;
	CrossHairInfo crossHairInfo;
	
	if (Verbosity() > 8) printf("Eyelink: Entering PsychEyelink_draw_image_line()\n");

	// Callbacks forcefully disabled by error-handling? Simply return with no-op, if so:
	if (0 == eyelinkDisplayCallbackFunc[0]) return;

	// width, line, totlines within valid range?
	if (width < 1 || width > eyewidth || line < 1 || line > eyeheight || totlines < 1 || totlines > eyeheight) {
		printf("EYELINK: WARNING! Eye camera image with invalid parameters received! (width = %i, line = %i, totlines = %i out of sane range %i x %i)!\n",
				width, line, totlines, eyewidth, eyeheight);
		printf("EYELINK: WARNING! Will try to clamp to valid values, but results may be junk.\n");
		width = eyewidth;
		line = (line < 1) ? 1 : line;
		line = (line > eyeheight) ? line : eyeheight;
		totlines = (totlines < 1) ? 1 : totlines;
		totlines = (totlines > eyeheight) ? totlines : eyeheight;
	}

	
	

	// Data structures properly initialized?
	if(eyeimage != NULL) {
		// Retrieve p as pointer to input pixel index color buffer:
		p = pixels;
		
		// Retrieve v0 as pointer to pixel row in output buffer:
		v0 = (unsigned int*) (( eyeimage + ( (totlines - line) * width * 4 ) ));
		
		// Copy one row of pixels from input- to output buffer:
		// This is a bit optimized, but we could do more if we're really bored with life ;-)
		for(i=0; i < width; i++) {
			// Decode pixel index value in *p via palette color lookup table and store to output buffer:
			*(v0++) = palmap32[*p++];
		}

		if (Verbosity() > 8) printf("Eyelink: PsychEyelink_draw_image_line(): Scanline %i received.\n", (int) line);

		// Premature wraparound?
		if (line < lastline) {
			// Premature wraparound due to too slow processing. Increase wrapcounter:
			wrapcount++;			
		}

		// More than some threshold?
		if (wrapcount > 10) {
			// Spill a warning?
			PsychGetAdjustedPrecisionTimerSeconds(&tnow);
			if (tnow - tlastwrap > 2.0) {
				// Last invocation longer than 2 seconds away:
				// Output some warning to console...
				if (Verbosity() > 1) {
					printf("Eyelink: Warning: Skipped videoframes from eye camera detected within last seconds (count=%i)\n", wrapcount);
					printf("Eyelink: Warning: Timing problems on your machine or network problems on tracker connection?!?\n\n");
				}

				// Update / Reset detector:
				tlastwrap = tnow;
				wrapcount = 0;
			}				
		}
		
		// Update skip detector:
		lastline = line;

		// Complete new eye image received?
		if (line == totlines) {
			// Yes. Our eyeimage buffer contains a new image.
			
			// Reset skip detector:
			lastline  = -1;
		
			crossHairInfo.w = eyewidth;
			crossHairInfo.h = eyeheight;
			crossHairInfo.drawLozenge = drawLozenge;
			crossHairInfo.drawLine = drawLine;
			crossHairInfo.getMouseState = mouseLoc?mouseLoc:getMouseState;
			crossHairInfo.userdata = eyeimage;
			
			eyelink_draw_cross_hair(&crossHairInfo);
			
			
			// Compute double-encoded Matlab/Octave compatible memory pointer to image buffer:
			teximage = PsychPtrToDouble((void*) eyeimage);

			
			
			
			// Ok, teximage is a memory pointer to our image buffer, encoded as a double.
			// Now we need to call our Matlab callback function which actually converts
			// the data in our internal image buffer into a PTB texture, then draws that
			// texture etc. to display the new eye camera image.
			if (Verbosity() > 6) printf("Eyelink: PsychEyelink_draw_image_line(): All %i Scanlines received. Calling Runtime!\n", (int) line);
			
			// Create a Matlab double matrix with 4 elements: 1st is command code '1'
			// 2nd is the double pointer, 3r//d is image width, 4th is image height:
			outputs[0]  = NULL;
			inputs[0]   = mxCreateDoubleMatrix(1, 4, mxREAL);
			callargs    = mxGetPr(inputs[0]);

			callargs[0] = 1; // 1 == Command code for "Show eye image".
			callargs[1] = teximage;
			callargs[2] = eyewidth;
			callargs[3] = eyeheight;

			// Call the runtime environment, on Matlab with the trap flag set, so control
			// returns to us on error, instead of returning to Matlab runtime system. Eyelink
			// runtime system doesn't like losing control and would crash otherwise!
			#if PSYCH_LANGUAGE == PSYCH_MATLAB
				mexSetTrapFlag(1);
			#endif
			
			rc = mexCallMATLAB(0, outputs, 1, inputs, eyelinkDisplayCallbackFunc);
			if(rc) {
				printf("EYELINK: WARNING! Failed to call eyelink camera image display callback function %s [rc=%i]!\n", eyelinkDisplayCallbackFunc, rc);
				printf("EYELINK: WARNING! Make sure that function is on your Matlab/Octave path and properly initialized.\n");
				printf("EYELINK: WARNING! May also be an error during execution of that function. Type ple at command prompt for error messages.\n");
				printf("EYELINK: WARNING! Auto-Disabling all callbacks to the runtime environment for safety reasons.\n");
				eyelinkDisplayCallbackFunc[0] = 0;
			}

			// Reset error handling to default on Matlab:
			#if PSYCH_LANGUAGE == PSYCH_MATLAB
				mexSetTrapFlag(0);
			#endif

			// Release our matrix again:
			mxDestroyArray(inputs[0]);
		}
	}
	
	
	
	
	// Done.
	return;
}

// PsychEyelink_set_image_palette() sets the color palette for decoding 1-byte color index
// values in an eyelink camera image into RGB8 1-byte-per-color-component color values:
static void ELCALLBACK PsychEyelink_set_image_palette(INT16 ncolors, byte r[], byte g[], byte b[])
{
	short i;
	
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_set_image_palette()\n");

	if (ncolors > 256) {
		printf("EYELINK: WARNING! Invalid color palette size %i (> 256 colors!) received from eyelink: Clamping to 256 colors.\n", (int) ncolors);
		ncolors = 256;
	}

	// Copy given r,g,b color arrays into internal remapping table palmap32:
	for(i=0; i < ncolors; i++) {
		// Format is ABGR
		palmap32[i] = 0xFF000000 | ((unsigned int) b[i] << 16) | ((unsigned int) g[i] << 8) | ((unsigned int) r[i]);
	}

	return;
}

static INT16  ELCALLBACK PsychEyelink_setup_cal_display(void)
{
	//nj added "hack" to disable flashing instructions in drift correction and to enable sending cal and val results
	int mode = -1;
	
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_setup_cal_display()\n");

	mode = eyelink_tracker_mode();
	
	if (mode == 1 || mode ==9 ) //EL_DRIFT_CORR_MODE)
		if (0xdeadbeef == PsychEyelinkCallRuntime(17, 0, 0, NULL)) {
			// Error condition. Return error to eyelink runtime:
			return(-1);
		}
	
	// Tell runtime to setup calibration display: Command code 7.
	if (0xdeadbeef == PsychEyelinkCallRuntime(7, 0, 0, NULL)) {
		// Error condition. Return error to eyelink runtime:
		return(-1);
	}

	// Return success:
	return(0);
}

static void ELCALLBACK   PsychEyelink_exit_cal_display(void)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_exit_cal_display()\n");

	// Tell runtime to exit calibration display: Command code 11.
	PsychEyelinkCallRuntime(11, 0, 0, NULL);
	return;
}

static void ELCALLBACK   PsychEyelink_clear_display(void)
{
	//NJ modified to add msg to call back 6 with cal and val result
	char strMessage[256];
	int result =-1;
	// Clear strings
	memset(strMessage, 0, sizeof(strMessage));
	
	
	result = eyelink_cal_message(strMessage);
	
	
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_clear_display()\n");

	// Tell runtime to clear display: Command code 6.
	PsychEyelinkCallRuntime(6, 0, 0, strMessage);//NULL);

		
	return;
}

static void ELCALLBACK   PsychEyelink_draw_cal_target(INT16 x, INT16 y)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_draw_cal_target(): x=%i y=%i.\n", (int) x, (int) y);

	// Tell runtime about where to draw calibration target: Command code 5.
	PsychEyelinkCallRuntime(5, (int) x, (int) y, NULL);

	return;
}

static void ELCALLBACK   PsychEyelink_erase_cal_target(void)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_erase_cal_target():\n");

	// Tell runtime about clear calibration target: Command code 10.
	PsychEyelinkCallRuntime(10, 0, 0, NULL);

	return;
}

static void ELCALLBACK   PsychEyelink_image_title(INT16 threshold, char *title)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_image_title(): threshold = %i : Title = %s\n", (int) threshold, title);

	//mexPrintf("C code: %s ... %d\n", title, threshold);
	//fflush(stdout);
	
	// Tell runtime about image title: Command code 4.
	PsychEyelinkCallRuntime(4, (int) threshold, 0, title);
	
	return;
}

#ifndef ELKEY_DOWN
#define ELKEY_DOWN 1 //temporary while we wait for sr-research's lib to get updated with this
#endif

static INT16 ELCALLBACK  PsychEyelink_get_input_key(InputEvent *keyinput)
{
	int ky = 0;
	double tnow;
	static double tlastquery = 0;
	const  double tmininterval = 0.1;	// Allow one query every 0.1 seconds.
	InputEvent *key_input = keyinput;

	// Throttling routine:
	// We don't want key queries to call out to the runtime too often, as this
	// creates a quite significant overhead, e.g., approx. 1 msec for a KbCheck
	// for a fast 2009'ish machine on OS/X!
	PsychGetAdjustedPrecisionTimerSeconds(&tnow);
	if (tnow - tlastquery < tmininterval) {
		// Last invocation less than tmininterval seconds away. Throttle this,
		// we just return "no key pressed".
		if (Verbosity() > 9) printf("Eyelink: In PsychEyelink_get_input_key(): Throttling...\n");
		return(0);
	}
	else {
		// Last invocation longer than tmininterval seconds away. Accept this
		// query and update timestamp:
		tlastquery = tnow;
	}
	
	if (Verbosity() > 7) printf("Eyelink: Entering PsychEyelink_get_input_key()\n");

	// Call runtime for keycode of pressed key (command code 2):
	if ((ky = PsychEyelinkCallRuntime(2, 0, 0, NULL)) == 0xdeadbeef) {
		// Error condition in runtime callback! Can't progress. We try to
		// shutdown the current eyelink runtime operation by sending a fake
		// keycode corresponding to the terminate key:
		if (Verbosity() > 0) printf("Eyelink: In PsychEyelink_get_input_key(): Error condition detected: Trying to send TERMINATE_KEY abort keycode!\n");
		ky = TERMINATE_KEY;
	}
	
	if (ky > 0) {
		// Fill Eyelinks InputEvent struct:
		memset(key_input, 0, sizeof(InputEvent));
		key_input->key.key = ky;
		key_input->key.state = ELKEY_DOWN; 
		key_input->key.type = KEYINPUT_EVENT;
		key_input->key.modifier = 0; //event.key.keysym.mod;
		key_input->key.unicode = 0;//event.key.keysym.unicode;
		
		// One key pressed:
		return(1);
	}
	else {
		// No key pressed:
		return(0);
	}
}

static void ELCALLBACK   PsychEyelink_alert_printf_hook(const char *msg)
{
	// Print error message to runtime console if error output is allowed:
	if (Verbosity() > 3) printf("Eyelink: Alert! Eyelink says: %s.\n\n", msg);
	
	// Tell runtime about alert condition: Command code 3.
	PsychEyelinkCallRuntime(3, 0, 0, (char*) msg);

	return;
}

static void ELCALLBACK   PsychEyelink_cal_target_beep_hook(void)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_cal_target_beep_hook():\n");

	PsychEyelinkCallRuntime(12, 0, 0, NULL);
	return;
}

static void ELCALLBACK   PsychEyelink_dc_target_beep_hook(void)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_dc_target_beep_hook():\n");

	PsychEyelinkCallRuntime(13, 0, 0, NULL);
	return;
}

static void ELCALLBACK	  PsychEyelink_cal_done_beep_hook(INT16 error)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_cal_done_beep_hook():\n");

	PsychEyelinkCallRuntime(14, (int) error, 0, NULL);
	return;
}

static void ELCALLBACK	  PsychEyelink_dc_done_beep_hook(INT16 error)
{
	if (Verbosity() > 5) printf("Eyelink: Entering PsychEyelink_dc_done_beep_hook():\n");

	PsychEyelinkCallRuntime(15, (int) error, 0, NULL);
	return;
}
