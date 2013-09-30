/*
	SCREENColorRange.c		

	AUTHORS:

		mario.kleiner at tuebingen.mpg.de 		mk 
  
	PLATFORMS:	
	
		All.

	HISTORY:

		12/05/06	mk		Created.

	DESCRIPTION:
  
		Set color scaling and clamping behaviour for a specific onscreen window.
		
	TO DO:  

*/


#include "Screen.h"

char fragmentTunnelSrc[] =
"/* Simple pass-through fragment shader: Just draws fragment at requested */ \n"
"/* position, but with the special unclampedFragColor instead of the FragColor: */ \n"
"\n"
"varying vec4 unclampedFragColor;\n"
"\n"
"void main()\n"
"{\n"
"    /* Simply copy input unclamped RGBA pixel color into output fragment color: */\n"
"    gl_FragColor = unclampedFragColor;\n"
"}\n\0";

char vertexTunnelSrc[] =
"/* Simple pass-through vertex shader: Emulates fixed function pipeline, but passes  */ \n"
"/* gl_MultiTexCoord0 as varying unclampedFragColor to circumvent vertex color       */ \n"
"/* clamping on gfx-hardware / OS combos that don't support unclamped operation:     */ \n"
"/* PTBs color handling is expected to pass the vertex color in gl_MultiTexCoord0    */ \n"
"/* for unclamped drawing for this reason. */ \n"
"\n"
"varying vec4 unclampedFragColor;\n"
"\n"
"void main()\n"
"{\n"
"    /* Simply copy input unclamped RGBA pixel color into output varying color: */\n"
"    unclampedFragColor = gl_MultiTexCoord0;\n"
"\n"
"    /* Set real fixed function vertex color to red -- helps debugging: */\n"
"    gl_FrontColor  = vec4(1.0, 0.0, 0.0, 1.0);\n"
"\n"
"    /* Output position is the same as fixed function pipeline: */\n"
"    gl_Position    = ftransform();\n"
"}\n\0";

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "[oldmaximumvalue, oldclampcolors, oldapplyToDoubleInputMakeTexture] = Screen('ColorRange', windowPtr [, maximumvalue][, clampcolors][, applyToDoubleInputMakeTexture]);";
static char synopsisString[] = 
	"Set or return the maximum color component value that PTB should allow for provided color "
	"values when drawing into a window 'windowPtr' or its associated Offscreen windows. "
	"'maximumvalue' is the optional new setting for the maximum allowed color component value: "
	"PTB expects the values of provided color components (red, green, blue, alpha or intensity) "
	"to be in the range 0 to maximumvalue. 0 is mapped to minimum output intensity, "
	"maximumvalue is mapped to maximum output intensity, values outside the range 0-maximumvalue "
	"are automatically saturated (clamped) at zero or maximumvalue if clamping is enabled. "
	"Initially color range clamping is enabled and the maximumvalue defaults to the biggest "
	"integral number displayable by your video hardware, e.g., 255 for a standard 8 bit "
	"per color component framebuffer as present on most consumer graphics hardware. \n"
	"A maximumvalue == 1.0 has special meaning: A maximumvalue of 1.0 will enable PTB to pass "
	"color values in OpenGL's native floating point color range of 0.0 to 1.0: This has two "
	"advantages: First, your color values are independent of display device depth, i.e. no need "
	"to rewrite your code when running it on higher resolution hardware. Second, PTB can skip any "
	"color range remapping operations - this can speed up drawing significantly in some cases. "
	"'clampcolors': By default, OpenGL clamps colors to the range 0-maximumvalue. Negative colors "
	"are clamped to zero, color values greater than maximumvalue are clamped to maximumvalue. If "
	"you set the optional flag 'clampcolors' to 0, PTB will disable color clamping on hardware that "
	"supports unclamped colors. This allows to pass arbitrary (even negative) floating point "
	"numbers as color values. This is useful in conjunction with special high dynamic range display "
	"devices and for certain image processing operations when using PTB's image processing pipeline. "
	"If your hardware/operating system does support shaders, but not unclamped colors, Screen will "
	"try to use a shader-based workaround to enable unclamped color processing despite missing "
	"hardware capabilities - This comes at some performance penalty. You can also force Screen "
	"to always use its own unclamped color implementation by setting 'clampcolors' to a value of "
	"-1. On some graphics hardware, this may increase the precision with which the color of "
	"drawn objects is handled, but again at some speed penalty.\n"
	"CAUTION: If you change the color range or clamping of an onscreen window with this function, the change will "
	"only affect textures and offscreen windows created *after* this function call, not ones created "
	"before. It's therefore recommended to execute this function immediately after creating an onscreen "
	"window to guarantee consistent behaviour of your code. Color values provided as uint8 arrays or as "
	"textures, e.g., from video capture, movies or Screen('MakeTexture') are not rescaled, as"
	"they are (expected) to be in a proper format for the given color depth already. However, the "
	"Screen('MakeTexture') command has an optional flag 'floatprecision' that allows you to pass image "
	"matrices unclamped and with either 16 bpc or 32 bpc floating point color precision if you want.\n"
	"Additionally you can force Screen('MakeTexture') to apply the 'maximumvalue' setting to regular "
	"textures which are provided as Matlab double type matrices by setting the optional parameter "
	"'applyToDoubleInputMakeTexture' to 1. This will allow an input range of double values between "
	"zero and 'maximumvalue' to make your code more consistent, but it will still limit the range of "
	"valid values to that range and will represent that range only with 8 bit for 256 different levels, "
	"ie., the description of the 'clampcolors' setting do not apply to such uint8 low-precision textures! "
	"For full precision unconstrained textures, you'll still need to set the 'floatprecision' flag accordingly, "
	"as such textures require more memory and processing resources.";

static char seeAlsoString[] = "OpenWindow OpenOffscreenWindow";	 

PsychError SCREENColorRange(void) 
{
	PsychWindowRecordType *windowRecord;
	double maxvalue, clampcolors, oldclampcolors;
	int applyToMakeTexture;
	GLboolean enabled, enabled1, enabled2, enabled3;
	GLuint tunnelShader = 0;

	//all subfunctions should have these two lines.  
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	PsychErrorExit(PsychCapNumInputArgs(4));     //The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1)); //The required number of inputs	
	PsychErrorExit(PsychCapNumOutputArgs(3));	 //The maximum number of outputs
        
	//get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
	if(!PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_user, "Tried to set color range on something else than a onscreen window");
            
	// Decode current setting:
	maxvalue = windowRecord->colorRange;
	clampcolors = (maxvalue > 0) ? 1 : 0;
	maxvalue = fabs(maxvalue);
	oldclampcolors = clampcolors;
	
	// Copy out optional return value:
	PsychCopyOutDoubleArg(1, FALSE, maxvalue);
	PsychCopyOutDoubleArg(2, FALSE, clampcolors);
	PsychCopyOutDoubleArg(3, FALSE, windowRecord->applyColorRangeToDoubleInputMakeTexture);

	// Get the optional new values:
	PsychCopyInDoubleArg(2, FALSE, &maxvalue);
	PsychCopyInDoubleArg(3, FALSE, &clampcolors);
	if (PsychCopyInIntegerArg(4, FALSE, &applyToMakeTexture)) {
		if ((applyToMakeTexture < 0) || (applyToMakeTexture > 1)) PsychErrorExitMsg(PsychError_user, "Tried to set invalid applyToDoubleInputMakeTexture flag (valid values are 0 and 1).");
		windowRecord->applyColorRangeToDoubleInputMakeTexture = applyToMakeTexture;
	}

	// Encode into one value:
	if (maxvalue<=0) PsychErrorExitMsg(PsychError_user, "Tried to set invalid color maximumvalue (negative or zero).");
		
	// Switch of clamping mode needed?
	if (oldclampcolors != clampcolors) {
		// Yes. Try to set new clamping behaviour:
		PsychSetGLContext(windowRecord);

		// Does graphics hardware/OS support clamping mode change via glClampColorARB and shall
		// we use it? A clampcolors setting of -1 would use our own shader based implementation
		// even if the hardware could do it -- This to test the precision of our approach vs.
		// hardware and to guarantee consistent results even if it means a performance hit.
		if (glClampColorARB && (clampcolors>=0)) {
			// Color clamping extension supported: Set new clamp mode in hardware via extension:
			
			enabled = (clampcolors > 0) ? GL_TRUE : GL_FALSE;
			glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, enabled);
			glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, enabled);
			glClampColorARB(GL_CLAMP_READ_COLOR_ARB, enabled);

			// Check if the clamp en-/disable worked:
			glGetBooleanv(GL_CLAMP_VERTEX_COLOR_ARB, &enabled1);
			glGetBooleanv(GL_CLAMP_FRAGMENT_COLOR_ARB, &enabled2);
			glGetBooleanv(GL_CLAMP_READ_COLOR_ARB, &enabled3);
			
			if ((clampcolors==0 && (enabled1 || enabled2 || enabled3)) || (clampcolors>0 && (!enabled1 || !enabled2 || !enabled3))) {
				if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Could not %s color value clamping as requested. Unsupported by your graphics hardware?\n", (clampcolors>0) ? "enable" : "disable");

				// Reset to old setting if the switch didn't work:
				clampcolors = oldclampcolors;
			}
			else if ((PsychPrefStateGet_Verbosity() > 3) && (clampcolors==0)) {
				printf("PTB-INFO: Disabled color clamping via hardware.\n");
			}
		}
		else {
			// Color clamping extensions unsupported, or user wants our own implementation: We need to use quite a bit of shader and cpu magic...
			if ((PsychPrefStateGet_Verbosity() > 3) && (clampcolors>=0)) printf("PTB-INFO: Switching and query of color clamping via glClampColorARB unsupported by your graphics hardware or operating system.\n");
			if ((PsychPrefStateGet_Verbosity() > 3) && (clampcolors< 0)) printf("PTB-INFO: Switching of color clamping via internal shader-based solution forcefully enabled by usercode.\n");

			// Not a big deal if clamping is supposed to be enabled, only if clamping shall be disabled:
			if ((clampcolors <= 0) && (windowRecord->unclampedDrawShader == 0)) {
				// Need to workaround by using a GLSL shader, but none yet created:
				
				// We create a GLSL vertex-fragment-shader combo that accepts unclamped 
				// colors as standard vertex attribute and "tunnels" that down to the 
				// raster-backends. Obviously only works with GLSL capable hardware and
				// only for internal drawing commands where we can bind the shader whenever needed:
				tunnelShader = PsychCreateGLSLProgram(fragmentTunnelSrc, vertexTunnelSrc, NULL);
				if (tunnelShader) {
					// Got a shader :-) -- Assign it as color clamping shader to onscreen windowRecord:
					windowRecord->unclampedDrawShader = tunnelShader;
					if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Color clamping disabled via internal shader-based solution.\n");
				}
				else {
					// Failed: This is a no-go:
					if (PsychPrefStateGet_Verbosity()>1) {
						printf("PTB-WARNING: Switching and query of color clamping via glClampColorARB unsupported by your graphics hardware.\n");
						printf("PTB-WARNING: Unfortunately the GLSL shader-based workaround didn't work either. Disable of color clamping not possible.\n");
					}

					// Reset to old setting if the switch didn't work:
					clampcolors = oldclampcolors;
				}
			}
			
			// Switch clamping on/off by shader assignment:
			if (clampcolors <= 0) {
				// Clamping off: Assign our unclamped drawing shader as default shader to be
				// used for drawing whenever no other special shader / override shader is bound:
				windowRecord->defaultDrawShader = windowRecord->unclampedDrawShader;
			}
			else {
				// Clamped drawing: Use fixed function pipeline:
				windowRecord->defaultDrawShader = 0;
			}
			// End of workaround code for missing glClampColorARB:
		}
	}
	
	// Encode maxcolor as well as new clamping mode (in sign):
	windowRecord->colorRange = (clampcolors>0) ? maxvalue : -maxvalue;
	
	// Successfully done:
	return(PsychError_none);
}
