/*
	PsychToolbox3/Source/Common/Screen/PsychImagingPipelineSupport.c
	
	PLATFORMS:	
	
		All. Well, all with sufficiently advanced graphics hardware...
				
	AUTHORS:
		
		Mario Kleiner           mk              mario.kleiner at tuebingen.mpg.de

	HISTORY:
	
		12/05/06	mk	Wrote it.
		
	DESCRIPTION:
	
		Infrastructure for all Screen imaging pipeline functions, i.e., hook callback functions and chains
		and GLSL based internal image processing pipeline.
		
		The level of support for PTB's imaging pipe strongly depends on the capabilities of the gfx-hardware,
		especially GLSL support, shader support and Framebuffer object support.
		
	NOTES:
	
	TO DO: 
	
		We could implement a less capable minimal fallback-path for gfx-hardware that doesn't support
		framebuffer objects, but does support rectangle textures. Such hw would be limited to the initial
		pipeline stage (image processing) and would only allow single-pass processing, no multi-pass algs
		due to the lack of bounce buffers. But it could still allow for slightly more precise textures,
		and simple shaders (if supported) or basic fixed-function processing, e.g., simple contrast
		corrections or geometric undistortion.

*/

#include "Screen.h"

static char texturePlanar1FragmentShaderSrc[] =
"\n"
" \n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image; \n"
"varying vec4 unclampedFragColor; \n"
"varying vec2 texNominalSize; \n"
" \n"
"void main() \n"
"{ \n"
"    vec4 texcolor; \n"
"    texcolor.rgb = vec3(texture2DRect(Image, gl_TexCoord[0].st).r); \n"
"    texcolor.a   = 1.0; \n"
"\n"
"    /* Multiply texcolor with incoming fragment color (GL_MODULATE emulation): */ \n"
"    /* Assign result as output fragment color: */ \n"
"    gl_FragColor = texcolor * unclampedFragColor; \n"
"} \n";

static char texturePlanar2FragmentShaderSrc[] =
"\n"
" \n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image; \n"
"varying vec4 unclampedFragColor; \n"
"varying vec2 texNominalSize; \n"
" \n"
"void main() \n"
"{ \n"
"    vec4 texcolor; \n"
"    texcolor.rgb = vec3(texture2DRect(Image, gl_TexCoord[0].st).r); \n"
"    texcolor.a   = texture2DRect(Image, gl_TexCoord[0].st + vec2(0.0, 1.0 * texNominalSize.y)).r; \n"
"\n"
"    /* Multiply texcolor with incoming fragment color (GL_MODULATE emulation): */ \n"
"    /* Assign result as output fragment color: */ \n"
"    gl_FragColor = texcolor * unclampedFragColor; \n"
"} \n";

static char texturePlanar3FragmentShaderSrc[] =
"\n"
" \n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image; \n"
"varying vec4 unclampedFragColor; \n"
"varying vec2 texNominalSize; \n"
" \n"
"void main() \n"
"{ \n"
"    vec4 texcolor; \n"
"    texcolor.r = texture2DRect(Image, gl_TexCoord[0].st).r; \n"
"    texcolor.g = texture2DRect(Image, gl_TexCoord[0].st + vec2(0.0, 1.0 * texNominalSize.y)).r; \n"
"    texcolor.b = texture2DRect(Image, gl_TexCoord[0].st + vec2(0.0, 2.0 * texNominalSize.y)).r; \n"
"    texcolor.a = 1.0; \n"
"\n"
"    /* Multiply texcolor with incoming fragment color (GL_MODULATE emulation): */ \n"
"    /* Assign result as output fragment color: */ \n"
"    gl_FragColor = texcolor * unclampedFragColor; \n"
"} \n";

static char texturePlanar4FragmentShaderSrc[] =
"\n"
" \n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image; \n"
"varying vec4 unclampedFragColor; \n"
"varying vec2 texNominalSize; \n"
" \n"
"void main() \n"
"{ \n"
"    vec4 texcolor; \n"
"    texcolor.r = texture2DRect(Image, gl_TexCoord[0].st).r; \n"
"    texcolor.g = texture2DRect(Image, gl_TexCoord[0].st + vec2(0.0, 1.0 * texNominalSize.y)).r; \n"
"    texcolor.b = texture2DRect(Image, gl_TexCoord[0].st + vec2(0.0, 2.0 * texNominalSize.y)).r; \n"
"    texcolor.a = texture2DRect(Image, gl_TexCoord[0].st + vec2(0.0, 3.0 * texNominalSize.y)).r; \n"
"\n"
"    /* Multiply texcolor with incoming fragment color (GL_MODULATE emulation): */ \n"
"    /* Assign result as output fragment color: */ \n"
"    gl_FragColor = texcolor * unclampedFragColor; \n"
"} \n";

/* Sampling and conversion shader from YUV-I420 planar format to standard RGBA8
 * format. Samples our YUV I420 planar luminance texture, builds yuv sample, then
 * performs color space conversion from yuv to rgb.
 *
 * This shader is based on BSD licensed example code from Peter Bengtsson, from
 * http://www.fourcc.org/source/YUV420P-OpenGL-GLSLang.c
 */
static char texturePlanarI420FragmentShaderSrc[] =
"/* YUV-I420 planar texture sampling fragment shader.               */ \n"
"/* Retrieves YUV sample from proper locations in planes.           */ \n"
"/* Converts YUV sample to RGB color triplet and applies            */ \n"
"/* GL_MODULATE texture function emulation before fragment output.  */ \n"
"\n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image; \n"
"varying vec4 unclampedFragColor; \n"
"varying vec2 texNominalSize; \n"
" \n"
"void main() \n"
"{ \n"
"    float r, g, b, y, u, v;\n"
"    float nx, ny;\n"
"\n"
"    nx = gl_TexCoord[0].x;\n"
"    ny = gl_TexCoord[0].y;\n"
"\n"
"    y = texture2DRect(Image, vec2(nx, ny)).r;\n"
"    ny = floor(ny * 0.5);\n"
"    nx = floor(nx * 0.5);\n"
"    if (mod(ny, 2.0) >= 0.5) {\n"
"        nx += texNominalSize.x * 0.5;\n"
"    }\n"
"\n"
"    ny = (ny - mod(ny, 2.0)) * 0.5;\n"
"    u = texture2DRect(Image, vec2(nx, ny + texNominalSize.y)).r; \n"
"    v = texture2DRect(Image, vec2(nx, ny + (1.25 * texNominalSize.y))).r; \n"
"\n"
"    y = 1.1643 * (y - 0.0625);\n"
"    u = u - 0.5;\n"
"    v = v - 0.5;\n"
"\n"
"    r = y + 1.5958 * v;\n"
"    g = y - 0.39173 * u - 0.81290 * v;\n"
"    b = y + 2.017 * u;\n"
"\n"
"    /* Multiply texcolor with incoming fragment color (GL_MODULATE emulation): */ \n"
"    /* Assign result as output fragment color: */ \n"
"    gl_FragColor = vec4(r, g, b, 1.0) * unclampedFragColor; \n"
"} \n";

/* Sampling and conversion shader from Y8-I800 planar format to standard RGBA8
 * format. Samples our Y8 I800 planar luminance texture, then performs color
 * space conversion from y to rgb.
 */
static char texturePlanarI800FragmentShaderSrc[] =
"/* Y8-I800 planar texture sampling fragment shader.                */ \n"
"/* Retrieves Y8 sample from proper location.                       */ \n"
"/* Converts Y8 sample to RGB color triplet and applies             */ \n"
"/* GL_MODULATE texture function emulation before fragment output.  */ \n"
"\n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image; \n"
"varying vec4 unclampedFragColor; \n"
" \n"
"void main() \n"
"{ \n"
"    float y = texture2DRect(Image, gl_TexCoord[0].xy).r;\n"
"    y = 1.1643 * (y - 0.0625);\n"
"\n"
"    /* Multiply texcolor with incoming fragment color (GL_MODULATE emulation): */ \n"
"    /* Assign result as output fragment color: */ \n"
"    gl_FragColor = vec4(y, y, y, 1.0) * unclampedFragColor; \n"
"} \n";

char texturePlanarVertexShaderSrc[] =
"/* Simple pass-through vertex shader: Emulates fixed function pipeline, but passes  */ \n"
"/* modulateColor as varying unclampedFragColor to circumvent vertex color       */ \n"
"/* clamping on gfx-hardware / OS combos that don't support unclamped operation:     */ \n"
"/* PTBs color handling is expected to pass the vertex color in modulateColor    */ \n"
"/* for unclamped drawing for this reason. */ \n"
"\n"
"varying vec4 unclampedFragColor;\n"
"varying vec2 texNominalSize;\n"
"attribute vec4 modulateColor;\n"
"attribute vec4 sizeAngleFilterMode;\n"
"\n"
"void main()\n"
"{\n"
"    /* Simply copy input unclamped RGBA pixel color into output varying color: */\n"
"    unclampedFragColor = modulateColor;\n"
"    texNominalSize = sizeAngleFilterMode.xy;\n"
"\n"
"    gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0;\n"
"\n"
"    /* Output position is the same as fixed function pipeline: */\n"
"    gl_Position    = ftransform();\n"
"}\n\0";

// Source code for a fragment shader that performs texture lookup and
// modulation, but taking the modulatecolor from 'unclampedFragColor'
// instead of standard interpolated fragment color. This is used when
// high-precision and unclamped fragment input colors are needed and
// the hw doesn't support that. It's a drop in replacement for fixed
// function pipe otherwise:
static char textureLookupFragmentShaderSrc[] =
"\n"
" \n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image; \n"
"varying vec4 unclampedFragColor; \n"
" \n"
"void main() \n"
"{ \n"
"    vec4 texcolor = texture2DRect(Image, gl_TexCoord[0].st); \n"
"    /* Multiply texcolor with incoming fragment color (GL_MODULATE emulation): */ \n"
"    /* Assign result as output fragment color: */ \n"
"    gl_FragColor = texcolor * unclampedFragColor; \n"
"} \n";

// Source code for a fragment shader that performs bilinear texture filtering.
// This shader is used as a drop-in replacement for GL's GL_LINEAR built-in
// texture filter, whenever that filter is not available: All pre ATI Radeon HD
// hardware and all pre GF6000 NVidia hardware can't filter float textures,
// GF6/7 series NVidia hardware can only filter 16bpc floats, not 32bpc floats.
static char textureBilinearFilterFragmentShaderSrc[] =
"\n"
" \n"
" \n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image; \n"
"varying vec4 unclampedFragColor; \n"
" \n"
"void main() \n"
"{ \n"
"    /* Get wanted texture coordinate for which we should filter: */ \n"
"    vec2 texinpos = (gl_TexCoord[0].st) - vec2(0.5, 0.5); \n"
"    /* Retrieve texel colors for 4 nearest neighbours: */ \n"
"    vec4 tl=texture2DRect(Image, floor(texinpos)); \n"
"    vec4 tr=texture2DRect(Image, floor(texinpos) + vec2(1.0, 0.0)); \n"
"    vec4 bl=texture2DRect(Image, floor(texinpos) + vec2(0.0, 1.0)); \n"
"    vec4 br=texture2DRect(Image, floor(texinpos) + vec2(1.0, 1.0)); \n"
"    /* Perform weighted linear interpolation -- bilinear interpolation of the 4: */ \n"
"    tl=mix(tl,tr,fract(texinpos.x)); \n"
"    bl=mix(bl,br,fract(texinpos.x)); \n"
"    vec4 texcolor = mix(tl, bl, fract(texinpos.y)); \n"
"    /* Multiply filtered texcolor with incoming fragment color (GL_MODULATE emulation): */ \n"
"    /* Assign result as output fragment color: */ \n"
"    gl_FragColor = texcolor * unclampedFragColor; \n"
"} \n";

char textureBilinearFilterVertexShaderSrc[] =
"/* Simple pass-through vertex shader: Emulates fixed function pipeline, but passes  */ \n"
"/* modulateColor as varying unclampedFragColor to circumvent vertex color       */ \n"
"/* clamping on gfx-hardware / OS combos that don't support unclamped operation:     */ \n"
"/* PTBs color handling is expected to pass the vertex color in modulateColor    */ \n"
"/* for unclamped drawing for this reason. */ \n"
"\n"
"varying vec4 unclampedFragColor;\n"
"attribute vec4 modulateColor;\n"
"\n"
"void main()\n"
"{\n"
"    /* Simply copy input unclamped RGBA pixel color into output varying color: */\n"
"    unclampedFragColor = modulateColor;\n"
"\n"
"    gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0;\n"
"\n"
"    /* Output position is the same as fixed function pipeline: */\n"
"    gl_Position    = ftransform();\n"
"}\n\0";



// Source code for our GLSL anaglyph stereo shader:
char anaglyphshadersrc[] = 
"/* Weight vector for conversion from RGB to Luminance, according to NTSC spec. */ \n"
" \n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform vec3 ColorToGrayWeights; \n"
"/* Bias to add to final product - The background color (normally (0,0,0)). */ \n"
"uniform vec3 ChannelBias; \n"
"/* Left image channel and right image channel: */ \n"
"uniform sampler2DRect Image1; \n"
"uniform sampler2DRect Image2; \n"
"uniform vec3 Gains1;\n"
"uniform vec3 Gains2;\n"
"\n"
"void main()\n"
"{\n"
"    /* Lookup RGBA pixel colors in left- and right buffer and convert to Luminance */\n"
"    vec3 incolor1 = texture2DRect(Image1, gl_TexCoord[0].st).rgb;\n"
"    float luminance1 = dot(incolor1, ColorToGrayWeights);\n"
"    vec3 incolor2 = texture2DRect(Image2, gl_TexCoord[0].st).rgb;\n"
"    float luminance2 = dot(incolor2, ColorToGrayWeights);\n"
"    /* Replicate in own RGBA tupel */\n"
"    vec3 channel1 = vec3(luminance1);\n"
"    vec3 channel2 = vec3(luminance2);\n"
"    /* Mask with per channel weights: */\n"
"    channel1 = channel1 * Gains1;\n"
"    channel2 = channel2 * Gains2;\n"
"    /* Add them up to form final output fragment: */\n"
"    gl_FragColor.rgb = channel1 + channel2 + ChannelBias;\n"
"    /* Alpha is forced to 1 - It does not matter anymore: */\n"
"    gl_FragColor.a = 1.0;\n"
"}\n\0";

char passthroughshadersrc[] =
" \n"
"#extension GL_ARB_texture_rectangle : enable \n"
" \n"
"uniform sampler2DRect Image1; \n"
"\n"
"void main()\n"
"{\n"
"    gl_FragColor.rgb = texture2DRect(Image1, gl_TexCoord[0].st).rgb;\n"
"    gl_FragColor.a = 1.0;\n"
"}\n\0";

// This array maps hook point name strings to indices. The symbolic constants in
// PsychImagingPipelineSupport.h define symbolic names for the indices for fast
// lookup by name:
#define MAX_HOOKNAME_LENGTH 40
#define MAX_HOOKSYNOPSIS_LENGTH 1024

char PsychHookPointNames[MAX_SCREEN_HOOKS][MAX_HOOKNAME_LENGTH] = {
	"CloseOnscreenWindowPreGLShutdown",
	"CloseOnscreenWindowPostGLShutdown",
	"UserspaceBufferDrawingFinished",
	"StereoLeftCompositingBlit",
	"StereoRightCompositingBlit",
	"StereoCompositingBlit",
	"PostCompositingBlit",
	"FinalOutputFormattingBlit",
	"UserspaceBufferDrawingPrepare",
	"IdentityBlitChain",
	"LeftFinalizerBlitChain",
	"RightFinalizerBlitChain",
	"UserDefinedBlit",
	"FinalOutputFormattingBlit0",
	"FinalOutputFormattingBlit1",
	"ScreenFlipImpliedOperations",
	"PreSwapbuffersOperations"
};

char PsychHookPointSynopsis[MAX_SCREEN_HOOKS][MAX_HOOKSYNOPSIS_LENGTH] = {
	"OpenGL based actions to be performed when an onscreen window is closed, e.g., teardown for special output devices.",
	"Non-graphics actions to be performed when an onscreen window is closed, e.g., teardown for special output devices.",
	"Operations to be performed after last drawing command, i.e. in Screen('Flip') or Screen('DrawingFinshed').",
	"Perform generic user-defined image processing on image content of left-eye (or mono) buffer.",
	"Perform generic user-defined image processing on image content of right-eye buffer.",
	"Internal(preinitialized): Compose left- and right-eye view into one combined image for all stereo modes except quad-buffered flip-frame stereo.",
	"Not yet used.",
	"Perform post-processing indifferent of stereo mode, e.g., special data formatting for devices like BrightSideHDR, Bits++, Video attenuators...",
	"Operations to be performed immediately after Screen('Flip') in order to prepare drawing commands of users script.",
	"Internal(preinitialized): Only for internal use. Only modify for debugging and testing of pipeline itself!",
	"Internal(preinitialized): Perform last time operation on left (or mono) channel, e.g., draw blue-sync lines.",
	"Internal(preinitialized): Perform last time operation on right channel, e.g., draw blue-sync lines.",
	"Defines a user defined image processing operation for the Screen('TransformTexture') command.",
	"Perform post-processing on 1st (left) output channel, e.g., special data formatting for devices like BrightSideHDR, Bits++, Video attenuators...",
	"Perform post-processing on 2nd (right) output channel, e.g., special data formatting for devices like BrightSideHDR, Bits++, Video attenuators...",
	"This is called after a bufferswap has truly completed by Screen('Flip'), Screen('AsyncFlipCheckEnd') or Screen('AsyncFlipEnd'). Use of OpenGL commands and Matlab commands is allowed here.",
	"Called before emitting the PsychOSFlipWindowBuffers() call, ie., less than 1 video refresh from flipdeadline away. No OpenGL ops allowed!"
};

/* PsychInitImagingPipelineDefaultsForWindowRecord()
 * Set "all-off" defaults in windowRecord. This is called during creation of a windowRecord.
 * It sets all imaging related fields to safe defaults.
 */
void PsychInitImagingPipelineDefaultsForWindowRecord(PsychWindowRecordType *windowRecord)
{
	int i;
	
	// Initialize everything to "all-off" default:
	for (i=0; i<MAX_SCREEN_HOOKS; i++) {
		windowRecord->HookChainEnabled[i]=FALSE;
		windowRecord->HookChain[i]=NULL;
	}
	
	// Disable all special framebuffer objects by default:
	windowRecord->drawBufferFBO[0]=-1;
	windowRecord->drawBufferFBO[1]=-1;
	windowRecord->inputBufferFBO[0]=-1;
	windowRecord->inputBufferFBO[1]=-1;
	windowRecord->processedDrawBufferFBO[0]=-1;
	windowRecord->processedDrawBufferFBO[1]=-1;
	windowRecord->processedDrawBufferFBO[2]=-1;
	windowRecord->preConversionFBO[0]=-1;
	windowRecord->preConversionFBO[1]=-1;
	windowRecord->preConversionFBO[2]=-1;
	windowRecord->finalizedFBO[0]=-1;
	windowRecord->finalizedFBO[1]=-1;
	windowRecord->fboCount = 0;

	// NULL-out fboTable:
	for (i=0; i<MAX_FBOTABLE_SLOTS; i++) windowRecord->fboTable[i] = NULL;
	
	// Setup mode switch in record to "all off":
	windowRecord->imagingMode = 0;

	return;
}

/*  PsychInitializeImagingPipeline()
 *
 *  Initialize imaging pipeline for windowRecord, applying the imagingmode flags. Called by Screen('OpenWindow').
 *
 *  This routine performs initial setup of an imaging pipeline for an onscreen window. It sets up reasonable
 *  default values in the windowRecord (imaging pipe is disabled by default if imagingmode is zero), based on
 *  the imagingmode flags and all the windowRecord and OpenGL settings.
 *
 *  1. FBO's are setup according to the requested imagingmode, stereomode and color depth of a window.
 *  2. Depending on stereo mode and imagingmode, some default GLSL shaders may get created and attached to
 *     some hook-chains for advanced stereo processing.
 */
void PsychInitializeImagingPipeline(PsychWindowRecordType *windowRecord, int imagingmode, int multiSample)
{	
	GLenum fboInternalFormat, finalizedFBOFormat;
	int newimagingmode = 0;
	int fbocount = 0;
	int winwidth, winheight;
    int clientwidth, clientheight;
	psych_bool needzbuffer, needoutputconversion, needimageprocessing, needseparatestreams, needfastbackingstore, targetisfinalFB;
	GLuint glsl;
	GLint redbits;
	float rg, gg, bg;	// Gains for color channels and color masking for anaglyph shader setup.
	char blittercfg[1000];

	// Processing ends here after minimal "all off" setup, if pipeline is disabled:
	if (imagingmode<=0) {
		imagingmode=0;
		return;
	}
		
	// Activate rendering context of this window:
	PsychSetGLContext(windowRecord);

	// Check if this system does support OpenGL framebuffer objects and rectangle textures:
	if (!(windowRecord->gfxcaps & kPsychGfxCapFBO)) {
		// Unsupported! This is a complete no-go :(
		printf("PTB-ERROR: Initialization of the built-in image processing pipeline failed. Your graphics hardware or graphics driver does not support\n");
		printf("PTB-ERROR: the required OpenGL framebuffer object extension. You may want to upgrade to the latest drivers or if that doesn't help, to a\n");
		printf("PTB-ERROR: more recent graphics card. You'll need at minimum a NVidia GeforceFX-5000 or a ATI Radeon 9600 or Intel GMA 950 for this to work.\n");
		printf("PTB-ERROR: See the www.psychtoolbox.org Wiki for recommendations. You can still use basic stereo support (with restricted performance and features)\n");
		printf("PTB-ERROR: by disabling the imaging pipeline (imagingmode = 0) but still selecting a stereomode in the 'OpenWindow' subfunction.\n");
		PsychErrorExitMsg(PsychError_user, "Imaging Pipeline setup: Sorry, your graphics card does not meet the minimum requirements for use of the imaging pipeline.");
	}

	// Another child protection:
	if ((windowRecord->windowType != kPsychDoubleBufferOnscreen) || PsychPrefStateGet_EmulateOldPTB()>0) {
		PsychErrorExitMsg(PsychError_user, "Imaging Pipeline setup: Sorry, imaging pipeline only supported on double buffered onscreen windows and if not in emulation mode for old PTB-2.\n");
	}

	// Specific setup of pipeline if real imaging ops are requested:
	
	// Setup mode switch in record:
	windowRecord->imagingMode = imagingmode;

	// Is this a request for fast Offscreen window support only? The special flag kPsychNeedFastOffscreenWindows,
	// if provided without any other flags, will not startup the full pipeline, but only allow creation of FBO-backed
	// offscreen windows and fast switching between them and the system framebuffer.
	if (imagingmode == kPsychNeedFastOffscreenWindows) {
		if (PsychPrefStateGet_Verbosity()>2) printf("PTB-INFO: Support for fast OffscreenWindows enabled.\n");
		fflush(NULL);
		return;
	}

	if (PsychPrefStateGet_Verbosity()>2) printf("PTB-INFO: Psychtoolbox imaging pipeline starting up for window with requested imagingmode %i ...\n", imagingmode);
	fflush(NULL);

	// Safe default:
	targetisfinalFB = FALSE;

    // Panel fitter requested? If so, framebuffer blit extension supported?
    if ((imagingmode & kPsychNeedGPUPanelFitter) && !(windowRecord->gfxcaps & kPsychGfxCapFBOBlit)) {
        // This is a no-go:
        printf("PTB-WARNING: You requested use of the panel-fitter via the 'clientRect' parameter of Screen('OpenWindow', ...);\n");
        printf("PTB-WARNING: but this graphics card or graphics driver does not support the required GL_EXT_framebuffer_blit extension.\n");
        printf("PTB-WARNING: Will use a lower performance fallback path to do it anyway. Some restrictions apply, check further status output.\n\n");
    }

	// Multisampled anti-aliasing requested?
	if (multiSample > 0) {
		// Yep. Supported by GPU?
		if (!(windowRecord->gfxcaps & kPsychGfxCapFBOMultisample)) {
			// No. We fall back to non-multisampled mode:
			multiSample = 0;
			windowRecord->multiSample = 0;
            
			// Tell user if warnings enabled:
			if (PsychPrefStateGet_Verbosity() > 1) {
				printf("PTB-WARNING: You requested stimulus anti-aliasing by multisampling by setting the multiSample parameter of Screen('OpenWindow', ...) to a non-zero value.\n");
				printf("PTB-WARNING: You also requested use of the imaging pipeline. Unfortunately, your combination of operating system, graphics hardware and driver does not\n");
				printf("PTB-WARNING: support simultaneous use of the imaging pipeline and multisampled anti-aliasing.\n");
				printf("PTB-WARNING: Will therefore continue without anti-aliasing...\n\n");
				printf("PTB-WARNING: A driver upgrade may resolve this issue. Users of MacOS-X need at least OS/X 10.5.2 Leopard for support on recent ATI hardware.\n\n");
			}
		}    // Panel scaling requested? If so we need support for scaled multisample resolve blits to satisfy needs of multisampling and scaling:
        else if ((imagingmode & kPsychNeedGPUPanelFitter) && !(windowRecord->gfxcaps & kPsychGfxCapFBOScaledResolveBlit)) {
            // Not supported by GPU. Disable multisampling to satisfy at least the requirement for panelscaling,
            // which is probably more important, as usercode usually only uses panel scaling to workaround serious
            // trouble with experimental setups, ie., it is more urgent:
			multiSample = 0;
			windowRecord->multiSample = 0;
			
			// Tell user if warnings enabled:
			if (PsychPrefStateGet_Verbosity() > 1) {
				printf("PTB-WARNING: You requested stimulus anti-aliasing by multisampling by setting the multiSample parameter of Screen('OpenWindow', ...) to a non-zero value.\n");
				printf("PTB-WARNING: You also requested use of the imaging pipeline and of the GPU panel-fitter / rescaler via the 'clientRect' argument.\n");
                printf("PTB-WARNING: Unfortunately, your combination of operating system, graphics hardware and driver does not support simultaneous\n");
				printf("PTB-WARNING: use of multisampled anti-aliasing and the panel-fitter. I assume your request for panel fitting is more important.\n");
				printf("PTB-WARNING: Will therefore continue without anti-aliasing to make the panel-fitter work.\n");
                printf("PTB-WARNING: You would need a graphics card, os or graphics driver that supports the GL_EXT_framebuffer_multisample_blit_scaled\n");
                printf("PTB-WARNING: extension to avoid this degradation of functionality.\n\n");
			}
        }
	}

	// Determine required precision for our framebuffer objects:

	// Start off with standard 8 bpc fixed point:
	fboInternalFormat = GL_RGBA8; windowRecord->bpc = 8;
	
	// Need 16 bpc fixed point precision?
	if (imagingmode & kPsychNeed16BPCFixed) { fboInternalFormat = ((windowRecord->gfxcaps & kPsychGfxCapSNTex16) ? GL_RGBA16_SNORM : GL_RGBA16); windowRecord->bpc = 16; }
	 
	// Need 16 bpc floating point precision?
	if (imagingmode & kPsychNeed16BPCFloat) { fboInternalFormat = GL_RGBA_FLOAT16_APPLE; windowRecord->bpc = 16; }
	
	// Need 32 bpc floating point precision?
	if (imagingmode & kPsychNeed32BPCFloat) { fboInternalFormat = GL_RGBA_FLOAT32_APPLE; windowRecord->bpc = 32; }

	// Want dynamic adaption of buffer precision?
	if (imagingmode & kPsychUse32BPCFloatAsap) {
		// Yes. If the gfx-hardware is capable of unrestricted hardware-accelerated 32 bpc
		// float framebuffer blending, we use 32bpc float for drawBufferFBOs, if not -> use 16 bpc.
		
		// Start off with 16 bpc for the stage 0 FBO's.
		fboInternalFormat = GL_RGBA_FLOAT16_APPLE; windowRecord->bpc = 16;

		// Blending on 32 bpc float FBO's supported? Upgrade to 32 bpc float for stage 0 if possible:
		if (windowRecord->gfxcaps & kPsychGfxCapFPBlend32) { fboInternalFormat = GL_RGBA_FLOAT32_APPLE; windowRecord->bpc = 32; }
	}

    // Floating point framebuffer on OpenGL-ES requested?
    if (PsychIsGLES(windowRecord) && (imagingmode & (kPsychNeed16BPCFloat | kPsychNeed32BPCFloat | kPsychUse32BPCFloatAsap))) {
        // Yes. We only support 32 bpc float framebuffers with alpha-blending. On less supportive hardware we fail:
        if (!(windowRecord->gfxcaps & kPsychGfxCapFPTex32) || !(windowRecord->gfxcaps & kPsychGfxCapFPFBO32)) {
            PsychErrorExitMsg(PsychError_user, "Sorry, the requested framebuffer color resolution of 32 bpc floating point is not supported by your graphics card. Game over.");
        }

        // Supported. Upgrade requested format to 32 bpc float, whatever it was before:
        fboInternalFormat = GL_RGBA_FLOAT32_APPLE; windowRecord->bpc = 32;
    }

    // Sanity check: Is a floating point framebuffer requested which requires floating point texture support and
    // the hardware doesn't support float textures? Bail early, if so.
    if (((windowRecord->bpc == 16) && !(windowRecord->gfxcaps & kPsychGfxCapFPTex16) && !(imagingmode & kPsychNeed16BPCFixed)) ||
        ((windowRecord->bpc == 32) && !(windowRecord->gfxcaps & kPsychGfxCapFPTex32))) {
        printf("PTB-ERROR: Your script requested a floating point resolution framebuffer with a resolution of more than 8 bits per color channel.\n");
        printf("PTB-ERROR: Your graphics hardware doesn't support floating point textures or framebuffers, so this is a no-go. Aborting...\n");
		PsychErrorExitMsg(PsychError_user, "Sorry, the requested framebuffer color resolution is not supported by your graphics card. Game over.");
    }

	if (PsychPrefStateGet_Verbosity()>2) {
		switch(fboInternalFormat) {
			case GL_RGBA8:
				printf("PTB-INFO: Will use 8 bits per color component framebuffer for stimulus drawing.\n");
			break;

			case GL_RGBA16:
				printf("PTB-INFO: Will use 16 bits per color component unsigned integer framebuffer for stimulus drawing. Alpha blending may not work.\n");
			break;

			case GL_RGBA16_SNORM:
				printf("PTB-INFO: Will use 15 bits per color component signed integer framebuffer for stimulus drawing. Alpha blending may not work.\n");
			break;

			case GL_RGBA_FLOAT16_APPLE:
				printf("PTB-INFO: Will use 16 bits per color component floating point framebuffer for stimulus drawing. ");
				if (windowRecord->gfxcaps & kPsychGfxCapFPBlend16) {
					printf("Alpha blending should work correctly.\n");
					if (imagingmode & kPsychUse32BPCFloatAsap) {
						printf("PTB-INFO: Can't use 32 bit precision for drawing because hardware doesn't support alpha-blending in 32 bpc.\n");
					}
				}
				else {
					printf("Alpha blending may not work on your system with this setup, but only for 8 bits per color component mode.\n");
				}
			break;

			case GL_RGBA_FLOAT32_APPLE:
				printf("PTB-INFO: Will use 32 bits per color component floating point framebuffer for stimulus drawing. ");
				if (windowRecord->gfxcaps & kPsychGfxCapFPBlend32) {
					printf("Alpha blending should work correctly.\n");
				}
				else {
					printf("Alpha blending may not work on your system with this setup, but only for lower precision modes.\n");
				}
			break;
		}
	}

	// Do we need additional depth buffer attachments?
	needzbuffer = (PsychPrefStateGet_3DGfx()>0) ? TRUE : FALSE;
	
	// Do we need separate streams for stereo? Only for OpenGL quad-buffered mode and dual-window stereo mode:
	needseparatestreams = (windowRecord->stereomode == kPsychOpenGLStereo || windowRecord->stereomode == kPsychDualWindowStereo || windowRecord->stereomode == kPsychFrameSequentialStereo) ? TRUE : FALSE;

	// Do we need some intermediate image processing?
	needimageprocessing= (imagingmode & kPsychNeedImageProcessing) ? TRUE : FALSE;

	// Do we need some final output formatting?
	needoutputconversion = (imagingmode & kPsychNeedOutputConversion) ? TRUE : FALSE;
	
	// Do we need fast backing store?
	needfastbackingstore = (imagingmode & kPsychNeedFastBackingStore) ? TRUE : FALSE;
	
	// Consolidate settings: Most settings imply kPsychNeedFastBackingStore.
	if (needoutputconversion || needimageprocessing || windowRecord->stereomode > 0 || fboInternalFormat!=GL_RGBA8) {
		imagingmode|=kPsychNeedFastBackingStore;
		needfastbackingstore = TRUE;
	}
	
	// Try to allocate and configure proper FBO's:
	fbocount = 0;
	
	// Define final default output buffers as system framebuffers: We create some pseudo-FBO's for these
	// which describe the system framebuffer (backbuffer). This is done to simplify pipeline design:

	// Allocate empty FBO info struct and assign it:
	winwidth=(int)PsychGetWidthFromRect(windowRecord->rect);
	winheight=(int)PsychGetHeightFromRect(windowRecord->rect);

	if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), 0, FALSE, winwidth, winheight, 0, 0)) {
		// Failed!
		PsychErrorExitMsg(PsychError_internal, "Imaging Pipeline setup: Could not setup stage 0 of imaging pipeline.");
	}

	// The pseudo-FBO initially contains a fboid of zero == system framebuffer, and empty (zero) attachments.
	// The up to now only useful information is the viewport geometry ie winwidth and winheight.

	// We use the same struct for both buffers, because in the end, there is only one backbuffer. Separate channels
	// with same mapping allow some interesting extensions in the future for additional stereo modes or snapshot
	// creation...
	windowRecord->finalizedFBO[0]=fbocount;
	windowRecord->finalizedFBO[1]=fbocount;
	fbocount++;

	// Compute format for a finalizedFBO which is not the backbuffer, but blits unmodified to the backbuffer:
	// This is by default always a standard 8bpc fixed point RGBA8 framebuffer without stencil- and z-buffers etc.
	// If bit depths of native backbuffer is more than 8 bits, we allocate a float32 FBO though. Should we need a
	// 32 bpc float FBO but the GPU doesn't support this, we try a 16 bit snorm FBO. A 16 bit snorm FBO has effective
	// 15 bits linear integer precision, which is enough for all currently existing hw framebuffers:
	// Query bit depths of native backbuffer: Assume red bits == green bits == blue bits == bit depths.
	glGetIntegerv(GL_RED_BITS, &redbits);
	// Decide on proper format:
	finalizedFBOFormat = (redbits <= 8) ? GL_RGBA8 : ((windowRecord->gfxcaps & kPsychGfxCapFPFBO32) ? GL_RGBA_FLOAT32_APPLE : GL_RGBA16_SNORM);

	if ((windowRecord->stereomode == kPsychDualWindowStereo) || (imagingmode & kPsychNeedDualWindowOutput)) {
		// Dual-window stereo or dual window output is a special case: This window contains the imaging pipeline for
		// both views, but its OpenGL context and framebuffer only represents the left-view channel.
		// The right-view channel is represented by a slave window and its associated context. We
		// create a framebuffer object and attach it to finalizedFBO[1]. This way, the final left view
		// stimulus image gets blitted directly into the system framebuffer for this window, but the
		// right view stimulus gets blitted into our real finalizedFBO[1]. PsychPreflipOperations() will
		// perform a last blitcopy operation at the end of pipeline processing, where it copies the content
		// of finalizedFBO[1] into the real system framebuffer for the onscreen window which represents the
		// user visible right view.
		//
		// In dual window output mode, we may only have one merged/composited stereo view or even only
		// a single monoscopic view, but we still distribute that view to both finalizedFBO's aka different
		// onscreen windows backbuffers, possibly with separate output formatting / postprocessing.
		if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), finalizedFBOFormat, FALSE, winwidth, winheight, 0, 0)) {
			// Failed!
			PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 0 of imaging pipeline for dual-window stereo.");
		}
		
		windowRecord->finalizedFBO[1]=fbocount;
		fbocount++;
	}

	if (windowRecord->stereomode == kPsychFrameSequentialStereo) {
		// Home-Grown frame-sequential stereo mode: Need one real finalizedFBO for each of the
		// two stereo streams:
		if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), finalizedFBOFormat, FALSE, winwidth, winheight, 0, 0)) {
			// Failed!
			PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 0 of imaging pipeline for frame-sequential stereo (left eye).");
		}
		
		windowRecord->finalizedFBO[0]=fbocount;
		fbocount++;

		if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), finalizedFBOFormat, FALSE, winwidth, winheight, 0, 0)) {
			// Failed!
			PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 0 of imaging pipeline for frame-sequential stereo (right eye).");
		}
		
		windowRecord->finalizedFBO[1]=fbocount;
		fbocount++;
	}

	// Now we preinit all further stages with the finalizedFBO assignment.
	if (needfastbackingstore) {
		// We need at least the 1st level drawBufferFBO's as rendertargets for all
		// user-space drawing, ie Screen 2D drawing functions, MOGL OpenGL rendering and
		// C-MEX OpenGL rendering plugins...
		
		// Define dimensions of 1st stage FBO:
		winwidth=(int)PsychGetWidthFromRect(windowRecord->rect);
		winheight=(int)PsychGetHeightFromRect(windowRecord->rect);

		// Adapt it for some stereo modes:
		if (windowRecord->specialflags & kPsychHalfWidthWindow) {
			// Special case for stereo: Only half the real window width:
			winwidth = winwidth / 2;
		}

		if (windowRecord->specialflags & kPsychTwiceWidthWindow) {
			// Special case: Twice the real window width:
			winwidth = winwidth * 2;
		}

		if (windowRecord->specialflags & kPsychHalfHeightWindow) {
			// Special case for stereo: Only half the real window height:
			winheight = winheight / 2;
		}

        // Special setup of FBO size for use with panel-fitter needed?
        if (imagingmode & kPsychNeedGPUPanelFitter) {
            // Panel-fitter: Use override values from windows clientrectangle:            
            clientwidth  = (int) PsychGetWidthFromRect(windowRecord->clientrect);
            clientheight = (int) PsychGetHeightFromRect(windowRecord->clientrect);
            
            if (PsychPrefStateGet_Verbosity() > 2) {
                printf("PTB-INFO: Enabling panel fitter. Rescaling from %i x %i pixels virtual size to real framebuffer resolution.\n", clientwidth, clientheight);
            }
        }
        else {
            // No panel-fitter, use calculated winwidth x winheight for 1st level drawBufferFBO's:
            clientwidth  = winwidth;
            clientheight = winheight;
        }

		// These FBO's may need a z-buffer or stencil buffer as well if 3D rendering is
		// enabled:
		if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, needzbuffer, clientwidth, clientheight, multiSample, 0)) {
			// Failed!
			PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 1 of imaging pipeline.");
		}
		
		if ((PsychPrefStateGet_Verbosity() > 2) && (windowRecord->fboTable[fbocount]->multisample > 0)) {
			printf("PTB-INFO: Created framebuffer for anti-aliasing with %i samples per pixel for use with imaging pipeline.\n", windowRecord->fboTable[fbocount]->multisample);
		}
		
		// Assign this FBO as drawBuffer for left-eye or mono channel:
		windowRecord->drawBufferFBO[0] = fbocount;
		fbocount++;
		
		// If we are in stereo mode, we'll need a 2nd buffer for the right-eye channel:
		if (windowRecord->stereomode > 0) {
			if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, needzbuffer, clientwidth, clientheight, multiSample, 0)) {
				// Failed!
				PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 1 of imaging pipeline.");
			}
			
			// Assign this FBO as drawBuffer for right-eye channel:
			windowRecord->drawBufferFBO[1] = fbocount;
			fbocount++;
		}
		
		// Windows with fast backing store always have 4 color channels RGBA, regardless what the
		// associated system framebuffer has:
		windowRecord->nrchannels = 4;
	}
	
	// Upgrade to 32 bpc float FBO's needed, starting with the 2nd stage of the pipe?
	// If so, we can now upgrade, because we won't need alpha-blending anymore in later stages:
	if (imagingmode & kPsychUse32BPCFloatAsap) fboInternalFormat = GL_RGBA_FLOAT32_APPLE;
	
	if (PsychPrefStateGet_Verbosity()>2) {
		switch(fboInternalFormat) {
			case GL_RGBA8:
				printf("PTB-INFO: Will use 8 bits per color component framebuffer for stimulus post-processing (if any).\n");
			break;

			case GL_RGBA16:
				printf("PTB-INFO: Will use 16 bits per color component unsigned integer framebuffer for stimulus post-processing (if any).\n");
			break;

			case GL_RGBA16_SNORM:
				printf("PTB-INFO: Will use 15 bits per color component signed integer framebuffer for stimulus post-processing (if any).\n");
			break;

			case GL_RGBA_FLOAT16_APPLE:
				printf("PTB-INFO: Will use 16 bits per color component floating point framebuffer for stimulus post-processing (if any).\n");
			break;

			case GL_RGBA_FLOAT32_APPLE:
				printf("PTB-INFO: Will use 32 bits per color component floating point framebuffer for stimulus post-processing (if any).\n");
			break;
		}
	}

	// Multisampling requested? Or panel-fitter active?
	if ((multiSample > 0) || (imagingmode & kPsychNeedGPUPanelFitter)) {
		// Multisampling or panel-fitting requested. Need to find out if we are an intermediate multisample-resolve / rescaler buffer
		// or if this is already the final destination and we can resolve and/or scale-blit directly into system framebuffer/
		// into finalizedFBO's:
	
		// The target of the drawBufferFBO's is already the final FB if not processing is needed. This is the case
		// if all of the following holds:
		// a) No image processing requested.
		// b) No stereo mode active, therefore no need for any kind of stereo compositing or merging.
		// c) No output conversion / final formatting needed.
		targetisfinalFB = ( !needimageprocessing && (windowRecord->stereomode == kPsychMonoscopic) && !needoutputconversion ) ? TRUE : FALSE;
	
		if (!targetisfinalFB) {
			// Yes. Setup real inputBuffers as multisample-resolve / scaler targets:
			if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, FALSE, winwidth, winheight, 0, 0)) {
				// Failed!
				PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 1 inputBufferFBO of imaging pipeline.");
			}
			
			// Assign this FBO as inputBufferFBO for left-eye or mono channel:
			windowRecord->inputBufferFBO[0] = fbocount;
			fbocount++;
		}
		else {
			// Nothing further to do! Just set us as final framebuffer:
			windowRecord->inputBufferFBO[0] = windowRecord->finalizedFBO[0];
		}
		
		// If we are in stereo mode, we'll need a 2nd buffer for the right-eye channel:
		if (windowRecord->stereomode > 0) {
			if (!targetisfinalFB) {
				if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, FALSE, winwidth, winheight, 0, 0)) {
					// Failed!
					PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 1 inputBufferFBO of imaging pipeline.");
				}
				
				// Assign this FBO as drawBuffer for right-eye channel:
				windowRecord->inputBufferFBO[1] = fbocount;
				fbocount++;
			}
			else {
				// Nothing further to do! Just set us as final framebuffer:
				windowRecord->inputBufferFBO[1] = windowRecord->finalizedFBO[1];
			}
		}
	}
	else {
		// No. Setup pass-through inputBuffers that do nothing: The "pointers" just point to / replicate the
		// assignment of the drawBufferFBOs so this is a zero-copy op:
		windowRecord->inputBufferFBO[0] = windowRecord->drawBufferFBO[0];
		windowRecord->inputBufferFBO[1] = windowRecord->drawBufferFBO[1];
	}
	
	// Do we need 2nd stage FBOs? We need them as targets for the processed data if support for misc image processing ops is requested.
	if (needimageprocessing) {
		// Need real FBO's as targets for image processing:

		// Define dimensions of 2nd stage FBO:
		winwidth=(int)PsychGetWidthFromRect(windowRecord->rect);
		winheight=(int)PsychGetHeightFromRect(windowRecord->rect);

		// Adapt it for some stereo modes:
		if (windowRecord->specialflags & kPsychHalfWidthWindow) {
			// Special case for stereo: Only half the real window width:
			winwidth = winwidth / 2;
		}

		if (windowRecord->specialflags & kPsychTwiceWidthWindow) {
			// Special case: Twice the real window width:
			winwidth = winwidth * 2;
		}

		if (windowRecord->specialflags & kPsychHalfHeightWindow) {
			// Special case for stereo: Only half the real window height:
			winheight = winheight / 2;
		}

		// Is the target of imageprocessing (our processedDrawBufferFBO) the final destination? This is true if there is no further need
		// for more rendering passes in later processing stages and we don't do any processing here with more than 2 rendering passes. In that
		// case we can directly output to the finalizedFBO without need for one more intermediate buffer -- faster!
		// In all other cases, we'll need an additional buffer. We also exclude quad-buffered stereo, because the image processing blit chains
		// cannot switch between left- and right backbuffer of the system framebuffer...
		targetisfinalFB = ( !(imagingmode & kPsychNeedMultiPass) && (windowRecord->stereomode == kPsychMonoscopic) && !needoutputconversion ) ? TRUE : FALSE;

		if (!targetisfinalFB) {
			// These FBO's don't need z- or stencil buffers anymore:
			if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, FALSE, winwidth, winheight, 0, 0)) {
				// Failed!
				PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 2 of imaging pipeline.");
			}			

			// Assign this FBO as processedDrawBuffer for left-eye or mono channel:
			windowRecord->processedDrawBufferFBO[0] = fbocount;
			fbocount++;			
		}
		else {
			// Can assign final destination:
			windowRecord->processedDrawBufferFBO[0] = windowRecord->finalizedFBO[0];
		}
				
		// If we are in stereo mode, we'll need a 2nd buffer for the right-eye channel:
		if (windowRecord->stereomode > 0) {
			if (!targetisfinalFB) {
				// These FBO's don't need z- or stencil buffers anymore:
				if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, FALSE, winwidth, winheight, 0, 0)) {
					// Failed!
					PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 2 of imaging pipeline.");
				}			
				
				// Assign this FBO as processedDrawBuffer for right-eye channel:
				windowRecord->processedDrawBufferFBO[1] = fbocount;
				fbocount++;
			}
			else {
				// Can assign final destination:
				windowRecord->processedDrawBufferFBO[1] = windowRecord->finalizedFBO[1];
			}
		}
		else {
			// Mono mode: No right-eye buffer:
			windowRecord->processedDrawBufferFBO[1] = -1;
		}
		
		// Allocate a bounce-buffer as well if multi-pass rendering is requested:
		if (imagingmode & kPsychNeedDualPass || imagingmode & kPsychNeedMultiPass) {
			if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, FALSE, winwidth, winheight, 0, 0)) {
				// Failed!
				PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 2 of imaging pipeline.");
			}
			
			// Assign this FBO as processedDrawBuffer for bounce buffer ops in multi-pass rendering:
			windowRecord->processedDrawBufferFBO[2] = fbocount;
			fbocount++;
		}
		else {
			// No need for bounce-buffers, only single-pass processing requested.
			windowRecord->processedDrawBufferFBO[2] = -1;
		}
	}
	else {
		// No image processing: Set 2nd stage FBO's to 1st stage FBO's:
		windowRecord->processedDrawBufferFBO[0] = windowRecord->inputBufferFBO[0];
		windowRecord->processedDrawBufferFBO[1] = windowRecord->inputBufferFBO[1];
		windowRecord->processedDrawBufferFBO[2] = -1;
	}
	
	// Stage 2 ready. Any need for real merged FBO's? We need a merged FBO if we are in stereo mode
	// and in need to merge output from the two views and to postprocess that output. In all other
	// cases there's no need for real merged FBO's and we do just a "pass-through" assignment.
	if ((windowRecord->stereomode > 0) && (!needseparatestreams) && (needoutputconversion)) {
		// Need real FBO's as targets for merger output.

		// Define dimensions of 3rd stage FBO:
		winwidth=(int)PsychGetWidthFromRect(windowRecord->rect);
		winheight=(int)PsychGetHeightFromRect(windowRecord->rect);

        // Must not take half width and half height flags into account,
        // but need to make sure we retain info in a double-width buffer
        // until we reach the final system framebuffer or final output fbo:
		if (windowRecord->specialflags & kPsychTwiceWidthWindow) {
			// Special case: Twice the real window width:
			winwidth = winwidth * 2;
		}

		// These FBO's don't need z- or stencil buffers anymore:
		if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, FALSE, winwidth, winheight, 0, 0)) {
			// Failed!
			PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 3 of imaging pipeline.");
		}
		
		// Assign this FBO for left-eye and right-eye channel: The FBO is shared accross channels...
		windowRecord->preConversionFBO[0] = fbocount;
		windowRecord->preConversionFBO[1] = fbocount;
		fbocount++;

		// Request bounce buffer:
		windowRecord->preConversionFBO[2] = -1000;
	}
	else {
		if ((windowRecord->stereomode > 0) && (!needseparatestreams)) {
			// Need to merge two streams, but don't need any output conversion on them. We
			// don't need an extra FBO for the merge results! We just write our merge results
			// into whatever the final framebuffer is - Could be another FBO if framebuffer
			// snapshotting is requested, but most likely its the pseudo-FBO of the system
			// backbuffer. Anyway, the proper ones are stored in finalizedFBO[]:
			windowRecord->preConversionFBO[0] = windowRecord->finalizedFBO[0];
			windowRecord->preConversionFBO[1] = windowRecord->finalizedFBO[1];
			// Request bounce buffer:
			windowRecord->preConversionFBO[2] = -1000;
		}
		else {		
			// No merge operation needed. Do we need output conversion?
			if (needoutputconversion) {
				// Output conversion needed. Set input for this stage to output of the
				// image processing.
				windowRecord->preConversionFBO[0] = windowRecord->processedDrawBufferFBO[0];
				windowRecord->preConversionFBO[1] = windowRecord->processedDrawBufferFBO[1];
				// Request bounce buffer:
				windowRecord->preConversionFBO[2] = -1000;
			}
			else {
				// No merge and no output conversion needed. In that case, PsychPreFlipOperations()
				// will behave as if output conversion is requested, but with the identity blit chain,
				// and merge stage is skipped, so we need to set the preConversionFBO's as if conversion
				// is done.
				windowRecord->preConversionFBO[0] = windowRecord->processedDrawBufferFBO[0];
				windowRecord->preConversionFBO[1] = windowRecord->processedDrawBufferFBO[1];
				// No bounce buffer needed:
				windowRecord->preConversionFBO[2] = -1;
			}
			
		}		
	}

	// Do we need a bounce buffer for merging and/or conversion?
	if (windowRecord->preConversionFBO[2] == -1000) {
		// Yes. We can reuse/share the bounce buffer of the image processing stage if
		// one exists and is of suitable size i.e. we're not in dual-view stereo - in
		// that case all buffers are of same size.
		if ((windowRecord->processedDrawBufferFBO[2]!=-1) && 
			!(windowRecord->stereomode==kPsychFreeFusionStereo || windowRecord->stereomode==kPsychFreeCrossFusionStereo)) {
			// Stage 1 bounce buffer is suitable for sharing, assign it:
			windowRecord->preConversionFBO[2] = windowRecord->processedDrawBufferFBO[2];
		}
		else {
			// We need a new, private bounce-buffer:
			if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, FALSE, winwidth, winheight, 0, 0)) {
				// Failed!
				PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 3 of imaging pipeline [1st bounce buffer].");
			}
			
			windowRecord->preConversionFBO[2] = fbocount;
			fbocount++;				
		}
		
		// In any case, we need a new private 2nd bounce buffer for the special case of the final processing chain:
		if (!PsychCreateFBO(&(windowRecord->fboTable[fbocount]), fboInternalFormat, FALSE, winwidth, winheight, 0, 0)) {
			// Failed!
			PsychErrorExitMsg(PsychError_system, "Imaging Pipeline setup: Could not setup stage 3 of imaging pipeline [2nd bounce buffer].");
		}
		
		windowRecord->preConversionFBO[3] = fbocount;
		fbocount++;
	}

	// If dualwindow output is requested and preConversionFBO[1] isn't assigned yet,
	// we set it to the same FBO as preConversionFBO[0], so the image data of our one
	// single image buffer is distributed to both output pipes for output conversion and
	// display:
	if ((imagingmode & kPsychNeedDualWindowOutput) && (windowRecord->preConversionFBO[1] == -1)) {
		windowRecord->preConversionFBO[1] = windowRecord->preConversionFBO[0];
	}

	// Setup imaging mode flags:
	newimagingmode = (needseparatestreams) ? kPsychNeedSeparateStreams : 0;
	if (!needseparatestreams && (windowRecord->stereomode > 0)) newimagingmode |= kPsychNeedStereoMergeOp;
	if (needfastbackingstore) newimagingmode |= kPsychNeedFastBackingStore;
	if (needoutputconversion) newimagingmode |= kPsychNeedOutputConversion;
	if (needimageprocessing)  newimagingmode |= kPsychNeedImageProcessing;
	if (imagingmode & kPsychNeed32BPCFloat) {
		newimagingmode |= kPsychNeed32BPCFloat;
	}
	else if (imagingmode & kPsychNeed16BPCFloat) {
		newimagingmode |= kPsychNeed16BPCFloat;
	}
	else if (imagingmode & kPsychNeed16BPCFixed) {
		newimagingmode |= kPsychNeed16BPCFixed;
	}
	if (imagingmode & kPsychNeedDualWindowOutput) newimagingmode |= kPsychNeedDualWindowOutput;
    if (imagingmode & kPsychNeedGPUPanelFitter) newimagingmode |= kPsychNeedGPUPanelFitter;
	
	// Set new final imaging mode and fbocount:
	windowRecord->imagingMode = newimagingmode;
	windowRecord->fboCount = fbocount;

	// The pipelines buffers and information flow are configured now...
	if (PsychPrefStateGet_Verbosity()>4) {
		printf("PTB-DEBUG: Buffer mappings follow...\n");
		printf("fboCount = %i\n", windowRecord->fboCount);
		printf("finalizedFBO = %i, %i\n", windowRecord->finalizedFBO[0], windowRecord->finalizedFBO[1]);
		printf("preConversionFBO = %i, %i, %i, %i\n", windowRecord->preConversionFBO[0], windowRecord->preConversionFBO[1], windowRecord->preConversionFBO[2], windowRecord->preConversionFBO[3]);
		printf("processedDrawBufferFBO = %i %i %i\n", windowRecord->processedDrawBufferFBO[0], windowRecord->processedDrawBufferFBO[1], windowRecord->processedDrawBufferFBO[2]);
		printf("inputBufferFBO = %i %i \n", windowRecord->inputBufferFBO[0], windowRecord->inputBufferFBO[1]);
		printf("drawBufferFBO = %i %i \n", windowRecord->drawBufferFBO[0], windowRecord->drawBufferFBO[1]);
		printf("-------------------------------------\n\n");
		fflush(NULL);
	}

	// Setup our default chain: This chain is executed if some stage of the imaging pipe is set up according
	// to imagingMode, but the corresponding hook-chain is empty or disabled. In that case we need to copy
	// the data for that stage from its input buffers to its output buffers via a simple blit operation.
	PsychPipelineAddBuiltinFunctionToHook(windowRecord, "IdentityBlitChain", "Builtin:IdentityBlit", INT_MAX, "");
	PsychPipelineEnableHook(windowRecord, "IdentityBlitChain");

	// Setup of GLSL stereo shaders for stereo modes that need some merging operations:
	// Quad-buffered stereo and mono mode don't need these...
	if (windowRecord->stereomode > kPsychOpenGLStereo) {
		// Merged stereo mode requested.
		glsl = 0;
		
		// Which mode?
		switch(windowRecord->stereomode) {
			// Anaglyph mode?
			case kPsychAnaglyphRGStereo:
			case kPsychAnaglyphGRStereo:
			case kPsychAnaglyphRBStereo:
			case kPsychAnaglyphBRStereo:
				// These share all code...
				
				// Create anaglyph shader and set proper defaults: These can be changed from the M-File if wanted.
				if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: Creating internal anaglyph stereo compositing shader...\n");
				glsl = PsychCreateGLSLProgram(anaglyphshadersrc, NULL, NULL);
				if (glsl) {
					// Bind it:
					glUseProgram(glsl);

					// Set channel to texture units assignments:
					glUniform1i(glGetUniformLocation(glsl, "Image1"), 0);
					glUniform1i(glGetUniformLocation(glsl, "Image2"), 1);
					
					// Set per-channel color gains: 0 masks the channel, >0 enables it. The values can be used
					// to compensate for differences in the color reproduction of different monitors to reduce
					// cross-talk / ghosting:
					
					// Left-eye channel (channel 1):
					rg = (float) ((windowRecord->stereomode==kPsychAnaglyphRGStereo || windowRecord->stereomode==kPsychAnaglyphRBStereo) ? 1.0 : 0.0);
					gg = (float) ((windowRecord->stereomode==kPsychAnaglyphGRStereo) ? 1.0 : 0.0);
					bg = (float) ((windowRecord->stereomode==kPsychAnaglyphBRStereo) ? 1.0 : 0.0);
					glUniform3f(glGetUniformLocation(glsl, "Gains1"), rg, gg, bg);

					// Right-eye channel (channel 2):
					rg = (float) ((windowRecord->stereomode==kPsychAnaglyphGRStereo || windowRecord->stereomode==kPsychAnaglyphBRStereo) ? 1.0 : 0.0);
					gg = (float) ((windowRecord->stereomode==kPsychAnaglyphRGStereo) ? 1.0 : 0.0);
					bg = (float) ((windowRecord->stereomode==kPsychAnaglyphRBStereo) ? 1.0 : 0.0);
					glUniform3f(glGetUniformLocation(glsl, "Gains2"), rg, gg, bg);
					
					// Define default weights for RGB -> Luminance conversion: We default to the standardized NTSC color weights.
					glUniform3f(glGetUniformLocation(glsl, "ColorToGrayWeights"), 0.299f, 0.587f, 0.114f);
					// Define background bias color to add: Normally zero for standard anaglyph:
					glUniform3f(glGetUniformLocation(glsl, "ChannelBias"), 0.0, 0.0, 0.0);

					// Unbind it, its ready!
					glUseProgram(0);

					if (glsl) {
						// Add shader to processing chain:
						PsychPipelineAddShaderToHook(windowRecord, "StereoCompositingBlit", "StereoCompositingShaderAnaglyph", INT_MAX, glsl, "", 0) ;
						
						// Enable stereo compositor:
						PsychPipelineEnableHook(windowRecord, "StereoCompositingBlit");		
					}
				}
				else {
					PsychErrorExitMsg(PsychError_user, "PTB-ERROR: Failed to create anaglyph stereo processing shader -- Anaglyph stereo won't work!\n");
				}
			break;
			
			case kPsychFreeFusionStereo:
			case kPsychFreeCrossFusionStereo:
				if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: Creating internal dualview stereo compositing shader...\n");
				
				glsl = PsychCreateGLSLProgram(passthroughshadersrc, NULL, NULL);
				if (glsl) {
					// Bind it:
					glUseProgram(glsl);
					// Set channel to texture units assignments:
					glUniform1i(glGetUniformLocation(glsl, "Image1"), (windowRecord->stereomode == kPsychFreeFusionStereo) ? 0 : 1);
					glUseProgram(0);
					
					// Add shader to processing chain:
					sprintf(blittercfg, "Builtin:IdentityBlit:Offset:%i:%i", 0, 0);
					PsychPipelineAddShaderToHook(windowRecord, "StereoCompositingBlit", "StereoCompositingShaderDualViewLeft", INT_MAX, glsl, blittercfg, 0);
				}
				else {
					PsychErrorExitMsg(PsychError_user, "PTB-ERROR: Failed to create left channel dualview stereo processing shader -- Dualview stereo won't work!\n");
				}
				
				glsl = PsychCreateGLSLProgram(passthroughshadersrc, NULL, NULL);
				if (glsl) {
					// Bind it:
					glUseProgram(glsl);
					// Set channel to texture units assignments:
					glUniform1i(glGetUniformLocation(glsl, "Image1"),  (windowRecord->stereomode == kPsychFreeFusionStereo) ? 1 : 0);
					glUseProgram(0);
					
					// Add shader to processing chain:
					sprintf(blittercfg, "Builtin:IdentityBlit:Offset:%i:%i", (int) PsychGetWidthFromRect(windowRecord->rect)/2, 0);
					PsychPipelineAddShaderToHook(windowRecord, "StereoCompositingBlit", "StereoCompositingShaderDualViewRight", INT_MAX, glsl, blittercfg, 0);
				}
				else {
					PsychErrorExitMsg(PsychError_user, "PTB-ERROR: Failed to create right channel dualview stereo processing shader -- Dualview stereo won't work!\n");
				}

				// Enable stereo compositor:
				PsychPipelineEnableHook(windowRecord, "StereoCompositingBlit");		
			break;
			
			case kPsychCompressedTLBRStereo:
			case kPsychCompressedTRBLStereo:
				if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: Creating internal vertical split stereo compositing shader...\n");
				
				glsl = PsychCreateGLSLProgram(passthroughshadersrc, NULL, NULL);
				if (glsl) {
					// Bind it:
					glUseProgram(glsl);
					// Set channel to texture units assignments:
					glUniform1i(glGetUniformLocation(glsl, "Image1"), (windowRecord->stereomode == kPsychCompressedTLBRStereo) ? 0 : 1);
					glUseProgram(0);
					
					// Add shader to processing chain:
					sprintf(blittercfg, "Builtin:IdentityBlit:Offset:%i:%i:Scaling:%f:%f", 0, 0, 1.0, 0.5);
					PsychPipelineAddShaderToHook(windowRecord, "StereoCompositingBlit", "StereoCompositingShaderCompressedTop", INT_MAX, glsl, blittercfg, 0);
				}
				else {
					PsychErrorExitMsg(PsychError_user, "PTB-ERROR: Failed to create left channel dualview stereo processing shader -- Dualview stereo won't work!\n");
				}
				
				glsl = PsychCreateGLSLProgram(passthroughshadersrc, NULL, NULL);
				if (glsl) {
					// Bind it:
					glUseProgram(glsl);
					// Set channel to texture units assignments:
					glUniform1i(glGetUniformLocation(glsl, "Image1"),  (windowRecord->stereomode == kPsychCompressedTLBRStereo) ? 1 : 0);
					glUseProgram(0);
					
					// Add shader to processing chain:
					sprintf(blittercfg, "Builtin:IdentityBlit:Offset:%i:%i:Scaling:%f:%f", 0, (int) PsychGetHeightFromRect(windowRecord->rect)/2, 1.0, 0.5);
					PsychPipelineAddShaderToHook(windowRecord, "StereoCompositingBlit", "StereoCompositingShaderCompressedBottom", INT_MAX, glsl, blittercfg, 0);
				}
				else {
					PsychErrorExitMsg(PsychError_user, "PTB-ERROR: Failed to create right channel dualview stereo processing shader -- Dualview stereo won't work!\n");
				}

				// Enable stereo compositor:
				PsychPipelineEnableHook(windowRecord, "StereoCompositingBlit");
			break;

			case kPsychOpenGLStereo:
			case kPsychFrameSequentialStereo:
				// Nothing to do for now: Setup of blue-line syncing is done in SCREENOpenWindow.c, because it also
				// applies to non-imaging mode...
			break;
			
			case kPsychDualWindowStereo:
				// Nothing to do for now.
			break;

			default:
				PsychErrorExitMsg(PsychError_internal, "Unknown stereo mode encountered! FIX SCREENOpenWindow.c to catch this at the appropriate place!\n");
		}
	}
	

	//PsychPipelineAddBuiltinFunctionToHook(windowRecord, "FinalOutputFormattingBlit", "Builtin:IdentityBlit", INT_MAX, "");
	//PsychPipelineEnableHook(windowRecord, "FinalOutputFormattingBlit");

	//PsychPipelineAddBuiltinFunctionToHook(windowRecord, "StereoLeftCompositingBlit", "Builtin:IdentityBlit", INT_MAX, "");
	//PsychPipelineAddBuiltinFunctionToHook(windowRecord, "StereoLeftCompositingBlit", "Builtin:FlipFBOs", INT_MAX, "");
	//PsychPipelineAddBuiltinFunctionToHook(windowRecord, "StereoLeftCompositingBlit", "Builtin:IdentityBlit", INT_MAX, "");
	//PsychPipelineEnableHook(windowRecord, "StereoLeftCompositingBlit");

	// Multisampling requested? If so, we need to enable it:
	if (multiSample > 0) {
		glEnable(GL_MULTISAMPLE);
        windowRecord->multiSample = multiSample;
	}

	// Perform a safe reset of current drawing target. This is a warm-start of PTB's drawing
	// engine, so the next drawing command will trigger binding the proper FBO of our pipeline.
	// Before this point (==OpenWindow time), all drawing was directly directed to the system
	// framebuffer - important for all the timing tests and calibrations to work correctly.
	PsychSetDrawingTarget((PsychWindowRecordType*) 0x1);

	// Well done.
	return;
}

/* PsychCreateGLSLProgram()
 *  Try to create GLSL shader from source strings and return handle to new shader.
 *  Returns the shader handle if it worked, 0 otherwise.
 *
 *  fragmentsrc - Source string for fragment shader. NULL if none needed.
 *  vertexsrc   - Source string for vertex shader. NULL if none needed.
 *  primitivesrc - Source string for primitive shader. NULL if none needed.
 *
 */
GLuint PsychCreateGLSLProgram(const char* fragmentsrc, const char* vertexsrc, const char* primitivesrc)
{
	GLuint glsl = 0;
	GLuint shader;
	GLint status;
	char errtxt[10000];
	
	// Reset error state:
	while (glGetError());
	
	// Supported at all on this hardware?
	if (!glewIsSupported("GL_ARB_shader_objects") || !glewIsSupported("GL_ARB_shading_language_100")) {
		printf("PTB-ERROR: Your graphics hardware does not support GLSL fragment shaders! Use of imaging pipeline with current settings impossible!\n");
		return(0);
	}
	
	// Create GLSL program object:
	glsl = glCreateProgram();
	
	// Fragment shader wanted?
	if (fragmentsrc) {
		if (PsychPrefStateGet_Verbosity()>4)  printf("PTB-INFO: Creating the following fragment shader, GLSL source code follows:\n\n%s\n\n", fragmentsrc);

		// Supported on this hardware?
		if (!glewIsSupported("GL_ARB_fragment_shader")) {
			printf("PTB-ERROR: Your graphics hardware does not support GLSL fragment shaders! Use of imaging pipeline with current settings impossible!\n");
			return(0);
		}
		
		// Create shader object:
		shader = glCreateShader(GL_FRAGMENT_SHADER);
		// Feed it with GLSL source code:
		glShaderSource(shader, 1, (const char**) &fragmentsrc, NULL);
		
		// Compile shader:
		glCompileShader(shader);
		
		glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
		if (status != GL_TRUE) {
			printf("PTB-ERROR: Shader compilation for builtin fragment shader failed:\n");
			glGetShaderInfoLog(shader, 9999, NULL, (GLchar*) &errtxt);
			printf("%s\n\n", errtxt);
			glDeleteShader(shader);
			glDeleteProgram(glsl);
			// Failed!
			while (glGetError());

			return(0);
		}
		
		// Attach it to program object:
		glAttachShader(glsl, shader);
	}

	// Vertex shader wanted?
	if (vertexsrc) {
		if (PsychPrefStateGet_Verbosity()>4)  printf("PTB-INFO: Creating the following vertex shader, GLSL source code follows:\n\n%s\n\n", vertexsrc);

		// Supported on this hardware?
		if (!glewIsSupported("GL_ARB_vertex_shader")) {
			printf("PTB-ERROR: Your graphics hardware does not support GLSL vertex shaders! Use of imaging pipeline with current settings impossible!\n");
			return(0);
		}

		// Create shader object:
		shader = glCreateShader(GL_VERTEX_SHADER);
		
		// Feed it with GLSL source code:
		glShaderSource(shader, 1, (const char**) &vertexsrc, NULL);
		
		// Compile shader:
		glCompileShader(shader);
		
		glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
		if (status != GL_TRUE) {
			printf("PTB-ERROR: Shader compilation for builtin vertex shader failed:\n");
			glGetShaderInfoLog(shader, 9999, NULL, (GLchar*) &errtxt);
			printf("%s\n\n", errtxt);
			glDeleteShader(shader);
			glDeleteProgram(glsl);
			// Failed!
			while (glGetError());
			return(0);
		}
		
		// Attach it to program object:
		glAttachShader(glsl, shader);
	}

	// Link into final program object:
	glLinkProgram(glsl);

	// Check link status:
	glGetProgramiv(glsl, GL_LINK_STATUS, &status);
	if (status != GL_TRUE) {
		printf("PTB-ERROR: Shader link operation for builtin glsl program failed:\n");
		glGetProgramInfoLog(glsl, 9999, NULL, (GLchar*) &errtxt);
		printf("Error output follows:\n\n%s\n\n", errtxt);
		glDeleteProgram(glsl);
		// Failed!
		while (glGetError());

		return(0);
	}
	
	while (glGetError());

	// Return new GLSL program object handle:
	return(glsl);
}

/* PsychCreateFBO()
 * Create OpenGL framebuffer object for internal rendering, setup PTB info struct for it.
 * This function creates an OpenGL framebuffer object, creates and attaches a texture of suitable size
 * (width x height) pixels with format fboInternalFormat as color buffer (color attachment 0). Optionally,
 * (if needzbuffer is true) it also creates and attaches suitable z-buffer and stencil-buffer attachments.
 * It checks for correct setup and then stores all relevant information in the PsychFBO struct, pointed by
 * fbo. On success it returns true, on failure it returns false.
 */
psych_bool PsychCreateFBO(PsychFBO** fbo, GLenum fboInternalFormat, psych_bool needzbuffer, int width, int height, int multisample, int specialFlags)
{
    GLenum texturetarget;
	GLenum fborc;
	GLint bpc, maxSamples;
	GLboolean isFloatBuffer;
	char fbodiag[1024];
    GLenum glerr;
    int twidth, theight;
    psych_bool multisampled_cb = FALSE;

	// Eat all GL errors:
	PsychTestForGLErrors();
	
    // Select type of texture target: Always GL_TEXTURE_2D on OpenGL-ES. Also, no multisample anti-aliasing
    // on ES:
    if (PsychIsGLES(NULL)) {
        specialFlags |= 1;
        if ((multisample > 0) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: Requested multisample anti-aliasing unsupported on OpenGL-ES hardware. Disabling.\n");
        multisample = 0;

        // Eat all errors:
        while(glGetError());
    }

    texturetarget = (specialFlags & 0x1) ? GL_TEXTURE_2D : GL_TEXTURE_RECTANGLE_EXT;

	// If fboInternalFormat!=1 then we need to allocate and assign a proper PsychFBO struct first:
	if (fboInternalFormat!=1) {
		*fbo = (PsychFBO*) malloc(sizeof(PsychFBO));
		if (*fbo == NULL) PsychErrorExitMsg(PsychError_outofMemory, "Out of system memory when trying to allocate PsychFBO struct!");

		// Start cleanly for error handling:
		(*fbo)->fboid = 0;
		(*fbo)->stexid = 0;
		(*fbo)->ztexid = 0;
		
		(*fbo)->width = width;
		(*fbo)->height = height;
		(*fbo)->multisample = multisample;
        (*fbo)->textarget = texturetarget;
		
		// fboInternalFormat == 0 --> Only allocate and assign, don't initialize FBO.
		if (fboInternalFormat==0) return(TRUE);
	}

	// Initialize target storage size to use size "as is", as a safe default:
	twidth = width;
	theight = height;

	// Standard path w/o multisampling? --> Setup and or bind texture as color attachment:
	if (multisample <= 0) {
		// Is there already a texture object defined for the color attachment?
		if (fboInternalFormat != (GLenum) 1) {
			// No, need to create one:
			
			// Build color buffer: Create a new texture handle for the color buffer attachment.
			glGenTextures(1, (GLuint*) &((*fbo)->coltexid));
			
			// Bind it:
			glBindTexture(texturetarget, (*fbo)->coltexid);
			
			// Create proper texture: Just allocate proper format, don't assign data.
            if (PsychIsGLES(NULL)) {
                // OpenGL-ES: Internal and external format must match:
                if (fboInternalFormat == GL_RGBA8) fboInternalFormat = GL_RGBA;
                if (fboInternalFormat == GL_RGB8) fboInternalFormat = GL_RGB;

                // Non-power-of-two textures supported on OES?
                if (!strstr(glGetString(GL_EXTENSIONS), "GL_OES_texture_npot")) {
                    // No: Find size of suitable pot texture (smallest power of two which is
                    // greater than or equal to the image size:
                    twidth=1;
                    while (twidth < width) twidth*=2;
                    theight=1;
                    while (theight < height) theight*=2;
                }
                else {
                    // Yes: Use size "as is":
                    twidth = width;
                    theight = height;
                }

                if (fboInternalFormat == GL_RGBA_FLOAT32_APPLE || fboInternalFormat == GL_RGB_FLOAT32_APPLE) {
                    // 32-bpc float path:
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In PsychCreateFBO: OpenGL-ES allocated 32-Bit float %s texture.\n",
                                                                  (fboInternalFormat == GL_RGBA_FLOAT32_APPLE) ? "RGBA32F" : "RGB32F");
                    glTexImage2D(texturetarget, 0, fboInternalFormat, twidth, theight, 0, (fboInternalFormat == GL_RGBA_FLOAT32_APPLE) ? GL_RGBA : GL_RGB, GL_FLOAT, NULL);
                }
                else {
                    // Classic 8-Bit integer path:
                    glTexImage2D(texturetarget, 0, fboInternalFormat, twidth, theight, 0, fboInternalFormat, GL_UNSIGNED_BYTE, NULL);
                    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: In PsychCreateFBO: OpenGL-ES allocated 8-Bit integer texture.\n");
                }
            }
            else {
                // Desktop OpenGL: Use size "as is":
                twidth = width;
                theight = height;
                glTexImage2D(texturetarget, 0, fboInternalFormat, twidth, theight, 0, GL_RGBA, GL_FLOAT, NULL);
            }
		}
		else {
			// Yes. Bind it as texture:
			glBindTexture(texturetarget, (*fbo)->coltexid);
		}
		
		if (glGetError()!=GL_NO_ERROR) {
			printf("PTB-ERROR: Failed to setup internal framebuffer objects color buffer attachment for imaging pipeline!\n");
			printf("PTB-ERROR: Most likely the requested size or colordepth of the window or texture is not supported by your graphics hardware.\n");
			return(FALSE);
		}
		
		// Setup texture wrapping behaviour to clamp, as other behaviours are
		// unsupported on many gfx-cards for rectangle textures:
		glTexParameterf(texturetarget,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
		glTexParameterf(texturetarget,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
		
		// Setup filtering for the textures - Use nearest neighbour by default, as floating
		// point filtering usually unsupported.
		glTexParameterf(texturetarget, GL_TEXTURE_MAG_FILTER,GL_NEAREST);
		glTexParameterf(texturetarget, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
		
		// Texture ready, unbind it.
		glBindTexture(texturetarget, 0);
        
        // Mark use of standard path:
        multisampled_cb = FALSE;
	}
	else {
		// Multisampled FBO: Setup a multisampled renderbuffer as color attachment;
		if (fboInternalFormat == (GLenum) 1) PsychErrorExitMsg(PsychError_internal, "Tried to setup a multisampled FBO, but fboInternalFormat was == 1. PTB implementation bug!");

		glGenRenderbuffersEXT(1, (GLuint*) &((*fbo)->coltexid));
		glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, (*fbo)->coltexid);

        // Query maximum supported number of samples for multi-sampling:
        glGetIntegerv(GL_MAX_SAMPLES_EXT, &maxSamples);
        
        // Clamp multisample level to maximum, warn user if she aimed to high:
        if (multisample > maxSamples) {
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-WARNING: Your graphics hardware does not support anti-aliasing with the requested minimum number of %i samples.\n", multisample);
                printf("PTB-WARNING: Will try requesting the theoretical maximum supported value of %i samples instead and see what i can get.\n", maxSamples);
            }
            
            multisample = maxSamples;
        }
        
        // Try creating multisampled renderbuffers. On failure, decrease requested sample count.
        // Hardware may deny requests if insufficient memory is available.
		do {
            if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychCreateFBO(): Trying to alloc multisample renderbuffer with %i samples.\n", multisample);
			glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, multisample--, fboInternalFormat, width, height);
		} while (((glerr = glGetError()) != GL_NO_ERROR) && (multisample >= 0));
		
        // Worked? Worst case we should have gotten at least a renderbuffer with multisample == 0,
        // ie., a non-multisampled buffer. If even that failed, then something's screwed, e.g.,
        // totally out of memory and we have to give up:
		if (glerr != GL_NO_ERROR) {
            if (PsychPrefStateGet_Verbosity() > 0) {
                printf("PTB-ERROR: Failed to setup internal framebuffer objects color buffer attachment as a multisampled renderbuffer for imaging pipeline!\n");
                if (glerr == GL_OUT_OF_MEMORY) {
                    printf("PTB-ERROR: Reason seems to be an out of graphics memory condition.\n");
                }
                else {
                    printf("PTB-ERROR: Reason for failure is unknown [OpenGL error: %s]\n", gluErrorString(glerr));
                }
            }
            
            glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
			return(FALSE);
		}

		// Got some renderbuffer. Query real number of samples for renderbuffer:
		glGetRenderbufferParameterivEXT(GL_RENDERBUFFER_EXT, GL_RENDERBUFFER_SAMPLES_EXT, &multisample);
        
		// Unbind, we're done with setup:
		glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);

        if ((multisample < (*fbo)->multisample) && (PsychPrefStateGet_Verbosity() > 1)) {
            printf("PTB-WARNING: Could only get %i samples instead of requested %i samples for multi-sampling from hardware.\n", multisample, (*fbo)->multisample);
        }
        
        // Assign effective multisample count:
		(*fbo)->multisample = multisample;

		if (PsychPrefStateGet_Verbosity() > 5) {
			printf("PTB-DEBUG: Created framebuffer object for anti-aliasing with %i samples per pixel for use with imaging pipeline.\n", (*fbo)->multisample);
		}
        
        // Mark use of multi-sampled renderbuffer path:
        multisampled_cb = TRUE;
	}
	
	// Create a new framebuffer object and bind it:
	glGenFramebuffersEXT(1, (GLuint*) &((*fbo)->fboid));
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, (*fbo)->fboid);
	
	if (!multisampled_cb) {
		// Attach the texture as color buffer zero:
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, texturetarget, (*fbo)->coltexid, 0);
	}
	else {
		// Attach the multi-sampled renderbuffer as color buffer zero:
		glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, (*fbo)->coltexid);
	}

	if (glGetError()!=GL_NO_ERROR) {
		printf("PTB-ERROR: Failed to attach internal framebuffer objects color buffer attachment for imaging pipeline!\n");
		printf("PTB-ERROR: Maybe the requested size & depth of the window or texture is not supported by your graphics hardware.\n");
		return(FALSE);
	}

	// Check for framebuffer completeness:
	fborc = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if (fborc!=GL_FRAMEBUFFER_COMPLETE_EXT && fborc!=0) {
		// Framebuffer incomplete!
		while(glGetError());
		printf("PTB-ERROR[Imaging pipeline]: Failed to setup color buffer attachment of internal FBO when trying to prepare drawing into a texture or window.\n");
		if (fborc==GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT) {
			printf("PTB-ERROR: Your graphics hardware does not support the provided or requested texture- or offscreen window format for drawing into it.\n");
			printf("PTB-ERROR: Most graphics cards do not support drawing into textures or offscreen windows which are not true-color, i.e. 1 layer pure\n");
			printf("PTB-ERROR: luminance or 2 layer luminance+alpha textures or offscreen windows won't work. Choose a 3-layer RGB or 4-layer RGBA texture\n");
			printf("PTB-ERROR: or offscreen window and retry.\n");
		}
		else if (fborc==GL_FRAMEBUFFER_UNSUPPORTED_EXT) {
			printf("PTB-ERROR: Your graphics hardware does not support the provided or requested texture- or offscreen window format for drawing into it.\n");
			printf("PTB-ERROR: Could be that the specific color depth of the texture or offscreen window is not supported for drawing on your hardware:\n");
			printf("PTB-ERROR: 8 bits per color component are supported on nearly all hardware, 16 bpc or 32 bpc floating point formats only on recent\n");
			printf("PTB-ERROR: hardware. 16 bpc fixed precision is not supported on any NVidia hardware, and on many systems only under restricted conditions.\n");
			printf("PTB-ERROR: Retry with the lowest acceptable (for your study) size and depth of the onscreen window or offscreen window.\n");
		}
		else {
			printf("PTB-ERROR: Exact reason for failure is unknown, most likely a Psychtoolbox bug, OpenGL-driver bug or unintented use. glCheckFramebufferStatus() returns code %i\n", fborc);
		}
		return(FALSE);
	}
	else if (fborc == 0) {
		// Checking command itself failed?!?
		if (PsychPrefStateGet_Verbosity() > 2) {
			// Warn the user:
			printf("PTB-WARNING: In setup of framebuffer object color attachment: glCheckFramebufferStatusEXT() malfunctioned. (glGetError reports: %s)\n", gluErrorString(glGetError()));
			printf("PTB-WARNING: Therefore can't determine if FBO setup worked or not. Will continue and hope for the best :(\n");
			printf("PTB-WARNING: This is most likely a graphics driver bug. You may want to update your graphics drivers, maybe it helps.\n");
		}
		while(glGetError());
	}
	
	// Do we need additional buffers for 3D rendering?
	if (needzbuffer) {
		// Yes. Try to setup and attach them: We use depth textures if they are supported and no MSAA is needed:
		if ((multisample <= 0) && (glewIsSupported("GL_ARB_depth_texture") || strstr(glGetString(GL_EXTENSIONS), "GL_OES_depth_texture"))) {
			// No multisampling requested on FBO and depth textures are supported. Use those to implement depth + stencil buffers:
			if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Trying to attach texture depth+stencil attachments to FBO...\n"); 
			
			// Create texture object for z-buffer (or z+stencil buffer) and set it up:
			glGenTextures(1, (GLuint*) &((*fbo)->ztexid));
			glBindTexture(texturetarget, (*fbo)->ztexid);
			
			// Setup texture wrapping behaviour to clamp, as other behaviours are
			// unsupported on many gfx-cards for rectangle textures:
            if (!PsychIsGLES(NULL)) {
                glTexParameterf(texturetarget,GL_TEXTURE_WRAP_S,GL_CLAMP);
                glTexParameterf(texturetarget,GL_TEXTURE_WRAP_T,GL_CLAMP);
			}
            else {
                glTexParameterf(texturetarget,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
                glTexParameterf(texturetarget,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
            }

			// Setup filtering for the textures - Use nearest neighbour by default, as floating
			// point filtering usually unsupported.
			glTexParameterf(texturetarget, GL_TEXTURE_MAG_FILTER,GL_NEAREST);
			glTexParameterf(texturetarget, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
			
			// Just to be safe...
			PsychTestForGLErrors();
			
			// Do we have support for combined 24 bit depth and 8 bit stencil buffer textures?
			if (glewIsSupported("GL_EXT_packed_depth_stencil") || (glewIsSupported("GL_NV_packed_depth_stencil") && glewIsSupported("GL_SGIX_depth_texture")) ||
                strstr(glGetString(GL_EXTENSIONS), "GL_OES_packed_depth_stencil")) {
				// Yes! Create combined depth and stencil texture:
				if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: packed_depth_stencil supported. Attaching combined 24 bit depth + 8 bit stencil texture...\n"); 
				
				// Create proper texture: Just allocate proper format, don't assign data.
                if (PsychIsGLES(NULL)) {
					glTexImage2D(texturetarget, 0, GL_DEPTH_STENCIL_EXT, twidth, theight, 0, GL_DEPTH_STENCIL_EXT, GL_UNSIGNED_INT_24_8_EXT, NULL);
                }
                else if (glewIsSupported("GL_EXT_packed_depth_stencil")) {
					glTexImage2D(texturetarget, 0, GL_DEPTH24_STENCIL8_EXT, twidth, theight, 0, GL_DEPTH_STENCIL_EXT, GL_UNSIGNED_INT_24_8_EXT, NULL);
				}
				else {
					// Ancient drivers with only NV extension support...
					glTexImage2D(texturetarget, 0, GL_DEPTH_COMPONENT24_SGIX, twidth, theight, 0, GL_DEPTH_STENCIL_EXT, GL_UNSIGNED_INT_24_8_EXT, NULL);
				}
				
				PsychTestForGLErrors();
				
				// Texture ready, unbind it.
				glBindTexture(texturetarget, 0);
				
				// Attach the texture as depth buffer...
				glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, texturetarget, (*fbo)->ztexid, 0);
				PsychTestForGLErrors();
				
				if (!(PsychPrefStateGet_ConserveVRAM() & kPsychDontAttachStencilToFBO)) {
					// ... and as stencil buffer ...
					glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_STENCIL_ATTACHMENT_EXT, texturetarget, (*fbo)->ztexid, 0);
					if (glGetError()) {
						// Attaching stencil buffer doesnt work :( We try to live without it...
						while(glGetError());
						if (PsychPrefStateGet_Verbosity()>1) {
							printf("PTB-WARNING: OpenGL stencil buffers not supported in imagingmode by your hardware. This won't affect Screen 2D drawing functions and won't affect\n");
							printf("PTB-WARNING: the majority of OpenGL (MOGL) 3D drawing code either, but OpenGL code that needs a stencil buffer will misbehave or fail in random ways!\n");
							printf("PTB-WARNING: If you need to use such code, you'll either have to disable the internal imaging pipeline, or carefully work-around this limitation by\n");
							printf("PTB-WARNING: proper modifications and testing of the affected code. Good luck... Alternatively, upgrade your graphics hardware or drivers.\n");
						}
					}
				}
				else {
					// Override: Do not attach stencil attachment!
					(*fbo)->stexid = 0;
					if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: PsychCreateFBO(): Won't attach a stencil buffer to FBO due to user override...\n"); 
				}
			}
			else {
				// Packed depth+stencil textures unsupported :( 
				// Allocate a single depth texture and attach it. Then see what we can do about the stencil attachment...
				if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: packed_depth_stencil unsupported. Attaching 24 bit depth texture and 8 bit stencil renderbuffer...\n"); 
				
				// Create proper texture: Just allocate proper format, don't assign data.
				glTexImage2D(texturetarget, 0, GL_DEPTH_COMPONENT, twidth, theight, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, NULL);
				
				// Depth texture ready, unbind it.
				glBindTexture(texturetarget, 0);
				
				// Attach the texture as depth buffer...
				glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, texturetarget, (*fbo)->ztexid, 0);
				
				if (!(PsychPrefStateGet_ConserveVRAM() & kPsychDontAttachStencilToFBO)) {
					// Create and attach renderbuffer as a stencil buffer of 8 bit depths:
					glGenRenderbuffersEXT(1, (GLuint*) &((*fbo)->stexid));
					glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, (*fbo)->stexid);
					glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_STENCIL_INDEX8_EXT, twidth, theight);
					glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
					glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_STENCIL_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, (*fbo)->stexid);
				}
				else {
					// Override: Do not attach stencil attachment!
					(*fbo)->stexid = 0;
					if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: PsychCreateFBO(): Pathway-2: Won't attach a stencil buffer to FBO due to user override...\n"); 
				}
				
				// See if we are framebuffer complete:
				fborc = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
				if (GL_FRAMEBUFFER_COMPLETE_EXT != fborc && 0 != fborc) {
					// Nope. Our trick doesn't work, this hardware won't let us attach a stencil buffer at all. Remove it
					// and live with a depth-buffer only setup.
					if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Stencil renderbuffer attachment failed! Detaching stencil buffer...\n"); 
					if (PsychPrefStateGet_Verbosity()>1) {
						printf("PTB-WARNING: OpenGL stencil buffers not supported in imagingmode by your hardware. This won't affect Screen 2D drawing functions and won't affect\n");
						printf("PTB-WARNING: the majority of OpenGL (MOGL) 3D drawing code either, but OpenGL code that needs a stencil buffer will misbehave or fail in random ways!\n");
						printf("PTB-WARNING: If you need to use such code, you'll either have to disable the internal imaging pipeline, or carefully work-around this limitation by\n");
						printf("PTB-WARNING: proper modifications and testing of the affected code. Good luck... Alternatively, upgrade your graphics hardware or drivers. According to specs,\n");
						printf("PTB-WARNING: all gfx-cards starting with GeForceFX 5200 on Windows and Linux and all cards on Intel-Macs except the Intel GMA cards should work, whereas\n");
						printf("PTB-WARNING: none of the PowerPC hardware is supported as of OS-X 10.4.9.\n"); 
					}
					
					glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_STENCIL_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, 0);
					glDeleteRenderbuffersEXT(1, (GLuint*) &((*fbo)->stexid));
					(*fbo)->stexid = 0;
				}
				else if (fborc == 0) {
					// Checking command itself failed?!?
					if (PsychPrefStateGet_Verbosity() > 2) {
						// Warn the user:
						printf("PTB-WARNING: In setup of framebuffer object stencil attachment: glCheckFramebufferStatusEXT() malfunctioned. (glGetError reports: %s)\n", gluErrorString(glGetError()));
						printf("PTB-WARNING: Therefore can't determine if FBO setup worked or not. Will continue and hope for the best :(\n");
						printf("PTB-WARNING: This is most likely a graphics driver bug. You may want to update your graphics drivers, maybe it helps.\n");
					}
					while(glGetError());
				}
			}
		}	// End of setup code for multiSample == 0.
		else {
			// Setup code of z- and stencilbuffer for multisampled mode. We must allocate these attachments as renderbuffers,
			// as they need the same sample count as the color buffers, and textures do not support multisampling.
			// We also use this path for non-multisampled renderbuffer setup if depth textures aren't available.
			if (PsychPrefStateGet_Verbosity()>4) {
				printf("PTB-DEBUG: Trying to attach %s renderbuffer depth+stencil attachments to FBO...\n", (multisample > 0) ? "multisample" : ""); 
			}

			glGenRenderbuffersEXT(1, (GLuint*) &((*fbo)->ztexid));
			glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, (*fbo)->ztexid);

            // Use of depth+stencil requested and packed depth+stencil possible?
            if (!(PsychPrefStateGet_ConserveVRAM() & kPsychDontAttachStencilToFBO) && 
                (glewIsSupported("GL_EXT_packed_depth_stencil") || strstr(glGetString(GL_EXTENSIONS), "GL_OES_packed_depth_stencil"))) {
                // Try a packed depth + stencil buffer:
                if (multisample > 0) {
                    glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, multisample, GL_DEPTH24_STENCIL8_EXT, twidth, theight);
                }
                else {
                    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH24_STENCIL8_EXT, twidth, theight);
                }

                (*fbo)->stexid = (*fbo)->ztexid;
            }
            else {
                // Depth buffer only:
                if (multisample > 0) {
                    glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, multisample, GL_DEPTH_COMPONENT24_ARB, twidth, theight);
                }
                else {
                    // OES doesn't guarantee 24 bit depth buffers, only 16 bit, so fallback to them on OES if 24 bits are unsupported:
                    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT,
                                             (!PsychIsGLES(NULL) || strstr(glGetString(GL_EXTENSIONS), "GL_OES_depth24")) ? GL_DEPTH_COMPONENT24_ARB : GL_DEPTH_COMPONENT16_ARB,
                                             twidth, theight);
                }

                (*fbo)->stexid = 0;
            }
            
			if (glGetError()!=GL_NO_ERROR) {
				printf("PTB-ERROR: Failed to setup internal framebuffer objects depths renderbuffer attachment as a %s renderbuffer for imaging pipeline!\n", (multisample > 0) ? "multisample" : "");
				printf("PTB-ERROR: Most likely the requested size, depth and multisampling level is not supported by your graphics hardware.\n");
				return(FALSE);
			}
			
            if (multisample > 0) {
                // Query real number of samples for renderbuffer:
                glGetRenderbufferParameterivEXT(GL_RENDERBUFFER_EXT, GL_RENDERBUFFER_SAMPLES_EXT, &multisample);
                if ((*fbo)->multisample != multisample) {
                    printf("PTB-ERROR: Failed to setup internal framebuffer objects depths renderbuffer attachment as a multisampled renderbuffer for imaging pipeline!\n");
                    printf("PTB-ERROR: Could not get the same number of samples for depths renderbuffer as for color buffer, which is a requirement.\n");
                    return(FALSE);
                }
			}

			// Unbind, we're done with setup:
			glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
			
            // Non-zero stexid marks that we've created a packed depth+stencil renderbuffer above:
            if ((*fbo)->stexid) {
                // Attach combined z + stencil buffer:
                glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER_EXT, (*fbo)->ztexid);                
            }
            else {
                // Attach z-buffer only:
                glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, (*fbo)->ztexid);
			}
            
			// Now try to attach stencil buffer:
			if (!((*fbo)->stexid) && !(PsychPrefStateGet_ConserveVRAM() & kPsychDontAttachStencilToFBO)) {
				// Create and attach renderbuffer as a stencil buffer of hopefully 8 bit depths:
				glGenRenderbuffersEXT(1, (GLuint*) &((*fbo)->stexid));
				glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, (*fbo)->stexid);

                // Try to get stencil buffer with matching sample count to depth and color buffers:
                if (multisample > 0) {
                    glRenderbufferStorageMultisampleEXT(GL_RENDERBUFFER_EXT, multisample, GL_STENCIL_INDEX8_EXT, twidth, theight);
                }
                else {
                    // OES does not guarantee 8 bit stencil buffers, only with extension. In fact, it does not even guarantee
                    // stencil buffers at all, so we aim for 8 bit, fallback to 4 bit and leave it to the error handling below to disable
                    // stencil buffers if the implementation doesn't support at least the 4 bit case:
                    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, (!PsychIsGLES(NULL) || strstr(glGetString(GL_EXTENSIONS), "GL_OES_stencil8")) ? GL_STENCIL_INDEX8_EXT : GL_STENCIL_INDEX4_EXT, twidth, theight);
                }

				if (glGetError()!=GL_NO_ERROR) {
					if (PsychPrefStateGet_Verbosity() > 2) {
                        printf("PTB-WARNING: Failed to setup internal framebuffer objects stencil renderbuffer attachment as a %s renderbuffer for imaging pipeline!\n", (multisample > 0) ? "multisample" : "");
                    }
                    
                    // Clean up:
                    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
					glDeleteRenderbuffersEXT(1, (GLuint*) &((*fbo)->stexid));
					(*fbo)->stexid = 0;
				}

                if (multisample > 0) {				
                    // Query real number of samples for renderbuffer:
                    glGetRenderbufferParameterivEXT(GL_RENDERBUFFER_EXT, GL_RENDERBUFFER_SAMPLES_EXT, &multisample);
                    if ((*fbo)->multisample != multisample) {
                        if (PsychPrefStateGet_Verbosity() > 2) {
                            printf("PTB-WARNING: Failed to setup internal framebuffer objects stencil renderbuffer attachment as a multisampled renderbuffer for imaging pipeline!\n");
                            printf("PTB-WARNING: Could not get the same number of samples for stencil renderbuffer as for color buffer, which is a requirement.\n");
                        }

                        // Clean up:
                        glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
                        glDeleteRenderbuffersEXT(1, (GLuint*) &((*fbo)->stexid));
                        (*fbo)->stexid = 0;
                    }
				}

                if ((*fbo)->stexid > 0) {
                    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
                    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_STENCIL_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, (*fbo)->stexid);
                    
                    // See if we are framebuffer complete:
                    fborc = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
                    if (GL_FRAMEBUFFER_COMPLETE_EXT != fborc && 0 != fborc) {
                        // Nope. Our trick doesn't work, this hardware won't let us attach a stencil buffer at all. Remove it
                        // and live with a depth-buffer only setup.
                        if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Stencil %s renderbuffer attachment failed! Detaching stencil buffer...\n", (multisample > 0) ? "multisample" : "");                        
                        glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_STENCIL_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, 0);
                        glDeleteRenderbuffersEXT(1, (GLuint*) &((*fbo)->stexid));
                        (*fbo)->stexid = 0;
                    }
                    else if (fborc == 0) {
                        // Checking command itself failed?!?
                        if (PsychPrefStateGet_Verbosity() > 1) {
                            // Warn the user:
                            printf("PTB-WARNING: In setup of framebuffer object %s stencil attachment: glCheckFramebufferStatusEXT() malfunctioned. (glGetError reports: %s)\n",
                                   (multisample > 0) ? "multisample" : "", gluErrorString(glGetError()));
                            printf("PTB-WARNING: Therefore can't determine if FBO setup worked or not. Will continue and hope for the best :(\n");
                            printf("PTB-WARNING: This is most likely a graphics driver bug. You may want to update your graphics drivers, maybe it helps.\n");
                        }
                        while(glGetError());
                    }
                }
                
                if (((*fbo)->stexid == 0) && (PsychPrefStateGet_Verbosity() > 1)) {
                    printf("PTB-WARNING: OpenGL %s stencil buffers not supported in imagingmode by your hardware. This won't affect Screen 2D drawing functions and won't affect\n",
                           (multisample > 0) ? "multisample" : "");
                    printf("PTB-WARNING: the majority of OpenGL (MOGL) 3D drawing code either, but OpenGL code that needs a stencil buffer will misbehave or fail in random ways!\n");
                    printf("PTB-WARNING: If you need to use such code, you'll either have to disable the internal imaging pipeline, or carefully work-around this limitation by\n");
                    printf("PTB-WARNING: proper modifications and testing of the affected code. Good luck... Alternatively, upgrade your graphics hardware or drivers.\n\n");
                }
			}
			else {
				// Override: Do not attach stencil attachment, or combined depth+stencil attached?
                if ((*fbo)->stexid) {
                    // Non-Zero stexid indicates we use a combined depth+stencil renderbuffer. Its handle
                    // is already stored in ztexid, so zero-out stexid to avoid redundant destruction on
                    // framebuffer cleanup later on:
                    (*fbo)->stexid = 0;
                } // Zero stexid indicates that no stencil buffer shall be used, according to usercode's preferences:
                else if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: PsychCreateFBO(): Pathway-3: Won't attach a stencil buffer to FBO due to user override...\n"); 
			}
			
			// Ok, all depths- and stencil- renderbuffers with same number of multisamples as color renderbuffer attached. Check for completeness will
			// happen further down the road...
		}
	}
	else {
		// Initialize additional buffers to zero:
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Only colorbuffer texture attached to FBO, no depth- or stencil buffers requested...\n"); 
		(*fbo)->stexid = 0;
		(*fbo)->ztexid = 0;
	}
	
	// Store dimensions:
	(*fbo)->width = width;
	(*fbo)->height = height;
	
	// Check for framebuffer completeness:
	fborc = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if (fborc!=GL_FRAMEBUFFER_COMPLETE_EXT && fborc!=0) {
		// Framebuffer incomplete!
		while(glGetError());
        
        switch(fborc) {
            case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
                sprintf(fbodiag, "Unsupported format");
            break;

            case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
                sprintf(fbodiag, "Framebuffer attachment incomplete.");
            break;
            
            case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE:
                sprintf(fbodiag, "Framebuffer attachment multisample incomplete.");
            break;

            case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
                sprintf(fbodiag, "Framebuffer attachments missing incomplete.");
            break;

            case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
                sprintf(fbodiag, "Framebuffer dimensions incomplete.");
            break;

            case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
                sprintf(fbodiag, "Framebuffer formats incomplete.");
            break;

            case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
                sprintf(fbodiag, "Framebuffer drawbuffer incomplete.");
            break;

            case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
                sprintf(fbodiag, "Framebuffer readbuffer incomplete.");
            break;

            default:
                sprintf(fbodiag, "Unknown error: glCheckFramebufferStatusEXT returns error code %i", fborc);                
        }
        
		printf("PTB-ERROR: Failed to setup internal framebuffer object for imaging pipeline [%s]! The most likely cause is that your hardware does not support\n", fbodiag);
		printf("PTB-ERROR: the required buffers at the given screen resolution (Additional 3D buffers for z- and stencil are %s).\n", (needzbuffer) ? "requested" : "disabled");
		printf("PTB-ERROR: You may want to retry with the lowest acceptable (for your study) display resolution or with 3D rendering support disabled.\n");
		return(FALSE);
	}
	else if (fborc == 0) {
		// Checking command itself failed?!?
		if (PsychPrefStateGet_Verbosity() > 2) {
			// Warn the user:
			printf("PTB-WARNING: In setup of framebuffer object: glCheckFramebufferStatusEXT() malfunctioned. (glGetError reports: %s)\n", gluErrorString(glGetError()));
			printf("PTB-WARNING: Therefore can't determine if FBO setup worked or not. Will continue and hope for the best :(\n");
			printf("PTB-WARNING: This is most likely a graphics driver bug. You may want to update your graphics drivers, maybe it helps.\n");
		}
		while(glGetError());
	}

	if (PsychPrefStateGet_Verbosity()>4) {
		// Output framebuffer properties:
		glGetIntegerv(GL_RED_BITS, &bpc);
		printf("PTB-DEBUG: FBO has %i bits precision per color component in ", bpc);
		if (glewIsSupported("GL_ARB_color_buffer_float") || glewIsSupported("GL_EXT_color_buffer_float")) {
			glGetBooleanv(GL_RGBA_FLOAT_MODE_ARB, &isFloatBuffer);
			if (isFloatBuffer) {
				printf("floating point format ");
			}
			else {
				printf("fixed point format ");
			}
		}
		else if (glewIsSupported("GL_APPLE_float_pixels")) { 
			glGetBooleanv(GL_COLOR_FLOAT_APPLE, &isFloatBuffer);
			if (isFloatBuffer) {
				printf("floating point format ");
			}
			else {
				printf("fixed point format ");
			}
		}
		else {
			isFloatBuffer = FALSE;
			printf("unknown (but likely fixed point) format ");
		}

		glGetIntegerv(GL_DEPTH_BITS, &bpc);
		printf("with  %i depths buffer bits ", bpc);
		glGetIntegerv(GL_STENCIL_BITS, &bpc);
		printf("and %i stencil buffer bits.\n", bpc);
	}

	// Unbind FBO:
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// Test all GL errors:
	PsychTestForGLErrors();

	// Well done.
	return(TRUE);
}

/* PsychCreateShadowFBOForTexture()
 * Check if provided PTB texture already has a PsychFBO attached. Do nothing if so.
 * If a FBO is missing, create one.
 *
 * If asRendertarget is FALSE, we only create the data structure, not a real OpenGL FBO,
 * so the texture is only suitable as image source for image processing.
 *
 * If asRendertraget is TRUE, we create a full blown FBO, so the texture can be used as
 * rendertarget.
 *
 */
void PsychCreateShadowFBOForTexture(PsychWindowRecordType *textureRecord, psych_bool asRendertarget, int forImagingmode)
{
	GLenum fboInternalFormat;
	
	// Do we already have a framebuffer object for this texture? All textures start off without one,
	// because most textures are just used for drawing them, not drawing *into* them. Therefore we
	// only create a full blown FBO on demand here.
	if (textureRecord->drawBufferFBO[0]==-1) {
		// No. This texture is used the first time as a drawing target.
		// Need to create a framebuffer object for it first.
		
		if (textureRecord->textureNumber > 0) {
			// Allocate and assign FBO object info structure PsychFBO:
			PsychCreateFBO(&(textureRecord->fboTable[0]), (GLenum) 0, (PsychPrefStateGet_3DGfx() > 0) ? TRUE : FALSE, (int) PsychGetWidthFromRect(textureRecord->rect), (int) PsychGetHeightFromRect(textureRecord->rect), 0, (PsychGetTextureTarget(textureRecord) == GL_TEXTURE_2D) ? 1 : 0);
			
			// Manually set up the color attachment texture id to our texture id:
			textureRecord->fboTable[0]->coltexid = textureRecord->textureNumber;
		}
		else {
			// No texture yet. Create suitable one for given imagingmode:

			// Start off with standard 8 bpc fixed point:
			fboInternalFormat = GL_RGBA8; textureRecord->bpc = 8;
			
			// Need 16 bpc fixed point precision?
			if (forImagingmode & kPsychNeed16BPCFixed) { fboInternalFormat = ((textureRecord->gfxcaps & kPsychGfxCapSNTex16) ? GL_RGBA16_SNORM : GL_RGBA16); textureRecord->bpc = 16; }
			
			// Need 16 bpc floating point precision?
			if (forImagingmode & kPsychNeed16BPCFloat) { fboInternalFormat = GL_RGBA_FLOAT16_APPLE; textureRecord->bpc = 16; }
			
			// Need 32 bpc floating point precision?
			if (forImagingmode & kPsychNeed32BPCFloat) { fboInternalFormat = GL_RGBA_FLOAT32_APPLE; textureRecord->bpc = 32; }
			
            // OpenGL-ES + float precision requested?
            if (PsychIsGLES(textureRecord) && (textureRecord->bpc > 8)) {
                // Upgrade to full 32bpc float -- we only support this, if at all:
                fboInternalFormat = GL_RGBA_FLOAT32_APPLE; textureRecord->bpc = 32;
            }

			PsychCreateFBO(&(textureRecord->fboTable[0]), fboInternalFormat, (PsychPrefStateGet_3DGfx() > 0) ? TRUE : FALSE, (int) PsychGetWidthFromRect(textureRecord->rect), (int) PsychGetHeightFromRect(textureRecord->rect), 0, (PsychGetTextureTarget(textureRecord) == GL_TEXTURE_2D) ? 1 : 0);
			
			// Manually set up the texture id from our color attachment texture id:
			textureRecord->textureNumber = textureRecord->fboTable[0]->coltexid;
		}
				
		// Worked. Set up remaining state:
		textureRecord->fboCount = 1;
		textureRecord->drawBufferFBO[0]=0;
	}
	
	// Does it need to be suitable as a rendertarget? If so, check if it is already, upgrade it if neccessary:
	if (asRendertarget && textureRecord->fboTable[0]->fboid==0) {
		// Initialize and setup real FBO object (optionally with z- and stencilbuffer) and attach the texture
		// as color attachment 0, aka main colorbuffer:				
		if (!PsychCreateFBO(&(textureRecord->fboTable[0]), (GLenum) 1, (PsychPrefStateGet_3DGfx() > 0) ? TRUE : FALSE, (int) PsychGetWidthFromRect(textureRecord->rect), (int) PsychGetHeightFromRect(textureRecord->rect), 0, (PsychGetTextureTarget(textureRecord) == GL_TEXTURE_2D) ? 1 : 0)) {
			// Failed!
			PsychErrorExitMsg(PsychError_internal, "Preparation of drawing into an offscreen window or texture failed when trying to create associated framebuffer object!");
		}
	}
	
	return;
}

/* PsychNormalizeTextureOrientation() - On demand texture reswapping.
 *
 * PTB supports multiple different ways of orienting and formatting textures,
 * optimized for different purposes. However, textures that one wants to draw
 * to as if they were Offscreen windows or textures to be fed into the image
 * processing pipeline (Screen('TransformTexture')) need to be in a standardized
 * upright, non-transposed format. This routine checks the orientation of a
 * texture and - if neccessary - transforms the texture from its current format
 * to the upright standard format. As a side effect, it also converts such textures
 * from Luminance or Luminance+Alpha formats into RGB or RGBA formats. This is
 * important, because only RGB(A) textures are suitable as FBO color buffer attachments.
 *
 */ 
void PsychNormalizeTextureOrientation(PsychWindowRecordType *sourceRecord)
{
	int tmpimagingmode;
	PsychFBO *fboptr;
	GLint fboInternalFormat;
	psych_bool needzbuffer, isplanar;
	int width, height;

    // Is this a planar encoding texture?
    isplanar = (PsychIsTexture(sourceRecord) && (sourceRecord->specialflags & kPsychPlanarTexture)) ? TRUE : FALSE;
    
	// The source texture sourceRecord could be in any of PTB's supported
	// internal texture orientations. It may be upright as an Offscreen window,
	// or flipped upside down as some textures from the video grabber or Quicktime,
	// or transposed, as textures from Matlab/Octave. However, handling all those
	// cases for image processing would be a debug and maintenance nightmare.
	// Therefore we check the format of the source texture and require it to be
	// a normal upright orientation. If this isn't the case, we perform a preprocessing
	// step to transform the texture into normalized orientation. We also perform a
	// preprocessing step on any CoreVideo texture from Quicktime. Although such a
	// texture may be properly oriented, it is of a non-renderable YUV color format, so
	// we need to recreate it in a RGB renderable format. Non-planar textures would also
    // wreak havoc if not converted into standard pixel-interleaved format:
	if (sourceRecord->textureOrientation != 2 || sourceRecord->targetSpecific.QuickTimeGLTexture != NULL || isplanar) {
		if (PsychPrefStateGet_Verbosity()>5) printf("PTB-DEBUG: In PsychNormalizeTextureOrientation(): Performing GPU renderswap or format conversion for source gl-texture %i --> ", sourceRecord->textureNumber);
		
		// Soft-reset drawing engine in a safe way:
		PsychSetDrawingTarget((PsychWindowRecordType*) 0x1);

		// Normalization needed. Create a suitable FBO as rendertarget:
		needzbuffer = FALSE;
		
		// First delete FBO of texture if one already exists:
		fboptr = sourceRecord->fboTable[0];
		if (fboptr!=NULL) { 			
			// Detach and delete color buffer texture:
			needzbuffer = (fboptr->ztexid) ? TRUE : FALSE;

			// Detach and delete depth buffer (and probably stencil buffer) texture, if any:
			if (fboptr->ztexid) {
                if (glIsTexture(fboptr->ztexid)) {
                    glDeleteTextures(1, &(fboptr->ztexid));
                }
                else {
                    glDeleteRenderbuffersEXT(1, &(fboptr->ztexid));
                }
            }

			// Detach and delete stencil renderbuffer, if a separate stencil buffer was needed:
			if (fboptr->stexid) glDeleteRenderbuffersEXT(1, &(fboptr->stexid));

			// Delete FBO itself:
			if (fboptr->fboid) glDeleteFramebuffersEXT(1, &(fboptr->fboid));
						
			// Delete PsychFBO struct associated with this FBO:
			free(fboptr); fboptr = NULL;
		}
		
		// Now recreate with proper format:

		// First need to know internal format of texture...
		glBindTexture(PsychGetTextureTarget(sourceRecord), sourceRecord->textureNumber);
		glGetTexLevelParameteriv(PsychGetTextureTarget(sourceRecord), 0, GL_TEXTURE_INTERNAL_FORMAT, &fboInternalFormat);

		// Need to query real size of underlying texture, not the logical size from sourceRecord->rect, otherwise we'd screw
		// up for padded textures (from Quicktime movie/vidcap) where the real texture is a bit bigger than its logical size.
		if (sourceRecord->textureOrientation > 1) {
			// Non-transposed textures, width and height are correct:
			glGetTexLevelParameteriv(PsychGetTextureTarget(sourceRecord), 0, GL_TEXTURE_WIDTH, &width);
			glGetTexLevelParameteriv(PsychGetTextureTarget(sourceRecord), 0, GL_TEXTURE_HEIGHT, &height);
		}
		else {
			// Transposed textures: Need to swap meaning of width and height:
			glGetTexLevelParameteriv(PsychGetTextureTarget(sourceRecord), 0, GL_TEXTURE_WIDTH, &height);
			glGetTexLevelParameteriv(PsychGetTextureTarget(sourceRecord), 0, GL_TEXTURE_HEIGHT, &width);
		}
		
		glBindTexture(PsychGetTextureTarget(sourceRecord), 0);

		// Override detected width and height for planar textures:
		if (isplanar) {
			// They are actually the net size as specified by their rect's:
			width  = (int) PsychGetWidthFromRect(sourceRecord->rect);
			height = (int) PsychGetHeightFromRect(sourceRecord->rect);
		}

        // Special mapping for OpenGL-ES:
        if (PsychIsGLES(sourceRecord)) {
            // Enforce RGBA8 mapping, so all code below proceeds as if this is the case:
            fboInternalFormat = GL_RGBA8;

            // OpenGL-ES + float precision requested?
            if ((sourceRecord->bpc > 8) || ((sourceRecord->depth / sourceRecord->nrchannels) >= 32)) {
                // Upgrade to full 32bpc float -- we only support this, if at all:
                fboInternalFormat = GL_RGBA_FLOAT32_APPLE;

                // This will force code into proper path below:
                sourceRecord->textureexternalformat = fboInternalFormat;
            }

			// Override width and height: Can't query size from OpenGL-ES due to lack of support,
			// get net width and height for recreation of FBO:
			width  = (int) PsychGetWidthFromRect(sourceRecord->rect);
			height = (int) PsychGetHeightFromRect(sourceRecord->rect);
        }

		// Renderable format? Pure luminance or luminance+alpha formats are not renderable on most hardware.

		// Upgrade 8 bit luminace to 8 bit RGBA:
		if (fboInternalFormat == GL_LUMINANCE8 || fboInternalFormat == GL_LUMINANCE8_ALPHA8 || sourceRecord->depth == 8) fboInternalFormat = GL_RGBA8;
		
		// Upgrade non-renderable floating point formats to their RGB or RGBA counterparts of matching precision:
		if (sourceRecord->nrchannels < 3 && fboInternalFormat != GL_RGBA8) {
			// Unsupported format for FBO rendertargets. Need to upgrade to something suitable...
			if (sourceRecord->textureexternalformat == GL_LUMINANCE) {
				// Upgrade luminance to RGB of matching precision:
				// printf("UPGRADING TO RGBFloat %i\n", (sourceRecord->textureinternalformat == GL_LUMINANCE_FLOAT16_APPLE) ? 0:1);
                if (sourceRecord->textureinternalformat == GL_LUMINANCE16_SNORM) {
                    fboInternalFormat = GL_RGB16_SNORM;
                } else {
                    fboInternalFormat = (sourceRecord->textureinternalformat == GL_LUMINANCE_FLOAT16_APPLE) ? GL_RGB_FLOAT16_APPLE : GL_RGB_FLOAT32_APPLE;
                }
			}
			else {
				// Upgrade luminance+alpha to RGBA of matching precision:
				// printf("UPGRADING TO RGBAFloat %i\n", (sourceRecord->textureinternalformat == GL_LUMINANCE_ALPHA_FLOAT16_APPLE) ? 0:1);
                if (sourceRecord->textureinternalformat == GL_LUMINANCE16_ALPHA16_SNORM) {
                    fboInternalFormat = GL_RGBA16_SNORM;
                } else {
                    fboInternalFormat = (sourceRecord->textureinternalformat == GL_LUMINANCE_ALPHA_FLOAT16_APPLE) ? GL_RGBA_FLOAT16_APPLE : GL_RGBA_FLOAT32_APPLE;
                }
			}
		}
		
		// Special case: Quicktime movie or video texture, created by CoreVideo in Apple specific YUV format.
		// This is a non-framebuffer renderable color format. Need to upgrade it to something safe:
		if (fboInternalFormat == GL_YCBCR_422_APPLE) fboInternalFormat = GL_RGBA8;
		
		// Special case: Planar texture encoded in a luminance texture other than LUMINANCE8. Need to
		// upgrade to a full RGBA format of sufficient precision:
		if (isplanar && (fboInternalFormat != GL_RGBA8)) {
			if (sourceRecord->textureinternalformat == GL_LUMINANCE16_SNORM) {
				fboInternalFormat = GL_RGBA16_SNORM;
			} else {
				fboInternalFormat = (sourceRecord->textureinternalformat == GL_LUMINANCE_FLOAT16_APPLE) ? GL_RGBA_FLOAT16_APPLE : GL_RGBA_FLOAT32_APPLE;
			}
		}
		
		// Now create proper FBO:
		if (!PsychCreateFBO(&(sourceRecord->fboTable[0]), (GLenum) fboInternalFormat, needzbuffer, width, height, 0, (PsychGetTextureTarget(sourceRecord) == GL_TEXTURE_2D) ? 1 : 0)) {
			PsychErrorExitMsg(PsychError_internal, "Failed to normalize texture orientation - Creation of framebuffer object failed!");
		}
		
		sourceRecord->drawBufferFBO[0] = 0;
		sourceRecord->fboCount = 1;
		
		tmpimagingmode = sourceRecord->imagingMode;
		sourceRecord->imagingMode = 1;

		// Set FBO of sourceRecord as rendertarget, including proper setup of render geometry:
		// We can't use PsychSetDrawingTarget() here, as we might get called by that function, i.e.
		// infinite recursion or other side effects if we tried to use it.
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, sourceRecord->fboTable[0]->fboid);
		PsychSetupView(sourceRecord, FALSE);
		// Reset MODELVIEW matrix, after backing it up...
		glPushMatrix();
		glLoadIdentity();
		
		// For planar textures we need to bind the planar -> interleaved conversion shader:
		if (isplanar) {
            if ((sourceRecord->textureFilterShader == 0) && (PsychPrefStateGet_Verbosity() > 1)) {
                printf("PTB-WARNING: Failed to normalize texture orientation and format: Conversion shader missing for this special texture!\n");
            }
            
            PsychSetShader(sourceRecord, -1 * sourceRecord->textureFilterShader);
		}
        
		// Now blit the old "disoriented" texture into the new FBO: The textureNumber of sourceRecord
		// references the old texture, the PsychFBO of sourceRecord defines the new texture...
		if (glIsEnabled(GL_BLEND)) {
			// Alpha blending enabled. Disable it, blit texture, reenable it:
			glDisable(GL_BLEND);
			PsychBlitTextureToDisplay(sourceRecord, sourceRecord, sourceRecord->rect, sourceRecord->rect, 0, 0, 1);
			glEnable(GL_BLEND);
		}
		else {
			// Alpha blending not enabled. Just blit it:
			PsychBlitTextureToDisplay(sourceRecord, sourceRecord, sourceRecord->rect, sourceRecord->rect, 0, 0, 1);
		}

		// Reset shader binding:
		if (isplanar) {
			PsychSetShader(sourceRecord, 0);
			
			// Now the texture has been turned into a regular pixel-interleaved texture,
			// it is no longer a planar texture, so we remove the planar->interleave
			// conversion shader and clear the planar texture flag:
			sourceRecord->textureFilterShader = 0;
			sourceRecord->specialflags &= ~kPsychPlanarTexture;
		}

		// Restore modelview matrix:
		glPopMatrix();
		
		sourceRecord->imagingMode = tmpimagingmode;
		PsychSetDrawingTarget(NULL);

		// At this point the color attachment of the sourceRecords FBO contains the properly oriented texture.
		// Delete the old texture, attach the FBO texture as new one:
		if (sourceRecord->targetSpecific.QuickTimeGLTexture != NULL) {
			// Special case: CoreVideo texture:
			PsychFreeMovieTexture(sourceRecord);
		}
		else {
			// Standard case:
			glDeleteTextures(1, &(sourceRecord->textureNumber));
		}
		
		// Assign new texture:
		sourceRecord->textureNumber = sourceRecord->fboTable[0]->coltexid;
		
		// Finally sourceRecord has the proper orientation:
		sourceRecord->textureOrientation = 2;
		
		// GPU renderswap finished.
		if (PsychPrefStateGet_Verbosity()>5) printf("%i.\n", sourceRecord->textureNumber);
	}

	return;
}

/* PsychShutdownImagingPipeline()
 * Shutdown imaging pipeline for a windowRecord and free all ressources associated with it.
 */
void PsychShutdownImagingPipeline(PsychWindowRecordType *windowRecord, psych_bool openglpart)
{
	int i;
	PsychFBO* fboptr;
	
	// Do OpenGL specific cleanup:
	if (openglpart) {
		// Yes. Mode specific cleanup:
		for (i=0; i<windowRecord->fboCount; i++) {
			// Delete i'th FBO, if any:
			fboptr = windowRecord->fboTable[i];
			if (fboptr!=NULL) { 
				// Delete all remaining references to this fbo:
				for (i=0; i<windowRecord->fboCount; i++) if (fboptr == windowRecord->fboTable[i]) windowRecord->fboTable[i] = NULL;
				
				// Detach and delete color buffer texture/renderbuffer:
				if (fboptr->coltexid) {
					if (glIsTexture(fboptr->coltexid)) {
						// Color buffer is a texture:
						glDeleteTextures(1, &(fboptr->coltexid));
					}
					else {
						// Color buffer is a renderbuffer:
						glDeleteRenderbuffersEXT(1, &(fboptr->coltexid));
					}
				}

				// Detach and delete depth buffer (and probably stencil buffer) texture, if any:
				if (fboptr->ztexid) {
					if (glIsTexture(fboptr->ztexid)) {
						// Depths buffer is a texture:
						glDeleteTextures(1, &(fboptr->ztexid));
					}
					else {
						// Depths buffer is a renderbuffer:
						glDeleteRenderbuffersEXT(1, &(fboptr->ztexid));
					}
				}

				// Detach and delete stencil renderbuffer, if a separate stencil buffer was needed:
				if (fboptr->stexid) glDeleteRenderbuffersEXT(1, &(fboptr->stexid));

				// Delete FBO itself:
				if (fboptr->fboid) glDeleteFramebuffersEXT(1, &(fboptr->fboid));
				
				// Delete PsychFBO struct associated with this FBO:
				free(fboptr); fboptr = NULL;
			}
		}
	} 

	// The following cleanup must only happen after OpenGL rendering context is already detached and
	// destroyed. It's part of phase-2 "post GL shutdown" of Screen('Close') and friends...
	if (!openglpart) {
		// Clear all hook chains:
		for (i=0; i<MAX_SCREEN_HOOKS; i++) {
			windowRecord->HookChainEnabled[i]=FALSE;
			PsychPipelineResetHook(windowRecord, PsychHookPointNames[i]);
		}
		
		// Global off:
		windowRecord->imagingMode=0;
	}
	
	// Cleanup done.
	return;
}

void PsychPipelineListAllHooks(PsychWindowRecordType *windowRecord)
{
	int i;
	printf("PTB-INFO: The Screen command currently provides the following hook functions:\n");
	printf("=============================================================================\n");
	for (i=0; i<MAX_SCREEN_HOOKS; i++) {
		printf("- %s : %s\n", PsychHookPointNames[i], PsychHookPointSynopsis[i]);
	}
	printf("=============================================================================\n\n");	
	fflush(NULL);
	// Well done.
	return;	
}

/* Map a hook name string to its corresponding hook chain index:
 * Returns -1 if such a hook name doesn't exist.
 */
int	PsychGetHookByName(const char* hookName)
{
	int i;
	for (i=0; i<MAX_SCREEN_HOOKS; i++) {
		if(strcmp(PsychHookPointNames[i], hookName)==0) break; 
	}
	return((i>=MAX_SCREEN_HOOKS) ? -1 : i);
}

/* Internal: PsychAddNewHookFunction()  - Add a new hook callback function to a hook-chain.
 * This helper function allocates a hook func struct, enqueues it into a hook chain and sets
 * all common struct fields to their proper values. Then it returns a pointer to the struct, so
 * the calling routine can fill the rest of the struct with information.
 *
 * windowRecord - The windowRecord of the window to attach to.
 * hookString   - String with the human-readable name of the hook-chain to attach to.
 * idString     - Arbitrary name string for this hook function (for debugging and finding it programmatically)
 * where        - Attachment point: 0 = Prepend to chain (Recommended). INT_MAX = Append to chain (Use with care!), or insert at position 'where'.
 * hookfunctype - Symbolic hook function type id, needed for interpreter to distinguish different types.
 *
 */
PsychHookFunction* PsychAddNewHookFunction(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, int hookfunctype)
{
	int hookidx, hookslotidx;
	PtrPsychHookFunction hookfunc, hookiter, *hookpreiter;
	
	// Lookup hook-chain idx for this name, if any:
	if ((hookidx=PsychGetHookByName(hookString))==-1) PsychErrorExitMsg(PsychError_user, "AddHook: Unknown (non-existent) hook name provided.");
	
	// Allocate a hook structure:
	hookfunc =	(PtrPsychHookFunction) calloc(1, sizeof(PsychHookFunction));
	if (hookfunc==NULL) PsychErrorExitMsg(PsychError_outofMemory, "Failed to allocate memory for new hook function.");
	
	// Enqueue at beginning or end of chain:
	if (where==0) {
		// Prepend it to existing chain:
		hookfunc->next = windowRecord->HookChain[hookidx];
		windowRecord->HookChain[hookidx] = hookfunc;
	}
	else {
		// Append it to existing chain:
		hookfunc->next = NULL;
		if (windowRecord->HookChain[hookidx]==NULL) {
			windowRecord->HookChain[hookidx]=hookfunc;
		}
		else {
			hookiter = windowRecord->HookChain[hookidx];
			hookslotidx = 0;
			hookpreiter = &(windowRecord->HookChain[hookidx]);
			
			while ((hookiter->next) && (hookslotidx!=where)) {
				hookpreiter = &(hookiter->next);
				hookiter = hookiter->next;
				hookslotidx++;
			}
			if (hookslotidx==where) {
				// Target slot index for insertion reached: Insert here.
				hookfunc->next = hookiter;
				*hookpreiter = hookfunc;
			}
			else {
				// End of chain reached: Append new slot, regardless of index:
				hookiter->next = hookfunc;
			}
		}
	}
	
	// New hookfunc struct is enqueued and zero initialized. Fill rest of its fields:
	hookfunc->idString = (idString) ? strdup(idString) : strdup("");
	hookfunc->hookfunctype = hookfunctype;
	
	// Return pointer to new hook slot:
	return(hookfunc);
}

/* PsychPipelibneDisableHook - Disable named hook chain. */
void PsychPipelineDisableHook(PsychWindowRecordType *windowRecord, const char* hookString)
{
	int hook=PsychGetHookByName(hookString);
	if (hook==-1) PsychErrorExitMsg(PsychError_user, "DisableHook: Unknown (non-existent) hook name provided.");
	windowRecord->HookChainEnabled[hook] = FALSE;
	return;
}

/* PsychPipelibneEnableHook - Enable named hook chain. */
void PsychPipelineEnableHook(PsychWindowRecordType *windowRecord, const char* hookString)
{
	int hook=PsychGetHookByName(hookString);
	if (hook==-1) PsychErrorExitMsg(PsychError_user, "EnableHook: Unknown (non-existent) hook name provided.");
	windowRecord->HookChainEnabled[hook] = TRUE;
	return;
}

/* PsychPipelineResetHook() - Reset named hook chain by deleting all assigned callback functions.
 * windowRecord - Window/Texture for which processing chain should be reset.
 * hookString   - Name string of hook chain to reset.
 */
void PsychPipelineResetHook(PsychWindowRecordType *windowRecord, const char* hookString)
{
	PtrPsychHookFunction hookfunc, hookiter;
	int hookidx=PsychGetHookByName(hookString);
	if (hookidx==-1) PsychErrorExitMsg(PsychError_user, "ResetHook: Unknown (non-existent) hook name provided.");
	hookiter = windowRecord->HookChain[hookidx]; 
	while(hookiter) {
			hookfunc = hookiter;
			hookiter = hookiter->next;
			// Delete all referenced memory:
			free(hookfunc->idString);
			free(hookfunc->pString1);
			// Delete hookfunc struct itself:
			free(hookfunc);
	}

	// Null-out hook chain:
	windowRecord->HookChain[hookidx]=NULL;
	return;
}

/* PsychPipelineAddShaderToHook()
 * Add a GLSL shader program object to a hook chain. The shader is executed when the corresponding
 * hook chain slot is executed, using the specified blitter and OpenGL context configuration and the
 * specified lut texture bound to unit 1.
 *
 * windowRecord - Execute for this window/texture.
 * hookString   - Attach to this named chain.
 * idString     - Arbitrary name string for identification (query) and debugging.
 * where        - Where to attach (0=Beginning, 1=End).
 * shaderid		- GLSL shader program object id.
 * blitterString - Config string to define the used blitter function and its config.
 * luttexid1	- Id of texture to be bound to first texture unit (unit 1).
 */
void PsychPipelineAddShaderToHook(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, unsigned int shaderid, const char* blitterString, unsigned int luttexid1)
{
	// Create and attach proper preinitialized hook function and return pointer to it for further initialization:
	PtrPsychHookFunction hookfunc = PsychAddNewHookFunction(windowRecord, hookString, idString, where, kPsychShaderFunc);
	// Init remaining fields:
	hookfunc->shaderid =  shaderid;
	hookfunc->pString1 =  (blitterString) ? strdup(blitterString) : strdup("");
	hookfunc->luttexid1 = luttexid1;
	return;
}

/* PsychPipelineAddCFunctionToHook()
 * Add a C callback function to a hook chain. The C callback function is executed when the corresponding
 * hook chain slot is executed, passing a set of parameters via a void* struct pointer. The set of parameters
 * is dependent on the exact hook chain, so users of this function must have knowledge of the PTB-3 source code
 * to know what to expect.
 *
 * windowRecord - Execute for this window/texture.
 * hookString   - Attach to this named chain.
 * idString     - Arbitrary name string for identification (query) and debugging.
 * where        - Where to attach (0=Beginning, 1=End).
 * procPtr		- A void* function pointer which will be cast to a proper function pointer during execution.
 */
void PsychPipelineAddCFunctionToHook(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, void* procPtr)
{
	// Create and attach proper preinitialized hook function and return pointer to it for further initialization:
	PtrPsychHookFunction hookfunc = PsychAddNewHookFunction(windowRecord, hookString, idString, where, kPsychCFunc);
	// Init remaining fields:
	hookfunc->cprocfunc =  procPtr;
	return;
}

/* PsychPipelineAddRuntimeFunctionToHook()
 * Add a runtime environment callback function to a hook chain. The function is executed when the corresponding
 * hook chain slot is executed, passing a set of parameters. The set of parameters depends on the exact hook
 * chain, so users of this function must have knowledge of the PTB-3 source code to know what to expect.
 *
 * The mechanism to execute a runtime function depends on the runtime environment. On Matlab and Octave, the
 * internal feval() functions are used to call a Matlab- or Octave function, either a builtin function or some
 * function defined as M-File or dynamically linked.
 *
 * windowRecord - Execute for this window/texture.
 * hookString   - Attach to this named chain.
 * idString     - Arbitrary name string for identification (query) and debugging.
 * where        - Where to attach (0=Beginning, 1=End).
 * evalString	- A function string to be passed to the runtime environment for evaluation during execution.
 */
void PsychPipelineAddRuntimeFunctionToHook(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, const char* evalString)
{
	// Create and attach proper preinitialized hook function and return pointer to it for further initialization:
	PtrPsychHookFunction hookfunc = PsychAddNewHookFunction(windowRecord, hookString, idString, where, kPsychMFunc);
	// Init remaining fields:
	hookfunc->pString1 =  (evalString) ? strdup(evalString) : strdup("");
	return;
}

/* PsychPipelineAddBuiltinFunctionToHook()
 * Add a builtin callback function to a hook chain. The function is executed when the corresponding
 * hook chain slot is executed, passing a set of parameters. The set of parameters depends on the exact hook
 * chain, so users of this function must have knowledge of the PTB-3 source code to know what to expect.
 *
 * windowRecord - Execute for this window/texture.
 * hookString   - Attach to this named chain.
 * idString     - This idString defines the builtin function to call.
 * where        - Where to attach (0=Beginning, 1=End).
 * configString	- A string with configuration parameters.
 */
void PsychPipelineAddBuiltinFunctionToHook(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, const char* configString)
{
	// Create and attach proper preinitialized hook function and return pointer to it for further initialization:
	PtrPsychHookFunction hookfunc = PsychAddNewHookFunction(windowRecord, hookString, idString, where, kPsychBuiltinFunc);
	// Init remaining fields:
	hookfunc->pString1 =  (configString) ? strdup(configString) : strdup("");
	return;
}

/* PsychPipelineDeleteHookSlot
 * Remove slot at index 'slotid', move all slots past this slot one slot up.
 */
void PsychPipelineDeleteHookSlot(PsychWindowRecordType *windowRecord, const char* hookString, int slotid)
{
	PtrPsychHookFunction hookfunc;
	PtrPsychHookFunction *prehookfunc;
	int idx;
	int hookidx=PsychGetHookByName(hookString);
	if (hookidx==-1) PsychErrorExitMsg(PsychError_user, "RemoveHook: Unknown (non-existent) hook name provided.");

	// Perform linear search until proper slot reached or proper name reached:
	idx=0;	
	hookfunc = windowRecord->HookChain[hookidx]; 
	prehookfunc = &(windowRecord->HookChain[hookidx]);
	while(hookfunc && idx<slotid) {
		prehookfunc = &(hookfunc->next);
		hookfunc = hookfunc->next;
		idx++;
	}

	// Anything found?
	if(hookfunc == NULL) return;
	
	// hookfunc is to be deleted.
	
	// Detach it from hookchain, update predecessors next pointer so it points to successor:
	*prehookfunc = hookfunc->next;
	
	// Detached. Delete hookfunc:
	free(hookfunc->pString1);
	free(hookfunc->idString);
	free(hookfunc);
	hookfunc = NULL;
	
	// Done.
	return;
}

/* PsychPipelineQueryHookSlot
 * Query properties of a specific hook slot in a specific hook chain of a specific window:
 * windowRecord - Query for this window/texture.
 * hookString   - Query this named chain.
 * idString     - This string defines the specific slot to query. Can contain an integral number, then the associated slot is
 *				  queried, or a idString (as assigned during creation), then a slot with that name is queried. Partial name
 *				  matches are also accepted to search for substrings...
 *
 * On successfull return, the following values are assigned, on unsuccessfull return (=-1), nothing is assigned:
 * insertString = The subcommand one would need to issue to (re-)add this slot at exactly the same place. E.g., 'InsertAt4MFunction' if slot is a M-Function at position 4.
 * idString = The name string of this slot *Read-Only*
 * blitterString = Config string for this slot.
 * doubleptr = Double encoded void* to the C-Callback function, if any.
 * shaderid = Double encoded shader handle for GLSL shader, if any.
 * luttexid = Double encoded lut texture handle for unit 1, if any.
 *
 * The return value contains the hook slot index where the function was found, or -1 if no matching function could be found. 
 */
int PsychPipelineQueryHookSlot(PsychWindowRecordType *windowRecord, const char* hookString, char** insertString, char** idString, char** blitterString, double* doubleptr, double* shaderid, double* luttexid1)
{
	PtrPsychHookFunction hookfunc;
	char myinsertString[100];
	char mytypeString[100];
	int targetidx, idx;
	int nrassigned = sscanf((*idString), "%i", &targetidx);
	int hookidx=PsychGetHookByName(hookString);
	if (hookidx==-1) PsychErrorExitMsg(PsychError_user, "QueryHook: Unknown (non-existent) hook name provided.");
	if (nrassigned != 1) targetidx=-1;
	idx=0;
	
	// Perform linear search until proper slot reached or proper name reached:
	hookfunc = windowRecord->HookChain[hookidx]; 
	while(hookfunc && ((targetidx>-1 && idx<targetidx) || (targetidx==-1 && strstr(hookfunc->idString, *idString)==NULL))) {
			hookfunc = hookfunc->next;
			idx++;
	}
	
	// If hookfunc is non-NULL then we found our slot:
	if (hookfunc==NULL) {
		*insertString = NULL;
		*idString = NULL;
		*blitterString = NULL;
		*doubleptr = 0;
		*shaderid = 0;
		*luttexid1 = 0;
		return(-1);
	}
	
	switch(hookfunc->hookfunctype) {
		case kPsychBuiltinFunc:
			sprintf(mytypeString, "Builtin");
		break;
		
		case kPsychMFunc:
			sprintf(mytypeString, "MFunction");
		break;

		case kPsychCFunc:
			sprintf(mytypeString, "CFunction");
		break;

		case kPsychShaderFunc:
			sprintf(mytypeString, "Shader");
		break;

		default:
			PsychErrorExitMsg(PsychError_internal, "Unknown id for hookfunction in PsychPipelineQueryHookSlot()! Update the code and recompile!!");
	}
	
	sprintf(myinsertString, "InsertAt%i%s", idx, mytypeString);
	*insertString = strdup(myinsertString);
	*idString = hookfunc->idString;
	*blitterString = hookfunc->pString1;
	*doubleptr = PsychPtrToDouble(hookfunc->cprocfunc);
	*shaderid = (double) hookfunc->shaderid;
	*luttexid1= (double) hookfunc->luttexid1;
	return(idx);
}

/* PsychPipelineDumpHook
 * Dump properties of a specific hook chain of a specific window in a human-readable format into
 * the console of the scripting environment:
 *
 * windowRecord - Query for this window/texture.
 * hookString   - Query this named chain.
 */
void PsychPipelineDumpHook(PsychWindowRecordType *windowRecord, const char* hookString)
{
	PtrPsychHookFunction hookfunc;
	int i=0;
	int hookidx=PsychGetHookByName(hookString);
	if (hookidx==-1) PsychErrorExitMsg(PsychError_user, "DumpHook: Unknown (non-existent) hook name provided.");
	
	hookfunc = windowRecord->HookChain[hookidx];
	printf("Hook chain %s is currently %s.\n", hookString, (windowRecord->HookChainEnabled[hookidx]) ? "enabled" : "disabled");
	if (hookfunc==NULL) {
		printf("No processing assigned to this hook-chain.\n");
	}
	else {
		printf("Following hook slots are assigned to this hook-chain:\n");
		printf("=====================================================\n");
	}
	
	while(hookfunc) {
		printf("Slot %i: Id='%s' : ", i, hookfunc->idString);
		switch(hookfunc->hookfunctype) {
			case kPsychShaderFunc:
				printf("GLSL-Shader      : id=%i , luttex1=%i , blitter=%s\n", hookfunc->shaderid, hookfunc->luttexid1, hookfunc->pString1);
			break;
			
			case kPsychCFunc:
				printf("C-Callback       : void*= %p\n", hookfunc->cprocfunc);
			break;

			case kPsychMFunc:
				printf("Runtime-Function : Evalstring= %s\n", hookfunc->pString1);
			break;
			
			case kPsychBuiltinFunc:
				printf("Builtin-Function : Name= %s\n", hookfunc->idString);
			break;
		}
		
		// Next one, if any:
		i++;
		hookfunc = hookfunc->next;
	}
	
	printf("=====================================================\n\n");
	fflush(NULL);
	return;
}

/* PsychPipelineDumpAllHooks
 * Dump current state of all hook-points for given window. See PsychPipelineDumpHook()
 * for more info.
 */
void PsychPipelineDumpAllHooks(PsychWindowRecordType *windowRecord)
{
	int i;
	for (i=0; i<MAX_SCREEN_HOOKS; i++) {
		PsychPipelineDumpHook(windowRecord, PsychHookPointNames[i]);
	}
	return;
}

psych_bool PsychIsHookChainOperational(PsychWindowRecordType *windowRecord, int hookid)
{
	// Child protection:
	if (hookid<0 || hookid>=MAX_SCREEN_HOOKS) PsychErrorExitMsg(PsychError_internal, "In PsychIsHookChainOperational: Was asked to check unknown (non-existent) hook chain with invalid id!");

	// Hook chain enabled for processing and contains at least one hook slot?
	if ((!windowRecord->HookChainEnabled[hookid]) || (windowRecord->HookChain[hookid] == NULL)) {
		// Chain is empty or disabled.
		return(FALSE);
	}
	
	// Chain operational:
	return(TRUE);
}

/* PsychPipelineExecuteHook()
 * Execute the full hook processing chain for a specific hook and a specific windowRecord.
 * This checks if the chain is enabled. If it isn't enabled, it skips processing.
 * If it is enabled, it iterates over the full chain, executes all assigned hook functions in order and uses the FBO's between minfbo and maxfbo
 * as pingpong buffers if neccessary.
 */
psych_bool PsychPipelineExecuteHook(PsychWindowRecordType *windowRecord, int hookId, void* hookUserData, void* hookBlitterFunction, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo)
{
	PtrPsychHookFunction hookfunc;
	int i=0;
	int pendingFBOpingpongs = 0;
	PsychFBO *mysrcfbo1, *mysrcfbo2, *mydstfbo, *mynxtfbo;
	PsychFBO **bouncefbo2;
	psych_bool gfxprocessing;
	GLint restorefboid = 0;
	psych_bool scissor_ignore = FALSE;
	psych_bool scissor_enabled = FALSE;
	int sciss_x, sciss_y, sciss_w, sciss_h;

	// Child protection:
	if (hookId<0 || hookId>=MAX_SCREEN_HOOKS) PsychErrorExitMsg(PsychError_internal, "In PsychPipelineExecuteHook: Was asked to execute unknown (non-existent) hook chain with invalid id!");

	// Hook chain enabled for processing and contains at least one hook slot?
	if ((!windowRecord->HookChainEnabled[hookId]) || (windowRecord->HookChain[hookId] == NULL)) {
		// Chain is empty or disabled.
		return(TRUE);
	}
	
	// Is this an image processing hook?
	gfxprocessing = (dstfbo!=NULL) ? TRUE : FALSE;

	// Get start of enabled chain:
	hookfunc = windowRecord->HookChain[hookId];

	// Count number of needed ping-pong FBO switches inside this chain:
	while(hookfunc) {
		// Pingpong command?
		if (hookfunc->hookfunctype == kPsychBuiltinFunc && strcmp(hookfunc->idString, "Builtin:FlipFBOs")==0) pendingFBOpingpongs++;
		// Process next hookfunc slot in chain, if any:
		hookfunc = hookfunc->next;
	}

	if (gfxprocessing) {
		// Prepare gfx-processing:
		
		// Slightly ugly, because it is a layering violation, but kind'a unavoidable,
		// as the final output formatting blit chain is a special case...
		// Assign bouncefbo2 if it is available and needed, otherwise assign dstfbo:
		// It is available if preConversionFBO[3]>=0. It is needed if pendingFBOpingpongs > 1
		// and this is one of the final output formatting blit chains. Other chains don't need
		// this special treatment with a 2nd bounce buffer:
		if ((windowRecord->preConversionFBO[3]>=0) && (pendingFBOpingpongs > 1) &&
			(hookId == kPsychFinalOutputFormattingBlit || hookId == kPsychFinalOutputFormattingBlit0 || hookId == kPsychFinalOutputFormattingBlit1)) {
				// Assign special 2nd bounce buffer:
				bouncefbo2 = &(windowRecord->fboTable[windowRecord->preConversionFBO[3]]);
		}
		else {
				// Standard case, dstfbo acts as final target and as 2nd bounce buffer if needed for multi-slot chains:
				bouncefbo2 = dstfbo;
		}

		// Backup scissoring state:
		scissor_enabled = glIsEnabled(GL_SCISSOR_TEST);
		
		// If this is a multi-pass chain we'll need a bounce buffer FBO:
		if ((pendingFBOpingpongs > 0 && bouncefbo == NULL) || (pendingFBOpingpongs > 1 && ((*bouncefbo2)->fboid == 0))) {
			printf("PTB-ERROR: Hook processing chain '%s' is a multi-pass processing chain with %i passes,\n", PsychHookPointNames[hookId], pendingFBOpingpongs + 1);
			printf("PTB-ERROR: but imaging pipeline is not configured for multi-pass processing! You need to supply the additional flag\n");
			printf("PTB-ERROR: %s as imagingmode to Screen('OpenWindow') to tell PTB about these requirements and then restart.\n\n",
					(pendingFBOpingpongs > 1) ? "kPsychNeedMultiPass" : "kPsychNeedDualPass");
			// Ok, abort...
			PsychErrorExitMsg(PsychError_user, "Insufficient pipeline configuration for processing. Adapt the 'imagingmode' flag according to my tips!");
		}
		
		if ((pendingFBOpingpongs % 2) == 0) {
			// Even number of ping-pongs needed in this chain. We stream from source fbo to
			// destination fbo in first pass.
			mysrcfbo1 = (srcfbo1) ? *srcfbo1 : NULL;
			mysrcfbo2 = (srcfbo2) ? *srcfbo2 : NULL;
			mydstfbo  = *bouncefbo2;
			mynxtfbo  = (bouncefbo) ? *bouncefbo : NULL;
		}
		else {
			// Odd number of ping-pongs needed. Initially stream from source to bouncefbo:
			mysrcfbo1 = (srcfbo1) ? *srcfbo1 : NULL;
			mysrcfbo2 = (srcfbo2) ? *srcfbo2 : NULL;
			mydstfbo  = (bouncefbo) ? *bouncefbo : NULL;
			mynxtfbo  = *bouncefbo2;
		}

		// Special case: If this is the last processing slot, aka pendingFBOpingpongs == 0,
		// then mydstfbo must be our real destination framebuffer:
		if (pendingFBOpingpongs == 0) mydstfbo = *dstfbo;

		// Enable associated GL context:
		PsychSetGLContext(windowRecord);
		
		// Save current FBO bindings for later restore on classic desktop OpenGL1/2:
		if (glBindFramebufferEXT && PsychIsGLClassic(windowRecord)) {
			glGetIntegerv(GL_FRAMEBUFFER_BINDING_EXT, &restorefboid);
		}
        else {
            // Unsupported on OES and next gen: TODO FIXME - Implement proper state tracking.
            restorefboid = 0;
        }
		
		// Setup initial source -> target binding:
		PsychPipelineSetupRenderFlow(mysrcfbo1, mysrcfbo2, mydstfbo, scissor_ignore);
	}

	// Reget start of enabled chain:
	hookfunc = windowRecord->HookChain[hookId];

	// Iterate over all slots:
	while(hookfunc) {
		// Debug output, if requested:
		if (PsychPrefStateGet_Verbosity()>4) {
			printf("Hookchain '%s' : Slot %i: Id='%s' : ", PsychHookPointNames[hookId], i, hookfunc->idString);
			switch(hookfunc->hookfunctype) {
				case kPsychShaderFunc:
					printf("GLSL-Shader      : id=%i , luttex1=%i , blitter=%s\n", hookfunc->shaderid, hookfunc->luttexid1, hookfunc->pString1);
					break;
					
				case kPsychCFunc:
					printf("C-Callback       : void*= %p\n", hookfunc->cprocfunc);
					break;
					
				case kPsychMFunc:
					printf("Runtime-Function : Evalstring= %s\n", hookfunc->pString1);
					break;
					
				case kPsychBuiltinFunc:
					printf("Builtin-Function : Name= %s : Params= %s\n", hookfunc->idString, hookfunc->pString1);
					break;
			}
		}
		
		// Is this a ping-pong command?
		if (hookfunc->hookfunctype == kPsychBuiltinFunc && strcmp(hookfunc->idString, "Builtin:FlipFBOs")==0) {
			// Ping pong buffer swap requested:
			pendingFBOpingpongs--;
			mysrcfbo1 = mydstfbo;
			mydstfbo  = mynxtfbo;
			if ((pendingFBOpingpongs % 2) == 0) {
				// Even number of ping-pongs remaining in this chain.
				mynxtfbo  = (bouncefbo) ? *bouncefbo : NULL;
			}
			else {
				// Odd number of ping-pongs remaining.
				mynxtfbo  = *bouncefbo2;
			}
			
			// Special case: If this is the last processing slot, aka pendingFBOpingpongs == 0,
			// then mydstfbo must be our real destination framebuffer:
			if (pendingFBOpingpongs == 0) mydstfbo = *dstfbo;
			
			if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: SWAPPING PING-PONG FBOS, %i swaps pending...\n", pendingFBOpingpongs);
			
			// Set new src -> dst binding:
			PsychPipelineSetupRenderFlow(mysrcfbo1, mysrcfbo2, mydstfbo, scissor_ignore);
		}
		else {
			// Restricted area processing?
			if (hookfunc->hookfunctype == kPsychBuiltinFunc && strstr(hookfunc->idString, "Builtin:RestrictToScissorROI")) {
				// Restrict pixel processing to specified region of interest ROI by setting
				// up a proper scissor rectangle and enabling scissor tests. The special
				// ROI (-1,-1,-1,-1) means: Disable scissor testing -> Unrestrict. 
				if (4!=sscanf(hookfunc->pString1, "%i:%i:%i:%i", &sciss_x, &sciss_y, &sciss_w, &sciss_h)) {
					if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: In PsychPipelineExecuteHook: Builtin:RestrictToScissorROI - Parameter parse error in string %s\n", hookfunc->idString);
					return(FALSE);
				}
				
				if (sciss_x==-1 && sciss_y==-1 && sciss_w==-1 && sciss_h==-1) {
					// Disable scissor tests:
					glDisable(GL_SCISSOR_TEST);
					scissor_ignore = FALSE;
				}
				else {
					// Setup and enable scissor test:
					glEnable(GL_SCISSOR_TEST);
					glScissor(sciss_x, sciss_y, sciss_w, sciss_h);
					// Make sure PsychSetupRenderFlow() ignores scissor setup:
					scissor_ignore = TRUE;
				}
			}
			else {
				// Normal hook function - Process this hook function:
				if (!PsychPipelineExecuteHookSlot(windowRecord, hookId, hookfunc, hookUserData, hookBlitterFunction, srcIsReadonly, allowFBOSwizzle, &mysrcfbo1, &mysrcfbo2, &mydstfbo, &mynxtfbo)) {
					// Failed!
					if (PsychPrefStateGet_Verbosity()>0) {
						printf("PTB-ERROR: Failed in processing of Hookchain '%s' : Slot %i: Id='%s'  --> Aborting chain processing. Set verbosity to 5 for extended debug output.\n", PsychHookPointNames[hookId], i, hookfunc->idString);
					}
					return(FALSE);
				}
			}
		}
		
		// Process next hookfunc slot in chain, if any:
		i++;
		hookfunc = hookfunc->next;
	}

	if (gfxprocessing) {
		// Disable renderflow:
		PsychPipelineSetupRenderFlow(NULL, NULL, NULL, scissor_ignore);
		
		// Restore old FBO bindings:
		if (glBindFramebufferEXT) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, restorefboid);
		
		// Restore scissor state:
		if (scissor_enabled) {
			glEnable(GL_SCISSOR_TEST);
		}
		else {
			glDisable(GL_SCISSOR_TEST);
		}
	}

	// Done.
	return(TRUE);
}

/* PsychPipelineProcessMacros()
 * Expand given cmdString, based on info in given windowRecord, return expanded string.
 * The expanded string is returned as a char* and the caller is responsible for freeing
 * the allocated memory.
 *
 * Returns 1 on succcess, 0 on error.
 */
int PsychPipelineProcessMacros(PsychWindowRecordType *windowRecord, char* cmdString)
{
	psych_bool repeatit = TRUE;
	char  varName[256];
	char* pCurrent		= cmdString;
	char* pToken		= NULL;
	int   i;
	
	double* dbltarget	= NULL;
	PsychGenericScriptType*	newvar = NULL;
	
	while ((repeatit) && (pCurrent) && (pCurrent[0]!=0)) {
		// Reset repeatit. Will be set by following code again if another parse-iteration is needed:
		repeatit	= FALSE;
		dbltarget	= NULL;
		newvar		= NULL;

		if (PsychPrefStateGet_Verbosity() > 10) printf("PTB-DEBUG: PsychPipelineProcessMacros: cmd = %s\n", pCurrent); 
		
		// Check current remaining command string for certain keywords to process:
		// -----------------------------------------------------------------------
		
		// Return of current gammatable requested?
		sprintf(varName, "IMAGINGPIPE_GAMMATABLE");
		if ((pToken = strstr(pCurrent, varName))) {
			// Wants us to assign current gammatable from Screen('LoadNormalizedGammatable', ..., 2) call
			// to the variable IMAGINGPIPE_GAMMATABLE:
			
			// Advance to position behind the varName:
			pCurrent = &pToken[strlen(varName)];
			repeatit = TRUE;
			
			// Any pending gammatable stored internally for update? If not, we skip further processing:
			if (windowRecord->inRedTable == NULL) {
				if (PsychPrefStateGet_Verbosity() > 10) printf("PTB-DEBUG: PsychPipelineProcessMacros: IMAGINGPIPE_GAMMATABLE variable assignment failed. No suitable CLUT set in Screen('LoadNormalizedGammaTable'). Skipped slot!\n");

				// Return error code zero to abort processing of this hook slot:
				return(0);
			}
			
			if (windowRecord->inTableSize < 1) {
				if (PsychPrefStateGet_Verbosity() > 10) printf("PTB-DEBUG: PsychPipelineProcessMacros: IMAGINGPIPE_GAMMATABLE variable assignment failed. CLUT has less than the required 1 slots. Skipped slot!\n");

				// Return error code zero to abort processing of this hook slot:
				return(0);
			}
			
			// Allocate runtime double matrix of sufficient size: 'newvar' is the handle of it,
			// dbltarget is a pointer to the target double matrix we need to fill with the clut:
			PsychAllocateNativeDoubleMat(windowRecord->inTableSize, 3, 1, &dbltarget, &newvar);
			
			for (i = 0; i < windowRecord->inTableSize; i++) *(dbltarget++) = (double) windowRecord->inRedTable[i];
			for (i = 0; i < windowRecord->inTableSize; i++) *(dbltarget++) = (double) windowRecord->inGreenTable[i];
			for (i = 0; i < windowRecord->inTableSize; i++) *(dbltarget++) = (double) windowRecord->inBlueTable[i];
			
			// Release the gamma table:
			free(windowRecord->inRedTable); windowRecord->inRedTable = NULL;
			free(windowRecord->inGreenTable); windowRecord->inGreenTable = NULL;
			free(windowRecord->inBlueTable); windowRecord->inBlueTable = NULL;
			windowRecord->inTableSize = 0;
			windowRecord->loadGammaTableOnNextFlip = 0;
			
			// Copy to base workspace:
			// N.B. caller workspace would be cleaner, but there's an incompatibility between Octave and
			// Matlab which causes trouble if you do so, so stick to "base" for sanity of mind!
			if (PsychRuntimePutVariable("base", varName, newvar)) {
				if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychPipelineProcessMacros: IMAGINGPIPE_GAMMATABLE variable assignment failed in runtime! Skipped slot!\n");
				return(0);
			}
			else {
				if (PsychPrefStateGet_Verbosity() > 10) printf("PTB-DEBUG: PsychPipelineProcessMacros: IMAGINGPIPE_GAMMATABLE variable assignment success in runtime!\n");
			}

			// Next parseloop iteration:
			continue;
		}

		// Return of current flipcount requested?
		sprintf(varName, "IMAGINGPIPE_FLIPCOUNT");
		if ((pToken = strstr(pCurrent, varName))) {
			// Wants us to assign current flipCount to the variable IMAGINGPIPE_FLIPCOUNT:
			
			// Advance to position behind the varName:
			pCurrent = &pToken[strlen(varName)];
			repeatit = TRUE;

			// Allocate runtime double matrix of scalar size: 'newvar' is the handle of it,
			// dbltarget is a pointer to the target double matrix we need to fill with the clut:
			PsychAllocateNativeDoubleMat(1, 1, 1, &dbltarget, &newvar);
			*dbltarget = (double) windowRecord->flipCount;

			// Copy to base workspace:
			// N.B. caller workspace would be cleaner, but there's an incompatibility between Octave and
			// Matlab which causes trouble if you do so, so stick to "base" for sanity of mind!
			if (PsychRuntimePutVariable("base", varName, newvar)) {
				if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychPipelineProcessMacros: IMAGINGPIPE_FLIPCOUNT variable assignment failed in runtime! Skipped slot!\n");
				return(0);
			}
			else {
				if (PsychPrefStateGet_Verbosity() > 10) printf("PTB-DEBUG: PsychPipelineProcessMacros: ASSIGNED %f\n", (float) *dbltarget); 
			}

			// Next parseloop iteration:
			continue;			
		}
		
		// Next parse loop iteration, if any...
	}	// Parse loop.
	
	// Successfully processed string:
	return(1);
}

/* PsychPipelineExecuteHookSlot()
 * Execute a single hookfunction slot in a hook chain for a specific window.
 */
psych_bool PsychPipelineExecuteHookSlot(PsychWindowRecordType *windowRecord, int hookId, PsychHookFunction* hookfunc, void* hookUserData, void* hookBlitterFunction, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo)
{
	psych_bool dispatched = FALSE;
	char* execString = NULL;

	// Dispatch by hook function type:
	switch(hookfunc->hookfunctype) {
		case kPsychShaderFunc:
			// Call a GLSL shader to do some image processing: We just execute the blitter, the shader gets assigned inside
			// this function.
			if (!PsychPipelineExecuteBlitter(windowRecord, hookfunc, hookUserData, hookBlitterFunction, srcIsReadonly, allowFBOSwizzle, srcfbo1, srcfbo2, dstfbo, bouncefbo)) {
				// Blitter failed!
				return(FALSE);
			}
			dispatched=TRUE;
		break;
			
		case kPsychCFunc:
			// Call a C callback function via the given memory function pointer:
			printf("TODO: EXECUTE -- C-Callback       : void*= %p\n", hookfunc->cprocfunc);
			dispatched=TRUE;
		break;
			
		case kPsychMFunc:
			// Call the eval() function of our scripting runtime environment to evaluate
			// function string pString1. Currently supported are Matlab & Octave, so this
			// can be the call string of an arbitrary Matlab/Octave builtin or M-Function.
			// Care has to be taken that the called functions do not invoke any Screen
			// subfunctions! Screen is not reentrant, so that would likely screw seriously!

			// Process certain macro keywords inside pString1 before calling the pString itself:
			if (PsychPipelineProcessMacros(windowRecord, hookfunc->pString1)) {
				// Call the runtime environment to process the command string:
				if (PsychRuntimeEvaluateString(hookfunc->pString1)) {
				}
			}
			dispatched=TRUE;
		break;
			
		case kPsychBuiltinFunc:
			// Dispatch to a builtin function:

			if (strcmp(hookfunc->idString, "Builtin:FlipFBOs")==0) { dispatched=TRUE; } // No op here. Done in upper layer...
			if (strstr(hookfunc->idString, "Builtin:RestrictToScissorROI")) { dispatched=TRUE; } // No op here. Done in upper layer...

			if (strstr(hookfunc->idString, "Builtin:IdentityBlit")) {
				// Perform the most simple blit operation: A simple one-to-one copy of input FBO to output FBO:
				if (!PsychPipelineExecuteBlitter(windowRecord, hookfunc, hookUserData, hookBlitterFunction, TRUE, FALSE, srcfbo1, NULL, dstfbo, NULL)) {
					// Blitter failed!
					return(FALSE);
				}
				dispatched=TRUE;
			}

			if (strcmp(hookfunc->idString, "Builtin:RenderClutBits++")==0) {
				// Compute the T-Lock encoded CLUT for Cambridge Research Bits++ system in Bits++ mode. The CLUT
				// is set via the standard Screen('LoadNormalizedGammaTable', ..., loadOnNextFlip) call by setting
				// loadOnNextFlip to a value of 2.
				if (!PsychPipelineBuiltinRenderClutBitsPlusPlus(windowRecord, hookfunc)) {
					// Operation failed!
					return(FALSE);
				}
				dispatched=TRUE;
			}

			if (strcmp(hookfunc->idString, "Builtin:RenderClutViaRuntime")==0) {
				// Pass the last CLUT that was set via the standard Screen('LoadNormalizedGammaTable', ..., loadOnNextFlip) call by setting
				// loadOnNextFlip to a value of 2 to the runtime environment.
				if (!PsychPipelineBuiltinRenderClutViaRuntime(windowRecord, hookfunc)) {
					// Operation failed!
					return(FALSE);
				}
				dispatched=TRUE;
			}
			
			if (strcmp(hookfunc->idString, "Builtin:RenderStereoSyncLine")==0) {
				// Draw a blue-line-sync sync line at the bottom of the current framebuffer. This is needed
				// to drive stereo shutter glasses with blueline-sync in quad-buffered frame-sequential stereo
				// mode.
				if (!PsychPipelineBuiltinRenderStereoSyncLine(windowRecord, hookId, hookfunc)) {
					// Operation failed!
					return(FALSE);
				}
				dispatched=TRUE;
			}

		break;
			
		default:
			PsychErrorExitMsg(PsychError_internal, "In PsychPipelineExecuteHookSlot: Was asked to execute unknown (non-existent) hook function type!");
	}
	
	if (!dispatched) {
		if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: Failed to dispatch hook slot - Unknown command or failure in command execution.\n");
		return(FALSE);
	}
	
	return(TRUE);
}

void PsychPipelineSetupRenderFlow(PsychFBO* srcfbo1, PsychFBO* srcfbo2, PsychFBO* dstfbo, psych_bool scissor_ignore)
{
	static int ow=0;
	static int oh=0;
	int w, h;

	// Select rendertarget:
	if (glBindFramebufferEXT) glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, (dstfbo) ? dstfbo->fboid : 0);

	// Assign color texture of srcfbo2, if any,  to texture unit 1:
	glActiveTextureARB(GL_TEXTURE1_ARB);
	if (srcfbo2) {
		// srcfbo2 is valid: Assign its color buffer texture:
		if (PsychPrefStateGet_Verbosity()>4) printf("TexUnit 1 reading from texid -- %i\n", srcfbo2->coltexid);
		glBindTexture(srcfbo2->textarget, srcfbo2->coltexid);
		// Set texture application mode to replace:
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		
		glTexParameteri(srcfbo2->textarget, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(srcfbo2->textarget, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(srcfbo2->textarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(srcfbo2->textarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_TEXTURE_RECTANGLE_EXT);
		glEnable(srcfbo2->textarget);
	}
	else {
		// srcfbo2 doesn't exist: Unbind and deactivate 2nd unit:
		if (PsychPrefStateGet_Verbosity() > 10) printf("TexUnit 1 not reading from srcfbo2\n");
		glBindTexture(GL_TEXTURE_RECTANGLE_EXT, 0);
		glDisable(GL_TEXTURE_RECTANGLE_EXT);
		glBindTexture(GL_TEXTURE_2D, 0);
		glDisable(GL_TEXTURE_2D);
	}

	// Assign color texture of srcfbo1 to texture unit 0:
	glActiveTextureARB(GL_TEXTURE0_ARB);
	if (srcfbo1) {
		// srcfbo1 is valid: Assign its color buffer texture:
		if (PsychPrefStateGet_Verbosity()>4) printf("TexUnit 0 reading from texid -- %i\n", srcfbo1->coltexid);
		glBindTexture(srcfbo1->textarget, srcfbo1->coltexid);

		// Set texture application mode to replace:
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		
		glTexParameteri(srcfbo1->textarget, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(srcfbo1->textarget, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(srcfbo1->textarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(srcfbo1->textarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_TEXTURE_RECTANGLE_EXT);
		glEnable(srcfbo1->textarget);
	}
	else {
		// srcfbo1 doesn't exist: Unbind and deactivate 1st unit:
		if (PsychPrefStateGet_Verbosity() > 10) printf("TexUnit 0 not reading from srcfbo1\n");
		glBindTexture(GL_TEXTURE_RECTANGLE_EXT, 0);
		glDisable(GL_TEXTURE_RECTANGLE_EXT);
		glBindTexture(GL_TEXTURE_2D, 0);
		glDisable(GL_TEXTURE_2D);
	}

	if (dstfbo) {
		// Setup viewport, scissor rectangle and projection matrix for orthonormal rendering into the
		// target FBO or system framebuffer:
		w = (int) dstfbo->width;
		h = (int) dstfbo->height;
		if (PsychPrefStateGet_Verbosity()>4) {
			if (dstfbo->fboid > 0) {
				printf("Blitting to Targettex = %i , w x h = %i %i\n", dstfbo->coltexid, w, h);
			} else {
				printf("Blitting to system framebuffer, w x h = %i %i\n", w, h);
			}
		}

		// Settings changed? We skip if not - state changes are expensive...
		if ((w!=ow || h!=oh) && PsychIsMasterThread()) {
			ow=w;
			oh=h;
			
			// Setup viewport and scissor for full FBO area:
			glViewport(0, 0, w, h);
			if (!scissor_ignore) glScissor(0, 0, w, h);
			
			// Setup projection matrix for a proper orthonormal projection for this framebuffer:
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
            
            if (!PsychIsGLES(NULL)) {
                gluOrtho2D(0, w, h, 0);
            }
            else {
                glOrthofOES(0, (float) w, (float) h, 0, -1, 1);
            }

			// Switch back to modelview matrix, but leave it unaltered:
			glMatrixMode(GL_MODELVIEW);
		}
	}
	else {
        if (PsychIsMasterThread()) {
            // Reset our cached settings:
            ow=0;
            oh=0;
        }
	}

	return;
}

psych_bool PsychPipelineExecuteBlitter(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc, void* hookUserData, void* hookBlitterFunction, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo)
{
	psych_bool rc = TRUE;
	PsychBlitterFunc blitterfnc = NULL;
	GLenum glerr;
	char*  pstrpos = NULL;
	int texunit, texid;
	
	// Select proper blitter function:
	
	// Initialize with master blitter function (if any). If none set,
	// this will init to NULL:
	blitterfnc = hookBlitterFunction;
	
	// Any special override blitter defined in parameter string?
	if (strstr(hookfunc->pString1, "Blitter:")) {
		// Yes. Which one?
		blitterfnc = NULL;
		
		// Standard blitter? This one does a one-to-one copy without special geometric transformations.
		if (strstr(hookfunc->pString1, "Blitter:IdentityBlit")) blitterfnc = &PsychBlitterIdentity; // Assign our standard one-to-one blitter.
		
		// Displaylist blitter? This one calls an externally setup OpenGL display list to perform complex geometric transformations:
		if (strstr(hookfunc->pString1, "Blitter:DisplayListBlit")) blitterfnc = &PsychBlitterDisplayList;

		// Blitter assigned?
		if (blitterfnc == NULL) {
			if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: Invalid (unknown) blitter specified in blitter string. Blit aborted.\n");
			return(FALSE);
		}
	}
	
	if (blitterfnc == NULL) {
		// No blitter set up to now. Assign the default blitter:
		blitterfnc = &PsychBlitterIdentity; // Assign our standard one-to-one blitter.
	}
	
	// TODO: Common setup code for texturing, filtering, alpha blending, z-test and such...
	
	// Setup code for 1D textures:
	pstrpos = hookfunc->pString1;
	while (pstrpos=strstr(pstrpos, "TEXTURE1D")) {
		if (2==sscanf(pstrpos, "TEXTURE1D(%i)=%i", &texunit, &texid)) {
			glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
			glEnable(GL_TEXTURE_1D);
			glBindTexture(GL_TEXTURE_1D, texid);
			if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Binding gltexid %i to GL_TEXTURE_1D target of texunit %i\n", texid, texunit);
		}
		pstrpos++;
	}

	// Setup code for 2D textures:
	pstrpos = hookfunc->pString1;
	while (pstrpos=strstr(pstrpos, "TEXTURE2D")) {
		if (2==sscanf(pstrpos, "TEXTURE2D(%i)=%i", &texunit, &texid)) {
			glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, texid);
			if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Binding gltexid %i to GL_TEXTURE_2D target of texunit %i\n", texid, texunit);
		}
		pstrpos++;
	}

	// Setup code for 2D rectangle textures:
	pstrpos = hookfunc->pString1;
	while (pstrpos=strstr(pstrpos, "TEXTURERECT2D")) {
		if (2==sscanf(pstrpos, "TEXTURERECT2D(%i)=%i", &texunit, &texid)) {
			glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
			glEnable(GL_TEXTURE_RECTANGLE_EXT);
			glBindTexture(GL_TEXTURE_RECTANGLE_EXT, texid);
			if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Binding gltexid %i to GL_TEXTURE_RECTANGLE_EXT target of texunit %i\n", texid, texunit);
		}
		pstrpos++;
	}

	// Setup code for 3D textures:
	pstrpos = hookfunc->pString1;
	while (pstrpos=strstr(pstrpos, "TEXTURE3D")) {
		if (2==sscanf(pstrpos, "TEXTURE3D(%i)=%i", &texunit, &texid)) {
			glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
			glEnable(GL_TEXTURE_3D);
			glBindTexture(GL_TEXTURE_3D, texid);
			if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: Binding gltexid %i to GL_TEXTURE_3D target of texunit %i\n", texid, texunit);
		}
		pstrpos++;
	}

	glActiveTextureARB(GL_TEXTURE0_ARB);
	
	// Need a shader for this blit op?
	if (hookfunc->shaderid) {
		// Setup shader, if any:
		if (!glUseProgram){
			if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: Blitter invocation failed: Blitter needs to attach GLSL shaders, but shaders are not supported on your hardware!\n");
			rc = FALSE;
		} else {
			// Attach shader:
			if (!(PsychPrefStateGet_ConserveVRAM() & kPsychAvoidCPUGPUSync)) while(glGetError());
			glUseProgram(hookfunc->shaderid);
			if (!(PsychPrefStateGet_ConserveVRAM() & kPsychAvoidCPUGPUSync) && ((glerr = glGetError())!=GL_NO_ERROR)) {
				if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: Blitter invocation failed: glUseProgram(%i) failed with error: %s\n", hookfunc->shaderid, gluErrorString(glerr));
				rc = FALSE;
			}
			
			#if PSYCH_SYSTEM == PSYCH_OSX
			// On OS-X we can query the OS if the bound shader is running on the GPU or if it is running in emulation mode on the CPU.
			// This is an expensive operation - it triggers OpenGL internal state revalidation. Only use for debugging and testing!
			if (PsychPrefStateGet_Verbosity()>5) {
				GLint vsgpu=0, fsgpu=0;
				CGLGetParameter(CGLGetCurrentContext(), kCGLCPGPUVertexProcessing, &vsgpu);
				CGLGetParameter(CGLGetCurrentContext(), kCGLCPGPUFragmentProcessing, &fsgpu);
				printf("PTB-DEBUG: Imaging pipeline GPU shading state: Vertex processing on %s : Fragment processing on %s.\n", (vsgpu) ? "GPU" : "CPU!!", (fsgpu) ? "GPU" : "CPU!!");
			}
			#endif
			
		}
	}
	
	// Execute blitter function:
	rc = (rc && blitterfnc(windowRecord, hookfunc, hookUserData, srcIsReadonly, allowFBOSwizzle, srcfbo1, srcfbo2, dstfbo, bouncefbo));
	if (!rc) {
		if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: Blitter invocation failed: OpenGL error state is: %s\n", gluErrorString(glGetError()));
		while(glGetError());
	}
	
	// TODO: Common teardown code for texturing, filtering and such...

	// Teardown code for 1D textures:
	pstrpos = hookfunc->pString1;
	while (pstrpos=strstr(pstrpos, "TEXTURE1D")) {
		if (2==sscanf(pstrpos, "TEXTURE1D(%i)=%i", &texunit, &texid)) {
			glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
			glBindTexture(GL_TEXTURE_1D, 0);
			glDisable(GL_TEXTURE_1D);
		}
		pstrpos++;
	}

	// Teardown code for 2D textures:
	pstrpos = hookfunc->pString1;
	while (pstrpos=strstr(pstrpos, "TEXTURE2D")) {
		if (2==sscanf(pstrpos, "TEXTURE2D(%i)=%i", &texunit, &texid)) {
			glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
			glBindTexture(GL_TEXTURE_2D, 0);
			glDisable(GL_TEXTURE_2D);
		}
		pstrpos++;
	}

	// Teardown code for 2D rectangle textures:
	pstrpos = hookfunc->pString1;
	while (pstrpos=strstr(pstrpos, "TEXTURERECT2D")) {
		if (2==sscanf(pstrpos, "TEXTURERECT2D(%i)=%i", &texunit, &texid)) {
			glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
			glBindTexture(GL_TEXTURE_RECTANGLE_EXT, 0);
			glDisable(GL_TEXTURE_RECTANGLE_EXT);
		}
		pstrpos++;
	}

	// Teardown code for 3D textures:
	pstrpos = hookfunc->pString1;
	while (pstrpos=strstr(pstrpos, "TEXTURE3D")) {
		if (2==sscanf(pstrpos, "TEXTURE3D(%i)=%i", &texunit, &texid)) {
			glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
			glBindTexture(GL_TEXTURE_3D, 0);
			glDisable(GL_TEXTURE_3D);
		}
		pstrpos++;
	}

	glActiveTextureARB(GL_TEXTURE0_ARB);

	// Reset shader assignment, if any:
	if ((hookfunc->shaderid) && glUseProgram) glUseProgram(0);

	// Return result code:
	return(rc);
}

/* PsychBlitterIdentity()  -- Default blitter.
 *
 * Identity blitter: Blits from srcfbo1 color attachment to dstfbo without geometric transformations or other extras.
 * This is the most common one for one-to-one copies or simple shader image processing. It gets automatically used
 * when no special (non-default) blitter is requested by core code or users blitter parameter string:
 */
psych_bool PsychBlitterIdentity(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc, void* hookUserData, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo)
{
	int w, h, x, y, wf, hf;
	float sx, sy, wt, ht;
	char* strp;
	psych_bool bilinearfiltering;

    // hookUserData, if non-NULL, can provide override parameter string:
    char* pString1 = (hookUserData) ? (char*) hookUserData : hookfunc->pString1;

	// Child protection:
	if (!(srcfbo1 && (*srcfbo1))) {
		PsychErrorExitMsg(PsychError_internal, "In PsychBlitterIdentity(): srcfbo1 is a NULL - Pointer!!!");
	}	

	// Query dimensions of viewport:
	w = (*srcfbo1)->width;
	h = (*srcfbo1)->height;
    
    // Same for texture coordinate space, depending on type of texture in use:
    wt = ((*srcfbo1)->textarget == GL_TEXTURE_2D) ? 1 : (float) w;
    ht = ((*srcfbo1)->textarget == GL_TEXTURE_2D) ? 1 : (float) h;

    // This pot-textures remapping mostly applies to OpenGL-ES 1.x hardware:
    if (((*srcfbo1)->textarget == GL_TEXTURE_2D) && !(windowRecord->gfxcaps & kPsychGfxCapNPOTTex)) {
        // Only power-of-two GL_TEXTURE_2D targets supported. Find real width of
        // FBO color buffer backing texture (wf, hf):
        wf = 1;
        while (wf < w) wf *= 2;
        hf = 1;
        while (hf < h) hf *= 2;

        // Remap texture coordinates relative to (wf, hf) to select a proper
        // subrectangle for the blit - the subrectangle with actual meaningful
        // framebuffer content:
        wt = (float) w / (float) wf;
        ht = (float) h / (float) hf;
    }

	// Check for override width x height parameter in the blitterString: An integral (w,h)
	// size the blit. This allows to blit a target quad with a size different from srcfbo1, without
	// scaling or filtering it. Mostly useful in conjunction with specific shaders.
	if (strp=strstr(pString1, "OvrSize:")) {
		// Parse and assign offset:
		if (sscanf(strp, "OvrSize:%i:%i", &w, &h)!=2) {
			PsychErrorExitMsg(PsychError_internal, "In PsychBlitterIdentity(): OvrSize: blit string parameter is invalid! Parse error...\n");
		}
	}

	// Bilinear filtering of srcfbo1 texture requested?
	if (strstr(pString1, "Bilinear")) {
		// Yes. Enable it.
		bilinearfiltering = TRUE;
		glTexParameteri((*srcfbo1)->textarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri((*srcfbo1)->textarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	}
	else {
		bilinearfiltering = FALSE;
	}

	// Check for offset parameter in the blitterString: An integral (x,y)
	// offset for the destination of the blit. This allows to blit the srcfbo1, without
	// scaling or filtering it, to a different start location than (0,0):
	x=y=0;
	if (strp=strstr(pString1, "Offset:")) {
		// Parse and assign offset:
		if (sscanf(strp, "Offset:%i:%i", &x, &y)!=2) {
			PsychErrorExitMsg(PsychError_internal, "In PsychBlitterIdentity(): Offset: blit string parameter is invalid! Parse error...\n");
		}
	}

	// Check for scaling parameter:
	sx = sy = 1.0;
	if (strp=strstr(pString1, "Scaling:")) {
		// Parse and assign offset:
		if (sscanf(strp, "Scaling:%f:%f", &sx, &sy)!=2) {
			PsychErrorExitMsg(PsychError_internal, "In PsychBlitterIdentity(): Scaling: blit string parameter is invalid! Parse error...\n");
		}
	}

	if (x!=0 || y!=0 || sx!=1.0 || sy!=1.0) {
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		
		// Apply global (x,y) offset:
		glTranslatef((float) x, (float) y, 0);
		
		// Apply scaling:
		glScalef(sx, sy, 1);
	}
	
	if (PsychPrefStateGet_Verbosity()>4) printf("PTB-DEBUG: PsychBlitterIdentity: Blitting x=%i y=%i sx=%f sy=%f w=%i h=%i\n", x, y, sx, sy, w, h);
	
    if (PsychIsGLClassic(windowRecord)) {
        // OpenGL-1/2: Do the blit, using a rectangular quad:
        glBegin(GL_QUADS);
	
        // Note the swapped y-coord for textures wrt. y-coord of vertex position!
        // Texture coordinate system has origin at bottom-left, y-axis pointing upward,
        // but PTB has framebuffer coordinate system with origin at top-left, with
        // y-axis pointing downward! Normally OpenGL would have origin always bottom-left,
        // but PTB has to use a different system (changed by special gluOrtho2D) transform),
        // because our 2D coordinate system needs to conform to the standards of the old
        // Psychtoolboxes and of typical windowing systems. -- A tribute to the past.

        // Upper left vertex in window
        glTexCoord2f(0, ht);
        glVertex2f(0, 0);

        // Lower left vertex in window
        glTexCoord2f(0, 0);
        glVertex2f(0, (float) h);

        // Lower right  vertex in window
        glTexCoord2f(wt, 0);
        glVertex2f((float) w, (float) h);

        // Upper right in window
        glTexCoord2f(wt, ht);
        glVertex2f((float) w, 0);

        glEnd();
    }
    else {
        // Other. Need to emulate immediate mode GL-QUADS via GL_TRIANGLE_STRIPs:
        // Also need to avoid immediate mode functions and use our own vertex arrray
        // on our local function call stack, because this blitter can be called from
        // async flipper threads for frame-sequential stereo emulation, so can't use
        // global convenience helpers (GLBEGIN, GLEND etc.), as they'd use shared data
        // This glverts array is local on our stack, thereby thread-local. Also our
        // thread has its own OpenGL context, so no danger there either:
        GLfloat glverts[4*4] = { 0,   0,  0, h,
                                 0,  ht,  0, 0,
                                 wt,  0,  w, h,
                                 wt, ht, w,  0};

        glVertexPointer(2, GL_FLOAT, 4 * sizeof(GLfloat), &glverts[2]);
        glTexCoordPointer(2, GL_FLOAT, 4 * sizeof(GLfloat), &glverts[0]);
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    }

	if (x!=0 || y!=0 || sx!=1.0 || sy!=1.0) {
		glPopMatrix();
	}

	if (bilinearfiltering) {
		// Disable filtering again:
		glTexParameteri((*srcfbo1)->textarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri((*srcfbo1)->textarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	}

	// Done.
	return(TRUE);
}

/* PsychBlitterDisplayList()  -- Displaylist blitter.
 *
 * Blits from srcfbo1 color attachment to dstfbo by calling a premade OpenGL display list.
 * Useful for application of geometric transformations, e.g., warping during blit. Typically
 * used for geometric display undistortion.
 */
psych_bool PsychBlitterDisplayList(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc, void* hookUserData, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo)
{
	int x, y;
	GLuint gllist;
	float sx, sy;
	char* strp;
	psych_bool bilinearfiltering;

    // Not available on non-classic OpenGL. Need to find some replacement for display lists at some point in time to make this work :(
    if (!PsychIsGLClassic(windowRecord)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychBlitterDisplayList(): Called, but function unsupported on this non-OpenGL-1/2 rendering context! Aborted.\n");
        return(FALSE);
	}

	// Child protection:
	if (!(srcfbo1 && (*srcfbo1))) {
		PsychErrorExitMsg(PsychError_internal, "In PsychBlitterDisplayList(): srcfbo1 is a NULL - Pointer!!!");
	}	

	// Query display list handle:
	if (strp=strstr(hookfunc->pString1, "Handle:")) {
		// Parse and assign offset:
		if (sscanf(strp, "Handle:%i", &gllist)!=1) {
			PsychErrorExitMsg(PsychError_internal, "In PsychBlitterDisplayList(): Handle: Parse error fetching display list handle!\n");
		}
	}
	else {
		PsychErrorExitMsg(PsychError_internal, "In PsychBlitterDisplayList(): No display list handle provided or parse-error fetching display list handle!\n");
	}
	
	// Handle valid?
	if (!glIsList(gllist)) PsychErrorExitMsg(PsychError_internal, "In PsychBlitterDisplayList(): Invalid display list handle provided!\n");
	
	// Bilinear filtering of srcfbo1 texture requested?
	if (strstr(hookfunc->pString1, "Bilinear")) {
		// Yes. Enable it.
		bilinearfiltering = TRUE;
		glTexParameteri((*srcfbo1)->textarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri((*srcfbo1)->textarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	}
	else {
		bilinearfiltering = FALSE;
	}

	// Check for offset parameter in the blitterString: An integral (x,y)
	// offset for the destination of the blit. This allows to blit the srcfbo1, without
	// scaling or filtering it, to a different start location than (0,0):
	x=y=0;
	if (strp=strstr(hookfunc->pString1, "Offset:")) {
		// Parse and assign offset:
		if (sscanf(strp, "Offset:%i:%i", &x, &y)!=2) {
			PsychErrorExitMsg(PsychError_internal, "In PsychBlitterDisplayList(): Offset: blit string parameter is invalid! Parse error...\n");
		}
	}

	// Check for scaling parameter:
	sx = sy = 1.0;
	if (strp=strstr(hookfunc->pString1, "Scaling:")) {
		// Parse and assign offset:
		if (sscanf(strp, "Scaling:%f:%f", &sx, &sy)!=2) {
			PsychErrorExitMsg(PsychError_internal, "In PsychBlitterDisplayList(): Scaling: blit string parameter is invalid! Parse error...\n");
		}
	}

	if (x!=0 || y!=0 || sx!=1.0 || sy!=1.0) {
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		
		// Apply global (x,y) offset:
		glTranslatef((float) x, (float) y, 0);
		
		// Apply scaling:
		glScalef(sx, sy, 1);
	}
	
	// Note the swapped y-coord for textures wrt. y-coord of vertex position!
	// Texture coordinate system has origin at bottom-left, y-axis pointing upward,
	// but PTB has framebuffer coordinate system with origin at top-left, with
	// y-axis pointing downward! Normally OpenGL would have origin always bottom-left,
	// but PTB has to use a different system (changed by special gluOrtho2D) transform),
	// because our 2D coordinate system needs to conform to the standards of the old
	// Psychtoolboxes and of typical windowing systems. -- A tribute to the past.
	
	// Call the display list: This will perform the blit operation.
	glCallList(gllist);
	
	if (x!=0 || y!=0 || sx!=1.0 || sy!=1.0) {
		glPopMatrix();
	}

	if (bilinearfiltering) {
		// Disable filtering again:
		glTexParameteri((*srcfbo1)->textarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri((*srcfbo1)->textarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	}

	// Done.
	return(TRUE);
}

/* PsychPipelineBuiltinRenderClutBitsPlusPlus - Encode Bits++ CLUT into framebuffer.
 * 
 * This builtin routine takes the current gamma table for this windowRecord, encodes it into a Bits++
 * compatible T-Lock CLUT and renders it into the framebuffer.
 */
psych_bool PsychPipelineBuiltinRenderClutBitsPlusPlus(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc)
{
	char* strp;
	const int bitshift = 16; // Bits++ expects 16 bit numbers, but ignores 2 least significant bits --> Effective 14 bit.
	int i, x, y;
	unsigned int r, g, b;
	double t1, t2;
	y=1;
	x=0;
	
	if (windowRecord->loadGammaTableOnNextFlip != 0 || windowRecord->inRedTable == NULL) {
		if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: Bits++ CLUT encoding failed. No suitable CLUT set in Screen('LoadNormalizedGammaTable'). Skipped!\n");
		return(FALSE);
	}
	
	if (windowRecord->inTableSize < 256) {
		if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: Bits++ CLUT encoding failed. CLUT has less than the required 256 entries. Skipped!\n");
		return(FALSE);
	}
	
	if (PsychPrefStateGet_Verbosity() > 4) {  
		glFinish();
		PsychGetAdjustedPrecisionTimerSeconds(&t1);
	}

	// Options provided?
	if (strlen(hookfunc->pString1)>0) {
		// Check for override vertical position for line. Default is first scanline of display.
		if (strp=strstr(hookfunc->pString1, "yPosition=")) {
			// Parse and assign offset:
			if (sscanf(strp, "yPosition=%i", &y)!=1) {
				PsychErrorExitMsg(PsychError_user, "builtin:RenderClutBitsPlusPlus: yPosition parameter for T-Lock line position is invalid! Parse error...\n");
			}
		}

		// Check for override horizontal position for line. Default is first pixel of scanline.
		if (strp=strstr(hookfunc->pString1, "xPosition=")) {
			// Parse and assign offset:
			if (sscanf(strp, "xPosition=%i", &x)!=1) {
				PsychErrorExitMsg(PsychError_user, "builtin:RenderClutBitsPlusPlus: xPosition parameter for T-Lock line position is invalid! Parse error...\n");
			}
		}
	}
	
	// Render CLUT as sequence of single points:	
	glPointSize(1);
	glBegin(GL_POINTS);
	
	// First the T-Lock unlock key:
	glColor3ub(36, 106, 133);
	glVertex2i(x++, y);
	glColor3ub(63, 136, 163);
	glVertex2i(x++, y);
	glColor3ub(8, 19, 138);
	glVertex2i(x++, y);
	glColor3ub(211, 25, 46);
	glVertex2i(x++, y);
	glColor3ub(3, 115, 164);
	glVertex2i(x++, y);
	glColor3ub(112, 68, 9);
	glVertex2i(x++, y);
	glColor3ub(56, 41, 49);
	glVertex2i(x++, y);
	glColor3ub(34, 159, 208);
	glVertex2i(x++, y);
	glColor3ub(0, 0, 0);
	glVertex2i(x++, y);
	glColor3ub(0, 0, 0);
	glVertex2i(x++, y);
	glColor3ub(0, 0, 0);
	glVertex2i(x++, y);
	glColor3ub(0, 0, 0);
	glVertex2i(x++, y);
	
	// Now the encoded CLUT: We encode 16 bit values in a high and a low pixel,
	// Bits++ throws away the two least significant bits - get 14 bit output resolution.
	for (i=0; i<256; i++) {
		// Convert 0.0 - 1.0 float value into 0 - 2^14 -1 integer range of Bits++
		r = (unsigned int)(windowRecord->inRedTable[i] * (float)((1 << bitshift) - 1) + 0.5f);
		g = (unsigned int)(windowRecord->inGreenTable[i] * (float)((1 << bitshift) - 1) + 0.5f);
		b = (unsigned int)(windowRecord->inBlueTable[i] * (float)((1 << bitshift) - 1) + 0.5f);
		
		// Pixel with high-byte of 16 bit value:
		glColor3ub((GLubyte) ((r >> 8) & 0xff), (GLubyte) ((g >> 8) & 0xff), (GLubyte) ((b >> 8) & 0xff));
		glVertex2i(x++, y);
		
		// Pixel with low-byte of 16 bit value:
		glColor3ub((GLubyte) (r & 0xff), (GLubyte) (g  & 0xff), (GLubyte) (b  & 0xff));
		glVertex2i(x++, y);
	}
	
	glEnd();
	
	if (PsychPrefStateGet_Verbosity() > 4) {  
		glFinish();
		PsychGetAdjustedPrecisionTimerSeconds(&t2);
		printf("PTB-DEBUG: Execution time of built-in Bits++ CLUT encoder was %lf ms.\n", (t2 - t1) * 1000.0f);
	}

	// Done.
	return(TRUE);
}

/* PsychPipelineBuiltinRenderClutViaRuntime() - Encode CLUT via callback to runtime environment.
 * 
 * This builtin routine takes the current gamma table for this windowRecord and calls back into
 * the runtime environment to execute some function and passes it the CLUT.
 */
psych_bool PsychPipelineBuiltinRenderClutViaRuntime(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc)
{
	char* strp;
	char* outcmd = NULL;
	int i, cmdlen;
	double t1, t2;
	
	// Be lazy: Only execute call if there is actually a pending CLUT for update:
	if (windowRecord->inRedTable == NULL) {
		if (PsychPrefStateGet_Verbosity()>10) printf("PTB-DEBUG: PsychPipelineBuiltinRenderClutViaRuntime: No new CLUT set via Screen('LoadNormalizedGammaTable'). Skipped!\n");
		return(TRUE);
	}
	
	if (windowRecord->inTableSize < 1) {
		if (PsychPrefStateGet_Verbosity()>0) printf("PTB-ERROR: PsychPipelineBuiltinRenderClutViaRuntime: CLUT encoding failed. CLUT has less than the required 1 entries. Skipped!\n");
		return(FALSE);
	}
	
	if (PsychPrefStateGet_Verbosity() > 4) {  
		PsychGetAdjustedPrecisionTimerSeconds(&t1);
	}

	cmdlen = (int) strlen(hookfunc->pString1);
	outcmd = (char*) calloc(cmdlen + 10 + (windowRecord->inTableSize * 3 * 10), sizeof(char));
	sprintf(outcmd, "%s [", hookfunc->pString1);
	strp = &outcmd[strlen(outcmd)];
	
	for (i = 0; i < windowRecord->inTableSize; i++) {
		sprintf(strp, "%06f %06f %06f ; ", windowRecord->inRedTable[i], windowRecord->inGreenTable[i], windowRecord->inBlueTable[i]);
		strp = &outcmd[strlen(outcmd)];		
		// printf("%06f %06f %06f ; ", windowRecord->inRedTable[i], windowRecord->inGreenTable[i], windowRecord->inBlueTable[i]);
	}

	strp = &outcmd[strlen(outcmd)];
	sprintf(strp, "]); ");
	
	// Release the gamma table:
	free(windowRecord->inRedTable); windowRecord->inRedTable = NULL;
	free(windowRecord->inGreenTable); windowRecord->inGreenTable = NULL;
	free(windowRecord->inBlueTable); windowRecord->inBlueTable = NULL;
	windowRecord->inTableSize = 0;
	windowRecord->loadGammaTableOnNextFlip = 0;
	
	// Execute callback into runtime:
	PsychRuntimeEvaluateString(outcmd);
	free(outcmd);
	
	if (PsychPrefStateGet_Verbosity() > 4) {  
		PsychGetAdjustedPrecisionTimerSeconds(&t2);
		printf("PTB-DEBUG: PsychPipelineBuiltinRenderClutViaRuntime: Execution time was %lf ms.\n", (t2 - t1) * 1000.0f);
	}

	// Done.
	return(TRUE);
}

/* PsychPipelineBuiltinRenderStereoSyncLine() -- Render sync trigger lines for quad-buffered stereo contexts.
 *
 * A builtin function to be called for drawing of blue-line-sync marker lines in quad-buffered stereo mode.
 */
psych_bool PsychPipelineBuiltinRenderStereoSyncLine(PsychWindowRecordType *windowRecord, int hookId, PsychHookFunction* hookfunc)
{
	GLenum draw_buffer;
	char* strp;
	float blackpoint, r, g, b;
	float fraction = 0.25;
	float w = (float) PsychGetWidthFromRect(windowRecord->rect);
	float h = (float) PsychGetHeightFromRect(windowRecord->rect);
	r=g=b=1.0;
	
	// We default to display height minus 1 for position of sync-line, instead of the lower most row
	// of the display: This is to account for a few display drivers that are off-by-one, so they would
	// actually clip the line outside display area if provided with correct coordinates (e.g., NVidia Geforce8600M on OS/X 10.4.10 and 10.5)!
	h = h - 1;
	
	// Options provided?
	if (strlen(hookfunc->pString1)>0) {
		// Check for override vertical position for sync line. Default is last scanline of display.
		if (strp=strstr(hookfunc->pString1, "yPosition=")) {
			// Parse and assign offset:
			if (sscanf(strp, "yPosition=%f", &h)!=1) {
				PsychErrorExitMsg(PsychError_user, "builtin:RenderStereoSyncLine: yPosition parameter for horizontal stereo blue-sync line position is invalid! Parse error...\n");
			}
		}

		// Check for override horizontal fraction for sync line. Default is 25% for left eye, 75% for right eye.
		if (strp=strstr(hookfunc->pString1, "hFraction=")) {
			// Parse and assign offset:
			if ((sscanf(strp, "hFraction=%f", &fraction)!=1) || (fraction < 0.0) || (fraction > 1.0)) {
				PsychErrorExitMsg(PsychError_user, "builtin:RenderStereoSyncLine: hFraction parameter for horizontal stereo blue-sync line length is invalid!\n");
			}
		}
		
		// Check for override color of sync-line. Default is white.
		if (strp=strstr(hookfunc->pString1, "Color=")) {
			// Parse and assign offset:
			if (sscanf(strp, "Color=%f %f %f", &r, &g, &b)!=3) {
				PsychErrorExitMsg(PsychError_user, "builtin:RenderStereoSyncLine: Color spec for stereo sync-line is invalid!\n");
			}
		}
	}
	
	// Query current target buffer:
	glGetIntegerv(GL_DRAW_BUFFER, (GLint*) &draw_buffer);

	// If a FBO is bound, we use the hookId of the originating chain to find out if this is left view or right view:
	if (draw_buffer == GL_COLOR_ATTACHMENT0_EXT) draw_buffer = ((hookId == kPsychLeftFinalizerBlit) ? GL_BACK_LEFT : GL_BACK_RIGHT);

	if (draw_buffer == GL_BACK_LEFT || draw_buffer == GL_FRONT_LEFT) {
		// Left stereo buffer:
		blackpoint = fraction;
	}
	else if (draw_buffer == GL_BACK_RIGHT || draw_buffer == GL_FRONT_RIGHT) {
		// Right stereo buffer:
		blackpoint = 1 - fraction;
	}
	else {
		// No stereo buffer! No stereo mode. This routine is a no-op...
		if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Stereo sync line renderer called on non-stereo framebuffer!?!\n");  
		return(TRUE);
	}
	
	// There must not be any image content below the sync-lines! Clear out everything to black:
	glColor3f(0, 0, 0);
	glBegin(GL_QUADS);
	glVertex2f(0, h-2);
	glVertex2f(w, h-2);
	glVertex2f(w, (float) PsychGetHeightFromRect(windowRecord->rect)+1);
	glVertex2f(0, (float) PsychGetHeightFromRect(windowRecord->rect)+1);
	glEnd();

	// Draw the sync-lines:
	glLineWidth(1);
	glBegin(GL_LINES);
	glColor3f(r, g, b);
	glVertex2f(0, h);
	glVertex2f(w*blackpoint, h);
	glColor3f(0, 0, 0);
	glVertex2f(w*blackpoint, h);
	glVertex2f(w, h);
	glEnd();
	
	return(TRUE);
}

/* PsychAssignHighPrecisionTextureShaders()
 *
 * Helper function, used by all functions that create textures.
 *
 * This function checks for a given texture 'textureRecord' and its parent 'windowRecord', if
 * special GLSL texture filter shaders (for bilinear filterin) and texture lookup shaders
 * (for nearest neighbour sampling) should be assigned, based on the properties of the
 * textures, parent window and provided request flags.
 *
 * If shaders are needed, assigns them. On first invocation, creates the neccessary
 * shaders.
 *
 * Conditions under which shaders are assigned:
 * a) Calling code absolutely wants it ('userRequest' > 0).
 * b) Color clamping is disabled via Screen('ColorRange') and gfx-hw is incapable of doing
 *    this internally and at high precision, so this works around lacking hw support/precision.
 * c) It is a floating point texture format that the hw can't handle with its built-in
 *    filtering facilities, so this is a workaround for such hw. (usefloatformat == 1 for 16bpc float, 2 == 32 bpc float).
 *
 * Returns true on success, false on error, but performs error/warning output itself.
 */
psych_bool PsychAssignHighPrecisionTextureShaders(PsychWindowRecordType* textureRecord, PsychWindowRecordType* windowRecord, int usefloatformat, int userRequest)
{
	// Detect if its a GL_TEXTURE_2D texture: Currently not handled by our shaders...
	unsigned int usepoweroftwo = (PsychGetTextureTarget(textureRecord) == GL_TEXTURE_2D) ? 1 : 0;
	
	// Remap windowRecord to its parent if any. We want the associated "toplevel" onscreen window,
	// because only that contains the required settings for colorrange, shaders and gfcaps in
	// a reliable way:
	windowRecord = PsychGetParentWindow(windowRecord);
	
	// Use a GLSL shader for texture mapping / filtering?
	// We use a GLSL shader (instead of the standard nearest-neighbour, or bilinear hardware samplers),
	// if any of these is true:
	// a) It is a floating point texture and the gfx-hardware is incapable of filtering it in hardware.
	// b) Color clamping is disabled via Screen('ColorRange') and the gfx-hardware is incapable of high precision vertex color interpolation / non-clamped mode.
	// c) The script wants us to use shaders via the 'userRequest' setting > 0.	
	if ((userRequest > 0) || ((windowRecord->colorRange < 0) && !(windowRecord->gfxcaps & kPsychGfxCapVCGood)) || ((usefloatformat == 1) && !(windowRecord->gfxcaps & kPsychGfxCapFPFilter16)) || ((usefloatformat == 2) && !(windowRecord->gfxcaps & kPsychGfxCapFPFilter32))) {
		// Need filtershaders at least for bilinear filtering:
		
		// Do we have bilinear filtershader already? Don't have this stuff for GL_TEXTURE_2D btw...
		if (windowRecord->textureFilterShader == 0 && !(usepoweroftwo & 1)) {
			// Nope. Need to create one:
			windowRecord->textureFilterShader = PsychCreateGLSLProgram(textureBilinearFilterFragmentShaderSrc, textureBilinearFilterVertexShaderSrc, NULL);
			if ((windowRecord->textureFilterShader == 0) && PsychPrefStateGet_Verbosity() > 1) {
				printf("PTB-WARNING: Created a floating point texture as requested, or manual filtering wanted, but was unable to create a float filter shader.\n");
				printf("PTB-WARNING: (Custom) Filtering - and therefore anti-aliasing - of this texture won't work or at least not at the requested precision.\n");
			}
			else {
				if (PsychPrefStateGet_Verbosity() > 3) {
					if (userRequest > 0) {
						printf("PTB-INFO: GLSL fragment filtershader created for custom high quality texture filtering.\n");
					}
					else if (windowRecord->colorRange < 0) {
						printf("PTB-INFO: GLSL fragment filtershader created for high quality texture filtering in high-precision unclamped color mode\n");
					}
					else if (usefloatformat > 0) {
						printf("PTB-INFO: %i bpc Floating point texture created. This gfx-hardware doesn't support automatic filtering of such\n", 16 * usefloatformat);
						printf("PTB-INFO: textures. A GLSL fragment filtershader was generated to work-around this.\n");
					}
				}
			}
		}
		else if ((usepoweroftwo & 1) && PsychPrefStateGet_Verbosity() > 1) {
			printf("PTB-WARNING: Created a power of two texture as requested, and you wanted high-precision texturing, but don't have a GLSL filter shader.\n");
			printf("PTB-WARNING: Filtering - and therefore anti-aliasing - and high-precision drawing of power of two textures may not work.\n");
		}
		
		// Do we need a simple texture lookup & modulate with vertex color shader for non-filtered texture
		// mapping? This is not needed for floating point textures per se (all hw can nearest neighbour sample them),
		// but it is needed if user wants it (userRequest flag > 0) or if unclamped high-precision drawing is
		// requested but the hw is not up to that task:
		if ((userRequest > 0) || ((windowRecord->colorRange < 0) && !(windowRecord->gfxcaps & kPsychGfxCapVCGood))) {
			// Yes, need it - Create a shader if one doesn't yet exist:
			if (windowRecord->textureLookupShader == 0 && !(usepoweroftwo & 1)) {
				// Create one:
				windowRecord->textureLookupShader = PsychCreateGLSLProgram(textureLookupFragmentShaderSrc, textureBilinearFilterVertexShaderSrc, NULL);
				if ((windowRecord->textureLookupShader == 0) && PsychPrefStateGet_Verbosity() > 1) {
					printf("PTB-WARNING: Failed to create a texture lookup shader. High precision texture drawing therefore won't work.\n");
				}				
			}
			else if ((usepoweroftwo & 1) && PsychPrefStateGet_Verbosity() > 1) {
				printf("PTB-WARNING: Created a power of two texture as requested, and you wanted high-precision texturing, but don't have a GLSL shader for that.\n");
				printf("PTB-WARNING: High-precision drawing of power of two textures therefore may not work.\n");
			}
			
		}
		
		// Assign our onscreen windows filtershader to this texture:
		textureRecord->textureFilterShader = windowRecord->textureFilterShader;
		textureRecord->textureLookupShader = windowRecord->textureLookupShader;
	}

	// Done.
	return(TRUE);
}

psych_bool PsychAssignPlanarTextureShaders(PsychWindowRecordType* textureRecord, PsychWindowRecordType* windowRecord, int channels)
{	
	// Remap windowRecord to its parent if any. We want the associated "toplevel" onscreen window,
	// because only that contains the required shaders and gfcaps in a reliable way:
	windowRecord = PsychGetParentWindow(windowRecord);

	// Do we have a planar texture shader for this channels-count already?
	if (windowRecord->texturePlanarShader[channels - 1] == 0) {
		// Nope. Need to create one:
		switch (channels) {
			case 1:
				windowRecord->texturePlanarShader[channels - 1] = PsychCreateGLSLProgram(texturePlanar1FragmentShaderSrc, texturePlanarVertexShaderSrc, NULL);
			break;

			case 2:
				windowRecord->texturePlanarShader[channels - 1] = PsychCreateGLSLProgram(texturePlanar2FragmentShaderSrc, texturePlanarVertexShaderSrc, NULL);
			break;

			case 3:
				windowRecord->texturePlanarShader[channels - 1] = PsychCreateGLSLProgram(texturePlanar3FragmentShaderSrc, texturePlanarVertexShaderSrc, NULL);
			break;

			case 4:
				windowRecord->texturePlanarShader[channels - 1] = PsychCreateGLSLProgram(texturePlanar4FragmentShaderSrc, texturePlanarVertexShaderSrc, NULL);
			break;

			default:
				printf("PTB-BUG: In PsychAssignPlanarTextureShaders() unknown channels count of %i !!\n", channels);
				return(FALSE);
		}
		
		if (windowRecord->texturePlanarShader[channels - 1] == 0) {
			printf("PTB-ERROR: Failed to create planar filter shader for planar texture (specialFlags == 4 in Screen('MakeTexture')). This will lead to a corrupted texture!\n");
			return(FALSE);
		}
	}

	// Assign our onscreen windows planarshader to this texture:
    // We don't support GL_TEXTURE_2D textures yet though, so only auto-assign shader for rectangle textures.
	if (textureRecord && !(PsychGetTextureTarget(textureRecord) == GL_TEXTURE_2D)) textureRecord->textureFilterShader = -1 * windowRecord->texturePlanarShader[channels - 1];

	// Done.
	return(TRUE);
}

psych_bool PsychAssignPlanarI420TextureShader(PsychWindowRecordType* textureRecord, PsychWindowRecordType* windowRecord)
{
    // Remap windowRecord to its parent if any. We want the associated "toplevel" onscreen window,
    // because only that contains the required shader and gfcaps in a reliable way:
    windowRecord = PsychGetParentWindow(windowRecord);

    // Do we already have a I420 planar sampling texture shader?
    if (windowRecord->textureI420PlanarShader == 0) {
        // Nope. Need to create one:
        windowRecord->textureI420PlanarShader = PsychCreateGLSLProgram(texturePlanarI420FragmentShaderSrc, texturePlanarVertexShaderSrc, NULL);

        if (windowRecord->textureI420PlanarShader == 0) {
            printf("PTB-ERROR: Failed to create planar YUV-I420 filter shader for video texture.\n");
            return(FALSE);
        }
    }

    // Assign our onscreen windows planar I420 shader to this texture:
    // We don't support GL_TEXTURE_2D textures yet though, so only auto-assign shader for rectangle textures.    
    if (textureRecord && !(PsychGetTextureTarget(textureRecord) == GL_TEXTURE_2D)) textureRecord->textureFilterShader = -1 * windowRecord->textureI420PlanarShader;

    // Done.
    return(TRUE);
}

psych_bool PsychAssignPlanarI800TextureShader(PsychWindowRecordType* textureRecord, PsychWindowRecordType* windowRecord)
{
    // Remap windowRecord to its parent if any. We want the associated "toplevel" onscreen window,
    // because only that contains the required shader and gfcaps in a reliable way:
    windowRecord = PsychGetParentWindow(windowRecord);
    
    // Do we already have a I800 planar sampling texture shader?
    if (windowRecord->textureI800PlanarShader == 0) {
        // Nope. Need to create one:
        windowRecord->textureI800PlanarShader = PsychCreateGLSLProgram(texturePlanarI800FragmentShaderSrc, texturePlanarVertexShaderSrc, NULL);
        
        if (windowRecord->textureI800PlanarShader == 0) {
            printf("PTB-ERROR: Failed to create planar Y8-I800 filter shader for video texture.\n");
            return(FALSE);
        }
    }
    
    // Assign our onscreen windows planar I800 shader to this texture:
    // We don't support GL_TEXTURE_2D textures yet though, so only auto-assign shader for rectangle textures.    
    if (textureRecord && !(PsychGetTextureTarget(textureRecord) == GL_TEXTURE_2D)) textureRecord->textureFilterShader = -1 * windowRecord->textureI800PlanarShader;
    
    // Done.
    return(TRUE);
}
