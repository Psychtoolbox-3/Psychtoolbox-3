/*
	PsychToolbox3/Source/Common/Screen/PsychImagingPipelineSupport.h
	
	PLATFORMS:	
	
		All.  
				
	AUTHORS:
	
		Mario Kleiner           mk              mario.kleiner at tuebingen.mpg.de

	HISTORY:
	
		12/05/06	mk	Wrote it.
		
	DESCRIPTION:
	
		Infrastructure for all Screen imaging pipeline functions, i.e., hook callback functions and chains
		and GLSL based internal image processing pipeline.
		
	NOTES:
	
	TO DO: 
	
*/

//include once
#ifndef PSYCH_IS_INCLUDED_PsychImagingPipelineSupport
#define PSYCH_IS_INCLUDED_PsychImagingPipelineSupport

#include "Screen.h"

// Definition of a pointer to a blitter function: See below for conforming blitter function prototypes:
typedef psych_bool (*PsychBlitterFunc)(PsychWindowRecordType*, PsychHookFunction*, void*, psych_bool, psych_bool, PsychFBO**, PsychFBO**, PsychFBO**, PsychFBO**);

// Symbolic names for Hook-Chains: The text strings for user-space and synopsis strings are in the PsychImagingPipelineSupport.c file.
typedef enum {
	kPsychCloseWindowPreGLShutdown =				0,	// Called when closing an onscreen window while context is still available.
	kPsychCloseWindowPostGLShutdown =				1,	// Ditto., but after OpenGL context shutdown -- only non-GL ops possible here.
	kPsychUserspaceBufferDrawingFinished =			2,	// Called as part of 'DrawingFinished' or 'Flip' to do generic after drawing ops.
	kPsychStereoLeftCompositingBlit =				3,	// Image processing on the mono channel or left channel (in stereomode).
	kPsychStereoRightCompositingBlit =				4,	// Ditto for right channel in stereomode. Can add generic image processing here.
	kPsychStereoCompositingBlit =					5,	// Performs merge operation from two stereo channels to single output in all modes but quad-buffered stereo.
	kPsychPostCompositingBlit =						6,	// Not used for now.
	kPsychFinalOutputFormattingBlit =				7,	// Performs unified final data conversion/image processing on both output views immediately before hitting framebuffer.
	kPsychUserspaceBufferDrawingPrepare =			8,	// Prepare transition to userspace after 'Flip' command.
	kPsychIdentityBlit=								9,	// Standard blit chain, used when nothing else available: Copy images from one buffer to next.
	kPsychLeftFinalizerBlit=						10,	// Very last (single-pass only!) operations on left- or mono channel, e.g., drawing stereo sync lines.
	kPsychRightFinalizerBlit=						11, // Same for right channel in stereo modes.
	kPsychUserDefinedBlit=							12, // User defined image processing, e.g., for Screen('TransformTexture').
	kPsychFinalOutputFormattingBlit0 =				13,	// Performs final data conversion/image processing on output view 0 immediately before hitting framebuffer.
	kPsychFinalOutputFormattingBlit1 =				14,	// Performs final data conversion/image processing on output view 1 immediately before hitting framebuffer.
	kPsychScreenFlipImpliedOperations =				15, // Performs synchronous operations on the masterthread at time of execution of Screen('Flip'), Screen('FlipCheckEnd') and Screen('AsyncFlipEnd') after a flip is truly finished.
	kPsychPreSwapbuffersOperations =				16, // Called before emitting the PsychOSFlipWindowBuffers() call, ie., less than 1 video refresh from flipdeadline away. No OpenGL ops allowed!
} PsychHookType;

// API for PTB core:
void	PsychInitImagingPipelineDefaultsForWindowRecord(PsychWindowRecordType *windowRecord);
void	PsychInitializeImagingPipeline(PsychWindowRecordType *windowRecord, int imagingmode, int multiSample);
void	PsychShutdownImagingPipeline(PsychWindowRecordType *windowRecord, psych_bool openglpart);
void	PsychPipelineListAllHooks(PsychWindowRecordType *windowRecord);
void	PsychPipelineDumpAllHooks(PsychWindowRecordType *windowRecord);
void	PsychPipelineDumpHook(PsychWindowRecordType *windowRecord, const char* hookString);
void	PsychPipelineDisableHook(PsychWindowRecordType *windowRecord, const char* hookString);
void	PsychPipelineEnableHook(PsychWindowRecordType *windowRecord, const char* hookString);
void	PsychPipelineResetHook(PsychWindowRecordType *windowRecord, const char* hookString);
int		PsychPipelineQueryHookSlot(PsychWindowRecordType *windowRecord, const char* hookString, char** insertString, char** idString, char** blitterString, double* doubleptr, double* shaderid, double* luttexid1);
void	PsychPipelineDeleteHookSlot(PsychWindowRecordType *windowRecord, const char* hookString, int slotid);
void	PsychPipelineAddBuiltinFunctionToHook(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, const char* configString);
void	PsychPipelineAddRuntimeFunctionToHook(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, const char* evalString);
void	PsychPipelineAddCFunctionToHook(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, void* procPtr);
void	PsychPipelineAddShaderToHook(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, unsigned int shaderid, const char* blitterString, unsigned int luttexid1);

psych_bool	PsychPipelineExecuteHook(PsychWindowRecordType *windowRecord, int hookId, void* hookUserData, void* hookBlitterFunction, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo);
psych_bool	PsychPipelineExecuteHookSlot(PsychWindowRecordType *windowRecord, int hookId, PsychHookFunction* hookfunc, void* hookUserData, void* hookBlitterFunction, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo);
int		PsychPipelineProcessMacros(PsychWindowRecordType *windowRecord, char* cmdString);

// Internal helper functions:

PsychHookFunction* PsychAddNewHookFunction(PsychWindowRecordType *windowRecord, const char* hookString, const char* idString, int where, int hookfunctype);
int		PsychGetHookByName(const char* hookName);

// Setup source -> rendertarget binding for next rendering pass:
void	PsychPipelineSetupRenderFlow(PsychFBO* srcfbo1, PsychFBO* srcfbo2, PsychFBO* dstfbo, psych_bool scissor_ignore);

// Create OpenGL framebuffer object for internal rendering, setup PTB info struct for it:
psych_bool PsychCreateFBO(PsychFBO** fbo, GLenum fboInternalFormat, psych_bool needzbuffer, int width, int height, int multisample);

// Check if provided PTB texture already has a PsychFBO attached. Do nothing if so. If a FBO is missing, create one:
void PsychCreateShadowFBOForTexture(PsychWindowRecordType *textureRecord, psych_bool asRendertarget, int forImagingmode);

// On demand normalization of texture orientation to standard upright format:
void PsychNormalizeTextureOrientation(PsychWindowRecordType *sourceRecord);

psych_bool PsychIsHookChainOperational(PsychWindowRecordType *windowRecord, int hookid);
psych_bool PsychPipelineExecuteBlitter(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc, void* hookUserData, void* hookBlitterFunction, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo);

// Try to create GLSL shader from source strings and return handle to new shader.
GLuint  PsychCreateGLSLProgram(const char* fragmentsrc, const char* vertexsrc, const char* primitivesrc);

// Assign special filter/lookup shaders to textures, e.g., in HDR mode, for float textures, etc...
psych_bool PsychAssignHighPrecisionTextureShaders(PsychWindowRecordType* textureRecord, PsychWindowRecordType* windowRecord, int usefloatformat, int userRequest);
psych_bool PsychAssignPlanarTextureShaders(PsychWindowRecordType* textureRecord, PsychWindowRecordType* windowRecord, int channels);
psych_bool PsychAssignPlanarI420TextureShader(PsychWindowRecordType* textureRecord, PsychWindowRecordType* windowRecord);
psych_bool PsychAssignPlanarI800TextureShader(PsychWindowRecordType* textureRecord, PsychWindowRecordType* windowRecord);

// Builtin functions:

// PsychPipelineBuiltinRenderClutBitsPlusPlus - Encode Bits++ CLUT into framebuffer.
psych_bool PsychPipelineBuiltinRenderClutBitsPlusPlus(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc);
// PsychPipelineBuiltinRenderStereoSyncLine - Rendering of blue-sync lines for stereo shutter glasses in quad-buffered stereo mode:
psych_bool PsychPipelineBuiltinRenderStereoSyncLine(PsychWindowRecordType *windowRecord, int hookId, PsychHookFunction* hookfunc);
// PsychPipelineBuiltinRenderClutViaRuntime - Encode CLUT via callback to runtime environment.
psych_bool PsychPipelineBuiltinRenderClutViaRuntime(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc);

// Blitter functions: Assignable to a function pointer of type PsychBlitterFunc:
// =============================================================================

// Identity blitter: Blits from srcfbo1 color attachment to dstfbo without geometric transformations or other extras.
// This is the most common one for one-to-one copies or simple shader image processing. It gets automatically used
// when no special (non-default) blitter is requested by core code or users blitter parameter string:
psych_bool PsychBlitterIdentity(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc, void* hookUserData, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo);
psych_bool PsychBlitterDisplayList(PsychWindowRecordType *windowRecord, PsychHookFunction* hookfunc, void* hookUserData, psych_bool srcIsReadonly, psych_bool allowFBOSwizzle, PsychFBO** srcfbo1, PsychFBO** srcfbo2, PsychFBO** dstfbo, PsychFBO** bouncefbo);

//end include once
#endif
