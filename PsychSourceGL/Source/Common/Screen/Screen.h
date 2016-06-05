/*
	Screen.h
	
	PLATFORMS:
			
		All.
			
	AUTHORS:
	
		Allen Ingling		awi		Allen.Ingling@nyu.edu
		Mario Kleiner		mk		mario.kleiner@tuebingen.mpg.de
		
	HISTORY:
	
		12/18/01		awi		wrote it.  Derived from joystick.h
		11/16/04		awi		Added  SCREENGlobalRect 
		1/25/05			awi		Merged in mk version which provides ScreenCloseAllWindows() declaration. Then restored SCREENGetTimeList declaration.  
		4/22/05          mk     Added new commands DrawLines, SelectStereoDrawBuffer and DrawingFinished
		5/09/05          mk     Added new command GetFlipInterval
		10/11/05         mk     Support for special movie playback added (include of PsychMovieSupport.h)
		01/13/08		 mk		Added warning text and more fixes for the totally broken #include chain issues...
 */

//begin include once 
#ifndef PSYCH_IS_INCLUDED_Screen
#define PSYCH_IS_INCLUDED_Screen

// CAUTION FIXME TODO: The whole include chain of PTB is pretty broken, containing
// lot's of circular dependencies. It was a well meant, but badly thought through
// design of the original developers, with lot's of issues creeping up with the
// evolution of PTB-3. The only clean fix would be a massive, painful redesign with
// lot's of care and lots of testing -- Nothing to be done quickly or hastily :-(
//
// For now this means that the exact order of #include statements here matters.
// A small change in order can break the whole build process in pretty nasty
// ways --> Thousands of compiler errors which are hard to track down!
// If you ever add a new #include here and then get severe and mysterious
// compile breakage, first check the order of includes here for circular
// dependency issues. Try to not change the order here unless absolutely needed
// and then only in the most minimal way that fixes your problems!!!!
//
// You may need to add pretty ugly hacks to fix'em then,
// See the comment below for PsychGraphicsHardwareHALSupport.h ...

//project includes includes for screen foundation
#include "Psych.h"
#include "PsychRects.h"
#include "ScreenTypes.h"

// Include POSIX Threading support on Unix systems (for Async Flips):
#if PSYCH_SYSTEM != PSYCH_WINDOWS
#include <pthread.h>
#endif

#ifdef PTB_USE_WAFFLE
#include <waffle.h>
#endif

// Include specifications of the GPU registers:
#include "PsychGraphicsCardRegisterSpecs.h"

#include "PsychScreenGlue.h"
#include "PsychWindowTextGlue.h"
#include "WindowBank.h"
// PsychGraphicsHardwareHALSupport *must* be behind "WindowBank" and "PsychWindowTextGlue.h"!!!
// Read the comments marked with CAUTION in WindowRecord.h and PsychGraphicsHardwareHALSupport.h
// for explanation.
#include "PsychGraphicsHardwareHALSupport.h"
#include "PsychWindowGlue.h"
#include "PsychWindowSupport.h"
#include "PsychMovieSupport.h"
#include "PsychTextureSupport.h"
#include "PsychAlphaBlending.h"
#include "PsychVideoCaptureSupport.h"
#include "PsychImagingPipelineSupport.h"
#include "PsychMovieWritingSupport.h"
#include "ScreenArguments.h"
#include "RegisterProject.h"
#include "WindowHelpers.h"

#if PSYCH_SYSTEM == PSYCH_OSX
#include "PsychFontGlue.h"
#include "FontInfo.h"
#endif

#include "ScreenPreferenceState.h"

//functions registered with the Psychtoolbox library 
PsychError ScreenExitFunction(void); 			//ScreenExit.c
PsychError PsychDisplayScreenSynopsis(void);	        //ScreenSynopsis.c

//internal screen functions
void InitializeSynopsis();				//ScreenExit.c
void ScreenCloseAllWindows();           //SCREENCloseAll.c

//PsychGLGlue.c
int		PsychConvertColorToDoubleVector(PsychColorType *color, PsychWindowRecordType *windowRecord, GLdouble *valueArray);
// int		PsychConvertColorAndColorSizeToDoubleVector(PsychColorType *color, int colorSize, GLdouble *valueArray);
void		PsychSetGLColor(PsychColorType *color, PsychWindowRecordType *windowRecord);
void		PsychSetupVertexColorArrays(PsychWindowRecordType *windowRecord, psych_bool enable, int mc, double* colors, unsigned char *bytecolors);
void		PsychSetArrayColor(PsychWindowRecordType *windowRecord, int i, int mc, double* colors, unsigned char *bytecolors);
void		PsychGLClear(PsychWindowRecordType *windowRecord);
void		PsychGLRect(PsychRectType psychRect);
char		*PsychGetGLErrorNameString(GLenum errorConstant);
#define		PsychTestForGLErrors()		PsychTestForGLErrorsC(__LINE__, __func__, __FILE__) 
void		PsychTestForGLErrorsC(int lineNum, const char *funcName, const char *fileName);
GLdouble	*PsychExtractQuadVertexFromRect(double *rect, int vertexNumber, GLdouble *vertex);
void		PsychPrepareRenderBatch(PsychWindowRecordType *windowRecord, int coords_pos, int* coords_count, double** xy, int colors_pos, int* colors_count, int* colorcomponent_count, double** colors, unsigned char** bytecolors, int sizes_pos, int* sizes_count, double** size, psych_bool usefloat);
void		PsychWaitPixelSyncToken(PsychWindowRecordType *windowRecord, psych_bool flushOnly);
psych_bool	PsychIsGLClassic(PsychWindowRecordType *windowRecord);
GLenum		PsychGLFloatType(PsychWindowRecordType *windowRecord);
#define PSYCHGLFLOAT PsychGLFloatType(windowRecord)
psych_bool	PsychIsGLES(PsychWindowRecordType *windowRecord);
#define PSYCHEXECNONGLES(x) if (!PsychIsGLES(windowRecord)) (x)

void PsychGLBegin(PsychWindowRecordType *windowRecord, GLenum primitive);
void PsychGLEnd(PsychWindowRecordType *windowRecord);
void PsychGLVertex4f(PsychWindowRecordType *windowRecord, float x, float y, float z, float w);
void PsychGLColor4f(PsychWindowRecordType *windowRecord, float r, float g, float b, float a);
void PsychGLTexCoord4f(PsychWindowRecordType *windowRecord, float s, float t, float u, float v);
void PsychGLRectd(PsychWindowRecordType *windowRecord, double x1, double y1, double x2, double y2);
void PsychDrawDisc(PsychWindowRecordType *windowRecord, float xc, float yc, float innerRadius, float outerRadius, int numSlices, float xScale, float yScale, float startAngle, float arcAngle);

#define GLBEGIN(p) PsychGLBegin(windowRecord, (p))
#define GLEND() PsychGLEnd(windowRecord)
#define GLVERTEX2f(x,y) PsychGLVertex4f(windowRecord, (x), (y), 0.0, 1.0)
#define GLVERTEX2d(x,y) PsychGLVertex4f(windowRecord, (float) (x), (float) (y), (float) 0.0, (float) 1.0)
#define GLRECTd(x1, y1, x2, y2) PsychGLRectd(windowRecord, (x1), (y1), (x2), (y2))
#define GLTEXCOORD2f(s,t) PsychGLTexCoord4f(windowRecord, (s), (t), 0.0, 1.0)

// Helper routines for vertically compressed stereo displays: Defined in SCREENSelectStereoDrawBuffer.c
int PsychSwitchCompressedStereoDrawBuffer(PsychWindowRecordType *windowRecord, int newbuffer);
void PsychComposeCompressedStereoBuffer(PsychWindowRecordType *windowRecord);

// Helper routines for text renderers:
void		PsychCleanupTextRenderer(PsychWindowRecordType* windowRecord);
psych_bool	PsychLoadTextRendererPlugin(PsychWindowRecordType* windowRecord);
void		PsychDrawCharText(PsychWindowRecordType* winRec, const char* textString, double* xp, double* yp, unsigned int yPositionIsBaseline, PsychColorType *textColor, PsychColorType *backgroundColor, PsychRectType* boundingbox);
PsychError	PsychDrawUnicodeText(PsychWindowRecordType* winRec, PsychRectType* boundingbox, unsigned int stringLengthChars, double* textUniDoubleString, double* xp, double* yp, double* theight, double* xAdvance, unsigned int yPositionIsBaseline, PsychColorType *textColor, PsychColorType *backgroundColor, int swapTextDirection);
PsychError	PsychOSDrawUnicodeText(PsychWindowRecordType* winRec, PsychRectType* boundingbox, unsigned int stringLengthChars, double* textUniDoubleString, double* xp, double* yp, unsigned int yPositionIsBaseline, PsychColorType *textColor, PsychColorType *backgroundColor);
psych_bool	PsychAllocInTextAsUnicode(int position, PsychArgRequirementType isRequired, int *textLength, double **unicodeText);
psych_bool	PsychSetUnicodeTextConversionLocale(char* mnewlocale);
const char* PsychGetUnicodeTextConversionLocale(void);

//functions implementing Screen subcommands. 
PsychError	SCREENNull(void);					
PsychError 	SCREENTestStructures(void);				
PsychError 	MODULEVersion(void);					
PsychError 	SCREENComputer(void);					
PsychError	SCREENScreens(void);					
PsychError 	SCREENPixelSize(void);					 
PsychError 	SCREENPixelSizes(void);					
PsychError 	SCREENNominalFramerate(void);			
PsychError 	SCREENOpenWindow(void);					
PsychError 	SCREENOpenOffscreenWindow(void);		
PsychError 	SCREENClose(void);						
PsychError	SCREENCloseAll(void);				
PsychError 	SCREENCopyWindow(void);				 
PsychError 	SCREENFlip(void);					
PsychError 	SCREENFillRect(void);					
PsychError	SCREENGetImage(void);					
PsychError 	SCREENPutImage(void);					
PsychError 	SCREENHideCursorHelper(void);					
PsychError 	SCREENShowCursorHelper(void);					
PsychError 	SCREENSetMouseHelper(void);				 
PsychError 	SCREENRect(void);					
PsychError 	SCREENWindowScreenNumber(void);			
PsychError 	SCREENWindows(void);					 
PsychError 	SCREENWindowKind(void);					
PsychError 	SCREENIsOffscreen(void);				
PsychError 	SCREENReadNormalizedGammaTable(void);	 
PsychError 	SCREENLoadNormalizedGammaTable(void);	
PsychError 	SCREENglPoint(void);					
PsychError 	SCREENgluDisk(void);					 
PsychError 	SCREENFillOval(void);					
PsychError 	SCREENFrameOval(void) ;					
PsychError 	SCREENTextModes(void);					
PsychError 	SCREENTextMode(void);					
PsychError 	SCREENTextSize(void);					
PsychError 	SCREENTextStyle(void);					
PsychError 	SCREENTextFont(void);					
PsychError      SCREENTextBounds(void);					
PsychError      SCREENTextTransform(void);
PsychError      SCREENDrawText(void);					
PsychError      SCREENTextColor(void);					
PsychError      SCREENPreference(void);					
PsychError      SCREENDrawTexture(void);			
PsychError      SCREENMakeTexture(void);			
PsychError      SCREENFrameRect(void);
PsychError      SCREENDrawLine(void);
PsychError      SCREENFillPoly(void);
PsychError      SCREENFramePoly(void);
PsychError      SCREENGlobalRect(void);
PsychError	SCREENDrawDots(void); 
PsychError	SCREENGetTimeList(void);
PsychError	SCREENClearTimeList(void);
PsychError	SCREENBlendFunction(void);
PsychError      SCREENWindowSize(void); 
PsychError	SCREENTextBackgroundColor(void); 
PsychError      SCREENLineStipple(void);
PsychError      SCREENSelectStereoDrawBuffer(void); 
PsychError      SCREENDrawingFinished(void); 
PsychError      SCREENDrawLines(void);
PsychError      SCREENGetFlipInterval(void);
PsychError      SCREENCloseMovie(void);
PsychError      SCREENOpenMovie(void);
PsychError      SCREENPlayMovie(void);
PsychError      SCREENSetMovieTimeIndex(void);
PsychError      SCREENGetMovieTimeIndex(void);
PsychError      SCREENGetMovieImage(void);
PsychError      SCREENglPushMatrix(void);
PsychError      SCREENglPopMatrix(void);
PsychError      SCREENglLoadIdentity(void);
PsychError      SCREENglTranslate(void);
PsychError      SCREENglScale(void);
PsychError      SCREENglRotate(void);
PsychError      SCREENPreloadTextures(void);
PsychError      SCREENFillArc(void);
PsychError      SCREENDrawArc(void);
PsychError      SCREENFrameArc(void);
PsychError      SCREENWaitBlanking(void);
PsychError      SCREENSetOpenGLTexture(void); 
PsychError      SCREENGetOpenGLTexture(void);
PsychError	SCREENVideoCaptureDevices(void);
PsychError      SCREENOpenVideoCapture(void); 
PsychError      SCREENCloseVideoCapture(void); 
PsychError      SCREENStartVideoCapture(void); 
PsychError      SCREENStopVideoCapture(void); 
PsychError      SCREENGetCapturedImage(void); 
PsychError      SCREENSetVideoCaptureParameter(void); 
PsychError      SCREENBeginOpenGL(void);
PsychError      SCREENEndOpenGL(void);
PsychError	SCREENGetOpenGLDrawMode(void);
PsychError      SCREENLoadCLUT(void);
PsychError      SCREENDisplaySize(void);
PsychError      SCREENSetOpenGLTextureFromMemPointer(void);
PsychError	SCREENColorRange(void);
PsychError	SCREENHookFunction(void);
PsychError	SCREENOpenProxy(void);
PsychError	SCREENTransformTexture(void);
PsychError	SCREENDrawTextures(void) ;
PsychError	SCREENGetWindowInfo(void);
PsychError	SCREENGetMouseHelper(void);
PsychError	SCREENResolution(void);
PsychError	SCREENResolutions(void);
PsychError	SCREENWaitUntilAsyncFlipCertain(void);
PsychError	SCREENCreateMovie(void);
PsychError	SCREENFinalizeMovie(void);
PsychError      SCREENAddAudioBufferToMovie(void);
PsychError      SCREENGetFlipInfo(void);
PsychError      SCREENConfigureDisplay(void);
PsychError      SCREENPanelFitter(void);
//PsychError SCREENSetGLSynchronous(void);		//SCREENSetGLSynchronous.c


//end include once
#endif

