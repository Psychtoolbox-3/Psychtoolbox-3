/*
	PsychToolbox2/Source/Common/Screen/SCREENPutImage.cpp		
  
	AUTHORS:
	
		Allen.Ingling@nyu.edu			awi 
		mario.kleiner@tuebingen.mpg.de	mk
  
	PLATFORMS:	
	
		All.
    
	HISTORY:

		02/06/03  	awi		Created.
		10/12/04	awi		In useString: moved commas to inside [].  Filled in synopsis
		1/14/05		awi		Commented out calls to enable alpha blending.
		2/25/05		awi		Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
		2/25/07		cgb		Now lives in the GLfloat world.
		2/01/08		mk		Add missing glPixelZoom(1,1); to avoid side effects.
		5/01/08		mk		Remove check for target windows number of layers - This works for all layer values as the GL handles this.
		3/20/11		mk		Make 64-bit clean.

	TO DO:
    
*/

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] =  "Screen('PutImage', windowPtr, imageArray [,rect])";
//                                             1          2            3
static char synopsisString[] =
	"Copy the matrix \"imageArray\" to a window, very slowly.\n"
    "This function is not (and will not ever get) supported on OpenGL-ES embedded graphics hardware.\n"
	"\"rect\" is in window coordinates. The whole image is copied to \"rect\", scaling if "
	"necessary. The rect default is the imageArray's rect, centered in the window. "
	"The orientation of the array in the window is identical to that of Matlab's "
	"numerical array displays in the Command Window. The first pixel is in the upper "
	"left, and the rows are horizontal. imageArray can be a pure luminance matrix, "
	"a RGB color matrix, or a RGBA color matrix with alpha channel. Its data type "
	"can be uint8 (faster) or double. Please note that this function is relatively "
	"slow and inflexible. Have a look at the 'MakeTexture' and 'DrawTexture' "
	"functions for a faster and more flexible way of drawing image matrices. ";
static char seeAlsoString[] = "GetImage OpenOffscreenWindow MakeTexture DrawTexture DrawTextures";
	
// Macro version of a function found in MiniBox.c.  Eliminates the unneeded overhead required by a function call.
// This improves speed by several milliseconds for medium to large images.
#define PSYCHINDEXELEMENTFROM3DARRAY(mDim, nDim, pDim, m, n, p) (p*mDim*nDim + n*mDim + m)

PsychError SCREENPutImage(void) 
{
	PsychRectType			windowRect, positionRect;
	int						ix, iy;
	size_t					matrixRedIndex, matrixGreenIndex, matrixBlueIndex, matrixAlphaIndex, matrixGrayIndex;
	int						inputM, inputN, inputP, positionRectWidth, positionRectHeight;
	size_t					pixelIndex = 0;
	PsychWindowRecordType	*windowRecord;
	unsigned char			*inputMatrixByte;
	double					*inputMatrixDouble;
	GLfloat					*pixelData;
	GLfloat					matrixGrayValue, matrixRedValue, matrixGreenValue, matrixBlueValue, matrixAlphaValue;
	PsychArgFormatType		inputMatrixType;
	GLfloat					xZoom = 1, yZoom = -1;
        
	// All sub functions should have these two lines.
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if (PsychIsGiveHelp()) {
		PsychGiveHelp();
		return PsychError_none;
	};

	// Cap the number of inputs.
	PsychErrorExit(PsychCapNumInputArgs(4));   //The maximum number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(0));  //The maximum number of outputs
        
	// Get the image matrix.
	inputMatrixType = PsychGetArgType(2);
	switch (inputMatrixType) {
		case PsychArgType_none:
		case PsychArgType_default:
			PsychErrorExitMsg(PsychError_user, "imageArray argument required");
			break;
		case PsychArgType_uint8:
			PsychAllocInUnsignedByteMatArg(2, TRUE, &inputM, &inputN, &inputP, &inputMatrixByte);
			break;
		case PsychArgType_double:
			PsychAllocInDoubleMatArg(2, TRUE, &inputM, &inputN, &inputP, &inputMatrixDouble);
			break;
		default:
			PsychErrorExitMsg(PsychError_user, "imageArray must be uint8 or double type");
			break;
	}

    if (inputP != 1 && inputP != 3 && inputP != 4) {
        PsychErrorExitMsg(PsychError_user, "Third dimension of image matrix must be 1, 3, or 4");
    }
        
	// Get the window and get the rect and stuff.
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);

    // A no-go on OES:
    if (PsychIsGLES(windowRecord)) {
        PsychErrorExitMsg(PsychError_unimplemented, "Sorry, Screen('PutImage') is not supported on OpenGL-ES embedded graphics hardware. Use 'MakeTexture' and 'DrawTexture' instead.");
    }

	PsychGetRectFromWindowRecord(windowRect, windowRecord);
	if (PsychCopyInRectArg(3, FALSE, positionRect)) {
		if (IsPsychRectEmpty(positionRect)) {
			return PsychError_none;
		}
		positionRectWidth  = (int) PsychGetWidthFromRect(positionRect);
		positionRectHeight = (int) PsychGetHeightFromRect(positionRect);
		if (positionRectWidth != inputN  || positionRectHeight != inputM) {
			// Calculate the zoom factor.
			xZoom = (GLfloat)   positionRectWidth  / (GLfloat) inputN;
			yZoom = -((GLfloat) positionRectHeight / (GLfloat) inputM);
		}
	}
	else {
	   positionRect[kPsychLeft] = 0;
	   positionRect[kPsychTop] = 0;
	   positionRect[kPsychRight] = inputN;
	   positionRect[kPsychBottom] = inputM;
	   PsychCenterRect(positionRect, windowRect, positionRect);
	}
        
	// Allocate memory to hold the pixel data that we'll later pass to OpenGL.
	pixelData = (GLfloat*) PsychMallocTemp(sizeof(GLfloat) * (size_t) inputN * (size_t) inputM * 4);
	
	// Loop through all rows and columns of the pixel data passed from Matlab, extract it,
	// and stick it into 'pixelData'.
	for (iy = 0; iy < inputM; iy++) {
		for (ix = 0; ix < inputN; ix++) {
			if (inputP == 1) { // Grayscale
							   // Extract the grayscale value.
				matrixGrayIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 1, (size_t) iy, (size_t) ix, 0);
				if (inputMatrixType == PsychArgType_uint8) {
					// If the color range is > 255, then force it to 255 for 8-bit values.
					matrixGrayValue = (GLfloat)inputMatrixByte[matrixGrayIndex];
					if (windowRecord->colorRange > 255) {
						matrixGrayValue /= (GLfloat)255;
					}
					else {
						matrixGrayValue /= (GLfloat)windowRecord->colorRange;
					}
				}
				else {
					matrixGrayValue = (GLfloat)(inputMatrixDouble[matrixGrayIndex] / windowRecord->colorRange);
				}
				
				// RGB will all be the same for grayscale.  We'll go ahead and fix alpha to the max value.
				pixelData[pixelIndex++] = matrixGrayValue; // R
				pixelData[pixelIndex++] = matrixGrayValue; // G
				pixelData[pixelIndex++] = matrixGrayValue; // B
				pixelData[pixelIndex++] = (GLfloat) 1.0; // A
			}
			else if (inputP == 3) { // RGB
				matrixRedIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 0);
				matrixGreenIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 1);
				matrixBlueIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 2);
				
				if (inputMatrixType == PsychArgType_uint8) {
					// If the color range is > 255, then force it to 255 for 8-bit values.
					matrixRedValue = (GLfloat)inputMatrixByte[matrixRedIndex];
					matrixGreenValue = (GLfloat)inputMatrixByte[matrixGreenIndex];
					matrixBlueValue = (GLfloat)inputMatrixByte[matrixBlueIndex];
					if (windowRecord->colorRange > 255) {
						matrixRedValue /= (GLfloat)255;
						matrixGreenValue /= (GLfloat)255;
						matrixBlueValue /= (GLfloat)255;
					}
					else {
						matrixRedValue /= (GLfloat)windowRecord->colorRange;
						matrixGreenValue /= (GLfloat)windowRecord->colorRange;
						matrixBlueValue /= (GLfloat)windowRecord->colorRange;
					}
				}
				else {
					matrixRedValue = (GLfloat)(inputMatrixDouble[matrixRedIndex] / windowRecord->colorRange);
					matrixGreenValue = (GLfloat)(inputMatrixDouble[matrixGreenIndex] / windowRecord->colorRange);
					matrixBlueValue = (GLfloat)(inputMatrixDouble[matrixBlueIndex] / windowRecord->colorRange);
				}
				
				pixelData[pixelIndex++] = matrixRedValue;
				pixelData[pixelIndex++] = matrixGreenValue;
				pixelData[pixelIndex++] = matrixBlueValue;
				pixelData[pixelIndex++] = (GLfloat)1.0;
			}
			else if (inputP == 4) { // RGBA
				matrixRedIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 0);
				matrixGreenIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 1);
				matrixBlueIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 2);
				matrixAlphaIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 3);
				
				if (inputMatrixType == PsychArgType_uint8) {
					// If the color range is > 255, then force it to 255 for 8-bit values.
					matrixRedValue = (GLfloat)inputMatrixByte[matrixRedIndex];
					matrixGreenValue = (GLfloat)inputMatrixByte[matrixGreenIndex];
					matrixBlueValue = (GLfloat)inputMatrixByte[matrixBlueIndex];
					matrixAlphaValue = (GLfloat)inputMatrixByte[matrixAlphaIndex];
					if (windowRecord->colorRange > 255) {
						matrixRedValue /= (GLfloat)255;
						matrixGreenValue /= (GLfloat)255;
						matrixBlueValue /= (GLfloat)255;
						matrixAlphaValue /= (GLfloat)255;
					}
					else {
						matrixRedValue /= (GLfloat)windowRecord->colorRange;
						matrixGreenValue /= (GLfloat)windowRecord->colorRange;
						matrixBlueValue /= (GLfloat)windowRecord->colorRange;
						matrixAlphaValue /= (GLfloat)windowRecord->colorRange;
					}
				}
				else {
					matrixRedValue = (GLfloat)(inputMatrixDouble[matrixRedIndex] / windowRecord->colorRange);
					matrixGreenValue = (GLfloat)(inputMatrixDouble[matrixGreenIndex] / (GLfloat)windowRecord->colorRange);
					matrixBlueValue = (GLfloat)(inputMatrixDouble[matrixBlueIndex] / (GLfloat)windowRecord->colorRange);
					matrixAlphaValue = (GLfloat)(inputMatrixDouble[matrixAlphaIndex] / (GLfloat)windowRecord->colorRange);
				}
				
				pixelData[pixelIndex++] = matrixRedValue;
				pixelData[pixelIndex++] = matrixGreenValue;
				pixelData[pixelIndex++] = matrixBlueValue;
				pixelData[pixelIndex++] = matrixAlphaValue;
			}
		} // for (iy = 0; iy < inputM; iy++)
	} // for (ix = 0; ix < inputN; ix++)
	
	// Enable this windowRecords framebuffer as current drawingtarget:
	PsychSetDrawingTarget(windowRecord);
	
	// Disable draw shader:
	PsychSetShader(windowRecord, 0);
	
	PsychUpdateAlphaBlendingFactorLazily(windowRecord);
	
	// Set the raster position so that we can draw starting at this location.
	glRasterPos2f((GLfloat)(positionRect[kPsychLeft]), (GLfloat)(positionRect[kPsychTop]));
	
	// Tell glDrawPixels to unpack the pixel array along GLfloat boundaries.
	glPixelStorei(GL_UNPACK_ALIGNMENT, (GLint)sizeof(GLfloat));
	
	// Dump the pixels onto the screen.
	glPixelZoom(xZoom, yZoom);
	glDrawPixels(inputN, inputM, GL_RGBA, GL_FLOAT, pixelData);
	glPixelZoom(1,1);
	
	PsychFlushGL(windowRecord);  // This does nothing if we are multi buffered, otherwise it glFlushes
	PsychTestForGLErrors();
	
	return PsychError_none;
}
