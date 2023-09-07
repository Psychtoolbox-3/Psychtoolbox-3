/*
        EyelinkImageTransfer.c

        AUTHORS:

        nuha@sr-research.com            nj
        li@sr-research.com              lj
        mario.kleiner.de@gmail.com      mk
	brian@sr-research.com		br

        PLATFORMS:

        All.

        HISTORY:

        10/10/10    nj      Created.
        05/24/11    mk      Cleanup.
        12/20/13    lj      Fixed EyelinkImageTransfer for crash due to unallocated memory in EYEBITMAP.
        07/07/23    mk      Allow passing in uint8 L/RGB/RGBA image matrices in addition to bmp file names.
	07/11/23    br      Alternate implementation for accepting image matricies, tested /w updated SR Research demos
*/

#include "PsychEyelink.h"
#include "bitmap.h"

// Macro version of a function found in MiniBox.c. Eliminates the unneeded overhead required by a function call.
// This improves speed by several milliseconds for medium to large images.
// Adapted from SCREENPutImage.c
#define PSYCHINDEXELEMENTFROM3DARRAY(mDim, nDim, pDim, m, n, p) (p*mDim*nDim + n*mDim + m)


// Help strings
static char useString[] =
"[status] = Eyelink('ImageTransfer', imagePath [, xPosition=0][, yPosition=0][, width=0][, height=0][, trackerXPosition=0][, trackerYPosition=0][, xferoptions=0]);\n"
"[status] = Eyelink('ImageTransfer', imageArray [, xPosition=0][, yPosition=0][, width=0][, height=0][, trackerXPosition=0][, trackerYPosition=0][, xferoptions=0]);";

static char synopsisString[] =
"This function transfers an image to the tracker PC as backdrop for gaze cursors.\n"
"ImageTransfer has two supported methods of operation:\n"
"\tEyelink('ImageTransfer', imagePath, ...)\n"
"\t\tWhere imagePath is a string containing path/filename to a .bmp image\n"
"\t\tthat must be a 24-bit or 32-bit bitmap (i.e. only .bmp format images).\n"
"\tEyelink('ImageTransfer', imageArray, ...)\n"
"\t\tWhere imageArray is image data passed in from the MATLAB workspace. e.g. from:\n"
"\t\t\timageArray = imread()\n"
"\t\tand supports image data that:\n"
"\t\t\t- is from any image format supported by imread()\n"
"\t\t\t- has a greyscale, RGB, or RGBA color space\n"
"\t\t\t- is uint8, single (16-bit/32-bit float), or double type representation\n"
"'width' and 'height' define image size as it will be displayed on the host overlay view.\n"
"\tIf width/height are set to zero, the original image height and width will be used.\n"
"'xferoptions' Transfer options set with bitwise OR of the following constants, determines how bitmap is processed:\n"
"\t0 Averaging combined pixels\n"
"\t1 Choosing darkest and keep thin dark lines\n"
"\t2 Choosing darkest and keep thin white lines and control how bitmap size is reduced to fit tracker display\n"
"\t4 Maximizes contrast for clearest image\n"
"\t8 Disables the dithering of the image\n"
"\t16 Converts the image to grayscale (works best for EyeLinkI,text,etc.)\n";

static char seeAlsoString[] = "";

PsychError EyelinkImageTransfer(void)
{
	int				iStatus = -1;
	int				xs = 0;
	int				ys = 0;
	int				width = 0;
	int				height = 0;
	int				xd = 0;
	int				yd = 0;
	int				i;
	int				xferoptions = 0;
	char			*filename;		
	int				pitch,j,size;
	int				pitch24, pitch32;
	EYEBITMAP		bmp;		
	BITMAPINFO		*BitmapInfo; /* Bitmap information */
	psych_uint8		*BitmapBits; /* Bitmap data */
	psych_bool		bitmapFromFile = TRUE;
	PsychArgFormatType 	inputMatrixType;
	psych_uint8			*inputMatrixByte;
    double				*inputMatrixDouble;
	float				*inputMatrixSingle;
	int					inputM, inputN, inputP, ix, iy;
	size_t 				matrixGrayIndex, matrixRedIndex, matrixGreenIndex, matrixBlueIndex, matrixAlphaIndex;
	psych_uint8 		matrixGrayValue, matrixRedValue, matrixGreenValue, matrixBlueValue, matrixAlphaValue;
	size_t				pixelIndex = 0;
	psych_uint8			*pixelData;
	
#if PSYCH_SYSTEM == PSYCH_WINDOWS
    BYTE* pPixels;
    BYTE* dest;
    BYTE* src;
    BYTE* pDest;
    BYTE* pReturn;
#else
    byte* pPixels;
    byte* dest;
    byte* src;
    byte* pDest;
    byte* pReturn;
#endif
    
	PsychPushHelp(useString, synopsisString, seeAlsoString);
    
	// Output help if asked
	if(PsychIsGiveHelp()) {
		PsychGiveHelp();
		return(PsychError_none);
	}
    
	// Check arguments
	PsychErrorExit(PsychCapNumInputArgs(8));
	PsychErrorExit(PsychRequireNumInputArgs(1));
	PsychErrorExit(PsychCapNumOutputArgs(1));
    
	// Verify eyelink is up and running
	EyelinkSystemIsConnected();
	EyelinkSystemIsInitialized();


	inputMatrixType = PsychGetArgType(1);
	switch (inputMatrixType) {
        case PsychArgType_none:
        case PsychArgType_default:
            PsychErrorExitMsg(PsychError_user, "image filename or imageArray required as first argument");
            break;
		case PsychArgType_char:
			PsychAllocInCharArg(1, kPsychArgRequired, &filename);
			bitmapFromFile = TRUE;
			break;
        case PsychArgType_uint8:
            PsychAllocInUnsignedByteMatArg(1, kPsychArgRequired, &inputM, &inputN, &inputP, &inputMatrixByte);
			bitmapFromFile = FALSE;
            break;
        case PsychArgType_double:
            PsychAllocInDoubleMatArg(1, kPsychArgRequired, &inputM, &inputN, &inputP, &inputMatrixDouble);
			bitmapFromFile = FALSE;
            break;
		case PsychArgType_single:
			PsychAllocInFloatMatArg(1, kPsychArgRequired, &inputM, &inputN, &inputP, &inputMatrixSingle);
			bitmapFromFile = FALSE;
			break;
        default:
			mexPrintf("first argument has type: %d\n", inputMatrixType);
            PsychErrorExitMsg(PsychError_user, "first argument must be char type if filename, or uint8/single/double type if imageArray");
            break;
	}

	PsychCopyInIntegerArg(2, FALSE, &xs);
	PsychCopyInIntegerArg(3, FALSE, &ys);
	PsychCopyInIntegerArg(4, FALSE, &width);
	PsychCopyInIntegerArg(5, FALSE, &height);
	PsychCopyInIntegerArg(6, FALSE, &xd);
	PsychCopyInIntegerArg(7, FALSE, &yd);
	PsychCopyInIntegerArg(8, FALSE, &xferoptions);
    
	if (bitmapFromFile) {
		BitmapBits = LoadDIBitmap(filename, &BitmapInfo);
		
		if (BitmapBits == NULL){		
			PsychErrorExitMsg(PsychError_user, "Loading Bitmap failed\n");
		}
		
		bmp.w = (int) BitmapInfo->bmiHeader.biWidth;
		bmp.h = (int) BitmapInfo->bmiHeader.biHeight;
		bmp.depth = (int) 32;
		bmp.pitch = (int) bmp.w * 4;
		
		pPixels = BitmapBits;
		
		size = (((int) BitmapInfo->bmiHeader.biBitCount * (int) BitmapInfo->bmiHeader.biWidth + 31) / 32) * 4 * bmp.h;  
		
		// invert
		pitch24 = ((abs(bmp.w) * 24 + 31) / 32) * 4; 
		pitch32 = abs(bmp.w) * 4;
		pDest   = (byte*) PsychMallocTemp(sizeof(byte) * pitch32 * abs(bmp.h));
		pReturn = (byte*) PsychMallocTemp(sizeof(byte) * pitch32 * abs(bmp.h));
		if (!pDest || !pReturn) PsychErrorExitMsg(PsychError_outofMemory, "Out of memory while trying to load bitmap!");
		
		memset(pDest,   0, pitch32 * abs(bmp.h));
		memset(pReturn, 0, pitch32 * abs(bmp.h));
		
		if ((int) BitmapInfo->bmiHeader.biBitCount == 24) {
			for( i=0; i<abs(bmp.h); ++i)
				for( j=0; j<abs(bmp.w); ++j) {
					pDest[i*pitch32+j*4+0] = BitmapBits[i*pitch24+j*3+0]; 
					pDest[i*pitch32+j*4+1] = BitmapBits[i*pitch24+j*3+1];
					pDest[i*pitch32+j*4+2] = BitmapBits[i*pitch24+j*3+2];
					pDest[i*pitch32+j*4+3] = 0;
				}
					
					if(bmp.w>0 && bmp.h>0) {
						pitch = bmp.w*4; // for 32bits img
						
						for(i=0; i<bmp.h; ++i) {
							dest = pReturn +((bmp.h-1-i)*(pitch));
							src = pDest +(i*pitch);
							for(j=0; j<pitch; ++j)
								dest[j]= src[j];
						}
					}
					else {
						memcpy(pReturn, pDest, abs(bmp.w)*4*abs(bmp.h));
					}
					bmp.pixels = pReturn;
		} else {
			if(bmp.w>0 && bmp.h>0) {
				pitch = bmp.w*4;	
				for(i=0; i<bmp.h; ++i)
				{
					dest = pReturn +((bmp.h-1-i)*(pitch));
					src = pPixels +(i*pitch);
					for(j=0; j<pitch; ++j)
						dest[j]= src[j];
				}
			} else {
				memcpy(pReturn, pPixels, abs(bmp.w)*4*abs(bmp.h));
			}
			
			bmp.pixels = pReturn; 
		}
		bmp.format = (EYEPIXELFORMAT* ) PsychMallocTemp(sizeof(EYEPIXELFORMAT));
		memset(bmp.format, 0, sizeof(EYEPIXELFORMAT));
		
		
		bmp.format->Rmask = 0x000000ff;
		bmp.format->Gmask = 0x0000ff00;
		bmp.format->Bmask = 0x00ff0000;
		bmp.format->Amask = 0xff000000;
	
		iStatus = el_bitmap_to_backdrop(&bmp, xs, ys, width, height, xd, yd, xferoptions);
		
		if (BitmapInfo) {
			free(BitmapInfo);
			if (BitmapBits) free(BitmapBits);
		}
		if(bmp.format)
		{		
			if(bmp.format) PsychFreeTemp(bmp.format);
		}

	} else {

		
		// Allocate memory to hold the pixel data that we'll later pass to OpenGL.
    	// pixelData = (GLfloat*) PsychMallocTemp(sizeof(GLfloat) * (size_t) inputN * (size_t) inputM * 4);
		pixelData = (psych_uint8*) PsychMallocTemp(sizeof(psych_uint8) * (size_t) inputN * (size_t) inputM * 4);

		// Loop through all rows and columns of the pixel data passed from Matlab, extract it,
		// and stick it into 'pixelData'.
		for (iy = 0; iy < inputM; iy++) {
			for (ix = 0; ix < inputN; ix++) {
				if (inputP == 1) { // Grayscale
					// Extract the grayscale value.
					matrixGrayIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 1, (size_t) iy, (size_t) ix, 0);
					if (inputMatrixType == PsychArgType_uint8) {
						matrixGrayValue = inputMatrixByte[matrixGrayIndex];
					}
					else if (inputMatrixType == PsychArgType_double) {
						matrixGrayValue = (psych_uint8)((inputMatrixDouble[matrixGrayIndex] * 255) + 0.5);
					}
					else if (inputMatrixType == PsychArgType_single) {
						matrixGrayValue = (psych_uint8)((inputMatrixSingle[matrixGrayIndex] * 255) + 0.5);
					}

					// RGB will all be the same for grayscale.  We'll go ahead and fix alpha to the max value.
					pixelData[pixelIndex++] = matrixGrayValue; // R
					pixelData[pixelIndex++] = matrixGrayValue; // G
					pixelData[pixelIndex++] = matrixGrayValue; // B
					pixelData[pixelIndex++] = (psych_uint8) 255; // A
				}
				else if (inputP == 3) { // RGB
					matrixRedIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 0);
					matrixGreenIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 1);
					matrixBlueIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 2);

					if (inputMatrixType == PsychArgType_uint8) {
						matrixRedValue = inputMatrixByte[matrixRedIndex];
						matrixGreenValue = inputMatrixByte[matrixGreenIndex];
						matrixBlueValue = inputMatrixByte[matrixBlueIndex];
					}
					else if (inputMatrixType == PsychArgType_double) {
						matrixRedValue = (psych_uint8)((inputMatrixDouble[matrixRedIndex] * 255) + 0.5);
						matrixGreenValue = (psych_uint8)((inputMatrixDouble[matrixGreenIndex] * 255) + 0.5);
						matrixBlueValue = (psych_uint8)((inputMatrixDouble[matrixBlueIndex] * 255) + 0.5);
					}
					else if (inputMatrixType == PsychArgType_single) {
						matrixRedValue = (psych_uint8)((inputMatrixSingle[matrixRedIndex] * 255) + 0.5);
						matrixGreenValue = (psych_uint8)((inputMatrixSingle[matrixGreenIndex] * 255) + 0.5);
						matrixBlueValue = (psych_uint8)((inputMatrixSingle[matrixBlueIndex] * 255) + 0.5);						
					}

					pixelData[pixelIndex++] = matrixBlueValue;
					pixelData[pixelIndex++] = matrixGreenValue;
					pixelData[pixelIndex++] = matrixRedValue;
					pixelData[pixelIndex++] = (psych_uint8)255;
				}
				else if (inputP == 4) { // RGBA
					matrixRedIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 0);
					matrixGreenIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 1);
					matrixBlueIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 2);
					matrixAlphaIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 3);

					if (inputMatrixType == PsychArgType_uint8) {
						matrixRedValue = inputMatrixByte[matrixRedIndex];
						matrixGreenValue = inputMatrixByte[matrixGreenIndex];
						matrixBlueValue = inputMatrixByte[matrixBlueIndex];
						matrixAlphaValue = inputMatrixByte[matrixAlphaIndex];
					}
					else if (inputMatrixType == PsychArgType_double) {
						matrixRedValue = (psych_uint8)((inputMatrixDouble[matrixRedIndex] * 255) + 0.5);
						matrixGreenValue = (psych_uint8)((inputMatrixDouble[matrixGreenIndex] * 255) + 0.5);
						matrixBlueValue = (psych_uint8)((inputMatrixDouble[matrixBlueIndex] * 255) + 0.5);
						matrixAlphaValue = (psych_uint8)((inputMatrixDouble[matrixAlphaIndex] * 255) + 0.5);
					}
					else if (inputMatrixType == PsychArgType_single) {
						matrixRedValue = (psych_uint8)((inputMatrixSingle[matrixRedIndex] * 255) + 0.5);
						matrixGreenValue = (psych_uint8)((inputMatrixSingle[matrixGreenIndex] * 255) + 0.5);
						matrixBlueValue = (psych_uint8)((inputMatrixSingle[matrixBlueIndex] * 255) + 0.5);
						matrixAlphaValue = (psych_uint8)((inputMatrixSingle[matrixAlphaIndex] * 255) + 0.5);					
					}

					pixelData[pixelIndex++] = matrixBlueValue;
					pixelData[pixelIndex++] = matrixGreenValue;
					pixelData[pixelIndex++] = matrixRedValue;
					pixelData[pixelIndex++] = matrixAlphaValue;
				}
			} // for (iy = 0; iy < inputM; iy++)
		} // for (ix = 0; ix < inputN; ix++)
		
		bmp.pixels = pixelData;

		bmp.w = inputN;
		bmp.h = inputM;
		bmp.depth = (int) 32;
		bmp.pitch = (int) bmp.w * 4;

		bmp.format = (EYEPIXELFORMAT* ) PsychMallocTemp(sizeof(EYEPIXELFORMAT));
		memset(bmp.format, 0, sizeof(EYEPIXELFORMAT));
		bmp.format->Rmask = 0x000000ff;
		bmp.format->Gmask = 0x0000ff00;
		bmp.format->Bmask = 0x00ff0000;
		bmp.format->Amask = 0xff000000;
	
		iStatus = el_bitmap_to_backdrop(&bmp, xs, ys, width, height, xd, yd, xferoptions);
		
		if(bmp.format)
		{		
			if(bmp.format) PsychFreeTemp(bmp.format);
		}
	}
    
	PsychCopyOutDoubleArg(1, FALSE, iStatus);
    
	return(PsychError_none);
}
