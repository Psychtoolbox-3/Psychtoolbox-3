/*
        EyelinkImageTransfer.c

        AUTHORS:

        nuha@sr-research.com            nj
        li@sr-research.com              lj
        mario.kleiner.de@gmail.com      mk

        PLATFORMS:

        All.

        HISTORY:

        10/10/10    nj      Created.
        05/24/11    mk      Cleanup.
        12/20/13    lj      Fixed EyelinkImageTransfer for crash due to unallocated memory in EYEBITMAP.
        07/07/23    mk      Allow passing in uint8 L/RGB/RGBA image matrices in addition to bmp file names.
*/

#include "PsychEyelink.h"
#include "bitmap.h"

// Macro version of a function found in MiniBox.c. Eliminates the unneeded overhead required by a function call.
// This improves speed by several milliseconds for medium to large images.
#define PSYCHINDEXELEMENTFROM3DARRAY(mDim, nDim, pDim, m, n, p) (p*mDim*nDim + n*mDim + m)

// Help strings
static char useString[] = "[status] = Eyelink('ImageTransfer', imagePathOrImageMatrix [, xPosition=0][, yPosition=0][, width=0][, height=0][, trackerXPosition=0][, trackerYPosition=0][, xferoptions=0]);";

static char synopsisString[] =
"This function transfers a 8 bpc image to the tracker computer for use as a backdrop for gaze cursors.\n"
"'imagePathOrImageMatrix' can be the filename of an uncompressed 8 bpc Microsoft BMP bitmap image file "
"with 3 (RGB) or 4 (RGBA) color channels. It also can be a uint8() image matrix with either 1 layer grayscale image, "
"or a 3 layer RGB color image, or a 4 layer RGBA color image with alpha channel. You can generate the matrix "
"in Octave/Matlab, via Psychtoolbox functions like Screen('GetImage'), or you can load an image file of a supported "
"image file format via the imread() function. This allows for way more flexibility than bmp files.\n"
"'width' and 'height' define output image size for display on the tracker computer. If set to zero, "
"full input image height and width will be used. Otherwise you could use 'xPosition' and 'yPosition' "
"to define a source region in the input image that should be displayed on the tracker computer, and "
"'trackerXPosition' and 'trackerYPosition' would define the start position of that subimage on the "
"display of the tracker computer.\n"
"'xferoptions' Transfer options set with bitwise OR of the following constants, determines how the image "
"is processed before display on the tracker computers operator display:\n"
" 0 Averaging combined pixels - the default\n"
" 1 Choosing darkest and keep thin dark lines\n"
" 2 Choosing darkest and keep thin white lines and control how image size is reduced to fit tracker display\n"
" 4 Maximizes contrast for clearest image\n"
" 8 Disables the dithering of the image\n"
" 16 Converts the image to grayscale (grayscale works best for EyeLink-I, text, etc.)\n";

static char seeAlsoString[] = "";

PsychError EyelinkImageTransfer(void)
{
    int         iStatus = -1;
    int         xs = 0;
    int         ys = 0;
    int         width = 0;
    int         height = 0;
    int         xd = 0;
    int         yd = 0;
    int         i;
    int         xferoptions = 0;
    char        *filename;
    int         pitch,j,size;
    int         pitch24, pitch32;
    EYEBITMAP   bmp;
    BITMAPINFO  *BitmapInfo; /* Bitmap information */
    psych_uint8 *BitmapBits; /* Bitmap data */

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
    if (PsychIsGiveHelp()) {
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

    // assign arguments
    PsychCopyInIntegerArg(2, FALSE, &xs);
    PsychCopyInIntegerArg(3, FALSE, &ys);
    PsychCopyInIntegerArg(4, FALSE, &width);
    PsychCopyInIntegerArg(5, FALSE, &height);
    PsychCopyInIntegerArg(6, FALSE, &xd);
    PsychCopyInIntegerArg(7, FALSE, &yd);
    PsychCopyInIntegerArg(8, FALSE, &xferoptions);

    // uint8 matrix passed in with image data?
    if (PsychGetArgType(1) == PsychArgType_uint8) {
        // Yes. Convert to an image memory buffer: Code adapted from Psychtoolbox SCREENPutImage.c implementation.
        // This implementation is low performance and quite limited, but so is the whole function.
        int         ix, iy, inputM, inputN, inputP;
        size_t      matrixRedIndex, matrixGreenIndex, matrixBlueIndex, matrixAlphaIndex, matrixGrayIndex;
        size_t      pixelIndex = 0;
        psych_uint8 *inputMatrixByte;
        psych_uint8 matrixGrayValue;

        // Get psych_uint8 pixel matrix:
        PsychAllocInUnsignedByteMatArg(1, TRUE, &inputM, &inputN, &inputP, &inputMatrixByte);
        if (inputP != 1 && inputP != 3 && inputP != 4)
            PsychErrorExitMsg(PsychError_user, "Third dimension of image matrix must be 1, 3, or 4.");

        if (inputM < 1 || inputN < 1)
            PsychErrorExitMsg(PsychError_user, "Image matrix must be at least 1 x 1 pixels in size.");

        // Fake suitable BitmapInfo for later code to consume:
        BitmapInfo = malloc(sizeof(BITMAPINFO));
        BitmapInfo->bmiHeader.biWidth =inputN;
        BitmapInfo->bmiHeader.biHeight = inputM;
        BitmapInfo->bmiHeader.biBitCount = 32;

        // Allocate memory to hold the pixel data that we'll later pass to common code below:
        BitmapBits = (psych_uint8*) malloc((size_t) inputN * (size_t) inputM * 4);

        // Loop through all rows and columns of the pixel data, extract it,
        // and stick it into 'BitmapBits'.
        for (iy = inputM-1; iy >= 0; iy--) { // decrement iy for correct vertical orientation
            for (ix = 0; ix < inputN; ix++) {
                if (inputP == 1) { // Grayscale8
                    // Extract the grayscale value.
                    matrixGrayIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 1, (size_t) iy, (size_t) ix, 0);
                    matrixGrayValue = inputMatrixByte[matrixGrayIndex];

                    // RGB will all be the same for grayscale. We fix alpha to the max opaque value.
                    BitmapBits[pixelIndex++] = matrixGrayValue; // B
                    BitmapBits[pixelIndex++] = matrixGrayValue; // G
                    BitmapBits[pixelIndex++] = matrixGrayValue; // R
                    BitmapBits[pixelIndex++] = 255;             // A
                }
                else if (inputP == 3) { // RGB8
                    matrixRedIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 0);
                    matrixGreenIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 1);
                    matrixBlueIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 3, (size_t) iy, (size_t) ix, 2);

                    // Ordering is BGRA - the only ordering actually handled by the Eyelink software, contrary to what
                    // the latter bmp.format-> color channel assignment mask below implies, according to private communication
                    // with Brian Richardson from SR-Research (thanks!). Iow. the color mask is completely ignored as of SDK
                    // version 2.1, and BGRA is always expected:
                    BitmapBits[pixelIndex++] = inputMatrixByte[matrixBlueIndex];        // B
                    BitmapBits[pixelIndex++] = inputMatrixByte[matrixGreenIndex];       // G
                    BitmapBits[pixelIndex++] = inputMatrixByte[matrixRedIndex];         // R
                    BitmapBits[pixelIndex++] = 255;                                     // A
                }
                else if (inputP == 4) { // RGBA8
                    matrixRedIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 0);
                    matrixGreenIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 1);
                    matrixBlueIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 2);
                    matrixAlphaIndex = PSYCHINDEXELEMENTFROM3DARRAY((size_t) inputM, (size_t) inputN, 4, (size_t) iy, (size_t) ix, 3);

                    // Ordering is BGRA - see above for explanation:
                    BitmapBits[pixelIndex++] = inputMatrixByte[matrixBlueIndex];        // B
                    BitmapBits[pixelIndex++] = inputMatrixByte[matrixGreenIndex];       // G
                    BitmapBits[pixelIndex++] = inputMatrixByte[matrixRedIndex];         // R
                    BitmapBits[pixelIndex++] = inputMatrixByte[matrixAlphaIndex];       // A
                }
            }
        }
    }
    else {
        // No. Assume that it is the filename of a BMP image file to load and use our own
        // primitive BMP file loader. Note this loader provices data in whatever ordering
        // the BMP image file contains, despite the fact that el_bitmap_to_backdrop() can
        // only handle BGR[A] ordering. According to Wikipedia, BGR or BGRA ordering is
        // the most common ordering in MS BMP image files, so this will work most of the
        // time, but there may be exotic BMP files which will result in swapped red and
        // blue color channels on the tracker host computers backdrop image:
        PsychAllocInCharArg(1, TRUE, &filename);
        BitmapBits = LoadDIBitmap(filename, &BitmapInfo);
    }

    if (BitmapBits == NULL) {
        PsychErrorExitMsg(PsychError_user, "Loading Bitmap failed\n");
    }

    bmp.w = (int) BitmapInfo->bmiHeader.biWidth;
    bmp.h = (int) BitmapInfo->bmiHeader.biHeight;
    bmp.depth = (int) 32;
    bmp.pitch = (int) bmp.w * 4;

    pPixels = BitmapBits;

    if (FALSE)
        printf("bmp.w %i bmp.h %i bibits %i :: ", bmp.w, bmp.h, (int) BitmapInfo->bmiHeader.biBitCount);

    // invert
    pitch24 = ((abs(bmp.w) * 24 + 31) / 32) * 4;
    pitch32 = abs(bmp.w) * 4;
    pDest   = (byte*) PsychMallocTemp(sizeof(byte) * pitch32 * abs(bmp.h));
    pReturn = (byte*) PsychMallocTemp(sizeof(byte) * pitch32 * abs(bmp.h));
    if (!pDest || !pReturn)
        PsychErrorExitMsg(PsychError_outofMemory, "Out of memory while trying to load bitmap!");

    memset(pDest,   0, pitch32 * abs(bmp.h));
    memset(pReturn, 0, pitch32 * abs(bmp.h));

    if ((int) BitmapInfo->bmiHeader.biBitCount == 24) {
        for (i = 0; i < abs(bmp.h); ++i)
            for (j = 0; j < abs(bmp.w); ++j) {
                pDest[i*pitch32+j*4+0] = BitmapBits[i*pitch24+j*3+0];
                pDest[i*pitch32+j*4+1] = BitmapBits[i*pitch24+j*3+1];
                pDest[i*pitch32+j*4+2] = BitmapBits[i*pitch24+j*3+2];
                pDest[i*pitch32+j*4+3] = 255;
            }

            if (bmp.w > 0 && bmp.h > 0) {
                pitch = bmp.w*4; // for 32bits img

                for (i = 0; i < bmp.h; ++i) {
                    dest = pReturn +((bmp.h-1-i)*(pitch));
                    src = pDest +(i*pitch);
                    for (j = 0; j < pitch; ++j)
                        dest[j]= src[j];
                }
            }
            else {
                memcpy(pReturn, pDest, abs(bmp.w)*4*abs(bmp.h));
            }
            bmp.pixels = pReturn;
    } else {
        if (bmp.w > 0 && bmp.h > 0) {
            pitch = bmp.w*4;
            for (i = 0; i < bmp.h; ++i)
            {
                dest = pReturn +((bmp.h-1-i)*(pitch));
                src = pPixels +(i*pitch);
                for (j = 0; j < pitch; ++j)
                    dest[j] = src[j];
            }
        } else {
            memcpy(pReturn, pPixels, abs(bmp.w)*4*abs(bmp.h));
        }

        bmp.pixels = pReturn;
    }

    bmp.format = (EYEPIXELFORMAT* ) PsychMallocTemp(sizeof(EYEPIXELFORMAT));
    memset(bmp.format, 0, sizeof(EYEPIXELFORMAT));

    // FIXME? Remove this mask assignment? According to SR-Research, as of v2.1,
    // this color mask is completely ignored, always assuming BGRA ordering, but
    // if this mask here would be used, it would be wrong, causing swapped red
    // and blue color channels.
    bmp.format->Rmask = 0x000000ff;
    bmp.format->Gmask = 0x0000ff00;
    bmp.format->Bmask = 0x00ff0000;
    bmp.format->Amask = 0xff000000;

    iStatus = el_bitmap_to_backdrop(&bmp, xs, ys, width, height, xd, yd, xferoptions);

    if (BitmapInfo) {
        free(BitmapInfo);
        if (BitmapBits)
            free(BitmapBits);
    }

    if(bmp.format)
        PsychFreeTemp(bmp.format);

    PsychCopyOutDoubleArg(1, FALSE, iStatus);

    if (FALSE) {
        for (i=0; i < 12; i++)
            printf("%i ", ((unsigned char*) bmp.pixels)[i]);
        printf("\n");
    }

    return(PsychError_none);
}
