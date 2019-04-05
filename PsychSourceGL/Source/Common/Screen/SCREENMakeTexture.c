/*
 *        SCREENMakeTexture.c
 *
 *        AUTHORS:
 *
 *                Allen Ingling        awi      Allen.Ingling@nyu.edu
 *                Mario Kleiner        mk       mario.kleiner.de@gmail.com
 *
 *        PLATFORMS:
 *
 *                All.
 *
 *        HISTORY:
 *
 *                6/25/04       awi     Created.
 *                1/4/05        mk      Performance optimizations, some bug-fixes.
 *                1/13/05       awi     Merged in Marios's changes from 1/4/05 into the master at Psychtoolbox.org.
 *                1/19/05       awi     Removed unused variables to eliminate compiler warnings.
 *                1/26/05       awi     Added StoreNowTime() calls.
 *                3/19/11       mk      Make 64-bit clean.
 *
 *        DESCRIPTION:
 *
 *                Conversion code for creating OpenGL textures from Matlab/Octave image matrices. Handles uint8 and
 *                double matrices as input. Converts into 8bpc textures by default, but also supports half_float and
 *                float 16 bpc, 32 bpc floating point textures on modern hardware.
 */

#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "textureIndex=Screen('MakeTexture', WindowIndex, imageMatrix [, optimizeForDrawAngle=0] [, specialFlags=0] [, floatprecision=0] [, textureOrientation=0] [, textureShader=0]);";
//                                                            1            2              3                          4                5                    6                        7
static char synopsisString[] =
"Convert the 2D or 3D matrix 'imageMatrix' into an OpenGL texture and return an index which may be passed to 'DrawTexture' to specify the texture.\n"
"In the OpenGL Psychtoolbox textures replace offscreen windows for fast drawing of images during animation.\n"
"The imageMatrix argument may consist of one monochrome plane (Luminance), LA planes, RGB planes, or RGBA planes where "
"A is alpha, the transparency of a pixel. Alpha values typically range between zero (=fully transparent) and 255 (=fully opaque).\n"
"The Screen('ColorRange') command affects the range of expected input values in 'imageMatrix' matrices of double precision type, "
"as does the optional 'floatprecision' flag discussed below.\n"
"You need to enable Alpha-Blending via Screen('BlendFunction',...) for the transparency values to have an effect.\n"
"The argument 'optimizeForDrawAngle' if provided, asks Psychtoolbox to optimize the texture for especially fast "
"drawing at the specified rotation angle. The default is 0 == Optimize for upright drawing. If 'specialFlags' is set "
"to 1 and the width and height of the imageMatrix are powers of two (e.g., 64 x 64, 256 x 256, 512 x 512, ...), or "
"your graphics card supports so called non-power-of-two textures, then the texture is created as an OpenGL texture "
"of type GL_TEXTURE_2D. Otherwise Psychtoolbox will try to "
"pick the most optimal format for fast drawing and low memory consumption. GL_TEXTURE_2D textures are especially useful "
"for animation of drifting gratings, for simple use with the OpenGL 3D graphics functions and for blurring. Use of "
"GL_TEXTURE_2D textures is currently not automatically compatible with use of specialFlags settings 2 or 4.\n"
"If 'specialFlags' is set to 2 then PTB will try to use its own high quality texture filtering algorithm for drawing "
"of bilinearly filtered textures instead of the hardwares built-in method. This only works on modern hardware with "
"fragment shader support and is slower than using the hardwares built in filtering, but it may provide higher precision "
"on some hardware. PTB automatically enables its own filter algorithm when used with floating point "
"textures or when Screen('ColorRange') is used to enable unclamped color processing: PTB will check if your hardware "
"is capable of unrestricted high precision color processing in that case. If your hardware can't guarantee high precision, "
"PTB will enable its own shader-based workarounds to provide higher precision at the cost of lower speed. \n"
"If 'specialFlags' is set to 4 then PTB tries to use an especially fast method of texture creation. This method can be "
"at least an order of magnitude faster on some systems. However, it only works on modern GPUs, only for certain maximum "
"image sizes, and with some restrictions, e.g., scrolling of textures or high-precision filtering may not work at all or "
"as well. Your mileage may vary, so only use this flag if you need extra speed and after verifying your stimuli still look "
"correct. The biggest speedup is expected for creation of standard 8 bit integer textures from uint8 input matrices, "
"e.g., images from imread(), but also for 8 bit integer Luminance+Alpha and RGB textures from double format input matrices.\n"
"A 'specialFlags' == 8 will prevent automatic mipmap-generation for GL_TEXTURE_2D textures.\n"
"A 'specialFlags' == 32 setting will prevent automatic closing of the texture if Screen('Close'); is called. Only "
"Screen('Close', textureIndex); would close the texture.\n"
"'floatprecision' defines the precision with which the texture should be stored and processed. Default value is zero, "
"which asks to store textures with 8 bit per color component precision, a suitable format for standard images read via "
"imread(). A non-zero value will store the textures color component values as floating point precision numbers, useful "
"for complex blending operations and calculations on the textures and for processing and display of high dynamic range "
"image textures, either on a LDR display via tone-mapping, or on a HDR display. If floatprecision is set to 1, the "
"texture gets stored in half_float format, i.e. 16 bit per color component - Suitable for most display purposes and "
"fast on recent gfx-hardware. A value of 2 asks for full 32 bit single precision float per color component. Useful for complex "
"computations and image processing, but can be extremely slow when texture filtering is used on any piece of graphics hardware "
"manufactured before the year 2007. If a value of 1 is provided, asking for 16 bit floating point textures, but the graphics "
"hardware does not support this, then PTB tries to allocate a 15 bit precision signed integer texture instead, assuming the "
"graphics hardware supports that. Such a texture is more precise than the 16 bit floating point texture it replaces, but can "
"not store values outside the range [-1.0; 1.0]. On OpenGL-ES hardware, a 32 bit floating point texture is selected instead.\n"
"'textureOrientation' This optional argument labels textures with a special orientation. "
"Normally (value 0) a Matlab matrix is passed in standard Matlab column-major dataformat. This is efficient for drawing of "
"textures but not for processing them via Screen('TransformTexture'). Therefore textures need to be transformed on demand "
"if used that way. This flag allows to short-cut the process: A setting of 1 will ask for immediate conversion into the "
"optimized format. A setting of 2 will tell PTB that the Matlab matrix has been already converted into optimal format, "
"so no further processing is needed. A value of 3 tells PTB that the texture is completely isotropic, with no real orientation, "
"therefore no conversion is required. This latter setting only makes sense for random noise textures or other textures generated "
"from a distribution with uncorrelated noise-like pixels, e.g., some power spectrum distribution.\n"
"'textureShader' - optional: If you provide the handle of an OpenGL GLSL shader program, then this shader program will be "
"executed (bound) during drawing of this texture via the Screen('DrawTexture',...); command -- The normal texture drawing "
"operation is replaced by your customized algorithm. This is useful for two purposes: a) Very basic on-the-fly image processing "
"on the texture. b) Procedural shading: Your texture matrix doesn't encode an image, but only per-pixel parameters as input "
"for some formula to compute the real image during drawing. E.g., instead of defining a gabor patch as image or other standard "
"stimulus, one could define it as a mathematical formula to be evaluated at draw-time. The Screen('SetOpenGLTexture') command "
"allows you to create purely virtual textures which only consist of such a shader and some virtual size, but don't have any "
"real data matrix associated with it -- all content is generated on the fly.\n";

static char seeAlsoString[] = "DrawTexture TransformTexture BlendFunction";

PsychError SCREENMakeTexture(void)
{
    size_t                      ix, iters;
    PsychWindowRecordType       *textureRecord;
    PsychWindowRecordType       *windowRecord;
    PsychRectType               rect;
    psych_bool                  isImageMatrixBytes, isImageMatrixDoubles;
    int                         numMatrixPlanes, xSize, ySize;
    unsigned char               *byteMatrix;
    double                      *doubleMatrix;
    GLuint                      *texturePointer;
    GLubyte                     *texturePointer_b;
    GLfloat                     *texturePointer_f;
    double                      *rp, *gp, *bp, *ap;
    GLubyte                     *rpb, *gpb, *bpb, *apb;
    int                         usepoweroftwo, usefloatformat, assume_texorientation, textureShader;
    double                      optimized_orientation;
    psych_bool                  bigendian;
    psych_bool                  planar_storage = FALSE;
    double                      scaled = 1.0;
    double                      offsetd;

    // Detect endianity (byte-order) of machine:
    ix=255;
    rpb=(GLubyte*) &ix;
    bigendian = ( *rpb == 255 ) ? FALSE : TRUE;
    ix = 0; rpb = NULL;

    if(PsychPrefStateGet_DebugMakeTexture())    //MARK #1
        StoreNowTime();

    //all subfunctions should have these two lines.
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //Get the window structure for the onscreen window.  It holds the onscreein GL context which we will need in the
    //final step when we copy the texture from system RAM onto the screen.
    PsychErrorExit(PsychCapNumInputArgs(7));
    PsychErrorExit(PsychRequireNumInputArgs(2));
    PsychErrorExit(PsychCapNumOutputArgs(1));

    //get the window record from the window record argument and get info from the window record
    PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);

    if((windowRecord->windowType!=kPsychDoubleBufferOnscreen) && (windowRecord->windowType!=kPsychSingleBufferOnscreen))
        PsychErrorExitMsg(PsychError_user, "MakeTexture called on something else than a onscreen window");

    // Compute offset to add when converting double input matrices into uint8 textures, so that
    // the uint8 texture data matches how the hardware converts floating point color values to
    // fixed-point integer color values. When we convert a double float matrix into uint8 textures,
    // we want to apply the same conversion behavior that the hw applies when drawing "floating point
    // color" primitives into a fixed point framebuffer. This is important, e.g., for a uint8 texture
    // created from double 0.5 luminance values to draw the same color value into a 8 bpc framebuffer,
    // as a framebuffer clear command, drawdots command etc. with normalized color 0.5 does.
    // See what happens if this goes wrong (as of Psychtoolbox 3.0.14) in forum message #23300.
    //
    // The original OpenGL 1.0 spec says hw should convert float -> fixed with round-to-nearest
    // behaviour (Section 2.12.9). The OpenGL 4.5 spec leaves it up to the GL implementation if it
    // wants to truncate, round, ceil up or down to the closest integer etc., but recommends
    // round-to-nearest (Section 2.3.5.2). In practice, some tested AMD and NVidia hw on Linux and macOS
    // truncates, whereas some tested Intel hw on Linux and macOS rounds to nearest, so we apparently
    // can't rely on defined behaviour and therefore detect what the given hw does during startup.
    offsetd = (windowRecord->gfxcaps & kPsychGfxCapFloatToIntRound) ? 0.5 : 0.0;

    if (windowRecord->applyColorRangeToDoubleInputMakeTexture == 1) {
        // Apply scaling to uint8 textures which are created from double input matrices, so that
        // the Screen('ColorRange', ...) affects those values just as it affects other color specs for other
        // drawing commands:
        scaled = 255.0 / fabs(windowRecord->colorRange);
    }
    else {
        //  Do not apply scaling based on colorRange to double input for uint8 textures:
        scaled = 1.0;
    }

    // Get optional texture orientation flag:
    assume_texorientation = 0;
    PsychCopyInIntegerArg(6, FALSE, &assume_texorientation);

    // Get optional texture shader handle:
    textureShader = 0;
    PsychCopyInIntegerArg(7, FALSE, &textureShader);

    //get the argument and sanity check it.
    isImageMatrixBytes=PsychAllocInUnsignedByteMatArg(2, kPsychArgAnything, &ySize, &xSize, &numMatrixPlanes, &byteMatrix);
    isImageMatrixDoubles=PsychAllocInDoubleMatArg(2, kPsychArgAnything, &ySize, &xSize, &numMatrixPlanes, &doubleMatrix);

    if (!(isImageMatrixBytes || isImageMatrixDoubles))
        PsychErrorExitMsg(PsychError_user, "Illegal argument type. Image matrices must be uint8 or double data type.");

    if (numMatrixPlanes < 1 || numMatrixPlanes > 4) {
        iters = 0; // Make compiler happy.
        PsychErrorExitMsg(PsychError_inputMatrixIllegalDimensionSize, "Specified image matrix exceeds maximum depth of 4 layers");
    }

    if (ySize<1 || xSize <1) {
        iters = 0; // Make compiler happy.
        PsychErrorExitMsg(PsychError_inputMatrixIllegalDimensionSize, "Specified image matrix must be at least 1 x 1 pixels in size");
    }

    // Is this a special image matrix which is already pre-transposed to fit our optimal format?
    if (assume_texorientation == 2) {
        // Yes. Swap xSize and ySize to take this into account:
        ix = (size_t) xSize;
        xSize = ySize;
        ySize = (int) ix;
        ix = 0;
    }

    // Build defining rect for this texture:
    PsychMakeRect(rect, 0, 0, xSize, ySize);

    // Copy in texture preferred draw orientation hint. We default to zero degrees, if
    // not provided.
    // This parameter is not yet used. It is silently ignorerd for now...
    optimized_orientation = 0;
    PsychCopyInDoubleArg(3, FALSE, &optimized_orientation);

    // Copy in special creation mode flag: It defaults to zero. If set to 1 then we
    // always create a power-of-two GL_TEXTURE_2D texture. This is useful if one wants
    // to create and use drifting gratings with no effort - texture wrapping is only available
    // for GL_TEXTURE_2D, not for non-pot types. It is also useful if the texture is to be
    // exported to external OpenGL code to simplify tex coords assignments.
    usepoweroftwo=0;
    PsychCopyInIntegerArg(4, FALSE, &usepoweroftwo);

    // Check if size constraints are fullfilled for power-of-two mode:
    // We relax this constraint if GPU supports non-power-of-two texture extension.
    if ((usepoweroftwo & 1) && !(windowRecord->gfxcaps & kPsychGfxCapNPOTTex)) {
        for(ix = 1; ix < (size_t) xSize; ix*=2);
        if (ix != (size_t) xSize) {
            PsychErrorExitMsg(PsychError_inputMatrixIllegalDimensionSize, "Power-of-two texture requested but width of imageMatrix is not a power of two!");
        }

        for(ix = 1; ix < (size_t) ySize; ix*=2);
        if (ix != (size_t) ySize) {
            PsychErrorExitMsg(PsychError_inputMatrixIllegalDimensionSize, "Power-of-two texture requested but height of imageMatrix is not a power of two!");
        }
    }

    // Check if creation of a floating point texture is requested? We default to non-floating point,
    // standard 8 bpc textures if this parameter is not provided.
    usefloatformat = 0;
    PsychCopyInIntegerArg(5, FALSE, &usefloatformat);
    if (usefloatformat<0 || usefloatformat>2)
        PsychErrorExitMsg(PsychError_user, "Invalid value for 'floatprecision' parameter provided! Valid values are 0 for 8bpc int, 1 for 16bpc float or 2 for 32bpc float.");

    if (usefloatformat && !isImageMatrixDoubles) {
        // Floating point texture requested. We only support this if our input is a double matrix, not
        // for uint8 matrices - converting them to float precision would be just a waste of ressources
        // without any benefit for precision.
        PsychErrorExitMsg(PsychError_user, "Creation of a floating point precision texture requested, but uint8 matrix provided! Only double matrices are acceptable for this mode.");
    }

    // Float texture on OpenGL-ES requested?
    if ((usefloatformat > 0) && PsychIsGLES(windowRecord)) {
        // 32-bpc float textures supported? We can't do 16-bpc so we fail if 32 bpc is unsupported:
        if (!(windowRecord->gfxcaps & kPsychGfxCapFPTex32))
            PsychErrorExitMsg(PsychError_user, "Creation of a floating point precision texture requested, but this is not supported by your hardware!");

        // Upgrade float format to 2 aka 32 bpc float, the only thing we can handle:
        usefloatformat = 2;
    }

    //Create a texture record.  Really just a window record adapted for textures.
    PsychCreateWindowRecord(&textureRecord);
    textureRecord->windowType=kPsychTexture;

    // MK: We need to assign the screen number of the onscreen-window, so PsychCreateTexture()
    // can query the size of the screen/onscreen-window...
    textureRecord->screenNumber=windowRecord->screenNumber;
    textureRecord->depth=32;
    PsychCopyRect(textureRecord->rect, rect);

    // Is texture storage in planar format explicitely requested by usercode? Do the gpu and its size
    // constraints on textures support planar storage for this image?
    // Can a proper planar -> interleaved remapping GLSL shader be generated and assigned for this texture?
    if ((usepoweroftwo == 4) && (numMatrixPlanes > 1) && (windowRecord->gfxcaps & kPsychGfxCapFBO) && !(PsychPrefStateGet_ConserveVRAM() & kPsychDontCacheTextures) &&
        (ySize * numMatrixPlanes <= windowRecord->maxTextureSize) && PsychAssignPlanarTextureShaders(textureRecord, windowRecord, numMatrixPlanes)) {
        // Yes: Use the planar texture storage fast-path.
        planar_storage = TRUE;
        if (PsychPrefStateGet_Verbosity() > 6)
            printf("PTB-DEBUG: Using planar storage for %i layer texture of size %i x %i texels.\n", numMatrixPlanes, xSize, ySize);
    }
    else {
        planar_storage = FALSE;
        if (PsychPrefStateGet_Verbosity() > 7) printf("PTB-DEBUG: Using standard storage for %i layer texture of size %i x %i texels.\n", numMatrixPlanes, xSize, ySize);
    }

    //Allocate the texture memory and copy the MATLAB matrix into the texture memory.
    if (usefloatformat || (planar_storage && !isImageMatrixBytes)) {
        // Allocate a double for each color component and pixel:
        textureRecord->textureMemorySizeBytes = sizeof(double) * (size_t) numMatrixPlanes * (size_t) xSize * (size_t) ySize;
    }
    else {
        // Allocate one byte per color component and pixel:
        textureRecord->textureMemorySizeBytes = (size_t) numMatrixPlanes * (size_t) xSize * (size_t) ySize;
    }

    // We allocate our own intermediate conversion buffer unless this is
    // creation of a single-layer luminance8 integer texture from a single
    // layer uint8 input matrix and client storage is disabled. In that case, we can use a zero-copy path:
    if ((isImageMatrixBytes && (numMatrixPlanes == 1) && (!usefloatformat) && !(PsychPrefStateGet_ConserveVRAM() & kPsychDontCacheTextures)) ||
        (isImageMatrixBytes && planar_storage)) {
        // Zero copy path:
        texturePointer = NULL;
        // Set usefloatformat = 0 to prevent false compiler warnings about iters
        // being used uninitialized:
        usefloatformat = 0;
    }
    else {
        // Allocate memory:
        if(PsychPrefStateGet_DebugMakeTexture()) StoreNowTime();
        textureRecord->textureMemory = malloc(textureRecord->textureMemorySizeBytes);
        if(PsychPrefStateGet_DebugMakeTexture()) StoreNowTime();
        texturePointer = textureRecord->textureMemory;
    }

    // Does script explicitely request usage of a GL_TEXTURE_2D texture?
    if (usepoweroftwo & 1) {
        // Enforce creation as a GL_TEXTURE_2D texture:
        textureRecord->texturetarget=GL_TEXTURE_2D;
    }

            // Now the conversion routines that convert Matlab/Octave matrices into memory
            // buffers suitable for OpenGL:
    if (planar_storage) {
        // Planar texture storage, backed by a LUMINANCE texture container:

        // Zero-Copy possible? Only for uint8 input -> uint8 output:
        if (texturePointer == NULL) {
            texturePointer = (GLuint*) byteMatrix;
            textureRecord->textureMemory = texturePointer;
            // Set size to zero, so PsychCreateTexture() does not free() our
            // input buffer:
            textureRecord->textureMemorySizeBytes = 0;

            // This is always a LUMINANCE8 texture, backing our planar uint8 texture:
            textureRecord->depth = 8 * numMatrixPlanes;
            textureRecord->textureexternaltype   = GL_UNSIGNED_BYTE;
            textureRecord->textureexternalformat = GL_LUMINANCE;
            textureRecord->textureinternalformat = GL_LUMINANCE8;
        }
        else {
            // Some cast operation needed from double input format.
            // We always cast from double to float, potentially with
            // normalization and/or checking of value range.
            textureRecord->textureexternalformat = GL_LUMINANCE;

            if (usefloatformat) {
                // Floating point or other high precision format:
                textureRecord->depth = ((usefloatformat == 1) ? 16 : 32) * numMatrixPlanes;
                textureRecord->textureexternaltype = GL_FLOAT;
                textureRecord->textureinternalformat = (usefloatformat == 1) ? GL_LUMINANCE_FLOAT16_APPLE : GL_LUMINANCE_FLOAT32_APPLE;

                // Override for missing floating point texture support: Try to use 16 bit fixed point signed normalized textures [-1.0 ; 1.0] resolved at 15 bits:
                if ((usefloatformat == 1) && !(windowRecord->gfxcaps & kPsychGfxCapFPTex16)) textureRecord->textureinternalformat = GL_LUMINANCE16_SNORM;

                // Perform copy with double -> float cast:
                iters = (size_t) xSize * (size_t) ySize * (size_t) numMatrixPlanes;
                texturePointer_f = (GLfloat*) texturePointer;
                for(ix = 0; ix < iters; ix++) {
                    *(texturePointer_f++) = (GLfloat) *(doubleMatrix++);
                }
                iters = (size_t) xSize * (size_t) ySize;
            }
            else {
                // 8 Bit format, but from double input matrix -> cast to uint8:
                textureRecord->depth = 8 * numMatrixPlanes;
                textureRecord->textureexternaltype = GL_UNSIGNED_BYTE;
                textureRecord->textureinternalformat = GL_LUMINANCE8;

                iters = (size_t) xSize * (size_t) ySize * (size_t) numMatrixPlanes;
                texturePointer_b = (GLubyte*) texturePointer;
                for(ix = 0; ix < iters; ix++) {
                    *(texturePointer_b++) = (GLubyte) (offsetd + scaled * *(doubleMatrix++));
                }
                iters = (size_t) xSize * (size_t) ySize;
            }
        }
    }
    else if (usefloatformat) {
        // Conversion routines for HDR 16 bpc or 32 bpc textures -- Slow path.
        // Our input is always double matrices...
        iters = (size_t) xSize * (size_t) ySize;

        // Our input buffer is always of GL_FLOAT precision:
        textureRecord->textureexternaltype = GL_FLOAT;
        texturePointer_f=(GLfloat*) texturePointer;

        if (numMatrixPlanes==1) {
            for(ix=0;ix<iters;ix++){
                *(texturePointer_f++)= (GLfloat) *(doubleMatrix++);
            }
            textureRecord->depth=(usefloatformat==1) ? 16 : 32;

            textureRecord->textureinternalformat = (usefloatformat==1) ? GL_LUMINANCE_FLOAT16_APPLE : GL_LUMINANCE_FLOAT32_APPLE;
            textureRecord->textureexternalformat = GL_LUMINANCE;

            // Override for missing floating point texture support: Try to use 16 bit fixed point signed normalized textures [-1.0 ; 1.0] resolved at 15 bits:
            if ((usefloatformat==1) && !(windowRecord->gfxcaps & kPsychGfxCapFPTex16)) textureRecord->textureinternalformat = GL_LUMINANCE16_SNORM;
        }

        if (numMatrixPlanes==2) {
            rp=(double*) ((size_t) doubleMatrix);
            ap=(double*) ((size_t) rp + (size_t) iters * sizeof(double));

            for(ix=0;ix<iters;ix++){
                *(texturePointer_f++)= (GLfloat) *(rp++);
                *(texturePointer_f++)= (GLfloat) *(ap++);
            }
            textureRecord->depth=(usefloatformat==1) ? 32 : 64;
            textureRecord->textureinternalformat = (usefloatformat==1) ? GL_LUMINANCE_ALPHA_FLOAT16_APPLE : GL_LUMINANCE_ALPHA_FLOAT32_APPLE;
            textureRecord->textureexternalformat = GL_LUMINANCE_ALPHA;

            // Override for missing floating point texture support: Try to use 16 bit fixed point signed normalized textures [-1.0 ; 1.0] resolved at 15 bits:
            if ((usefloatformat==1) && !(windowRecord->gfxcaps & kPsychGfxCapFPTex16)) textureRecord->textureinternalformat = GL_LUMINANCE16_ALPHA16_SNORM;
        }

        if (numMatrixPlanes==3) {
            rp=(double*) ((size_t) doubleMatrix);
            gp=(double*) ((size_t) rp + (size_t) iters * sizeof(double));
            bp=(double*) ((size_t) gp + (size_t) iters * sizeof(double));

            for(ix=0;ix<iters;ix++){
                *(texturePointer_f++)= (GLfloat) *(rp++);
                *(texturePointer_f++)= (GLfloat) *(gp++);
                *(texturePointer_f++)= (GLfloat) *(bp++);
            }
            textureRecord->depth=(usefloatformat==1) ? 48 : 96;
            textureRecord->textureinternalformat = (usefloatformat==1) ? GL_RGB_FLOAT16_APPLE : GL_RGB_FLOAT32_APPLE;
            textureRecord->textureexternalformat = GL_RGB;

            // Override for missing floating point texture support: Try to use 16 bit fixed point signed normalized textures [-1.0 ; 1.0] resolved at 15 bits:
            if ((usefloatformat==1) && !(windowRecord->gfxcaps & kPsychGfxCapFPTex16)) textureRecord->textureinternalformat = GL_RGB16_SNORM;
        }

        if (numMatrixPlanes==4) {
            rp=(double*) ((size_t) doubleMatrix);
            gp=(double*) ((size_t) rp + (size_t) iters * sizeof(double));
            bp=(double*) ((size_t) gp + (size_t) iters * sizeof(double));
            ap=(double*) ((size_t) bp + (size_t) iters * sizeof(double));

            for(ix=0;ix<iters;ix++){
                *(texturePointer_f++)= (GLfloat) *(rp++);
                *(texturePointer_f++)= (GLfloat) *(gp++);
                *(texturePointer_f++)= (GLfloat) *(bp++);
                *(texturePointer_f++)= (GLfloat) *(ap++);
            }
            textureRecord->depth=(usefloatformat==1) ? 64 : 128;
            textureRecord->textureinternalformat = (usefloatformat==1) ? GL_RGBA_FLOAT16_APPLE : GL_RGBA_FLOAT32_APPLE;
            textureRecord->textureexternalformat = GL_RGBA;

            // Override for missing floating point texture support: Try to use 16 bit fixed point signed normalized textures [-1.0 ; 1.0] resolved at 15 bits:
            if ((usefloatformat==1) && !(windowRecord->gfxcaps & kPsychGfxCapFPTex16)) textureRecord->textureinternalformat = GL_RGBA16_SNORM;
        }
        // End of HDR conversion code...
    }
    else {
        // Standard LDR texture 8 bpc conversion routines -- Fast path.
        iters = (size_t) xSize * (size_t) ySize;

        // Improved implementation: Takes 13 ms on a 800x800 texture...
        if (isImageMatrixDoubles && numMatrixPlanes==1){
            texturePointer_b=(GLubyte*) texturePointer;
            for(ix=0;ix<iters;ix++){
                *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(doubleMatrix++));
            }
            textureRecord->depth=8;
        }

        // Improved version: Takes 3 ms on a 800x800 texture...
        // NB: Implementing memcpy manually by a for-loop takes 10 ms! This is a huge difference.
        // -> That's because memcpy on MacOS-X is implemented with hand-coded, highly tuned Assembler code for PowerPC.
        // -> It's always wise to use system-routines if available, instead of coding it by yourself!
        if (isImageMatrixBytes && numMatrixPlanes==1) {
            if (texturePointer) {
                // Need to do a copy. Use optimized memcpy():
                memcpy((void*) texturePointer, (void*) byteMatrix, iters);

                //texturePointer_b=(GLubyte*) texturePointer;
                //for(ix=0;ix<iters;ix++){
                //    *(texturePointer_b++) = *(byteMatrix++);
                //}
            }
            else {
                // Zero-Copy path. Just pass a pointer to our input matrix:
                texturePointer = (GLuint*) byteMatrix;
                textureRecord->textureMemory = texturePointer;
                // Set size to zero, so PsychCreateTexture() does not free() our
                // input buffer:
                textureRecord->textureMemorySizeBytes = 0;
            }

            textureRecord->depth=8;
        }

        // New version: Takes 33 ms on a 800x800 texture...
        if (isImageMatrixDoubles && numMatrixPlanes==2){
            texturePointer_b=(GLubyte*) texturePointer;
            rp=(double*) ((size_t) doubleMatrix);
            ap=(double*) ((size_t) rp + (size_t) iters * sizeof(double));
            for(ix=0;ix<iters;ix++){
                *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(rp++));
                *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(ap++));
            }
            textureRecord->depth=16;
        }

        // New version: Takes 20 ms on a 800x800 texture...
        if (isImageMatrixBytes && numMatrixPlanes==2){
            texturePointer_b=(GLubyte*) texturePointer;
            rpb=(GLubyte*) ((size_t) byteMatrix);
            apb=(GLubyte*) ((size_t) rpb + (size_t) iters);
            for(ix=0;ix<iters;ix++){
                *(texturePointer_b++)= *(rpb++);
                *(texturePointer_b++)= *(apb++);
            }
            textureRecord->depth=16;
        }

        // Improved version: Takes 43 ms on a 800x800 texture...
        if (isImageMatrixDoubles && numMatrixPlanes==3){
            texturePointer_b=(GLubyte*) texturePointer;
            rp=(double*) ((size_t) doubleMatrix);
            gp=(double*) ((size_t) rp + (size_t) iters * sizeof(double));
            bp=(double*) ((size_t) gp + (size_t) iters * sizeof(double));
            for(ix=0;ix<iters;ix++){
                *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(rp++));
                *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(gp++));
                *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(bp++));
            }
            textureRecord->depth=24;
        }

        // Improved version: Takes 25 ms on a 800x800 texture...
        if (isImageMatrixBytes && numMatrixPlanes==3){
            texturePointer_b=(GLubyte*) texturePointer;
            rpb=(GLubyte*) ((size_t) byteMatrix);
            gpb=(GLubyte*) ((size_t) rpb + (size_t) iters);
            bpb=(GLubyte*) ((size_t) gpb + (size_t) iters);
            for(ix=0;ix<iters;ix++){
                *(texturePointer_b++)= *(rpb++);
                *(texturePointer_b++)= *(gpb++);
                *(texturePointer_b++)= *(bpb++);
            }
            textureRecord->depth=24;
        }

        // Improved version: Takes 55 ms on a 800x800 texture...
        if (isImageMatrixDoubles && numMatrixPlanes==4){
            texturePointer_b=(GLubyte*) texturePointer;
            rp=(double*) ((size_t) doubleMatrix);
            gp=(double*) ((size_t) rp + (size_t) iters * sizeof(double));
            bp=(double*) ((size_t) gp + (size_t) iters * sizeof(double));
            ap=(double*) ((size_t) bp + (size_t) iters * sizeof(double));
            if (bigendian) {
                // Code for big-endian machines like PowerPC:
                for(ix=0;ix<iters;ix++){
                    *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(ap++));
                    *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(rp++));
                    *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(gp++));
                    *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(bp++));
                }
            }
            else {
                // Code for little-endian machines like Intel Pentium:
                for(ix=0;ix<iters;ix++){
                    *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(bp++));
                    *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(gp++));
                    *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(rp++));
                    *(texturePointer_b++)= (GLubyte) (offsetd + scaled * *(ap++));
                }
            }

            textureRecord->depth=32;
        }

        // Improved version: Takes 33 ms on a 800x800 texture...
        if (isImageMatrixBytes && numMatrixPlanes==4){
            texturePointer_b=(GLubyte*) texturePointer;
            rpb=(GLubyte*) ((size_t) byteMatrix);
            gpb=(GLubyte*) ((size_t) rpb + (size_t) iters);
            bpb=(GLubyte*) ((size_t) gpb + (size_t) iters);
            apb=(GLubyte*) ((size_t) bpb + (size_t) iters);
            if (bigendian) {
                // Code for big-endian machines like PowerPC:
                for(ix=0;ix<iters;ix++){
                    *(texturePointer_b++)= *(apb++);
                    *(texturePointer_b++)= *(rpb++);
                    *(texturePointer_b++)= *(gpb++);
                    *(texturePointer_b++)= *(bpb++);
                }
            }
            else {
                // Code for little-endian machines like Intel Pentium:
                for(ix=0;ix<iters;ix++){
                    *(texturePointer_b++)= *(bpb++);
                    *(texturePointer_b++)= *(gpb++);
                    *(texturePointer_b++)= *(rpb++);
                    *(texturePointer_b++)= *(apb++);
                }
            }

            textureRecord->depth=32;
        }
    } // End of 8 bpc texture conversion code (fast-path for LDR textures)

    // Override for missing floating point texture support?
    if ((usefloatformat==1) && !(windowRecord->gfxcaps & kPsychGfxCapFPTex16)) {
        // Override enabled. Instead of a 16bpc float texture with 11 bit linear precision in the
        // range [-1.0 ; 1.0], we use a 16 bit signed normalized texture with a normalized value
        // range of [-1.0; 1.0], encoded with 1 bit sign and 15 bit magnitude. These textures have
        // an effective linear precision of 15 bits - better than 16 bpc float - but they are restricted
        // to a value range of [-1.0 ; 1.0], as opposed to 16 bpc float textures. Tell user about this
        // replacement at high verbosity levels:
        if (PsychPrefStateGet_Verbosity() > 4)
            printf("PTB-INFO:MakeTexture: Code requested 16 bpc float texture, but this is unsupported. Trying to use 16 bit snorm texture instead.\n");

        // Signed normalized textures supported? Otherwise we bail...
        if (!(windowRecord->gfxcaps & kPsychGfxCapSNTex16)) {
            printf("PTB-ERROR:MakeTexture: Code requested 16 bpc floating point texture, but this is unsupported by this graphics card.\n");
            printf("PTB-ERROR:MakeTexture: Tried to use 16 bit snorm texture instead, but failed as this is unsupported as well.\n");
            PsychErrorExitMsg(PsychError_user, "Creation of 15 bit linear precision signed normalized texture failed. Not supported by your graphics hardware!");
        }

        // Check value range of pixels. This will not work for out of [-1; 1] range values.
        texturePointer_f=(GLfloat*) texturePointer;
        iters = iters * (size_t) numMatrixPlanes;
        for (ix=0; ix<iters; ix++, texturePointer_f++) {
            if(fabs((double) *texturePointer_f) > 1.0) {
                // Game over!
                printf("PTB-ERROR:MakeTexture: Code requested 16 bpc floating point texture, but this is unsupported by this graphics card.\n");
                printf("PTB-ERROR:MakeTexture: Tried to use 16 bit snorm texture instead, but failed because some pixels are outside the\n");
                printf("PTB-ERROR:MakeTexture: representable range -1.0 to 1.0 for this texture type. Change your code or update your graphics hardware.\n");
                PsychErrorExitMsg(PsychError_user, "Creation of 15 bit linear precision signed normalized texture failed due to out of [-1 ; +1] range pixel values!");
            }
        }
    }

    // This is a special workaround for bugs in FLOAT16 texture creation on Mac OS/X 10.4.x and 10.5.x.
    // The OpenGL fails to properly flush very small values (< 1e-9) to zero when creating a FLOAT16
    // type texture. Instead it seems to initialize with trash data, corrupting the texture.
    // Therefore, if FLOAT16 texture creation is requested, we loop over the whole input buffer and
    // set all values with magnitude smaller than 1e-9 to zero. Better safe than sorry...
    if ((usefloatformat==1) && (windowRecord->gfxcaps & kPsychGfxCapFPTex16)) {
        texturePointer_f=(GLfloat*) texturePointer;
        iters = iters * (size_t) numMatrixPlanes;
        for(ix=0; ix<iters; ix++, texturePointer_f++) if(fabs((double) *texturePointer_f) < 1e-9) { *texturePointer_f = 0.0; }
    }

    // On OpenGL-ES, 32 bpc floating point textures are selected via the GL_FLOAT type specifier, and
    // internal format must be == external format == not defining resolution. External format is already
    // properly set by common desktop/es HDR setup code, as is type spec, so we just need to make sure that
    // internal format is consistent with external one:
    // MK NOPE: Seems i misread the spec and float textures are treated identical to desktop OpenGL, so this
    // is probably not needed, at least not on NVidia. Leave it here in case this is a NVidia peculiarity and
    // my old interpretation was actually correct for any other hardware.
    //    if ((usefloatformat == 2) && PsychIsGLES(windowRecord)) textureRecord->textureinternalformat = textureRecord->textureexternalformat;

    // The memory buffer now contains our texture data in a format ready to submit to OpenGL.

    // Assign parent window and copy its inheritable properties:
    PsychAssignParentWindow(textureRecord, windowRecord);

    // Texture orientation is zero aka transposed aka non-renderswapped.
    textureRecord->textureOrientation = ((assume_texorientation != 2) && (assume_texorientation != 3)) ? 0 : 2;

    // This is our best guess about the number of image channels:
    textureRecord->nrchannels = numMatrixPlanes;

    if (planar_storage) {
        // Setup special rect to fake PsychCreateTexture() into creating a luminance
        // texture numMatrixPlanes times the height (in rows) of the texture, to store the
        // numMatrixPlanes layers concatenated to each other.
        if (textureRecord->textureOrientation == 0) {
            // Normal case: Transposed storage.
            PsychMakeRect(&(textureRecord->rect[0]), 0, 0, xSize * numMatrixPlanes, ySize);
        }
        else {
            // Special case: Non-transposed or isotropic storage:
            PsychMakeRect(&(textureRecord->rect[0]), 0, 0, xSize, ySize * numMatrixPlanes);
        }

        // Create planar texture:
        PsychCreateTexture(textureRecord);

        // Restore rect and clientrect of texture to effective size:
        PsychMakeRect(&(textureRecord->rect[0]), 0, 0, xSize, ySize);
        PsychCopyRect(textureRecord->clientrect, textureRecord->rect);

        textureRecord->specialflags = kPsychPlanarTexture;
    }
    else {
        // Let's create and bind a new texture object and fill it with our new texture data.
        PsychCreateTexture(textureRecord);

        // Assign GLSL filter-/lookup-shaders if needed:
        PsychAssignHighPrecisionTextureShaders(textureRecord, windowRecord, usefloatformat, (usepoweroftwo & 2) ? 1 : 0);
    }

    // User specified override shader for this texture provided? This is useful for
    // basic image processing and procedural texture shading:
    if (textureShader!=0) {
        // Assign provided shader as filtershader to this texture: We negate it so
        // that the texture blitter routines know this is a custom shader, not our
        // built in filter shader:
        textureRecord->textureFilterShader = -1 * textureShader;
    }

    // Texture ready. Mark it valid and return handle to userspace:
    PsychSetWindowRecordValid(textureRecord);
    PsychCopyOutDoubleArg(1, FALSE, textureRecord->windowIndex);

    // Swapping the texture to upright orientation requested?
    if (assume_texorientation == 1) {
        // Transform sourceRecord source texture into a normalized, upright texture if it isn't already in
        // that format. We require this standard orientation for simplified shader design.
        PsychSetShader(windowRecord, 0);
        PsychNormalizeTextureOrientation(textureRecord);
    }

    // Shall the texture be finally declared "normally oriented"?
    // This is either due to explicit renderswapping if assume_textureorientation == 1,
    // or because it was already pretransposed in Matlab/Octave if assume_textureorientation == 2,
    // or because user space tells us the texture is isotropic if assume_textureorientation == 3.
    if (assume_texorientation > 0) {
        // Yes. Label it as such:
        textureRecord->textureOrientation = 2;
    }

    // specialFlags setting 8? Disable auto-mipmap generation:
    if (usepoweroftwo & 0x8) textureRecord->specialflags |= kPsychDontAutoGenMipMaps;

    // A specialFlags setting of 32? Protect texture against deletion via Screen('Close') without providing a explicit handle:
    if (usepoweroftwo & 32) textureRecord->specialflags |= kPsychDontDeleteOnClose;

    if (PsychPrefStateGet_DebugMakeTexture())     //MARK #4
        StoreNowTime();

    return(PsychError_none);
}
