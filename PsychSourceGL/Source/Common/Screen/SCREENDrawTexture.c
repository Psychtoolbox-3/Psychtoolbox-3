/*
 *    SCREENDrawTexture.c
 *
 *    AUTHORS:
 *
 *    Allen.Ingling@nyu.edu               awi
 *    mario.kleiner.de@gmail.com          mk
 *
 *    PLATFORMS:
 *
 *    All
 *
 *    HISTORY:
 *
 *        6/25/04         awi     Created.
 *        10/12/04        awi     In useString: changed "SCREEN" to "Screen", and moved commas to inside [].
 *        1/25/05         awi     Added updated provided by mk
 *                                - Screen-sized textures work correclty
 *                                - Textures larger than the screen are rejected
 *                                - Fixed sourerect parameter
 *        2/25/05         awi     Added call to PsychSetGLContext
 *                                Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
 *        4/22/05         mk      Commented out glFinish() - This would be a performance killer!
 *        5/13/05         mk      Support for rotated drawing of textures.
 *        7/23/05         mk      New options filterMode and globalAlpha.
 *        9/30/05         mk      Remove size check for texturesize <= windowsize. This restriction doesn't apply anymore for new texture mapping code.
 *
 */

#include "Screen.h"

static char seeAlsoString[] = "MakeTexture DrawTexture DrawTextures";

PsychError SCREENDrawTexture(void)
{
    static char synopsisString[] =
    "Draw the texture specified via 'texturePointer' into the target window specified via 'windowPointer'. "
    "In the the OS X Psychtoolbox textures replace offscreen windows for fast drawing of images during animation."
    "'sourceRect' specifies a rectangular subpart of the texture to be drawn (Defaults to full texture). "
    "'destinationRect' defines the rectangular subpart of the window where the texture should be drawn. This defaults"
    "to centered on the screen. "
    "'rotationAngle' Specifies a rotation angle in degree for rotated drawing of the texture (Defaults to 0 deg. = upright). "
    "'filterMode' How to compute the pixel color values when the texture is drawn magnified, minified or drawn shifted, e.g., "
    "if sourceRect and destinationRect do not have the same size or if sourceRect specifies fractional pixel values. 0 = Nearest "
    "neighbour filtering, 1 = Bilinear filtering - this is the default. Values 2 or 3 select use of OpenGL mip-mapping for improved "
    "quality: 2 = Bilinear filtering for nearest mipmap level, 3 = Trilinear filtering across mipmap levels, 4 = Nearest neighbour "
    "filtering for nearest mipmap level, 5 = nearest neighbour filtering with linear interpolation between mipmap levels. Mipmap filtering is "
    "only supported for GL_TEXTURE_2D textures (see description of 'specialFlags' flag 1 below). A negative filterMode value will "
    "also use mip-mapping for fast drawing of blurred textures if the GL_TEXTURE_2D format is used: Mip-maps are essentially image "
    "resolution pyramids, the filterMode value selects a specific layer in that pyramid. A value of -1 draws the highest resolution "
    "layer, a value of -2 draws a half-resolution layer, a value of -3 draws a quarter resolution layer and so on. Each layer has "
    "half the resolution of the preceeding layer. This allows for very fast drawing of blurred or low-pass filtered images, e.g., for "
    "gaze-contingent displays. However, the filter function for downsampling is system dependent and may vary across graphics cards, "
    "although a box-filter is the most common type. If you need a well defined filter function, use a custom written GLSL shader "
    "instead, so you have full control over the mathematical properties of the downsampling function. This would incur a "
    "performance penalty.\n"
    "'globalAlpha' A global alpha transparency value to apply "
    "to the whole texture for blending. Range is 0 = fully transparent to 1 = fully opaque, defaults to one. If both, an alpha-channel "
    "and globalAlpha are provided, then the final alpha is the product of both values. 'modulateColor', if provided, overrides the "
    "'globalAlpha' value. If 'modulateColor' is specified, the 'globalAlpha' value will be ignored. 'modulateColor' will be a global "
    "color that gets applied to the texture as a whole, i.e., it modulates each color channel. E.g., modulateColor = [128 255 0] would "
    "leave the green- and alpha-channel untouched, but it would multiply the blue channel with 0 - set it to zero blue intensity, and "
    "it would multiply each texel in the red channel by 128/255 - reduce its intensity to 50%. The most interesting application of "
    "'modulateColor' is drawing of arbitrary complex shapes of selectable color: Simply generate an all-white luminance texture of "
    "arbitrary shape, possibly with alpha channel, then draw it with 'modulateColor' set to the wanted color and global alpha value.\n"
    "'textureShader' (optional): If you provide a valid handle of a GLSL shader, this shader will be applied to the texture during "
    "drawing. If the texture already has a shader assigned (via Screen('MakeTexture') or automatically by PTB for some reason), then "
    "the shader provided here as 'textureShader' will silently override the shader assigned earlier. Application of shaders this way "
    "is mostly useful for application of simple single-pass image processing operations to a texture, e.g., a simple blur or a "
    "deinterlacing operation for a video texture. If you intend to use this texture multiple times or if you need more complex image "
    "processing, e.g., multi-pass operations, better use the Screen('TransformTexture') command. It allows for complex operations to "
    "be applied and is more flexible.\n"
    "'specialFlags' optional argument: Allows to pass a couple of special flags to influence the drawing. The flags can be combined "
    "by mor() ing them together. A value of kPsychUseTextureMatrixForRotation will use a different mode of operation for drawing of "
    "rotated textures, where the drawn 'dstRect' texture rectangle is always upright, but texels are retrieved at rotated positions, "
    "as if the 'srcRect' rectangle would be rotated. If you set a value of kPsychDontDoRotation then the rotation angle will not be "
    "used to rotate the texture. Instead it will be passed to a bount texture shader (if any), which is free to interpret the "
    "'rotationAngle' parameters is it wants - e.g., to implement custom texture rotation."
    "\n\n"
    "'auxParameters' optional argument: If this is set as a vector with at least 4 components, and a multiple of four components, "
    "then these values are passed to a shader (if any is bound) as 'auxParameter0....n'. The current implementation supports at "
    "most 32 values per draw call. This is mostly useful when drawing procedural textures if one needs to pass more additional "
    "parameters to define the texture than can fit into other parameter fields. See 'help ProceduralShadingAPI' for more info. "
    "\n\n"
    "If you want to draw many textures to the same onscreen- or offscreen window, use the function Screen('DrawTextures'). "
    "It accepts the same arguments as this function, but is optimized to draw many textures in one call.";

    // If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
    static char useString[] = "Screen('DrawTexture', windowPointer, texturePointer [,sourceRect] [,destinationRect] [,rotationAngle] [, filterMode] [, globalAlpha] [, modulateColor] [, textureShader] [, specialFlags] [, auxParameters]);";
    //                                               1              2                3             4                  5                 6              7               8                 9                 10               11

    PsychWindowRecordType *source, *target;
    PsychRectType sourceRect, targetRect, tempRect;
    double rotationAngle = 0;   // Default rotation angle is zero deg. = upright.
    int filterMode = 1;         // Default filter mode is bilinear filtering.
    double globalAlpha = 1.0;   // Default global alpha is 1 == no effect.
    PsychColorType color;
    int textureShader, backupShader;
    double* auxParameters;
    int numAuxParams, m, n, p;
    psych_bool isclassic;
    int specialFlags = 0;

    //all subfunctions should have these two lines.
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) { PsychGiveHelp();return(PsychError_none); };

    //Get the window structure for the onscreen window.  It holds the onscreein GL context which we will need in the
    //final step when we copy the texture from system RAM onto the screen.
    PsychErrorExit(PsychCapNumInputArgs(11));
    PsychErrorExit(PsychRequireNumInputArgs(2));
    PsychErrorExit(PsychCapNumOutputArgs(0));

    //Read in arguments
    PsychAllocInWindowRecordArg(1, kPsychArgRequired, &target);
    PsychAllocInWindowRecordArg(2, kPsychArgRequired, &source);
    if(source->windowType!=kPsychTexture) {
        PsychErrorExitMsg(PsychError_user, "The first argument supplied was a window pointer, not a texture pointer");
    }

    // Classic OpenGL-1/2?
    isclassic = PsychIsGLClassic(target);

    PsychCopyRect(sourceRect,source->clientrect);
    PsychCopyInRectArg(3, kPsychArgOptional, sourceRect);
    if (IsPsychRectEmpty(sourceRect)) return(PsychError_none);

    PsychCopyRect(tempRect, target->clientrect);
    PsychCenterRectInRect(sourceRect, tempRect, targetRect);
    PsychCopyInRectArg(4, kPsychArgOptional, targetRect);
    if (IsPsychRectEmpty(targetRect)) return(PsychError_none);

    PsychCopyInDoubleArg(5, kPsychArgOptional, &rotationAngle);
    PsychCopyInIntegerArg(6, kPsychArgOptional, &filterMode);
    if (filterMode > 5) {
        PsychErrorExitMsg(PsychError_user, "filterMode needs to be negative for a specific blur level, or at most 5 for other modes.");
    }

    // Copy in optional 'globalAlpha': We don't put restrictions on its valid range
    // anymore - That made sense for pure fixed function LDR rendering, but no longer
    // for HDR rendering or procedural shading.
    PsychCopyInDoubleArg(7, kPsychArgOptional, &globalAlpha);

    PsychSetDrawingTarget(target);
    PsychUpdateAlphaBlendingFactorLazily(target);

    if(PsychCopyInColorArg(8, kPsychArgOptional, &color)) {
        // set globalAlpha to DBL_MAX to signal that PsychBlitTexture() shouldn't
        // use this parameter and not set any modulate color, cause we do it.
        globalAlpha = DBL_MAX;

        // Setup global vertex color as modulate color for texture drawing:
        PsychCoerceColorMode(&color);

        // This call stores unclamped color in target->currentColor, as needed
        // if color is to be processed by some bound shader (procedural or filtershader)
        // inside PsychBlitTextureToDisplay():
        PsychConvertColorToDoubleVector(&color, target, (GLdouble*) &(target->currentColor));

        // Submit the same color to fixed function pipe attribute as well, in case no
        // shader is bound, or shader pulls from standard color attribute (we can't know yet):
        if (isclassic) {
            glColor4dv(target->currentColor);
        }
        else {
            PsychGLColor4f(target, (float) target->currentColor[0], (float) target->currentColor[1], (float) target->currentColor[2], (float) target->currentColor[3]);
        }
    }

    // Assign optional override texture shader, if any provided:
    textureShader = -1;
    PsychCopyInIntegerArg(9, kPsychArgOptional, &textureShader);

    // Assign any other optional special flags:
    PsychCopyInIntegerArg(10, kPsychArgOptional, &specialFlags);

    // Set rotation mode flag for texture matrix rotation if secialFlags is set accordingly:
    if (specialFlags & kPsychUseTextureMatrixForRotation) source->specialflags|=kPsychUseTextureMatrixForRotation;
    // Set rotation mode flag for no fixed function pipeline rotation if secialFlags is set accordingly:
    if (specialFlags & kPsychDontDoRotation) source->specialflags|=kPsychDontDoRotation;

    // Optional auxParameters:
    auxParameters = NULL;
    m=n=p=0;
    if (PsychAllocInDoubleMatArg(11, kPsychArgOptional, &m, &n, &p, &auxParameters)) {
        if ((p!=1) || (m * n < 4) || (((m*n) % 4)!=0)) PsychErrorExitMsg(PsychError_user, "The 11th argument must be a vector of 'auxParameter' values with a multiple of 4 components.");
    }
    numAuxParams = m*n;
    target->auxShaderParamsCount = numAuxParams;

    // Pass auxParameters for current primitive in the auxShaderParams field.
    if (numAuxParams > 0) {
        target->auxShaderParams = auxParameters;
    }
    else {
        target->auxShaderParams = NULL;
    }

    if (textureShader > -1) {
        backupShader = source->textureFilterShader;
        source->textureFilterShader = -1 * textureShader;
        PsychBlitTextureToDisplay(source, target, sourceRect, targetRect, rotationAngle, filterMode, globalAlpha);
        source->textureFilterShader = backupShader;
    }
    else {
        PsychBlitTextureToDisplay(source, target, sourceRect, targetRect, rotationAngle, filterMode, globalAlpha);
    }

    // Reset rotation mode flag:
    source->specialflags &= ~(kPsychUseTextureMatrixForRotation | kPsychDontDoRotation);

    target->auxShaderParams = NULL;
    target->auxShaderParamsCount = 0;

    // Mark end of drawing op. This is needed for single buffered drawing:
    PsychFlushGL(target);

    return(PsychError_none);

}

// Batch-drawing version of DrawTexture:
PsychError SCREENDrawTextures(void)
{
    // If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c 1 2 3 4 5 6 7 8
    static char useString[] = "Screen('DrawTextures', windowPointer, texturePointer(s) [, sourceRect(s)] [, destinationRect(s)] [, rotationAngle(s)] [, filterMode(s)] [, globalAlpha(s)] [, modulateColor(s)] [, textureShader] [, specialFlags] [, auxParameters]);";
    //                                                1              2                    3                 4                      5                    6                 7                  8                    9                 10               11

    static char synopsisString[] = "Draw many textures at once, either one texture to many locations or many textures.\n"
    "This function accepts the same parameters as Screen('DrawTexture'), but it is optimized for drawing many textures. "
    "You can leave out each argument, a default setting will be used in that case, provide it once to apply it to all "
    "drawn items, or provide a vector or matrix with a individual setting for each drawn item. If you provide multiple "
    "settings per argument, then the number must match between all arguments.\n\n"
    "Examples:\n"
    "a) One texture drawn to different locations at different orientations: Provide one texture handle for the texturePointer, "
    "a 4 row by n columns matrix for 'destinationRect' to provide target rectangles for n locations, provide a n component "
    "vector of 'rotationAngles' for the n different orientations of the n drawn texture patches.\n"
    "b) n textures drawn to n different locations: Same as a) but provide a n component vector of 'texturePointers' one for "
    "each texture to be drawn to one of n locations at n angles.\n";

    PsychWindowRecordType *source, *target;
    PsychRectType sourceRect, targetRect, tempRect;
    PsychColorType color;
    double *dstRects, *srcRects, *colors, *penSizes, *globalAlphas, *filterModes, *rotationAngles;
    unsigned char *bytecolors;
    int numTexs, numdstRects, numsrcRects, i, nc, mc, nrsize, m, n, p, numAngles, numFilterModes, numAlphas, numRef;
    double* texids;
    double rotationAngle, globalAlpha, filterMode;
    double* auxParameters;
    int numAuxParams, numAuxComponents;
    psych_bool isclassic;
    int textureShader, backupShader = 0;
    int specialFlags = 0;
    psych_bool batchIt = FALSE;

    //all subfunctions should have these two lines.
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){ PsychGiveHelp();return(PsychError_none); };

    //Get the window structure for the onscreen window.  It holds the onscreen GL context which we will need in the
    //final step when we copy the texture from system RAM onto the screen.
    PsychErrorExit(PsychCapNumInputArgs(11));
    PsychErrorExit(PsychRequireNumInputArgs(2));
    PsychErrorExit(PsychCapNumOutputArgs(0));

    // The target window is a fixed parameter:
    PsychAllocInWindowRecordArg(1, kPsychArgRequired, &target);

    // Classic OpenGL-1/2?
    isclassic = PsychIsGLClassic(target);

    // First get all source texture handles:
    PsychAllocInDoubleMatArg(2, kPsychArgRequired, &m, &n, &p, &texids);
    if ((p!=1) || (m>1 && n!=1) || (n>1 && m!=1)) PsychErrorExitMsg(PsychError_user, "The second argument must be either a row- or columnvector of valid texture handles.");

    // This is the number of texture handles:
    numTexs = m * n;

    // Only one texture?
    if (numTexs == 1) {
        // Yes. Allocate it in the conventional way:
        PsychAllocInWindowRecordArg(2, kPsychArgRequired, &source);
        if(source->windowType!=kPsychTexture) {
            PsychErrorExitMsg(PsychError_user, "The second argument supplied was not a texture handle!");
        }
    }

    // Query, allocate and copy in all vectors...
    numdstRects = 4;
    nrsize = 0;
    colors = NULL;
    bytecolors = NULL;
    penSizes = NULL;

    // The negative position -4 means: dstRects coords are expected at position 4, but they are optional.
    // NULL means - don't want a size's vector.
    PsychPrepareRenderBatch(target, -4, &numdstRects, &dstRects, 8, &nc, &mc, &colors, &bytecolors, 5, &nrsize, &penSizes, FALSE);

    // At this point, target is set up as target window, i.e. its GL-Context is active, it is set as drawing target,
    // alpha blending is set up according to Screen('BlendFunction'), and the drawing color is set if it is a singular one.
    if (nc <= 1) {
        // Only one - or no - color provided. One or none?
        if(PsychCopyInColorArg(8, kPsychArgOptional, &color)) {
            // One global modulate color provided:

            // Setup global vertex color as modulate color for texture drawing:
            PsychCoerceColorMode(&color);
            PsychSetGLColor(&color, target);
        }
        else {
            // No modulateColor provided: Don't use this parameter:
            nc = 0;
        }
    }

    // Try to get source rects:
    m=n=p=0;
    if (PsychAllocInDoubleMatArg(3, kPsychArgOptional, &m, &n, &p, &srcRects)) {
        if ((p!=1) || (m!=1 && m!=4)) PsychErrorExitMsg(PsychError_user, "The third argument must be either empty, or a single srcRect 4 component row vector, or a 4 row by n column matrix with srcRects for all objects to draw, not a 3D matrix!");
        // Ok, its a one row or four row matrix:
        if (m==4) {
            // Potentially multiple source rects provided:
            numsrcRects = n;
        }
        else {
            // Its a one row vector: This is either a single srcRect for all textures, or something invalid:
            if (n!=4) PsychErrorExitMsg(PsychError_user, "The third argument must be either empty, or a single srcRect 4 component row vector, or a 4 row by n column matrix with srcRects for all objects to draw!");
            // Single srcRect provided:
            numsrcRects = 1;
        }
    }
    else {
        // No srcRects provided:
        numsrcRects = 0;
    }

    // Optional rotation angles:
    m=n=p=0;
    if (PsychAllocInDoubleMatArg(5, kPsychArgOptional, &m, &n, &p, &rotationAngles)) {
        if ((p!=1) || (m>1 && n!=1) || (n>1 && m!=1)) PsychErrorExitMsg(PsychError_user, "The fifth argument must be either a row- or columnvector of rotation angles.");
    }
    numAngles = m * n;

    // Default to 0 degree rotation -- upright drawing:
    rotationAngle = (numAngles == 1) ? rotationAngles[0] : 0.0;

    // Optional filter modes:
    m=n=p=0;
    if (PsychAllocInDoubleMatArg(6, kPsychArgOptional, &m, &n, &p, &filterModes)) {
        if ((p!=1) || (m>1 && n!=1) || (n>1 && m!=1)) PsychErrorExitMsg(PsychError_user, "The sixth argument must be either a row- or columnvector of filterModes.");
    }
    numFilterModes = m * n;

    // Default to bilinear filtering:
    filterMode = (numFilterModes == 1) ? filterModes[0] : 1;

    // Optional globalAlphas:
    m=n=p=0;
    if (PsychAllocInDoubleMatArg(7, kPsychArgOptional, &m, &n, &p, &globalAlphas)) {
        if ((p!=1) || (m>1 && n!=1) || (n>1 && m!=1)) PsychErrorExitMsg(PsychError_user, "The seventh argument must be either a row- or columnvector of globalAlpha values.");
    }
    numAlphas = m * n;
    globalAlpha = (numAlphas == 1) ? globalAlphas[0] : 1.0;

    // Optional auxParameters:
    auxParameters = NULL;
    m=n=p=0;
    if (PsychAllocInDoubleMatArg(11, kPsychArgOptional, &m, &n, &p, &auxParameters)) {
        if ((p!=1) || (m < 4) || ((m % 4) !=0)|| (n < 1)) PsychErrorExitMsg(PsychError_user, "The 11th argument must be a column vector or matrix of 'auxParameter' values with at least 4 components and component count a multiple of four.");
    }
    numAuxParams = n;
    numAuxComponents = m;

    // Check for consistency: Each parameter must be either not present, present once,
    // or present as many times as all other multi-parameters:
    numRef = (numsrcRects > numdstRects) ? numsrcRects : numdstRects;
    numRef = (numRef > numTexs) ? numRef : numTexs;
    numRef = (numRef > nc) ? numRef : nc;
    numRef = (numRef > numAlphas) ? numRef : numAlphas;
    numRef = (numRef > numFilterModes) ? numRef : numFilterModes;
    numRef = (numRef > numAngles) ? numRef : numAngles;
    numRef = (numRef > numAuxParams) ? numRef : numAuxParams;

    if (numTexs > 1 && numTexs != numRef) {
        printf("PTB-ERROR: Number of provided texture handles %i doesn't match number of other primitives %i!\n", numTexs, numRef);
        PsychErrorExitMsg(PsychError_user, "Inconsistent number of arguments provided to Screen('DrawTextures').");
    }

    if (numsrcRects > 1 && numsrcRects != numRef) {
        printf("PTB-ERROR: Number of provided source rectangles %i doesn't match number of other primitives %i!\n", numsrcRects, numRef);
        PsychErrorExitMsg(PsychError_user, "Inconsistent number of arguments provided to Screen('DrawTextures').");
    }

    if (numdstRects > 1 && numdstRects != numRef) {
        printf("PTB-ERROR: Number of provided destination rectangles %i doesn't match number of other primitives %i!\n", numdstRects, numRef);
        PsychErrorExitMsg(PsychError_user, "Inconsistent number of arguments provided to Screen('DrawTextures').");
    }

    if (numAngles > 1 && numAngles != numRef) {
        printf("PTB-ERROR: Number of provided rotation angles %i doesn't match number of other primitives %i!\n", numAngles, numRef);
        PsychErrorExitMsg(PsychError_user, "Inconsistent number of arguments provided to Screen('DrawTextures').");
    }

    if (numAlphas > 1 && numAlphas != numRef) {
        printf("PTB-ERROR: Number of provided global alpha values %i doesn't match number of other primitives %i!\n", numAlphas, numRef);
        PsychErrorExitMsg(PsychError_user, "Inconsistent number of arguments provided to Screen('DrawTextures').");
    }

    if (numFilterModes > 1 && numFilterModes != numRef) {
        printf("PTB-ERROR: Number of provided filtermode values %i doesn't match number of other primitives %i!\n", numFilterModes, numRef);
        PsychErrorExitMsg(PsychError_user, "Inconsistent number of arguments provided to Screen('DrawTextures').");
    }

    if (nc > 1 && nc != numRef) {
        printf("PTB-ERROR: Number of provided modulateColors %i doesn't match number of other primitives %i!\n", nc, numRef);
        PsychErrorExitMsg(PsychError_user, "Inconsistent number of arguments provided to Screen('DrawTextures').");
    }

    if (numAuxParams > 1 && numAuxParams != numRef) {
        printf("PTB-ERROR: Number of provided 'auxParameter' column vectors %i doesn't match number of other primitives %i!\n", numAuxParams, numRef);
        PsychErrorExitMsg(PsychError_user, "Inconsistent number of arguments provided to Screen('DrawTextures').");
    }

    // Assign optional override texture shader, if any provided:
    textureShader = -1;
    PsychCopyInIntegerArg(9, kPsychArgOptional, &textureShader);

    // Assign any other optional special flags:
    PsychCopyInIntegerArg(10, kPsychArgOptional, &specialFlags);

    // Check if efficient batch drawing is possible at the GL level:
    if (isclassic && (numTexs == 1) && (numFilterModes <= 1)) {
        batchIt = TRUE;
    }
    else {
        batchIt = FALSE;
    }

    if (PsychPrefStateGet_Verbosity() > 5)
        printf("PTB-DEBUG: DrawTextures optimized batch submit: %i\n", (int) batchIt);

    if (batchIt) {
        // Signal start of new batch with numRef drawn textures, all sourced from source and
        // drawn into windo target with filterMode:
        PsychBatchBlitTexturesToDisplay(0, numRef, source, target, NULL, NULL, 0, (int) filterMode, 1.0);
    }

    // Texture blitting loop:
    for (i=0; i < numRef; i++) {
        // Draw i'th texture:

        // Check if more than one texture provided. If not then the one single texture has been
        // setup already above:
        if (numTexs > 1) {
            // More than one texture handle provided: Need to allocate i'th one in:
            if(!IsWindowIndex((PsychWindowIndexType) texids[i])) {
                printf("PTB-ERROR: %i th entry in texture handle vector is not a valid handle!\n", i + 1);
                PsychErrorExitMsg(PsychError_user, "Invalid texture handle provided to Screen('DrawTextures').");
            }

            // Get it:
            FindWindowRecord((PsychWindowIndexType) texids[i], &source);
            if(source->windowType!=kPsychTexture) {
                printf("PTB-ERROR: %i th entry in texture handle vector is not a valid handle!\n", i + 1);
                PsychErrorExitMsg(PsychError_user, "The second argument supplied was not a texture handle!");
            }
        }

        // Source rectangle provided?
        if (numsrcRects > 1) {
            // Get i'th source rectangle:
            PsychCopyRect(sourceRect, &(srcRects[i*4]));
        } else if (numsrcRects == 1) {
            // Single source rect provided - get it:
            PsychCopyRect(sourceRect, &(srcRects[0]));
        } else {
            // No source rect provided: Take rectangle of current texture as srcRect:
            PsychCopyRect(sourceRect,source->clientrect);
        }

        // Skip this texture if sourceRect is an empty rect:
        if (IsPsychRectEmpty(sourceRect)) continue;

        // Destination rectangle provided?
        if (numdstRects > 1) {
            // Get i'th destination rectangle:
            PsychCopyRect(targetRect, &(dstRects[i*4]));
        } else if (numdstRects == 1) {
            // Single destination rect provided - get it:
            PsychCopyRect(targetRect, &(dstRects[0]));
        } else {
            // No destination rect provided: Center the current sourceRect in the current
            // target window and use that as destination:
            PsychCopyRect(tempRect, target->clientrect);
            PsychCenterRectInRect(sourceRect, tempRect, targetRect);
        }

        // Skip this texture if targetRect is an empty rect:
        if (IsPsychRectEmpty(targetRect)) continue;

        if (numAngles > 1) rotationAngle = rotationAngles[i];
        if (numFilterModes > 1) filterMode = filterModes[i];
        if (numAlphas > 1) globalAlpha = globalAlphas[i];

        // Disable alpha if modulateColor active:
        if (nc > 0) globalAlpha = DBL_MAX;

        // Pass auxParameters for current primitive in the auxShaderParams field.
        target->auxShaderParamsCount = numAuxComponents;
        if (numAuxParams > 0) {
            if (numAuxParams == 1) {
                target->auxShaderParams = auxParameters;
            }
            else {
                target->auxShaderParams = &(auxParameters[i * numAuxComponents]);
            }
        }
        else {
            target->auxShaderParams = NULL;
        }

        // Multiple modulateColors provided?
        if (nc > 1) {
            if (isclassic) {
                // Yes. Set it up as current vertex color: We submit to internal currentColor for
                // shader based color processing and via glColorXXX() for fixed pipe processing:
                if (mc==3) {
                    if (colors) {
                        // RGB double:
                        target->currentColor[0]=colors[i*3 + 0];
                        target->currentColor[1]=colors[i*3 + 1];
                        target->currentColor[2]=colors[i*3 + 2];
                        target->currentColor[3]=1.0;
                    }
                    else {
                        // RGB uint8:
                        target->currentColor[0]=((double) bytecolors[i*3 + 0] / 255.0);
                        target->currentColor[1]=((double) bytecolors[i*3 + 1] / 255.0);
                        target->currentColor[2]=((double) bytecolors[i*3 + 2] / 255.0);
                        target->currentColor[3]=1.0;
                    }
                }
                else {
                    if (colors) {
                        // RGBA double:
                        target->currentColor[0]=colors[i*4 + 0];
                        target->currentColor[1]=colors[i*4 + 1];
                        target->currentColor[2]=colors[i*4 + 2];
                        target->currentColor[3]=colors[i*4 + 3];
                    }
                    else {
                        // RGBA uint8:
                        target->currentColor[0]=((double) bytecolors[i*4 + 0] / 255.0);
                        target->currentColor[1]=((double) bytecolors[i*4 + 1] / 255.0);
                        target->currentColor[2]=((double) bytecolors[i*4 + 2] / 255.0);
                        target->currentColor[3]=((double) bytecolors[i*4 + 3] / 255.0);
                    }
                }
            }
            else {
                PsychSetArrayColor(target, i, mc, colors, bytecolors);
            }
        }

        // Ok, everything assigned. Check parameters:
        if (filterMode > 5) {
            PsychErrorExitMsg(PsychError_user, "filterMode needs to be negative for a specific blur level, or at most 5 for other modes.");
        }

        // Set rotation mode flag for texture matrix rotation if secialFlags is set accordingly:
        if (specialFlags & kPsychUseTextureMatrixForRotation) source->specialflags|=kPsychUseTextureMatrixForRotation;
        if (specialFlags & kPsychDontDoRotation) source->specialflags|=kPsychDontDoRotation;

        if (textureShader > -1) {
            backupShader = source->textureFilterShader;
            source->textureFilterShader = -1 * textureShader;
        }

        if (batchIt) {
            // Add current element to the batch to be drawn:
            PsychBatchBlitTexturesToDisplay(2, numRef, source, target, sourceRect, targetRect, rotationAngle, (int) filterMode, globalAlpha);
        }
        else {
            // Perform blit operation for i'th texture, either with or without an override texture shader applied:
            glColor4dv(target->currentColor);
            PsychBlitTextureToDisplay(source, target, sourceRect, targetRect, rotationAngle, (int) filterMode, globalAlpha);
        }

        if (textureShader > -1)
            source->textureFilterShader = backupShader;

        // Reset rotation mode flag:
        source->specialflags &= ~(kPsychUseTextureMatrixForRotation | kPsychDontDoRotation);

        // Next one...
    }

    target->auxShaderParams = NULL;
    target->auxShaderParamsCount = 0;

    if (batchIt) {
        // Finalize batch drawing:
        PsychBatchBlitTexturesToDisplay(1, numRef, source, target, NULL, NULL, 0, (int) filterMode, 1.0);
    }

    // Mark end of drawing op. This is needed for single buffered drawing:
    PsychFlushGL(target);

    return(PsychError_none);
}
