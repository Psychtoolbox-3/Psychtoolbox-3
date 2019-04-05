/*
    PsychToolbox3/Source/Common/Screen/PsychGLGlue.c

    PLATFORMS: All

    AUTHORS:
    Allen Ingling   awi     Allen.Ingling@nyu.edu
    Mario Kleiner   mk      mario.kleiner.de@gmail.com

    HISTORY:
    09/09/02        awi     wrote it.

    DESCRIPTION:

    Functions to convert between Psych toolbox data types and GL data types.
*/

#include "Screen.h"

// Cached global glApiType value for cases where no windowRecord is available:
static int global_glApiType = 0;

/*
    PsychConvertColorToDoubleVector()

    Accept a color structure and a screen depth and return either three or four double values in the range between
    0-1 which specify r, g, b and optinally alpha values.

    The value array argument should be be four elements long.

*/
int PsychConvertColorToDoubleVector(PsychColorType *color, PsychWindowRecordType *windowRecord, GLdouble *valueArray)
{
    GLdouble deno;

    // Read denominator from windowRecord. Need to get rid of the sign, because it
    // encodes if we have color clamping enabled or not:
    deno = fabs(windowRecord->colorRange);

    switch(color->mode){
        case kPsychIndexColor:
            valueArray[0]=color->value.index.i/deno;
            return(1);
        case kPsychRGBColor:
            valueArray[0]=color->value.rgb.r/deno;
            valueArray[1]=color->value.rgb.g/deno;
            valueArray[2]=color->value.rgb.b/deno;
            valueArray[3]=1.0;	// Needed to init the valueArray completely.
            return(3);
        case kPsychRGBAColor:
            valueArray[0]=color->value.rgba.r/deno;
            valueArray[1]=color->value.rgba.g/deno;
            valueArray[2]=color->value.rgba.b/deno;
            valueArray[3]=(color->value.rgba.a == DBL_MAX) ? 1.0 : color->value.rgba.a/deno;
            return(4);
        case kPsychUnknownColor:
            PsychErrorExitMsg(PsychError_internal,"Unspecified display mode");
    }
    PsychErrorExitMsg(PsychError_internal,"Unknown display mode");
    return(0); //make the compiler happy.
}

// Define submission command for submitting single unclamped colors to drawshader.
// For now, we use the first (primary) 4D texture coordinate, as this is a predefined
// attribute:
#define HDRglColor4dv(v) glTexCoord4dv((v))

/*
 *   PsychSetGLColor()
 *
 *   Accept a Psych color structure and a depth value and call the appropriate variant of glColor.
 */
void PsychSetGLColor(PsychColorType *color, PsychWindowRecordType *windowRecord)
{
    int numVals;

    numVals=PsychConvertColorToDoubleVector(color, windowRecord, (GLdouble*) &(windowRecord->currentColor));
    if (numVals < 3 || numVals > 4) PsychErrorExitMsg(PsychError_internal, "Palette mode not yet implemented or illegal color specifier.");

    // Set the color in GL:
    if (windowRecord->defaultDrawShader) {
        // Drawshader color submission:
        HDRglColor4dv(windowRecord->currentColor);
    }
    else {
        // Fixed function pipe:
        if (PsychIsGLClassic(windowRecord)) {
            // OpenGL-1/2:
            glColor4dv(windowRecord->currentColor);
        }
        else {
            // OpenGL-ES 1.x: glColor4f() is only available function.
            glColor4f((float) windowRecord->currentColor[0], (float) windowRecord->currentColor[1], (float) windowRecord->currentColor[2], (float) windowRecord->currentColor[3]);
            PsychGLColor4f(windowRecord, (float) windowRecord->currentColor[0], (float) windowRecord->currentColor[1], (float) windowRecord->currentColor[2], (float) windowRecord->currentColor[3]);
        }
    }
}

/* PsychSetupVertexColorArrays()
 * Helper routine, called from the different batch drawing functions of Screen():
 */
void PsychSetupVertexColorArrays(PsychWindowRecordType *windowRecord, psych_bool enable, int mc, double* colors, unsigned char *bytecolors)
{
    if (enable) {
        // Enable and setup whatever's used:
        if (windowRecord->defaultDrawShader) {
            // Shader based unclamped path:
            if (colors) glTexCoordPointer(mc, PSYCHGLFLOAT, 0, colors);

            // Can't support uint8 datatype for this vertex attribute :-(
            if (bytecolors) PsychErrorExitMsg(PsychError_user, "Sorry, this function can't accept matrices of uint8 type for colors\nif color clamping is disabled or high precision mode active.\n Use the double() operator to convert to double matrix.");

            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
            glColorPointer(4, PSYCHGLFLOAT, 0, NULL);
        }
        else {
            // Standard path:
            if (colors)     glColorPointer(mc, PSYCHGLFLOAT, 0, colors);
            if (bytecolors) glColorPointer(mc, GL_UNSIGNED_BYTE, 0, bytecolors);

            glEnableClientState(GL_COLOR_ARRAY);
            glTexCoordPointer(4, PSYCHGLFLOAT, 0, NULL);
        }
    }
    else {
        // Disable whatever's used:
        if (windowRecord->defaultDrawShader) {
            // Shader based unclamped path:
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        }
        else {
            // Standard path:
            glDisableClientState(GL_COLOR_ARRAY);
        }

        glColorPointer(4, PSYCHGLFLOAT, 0, NULL);
        glTexCoordPointer(4, PSYCHGLFLOAT, 0, NULL);
    }
}

/* PsychSetArrayColor()
 * Helper routine, called from the different batch drawing functions of Screen():
 */
static GLdouble currentColor[4];
void PsychSetArrayColor(PsychWindowRecordType *windowRecord, int i, int mc, double* colors, unsigned char *bytecolors)
{
    psych_bool isgles = !PsychIsGLClassic(windowRecord);

    if ((windowRecord->defaultDrawShader) || isgles) {
        // Draw shader assigned or OpenGL-ES in use. Need to feed color values into high-precision
        // alternative channel for unclamped, high-precision color handling:
        if (mc==3) {
            i = i * 3;
            if (colors) {
                // RGB double:
                currentColor[0]=colors[i++];
                currentColor[1]=colors[i++];
                currentColor[2]=colors[i++];
                currentColor[3]=1.0;
            }
            else {
                // RGB uint8:
                currentColor[0]=((double) bytecolors[i++] / 255.0);
                currentColor[1]=((double) bytecolors[i++] / 255.0);
                currentColor[2]=((double) bytecolors[i++] / 255.0);
                currentColor[3]=1.0;
            }
        }
        else {
            i = i * 4;
            if (colors) {
                // RGBA double:
                currentColor[0]=colors[i++];
                currentColor[1]=colors[i++];
                currentColor[2]=colors[i++];
                currentColor[3]=colors[i++];
            }
            else {
                // RGBA uint8:
                currentColor[0]=((double) bytecolors[i++] / 255.0);
                currentColor[1]=((double) bytecolors[i++] / 255.0);
                currentColor[2]=((double) bytecolors[i++] / 255.0);
                currentColor[3]=((double) bytecolors[i++] / 255.0);
            }
        }

        if (isgles) {
            // GLES can only do glColor4f(), nothing else:
            glColor4f((float) currentColor[0], (float) currentColor[1], (float) currentColor[2], (float) currentColor[3]);
        }
        else {
            HDRglColor4dv(currentColor);
        }
    }
    else {
        // Standard fixed-function pipeline assigned: Feed into standard glColorXXX() calls:
        if (mc==3) {
            if (colors) {
                // RGB double:
                glColor3dv(&(colors[i*3]));
            }
            else {
                // RGB uint8:
                glColor3ubv(&(bytecolors[i*3]));
            }
        }
        else {
            if (colors) {
                // RGBA double:
                glColor4dv(&(colors[i*4]));
            }
            else {
                // RGBA uint8:
                glColor4ubv(&(bytecolors[i*4]));
            }
        }
    }

    return;
}

/* PsychGLClear()
 *
 * Helper around glClearColor() and glClear() - takes special issues
 * caused by HDR framebuffer support into account.
 */
void PsychGLClear(PsychWindowRecordType *windowRecord)
{
    int oldShader, nowShader;

    // Unclamped/High-precision color mode enabled via GLSL shaders?
    if ((windowRecord->defaultDrawShader != 0) && (windowRecord->defaultDrawShader == windowRecord->unclampedDrawShader)) {
        // Yes. Can't use standard clear, but need to clear via drawing a full-window rect with
        // clear color:

        // Query currently bound shader:
        oldShader = PsychGetCurrentShader(windowRecord);

        // Assign hdr draw shader:
        nowShader = PsychSetShader(windowRecord, -1);

        // Assign HDR clear color for window:
        HDRglColor4dv(windowRecord->clearColor);

        // Draw a fullscreen rect in the clear color, make sure
        // no alpha blending is active:
        if (glIsEnabled(GL_BLEND)) {
            glDisable(GL_BLEND);
            PsychGLRect(windowRecord->rect);
            glEnable(GL_BLEND);
        }
        else {
            PsychGLRect(windowRecord->rect);
        }

        // Revert to old shader binding:
        if (nowShader != oldShader) PsychSetShader(windowRecord, oldShader);
    }
    else {
        // Standard clear path: Can use OpenGL's fast color buffer clear:
        glClearColor((GLclampf) windowRecord->clearColor[0], (GLclampf) windowRecord->clearColor[1], (GLclampf) windowRecord->clearColor[2], (GLclampf) windowRecord->clearColor[3]);
        glClear(GL_COLOR_BUFFER_BIT);
    }

    return;
}

/*
 * PsychGLRect()
 */
void PsychGLRect(PsychRectType psychRect)
{
    // Call our glRectd helper which makes sure glRectd() works
    // not only on classic desktop GL, but also on latest desktop
    // OpenGL and OpenGL-ES. Passing a NULL windowRecord will make
    // it use the cached api type spec:
    PsychGLRectd(NULL,
                (GLdouble)(psychRect[kPsychLeft]),
                (GLdouble)(psychRect[kPsychTop]),
                (GLdouble)(psychRect[kPsychRight]),
                (GLdouble)(psychRect[kPsychBottom]));
}

char *PsychGetGLErrorNameString(GLenum errorConstant)
{
    static char GL_NO_ERROR_str[] = "GL_NO_ERROR";
    static char GL_INVALID_ENUM_str[] = "GL_INVALID_ENUM";
    static char GL_INVALID_VALUE_str[] = "GL_INVALID_VALUE";
    static char GL_INVALID_OPERATION_str[] = "GL_INVALID_OPERATION";
    static char GL_STACK_OVERFLOW_str[] = "GL_STACK_OVERFLOW";
    static char GL_STACK_UNDERFLOW_str[] = "GL_STACK_UNDERFLOW";
    static char GL_OUT_OF_MEMORY_str[] = "GL_OUT_OF_MEMORY";
    static char GL_TABLE_TOO_LARGE_str[] = "GL_TABLE_TOO_LARGE";
    static char GL_INVALID_FRAMEBUFFER_OPERATION_EXT_str[] = "GL_INVALID_FRAMEBUFFER_OPERATION_EXT";
    static char unrecognized_error_str[] = "unrecognized GL error constant";

    switch(errorConstant){
        case GL_NO_ERROR: return(GL_NO_ERROR_str);
        case GL_INVALID_ENUM: return(GL_INVALID_ENUM_str);
        case GL_INVALID_VALUE: return(GL_INVALID_VALUE_str);
        case GL_INVALID_OPERATION: return(GL_INVALID_OPERATION_str);
        case GL_STACK_OVERFLOW: return(GL_STACK_OVERFLOW_str);
        case GL_STACK_UNDERFLOW: return(GL_STACK_UNDERFLOW_str);
        case GL_OUT_OF_MEMORY: return(GL_OUT_OF_MEMORY_str);
        case GL_TABLE_TOO_LARGE: return(GL_TABLE_TOO_LARGE_str);
        case GL_INVALID_FRAMEBUFFER_OPERATION_EXT: return(GL_INVALID_FRAMEBUFFER_OPERATION_EXT_str);
    }
    return(unrecognized_error_str);
}

/*
 *  PsychGetGLErrorListString()
 */
psych_bool PsychGetGLErrorListString(const char **errorListStr)
{
    #define MAX_GL_ERROR_LIST_LEN               2048
    #define MAX_GL_ERROR_LIST_DELTA_LEN         256
    static char errorListString[MAX_GL_ERROR_LIST_LEN];
    char        *errorNameStr;
    int         currentIndex, deltaStrLen, nextCurrentIndex;
    GLenum      glError;
    psych_bool  isError=FALSE;

    // Skip this routine with "no error" return status, if kPsychAvoidCPUGPUSync
    // is set as conserveVRAM setting by usercode:
    *errorListStr = NULL;
    if (PsychPrefStateGet_ConserveVRAM() & kPsychAvoidCPUGPUSync) return(FALSE);

    currentIndex=0;
    for(glError=glGetError(); glError!=GL_NO_ERROR; glError=glGetError()){
        errorNameStr=PsychGetGLErrorNameString(glError);
        deltaStrLen = (int) strlen(errorNameStr)+2;  //2 chars: comma and space
        nextCurrentIndex=currentIndex+deltaStrLen;
        if (nextCurrentIndex >= MAX_GL_ERROR_LIST_LEN)
            PsychErrorExitMsg(PsychError_internal,"string memory overflow");
        if (isError)
            sprintf(&(errorListString[currentIndex]), " ,%s", errorNameStr);
        else
            sprintf(&(errorListString[currentIndex]), "%s", errorNameStr);
        currentIndex=nextCurrentIndex;
        isError=TRUE;
    }

    if(isError)
        *errorListStr=errorListString;
    else
        *errorListStr=NULL;
    return(isError);
}

void PsychTestForGLErrorsC(int lineNum, const char *funcName, const char *fileName)
{
    psych_bool  isError;
    const char  *glErrorListString;

    isError=PsychGetGLErrorListString(&glErrorListString);
    if (isError)
        PsychErrorExitC(PsychError_OpenGL,
                        glErrorListString,
                        lineNum,
                        funcName,
                        fileName);
}

/*
        PsychExtractQuadVertexFromRect()

        Return one of the four vertices define by a Psych rect in a 2-element array of GLdoubles.
        Vertices are numbered from the top left corner (0) clockwise to the bottom left corner (3).
*/
GLdouble *PsychExtractQuadVertexFromRect(double *rect, int vertexNumber, GLdouble *vertex)
{
    switch(vertexNumber){
        case 0:
            vertex[0]=(GLdouble)rect[0];
            vertex[1]=(GLdouble)rect[1];
            break;
        case 1:
            vertex[0]=(GLdouble)rect[2];
            vertex[1]=(GLdouble)rect[1];
            break;
        case 2:
            vertex[0]=(GLdouble)rect[2];
            vertex[1]=(GLdouble)rect[3];
            break;
        case 3:
            vertex[0]=(GLdouble)rect[0];
            vertex[1]=(GLdouble)rect[3];
            break;
        default:
            PsychErrorExitMsg(PsychError_internal, "Illegal vertex value");
    }
    return(vertex);
}

/* PsychPrepareRenderBatch()
 *
 * Perform setup for a batch of render requests for a specific primitive. Some 2D Screen
 * drawing commands allow to specify a list of primitives to draw instead of only a single
 * one. E.g. 'DrawDots' allows to draw thousands of dots with one single DrawDots command.
 * This helper routine is called by such batch-capable commands. It checks which input arguments
 * are provided and if it's a single one or multiple ones. It sets up the rendering pipe accordingly,
 * performing required conversion steps. The actual drawing routine just needs to perform primitive
 * specific code.
 */
void PsychPrepareRenderBatch(PsychWindowRecordType *windowRecord, int coords_pos, int* coords_count, double** xy, int colors_pos, int* colors_count, int* colorcomponent_count, double** colors, unsigned char** bytecolors, int sizes_pos, int* sizes_count, double** size, psych_bool usefloat)
{
    PsychColorType      color;
    int                 m,n,p,mc,nc,pc;
    int                 i, nrpoints, nrsize;
    psych_bool          isArgThere, isdoublecolors, isuint8colors, usecolorvector, needxy;
    double              *tmpcolors, *pcolors, *tcolors;
    double              convfactor, whiteValue;
    float               *pcolorsf, *tcolorsf;
    float               convfactorf;
    unsigned char       *pcolorsb, *tcolorsb;

    needxy = (coords_pos > 0) ? TRUE: FALSE;
    coords_pos = abs(coords_pos);
    colors_pos = abs(colors_pos);
    sizes_pos = abs(sizes_pos);

    // Get mandatory or optional xy coordinates argument
    isArgThere = PsychIsArgPresent(PsychArgIn, coords_pos);
    if(!isArgThere && needxy) {
        PsychErrorExitMsg(PsychError_user, "No position argument supplied");
    }

    if (isArgThere) {
        if (usefloat) {
            PsychAllocInFloatMatArg(coords_pos, TRUE, &m, &n, &p, (float**) xy);
        }
        else {
            PsychAllocInDoubleMatArg(coords_pos, TRUE, &m, &n, &p, xy);
        }

        if (p!=1 || (m!=*coords_count && (m*n)!=*coords_count)) {
            printf("PTB-ERROR: Coordinates must be a %i tuple or a %i rows vector.\n", *coords_count, *coords_count);
            PsychErrorExitMsg(PsychError_user, "Invalid format for coordinate specification.");
        }

        if (m!=1) {
            nrpoints=n;
            *coords_count = n;
        }
        else {
            // Special case: 1 row vector provided for single argument.
            nrpoints=1;
            *coords_count = 1;
        }
    }
    else {
        nrpoints = 0;
        *coords_count = 0;
    }

    if (size) {
        // Get optional size argument
        isArgThere = PsychIsArgPresent(PsychArgIn, sizes_pos);
        if (!isArgThere) {
                // No size provided: Use a default size of 1.0:
                *size = (double *) PsychMallocTemp(sizeof(double));
                if (usefloat)
                    **((float**) size) = 1;
                else
                    *size[0] = 1;
                nrsize=1;
        } else {
            if (usefloat) {
                PsychAllocInFloatMatArg(sizes_pos, TRUE, &m, &n, &p, (float**) size);
            }
            else {
                PsychAllocInDoubleMatArg(sizes_pos, TRUE, &m, &n, &p, size);
            }

            if (p!=1) PsychErrorExitMsg(PsychError_user, "Size must be a scalar or a vector with one column or row");
            nrsize=m*n;
            if (nrsize!=nrpoints && nrsize!=1 && *sizes_count!=1) PsychErrorExitMsg(PsychError_user, "Size vector must contain one size value per item.");
        }

        *sizes_count = nrsize;
    }

    // Check if color argument is provided:
    isArgThere = PsychIsArgPresent(PsychArgIn, colors_pos);
    if(!isArgThere) {
        // No color argument provided - Use defaults:
        whiteValue=PsychGetWhiteValueFromWindow(windowRecord);
        PsychLoadColorStruct(&color, kPsychIndexColor, whiteValue ); //index mode will coerce to any other.
        usecolorvector=false;
    }
    else {
        // Some color argument provided. Check first, if it's a valid color vector:
        if (usefloat) {
            isdoublecolors = PsychAllocInFloatMatArg(colors_pos, kPsychArgAnything, &mc, &nc, &pc, (float**) colors);
        }
        else {
            isdoublecolors = PsychAllocInDoubleMatArg(colors_pos, kPsychArgAnything, &mc, &nc, &pc, colors);
        }

        isuint8colors  = PsychAllocInUnsignedByteMatArg(colors_pos, kPsychArgAnything, &mc, &nc, &pc, bytecolors);

        // Do we have a color vector, aka one element per vertex?
        if ((isdoublecolors || isuint8colors) && pc==1 && mc!=1 && nc==nrpoints && nrpoints>1) {
            // Looks like we might have a color vector... ... Double-check it:
            if (mc!=3 && mc!=4) PsychErrorExitMsg(PsychError_user, "Color vector must be a 3 or 4 row vector");

            // Yes. colors is a valid pointer to it.
            usecolorvector=true;

            if (isdoublecolors) {
                if (fabs(windowRecord->colorRange)!=1) {
                    // We have to loop through the vector and divide all values by windowRecord->colorRange, so the input values
                    // 0-colorRange get mapped to the range 0.0-1.0, as OpenGL expects values in range 0-1 when
                    // a color vector is passed in Double- or Float format.
                    // This is inefficient, as it burns some cpu-cycles, but necessary to keep color
                    // specifications consistent in the PTB - API.
                    if (usefloat) {
                        // OpenGL-ES 1 code path: Color arrays must have 4 component RGBA spec in
                        // single precision float format. Make it so:
                        convfactorf = (float) (1.0 / fabs(windowRecord->colorRange));
                        tmpcolors = PsychMallocTemp(sizeof(float) * nc * 4);
                        pcolorsf = (float*) *colors;
                        tcolorsf = (float*) tmpcolors;
                        if (mc == 4) {
                            // 4 channel RGBA input: Normalize.
                            for (i=0; i<(nc*mc); i++) {
                                *(tcolorsf++) = *(pcolorsf++) * convfactorf;
                            }
                        }
                        else {
                            // 3 channel RGB input: Normalize RGB, add a 1.0 alpha channel:
                            for (i=0; i < nc; i++) {
                                *(tcolorsf++) = *(pcolorsf++) * convfactorf;
                                *(tcolorsf++) = *(pcolorsf++) * convfactorf;
                                *(tcolorsf++) = *(pcolorsf++) * convfactorf;
                                *(tcolorsf++) = 1.0;
                            }
                        }
                    }
                    else {
                        // Desktop OpenGL:
                        convfactor = 1.0 / fabs(windowRecord->colorRange);
                        tmpcolors=PsychMallocTemp(sizeof(double) * nc * mc);
                        pcolors = *colors;
                        tcolors = tmpcolors;
                        for (i=0; i<(nc*mc); i++) {
                            *(tcolors++)=(*pcolors++) * convfactor;
                        }
                    }
                }
                else {
                    // colorRange is == 1 --> No remapping needed as colors are already in proper range.
                    if (usefloat && (mc == 3)) {
                        // OpenGL-ES 1 and only 3 channel RGB input: Extend to RGBA:
                        tmpcolors = PsychMallocTemp(sizeof(float) * nc * 4);
                        pcolorsf = (float*) *colors;
                        tcolorsf = (float*) tmpcolors;

                        // 3 channel RGB input: Add a 1.0 alpha channel:
                        for (i = 0; i < nc; i++) {
                            *(tcolorsf++) = *(pcolorsf++);
                            *(tcolorsf++) = *(pcolorsf++);
                            *(tcolorsf++) = *(pcolorsf++);
                            *(tcolorsf++) = 1.0;
                        }
                    }
                    else {
                        // Desktop OpenGL: Just setup pointer to our unaltered input color vector:
                        tmpcolors=*colors;
                    }
                }

                *colors = tmpcolors;
            }
            else {
                // Color vector in uint8 format. Nothing to do, unless this is OpenGL-ES 1 and
                // input is only RGB instead of required RGBA:
                if (usefloat && (mc == 3)) {
                    // OpenGL-ES 1 and only 3 channel RGB input: Extend to RGBA:
                    tmpcolors = PsychMallocTemp(nc * 4);
                    pcolorsb = *bytecolors;
                    tcolorsb = (unsigned char*) tmpcolors;

                    // 3 channel RGB input: Add a 1.0 aka 255 alpha channel:
                    for (i = 0; i < nc; i++) {
                        *(tcolorsb++) = *(pcolorsb++);
                        *(tcolorsb++) = *(pcolorsb++);
                        *(tcolorsb++) = *(pcolorsb++);
                        *(tcolorsb++) = 255;
                    }

                    *bytecolors = (unsigned char*) tmpcolors;
                }
            }
        }
        else {
            // No color vector provided: Check for a single valid color triplet or quadruple:
            usecolorvector=false;
            isArgThere=PsychCopyInColorArg(colors_pos, TRUE, &color);
        }
    }

    // OpenGL-ES 1 always has 4 component RGBA color vectors:
    if (usefloat) mc = 4;

    // Enable this windowRecords framebuffer as current drawingtarget:
    PsychSetDrawingTarget(windowRecord);

    // Setup default drawshader:
    PsychSetShader(windowRecord, -1);

    // Setup alpha blending properly:
    PsychUpdateAlphaBlendingFactorLazily(windowRecord);

    // Setup common color for all objects if no color vector has been provided:
    if (!usecolorvector) {
        PsychCoerceColorMode(&color);
        PsychSetGLColor(&color, windowRecord);
        *colors_count = 1;
    }
    else {
        *colors_count = nc;
    }
    *colorcomponent_count = mc;

    return;
}

/* Emit a single pixel in top-left corner of window and wait for its rendering
* to complete. Our classic trick to wait for double-buffer swap completion on
* systems where we don't have better system-provided timestamping and syncing
* methods. This needs different implementations on classic OpenGL vs. non-
* immediate mode OpenGL.
*
* A 'flushOnly' flag of TRUE will only flush, not wait.
*
*/
void PsychWaitPixelSyncToken(PsychWindowRecordType *windowRecord, psych_bool flushOnly)
{
    // Classic desktop OpenGL in use?
    if (PsychIsGLClassic(windowRecord)) {
        // Yes. Use our classic fixed-function immediate mode method:
        glBegin(GL_POINTS);
        glColor4f(0, 0, 0, 0);
        glVertex2i(10, 10);
        glEnd();
    }
    else {
        // No. Avoid immediate mode functions, they won't work:
        GLfloat glverts[2] = { 10, 10 };
        glVertexPointer(2, GL_FLOAT, 0, glverts);
        glEnableClientState(GL_VERTEX_ARRAY);
        glDrawArrays(GL_POINTS, 0, 1);
        glDisableClientState(GL_VERTEX_ARRAY);
    }

    // flushOnly flag - Don't wait for write to happen, just flush it:
    if (!flushOnly) {
        // Wait for write completion - Used for sync and timestamping:
        glFinish();
    }
    else {
        // Only flush:
        glFlush();
    }

    if (flushOnly && (PsychPrefStateGet_Verbosity() > 15)) printf("PTB-DEBUG: PixelSyncToken write + glFlush().\n");
}

GLenum PsychGLFloatType(PsychWindowRecordType *windowRecord)
{
    // On OpenGL-ES we only have GL_FLOAT data type, not GL_DOUBLE, so we need
    // to convert all double input data into float data for the calling routine:
    return((PsychIsGLES(windowRecord)) ? GL_FLOAT : GL_DOUBLE);
}

psych_bool PsychIsGLES(PsychWindowRecordType *windowRecord)
{
    // No windowRecord (NULL) given? Use cached global api type:
    if (!windowRecord) return((global_glApiType > 0) ? TRUE : FALSE);

    // windowRecord given. Cache its api type as global api type:
    global_glApiType = windowRecord->glApiType;

    // And return decision based on its api type:
    return((windowRecord->glApiType > 0) ? TRUE : FALSE);
}

psych_bool PsychIsGLClassic(PsychWindowRecordType *windowRecord)
{
    // No windowRecord (NULL) given? Use cached global api type:
    if (!windowRecord) return((global_glApiType == 0) ? TRUE : FALSE);

    // windowRecord given. Cache its api type as global api type:
    global_glApiType = windowRecord->glApiType;

    // And return decision based on its api type:
    return((windowRecord->glApiType == 0) ? TRUE : FALSE);
}

#define PSYCH_MAX_IMMSIM_VERTEX_COMPONENTS 4*3*16384
static unsigned int gl_buffer_index = 0;
static GLenum gl_buffer_primitivetype;
static float gl_buffer[PSYCH_MAX_IMMSIM_VERTEX_COMPONENTS];
static float currentTexCoord[4] = { 0, 0, 0, 1 };

void PsychGLBegin(PsychWindowRecordType *windowRecord, GLenum primitive)
{
    // Classic OpenGL? Use glBegin() and be done:
    if (PsychIsGLClassic(windowRecord)) {
        glBegin(primitive);
        return;
    }

    // OpenGL 3+ or OpenGL-ES. Need to emulate immediate mode stuff:
    if (gl_buffer_index > 0) PsychErrorExitMsg(PsychError_internal, "PsychGLBegin() called while already inside PsychGLBegin() block!");

    // Assign primitive for this Begin->End block:
    gl_buffer_primitivetype = primitive;

    return;
}

void PsychGLEnd(PsychWindowRecordType *windowRecord)
{
    // Classic OpenGL? Use glEnd() and be done:
    if (PsychIsGLClassic(windowRecord)) {
        glEnd();
        return;
    }

    // OpenGL 3+ or OpenGL-ES. Need to emulate immediate mode stuff:

    // Work to do?
    if (gl_buffer_index > 0) {
        // Yes: Submit draw call:
        glVertexPointer(4, GL_FLOAT, 4 * 3 * sizeof(float), &gl_buffer[0 * 4]);
        glColorPointer(4,  GL_FLOAT, 4 * 3 * sizeof(float), &gl_buffer[1 * 4]);
        glTexCoordPointer(4,  GL_FLOAT, 4 * 3 * sizeof(float), &gl_buffer[2 * 4]);

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);

        glDrawArrays(gl_buffer_primitivetype, 0, gl_buffer_index / (4 * 3));

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    }

    // Reset primitive for this Begin->End block:
    gl_buffer_primitivetype = (GLenum) 0xffff;
    gl_buffer_index = 0;

    return;
}

void PsychGLVertex4f(PsychWindowRecordType *windowRecord, float x, float y, float z, float w)
{
    // Classic OpenGL? Use glEnd() and be done:
    if (PsychIsGLClassic(windowRecord)) {
        glVertex4f(x,y,z,w);
        return;
    }

    // Make sure we don't overflow our internal static vertex array bufer:
    if (gl_buffer_index > PSYCH_MAX_IMMSIM_VERTEX_COMPONENTS - 12) {
        printf("PTB-ERROR: Overflow of internal immediate mode rendering buffer! Limit of %i components reached.\n", PSYCH_MAX_IMMSIM_VERTEX_COMPONENTS);
        printf("PTB-ERROR: You must reduce the workload, or your stimulus image will be damaged. Likely culprit is Screen('Fill/FrameOval') with too high level of detail!\n\n");
        return;
    }

    // Store (x,y,z,w) vertex pos first:
    gl_buffer[gl_buffer_index++] = x;
    gl_buffer[gl_buffer_index++] = y;
    gl_buffer[gl_buffer_index++] = z;
    gl_buffer[gl_buffer_index++] = w;

    // Then current color (r,g,b,a):
    gl_buffer[gl_buffer_index++] = (float) currentColor[0];
    gl_buffer[gl_buffer_index++] = (float) currentColor[1];
    gl_buffer[gl_buffer_index++] = (float) currentColor[2];
    gl_buffer[gl_buffer_index++] = (float) currentColor[3];

    // Then current texture coordinates (s,t,u,v):
    gl_buffer[gl_buffer_index++] = currentTexCoord[0];
    gl_buffer[gl_buffer_index++] = currentTexCoord[1];
    gl_buffer[gl_buffer_index++] = currentTexCoord[2];
    gl_buffer[gl_buffer_index++] = currentTexCoord[3];

    return;
}

void PsychGLColor4f(PsychWindowRecordType *windowRecord, float r, float g, float b, float a)
{
    double colors[4];
    colors[0] = r;
    colors[1] = g;
    colors[2] = b;
    colors[3] = a;
    PsychSetArrayColor(windowRecord, 0, 4, colors, NULL);
}

void PsychGLTexCoord4f(PsychWindowRecordType *windowRecord, float s, float t, float u, float v)
{
    // Classic OpenGL? Use glEnd() and be done:
    if (PsychIsGLClassic(windowRecord)) {
        glTexCoord4f(s, t, u, v);
        return;
    }

    currentTexCoord[0] = s;
    currentTexCoord[1] = t;
    currentTexCoord[2] = u;
    currentTexCoord[3] = v;
}

void PsychGLRectd(PsychWindowRecordType *windowRecord, double x1, double y1, double x2, double y2)
{
    // Classic OpenGL? Use glRectd() and be done:
    if (PsychIsGLClassic(windowRecord)) {
        glRectd(x1, y1, x2, y2);
        return;
    }

    // Want to emulate a GL_QUAD draw: We do it as a triangle strip:
    GLBEGIN(GL_TRIANGLE_STRIP);
    GLVERTEX2d(x1, y1);
    GLVERTEX2d(x2, y1);
    GLVERTEX2d(x1, y2);
    GLVERTEX2d(x2, y2);
    GLEND();
}

/* A rather lame implementation of gluDisc() - OpenGL-ES compatible, but as inefficient as
* it gets. Good enough for a start though...
*
* Thinking about this was avoided by adapting freely available sample code from Thomas Visser:
* http://www.cocos2d-iphone.org/forum/topic/2207
*
*/
void PsychDrawDisc(PsychWindowRecordType *windowRecord, float xc, float yc, float innerRadius, float outerRadius, int numSlices, float xScale, float yScale, float startAngle, float arcAngle)
{
    float Rads, outerX, outerY, innerX, innerY;
        int i;
        #ifndef M_PI
        #define M_PI 3.141592654
        #endif

    // Sweep clock-wise over arcAngle degrees, split up into 'numSlices' steps:
    float step = (float) ((arcAngle * 2 * M_PI / 360) / (float) numSlices);

    // Start at startAngle degrees, where 0 deg. = Upward, 90 deg. = Rightward.
    startAngle = (float) (M_PI/2 - startAngle * 2 * M_PI / 360);

    GLBEGIN(GL_TRIANGLE_STRIP);
        for (i = 0; i <= numSlices; i++) {
            // calculating the current vertice on the outer side of the segment
            Rads = startAngle - i * step;
            outerX = xScale * outerRadius * (float) cos( Rads );
            outerY = yScale * outerRadius * (float) sin( Rads );
            GLVERTEX2f(xc + outerX, yc + outerY);

            // calculating the current vertice on the inner side of the segment
            innerX = xScale * innerRadius * (float) cos( Rads );
            innerY = yScale * innerRadius * (float) sin( Rads );
            GLVERTEX2f(xc + innerX, yc + innerY);
        }
    GLEND();

    return;
}
