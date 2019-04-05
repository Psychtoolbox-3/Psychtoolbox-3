/*
        SCREENDrawDots.c

        AUTHORS:

            kas@princeton.edu               kas     Keith Schneider
            fcalabro@bu.edu                 fjc     Finnegan Calabro
            mario.kleiner.de@gmail.com      mk      Mario Kleiner

        PLATFORMS:

            All.

        HISTORY:

            mm/dd/yy

            12/13/04        kas     Created
            12/17/04        fjc     added antialiased points
            12/17/04        kas     added dot_type flag
            1/11/04         awi     Merged into pychtoolbox.org distribution.
            1/01/05         mk      Performance optimizations and small bug-fixes
                                    - Fixes some minor bugs (or possible bugs, they would only
                                    trigger in some usage cases): GL-Calls were missing to undo
                                    parameter changes like glPointSize,
                                    glEnable(GL_POINT_SMOOTH). This can affect later
                                    invocations of DrawDots, e.g., Call DrawDots with point
                                    smoothing on --> Smoothing gets enabled, after that it
                                    would not be possible to disable point smoothing again -->
                                    any following DrawDots call would always render smooth
                                    points, regardless of setting.

                                    - Uses an optimized way of drawing the dots: Instead of
                                    looping over the 2-column vector of dots in a for-loop and
                                    submitting each Point via glVertex2d(), it now *requires* a
                                    2-row vector of dots as input and then renders all dots
                                    with one function call to glDrawArrays(). This is the most
                                    efficient way of passing a huge number of primitives
                                    (Points, rectangles, ...) to the OpenGL-System, avoiding
                                    for-loop and function call overhead and allowing possible
                                    optimizations by the graphics driver.
            1/12/05         awi     Merged mk's improvements into psychtoolbox.org distribution.
            1/30/05         mk      Bug-fixes: Adding PsychMallocTemp()'s and a few checks to avoid memory corruption
                                    crashes.
            2/25/05         awi     Added call to PsychUpdateAlphaBlendingFactorLazily().  Drawing now obeys settings by Screen('BlendFunction').
            3/22/05         mk      Added possibility to spec vectors with individual color and size spec per dot.
            4/29/05         mk      Bugfix for color vectors: They should also take values in range 0-255 instead of 0.0-1.0.
            11/14/06        mk      We now also accept color vectors in uint8 format and pass them directly for higher efficiency.
*/

#include "Screen.h"

static char PointSmoothFragmentShaderSrc[] =
"\n"
"uniform int drawRoundDots;\n"
"varying vec4 unclampedFragColor;\n"
"varying float pointSize;\n"
"\n"
"void main()\n"
"{\n"
"    /* Non-round, aliased square dots requested? */\n"
"    if (drawRoundDots == 0) {\n"
"       /* Yes. Simply passthrough unclamped color and be done: */\n"
"       gl_FragColor = unclampedFragColor;\n"
"       return;\n"
"    }\n"
"\n"
"    /* Passthrough RGB color values: */\n"
"    gl_FragColor.rgb = unclampedFragColor.rgb;\n"
"\n"
"    /* Adapt alpha value dependent on relative radius of the fragment within a dot:   */\n"
"    /* This for point smoothing on GPU's that don't support this themselves.          */\n"
"    /* Points on the border of the dot, at [radius - 0.5 ; radius + 0.5] pixels, will */\n"
"    /* get their alpha value reduced from 1.0 * alpha to 0.0, so they completely      */\n"
"    /* disappear over a distance of 1 pixel distance unit. The - 1.0 subtraction      */\n"
"    /* in clamp() accounts for the fact that pointSize is 2.0 pixels larger than user */\n"
"    /* code requested. This 1 pixel padding around true size avoids cutoff artifacts. */\n"
"    float r = length(gl_TexCoord[1].st - vec2(0.5, 0.5)) * pointSize;\n"
"    r = 1.0 - clamp(r - (0.5 * pointSize - 1.0 - 0.5), 0.0, 1.0);\n"
"    gl_FragColor.a = unclampedFragColor.a * r;\n"
"    if (r <= 0.0)\n"
"        discard;\n"
"}\n\0";

char PointSmoothVertexShaderSrc[] =
"/* Vertex shader: Emulates fixed function pipeline, but in HDR color mode passes    */ \n"
"/* gl_MultiTexCoord0 as varying unclampedFragColor to circumvent vertex color       */ \n"
"/* clamping on gfx-hardware / OS combos that don't support unclamped operation:     */ \n"
"/* PTBs color handling is expected to pass the vertex color in gl_MultiTexCoord0    */ \n"
"/* for unclamped drawing for this reason in unclamped color mode. gl_MultiTexCoord2 */ \n"
"/* delivers individual point size (diameter) information for each point.            */ \n"
"\n"
"uniform int useUnclampedFragColor;\n"
"varying float pointSize;\n"
"varying vec4 unclampedFragColor;\n"
"\n"
"void main()\n"
"{\n"
"    if (useUnclampedFragColor > 0) {\n"
"       /* Simply copy input unclamped RGBA pixel color into output varying color: */\n"
"       unclampedFragColor = gl_MultiTexCoord0;\n"
"    }\n"
"    else {\n"
"       /* Simply copy regular RGBA pixel color into output varying color: */\n"
"       unclampedFragColor = gl_Color;\n"
"    }\n"
"\n"
"    /* Output position is the same as fixed function pipeline: */\n"
"    gl_Position = ftransform();\n"
"\n"
"    /* Point size comes via texture coordinate set 2: Make diameter 2 pixels bigger    */\n"
"    /* than requested, to have some 1 pixel security margin around the dot, to avoid   */\n"
"    /* cutoff artifacts for the rendered point-sprite quad. Compensate in frag-shader. */\n"
"    pointSize = gl_MultiTexCoord2[0] + 2.0;\n"
"    gl_PointSize = pointSize;\n"
"}\n\0";

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "[minSmoothPointSize, maxSmoothPointSize, minAliasedPointSize, maxAliasedPointSize] = Screen('DrawDots', windowPtr, xy [,size] [,color] [,center] [,dot_type][, lenient]);";
//                          1                   2                   3                    4                                         1          2    3       4        5         6           7
static char synopsisString[] =
"Quickly draw an array of dots.  "
"\"xy\" is a two-row vector containing the x and y coordinates of the dot centers, "
"relative to \"center\" (default center is [0 0]).\n"
"\"size\" is the diameter of each dot in pixels (default is 1). "
"Instead of a common size for all dots you can also provide a "
"vector which defines a different dot size for each dot. Different graphics cards do "
"have different limits on the maximum size of dots, 10 or 63 are typical limits.\n"
"\"color\" is the the clut index (scalar or [r g b a] vector) "
"that you want to poke into each dot pixel (default is black).  "
"Instead of a single \"color\" you can also provide a 3 or 4 row vector,"
"which specifies an individual RGB or RGBA color for each corresponding point.\n"
"\"dot_type\" is a flag that determines what type of dot is drawn: "
"0 (default) and 4 draw square dots, whereas 1, 2 and 3 draw round dots (circles) with anti-aliasing: "
"1 favors performance, 2 tries to use high-quality anti-aliasing, if supported by your hardware. "
"3 Uses a builtin shader-based implementation. "
"dot_type 1 and 2 may not be supported by all graphics cards and drivers. On some systems "
"Screen() will then automatically switch to dot_type 3 - our own implementation - in such a case. "
"If you use round dot_type 1, 2 or 3 you'll also need to set a proper blending mode with the "
"Screen('BlendFunction') command, e.g., GL_SRC_ALPHA + GL_ONE_MINUS_SRC_ALPHA. A dot_type of 4 will "
"draw square dots like dot_type 0, but may be faster when drawing lots of dots of different sizes by "
"use of an efficient shader based path.\n"
"\"lenient\" If set to 1, will not check the sizes of dots for validity, so you can try requesting "
"sizes bigger than what the hardware claims to support.\n\n"
"The optional return arguments [minSmoothPointSize, maxSmoothPointSize, minAliasedPointSize, maxAliasedPointSize] "
"allow you to query the minimum and maximum allowed 'size' for smooth anti-aliased dots (dot_type 1,2,3) and for "
"non anti-aliased square dots (dot_type 0 and 4). Calling [...] = Screen('DrawDots', windowPtr) will only query "
"these point size limits without drawing any dots.\n";

static char seeAlsoString[] = "BlendFunction";

PsychError SCREENDrawDots(void)
{
    PsychWindowRecordType                   *windowRecord, *parentWindowRecord;
    int                                     m,n,p,mc,nc,idot_type;
    int                                     i, nrpoints, nrsize;
    psych_bool                              isArgThere, usecolorvector;
    double                                  *xy, *size, *center, *dot_type, *colors;
    float                                   *sizef;
    unsigned char                           *bytecolors;
    GLfloat                                 pointsizerange[2];
    psych_bool                              lenient = FALSE;
    psych_bool                              usePointSizeArray = FALSE;
    static psych_bool                       nocando = FALSE;
    int                                     oldverbosity;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    // Check for superfluous arguments
    PsychErrorExit(PsychCapNumInputArgs(7));   //The maximum number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(4));  //The maximum number of outputs

    // Get the window record from the window record argument and get info from the window record
    PsychAllocInWindowRecordArg(1, kPsychArgRequired, &windowRecord);

    // Get dot_type argument, if any, as it is already needed for a pure point size range query below:
    isArgThere = PsychIsArgPresent(PsychArgIn, 6);
    if(!isArgThere){
        idot_type = 0;
    } else {
        PsychAllocInDoubleMatArg(6, TRUE, &m, &n, &p, &dot_type);
        if(p != 1 || n != 1 || m != 1 || (dot_type[0] < 0 || dot_type[0] > 4))
            PsychErrorExitMsg(PsychError_user, "dot_type must be 0, 1, 2, 3 or 4");
        idot_type = (int) dot_type[0];
    }

    // Query for supported point size range?
    if (PsychGetNumOutputArgs() > 0) {
        PsychSetDrawingTarget(windowRecord);

        // Always query and return aliased range:
        glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, (GLfloat*) &pointsizerange);
        PsychCopyOutDoubleArg(3, FALSE, (double) pointsizerange[0]);
        PsychCopyOutDoubleArg(4, FALSE, (double) pointsizerange[1]);

        // If driver supports smooth points and usercode doesn't specify a dot type (idot_type 0)
        // or does not request shader + point-sprite based drawing then return smooth point
        // size range as "smooth point size range" - query and assign it. Otherwise, ie., code
        // explicitely wants to use a shader (idot_type >= 3) or has to use one, we will use
        // point-sprites and that means the GL_ALIASED_POINT_SIZE_RANGE limits apply also to
        // our shader based smooth dots, so return those:
        if ((windowRecord->gfxcaps & kPsychGfxCapSmoothPrimitives) && (idot_type < 3))
            glGetFloatv(GL_POINT_SIZE_RANGE, (GLfloat*) &pointsizerange);

        // Whatever the final choice for smooth dots is, return its limits:
        PsychCopyOutDoubleArg(1, FALSE, (double) pointsizerange[0]);
        PsychCopyOutDoubleArg(2, FALSE, (double) pointsizerange[1]);

        // If this was only a query then we are done:
        if (PsychGetNumInputArgs() < 2)
            return(PsychError_none);
    }

    // Query, allocate and copy in all vectors...
    nrpoints = 2;
    nrsize = 0;
    colors = NULL;
    bytecolors = NULL;

    PsychPrepareRenderBatch(windowRecord, 2, &nrpoints, &xy, 4, &nc, &mc, &colors, &bytecolors, 3, &nrsize, &size, (GL_FLOAT == PsychGLFloatType(windowRecord)));
    usecolorvector = (nc>1) ? TRUE:FALSE;

    // Assign sizef as float-type array of sizes, if float mode active, NULL otherwise:
    sizef = (GL_FLOAT == PsychGLFloatType(windowRecord)) ? (float*) size : NULL;

    // Get center argument
    isArgThere = PsychIsArgPresent(PsychArgIn, 5);
    if(!isArgThere){
        center = (double *) PsychMallocTemp(2 * sizeof(double));
        center[0] = 0;
        center[1] = 0;
    } else {
        PsychAllocInDoubleMatArg(5, TRUE, &m, &n, &p, &center);
        if(p!=1 || n!=2 || m!=1) PsychErrorExitMsg(PsychError_user, "center must be a 1-by-2 vector");
    }

    // Turn on antialiasing to draw circles? Or idot_type 4 for shader based square dots?
    if (idot_type) {
        // Smooth point rendering supported by gfx-driver and hardware? And user does not request our own stuff?
        if ((idot_type == 3) || (idot_type == 4) || !(windowRecord->gfxcaps & kPsychGfxCapSmoothPrimitives)) {
            // No. Need to roll our own shader + point sprite solution.
            if (!windowRecord->smoothPointShader && !nocando) {
                parentWindowRecord = PsychGetParentWindow(windowRecord);
                if (!parentWindowRecord->smoothPointShader) {
                    // Build and assign shader to parent window, but allow this to silently fail:
                    oldverbosity = PsychPrefStateGet_Verbosity();
                    PsychPrefStateSet_Verbosity(0);
                    parentWindowRecord->smoothPointShader = PsychCreateGLSLProgram(PointSmoothFragmentShaderSrc, PointSmoothVertexShaderSrc, NULL);
                    PsychPrefStateSet_Verbosity(oldverbosity);
                }

                if (parentWindowRecord->smoothPointShader) {
                    // Got one compiled - assign it for use:
                    windowRecord->smoothPointShader = parentWindowRecord->smoothPointShader;
                }
                else {
                    // Failed. Record this failure so we can avoid retrying at next DrawDots invocation:
                    nocando = TRUE;
                }
            }

            if (windowRecord->smoothPointShader) {
                // Activate point smooth shader, and point sprite operation on texunit 1 for coordinates on set 1:
                PsychSetShader(windowRecord, windowRecord->smoothPointShader);
                glActiveTexture(GL_TEXTURE1);
                glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
                glActiveTexture(GL_TEXTURE0);
                glEnable(GL_POINT_SPRITE);

                // Tell shader from where to get its color information: Unclamped high precision colors from texture coordinate set 0, or regular colors from vertex color attribute?
                glUniform1i(glGetUniformLocation(windowRecord->smoothPointShader, "useUnclampedFragColor"), (windowRecord->defaultDrawShader) ? 1 : 0);

                // Tell shader if it should shade smooth round dots, or square dots:
                glUniform1i(glGetUniformLocation(windowRecord->smoothPointShader, "drawRoundDots"), (idot_type != 4) ? 1 : 0);

                // Tell shader about current point size in pointSize uniform:
                glEnable(GL_PROGRAM_POINT_SIZE);
                usePointSizeArray = TRUE;
            }
            else if (idot_type != 4) {
                // Game over for round dot drawing:
                PsychErrorExitMsg(PsychError_user, "Point smoothing unsupported on your system and our shader based implementation failed as well in Screen('DrawDots').");
            }
            else {
                // Type 4 requested but unsupported. Fallback to type 0, which is the same, just slower:
                idot_type = 0;
            }

            // Request square dots, without anti-aliasing: Better compatibility with
            // shader + point sprite operation, and needed for idot_type 0 fallback.
            glDisable(GL_POINT_SMOOTH);
            glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, (GLfloat*) &pointsizerange);
        }
        else {
            // User wants hw anti-aliased round smooth dots (idot_type = 1 or 2) and
            // hardware + driver support this. Request smooth points from hardware:
            glEnable(GL_POINT_SMOOTH);
            glGetFloatv(GL_POINT_SIZE_RANGE, (GLfloat*) &pointsizerange);

            // A dot type of 2 requests highest quality point smoothing:
            glHint(GL_POINT_SMOOTH_HINT, (idot_type > 1) ? GL_NICEST : GL_DONT_CARE);
        }
    }
    else {
        glDisable(GL_POINT_SMOOTH);
        glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, (GLfloat*) &pointsizerange);
    }

    // Does ES-GPU only support a fixed point diameter of 1 pixel?
    if ((pointsizerange[1] <= 1) && PsychIsGLES(windowRecord)) {
        // Yes. Not much point bailing on this, as it should be easily visible
        // during testing of a studies code on a OpenGL-ES device.
        lenient = TRUE;
    }

    // Accept optional 'lenient' flag, if provided:
    PsychCopyInFlagArg(7, FALSE, &lenient);

    // Set size of a single dot:
    if (!lenient && ((sizef && (sizef[0] > pointsizerange[1] || sizef[0] < pointsizerange[0])) ||
        (!sizef && (size[0] > pointsizerange[1] || size[0] < pointsizerange[0])))) {
        printf("PTB-ERROR: You requested a point size of %f units, which is not in the range (%f to %f) supported by your graphics hardware.\n",
                (sizef) ? sizef[0] : size[0], pointsizerange[0], pointsizerange[1]);
        PsychErrorExitMsg(PsychError_user, "Unsupported point size requested in Screen('DrawDots').");
    }

    // Setup initial common point size for all points:
    if (!usePointSizeArray) glPointSize((sizef) ? sizef[0] : (float) size[0]);
    if (usePointSizeArray) glMultiTexCoord1f(GL_TEXTURE2, (sizef) ? sizef[0] : (float) size[0]);

    // Setup modelview matrix to perform translation by 'center':
    glMatrixMode(GL_MODELVIEW);

    // Make a backup copy of the matrix:
    glPushMatrix();

    // Apply a global translation of (center(x,y)) pixels to all following points:
    glTranslatef((float) center[0], (float) center[1], 0);

    // Render the array of 2D-Points - Efficient version:
    // This command sequence allows fast processing of whole arrays
    // of vertices (or points, in this case). It saves the call overhead
    // associated with the original implementation below and is potentially
    // optimized in specific OpenGL implementations.

    // Pass a pointer to the start of the point-coordinate array:
    glVertexPointer(2, PSYCHGLFLOAT, 0, &xy[0]);

    // Enable fast rendering of arrays:
    glEnableClientState(GL_VERTEX_ARRAY);

    if (usecolorvector) {
        PsychSetupVertexColorArrays(windowRecord, TRUE, mc, colors, bytecolors);
    }

    // Render all n points, starting at point 0, render them as POINTS:
    if ((nrsize == 1) || usePointSizeArray) {
        // Only one common size provided, or efficient shader based
        // path in use. We can use the fast path of only submitting
        // one glDrawArrays call to draw all GL_POINTS. For a single
        // common size, no further setup is needed.
        if (nrsize > 1) {
            // Individual size for each dot provided. Setup texture unit 2
            // with a 1D texcoord array that stores per point size info in
            // texture coordinate set 2. But first validate point sizes:
            for (i = 0; i < nrpoints; i++) {
                if (!lenient && ((sizef && (sizef[i] > pointsizerange[1] || sizef[i] < pointsizerange[0])) ||
                    (!sizef && (size[i] > pointsizerange[1] || size[i] < pointsizerange[0])))) {
                    printf("PTB-ERROR: You requested a point size of %f units, which is not in the range (%f to %f) supported by your graphics hardware.\n",
                           (sizef) ? sizef[i] : size[i], pointsizerange[0], pointsizerange[1]);
                    PsychErrorExitMsg(PsychError_user, "Unsupported point size requested in Screen('DrawDots').");
                }
            }

            // Do we need the GL_FLOAT data glTexCoordPointer(1, ...) workaround?
            // See explanation in PsychWindowSupport.c: PsychDetectAndAssignGfxCapabilities():
            if (windowRecord->specialflags & kPsychNeedVBODouble12Workaround) {
                // Yes: Convert double size buffer to a float sizef buffer and then
                // submit that as GL_FLOAT:
                sizef = (float*) PsychMallocTemp(nrpoints * sizeof(float));
                for (i = 0; i < nrpoints; i++)
                    sizef[i] = (float) size[i];
            }

            // Sizes are fine, setup texunit 2:
            glClientActiveTexture(GL_TEXTURE2);
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
            glTexCoordPointer(1, (sizef) ? GL_FLOAT : GL_DOUBLE, 0, (sizef) ? (const GLvoid*) sizef : (const GLvoid*) size);
        }

        // Draw all points:
        glDrawArrays(GL_POINTS, 0, nrpoints);

        if (nrsize > 1) {
            // Individual size for each dot provided. Reset texture unit 2:
            glTexCoordPointer(1, (sizef) ? GL_FLOAT : GL_DOUBLE, 0, (const GLvoid*) NULL);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);

            // Back to default texunit 0:
            glClientActiveTexture(GL_TEXTURE0);
        }
    }
    else {
        // Different size for each dot provided and we can't use our shader based implementation:
        // We have to do One GL - call per dot:
        for (i = 0; i < nrpoints; i++) {
            if (!lenient && ((sizef && (sizef[i] > pointsizerange[1] || sizef[i] < pointsizerange[0])) ||
                (!sizef && (size[i] > pointsizerange[1] || size[i] < pointsizerange[0])))) {
                printf("PTB-ERROR: You requested a point size of %f units, which is not in the range (%f to %f) supported by your graphics hardware.\n",
                        (sizef) ? sizef[i] : size[i], pointsizerange[0], pointsizerange[1]);
                PsychErrorExitMsg(PsychError_user, "Unsupported point size requested in Screen('DrawDots').");
            }

            // Setup point size for this point:
            glPointSize((sizef) ? sizef[i] : (float) size[i]);

            // Render point:
            glDrawArrays(GL_POINTS, i, 1);
        }
    }

    // Disable fast rendering of arrays:
    glDisableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, PSYCHGLFLOAT, 0, NULL);

    if (usecolorvector) PsychSetupVertexColorArrays(windowRecord, FALSE, 0, NULL, NULL);

    // Restore old matrix from backup copy, undoing the global translation:
    glPopMatrix();

    // turn off antialiasing again
    if (idot_type) {
        glDisable(GL_POINT_SMOOTH);

        if (windowRecord->smoothPointShader) {
            // Deactivate point smooth shader and point sprite operation on texunit 1:
            PsychSetShader(windowRecord, 0);
            glActiveTexture(GL_TEXTURE1);
            glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_FALSE);
            glActiveTexture(GL_TEXTURE0);
            glDisable(GL_POINT_SPRITE);
            glDisable(GL_PROGRAM_POINT_SIZE);
        }
    }

    // Reset pointsize to 1.0
    glPointSize(1);

    // Mark end of drawing op. This is needed for single buffered drawing:
    PsychFlushGL(windowRecord);

    //All psychfunctions require this.
    return(PsychError_none);
}
