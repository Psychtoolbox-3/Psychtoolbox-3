/*
 Common/Screen/SCREENSetOpenGLTextureFromMemPointer.c
 
 AUTHORS:
 mario.kleiner at tuebingen.mpg.de  mk
 
 PLATFORMS:	All
 
 
 HISTORY:
 06/23/06  mk		Created.
 
 DESCRIPTION:
 Creates or updates a Psychtoolbox texture from/with the pixel data that is
 passed as a memory pointer. External C-code can create image data and store
 it in a standard C memory array or buffer. Then it can use this function to
 pass a memory pointer (void*) to PTB. The pointer has to point to the start
 of the memory buffer. PTB will create or update a texture from the data.
 */


#include "Screen.h"

static char useString[]="[textureHandle rect] = Screen('SetOpenGLTextureFromMemPointer', windowPtr, textureHandle, imagePtr, width, height, depth [, upsidedown][, target][, glinternalformat][, gltype][, extdataformat][, specialFlags]);";
static char synopsisString[] = 
"DANGEROUS! C-PROGRAMMING EXPERTS ONLY! OTHERS STAY AWAY! "
"Convert raw image data, provided as a pointer to a system memory buffer, into a Psychtoolbox texture. "
"CAUTION: Providing wrong arguments to this subfunction will definitely crash Psychtoolbox and Matlab! "
"\"windowPtr\" is the handle of the onscreen window to which the texture should be attached. "
"'textureHandle' is either the Psychtoolbox handle for an existing PTB texture that should be updated "
"with the new raw pixel data, or the special value zero or [] if a completely new PTB texture should "
"be created from the raw pixel data. 'imagePtr' is a C programming language (void*) memory pointer that "
"points to the start of a memory buffer with the raw image data. The pointer needs to be encoded in a "
"special way. Read the source code of Memorybuffer2Texture.c and Memorybuffer2TextureDemo.m for minimal "
"examples on how to do this properly. Memorybuffer2Texture.cc demonstrates the same for GNU/Octave. "
"The optional parameter 'target' is the type of texture object to create, e.g., GL_TEXTURE_2D. Normally "
"you'll just leave it out, so PTB can choose the optimal texture format for a system. 'width' and 'height' are the "
"width and height of the image. 'depth' is the depth of the image: 1 = Input is a greyscale image with "
"1 Byte per pixel. 3 = Input is a RGB truecolor image with 3 Bytes per pixel. 2 and 4 are like "
"1 and 3, but with one byte for an alpha (transparency) channel added. The optional flag 'upsidedown' can "
"be set to one (default is zero) to create textures upside-down - inverted in vertical direction. "
"'glinternalformat', 'gltype' and 'extdataformat' are all optional. If you provide one of them, you'll need "
"to provide all of them. In that case, PTB will create an OpenGL texture with exactly the specified internal "
"format, assuming input data that is of type 'gltype' and in numeric data format 'extdataformat'. You are "
"completely responsible for passing properly formatted data. This is mostly useful for injecting high dynamic "
"range texture images and other exotic texture formats.\n"
"'specialFlags' Special optional texture flags, see help for Screen('MakeTexture').\n"
"The function returns the textureHandle of the PTB texture and its defining rectangle. "
"This routine allows external C-Code (Mex- or Oct-Files) to inject raw image data as a texture into PTB. ";

static char seeAlsoString[] = "SetOpenGLTexture MakeTexture ";

PsychError SCREENSetOpenGLTextureFromMemPointer(void) 
{
    PsychWindowRecordType *windowRecord, *textureRecord;
    int w, h, d, testarg, upsidedown, glinternalformat, glexternaltype, glexternalformat, usefloatformat;
    int specialFlags = 0;
    double doubleMemPtr;
    GLenum target = 0;
    w=h=d=-1;
    doubleMemPtr = 0;
    upsidedown = 0;

    //all subfunctions should have these two lines.  
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    
    PsychErrorExit(PsychCapNumInputArgs(12));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(5)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(2));    // The maximum number of outputs
    
    // Get the window record from the window record argument and get info from the window record
    PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
    
    // Get the texture record from the texture record argument.
    // Check if either none ( [] or '' ) or the special value zero was
    // provided as Psychtoolbox textureHandle. In that case, we create a new
    // empty texture record instead of reusing an existing one.
    testarg=0;
    PsychCopyInIntegerArg(2, FALSE, &testarg);
    if (testarg==0) {
	// No valid textureHandle provided. Create a new empty textureRecord.
        PsychCreateWindowRecord(&textureRecord);
        textureRecord->windowType=kPsychTexture;
        textureRecord->screenNumber = windowRecord->screenNumber;

	// Assign parent window and copy its inheritable properties:
	PsychAssignParentWindow(textureRecord, windowRecord);

        // Mark it valid and return handle to userspace:
        PsychSetWindowRecordValid(textureRecord);
    }
    else {
        // None of the special values provided. We assume its a handle to a valid
        // and existing PTB texture and try to retrieve the textureRecord:
        PsychAllocInWindowRecordArg(2, TRUE, &textureRecord);
    }
    
    // Is it  a textureRecord?
    if (!PsychIsTexture(textureRecord)) {
        PsychErrorExitMsg(PsychError_user, "You tried to set texture information on something else than a texture!");
    }
    
    // Query double-encoded memory pointer:
    PsychCopyInDoubleArg(3, TRUE, &doubleMemPtr);
    
    // Query width:
    PsychCopyInIntegerArg(4, TRUE, &w);

    // Query height:
    PsychCopyInIntegerArg(5, TRUE, &h);

    // Query depth:
    PsychCopyInIntegerArg(6, TRUE, &d);

    // Query (optional) upsidedown - flag:
    PsychCopyInIntegerArg(7, FALSE, &upsidedown);
    
    // Query (optional) OpenGL texture target:
    PsychCopyInIntegerArg(8, FALSE, (int*) &target);
 
    // Query (optional) full format spec:
    glinternalformat = 0;
    PsychCopyInIntegerArg(9, FALSE, &glinternalformat);
    if (glinternalformat>0) {
      // Ok copy the (now non-optional) remaining format spec:
      PsychCopyInIntegerArg(10, TRUE, &glexternalformat);
      PsychCopyInIntegerArg(11, TRUE, &glexternaltype);      
    }

    // Query optional specialFlags:
    PsychCopyInIntegerArg(12, FALSE, &specialFlags);
    
    // Safety checks:
    if (doubleMemPtr == 0) {
        PsychErrorExitMsg(PsychError_user, "You tried to set invalid (NULL) imagePtr.");
    }

    if (w<=0) {
        PsychErrorExitMsg(PsychError_user, "You tried to set invalid (negative) texture width.");
    }

    if (h<=0) {
        PsychErrorExitMsg(PsychError_user, "You tried to set invalid (negative) texture height.");
    }
    
    if (d<=0) {
        PsychErrorExitMsg(PsychError_user, "You tried to set invalid (negative) texture depth.");
    }
    
    if (d>4) {
        PsychErrorExitMsg(PsychError_user, "You tried to set invalid (greater than four) texture depth.");
    }

    if (target!=0 && target!=GL_TEXTURE_RECTANGLE_EXT && target!=GL_TEXTURE_2D) {
        PsychErrorExitMsg(PsychError_user, "You tried to set invalid texture target.");
    }

    // Activate OpenGL rendering context of windowRecord:
    PsychSetGLContext(windowRecord);

    PsychTestForGLErrors();

    // Ok, setup texture record for texture:
    textureRecord->depth = d * 8;
    textureRecord->nrchannels = d;
    PsychMakeRect(textureRecord->rect, 0, 0, w, h);

    // Client rect of a texture is always == rect of it:
    PsychCopyRect(textureRecord->clientrect, textureRecord->rect);

    // Override texture target, if one was provided:
    if (target!=0) textureRecord->texturetarget = target;

    // Orientation is normally set to 2 - like an upright Offscreen window texture.
    // If upsidedown flag is set, then we do 3 - an upside down Offscreen window texture.
    textureRecord->textureOrientation = (upsidedown>0) ? 3 : 2;
    
    if (glinternalformat!=0) {
      textureRecord->textureinternalformat = glinternalformat;
      textureRecord->textureexternalformat = glexternalformat;
      textureRecord->textureexternaltype = glexternaltype;
    }

    // Setting memsize to zero prevents unwanted free() operation in PsychDeleteTexture...
    textureRecord->textureMemorySizeBytes = 0;

    // This will retrieve an OpenGL compatible pointer to the raw pixel data and assign it to our texmemptr:
    textureRecord->textureMemory = (GLuint*) PsychDoubleToPtr(doubleMemPtr);
    // printf("InTexPtr %p , %.20e", PsychDoubleToPtr(doubleMemPtr), doubleMemPtr);

    // Let PsychCreateTexture() do the rest of the job of creating, setting up and
    // filling an OpenGL texture with memory buffers image content:
    PsychCreateTexture(textureRecord);

    // Assign GLSL filter-/lookup-shaders if needed: usefloatformat is queried.
    // The 'userRequest' flag is set to zero for now.
    target = PsychGetTextureTarget(textureRecord);
    glBindTexture(target, textureRecord->textureNumber);
    glGetTexLevelParameteriv(target, 0, GL_TEXTURE_RED_SIZE, (GLint*) &d);
    if (d <= 0) glGetTexLevelParameteriv(target, 0, GL_TEXTURE_LUMINANCE_SIZE, (GLint*) &d);
    glBindTexture(target, 0);

    usefloatformat = 0;
    if (d == 16) usefloatformat = 1;
    if (d >= 32) usefloatformat = 2;

    PsychAssignHighPrecisionTextureShaders(textureRecord, windowRecord, usefloatformat, (specialFlags & 2) ? 1 : 0);

    // specialFlags setting 8? Disable auto-mipmap generation:
    if (specialFlags & 0x8) textureRecord->specialflags |= kPsychDontAutoGenMipMaps;    

    // A specialFlags setting of 32? Protect texture against deletion via Screen('Close') without providing a explicit handle:
    if (specialFlags & 32) textureRecord->specialflags |= kPsychDontDeleteOnClose;    

    // Return new (or old) PTB textureHandle for this texture:
    PsychCopyOutDoubleArg(1, FALSE, textureRecord->windowIndex);
    PsychCopyOutRectArg(2, FALSE, textureRecord->rect);

    // Done.
    return(PsychError_none);
}
