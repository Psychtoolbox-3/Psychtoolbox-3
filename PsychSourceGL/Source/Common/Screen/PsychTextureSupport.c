/*
 *    PsychTextureSupport.c
 *
 *    PLATFORMS:
 *
 *        This is the OS independent version.
 *
 *    AUTHORS:
 *
 *        Allen Ingling         awi     Allen.Ingling@nyu.edu
 *        Mario Kleiner         mk      mario.kleiner.de@gmail.com
 *
 *    HISTORY:
 *
 *        3/9/04        awi     Wrote it
 *        7/27/04       awi     Removed commented-out code.
 *        1/4/05        mk      Performance optimizations, some bug-fixes.
 *        1/13/05       awi     Merged in Mario's changes of 1/4/05 into the psychtoolbox.org master.
 *        1/23/05       awi     Merged mk's update which adds glEnable(GL_TEXTURE_RECTANGLE_EXT) call.  Fixes DrawText and DrawTexture
 *                              interference reported by Frans Cornelissen and others;
 *        1/30/05       mk      Small bugfix: Now specifying sourceRect in DrawTexture really works as expected.
 *        5/13/05       mk      Support for rotated drawing of textures. Switched to bilinear texture filtering instead of nearest neighbour.
 *        6/07/05       mk      Change definition of texture coords from GLint to GLflot -> Allows subpixel accurate texturing.
 *        7/23/05       mk      New options filterMode and globalAlpha for PsychBlitTextureToDisplay().
 *        8/11/05       mk      New texture handling: We now use real OpenGL textures instead of our old "pseudo-texture"
 *                              implementation. This is *way* faster, e.g., drawing times decrease from 35 ms to 3 ms for big
 *                              textures. Some experimental optimizations are implemented but not yet enabled...
 *        10/11/05      mk      Support for special Quicktime movie textures added.
 *        01/02/05      mk      Moved from OSX folder to Common folder. Contains nearly only shared code.
 *        3/07/06       awi     Print warnings conditionally according to PsychPrefStateGet_SuppressAllWarnings().
 *
 *    DESCRIPTION:
 *
 *        Psychtoolbox functions for dealing with textures.
 *
 */

#include "Screen.h"

// If set to true, then the apple client storage extensions are used: I doubt that they have any
// advantage for the current way PTB is used, but it can be useful to conserve VRAM on very
// low-mem gfx cards if Screen('Preference', 'ConserveVRAM') is set appropriately.
static psych_bool clientstorage = FALSE;

// This stores the texture format/mode to use: We autodetect available types at first
// invocation of PsychCreateTexture()... We try to use GL_EXT_TEXTURE_RECTANGLE_2D textures for
// higher speed/efficiency and lower memory consumption. If that fails, we try
// vendor specifics like GL_NV_TEXTURE_RECTANGLE, if everything fails, we resort to
// GL_TEXTURE_2D... This switch defines the global mode for the texture mapping engine...
static GLenum  texturetarget = 0;

// A rough guess of how much memory is currently consumed by textures... Can be grossly wrong,
// only used if texture creation failed and out-of-memory is a likely suspect.
static size_t texmemguesstimate = 0;

void PsychDetectTextureTarget(PsychWindowRecordType *win)
{
    // First time invocation?
    if (texturetarget==0) {
        // Yes. Need to auto-detect texturetarget to use. This routine is called with
        // the OpenGL context for the 'win' already attached, from PsychOpenOnscreenWindow().
        if ((strstr((char*) glGetString(GL_EXTENSIONS), "GL_EXT_texture_rectangle") || strstr((char*) glGetString(GL_EXTENSIONS), "GL_ARB_texture_rectangle")) && GL_TEXTURE_RECTANGLE_EXT != GL_TEXTURE_2D) {
            // Great! GL_TEXTURE_RECTANGLE_EXT is available! Use it.
            texturetarget = GL_TEXTURE_RECTANGLE_EXT;
            if(PsychPrefStateGet_Verbosity()>3)
                printf("PTB-INFO: Using OpenGL GL_TEXTURE_RECTANGLE_EXT extension for efficient high-performance texture mapping...\n");
        }
        else if (strstr((char*) glGetString(GL_EXTENSIONS), "GL_NV_texture_rectangle") && GL_TEXTURE_RECTANGLE_NV != GL_TEXTURE_2D) {
            // Try NVidia specific texture rectangle extension:
            texturetarget = GL_TEXTURE_RECTANGLE_NV;
            if(PsychPrefStateGet_Verbosity()>3)
                printf("PTB-INFO: Using NVidia's GL_TEXTURE_RECTANGLE_NV extension for efficient high-performance texture mapping...\n");
        }
        else {
            // No texture rectangle extension available :(
            // We fall back to standard power-of-two or non-power-of-two normalized textures...
            texturetarget = GL_TEXTURE_2D;

            // Warn only about lack of rectangle textures on desktop OpenGL, not on OpenGL-ES, as ES does not
            // support rectangle textures at all, so no point about making noise about it. At least OpenGL-ES-2
            // usually supports non-power-of-two textures, so we can work around many of the limitations there,
            // but no remedy for ES-1:
            if ((PsychPrefStateGet_Verbosity() > 1) && !PsychIsGLES(win)){
                printf("\nPTB-WARNING: Your graphics hardware & driver doesn't support OpenGL rectangle textures.\n");
                printf("PTB-WARNING: This won't affect the correctness or visual accuracy of image drawing, but it can significantly\n");
                printf("PTB-WARNING: degrade performance/speed and increase memory consumption of images by up to a factor of 4!\n");
                printf("PTB-WARNING: If you use a lot of image stimuli (DrawTexture, Offscreen windows, Stereo display, movies)\n");
                printf("PTB-WARNING: and you are unhappy with the performance, then please upgrade your graphics driver and possibly\n");
                printf("PTB-WARNING: your gfx hardware if you need higher performance...\n");
            }
        }
    }
    return;
}

void PsychInitWindowRecordTextureFields(PsychWindowRecordType *win)
{
    win->textureMemory=NULL;
    win->textureNumber=0;
    win->textureMemorySizeBytes=0;

    // Setup initial texture orientation: 0 = Transposed texture == Format of Matlab image matrices.
    // This number defines how the height and width of a texture need to be interpreted and how
    // texture coordinates are assigned in PsychBlitTextureToDisplay().
    win->textureOrientation=0;

    // Set to default minus 1-value (== disabled): Meaning of this field is specific for source of texture. It should
    // somehow identify the cache data structure for textures of specifif origin in a unique way.
    // If this is a cached texture for use by the PsychMovieSupport subsystem then this points to
    // the movieRecord of the movie which is associated with this texture...
    win->texturecache_slot=-1;

    // Explicit storage of the type of texture target for this texture: Zero means - Autodetect.
    win->texturetarget=0;

    // We do not have a preset for texture representation by default:
    win->textureinternalformat=0;
    win->textureexternalformat=0;
    win->textureexternaltype=0;

    // Optional GLSL filtershader: This defaults to zero i.e. no such thing. SCREENMakeTexture.c will
    // initialize this on demand: If a floating point texture is created on a piece of gfx-hardware that
    // doesn't support float-texture filtering. Then a shader is created which can reimplement that for
    // float-textures and the handle is stored in the onscreen window record and in the corresponding
    // texture. The blitters can then optionally bind that shader if filtering is requested.
    win->textureFilterShader=0;

    // Same for nearest neighbour lookup -- needed for unclamped high precision drawing.
    win->textureLookupShader=0;

    // Input pixel data for textures is not aligned to some x-byte boundary by default. A non-zero
    // setting will be used for the GL_UNPACK_ALIGNMENT setting in PsychCreateTexture() and friends
    // to optimize texture upload:
    win->textureByteAligned=0;
}

void PsychCreateTexture(PsychWindowRecordType *win)
{
    #if PSYCH_SYSTEM == PSYCH_OSX
    GLenum textureHint;
    #endif
    GLenum texturetarget, oldtexturetarget = GL_TEXTURE_RECTANGLE_EXT;
    double sourceWidth, sourceHeight;
    GLint glinternalFormat = 0, gl_realinternalformat = 0;
    static GLint gl_lastrequestedinternalFormat = 0;
    GLint gl_rbits=0, gl_gbits=0, gl_bbits=0, gl_abits=0, gl_lbits=0;
    int twidth, theight, pass, texcount;
    void* texmemptr;
    psych_bool recycle = FALSE, avoidCPUGPUSync;
    GLenum glerr;
    int verbosity;

    verbosity = PsychPrefStateGet_Verbosity();

    // Check if any calls that can cause a CPU<->GPU sync should be avoided at (nearly all) costs.
    // Avoiding CPU<->GPU sync is a robustness vs. performance tradeoff: Performance is potentially
    // significantly increased, but the amount of error checking and error handling is drastically
    // reduced --> Fast, but potential errors or malfunctions don't get catched and diagnosed in a
    // useful way, instead silent failure occurs.
    avoidCPUGPUSync = (PsychPrefStateGet_ConserveVRAM() & kPsychAvoidCPUGPUSync) ? TRUE : FALSE;

    // Enable the proper OpenGL rendering context for the window associated with this texture:
    PsychSetGLContext(win);

    // Make sure we don't have any dangling GL errors from other operations...
    if (!avoidCPUGPUSync || (verbosity > 10)) PsychTestForGLErrors();

    // Setup texture-target if not already done:
    PsychDetectTextureTarget(win);

    // Assign proper texturetarget for creation:
    texturetarget = PsychGetTextureTarget(win);

    // Check if user requested explicit use of clientstorage + Use of System RAM for
    // storage of textures instead of VRAM caching in order to conserve VRAM memory on
    // low-mem gfx-cards. Enable clientstorage, if so...
    clientstorage = (PsychPrefStateGet_ConserveVRAM() & kPsychDontCacheTextures) ? TRUE : FALSE;

    // Create a unique texture handle for this texture:
    // If the texture already has a handle assigned then this means that we shouldn't
    // create and setup a new OpenGL texture from scratch, but bind and recycle the
    // given texture object. Just bind it and update its content via glTexSubImage()... calls.
    // Updating textures is potentially faster than recreating them -> movie playback
    // and the Videocapture code et al. will benefit from this...
    if (win->textureNumber == 0) {
        glGenTextures(1, &win->textureNumber);
        recycle = FALSE;
        //printf("CREATING NEW TEX %i\n", win->textureNumber);
    }
    else {
        recycle = TRUE;
        //printf("RECYCLING TEX %i\n", win->textureNumber);
    }

    // Setup texturing:
    glDisable(GL_TEXTURE_2D);
    glEnable(texturetarget);

    // Create & bind a new OpenGL texture object and attach it to our new texhandle:
    glBindTexture(texturetarget, win->textureNumber);

    // Setup texture parameters like optimization, storage format et al.

    // Choose the texture acceleration extension out of GL_STORAGE_PRIVATE_APPLE, GL_STORAGE_CACHED_APPLE, GL_STORAGE_SHARED_APPLE
    // We normally use CACHED storage for caching textures in gfx-cards VRAM for high-perf drawing,
    // but if user explicitely requests client storage for saving VRAM memory, we do so and
    // use SHARED storage in system RAM --> Slower but saves VRAM memory.
    #if PSYCH_SYSTEM == PSYCH_OSX
        textureHint= (clientstorage) ? GL_STORAGE_SHARED_APPLE : GL_STORAGE_CACHED_APPLE;
        glTexParameteri(texturetarget, GL_TEXTURE_STORAGE_HINT_APPLE , textureHint);
        glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, (clientstorage) ? GL_TRUE : GL_FALSE);
    #endif

    // Not using GL_STORAGE_SHARED_APPLE provided increased reliability of timing and significantly shorter rendering times
    // when testing with a G5-Mac with 1.6Ghz CPU and 256 MB RAM, MacOS-X 10.3.7, using 12 textures of 800x800 pixels each.
    // Rendering times with this code are around 6 msecs, while original PTB 1.0.40 code took 17 ms on average...
    // Can't test this on other machines like G4 or machines with more RAM...
    // -> Sometimes, GL_STORAGE_SHARED_APPLE is faster, but only if texture width and height are divisable by 16 and
    // all used memory is page-aligned and a couple other conditions are met. So... Sometimes you are 20% faster, but most of
    // the time you are 2 to 3 times slower than without this extensions...
    // The "Principle of least surprise" would suggest to disable the extension, because the end-user doesn't
    // expect sudden and random changes in performance of his PTB scripts.
    // Alternatively one could code up different path's depending on if the preconditions are met or not...
    //
    // We could reenable the extension, if wanted, but then the MakeTexture code needs to be modified in a way
    // that will slow down MakeTexture a bit. It's a tradeoff between speed of DrawTexture and speed of MakeTexture.
    //
    // BTW -> Does disabling the extension solve "severe tearing bug" reported in Forum message 3007?
    // Explanation: GL_STORAGE_SHARED_APPLE enables texture fetches over AGP bus via DMA operations and
    // should increase performance. But DMA only triggers when texture width is divisible by 8, in all other
    // cases it's disabled. Bug in message 3007 only happens when texture width is divisible by 8. Could this
    // be a bug in the G4 Laptops graphics hardware (DMA-Engine) or in its OpenGL driver???
    // Would be interesting to find out...

    // Setting GL_UNPACK_ALIGNMENT == 1 fixes a bug, where textures are drawn incorrectly, if their
    // width or height is not divisible by 4.
    if (win->textureByteAligned > 1) {
        // Aligned to some byte boundary: Set alignment to allow driver to optimize texture upload:
        glPixelStorei(GL_UNPACK_ALIGNMENT, win->textureByteAligned);
    }
    else {
        // No alignment (1) or alignment unknown (0): Set to safe default of 1 byte:
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    }

    // The texture object is ready for use: Assign it our texture data:

    // Definition of width and height is swapped due to texture rotation trick, see comments in PsychBlit.....
    if (win->textureOrientation==0 || win->textureOrientation==1) {
        // Transposed case: Optimized for fast MakeTexture from Matlab image matrix.
        // This is true for all calls from MakeTexure, except ones with the textureOrientation flag set to 2:
        sourceHeight=PsychGetWidthFromRect(win->rect);
        sourceWidth=PsychGetHeightFromRect(win->rect);
    }
    else {
        // Non-transposed upright case: This is used for textures created by 'OpenOffscreenWindow'
        // One can directly draw to these textures as rendertargets aka OpenGL framebuffer objects...
        sourceHeight=PsychGetHeightFromRect(win->rect);
        sourceWidth=PsychGetWidthFromRect(win->rect);
    }

    // Desktop-GL only:
    if (!PsychIsGLES(win)) glPixelStorei(GL_UNPACK_ROW_LENGTH, (int) sourceWidth);

    // We used to have different cases for Luminance, Luminance+Alpha, RGB, RGBA.
    // This way we saved texture memory for the source->textureMemory -- Arrays, as well as copy-time
    // in MakeTexture - In theory...
    // Reality is: We always use GL_RGBA8 as internal format, except for pure luminance textures.
    // This obviously wastes storage space for LA and RGB textures, but it is the only mode that is
    // well supported (=fast) on all common gfx-hardware. Only the very latest models of NVidia and ATI
    // are capable of handling the other formats natively in hardware :-(
    if ((texturetarget == GL_TEXTURE_2D) && !(win->gfxcaps & kPsychGfxCapNPOTTex)) {
        // This hardware doesn't support non-power-of-two GL_TEXTURE_2D textures. We create and use power of two
        // textures to emulate rectangle textures...

        // Compute smallest power of two dimension that fits the texture.
        twidth=1;
        while (twidth<sourceWidth) twidth*=2;
        theight=1;
        while (theight<sourceHeight) theight*=2;
        // First we only use glTexImage2D with NULL data pointer to create a properly sized empty texture:
        texmemptr=NULL;
    }
    else {
        // Hardware supports rectangular textures: Use texture as-is:
        twidth  = (int) sourceWidth;
        theight = (int) sourceHeight;
        texmemptr=win->textureMemory;
    }

    // We only execute this pass for really new textures, not for recycled ones:
    if (!recycle) {
        // This is a two-pass procedure. First we check with a proxy-texture if texture
        // creation will succeed without trouble, then we either fail in case of error,
        // or we do the real thing and create the texture:
        for (pass=0; pass < 2; pass++) {
            // Skip first probing pass if usercode wants it, or if we're on OpenGL-ES, as
            // ES does not support proxy textures, so we can't do the probing:
            if ((PsychIsGLES(win) || (avoidCPUGPUSync && !(verbosity > 10))) && (pass == 0)) {
                // Usercode wants us to skip dual-path texture creation in order to save
                // some time: Skip pass 0...
                oldtexturetarget = texturetarget;
                continue;
            }

            if (pass == 0) {
                // Prepare proxy-pass:
                oldtexturetarget = texturetarget;
                texturetarget = (texturetarget == GL_TEXTURE_2D) ? GL_PROXY_TEXTURE_2D : GL_PROXY_TEXTURE_RECTANGLE_ARB;
            }
            else {
                // Restore real texture target from saved one in pass 1:
                texturetarget = oldtexturetarget;
            }

            if (win->textureinternalformat == 0) {
                // Standard path: Derive texture format and such from requested pixeldepth:
                switch(win->depth) {
                    case 8:
                        glinternalFormat = PsychIsGLES(win) ? GL_LUMINANCE : GL_LUMINANCE8;
                        glTexImage2D(texturetarget, 0, glinternalFormat, (GLsizei) twidth, (GLsizei) theight, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, texmemptr);
                        break;

                    case 16:
                        //glinternalFormat=GL_LUMINANCE8_ALPHA8;
                        glinternalFormat = PsychIsGLES(win) ? GL_LUMINANCE_ALPHA : GL_RGBA8;
                        glTexImage2D(texturetarget, 0, glinternalFormat, (GLsizei) twidth, (GLsizei) theight, 0, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, texmemptr);
                        break;

                    case 24:
                        //glinternalFormat=GL_RGB8;
                        glinternalFormat = PsychIsGLES(win) ? GL_RGB : GL_RGBA8;
                        glTexImage2D(texturetarget, 0, glinternalFormat, (GLsizei) twidth, (GLsizei) theight, 0, GL_RGB, GL_UNSIGNED_BYTE, texmemptr);
                        break;

                    case 32:
                        if (PsychIsGLES(win)) {
                            // GLES is much more restricted:
                            if (strstr((const char*) glGetString(GL_EXTENSIONS), "GL_EXT_texture_format_BGRA8888")) {
                                glinternalFormat = GL_BGRA_EXT;
                                glTexImage2D(texturetarget, 0, glinternalFormat, (GLsizei) twidth, (GLsizei) theight, 0, GL_BGRA_EXT, GL_UNSIGNED_BYTE, texmemptr);
                            }
                            else {
                                glinternalFormat = GL_RGBA;
                                glTexImage2D(texturetarget, 0, glinternalFormat, (GLsizei) twidth, (GLsizei) theight, 0, GL_RGBA, GL_UNSIGNED_BYTE, texmemptr);
                            }
                        }
                        else {
                            // Classic path with optimizations:
                            glinternalFormat = GL_RGBA8;
                            glTexImage2D(texturetarget, 0, glinternalFormat, (GLsizei) twidth, (GLsizei) theight, 0, GL_BGRA, ((win->gfxcaps & kPsychGfxCapNeedsUnsignedByteRGBATextureUpload) ? GL_UNSIGNED_BYTE : GL_UNSIGNED_INT_8_8_8_8_REV), texmemptr);
                        }
                        break;
                }
            }
            else {
                // Requested internal format and external data representation are explicitely requested: Use it.
                glTexImage2D(texturetarget, 0, win->textureinternalformat, (GLsizei) twidth, (GLsizei) theight, 0, win->textureexternalformat,
                             win->textureexternaltype, texmemptr);
                glinternalFormat = win->textureinternalformat;
            }

            if (!PsychIsGLES(win)) {
                // Only query real color depths per channel if either a special internal format was requested, ie., we
                // need to query the real depths because we don't know it, or if verbosity is very high or avoidance
                // of synchronizing calls is not disabled.
                if ((!avoidCPUGPUSync || (verbosity > 10)) || ((glinternalFormat!=GL_RGBA8) && (glinternalFormat!=GL_LUMINANCE8))) {
                    // Request sizes of created (proxy-)texture:
                    glGetTexLevelParameteriv(texturetarget, 0, GL_TEXTURE_RED_SIZE, &gl_rbits);
                    glGetTexLevelParameteriv(texturetarget, 0, GL_TEXTURE_GREEN_SIZE, &gl_gbits);
                    glGetTexLevelParameteriv(texturetarget, 0, GL_TEXTURE_BLUE_SIZE, &gl_bbits);
                    glGetTexLevelParameteriv(texturetarget, 0, GL_TEXTURE_ALPHA_SIZE, &gl_abits);
                    glGetTexLevelParameteriv(texturetarget, 0, GL_TEXTURE_LUMINANCE_SIZE, &gl_lbits);

                    // Special override for YCBCR textures which return all component bits as zero:
                    if (glinternalFormat == GL_YCBCR_MESA || win->textureexternalformat == GL_YCBCR_422_APPLE) gl_rbits = 8;

                    // Store override per-component bit-depths:
                    win->bpc = (int) ((gl_rbits > gl_lbits) ? gl_rbits : gl_lbits);
                }
                else {
                    // Null out to catch implementation error in our own code, if any:
                    gl_rbits = gl_gbits = gl_bbits = gl_abits = gl_lbits = 0;

                    // Assume 8 bpc. This is the default anyway, just here to make it explicit:
                    win->bpc = 8;
                }
            }
            else {
                // OpenGL-ES: Can't query bit depth of texture. Assume 8 (integer) or 32 (float) bits per channel:
                win->bpc = (win->textureexternaltype == GL_FLOAT) ? 32 : 8;
                gl_rbits = gl_gbits = gl_bbits = gl_abits = gl_lbits = win->bpc;
            }

            if (PsychIsGLES(win)) {
                // Our error checking for ES is very limited:
                if ((!avoidCPUGPUSync || (verbosity > 10)) && ((glerr = glGetError()) !=0 )) {
                    glBindTexture(texturetarget, 0);
                    win->textureNumber = 0;
                    if (!clientstorage) {
                        if (win->textureMemory && (win->textureMemorySizeBytes > 0)) free(win->textureMemory);
                        win->textureMemory=NULL;
                        win->textureMemorySizeBytes=0;
                    }
                    while(glGetError()) {};

                    if (verbosity > 0) {
                        printf("\n\nPTB-ERROR: Texture creation failed or malfunctioned for a texture of requested size w x h = %i x %i texels\n", twidth, theight);
                        printf("PTB-ERROR: and at least %i bytes VRAM memory consumption per texel.\n", win->depth / 8);
                        if (glerr == GL_INVALID_VALUE) {
                            glGetIntegerv(GL_MAX_TEXTURE_SIZE, &gl_rbits);
                            if (gl_rbits < twidth || gl_rbits < theight) {
                                // Hard hardware limit exceeded:
                                printf("PTB-ERROR: Your image or texture exceeds the maximum width and/or height of %i texels supported by your graphics hardware.\n", gl_rbits);
                                printf("PTB-ERROR: You'll have to either reduce the size of your images below that limit, or upgrade your hardware.\n\n");
                            }
                            else {
                                printf("PTB-ERROR: The hardware likely does not support the requested texture storage format.\n");
                            }
                        }
                        else {
                            printf("PTB-ERROR: The gfx-hardware reports the following error: %s\n", gluErrorString(glerr));
                        }
                    }
                    PsychErrorExitMsg(PsychError_user, "Texture creation failed, most likely due to unsupported precision or insufficient free memory.");
                }
            }

            // Checking for desktop OpenGL - much more elaborate...
            // Sanity check: A sum of zero over all texture channels would indicate that texture
            // creation failed, most likely due to an out of memory condition in the graphics
            // hardwares VRAM:
            else if ((!avoidCPUGPUSync || (verbosity > 10)) && (((glerr = glGetError())!=0) || (gl_rbits + gl_gbits + gl_bbits + gl_abits + gl_lbits == 0))) {
                // Texture creation failed or malfunctioned!
                if (verbosity > 0) {
                    // Abort with error:

                    // First eat up any pending GL errors to make sure our shutdown path doesn't fail:
                    while(glGetError());

                    // Free all ressources already allocated for this failed texture creation request:
                    glBindTexture(texturetarget, 0);
                    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
                    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
                    glDeleteTextures(1, &win->textureNumber);
                    win->textureNumber = 0;
                    if (!clientstorage) {
                        if (win->textureMemory && (win->textureMemorySizeBytes > 0)) free(win->textureMemory);
                        win->textureMemory=NULL;
                        win->textureMemorySizeBytes=0;
                    }

                    while(glGetError());

                    if (glerr == GL_INVALID_VALUE || (gl_rbits + gl_gbits + gl_bbits + gl_abits + gl_lbits == 0)) {
                        // Most likely texture too big for implementation or out of memory condition in VRAM or unsupported format:
                        printf("\n\nPTB-ERROR: Texture creation failed or malfunctioned for a texture of requested size w x h = %i x %i texels\n", twidth, theight);
                        printf("PTB-ERROR: and at least %i bytes VRAM memory consumption per texel.\n", (glinternalFormat==GL_RGBA8) ? 4 : win->depth / 8);

                        // Query maximum size of textures supported by hardware:
                        glGetIntegerv((texturetarget == GL_TEXTURE_2D || texturetarget == GL_PROXY_TEXTURE_2D) ? GL_MAX_TEXTURE_SIZE : GL_MAX_RECTANGLE_TEXTURE_SIZE_ARB, &gl_rbits);
                        if (gl_rbits < twidth || gl_rbits < theight) {
                            // Hard hardware limit exceeded:
                            printf("PTB-ERROR: Your image or texture exceeds the maximum width and/or height of %i texels supported by your graphics hardware.\n", gl_rbits);
                            printf("PTB-ERROR: You'll have to either reduce the size of your images below that limit, or upgrade your hardware.\n\n");
                        }
                        else {
                            // Either out-of-memory in VRAM for such large textures, or unsupported format/precision:
                            if (glinternalFormat!=GL_LUMINANCE8 && glinternalFormat!=GL_RGBA8 && !glewIsSupported("GL_APPLE_float_pixels") && !glewIsSupported("GL_ATI_texture_float") && !glewIsSupported("GL_ARB_texture_float")) {
                                // Requested format is not one of the 8bpc fixed-point LDR formats, but a HDR format which
                                // doesn't seem to be supported by the GL implementation:
                                printf("PTB-ERROR: The image is of a high precision (HDR) format, not of a standard 8bpc (LDR) format.\n");
                                printf("PTB-ERROR: Seems that such texture formats are not supported by your graphics hardware. You'll need to\n");
                                printf("PTB-ERROR: update your graphics driver and probably your graphics hardware to make use of HDR textures.\n");
                            }
                            else {
                                if (glinternalFormat!=GL_LUMINANCE8 && glinternalFormat!=GL_RGBA8) {
                                    printf("PTB-ERROR: The image is of a high precision (HDR) format, not of a standard 8bpc (LDR) format.\n");
                                    printf("PTB-ERROR: Such HDR textures have very high VRAM memory demands.\n");
                                }

                                // Check if a general out-of-memory condition is the likely culprit, due to wrong texture management
                                // on the users side:
                                texcount = PsychRessourceCheckAndReminder(FALSE);
                                printf("PTB-ERROR: Currently there are already %i textures, offscreen windows, movies or proxies open.\n", texcount);
                                printf("PTB-ERROR: All these objects consume system memory and could lead to ressource shortage.\n");
                                printf("PTB-ERROR: My current (rough and probably way too low) estimate is that at least %f MB of memory are\n", (float) texmemguesstimate / 1024 / 1024);
                                printf("PTB-ERROR: consumed for textures, offscreen windows and similar objects.\n");
                                if (texcount > 100) {
                                    printf("PTB-ERROR: The count is above one hundred objects. Could it be that you forgot to dispose no longer\n");
                                    printf("PTB-ERROR: needed objects from previous experiment trials (missing Screen('Close' [, texturePtr]) or Screen('CloseMovie', moviePtr))??\n");
                                }

                                if (texmemguesstimate > 100 * 1024 * 1024) {
                                    printf("PTB-ERROR: At least 100 MB memory consumed for textures, probably much more. Could it be that you forgot to dispose no longer\n");
                                    printf("PTB-ERROR: needed textures from previous experiment trials (missing Screen('Close' [, texturePtr]))??\n");
                                }

                                printf("PTB-ERROR: Another cause of failure could be that your graphics hardware doesn't have sufficient amounts of\n");
                                printf("PTB-ERROR: free VRAM memory. Try to reduce the precision and/or size of your texture image to the lowest\n");
                                printf("PTB-ERROR: acceptable setting for your purpose.\n");
                                printf("PTB-ERROR: Read the online help for Screen MakeTexture? or Screen OpenOffscreenWindow? for information\n");
                                printf("PTB-ERROR: about how to reduce the number of color channels and/or precision of the texture.\n");
                                printf("PTB-ERROR: It may also help to reduce general VRAM memory consumption:\n");
                                printf("PTB-ERROR: Reduce the display resolution to the smallest acceptable resolution or disable the 2nd\n");
                                printf("PTB-ERROR: display of a dual-display setup if it isn't strictly needed for your study.\n");
                                printf("PTB-ERROR: Other than that you could upgrade to more powerful graphics hardware (more VRAM) or try to split\n");
                                printf("PTB-ERROR: your oversized image into multiple separate subimages and draw them sequentially to the screen.\n");
                                printf("PTB-ERROR: That way you'll be able to draw larger images, albeit at drastically lowered performance and more\n");
                                printf("PTB-ERROR: coding hazzle.\n\n");
                            }
                        }

                        PsychErrorExitMsg(PsychError_user, "Texture creation failed, most likely due to unsupported precision or insufficient free memory.");
                    }
                    else {
                        if ((glerr == GL_INVALID_ENUM) && (glinternalFormat == GL_YCBCR_MESA || win->textureexternalformat == GL_YCBCR_422_APPLE)) {
                            printf("PTB-ERROR: Tried to use a memory-optimized YCBYCR texture, but your GPU + graphics driver doesn't support this.\n");
                            printf("PTB-ERROR: You'll need to update your graphics driver and/or upgrade your graphics card, or change your script code\n");
                            printf("PTB-ERROR: to avoid this special unsupported texture format.\n");
                            PsychErrorExitMsg(PsychError_user, "Texture creation failed, most likely due to use of an unsupported texture format.");
                        }

                        // Some other error:
                        printf("\n\nPTB-ERROR: Texture creation failed! OpenGL reported the following error condition: %s.\n", gluErrorString(glerr));
                        PsychErrorExitMsg(PsychError_user, "Texture creation failed for unknown reason. You may want to contact the Psychtoolbox forum for help.");
                    }
                }
            }    // End of error checking...
        }  // End of dual-pass texture creation (check + create).

        // Accounting... ...this is only a rough guesstimate:
        win->surfaceSizeBytes = ((size_t) ((glinternalFormat==GL_RGBA8) ? 4 : win->depth / 8)) * (size_t) twidth * (size_t) theight;
        texmemguesstimate+= win->surfaceSizeBytes;
    }  // End of new texture creation.

    // Stage 2: If it is a 2D texture or a recycled texture, fill it with content via glTexSubImage2D:
    if (texturetarget == GL_TEXTURE_2D || recycle) {
        // Special setup code for pot2 textures: Fill the empty power of two texture object with content:
        // We only fill a subrectangle (of sourceWidth x sourceHeight size) with our images content. The
        // unused border contains all zero == black.
        // The same path is used for efficient refilling existing textures that are to be recycled:
        if (win->textureinternalformat==0) {
            // Standard path: Derive texture format and such from requested pixeldepth:
            switch(win->depth) {
                case 8:
                    glTexSubImage2D(texturetarget, 0, 0, 0, (GLsizei)sourceWidth, (GLsizei)sourceHeight, GL_LUMINANCE, GL_UNSIGNED_BYTE, win->textureMemory);
                    break;

                case 16:
                    glTexSubImage2D(texturetarget, 0, 0, 0, (GLsizei)sourceWidth, (GLsizei)sourceHeight, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, win->textureMemory);
                    break;

                case 24:
                    glTexSubImage2D(texturetarget, 0, 0, 0, (GLsizei)sourceWidth, (GLsizei)sourceHeight, GL_RGB, GL_UNSIGNED_BYTE, win->textureMemory);
                    break;

                case 32:
                    if (PsychIsGLES(win)) {
                        // GLES is much more restricted:
                        if (strstr((const char*) glGetString(GL_EXTENSIONS), "GL_EXT_texture_format_BGRA8888")) {
                            glTexSubImage2D(texturetarget, 0, 0, 0, (GLsizei)sourceWidth, (GLsizei)sourceHeight, GL_BGRA_EXT, GL_UNSIGNED_BYTE, win->textureMemory);
                        }
                        else {
                            glTexSubImage2D(texturetarget, 0, 0, 0, (GLsizei)sourceWidth, (GLsizei)sourceHeight, GL_RGBA, GL_UNSIGNED_BYTE, win->textureMemory);
                        }
                    }
                    else {
                        // Classic path:
                        glTexSubImage2D(texturetarget, 0, 0, 0, (GLsizei)sourceWidth, (GLsizei)sourceHeight, GL_BGRA, ((win->gfxcaps & kPsychGfxCapNeedsUnsignedByteRGBATextureUpload) ? GL_UNSIGNED_BYTE : GL_UNSIGNED_INT_8_8_8_8_REV), win->textureMemory);
                    }
                    break;
            }
        }
        else {
            // Requested internal format and external data representation are explicitely requested: Use it.
            glTexSubImage2D(texturetarget, 0, 0, 0, (GLsizei)sourceWidth, (GLsizei)sourceHeight, win->textureexternalformat, win->textureexternaltype, win->textureMemory);
            glinternalFormat = win->textureinternalformat;
        }
    }

    if (!PsychIsGLES(win)) {
        // New internal format requested?
        if ((!avoidCPUGPUSync || (verbosity > 10)) && (!recycle && (gl_lastrequestedinternalFormat != glinternalFormat))) {
            // Seems so...
            gl_lastrequestedinternalFormat = glinternalFormat;

            // Query real internal format and params...
            glGetTexLevelParameteriv(texturetarget, 0, GL_TEXTURE_INTERNAL_FORMAT, &gl_realinternalformat);
            // If there is a mismatch between wish and reality, report it:
            if (false || gl_realinternalformat != glinternalFormat) {
                // Mismatch between requested format and format that the OpenGL has chosen:
                printf("PTB-WARNING: In glTexImage2D: Mismatch between requested and real format: depth=%i, fcode=x%x\n", win->depth, gl_realinternalformat);
                printf("PTB-WARNING: Requested size = %i bits, real size = %i bits.\n", win->depth, gl_rbits + gl_gbits + gl_bbits + gl_abits + gl_lbits);
                printf("PTB-WARNING: This could mean that something went wrong when creating the texture!\n");
                fflush(NULL);
            }
        }
    }

    // Free system RAM backing memory buffer, if client storage extensions are not used for this texture:
    if (!clientstorage) {
        if (win->textureMemory && (win->textureMemorySizeBytes > 0)) free(win->textureMemory);
        win->textureMemory=NULL;
        win->textureMemorySizeBytes=0;
    }

    // Texture object ready for future use. Unbind it:
    glBindTexture(texturetarget, 0);
    glDisable(texturetarget);

    // Reset pixel storage parameter:
    if (!PsychIsGLES(win)) glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // Client rect of a texture is always == rect of it:
    PsychCopyRect(win->clientrect, win->rect);

    // Finished!
    return;
}

/*
 *    PsychFreeTextureForWindowRecord()
 *
 *    Accepts a window record for either a texture or an onscreen window and deallocte texture resources.
 *
 */
void PsychFreeTextureForWindowRecord(PsychWindowRecordType *win)
{
    // Destroy OpenGL texture object for windows that have one:
    if((win->windowType==kPsychSingleBufferOnscreen || win->windowType==kPsychDoubleBufferOnscreen || win->windowType==kPsychTexture) &&
        (win->targetSpecific.contextObject)) {
        // Activate associated OpenGL context:
        PsychSetGLContext(win);

        // PsychTestForGLErrors() is a GPU-CPU synchronization point, so in order to keep good
        // parallelism, we only do it at verbosity levels of 5 and greater.
        if (PsychPrefStateGet_Verbosity() > 4) PsychTestForGLErrors();

        // Call special texture release routine for Movie textures: This routine will
        // check if 'win' is a movie texture and perform the necessary cleanup work, if so:
        PsychFreeMovieTexture(win);

        // If we use client-storage textures, we need to wait for completion of texture operations on the
        // to-be-released client texture buffer before deleting it and freeing the RAM backing buffers. Waiting for
        // completion is done via FinishObjectApple...
        // We need to use glFinish() here. FinishObjectApple would be better (more async operations) but it doesn't
        // work for some strange reason :(
        if ((win->textureMemory) && (win->textureNumber > 0)) glFinish(); // FinishObjectAPPLE(GL_TEXTURE_2D, win->textureNumber);

        // Perform standard OpenGL texture cleanup if needed:
        if (win->textureNumber != 0) {
            glDeleteTextures(1, &win->textureNumber);

            // Accounting... ...this is only a rough guesstimate:
            texmemguesstimate-= win->surfaceSizeBytes;
        }

        // PsychTestForGLErrors() is a GPU-CPU synchronization point, so in order to keep good
        // parallelism, we only do it at verbosity levels of 5 and greater.
        if (PsychPrefStateGet_Verbosity() > 4) PsychTestForGLErrors();
    }

    // Free system RAM backing memory buffer, if any:
    if (win->textureMemory) free(win->textureMemory);
    win->textureMemory=NULL;
    win->textureMemorySizeBytes=0;
    win->textureNumber=0;

    return;
}


void PsychBlitTextureToDisplay(PsychWindowRecordType *source, PsychWindowRecordType *target, double *sourceRect, double *targetRect,
                               double rotationAngle, int filterMode, double globalAlpha)
{
    int tWidth = 0, tHeight = 0;
    GLdouble sourceWidth, sourceHeight;
    GLdouble sourceX, sourceY, sourceXEnd, sourceYEnd;
    double transX, transY;
    GLenum texturetarget;
    GLint attrib;
    GLuint shader = 0;

    // Enable targets framebuffer as current drawingtarget, except if this is a
    // blit operation from a window into itself and the imaging pipe is on:
    if ((source != target) || (target->imagingMode==0)) {
        PsychSetDrawingTarget(target);
    }
    else {
        // Activate rendering context of target window without changing drawing target:
        PsychSetGLContext(target);
    }

    // Setup texture-target if not already done:
    PsychDetectTextureTarget(target);

    // Query target for this specific texture:
    texturetarget = PsychGetTextureTarget(source);

    //printf("%i\n", source->textureOrientation);

    // This code allows the application of sourceRect, as it is meant to be:
    // CAUTION: This calculation with sourceHeight - xxxx  depends on if GPU texture swapping
    // is on or off!!!!
    // 0 == Transposed as from Matlab image array. 2 == Offscreen window in normal orientation.
    if (source->textureOrientation == 2) {
        sourceHeight=PsychGetHeightFromRect(source->rect);
        sourceWidth=PsychGetWidthFromRect(source->rect);

        sourceX=sourceRect[kPsychLeft];
        sourceY=sourceHeight - sourceRect[kPsychBottom];
        sourceXEnd=sourceRect[kPsychRight];
        sourceYEnd=sourceHeight - sourceRect[kPsychTop];
    }
    else {
        sourceHeight=PsychGetWidthFromRect(source->rect);
        sourceWidth=PsychGetHeightFromRect(source->rect);
        sourceX=sourceRect[kPsychTop];
        sourceY=sourceRect[kPsychLeft];
        sourceXEnd=sourceRect[kPsychBottom];
        sourceYEnd=sourceRect[kPsychRight];
    }

    // Overrides for special cases: Upside-down texture.
    if (source->textureOrientation == 3) {
        sourceHeight=PsychGetHeightFromRect(source->rect);
        sourceWidth=PsychGetWidthFromRect(source->rect);
        sourceX=sourceRect[kPsychLeft];
        sourceY=sourceRect[kPsychBottom];
        sourceXEnd=sourceRect[kPsychRight];
        sourceYEnd=sourceRect[kPsychTop];
    }

    // This case can happen with some QT movies, they are upside down in an unusual way:
    if (source->textureOrientation == 4) {
        sourceHeight=PsychGetHeightFromRect(source->rect);
        sourceWidth=PsychGetWidthFromRect(source->rect);
        sourceX=sourceRect[kPsychLeft];
        sourceY=sourceHeight - sourceRect[kPsychBottom];
        sourceXEnd=sourceRect[kPsychRight];
        sourceYEnd=sourceHeight - sourceRect[kPsychTop];
    }

    // Special case handling for GL_TEXTURE_2D textures. We need to map the
    // absolute texture coordinates (in pixels) to the interval 0.0 - 1.0 where
    // 1.0 == full extent of power of two texture...
    if (texturetarget == GL_TEXTURE_2D) {
        // NPOT supported?
        if (!(source->gfxcaps & kPsychGfxCapNPOTTex)) {
            // No: Find size of real underlying texture (smallest power of two which is
            // greater than or equal to the image size:
            tWidth=1;
            while (tWidth < (int) sourceWidth) tWidth*=2;
            tHeight=1;
            while (tHeight < (int) sourceHeight) tHeight*=2;
        }
        else {
            // Yes:
            tWidth = (int) sourceWidth;
            tHeight = (int) sourceHeight;
        }

        // Remap texcoords into 0-1 subrange: We subtract 0.5 pixel-units before
        // mapping to accomodate for roundoff-error in the power-of-two gfx
        // hardware...
        // For a good intro into the issue of texture border seams, due to interpolation
        // problems at texture borders, see:
        // http://home.planet.nl/~monstrous/skybox.html

        sourceXEnd-=0.5f;
        sourceYEnd-=0.5f;

        // Remap:
        sourceX=sourceX / (double) tWidth;
        sourceXEnd=sourceXEnd / (double) tWidth;
        sourceY=sourceY / (double) tHeight;
        sourceYEnd=sourceYEnd / (double) tHeight;
    }

    // MK: We need to reenable the proper texturing mode. This fixes bug reported in Forum message 3055,
    // because SCREENDrawText glDisable'd GL_TEXTURE_RECTANGLE_EXT, without this routine reenabling it.
    glDisable(GL_TEXTURE_2D);

    // Only enable actual texture hardware if a real texture is provided.
    // In the case of no real texture, we don't bind a real texture, don't
    // enable texture mapping and just blit the quad, with interpolated
    // texture coordinates set up for purely procedural shading.
    if (source->textureNumber > 0) {
        glEnable(texturetarget);
        glBindTexture(texturetarget, source->textureNumber);
    }

    // Use of OpenGL mip-mapping requested? And automatic mipmap generation wanted - aka not forbidden?
    if (!(source->specialflags & kPsychDontAutoGenMipMaps) && (filterMode < 0 || filterMode > 1)) {
        // Yes: Automatically build a mip-map pyramid.
        if (texturetarget != GL_TEXTURE_2D) PsychErrorExitMsg(PsychError_user, "You asked me to use mip-mapped texture filtering on a texture that is not of GL_TEXTURE_2D type! Unsupported.");

        if (NULL != glGenerateMipmapEXT) {
            GLint mipattrib = 0;

            // Select highest quality downsampling method:
            glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);

            // Is automatic mipmap generation already enabled for this texture?
            glGetTexParameteriv(texturetarget, GL_GENERATE_MIPMAP, &mipattrib);

            // Need to manually trigger regen if automatic mode not yet enabled, or
            // if the "dirty" flag is set, because some render-to-texture activity has
            // happened since last drawing this texture, which is not covered by the auto-update:
            if (source->needsViewportSetup || !mipattrib) {
                // No: We trigger hardware-accelerated mipmap generation manually for this draw call:
                glGenerateMipmapEXT(texturetarget);

                // Enable automatic mipmap generation for future updates to this texture object. This
                // will automatically trigger regen if new image content is uploaded into this texture
                // object:
                glTexParameteri(texturetarget, GL_GENERATE_MIPMAP, GL_TRUE);

                // Clear "dirty" flag:
                source->needsViewportSetup = FALSE;
            }
        }
        else if (PsychPrefStateGet_Verbosity() > 1) {
            printf("PTB-WARNING: Was asked to draw a texture with mip-mapping, but automatic mipmap generation unsupported by this system! Check your stimulus!\n");
        }
    }

    // Linear filtering on non-capable hardware via shader emulation?
    if ((filterMode != 0) && (source->textureFilterShader > 0)) {
        // Yes. Bind the shader:
        shader = source->textureFilterShader;
        if (0 == PsychSetShader(target, shader)) PsychErrorExitMsg(PsychError_user, "Tried to use a bilinear texture filter shader, but your hardware doesn't support GLSL shaders.");

        if (filterMode < 0 || filterMode > 1) {
            // Some mip-mapped filtermode. Choose nearest neighbour sampling within mipmap levels, so shader can decide about sample locations itself.
            // In filterMode 2 choose the nearest mipmap, in others interpolate linearly between two nearest mipmap levels:
            glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, (filterMode == 2) ? GL_NEAREST_MIPMAP_NEAREST : GL_NEAREST_MIPMAP_LINEAR);
            glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        }
        else {
            // No mip-mapping: Switch hardware samplers into nearest neighbour mode so we don't get any interference:
            glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        }

        // Don't restrict mipmap-levels for sampling, reset to initial system defaults:
        // This makes even sense for negative filterMode arguments, because the filterMode
        // parameter is passed as an attribute to the filtershader, so the shader itself can
        // decide how to implement a specific blur level on its own, unrestricted by us:
        if ((texturetarget == GL_TEXTURE_2D) && !PsychIsGLES(source)) {
            glTexParameteri(texturetarget, GL_TEXTURE_BASE_LEVEL, 0);
            glTexParameteri(texturetarget, GL_TEXTURE_MAX_LEVEL,  1000);
        }
    }
    else {
        // Standard hardware texture sampling/filtering: Select filter-mode for texturing:
        if (filterMode >= 0) {
            // Select specific hardware sampling strategy:
            switch (filterMode) {
                case 0: // Nearest-Neighbour filtering:
                    glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
                    glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                    break;

                case 1: // Bilinear filtering:
                    glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
                    glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                    break;

                case 2: // Linear filtering with nearest neighbour mipmapping:
                    glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
                    glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                    break;

                case 3: // Linear filtering with linear mipmapping --> This is full trilinear filtering:
                    glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
                    glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                    break;

                case 4: // Nearest-Neighbour filtering with nearest neighbour mipmapping:
                    glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
                    glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                    break;

                case 5: // Nearest-Neighbour filtering with linear mipmapping:
                    glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_LINEAR);
                    glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                    break;
            }

            // Don't restrict mipmap-levels for sampling, reset to initial system defaults:
            if ((texturetarget == GL_TEXTURE_2D) && !PsychIsGLES(source)) {
                glTexParameteri(texturetarget, GL_TEXTURE_BASE_LEVEL, 0);
                glTexParameteri(texturetarget, GL_TEXTURE_MAX_LEVEL,  1000);
            }
        }
        else {
            // Negative filterMode: This is mostly meant for fast drawing of blurred (low-pass filtered) textures
            // by selecting a specific integral mip-level:
            glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
            glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            // A negative filterMode means to select a specific mip-level in the
            // mipmap pyramid, according to the filterMode, starting with mip level 0, i.e,
            // full resolution for a value of -1, then level 1 aka half-resolution for a value
            // of -2 etc.:
            if ((texturetarget == GL_TEXTURE_2D) && !PsychIsGLES(source)) {
                glTexParameteri(texturetarget, GL_TEXTURE_BASE_LEVEL, (-1 * filterMode) - 1);
                glTexParameteri(texturetarget, GL_TEXTURE_MAX_LEVEL,  (-1 * filterMode) - 1);
            }
        }

        // Optional texture lookup shader set up (in Screen('MakeTexture') or due to disabled color clamping...)
        if (source->textureLookupShader > 0) {
            shader = source->textureLookupShader;
            if (0 == PsychSetShader(target, shader)) PsychErrorExitMsg(PsychError_user, "Tried to use a texture lookup shader, but your hardware doesn't support GLSL shaders.");
        }
    }

    // Any automatic shader assigned yet?
    if (shader > 0) {
        // In case our texture (filter)/(lookup) shader also requests/defines a 'modulateColor'
        // attribute in its vertex shader part, this attribute is assigned the
        // unclamped RGBA 'modulateColor' after normalization via the colorrange
        // value of Screen('ColorRange'), or the unclamped globalAlpha value:
        if ((attrib = glGetAttribLocationARB(shader, "modulateColor")) >= 0) {
            if(globalAlpha == DBL_MAX) {
                // globalAlpha disabled: Pass the 'modulateColor' vector:
                glVertexAttrib4dvARB(attrib, target->currentColor);
            }
            else {
                // modulateColor disabled: Pass (1,1,1) as RGB color and globalAlpha as alpha:
                glVertexAttrib4fARB(attrib, 1.0, 1.0, 1.0, (GLfloat) globalAlpha);
            }
        }
    }

    // Setup texture wrap-mode: We usually default to clamping - the best we can do
    // for the rectangle textures we usually use. Special case is the intentional
    // use of power-of-two textures with a real power-of-two size. In that case we
    // enable wrapping mode to allow for scrolling effects -- useful for drifting
    // gratings.
    if (texturetarget==GL_TEXTURE_2D && tWidth==sourceWidth && tHeight==sourceHeight) {
        // Special case: Scrollable real power-of-two textures. Enable wrapping.
        glTexParameteri(texturetarget, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(texturetarget, GL_TEXTURE_WRAP_T, GL_REPEAT);
    }
    else {
        // Default: Clamp to edge.
        glTexParameteri(texturetarget, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(texturetarget, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    }

    // We use GL_MODULATE texture application mode together with the special rectangle color
    // (1,1,1,globalAlpha) -- This way, the alpha blending value is the product of the alpha-
    // value of each texel and the globalAlpha value. --> Can apply global alpha value for
    // global blending without need for a texture alpha-channel...
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    // A globalAlpha of DBL_MAX means: Don't set vertex color here, higher-level code
    // has done it already. Used in SCREENDrawTexture for a global override color...
    if (globalAlpha != DBL_MAX) {
        if (PsychIsGLClassic(source)) {
            glColor4f(1, 1, 1, (GLfloat) globalAlpha);
        }
        else PsychGLColor4f(source, 1, 1, 1, (GLfloat) globalAlpha);
    }

    // Apply a rotation transform for rotated drawing, either to modelview-,
    // or texture matrix.
    if ((rotationAngle != 0.0) && !(source->specialflags & kPsychDontDoRotation)) {
        if (!(source->specialflags & kPsychUseTextureMatrixForRotation)) {
            // Standard case: Transform quad -> Modelview matrix.
            glMatrixMode(GL_MODELVIEW);
            transX=(targetRect[kPsychRight] + targetRect[kPsychLeft]) * 0.5;
            transY=(targetRect[kPsychTop] + targetRect[kPsychBottom]) * 0.5;
        }
        else {
            // Transform texture coordinates -> Texture matrix.
            glMatrixMode(GL_TEXTURE);
            transX=(sourceX + sourceXEnd) * 0.5;
            transY=(sourceY + sourceYEnd) * 0.5;
        }

        glPushMatrix();
        glTranslatef((float) +transX, (float) +transY, 0);
        glRotatef((float) rotationAngle, 0, 0, 1);
        glTranslatef((float) -transX, (float) -transY, 0);
        // Rotation transform ready...
    }

    // Support for basic shading during texture blitting: Useful for very simple
    // single-pass isotropic image processing and for procedural texture mapping:
    if (source->textureFilterShader < 0) {
        // User supplied texture shader for either some on-the-fly texture image processing,
        // or for procedural texture shading/on-the-fly texture synthesis. These can be
        // assigned in Screen('MakeTexture') for procedural texture shading or via a
        // optional shader handle to Screen('DrawTexture');
        shader = -1 * source->textureFilterShader;
        if (0 == PsychSetShader(target, shader)) PsychErrorExitMsg(PsychError_user, "Tried to use a user defined texture shader or procedural texture, but your hardware doesn't support GLSL shaders.");

        // Parameter transfer for advanced procedural shading:
        // We encode all parameters about the blit operation into additional
        // vertex attributes so a complex shader can derive useful information.

        // 'srcRect' parameter: The glTexCoord() calls below encode texture coordinates
        // - and thereby the corners of 'srcRect' - into each vertex, however this
        // info gets potentially transformed by the texture matrix, also each vertex
        // only sees one corner of the srcRect: Therefore we encode srcrect = [left top right bottom]
        // on demand:
        if ((attrib = glGetAttribLocationARB(shader, "srcRect")) >= 0) glVertexAttrib4fARB(attrib, (GLfloat) sourceRect[kPsychLeft], (GLfloat) sourceRect[kPsychTop], (GLfloat) sourceRect[kPsychRight], (GLfloat) sourceRect[kPsychBottom]);

        // 'dstRect' parameter: The glVertex() calls below encode target pixel coordinates
        // - and thereby the corners of 'dstRect' - into each vertex, however this
        // info gets potentially transformed by the modelview/proj. matrix, also each vertex
        // only sees one corner of the dstRect: Therefore we encode dstrect = [left top right bottom]
        // on demand:
        if ((attrib = glGetAttribLocationARB(shader, "dstRect")) >= 0) glVertexAttrib4fARB(attrib, (GLfloat) targetRect[kPsychLeft], (GLfloat) targetRect[kPsychTop], (GLfloat) targetRect[kPsychRight], (GLfloat) targetRect[kPsychBottom]);

        // 'sizeAngleFilterMode' - if requested - encodes texture width in .x component, height in .y
        // requested rotationAngle in .z and the 'filterMode' flags in .w:
        if ((attrib = glGetAttribLocationARB(shader, "sizeAngleFilterMode")) >= 0) glVertexAttrib4fARB(attrib, (GLfloat) sourceWidth, (GLfloat) sourceHeight, (GLfloat) rotationAngle, (GLfloat) filterMode);

        // 'modulateColor' - if requested - encodes the RGBA 'modulateColor' after normalization
        // via the colorrange value of Screen('ColorRange').
        if ((attrib = glGetAttribLocationARB(shader, "modulateColor")) >= 0) {
            if(globalAlpha == DBL_MAX) {
                // globalAlpha disabled: Pass the 'modulateColor' vector:
                glVertexAttrib4dvARB(attrib, target->currentColor);
            }
            else {
                // modulateColor disabled: Pass (1,1,1) as RGB color and globalAlpha as alpha:
                glVertexAttrib4fARB(attrib, 1.0, 1.0, 1.0, (GLfloat) globalAlpha);
            }
        }

        // 'auxParameters0' is the first for components (rows) of the 'auxParameters' argument
        // of Screen('DrawTexture(s)') - if such an argument was spec'd:
        if (target->auxShaderParams) {
            if ((target->auxShaderParamsCount >=4) && ((attrib = glGetAttribLocationARB(shader, "auxParameters0")) >= 0)) glVertexAttrib4dvARB(attrib, target->auxShaderParams);
            if ((target->auxShaderParamsCount >=8) && ((attrib = glGetAttribLocationARB(shader, "auxParameters1")) >= 0)) glVertexAttrib4dvARB(attrib, &(target->auxShaderParams[4]));
            if ((target->auxShaderParamsCount >=12) && ((attrib = glGetAttribLocationARB(shader, "auxParameters2")) >= 0)) glVertexAttrib4dvARB(attrib, &(target->auxShaderParams[8]));
            if ((target->auxShaderParamsCount >=16) && ((attrib = glGetAttribLocationARB(shader, "auxParameters3")) >= 0)) glVertexAttrib4dvARB(attrib, &(target->auxShaderParams[12]));
            if ((target->auxShaderParamsCount >=20) && ((attrib = glGetAttribLocationARB(shader, "auxParameters4")) >= 0)) glVertexAttrib4dvARB(attrib, &(target->auxShaderParams[16]));
            if ((target->auxShaderParamsCount >=24) && ((attrib = glGetAttribLocationARB(shader, "auxParameters5")) >= 0)) glVertexAttrib4dvARB(attrib, &(target->auxShaderParams[20]));
            if ((target->auxShaderParamsCount >=28) && ((attrib = glGetAttribLocationARB(shader, "auxParameters6")) >= 0)) glVertexAttrib4dvARB(attrib, &(target->auxShaderParams[24]));
            if ((target->auxShaderParamsCount >=32) && ((attrib = glGetAttribLocationARB(shader, "auxParameters7")) >= 0)) glVertexAttrib4dvARB(attrib, &(target->auxShaderParams[28]));
        }
    }

    // Test for standard case: No shader requested for this texture. In that case we make sure that really no shader is bound.
    if (shader == 0) {
        PsychSetShader(target, 0);
    }
    else {
        #if PSYCH_SYSTEM == PSYCH_OSX
            // On OS-X we can query the OS if the bound shader is running on the GPU or if it is running in emulation mode on the CPU.
            // This is an expensive operation - it triggers OpenGL internal state revalidation. Only use for debugging and testing!
            if (PsychPrefStateGet_Verbosity() > 10) {
                GLint vsgpu=0, fsgpu=0;
                CGLGetParameter(CGLGetCurrentContext(), kCGLCPGPUVertexProcessing, &vsgpu);
                CGLGetParameter(CGLGetCurrentContext(), kCGLCPGPUFragmentProcessing, &fsgpu);
                printf("PTB-DEBUG: In Screen('DrawTexture') aka PsychBlitTextureToDisplay():  GPU shading state: Vertex processing on %s : Fragment processing on %s.\n", (vsgpu) ? "GPU" : "CPU!!", (fsgpu) ? "GPU" : "CPU!!");
            }
        #endif
    }

    if (PsychIsGLClassic(source)) {
        // Classic OpenGL-1/2:

        // matchups for inverted Y coordinate frame (which is inverted?)
        // MK: Texture coordinate assignments have been changed.
        // Explanation: Matlab stores matrices in column-major order, but OpenGL requires
        // textures in row-major order. The old implementation of AWI performed row-column
        // swapping in MakeTexture via C-Code on the CPU. This makes copy-loop implementation
        // complex and creates "Cash trashing" effects on the processor. --> slow MakeTexture performance.
        // Now we store the textures as provided by Matlab, simplifying MakeTexture's implementation,
        // and let the Graphics hardware do the job of "swapping" during rendering, by drawing the texture
        // in some rotated and mirrored order. This is way faster, as the GPU is optimized for such things...
        glBegin(GL_QUADS);
        // Coordinate assignments depend on internal texture orientation...
        if (source->textureOrientation == 2 ||
            source->textureOrientation == 3 || source->textureOrientation == 4) {
            // Use "normal" coordinate assignments, so that the rotation == 0 deg. case
            // is the fastest case --> Most common orientation has highest performance.
            //lower left
            glTexCoord2f((GLfloat)sourceX, (GLfloat)sourceYEnd);
            glVertex2f((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychTop]));        //upper left vertex in window

            //upper left
            glTexCoord2f((GLfloat)sourceX, (GLfloat)sourceY);
            glVertex2f((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychBottom]));     //lower left vertex in window

            //upper right
            glTexCoord2f((GLfloat)sourceXEnd, (GLfloat)sourceY);
            glVertex2f((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychBottom]) );   //lower right  vertex in window

            //lower right
            glTexCoord2f((GLfloat)sourceXEnd, (GLfloat)sourceYEnd);
            glVertex2f((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychTop]));       //upper right in window
        }
        else {
            // Use swapped texture coordinates....
            //lower left
            glTexCoord2f((GLfloat)sourceX, (GLfloat)sourceY);                                       //lower left vertex in  window
            glVertex2f((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychTop]));        //upper left vertex in window

            //upper left
            glTexCoord2f((GLfloat)sourceXEnd, (GLfloat)sourceY);                                    //upper left vertex in texture
            glVertex2f((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychBottom]));     //lower left vertex in window

            //upper right
            glTexCoord2f((GLfloat)sourceXEnd, (GLfloat)sourceYEnd);                                 //upper right vertex in texture
            glVertex2f((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychBottom]) );   //lower right  vertex in window

            //lower right
            glTexCoord2f((GLfloat)sourceX, (GLfloat)sourceYEnd);                                    //lower right in texture
            glVertex2f((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychTop]));       //upper right in window
        }

        glEnd();
    }
    else {
        // Non-classic path (OpenGL-3/4/... and OpenGL-ES). No GL_QUADS and no immediate mode,
        // need to do it with triangles and vertex arrays:

        // Define "windowRecord", so our convenience macros work:
        PsychWindowRecordType *windowRecord = source;

        GLBEGIN(GL_TRIANGLE_STRIP);
        // Coordinate assignments depend on internal texture orientation...
        if (source->textureOrientation == 2 ||
            source->textureOrientation == 3 || source->textureOrientation == 4) {
            // Use "normal" coordinate assignments, so that the rotation == 0 deg. case
            // is the fastest case --> Most common orientation has highest performance.
            //upper left
            GLTEXCOORD2f((GLfloat)sourceX, (GLfloat)sourceY);
            GLVERTEX2f((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychBottom]));     //lower left vertex in window

            //lower left
            GLTEXCOORD2f((GLfloat)sourceX, (GLfloat)sourceYEnd);
            GLVERTEX2f((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychTop]));        //upper left vertex in window

            //upper right
            GLTEXCOORD2f((GLfloat)sourceXEnd, (GLfloat)sourceY);
            GLVERTEX2f((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychBottom]) );   //lower right  vertex in window

            //lower right
            GLTEXCOORD2f((GLfloat)sourceXEnd, (GLfloat)sourceYEnd);
            GLVERTEX2f((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychTop]));       //upper right in window
        }
        else {
            // Use swapped texture coordinates....

            //upper left
            GLTEXCOORD2f((GLfloat)sourceXEnd, (GLfloat)sourceY);                                    //upper left vertex in texture
            GLVERTEX2f((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychBottom]));     //lower left vertex in window

            //lower left
            GLTEXCOORD2f((GLfloat)sourceX, (GLfloat)sourceY);                                       //lower left vertex in  window
            GLVERTEX2f((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychTop]));        //upper left vertex in window

            //upper right
            GLTEXCOORD2f((GLfloat)sourceXEnd, (GLfloat)sourceYEnd);                                 //upper right vertex in texture
            GLVERTEX2f((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychBottom]) );   //lower right  vertex in window

            //lower right
            GLTEXCOORD2f((GLfloat)sourceX, (GLfloat)sourceYEnd);                                    //lower right in texture
            GLVERTEX2f((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychTop]));       //upper right in window
        }

        GLEND();
    }

    // Undo rotation transform, if any...
    if ((rotationAngle != 0.0) && !(source->specialflags & kPsychDontDoRotation)) {
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
    }

    // Only disable texture mapping if we actually enabled it.
    if (source->textureNumber > 0) {
        // Reset filters to nearest: This is important in case this texture
        // is used as color buffer attachment of a FBO, because using the
        // FBO would fail in puzzling ways if filtermode!=GL_NEAREST.
        glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        // Don't restrict mipmap-levels for sampling, reset to initial system defaults:
        if  ((texturetarget == GL_TEXTURE_2D) && !PsychIsGLES(source)) {
            glTexParameteri(texturetarget, GL_TEXTURE_BASE_LEVEL, 0);
            glTexParameteri(texturetarget, GL_TEXTURE_MAX_LEVEL,  1000);
        }

        // Unbind texture:
        glBindTexture(texturetarget, 0);
        glDisable(texturetarget);
    }

    /* Dead and disabled: Left here for documentation...
     *    if ((filterMode > 0 && source->textureFilterShader > 0) || (source->textureFilterShader < 0)) {
     *        // Filtershader was used. Unbind it:
     *        // MK: Not really: We used to do this, but have a new policy starting 9th March 2008. All rendering
     *        // functions have to use PsychSetShader() or glUseProgram() calls to define their wanted shader
     *        // binding. Therefore no need to revert to a zero binding at end of a routine. This avoids lots
     *        // of redundant shader switches in batch drawing routines, e.g., if this routine is called via
     *        // Screen('DrawTextures'), especially with procedural textures bound.
     *        // glUseProgram(0);
     *    }
     */

    // Finished!
    return;
}

/* PsychGetTextureTarget
 * Returns GLenum with the texture target used for all PTB operations.
 * This way, external code can bind the correct target for given hardware.
 */
GLenum PsychGetTextureTarget(PsychWindowRecordType *win)
{
    (void) win;

    // If texturetaget field for this texture isn't yet initialized, then
    // init it now from our global setting:
    if (win->texturetarget == 0) {
        // Setup texture-target if not already done:
        PsychDetectTextureTarget(win);

        win->texturetarget = texturetarget;
    }

    // Return texturetarget for this window:
    return(win->texturetarget);
}

void PsychMapTexCoord(PsychWindowRecordType *tex, double* tx, double* ty)
{
    GLdouble sourceWidth, sourceHeight, tWidth, tHeight;
    GLdouble sourceX, sourceY;
    GLenum texturetarget;

    if (!PsychIsTexture(tex)) {
        PsychErrorExitMsg(PsychError_user, "Tried to map (x,y) texel position to texture coordinate for something else than a texture!");
    }

    if (tx==NULL || ty==NULL) PsychErrorExitMsg(PsychError_internal, "NULL-Ptr passed as tx or ty into PsychMapTexCoord()!!!");

    // Perform mapping: This mostly duplicates code in PsychBlitTextureToDisplay():

    // Setup texture-target if not already done:
    PsychDetectTextureTarget(tex);

    // Assign proper texturetarget for mapping:
    texturetarget = PsychGetTextureTarget(tex);

    // Basic mapping for rectangular textures:
    // 0 == Transposed as from Matlab image array. 2 == Offscreen window in normal orientation.
    if (tex->textureOrientation == 2) {
        sourceHeight=PsychGetHeightFromRect(tex->rect);
        sourceWidth=PsychGetWidthFromRect(tex->rect);

        sourceX=*tx;
        sourceY=sourceHeight - *ty;
    }
    else {
        // Transposed texture from Matlab:
        sourceHeight=PsychGetWidthFromRect(tex->rect);
        sourceWidth=PsychGetHeightFromRect(tex->rect);
        sourceX=*tx;
        sourceY=*ty;
    }

    // Override for special case: Upside-Down texture
    if (tex->textureOrientation == 3) {
        sourceHeight=PsychGetHeightFromRect(tex->rect);
        sourceWidth=PsychGetWidthFromRect(tex->rect);
        sourceX=*tx;
        sourceY=*ty;
    }

    // Special case handling for GL_TEXTURE_2D textures. We need to map the
    // absolute texture coordinates (in pixels) to the interval 0.0 - 1.0 where
    // 1.0 == full extent of power of two texture...
    if (texturetarget==GL_TEXTURE_2D) {
        // NPOT supported?
        if (!(tex->gfxcaps & kPsychGfxCapNPOTTex)) {
            // No: Find size of real underlying texture (smallest power of two which is
            // greater than or equal to the image size:
            tWidth=1;
            while (tWidth < sourceWidth) tWidth*=2;
            tHeight=1;
            while (tHeight < sourceHeight) tHeight*=2;
        }
        else {
            // Yes:
            tWidth = sourceWidth;
            tHeight = sourceHeight;
        }

        // Remap texcoords into 0-1 subrange: We subtract 0.5 pixel-units before
        // mapping to accomodate for roundoff-error in the power-of-two gfx
        // hardware...
        // For a good intro into the issue of texture border seams, due to interpolation
        // problems at texture borders, see:
        // http://home.planet.nl/~monstrous/skybox.html
        sourceX-=0.5f;
        sourceY-=0.5f;

        // Remap:
        sourceX=sourceX / tWidth;
        sourceY=sourceY / tHeight;
    }

    // Return mapped coords:
    *tx = sourceX;
    *ty = sourceY;

    // Done.
    return;
}

static float transX, transY;
static float crt, srt;
static unsigned int useXForm = 0;

static inline void glVertexXform(GLfloat x, GLfloat y)
{
    GLfloat xo, yo;

    if (useXForm == 1) {
        x = x - transX;
        y = y - transY;

        xo = crt * x - srt * y;
        yo = srt * x + crt * y;

        xo = xo + transX;
        yo = yo + transY;
    }
    else {
        xo = x;
        yo = y;
    }

    glVertex2f(xo, yo);
}

static inline void glTexCoordXform(GLfloat x, GLfloat y)
{
    GLfloat xo, yo;

    if (useXForm == 2) {
        x = x - transX;
        y = y - transY;

        xo = crt * x - srt * y;
        yo = srt * x + crt * y;

        xo = xo + transX;
        yo = yo + transY;
    }
    else {
        xo = x;
        yo = y;
    }

    glTexCoord2f(xo, yo);
}

void PsychBatchBlitTexturesToDisplay(unsigned int opMode, unsigned int count, PsychWindowRecordType *source, PsychWindowRecordType *target, double *sourceRect, double *targetRect,
                                     double rotationAngle, int filterMode, double globalAlpha)
{
    static unsigned int index = 0;
    static float *vertices, *colors, *texcoords;

    static GLint attribs[11];
    static GLenum texturetarget;
    static GLint textureNumber = -1;
    static GLint mattrib = -1;
    static GLuint shader = 0;
    static int tWidth = 0, tHeight = 0;
    static double oldRotationAngle;
    static GLdouble sourceWidth, sourceHeight;
    GLdouble sourceX, sourceY, sourceXEnd, sourceYEnd;

    if (opMode == 0) {
        // Start new batch:
        if (vertices || colors || texcoords)
            PsychErrorExitMsg(PsychError_internal, "Non-NULL arrays at start of batch for opMode 0!\n");

        index = 0;
        //vertices = malloc(count * 2 *sizeof(float));
        //colors = malloc(count * 4 * sizeof(float));
        //texcoords = malloc(count * 2 * sizeof(float));

        // Enable targets framebuffer as current drawingtarget, except if this is a
        // blit operation from a window into itself and the imaging pipe is on:
        if ((source != target) || (target->imagingMode==0)) {
            PsychSetDrawingTarget(target);
        }
        else {
            // Activate rendering context of target window without changing drawing target:
            PsychSetGLContext(target);
        }

        // Disable Transform:
        useXForm = 0;
        oldRotationAngle = 0;

        shader = 0;
        return;
    }

    if (opMode == 1) {
        // Finalize this batch:

        // DRAW DRAW DRAW DRAW!

        glEnd();

        // Disable Transform:
        useXForm = 0;
        oldRotationAngle = 0;

        // Only disable texture mapping if we actually enabled it.
        if (textureNumber > 0) {
            // Reset filters to nearest: This is important in case this texture
            // is used as color buffer attachment of a FBO, because using the
            // FBO would fail in puzzling ways if filtermode!=GL_NEAREST.
            glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

            // Don't restrict mipmap-levels for sampling, reset to initial system defaults:
            if  ((texturetarget == GL_TEXTURE_2D) && !PsychIsGLES(source)) {
                glTexParameteri(texturetarget, GL_TEXTURE_BASE_LEVEL, 0);
                glTexParameteri(texturetarget, GL_TEXTURE_MAX_LEVEL,  1000);
            }

            // Unbind texture:
            glBindTexture(texturetarget, 0);
            glDisable(texturetarget);
        }

        //free(vertices);
        vertices = NULL;
        //free(colors);
        colors = NULL;
        //free(texcoords);
        texcoords = NULL;

        return;
    }

    // opMode 2: Add a new texture to buffers:

    // First element to draw? Need some more setup from information derived from
    // first item:
    if (index == 0) {
        // Setup texture-target if not already done:
        PsychDetectTextureTarget(target);

        // Query target for this specific texture:
        texturetarget = PsychGetTextureTarget(source);

        // 0 == Transposed as from Matlab image array. 2 == Offscreen window in normal orientation.
        if (source->textureOrientation == 2) {
            sourceHeight=PsychGetHeightFromRect(source->rect);
            sourceWidth=PsychGetWidthFromRect(source->rect);
        }
        else {
            sourceHeight=PsychGetWidthFromRect(source->rect);
            sourceWidth=PsychGetHeightFromRect(source->rect);
        }

        // Overrides for special cases: Upside-down texture.
        if (source->textureOrientation == 3) {
            sourceHeight=PsychGetHeightFromRect(source->rect);
            sourceWidth=PsychGetWidthFromRect(source->rect);
        }

        // This case can happen with some QT movies, they are upside down in an unusual way:
        if (source->textureOrientation == 4) {
            sourceHeight=PsychGetHeightFromRect(source->rect);
            sourceWidth=PsychGetWidthFromRect(source->rect);
        }

        // Special case handling for GL_TEXTURE_2D textures. We need to map the
        // absolute texture coordinates (in pixels) to the interval 0.0 - 1.0 where
        // 1.0 == full extent of power of two texture...
        if (texturetarget == GL_TEXTURE_2D) {
            // NPOT supported?
            if (!(source->gfxcaps & kPsychGfxCapNPOTTex)) {
                // No: Find size of real underlying texture (smallest power of two which is
                // greater than or equal to the image size:
                tWidth=1;
                while (tWidth < (int) sourceWidth) tWidth*=2;
                tHeight=1;
                while (tHeight < (int) sourceHeight) tHeight*=2;
            }
            else {
                // Yes:
                tWidth = (int) sourceWidth;
                tHeight = (int) sourceHeight;
            }
        }

        // Only enable actual texture hardware if a real texture is provided.
        // In the case of no real texture, we don't bind a real texture, don't
        // enable texture mapping and just blit the quad, with interpolated
        // texture coordinates set up for purely procedural shading.
        glDisable(GL_TEXTURE_2D);
        if (source->textureNumber > 0) {
            glEnable(texturetarget);
            glBindTexture(texturetarget, source->textureNumber);
        }

        // Use of OpenGL mip-mapping requested? And automatic mipmap generation wanted - aka not forbidden?
        if (!(source->specialflags & kPsychDontAutoGenMipMaps) && (filterMode < 0 || filterMode > 1)) {
            // Yes: Automatically build a mip-map pyramid.
            if (texturetarget != GL_TEXTURE_2D) PsychErrorExitMsg(PsychError_user, "You asked me to use mip-mapped texture filtering on a texture that is not of GL_TEXTURE_2D type! Unsupported.");

            if (NULL != glGenerateMipmapEXT) {
                GLint mipattrib = 0;

                // Select highest quality downsampling method:
                glHint(GL_GENERATE_MIPMAP_HINT, GL_NICEST);

                // Is automatic mipmap generation already enabled for this texture?
                glGetTexParameteriv(texturetarget, GL_GENERATE_MIPMAP, &mipattrib);

                // Need to manually trigger regen if automatic mode not yet enabled, or
                // if the "dirty" flag is set, because some render-to-texture activity has
                // happened since last drawing this texture, which is not covered by the auto-update:
                if (source->needsViewportSetup || !mipattrib) {
                    // No: We trigger hardware-accelerated mipmap generation manually for this draw call:
                    glGenerateMipmapEXT(texturetarget);

                    // Enable automatic mipmap generation for future updates to this texture object. This
                    // will automatically trigger regen if new image content is uploaded into this texture
                    // object:
                    glTexParameteri(texturetarget, GL_GENERATE_MIPMAP, GL_TRUE);

                    // Clear "dirty" flag:
                    source->needsViewportSetup = FALSE;
                }
            }
            else if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-WARNING: Was asked to draw a texture with mip-mapping, but automatic mipmap generation unsupported by this system! Check your stimulus!\n");
            }
        }

        // Linear filtering on non-capable hardware via shader emulation?
        if ((filterMode != 0) && (source->textureFilterShader > 0)) {
            // Yes. Bind the shader:
            shader = source->textureFilterShader;
            if (0 == PsychSetShader(target, shader)) PsychErrorExitMsg(PsychError_user, "Tried to use a bilinear texture filter shader, but your hardware doesn't support GLSL shaders.");

            if (filterMode < 0 || filterMode > 1) {
                // Some mip-mapped filtermode. Choose nearest neighbour sampling within mipmap levels, so shader can decide about sample locations itself.
                // In filterMode 2 choose the nearest mipmap, in others interpolate linearly between two nearest mipmap levels:
                glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, (filterMode == 2) ? GL_NEAREST_MIPMAP_NEAREST : GL_NEAREST_MIPMAP_LINEAR);
                glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            }
            else {
                // No mip-mapping: Switch hardware samplers into nearest neighbour mode so we don't get any interference:
                glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
                glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            }

            // Don't restrict mipmap-levels for sampling, reset to initial system defaults:
            // This makes even sense for negative filterMode arguments, because the filterMode
            // parameter is passed as an attribute to the filtershader, so the shader itself can
            // decide how to implement a specific blur level on its own, unrestricted by us:
            if ((texturetarget == GL_TEXTURE_2D) && !PsychIsGLES(source)) {
                glTexParameteri(texturetarget, GL_TEXTURE_BASE_LEVEL, 0);
                glTexParameteri(texturetarget, GL_TEXTURE_MAX_LEVEL,  1000);
            }
        }
        else {
            // Standard hardware texture sampling/filtering: Select filter-mode for texturing:
            if (filterMode >= 0) {
                // Select specific hardware sampling strategy:
                switch (filterMode) {
                    case 0: // Nearest-Neighbour filtering:
                        glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
                        glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                        break;

                    case 1: // Bilinear filtering:
                        glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
                        glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                        break;

                    case 2: // Linear filtering with nearest neighbour mipmapping:
                        glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
                        glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                        break;

                    case 3: // Linear filtering with linear mipmapping --> This is full trilinear filtering:
                        glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
                        glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                        break;

                    case 4: // Nearest-Neighbour filtering with nearest neighbour mipmapping:
                        glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
                        glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                        break;

                    case 5: // Nearest-Neighbour filtering with linear mipmapping:
                        glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_LINEAR);
                        glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                        break;
                }

                // Don't restrict mipmap-levels for sampling, reset to initial system defaults:
                if ((texturetarget == GL_TEXTURE_2D) && !PsychIsGLES(source)) {
                    glTexParameteri(texturetarget, GL_TEXTURE_BASE_LEVEL, 0);
                    glTexParameteri(texturetarget, GL_TEXTURE_MAX_LEVEL,  1000);
                }
            }
            else {
                // Negative filterMode: This is mostly meant for fast drawing of blurred (low-pass filtered) textures
                // by selecting a specific integral mip-level:
                glTexParameteri(texturetarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
                glTexParameteri(texturetarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

                // A negative filterMode means to select a specific mip-level in the
                // mipmap pyramid, according to the filterMode, starting with mip level 0, i.e,
                // full resolution for a value of -1, then level 1 aka half-resolution for a value
                // of -2 etc.:
                if ((texturetarget == GL_TEXTURE_2D) && !PsychIsGLES(source)) {
                    glTexParameteri(texturetarget, GL_TEXTURE_BASE_LEVEL, (-1 * filterMode) - 1);
                    glTexParameteri(texturetarget, GL_TEXTURE_MAX_LEVEL,  (-1 * filterMode) - 1);
                }
            }

            // Optional texture lookup shader set up (in Screen('MakeTexture') or due to disabled color clamping...)
            if (source->textureLookupShader > 0) {
                shader = source->textureLookupShader;
                if (0 == PsychSetShader(target, shader)) PsychErrorExitMsg(PsychError_user, "Tried to use a texture lookup shader, but your hardware doesn't support GLSL shaders.");
            }
        }

        // Support for basic shading during texture blitting: Useful for very simple
        // single-pass isotropic image processing and for procedural texture mapping:
        if (source->textureFilterShader < 0) {
            // User supplied texture shader for either some on-the-fly texture image processing,
            // or for procedural texture shading/on-the-fly texture synthesis. These can be
            // assigned in Screen('MakeTexture') for procedural texture shading or via a
            // optional shader handle to Screen('DrawTexture');
            shader = -1 * source->textureFilterShader;
            if (0 == PsychSetShader(target, shader)) PsychErrorExitMsg(PsychError_user, "Tried to use a user defined texture shader or procedural texture, but your hardware doesn't support GLSL shaders.");

            // Query and cache all vertex attrib locations - Can't do within glBegin/End():
            attribs[0] = glGetAttribLocationARB(shader, "srcRect");
            attribs[1] = glGetAttribLocationARB(shader, "dstRect");
            attribs[2] = glGetAttribLocationARB(shader, "sizeAngleFilterMode");
            attribs[3] = glGetAttribLocationARB(shader, "auxParameters0");
            attribs[4] = glGetAttribLocationARB(shader, "auxParameters1");
            attribs[5] = glGetAttribLocationARB(shader, "auxParameters2");
            attribs[6] = glGetAttribLocationARB(shader, "auxParameters3");
            attribs[7] = glGetAttribLocationARB(shader, "auxParameters4");
            attribs[8] = glGetAttribLocationARB(shader, "auxParameters5");
            attribs[9] = glGetAttribLocationARB(shader, "auxParameters6");
            attribs[10] = glGetAttribLocationARB(shader, "auxParameters7");
        }

        if (shader > 0) {
            // In case our texture (filter)/(lookup) shader also requests/defines a 'modulateColor'
            // attribute in its vertex shader part, this attribute is assigned the
            // unclamped RGBA 'modulateColor' after normalization via the colorrange
            // value of Screen('ColorRange'), or the unclamped globalAlpha value:
            mattrib = glGetAttribLocationARB(shader, "modulateColor");
        }
        else {
            // Not needed:
            mattrib = -1;
        }

        // Setup texture wrap-mode: We usually default to clamping - the best we can do
        // for the rectangle textures we usually use. Special case is the intentional
        // use of power-of-two textures with a real power-of-two size. In that case we
        // enable wrapping mode to allow for scrolling effects -- useful for drifting
        // gratings.
        if (texturetarget==GL_TEXTURE_2D && tWidth==sourceWidth && tHeight==sourceHeight) {
            // Special case: Scrollable real power-of-two textures. Enable wrapping.
            glTexParameteri(texturetarget, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(texturetarget, GL_TEXTURE_WRAP_T, GL_REPEAT);
        }
        else {
            // Default: Clamp to edge.
            glTexParameteri(texturetarget, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(texturetarget, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        }

        // We use GL_MODULATE texture application mode together with the special rectangle color
        // (1,1,1,globalAlpha) -- This way, the alpha blending value is the product of the alpha-
        // value of each texel and the globalAlpha value. --> Can apply global alpha value for
        // global blending without need for a texture alpha-channel...
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

        // Test for standard case: No shader requested for this texture. In that case we make sure that really no shader is bound.
        if (shader == 0) {
            PsychSetShader(target, 0);
        }
        else {
            #if PSYCH_SYSTEM == PSYCH_OSX
            // On OS-X we can query the OS if the bound shader is running on the GPU or if it is running in emulation mode on the CPU.
            // This is an expensive operation - it triggers OpenGL internal state revalidation. Only use for debugging and testing!
            if (PsychPrefStateGet_Verbosity() > 10) {
                GLint vsgpu=0, fsgpu=0;
                CGLGetParameter(CGLGetCurrentContext(), kCGLCPGPUVertexProcessing, &vsgpu);
                CGLGetParameter(CGLGetCurrentContext(), kCGLCPGPUFragmentProcessing, &fsgpu);
                printf("PTB-DEBUG: In Screen('DrawTextures') aka PsychBatchBlitTexturesToDisplay():  GPU shading state: Vertex processing on %s : Fragment processing on %s.\n", (vsgpu) ? "GPU" : "CPU!!", (fsgpu) ? "GPU" : "CPU!!");
            }
            #endif
        }

        textureNumber = source->textureNumber;

        glBegin(GL_QUADS);

        // End of prep for first texture quad.
    }

    // 0 == Transposed as from Matlab image array. 2 == Offscreen window in normal orientation.
    if (source->textureOrientation == 2) {
        sourceX=sourceRect[kPsychLeft];
        sourceY=sourceHeight - sourceRect[kPsychBottom];
        sourceXEnd=sourceRect[kPsychRight];
        sourceYEnd=sourceHeight - sourceRect[kPsychTop];
    }
    else {
        sourceX=sourceRect[kPsychTop];
        sourceY=sourceRect[kPsychLeft];
        sourceXEnd=sourceRect[kPsychBottom];
        sourceYEnd=sourceRect[kPsychRight];
    }

    // Overrides for special cases: Upside-down texture.
    if (source->textureOrientation == 3) {
        sourceX=sourceRect[kPsychLeft];
        sourceY=sourceRect[kPsychBottom];
        sourceXEnd=sourceRect[kPsychRight];
        sourceYEnd=sourceRect[kPsychTop];
    }

    // This case can happen with some QT movies, they are upside down in an unusual way:
    if (source->textureOrientation == 4) {
        sourceX=sourceRect[kPsychLeft];
        sourceY=sourceHeight - sourceRect[kPsychBottom];
        sourceXEnd=sourceRect[kPsychRight];
        sourceYEnd=sourceHeight - sourceRect[kPsychTop];
    }

    // Special case handling for GL_TEXTURE_2D textures. We need to map the
    // absolute texture coordinates (in pixels) to the interval 0.0 - 1.0 where
    // 1.0 == full extent of power of two texture...
    if (texturetarget == GL_TEXTURE_2D) {
        // Remap texcoords into 0-1 subrange: We subtract 0.5 pixel-units before
        // mapping to accomodate for roundoff-error in the power-of-two gfx
        // hardware...
        // For a good intro into the issue of texture border seams, due to interpolation
        // problems at texture borders, see:
        // http://home.planet.nl/~monstrous/skybox.html
        sourceXEnd-=0.5f;
        sourceYEnd-=0.5f;

        // Remap:
        sourceX=sourceX / (double) tWidth;
        sourceXEnd=sourceXEnd / (double) tWidth;
        sourceY=sourceY / (double) tHeight;
        sourceYEnd=sourceYEnd / (double) tHeight;
    }

    // Any automatic shader assigned yet?
    if (shader > 0 && mattrib >= 0) {
        if (globalAlpha == DBL_MAX) {
            // globalAlpha disabled: Pass the 'modulateColor' vector:
            glVertexAttrib4dvARB(mattrib, target->currentColor);
        }
        else {
            // modulateColor disabled: Pass (1,1,1) as RGB color and globalAlpha as alpha:
            glVertexAttrib4fARB(mattrib, 1.0, 1.0, 1.0, (GLfloat) globalAlpha);
        }
    }

    // Fixed function pipeline color assignment:
    if (globalAlpha == DBL_MAX) {
        glColor4dv(target->currentColor);
    }
    else {
        glColor4f(1, 1, 1, (GLfloat) globalAlpha);
    }

    if ((rotationAngle != 0) && !(source->specialflags & kPsychDontDoRotation)) {
        // Apply a rotation transform for rotated drawing, either to modelview-,
        // or texture matrix.
        if (!(source->specialflags & kPsychUseTextureMatrixForRotation)) {
            // Standard case: Transform quad vertex coordinates:
            useXForm = 1;
            transX = (float) ((targetRect[kPsychRight] + targetRect[kPsychLeft]) * 0.5);
            transY = (float) ((targetRect[kPsychTop] + targetRect[kPsychBottom]) * 0.5);
        }
        else {
            // Transform texture coordinates:
            useXForm = 2;
            transX = (float) ((sourceX + sourceXEnd) * 0.5);
            transY = (float) ((sourceY + sourceYEnd) * 0.5);
        }

        // Avoid redundant updates of rotation transform 2D "matrix":
        if (rotationAngle != oldRotationAngle) {
            double rotAngleRad;

            // New angle. Update rotation matrix coefficients:
            rotAngleRad = rotationAngle * M_PI / 180.0;
            crt = (float) cos(rotAngleRad);
            srt = (float) sin(rotAngleRad);
        }
    }
    else {
        // Disable transform:
        useXForm = 0;
    }

    oldRotationAngle = rotationAngle;

    // Support for basic shading during texture blitting: Useful for very simple
    // single-pass isotropic image processing and for procedural texture mapping:
    if (source->textureFilterShader < 0) {
        // User supplied texture shader for either some on-the-fly texture image processing,
        // or for procedural texture shading/on-the-fly texture synthesis already assigned
        // at first quad submission above.
        //
        // Now need parameter transfer for advanced procedural shading:
        // We encode all parameters about the blit operation into additional
        // vertex attributes so a complex shader can derive useful information.

        // 'srcRect' parameter: The glTexCoord() calls below encode texture coordinates
        // - and thereby the corners of 'srcRect' - into each vertex, however this
        // info gets potentially transformed by the texture matrix, also each vertex
        // only sees one corner of the srcRect: Therefore we encode srcrect = [left top right bottom]
        // on demand:
        if (attribs[0] >= 0) glVertexAttrib4fARB(attribs[0], (GLfloat) sourceRect[kPsychLeft], (GLfloat) sourceRect[kPsychTop], (GLfloat) sourceRect[kPsychRight], (GLfloat) sourceRect[kPsychBottom]);

        // 'dstRect' parameter: The glVertex() calls below encode target pixel coordinates
        // - and thereby the corners of 'dstRect' - into each vertex, however this
        // info gets potentially transformed by the modelview/proj. matrix, also each vertex
        // only sees one corner of the dstRect: Therefore we encode dstrect = [left top right bottom]
        // on demand:
        if (attribs[1] >= 0) glVertexAttrib4fARB(attribs[1], (GLfloat) targetRect[kPsychLeft], (GLfloat) targetRect[kPsychTop], (GLfloat) targetRect[kPsychRight], (GLfloat) targetRect[kPsychBottom]);

        // 'sizeAngleFilterMode' - if requested - encodes texture width in .x component, height in .y
        // requested rotationAngle in .z and the 'filterMode' flags in .w:
        if (attribs[2] >= 0) glVertexAttrib4fARB(attribs[2], (GLfloat) sourceWidth, (GLfloat) sourceHeight, (GLfloat) rotationAngle, (GLfloat) filterMode);

        // 'auxParameters0' is the first for components (rows) of the 'auxParameters' argument
        // of Screen('DrawTexture(s)') - if such an argument was spec'd:
        if (target->auxShaderParams) {
            if ((target->auxShaderParamsCount >=4) && (attribs[3] >= 0)) glVertexAttrib4dvARB(attribs[3], target->auxShaderParams);
            if ((target->auxShaderParamsCount >=8) && (attribs[4] >= 0)) glVertexAttrib4dvARB(attribs[4], &(target->auxShaderParams[4]));
            if ((target->auxShaderParamsCount >=12) && (attribs[5] >= 0)) glVertexAttrib4dvARB(attribs[5], &(target->auxShaderParams[8]));
            if ((target->auxShaderParamsCount >=16) && (attribs[6] >= 0)) glVertexAttrib4dvARB(attribs[6], &(target->auxShaderParams[12]));
            if ((target->auxShaderParamsCount >=20) && (attribs[7] >= 0)) glVertexAttrib4dvARB(attribs[7], &(target->auxShaderParams[16]));
            if ((target->auxShaderParamsCount >=24) && (attribs[8] >= 0)) glVertexAttrib4dvARB(attribs[8], &(target->auxShaderParams[20]));
            if ((target->auxShaderParamsCount >=28) && (attribs[9] >= 0)) glVertexAttrib4dvARB(attribs[9], &(target->auxShaderParams[24]));
            if ((target->auxShaderParamsCount >=32) && (attribs[10] >= 0)) glVertexAttrib4dvARB(attribs[10], &(target->auxShaderParams[28]));
        }
    }

    // Coordinate assignments depend on internal texture orientation...
    if (source->textureOrientation == 2 ||
        source->textureOrientation == 3 || source->textureOrientation == 4) {
        // Use "normal" coordinate assignments, so that the rotation == 0 deg. case
        // is the fastest case --> Most common orientation has highest performance.
        //lower left
        glTexCoordXform((GLfloat)sourceX, (GLfloat)sourceYEnd);
        glVertexXform((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychTop]));        //upper left vertex in window

        //upper left
        glTexCoordXform((GLfloat)sourceX, (GLfloat)sourceY);
        glVertexXform((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychBottom]));     //lower left vertex in window

        //upper right
        glTexCoordXform((GLfloat)sourceXEnd, (GLfloat)sourceY);
        glVertexXform((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychBottom]) );   //lower right  vertex in window

        //lower right
        glTexCoordXform((GLfloat)sourceXEnd, (GLfloat)sourceYEnd);
        glVertexXform((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychTop]));       //upper right in window
    }
    else {
        // Use swapped texture coordinates....
        //lower left
        glTexCoordXform((GLfloat)sourceX, (GLfloat)sourceY);                                       //lower left vertex in  window
        glVertexXform((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychTop]));        //upper left vertex in window

        //upper left
        glTexCoordXform((GLfloat)sourceXEnd, (GLfloat)sourceY);                                    //upper left vertex in texture
        glVertexXform((GLfloat)(targetRect[kPsychLeft]), (GLfloat)(targetRect[kPsychBottom]));     //lower left vertex in window

        //upper right
        glTexCoordXform((GLfloat)sourceXEnd, (GLfloat)sourceYEnd);                                 //upper right vertex in texture
        glVertexXform((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychBottom]) );   //lower right  vertex in window

        //lower right
        glTexCoordXform((GLfloat)sourceX, (GLfloat)sourceYEnd);                                    //lower right in texture
        glVertexXform((GLfloat)(targetRect[kPsychRight]), (GLfloat)(targetRect[kPsychTop]));       //upper right in window
    }

    index++;

    // Finished!
    return;
}
