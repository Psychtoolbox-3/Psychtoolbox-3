/*
 *
 * moglcore.cc -- MATLAB MEX and GNU/OCTAVE OCT file interface to OpenGL.
 *
 * 08-May-2005 -- created (RFM)
 * 08-Dec-2005 -- reworked into direct interface to gl, glu, and glm functions (RFM)
 * 05-Mar-2006 -- reworked to make inclusion of glm optional (for Psychtoolbox) (MK)
 * 20-Mar-2006 -- Included support for GLEW lib for auto-detection of OpenGL extensions. (MK)
 * 15-Apr-2006 -- Dynamic rebinding of OpenGL-2 core functions to ARB extensions. (MK)
 * 16-Apr-2006 -- Built-in error detection and handling via the debuglevel parameter. (MK)
 * 16-May-2006 -- Implementation of a buffer memory manager for dynamic buffer allocation and such.
 *                This is needed for commands like glFeedbackBuffer() to work properly.
 * 19-Jun-2006 -- Implement support for GNU/Octave (MK).
 * 24-Mar-2011 -- Make 64-bit clean (MK).
 * 27-Mar-2011 -- Remove obsolete and totally bitrotten Octave-2 support (MK).
 * 03-Apr-2011 -- Allow to receive pointers encoded in double's, uint32 or uint64. Adapt dynamically (MK).
 *
 */

#include "mogltypes.h"

/* Build and include support for glm if BUILD_GLM is defined.
   Otherwise, only build OpenGL wrappers.
*/

#ifdef BUILD_GLM
extern int glm_map_count;
extern cmdhandler glm_map[];
#endif

extern int gl_manual_map_count, gl_auto_map_count;
extern cmdhandler gl_manual_map[], gl_auto_map[];

// Pointer to the start of our list of allocated temporary memory buffers:
static size_t* PsychTempMemHead[4] = { NULL , NULL , NULL , NULL };

// Total count of allocated memory in Bytes:
static size_t totalTempMemAllocated[4] = { 0 , 0 , 0 , 0 };

// Flag that signals first real invocation of moglcore:
static int firsttime = 1;

// Debuglevel: Defines if moglcore should check for errors and how to respond to them.
static int debuglevel = 1;

// Dummymode: If set to a value greater than zero, moglcore does not really execute
// commands anymore. 
static int dummymode = 0;

// This flag is > 0 if we are in a section started by glBegin() in that case we are
// not allowed to call glGetError() as that itself is an error, until the corresponding
// glEnd gets executed.
static int glBeginLevel = 0;

// command string
#define CMDLEN 64
char cmd[CMDLEN];

// binary search routine
int binsearch(cmdhandler *map, int mapsize, char *str);

// error handler
void mogl_usageerr();

// Dynamic extension rebinding: This is in an extra file mogl_rebinder.c and defined
// with explicit C linkage. Reason: For the Octave version, we need to compile all
// files (except glew.c) as C++, not C, otherwise won't work. The C++ compiler is
// more picky about function pointer assignments and will abort for hacks like in
// mogl_rebind...  -> Move this function into separate file, compile that file as
// C file, make the compiler happy.

// Need different cases for C++ compilation vs. C compilation: On OS/X we do C
// compiles for Matlab builds, but C++ compiles for Octave builds. On Windows
// (i think) we do pure C builds.
#ifdef PTBOCTAVE
extern "C" void mogl_rebindARBExtensionsToCore(void);
#else
extern void mogl_rebindARBExtensionsToCore(void);
#endif

// Automatic checking and handling for glError's and GLSL errors.
void mogl_checkerrors(const char* cmd, const mxArray *prhs[]);

void mexExitFunction(void)
{
  // Release all memory in bufferlist 1 - The one that usually
  // persists over calls to moglcore.
  PsychFreeAllTempMemory(1);
  
  // Release all memory for persistent GLU tesselator memory list:
  PsychFreeAllTempMemory(2);
  PsychFreeAllTempMemory(3);
  
  firsttime = 1;
}

/* This is the main entry point for Matlab or Octave. It gets called by
   Matlab or Octave, handles first-time initialization, error handling
   and subfunction dispatching.
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Start of dispatcher:
    int i;
    GLenum err;
    
    // FreeGlut must be initialized, otherwise it will emergency abort the whole application.
    // These variables are needed for it:
    int noargs = 1;
    char dummyarg[] = "ptbmoglcore";
    char *dummyargp = &dummyarg[0];

    // see whether there's a string command
    if(nrhs<1 || !mxIsChar(prhs[0])) mogl_usageerr();
    
    // get string command
    mxGetString(prhs[0],cmd,CMDLEN);

    // Special case. If we're called with the special command "PREINIT", then
    // we return immediately. moglcore('PREINIT') is called by ptbmoglinit.m
    // on M$-Windows in order to preload the moglcore Mex-file into Matlab while
    // the current working directory is set to ..MOGL/core/ . This way, the dynamic
    // linker can find our own local version of glut32.dll and link it against moglcore.
    // Without this trick, we would need to install glut32.dll into the Windows system
    // folder which requires admin privileges and makes installation of Psychtoolbox
    // more complicated on M$-Windows...
    if (strcmp(cmd, "PREINIT")==0) {
        glBeginLevel=0;
        goto moglreturn;
    }
    
    if (strcmp(cmd, "DUMMYMODE")==0) {
        if (nrhs<2 || mxGetScalar(prhs[1])<0) {
            mexErrMsgTxt("MOGL-ERROR: No dummy mode level or invalid level (<0) given for subcommand DUMMYMODE!");
        }
        dummymode = (int) mxGetScalar(prhs[1]);
	if (dummymode>0) printf("MOGL-INFO: Switched to dummy mode level %i.\n", dummymode);
        goto moglreturn;
    }

    // Special command to set MOGL debug level:
    // debuglevel = 0 --> Shut up in any case, leave error-handling to higher-level code.
    // debuglevel > 0 --> Check for OpenGL error conditions.
    // debuglevel > 0 --> Output glError()'s in clear-text and abort. Output GLSL errors and abort.
    // debuglevel > 1 --> Output GLSL diagnostic messages as well.
    // debuglevel > 2 --> Be very verbose!
    if (strcmp(cmd, "DEBUGLEVEL")==0) {
        if (nrhs<2 || mxGetScalar(prhs[1])<0) {
            mexErrMsgTxt("MOGL-ERROR: No debuglevel or invalid debuglevel (<0) given for subcommand DEBUGLEVEL!");
        }

        debuglevel = (int) mxGetScalar(prhs[1]);
        goto moglreturn;
    }
    
    // Special cleanup subcommand needed for GNU/Octave: See explanation below in firstTime init.
    if (strcmp(cmd, "JettisonModuleHelper")==0) {
      goto moglreturn;
    }

    // Abort here if dummymode >= 10: Input arg. processing run, but no real
    // command parsing and processing;
    if (dummymode >= 10) {
      printf("MOGL-INFO: Dummy call to %s() - Ignored in dummy mode %i ...\n", cmd, dummymode);
      goto moglreturn;
    }

    #ifdef BUILD_GLM
    // GLM module is included and supported in moglcore: This is necessary if
    // one wants to use MOGL independent from Psychtoolbox. GLM is only supported
    // on MacOS-X, not on Linux or Windows...
    
    // We execute glm-commands without performing GLEW first-time initialization,
    // because to execute glewinit() we need a valid OpenGL context. This context is
    // either created by Psychtoolbox or by glm. Therefore glm-commands must be able
    // to execute before glewinit() happened.
    
    // look for command in glm command map
    if( (i=binsearch(glm_map,glm_map_count,cmd))>=0 ) {
        glm_map[i].cmdfn(nlhs,plhs,nrhs-1,prhs+1);
        goto moglreturn;
    }

    #endif

    // Is this the first invocation of moglcore?
    if (firsttime) {
        // Yes. Initialize GLEW, the GL Extension Wrangler Library. This will
        // auto-detect and dynamically link/bind all core OpenGL functionality
        // as well as all possible OpenGL extensions on OS-X, Linux and Windows.
        err = GLEW_OK;
        #ifdef PTB_USE_WAFFLE
        // Linux is special: If we use the Waffle backend for display system binding, then our display backend
        // may be something else than GLX (e.g., X11/EGL, Wayland/EGL, GBM/EGL, ANDROID/EGL etc.), in which case
        // glewInit() would not work and would crash hard. Detect if we're on classic Linux or Linux with X11/GLX.
        // If so, execute glewInit(), otherwise call glewContextInit() - a routine which skips GLX specific setup,
        // therefore only initializes the non-GLX parts. We need a hacked glew.c for this function to be available,
        // as original upstream GLEW makes that function private (static):
        if (!getenv("PSYCH_USE_DISPLAY_BACKEND") || strstr(getenv("PSYCH_USE_DISPLAY_BACKEND"), "glx")) {
            // Classic backend or GLX backend: The full show.
            err = glewInit();
        }
        else {
            // Non-GLX backend, probably EGL: Reduced show.
            err = glewContextInit();
        }
        #else
            // Other os'es, or Linux without Waffle backend: Always init GLEW:
            err = glewInit();
        #endif

        if (GLEW_OK != err) {
            // Failed! Something is seriously wrong - We have to abort :(
            printf("MOGL: Failed to initialize! Probably you called an OpenGL command *before* opening an onscreen window?!?\n");
            printf("GLEW reported the following error: %s\n", glewGetErrorString(err)); fflush(NULL);
            goto moglreturn;
        }

        // Success. Ready to go...
		if (debuglevel > 1) {
			printf("MOGL - OpenGL for Matlab & GNU/Octave initialized. MOGL is (c) 2006-2013 Richard F. Murray & Mario Kleiner, licensed to you under MIT license.\n");
            #ifdef WINDOWS
			printf("On MS-Windows, we make use of the freeglut library, which is Copyright (c) 1999-2000 Pawel W. Olszta, licensed under compatible MIT/X11 license.\n");
            printf("The precompiled Windows binary DLL's have been kindly provided by http://www.transmissionzero.co.uk/software/freeglut-devel/ -- Thanks!\n");
            #endif
			printf("See file 'License.txt' in the Psychtoolbox root folder for the exact licensing conditions.\n");
		}
        fflush(NULL);
        
        // Perform dynamic rebinding of ARB extensions to core functions, if necessary:
        mogl_rebindARBExtensionsToCore();
        
		#ifdef FREEGLUT
		// FreeGlut must be initialized, otherwise it will emergency abort the whole application.
		// However, we skip init if we're on a setup without GLX display backend, as this would
		// abort us due to lack of GLX. On non-GLX we simply can't use FreeGlut at all.
		if (!getenv("PSYCH_USE_DISPLAY_BACKEND") || strstr(getenv("PSYCH_USE_DISPLAY_BACKEND"), "glx")) {
			// GLX display backend - Init and use FreeGlut:
			glutInit( &noargs, &dummyargp);
		}
		#endif

        // Running on a OpenGL-ES rendering api under Linux?
        if (getenv("PSYCH_USE_GFX_BACKEND") && strstr(getenv("PSYCH_USE_GFX_BACKEND"), "gles")) {
            // Yes. We emulate some immediate mode rendering commands, which aren't available
            // in OpenGL Embedded Subset at all, via "our" own emulation code. This code emulates
            // immediate mode on top of client vertex arrays and batch submission.
            if (debuglevel > 1) {
                printf("OpenGL-ES rendering API active: Emulating immediate mode rendering via David Petrie's ftglesGlue emulation code.\n");
            }
        }

		// Register exit-handler: When flushing the mex-file, we free all allocated buffer memory:
		mexAtExit(&mexExitFunction);
		
        // Done with first time initialization:
        firsttime = 0;
    }   
	
    // If glBeginLevel >  1 then most probably the script was aborted after execution of glBegin() but
    // before execution of glEnd(). In that case, we reset the level to zero.
    if (glBeginLevel > 1) glBeginLevel = 0;

    // Reset OpenGL error state so we can be sure that any of our error queries really
    // relate to errors caused by us:
    if (glBeginLevel == 0 && debuglevel > 0 && (strstr(cmd, "glGetError")==NULL)) glGetError();
        
    // look for command in manual command map
    if( (i=binsearch(gl_manual_map,gl_manual_map_count,cmd))>=0 ) {
        gl_manual_map[i].cmdfn(nlhs,plhs,nrhs-1,(const mxArray**) prhs+1);
        if (debuglevel > 0) mogl_checkerrors(cmd, prhs);
        goto moglreturn;
    }
    
    // look for command in auto command map
    if( (i=binsearch(gl_auto_map,gl_auto_map_count,cmd))>=0 ) {
        gl_auto_map[i].cmdfn(nlhs,plhs,nrhs-1,(const mxArray**) prhs+1);
        if (debuglevel > 0) mogl_checkerrors(cmd, prhs);
        goto moglreturn;
    }
    
    // no match
    mogl_usageerr();
    
    // moglreturn: Is the exit path of mogl. All execution ends at this point.
 moglreturn:
    return;
}

// do binary search in a command map for a command string
int binsearch(cmdhandler *map, int mapsize, char *str) {
    int m=0,n=mapsize-1,i,k,count=0;
    while( m<=n && count < 100) {
      count++;
        i=(int)((m+n)/2);
        k=strcmp(str,map[i].cmdstr);
        if( k==0 )
            return( i );
        else if( k<0 )
            n=i-1;
        else
            m=i+1;
    }
    return( -1 );
}

// error handler
void mogl_usageerr() {
    glBeginLevel = 0;
    printf("MOGL-Command: %s\n", cmd);
    mexErrMsgTxt("invalid moglcore command");
}

// Error handler for unsupported core OpenGL functions or extensions.
// This handler gets called by the subroutines in gl_auto.c and gl_manual.c
// whenever a gl function is not bound == not supported by current OS/driver/gfx-hardware.
// As we use the GLEW library to dynamically detect and bind all OpenGL functions, we can
// easily check if a function is supported, e.g.,
// if (glCreateShader == NULL) mogl_glunsupported("glCreateShader");
//
void mogl_glunsupported(const char* fname)
{
    char errtxt[1000];
    sprintf(errtxt, "MOGL-Error: Your Matlab code tried to call the OpenGL function %s(), which is not supported\n"
                    "MOGL-Error: by your combination of graphics hardware + graphics device driver.\n"
                    "MOGL-Error: You'll have to download+install the latest gfx-drivers for your gfx-hardware\n"
                    "MOGL-Error: or upgrade your gfx-hardware with more recent one to use this function. Aborted.\n\n", fname);

    glBeginLevel = 0;
    
    // Exit to Matlab prompt with error message:
    mexErrMsgTxt(errtxt);
}

// Function prototype for exit with printf(), like mogl_glunsupported...
void mogl_printfexit(const char* str)
{
    // Reset recursion level for glBegin() / glEnd() statements:
    glBeginLevel = 0;
    
    // Exit to Matlab prompt with error message in 'str':
    mexErrMsgTxt(str);
}

void mogl_checkerrors(const char* cmd, const mxArray *prhs[]) 
{
    char errtxt[10000];
    GLint err, status, handle;
    
    // Reject no-op calls:
    if (debuglevel<=0) return;
    
    // Accounting for glBegin and glEnd calls...
    if (strcmp(cmd, "glBegin")==0) glBeginLevel++;
    if (strcmp(cmd, "glEnd")==0 && glBeginLevel > 0) glBeginLevel--;
    
    // Skip error-checking if we are in the middle of a glBegin(), glEnd() pair.
    // Calling glGetError() in between these two commands itself is an error!
    if (glBeginLevel > 0) return;
    
    // Check for glErrors():
    if ((err=glGetError())>0) {
        // Last command caused an OpenGL error condition: Report it and abort.
        sprintf(errtxt, "MOGL-Error: Your OpenGL command %s() caused the following OpenGL error: %s. Aborted.\n", cmd, (const char*) gluErrorString(err));
        // Exit to Matlab prompt with error message:
        glBeginLevel = 0;
        mexErrMsgTxt(errtxt);
    }
    
    // No OpenGL core system errors. Check if a special GLSL command was executed:
    if (strcmp(cmd, "glCompileShader")==0) {
        // A GLSL shader got just compiled. Check its compile status...
        handle = (GLuint) mxGetScalar(prhs[1]);
        glGetShaderiv(handle, GL_COMPILE_STATUS, &status);
        if (!status) printf("MOGL-ERROR: Compilation of the GLSL shader object %i via glCompileShader(%i) failed!\n", handle, handle);
        if (debuglevel>1 || !status) {
            // Output shader info-log:
            glGetShaderInfoLog(handle, 9999, NULL, (GLchar*) &errtxt);
            printf("The shader info log for shader %i tells us the following:\n", handle);
            printf("%s \n\n", errtxt);
        }

        // Exit to Matlab prompt with error message if an error happened.
        if (!status) mexErrMsgTxt("Shader compilation failed!"); else return;
    }
        
    if (strcmp(cmd, "glLinkProgram")==0) {
        // A GLSL shader got just compiled. Check its compile status...
        handle = (GLuint) mxGetScalar(prhs[1]);
        glGetProgramiv(handle, GL_LINK_STATUS, &status);
        if (!status) printf("MOGL-ERROR: Linking of the GLSL shader program %i via glLinkProgram(%i) failed!\n", handle, handle);
        if (debuglevel>1 || !status) {
            // Output shader info-log:
            glGetProgramInfoLog(handle, 9999, NULL, (GLchar*) &errtxt);
            printf("The program info log for program %i tells us the following:\n", handle);
            printf("%s \n\n", errtxt);
        }
        
        if (debuglevel>1) {
            // Output shader info-log after program validation:
            glValidateProgram(handle);
            glGetProgramInfoLog(handle, 9999, NULL, (GLchar*) &errtxt);
            printf("The program info log for program %i tells us the following after calling glValidateProgram():\n", handle);
            printf("%s \n\n", errtxt);
        }

        // Exit to Matlab prompt with error message if an error happened.
        if (!status) mexErrMsgTxt("GLSL link operation failed!"); else return;
    }

    return;
}

// Our memory buffer allocator, adapted from Psychtoolboxs PsychMemory.c
// allocator:

#define PTBTEMPMEMDEC(n,m) totalTempMemAllocated[(m)] -=(n)

// Enqueues a new record into our linked list of temp. memory buffers.
// Returns the memory pointer to be passed to rest of Psychtoolbox.
void* PsychEnqueueTempMemory(void* p, size_t n, int mlist)
{
  // Add current buffer-head ptr as next-pointer to our new buffer:
  *((size_t*) p) = (size_t) PsychTempMemHead[mlist];

  // Set our buffer as new head of list:
  PsychTempMemHead[mlist] = (size_t*) p;

  // Add allocated buffer size as 2nd element:
  p = (unsigned char*) p + sizeof(PsychTempMemHead[mlist]);
  *((size_t*) p) = n;

  // Accounting:
  totalTempMemAllocated[mlist] += n;

  // Increment p again to get real start of user-visible buffer:
  p = (unsigned char*) p + sizeof(n);

  // MK TODO FIXME: This can print bogus values if bytecount > INT_MAX ...
  if (debuglevel > 1) mexPrintf("MOGL: Memlist %i : Allocated new buffer %p of %i Bytes,  new total = %i.\n", mlist, p, n, totalTempMemAllocated[mlist]); fflush(NULL);

  // Return ptr:
  return(p);
}

void *PsychCallocTemp(size_t n, size_t size, int mlist)
{
  void *ret;
  // MK: This could create an overflow if product n * size is
  // bigger than length of size_t --> Only
  // happens if more than 4 GB of RAM are allocated at once on
  // a 32-bit system.
  // --> Improbable for PTB, unless someones trying a buffer
  // overflow attack -- PTB would lose there badly anyway...
  size_t realsize = n * size + sizeof(void*) + sizeof(realsize);

  // realsize has extra bytes allocated for our little header...  
  if(NULL==(ret=calloc((size_t) 1, realsize))) {
    mexErrMsgTxt("MOGL-FATAL ERROR: Out of memory in PsychCallocTemp!\n");
  }

  // Need to enqueue memory buffer...
  return(PsychEnqueueTempMemory(ret, realsize, mlist));
}

void *PsychMallocTemp(size_t n, int mlist)
{
  void *ret;

  // Allocate some extra bytes for our little header...
  n = n + sizeof(void*) + sizeof(n);
  if(NULL==(ret=malloc(n))){
    mexErrMsgTxt("MOGL-FATAL ERROR: Out of memory in PsychMallocTemp!\n");
  }

  // Need to enqueue memory buffer...
  return(PsychEnqueueTempMemory(ret, n, mlist));
}

// Free a single spec'd temp memory buffer.
// TODO Note: The current implementation of our allocator
// uses a single-linked list, which has O(1) cost for
// allocating memory (Optimal!) and O(n) cost for freeing
// all allocated memory (Optimal!), but it has up to
// O(n) cost for deleting a single memory buffer as well,
// be n the total number of allocated buffers. This is
// worst-case upper bound. If PsychFreeTemp() is used a
// lot on long buffer lists, this will incur significant
// overhead! A better implementation would use a double-
// linked list or even a binary tree or hash structure,
// but for now this has to be good enough(TM).
void PsychFreeTemp(void* ptr, int mlist)
{
  void* ptrbackup = ptr;
  size_t* next = PsychTempMemHead[mlist];
  size_t* prevptr = NULL;

  if (ptr == NULL) return;
 
  // Convert ptb supplied pointer ptr into real start
  // of our buffer, including our header:
  ptr = (unsigned char*) ptr - sizeof((unsigned char*) ptr) - sizeof(size_t);
  if (ptr == NULL) return;

  if (PsychTempMemHead[mlist] == ptr) {
    // Special case: ptr is first buffer in queue. Dequeue:
    PsychTempMemHead[mlist] = (size_t*) *(PsychTempMemHead[mlist]);

    // Some accounting:
    PTBTEMPMEMDEC(((size_t*)ptr)[1], mlist);
    if (debuglevel > 1) mexPrintf("MOGL: Memlist %i : Freed buffer at %p, new total = %i.\n", mlist, ptrbackup, totalTempMemAllocated[mlist]); fflush(NULL);

    // Release it:
    free(ptr);

    return;
  }

  // ptr valid and not first buffer in queue.
  // Walk the whole buffer list until we encounter our buffer:
  while (next != NULL && next!=ptr) {
    prevptr = next;
    next = (size_t*) *next;
  }

  // Done with search loop. Did we find our buffer?
  if (next == ptr) {
    // Found it! Set next-ptr of previous buffer to next-ptr
    // of this buffer to dequeue from list:
    *prevptr = *next;

    // Some accounting:
    PTBTEMPMEMDEC(next[1], mlist);
    if (debuglevel > 1) mexPrintf("MOGL: Memlist %i: Freed buffer at %p, new total = %i.\n", mlist, ptrbackup, totalTempMemAllocated[mlist]); fflush(NULL);
    
    // Release:
    free(ptr);

    // Done.
    return;
  }

  // Oops.: Did not find matching buffer to pointer --> Trouble!
  printf("MOGL-BUG: In PsychFreeTemp: Tried to free non-existent temporary membuffer %p!!! Ignored.\n", ptrbackup);
  fflush(NULL);
  return;
}

// Master cleanup routine: Frees all allocated memory:
void PsychFreeAllTempMemory(int mlist)
{
  size_t* p = NULL;
  size_t* next = PsychTempMemHead[mlist];

  // Walk our whole buffer list and release all buffers on it:
  while (next != NULL) {
    // next points to current buffer to release. Make a copy of
    // next:
    p = next;

    // Update next to point to the next buffer to release:
    next = (size_t*) *p;

    // Some accounting:
    PTBTEMPMEMDEC(p[1], mlist);

    // Release buffer p:
    free(p);
    
    // We're done with this buffer, next points to next one to release
    // or is NULL if all released...
  }

  // Done. NULL-out the list start ptr:
  PsychTempMemHead[mlist] = NULL;

  // Sanity check:
  if (totalTempMemAllocated[mlist] != 0) {
    printf("MOGL-CRITICAL BUG: Inconsistency detected in temporary memory allocator!\n");
    printf("MOGL-CRITICAL BUG: totalTempMemAllocated[%i] = %i after PsychFreeAllTempMemory(%i)!!!!\n",
	   mlist, totalTempMemAllocated[mlist], mlist);
    fflush(NULL);

    // Reset to defined state.
    totalTempMemAllocated[mlist] = 0;
  }

  if (debuglevel > 1) printf("MOGL: Memlist %i : Freed all internal memory buffers.\n", mlist); fflush(NULL);

  return;
}

// Convert a double value (which encodes a memory address) into a ptr:
void*  PsychDoubleToPtr(volatile double dptr)
{
  volatile psych_uint64* iptr = (psych_uint64*) &dptr;
  volatile psych_uint64 ival = *iptr;
  return((void*) ((size_t) ival));
}

// Convert a memory address pointer into a double value:
double PsychPtrToDouble(void* ptr)
{
  volatile psych_uint64 ival = (psych_uint64) ((size_t) ptr);
  volatile double* dptr = (double*) &ival;
  volatile double outval = *dptr;
  return(outval);
}

// Return the size of the buffer pointed to by ptr in bytes.
// CAUTION: Only works with buffers allocated via PsychMallocTemp()
// or PsychCallocTemp(). Will segfault, crash & burn with other pointers!
// The routine returns the net-size of the buffer (the size useable by
// code), not the allocated size (which would be a few additional bytes
// for the buffer-header).
size_t PsychGetBufferSizeForPtr(void* ptr)
{
  size_t mysize;

  // Decrement pointer to let it point to our size field in the buffer-header:
  ptr = (unsigned char*) ptr - sizeof(size_t);

  // Retrieve size value in header, subtract size of header itself:
  mysize = (*((size_t*) ptr)) - sizeof(void*) - sizeof(size_t);

  // Return it:
  return((size_t) mysize);
}

// Mapping of scalar buffer offset value (in units of bytes) to an
// equivalent memory void*. Handles doubles, uint32 and uint64:
void* moglScalarToPtrOffset(const mxArray *m) {
	if (mxIsDouble(m)) return((void*) (size_t) mxGetScalar(m));
	if (mxIsUint32(m)) return((void*) (size_t) (((unsigned int*) mxGetData(m))[0]));

    #ifndef MATLABR11
    if (mxIsUint64(m)) return((void*) (size_t) (((psych_uint64*) mxGetData(m))[0]));
	#endif
    
	// Invalid input type - Error abort:
    glBeginLevel = 0;
    printf("MOGL-Command: %s\n", cmd);
	mexErrMsgTxt("Provided pointer or buffer offset argument is of invalid type - Not double, uint32 or uint64!");
	return(NULL);
}
