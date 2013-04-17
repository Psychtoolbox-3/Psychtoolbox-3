
/*
 * gl_manual.c
 *
 * 19-Dec-2005 -- created (RFM)
 * 24-Mar-2011 -- Make 64-bit clean (MK).
 *
 */

#include "mogltypes.h"
#include "ftglesGlue.h"

void gl_getbufferpointerv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
	glGetBufferPointerv((GLenum)mxGetScalar(prhs[0]),
        (GLenum)mxGetScalar(prhs[1]),
        (GLvoid **)mxGetData(prhs[2]));

}

void gl_getpointerv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
	glGetPointerv((GLenum)mxGetScalar(prhs[0]),
		(GLvoid **)mxGetData(prhs[1]));

}

void gl_getstring( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    plhs[0]=mxCreateString((const char *)glGetString((GLenum)mxGetScalar(prhs[0])));

}

void gl_getvertexattribpointerv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	glGetVertexAttribPointerv((GLuint)mxGetScalar(prhs[0]),
        (GLenum)mxGetScalar(prhs[1]),
        (GLvoid **)mxGetData(prhs[2]));

}

void glu_errorstring( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    plhs[0]=mxCreateString((const char *)gluErrorString((GLenum)mxGetScalar(prhs[0])));

}

void glu_getstring( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    plhs[0]=mxCreateString((const char *)gluGetString((GLenum)mxGetScalar(prhs[0])));

}

void gl_samplepass( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    // MK: This function is a zombie. It was specified in some early draft of the
	// multisample extension and then removed from the spec. Unfortunately its definition
	// is still floating around in some sdk files on the internet, although the function
	// is not implemented in any OpenGL driver. We define the function but error out
	// if someone tries to use it:
	mogl_glunsupported("glSamplePass");
}

void gl_shadersource( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    char* sourcestring;
    char** srcstrings;
    GLuint handle;
    size_t i, count, savedlength;
    count = 0;
    
    // Ok, glShaderSource needs a list of one-line strings for the single lines
    // of code in a shader program. We take a single big string of newline-
    // separated lines of code and break it up into a list of strings.
    
    // Retrieve handle to shader:
    handle = (GLuint) mxGetScalar(prhs[0]);
    // Retrieve shader source code string:
    sourcestring = mxArrayToString(prhs[1]);
    // Count number of lines in string:
    savedlength = strlen(sourcestring);
    for (i=0; i<savedlength; i++) if(sourcestring[i]=='\n') count++;
    count++;
    // printf("COUNT %i\n", count); fflush(NULL);
    // Allocate char* array of proper capacity:
    srcstrings=(char**) malloc(count * sizeof(char*));
    // Initialize array of char-ptrs and setup string:
    count=0;
    srcstrings[0]=(char*) sourcestring;
    count++;
    for (i=0; i<savedlength; i++) if(sourcestring[i]=='\n') {
        // NULL-out the end-of-line terminator to create a null-terminated piece of
        // substring:
        sourcestring[i]=0;
        // Setup a new char* that points behind the sourcestring[i]:
        i++;
        srcstrings[count]=(char*) &(sourcestring[i]);
        count++;
    }

    if (mxGetScalar(prhs[2])>0) {
        printf("\n\n");
        for(i=0; i<count; i++) printf("Shader Line %i: %s\n", (int) i, srcstrings[i]);
        printf("\n\n");
        fflush(NULL);
        // Free the sourcestring:
        mxFree(sourcestring);
        // Free our array:
        free(srcstrings);
        // Abort with error:
        mogl_printfexit("MOGL-Info: In call to glShaderSource(): Called with debug flag set to non-zero.\nMOGL-Info: Therefore i've dumped the shader source code to the console and will now exit.");
    }
    
    if (NULL == glShaderSource) {        
        // Free the sourcestring:
        mxFree(sourcestring);
        // Free our array:
        free(srcstrings);
        // Abort with error:
        mogl_glunsupported("glShaderSource");
    }
    
    // Ok, now srcstrings should be an array of count char*'s to single line, null-terminated
    // strings, suitable for glShaderSource. Call it.
    glGetError();
    glShaderSource(handle, (GLsizei) count, (const char**) srcstrings, NULL);
    
    // Free the sourcestring:
    mxFree(sourcestring);
    // Free our array:
    free(srcstrings);
    
    // Done.
    return;
}

void gl_feedbackbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

  // Retrieve memory ptr from double argument:
  GLfloat* ptr = (GLfloat*) PsychDoubleToPtr((GLdouble) mxGetScalar(prhs[2]));

  if (NULL == glFeedbackBuffer) mogl_glunsupported("glFeedbackBuffer");

  glFeedbackBuffer((GLsizei)mxGetScalar(prhs[0]),
		   (GLenum)mxGetScalar(prhs[1]),
		   ptr);
}

void gl_selectbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

  // Retrieve memory ptr from double argument:
  GLuint* ptr = (GLuint*) PsychDoubleToPtr((GLdouble) mxGetScalar(prhs[1]));

  if (NULL == glSelectBuffer) mogl_glunsupported("glSelectBuffer");
  glSelectBuffer((GLsizei)mxGetScalar(prhs[0]), ptr);
}

void gl_teximage2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexImage2D) mogl_glunsupported("glTexImage2D");
	glTexImage2D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLenum)mxGetScalar(prhs[6]),
		(GLenum)mxGetScalar(prhs[7]),
        (const GLvoid*) ((mxGetM(prhs[8]) * mxGetN(prhs[8]) > 1) ? mxGetData(prhs[8]) : moglScalarToPtrOffset(prhs[8])));
}

void moglmalloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  // Allocate a memory buffer of prhs[0] bytes size. ptr points to start of buffer:
  void* ptr = PsychMallocTemp((size_t) mxGetScalar(prhs[0]), 1);

  // Convert ptr into a double value and assign it as first return argument:
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  *((GLdouble*) mxGetPr(plhs[0])) = PsychPtrToDouble(ptr);
}

void moglcalloc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  // Allocate a memory buffer of prhs[0] bytes size. ptr points to start of buffer:
  void* ptr = PsychCallocTemp((size_t) mxGetScalar(prhs[0]), (size_t) mxGetScalar(prhs[1]), 1);

  // Convert ptr into a double value and assign it as first return argument:
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  *((GLdouble*) mxGetPr(plhs[0])) = PsychPtrToDouble(ptr);
}

void moglfree(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  // Retrieve ptr to membuffer, encoded as double and convert it into a void*
  void* ptr = PsychDoubleToPtr((GLdouble) mxGetScalar(prhs[0]));

  // Free memory buffer:
  PsychFreeTemp(ptr, 1);
}

void moglfreeall(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  PsychFreeAllTempMemory(1);
}

void moglcopybuffertomatrix(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  int type;
  size_t dims;
  mwSize outdims;
  GLenum mattype;
  GLfloat* dst;

  // Retrieve ptr to membuffer, encoded as double and convert it into a void*
  void* src = PsychDoubleToPtr((GLdouble) mxGetScalar(prhs[0]));

  // Retrieve size of buffer pointed to by src:
  size_t n = PsychGetBufferSizeForPtr(src);

  // Retrieve max number of bytes to copy:
  size_t nmax = (size_t) mxGetScalar(prhs[2]);

  if (nmax < n) n = nmax;

  // Retrieve type of matrix to create:
  mattype = (GLenum) mxGetScalar(prhs[1]);

  switch(mattype)
    {
    case GL_DOUBLE:
      type = mxDOUBLE_CLASS;
      dims = n / sizeof(double);
      if (dims * sizeof(double) < n) dims++;
      break;
    case GL_FLOAT:
      type = mxSINGLE_CLASS;
      dims = n / sizeof(float);
      if (dims * sizeof(float) < n) dims++;
      break;
    case GL_UNSIGNED_INT:
      type = mxUINT32_CLASS;
      dims = n / sizeof(unsigned int);
      if (dims * sizeof(unsigned int) < n) dims++;
      break;
    case GL_UNSIGNED_BYTE:
      type = mxUINT8_CLASS;
      dims = n / sizeof(unsigned char);
      if (dims * sizeof(unsigned char) < n) dims++;
      break;
    default:
	  type = 0;
      mexErrMsgTxt("MOGL-ERROR: Unknown matrix type requested in moglgetbuffer()! Ignored.");
    }

  if ((sizeof(outdims) < sizeof(dims)) && (dims > INT_MAX)) {
      mexErrMsgTxt("MOGL-ERROR: In moglgetbuffer(): Returned matrix is too big (more than 2^31 - 1 elements) for your version of Matlab or Octave! Aborted.");
  }
  
  // Allocate the beast:
  outdims = (mwSize) dims;
  plhs[0] = mxCreateNumericArray(1, &outdims, type, mxREAL);

  // Retrieve pointer to output matrix:
  dst = (GLfloat*) mxGetData(plhs[0]);

  // Do the copy:
  memcpy(dst, src, n);
}

void moglcopymatrixtobuffer(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  // Retrieve ptr to membuffer, encoded as double and convert it into a void*
  void* dst = PsychDoubleToPtr((GLdouble) mxGetScalar(prhs[1]));

  // Retrieve pointer to input matrix:
  GLfloat* src = (GLfloat*) mxGetData(prhs[0]);

  // Retrieve size of buffer pointed to by dst:
  size_t nmax = PsychGetBufferSizeForPtr(dst);

  // Set final size of data to copy:
  size_t nin = (size_t) mxGetScalar(prhs[2]);
  if (nin > nmax) nin = nmax;

  // Do the copy:
  memcpy(dst, src, nin);
}

void gl_bufferdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBufferData) mogl_glunsupported("glBufferData");
	glBufferData((GLenum)mxGetScalar(prhs[0]),
		(GLsizeiptr)mxGetScalar(prhs[1]),
		(const GLvoid*) ((mxGetM(prhs[2]) * mxGetN(prhs[2]) > 1) ? mxGetData(prhs[2]) : moglScalarToPtrOffset(prhs[2])),
		(GLenum)mxGetScalar(prhs[3]));

}

void gl_readpixels( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glReadPixels) mogl_glunsupported("glReadPixels");
	glReadPixels((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(GLvoid*) ((mxGetM(prhs[6]) * mxGetN(prhs[6]) > 1) ? mxGetData(prhs[6]) :  moglScalarToPtrOffset(prhs[6])));

}

void gl_vertexpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexPointer) mogl_glunsupported("glVertexPointer");
	glVertexPointer((GLint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
	        (const GLvoid*) ((mxGetM(prhs[3]) * mxGetN(prhs[3]) > 1) ? mxGetData(prhs[3]) :  moglScalarToPtrOffset(prhs[3])));

}

void gl_normalpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormalPointer) mogl_glunsupported("glNormalPointer");
	glNormalPointer((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLvoid*) ((mxGetM(prhs[2]) * mxGetN(prhs[2]) > 1) ? mxGetData(prhs[2]) :  moglScalarToPtrOffset(prhs[2])));

}

void gl_texcoordpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordPointer) mogl_glunsupported("glTexCoordPointer");
	glTexCoordPointer((GLint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLvoid*) ((mxGetM(prhs[3]) * mxGetN(prhs[3]) > 1) ? mxGetData(prhs[3]) :  moglScalarToPtrOffset(prhs[3])));

}

void gl_colorpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorPointer) mogl_glunsupported("glColorPointer");
	glColorPointer((GLint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLvoid*) ((mxGetM(prhs[3]) * mxGetN(prhs[3]) > 1) ? mxGetData(prhs[3]) :  moglScalarToPtrOffset(prhs[3])));

}

void gl_vertexattribpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribPointer) mogl_glunsupported("glVertexAttribPointer");
	glVertexAttribPointer((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(const GLvoid*) ((mxGetM(prhs[5]) * mxGetN(prhs[5]) > 1) ? mxGetData(prhs[5]) :  moglScalarToPtrOffset(prhs[5])));

}

void gl_drawelements( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawElements) mogl_glunsupported("glDrawElements");
	glDrawElements((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(const GLvoid*) ((mxGetM(prhs[3]) * mxGetN(prhs[3]) > 1) ? mxGetData(prhs[3]) :  moglScalarToPtrOffset(prhs[3])));

}

void gl_drawrangeelements( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawRangeElements) mogl_glunsupported("glDrawRangeElements");
	glDrawRangeElements((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(const GLvoid*) ((mxGetM(prhs[5]) * mxGetN(prhs[5]) > 1) ? mxGetData(prhs[5]) :  moglScalarToPtrOffset(prhs[5])));

}

void glu_newtess( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    mogl_tess_struct* mytess;
	if (NULL == gluNewTess) mogl_glunsupported("gluNewTess");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);

    // Create our own virtual tesselator struct:
    mytess = (mogl_tess_struct*) PsychCallocTemp(1, sizeof(mogl_tess_struct), 2);
    
    // Alloc a gluTesselator and store its ptr inside our struct:
    mytess->glutesselator = gluNewTess();
    
    // Don't need to init other parts of struct, as calloc() inits them to
    // all-zero, which is our wanted default.
    
    // Assign ptr to our own struct as returned tesselator pointer:
	*(double*) mxGetData(plhs[0]) = PsychPtrToDouble((void*) mytess);

}

void glu_deletetess( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    // Define mytess struct pointer for this tesselator:
    MOGLDEFMYTESS;
	
	if (NULL == gluDeleteTess) mogl_glunsupported("gluDeleteTess");
        
    // Destroy GLUtesselator:
	gluDeleteTess((GLUtesselator*) mytess->glutesselator);
        
    // Release destructBuffers:
    PsychFreeAllTempMemory(3);    
    
    // Release our own data structure:
    PsychFreeTemp((void*) mytess, 2);

}

void glu_tessbegincontour( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    MOGLDEFMYTESS;
	if (NULL == gluTessBeginContour) mogl_glunsupported("gluTessBeginContour");

	gluTessBeginContour((GLUtesselator*) mytess->glutesselator);

}

void CALLCONV mogl_GLU_TESS_BEGIN_DATA(GLenum type, void* polygondata)
{
    mogl_tess_struct* mytess = (mogl_tess_struct*) polygondata;
    mxArray* prhs[2];
    //mexPrintf("BEGINCB\n");
    
    // Assign type as 1st argument of double type:
    prhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    *(double*) mxGetData(prhs[0]) = (double) type;
    
    if (mytess->userData) {
        // Assign user-provided polygondata as double:
		prhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		*(double*) mxGetData(prhs[1]) = mytess->polygondata;
        mexCallMATLAB(0, NULL, 2, prhs, mytess->nGLU_TESS_BEGIN_DATA);
    }
    else {
        // Don't assign polygondata:
        mexCallMATLAB(0, NULL, 1, prhs, mytess->nGLU_TESS_BEGIN);
    }
}

void CALLCONV mogl_GLU_TESS_END_DATA(void* polygondata)
{
    mogl_tess_struct* mytess = (mogl_tess_struct*) polygondata;
    mxArray* prhs[1];
    //mexPrintf("ENDCB\n");
        
    if (mytess->userData) {
        // Assign user-provided polygondata as double:
		prhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		*(double*) mxGetData(prhs[0]) = mytess->polygondata;
        mexCallMATLAB(0, NULL, 1, prhs, mytess->nGLU_TESS_END_DATA);
    }
    else {
        // Don't assign polygondata:
        mexCallMATLAB(0, NULL, 0, NULL, mytess->nGLU_TESS_END);
    }
}

void CALLCONV mogl_GLU_TESS_EDGE_FLAG_DATA(GLboolean flag, void* polygondata)
{
    mogl_tess_struct* mytess = (mogl_tess_struct*) polygondata;
    mxArray* prhs[2];
    //mexPrintf("EDGECB\n");
    
    // Assign flag as 1st argument of double type:
    prhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    *(double*) mxGetData(prhs[0]) = (double) flag;
    
    if (mytess->userData) {
        // Assign user-provided polygondata as double:
		prhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		*(double*) mxGetData(prhs[1]) = mytess->polygondata;
        mexCallMATLAB(0, NULL, 2, prhs, mytess->nGLU_TESS_EDGE_FLAG_DATA);
    }
    else {
        // Don't assign polygondata:
        mexCallMATLAB(0, NULL, 1, prhs, mytess->nGLU_TESS_EDGE_FLAG);
    }
}

void CALLCONV mogl_GLU_TESS_VERTEX_DATA(void* vertex_data, void* polygondata)
{
    mogl_tess_struct* mytess = (mogl_tess_struct*) polygondata;
    mxArray* prhs[2];
    double*  ddat;
    //mexPrintf("VERTEXCB\n");


    prhs[0] = mxCreateNumericMatrix(mytess->nrElements, 1, mxDOUBLE_CLASS, mxREAL);
    ddat = (double*) mxGetData(prhs[0]);
    memcpy(ddat, vertex_data, sizeof(double) * mytess->nrElements);
        
    if (mytess->userData) {
        // Assign user-provided polygondata as double:
		prhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		*(double*) mxGetData(prhs[1]) = mytess->polygondata;
        mexCallMATLAB(0, NULL, 2, prhs, mytess->nGLU_TESS_VERTEX_DATA);
    }
    else {
        // Don't assign polygondata:
        // mexPrintf("IN VERTEXCALLBACK: %p\n", vertex_data);
        mexCallMATLAB(0, NULL, 1, prhs, mytess->nGLU_TESS_VERTEX);
    }
}

void CALLCONV mogl_GLU_TESS_ERROR_DATA(GLenum type, void* polygondata)
{
    mogl_tess_struct* mytess = (mogl_tess_struct*) polygondata;
    mxArray* prhs[2];
    //mexPrintf("ERRORCB\n");

    // Assign type as 1st argument of double type:
    prhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    *(double*) mxGetData(prhs[0]) = (double) type;
    
    if (mytess->userData) {
        // Assign user-provided polygondata as double:
		prhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		*(double*) mxGetData(prhs[1]) = mytess->polygondata;
        mexCallMATLAB(0, NULL, 2, prhs, mytess->nGLU_TESS_ERROR_DATA);
    }
    else {
        // Don't assign polygondata:
        mexCallMATLAB(0, NULL, 1, prhs, mytess->nGLU_TESS_ERROR);
    }
}

void* mogl_enqueueVertex(mogl_tess_struct* mytess, mxArray* vdat)
{
    void* dst;
    double* newdestructBuffer;
	
    mytess->nrElements = mxGetNumberOfElements(vdat);
    
    if (mytess->destructCount >= mytess->destructSize) {
        // Keep track of biggest buffersize for this tesselator so far.
        // Init to biggest size so far at initial allocation. We want to
        // grow the buffer quickly to a sufficient capacity to reduce alloc
        // overhead and memory fragmentation:
        if (mytess->destructSize < mytess->maxdestructSize) {
            mytess->destructSize = mytess->maxdestructSize;
        }
        else {
            mytess->destructSize += (1000 * mytess->nrElements);
            mytess->maxdestructSize = mytess->destructSize;
        }
        
        // Alloc:
        // mexPrintf("REALLOC VBUFFER of size %i elements.\n", mytess->destructSize);
        newdestructBuffer = (double*) PsychMallocTemp(sizeof(double) * mytess->destructSize, 3);
        if (newdestructBuffer) {
            mytess->destructBuffer = newdestructBuffer;
            mytess->destructCount  = 0;
        }
        else {
            mytess->destructBuffer = NULL;
            mytess->destructSize  = 0;
            mytess->destructCount = 0;
            PsychFreeAllTempMemory(3);
            mexErrMsgTxt("MOGL-ERROR: Out of memory error while processing gluTessCallback() or gluTessVertex()! Aborting!");
        }
    }
    
    dst = (void*) &(mytess->destructBuffer[mytess->destructCount]);
    memcpy(dst, mxGetData(vdat), sizeof(double) * mytess->nrElements);
    mytess->destructCount += mytess->nrElements;
    
    return(dst);
}

void CALLCONV mogl_GLU_TESS_COMBINE_DATA(GLdouble coords[3], void *vertex_data[4],
                                GLfloat weight[4], void **outData, void* polygondata)
{
    mxArray* prhs[4];
    mxArray* retData[1];
    double*  ddat;

    mogl_tess_struct* mytess = (mogl_tess_struct*) polygondata;

    //mexPrintf("COMBINERCB\n");

    // Assign coords as 1st argument:
    prhs[0] = mxCreateNumericMatrix(3, 1, mxDOUBLE_CLASS, mxREAL);
    ddat = (double*) mxGetData(prhs[0]);
    ddat[0] = coords[0];
    ddat[1] = coords[1];
    ddat[2] = coords[2];

    //mexPrintf("COMBINERCB2: %i , %p , %p , %p, %p\n", mytess->nrElements, vertex_data[0], vertex_data[1],vertex_data[2],vertex_data[3]);
    //mexPrintf("COMBINERCB2: %i , %f , %f , %f, %f\n", mytess->nrElements, weight[0], weight[1],weight[2],weight[3]);

    // Assign vertex_data[0-3] as 2nd argument:
    prhs[1] = mxCreateNumericMatrix(mytess->nrElements, 4, mxDOUBLE_CLASS, mxREAL);
    ddat = (double*) mxGetData(prhs[1]);
    memset(ddat, 0, sizeof(double) * 4 * mytess->nrElements);

    if (weight[0] && vertex_data[0])
        memcpy(&(ddat[0]), vertex_data[0], sizeof(double) * mytess->nrElements);

    ddat = &(ddat[mytess->nrElements]);

    if (weight[1] && vertex_data[1])
        memcpy(&(ddat[0]), vertex_data[1], sizeof(double) * mytess->nrElements);

    ddat = &(ddat[mytess->nrElements]);

    if (weight[2] && vertex_data[2])
        memcpy(&(ddat[0]), vertex_data[2], sizeof(double) * mytess->nrElements);

    ddat = &(ddat[mytess->nrElements]);

    if (weight[3] && vertex_data[3])
        memcpy(&(ddat[0]), vertex_data[3], sizeof(double) * mytess->nrElements);

    // mexPrintf("COMBINERCB3\n");


    // Assign weight as 3rd argument:
    prhs[2] = mxCreateNumericMatrix(4, 1, mxDOUBLE_CLASS, mxREAL);
    ddat = (double*) mxGetData(prhs[2]);
    ddat[0] = weight[0];
    ddat[1] = weight[1];
    ddat[2] = weight[2];
    ddat[3] = weight[3];
    //mexPrintf("COMBINERCB4\n");
    
    if (mytess->userData) {
        // Assign user-provided polygondata as double:
		prhs[3] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		*(double*) mxGetData(prhs[3]) = mytess->polygondata;
        mexCallMATLAB(1, retData, 4, prhs, mytess->nGLU_TESS_COMBINE_DATA);
    }
    else {
        // Don't assign polygondata:
        mexCallMATLAB(1, retData, 3, prhs, mytess->nGLU_TESS_COMBINE);
    }
    
    //mexPrintf("COMBINERCB5\n");

    // Assign pointer to new created databuffer location with new vertex data:
    *outData = mogl_enqueueVertex(mytess, retData[0]);
    //mexPrintf("COMBINERCBEND\n");
    
}

#ifndef MOGLSETTESSCALLBACK
#define MOGLSETTESSCALLBACK(which)                                                                                          \
{                                                                                                                           \
    if (mytess->n##which##_DATA[0]) {                                                                                       \
        gluTessCallback((GLUtesselator*) mytess->glutesselator, which##_DATA, mogl_##which##_DATA);                         \
        mytess->userData = 1;                                                                                               \
    }                                                                                                                       \
     else if (mytess->n##which[0]) {                                                                                        \
         gluTessCallback((GLUtesselator*) mytess->glutesselator, which##_DATA, mogl_##which##_DATA);                        \
         mytess->userData = 0;                                                                                              \
     }                                                                                                                      \
     else {                                                                                                                 \
         gluTessCallback((GLUtesselator*) mytess->glutesselator, which##_DATA, NULL );                                      \
     }                                                                                                                      \
}
#endif

void glu_tessbeginpolygon( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    MOGLDEFMYTESS;

	if (NULL == gluTessBeginPolygon) mogl_glunsupported("gluTessBeginPolygon");

    // Setup all tesselator callbacks for this polygon:
    MOGLSETTESSCALLBACK(GLU_TESS_BEGIN);
    MOGLSETTESSCALLBACK(GLU_TESS_EDGE_FLAG);
    MOGLSETTESSCALLBACK(GLU_TESS_VERTEX);
    MOGLSETTESSCALLBACK(GLU_TESS_END);
    MOGLSETTESSCALLBACK(GLU_TESS_COMBINE);
    MOGLSETTESSCALLBACK(GLU_TESS_ERROR);
        
    // Store user-provided polygon data pointer internally:
    mytess->polygondata = mxGetScalar(prhs[1]);
    
    // Pass pointer to our own struct as polygondata to the real tesselator:
	gluTessBeginPolygon((GLUtesselator*) mytess->glutesselator, (GLvoid*) mytess);

}

void glu_tessendcontour( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    MOGLDEFMYTESS;    
	if (NULL == gluTessEndContour) mogl_glunsupported("gluTessEndContour");

    gluTessEndContour((GLUtesselator*) mytess->glutesselator);

}

void glu_tessendpolygon( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    MOGLDEFMYTESS;    
	if (NULL == gluTessEndPolygon) mogl_glunsupported("gluTessEndPolygon");

    gluTessEndPolygon((GLUtesselator*) mytess->glutesselator);
    mytess->destructBuffer = NULL;
    mytess->destructSize  = 0;
    mytess->destructCount = 0;
    PsychFreeAllTempMemory(3);    
}

void glu_tessnormal( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    MOGLDEFMYTESS;
	if (NULL == gluTessNormal) mogl_glunsupported("gluTessNormal");

    gluTessNormal((GLUtesselator*) mytess->glutesselator,
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void glu_tessproperty( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    MOGLDEFMYTESS;    
	if (NULL == gluTessProperty) mogl_glunsupported("gluTessProperty");

    gluTessProperty((GLUtesselator*) mytess->glutesselator,
		(GLenum)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void glu_tessvertex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    void* vdat;    
    MOGLDEFMYTESS;    
    vdat = mogl_enqueueVertex(mytess, (mxArray*) prhs[2]);
    
	if (NULL == gluTessVertex) mogl_glunsupported("gluTessVertex");

    gluTessVertex((GLUtesselator*) mytess->glutesselator,
		(GLdouble*)mxGetData(prhs[1]),
		(GLvoid*) vdat);
    
    // mexPrintf("IN PPUSHVERTEX: %p\n", vdat);
}

void glu_tesscallback( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    GLenum which;
    MOGLDEFMYTESS;
	if (NULL == gluTessCallback) mogl_glunsupported("gluTessCallback");

    which = (GLenum)mxGetScalar(prhs[1]);
    switch(which) {
        case GLU_TESS_BEGIN:
            mxGetString(prhs[2], mytess->nGLU_TESS_BEGIN, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_BEGIN_DATA:
            mxGetString(prhs[2], mytess->nGLU_TESS_BEGIN_DATA, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_EDGE_FLAG:
            mxGetString(prhs[2], mytess->nGLU_TESS_EDGE_FLAG, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_EDGE_FLAG_DATA:
            mxGetString(prhs[2], mytess->nGLU_TESS_EDGE_FLAG_DATA, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_VERTEX:
            mxGetString(prhs[2], mytess->nGLU_TESS_VERTEX, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_VERTEX_DATA:
            mxGetString(prhs[2], mytess->nGLU_TESS_VERTEX_DATA, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_END:
            mxGetString(prhs[2], mytess->nGLU_TESS_END, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_END_DATA:
            mxGetString(prhs[2], mytess->nGLU_TESS_END_DATA, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_COMBINE:
            mxGetString(prhs[2], mytess->nGLU_TESS_COMBINE, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_COMBINE_DATA:
            mxGetString(prhs[2], mytess->nGLU_TESS_COMBINE_DATA, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_ERROR:
            mxGetString(prhs[2], mytess->nGLU_TESS_ERROR, MAX_TESSCBNAME);
            break;
            
        case GLU_TESS_ERROR_DATA:
            mxGetString(prhs[2], mytess->nGLU_TESS_ERROR_DATA, MAX_TESSCBNAME);
            break;
            
        default:
            mexErrMsgTxt("MOGL-ERROR: Unknown 'which' parameter given to gluTessCallback()!");
    }
    
}

// GLES emulation wrapper wrappers:
void gles_color4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == ftglColor4f) mogl_glunsupported("ftglColor4f");
	ftglColor4f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));
}

void gles_texcoord2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == ftglTexCoord2f) mogl_glunsupported("ftglTexCoord2f");
	ftglTexCoord2f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));
}

void gles_vertex2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == ftglVertex2f) mogl_glunsupported("ftglVertex2f");
	ftglVertex2f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));
}

void gles_vertex3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == ftglVertex3f) mogl_glunsupported("ftglVertex3f");
	ftglVertex3f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));
}

void gles_begin( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == ftglBegin) mogl_glunsupported("ftglBegin");
	ftglBegin((GLenum)mxGetScalar(prhs[0]));
}

void gles_end( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == ftglEnd) mogl_glunsupported("ftglEnd");
	ftglEnd();
}


// command map:  moglcore string commands and functions that handle them
// *** it's important that this list be kept in alphabetical order, 
//     and that gl_manual_map_count be updated
//     for each new entry ***
int gl_manual_map_count=42;
cmdhandler gl_manual_map[] = {
{ "ftglBegin",                      gles_begin                          },
{ "ftglColor4f",                    gles_color4f                        },
{ "ftglEnd",                        gles_end                            },
{ "ftglTexCoord2f",                 gles_texcoord2f                     },
{ "ftglVertex2f",                   gles_vertex2f                       },
{ "ftglVertex3f",                   gles_vertex3f                       },
{ "glBufferData",                   gl_bufferdata                       },
{ "glColorPointer",                 gl_colorpointer                     },
{ "glDrawElements",                 gl_drawelements                     },
{ "glDrawRangeElements",            gl_drawrangeelements                },
{ "glFeedbackBuffer",               gl_feedbackbuffer                   },
{ "glGetBufferPointerv",            gl_getbufferpointerv                },
{ "glGetPointerv",                  gl_getpointerv                      },
{ "glGetString",                    gl_getstring                        },
{ "glGetVertexAttribPointerv",      gl_getvertexattribpointerv          },
{ "glNormalPointer",                gl_normalpointer                    },
{ "glReadPixels",                   gl_readpixels                       },
{ "glSamplePass",                   gl_samplepass                       },
{ "glSelectBuffer",                 gl_selectbuffer                     },
{ "glShaderSource",                 gl_shadersource                     },
{ "glTexCoordPointer",              gl_texcoordpointer                  },
{ "glTexImage2D",                   gl_teximage2d                       },
{ "glVertexAttribPointer",          gl_vertexattribpointer              },
{ "glVertexPointer",                gl_vertexpointer                    },
{ "gluDeleteTess",                  glu_deletetess                      },
{ "gluErrorString",                 glu_errorstring                     },
{ "gluGetString",                   glu_getstring                       },
{ "gluNewTess",                     glu_newtess                         },
{ "gluTessBeginContour",            glu_tessbegincontour                },
{ "gluTessBeginPolygon",            glu_tessbeginpolygon                },
{ "gluTessCallback",                glu_tesscallback                    },
{ "gluTessEndContour",              glu_tessendcontour                  },
{ "gluTessEndPolygon",              glu_tessendpolygon                  },
{ "gluTessNormal",                  glu_tessnormal                      },
{ "gluTessProperty",                glu_tessproperty                    },
{ "gluTessVertex",                  glu_tessvertex                      },
{ "moglcalloc",                     moglcalloc                          },
{ "moglcopybuffertomatrix",         moglcopybuffertomatrix              },
{ "moglcopymatrixtobuffer",         moglcopymatrixtobuffer              },
{ "moglfree",                       moglfree                            },
{ "moglfreeall",                    moglfreeall                         },
{ "moglmalloc",                     moglmalloc                          }};
