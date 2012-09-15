
/*
 * gl_auto.c
 *
 * 30-Aug-2012 -- created (moglgen)
 *
 */

#include "mogltypes.h"

void gl_accum( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glAccum) mogl_glunsupported("glAccum");
	glAccum((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_alphafunc( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glAlphaFunc) mogl_glunsupported("glAlphaFunc");
	glAlphaFunc((GLenum)mxGetScalar(prhs[0]),
		(GLclampf)mxGetScalar(prhs[1]));

}

void gl_aretexturesresident( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glAreTexturesResident) mogl_glunsupported("glAreTexturesResident");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glAreTexturesResident((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]),
		(GLboolean*)mxGetData(prhs[2]));

}

void gl_arrayelement( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glArrayElement) mogl_glunsupported("glArrayElement");
	glArrayElement((GLint)mxGetScalar(prhs[0]));

}

void gl_begin( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBegin) mogl_glunsupported("glBegin");
	glBegin((GLenum)mxGetScalar(prhs[0]));

}

void gl_bindtexture( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindTexture) mogl_glunsupported("glBindTexture");
	glBindTexture((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_bitmap( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBitmap) mogl_glunsupported("glBitmap");
	glBitmap((GLsizei)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]),
		(GLfloat)mxGetScalar(prhs[4]),
		(GLfloat)mxGetScalar(prhs[5]),
		(const GLubyte*)mxGetData(prhs[6]));

}

void gl_blendcolor( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendColor) mogl_glunsupported("glBlendColor");
	glBlendColor((GLclampf)mxGetScalar(prhs[0]),
		(GLclampf)mxGetScalar(prhs[1]),
		(GLclampf)mxGetScalar(prhs[2]),
		(GLclampf)mxGetScalar(prhs[3]));

}

void gl_blendequation( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendEquation) mogl_glunsupported("glBlendEquation");
	glBlendEquation((GLenum)mxGetScalar(prhs[0]));

}

void gl_blendequationseparate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendEquationSeparate) mogl_glunsupported("glBlendEquationSeparate");
	glBlendEquationSeparate((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]));

}

void gl_blendfunc( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendFunc) mogl_glunsupported("glBlendFunc");
	glBlendFunc((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]));

}

void gl_calllist( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCallList) mogl_glunsupported("glCallList");
	glCallList((GLuint)mxGetScalar(prhs[0]));

}

void gl_calllists( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCallLists) mogl_glunsupported("glCallLists");
	glCallLists((GLsizei)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLvoid*)mxGetData(prhs[2]));

}

void gl_clear( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClear) mogl_glunsupported("glClear");
	glClear((GLbitfield)mxGetScalar(prhs[0]));

}

void gl_clearaccum( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearAccum) mogl_glunsupported("glClearAccum");
	glClearAccum((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_clearcolor( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearColor) mogl_glunsupported("glClearColor");
	glClearColor((GLclampf)mxGetScalar(prhs[0]),
		(GLclampf)mxGetScalar(prhs[1]),
		(GLclampf)mxGetScalar(prhs[2]),
		(GLclampf)mxGetScalar(prhs[3]));

}

void gl_cleardepth( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearDepth) mogl_glunsupported("glClearDepth");
	glClearDepth((GLclampd)mxGetScalar(prhs[0]));

}

void gl_clearindex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearIndex) mogl_glunsupported("glClearIndex");
	glClearIndex((GLfloat)mxGetScalar(prhs[0]));

}

void gl_clearstencil( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearStencil) mogl_glunsupported("glClearStencil");
	glClearStencil((GLint)mxGetScalar(prhs[0]));

}

void gl_clipplane( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClipPlane) mogl_glunsupported("glClipPlane");
	glClipPlane((GLenum)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_color3b( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3b) mogl_glunsupported("glColor3b");
	glColor3b((GLbyte)mxGetScalar(prhs[0]),
		(GLbyte)mxGetScalar(prhs[1]),
		(GLbyte)mxGetScalar(prhs[2]));

}

void gl_color3bv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3bv) mogl_glunsupported("glColor3bv");
	glColor3bv((const GLbyte*)mxGetData(prhs[0]));

}

void gl_color3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3d) mogl_glunsupported("glColor3d");
	glColor3d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_color3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3dv) mogl_glunsupported("glColor3dv");
	glColor3dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_color3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3f) mogl_glunsupported("glColor3f");
	glColor3f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_color3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3fv) mogl_glunsupported("glColor3fv");
	glColor3fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_color3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3i) mogl_glunsupported("glColor3i");
	glColor3i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_color3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3iv) mogl_glunsupported("glColor3iv");
	glColor3iv((const GLint*)mxGetData(prhs[0]));

}

void gl_color3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3s) mogl_glunsupported("glColor3s");
	glColor3s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_color3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3sv) mogl_glunsupported("glColor3sv");
	glColor3sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_color3ub( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3ub) mogl_glunsupported("glColor3ub");
	glColor3ub((GLubyte)mxGetScalar(prhs[0]),
		(GLubyte)mxGetScalar(prhs[1]),
		(GLubyte)mxGetScalar(prhs[2]));

}

void gl_color3ubv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3ubv) mogl_glunsupported("glColor3ubv");
	glColor3ubv((const GLubyte*)mxGetData(prhs[0]));

}

void gl_color3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3ui) mogl_glunsupported("glColor3ui");
	glColor3ui((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_color3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3uiv) mogl_glunsupported("glColor3uiv");
	glColor3uiv((const GLuint*)mxGetData(prhs[0]));

}

void gl_color3us( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3us) mogl_glunsupported("glColor3us");
	glColor3us((GLushort)mxGetScalar(prhs[0]),
		(GLushort)mxGetScalar(prhs[1]),
		(GLushort)mxGetScalar(prhs[2]));

}

void gl_color3usv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor3usv) mogl_glunsupported("glColor3usv");
	glColor3usv((const GLushort*)mxGetData(prhs[0]));

}

void gl_color4b( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4b) mogl_glunsupported("glColor4b");
	glColor4b((GLbyte)mxGetScalar(prhs[0]),
		(GLbyte)mxGetScalar(prhs[1]),
		(GLbyte)mxGetScalar(prhs[2]),
		(GLbyte)mxGetScalar(prhs[3]));

}

void gl_color4bv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4bv) mogl_glunsupported("glColor4bv");
	glColor4bv((const GLbyte*)mxGetData(prhs[0]));

}

void gl_color4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4d) mogl_glunsupported("glColor4d");
	glColor4d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_color4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4dv) mogl_glunsupported("glColor4dv");
	glColor4dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_color4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4f) mogl_glunsupported("glColor4f");
	glColor4f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_color4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4fv) mogl_glunsupported("glColor4fv");
	glColor4fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_color4i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4i) mogl_glunsupported("glColor4i");
	glColor4i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_color4iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4iv) mogl_glunsupported("glColor4iv");
	glColor4iv((const GLint*)mxGetData(prhs[0]));

}

void gl_color4s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4s) mogl_glunsupported("glColor4s");
	glColor4s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]));

}

void gl_color4sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4sv) mogl_glunsupported("glColor4sv");
	glColor4sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_color4ub( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4ub) mogl_glunsupported("glColor4ub");
	glColor4ub((GLubyte)mxGetScalar(prhs[0]),
		(GLubyte)mxGetScalar(prhs[1]),
		(GLubyte)mxGetScalar(prhs[2]),
		(GLubyte)mxGetScalar(prhs[3]));

}

void gl_color4ubv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4ubv) mogl_glunsupported("glColor4ubv");
	glColor4ubv((const GLubyte*)mxGetData(prhs[0]));

}

void gl_color4ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4ui) mogl_glunsupported("glColor4ui");
	glColor4ui((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_color4uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4uiv) mogl_glunsupported("glColor4uiv");
	glColor4uiv((const GLuint*)mxGetData(prhs[0]));

}

void gl_color4us( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4us) mogl_glunsupported("glColor4us");
	glColor4us((GLushort)mxGetScalar(prhs[0]),
		(GLushort)mxGetScalar(prhs[1]),
		(GLushort)mxGetScalar(prhs[2]),
		(GLushort)mxGetScalar(prhs[3]));

}

void gl_color4usv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColor4usv) mogl_glunsupported("glColor4usv");
	glColor4usv((const GLushort*)mxGetData(prhs[0]));

}

void gl_colormask( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorMask) mogl_glunsupported("glColorMask");
	glColorMask((GLboolean)mxGetScalar(prhs[0]),
		(GLboolean)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]));

}

void gl_colormaterial( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorMaterial) mogl_glunsupported("glColorMaterial");
	glColorMaterial((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]));

}

void gl_colorsubtable( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorSubTable) mogl_glunsupported("glColorSubTable");
	glColorSubTable((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(const GLvoid*)mxGetData(prhs[5]));

}

void gl_colortable( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorTable) mogl_glunsupported("glColorTable");
	glColorTable((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(const GLvoid*)mxGetData(prhs[5]));

}

void gl_colortableparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorTableParameterfv) mogl_glunsupported("glColorTableParameterfv");
	glColorTableParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_colortableparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorTableParameteriv) mogl_glunsupported("glColorTableParameteriv");
	glColorTableParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_convolutionfilter1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glConvolutionFilter1D) mogl_glunsupported("glConvolutionFilter1D");
	glConvolutionFilter1D((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(const GLvoid*)mxGetData(prhs[5]));

}

void gl_convolutionfilter2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glConvolutionFilter2D) mogl_glunsupported("glConvolutionFilter2D");
	glConvolutionFilter2D((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(const GLvoid*)mxGetData(prhs[6]));

}

void gl_convolutionparameterf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glConvolutionParameterf) mogl_glunsupported("glConvolutionParameterf");
	glConvolutionParameterf((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_convolutionparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glConvolutionParameterfv) mogl_glunsupported("glConvolutionParameterfv");
	glConvolutionParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_convolutionparameteri( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glConvolutionParameteri) mogl_glunsupported("glConvolutionParameteri");
	glConvolutionParameteri((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_convolutionparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glConvolutionParameteriv) mogl_glunsupported("glConvolutionParameteriv");
	glConvolutionParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_copycolorsubtable( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyColorSubTable) mogl_glunsupported("glCopyColorSubTable");
	glCopyColorSubTable((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_copycolortable( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyColorTable) mogl_glunsupported("glCopyColorTable");
	glCopyColorTable((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_copyconvolutionfilter1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyConvolutionFilter1D) mogl_glunsupported("glCopyConvolutionFilter1D");
	glCopyConvolutionFilter1D((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_copyconvolutionfilter2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyConvolutionFilter2D) mogl_glunsupported("glCopyConvolutionFilter2D");
	glCopyConvolutionFilter2D((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]));

}

void gl_copypixels( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyPixels) mogl_glunsupported("glCopyPixels");
	glCopyPixels((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]));

}

void gl_copyteximage1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyTexImage1D) mogl_glunsupported("glCopyTexImage1D");
	glCopyTexImage1D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLint)mxGetScalar(prhs[6]));

}

void gl_copyteximage2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyTexImage2D) mogl_glunsupported("glCopyTexImage2D");
	glCopyTexImage2D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(GLint)mxGetScalar(prhs[7]));

}

void gl_copytexsubimage1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyTexSubImage1D) mogl_glunsupported("glCopyTexSubImage1D");
	glCopyTexSubImage1D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]));

}

void gl_copytexsubimage2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyTexSubImage2D) mogl_glunsupported("glCopyTexSubImage2D");
	glCopyTexSubImage2D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(GLsizei)mxGetScalar(prhs[7]));

}

void gl_copytexsubimage3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyTexSubImage3D) mogl_glunsupported("glCopyTexSubImage3D");
	glCopyTexSubImage3D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLint)mxGetScalar(prhs[6]),
		(GLsizei)mxGetScalar(prhs[7]),
		(GLsizei)mxGetScalar(prhs[8]));

}

void gl_cullface( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCullFace) mogl_glunsupported("glCullFace");
	glCullFace((GLenum)mxGetScalar(prhs[0]));

}

void gl_deletelists( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteLists) mogl_glunsupported("glDeleteLists");
	glDeleteLists((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]));

}

void gl_deletetextures( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteTextures) mogl_glunsupported("glDeleteTextures");
	glDeleteTextures((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_depthfunc( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDepthFunc) mogl_glunsupported("glDepthFunc");
	glDepthFunc((GLenum)mxGetScalar(prhs[0]));

}

void gl_depthmask( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDepthMask) mogl_glunsupported("glDepthMask");
	glDepthMask((GLboolean)mxGetScalar(prhs[0]));

}

void gl_depthrange( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDepthRange) mogl_glunsupported("glDepthRange");
	glDepthRange((GLclampd)mxGetScalar(prhs[0]),
		(GLclampd)mxGetScalar(prhs[1]));

}

void gl_disable( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDisable) mogl_glunsupported("glDisable");
	glDisable((GLenum)mxGetScalar(prhs[0]));

}

void gl_disableclientstate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDisableClientState) mogl_glunsupported("glDisableClientState");
	glDisableClientState((GLenum)mxGetScalar(prhs[0]));

}

void gl_drawarrays( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawArrays) mogl_glunsupported("glDrawArrays");
	glDrawArrays((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]));

}

void gl_drawbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawBuffer) mogl_glunsupported("glDrawBuffer");
	glDrawBuffer((GLenum)mxGetScalar(prhs[0]));

}

void gl_drawpixels( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawPixels) mogl_glunsupported("glDrawPixels");
	glDrawPixels((GLsizei)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(const GLvoid*)mxGetData(prhs[4]));

}

void gl_edgeflag( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEdgeFlag) mogl_glunsupported("glEdgeFlag");
	glEdgeFlag((GLboolean)mxGetScalar(prhs[0]));

}

void gl_edgeflagpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEdgeFlagPointer) mogl_glunsupported("glEdgeFlagPointer");
	glEdgeFlagPointer((GLsizei)mxGetScalar(prhs[0]),
		(const GLvoid*)mxGetData(prhs[1]));

}

void gl_edgeflagv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEdgeFlagv) mogl_glunsupported("glEdgeFlagv");
	glEdgeFlagv((const GLboolean*)mxGetData(prhs[0]));

}

void gl_enable( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEnable) mogl_glunsupported("glEnable");
	glEnable((GLenum)mxGetScalar(prhs[0]));

}

void gl_enableclientstate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEnableClientState) mogl_glunsupported("glEnableClientState");
	glEnableClientState((GLenum)mxGetScalar(prhs[0]));

}

void gl_end( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEnd) mogl_glunsupported("glEnd");
	glEnd();

}

void gl_endlist( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEndList) mogl_glunsupported("glEndList");
	glEndList();

}

void gl_evalcoord1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalCoord1d) mogl_glunsupported("glEvalCoord1d");
	glEvalCoord1d((GLdouble)mxGetScalar(prhs[0]));

}

void gl_evalcoord1dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalCoord1dv) mogl_glunsupported("glEvalCoord1dv");
	glEvalCoord1dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_evalcoord1f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalCoord1f) mogl_glunsupported("glEvalCoord1f");
	glEvalCoord1f((GLfloat)mxGetScalar(prhs[0]));

}

void gl_evalcoord1fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalCoord1fv) mogl_glunsupported("glEvalCoord1fv");
	glEvalCoord1fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_evalcoord2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalCoord2d) mogl_glunsupported("glEvalCoord2d");
	glEvalCoord2d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_evalcoord2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalCoord2dv) mogl_glunsupported("glEvalCoord2dv");
	glEvalCoord2dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_evalcoord2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalCoord2f) mogl_glunsupported("glEvalCoord2f");
	glEvalCoord2f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_evalcoord2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalCoord2fv) mogl_glunsupported("glEvalCoord2fv");
	glEvalCoord2fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_evalmesh1( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalMesh1) mogl_glunsupported("glEvalMesh1");
	glEvalMesh1((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_evalmesh2( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalMesh2) mogl_glunsupported("glEvalMesh2");
	glEvalMesh2((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_evalpoint1( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalPoint1) mogl_glunsupported("glEvalPoint1");
	glEvalPoint1((GLint)mxGetScalar(prhs[0]));

}

void gl_evalpoint2( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEvalPoint2) mogl_glunsupported("glEvalPoint2");
	glEvalPoint2((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_finish( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFinish) mogl_glunsupported("glFinish");
	glFinish();

}

void gl_flush( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFlush) mogl_glunsupported("glFlush");
	glFlush();

}

void gl_fogf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogf) mogl_glunsupported("glFogf");
	glFogf((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_fogfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogfv) mogl_glunsupported("glFogfv");
	glFogfv((GLenum)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_fogi( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogi) mogl_glunsupported("glFogi");
	glFogi((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_fogiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogiv) mogl_glunsupported("glFogiv");
	glFogiv((GLenum)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_frontface( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFrontFace) mogl_glunsupported("glFrontFace");
	glFrontFace((GLenum)mxGetScalar(prhs[0]));

}

void gl_frustum( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFrustum) mogl_glunsupported("glFrustum");
	glFrustum((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]),
		(GLdouble)mxGetScalar(prhs[5]));

}

void gl_genlists( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenLists) mogl_glunsupported("glGenLists");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGenLists((GLsizei)mxGetScalar(prhs[0]));

}

void gl_gentextures( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenTextures) mogl_glunsupported("glGenTextures");
	glGenTextures((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_getbooleanv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetBooleanv) mogl_glunsupported("glGetBooleanv");
	glGetBooleanv((GLenum)mxGetScalar(prhs[0]),
		(GLboolean*)mxGetData(prhs[1]));

}

void gl_getclipplane( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetClipPlane) mogl_glunsupported("glGetClipPlane");
	glGetClipPlane((GLenum)mxGetScalar(prhs[0]),
		(GLdouble*)mxGetData(prhs[1]));

}

void gl_getcolortable( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetColorTable) mogl_glunsupported("glGetColorTable");
	glGetColorTable((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLvoid*)mxGetData(prhs[3]));

}

void gl_getcolortableparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetColorTableParameterfv) mogl_glunsupported("glGetColorTableParameterfv");
	glGetColorTableParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getcolortableparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetColorTableParameteriv) mogl_glunsupported("glGetColorTableParameteriv");
	glGetColorTableParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getconvolutionfilter( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetConvolutionFilter) mogl_glunsupported("glGetConvolutionFilter");
	glGetConvolutionFilter((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLvoid*)mxGetData(prhs[3]));

}

void gl_getconvolutionparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetConvolutionParameterfv) mogl_glunsupported("glGetConvolutionParameterfv");
	glGetConvolutionParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getconvolutionparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetConvolutionParameteriv) mogl_glunsupported("glGetConvolutionParameteriv");
	glGetConvolutionParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getdoublev( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetDoublev) mogl_glunsupported("glGetDoublev");
	glGetDoublev((GLenum)mxGetScalar(prhs[0]),
		(GLdouble*)mxGetData(prhs[1]));

}

void gl_geterror( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetError) mogl_glunsupported("glGetError");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetError();

}

void gl_getfloatv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetFloatv) mogl_glunsupported("glGetFloatv");
	glGetFloatv((GLenum)mxGetScalar(prhs[0]),
		(GLfloat*)mxGetData(prhs[1]));

}

void gl_gethistogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetHistogram) mogl_glunsupported("glGetHistogram");
	glGetHistogram((GLenum)mxGetScalar(prhs[0]),
		(GLboolean)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLvoid*)mxGetData(prhs[4]));

}

void gl_gethistogramparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetHistogramParameterfv) mogl_glunsupported("glGetHistogramParameterfv");
	glGetHistogramParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_gethistogramparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetHistogramParameteriv) mogl_glunsupported("glGetHistogramParameteriv");
	glGetHistogramParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getintegerv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetIntegerv) mogl_glunsupported("glGetIntegerv");
	glGetIntegerv((GLenum)mxGetScalar(prhs[0]),
		(GLint*)mxGetData(prhs[1]));

}

void gl_getlightfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetLightfv) mogl_glunsupported("glGetLightfv");
	glGetLightfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getlightiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetLightiv) mogl_glunsupported("glGetLightiv");
	glGetLightiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getmapdv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMapdv) mogl_glunsupported("glGetMapdv");
	glGetMapdv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLdouble*)mxGetData(prhs[2]));

}

void gl_getmapfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMapfv) mogl_glunsupported("glGetMapfv");
	glGetMapfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getmapiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMapiv) mogl_glunsupported("glGetMapiv");
	glGetMapiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getmaterialfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMaterialfv) mogl_glunsupported("glGetMaterialfv");
	glGetMaterialfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getmaterialiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMaterialiv) mogl_glunsupported("glGetMaterialiv");
	glGetMaterialiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getminmax( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMinmax) mogl_glunsupported("glGetMinmax");
	glGetMinmax((GLenum)mxGetScalar(prhs[0]),
		(GLboolean)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLvoid*)mxGetData(prhs[4]));

}

void gl_getminmaxparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMinmaxParameterfv) mogl_glunsupported("glGetMinmaxParameterfv");
	glGetMinmaxParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getminmaxparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMinmaxParameteriv) mogl_glunsupported("glGetMinmaxParameteriv");
	glGetMinmaxParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getpixelmapfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetPixelMapfv) mogl_glunsupported("glGetPixelMapfv");
	glGetPixelMapfv((GLenum)mxGetScalar(prhs[0]),
		(GLfloat*)mxGetData(prhs[1]));

}

void gl_getpixelmapuiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetPixelMapuiv) mogl_glunsupported("glGetPixelMapuiv");
	glGetPixelMapuiv((GLenum)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_getpixelmapusv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetPixelMapusv) mogl_glunsupported("glGetPixelMapusv");
	glGetPixelMapusv((GLenum)mxGetScalar(prhs[0]),
		(GLushort*)mxGetData(prhs[1]));

}

void gl_getpolygonstipple( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetPolygonStipple) mogl_glunsupported("glGetPolygonStipple");
	glGetPolygonStipple((GLubyte*)mxGetData(prhs[0]));

}

void gl_getseparablefilter( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetSeparableFilter) mogl_glunsupported("glGetSeparableFilter");
	glGetSeparableFilter((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLvoid*)mxGetData(prhs[3]),
		(GLvoid*)mxGetData(prhs[4]),
		(GLvoid*)mxGetData(prhs[5]));

}

void gl_gettexenvfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexEnvfv) mogl_glunsupported("glGetTexEnvfv");
	glGetTexEnvfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_gettexenviv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexEnviv) mogl_glunsupported("glGetTexEnviv");
	glGetTexEnviv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_gettexgendv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexGendv) mogl_glunsupported("glGetTexGendv");
	glGetTexGendv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLdouble*)mxGetData(prhs[2]));

}

void gl_gettexgenfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexGenfv) mogl_glunsupported("glGetTexGenfv");
	glGetTexGenfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_gettexgeniv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexGeniv) mogl_glunsupported("glGetTexGeniv");
	glGetTexGeniv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getteximage( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexImage) mogl_glunsupported("glGetTexImage");
	glGetTexImage((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLvoid*)mxGetData(prhs[4]));

}

void gl_gettexlevelparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexLevelParameterfv) mogl_glunsupported("glGetTexLevelParameterfv");
	glGetTexLevelParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLfloat*)mxGetData(prhs[3]));

}

void gl_gettexlevelparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexLevelParameteriv) mogl_glunsupported("glGetTexLevelParameteriv");
	glGetTexLevelParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_gettexparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexParameterfv) mogl_glunsupported("glGetTexParameterfv");
	glGetTexParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_gettexparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetTexParameteriv) mogl_glunsupported("glGetTexParameteriv");
	glGetTexParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_hint( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glHint) mogl_glunsupported("glHint");
	glHint((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]));

}

void gl_histogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glHistogram) mogl_glunsupported("glHistogram");
	glHistogram((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]));

}

void gl_indexmask( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexMask) mogl_glunsupported("glIndexMask");
	glIndexMask((GLuint)mxGetScalar(prhs[0]));

}

void gl_indexpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexPointer) mogl_glunsupported("glIndexPointer");
	glIndexPointer((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLvoid*)mxGetData(prhs[2]));

}

void gl_indexd( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexd) mogl_glunsupported("glIndexd");
	glIndexd((GLdouble)mxGetScalar(prhs[0]));

}

void gl_indexdv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexdv) mogl_glunsupported("glIndexdv");
	glIndexdv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_indexf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexf) mogl_glunsupported("glIndexf");
	glIndexf((GLfloat)mxGetScalar(prhs[0]));

}

void gl_indexfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexfv) mogl_glunsupported("glIndexfv");
	glIndexfv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_indexi( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexi) mogl_glunsupported("glIndexi");
	glIndexi((GLint)mxGetScalar(prhs[0]));

}

void gl_indexiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexiv) mogl_glunsupported("glIndexiv");
	glIndexiv((const GLint*)mxGetData(prhs[0]));

}

void gl_indexs( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexs) mogl_glunsupported("glIndexs");
	glIndexs((GLshort)mxGetScalar(prhs[0]));

}

void gl_indexsv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexsv) mogl_glunsupported("glIndexsv");
	glIndexsv((const GLshort*)mxGetData(prhs[0]));

}

void gl_indexub( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexub) mogl_glunsupported("glIndexub");
	glIndexub((GLubyte)mxGetScalar(prhs[0]));

}

void gl_indexubv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIndexubv) mogl_glunsupported("glIndexubv");
	glIndexubv((const GLubyte*)mxGetData(prhs[0]));

}

void gl_initnames( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glInitNames) mogl_glunsupported("glInitNames");
	glInitNames();

}

void gl_interleavedarrays( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glInterleavedArrays) mogl_glunsupported("glInterleavedArrays");
	glInterleavedArrays((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLvoid*)mxGetData(prhs[2]));

}

void gl_isenabled( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsEnabled) mogl_glunsupported("glIsEnabled");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsEnabled((GLenum)mxGetScalar(prhs[0]));

}

void gl_islist( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsList) mogl_glunsupported("glIsList");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsList((GLuint)mxGetScalar(prhs[0]));

}

void gl_istexture( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsTexture) mogl_glunsupported("glIsTexture");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsTexture((GLuint)mxGetScalar(prhs[0]));

}

void gl_lightmodelf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLightModelf) mogl_glunsupported("glLightModelf");
	glLightModelf((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_lightmodelfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLightModelfv) mogl_glunsupported("glLightModelfv");
	glLightModelfv((GLenum)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_lightmodeli( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLightModeli) mogl_glunsupported("glLightModeli");
	glLightModeli((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_lightmodeliv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLightModeliv) mogl_glunsupported("glLightModeliv");
	glLightModeliv((GLenum)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_lightf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLightf) mogl_glunsupported("glLightf");
	glLightf((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_lightfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLightfv) mogl_glunsupported("glLightfv");
	glLightfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_lighti( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLighti) mogl_glunsupported("glLighti");
	glLighti((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_lightiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLightiv) mogl_glunsupported("glLightiv");
	glLightiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_linestipple( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLineStipple) mogl_glunsupported("glLineStipple");
	glLineStipple((GLint)mxGetScalar(prhs[0]),
		(GLushort)mxGetScalar(prhs[1]));

}

void gl_linewidth( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLineWidth) mogl_glunsupported("glLineWidth");
	glLineWidth((GLfloat)mxGetScalar(prhs[0]));

}

void gl_listbase( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glListBase) mogl_glunsupported("glListBase");
	glListBase((GLuint)mxGetScalar(prhs[0]));

}

void gl_loadidentity( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLoadIdentity) mogl_glunsupported("glLoadIdentity");
	glLoadIdentity();

}

void gl_loadmatrixd( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLoadMatrixd) mogl_glunsupported("glLoadMatrixd");
	glLoadMatrixd((const GLdouble*)mxGetData(prhs[0]));

}

void gl_loadmatrixf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLoadMatrixf) mogl_glunsupported("glLoadMatrixf");
	glLoadMatrixf((const GLfloat*)mxGetData(prhs[0]));

}

void gl_loadname( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLoadName) mogl_glunsupported("glLoadName");
	glLoadName((GLuint)mxGetScalar(prhs[0]));

}

void gl_logicop( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLogicOp) mogl_glunsupported("glLogicOp");
	glLogicOp((GLenum)mxGetScalar(prhs[0]));

}

void gl_map1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMap1d) mogl_glunsupported("glMap1d");
	glMap1d((GLenum)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(const GLdouble*)mxGetData(prhs[5]));

}

void gl_map1f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMap1f) mogl_glunsupported("glMap1f");
	glMap1f((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(const GLfloat*)mxGetData(prhs[5]));

}

void gl_map2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMap2d) mogl_glunsupported("glMap2d");
	glMap2d((GLenum)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLdouble)mxGetScalar(prhs[5]),
		(GLdouble)mxGetScalar(prhs[6]),
		(GLint)mxGetScalar(prhs[7]),
		(GLint)mxGetScalar(prhs[8]),
		(const GLdouble*)mxGetData(prhs[9]));

}

void gl_map2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMap2f) mogl_glunsupported("glMap2f");
	glMap2f((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLfloat)mxGetScalar(prhs[5]),
		(GLfloat)mxGetScalar(prhs[6]),
		(GLint)mxGetScalar(prhs[7]),
		(GLint)mxGetScalar(prhs[8]),
		(const GLfloat*)mxGetData(prhs[9]));

}

void gl_mapgrid1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMapGrid1d) mogl_glunsupported("glMapGrid1d");
	glMapGrid1d((GLint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_mapgrid1f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMapGrid1f) mogl_glunsupported("glMapGrid1f");
	glMapGrid1f((GLint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_mapgrid2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMapGrid2d) mogl_glunsupported("glMapGrid2d");
	glMapGrid2d((GLint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]),
		(GLdouble)mxGetScalar(prhs[5]));

}

void gl_mapgrid2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMapGrid2f) mogl_glunsupported("glMapGrid2f");
	glMapGrid2f((GLint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLfloat)mxGetScalar(prhs[4]),
		(GLfloat)mxGetScalar(prhs[5]));

}

void gl_materialf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMaterialf) mogl_glunsupported("glMaterialf");
	glMaterialf((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_materialfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMaterialfv) mogl_glunsupported("glMaterialfv");
	glMaterialfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_materiali( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMateriali) mogl_glunsupported("glMateriali");
	glMateriali((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_materialiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMaterialiv) mogl_glunsupported("glMaterialiv");
	glMaterialiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_matrixmode( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMatrixMode) mogl_glunsupported("glMatrixMode");
	glMatrixMode((GLenum)mxGetScalar(prhs[0]));

}

void gl_minmax( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMinmax) mogl_glunsupported("glMinmax");
	glMinmax((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]));

}

void gl_multmatrixd( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultMatrixd) mogl_glunsupported("glMultMatrixd");
	glMultMatrixd((const GLdouble*)mxGetData(prhs[0]));

}

void gl_multmatrixf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultMatrixf) mogl_glunsupported("glMultMatrixf");
	glMultMatrixf((const GLfloat*)mxGetData(prhs[0]));

}

void gl_newlist( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNewList) mogl_glunsupported("glNewList");
	glNewList((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]));

}

void gl_normal3b( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3b) mogl_glunsupported("glNormal3b");
	glNormal3b((GLbyte)mxGetScalar(prhs[0]),
		(GLbyte)mxGetScalar(prhs[1]),
		(GLbyte)mxGetScalar(prhs[2]));

}

void gl_normal3bv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3bv) mogl_glunsupported("glNormal3bv");
	glNormal3bv((const GLbyte*)mxGetData(prhs[0]));

}

void gl_normal3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3d) mogl_glunsupported("glNormal3d");
	glNormal3d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_normal3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3dv) mogl_glunsupported("glNormal3dv");
	glNormal3dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_normal3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3f) mogl_glunsupported("glNormal3f");
	glNormal3f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_normal3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3fv) mogl_glunsupported("glNormal3fv");
	glNormal3fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_normal3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3i) mogl_glunsupported("glNormal3i");
	glNormal3i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_normal3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3iv) mogl_glunsupported("glNormal3iv");
	glNormal3iv((const GLint*)mxGetData(prhs[0]));

}

void gl_normal3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3s) mogl_glunsupported("glNormal3s");
	glNormal3s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_normal3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormal3sv) mogl_glunsupported("glNormal3sv");
	glNormal3sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_ortho( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glOrtho) mogl_glunsupported("glOrtho");
	glOrtho((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]),
		(GLdouble)mxGetScalar(prhs[5]));

}

void gl_passthrough( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPassThrough) mogl_glunsupported("glPassThrough");
	glPassThrough((GLfloat)mxGetScalar(prhs[0]));

}

void gl_pixelmapfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPixelMapfv) mogl_glunsupported("glPixelMapfv");
	glPixelMapfv((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_pixelmapuiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPixelMapuiv) mogl_glunsupported("glPixelMapuiv");
	glPixelMapuiv((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(const GLuint*)mxGetData(prhs[2]));

}

void gl_pixelmapusv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPixelMapusv) mogl_glunsupported("glPixelMapusv");
	glPixelMapusv((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(const GLushort*)mxGetData(prhs[2]));

}

void gl_pixelstoref( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPixelStoref) mogl_glunsupported("glPixelStoref");
	glPixelStoref((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_pixelstorei( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPixelStorei) mogl_glunsupported("glPixelStorei");
	glPixelStorei((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_pixeltransferf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPixelTransferf) mogl_glunsupported("glPixelTransferf");
	glPixelTransferf((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_pixeltransferi( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPixelTransferi) mogl_glunsupported("glPixelTransferi");
	glPixelTransferi((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_pixelzoom( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPixelZoom) mogl_glunsupported("glPixelZoom");
	glPixelZoom((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_pointsize( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPointSize) mogl_glunsupported("glPointSize");
	glPointSize((GLfloat)mxGetScalar(prhs[0]));

}

void gl_polygonmode( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPolygonMode) mogl_glunsupported("glPolygonMode");
	glPolygonMode((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]));

}

void gl_polygonoffset( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPolygonOffset) mogl_glunsupported("glPolygonOffset");
	glPolygonOffset((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_polygonstipple( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPolygonStipple) mogl_glunsupported("glPolygonStipple");
	glPolygonStipple((const GLubyte*)mxGetData(prhs[0]));

}

void gl_popattrib( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPopAttrib) mogl_glunsupported("glPopAttrib");
	glPopAttrib();

}

void gl_popclientattrib( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPopClientAttrib) mogl_glunsupported("glPopClientAttrib");
	glPopClientAttrib();

}

void gl_popmatrix( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPopMatrix) mogl_glunsupported("glPopMatrix");
	glPopMatrix();

}

void gl_popname( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPopName) mogl_glunsupported("glPopName");
	glPopName();

}

void gl_prioritizetextures( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPrioritizeTextures) mogl_glunsupported("glPrioritizeTextures");
	glPrioritizeTextures((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]),
		(const GLclampf*)mxGetData(prhs[2]));

}

void gl_pushattrib( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPushAttrib) mogl_glunsupported("glPushAttrib");
	glPushAttrib((GLbitfield)mxGetScalar(prhs[0]));

}

void gl_pushclientattrib( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPushClientAttrib) mogl_glunsupported("glPushClientAttrib");
	glPushClientAttrib((GLbitfield)mxGetScalar(prhs[0]));

}

void gl_pushmatrix( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPushMatrix) mogl_glunsupported("glPushMatrix");
	glPushMatrix();

}

void gl_pushname( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPushName) mogl_glunsupported("glPushName");
	glPushName((GLuint)mxGetScalar(prhs[0]));

}

void gl_rasterpos2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos2d) mogl_glunsupported("glRasterPos2d");
	glRasterPos2d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_rasterpos2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos2dv) mogl_glunsupported("glRasterPos2dv");
	glRasterPos2dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_rasterpos2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos2f) mogl_glunsupported("glRasterPos2f");
	glRasterPos2f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_rasterpos2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos2fv) mogl_glunsupported("glRasterPos2fv");
	glRasterPos2fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_rasterpos2i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos2i) mogl_glunsupported("glRasterPos2i");
	glRasterPos2i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_rasterpos2iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos2iv) mogl_glunsupported("glRasterPos2iv");
	glRasterPos2iv((const GLint*)mxGetData(prhs[0]));

}

void gl_rasterpos2s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos2s) mogl_glunsupported("glRasterPos2s");
	glRasterPos2s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]));

}

void gl_rasterpos2sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos2sv) mogl_glunsupported("glRasterPos2sv");
	glRasterPos2sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_rasterpos3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos3d) mogl_glunsupported("glRasterPos3d");
	glRasterPos3d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_rasterpos3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos3dv) mogl_glunsupported("glRasterPos3dv");
	glRasterPos3dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_rasterpos3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos3f) mogl_glunsupported("glRasterPos3f");
	glRasterPos3f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_rasterpos3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos3fv) mogl_glunsupported("glRasterPos3fv");
	glRasterPos3fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_rasterpos3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos3i) mogl_glunsupported("glRasterPos3i");
	glRasterPos3i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_rasterpos3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos3iv) mogl_glunsupported("glRasterPos3iv");
	glRasterPos3iv((const GLint*)mxGetData(prhs[0]));

}

void gl_rasterpos3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos3s) mogl_glunsupported("glRasterPos3s");
	glRasterPos3s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_rasterpos3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos3sv) mogl_glunsupported("glRasterPos3sv");
	glRasterPos3sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_rasterpos4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos4d) mogl_glunsupported("glRasterPos4d");
	glRasterPos4d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_rasterpos4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos4dv) mogl_glunsupported("glRasterPos4dv");
	glRasterPos4dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_rasterpos4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos4f) mogl_glunsupported("glRasterPos4f");
	glRasterPos4f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_rasterpos4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos4fv) mogl_glunsupported("glRasterPos4fv");
	glRasterPos4fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_rasterpos4i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos4i) mogl_glunsupported("glRasterPos4i");
	glRasterPos4i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_rasterpos4iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos4iv) mogl_glunsupported("glRasterPos4iv");
	glRasterPos4iv((const GLint*)mxGetData(prhs[0]));

}

void gl_rasterpos4s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos4s) mogl_glunsupported("glRasterPos4s");
	glRasterPos4s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]));

}

void gl_rasterpos4sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRasterPos4sv) mogl_glunsupported("glRasterPos4sv");
	glRasterPos4sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_readbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glReadBuffer) mogl_glunsupported("glReadBuffer");
	glReadBuffer((GLenum)mxGetScalar(prhs[0]));

}

void gl_rectd( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRectd) mogl_glunsupported("glRectd");
	glRectd((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_rectdv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRectdv) mogl_glunsupported("glRectdv");
	glRectdv((const GLdouble*)mxGetData(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_rectf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRectf) mogl_glunsupported("glRectf");
	glRectf((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_rectfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRectfv) mogl_glunsupported("glRectfv");
	glRectfv((const GLfloat*)mxGetData(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_recti( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRecti) mogl_glunsupported("glRecti");
	glRecti((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_rectiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRectiv) mogl_glunsupported("glRectiv");
	glRectiv((const GLint*)mxGetData(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_rects( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRects) mogl_glunsupported("glRects");
	glRects((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]));

}

void gl_rectsv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRectsv) mogl_glunsupported("glRectsv");
	glRectsv((const GLshort*)mxGetData(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_rendermode( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRenderMode) mogl_glunsupported("glRenderMode");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glRenderMode((GLenum)mxGetScalar(prhs[0]));

}

void gl_resethistogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glResetHistogram) mogl_glunsupported("glResetHistogram");
	glResetHistogram((GLenum)mxGetScalar(prhs[0]));

}

void gl_resetminmax( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glResetMinmax) mogl_glunsupported("glResetMinmax");
	glResetMinmax((GLenum)mxGetScalar(prhs[0]));

}

void gl_rotated( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRotated) mogl_glunsupported("glRotated");
	glRotated((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_rotatef( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRotatef) mogl_glunsupported("glRotatef");
	glRotatef((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_scaled( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glScaled) mogl_glunsupported("glScaled");
	glScaled((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_scalef( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glScalef) mogl_glunsupported("glScalef");
	glScalef((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_scissor( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glScissor) mogl_glunsupported("glScissor");
	glScissor((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_separablefilter2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSeparableFilter2D) mogl_glunsupported("glSeparableFilter2D");
	glSeparableFilter2D((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(const GLvoid*)mxGetData(prhs[6]),
		(const GLvoid*)mxGetData(prhs[7]));

}

void gl_shademodel( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glShadeModel) mogl_glunsupported("glShadeModel");
	glShadeModel((GLenum)mxGetScalar(prhs[0]));

}

void gl_stencilfunc( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glStencilFunc) mogl_glunsupported("glStencilFunc");
	glStencilFunc((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_stencilmask( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glStencilMask) mogl_glunsupported("glStencilMask");
	glStencilMask((GLuint)mxGetScalar(prhs[0]));

}

void gl_stencilop( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glStencilOp) mogl_glunsupported("glStencilOp");
	glStencilOp((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]));

}

void gl_texcoord1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord1d) mogl_glunsupported("glTexCoord1d");
	glTexCoord1d((GLdouble)mxGetScalar(prhs[0]));

}

void gl_texcoord1dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord1dv) mogl_glunsupported("glTexCoord1dv");
	glTexCoord1dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_texcoord1f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord1f) mogl_glunsupported("glTexCoord1f");
	glTexCoord1f((GLfloat)mxGetScalar(prhs[0]));

}

void gl_texcoord1fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord1fv) mogl_glunsupported("glTexCoord1fv");
	glTexCoord1fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_texcoord1i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord1i) mogl_glunsupported("glTexCoord1i");
	glTexCoord1i((GLint)mxGetScalar(prhs[0]));

}

void gl_texcoord1iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord1iv) mogl_glunsupported("glTexCoord1iv");
	glTexCoord1iv((const GLint*)mxGetData(prhs[0]));

}

void gl_texcoord1s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord1s) mogl_glunsupported("glTexCoord1s");
	glTexCoord1s((GLshort)mxGetScalar(prhs[0]));

}

void gl_texcoord1sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord1sv) mogl_glunsupported("glTexCoord1sv");
	glTexCoord1sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_texcoord2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord2d) mogl_glunsupported("glTexCoord2d");
	glTexCoord2d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_texcoord2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord2dv) mogl_glunsupported("glTexCoord2dv");
	glTexCoord2dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_texcoord2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord2f) mogl_glunsupported("glTexCoord2f");
	glTexCoord2f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_texcoord2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord2fv) mogl_glunsupported("glTexCoord2fv");
	glTexCoord2fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_texcoord2i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord2i) mogl_glunsupported("glTexCoord2i");
	glTexCoord2i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_texcoord2iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord2iv) mogl_glunsupported("glTexCoord2iv");
	glTexCoord2iv((const GLint*)mxGetData(prhs[0]));

}

void gl_texcoord2s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord2s) mogl_glunsupported("glTexCoord2s");
	glTexCoord2s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]));

}

void gl_texcoord2sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord2sv) mogl_glunsupported("glTexCoord2sv");
	glTexCoord2sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_texcoord3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord3d) mogl_glunsupported("glTexCoord3d");
	glTexCoord3d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_texcoord3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord3dv) mogl_glunsupported("glTexCoord3dv");
	glTexCoord3dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_texcoord3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord3f) mogl_glunsupported("glTexCoord3f");
	glTexCoord3f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_texcoord3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord3fv) mogl_glunsupported("glTexCoord3fv");
	glTexCoord3fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_texcoord3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord3i) mogl_glunsupported("glTexCoord3i");
	glTexCoord3i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_texcoord3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord3iv) mogl_glunsupported("glTexCoord3iv");
	glTexCoord3iv((const GLint*)mxGetData(prhs[0]));

}

void gl_texcoord3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord3s) mogl_glunsupported("glTexCoord3s");
	glTexCoord3s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_texcoord3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord3sv) mogl_glunsupported("glTexCoord3sv");
	glTexCoord3sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_texcoord4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord4d) mogl_glunsupported("glTexCoord4d");
	glTexCoord4d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_texcoord4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord4dv) mogl_glunsupported("glTexCoord4dv");
	glTexCoord4dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_texcoord4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord4f) mogl_glunsupported("glTexCoord4f");
	glTexCoord4f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_texcoord4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord4fv) mogl_glunsupported("glTexCoord4fv");
	glTexCoord4fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_texcoord4i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord4i) mogl_glunsupported("glTexCoord4i");
	glTexCoord4i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_texcoord4iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord4iv) mogl_glunsupported("glTexCoord4iv");
	glTexCoord4iv((const GLint*)mxGetData(prhs[0]));

}

void gl_texcoord4s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord4s) mogl_glunsupported("glTexCoord4s");
	glTexCoord4s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]));

}

void gl_texcoord4sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoord4sv) mogl_glunsupported("glTexCoord4sv");
	glTexCoord4sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_texenvf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexEnvf) mogl_glunsupported("glTexEnvf");
	glTexEnvf((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_texenvfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexEnvfv) mogl_glunsupported("glTexEnvfv");
	glTexEnvfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_texenvi( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexEnvi) mogl_glunsupported("glTexEnvi");
	glTexEnvi((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_texenviv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexEnviv) mogl_glunsupported("glTexEnviv");
	glTexEnviv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_texgend( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexGend) mogl_glunsupported("glTexGend");
	glTexGend((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_texgendv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexGendv) mogl_glunsupported("glTexGendv");
	glTexGendv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLdouble*)mxGetData(prhs[2]));

}

void gl_texgenf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexGenf) mogl_glunsupported("glTexGenf");
	glTexGenf((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_texgenfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexGenfv) mogl_glunsupported("glTexGenfv");
	glTexGenfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_texgeni( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexGeni) mogl_glunsupported("glTexGeni");
	glTexGeni((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_texgeniv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexGeniv) mogl_glunsupported("glTexGeniv");
	glTexGeniv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_teximage1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexImage1D) mogl_glunsupported("glTexImage1D");
	glTexImage1D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(GLenum)mxGetScalar(prhs[6]),
		(const GLvoid*)mxGetData(prhs[7]));

}

void gl_teximage3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexImage3D) mogl_glunsupported("glTexImage3D");
	glTexImage3D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLint)mxGetScalar(prhs[6]),
		(GLenum)mxGetScalar(prhs[7]),
		(GLenum)mxGetScalar(prhs[8]),
		(const GLvoid*)mxGetData(prhs[9]));

}

void gl_texparameterf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexParameterf) mogl_glunsupported("glTexParameterf");
	glTexParameterf((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_texparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexParameterfv) mogl_glunsupported("glTexParameterfv");
	glTexParameterfv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_texparameteri( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexParameteri) mogl_glunsupported("glTexParameteri");
	glTexParameteri((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_texparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexParameteriv) mogl_glunsupported("glTexParameteriv");
	glTexParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_texsubimage1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexSubImage1D) mogl_glunsupported("glTexSubImage1D");
	glTexSubImage1D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(const GLvoid*)mxGetData(prhs[6]));

}

void gl_texsubimage2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexSubImage2D) mogl_glunsupported("glTexSubImage2D");
	glTexSubImage2D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLenum)mxGetScalar(prhs[6]),
		(GLenum)mxGetScalar(prhs[7]),
		(const GLvoid*)mxGetData(prhs[8]));

}

void gl_texsubimage3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexSubImage3D) mogl_glunsupported("glTexSubImage3D");
	glTexSubImage3D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(GLsizei)mxGetScalar(prhs[7]),
		(GLenum)mxGetScalar(prhs[8]),
		(GLenum)mxGetScalar(prhs[9]),
		(const GLvoid*)mxGetData(prhs[10]));

}

void gl_translated( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTranslated) mogl_glunsupported("glTranslated");
	glTranslated((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_translatef( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTranslatef) mogl_glunsupported("glTranslatef");
	glTranslatef((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_vertex2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex2d) mogl_glunsupported("glVertex2d");
	glVertex2d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_vertex2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex2dv) mogl_glunsupported("glVertex2dv");
	glVertex2dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_vertex2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex2f) mogl_glunsupported("glVertex2f");
	glVertex2f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_vertex2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex2fv) mogl_glunsupported("glVertex2fv");
	glVertex2fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_vertex2i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex2i) mogl_glunsupported("glVertex2i");
	glVertex2i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_vertex2iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex2iv) mogl_glunsupported("glVertex2iv");
	glVertex2iv((const GLint*)mxGetData(prhs[0]));

}

void gl_vertex2s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex2s) mogl_glunsupported("glVertex2s");
	glVertex2s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]));

}

void gl_vertex2sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex2sv) mogl_glunsupported("glVertex2sv");
	glVertex2sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_vertex3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex3d) mogl_glunsupported("glVertex3d");
	glVertex3d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_vertex3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex3dv) mogl_glunsupported("glVertex3dv");
	glVertex3dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_vertex3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex3f) mogl_glunsupported("glVertex3f");
	glVertex3f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_vertex3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex3fv) mogl_glunsupported("glVertex3fv");
	glVertex3fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_vertex3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex3i) mogl_glunsupported("glVertex3i");
	glVertex3i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_vertex3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex3iv) mogl_glunsupported("glVertex3iv");
	glVertex3iv((const GLint*)mxGetData(prhs[0]));

}

void gl_vertex3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex3s) mogl_glunsupported("glVertex3s");
	glVertex3s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_vertex3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex3sv) mogl_glunsupported("glVertex3sv");
	glVertex3sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_vertex4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex4d) mogl_glunsupported("glVertex4d");
	glVertex4d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_vertex4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex4dv) mogl_glunsupported("glVertex4dv");
	glVertex4dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_vertex4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex4f) mogl_glunsupported("glVertex4f");
	glVertex4f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_vertex4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex4fv) mogl_glunsupported("glVertex4fv");
	glVertex4fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_vertex4i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex4i) mogl_glunsupported("glVertex4i");
	glVertex4i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_vertex4iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex4iv) mogl_glunsupported("glVertex4iv");
	glVertex4iv((const GLint*)mxGetData(prhs[0]));

}

void gl_vertex4s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex4s) mogl_glunsupported("glVertex4s");
	glVertex4s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]));

}

void gl_vertex4sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertex4sv) mogl_glunsupported("glVertex4sv");
	glVertex4sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_viewport( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glViewport) mogl_glunsupported("glViewport");
	glViewport((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_samplecoverage( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSampleCoverage) mogl_glunsupported("glSampleCoverage");
	glSampleCoverage((GLclampf)mxGetScalar(prhs[0]),
		(GLboolean)mxGetScalar(prhs[1]));

}

void gl_loadtransposematrixf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLoadTransposeMatrixf) mogl_glunsupported("glLoadTransposeMatrixf");
	glLoadTransposeMatrixf((const GLfloat*)mxGetData(prhs[0]));

}

void gl_loadtransposematrixd( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLoadTransposeMatrixd) mogl_glunsupported("glLoadTransposeMatrixd");
	glLoadTransposeMatrixd((const GLdouble*)mxGetData(prhs[0]));

}

void gl_multtransposematrixf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultTransposeMatrixf) mogl_glunsupported("glMultTransposeMatrixf");
	glMultTransposeMatrixf((const GLfloat*)mxGetData(prhs[0]));

}

void gl_multtransposematrixd( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultTransposeMatrixd) mogl_glunsupported("glMultTransposeMatrixd");
	glMultTransposeMatrixd((const GLdouble*)mxGetData(prhs[0]));

}

void gl_compressedteximage3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCompressedTexImage3D) mogl_glunsupported("glCompressedTexImage3D");
	glCompressedTexImage3D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLint)mxGetScalar(prhs[6]),
		(GLsizei)mxGetScalar(prhs[7]),
		(const GLvoid*)mxGetData(prhs[8]));

}

void gl_compressedteximage2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCompressedTexImage2D) mogl_glunsupported("glCompressedTexImage2D");
	glCompressedTexImage2D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(const GLvoid*)mxGetData(prhs[7]));

}

void gl_compressedteximage1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCompressedTexImage1D) mogl_glunsupported("glCompressedTexImage1D");
	glCompressedTexImage1D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(const GLvoid*)mxGetData(prhs[6]));

}

void gl_compressedtexsubimage3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCompressedTexSubImage3D) mogl_glunsupported("glCompressedTexSubImage3D");
	glCompressedTexSubImage3D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(GLsizei)mxGetScalar(prhs[7]),
		(GLenum)mxGetScalar(prhs[8]),
		(GLsizei)mxGetScalar(prhs[9]),
		(const GLvoid*)mxGetData(prhs[10]));

}

void gl_compressedtexsubimage2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCompressedTexSubImage2D) mogl_glunsupported("glCompressedTexSubImage2D");
	glCompressedTexSubImage2D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLenum)mxGetScalar(prhs[6]),
		(GLsizei)mxGetScalar(prhs[7]),
		(const GLvoid*)mxGetData(prhs[8]));

}

void gl_compressedtexsubimage1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCompressedTexSubImage1D) mogl_glunsupported("glCompressedTexSubImage1D");
	glCompressedTexSubImage1D((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(const GLvoid*)mxGetData(prhs[6]));

}

void gl_getcompressedteximage( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetCompressedTexImage) mogl_glunsupported("glGetCompressedTexImage");
	glGetCompressedTexImage((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLvoid*)mxGetData(prhs[2]));

}

void gl_activetexture( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glActiveTexture) mogl_glunsupported("glActiveTexture");
	glActiveTexture((GLenum)mxGetScalar(prhs[0]));

}

void gl_clientactivetexture( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClientActiveTexture) mogl_glunsupported("glClientActiveTexture");
	glClientActiveTexture((GLenum)mxGetScalar(prhs[0]));

}

void gl_multitexcoord1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord1d) mogl_glunsupported("glMultiTexCoord1d");
	glMultiTexCoord1d((GLenum)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_multitexcoord1dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord1dv) mogl_glunsupported("glMultiTexCoord1dv");
	glMultiTexCoord1dv((GLenum)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_multitexcoord1f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord1f) mogl_glunsupported("glMultiTexCoord1f");
	glMultiTexCoord1f((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_multitexcoord1fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord1fv) mogl_glunsupported("glMultiTexCoord1fv");
	glMultiTexCoord1fv((GLenum)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_multitexcoord1i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord1i) mogl_glunsupported("glMultiTexCoord1i");
	glMultiTexCoord1i((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_multitexcoord1iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord1iv) mogl_glunsupported("glMultiTexCoord1iv");
	glMultiTexCoord1iv((GLenum)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_multitexcoord1s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord1s) mogl_glunsupported("glMultiTexCoord1s");
	glMultiTexCoord1s((GLenum)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]));

}

void gl_multitexcoord1sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord1sv) mogl_glunsupported("glMultiTexCoord1sv");
	glMultiTexCoord1sv((GLenum)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_multitexcoord2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord2d) mogl_glunsupported("glMultiTexCoord2d");
	glMultiTexCoord2d((GLenum)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_multitexcoord2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord2dv) mogl_glunsupported("glMultiTexCoord2dv");
	glMultiTexCoord2dv((GLenum)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_multitexcoord2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord2f) mogl_glunsupported("glMultiTexCoord2f");
	glMultiTexCoord2f((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_multitexcoord2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord2fv) mogl_glunsupported("glMultiTexCoord2fv");
	glMultiTexCoord2fv((GLenum)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_multitexcoord2i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord2i) mogl_glunsupported("glMultiTexCoord2i");
	glMultiTexCoord2i((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_multitexcoord2iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord2iv) mogl_glunsupported("glMultiTexCoord2iv");
	glMultiTexCoord2iv((GLenum)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_multitexcoord2s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord2s) mogl_glunsupported("glMultiTexCoord2s");
	glMultiTexCoord2s((GLenum)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_multitexcoord2sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord2sv) mogl_glunsupported("glMultiTexCoord2sv");
	glMultiTexCoord2sv((GLenum)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_multitexcoord3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord3d) mogl_glunsupported("glMultiTexCoord3d");
	glMultiTexCoord3d((GLenum)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_multitexcoord3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord3dv) mogl_glunsupported("glMultiTexCoord3dv");
	glMultiTexCoord3dv((GLenum)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_multitexcoord3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord3f) mogl_glunsupported("glMultiTexCoord3f");
	glMultiTexCoord3f((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_multitexcoord3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord3fv) mogl_glunsupported("glMultiTexCoord3fv");
	glMultiTexCoord3fv((GLenum)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_multitexcoord3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord3i) mogl_glunsupported("glMultiTexCoord3i");
	glMultiTexCoord3i((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_multitexcoord3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord3iv) mogl_glunsupported("glMultiTexCoord3iv");
	glMultiTexCoord3iv((GLenum)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_multitexcoord3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord3s) mogl_glunsupported("glMultiTexCoord3s");
	glMultiTexCoord3s((GLenum)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]));

}

void gl_multitexcoord3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord3sv) mogl_glunsupported("glMultiTexCoord3sv");
	glMultiTexCoord3sv((GLenum)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_multitexcoord4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord4d) mogl_glunsupported("glMultiTexCoord4d");
	glMultiTexCoord4d((GLenum)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]));

}

void gl_multitexcoord4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord4dv) mogl_glunsupported("glMultiTexCoord4dv");
	glMultiTexCoord4dv((GLenum)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_multitexcoord4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord4f) mogl_glunsupported("glMultiTexCoord4f");
	glMultiTexCoord4f((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]),
		(GLfloat)mxGetScalar(prhs[4]));

}

void gl_multitexcoord4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord4fv) mogl_glunsupported("glMultiTexCoord4fv");
	glMultiTexCoord4fv((GLenum)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_multitexcoord4i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord4i) mogl_glunsupported("glMultiTexCoord4i");
	glMultiTexCoord4i((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_multitexcoord4iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord4iv) mogl_glunsupported("glMultiTexCoord4iv");
	glMultiTexCoord4iv((GLenum)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_multitexcoord4s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord4s) mogl_glunsupported("glMultiTexCoord4s");
	glMultiTexCoord4s((GLenum)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]),
		(GLshort)mxGetScalar(prhs[4]));

}

void gl_multitexcoord4sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoord4sv) mogl_glunsupported("glMultiTexCoord4sv");
	glMultiTexCoord4sv((GLenum)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_fogcoordf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogCoordf) mogl_glunsupported("glFogCoordf");
	glFogCoordf((GLfloat)mxGetScalar(prhs[0]));

}

void gl_fogcoordfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogCoordfv) mogl_glunsupported("glFogCoordfv");
	glFogCoordfv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_fogcoordd( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogCoordd) mogl_glunsupported("glFogCoordd");
	glFogCoordd((GLdouble)mxGetScalar(prhs[0]));

}

void gl_fogcoorddv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogCoorddv) mogl_glunsupported("glFogCoorddv");
	glFogCoorddv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_fogcoordpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFogCoordPointer) mogl_glunsupported("glFogCoordPointer");
	glFogCoordPointer((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLvoid*)mxGetData(prhs[2]));

}

void gl_secondarycolor3b( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3b) mogl_glunsupported("glSecondaryColor3b");
	glSecondaryColor3b((GLbyte)mxGetScalar(prhs[0]),
		(GLbyte)mxGetScalar(prhs[1]),
		(GLbyte)mxGetScalar(prhs[2]));

}

void gl_secondarycolor3bv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3bv) mogl_glunsupported("glSecondaryColor3bv");
	glSecondaryColor3bv((const GLbyte*)mxGetData(prhs[0]));

}

void gl_secondarycolor3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3d) mogl_glunsupported("glSecondaryColor3d");
	glSecondaryColor3d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_secondarycolor3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3dv) mogl_glunsupported("glSecondaryColor3dv");
	glSecondaryColor3dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_secondarycolor3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3f) mogl_glunsupported("glSecondaryColor3f");
	glSecondaryColor3f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_secondarycolor3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3fv) mogl_glunsupported("glSecondaryColor3fv");
	glSecondaryColor3fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_secondarycolor3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3i) mogl_glunsupported("glSecondaryColor3i");
	glSecondaryColor3i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_secondarycolor3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3iv) mogl_glunsupported("glSecondaryColor3iv");
	glSecondaryColor3iv((const GLint*)mxGetData(prhs[0]));

}

void gl_secondarycolor3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3s) mogl_glunsupported("glSecondaryColor3s");
	glSecondaryColor3s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_secondarycolor3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3sv) mogl_glunsupported("glSecondaryColor3sv");
	glSecondaryColor3sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_secondarycolor3ub( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3ub) mogl_glunsupported("glSecondaryColor3ub");
	glSecondaryColor3ub((GLubyte)mxGetScalar(prhs[0]),
		(GLubyte)mxGetScalar(prhs[1]),
		(GLubyte)mxGetScalar(prhs[2]));

}

void gl_secondarycolor3ubv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3ubv) mogl_glunsupported("glSecondaryColor3ubv");
	glSecondaryColor3ubv((const GLubyte*)mxGetData(prhs[0]));

}

void gl_secondarycolor3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3ui) mogl_glunsupported("glSecondaryColor3ui");
	glSecondaryColor3ui((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_secondarycolor3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3uiv) mogl_glunsupported("glSecondaryColor3uiv");
	glSecondaryColor3uiv((const GLuint*)mxGetData(prhs[0]));

}

void gl_secondarycolor3us( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3us) mogl_glunsupported("glSecondaryColor3us");
	glSecondaryColor3us((GLushort)mxGetScalar(prhs[0]),
		(GLushort)mxGetScalar(prhs[1]),
		(GLushort)mxGetScalar(prhs[2]));

}

void gl_secondarycolor3usv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColor3usv) mogl_glunsupported("glSecondaryColor3usv");
	glSecondaryColor3usv((const GLushort*)mxGetData(prhs[0]));

}

void gl_secondarycolorpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColorPointer) mogl_glunsupported("glSecondaryColorPointer");
	glSecondaryColorPointer((GLint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLvoid*)mxGetData(prhs[3]));

}

void gl_pointparameterf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPointParameterf) mogl_glunsupported("glPointParameterf");
	glPointParameterf((GLenum)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_pointparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPointParameterfv) mogl_glunsupported("glPointParameterfv");
	glPointParameterfv((GLenum)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_pointparameteri( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPointParameteri) mogl_glunsupported("glPointParameteri");
	glPointParameteri((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_pointparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPointParameteriv) mogl_glunsupported("glPointParameteriv");
	glPointParameteriv((GLenum)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_blendfuncseparate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendFuncSeparate) mogl_glunsupported("glBlendFuncSeparate");
	glBlendFuncSeparate((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]));

}

void gl_multidrawarrays( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiDrawArrays) mogl_glunsupported("glMultiDrawArrays");
	glMultiDrawArrays((GLenum)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]),
		(const GLsizei*)mxGetData(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_windowpos2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos2d) mogl_glunsupported("glWindowPos2d");
	glWindowPos2d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_windowpos2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos2dv) mogl_glunsupported("glWindowPos2dv");
	glWindowPos2dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_windowpos2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos2f) mogl_glunsupported("glWindowPos2f");
	glWindowPos2f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_windowpos2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos2fv) mogl_glunsupported("glWindowPos2fv");
	glWindowPos2fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_windowpos2i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos2i) mogl_glunsupported("glWindowPos2i");
	glWindowPos2i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_windowpos2iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos2iv) mogl_glunsupported("glWindowPos2iv");
	glWindowPos2iv((const GLint*)mxGetData(prhs[0]));

}

void gl_windowpos2s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos2s) mogl_glunsupported("glWindowPos2s");
	glWindowPos2s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]));

}

void gl_windowpos2sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos2sv) mogl_glunsupported("glWindowPos2sv");
	glWindowPos2sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_windowpos3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos3d) mogl_glunsupported("glWindowPos3d");
	glWindowPos3d((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_windowpos3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos3dv) mogl_glunsupported("glWindowPos3dv");
	glWindowPos3dv((const GLdouble*)mxGetData(prhs[0]));

}

void gl_windowpos3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos3f) mogl_glunsupported("glWindowPos3f");
	glWindowPos3f((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_windowpos3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos3fv) mogl_glunsupported("glWindowPos3fv");
	glWindowPos3fv((const GLfloat*)mxGetData(prhs[0]));

}

void gl_windowpos3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos3i) mogl_glunsupported("glWindowPos3i");
	glWindowPos3i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_windowpos3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos3iv) mogl_glunsupported("glWindowPos3iv");
	glWindowPos3iv((const GLint*)mxGetData(prhs[0]));

}

void gl_windowpos3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos3s) mogl_glunsupported("glWindowPos3s");
	glWindowPos3s((GLshort)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_windowpos3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWindowPos3sv) mogl_glunsupported("glWindowPos3sv");
	glWindowPos3sv((const GLshort*)mxGetData(prhs[0]));

}

void gl_genqueries( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenQueries) mogl_glunsupported("glGenQueries");
	glGenQueries((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_deletequeries( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteQueries) mogl_glunsupported("glDeleteQueries");
	glDeleteQueries((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_isquery( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsQuery) mogl_glunsupported("glIsQuery");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsQuery((GLuint)mxGetScalar(prhs[0]));

}

void gl_beginquery( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBeginQuery) mogl_glunsupported("glBeginQuery");
	glBeginQuery((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_endquery( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEndQuery) mogl_glunsupported("glEndQuery");
	glEndQuery((GLenum)mxGetScalar(prhs[0]));

}

void gl_getqueryiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetQueryiv) mogl_glunsupported("glGetQueryiv");
	glGetQueryiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getqueryobjectiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetQueryObjectiv) mogl_glunsupported("glGetQueryObjectiv");
	glGetQueryObjectiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getqueryobjectuiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetQueryObjectuiv) mogl_glunsupported("glGetQueryObjectuiv");
	glGetQueryObjectuiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint*)mxGetData(prhs[2]));

}

void gl_bindbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindBuffer) mogl_glunsupported("glBindBuffer");
	glBindBuffer((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_deletebuffers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteBuffers) mogl_glunsupported("glDeleteBuffers");
	glDeleteBuffers((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_genbuffers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenBuffers) mogl_glunsupported("glGenBuffers");
	glGenBuffers((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_isbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsBuffer) mogl_glunsupported("glIsBuffer");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsBuffer((GLuint)mxGetScalar(prhs[0]));

}

void gl_buffersubdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBufferSubData) mogl_glunsupported("glBufferSubData");
	glBufferSubData((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLvoid*)mxGetData(prhs[3]));

}

void gl_getbuffersubdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetBufferSubData) mogl_glunsupported("glGetBufferSubData");
	glGetBufferSubData((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLvoid*)mxGetData(prhs[3]));

}

void gl_mapbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMapBuffer) mogl_glunsupported("glMapBuffer");
	plhs[0]=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	*(double *)mxGetData(plhs[0])=PsychPtrToDouble((void*) glMapBuffer((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1])));

}

void gl_unmapbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUnmapBuffer) mogl_glunsupported("glUnmapBuffer");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glUnmapBuffer((GLenum)mxGetScalar(prhs[0]));

}

void gl_getbufferparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetBufferParameteriv) mogl_glunsupported("glGetBufferParameteriv");
	glGetBufferParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_drawbuffers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawBuffers) mogl_glunsupported("glDrawBuffers");
	glDrawBuffers((GLsizei)mxGetScalar(prhs[0]),
		(const GLenum*)mxGetData(prhs[1]));

}

void gl_vertexattrib1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib1d) mogl_glunsupported("glVertexAttrib1d");
	glVertexAttrib1d((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_vertexattrib1dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib1dv) mogl_glunsupported("glVertexAttrib1dv");
	glVertexAttrib1dv((GLuint)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_vertexattrib1f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib1f) mogl_glunsupported("glVertexAttrib1f");
	glVertexAttrib1f((GLuint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_vertexattrib1fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib1fv) mogl_glunsupported("glVertexAttrib1fv");
	glVertexAttrib1fv((GLuint)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_vertexattrib1s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib1s) mogl_glunsupported("glVertexAttrib1s");
	glVertexAttrib1s((GLuint)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]));

}

void gl_vertexattrib1sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib1sv) mogl_glunsupported("glVertexAttrib1sv");
	glVertexAttrib1sv((GLuint)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_vertexattrib2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib2d) mogl_glunsupported("glVertexAttrib2d");
	glVertexAttrib2d((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_vertexattrib2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib2dv) mogl_glunsupported("glVertexAttrib2dv");
	glVertexAttrib2dv((GLuint)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_vertexattrib2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib2f) mogl_glunsupported("glVertexAttrib2f");
	glVertexAttrib2f((GLuint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_vertexattrib2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib2fv) mogl_glunsupported("glVertexAttrib2fv");
	glVertexAttrib2fv((GLuint)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_vertexattrib2s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib2s) mogl_glunsupported("glVertexAttrib2s");
	glVertexAttrib2s((GLuint)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]));

}

void gl_vertexattrib2sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib2sv) mogl_glunsupported("glVertexAttrib2sv");
	glVertexAttrib2sv((GLuint)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_vertexattrib3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib3d) mogl_glunsupported("glVertexAttrib3d");
	glVertexAttrib3d((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_vertexattrib3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib3dv) mogl_glunsupported("glVertexAttrib3dv");
	glVertexAttrib3dv((GLuint)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_vertexattrib3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib3f) mogl_glunsupported("glVertexAttrib3f");
	glVertexAttrib3f((GLuint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_vertexattrib3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib3fv) mogl_glunsupported("glVertexAttrib3fv");
	glVertexAttrib3fv((GLuint)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_vertexattrib3s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib3s) mogl_glunsupported("glVertexAttrib3s");
	glVertexAttrib3s((GLuint)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]));

}

void gl_vertexattrib3sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib3sv) mogl_glunsupported("glVertexAttrib3sv");
	glVertexAttrib3sv((GLuint)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_vertexattrib4nbv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4Nbv) mogl_glunsupported("glVertexAttrib4Nbv");
	glVertexAttrib4Nbv((GLuint)mxGetScalar(prhs[0]),
		(const GLbyte*)mxGetData(prhs[1]));

}

void gl_vertexattrib4niv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4Niv) mogl_glunsupported("glVertexAttrib4Niv");
	glVertexAttrib4Niv((GLuint)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_vertexattrib4nsv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4Nsv) mogl_glunsupported("glVertexAttrib4Nsv");
	glVertexAttrib4Nsv((GLuint)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_vertexattrib4nub( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4Nub) mogl_glunsupported("glVertexAttrib4Nub");
	glVertexAttrib4Nub((GLuint)mxGetScalar(prhs[0]),
		(GLubyte)mxGetScalar(prhs[1]),
		(GLubyte)mxGetScalar(prhs[2]),
		(GLubyte)mxGetScalar(prhs[3]),
		(GLubyte)mxGetScalar(prhs[4]));

}

void gl_vertexattrib4nubv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4Nubv) mogl_glunsupported("glVertexAttrib4Nubv");
	glVertexAttrib4Nubv((GLuint)mxGetScalar(prhs[0]),
		(const GLubyte*)mxGetData(prhs[1]));

}

void gl_vertexattrib4nuiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4Nuiv) mogl_glunsupported("glVertexAttrib4Nuiv");
	glVertexAttrib4Nuiv((GLuint)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_vertexattrib4nusv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4Nusv) mogl_glunsupported("glVertexAttrib4Nusv");
	glVertexAttrib4Nusv((GLuint)mxGetScalar(prhs[0]),
		(const GLushort*)mxGetData(prhs[1]));

}

void gl_vertexattrib4bv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4bv) mogl_glunsupported("glVertexAttrib4bv");
	glVertexAttrib4bv((GLuint)mxGetScalar(prhs[0]),
		(const GLbyte*)mxGetData(prhs[1]));

}

void gl_vertexattrib4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4d) mogl_glunsupported("glVertexAttrib4d");
	glVertexAttrib4d((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]));

}

void gl_vertexattrib4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4dv) mogl_glunsupported("glVertexAttrib4dv");
	glVertexAttrib4dv((GLuint)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_vertexattrib4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4f) mogl_glunsupported("glVertexAttrib4f");
	glVertexAttrib4f((GLuint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]),
		(GLfloat)mxGetScalar(prhs[4]));

}

void gl_vertexattrib4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4fv) mogl_glunsupported("glVertexAttrib4fv");
	glVertexAttrib4fv((GLuint)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_vertexattrib4iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4iv) mogl_glunsupported("glVertexAttrib4iv");
	glVertexAttrib4iv((GLuint)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_vertexattrib4s( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4s) mogl_glunsupported("glVertexAttrib4s");
	glVertexAttrib4s((GLuint)mxGetScalar(prhs[0]),
		(GLshort)mxGetScalar(prhs[1]),
		(GLshort)mxGetScalar(prhs[2]),
		(GLshort)mxGetScalar(prhs[3]),
		(GLshort)mxGetScalar(prhs[4]));

}

void gl_vertexattrib4sv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4sv) mogl_glunsupported("glVertexAttrib4sv");
	glVertexAttrib4sv((GLuint)mxGetScalar(prhs[0]),
		(const GLshort*)mxGetData(prhs[1]));

}

void gl_vertexattrib4ubv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4ubv) mogl_glunsupported("glVertexAttrib4ubv");
	glVertexAttrib4ubv((GLuint)mxGetScalar(prhs[0]),
		(const GLubyte*)mxGetData(prhs[1]));

}

void gl_vertexattrib4uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4uiv) mogl_glunsupported("glVertexAttrib4uiv");
	glVertexAttrib4uiv((GLuint)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_vertexattrib4usv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttrib4usv) mogl_glunsupported("glVertexAttrib4usv");
	glVertexAttrib4usv((GLuint)mxGetScalar(prhs[0]),
		(const GLushort*)mxGetData(prhs[1]));

}

void gl_enablevertexattribarray( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEnableVertexAttribArray) mogl_glunsupported("glEnableVertexAttribArray");
	glEnableVertexAttribArray((GLuint)mxGetScalar(prhs[0]));

}

void gl_disablevertexattribarray( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDisableVertexAttribArray) mogl_glunsupported("glDisableVertexAttribArray");
	glDisableVertexAttribArray((GLuint)mxGetScalar(prhs[0]));

}

void gl_getvertexattribdv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetVertexAttribdv) mogl_glunsupported("glGetVertexAttribdv");
	glGetVertexAttribdv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLdouble*)mxGetData(prhs[2]));

}

void gl_getvertexattribfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetVertexAttribfv) mogl_glunsupported("glGetVertexAttribfv");
	glGetVertexAttribfv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getvertexattribiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetVertexAttribiv) mogl_glunsupported("glGetVertexAttribiv");
	glGetVertexAttribiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_deleteshader( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteShader) mogl_glunsupported("glDeleteShader");
	glDeleteShader((GLuint)mxGetScalar(prhs[0]));

}

void gl_detachshader( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDetachShader) mogl_glunsupported("glDetachShader");
	glDetachShader((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_createshader( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCreateShader) mogl_glunsupported("glCreateShader");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glCreateShader((GLenum)mxGetScalar(prhs[0]));

}

void gl_compileshader( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCompileShader) mogl_glunsupported("glCompileShader");
	glCompileShader((GLuint)mxGetScalar(prhs[0]));

}

void gl_createprogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCreateProgram) mogl_glunsupported("glCreateProgram");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glCreateProgram();

}

void gl_attachshader( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glAttachShader) mogl_glunsupported("glAttachShader");
	glAttachShader((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_linkprogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glLinkProgram) mogl_glunsupported("glLinkProgram");
	glLinkProgram((GLuint)mxGetScalar(prhs[0]));

}

void gl_useprogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUseProgram) mogl_glunsupported("glUseProgram");
	glUseProgram((GLuint)mxGetScalar(prhs[0]));

}

void gl_deleteprogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteProgram) mogl_glunsupported("glDeleteProgram");
	glDeleteProgram((GLuint)mxGetScalar(prhs[0]));

}

void gl_validateprogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glValidateProgram) mogl_glunsupported("glValidateProgram");
	glValidateProgram((GLuint)mxGetScalar(prhs[0]));

}

void gl_uniform1f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform1f) mogl_glunsupported("glUniform1f");
	glUniform1f((GLint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_uniform2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform2f) mogl_glunsupported("glUniform2f");
	glUniform2f((GLint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_uniform3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform3f) mogl_glunsupported("glUniform3f");
	glUniform3f((GLint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_uniform4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform4f) mogl_glunsupported("glUniform4f");
	glUniform4f((GLint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]),
		(GLfloat)mxGetScalar(prhs[4]));

}

void gl_uniform1i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform1i) mogl_glunsupported("glUniform1i");
	glUniform1i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_uniform2i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform2i) mogl_glunsupported("glUniform2i");
	glUniform2i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_uniform3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform3i) mogl_glunsupported("glUniform3i");
	glUniform3i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_uniform4i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform4i) mogl_glunsupported("glUniform4i");
	glUniform4i((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_uniform1fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform1fv) mogl_glunsupported("glUniform1fv");
	glUniform1fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_uniform2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform2fv) mogl_glunsupported("glUniform2fv");
	glUniform2fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_uniform3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform3fv) mogl_glunsupported("glUniform3fv");
	glUniform3fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_uniform4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform4fv) mogl_glunsupported("glUniform4fv");
	glUniform4fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_uniform1iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform1iv) mogl_glunsupported("glUniform1iv");
	glUniform1iv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_uniform2iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform2iv) mogl_glunsupported("glUniform2iv");
	glUniform2iv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_uniform3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform3iv) mogl_glunsupported("glUniform3iv");
	glUniform3iv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_uniform4iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform4iv) mogl_glunsupported("glUniform4iv");
	glUniform4iv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_uniformmatrix2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix2fv) mogl_glunsupported("glUniformMatrix2fv");
	glUniformMatrix2fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_uniformmatrix3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix3fv) mogl_glunsupported("glUniformMatrix3fv");
	glUniformMatrix3fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_uniformmatrix4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix4fv) mogl_glunsupported("glUniformMatrix4fv");
	glUniformMatrix4fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_isshader( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsShader) mogl_glunsupported("glIsShader");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsShader((GLuint)mxGetScalar(prhs[0]));

}

void gl_isprogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsProgram) mogl_glunsupported("glIsProgram");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsProgram((GLuint)mxGetScalar(prhs[0]));

}

void gl_getshaderiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetShaderiv) mogl_glunsupported("glGetShaderiv");
	glGetShaderiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getprogramiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramiv) mogl_glunsupported("glGetProgramiv");
	glGetProgramiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getattachedshaders( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetAttachedShaders) mogl_glunsupported("glGetAttachedShaders");
	glGetAttachedShaders((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei*)mxGetData(prhs[2]),
		(GLuint*)mxGetData(prhs[3]));

}

void gl_getshaderinfolog( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetShaderInfoLog) mogl_glunsupported("glGetShaderInfoLog");
	glGetShaderInfoLog((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei*)mxGetData(prhs[2]),
		(GLchar*)mxGetData(prhs[3]));

}

void gl_getprograminfolog( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramInfoLog) mogl_glunsupported("glGetProgramInfoLog");
	glGetProgramInfoLog((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei*)mxGetData(prhs[2]),
		(GLchar*)mxGetData(prhs[3]));

}

void gl_getuniformlocation( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetUniformLocation) mogl_glunsupported("glGetUniformLocation");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetUniformLocation((GLuint)mxGetScalar(prhs[0]),
		(const GLchar*)mxGetData(prhs[1]));

}

void gl_getactiveuniform( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveUniform) mogl_glunsupported("glGetActiveUniform");
	glGetActiveUniform((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei*)mxGetData(prhs[3]),
		(GLint*)mxGetData(prhs[4]),
		(GLenum*)mxGetData(prhs[5]),
		(GLchar*)mxGetData(prhs[6]));

}

void gl_getuniformfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetUniformfv) mogl_glunsupported("glGetUniformfv");
	glGetUniformfv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getuniformiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetUniformiv) mogl_glunsupported("glGetUniformiv");
	glGetUniformiv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getshadersource( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetShaderSource) mogl_glunsupported("glGetShaderSource");
	glGetShaderSource((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei*)mxGetData(prhs[2]),
		(GLchar*)mxGetData(prhs[3]));

}

void gl_bindattriblocation( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindAttribLocation) mogl_glunsupported("glBindAttribLocation");
	glBindAttribLocation((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]));

}

void gl_getactiveattrib( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveAttrib) mogl_glunsupported("glGetActiveAttrib");
	glGetActiveAttrib((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei*)mxGetData(prhs[3]),
		(GLint*)mxGetData(prhs[4]),
		(GLenum*)mxGetData(prhs[5]),
		(GLchar*)mxGetData(prhs[6]));

}

void gl_getattriblocation( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetAttribLocation) mogl_glunsupported("glGetAttribLocation");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetAttribLocation((GLuint)mxGetScalar(prhs[0]),
		(const GLchar*)mxGetData(prhs[1]));

}

void gl_stencilfuncseparate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glStencilFuncSeparate) mogl_glunsupported("glStencilFuncSeparate");
	glStencilFuncSeparate((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_stencilopseparate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glStencilOpSeparate) mogl_glunsupported("glStencilOpSeparate");
	glStencilOpSeparate((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]));

}

void gl_stencilmaskseparate( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glStencilMaskSeparate) mogl_glunsupported("glStencilMaskSeparate");
	glStencilMaskSeparate((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_uniformmatrix2x3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix2x3fv) mogl_glunsupported("glUniformMatrix2x3fv");
	glUniformMatrix2x3fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_uniformmatrix3x2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix3x2fv) mogl_glunsupported("glUniformMatrix3x2fv");
	glUniformMatrix3x2fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_uniformmatrix2x4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix2x4fv) mogl_glunsupported("glUniformMatrix2x4fv");
	glUniformMatrix2x4fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_uniformmatrix4x2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix4x2fv) mogl_glunsupported("glUniformMatrix4x2fv");
	glUniformMatrix4x2fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_uniformmatrix3x4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix3x4fv) mogl_glunsupported("glUniformMatrix3x4fv");
	glUniformMatrix3x4fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_uniformmatrix4x3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix4x3fv) mogl_glunsupported("glUniformMatrix4x3fv");
	glUniformMatrix4x3fv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void glu_begincurve( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBeginCurve) mogl_glunsupported("gluBeginCurve");
	gluBeginCurve((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void glu_beginsurface( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBeginSurface) mogl_glunsupported("gluBeginSurface");
	gluBeginSurface((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void glu_begintrim( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBeginTrim) mogl_glunsupported("gluBeginTrim");
	gluBeginTrim((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void glu_build1dmipmaplevels( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBuild1DMipmapLevels) mogl_glunsupported("gluBuild1DMipmapLevels");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluBuild1DMipmapLevels((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLint)mxGetScalar(prhs[6]),
		(GLint)mxGetScalar(prhs[7]),
		(const void*)mxGetData(prhs[8]));

}

void glu_build1dmipmaps( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBuild1DMipmaps) mogl_glunsupported("gluBuild1DMipmaps");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluBuild1DMipmaps((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(const void*)mxGetData(prhs[5]));

}

void glu_build2dmipmaplevels( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBuild2DMipmapLevels) mogl_glunsupported("gluBuild2DMipmapLevels");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluBuild2DMipmapLevels((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(GLint)mxGetScalar(prhs[6]),
		(GLint)mxGetScalar(prhs[7]),
		(GLint)mxGetScalar(prhs[8]),
		(const void*)mxGetData(prhs[9]));

}

void glu_build2dmipmaps( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBuild2DMipmaps) mogl_glunsupported("gluBuild2DMipmaps");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluBuild2DMipmaps((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(const void*)mxGetData(prhs[6]));

}

void glu_build3dmipmaplevels( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBuild3DMipmapLevels) mogl_glunsupported("gluBuild3DMipmapLevels");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluBuild3DMipmapLevels((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(GLenum)mxGetScalar(prhs[6]),
		(GLint)mxGetScalar(prhs[7]),
		(GLint)mxGetScalar(prhs[8]),
		(GLint)mxGetScalar(prhs[9]),
		(const void*)mxGetData(prhs[10]));

}

void glu_build3dmipmaps( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluBuild3DMipmaps) mogl_glunsupported("gluBuild3DMipmaps");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluBuild3DMipmaps((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(GLenum)mxGetScalar(prhs[6]),
		(const void*)mxGetData(prhs[7]));

}

void glu_checkextension( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluCheckExtension) mogl_glunsupported("gluCheckExtension");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluCheckExtension((const GLubyte*)mxGetData(prhs[0]),
		(const GLubyte*)mxGetData(prhs[1]));

}

void glu_cylinder( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluCylinder) mogl_glunsupported("gluCylinder");
	gluCylinder((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]));

}

void glu_deletenurbsrenderer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluDeleteNurbsRenderer) mogl_glunsupported("gluDeleteNurbsRenderer");
	gluDeleteNurbsRenderer((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void glu_deletequadric( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluDeleteQuadric) mogl_glunsupported("gluDeleteQuadric");
	gluDeleteQuadric((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void glu_disk( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluDisk) mogl_glunsupported("gluDisk");
	gluDisk((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void glu_endcurve( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluEndCurve) mogl_glunsupported("gluEndCurve");
	gluEndCurve((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void glu_endsurface( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluEndSurface) mogl_glunsupported("gluEndSurface");
	gluEndSurface((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void glu_endtrim( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluEndTrim) mogl_glunsupported("gluEndTrim");
	gluEndTrim((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void glu_getnurbsproperty( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluGetNurbsProperty) mogl_glunsupported("gluGetNurbsProperty");
	gluGetNurbsProperty((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void glu_gettessproperty( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluGetTessProperty) mogl_glunsupported("gluGetTessProperty");
	gluGetTessProperty((GLUtesselator*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLenum)mxGetScalar(prhs[1]),
		(GLdouble*)mxGetData(prhs[2]));

}

void glu_loadsamplingmatrices( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluLoadSamplingMatrices) mogl_glunsupported("gluLoadSamplingMatrices");
	gluLoadSamplingMatrices((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(const GLfloat*)mxGetData(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]),
		(const GLint*)mxGetData(prhs[3]));

}

void glu_lookat( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluLookAt) mogl_glunsupported("gluLookAt");
	gluLookAt((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]),
		(GLdouble)mxGetScalar(prhs[5]),
		(GLdouble)mxGetScalar(prhs[6]),
		(GLdouble)mxGetScalar(prhs[7]),
		(GLdouble)mxGetScalar(prhs[8]));

}

void glu_newnurbsrenderer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluNewNurbsRenderer) mogl_glunsupported("gluNewNurbsRenderer");
	plhs[0]=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	*(double *)mxGetData(plhs[0])=PsychPtrToDouble((void*) gluNewNurbsRenderer());

}

void glu_newquadric( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluNewQuadric) mogl_glunsupported("gluNewQuadric");
	plhs[0]=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	*(double *)mxGetData(plhs[0])=PsychPtrToDouble((void*) gluNewQuadric());

}

void glu_nurbscurve( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluNurbsCurve) mogl_glunsupported("gluNurbsCurve");
	gluNurbsCurve((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLint)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLfloat*)mxGetData(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLenum)mxGetScalar(prhs[6]));

}

void glu_nurbsproperty( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluNurbsProperty) mogl_glunsupported("gluNurbsProperty");
	gluNurbsProperty((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void glu_nurbssurface( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluNurbsSurface) mogl_glunsupported("gluNurbsSurface");
	gluNurbsSurface((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLint)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLfloat*)mxGetData(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLint)mxGetScalar(prhs[6]),
		(GLfloat*)mxGetData(prhs[7]),
		(GLint)mxGetScalar(prhs[8]),
		(GLint)mxGetScalar(prhs[9]),
		(GLenum)mxGetScalar(prhs[10]));

}

void glu_ortho2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluOrtho2D) mogl_glunsupported("gluOrtho2D");
	gluOrtho2D((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void glu_partialdisk( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluPartialDisk) mogl_glunsupported("gluPartialDisk");
	gluPartialDisk((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLdouble)mxGetScalar(prhs[5]),
		(GLdouble)mxGetScalar(prhs[6]));

}

void glu_perspective( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluPerspective) mogl_glunsupported("gluPerspective");
	gluPerspective((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void glu_pickmatrix( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluPickMatrix) mogl_glunsupported("gluPickMatrix");
	gluPickMatrix((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLint*)mxGetData(prhs[4]));

}

void glu_project( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluProject) mogl_glunsupported("gluProject");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluProject((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]),
		(const GLint*)mxGetData(prhs[5]),
		(GLdouble*)mxGetData(prhs[6]),
		(GLdouble*)mxGetData(prhs[7]),
		(GLdouble*)mxGetData(prhs[8]));

}

void glu_pwlcurve( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluPwlCurve) mogl_glunsupported("gluPwlCurve");
	gluPwlCurve((GLUnurbs*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLint)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]));

}

void glu_quadricdrawstyle( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluQuadricDrawStyle) mogl_glunsupported("gluQuadricDrawStyle");
	gluQuadricDrawStyle((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLenum)mxGetScalar(prhs[1]));

}

void glu_quadricnormals( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluQuadricNormals) mogl_glunsupported("gluQuadricNormals");
	gluQuadricNormals((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLenum)mxGetScalar(prhs[1]));

}

void glu_quadricorientation( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluQuadricOrientation) mogl_glunsupported("gluQuadricOrientation");
	gluQuadricOrientation((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLenum)mxGetScalar(prhs[1]));

}

void glu_quadrictexture( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluQuadricTexture) mogl_glunsupported("gluQuadricTexture");
	gluQuadricTexture((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLboolean)mxGetScalar(prhs[1]));

}

void glu_scaleimage( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluScaleImage) mogl_glunsupported("gluScaleImage");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluScaleImage((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(const void*)mxGetData(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(GLenum)mxGetScalar(prhs[7]),
		(GLvoid*)mxGetData(prhs[8]));

}

void glu_sphere( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluSphere) mogl_glunsupported("gluSphere");
	gluSphere((GLUquadric*) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void glu_unproject( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluUnProject) mogl_glunsupported("gluUnProject");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluUnProject((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]),
		(const GLint*)mxGetData(prhs[5]),
		(GLdouble*)mxGetData(prhs[6]),
		(GLdouble*)mxGetData(prhs[7]),
		(GLdouble*)mxGetData(prhs[8]));

}

void glu_unproject4( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == gluUnProject4) mogl_glunsupported("gluUnProject4");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)gluUnProject4((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]),
		(const GLdouble*)mxGetData(prhs[5]),
		(const GLint*)mxGetData(prhs[6]),
		(GLdouble)mxGetScalar(prhs[7]),
		(GLdouble)mxGetScalar(prhs[8]),
		(GLdouble*)mxGetData(prhs[9]),
		(GLdouble*)mxGetData(prhs[10]),
		(GLdouble*)mxGetData(prhs[11]),
		(GLdouble*)mxGetData(prhs[12]));

}

void glut_wiresphere( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireSphere) mogl_glunsupported("glutWireSphere");
	glutWireSphere((GLdouble)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void glut_solidsphere( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidSphere) mogl_glunsupported("glutSolidSphere");
	glutSolidSphere((GLdouble)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void glut_wirecone( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireCone) mogl_glunsupported("glutWireCone");
	glutWireCone((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void glut_solidcone( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidCone) mogl_glunsupported("glutSolidCone");
	glutSolidCone((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void glut_wirecube( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireCube) mogl_glunsupported("glutWireCube");
	glutWireCube((GLdouble)mxGetScalar(prhs[0]));

}

void glut_solidcube( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidCube) mogl_glunsupported("glutSolidCube");
	glutSolidCube((GLdouble)mxGetScalar(prhs[0]));

}

void glut_wiretorus( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireTorus) mogl_glunsupported("glutWireTorus");
	glutWireTorus((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void glut_solidtorus( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidTorus) mogl_glunsupported("glutSolidTorus");
	glutSolidTorus((GLdouble)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void glut_wiredodecahedron( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireDodecahedron) mogl_glunsupported("glutWireDodecahedron");
	glutWireDodecahedron();

}

void glut_soliddodecahedron( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidDodecahedron) mogl_glunsupported("glutSolidDodecahedron");
	glutSolidDodecahedron();

}

void glut_wireteapot( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireTeapot) mogl_glunsupported("glutWireTeapot");
	glutWireTeapot((GLdouble)mxGetScalar(prhs[0]));

}

void glut_solidteapot( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidTeapot) mogl_glunsupported("glutSolidTeapot");
	glutSolidTeapot((GLdouble)mxGetScalar(prhs[0]));

}

void glut_wireoctahedron( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireOctahedron) mogl_glunsupported("glutWireOctahedron");
	glutWireOctahedron();

}

void glut_solidoctahedron( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidOctahedron) mogl_glunsupported("glutSolidOctahedron");
	glutSolidOctahedron();

}

void glut_wiretetrahedron( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireTetrahedron) mogl_glunsupported("glutWireTetrahedron");
	glutWireTetrahedron();

}

void glut_solidtetrahedron( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidTetrahedron) mogl_glunsupported("glutSolidTetrahedron");
	glutSolidTetrahedron();

}

void glut_wireicosahedron( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutWireIcosahedron) mogl_glunsupported("glutWireIcosahedron");
	glutWireIcosahedron();

}

void glut_solidicosahedron( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glutSolidIcosahedron) mogl_glunsupported("glutSolidIcosahedron");
	glutSolidIcosahedron();

}

void gl_isrenderbufferext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsRenderbufferEXT) mogl_glunsupported("glIsRenderbufferEXT");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsRenderbufferEXT((GLuint)mxGetScalar(prhs[0]));

}

void gl_bindrenderbufferext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindRenderbufferEXT) mogl_glunsupported("glBindRenderbufferEXT");
	glBindRenderbufferEXT((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_deleterenderbuffersext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteRenderbuffersEXT) mogl_glunsupported("glDeleteRenderbuffersEXT");
	glDeleteRenderbuffersEXT((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_genrenderbuffersext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenRenderbuffersEXT) mogl_glunsupported("glGenRenderbuffersEXT");
	glGenRenderbuffersEXT((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_renderbufferstorageext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRenderbufferStorageEXT) mogl_glunsupported("glRenderbufferStorageEXT");
	glRenderbufferStorageEXT((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_getrenderbufferparameterivext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetRenderbufferParameterivEXT) mogl_glunsupported("glGetRenderbufferParameterivEXT");
	glGetRenderbufferParameterivEXT((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_isframebufferext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsFramebufferEXT) mogl_glunsupported("glIsFramebufferEXT");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsFramebufferEXT((GLuint)mxGetScalar(prhs[0]));

}

void gl_bindframebufferext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindFramebufferEXT) mogl_glunsupported("glBindFramebufferEXT");
	glBindFramebufferEXT((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_deleteframebuffersext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteFramebuffersEXT) mogl_glunsupported("glDeleteFramebuffersEXT");
	glDeleteFramebuffersEXT((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_genframebuffersext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenFramebuffersEXT) mogl_glunsupported("glGenFramebuffersEXT");
	glGenFramebuffersEXT((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_checkframebufferstatusext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCheckFramebufferStatusEXT) mogl_glunsupported("glCheckFramebufferStatusEXT");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glCheckFramebufferStatusEXT((GLenum)mxGetScalar(prhs[0]));

}

void gl_framebuffertexture1dext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTexture1DEXT) mogl_glunsupported("glFramebufferTexture1DEXT");
	glFramebufferTexture1DEXT((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_framebuffertexture2dext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTexture2DEXT) mogl_glunsupported("glFramebufferTexture2DEXT");
	glFramebufferTexture2DEXT((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_framebuffertexture3dext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTexture3DEXT) mogl_glunsupported("glFramebufferTexture3DEXT");
	glFramebufferTexture3DEXT((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]));

}

void gl_framebufferrenderbufferext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferRenderbufferEXT) mogl_glunsupported("glFramebufferRenderbufferEXT");
	glFramebufferRenderbufferEXT((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_getframebufferattachmentparameterivext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetFramebufferAttachmentParameterivEXT) mogl_glunsupported("glGetFramebufferAttachmentParameterivEXT");
	glGetFramebufferAttachmentParameterivEXT((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_generatemipmapext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenerateMipmapEXT) mogl_glunsupported("glGenerateMipmapEXT");
	glGenerateMipmapEXT((GLenum)mxGetScalar(prhs[0]));

}

void gl_samplemasksgis( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSampleMaskSGIS) mogl_glunsupported("glSampleMaskSGIS");
	glSampleMaskSGIS((GLclampf)mxGetScalar(prhs[0]),
		(GLboolean)mxGetScalar(prhs[1]));

}

void gl_samplepatternsgis( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSamplePatternSGIS) mogl_glunsupported("glSamplePatternSGIS");
	glSamplePatternSGIS((GLenum)mxGetScalar(prhs[0]));

}

void gl_drawarraysinstancedarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawArraysInstancedARB) mogl_glunsupported("glDrawArraysInstancedARB");
	glDrawArraysInstancedARB((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_drawelementsinstancedarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawElementsInstancedARB) mogl_glunsupported("glDrawElementsInstancedARB");
	glDrawElementsInstancedARB((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(const GLvoid*)mxGetData(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_isrenderbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsRenderbuffer) mogl_glunsupported("glIsRenderbuffer");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsRenderbuffer((GLuint)mxGetScalar(prhs[0]));

}

void gl_bindrenderbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindRenderbuffer) mogl_glunsupported("glBindRenderbuffer");
	glBindRenderbuffer((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_deleterenderbuffers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteRenderbuffers) mogl_glunsupported("glDeleteRenderbuffers");
	glDeleteRenderbuffers((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_genrenderbuffers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenRenderbuffers) mogl_glunsupported("glGenRenderbuffers");
	glGenRenderbuffers((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_renderbufferstorage( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRenderbufferStorage) mogl_glunsupported("glRenderbufferStorage");
	glRenderbufferStorage((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_getrenderbufferparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetRenderbufferParameteriv) mogl_glunsupported("glGetRenderbufferParameteriv");
	glGetRenderbufferParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_isframebuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsFramebuffer) mogl_glunsupported("glIsFramebuffer");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsFramebuffer((GLuint)mxGetScalar(prhs[0]));

}

void gl_bindframebuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindFramebuffer) mogl_glunsupported("glBindFramebuffer");
	glBindFramebuffer((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_deleteframebuffers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteFramebuffers) mogl_glunsupported("glDeleteFramebuffers");
	glDeleteFramebuffers((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_genframebuffers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenFramebuffers) mogl_glunsupported("glGenFramebuffers");
	glGenFramebuffers((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_checkframebufferstatus( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCheckFramebufferStatus) mogl_glunsupported("glCheckFramebufferStatus");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glCheckFramebufferStatus((GLenum)mxGetScalar(prhs[0]));

}

void gl_framebuffertexture1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTexture1D) mogl_glunsupported("glFramebufferTexture1D");
	glFramebufferTexture1D((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_framebuffertexture2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTexture2D) mogl_glunsupported("glFramebufferTexture2D");
	glFramebufferTexture2D((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_framebuffertexture3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTexture3D) mogl_glunsupported("glFramebufferTexture3D");
	glFramebufferTexture3D((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]));

}

void gl_framebufferrenderbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferRenderbuffer) mogl_glunsupported("glFramebufferRenderbuffer");
	glFramebufferRenderbuffer((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_getframebufferattachmentparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetFramebufferAttachmentParameteriv) mogl_glunsupported("glGetFramebufferAttachmentParameteriv");
	glGetFramebufferAttachmentParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_generatemipmap( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenerateMipmap) mogl_glunsupported("glGenerateMipmap");
	glGenerateMipmap((GLenum)mxGetScalar(prhs[0]));

}

void gl_blitframebuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlitFramebuffer) mogl_glunsupported("glBlitFramebuffer");
	glBlitFramebuffer((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLint)mxGetScalar(prhs[6]),
		(GLint)mxGetScalar(prhs[7]),
		(GLbitfield)mxGetScalar(prhs[8]),
		(GLenum)mxGetScalar(prhs[9]));

}

void gl_renderbufferstoragemultisample( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glRenderbufferStorageMultisample) mogl_glunsupported("glRenderbufferStorageMultisample");
	glRenderbufferStorageMultisample((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_framebuffertexturelayer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTextureLayer) mogl_glunsupported("glFramebufferTextureLayer");
	glFramebufferTextureLayer((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_programparameteriarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramParameteriARB) mogl_glunsupported("glProgramParameteriARB");
	glProgramParameteriARB((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_framebuffertexturearb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTextureARB) mogl_glunsupported("glFramebufferTextureARB");
	glFramebufferTextureARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_framebuffertexturelayerarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTextureLayerARB) mogl_glunsupported("glFramebufferTextureLayerARB");
	glFramebufferTextureLayerARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_framebuffertexturefacearb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferTextureFaceARB) mogl_glunsupported("glFramebufferTextureFaceARB");
	glFramebufferTextureFaceARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]));

}

void gl_vertexattribdivisorarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribDivisorARB) mogl_glunsupported("glVertexAttribDivisorARB");
	glVertexAttribDivisorARB((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_mapbufferrange( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMapBufferRange) mogl_glunsupported("glMapBufferRange");
	plhs[0]=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	*(double *)mxGetData(plhs[0])=PsychPtrToDouble((void*) glMapBufferRange((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLbitfield)mxGetScalar(prhs[3])));

}

void gl_flushmappedbufferrange( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFlushMappedBufferRange) mogl_glunsupported("glFlushMappedBufferRange");
	glFlushMappedBufferRange((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]));

}

void gl_texbufferarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexBufferARB) mogl_glunsupported("glTexBufferARB");
	glTexBufferARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_bindvertexarray( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindVertexArray) mogl_glunsupported("glBindVertexArray");
	glBindVertexArray((GLuint)mxGetScalar(prhs[0]));

}

void gl_deletevertexarrays( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteVertexArrays) mogl_glunsupported("glDeleteVertexArrays");
	glDeleteVertexArrays((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_genvertexarrays( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenVertexArrays) mogl_glunsupported("glGenVertexArrays");
	glGenVertexArrays((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_isvertexarray( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsVertexArray) mogl_glunsupported("glIsVertexArray");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsVertexArray((GLuint)mxGetScalar(prhs[0]));

}

void gl_getuniformindices( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetUniformIndices) mogl_glunsupported("glGetUniformIndices");
	glGetUniformIndices((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]),
		(GLuint*)mxGetData(prhs[3]));

}

void gl_getactiveuniformsiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveUniformsiv) mogl_glunsupported("glGetActiveUniformsiv");
	glGetActiveUniformsiv((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLuint*)mxGetData(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLint*)mxGetData(prhs[4]));

}

void gl_getactiveuniformname( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveUniformName) mogl_glunsupported("glGetActiveUniformName");
	glGetActiveUniformName((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei*)mxGetData(prhs[3]),
		(GLchar*)mxGetData(prhs[4]));

}

void gl_getuniformblockindex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetUniformBlockIndex) mogl_glunsupported("glGetUniformBlockIndex");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetUniformBlockIndex((GLuint)mxGetScalar(prhs[0]),
		(const GLchar*)mxGetData(prhs[1]));

}

void gl_getactiveuniformblockiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveUniformBlockiv) mogl_glunsupported("glGetActiveUniformBlockiv");
	glGetActiveUniformBlockiv((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_getactiveuniformblockname( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveUniformBlockName) mogl_glunsupported("glGetActiveUniformBlockName");
	glGetActiveUniformBlockName((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei*)mxGetData(prhs[3]),
		(GLchar*)mxGetData(prhs[4]));

}

void gl_uniformblockbinding( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformBlockBinding) mogl_glunsupported("glUniformBlockBinding");
	glUniformBlockBinding((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_copybuffersubdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyBufferSubData) mogl_glunsupported("glCopyBufferSubData");
	glCopyBufferSubData((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_drawelementsbasevertex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawElementsBaseVertex) mogl_glunsupported("glDrawElementsBaseVertex");
	glDrawElementsBaseVertex((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(const GLvoid*)mxGetData(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_drawrangeelementsbasevertex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawRangeElementsBaseVertex) mogl_glunsupported("glDrawRangeElementsBaseVertex");
	glDrawRangeElementsBaseVertex((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(const GLvoid*)mxGetData(prhs[5]),
		(GLint)mxGetScalar(prhs[6]));

}

void gl_drawelementsinstancedbasevertex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawElementsInstancedBaseVertex) mogl_glunsupported("glDrawElementsInstancedBaseVertex");
	glDrawElementsInstancedBaseVertex((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(const GLvoid*)mxGetData(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]));

}

void gl_multidrawelementsbasevertex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiDrawElementsBaseVertex) mogl_glunsupported("glMultiDrawElementsBaseVertex");
	glMultiDrawElementsBaseVertex((GLenum)mxGetScalar(prhs[0]),
		(const GLsizei*)mxGetData(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(const GLvoid*)mxGetData(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(const GLint*)mxGetData(prhs[5]));

}

void gl_provokingvertex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProvokingVertex) mogl_glunsupported("glProvokingVertex");
	glProvokingVertex((GLenum)mxGetScalar(prhs[0]));

}

void gl_fencesync( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFenceSync) mogl_glunsupported("glFenceSync");
	plhs[0]=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	*(double *)mxGetData(plhs[0])=PsychPtrToDouble((void*) glFenceSync((GLenum)mxGetScalar(prhs[0]),
		(GLbitfield)mxGetScalar(prhs[1])));

}

void gl_issync( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsSync) mogl_glunsupported("glIsSync");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsSync((GLsync) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void gl_deletesync( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteSync) mogl_glunsupported("glDeleteSync");
	glDeleteSync((GLsync) PsychDoubleToPtr(mxGetScalar(prhs[0])));

}

void gl_clientwaitsync( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClientWaitSync) mogl_glunsupported("glClientWaitSync");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glClientWaitSync((GLsync) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLbitfield)mxGetScalar(prhs[1]),
		(GLuint64) *((GLuint64*) mxGetData(prhs[2])));

}

void gl_waitsync( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glWaitSync) mogl_glunsupported("glWaitSync");
	glWaitSync((GLsync) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLbitfield)mxGetScalar(prhs[1]),
		(GLuint64) *((GLuint64*) mxGetData(prhs[2])));

}

void gl_getinteger64v( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetInteger64v) mogl_glunsupported("glGetInteger64v");
	glGetInteger64v((GLenum)mxGetScalar(prhs[0]),
		(GLint64*)mxGetData(prhs[1]));

}

void gl_getsynciv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetSynciv) mogl_glunsupported("glGetSynciv");
	glGetSynciv((GLsync) PsychDoubleToPtr(mxGetScalar(prhs[0])),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei*)mxGetData(prhs[3]),
		(GLint*)mxGetData(prhs[4]));

}

void gl_teximage2dmultisample( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexImage2DMultisample) mogl_glunsupported("glTexImage2DMultisample");
	glTexImage2DMultisample((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLboolean)mxGetScalar(prhs[5]));

}

void gl_teximage3dmultisample( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexImage3DMultisample) mogl_glunsupported("glTexImage3DMultisample");
	glTexImage3DMultisample((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLboolean)mxGetScalar(prhs[6]));

}

void gl_getmultisamplefv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetMultisamplefv) mogl_glunsupported("glGetMultisamplefv");
	glGetMultisamplefv((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_samplemaski( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSampleMaski) mogl_glunsupported("glSampleMaski");
	glSampleMaski((GLuint)mxGetScalar(prhs[0]),
		(GLbitfield)mxGetScalar(prhs[1]));

}

void gl_blendequationiarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendEquationiARB) mogl_glunsupported("glBlendEquationiARB");
	glBlendEquationiARB((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]));

}

void gl_blendequationseparateiarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendEquationSeparateiARB) mogl_glunsupported("glBlendEquationSeparateiARB");
	glBlendEquationSeparateiARB((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]));

}

void gl_blendfunciarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendFunciARB) mogl_glunsupported("glBlendFunciARB");
	glBlendFunciARB((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]));

}

void gl_blendfuncseparateiarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBlendFuncSeparateiARB) mogl_glunsupported("glBlendFuncSeparateiARB");
	glBlendFuncSeparateiARB((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]));

}

void gl_minsampleshadingarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMinSampleShadingARB) mogl_glunsupported("glMinSampleShadingARB");
	glMinSampleShadingARB((GLfloat)mxGetScalar(prhs[0]));

}

void gl_namedstringarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNamedStringARB) mogl_glunsupported("glNamedStringARB");
	glNamedStringARB((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(const GLchar*)mxGetData(prhs[4]));

}

void gl_deletenamedstringarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteNamedStringARB) mogl_glunsupported("glDeleteNamedStringARB");
	glDeleteNamedStringARB((GLint)mxGetScalar(prhs[0]),
		(const GLchar*)mxGetData(prhs[1]));

}

void gl_isnamedstringarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsNamedStringARB) mogl_glunsupported("glIsNamedStringARB");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsNamedStringARB((GLint)mxGetScalar(prhs[0]),
		(const GLchar*)mxGetData(prhs[1]));

}

void gl_getnamedstringarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetNamedStringARB) mogl_glunsupported("glGetNamedStringARB");
	glGetNamedStringARB((GLint)mxGetScalar(prhs[0]),
		(const GLchar*)mxGetData(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]),
		(GLchar*)mxGetData(prhs[4]));

}

void gl_getnamedstringivarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetNamedStringivARB) mogl_glunsupported("glGetNamedStringivARB");
	glGetNamedStringivARB((GLint)mxGetScalar(prhs[0]),
		(const GLchar*)mxGetData(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_bindfragdatalocationindexed( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindFragDataLocationIndexed) mogl_glunsupported("glBindFragDataLocationIndexed");
	glBindFragDataLocationIndexed((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(const GLchar*)mxGetData(prhs[3]));

}

void gl_getfragdataindex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetFragDataIndex) mogl_glunsupported("glGetFragDataIndex");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetFragDataIndex((GLuint)mxGetScalar(prhs[0]),
		(const GLchar*)mxGetData(prhs[1]));

}

void gl_gensamplers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenSamplers) mogl_glunsupported("glGenSamplers");
	glGenSamplers((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_deletesamplers( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteSamplers) mogl_glunsupported("glDeleteSamplers");
	glDeleteSamplers((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_issampler( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsSampler) mogl_glunsupported("glIsSampler");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsSampler((GLuint)mxGetScalar(prhs[0]));

}

void gl_bindsampler( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindSampler) mogl_glunsupported("glBindSampler");
	glBindSampler((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_samplerparameteri( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSamplerParameteri) mogl_glunsupported("glSamplerParameteri");
	glSamplerParameteri((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_samplerparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSamplerParameteriv) mogl_glunsupported("glSamplerParameteriv");
	glSamplerParameteriv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_samplerparameterf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSamplerParameterf) mogl_glunsupported("glSamplerParameterf");
	glSamplerParameterf((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_samplerparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSamplerParameterfv) mogl_glunsupported("glSamplerParameterfv");
	glSamplerParameterfv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_samplerparameteriiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSamplerParameterIiv) mogl_glunsupported("glSamplerParameterIiv");
	glSamplerParameterIiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_samplerparameteriuiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSamplerParameterIuiv) mogl_glunsupported("glSamplerParameterIuiv");
	glSamplerParameterIuiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLuint*)mxGetData(prhs[2]));

}

void gl_getsamplerparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetSamplerParameteriv) mogl_glunsupported("glGetSamplerParameteriv");
	glGetSamplerParameteriv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getsamplerparameteriiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetSamplerParameterIiv) mogl_glunsupported("glGetSamplerParameterIiv");
	glGetSamplerParameterIiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getsamplerparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetSamplerParameterfv) mogl_glunsupported("glGetSamplerParameterfv");
	glGetSamplerParameterfv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getsamplerparameteriuiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetSamplerParameterIuiv) mogl_glunsupported("glGetSamplerParameterIuiv");
	glGetSamplerParameterIuiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint*)mxGetData(prhs[2]));

}

void gl_querycounter( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glQueryCounter) mogl_glunsupported("glQueryCounter");
	glQueryCounter((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]));

}

void gl_getqueryobjecti64v( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetQueryObjecti64v) mogl_glunsupported("glGetQueryObjecti64v");
	glGetQueryObjecti64v((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint64*)mxGetData(prhs[2]));

}

void gl_getqueryobjectui64v( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetQueryObjectui64v) mogl_glunsupported("glGetQueryObjectui64v");
	glGetQueryObjectui64v((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint64*)mxGetData(prhs[2]));

}

void gl_vertexp2ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexP2ui) mogl_glunsupported("glVertexP2ui");
	glVertexP2ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_vertexp2uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexP2uiv) mogl_glunsupported("glVertexP2uiv");
	glVertexP2uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_vertexp3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexP3ui) mogl_glunsupported("glVertexP3ui");
	glVertexP3ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_vertexp3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexP3uiv) mogl_glunsupported("glVertexP3uiv");
	glVertexP3uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_vertexp4ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexP4ui) mogl_glunsupported("glVertexP4ui");
	glVertexP4ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_vertexp4uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexP4uiv) mogl_glunsupported("glVertexP4uiv");
	glVertexP4uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_texcoordp1ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordP1ui) mogl_glunsupported("glTexCoordP1ui");
	glTexCoordP1ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_texcoordp1uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordP1uiv) mogl_glunsupported("glTexCoordP1uiv");
	glTexCoordP1uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_texcoordp2ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordP2ui) mogl_glunsupported("glTexCoordP2ui");
	glTexCoordP2ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_texcoordp2uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordP2uiv) mogl_glunsupported("glTexCoordP2uiv");
	glTexCoordP2uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_texcoordp3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordP3ui) mogl_glunsupported("glTexCoordP3ui");
	glTexCoordP3ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_texcoordp3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordP3uiv) mogl_glunsupported("glTexCoordP3uiv");
	glTexCoordP3uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_texcoordp4ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordP4ui) mogl_glunsupported("glTexCoordP4ui");
	glTexCoordP4ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_texcoordp4uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexCoordP4uiv) mogl_glunsupported("glTexCoordP4uiv");
	glTexCoordP4uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_multitexcoordp1ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoordP1ui) mogl_glunsupported("glMultiTexCoordP1ui");
	glMultiTexCoordP1ui((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_multitexcoordp1uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoordP1uiv) mogl_glunsupported("glMultiTexCoordP1uiv");
	glMultiTexCoordP1uiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLuint*)mxGetData(prhs[2]));

}

void gl_multitexcoordp2ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoordP2ui) mogl_glunsupported("glMultiTexCoordP2ui");
	glMultiTexCoordP2ui((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_multitexcoordp2uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoordP2uiv) mogl_glunsupported("glMultiTexCoordP2uiv");
	glMultiTexCoordP2uiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLuint*)mxGetData(prhs[2]));

}

void gl_multitexcoordp3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoordP3ui) mogl_glunsupported("glMultiTexCoordP3ui");
	glMultiTexCoordP3ui((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_multitexcoordp3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoordP3uiv) mogl_glunsupported("glMultiTexCoordP3uiv");
	glMultiTexCoordP3uiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLuint*)mxGetData(prhs[2]));

}

void gl_multitexcoordp4ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoordP4ui) mogl_glunsupported("glMultiTexCoordP4ui");
	glMultiTexCoordP4ui((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_multitexcoordp4uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiTexCoordP4uiv) mogl_glunsupported("glMultiTexCoordP4uiv");
	glMultiTexCoordP4uiv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLuint*)mxGetData(prhs[2]));

}

void gl_normalp3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormalP3ui) mogl_glunsupported("glNormalP3ui");
	glNormalP3ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_normalp3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNormalP3uiv) mogl_glunsupported("glNormalP3uiv");
	glNormalP3uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_colorp3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorP3ui) mogl_glunsupported("glColorP3ui");
	glColorP3ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_colorp3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorP3uiv) mogl_glunsupported("glColorP3uiv");
	glColorP3uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_colorp4ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorP4ui) mogl_glunsupported("glColorP4ui");
	glColorP4ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_colorp4uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glColorP4uiv) mogl_glunsupported("glColorP4uiv");
	glColorP4uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_secondarycolorp3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColorP3ui) mogl_glunsupported("glSecondaryColorP3ui");
	glSecondaryColorP3ui((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_secondarycolorp3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glSecondaryColorP3uiv) mogl_glunsupported("glSecondaryColorP3uiv");
	glSecondaryColorP3uiv((GLenum)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_vertexattribp1ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribP1ui) mogl_glunsupported("glVertexAttribP1ui");
	glVertexAttribP1ui((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_vertexattribp1uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribP1uiv) mogl_glunsupported("glVertexAttribP1uiv");
	glVertexAttribP1uiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLuint*)mxGetData(prhs[3]));

}

void gl_vertexattribp2ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribP2ui) mogl_glunsupported("glVertexAttribP2ui");
	glVertexAttribP2ui((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_vertexattribp2uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribP2uiv) mogl_glunsupported("glVertexAttribP2uiv");
	glVertexAttribP2uiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLuint*)mxGetData(prhs[3]));

}

void gl_vertexattribp3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribP3ui) mogl_glunsupported("glVertexAttribP3ui");
	glVertexAttribP3ui((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_vertexattribp3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribP3uiv) mogl_glunsupported("glVertexAttribP3uiv");
	glVertexAttribP3uiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLuint*)mxGetData(prhs[3]));

}

void gl_vertexattribp4ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribP4ui) mogl_glunsupported("glVertexAttribP4ui");
	glVertexAttribP4ui((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_vertexattribp4uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribP4uiv) mogl_glunsupported("glVertexAttribP4uiv");
	glVertexAttribP4uiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLuint*)mxGetData(prhs[3]));

}

void gl_drawarraysindirect( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawArraysIndirect) mogl_glunsupported("glDrawArraysIndirect");
	glDrawArraysIndirect((GLenum)mxGetScalar(prhs[0]),
		(const GLvoid*)mxGetData(prhs[1]));

}

void gl_drawelementsindirect( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawElementsIndirect) mogl_glunsupported("glDrawElementsIndirect");
	glDrawElementsIndirect((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLvoid*)mxGetData(prhs[2]));

}

void gl_uniform1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform1d) mogl_glunsupported("glUniform1d");
	glUniform1d((GLint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_uniform2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform2d) mogl_glunsupported("glUniform2d");
	glUniform2d((GLint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_uniform3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform3d) mogl_glunsupported("glUniform3d");
	glUniform3d((GLint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_uniform4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform4d) mogl_glunsupported("glUniform4d");
	glUniform4d((GLint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]));

}

void gl_uniform1dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform1dv) mogl_glunsupported("glUniform1dv");
	glUniform1dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLdouble*)mxGetData(prhs[2]));

}

void gl_uniform2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform2dv) mogl_glunsupported("glUniform2dv");
	glUniform2dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLdouble*)mxGetData(prhs[2]));

}

void gl_uniform3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform3dv) mogl_glunsupported("glUniform3dv");
	glUniform3dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLdouble*)mxGetData(prhs[2]));

}

void gl_uniform4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniform4dv) mogl_glunsupported("glUniform4dv");
	glUniform4dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLdouble*)mxGetData(prhs[2]));

}

void gl_uniformmatrix2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix2dv) mogl_glunsupported("glUniformMatrix2dv");
	glUniformMatrix2dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_uniformmatrix3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix3dv) mogl_glunsupported("glUniformMatrix3dv");
	glUniformMatrix3dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_uniformmatrix4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix4dv) mogl_glunsupported("glUniformMatrix4dv");
	glUniformMatrix4dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_uniformmatrix2x3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix2x3dv) mogl_glunsupported("glUniformMatrix2x3dv");
	glUniformMatrix2x3dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_uniformmatrix2x4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix2x4dv) mogl_glunsupported("glUniformMatrix2x4dv");
	glUniformMatrix2x4dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_uniformmatrix3x2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix3x2dv) mogl_glunsupported("glUniformMatrix3x2dv");
	glUniformMatrix3x2dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_uniformmatrix3x4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix3x4dv) mogl_glunsupported("glUniformMatrix3x4dv");
	glUniformMatrix3x4dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_uniformmatrix4x2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix4x2dv) mogl_glunsupported("glUniformMatrix4x2dv");
	glUniformMatrix4x2dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_uniformmatrix4x3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformMatrix4x3dv) mogl_glunsupported("glUniformMatrix4x3dv");
	glUniformMatrix4x3dv((GLint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLboolean)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_getuniformdv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetUniformdv) mogl_glunsupported("glGetUniformdv");
	glGetUniformdv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLdouble*)mxGetData(prhs[2]));

}

void gl_getsubroutineuniformlocation( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetSubroutineUniformLocation) mogl_glunsupported("glGetSubroutineUniformLocation");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetSubroutineUniformLocation((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]));

}

void gl_getsubroutineindex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetSubroutineIndex) mogl_glunsupported("glGetSubroutineIndex");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetSubroutineIndex((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]));

}

void gl_getactivesubroutineuniformiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveSubroutineUniformiv) mogl_glunsupported("glGetActiveSubroutineUniformiv");
	glGetActiveSubroutineUniformiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLint*)mxGetData(prhs[4]));

}

void gl_getactivesubroutineuniformname( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveSubroutineUniformName) mogl_glunsupported("glGetActiveSubroutineUniformName");
	glGetActiveSubroutineUniformName((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei*)mxGetData(prhs[4]),
		(GLchar*)mxGetData(prhs[5]));

}

void gl_getactivesubroutinename( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveSubroutineName) mogl_glunsupported("glGetActiveSubroutineName");
	glGetActiveSubroutineName((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei*)mxGetData(prhs[4]),
		(GLchar*)mxGetData(prhs[5]));

}

void gl_uniformsubroutinesuiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUniformSubroutinesuiv) mogl_glunsupported("glUniformSubroutinesuiv");
	glUniformSubroutinesuiv((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLuint*)mxGetData(prhs[2]));

}

void gl_getuniformsubroutineuiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetUniformSubroutineuiv) mogl_glunsupported("glGetUniformSubroutineuiv");
	glGetUniformSubroutineuiv((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLuint*)mxGetData(prhs[2]));

}

void gl_getprogramstageiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramStageiv) mogl_glunsupported("glGetProgramStageiv");
	glGetProgramStageiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_patchparameteri( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPatchParameteri) mogl_glunsupported("glPatchParameteri");
	glPatchParameteri((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_patchparameterfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPatchParameterfv) mogl_glunsupported("glPatchParameterfv");
	glPatchParameterfv((GLenum)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_bindtransformfeedback( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindTransformFeedback) mogl_glunsupported("glBindTransformFeedback");
	glBindTransformFeedback((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_deletetransformfeedbacks( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteTransformFeedbacks) mogl_glunsupported("glDeleteTransformFeedbacks");
	glDeleteTransformFeedbacks((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_gentransformfeedbacks( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenTransformFeedbacks) mogl_glunsupported("glGenTransformFeedbacks");
	glGenTransformFeedbacks((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_istransformfeedback( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsTransformFeedback) mogl_glunsupported("glIsTransformFeedback");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsTransformFeedback((GLuint)mxGetScalar(prhs[0]));

}

void gl_pausetransformfeedback( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPauseTransformFeedback) mogl_glunsupported("glPauseTransformFeedback");
	glPauseTransformFeedback();

}

void gl_resumetransformfeedback( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glResumeTransformFeedback) mogl_glunsupported("glResumeTransformFeedback");
	glResumeTransformFeedback();

}

void gl_drawtransformfeedback( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawTransformFeedback) mogl_glunsupported("glDrawTransformFeedback");
	glDrawTransformFeedback((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_drawtransformfeedbackstream( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawTransformFeedbackStream) mogl_glunsupported("glDrawTransformFeedbackStream");
	glDrawTransformFeedbackStream((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_beginqueryindexed( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBeginQueryIndexed) mogl_glunsupported("glBeginQueryIndexed");
	glBeginQueryIndexed((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_endqueryindexed( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glEndQueryIndexed) mogl_glunsupported("glEndQueryIndexed");
	glEndQueryIndexed((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_getqueryindexediv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetQueryIndexediv) mogl_glunsupported("glGetQueryIndexediv");
	glGetQueryIndexediv((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_releaseshadercompiler( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glReleaseShaderCompiler) mogl_glunsupported("glReleaseShaderCompiler");
	glReleaseShaderCompiler();

}

void gl_shaderbinary( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glShaderBinary) mogl_glunsupported("glShaderBinary");
	glShaderBinary((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(const GLvoid*)mxGetData(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_getshaderprecisionformat( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetShaderPrecisionFormat) mogl_glunsupported("glGetShaderPrecisionFormat");
	glGetShaderPrecisionFormat((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_depthrangef( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDepthRangef) mogl_glunsupported("glDepthRangef");
	glDepthRangef((GLfloat)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]));

}

void gl_cleardepthf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearDepthf) mogl_glunsupported("glClearDepthf");
	glClearDepthf((GLfloat)mxGetScalar(prhs[0]));

}

void gl_getprogrambinary( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramBinary) mogl_glunsupported("glGetProgramBinary");
	glGetProgramBinary((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei*)mxGetData(prhs[2]),
		(GLenum*)mxGetData(prhs[3]),
		(GLvoid*)mxGetData(prhs[4]));

}

void gl_programbinary( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramBinary) mogl_glunsupported("glProgramBinary");
	glProgramBinary((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLvoid*)mxGetData(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_programparameteri( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramParameteri) mogl_glunsupported("glProgramParameteri");
	glProgramParameteri((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_useprogramstages( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glUseProgramStages) mogl_glunsupported("glUseProgramStages");
	glUseProgramStages((GLuint)mxGetScalar(prhs[0]),
		(GLbitfield)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_activeshaderprogram( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glActiveShaderProgram) mogl_glunsupported("glActiveShaderProgram");
	glActiveShaderProgram((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_createshaderprogramv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCreateShaderProgramv) mogl_glunsupported("glCreateShaderProgramv");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glCreateShaderProgramv((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]));

}

void gl_bindprogrampipeline( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindProgramPipeline) mogl_glunsupported("glBindProgramPipeline");
	glBindProgramPipeline((GLuint)mxGetScalar(prhs[0]));

}

void gl_deleteprogrampipelines( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDeleteProgramPipelines) mogl_glunsupported("glDeleteProgramPipelines");
	glDeleteProgramPipelines((GLsizei)mxGetScalar(prhs[0]),
		(const GLuint*)mxGetData(prhs[1]));

}

void gl_genprogrampipelines( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGenProgramPipelines) mogl_glunsupported("glGenProgramPipelines");
	glGenProgramPipelines((GLsizei)mxGetScalar(prhs[0]),
		(GLuint*)mxGetData(prhs[1]));

}

void gl_isprogrampipeline( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glIsProgramPipeline) mogl_glunsupported("glIsProgramPipeline");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glIsProgramPipeline((GLuint)mxGetScalar(prhs[0]));

}

void gl_getprogrampipelineiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramPipelineiv) mogl_glunsupported("glGetProgramPipelineiv");
	glGetProgramPipelineiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_programuniform1i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform1i) mogl_glunsupported("glProgramUniform1i");
	glProgramUniform1i((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_programuniform1iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform1iv) mogl_glunsupported("glProgramUniform1iv");
	glProgramUniform1iv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLint*)mxGetData(prhs[3]));

}

void gl_programuniform1f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform1f) mogl_glunsupported("glProgramUniform1f");
	glProgramUniform1f((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]));

}

void gl_programuniform1fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform1fv) mogl_glunsupported("glProgramUniform1fv");
	glProgramUniform1fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_programuniform1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform1d) mogl_glunsupported("glProgramUniform1d");
	glProgramUniform1d((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_programuniform1dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform1dv) mogl_glunsupported("glProgramUniform1dv");
	glProgramUniform1dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_programuniform1ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform1ui) mogl_glunsupported("glProgramUniform1ui");
	glProgramUniform1ui((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_programuniform1uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform1uiv) mogl_glunsupported("glProgramUniform1uiv");
	glProgramUniform1uiv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLuint*)mxGetData(prhs[3]));

}

void gl_programuniform2i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform2i) mogl_glunsupported("glProgramUniform2i");
	glProgramUniform2i((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]));

}

void gl_programuniform2iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform2iv) mogl_glunsupported("glProgramUniform2iv");
	glProgramUniform2iv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLint*)mxGetData(prhs[3]));

}

void gl_programuniform2f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform2f) mogl_glunsupported("glProgramUniform2f");
	glProgramUniform2f((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]));

}

void gl_programuniform2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform2fv) mogl_glunsupported("glProgramUniform2fv");
	glProgramUniform2fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_programuniform2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform2d) mogl_glunsupported("glProgramUniform2d");
	glProgramUniform2d((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_programuniform2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform2dv) mogl_glunsupported("glProgramUniform2dv");
	glProgramUniform2dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_programuniform2ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform2ui) mogl_glunsupported("glProgramUniform2ui");
	glProgramUniform2ui((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_programuniform2uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform2uiv) mogl_glunsupported("glProgramUniform2uiv");
	glProgramUniform2uiv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLuint*)mxGetData(prhs[3]));

}

void gl_programuniform3i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform3i) mogl_glunsupported("glProgramUniform3i");
	glProgramUniform3i((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]));

}

void gl_programuniform3iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform3iv) mogl_glunsupported("glProgramUniform3iv");
	glProgramUniform3iv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLint*)mxGetData(prhs[3]));

}

void gl_programuniform3f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform3f) mogl_glunsupported("glProgramUniform3f");
	glProgramUniform3f((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]),
		(GLfloat)mxGetScalar(prhs[4]));

}

void gl_programuniform3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform3fv) mogl_glunsupported("glProgramUniform3fv");
	glProgramUniform3fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_programuniform3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform3d) mogl_glunsupported("glProgramUniform3d");
	glProgramUniform3d((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]));

}

void gl_programuniform3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform3dv) mogl_glunsupported("glProgramUniform3dv");
	glProgramUniform3dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_programuniform3ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform3ui) mogl_glunsupported("glProgramUniform3ui");
	glProgramUniform3ui((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLuint)mxGetScalar(prhs[4]));

}

void gl_programuniform3uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform3uiv) mogl_glunsupported("glProgramUniform3uiv");
	glProgramUniform3uiv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLuint*)mxGetData(prhs[3]));

}

void gl_programuniform4i( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform4i) mogl_glunsupported("glProgramUniform4i");
	glProgramUniform4i((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]));

}

void gl_programuniform4iv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform4iv) mogl_glunsupported("glProgramUniform4iv");
	glProgramUniform4iv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLint*)mxGetData(prhs[3]));

}

void gl_programuniform4f( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform4f) mogl_glunsupported("glProgramUniform4f");
	glProgramUniform4f((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]),
		(GLfloat)mxGetScalar(prhs[4]),
		(GLfloat)mxGetScalar(prhs[5]));

}

void gl_programuniform4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform4fv) mogl_glunsupported("glProgramUniform4fv");
	glProgramUniform4fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLfloat*)mxGetData(prhs[3]));

}

void gl_programuniform4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform4d) mogl_glunsupported("glProgramUniform4d");
	glProgramUniform4d((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]),
		(GLdouble)mxGetScalar(prhs[5]));

}

void gl_programuniform4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform4dv) mogl_glunsupported("glProgramUniform4dv");
	glProgramUniform4dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLdouble*)mxGetData(prhs[3]));

}

void gl_programuniform4ui( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform4ui) mogl_glunsupported("glProgramUniform4ui");
	glProgramUniform4ui((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLuint)mxGetScalar(prhs[4]),
		(GLuint)mxGetScalar(prhs[5]));

}

void gl_programuniform4uiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniform4uiv) mogl_glunsupported("glProgramUniform4uiv");
	glProgramUniform4uiv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLuint*)mxGetData(prhs[3]));

}

void gl_programuniformmatrix2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix2fv) mogl_glunsupported("glProgramUniformMatrix2fv");
	glProgramUniformMatrix2fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix3fv) mogl_glunsupported("glProgramUniformMatrix3fv");
	glProgramUniformMatrix3fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix4fv) mogl_glunsupported("glProgramUniformMatrix4fv");
	glProgramUniformMatrix4fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix2dv) mogl_glunsupported("glProgramUniformMatrix2dv");
	glProgramUniformMatrix2dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix3dv) mogl_glunsupported("glProgramUniformMatrix3dv");
	glProgramUniformMatrix3dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix4dv) mogl_glunsupported("glProgramUniformMatrix4dv");
	glProgramUniformMatrix4dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix2x3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix2x3fv) mogl_glunsupported("glProgramUniformMatrix2x3fv");
	glProgramUniformMatrix2x3fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix3x2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix3x2fv) mogl_glunsupported("glProgramUniformMatrix3x2fv");
	glProgramUniformMatrix3x2fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix2x4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix2x4fv) mogl_glunsupported("glProgramUniformMatrix2x4fv");
	glProgramUniformMatrix2x4fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix4x2fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix4x2fv) mogl_glunsupported("glProgramUniformMatrix4x2fv");
	glProgramUniformMatrix4x2fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix3x4fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix3x4fv) mogl_glunsupported("glProgramUniformMatrix3x4fv");
	glProgramUniformMatrix3x4fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix4x3fv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix4x3fv) mogl_glunsupported("glProgramUniformMatrix4x3fv");
	glProgramUniformMatrix4x3fv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLfloat*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix2x3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix2x3dv) mogl_glunsupported("glProgramUniformMatrix2x3dv");
	glProgramUniformMatrix2x3dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix3x2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix3x2dv) mogl_glunsupported("glProgramUniformMatrix3x2dv");
	glProgramUniformMatrix3x2dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix2x4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix2x4dv) mogl_glunsupported("glProgramUniformMatrix2x4dv");
	glProgramUniformMatrix2x4dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix4x2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix4x2dv) mogl_glunsupported("glProgramUniformMatrix4x2dv");
	glProgramUniformMatrix4x2dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix3x4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix3x4dv) mogl_glunsupported("glProgramUniformMatrix3x4dv");
	glProgramUniformMatrix3x4dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_programuniformmatrix4x3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glProgramUniformMatrix4x3dv) mogl_glunsupported("glProgramUniformMatrix4x3dv");
	glProgramUniformMatrix4x3dv((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(const GLdouble*)mxGetData(prhs[4]));

}

void gl_validateprogrampipeline( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glValidateProgramPipeline) mogl_glunsupported("glValidateProgramPipeline");
	glValidateProgramPipeline((GLuint)mxGetScalar(prhs[0]));

}

void gl_getprogrampipelineinfolog( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramPipelineInfoLog) mogl_glunsupported("glGetProgramPipelineInfoLog");
	glGetProgramPipelineInfoLog((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei*)mxGetData(prhs[2]),
		(GLchar*)mxGetData(prhs[3]));

}

void gl_vertexattribl1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribL1d) mogl_glunsupported("glVertexAttribL1d");
	glVertexAttribL1d((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]));

}

void gl_vertexattribl2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribL2d) mogl_glunsupported("glVertexAttribL2d");
	glVertexAttribL2d((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_vertexattribl3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribL3d) mogl_glunsupported("glVertexAttribL3d");
	glVertexAttribL3d((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]));

}

void gl_vertexattribl4d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribL4d) mogl_glunsupported("glVertexAttribL4d");
	glVertexAttribL4d((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]),
		(GLdouble)mxGetScalar(prhs[3]),
		(GLdouble)mxGetScalar(prhs[4]));

}

void gl_vertexattribl1dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribL1dv) mogl_glunsupported("glVertexAttribL1dv");
	glVertexAttribL1dv((GLuint)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_vertexattribl2dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribL2dv) mogl_glunsupported("glVertexAttribL2dv");
	glVertexAttribL2dv((GLuint)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_vertexattribl3dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribL3dv) mogl_glunsupported("glVertexAttribL3dv");
	glVertexAttribL3dv((GLuint)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_vertexattribl4dv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribL4dv) mogl_glunsupported("glVertexAttribL4dv");
	glVertexAttribL4dv((GLuint)mxGetScalar(prhs[0]),
		(const GLdouble*)mxGetData(prhs[1]));

}

void gl_vertexattriblpointer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribLPointer) mogl_glunsupported("glVertexAttribLPointer");
	glVertexAttribLPointer((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(const GLvoid*)mxGetData(prhs[4]));

}

void gl_getvertexattribldv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetVertexAttribLdv) mogl_glunsupported("glGetVertexAttribLdv");
	glGetVertexAttribLdv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLdouble*)mxGetData(prhs[2]));

}

void gl_viewportarrayv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glViewportArrayv) mogl_glunsupported("glViewportArrayv");
	glViewportArrayv((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLfloat*)mxGetData(prhs[2]));

}

void gl_viewportindexedf( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glViewportIndexedf) mogl_glunsupported("glViewportIndexedf");
	glViewportIndexedf((GLuint)mxGetScalar(prhs[0]),
		(GLfloat)mxGetScalar(prhs[1]),
		(GLfloat)mxGetScalar(prhs[2]),
		(GLfloat)mxGetScalar(prhs[3]),
		(GLfloat)mxGetScalar(prhs[4]));

}

void gl_viewportindexedfv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glViewportIndexedfv) mogl_glunsupported("glViewportIndexedfv");
	glViewportIndexedfv((GLuint)mxGetScalar(prhs[0]),
		(const GLfloat*)mxGetData(prhs[1]));

}

void gl_scissorarrayv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glScissorArrayv) mogl_glunsupported("glScissorArrayv");
	glScissorArrayv((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLint*)mxGetData(prhs[2]));

}

void gl_scissorindexed( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glScissorIndexed) mogl_glunsupported("glScissorIndexed");
	glScissorIndexed((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_scissorindexedv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glScissorIndexedv) mogl_glunsupported("glScissorIndexedv");
	glScissorIndexedv((GLuint)mxGetScalar(prhs[0]),
		(const GLint*)mxGetData(prhs[1]));

}

void gl_depthrangearrayv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDepthRangeArrayv) mogl_glunsupported("glDepthRangeArrayv");
	glDepthRangeArrayv((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLdouble*)mxGetData(prhs[2]));

}

void gl_depthrangeindexed( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDepthRangeIndexed) mogl_glunsupported("glDepthRangeIndexed");
	glDepthRangeIndexed((GLuint)mxGetScalar(prhs[0]),
		(GLdouble)mxGetScalar(prhs[1]),
		(GLdouble)mxGetScalar(prhs[2]));

}

void gl_getfloati_v( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetFloati_v) mogl_glunsupported("glGetFloati_v");
	glGetFloati_v((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getdoublei_v( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetDoublei_v) mogl_glunsupported("glGetDoublei_v");
	glGetDoublei_v((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLdouble*)mxGetData(prhs[2]));

}

void gl_debugmessagecontrolarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDebugMessageControlARB) mogl_glunsupported("glDebugMessageControlARB");
	glDebugMessageControlARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(const GLuint*)mxGetData(prhs[4]),
		(GLboolean)mxGetScalar(prhs[5]));

}

void gl_debugmessageinsertarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDebugMessageInsertARB) mogl_glunsupported("glDebugMessageInsertARB");
	glDebugMessageInsertARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(const GLchar*)mxGetData(prhs[5]));

}

void gl_getdebugmessagelogarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetDebugMessageLogARB) mogl_glunsupported("glGetDebugMessageLogARB");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetDebugMessageLogARB((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum*)mxGetData(prhs[2]),
		(GLenum*)mxGetData(prhs[3]),
		(GLuint*)mxGetData(prhs[4]),
		(GLenum*)mxGetData(prhs[5]),
		(GLsizei*)mxGetData(prhs[6]),
		(GLchar*)mxGetData(prhs[7]));

}

void gl_getgraphicsresetstatusarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetGraphicsResetStatusARB) mogl_glunsupported("glGetGraphicsResetStatusARB");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetGraphicsResetStatusARB();

}

void gl_getnmapdvarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnMapdvARB) mogl_glunsupported("glGetnMapdvARB");
	glGetnMapdvARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLdouble*)mxGetData(prhs[3]));

}

void gl_getnmapfvarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnMapfvARB) mogl_glunsupported("glGetnMapfvARB");
	glGetnMapfvARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLfloat*)mxGetData(prhs[3]));

}

void gl_getnmapivarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnMapivARB) mogl_glunsupported("glGetnMapivARB");
	glGetnMapivARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_getnpixelmapfvarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnPixelMapfvARB) mogl_glunsupported("glGetnPixelMapfvARB");
	glGetnPixelMapfvARB((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLfloat*)mxGetData(prhs[2]));

}

void gl_getnpixelmapuivarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnPixelMapuivARB) mogl_glunsupported("glGetnPixelMapuivARB");
	glGetnPixelMapuivARB((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLuint*)mxGetData(prhs[2]));

}

void gl_getnpixelmapusvarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnPixelMapusvARB) mogl_glunsupported("glGetnPixelMapusvARB");
	glGetnPixelMapusvARB((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLushort*)mxGetData(prhs[2]));

}

void gl_getnpolygonstipplearb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnPolygonStippleARB) mogl_glunsupported("glGetnPolygonStippleARB");
	glGetnPolygonStippleARB((GLsizei)mxGetScalar(prhs[0]),
		(GLubyte*)mxGetData(prhs[1]));

}

void gl_getncolortablearb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnColorTableARB) mogl_glunsupported("glGetnColorTableARB");
	glGetnColorTableARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLvoid*)mxGetData(prhs[4]));

}

void gl_getnconvolutionfilterarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnConvolutionFilterARB) mogl_glunsupported("glGetnConvolutionFilterARB");
	glGetnConvolutionFilterARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLvoid*)mxGetData(prhs[4]));

}

void gl_getnseparablefilterarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnSeparableFilterARB) mogl_glunsupported("glGetnSeparableFilterARB");
	glGetnSeparableFilterARB((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLvoid*)mxGetData(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLvoid*)mxGetData(prhs[6]),
		(GLvoid*)mxGetData(prhs[7]));

}

void gl_getnhistogramarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnHistogramARB) mogl_glunsupported("glGetnHistogramARB");
	glGetnHistogramARB((GLenum)mxGetScalar(prhs[0]),
		(GLboolean)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLvoid*)mxGetData(prhs[5]));

}

void gl_getnminmaxarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnMinmaxARB) mogl_glunsupported("glGetnMinmaxARB");
	glGetnMinmaxARB((GLenum)mxGetScalar(prhs[0]),
		(GLboolean)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLvoid*)mxGetData(prhs[5]));

}

void gl_getnteximagearb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnTexImageARB) mogl_glunsupported("glGetnTexImageARB");
	glGetnTexImageARB((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLvoid*)mxGetData(prhs[5]));

}

void gl_readnpixelsarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glReadnPixelsARB) mogl_glunsupported("glReadnPixelsARB");
	glReadnPixelsARB((GLint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(GLvoid*)mxGetData(prhs[7]));

}

void gl_getncompressedteximagearb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnCompressedTexImageARB) mogl_glunsupported("glGetnCompressedTexImageARB");
	glGetnCompressedTexImageARB((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLvoid*)mxGetData(prhs[3]));

}

void gl_getnuniformfvarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnUniformfvARB) mogl_glunsupported("glGetnUniformfvARB");
	glGetnUniformfvARB((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLfloat*)mxGetData(prhs[3]));

}

void gl_getnuniformivarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnUniformivARB) mogl_glunsupported("glGetnUniformivARB");
	glGetnUniformivARB((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_getnuniformuivarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnUniformuivARB) mogl_glunsupported("glGetnUniformuivARB");
	glGetnUniformuivARB((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLuint*)mxGetData(prhs[3]));

}

void gl_getnuniformdvarb( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetnUniformdvARB) mogl_glunsupported("glGetnUniformdvARB");
	glGetnUniformdvARB((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLdouble*)mxGetData(prhs[3]));

}

void gl_drawarraysinstancedbaseinstance( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawArraysInstancedBaseInstance) mogl_glunsupported("glDrawArraysInstancedBaseInstance");
	glDrawArraysInstancedBaseInstance((GLenum)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLuint)mxGetScalar(prhs[4]));

}

void gl_drawelementsinstancedbaseinstance( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawElementsInstancedBaseInstance) mogl_glunsupported("glDrawElementsInstancedBaseInstance");
	glDrawElementsInstancedBaseInstance((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(const void*)mxGetData(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLuint)mxGetScalar(prhs[5]));

}

void gl_drawelementsinstancedbasevertexbaseinstance( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawElementsInstancedBaseVertexBaseInstance) mogl_glunsupported("glDrawElementsInstancedBaseVertexBaseInstance");
	glDrawElementsInstancedBaseVertexBaseInstance((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(const void*)mxGetData(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLuint)mxGetScalar(prhs[6]));

}

void gl_drawtransformfeedbackinstanced( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawTransformFeedbackInstanced) mogl_glunsupported("glDrawTransformFeedbackInstanced");
	glDrawTransformFeedbackInstanced((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]));

}

void gl_drawtransformfeedbackstreaminstanced( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDrawTransformFeedbackStreamInstanced) mogl_glunsupported("glDrawTransformFeedbackStreamInstanced");
	glDrawTransformFeedbackStreamInstanced((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_getinternalformativ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetInternalformativ) mogl_glunsupported("glGetInternalformativ");
	glGetInternalformativ((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLint*)mxGetData(prhs[4]));

}

void gl_getactiveatomiccounterbufferiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetActiveAtomicCounterBufferiv) mogl_glunsupported("glGetActiveAtomicCounterBufferiv");
	glGetActiveAtomicCounterBufferiv((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_bindimagetexture( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindImageTexture) mogl_glunsupported("glBindImageTexture");
	glBindImageTexture((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(GLenum)mxGetScalar(prhs[6]));

}

void gl_memorybarrier( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMemoryBarrier) mogl_glunsupported("glMemoryBarrier");
	glMemoryBarrier((GLbitfield)mxGetScalar(prhs[0]));

}

void gl_texstorage1d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexStorage1D) mogl_glunsupported("glTexStorage1D");
	glTexStorage1D((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_texstorage2d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexStorage2D) mogl_glunsupported("glTexStorage2D");
	glTexStorage2D((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_texstorage3d( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexStorage3D) mogl_glunsupported("glTexStorage3D");
	glTexStorage3D((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]));

}

void gl_texturestorage1dext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTextureStorage1DEXT) mogl_glunsupported("glTextureStorage1DEXT");
	glTextureStorage1DEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_texturestorage2dext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTextureStorage2DEXT) mogl_glunsupported("glTextureStorage2DEXT");
	glTextureStorage2DEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]));

}

void gl_texturestorage3dext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTextureStorage3DEXT) mogl_glunsupported("glTextureStorage3DEXT");
	glTextureStorage3DEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]));

}

void gl_debugmessagecontrol( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDebugMessageControl) mogl_glunsupported("glDebugMessageControl");
	glDebugMessageControl((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(const GLuint*)mxGetData(prhs[4]),
		(GLboolean)mxGetScalar(prhs[5]));

}

void gl_debugmessageinsert( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDebugMessageInsert) mogl_glunsupported("glDebugMessageInsert");
	glDebugMessageInsert((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(const GLchar*)mxGetData(prhs[5]));

}

void gl_getdebugmessagelog( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetDebugMessageLog) mogl_glunsupported("glGetDebugMessageLog");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetDebugMessageLog((GLuint)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum*)mxGetData(prhs[2]),
		(GLenum*)mxGetData(prhs[3]),
		(GLuint*)mxGetData(prhs[4]),
		(GLenum*)mxGetData(prhs[5]),
		(GLsizei*)mxGetData(prhs[6]),
		(GLchar*)mxGetData(prhs[7]));

}

void gl_pushdebuggroup( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glPushDebugGroup) mogl_glunsupported("glPushDebugGroup");
	glPushDebugGroup((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLchar*)mxGetData(prhs[3]));

}

void gl_objectlabel( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glObjectLabel) mogl_glunsupported("glObjectLabel");
	glObjectLabel((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(const GLchar*)mxGetData(prhs[3]));

}

void gl_getobjectlabel( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetObjectLabel) mogl_glunsupported("glGetObjectLabel");
	glGetObjectLabel((GLenum)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei*)mxGetData(prhs[3]),
		(GLchar*)mxGetData(prhs[4]));

}

void gl_objectptrlabel( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glObjectPtrLabel) mogl_glunsupported("glObjectPtrLabel");
	glObjectPtrLabel((const void*)mxGetData(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]));

}

void gl_getobjectptrlabel( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetObjectPtrLabel) mogl_glunsupported("glGetObjectPtrLabel");
	glGetObjectPtrLabel((const void*)mxGetData(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLsizei*)mxGetData(prhs[2]),
		(GLchar*)mxGetData(prhs[3]));

}

void gl_clearbufferdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearBufferData) mogl_glunsupported("glClearBufferData");
	glClearBufferData((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(const void*)mxGetData(prhs[4]));

}

void gl_clearbuffersubdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearBufferSubData) mogl_glunsupported("glClearBufferSubData");
	glClearBufferSubData((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLenum)mxGetScalar(prhs[4]),
		(GLenum)mxGetScalar(prhs[5]),
		(const void*)mxGetData(prhs[6]));

}

void gl_clearnamedbufferdataext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearNamedBufferDataEXT) mogl_glunsupported("glClearNamedBufferDataEXT");
	glClearNamedBufferDataEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(const void*)mxGetData(prhs[4]));

}

void gl_clearnamedbuffersubdataext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glClearNamedBufferSubDataEXT) mogl_glunsupported("glClearNamedBufferSubDataEXT");
	glClearNamedBufferSubDataEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(const void*)mxGetData(prhs[6]));

}

void gl_dispatchcompute( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDispatchCompute) mogl_glunsupported("glDispatchCompute");
	glDispatchCompute((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_dispatchcomputeindirect( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glDispatchComputeIndirect) mogl_glunsupported("glDispatchComputeIndirect");
	glDispatchComputeIndirect((GLint)mxGetScalar(prhs[0]));

}

void gl_copyimagesubdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glCopyImageSubData) mogl_glunsupported("glCopyImageSubData");
	glCopyImageSubData((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLint)mxGetScalar(prhs[5]),
		(GLuint)mxGetScalar(prhs[6]),
		(GLenum)mxGetScalar(prhs[7]),
		(GLint)mxGetScalar(prhs[8]),
		(GLint)mxGetScalar(prhs[9]),
		(GLint)mxGetScalar(prhs[10]),
		(GLint)mxGetScalar(prhs[11]),
		(GLsizei)mxGetScalar(prhs[12]),
		(GLsizei)mxGetScalar(prhs[13]),
		(GLsizei)mxGetScalar(prhs[14]));

}

void gl_textureview( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTextureView) mogl_glunsupported("glTextureView");
	glTextureView((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLuint)mxGetScalar(prhs[4]),
		(GLuint)mxGetScalar(prhs[5]),
		(GLuint)mxGetScalar(prhs[6]),
		(GLuint)mxGetScalar(prhs[7]));

}

void gl_bindvertexbuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glBindVertexBuffer) mogl_glunsupported("glBindVertexBuffer");
	glBindVertexBuffer((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_vertexattribformat( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribFormat) mogl_glunsupported("glVertexAttribFormat");
	glVertexAttribFormat((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLboolean)mxGetScalar(prhs[3]),
		(GLuint)mxGetScalar(prhs[4]));

}

void gl_vertexattribiformat( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribIFormat) mogl_glunsupported("glVertexAttribIFormat");
	glVertexAttribIFormat((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_vertexattriblformat( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribLFormat) mogl_glunsupported("glVertexAttribLFormat");
	glVertexAttribLFormat((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]));

}

void gl_vertexattribbinding( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexAttribBinding) mogl_glunsupported("glVertexAttribBinding");
	glVertexAttribBinding((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_vertexbindingdivisor( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glVertexBindingDivisor) mogl_glunsupported("glVertexBindingDivisor");
	glVertexBindingDivisor((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]));

}

void gl_framebufferparameteri( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glFramebufferParameteri) mogl_glunsupported("glFramebufferParameteri");
	glFramebufferParameteri((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_getframebufferparameteriv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetFramebufferParameteriv) mogl_glunsupported("glGetFramebufferParameteriv");
	glGetFramebufferParameteriv((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_namedframebufferparameteriext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glNamedFramebufferParameteriEXT) mogl_glunsupported("glNamedFramebufferParameteriEXT");
	glNamedFramebufferParameteriEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]));

}

void gl_getnamedframebufferparameterivext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetNamedFramebufferParameterivEXT) mogl_glunsupported("glGetNamedFramebufferParameterivEXT");
	glGetNamedFramebufferParameterivEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLint*)mxGetData(prhs[2]));

}

void gl_getinternalformati64v( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetInternalformati64v) mogl_glunsupported("glGetInternalformati64v");
	glGetInternalformati64v((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLint64*)mxGetData(prhs[4]));

}

void gl_invalidatetexsubimage( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glInvalidateTexSubImage) mogl_glunsupported("glInvalidateTexSubImage");
	glInvalidateTexSubImage((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(GLsizei)mxGetScalar(prhs[7]));

}

void gl_invalidateteximage( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glInvalidateTexImage) mogl_glunsupported("glInvalidateTexImage");
	glInvalidateTexImage((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]));

}

void gl_invalidatebuffersubdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glInvalidateBufferSubData) mogl_glunsupported("glInvalidateBufferSubData");
	glInvalidateBufferSubData((GLuint)mxGetScalar(prhs[0]),
		(GLint)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]));

}

void gl_invalidatebufferdata( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glInvalidateBufferData) mogl_glunsupported("glInvalidateBufferData");
	glInvalidateBufferData((GLuint)mxGetScalar(prhs[0]));

}

void gl_invalidateframebuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glInvalidateFramebuffer) mogl_glunsupported("glInvalidateFramebuffer");
	glInvalidateFramebuffer((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLenum*)mxGetData(prhs[2]));

}

void gl_invalidatesubframebuffer( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glInvalidateSubFramebuffer) mogl_glunsupported("glInvalidateSubFramebuffer");
	glInvalidateSubFramebuffer((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(const GLenum*)mxGetData(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]));

}

void gl_multidrawarraysindirect( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiDrawArraysIndirect) mogl_glunsupported("glMultiDrawArraysIndirect");
	glMultiDrawArraysIndirect((GLenum)mxGetScalar(prhs[0]),
		(const void*)mxGetData(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]));

}

void gl_multidrawelementsindirect( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glMultiDrawElementsIndirect) mogl_glunsupported("glMultiDrawElementsIndirect");
	glMultiDrawElementsIndirect((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const void*)mxGetData(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_getprograminterfaceiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramInterfaceiv) mogl_glunsupported("glGetProgramInterfaceiv");
	glGetProgramInterfaceiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLint*)mxGetData(prhs[3]));

}

void gl_getprogramresourceindex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramResourceIndex) mogl_glunsupported("glGetProgramResourceIndex");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetProgramResourceIndex((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]));

}

void gl_getprogramresourcename( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramResourceName) mogl_glunsupported("glGetProgramResourceName");
	glGetProgramResourceName((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei*)mxGetData(prhs[4]),
		(GLchar*)mxGetData(prhs[5]));

}

void gl_getprogramresourceiv( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramResourceiv) mogl_glunsupported("glGetProgramResourceiv");
	glGetProgramResourceiv((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(const GLenum*)mxGetData(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei*)mxGetData(prhs[6]),
		(GLint*)mxGetData(prhs[7]));

}

void gl_getprogramresourcelocation( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramResourceLocation) mogl_glunsupported("glGetProgramResourceLocation");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetProgramResourceLocation((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]));

}

void gl_getprogramresourcelocationindex( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glGetProgramResourceLocationIndex) mogl_glunsupported("glGetProgramResourceLocationIndex");
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[0])=(double)glGetProgramResourceLocationIndex((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(const GLchar*)mxGetData(prhs[2]));

}

void gl_shaderstorageblockbinding( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glShaderStorageBlockBinding) mogl_glunsupported("glShaderStorageBlockBinding");
	glShaderStorageBlockBinding((GLuint)mxGetScalar(prhs[0]),
		(GLuint)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]));

}

void gl_texbufferrange( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexBufferRange) mogl_glunsupported("glTexBufferRange");
	glTexBufferRange((GLenum)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLuint)mxGetScalar(prhs[2]),
		(GLint)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]));

}

void gl_texturebufferrangeext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTextureBufferRangeEXT) mogl_glunsupported("glTextureBufferRangeEXT");
	glTextureBufferRangeEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLuint)mxGetScalar(prhs[3]),
		(GLint)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]));

}

void gl_texstorage2dmultisample( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexStorage2DMultisample) mogl_glunsupported("glTexStorage2DMultisample");
	glTexStorage2DMultisample((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLboolean)mxGetScalar(prhs[5]));

}

void gl_texstorage3dmultisample( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTexStorage3DMultisample) mogl_glunsupported("glTexStorage3DMultisample");
	glTexStorage3DMultisample((GLenum)mxGetScalar(prhs[0]),
		(GLsizei)mxGetScalar(prhs[1]),
		(GLenum)mxGetScalar(prhs[2]),
		(GLsizei)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLboolean)mxGetScalar(prhs[6]));

}

void gl_texturestorage2dmultisampleext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTextureStorage2DMultisampleEXT) mogl_glunsupported("glTextureStorage2DMultisampleEXT");
	glTextureStorage2DMultisampleEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLboolean)mxGetScalar(prhs[6]));

}

void gl_texturestorage3dmultisampleext( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	if (NULL == glTextureStorage3DMultisampleEXT) mogl_glunsupported("glTextureStorage3DMultisampleEXT");
	glTextureStorage3DMultisampleEXT((GLuint)mxGetScalar(prhs[0]),
		(GLenum)mxGetScalar(prhs[1]),
		(GLsizei)mxGetScalar(prhs[2]),
		(GLenum)mxGetScalar(prhs[3]),
		(GLsizei)mxGetScalar(prhs[4]),
		(GLsizei)mxGetScalar(prhs[5]),
		(GLsizei)mxGetScalar(prhs[6]),
		(GLboolean)mxGetScalar(prhs[7]));

}

int gl_auto_map_count=982;
cmdhandler gl_auto_map[] = {
{ "glAccum",                         gl_accum                            },
{ "glActiveShaderProgram",           gl_activeshaderprogram              },
{ "glActiveTexture",                 gl_activetexture                    },
{ "glAlphaFunc",                     gl_alphafunc                        },
{ "glAreTexturesResident",           gl_aretexturesresident              },
{ "glArrayElement",                  gl_arrayelement                     },
{ "glAttachShader",                  gl_attachshader                     },
{ "glBegin",                         gl_begin                            },
{ "glBeginQuery",                    gl_beginquery                       },
{ "glBeginQueryIndexed",             gl_beginqueryindexed                },
{ "glBindAttribLocation",            gl_bindattriblocation               },
{ "glBindBuffer",                    gl_bindbuffer                       },
{ "glBindFragDataLocationIndexed",   gl_bindfragdatalocationindexed      },
{ "glBindFramebuffer",               gl_bindframebuffer                  },
{ "glBindFramebufferEXT",            gl_bindframebufferext               },
{ "glBindImageTexture",              gl_bindimagetexture                 },
{ "glBindProgramPipeline",           gl_bindprogrampipeline              },
{ "glBindRenderbuffer",              gl_bindrenderbuffer                 },
{ "glBindRenderbufferEXT",           gl_bindrenderbufferext              },
{ "glBindSampler",                   gl_bindsampler                      },
{ "glBindTexture",                   gl_bindtexture                      },
{ "glBindTransformFeedback",         gl_bindtransformfeedback            },
{ "glBindVertexArray",               gl_bindvertexarray                  },
{ "glBindVertexBuffer",              gl_bindvertexbuffer                 },
{ "glBitmap",                        gl_bitmap                           },
{ "glBlendColor",                    gl_blendcolor                       },
{ "glBlendEquation",                 gl_blendequation                    },
{ "glBlendEquationSeparate",         gl_blendequationseparate            },
{ "glBlendEquationSeparateiARB",     gl_blendequationseparateiarb        },
{ "glBlendEquationiARB",             gl_blendequationiarb                },
{ "glBlendFunc",                     gl_blendfunc                        },
{ "glBlendFuncSeparate",             gl_blendfuncseparate                },
{ "glBlendFuncSeparateiARB",         gl_blendfuncseparateiarb            },
{ "glBlendFunciARB",                 gl_blendfunciarb                    },
{ "glBlitFramebuffer",               gl_blitframebuffer                  },
{ "glBufferSubData",                 gl_buffersubdata                    },
{ "glCallList",                      gl_calllist                         },
{ "glCallLists",                     gl_calllists                        },
{ "glCheckFramebufferStatus",        gl_checkframebufferstatus           },
{ "glCheckFramebufferStatusEXT",     gl_checkframebufferstatusext        },
{ "glClear",                         gl_clear                            },
{ "glClearAccum",                    gl_clearaccum                       },
{ "glClearBufferData",               gl_clearbufferdata                  },
{ "glClearBufferSubData",            gl_clearbuffersubdata               },
{ "glClearColor",                    gl_clearcolor                       },
{ "glClearDepth",                    gl_cleardepth                       },
{ "glClearDepthf",                   gl_cleardepthf                      },
{ "glClearIndex",                    gl_clearindex                       },
{ "glClearNamedBufferDataEXT",       gl_clearnamedbufferdataext          },
{ "glClearNamedBufferSubDataEXT",    gl_clearnamedbuffersubdataext       },
{ "glClearStencil",                  gl_clearstencil                     },
{ "glClientActiveTexture",           gl_clientactivetexture              },
{ "glClientWaitSync",                gl_clientwaitsync                   },
{ "glClipPlane",                     gl_clipplane                        },
{ "glColor3b",                       gl_color3b                          },
{ "glColor3bv",                      gl_color3bv                         },
{ "glColor3d",                       gl_color3d                          },
{ "glColor3dv",                      gl_color3dv                         },
{ "glColor3f",                       gl_color3f                          },
{ "glColor3fv",                      gl_color3fv                         },
{ "glColor3i",                       gl_color3i                          },
{ "glColor3iv",                      gl_color3iv                         },
{ "glColor3s",                       gl_color3s                          },
{ "glColor3sv",                      gl_color3sv                         },
{ "glColor3ub",                      gl_color3ub                         },
{ "glColor3ubv",                     gl_color3ubv                        },
{ "glColor3ui",                      gl_color3ui                         },
{ "glColor3uiv",                     gl_color3uiv                        },
{ "glColor3us",                      gl_color3us                         },
{ "glColor3usv",                     gl_color3usv                        },
{ "glColor4b",                       gl_color4b                          },
{ "glColor4bv",                      gl_color4bv                         },
{ "glColor4d",                       gl_color4d                          },
{ "glColor4dv",                      gl_color4dv                         },
{ "glColor4f",                       gl_color4f                          },
{ "glColor4fv",                      gl_color4fv                         },
{ "glColor4i",                       gl_color4i                          },
{ "glColor4iv",                      gl_color4iv                         },
{ "glColor4s",                       gl_color4s                          },
{ "glColor4sv",                      gl_color4sv                         },
{ "glColor4ub",                      gl_color4ub                         },
{ "glColor4ubv",                     gl_color4ubv                        },
{ "glColor4ui",                      gl_color4ui                         },
{ "glColor4uiv",                     gl_color4uiv                        },
{ "glColor4us",                      gl_color4us                         },
{ "glColor4usv",                     gl_color4usv                        },
{ "glColorMask",                     gl_colormask                        },
{ "glColorMaterial",                 gl_colormaterial                    },
{ "glColorP3ui",                     gl_colorp3ui                        },
{ "glColorP3uiv",                    gl_colorp3uiv                       },
{ "glColorP4ui",                     gl_colorp4ui                        },
{ "glColorP4uiv",                    gl_colorp4uiv                       },
{ "glColorSubTable",                 gl_colorsubtable                    },
{ "glColorTable",                    gl_colortable                       },
{ "glColorTableParameterfv",         gl_colortableparameterfv            },
{ "glColorTableParameteriv",         gl_colortableparameteriv            },
{ "glCompileShader",                 gl_compileshader                    },
{ "glCompressedTexImage1D",          gl_compressedteximage1d             },
{ "glCompressedTexImage2D",          gl_compressedteximage2d             },
{ "glCompressedTexImage3D",          gl_compressedteximage3d             },
{ "glCompressedTexSubImage1D",       gl_compressedtexsubimage1d          },
{ "glCompressedTexSubImage2D",       gl_compressedtexsubimage2d          },
{ "glCompressedTexSubImage3D",       gl_compressedtexsubimage3d          },
{ "glConvolutionFilter1D",           gl_convolutionfilter1d              },
{ "glConvolutionFilter2D",           gl_convolutionfilter2d              },
{ "glConvolutionParameterf",         gl_convolutionparameterf            },
{ "glConvolutionParameterfv",        gl_convolutionparameterfv           },
{ "glConvolutionParameteri",         gl_convolutionparameteri            },
{ "glConvolutionParameteriv",        gl_convolutionparameteriv           },
{ "glCopyBufferSubData",             gl_copybuffersubdata                },
{ "glCopyColorSubTable",             gl_copycolorsubtable                },
{ "glCopyColorTable",                gl_copycolortable                   },
{ "glCopyConvolutionFilter1D",       gl_copyconvolutionfilter1d          },
{ "glCopyConvolutionFilter2D",       gl_copyconvolutionfilter2d          },
{ "glCopyImageSubData",              gl_copyimagesubdata                 },
{ "glCopyPixels",                    gl_copypixels                       },
{ "glCopyTexImage1D",                gl_copyteximage1d                   },
{ "glCopyTexImage2D",                gl_copyteximage2d                   },
{ "glCopyTexSubImage1D",             gl_copytexsubimage1d                },
{ "glCopyTexSubImage2D",             gl_copytexsubimage2d                },
{ "glCopyTexSubImage3D",             gl_copytexsubimage3d                },
{ "glCreateProgram",                 gl_createprogram                    },
{ "glCreateShader",                  gl_createshader                     },
{ "glCreateShaderProgramv",          gl_createshaderprogramv             },
{ "glCullFace",                      gl_cullface                         },
{ "glDebugMessageControl",           gl_debugmessagecontrol              },
{ "glDebugMessageControlARB",        gl_debugmessagecontrolarb           },
{ "glDebugMessageInsert",            gl_debugmessageinsert               },
{ "glDebugMessageInsertARB",         gl_debugmessageinsertarb            },
{ "glDeleteBuffers",                 gl_deletebuffers                    },
{ "glDeleteFramebuffers",            gl_deleteframebuffers               },
{ "glDeleteFramebuffersEXT",         gl_deleteframebuffersext            },
{ "glDeleteLists",                   gl_deletelists                      },
{ "glDeleteNamedStringARB",          gl_deletenamedstringarb             },
{ "glDeleteProgram",                 gl_deleteprogram                    },
{ "glDeleteProgramPipelines",        gl_deleteprogrampipelines           },
{ "glDeleteQueries",                 gl_deletequeries                    },
{ "glDeleteRenderbuffers",           gl_deleterenderbuffers              },
{ "glDeleteRenderbuffersEXT",        gl_deleterenderbuffersext           },
{ "glDeleteSamplers",                gl_deletesamplers                   },
{ "glDeleteShader",                  gl_deleteshader                     },
{ "glDeleteSync",                    gl_deletesync                       },
{ "glDeleteTextures",                gl_deletetextures                   },
{ "glDeleteTransformFeedbacks",      gl_deletetransformfeedbacks         },
{ "glDeleteVertexArrays",            gl_deletevertexarrays               },
{ "glDepthFunc",                     gl_depthfunc                        },
{ "glDepthMask",                     gl_depthmask                        },
{ "glDepthRange",                    gl_depthrange                       },
{ "glDepthRangeArrayv",              gl_depthrangearrayv                 },
{ "glDepthRangeIndexed",             gl_depthrangeindexed                },
{ "glDepthRangef",                   gl_depthrangef                      },
{ "glDetachShader",                  gl_detachshader                     },
{ "glDisable",                       gl_disable                          },
{ "glDisableClientState",            gl_disableclientstate               },
{ "glDisableVertexAttribArray",      gl_disablevertexattribarray         },
{ "glDispatchCompute",               gl_dispatchcompute                  },
{ "glDispatchComputeIndirect",       gl_dispatchcomputeindirect          },
{ "glDrawArrays",                    gl_drawarrays                       },
{ "glDrawArraysIndirect",            gl_drawarraysindirect               },
{ "glDrawArraysInstancedARB",        gl_drawarraysinstancedarb           },
{ "glDrawArraysInstancedBaseInstance",gl_drawarraysinstancedbaseinstance  },
{ "glDrawBuffer",                    gl_drawbuffer                       },
{ "glDrawBuffers",                   gl_drawbuffers                      },
{ "glDrawElementsBaseVertex",        gl_drawelementsbasevertex           },
{ "glDrawElementsIndirect",          gl_drawelementsindirect             },
{ "glDrawElementsInstancedARB",      gl_drawelementsinstancedarb         },
{ "glDrawElementsInstancedBaseInstance",gl_drawelementsinstancedbaseinstance },
{ "glDrawElementsInstancedBaseVertex",gl_drawelementsinstancedbasevertex  },
{ "glDrawElementsInstancedBaseVertexBaseInstance",gl_drawelementsinstancedbasevertexbaseinstance },
{ "glDrawPixels",                    gl_drawpixels                       },
{ "glDrawRangeElementsBaseVertex",   gl_drawrangeelementsbasevertex      },
{ "glDrawTransformFeedback",         gl_drawtransformfeedback            },
{ "glDrawTransformFeedbackInstanced",gl_drawtransformfeedbackinstanced   },
{ "glDrawTransformFeedbackStream",   gl_drawtransformfeedbackstream      },
{ "glDrawTransformFeedbackStreamInstanced",gl_drawtransformfeedbackstreaminstanced },
{ "glEdgeFlag",                      gl_edgeflag                         },
{ "glEdgeFlagPointer",               gl_edgeflagpointer                  },
{ "glEdgeFlagv",                     gl_edgeflagv                        },
{ "glEnable",                        gl_enable                           },
{ "glEnableClientState",             gl_enableclientstate                },
{ "glEnableVertexAttribArray",       gl_enablevertexattribarray          },
{ "glEnd",                           gl_end                              },
{ "glEndList",                       gl_endlist                          },
{ "glEndQuery",                      gl_endquery                         },
{ "glEndQueryIndexed",               gl_endqueryindexed                  },
{ "glEvalCoord1d",                   gl_evalcoord1d                      },
{ "glEvalCoord1dv",                  gl_evalcoord1dv                     },
{ "glEvalCoord1f",                   gl_evalcoord1f                      },
{ "glEvalCoord1fv",                  gl_evalcoord1fv                     },
{ "glEvalCoord2d",                   gl_evalcoord2d                      },
{ "glEvalCoord2dv",                  gl_evalcoord2dv                     },
{ "glEvalCoord2f",                   gl_evalcoord2f                      },
{ "glEvalCoord2fv",                  gl_evalcoord2fv                     },
{ "glEvalMesh1",                     gl_evalmesh1                        },
{ "glEvalMesh2",                     gl_evalmesh2                        },
{ "glEvalPoint1",                    gl_evalpoint1                       },
{ "glEvalPoint2",                    gl_evalpoint2                       },
{ "glFenceSync",                     gl_fencesync                        },
{ "glFinish",                        gl_finish                           },
{ "glFlush",                         gl_flush                            },
{ "glFlushMappedBufferRange",        gl_flushmappedbufferrange           },
{ "glFogCoordPointer",               gl_fogcoordpointer                  },
{ "glFogCoordd",                     gl_fogcoordd                        },
{ "glFogCoorddv",                    gl_fogcoorddv                       },
{ "glFogCoordf",                     gl_fogcoordf                        },
{ "glFogCoordfv",                    gl_fogcoordfv                       },
{ "glFogf",                          gl_fogf                             },
{ "glFogfv",                         gl_fogfv                            },
{ "glFogi",                          gl_fogi                             },
{ "glFogiv",                         gl_fogiv                            },
{ "glFramebufferParameteri",         gl_framebufferparameteri            },
{ "glFramebufferRenderbuffer",       gl_framebufferrenderbuffer          },
{ "glFramebufferRenderbufferEXT",    gl_framebufferrenderbufferext       },
{ "glFramebufferTexture1D",          gl_framebuffertexture1d             },
{ "glFramebufferTexture1DEXT",       gl_framebuffertexture1dext          },
{ "glFramebufferTexture2D",          gl_framebuffertexture2d             },
{ "glFramebufferTexture2DEXT",       gl_framebuffertexture2dext          },
{ "glFramebufferTexture3D",          gl_framebuffertexture3d             },
{ "glFramebufferTexture3DEXT",       gl_framebuffertexture3dext          },
{ "glFramebufferTextureARB",         gl_framebuffertexturearb            },
{ "glFramebufferTextureFaceARB",     gl_framebuffertexturefacearb        },
{ "glFramebufferTextureLayer",       gl_framebuffertexturelayer          },
{ "glFramebufferTextureLayerARB",    gl_framebuffertexturelayerarb       },
{ "glFrontFace",                     gl_frontface                        },
{ "glFrustum",                       gl_frustum                          },
{ "glGenBuffers",                    gl_genbuffers                       },
{ "glGenFramebuffers",               gl_genframebuffers                  },
{ "glGenFramebuffersEXT",            gl_genframebuffersext               },
{ "glGenLists",                      gl_genlists                         },
{ "glGenProgramPipelines",           gl_genprogrampipelines              },
{ "glGenQueries",                    gl_genqueries                       },
{ "glGenRenderbuffers",              gl_genrenderbuffers                 },
{ "glGenRenderbuffersEXT",           gl_genrenderbuffersext              },
{ "glGenSamplers",                   gl_gensamplers                      },
{ "glGenTextures",                   gl_gentextures                      },
{ "glGenTransformFeedbacks",         gl_gentransformfeedbacks            },
{ "glGenVertexArrays",               gl_genvertexarrays                  },
{ "glGenerateMipmap",                gl_generatemipmap                   },
{ "glGenerateMipmapEXT",             gl_generatemipmapext                },
{ "glGetActiveAtomicCounterBufferiv",gl_getactiveatomiccounterbufferiv   },
{ "glGetActiveAttrib",               gl_getactiveattrib                  },
{ "glGetActiveSubroutineName",       gl_getactivesubroutinename          },
{ "glGetActiveSubroutineUniformName",gl_getactivesubroutineuniformname   },
{ "glGetActiveSubroutineUniformiv",  gl_getactivesubroutineuniformiv     },
{ "glGetActiveUniform",              gl_getactiveuniform                 },
{ "glGetActiveUniformBlockName",     gl_getactiveuniformblockname        },
{ "glGetActiveUniformBlockiv",       gl_getactiveuniformblockiv          },
{ "glGetActiveUniformName",          gl_getactiveuniformname             },
{ "glGetActiveUniformsiv",           gl_getactiveuniformsiv              },
{ "glGetAttachedShaders",            gl_getattachedshaders               },
{ "glGetAttribLocation",             gl_getattriblocation                },
{ "glGetBooleanv",                   gl_getbooleanv                      },
{ "glGetBufferParameteriv",          gl_getbufferparameteriv             },
{ "glGetBufferSubData",              gl_getbuffersubdata                 },
{ "glGetClipPlane",                  gl_getclipplane                     },
{ "glGetColorTable",                 gl_getcolortable                    },
{ "glGetColorTableParameterfv",      gl_getcolortableparameterfv         },
{ "glGetColorTableParameteriv",      gl_getcolortableparameteriv         },
{ "glGetCompressedTexImage",         gl_getcompressedteximage            },
{ "glGetConvolutionFilter",          gl_getconvolutionfilter             },
{ "glGetConvolutionParameterfv",     gl_getconvolutionparameterfv        },
{ "glGetConvolutionParameteriv",     gl_getconvolutionparameteriv        },
{ "glGetDebugMessageLog",            gl_getdebugmessagelog               },
{ "glGetDebugMessageLogARB",         gl_getdebugmessagelogarb            },
{ "glGetDoublei_v",                  gl_getdoublei_v                     },
{ "glGetDoublev",                    gl_getdoublev                       },
{ "glGetError",                      gl_geterror                         },
{ "glGetFloati_v",                   gl_getfloati_v                      },
{ "glGetFloatv",                     gl_getfloatv                        },
{ "glGetFragDataIndex",              gl_getfragdataindex                 },
{ "glGetFramebufferAttachmentParameteriv",gl_getframebufferattachmentparameteriv },
{ "glGetFramebufferAttachmentParameterivEXT",gl_getframebufferattachmentparameterivext },
{ "glGetFramebufferParameteriv",     gl_getframebufferparameteriv        },
{ "glGetGraphicsResetStatusARB",     gl_getgraphicsresetstatusarb        },
{ "glGetHistogram",                  gl_gethistogram                     },
{ "glGetHistogramParameterfv",       gl_gethistogramparameterfv          },
{ "glGetHistogramParameteriv",       gl_gethistogramparameteriv          },
{ "glGetInteger64v",                 gl_getinteger64v                    },
{ "glGetIntegerv",                   gl_getintegerv                      },
{ "glGetInternalformati64v",         gl_getinternalformati64v            },
{ "glGetInternalformativ",           gl_getinternalformativ              },
{ "glGetLightfv",                    gl_getlightfv                       },
{ "glGetLightiv",                    gl_getlightiv                       },
{ "glGetMapdv",                      gl_getmapdv                         },
{ "glGetMapfv",                      gl_getmapfv                         },
{ "glGetMapiv",                      gl_getmapiv                         },
{ "glGetMaterialfv",                 gl_getmaterialfv                    },
{ "glGetMaterialiv",                 gl_getmaterialiv                    },
{ "glGetMinmax",                     gl_getminmax                        },
{ "glGetMinmaxParameterfv",          gl_getminmaxparameterfv             },
{ "glGetMinmaxParameteriv",          gl_getminmaxparameteriv             },
{ "glGetMultisamplefv",              gl_getmultisamplefv                 },
{ "glGetNamedFramebufferParameterivEXT",gl_getnamedframebufferparameterivext },
{ "glGetNamedStringARB",             gl_getnamedstringarb                },
{ "glGetNamedStringivARB",           gl_getnamedstringivarb              },
{ "glGetObjectLabel",                gl_getobjectlabel                   },
{ "glGetObjectPtrLabel",             gl_getobjectptrlabel                },
{ "glGetPixelMapfv",                 gl_getpixelmapfv                    },
{ "glGetPixelMapuiv",                gl_getpixelmapuiv                   },
{ "glGetPixelMapusv",                gl_getpixelmapusv                   },
{ "glGetPolygonStipple",             gl_getpolygonstipple                },
{ "glGetProgramBinary",              gl_getprogrambinary                 },
{ "glGetProgramInfoLog",             gl_getprograminfolog                },
{ "glGetProgramInterfaceiv",         gl_getprograminterfaceiv            },
{ "glGetProgramPipelineInfoLog",     gl_getprogrampipelineinfolog        },
{ "glGetProgramPipelineiv",          gl_getprogrampipelineiv             },
{ "glGetProgramResourceIndex",       gl_getprogramresourceindex          },
{ "glGetProgramResourceLocation",    gl_getprogramresourcelocation       },
{ "glGetProgramResourceLocationIndex",gl_getprogramresourcelocationindex  },
{ "glGetProgramResourceName",        gl_getprogramresourcename           },
{ "glGetProgramResourceiv",          gl_getprogramresourceiv             },
{ "glGetProgramStageiv",             gl_getprogramstageiv                },
{ "glGetProgramiv",                  gl_getprogramiv                     },
{ "glGetQueryIndexediv",             gl_getqueryindexediv                },
{ "glGetQueryObjecti64v",            gl_getqueryobjecti64v               },
{ "glGetQueryObjectiv",              gl_getqueryobjectiv                 },
{ "glGetQueryObjectui64v",           gl_getqueryobjectui64v              },
{ "glGetQueryObjectuiv",             gl_getqueryobjectuiv                },
{ "glGetQueryiv",                    gl_getqueryiv                       },
{ "glGetRenderbufferParameteriv",    gl_getrenderbufferparameteriv       },
{ "glGetRenderbufferParameterivEXT", gl_getrenderbufferparameterivext    },
{ "glGetSamplerParameterIiv",        gl_getsamplerparameteriiv           },
{ "glGetSamplerParameterIuiv",       gl_getsamplerparameteriuiv          },
{ "glGetSamplerParameterfv",         gl_getsamplerparameterfv            },
{ "glGetSamplerParameteriv",         gl_getsamplerparameteriv            },
{ "glGetSeparableFilter",            gl_getseparablefilter               },
{ "glGetShaderInfoLog",              gl_getshaderinfolog                 },
{ "glGetShaderPrecisionFormat",      gl_getshaderprecisionformat         },
{ "glGetShaderSource",               gl_getshadersource                  },
{ "glGetShaderiv",                   gl_getshaderiv                      },
{ "glGetSubroutineIndex",            gl_getsubroutineindex               },
{ "glGetSubroutineUniformLocation",  gl_getsubroutineuniformlocation     },
{ "glGetSynciv",                     gl_getsynciv                        },
{ "glGetTexEnvfv",                   gl_gettexenvfv                      },
{ "glGetTexEnviv",                   gl_gettexenviv                      },
{ "glGetTexGendv",                   gl_gettexgendv                      },
{ "glGetTexGenfv",                   gl_gettexgenfv                      },
{ "glGetTexGeniv",                   gl_gettexgeniv                      },
{ "glGetTexImage",                   gl_getteximage                      },
{ "glGetTexLevelParameterfv",        gl_gettexlevelparameterfv           },
{ "glGetTexLevelParameteriv",        gl_gettexlevelparameteriv           },
{ "glGetTexParameterfv",             gl_gettexparameterfv                },
{ "glGetTexParameteriv",             gl_gettexparameteriv                },
{ "glGetUniformBlockIndex",          gl_getuniformblockindex             },
{ "glGetUniformIndices",             gl_getuniformindices                },
{ "glGetUniformLocation",            gl_getuniformlocation               },
{ "glGetUniformSubroutineuiv",       gl_getuniformsubroutineuiv          },
{ "glGetUniformdv",                  gl_getuniformdv                     },
{ "glGetUniformfv",                  gl_getuniformfv                     },
{ "glGetUniformiv",                  gl_getuniformiv                     },
{ "glGetVertexAttribLdv",            gl_getvertexattribldv               },
{ "glGetVertexAttribdv",             gl_getvertexattribdv                },
{ "glGetVertexAttribfv",             gl_getvertexattribfv                },
{ "glGetVertexAttribiv",             gl_getvertexattribiv                },
{ "glGetnColorTableARB",             gl_getncolortablearb                },
{ "glGetnCompressedTexImageARB",     gl_getncompressedteximagearb        },
{ "glGetnConvolutionFilterARB",      gl_getnconvolutionfilterarb         },
{ "glGetnHistogramARB",              gl_getnhistogramarb                 },
{ "glGetnMapdvARB",                  gl_getnmapdvarb                     },
{ "glGetnMapfvARB",                  gl_getnmapfvarb                     },
{ "glGetnMapivARB",                  gl_getnmapivarb                     },
{ "glGetnMinmaxARB",                 gl_getnminmaxarb                    },
{ "glGetnPixelMapfvARB",             gl_getnpixelmapfvarb                },
{ "glGetnPixelMapuivARB",            gl_getnpixelmapuivarb               },
{ "glGetnPixelMapusvARB",            gl_getnpixelmapusvarb               },
{ "glGetnPolygonStippleARB",         gl_getnpolygonstipplearb            },
{ "glGetnSeparableFilterARB",        gl_getnseparablefilterarb           },
{ "glGetnTexImageARB",               gl_getnteximagearb                  },
{ "glGetnUniformdvARB",              gl_getnuniformdvarb                 },
{ "glGetnUniformfvARB",              gl_getnuniformfvarb                 },
{ "glGetnUniformivARB",              gl_getnuniformivarb                 },
{ "glGetnUniformuivARB",             gl_getnuniformuivarb                },
{ "glHint",                          gl_hint                             },
{ "glHistogram",                     gl_histogram                        },
{ "glIndexMask",                     gl_indexmask                        },
{ "glIndexPointer",                  gl_indexpointer                     },
{ "glIndexd",                        gl_indexd                           },
{ "glIndexdv",                       gl_indexdv                          },
{ "glIndexf",                        gl_indexf                           },
{ "glIndexfv",                       gl_indexfv                          },
{ "glIndexi",                        gl_indexi                           },
{ "glIndexiv",                       gl_indexiv                          },
{ "glIndexs",                        gl_indexs                           },
{ "glIndexsv",                       gl_indexsv                          },
{ "glIndexub",                       gl_indexub                          },
{ "glIndexubv",                      gl_indexubv                         },
{ "glInitNames",                     gl_initnames                        },
{ "glInterleavedArrays",             gl_interleavedarrays                },
{ "glInvalidateBufferData",          gl_invalidatebufferdata             },
{ "glInvalidateBufferSubData",       gl_invalidatebuffersubdata          },
{ "glInvalidateFramebuffer",         gl_invalidateframebuffer            },
{ "glInvalidateSubFramebuffer",      gl_invalidatesubframebuffer         },
{ "glInvalidateTexImage",            gl_invalidateteximage               },
{ "glInvalidateTexSubImage",         gl_invalidatetexsubimage            },
{ "glIsBuffer",                      gl_isbuffer                         },
{ "glIsEnabled",                     gl_isenabled                        },
{ "glIsFramebuffer",                 gl_isframebuffer                    },
{ "glIsFramebufferEXT",              gl_isframebufferext                 },
{ "glIsList",                        gl_islist                           },
{ "glIsNamedStringARB",              gl_isnamedstringarb                 },
{ "glIsProgram",                     gl_isprogram                        },
{ "glIsProgramPipeline",             gl_isprogrampipeline                },
{ "glIsQuery",                       gl_isquery                          },
{ "glIsRenderbuffer",                gl_isrenderbuffer                   },
{ "glIsRenderbufferEXT",             gl_isrenderbufferext                },
{ "glIsSampler",                     gl_issampler                        },
{ "glIsShader",                      gl_isshader                         },
{ "glIsSync",                        gl_issync                           },
{ "glIsTexture",                     gl_istexture                        },
{ "glIsTransformFeedback",           gl_istransformfeedback              },
{ "glIsVertexArray",                 gl_isvertexarray                    },
{ "glLightModelf",                   gl_lightmodelf                      },
{ "glLightModelfv",                  gl_lightmodelfv                     },
{ "glLightModeli",                   gl_lightmodeli                      },
{ "glLightModeliv",                  gl_lightmodeliv                     },
{ "glLightf",                        gl_lightf                           },
{ "glLightfv",                       gl_lightfv                          },
{ "glLighti",                        gl_lighti                           },
{ "glLightiv",                       gl_lightiv                          },
{ "glLineStipple",                   gl_linestipple                      },
{ "glLineWidth",                     gl_linewidth                        },
{ "glLinkProgram",                   gl_linkprogram                      },
{ "glListBase",                      gl_listbase                         },
{ "glLoadIdentity",                  gl_loadidentity                     },
{ "glLoadMatrixd",                   gl_loadmatrixd                      },
{ "glLoadMatrixf",                   gl_loadmatrixf                      },
{ "glLoadName",                      gl_loadname                         },
{ "glLoadTransposeMatrixd",          gl_loadtransposematrixd             },
{ "glLoadTransposeMatrixf",          gl_loadtransposematrixf             },
{ "glLogicOp",                       gl_logicop                          },
{ "glMap1d",                         gl_map1d                            },
{ "glMap1f",                         gl_map1f                            },
{ "glMap2d",                         gl_map2d                            },
{ "glMap2f",                         gl_map2f                            },
{ "glMapBuffer",                     gl_mapbuffer                        },
{ "glMapBufferRange",                gl_mapbufferrange                   },
{ "glMapGrid1d",                     gl_mapgrid1d                        },
{ "glMapGrid1f",                     gl_mapgrid1f                        },
{ "glMapGrid2d",                     gl_mapgrid2d                        },
{ "glMapGrid2f",                     gl_mapgrid2f                        },
{ "glMaterialf",                     gl_materialf                        },
{ "glMaterialfv",                    gl_materialfv                       },
{ "glMateriali",                     gl_materiali                        },
{ "glMaterialiv",                    gl_materialiv                       },
{ "glMatrixMode",                    gl_matrixmode                       },
{ "glMemoryBarrier",                 gl_memorybarrier                    },
{ "glMinSampleShadingARB",           gl_minsampleshadingarb              },
{ "glMinmax",                        gl_minmax                           },
{ "glMultMatrixd",                   gl_multmatrixd                      },
{ "glMultMatrixf",                   gl_multmatrixf                      },
{ "glMultTransposeMatrixd",          gl_multtransposematrixd             },
{ "glMultTransposeMatrixf",          gl_multtransposematrixf             },
{ "glMultiDrawArrays",               gl_multidrawarrays                  },
{ "glMultiDrawArraysIndirect",       gl_multidrawarraysindirect          },
{ "glMultiDrawElementsBaseVertex",   gl_multidrawelementsbasevertex      },
{ "glMultiDrawElementsIndirect",     gl_multidrawelementsindirect        },
{ "glMultiTexCoord1d",               gl_multitexcoord1d                  },
{ "glMultiTexCoord1dv",              gl_multitexcoord1dv                 },
{ "glMultiTexCoord1f",               gl_multitexcoord1f                  },
{ "glMultiTexCoord1fv",              gl_multitexcoord1fv                 },
{ "glMultiTexCoord1i",               gl_multitexcoord1i                  },
{ "glMultiTexCoord1iv",              gl_multitexcoord1iv                 },
{ "glMultiTexCoord1s",               gl_multitexcoord1s                  },
{ "glMultiTexCoord1sv",              gl_multitexcoord1sv                 },
{ "glMultiTexCoord2d",               gl_multitexcoord2d                  },
{ "glMultiTexCoord2dv",              gl_multitexcoord2dv                 },
{ "glMultiTexCoord2f",               gl_multitexcoord2f                  },
{ "glMultiTexCoord2fv",              gl_multitexcoord2fv                 },
{ "glMultiTexCoord2i",               gl_multitexcoord2i                  },
{ "glMultiTexCoord2iv",              gl_multitexcoord2iv                 },
{ "glMultiTexCoord2s",               gl_multitexcoord2s                  },
{ "glMultiTexCoord2sv",              gl_multitexcoord2sv                 },
{ "glMultiTexCoord3d",               gl_multitexcoord3d                  },
{ "glMultiTexCoord3dv",              gl_multitexcoord3dv                 },
{ "glMultiTexCoord3f",               gl_multitexcoord3f                  },
{ "glMultiTexCoord3fv",              gl_multitexcoord3fv                 },
{ "glMultiTexCoord3i",               gl_multitexcoord3i                  },
{ "glMultiTexCoord3iv",              gl_multitexcoord3iv                 },
{ "glMultiTexCoord3s",               gl_multitexcoord3s                  },
{ "glMultiTexCoord3sv",              gl_multitexcoord3sv                 },
{ "glMultiTexCoord4d",               gl_multitexcoord4d                  },
{ "glMultiTexCoord4dv",              gl_multitexcoord4dv                 },
{ "glMultiTexCoord4f",               gl_multitexcoord4f                  },
{ "glMultiTexCoord4fv",              gl_multitexcoord4fv                 },
{ "glMultiTexCoord4i",               gl_multitexcoord4i                  },
{ "glMultiTexCoord4iv",              gl_multitexcoord4iv                 },
{ "glMultiTexCoord4s",               gl_multitexcoord4s                  },
{ "glMultiTexCoord4sv",              gl_multitexcoord4sv                 },
{ "glMultiTexCoordP1ui",             gl_multitexcoordp1ui                },
{ "glMultiTexCoordP1uiv",            gl_multitexcoordp1uiv               },
{ "glMultiTexCoordP2ui",             gl_multitexcoordp2ui                },
{ "glMultiTexCoordP2uiv",            gl_multitexcoordp2uiv               },
{ "glMultiTexCoordP3ui",             gl_multitexcoordp3ui                },
{ "glMultiTexCoordP3uiv",            gl_multitexcoordp3uiv               },
{ "glMultiTexCoordP4ui",             gl_multitexcoordp4ui                },
{ "glMultiTexCoordP4uiv",            gl_multitexcoordp4uiv               },
{ "glNamedFramebufferParameteriEXT", gl_namedframebufferparameteriext    },
{ "glNamedStringARB",                gl_namedstringarb                   },
{ "glNewList",                       gl_newlist                          },
{ "glNormal3b",                      gl_normal3b                         },
{ "glNormal3bv",                     gl_normal3bv                        },
{ "glNormal3d",                      gl_normal3d                         },
{ "glNormal3dv",                     gl_normal3dv                        },
{ "glNormal3f",                      gl_normal3f                         },
{ "glNormal3fv",                     gl_normal3fv                        },
{ "glNormal3i",                      gl_normal3i                         },
{ "glNormal3iv",                     gl_normal3iv                        },
{ "glNormal3s",                      gl_normal3s                         },
{ "glNormal3sv",                     gl_normal3sv                        },
{ "glNormalP3ui",                    gl_normalp3ui                       },
{ "glNormalP3uiv",                   gl_normalp3uiv                      },
{ "glObjectLabel",                   gl_objectlabel                      },
{ "glObjectPtrLabel",                gl_objectptrlabel                   },
{ "glOrtho",                         gl_ortho                            },
{ "glPassThrough",                   gl_passthrough                      },
{ "glPatchParameterfv",              gl_patchparameterfv                 },
{ "glPatchParameteri",               gl_patchparameteri                  },
{ "glPauseTransformFeedback",        gl_pausetransformfeedback           },
{ "glPixelMapfv",                    gl_pixelmapfv                       },
{ "glPixelMapuiv",                   gl_pixelmapuiv                      },
{ "glPixelMapusv",                   gl_pixelmapusv                      },
{ "glPixelStoref",                   gl_pixelstoref                      },
{ "glPixelStorei",                   gl_pixelstorei                      },
{ "glPixelTransferf",                gl_pixeltransferf                   },
{ "glPixelTransferi",                gl_pixeltransferi                   },
{ "glPixelZoom",                     gl_pixelzoom                        },
{ "glPointParameterf",               gl_pointparameterf                  },
{ "glPointParameterfv",              gl_pointparameterfv                 },
{ "glPointParameteri",               gl_pointparameteri                  },
{ "glPointParameteriv",              gl_pointparameteriv                 },
{ "glPointSize",                     gl_pointsize                        },
{ "glPolygonMode",                   gl_polygonmode                      },
{ "glPolygonOffset",                 gl_polygonoffset                    },
{ "glPolygonStipple",                gl_polygonstipple                   },
{ "glPopAttrib",                     gl_popattrib                        },
{ "glPopClientAttrib",               gl_popclientattrib                  },
{ "glPopMatrix",                     gl_popmatrix                        },
{ "glPopName",                       gl_popname                          },
{ "glPrioritizeTextures",            gl_prioritizetextures               },
{ "glProgramBinary",                 gl_programbinary                    },
{ "glProgramParameteri",             gl_programparameteri                },
{ "glProgramParameteriARB",          gl_programparameteriarb             },
{ "glProgramUniform1d",              gl_programuniform1d                 },
{ "glProgramUniform1dv",             gl_programuniform1dv                },
{ "glProgramUniform1f",              gl_programuniform1f                 },
{ "glProgramUniform1fv",             gl_programuniform1fv                },
{ "glProgramUniform1i",              gl_programuniform1i                 },
{ "glProgramUniform1iv",             gl_programuniform1iv                },
{ "glProgramUniform1ui",             gl_programuniform1ui                },
{ "glProgramUniform1uiv",            gl_programuniform1uiv               },
{ "glProgramUniform2d",              gl_programuniform2d                 },
{ "glProgramUniform2dv",             gl_programuniform2dv                },
{ "glProgramUniform2f",              gl_programuniform2f                 },
{ "glProgramUniform2fv",             gl_programuniform2fv                },
{ "glProgramUniform2i",              gl_programuniform2i                 },
{ "glProgramUniform2iv",             gl_programuniform2iv                },
{ "glProgramUniform2ui",             gl_programuniform2ui                },
{ "glProgramUniform2uiv",            gl_programuniform2uiv               },
{ "glProgramUniform3d",              gl_programuniform3d                 },
{ "glProgramUniform3dv",             gl_programuniform3dv                },
{ "glProgramUniform3f",              gl_programuniform3f                 },
{ "glProgramUniform3fv",             gl_programuniform3fv                },
{ "glProgramUniform3i",              gl_programuniform3i                 },
{ "glProgramUniform3iv",             gl_programuniform3iv                },
{ "glProgramUniform3ui",             gl_programuniform3ui                },
{ "glProgramUniform3uiv",            gl_programuniform3uiv               },
{ "glProgramUniform4d",              gl_programuniform4d                 },
{ "glProgramUniform4dv",             gl_programuniform4dv                },
{ "glProgramUniform4f",              gl_programuniform4f                 },
{ "glProgramUniform4fv",             gl_programuniform4fv                },
{ "glProgramUniform4i",              gl_programuniform4i                 },
{ "glProgramUniform4iv",             gl_programuniform4iv                },
{ "glProgramUniform4ui",             gl_programuniform4ui                },
{ "glProgramUniform4uiv",            gl_programuniform4uiv               },
{ "glProgramUniformMatrix2dv",       gl_programuniformmatrix2dv          },
{ "glProgramUniformMatrix2fv",       gl_programuniformmatrix2fv          },
{ "glProgramUniformMatrix2x3dv",     gl_programuniformmatrix2x3dv        },
{ "glProgramUniformMatrix2x3fv",     gl_programuniformmatrix2x3fv        },
{ "glProgramUniformMatrix2x4dv",     gl_programuniformmatrix2x4dv        },
{ "glProgramUniformMatrix2x4fv",     gl_programuniformmatrix2x4fv        },
{ "glProgramUniformMatrix3dv",       gl_programuniformmatrix3dv          },
{ "glProgramUniformMatrix3fv",       gl_programuniformmatrix3fv          },
{ "glProgramUniformMatrix3x2dv",     gl_programuniformmatrix3x2dv        },
{ "glProgramUniformMatrix3x2fv",     gl_programuniformmatrix3x2fv        },
{ "glProgramUniformMatrix3x4dv",     gl_programuniformmatrix3x4dv        },
{ "glProgramUniformMatrix3x4fv",     gl_programuniformmatrix3x4fv        },
{ "glProgramUniformMatrix4dv",       gl_programuniformmatrix4dv          },
{ "glProgramUniformMatrix4fv",       gl_programuniformmatrix4fv          },
{ "glProgramUniformMatrix4x2dv",     gl_programuniformmatrix4x2dv        },
{ "glProgramUniformMatrix4x2fv",     gl_programuniformmatrix4x2fv        },
{ "glProgramUniformMatrix4x3dv",     gl_programuniformmatrix4x3dv        },
{ "glProgramUniformMatrix4x3fv",     gl_programuniformmatrix4x3fv        },
{ "glProvokingVertex",               gl_provokingvertex                  },
{ "glPushAttrib",                    gl_pushattrib                       },
{ "glPushClientAttrib",              gl_pushclientattrib                 },
{ "glPushDebugGroup",                gl_pushdebuggroup                   },
{ "glPushMatrix",                    gl_pushmatrix                       },
{ "glPushName",                      gl_pushname                         },
{ "glQueryCounter",                  gl_querycounter                     },
{ "glRasterPos2d",                   gl_rasterpos2d                      },
{ "glRasterPos2dv",                  gl_rasterpos2dv                     },
{ "glRasterPos2f",                   gl_rasterpos2f                      },
{ "glRasterPos2fv",                  gl_rasterpos2fv                     },
{ "glRasterPos2i",                   gl_rasterpos2i                      },
{ "glRasterPos2iv",                  gl_rasterpos2iv                     },
{ "glRasterPos2s",                   gl_rasterpos2s                      },
{ "glRasterPos2sv",                  gl_rasterpos2sv                     },
{ "glRasterPos3d",                   gl_rasterpos3d                      },
{ "glRasterPos3dv",                  gl_rasterpos3dv                     },
{ "glRasterPos3f",                   gl_rasterpos3f                      },
{ "glRasterPos3fv",                  gl_rasterpos3fv                     },
{ "glRasterPos3i",                   gl_rasterpos3i                      },
{ "glRasterPos3iv",                  gl_rasterpos3iv                     },
{ "glRasterPos3s",                   gl_rasterpos3s                      },
{ "glRasterPos3sv",                  gl_rasterpos3sv                     },
{ "glRasterPos4d",                   gl_rasterpos4d                      },
{ "glRasterPos4dv",                  gl_rasterpos4dv                     },
{ "glRasterPos4f",                   gl_rasterpos4f                      },
{ "glRasterPos4fv",                  gl_rasterpos4fv                     },
{ "glRasterPos4i",                   gl_rasterpos4i                      },
{ "glRasterPos4iv",                  gl_rasterpos4iv                     },
{ "glRasterPos4s",                   gl_rasterpos4s                      },
{ "glRasterPos4sv",                  gl_rasterpos4sv                     },
{ "glReadBuffer",                    gl_readbuffer                       },
{ "glReadnPixelsARB",                gl_readnpixelsarb                   },
{ "glRectd",                         gl_rectd                            },
{ "glRectdv",                        gl_rectdv                           },
{ "glRectf",                         gl_rectf                            },
{ "glRectfv",                        gl_rectfv                           },
{ "glRecti",                         gl_recti                            },
{ "glRectiv",                        gl_rectiv                           },
{ "glRects",                         gl_rects                            },
{ "glRectsv",                        gl_rectsv                           },
{ "glReleaseShaderCompiler",         gl_releaseshadercompiler            },
{ "glRenderMode",                    gl_rendermode                       },
{ "glRenderbufferStorage",           gl_renderbufferstorage              },
{ "glRenderbufferStorageEXT",        gl_renderbufferstorageext           },
{ "glRenderbufferStorageMultisample",gl_renderbufferstoragemultisample   },
{ "glResetHistogram",                gl_resethistogram                   },
{ "glResetMinmax",                   gl_resetminmax                      },
{ "glResumeTransformFeedback",       gl_resumetransformfeedback          },
{ "glRotated",                       gl_rotated                          },
{ "glRotatef",                       gl_rotatef                          },
{ "glSampleCoverage",                gl_samplecoverage                   },
{ "glSampleMaskSGIS",                gl_samplemasksgis                   },
{ "glSampleMaski",                   gl_samplemaski                      },
{ "glSamplePatternSGIS",             gl_samplepatternsgis                },
{ "glSamplerParameterIiv",           gl_samplerparameteriiv              },
{ "glSamplerParameterIuiv",          gl_samplerparameteriuiv             },
{ "glSamplerParameterf",             gl_samplerparameterf                },
{ "glSamplerParameterfv",            gl_samplerparameterfv               },
{ "glSamplerParameteri",             gl_samplerparameteri                },
{ "glSamplerParameteriv",            gl_samplerparameteriv               },
{ "glScaled",                        gl_scaled                           },
{ "glScalef",                        gl_scalef                           },
{ "glScissor",                       gl_scissor                          },
{ "glScissorArrayv",                 gl_scissorarrayv                    },
{ "glScissorIndexed",                gl_scissorindexed                   },
{ "glScissorIndexedv",               gl_scissorindexedv                  },
{ "glSecondaryColor3b",              gl_secondarycolor3b                 },
{ "glSecondaryColor3bv",             gl_secondarycolor3bv                },
{ "glSecondaryColor3d",              gl_secondarycolor3d                 },
{ "glSecondaryColor3dv",             gl_secondarycolor3dv                },
{ "glSecondaryColor3f",              gl_secondarycolor3f                 },
{ "glSecondaryColor3fv",             gl_secondarycolor3fv                },
{ "glSecondaryColor3i",              gl_secondarycolor3i                 },
{ "glSecondaryColor3iv",             gl_secondarycolor3iv                },
{ "glSecondaryColor3s",              gl_secondarycolor3s                 },
{ "glSecondaryColor3sv",             gl_secondarycolor3sv                },
{ "glSecondaryColor3ub",             gl_secondarycolor3ub                },
{ "glSecondaryColor3ubv",            gl_secondarycolor3ubv               },
{ "glSecondaryColor3ui",             gl_secondarycolor3ui                },
{ "glSecondaryColor3uiv",            gl_secondarycolor3uiv               },
{ "glSecondaryColor3us",             gl_secondarycolor3us                },
{ "glSecondaryColor3usv",            gl_secondarycolor3usv               },
{ "glSecondaryColorP3ui",            gl_secondarycolorp3ui               },
{ "glSecondaryColorP3uiv",           gl_secondarycolorp3uiv              },
{ "glSecondaryColorPointer",         gl_secondarycolorpointer            },
{ "glSeparableFilter2D",             gl_separablefilter2d                },
{ "glShadeModel",                    gl_shademodel                       },
{ "glShaderBinary",                  gl_shaderbinary                     },
{ "glShaderStorageBlockBinding",     gl_shaderstorageblockbinding        },
{ "glStencilFunc",                   gl_stencilfunc                      },
{ "glStencilFuncSeparate",           gl_stencilfuncseparate              },
{ "glStencilMask",                   gl_stencilmask                      },
{ "glStencilMaskSeparate",           gl_stencilmaskseparate              },
{ "glStencilOp",                     gl_stencilop                        },
{ "glStencilOpSeparate",             gl_stencilopseparate                },
{ "glTexBufferARB",                  gl_texbufferarb                     },
{ "glTexBufferRange",                gl_texbufferrange                   },
{ "glTexCoord1d",                    gl_texcoord1d                       },
{ "glTexCoord1dv",                   gl_texcoord1dv                      },
{ "glTexCoord1f",                    gl_texcoord1f                       },
{ "glTexCoord1fv",                   gl_texcoord1fv                      },
{ "glTexCoord1i",                    gl_texcoord1i                       },
{ "glTexCoord1iv",                   gl_texcoord1iv                      },
{ "glTexCoord1s",                    gl_texcoord1s                       },
{ "glTexCoord1sv",                   gl_texcoord1sv                      },
{ "glTexCoord2d",                    gl_texcoord2d                       },
{ "glTexCoord2dv",                   gl_texcoord2dv                      },
{ "glTexCoord2f",                    gl_texcoord2f                       },
{ "glTexCoord2fv",                   gl_texcoord2fv                      },
{ "glTexCoord2i",                    gl_texcoord2i                       },
{ "glTexCoord2iv",                   gl_texcoord2iv                      },
{ "glTexCoord2s",                    gl_texcoord2s                       },
{ "glTexCoord2sv",                   gl_texcoord2sv                      },
{ "glTexCoord3d",                    gl_texcoord3d                       },
{ "glTexCoord3dv",                   gl_texcoord3dv                      },
{ "glTexCoord3f",                    gl_texcoord3f                       },
{ "glTexCoord3fv",                   gl_texcoord3fv                      },
{ "glTexCoord3i",                    gl_texcoord3i                       },
{ "glTexCoord3iv",                   gl_texcoord3iv                      },
{ "glTexCoord3s",                    gl_texcoord3s                       },
{ "glTexCoord3sv",                   gl_texcoord3sv                      },
{ "glTexCoord4d",                    gl_texcoord4d                       },
{ "glTexCoord4dv",                   gl_texcoord4dv                      },
{ "glTexCoord4f",                    gl_texcoord4f                       },
{ "glTexCoord4fv",                   gl_texcoord4fv                      },
{ "glTexCoord4i",                    gl_texcoord4i                       },
{ "glTexCoord4iv",                   gl_texcoord4iv                      },
{ "glTexCoord4s",                    gl_texcoord4s                       },
{ "glTexCoord4sv",                   gl_texcoord4sv                      },
{ "glTexCoordP1ui",                  gl_texcoordp1ui                     },
{ "glTexCoordP1uiv",                 gl_texcoordp1uiv                    },
{ "glTexCoordP2ui",                  gl_texcoordp2ui                     },
{ "glTexCoordP2uiv",                 gl_texcoordp2uiv                    },
{ "glTexCoordP3ui",                  gl_texcoordp3ui                     },
{ "glTexCoordP3uiv",                 gl_texcoordp3uiv                    },
{ "glTexCoordP4ui",                  gl_texcoordp4ui                     },
{ "glTexCoordP4uiv",                 gl_texcoordp4uiv                    },
{ "glTexEnvf",                       gl_texenvf                          },
{ "glTexEnvfv",                      gl_texenvfv                         },
{ "glTexEnvi",                       gl_texenvi                          },
{ "glTexEnviv",                      gl_texenviv                         },
{ "glTexGend",                       gl_texgend                          },
{ "glTexGendv",                      gl_texgendv                         },
{ "glTexGenf",                       gl_texgenf                          },
{ "glTexGenfv",                      gl_texgenfv                         },
{ "glTexGeni",                       gl_texgeni                          },
{ "glTexGeniv",                      gl_texgeniv                         },
{ "glTexImage1D",                    gl_teximage1d                       },
{ "glTexImage2DMultisample",         gl_teximage2dmultisample            },
{ "glTexImage3D",                    gl_teximage3d                       },
{ "glTexImage3DMultisample",         gl_teximage3dmultisample            },
{ "glTexParameterf",                 gl_texparameterf                    },
{ "glTexParameterfv",                gl_texparameterfv                   },
{ "glTexParameteri",                 gl_texparameteri                    },
{ "glTexParameteriv",                gl_texparameteriv                   },
{ "glTexStorage1D",                  gl_texstorage1d                     },
{ "glTexStorage2D",                  gl_texstorage2d                     },
{ "glTexStorage2DMultisample",       gl_texstorage2dmultisample          },
{ "glTexStorage3D",                  gl_texstorage3d                     },
{ "glTexStorage3DMultisample",       gl_texstorage3dmultisample          },
{ "glTexSubImage1D",                 gl_texsubimage1d                    },
{ "glTexSubImage2D",                 gl_texsubimage2d                    },
{ "glTexSubImage3D",                 gl_texsubimage3d                    },
{ "glTextureBufferRangeEXT",         gl_texturebufferrangeext            },
{ "glTextureStorage1DEXT",           gl_texturestorage1dext              },
{ "glTextureStorage2DEXT",           gl_texturestorage2dext              },
{ "glTextureStorage2DMultisampleEXT",gl_texturestorage2dmultisampleext   },
{ "glTextureStorage3DEXT",           gl_texturestorage3dext              },
{ "glTextureStorage3DMultisampleEXT",gl_texturestorage3dmultisampleext   },
{ "glTextureView",                   gl_textureview                      },
{ "glTranslated",                    gl_translated                       },
{ "glTranslatef",                    gl_translatef                       },
{ "glUniform1d",                     gl_uniform1d                        },
{ "glUniform1dv",                    gl_uniform1dv                       },
{ "glUniform1f",                     gl_uniform1f                        },
{ "glUniform1fv",                    gl_uniform1fv                       },
{ "glUniform1i",                     gl_uniform1i                        },
{ "glUniform1iv",                    gl_uniform1iv                       },
{ "glUniform2d",                     gl_uniform2d                        },
{ "glUniform2dv",                    gl_uniform2dv                       },
{ "glUniform2f",                     gl_uniform2f                        },
{ "glUniform2fv",                    gl_uniform2fv                       },
{ "glUniform2i",                     gl_uniform2i                        },
{ "glUniform2iv",                    gl_uniform2iv                       },
{ "glUniform3d",                     gl_uniform3d                        },
{ "glUniform3dv",                    gl_uniform3dv                       },
{ "glUniform3f",                     gl_uniform3f                        },
{ "glUniform3fv",                    gl_uniform3fv                       },
{ "glUniform3i",                     gl_uniform3i                        },
{ "glUniform3iv",                    gl_uniform3iv                       },
{ "glUniform4d",                     gl_uniform4d                        },
{ "glUniform4dv",                    gl_uniform4dv                       },
{ "glUniform4f",                     gl_uniform4f                        },
{ "glUniform4fv",                    gl_uniform4fv                       },
{ "glUniform4i",                     gl_uniform4i                        },
{ "glUniform4iv",                    gl_uniform4iv                       },
{ "glUniformBlockBinding",           gl_uniformblockbinding              },
{ "glUniformMatrix2dv",              gl_uniformmatrix2dv                 },
{ "glUniformMatrix2fv",              gl_uniformmatrix2fv                 },
{ "glUniformMatrix2x3dv",            gl_uniformmatrix2x3dv               },
{ "glUniformMatrix2x3fv",            gl_uniformmatrix2x3fv               },
{ "glUniformMatrix2x4dv",            gl_uniformmatrix2x4dv               },
{ "glUniformMatrix2x4fv",            gl_uniformmatrix2x4fv               },
{ "glUniformMatrix3dv",              gl_uniformmatrix3dv                 },
{ "glUniformMatrix3fv",              gl_uniformmatrix3fv                 },
{ "glUniformMatrix3x2dv",            gl_uniformmatrix3x2dv               },
{ "glUniformMatrix3x2fv",            gl_uniformmatrix3x2fv               },
{ "glUniformMatrix3x4dv",            gl_uniformmatrix3x4dv               },
{ "glUniformMatrix3x4fv",            gl_uniformmatrix3x4fv               },
{ "glUniformMatrix4dv",              gl_uniformmatrix4dv                 },
{ "glUniformMatrix4fv",              gl_uniformmatrix4fv                 },
{ "glUniformMatrix4x2dv",            gl_uniformmatrix4x2dv               },
{ "glUniformMatrix4x2fv",            gl_uniformmatrix4x2fv               },
{ "glUniformMatrix4x3dv",            gl_uniformmatrix4x3dv               },
{ "glUniformMatrix4x3fv",            gl_uniformmatrix4x3fv               },
{ "glUniformSubroutinesuiv",         gl_uniformsubroutinesuiv            },
{ "glUnmapBuffer",                   gl_unmapbuffer                      },
{ "glUseProgram",                    gl_useprogram                       },
{ "glUseProgramStages",              gl_useprogramstages                 },
{ "glValidateProgram",               gl_validateprogram                  },
{ "glValidateProgramPipeline",       gl_validateprogrampipeline          },
{ "glVertex2d",                      gl_vertex2d                         },
{ "glVertex2dv",                     gl_vertex2dv                        },
{ "glVertex2f",                      gl_vertex2f                         },
{ "glVertex2fv",                     gl_vertex2fv                        },
{ "glVertex2i",                      gl_vertex2i                         },
{ "glVertex2iv",                     gl_vertex2iv                        },
{ "glVertex2s",                      gl_vertex2s                         },
{ "glVertex2sv",                     gl_vertex2sv                        },
{ "glVertex3d",                      gl_vertex3d                         },
{ "glVertex3dv",                     gl_vertex3dv                        },
{ "glVertex3f",                      gl_vertex3f                         },
{ "glVertex3fv",                     gl_vertex3fv                        },
{ "glVertex3i",                      gl_vertex3i                         },
{ "glVertex3iv",                     gl_vertex3iv                        },
{ "glVertex3s",                      gl_vertex3s                         },
{ "glVertex3sv",                     gl_vertex3sv                        },
{ "glVertex4d",                      gl_vertex4d                         },
{ "glVertex4dv",                     gl_vertex4dv                        },
{ "glVertex4f",                      gl_vertex4f                         },
{ "glVertex4fv",                     gl_vertex4fv                        },
{ "glVertex4i",                      gl_vertex4i                         },
{ "glVertex4iv",                     gl_vertex4iv                        },
{ "glVertex4s",                      gl_vertex4s                         },
{ "glVertex4sv",                     gl_vertex4sv                        },
{ "glVertexAttrib1d",                gl_vertexattrib1d                   },
{ "glVertexAttrib1dv",               gl_vertexattrib1dv                  },
{ "glVertexAttrib1f",                gl_vertexattrib1f                   },
{ "glVertexAttrib1fv",               gl_vertexattrib1fv                  },
{ "glVertexAttrib1s",                gl_vertexattrib1s                   },
{ "glVertexAttrib1sv",               gl_vertexattrib1sv                  },
{ "glVertexAttrib2d",                gl_vertexattrib2d                   },
{ "glVertexAttrib2dv",               gl_vertexattrib2dv                  },
{ "glVertexAttrib2f",                gl_vertexattrib2f                   },
{ "glVertexAttrib2fv",               gl_vertexattrib2fv                  },
{ "glVertexAttrib2s",                gl_vertexattrib2s                   },
{ "glVertexAttrib2sv",               gl_vertexattrib2sv                  },
{ "glVertexAttrib3d",                gl_vertexattrib3d                   },
{ "glVertexAttrib3dv",               gl_vertexattrib3dv                  },
{ "glVertexAttrib3f",                gl_vertexattrib3f                   },
{ "glVertexAttrib3fv",               gl_vertexattrib3fv                  },
{ "glVertexAttrib3s",                gl_vertexattrib3s                   },
{ "glVertexAttrib3sv",               gl_vertexattrib3sv                  },
{ "glVertexAttrib4Nbv",              gl_vertexattrib4nbv                 },
{ "glVertexAttrib4Niv",              gl_vertexattrib4niv                 },
{ "glVertexAttrib4Nsv",              gl_vertexattrib4nsv                 },
{ "glVertexAttrib4Nub",              gl_vertexattrib4nub                 },
{ "glVertexAttrib4Nubv",             gl_vertexattrib4nubv                },
{ "glVertexAttrib4Nuiv",             gl_vertexattrib4nuiv                },
{ "glVertexAttrib4Nusv",             gl_vertexattrib4nusv                },
{ "glVertexAttrib4bv",               gl_vertexattrib4bv                  },
{ "glVertexAttrib4d",                gl_vertexattrib4d                   },
{ "glVertexAttrib4dv",               gl_vertexattrib4dv                  },
{ "glVertexAttrib4f",                gl_vertexattrib4f                   },
{ "glVertexAttrib4fv",               gl_vertexattrib4fv                  },
{ "glVertexAttrib4iv",               gl_vertexattrib4iv                  },
{ "glVertexAttrib4s",                gl_vertexattrib4s                   },
{ "glVertexAttrib4sv",               gl_vertexattrib4sv                  },
{ "glVertexAttrib4ubv",              gl_vertexattrib4ubv                 },
{ "glVertexAttrib4uiv",              gl_vertexattrib4uiv                 },
{ "glVertexAttrib4usv",              gl_vertexattrib4usv                 },
{ "glVertexAttribBinding",           gl_vertexattribbinding              },
{ "glVertexAttribDivisorARB",        gl_vertexattribdivisorarb           },
{ "glVertexAttribFormat",            gl_vertexattribformat               },
{ "glVertexAttribIFormat",           gl_vertexattribiformat              },
{ "glVertexAttribL1d",               gl_vertexattribl1d                  },
{ "glVertexAttribL1dv",              gl_vertexattribl1dv                 },
{ "glVertexAttribL2d",               gl_vertexattribl2d                  },
{ "glVertexAttribL2dv",              gl_vertexattribl2dv                 },
{ "glVertexAttribL3d",               gl_vertexattribl3d                  },
{ "glVertexAttribL3dv",              gl_vertexattribl3dv                 },
{ "glVertexAttribL4d",               gl_vertexattribl4d                  },
{ "glVertexAttribL4dv",              gl_vertexattribl4dv                 },
{ "glVertexAttribLFormat",           gl_vertexattriblformat              },
{ "glVertexAttribLPointer",          gl_vertexattriblpointer             },
{ "glVertexAttribP1ui",              gl_vertexattribp1ui                 },
{ "glVertexAttribP1uiv",             gl_vertexattribp1uiv                },
{ "glVertexAttribP2ui",              gl_vertexattribp2ui                 },
{ "glVertexAttribP2uiv",             gl_vertexattribp2uiv                },
{ "glVertexAttribP3ui",              gl_vertexattribp3ui                 },
{ "glVertexAttribP3uiv",             gl_vertexattribp3uiv                },
{ "glVertexAttribP4ui",              gl_vertexattribp4ui                 },
{ "glVertexAttribP4uiv",             gl_vertexattribp4uiv                },
{ "glVertexBindingDivisor",          gl_vertexbindingdivisor             },
{ "glVertexP2ui",                    gl_vertexp2ui                       },
{ "glVertexP2uiv",                   gl_vertexp2uiv                      },
{ "glVertexP3ui",                    gl_vertexp3ui                       },
{ "glVertexP3uiv",                   gl_vertexp3uiv                      },
{ "glVertexP4ui",                    gl_vertexp4ui                       },
{ "glVertexP4uiv",                   gl_vertexp4uiv                      },
{ "glViewport",                      gl_viewport                         },
{ "glViewportArrayv",                gl_viewportarrayv                   },
{ "glViewportIndexedf",              gl_viewportindexedf                 },
{ "glViewportIndexedfv",             gl_viewportindexedfv                },
{ "glWaitSync",                      gl_waitsync                         },
{ "glWindowPos2d",                   gl_windowpos2d                      },
{ "glWindowPos2dv",                  gl_windowpos2dv                     },
{ "glWindowPos2f",                   gl_windowpos2f                      },
{ "glWindowPos2fv",                  gl_windowpos2fv                     },
{ "glWindowPos2i",                   gl_windowpos2i                      },
{ "glWindowPos2iv",                  gl_windowpos2iv                     },
{ "glWindowPos2s",                   gl_windowpos2s                      },
{ "glWindowPos2sv",                  gl_windowpos2sv                     },
{ "glWindowPos3d",                   gl_windowpos3d                      },
{ "glWindowPos3dv",                  gl_windowpos3dv                     },
{ "glWindowPos3f",                   gl_windowpos3f                      },
{ "glWindowPos3fv",                  gl_windowpos3fv                     },
{ "glWindowPos3i",                   gl_windowpos3i                      },
{ "glWindowPos3iv",                  gl_windowpos3iv                     },
{ "glWindowPos3s",                   gl_windowpos3s                      },
{ "glWindowPos3sv",                  gl_windowpos3sv                     },
{ "gluBeginCurve",                   glu_begincurve                      },
{ "gluBeginSurface",                 glu_beginsurface                    },
{ "gluBeginTrim",                    glu_begintrim                       },
{ "gluBuild1DMipmapLevels",          glu_build1dmipmaplevels             },
{ "gluBuild1DMipmaps",               glu_build1dmipmaps                  },
{ "gluBuild2DMipmapLevels",          glu_build2dmipmaplevels             },
{ "gluBuild2DMipmaps",               glu_build2dmipmaps                  },
{ "gluBuild3DMipmapLevels",          glu_build3dmipmaplevels             },
{ "gluBuild3DMipmaps",               glu_build3dmipmaps                  },
{ "gluCheckExtension",               glu_checkextension                  },
{ "gluCylinder",                     glu_cylinder                        },
{ "gluDeleteNurbsRenderer",          glu_deletenurbsrenderer             },
{ "gluDeleteQuadric",                glu_deletequadric                   },
{ "gluDisk",                         glu_disk                            },
{ "gluEndCurve",                     glu_endcurve                        },
{ "gluEndSurface",                   glu_endsurface                      },
{ "gluEndTrim",                      glu_endtrim                         },
{ "gluGetNurbsProperty",             glu_getnurbsproperty                },
{ "gluGetTessProperty",              glu_gettessproperty                 },
{ "gluLoadSamplingMatrices",         glu_loadsamplingmatrices            },
{ "gluLookAt",                       glu_lookat                          },
{ "gluNewNurbsRenderer",             glu_newnurbsrenderer                },
{ "gluNewQuadric",                   glu_newquadric                      },
{ "gluNurbsCurve",                   glu_nurbscurve                      },
{ "gluNurbsProperty",                glu_nurbsproperty                   },
{ "gluNurbsSurface",                 glu_nurbssurface                    },
{ "gluOrtho2D",                      glu_ortho2d                         },
{ "gluPartialDisk",                  glu_partialdisk                     },
{ "gluPerspective",                  glu_perspective                     },
{ "gluPickMatrix",                   glu_pickmatrix                      },
{ "gluProject",                      glu_project                         },
{ "gluPwlCurve",                     glu_pwlcurve                        },
{ "gluQuadricDrawStyle",             glu_quadricdrawstyle                },
{ "gluQuadricNormals",               glu_quadricnormals                  },
{ "gluQuadricOrientation",           glu_quadricorientation              },
{ "gluQuadricTexture",               glu_quadrictexture                  },
{ "gluScaleImage",                   glu_scaleimage                      },
{ "gluSphere",                       glu_sphere                          },
{ "gluUnProject",                    glu_unproject                       },
{ "gluUnProject4",                   glu_unproject4                      },
{ "glutSolidCone",                   glut_solidcone                      },
{ "glutSolidCube",                   glut_solidcube                      },
{ "glutSolidDodecahedron",           glut_soliddodecahedron              },
{ "glutSolidIcosahedron",            glut_solidicosahedron               },
{ "glutSolidOctahedron",             glut_solidoctahedron                },
{ "glutSolidSphere",                 glut_solidsphere                    },
{ "glutSolidTeapot",                 glut_solidteapot                    },
{ "glutSolidTetrahedron",            glut_solidtetrahedron               },
{ "glutSolidTorus",                  glut_solidtorus                     },
{ "glutWireCone",                    glut_wirecone                       },
{ "glutWireCube",                    glut_wirecube                       },
{ "glutWireDodecahedron",            glut_wiredodecahedron               },
{ "glutWireIcosahedron",             glut_wireicosahedron                },
{ "glutWireOctahedron",              glut_wireoctahedron                 },
{ "glutWireSphere",                  glut_wiresphere                     },
{ "glutWireTeapot",                  glut_wireteapot                     },
{ "glutWireTetrahedron",             glut_wiretetrahedron                },
{ "glutWireTorus",                   glut_wiretorus                      }};
