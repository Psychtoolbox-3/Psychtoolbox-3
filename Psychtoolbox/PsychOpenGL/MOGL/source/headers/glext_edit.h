/* NOTE: This is a drastically stripped-down version of the official glext.h
 * file from http://www.opengl.org/registry/api/glext.h
 *
 * Only the subset of function prototype definitions is extracted from the
 * original file, which we actually use -- Plus a bit of context, e.g., all
 * preprocessor directives (#define, #ifdef etc.) and all 'typedef' lines
 * and comments are ignored by our parser. We just leave this dead lines
 * here to simplify diff's to the original file etc.
 *
 * Adding new definitions of OpenGL extensions for use with PTB+MOGL and
 * generating code would be as easy as copying a larger subset of function
 * definitions from the original file - and potentially extending the parser
 * to handle more complex function definitions, e.g., we currently don't handle
 * functions whose pointer arguments contain double-indirections.
 *
 * Updated by Mario Kleiner at 1st April 2012.
 * Updated by Mario Kleiner at 28. August 2012 up to GL_ARB_texture_storage_multisample
 *
 */

/*
** Copyright (c) 2007-2012 The Khronos Group Inc.
** 
** Permission is hereby granted, free of charge, to any person obtaining a
** copy of this software and/or associated documentation files (the
** "Materials"), to deal in the Materials without restriction, including
** without limitation the rights to use, copy, modify, merge, publish,
** distribute, sublicense, and/or sell copies of the Materials, and to
** permit persons to whom the Materials are furnished to do so, subject to
** the following conditions:
** 
** The above copyright notice and this permission notice shall be included
** in all copies or substantial portions of the Materials.
** 
** THE MATERIALS ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
** MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
** TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
** MATERIALS OR THE USE OR OTHER DEALINGS IN THE MATERIALS.
*/

/* Header file version number, required by OpenGL ABI for Linux */
/* glext.h last updated $Date: 2012-03-28 22:34:23 -0700 (Wed, 28 Mar 2012) $ */
/* glext.h last updated $Date: 2012-08-13 16:18:01 -0700 (Mon, 13 Aug 2012) $ */

/* Current version at http://www.opengl.org/registry/ */
#define GL_GLEXT_VERSION 78
/* Function declaration macros - to move into glplatform.h */

#ifndef GL_EXT_framebuffer_object
#define GL_EXT_framebuffer_object 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI GLboolean APIENTRY glIsRenderbufferEXT (GLuint renderbuffer);
GLAPI void APIENTRY glBindRenderbufferEXT (GLenum target, GLuint renderbuffer);
GLAPI void APIENTRY glDeleteRenderbuffersEXT (GLsizei n, const GLuint *renderbuffers);
GLAPI void APIENTRY glGenRenderbuffersEXT (GLsizei n, GLuint *renderbuffers);
GLAPI void APIENTRY glRenderbufferStorageEXT (GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
GLAPI void APIENTRY glGetRenderbufferParameterivEXT (GLenum target, GLenum pname, GLint *params);
GLAPI GLboolean APIENTRY glIsFramebufferEXT (GLuint framebuffer);
GLAPI void APIENTRY glBindFramebufferEXT (GLenum target, GLuint framebuffer);
GLAPI void APIENTRY glDeleteFramebuffersEXT (GLsizei n, const GLuint *framebuffers);
GLAPI void APIENTRY glGenFramebuffersEXT (GLsizei n, GLuint *framebuffers);
GLAPI GLenum APIENTRY glCheckFramebufferStatusEXT (GLenum target);
GLAPI void APIENTRY glFramebufferTexture1DEXT (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
GLAPI void APIENTRY glFramebufferTexture2DEXT (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
GLAPI void APIENTRY glFramebufferTexture3DEXT (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset);
GLAPI void APIENTRY glFramebufferRenderbufferEXT (GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);
GLAPI void APIENTRY glGetFramebufferAttachmentParameterivEXT (GLenum target, GLenum attachment, GLenum pname, GLint *params);
GLAPI void APIENTRY glGenerateMipmapEXT (GLenum target);
#endif /* GL_GLEXT_PROTOTYPES */
typedef GLboolean (APIENTRYP PFNGLISRENDERBUFFEREXTPROC) (GLuint renderbuffer);
typedef void (APIENTRYP PFNGLBINDRENDERBUFFEREXTPROC) (GLenum target, GLuint renderbuffer);
typedef void (APIENTRYP PFNGLDELETERENDERBUFFERSEXTPROC) (GLsizei n, const GLuint *renderbuffers);
typedef void (APIENTRYP PFNGLGENRENDERBUFFERSEXTPROC) (GLsizei n, GLuint *renderbuffers);
typedef void (APIENTRYP PFNGLRENDERBUFFERSTORAGEEXTPROC) (GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
typedef void (APIENTRYP PFNGLGETRENDERBUFFERPARAMETERIVEXTPROC) (GLenum target, GLenum pname, GLint *params);
typedef GLboolean (APIENTRYP PFNGLISFRAMEBUFFEREXTPROC) (GLuint framebuffer);
typedef void (APIENTRYP PFNGLBINDFRAMEBUFFEREXTPROC) (GLenum target, GLuint framebuffer);
typedef void (APIENTRYP PFNGLDELETEFRAMEBUFFERSEXTPROC) (GLsizei n, const GLuint *framebuffers);
typedef void (APIENTRYP PFNGLGENFRAMEBUFFERSEXTPROC) (GLsizei n, GLuint *framebuffers);
typedef GLenum (APIENTRYP PFNGLCHECKFRAMEBUFFERSTATUSEXTPROC) (GLenum target);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE1DEXTPROC) (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE2DEXTPROC) (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE3DEXTPROC) (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset);
typedef void (APIENTRYP PFNGLFRAMEBUFFERRENDERBUFFEREXTPROC) (GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);
typedef void (APIENTRYP PFNGLGETFRAMEBUFFERATTACHMENTPARAMETERIVEXTPROC) (GLenum target, GLenum attachment, GLenum pname, GLint *params);
typedef void (APIENTRYP PFNGLGENERATEMIPMAPEXTPROC) (GLenum target);
#endif

#ifndef GL_SGIS_multisample
#define GL_SGIS_multisample 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glSampleMaskSGIS (GLclampf value, GLboolean invert);
GLAPI void APIENTRY glSamplePatternSGIS (GLenum pattern);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLSAMPLEMASKSGISPROC) (GLclampf value, GLboolean invert);
typedef void (APIENTRYP PFNGLSAMPLEPATTERNSGISPROC) (GLenum pattern);
#endif

/* NEW SINCE 28-August-2012 with glext.h update to GL_GLEXT_VERSION 84: */

#ifndef GL_ARB_draw_instanced
#define GL_ARB_draw_instanced 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDrawArraysInstancedARB (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
GLAPI void APIENTRY glDrawElementsInstancedARB (GLenum mode, GLsizei count, GLenum type, const GLvoid *indices, GLsizei primcount);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
typedef void (APIENTRYP PFNGLDRAWELEMENTSINSTANCEDARBPROC) (GLenum mode, GLsizei count, GLenum type, const GLvoid *indices, GLsizei primcount);
#endif

#ifndef GL_ARB_framebuffer_object
#define GL_ARB_framebuffer_object 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI GLboolean APIENTRY glIsRenderbuffer (GLuint renderbuffer);
GLAPI void APIENTRY glBindRenderbuffer (GLenum target, GLuint renderbuffer);
GLAPI void APIENTRY glDeleteRenderbuffers (GLsizei n, const GLuint *renderbuffers);
GLAPI void APIENTRY glGenRenderbuffers (GLsizei n, GLuint *renderbuffers);
GLAPI void APIENTRY glRenderbufferStorage (GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
GLAPI void APIENTRY glGetRenderbufferParameteriv (GLenum target, GLenum pname, GLint *params);
GLAPI GLboolean APIENTRY glIsFramebuffer (GLuint framebuffer);
GLAPI void APIENTRY glBindFramebuffer (GLenum target, GLuint framebuffer);
GLAPI void APIENTRY glDeleteFramebuffers (GLsizei n, const GLuint *framebuffers);
GLAPI void APIENTRY glGenFramebuffers (GLsizei n, GLuint *framebuffers);
GLAPI GLenum APIENTRY glCheckFramebufferStatus (GLenum target);
GLAPI void APIENTRY glFramebufferTexture1D (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
GLAPI void APIENTRY glFramebufferTexture2D (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
GLAPI void APIENTRY glFramebufferTexture3D (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset);
GLAPI void APIENTRY glFramebufferRenderbuffer (GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);
GLAPI void APIENTRY glGetFramebufferAttachmentParameteriv (GLenum target, GLenum attachment, GLenum pname, GLint *params);
GLAPI void APIENTRY glGenerateMipmap (GLenum target);
GLAPI void APIENTRY glBlitFramebuffer (GLint srcX0, GLint srcY0, GLint srcX1, GLint srcY1, GLint dstX0, GLint dstY0, GLint dstX1, GLint dstY1, GLbitfield mask, GLenum filter);
GLAPI void APIENTRY glRenderbufferStorageMultisample (GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height);
GLAPI void APIENTRY glFramebufferTextureLayer (GLenum target, GLenum attachment, GLuint texture, GLint level, GLint layer);
#endif /* GL_GLEXT_PROTOTYPES */
typedef GLboolean (APIENTRYP PFNGLISRENDERBUFFERPROC) (GLuint renderbuffer);
typedef void (APIENTRYP PFNGLBINDRENDERBUFFERPROC) (GLenum target, GLuint renderbuffer);
typedef void (APIENTRYP PFNGLDELETERENDERBUFFERSPROC) (GLsizei n, const GLuint *renderbuffers);
typedef void (APIENTRYP PFNGLGENRENDERBUFFERSPROC) (GLsizei n, GLuint *renderbuffers);
typedef void (APIENTRYP PFNGLRENDERBUFFERSTORAGEPROC) (GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
typedef void (APIENTRYP PFNGLGETRENDERBUFFERPARAMETERIVPROC) (GLenum target, GLenum pname, GLint *params);
typedef GLboolean (APIENTRYP PFNGLISFRAMEBUFFERPROC) (GLuint framebuffer);
typedef void (APIENTRYP PFNGLBINDFRAMEBUFFERPROC) (GLenum target, GLuint framebuffer);
typedef void (APIENTRYP PFNGLDELETEFRAMEBUFFERSPROC) (GLsizei n, const GLuint *framebuffers);
typedef void (APIENTRYP PFNGLGENFRAMEBUFFERSPROC) (GLsizei n, GLuint *framebuffers);
typedef GLenum (APIENTRYP PFNGLCHECKFRAMEBUFFERSTATUSPROC) (GLenum target);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE1DPROC) (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE2DPROC) (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURE3DPROC) (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset);
typedef void (APIENTRYP PFNGLFRAMEBUFFERRENDERBUFFERPROC) (GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);
typedef void (APIENTRYP PFNGLGETFRAMEBUFFERATTACHMENTPARAMETERIVPROC) (GLenum target, GLenum attachment, GLenum pname, GLint *params);
typedef void (APIENTRYP PFNGLGENERATEMIPMAPPROC) (GLenum target);
typedef void (APIENTRYP PFNGLBLITFRAMEBUFFERPROC) (GLint srcX0, GLint srcY0, GLint srcX1, GLint srcY1, GLint dstX0, GLint dstY0, GLint dstX1, GLint dstY1, GLbitfield mask, GLenum filter);
typedef void (APIENTRYP PFNGLRENDERBUFFERSTORAGEMULTISAMPLEPROC) (GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURELAYERPROC) (GLenum target, GLenum attachment, GLuint texture, GLint level, GLint layer);
#endif

#ifndef GL_ARB_framebuffer_sRGB
#define GL_ARB_framebuffer_sRGB 1
#endif

#ifndef GL_ARB_geometry_shader4
#define GL_ARB_geometry_shader4 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glProgramParameteriARB (GLuint program, GLenum pname, GLint value);
GLAPI void APIENTRY glFramebufferTextureARB (GLenum target, GLenum attachment, GLuint texture, GLint level);
GLAPI void APIENTRY glFramebufferTextureLayerARB (GLenum target, GLenum attachment, GLuint texture, GLint level, GLint layer);
GLAPI void APIENTRY glFramebufferTextureFaceARB (GLenum target, GLenum attachment, GLuint texture, GLint level, GLenum face);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLPROGRAMPARAMETERIARBPROC) (GLuint program, GLenum pname, GLint value);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTUREARBPROC) (GLenum target, GLenum attachment, GLuint texture, GLint level);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTURELAYERARBPROC) (GLenum target, GLenum attachment, GLuint texture, GLint level, GLint layer);
typedef void (APIENTRYP PFNGLFRAMEBUFFERTEXTUREFACEARBPROC) (GLenum target, GLenum attachment, GLuint texture, GLint level, GLenum face);
#endif

#ifndef GL_ARB_half_float_vertex
#define GL_ARB_half_float_vertex 1
#endif

#ifndef GL_ARB_instanced_arrays
#define GL_ARB_instanced_arrays 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glVertexAttribDivisorARB (GLuint index, GLuint divisor);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
#endif

#ifndef GL_ARB_map_buffer_range
#define GL_ARB_map_buffer_range 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI GLvoid* APIENTRY glMapBufferRange (GLenum target, GLintptr offset, GLsizeiptr length, GLbitfield access);
GLAPI void APIENTRY glFlushMappedBufferRange (GLenum target, GLintptr offset, GLsizeiptr length);
#endif /* GL_GLEXT_PROTOTYPES */
typedef GLvoid* (APIENTRYP PFNGLMAPBUFFERRANGEPROC) (GLenum target, GLintptr offset, GLsizeiptr length, GLbitfield access);
typedef void (APIENTRYP PFNGLFLUSHMAPPEDBUFFERRANGEPROC) (GLenum target, GLintptr offset, GLsizeiptr length);
#endif

#ifndef GL_ARB_texture_buffer_object
#define GL_ARB_texture_buffer_object 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glTexBufferARB (GLenum target, GLenum internalformat, GLuint buffer);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLTEXBUFFERARBPROC) (GLenum target, GLenum internalformat, GLuint buffer);
#endif

#ifndef GL_ARB_texture_compression_rgtc
#define GL_ARB_texture_compression_rgtc 1
#endif

#ifndef GL_ARB_texture_rg
#define GL_ARB_texture_rg 1
#endif

#ifndef GL_ARB_vertex_array_object
#define GL_ARB_vertex_array_object 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glBindVertexArray (GLuint array);
GLAPI void APIENTRY glDeleteVertexArrays (GLsizei n, const GLuint *arrays);
GLAPI void APIENTRY glGenVertexArrays (GLsizei n, GLuint *arrays);
GLAPI GLboolean APIENTRY glIsVertexArray (GLuint array);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLBINDVERTEXARRAYPROC) (GLuint array);
typedef void (APIENTRYP PFNGLDELETEVERTEXARRAYSPROC) (GLsizei n, const GLuint *arrays);
typedef void (APIENTRYP PFNGLGENVERTEXARRAYSPROC) (GLsizei n, GLuint *arrays);
typedef GLboolean (APIENTRYP PFNGLISVERTEXARRAYPROC) (GLuint array);
#endif

#ifndef GL_ARB_uniform_buffer_object
#define GL_ARB_uniform_buffer_object 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glGetUniformIndices (GLuint program, GLsizei uniformCount, const GLchar* const *uniformNames, GLuint *uniformIndices);
GLAPI void APIENTRY glGetActiveUniformsiv (GLuint program, GLsizei uniformCount, const GLuint *uniformIndices, GLenum pname, GLint *params);
GLAPI void APIENTRY glGetActiveUniformName (GLuint program, GLuint uniformIndex, GLsizei bufSize, GLsizei *length, GLchar *uniformName);
GLAPI GLuint APIENTRY glGetUniformBlockIndex (GLuint program, const GLchar *uniformBlockName);
GLAPI void APIENTRY glGetActiveUniformBlockiv (GLuint program, GLuint uniformBlockIndex, GLenum pname, GLint *params);
GLAPI void APIENTRY glGetActiveUniformBlockName (GLuint program, GLuint uniformBlockIndex, GLsizei bufSize, GLsizei *length, GLchar *uniformBlockName);
GLAPI void APIENTRY glUniformBlockBinding (GLuint program, GLuint uniformBlockIndex, GLuint uniformBlockBinding);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLGETUNIFORMINDICESPROC) (GLuint program, GLsizei uniformCount, const GLchar* const *uniformNames, GLuint *uniformIndices);
typedef void (APIENTRYP PFNGLGETACTIVEUNIFORMSIVPROC) (GLuint program, GLsizei uniformCount, const GLuint *uniformIndices, GLenum pname, GLint *params);
typedef void (APIENTRYP PFNGLGETACTIVEUNIFORMNAMEPROC) (GLuint program, GLuint uniformIndex, GLsizei bufSize, GLsizei *length, GLchar *uniformName);
typedef GLuint (APIENTRYP PFNGLGETUNIFORMBLOCKINDEXPROC) (GLuint program, const GLchar *uniformBlockName);
typedef void (APIENTRYP PFNGLGETACTIVEUNIFORMBLOCKIVPROC) (GLuint program, GLuint uniformBlockIndex, GLenum pname, GLint *params);
typedef void (APIENTRYP PFNGLGETACTIVEUNIFORMBLOCKNAMEPROC) (GLuint program, GLuint uniformBlockIndex, GLsizei bufSize, GLsizei *length, GLchar *uniformBlockName);
typedef void (APIENTRYP PFNGLUNIFORMBLOCKBINDINGPROC) (GLuint program, GLuint uniformBlockIndex, GLuint uniformBlockBinding);
#endif

#ifndef GL_ARB_compatibility
#define GL_ARB_compatibility 1
#endif

#ifndef GL_ARB_copy_buffer
#define GL_ARB_copy_buffer 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glCopyBufferSubData (GLenum readTarget, GLenum writeTarget, GLintptr readOffset, GLintptr writeOffset, GLsizeiptr size);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLCOPYBUFFERSUBDATAPROC) (GLenum readTarget, GLenum writeTarget, GLintptr readOffset, GLintptr writeOffset, GLsizeiptr size);
#endif

#ifndef GL_ARB_shader_texture_lod
#define GL_ARB_shader_texture_lod 1
#endif

#ifndef GL_ARB_depth_clamp
#define GL_ARB_depth_clamp 1
#endif

#ifndef GL_ARB_draw_elements_base_vertex
#define GL_ARB_draw_elements_base_vertex 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDrawElementsBaseVertex (GLenum mode, GLsizei count, GLenum type, const GLvoid *indices, GLint basevertex);
GLAPI void APIENTRY glDrawRangeElementsBaseVertex (GLenum mode, GLuint start, GLuint end, GLsizei count, GLenum type, const GLvoid *indices, GLint basevertex);
GLAPI void APIENTRY glDrawElementsInstancedBaseVertex (GLenum mode, GLsizei count, GLenum type, const GLvoid *indices, GLsizei instancecount, GLint basevertex);
GLAPI void APIENTRY glMultiDrawElementsBaseVertex (GLenum mode, const GLsizei *count, GLenum type, const GLvoid* const *indices, GLsizei drawcount, const GLint *basevertex);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDRAWELEMENTSBASEVERTEXPROC) (GLenum mode, GLsizei count, GLenum type, const GLvoid *indices, GLint basevertex);
typedef void (APIENTRYP PFNGLDRAWRANGEELEMENTSBASEVERTEXPROC) (GLenum mode, GLuint start, GLuint end, GLsizei count, GLenum type, const GLvoid *indices, GLint basevertex);
typedef void (APIENTRYP PFNGLDRAWELEMENTSINSTANCEDBASEVERTEXPROC) (GLenum mode, GLsizei count, GLenum type, const GLvoid *indices, GLsizei instancecount, GLint basevertex);
typedef void (APIENTRYP PFNGLMULTIDRAWELEMENTSBASEVERTEXPROC) (GLenum mode, const GLsizei *count, GLenum type, const GLvoid* const *indices, GLsizei drawcount, const GLint *basevertex);
#endif

#ifndef GL_ARB_fragment_coord_conventions
#define GL_ARB_fragment_coord_conventions 1
#endif

#ifndef GL_ARB_provoking_vertex
#define GL_ARB_provoking_vertex 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glProvokingVertex (GLenum mode);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLPROVOKINGVERTEXPROC) (GLenum mode);
#endif

#ifndef GL_ARB_seamless_cube_map
#define GL_ARB_seamless_cube_map 1
#endif

#ifndef GL_ARB_sync
#define GL_ARB_sync 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI GLsync APIENTRY glFenceSync (GLenum condition, GLbitfield flags);
GLAPI GLboolean APIENTRY glIsSync (GLsync sync);
GLAPI void APIENTRY glDeleteSync (GLsync sync);
GLAPI GLenum APIENTRY glClientWaitSync (GLsync sync, GLbitfield flags, GLuint64 timeout);
GLAPI void APIENTRY glWaitSync (GLsync sync, GLbitfield flags, GLuint64 timeout);
GLAPI void APIENTRY glGetInteger64v (GLenum pname, GLint64 *params);
GLAPI void APIENTRY glGetSynciv (GLsync sync, GLenum pname, GLsizei bufSize, GLsizei *length, GLint *values);
#endif /* GL_GLEXT_PROTOTYPES */
typedef GLsync (APIENTRYP PFNGLFENCESYNCPROC) (GLenum condition, GLbitfield flags);
typedef GLboolean (APIENTRYP PFNGLISSYNCPROC) (GLsync sync);
typedef void (APIENTRYP PFNGLDELETESYNCPROC) (GLsync sync);
typedef GLenum (APIENTRYP PFNGLCLIENTWAITSYNCPROC) (GLsync sync, GLbitfield flags, GLuint64 timeout);
typedef void (APIENTRYP PFNGLWAITSYNCPROC) (GLsync sync, GLbitfield flags, GLuint64 timeout);
typedef void (APIENTRYP PFNGLGETINTEGER64VPROC) (GLenum pname, GLint64 *params);
typedef void (APIENTRYP PFNGLGETSYNCIVPROC) (GLsync sync, GLenum pname, GLsizei bufSize, GLsizei *length, GLint *values);
#endif

#ifndef GL_ARB_texture_multisample
#define GL_ARB_texture_multisample 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glTexImage2DMultisample (GLenum target, GLsizei samples, GLint internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations);
GLAPI void APIENTRY glTexImage3DMultisample (GLenum target, GLsizei samples, GLint internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations);
GLAPI void APIENTRY glGetMultisamplefv (GLenum pname, GLuint index, GLfloat *val);
GLAPI void APIENTRY glSampleMaski (GLuint index, GLbitfield mask);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLTEXIMAGE2DMULTISAMPLEPROC) (GLenum target, GLsizei samples, GLint internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations);
typedef void (APIENTRYP PFNGLTEXIMAGE3DMULTISAMPLEPROC) (GLenum target, GLsizei samples, GLint internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations);
typedef void (APIENTRYP PFNGLGETMULTISAMPLEFVPROC) (GLenum pname, GLuint index, GLfloat *val);
typedef void (APIENTRYP PFNGLSAMPLEMASKIPROC) (GLuint index, GLbitfield mask);
#endif

#ifndef GL_ARB_vertex_array_bgra
#define GL_ARB_vertex_array_bgra 1
#endif

#ifndef GL_ARB_draw_buffers_blend
#define GL_ARB_draw_buffers_blend 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glBlendEquationiARB (GLuint buf, GLenum mode);
GLAPI void APIENTRY glBlendEquationSeparateiARB (GLuint buf, GLenum modeRGB, GLenum modeAlpha);
GLAPI void APIENTRY glBlendFunciARB (GLuint buf, GLenum src, GLenum dst);
GLAPI void APIENTRY glBlendFuncSeparateiARB (GLuint buf, GLenum srcRGB, GLenum dstRGB, GLenum srcAlpha, GLenum dstAlpha);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLBLENDEQUATIONIARBPROC) (GLuint buf, GLenum mode);
typedef void (APIENTRYP PFNGLBLENDEQUATIONSEPARATEIARBPROC) (GLuint buf, GLenum modeRGB, GLenum modeAlpha);
typedef void (APIENTRYP PFNGLBLENDFUNCIARBPROC) (GLuint buf, GLenum src, GLenum dst);
typedef void (APIENTRYP PFNGLBLENDFUNCSEPARATEIARBPROC) (GLuint buf, GLenum srcRGB, GLenum dstRGB, GLenum srcAlpha, GLenum dstAlpha);
#endif

#ifndef GL_ARB_sample_shading
#define GL_ARB_sample_shading 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glMinSampleShadingARB (GLfloat value);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLMINSAMPLESHADINGARBPROC) (GLfloat value);
#endif

#ifndef GL_ARB_texture_cube_map_array
#define GL_ARB_texture_cube_map_array 1
#endif

#ifndef GL_ARB_texture_gather
#define GL_ARB_texture_gather 1
#endif

#ifndef GL_ARB_texture_query_lod
#define GL_ARB_texture_query_lod 1
#endif

#ifndef GL_ARB_shading_language_include
#define GL_ARB_shading_language_include 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glNamedStringARB (GLenum type, GLint namelen, const GLchar *name, GLint stringlen, const GLchar *string);
GLAPI void APIENTRY glDeleteNamedStringARB (GLint namelen, const GLchar *name);
GLAPI void APIENTRY glCompileShaderIncludeARB (GLuint shader, GLsizei count, const GLchar* *path, const GLint *length);
GLAPI GLboolean APIENTRY glIsNamedStringARB (GLint namelen, const GLchar *name);
GLAPI void APIENTRY glGetNamedStringARB (GLint namelen, const GLchar *name, GLsizei bufSize, GLint *stringlen, GLchar *string);
GLAPI void APIENTRY glGetNamedStringivARB (GLint namelen, const GLchar *name, GLenum pname, GLint *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLNAMEDSTRINGARBPROC) (GLenum type, GLint namelen, const GLchar *name, GLint stringlen, const GLchar *string);
typedef void (APIENTRYP PFNGLDELETENAMEDSTRINGARBPROC) (GLint namelen, const GLchar *name);
typedef void (APIENTRYP PFNGLCOMPILESHADERINCLUDEARBPROC) (GLuint shader, GLsizei count, const GLchar* *path, const GLint *length);
typedef GLboolean (APIENTRYP PFNGLISNAMEDSTRINGARBPROC) (GLint namelen, const GLchar *name);
typedef void (APIENTRYP PFNGLGETNAMEDSTRINGARBPROC) (GLint namelen, const GLchar *name, GLsizei bufSize, GLint *stringlen, GLchar *string);
typedef void (APIENTRYP PFNGLGETNAMEDSTRINGIVARBPROC) (GLint namelen, const GLchar *name, GLenum pname, GLint *params);
#endif

#ifndef GL_ARB_texture_compression_bptc
#define GL_ARB_texture_compression_bptc 1
#endif

#ifndef GL_ARB_blend_func_extended
#define GL_ARB_blend_func_extended 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glBindFragDataLocationIndexed (GLuint program, GLuint colorNumber, GLuint index, const GLchar *name);
GLAPI GLint APIENTRY glGetFragDataIndex (GLuint program, const GLchar *name);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLBINDFRAGDATALOCATIONINDEXEDPROC) (GLuint program, GLuint colorNumber, GLuint index, const GLchar *name);
typedef GLint (APIENTRYP PFNGLGETFRAGDATAINDEXPROC) (GLuint program, const GLchar *name);
#endif

#ifndef GL_ARB_explicit_attrib_location
#define GL_ARB_explicit_attrib_location 1
#endif

#ifndef GL_ARB_occlusion_query2
#define GL_ARB_occlusion_query2 1
#endif

#ifndef GL_ARB_sampler_objects
#define GL_ARB_sampler_objects 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glGenSamplers (GLsizei count, GLuint *samplers);
GLAPI void APIENTRY glDeleteSamplers (GLsizei count, const GLuint *samplers);
GLAPI GLboolean APIENTRY glIsSampler (GLuint sampler);
GLAPI void APIENTRY glBindSampler (GLuint unit, GLuint sampler);
GLAPI void APIENTRY glSamplerParameteri (GLuint sampler, GLenum pname, GLint param);
GLAPI void APIENTRY glSamplerParameteriv (GLuint sampler, GLenum pname, const GLint *param);
GLAPI void APIENTRY glSamplerParameterf (GLuint sampler, GLenum pname, GLfloat param);
GLAPI void APIENTRY glSamplerParameterfv (GLuint sampler, GLenum pname, const GLfloat *param);
GLAPI void APIENTRY glSamplerParameterIiv (GLuint sampler, GLenum pname, const GLint *param);
GLAPI void APIENTRY glSamplerParameterIuiv (GLuint sampler, GLenum pname, const GLuint *param);
GLAPI void APIENTRY glGetSamplerParameteriv (GLuint sampler, GLenum pname, GLint *params);
GLAPI void APIENTRY glGetSamplerParameterIiv (GLuint sampler, GLenum pname, GLint *params);
GLAPI void APIENTRY glGetSamplerParameterfv (GLuint sampler, GLenum pname, GLfloat *params);
GLAPI void APIENTRY glGetSamplerParameterIuiv (GLuint sampler, GLenum pname, GLuint *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLGENSAMPLERSPROC) (GLsizei count, GLuint *samplers);
typedef void (APIENTRYP PFNGLDELETESAMPLERSPROC) (GLsizei count, const GLuint *samplers);
typedef GLboolean (APIENTRYP PFNGLISSAMPLERPROC) (GLuint sampler);
typedef void (APIENTRYP PFNGLBINDSAMPLERPROC) (GLuint unit, GLuint sampler);
typedef void (APIENTRYP PFNGLSAMPLERPARAMETERIPROC) (GLuint sampler, GLenum pname, GLint param);
typedef void (APIENTRYP PFNGLSAMPLERPARAMETERIVPROC) (GLuint sampler, GLenum pname, const GLint *param);
typedef void (APIENTRYP PFNGLSAMPLERPARAMETERFPROC) (GLuint sampler, GLenum pname, GLfloat param);
typedef void (APIENTRYP PFNGLSAMPLERPARAMETERFVPROC) (GLuint sampler, GLenum pname, const GLfloat *param);
typedef void (APIENTRYP PFNGLSAMPLERPARAMETERIIVPROC) (GLuint sampler, GLenum pname, const GLint *param);
typedef void (APIENTRYP PFNGLSAMPLERPARAMETERIUIVPROC) (GLuint sampler, GLenum pname, const GLuint *param);
typedef void (APIENTRYP PFNGLGETSAMPLERPARAMETERIVPROC) (GLuint sampler, GLenum pname, GLint *params);
typedef void (APIENTRYP PFNGLGETSAMPLERPARAMETERIIVPROC) (GLuint sampler, GLenum pname, GLint *params);
typedef void (APIENTRYP PFNGLGETSAMPLERPARAMETERFVPROC) (GLuint sampler, GLenum pname, GLfloat *params);
typedef void (APIENTRYP PFNGLGETSAMPLERPARAMETERIUIVPROC) (GLuint sampler, GLenum pname, GLuint *params);
#endif

#ifndef GL_ARB_shader_bit_encoding
#define GL_ARB_shader_bit_encoding 1
#endif

#ifndef GL_ARB_texture_rgb10_a2ui
#define GL_ARB_texture_rgb10_a2ui 1
#endif

#ifndef GL_ARB_texture_swizzle
#define GL_ARB_texture_swizzle 1
#endif

#ifndef GL_ARB_timer_query
#define GL_ARB_timer_query 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glQueryCounter (GLuint id, GLenum target);
GLAPI void APIENTRY glGetQueryObjecti64v (GLuint id, GLenum pname, GLint64 *params);
GLAPI void APIENTRY glGetQueryObjectui64v (GLuint id, GLenum pname, GLuint64 *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLQUERYCOUNTERPROC) (GLuint id, GLenum target);
typedef void (APIENTRYP PFNGLGETQUERYOBJECTI64VPROC) (GLuint id, GLenum pname, GLint64 *params);
typedef void (APIENTRYP PFNGLGETQUERYOBJECTUI64VPROC) (GLuint id, GLenum pname, GLuint64 *params);
#endif

#ifndef GL_ARB_vertex_type_2_10_10_10_rev
#define GL_ARB_vertex_type_2_10_10_10_rev 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glVertexP2ui (GLenum type, GLuint value);
GLAPI void APIENTRY glVertexP2uiv (GLenum type, const GLuint *value);
GLAPI void APIENTRY glVertexP3ui (GLenum type, GLuint value);
GLAPI void APIENTRY glVertexP3uiv (GLenum type, const GLuint *value);
GLAPI void APIENTRY glVertexP4ui (GLenum type, GLuint value);
GLAPI void APIENTRY glVertexP4uiv (GLenum type, const GLuint *value);
GLAPI void APIENTRY glTexCoordP1ui (GLenum type, GLuint coords);
GLAPI void APIENTRY glTexCoordP1uiv (GLenum type, const GLuint *coords);
GLAPI void APIENTRY glTexCoordP2ui (GLenum type, GLuint coords);
GLAPI void APIENTRY glTexCoordP2uiv (GLenum type, const GLuint *coords);
GLAPI void APIENTRY glTexCoordP3ui (GLenum type, GLuint coords);
GLAPI void APIENTRY glTexCoordP3uiv (GLenum type, const GLuint *coords);
GLAPI void APIENTRY glTexCoordP4ui (GLenum type, GLuint coords);
GLAPI void APIENTRY glTexCoordP4uiv (GLenum type, const GLuint *coords);
GLAPI void APIENTRY glMultiTexCoordP1ui (GLenum texture, GLenum type, GLuint coords);
GLAPI void APIENTRY glMultiTexCoordP1uiv (GLenum texture, GLenum type, const GLuint *coords);
GLAPI void APIENTRY glMultiTexCoordP2ui (GLenum texture, GLenum type, GLuint coords);
GLAPI void APIENTRY glMultiTexCoordP2uiv (GLenum texture, GLenum type, const GLuint *coords);
GLAPI void APIENTRY glMultiTexCoordP3ui (GLenum texture, GLenum type, GLuint coords);
GLAPI void APIENTRY glMultiTexCoordP3uiv (GLenum texture, GLenum type, const GLuint *coords);
GLAPI void APIENTRY glMultiTexCoordP4ui (GLenum texture, GLenum type, GLuint coords);
GLAPI void APIENTRY glMultiTexCoordP4uiv (GLenum texture, GLenum type, const GLuint *coords);
GLAPI void APIENTRY glNormalP3ui (GLenum type, GLuint coords);
GLAPI void APIENTRY glNormalP3uiv (GLenum type, const GLuint *coords);
GLAPI void APIENTRY glColorP3ui (GLenum type, GLuint color);
GLAPI void APIENTRY glColorP3uiv (GLenum type, const GLuint *color);
GLAPI void APIENTRY glColorP4ui (GLenum type, GLuint color);
GLAPI void APIENTRY glColorP4uiv (GLenum type, const GLuint *color);
GLAPI void APIENTRY glSecondaryColorP3ui (GLenum type, GLuint color);
GLAPI void APIENTRY glSecondaryColorP3uiv (GLenum type, const GLuint *color);
GLAPI void APIENTRY glVertexAttribP1ui (GLuint index, GLenum type, GLboolean normalized, GLuint value);
GLAPI void APIENTRY glVertexAttribP1uiv (GLuint index, GLenum type, GLboolean normalized, const GLuint *value);
GLAPI void APIENTRY glVertexAttribP2ui (GLuint index, GLenum type, GLboolean normalized, GLuint value);
GLAPI void APIENTRY glVertexAttribP2uiv (GLuint index, GLenum type, GLboolean normalized, const GLuint *value);
GLAPI void APIENTRY glVertexAttribP3ui (GLuint index, GLenum type, GLboolean normalized, GLuint value);
GLAPI void APIENTRY glVertexAttribP3uiv (GLuint index, GLenum type, GLboolean normalized, const GLuint *value);
GLAPI void APIENTRY glVertexAttribP4ui (GLuint index, GLenum type, GLboolean normalized, GLuint value);
GLAPI void APIENTRY glVertexAttribP4uiv (GLuint index, GLenum type, GLboolean normalized, const GLuint *value);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLVERTEXP2UIPROC) (GLenum type, GLuint value);
typedef void (APIENTRYP PFNGLVERTEXP2UIVPROC) (GLenum type, const GLuint *value);
typedef void (APIENTRYP PFNGLVERTEXP3UIPROC) (GLenum type, GLuint value);
typedef void (APIENTRYP PFNGLVERTEXP3UIVPROC) (GLenum type, const GLuint *value);
typedef void (APIENTRYP PFNGLVERTEXP4UIPROC) (GLenum type, GLuint value);
typedef void (APIENTRYP PFNGLVERTEXP4UIVPROC) (GLenum type, const GLuint *value);
typedef void (APIENTRYP PFNGLTEXCOORDP1UIPROC) (GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLTEXCOORDP1UIVPROC) (GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLTEXCOORDP2UIPROC) (GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLTEXCOORDP2UIVPROC) (GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLTEXCOORDP3UIPROC) (GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLTEXCOORDP3UIVPROC) (GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLTEXCOORDP4UIPROC) (GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLTEXCOORDP4UIVPROC) (GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLMULTITEXCOORDP1UIPROC) (GLenum texture, GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLMULTITEXCOORDP1UIVPROC) (GLenum texture, GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLMULTITEXCOORDP2UIPROC) (GLenum texture, GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLMULTITEXCOORDP2UIVPROC) (GLenum texture, GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLMULTITEXCOORDP3UIPROC) (GLenum texture, GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLMULTITEXCOORDP3UIVPROC) (GLenum texture, GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLMULTITEXCOORDP4UIPROC) (GLenum texture, GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLMULTITEXCOORDP4UIVPROC) (GLenum texture, GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLNORMALP3UIPROC) (GLenum type, GLuint coords);
typedef void (APIENTRYP PFNGLNORMALP3UIVPROC) (GLenum type, const GLuint *coords);
typedef void (APIENTRYP PFNGLCOLORP3UIPROC) (GLenum type, GLuint color);
typedef void (APIENTRYP PFNGLCOLORP3UIVPROC) (GLenum type, const GLuint *color);
typedef void (APIENTRYP PFNGLCOLORP4UIPROC) (GLenum type, GLuint color);
typedef void (APIENTRYP PFNGLCOLORP4UIVPROC) (GLenum type, const GLuint *color);
typedef void (APIENTRYP PFNGLSECONDARYCOLORP3UIPROC) (GLenum type, GLuint color);
typedef void (APIENTRYP PFNGLSECONDARYCOLORP3UIVPROC) (GLenum type, const GLuint *color);
typedef void (APIENTRYP PFNGLVERTEXATTRIBP1UIPROC) (GLuint index, GLenum type, GLboolean normalized, GLuint value);
typedef void (APIENTRYP PFNGLVERTEXATTRIBP1UIVPROC) (GLuint index, GLenum type, GLboolean normalized, const GLuint *value);
typedef void (APIENTRYP PFNGLVERTEXATTRIBP2UIPROC) (GLuint index, GLenum type, GLboolean normalized, GLuint value);
typedef void (APIENTRYP PFNGLVERTEXATTRIBP2UIVPROC) (GLuint index, GLenum type, GLboolean normalized, const GLuint *value);
typedef void (APIENTRYP PFNGLVERTEXATTRIBP3UIPROC) (GLuint index, GLenum type, GLboolean normalized, GLuint value);
typedef void (APIENTRYP PFNGLVERTEXATTRIBP3UIVPROC) (GLuint index, GLenum type, GLboolean normalized, const GLuint *value);
typedef void (APIENTRYP PFNGLVERTEXATTRIBP4UIPROC) (GLuint index, GLenum type, GLboolean normalized, GLuint value);
typedef void (APIENTRYP PFNGLVERTEXATTRIBP4UIVPROC) (GLuint index, GLenum type, GLboolean normalized, const GLuint *value);
#endif

#ifndef GL_ARB_draw_indirect
#define GL_ARB_draw_indirect 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDrawArraysIndirect (GLenum mode, const GLvoid *indirect);
GLAPI void APIENTRY glDrawElementsIndirect (GLenum mode, GLenum type, const GLvoid *indirect);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDRAWARRAYSINDIRECTPROC) (GLenum mode, const GLvoid *indirect);
typedef void (APIENTRYP PFNGLDRAWELEMENTSINDIRECTPROC) (GLenum mode, GLenum type, const GLvoid *indirect);
#endif

#ifndef GL_ARB_gpu_shader5
#define GL_ARB_gpu_shader5 1
#endif

#ifndef GL_ARB_gpu_shader_fp64
#define GL_ARB_gpu_shader_fp64 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glUniform1d (GLint location, GLdouble x);
GLAPI void APIENTRY glUniform2d (GLint location, GLdouble x, GLdouble y);
GLAPI void APIENTRY glUniform3d (GLint location, GLdouble x, GLdouble y, GLdouble z);
GLAPI void APIENTRY glUniform4d (GLint location, GLdouble x, GLdouble y, GLdouble z, GLdouble w);
GLAPI void APIENTRY glUniform1dv (GLint location, GLsizei count, const GLdouble *value);
GLAPI void APIENTRY glUniform2dv (GLint location, GLsizei count, const GLdouble *value);
GLAPI void APIENTRY glUniform3dv (GLint location, GLsizei count, const GLdouble *value);
GLAPI void APIENTRY glUniform4dv (GLint location, GLsizei count, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix2dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix3dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix4dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix2x3dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix2x4dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix3x2dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix3x4dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix4x2dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glUniformMatrix4x3dv (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glGetUniformdv (GLuint program, GLint location, GLdouble *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLUNIFORM1DPROC) (GLint location, GLdouble x);
typedef void (APIENTRYP PFNGLUNIFORM2DPROC) (GLint location, GLdouble x, GLdouble y);
typedef void (APIENTRYP PFNGLUNIFORM3DPROC) (GLint location, GLdouble x, GLdouble y, GLdouble z);
typedef void (APIENTRYP PFNGLUNIFORM4DPROC) (GLint location, GLdouble x, GLdouble y, GLdouble z, GLdouble w);
typedef void (APIENTRYP PFNGLUNIFORM1DVPROC) (GLint location, GLsizei count, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORM2DVPROC) (GLint location, GLsizei count, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORM3DVPROC) (GLint location, GLsizei count, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORM4DVPROC) (GLint location, GLsizei count, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX2DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX3DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX4DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX2X3DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX2X4DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX3X2DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX3X4DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX4X2DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLUNIFORMMATRIX4X3DVPROC) (GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLGETUNIFORMDVPROC) (GLuint program, GLint location, GLdouble *params);
#endif

#ifndef GL_ARB_shader_subroutine
#define GL_ARB_shader_subroutine 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI GLint APIENTRY glGetSubroutineUniformLocation (GLuint program, GLenum shadertype, const GLchar *name);
GLAPI GLuint APIENTRY glGetSubroutineIndex (GLuint program, GLenum shadertype, const GLchar *name);
GLAPI void APIENTRY glGetActiveSubroutineUniformiv (GLuint program, GLenum shadertype, GLuint index, GLenum pname, GLint *values);
GLAPI void APIENTRY glGetActiveSubroutineUniformName (GLuint program, GLenum shadertype, GLuint index, GLsizei bufsize, GLsizei *length, GLchar *name);
GLAPI void APIENTRY glGetActiveSubroutineName (GLuint program, GLenum shadertype, GLuint index, GLsizei bufsize, GLsizei *length, GLchar *name);
GLAPI void APIENTRY glUniformSubroutinesuiv (GLenum shadertype, GLsizei count, const GLuint *indices);
GLAPI void APIENTRY glGetUniformSubroutineuiv (GLenum shadertype, GLint location, GLuint *params);
GLAPI void APIENTRY glGetProgramStageiv (GLuint program, GLenum shadertype, GLenum pname, GLint *values);
#endif /* GL_GLEXT_PROTOTYPES */
typedef GLint (APIENTRYP PFNGLGETSUBROUTINEUNIFORMLOCATIONPROC) (GLuint program, GLenum shadertype, const GLchar *name);
typedef GLuint (APIENTRYP PFNGLGETSUBROUTINEINDEXPROC) (GLuint program, GLenum shadertype, const GLchar *name);
typedef void (APIENTRYP PFNGLGETACTIVESUBROUTINEUNIFORMIVPROC) (GLuint program, GLenum shadertype, GLuint index, GLenum pname, GLint *values);
typedef void (APIENTRYP PFNGLGETACTIVESUBROUTINEUNIFORMNAMEPROC) (GLuint program, GLenum shadertype, GLuint index, GLsizei bufsize, GLsizei *length, GLchar *name);
typedef void (APIENTRYP PFNGLGETACTIVESUBROUTINENAMEPROC) (GLuint program, GLenum shadertype, GLuint index, GLsizei bufsize, GLsizei *length, GLchar *name);
typedef void (APIENTRYP PFNGLUNIFORMSUBROUTINESUIVPROC) (GLenum shadertype, GLsizei count, const GLuint *indices);
typedef void (APIENTRYP PFNGLGETUNIFORMSUBROUTINEUIVPROC) (GLenum shadertype, GLint location, GLuint *params);
typedef void (APIENTRYP PFNGLGETPROGRAMSTAGEIVPROC) (GLuint program, GLenum shadertype, GLenum pname, GLint *values);
#endif

#ifndef GL_ARB_tessellation_shader
#define GL_ARB_tessellation_shader 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glPatchParameteri (GLenum pname, GLint value);
GLAPI void APIENTRY glPatchParameterfv (GLenum pname, const GLfloat *values);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLPATCHPARAMETERIPROC) (GLenum pname, GLint value);
typedef void (APIENTRYP PFNGLPATCHPARAMETERFVPROC) (GLenum pname, const GLfloat *values);
#endif

#ifndef GL_ARB_texture_buffer_object_rgb32
#define GL_ARB_texture_buffer_object_rgb32 1
#endif

#ifndef GL_ARB_transform_feedback2
#define GL_ARB_transform_feedback2 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glBindTransformFeedback (GLenum target, GLuint id);
GLAPI void APIENTRY glDeleteTransformFeedbacks (GLsizei n, const GLuint *ids);
GLAPI void APIENTRY glGenTransformFeedbacks (GLsizei n, GLuint *ids);
GLAPI GLboolean APIENTRY glIsTransformFeedback (GLuint id);
GLAPI void APIENTRY glPauseTransformFeedback (void);
GLAPI void APIENTRY glResumeTransformFeedback (void);
GLAPI void APIENTRY glDrawTransformFeedback (GLenum mode, GLuint id);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLBINDTRANSFORMFEEDBACKPROC) (GLenum target, GLuint id);
typedef void (APIENTRYP PFNGLDELETETRANSFORMFEEDBACKSPROC) (GLsizei n, const GLuint *ids);
typedef void (APIENTRYP PFNGLGENTRANSFORMFEEDBACKSPROC) (GLsizei n, GLuint *ids);
typedef GLboolean (APIENTRYP PFNGLISTRANSFORMFEEDBACKPROC) (GLuint id);
typedef void (APIENTRYP PFNGLPAUSETRANSFORMFEEDBACKPROC) (void);
typedef void (APIENTRYP PFNGLRESUMETRANSFORMFEEDBACKPROC) (void);
typedef void (APIENTRYP PFNGLDRAWTRANSFORMFEEDBACKPROC) (GLenum mode, GLuint id);
#endif

#ifndef GL_ARB_transform_feedback3
#define GL_ARB_transform_feedback3 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDrawTransformFeedbackStream (GLenum mode, GLuint id, GLuint stream);
GLAPI void APIENTRY glBeginQueryIndexed (GLenum target, GLuint index, GLuint id);
GLAPI void APIENTRY glEndQueryIndexed (GLenum target, GLuint index);
GLAPI void APIENTRY glGetQueryIndexediv (GLenum target, GLuint index, GLenum pname, GLint *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDRAWTRANSFORMFEEDBACKSTREAMPROC) (GLenum mode, GLuint id, GLuint stream);
typedef void (APIENTRYP PFNGLBEGINQUERYINDEXEDPROC) (GLenum target, GLuint index, GLuint id);
typedef void (APIENTRYP PFNGLENDQUERYINDEXEDPROC) (GLenum target, GLuint index);
typedef void (APIENTRYP PFNGLGETQUERYINDEXEDIVPROC) (GLenum target, GLuint index, GLenum pname, GLint *params);
#endif

#ifndef GL_ARB_ES2_compatibility
#define GL_ARB_ES2_compatibility 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glReleaseShaderCompiler (void);
GLAPI void APIENTRY glShaderBinary (GLsizei count, const GLuint *shaders, GLenum binaryformat, const GLvoid *binary, GLsizei length);
GLAPI void APIENTRY glGetShaderPrecisionFormat (GLenum shadertype, GLenum precisiontype, GLint *range, GLint *precision);
GLAPI void APIENTRY glDepthRangef (GLfloat n, GLfloat f);
GLAPI void APIENTRY glClearDepthf (GLfloat d);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLRELEASESHADERCOMPILERPROC) (void);
typedef void (APIENTRYP PFNGLSHADERBINARYPROC) (GLsizei count, const GLuint *shaders, GLenum binaryformat, const GLvoid *binary, GLsizei length);
typedef void (APIENTRYP PFNGLGETSHADERPRECISIONFORMATPROC) (GLenum shadertype, GLenum precisiontype, GLint *range, GLint *precision);
typedef void (APIENTRYP PFNGLDEPTHRANGEFPROC) (GLfloat n, GLfloat f);
typedef void (APIENTRYP PFNGLCLEARDEPTHFPROC) (GLfloat d);
#endif

#ifndef GL_ARB_get_program_binary
#define GL_ARB_get_program_binary 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glGetProgramBinary (GLuint program, GLsizei bufSize, GLsizei *length, GLenum *binaryFormat, GLvoid *binary);
GLAPI void APIENTRY glProgramBinary (GLuint program, GLenum binaryFormat, const GLvoid *binary, GLsizei length);
GLAPI void APIENTRY glProgramParameteri (GLuint program, GLenum pname, GLint value);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLGETPROGRAMBINARYPROC) (GLuint program, GLsizei bufSize, GLsizei *length, GLenum *binaryFormat, GLvoid *binary);
typedef void (APIENTRYP PFNGLPROGRAMBINARYPROC) (GLuint program, GLenum binaryFormat, const GLvoid *binary, GLsizei length);
typedef void (APIENTRYP PFNGLPROGRAMPARAMETERIPROC) (GLuint program, GLenum pname, GLint value);
#endif

#ifndef GL_ARB_separate_shader_objects
#define GL_ARB_separate_shader_objects 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glUseProgramStages (GLuint pipeline, GLbitfield stages, GLuint program);
GLAPI void APIENTRY glActiveShaderProgram (GLuint pipeline, GLuint program);
GLAPI GLuint APIENTRY glCreateShaderProgramv (GLenum type, GLsizei count, const GLchar* const *strings);
GLAPI void APIENTRY glBindProgramPipeline (GLuint pipeline);
GLAPI void APIENTRY glDeleteProgramPipelines (GLsizei n, const GLuint *pipelines);
GLAPI void APIENTRY glGenProgramPipelines (GLsizei n, GLuint *pipelines);
GLAPI GLboolean APIENTRY glIsProgramPipeline (GLuint pipeline);
GLAPI void APIENTRY glGetProgramPipelineiv (GLuint pipeline, GLenum pname, GLint *params);
GLAPI void APIENTRY glProgramUniform1i (GLuint program, GLint location, GLint v0);
GLAPI void APIENTRY glProgramUniform1iv (GLuint program, GLint location, GLsizei count, const GLint *value);
GLAPI void APIENTRY glProgramUniform1f (GLuint program, GLint location, GLfloat v0);
GLAPI void APIENTRY glProgramUniform1fv (GLuint program, GLint location, GLsizei count, const GLfloat *value);
GLAPI void APIENTRY glProgramUniform1d (GLuint program, GLint location, GLdouble v0);
GLAPI void APIENTRY glProgramUniform1dv (GLuint program, GLint location, GLsizei count, const GLdouble *value);
GLAPI void APIENTRY glProgramUniform1ui (GLuint program, GLint location, GLuint v0);
GLAPI void APIENTRY glProgramUniform1uiv (GLuint program, GLint location, GLsizei count, const GLuint *value);
GLAPI void APIENTRY glProgramUniform2i (GLuint program, GLint location, GLint v0, GLint v1);
GLAPI void APIENTRY glProgramUniform2iv (GLuint program, GLint location, GLsizei count, const GLint *value);
GLAPI void APIENTRY glProgramUniform2f (GLuint program, GLint location, GLfloat v0, GLfloat v1);
GLAPI void APIENTRY glProgramUniform2fv (GLuint program, GLint location, GLsizei count, const GLfloat *value);
GLAPI void APIENTRY glProgramUniform2d (GLuint program, GLint location, GLdouble v0, GLdouble v1);
GLAPI void APIENTRY glProgramUniform2dv (GLuint program, GLint location, GLsizei count, const GLdouble *value);
GLAPI void APIENTRY glProgramUniform2ui (GLuint program, GLint location, GLuint v0, GLuint v1);
GLAPI void APIENTRY glProgramUniform2uiv (GLuint program, GLint location, GLsizei count, const GLuint *value);
GLAPI void APIENTRY glProgramUniform3i (GLuint program, GLint location, GLint v0, GLint v1, GLint v2);
GLAPI void APIENTRY glProgramUniform3iv (GLuint program, GLint location, GLsizei count, const GLint *value);
GLAPI void APIENTRY glProgramUniform3f (GLuint program, GLint location, GLfloat v0, GLfloat v1, GLfloat v2);
GLAPI void APIENTRY glProgramUniform3fv (GLuint program, GLint location, GLsizei count, const GLfloat *value);
GLAPI void APIENTRY glProgramUniform3d (GLuint program, GLint location, GLdouble v0, GLdouble v1, GLdouble v2);
GLAPI void APIENTRY glProgramUniform3dv (GLuint program, GLint location, GLsizei count, const GLdouble *value);
GLAPI void APIENTRY glProgramUniform3ui (GLuint program, GLint location, GLuint v0, GLuint v1, GLuint v2);
GLAPI void APIENTRY glProgramUniform3uiv (GLuint program, GLint location, GLsizei count, const GLuint *value);
GLAPI void APIENTRY glProgramUniform4i (GLuint program, GLint location, GLint v0, GLint v1, GLint v2, GLint v3);
GLAPI void APIENTRY glProgramUniform4iv (GLuint program, GLint location, GLsizei count, const GLint *value);
GLAPI void APIENTRY glProgramUniform4f (GLuint program, GLint location, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3);
GLAPI void APIENTRY glProgramUniform4fv (GLuint program, GLint location, GLsizei count, const GLfloat *value);
GLAPI void APIENTRY glProgramUniform4d (GLuint program, GLint location, GLdouble v0, GLdouble v1, GLdouble v2, GLdouble v3);
GLAPI void APIENTRY glProgramUniform4dv (GLuint program, GLint location, GLsizei count, const GLdouble *value);
GLAPI void APIENTRY glProgramUniform4ui (GLuint program, GLint location, GLuint v0, GLuint v1, GLuint v2, GLuint v3);
GLAPI void APIENTRY glProgramUniform4uiv (GLuint program, GLint location, GLsizei count, const GLuint *value);
GLAPI void APIENTRY glProgramUniformMatrix2fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix3fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix4fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix2dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glProgramUniformMatrix3dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glProgramUniformMatrix4dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glProgramUniformMatrix2x3fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix3x2fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix2x4fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix4x2fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix3x4fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix4x3fv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
GLAPI void APIENTRY glProgramUniformMatrix2x3dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glProgramUniformMatrix3x2dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glProgramUniformMatrix2x4dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glProgramUniformMatrix4x2dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glProgramUniformMatrix3x4dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glProgramUniformMatrix4x3dv (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
GLAPI void APIENTRY glValidateProgramPipeline (GLuint pipeline);
GLAPI void APIENTRY glGetProgramPipelineInfoLog (GLuint pipeline, GLsizei bufSize, GLsizei *length, GLchar *infoLog);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLUSEPROGRAMSTAGESPROC) (GLuint pipeline, GLbitfield stages, GLuint program);
typedef void (APIENTRYP PFNGLACTIVESHADERPROGRAMPROC) (GLuint pipeline, GLuint program);
typedef GLuint (APIENTRYP PFNGLCREATESHADERPROGRAMVPROC) (GLenum type, GLsizei count, const GLchar* const *strings);
typedef void (APIENTRYP PFNGLBINDPROGRAMPIPELINEPROC) (GLuint pipeline);
typedef void (APIENTRYP PFNGLDELETEPROGRAMPIPELINESPROC) (GLsizei n, const GLuint *pipelines);
typedef void (APIENTRYP PFNGLGENPROGRAMPIPELINESPROC) (GLsizei n, GLuint *pipelines);
typedef GLboolean (APIENTRYP PFNGLISPROGRAMPIPELINEPROC) (GLuint pipeline);
typedef void (APIENTRYP PFNGLGETPROGRAMPIPELINEIVPROC) (GLuint pipeline, GLenum pname, GLint *params);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM1IPROC) (GLuint program, GLint location, GLint v0);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM1IVPROC) (GLuint program, GLint location, GLsizei count, const GLint *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM1FPROC) (GLuint program, GLint location, GLfloat v0);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM1FVPROC) (GLuint program, GLint location, GLsizei count, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM1DPROC) (GLuint program, GLint location, GLdouble v0);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM1DVPROC) (GLuint program, GLint location, GLsizei count, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM1UIPROC) (GLuint program, GLint location, GLuint v0);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM1UIVPROC) (GLuint program, GLint location, GLsizei count, const GLuint *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM2IPROC) (GLuint program, GLint location, GLint v0, GLint v1);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM2IVPROC) (GLuint program, GLint location, GLsizei count, const GLint *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM2FPROC) (GLuint program, GLint location, GLfloat v0, GLfloat v1);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM2FVPROC) (GLuint program, GLint location, GLsizei count, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM2DPROC) (GLuint program, GLint location, GLdouble v0, GLdouble v1);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM2DVPROC) (GLuint program, GLint location, GLsizei count, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM2UIPROC) (GLuint program, GLint location, GLuint v0, GLuint v1);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM2UIVPROC) (GLuint program, GLint location, GLsizei count, const GLuint *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM3IPROC) (GLuint program, GLint location, GLint v0, GLint v1, GLint v2);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM3IVPROC) (GLuint program, GLint location, GLsizei count, const GLint *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM3FPROC) (GLuint program, GLint location, GLfloat v0, GLfloat v1, GLfloat v2);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM3FVPROC) (GLuint program, GLint location, GLsizei count, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM3DPROC) (GLuint program, GLint location, GLdouble v0, GLdouble v1, GLdouble v2);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM3DVPROC) (GLuint program, GLint location, GLsizei count, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM3UIPROC) (GLuint program, GLint location, GLuint v0, GLuint v1, GLuint v2);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM3UIVPROC) (GLuint program, GLint location, GLsizei count, const GLuint *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM4IPROC) (GLuint program, GLint location, GLint v0, GLint v1, GLint v2, GLint v3);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM4IVPROC) (GLuint program, GLint location, GLsizei count, const GLint *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM4FPROC) (GLuint program, GLint location, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM4FVPROC) (GLuint program, GLint location, GLsizei count, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM4DPROC) (GLuint program, GLint location, GLdouble v0, GLdouble v1, GLdouble v2, GLdouble v3);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM4DVPROC) (GLuint program, GLint location, GLsizei count, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM4UIPROC) (GLuint program, GLint location, GLuint v0, GLuint v1, GLuint v2, GLuint v3);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORM4UIVPROC) (GLuint program, GLint location, GLsizei count, const GLuint *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX2FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX3FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX4FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX2DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX3DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX4DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX2X3FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX3X2FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX2X4FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX4X2FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX3X4FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX4X3FVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX2X3DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX3X2DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX2X4DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX4X2DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX3X4DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLPROGRAMUNIFORMMATRIX4X3DVPROC) (GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble *value);
typedef void (APIENTRYP PFNGLVALIDATEPROGRAMPIPELINEPROC) (GLuint pipeline);
typedef void (APIENTRYP PFNGLGETPROGRAMPIPELINEINFOLOGPROC) (GLuint pipeline, GLsizei bufSize, GLsizei *length, GLchar *infoLog);
#endif

#ifndef GL_ARB_vertex_attrib_64bit
#define GL_ARB_vertex_attrib_64bit 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glVertexAttribL1d (GLuint index, GLdouble x);
GLAPI void APIENTRY glVertexAttribL2d (GLuint index, GLdouble x, GLdouble y);
GLAPI void APIENTRY glVertexAttribL3d (GLuint index, GLdouble x, GLdouble y, GLdouble z);
GLAPI void APIENTRY glVertexAttribL4d (GLuint index, GLdouble x, GLdouble y, GLdouble z, GLdouble w);
GLAPI void APIENTRY glVertexAttribL1dv (GLuint index, const GLdouble *v);
GLAPI void APIENTRY glVertexAttribL2dv (GLuint index, const GLdouble *v);
GLAPI void APIENTRY glVertexAttribL3dv (GLuint index, const GLdouble *v);
GLAPI void APIENTRY glVertexAttribL4dv (GLuint index, const GLdouble *v);
GLAPI void APIENTRY glVertexAttribLPointer (GLuint index, GLint size, GLenum type, GLsizei stride, const GLvoid *pointer);
GLAPI void APIENTRY glGetVertexAttribLdv (GLuint index, GLenum pname, GLdouble *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLVERTEXATTRIBL1DPROC) (GLuint index, GLdouble x);
typedef void (APIENTRYP PFNGLVERTEXATTRIBL2DPROC) (GLuint index, GLdouble x, GLdouble y);
typedef void (APIENTRYP PFNGLVERTEXATTRIBL3DPROC) (GLuint index, GLdouble x, GLdouble y, GLdouble z);
typedef void (APIENTRYP PFNGLVERTEXATTRIBL4DPROC) (GLuint index, GLdouble x, GLdouble y, GLdouble z, GLdouble w);
typedef void (APIENTRYP PFNGLVERTEXATTRIBL1DVPROC) (GLuint index, const GLdouble *v);
typedef void (APIENTRYP PFNGLVERTEXATTRIBL2DVPROC) (GLuint index, const GLdouble *v);
typedef void (APIENTRYP PFNGLVERTEXATTRIBL3DVPROC) (GLuint index, const GLdouble *v);
typedef void (APIENTRYP PFNGLVERTEXATTRIBL4DVPROC) (GLuint index, const GLdouble *v);
typedef void (APIENTRYP PFNGLVERTEXATTRIBLPOINTERPROC) (GLuint index, GLint size, GLenum type, GLsizei stride, const GLvoid *pointer);
typedef void (APIENTRYP PFNGLGETVERTEXATTRIBLDVPROC) (GLuint index, GLenum pname, GLdouble *params);
#endif

#ifndef GL_ARB_viewport_array
#define GL_ARB_viewport_array 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glViewportArrayv (GLuint first, GLsizei count, const GLfloat *v);
GLAPI void APIENTRY glViewportIndexedf (GLuint index, GLfloat x, GLfloat y, GLfloat w, GLfloat h);
GLAPI void APIENTRY glViewportIndexedfv (GLuint index, const GLfloat *v);
GLAPI void APIENTRY glScissorArrayv (GLuint first, GLsizei count, const GLint *v);
GLAPI void APIENTRY glScissorIndexed (GLuint index, GLint left, GLint bottom, GLsizei width, GLsizei height);
GLAPI void APIENTRY glScissorIndexedv (GLuint index, const GLint *v);
GLAPI void APIENTRY glDepthRangeArrayv (GLuint first, GLsizei count, const GLdouble *v);
GLAPI void APIENTRY glDepthRangeIndexed (GLuint index, GLdouble n, GLdouble f);
GLAPI void APIENTRY glGetFloati_v (GLenum target, GLuint index, GLfloat *data);
GLAPI void APIENTRY glGetDoublei_v (GLenum target, GLuint index, GLdouble *data);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLVIEWPORTARRAYVPROC) (GLuint first, GLsizei count, const GLfloat *v);
typedef void (APIENTRYP PFNGLVIEWPORTINDEXEDFPROC) (GLuint index, GLfloat x, GLfloat y, GLfloat w, GLfloat h);
typedef void (APIENTRYP PFNGLVIEWPORTINDEXEDFVPROC) (GLuint index, const GLfloat *v);
typedef void (APIENTRYP PFNGLSCISSORARRAYVPROC) (GLuint first, GLsizei count, const GLint *v);
typedef void (APIENTRYP PFNGLSCISSORINDEXEDPROC) (GLuint index, GLint left, GLint bottom, GLsizei width, GLsizei height);
typedef void (APIENTRYP PFNGLSCISSORINDEXEDVPROC) (GLuint index, const GLint *v);
typedef void (APIENTRYP PFNGLDEPTHRANGEARRAYVPROC) (GLuint first, GLsizei count, const GLdouble *v);
typedef void (APIENTRYP PFNGLDEPTHRANGEINDEXEDPROC) (GLuint index, GLdouble n, GLdouble f);
typedef void (APIENTRYP PFNGLGETFLOATI_VPROC) (GLenum target, GLuint index, GLfloat *data);
typedef void (APIENTRYP PFNGLGETDOUBLEI_VPROC) (GLenum target, GLuint index, GLdouble *data);
#endif

#ifndef GL_ARB_cl_event
#define GL_ARB_cl_event 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI GLsync APIENTRY glCreateSyncFromCLeventARB (struct _cl_context * context, struct _cl_event * event, GLbitfield flags);
#endif /* GL_GLEXT_PROTOTYPES */
typedef GLsync (APIENTRYP PFNGLCREATESYNCFROMCLEVENTARBPROC) (struct _cl_context * context, struct _cl_event * event, GLbitfield flags);
#endif

#ifndef GL_ARB_debug_output
#define GL_ARB_debug_output 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDebugMessageControlARB (GLenum source, GLenum type, GLenum severity, GLsizei count, const GLuint *ids, GLboolean enabled);
GLAPI void APIENTRY glDebugMessageInsertARB (GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *buf);
GLAPI void APIENTRY glDebugMessageCallbackARB (GLDEBUGPROCARB callback, const GLvoid *userParam);
GLAPI GLuint APIENTRY glGetDebugMessageLogARB (GLuint count, GLsizei bufsize, GLenum *sources, GLenum *types, GLuint *ids, GLenum *severities, GLsizei *lengths, GLchar *messageLog);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDEBUGMESSAGECONTROLARBPROC) (GLenum source, GLenum type, GLenum severity, GLsizei count, const GLuint *ids, GLboolean enabled);
typedef void (APIENTRYP PFNGLDEBUGMESSAGEINSERTARBPROC) (GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *buf);
typedef void (APIENTRYP PFNGLDEBUGMESSAGECALLBACKARBPROC) (GLDEBUGPROCARB callback, const GLvoid *userParam);
typedef GLuint (APIENTRYP PFNGLGETDEBUGMESSAGELOGARBPROC) (GLuint count, GLsizei bufsize, GLenum *sources, GLenum *types, GLuint *ids, GLenum *severities, GLsizei *lengths, GLchar *messageLog);
#endif

#ifndef GL_ARB_robustness
#define GL_ARB_robustness 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI GLenum APIENTRY glGetGraphicsResetStatusARB (void);
GLAPI void APIENTRY glGetnMapdvARB (GLenum target, GLenum query, GLsizei bufSize, GLdouble *v);
GLAPI void APIENTRY glGetnMapfvARB (GLenum target, GLenum query, GLsizei bufSize, GLfloat *v);
GLAPI void APIENTRY glGetnMapivARB (GLenum target, GLenum query, GLsizei bufSize, GLint *v);
GLAPI void APIENTRY glGetnPixelMapfvARB (GLenum map, GLsizei bufSize, GLfloat *values);
GLAPI void APIENTRY glGetnPixelMapuivARB (GLenum map, GLsizei bufSize, GLuint *values);
GLAPI void APIENTRY glGetnPixelMapusvARB (GLenum map, GLsizei bufSize, GLushort *values);
GLAPI void APIENTRY glGetnPolygonStippleARB (GLsizei bufSize, GLubyte *pattern);
GLAPI void APIENTRY glGetnColorTableARB (GLenum target, GLenum format, GLenum type, GLsizei bufSize, GLvoid *table);
GLAPI void APIENTRY glGetnConvolutionFilterARB (GLenum target, GLenum format, GLenum type, GLsizei bufSize, GLvoid *image);
GLAPI void APIENTRY glGetnSeparableFilterARB (GLenum target, GLenum format, GLenum type, GLsizei rowBufSize, GLvoid *row, GLsizei columnBufSize, GLvoid *column, GLvoid *span);
GLAPI void APIENTRY glGetnHistogramARB (GLenum target, GLboolean reset, GLenum format, GLenum type, GLsizei bufSize, GLvoid *values);
GLAPI void APIENTRY glGetnMinmaxARB (GLenum target, GLboolean reset, GLenum format, GLenum type, GLsizei bufSize, GLvoid *values);
GLAPI void APIENTRY glGetnTexImageARB (GLenum target, GLint level, GLenum format, GLenum type, GLsizei bufSize, GLvoid *img);
GLAPI void APIENTRY glReadnPixelsARB (GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLsizei bufSize, GLvoid *data);
GLAPI void APIENTRY glGetnCompressedTexImageARB (GLenum target, GLint lod, GLsizei bufSize, GLvoid *img);
GLAPI void APIENTRY glGetnUniformfvARB (GLuint program, GLint location, GLsizei bufSize, GLfloat *params);
GLAPI void APIENTRY glGetnUniformivARB (GLuint program, GLint location, GLsizei bufSize, GLint *params);
GLAPI void APIENTRY glGetnUniformuivARB (GLuint program, GLint location, GLsizei bufSize, GLuint *params);
GLAPI void APIENTRY glGetnUniformdvARB (GLuint program, GLint location, GLsizei bufSize, GLdouble *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef GLenum (APIENTRYP PFNGLGETGRAPHICSRESETSTATUSARBPROC) (void);
typedef void (APIENTRYP PFNGLGETNMAPDVARBPROC) (GLenum target, GLenum query, GLsizei bufSize, GLdouble *v);
typedef void (APIENTRYP PFNGLGETNMAPFVARBPROC) (GLenum target, GLenum query, GLsizei bufSize, GLfloat *v);
typedef void (APIENTRYP PFNGLGETNMAPIVARBPROC) (GLenum target, GLenum query, GLsizei bufSize, GLint *v);
typedef void (APIENTRYP PFNGLGETNPIXELMAPFVARBPROC) (GLenum map, GLsizei bufSize, GLfloat *values);
typedef void (APIENTRYP PFNGLGETNPIXELMAPUIVARBPROC) (GLenum map, GLsizei bufSize, GLuint *values);
typedef void (APIENTRYP PFNGLGETNPIXELMAPUSVARBPROC) (GLenum map, GLsizei bufSize, GLushort *values);
typedef void (APIENTRYP PFNGLGETNPOLYGONSTIPPLEARBPROC) (GLsizei bufSize, GLubyte *pattern);
typedef void (APIENTRYP PFNGLGETNCOLORTABLEARBPROC) (GLenum target, GLenum format, GLenum type, GLsizei bufSize, GLvoid *table);
typedef void (APIENTRYP PFNGLGETNCONVOLUTIONFILTERARBPROC) (GLenum target, GLenum format, GLenum type, GLsizei bufSize, GLvoid *image);
typedef void (APIENTRYP PFNGLGETNSEPARABLEFILTERARBPROC) (GLenum target, GLenum format, GLenum type, GLsizei rowBufSize, GLvoid *row, GLsizei columnBufSize, GLvoid *column, GLvoid *span);
typedef void (APIENTRYP PFNGLGETNHISTOGRAMARBPROC) (GLenum target, GLboolean reset, GLenum format, GLenum type, GLsizei bufSize, GLvoid *values);
typedef void (APIENTRYP PFNGLGETNMINMAXARBPROC) (GLenum target, GLboolean reset, GLenum format, GLenum type, GLsizei bufSize, GLvoid *values);
typedef void (APIENTRYP PFNGLGETNTEXIMAGEARBPROC) (GLenum target, GLint level, GLenum format, GLenum type, GLsizei bufSize, GLvoid *img);
typedef void (APIENTRYP PFNGLREADNPIXELSARBPROC) (GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLsizei bufSize, GLvoid *data);
typedef void (APIENTRYP PFNGLGETNCOMPRESSEDTEXIMAGEARBPROC) (GLenum target, GLint lod, GLsizei bufSize, GLvoid *img);
typedef void (APIENTRYP PFNGLGETNUNIFORMFVARBPROC) (GLuint program, GLint location, GLsizei bufSize, GLfloat *params);
typedef void (APIENTRYP PFNGLGETNUNIFORMIVARBPROC) (GLuint program, GLint location, GLsizei bufSize, GLint *params);
typedef void (APIENTRYP PFNGLGETNUNIFORMUIVARBPROC) (GLuint program, GLint location, GLsizei bufSize, GLuint *params);
typedef void (APIENTRYP PFNGLGETNUNIFORMDVARBPROC) (GLuint program, GLint location, GLsizei bufSize, GLdouble *params);
#endif

#ifndef GL_ARB_shader_stencil_export
#define GL_ARB_shader_stencil_export 1
#endif

#ifndef GL_ARB_base_instance
#define GL_ARB_base_instance 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDrawArraysInstancedBaseInstance (GLenum mode, GLint first, GLsizei count, GLsizei instancecount, GLuint baseinstance);
GLAPI void APIENTRY glDrawElementsInstancedBaseInstance (GLenum mode, GLsizei count, GLenum type, const void *indices, GLsizei instancecount, GLuint baseinstance);
GLAPI void APIENTRY glDrawElementsInstancedBaseVertexBaseInstance (GLenum mode, GLsizei count, GLenum type, const void *indices, GLsizei instancecount, GLint basevertex, GLuint baseinstance);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDBASEINSTANCEPROC) (GLenum mode, GLint first, GLsizei count, GLsizei instancecount, GLuint baseinstance);
typedef void (APIENTRYP PFNGLDRAWELEMENTSINSTANCEDBASEINSTANCEPROC) (GLenum mode, GLsizei count, GLenum type, const void *indices, GLsizei instancecount, GLuint baseinstance);
typedef void (APIENTRYP PFNGLDRAWELEMENTSINSTANCEDBASEVERTEXBASEINSTANCEPROC) (GLenum mode, GLsizei count, GLenum type, const void *indices, GLsizei instancecount, GLint basevertex, GLuint baseinstance);
#endif

#ifndef GL_ARB_shading_language_420pack
#define GL_ARB_shading_language_420pack 1
#endif

#ifndef GL_ARB_transform_feedback_instanced
#define GL_ARB_transform_feedback_instanced 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDrawTransformFeedbackInstanced (GLenum mode, GLuint id, GLsizei instancecount);
GLAPI void APIENTRY glDrawTransformFeedbackStreamInstanced (GLenum mode, GLuint id, GLuint stream, GLsizei instancecount);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDRAWTRANSFORMFEEDBACKINSTANCEDPROC) (GLenum mode, GLuint id, GLsizei instancecount);
typedef void (APIENTRYP PFNGLDRAWTRANSFORMFEEDBACKSTREAMINSTANCEDPROC) (GLenum mode, GLuint id, GLuint stream, GLsizei instancecount);
#endif

#ifndef GL_ARB_compressed_texture_pixel_storage
#define GL_ARB_compressed_texture_pixel_storage 1
#endif

#ifndef GL_ARB_conservative_depth
#define GL_ARB_conservative_depth 1
#endif

#ifndef GL_ARB_internalformat_query
#define GL_ARB_internalformat_query 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glGetInternalformativ (GLenum target, GLenum internalformat, GLenum pname, GLsizei bufSize, GLint *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLGETINTERNALFORMATIVPROC) (GLenum target, GLenum internalformat, GLenum pname, GLsizei bufSize, GLint *params);
#endif

#ifndef GL_ARB_map_buffer_alignment
#define GL_ARB_map_buffer_alignment 1
#endif

#ifndef GL_ARB_shader_atomic_counters
#define GL_ARB_shader_atomic_counters 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glGetActiveAtomicCounterBufferiv (GLuint program, GLuint bufferIndex, GLenum pname, GLint *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLGETACTIVEATOMICCOUNTERBUFFERIVPROC) (GLuint program, GLuint bufferIndex, GLenum pname, GLint *params);
#endif

#ifndef GL_ARB_shader_image_load_store
#define GL_ARB_shader_image_load_store 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glBindImageTexture (GLuint unit, GLuint texture, GLint level, GLboolean layered, GLint layer, GLenum access, GLenum format);
GLAPI void APIENTRY glMemoryBarrier (GLbitfield barriers);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLBINDIMAGETEXTUREPROC) (GLuint unit, GLuint texture, GLint level, GLboolean layered, GLint layer, GLenum access, GLenum format);
typedef void (APIENTRYP PFNGLMEMORYBARRIERPROC) (GLbitfield barriers);
#endif

#ifndef GL_ARB_shading_language_packing
#define GL_ARB_shading_language_packing 1
#endif

#ifndef GL_ARB_texture_storage
#define GL_ARB_texture_storage 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glTexStorage1D (GLenum target, GLsizei levels, GLenum internalformat, GLsizei width);
GLAPI void APIENTRY glTexStorage2D (GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height);
GLAPI void APIENTRY glTexStorage3D (GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth);
GLAPI void APIENTRY glTextureStorage1DEXT (GLuint texture, GLenum target, GLsizei levels, GLenum internalformat, GLsizei width);
GLAPI void APIENTRY glTextureStorage2DEXT (GLuint texture, GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height);
GLAPI void APIENTRY glTextureStorage3DEXT (GLuint texture, GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLTEXSTORAGE1DPROC) (GLenum target, GLsizei levels, GLenum internalformat, GLsizei width);
typedef void (APIENTRYP PFNGLTEXSTORAGE2DPROC) (GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height);
typedef void (APIENTRYP PFNGLTEXSTORAGE3DPROC) (GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth);
typedef void (APIENTRYP PFNGLTEXTURESTORAGE1DEXTPROC) (GLuint texture, GLenum target, GLsizei levels, GLenum internalformat, GLsizei width);
typedef void (APIENTRYP PFNGLTEXTURESTORAGE2DEXTPROC) (GLuint texture, GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height);
typedef void (APIENTRYP PFNGLTEXTURESTORAGE3DEXTPROC) (GLuint texture, GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth);
#endif

#ifndef GL_KHR_texture_compression_astc_ldr
#define GL_KHR_texture_compression_astc_ldr 1
#endif

#ifndef GL_KHR_debug
#define GL_KHR_debug 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDebugMessageControl (GLenum source, GLenum type, GLenum severity, GLsizei count, const GLuint *ids, GLboolean enabled);
GLAPI void APIENTRY glDebugMessageInsert (GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *buf);
GLAPI void APIENTRY glDebugMessageCallback (GLDEBUGPROC callback, const void *userParam);
GLAPI GLuint APIENTRY glGetDebugMessageLog (GLuint count, GLsizei bufsize, GLenum *sources, GLenum *types, GLuint *ids, GLenum *severities, GLsizei *lengths, GLchar *messageLog);
GLAPI void APIENTRY glPushDebugGroup (GLenum source, GLuint id, GLsizei length, const GLchar *message);
GLAPI void APIENTRY glPopDebugGroup (void);
GLAPI void APIENTRY glObjectLabel (GLenum identifier, GLuint name, GLsizei length, const GLchar *label);
GLAPI void APIENTRY glGetObjectLabel (GLenum identifier, GLuint name, GLsizei bufSize, GLsizei *length, GLchar *label);
GLAPI void APIENTRY glObjectPtrLabel (const void *ptr, GLsizei length, const GLchar *label);
GLAPI void APIENTRY glGetObjectPtrLabel (const void *ptr, GLsizei bufSize, GLsizei *length, GLchar *label);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDEBUGMESSAGECONTROLPROC) (GLenum source, GLenum type, GLenum severity, GLsizei count, const GLuint *ids, GLboolean enabled);
typedef void (APIENTRYP PFNGLDEBUGMESSAGEINSERTPROC) (GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar *buf);
typedef void (APIENTRYP PFNGLDEBUGMESSAGECALLBACKPROC) (GLDEBUGPROC callback, const void *userParam);
typedef GLuint (APIENTRYP PFNGLGETDEBUGMESSAGELOGPROC) (GLuint count, GLsizei bufsize, GLenum *sources, GLenum *types, GLuint *ids, GLenum *severities, GLsizei *lengths, GLchar *messageLog);
typedef void (APIENTRYP PFNGLPUSHDEBUGGROUPPROC) (GLenum source, GLuint id, GLsizei length, const GLchar *message);
typedef void (APIENTRYP PFNGLPOPDEBUGGROUPPROC) (void);
typedef void (APIENTRYP PFNGLOBJECTLABELPROC) (GLenum identifier, GLuint name, GLsizei length, const GLchar *label);
typedef void (APIENTRYP PFNGLGETOBJECTLABELPROC) (GLenum identifier, GLuint name, GLsizei bufSize, GLsizei *length, GLchar *label);
typedef void (APIENTRYP PFNGLOBJECTPTRLABELPROC) (const void *ptr, GLsizei length, const GLchar *label);
typedef void (APIENTRYP PFNGLGETOBJECTPTRLABELPROC) (const void *ptr, GLsizei bufSize, GLsizei *length, GLchar *label);
#endif

#ifndef GL_ARB_arrays_of_arrays
#define GL_ARB_arrays_of_arrays 1
#endif

#ifndef GL_ARB_clear_buffer_object
#define GL_ARB_clear_buffer_object 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glClearBufferData (GLenum target, GLenum internalformat, GLenum format, GLenum type, const void *data);
GLAPI void APIENTRY glClearBufferSubData (GLenum target, GLenum internalformat, GLintptr offset, GLsizeiptr size, GLenum format, GLenum type, const void *data);
GLAPI void APIENTRY glClearNamedBufferDataEXT (GLuint buffer, GLenum internalformat, GLenum format, GLenum type, const void *data);
GLAPI void APIENTRY glClearNamedBufferSubDataEXT (GLuint buffer, GLenum internalformat, GLenum format, GLenum type, GLsizeiptr offset, GLsizeiptr size, const void *data);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLCLEARBUFFERDATAPROC) (GLenum target, GLenum internalformat, GLenum format, GLenum type, const void *data);
typedef void (APIENTRYP PFNGLCLEARBUFFERSUBDATAPROC) (GLenum target, GLenum internalformat, GLintptr offset, GLsizeiptr size, GLenum format, GLenum type, const void *data);
typedef void (APIENTRYP PFNGLCLEARNAMEDBUFFERDATAEXTPROC) (GLuint buffer, GLenum internalformat, GLenum format, GLenum type, const void *data);
typedef void (APIENTRYP PFNGLCLEARNAMEDBUFFERSUBDATAEXTPROC) (GLuint buffer, GLenum internalformat, GLenum format, GLenum type, GLsizeiptr offset, GLsizeiptr size, const void *data);
#endif

#ifndef GL_ARB_compute_shader
#define GL_ARB_compute_shader 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glDispatchCompute (GLuint num_groups_x, GLuint num_groups_y, GLuint num_groups_z);
GLAPI void APIENTRY glDispatchComputeIndirect (GLintptr indirect);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLDISPATCHCOMPUTEPROC) (GLuint num_groups_x, GLuint num_groups_y, GLuint num_groups_z);
typedef void (APIENTRYP PFNGLDISPATCHCOMPUTEINDIRECTPROC) (GLintptr indirect);
#endif

#ifndef GL_ARB_copy_image
#define GL_ARB_copy_image 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glCopyImageSubData (GLuint srcName, GLenum srcTarget, GLint srcLevel, GLint srcX, GLint srcY, GLint srcZ, GLuint dstName, GLenum dstTarget, GLint dstLevel, GLint dstX, GLint dstY, GLint dstZ, GLsizei srcWidth, GLsizei srcHeight, GLsizei srcDepth);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLCOPYIMAGESUBDATAPROC) (GLuint srcName, GLenum srcTarget, GLint srcLevel, GLint srcX, GLint srcY, GLint srcZ, GLuint dstName, GLenum dstTarget, GLint dstLevel, GLint dstX, GLint dstY, GLint dstZ, GLsizei srcWidth, GLsizei srcHeight, GLsizei srcDepth);
#endif

#ifndef GL_ARB_texture_view
#define GL_ARB_texture_view 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glTextureView (GLuint texture, GLenum target, GLuint origtexture, GLenum internalformat, GLuint minlevel, GLuint numlevels, GLuint minlayer, GLuint numlayers);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLTEXTUREVIEWPROC) (GLuint texture, GLenum target, GLuint origtexture, GLenum internalformat, GLuint minlevel, GLuint numlevels, GLuint minlayer, GLuint numlayers);
#endif

#ifndef GL_ARB_vertex_attrib_binding
#define GL_ARB_vertex_attrib_binding 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glBindVertexBuffer (GLuint bindingindex, GLuint buffer, GLintptr offset, GLsizei stride);
GLAPI void APIENTRY glVertexAttribFormat (GLuint attribindex, GLint size, GLenum type, GLboolean normalized, GLuint relativeoffset);
GLAPI void APIENTRY glVertexAttribIFormat (GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset);
GLAPI void APIENTRY glVertexAttribLFormat (GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset);
GLAPI void APIENTRY glVertexAttribBinding (GLuint attribindex, GLuint bindingindex);
GLAPI void APIENTRY glVertexBindingDivisor (GLuint bindingindex, GLuint divisor);
GLAPI void APIENTRY glVertexArrayBindVertexBufferEXT (GLuint vaobj, GLuint bindingindex, GLuint buffer, GLintptr offset, GLsizei stride);
GLAPI void APIENTRY glVertexArrayVertexAttribFormatEXT (GLuint vaobj, GLuint attribindex, GLint size, GLenum type, GLboolean normalized, GLuint relativeoffset);
GLAPI void APIENTRY glVertexArrayVertexAttribIFormatEXT (GLuint vaobj, GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset);
GLAPI void APIENTRY glVertexArrayVertexAttribLFormatEXT (GLuint vaobj, GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset);
GLAPI void APIENTRY glVertexArrayVertexAttribBindingEXT (GLuint vaobj, GLuint attribindex, GLuint bindingindex);
GLAPI void APIENTRY glVertexArrayVertexBindingDivisorEXT (GLuint vaobj, GLuint bindingindex, GLuint divisor);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLBINDVERTEXBUFFERPROC) (GLuint bindingindex, GLuint buffer, GLintptr offset, GLsizei stride);
typedef void (APIENTRYP PFNGLVERTEXATTRIBFORMATPROC) (GLuint attribindex, GLint size, GLenum type, GLboolean normalized, GLuint relativeoffset);
typedef void (APIENTRYP PFNGLVERTEXATTRIBIFORMATPROC) (GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset);
typedef void (APIENTRYP PFNGLVERTEXATTRIBLFORMATPROC) (GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset);
typedef void (APIENTRYP PFNGLVERTEXATTRIBBINDINGPROC) (GLuint attribindex, GLuint bindingindex);
typedef void (APIENTRYP PFNGLVERTEXBINDINGDIVISORPROC) (GLuint bindingindex, GLuint divisor);
typedef void (APIENTRYP PFNGLVERTEXARRAYBINDVERTEXBUFFEREXTPROC) (GLuint vaobj, GLuint bindingindex, GLuint buffer, GLintptr offset, GLsizei stride);
typedef void (APIENTRYP PFNGLVERTEXARRAYVERTEXATTRIBFORMATEXTPROC) (GLuint vaobj, GLuint attribindex, GLint size, GLenum type, GLboolean normalized, GLuint relativeoffset);
typedef void (APIENTRYP PFNGLVERTEXARRAYVERTEXATTRIBIFORMATEXTPROC) (GLuint vaobj, GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset);
typedef void (APIENTRYP PFNGLVERTEXARRAYVERTEXATTRIBLFORMATEXTPROC) (GLuint vaobj, GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset);
typedef void (APIENTRYP PFNGLVERTEXARRAYVERTEXATTRIBBINDINGEXTPROC) (GLuint vaobj, GLuint attribindex, GLuint bindingindex);
typedef void (APIENTRYP PFNGLVERTEXARRAYVERTEXBINDINGDIVISOREXTPROC) (GLuint vaobj, GLuint bindingindex, GLuint divisor);
#endif

#ifndef GL_ARB_robustness_isolation
#define GL_ARB_robustness_isolation 1
#endif

#ifndef GL_ARB_ES3_compatibility
#define GL_ARB_ES3_compatibility 1
#endif

#ifndef GL_ARB_explicit_uniform_location
#define GL_ARB_explicit_uniform_location 1
#endif

#ifndef GL_ARB_fragment_layer_viewport
#define GL_ARB_fragment_layer_viewport 1
#endif

#ifndef GL_ARB_framebuffer_no_attachments
#define GL_ARB_framebuffer_no_attachments 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glFramebufferParameteri (GLenum target, GLenum pname, GLint param);
GLAPI void APIENTRY glGetFramebufferParameteriv (GLenum target, GLenum pname, GLint *params);
GLAPI void APIENTRY glNamedFramebufferParameteriEXT (GLuint framebuffer, GLenum pname, GLint param);
GLAPI void APIENTRY glGetNamedFramebufferParameterivEXT (GLuint framebuffer, GLenum pname, GLint *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLFRAMEBUFFERPARAMETERIPROC) (GLenum target, GLenum pname, GLint param);
typedef void (APIENTRYP PFNGLGETFRAMEBUFFERPARAMETERIVPROC) (GLenum target, GLenum pname, GLint *params);
typedef void (APIENTRYP PFNGLNAMEDFRAMEBUFFERPARAMETERIEXTPROC) (GLuint framebuffer, GLenum pname, GLint param);
typedef void (APIENTRYP PFNGLGETNAMEDFRAMEBUFFERPARAMETERIVEXTPROC) (GLuint framebuffer, GLenum pname, GLint *params);
#endif

#ifndef GL_ARB_internalformat_query2
#define GL_ARB_internalformat_query2 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glGetInternalformati64v (GLenum target, GLenum internalformat, GLenum pname, GLsizei bufSize, GLint64 *params);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLGETINTERNALFORMATI64VPROC) (GLenum target, GLenum internalformat, GLenum pname, GLsizei bufSize, GLint64 *params);
#endif

#ifndef GL_ARB_invalidate_subdata
#define GL_ARB_invalidate_subdata 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glInvalidateTexSubImage (GLuint texture, GLint level, GLint xoffset, GLint yoffset, GLint zoffset, GLsizei width, GLsizei height, GLsizei depth);
GLAPI void APIENTRY glInvalidateTexImage (GLuint texture, GLint level);
GLAPI void APIENTRY glInvalidateBufferSubData (GLuint buffer, GLintptr offset, GLsizeiptr length);
GLAPI void APIENTRY glInvalidateBufferData (GLuint buffer);
GLAPI void APIENTRY glInvalidateFramebuffer (GLenum target, GLsizei numAttachments, const GLenum *attachments);
GLAPI void APIENTRY glInvalidateSubFramebuffer (GLenum target, GLsizei numAttachments, const GLenum *attachments, GLint x, GLint y, GLsizei width, GLsizei height);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLINVALIDATETEXSUBIMAGEPROC) (GLuint texture, GLint level, GLint xoffset, GLint yoffset, GLint zoffset, GLsizei width, GLsizei height, GLsizei depth);
typedef void (APIENTRYP PFNGLINVALIDATETEXIMAGEPROC) (GLuint texture, GLint level);
typedef void (APIENTRYP PFNGLINVALIDATEBUFFERSUBDATAPROC) (GLuint buffer, GLintptr offset, GLsizeiptr length);
typedef void (APIENTRYP PFNGLINVALIDATEBUFFERDATAPROC) (GLuint buffer);
typedef void (APIENTRYP PFNGLINVALIDATEFRAMEBUFFERPROC) (GLenum target, GLsizei numAttachments, const GLenum *attachments);
typedef void (APIENTRYP PFNGLINVALIDATESUBFRAMEBUFFERPROC) (GLenum target, GLsizei numAttachments, const GLenum *attachments, GLint x, GLint y, GLsizei width, GLsizei height);
#endif

#ifndef GL_ARB_multi_draw_indirect
#define GL_ARB_multi_draw_indirect 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glMultiDrawArraysIndirect (GLenum mode, const void *indirect, GLsizei drawcount, GLsizei stride);
GLAPI void APIENTRY glMultiDrawElementsIndirect (GLenum mode, GLenum type, const void *indirect, GLsizei drawcount, GLsizei stride);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLMULTIDRAWARRAYSINDIRECTPROC) (GLenum mode, const void *indirect, GLsizei drawcount, GLsizei stride);
typedef void (APIENTRYP PFNGLMULTIDRAWELEMENTSINDIRECTPROC) (GLenum mode, GLenum type, const void *indirect, GLsizei drawcount, GLsizei stride);
#endif

#ifndef GL_ARB_program_interface_query
#define GL_ARB_program_interface_query 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glGetProgramInterfaceiv (GLuint program, GLenum programInterface, GLenum pname, GLint *params);
GLAPI GLuint APIENTRY glGetProgramResourceIndex (GLuint program, GLenum programInterface, const GLchar *name);
GLAPI void APIENTRY glGetProgramResourceName (GLuint program, GLenum programInterface, GLuint index, GLsizei bufSize, GLsizei *length, GLchar *name);
GLAPI void APIENTRY glGetProgramResourceiv (GLuint program, GLenum programInterface, GLuint index, GLsizei propCount, const GLenum *props, GLsizei bufSize, GLsizei *length, GLint *params);
GLAPI GLint APIENTRY glGetProgramResourceLocation (GLuint program, GLenum programInterface, const GLchar *name);
GLAPI GLint APIENTRY glGetProgramResourceLocationIndex (GLuint program, GLenum programInterface, const GLchar *name);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLGETPROGRAMINTERFACEIVPROC) (GLuint program, GLenum programInterface, GLenum pname, GLint *params);
typedef GLuint (APIENTRYP PFNGLGETPROGRAMRESOURCEINDEXPROC) (GLuint program, GLenum programInterface, const GLchar *name);
typedef void (APIENTRYP PFNGLGETPROGRAMRESOURCENAMEPROC) (GLuint program, GLenum programInterface, GLuint index, GLsizei bufSize, GLsizei *length, GLchar *name);
typedef void (APIENTRYP PFNGLGETPROGRAMRESOURCEIVPROC) (GLuint program, GLenum programInterface, GLuint index, GLsizei propCount, const GLenum *props, GLsizei bufSize, GLsizei *length, GLint *params);
typedef GLint (APIENTRYP PFNGLGETPROGRAMRESOURCELOCATIONPROC) (GLuint program, GLenum programInterface, const GLchar *name);
typedef GLint (APIENTRYP PFNGLGETPROGRAMRESOURCELOCATIONINDEXPROC) (GLuint program, GLenum programInterface, const GLchar *name);
#endif

#ifndef GL_ARB_robust_buffer_access_behavior
#define GL_ARB_robust_buffer_access_behavior 1
#endif

#ifndef GL_ARB_shader_image_size
#define GL_ARB_shader_image_size 1
#endif

#ifndef GL_ARB_shader_storage_buffer_object
#define GL_ARB_shader_storage_buffer_object 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glShaderStorageBlockBinding (GLuint program, GLuint storageBlockIndex, GLuint storageBlockBinding);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLSHADERSTORAGEBLOCKBINDINGPROC) (GLuint program, GLuint storageBlockIndex, GLuint storageBlockBinding);
#endif

#ifndef GL_ARB_stencil_texturing
#define GL_ARB_stencil_texturing 1
#endif

#ifndef GL_ARB_texture_buffer_range
#define GL_ARB_texture_buffer_range 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glTexBufferRange (GLenum target, GLenum internalformat, GLuint buffer, GLintptr offset, GLsizeiptr size);
GLAPI void APIENTRY glTextureBufferRangeEXT (GLuint texture, GLenum target, GLenum internalformat, GLuint buffer, GLintptr offset, GLsizeiptr size);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLTEXBUFFERRANGEPROC) (GLenum target, GLenum internalformat, GLuint buffer, GLintptr offset, GLsizeiptr size);
typedef void (APIENTRYP PFNGLTEXTUREBUFFERRANGEEXTPROC) (GLuint texture, GLenum target, GLenum internalformat, GLuint buffer, GLintptr offset, GLsizeiptr size);
#endif

#ifndef GL_ARB_texture_query_levels
#define GL_ARB_texture_query_levels 1
#endif

#ifndef GL_ARB_texture_storage_multisample
#define GL_ARB_texture_storage_multisample 1
#ifdef GL_GLEXT_PROTOTYPES
GLAPI void APIENTRY glTexStorage2DMultisample (GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations);
GLAPI void APIENTRY glTexStorage3DMultisample (GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations);
GLAPI void APIENTRY glTextureStorage2DMultisampleEXT (GLuint texture, GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations);
GLAPI void APIENTRY glTextureStorage3DMultisampleEXT (GLuint texture, GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations);
#endif /* GL_GLEXT_PROTOTYPES */
typedef void (APIENTRYP PFNGLTEXSTORAGE2DMULTISAMPLEPROC) (GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations);
typedef void (APIENTRYP PFNGLTEXSTORAGE3DMULTISAMPLEPROC) (GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations);
typedef void (APIENTRYP PFNGLTEXTURESTORAGE2DMULTISAMPLEEXTPROC) (GLuint texture, GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations);
typedef void (APIENTRYP PFNGLTEXTURESTORAGE3DMULTISAMPLEEXTPROC) (GLuint texture, GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations);
#endif
