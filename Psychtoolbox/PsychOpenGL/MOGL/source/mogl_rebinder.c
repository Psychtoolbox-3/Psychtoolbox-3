
/*
 * mogl_rebinder.c -- Dynamic extension rebinding.
 *
 */

#include "glew.h"

// Dynamic extension rebinding: This is in an extra file mogl_rebinder.c and defined
// with explicit C linkage. Reason: For the Octave version, we need to compile all
// files (except glew.c) as C++, not C, otherwise won't work. The C++ compiler is
// more picky about function pointer assignments and will abort for hacks like in
// mogl_rebind...  -> Move this function into separate file, compile that file as
// C file, make the compiler happy.

// Dynamic rebinding of ARB extensions to core routines:
// This is a trick to get GLSL working on current OS-X (10.4.4). MacOS-X supports the OpenGL
// shading language on all graphics cards as an ARB extension. But as OS-X only supports
// OpenGL versions < 2.0 as of now, the functionality is not available as core functions, but
// only as their ARB counterparts. e.g., glCreateProgram() is always a NULL-Ptr on OS-X, but
// glCreateProgramObjectARB() is supported with exactly the same syntax and behaviour. By
// binding glCreateProgram as glCreateProgramObjectARB, we allow users to write Matlab code
// that uses glCreateProgram -- which is cleaner code than using glCreateProgramObjectARB,
// and it saves us from parsing tons of additional redundant function definitions anc code
// generation...
// In this function, we try to detect such OS dependent quirks and try to work around them...
#ifdef PTBOCTAVE
extern "C" void mogl_rebindARBExtensionsToCore(void)
#else
void mogl_rebindARBExtensionsToCore(void)
#endif
{   
#ifndef NULL
#define NULL 0
#endif
    // Remap unsupported OpenGL 2.0 core functions for GLSL to supported ARB extension counterparts:
    if (NULL == glCreateProgram) glCreateProgram = glCreateProgramObjectARB;
    if (NULL == glCreateShader) glCreateShader = glCreateShaderObjectARB;
    if (NULL == glShaderSource) glShaderSource = glShaderSourceARB;
    if (NULL == glCompileShader) glCompileShader = glCompileShaderARB;
    if (NULL == glAttachShader) glAttachShader = glAttachObjectARB;
    if (NULL == glLinkProgram) glLinkProgram = glLinkProgramARB;
    if (NULL == glUseProgram) glUseProgram = glUseProgramObjectARB;
    if (NULL == glGetAttribLocation) glGetAttribLocation = glGetAttribLocationARB;
#ifndef WINR2007a
    // Windows + Octave or any other OS:
    if (NULL == glGetUniformLocation) glGetUniformLocation = (GLint (*)(GLint, const GLchar*)) glGetUniformLocationARB;
#else
    // Windows + Matlab R2007a or later:
    if (NULL == glGetUniformLocation) glGetUniformLocation = glGetUniformLocationARB;
#endif
    if (NULL == glGetUniformfv) glGetUniformfv = glGetUniformfvARB;
    if (NULL == glUniform1f) glUniform1f = glUniform1fARB;
    if (NULL == glUniform2f) glUniform2f = glUniform2fARB;
    if (NULL == glUniform3f) glUniform3f = glUniform3fARB;
    if (NULL == glUniform4f) glUniform4f = glUniform4fARB;
    if (NULL == glUniform1fv) glUniform1fv = glUniform1fvARB;
    if (NULL == glUniform2fv) glUniform2fv = glUniform2fvARB;
    if (NULL == glUniform3fv) glUniform3fv = glUniform3fvARB;
    if (NULL == glUniform4fv) glUniform4fv = glUniform4fvARB;
    if (NULL == glUniform1i) glUniform1i = glUniform1iARB;
    if (NULL == glUniform2i) glUniform2i = glUniform2iARB;
    if (NULL == glUniform3i) glUniform3i = glUniform3iARB;
    if (NULL == glUniform4i) glUniform4i = glUniform4iARB;
    if (NULL == glUniform1iv) glUniform1iv = glUniform1ivARB;
    if (NULL == glUniform2iv) glUniform2iv = glUniform2ivARB;
    if (NULL == glUniform3iv) glUniform3iv = glUniform3ivARB;
    if (NULL == glUniform4iv) glUniform4iv = glUniform4ivARB;
    if (NULL == glUniformMatrix2fv) glUniformMatrix2fv = glUniformMatrix2fvARB;
    if (NULL == glUniformMatrix3fv) glUniformMatrix3fv = glUniformMatrix3fvARB;
    if (NULL == glUniformMatrix4fv) glUniformMatrix4fv = glUniformMatrix4fvARB;
    if (NULL == glGetShaderiv) glGetShaderiv = glGetObjectParameterivARB;
    if (NULL == glGetProgramiv) glGetProgramiv = glGetObjectParameterivARB;
    if (NULL == glGetShaderInfoLog) glGetShaderInfoLog = glGetInfoLogARB;
    if (NULL == glGetProgramInfoLog) glGetProgramInfoLog = glGetInfoLogARB;
    if (NULL == glValidateProgram) glValidateProgram = glValidateProgramARB;
    if (NULL == glIsProgram) glIsProgram = glIsProgramARB;
    
    // ARB_vertex_buffer_object:
    if (NULL == glGenBuffers) glGenBuffers = glGenBuffersARB;
    if (NULL == glDeleteBuffers) glDeleteBuffers = glDeleteBuffersARB;
    if (NULL == glBindBuffer) glBindBuffer = glBindBufferARB;
    if (NULL == glBufferData) glBufferData = glBufferDataARB;
    if (NULL == glBufferSubData) glBufferSubData = glBufferSubDataARB;
    if (NULL == glGetBufferSubData) glGetBufferSubData = glGetBufferSubDataARB;
    if (NULL == glIsBuffer) glIsBuffer = glIsBufferARB;
    if (NULL == glMapBuffer) glMapBuffer = glMapBufferARB;
    if (NULL == glUnmapBuffer) glUnmapBuffer = glUnmapBufferARB;
    if (NULL == glGetBufferPointerv) glGetBufferPointerv = glGetBufferPointervARB;
    if (NULL == glGetBufferParameteriv) glGetBufferParameteriv = glGetBufferParameterivARB;

    // Misc other stuff to remap...
    if (NULL == glDrawRangeElements) glDrawRangeElements = glDrawRangeElementsEXT;

    return;
}
