/*
 
 Copyright (c) 2010 David Petrie
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 
 */

#ifndef FTGL_ES_GLUE_H
#define FTGL_ES_GLUE_H

#ifdef __linux__
#include <GLES/gl.h>
#include <GLES/glext.h>
#else 
#include <OpenGLES/ES1/gl.h>
#include <OpenGLES/ES1/glext.h>
#endif

#include <stdio.h>
#include <assert.h>

#define GL_QUADS 888

#ifdef __cplusplus
extern "C" {
#endif
	
	extern GLvoid ftglBegin( GLenum prim );
	
	extern GLvoid ftglVertex3f( float x, float y, float z );
	
	extern GLvoid ftglVertex2f( float x, float y);
	
	extern GLvoid ftglColor4f( GLfloat r, GLfloat g, GLfloat b, GLfloat a );
		
	extern GLvoid ftglTexCoord2f( GLfloat s, GLfloat t );
	
	extern GLvoid ftglEnd();
	
	extern GLvoid ftglError(const char *source);
	
#ifdef __cplusplus
}
#endif



#endif
