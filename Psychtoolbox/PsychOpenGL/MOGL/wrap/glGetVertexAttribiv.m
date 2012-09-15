function params = glGetVertexAttribiv( index, pname )

% glGetVertexAttribiv  Interface to OpenGL function glGetVertexAttribiv
%
% usage:  params = glGetVertexAttribiv( index, pname )
%
% C function:  void glGetVertexAttribiv(GLuint index, GLenum pname, GLint* params)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = int32(repmat(intmax('int32'),[ 32 1 ]));
moglcore( 'glGetVertexAttribiv', index, pname, params );
params = params(find(params ~= intmax('int32'))); %#ok<FNDSB>
return
