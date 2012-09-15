function params = glGetVertexAttribdv( index, pname )

% glGetVertexAttribdv  Interface to OpenGL function glGetVertexAttribdv
%
% usage:  params = glGetVertexAttribdv( index, pname )
%
% C function:  void glGetVertexAttribdv(GLuint index, GLenum pname, GLdouble* params)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = double(repmat(NaN,[ 32 1 ]));
moglcore( 'glGetVertexAttribdv', index, pname, params );
params = params(find(~isnan(params))); %#ok<FNDSB>

return
