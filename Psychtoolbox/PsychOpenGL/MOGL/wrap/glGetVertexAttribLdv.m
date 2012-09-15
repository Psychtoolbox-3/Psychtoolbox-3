function params = glGetVertexAttribLdv( index, pname )

% glGetVertexAttribLdv  Interface to OpenGL function glGetVertexAttribLdv
%
% usage:  params = glGetVertexAttribLdv( index, pname )
%
% C function:  void glGetVertexAttribLdv(GLuint index, GLenum pname, GLdouble* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = double(repmat(NaN,[ 32 1 ]));
moglcore( 'glGetVertexAttribLdv', index, pname, params );
params = params(find(~isnan(params))); %#ok<FNDSB>

return
