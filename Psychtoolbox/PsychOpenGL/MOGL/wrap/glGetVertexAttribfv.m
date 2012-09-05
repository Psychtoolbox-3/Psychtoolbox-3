function params = glGetVertexAttribfv( index, pname )

% glGetVertexAttribfv  Interface to OpenGL function glGetVertexAttribfv
%
% usage:  params = glGetVertexAttribfv( index, pname )
%
% C function:  void glGetVertexAttribfv(GLuint index, GLenum pname, GLfloat* params)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = single(repmat(NaN,[ 32 1 ]));
moglcore( 'glGetVertexAttribfv', index, pname, params );
params = double(params);
params = params(find(~isnan(params))); %#ok<FNDSB>

return
