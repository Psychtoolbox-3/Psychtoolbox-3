function glVertexAttribFormat( attribindex, size, type, normalized, relativeoffset )

% glVertexAttribFormat  Interface to OpenGL function glVertexAttribFormat
%
% usage:  glVertexAttribFormat( attribindex, size, type, normalized, relativeoffset )
%
% C function:  void glVertexAttribFormat(GLuint attribindex, GLint size, GLenum type, GLboolean normalized, GLuint relativeoffset)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribFormat', attribindex, size, type, normalized, relativeoffset );

return
