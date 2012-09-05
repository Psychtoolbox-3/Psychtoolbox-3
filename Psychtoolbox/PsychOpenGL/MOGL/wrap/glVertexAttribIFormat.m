function glVertexAttribIFormat( attribindex, size, type, relativeoffset )

% glVertexAttribIFormat  Interface to OpenGL function glVertexAttribIFormat
%
% usage:  glVertexAttribIFormat( attribindex, size, type, relativeoffset )
%
% C function:  void glVertexAttribIFormat(GLuint attribindex, GLint size, GLenum type, GLuint relativeoffset)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribIFormat', attribindex, size, type, relativeoffset );

return
