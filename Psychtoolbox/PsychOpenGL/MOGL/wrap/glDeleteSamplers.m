function glDeleteSamplers( count, samplers )

% glDeleteSamplers  Interface to OpenGL function glDeleteSamplers
%
% usage:  glDeleteSamplers( count, samplers )
%
% C function:  void glDeleteSamplers(GLsizei count, const GLuint* samplers)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glDeleteSamplers', count, uint32(samplers) );

return
