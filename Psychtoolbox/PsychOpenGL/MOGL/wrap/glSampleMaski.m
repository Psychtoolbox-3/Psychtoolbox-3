function glSampleMaski( index, mask )

% glSampleMaski  Interface to OpenGL function glSampleMaski
%
% usage:  glSampleMaski( index, mask )
%
% C function:  void glSampleMaski(GLuint index, GLbitfield mask)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glSampleMaski', index, mask );

return
