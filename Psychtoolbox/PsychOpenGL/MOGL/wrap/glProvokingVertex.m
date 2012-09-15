function glProvokingVertex( mode )

% glProvokingVertex  Interface to OpenGL function glProvokingVertex
%
% usage:  glProvokingVertex( mode )
%
% C function:  void glProvokingVertex(GLenum mode)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glProvokingVertex', mode );

return
