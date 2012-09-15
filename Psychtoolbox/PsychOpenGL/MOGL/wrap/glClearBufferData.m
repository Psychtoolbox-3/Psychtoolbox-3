function glClearBufferData( target, internalformat, format, type, data )

% glClearBufferData  Interface to OpenGL function glClearBufferData
%
% usage:  glClearBufferData( target, internalformat, format, type, data )
%
% C function:  void glClearBufferData(GLenum target, GLenum internalformat, GLenum format, GLenum type, const void* data)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glClearBufferData', target, internalformat, format, type, data );

return
