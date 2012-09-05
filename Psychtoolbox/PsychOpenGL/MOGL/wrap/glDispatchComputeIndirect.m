function glDispatchComputeIndirect( ptr )

% glDispatchComputeIndirect  Interface to OpenGL function glDispatchComputeIndirect
%
% usage:  glDispatchComputeIndirect( ptr )
%
% C function:  void glDispatchComputeIndirect(GLint ptr)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glDispatchComputeIndirect', ptr );

return
