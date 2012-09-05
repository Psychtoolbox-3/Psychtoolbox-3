function glDrawArraysInstancedBaseInstance( mode, first, count, instancecount, baseinstance )

% glDrawArraysInstancedBaseInstance  Interface to OpenGL function glDrawArraysInstancedBaseInstance
%
% usage:  glDrawArraysInstancedBaseInstance( mode, first, count, instancecount, baseinstance )
%
% C function:  void glDrawArraysInstancedBaseInstance(GLenum mode, GLint first, GLsizei count, GLsizei instancecount, GLuint baseinstance)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glDrawArraysInstancedBaseInstance', mode, first, count, instancecount, baseinstance );

return
