function params = glGetnUniformfvARB( program, location, bufSize )

% glGetnUniformfvARB  Interface to OpenGL function glGetnUniformfvARB
%
% usage:  params = glGetnUniformfvARB( program, location, bufSize )
%
% C function:  void glGetnUniformfvARB(GLuint program, GLint location, GLsizei bufSize, GLfloat* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---

if nargin~=3,
    error('invalid number of arguments');
end

params = single(0);

moglcore( 'glGetnUniformfvARB', program, location, bufSize, params );

return
