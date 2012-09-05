function params = glGetnUniformuivARB( program, location, bufSize )

% glGetnUniformuivARB  Interface to OpenGL function glGetnUniformuivARB
%
% usage:  params = glGetnUniformuivARB( program, location, bufSize )
%
% C function:  void glGetnUniformuivARB(GLuint program, GLint location, GLsizei bufSize, GLuint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---

if nargin~=3,
    error('invalid number of arguments');
end

params = uint32(0);

moglcore( 'glGetnUniformuivARB', program, location, bufSize, params );

return
