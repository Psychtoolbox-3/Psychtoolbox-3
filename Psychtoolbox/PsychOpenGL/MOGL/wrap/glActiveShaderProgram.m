function glActiveShaderProgram( pipeline, program )

% glActiveShaderProgram  Interface to OpenGL function glActiveShaderProgram
%
% usage:  glActiveShaderProgram( pipeline, program )
%
% C function:  void glActiveShaderProgram(GLuint pipeline, GLuint program)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glActiveShaderProgram', pipeline, program );

return
