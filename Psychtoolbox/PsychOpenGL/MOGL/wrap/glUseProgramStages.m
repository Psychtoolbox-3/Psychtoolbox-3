function glUseProgramStages( pipeline, stages, program )

% glUseProgramStages  Interface to OpenGL function glUseProgramStages
%
% usage:  glUseProgramStages( pipeline, stages, program )
%
% C function:  void glUseProgramStages(GLuint pipeline, GLbitfield stages, GLuint program)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glUseProgramStages', pipeline, stages, program );

return
