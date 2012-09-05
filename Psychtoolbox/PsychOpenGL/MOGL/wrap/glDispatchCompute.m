function glDispatchCompute( num_groups_x, num_groups_y, num_groups_z )

% glDispatchCompute  Interface to OpenGL function glDispatchCompute
%
% usage:  glDispatchCompute( num_groups_x, num_groups_y, num_groups_z )
%
% C function:  void glDispatchCompute(GLuint num_groups_x, GLuint num_groups_y, GLuint num_groups_z)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glDispatchCompute', num_groups_x, num_groups_y, num_groups_z );

return
