function glViewportIndexedfv( index, v )

% glViewportIndexedfv  Interface to OpenGL function glViewportIndexedfv
%
% usage:  glViewportIndexedfv( index, v )
%
% C function:  void glViewportIndexedfv(GLuint index, const GLfloat* v)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glViewportIndexedfv', index, single(v) );

return
