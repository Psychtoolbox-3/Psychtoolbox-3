function params = glGetUniformdv( program, location )

% glGetUniformdv  Interface to OpenGL function glGetUniformdv
%
% usage:  params = glGetUniformdv( program, location )
%
% C function:  void glGetUniformdv(GLuint program, GLint location, GLdouble* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---
% ---skip---

if nargin~=2,
    error('invalid number of arguments');
end

% Allocate enough space of float double datatype for a sufficient
% number of arguments (worst case): I believe the biggest return could be a
% 4 x 4 matrix, i.e. 16 elements, but i'm not sure about arrays...
% FIXME: Verify this assumption!
params = double(repmat(NaN,[ 16 1 ]));

moglcore( 'glGetUniformdv', program, location, params );
params = params(find(~isnan(params))); %#ok<FNDSB>
