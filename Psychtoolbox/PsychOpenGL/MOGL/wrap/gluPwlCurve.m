function data = gluPwlCurve( nurb, count, stride, type )

% gluPwlCurve  Interface to OpenGL function gluPwlCurve
%
% usage:  data = gluPwlCurve( nurb, count, stride, type )
%
% C function:  void gluPwlCurve(GLUnurbs* nurb, GLint count, GLfloat* data, GLint stride, GLenum type)

% 24-Jan-2006 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

if ~isa(nurb,'double')
	error('argument ''nurb'' must be a pointer coded as type double');
end

moglcore( 'gluPwlCurve', nurb, count, single(data), stride, type );

return
