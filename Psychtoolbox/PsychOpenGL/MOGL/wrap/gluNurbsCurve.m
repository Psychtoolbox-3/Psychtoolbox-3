function gluNurbsCurve( nurb, knotCount, knots, stride, control, order, type )

% gluNurbsCurve  Interface to OpenGL function gluNurbsCurve
%
% usage:  gluNurbsCurve( nurb, knotCount, knots, stride, control, order, type )
%
% C function:  void gluNurbsCurve(GLUnurbs* nurb, GLint knotCount, GLfloat* knots, GLint stride, GLfloat* control, GLint order, GLenum type)

% 24-Jan-2006 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=7,
    error('invalid number of arguments');
end

if ~isa(nurb,'double')
	error('argument ''nurb'' must be a pointer coded as type double');
end

moglcore( 'gluNurbsCurve', nurb, knotCount, single(knots), stride, single(control), order, type );

return
