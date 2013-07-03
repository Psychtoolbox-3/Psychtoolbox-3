function gluNurbsSurface( nurb, sKnotCount, sKnots, tKnotCount, tKnots, sStride, tStride, control, sOrder, tOrder, type )

% gluNurbsSurface  Interface to OpenGL function gluNurbsSurface
%
% usage:  gluNurbsSurface( nurb, sKnotCount, sKnots, tKnotCount, tKnots, sStride, tStride, control, sOrder, tOrder, type )
%
% C function:  void gluNurbsSurface(GLUnurbs* nurb, GLint sKnotCount, GLfloat* sKnots, GLint tKnotCount, GLfloat* tKnots, GLint sStride, GLint tStride, GLfloat* control, GLint sOrder, GLint tOrder, GLenum type)

% 24-Jan-2006 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=11,
    error('invalid number of arguments');
end

if ~isa(nurb,'double')
	error('argument ''nurb'' must be a pointer coded as type double');
end

moglcore( 'gluNurbsSurface', nurb, sKnotCount, single(sKnots), tKnotCount, single(tKnots), sStride, tStride, single(control), sOrder, tOrder, type );

return
