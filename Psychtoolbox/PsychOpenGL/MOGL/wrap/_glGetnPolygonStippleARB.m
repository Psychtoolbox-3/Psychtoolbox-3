function pattern = glGetnPolygonStippleARB( bufSize )

% glGetnPolygonStippleARB  Interface to OpenGL function glGetnPolygonStippleARB
%
% usage:  pattern = glGetnPolygonStippleARB( bufSize )
%
% C function:  void glGetnPolygonStippleARB(GLsizei bufSize, GLubyte* pattern)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---

if nargin~=1,
    error('invalid number of arguments');
end

pattern = uint8(0);

moglcore( 'glGetnPolygonStippleARB', bufSize, pattern );

return
