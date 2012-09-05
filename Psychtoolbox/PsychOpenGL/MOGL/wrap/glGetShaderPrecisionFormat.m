function [ range, precision ] = glGetShaderPrecisionFormat( shadertype, precisiontype )

% glGetShaderPrecisionFormat  Interface to OpenGL function glGetShaderPrecisionFormat
%
% usage:  [ range, precision ] = glGetShaderPrecisionFormat( shadertype, precisiontype )
%
% C function:  void glGetShaderPrecisionFormat(GLenum shadertype, GLenum precisiontype, GLint* range, GLint* precision)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

range = int32(zeros(1,2));
precision = int32(0);

moglcore( 'glGetShaderPrecisionFormat', shadertype, precisiontype, range, precision );

return
