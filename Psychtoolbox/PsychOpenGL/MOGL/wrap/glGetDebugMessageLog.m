function [ r, sources, types, ids, severities, lengths, messageLog ] = glGetDebugMessageLog( count, bufsize )

% glGetDebugMessageLog  Interface to OpenGL function glGetDebugMessageLog
%
% usage:  [ r, sources, types, ids, severities, lengths, messageLog ] = glGetDebugMessageLog( count, bufsize )
%
% C function:  GLuint glGetDebugMessageLog(GLuint count, GLsizei bufsize, GLenum* sources, GLenum* types, GLuint* ids, GLenum* severities, GLsizei* lengths, GLchar* messageLog)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

sources = uint32(zeros(1,count));
types = uint32(zeros(1,count));
ids = uint32(zeros(1,count));
severities = uint32(zeros(1,count));
lengths = int32(zeros(1,count));
messageLog = uint8(zeros(1,bufsize));

r = moglcore( 'glGetDebugMessageLog', count, bufsize, sources, types, ids, severities, lengths, messageLog );
messageLog = char(messageLog);

return
