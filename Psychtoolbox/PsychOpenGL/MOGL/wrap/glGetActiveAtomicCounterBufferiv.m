function params = glGetActiveAtomicCounterBufferiv( program, bufferIndex, pname )

% glGetActiveAtomicCounterBufferiv  Interface to OpenGL function glGetActiveAtomicCounterBufferiv
%
% usage:  params = glGetActiveAtomicCounterBufferiv( program, bufferIndex, pname )
%
% C function:  void glGetActiveAtomicCounterBufferiv(GLuint program, GLuint bufferIndex, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

global GL;

if nargin~=3,
    error('invalid number of arguments');
end

% GL.ATOMIC_COUNTER_BUFFER_ACTIVE_ATOMIC_COUNTER_INDICES has a variable
% number of return indices, which we need to query to prealloc properly:
if pname == GL.ATOMIC_COUNTER_BUFFER_ACTIVE_ATOMIC_COUNTER_INDICES
    params = int32(zeros(1, int32(glGetActiveAtomicCounterBufferiv(program, bufferIndex, GL.ATOMIC_COUNTER_BUFFER_ACTIVE_ATOMIC_COUNTERS))));
else
    % All other pname's return exactly one scalar value:
    params = int32(0);
end

moglcore( 'glGetActiveAtomicCounterBufferiv', program, bufferIndex, pname, params );

return
