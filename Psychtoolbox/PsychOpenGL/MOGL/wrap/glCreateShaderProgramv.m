function program = glCreateShaderProgramv( type, count, strings, debug ) %#ok<INUSL>
% glCreateShaderProgramv  Interface to OpenGL function glCreateShaderProgramv
%
% usage:  r = glCreateShaderProgramv( type, count, strings )
%
% You have to pass a single character string 'strings' that
% contains the ASCII text of the shaders source code. Each line in the
% string needs to be terminated by a '\n' aka ASCII code 10 aka newline.
% This string will be split up into the array of strings as expected by
% the C function glCreateShaderProgramv.
%
% The 'count' parameter is ignored.
% 'type' is the type of shader.
%
% A simple way to read a shader from a standard text file, e.g., for a
% vertex shader is as follows:
% fid=fopen('MyShader.txt', 'r');
% shadersource=fread(fid);
% fclose(fid);
% shader = glCreateShaderProgramv(GL.VERTEX_SHADER, 1, shadersource);
%
% C function:  GLuint glCreateShaderProgramv(GLenum type, GLsizei count, const GLchar* const)

% 1-Sep-2012 -- Manually written by MK, based on the official spec.

% ---protected---
% ---skip---

global GL;

if nargin~=3 && nargin~=4
    error('invalid number of arguments');
end

% Special debug flag: If provided and set to non-zero, then glShaderSource is
% not executed, but only the parser inside moglcore is exercised and the result
% of parsing is output to the Matlab command window.
if nargin < 4
    debug = 0;
end

shader = glCreateShader(type);
if shader == 0
    program = 0;
    return;
end

glShaderSource(shader, strings, debug);
glCompileShader(shader);
program = glCreateProgram();
if program
    compiled = glGetShaderiv(shader, GL.COMPILE_STATUS);
    glProgramParameteri(program, GL.PROGRAM_SEPARABLE, GL.TRUE);
    if compiled
        glAttachShader(program, shader);
        glLinkProgram(program);
        glDetachShader(program, shader); 
    end
end

glDeleteShader(shader);

return;
