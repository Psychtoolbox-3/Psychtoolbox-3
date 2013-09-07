function handle = LoadShaderFromFile(filename, shadertype, debug)
% handle = LoadShaderFromFile(filename [, shadertype] [, debug])
%
% Loads a GLSL OpenGL shader definition from textfile 'filename' and
% creates an OpenGL GLSL shader of type 'shadertype'. Returns a handle to
% the new shader. If shadertype is omitted, the type of the shader is
% derived from the filename extension: .vert == GL_VERTEX_SHADER, .frag ==
% GL_FRAGMENT_SHADER, .geom == GL_GEOMETRY_SHADER, .tesscontrol =
% GL_TESS_CONTROL_SHADER, .tesseval = GL_TESS_EVALUATION_SHADER. The
% optional 'debug' flag allows to enable debugging output.
%
% On successfull load & creation, a 'handle' to the new shader is returned.
%

% 29-Mar-2006 written by MK
% 29-Jul-2009 Bugfixes by MK: debug must default to 1, not 2!
% 06-Sep-2013 Add geometry shader and tesselation shader support. (MK)

global GL;

if isempty(GL)
    InitializeMatlabOpenGL([],[],1);
end;

% Make sure we run on a GLSL capable system.
AssertGLSL;

if nargin < 3
    debug = 1;
end;

if nargin < 1 
    filename = [];
end

if isempty(filename)
    % Set default shader file name:
    error('LoadShaderFromFile() called without any shader filename!');
end;


% Add path to our own shader directory if no path given:
if isempty(fileparts(filename))
    % No path to the shader defined as part of the filename. Add our
    % default path to our standard shader directory:
    filename = [ PsychtoolboxRoot 'PsychOpenGL/PsychGLSLShaders/' filename ];
end;

if nargin < 2
    shadertype = [];
end

if isempty(shadertype)
    % Try to autodetect shadertype from file extension:
    if ~isempty(strfind(filename, '.frag'))
        shadertype = GL.FRAGMENT_SHADER;
        if debug>0, fprintf('Building a fragment shader:'); end;
    elseif ~isempty(strfind(filename, '.vert'))
        shadertype = GL.VERTEX_SHADER;
        if debug>0, fprintf('Building a vertex shader:'); end;
    elseif ~isempty(strfind(filename, '.geom'))
        shadertype = GL.GEOMETRY_SHADER;
        if debug>0, fprintf('Building a geometry shader:'); end;
    elseif ~isempty(strfind(filename, '.tesscontrol'))
        shadertype = GL.TESS_CONTROL_SHADER;
        if debug>0, fprintf('Building a tesselation control shader:'); end;        
    elseif ~isempty(strfind(filename, '.tesseval'))
        shadertype = GL.TESS_EVALUATION_SHADER;
        if debug>0, fprintf('Building a tesselation evaluation shader:'); end;        
    else
        error('No shadertype provided and could not derive type from filename extension!');
    end;
end;

if debug>0
    fprintf('Reading shader from file %s ...\n', filename);
end;

% Read shader source code from file:
[fid errmsg]=fopen(filename, 'rt');
if fid<0
    error('Could not open shader definition file %s [%s]!', filename, errmsg);
end;
shadersrc = fread(fid);
fclose(fid);

% Create shader, assign sourcecode and compile it:
handle = glCreateShader(shadertype);
if handle <= 0
    fprintf('The handle created by glCreateShader is %i -- An invalid handle!\n', handle);
    error('LoadShaderFromFile: glCreateShader failed to create a valid shader! Something is wrong with your graphics drivers!');
end

if debug > 1
    glShaderSource(handle, shadersrc, debug);
else
    glShaderSource(handle, shadersrc);
end;

glCompileShader(handle);

% Done.
return;
