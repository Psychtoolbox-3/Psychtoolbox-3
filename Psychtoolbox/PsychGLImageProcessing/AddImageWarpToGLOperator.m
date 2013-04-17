function [warpmap, warpmapgltexid] = AddImageWarpToGLOperator(gloperator, warprectangle)
% [warpmap warpmapgltexid] = AddImageWarpToGLOperator(gloperator [, warprectangle])
%
% Add a generic image warp (geometric distortion) operation to a
% GLOperator. 

global GL;

if isempty(GL)
    InitializeMatlabOpenGL([], [], 1);
end

% Check input arguments for validity and assign defaults:
if nargin < 1 || isempty(gloperator) || Screen('WindowKind', gloperator)~=4
    error('You must provide the handle of a valid GL imaging operator ''gloperator''!');
end

% Make sure gloperator is enabled for imaging operations:
Screen('HookFunction', gloperator, 'ImagingMode', mor(kPsychNeedFastBackingStore, Screen('HookFunction', gloperator, 'ImagingMode')));

% Check if gloperator is of at least 16 bpc float precision, upgrade it to
% 16 bpc float if it isn't already:
if bitand(Screen('HookFunction', gloperator, 'ImagingMode'), mor(kPsychNeed16BPCFloat, kPsychNeed32BPCFloat)) == 0
    % Not yet set:
    Screen('HookFunction', gloperator, 'ImagingMode', mor(kPsychNeed16BPCFloat, Screen('HookFunction', gloperator, 'ImagingMode')));
end

% Size of warp map and to be warped area provided?
if nargin >= 2 && ~isempty(warprectangle)
    % Create an offscreen window of corresponding size: This will define the
    % warp-map, the 2D distortion field to apply to source images or
    % framebuffers:
    warpmap = Screen('OpenOffscreenWindow', gloperator, [0 0 0 0], warprectangle, 64, 32);

    % Assign its OpenGL texture handle to our return argument. This is not
    % needed for Screen('TransformTexture'), but it is needed for application
    % in the imaging pipeline for fullscreen stimulus post-processing:
    warpmapgltexid = Screen('GetOpenGLTexture', gloperator, warpmap);
else
    % No rectangle provided: Don't create a warpmap yet:
    warpmap = [];
    warpmapgltexid = [];
end

% Load GLSL warp shader:
warpshader = LoadGLSLProgramFromFiles('TextureWarpShader.frag.txt');

% Setup texture unit mappings for shader:
glUseProgram(warpshader);

% Bind input image to unit 0:
glUniform1i(glGetUniformLocation(warpshader, 'Image'), 0);
% Bind input image to unit 1:
glUniform1i(glGetUniformLocation(warpshader, 'WarpMap'), 1);

glUseProgram(0);

% Add it to gloperator:
AddToGLOperator(gloperator, 'TextureWarpShader', warpshader);

% Ready!
return;
