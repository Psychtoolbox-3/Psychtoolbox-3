function AddImageUndistortionToGLOperator(gloperator, exampleImage, calibrationStructureOrFile, showCalibOutput, varargin)
% AddImageUndistortionToGLOperator(gloperator, exampleImage, calibrationStructureOrFile [, showCalibOutput=0][, varargin])
%
% Add a geometric undistortion operation to a given imaging pipeline image
% processing operator 'gloperator'. 'gloperator' is the operator to add to.
% 'exampleImage' is a texture- or offscreen window handle to a texture or
% offscreen window which has exactly the color depths and size of the input
% images you want to geometrically undistort (and scale) later on. The
% created operator will be adapted to only work correctly on input images
% of that size! 'calibrationStructureOrFile' Either the file name to a
% geometric calibration file, as created by, e.g.,
% DisplayUndistortionBezier.m or DisplayUndistortionBVL.m, or a struct with
% the neccessary information as created, e.g., by CreateDisplayWarp().
% 'showCalibOutput' optional flag: If set to non-zero value, the routine
% will plot some debug output to the console or into figure windows.
% 'varargin' may contain optional parameters that will be passed to the
% routine CreateDisplayWarp() as optional arguments (see help
% CreateDisplayWarp).
%
% This routine can be used if you have geometrically distorted images from
% some source (movie file, video capture, etc.), given as textures or
% offscreen windows. It allows to define an operator that applies some
% geometric correction to that images and returns corrected versions of
% that images.
%
% A typical way of using this:
%
% a) Create an "undistortion definition file" by use of the interactive
% display calibration routines, e.g., DisplayUndistortionBezier.m or
% DisplayUndistortionBVL.m, and save it somewhere.
%
% b) In your script, call this routine, passing the filename of that
% calibration file, and an 'exampleImage' of the size and format of the
% input images to undistort, to create a suitable gloperator.
%
% c) In your code, apply the operator to distorted images by use of the
% Screen('TransformTexture') function, e.g.:
%
% correctedImage = Screen('TransformTexture', distortedImage, gloperator);
%
% An example of this procedure can be seen in the demo
% "ImageUndistortionDemo".
%

% History:
% 7/27/8 mk Written.

% Global GL const struct:
global GL;

% Initialize GL struct if not done already by calling code:
if isempty(GL)
    % Only initialize in 2D mode, ie. noswitchTo3D == 1:
    InitializeMatlabOpenGL([], [], 1);
end

% Check input arguments for validity and assign defaults:
if nargin < 1 || isempty(gloperator) || Screen('WindowKind', gloperator)~=4
    error('You must provide the handle of a valid GL imaging operator ''gloperator''!');
end

if nargin < 2 || isempty(exampleImage)
    error('You must provide the handle to an ''exampleImage'' as 2nd argument!');
end

if nargin < 3 || isempty(calibrationStructureOrFile)
    error('You must provide the ''calibrationStructureOrFile'' parameter as 3rd argument!');
end

if nargin < 4 || isempty(showCalibOutput)
    showCalibOutput = 0;
end

% Switch to gloperators OpenGL context:
Screen('GetWindowInfo', gloperator);

% Make sure gloperator is enabled for imaging operations:
Screen('HookFunction', gloperator, 'ImagingMode', mor(kPsychNeedFastBackingStore, Screen('HookFunction', gloperator, 'ImagingMode')));

if ischar(calibrationStructureOrFile)
    % Setup calibration: 'calibrationStructureOrFile' is the path to the calibration file
    % that does undistortion and scaling: varargin{1:2} is resolution of warp-mesh.
    % Higher numbers == finer undistortion but higher draw time.
    % 'exampleImage' defines size and color depths of the input image.
    calib = CreateDisplayWarp(exampleImage, calibrationStructureOrFile, showCalibOutput, varargin{:});
else
    if ~isstruct(calibrationStructureOrFile)
        error('Parameter ''calibrationStructureOrFile'' is neither a calibration structure, nor the filename of a calibration file!');
    else
        calib = calibrationStructureOrFile;
    end
end

% Setup blitter configuration for display list blitter and optional custom
% filtering shader:
glsl = calib.glsl;
gld = calib.gld;

if ~isempty(glsl)
    % Shader assigned for custom filtering:
    blittercfg = sprintf('Blitter:DisplayListBlit:Handle:%i', gld);
else
    % No shader: Select standard bilinear filtering...
    blittercfg = sprintf('Blitter:DisplayListBlit:Handle:%i:Bilinear', gld);
    % ...and assign a dummy zero shader, aka the fixed-function pipeline:
    glsl = 0;
end

% Add it to gloperator:
AddToGLOperator(gloperator, 'ImageUndistortion', glsl, blittercfg);

% Ready!
return;
