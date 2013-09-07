function PsychDrawSprites2D(windowPtr, spriteTex, xy, spriteScale, spriteAngle, spriteColor, center2D, filterMode, spriteShader)
% Draw little texture in multiple places, with a similar syntax as
% Screen('DrawDots').
%
% Usage: PsychDrawSprites2D(windowPtr, spriteTex, xy [, spriteScale=1][, spriteAngle=0][, spriteColor=white][, center2D=[0,0]][, spriteShader]);
%
% This is basically a 'DrawDots' work-alike, just that it doesn't draw nice
% round dots, but many copies of nice little texture images, each
% individually colored, scaled, positioned, shaded and rotated.
%
% It is a convenience wrapper around the Screen('DrawTextures') command, a
% more low-level, more flexible command for drawing of many textures at
% once.
%
% 'windowPtr' target window. 'spriteTex' texture to draw. 'xy' 2-row matrix
% with x,y locations of the centers of the drawn textures. 'center2D' adds
% an optional offset to these locations. 'spriteScale' global or per-sprite
% scaling factor. 'spriteAngle' global or per-sprite rotation angle.
% 'spriteColor' global or per sprite 1-, 3- or 4-row color matrix.
% 'spriteShader' a shader handle for shaded drawing.
%
% See Screen DrawTextures? for more help on the options and Screen
% DrawDots? for similar syntax. See DotDemo.m or DotDemo(1) for example
% usage.
%

% History:
% 18.04.2010  mk  Wrote it.
% 05.09.2013  mk  Use WhiteIndex() instead of hard-coded 255 color value.

if nargin < 1 || isempty(windowPtr)
    error('"windowPtr" window handle missing! This is required!');
end

if nargin < 2 || isempty(spriteTex)
    error('"spriteTex" handle missing! This is required!');
end

% Unpack call arguments:
if nargin < 3
    error('Required dotposition vector "xy" is missing!');
end

if isempty(xy)
    % xy position matrix is empty! Nothing to do for us:
    return;
end

% Must be a 2D matrix:
if ndims(xy)~=2
    error('"xy" position argument is not a 2D matrix! This is required!');
end

if size(xy,1) == 1
    xy = xy';
end

% Want xy matrix as a 2 row,  x-column matrix:
if (size(xy, 1) ~= 2)
    error('"xy" dot position argument is not a 2-row matrix! This is required!');
end

% Number of dots:
nrdots = size(xy, 2);

if nargin < 4 || isempty(spriteScale)
    spriteScale = 1;
end

nsizes = length(spriteScale);
if ~isvector(spriteScale) || (nsizes~=1 && nsizes~=nrdots)
    error('"spriteScale" argument must be a vector with same number of elements as textures to draw, or a single scalar value!');
end

if min(spriteScale) <= 0
    error('"spriteScale" argument contains negative or zero scale factors!');
end

if nargin < 5
    spriteAngle = [];
end

nangles = length(spriteAngle);
if ~isempty(spriteAngle)
    if ~isvector(spriteAngle) || (nangles~=1 && nangles~=nrdots)
        error('"spriteAngle" argument must be a vector with same number of elements as textures to draw, or a single scalar value!');
    end
end

if nargin < 6
    spriteColor = [];
end

if ~isempty(spriteColor)
    % Want single spriteColor vector as a 1-4 row, 1 column vector:
    if (size(spriteColor, 1) == 1) && (ismember(size(spriteColor, 2), [1,3,4]))
        spriteColor = transpose(spriteColor);
    end

    ncolors = size(spriteColor, 2);
    ncolcomps = size(spriteColor, 1);
    if  ~ismember(ncolcomps, [1,3,4]) || (ncolors~=1 && ncolors~=nrdots)
        error('"spriteColor" must be a matrix with 3 or 4 rows and at least 1 column, or as many columns as sprites to draw!');
    end
else
    ncolors = 0; %#ok<NASGU>
    spriteColor = WhiteIndex(windowPtr);
end

if nargin < 8
    filterMode = [];
end

if nargin < 9
    spriteShader = [];
end

nshaders = length(spriteShader);
if ~isempty(spriteShader)
    if ~isvector(spriteShader) || (nshaders~=1 && nshaders~=nrdots)
        error('"spriteShader" argument must be a vector with same number of elements as textures to draw, or a single scalar value!');
    end
end

% 'center2D' argument specified?
if nargin < 7
    % Default to "no center set":
    center2D = [];
end

if ~isempty(center2D)
    % Center valid?
    if length(center2D) ~= 2
        error('center2D argument must be a 2-component vector with [x,y] center position.');
    end

    % Yes: Add its offsets to all components of xy vector:
    xy = xy + repmat(center2D(1:2)', 1, size(xy,2));
end

% Build dstRects from xy and size, srcRects is always the default []:
srcRect = Screen('Rect', spriteTex);
brects = srcRect' * spriteScale;

% Position the final dstRects for each texture centered on the
% target xy positions:
dstRects = CenterRectOnPointd(brects', xy(1,:)', xy(2,:)')';

% Use DrawTextures to batch-draw all textures at their proper location with
% the proper texture objects and shaders:
Screen('DrawTextures', windowPtr, spriteTex, [], dstRects, spriteAngle, filterMode, [], spriteColor, spriteShader);

% Done.
return;
