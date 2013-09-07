function ScreenDrawDots(windowPtr, varargin)
% Workaround-Wrapper for the Screen('DrawDots') function.
%
% Usage: ScreenDrawDots(windowPtr, xy [, dotdiameter=1][, dotcolor=white][, center2D][, dot_type=1]);
%
% This function is the equivalent of the Screen('DrawDots') subfunction
% for fast drawing of 2D dots. It has the same parameters as that function,
% so it can be used as a drop-in replacement in many cases.
%
% On most systems, this function will simply call Screen('DrawDots',...);
% passing it all parameters you provided, so you'll get exactly that
% functionality, minus a little bit of performance loss due to the extra
% call overhead.
%
% On severely broken systems (like the piece of junk that Apple calls Snow
% Leopard 10.6.3), where the iPhone company managed to screw up their OpenGL
% implementation so badly that even Screen('DrawDots') doesn't work, it
% will try to use a different method to emulate 'DrawDots' behaviour as
% closely as possible. This is a hack! It will incur significant
% performance loss and it will provide an imperfect emulation which
% should work ok for simple cases of usage of 'DrawDots', but may not work
% perfectly (or at all) for complex usage or with high precision framebuffers.
%
% For syntax and further help, see help of "Screen DrawDots?"
%
% To replace direct calls to Screen('DrawDots',...); with calls to this
% function, simply use text search & replace of your text editor to search
% for:
%       Screen('DrawDots'
%
% and set the replace text to:
%
%       ScreenDrawDots(
%
%
% You can enable the workaround (across sessions) via a call:
% clear all; ScreenDrawDots(1); on the Matlab/Octave command line.
%
% You can disable the workaround (across sessions) via a call:
% clear all; ScreenDrawDots(0); on the Matlab/Octave command line.
%

% History:
% 03.04.2010  mk  Written to compensate for severe OS/X 10.6.3 OpenGL bugs.
% 05.09.2013  mk  Use WhiteIndex() instead of hard-coded 255 color value.

persistent needWorkaround
persistent spriteTextures
persistent maxdiameter

if nargin < 1 || isempty(windowPtr)
    error('"windowPtr" window handle missing! This is required!');
end

if isempty(needWorkaround)
    
    markerfile = [PsychtoolboxConfigDir('Workarounds') 'ptbDrawDotsUseWorkaround.txt'];

    if windowPtr == 1
        % Special windowPtr '1': Create markerfile to enable workaround:
        fid = fopen(markerfile, 'wt');
        fclose(fid);
        fprintf('ScreenDrawDots: Workaround enabled for this and future sessions.\n');
        return;
    end

    if windowPtr == 0
        % Special windowPtr '0': Delete markerfile to disable workaround:
        delete(markerfile);
        fprintf('ScreenDrawDots: Workaround disabled for this and future sessions.\n');
        return;
    end

    if exist(markerfile, 'file')
        % Broken. Setup replacement:
        needWorkaround = 1;
        
        fprintf('\n\nPTB-WARNING: Found marker file %s\n', markerfile);
        fprintf('PTB-WARNING: This indicates a broken Screen(''DrawText'') implementation\n');
        fprintf('PTB-WARNING: in your operating system! Will enable a slower workaround.\n');
        fprintf('PTB-WARNING: Call this function as ScreenDrawDots(0) once, if you want to disable the workaround.\n\n');
        fprintf('PTB-WARNING: You must call "clear all" before each run of your script, or\n');
        fprintf('PTB-WARNING: the workaround will abort with an error about invalid texture handles.\n');
        fprintf('PTB-WARNING: You will also see warnings at the end of each session about at least 64\n');
        fprintf('PTB-WARNING: textures being open at the end of your script...\n\n');

        % Maximum covered diameter:
        maxdiameter = 64;

        % Onetime init: Create proper textures, one for each dot size:
        spriteTextures = BuildDotTextures(windowPtr, maxdiameter);
    else
        % Everything fine, use regular Screen('DrawDots') subfunction:
        needWorkaround = 0;
    end
end

% Workaround needed?
if ~needWorkaround
    % No: Just dispatch to Screen:
    Screen('DrawDots', windowPtr, varargin{:});
    return;
end

% Workaround is needed :-( -- Unpack call arguments:
if length(varargin) < 1
    error('Required dotposition vector "xy" is missing!');
else
    xy = varargin{1};
end

if isempty(xy)
    % xy dot position matrix is empty! Nothing to do for us:
    return;
end

% Must be a 2D matrix:
if ndims(xy)~=2 %#ok<ISMAT>
    error('"xy" dot position argument is not a 2D matrix! This is required!');
end

if size(xy,1) == 1
    xy = xy';
end

% Want single xy vector as a 3 or 4 row, 1 column vector:
if (size(xy, 1) ~= 2)
    error('"xy" dot position argument is not a 2-row matrix! This is required!');
end

% Number of dots:
nrdots = size(xy, 2);

if length(varargin) < 2 || isempty(varargin{2})
    dotdiameter = 1;
else
    dotdiameter = round(varargin{2});
end

nsizes = length(dotdiameter);
if ~isvector(dotdiameter) || (nsizes~=1 && nsizes~=nrdots)
    error('"size" argument must be a vector with same number of elements as dots to draw, or a single scalar value!');
end

if min(dotdiameter) < 1 || max(dotdiameter) > maxdiameter
    error('"size" argument contains dot sizes smaller than one or greater than %f, which is not supported!', maxdiameter);
end

if length(varargin) < 3
    dotcolor = [];
else
    dotcolor = varargin{3};
end

if ~isempty(dotcolor)
    % Want single dotcolor vector as a 1-4 row, 1 column vector:
    if (size(dotcolor, 1) == 1) && (ismember(size(dotcolor, 2), [1,3,4]))
        dotcolor = transpose(dotcolor);
    end

    ncolors = size(dotcolor, 2);
    ncolcomps = size(dotcolor, 1);
    if  ~ismember(ncolcomps, [1,3,4]) || (ncolors~=1 && ncolors~=nrdots)
        error('"dotcolor" must be a matrix with 3 or 4 rows and at least 1 column, or as many columns as dots to draw!');
    end
else
    ncolors = 0; %#ok<NASGU>
    dotcolor = WhiteIndex(windowPtr);
end

% 'center2D' argument specified?
if length(varargin) < 4
    % Default to "no center set":
    center2D = [];
else
    center2D = varargin{4};
end

if ~isempty(center2D)
    % Center valid?
    if length(center2D) ~= 2
        error('center2D argument must be a 2-component vector with [x,y] center position.');
    end

    % Yes: Add its offsets to all components of xy vector:
    xy = xy + repmat(center2D(1:2)', 1, size(xy,2));
end

% 'dot_type' argument given?
if length(varargin) < 5 || isempty(varargin{5})
    % Default to "no point smoothing set":
    dot_type = 0;
else
    dot_type = varargin{5};
end

% Point smoothing wanted?
if dot_type > 0
    % Point smoothing: Use texture mapping to draw precomputed little nice
    % anti-aliased dot images:
    
    % Map dotdiameters to texture handles for textures for dotsize:
    usetex = spriteTextures(dotdiameter);

    % Build dstRects from xy and size, srcRects is always the default []:
    % We extend the dstRects by 2 because the textures have a single-pixel
    % border on each side to avoid sampling artifacts:
    if nsizes == 1
        brects = [0, 0, dotdiameter+2, dotdiameter+2]';
    else
        brects = zeros(4, nsizes);
        for i=1:nsizes
            brects(1:4, i) = [0 , 0 , dotdiameter(i)+2 , dotdiameter(i)+2];
        end
    end

    % Position the final dstRects for each dottexture centered on the
    % target xy dot positions:
    dstRects = CenterRectOnPointd(brects', xy(1,:)', xy(2,:)')';

    % Use DrawTextures to batch-draw all dot textures at their proper
    % location with the proper texture object:
    Screen('DrawTextures', windowPtr, usetex, [], dstRects, [], 1, [], dotcolor);
else
    % Square dots - just draw rectangles:

    % Build dstRects from xy and size, srcRects is always the default []:
    if nsizes == 1
        brects = [0, 0, dotdiameter, dotdiameter]';
    else
        brects = zeros(4, nsizes);
        for i=1:nsizes
            brects(1:4, i) = [0 , 0 , dotdiameter(i) , dotdiameter(i)];
        end
    end

    % Position the final dstRects for each dottexture centered on the
    % target xy dot positions:
    dstRects = CenterRectOnPointd(brects', xy(1,:)', xy(2,:)')';
    
    % Use FillRect for square dots:
    Screen('FillRect', windowPtr, dotcolor, dstRects);
end

% Done.
return;
end

% Create textures for dots of all integral diameters from 1 to maxms:
function tex = BuildDotTextures(w, maxms)
    % Prealloc texture vector:
    tex = zeros(1, maxms);
    
    % Get white value:
    white = WhiteIndex(w);

    % Create images for different diameters and corresponding textures:
    for ms = 1:maxms
        % Two-Layer Luminance+Alpha image: 1st Layer - Luminance - is
        % always a full white white. We extend the dotimage by 2 pixels in
        % each direction to have some safety margin around the dot:
        dotimage = ones(ms+2, ms+2, 2) * white;

        % 2nd alpha layer starts with zero:
        dotimage(:,:,2) = 0;
        
        % Iterate over whole image in quarter pixel steps. Sample each
        % subpixel at its center, check if it is inside the circle or
        % outside. If inside, increment "subpixel inside count" of
        % corresponding parent-pixel (tx,ty) by 1. This basically
        % partitions each pixel into a 4-by-4 subpixel grid and computes
        % the contribution of each of the 16 subpixels (inside or outside
        % mathematical circle radius) and integrates over them to find
        % final pixel color:
        
        % cp is the mathematical center of the sampling grid:
        cp  = (ms+2)/2;
        cpf = (ms/2);
        for x=0.125:0.25:ms+1-0.125
            for y=0.125:0.25:ms+1-0.125
                d = sqrt(((x-cp).^2) + ((y-cp).^2));
                if d < cpf
                    tx = floor(x)+1;
                    ty = floor(y)+1;
                    dotimage(tx, ty, 2) = dotimage(tx, ty, 2) + 1;
                end
            end
        end

        % Normalize the score between 0 and 16 for each pixel to 0.0 - 1.0,
        % then map to alpha range 0 - white to build final alpha-layer:
        dotimage(:,:,2) = dotimage(:,:,2) / 16 * white;

        if 0
            imagesc(dotimage(:,:,2));
            drawnow;
            keyboard;
            KbStrokeWait;
        end
        
        % Build texture from it and assign it to proper index in vector:
        tex(ms) = Screen('MakeTexture', w, dotimage, [], 32);
        
        % Build next diameter's texture...
    end
    
    % Done, return vector of texture handles for all possible dot sizes.
    return;
end
