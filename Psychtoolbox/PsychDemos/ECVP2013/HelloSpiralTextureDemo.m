% SpiralTextureDemo

% Written by Peter Scarfe.

PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define grey and white
white = WhiteIndex(screenNumber);
grey = GrayIndex(screenNumber);
inc = white - grey;

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Define a simple spiral texture by defining X and Y coordinates with the
% meshgrid command, converting these to polar coordinates and finally
% defining the spiral texture
[x, y] = meshgrid(-120:1:120, -120:1:120);
[th, r] = cart2pol(x, y);
spiral = grey + inc .* cos(r / 5 + th * 5);

% Make our spiral texure into a screen texture for drawing
spiralTexture = Screen('MakeTexture', window, spiral)

% We are going to draw four textures to show how a black and white texture
% can be color modulated upon drawing.
yPos = yCenter;
xPos = linspace(screenXpixels * 0.2, screenXpixels * 0.8, 4);

% Define the destination rectangles for our spiral textures. For this demo
% these will be the same size as our actual texture, but this doesn't have
% to be the case. See: ScaleSpiralTextureDemo and CheckerboardTextureDemo.
[s1, s2] = size(x);
baseRect = [0 0 s1 s2];
dstRects = nan(4, 4);
for i = 1:4
    dstRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos);
end

% Color Modulation
colorMod = [1 1 1; 1 0 0; 0 1 0; 0 0 1]';

% Batch Draw all of the texures to screen
Screen('DrawTextures', window, spiralTexture, [], dstRects, [], [], [], colorMod);

% Flip to the screen
Screen('Flip', window);

% Wait for a key press
KbStrokeWait;

% Clear the screen
sca;
