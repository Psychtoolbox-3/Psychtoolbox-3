% HelloAlphaDemo -- Minimal demonstration of alpha blending.
%
% Derived from code written by Peter Scarfe.
%

PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Get value of color white:
white = WhiteIndex(screenNumber)

% Open an on screen window with a white background color:
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, white);

% Set blend function for alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Make a base Rect of 300 by 300 pixels
baseRect = [0 0 300 300];

% Screen X positions of our three rectangles:
squareXpos = [xCenter - 200 xCenter + 200 xCenter];
squareYpos = [yCenter yCenter yCenter + 200];
numSqaures = length(squareXpos);

% Set the colors to Red, Green and Blue, with the fourth value being the
% "alpha" value. This also takes a value between 0 and 255 just like a
% normal colour, however now 0 = totally transparent and 255 = totally
% opaque. Our RGB triplets are now RGBA values.
allColors = [1 0 0 1; 0 1 0 1; 0 0 1 0.5]';

% Make our rectangle coordinates
allRects = nan(4, 3);
for i = 1:numSqaures
    allRects(:, i) = CenterRectOnPointd(baseRect, squareXpos(i), squareYpos(i));
end

% Draw the rect to the screen
Screen('FillRect', window, allColors, allRects);

% Flip to the screen
Screen('Flip', window);

% Wait for a key press
KbStrokeWait;

% Clear the screen
sca;
