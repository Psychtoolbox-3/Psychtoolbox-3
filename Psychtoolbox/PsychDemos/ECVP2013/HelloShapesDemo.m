% HelloShapesDemo - Draw some basic shapes in different colors.

% Make sure Psychtoolbox is properly installed and set it up for use at
% feature level 2: Normalized 0-1 color range and unified key mapping.
PsychDefaultSetup(2);

% Select the display screen on which to show our window:
screenid = max(Screen('Screens'));

% Open a window on display screen 'screenid'. We choose a 50% gray
% background by specifying a luminance level of 0.5:
win = PsychImaging('Openwindow', screenid, 0.5);

% Set a big text size:
Screen('TextSize', win, 48);

% Draw 'Hello World!' into the center of the screen, and in a new line,
% as marked by the '\n' new line marker, we also draw 'Press a key':
DrawFormattedText(win, 'Hello World!\nPress a key', 'center', 'center');

% Show it on the display:
Screen('Flip', win);

% Wait for a keystroke on the keyboard:
KbStrokeWait;

% Draw some basic shapes:

% Draw some framed and filled rects. We define the bounding rectangles
% [left, top, right, bottom] of the rectangles:
r1=[0 0 100 100];

% Helper functions (help PsychRects) allow to manipulate rectangles,
% e.g., shift a rect 'r1' around by an x,y offset:
r2=OffsetRect(r1, 200,200);
r3=OffsetRect(r1, 250,250);
r4=OffsetRect(r1, 300,300);

% Draw a rectangle contour in color [red, green, blue] = [1, 0, 0] at
% location r2, with a line width of 1 pixel:
Screen('FrameRect',win, [1, 0, 0], r2, 1);

% Another one in green at location r3, 4 pixels wide border:
Screen('FrameRect',win, [0, 1, 0], r3, 4);

% And in blue at location r4 with 9 pixels border:
Screen('FrameRect',win, [0, 0, 1], r4, 9);

% And a yellow one at r4, 1 pixel thin:
Screen('FrameRect',win, [1, 1, 0], r4, 1);

% And a filled rectangle at location r5, with a orange color:
r5=OffsetRect(r1, 100,300);
Screen('FillRect', win, [0.8, 0.5, 0.0], r5);

% And some oval:
r6=OffsetRect([0 0 160 80], 10,150);
Screen('FillOval', win, [0.0, 0.5, 0.8], r6);

% Draw a white (=1.0 intensity) line from (100,100) to (400, 400), 3 pixels wide:
Screen('DrawLine', win, 1, 100,100, 400,400, 3);

% Draw some filled and framed polygons:
polyCenterX=500;
polyCenterY=300;
numPoints=15;
polyRadius=50;

% Place the 'numPoints' points at random locations on a circle with
% radius 'polyRadius' and center [polyCenterX, polyCenterY]. For this
% we choose random angles on the circle and calculate the [x,y] locations
% of the points via triogonometric functions. 'polyPoints' is the list of [x,y]
% points which define the polygon:
ClockRandSeed;
angles=sort(rand(1,numPoints) * 2 * pi);
polyPoints=round([ [cos(angles)*polyRadius+polyCenterX]', [sin(angles)*polyRadius+polyCenterY]'])

% Draw a filled polygon:
Screen('FillPoly', win, [0.5, 0.5, 1], polyPoints);

% Draw polygon contour at different vertical location:
polyCenterY=150;
polyPoints=round([ [cos(angles)*polyRadius+polyCenterX]', [sin(angles)*polyRadius+polyCenterY]']);

% Draw only polygon contour:
Screen('FramePoly', win, [0.0, 0.5, 0.5], polyPoints,2);

% Draw instructions, but now at starting location [x,y] = [10, 40], in blue color:
DrawFormattedText(win, 'Hello World!\nPress a key', 10, 40, [0 0 1]);

% Show it on the display:
Screen('Flip', win);

% Wait for a keystroke on the keyboard:
KbStrokeWait;

% Done. Close the window and clean up.
sca;
