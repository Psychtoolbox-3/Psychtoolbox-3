% HelloWorld - The most minimalistic example of a visual stimulus.

% Make sure Psychtoolbox is properly installed and set it up for use at
% feature level 2: Normalized 0-1 color range and unified key mapping.
% This allows to specify all color or intensity values in a normalized
% range between 0.0 (for 0% output intensity) and 1.0 (for 100% output
% intensity).
PsychDefaultSetup(2);

% Select the display screen on which to show our window:
% First we get a list of all display screens and store it in the vector
% 'screens':
screens = Screen('Screens')

% 'screens' contains a list of all available display screens, numbered
% from 0 to n. 0 usually corresponds to the first display on a setup,
% or the internal flat panel of a Laptop. We select the max()imum screen
% number as 'screenid' to show our window.
screenid = max(screens)

% Open a window on display screen 'screenid', with its default background
% color, which is white:
win = PsychImaging('Openwindow', screenid)

% Set a big text size:
Screen('TextSize', win, 48);

% Draw 'Hello World!' into the center of the screen:
DrawFormattedText(win, 'Hello World!', 'center', 'center');

% Show it on the display:
Screen('Flip', win);

% Wait for a keystroke on the keyboard:
KbStrokeWait;

% Done. Close the window and clean up.
sca;
