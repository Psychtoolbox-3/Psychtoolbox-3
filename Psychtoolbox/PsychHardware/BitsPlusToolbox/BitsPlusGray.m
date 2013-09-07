function BitsPlusGray(whichScreen)
% BitsPlusGray([screen=max])
%
% Put a gray ramp into the Bits++ box internal CLUT.
%
% 18/04/05  ejw  Rewritten to run with OSX version of Psychtoolbox

% Define screen:
if nargin < 1 || isempty(whichScreen)
    whichScreen=max(Screen('Screens'));
end

gray=GrayIndex(whichScreen);

% Open up the a window on the screen.
% This is done by calling the Screen function of the Psychophysics
% Toolbox. We need to open the onscreen window first, so that
% the offscreen storage allocated subsequently is properly
% matched to the onscreen frame buffer.
window = Screen('OpenWindow', whichScreen, gray);

% THE FOLLOWING STEP IS IMPORTANT.
% make sure the graphics card LUT is set to a linear ramp
% (else the encoded data will not be recognised by Bits++).
LoadIdentityClut(window);

% =================================================================
% CODE NEEDED HERE !
% "linear_lut" should be replaced here with one giving the inverse
% characteristic of the monitor.
% =================================================================
% restore the Bits++ LUT to a linear ramp
linear_lut = repmat(round(linspace(0, 2^16 -1, 256))', 1, 3);
BitsPlusSetClut(window, linear_lut);

% draw a gray background on front and back buffers to clear out any old LUTs
Screen('FillRect',window, gray);
Screen('Flip', window);
Screen('FillRect',window, gray);
Screen('Flip', window);

% Close the window.
Screen('CloseAll');
