function result = bvlAcceptCalibGuiInput(timestamp)
% result = bvlAcceptCalibGuiInput([timestamp])
%
% optional timestamp set's the new timestamp to compare against with later
% calls when the timestamp is not passed in. timestamp is in seconds.
%
% Test's whether enough time has passed since we finished calibrating to
% accept input to the Gui Widgets.  This is here to get around flaws in
% Matlab's event queue.  Mouse clicks get stored up and flood the widgets
% when the foreground window closes.

%
% HACK DETAILS:  There are bugs in the handling of the event queue in Matlab.
% Mouse events get stored up (even when handled) and pass on to the dialog
% in the background.  So while we are calibrating in fullscreen mode, any
% mouse event (like a left-click) that happens over the pixel-space of a
% button in this GUI, which is IN THE BACKGROUND, is stored up and all
% passed to the GUI Widget's when the fullscreen foreground window closes.
% This is a hack, I wish I didn't have to do this, but after a full day
% of looking, I cannot find a better solution.
% We'll create a global variable that will store the timestamp of when
% we finished making any calibration edits.  Then in our GUI Callback
% functions, only accept input if the current timestamp is far enough
% past the "last calibration adjustment".  Create a function to test
% this and wrap it's functionality.
%
% Created 2007-05-15 - cburns

global gTimeOfLastCalibAdjust

inputTimeBuffer = 2.0;

if nargin == 1
    gTimeOfLastCalibAdjust = timestamp;
    result = false;
else
    currTime = GetSecs;
    deltaTime = currTime - gTimeOfLastCalibAdjust;
    
    %fprintf('bvlAcceptCalibGuiInput:\n');
    %fprintf('\t deltaTime(%f) = currTime(%f) - gTimeOfLastCalibAdjust(%f)\n', deltaTime, currTime, gTimeOfLastCalibAdjust);

    if deltaTime > inputTimeBuffer
        result = true;
    else
        result = false;
    end
end
