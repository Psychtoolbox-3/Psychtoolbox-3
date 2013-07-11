function HIDIntervalTest
% HIDIntervalTest -- Illustrate human interface device sampling interval.
%
% The routine queries first the keyboard, later on the mouse buttons for
% changes in key press state or mouse button press state, recording the
% timestamp of any change. This is done first for 10 seconds on the
% keyboard, then for 10 seconds on the mouse.
%
% Time difference between consecutive "state-change" samples is computed
% and plotted as a histogram.
%
% At least on USB-HID devices - keyboards and mice connected via USB - you
% should see some characteristic distribution in the histograms. You'll see
% that clusters of samples are spaced in even time intervals, typically
% separated by multiples of 8-10 msecs. This indicates that operating
% systems only update mouse and keyboard state at a sampling frequency of
% about 100 - 125 Hz and therefore measured mouse or keyboard response
% times will be quantized with at least that granularity.
%
% There are more complex interactions, especially inside keyboards, that
% will increase latency and uncertainty about the timing of key presses
% beyond the 8 msecs granularity of the USB-HID scan cycle.
%
% Be very cautious when using normal keyboards for reaction time
% measurements, or better use a dedicated response box.

close all;

Priority(MaxPriority(0));

fprintf('\n\nTesting keyboard polling interval... Type fast for 10 seconds!!\n');
ListenChar(2);
[secs, oldKeyCode] = KbReleaseWait;

tend = secs + 10;
tsample = [];
i=0;

while GetSecs < tend
    [downstate tchange keyCode] = KbCheck; %#ok<ASGLU>
    
    if ~isequal(keyCode, oldKeyCode)
        oldKeyCode = keyCode;
        i=i+1;
        tsample(i) = tchange; %#ok<AGROW>
    end
end

ListenChar(0);
Beeper;
KbReleaseWait;

if i == 0
    fprintf('\n\nYou did not press any keys as requested. Very uncooperative! Skipping.\n');
else
    nbins = ceil(max(diff(tsample)*1000));
    hist(diff(tsample)*1000, nbins);
    title('Distribution of delta between detected keyboard events - msecs:');
end

fprintf('\n\nTesting mouse polling interval... Press mouse buttons wildly for 10 seconds!!\n');
tsample = [];
i = 0;

tend = GetSecs + 10;
[x y button] = GetMouse;
oldstate = [x,y,button];

while GetSecs < tend
    [x y button] = GetMouse;
    newstate = [x,y,button];
    
    if ~isequal(newstate, oldstate)
        oldstate = newstate;
        i=i+1;
        tsample(i) = GetSecs; %#ok<AGROW>
    end
end

Beeper;

if i == 0
    fprintf('\n\nYou did not press any mouse buttons as requested. Very uncooperative! Skipping.\n');
else
    figure;
    xxx=diff(tsample)*1000;
    nbins = ceil(max(xxx));
    hist(xxx, nbins);
    title('Distribution of delta between detected mouse events - msecs:');
end

Priority(0);

return;
