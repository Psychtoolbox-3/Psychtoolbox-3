function OMLBasicTest(screenid)
% OMLBasicTest([screenid=max]) - Test basic correctness of OpenML timestamping.
%
% Performs a sequence of 300 flips, acquires timestamps of
% Flip completion according to OML_sync_control timestamping,
% as well as a post-swap GetSecs timestamp and a timestamp of
% the last vblank. If OpenML timestamping works correctly, then
% all three timestamps should be almost identical, ie., the vblank
% and flip timestamps should be identical, the GetSecs timestamp
% only minimally later than the flip timestamp, depending on
% system scheduling noise and load.
%
% If the timestamps (minus a few outliers) disagree, then
% something is likely broken in OpenML flip timestamping.
%

global t;
global drt;
global dsw;

AssertOpenGL;

Screen('Preference','VBLTimestampingmode', 4);

if nargin < 1 || isempty(screenid)
    screenid=max(Screen('Screens'));
end

win = Screen('OpenWindow', screenid);
ifi = Screen('GetFlipInterval', win);

GetSecs;
WaitSecs(0);
Priority(MaxPriority(win));

n=300;
t=zeros(3,n);

for i=1:n
	t(1,i) = Screen('Flip', win);
	t(3,i) = GetSecs;
	winfo = Screen('GetWindowInfo', win);
	t(2,i) = winfo.LastVBLTime;
end

Priority(0);
sca;

drt = 1000 * (t(3,:) - t(1,:));
dsw = 1000 * (t(2,:) - t(1,:));

mrt = max(abs(drt));
msw = max(abs(dsw));

fprintf('\n\n');
fprintf('Maximum deviation between GetSecs and stimulus onset timestamp: %f msecs. --> ', mrt);
if mrt < ifi * 1000 / 2
    fprintf('Decent.\n');
else
    fprintf('Troublesome!\n');
end

fprintf('Maximum deviation between VBlank and stimulus onset timestamp: %f msecs. --> ', msw);
if msw < ifi * 1000 / 2
    fprintf('Decent.\n');
else
    fprintf('Troublesome!\n');
end

fprintf('See plots for details...\n');

close all;
plot(drt);
title('GetSecs - swapcomplete [msecs]:');

figure;
plot(dsw);
title('vblank - swapcomplete [msecs]:');

return;
