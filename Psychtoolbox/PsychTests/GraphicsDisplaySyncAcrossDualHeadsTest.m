function GraphicsDisplaySyncAcrossDualHeadsTest(screenids, nrtrials, driftsync)
% GraphicsDisplaySyncAcrossDualHeadsTest([screenids] [, nrtrials=6000] [, driftsync=0])
%
% Test synchronizity between the scanout/refresh cycles of two different
% display heads 'screenids(1)' and 'screenids(2)'. If 'screenids' is
% omitted, we test the two heads with maximal id. 'nrtrial' sample passes
% with a sampling interval of roughly 1 msecs are conducted. 'nrtrials'
% defaults to 6000 samples.
%
% Each sample consists of querying the current rasterbeam position and
% timestamp of last VBL interrupt on each display head. At the end, the
% samples of both heads are plotted against each other and compared.
%
% The help text in the Matlab prompt tells you how to interpret the plots.

% History:
% 12/11/07 Written (MK).

AssertOpenGL;

if nargin < 1
    screenids = [];
end

if isempty(screenids)
    screenids = Screen('Screens');
    if length(screenids)<2
        % Only single head - Nothing to test?!?
        error('Only one display head attached to graphics hardware. How am i supposed to test sync across heads?!?');
    end
    
    if length(screenids)>2
        % If more than two heads available, choose the two heads with
        % maximum id:
        screenids = screenids(end-1:end);
    end
end

if length(screenids) > 2
    error('"screenids" parameter contains indices of more than two screens?!? This is a dual-head sync test!');
end

if nargin < 2
    nrtrials = [];
end

if isempty(nrtrials)
    nrtrials = 6000;
end

if nargin < 3
    driftsync = [];
end

if isempty(driftsync)
    driftsync = 0;
end

driftsync = driftsync * 2;

testvblirqs = 1;

% Open windows:
w(1) = Screen('OpenWindow', screenids(1), 0, []);
w(2) = Screen('OpenWindow', screenids(2), 0, []);

if driftsync > 0
    WaitSecs(2);

    % Query current display settings for 1st display:
    oldResolution(1)=Screen('Resolution', screenids(1));

    % Switch 1st display to 60 Hz refresh rate:
    if oldResolution(1).hz ~=0
        Screen('Resolution', screenids(1), oldResolution(1).width, oldResolution(1).height, 60, oldResolution(1).pixelSize, 1);
    end
    
    % Query current display settings for 2nd display:
    oldResolution(2)=Screen('Resolution', screenids(2));

    % Switch 2nd display to 75 Hz refresh rate:
    Screen('Resolution', screenids(2), oldResolution(2).width, oldResolution(2).height, 75, oldResolution(2).pixelSize, 1);
    
    WaitSecs(2);
end

% Some beauty...
Screen('Flip', w(1), 0, 2);
Screen('Flip', w(2), 0, 2);
Screen('FillRect', w(1), 255);
Screen('FillRect', w(2), 255);
Screen('Flip', w(1), 0, 2);
Screen('Flip', w(2), 0, 2);

if driftsync > 0
    nrtrials = 100000;
end

timinginfo = zeros(2, 2, nrtrials);
syncedattrial = -1;

% Test loop one:
for i=1:nrtrials
    % Sample in roughly 1 msec intervals:
    WaitSecs(0.001);
    
    % Sample in close sync:
    for j=1:2
        winfo(j) = Screen('GetWindowInfo', w(j)); %#ok<AGROW>
    end
    
    for j=1:2
        timinginfo(1, j,i) = winfo(j).Beamposition;
        timinginfo(2, j,i) = winfo(j).LastVBLTime;
    end
    
    if driftsync > 1
        if abs(winfo(1).LastVBLTime - winfo(2).LastVBLTime) < 0.0005
            % Switch display 2 also to 60 Hz refresh:
            Screen('Resolution', screenids(2), oldResolution(2).width, oldResolution(2).height, 60, oldResolution(2).pixelSize, 1);
            driftsync = 1;
            syncedattrial = i %#ok<NOPRT>
            drawnow;
            Beeper;
        end
    end
    
    if (driftsync == 1) && (i > syncedattrial + 6000)
        stoppedattrial = i %#ok<NASGU,NOPRT>
        drawnow;
        break;
    end
end

timinginfo = timinginfo(:,:,1:i);
nrtrials = i;

%Beeper(800);
Screen('CloseAll');
%Beeper(1200);

if driftsync
    % Restore original display settings:
    Screen('Resolution', screenids(1), oldResolution(1).width, oldResolution(1).height, oldResolution(1).hz, oldResolution(1).pixelSize, 0);
    Screen('Resolution', screenids(2), oldResolution(2).width, oldResolution(2).height, oldResolution(2).hz, oldResolution(2).pixelSize, 0);
end

fprintf('What you should see in the following plots for well synchronized displays:\n');
fprintf('The "drift" plots should either show a nice zero flat line or a flat line with\n');
fprintf('a small constant offset. Of course a bit of kind of gaussion noisyness has to be\n');
fprintf('expected in the system. What you should not see is some clear trend in one direction\n');
fprintf('or a sawtooth-like pattern. That would indicate drift between the display heads due\n');
fprintf('to missing sync across the heads.\n\n');

% Some plotting:
close all;
if (range(timinginfo(1, 1, :)) > 0) && (range(timinginfo(1, 2, :)) > 0)
    plot(1:nrtrials, squeeze(timinginfo(1,1,:)), 'r', 1:nrtrials, squeeze(timinginfo(1,2,:)), 'b');
    title('Beampositions for both displays: Red = Head 1, Blue=Head 2');

    figure;
    plot(1:nrtrials, squeeze(timinginfo(1,1,:)) - squeeze(timinginfo(1,2,:)));
    title('Beamposition drift between displays:');
end

if testvblirqs
    figure;
    plot(1:nrtrials, squeeze(timinginfo(2,1,:)), 'r', 1:nrtrials, squeeze(timinginfo(2,2,:)), 'b');
    title('VBL Interrupt timestamps for both displays: Red = Head 1, Blue=Head 2');

    figure;
    plot(1:nrtrials, squeeze(timinginfo(2,1,:)) - squeeze(timinginfo(2,2,:)));
    title('VBL Interrupt timestamp drift between displays:');
end

return;
