function res = CV1Test(useRTbox)
% res = CV1Test([useRTbox=0]) - A timing test script for HMDs by use of a photometer.
%
% Needs the RTBox, and a photo-diode or such, e.g., a ColorCal-II,
% connected to the TTL trigger input.
%
% Atm., we still measure a discrepancy of about 75 msecs between what
% 'Flip' reports according to our drivers timestamping, and what the
% photometer based timestamping reports. That's still not good enough, and
% i'm out of ideas on how to improve this.
%
% Onset scheduling also becomes erratic for short waitframes intervals
% between white flashes.
%

if nargin < 1 || isempty(useRTbox)
    useRTbox = 0;
end

waitframes = 90;

% Setup unified keymapping and unit color range:
PsychDefaultSetup(2);

% Select screen with highest id as Oculus output display:
screenid = max(Screen('Screens'));
%PsychDebugWindowConfiguration([],0.75);
% Open our fullscreen onscreen window with black background clear color:
PsychImaging('PrepareConfiguration');
% Setup the HMD to act as a regular "monoscopic" display monitor
% by displaying the same image to both eyes:
% PsychDebugWindowConfiguration; hmd = PsychVRHMD('AutoSetupHMD', 'Monoscopic', 'HUD=0 DebugDisplay');
hmd = PsychVRHMD('AutoSetupHMD', 'Monoscopic');
PsychOculusVR1('Verbosity', 3);
Screen('Preference','Verbosity', 3);

[win, rect] = PsychImaging('OpenWindow', screenid);
ifi = Screen('GetFlipInterval', win);

% Render one view for each eye in stereoscopic mode, in an animation loop:
res.vbl = [];
res.failFlag = [];
res.tBase = [];
res.measuredTime = [];

if useRTbox
    rtbox = PsychRTBox('Open'); %, 'COM5');

    % Query and print all box settings inside the returned struct 'boxinfo':
    res.boxinfo = PsychRTBox('BoxInfo', rtbox);
    disp(res.boxinfo);

    % Enable photo-diode and TTL trigger input of box, and only those:
    PsychRTBox('Disable', rtbox, 'all');
    PsychRTBox('Enable', rtbox, 'pulse');
    WaitSecs(1);

    % Clear receive buffers to start clean:
    PsychRTBox('Stop', rtbox);
    PsychRTBox('Clear', rtbox);
    PsychRTBox('Start', rtbox);

    if ~IsWin
        % Hack: Enable async background reads to speedup box operations:
        IOPort ('ConfigureSerialport', res.boxinfo.handle, 'BlockingBackgroundRead=1 StartBackgroundRead=1');
    end
end

while ~KbCheck
    Screen('FillRect', win, 0);
    tBase = Screen('Flip', win);
    %Screen('Flip', win);
    %tBase = WaitSecs(0.5);
    res.tBase(end+1) = tBase;
    Screen('FillRect', win, 1);
    res.vbl(end+1) = Screen('Flip', win, tBase + (waitframes - 0.5) * ifi);

    % Measure real onset time:
    if useRTbox
        % Fetch sample immediately to preserve correspondence:
        [time, event, mytstamp] = PsychRTBox('GetSecs', rtbox)
        if isempty(mytstamp)
            % Failed within expected time window. This probably due to
            % tearing artifacts or GPU malfunction. Mark it as "tearing"
            % and retry for 1 full more video refresh:
            res.failFlag(end+1) = 1;
            [time, event, mytstamp] = PsychRTBox('GetSecs', rtbox);
            if isempty(mytstamp)
                % Ok, this is fucked up. No way to recover :-(
                res.failFlag(end) = 2;
                res.measuredTime(end+1) = nan;
            else
                % Got something:
                res.measuredTime(end+1) = max(mytstamp);
            end
        else
            % Success!
            res.failFlag(end+1) = 0;
            res.measuredTime(end+1) = max(mytstamp);
        end

        if ~isempty(time)
            fprintf('DT Flip %f msecs. Box uncorrected %f msecs.\n', 1000 * (res.vbl(end) - res.tBase(end)), 1000 * (max(time) - res.tBase(end)));
        end
    else
        fprintf('DT Flip %f msecs.\n', 1000 * (res.vbl(end) - res.tBase(end)));
    end
end

% Backup save for safety:
save('OculusTimingResults.mat', 'res', '-V6');

sca;

close all;
figure;

if useRTbox
    PsychRTBox('Stop', rtbox);
    PsychRTBox('Clear', rtbox);

    if ~IsWin
        % Hack: Disable async background reads:
        fprintf('Stopping background read op on box...\n');
        IOPort ('ConfigureSerialport', res.boxinfo.handle, 'StopBackgroundRead');
        IOPort ('ConfigureSerialport', res.boxinfo.handle, 'BlockingBackgroundRead=0');
        fprintf('...done, now remapping timestamps.\n');
    end

    % Primary save for safety:
    save('OculusTimingResults.mat', 'res', '-V6');

    % Remap box timestamps to GetSecs timestamps:
    res.measuredTime = PsychRTBox('BoxsecsToGetsecs', rtbox, res.measuredTime);
    fprintf('...done, saving backup copy of data, then closing box.\n');

    plot(1:length(res.vbl), 1000 * (res.vbl - res.tBase), 1:length(res.measuredTime), 1000 * (res.measuredTime - res.tBase));

    % Close connection to box:
    PsychRTBox('CloseAll');

    dT = res.measuredTime - res.vbl;
    dT = dT(~isnan(dT));
    dT = 1000 * mean(dT);
    ifi = ifi * 1000;
    fprintf('Mean difference Measurement - Flip: %f msecs [-half refresh: %f msecs], frames %f\n', dT, dT - ifi / 2, (dT - ifi / 2) / ifi);
else
    plot(1:length(res.vbl), 1000 * (res.vbl - res.tBase));
end
title('Corrected data [msecs]:');
