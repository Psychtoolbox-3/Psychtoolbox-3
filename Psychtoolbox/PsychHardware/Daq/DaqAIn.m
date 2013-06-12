function v=DaqAIn(daq,channel,range,UnCal)
% v=DaqAIn(DeviceIndex,channel,range,[DoNotCalibrate])
% Analog input.
% USB-1208FS:
% "v" is the measured voltage, in Volts, or NaN if no data were received.
% "DeviceIndex" is a small integer, the array index specifying which HID
%     device in the array returned by PsychHID('Devices') is interface 0
%     of the desired Daq device.
% "channel" (0 to 15) selects any of various single-ended or differential
%     measurements.
%     "channel" Measurement
%      0        0-1 (differential)
%      1        2-3 (differential)
%      2        4-5 (differential)
%      3        6-7 (differential)
%      4        1-0 (differential)
%      5        3-2 (differential)
%      6        5-4 (differential)
%      7        7-6 (differential)
%      8          0 (single-ended)
%      9          1 (single-ended)
%     10          2 (single-ended)
%     11          3 (single-ended)
%     12          4 (single-ended)
%     13          5 (single-ended)
%     14          6 (single-ended)
%     15          7 (single-ended)
% "range" (0 to 7) sets the gain (hence voltage range):
%     for single-ended measurements (channels 8-15), range is always +/- 10 V,
%     and attempts to set ranges other than 0 usually result in saturation, so
%     setting range in this function does not do anything.  For differential
%     measurements (channels 0-7), the map between range as input to this
%     function and the output generated is:
%     0 for Gain 1x (+/-20 V),   1 for Gain 2x (+/-10 V),
%     2 for Gain 4x (+/-5 V),    3 for Gain 5x (+/-4 V),
%     4 for Gain 8x (+/-2.5 V),  5 for Gain 10x (+/-2 V),
%     6 for Gain 16x (+/-1.25 V),  7 for Gain 20x (+/-1 V).
% "DoNotCalibrate" for 1208FS has no effect
%
% USB-1608FS:
% only single-ended inputs are defined, so "channel" is restricted to integers
% from 0 to 7.  Ignore the manual which claims that gains are restricted to only
% four possible values.  Thank the linux people for discovering:
%      0 for Gain 1x (+/- 10 V),      1 for Gain 2x (+/-5 V),
%      2 for Gain 4x (+/- 2.5 V),     3 for Gain 5x (+/- 2 V),
%      4 for Gain 8x (+/- 1.25 V),    5 for Gain 10x (+/- 1V),
%      6 for Gain 16x (+/- 0.625 V),  7 for Gain 32x (+/- 0.3125 V)
%
% "DoNotCalibrate" (0 or 1) if non-zero, function does not use information
% acquired from calibration measurements (see DaqCalibrateAIn).  This should
% probably always be 0 when function is called by a user (hence it defaults to
% 0).  Flag was created for the purpose of acquiring calibration data.
% Otherwise subsequent calibration measurements would be fits of fits.
%
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest,
% DaqDeviceIndex, DaqDIn, DaqDOut, DaqAIn, DaqAOut, DaqAInScan,DaqAOutScan.
%
% 4/15/05 dgp Wrote it.
% 1/21/07 asg changed normalization value from 2^16 to 2^15 to account for 16th bit as a sign bit (not data bit)
%             and modified the "range" value help.
% 6/17/07 mk  Add proper sign handling for negative voltages.
% 12/2x/07-1/x/08  mpr   modified to work with USB-1608FS and changed some
%                           terminology; particularly did away with "sign" name
%                           conflict ("sign" is a Matlab function)
% 6/07/13 mk  Try to make it more robust: Retry on no-date received, proper
%             error handling with error abort on error instead of silent
%             failure returning NaN. Cleanup.

% Perform internal caching of list of HID devices to speedup call:
persistent devices;
if isempty(devices)
    devices = PsychHIDDAQS;
end

if nargin < 4 || isempty(UnCal)
    UnCal=0;
end

if strcmp(devices(daq).product(5:6),'16')
    Is1608=1;
    MaxChannelID = 7;
else
    Is1608=0;
    MaxChannelID = 15;
    % ignore range input for single-ended measurements
    if channel > 7
        range=0;
    end
end

% Receive and flush all stale reports from device or in receive queue:
err = PsychHID('ReceiveReports',daq);
if err.n
    error('DaqAIn ReceiveReports error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

err = PsychHID('ReceiveReportsStop',daq);
if err.n
    error('DaqAIn ReceiveReportsStop error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

[reports,err] = PsychHID('GiveMeReports',daq);
if err.n
    error('DaqAIn GiveMeReports error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

if ~ismember(channel,0:MaxChannelID)
    error('DaqAIn: "channel" must be an integer 0 to %d.',MaxChannelID);
end

if ~ismember(range,0:7)
    error('DaqAIn: "range" must be an integer 0 to 7.');
end

% Send analog measurement request to device:
err = PsychHID('SetReport',daq,2,16,uint8([16 channel range]));
if err.n
    error('DaqAIn SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

% Enable reception of data from device and retry reception until a valid
% non-empty report of 3 Bytes size arrives:
report = [];
while isempty(report) || length(report)~=3
    err = PsychHID('ReceiveReports',daq);
    if err.n
        error('DaqAIn ReceiveReports-II error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
    end
    [report,err] = PsychHID('GetReport',daq,1,16,3); % Get report
    if err.n
        error('DaqAIn GetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
    end
end

% Ok, got a non-empty report of length 3 as expected. Stop reception:
err = PsychHID('ReceiveReportsStop',daq);
if err.n
    error('DaqAIn ReceiveReportsStop-II error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

if Is1608
    vmax = [10 5 2.5 2 1.25 1 0.625 0.3125];

    RawValue = double(report(3))*256 + double(report(2));
    v = vmax(range+1)*(RawValue/32768-1);
    if UnCal
        return;
    else
        DaqPrefsDir = DaqtoolboxConfigDir;
        if exist([DaqPrefsDir filesep 'DaqPrefs.mat'],'file')
            DaqVars = load([DaqPrefsDir filesep 'DaqPrefs']);
            if isfield(DaqVars,'CalData')
                CalData = DaqVars.CalData;
                GoodIndices = find(CalData(:,1) == channel & CalData(:,2) == range);
                if ~isempty(GoodIndices)
                    TheDays = CalData(GoodIndices,end);
                    ThisDay = datenum(date);
                    [DaysSinceLast,BestIndex] = min(ThisDay-TheDays);
                    AllThatDay = find(TheDays == TheDays(BestIndex));
                    MostRecentPolyFit = CalData(GoodIndices(AllThatDay(end)),3:5);

                    if DaysSinceLast > 30
                        warning('Psychtoolbox:Daq:outdatedCalibration', 'Calibration of this channel has not been performed since %s!!', ...
                            datestr(TheDays(BestIndex)));
                    end

                    v = polyval(MostRecentPolyFit,v);
                    return;
                end
            end % if isfield(DaqVars,'CalData')
        end % if exist([DaqPrefsDir filesep 'DaqPrefs.mat'],'file')
    end % if UnCal

    warning('Psychtoolbox:Daq:missingCalibration', ['It looks like this channel has not yet been calibrated. In my\n' ...
        'tests, uncalibrated values could be off by as much as 15%%!']);
else
    % Mapping table value -> voltage for differential gains:
    vmax=[20 10 5 4 2.5 2 1.25 1];

    RawReturn = double(report(2:3))*[1; 256];
    if channel < 8
        % combined two-bytes of report make a 2's complement 12-bit value
        DigitalValue = bitshift(RawReturn,-4);
        if bitget(DigitalValue,12)
            DigitalValue = -bitcmp(DigitalValue,12)-1;
        end
    else
        % range needs to be zero above during call to PsychHID('SetReport',... but
        % must be 1 to get scale below.
        range=1;
        % combined two-bytes of report make a 2's complement 11-bit value
        if RawReturn > 32752
            DigitalValue = -2048;
        elseif RawReturn > 32736
            DigitalValue = 2047;
        else
            DigitalValue = bitand(4095,bitshift(RawReturn,-3))-2048;
        end
    end

    v=vmax(range+1)*DigitalValue/2047;
end

return;
