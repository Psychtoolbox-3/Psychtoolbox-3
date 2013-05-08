function err=DaqALoadQueue(daq,channel,range)
% err=DaqALoadQueue(DeviceIndex,channel,range)
% USB-1208FS: Load the channel/gain queue. The USB-1208FS can scan an
% arbitrary sequence of analog input channels, each with an arbitrary gain
% setting. DaqALoadQueue allows you to specify that sequence. 
% "DeviceIndex" is a small integer, the array index specifying which HID
%     device in the array returned by PsychHID('Devices') is interface 0
%     of the desired USB-1208FS box.
% "channel" is a vector of length 1 to 8; each value (0 to 15) selects any of
%     various single-ended or differential measurements.
%  "channel"    Measurement
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
% "range" is a vector of the same length, with values of 0 to 7, specifying
%     the desired gain (and voltage range) for the corresponding channel...
%     unless we're talking channels greater than 7.  Single-ended inputs always
%     have a range of +/- 10 V, so if you try to set the range to anything other
%     than zero for a single-ended measurement, your setting will be ignored.
%     For differential measurements (channels 0:7), the mapping between the
%     values passed to this function and the gain (and actual range) are:
%     0 for 1x (+/-20 V),    1 for 2x (+/-10 V),
%     2 for 4x (+/-5 V),     3 for 5x (+/-4 V),
%     4 for 8x (+/-2.5 V),   5 for 10x (+/-2 V),
%     6 for 16x (+/-1.25 V), 7 for 20x (+/-1 V).
%
% USB-1608FS: 
% "DeviceIndex" has the same meaning here as for the USB-1208FS.
% "channel" is a bit different...
% Only single ended inputs are defined for the 1608FS, so channels range 
% only from 0 through 7.  Also, when the gain queue is loaded, it is 
% supposed to be loaded for all 8 channels -- AInScan does not allow you to
% use the gain queue as a way to determine which channels to check; you can
% *NOT* scan an arbitrary sequence of channels in a 1608 box like you can
% with a 1208 (and I believe 1408) box.  There also appears to be no way to
% poll any of the device to find out the current state of the gain queue.  
% All of this means choices have to be made.  To keep this function as 
% similar as possible for operations of both types of device, user is still
% allowed to pass an arbitrary set of channel/gain pairs (as long as there 
% no repeats), and the function will do its best not to touch the gains of 
% other channels.  To be safe, you should always specify the gains for all 
% 8 channels, but I figure you will only read from the channels whose gains
% you explicitly set...  Even though that behavior does not happen 
% automatically through calls to this function as it may with the 1208FS.
% However, since we cannot poll the device to find out the current gain
% queue settings, we're going to try to use a Preferences file.  I suspect
% this will generally work, but if multiple accounts use the device or if
% the device is sometimes attached to a different machine, then this
% function may cause some gains to get changed even if not specified by the
% user.  Did I mention that you should always specify the gains for all 8
% channels?  Anyhoo...
% "range" has the same conceptual meaning for the 1608 as it does for the 
% 1208, but the numbers are different:
%
%     0 for 1x (+/- 10 V),      1 for 2x (+/- 5 V),
%     2 for 4x (+/- 2.5 V),     3 for 5x (+/- 2 V),
%     4 for 8x (+/- 1.25 V),    5 for 10x (+/- 1 V),
%     6 for 16x (+/- 0.625 V),  7 for 32x (+/- 0.3125 V).
% 
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHidTest, DaqAIn, 
% DaqAInScan, DaqAInScanBegin.
%
% 4/15/05 dgp Wrote it.
% 1/8/08  mpr created behavior appropriate for 1608

devices = PsychHIDDAQS;

if strcmp(devices(daq).product(5:6),'16')
  Is1608 = 1;
  MaxChannelID = 7;
else
  Is1608 = 0;
  MaxChannelID = 15;
end

if length(channel)~=length(range)
    error('Length of "channel" and "range" vectors must be equal.');
end
if length(channel)~=numel(channel) || length(range)~=numel(range)
    error('"channel" and "range" must be vectors.');
end
if any(~ismember(channel,0:MaxChannelID)) || any(~ismember(range,0:7))
    error(sprintf('Each "channel" value must be in 0:%d and each range value must be in 0:7.',MaxChannelID));
end

if Is1608
  if length(channel) ~= length(unique(channel))
    error('You cannot specify the same input channel multiple times!');
  end
  DaqPrefsDir = DaqtoolboxConfigDir;
  PrefsExist = exist([DaqPrefsDir filesep 'DaqPrefs.mat'],'file');
  if PrefsExist
    DaqVars=load([DaqPrefsDir filesep 'DaqPrefs']);
    if isfield(DaqVars,'OldGains')
      OldGains = DaqVars.OldGains;
    else
      PrefsExist=0;
    end
  end % if PrefsExist
    
  LeftOut = setdiff(0:7,channel);
  range(channel+1) = range;
  if PrefsExist
    range(LeftOut+1) = OldGains(LeftOut+1);
  else
    range(LeftOut+1) = zeros(1,length(LeftOut));
  end
  report=uint8([19 range(:)']);
  err=PsychHID('SetReport',daq,2,19,report); % ALoadQueue
  if err.n
    fprintf('DaqALoadQueue error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
  else
    DaqVars.OldGains = range;
    save([DaqPrefsDir filesep 'DaqPrefs'],'-Struct','DaqVars');
  end
else % if Is1608
  % ignore range inputs for single-ended channels
  range(find(channel > 7)) = zeros(size(find(channel > 7)));
  
  report=zeros(1,2+2*length(channel));
  report(1)=19;
  report(2)=length(channel);
  for i=1:length(channel)
      report(2*i+1)=channel(i);
      report(2*i+2)=range(i);
  end
  err=PsychHID('SetReport',daq,2,19,uint8(report)); % ALoadQueue
  if err.n
      fprintf('DaqALoadQueue error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
  end
end % if Is1608; else

return;
