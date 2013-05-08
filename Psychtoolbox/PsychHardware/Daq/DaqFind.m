function daq = DaqFind
% Syntax: DeviceIndex = DaqFind
%
% Purpose: To allow various Daq functions to run without explicitly being told
%          the index of the daq in PsychHID's enumeration.  Works iff there
%          is one device.
%
% History: 12/10/07   mpr   consolidated calls from other code
%           1/23/08   mpr   added second attempt option if No daq found
%
% see also DaqDeviceIndex and DaqReset

persistent BeenToDaqFind;

daq = DaqDeviceIndex;

if isempty(daq)
  if isempty(BeenToDaqFind)
    TryAgain=TwoStateQuery('Did not find a device, should I clear PsychHID and try again?');
    if TryAgain > 0
      BeenToDaqFind = 1;
      clear PsychHID;
      clear PsychHIDDAQS;
      daq = DaqFind;
      if isempty(daq)
        error('Did not find daq after clearance; are you sure it is connected?');
      end  
    else
      error('No daq listed in current PsychHID enumeration.');
    end
  else
    error('No daq found');
  end
elseif length(daq) > 1
  error('More than one daq found.  You must specify one of them.');
end

return;
