function err=DaqCInit(daq)
% err=DaqCInit(DeviceIndex)
% USB-1208FS: Initialize counter to zero.
% "DeviceIndex" is a small integer, the array index specifying which HID
%         device in the array returned by PsychHID('Devices') is interface 0
%         of the desired USB-1208FS box.
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest.
%
% 4/15/05 dgp Wrote it.
% 12/20/07  mpr   tested it on 1608FS and then made input argument optional
% 1/10/08   mpr   swept through trying to improve consistency across Daq
%                     functions
% 5/22/08   mk  Add (untested!) support for USB-1024LS box. 
% 5/23/08   mk  Add caching for HID device list. 

% Perform internal caching of list of HID devices in 'TheDevices'
% to speedup call:
persistent TheDevices;
if isempty(TheDevices)
    TheDevices = PsychHIDDAQS;
end

if ~nargin || isempty(daq)
  daq=DaqFind;
end

% Default reportId for 1x08FS devices is 32:
reportId = 32;
TheReport = uint8(32);

% Is this a USB-1024LS or similar?
if ismember(TheDevices(daq).productID, [ 118, 127 ])
    % Yes. Need different reportId and report:
    reportId = 0;
    TheReport=uint8([5 0 0 0 0 0 0 0]);    
end

err=PsychHID('SetReport',daq,2,reportId,TheReport); % CInit
if err.n
    fprintf('CInit SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

return;
