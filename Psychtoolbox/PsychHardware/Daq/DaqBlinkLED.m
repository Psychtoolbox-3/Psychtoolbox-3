function err=DaqBlinkLED(daq)
% err=DaqBlinkLED(DeviceIndex)
% USB-1208FS: blink LED.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0
%       of the desired USB-1208FS box.
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest.
%
% 4/15/05 dgp Wrote it.
% 11/13/07  mpr Tested it on USB-1608FS called by Matlab 2007b from a Mac
%                   Pro running Leopard.  Worked with no changes! 
% 1/10/08   mpr worked to get improved internal consistency (changed
%                   "device" to "daq", fixed "TestDaq" and "TestPsychHid"
% 5/22/08   mk  Add (untested!) support for USB-1024LS box. 
% 5/23/08   mk  Add caching for HID device list. 

% Perform internal caching of list of HID devices in 'TheDevices'
% to speedup call:
persistent TheDevices;
if isempty(TheDevices)
    TheDevices = PsychHIDDAQS;
end

% Default reportId for 1x08FS devices is 64:
reportId = 64;
TheReport = uint8(64);

% Is this a USB-1024LS or similar?
if ismember(TheDevices(daq).productID, [ 118, 127 ])
    % Yes. Need different reportId and report:
    reportId = 0;
    TheReport=uint8([11 0 0 0 0 0 0 0]);    
end

% DaqBlinkLED
err=PsychHID('SetReport',daq,2,reportId,TheReport); % Blink LED
if err.n
    fprintf('DaqBlinkLED SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

return;
