function err=DaqAOutStop(daq)
% err=DaqAOutStop(DeviceIndex)
% USB-1208FS stop analog output scan.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0
%       of the desired USB-1208FS box.
% See also Daq, DaqFunctions, DaqPins, DaqAOutScan, DaqTest, 
% PsychHIDTest.

% 4/15/05 dgp Wrote it.
% 1/10/08 mpr glanced at it; made same changes as in DaqAOut.

AllHIDDevices = PsychHIDDAQS;
if strcmp(AllHIDDevices(daq).product(5:6),'16')
  error('It looks like you are trying to run DaqAOutStop from a USB-1608FS, but such devices have no analog outputs.');
end

err=PsychHID('SetReport',daq,2,22,uint8(22)); % AOutStop
if err.n
    fprintf('DaqAOutStop SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

return;
