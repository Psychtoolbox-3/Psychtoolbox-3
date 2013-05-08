function err=DaqSetCal(daq,on)
% err=DaqSetCal(DeviceIndex,on)
% USB-1208FS set CAL output.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0
%       of the desired USB-1208FS box.
% "on" is 1 for on (+2.5 V), and 0 for off (0 V).
% USB-1608FS:
% arguments and logic are the same except that "on" can take more values to
% produce a wider variety of results:
%
% on                    output (Volts)
% ---                   ------
%  0                       0 
%  1                     0.625
%  2                      1.25
%  3                      2.5
%  4                        5
%
% The calibration terminal is pin 17 on the USB-1608FS.
%
% See also Daq, DaqFunctions, DaqPins, ,DaqCalibrateAIn, DaqTest, PsychHIDTest.

% 4/15/05 dgp Wrote it.
% 1/10/08 mpr updated it for USB-1608FS

AllHIDDevices = PsychHIDDAQS;
if ~isempty(strfind(AllHIDDevices(daq).product,'1608'))
	AcceptableStates = 0:4;
else
	AcceptableStates = 0:1;
end
  
if ~ismember(on,AcceptableStates)
	error(sprintf('"on" must be in the range of 0 to %d.',max(AcceptableStates)));
end
err=PsychHID('SetReport',daq,2,69,uint8([69 on])); % SetCal
if err.n
    fprintf('DaqSetCal SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

return;
