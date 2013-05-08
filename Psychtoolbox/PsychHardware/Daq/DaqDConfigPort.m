function err=DaqDConfigPort(daq,port,direction)
% err=DaqDConfigPort(DeviceIndex,port,direction)
% USB-1208FS: Configure digital port. This command sets the direction of
% the DIO port to input or output.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0
%       of the desired USB-1208FS box.
% "port" 0 = port A, 1 = port B
% "direction" 0 = output, 1 = input
%
% USB-1608FS: has only one port, so second argument is meaningless.  If three
% arguments are passed and device is a 1608FS, second argument is ignored.
% Function could be made more efficient if a separate one were written for
% 1608FS, but I'm guessing this function will never be used at a time when 
% maximal speed is an overarching requirement.  
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHidTest.
%
% USB-1024LS: Has three ports, probably numbered: 1 = port A, 4 = port B,
% 10 = port C (the sum of: 8 = portC low, 2 = portC high).
%
% 4/15/05   dgp Wrote it.
% 12/18/07  mpr extended for 1608
% 1/11/08   mpr swept through trying to improve consistency across daq
%                 functions
% 5/22/08   mk  Add (untested!) support for USB-1024LS box. 
% 5/23/08   mk  Add caching for HID device list. 
% 8/12/2010 sdv Fixed error when other HID devices have short names

% Perform internal caching of list of HID devices in 'TheDevices'
% to speedup call:
persistent TheDevices;
if isempty(TheDevices)
    TheDevices = PsychHIDDAQS;
end

% Default reportId for 1x08FS devices is 1:
reportId = 1;

Is1608=0; % assume no 1608, change this if there is
if length(TheDevices(daq).product)>=6,  % check for sufficient length
    if strcmp(TheDevices(daq).product(5:6),'16')
      Is1608=1;
      if nargin == 2
        direction = port;
      end
    
    end
end

% Denis(?) commented this out for some reason...
% if ~ismember(port,0:1)
% 	error('"port" must be 0 or 1.');
% end
if ~ismember(direction,0:1)
	error('"direction" must be 0 (out) or 1 (in).');
end
if Is1608
  report = uint8([1 direction]);
else
  report=uint8([1 port direction]);
end

% Is this a USB-1024LS or similar?
if ismember(TheDevices(daq).productID, [ 118, 127 ])
    % Yes. Need different reportId and report:
    reportId = 0;
    report=uint8([13 port direction 0 0 0 0 0]);    
end

err=PsychHID('SetReport',daq,2,reportId,report); % Configure digital port.
if err.n
    fprintf('DaqDConfigPort error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end
return
