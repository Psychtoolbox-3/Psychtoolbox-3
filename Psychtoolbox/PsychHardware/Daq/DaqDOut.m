function err=DaqDOut(daq,port,data)
% err=DaqDOut(DeviceIndex,port,data)
% USB-1208FS: Write digital port. This command writes data to the DIO port
% bits that are configured as outputs.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0
%       of the desired USB-1208FS box.
% "port" 0 = port A, 1 = port B
% "data" 8-bit value, 0 to 255.
%
% USB-1608FS: Only has one DIO port, so second argument is irrelevant and will 
% be ignored if passed.  Function will check to determine if device is a 1608,  
% and act accordingly.
%
% USB-1024LS: Has three ports, probably numbered: 1 = port A, 4 = port B,
% 10 = port C (the sum of: 8 = portC low, 2 = portC high). Maybe you'll
% need to set port C separately in two calls for low- and high- part.
%
%
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest.
% DaqDeviceIndex, DaqDIn, DaqDOut, DaqAIn, DaqAOut, DaqAInScan,DaqAOutScan.
%
% 4/15/05 dgp Wrote it.
% 12/18/07  mpr   updated for use with USB-1608FS
% 1/11/08   mpr   swept through trying to improve consistency across daq
%                   functions
% 5/22/08   mk  Add (untested!) support for USB-1024LS box. 
% 5/23/08   mk  Add caching for HID device list. 

% Perform internal caching of list of HID devices in 'TheDevices'
% to speedup call:
persistent TheDevices;
if isempty(TheDevices)
    TheDevices = PsychHIDDAQS;
end

% Denis(?) commented these out.  Don't know why.  For speed, I would guess.
% If you choose to comment them back in, then structure of code should
% change.  Logic should remain that two arguments are expected for 1608FS
% and three arguments for other device types.
%
% if ~ismember(port,0:1)
% 	error('"port" must be 0 or 1.');
% end
% if ~ismember(data,0:255)
% 	error('"data" must be in range 0:255.');
% end

if strcmp(TheDevices(daq).product(5:6),'16')
  if nargin == 2
    data = port;
  end
  TheReport = uint8([0 data]);
  % Default reportId for 1x08FS devices is 4:
  reportId = 4;
else
  TheReport = uint8([0 port data]);
  % Default reportId for 1x08FS devices is 4:
  reportId = 4;
end

% Is this a USB-1024LS or similar?
if ismember(TheDevices(daq).productID, [ 118, 127 ])
    % Yes. Need different reportId and report:
    reportId = 0;
    TheReport=uint8([1 port data 0 0 0 0 0]);    
end

err=PsychHID('SetReport',daq,2,reportId,TheReport); % DOut
if err.n
    fprintf('DaqDOut error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

return;
