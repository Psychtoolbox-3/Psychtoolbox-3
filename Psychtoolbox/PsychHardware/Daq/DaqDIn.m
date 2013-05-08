function data=DaqDIn(daq,NumberOfPorts, port)
% data=DaqDIn([DeviceIndex],[NumberOfPorts])
% USB-1208FS: Read digital ports. This command reads the current state of
% the DIO ports.  The return value will be the value seen at the port pins.
% data(1) is the 8-bit value read from port A.
% data(2) is the 8-bit value read from port B.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0
%       of the desired USB-1208FS box.
%
% USB-1608FS: There is only one port, so only one meaningful return value.  
% First argument is optional if you have only one DAQ device, but code will
% run faster if DeviceIndex is specified.  Second argument is optional; if
% it is not specified, software will attempt to determine if it should be 1
% or 2 based on the type of device.  Specifying it explicitly speeds up
% code about 1 to 2  milliseconds per call on a dual core MacI with 2.67
% GHz processors. -- mpr 
%
% USB-1024LS: Has three ports, probably numbered: 1 = port A, 4 = port B,
% 10 = port C (the sum of: 8 = portC low, 2 = portC high). Maybe you'll
% need to read port C separately in two calls for low- and high- part.
%
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest.
% DaqDeviceIndex, DaqDOut, DaqAIn, DaqAOut, DaqAInScan,DaqAOutScan.
%
% 4/15/05 dgp Wrote it.
% 12/18/07  mpr   added second argument and made first optional
% 1/11/08   mpr   swept through trying to improve consistency across daq
%                     functions
% 5/22/08   mk  Add (untested!) support for USB-1024LS box. 
% 5/23/08   mk  Add caching for HID device list. 

% Perform internal caching of list of HID devices in 'TheDevices'
% to speedup call:
persistent TheDevices;
if isempty(TheDevices)
    TheDevices = PsychHIDDAQS;
end

% Default reportId for 1x08FS devices is 3:
reportId = 3;
TheReport = uint8(0);

if ~nargin || isempty(daq)
  daq=DaqFind;
end

if nargin < 2 || isempty(NumberOfPorts)
  if strcmp(TheDevices(daq).product(1:6),'USB-16')
    NumberOfPorts = 1;
  else
    NumberOfPorts = 2;
  end
end

% Is this a USB-1024LS or similar?
if ismember(TheDevices(daq).productID, [ 118, 127 ])
    % Yes. Need different reportId and report:
    Is1024LS=1;
    reportId = 0;
    NumberOfPorts = 0; % Hack: Want exactly 1 byte from device...
    TheReport=uint8([0 port 0 0 0 0 0 0]);    
else
    Is1024LS=0;
end

PsychHID('ReceiveReports',daq);
PsychHID('GiveMeReports',daq);

err=PsychHID('SetReport',daq,2,reportId, TheReport); % DIn: Read digital ports

if err.n
    fprintf('DaqDIn SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end
[report,err]=PsychHID('GetReport',daq,1,reportId,NumberOfPorts+1); % Get report
if err.n
    fprintf('DaqDIn GetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

if ~Is1024LS
    % Default 1x08 devices:
    if length(report)==NumberOfPorts+1 && report(1)==reportId
        % 2 bytes for 1208FS or 1408FS, but only 1 byte for 1608FS
        data=double(report(2:end));
    else
        data=[];
    end
else
    % USB-1024LS: Only 1 byte with the port value:
    data = double(report(1));
end

PsychHID('ReceiveReportsStop',daq);

return;
