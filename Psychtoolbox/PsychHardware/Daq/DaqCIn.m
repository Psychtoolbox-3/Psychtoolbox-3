function count=DaqCIn(daq)
% count=DaqCIn(DeviceIndex)
% USB-1208FS: Read counter.
% This function reads the 32-bit event counter in the device. This counter
% tallies the transitions of an external input attached to the CTR pin on
% the screw terminal of the device.
% "count" is the 32-bit value.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0
%       of the desired USB-1208FS box.
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest.
%
% 4/15/05 dgp Wrote it.
% 12/20/07  mpr   tested it on a 1608FS and made input argument optional
% 1/13/08   mpr   swept through and tried to make terminology consistent
%                     with that of other daq functions
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

% Default reportId for 1x08FS devices is 33:
reportId = 33;
TheReport = uint8(33);
replength = 5;

% Is this a USB-1024LS or similar?
if ismember(TheDevices(daq).productID, [ 118, 127 ])
    % Yes. Need different reportId and report:
    reportId = 0;
    TheReport=uint8([4 0 0 0 0 0 0 0]);
    replength = 4;
    Is1024LS=1;
else
    Is1024LS=0;
end

err=PsychHID('ReceiveReports',daq);
err=PsychHID('ReceiveReportsStop',daq);
[reports,err]=PsychHID('GiveMeReports',daq);
err=PsychHID('SetReport',daq,2,reportId,TheReport); % CIn
if err.n
    fprintf('DaqCIn SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end
[report,err]=PsychHID('GetReport',daq,1,reportId,replength); % get report
if err.n
    fprintf('DaqCIn GetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end

if ~Is1024LS
    if length(report)==5
        % 4 bytes
        report=double(report);
        count=report(2)+256*(report(3)+256*(report(4)+256*report(5)));
    else
        count=[];
    end
else
    if length(report)==4
        % 4 bytes
        report=double(report);
        count=report(1)+256*(report(2)+256*(report(3)+256*report(4)));
    else
        count=[];
    end
end

err=PsychHID('ReceiveReportsStop',daq);

return;
