function data=DaqGetAll(daq)
% data=DaqGetAll(DeviceIndex)
% USB-1208FS: Retrieve all analog and digital input values. This command
% reads the value from all analog input channels and digital I/Os.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0 
%       of the desired USB-1208FS box.
% data.analog is a 1x16 array of the values read from the analog input
%       channels.
% data.digital(1) is the value (0 to 255) read from digital port A.
% data.digital(2) is the value (0 to 255) read from digital port B.
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest.

% 4/15/05 dgp Wrote it.
% 12/27/07  mpr   added some commentary; note that so far as I know, this
%                  function does *not* work with the USB-1608FS.  I tried 
%                  several variations of different parts but gave up trying
%                  before I was able to figure out why.
% 1/13/08   mpr   swept through and tried to make terminology conistent
%                   with that of other daq functions

AllHIDDevices = PsychHIDDAQS;
if ~isempty(strfind(AllHIDDevices(daq).product,'1608'))
  Is1608=1;
else
  Is1608=0;
end

ReceiveReportsErr=PsychHID('ReceiveReports',daq);
ReceiveReportsStopErr=PsychHID('ReceiveReportsStop',daq);
[TheReports,GiveMeReportsErr] = PsychHID('GiveMeReports',daq);
% err=PsychHID('SetReport',daq,2,70,uint8(70)); % GetAll
SetReportErr=PsychHID('SetReport',daq,2,70,uint8(0)); % GetAll
if SetReportErr.n
  fprintf('GetAll SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end
% Don't ask me... Denis is responsible for this next conditional statement --
% mpr
if 0
  % The last value passed in the argument list in this next statement is
  % almost certainly wrong for the 1608FS, but since this code is effectively
  % commented out, I don't see any upside to fixing it.
  [report,err]=PsychHID('GetReport',daq,1,70,35);
  if err.n
    fprintf('GetAll GetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
    data=[];
  end
else
  ReceiveReports2Err=PsychHID('ReceiveReports',daq);
  [TheReports2,GiveMeReports2Err]=PsychHID('GiveMeReports',daq);
  if GiveMeReports2Err.n
    fprintf('GetAll GiveMeReports error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
  end
  report=[];
  for i=1:length(TheReports2)
    if TheReports2(i).report(1)==70
      report=TheReports2(1).report;
    end
  end
end
if ~isempty(report) && report(1)==70
  % 35 byte report
  data.analog=zeros(1,16);
  for i=1:16
    data.analog(i)=report(i*2)+256*report(i*2+1);
    % if sign bit is set, make it negative
    if data.analog(i)>(2^16-1)
      data.analog(i)=data.analog(i)-2^16;
    end
  end
  data.digital=double(report(34:35));
else  
  data=[];
  fprintf(['\nIf you are using a 1208FS, you should consider asking Denis Pelli if he\n' ...
           'ever got this function working on that device.  I never got anything sensible\n' ...
           'from it on a 1608FS, but I eventually decided I would never use the function,\n' ...
           'so I gave up trying.  The "if 0" section of code looks like something Denis\n' ...
           'tried as he worked on this function.  Which suggests he might not have gotten it\n' ...
           'working on the 1208FS either.  If you find anything wrong in what I just wrote\n' ...
           'or if you can fix this function for any device, please edit my comments and\n' ...
           'update the PsychToolbox version. -- Mickey P. Rowe, 27-Dec-2007\n\n']);
  if Is1608
    for k=(daq-6):(daq-1)
      if strcmp(AllHIDDevices(k).serialNumber,AllHIDDevices(daq).serialNumber)
        err=PsychHID('ReceiveReports',k);
        [ThisReport,ThisErr] = PsychHID('GiveMeReports',k);
        if ~isempty(ThisReport)
          fprintf('Hold the presses... I got a report from DeviceIndex %d!\n',k);
        end
        tmperr=PsychHID('ReceiveReportsStop',k);
      end
    end
  end
end
ReceiveReportsStop2Err=PsychHID('ReceiveReportsStop',daq);
return

