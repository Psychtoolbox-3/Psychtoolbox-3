% DaqWaitButtonDumbTest1408FS.m
%
% A quick and dirty demo on how to poll the digital inputs
% of a USB-1408FS as fast as possible, without overhead.
%
% This is essentially a stripped down version of DaqDIn() with
% some polling loops wrapped around and a options.secs value of
% zero for essentially polling with no timeout.
%

daq = DaqFind;
reportId = 3;
TheReport = uint8(0);
NumberOfPorts = 2;

if IsWin
 % Windows needs some minimal polling time:
 options.secs = 0.001;
else
 options.secs = 0.000;
end

PsychHID('ReceiveReportsStop',daq);
PsychHID('GiveMeReports',daq);
PsychHID('ReceiveReports',daq, options);

while 1
% Emit query to device:
PsychHID('SetReport',daq,2,reportId, TheReport);

% Wait for result from device:
inreport = [];
while isempty(inreport)
% Yield some minimum amount of time to other processes if you want to be
WaitSecs('YieldSecs', 0);
% New data available?
inreport = PsychHID('GetReport',daq,1,reportId,NumberOfPorts+1);
end

% inreport contains latest button query result: Button pressed?
if inreport(3) ~= 255
% yes: Done!
break;
end

% No. Repeat sampling...
end

fprintf('Change detected! Bye.\n');
