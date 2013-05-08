function params=DaqAOutScan(daq,v,options)
% params=DaqAOutScan(DeviceIndex,v,options)
% USB-1208FS: Analog output scan. Produce sampled analog output voltage
% waveforms on one or two channels. This command sends the values in "v"
% (one column per channel) to the specified range of (one or two) output
% channels.
% 
% "params.fActual" is the actual sampling frequency, sample/channel/s.
%       It is as close as possible to the requested sampling frequency
%       options.f.
% "params.countActual" is the actual total sample/channel sent.
% "params.start" is when the first report ("v" data) was sent to the device.
% "params.end" is when the final report was sent to device.
% "DeviceIndex" is a small integer, the array index specifying which HID
%       device in the array returned by PsychHID('Devices') is interface 0 
%       of the desired USB-1208FS box. 
% "v" is a matrix, with one column per channel. (If you're using only one  
%       channel, then you're allowed to send the vector as either a column
%       or a row.) Each value is a double, in the range 0 to 1, which will
%       produce an output voltage in the range 0 to 4.095 V.
% "options.FirstChannel" is the first channel of the scan (0 - 1). (formerly
%       "options.lowChannel" -- that terminology is deprecated.)
% "options.LastChannel" is the last channel of the scan (0 - 1), and must 
%       be greater than or equal to options.FirstChannel. The values
%       options.FirstChannel and options.LastChannel specify the channel
%       range for the scan. (Formerly "options.highChannel" -- that terminology
%       is deprecated.)
% "options.f" is the desired sampling frequency, sample/channel/s, in the
%       range 0.596/c to 10e6/c Hz, where c is the number of channels to be
%       scanned.
% "options.trigger" is 1 = use external trigger. The external trigger 
%       may be used to start the scan synchronously. If the bit is set, the
%       device will wait until the appropriate trigger edge is detected,
%       then begin emitting samples at the specified rate. 
% "options.getReports" is 1 (default) = wait to receive report from
%       USB-1208FS before sending each block of data in "v". 0 = send
%       blocks of data at what should be the right rate.
% "options.print" is 1 = enable diagnostic printing; 0 (default) no 
%       diagnostics.
% 
% LIMITATIONS: As of 17 April 2005, DaqAoutScan works fine at up to 2000
% sample/s with one channel, or up to 200 sample/channel/s with two
% channels. (It is surprising that these two limits are not exactly a
% factor of 2 apart.) When I try to go faster than that I encounter various
% problems. Sometimes I get "data underrun" errors and sometimes the
% PsychHID SetReport call hangs up forever, never returning, making it
% necessary to quit MATLAB. I don't understand why these things happen;
% they shouldn't. I am trying to get help from Measurement Computing and
% Apple to understand and overcome this problems. However, in the meantime,
% the attained rates are enough for my immediate needs, and hopefully that
% will be true for many other users. 
% 
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest,
% DaqDeviceIndex, DaqDIn, DaqDOut, DaqAIn, DaqAOut, DaqAInScan.

% 4/15/05 dgp Wrote it.
% 4/15/05 dgp Merged several arguments into the "options" struct.
% 1/10/08 mpr glanced at it; made same changes as in DaqAOut.
% 1/14/08 mpr changed lowChannel/highChannel terminology to
%               FirstChannel/LastChannel a la DaqAInScan
% 3/15/09 mk  timer_preload calculation changed according to bug report and bugfix
%             suggested by Peter Meilstrup in forum message 9221. There was an
%             off-by-one bug present...

AllHIDDevices = PsychHIDDAQS;
if strcmp(AllHIDDevices(daq).product(5:6),'16')
  error('It looks like you are trying to run DaqAOutScan from a USB-1608FS, but such devices have no analog outputs.');
end

if isfield(options,'lowChannel')
  if isfield(options,'FirstChannel')
    if length(options.FirstChannel) ~= length(options.lowChannel) || ~all(options.FirstChannel == options.lowChannel)
      error('"options.lowChannel" is deprecated; new name is options.FirstChannel.  Do not try to use both inconsistently!');
    end
  else
    options.FirstChannel = options.lowChannel;
  end
end

if isfield(options,'highChannel')
  if isfield(options,'LastChannel')
    if length(options.LastChannel) ~= length(options.highChannel) || ~all(options.LastChannel == options.highChannel)
      error('"options.highChannel" is deprecated; new name is options.LastChannel.  Do not try to use both inconsistently!');
    end
  else
    options.LastChannel = options.highChannel;
  end
end

if ~isfield(options,'trigger')
    options.trigger=0;
end
if ~isfield(options,'getReports')
    options.getReports=1;
end
if ~isfield(options,'FirstChannel')
    options.FirstChannel=0;
end
if ~isfield(options,'LastChannel')
    options.LastChannel=0;
end
if ~isfield(options,'f')
    options.f=200;
end
if ~isfield(options,'print')
    options.print=0;
end
if ~ismember(options.LastChannel,0:1) || ~ismember(options.FirstChannel,0:1) 
    error('options.FirstChannel and options.LastChannel must each be 0 or 1.');
end
if options.LastChannel<options.FirstChannel
    error('options.LastChannel must be greater than or equal to options.FirstChannel.');
end
if ~ismember(options.trigger,0:1)
    error('If specified, options.trigger must be 0 or 1.');
end
c=1+options.LastChannel-options.FirstChannel;
if c==1 && size(v,1)==1 && size(v,2)>1
    v=v';
end
if c~=size(v,2)
    error('Number of columns in v must equal the number of channels in the scan.');
end
count=size(v,1); % sample/channel
if isinf(count)
    options.counted=0;
else
    if count<0 || count>intmax('uint32');
        error('count is out of 32-bit range, yet not INF.');
    end
    options.counted=1;
end
timer_prescale=ceil(log2(10e6/65535/options.f)); % Use smallest integer timer_prescale consistent with options.f.
timer_prescale=max(0,min(8,timer_prescale)); % Restrict to available range.

% 3/15/09 mk  timer_preload calculation changed according to bug report and bugfix
%             suggested by Peter Meilstrup in forum message 9221. There was an
%             off-by-one bug present...
timer_preload=round(10e6/2^timer_prescale/options.f) - 1;
timer_preload=max(0,min(65535,timer_preload)); % Restrict to available range.
params.fActual=10e6/2^timer_prescale/(timer_preload+1);

if abs(log10(params.fActual/options.f))>log10(1.1)
    error('Nearest attainable sampling frequency %.4f kHz is too far from requested %.4f kHz.',...
        params.fActual/1000,options.f/1000);
end
% v
v=v'; % transpose, yielding one row per channel.
v=v(:); % interleave the channels in one long vector.
v=max(0,min(1,v)); % clip
v=round(65535*v); % scale up to unsigned 16-bit integer
v8=uint8([bitand(v,255);bitshift(v,-8)]);
v8=v8(:);
% report
report=uint8(zeros(1,11));
report(1)=21;
report(2)=options.FirstChannel;  % the first channel of the scan (0 - 1)
report(3)=options.LastChannel; % the last channel of the scan (0 - 1)
ct=count;
report(4)=bitand(ct,255); % count (4 bytes). The number of scans to perform.
ct=bitshift(ct,-8);
report(5)=bitand(ct,255); % count (4 bytes). 
ct=bitshift(ct,-8);
report(6)=bitand(ct,255); % count (4 bytes). 
ct=bitshift(ct,-8);
report(7)=bitand(ct,255); % count (4 bytes). 
report(8)=timer_prescale; % timer_prescale
report(9)=bitand(timer_preload,255);   % timer_preload (2 bytes), low byte
report(10)=bitshift(timer_preload,-8); % timer_preload (2 bytes), hi byte
report(11)=options.counted+2*options.trigger; % options; bit field that controls various options
% 	bit 0: 	1 = single execution, 0 = continuous execution
% 	bit 1: 	1 = use external trigger
% 	bits 2-7: 	not used
if options.print
    report
end
err=DaqAOutStop(daq); % It might already be scanning. Stop it.
err=PsychHID('SetReport',daq,2,21,report); % AOutScan
if err.n
    fprintf('DaqAOutScan SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end
startingTime=GetSecs;
quittingTime=startingTime+count/params.fActual+0.2;
nextTime=startingTime;
samplesSent=0;
first=1;
reports=[];
params.start=[];
params.end=[];
params.numberOfReportsReceived=0;
params.numberOfReportsSent=0;
while samplesSent<count*c && GetSecs<quittingTime
    if ~first
        samplesSendEnd=min(samplesSent+32,count*c);
        WaitSecs(0.001);
        if options.print
            fprintf('%.3f s. Sending %d samples %d:%d\n',...
                GetSecs-startingTime,samplesSendEnd-samplesSent,samplesSent+1,samplesSendEnd);
        end
        err=PsychHID('SetReport',daq-1,2,0,v8(2*samplesSent+1:2*samplesSendEnd)); % send to Interface 1
    %     err=PsychHID('SetReport',daq-1,2,0,uint16(v(samplesSent+1:samplesSendEnd))); % send to Interface 1
        if options.print
            fprintf('%.3f s. Sent.\n',GetSecs-startingTime);
        end
        if err.n
            fprintf('DaqAOutScan interface 1, SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
        end
        params.numberOfReportsSent=params.numberOfReportsSent+1;
        nextTime=startingTime+samplesSent/c/options.f;
        samplesSent=samplesSendEnd;
    %     nextTime=startingTime+samplesSent/c/options.f;
        if samplesSent==count*c
            break;
        end
        if isempty(params.start)
            params.start=GetSecs;
        end
        params.end=GetSecs;
    end
    report=[];
    while isempty(report) && GetSecs<quittingTime && (options.getReports || GetSecs<nextTime)
        % Wait for AOutScan report from the USB-1208FS, requesting data.
        WaitSecs(0.001);
        if 1
            % hex2dec('15') is 21.
            [report,err]=PsychHID('GetReport',daq,1,21,2);
            if err.n
                fprintf('DaqAOutScan GetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
            end
        else
            err=PsychHID('ReceiveReports',daq);
            if err.n
                fprintf('DaqAOutScan ReceiveReports error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
            end
            [rr,err]=PsychHID('GiveMeReports',daq);
            if err.n
                fprintf('DaqAOutScan GiveMeReports error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
            end
            reports=[reports rr];
            report=[];
            for i=1:length(reports)
                if length(reports(i).report)>0
                    report=reports(i).report;
                    break;
                end
            end
            if i<length(reports)
                reports=reports(i+1:end);
            else
                reports=[];
            end
        end
    end
    if ~isempty(report) && report(1)==21
        params.numberOfReportsReceived=params.numberOfReportsReceived+1;
        if options.print
            if first
                fprintf('0.000 s. Got first data request %.3f s after sending AOutScan command. report [%d %d]\n',GetSecs-startingTime,report);
            else
                fprintf('%.3f s. Got data request. report [%d %d]\n',GetSecs-startingTime,report);
            end
        end
        if first
            startingTime=GetSecs;
            quittingTime=startingTime+count/params.fActual+0.2;
        end
        switch report(2)
            case 0;
                % ok
            case 1;
                fprintf('%.3f s. DaqAOutScan USB-1208FS "data underrun". Scan terminated. report [%d %d]. \n',GetSecs-startingTime,report);
                break;
            otherwise;
                fprintf('%.3f s. DaqAOutScan USB-1208FS "unknown error". Continuing. report [%d %d]. \n',GetSecs-startingTime,report);
        end
    end
    first=0;
end
% fprintf('samplesSent %d, count*c %d GetSecs %.3f s, quittingTime %.3f s, params.end %.3f s.\n',...
%     samplesSent, count*c,GetSecs-startingTime,quittingTime-startingTime, params.end-startingTime);

params.countActual=samplesSent/c;

% Check to make sure that no extra data requests come after we send our
% last block of data. Receiving a data request now would indicate that we
% have not sent enough data, which would suggest a bug in the code above.
if options.print
    % Check for extra AOutScan reports from the USB-1208FS, requesting data.
    WaitSecs(0.001);
    if 0
        [report,err]=PsychHID('GetReport',daq,1,21,2);
        if err.n
            fprintf('DaqAOutScan GetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
        end
    else
        err=PsychHID('ReceiveReports',daq);
        if err.n
            fprintf('DaqAOutScan ReceiveReports error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
        end
        [reports,err]=PsychHID('GiveMeReports',daq);
        if err.n
            fprintf('DaqAOutScan GiveMeReports error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
        end
        for i=1:length(reports)
            report=reports(i).report;
            if ~isempty(report) && report(1)==21
                fprintf('%.3f s. DaqAOutScan received extra data request afterwards. report [%d %d]. \n',GetSecs-startingTime,report);
            end
        end
    end
end
err=DaqAOutStop(daq); 
err=PsychHID('ReceiveReports',daq);
PsychHID('GiveMeReports',daq); % flush any remaining input reports.

return;
