function [data,params]=DaqAInScan(daq,options)
% [data,params]=DaqAInScan(DeviceIndex,options)
%
% ****************************************************************************
% *                                                                          *
% * For USB-1608FS, make sure you look at the bottom of this help as there   *
% * are significant differences between the function of that device and of   *
% * the USB-1208FS.  "channel" means something a bit different, and the      *
% * recommendation for that device is oppostite of that for the USB-1208FS   *
% * with respect to using options.FIrstchannel and options.Lastchannel.  The *
% * differences indicating need for that change in recommendation are        *
% * unfortunate as the two different methods for choosing which channels to  *
% * scan were the most confusing things in this function.  But the           *
% * USB-1208FS can be told to select particular channels by utilizing the    *
% * gain queue, and the USB-1608FS does not allow that behavior.  See the    *
% * help for DaqALoadQueue for more insight.  -- mpr                         *
% *                                                                          *
% ****************************************************************************
%
% USB-1208FS: Analog input scan. Clocked analog to digital conversion. This
% command samples several analog input channels at a specified rate and
% sends the readings to the host. Use DaqAInScan when you want to complete
% the whole operation in one call, tying up the computer until the sampling
% is done. Use DaqAInScanBegin, DaqAInScanContinue, and DaqAInScanEnd,
% instead, when you want to use the computer during sampling.
%
% USE CHANNEL AND RANGE: You can specify the channels to be scanned in two
% different ways, "FirstChannel" and "LastChannel", or as an arbitrary
% "channel" list. We recommend the latter because it's much more general.
% With "channel" you also provide a corresponding list of the desired gain
% "range". If you use FirstChannel and LastChannel, then the sampling will
% use the existing gains that were last set on those channels (by AIn or
% ALoadQueue). However, DaqAInScan doesn't know what that gain range was,
% so it arbitrarily assumes a gain range of 3 in computing the scale factor
% to convert your readings into voltages. Thus every reason recommends that
% you supply channel and range when you call DaqAInScan.
%
% "data" is an NxM matrix, with one column per channel. Each reading is a
%     double, based on a 16-bit value in the report. In the USB-1208FS
%     only the upper 12 of those 16 bits are significant.
% "params.fActual" is the actual sampling frequency, sample/channel/s. It
%     is as close as possible to the requested sampling frequency
%     options.f.
% "params.times" is the times (GetSecs) of receipt by PsychHID of the
%     reports in data.
% "DeviceIndex" is a small integer, the array index specifying which HID
%     device in the array returned by PsychHID('Devices') is interface 0
%     of the desired USB-1208FS box.
% "options.channel" is a vector of length 1 to 8. Each value (0 to 15)
%     selects any of various single-ended or differential measurements
%     using the eight digitizers. By default, DaqAInScan sends this for you
%     to the device by calling DaqALoadQueue before issuing the AInScan
%     command.
%     "channel" Measurement
%      0        0-1 (differential)
%      1        2-3 (differential)
%      2        4-5 (differential)
%      3        6-7 (differential)
%      4        1-0 (differential)
%      5        3-2 (differential)
%      6        5-4 (differential)
%      7        7-6 (differential)
%      8          0 (single-ended)
%      9          1 (single-ended)
%     10          2 (single-ended)
%     11          3 (single-ended)
%     12          4 (single-ended)
%     13          5 (single-ended)
%     14          6 (single-ended)
%     15          7 (single-ended)
%
% "options.range" is a vector of the same length as options.channel, with
%     values of 0 to 7, specifying the desired gain (and voltage range) for
%     the corresponding channel. By default, DaqAInScan sends this for you
%     to the device by calling DaqALoadQueue before issuing the AInScan
%     command. When options.range is not specified DaqAInScan assumes a
%     value of 3 in computing the scale factor applied to the results to
%     convert them into volts unless you are making a single-ended measurement.
%     For single-ended measures, the range is always +/- 10 V, so if you pass a
%     channel higher than 7, any range values you pass for that channel will be
%     ignored.
%
%     0 for Gain 1x (+-20 V),     1 for Gain 2x (+-10 V),
%     2 for Gain 4x (+-5 V),      3 for Gain 5x (+-4 V),
%     4 for Gain 8x (+-2.5 V),    5 for Gain 10x (+-2 V),
%     6 for Gain 16x (+-1.25 V),  7 for Gain 20x (+-1 V).
%
% "options.FirstChannel" (0 to 15) is the first channel of the scan. NOT
%     RECOMMENDED; use options.channel instead. If you specify FirstChannel
%     and LastChannel then don't specify channel and range, above.  (formerly
%     "options.lowChannel" -- that terminology is deprecated.)
% "options.LastChannel" (0 to 15) is the last channel of the scan. The
%     values FirstChannel and LastChannel specify the channel range for the
%     scan. If FirstChannel is greater than LastChannel, the scan will wrap
%     (i.e. if FirstChannel is 14 and LastChannel is 1, the scan will include
%     channels 14, 15, 0, and 1.) NOT RECOMMENDED; use options.channel
%     instead. If you specify FirstChannel and LastChannel then don't specify
%     channel and range, above. (formerly "options.highChannel" -- that
%     terminology is deprecated.)
% "options.count" is the desired number of samples per channel, either in
%     the 32-bit range 1 to 4e+09, or INF. The value INF invokes the
%     USB-1208FS's continuous data collection mode, which runs indefinitely.
% "options.ReleaseTime" effective only when options.continue (see below) is set
%     to 1 and options.count is Inf.  Setting options.ReleaseTime to some
%     reasonable time in the future and then calling DaqAInScanContinue will
%     allow data to be offloaded from the device's FIFO.  If you select the
%     option appropriately and call DaqAInScanContinue frequently enough, you
%     will not lose data.  Otherwise the FIFO fills, and you get jack diddly
%     squat.
% "options.f" is desired sampling frequency, sample/channel/s, in the range
%     0.596/c to 10e6/c Hz, where c is the number of channels being
%     scanned.
% "options.immediate" is 1 = immediate-transfer mode, 0 (default) =
%     block-transfer mode. At low sampling rates, immediate-transfer mode
%     will reduce the delay in receiving the sampled data. In
%     immediate-transfer mode, each 16-bit sample is sent immediately (a
%     2-byte report), rather than waiting for the buffer to fill (62
%     bytes). This mode should not be used if the aggregate sampling rate
%     is greater than 2,000 samples per second in order to avoid data loss.
%     The difference between the two modes will be more apparent if you
%     enable options.print, allowing you to see that you get one-sample
%     reports (2 bytes) in immediate mode, and less-frequent 31-sample
%     reports (62 bytes) in block-transfer mode.
% "options.trigger" is 1 = use external trigger; 0 (default).
%     The external trigger may be used to start data collection
%     synchronously. If the bit is set, the device will wait until the
%     appropriate trigger edge is detected, then begin sampling data at
%     the specified rate.  No messages will be sent until the trigger is
%     detected.
% "options.retrigger" is 1 = retrigger mode, 0 (default) = normal trigger.
%     The retrigger mode option is only used if trigger is used.  This
%     option will cause the trigger to be rearmed after options.count
%     samples are acquired if in continuous mode.
% "options.sendChannelRange is 1 (default) to ask that the options.channel
%     and options.range lists (if provided) be sent to the device by
%     calling DaqALoadQueue before the sending the AInScan command, or 0 to
%     not send them. Since those lists persist in the device you can save
%     time by skipping this transmission after you've already sent them
%     once (either through DaqAinScan or by explicitly calling
%     DaqALoadQueue yourself). Even if you set this option to 0 you should
%     still supply the options.channel and options.range arguments to
%     DaqAInScan because it needs to know how many channels you are using
%     in order to set up the count and frequency parameters in the AInScan
%     command to the USB-1208FS. This default is 0 if you supply
%     FirstChannel and LastChannel.
% "options.queue" is 1 = use channel gain queue (i.e. "channel" and
%     "range"), 0 = use options.FirstChannel and options.LastChannel
%     arguments. Don't bother setting this. It will be set for you, based
%     on which parameters you provide.
% "options.print" is 1 = enable diagnostic printing of the reports; 0
%     (default) no diagnostics.
% "options.begin" is 1 (default) to begin a new scan.
% "options.continue" is 1 (default) to receive reports (inside PsychHID).
% "options.end" is 1 (default) to wait until done and return the result.
%
% LIMITATION: The literature from Measurement Computing mentions a speed of
% 50 kHz, but that's a theoretical limit. As of 17 April 2005, DaqAInScan
% achieves 2000/c sample/channel/s, where c is the number of channels being
% sampled.
%
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest,
% DaqAInScanBegin, DaqAInScanContinue, DaqAInScanEnd,
% DaqDeviceIndex, DaqDIn, DaqDOut, DaqAIn, DaqAOut, DaqAInScan,DaqAOutScan.
%
% In some ways the USB-1608FS is less functional than the 1208FS.  For instance,
% where the 1208FS allows you to sample from the different channels in arbitrary
% order, the 1608FS only allows you to sample from sequential channels.  Hence,
% I think the only way to get this function to wrap for the 1608FS would be if
% you sampled from *all* channels and threw away the data not requested by the
% user.  Since that would entail a potentially mystifying performance hit, I
% instead opt to tell user that they cannot set FirstChannel higher than
% LastChannel.  What follows are other differences in options...  If an option
% is not specified below, then it behaves the same for the 1608FS as for the
% 1208FS.
%
% "options.channel" will be ignored (though you will get a warning).  Do *NOT*
%     even consider using this method with a 1608FS as it would be too confusing
%     to try to set this up... you could not poll an arbitrary sequence, and if
%     you want non-consecutive channels... see previous paragraph.
%
% "options.range" if specified, is either a scalar or a vector that must be of
%     the same length as the vector: options.FirstChannel:options.LastChannel.
%     Values must be in the set 0:7, specifying the desired gain (and hence
%     voltage range) for the corresponding channel as in the 1208FS except that
%     the mapping differs.  For the USB-1608FS, the values mean:
%
%     0 for Gain 1x (+/- 10 V),     1 for Gain 2x (+/- 5 V),
%     2 for Gain 4x (+/- 2.5 V),    3 for Gain 5x (+/- 2 V),
%     4 for Gain 8x (+/- 1.25 V),   5 for Gain 10x (+/- 1 V),
%     6 for Gain 16x (+/- 0.625 V), 7 for Gain 32x (+/- 0.3125 V)
%
% If options.range is passed, this function will call DaqALoadQueue for you.  If
% it is not, then it is assumed you do not wish to change the gains, and that
% function is not called.  Since we cannot read the gains, a preferences file
% will be created or edited for you when you call DaqALoadQueue.
%
% "options.FirstChannel" (0 to 7) is the first channel of the scan.  If you do not
%     pass this option, it is assumed to be 0.
% "options.LastChannel" (options.Firstchannel to 7) is the last channel of the
%     scan.  If you do not specify this option, it is assumed you want only to
%     sample from options.Firstchannel.
% "options.retrigger" is apparently not an option for the 1608FS.  I could be
%     wrong about whether it *could* be an option: there is a potential option
%     called "external sync", and that might be the same as retrigger.  However,
%     even if that is correct, I have not implemented it.  (see note below on
%     numerical definitions of options)
% "options.sendChannelRange" has no meaning here.  I will use preference file
%     instead.
% "options.queue" again, not an option.  Don't try to use it.
%
% These options do not exist (as such) for the 1208FS.  They all default to 0.
%
% "options.burst" if 1, acquire data in burst mode.  It is not clear to me what
%     differences are entailed by setting this to 1 compared to just setting
%     "options.immediate" to 1, but if nothing else, data from multiple channels
%     can be acquired faster; I've not tested it, but it may be that in burst
%     mode the device acquires data until its FIFO is full and may ignore
%     efforts to communicate with it until that occurs.
% "options.ExternSync" if 1, use external sync signal.  I believe this is mainly
%     for using to daqs to acquire data simultaneously.  Connect the sync
%     terminals of the two devices, configure one to output a sync signal and
%     set this option to 1 so that this device listens for that sync signal and
%     times its data acquisition accordingly.
% "options.debug" I have no idea what this does... "debug" doesn't even appear
%     in the user manual for the device.
%
% To use this function, you don't need the following information.  However, if
% you want to modify it, this information may prove useful to you.  The options
% that are passed through PsychHID are combined into a single 8-bit number.  The
% settings for the bits (that I know of) for the 1208FS are:
%    1 (0x1) Execution Mode (1=counted, 0=continuous)
%    2 (0x2) Transfer mode (1=immediate, 0=blocked)
%    3 (0x4) Trigger mode (1=external, 0=internal)
%    5 (0x10) Gain Queue mode (1=use stored queue, 0=use channels passed)
%    6 (0x20) Retrigger mode (1=reset trigger, 0=don't)
%
% and for the 1608FS:
%    1 (0x1) Execution Mode (1=counted, 0=continuous)
%    2 (0x2) Burst Mode (1=burst I/O, 0=normal I/O)
%    3 (0x4) Transfer mode (1=immediate, 0=blocked)
%    4 (0x8) Trigger mode (1=external, 0=internal)
%    5 (0x10) External sync (1=use external, 0=don't); as noted above, this may
%              be the same as retrigger mode
%    6 (0x20) Debug mode (1=debug, 0=non-debug) (I don't know what this does)
%
% 4/15/05 dgp Wrote it.
% 4/15/05 dgp Merged several of the arguments into the "options" struct.
% 4/25/05 dgp Added options.channel, options.range, and options.sendChannelRange.
% 4/27/05 dgp Fixed immediate-transfer mode once I realized that the report
%             contains only one sample, not one sample per channel.
% 8/26/05 dgp Fixed bug in extraction of serial number, as suggested by Steve
%             Van Hooser, vanhooser@neuro.duke.edu
% 8/26/05 dgp Incorporated bug fix for compatibility with Mac OS X Tiger
%             suggested by Maria Mckinley <parody@u.washington.edu>. The reported
%             number of outputs of the USB-1208FS has changed in Tiger.
%             http://groups.yahoo.com/group/psychtoolbox/message/3614
% 1/02/07 mk  Add 'persistent dinc' to line 245. This apparently fixes some
%             bug, bugfix proposed by Florian Stendel.
% 1/x/08  mpr modified for use of USB-1608FS, added options.ReleaseTime; changed
%               low and high Channels to First and Last because previous
%               terminology seemed too easily confused with high and low bytes
%               or high and low channels in 1208FS pin out diagrams.
% 3/14/08 mpr added warning when data are discarded
% 3/15/09 mk  timer_preload calculation changed according to bug report and bugfix
%             suggested by Peter Meilstrup in forum message 9221. There was an
%             off-by-one bug present...
% 3/24/12 mk  Add handling for options.livedata -- retrieval of data while
%             DAQ is running.

% These USB-1280FS parameters are computed from the user-supplied
% arguments.
% "timer_prescale" chooses the timer prescaler (0 - 8).
% 	0 = 1:1 prescale
% 	1 = 1:2 prescale
% 	2 = 1:4 prescale
% 	3 = 1:8 prescale
% 	4 = 1:16 prescale
% 	5 = 1:32 prescale
% 	6 = 1:64 prescale
% 	7 = 1:128 prescale
% 	8 = 1:256 prescale
% "timer_preload" is the 16-bit timer preload value.
% "options" controls various options
% 	bit 0: 	1 = single execution, 0 = continuous execution
% 	bit 1:	1 = immediate-transfer mode, 0 = block transfer mode
% 	bit 2: 	1 = use external trigger
% 	bit 3:	not used
% 	bit 4:	1 = use channel gain queue, 0 = use channel parameters specified
% 	bit 5:	1 = retrigger mode, 0 = normal trigger
% 	bits 6-7: 	not used
%
% The sample rate is set by the internal 16-bit incrementing timer running
% at a base rate of 10 MHz.  The timer is controlled by timer_prescale and
% timer_preload.  These values are only used if the device has been set to
% master the SYNC pin with the SetSync command. The timer will provide an
% internal interrupt when its value rolls over.  The timer is preloaded
% with the value specified in timer_preload, and counts down.  This allows
% for a lowest rate of 0.596 Hz (1:256 prescale, preload = 65535).  It is
% preferable to keep the prescaler to the lowest value that will achieve
% the desired rate.

% MPR changed number of outputs so that this function would work with a
% USB-1608FS under Leopard.

% Perform internal caching of list of HID devices to speedup call:
persistent AllHIDDevices;
if isempty(AllHIDDevices)
    AllHIDDevices = PsychHIDDAQS;
end

if strcmp(AllHIDDevices(daq).product(5:6),'16')
    Is1608 = 1;
    
    % These options differ from those of 1208
    if ~isfield(options,'burst')
        options.burst = 0;
    end
    if ~isfield(options,'ExternSync')
        options.ExternSync = 0;
    end
    if ~isfield(options,'debug')
        options.debug = 0;
    end
else
    Is1608 = 0;
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


% Fill in defaults for any unspecified options.
% Parameters of the AInScan command sent to the USB-1x08FS.
if ~isfield(options,'channel')
    options.channel=[];
elseif Is1608
    warning('Psych:Daq:UnusedOption', '"channel" is not a useable option for USB-1608FS and will be ignored.');
end
if ~isfield(options,'range')
    options.range=[];
end
if ~Is1608 && length(options.channel)~=length(options.range)
    error('"options.channel" and "options.range" vectors must be of the same length.');
end
if Is1608
    channelRangeOk = 0;
else
    channelRangeOk=~isempty(options.channel) && ~isempty(options.range);
end
if ~isfield(options,'FirstChannel')
    if Is1608
        options.FirstChannel = 0;
    else
        options.FirstChannel=[];
    end
end
if ~isfield(options,'LastChannel')
    if Is1608
        options.LastChannel = options.FirstChannel;
    else
        options.LastChannel=[];
    end
end
FirstLastOk=~isempty(options.FirstChannel) && ~isempty(options.LastChannel);
if FirstLastOk==channelRangeOk
    error('Please specify either options.channel and options.range OR options.FirstChannel and options.LastChannel.');
end
if length(options.FirstChannel)>1 || length(options.LastChannel)>1
    error('options.FirstChannel and options.LastChannel, if specified, must be scalars.');
end
if ~isfield(options,'count')
    options.count=1000;
end
options.count=round(options.count); % make it an integer.
if ~isfield(options,'f')
    options.f=1000;
end
if ~isfield(options,'immediate')
    options.immediate=0;
end
if ~isfield(options,'trigger')
    options.trigger=0;
end
if ~isfield(options,'queue')
    options.queue=channelRangeOk;
elseif Is1608
    options.queue = 0;
    warning('Psych:Daq:UnusedOption', '"queue" option has no meaning for USB-1608FS, so will be ignored');
end
if options.queue~=channelRangeOk
    error('options.queue setting inconsistent with other parameters. Let us set it for you.');
end
if ~isfield(options,'retrigger')
    options.retrigger=0;
end
% User options.
if ~isfield(options,'sendChannelRange')
    if ~Is1608
        options.sendChannelRange=channelRangeOk;
    end
elseif Is1608
    warning('Psych:Daq:UnusedOption', '"sendChannelRange" has no meaning for USB-1608FS, so will be ignored.');
end
if ~Is1608 && (options.sendChannelRange && ~channelRangeOk)
    error('You are not supplying channel and range, so you can''t send them. Omit options.sendChannelRange.');
end
if ~isfield(options,'print')
    options.print=0;
end
if ~isfield(options,'begin')
    options.begin=1;
end
if ~isfield(options,'continue')
    options.continue=1;
end
if ~isfield(options,'end')
    options.end=1;
end

% This function may be called from any of DaqAInScanBegin, DaqAInScanContinue,
% or DaqAInScanEnd.  Keeping these persistent means they will remain even when
% options.begin = 0
persistent start;
persistent IndexRange;

if options.begin
    % It might be running, so stop it.
    % hex2dec('12') is 18.
    err=PsychHID('SetReport',daq,2,18,uint8(0)); % AInStop
    if err.n
        fprintf('AInStop SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
    end
    % The user supplies us only the device index "daq" corresponding to
    % interface 0. The reports containing the samples arrive on interfaces
    % 1,2,3. As returned by PsychHID('Devices'), interface "i" is at
    % daq-i, and we proceed on that basis after doing a quick check to
    % confirm our assumption. However, to be platform-independent, it would
    % be better to actually find all four device interfaces and confirm
    % their interface numbers. USB Probe reports interface numbers, but, as
    % far as I can tell, Apple's HID Explorer and our PsychHID do not
    % return this information. However, PsychHID('Devices') does report the
    % number of outputs: the USB-1208FS interface 0 has 229 (pre-Tiger) 70 (Tiger)
    % outputs, interface 1 has 65 (pre-Tiger) 1 (Tiger) outputs, and interfaces 2
    % and 3 have zero outputs.  I have no idea why the number of outputs changed
    % with the arrival of Mac OS X Tiger. -- author unknown
    %
    % The USB-1608FS interface 0 has only 66 outputs under Leopard.  It should
    % have seven interfaces.  I have frequently found that PsychHID finds only six
    % or (rarely) five interfaces.  I have not figured out why, but usually
    % subsequent enumerations find all seven.  I have tried to take care of this
    % in DaqFind.  Nevertheless, to be safe I will add a check here to make
    % certain the interfaces are consecutive, and that there are seven of them.
    % -- mpr
    
    if Is1608
        SN = AllHIDDevices(daq).product;
        
        % Fill empty .product fields with filler, otherwise the
        % following strvcat deletes them!
        for kk = 1:length(AllHIDDevices)
            if isempty(AllHIDDevices(kk).product)
                AllHIDDevices(kk).product = 'xoxo';
            end
        end
        
        AllSNs = strvcat(AllHIDDevices.product);
        InterfaceInds = strmatch(SN,AllSNs);
        if length(InterfaceInds) ~= 7 || ~all(InterfaceInds' == (daq-6):daq)
            % Horrible hack for the horrible 64-Bit OSX:
            if ~IsOSX(1)
                error('Not all interfaces found.  Run "help DaqReset" for suggestions.');
            else
                warning('Not all 7 interfaces found. Will fake most common interface config and hope for the best. Run "help DaqReset" for suggestions.');
            end
        end
        IndexRange = -1:-1:-6;
    else
        % Find all other interfaces of device 'daq', by looking for devices
        % with the same serial number:
        SN = AllHIDDevices(daq).product;

        % Fill empty .product fields with filler, otherwise the
        % following strvcat deletes them!
        for kk = 1:length(AllHIDDevices)
            if isempty(AllHIDDevices(kk).product)
                AllHIDDevices(kk).product = 'xoxo';
            end
        end

        AllSNs = strvcat(AllHIDDevices.product);
        InterfaceInds = transpose(strmatch(SN,AllSNs));
        if length(InterfaceInds) ~= 4
            % Horrible hack for the horrible 64-Bit OSX:
            if ~IsOSX(1)
                error('Not all interfaces found.  Run "help DaqReset" for suggestions.');
            else
                warning('Not all 4 interfaces found. Will fake most common interface config and hope for the best. Run "help DaqReset" for suggestions.');
            end
        end
        
        % Throw out the primary interface with index 'daq':
        InterfaceInds = InterfaceInds(find(InterfaceInds ~= daq));
        
        % Convert to indices/range relative to 'daq', as needed later on:
        IndexRange = InterfaceInds - daq;
        
        % Horrible hack for the horrible 64-Bit OSX:
        if IsOSX(1)
          % Hardcode index range, in the hope it helps that brain-dead os:
          IndexRange = -1:-1:-3;
        end
    end

    % Flush any stale reports.
    for d=IndexRange % Interfaces 1,2,3 (1208FS) or 1:6 (1608FS)
        err=PsychHID('ReceiveReports',daq+d);
    end
    for d=IndexRange % Interfaces 1,2,3 (1208FS) or 1:6 (1608FS)
        [reports,err]=PsychHID('GiveMeReports',daq+d);
        if ~isempty(reports) && options.print
            fprintf('Flushing %d stale reports from DeviceIndex %d.\n',length(reports),daq+d);
        end
    end
end
err.n = 0;
if Is1608
    % How many channels?
    c = options.LastChannel-options.FirstChannel+1;
    if length(options.range) ~= c && ~ismember(numel(options.range),0:1)
        error('Specified length of range is %d and number of channels is %d',length(options.range),options.LastChannel-options.FirstChannel+1);
    end
    if isempty(options.range)
        PrefsNotFound = 1;
        DaqPrefsDir = DaqtoolboxConfigDir;
        if exist([DaqPrefsDir filesep 'DaqPrefs.mat'],'file')
            DaqVars = load([DaqPrefsDir filesep 'DaqPrefs']);
            if isfield(DaqVars,'OldGains')
                options.range = DaqVars.OldGains((options.FirstChannel:options.LastChannel)+1);
                PrefsNotFound = 0;
            end
        end
        if PrefsNotFound
            fprintf(['\n\nDanger, Will Robinson!!  Danger!!\n\n' ...
                'Input options did not specify the ranges, and I could not find the values in a\n' ...
                'Preference file.  I will assume that gains are all 1 (range +/- 10 V), but if that\n' ...
                'is not what your daq thinks the gains are, then your data will be off by a scale\n' ...
                'factor or worse because calibration may be quite off too.\n\n' ...
                'I strongly recommend that you run DaqALoadQueue or use options.range in your next\n' ...
                'call to DaqAInScan.\n\n']);
            options.range = zeros(1,c);
        end
    else % if isempty(options.range)
        if numel(options.range) == 1
            options.range = options.range*ones(1,c);
        end
        err=DaqALoadQueue(daq,options.FirstChannel:options.LastChannel,options.range);
    end % if isempty(otptions.range); else
else % if Is1608
    if options.sendChannelRange
        % for single-ended measurements, there is only one value for range
        options.range(find(options.channel > 7)) = zeros(size(find(options.channel > 7)));
        err=DaqALoadQueue(daq,options.channel,options.range);
    end
    % How many channels?
    if channelRangeOk
        c=length(options.channel);
        % will need channel to be defined to check for single-ended measurements
        % below in order to ge scale correct.
        channel=options.channel;
    else
        if ~ismember(options.LastChannel,0:15) || ~ismember(options.FirstChannel,0:15)
            error('options.FirstChannel and options.LastChannel must each be in the range 0:15.');
        end
        if options.FirstChannel<=options.LastChannel
            channel=options.FirstChannel:options.LastChannel;
        else
            channel=[options.FirstChannel:15 0:options.LastChannel]; % wrap around.
        end
        c=length(channel);
    end % if channelRangeOk; else
end % if Is1608; else
if err.n
    fprintf('DaqALoadQueue error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end
% Convert from user to USB-1208FS usage. "options.f" and "options.count" refer to
% samples per channel.
% USB-1208FS: pmd.f and pmd.count refer to
% samples.
pmd.f=c*options.f;
pmd.count=c*options.count;
if isinf(pmd.count)
    options.counted=0;
else
    if pmd.count<0 || pmd.count>intmax('uint32');
        error('pmd.count is out of 32-bit range, yet not INF.');
    end
    options.counted=1;
end

timer_prescale=ceil(log2(10e6/65535/pmd.f)); % Use smallest integer timer_prescale consistent with pmd.f.
timer_prescale=max(0,min(8,timer_prescale)); % Restrict to available range.

% MK: timer_preload calculation changed according to bug report and bugfix
% suggested by Peter Meilstrup in forum message 9221. There was an
% off-by-one bug present...
timer_preload=round(10e6/2^timer_prescale/pmd.f) - 1;
timer_preload=max(0,min(65535,timer_preload)); % Restrict to available range.
pmd.fActual=10e6/2^timer_prescale/(timer_preload+1);

% Convert from USB-1208FS to user usage.
% "options.f" and "options.count" refer to samples/channel/s and
% samples/channel.
params.fActual=pmd.fActual/c;
if options.begin
    if abs(log10(params.fActual/options.f))>log10(1.1)
        error('Nearest attainable sampling frequency %.4f kHz is too far from requested %.4f kHz.',params.fActual/1000,options.f/1000);
    end
    report=uint8(zeros(1,11));
    report(1)=17;
    if Is1608 || ~options.queue
        report(2) = options.FirstChannel;
        report(3) = options.LastChannel;
    end
    if isfinite(pmd.count)
        countsave=pmd.count;
        report(4)=bitand(pmd.count,255); % count
        pmd.count=bitshift(pmd.count,-8);
        report(5)=bitand(pmd.count,255);
        pmd.count=bitshift(pmd.count,-8);
        report(6)=bitand(pmd.count,255);
        pmd.count=bitshift(pmd.count,-8);
        report(7)=bitand(pmd.count,255);
        pmd.count=countsave;
    end
    report(8)=timer_prescale; % timer_prescale
    report(9)=bitand(timer_preload,255); % timer_preload
    report(10)=bitshift(timer_preload,-8);
    if Is1608
        report(11) = options.counted+2*options.burst+4*options.immediate+8*options.trigger+16*options.ExternSync+32*options.debug;
    else
        report(11)=options.counted+2*options.immediate+4*options.trigger+16*options.queue+32*options.retrigger;
    end
    % Tell the USB1208FS to start sampling and sending reports.
    % hex2dec('11') is 17.
    err=PsychHID('SetReport',daq,2,17,report); % AInScan
    start=GetSecs;
    if err.n
        fprintf('AInScan SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
    end
end

if options.continue
    % AInScan responses from the USB-1208FS will be sent sequentially from
    % interrupt IN endpoints in interfaces 1-3 as input reports with a
    % report ID of 0. The last two bytes are a sequential index,
    % "reports.serial".
    if isfinite(options.count)
        % Perform exactly one pass per interface:
        for d=IndexRange % Interfaces 1,2,3 (1208FS) or 1:6 (1608FS)
            % Tell PsychHID to receive reports, storing them internally.
            err=PsychHID('ReceiveReports',daq+d,options);
            if err.n
                fprintf('AInScan DeviceIndex %d, ReceiveReports error 0x%s. %s: %s\n',daq+d,hexstr(err.n),err.name,err.description);
            end
        end
    else % isfinite(options.count)
        if isfield(options,'ReleaseTime')
            while GetSecs < options.ReleaseTime
                for d=IndexRange
                    err=PsychHID('ReceiveReports',daq+d,options);
                    if err.n
                        fprintf('AInScan DeviceIndex %d, ReceiveReports error 0x%s. %s: %s\n',daq+d,hexstr(err.n),err.name,err.description);
                    end
                end
            end
        end
    end
end

% Large parts of code are shared between end of aquisition with return of
% all captured data, and retrieval of livedata without stop of aquisition:
if options.end || (isfield(options, 'livedata') && options.livedata)
    % Get and combine reports from all interfaces.
    r=[];
    params.times=[];
    for d=IndexRange
        [reports,err]=PsychHID('GiveMeReports',daq+d);
        if err.n
            fprintf('AInScan DeviceIndex %d, GiveMeReports error 0x%s. %s: %s\n',daq+d,hexstr(err.n),err.name,err.description);
        end
        r=[r reports];
        
        % Only stop reception if this is really end of aquisition:
        if options.end
            err=PsychHID('ReceiveReportsStop',daq+d);
            if err.n
                fprintf('AInScan DeviceIndex %d, ReceiveReportsStop error 0x%s. %s: %s\n',daq+d,hexstr(err.n),err.name,err.description);
            end
        end
    end
    
    % Stop DAQ and flush receive buffers if this is really end of aquisition:
    if options.end
        % hex2dec('12') is 18.
        err=PsychHID('SetReport',daq,2,18,uint8(0)); % AInStop
        if err.n
            fprintf('AInStop SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
        end
        for d=IndexRange
            err=PsychHID('ReceiveReports',daq+d);
            if err.n
                fprintf('AInScan DeviceIndex %d, ReceiveReports error 0x%s. %s: %s\n',daq+d,hexstr(err.n),err.name,err.description);
            end
            [reports,err]=PsychHID('GiveMeReports',daq+d);
            if err.n
                fprintf('AInScan DeviceIndex %d, GiveMeReports error 0x%s. %s: %s\n',daq+d,hexstr(err.n),err.name,err.description);
            end
            if isfinite(options.count)
                if ~isempty(reports)
                    fprintf('WARNING: received %d post-deadline reports from DeviceIndex %d.\n',length(reports),daq+d);
                end
            else
                r=[r reports]; %#ok<*AGROW>
            end
            err=PsychHID('ReceiveReportsStop',daq+d);
            if err.n
                fprintf('AInScan DeviceIndex %d, ReceiveReportsStop error 0x%s. %s: %s\n',daq+d,hexstr(err.n),err.name,err.description);
            end
        end
    end
    
    % Post processing of raw report data:
    reports=r;
    if isempty(reports)
        data=[];
        fprintf('Nothing received.\n');
        return
    end
    % Extract the report serial number provided by the USB-1208FS
    for k=1:length(reports)
        reports(k).serial=double(reports(k).report(63))+double(reports(k).report(64))*256;
    end
    
    % Sort by serial number.
    [r,ri]=sort([reports.serial]);
    reports=reports(ri);
    
    InitReports = reports;
    
    % Keep only the consecutively numbered reports, from 0 on. Discard
    % everything after a gap.  -- author presumably dgp...
    %
    % This turned out to be a real problem for me; not sure why, but for some
    % reason a report or two can get lost even though the FIFO buffer is not full
    % and you have not reached the limit on the number of reports.  I was losing
    % data hand over fist with my USB-1608; not sure what conditions cause it, but
    % here's some of what I can say: I found I cannot use DaqAInScanContinue,
    % because it takes too long to send "ReceiveReports" for every interface, so
    % when timing is an issue, I try to send just enough calls to
    % "PsychHID('ReceiveReports',daq+k);" to prevent the FIFO from overflowing.
    % For "k" I use an index of the form k=rem(k,NumberOfInterfaces)+1; so that
    % each successive call asks for the next interface in order.  But for some
    % reason -- I think particularly if the computer spends a bit too much time
    % doing other things -- a report or two gets dropped.  I've added a check for
    % that.  Previously data were discarded here without warning!
    ri=[reports.serial]==0:length(reports)-1;
    reports=reports(ri);
    
    if length(InitReports) ~= length(reports)
        if length(InitReports)-length(reports) == 1
            diffindex = diff(r);
            % I think the index math is such that this will always be
            % length(reports)-1, but to be safe:
            BadOne = r(find(diffindex > 1)); %#ok<*FNDSB>
            fprintf('\n1 report discarded (%d reports received, but report(s) missed after report %d).\n', ...
                length(InitReports),BadOne);
        else
            % This warns user if reports are thrown away; if you get this warning
            % unacceptably often, I suggest you first try to modify your code to put
            % in more frequent calls to PsychHID('ReceiveReports, ...)  If you can't
            % do that or you can and it doesn't work, then you may be faced with the
            % hard choice of deciding how to handle the problem...  One thing you
            % might do is modify this code so that it doesn't just throw away the
            % reports obtained after a miss.  In my case, I frequently found that I'd
            % lose more than 9000 reports just because 1 report was missed after the
            % 500th (or so) report.  I got the problem to go away by reducing the time
            % the computer spent doing other things between ReceiveReports calls, but
            % I can imagine wanting to make due with all the data I managed to get
            % from the device rather than just the first sequential reports.  So you
            % could rewrite the above code so that instead of throwing away data
            % acquired after a missed report, the function returns all the data it
            % gets back along with a vector specifying the indices.  Difficulties
            % you'll face: meaning of sequence numbers differ according to value of
            % options.immediate, fact that params.times isn't helpful (time that
            % reports are received from device not monotonically related to time that
            % data within said reports are acquired).
            fprintf('\n%d reports discarded! (out of %d reports received; %d reports missed)\n', ...
                length(InitReports)-length(reports),length(InitReports),InitReports(end).serial+1-length(reports));
        end
        %% For Diagnostic purposes:
        %
        % [ErrMsg,HomeDir] = unix('echo $HOME');
        %% trim trailing carriage return
        % HomeDir(end) = [];
        % FileNo = 0;
        % while exist([HomeDir 'Desktop/DiscardedReports' int2str(FileNo) '.mat'],'file')
        %   FileNo = FileNo+1;
        % end
        % save([HomeDir 'Desktop/DiscardedReports' int2str(FileNo)]);
    end
    % Diagnostic print of the raw reports
    if options.immediate
        % report length
        bytes=2;
    else
        bytes=62;
    end
    if options.print
        fprintf('Got %d reports.\n',length(reports));
        fprintf('Serial number, Device, Time received (since AInScan command), Report rate, Report as 16-bit words.\n');
        for k=1:length(reports)
            fprintf('Report %3d, DeviceIndex %d, %4.0f ms, %4.0f report/s.',...
                reports(k).serial,reports(k).device,1000*(reports(k).time-start),reports(k).serial/(reports(k).time-start));
            r=double(reports(k).report);
            r=r(1:bytes);
            fprintf(' %5d',r(1:2:end)+256*r(2:2:end));
            fprintf('\n');
        end
    end
    % Return times.
    params.times=[reports.time];
    % Extract the data from the reports.
    data=[];
    for k=1:length(reports)
        data=[data reports(k).report(1:bytes)];
    end
    
    % Combine two bytes for each reading.
    data = double(data(1:2:end))+double(data(2:2:end))*256;
    % Discard any extra 16-bit words at the end of the last report.
    if length(data)>c*options.count
        if 2*(length(data)-options.count*c)>60
            fprintf('Trimming off an extra %.1f sample/channel, %.0f bytes.\n',length(data)/c-options.count,2*(length(data)-options.count*c));
        end
        data=data(1:c*options.count);
    end
    if c*options.count>length(data)
        if isfinite(options.count)
            fprintf('Missing %.1f sample/channel, %.0f bytes.\n',options.count-length(data)/c,2*(options.count*c-length(data)));
        end
        options.count=floor(length(data)/c);
        data=data(1:c*options.count);
    end
    data=reshape(data,c,options.count)';
    
    range=ones([1 size(data,2)]);
    if ~isempty(options.range)
        if length(options.range)~=size(data,2)
            error('Wrong length of "options.range" vector.');
        end
        range(:)=options.range(:);
    else
        range=3*range;
    end
    
    if Is1608
        vmax=[10 5 2.5 2 1.25 1 0.625 0.3125];
        data = data/32768-1;
    else
        vmax=[20 10 5 4 2.5 2 1.25 1];
        
        DiffChannels = find(channel < 8);
        SE_Channels = find(channel > 7);
        
        data(:,DiffChannels) = bitshift(data(:,DiffChannels),-4);
        
		ix = data(:,DiffChannels) > 2048;
        data(ix) = -bitcmp(data(ix),12)-1;
		
        SE_Data = data(:,SE_Channels);
        
        OverflowInds = find(SE_Data > 32752);
        UnderflowInds = find(SE_Data > 32736);
        OKInds = find(SE_Data < 32737);
        
        SE_Data(OverflowInds) = -2048*ones(size(OverflowInds));
        SE_Data(UnderflowInds) = 2047*ones(size(UnderflowInds));
        SE_Data(OKInds) = bitand(4095,bitshift(SE_Data(OKInds),-3))-2048;
        
        data(:,SE_Channels) = SE_Data;
        
        data=data/2047;
        
        % for PsychHID calls, single-ended measurements must have range set to 0,
        % but
        % for vmax determination range must be 1 because scale is +/- 10 V
        range(SE_Channels) = ones(length(SE_Channels),1);
    end
    
    range=vmax(range+1);
    data=repmat(range,size(data,1),1).*data;
    
    if Is1608
        TheChannels = options.FirstChannel:options.LastChannel;
        DataUncalibrated = ones(1,c);
        DaqPrefsDir = DaqtoolboxConfigDir;
        if exist([DaqPrefsDir filesep 'DaqPrefs.mat'],'file')
            DaqVars = load([DaqPrefsDir filesep 'DaqPrefs']);
            if isfield(DaqVars,'CalData')
                CalData = DaqVars.CalData;
                for k=1:length(TheChannels)
                    GoodIndices = find(CalData(:,1) == TheChannels(k) & CalData(:,2) == options.range(k));
                    if ~isempty(GoodIndices)
                        TheDays = CalData(GoodIndices,end);
                        ThisDay = datenum(date);
                        [DaysSinceLast,BestIndex] = min(ThisDay-TheDays);
                        AllThatDay = find(TheDays == TheDays(BestIndex));
                        MostRecentPolyFit = CalData(GoodIndices(AllThatDay(end)),3:5);
                        
                        if DaysSinceLast > 30
                            warning('Psych:Daq:CalibOutdated', 'Calibration of channel %d has not been performed since %s!!',TheChannels(k),datestr(TheDays(BestIndex)));
                        end
                        
                        data(:,k) = polyval(MostRecentPolyFit,data(:,k));
                        DataUncalibrated(k) = 0;
                    end
                end
            end % if isfield(DaqVars,'CalData')
        end % if exist([Preferences file])
        if any(DataUncalibrated)
            warning('Psych:Daq:CalibMissing', 'It looks like some of your channels have not been calibrated.');
            TheIndices = find(DataUncalibrated);
            if length(TheIndices) == 1
                fprintf('I did not find calibration data for channel %d.\n\n',TheChannels(find(DataUncalibrated)));
            else
                fprintf('I did not find calibration data for channels ');
                for k=1:(length(TheIndices)-1)
                    fprintf('%d, ',TheIndices(k));
                end
                fprintf('and %d.\n\n',TheIndices(end));
            end
            fprintf('In my tests, uncalibrated values could be off by as much as 15%%!\n\n');
        end
    end % if Is1608
else % if options.end
    data=[];
    params.times=[];
end % if options.end; else

return;
