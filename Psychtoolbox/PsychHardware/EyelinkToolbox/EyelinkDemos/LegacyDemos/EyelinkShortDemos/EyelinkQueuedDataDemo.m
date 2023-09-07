function EyelinkQueuedDataDemo
% Demos use of Eyelink('GetQueuedData')
% records for a while, then plots samples, events, and execution
% times of GetQueuedData and some non-preallocated storage including
% data-conversion to take less space

% History:
% ??/??/??   Written by Erik Flister.
% 20/04/15   Removed intwarning() - no longer suppored by Matlab (MK).

format long g
clear
clear classes
close all
clc

KbName('UnifyKeyNames')

if (Eyelink('initialize') ~= 0)
    error('could not init connection to Eyelink')
    return;
end;

try
    Eyelink('Verbosity', 0);
    oldPriority = Priority(MaxPriority('KbCheck'));

    el=EyelinkInitDefaults();
    el.LOSTDATAEVENT=hex2dec('3F');

    eyelink('command', 'active_eye = RIGHT');
    
    % configure eyelink to send raw data
    status=Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER');
    if status~=0
        status
        error('link_sample_data error')
    end

    status=Eyelink('command','file_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER');
    if status~=0
        status
        error('file_sample_data error')
    end

    status=Eyelink('command','inputword_is_window = ON');
    if status~=0
        status
        error('inputword_is_window error')
    end

    edfFile='demo.edf';
    status=Eyelink('openfile',edfFile);
    if status~=0
        status
        error('openfile error')
    end

    [version, versionString]  = Eyelink('GetTrackerVersion')

    status=Eyelink('startrecording');
    if status~=0
        status
        error('startrecording error')
    end

    WaitSecs(0.1);
    status=Eyelink('message','SYNCTIME');
    if status~=0
        status
        error('message error')
    end

    screenFlipSimulationTime=0;
    maxEvents=100;
    maxSamples=2000*10*60;

    samples=nan(48,maxSamples);
    events=nan(30,maxEvents);
    sampleNum=1;
    eventNum=1;

    times=[];
    convertTimes=[];
    unsafeConvertTimes=[];
    smalls={};

    fprintf('hit spacebar to stop\n')

    stopkey=KbName('space');

    eye_used = Eyelink('EyeAvailable');

    switch eye_used
        case el.BINOCULAR
            error('tracker indicates binocular')
        case el.LEFT_EYE
            error('tracker indicates left eye')
        case el.RIGHT_EYE
            disp('tracker indicates right eye')
        case -1
            error('eyeavailable returned -1')
        otherwise
            eye_used
            error('unexpected result from eyeavailable')
    end

    recordingStartTime=GetSecs();

    while true
        err=Eyelink('checkrecording');
        if(err~=0)
            err
            error('checkrecording problem')
        end

        drained = false;
        while ~drained
            start=GetSecs;
            [samplesIn, eventsIn, drained] = Eyelink('GetQueuedData', eye_used);
            times(end+1)=GetSecs-start;

            if ~isempty(samplesIn)
                samples(:,sampleNum:sampleNum+size(samplesIn,2)-1)=samplesIn;
                sampleNum=sampleNum+size(samplesIn,2);
                
                lost=sum(samplesIn(2,:)==el.LOSTDATAEVENT);
                if lost
                    fprintf('got %d losts!\n',lost)
                end

                start=GetSecs;
                smalls{1,end+1}=samplesIn(2,:)==el.LOSTDATAEVENT;
                smalls{2,end}=uint16(samplesIn(20,:));
                smalls{3,end}=uint32(samplesIn([36:43 46],:));
                convertTimes(end+1)=GetSecs-start;

                start=GetSecs;
                smalls{1,end+1}=samplesIn(2,:)==el.LOSTDATAEVENT;
                smalls{2,end}=uint16(samplesIn(20,:));
                smalls{3,end}=uint32(samplesIn([36:43 46],:));
                unsafeConvertTimes(end+1)=GetSecs-start;
            end
            if ~isempty(eventsIn)
                events(:,eventNum:eventNum+size(eventsIn,2)-1)=eventsIn;
                eventNum=eventNum+size(eventsIn,2);
            end

            if false && (size(samples,2)>0 || size(events,2)>0)
                fprintf('got %d samples, %d events, and drained is %d\n',size(samples,2),size(events,2),drained)
            end

            if ~drained
                fprintf('got a drain failure!\n')
            end
        end

        if sampleNum>=maxSamples || eventNum>=maxEvents
            fprintf('stopping cuz filled sample or event allocation\n')
            break;
        end

        [keyIsDown,secs,keyCode] = KbCheck;
        if keyCode(stopkey)
            fprintf('stopping cuz got stopkey')
            break;
        end

        if screenFlipSimulationTime
            WaitSecs(screenFlipSimulationTime);
        end
    end
    
    figure
    subplot(3,1,1)
    plot(samples')
    title('samples')
    xlabel('sampleNum')

    subplot(3,1,2)
    plot(events')
    title('events')
    xlabel('eventNum')

    subplot(3,1,3)
    plot(times,'b');
    hold on
    plot(convertTimes,'r')
    plot(unsafeConvertTimes,'g')
    title('times')
    legend({'GetQueuedData','type conversions','unsafe type conversions'})
    xlabel('call num')

    WaitSecs(0.1);
    cleanup;
catch
    cleanup;
end

function cleanup
        chk=Eyelink('checkrecording');
        if chk~=0
            disp('problem: wasn''t recording but should have been')
        end
        Eyelink('stoprecording');
        ShowCursor;
        Priority(oldPriority);
        status=Eyelink('closefile');
        if status ~=0
            disp(sprintf('closefile error, status: %d',status))
        end
        status=Eyelink('ReceiveFile',edfFile,pwd,1);
        if status~=0
            fprintf('problem: ReceiveFile status: %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        else
            disp('unknown where data file went')
        end
        Eyelink('shutdown');
    end
end
