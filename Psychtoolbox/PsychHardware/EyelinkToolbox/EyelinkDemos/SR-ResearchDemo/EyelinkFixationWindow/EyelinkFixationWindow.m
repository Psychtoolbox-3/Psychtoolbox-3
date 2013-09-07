function EyelinkFixationWindow
%
% ___________________________________________________________________
%
% Demo implemetation of detecting eye in fixation window on display
% computer through samples from the tracker
%
% stimulus: A fixation dot
%
% This task demonstrates stimuli presentation with Eyelink and Data Viewer
% integration.
%
% press space to end experiment
% ___________________________________________________________________

% HISTORY
% mm/dd/yy
%
% 01/28/11  NJ  created

PsychDefaultSetup(1);

% trial defaults
dummymode=0;

if ~IsOctave
    commandwindow;
else
    more off;
end
stopkey=KbName('space');
firstRun = 1;
infix=0;
dotSize = 10;
fixWinSize = 100;
imageList = {'HappyFace.jpg' 'SadFace.jpg'};

%[versionString, versionStructure]=PsychtoolboxVersion;

try
    
    %%%%%%%%%%
    % STEP 1 %
    %%%%%%%%%%
    
    % Added a dialog box to set your own EDF file name before opening
    % experiment graphics. Make sure the entered EDF file name is 1 to 8
    % characters in length and only numbers or letters are allowed.
    if IsOctave
        edfFile = 'DEMO';
    else 
        prompt = {'Enter tracker EDF file name (1 to 8 letters or numbers)'};
        dlg_title = 'Create EDF file';
        num_lines= 1;
        def     = {'DEMO'};
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        edfFile = answer{1};
        fprintf('EDFFile: %s\n', edfFile );
    end
    %%%%%%%%%%
    % STEP 2 %
    %%%%%%%%%%
    
    % Open a graphics window on the main screen
    % using the PsychToolbox's Screen function.
    screenNumber=max(Screen('Screens'));
    [window, wRect]=Screen('OpenWindow', screenNumber, 0,[],32,2);
    Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [winWidth, winHeight] = WindowSize(window);
    
    % create fixation dot and fixation window rectangles
    fixationDot = [-dotSize -dotSize dotSize dotSize];
    fixationDot = CenterRect(fixationDot, wRect);    
    fixationWindow = [-fixWinSize -fixWinSize fixWinSize fixWinSize];
    fixationWindow = CenterRect(fixationWindow, wRect);
    
    % make textures 
    happy = Screen('MakeTexture',window, imread(char(imageList(1))));
    sad = Screen('MakeTexture',window, imread(char(imageList(2))));
    
    %%%%%%%%%%
    % STEP 3 %
    %%%%%%%%%%
    
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    
    el=EyelinkInitDefaults(window);
    
    % We are changing calibration to a black background with white targets,
    % no sound and smaller targets
    el.backgroundcolour = BlackIndex(el.window);
    el.msgfontcolour  = WhiteIndex(el.window);
    el.imgtitlecolour = WhiteIndex(el.window);
    el.targetbeep = 0;
    el.calibrationtargetcolour = WhiteIndex(el.window);
    
    % for lower resolutions you might have to play around with these values
    % a little. If you would like to draw larger targets on lower res
    % settings please edit PsychEyelinkDispatchCallback.m and see comments
    % in the EyelinkDrawCalibrationTarget function
    el.calibrationtargetsize= 1;
    el.calibrationtargetwidth=0.5;
    % call this function for changes to the calibration structure to take
    % affect
    EyelinkUpdateDefaults(el);
    
    %%%%%%%%%%
    % STEP 4 %
    %%%%%%%%%%
    
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end
      
    % open file to record data to
    i = Eyelink('Openfile', edfFile);
    if i~=0
        fprintf('Cannot create EDF file ''%s'' ', edffilename);
        cleanup;
        return;
    end
    
    % make sure we're still connected.
    if Eyelink('IsConnected')~=1 && ~dummymode
        cleanup;
        return;
    end;
        
    %%%%%%%%%%
    % STEP 5 %
    %%%%%%%%%%
    
    % SET UP TRACKER CONFIGURATION
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    Eyelink('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox demo-experiment''');
   
    % This command is crucial to map the gaze positions from the tracker to
    % screen pixel positions to determine fixation
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, winWidth-1, winHeight-1);
    
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, winWidth-1, winHeight-1);
    % set calibration type.
    Eyelink('command', 'calibration_type = HV9');
    Eyelink('command', 'generate_default_targets = YES');
    % set parser (conservative saccade thresholds)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    % set EDF file contents
        % 5.1 retrieve tracker version and tracker software version
    [v,vs] = Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    vsn = regexp(vs,'\d','match');
    
    if v ==3 && str2double(vsn{1}) == 4 % if EL 1000 and tracker version 4.xx
        
        % remote mode possible add HTARGET ( head target)
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
    else
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    
    % calibration/drift correction target
    Eyelink('command', 'button_function 5 "accept_target_fixation"');
    
    %%%%%%%%%%
    % STEP 6 %
    %%%%%%%%%%
    
    if ~dummymode
        % Hide the mouse cursor and setup the eye calibration window
        Screen('HideCursorHelper', window);
    end
    % enter Eyetracker camera setup mode, calibration and validation
    EyelinkDoTrackerSetup(el);
    
    
    %%%%%%%%%
    % STEP 7%
    %%%%%%%%%
    
    % Now starts running individual trials
    % You can keep the rest of the code except for the implementation
    % of graphics and event monitoring
    % Each trial should have a pair of "StartRecording" and "StopRecording"
    % calls as well integration messages to the data file (message to mark
    % the time of critical events and the image/interest area/condition
    % information for the trial)
    trial = 1;
    numTrials = 3;
    index = 1;
    
    % repeat until we have 3 sucessful trials
    while trial <= numTrials
        
        % wait a second between trials
        WaitSecs(1);
        % STEP 7.1
        % Sending a 'TRIALID' message to mark the start of a trial in Data
        % Viewer.  This is different than the start of recording message
        % START that is logged when the trial recording begins. The viewer
        % will not parse any messages, events, or samples, that exist in
        % the data file prior to this message.       
        Eyelink('Message', 'TRIALID %d', trial);
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d/%d"', trial,numTrials);
        Eyelink('Command', 'set_idle_mode');
        % clear tracker display and draw box at center
        Eyelink('Command', 'clear_screen %d', 0);
        % draw fixation and fixation window shapes on host PC
        Eyelink('command', 'draw_cross %d %d 15', winWidth/2,winHeight/2);
        Eyelink('command', 'draw_box %d %d %d %d 15', fixationWindow(1), fixationWindow(2), fixationWindow(3), fixationWindow(4));
            
        % STEP 7.2
        % Do a drift correction at the beginning of each trial
        % Performing drift correction (checking) is optional for
        % EyeLink 1000 eye trackers. Drift correcting at different
        % locations x and y depending on where the ball will start
        % we change the location of the drift correction to match that of
        % the target start position
        EyelinkDoDriftCorrection(el);
        
        % STEP 7.3
        % start recording eye position (preceded by a short pause so that
        % the tracker can finish the mode transition)
        % The paramerters for the 'StartRecording' call controls the
        % file_samples, file_events, link_samples, link_events availability
       
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.05);
        Eyelink('StartRecording');
        eye_used = Eyelink('EyeAvailable'); % get eye that's tracked  
        % returns 0 (LEFT_EYE), 1 (RIGHT_EYE) or 2 (BINOCULAR) depending on what data is
        if eye_used == 2
            eye_used = 1; % use the right_eye data
        end
        % record a few samples before we actually start displaying
        % otherwise you may lose a few msec of data
        WaitSecs(0.1);
        
        % STEP 7.4
        % Prepare and show the screen.
        Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('FillRect', window, el.backgroundcolour);
        Screen('FillOval', window,round([rand*255 rand*255 rand*255]), fixationDot);
        Screen('Flip',window);
        Eyelink('Message', 'SYNCTIME');
        
        % get screen image from the first display to use as Data Viewer
        % trial overlay. Note this call is very slow and will affect your
        % timing for the first screen blanking
        if firstRun
            imageArray = Screen('GetImage', window);
            firstRun =0;
        end
                
        % set fixation display to be randomly chose between 650 and 1500
        fixateTime = GetSecs + round(650 + (1500-650).*rand)/1000 + 200/1000;
        graceTime = GetSecs + 200/1000;
        while GetSecs < fixateTime
            
            if dummymode==0
                error=Eyelink('CheckRecording');
                if(error~=0)
                    break;
                end
                
                if Eyelink( 'NewFloatSampleAvailable') > 0
                    % get the sample in the form of an event structure
                    evt = Eyelink( 'NewestFloatSample');
                    evt.gx
                    evt.gy
                    if eye_used ~= -1 % do we know which eye to use yet?
                        % if we do, get current gaze position from sample
                        x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                        y = evt.gy(eye_used+1);
                        % do we have valid data and is the pupil visible?
                        if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                            mx=x;
                            my=y;
                        end
                    end
                end
            else
                
                % Query current mouse cursor position (our "pseudo-eyetracker") -
                % (mx,my) is our gaze position.
                [mx, my]=GetMouse(window); %#ok<*NASGU>
           
            end
            if infixationWindow(mx,my) && ~infix
                
                Eyelink('Message', 'Fixation Start');
                Beeper(el.calibration_success_beep(1), el.calibration_success_beep(2), el.calibration_success_beep(3));
                infix = 1;
            elseif ~infixationWindow(mx,my) && infix && GetSecs > graceTime
                
                Screen('DrawTexture', window, sad);
                Screen('Flip',window);
                disp('broke fix');
                Eyelink('Message', 'Fixation broke or grace time ended');
                Beeper(el.calibration_failed_beep(1), el.calibration_failed_beep(2), el.calibration_failed_beep(3));
                infix = 0;
                break;
            end
            
            [keyIsDown,secs,keyCode] = KbCheck; %#ok<*ASGLU>
            % if spacebar was pressed stop display
            if keyCode(stopkey )
                sprintf('Space pressed, exiting trial\n');
                EyeLink('Message', 'Key pressed');
                break;
            end
            
        end
        
        if infix
            Screen('DrawTexture', window, happy);
            Screen('Flip',window);
            sprintf('Trial completed. Trial %d of %d\n', trial, numTrials);
            trial = trial + 1;
            WaitSecs(1);
        end        
       
        
        % STEP 7.5
        % add 100 msec of data to catch final events and blank display
        WaitSecs(0.1);
        Eyelink('StopRecording');

        index = index + 1;
        
        Screen('FillRect', window, el.backgroundcolour);
        Screen('Flip', window);
        
        imwrite(imageArray, 'imgfile.jpg', 'jpg');
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', 'imgfile.jpg', winWidth/2, winHeight/2);
        
        % STEP 7.6    
        % Send out necessary integration messages for data analysis
        % Send out interest area information for the trial
        % See "Protocol for EyeLink Data to Viewer Integration-> Interest
        % Area Commands" section of the EyeLink Data Viewer User Manual
        % IMPORTANT! Don't send too many messages in a very short period of
        % time or the EyeLink tracker may not be ablwWe to write them all
        % to the EDF file.
        % Consider adding a short delay every few messages.
        WaitSecs(0.001);
        Eyelink('Message', '!V IAREA ELLIPSE %d %d %d %d %d %s', 1, winWidth/2-dotSize, winHeight/2-dotSize, winWidth/2+dotSize, winHeight/2+dotSize,'center');
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 2, winWidth/2-fixWinSize, winHeight/2-fixWinSize, winWidth/2+fixWinSize, winHeight/2+fixWinSize,'centerWin');
        
        
        % Send messages to report trial condition information
        % Each message may be a pair of trial condition variable and its
        % corresponding value follwing the '!V TRIAL_VAR' token message
        % See "Protocol for EyeLink Data to Viewer Integration-> Trial
        % Message Commands" section of the EyeLink Data Viewer User Manual
        WaitSecs(0.001);
        Eyelink('Message', '!V TRIAL_VAR index %d', index);
        Eyelink('Message', '!V TRIAL_VAR imgfile %s', 'imgfile.jpg');
        if infix
            Eyelink('Message', '!V TRIAL_VAR trialOutcome %s', 'succesful');
        else
            Eyelink('Message', '!V TRIAL_VAR trialOutcome %s', 'recycled');
        end
        % STEP 9
        % Sending a 'TRIAL_RESULT' message to mark the end of a trial in
        % Data Viewer. This is different than the end of recording message
        % END that is logged when the trial recording ends. The viewer will
        % not parse any messages, events, or samples that exist in the data
        % file after this message.
        Eyelink('Message', 'TRIAL_RESULT 0');
        
    end
    
     %%%%%%%%%%
     % STEP 8 %
     %%%%%%%%%%
 
    % End of Experiment; close the file first
    % close graphics window, close data file and shut down tracker
    Screen('CloseAll');
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    
    % download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch %#ok<*CTCH>
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end
    
   
    %%%%%%%%%%
    % STEP 9 %
    %%%%%%%%%%
    
    % run cleanup function (close the eye tracker and window).
    cleanup;
    
catch
    cleanup;
    fprintf('%s: some error occured\n', mfilename);
    psychrethrow(lasterror); %#ok<*LERR>
    
end

    function cleanup
        % Shutdown Eyelink:
        Eyelink('Shutdown');
        Screen('CloseAll');
    end

    function fix = infixationWindow(mx,my)
        % determine if gx and gy are within fixation window
        fix = mx > fixationWindow(1) &&  mx <  fixationWindow(3) && ...
            my > fixationWindow(2) && my < fixationWindow(4) ;
    end
end


