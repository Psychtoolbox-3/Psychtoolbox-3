function AntiSaccade(placeHolderFlag,gapManipulation)
%
% ___________________________________________________________________
%
% AntiSaccade(placeHolderFlag,gapManipulation)
% placeHolderFlag = display place holders. Default 1
% gapManipulation = time gap. Default 1
%
% Demo implemetation of prosaccade and antisaccade task.
%
% stimulus: place holders, random duration of the initial fixation, and gap manipulation before the target presentation.
%
% This task demonstrates stimuli presentation with Eyelink and Data Viewer
% integration.
%
% ___________________________________________________________________

% HISTORY
% mm/dd/yy
%
% 01/28/11  NJ  created

PsychDefaultSetup(1);

if nargin < 1
    placeHolderFlag  = true;
end
if nargin < 2
    gapManipulation = true;
end

dummymode = 0;
TARGET_TIMER = 1250;
GAP_TIMER = 200;

task = zeros(4, 32);

task(1:16,1)   = 1;   % prosccade
task(1:8,2)    = 1;   % left
task(1:4,3)    = 5;   % amplitude
task(5:8,3)   = 10;
task(9:16,2)   = 2;   % Right
task(9:12,3)   = 5;   % amplitude
task(13:16,3) = 10;

task(17:32,1)  = 2;   % antisaccade
task(17:24,2)  = 1;   % left
task(17:20,3)  = 5;   % amplitude
task(21:24,3) = 10;
task(25:32,2)  = 2;   % Right
task(25:28,3)  = 5;   % amplitude
task(29:32,3) = 10;

saccade = '';
feedbackpos = '';

instructionText = [' In this experiment, you are going to see a series of events happening in each trial.\n ' ...
    ' A colored fixation box will be presented at the beginning of the trial ,\n'...
    ' followed by a white box appearing on either side of the screen.\n\n' ...
    '* If the color of the initial fixation box is green,\n'...
    '     you should look at the white box.\n' ...
    '* If the color of the initial fixation box is red, \n'...
    '     you should look in the opposite direction from the white box with same amount of\n'...
    '     distance to the center of the screen.\n'...
    '* Please try to perform the task as quickly and as precisely as possible.\n\n ' ...
    'Press any key to continue.'];


% Pixels Per Degree in X and Y
PPD_X = 30;
PPD_Y = 30;

if ~isoctave
    commandwindow;
else
    more off;
end

try
    
    %%%%%%%%%%
    % STEP 1 %
    %%%%%%%%%%
    
    % Added a dialog box to set your own EDF file name before opening
    % experiment graphics. Make sure the entered EDF file name is 1 to 8
    % characters in length and only numbers or letters are allowed.
    % Note: Octave does not support GUIs. replace lines below with
    % %edfFile= 'DEMO.EDF'
    
    if isoctave
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
    [window, wRect]=Screen('OpenWindow', screenNumber, 0,[],32,2); %#ok<*NASGU>
    Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [winWidth, winHeight] = WindowSize(window);
    
    % Select specific text font, style and size:
    Screen('TextFont',window, 'Courier New');
    Screen('TextSize',window, round(winWidth * 0.015));
    Screen('TextStyle', window, 1+2);
    
    DrawFormattedText(window, instructionText, 100,200, [255 255 255]);
    Screen('Flip',window);
    
    KbWait;
    while KbCheck; end;
    
    
    %%%%%%%%%%
    % STEP 3 %
    %%%%%%%%%%
    
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    % make necessary changes to calibration structure parameters and pass
    % it to EyelinkUpdateDefaults for changes to take affect
    
    el=EyelinkInitDefaults(window);
    
    % We are changing calibration to a black background with white targets,
    % no sound and smaller targets
    el.backgroundcolour = BlackIndex(el.window);
    el.msgfontcolour  = WhiteIndex(el.window);
    el.imgtitlecolour = WhiteIndex(el.window);
    el.targetbeep = 0;
    el.calibrationtargetcolour= WhiteIndex(el.window);
    % for lower resolutions you might have to play around with these values
    % a little. If you would like to draw larger targets on lower res
    % settings please edit PsychEyelinkDispatchCallback.m and see comments
    % in the EyelinkDrawCalibrationTarget function
    el.calibrationtargetsize= 1;
    el.calibrationtargetwidth=0.5;
    
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
    res = Eyelink('Openfile', edfFile);
    if res~=0
        fprintf('Cannot create EDF file ''%s'' ', edffilename);
        cleanup;
        return;
    end
    
    % make sure we're still connected.
    if Eyelink('IsConnected')~=1 && ~dummymode
        cleanup;
        return;
    end
    
    
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
    
    if v == 3 && str2double(vsn{1}) == 4 % if EL 1000 and tracker version 4.xx
        
        % remote mode possible add HTARGET ( head target)
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
    else
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    % allow to use the big button on the eyelink gamepad to accept the
    % calibration/drift correction target
    Eyelink('command', 'button_function 5 "accept_target_fixation"');
    
    %%%%%%%%%%
    % STEP 6 %
    %%%%%%%%%%
    
    % Hide the mouse cursor and Calibrate the eye tracker
    
    Screen('HideCursorHelper', window);
    % enter Eyetracker camera setup mode, calibration and validation
    EyelinkDoTrackerSetup(el);
    
    
    % setup a random permutation of the trials
    order = randperm(32);
    dots = zeros(2,3);
    dots(2,:) = winHeight/2;
    dots(1,1) = winWidth/2;
    color = zeros(3,3);
    if placeHolderFlag
        color(:,2:3) = [255 255 255 ; 255 255 255]';
    else
        color(:,2:3) = [0 0 0 ; 0 0 0]';
    end
    feedback = [-15 -15 15 15];
    
    %%%%%%%%%%
    % STEP 7 %
    %%%%%%%%%%
    
    % Now starts running individual trials;
    % You can keep the rest of the code except for the implementation
    % of graphics and event monitoring
    % Each trial should have a pair of "StartRecording" and "StopRecording"
    % calls as well integration messages to the data file (message to mark
    % the time of critical events and the image/interest area/condition
    % information for the trial)
    
    
    for i=1:32
        
        % determine current trial type and send msg to edf
        perm = task(order(i),:);
        Eyelink('Message', 'permutation number %d', order(i));
        if perm(1) == 1
            Eyelink('Message', 'Prosaccade');
            color(:,1) = [0 255 0]';
        else
            Eyelink('Message', 'Antisaccade');
            color(:,1) = [255 0 0]';
        end
        if perm(2) == 1
            feedbackpos = 'left';
            Eyelink('Message', 'Left');
        else
            feedbackpos = 'right';
            Eyelink('Message', 'Right');
        end
        Eyelink('Message', 'Amplitude %d', perm(3));
        
        dots(1,2) = winWidth/2 + perm(3) * PPD_X; % right
        dots(1,3) = winWidth/2 - perm(3) * PPD_X; % left
        
        if perm(2) == 1
            feedback = CenterRect(feedback, [dots(1,3)-5 dots(2,3)-5     dots(1,3)+5 dots(2,3)+5]);
        else
            feedback = CenterRect(feedback, [dots(1,2)-5 dots(2,2)-5     dots(1,2)+5 dots(2,2)+5]);
        end
        
        % STEP 7.1
        % Sending a 'TRIALID' message to mark the start of a trial in Data
        % Viewer.  This is different than the start of recording message
        % START that is logged when the trial recording begins. The viewer
        % will not parse any messages, events, or samples, that exist in
        % the data file prior to this message.
        Eyelink('Message', 'TRIALID %d', i);
        
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d/%d"', i, 32);
        % Before recording, we place reference graphics on the host display
        % Must be offline to draw to EyeLink screen
        Eyelink('Command', 'set_idle_mode');
        % clear tracker display and draw box at center
        Eyelink('Command', 'clear_screen 0');
        
        
        if perm(1) == 1
            % prosaccade
            
            saccade = 'pro';
            
            if perm(2) == 1
                rect1 = [dots(1,3)-60, dots(2,3)-60, dots(1,3)+60, dots(2,3) + 60];
                rect2 =[dots(1,2)-60, dots(2,2)-60, dots(1,2)+60, dots(2,2) + 60];
            else
                rect1 = [dots(1,2)-60, dots(2,2)-60, dots(1,2)+60, dots(2,2) + 60];
                rect2 =[dots(1,3)-60, dots(2,3)-60, dots(1,3)+60, dots(2,3) + 60];
            end
        else
            % anti saccade
            
            saccade = 'anti';
            
            if perm(2) == 2
                rect1 = [dots(1,3)-60, dots(2,3)-60, dots(1,3)+60, dots(2,3) + 60];
                rect2 =[dots(1,2)-60, dots(2,2)-60, dots(1,2)+60, dots(2,2) + 60];
            else
                rect1 = [dots(1,2)-60, dots(2,2)-60, dots(1,2)+60, dots(2,2) + 60];
                rect2 =[dots(1,3)-60, dots(2,3)-60, dots(1,3)+60, dots(2,3) + 60];
            end
        end
        
        % draw shapes to host pc
        Eyelink('command', 'draw_filled_box %d %d %d %d 2' ,rect1(1),rect1(2),rect1(3),rect1(4));
        Eyelink('command', 'draw_filled_box %d %d %d %d 4' ,rect2(1),rect2(2),rect2(3),rect2(4));
        
        % STEP 7.2
        % Do a drift correction at the beginning of each trial
        % Performing drift correction (checking) is optional for
        % EyeLink 1000 eye trackers.
        EyelinkDoDriftCorrection(el);
        
        % STEP 7.3
        % start recording eye position (preceded by a short pause so that
        % the tracker can finish the mode transition)
        % The paramerters for the 'StartRecording' call controls the
        % file_samples, file_events, link_samples, link_events availability
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.05);
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        % otherwise you may lose a few msec of data
        WaitSecs(0.1);
        
        % random fixation time between 800 and 1200 ms
        fixationTime = GetSecs + ((800 + (1200-800) * rand)/1000);
        while GetSecs < fixationTime
            
            % STEP 7.4
            % Prepare and show the screen
            Screen('FillRect', window, el.backgroundcolour);
            Screen('DrawDots',window, dots,10, color);
            Screen('Flip', window);
            Eyelink('Message', 'FIXATION_DISPLAY');
        end
        
        if gapManipulation
            gapTime = GetSecs + GAP_TIMER/1000;
            while GetSecs < gapTime
                Screen('FillRect', window, el.backgroundcolour);
                Screen('Flip', window);
                Eyelink('Message', 'GAP_DISPLAY');
            end
        end
        
        targetTime = GetSecs + TARGET_TIMER/1000;
        firstloop = 1;
        imgfile = sprintf('img%d%d.jpg' , perm(3),perm(2));
        while GetSecs < targetTime
            Screen('FillRect', window, el.backgroundcolour);
            Screen('DrawDots',window, dots(:,2:3),10, color(:,2:3) );
            Screen('FrameRect',window,255, feedback);
            Screen('Flip', window);
            % mark zero-plot time in data file
            Eyelink('Message', 'TARGET_DISPLAY');
            % this is done so that you do not write images everytime you might
            % want to delete the image files if you change the trial
            % parameters. Another option is hand drawing iA in Data Viewer
            if firstloop && ~exist(imgfile,'file')
                imageArray = Screen('GetImage', window);
                imwrite(imageArray,imgfile);
                firstloop = 0;
            end
        end
        
        
        % STEP 7.6
        % add 100 msec of data to catch final events and blank display
        WaitSecs(0.1);
        Eyelink('StopRecording');
        Screen('FillRect', window, el.backgroundcolour);
        Screen('Flip', window);
        
        % STEP 7.7
        % Send out necessary integration messages for data analysis
        % See "Protocol for EyeLink Data to Viewer Integration-> Interest
        % Area Commands" section of the EyeLink Data Viewer User Manual
        % IMPORTANT! Don't send too many messages in a very short period of
        % time or the EyeLink tracker may not be able to write them all
        % to the EDF file.
        % Consider adding a short delay every few messages.
        WaitSecs(0.001);
        % Send an integration message so that an image can be loaded as
        % overlay backgound when performing Data Viewer analysis.  This
        % message can be placed anywhere within the scope of a trial (i.e.,
        % after the 'TRIALID' message and before 'TRIAL_RESULT')
        % See "Protocol for EyeLink Data to Viewer Integration -> Image
        % Commands" section of the EyeLink Data Viewer User Manual.
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', imgfile,  winWidth/2, winHeight/2);
        
        % interest areas
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, rect1(1), rect1(2), rect1(3), rect1(4),'target');
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 2, rect2(1), rect2(2), rect2(3), rect2(4),'distractor');
        % Send messages to report trial condition information
        % Each message may be a pair of trial condition variable and its
        % corresponding value follwing the '!V TRIAL_VAR' token message
        % See "Protocol for EyeLink Data to Viewer Integration-> Trial
        % Message Commands" section of the EyeLink Data Viewer User Manual
        WaitSecs(0.001);
        
        Eyelink('Message', '!V TRIAL_VAR index %d', i);
        Eyelink('Message', '!V TRIAL_VAR amplitude %d', perm(3)); % 5 or 10
        Eyelink('Message', '!V TRIAL_VAR saccade %s', saccade);   % pro or anti
        Eyelink('Message', '!V TRIAL_VAR feedback %s', feedbackpos); % left or right?
        % STEP 7.8
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
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    
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
    psychrethrow(lasterror); %#ok<LERR>
end

    function cleanup
        % Shutdown Eyelink:
        Eyelink('Shutdown');
        Screen('CloseAll');
    end

end
