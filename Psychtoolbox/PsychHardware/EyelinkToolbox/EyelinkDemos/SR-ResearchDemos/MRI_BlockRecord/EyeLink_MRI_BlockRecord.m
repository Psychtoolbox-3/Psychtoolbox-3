function EyeLink_MRI_BlockRecord(screenNumber)
% Simple MRI demo with EyeLink integration.
% 6 trials are presented in 2 blocks of 3 trials. Trial duration is 5.5s during which a 4s stimulus is presented.
% A block starts with a drift-check followed by presentation of central crosshairs. Eye movements are recorded while 
% waiting for an MRI trigger (keyboard key 't' in this demo). The stimulus is presented when trigger is received.
% A fixed ITI is maintained by presenting crosshairs between each 4s stimulus. Eye movements are recorded throughout
% an entire block rather than on a trial-by-trial basis. 
%
% In STEP 5 it is shown how to:
% - shrink the spread of the calibration/validation targets so they are all visible if the MRI bore blocks part of the screen
% - apply an optional online drift correction (see EyeLink 1000 Plus User Manual section 3.11.2)
%
% Usage:
% Eyelink_MRI_BlockRecord(screenNumber)
%
% screenNumber is an optional parameter which can be used to pass a specific value to PsychImaging('OpenWindow', ...)
% If screenNumber is not specified, or if isempty(screenNumber) then the default:
% screenNumber = max(Screen('Screens'));
% will be used.

% Initialize PsychSound for calibration/validation audio feedback
InitializePsychSound();

%% STEP 1: PROVIDE SOME SESSION PARAMETERS

stimDur = 4.0; % stimulus duration in seconds
trialDur = 5.5; % trial duration in seconds

% Bring the Command Window to the front if it is already open
if ~IsOctave; commandwindow; end
% Use default screenNumber if none specified
if (nargin < 1)
    screenNumber = [];
end

try
    %% STEP 2: INITIALIZE EYELINK CONNECTION; OPEN EDF FILE; GET EYELINK TRACKER VERSION
    
    % Initialize EyeLink connection (dummymode = 0) or run in "Dummy Mode" without an EyeLink connection (dummymode = 1);
    dummymode = 0;
    EyelinkInit(dummymode); % Initialize EyeLink connection
    status = Eyelink('IsConnected');
    if status < 1 % If EyeLink not connected
        dummymode = 1; 
    end
    
    % Open dialog box for EyeLink Data file name entry. File name up to 8 characters
    prompt = {'Enter EDF file name (up to 8 characters)'};
    dlg_title = 'Create EDF file';
    def = {'demo'}; % Create a default edf file name
    answer = inputdlg(prompt, dlg_title, 1, def); % Prompt for new EDF file name    
    % Print some text in Matlab's Command Window if a file name has not been entered
    if  isempty(answer)
        fprintf('Session cancelled by user\n')
        cleanup; % Abort experiment (see cleanup function below)
        return
    end    
    edfFile = answer{1}; % Save file name to a variable    
    % Print some text in Matlab's Command Window if file name is longer than 8 characters
    if length(edfFile) > 8
        fprintf('Filename needs to be no more than 8 characters long (letters, numbers and underscores only)\n');
        cleanup; % Abort experiment (see cleanup function below)
        return
    end
    
    % Open an EDF file and name it
    failOpen = Eyelink('OpenFile', edfFile);
    if failOpen ~= 0 % Abort if it fails to open
        fprintf('Cannot create EDF file %s', edfFile); % Print some text in Matlab's Command Window
        cleanup; %see cleanup function below
        return
    end
    
    % Get EyeLink tracker and software version
    % <ver> returns 0 if not connected
    % <versionstring> returns 'EYELINK I', 'EYELINK II x.xx', 'EYELINK CL x.xx' where 'x.xx' is the software version
    ELsoftwareVersion = 0; % Default EyeLink version in dummy mode
    [ver, versionstring] = Eyelink('GetTrackerVersion');
    if dummymode == 0 % If connected to EyeLink
        % Extract software version number. 
        [~, vnumcell] = regexp(versionstring,'.*?(\d)\.\d*?','Match','Tokens'); % Extract EL version before decimal point
        ELsoftwareVersion = str2double(vnumcell{1}{1}); % Returns 1 for EyeLink I, 2 for EyeLink II, 3/4 for EyeLink 1K, 5 for EyeLink 1KPlus, 6 for Portable Duo              
        % Print some text in Matlab's Command Window
        fprintf('Running experiment on %s version %d\n', versionstring, ver );
    end
    % Add a line of text in the EDF file to identify the current experimemt name and session. This is optional.
    % If your text starts with "RECORDED BY " it will be available in DataViewer's Inspector window by clicking
    % the EDF session node in the top panel and looking for the "Recorded By:" field in the bottom panel of the Inspector.
    preambleText = sprintf('RECORDED BY Psychtoolbox demo %s session name: %s', mfilename, edfFile);
    Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);
    
    
    %% STEP 3: SELECT AVAILABLE SAMPLE/EVENT DATA
    % See EyeLinkProgrammers Guide manual > Useful EyeLink Commands > File Data Control & Link Data Control
    
    % Select which events are saved in the EDF file. Include everything just in case
    Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    % Select which events are available online for gaze-contingent experiments. Include everything just in case
    Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON,FIXUPDATE,INPUT');
    % Select which sample data is saved in EDF file or available online. Include everything just in case
    if ELsoftwareVersion > 3  % Check tracker version and include 'HTARGET' to save head target sticker data for supported eye trackers
        Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT');
        Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
    else
        Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT');
        Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    
    
    %% STEP 4: OPEN GRAPHICS WINDOW
    
    % Open experiment graphics on the specified screen
    if isempty(screenNumber)
        screenNumber = max(Screen('Screens')); % Use default screen if none specified
    end
    PsychDefaultSetup(2);
    window = PsychImaging('OpenWindow', screenNumber, GrayIndex(screenNumber)); % Open graphics window
    Screen('Flip', window);
    
    ifi = Screen('GetFlipInterval', window); % Return an estimate of the monitor flip interval 

    % Get max color value for rescaling  to RGB for Host PC & Data Viewer integration
    colorMaxVal = Screen('ColorRange', window);
    % Return width and height of the graphics window/screen in pixels
    [width, height] = Screen('WindowSize', window);
    
    
    %% STEP 5: SET CALIBRATION SCREEN COLOURS; PROVIDE WINDOW SIZE TO EYELINK HOST & DATAVIEWER; SET CALIBRATION PARAMETERS; CALIBRATE
    
    % Provide EyeLink with some defaults, which are returned in the structure "el".
    el = EyelinkInitDefaults(window);
    % set calibration/validation/drift-check(or drift-correct) size as well as background and target colors. 
    % It is important that this background colour is similar to that of the stimuli to prevent large luminance-based 
    % pupil size changes (which can cause a drift in the eye movement data)
    el.calibrationtargetsize = 3;% Outer target size as percentage of the screen
    el.calibrationtargetwidth = 0.7;% Inner target size as percentage of the screen
    el.backgroundcolour = repmat(GrayIndex(screenNumber),1,3);
    el.calibrationtargetcolour = repmat(BlackIndex(screenNumber),1,3);
    % set "Camera Setup" instructions text colour so it is different from background colour
    el.msgfontcolour = repmat(BlackIndex(screenNumber),1,3);
    % You must call this function to apply the changes made to the el structure above
    EyelinkUpdateDefaults(el);
    
    % Set display coordinates for EyeLink data by entering left, top, right and bottom coordinates in screen pixels
    Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
    % Write DISPLAY_COORDS message to EDF file: sets display coordinates in DataViewer
    % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Pre-trial Message Commands
    Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);    
    % Set number of calibration/validation dots and spread: horizontal-only(H) or horizontal-vertical(HV) as H3, HV3, HV5, HV9 or HV13
    Eyelink('Command', 'calibration_type = HV9'); % horizontal-vertical 9-points
    
    % Optional: shrink the spread of the calibration/validation targets <x, y display proportion>
    % if default outermost targets are not all visible in the bore.
    % Default spread is 0.88, 0.83 (88% of the display horizontally and 83% vertically)
    Eyelink('command', 'calibration_area_proportion 0.88 0.83');
    Eyelink('command', 'validation_area_proportion 0.88 0.83');
    
    %-------------------------------------------------------------------------------------------
    % Optional: online drift correction. See section 3.11.2 in the EyeLink 1000 / EyeLink 1000 Plus User Manual    
    %   % Online drift correction to mouse-click position:
    %     Eyelink('Command', 'driftcorrect_cr_disable = OFF');
    %     Eyelink('Command', 'normal_click_dcorr = ON');
    %
    %   % Online drift correction to a fixed location:
    %     Eyelink('Command', 'driftcorrect_cr_disable = OFF');
    %     Eyelink('Command', 'online_dcorr_refposn 512,384');
    %     Eyelink('Command', 'online_dcorr_button = ON');
    %     Eyelink('Command', 'normal_click_dcorr = OFF');
    %-------------------------------------------------------------------------------------------
    
    % Hide mouse cursor
    HideCursor(screenNumber);
    % Start listening for keyboard input. Suppress keypresses to Matlab windows.
    ListenChar(-1);
    Eyelink('Command', 'clear_screen 0'); % Clear Host PC display from any previus drawing
    % Put EyeLink Host PC in Camera Setup mode for participant setup/calibration
    EyelinkDoTrackerSetup(el);
    
    
    %% STEP 6: START BLOCK LOOP; DRAW FEEDBACK GRAPHICS/TEXT ON HOST PC; DRIFT-CHECK; START RECORDING; DRAW CROSSHAIRS ON SCREEN 
    
    % Provide list of images for each block (3 images, 2 blocks)
    imgList = {{'img1.jpg', 'img2.jpg', 'img3.jpg'}, {'img1.jpg', 'img2.jpg', 'img3.jpg'}};
    trialCount = 0; % Trial counter
    for iBlock = 1:length(imgList)
        
        % Put tracker in idle/offline mode before Host PC feedback graphics drawing
        Eyelink('SetOfflineMode');
        
        % Optional: draw feedback graphics on Host PC interface
        % See section 25.7 'Drawing Commands' in the EyeLink Programmers Guide manual
        imgSize = [1024, 768]; % Use stimulus image size for drawing feedback graphics
        Eyelink('Command', 'clear_screen 0'); %Clear Host screen to black
        Eyelink('Command', 'draw_box %d %d %d %d 15', round(width/2-imgSize(1)/2), round(height/2-imgSize(2)/2), round(width/2+imgSize(1)/2), round(height/2+imgSize(2)/2));
        Eyelink('Command', 'draw_line %d %d %d %d 15', width/2, 1, width/2, height);
        Eyelink('Command', 'draw_line %d %d %d %d 15', 1, height/2, width, height/2);
        % Supply the block number as a line of text on Host PC screen
        Eyelink('Command', 'record_status_message "BLOCK %d/%d"', iBlock, length(imgList));
        
        % Perform a drift check/correction. EyeLink 1000 and 1000 Plus perform a drift-check by default
        % Optionally provide x y target location, otherwise target is presented at screen centre
        EyelinkDoDriftCorrection(el, round(width/2), round(height/2));
        
        % Write TRIALID message to EDF file: marks the start of first trial for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial.
        % TRIALID before StartRecording prevents extra initial trial in DataViewer when recording continuously
        trialCount = trialCount + 1; % Add 1 to trial counter
        Eyelink('Message', 'TRIALID %d', trialCount);
        
        % Put tracker in idle/offline mode before recording. Eyelink('SetOfflineMode') is recommended 
        % however if Eyelink('Command', 'set_idle_mode') is used allow 50ms before recording as shown in the commented code:        
        % Eyelink('Command', 'set_idle_mode');% Put tracker in idle/offline mode before recording
        % WaitSecs(0.05); % Allow some time for transition           
        Eyelink('SetOfflineMode');% Put tracker in idle/offline mode before recording
        Eyelink('StartRecording'); % Start tracker recording
        WaitSecs(0.1); % Allow some time to record a few samples before presenting first stimulus
        
        % Draw crosshairs centrally on the screen and present some text
        Screen('FillRect', window, el.backgroundcolour);
        Screen('DrawLine', window, 0, round(width/2-20), round(height/2), round(width/2+20), round(height/2), 5);
        Screen('DrawLine', window, 0, round(width/2), round(height/2-20), round(width/2), round(height/2+20), 5);
        Screen('DrawText', window, 'Waiting for MRI trigger (t)', 5, height-35, 0); % Prepare text
        Screen('Flip', window);
        
        
        %% STEP 7: TRIAL LOOP
        
        for i = 1:length(imgList{iBlock})
            
            % STEP 7.1: PREBUILD STIMULUS AND SAVE .BMP FOR DATAVIEWER
            
            % Prepare grey background on backbuffer
            Screen('FillRect', window, el.backgroundcolour);
            % Use 'drawBuffer' to copy unprocessed backbuffer images without additional processing. Prevents image size info issues on Retina displays
            backgroundArray = Screen('GetImage', window, [], 'drawBuffer'); % Copy unprocessed backbuffer
            backgroundTexture = Screen('MakeTexture', window, backgroundArray); % Convert background to texture so it is ready for drawing later on            
            % Prepare image on backbuffer
            imgName = char(imgList{iBlock}{i});% Get image file name for current trial
            imgInfo = imfinfo(imgName); % Get image file info
            imgData = imread(imgName); % Read image from file
            imgTexture = Screen('MakeTexture',window, imgData); % Convert image file to texture
            Screen('DrawTexture', window, imgTexture); % Prepare image texture on backbuffer           
            % Save complete backbuffer as trial*.bmp to be used as stimulus & DataViewer backdrop
            stimName = ['blk' num2str(iBlock) 'trial' num2str(i) '.bmp']; % Prepare stimulus file name
            stimArray = Screen('GetImage', window, [], 'drawBuffer'); % Copy backbuffer to be used as stimulus
            imwrite(stimArray, stimName); % Save .bmp stimulus file in experment folder           
            % Convert stimulus to texture so it is ready for drawing later on
            stimTexture = Screen('MakeTexture', window, stimArray); % Convert to texture
            
            % STEP 7.2: START TRIAL AND DRAW CROSSHAIRS FOR DATAVIEWER
            
            if i ~= 1 % If not first trial in a block
                % Write TRIALID message to EDF file: marks the start of a trial for DataViewer
                % TRIALID message for first trial in the block was written before recording began
                % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
                trialCount = trialCount + 1; % Add 1 to trial counter
                Eyelink('Message', 'TRIALID %d', trialCount);
            end
            
            % Write messages to EDF: prepare backdrop and draw central crosshairs for DataViewer
            % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
            Eyelink('Message', '!V CLEAR %d %d %d', round(el.backgroundcolour(1)/colorMaxVal*255), round(el.backgroundcolour(2)/colorMaxVal*255), round(el.backgroundcolour(3)/colorMaxVal*255));
            Eyelink('Message', '!V DRAWLINE 0 0 0 %d %d %d %d', round(width/2-20), round(height/2), round(width/2+20), round(height/2));
            Eyelink('Message', '!V DRAWLINE 0 0 0 %d %d %d %d', round(width/2), round(height/2-20), round(width/2), round(height/2+20));
            
            % STEP 7.3: WAIT FOR MRI TRIGGER AND PRESENT STIMULUS; CREATE DATAVIEWER STIMULUS BACKDROP AND INTEREST AREA
            
            % Present initial trial image
            Screen('DrawTexture', window, stimTexture); % Prepare stimulus texture on backbuffer            
            % Wait for MRI trigger (keyboard key "t" in this demo)
            if i == 1 % If first trial in a block
                KbReleaseWait; % Wait until user releases keys on keyboard:
                keyTrigger = KbName('t');% Identify key code for keyboard key "t"
                while 1
                    [~, ~, keyCode] = KbCheck;
                    if keyCode(keyTrigger)
                        % Write message to EDF file to mark the time when the trigger is received
                        Eyelink('Message', 'TRIGGER_RECEIVED');
                        blockOnset = GetSecs; % Block onset time
                        vbl = Screen('Flip', window); % Present stimulus
                        break;
                    end
                end
            else % All subsequent trials in block
                vbl = Screen('Flip', window, blockOnset + (trialDur*(i-1)) - 0.5*ifi); % Present stimulus. Allow half flip interval for precise flip timing)
            end
            % Write message to EDF file to mark the start time of stimulus presentation.
            Eyelink('Message', 'STIM_ONSET');
            % Write !V IMGLOAD message to EDF file: creates backdrop image for DataViewer
            % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Image Commands
            Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', stimName, width/2, height/2);
            % Write !V IAREA message to EDF file: creates interest area around image in DataViewer
            % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Interest Area Commands
            Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, round(width/2-imgInfo.Width/2), round(height/2-imgInfo.Height/2), round(width/2+imgInfo.Width/2), round(height/2+imgInfo.Height/2),'IMAGE_IA');
            
            % STEP 7.4: WAIT FOR STIMULUS TIMER; CLEAR SCREEN AND SHOW CROSSHAIRS
            
            % Prepare crosshairs for end of trial
            Screen('DrawTexture', window, backgroundTexture); % Prepare background texture on backbuffer
            Screen('DrawLine', window, 0, round(width/2-20), round(height/2), round(width/2+20), round(height/2), 5);
            Screen('DrawLine', window, 0, round(width/2), round(height/2-20), round(width/2), round(height/2+20), 5);
                % Check that eye tracker is  still recording. Otherwise close and transfer copy of EDF file to Display PC
                err = Eyelink('CheckRecording');
                if(err ~= 0)
                    fprintf('EyeLink Recording stopped!\n');
                    % Transfer a copy of the EDF file to Display PC
                    Eyelink('SetOfflineMode');% Put tracker in idle/offline mode
                    Eyelink('CloseFile'); % Close EDF file on Host PC
                    Eyelink('Command', 'clear_screen 0'); % Clear trial image on Host PC at the end of the experiment
                    WaitSecs(0.1); % Allow some time for screen drawing
                    % Transfer a copy of the EDF file to Display PC
                    transferFile; % See transferFile function below
                    cleanup; % Abort experiment (see cleanup function below)
                    return
                end
                Screen('Flip', window, stimDur + vbl - 0.5*ifi); % Present crosshairs. Allow half flip interval for precise flip timing)
                % Write message to EDF file to mark time when blank screen is presented
                Eyelink('Message', 'BLANK_SCREEN');
                
            % STEP 7.5: CREATE VARIABLES FOR DATAVIEWER; WAIT FOR END OF LAST TRIAL IN BLOCK 
            
            if i == length(imgList{iBlock}) % If last trial in block
                % Clear screen and draw crosshairs for DataViewer
                % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
                Eyelink('Message', '!V CLEAR %d %d %d', round(el.backgroundcolour(1)/colorMaxVal*255), round(el.backgroundcolour(2)/colorMaxVal*255), round(el.backgroundcolour(3)/colorMaxVal*255));
                Eyelink('Message', '!V DRAWLINE 0 0 0 %d %d %d %d', round(width/2-20), round(height/2), round(width/2+20), round(height/2));
                Eyelink('Message', '!V DRAWLINE 0 0 0 %d %d %d %d', round(width/2), round(height/2-20), round(width/2), round(height/2+20));
                while 1
                    % Check that eye tracker is  still recording. Otherwise close and transfer copy of EDF file to Display PC
                    err = Eyelink('CheckRecording');
                    if(err ~= 0)
                        fprintf('EyeLink Recording stopped!\n');
                        % Transfer a copy of the EDF file to Display PC
                        Eyelink('SetOfflineMode');% Put tracker in idle/offline mode
                        Eyelink('CloseFile'); % Close EDF file on Host PC
                        Eyelink('Command', 'clear_screen 0'); % Clear trial image on Host PC at the end of the experiment
                        WaitSecs(0.1); % Allow some time for screen drawing
                        % Transfer a copy of the EDF file to Display PC
                        transferFile; % See transferFile function below
                        cleanup; % Abort experiment (see cleanup function below)
                        return
                    end
                    % Wait for end of last trial in block (duration relative to MRI trigger time)
                    if GetSecs-blockOnset >= (trialDur*(i))
                        % Write message to EDF file to mark time when block has ended
                        Eyelink('Message', 'BLOCK_END');
                        break;
                    end
                end
            end
            
            % Write !V TRIAL_VAR messages to EDF file: creates trial variables in DataViewer
            % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Trial Message Commands
            Eyelink('Message', '!V TRIAL_VAR block %d', iBlock); % Block iteration
            Eyelink('Message', '!V TRIAL_VAR iteration %d', i); % Trial iteration
            Eyelink('Message', '!V TRIAL_VAR image %s', imgName); % Image name            
            % Write TRIAL_RESULT message to EDF file: marks the end of a trial for DataViewer
            % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
            Eyelink('Message', 'TRIAL_RESULT 0');
            WaitSecs(0.01); % Allow some time before ending the trial
            
            % Clear Screen() textures that were initialized for each trial iteration
            Screen('Close', backgroundTexture);
            Screen('Close', imgTexture);
            Screen('Close', stimTexture);
        end % End trial loop
        
        % STEP 7.6: STOP TRACKER RECORDING AT THE END OF EACH BLOCK
        
        WaitSecs(0.1); % Add 100 msec of data to catch final events before stopping
        Eyelink('StopRecording'); % Stop tracker recording
    end % End block looop
    
    
    %% STEP 8: CLOSE EDF FILE. TRANSFER EDF COPY TO DISPLAY PC. CLOSE EYELINK CONNECTION. FINISH UP
    
    % Put tracker in idle/offline mode before closing file. Eyelink('SetOfflineMode') is recommended.
    % However if Eyelink('Command', 'set_idle_mode') is used, allow 50ms before closing the file as shown in the commented code:
    % Eyelink('Command', 'set_idle_mode');% Put tracker in idle/offline mode
    % WaitSecs(0.05); % Allow some time for transition
    Eyelink('SetOfflineMode'); % Put tracker in idle/offline mode
    Eyelink('Command', 'clear_screen 0'); % Clear Host PC backdrop graphics at the end of the experiment
    WaitSecs(0.5); % Allow some time before closing and transferring file    
    Eyelink('CloseFile'); % Close EDF file on Host PC       
    % Transfer a copy of the EDF file to Display PC
    transferFile; % See transferFile function below
catch % If syntax error is detected
    cleanup;
    % Print error message and line number in Matlab's Command Window
    psychrethrow(psychlasterror);
end

% Cleanup function used throughout the script above
    function cleanup
        try
            Screen('CloseAll'); % Close window if it is open
        end
        Eyelink('Shutdown'); % Close EyeLink connection
        ListenChar(0); % Restore keyboard output to Matlab
        ShowCursor; % Restore mouse cursor
        if ~IsOctave; commandwindow; end % Bring Command Window to front
    end

% Function for transferring copy of EDF file to the experiment folder on Display PC.
% Allows for optional destination path which is different from experiment folder
    function transferFile
        try
            if dummymode ==0 % If connected to EyeLink
                % Show 'Receiving data file...' text until file transfer is complete
                Screen('FillRect', window, el.backgroundcolour); % Prepare background on backbuffer
                Screen('DrawText', window, 'Receiving data file...', 5, height-35, 0); % Prepare text
                Screen('Flip', window); % Present text
                fprintf('Receiving data file ''%s.edf''\n', edfFile); % Print some text in Matlab's Command Window
                
                % Transfer EDF file to Host PC
                % [status =] Eyelink('ReceiveFile',['src'], ['dest'], ['dest_is_path'])
                status = Eyelink('ReceiveFile');
                
                % Check if EDF file has been transferred successfully and print file size in Matlab's Command Window
                if status > 0
                    fprintf('EDF file size: %.1f KB\n', status/1024); % Divide file size by 1024 to convert bytes to KB
                end
                % Print transferred EDF file path in Matlab's Command Window
                fprintf('Data file ''%s.edf'' can be found in ''%s''\n', edfFile, pwd);
            else
                fprintf('No EDF file saved in Dummy mode\n');
            end
            cleanup;
        catch % Catch a file-transfer error and print some text in Matlab's Command Window
            fprintf('Problem receiving data file ''%s''\n', edfFile);
            cleanup;
            psychrethrow(psychlasterror);
        end
    end
end