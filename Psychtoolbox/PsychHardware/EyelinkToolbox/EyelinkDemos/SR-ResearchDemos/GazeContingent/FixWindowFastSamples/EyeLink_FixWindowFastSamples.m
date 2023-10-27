function EyeLink_FixWindowFastSamples(screenNumber)
% EyeLink gaze-contingent demo showing how to retrieve fast gaze samples online.
% In each trial central crosshairs are shown until gaze is detected continuously within a central 
% square window for 500ms or until the space bar is pressed. An image is 
% then presented until the space bar is pressed to end the trial.
%
% Usage:
% Eyelink_FixWindowFastSamples(screenNumber)
%
% screenNumber is an optional parameter which can be used to pass a specific value to PsychImaging('OpenWindow', ...)
% If screenNumber is not specified, or if isempty(screenNumber) then the default:
% screenNumber = max(Screen('Screens'));
% will be used.
%
% The demo checks if a new sample is available online via the link. 
% This is equivalent to eyeLink_newest_float_sample() in C API. 
% See EyeLink Programmers Guide manual > Function Lists > Message and Command Sending/Receiving > Functions
% It allows access to the following sample properties:
%
%       time (sample time)
%       type (SAMPLE=200)
%       gx ([left gaze x, right gaze x])
%       gy ([left gaze y, right gaze y])
%       pa ([lef eye pupil size, right eye pupil size])
%       rx (x 'pixel per degree' value)
%       ry (y 'pixel per degree' value)
%       buttons (button state and changes)
%       hdata (contains a list of 8 fields. Only the first 4 values are important:
%             [uncalibrated target sticker x, uncalibrated target sticker y, target sticker distance in mm, target flags)

% Initialize PsychSound for calibration/validation audio feedback
InitializePsychSound();

% Some initial parameters:
fixWinSize = 100; % Width and Height of square fixation window [in pixels]
fixateTime = 500; % Duration of gaze inside fixation window required before stimulus presentation [ms]

% Bring the Command Window to the front if it is already open
if ~IsOctave; commandwindow; end

% Use default screenNumber if none specified
if (nargin < 1)
    screenNumber = [];
end

try
    %% STEP 1: INITIALIZE EYELINK CONNECTION; OPEN EDF FILE; GET EYELINK TRACKER VERSION
    
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
    
    
    %% STEP 2: SELECT AVAILABLE SAMPLE/EVENT DATA
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
    
    
    %% STEP 3: OPEN GRAPHICS WINDOW
    
    % Open experiment graphics on the specified screen
    if isempty(screenNumber)
        screenNumber = max(Screen('Screens')); % Use default screen if none specified
    end
    PsychDefaultSetup(2);
    [window, wRect] = PsychImaging('OpenWindow', screenNumber, GrayIndex(screenNumber)); % Open graphics window
    Screen('Flip', window);
    
    % Get max color value for rescaling  to RGB for Host PC & Data Viewer integration
    colorMaxVal = Screen('ColorRange', window);
    % Return width and height of the graphics window/screen in pixels
    [width, height] = Screen('WindowSize', window);
    
    
    %% STEP 4: SET CALIBRATION SCREEN COLOURS; PROVIDE WINDOW SIZE TO EYELINK HOST & DATAVIEWER; SET CALIBRATION PARAMETERS; CALIBRATE
    
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
    % Allow a supported EyeLink Host PC button box to accept calibration or drift-check/correction targets via button 5
    Eyelink('Command', 'button_function 5 "accept_target_fixation"');
    % Hide mouse cursor
    HideCursor(screenNumber);
    % Start listening for keyboard input. Suppress keypresses to Matlab windows.
    ListenChar(-1);
    Eyelink('Command', 'clear_screen 0'); % Clear Host PC display from any previus drawing
    % Put EyeLink Host PC in Camera Setup mode for participant setup/calibration
    EyelinkDoTrackerSetup(el);
    
    
    %% STEP 5: TRIAL LOOP.
    
    imgList = {'img1.jpg' 'img2.jpg'};% Provide image list for 2 trials
    spaceBar = KbName('space');% Identify keyboard key code for space bar to end each trial later on
    % Create central square fixation window
    fixationWindow = [-fixWinSize -fixWinSize fixWinSize fixWinSize];
    fixationWindow = CenterRect(fixationWindow, wRect);
            
    for i = 1:length(imgList) % Trial loop
        
        % STEP 5.1: PREBUILD STIMULUS (GREY BACKGROUND + IMAGE + TEXT)
        
        % Prepare grey background on backbuffer
        Screen('FillRect', window, el.backgroundcolour);
        % Use 'drawBuffer' to copy unprocessed backbuffer images without additional processing. Prevents image size info issues on Retina displays
        backgroundArray = Screen('GetImage', window, [], 'drawBuffer'); % Copy unprocessed backbuffer
        backgroundTexture = Screen('MakeTexture', window, backgroundArray); % Convert background to texture so it is ready for drawing later on        
        % Prepare image on backbuffer
        imgName = char(imgList(i)); % Get image file name for current trial
        imgInfo = imfinfo(imgName); % Get image file info
        imgData = imread(imgName); % Read image from file
        imgTexture = Screen('MakeTexture',window, imgData); % Convert image file to texture
        Screen('DrawTexture', window, imgTexture); % Prepare image texture on backbuffer        
        % Prepare text on backbuffer
        Screen('TextSize', window, 30); % Specify text size
        Screen('DrawText', window, 'Press space to end trial', 5, height-35, 0); % Prepare text on backbuffer        
        % Save complete backbuffer as trial*.bmp to be used as stimulus and as Host PC & DataViewer backdrop
        stimName = ['trial' num2str(i) '.bmp']; % Prepare stimulus file name
        stimArray = Screen('GetImage', window, [], 'drawBuffer'); % Copy backbuffer to be used as stimulus
        imwrite(stimArray, stimName); % Save .bmp stimulus file in experment folder        
        % Convert stimulus to texture so it is ready for drawing later on
        stimInfo = imfinfo(stimName); % Get stimulus info
        stimTexture = Screen('MakeTexture', window, stimArray); % Convert to texture
                
        % STEP 5.2: START TRIAL; SHOW TRIAL INFO ON HOST PC; SHOW BACKDROP IMAGE AND/OR DRAW FEEDBACK GRAPHICS ON HOST PC; DRIFT-CHECK/CORRECTION
        
        % Write TRIALID message to EDF file: marks the start of a trial for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
        Eyelink('Message', 'TRIALID %d', i);
       
        % Write !V CLEAR message to EDF file: creates blank backdrop for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
        Eyelink('Message', '!V CLEAR %d %d %d', round(el.backgroundcolour(1)/colorMaxVal*255), round(el.backgroundcolour(2)/colorMaxVal*255), round(el.backgroundcolour(3)/colorMaxVal*255));
        
        % Supply the trial number as a line of text on Host PC screen
        Eyelink('Command', 'record_status_message "TRIAL %d/%d"', i, length(imgList));
        
        % Draw graphics on the EyeLink Host PC display. See COMMANDS.INI in the Host PC's exe folder for a list of commands
        Eyelink('SetOfflineMode'); % Put tracker in idle/offline mode before drawing Host PC graphics and before recording        
        Eyelink('Command', 'clear_screen 0'); % Clear Host PC display from any previus drawing
        % Optional: Send an image to the Host PC to be displayed as the backdrop image over which 
        % the gaze-cursor is overlayed during trial recordings.
        % See Eyelink('ImageTransfer?') for information about supported syntax and compatible image formats.
        % Below, we use the new option to pass image data from imread() as the imageArray parameter, which
        % enables the use of many image formats.
        % [status] = Eyelink('ImageTransfer', imageArray, xs, ys, width, height, xd, yd, options);
        % xs, ys: top-left corner of the region to be transferred within the source image
        % width, height: size of region to be transferred within the source image (note, values of 0 will include the entire width/height)
        % xd, yd: location (top-left) where image region to be transferred will be presented on the Host PC
        % This image transfer function works for non-resized image presentation only. If you need to resize images and use this function please resize
        % the original image files beforehand
        transferStatus = Eyelink('ImageTransfer', stimArray, 0, 0, 0, 0, 0, 0);
        if dummymode == 0 && transferStatus ~= 0 % If connected to EyeLink and image transfer fails
            fprintf('Image transfer Failed\n'); % Print some text in Matlab's Command Window
        end      
        % Optional: draw feedback box and lines on Host PC interface instead of (or on top of) backdrop image.
        % See section 25.7 'Drawing Commands' in the EyeLink Programmers Guide manual
        Eyelink('Command', 'draw_box %d %d %d %d 15', fixationWindow(1), fixationWindow(2), fixationWindow(3), fixationWindow(4)); % Fixation window
        Eyelink('Command', 'draw_cross %d %d 15 ', width/2, height/2); % Central crosshairs
        
        % Perform a drift check/correction.
        % Optionally provide x y target location, otherwise target is presented on screen centre
        EyelinkDoDriftCorrection(el, round(width/2), round(height/2));
                
        %STEP 5.3: START RECORDING
        
        % Put tracker in idle/offline mode before recording. Eyelink('SetOfflineMode') is recommended 
        % however if Eyelink('Command', 'set_idle_mode') is used allow 50ms before recording as shown in the commented code:        
        % Eyelink('Command', 'set_idle_mode');% Put tracker in idle/offline mode before recording
        % WaitSecs(0.05); % Allow some time for transition            
        Eyelink('SetOfflineMode');% Put tracker in idle/offline mode before recording
        Eyelink('StartRecording'); % Start tracker recording
        WaitSecs(0.1); % Allow some time to record a few samples before presenting first stimulus
        
        % Check which eye is available online. Returns 0 (left), 1 (right) or 2 (binocular)
        eyeUsed = Eyelink('EyeAvailable');
        % Get events from right eye if binocular
        if eyeUsed == 2
            eyeUsed = 1;
        end
               
        % STEP 5.4: PRESENT CROSSHAIRS; WAIT FOR GAZE INSIDE WINDOW OR FOR KEYPRESS
              
        % Present central crosshairs on a grey background
        Screen('DrawTexture', window, backgroundTexture); % Prepare background texture on backbuffer
        Screen('DrawLine', window, 0, round(width/2-20), round(height/2), round(width/2+20), round(height/2), 5);
        Screen('DrawLine', window, 0, round(width/2), round(height/2-20), round(width/2), round(height/2+20), 5);
        [~, gazeWinStart] = Screen('Flip', window); % Present crosshairs. Start timer for fixation window
        % Write message to EDF file to mark the crosshairs presentation time.
        Eyelink('Message', 'CROSSHAIRS');
        % Write messages to EDF: prepare backdrop and draw central crosshairs in DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
        Eyelink('Message', '!V CLEAR %d %d %d', round(el.backgroundcolour(1)/colorMaxVal*255), round(el.backgroundcolour(2)/colorMaxVal*255), round(el.backgroundcolour(3)/colorMaxVal*255));
        Eyelink('Message', '!V DRAWLINE 0 0 0 %d %d %d %d', round(width/2-20), round(height/2), round(width/2+20), round(height/2));
        Eyelink('Message', '!V DRAWLINE 0 0 0 %d %d %d %d', round(width/2), round(height/2-20), round(width/2), round(height/2+20));
        % Write !V IAREA message to EDF file: creates fixation window interest areas in DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Interest Area Commands
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, fixationWindow(1), fixationWindow(2), fixationWindow(3), fixationWindow(4),'FIXWINDOW_IA');                     
        fixWinComplete = 'yes'; % Reset variable for gaze maintained inside fixation window successfully
        while 1 % loop until error or space bar press
            % Check tracker is  still recording, otherwise close and transfer copy of EDF file to Display PC
            err = Eyelink('CheckRecording');
            if(err ~= 0)
                fprintf('EyeLink Recording stopped!\n');
                % Transfer a copy of the EDF file to Display PC
                Eyelink('SetOfflineMode'); % Put tracker in idle/offline mode
                Eyelink('CloseFile'); % Close EDF file on Host PC
                Eyelink('Command', 'clear_screen 0'); % Clear trial image on Host PC at the end of the experiment
                WaitSecs(0.1); % Allow some time for screen drawing
                % Transfer a copy of the EDF file to Display PC
                transferFile; % See transferFile function below
                cleanup; % Abort experiment (see cleanup function below)
                return
            end
            % Check if a new sample is available online via the link.
            % This is equivalent to eyeLink_newest_float_sample() in C API. See EyeLink Programmers Guide manual > Function Lists > Message and Command Sending/Receiving > Functions
            if Eyelink('NewFloatSampleAvailable') > 0
                % Get sample data in a Matlab structure
                % This is equivalent to eyeLink_newest_float_sample() in C API. See EyeLink Programmers Guide manual > Function Lists > Message and Command Sending/Receiving > Functions
                evt = Eyelink('NewestFloatSample');
                % Save sample properties as variables. See EyeLink Programmers Guide manual > Data Structures > FSAMPLE
                x_gaze = evt.gx(eyeUsed+1); % [left eye gaze x, right eye gaze x] +1 as we're accessing an array
                y_gaze = evt.gy(eyeUsed+1); % [left eye gaze y, right eye gaze y]
                if inFixWindow(x_gaze,y_gaze) % If gaze sample is within fixation window (see inFixWindow function below)
                    if (GetSecs - gazeWinStart)*1000 >= fixateTime % If gaze duration >= minimum fixation window time (fxateTime)
                        break; % break while loop to show stimulus
                    end
                elseif ~inFixWindow(x_gaze,y_gaze) % If gaze sample is not within fixation window
                    gazeWinStart = GetSecs; % Reset fixation window timer
                end
            end
            % Wait for space bar to end crosshairs if participant is unable to maintain gaze inside window for duration 'fixateTime'
            [~, ~, keyCode] = KbCheck;
            if keyCode(spaceBar)
                % Write message to EDF file to mark the space bar press time
                Eyelink('Message', 'FIXATION_KEY_PRESSED');
                fixWinComplete = 'no'; % Update variable for gaze not maintained inside window
                break; % break while loop to show stimulus
            end           
        end % End of gaze-checking while loop
               
        % STEP 5.5: PRESENT STIMULUS; CREATE DATAVIEWER BACKDROP AND INTEREST AREA
        
        % Present initial trial image
        Screen('DrawTexture', window, stimTexture); % Prepare stimulus texture on backbuffer
        [~, RtStart] = Screen('Flip', window); % Present stimulus
        % Write message to EDF file to mark the start time of stimulus presentation.
        Eyelink('Message', 'STIM_ONSET');        
        % Write !V IMGLOAD message to EDF file: creates backdrop image for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Image Commands
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', stimName, width/2, height/2);        
        % Write !V IAREA message to EDF file: creates image interest area in DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Interest Area Commands
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 2, round(width/2-imgInfo.Width/2), round(height/2-imgInfo.Height/2), round(width/2+imgInfo.Width/2), round(height/2+imgInfo.Height/2),'IMAGE_IA');       
        
        % STEP 5.6: WAIT FOR KEYPRESS; SHOW BLANK SCREEN; STOP RECORDING
        
        KbReleaseWait; % Wait until space bar release if pressed in prevous while loop        
        while 1 % loop until error or space bar press            
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
            % End trial if space bar is pressed
            [~, RtEnd, keyCode] = KbCheck;
            if keyCode(spaceBar)
                % Write message to EDF file to mark the space bar press time
                Eyelink('Message', 'KEY_PRESSED');
                reactionTime = round((RtEnd - RtStart)*1000); % Calculate RT [ms] from stimulus onset
                break;
            end           
        end % End of while loop               
        % Draw blank screen at end of trial
        Screen('DrawTexture', window, backgroundTexture); % Prepare background texture on backbuffer
        Screen('Flip', window); % Present blank screen
        % Write message to EDF file to mark time when blank screen is presented
        Eyelink('Message', 'BLANK_SCREEN');
        % Write !V CLEAR message to EDF file: creates blank backdrop for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
        Eyelink('Message', '!V CLEAR %d %d %d', round(el.backgroundcolour(1)/colorMaxVal*255), round(el.backgroundcolour(2)/colorMaxVal*255), round(el.backgroundcolour(3)/colorMaxVal*255));
        
        % Stop recording eye movements at the end of each trial
        WaitSecs(0.1); % Add 100 msec of data to catch final events before stopping
        Eyelink('StopRecording'); % Stop tracker recording
             
        % STEP 5.7: CREATE VARIABLES FOR DATAVIEWER; END TRIAL
        
        % Write !V TRIAL_VAR messages to EDF file: creates trial variables in DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Trial Message Commands
        Eyelink('Message', '!V TRIAL_VAR iteration %d', i); % Trial iteration
        Eyelink('Message', '!V TRIAL_VAR image %s', imgName); % Image name
        WaitSecs(0.001); % Allow some time between messages. Some messages can be lost if too many are written at the same time
        Eyelink('Message', '!V TRIAL_VAR fix_completed %s', fixWinComplete); % Was gaze maintained inside fixation window successfully (yes/no)?        
        Eyelink('Message', '!V TRIAL_VAR rt %d', reactionTime); % Key press RT [ms] from stimulus onset                
        % Write TRIAL_RESULT message to EDF file: marks the end of a trial for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
        Eyelink('Message', 'TRIAL_RESULT 0');
        WaitSecs(0.01); % Allow some time before ending the trial
        
        % Clear Screen() textures that were initialized for each trial iteration
        Screen('Close', backgroundTexture);
        Screen('Close', imgTexture);
        Screen('Close', stimTexture);       
    end % End trial loop
    
    
    %% STEP 6: CLOSE EDF FILE. TRANSFER EDF COPY TO DISPLAY PC. CLOSE EYELINK CONNECTION. FINISH UP
    
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

% Function that determines if gaze x y coordinates are within fixation window
    function fix = inFixWindow(mx,my)        
        fix = mx > fixationWindow(1) &&  mx <  fixationWindow(3) && ...
            my > fixationWindow(2) && my < fixationWindow(4) ;
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
