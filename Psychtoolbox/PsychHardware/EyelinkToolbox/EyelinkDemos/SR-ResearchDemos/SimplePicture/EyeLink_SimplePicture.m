function EyeLink_SimplePicture(screenNumber)
% A simple EyeLink integration demo that records eye movements passively 
% while an image is presented on the screen. Each trial ends when the
% space bar or a button is pressed.
%
% Usage:
% Eyelink_SimplePicture(screenNumber)
%
% screenNumber is an optional parameter which can be used to pass a specific value to Screen('OpenWindow', ...)
% If screenNumber is not specified, or if isempty(screenNumber) then the default:
% screenNumber = max(Screen('Screens'));
% will be used.

% Bring the Command Window to the front if it is already open
if ~IsOctave; commandwindow; end

% Initialize PsychSound for calibration/validation audio feedback
InitializePsychSound();

% Use default screenNumber if none specified
if (nargin < 1)
    screenNumber = [];
end

try
    %% STEP 1: INITIALIZE EYELINK CONNECTION; OPEN EDF FILE; GET EYELINK TRACKER VERSION
    
    % Initialize EyeLink connection (dummymode = 0) or run in "Dummy Mode" without an EyeLink connection (dummymode = 1);
    dummymode = 0;
    
    % Optional: Set IP address of eyelink tracker computer to connect to.
    % Call this before initializing an EyeLink connection if you want to use a non-default IP address for the Host PC.
    %Eyelink('SetAddress', '10.10.10.240');
    
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
        error('Session cancelled by user'); % Abort experiment (see cleanup function below)
    end    
    edfFile = answer{1}; % Save file name to a variable    
    % Print some text in Matlab's Command Window if file name is longer than 8 characters
    if length(edfFile) > 8
        fprintf('Filename needs to be no more than 8 characters long (letters, numbers and underscores only)\n');
        error('Filename needs to be no more than 8 characters long (letters, numbers and underscores only)');
    end
 
    % Open an EDF file and name it
    failOpen = Eyelink('OpenFile', edfFile);
    if failOpen ~= 0 % Abort if it fails to open
        fprintf('Cannot create EDF file %s', edfFile); % Print some text in Matlab's Command Window
        error('Cannot create EDF file %s', edfFile); % Print some text in Matlab's Command Window
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
    window = Screen('OpenWindow', screenNumber, [128 128 128]); % Open graphics window
    Screen('Flip', window);
    % Return width and height of the graphics window/screen in pixels
    [width, height] = Screen('WindowSize', window);
    
    
    %% STEP 4: SET CALIBRATION SCREEN COLOURS/SOUNDS; PROVIDE WINDOW SIZE TO EYELINK HOST & DATAVIEWER; SET CALIBRATION PARAMETERS; CALIBRATE
    
    % Provide EyeLink with some defaults, which are returned in the structure "el".
    el = EyelinkInitDefaults(window);
    % set calibration/validation/drift-check(or drift-correct) size as well as background and target colors. 
    % It is important that this background colour is similar to that of the stimuli to prevent large luminance-based 
    % pupil size changes (which can cause a drift in the eye movement data)
    el.calibrationtargetsize = 3;% Outer target size as percentage of the screen
    el.calibrationtargetwidth = 0.7;% Inner target size as percentage of the screen
    el.backgroundcolour = [128 128 128];% RGB grey
    el.calibrationtargetcolour = [0 0 0];% RGB black
    % set "Camera Setup" instructions text colour so it is different from background colour
    el.msgfontcolour = [0 0 0];% RGB black
        
    % Use an image file instead of the default calibration bull's eye targets. 
    % Commenting out the following two lines will use default targets:
    el.calTargetType = 'image';
    el.calImageTargetFilename = [pwd '/' 'fixTarget.jpg'];
    
    % Set calibration beeps (0 = sound off, 1 = sound on)
    el.targetbeep = 1;  % sound a beep when a target is presented
    el.feedbackbeep = 1;  % sound a beep after calibration or drift check/correction
    
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
    
    spaceBar = KbName('space');% Identify keyboard key code for spacebar to end each trial later on
    imgList = {'img1.jpg' 'img2.jpg'};% Provide image list for 2 trials
    for i = 1:length(imgList)
        
        % STEP 5.1: START TRIAL; SHOW TRIAL INFO ON HOST PC; SHOW BACKDROP IMAGE AND/OR DRAW FEEDBACK GRAPHICS ON HOST PC; DRIFT-CHECK/CORRECTION
        
        % Write TRIALID message to EDF file: marks the start of a trial for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
        Eyelink('Message', 'TRIALID %d', i);
        % Write !V CLEAR message to EDF file: creates blank backdrop for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
        Eyelink('Message', '!V CLEAR %d %d %d', el.backgroundcolour(1), el.backgroundcolour(2), el.backgroundcolour(3));
        % Supply the trial number as a line of text on Host PC screen
        Eyelink('Command', 'record_status_message "TRIAL %d/%d"', i, length(imgList));
        
        % Draw graphics on the EyeLink Host PC display. See COMMANDS.INI in the Host PC's exe folder for a list of commands
        Eyelink('SetOfflineMode');% Put tracker in idle/offline mode before drawing Host PC graphics and before recording        
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
        imgName = char(imgList(i)); % Get image file name for current trial
        imgInfo = imfinfo(imgName); % Get image file info
        imgData = imread(imgName); % Get image file data
        transferStatus = Eyelink('ImageTransfer', imgData, 0, 0, 0, 0, round(width/2-imgInfo.Width/2), round(height/2-imgInfo.Height/2));

        if dummymode == 0 && transferStatus ~= 0 % If connected to EyeLink and image transfer fails
            fprintf('Image transfer Failed\n'); % Print some text in Matlab's Command Window
        end
        
        % Optional: draw feedback box and lines on Host PC interface instead of (or on top of) backdrop image.
        % See section 25.7 'Drawing Commands' in the EyeLink Programmers Guide manual
        Eyelink('Command', 'draw_box %d %d %d %d 15', round(width/2-imgInfo.Width/2), round(height/2-imgInfo.Height/2), round(width/2+imgInfo.Width/2), round(height/2+imgInfo.Height/2));
        
        % Perform a drift check/correction.
        % Optionally provide x y target location, otherwise target is presented on screen centre
        EyelinkDoDriftCorrection(el, round(width/2), round(height/2));
        
        %STEP 5.2: START RECORDING
        
        % Put tracker in idle/offline mode before recording. Eyelink('SetOfflineMode') is recommended 
        % however if Eyelink('Command', 'set_idle_mode') is used allow 50ms before recording as shown in the commented code:        
        % Eyelink('Command', 'set_idle_mode');% Put tracker in idle/offline mode before recording
        % WaitSecs(0.05); % Allow some time for transition        
        Eyelink('SetOfflineMode');% Put tracker in idle/offline mode before recording
        Eyelink('StartRecording'); % Start tracker recording
        WaitSecs(0.1); % Allow some time to record a few samples before presenting first stimulus
        
        % STEP 5.3: PRESENT STIMULUS; CREATE DATAVIEWER BACKDROP AND INTEREST AREA
        
        % Prepare and present stimulus                       
        Screen('FillRect', window, el.backgroundcolour);% Prepare grey background on backbuffer
        imgTexture = Screen('MakeTexture',window, imgData); % Convert image file to texture
        Screen('DrawTexture', window, imgTexture); % Prepare image texture on backbuffer        
        Screen('TextSize', window, 30); % Specify text size
        Screen('DrawText', window, 'Press space or button to end trial', 5, height-35, 0); % Prepare text on backbuffer  
        [~, RtStart] = Screen('Flip', window); % Present stimulus
        % Write message to EDF file to mark the start time of stimulus presentation.
        Eyelink('Message', 'STIM_ONSET');        
        % Write !V IMGLOAD message to EDF file: provides instructions for DataViewer so it will show trial stimulus as backdrop
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Image Commands
        Eyelink('Message', '!V IMGLOAD CENTER %s %d %d', imgName, width/2, height/2);        
        % Write !V IAREA message to EDF file: creates interest area around image in DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Interest Area Commands
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, round(width/2-imgInfo.Width/2), round(height/2-imgInfo.Height/2), round(width/2+imgInfo.Width/2), round(height/2+imgInfo.Height/2),'IMAGE_IA');
        
        % STEP 5.4: WAIT FOR KEYPRESS/BUTTON; SHOW BLANK SCREEN; STOP RECORDING
        
        while 1 % loop until error, space bar or button press
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
                transferFile; % See transferFile function below)
                error('EyeLink is not in record mode when it should be. Unknown error. EDF transferred from Host PC to Display PC, please check its integrity.');
            end
            % End trial if spacebar is pressed
            [~, RtEnd, keyCode] = KbCheck;
            if keyCode(spaceBar)
                % Write message to EDF file to mark the spacebar press time
                Eyelink('Message', 'KEY_PRESSED');
                reactionTime = round((RtEnd-RtStart)*1000); % Calculate RT from stimulus onset
                break; % Exit while loop
            end
            % End trial if button 5 on a supported Host PC button box is pressed
            % Use (button number * -1) + 1 to determine bitshift value
            % (e.g., button 5 should use bitshift value of -4)
            buttonResult = Eyelink('ButtonStates');
            if buttonResult
                if bitshift(buttonResult, -4) == 1
                    % Write message to EDF file to mark the button press time
                    Eyelink('Message', 'BUTTON_PRESSED');
                    reactionTime = round((GetSecs-RtStart)*1000); % Calculate RT from stimulus onset
                    break; % Exit while loop
                end
            end
        end % End of while loop
        
        % Draw blank screen at end of trial
        Screen('FillRect', window, el.backgroundcolour); % Prepare grey background on backbuffer
        Screen('Flip', window); % Present blank screen
        % Write message to EDF file to mark time when blank screen is presented
        Eyelink('Message', 'BLANK_SCREEN');
        % Write !V CLEAR message to EDF file: creates blank backdrop for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
        Eyelink('Message', '!V CLEAR %d %d %d', el.backgroundcolour(1), el.backgroundcolour(2), el.backgroundcolour(3));
        
        % Stop recording eye movements at the end of each trial
        WaitSecs(0.1); % Add 100 msec of data to catch final events before stopping
        Eyelink('StopRecording'); % Stop tracker recording
        
        % STEP 5.5: CREATE VARIABLES FOR DATAVIEWER; END TRIAL
        
        % Write !V TRIAL_VAR messages to EDF file: creates trial variables in DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Trial Message Commands
        Eyelink('Message', '!V TRIAL_VAR iteration %d', i); % Trial iteration
        Eyelink('Message', '!V TRIAL_VAR image %s', imgName); % Image name
        WaitSecs(0.001); % Allow some time between messages. Some messages can be lost if too many are written at the same time
        Eyelink('Message', '!V TRIAL_VAR rt %d', reactionTime); % Reaction time
        % Write TRIAL_RESULT message to EDF file: marks the end of a trial for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
        Eyelink('Message', 'TRIAL_RESULT 0');
        WaitSecs(0.01); % Allow some time before ending the trial
        
        % Clear Screen() textures that were initialized for each trial iteration
        Screen('Close', imgTexture);
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
    % Print error message and line number in Matlab's Command Window
    psychrethrow(psychlasterror);
end
cleanup;

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
                % Optionally uncomment below to change edf file name when a copy is transferred to the Display PC
                % % If <src> is omitted, tracker will send last opened data file.
                % % If <dest> is omitted, creates local file with source file name.
                % % Else, creates file using <dest> as name.  If <dest_is_path> is supplied and non-zero
                % % uses source file name but adds <dest> as directory path.
                % newName = ['Test_',char(datetime('now','TimeZone','local','Format','y_M_d_HH_mm')),'.edf'];                
                % status = Eyelink('ReceiveFile', [], newName, 0);
                
                % Check if EDF file has been transferred successfully and print file size in Matlab's Command Window
                if status > 0
                    fprintf('EDF file size: %.1f KB\n', status/1024); % Divide file size by 1024 to convert bytes to KB
                end
                % Print transferred EDF file path in Matlab's Command Window
                fprintf('Data file ''%s.edf'' can be found in ''%s''\n', edfFile, pwd);
            else
                fprintf('No EDF file saved in Dummy mode\n');
            end
        catch % Catch a file-transfer error and print some text in Matlab's Command Window
            fprintf('Problem receiving data file ''%s''\n', edfFile);
            psychrethrow(psychlasterror);
        end
    end
end