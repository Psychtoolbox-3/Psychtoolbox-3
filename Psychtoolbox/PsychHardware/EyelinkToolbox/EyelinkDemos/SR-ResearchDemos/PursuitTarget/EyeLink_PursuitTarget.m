function EyeLink_PursuitTarget(screenNumber)
% A smooth pursuit EyeLink integration demo that records eye movements
% while a target moves sinusoidally across the screen. Each trial ends after 5s.
%
% Illustrates how to:
%   - change the drift-check/correction target location before each trial
%   - create a moving target for Data Viewer's Play Back Animation view
%   - create dynamic target location for Data Viewer's Temporal Graph view and sample reports
%   - create target dynamic interest areas for Data Viewer
%
% Usage:
% Eyelink_PursuitTarget(screenNumber)
%
% screenNumber is an optional parameter which can be used to pass a specific value to PsychImaging('OpenWindow', ...)
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
    window = PsychImaging('OpenWindow', screenNumber, GrayIndex(screenNumber)); % Open graphics window
    Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
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
    
    % Specify some initial parameters:
    targetSize = 10; % radius of pursuit target [pixels]
    IAsize = 60; % radius of dynamic IA ellipse [pixels]
    targetDuration = 5000; % Trial duration [ms]
    sinePlotX = round(width/2); % X centre of target sinusoidal pattern (screen centre by default)
    sinePlotY = round(height/2); % Y centre of target sinusoidal pattern (screen centre by default)
    amplitudeX = round(width/3); % X amplitude of sinusoidal pattern [pixels]
    amplitudeY = round(height/3); % Y amplitude of sinusoidal pattern [pixels]
    % Assign a label for each trial. Used later on for the Data Viewer variable "type"
    type = {'HorizontalSlow' 'HorizontalFast' 'VerticalSlow' 'VerticalFast' 'EllipticSlow' 'EllipticFast'};
    % Assign the freq_x freq_y phase_x phase_y for each trial's sinusoidal pattern
    trials = [[0.2;0;270;0] [0.5;0;90;0] [0;0.2;0;0] [0;0.3;0;180] [0.2;0.2;270;180] [0.3;0.3;270;0]];
    
    for i = 1:length(trials)
        
        % STEP 5.1: OPEN DYNAMIC IA TEXT FILE; PREPARE TARGET SINUSOIDAL PATTERN
        
        frameNo = 0; % Reset frame counter variable at each trial        
        % Open text file for each trial. This will be used later on to write each instance of a dynamic interest area
        % See DataViewer manual section: Working with Events, Samples and Interest Areas > Interest Areas
        IAfileName = [edfFile '_' 'IA_' num2str(i) '.ias'];
        f = fopen(IAfileName, 'w'); % Open text file before saving dynamic interest area data info in it
        
        % Prepare sinusoidal pattern and target
        phaseX = (trials(3,i)/360 + ((0)) * trials(1,i)) * 2*pi; % Start phase x in radians
        phaseY = (trials(4,i)/360 + ((0)) * trials(2,i))* 2*pi; % Start phase y in radians
        x =  sinePlotX + amplitudeX * sin(phaseX); % Sine pattern x
        y =  sinePlotY + amplitudeY * sin(phaseY); % Sine pattern y
        targ([1 3]) = [x-targetSize x+targetSize]; % Target x
        targ([2 4]) =  [y-targetSize y+targetSize]; % Target y
        
        % STEP 5.2: START TRIAL; SHOW TRIAL INFO ON HOST PC; DRAW FEEDBACK GRAPHICS ON HOST PC; DRIFT-CHECK/CORRECTION
        
        % Write TRIALID message to EDF file: marks the start of a trial for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
        Eyelink('Message', 'TRIALID %d', i);
        % Write !V CLEAR message to EDF file: creates blank backdrop for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
        Eyelink('Message', '!V CLEAR %d %d %d', round(el.backgroundcolour(1)/colorMaxVal*255), round(el.backgroundcolour(2)/colorMaxVal*255), round(el.backgroundcolour(3)/colorMaxVal*255));
        % Supply the trial number and type as a line of text on Host PC screen
        Eyelink('Command', 'record_status_message "TRIAL %d/%d %s"', i, length(trials), char(type(i)));
        
        % Draw target trajectory lines on the EyeLink Host PC display.
        % See section 25.7 'Drawing Commands' in the EyeLink Programmers Guide manual
        Eyelink('SetOfflineMode');% Put tracker in idle/offline mode before drawing Host PC graphics and before recording
        Eyelink('Command', 'clear_screen 0'); % Clear Host PC display from any previus drawing
        Eyelink('Command', 'draw_line %d %d %d %d 15', (width/2)-amplitudeX, height/2, (width/2)+amplitudeX, height/2);
        Eyelink('Command', 'draw_line %d %d %d %d 15', (width/2), height/2-amplitudeY, (width/2), height/2+amplitudeY);
        WaitSecs(0.1);  % Allow some time for drawing
        
        % Perform a drift check/correction.
        % Present the drift-check/correction target at each trial's start x y pursuit target location
        EyelinkDoDriftCorrection(el, round(x), round(y));
        
        %STEP 5.3: START RECORDING
        
        % Put tracker in idle/offline mode before recording. Eyelink('SetOfflineMode') is recommended 
        % however if Eyelink('Command', 'set_idle_mode') is used allow 50ms before recording as shown in the commented code:        
        % Eyelink('Command', 'set_idle_mode');% Put tracker in idle/offline mode before recording
        % WaitSecs(0.05); % Allow some time for transition            
        Eyelink('SetOfflineMode');% Put tracker in idle/offline mode before recording
        Eyelink('StartRecording'); % Start tracker recording
        WaitSecs(0.1); % Allow some time to record a few samples before presenting first stimulus
        
        % STEP 5.4: PRESENT TARGET; CREATE PURSUIT TARGET AND DYNAMIC IA FOR DATAVIEWER; STOP RECORDING 
        
        while 1 % loop until error, or timeout
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
            % Prepare and show the screen.
            % Enable alpha blending for drawing of smooth points
            Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('FillRect', window, el.backgroundcolour);
            Screen('FillOval', window,[255 0 0], targ);
            [~, stTime] = Screen('Flip', window); % Present stimulus
            frameNo = frameNo + 1; % Update frame counter
            if frameNo == 1 % If first frame...
                % Write message to EDF file to mark the start time of stimulus presentation
                Eyelink('Message', 'TARGET_ONSET');
                % Write !V IAREA FILE message to EDF file: allows the correct dynamic IA file to be used
                % by Data Viewer. Start of dynamic IA will be synced to this message
                % See DataViewer manual section: Working with Events, Samples and Interest Areas > Interest Areas
                Eyelink('Message', '!V IAREA FILE %s', IAfileName);
                stStart = stTime; % Pursuit target start time
            else
                % Write an instance of dynamic IA in the text file for each frame
                % See DataViewer manual section: Working with Events, Samples and Interest Areas > Interest Areas
                fprintf(f,'%d %d ELLIPSE 1 %d %d %d %d TARGET\n', round((previousTime-stStart)*1000)*-1, round(((stTime-stStart)*1000)-1)*-1,...
                    round(previousX)-IAsize, round(previousY)-IAsize, round(previousX)+IAsize, round(previousY)+IAsize);
            end
            
            % Write !V TARGET_POS message to EDF file: allows target location at each frame to be available in Data Viewer's
            % Temporal Graph view and sample reports.
            % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Target Position Commands
            Eyelink('Message', '!V TARGET_POS TARG1 (%d, %d) 1 0', round(x), round(y));
            % Write !V CLEAR message to EDF file: creates blank backdrop for DataViewer before drawing target location in
            % Data Viewer's Play Back Animation view.
            % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
            Eyelink('Message', '!V CLEAR %d %d %d', round(el.backgroundcolour(1)/colorMaxVal*255), round(el.backgroundcolour(2)/colorMaxVal*255), round(el.backgroundcolour(3)/colorMaxVal*255));
            % Write !V FIXPOINT message to EDF file: draws the new target location in DataViewer's Play Back Animation View
            % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing Commands
            Eyelink('Message', '!V FIXPOINT %d %d %d %d %d %d %d %d %d %d', 255, 0, 0, 255, 0, 0, round(x), round(y), targetSize*2, targetSize*2);
            
            % Keep some target information available for next frame. Used in creating an instance of dynamic IA at every frame
            previousTime = stTime;
            previousX = x;
            previousY = y;
            
            % Prepare next frame's target location
            phaseX = (trials(3,i)/360 + ((GetSecs-stStart)) * trials(1,i))*2*pi; % Current phase x in radians
            phaseY = (trials(4,i)/360 + ((GetSecs-stStart)) * trials(2,i))*2*pi; % Current phase y in radians
            x =  sinePlotX + amplitudeX * sin(phaseX); % Sine pattern x
            y =  sinePlotY + amplitudeY * sin(phaseY); % Sine pattern y
            targ([1 3]) = [x-targetSize x+targetSize]; % Target x
            targ([2 4]) = [y-targetSize y+targetSize]; % Target y
            
            % Break loop when target duration reached
            if GetSecs-stStart >= targetDuration/1000
                break
            end
        end
        
        % Draw blank screen at end of trial
        Screen('FillRect', window, el.backgroundcolour);
        [~, blankTime] = Screen('Flip', window); % Present blank screen
        % Write message to EDF file to mark time when blank screen is presented
        Eyelink('Message', 'BLANK_SCREEN');
        % Write !V CLEAR message to EDF file: creates blank backdrop for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Simple Drawing
        Eyelink('Message', '!V CLEAR %d %d %d', round(el.backgroundcolour(1)/colorMaxVal*255), round(el.backgroundcolour(2)/colorMaxVal*255), round(el.backgroundcolour(3)/colorMaxVal*255));
        
        % Finish writing last instance of dynamic IA in the text file for each trial
        % See DataViewer manual section: Working with Events, Samples and Interest Areas > Interest Areas
        fprintf(f,'%d %d ELLIPSE 1 %d %d %d %d TARGET\n', round((stTime-stStart)*1000)*-1, round(((blankTime-stStart)*1000))*-1,...
            round(previousX)-IAsize, round(previousY)-IAsize, round(x)+IAsize, round(y)+IAsize);
        
        % Stop recording eye movements at the end of each trial
        WaitSecs(0.1); % Add 100 msec of data to catch final events before stopping
        Eyelink('StopRecording'); % Stop tracker recording
        
        % STEP 5.5: CREATE VARIABLES FOR DATAVIEWER; END TRIAL; CLOSE DYNAMIC IA FILE
        
        % Write !V TRIAL_VAR messages to EDF file: creates trial variables in DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Trial Message Commands
        Eyelink('Message', '!V TRIAL_VAR iteration %d', i); % Trial iteration
        Eyelink('Message', '!V TRIAL_VAR type %s', char(type(i))); % Image name
        Eyelink('Message', '!V TRIAL_VAR frequency_x %s', num2str(trials(1,i)));
        Eyelink('Message', '!V TRIAL_VAR frequency_y %s', num2str(trials(2,i)));
        WaitSecs(0.001); % Allow some time between messages. Some messages can be lost if too many are written at the same time        
        Eyelink('Message', '!V TRIAL_VAR amplitude_x %d', round(amplitudeX));
        Eyelink('Message', '!V TRIAL_VAR amplitude_y %d', round(amplitudeY));
        Eyelink('Message', '!V TRIAL_VAR phase_x %d', trials(3,i));
        Eyelink('Message', '!V TRIAL_VAR phase_y %d', trials(4,i));        

        
        % Write TRIAL_RESULT message to EDF file: marks the end of a trial for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
        Eyelink('Message', 'TRIAL_RESULT 0'); 
        WaitSecs(0.01); % Allow some time before ending the trial
        fclose(f); % Close dynamic IA file
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