function EyelinkExample

% Short MATLAB example program that uses the Eyelink and Psychophysics
% Toolboxes to create a real-time gaze-dependent display.
% This is the example as shown in the EyelinkToolbox article in BRMIC
% Cornelissen, Peters and Palmer 2002), but updated to use new routines
% and functionality.
%
% History
% ~2006     fwc    created it, to use updated functions
% 15-06-10  fwc    updated to enable eye image display
% 17-06-10  fwc    made colour of the gaze dot change, just for fun

clear all;
commandwindow;

PsychDefaultSetup(1);

try
    
    fprintf('EyelinkToolbox Example\n\n\t');
    dummymode=0;       % set to 1 to initialize in dummymode (rather pointless for this example though)
    
    % STEP 1
    % Open a graphics window on the main screen
    % using the PsychToolbox's Screen function.
    screenNumber=max(Screen('Screens'));
    window=Screen('OpenWindow', screenNumber);
     
    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    el=EyelinkInitDefaults(window);
    % Disable key output to Matlab window:
    ListenChar(2);
        
    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode, 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end
    
    [v vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
    
    % open file to record data to
    edfFile='demo.edf';
    Eyelink('Openfile', edfFile);
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);

    % do a final check of calibration using driftcorrection
    EyelinkDoDriftCorrection(el);
    
    % STEP 5
    % start recording eye position
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.1);
    % mark zero-plot time in data file
    Eyelink('Message', 'SYNCTIME');
    stopkey=KbName('space');
    eye_used = -1;

    Screen('FillRect', el.window, el.backgroundcolour);
    Screen('TextFont', el.window, el.msgfont);
    Screen('TextSize', el.window, el.msgfontsize);
    [width, height]=Screen('WindowSize', el.window);
    message='Press space to stop.';
    Screen('DrawText', el.window, message, 200, height-el.msgfontsize-20, el.msgfontcolour);
    Screen('Flip',  el.window, [], 1);
    
    % STEP 6
    % show gaze-dependent display
    while 1 % loop till error or space bar is pressed
        % Check recording status, stop display if error
        error=Eyelink('CheckRecording');
        if(error~=0)
            break;
        end
        % check for keyboard press
        [keyIsDown, secs, keyCode] = KbCheck;
        % if spacebar was pressed stop display
        if keyCode(stopkey)
            break;
        end
        % check for presence of a new sample update
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            evt = Eyelink( 'NewestFloatSample');
            if eye_used ~= -1 % do we know which eye to use yet?
                % if we do, get current gaze position from sample
                x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used+1);
                % do we have valid data and is the pupil visible?
                if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                    % if data is valid, draw a circle on the screen at current gaze position
                    % using PsychToolbox's Screen function
                    gazeRect=[ x-9 y-9 x+10 y+10];
                    colour=round(rand(3,1)*255); % coloured dot
                    Screen('FillOval', window, colour, gazeRect);
                    Screen('Flip',  el.window, [], 1); % don't erase
                else
                    % if data is invalid (e.g. during a blink), clear display
                    Screen('FillRect', window, el.backgroundcolour);
                    Screen('DrawText', window, message, 200, height-el.msgfontsize-20, el.msgfontcolour);
                    Screen('Flip',  el.window, [], 1); % don't erase
                end
            else % if we don't, first find eye that's being tracked
                eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                if eye_used == el.BINOCULAR; % if both eyes are tracked
                    eye_used = el.LEFT_EYE; % use left eye
                end
            end
        end % if sample available
    end % main loop
    % wait a while to record a few more samples
    WaitSecs(0.1);
    
    % STEP 7
    % finish up: stop recording eye-movements,
    % close graphics window, close data file and shut down tracker
    Eyelink('StopRecording');
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
    catch rdf
        fprintf('Problem receiving data file ''%s''\n', edfFile );
        rdf;
    end
    
    cleanup;
    
catch myerr
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    cleanup;
    commandwindow;
    myerr;
    myerr.message
    myerr.stack.line

end %try..catch.


% Cleanup routine:
function cleanup
% Shutdown Eyelink:
Eyelink('Shutdown');

% Close window:
sca;

% Restore keyboard output to Matlab:
ListenChar(0);
