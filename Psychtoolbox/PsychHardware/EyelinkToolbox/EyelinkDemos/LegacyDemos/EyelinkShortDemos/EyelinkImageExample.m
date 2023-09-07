function result=EyelinkImageExample

% Short MATLAB example program that uses the Eyelink and Psychophysics
% Toolboxes to measure viewing behaviour for an image.
%
% History
% 15-06-10  fwc     created it, based on EyelinkExample.m
try
    
    fprintf('EyelinkToolbox Image View Example\n\n\t');
    
    dummymode=0;       % set to 1 to initialize in dummymode

    % STEP 1
    % Open a graphics window on the main screen
    screenNumber=max(Screen('Screens'));
    [window, wRect]=Screen('OpenWindow', screenNumber);
    
    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    el=EyelinkInitDefaults(window);
   
    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end
    
    [v vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % open file for recording data
    edfFile='demo.edf';
    Eyelink('Openfile', edfFile);
 
    % STEP 4
    % Do setup and calibrate the eye tracker
    EyelinkDoTrackerSetup(el);

    % do a final check of calibration using driftcorrection
    % You have to hit esc before return.
    EyelinkDoDriftCorrection(el);
    
    % STEP 5
    % Start recording eye position
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.5);

    % STEP 6
    % Show image on display
    basepath = [ PsychtoolboxRoot 'PsychDemos' filesep ];
    myimg = [basepath 'konijntjes1024x768.jpg'];

    imdata=imread(myimg);
    imtex=Screen('MakeTexture', el.window, imdata);
    
    Screen('FillRect', el.window, el.backgroundcolour);
    Screen('DrawTexture', el.window, imtex, [], wRect);  % fill screen with image
    
    % Some useful info for user about what's happening.
    [width, height]=Screen('WindowSize', el.window);

    Screen('TextFont', el.window, el.msgfont);
    Screen('TextSize', el.window, el.msgfontsize);
    Screen('DrawText', el.window, 'Just look around for a while.', 200, height-el.msgfontsize-20, el.msgfontcolour);

    % Show result on screen:
    Screen('Flip', el.window)
      
    % mark zero-plot time in data file
    Eyelink('Message', 'SYNCTIME');
    % wait a while to record a bunch of samples  
    WaitSecs(3);
    
    % STEP 7 remove image
    Screen('FillRect', el.window, el.backgroundcolour);
    Screen('DrawText', el.window, 'Done', width/2, height/2, el.msgfontcolour);
    Screen('Flip', el.window);
    % mark image removal time in data file
    Eyelink('Message', 'ENDTIME');
    WaitSecs(0.5);
    Eyelink('Message', 'TRIAL_END');
    
    % STEP 8
    % finish up: stop recording eye-movements, 
    % close graphics window, close data file and shut down tracker
    Eyelink('StopRecording');
    Eyelink('CloseFile');    
    cleanup;
    
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if it's open.
    cleanup;
    psychrethrow(psychlasterror);
end %try..catch.


% Cleanup routine:
function cleanup
    % Shutdown Eyelink:
    Eyelink('Shutdown');

    % Close window:
    sca;
    commandwindow;
    % Restore keyboard output to Matlab:
    ListenChar(0);
