function OldNewRecogExp(subNo,hand)
% OldNewRecogExp(subNo,hand);
%
% Example of an old/new recognition experiment.
%
% This is an example of a simple old/new recognition experiment. It is split
% into two phases, a study phase and a test phase:
%
% Study phase:
%
% In the study phase, the subject is presented with 3 images of objects in
% randomized presentation order and has to learn/memorize them. Each image
% is presented for a 'duration' of 2 seconds, then the image disappears and
% the subject has to advance to the next presentation by pressing the 'n'
% key on the keyboard. The list of 'study' objects is read from the file
% 'studylist.txt'.
%
% Test phase:
%
% In the test phase, the subject is presented with test images defined in
% the file 'testlist.txt', again in randomized order. Test images are
% presented for a 'duration' of 0.5 seconds or until the subjects responds
% with a keypress. The subject has to press one of two keys, telling if the
% test image is an "old" image - previously presented in the study phase,
% or a "new" image - not seen in the study phase. The keys used are 'c' and
% 'm', the mapping of keys to response ("old" or "new") is selected by the
% input argument "hand" -- allowing to balance for handedness of subjects /
% response bias.
%
% At the end of a session (when all images in the 'testlist.txt' have been
% presented and tested), the results - both response of the subject and its
% reaction time - are stored to a file 'OldNewRecogExp_xx.dat', where 'xx'
% is the subject number given as input argument "subNo" to this script.
%
% Input parameters:
%
% subNo    subject number; use subNo>99 to skip check for existing file
% hand     response mapping for test phase
%
% e.g.: OldNewRecogExp(99,1);
%
% Example DESIGN:
%
% STUDY PHASE: study 3 objects
% TEST  PHASE: shown 6 objects, decide whether the object is old or new
%
% This script demonstrates:
%
%    - reading from files to get condition on each trial
%    - randomizing conditions (for the study and test phase)
%    - showing image/collecting response (response time, accuracy)
%    - writing data to file "OldNewRecogExp_<subNo>.dat
%
% Please refer to other included demos for new functions of Psychtoolbox-3
% vs. the old Psychtoolbox-2.
%
% Other than a few changes how to call the Screen function,
% and the try ... catch statements, this code replicates a
% typical OS9 experiment
%
% NOTE to previous MacOS-9 users: OSX is case sensitive!!!
%
% History:
% 
% 05/24/05 Quoc Vuong, PhD, University of NewCastle wrote and contributed
% it, as an example for usage of Psychtoolbox-3, Version 1.0.6.
% 03/01/08 Mario Kleiner modified the code to make use of new functionality
% added in Psychtoolbox 3.0.8.

%%%%%%%%%%%%%%%%%%%%%%%%%
% any preliminary stuff
%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Matlab/Octave window:
clc;

% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

% Check if all needed parameters given:
if nargin < 2
    error('Must provide required input parameters "subNo" and "hand"!');
end

% Reseed the random-number generator for each expt.
rand('state',sum(100*clock));

% Make sure keyboard mapping is the same on all supported operating systems
% Apple MacOS/X, MS-Windows and GNU/Linux:
KbName('UnifyKeyNames');

% Init keyboard responses (caps doesn't matter)
advancestudytrial=KbName('n');

% Use input variable "hand" to determine response mapping for this session.
if (hand==1)
    oldresp=KbName('c'); % "old" response via key 'c'
    newresp=KbName('m'); % "new" response via key 'm'
else
    oldresp=KbName('m'); % Keys are switched in this case.
    newresp=KbName('c');
end

%%%%%%%%%%%%%%%%%%%%%%
% file handling
%%%%%%%%%%%%%%%%%%%%%%

% Define filenames of input files and result file:
datafilename = strcat('OldNewRecogExp_',num2str(subNo),'.dat'); % name of data file to write to
studyfilename = 'studylist.txt';                           % study list
testfilename  = 'testlist.txt';                            % test list

% check for existing result file to prevent accidentally overwriting
% files from a previous subject/session (except for subject numbers > 99):
if subNo<99 && fopen(datafilename, 'rt')~=-1
    fclose('all');
    error('Result data file already exists! Choose a different subject number.');
else
    datafilepointer = fopen(datafilename,'wt'); % open ASCII file for writing
end

%%%%%%%%%%%%%%%%%%%%%%
% experiment
%%%%%%%%%%%%%%%%%%%%%%

% Embed core of code in try ... catch statement. If anything goes wrong
% inside the 'try' block (Matlab error), the 'catch' block is executed to
% clean up, save results, close the onscreen window etc.
try
    % Get screenNumber of stimulation display. We choose the display with
    % the maximum index, which is usually the right one, e.g., the external
    % display on a Laptop:
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    % Hide the mouse cursor:
    HideCursor;
    
    % Returns as default the mean gray value of screen:
    gray=GrayIndex(screenNumber); 

    % Open a double buffered fullscreen window on the stimulation screen
    % 'screenNumber' and choose/draw a gray background. 'w' is the handle
    % used to direct all drawing commands to that window - the "Name" of
    % the window. 'wRect' is a rectangle defining the size of the window.
    % See "help PsychRects" for help on such rectangles and useful helper
    % functions:
    [w, wRect]=Screen('OpenWindow',screenNumber, gray);
    
    % Set text size (Most Screen functions must be called after
    % opening an onscreen window, as they only take window handles 'w' as
    % input:
    Screen('TextSize', w, 32);
    
    % Do dummy calls to GetSecs, WaitSecs, KbCheck to make sure
    % they are loaded and ready when we need them - without delays
    % in the wrong moment:
    KbCheck;
    WaitSecs(0.1);
    GetSecs;

    % Set priority for script execution to realtime priority:
    priorityLevel=MaxPriority(w);
    Priority(priorityLevel);
    
    % run through study and test phase
    for phase=1:2 % 1 is study phase, 2 is test phase                

        % Setup experiment variables etc. depending on phase:
        if phase==1 % study phase
            
            % define variables for current phase
            phaselabel='study';
            duration=2.000; % Duration of study image presentation in secs.
            trialfilename=studyfilename;
            
            message = 'study phase ...\nstudy each picture ... press _n_ when it disappears ...\n... press mouse button to begin ...';
%             Screen('DrawText', w, , 10, 10, 255);
%             Screen('DrawText', w, , 10, 40, 255);
%             Screen('DrawText', w, , 10, 70, 255);
        else        % test phase
            
            % define variables
            phaselabel='test';
            duration=0.500;  %sec
            trialfilename=testfilename;
            
            % write message to subject
            str=sprintf('Press _%s_ for OLD and _%s_ for NEW\n',KbName(oldresp),KbName(newresp));
            message = ['test phase ...\n' str '... press mouse button to begin ...'];
%             Screen('DrawText', w, 'test phase ...', 10, 10, 255);
%             Screen('DrawText', w, str, 10, 40, 255);
%             Screen('DrawText', w, '... press mouse to begin ...', 10, 70,
%             255);
        end

        % Write instruction message for subject, nicely centered in the
        % middle of the display, in white color. As usual, the special
        % character '\n' introduces a line-break:
        DrawFormattedText(w, message, 'center', 'center', WhiteIndex(w));

        % Update the display to show the instruction text:
        Screen('Flip', w);
        
        % Wait for mouse click:
        GetClicks(w);
                
        % Clear screen to background color (our 'gray' as set at the
        % beginning):
        Screen('Flip', w);
        
        % Wait a second before starting trial
        WaitSecs(1.000);
        
        % read list of conditions/stimulus images -- textread() is a matlab function
        % objnumber  arbitrary number of stimulus
        % objname    stimulus filename
        % objtype    1=old stimulus, 2=new stimulus
        %            for study list, stimulus coded as "old"
        [ objnumber, objname, objtype ] = textread(trialfilename,'%d %s %d');
        
        % Randomize order of list
        ntrials=length(objnumber);         % get number of trials
        randomorder=randperm(ntrials);     % randperm() is a matlab function
        objnumber=objnumber(randomorder);  % need to randomize each list!
        objname=objname(randomorder);      %
        objtype=objtype(randomorder);      %
        
        % loop through trials
        for trial=1:ntrials
                        
            % wait a bit between trials
            WaitSecs(0.500);
            
            % initialize KbCheck and variables to make sure they're
            % properly initialized/allocted by Matlab - this to avoid time
            % delays in the critical reaction time measurement part of the
            % script:
            [KeyIsDown, endrt, KeyCode]=KbCheck;
            
            % read stimulus image into matlab matrix 'imdata':
            stimfilename=strcat('stims/',char(objname(trial))); % assume stims are in subfolder "stims"
            imdata=imread(char(stimfilename));
            
            % make texture image out of image matrix 'imdata'
            tex=Screen('MakeTexture', w, imdata);
            
            % Draw texture image to backbuffer. It will be automatically
            % centered in the middle of the display if you don't specify a
            % different destination:
            Screen('DrawTexture', w, tex);
            
            % Show stimulus on screen at next possible display refresh cycle,
            % and record stimulus onset time in 'startrt':
            [VBLTimestamp startrt]=Screen('Flip', w);
            
            % while loop to show stimulus until subjects response or until
            % "duration" seconds elapsed.
            while (GetSecs - startrt)<=duration
                % poll for a resp
                % during test phase, subjects can response
                % before stimulus terminates
                if ( phase==2 ) % if test phase
                    if ( KeyCode(oldresp)==1 || KeyCode(newresp)==1 )
                        break;
                    end
                    [KeyIsDown, endrt, KeyCode]=KbCheck;
                end

                % Wait 1 ms before checking the keyboard again to prevent
                % overload of the machine at elevated Priority():
                WaitSecs(0.001);
            end

            % Clear screen to background color after fixed 'duration'
            % or after subjects response (on test phase)
            Screen('Flip', w);
            
            % loop until valid key is pressed
            % if a response is made already, then this loop will be skipped
            if ( phase==1 ) % study phase
                while (KeyCode(advancestudytrial)==0)
                    [KeyIsDown, endrt, KeyCode]=KbCheck;
                    WaitSecs(0.001);
                end
            end
            
            if ( phase==2 ) % test phase
                while ( KeyCode(oldresp)==0 && KeyCode(newresp)==0 )
                    [KeyIsDown, endrt, KeyCode]=KbCheck;
                     WaitSecs(0.001);
                end
            end
            
            % compute response time
            rt=round(1000*(endrt-startrt));
            
            % compute accuracy
            if (phase==1 ) % study phase
                ac=1;
            else           % test phase
                ac=0;
                % code correct if old-response with old stimulus,
                % new-response with new stimulus, or study phase
                if ( (KeyCode(oldresp)==1 && objtype(trial)==1) || (KeyCode(newresp)==1 && objtype(trial)==2) )
                    ac=1;
                end
            end
            
            resp=KbName(KeyCode); % get key pressed by subject
            
            % Write trial result to file:
            fprintf(datafilepointer,'%i %i %s %i %s %i %s %i %i %i\n', ...
                subNo, ...
                hand, ...
                phaselabel, ...
                trial, ...
                resp, ...
                objnumber(trial), ...
                char(objname(trial)), ...
                objtype(trial), ...
                ac, ...
                rt);
            
        end % for trial loop
    end % phase loop
    
    % Cleanup at end of experiment - Close window, show mouse cursor, close
    % result file, switch Matlab/Octave back to priority 0 -- normal
    % priority:
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    
    % End of experiment:
    return;
catch
    % catch error: This is executed in case something goes wrong in the
    % 'try' part due to programming error etc.:
    
    % Do same cleanup as at the end of a regular session...
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
end % try ... catch %

% only 352 lines of codes ...  :)
