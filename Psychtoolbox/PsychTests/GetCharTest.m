function GetCharTest(startAt)
% GetCharTest([startAt])
%
% Test the Psychtoolbox function "GetChar".
%
% For platforms such as OS X where GetCharTest is divided into subtests,
% begin at subtest "startAt" if it is in the allowable range.  Otherwise
% ignore the argument and begin with the first subtest.    
%
% GetCharTest tests GetChar for likely failure modes. It also tests associated
% functions for character events, such as CharAvail, ListenChar and FlushEvents. 

% HISTORY
% 5/31/02 ds   Written by Dan Shima <daniel.shima@vanderbilt.edu>
% 6/1/02  dgp  Polished.
% 6/3/02  dgp  Added 4, above.
% 3/3/06  awi  Added OS X section.
% 6/15/06 awi  Added FlushEvents call after opening onscreen window.  
%               Erases panic keypresses made during l o n g time it takes to open the onscreen window.   
%              Added WaitSecs and FlushEvents to onscreen window test to block keypresses carried over 
%               from previous test.  
% 10/14/06 dhb  Removed OS 9 branch.
% 

% NOTES
%
% 3/4/06 - awi Turned off video synch testing in screen preferences during
% GetCharTest.  For purposes of testing GetChar, we are not concerned with
% video synching. Generaly, PsychTests should be orthogonal. So do not fail the
% GetChar test because we fail the synchtest.
%
% 3/9/06 - awi In section 5 of the test, a better way to break lines would
% be to set breaks at fixed pixel widths, not fixed character widths. We
% use proportially spaced fonts, so fixed-character-width lines give ragged
% line breaks.  

% We test for some common failings of previous implementations of
% GetChar.  The labeling of each testing step is somewhat
% inaccurate because
fprintf('\nGetCharTest tests the Psychtoolbox function "GetChar".\n\n');
query=1;
% Enable character key press listening. We use mode 2, which will
% suppress any output to Matlabs windows.
try
    ListenChar(2);
    FlushEvents('keyDown');
    while query
        fprintf('Proceed with the test [Y/N]? ');
        response=GetChar;
        fprintf([upper(response) '\n']);
        yesFlag=(strcmp(upper(response), 'Y'));
        noFlag=(strcmp(upper(response), 'N'));
        if ~(yesFlag || noFlag)
            fprintf('\nI did not understand your response.  Please press the "Y" or "N" key on your keyboard.\n');
        else
            query=0;
        end
    end
    if ~yesFlag
        fprintf('You chose not to run the test. Bye!\n');
    else
        fprintf('\n');
        fprintf('Proceeding with GetChar test...\n');
        % TEST 1 of 6: CharAvail
        if nargin==0 || (nargin==1 && startAt == 1)
            fprintf('\n');
            fprintf('GetChar test 1 of 7: CharAvail\n');
            fprintf('  CharAvail reports if a keyboard key has been pressed, but without removing keys from the GetChar queue.\n');
            fprintf('  We will test CharAvail by flushing all characters from the GetChar queue, looping until CharAvail reports\n');
            fprintf('  a keypress and then reading that keypress using GetChar.\n');
            fprintf('     Calling "FlushEvents(''keyDown'')" to remove queued keypresses...');
            FlushEvents('keyDown');
            fprintf(' Done.\n');
            fprintf('     Calling "CharAvail" within a loop.  Press a character key to continue...');
            gotChar=0;
            while ~gotChar
                gotChar=CharAvail;
            end
            fprintf(' Detected a character using CharAvail.\n');
            fprintf('     Reading the character with GetChar...');
            c=GetChar;
            fprintf(['You pressed: ' c '\n']);
        end
        % TEST 2 of 6: FlushEvents
        if nargin==0 || (nargin==1 && startAt <= 2)
            fprintf('\n');
            fprintf('GetChar test 2 of 7: FlushEvents\n');
            fprintf('  FlushEvents removes all characters from the GetChar queue.\n');
            fprintf('  We will test FlushEvents by gathering keypresses, flushing them, then checking to see if the queue is empty\n');
            fprintf('  by using CharAvail.\n');
            fprintf('     Calling "FlushEvents(''keyDown'')" to empty the queue before beginning...');
            FlushEvents('keyDown');
            fprintf(' Done.\n');
            fprintf('     Testing that the queue is empty by calling CharAvail... ');
            gotChar=CharAvail;
            if gotChar
                fprintf('\n');
                error('CharAvail detected a character after calling "FlushEvents".');
            else
                fprintf(' Passed.\n');
            end
            fprintf('     Gathering a keypress to occupy the queue:  Please type a character key on your keyboard... ');
            gotChar=0;
            while ~gotChar
                gotChar=CharAvail;
            end
            fprintf(' Done.\n');
            fprintf('     Flushing the keyboard character from the queue by calling "FlushEvents(''keyDown'')"... ');
            FlushEvents('keyDown');
            fprintf(' Done.\n');
            gotChar=CharAvail;
            fprintf('     Testing that the queue is empty by calling CharAvail... ');
            if gotChar
                fprintf('\n');
                error('CharAvail detected a character after calling "FlushEvents".');
            else
                fprintf(' Passed.\n');
            end
        end
        % TEST 3 of 6: ListenChar
        if nargin==0 || (nargin==1 && startAt <= 3)
            fprintf('\n');
            fprintf('GetChar test 3 of 7: ListenChar\n');
            fprintf('  ListenChar directs keyboard input to the GetChar queue.  We will test ListenChar by first.\n');
            fprintf('  directing keyboard input to the MATLAB command window, then redirecting it to the GetChar queue\n');
            fprintf('  using ListenChar, then checking for a character.\n');
            fprintf('     Calling "FlushEvents(''keyDown'')" to empty the queue before beginning...');
            FlushEvents('keyDown');
            fprintf(' Done.\n');
            % Enable keyboard dispatch to Matlab again:
            ListenChar(1);
            fprintf('     Directing keyboard input into the MATLAB Command Window by using the "input" command.\n');
            input('     Hit the <RETURN> key to continue.', 's');
            fprintf('     Directing keyboard input into the GetChar queue by calling "ListenChar"... ');
            ListenChar(2);
            fprintf(' Done.\n');
            fprintf('     Testing that GetChar receives keys.  Press any character key on your keyboard... ');
            FlushEvents('keyDown');
            c=GetChar;
            fprintf(['You pressed: ' c '\n']);
        end
        % TEST 4 of 6: Typing
        if nargin==0 || (nargin==1 && startAt <= 4)
            fprintf('\n');
            fprintf('GetChar test 4 of 7: Typing Test\n');
            fprintf('  Getchar should reliably return characters.  Type and display characters to check for errors.\n');
            fprintf('      Entering the GetChar loop, Use the "Q" key to quit...\n');
            fprintf('\n      ');
            q=0;
            lineAcc=0;
            while ~q
                c=GetChar;
                fprintf(c);
                lineAcc=lineAcc+1;
                if lineAcc == 80
                    fprintf('\n      ');
                    lineAcc=0;
                end
                q=streq('Q',upper(c));
            end
            fprintf('\n');
        end
        % TEST 5 of 6: Compatibility with Fullscreen Windows
        closeSecondWindow=0;
        if nargin==0 || (nargin==1 && startAt <= 5)
            fprintf('\n');
            fprintf('GetChar test 5 of 7: Interoperability with Fullscreen Windows\n');
            fprintf('  Test Using GetChar at the same time fullscreen windows are open.  We will open fullscreen.\n');
            fprintf('  Windows on all screens and then use GetChar to gather and display kepresses.\n');
            % Pause...
            WaitSecs(2);
            % ...Then block any inadvertant carry over from the previous test.  Too much fast typying on that and we can breeze
            % right through this test.
            FlushEvents('keyDown');
            fprintf('      Press any character key to continue.\n');
            GetChar;
            try
                skipTestFlagOld=Screen('Preference','SkipSyncTests');
                suppressAllWarningsFlagOld=Screen('Preference','SuppressAllWarnings');
                Screen('Preference','SkipSyncTests', 2);
                % Warning: Only use SuppressAllWarnings in tests independent of those warnings.
                Screen('Preference','SuppressAllWarnings', 1);
                screenNumbers= Screen('Screens');
                for i=1:length(screenNumbers)
                    screenWindow(i)= Screen('OpenWindow', screenNumbers(i));
                    Screen('FillRect', screenWindow(i), WhiteIndex(screenNumbers(i)));
                    Screen('Flip', screenWindow(i));
                    Screen('FillRect', screenWindow(i), WhiteIndex(screenNumbers(i)));
                end
                if length(screenNumbers) == 2 && closeSecondWindow
                    Screen('Close', screenWindow(2));
                end
                % Draw instructions to the screen.  We will do all our work
                % on Screen 0.  We have opened other onscreen windows only to
                % demonstrate that GetChar works then.
                if ~IsLinux
                    Screen('TextFont', screenWindow(1), 'Times');
                    Screen('TextSize', screenWindow(1), 40);
                end

                Screen('TextColor', screenWindow(1), [0 0 0]);
                Screen('DrawText', screenWindow(1), 'Hello! Welcome to the Fullscreen Window section of GetCharTest.', 100, 100);
                Screen('DrawText', screenWindow(1), 'Press any character key to continue.', 100, 150);
                Screen('Flip', screenWindow(1));
                x=GetChar;
                Screen('FillRect', screenWindow(1), WhiteIndex(screenNumbers(1)));
                Screen('DrawText', screenWindow(1), 'This is the typing test.', 100, 100);
                Screen('DrawText', screenWindow(1), 'Please type on your keyboard. Use the "Q" key to quit.', 100, 150);
                Screen('Flip', screenWindow(1));
                FlushEvents('keyDown');
                charBuffer='';
                stopLoop=0;
                lineBuffer={};
                totalChars=0;
                lineBufferIndex=0;
                startingYPosition=200;
                textHeight=0;
                textYPosition=startingYPosition;
                Screen('TextSize', screenWindow(1), 20);
                charsPerLine=100;
                totalLines=5;
                while ~stopLoop
                    % Get the character
                    c=GetChar;
                    totalChars=totalChars+1;
                    % test for doneness
                    if streq(upper(c), 'Q') || lineBufferIndex > totalLines
                        break;
                    end
                    % break lines
                    charBuffer=[charBuffer c];
                    if length(charBuffer) == charsPerLine
                        lineBufferIndex= lineBufferIndex+1;
                        lineBuffer{lineBufferIndex}= charBuffer;
                        charBuffer='';
                    end
                    % Draw the characters...
                    Screen('FillRect', screenWindow(1), WhiteIndex(screenNumbers(1)));
                    % First draw the instructions
                    %Screen('TextSize', screenWindow(1), 40);
                    Screen('DrawText', screenWindow(1), 'This is the typing test.', 100, 100);
                    Screen('DrawText', screenWindow(1), 'Please type on your keyboard. Use the "Q" key to quit.', 100, 150);
                    Screen('TextSize', screenWindow(1), 20);
                    % then display full lines...
                    textYPosition=startingYPosition;
                    for i=1:lineBufferIndex
                        Screen('DrawText', screenWindow(1), lineBuffer{i}, 100, textYPosition);
                        if  textHeight==0
                            textRect= Screen('TextBounds', screenWindow(1), lineBuffer{i});
                            textHeight= RectHeight(textRect);
                        end
                        textYPosition= floor(textYPosition + textHeight + 0.2 * textHeight);
                    end
                    % then display remainder lines.
                    if length(charBuffer) > 0;
                        %textYPosition= floor(startingYPosition + (i+1) * (textHeight + 0.2 * textHeight));
                        Screen('DrawText', screenWindow(1), charBuffer, 100, textYPosition);
                    end
                    Screen('Flip',  screenWindow(1));
                end
                % Draw the characters...
                Screen('FillRect', screenWindow(1), WhiteIndex(screenNumbers(1)));
                % First draw the instructions
                if ~IsLinux
                    Screen('TextSize', screenWindow(1), 40);
                end
                Screen('DrawText', screenWindow(1), 'End of Typing Test.', 100, 100);
                Screen('TextSize', screenWindow(1), 20);
                Screen('DrawText', screenWindow(1), 'Either you ended the test by typing Q, or you reached the limit of the test.', 100, startingYPosition);
                Screen('Flip', screenWindow(1));
                WaitSecs(4); % In rapid typing the preceeding message vanishes without a pause here.
                FlushEvents('keyDown');
                Screen('DrawText', screenWindow(1), 'Press any character key to continue.', 100, startingYPosition);
                Screen('Flip', screenWindow(1));
                GetChar;
                Screen('CloseAll');
            catch
                Screen('CloseAll');
                Screen('Preference','SkipSyncTests', skipTestFlagOld);
                Screen('Preference','SuppressAllWarnings', suppressAllWarningsFlagOld);
                error(lasterr);
            end
            Screen('Preference','SkipSyncTests', skipTestFlagOld);
            Screen('Preference','SuppressAllWarnings', suppressAllWarningsFlagOld);
            FlushEvents('keyDown');
        end
        % TEST 6 of 6: GetChar during WaitSecs
        if nargin==0 || (nargin==1 && startAt <= 6)
            fprintf('\n');
            fprintf('GetChar test 6 of 7: GetChar during WaitSecs\n');
            fprintf('  Test Using GetChar to gather keypresses during WaitSecs. Unlike KbCheck, GetChar should.\n');
            fprintf('  detect keypresses in the background during WaitSecs.\n');
            fprintf('\n');
            fprintf('  GetChar test is about to call WaitSecs. During execution of WaitSecs, you should type on\n');
            fprintf('  your keyboard.  Press any character key when you are read to begin the test.  \n');
            GetChar;
            fprintf('\n');
            fprintf('  You have ten seconds to type character keys on your keyboard during WaitSecs. Please begin\n');
            fprintf('  typing now!\n');
            FlushEvents('keyDown');
            ListenChar(2);
            WaitSecs(10);
            tChars=[];
            while CharAvail
                tChars=[tChars GetChar];
            end
            FlushEvents('keyDown');
            fprintf(['  You typed ' int2str(length(tChars)) ' characters during WaitSecs\n']);
            fprintf('  The characters which you typed are:\n');
            fprintf(['  ' tChars '\n\n\n']);
        end
        % TEST 7 of 7: Comparing GetChar timing with KbWait timing: Only
        % makes sense with Java-based GetChar, skip test in -nojvm mode.
        if (nargin==0 || (nargin==1 && startAt <= 7))
            fprintf('GetChar test 7 of 7: Comparison of keypress timing reported by KbWait and GetChar\n');
            fprintf('  Timing test: The test will wait five times for a keypress of you. Each time, waiting\n');
            fprintf('  is performed via KbWait. Then the typed character is retrieved via GetChar as well and\n');
            fprintf('  the time stamps reported by KbWait and GetChar are compared to each other. The difference\n');
            fprintf('  between both values will give you a feeling for how much the reported timing is off between\n');
            fprintf('  KbWait and GetChar.\n\n');
            for i=1:5
                fprintf('  Press a character key (a,b,c,d,e,...): ');
                while KbCheck; end;
                FlushEvents('keyDown');
                t1=KbWait;
                [ch, when]=GetChar;
                if isempty(when)
                    break;
                end
                tdelta = 1000 * (t1 - when.secs);
                fprintf('  Pressed %s. Time delta between KbWait and GetChar is %f milliseconds.\n', ch, tdelta);
            end

            fprintf('\n\n  Thanks for typing!\n\n');
        end

        fprintf('\n');
        fprintf('The GetCharTest has reached the end.  Thank you for testing.\n');
        fprintf('\n');
    end
catch
    % This empty catch block makes sure that the following
    % ListenChar(0) command is called in case of error, so Matlab isn't
    % left with a dead keyboard.
end

% Disable listening and flush queue, reenable keyboard dispatching
% to Matlabs windows.
ListenChar(0);

return;

