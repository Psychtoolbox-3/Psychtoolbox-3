function EyelinkEventExample
% Short MATLAB example program to demonstrate the use of events
% with the Eyelink and Psychophysics Toolboxes. When subject gazes at an object
% it changes color. When gaze moves away, it toggles back
%

% history
% 28-10-02	fwc	created it, based on eyelinkexample
% 27-11-02	fwc	changed dialog for dummy mode question, fixed dodriftcorrection bug
% 07-07-10  fwc adapted to use new psychtooblox/eyelink functions, finally
%           got toggling correct too ;-)
%
dummymode=0;       % set to 1 to run in dummymode (using mouse as pseudo-eyetracker)
showboxes=1;        % set to 1 to show boxes within which you have to fixate in order to toggle a number

try
    
    fprintf('EyelinkToolbox Event Example\n\n\t');
    
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
    
    % make sure that we get event data from the Eyelink
%         Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
    Eyelink('command', 'link_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY');
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,BLINK,SACCADE,BUTTON');
    
    % open file to record data to
    edfFile='demo.edf';
    Eyelink('Openfile', edfFile);
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    
    % STEP 5
    % draw a stimulus in a set of buffers
    [buffer, altbuffer, object]=DrawStimulus(window, showboxes);
    
    Screen('TextFont', buffer, el.msgfont);
    Screen('TextSize', buffer, el.msgfontsize);
    Screen('TextFont', altbuffer, el.msgfont);
    Screen('TextSize', altbuffer, el.msgfontsize);
    [width, height]=Screen('WindowSize', window);
    message='Press space to stop.';
    Screen('DrawText', buffer, message, 200, height-el.msgfontsize-20, el.msgfontcolour);
    Screen('DrawText', altbuffer, message, 200, height-el.msgfontsize-20, el.msgfontcolour);

    
    % STEP 6
    % do a final check of calibration using driftcorrection
    success=EyelinkDoDriftCorrection(el);
    if success~=1
        cleanup;
        return;
    end
    
    % STEP 7
    % some final preparations
    stopkey=KbName('space');
    lastchoice=-1;
    choice=-1;
    doflip=0;
    % put buffer on screen
    Screen('CopyWindow', buffer, window);
    Screen('Flip',  window, [], 1); % don't erase buffer
    
    % start recording eye position
    Eyelink('startrecording');
    % record a few samples before we actually start displaying
    waitsecs(0.1);
    % mark zero-plot time in data file
    Eyelink('message', 'SYNCTIME');
    
    % STEP 8
    % adapt display based on END-SACCADE events
    while 1 % loop till error or space bar is pressed
        error=Eyelink('checkrecording');        % Check recording status, stop display if error
        if(error~=0)
            break;
        end
        
        [keyIsDown,secs,keyCode] = KbCheck; % check for keyboard press
        if keyCode(stopkey) % if spacebar was pressed stop display
            break;
        end
        
        % check for endsaccade events
        if Eyelink('isconnected') == el.dummyconnected % in dummy mode use mousecoordinates
            [x,y,button] = GetMouse(window);
            evt.type=el.ENDSACC;
            evt.genx=x;
            evt.geny=y;
            evtype=el.ENDSACC;
        else % check for events
            evtype=Eyelink('getnextdatatype');
        end
        if evtype==el.ENDSACC		% if the subject finished a saccade check if it fell on an object
            if Eyelink('isconnected') == el.connected % if we're really measuring eye-movements
                evt = Eyelink('getfloatdata', evtype); % get data
            end
            % check if saccade landed on an object
            choice=-1;
            noobject=0;
            i=1;
            while 1
                if 1==IsInRect(evt.genx,evt.geny, object(i).rect )
                    choice=i;
                    break;
                end
                i=i+1;
                if i>length(object)
                    noobject=1;
                    break;
                end
            end
            if lastchoice>0 && (choice~=lastchoice || noobject==1) % toggle object color
                if object(lastchoice).on==1 % restore screen
                    Screen('CopyWindow', buffer, window, object(lastchoice).rect, object(lastchoice).rect);
                    object(lastchoice).on=0;
                    lastchoice=-1;
                    doflip=1;
                end
            end
            if choice>0 && choice~=lastchoice % toggle object color
                if object(choice).on==0 % toggle object on screen
                    Screen('CopyWindow', altbuffer, window, object(choice).rect, object(choice).rect);
                    object(choice).on=1;
                    doflip=1;
                end
                lastchoice=choice;
            end
            if doflip==1
                Screen('Flip',  window, [], 1);
                doflip=0;
            end
        end % saccade?
    end % main loop
    
    % STEP 9
    % finish demo
    waitsecs(0.1);	% wait a while to record a few more samples 
    cleanup;
    
catch myerr
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    cleanup;
    psychrethrow(psychlasterror);
end %try..catch.


% Cleanup routine:
function cleanup

% finish up: stop recording eye-movements,
% close graphics window, close data file and shut down tracker
% Eyelink('Stoprecording');
% Eyelink('CloseFile');
Eyelink('Shutdown');

% Close window:
sca;

% Restore keyboard output to Matlab:
ListenChar(0);

function [buffer, altbuffer, object]=DrawStimulus(window, showboxes)

if ~exist('showboxes', 'var') || isempty(showboxes)
    showboxes=0;
end
white=WhiteIndex(window);
black=BlackIndex(window);
gray=GrayIndex(window);
Screen(window, 'FillRect',gray);
[w,h]=WindowSize(window);
[x,y]=WindowCenter(window);
fontsize=round(w/100*5);
oldFont=Screen(window,'TextFont','Arial');
oldFontSize=Screen(window,'TextSize',fontsize);

buffer=Screen(window, 'OpenOffscreenWindow', gray);
altbuffer=Screen(window, 'OpenOffscreenWindow', gray);

cols=5;
rows=4;

Screen(buffer,'TextFont','Arial');
Screen(altbuffer,'TextFont','Arial');
Screen(buffer,'TextSize',fontsize);
Screen(altbuffer,'TextSize',fontsize);
k=1;
xdist=round(w/(cols));
ydist=round(h/(rows));
x0=x-((cols-1)/2)*xdist;
y0=y-((rows-1)/2)*ydist;
for i=1:rows
    for j=1:cols
        text=num2str(k);
        xpos=round(x0+(j-1)*xdist-fontsize/2);
        ypos=round(y0+(i-1)*ydist-fontsize/2);
        Screen(buffer,'DrawText',text,xpos,ypos,white);
        Screen(altbuffer,'DrawText',text,xpos,ypos,black);
        rect= Screen('TextBounds', window, text, xpos,ypos ); %[,yPositionIsBaseline] [,swapTextDirection]);
        object(k).rect=OffsetRect(rect, xpos, ypos);
        if showboxes==1
            Screen('FrameRect', buffer, black, object(k).rect);
            Screen('FrameRect', altbuffer, white, object(k).rect);
        end
        object(k).on=0;
        k=k+1;
    end
end

Screen(window,'TextFont',oldFont);
Screen(window,'TextSize',oldFontSize);
