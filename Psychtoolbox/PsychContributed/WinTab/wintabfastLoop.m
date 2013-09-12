function fastLoop
% function fastLoop
% 
% Opens a PTB Screen, attaches a tablet to the associated screen pointer wPtr and dumps out the contents of an event queue.
%
%
% pkt = WinTabMex(5) reads the latest data point out of a tablet's event queue. This queue is a buffer that can hold around 50 data points.
% This queue begins filling up after a call to WinTabMex(2), which empties the queue, and can be used to record movement data during 
% stimulus presentation in a 'fast loop' (see fastLoop.m). 
%
% pkt is a 8x1 column vector
% 
% tabletTestData.mat:
%           xPos                = pkt(1), x axis position (tablet coordinates)
%           yPos                = pkt(2), y axis position (tablet coordinates)
%           zPos                = pkt(3), z axis position (tablet coordinates)
%           buttonState         = pkt(4), encoded button state (works a little erratically)
%           serialNumber        = pkt(5)
%           tabletTimeStamp     = pkt(6), time in ms from the tablet
%           penStatus           = pkt(7), signals various events (eg penOutOfRange)
%           penChange           = pkt(8), flags what has changed since the last sample
%
%pkt(1:6) are straightforward, pkt(7:8) need some work on figuring out what it all means
%
%Andrew D. Wilson, 2009 (adwkiwi@gmail.com)

KbName('UnifyKeyNames');
spaceKeyID = KbName('space');

data = [];

try
    [wPtr rect] = Screen('OpenWindow', 0);
    WinTabMex(0, wPtr); %Initialize tablet driver, connect it to 'wPtr'   

    
    %PTB things
    Priority(1);
    ListenChar(2);
    HideCursor;
    black = BlackIndex(wPtr);
    white = WhiteIndex(wPtr);
    
    %*******************************************************************************************
    %Loads a window and waits for input to start recording - just for getting set up and ready
    instructions = 'Place the stylus on the tablet and press the space bar to begin recording';
    [instructionsX instructionsY] = centreText(wPtr, instructions, 15);
    Screen('DrawText', wPtr, instructions, instructionsX, instructionsY, black); Screen(wPtr, 'Flip');
    %Waits for space bar
    [keyIsDown, secs, keyCode] = KbCheck;
    while keyCode(spaceKeyID)~=1
        [keyIsDown, secs, keyCode] = KbCheck;
    end    
    %*******************************************************************************************
   
    WinTabMex(2); %Empties the packet queue in preparation for collecting actual data
    information = 'Event queue is filling up...you will be able to retrive up to 50 events (500ms @ 100Hz)';
    [informationX informationY] = centreText(wPtr, information, 15);
    Screen('DrawText', wPtr, information, informationX, informationY, black); 
    Screen(wPtr, 'Flip');
    WaitSecs(3);
    
    %Fast loop to dump all data from the queue into a variable. Each call to WintabMex(5) pulls a data point off the top of the queue, and
    %this loops runs as fast as Matlab can execute it until the queue is drained
    %
    %NB: This loop to drain the queue runs after the data is collected
    while 1
        pkt = WinTabMex(5);

        if isempty(pkt) %ie the queue is empty
            data = [data zeros(9,1)]; %To index that a break got called
            break;
        else
            data = [data pkt];
        end

        status = uint32(pkt(7));
        if bitget(status, 1)
            disp('Pen is outside the tablet');
        end
        if bitget(status, 2)
            disp('Queue overflow - packets getting lost');
        end
    end

    data = data' %Do it this way to suit the way Matlab memory preferences treat rows and columns, etc
    samplingRate = 1000 / mean(diff(data(1:end-1,6)));
    disp(['Your sampling rate appears to be ', int2str(samplingRate), 'Hz.']);

    endInstructions = 'Queue drained...pull the stylus away from the tablet';
    [endX endY] = centreText(wPtr, endInstructions, 15);
    Screen('DrawText', wPtr, endInstructions, endX, endY, black); Screen(wPtr, 'Flip');
    WaitSecs(3);
    
    WinTabMex(3); % Stop/Pause data acquisition.
    WinTabMex(1); % Shutdown driver.
    ListenChar(1);  ShowCursor;
    Screen('CloseAll');
catch
    disp('Quit with error (wintabfastLoop.m)');
    WinTabMex(3); WinTabMex(1);
    ListenChar(1);  ShowCursor;
    Screen('CloseAll');
    
    tabletTestError = lasterror;
    message = tabletTestError(1).message
    identifier = tabletTestError(1).identifier
    stack = tabletTestError(1).stack;
    file = stack.file
    line = stack.line
end

function [x y] = centreText(wPtr, text, preferredFontSize)
%function [x y] = centreText(wPtr, text, preferredFontSize)
%
%For use with the Psychtoolbox. 
%Takes a window pointer & a string and returns the x,y coordinates that will centre the text if fed to Screen('DrawText')
%
%ARGUMENTS:
%wPtr = the window pointer returned by a call to Screen('OpenWindow')
%text = a string to be sent to the screen
%preferredFontSize = the maximum sized font you would like to have. If the text is too big, centreText will cycle through 
%   text sizes until it works and set it to the biggest that fits
%
%RETURN VALUES:
%[x y] = a vector with the (x,y) coordinates to centre that text.
%
%NB Requires that you have used Screen('OpenWindow')
%
%Andrew D. Wilson (v1.0.3; 11 April 2008)

x=-1; y=-1; %Allows loop to run once

while x<0
    Screen('TextSize', wPtr, preferredFontSize);

    [normBoundsRect, offsetBoundsRect]= Screen('TextBounds', wPtr, text);
    rect = Screen('Rect', wPtr);

    windowCentre = [rect(3)/2 rect(4)/2];
    textCentre = [normBoundsRect(3)/2 normBoundsRect(4)/2];

    x = windowCentre(1) - textCentre(1);
    y = windowCentre(2) - textCentre(2);
    
    if x < 0 || y < 0 %ie if the text ends up being drawn offscreen
        preferredFontSize = preferredFontSize-1;
    end
end

return;
    
