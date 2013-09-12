function slowLoop(samplingRate)
% function slowLoop(samplingRate)
% 
% Opens a PTB Screen, attaches a tablet to the associated screen pointer wPtr and records data from a graphics tablet. The data is sorted 
% and output to a .mat file and an Excel file.
%
% INPUTS: 
% samplingRate: sampling rate in Hz for your tablet. Establish this using fastLoop.m
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
%           serialNumber        = pkt(5), an index of which data point this is since the last call to WinTabMex(2)
%           tabletTimeStamp     = pkt(6), time in ms from the tablet
%           penStatus           = pkt(7), signals various events (eg penOutOfRange)
%           penChange           = pkt(8), flags what has changed since the last sample
%           getsecTimeStamp     = the time at which the data was collected from the start of the trial, using PTB's GetSecs function
%           pktData             = matrix of data, compiled into columns for writing to Excel
%
%pkt(1:6) are straightforward, pkt(7:8) need some work on figuring out what it all means
%
%Andrew D. Wilson, 2009 (adwkiwi@gmail.com)

trialLength = 5;
deltaT = 1/samplingRate;

%Set up a variable to store data
pktData = [];

KbName('UnifyKeyNames');
spaceKeyID = KbName('space');

try
    wPtr = Screen('OpenWindow', 0); 
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
    
    record = 'The tablet is currently recording your movements...';
    [recordX recordY] = centreText(wPtr, record, 15);
    Screen('DrawText', wPtr, record, recordX, recordY, black); Screen(wPtr, 'Flip');

    WinTabMex(2);       %Empties the packet queue in preparation for collecting actual data
                        %Call this immediately before beginning any slow loop
            
    %This loop runs for trialLength seconds. 
    start = GetSecs;
    stop  = start + trialLength;
    while GetSecs<stop
        loopStart = GetSecs;
      
        %This loop runs for deltaT or until it successfully retrieves some data from the queue
        while 1  %Note this loop MUST be broken manually, as 'while 1' always returns TRUE
            pkt = WinTabMex(5);
            
            %This check breaks the loop if data is recovered from the queue before deltaT is up
            if ~isempty(pkt)
                break
            end
            
            %This check breaks the loop after deltaT if pkt was always empty
            if GetSecs>(loopStart+deltaT)
                pkt = zeros(9,1); %Dummy data representing a missed data point
                break;
            end
        end
        pkt = [pkt; (GetSecs - start)];
        pktData = [pktData pkt];
        
        %Waits to end of deltaT if need be
        if GetSecs<(loopStart+deltaT)
            WaitSecs('UntilTime', loopStart+deltaT);
        end
    end
    pktData = pktData';  %Assemble the data and then transpose to arrange data in columns because of Matlab memory preferences
    
    WinTabMex(3); % Stop/Pause data acquisition.
    Priority(0);

    %Sorts data and outputs to a save file
    if ~isempty(pktData)
        xPos = pktData(:,1);
        yPos = pktData(:,2);
        zPos = pktData(:,3);
        buttonState = pktData(:,4);
        serialNumber = pktData(:,5);
        tabletTimeStamp = pktData(:,6);
        penStatus = uint32(pktData(:,7));
        penChange = uint32(pktData(:,8));
        getsecTimeStamp = pktData(:,9);
        
        save tabletTestData pktData xPos yPos zPos buttonState serialNumber tabletTimeStamp penStatus penChange getsecTimeStamp
        
        xlswrite('tabletTestData', [{'x Position'} {'y Position'} {'z Position'} {'Button State'} {'Serial Number'} {'Tablet Timestamp'},...
            {'Pen Status'} {'Pen Change'} {'PTB Timestamp (GetSecs)'}], 'Sheet1');
        xlswrite('tabletTestData', pktData, 'Sheet1', 'A2');
        
        disp('Your data are in two files: tabletTestData.mat and tabletTestData.xls');
    else
        disp('Error: variable ''pktData'' was empty. (tabletTest.m)');
    end
    
    %Closing out stuff
    WinTabMex(1); % Shutdown driver.
    ListenChar(1);  ShowCursor;
    Screen('CloseAll');    
catch
    disp('Quit with error (tabletTest.m)');
    
    %Close tablet things
	WinTabMex(3); WinTabMex(1); %Close out data acquisition and the handle on the driver
    
    %Close PTB things
    ListenChar(1); ShowCursor; Priority(0); Screen('CloseAll');
    
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
