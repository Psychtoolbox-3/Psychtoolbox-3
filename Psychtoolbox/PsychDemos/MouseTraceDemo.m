% MouseTraceDemo
% ___________________________________________________________________
% 
% Draw a curve with the mouse.
% _______________________________________________________________________
%
% See also: PsychDemos, MouseTraceDemo2, GetMouse.

% HISTORY
% 8/12/97  dhb, wtf  Wrote it.
% 8/13/97  dhb		 Small modifications.
% 8/15/97  dgp	     Drag mouse instead of clicking before and after.
% 8/07/01  awi       Added font conditional, changed "button" to "buttons"
%                    and indexed the mouse button result for Windows.
% 4/11/01  awi		 Cosmetic editing of comments.
% 4/13/02  dgp       Use Arial, no need for conditional.
% 11/18/04 awi       Modified to make it work on OS X and renamed to MouseTraceDemoOSX.  
%                    Changed the open command to specify double
%                        buffers and 32 depth because we don't yet support
%                        8-bit depth on OS X.  
%                    Double buffer the display;  
%                        the path does not accumulate in the window so
%                        render the entire path on every frame, not just the
%                        path delta since the last segment was drawn.  
%                    Added flip command because we are double
%                        buffered.
%                    Added try..catch to close onscreen window in
%                        event of failure.  
%                    Open on the highest-numbered display and pass
%                        GetMouse the window pointer.  Not necessary for OS
%                        X, but often desirable to use secondary display if
%                        available.                
% 11/19/06 dhb       Remove OSX from name.

try
    % Open up a window on the screen and clear it.
    whichScreen = max(Screen('Screens'));
    [theWindow,theRect] = Screen(whichScreen,'OpenWindow',0,[],[],2);

    % Move the cursor to the center of the screen
    theX = theRect(RectRight)/2;
    theY = theRect(RectBottom)/2;
    SetMouse(theX,theY);
    ShowCursor(0);

    % Wait for a click and hide the cursor
    Screen(theWindow,'FillRect',0);
    %Screen(theWindow,'TextFont','Arial');
    %Screen(theWindow,'TextSize',18);
    Screen(theWindow,'DrawText','Drag mouse (i.e. hold button down) to draw',50,50,255);
    Screen('Flip', theWindow);
    while (1)
        [x,y,buttons] = GetMouse(whichScreen);
        if buttons(1) || KbCheck
          break;
        end
    end
    Screen(theWindow,'DrawText','Drag mouse (i.e. hold button down) to draw',50,50,0);
    %Screen(theWindow,'DrawText','Release button to finish',50,50,255);

    % Loop and track the mouse, drawing the contour
    [theX,theY] = GetMouse(whichScreen);
    thePoints = [theX theY];
    Screen(theWindow,'DrawLine',255,theX,theY,theX,theY);
    Screen('Flip', theWindow);
    sampleTime = 0.01;
    startTime = GetSecs;
    nextTime = startTime+sampleTime;
    while (1)
        [x,y,buttons] = GetMouse(whichScreen);	
        if ~buttons(1)
            break;
        end
        if (x ~= theX || y ~= theY)
            [numPoints, two]=size(thePoints);
            for i= 1:numPoints-1
                Screen(theWindow,'DrawLine',128,thePoints(i,1),thePoints(i,2),thePoints(i+1,1),thePoints(i+1,2));
            end
            Screen('Flip', theWindow);
            theX = x; theY = y;
        end
        if (GetSecs > nextTime)
            thePoints = [thePoints ; x y];
            nextTime = nextTime+sampleTime;
        end
    end

    % Close up
    Screen(theWindow,'DrawText','Click mouse to finish',50,50,255);
    ShowCursor;
    Screen(theWindow,'Close');

    % Plot the contour in a Matlab figure
    plot(thePoints(:,1),theRect(RectBottom)-thePoints(:,2));
catch
    Screen('CloseAll');
    ShowCursor;
    psychrethrow(psychlasterror);
end %try..catch..
