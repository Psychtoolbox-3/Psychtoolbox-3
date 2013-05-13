% MouseTraceDemo2
%
% ___________________________________________________________________
% 
% Draw a curve with the mouse. Same as MouseTraceDemo, but asks
% Screen('Flip') to not clear the framebuffer after flip. This way,
% we don't need to redraw the whole mousetrace in each frame.
% _________________________________________________________________________
%
% See also: PsychDemos, MouseTraceDemo, GetMouse.
%
% HISTORY
%                       
% 4/23/05  mk       Derived from MouseTraceDemoOSX: Uses new "Don't clear" mode of Flip.
%                   to avoid redrawing the whole past mousetrace after each
%                   Flip --> Faster.

try
    % Open up a window on the screen and clear it.
    whichScreen = max(Screen('Screens'));
    [theWindow,theRect] = Screen(whichScreen,'OpenWindow',0);

    % Move the cursor to the center of the screen
    theX = theRect(RectRight)/2;
    theY = theRect(RectBottom)/2;
    SetMouse(theX,theY);

    % Wait for a click and hide the cursor
    Screen(theWindow,'FillRect',0);
    Screen(theWindow,'TextSize',24);
    Screen(theWindow,'DrawText','Drag mouse (i.e. hold button down) to draw',50,50,255);
    Screen('Flip', theWindow);
    while (1)
        [x,y,buttons] = GetMouse(theWindow);
        if buttons(1)
          break;
        end
    end
    Screen(theWindow,'DrawText','Release button to finish',50,50,255);

    HideCursor;

    % Loop and track the mouse, drawing the contour
    [theX,theY] = GetMouse(theWindow);
    thePoints = [theX theY];
    Screen(theWindow,'DrawLine',255,theX,theY,theX,theY);
    % Set the 'dontclear' flag of Flip to 1 to prevent erasing the
    % frame-buffer:
    Screen('Flip', theWindow, 0, 1);
    while (1)
        [x,y,buttons] = GetMouse(theWindow);	
        if ~buttons(1)
            break;
        end
        if (x ~= theX || y ~= theY)
            thePoints = [thePoints ; x y]; %#ok<AGROW>
            [numPoints, two]=size(thePoints);
            % Only draw the most recent line segment: This is possible,
            % because...
            Screen(theWindow,'DrawLine',128,thePoints(numPoints-1,1),thePoints(numPoints-1,2),thePoints(numPoints,1),thePoints(numPoints,2));
            % ...we ask Flip to not clear the framebuffer after flipping:
            Screen('Flip', theWindow, 0, 1);
            theX = x; theY = y;
        end
    end

    % Close up
    Screen(theWindow,'DrawText','Click mouse to finish',50,50,255);
    ShowCursor;
    Screen(theWindow,'Close');

    % Plot the contour in a Matlab figure
    plot(thePoints(:,1),theRect(RectBottom)-thePoints(:,2));
catch
    Screen('CloseAll')
    ShowCursor;
    psychrethrow(psychlasterror);
end %try..catch..
