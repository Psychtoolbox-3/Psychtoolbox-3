function [strInputName, x, y] = bvlWaitForInput(scal,fSeconds, fBufferWindow)
% [strInputName, x, y] = bvlWaitForInput(fSeconds [, fBufferWindow])
%
% Wait for a key press or mouse input, then return the input name and mouse position.
% strInputName is a string identifier of the key. (Defined in KbName)
% x, y are the pixel coordinates of the mouse
%
% Returns after fSeconds, if no key is pressed.
% If fSeconds is not specified, returns after 1.0 seconds
%
% fBufferWindow allows you to specify the time buffer (window) between
% key/mouse events that will be accepted. This way we don't register the
% key-down and key-up of the same key as two events.  Nor do we register
% one key press as multiple events.  If this parameter is not passed in,
% the default is 0.250 (250 milliseconds).  Reduce this if you need
% constant feedback, for instance, in tracing out a rectangle when
% selecting points in the calibration, I use 0.010 seconds.
%
% Some example strInputNames:
%     LeftArrow
%     RightArrow
%     UpArrow
%     DownArrow
%     ESCAPE
%     LeftMouse
%     MiddleMouse
%     RightMouse
%
% 2007-05-15  - cburns - update to fix bugs in events not being handled.
% 2007-05-10  - cburns - Rewrote for include mouse input and update
%                        for Psychtoolbox 3.x install.
%
% 2006-01-31  - cburns - Bankslab, UC Berkeley
%

global gTimeOfPrevKeyEvent

% bv180108 Try to find out what the current display window is. If there are
% more than 1 active windows, this method won;t work and we have to proceed
% without information about the presentation window (indicated by an empty
% variable windowlist)


    
    


if isempty(gTimeOfPrevKeyEvent)
    % initialize the global timestamp otherwise we'll never enter the
    % event-handling blocks below and set a new timestamp
    gTimeOfPrevKeyEvent = GetSecs;
end

if (nargin < 2)
    fSeconds = 1.0;
    windowlist=0;
else
    windowlist = scal.windowPtr;
end

if (nargin < 3)
    fBufferWindow = 0.250;
end

% initialize return variables
strInputName = '';
x = 0;
y = 0;

fCurrTime = GetSecs;
fEndTime = fCurrTime + fSeconds;

while (fCurrTime < fEndTime)

    % Check state for entire keyboard
    % iKeyCode is an array of size 256
    [iKeyIsDown, fSecs, iKeyboardCodes] = KbCheck(-1);

    % Handle keyboard input
    if iKeyIsDown
        % We want to allow people to hold down keys for repeated behavior
        % but we need to dampen the speed of the input.  So grab the events
        % as they come in but only process them ever 250ms.
        % A hacked way of flushing the event queue.
        % This is a timed version of the "while KbCheck; end" used in the
        % KbDemo.m in the PTB.
        if fSecs - gTimeOfPrevKeyEvent > fBufferWindow
            iKeyIndex = find(iKeyboardCodes);
            % For Debugging:
            %KbName(iKeyIndex)

            % Get the keyboard name from PTB.
            strInputName = KbName(iKeyIndex);

            % BUGGY:  On the same WindowsXP development machine, using Matlab
            % 7.0.1 and PTB 3.0.8, I received both the Windows keynames and the
            % OSX keynames.... I don't know why.  Map them to one here so we
            % can count on the same from the calling function.
            if strcmp(strInputName, 'left')
                strInputName = 'LeftArrow';
            elseif strcmp(strInputName, 'right')
                strInputName = 'RightArrow';
            elseif strcmp(strInputName, 'up')
                strInputName = 'UpArrow';
            elseif strcmp(strInputName, 'down')
                strInputName = 'DownArrow';
            elseif strcmp(strInputName, 'esc')
                strInputName = 'ESCAPE';
            end

            % event handled break out of loop
            gTimeOfPrevKeyEvent = fSecs;
            break;
        end
    end

    % No keyboard input, Check mouse state. 
    % bv180108 If there is information about
    % the presentation window, then we will retrieve the mouse position
    % relative to that window. 
    if isempty(windowlist),
        [xMouse, yMouse, buttons] = GetMouse;
    else 
        [xMouse, yMouse, buttons] = GetMouse(windowlist);     
    end
    x = xMouse;
    y = yMouse;

    % Handle mouse input
    if buttons(1) || buttons(2) || buttons(3)
        % Just as for the keyboard input, only accept new input every
        % 250ms.  This allows someone to hold down the mouse button and
        % cycle through items at a reasonable pace
        if fSecs - gTimeOfPrevKeyEvent > fBufferWindow
            if IsWin || IsLinux
                if buttons(1)
                    strInputName = 'LeftMouse';
                elseif buttons(2)
                    strInputName = 'MiddleMouse';
                elseif buttons(3)
                    strInputName = 'RightMouse';
                end
            end
            
            if IsOSX
                if buttons(1)
                    strInputName = 'LeftMouse';
                elseif buttons(3)
                    strInputName = 'MiddleMouse';
                elseif buttons(2)
                    strInputName = 'RightMouse';
                end
            end
            
            % event handled break out of loop
            gTimeOfPrevKeyEvent = fSecs;
            break;
        end
    end

    % No need to check more then every 5 millisecond
    WaitSecs('YieldSecs', 0.005);
    fCurrTime = GetSecs;
end
