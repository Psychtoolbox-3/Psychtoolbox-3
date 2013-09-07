function white=WhiteIndex(windowPtrOrScreenNumber)
% color=WhiteIndex(windowPtrOrScreenNumber)
% Returns the intensity value to produce white at the current screen depth,
% assuming a standard color lookup table for that depth. E.g.
%
% white=WhiteIndex(w);
% Screen(w,'FillRect',white);
% 
% See BlackIndex.
% 

% 3/10/98	dgp Wrote it.
% 3/30/98	dgp Consider only one channel, even for 16 and 32 bit modes.
% 3/8/2000  emw Added Platform Conditionals
% 3/8/2000	dgp Fixed platform conditionals
% 3/30/2004 awi Added OS X case. For now OS X only supports true-color mode, so
%               WhiteIndex behavior on OS X will have to change when we add
%               more depth modes.
% 1/29/05   dgp Cosmetic.
% 03/1/08   mk  Adapted to the much more flexible scheme of PTB-3.
% 08/24/13  mk  Select 1.0 as default white index if PsychDefaultSetup(2+) was used.

% Default colormode to use: 0 = clamped, 0-255 range. 1 = unclamped 0-1 range.
global psych_default_colormode;

if nargin~=1
	error('Usage: color=WhiteIndex(windowPtrOrScreenNumber)');
end

% Is a default colormode specified via psych_default_colormode variable and the level
% is at least 1? If so, future created onscreen windows will have a [0;1] colorrange
% without clamping by default.
if ~isempty(psych_default_colormode) && (psych_default_colormode >= 1)
    % 0-1 normalized range preset. A default white index of 1.0 is a reasonable assumption,
    % as PsychImaging('Openwindow') would select 1.0 as maximum value:
    defaultwhite = 1;
else
    % No preset: Assume a traditional default of 255:
    defaultwhite = 255;
end

% Screen number given?
if ~isempty(find(Screen('Screens')==windowPtrOrScreenNumber))
    % Need to find corresponding onscreen window:
    windows = Screen('Windows');
    
    if isempty(windows)
        % No open windows associated with this screen. Just return the
        % default value of "defaultwhite", our default maximum pixel color component
        % value, which is valid irrespective of the actual pixel depths of
        % the screen as OpenGL takes care of such things / is invariant to
        % them:
        white = defaultwhite;
        return;
    end
    
    % At least one onscreen window open: Find the one with this screen as
    % parent:
    win = [];
    for i=windows
        if (Screen('WindowKind', i) == 1)
            % It's an onscreen window. Associated with this screen?
            if windowPtrOrScreenNumber == Screen('WindowScreenNumber', i)
                % This is it:
                win = i;
                break;
            end
        end
    end
    
    if isempty(win)
        % No onscreen window on this screen. Return default "defaultwhite":
        white = defaultwhite;
        return;
    end
    
    % Onscreen window id assigned to 'win'. Leave the rest of the job to
    % common code below...
else
    % No screen number given. Window number given?
    if isempty(find(Screen('Windows')==windowPtrOrScreenNumber))
        % No window number either. This is an invalid index:
        error('Provided "windowPtrOrScreenNumber" is neither a valid screen, nor window!');
    end

    % Its a window: Assign it to 'win'
    win = windowPtrOrScreenNumber;
end

% If we reach this point then we have the window handle of the relevant
% window to query in 'win'. use Screen('ColorRange') to query its maximum
% color value. By default this will be again "255" - the maximum for a 8bpc
% standard framebuffer. However, when used with special HDR
% devices/framebuffers or some specific setup was done via
% Screen('ColorRange'), the maximum value corresponding to
% white may be any positive number:
white = Screen('ColorRange', win);

% Done.
return;
