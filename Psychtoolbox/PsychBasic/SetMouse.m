function SetMouse(x,y,windowPtrOrScreenNumber, mouseid)
% SetMouse(x, y [, windowPtrOrScreenNumber][, mouseid])
% 
% Position the mouse cursor on the screen.
%
% The cursor position (x,y) is "local", i.e. relative to the origin of the
% window or screen, if supplied. Otherwise it's "global", i.e. relative to
% the origin of the main screen (the one with the menu bar).
%
% On Linux, the optional 'mouseid' parameter allows to select which
% of potentially multiple cursors should be repositioned. On OS/X and
% Windows this parameter is silently ignored.
% 
% OS 9: ___________________________________________________________________
%
% After calling SetMouse, there's a short period (until the next tick?)
% during which GetMouse reports the old position, after which it reports
% the new position. To be conservative, you can code as follows:
%  while 1
%    SetMouse(theX,theY);
%    [checkX,checkY] = GetMouse;
%    if (checkX==theX) && (checkY==theY)
%	   break;
%    end
%  end
%
% Cursor updating is usually carried out by a System task that runs once
% per tick, so it's likely that SetMouse doesn't take effect until the
% next tick. Instead of the elaborate check suggested above, it might
% be enough, after calling SetMouse, to simply WaitTicks(1) before calling
% GetMouse to be sure of getting the new position.
%
% OS-X, Linux & WINDOWS:___________________________________________________
%
% Psychtoolbox will accept the optional windowPtrOrScreenNumber
% argument and check it for validity.  However, supplying the argument will
% not  influence the position of the mouse cursor.  The cursor is always
% positioned  in absolute coordinates on the main screen.  It does not
% matter that the  cursor is always positioned in absolute coordinates,
% because in Windows  onscreen windows are always the full size of the
% display, absolute and relative coordinates are always the same.  It would
% be good if SetMouse could position the mouse cursor on secondary
% displays, but we dont' support  multiple displays yet in Windows.
%
% The delay between a call to SetMouse and when GetMouse will report the
% new mouse cursor position is not known.  GetMouse seems to report the new
% position immediately, but we have no guarantee that it always will.
%
% _________________________________________________________________________
%
% See Also: GetMouse, GetClicks

% 6/7/96    dhb     Wrote it.
% 8/23/96   dhb     Added support for windowInfo argument.
% 3/23/97   dgp     Updated.
% 8/14/97   dhb     Added comment about delay.
% 8/15/97   dgp     Suggest WaitTicks(1).
% 4/24/01   awi     Added WINDOWS section.
% 6/10/01   awi     Added See Also.  
% 4/14/03   awi     ****** OS X-specific fork from the OS 9 version *******
%                   Added call to Screen('PositionCursor'...) for OS X.
% 10/12/04  awi     Cosmetic changes to help.  This file should be modified 
%                   after the great mouse shift to state SetMouse is depricated.
% 11/18/04  awi     Renamed "PositionCursor" to "PositionCursorHelper".
% 02/21/06  mk      Added Linux support.
% 06/17/06  mk      Added Windows support.

% SetMouse.m wraps the Screen('PositionCursor',..) call to emulate the old SetMouse.mex

if nargin < 2
   error('SetMouse requires x and y positions');
end

if nargin < 3
   windowPtrOrScreenNumber = 0;
end

if nargin < 4
   mouseid = [];
end

Screen('SetMouseHelper', windowPtrOrScreenNumber, x, y, mouseid);
