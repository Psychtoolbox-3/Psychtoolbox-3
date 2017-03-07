function SetMouse(x,y,windowPtrOrScreenNumber, mouseid, detachFromMouse)
% SetMouse(x, y [, windowPtrOrScreenNumber][, mouseid][, detachFromMouse=0])
% 
% Position the mouse cursor on the screen.
%
% The cursor position (x,y) is "local" to the screen, i.e. relative to the
% origin of the screen if a screen number is supplied, or relative to the
% origin of a screen on which a supplied onscreen window is displayed.
% Otherwise it's "global", i.e. relative to the origin of the main screen
% (Screen 0).
%
% On Linux with X11, the optional 'mouseid' parameter allows to select
% which of potentially multiple cursors should be repositioned. On OS/X and
% Windows this parameter is silently ignored.
%
% On OSX, the optional 'detachFromMouse' parameter, if set to 1 instead of
% its default value of zero, will detach mouse cursor movements from mouse
% movements. The cursor will be frozen in place and can only be moved via
% SetMouse(), but no longer via mouse movements. On Linux and Windows this
% parameter is currently silently ignored.
%
% On Linux with the Wayland backend, this function does nothing.
%
% The delay between a call to SetMouse and when GetMouse will report the
% new mouse cursor position is not known. GetMouse seems to report the new
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
% 11/04/14  mk      round() x,y coords for integral coordinates to avoid error.
% 04/18/15  mk      Update help text - local coordinates now should also
%                   work on MS-Windows, not only Linux and OSX.
% 07/20/15  mk      Add support for 'detachFromMouse', on OSX for now.

% SetMouse.m wraps the Screen('PositionCursor',..) call to emulate the old SetMouse.mex

if nargin < 2
  error('SetMouse requires x and y positions');
end

if nargin < 3 || isempty(windowPtrOrScreenNumber)
  windowPtrOrScreenNumber = 0;
end

if nargin < 4
  mouseid = [];
end

if nargin < 5 || isempty(detachFromMouse)
  detachFromMouse = 0;
end
if isunix && ~ismac
    Screen('SetMouseHelper', windowPtrOrScreenNumber, round(x), round(y), mouseid);
else
    Screen('SetMouseHelper', windowPtrOrScreenNumber, round(x), round(y), mouseid, detachFromMouse);
end
