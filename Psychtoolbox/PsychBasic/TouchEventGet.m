function [event, nremaining] = TouchEventGet(deviceIndex, windowHandle, maxWaitTimeSecs)
% [event, nremaining] = TouchEventGet(deviceIndex, windowHandle [, maxWaitTimeSecs=0])
%
% Return oldest pending event, if any, in return argument 'event', and the
% remaining number of recorded events in the event buffer of a touch
% queue in the return argument 'nremaining'.
%
% TouchEventGet() will wait up to 'maxWaitTimeSecs' seconds for at least one
% event to show up before it gives up. By default, if 'maxWaitTimeSecs' is 0
% or empty or omitted, it doesn't wait but just gives up if there aren't any
% events queued at time of invocation. In that case it returns an empty 'event'.
%
% 'event' is either empty if there aren't any events available, or it is a
% struct with various information about the touch event. The returned event
% struct currently contains at least the following useful fields:
%
% 'Type' = 0 for button presses or releases if the touch device also has
%          physical or logical buttons. 'Pressed' will tell if this is a
%          button press (1) or release (0) event.
%
%        = 1 single touch point move: Only happens on some mousepads if mouse
%          emulation kicks in. 'X', 'mappedX' and 'Y' / 'mappedY' will tell the
%          touch point / emulated mouse cursor position.
%
%        = 2 New touch: 'X/NormX/mappedX' and 'Y/NormY/mappedY' tell start location.
%            Happens, e.g., when a finger touches the touch surface. The 'Keycode'
%            contains a numeric id which uniquely identifies this touch point while
%            it is active, e.g., while the finger stays on the surface.
%
%        = 3 Touch moved: 'X/NormX/mappedX' and 'Y/NormY/mappedY' tell new location.
%            E.g., when a finger moves over the touch surface. 'Keycode' allows you
%            to know which of possible multiple fingers or tools is moving.
%
%        = 4 Touch finished: 'X/NormX/mappedX' and 'Y/NormY/mappedY' tell final location.
%            Happens when the finger or tool is lifted from the touch surface. The
%            'Keycode' tells you which finger or tool was lifted from the surface.
%
%        = 5 Touch data lost: Some other application or the GUI took over our
%            touch input, and cut us off from it, so the recorded touch data
%            sequence is incomplete. You should wait a second, then call
%            TouchEventFlush() to discard all remaining events in the queue,
%            then mark your trial invalid and take some corrective action,
%            maybe asking the experimenter to disable other running touch
%            applications or GUI services, e.g., touch gesture recognition of
%            the system user interface, so they can no longer interfere with
%            experiment data collection.
%
% 'X'    = x-position in units of screen pixels, with fractional (sub-pixel) resolution.
% 'Y'    = y-position in units of screen pixels, with fractional (sub-pixel) resolution.
%
% The origin (0,0) for 'X' and 'Y' coordinates is the top-left corner of the screen
% on MS-Windows and on Linux with classic X11 X-Window display system. On Linux with
% the next-generation Wayland display system, (0,0) refers to the top-left corner of
% the associated onscreen window for the touch event. Use 'mappedX' and 'mappedY'
% below if you always want onscreen window relative coordinates, transformed to the
% top-left corner of the provided 'windowHandle'.
%
% 'NormX' = Normalized x-position in range 0 - 1. Constant 0 if unavailable.
% 'NormY' = Normalized y-position in range 0 - 1. Constant 0 if unavailable.
%
% 'MappedX' = x-position relative to provided onscreen window with 'windowHandle'.
% 'MappedY' = y-position relative to provided onscreen window with 'windowHandle'.
%
% 'Time' = The GetSecs time when the event was received.
%
% 'Keycode' = The unique numeric id key of this touch point. A specific physical
%             touch, e.g., one specific finger touching a touchscreen, will get
%             a unique number assigned, which will persist while the finger rests
%             or moves on the touch surface, until it is lifted off the surface,
%             ie. the unique number is assigned in a Type=2 event, kept throughout
%             Type=3 events, and last used in Type=4 events, before it goes the
%             way of all mortal things. If you'd put the same finger down onto the
%             touch surface again, it would get a new unique number assigned, as
%             the touch hardware doesn't know it is the same finger.
%
% 'Pressed' = 1 while a touch point is active (finger touches screen), 0 when the
%             touch point is removed, e.g., finger lifted from screen.
%
% 'Motion'  = 1 while a touch point (=finger) is moving over the surface, 0 while
%             it is resting.
%

% NOTE:
%
% Linux/X11 uses absolute touch device coordinates, which usually resolve much
% finer than display pixels, and can be mapped into absolute screen pixels, or
% window relative pixels with ease.
%
% Linux/Wayland only provides surface local coordinates, iow. the absolute pixels
% location on the screen is unknown, unless a fullsreen window is used where
% window == screen.
%
% Windows uses absolute screen coordinates in 1/100th pixels, so can be remapped
% with ease.
%
% macOS doesn't know what a touchscreen is.
%
% So X and Y is not fully portable, but mappedX and mappedY or normX and normY
% can be computed in a reasonably meaningful way on all systems.
%

% 1-Oct-2017  mk  Written.
% 7-Aug-2018  mk  Change .MappedX/Y to derive from normalized touch surface
%                 coordinates and screen geometry, as that does cause less
%                 hassle on multi-x-screen setups.

if nargin < 1 || isempty(deviceIndex)
  error('Required deviceIndex missing.');
end

if nargin < 2 || isempty(windowHandle)
  error('Required windowHandle missing.');
end

if nargin < 3
  maxWaitTimeSecs = [];
end

% Go get it, go-getter!
[event, nremaining] = PsychHID('KbQueueGetEvent', deviceIndex, maxWaitTimeSecs);

if ~isempty(event)
  if event.Type == 5
    % Sequence fail event!
    fprintf('TouchEventGet: WARNING! Data loss in touch sequence detected!\n');
  else
    % Normal event, map also to window relative coordinates:
    [w, h] = Screen('Windowsize', Screen('WindowScreenNumber', windowHandle));
    winRect = Screen('GlobalRect', windowHandle);
    event.MappedX = (event.NormX * w) - winRect(1);
    event.MappedY = (event.NormY * h) - winRect(2);
  end
end

return;
