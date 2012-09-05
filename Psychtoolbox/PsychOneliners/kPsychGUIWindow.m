function rc = kPsychGUIWindow
% kPsychGUIWindow -- Create onscreen windows with behaviour of normal GUI windows.
%
% This flag can be passed to the optional 'specialFlags' parameter of
% Screen('OpenWindow', ...) or PsychImaging('OpenWindow', ...).
%
% It will cause the onscreen window to be a "regular" window that mostly
% behaves like typical GUI windows on your system. The window will have a
% titlebar and title, a border and other decorations. It will have buttons
% and handles to allow it to be moved around, resized, minimized or
% maximized, hidden and so on. Functions like Screen('Rect'),
% Screen('GlobalRect') and Screen('WindowSize') will report the true size
% and position of the window after it has been resized or moved around. The
% GetMouse() function will optionally report if the window has keyboard
% input focus because it is the active foreground window.
%
% Window stacking order, transparency and other window manager interactions
% should mostly behave as with other application windows.
%
% Please note that timing precision and timestamp precision for visual
% stimulus onset for this mode will not be guaranteed. Performance may be
% reduced. Other limitations may apply.
%
% GUI window mode is a "best effort" behaviour, as Psychtoolbox is not
% really designed to be - or exactly behave - like a regular GUI toolkit.
%

% This is the numeric constant for GUI window mode:
rc = 32;

return;
