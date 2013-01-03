function PsychDebugWindowConfiguration(opaqueForHID, opacity)
% Switch PTB's onscreen window into a display mode suitable for easy debugging on single-screen setups.
%
% This function allows to setup Screen onscreen windows to be
% half-transparent, so one can simultaneously see the stimulus display and
% the Matlab window and other GUI windows.
%
% Usage: PsychDebugWindowConfiguration([opaqueForHID=0][, opacity=0.5])
%
% To enable: Call PsychDebugWindowConfiguration before the Screen('OpenWindow',...) call!
% To disable: Type "clear Screen"
%
% The optional parameter 'opaqueForHID' if set to a non-zero value will
% disallow mouse clicks and other mouse actions to "get through" to the
% GUI, ie., it will make the onscreen window opaque to the mouse pointer.
%
% A setting of -1 will disable the debug mode again.
%
% The optional parameter 'opacity' controls how opaque the onscreen window
% is, in a range of 0.0 to 1.0 for 0% to 100% opacity. By default the
% window is 50% opaque (or 50% transparent if you like).
%
% Stimulus onset timing and timestamping will be inaccurate in this mode
% and graphics performance will be reduced! Don't use for timing tests or
% during real experiment sessions!
%
% This feature will only work reliably - or at all - if your operating
% system is running with a compositing window manager installed and
% enabled. This is the case for Windows Vista, Windows-7 and later MS
% operating systems, MacOS/X, and GNU/Linux distributions that have the
% Compiz window manager installed and enabled, e.g., Ubuntu-9.1 and later.
%
% Keyboard and mouse input may not work as expected under all conditions,
% i.e., it may by impaired in either Psychtoolbox, or for the other running
% applications. Good luck!
%

% History:
% 30.07.2009 mk  Written.
% 15.11.2009 mk  Now also for Windows and Linux.

if nargin < 1
    opaqueForHID = [];
end

if isempty(opaqueForHID)
    opaqueForHID = 0;
end

if nargin < 2
    opacity=0.5;
end

if IsOSX
    % Disable high precision timestamping:
    Screen('Preference', 'VBLTimestampingMode', -1);
    
    % Use AGL + Carbon, even for fullscreen windows:
    oldconserve = Screen('Preference', 'ConserveVRAM');
    Screen('Preference', 'ConserveVRAM', bitor(oldconserve, 16384));
end

% Skip sync tests:
Screen('Preference', 'SkipSyncTests', 2);

% Map range 0.0 - 1.0 to 0 - 499:
opacity = floor(opacity * 499);
opacity = max(min(opacity, 499), 0);

if opaqueForHID
    % Set windows to be transparent, but not for mouse and keyboard:
    Screen('Preference', 'WindowShieldingLevel', 1500 + opacity);
else
    % Set windows to be transparent, also for mouse and keyboard:
    Screen('Preference', 'WindowShieldingLevel', 1000 + opacity);
end

if opaqueForHID == -1
    % Set windows to be normal, i.e., completely opaque:
    Screen('Preference', 'WindowShieldingLevel', 2000);
end

return;
