function [secs, keyCode, deltaSecs] = KbWait(deviceNumber, forWhat, untilTime)
% [secs, keyCode, deltaSecs] = KbWait([deviceNumber][, forWhat=0][, untilTime=inf])
%
% Waits until any key is down and optionally returns the time in seconds
% and the keyCode vector of keyboard states, just as KbCheck would do. Also
% allows to wait for release of all keys or for single keystrokes, see
% below.
%
% If the optional parameter 'untilTime' is provided, KbWait will only wait
% until that time and then return regardless if anything happened on the
% keyboard or not.
%
% CAUTION: KbWait periodically checks the keyboard. After each failed check
% (ie. no change in keyboard state) it will wait for 5 msecs before the
% next check. This is done to reduce the load on your system, and it is
% important to do so. However if you want to measure reaction times this is
% clearly not what you want, as it adds up to 5 msecs extra uncertainty to
% all measurements!
%
% If you have trouble with KbWait always returning immediately, this could
% be due to "stuck keys". See "help DisableKeysForKbCheck" on how to work
% around this problem. See also "help RestrictKeysForKbCheck".
%
% GetChar and CharAvail are character oriented (and slow), whereas KbCheck
% and KbWait are keypress oriented (and fast).
%
% Using KbWait from the MATLAB command line: When you type "KbWait" at the
% prompt and hit the enter/return key to execute that command, then KbWait
% will detect the enter/return key press and return immediatly.  If you
% want to test KbWait from the command line, then try this:
%
%  WaitSecs(0.2);KbWait
%
% KbWait can also wait for releasing of keys instead of pressing of keys
% if you set the optional 2nd argument 'forWhat' to 1.
%
% If you want to wait for a single keystroke, set the 'forWhat' value to 2.
% KbWait will then first wait until all keys are released, then for the
% first keypress, then it will return. The above example could be realized
% via:
%
%  KbWait([], 2);
%
% If you would set 'forWhat' to 3 then it would wait for releasing the key
% after pressing it againg, ie. waitForAllKeysReleased -> waitForKeypress
% -> waitForAllKeysReleased -> Return [secs, keyCode] of the key press.
%
%
% OSX and Linux: __________________________________________________________
%
% KbWait uses the PsychHID function, a general purpose function for
% reading from the Human Interface Device (HID) class of USB devices.
%
% KbWait tests the first USB-HID keyboard device by default. Optionally
% you can pass in a 'deviceNumber' to test a different keyboard if multiple
% keyboards are connected to your machine.  If deviceNumber is -1, all
% keyboard devices will be checked.  If deviceNumber is -2, all keypad
% devices (if any) will be checked. If deviceNumber is -3, all keyboard and
% keypad devices will be checked. The device numbers to be checked are
% determined only on the first call to the function.  If these numbers
% change, the function can be reset using "clear KbWait".
%
% As a little bonus, KbWait can also query other HID human input devices
% which have keys or buttons as if they were keyboards. If you pass in the
% deviceIndex of a mouse (GetMouseIndices will provide with them), it will
% treat mouse button state as keyboard state. Similar behaviour usually
% works with Joysticks, Gamepads and other input controllers.
%
% _________________________________________________________________________
%
% See also: KbCheck, KbStrokeWait, KbPressWait, KbReleaseWait, GetChar, CharAvail, KbDemo.

% 3/6/97    dhb  Wrote it.
% 8/2/97    dgp  Explain difference between key and character. See KbCheck.
% 9/06/03   awi  ****** OS X-specific fork from the OS 9 version *******
%                  Added OS X conditional.
% 7/12/04   awi  Cosmetic.  OS 9 Section. Uses IsOSX.
% 4/11/05   awi  Added to help note about testing kbWait from command line.
% 11/29/05  mk   Fixed really stupid bug: deviceNumber wasn't queried!
% 02/22/06  mk   Modified for Linux: Currently a hack.
% 10/24/06  mk   Replaced by a generic implementation that just uses KbCheck
%                in a while loop. This way we directly benefit from KbChecks
%                improvements.
% 3/15/07   kas  Added in option to poll all keyboard devices by passing
%                deviceNumber == -1
%
% 3/03/08   mk   Added option 'forWhat' to optionally wait for key release
%                or isolated keystrokes, and optional return argument 'keyCode'
%                to return keyCode vector, just as KbCheck does.
%
% 12/27/09  mk   Remove all the redundant code for 'deviceNumber' specific
%                behaviour. This is already covered by code in KbCheck!
%                This also fixes a bug reported in forum message 10468
%                where KbReleaseWait(-1) didn't wait for all keys on all
%                keyboards to be released.
% 12/18/09  rpw  Added documentation about keypad devices on OS/X.

% Time (in seconds) to wait between "failed" checks, in order to not
% overload the system in realtime mode. 5 msecs seems to be an ok value...
yieldInterval = 0.005;

if nargin < 2
    forWhat = 0;
end

if nargin == 0
    deviceNumber = [];
end

if nargin < 3
    untilTime = inf;
end

if isempty(untilTime)
    untilTime = inf;
end

% Wait for keystroke?
if (forWhat == 2) || (forWhat == 3)
    % Wait for keystroke, ie., first make sure all keys are released, then
    % wait for a keypress:
    
    % Wait for key release. we know we have deviceNumber valid here:
    KbWait(deviceNumber, 1, untilTime);
    
    if forWhat == 2
        % Now just go on with forWhat = 0, ie., wait for keypress:
        forWhat = 0;
    else
        % Wait for keypress:
        [secs, keyCode, deltaSecs] = KbWait(deviceNumber, 0, untilTime);
        
        % Wait for key release. we know we have deviceNumber valid here:
        KbWait(deviceNumber, 1, untilTime);

        return;
    end
end

secs = -inf;
while secs < untilTime
    [isDown, secs, keyCode, deltaSecs] = KbCheck(deviceNumber);
    if (isDown == ~forWhat) || (secs >= untilTime)
        return;
    end

    % Wait for yieldInterval to prevent system overload.
    secs = WaitSecs('YieldSecs', yieldInterval);
end
