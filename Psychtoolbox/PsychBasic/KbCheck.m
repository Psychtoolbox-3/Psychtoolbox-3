function [keyIsDown,secs, keyCode, deltaSecs] = KbCheck(deviceNumber, unusedUntilTime, varargin)
% [keyIsDown, secs, keyCode, deltaSecs] = KbCheck([deviceNumber])
% 
% Return keyboard status (keyIsDown), time (secs) of the status check, and
% keyboard scan code (keyCode).
% 
%    keyIsDown      1 if any key, including modifiers such as <shift>,
%                   <control> or <caps lock> is down. 0 otherwise.
% 
%    secs           Time of keypress as returned by GetSecs.
% 
%    keyCode        A 256-element logical array.  Each bit
%                   within the logical array represents one keyboard key. 
%                   If a key is pressed, its bit is set, othewise the bit 
%                   is clear. To convert a keyCode to a vector of key  
%                   numbers use FIND(keyCode). To find a key's keyNumber 
%                   use KbName or KbDemo.
% 
%    deltaSecs      Time in seconds since this KbCheck query and the most
%                   recent previous query (if any). This value is in some
%                   sense a confidence interval, e.g., for reaction time
%                   measurements. If KbCheck returns the information that a
%                   key is pressed by the subject, then the subject could
%                   have pressed the key down anytime between this
%                   invocation of KbCheck at time 'secs' and the most
%                   recent previous invocation. Therefore, 'deltaSecs'
%                   tells you about the interval in which depression of the
%                   key(s) might have happened: [secs - deltaSecs; secs].
%                   for practical purpose this means that "measured" RT's
%                   can't be more accurate than 'deltaSecs' seconds - the
%                   interval between the two most recent keyboard checks.
%                   Please note however, that standard computer keyboards
%                   can incur additional delays and timing uncertainty of
%                   up to 50 msecs, so the real uncertainty can be higher
%                   than 'deltaSecs' -- 'deltaSecs' is just a lower bound!
%
% KbCheck and KbWait determine whether any key is down now, including the
% meta keys: <caps lock>, <shift>, <command>, <control>, and <option>. The
% only key not reported is the start key (triangle) used to power on your
% computer.
% 
% Some users of Laptops experienced the problem of "stuck keys": Some keys
% are always reported as "down", so KbWait returns immediately and KbCheck
% always reports keyIsDown == 1. This is often due to special function keys.
% These keys or system functionality are assigned vendor specific
% key codes, e.g., the status of the Laptop lid (opened/closed) could be
% reported by some special keycode. Whenever the Laptop lid is open, this key
% will be reported as pressed. You can work around this problem by passing
% a list of keycodes to be ignored by KbCheck and KbWait. See
% "help DisableKeysForKbCheck" on how to do this.
%
% Keys pressed by the subject often show up in the Matlab command window as
% well, cluttering that window with useless character junk. You can prevent
% this from happening by disabling keyboard input to Matlab: Add a
% ListenChar(2); command at the beginning of your script and a
% ListenChar(0); to the end of your script to enable/disable transmission of
% keypresses to Matlab. If your script should abort and your keyboard is
% dead, press CTRL+C to reenable keyboard input -- It is the same as
% ListenChar(0). See 'help ListenChar' for more info.
%
% GetChar and CharAvail are character-oriented (and slow), whereas KbCheck
% and KbWait are keypress-oriented (and fast). If only a meta key was hit,
% KbCheck will return true, because a key was pressed, but CharAvail will
% return false, because no character was generated. See GetChar.
% 
% KbCheck and KbWait are MEX files, which take time to load when they're
% first called. They'll then stay loaded until you flush them (e.g. by
% changing directory or calling CLEAR MEX).
%
% OSX, Linux, and Windows with Octave or Matlab R2007a and later: _________
%
% KbCheck uses the PsychHID function, a general purpose function for
% reading from the Human Interface Device (HID) class of USB devices.
%
% KbCheck queries the first USB-HID keyboard device by default. Optionally,
% when multiple keyboards are attached to your machine, you can pass in a
% 'deviceNumber':  When 'deviceNumber' is -1, KbCheck will query all
% keyboard devices and return their "merged state" - The 'keyCode' vector
% will represent the state of all keys of all keyboards, and the
% 'keyIsDown' flag will be equal to one if at least one key on any of the
% keyboards is pressed. When 'deviceNumber' is -2, KbCheck will query all
% keypad devices (if any) and return their "merged state", and when
% 'deviceNumber' is -3, KbCheck will query all keyboard and keypad devices
% and return their "merged state". When 'deviceNumber' is greater than 0, it
% will query only the specified HID keyboard device corresponding to that
% 'deviceNumber'. The function GetKeyboardIndices() allows to query the
% device numbers of all attached keyboards, or keyboards matching specific
% criteria, and the function GetKeypadIndices() allows the same for keypads.
%
% On MS-Windows 2000 and earlier, KbCheck can address individual keyboards.
% On Windows-XP and later, it can't.
% 
% As a little bonus, KbCheck can also query other HID human input devices
% which have keys or buttons as if they were keyboards. If you pass in the
% deviceIndex of a mouse (GetMouseIndices will provide with them), it will
% report mouse button state as keyboard state. Similar behaviour usually
% works with Joysticks, Gamepads and other input controllers.
% _________________________________________________________________________
% 
% See also: FlushEvents, KbName, KbDemo, KbWait, GetChar, CharAvail.

% TO DO:
%
%  - Mention that on USB systems there the USB bus is sampled at 100 Hz.
%  - We could augment this to to accept an optional keyboard device number. 

% 3/6/97  dhb  Wrote it.
% 8/2/97  dgp  Explain difference between key and character.
% 1/28/98 dgp  Explain CapsLock.
% 2/4/98  awi  Explain keyCode.
% 2/13/98 awi  Changed keyCode to logical array, pointers to KbDecode, KbExplore.
% 2/19/98 dgp  Shortened by moving some text to GetChar.m.
% 3/15/99 xmz  Added comment for Windows version.
% 6/23/00 awi  Added paragraph contrasting queuing of GetChar and KbCheck.
% 7/7/00  dgp  Cosmetic.
% 6/17/02 awi  ****** OS X-specific fork from the OS 9 version *******
%                Added conditional invocation of PsychHID on OSX
% 7/12/04 awi  Cosmetic.  Separted platform-specific help. Use IsOSX now. 
% 10/4/05 awi Note here cosmetic changes by dbp on unknown date between 7/12/04 and 10/4/05.  
% 10/24/06 mk Windows and Linux implementation: Use built-in helper code in Screen.
% 10/24/06 mk Add code for disabling "stuck keys".
% 6/13/08 abl Option for OS X to poll all keyboard devices by passing deviceNumber == -1, \
%             based on kas's modification of KbWait.
% 11/16/8  mk Allow to use "white-list" of enabled keys in
%             ptb_kbcheck_enabledKeys. This is set via
%             RestrictKeysForKbCheck, and passed to PsychHID in order to
%             restrict keyboard scans to a subset of enabled keys. This
%             provides a significant speedup for KbChecks if used properly.
%             On non-OS/X, this is emulated in software and does not
%             provide any speedup.
% 03/03/9  mk Bugfix for "white-list" code on old Matlab releases. Need to
%             cast to double and back to uint8, as old Matlabs don't
%             support .* operation on uint8 class arrays.
% 12/18/09 rpw Added support for polling keypads on OSX via deviceNumber of -2 or -3
% 01/07/10 mk Code refactoring: Unified check-code for all deviceIndex
%             values.

% ptb_kbcheck_disabledKeys is a vector of keyboard scancodes. It allows
% to define keys which should never be reported as 'down', i.e. disabled
% keys. The vector is empty by default. If you have special quirky hardware,
% e.g., some Laptop keyboards, that reports some keys as 'always down', you
% can work around this 'stuck keys' by defining them in the ptb_kbcheck_disabledKeys
% vector.
global ptb_kbcheck_disabledKeys;

% ptb_kbcheck_enabledKeys is a white-list of keys to query. If set to
% non-empty 
global ptb_kbcheck_enabledKeys;

% Store timestamp of previous KbCheck:
persistent oldSecs;

% Cache operating system type to speed up the code below:
persistent macosx;
% ...and all keyboard indices as well:
persistent kbs kps;
persistent keyboardsdetected;

if isempty(macosx)
    % First time invocation: Query and cache type of OS:
    macosx = IsOSX;
    
    % Set initial oldSecs to minus infinity: No such query before...
    oldSecs = -inf;
    
    % Init ptb_kbcheck_enabledKeys to empty, if it hasn't been set
    % externally already:
    if ~exist('ptb_kbcheck_enabledKeys', 'var')
        ptb_kbcheck_enabledKeys = [];
    end
end

if nargin < 1
    deviceNumber = [];
end

if ~IsWin || (IsWin && ~isempty(deviceNumber))
    if ~isempty(deviceNumber)
        % All attached keyboards already detected?
        if isempty(keyboardsdetected)
            % No. Do it now:
            % Query indices of all attached keyboards, in case we need'em:
            kbs=GetKeyboardIndices;
            kps=GetKeypadIndices;
            keyboardsdetected = 1;
        end
        
        % Select keyboard(s):
        if deviceNumber==-1
            % Query all attached keyboards:
            keyt = kbs;
        elseif deviceNumber==-2
            % Query all attached keypads:
            keyt = kps; 
        elseif deviceNumber==-3
            % Query all attached keyboards and keypads:
            keyt = [kbs kps]; 
        else
            % Query a specific keyboard device number:
            keyt = deviceNumber;
        end

        if ~isempty(keyt)
            % Check all devices in vector keyt and merge their state:
            keyIsDown=0; keyCode=zeros(1,256);  % preallocate these variables
            for i=keyt
                [DeviceKeyIsDown, secs, DeviceKeyCode]= PsychHID('KbCheck', i, ptb_kbcheck_enabledKeys);
                keyIsDown = keyIsDown || DeviceKeyIsDown;
                keyCode = keyCode | DeviceKeyCode;
            end
        else
            [keyIsDown, secs, keyCode]= PsychHID('KbCheck', [], ptb_kbcheck_enabledKeys);
        end   
    else
        % Query primary keyboard:
        [keyIsDown, secs, keyCode]= PsychHID('KbCheck', [], ptb_kbcheck_enabledKeys);
    end
else
   % We use the built-in KbCheck facility of Screen on MS-Windows
   % for KbChecks if the usercode didn't specify any 'deviceIndex', so
   % the user gets a good "works out of the box" experience for the default
   % use case without the need to install the libusb-1.0.dll in the system,
   % which would be required for PsychHID on Windows to work.
   [keyIsDown,secs, keyCode]= Screen('GetMouseHelper', -1);
end

% Compute time delta since previous keyboard query, and update internal
% cached value:
deltaSecs = secs - oldSecs;
oldSecs = secs;

% Only need to apply ptb_kbcheck_enabledKeys manually on non-OS/X systems,
% as this is done internally in PsychHID('KbCheck') on OS/X:
if ~macosx && ~isempty(ptb_kbcheck_enabledKeys)
    % Mask all keys with the enabled keys:
    keyCode = uint8(double(keyCode) .* ptb_kbcheck_enabledKeys);

    % Reevaluate global key down state:
    keyIsDown = any(keyCode);
end

% Any dead keys defined?
if ~isempty(ptb_kbcheck_disabledKeys)
   % Yes. Disable all dead keys - force them to 'not pressed':
   keyCode(ptb_kbcheck_disabledKeys)=0;
   % Reevaluate global key down state:
   keyIsDown = any(keyCode);
end

% Must be double format for some client routines:
keyCode = double(keyCode);
return;
