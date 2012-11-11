function [ch, when] = GetChar(getExtendedData, getRawCode)
% [ch, when] = GetChar([getExtendedData], [getRawCode])
% 
% Wait for a typed character and return it.  If a character was typed
% before calling GetChar then GetChar will return that character immediatly.
% Characters flushed by FlushEvents are all ignored by GetChar. Characters
% are returned in the first return argument "ch".
%
% Please read the 'help ListenChar' carefully to understand various
% limitations and caveats of this function, and to learn about - often
% better - alternatives.
% 
% CAUTION: Do not rely on the keypress timestamps returned by GetChar
% without fully reading and understanding this help text. Run your own
% timing tests on GetChar and KbCheck to verify that the timing is good
% enough and avoid GetChar for timed keypresses if possible at all. Use
% KbWait and KbCheck instead.
%
% The main purpose of GetChar is to reliably collect keyboard input in the
% background while your experiment script is occupied with performing other
% operations, e.g., Matlab computations, sound output or visual stimulus
% drawing. After an initial call to ListenChar, the operating system will
% record all keyboard input into an internal queue. GetChar removes
% characters from that queue, one character per invocation of GetChar. You
% can empty that queue any time by calling FlushEvents('keyDown').
%
% If you want to check the current state of the keyboard, e.g., for
% triggering immediate actions in response to a key press, waiting for a
% subjects response, synchronizing to keytriggers (e.g., fMRI machines) or
% if you require high timing precision then use KbCheck instead of GetChar.
%
% GetChar should work on all platforms, but its specific functionality,
% beyond simply returning typed characters, will vary depending on OS type
% and version, if you use Matlab or Octave, and if you use Matlab with or
% without Java based GUI active. For portability it is therefore best to
% ignore all info returned beyond the character code.
%
% "when" is a struct. It (used to) return the time of the keypress, the "adb"
% address of the input device, and the state of all the modifier keys
% (shift, control, command, option, alphaLock) and the mouse button.
% "when.secs" is an estimate, of what GetSecs would have been. Since it's
% derived sometime from a timebase different from the timebase of GetSecs,
% times returned by GetSecs are not directly comparable to when.secs.
%
% By setting getExtendedData to 0, all extended timing/modifier information
% will not be collected and "when" will be returned empty.  This speeds up
% calls to this function. If ommitted or set to 1, the "when" data
% structure is filled.  getRawCode set to 1 will set "ch" to be the integer
% ACII code of the available character.  If ommitted or set to 0, "ch" will
% be in char format. When running under Linux in "matlab -nojvm" mode or on
% Octave, "when" will be returned empty. When running on any other
% operating system under Octave or in "matlab -nojvm" mode, or on Windows
% Vista and later versions of the Windows OS, when will only contain a
% valid timestamp, but all other fields will be meaningless.
%
% GetChar and CharAvail are character-oriented (and slow), whereas KbCheck
% and KbWait are keypress-oriented (and fast). If only a meta key (like
% <option> or <shift>) was hit, KbCheck will return true, because a key was
% pressed, but CharAvail will return false, because no character was
% generated. See KbCheck.
% 
% CharAvail and GetChar use the system event queue to retrieve the character
% generated, not the raw key press(es) per se. If the user presses "a",
% GetChar returns 'a', but if the user presses option-e followed by "a",
% this selects an accented a, "?", which is treated by GetChar as a single
% character, even though it took the user three keypresses (counting the
% option key) to produce it.
% 
% There can be some delay between when the key is pressed and when CharAvail
% or GetChar detects it, due to internal processing overhead in Matlabs Java
% implementation. GetChar internally collects timestamps in the timebase
% used by Matlabs Java implementation, whereas other Psychtoolbox timing functions
% (GetSecs, Screen('Flip'), KbCheck, KbWait, ...) use time reported by some
% high precision system timer. The "when.secs" time reported by GetChar is
% converted from Java timebase to Psychtoolboxs timebase. Due to conversion
% errors mostly out of our control, the reported values can be off by
% multiple dozen or even hundreds of milliseconds from what KbWait, KbCheck
% or GetSecs would report. Example: A high-end Pentium-4 3.2 Ghz system
% running Windows-XP has been measured to be off by 40 to 70 milliseconds.
%
% Some Java implementations are also known to have problems/bugs in
% timestamping keyboard presses properly and each Matlab version on each
% operating system is bundled with a different Java version, so some Matlab
% versions may be reliable with respect to GetChars timing, whereas others
% are not.
%
% ---> If precise timing of the keypress is important, use KbCheck or
% KbWait or KbQueueXXX functions or KbEventGet for consistent results!
%
% OS X / Windows-XP / Linux with Matlab and Java enabled: _________________
%
% JAVA PATH: The GetChar implementation for Matlab is based on Java.
% Therefore, the Psychtoolbox subfolder PsychJava must be added to Matlabs
% static classpath. Normally this is done by the Psychtoolbox installer by
% editing the Matlab file "classpath.txt" (enter which('classpath.txt') to
% find the location of that file). If the installer fails to edit the file
% properly, you'll need to perform that step manually by following the
% instructions of the installer. See 'help PsychJavaTrouble' for more infos
% on this.
%
% KEYSTROKES IN THE BACKGROUND: To detect keypresses made before the
% GetChar call, you must have called "ListenChar" earlier.  ListenChar
% redirects keystrokes to the GetChar queue. Calling ListenChar at the
% beginning of your script should cause GetChar to behave identically
% to OS 9 GetChar with respect to background keystroke collection.
%
% KEYSTROKES IN MATLAB WINDOW: By default, all keystrokes are also sent to
% Matlabs window, generating some ugly clutter. You can suppress this by
% calling ListenChar(2), so your Matlab console stays nice and clean. Don't
% forget to call ListenChar(1) or ListenChar(0) though before the end of
% your script. If Matlab returns to its command prompt without reenabling
% keyboard input via ListenChar(0) or ListenChar(1), Matlab will be left
% with a dead keyboard until you press the CTRL+C key combo. This silencing
% of clutter does currently not work in matlab -nojvm mode, or if you use
% GNU/Octave instead of Matlab.
%
% OTHER "when" RETURN ARGUMENT FIELDS: Owing to differences in what
% accessory information the underlying operating systems provides about
% keystrokes, "when' return argument fields differs between operating systems.
% GetChar sets fields for which it returns no value to the value 'Nan'.  
%
% See also: ListenChar, CharAvail, FlushEvents, TestGetChar, KbCheck,
% KbWait

% 5/7/96  dgp	Wrote this help file.
% 1/22/97 dhb	Added comment and pointer to TIMER routines.
% 3/6/97  dhb	References to KbWait, KbCheck.
% 7/23/97 dgp	It's a character not a keypress.
% 8/2/97  dgp	Explain difference between key and character. See KbCheck.
% 2/7/98  dgp	Streamlined. Eliminated call to GetKey, since it's now GetChar.mex.
% 3/24/98 dgp	Explain backgrounding and meta keys. Don't mention obsolete GetKey and KbHit.
% 3/15/99 xmz	Put in some comment for Windows version.
% 3/19/99 dgp	Update explanation of backgrounding. 
% 3/28/99 dgp	Show how to turn off backgrounding. 
% 8/19/00 dgp	Cosmetic. 
% 4/23/02 dgp   Fixed GetChar.mex to always quit on command-period.
% 4/27/02 dgp	Added optional second return argument.
% 6/1/02  dgp	Mention Tick0Secs.
% 9/21/02 dgp   Added address field to "when", using C code donated by Tom Busey.
% 7/12/04 awi   ****** OS X-specific fork from the OS 9 version *******
%                   Expanded on explantion in the first paragraph about when
%                   GetChar waits, when it returns immediatly, what
%                   it ignores.  Retains OS 9-specific comments.    
% 1/27/04 awi   Issue an error when calling GetChar and suggest KbWait. 
% 1/19/05 awi   Implemented GetChar on OS X.  Added AssertMex for OS 9 and OS X conditional block.
% 7/20/05 awi   Wrote OS X documentation section.
% 2/2/06  awi   Tested to see if this works when the MATLAB text editing
%               window is minimized. It does not.
% 2/22/06 awi  Commented out Cocoa wrapper and wrote Java wrapper.
% 3/28/06 awi  Detect buffer overflow.
%              Handle new double value from .getChar(), was char type.
%              Changed "char" return value to "ch" to avoid name conflict with
%               built-in MATLAB function "char" 
% 6/15/06 awi  Added a second return argument.
%              Updated built-in help for the Java implementation.
% 6/15/06 awi  Added break on CTRL-C
%              Added TO DO section and item to detect genuine KeyDown
%              events.
% 6/20/06 awi  Use AddPsychJavaPath instead of AssertGetCharJava.
% 8/16/06 cgb  Now using the new GetChar system which taps straight into
%              the java keypress dispatcher.
% 9/18/06  mk  GetChar now works on all Matlabs (OS-X, Windows) in JVM
%              mode. In -nojvm mode on Windows, it falls back to the old
%              Windows DLL ...
%
% 05/31/09 mk  Add support for Octave and Matlab in noJVM mode.
% 10/22/12 mk  Remove support for legacy Matlab R11 GetCharNoJVM.dll.
% 10/22/12 mk  Add support for KbQueue-Based implementation.

% NOTES:
%
% The second return argument from OS 9 PTB looks like this:
%     ticks: 5760808
%     secs: 1.4681e+05
%     address: 2
%     mouseButton: 0
%     alphaLock: 0
%     commandKey: 0
%     controlKey: 0
%     optionKey: 0
%     shiftKey: 0
% 

global OSX_JAVA_GETCHAR;

% If no command line argument was passed we'll assume that the user only
% wants to get character data and timing/modifier data.
if nargin == 0
    getExtendedData = 1;
    getRawCode = 0;
elseif nargin == 1
    getRawCode = 0;
end

% Is this Matlab? Is the JVM running? Isn't this Windows Vista or later?
if psychusejava('desktop') && ~IsWinVista
    % Java virtual machine and AWT and Desktop are running. Use our Java based
    % GetChar.

    % Make sure that the GetCharJava class is loaded and registered with
    % the java focus manager.
    if isempty(OSX_JAVA_GETCHAR)
        try
            OSX_JAVA_GETCHAR = AssignGetCharJava;
        catch %#ok<*CTCH>
            error('Could not load Java class GetCharJava! Read ''help PsychJavaTrouble'' for help.');
        end
        OSX_JAVA_GETCHAR.register;
    end

    % Make sure the Matlab window has keyboard focus:
    if exist('commandwindow') %#ok<EXIST>
        % Call builtin implementation:
        commandwindow;
    end

    % Loop until we receive character input.
    keepChecking = 1;
    while keepChecking
        % Check to see if a character is available, and stop looking if
        % we've found one.
        charValue = OSX_JAVA_GETCHAR.getChar;
        keepChecking = charValue == 0;
        if keepChecking
            drawnow;
        end
    end

    % Throw an error if we've exceeded the buffer size.
    if charValue == -1
        % Reenable keystroke dispatching to Matlab to leave us with a
        % functional Matlab console.
        OSX_JAVA_GETCHAR.setRedispatchFlag(0);
        error('GetChar buffer overflow. Use "FlushEvents" to clear error');
    end

    % Get the typed character.
    if getRawCode
        ch = charValue;
    else
        ch = char(charValue);
    end

    % Only fill up the 'when' data stucture if extended data was requested.
    if getExtendedData
        when.address=nan;
        when.mouseButton=nan;
        when.alphaLock=nan;
        modifiers = OSX_JAVA_GETCHAR.getModifiers;
        when.commandKey = modifiers(1);
        when.controlKey = modifiers(2);
        when.optionKey = modifiers(3);
        when.shiftKey = modifiers(4);
        rawEventTimeMs = OSX_JAVA_GETCHAR.getEventTime;  % result is in units of ms.
        when.ticks = nan;
        when.secs = JavaTimeToGetSecs(rawEventTimeMs, -1);
    else
        when = [];
    end

    return;
end

% Running either on Octave, or on Matlab in No JVM mode or on MS-Vista+:

% If we are on Linux and the keyboard queue is already in use by usercode,
% we can fall back to 'GetMouseHelper' low-level terminal tty magic. The
% only downside is that typed characters will spill into the console, ie.,
% ListenChar(2) suppression is unsupported:
if IsLinux && KbQueueReserve(3, 2, [])
    % Loop until we receive character input.
    keepChecking = 1;
    while keepChecking
        % Check to see if a character is available, and stop looking if
        % we've found one.
        
        % KeyboardHelper with command code 15 delivers
        % id of currently pending characters on stdin:
        charValue = PsychHID('KeyboardHelper', -15);
        keepChecking = charValue == 0;
        if keepChecking
            drawnow;
            if exist('fflush') %#ok<EXIST>
                builtin('fflush', 1);
            end
        end
    end
    
    % Throw an error if we've exceeded the buffer size.
    if charValue == -1
        % Reenable keystroke display to leave us with a
        % functional console.
        PsychHID('KeyboardHelper', -11);
        error('GetChar buffer overflow. Use "FlushEvents" to clear error');
    end

    % No extended data in this mode:
    when = [];
else
    % Use keyboard queue:
    
    % Only need to reserve/create/start queue if we don't have it
    % already:
    if ~KbQueueReserve(3, 1, [])
        % LoadPsychHID is needed on MS-Windows. It no-ops if called redundantly:
        LoadPsychHID;
        
        % Try to reserve default keyboard queue for our exclusive use:
        if ~KbQueueReserve(1, 1, [])
            % Ok, we have to abort here. While the same issue is only worth
            % a warning for CharAvail and FlushEvents, it is game over for
            % GetChar, as we must not touch a user managed kbqueue, and we
            % can't provide any sensible behaviour if we can't do that:
            error('Keyboard queue for default keyboard device already in use by KbQueue/KbEvent functions et al. Use of ListenChar/GetChar etc. and keyboard queues is mutually exclusive!');
        end
        
        % Got it. Allocate and start it:
        PsychHID('KbQueueCreate');
        PsychHID('KbQueueStart');

        if (IsOSX(1) || (IsOctave && IsGUI))
            % Enable keystroke redirection via kbyqueue and pty to bypass
            % blockade of onscreen windows:
            PsychHID('KeyboardHelper', -14);
        end        
    end
    
    % Queue is running: Poll it for new events we're interested in:
    % Loop until we receive character input.
    keepChecking = 1;
    while keepChecking
        % Check to see if a character is available, and stop looking if
        % we've found one.        
        event = PsychHID('KbQueueGetEvent', [], 0.1);
        if ~isempty(event) && event.Pressed && (event.CookedKey > 0)
            charValue = event.CookedKey;
            keepChecking = 0;
        else
            charValue = 0;
            keepChecking = 1;
        end
        
        if keepChecking
            drawnow;
            if exist('fflush') %#ok<EXIST>
                builtin('fflush', 1);
            end
        end
    end
    
    % Only fill up the 'when' data stucture if extended data was requested.
    if getExtendedData
        when.address=nan;
        when.mouseButton=nan;
        when.alphaLock=nan;
        modifiers = [nan, nan, nan, nan];
        when.commandKey = modifiers(1);
        when.controlKey = modifiers(2);
        when.optionKey = modifiers(3);
        when.shiftKey = modifiers(4);
        when.ticks = nan;
        when.secs = event.Time;
    else
        when = [];
    end
end

% Get the typed character.
if getRawCode
    ch = charValue;
else
    ch = char(charValue);
end

return;
