function ListenChar(listenFlag)
% function ListenChar([listenFlag])
%
% Tell the Psychtoolbox function "GetChar" to start or stop listening for
% keyboard input. By default ListenChar enables listening when no argument is
% supplied. Passing 0 will turn off character listening and reset the
% buffer which holds the captured characters. Passing a value of 1 or not
% passing any value will enable listening. Passing a value of 2 will enable
% listening, additionally any output of keypresses to Matlabs or Octaves
% windows is suppressed. Use this with care, if your script aborts with an
% error, Matlab or Octave may be left with a dead keyboard until the user
% presses CTRL+C to re-enable keyboard input. 'listenFlag' 2 is silently
% ignored on Matlab in -nojvm mode under MS-Windows.
%
% Passing the listenFlag -1 will only suppress keypresses into Matlabs or
% Octaves command window, but not collect characters for use with CharAvail
% or GetChar. This allows concurrent use of ListenChar for character suppression
% while at the same time using keyboard queues for character and key input.
% Reenabling keypresses into Matlab and Octave is only possible via a call
% with listenFlag 0 if suppression was enabled with listenFlag -1. Use of
% listenFlag -1 does not block keyboard queues, so if you want to get keyboard
% input, the better way might be to use keyboard queues (see KbEventAvail,
% KbEventGet and KbEventFlush) for input and ListenChar(-1) and ListenChar(0)
% to suppress characters into the Matlab/Octave command window or console.
%
% This function isn't entirely necessary to turn on listening as calling
% GetChar, CharAvail, or FlushEvents will trigger listening on. However,
% it is the only method by which to disable listening or switch between
% suppression of keyboard input to Matlab and unsuppressed mode.
%
% Please note that the commands ListenChar, CharAvail and GetChar are
% subject to various constraints and limitations, depending on the
% operating system you use, if you use a Matlab with/without Java based GUI,
% or if you use Octave, if you have Screen() onscreen windows open or not, or
% if you use KbQueueXXX functions in parallel or not. Therefore use of these
% functions can be troublesome for any but the most simple usages. Use of
% KbCheck, KbWait, KbStroke/Press/ReleaseWait is often simpler if you are
% just after keyboard input. Use of KbQueue functions, e.g., KbQueueCheck,
% KbQueueWait, KbTriggerWait is better suited for background keystroke
% collection. Use of KbEventAvail and KbEventGet is often more convenient,
% more flexible and subject to less restrictions and gotchas than use of
% GetChar et al.
%
% Some of the restrictions and caveats:
%
% 1. Works very well with Matlab releases older than R2025a and their Java
% based GUIs enabled on Linux (but only on desktop GUI's other than KDE),
% and on macOS.
%
% 2. When used on Microsoft Windows with old Matlab versions with Java GUI,
% you cannot use any KbQueue functions at the same time, ie.,
% KbQueueCreate/Start/Stop/Check/Wait as well as KbWaitTrigger,
% KbEventFlush, KbEventAvail, and KbEventGet are off limits after any call
% to ListenChar, ListenChar(1), ListenChar(2), FlushEvents, CharAvail or
% GetChar. You would need to call ListenChar(0) before you could call
% KbQueueCreate and then use those functions. Vice versa, after a call to
% KbQueueCreate, CharAvail, FlushEvents, and GetChar are dysfunctional and
% ListenChar may be limited. You need to call KbQueueRelease before you can
% use them again. Use of other devices, e.g., mouse or joystick, is not
% prohibited during use of GetChar et al.
%
% 3. If you use Matlab in "matlab -nojvm" or "matlab -nodesktop" mode without
% its GUI, or if you use a Matlab of version R2025a or later, or you use
% GNU/Octave instead of Matlab, the same restrictions as in 2. apply -
% no parallel use of the default keyboards KbQueue. KbQueues can be used
% for other input devices and on Linux and macOS for keyboards other than the
% default keyboard.
%
% The only feature that works in parallel with KbQueues on the default keyboard
% is the suppression of spilling of keystroke characters into the Matlab or
% Octave window during ListenChar(2) - at least on Linux and macOS with old
% Matlab versions before R2025a and their Java based GUIs. On Windows this
% can't be prevented at all in "matlab -nojvm" mode. However, if you switch
% to ListenChar(2) mode, you cannot break out of it by pressing CTRL+C on
% Linux if the keyboard queue that is in parallel use didn't get
% KbQueueStart called, ie. if it is stopped. On macOS with a stopped
% Keyboard queue, neither CTRL+C nor stopping a runaway script works.
%
% 4. On Linux, as a exception, some GetChar, CharAvail functionality may
% still work in case 3. under certain conditions, e.g., if you don't use
% ListenChar(2) and your Matlab/Octave is not in interactive mode.
%
% Also, GetChar can only collect keystrokes from multiple connected
% keyboards in case 1. In all other cases, it can only collect keystrokes,
% or respond to press of CTRL+C, for the default keyboard device. It will
% ignore other connected keyboards.
%
% Another limitation in cases 2 and 3 is that international keyboards with non-US
% layout and non-ASCII/Latin-1 characters may need special treatment and Octave
% may have some trouble processing such characters. See "help KbEventGet" for a
% detailed explanation.
%
% Basically: Mixing GetChar et al. and modern KbQueue functions is usually
% not advisable, or if needed, great care must be taken to sidestep all the
% mentioned limitations. Also the KbQueue functions usually have better
% timing precision and allow to flexibly address multiple keyboards
% separately at least on Linux and macOS.
%
%
% For further explanation see help for "GetChar".  
%
% _________________________________________________________________________
%
% See also: GetChar

% HISTORY
%
% 7/19/05  awi   Wrote it.
% 6/20/06  awi   Use AddPsychJavaPath instead of AssertGetCharJava.
% 8/31/06  cgb   Works with the new character listening system.
% 9/19/06  mk    Modified to work on all Java enabled Matlabs and be a no-op
%                in all other configurations.
% 10/13/06 mk    Support for setting the redispatch-mode of GetChar and
%                friends.
% 05/31/09 mk    Add support for Octave and Matlab in noJVM mode.
% 01/31/16 mk    Add support for listenFlag -1 for only blocking input.
% 06/20/19 mk    Try to protect against KDE focus stealing nastiness via kbqueues.
% 06/22/25 mk    Make baseline compatible with Matlab R2025a+ non-Java GUI.

global OSX_JAVA_GETCHAR; %#ok<GVMIS>
persistent keyboard_blocked;
persistent isjavadesktop;
 
if nargin == 0
    listenFlag = 1;
elseif nargin > 1
    error('Too many arguments to ListenChar!  See "help ListenChar" for more information');
end

if ~ismember(listenFlag, [-1,0,1,2])
    error('Invalid listenFlag provided!  See "help ListenChar" for more information');
end

if isempty(keyboard_blocked)
    keyboard_blocked = 0;
end

% Is this old Matlab with Java based GUI? Only check this once because psychusejava is a slow command.
if isempty(isjavadesktop)
    isjavadesktop = psychusejava('desktop');
end

if isjavadesktop
    % Java GUI enabled on Matlab. There's work to do.

    % Make sure that the GetCharJava class is loaded.
    if isempty(OSX_JAVA_GETCHAR)
        try
            OSX_JAVA_GETCHAR = AssignGetCharJava;
        catch %#ok<*CTCH>
            error('Could not load Java class GetCharJava! Read ''help PsychJavaTrouble'' for help.');
        end
    end

    if listenFlag ~= 0
        % Start listening for characters.
        OSX_JAVA_GETCHAR.register;

        % Make sure the Matlab window has keyboard focus:
        if ~IsWin && exist('commandwindow') %#ok<EXIST>
            % Call builtin implementation:
            commandwindow;
            drawnow;
        end

        % Should we block output of characters to Matlab?
        if listenFlag > 1 || listenFlag == -1
            % Disable redispatching:
            OSX_JAVA_GETCHAR.setRedispatchFlag(1);
        else
            % Enable redispatching: This is the startup default.
            OSX_JAVA_GETCHAR.setRedispatchFlag(0);
        end
    else
        % Stop listening for characters and clear the buffer.
        OSX_JAVA_GETCHAR.unregister;
        OSX_JAVA_GETCHAR.clear;
        % Enable redispatching:
        OSX_JAVA_GETCHAR.setRedispatchFlag(0);
    end

    % On non-Vista we're done. On Vista and later, ie. on all versions we
    % still support, we fall-through to the fallback path below, as Java
    % based GetChar() is only useful to suppress character output to the
    % Matlab command window, aka clutter prevention, not for actually
    % recording key strokes. If we are running on Linux with the KDE
    % desktop GUI, we also need to use non-Java fallbacks for keystroke
    % recording, as KDE's window manager has the nasty habit of removing
    % keyboard input focus from the Matlab window, as soon as the onscreen
    % window opens, so Java based GetChar doesn't get input.
    if ~IsWin && isempty(getenv('KDE_FULL_SESSION'))
        return;
    end

    % Windows-10 or later with old Matlabs Java based GUI running.

    % If only the keyboard was blocked via listenFlag -1 before and
    % this is a unblock request via listenFlag 0 then we are done, as
    % setRedispatchFlag() above has already done all the work.
    if listenFlag == 0 && keyboard_blocked
        keyboard_blocked = 0;
        return;
    end

    % If just keyboard blocking was requested via listenFlag -1 then
    % setRedispatchFlag() above has done all the work and we just need
    % to note this in keyboard_blocked, other than that we are done.
    if listenFlag == -1
        keyboard_blocked = 1;
        return;
    end
end

% Running either on Octave, or on Matlab command line, or on Matlab R2025a+, or on MS-Windows:

% Does the user only want to block keyboard input from spilling into the
% command window / console, but not use GetChar et al.?
if listenFlag == -1
    % Yes. No need for use of keyboard queues, only use low level tricks
    % to prevent keyboard input on Octave or on Matlab nojvm.
    keyboard_blocked = 1;

    % LoadPsychHID is needed on MS-Windows. It no-ops if called redundantly:
    LoadPsychHID;

    % Disable character forwarding to console:
    PsychHID('KeyboardHelper', -12);

    return;
end

% Does the user only want to unblock keyboard input from spilling into the
% command window / console after blocking it before with listenFlag -1?
if listenFlag == 0 && keyboard_blocked
    % Yes. Just re-enable keyboard input
    keyboard_blocked = 0;

    % LoadPsychHID is needed on MS-Windows. It no-ops if called redundantly:
    LoadPsychHID;

    % Re-enable character forwarding to console:
    PsychHID('KeyboardHelper', -10);

    return;
end

if keyboard_blocked
    % Reset keyboard only blocked state on other transitions:
    keyboard_blocked = 0;
end

% On all systems we prefer to (ab)use keyboard queues. This allows
% character suppression via ListenChar(2) to work at least on macOS and
% Linux with old Matlab pre R2025a and provides high robustness against
% keyboard focus changes. If we can't get the relevant keyboard queue on
% macOS or Windows at this point, we have to fail. However, if we are on
% Linux and the keyboard queue is already in use by usercode, we can fall
% back to 'GetMouseHelper' low-level terminal tty magic. The only downside
% is that typed characters will spill into the console, ie., ListenChar(2)
% suppression is unsupported:
if ~IsLinux || ~KbQueueReserve(3, 2, [])
    % We can use the default keyboard's keyboard queue - Good:

    % LoadPsychHID is needed on MS-Windows. It no-ops if called redundantly:
    LoadPsychHID;

    if listenFlag > 0
        % Only need to reserve/create/start queue if we don't have it
        % already:
        if ~KbQueueReserve(3, 1, [])
            % Try to reserve default keyboard queue for our exclusive use:
            if ~KbQueueReserve(1, 1, [])
                % This is non-fatal, only worth a warning:
                if IsOSX
                    % macOS:
                    warning('PTB3:KbQueueBusy', 'Keyboard queue for default keyboard device already in use by KbQueue/KbEvent functions et al. Use of ListenChar(2) may work for keystroke suppression, but GetChar() etc. will not work.\n');
                else
                    % MS-Windows:
                    warning('PTB3:KbQueueBusy', 'Keyboard queue for default keyboard device already in use by KbQueue/KbEvent functions et al. Use of ListenChar/GetChar/CharAvail/FlushEvents etc. and keyboard queues is mutually exclusive!');
                end

                % We fall through to KeyboardHelper to enable input
                % redirection on macOS. While our CharAvail() and
                % GetChar() are lost causes, input redirection and CTRL+C
                % can work if usercode has called KbQueueStart, as the
                % users kbqueue-thread gives us a free-ride for our
                % purpose.
            else
                % Got it. Allocate and start it:
                PsychHID('KbQueueCreate');
                PsychHID('KbQueueStart');
            end
        end
    else
        % Does default keyboard queue belong to us?
        if KbQueueReserve(3, 1, [])
            % Yes. Stop and release it:
            PsychHID('KbQueueStop');
            PsychHID('KbQueueRelease');
            KbQueueReserve(2, 1, []);
        end
    end

    if (listenFlag > 1) && (~IsOSX || ~IsOctave || IsGUI)
        % Disable character forwarding to console:
        PsychHID('KeyboardHelper', -12);
    elseif (listenFlag == 1) && (IsOctave && IsGUI)
        % Enable character forwarding to the runtime/console.
        % This is special: We receive our characters via the KbQueues event
        % buffer. At the same time, the runtime receives characters via
        % stdin, which are fed by our Kbqueue and a special unix pipe:
        PsychHID('KeyboardHelper', -11);
    else
        % Enable character forwarding to console,
        % disable it for us, as we use keyboard
        % queues, not tty magic:
        PsychHID('KeyboardHelper', -10);
    end

    return;
end

% This fallback code is only executed on Linux as a last resort. It uses
% low-level tty magic to get some characters from the stdin stream of the
% controlling tty:
if listenFlag > 1
    % Disable character forwarding to console: This will also prevent us
    % from seeing any characters, as we can only see what the runtime aka
    % stdin/tty sees - which is nothing. Can't use kbqueues to get the data
    % as usual:
    PsychHID('KeyboardHelper', -12);
    warning('PTB3:KbQueueBusy', 'Keyboard queue for default keyboard device already in use by KbQueue/KbEvent functions et al.\nUse of ListenChar(2) may work for keystroke suppression, but GetChar() etc. will not work!\n');    
elseif listenFlag == 1
    % Enable character forwarding to the runtime/console. Runtime gets
    % keystrokes and CTRL+C, but we will only get data via stdin/tty if the
    % runtime doesn't steal characters from us, as we can't use kbqueues
    % here:
    PsychHID('KeyboardHelper', -11);
    warning('PTB3:KbQueueBusy', 'Keyboard queue for default keyboard device already in use by KbQueue/KbEvent functions et al.\nUse of ListenChar(2) may work for keystroke suppression, but CharAvail()/GetChar() will not work reliably under all conditions!\n');
else
    % Enable character forwarding to console. Runtime gets keystrokes and
    % CTRL+C, but we will only get data via stdin/tty if the runtime
    % doesn't steal characters from us, as we can't use kbqueues here:
    PsychHID('KeyboardHelper', -10);
    warning('PTB3:KbQueueBusy', 'Keyboard queue for default keyboard device already in use by KbQueue/KbEvent functions et al.\nUse of ListenChar(2) may work for keystroke suppression, but CharAvail()/GetChar() will not work reliably under all conditions!\n');
end

return;
