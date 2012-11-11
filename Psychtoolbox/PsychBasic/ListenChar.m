function ListenChar(listenFlag)
% function ListenChar([listenFlag])
%
% Tell the Psychtoolbox function "GetChar" to start or stop listening for
% keyboard input.  By default ListenChar listens when no argument is
% supplied.  Passing 0 will turn off character listening and reset the
% buffer which holds the captured characters. Passing a value of 1 or not
% passing any value will enable listening. Passing a value of 2 will enable
% listening, additionally any output of keypresses to Matlabs or Octaves
% windows is suppressed. Use this with care, if your script aborts with an
% error, Matlab or Octave may be left with a dead keyboard until the user
% presses CTRL+C to reenable keyboard input. 'listenFlag' 2 is silently
% ignored on Matlab in -nojvm mode under MS-Windows.
%
% This function isn't entirely necessary to turn on listening as calling
% GetChar, CharAvail, or FlushEvents will trigger listening on. However,
% it is the only method by which to disable listening or switch between
% suppression of keyboard input to Matlab and unsuppressed mode.
%
% Please note that the commands ListenChar, CharAvail and GetChar are
% subject to various constraints and limitations, depending on the
% operating system you use, if you use Matlab in Java or -nojvm mode, if
% you use Octave, if you have Screen() onscreen windows open or not, or if
% you use KbQueueXXX functions in parallel or not. Therefore use of these
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
% 1. Works very well with Matlab and its Java based GUI enabled on Linux
% and MacOSX, as well as on WindowsXP.
%
% 2. When used on Windows Vista or later (Vista, Windows-7, Windows-8, ...)
% with Matlab's Java GUI, you cannot use any KbQueue functions at the same
% time, ie., KbQueueCreate/Start/Stop/Check/Wait as well as KbWaitTrigger,
% KbEventFlush, KbEventAvail, and KbEventGet are off limits after any call
% to ListenChar, ListenChar(1), ListenChar(2), FlushEvents, CharAvail or
% GetChar. You would need to call ListenChar(0) before you could call
% KbQueueCreate and then use those functions. Vice versa, after a call to
% KbQueueCreate, CharAvail, FlushEvents, and GetChar are dysfunctional and
% ListenChar may be limited. You need to call KbQueueRelease before you can
% use them again. This is generally true on MacOSX, and true for use of the
% default keyboard (the keyboard that GetChar et al. use) on MS-Windows and
% Linux. Use of other keyboards than the default keyboard, or of other
% devices, e.g., mouse or joystick, is not prohibited during use of GetChar
% et al.
%
% 3. If you use Matlab in "matlab -nojvm" mode without its GUI, or if you
% use GNU/Octave instead of Matlab, the same restrictions as in 2. apply -
% no parallel use of the default keyboards KbQueue - or any KbQueue on OSX
% with GetChar et al. The only feature that works in parallel with KbQueues
% is the suppression of spilling of keystroke characters into the Matlab or
% Octave window during ListenChar(2) - at least on Linux and OSX, on
% Windows this can't be prevented at all in "matlab -nojvm" mode. However,
% if you switch to ListenChar(2) mode, you cannot break out of it by
% pressing CTRL+C on Linux if the keyboard queue that is in parallel use
% didn't get KbQueueStart called, ie. if it is stopped. On OSX with a
% stopped Keyboard queue, neither CTRL+C nor stopping a runaway script
% works.
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
% Basically: Mixing GetChar et al. and modern KbQueue functions is usually
% not advisable, or if needed, great care must be taken to sidestep all the
% mentioned limitations. Also the KbQueue functions usually have better
% timing precision and allow to flexibly address multiple keyboards
% separately at least on Linux.
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

global OSX_JAVA_GETCHAR;
 
if nargin == 0
    listenFlag = 1;
elseif nargin > 1
    error('Too many arguments to ListenChar!  See "help ListenChar" for more information');
end

if ~ismember(listenFlag, [0,1,2])
    error('Invalid listenFlag provided!  See "help ListenChar" for more information');
end


% Is this Matlab? Is the JVM running?
if psychusejava('desktop')
    % Java enabled on Matlab. There's work to do.

    % Make sure that the GetCharJava class is loaded.
    if isempty(OSX_JAVA_GETCHAR)
        try
            OSX_JAVA_GETCHAR = AssignGetCharJava;
        catch %#ok<*CTCH>
            error('Could not load Java class GetCharJava! Read ''help PsychJavaTrouble'' for help.');
        end
    end

    if listenFlag
        % Start listening for characters.
        OSX_JAVA_GETCHAR.register;

        % Make sure the Matlab window has keyboard focus:
        if ~IsWinVista && exist('commandwindow') %#ok<EXIST>
            % Call builtin implementation:
            commandwindow;
            drawnow;
        end

        % Should we block output of characters to Matlab?
        if listenFlag > 1
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

    % On non-Vista we're done. On Vista and later, we fall-through to the
    % fallback path below, as Java based GetChar() is only useful to
    % suppress character output to the Matlab command window, aka clutter
    % prevention, not for actually recording key strokes.
    if ~IsWinVista
        return;
    end
end

% Running either on Octave, or on Matlab in No JVM mode, or on a MS-Vista
% system or later.

% On all systems we prefer to (ab)use keyboard queues. This allows character
% suppression via ListenChar(2) to work at least on OSX and Linux and provides
% high robustness against keyboard focus changes. If we can't get the relevant
% keyboard queue on OSX or Windows at this point, we have to fail. However,
% if we are on Linux and the keyboard queue is already in use by usercode,
% we can fall back to 'GetMouseHelper' low-level terminal tty magic. The
% only downside is that typed characters will spill into the console, ie.,
% ListenChar(2) suppression is unsupported:
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
                if IsOSX(1)
                    % OSX:
                    warning('PTB3:KbQueueBusy', 'Keyboard queue for default keyboard device already in use by KbQueue/KbEvent functions et al. Use of ListenChar(2) may work for keystroke suppression, but GetChar() etc. will not work.\n');
                else
                    % 32-Bit OSX, or MS-Windows:
                    warning('PTB3:KbQueueBusy', 'Keyboard queue for default keyboard device already in use by KbQueue/KbEvent functions et al. Use of ListenChar/GetChar etc. and keyboard queues is mutually exclusive!');
                end
                
                % We fall through to KeyboardHelper to enable input
                % redirection on 64-Bit OSX. While our CharAvail() and
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

    if listenFlag > 1
        % Disable character forwarding to console:
        PsychHID('KeyboardHelper', -12);
    elseif (listenFlag == 1) && (IsOSX(1) || (IsOctave && IsGUI))
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
