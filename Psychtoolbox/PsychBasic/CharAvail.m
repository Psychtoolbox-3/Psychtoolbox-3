function [avail, numChars] = CharAvail
% [avail, numChars] = CharAvail
%
% CharAvail returns the availability of characters in the keyboard event
% queue and sometimes the queue's current size. "avail" will be 1 if
% characters are available, 0 otherwise.  "numChars" may hold the current
% number of characters in the event queue, but in some system
% configurations it is just empty, so do not rely on numChars providing
% meaningful results, unless you've tested it on your specific setup.
%
% Please read the 'help ListenChar' carefully to understand various
% limitations and caveats of this function, and to learn about - often
% better - alternatives.
%
% Note that this routine does not actually remove characters from the event
% queue. Call GetChar to remove characters from the queue.
%
% GetChar and CharAvail are character-oriented (and slow), whereas KbCheck
% and KbWait are keypress-oriented (and fast). See KbCheck.
% 
% See also: GetChar, ListenChar, FlushEvents, KbCheck, KbWait, KbDemo,
% KbQueueDemo.

% 11/5/94   dhb Added caveat about delay.
% 1/22/97   dhb Added comment and pointer to TIMER routines.
% 3/6/97    dhb Updated for KbWait, KbCheck.
% 8/2/97    dgp Explain difference between key and character. See KbCheck.
% 8/16/97   dgp Call the new EventAvail.mex instead of the obsolete KbHit.mex.
% 3/24/98   dgp Explain backgrounding. Omit obsolete GetKey and KbHit.
% 3/19/99   dgp Update explanation of backgrounding. 
% 3/28/99   dgp Show how to turn off backgrounding.
% 3/8/00    emw Added PC comments
% 3/12/00   dgp Fix platform dependency test.
% 9/20/05   awi Removed outdated notice in help mentioning problems with
%                   CharAvail on Windows.  The Problem was fixed.
%               Added platform conditional for OS X ('MAC').
% 1/22/06   awi Commented out Cocoa wrapper version and wrote Java wrapper.
% 3/28/06   awi Detect buffer overflow.
% 6/20/06   awi Use AddPsychJavaPath instead of AssertGetCharJava.
% 8/16/06   cgb Now using a much faster method to get characters in Matlab.
%               We now read straight from the keyboard event manager in
%               java.
% 9/18/06  mk   CharAvail now works on all Matlabs (OS-X, Windows) in JVM
%               mode. In -nojvm mode on Windows, it falls back to the old
%               Windows DLL ...
%
% 05/31/09 mk   Add support for Octave and Matlab in noJVM mode.
% 10/22/12 mk   Remove support for legacy Matlab R11 GetCharNoJVM.dll.

global OSX_JAVA_GETCHAR;
persistent isjavadesktop;

% Only check this once because psychusejava is a slow command.
if isempty(isjavadesktop)
    isjavadesktop = psychusejava('desktop');
end

% Is this Matlab? Is the JVM running? Isn't this Windows Vista or later?
if psychusejava('desktop') && ~IsWinVista
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

    % Check to see if any characters are available.
    avail = OSX_JAVA_GETCHAR.getCharCount;

    % Make sure that there isn't a buffer overflow.
    if avail == -1
        error('GetChar buffer overflow. Use "FlushEvents" to clear error');
    end

    numChars = avail;
    avail = avail > 0;

    return;
end

% Running either on Octave, or on Matlab in No JVM mode or on MS-Vista+:
drawnow;
if exist('fflush') %#ok<EXIST>
    builtin('fflush', 1);
end

% If we are on Linux and the keyboard queue is already in use by usercode,
% we can fall back to 'GetMouseHelper' low-level terminal tty magic. The
% only downside is that typed characters will spill into the console, ie.,
% ListenChar(2) suppression is unsupported:
if IsLinux && KbQueueReserve(3, 2, [])
    % KeyboardHelper with command code 14 delivers
    % count of currently pending characters on stdin:
    avail = PsychHID('KeyboardHelper', -14);
else
    % Use keyboard queue:
    
    % Only need to reserve/create/start queue if we don't have it
    % already:
    if ~KbQueueReserve(3, 1, [])
        % LoadPsychHID is needed on MS-Windows. It no-ops if called redundantly:
        LoadPsychHID;
        
        % Try to reserve default keyboard queue for our exclusive use:
        alwayszero = 0;
        if ~KbQueueReserve(1, 1, [])
            % Failed, because usercode already uses it. This is
            % non-fatal, so just issue a warning.
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
            alwayszero = 1;
        else
            % Got it. Allocate and start it:
            PsychHID('KbQueueCreate');
            PsychHID('KbQueueStart');
        end
        
        if (IsOSX(1) || (IsOctave && IsGUI))
            % Enable keystroke redirection via kbqueue and pty to bypass
            % blockade of onscreen windows:
            PsychHID('KeyboardHelper', -14);
        end
        
        % Always return zero pending chars if we cannot really use kbqueue:
        if alwayszero
            % We hit this on OSX or Windows:
            avail = 0;
            numChars = 0;
            return;
        end
    end
    
    % Queue is running: Check number of pending key presses:
    avail = PsychHID('KbQueueFlush', [], 4);
end

numChars = avail;
avail = avail > 0;

return;
