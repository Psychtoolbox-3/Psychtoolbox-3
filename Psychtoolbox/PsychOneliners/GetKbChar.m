function [ch, when] = GetKbChar(varargin)
% [ch, when] = GetKbChar([deviceIndex][, untilTime=inf][, optional KbCheck arguments...]);
%
% GetKbChar reimplements basic functionality of GetChar() by use of KbCheck
% and KbPressWait. It accepts optionally all arguments that KbCheck accepts
% and passes those arguments to KbCheck, e.g., a keyboard index in order to
% only query a specific keyboard for input.
%
% GetKbChar also returns with empty return values if the optional deadline
% 'untilTime' is reached.
%
% The function only recognizes standard alpha-numeric keys, i.e., letters
% and numbers, and a few special symbols like the ones on top of the
% numeric keys. It only recognizes the delete, space and return keys as
% special function keys, not other keys like Function keys, CTRL, ALT or
% cursor keys. It always assumes a US keyboard mapping.
%
% It polls the keyboard, so may miss very brief keystrokes and doesn't use
% the keyboard queue.
%
% Use this function if you need a GetChar like interface for simple string
% and number input in situations where GetChar doesn't work reliably, e.g.,
% on some Octave configurations, with Matlab in -nojvm mode or on
% MS-Windows Vista or Windows-7.
%

% History:
% 22.10.2010  mk  Wrote it.

persistent shiftkeys;

if isempty(shiftkeys)
    shiftkeys = [KbName('LeftShift'), KbName('RightShift')];
end

% Get keypress, KbCheck style:
when = KbPressWait(varargin{:});
keycode = zeros(1,256);
down = 1;
secs = when;

if length(varargin) < 2
    untilTime = inf;
else
    untilTime = varargin{2};
    if isempty(untilTime)
        untilTime = inf;
    end
end

while down && (secs < untilTime)
    [down, secs, keyincode] = KbCheck(varargin{:});
    if down
        keycode = keycode + keyincode;
        secs = WaitSecs('YieldSecs', 0.001);
    end
end

% Force keycode to 0 or 1:
keycode(keycode > 0) = 1;

% Shift pressed?
if any(keycode(shiftkeys))
    shift = 2;
else
    shift = 1;
end

% Remove shift keys:
keycode(shiftkeys) = 0;

% Translate to ascii style:
ch = KbName(keycode);

% If multiple keys pressed, only use 1st one:
if iscell(ch)
    ch = ch{1};
end

% Decode 1st or 2nd char, depending if shift key was pressed:
if length(ch) == 1
    if shift > 1 && ismember(ch, 'abcdefghijklmnopqrstuvwxyz')
        ch = upper(ch);
    end
elseif length(ch) == 2
    ch = ch(shift);
elseif length(ch) > 2
    if strcmpi(ch, 'Return')
        ch = char(13);
    end
    if strcmpi(ch, 'ENTER')
        ch = char(13);
    end
    if strcmpi(ch, 'space')
        ch = char(32);
    end
    if strcmpi(ch, 'DELETE') || strcmpi(ch, 'BackSpace')
        ch = char(8);
    end
    
    % Catch otherwise unhandled special keys:
    if length(ch) > 1
        % Call ourselves recursively, thereby discarding this unmanageable
        % result.
        fprintf('GetKbChar: Warning: Received keypress for key %s. Don''t really know what to do with it.\n', ch);
        fprintf('Maybe you should check for stuck or invalid keys?\n');
        [ch, when] = GetKbChar(varargin{:});
        return;
    end
end

% Done.
return;
