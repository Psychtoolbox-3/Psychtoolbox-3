function ret = DisableKeysForKbCheck(disablekeys)
% olddisabledkeys = DisableKeysForKbCheck([disablekeys])
%
% Specify a vector of keycodes for keys which should be
% ignored by KbCheck and KbWait. This is useful to
% disable specific keys, e.g., if they are stuck. The function
% returns the old set of disabled keys.
%
% Example: To disable the keys with keycodes 4, 6 and 7, do
% DisableKeysForKbCheck([4, 6, 7]);
%
% Caution: This setting is reset to "empty" during a "clear all" command!
%
% Background info:
%
% Some users of Laptops experienced the problem of "stuck keys": Some keys
% are always reported as "down", so KbWait returns immediately and KbCheck
% always reports keyIsDown == 1. This is often due to special function keys.
% These keys or system functionality are assigned vendor specific
% key codes, e.g., the status of the Laptop lid (opened/closed) could be
% reported by some special keycode. Whenever the Laptop lid is open, this key
% will be reported as pressed. You can work around this problem by passing
% a list of keycodes to be ignored by KbCheck and KbWait.
% _________________________________________________________________________
% 
% See also: FlushEvents, KbName, KbDemo, KbWait, KbCheck, GetChar, CharAvail.

% History:
% 10/24/06 Written (MK).

% This global variable is used to pass the vector of disabled keys
% to the KbCheck command:
global ptb_kbcheck_disabledKeys;

% Return old vector:
ret = ptb_kbcheck_disabledKeys;

if nargin>0
   % Set new vector:
   ptb_kbcheck_disabledKeys = disablekeys;
end

return;
