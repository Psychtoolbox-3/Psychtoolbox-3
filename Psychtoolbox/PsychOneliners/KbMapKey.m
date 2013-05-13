function q = KbMapKey(KeyNums,keyCode)
% q = KbMapKey(KeyNums,keyCode)
% Checks if any of specified keys is depressed
% input: vector of key numbers, keystate boolean vector as returned by e.g.
%        KbCheck.
% output: boolean vector specifying true for every key when its depressed
%
% SEE ALSO: KbCheck

% 2008     DN  Wrote it.
% 06-02-09 DN  A vector is now accepted as input
% 08-02-09 DN  Mario's advice: changed to get keyCode (returned by, e.g.
%              KbCheck) as input instead of call KbCheck itself -> more
%              flexible

q = ismember(KeyNums,find(keyCode));
