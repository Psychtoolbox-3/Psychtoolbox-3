function resultFlag = IsOSX(is64)

% resultFlag = IsOSX([is64=0])
%
% Returns true if the operating system is Mac OS X.
% If the optional 'is64' flag is set to one, returns
% true if the runtime is 64 bit and on OSX.
%
% See also: IsWin, IsLinux, OSName, computer

% HISTORY
% ??/??/?? awi Wrote it.
% 6/30/06  awi Fixed help section.
% 4/30/12  mk  Support 64-Bit OS/X.

persistent rc;
persistent rc64;

% check input
if nargin < 1 || isempty(is64)
     is64 = 0;
end


if isempty(rc)
     rc = streq(computer,'MAC') || streq(computer, 'MACI') || streq(computer, 'MACI64') || ~isempty(strfind(computer, 'apple-darwin'));
end

if isempty(rc64)
     rc64 = rc && ~isempty(strfind(computer, '64'));
end

if is64 == 0
     resultFlag = rc;
else
     resultFlag = rc64;
end
