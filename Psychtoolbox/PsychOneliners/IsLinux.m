function resultFlag = IsLinux(is64)

% resultFlag = IsLinux([is64=0])
%
% Returns true if the operating system is Linux.
% If optional 'is64' flag is set to one, returns
% true if the runtime is 64 bit and on Linux.
% 
% See also: IsOSX, IsWin, IsOS9, OSName, computer

% HISTORY
% ??/??/?? ??? Wrote it.
% 6/30/06  awi Fixed help section.  
% 6/18/11   mk Handle query if 64 bit runtime as well.

persistent rc;
persistent rc64;

if isempty(rc)
     rc= streq(computer,'GLNX86') | streq(computer,'GLNXA64') | ~isempty(strfind(computer, 'linux-gnu'));
end

if isempty(rc64)
     rc64 = rc & ~isempty(strfind(computer, '64'));
end

if nargin < 1
     resultFlag = rc;
     return;
end

if isempty(is64)
     resultFlag = rc;
     return;
end

if is64 == 0
     resultFlag = rc;
else
     resultFlag = rc64;
end
