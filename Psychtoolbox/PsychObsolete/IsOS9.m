function resultFlag = IsOS9

% resultFlag = IsOS9
%
% Returns true if the operating system is Mac OS 9.
% 
% See also: IsOSX, IsWin, IsLinux, OSName, computer

% HISTORY
% ??/??/?? awi Wrote it.
% 6/30/06  awi Fixed help section.  

persistent rc;

if isempty(rc)
     rc= streq(computer,'MAC2');
end;
resultFlag=rc;
