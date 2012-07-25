function osNameStr = OSName
% sysNameStr = OSName
%
% Return the convential English-language name for your operating system (OS).
% OSName is useful in constructing error message strings which refer to
% particular operating systems. System name strings as returned by the MATLAB
% command "computer" are unsuitable for this purpose because they are
% abbreviations, not names.  
% 
% Currently possible returned namestrings are 'Windows', 'Linux' or 'OSX'.
%
% see also: computer, IsOSX, IsWin, IsLinux, MacModelName, DescribeComputer

% HISTORY
%
% 3/5/06  awi  Wrote it.  For use in ListenChar.
% 9/20/09  mk  Updated.

if IsWin
    osNameStr='Windows';
elseif IsLinux
    osNameStr='Linux';
elseif IsOSX
    osNameStr='OSX';
else
    osNameStr='Unknown OS';
end

return;
