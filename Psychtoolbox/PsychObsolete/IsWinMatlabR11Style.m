function resultFlag = IsWinMatlabR11Style
% resultFlag = IsWinMatlabR11Style
%
% Return true if Psychtoolbox is running on Microsoft Windows,
% with a Matlab version older than R2007a, i.e. using MEX files that
% are compiled against Matlab R11.
% 
% See also: IsOSX, IsLinux, IsOS9, IsWin, OSName, computer

% HISTORY
% 1/03/10 mk    Wrote it.

persistent rc;

if isempty(rc)
    % Running on Windows, but not on Octave (ie., on Matlab) and the Screen
    % mex file has a .dll extension? If so we're on a R11 based system:
    if IsWin & ~IsOctave & ~isempty(findstr(which('Screen'), 'Screen.dll')) %#ok<AND2>
        rc = 1;
    else
        rc = 0;
    end
end;
resultFlag = rc;
