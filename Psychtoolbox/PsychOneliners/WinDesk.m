function WinDesk
% minimizes all windows, thus shows windows desktop.
% Not the same as Win+D, the same as Win+M
% All open windows are minimized except (those with) open dialogue boxes

% DN 2008-06    Wrote it
% DN 2008-12-06 Added test if on Windows

if IsWin && ~IsOctave % from PTB - test if windows as we use a COM Automation server (Octave doesn't have those either)
    shell = actxserver('Shell.Application');
    shell.MinimizeAll;
    shell.delete;         % release interface
end
