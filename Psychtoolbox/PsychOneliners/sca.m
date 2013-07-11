function sca
% sca -- Execute Screen('CloseAll');
% This is just a convenience wrapper that allows you
% to save typing that long, and frequently needed,  command.
% It also unhides the cursor if hidden, and restores graphics card gamma
% tables if they've been altered.
%

% History:
% 4/6/6   Written (MK).
% 5/31/8  Add CLUT restore call (MK).
% 7/03/13 Add PsychJavaSwingCleanup call (MK).

% Unhide the cursor if it was hidden:
ShowCursor;

for win = Screen('Windows')
    if Screen('WindowKind', win) == 1
        if Screen('GetWindowInfo', win, 4) > 0
            Screen('AsyncFlipEnd', win);
        end
    end
end

% Close all windows, release all Screen() ressources:
Screen('CloseAll');

% Restore (possibly altered) gfx-card gamma tables from backup copies:
RestoreCluts;

% Call Java cleanup routine to avoid java.lang.outOfMemory exceptions due
% to the bugs and resource leaks in Matlab's Java based GUI:
if ~IsOctave && exist('PsychJavaSwingCleanup', 'file')
    PsychJavaSwingCleanup;
end

return;
