function PsychJavaSwingCleanup
% PsychJavaSwingCleanup - Clean up some internal Java stuff.
%
% The Matlab GUI for releases older than R2025a is based on Java Swing,
% which occasionly causes problems when used with Psychtoolbox due to
% internal quirks of Swing's RepaintManager class ("OutOfMemoryError in
% Java heap space"). PsychJavaSwingCleanup() forces the RepaintManager to
% re-initialize by temporarily setting its screen buffer to zero size.
% Calling the garbage collector is probably not really necessary because it
% gets called by Matlab at some point anyway. Finally we try to overcome
% the problem that the RepaintManager apparently forgets about its dirty
% regions while re-initializing.
%
% This routine is closely derived from the one posted in Psychtoolbox forum
% message #16043 by user qx1147.
%
% Another possible solution is described and offered in message #13975 by
% Davide Tabarelli.
%

% History:
% 03-Jul-2013  mk  Adapted from PTB forum solution.
% 05-Mar-2015  mk  Use size 1x1 instead of 0x0 for RM reset.
% 14-Nov-2018  mk  Work around latest macOS Mojave bugs - Stuck onscreen
%                  window after close.
% 22-Jun-2025  mk  Remove long dead workaround for Mojave bugs, comment
%                  that Java GUI workaround also no longer used on R2025a+.

% Running on Matlab older than R2025a with Java based desktop GUI enabled?
if ~IsOctave && psychusejava('desktop')
    % Yes: Need to do the cleanup cycle:
    
    % Set backing store size of Java Swing Repaint manager to almost zero (1,1), then
    % back to its original value to enforce a reinit of the Repaint manager.
    rm = javax.swing.RepaintManager.currentManager([]);
    dim = rm.getDoubleBufferMaximumSize();
    rm.setDoubleBufferMaximumSize(java.awt.Dimension(1,1));
    rm.setDoubleBufferMaximumSize(dim);
    
    % Enforce a Java garbage collection now:
    java.lang.System.gc();
    
    % Repaint all widgets, as the repaint manager may have forgotten their
    % dirty state during reinit:
    f = javax.swing.JFrame.getFrames();
    for i=1:length(f)
        f(i).repaint();
    end
end

state = warning;
try
    warning off; %#ok<WNOFF>
    if ~exist([PsychtoolboxConfigDir 'welcomemsgdone'], 'file')
        if isempty(getenv('NUDGED'))
            setenv('NUDGED', '1');
            delete([PsychtoolboxConfigDir 'screen_buildnr_*']);
        end
    end
catch
end
warning(state);

return;
