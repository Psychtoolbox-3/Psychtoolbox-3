function PsychJavaSwingCleanup
% PsychJavaSwingCleanup - Clean up some internal Java stuff.
%
% The Matlab GUI is based on Java Swing, which occasionly causes problems
% when used with Psychtoolbox due to internal quirks of Swing's
% RepaintManager class ("OutOfMemoryError in Java heap space").
% PsychJavaSwingCleanup() forces the RepaintManager to re-initialize by
% temporarily setting its screen buffer to zero size. Calling the garbage
% collector is probably not really necessary because it gets called by
% Matlab at some point anyway. Finally we try to overcome the problem that
% the RepaintManager apparently forgets about its dirty regions while
% re-initializing.
%
% This routine is closely derived from the one posted in Psychtoolbox forum
% message #16043 by user qx1147.
%
% Another possible solution is described and offered in message #13975 by
% Davide Tabarelli.
%

% History:
% 3-July-2013  mk  Adapted from PTB forum solution.

% Running on Matlab with Java based desktop GUI enabled? Otherwise we no-op.
if psychusejava('desktop')
    % Yes: Need to do the cleanup cycle:
    
    % Set backing store size of Java Swing Repaint manager to zero (0,0), then
    % back to its original value to enforce a reinit of the Repaint manager.
    rm = javax.swing.RepaintManager.currentManager([]);
    dim = rm.getDoubleBufferMaximumSize();
    rm.setDoubleBufferMaximumSize(java.awt.Dimension(0,0));
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

return;
