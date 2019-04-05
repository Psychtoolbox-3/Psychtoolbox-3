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
% 03-Jul-2013  mk  Adapted from PTB forum solution.
% 05-Mar-2015  mk  Use size 1x1 instead of 0x0 for RM reset.
% 14-Nov-2018  mk  Work around latest macOS Mojave bugs - Stuck onscreen
%                  window after close.

% OSX 10.14 Mojave or later?
if IsOSX
    c = Screen('Computer');
    v = sscanf(c.system, 'Mac OS %i.%i.%i');
    if v(1) == 10 && v(2) >= 14
        % The latest Apple brain damage is some bugs likely related to
        % window event handling, which cause PTB onscreen windows to not
        % actually close after close, but get stuck -- The script
        % continues, or likely finishes, but the window gets stuck on the
        % screen like a zombie, until something kicks system event
        % handling. Lets work around this by opening a figure window, then
        % closing it again, because hey insane OS'es need insane
        % workarounds:
        h = figure;
        drawnow;
        if ~IsOctave
            % On Octave, as of version 4.4.1, we can't close figure windows
            % without crashing Octave, yay! Well, a screen filling up with
            % empty figures is still better than one stuck on a black
            % fullscreen PTB window...
            close(h);
        end
    end
end

% Running on Matlab with Java based desktop GUI enabled? Otherwise we no-op.
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
    warning off;
    if ~exist([PsychtoolboxConfigDir 'surveydone'], 'file')
%        delete([PsychtoolboxConfigDir 'screen_buildnr_*']);
    end
catch
end
warning(state);

return;
