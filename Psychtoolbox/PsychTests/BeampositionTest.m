function BeampositionTest(n, screenid, synced)
% BeampositionTest - Test GPU scanout position ("beamposition") queries.
%
% Usage: BeampositionTest([n=10000][, screenid=max][, synced=0])
%
% Opens a onscreen window on screen 'screenid' (default = maximum screenid)
% and runs a sample loop, performing 'n' (default 10000) beamposition
% queries, either as fast as it can if 'synced' == 0, or it performs
% Screen('Flip') and takes one sample after each flip, ie., in sync with
% vertical retrace.
%
% Then it plots its samples in multiple figures:
% Beamposition vs. sample count, Beamposition vs. time, a histogram of all
% sampled beampositions, and the duration of each query in msecs.
%
% Here's what you should see on a properly working system:
%
% Beamposition against time or samplecount should show a periodic sawtooth pattern,
% where scanout starts at zero and increments up to the value VBL Endline,
% as reported by PTB's Matlab output. The histogram should show an even
% distribution (a flat line) across all bins from bin 0 to bin "VBL
% Endline". Query times should be in the low microsecond range, e.g.,
% not much more than 20-30 usecs on a modern machine on average.
%
% Other things you might see:
%
% Periodic spikes in the query times and possibly some gaps, or total loss
% of signal in the histogram at values higher than "VBL Startline" and or a
% spike in the histogram for beamposition bin zero: Your graphics driver or
% hardware has a small bug for which PTB can compensate, at the cost of
% higher query times and thereby slightly degraded realtime behaviour of
% your system. This is usually accompanied by some warnings wrt. slightly
% broken beamposition queries.
%
% No spikes in query times but no values abouve "VBL Startline" at all in
% the histograms and plots: A slightly broken graphics driver, but PTB
% couldn't compensate for the bug because it didn't detect the problem.
% Manual ways to enable the workaround are explained in 'help
% ConserveVRAMSettings', section 'kPsychUseBeampositionQueryWorkaround'.
%
% Beampositions implausibly high or low for your display. PTB may query the
% wrong display on a multi-display setup. See 'help DisplayOutputMappings'
% for remedies.
%
% A flat-line of all-zero beampositions: Beamposition queries are
% unsupported, disabled or broken on your setup. If you know they should
% work, then PTB may query the wrong display on a multi-display setup. See
% 'help DisplayOutputMappings' for remedies. Another possible fix is
% installing and loading the PsychtoolboxKernelDriver on OSX ('help
% PsychtoolboxKernelDriver') or running PsychLinuxConfiguration on Linux
% and rebooting the machine once to enable our own beamposition query
% mechanism. On Linux this is not needed when using the free and
% open-source graphics drivers, as these have a builtin timestamping
% mechanism which is even more robust and precise than beampositionq
% queries. On Linux with NVidia or AMD proprietary binary drivers this is
% always needed once for a given machine. On OSX it is increasingly needed
% due to total lack of OS support on Intel and AMD graphics cards and
% increasingly broken support for NVidia graphics cards. On OSX with NVidia
% graphics cards you additionally need to apply the measures explained in
% 'help ConserveVRAMSettings' - the section about kPsychDontUseNativeBeamposQuery.
%
% A flat-line not at zero, but at some fixed value. See previous paragraph:
% Install the kernel drivers or setup Linux properly, apply the kPsychDontUseNativeBeamposQuery
% fix on OSX, possibly apply measures from 'help DisplayOutputMappings'.
%
% Weird looking non-sawtooth patterns: This is a totally broken graphics
% driver, often to be found on Apple machines running OSX with Apple Retina
% displays or other LCD flat panels. See previous paragraph for fixes.
%
% -> Most fixes involve switching to our own implementation of beamposition
% queries instead of using broken or missing OS implementations, and then
% some fiddling around.
%
% -> On MS-Windows there's usually no way to fix broken queries beyond the
% automatic fixes applied by Screen() if it detects a problem.
%

% History:
% ??/??/????  mk     Written.
% 10.01.2013  iandol Cleanup of code, fixes, beautification of plots.
%

AssertOpenGL;

if nargin < 1 || isempty(n)
    n = 10000;
end

if nargin < 2 || isempty(screenid)
    screenid = max(Screen('Screens'));
end
fprintf('Testing screen #: %g\n',screenid);

if nargin < 3 || isempty(synced)
    synced = 0;
end

w=Screen('OpenWindow', screenid, 0);
Screen('Flip',w);
beampos = zeros(n, 1);
t = zeros(n, 1);
for i=1:n
    if synced
        Screen('Flip', w);
    end
    beampos(i) = Screen('GetWindowInfo', w, 1);
    t(i) = GetSecs;
end

% Extract VBLEndline for number of bins in histogram plot:
winfo = Screen('GetWindowInfo', w);
VBLEndline = winfo.VBLEndline;

% Skip this test if VBLEndline == -1 aka beamposition queries broken or
% disabled. Otherwise we would hang infinitely here:
if VBLEndline > -1
    while Screen('GetWindowInfo', w, 1) < 500; end;
    ta=GetSecs;
    while Screen('GetWindowInfo', w, 1) < 515; end;
    te=GetSecs;
    
    tel = (te - ta) * 1000
else
    fprintf('WARNING: Beamposition queries broken or disabled or unsupported on your setup!\n');
end

Screen('CloseAll');

figure;
subplot(2,2,1)
plot(beampos,'k-');
axis tight
title(['BeamPosition: Screen ' num2str(screenid)]); ylabel('Scanlines'); xlabel('Loop #')
subplot(2,2,2)
plot(diff(t)*1000,'k.');
axis tight
title('TimeStamp Deltas'); ylabel('Time (ms)'); xlabel('Loop #')
subplot(2,2,3)
hist(beampos, VBLEndline + 10);
axis tight
title('Histogram of BeamPosition'); ylabel('Number'); xlabel('Scanline #')
subplot(2,2,4)
plot(t, beampos,'k-');
axis tight
title('BeamPosition against Timestamp'); ylabel('Scanlines'); xlabel('Time (s)')
return;
