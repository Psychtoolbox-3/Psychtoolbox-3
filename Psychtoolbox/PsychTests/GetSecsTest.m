function GetSecsTest(n)
% GetSecsTest([n=100000])
%
% This test is meant for Microsoft Windows only!
%
% Performs a reliability test of your systems timing hardware. This script
% tries to find out if your systems clock works correctly, ie., if
% GetSecs(), WaitSecs(), Screen('Flip') and the PsychPortAudio functions
% for timed stimulus onset and clock queries will work correctly.
%
% The optional parameter 'n' specifies how many samples to use for each
% test run. A larger number means longer runtime, but higher accuracy of
% the results. Default is 100000 samples.
%
% The test shows a lot of different data plots, but you should focus on the
% results printed to the Matlab command window. If they sound ok then you
% can probably ignore the plots.
%
% The test exercises your timing hardware: It takes 'n' time samples in a
% loop which repeatedly asks GetSecs for the time. GetSecs is used in a
% special debug-mode which also provides time readings from an alternative
% hardware clock. After sampling, the script compares timestamps from the
% normally used high precision timer against samples taken from the low
% precision reference timer. In an ideal world on a well working system,
% both timers should deliver roughly the same time readings. A fixed offset
% between the two clocks is normal - they both count time from system
% bootup or reset, but they are started with a slight offset at system
% start. A slowly drifting offset is not uncommon on lower quality
% hardware, this is considered ok as long as the drift is small. Important
% is that time is monotonically increasing on both clocks and that the
% ratio of elapsed time in both clocks is nearly the same: Dividing elapsed
% time as reported by the high-res timer by the elapsed time reported by
% the low-res timer should yield a ratio of very close to 1.0. The test
% also checks for regular low-res timer ticks, if pause()'ing Matlabs
% execution for multiple seconds affects the clocks in some worrying
% manner (e.g., time slowing down) or if quick switching between multiple
% cpu cores of a multi-core machine causes time inconsistencies, e.g., time
% going backwards.
%
% The test first runs a block of trials with PTBs timing in "normal mode",
% where PTB applies a couple of fixes for broken hardware.
%
% Then it disables those fixes and retests, to assess how earlier PTB
% versions worked on your system - when no such special workarounds were
% available.
%
% Please note that this test is not 100% water-proof, it could
% theoretically have 'false negatives', as many problems with the clocks
% are associated with the behaviour of system power management, which
% itself is influenced by all kind of factors, e.g., what other
% applications are running in parallel to Matlab, if there is any network
% traffic, what kind of hardware is connected to the system...
%
% We will try to improve the test as we learn more about possible causes of
% trouble.
%
% Starting with PTB releases as of 'beta' from 26th November 2007, PTB also
% performs a couple of runtime checks to spot more timer problems. These
% checks are performed all the time while your scripts are executing, so if
% you see some "CRITICAL-WARNING" messages about timer problems at the
% Matlab prompt while your scripts are executing, better take them serious!
%
% For more (up to date) information about system configurations that might
% suffer from clock problems, background information and troubleshooting
% tips, visit the Psychtoolbox Wiki's FAQ section, specifically:
%
% http://psychtoolbox.org/wikka.php?wakka=FaqGetSecsTestFails
%
% Btw. Currently there are no known problems with timers and clocks on PC
% hardware running recent Linux distributions or on any Apple Macintosh
% computers running MacOS/X ;-)  -- That's why first versions of this test
% are only targeted at MS-Windows.

% History:
% 11/26/2007 Written (MK).

if IsWin
    hardfail = 0;
    softfail = 0;

    % Set default number of samples to 100000:
    if nargin < 1
        n = [];
    end

    if isempty(n)
        n = 100000;
    end

    % Close all plot figure windows:
    close all;
    drawnow;
    clear mex;
    drawnow;

    % Read the boot.ini if possible:
    fid = fopen('C:\boot.ini');
    if fid==-1
        fprintf('Failed to read C:\\boot.ini - not present or inaccessible.\n\n\n');
    else
        bootini = transpose(char(fread(fid)));
        fclose(fid);
        fprintf('Startup config file C:\\boot.ini looks like this:\n');
        disp(bootini);

        if ~isempty(findstr(bootini, '/usepmtimer'))
            fprintf('Seems that your system vendor applied some proper fixes for timing on your system... - Good! ;-)\n\n\n');
        else
            fprintf('No special configuration options found in your startup file. Either not neccessary, or your\n');
            fprintf('system vendor is a moron. We''ll see...\n\n\n');
        end
    end

    % Perform cold init of GetSecs:
    GetSecs;

    Priority(MaxPriority(0));

    fprintf('The first block of tests tests timing precision with the new timing fixes applied...\n\n');
    fprintf('Testing precision and tick rate of the low precision reference timer and high precision timer...\n');
    % Take n samples of all timers in a measurement loop:
    ticks = zeros(1,n);
    raw = zeros(1,n);
    cooked = zeros(1,n);
    for i=1:n
        % The -1 argument tells GetSecs to return unfiltered raw values as
        % well:
        [cooked(i) ticks(i) raw(i)] = GetSecs(-1);
    end

    % Compute increment between successive low-res timer ticks:
    tickinc = diff(ticks);
    plot(tickinc * 1000);
    title('Delta of low-res timer (msecs) between consecutive samples:');
    figure;
    hist(tickinc * 1000, 50);
    title('Distribution of low-res timer deltas (msecs):');

    % Find average duration of a timer tick: This should be 1 msec...
    avgtickdelta = median(tickinc(find(tickinc > 0)));
    fprintf('The average duration of a low-res timer tick is %f milliseconds.\n', avgtickdelta * 1000);
    if abs(avgtickdelta - 0.001) < 0.0001
        fprintf('The low-res timer runs at 1 msec resolution. It is suitable as a fallback-timesource in case\n')
        fprintf('of trouble with the high-res timer. Good!\n\n');
    else
        fprintf('The low-res timer did not switch to 1 msec resolution as requested :( - Not good!\n');
        fprintf('* PTB has to use a much more inefficient implementation of WaitSecs for timed operations.\n');
        fprintf('* The likelihood of defective power management drivers interfering with timing is higher.\n');
        fprintf('* Should PTB need to disable the high res timer, then the low-res timer will only provide inadequate\n')
        fprintf('  resolution for timed stimulus presentation and timing.\n\n');
    end

    overrunticks = length(find(tickinc > avgtickdelta + 0.0001));
    fprintf('The low-res timer showed timing spikes a total of %i times.', overrunticks);

    if overrunticks > 0
        fprintf(' Maximum spike %f msecs.\n\n', max(tickinc)*1000);
        fprintf('The low-res timer seems to be a bit shaky. That will only be a problem if it shows large spikes often,\n');
        fprintf('because it could cause the online correctness tests to fail (false positive) and trigger an unwanted\n');
        fprintf('switch from the high-res timer to the unreliable low-res timer.\n');
        fprintf('Small occasional spikes are not uncommon on such a deficient system as Windows. Many spikes could\n');
        fprintf('indicate that some of your hardware devices or its drivers don''t operate optimally or the system\n');
        fprintf('is severly overloaded in some way.\n\n');
    else
        fprintf('\n\n');
    end

    figure;
    plot(ticks, raw);
    title('Elapsed time (msecs) of low-res (x-axis) vs. high-res (y-axis) timer:');

    figure;
    plot((raw - ticks)*1000);
    title('Lag of high-res timer vs. low-res timer, plotted against sample count:');
    avglag = mean(raw-ticks)*1000;
    lag1 = (raw(1) - ticks(1)) * 1000;
    lagn = (raw(n) - ticks(n)) * 1000;

    % Most important: The ratio between low-res and high-res timer:
    ratio = (raw(n)-raw(1))/(ticks(n)-ticks(1));
    fprintf('lag1 = %15.6f , lag%i = %15.6f - Not too important...\n\n', lag1, n, lagn);
    fprintf('Ratio of elapsed time in high-res vs. low-res clock: %15.10f\n', ratio);

    if abs(ratio - 1.0) > 0.05
        fprintf('\n\nTimers DISAGREE by more than 5 percent! --> This likely indicates an UNRELIABLE high-precision timer!!\n\n');
        hardfail = 1;
    end

    if min(diff(raw)) < 0
        fprintf('TIME WAS RUNNING BACKWARDS!!! The TSC timers are not synchronized across cpu cores and\n');
        fprintf('PTB''s countermeasures were unable to fix this problem! Your system is highly unreliable!\n\n');
        hardfail = 1;
    end

    % Test 3: Let CPU go to sleep for a few seconds and see if there's drift:
    % We use pause() here, because WaitSecs would take counter-measures,
    % rendering our little test non-diagnostic.
    fprintf('\nGoing to sleep for 10 seconds, then checking if this idle time did something bad...\n');
    [cooked1 lowres1 highres1] = GetSecs(-1);
    pause(10);
    [cooked2 lowres2 highres2] = GetSecs(-1);
    fprintf('Elapsed time should be 10 seconds: lowres timer says %f secs, highres timer says %f secs.\n', lowres2-lowres1, highres2-highres1);
    deltaratio = abs(((lowres2-lowres1) / (highres2-highres1)) - 1.0);
    if deltaratio > 0.05
        fprintf('\n\nTimers DISAGREE by more than 5 percent! --> This likely indicates an UNRELIABLE high-precision timer (delta %f %%)!\n\n', deltaratio * 100);
        hardfail = 1;
    end

    % Give a summary of behaviour with new PsychTimeGlue implementation:
    if hardfail > 0
        fprintf('YOUR SYSTEMS TIMING IS BROKEN AND UNRELIABLE!!!\n');
        fprintf('PTB''s built in counter-measures are not an effective fix on your system.\n');
        fprintf('Read in the online help of this test script and on the PTB Wiki about measures\n');
        fprintf('that may allow to fix your systems timing. Apply them, then retest.\n\n');
    else
        fprintf('So far your system behaved well. This could mean it is well behaved.\n');
        fprintf('it could also mean that it is broken, but PTB''s built in work-arounds provide\n');
        fprintf('an effective fix for the brokeness.\n');
        fprintf('Unfortunately it could also mean that the test is unable to detect real problems (false negative) :(\n\n');
    end

    fprintf('\n\nThe second block of tests tries to assess how your system would behave without\n');
    fprintf('PTB''s special builtin work-arounds enabled. This is a hint on how the system behaved\n');
    fprintf('with Psychtoolbox versions prior to 23rd November 2007, and how it would likely behave\n');
    fprintf('with other psychophysics toolkits...\n\n');

    % Disable PTB's fixes:

    % This will disable thread affinity to core zero and allow scheduling to
    % cores zero and one - the first two cores in a multi-core system. It will
    % not do anything on single-core systems. If there is no sync between
    % different cores, this should facilitate detection of the problem:
    GetSecs(1+2);

    % This will set the low-res timer and scheduling rate back to "factory
    % defaults", ie., around 15-16 msecs intervals instead of 1 msec intervals.
    % Given the much lower timer IRQ load, this increases the likelihood of the
    % system entering lower power modes (P-States or C-States), thereby
    % increasing the effects of possible broken power management drivers:
    for j=1:10
        GetSecs(-2);
    end
    pause(1);

    % Take n samples of all timers in a measurement loop:
    ticks = zeros(1,n);
    raw = zeros(1,n);
    cooked = zeros(1,n);
    for i=1:n
        % The -1 argument tells GetSecs to return unfiltered raw values as
        % well:
        [cooked(i) ticks(i) raw(i)] = GetSecs(-1);
    end

    % Compute increment between successive low-res timer ticks:
    tickinc = diff(ticks);
    figure;
    plot(tickinc * 1000);
    title('Block II: Delta of low-res timer (msecs) between consecutive samples:');
    figure;
    hist(tickinc * 1000, 50);
    title('Block II: Distribution of low-res timer deltas (msecs):');

    % Find average duration of a timer tick: This should be 1 msec...
    avgtickdelta = median(tickinc(find(tickinc > 0)));
    fprintf('The average duration of a low-res timer tick is %f milliseconds.\n', avgtickdelta * 1000);
    if abs(avgtickdelta - 0.001) < 0.0001
        fprintf('The low-res timer still runs at 1 msec resolution. That means that our test is\n');
        fprintf('not diagnostic. Please exit and restart Matlab, then retry. Make sure that any\n');
        fprintf('kind of multi-media applications are closed while the test is running\n\n');
        error('Test failed, because could not disable 1khZ timer Interrupts.');
    end

    overrunticks = length(find(tickinc > avgtickdelta + 0.0001));
    fprintf('The low-res timer showed timing spikes a total of %i times.', overrunticks);

    if overrunticks > 0
        fprintf(' Maximum spike %f msecs.\n\n', max(tickinc)*1000);
        fprintf('The low-res timer seems to be a bit shaky. That will only be a problem if it shows large spikes often,\n');
        fprintf('because it could cause the online correctness tests to fail (false positive) and trigger an unwanted\n');
        fprintf('switch from the high-res timer to the unreliable low-res timer.\n');
        fprintf('Small occasional spikes are not uncommon on such a deficient system as Windows. Many spikes could\n');
        fprintf('indicate that some of your hardware devices or its drivers don''t operate optimally or the system\n');
        fprintf('is severly overloaded in some way.\n\n');
    else
        fprintf('\n\n');
    end

    figure;
    plot(ticks, raw);
    title('Block II: Elapsed time (msecs) of low-res (x-axis) vs. high-res (y-axis) timer:');

    figure;
    plot((raw - ticks)*1000);
    title('Block II: Lag of high-res timer vs. low-res timer, plotted against sample count:');
    avglag = mean(raw-ticks)*1000;
    lag1 = (raw(1) - ticks(1)) * 1000;
    lagn = (raw(n) - ticks(n)) * 1000;

    % Most important: The ratio between low-res and high-res timer:
    ratio = (raw(n)-raw(1))/(ticks(n)-ticks(1));
    fprintf('lag1 = %15.6f , lag%i = %15.6f - Not too important...\n\n', lag1, n, lagn);
    fprintf('Ratio of elapsed time in high-res vs. low-res clock: %15.10f\n', ratio);

    if abs(ratio - 1.0) > 0.05
        fprintf('\n\nTimers DISAGREE by more than 5 percent! --> This likely indicates an UNRELIABLE high-precision timer!!\n\n');
        softfail = 1;
    end

    if min(diff(raw)) < 0
        fprintf('TIME WAS RUNNING BACKWARDS!!! The TSC timers are not synchronized across cpu cores!\n');
        softfail = 1;
    end

    % Test 3: Let CPU go to sleep for a few seconds and see if there's drift:
    % We use pause() here, because WaitSecs would take counter-measures,
    % rendering our little test non-diagnostic.
    fprintf('\nGoing to sleep for 10 seconds, then checking if this idle time did something bad...\n');
    [cooked1 lowres1 highres1] = GetSecs(-1);
    pause(10);
    [cooked2 lowres2 highres2] = GetSecs(-1);
    fprintf('Elapsed time should be 10 seconds: lowres timer says %f secs, highres timer says %f secs.\n', lowres2-lowres1, highres2-highres1);
    deltaratio = abs(((lowres2-lowres1) / (highres2-highres1)) - 1.0);
    if deltaratio > 0.05
        fprintf('\n\nTimers DISAGREE by more than 5 percent! --> This likely indicates an UNRELIABLE high-precision timer (delta %f %%)!\n\n', deltaratio * 100);
        softfail = 1;
    end

    fprintf('\nTesting TSC sync across cpu cores on multi-core machines during enforced cpu hopping.\n');
    % Test 4: Multicore hopping stress test:
    timewarp = 0;
    told = GetSecs;
    for i=1:n
        tnew = GetSecs(mod(i,2)+1);
        if tnew < told
            fprintf('Time warps (time going backwards!) detected when switching between CPU cores!!! (delta = %f secs)\n', tnew - told);
            timewarp = timewarp + 1;
        end
        told = tnew;
    end

    if timewarp > 0
        fprintf('TIME WAS RUNNING BACKWARDS %i TIMES! Synchronization across cpus is broken on your setup despite PTB taking counter-measures!\n', timewarp);
        softfail = 1;
    end

    % Give a summary of behaviour without new PsychTimeGlue implementation:
    if softfail > 0
        fprintf('YOUR SYSTEMS TIMING IS BROKEN AND UNRELIABLE!!!\n');
        fprintf('Read in the online help of this test script and on the PTB Wiki about measures\n');
        fprintf('that may allow to fix your systems timing. Apply them, then retest.\n\n');
    end

    % Give a final conclusion:
    fprintf('\n\nPRELIMINARY CONCLUSION (Disclaimer: No test is perfect, this one isn''t an exception)\n\n');
    if hardfail > 0
        fprintf('Your system shows severe timing problems and PTB was unable to fix them\n');
        fprintf('with its built-in workarounds. Likely, PTB will detect this problems at runtime\n');
        fprintf('and switch to a low resolution backup time source, so it is at least useable for\n');
        fprintf('studies with low requirements on timing precision.\n');
        fprintf('However, if PTB should fail to detect the problems at runtime, you are in trouble.\n');
        fprintf('Read the help text of this script and on our Wiki how to solve this properly.\n\n');
    else
        fprintf('The test could not find indications of timing problems on your system when running\n');
        fprintf('with the current PTB. However, there is some likelihood of the test reporting false\n');
        fprintf('negatives.\n\n');
    end
    
    if hardfail == 0
        if softfail > 0
            fprintf('When running the timing tests without PTB''s automatic workarounds,\n');
            fprintf('your system showed severe timing problems!\n');
            fprintf('Studies run with older versions of PTB or other toolkits may have been\n');
            fprintf('affected by false time measurements. You may want to run some of your old\n');
            fprintf('scripts with the new PTB to see if you get different results.\n\n');
            fprintf('Read the help text of this script and on our Wiki for more background information.\n\n');
        else
            fprintf('So far your system behaved well, even without PTB''s fixes applied. This could mean it is well behaved.\n');
            fprintf('Unfortunately it could also mean that the test is unable to detect real problems (false negative) :(\n\n');
        end
    end

    Priority(0);
    fprintf('Test finished.\n');
end

if ~IsWin
    error('GetSecsTest is not yet capable of meaningful operation on non-Windows systems.');
end

return;
