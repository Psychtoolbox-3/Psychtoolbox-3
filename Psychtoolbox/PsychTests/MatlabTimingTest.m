% MatlabTimingTest
% 
% Assign the main MATLAB thread "time constraint" priority status and 
% run a tight loop which records with high preceision the time at every 
% pass through the loop.  Calculate the time intervals between succussive
% passes through the loop.  
%
% Granting MATLAB "time constraint" priority prevents any other thread on the
% system from preempting the main MATLAB thread. Other threads recieve
% CPU time only at the discretion of the main MATLAB thread which executes
% your MATLAB scripts and functions.  Therefore, the delay between
% the nth and the (n+1)th timing loop is not caused by other threads
% preempting MATLAB, but instead by some activity of MATLAB itself or the script which it executes.
%
% MatlabTimingTest samples the time using GetSecsMex, a custom mex file
% which calls the native OSX Core Audio function AudioGetCurrentHostTime()
% The precision depends on the CPU and clocks speed, on a 1GHZ G4 clock tick 
% period is 30 nanoseconds)
%
% AUTHORS: Allen Ingling
% SEE ALSO: GetSecsMex

% HISTORY:
% 04/09/03      awi Wrote it. 
% 04/14/03      awi rewrote it after changing Priority to use mach threads.
% 8/13/03       awi Changed name to "TestMATLABTimingOSX.m and re-wrote for
%                   timing demo package mex functions.  Clean up
%                   documentaion and explained the script in detail for the
%                   benefit of Mathworks.
% 9/24/03       awi Fixed a dependencies on obsolete files,
%                   detected when packaging the timing tests for release.
% 11/04/03      awi Added axis labels. 
% 1/29/05       dgp Cosmetic. Changed "Seconds" to "Secs".


% Setting "time constraint" priority settings and blocking duration:
%
% periodSecs= 1/1000;                                  
% computationSecs = periodSecs/10;
% constraintSecs= computationSecs;              
% preemptibleFlag=1;
% blockingIntervalSecs=periodSecs;
%
% Out of every millisecond, we guarantee the MATLAB process up to 100 microseconds
% of CPU time.  That is, MATLAB will lay claim to up to 1/10 of the total CPU time.
% Gurantee that MATLAB's 100 microseconds of CPU time fall within a
% 100 microsecond window, from the start of computation to the end of
% computation.  Within that window, allow MATLAB to be interrupted by other
% threads.
%
% Note that we do not know in advance of the timing loop how much CPU time we will
% actually need, that depends on how much CPU time MATLAB expends 
% in processing a scripted loop pass. We can however retroactively validate our
% choice of values by examing timing results:
% 
%  blockingInterval + actualMATLABCPUUsage = loopPassDuration
% 
% therefore:
%
%   loopPassDuration - blockingInterval = actualMATLABCPUUsage
%
% we know everything on the left hand side can find actualMATLABCPUUsage.  If  we have chosen 
% parameters correctly then:
%
% actualMATLABCPUUsage <= computationSecs
% 
% and, 
%
% actualMATLABCPUUsage/(actualMATLABCPUUsage + blockingInterval) < computationSecs/periodSecs)
% 
% unknown quantity actualMATLABCPUUsage is found:
%
% loopPassDuration = actualMATLABCPUUsage + blockingInterval
% actualMATLABCPUUsage = loopPassDuration - blockingInterval
% 
% substituting into the inequality we get:
%
% (loopPassDuration - blockingInterval)/loopPassDuration < computationSecs/periodSecs)
%
% Note that the above is sufficient when loopPassDuration < periodSecs.
% If loopPassDuration > periodSecs then we must also have
% actualMATLABCPUUsage <= computationSecs.

clear all;  

% designate an index for the block of test parameters to be used. 
tci = 2;                    %tci stands for  test condition index.  
fixedAxisY=0;
fixedAxisValueY=0.020;

% Test condition 1 shows delays at 30-second intervals.
%
%
% designate an index for this block of test parameters. 
tcb= 1;
%set flags which specify various test conditions
tc(tcb).enableTCPriority=1;
tc(tcb).useFlip=0;              %if set use Screen('flip') command instead of instead of BlockSecs command to block during timing loop.
tc(tcb).useMexOnly=1;        % GetSecsMex is a bare-bones mex file. The alternative, GetSecs, compiles in an extra abstraction layer.  
%specify the 'time constraint' parameters
tc(tcb).periodSecs= 1/1000;                                  
tc(tcb).computationSecs = tc(tcb).periodSecs/10 * 5 ;
tc(tcb).constraintSecs= tc(tcb).computationSecs;              
tc(tcb).preemptibleFlag=0;
tc(tcb).blockingIntervalSecs=tc(tcb).periodSecs - tc(tcb).computationSecs;
% specify for how many seconds we want to run out timing loop
tc(tcb).testDurationSecs = 2*60+1;
%specify a factor by which to mulitply the memory allocation estimate, for
%saftey margin.
tc(tcb).allocationMarginFactor=2.0;
tc(tcb).maxDisplayExcessComputationDurations=100;
% specify values describing period glithes.  As part of the analysis after
% the timing test we search for regular glitches described by these
% paramters.
tc(tcb).interGlitchIntevalSecs = 30;
tc(tcb).interGlitchIntervalJitter = 0.1;        % tolerance of interGlitchIntevalSecs value.

% Test condition 2 shows glitches at 60 second intervals which increase in
% magnitude over the duration of the test period.  
%
% designate an index for this block of test parameters. 
tcb= 2;
%set flags which specify various test conditions
tc(tcb).enableTCPriority=1;
tc(tcb).useFlip=1;      %if set use Screen('flip') command instead of instead of BlockSecs command to block during timing loop.
tc(tcb).useMexOnly=1;        % GetSecsMex is a bare-bones mex file. The alternative, GetSecs, compiles in an extra abstraction layer.  
%specify the 'time constraint' parameters
tc(tcb).periodSecs= 1/60;                                  
tc(tcb).computationSecs = tc(tcb).periodSecs/10 * 4;
tc(tcb).constraintSecs= tc(tcb).computationSecs;              
tc(tcb).preemptibleFlag=0;
tc(tcb).blockingIntervalSecs=tc(tcb).periodSecs - tc(tcb).computationSecs;
% specify for how many seconds we want to run out timing loop
tc(tcb).testDurationSecs = 4*60+1;
%specify a factor by which to mulitply the memory allocation estimate, for
%saftey margin.
tc(tcb).allocationMarginFactor=1.25;
tc(tcb).maxDisplayExcessComputationDurations=100;
% specify values describing period glithes.  As part of the analysis after
% the timing test we search for regular glitches described by these
% paramters.
tc(tcb).interGlitchIntevalSecs = 30;
tc(tcb).interGlitchIntervalJitter = 0.1;        % tolerance of interGlitchIntevalSecs value.
 
% Test condition 3 shows ?
%
% designate an index for this block of test parameters. 
tcb= 3;
%set flags which specify various test conditions
tc(tcb).enableTCPriority=1;
tc(tcb).useFlip=0;      %if set use Screen('flip') command instead of instead of BlockSecs command to block during timing loop.
tc(tcb).useMexOnly=1;        % GetSecsMex is a bare-bones mex file. The alternative, GetSecs, compiles in an extra abstraction layer.  
%specify the 'time constraint' parameters
tc(tcb).periodSecs= 1/100;                                  
tc(tcb).computationSecs = tc(tcb).periodSecs/10 * 1 ;
tc(tcb).constraintSecs= tc(tcb).computationSecs;              
tc(tcb).preemptibleFlag=0;
tc(tcb).blockingIntervalSecs=tc(tcb).periodSecs - tc(tcb).computationSecs;
% specify for how many seconds we want to run out timing loop
tc(tcb).testDurationSecs = 2*60+1;
%specify a factor by which to mulitply the memory allocation estimate, for
%saftey margin.
tc(tcb).allocationMarginFactor=1.25;
tc(tcb).maxDisplayExcessComputationDurations=100;
% specify values describing period glithes.  As part of the analysis after
% the timing test we search for regular glitches described by these
% parameters.
tc(tcb).interGlitchIntevalSecs = 30;
tc(tcb).interGlitchIntervalJitter = 0.1;        % tolerance of interGlitchIntevalSecs value.

% Test condition 4 shows ?
%
% designate an index for this block of test parameters. 
tcb= 4;
%set flags which specify various test conditions
tc(tcb).enableTCPriority=1;
tc(tcb).useFlip=0;      %if set use Screen('flip') command instead of instead of BlockSecs command to block during timing loop.
tc(tcb).useMexOnly=1;        % GetSecsMex is a bare-bones mex file. The alternative, GetSecs, compiles in an extra abstraction layer.  
%specify the 'time constraint' parameters
tc(tcb).periodSecs= 1/100;                                  
tc(tcb).computationSecs = tc(tcb).periodSecs/10 * 3 ;
tc(tcb).constraintSecs= tc(tcb).computationSecs;              
tc(tcb).preemptibleFlag=0;
tc(tcb).blockingIntervalSecs=tc(tcb).periodSecs - tc(tcb).computationSecs;
% specify for how many seconds we want to run out timing loop
tc(tcb).testDurationSecs = 2*60+1;
%specify a factor by which to mulitply the memory allocation estimate, for
%saftey margin.
tc(tcb).allocationMarginFactor=1.25;
tc(tcb).maxDisplayExcessComputationDurations=100;
% specify values describing period glithes.  As part of the analysis after
% the timing test we search for regular glitches described by these
% parameters.
tc(tcb).interGlitchIntevalSecs = 30;
tc(tcb).interGlitchIntervalJitter = 0.1;        % tolerance of interGlitchIntevalSecs value.


% Test condition 5 is like 1 except for without TC priority. 
%
%
% designate an index for this block of test parameters. 
tcb= 5;
%set flags which specify various test conditions
tc(tcb).enableTCPriority=0;
tc(tcb).useFlip=0;              %if set use Screen('flip') command instead of instead of BlockSecs command to block during timing loop.
tc(tcb).useMexOnly=1;        % GetSecsMex is a bare-bones mex file. The alternative, GetSecs, compiles in an extra abstraction layer.  
%specify the 'time constraint' parameters
tc(tcb).periodSecs= 1/1000;                                  
tc(tcb).computationSecs = tc(tcb).periodSecs/10 * 5 ;
tc(tcb).constraintSecs= tc(tcb).computationSecs;              
tc(tcb).preemptibleFlag=0;
tc(tcb).blockingIntervalSecs=tc(tcb).periodSecs - tc(tcb).computationSecs;
% specify for how many seconds we want to run out timing loop
tc(tcb).testDurationSecs = 2*60+1;
%specify a factor by which to mulitply the memory allocation estimate, for
%saftey margin.
tc(tcb).allocationMarginFactor=1.25;
tc(tcb).maxDisplayExcessComputationDurations=100;
% specify values describing period glithes.  As part of the analysis after
% the timing test we search for regular glitches described by these
% paramters.
tc(tcb).interGlitchIntevalSecs = 30;
tc(tcb).interGlitchIntervalJitter = 0.1;        % tolerance of interGlitchIntevalSecs value.


% Test condition 6 was written for use with Thread Viewer. 
% We slow down the sample rate 
%
%
% designate an index for this block of test parameters. 
tcb= 6;
%set flags which specify various test conditions
tc(tcb).enableTCPriority=1;
tc(tcb).useFlip=0;              %if set use Screen('flip') command instead of instead of BlockSecs command to block during timing loop.
tc(tcb).useMexOnly=1;        % GetSecsMex is a bare-bones mex file. The alternative, GetSecs, compiles in an extra abstraction layer.  
%specify the 'time constraint' parameters
tc(tcb).periodSecs= 0.010;                                  
tc(tcb).computationSecs = 0.001;
tc(tcb).constraintSecs= 0.001;              
tc(tcb).preemptibleFlag=0;
tc(tcb).blockingIntervalSecs= 0.009;
% specify for how many seconds we want to run out timing loop
tc(tcb).testDurationSecs = 60;
%specify a factor by which to mulitply the memory allocation estimate, for
%saftey margin.
tc(tcb).allocationMarginFactor=2.0;
tc(tcb).maxDisplayExcessComputationDurations=100;
% specify values describing period glithes.  As part of the analysis after
% the timing test we search for regular glitches described by these
% paramters.
tc(tcb).interGlitchIntevalSecs = 30;
tc(tcb).interGlitchIntervalJitter = 0.1;        % tolerance of interGlitchIntevalSecs value.


clear tcb;      %prevent  this guy from sneeking into the scriptage below.
% _________________________________________________________________________


%if we are using screen buffer flips to block during the timing loop
%instead of BlockSecs then open a window
if tc(tci).useFlip
    sNumber=max(Screen('Screens'));
    w=Screen('OpenWindow', sNumber, [],[],[], 2);
end


%we should pre-allocate the vector used to store time samples so to 
% avoid delays during the loop caused by memory allocation and garbage
% collection.  Estimate how many elements we need to allocate based on 
% a short trial loop
preTrialLoopDurationSecs=5;
safteyFactor=1.5;     
fprintf('Running %d second pre-trial loop...\n', preTrialLoopDurationSecs); 
preTrialNumLoops=0;
if tc(tci).enableTCPriority
    MachSetTimeConstraintPriority(tc(tci).periodSecs, tc(tci).computationSecs, tc(tci).constraintSecs, tc(tci).preemptibleFlag);
end
fprintf(['Priority flavor set to ' MachGetPriorityFlavor '\n']);
t=GetSecsMex;
tStart=GetSecsMex;
endTime=tStart+preTrialLoopDurationSecs;
t=tStart;
if ~ tc(tci).useFlip
    if tc(tci).useMexOnly;       
	    while t<endTime
            t=GetSecsMex;
            preTrialNumLoops= preTrialNumLoops+1;
            SleepSecsMex(tc(tci).blockingIntervalSecs);   
        end
    else
	    while t<endTime
            t=GetSecs;
            preTrialNumLoops= preTrialNumLoops+1;
            SleepSecsMex(tc(tci).blockingIntervalSecs);   
        end
    end 
else
    if tc(tci).useMexOnly;       
        while t<endTime
            t=GetSecsMex;
            preTrialNumLoops= preTrialNumLoops+1;
            Screen('Flip',w);   
        end
    else 
        while t<endTime
            t=GetSecs;
            preTrialNumLoops= preTrialNumLoops+1;
            Screen('Flip',w);   
        end
    end
end
   
preTrialMeasuredDurationSecs=GetSecsMex-tStart;
if tc(tci).enableTCPriority
    MachSetStandardPriority;
end
fprintf(['Priority flavor set to ' MachGetPriorityFlavor '\n']);
fprintf('pre-trial loop complete\n\n'); 
numEstimatedTimingLoops=round((tc(tci).testDurationSecs/preTrialMeasuredDurationSecs) * preTrialNumLoops * tc(tci).allocationMarginFactor);   
tVec=1:numEstimatedTimingLoops;

%run the timing loop
fprintf('Running %d second timing loop\n', tc(tci).testDurationSecs);
if tc(tci).enableTCPriority
    MachSetTimeConstraintPriority(tc(tci).periodSecs, tc(tci).computationSecs, tc(tci).constraintSecs, tc(tci).preemptibleFlag);
end
fprintf(['Priority flavor set to ' MachGetPriorityFlavor '\n']);
i=1;
tVec(1)=GetSecsMex;
endTime=tVec(1)+tc(tci).testDurationSecs;
if ~ tc(tci).useFlip
    if tc(tci).useMexOnly;       
		while tVec(i)<endTime
            i=i+1;
            tVec(i)=GetSecsMex; 
            SleepSecsMex(tc(tci).blockingIntervalSecs);
		end
    else
		while tVec(i)<endTime
            i=i+1;
            tVec(i)=GetSecs; 
            SleepSecsMex(tc(tci).blockingIntervalSecs);
		end
    end
else
    if tc(tci).useMexOnly;       
		while tVec(i)<endTime
            i=i+1;
            tVec(i)=GetSecsMex; 
            Screen('Flip',w);   
		end
    else
		while tVec(i)<endTime
            i=i+1;
            tVec(i)=GetSecs 
            Screen('Flip',w);   
		end
    end
end

if tc(tci).enableTCPriority
    MachSetStandardPriority;
end
fprintf(['Priority flavor set to ' MachGetPriorityFlavor '\n']);
numTrialLoops=i;
trialMeasuredDuration=tVec(i)-tVec(2);
fprintf('Timing loop complete\n\n'); 

if numEstimatedTimingLoops>=numTrialLoops
    fprintf('The actual number of timing loops was %d loops in %d seconds.\n', numTrialLoops, trialMeasuredDuration);
    fprintf('We preallocated enough memory for %d timing loops, a saftey factor of %f.\n', numEstimatedTimingLoops, numEstimatedTimingLoops/numTrialLoops);
else
    fprintf('WARNING:\n');
    fprintf('\tfailed to preallocate sufficient memory for timing results.  Try increasing the variable "marginFactor" and try again\n');
    fprintf('The actual number of timing loops was %d loops in %d seconds.\n', numTrialLoops, trialMeasuredDuration);
    fprintf('We preallocated enough memory for %d timing loops, a saftey factor of %f\n.', numEstimatedTimingLoops, numEstimatedTimingLoops/numTrialLoops);
end
tVec=tVec(2:i); 
tDiffVec=diff(tVec);
tDiffVecSampleTimes=tVec(1:end-1)-tVec(1);
plot(tDiffVecSampleTimes, tDiffVec, 'b');
xlabel('loop start time (seconds)');
ylabel('loop pass time (seconds)');
if fixedAxisY
    axis([-5, max(tDiffVecSampleTimes)+5, 0, fixedAxisValueY]);
else
    axis([-5, max(tDiffVecSampleTimes)+5, 0, 1.25 * max(tDiffVec)]);
end
minLoopDuration=min(tDiffVec);
maxLoopDuration=max(tDiffVec);
medianLoopDuration=median(tDiffVec);
fprintf('The shortest loop was: %f seconds\n', minLoopDuration);
fprintf('The longest loop was: %f seconds\n', maxLoopDuration);
fprintf('The median loop was: %f seconds\n', medianLoopDuration);

%plot horizontal lines across the graph marking "time constraint"
%parameters "period" and the blocking interval.
hold on;
periodSecsLineX=[0 tc(tci).testDurationSecs];
periodSecsLineY=[tc(tci).periodSecs tc(tci).periodSecs ];
plot(periodSecsLineX, periodSecsLineY, 'g');
blockingSecsLineX=periodSecsLineX;
blockingSecsLineY=[tc(tci).blockingIntervalSecs tc(tci).blockingIntervalSecs ];
plot(blockingSecsLineX, blockingSecsLineY, 'r');

% Find the n longest glitches.  n is set with the assumption that
% the total number of glitches is predicted by their falling
% at 30-second intervals.  If that assumption is wrong then 
% we will not find all of the long glitches
numExpectedGlitches=floor(tc(tci).testDurationSecs/tc(tci).interGlitchIntevalSecs);
[sortedtDiffVec, sortedtDiffVecIndices]=sort(tDiffVec);
sortedtDiffVecRev=fliplr(sortedtDiffVec);
sortedtDiffVecIndicesRev=fliplr(sortedtDiffVecIndices);

nLongestGlitches= sortedtDiffVecRev(1:numExpectedGlitches);
nLongestGlitchesIndices=sortedtDiffVecIndicesRev(1:numExpectedGlitches);
nLongestGlitchesTimestamps=tDiffVecSampleTimes(nLongestGlitchesIndices);
[nLongestGlitchesTimestampsSequenced, nLongestGlitchesTimestampsSequencedIndices]=sort(nLongestGlitchesTimestamps);
interGlitchIntervals=diff(nLongestGlitchesTimestampsSequenced);
for i=1:numExpectedGlitches
    plot(nLongestGlitchesTimestamps(i), nLongestGlitches(i), 'rx');
end
hold off;

fprintf([int2str(numExpectedGlitches) ' glitches predicted in ' num2str(tc(tci).testDurationSecs) ' second test interval, assuming 30-second interval between glitches\n']);
fprintf(['The ' int2str(numExpectedGlitches) ' longest loop delays occured at times and intervals:\n']);
for i=1:numExpectedGlitches
    fprintf([ '\t' num2str(nLongestGlitchesTimestampsSequenced(i)) ' s\n']);
    if  i<numExpectedGlitches
        fprintf([ '\t\tdelta=' num2str(interGlitchIntervals(i)) ' s\n']);
    end
end

%identify those glitches among the longest 4 which fall at 30-second
%intervals.  First find the intervals of 30-second duration, then find
%the samples which mark the end points of those intervals.
synchedIntervalDoubleIndices=find( (interGlitchIntervals < (tc(tci).interGlitchIntevalSecs + tc(tci).interGlitchIntervalJitter)) & (interGlitchIntervals > (tc(tci).interGlitchIntevalSecs - tc(tci).interGlitchIntervalJitter)));
if ~isempty(synchedIntervalDoubleIndices) 
    synchedIntervals=interGlitchIntervals(synchedIntervalDoubleIndices);
    intervalsStartGlitchesTimestamps=nLongestGlitchesTimestampsSequenced(1:end-1);
    intervalsEndGlitchesTimestamps=nLongestGlitchesTimestampsSequenced(2:end);
    synchedGlitchesTimestamps=unique([intervalsStartGlitchesTimestamps(synchedIntervalDoubleIndices) intervalsEndGlitchesTimestamps(synchedIntervalDoubleIndices)]);
else
    synchedGlitchesTimestamps=[];
end
 
    


% find the next longest glitch.  
nextLongestDelay=sortedtDiffVecRev(numExpectedGlitches+1); 
nextLongestDelayIndex=sortedtDiffVecIndicesRev(numExpectedGlitches+1);
nextLongestDelayTimestamp=tDiffVecSampleTimes(nextLongestDelayIndex);

fprintf(['The ' int2str(numExpectedGlitches+1) 'th longest loop delay, ']);
fprintf(['at time ' num2str(nextLongestDelayTimestamp) ' seconds, was ' num2str(nextLongestDelay) ' seconds, ']);
fprintf([num2str(nLongestGlitches(numExpectedGlitches)/nextLongestDelay) ' times smaller than the ' int2str(numExpectedGlitches) 'th longest delay\n']);

% Check to see that we always blocked for the specified period. 
underBlockIndices=find(tDiffVec < tc(tci).blockingIntervalSecs);
underBlocks=tDiffVec(underBlockIndices);
if ~isempty(underBlocks)
    fprintf('Detected loop durations shorter than blocking period.  The durations and times are:\n');
    for i=1:length(underBlocks);
        fprintf(['\t' num2str(underBlocks(i)) '\t' num2str(tDiffVecSampleTimes(underBlockIndices(i))) '\n']);
    end
end

% Check to see when we exceeded the specified "computation" CPU time
% allowance specified when assigning 'time constraint' priority. In the
% case that we use Screen('Flip') to block then the blocking interval is
% the video frame period which we do not know. So when using flip  we
% assume the blocking interval to be the median of the timing loop
% durations, which is actually an overestimate.  However, the error should
% be small in comparison to the  blocking interval itself because flip is
% fast.  
if tc(tci).useFlip
    measuredComputationSecs=tDiffVec-medianLoopDuration; 
else
    measuredComputationSecs=tDiffVec-tc(tci).blockingIntervalSecs;
end
excessComputationSecsIndices=find(measuredComputationSecs > tc(tci).computationSecs);
excessComputationSecs=measuredComputationSecs(excessComputationSecsIndices);
excessComputationTotalSecs=tDiffVec(excessComputationSecsIndices);
exessComputationSecsTimes=tDiffVecSampleTimes(excessComputationSecsIndices);
numExcessComputationLoops=length(excessComputationSecs);
fprintf(['The timing loop computation time exceeded the allocated computation time on ' int2str(numExcessComputationLoops) ' loops.\n']);
if numExcessComputationLoops > 0   
	fprintf('The loop durations, computation times, and timestamps are listed below.\n');
    if numExcessComputationLoops > tc(tci).maxDisplayExcessComputationDurations 
        fprintf(['(only the first ' int2str(tc(tci).maxDisplayExcessComputationDurations) ' are dislayed)\n']);
    end
	for i = 1:min([numExcessComputationLoops, tc(tci).maxDisplayExcessComputationDurations])
        fprintf('\t%1.5f\t%1.5f\t%2.5f\n', excessComputationTotalSecs(i),  excessComputationSecs(i), exessComputationSecsTimes(i) );
	end
end

% Check for compliance with:
% (loopPassDuration - blockingInterval)/loopPassDuration < computationSecs/periodSecs)
% For more info see notes at top of the file about this.
% For the tests which use flip instead of BlockSecs we do not actually
% know the blocking interval. This section should be modified accordingly. 
computationRatioLimit = tc(tci).computationSecs / tc(tci).periodSecs;
computationRatios = (tDiffVec - tc(tci).blockingIntervalSecs) / medianLoopDuration;  
computationRatioViolationIndices=find(computationRatios > computationRatioLimit); 
computationRatioViolationDurations=tDiffVec(computationRatioViolationIndices);
computationRatioViolationTimestamps=tDiffVecSampleTimes(computationRatioViolationIndices);
numComputationRatioViolations=length(computationRatioViolationIndices);

% Check to see if we exceed the "computation" CPU time allowance before the
% first glitch.  One possible explanation for the glitches is that the 
% Mach Kernel revokes realtime status at 30-second intervals as penalty 
% for the main MATLAB thread exceeding "computation" allowance.  However,
% if the thread has not exceeded that allowance before the first
% synchronized glitch, this strongly suggests that in fact the source
% of the glitch is MATLAB.
if ~isempty(synchedGlitchesTimestamps)
    earlyComputationViolations=exessComputationSecsTimes(exessComputationSecsTimes < min(synchedGlitchesTimestamps));
    if isempty(earlyComputationViolations)
        fprintf('The MATLAB thread did not exceed the CPU computation allowance before the first synchronized glitch.\n');
    else
        fprintf('The MATLAB thread exceeded the CPU computation allowance before the first synchronized glitch.\n');
        fprintf(['The violations occured at times: ' num2str(earlyComputationViolations) '\n']);
    end
end


if tc(tci).useFlip
    Screen('CloseAll');
end

% save results to files, including the plot.  We append to the existing
% file of the same name to avoid overwriting previous test results.

% identify which of the delays are those which fall at 30-second intervals,
% count them to see if we have fallen short.  
% Check wich of those corrspond to excess computations, list those which do
% and those which do not, and state whether we violated excess computation
% before the delays

% Check the theory that we only get under times for blocking if we allow
% preemption
