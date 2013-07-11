function oldPriority=Priority(newPriority)
% oldPriority=Priority([newPriority])
%
% Query and optionally set the execution priority of your script to
% 'newPriority' and switch from non-realtime to realtime scheduling (for
% 'newPriority' values greater than zero).
%
% Higher priority levels will reduce the chance of other running
% applications and system processes interfering with the execution timing
% of your script, or reduce the severity of interference. However, various
% factors influence the timing behaviour of your computer and its hardware,
% so while use of Priority() is one step to ensure good timing, it is not
% the only required step. The PsychDocumentation subfolder of your
% distribution as well as the Psychtoolbox Wiki contain various documents
% with tips on how to improve system timing behaviour. All else being
% equal, the Linux operating system should provide best possible timing
% behaviour out of the box, with various possibilities to tweak and tune the
% system for even better timing behaviour. Mac OSX is the second best
% choice for precise timing.
%
% OS X: ___________________________________________________________________
%
% Unlike on OS 9 and Windows, On OS X setting a high priority will not
% disable important parts of your computer such as the keybaord and hard
% drive. Priority calls the Psychtoolbox mex function MachSetPriority.  For
% simple, OS-neutral control of process priority on OS X use Priority. For
% complicated, fine-grained, OS-specific control of process priority use
% MachSetPriority.    
%
% At priority settings greater than 0, your script must limit MATLAB's  use
% of CPU time to stay within the limits specified by the priority value.
% The consequence of failing to restrict the CPU time to the limits
% specified by the priority values is that OS X will demote MATLAB back to
% priority 0. You can check to see if this has happened by calling
% Priority with no arguments.
%
% Within a script there are two ways to limit MATLAB's use of CPU time:
%
%   1. Call "WaitSecs" 
%   2. Call Screen('Flip',...)
% 
% Both calls will sleep the main MATLAB thread, surrendering CPU time to
% other threads on the system until the MATLAB thread is awakened.
% WaitSecs sleeps the MATLAB thread for a specified period.
% Screen('Flip',...) sleeps the MATLAB thread until the next video blanking
% interval.
%
% The priority value corresponds to the proportion of the video frame
% period  during which MATLAB is guranteed 100 percent of CPU time:
%
%   GuaranteedCPUtimePerFramePeriod = priority/10 * framePeriod
%
% If your computer has multiple video displays running at different frame
% rates then Priority will choose the shortest frame period.
%
% The consequence of a thread exceeding the specified limits on CPU usage
% are that the MACH kernel will revoke its real-time priority status.
% However, the threshold at which the MACH kernel decices to revoke
% priority is not known.  We have observed that at 100% CPU usage the
% Jaguar Mach kernel will revoke priority after 2.5 seconds, indepent of
% the CPU time requested when setting priority.  
%
% Note that if you change the video frame rate after setting priority then
% the framePeriod which Priority used when setting the priority level will
% no longer match the true frame period. Therefore you should change the
% priority level after changing the video mode and not before. 
%
%
% WINDOWS: ________________________________________________________________
% 
% For the Windows version of Priority (and Rush), the priority levels set
% are  "process priority levels". There are 3 priority levels available,
% levels 0, 1, and 2. Level 0 is "normal priority level", level 1 is "high
% priority level", and level 2 is "real time priority level". Combined with
% thread priority levels, they determine the absolute priority level of the
% Psychtoolbox threads. Threads are executed in a "round robin" fashion on
% Windows, with the  lower priority threads getting cpu time slice only
% when no higher priority thread is ready to execute. Currently, no tests
% have been done to see what tasks are pre-empted by setting the Matlab
% process to real-time priority. It does seem to block keyboard input,
% though, so for example if you have a clut animation going on at priority
% level 2, then the force-quit key combo (Ctrl-Alt-Delete) does not  work.
% However, the keyboard inputs are still sent to the message queue, so
% GetChar or GetClicks still work if they are also called at priority level
% 2.
%
% Typically you will not want to choose a higher priority than 1 unless you
% know exactly what you're doing.
% 
% LINUX: __________________________________________________________________
%
% To enable use of Priority(), you must run the script
% PsychLinuxConfiguration at least once and follow its instructions.
% Normally the script will get executed automatically by our installer if
% you got Psychtoolbox via DownloadPsychtoolbox or UpdatePsychtoolbox.
% However, if you got Psychtoolbox directly from the package manager of
% your Linux distribution or from the NeuroDebian repositories, it may be
% neccessary to execute the PsychLinuxConfiguration function manually.
%
% GNU/Linux supports priority levels between 0 and 99. Zero means standard
% non-realtime timesharing operation -- Play fair with all other
% applications. Increasing values greater than zero mean realtime mode: The
% Matlab/Octave process is scheduled in realtime mode with a priority
% corresponding to the given number - higher means better, but also more
% likelihood of interfering with system processes. Try to stick to a level
% of 1 unless you have reason to go higher. Level one is usually
% sufficient, unless you run other realtime applications on the same
% machine in parallel and want to prioritize Psychtoolbox relative to those
% applications.
%
% In realtime mode, PTB will also try to lock down all of Matlab/Octaves
% memory into physical RAM, preventing it from being paged out to disk by
% the virtual memory manager. If it works, it's great! However, the amount
% of lockable memory is restricted to max. 50% of installed RAM memory on
% some older Linux setups, so if Matlab/Octave/your experiment would need
% more than 50% of available system memory, this will fail. PTB will output
% an informative warning in this case, but continue otherwise unaffected.
% Realtime scheduling will be still effective, you'll just lose the bonus
% of memory locking.
%
% On old single processor computers, be careful not to create any
% uninterruptible infinite loops in your code when running realtime,
% otherwise your system may lock up, requiring a hard reboot!
% _________________________________________________________________________
%
% see also OS X:    Rush
% see also Windows: Rush
% see also OS 9:    Rush, ScreenTest, RushTest, LoopTest

% NOTES:
% priorityStruct looks like this.   
%
% threadID:
% flavor: 'THREAD_TIME_CONSTRAINT_POLICY'
% policy:
%     period: 
%     computation: 
%     constraint: 
%     preemptible: 
% policySize:
% policyFillSize:
% getDefault:
% isDefault:

% HISTORY:
% 4/13/97   dgp  Wrote it
% 4/22/97   dgp  Updated
% 5/1/97    dgp  Updated
% 5/31/97   dgp  Updated
% 6/3/97    dgp  Updated
% 3/15/99   xmz  Added comments for Windows version.
% 2/4/00    dgp  Updated for Mac OS 9.
% 6/13/02   awi  Cosmetic 
% 12/11/03  awi  Wrote documentation and scripting for OS X
% 7/16/03   awi  Merged together documetation for OS X and OS 9 and
%                Windows.
% 8/24/05   awi  Kill and start update under standard, not time constraint
%                priority.
%                When raising priority, issue error and return if we fail
%                to kill the update process.  Do not set time contstraint
%                priority.   
% 2/21/06   mk   Updated for Linux.
% 5/29/06   mk   Merged a fix for properly restarting update process.
%                Bug found & fix provided by Michael Ross.
% 8/02/06   mk   Bugfix by Michael Ross merged: oldPriority was sometimes
%                reported with some roundoff error.
% 8/02/06   mk   We now detect if we're running on a OS-X 10.4.7 system or
%                later. We do not kill and restart the system update process
%                anymore in that case as it isn't necessary anymore.
% 5/15/07   mk   Priority.dll on Windoze replaced by code in this M-File.
% 6/01/09   mk   Enable Priority() support for OS/X + Octave-3.
% 7/10/13   mk   Remove handling of update process, as we no longer support
%                OSX versions < 10.6, so no need to kill/restart updated.

if ~isempty(getenv('PSYCH_IN_VM'))
   % Running inside Virtual machine. We no-op. Everything else has potential
   % for bogging the VM down to a standstill:
   oldPriority = 0; % Dummy return - No realtime sched active.
   return;
end

if IsLinux
   % Linux: We do not use a separate MEX file anymore. Instead we use a
   % built-in helper subroutine of Screen(), accessed via the special code -5
   % when calling the 'GetMouseHelper' subfunction.
   if nargin < 1
      newPriority = [];
   end

   if isempty(newPriority)
      % No priority provided: Just return the current level:
      oldPriority = Screen('GetMouseHelper', -5);
   else
      % New priority provided: Query and return old level, set new one:
      if newPriority < 0 || newPriority > 99
         error('Invalid Priority level specified! Not one between 0 and 99.');
      end

      oldPriority = Screen('GetMouseHelper', -5, newPriority);
   end

   % Done for Linux:
   return;
end;

if IsOSX
    %Get the current settings.
    [flavorNameString, priorityStruct] = MachGetPriorityFlavor;
    if strcmp('THREAD_STANDARD_POLICY', flavorNameString)
        oldPriority=0;
    else % strcmp('THREAD_TIME_CONSTRAINT_POLICY', flavorNameString)
        %Values in priority struct returned by MachGetPriorityFlavor are in ticks but should be in seconds
        %it does not matter here because we are only concerned with the ratio
        if priorityStruct.policy.period == 0 || priorityStruct.policy.computation == 0
            %this is an illegitimate setting, so restore to standard priority
            %and go from there.
            MachSetStandardPriority;
            oldPriority=0;
        else
            oldPriorityRatio= priorityStruct.policy.computation / priorityStruct.policy.period;
            % Make sure we never exceed an oldPriority of 9, even in case
            % of roundoff-error. Bug found and fixed by Michael Ross.
            oldPriority=min(oldPriorityRatio * 10, 9);
        end
    end

    %if no new setting is given then return
    if nargin==0
        return
    end
    
    %bounds check the input argument
    if newPriority > 9
        error('"newPriority" value exceeds maximum allowable value of 9');
    elseif newPriority < 0
        error('"newPriority" value is less than the minimum allowable value of 0');
    end
    
    %if the priority level calls for time constraint priority...
    if newPriority > 0
        % Find the frame periods.  FrameRate returns the nominal frame rate, which
        % for LCD displays is 0.  We assume 60 in that case.
        defaultFrameRate=60;    %Hz
        screenNumbers=Screen('Screens');
        for i = 1:length(screenNumbers)
            frameRates(i)=Screen('FrameRate', screenNumbers(i)); %#ok<AGROW>
        end
        [zeroRates, zeroRateIndices]=find(frameRates==0); %#ok<*ASGLU>
        frameRates(zeroRateIndices)=defaultFrameRate;
        framePeriods=1./frameRates;
        
        % MachSetTimeConstraintPriority(periodSecs,computationSecs, constraintSecs, preemptibleFlag) 
        periodSecs=min(framePeriods);
        computationSecs = newPriority/10 * periodSecs;
        constraintSecs=computationSecs;
        preemptibleFlag=1;
        
        
        %set the prioirity
        MachSetTimeConstraintPriority(periodSecs,computationSecs, constraintSecs, preemptibleFlag);
        
        %if the priority level calls for standard priority then ..
    else    %priority==0
        %restore standard priority
        MachSetStandardPriority;
    end
end

% Microsoft Windows?
if IsWin
   % Yes. We do not use a separate MEX file anymore. Instead we use a
   % built-in helper subroutine of Screen(), accessed via the special code -3
   % when calling the 'GetMouseHelper' subfunction.
   if nargin < 1
      newPriority = [];
   end
   
   if isempty(newPriority)
      % No priority provided: Just return the current level:
      oldPriority = Screen('GetMouseHelper', -3);
   else
      % New priority provided: Query and return old level, set new one:
      if newPriority<0 || newPriority>2
         error('Invalid Priority level specified! Not one of 0, 1 or 2.');
      end
      
      oldPriority = Screen('GetMouseHelper', -3, newPriority);
   end

   % Done for Windows...
end;
