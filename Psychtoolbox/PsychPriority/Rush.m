function Rush(rushedCode, priorityLevel)
% Rush(rushedCode, priorityLevel)
%
% Note: This function is not needed anymore. Use Priority() instead to
% simplify your life. The function is only left for backward compatibility
% to keep old code running.
%
% Rush runs a critical bit of your Matlab code with minimal interruption by
% other tasks. The first argument is a string containing Matlab code to be
% passed to EVAL. Within the string, you can have multiple statements
% separated by ";" or ",".
% 
% The optional "priorityLevel" argument specifies how much to block
% interrupt tasks. Use MaxPriority to determine the highest priority that
% allows normal operation of the functions you use, e.g. SND and SCREEN
% 'WaitBlanking'. We suggest you always call MaxPriority rather than hard
% coding any particular priorityLevel, so that your program will gracefully
% adapt to run optimally on any computer. Here's a typical use:
% 
% 	Screen('Screens');	% Make sure all functions (SCREEN.mex) are in memory.
% 	i=0;				% Allocate all variables.
% 	loop={
% 		'for i=1:100;'
% 			'Screen(window,''WaitBlanking'');'
% 			'Screen(''CopyWindow'',w(i),window);'
% 		'end;'
% 	};
% 	priorityLevel=MaxPriority(window,'WaitBlanking');
% 	Rush(loop,priorityLevel);
% 
% Allowable 'priorityLevel' settings are described in "help Priority".
%
% see also: Priority

% HISTORY
% 4/21/97   dgp Wrote it.
% 4/22/97   dgp Added alternate calling form.
% 4/23/97   dgp Updated.
% 4/30/97   dgp Updated to incorporate results of running RushTest.
% 5/6/97    dgp Explain how to choose a priorityLevel.
% 5/30/97   dgp Explain WaitForVBLInterrupt Preference.
% 6/22/97   dgp WaitBlanking calls cscSetEntries, not cscGetEntries.
% 2/1/98    dgp MaxPriority.
% 2/8/98    dgp add priorityLevel 0.5.
% 3/24/98   dgp eliminate priorityLevel -1.
% 7/16/98   dgp demonstrate use of cell array in example.
% 7/22/98   dgp Added warnings about loading functions and Backgrounding.
% 3/15/99   xmz Added comment for Windows version.
% 3/18/99   dgp Added warnings about VM. Document dynamic control of backgrounding.
% 2/4/00    dgp Updated for Mac OS 9.
% 12/16/03  awi Wrote OS X Documentation and merged with OS 9 + Windows Rush.m.
% 7/16/04   awi Cosmetic.  Restored HISTORY.
% 5/22/07   mk  Reimplemented Rush() in this script. Long-term we wanna get
%               rid of the mex files.
% 1/26/08   mk  Small fix: Use evalin('caller', ...) instead of eval(...),
%               so variables in the scope of calling M-Function are available
%               to Rushed code.

if nargin < 2
    error('You must specify both, the string with code to rush as well as the rush priority level!');
end

try
    Priority(priorityLevel);
    evalin('caller', rushedCode);
    Priority(0);
catch %#ok<CTCH>
    Priority(0);
    psychrethrow(psychlasterror);
end
