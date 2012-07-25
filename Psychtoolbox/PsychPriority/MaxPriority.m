function priorityLevel=MaxPriority(varargin);
% priorityLevel=MaxPriority([windowPtrOrScreenNumber],['WaitBlanking'],['PeekBlanking'],...
% 							['BlankingInterrupt'],['SetClut'],['ClutMovie'],...
% 							['SND'],['sound'],['speak'],...
% 							['GetSecs'],['WaitSecs'],['cputime'],...
% 							['KbCheck'],['KbWait'],['CharAvail'],['GetChar'],...
% 							['EventAvail'],['GetClicks'],['GetMouse'],['GetTicks'])
%
% MaxPriority.m receives a list of one or more function names, in any
% order, and returns the maximum priorityLevel that will allow all the
% named functions to work normally on this computer. Use MaxPriority
% before calling RUSH, to select the highest priorityLevel that's
% compatible with all the functions that you're rushing.
% 
% The name matching ignores case.
%
% OS X: ___________________________________________________________________
%
% On OS X all priority levels are safe for all functions. MaxPriority
% always returns 9, the highest priority level.  
%
% To preserve compatibility with other platforms we recommend using
% MaxPriority in your script on OS X, instead of the constant 9.
%
% WINDOWS: ________________________________________________________________
%
% Priority levels returned by MaxPriority are 0, 1 and 2.
% Although use of priority levels > 1 is possible and allowed by MaxPriority
% if you don't try to acquire input from keyboard or mouse, it is discouraged
% to use levels > 1 as this can interfere with execution of a lot of
% important system processes and severely reduce the stability of
% Windows execution.
%
% LINUX: __________________________________________________________________
%
% MaxPriority always returns 1, although levels of up to 99 are possible.
% We recommend to sticking to the lowest level, unless some tweaking for a
% specific setup or situation is necessary.
%
% _________________________________________________________________________
%
% See RUSH, Priority, MovieDemo, ScreenTest, SCREEN Preference MaxPriorityForBlankingInterrupt,
% SCREEN Preference SetClutPunchesBlankingClock.

% HISTORY
% 2/1/98 dgp wrote it.
% 2/8/98 dgp Added priorityLevel 0.5.
% 2/28/98 dgp Mention GetSecs feature of PeekBlanking.
% 2/28/98 dgp Reduced priority for WaitBlanking and PeekBlanking from 1 to 0.5
% 3/13/98 dgp Use the new SCREEN preference MaxPriorityForBlankingInterrupt.
% 3/24/98 dgp Use SCREEN preference WaitBlankingAlwaysCallsSetClut.
% 3/24/98 dgp Add 'SetClut'. Enhance 'SND' case to check for 6100.
% 7/25/98 dgp Updated for new SCREEN Preference names.
% 8/1/98 dgp Make smarter analysis of PeekBlanking.
% 3/18/99 dgp Check for VM and RAM Doubler.
% 3/23/99 dgp Cosmetic.
% 2/4/00 dgp Updated GetSecs test to directly test for UpTime trap instead of testing for PCI.
% 3/13/01 dgp Updated to accept struct from SCREEN('Computer'), as suggested by Harriet Allen.
% 2/12/02 dgp Fixed logic for 'WaitBlanking' to allow high priority if SetClutDriverWaitsForBlanking. 
%             Eliminated unnecessary variables "analyzePeekBlanking" and "peekBlankingNeedsInterrupt".
% 4/3/02  awi Added section for Windows.  Priorities are based on the comment at the bottom of the 
%             priority help file.
% 4/9/02 dgp Give VM warning only once.
% 6/20/02  awi Decrease Win priorities for GetMouse and GetClicks from to 2 to 1. 
% 6/21/02 dgp Deleted obsolete ClutMovie option. Use new MaximumSetClutPriority pref.
% 6/22/02 dgp On second thought, DON'T use the new MaximumSetClutPriority pref. That's already
%             taken care of by the C code surrounding the SetClut call, so we don't need
% 			  to worry about it here. The rest of the MATLAB loop can usefully run at
% 			  higher priority.
% 7/29/02 awi fixed comment of 6/20/02. 
% 7/17/04 awi Added OS X condition.  Partitioned help platformwise.
% 3/12/05 dgp Changed "strcmp" to "streq".
% 10/10/05 awi Noted changes by dgp on 3/12/05 
% 12 /31/05 mk Bugfix for windows part. Never worked when passing in a windowPtrOrScreenNumber
%              instead of keywords. Even if passing in keywords, the priority matching gave
%              wrong results (possibly too high priority) when multiple keywords were supplied.

% NOTES
% 7/17/04 awi
% The OS X section does not do thorough validation of arguments.  It does do:  
%  ¥If the window pointer or screen number argument is present it checks that it is
%  valid.  
%  ¥If a command is included, it checks that it is a valid command
%  for MacPriority.  
%  ¥It checks that at least one argument was passed.
%
% It does not check other requirements, for example some command arguments
% require the screen number or window pointer.
%
% 7/17/04 awi
% If we can figure out which MATLAB routines cause prolonged kernel locks
% then we could have MaxPriority on OS X return 0 for those.  See Mario
% Kleiner's list of Posix functions which cause Kernel locking and ask
% Scott French when MATLAB calls those.  Sigsetjmp is embedded in MATLAB
% run loop and that just has to go.  But other POSIX functions called only
% from some particular MATLAB functions could be avoided by using
% MaxPriority.  
%
% Also, since we can uncover the name of a function from within itself, we
% could improve MaxPriority by having it digest the script from which it is
% called, filter out comments, and search for functions which restrict
% priority level.  Thats not going to to be foolproof because quoted text
% is ambiguous, may be rushed or not.  The problem remains that you have to
% remember what in your script you have to pass to MaxPriority.  How about
% MaxPriorityScan. Returns keywords in script which are candidates for
% MaxPriority. Can be used as a programming tool from the command line, or
% embedded in a script and the results passed right into MacPriority.  
% 

%_______________________________________________________________________________________________________________________________
%The Windows Part

if IsWin
	if nargin<1
		error(['Usage: priorityLevel=MaxPriority([windowPtrOrScreenNumber],[''WaitBlanking''],[''PeekBlanking''],...' ...
		char(13) '                        [''BlankingInterrupt''],[''SetClut''],...'...
		char(13) '                        [''SND''],[''sound''],[''speak''],...'...
		char(13) '                        [''GetSecs''],[''WaitSecs''],[''cputime''],...'...
		char(13) '                        [''KbCheck''],[''KbWait''],[''CharAvail''],[''GetChar''],...'...
		char(13) '                        [''EventAvail''],[''GetClicks''],[''GetMouse''],[''GetTicks''])']);
	end
   
   match=0;
   priorityLevel = 2;    
	for i=1:nargin
		if ~ischar(varargin{i}) & ~isnumeric(varargin{i})
			error([ 'argument ' num2str(i) ' is of wrong type']);
		end
		if ischar(varargin{i})
			name=lower(varargin{i});
            %****cases which change the priority****
			if streq(name,lower('CharAvail'))
				priorityLevel=min(priorityLevel,1);
				match=1;
			end
            if streq(name,lower('GetChar'))
				priorityLevel=min(priorityLevel,1);
				match=1;
			end
			if streq(name,lower('GetMouse'))
				priorityLevel=min(priorityLevel,1);
				match=1;
			end
			if streq(name,lower('GetClicks'))
				priorityLevel=min(priorityLevel,1);
				match=1;
			end
            %****cases which are just here to check for valid arguments****
            if streq(name,lower('SND'))
				priorityLevel=min(priorityLevel,0);
				match=1;
			end
			if streq(name,lower('sound'))
				priorityLevel=min(priorityLevel,0);
				match=1;
			end
			if streq(name,lower('speak'))
				priorityLevel=min(priorityLevel,0);
				match=1;
			end
			if streq(name,lower('CopyWindow'))
				match=1;
			end
			if streq(name,lower('WaitBlanking')) | streq(name,lower('WaitVBL'))
				match=1;
			end
			if streq(name,lower('PeekBlanking')) | streq(name,lower('PeekVBL'))
				match=1;
			end
			if streq(name,lower('BlankingInterrupt'))
				match=1;
			end
			if streq(name,lower('SetClut'))
				match=1;
			end
			if streq(name,lower('GetSecs')) | streq(name,lower('WaitSecs'))
				match=1;
			end
			if streq(name,'cputime')
				match=1;
			end
			if streq(name,lower('EventAvail'))
				match=1;
			end
			if streq(name,lower('GetTicks'))
				match=1;
			end
			if streq(name,lower('KbCheck'))
				priorityLevel=min(priorityLevel,1);
				match=1;
			end
			if streq(name,lower('KbWait'))
				priorityLevel=min(priorityLevel,1);
				match=1;
			end
			if streq(name,lower('fopen'))
				match=1;
			end
			if streq(name,lower('fclose'))
				match=1;
			end
			if streq(name,lower('fprintf'))
				match=1;
			end
        end
		if isnumeric(varargin{i})
			w=varargin{i};
         % If a generic windowPtrOrScreenNumber argument is supplied, instead of the
         % more specific keywords, we lower priority to 1. This may be a bit
         % strict, but better safe than sorry...
         priorityLevel=min(priorityLevel,1);
         match=1;
		end
		if ~match
			error(['Unknown function ''' varargin{i} '''']);
		end
	end
%_______________________________________________________________________________________________________________________________
% OSX     
elseif IsOSX
        
    %validate the arguments and return 9.  
	if nargin<1
		error(['Usage: priorityLevel=MaxPriority([windowPtrOrScreenNumber],[''WaitBlanking''],[''PeekBlanking''],...' ...
		char(13) '                        [''BlankingInterrupt''],[''SetClut''],...'...
		char(13) '                        [''SND''],[''sound''],[''speak''],...'...
		char(13) '                        [''GetSecs''],[''WaitSecs''],[''cputime''],...'...
		char(13) '                        [''KbCheck''],[''KbWait''],[''CharAvail''],[''GetChar''],...'...
		char(13) '                        [''EventAvail''],[''GetClicks''],[''GetMouse''],[''GetTicks''])']);
    end %if
    %some of these things only exist on OS 9, accepting them as arguments
    %to MacPriority on OS X seems harmless and may be a minor convenience.  
    validFunctionArgsOS9=upper({'WaitBlanking','PeekBlanking','BlankingInterrupt','SetClut','SND','sound', 'speak', 'GetSecs', 'WaitSecs', ...
                    'cputime','KbCheck', 'KbWait', 'CharAvail','GetChar', 'EventAvail', 'GetClicks', 'GetMouse','GetTicks'});
    validWindowsAndScreens=union(Screen('Windows'), Screen('Screens'));
    for i=1:nargin
        %its a window pointer or screen number
		if isnumeric(varargin{i})
            if ~ismember(varargin{i}, validWindowsAndScreens);
                error('Invalid window pointer or screen number');
            end
        %its a function name
        elseif ischar(varargin{i})
            if ~ismember(upper(varargin{i}), validFunctionArgsOS9)
                error(['Unknown function ''' varargin{i} '''']);
            end
        %its something else
        else
			error([ 'argument ' num2str(i) ' is of wrong type']);
        end %else
    end %for
    priorityLevel=9;

%_______________________________________________________________________________________________________________________________
elseif IsLinux
    % We only differentiante between 0 = Normal priority, 1-99 = Realtime Priority mlocked Realtime.
    % We suggest 1 as a good level. The user is free to take higher numbers of course...
    priorityLevel=1;
else
    error('OS not supported by Psychtoolbox MaxPriority.m');
end
