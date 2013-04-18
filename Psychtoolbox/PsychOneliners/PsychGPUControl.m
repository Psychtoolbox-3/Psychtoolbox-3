function rc = PsychGPUControl(cmd, varargin)
% rc = PsychGPUControl(cmd, arg); -- Control low-level GPU settings.
%
% PsychGPUControl calls into external helper tools to change certain
% low-level operating settings of your systems graphics card (GPU).
%
% Not all operating systems and GPU's support this. The function will
% do nothing on unsupported OS/GPU combos.
%
% Currently OS/X doesn't support this function at all, and on MS-Windows
% and GNU/Linux, only recent ATI GPU's with recent drivers do support it.
%
% All subfunctions return an optional 'rc' return code of zero on success,
% non-zero on error or if the feature is unsupported.
%
% Subfunctions and their syntax & meaning:
%
% rc = PsychGPUControl('SetDitheringEnabled', enableFlag);
% - Depending on the setting of 'enableFlag', either enable (=1) or
% disable (=0) display color dithering on all connected displays.
%
% Under normal circumstances, the GPU should decide itself if dithering
% should be used or not. This function allows you to override the GPU's
% automatic choice.
%
%
% rc = PsychGPUControl('SetGPUPerformance', gpuPerformance);
% - Select the performance state of the GPU. 'gpuPerformance' can be set to
% 0 if the GPU shall automatically adjust its performance and power-
% consumption, or to one of 10 fixed levels between 1 and 10, where 1 means
% the lowest performance - and power consumption, whereas 10 means the
% highest performance - and maximum power consumption.
%
% If in doubt, choose 10 for best accuracy of visual stimulus onset timing,
% 0 for non-critical activities to leave the decision up to the graphics
% driver and GPU.
%
%
% rc = PsychGPUControl('FullScreenWindowDisablesCompositor', flag [, screenIds]);
% - Select if desktop composition should be disabled for displays where
% a Psychtoolbox fullscreen onscreen window is displayed. 'flag' == 1 means
% to disable composition for fullscreen windows, 0 means to enable composition
% for fullscreen windows. You usually want composition to be disabled, as this
% is currently the only way to get decent timing and precise visual stimulus
% onset timestamping. The optional vector of 'screenIds' selects which screens
% should be affected by the change. If left out or set to [], all detected
% screens will be changed.
%
% 

% History:
%  3.01.2010  mk  Written.
% 19.04.2010  mk  Add quotes around path to command to protect against
%                 blanks in path to executable.
% 16.01.2011  mk  Add function to control desktop composition on Linux with
%                 Compiz.
% 18.04.2013  mk  Add use of 64-Bit ATIRadeonperf_Linux64 exe on 64-Bit Linux.
%

if nargin < 1
	error('Subfunction command argument missing!');
end

if strcmpi(cmd, 'SetDitheringEnabled')
	if isempty(varargin)
		error('SetDitheringEnabled: enableFlag argument missing!');
	end

	enable = varargin{1};
	if ~ismember(enable, 0:1)
		error('SetDitheringEnabled: Invalid enableFlag argument, not 0 or 1!');
	end

	% Command code 1 means: Control ditherstate, according to 2nd arg 0 or 1 == disable, enable.
	cmdpostfix = sprintf(' 1 %i', enable);
	rc = executeRadeoncmd(cmdpostfix);
	return;
end

if strcmpi(cmd, 'SetGPUPerformance')
	if isempty(varargin)
		error('SetGPUPerformance: gpuPerformance argument missing!');
	end

	gpuperf = varargin{1};
	if ~ismember(gpuperf, 0:10)
		error('SetGPUPerformance: Invalid gpuPerformance argument, not an integer in range 0 - 10!');
	end

	% Map range 1 to 5 to "minimum performance" on ATI GPU's:
	if gpuperf > 0 && gpuperf <= 5
		perfflag = 2;
	end

	% Map range 6 to 10 to "maximum performance" on ATI GPU's:
	if gpuperf > 5 && gpuperf <= 10
		perfflag = 1;
	end

	if gpuperf == 0
		perfflag = 0;
	end

	% Command code 0 means: Control performance state: 0 = AUTO, 1 = MAX, 2 = MIN.
	cmdpostfix = sprintf(' 0 %i', perfflag);
	rc = executeRadeoncmd(cmdpostfix);
	return;
end

if strcmpi(cmd, 'FullScreenWindowDisablesCompositor')
	if isempty(varargin)
		error('FullScreenWindowDisablesCompositor: flag argument missing!');
	end

	compositorOff = varargin{1};
	if ~ismember(compositorOff, 0:1)
		error('FullScreenWindowsDisableCompositor: Invalid flag argument, not an integer of value 0 or 1!');
	end

	if length(varargin) < 2 || isempty(varargin{2})
		screenIds = Screen('Screens');
	else
		screenIds = varargin{2};
	end

	% Which OS?
	if ~IsLinux
		% Nothing to do on non-Linux as compositor handling is
		% implemented in Screen internally:
		rc = 1;
		return;
	end

	% We only know how to do this for Compiz, so we try that. Settings are persistent
	% across sessions and take effect immediately:
	rc = [];
	for screenId=screenIds
		if compositorOff
			% Enable un-redirection: Fullscreen windows aren't subject to treatment by compositor,
			% but can do (e.g. page-flipping) whatever they want:
			newstate = 'dis';
			rc(end+1) = system(sprintf('gconftool-2 -s --type bool /apps/compiz/general/screen%i/options/unredirect_fullscreen_windows true', screenId));
			rc(end+1) = system(sprintf('gconftool-2 -s --type bool /apps/compiz-1/plugins/composite/screen%i/options/unredirect_fullscreen_windows true', screenId));
			rc(end+1) = system(sprintf('dconf write /org/compiz/profiles/unity/plugins/composite/unredirect-fullscreen-windows true'));
		else
			% Disable un-redirection: Fullscreen windows get composited as all other windows:
			newstate = 'en';
			rc(end+1) = system(sprintf('gconftool-2 -s --type bool /apps/compiz/general/screen%i/options/unredirect_fullscreen_windows false', screenId));
			rc(end+1) = system(sprintf('gconftool-2 -s --type bool /apps/compiz-1/plugins/composite/screen%i/options/unredirect_fullscreen_windows false', screenId));
			rc(end+1) = system(sprintf('dconf write /org/compiz/profiles/unity/plugins/composite/unredirect-fullscreen-windows false'));
		end

		if ~all(rc)
			rc = 0;
			fprintf('PsychGPUControl:FullScreenWindowDisablesCompositor: Desktop composition for fullscreen windows on screen %i %sabled.\n', screenId, newstate);
		else
			rc = 1;
			fprintf('PsychGPUControl:FullScreenWindowDisablesCompositor: FAILED to %sable desktop composition for fullscreen windows on screen %i!\n', newstate, screenId);
			fprintf('This can cause visual onset timing problems! See ''help SyncTrouble'' - the Linux specific subsection for tips.\n');
		end
	end
	return;
end

error('Invalid subfunction provided. Read the help for valid commands!');
return; %#ok<UNRCH>
end

function rc = executeRadeoncmd(cmdpostfix)
    if IsOSX
        % A no-op on OS/X, as this is not supported at all.
        rc = 1;
        return;
    end

    if IsLinux
        % For 32-Bit Linux:
        cmdprefix = '/PsychContributed/ATIRadeonperf_Linux ';
    end

    if IsLinux(1)
        % For 64-Bit Linux: We have a dedicated 64-Bit executable, so we
        % do not require installation of 32-Bit compatibility libraries on
        % 64-Bit Linux systems just for our little exe here, as that would
        % be wasteful:
        cmdprefix = '/PsychContributed/ATIRadeonperf_Linux64 ';
    end

    if IsWin
        % For 32-Bit or 64-Bit Windows we always use a 32-Bit executable:
        cmdprefix = '/PsychContributed/ATIRadeonperf_Windows ';
    end

    % Create quoted version of path, so blanks in path are handled properly:
    doCmd = strcat('"', [PsychtoolboxRoot cmdprefix] ,'"');

    % Call final command, return its return status code:
    [rc, msg] = system([doCmd cmdpostfix]);

    % Code has it backwards 1 = success, 0 = failure. Remap to our
    % convention:
    rc = 1 - rc;
    
    % Output potential status or error messages, unless it is a message
    % that signals we are not executing on a AMD/ATI GPU with Catalyst
    % driver, ie., that the whole thing was a no-op:
    if isempty(strfind(msg, 'ADL library not found!'))
        disp(msg);
    end

    return;
end
