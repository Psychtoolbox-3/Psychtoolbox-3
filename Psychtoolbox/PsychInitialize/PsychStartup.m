function PsychStartup
% PsychStartup -- Perform Psychtoolbox related setup at Matlab/Octave startup.
%
% This performs setup of Matlab or Octave at session startup for
% Psychtoolbox.
%
% On MS-Windows with 64-Bit Matlab it currently detects if the GStreamer
% SDK runtime is installed, as this is required for Screen() to work. It
% performs GStreamer setup, or outputs a warning if the runtime is missing.
%
% This function is normally called from the startup.m Matlab startup file.
%

% History:
% 14.09.2012  mk  Written.
% 14.01.2013  mk  Make path detection more robust.

% Try-Catch protect the function, so Matlab startup won't fail due to
% errors in this function:
try
    % Setup for 64-Bit Matlab on MS-Windows:
    if IsWin(1) && ~IsOctave
        % Need to assign a proper install location for the
        % www.gstreamer.com GStreamer-SDK runtime libraries, otherwise
        % loading of the Screen() MEX file would fail due to unresolved
        % link dependencies:
        
        % Find path to SDK-Root folder: Should be defined in env variable
        % by installer:
        sdkroot = getenv('GSTREAMER_SDK_ROOT_X86_64');
        if ~exist(sdkroot, 'dir')
            % Env variable points to non-existent SDK dir. How peculiar?
            % Invalidate invalid sdkroot, so fallback code can run:
            fprintf('PsychStartup: Environment variable GSTREAMER_SDK_ROOT_X86_64 points to non-existent SDK folder?!?\n');
            fprintf('PsychStartup: The missing or inaccessible path to GStreamer is: %s\n', sdkroot);
            fprintf('PsychStartup: Something is botched. Trying various common locations for the SDK to keep going.\n');
            sdkroot = [];
        end

        % Probe standard install location on drives C,D,E,F,G:
        if isempty(sdkroot) && exist('C:\gstreamer-sdk\0.10\', 'dir')
            sdkroot = 'C:\gstreamer-sdk\0.10\x86_64\';
        end

        if isempty(sdkroot) && exist('D:\gstreamer-sdk\0.10\', 'dir')
            sdkroot = 'D:\gstreamer-sdk\0.10\x86_64\';
        end

        if isempty(sdkroot) && exist('E:\gstreamer-sdk\0.10\', 'dir')
            sdkroot = 'E:\gstreamer-sdk\0.10\x86_64\';
        end

        if isempty(sdkroot) && exist('F:\gstreamer-sdk\0.10\', 'dir')
            sdkroot = 'F:\gstreamer-sdk\0.10\x86_64\';
        end

        if isempty(sdkroot) && exist('G:\gstreamer-sdk\0.10\', 'dir')
            sdkroot = 'G:\gstreamer-sdk\0.10\x86_64\';
        end
        
        % Probe install locations of legacy SDK's:
        if isempty(sdkroot) && exist('C:\gstreamer-sdk\2012.5\', 'dir')
            sdkroot = 'C:\gstreamer-sdk\2012.5\x86_64\';
        end
        
        if isempty(sdkroot) && exist('C:\gstreamer-sdk\2012.7\', 'dir')
            sdkroot = 'C:\gstreamer-sdk\2012.7\x86_64\';
        end
        
        if isempty(sdkroot) && exist('D:\gstreamer-sdk\2012.5\', 'dir')
            sdkroot = 'D:\gstreamer-sdk\2012.5\x86_64\';
        end
        
        if isempty(sdkroot) && exist('D:\gstreamer-sdk\2012.7\', 'dir')
            sdkroot = 'D:\gstreamer-sdk\2012.7\x86_64\';
        end
        
        if isempty(sdkroot)
            fprintf('\nPsychStartup: Path to GStreamer runtime is undefined! This probably means that\n');
            fprintf('PsychStartup: the 64-Bit GStreamer SDK from www.gstreamer.com is not installed.\n');
            fprintf('PsychStartup: The Psychtoolbox Screen() mex file will fail to work until you fix\n');
            fprintf('PsychStartup: this! Read ''help GStreamer'' for instructions.\n\n');
        else
            sdkroot = [sdkroot 'bin'];
            
            fprintf('\nPsychStartup: Adding path of installed GStreamer runtime to library path. [%s]\n', sdkroot);
            
            % Get current path:
            path = getenv('PATH');
            
            % Prepend sdkroot to path:
            newpath = [sdkroot ';' path];
            setenv('PATH', newpath);
        end
    end
catch %#ok<*CTCH>
    fprintf('PTB-WARNING: Call to PsychStartup() failed!!\n');
end

end
