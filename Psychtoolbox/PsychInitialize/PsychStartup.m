function PsychStartup
% PsychStartup -- Perform Psychtoolbox related setup at Matlab/Octave startup.
%
% This performs setup of Matlab or Octave at session startup for
% Psychtoolbox.
%
% On MS-Windows, it currently detects if the GStreamer 1.16+ runtime is
% installed, as this is required for Screen() multi-media functions to
% work. It performs GStreamer setup, or outputs a warning if the runtime is
% missing.
%
% This function is normally called from the startup.m Matlab startup file,
% or from the .octaverc startup file on GNU/Octave.
%

% History:
% 14.09.2012  mk  Written.
% 14.01.2013  mk  Make path detection more robust.
% 12.09.2013  mk  Also apply GStreamer-SDK setup to 32-Bit Matlab on Windows.
% 26.08.2014  mk  Switch the 64-Bit setup to the GStreamer-1.4.0+ runtime.
% 17.10.2015  mk  Add warning about complete failure of Screen() for Octave on Windows.
% 30.10.2015  mk  Setenv PSYCH_GSTREAMER_SDK_ROOT to help GStreamer on Octave-4 for Windows.
% 14.01.2016  mk  Make detection more robust, delete dead code for GStreamer 0.10.
% 07.08.2019  mk  Update for Psychtoolbox 3.0.16 - GStreamer 1.16+ MSVC needed.

% Try-Catch protect the function, so Matlab startup won't fail due to
% errors in this function:
try
    % Setup for MS-Windows:
    if IsWin
        % Need to assign a proper install location for the
        % www.gstreamer.net GStreamer-1.16.0+ MSVC runtime libraries, otherwise
        % use of GStreamer based functions would fail due to unresolved
        % link dependencies:

        % Must omit GStreamer runtime path on Octave in CLI mode, to avoid
        % build trouble due to lib clashes between Octave and GStreamer:
        useGStreamer = ~IsOctave || ~any(strfind(cell2mat (argv'), '--no-gui'));

        % Find path to SDK-Root folder: Should be defined in env variable
        % by installer:
        if Is64Bit
            % Use 64-Bit GStreamer runtime for 64-Bit Matlab/Octave on Windows:
            sdkroot = getenv('GSTREAMER_1_0_ROOT_X86_64');
            suffix = 'x86_64\';
        else
            % Use 32-Bit GStreamer runtime atm. for 32-Bit Octave-4.0 on Windows:
            sdkroot = getenv('GSTREAMER_1_0_ROOT_X86');
            suffix = 'x86\';
        end
        if isempty(sdkroot) %try alternative Windows enviroment variable name
            if Is64Bit
                % Use 64-Bit GStreamer runtime for 64-Bit Matlab/Octave on Windows:
                sdkroot = getenv('GSTREAMER_1_0_ROOT_MSVC_X86_64');
            else
                % Use 32-Bit GStreamer runtime atm. for 32-Bit Octave-4.0 on Windows:
                sdkroot = getenv('GSTREAMER_1_0_ROOT_MSVC_X86');
            end
        end
        if isempty(sdkroot)
            if Is64Bit
                fprintf('PsychStartup: Environment variable GSTREAMER_1_0_ROOT_X86_64 is undefined.\n');
            else
                fprintf('PsychStartup: Environment variable GSTREAMER_1_0_ROOT_X86 is undefined.\n');
            end
            fprintf('PsychStartup: Either GStreamer-1.16 MSVC is not installed at all, or if it is installed then something\n');
            fprintf('PsychStartup: is botched. Trying various common locations for the GStreamer runtime to keep going.\n');
        else
            if ~exist(sdkroot, 'dir')
                % Env variable points to non-existent SDK dir. How peculiar?
                % Invalidate invalid sdkroot, so fallback code can run:
                if Is64Bit
                    fprintf('PsychStartup: Environment variable GSTREAMER_1_0_ROOT_X86_64 points to non-existent folder?!?\n');
                else
                    fprintf('PsychStartup: Environment variable GSTREAMER_1_0_ROOT_X86 points to non-existent folder?!?\n');
                end

                fprintf('PsychStartup: The missing or inaccessible path to GStreamer is: %s\n', sdkroot);
                fprintf('PsychStartup: Something is botched. Trying various common locations for the SDK to keep going.\n');
                sdkroot = [];
            end
        end

        % Probe standard install location on drives C,D,E,F,G:
        if isempty(sdkroot) && exist('C:\gstreamer\1.0\', 'dir')
            sdkroot = ['C:\gstreamer\1.0\' suffix];
        end

        if isempty(sdkroot) && exist('D:\gstreamer\1.0\', 'dir')
            sdkroot = ['D:\gstreamer\1.0\' suffix];
        end

        if isempty(sdkroot) && exist('E:\gstreamer\1.0\', 'dir')
            sdkroot = ['E:\gstreamer\1.0\' suffix];
        end

        if isempty(sdkroot) && exist('F:\gstreamer\1.0\', 'dir')
            sdkroot = ['F:\gstreamer\1.0\' suffix];
        end

        if isempty(sdkroot) && exist('G:\gstreamer\1.0\', 'dir')
            sdkroot = ['G:\gstreamer\1.0\' suffix];
        end

        if isempty(sdkroot) || ~exist(sdkroot, 'dir')
            fprintf('\nPsychStartup: Path to GStreamer runtime is undefined! This probably means that\n');
            fprintf('PsychStartup: the 64-Bit GStreamer 1.16+ MSVC runtime from www.gstreamer.net is not installed.\n');

            % PTB 3.0.16 - GStreamer always needed, also on Matlab.
            if 1 || IsOctave
                fprintf('PsychStartup: The Psychtoolbox Screen() function will not work at all until you fix\n');
            else
                fprintf('PsychStartup: The Psychtoolbox Screen() multimedia functions and the new ''DrawText''\n');
                fprintf('PsychStartup: high quality text renderer will fail to work until you fix\n');
            end

            fprintf('PsychStartup: this! Read ''help GStreamer'' for instructions.\n\n');
        else
            sdkroot = [sdkroot 'bin'];

            % Get current path:
            path = getenv('PATH');

            if useGStreamer
                % Matlab, or Octave in GUI mode: Prepend sdkroot to path:
                newpath = [sdkroot ';' path];

                % Check if we have the right flavor of GStreamer, the MSVC variant:
                addpath(sdkroot);
                if exist('libgstreamer-1.0-0.dll', 'file') && ~exist('gstreamer-1.0-0.dll', 'file')
                    % Wrong type, the MinGW build instead of the MSVC build!
                    fprintf('\n\n');
                    fprintf('PsychStartup: WRONG variant of the GStreamer runtimes installed! This is\n');
                    fprintf('PsychStartup: the MinGW variant, but Psychtoolbox needs the MSVC variant!\n');
                    fprintf('PsychStartup: The Psychtoolbox Screen() function will NOT WORK AT ALL until\n');
                    fprintf('PsychStartup: you fix this! Read ''help GStreamer'' for instructions.\n\n');
                    warning('WRONG TYPE OF GStreamer packages INSTALLED! WE NEED THE MSVC variant, this is the MINGW variant!');
                end
                rmpath(sdkroot);

                fprintf('\nPsychStartup: Adding path of installed GStreamer runtime to library path. [%s]\n', sdkroot);
            else
                % Octave CLI: Do not add GStreamer runtime, so we can build
                % our mex files without interference of GStreamer's libstdc++.6.dll:
                newpath = path;
                fprintf('\nPsychStartup: Use of GStreamer disabled for Octave CLI session.');
                fprintf('\nPsychStartup: Screen() mex file will not work.\n');
            end

            % For Octave-5.1 on Windows, also need to prepend path to portaudio dll
            % to library path, otherwise PsychPortAudio won't load:
            driverloadpath = [PsychtoolboxRoot 'PsychSound'];
            newpath = [driverloadpath ';' newpath];

            % For moglcore we need to prepend the path to freeglut.dll:
            if Is64Bit
                % Need 64-Bit freeglut.dll:
                driverloadpath = [PsychtoolboxRoot 'PsychOpenGL/MOGL/core/x64'];
            else
                % Need 32-Bit freeglut.dll:
                driverloadpath = [PsychtoolboxRoot 'PsychOpenGL/MOGL/core'];
            end
            newpath = [driverloadpath ';' newpath];

            % For PsychHID we need to prepend the path to libusb-1.0.dll:
            if Is64Bit
                % 64-Bit version of libusb.dll
                driverloadpath = [PsychtoolboxRoot 'PsychContributed' filesep 'x64' filesep];
            else
                % 32-Bit version of libusb.dll
                driverloadpath = [PsychtoolboxRoot 'PsychContributed' filesep ];
            end
            newpath = [driverloadpath ';' newpath];

            setenv('PATH', newpath);

            % Also store sdkroot in a separate environment variable, to be used
            % by Screen for Octave internally:
            setenv('PSYCH_GSTREAMER_SDK_ROOT', sdkroot);
        end
    end

catch %#ok<*CTCH>
    fprintf('PTB-WARNING: Call to PsychStartup() failed!!\n');
end

end
