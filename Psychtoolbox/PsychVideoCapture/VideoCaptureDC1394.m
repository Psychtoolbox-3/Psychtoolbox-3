% VideoCaptureDC1394 - Setup instructions for IEEE-1394 firewire video
% capture.
%
% Psychtoolbox on GNU/Linux and MacOSX supports a special video capture
% engine for IIDC compliant machine vision cameras connected via IEEE-1394
% firewire bus or via USB bus and Firewire-over-USB protocol.
%
% This engine with engineID 1, as selectable via ...
% Screen('Preference', 'DefaultVideocaptureEngine', 1)
% ... provides especially precise low-level control over the cameras, high
% framerates, low-latency and high timing precision.
%
% The engine is based on the open-source library libdc1394. For more info
% about it see 'help VideoCapture' in the section about professional class
% video capture engines.
%
% INSTALLATION:
%
% On Linux, usually no installation is required if you use the Psychtoolbox
% provided by Debian or Ubuntu via the Debian main package archive or the
% NeuroDebian project. If you have a different Linux distribution, you'll
% need to install the "libdc1394" package via your distributions software
% manager, e.g., via "sudo apt-get install libdc1394" on a Debian based
% system.
%
% The software manager on Linux will keep your system up to date with the
% latest stable version of the library.
%
% On MacOS/X for 32-Bit Matlab & Octave: Screen() is currently statically
% linked against a version of the library. If you want access to its
% corresponding source-code, read "help UseTheSource" on how to find the
% PTB C-Source code. The PsychSourceGL/Cohorts/libDC1394/ subfolder
% contains the source in libdc1394-2.0.1.tar.gz at the moment of this
% writing.
%
% On MacOS/X for 64-Bit Matlab & Octave: Screen() links dynamically against
% a system-installed .dylib version of the library, if such a version is
% installed on your system. A precompiled library can be found in the
% Psychtoolbox/PsychVideoCapture/ subfolder as libdc1394.22.dylib. Copy it
% into the /usr/local/lib/ folder via executing this function
% 'VideoCaptureDC1394' and entering your administrator password on request.
% The corresponding source code of this library can be found under
% PsychSourceGL/Cohorts/libDC1394/ as libdc1394-2.2.0.tar at the moment of
% this writing. However, you can download the most recent copies of the
% libraries from libDC's sourceforge website anytime and install it instead
% to get access to the latest features and bugfixes.
%

% History:
% 6.05.2012  mk  Written.

% Perform installation on 64-Bit OS/X:
if IsOSX(1)
    % Check if /usr/local/lib/libdc1394... file already exists:
    needinstall = 0;
    if ~exist('/usr/local/lib/libdc1394.dylib', 'file')
        % No: Needs to be installed.
        needinstall = 1;
        fprintf('The libdc1394.dylib file for Psychtoolbox is not installed on your system.\n');
        answer = input('Should i install it? [y/n] : ', 's');
    else
        fprintf('LibDC1394 already installed. Nothing to do.\n');
        return;
    end
    
    if needinstall && answer == 'y'
        fprintf('I will copy my most recent libDC1394 file to your system. Please enter\n');
        fprintf('now your system administrator password. You may not see any feedback.\n');
        drawnow;
        
        cmd = sprintf('sudo cp %s/PsychVideoCapture/libdc1394.22.dylib /usr/local/lib/libdc1394.22.dylib', PsychtoolboxRoot);
        [rc, msg] = system(cmd);
        if rc == 0
            cmd = sprintf('sudo ln -s /usr/local/lib/libdc1394.22.dylib /usr/local/lib/libdc1394.dylib');
            [rc, msg] = system(cmd);
        end
        
        if rc == 0
            fprintf('Success! You may need to restart Matlab or Octave for some changes to take effect.\n');
        else
            fprintf('Failed! The error message was: %s\n', msg);
        end
    end
    
    return;
end
