% InstallKinect - Kinect driver installation instructions.
%
% MS-Windows:
% -----------
%
% 1. Unzip the file Psychtoolbox/PsychContributed/Kinect-v16-WindowsDrivers.zip
%    into a temporary folder, e.g., C:\tmp\KinectDrivers.
%
% 2. Plug in your Kinect box into a USB port via the interface cable.
%
% 3. The Windows hardware setup assistant will start, tell you about new
%    detected hardware and ask you for drivers. Say "No" to the
%    "automatically search for drivers" option, but select to "provide your
%    own drivers". In the next tab, select "browse for drivers", and then
%    select the temporary folder with the unzipped file (from step 1) as
%    driver folder. Press ok.
%
% 4. The driver will be installed from the zip file, the device manager
%    will notifiy you of the new device, then prompt you for installation
%    of another device. Repeat the same procedure from step 3 until no more
%    devices need to be installed. This procedure will repeat three to four
%    times until all drivers are installed, as the Kinect shows up as
%    multiple separate devices (Video camera, Audio soundcard and
%    Kinect motor).
%
% 5. The setup assistant will tell you that your new device is fully
%    operational. Quit the assistant.
%
% 6. Unplug the Kinect, wait a few seconds, plug it in again, just to be
%    sure it is correctly detected.
%
% 7. Now you can start Matlab or Octave and use the Kinect! Try KinectDemo
%    and Kinect3DDemo for a start, then delve into your own Kinect adventures.
%
% The current Kinect low level drivers are still early prototypes, so
% expect occassional bugs or weird behaviour.
%
%
% GNU/Linux:
% ----------
%
% If you use the Psychtoolbox distribution provided by the NeuroDebian
% project (http://neuro.debian.net) there's nothing to do. It should "just
% work(tm)", well almost: Skip to step 2. Otherwise the following step 1
% is required:
%
% 1. If you have Ubuntu Linux 10.04 LTS or later installed, open a terminal
% window and type this sequence of commands, providing your administrator
% password as requested. (Same procedure for Debian 4.0 or later)
%
%    a) Add the NeuroDebian repository to your software sources, as described
%       at http://neuro.debian.net#how-to-use-this-repository
%    b) sudo apt-get update
%    c) sudo apt-get install freenect
%    d) sudo adduser YOURNAME video
%       --> (YOURNAME) is your user login name!
%    e) Log out and Log in again.
%
% For non-Debian or non-Ubuntu Linux distributions, you'll need to install
% a version of libfreenect that is compatible to version 0.1.2 via whatever
% means your system provides to do this.
%
% 2. Kinect is now useable from within Matlab or Octave. Well almost.
% Systems with Linux kernel version 3.0 or later can use the video camera
% and microphones of the Kinect as regular sound and video devices, e.g.,
% for use by the Psychtoolbox videocapture and recording functions or other
% video apps (Skype, etc.). This however blocks use of the Kinect by our
% PsychKinect() driver. If you want the Kinect as 3D depths camera with
% our driver or other Kinect-specific software, you need to disable the
% standard Linux kinect driver "gspca_kinect" by black-listing it. On
% Ubuntu Linux (and probably most other distributions) you can do this
% by copying the file linux_blacklist_kinectvideo from Psychtoolboxs
% PsychContributed folder to Linux /etc/modprobe.d/ directory as a root
% user. This is most simply done by executing this function InstallKinect,
% and blindly entering your password while logged in as a user with
% administrator rights (as the script calls the sudo command).
%
% 3. After this procedure, the Kinect should be fully useable by Psychtoolbox.
%
%
% Mac OS/X:
% ---------
%
% See the OS/X section at http://openkinect.org/wiki/Getting_Started
% You will need libfreenect version 0.1.2 or compatible for this to work.
%
% The easiest way to get these libraries on OSX is via Homebrew.
% Get it at: http://mxcl.github.com/homebrew
%
% Once Homebrew is installed, one first needs to "brew install automake",
% and "brew install autoconfig", and potentially cmake before one can
% "brew install libfreenect".
%     
% PsychKinectCore links dynamically against those two libraries. We don't
% distribute them for now, as that would require us to distribute the
% corresponding source code of libusb-1.0.0 as well due to LGPLv2
% requirements.
%     
%
% CAVEATS:
% --------
%
% This is still early prototype software, expect bugs, bumps and hickups.
%
% The Kinect driver has been successfully tested with "Microsoft XBox
% Kinect". This version doesn't yet work with "Microsoft Kinect for
% Windows", due to lack of support in the libfreenect-0.1.2 library.
%
% Further installation instructions for other systems can be found at
% http://openkinect.org/wiki/Getting_Started
%
%
% Licenses: The driver consists of multiple components, which are licensed
% under different free software / open source licenses. The drivers are
% developed and licensed to you by their respective developers, the members
% of the OpenKinect community. See the respective web sites and licenses
% for copyright, authors, credits etc.:
%
% libusb - The underlying USB communication library is licensed under LGPL.
% libfreenect - The Kinect driver is dual-licensed (at your option) under
% GPLv2 or Apache license.
%
% The source code of the Windows version of the Kinect driver and libusb
% can be found after checking out Psychtoolbox's C source code from the SVN
% repository (see "help UseTheSource" for instructions) as a zip file
% under:
%
% trunk/PsychSourceGL/Cohorts/Kinect-v16-withsource.zip
%
% The source code for libfreenect (Unix aka Linux / MacOS/X) can be found
% under:
%
% web: http://openkinect.org/wiki/Contributing_Code
%
% The source code of libusb-1.0 can be found at:
%
% web: http://www.libusb.org/
%

cmd = ['sudo cp ' PsychtoolboxRoot '/PsychContributed/linux_blacklist_kinectvideo /etc/modprobe.d/'];
fprintf('Will copy the Kinect video driver blacklist file to your system. This will require\n');
fprintf('administrator root permission. Please enter your admin password now. You will likely\n');
fprintf('not see any visual feedback until you pressed ENTER and entered a valid password.\n');
drawnow;
[rc, msg] = system(cmd);
if rc == 0
    fprintf('Success! Now disabling the Kinect video driver. If this does not complete within a second,\n');
    fprintf('you may need to blindly type your password again + ENTER, but usually it just completes.\n');
    drawnow;
    [rc, msg] = system('sudo rmmod gspca_kinect');
    if rc == 0
        fprintf('Success! Your Kinect should now be useable by Psychtoolbox PsychKinect driver.\n');
    else
        fprintf('Failed! Maybe retry? Other than that, unplug your Kinect and reboot your machine to make it work.\n');
        fprintf('Reported error was: %s\n', msg);
    end
else
    fprintf('Failed! Maybe retry? Or ask a system administrator for help.\n');
    fprintf('Reported error was: %s\n', msg);
end

fprintf('If you want to enable use of Kinect as a regular webcam again, type the\n');
fprintf('following: sudo rm /etc/modprobe.d/linux_blacklist_kinectvideo\n');
fprintf('After that, unplug and replug your Kinect, or maybe reboot your machine.\n\n');

return;
