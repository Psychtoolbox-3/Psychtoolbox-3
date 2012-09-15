% Psychtoolbox:PsychContributed
%
% Scripts and software contributed by Psychtoolbox users and
% third-party tools and software that are used in Psychtoolbox. These
% third-party tools are covered by different licenses which allow distribution
% with Psychtoolbox, but are not neccessarily free-software or open-source.
% See the descriptions below or the help files of the functions for
% exact licensing conditions of individual contributions.
%
% help Psychtoolbox % For an overview, triple-click me & hit enter.
%
%
% linux_blacklist_kinectvideo  Blacklist file for Linux, to be copied into the
%                              /etc/modprobe.conf/ directory. Disables the Linux
%                              standard video driver for Microsoft Kinect box. Allows
%                              use of the Kinect with PsychKinect() driver et al.
%                              See 'help InstallKinect'.
%
% ATIRadeonperf                ATIRadeonperf_Linux and ATIRadeonperf_Windows.exe are
%                              compiled Linux and Windows executables that are called
%                              as helpers from the function PsychGPUControl() via a
%                              system() call, ie., they execute in the system shell as
%                              separate processes.
%
%                              These binaries perform the actual GPU setup work on behalf
%                              of PsychGPUControl() using ATI/AMD proprietary low-level
%                              API's. The API's are available to registered developers from
%                              AMD/ATI as "ADL SDK" without royalty fees. While the SDK
%                              license does allow redistribution of executables using the
%                              SDK, it doesn't allow redistribution in source code form.
%                              For that reason, these tools are only made available as
%                              precompiled binaries for your use, not under any open-source
%                              or free-software license. The tools were compiled by Mario
%                              Kleiner, so contact him if you have further questions about
%                              the tools.
%
% Kinect-v16-WindowsDrivers.zip Zip file containing the Microsoft Windows
%                               drivers and installer that are needed to
%                               use the Microsoft XBOX Kinect under
%                               Windows. Type "help PsychKinect" for more
%                               infos. The PsychSourceGL/Cohorts/Kinect-v16-withsource.zip
%                               file contains corresponding source code.
%
% libusb-1.0.dll               This is the 32 bit Windows experimental version of
%                              libusb-1.0 from ...
%                              http://libusb.org/wiki/windows_backend
%                              The dll is used by PsychHID on MS-Windows.
%                              It is licensed under LGPL v2+. You can get a
%                              more recent version of this dll from the
%                              website mentioned above if you like.
%
%                              The corresponding source code can be found
%                              in the PsychSourceGL/Cohorts/libusb1-win32/
%                              subfolder after you have downloaded the full
%                              ptb source (help UseTheSource).
%
% x64/libusb-1.0.dll           This is the 64-Bit windows version for
%                              64-Bit Matlab.
%
%
% simplepsychtoolboxsetup.sh   Configure OS X priority.  Call from Unix shell, not MATLAB.
%
% macid.exe                    Microsoft Windows command line utility for retrieving the MAC
%                              ethernet address of the primary network
%                              adapter. This free software GPL utility is
%                              downloaded from www.mirkes.de Windows batch
%                              tools. Author: Markus Stephany
%                              URL:
%                              http://www.mirkes.de/en/freeware/batch.php
%
% macidpascalsource/           This folder contains the Delphi-7 Pascal source code
%                              of macid.exe, as required by macid's GPL license.
%
% ple                          Prints useful debugging output, including a stack trace in case
%                              a script exits via rethrow() or psychrethrow(). See help ple for
%                              copyright info.
%
% read_rle_rgbe                Reader routine for RLE encoded RGBE high dynamic range images.
%                              See help read_rle_rgbe for copyright & authorship info. This is
%                              a helper routine for ReadHDR - the recommended function for user
%                              scripts.
%
% vcredist_x86.exe             Microsoft Visual C++ 2010 runtime 32-Bit
%                              redistributable installer. You must use this
%                              to install the Visual C runtime if you want
%                              to use the Psychtoolbox for 32-Bit Matlab on
%                              Microsoft Windows.
%
% vcredist_x64.exe             Microsoft Visual C++ 2010 runtime 64-Bit
%                              redistributable installer. You must use this
%                              to install the Visual C runtime if you want
%                              to use the Psychtoolbox for 64-Bit Matlab on
%                              Microsoft Windows.
%
% WinJoystickMex               WinJoystickMex is a very simple/rudimentary
%                              joystick driver for basic query of position
%                              and button state of joysticks under
%                              MS-Windows, written by Mario Kleiner. The
%                              code is contributed under GPLv2.
%
% WinTab/                      This folder contains the (unsupported!)
%                              WinTabMex MEX file driver (and its C source
%                              code) for controlling touch/digitizer
%                              tablets via the WinTab API on Microsoft
%                              Windows operating systems. It also contains
%                              some basic usage instructions and two M-File
%                              demos, contributed by Andrew Wilson.
%
