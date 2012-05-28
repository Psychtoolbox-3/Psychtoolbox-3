% PsychtoolboxKernelDriver - A low level support driver for Apple OS-X
%
% The PsychtoolboxKernelDriver (PKD) is a MacOS-X kernel extension (a
% kext). It currently should fully work with AMD/ATI Radeon graphics cards
% of the X1000, HD2000, HD3000 and HD4000 series. It works with the most
% recent HD-5000 series chips and later, although some functionality (10
% Bit native framebuffer support) is not yet available. On NVidia cards,
% beamposition queries are supported.
%
% The driver needs to be manually loaded by a user with administrator
% privileges and provides a few special services to PTB-3, ie., PTB's
% functionality is extended/enhanced if it detects such a driver at startup
% time:
%
% * Beamposition queries (See help BeampositionQueries) allow for
% especially robust and accurate stimulus onset timestamping. They are
% supported by OS-X on all PowerPC systems and on Intel based systems with
% NVidia graphics cards, but not on Intel systems with AMD/ATI graphics. The
% driver restores this functionality.
%
% * Multihead display sync: The driver allows to synchronize the display
% refresh cycles of the displays connected to a multihead graphics card
% to allow for high quality tear-free binocular and stereo stimulation of
% subjects.
%
% * Use of 10 bit per color component framebuffers: The driver enables two
% extra bits of color output precision per color channel on your card,
% allowing for 1 billion shades of different colors instead of the 16.8
% million colors available without the driver.
%
% * The driver implements workarounds to fix some problems caused by
% graphics driver and operating system bugs when the graphics card is used
% with high color/luminance precision display devices like the CRS Bits+ or
% Bits# boxes, or the VPixx Inc. DataPixx and ViewPixx devices, and similar
% equipment.
%
%
% The driver only works if a *single* graphics card (single-head or
% dual/multi-head) is active in your machine. It may not work reliably if
% more than one card is installed and or active. Specifically, we never
% tested how (well) it behaves in the presence of multiple active graphics
% cards.
%
% As this driver accesses the hardware at a very low level, bypassing the
% whole operating system, its graphics subsystem and the drivers from ATI,
% there is some potential for things to go wrong. Although our testing so
% far didn't show any such problems, it may happen on your system. That is
% why this driver is an *experimental* feature and why you need to have
% administrator privileges to install and load the driver the first time.
%
% How to install (one time setup):
% --------------------------------
%
% There are two versions of the driver, one for 32-Bit OSX kernels, one for
% 64-Bit kernels. You can find out which one you need by typing:
% !uname -m
%
% If the output says "x86_64" you'll need the 64-Bit driver, otherwise the
% 32-Bit driver.
%
% For 32-Bit OSX kernel:
% ----------------------
%
% Unzip the file
% Psychtoolbox/PsychHardware/PsychtoolboxKernelDriver32Bit.zip, then copy
% the unpacked driver PsychtoolboxKernelDriver.kext into the system folder
% /System/Library/Extensions/
% You must do this from the terminal via:
%
% sudo cp -R PsychtoolboxKernelDriver.kext /System/Library/Extensions/
%
% assuming the unzipped driver is in your current working directory.
%
% For 64-Bit OSX kernel:
% ----------------------
%
% Unzip the file
% Psychtoolbox/PsychHardware/PsychtoolboxKernelDriver64Bit.zip, then copy
% the unpacked driver PsychtoolboxKernelDriver.kext into the system folder
% /System/Library/Extensions/
% You must do this from the terminal via:
%
% sudo cp -R PsychtoolboxKernelDriver.kext /System/Library/Extensions/
%
% assuming the unzipped driver is in your current working directory.
%
%
% On modern OSX systems, this will automatically load the driver after a
% few seconds, and also after each system restart. On older systems you
% have to load the driver manually - on new systems you could do this for
% debugging purposes: How to enable (each time after restarting your system):
%
%
% In a terminal type:
% sudo kextload /System/Library/Extensions/PsychtoolboxKernelDriver.kext
%
% Both sudo commands will ask you for your password.
%
% If everything went well, Psychtoolbox will report availability (and its
% use) of the driver the first time the Screen() mex file gets loaded and
% used. If it doesn't happen immediately, type "clear Screen" to force a
% reload of Screen.
%
