libdc1394 libraries:

LibDC1394 is free software, licensed under LGPL, therefore
compatible with PTB's license and included here.

The official website of libdc1394 is currently:

http://sourceforge.net/projects/libdc1394/

=> Screen makes use of libdc1394 for its Firewire video capture engine
   under GNU/Linux and MacOS/X.

How to rebuild:

On OS/X you'd follow the following instructions. On Linux you'll usually
just download the latest development libraries with your distributions
package management system.

The libraries bundled in PsychSourceGL are the ones the official PTB
distro is linked against, but you can of course download more recent
versions of the libraries from their official home anytime and use
these.

1. Unpack the PsychSourceGL/Cohorts/libDC1394/libdc1394...tar.gz
   file into a local directory.

2. Follow the installation instructions inside the directory to
   build and install the library on your system.

3. By default, we install the library build products to
   /usr/local/lib and link statically against /usr/local/lib/libdc1394.a
   and have the header files in /usr/local/include

   -> Change the XCode project settings or other makefile settings to
      adapt for different paths or dynamic linking.

4. Rebuild Screen.

-> Done.

Note: On Linux we always link dynamically against the system installed lidc1394, this way the Linux package management system can keep the library up to date with the latest features, enhancements and bug-fixes.

On OSX for 64-Bit pub we also weak-link dynamically. We use libdc v-2.2.0, which is contained in the libdc1394-2.2.0.tar file. A precompiled v2.2.0 64-Bit .dylib is contained inside the Psychtoolbox/PsychVideoCapture/ subfolder.


In PsychSourceGL/Cohorts/ARToolkit:

For OS/X builds, the precompiled static link library libdc1394.a (.a's) is
stored inside PsychSourceGL/Projects/MacOSX/PsychToolbox/libDC1394/
This one is just stored here for convenience to save you the time for
recompiling yourself.

==> All content can be derived from the included .tar.gz file!
