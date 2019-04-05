As of Psychtoolbox 3.0.15, this folder mostly contains files of historic value.

The following files are meaningful and current as of PTB 3.0.15:

libportaudio_osx_64.a      = Static link library for building PsychPortaudio for 64-Bit macOS with static linking.
portaudio_x64.lib          = Import library used for building PsychPortaudio.mexw64 for 64-Bit Matlab on Windows.
MinGW64/portaudio_x64.lib  = Import library used for building PsychPortaudio.mex for 64-Bit Octave-4.4+ on Windows.

These files are built from upstream Portaudio V19.6.0+ (https://portaudio.com).
However, for the Windows variants, any trace of ASIO support code has been removed from the
upstream source before building, to have a completely ASIO free build.

The corresponding ASIO-less portaudio_x64.dll Windows runtime library is in the Psychtoolbox/PsychSound folder as usual.
It is build with MME, WDMKS and WASAPI backends, and without DirectSound or ASIO backends.

On the Linux side, we now just dynamically build and link against the distro installed Portaudio shared libraries,
which are usually of version 19.6.0+ on recent distros, and on older distros should be close enough to the 19.6.0
state to be useful for us.

All other files in this folder were only used for Psychtoolbox 3.0.14 and earlier, containing heavily
patched/modified/enhanced/bug-fixed source files of a very old Portaudio V19 devel.
We keep them here for reference for now, but will probably delete them soon, maybe once
i've decided if i want to forward-port our direct input monitoring support for macOS (but probably not).

History:

portaudio_unpatched_except4OSX.zip is the original portaudio version which we use as

basis for our own enhancements. Well, almost - the patches for OS/X Coreaudio are 
partially
applied to this codebase. 
To build portaudio from scratch, you'll need to unpack this zip file,
then replace 
the files inside the unzipped source tree by the patched files in this folder, 
overwriting the originals. For Windows+ASIO you'll need to download the freely available

Steinberg ASIO SDK from their website and copy its files into the proper portaudio 
folders
and set the proper build flags, as described in portaudio's documentation 
Finally, build the thing.

This zip file is just here to document our 'baseline'.
