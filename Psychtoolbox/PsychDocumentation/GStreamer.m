% GStreamer - Installation instructions for the GStreamer media framework.
%
% Psychtoolbox for all versions of GNU/Linux and MS-Windows, as well as
% 64 bit versions of Mac OS/X uses the GStreamer multi-media framework for
% all multi-media related operations. For 32 bit versions of Matlab on
% OS/X, use of GStreamer is still optional, Apple Quicktime is the default
% on such setups, but Quicktime is a deprecated technology in Psychtoolbox,
% so we encourage you to use GStreamer on all system configurations.
%
% All movie playback, movie creation, video capture and video recording
% operations are based on GStreamer. These functions won't work without a
% working GStreamer installation on your system.
%
% You will need at least version 0.10.24 of GStreamer, but we recommend to
% use the latest available versions of the 0.10.x series, as not all
% Psychtoolbox features are fully supported by version 0.10.24. Especially
% video capture, video recording and movie creation may be limited in
% functionality or performance if you don't use a recent version of
% GStreamer.
%
% Currently, the 32-bit version of Psychtoolbox for MacOS/X uses Apple's
% Quicktime for those operations by default, so GStreamer is not strictly
% needed there. You can however use GStreamer (once installed) instead of
% Quicktime if you add the following command to the top of your experiment
% script, before any other Screen() command:
%
% This uses GStreamer instead of Quicktime for movie playback:
% Screen('Preference', 'OverrideMultimediaEngine', 1);
%
% This uses GStreamer instead of Quicktime for video capture or video
% recording:
%
% Screen('Preference', 'DefaultVideocaptureEngine', 3);
%
% So what's the deal with GStreamer?
%
% GStreamer (Webpage web: http://www.gstreamer.net ) is a free software /
% open-source multimedia framework. The summary on its webpage summarizes
% it as follows:
%
% "GStreamer is a library for constructing graphs of media-handling
% components. The applications it supports range from simple Ogg/Vorbis
% playback, audio/video streaming to complex audio (mixing) and video
% (non-linear editing) processing.
%
% Applications can take advantage of advances in codec and filter
% technology transparently. Developers can add new codecs and filters by
% writing a simple plugin with a clean, generic interface.
%
% License: GStreamer is released under the LGPL. The 0.10 series is API and
% ABI stable."
%
% In summary, GStreamer is for the free software world what Quicktime is
% for Apple's "walled garden". It is used a the standard media framework
% for many (most?) Linux desktop systems, netbooks, and MeeGo, Web-OS (and
% also Android) based tablets and Smartphones. It is also used on many 
% Linux based embedded devices.
%
% Installation instructions:
%
% GNU/Linux: Any half-way recent Linux distribution will include support
% for GStreamer in its package management system, so you can easily install
% it via the software management tool of your system.
%
% You will need at least the core GStreamer runtime and the gstreamer-base set
% of plugins. Modern desktop distros usually already ship with these. A simple
% way to make sure GStreamer is available is to install the "totem" multimedia
% player. This will pull in the needed GStreamer packages as a dependency.
%
% E.g. via 
%
%    sudo apt-get install totem
%
% for Debian distros and derivatives like Ubuntu.
%
% You may need to install some restricted plugin packages to play back all
% common audio and video file formats like MP3 and MP4. Video or movie
% recording with high quality (DivX, H.264) will also require recent versions
% of additional plugin packages containing support for these formats. These
% are usually not installed by default due to licensing and patent clauses in
% place for some territories. You may want to specifically add them to your
% system depending on your format needs.
%
% On Debian 
%
%     sudo apt-get install gstreamer0.10-plugins-bad-multiverse gstreamer0.10-plugins-ugly-multiverse
%
% MS-Windows:
%
% We currently provide GStreamer support on Windows for 32-bit/64-bit Matlab.
%
% In fact, a fully working GStreamer installation is mandatory for
% Psychtoolbox to work on 64-Bit Matlab for Windows, otherwise the Screen()
% function will fail to work at all, regardless if you use multi-media
% functions or not. For 32-Bit Matlab, GStreamer is optional, as long as
% you don't need to use any movie playback/writing or video
% capture/recording.
%
% For 32-Bit Matlab on Windows, you can download and run a standard
% MS-Windows setup.exe installer from:
%
% web: http://code.google.com/p/ossbuild/
%
% For best results, you'd currently (beginning June 2011) choose this
% one or any later version: "GStreamer-WinBuilds-GPL-x86-Beta04-0.10.7.msi"
% The current download link is:
%
% <http://ossbuild.googlecode.com/files/GStreamer-WinBuilds-GPL-x86-Beta04-0.10.7.msi>
%
% If you are a fan of stable versions with more limited functionality, we
% recommend this one: "GStreamer WinBuilds v0.10.6, LGPL (x86)"
%
% However, some features may be disabled or limited with such an older
% version, e.g., playback performance may be reduced, certain video formats
% may not play at all, etc.
%
% For 64-Bit Matlab on Windows, install the 64-Bit GStreamer runtime from
% <http://www.gstreamer.com>. This runtime has higher performance for movie
% playback, but does lack support for video capture/recording as of
% September 2012.
%
%
% Mac OS/X:
%
% For versions of Mac OS/X 10.6 "Snow Leopard" and later, you can download
% and install GStreamer via an easy to use installer from here:
%
% <http://www.gstreamer.com>
%
% This is the most convenient installation method and provides
% high-performance support for movie playback and movie writing (although
% with some limitations on available codecs for movie writing). If you also
% need to capture or record live video from cameras or other live video
% sources, you will need to install GStreamer via the Homebrew package
% manager instead:
%
% Get Homebrew at ...
%
% http://mxcl.github.com/homebrew
%
% ... then install all components via:
%
% brew install gstreamer gst-plugins-base gst-plugins-good gst-plugins-bad gst-plugins-ugly gst-ffmpeg
%
% After a couple minutes of download and compile time, you'll have a fully
% functional GStreamer compiled from source - assuming everything goes
% well, that is.
%
