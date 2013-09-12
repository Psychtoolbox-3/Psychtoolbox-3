% GStreamer - Installation instructions for the GStreamer media framework.
%
% Psychtoolbox uses the GStreamer multi-media framework for all multi-media
% related operations.
%
% All movie playback, movie creation, video capture and video recording
% operations are based on GStreamer. These functions won't work without a
% working GStreamer installation on your system.
%
% You will need at least version 0.10.24 of GStreamer, but we recommend to
% use the latest available versions of the version 0.10.x series, as not
% all Psychtoolbox features are fully supported by version 0.10.24.
% Especially video capture, video recording and movie creation may be
% limited in functionality or performance if you don't use a recent version
% of GStreamer.
%
% GStreamer is a free software / open-source multimedia framework. The
% summary on its webpage tells us:
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
% for many (most?) Linux desktop systems. It is also used on many Linux
% based embedded devices.
%
%
% Installation instructions:
%
%
% GNU/Linux Debian or Ubuntu distributions: If you install Psychtoolbox via
% the package management system of your distribution or via NeuroDebian
% then you don't need to do anything extra.
%
% GNU/Linux other: Any half-way recent Linux distribution will include
% support for GStreamer in its package management system, so you can easily
% install it via the software management tool of your system.
%
% You will need at least the core GStreamer runtime and the gstreamer-base
% set of plugins, but we recommend to install all available plugins for
% full functionality and performance. Modern desktop distros usually
% already ship with these.
%
% You may need to install some restricted plugin packages to play back all
% common audio and video file formats like MP3 and MP4. Video or movie
% recording with high quality (DivX, H.264) will also require recent
% versions of additional plugin packages which contain support for these
% formats. These are usually not installed by default due to licensing and
% patent clauses in place for some territories. You may want to
% specifically add them to your system depending on your format needs.
%
% E.g., on Debian:
%
%     sudo apt-get install gstreamer0.10-plugins-bad-multiverse
%     gstreamer0.10-plugins-ugly-multiverse
%
%
% MS-Windows and Apple MacOSX:
%
% You must install GStreamer if you want to use multi-media functions! If
% you don't intend to use multi-media functions then installation of
% GStreamer is optional.
%
% Download and install the latest GStreamer runtimes from:
%
% <http://www.gstreamer.com>
%
% When the installer asks you which components it should install, go for a
% "full installation" or "complete installation", or in case of a "custom
% install" select all components manually in the checklist if you want
% support for all video formats and all functionality. Without this, many
% popular video formats like H264 video will not play at all, or video
% recording / video capture and similar functions will not work. In fact,
% even our own demos, e.g., SimpleMovieDemo will fail if you don't have all
% codecs installed!

% On Apple MacOSX, optionally you can also install GStreamer version 0.10
% via the Homebrew package manager:
%
% Get Homebrew at ...
%
% http://mxcl.github.com/homebrew
%
% ... then install all components via:
%
% brew install gstreamer gst-plugins-base gst-plugins-good gst-plugins-bad
% gst-plugins-ugly gst-ffmpeg
%
% After a couple minutes of download and compile time, you'll have a fully
% functional GStreamer compiled from source - assuming everything goes
% well, that is.
%
