% DisplayOutputMappings -- How they work, how to resolve problems.
%
% You probably read this text, because Screen() gave you some error or
% warning message about problems with beamposition queries or other
% low-level functions and directed you to this document for
% troubleshooting.
%
% So, what happened here?
%
% For various low-level operations, e.g., beamposition queries for
% timestamping, control of high precision framebuffers and high precision
% display devices, and for certain stereo display modes, Screen() needs to
% know which video output connector on your graphics card is connected to
% which electronic display scanout engine (a so called "crtc") in the
% graphics card, because it needs to exercise low-level control over the
% crtc associated with a specific display. The actual wiring between crtc's
% and display connectors however is not fixed on modern graphics cards, but
% it is flexible. The wiring is controlled by electronic programmable
% switches, so called multiplexers ("MUXers"), which are configured by the
% operating system and graphics card driver. Depending on the specific
% hardware, operating system and display setup, the wiring can change
% dynamically.
%
% You can find more technical details at the following links if you are
% interested: <http://www.botchco.com/agd5f/?p=51> for AMD/ATI hardware,
% and <http://virtuousgeek.org/blog/index.php/jbarnes/2011/01/> for Intel
% hardware.
%
%
% Sometimes, due to operating system or graphics driver bugs, or system
% misconfiguration, the operating system "lies" to Screen() about the true
% wiring. On some systems, the operating system can't tell Screen() about the
% wiring and Screen() has to make an educated guess based on some
% heuristics - a guess which can go wrong.
%
% Long story short, Screen() can be wrong about which crtc to access for a
% given display, causing the kind of malfunctions, warnings and errors
% messages that brought you to this help text. To fix such problems you
% need to help Screen() in one of multiple different ways:
%
% 1. On MS-Windows, try to update your graphics driver in the hope that
%    this fixes the problem. Don't forget to reboot!
%
% 2. If this doesn't help, or if you are on Linux or OS/X, try to replug
%    the displays into different video output connectors on your computer.
%    E.g., on a single monitor setup, plug the monitor into the other video
%    output connector. On a dual-display setup, exchange which monitor is
%    plugged into which connector etc. This replugging will make reality
%    match the expectations of Screen(). Don't forget to restart Matlab or
%    Octave after the change.
%
% 3. If 2. doesn't help or is infeasible or problematic, you can also tell
%    Screen() about the true wiring by adding the command
%    Screen('Preference', 'ScreenToHead', screen, head, crtc[, rank]); to the
%    top of your script, before other Screen() commmands:
%
%    On OS/X or Windows, Screen('Preference', 'ScreenToHead', 0, 1, 1);
%    would tell Screen() that the Psychtoolbox screen with screenid 0 is
%    not connected to video output 0 and crtc 0 (as would be the default),
%    but to video output 1 and crtc 1.
%
%    On Linux the same logic applies. However, on Linux, multiple video
%    outputs (and thereby display monitors) can be connected to one single
%    Psychtoolbox screenid (aka X-Window system X-Screen) to allow for more
%    flexibility than on the other systems. For this reason an additional
%    'rank' parameter controls which of multiple possible outputs per
%    screen is remapped. The default 'rank' of 0 refers to the primary
%    display output, the one which is used for stimulus onset timestamping
%    or framerate queries. It may therefore be neccessary to play with the
%    'rank' parameter as well on multi-display setups with multiple
%    monitors per Psychtoolbox screen.
%
% 4. If you have a multi-GPU setup, ie., multiple graphics cards installed
%    and active at the same time, then Screen() low-level functions may not
%    work. More precisely, they will not work at all on MS-Windows. On
%    Apple OS/X or on Linux you can make them work on exactly one GPU, the
%    other GPU's will be ignored. By default, the first GPU in the system
%    is chosen. You can override the choice by executing the command
%    PsychTweak('UseGPUIndex', gpuidx); *before* other calls to Screen() to
%    instruct Screen() to use the GPU with index gpuidx. If you call the
%    command after Screen() has been called already, it will usually not
%    pick up your new setting! E.g., PsychTweak('UseGPUIndex', 1); would
%    instruct Screen() to use the GPU with index 1, which is the 2nd GPU in
%    the system, as numbering starts with zero for the 1st GPU.
%
%    If you use PsychTweak() to select a different GPU than the first one,
%    you often also need to perform the actions explained in step 3 to
%    manually define which GPU output maps to which display connector, as
%    the heuristic will most often result in wrong assignments for other
%    GPU's than the first one (aka GPU zero).
%
%    Also remember that all your visual stimulation (== "experiment")
%    displays should be connected to the chosen GPU, not one of the other
%    GPUs! Other GPUs should usually only be used for non-Psychtoolbox
%    purposes, e.g., operator displays like the Matlab or Octave window,
%    etc. While it is possible to open Screen() onscreen windows on the
%    other GPUs as well, visual stimulus onset timestamping and other
%    special functions won't work there and you'll likely see a bunch of
%    timing error warnings ("Synchronization trouble") for those windows,
%    which indicate that visual timing on those windows will be broken or
%    highly unreliable.
%
%    For an example case of such a multi-gpu setup under Linux, and the
%    procedure that resolved problems and made timestamping work, see
%    Psychtoolbox forum message #16208:
%
%    http://tech.groups.yahoo.com/group/psychtoolbox/message/16208
%
% Good luck with troubleshooting!
%
