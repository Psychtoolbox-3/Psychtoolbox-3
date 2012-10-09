% SyncTrouble -- Causes and solutions for synchronization problems.
%
% You most probably arrived at this help page because Psychtoolbox
% aborted with "SYNCHRONIZATION FAILURE" or a similar error message and
% asked you to read this page.
%
% BACKGROUND: Why proper synchronization to retrace is important.
%
% When executing Screen('OpenWindow'), Psychtoolbox executes different
% calibration routines to find out if back- and frontbuffer swaps
% (what Screen('Flip') does) are properly synchronized to the vertical
% retrace signal (also known as VBL) of your display. At the same time,
% it measures the real monitor video refresh interval - the elapsed time
% between two VBL signals. It is crucial for flicker free, tear free,
% properly timed visual stimulus presentation that buffer swaps only
% happen during the VBL period of the display. The VBL (vertical blank) is
% the small gap in time that occurs when the display has updated its last scanline
% and before it starts redrawing its display surface starting at the first
% scanline again. This small gap is a neccessity for CRT displays and it
% is preserved for compatibility reasons or other technical reasons on flat
% panels and video beamers. After issuing the Screen('Flip') command, the
% graphics hardware finalizes all pending drawing- and imageprocessing
% operations in the backbuffer of an onscreen window to make sure that the
% final stimulus image is ready in the backbuffer for presentation. Then
% it waits for onset of the next VBL interval before flipp'ing the back-
% and frontbuffer of the onscreen window: The previous backbuffer with your
% newly drawn stimulus becomes the frontbuffer, so it will get scanned out
% and displayed to the subject, starting with the next refresh cycle of your
% display device and on all consecutive refresh cycles, until you draw a new
% stimulus to the onscreen window and update the display again via Screen('Flip').
%
% On a properly working system, this double buffer swap happens in less than
% a microsecond, synchronized to VBL onset with an accuracy of better than a
% microsecond. All change of visual content therefore only happens during
% the VBL period when the display is not updating, thereby avoiding any kind of
% visual flicker or tearing that would be caused by a mixup of an old stimulus and
% a new (incompletely drawn) stimulus when changing image content during the
% scanout cycle of the display. The exact point in time when this
% buffer swap happened, is returned as timestamp by the Screen('Flip') command.
% It is the most well defined timestamp of visual stimulus onset, and it allows to
% define stimulus onset of future stims relative to this accurate baseline,
% using the 'when' argument of Screen('Flip').
%
% Without proper synchronization, you would see very strong visual flicker and
% tearing artifacts in animated (movie / moving) stimuli, you would not have any
% well defined stimulus onset for sequences of static stimuli or rapid stimulus
% presentation, and no means of synchronizing visual stimulus presentation to any
% external stimulation- or acquisition devices like fMRI, EEG, sound, ... You also
% would not have any accurate way of getting a stimulus onset timestamp.
%
% However, if you have very special needs, you can disable either Matlabs / Octaves
% synchronization of execution to the vertical retrace or you can disable synchronization
% of stimulus onset to the vertical retrace completely by setting the 'dontsync' flag
% of Screen('Flip') accordingly.
%
% For more infos about tearing, see Wikipedia articles about "Tearing", "Double buffering",
% "Vertical Synchronization" and the info pages on www.psychtoolbox.org
%
% TESTS: How Psychtoolbox tests for proper synchronization to retrace.
%
% After opening an onscreen window, Psychtoolbox executes a measurement loop,
% where it issues Screen('Flip') commands and measures the time elapsed between
% execution of two consecutive flip commands, getting one refresh sample per
% loop iteration. Each sample is checked for validity: Duration must be longer than
% 4 milliseconds and shorter than 40 milliseconds, because we assume that none of
% the available display devices updates slower than 25 Hz or faster than 250 Hz. Each
% sample is also tested against the expected value provided by the operating system, e.g.,
% if the operating system reports a nominal refresh rate of 100 Hz, then a sample should have
% a duration of roughly 1000 ms / 100 Hz == 10 milliseconds. We accept any sample in a
% range of +/- 20% around this expected value as valid, because timing jitter present in
% any computer system can cause some deviation from the expected value. Samples that don't
% pass this basic test are rejected. Valid samples are used to update a mean value, standard
% deviation of the mean is also calculated: The measurement loop ends when at least 50 valid
% samples have been taken and the standard deviation from the mean is less than 200 microseconds.
% If it is not possible to satisfy this criteria during a five second measurement interval, then the
% calibration is aborted and repeated for up to three times. Failure to get a valid measurement
% during up to three calibration runs is indicating massive timing problems or the inability
% of the gfx-hardware to properly synchronize buffer swaps to the vertical retrace. This leads
% to abortion with the "SYNCHRONIZATION FAILURE" error message.
% Assuming that this calibration loop did provide a valid mean measurement of monitor refresh,
% the value is checked against the value reported by the operating system and - on MacOS-X - against
% the result of an independent measurement loop that uses direct queries of rasterbeam positions
% to measure the monitor refresh interval. Only if all available measurements yield similar results,
% the test is finally rated as PASSED, Psychtoolbox continues execution and the computed monitor
% refresh interval is used internally for all built-in timing checks and for properly timed
% stimulus presentation.
%
% REASONS FOR FAILING THE SYNC TESTS AND HOW TO FIX THEM:
%
% There are multiple classes of possible causes for sync failure. Work down this
% list of causes and solutions down until your problem is resolved or you
% hit the bottom of the list:
%
% 1. Wrong configuration settings: This usually only affects MS-Windows
% systems, where the display settings control panel for your graphics card
% allows to customize a couple of graphics driver parameters. Some of these
% settings can cause sync failure if they are wrong:
%
% -> Make sure the "Synchronize bufferswaps to the vertical retrace" option
% is set to "Application controlled" or "Application controlled, default to
% on". The wording of the option differs between different graphics cards,
% search for something like that. Examples of other names: "Wait for
% vertical sync", "Wait for vertical refresh" ...
% If this setting is forced to off and *not* application controlled, then
% the sync tests will fail because the hardware doesn't synchronize its
% image onset (bufferswap) to the video refresh cycle of your display.
%
% -> Make sure the "Triple buffering" setting is off, or if you can select
% some "Multibuffering" setting, that it is set to "double buffering" or
% "wait for 1 video refresh" or "swap every refresh". This option may not
% exist, but if it does, any other setting will cause the sync tests to
% possibly succeed, but later stimulus onset timestamping to fail with
% errors.
%
% -> If there is an option "Buffer swap mode" or "Bufferswap strategy", it
% should be set to "Auto select" or "Page flipping" or "Exchange buffers".
% The so called "Copy buffers" or "Blitting" option would result in lower
% performance and inaccurate timing.
%
% -> On dual/multi display setups MS-Windows allows you to assign one
% monitor the role of the "primary monitor" or "primary display". It is
% important that the display device which you use for stimulus presentation
% is the "primary display", otherwise random things may go wrong wrt. sync
% tests and timing.
%
% -> If you have the choice to set your multi-monitor configuration to
% either "dual display mode"/"dual display performance mode"/"separate
% displays" or instead to "extended desktop mode" or "horizontal spanning",
% you should choose "extended desktop mode" or "horizontal spanning" modes
% for best timing and stimulus quality. Please note that this choice
% doesn't exist anymore on Windows-Vista and later.
%
% -> On all operating systems in dual display or multi display mode it is
% important that you configure both displays for exactly the same color
% depths, resolution and refresh rate. If there is some option you can
% choose for "genlocked modes" or "genlocked modes only", choose or enable
% that one. Failing to configure dual display setups like this will cause
% massive timing problems or tearing artifacts on one of the display if you
% do dual display stimulation. It may also cause failures in timetamping.
%
%
% 2. Temporary timing glitches or system malfunction: It may help to
% restart Matlab/Octave, or to reboot your machine. Sometimes this resolves
% intermittent problems on your system, especially after the system was
% running without reboot for a long time, on high load, or if display
% settings or display configuration has been changed frequently.
%
%
% 3. Driver bugs: Many graphics card device drivers have bugs that cause
% synchronization to fail. If none of the above steps resolves your
% problems, check the website of your computer vendor or graphics card
% vendor, or use the "Check for driver updates" function of some operating
% systems to find out if new, more recent graphics drivers have been
% released for your graphics card. If so, update to them. A tremendeously
% large number of problems can be resolved by a simple driver update!
%
% 4. Driver/Hardware limitations:
%
% Many systems can't provide reliable research grade timing if you don't
% display your stimuli in fullscreen windows, but use windowed mode
% instead. This can lead to sync failures, problems with timestamping and
% other performance problems. Only use non-fullscreen windows for
% development, debugging and leisure, not for running your studies!
%
% Some systems have serious problems if more than one graphics card is
% connected and enabled on the computer. They only work well in
% single-display mode or dual display mode from a single dual-output
% graphics card.
%
% Microsoft Windows Vista and Windows-7 may provide poor performance on
% dual display setups if you present on both displays simultaneously,
% although your mileage may vary widely depending on exact setup.
%
% On Vista and Windows-7, you may run into drastic timing and performance
% problems if the stimulus presentation window loses the "keyboard focus".
% The window with the "keyboard focus" is the one which is in the
% foreground (in front of all other windows), has its titlebar highlighted
% instead of shaded (assuming it has a titlebar) and receives all keyboard
% input, i.e., key presses. Therefore we assign "keyboard focus" to our
% onscreen windows automatically. However, if the user clicks into windows
% other than our window with the mouse, or onto the desktop background, or
% uses key combos like ALT+TAB or Windows+TAB to switch between windows,
% our window will "lose" the keyboard focus and severe timing and
% performance problems may occur. Obviously if any window on the screen is
% highlighted, this means it *has stolen* the keyboard focus from our
% window. This weird keyboard focus problem is an unfortunate design
% decision (or rather design flaw) of the Windows Vista/Win-7 graphics
% subsystem. There isn't anything we or the graphics cards vendors could do
% about it, so you'll have to accept it and work-around it. Of course this
% becomes mostly a problem on dual-display setups where one display shows
% the desktop and GUI, so avoid such configurations if you can.
%
% Further examples:
%
% On most systems, graphics adapters are only capable of properly timed
% bufferswaps and proper visual stimulus onset timestamps syncing to
% retrace if the onscreen window is a full-screen window. Synchronization
% fails if the onscreen window only covers part of the screen (i.e., when
% providing a 'rect' argument to Screen('OpenWindow') other than the
% default full-screen rect). Solution is to only use fullscreen windows for
% stimulus presentation. On Windows, Linux and MacOS-X, some graphics cards
% (e.g., many - if not all - mobile graphics cards built into Laptops) are
% only capable of synchronizing to the retrace of one display. On a single
% display setup, this will simply work. On a dual display setup, e.g.,
% Laptop connected to external video beamer or CRT, the driver/hardware can
% sync to the wrong output device. A simple, although inconvenient solution
% is to disable the internal flat panel of a Laptop while running your
% study, so the hardware is guaranteed to sync to the external display.
% Depending on the hardware it may also help to try dual display output
% with either: Non-mirror mode, running both displays at different refresh
% rates, mirror mode running both displays at different rates, mirror mode
% running both displays at exactly the same resolution, color depth and
% refresh rate. You'll have to try, as it has been found to be highly
% dependent on hardware, driver and operating system, which combinations
% work and which don't.
%
% 5. Graphics system overload: If you ask too much from your poor graphics
% hardware, the system may enter a state where the electronics is not
% capable of performing drawing operations in hardware, either because it
% runs out of video memory ressources, or because it is lacking the
% neccessary features. In that case, some drivers (e.g., on Microsoft
% Windows or MacOS-X) may activate a software rendering fallback-path: The
% graphics engine is switched off, all rendering is performed by slow
% software in system memory on the cpu and the final image is copied to the
% onscreen framebuffer. While this produces visually correct stimuli,
% presentation timing is completely screwed and not synchronized to the
% monitors refresh at all. On Microsoft Windows, Psychtoolbox will detect
% this case and output some warnings to the Matlab window.
%
% Possible causes of such an overload: Running with Anti-Aliasing enabled
% at a setting that is too high for the given screen resolution (see 'help
% AntiAliasing'), or running at a display resolution that is too high,
% given the amount of video memory installed on your graphics adapter.
% There may be other cases, although we didn't encounter any of them up to
% now. The same could happen if you run a dual display setup that is not
% switched to mirror-mode (or clone mode), so you take up twice the amount
% of video memory for two separate framebuffers.
%
% Troubleshooting: Try lower display resolutions and multisampling levels,
% switch dual display setups into mirror-mode if possible, or buy a
% graphics adapter with more onboard memory.
%
% 6. General system overload: If you run too many applications on your
% system in parallel to your Psychtoolbox+Matlab/Octave session, then these
% applications may cause significant timing jitter in your system, so the
% execution of Psychtoolbox - and its measurement loops - becomes
% non-deterministic up to the point of being unuseable.
%
% Troubleshooting: Quit and disable all applications and services not
% needed for your study, then retry. The usual suspects are: Virus
% scanners, applications accessing the network or the harddiscs,
% applications like iTunes, system software update...
%
% 7. Bad drivers or hardware in your system that interferes with general
% system timing: This is difficult to diagnose. At least on MS-Windows, you
% can download a free tool "dpclat.exe" from the internet (Use Google to
% find it). If you run it, it will tell you if there are potential problems
% with your systems timing and give hints on how to resolve them.
%
% 8. Other: Search the FAQ pages on the www.psychtoolbox.org Wiki and (via
% Google search) the Psychtoolbox forum for other problems and solutions.
%
% 9. If everything else fails, post on the forum for help, but read our
% instructions on how to ask questions on the forum properly. You can find
% these instructions on the "Forum" and "Bugs" pages of our Wiki. If we
% find that you didn't read the instructions and you're basically wasting
% our time due to your omissions, we will simply ignore your request for
% help.
%
%
% HOW TO OVERRIDE THE SYNC TESTS:
%
% That all said, there may be occassions where you do not care about
% perfect sync to retrace or millisecond accurate stimulus presentation
% timing, but you do care about listening to iTunes or getting your
% stimulus running quickly, e.g., during development and debugging of your
% experiment or when showing a quick & dirty online demo of your stimulus
% during a presentation. In these situations you can add the command
% Screen('Preference','SkipSyncTests', 1); at the top of your script,
% before the first call to Screen('OpenWindow'). This will shorten the
% maximum duration of the sync tests to 3 seconds worst case and it will
% force Psychtoolbox to continue with execution of your script, even if the
% sync tests failed completely. Psychtoolbox will still print error
% messages to the Matlab/Octave command window and it will nag about the
% issue by showing the red flashing warning sign for one second. You can
% disable all visual alerts via Screen('Preference','VisualDebugLevel', 0);
% You can disable all output to the command window via Screen('Preference',
% 'SuppressAllWarnings', 1);
%
% If your graphics system basically works, but your computer has just very
% noisy timing you can adjust the threshold settings we use for our tests
% via the setting:
%
% Screen('Preference','SyncTestSettings' [, maxStddev=0.001 secs][, minSamples=50][, maxDeviation=0.1][, maxDuration=5 secs]);
%
% 'maxStddev' selects the amount of tolerable noisyness, the standard
% deviation of measured timing samples from the computed mean. We default
% to 0.001, ie., 1 msec.
%
% 'minSamples' controls the minimum amount of valid measurements to be
% taken for successfull tests: We require at least 50 valid samples by
% default.
%
% 'maxDeviation' sets a tolerance threshold for the maximum percentual
% deviation of the measured video refresh interval duration from the
% duration suggested by the operating system (the nominal value). Our
% default setting of 0.1 allows for +/- 10% of tolerance between
% measurement and expectation before we fail our tests.
%
% 'maxDuration' Controls the maximum duration of a single test run in
% seconds. We default to 5 seconds per run, with 3 repetitions if
% neccessary. A well working system will complete the tests in less than 1
% second though.
%
% Empirically we've found that especially Microsoft Windows Vista and
% Windows-7 may need some tweaking of these parameters, as some of those
% setups do have rather noisy timing.
%
%
%
% MORE WAYS TO TEST:
%
% The script VBLSyncTest() allows you to assess the timing of Psychtoolbox
% on your specific setup in a variety of conditions. It expects many
% parameters and displays a couple of plots at the end, so there is no way
% around reading the 'help VBLSyncTest' if you want to use it.
%
% The script PerceptualVBLSyncTest() shows some flickering stimulus on the
% screen and allows you to assess visually, if synchronization works
% properly.
%
% Both tests are for the really cautious: The built-in test of
% Screen('OpenWindow') should be able to catch about 99% of all conceivable
% synchronization problems.
%
% MORE READING: See the help for 'help MirrorMode' and 'help
% BeampositionQueries' for more info about display issues.
%
%
% LINUX specific tips:
%
% 1. Just as on all other operating systems, timed visual stimulus onset
% and stimulus onset timestamping is not reliable on regular (non
% fullscreen) windows, windows that don't cover the complete desktop of a
% Psychtoolbox screen (also known as X-Screen) or transparent windows,
% e.g., when the PsychDebugWindowConfiguration() command was used. Use of
% rotated display output (90, 180, 270 degrees etc.) will also prevent
% proper timing.
%
% Your windows must be non-transparent, decoration/borderless, fullscreen
% and cover a complete X-Window system screen. On a multi-display setup
% that means that either your window must cover all connected displays, or
% you need to setup separate X-Screens in the graphics driver control panel
% GUI or the /etc/X11/xorg.conf file for different displays, so that all
% stimulus displays are grouped in one (or multiple) X-Screen which are
% fully covered by your PTB onscreen window, and other displays, e.g.,
% operator GUI displays, are grouped into a different X-Screen.
%
% 2. A major source of timing trouble can be 3D desktop compositors. Either
% use a user interface that doesn't employ a desktop compositor, or
% configure your desktop compositor of choice to
% "unredirect_fullscreen_windows", which will disable the compositor for
% Psychtoolbox fullscreen windows. With most versions of compiz, the
% command PsychGPUControl('FullScreenWindowDisablesCompositor', 1); can do
% this setup step for you. If you use a different desktop compositor, and
% on a few versions of compiz, the command doesn't yet work, so you need to
% check the manuals/help of your system on how to enable the option
% "unredirect_fullscreen_windows" manually. If in doubt, just use a desktop
% session without 3D compositor for running the actual data collection of
% your studies. Examples of desktops which use a compositor: Ubuntu Unity,
% GNOME-3, GNOME-2, KWin. Examples which don't use a compositor: GNOME-2
% classic, Mate desktop, XFCE.
%
% If you use the NVidia or AMD proprietary graphics drivers, frequent
% synchronization failures, or other sync related warnings, or unsteady,
% irregularly timed flicker during PerceptualVBLSyncTest indicate that a
% desktop compositor is in use. On Linux with the free and open-source
% graphics drivers radeon, nouveau or intel, PTB will output warnings about
% non-pageflipped flips in such a case.
%
% 3. Another reason for timestamping problems can be the use of
% triple-buffering if your system has an X-Server older than version 1.12,
% or if you have an old system with an outdated graphics driver.
% Psychtoolbox would warn you at least when the free graphics drivers are
% in use, about some error in "PsychOSGetSwapCompletionTimestamp" and some
% system configuration problems.
% On such setups, triple-buffering can be disabled with driver specific
% options in xorg.conf:
%
% On Intel graphics drivers, add the options:
%
% Option "TripleBuffer"    "off"
% Option "SwapbuffersWait"  "on"
%
% On the nouveau driver for NVidia cards, add:
% Option "SwapLimit" "1"
% Option "GLXVBlank" "on"
%
% On Linux with the open-source intel, nouveau and radeon graphics drivers,
% the script OMLBasicTest allows some additional correctness checks.
%

% History:
% 17.06.2006 written (MK).
% 09.10.2012 Add Linux-specific section (MK).
