% ConserveVRAMSettings: Workaround for flawed hardware and drivers
%
% The command Screen('Preference', 'ConserveVRAM', mode); can be used to
% enable a couple of special work-arounds inside Screen to work around
% broken operating systems, graphics drivers or graphics hardware, or to
% work around ressource limitations of graphics hardware.
%
% You define the requested workaround by setting the parameter 'mode' to a
% sum of the following values:
%
% Allowerd summands (flags) for 'mode' and their effect:
%
% 1 == kPsychDisableAUXBuffers: A setting of 1 asks Psychtoolbox to not
% allocate any OpenGL AUXiliary buffers when opening a new onscreen window.
% AUX buffers are only needed if you want to run the Screen('Flip') command
% with the optional argument 'dontclear = 1' and you are not using the
% imaging pipeline, or if you want to use stereomode 2 or 3 without using
% the imaging pipeline. If you do use the imaging pipeline or don't use any
% of the above, there's no need for AUX buffers.
%
% This setting is mostly meant to save a bit of VRAM on graphics hardware
% that only has very small amounts of VRAM, e.g., only 16 MB or 8 MB VRAM.
%
%
% 2 == kPsychDontCacheTextures: A setting of 2 asks Psychtoolbox not to
% cache used textures in the graphics hardware local VRAM memory. This will
% save some VRAM memory at the expense of lower drawing performance. Only
% useful on gfx-hardware with low amounts of VRAM and only works on
% MacOS/X. The flag is silently ignored on Windows and Linux.
%
%
% 4 == kPsychOverrideWglChoosePixelformat: This is a workaround for broken
% MS-Windows graphics drivers: Ask Screen to not use the
% wglChoosePixelFormat() command when creating a new onscreen window. This
% can prevent crashes on such broken setups, but it will disable OpenGL
% multisampling for anti-aliasing, ie., the 'multisample' parameter of
% Screen('OpenWindow') will be ignored. In the future, other special
% capabilities will be disabled as well.
%
%
% 8 == kPsychDisableContextIsolation: This is a workaround for broken
% MS-Windows graphics drivers: Do not create separate isolated OpenGL
% rendering contexts for Screen and MOGL when using low level OpenGL 3D
% graphics commands with OpenGL for Matlab. This prevents crashes on broken
% setups, but debugging of your own 3D code may become much harder. Its
% better to upgrade to the latest fixed drivers.
% Before you try this setting 8, first try if the setting 256 (see below)
% fixes the problem for you. That is a softer approach - If it works for
% you then you won't lose any important functionality!
%
%
% 16 == kPsychDontAttachStencilToFBO: Do not attach stencil buffer
% attachments to OpenGL framebuffer objects when using OpenGL 3D graphics
% in conjunction with the Psychtoolbox imaging pipeline. This is again a
% workaround for some broken MS-Windows graphics drivers to make the 3D +
% imaging combo work at least when no stencil buffer is needed.
%
%
% 32 == kPsychDontShareContextRessources: Do not share ressources between
% different onscreen windows. Usually you want PTB to share all ressources
% like offscreen windows, textures and GLSL shaders among all open onscreen
% windows. If that causes trouble for some weird reason, you can prevent
% automatic sharing with this flag.
%
%
% 64 == kPsychUseSoftwareRenderer: Request use of a software implemented
% renderer instead of the GPU hardware renderer. This request is silently
% ignored if your platform doesn't support software rendering. Currently
% only MacOS/X 10.4 and later in windowed mode (i.e. not fullscreen)
% supports this via the Apple floating point renderer. Mostly useful for
% testing and debugging of scripts that need floating point support on
% hardware that doesn't support this. Not generally useful for production
% use.
%
% On MS-Windows this allows to use the Microsoft GDI software renderer,
% ie., to not abort if that renderer is detected. Normally Screen() would
% abort when detecting the GDI renderer.
%
%
% 128 == kPsychEnforceForegroundWindow: Request application of the Windows
% GDI calls SetForegroundWindow() and SetFocus() on each created onscreen
% window on MS-Windows. This may improve reliabilty of onscreen windows
% staying in front of all other windows, but is incompatible with the use
% of GetChar, CharAvail and ListenChar, so it must be requested with this
% flag. These calls are unfortunately absolutely crucial on MS-Vista and
% later operating systems to guarantee artifact free (tear-free) and timing
% accuratestimulus onset and robust and accurate stimulus onset
% timestamping. Therefore they are automatically applied to all fullscreen
% windows on Windows Vista and later operating systems. See the option
% kPsychPreventForegroundWindow to forcefully disable/prevent use of these
% options, if use of GetChar() et al. is more important than artifact free
% stimulus presentation.
%
%
% 256 == kPsychUseWindowsContextSharingWorkaround1
% On MS-Windows, skip a few not too essential setup steps when creating a
% userspace OpenGL rendering context for 3D mode. This is a "soft" version
% of kPsychDisableContextIsolation -- Less intrusive as it doesn't disable
% context isolation completely, but only a subset. May be able to
% work-around an NVidia driver bug reported in March 2008 on GF8xxx series.
%
%
% 512 == kPsychAvoidCPUGPUSync: Avoid any internal calls (if possible) that
% could cause a synchronization of the CPU and GPU. Synchronization is a
% potentially expensive operation that can degrade performance in certain
% circumstances. Its often needed for error checking. Setting this flag may
% give you a speedup on certain operations, but at the cost of reduced
% error checking and error handling: Error conditions detected otherwise
% may silently slip through and cause mysterious malfunctions or stimulus
% corruption without PTB noticing this or providing any troubleshooting
% tips. The usefulness of this flag highly depends on your graphics
% hardware, driver and operating system. It may give a large speedup, or no
% speedup at all, but it will always reduce robustness!
%
%
% 1024 == kPsychTextureUploadFormatOverride
% Tell PTB to use the opposite texture format of what its auto-detection
% thinks is optimal. Screen contains code to auto-detect certain type of
% graphics chips with broken drivers and tries to work-around them by
% chosing different parameters for fast texture creation in certains
% circumstances. In case those vendors should ever fix their drivers and
% thereby the built-in workaround becoming invalid, this allows to override
% PTB's choice. This is mostly to work around broken ATI drivers on
% MS-Windows which cause miserable texture creation performance with the
% standard optimized settings.
%
%
% 2048 == kPsychAvoidFramebufferBlitIfPossible
% Tell PTB to not use the EXT_framebuffer_blit extension if a lower-speed
% workaround solution exists. This will mostly affect the operation of
% Screen('CopyWindow') when the imaging pipeline is active. Normally a more
% flexible, capable, faster method would be used, unless you set this flag
% to fall back to the old solution.
%
%
% 4096 == kPsychUseBeampositionQueryWorkaround
% Tell PTB to always use the workaround for broken beamposition queries in
% VBL on MS-Windows, even if the automatic startup test does not detect any
% problems. This for rare cases where the test fails to detect broken
% setups.
%
%
% 8192 == kPsychUseAGLForFullscreenWindows
% Tell PTB on Mac OS/X to always use the Cocoa/NSOpenGL API for OpenGL
% system setup, even if the requested onscreen window is a fullscreen
% window. Normally PTB would use the CGL API for fullscreen windows for
% higher efficiency. This setting is automatically implicitely applied if
% PTB is running on OSX version 10.8 "Mountain Lion" or later, to work
% around various hilarious graphics driver bugs.
%
%
% 16384 == kPsychUseCompositorForFullscreenWindows
% Tell PTB to use a compositing window manager for stimulus display if such
% a desktop compositor is supported on your operating system. Currently
% this flags affects operation on MacOS/X and on Microsoft Windows Vista
% and later versions of Windows.
%
% On Windows Vista and Windows-7, it will enforce use of the Windows Aero
% desktop compositor (aka DWM or Desktop Window Manager). Accuracy and
% reliability of visual stimulus onset timing and the accuracy and
% reliability of stimulus onset timestamps will be greatly reduced in this
% mode - up to the point of being completely useless for timed stimulus
% presentation. The mode may however be useful for debugging and code
% development as a convenience feature. By default, PTB disables the DWM as
% soon as a fullscreen window is opened, unless the
% PsychDebugWindowConfiguration() function was used to switch to debug
% mode.
%
% On Mac OS/X this will cause PTB to always use the Cocoa API for OpenGL
% system setup, and to always use Quartz composited regular windows instead
% of fullscreen contexts, even if the requested onscreen window is a
% fullscreen window. Normally PTB would use the CGL API for fullscreen
% windows, but this was broken for dual-display operations on some systems
% running OS/X Leopard 10.5 with all modern NVidia cards. Stimulus onset
% timing, animations and timestamping precision will be horrible in this
% mode. Only suitable for presentation of mostly static stimuli with no
% requirements for frame-accurate timing.
%
%
% 32768 == kPsychBusyWaitForVBLBeforeBufferSwapRequest
% If Screen('Flip') in sync with vertical retrace is requested and
% beamposition queries are supported, use a busy-waiting, high cpu load
% spin-wait loop to wait for onset of vertical blank interval (VBL) before
% submitting doublebuffer swaprequests to the GPU. This is meant as a
% last-ressort workaround for GPU's with severely broken sync-to-VBL
% support. The only known current example is the Apple Leopard operatings
% system when used with NVidia Geforce 8000 GPU's or later in
% frame-sequential stereo mode. This will create a very high cpu load and
% may have negative side effects on system timing. Use as last resort!
%
%
% 65536 (= 2^16) == kPsychDontUseNativeBeamposQuery
% Do not use operating system native beamposition queries, but try to use
% own mechanism, or none at all. This to work around bugs in OS native
% beamposition query mechanisms, e.g., Leopard 10.5.7 + ATI GPU's.
%
%
% 131072 (= 2^17) == kPsychDisableAeroDWM 
% Disable the Aero DWM desktop composition manager on Windows Vista and
% later. By default, Psychtoolbox will try to keep the WDM running,
% actually enforce use of it, as this provides better timing behaviour. If
% you find otherwise on your setup or have special needs, specify this flag
% to force the DWM off.
%
%
% 262144 (= 2^18) == kPsychPreventForegroundWindow
% Prevent calls to the Windows GDI functions SetForegroundWindow() and
% SetFocus() on each created fullscreen onscreen window on MS-Windows.
% These calls would prevent proper use of GetChar(), but are needed on
% MS-Vista and later for proper visual stimulus onset timing. With this
% flag you can make a conscious decision between proper stimulus display
% and use of GetChar.
%
%
% 524288 (= 2^19) == kPsychDisableOpenMLScheduling
% Disable use of OpenML scheduling for Screen('Flip') bufferswaps. OpenML
% is currently supported on some recent versions of GNU/Linux with certain
% graphics cards and drivers (e.g., Free graphics stack on Ubuntu 10.10 and
% later with modern Intel and ATI/AMD GPU's and XOrg Servers 1.8.2, 1.9.x and
% later, Linux kernel 2.6.35 and later). PTB will use OpenML for scheduling
% and timestamping of visual stimulus onset if it detects a Linux system
% with support for OpenML. The kPsychDisableOpenMLScheduling flag will
% forcefully disable use of OpenML, e.g., for debugging/testing purpose.
%
%
% 1048576 (= 2^20) == kPsychBypassLUTFor10BitFramebuffer
% If a 30 bpp, 10 bpc native 10 bit framebuffer is requested and
% Psychtoolbox is executing on Linux (as superuser) or OS/X (with the
% PsychtoolboxKernelDriver loaded), then apply the 10 bit LUT bypass enable
% hack even on graphics cards where this should not be required, e.g., on
% the FireGL or FirePro cards from ATI/AMD. This as a workaround for broken
% ATI/AMD graphics drivers which are able to configure a 10 bit framebuffer
% and scanout, but fail to setup the LUT's properly.
%
%
% 2097152 (= 2^21) == kPsychEnforce10BitFramebufferHack
% Use 10 bpc framebuffer hack even if PTB thinks it is not needed or
% appropriate. This implies kPsychBypassLUTFor10BitFramebuffer, and the
% same conditions must be met for it to possibly work. This can be used to
% enforce 10 bpc mode on FireGL/FirePro GPU's with broken drivers. It
% overrides all our safeguards and may end in a hard machine crash is used
% on the wrong setup. Try to use the regular way of enabling this in ATI's
% Catalyst Control Center application if possible.
%
%
% 4194304 (= 2^22) == kPsychIgnoreNominalFramerate
% Do not use the nominal video refresh rate of a screen as reported by the
% operating system for internal calibrations and tests. Return zero instead
% of this rate in calls to Screen('Framerate') or Screen('NominalFramerate').
% This to work around broken or problematic video refresh reporting mechanisms.
%
%
% 2^23 == kPsychUseOldStyleAsyncFlips
% Do not use the enhanced Screen('AsyncFlipBegin') implementation which
% allows for more parallelism between your code and pending async flips.
% There is no reason to use this flag except for benchmarking by PTB
% developers. Or, of course, if you should happen to have an operating
% system + graphics driver + GPU combo which performs much worse with the
% new method than with the old one.
%
%
% 2^24 == kPsychDontAutoEnableImagingPipeline
% Do not automatically enable support for fast offscreen windows on
% graphics cards (GPU's) that support this. Do not automatically enable the
% full Psychtoolbox image processing pipeline on supported GPU's for stereo
% display modes which would benefit from it.
%
% On Psychtoolbox versions released before the year 2012, you needed to
% enable those two features manually. 2012+ PTB's will auto-enable on
% suitable hardware if this promises improvements in performance,
% flexibility or robustness. In the unlikely case you should encounter
% problems with this behaviour due to graphics driver or operating system
% bugs, you can revert to the old opt-in behaviour: The PsychImaging()
% command then allows you to enable those features manually and separately,
% just as on pre 2012 PTB's.
%
%
% 2^25 == kPsychOldStyleOverrideRedirect
% Use the old-school method of setting the override_redirect property
% of X11 Windows on Linux, as it was done until end 2012. This just as
% a safe-guard in case somebody runs a very exotic or old Linux setup,
% where the new method doesn't work.
%
%
% 2^26 == kPsychForceUseNativeBeamposQuery
% Force use of the OSX builtin beamposition query mechanism, even if the
% PsychtoolboxKernelDriver is installed and would provide better results.
% This is mostly for internal testing and benchmarking of Psychtoolbox. As
% of April 2013 and MacOSX 10.8 "Mountain Lion", our own implementation is
% superior to Apple's broken implementation, so we default to use of our
% own implementation whenever possible. The OSX native implementation is
% only supported on NVidia gpu's anymore, scheduled for removal in a future
% OSX version (marked as deprecated in Apple developer documentation!),
% always slightly broken and more noisy, imprecise and higher overhead than
% our own implementation. At least on OSX 10.8 with GeForce-600 gpu's, the
% OSX implementation is completely broken for anything but analog VGA CRT
% displays, ie., it delivers bogus results for any kind of internal or
% external digital lcd flat panel!
%
%
% 2^27 == kPsychForceOpenMLTSWorkaround
% Force use of OpenML swap completion timestamp workaround. This to
% workaround certain potential timestamping bugs in some graphics drivers.
% Currently no such bugs exist, so this option is just to future-proof
% your Psychtoolbox against potential bugs in future operating systems.
% Bugs were present in Linux versions 3.13 - 3.15 for a short period of
% time between April and July 2014 which made this workaround neccessary.
% However, the workaround is automatically enabled on such Linux versions
% without the need for this conservevram setting. The relevant Linux bugs
% have been fixed by mid-July 2014, ie., at the time of this writing, in all
% versions of the Linux kernel (Linux >= 3.13.11.5+, 3.14.12+, 3.15.5+, 3.16+).
% Nonetheless we leave this manual way to enable the workaround, should the
% need ever arise again in the future.
%
%
% --> It's always better to update your graphics drivers with fixed
% versions or buy proper hardware than using these workarounds. They are
% meant as a last ressort, e.g., if you need to get something going quickly
% or can't get access to bug-fixed drivers.
%

% Need this extra line, separated by a blank line for help to work on Octave - weird.
