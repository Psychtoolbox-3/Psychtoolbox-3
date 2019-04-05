% RetinaDisplay - How to deal with HiDPI "Retina" displays.
%
% On Linux with the standard X11 display system, no special measures are
% needed and HiDPI "Retina" displays should just work like regular displays
% as far as Psychtoolbox is concerned.
%
% On Apple OSX, Psychtoolbox by default will emulate a standard resolution
% display when used on a Retina display to keep old user code working which
% was written before the era of Retina displays, or with Psychtoolbox
% 3.0.11 or earlier versions. If you want to make full use of the native
% high resolution of your Retina display, use the imaging pipeline task
% PsychImaging('AddTask', 'General', 'UseRetinaResolution');. This will
% expose the full resolution but also require more graphics performance
% from your graphics card.
%
% On MS-Windows the situation can be difficult, depending on the Windows
% version and Matlab version in use, and potentially on what other
% applications are running and their status with respect to HiDPI Retina
% displays:
%
% On Windows XP and earlier versions, stuff will just work, but your
% application windows will look tiny when displayed on a Retina display.
%
% On Windows Vista and Windows 7, Psychtoolbox disables the DWM (Desktop
% Window Manager) by default for fullscreen onscreen windows. This will
% ensure optimal display timing precision and avoid Retina display trouble,
% essentially turning your system into a Windows XP system. DPI aware
% applications will continue to display correctly in such a configuration,
% but non DPI aware apps may appear tiny and difficult to use on Retina
% displays. All versions of Matlab older than R2015b are not DPI aware and
% will likely suffer from this problems.
%
% On Windows 8 and later, Psychtoolbox can not disable the DWM. If you
% display a Psychtoolbox window on a standard resolution display, stuff
% will probably work. If you try to display on a HiDPI Retina display then
% you will need to make sure that you either use GNU/Octave, or you use a
% Matlab version that is DPI aware, ie., at least R2015b. You will probably
% also have to make sure that no other non-DPI aware application is running
% or at least not displaying on your stimulus monitor, as such applications
% would trigger use of the DWM and thereby cause stimulus timing problems.
%
% Additionally, on Windows 8 you must make sure that the stimulation
% display is the designated "primary display", so Windows 8 DPI settings,
% which always apply to the primary display, thereby apply to your stimulus
% display. Applications displayed on other than the primary display will be
% subject to DWM scaling if those other displays do not have the same DPI
% setting as the primary display. A running Matlab or Octave would trigger
% harmful use of the DWM on such a non-matching display. In general it is a
% good idea for timing precision to make your stimulus display the
% designated primary display (this can be set in the display settings
% control display, maybe via a "Make this my main display" checkbox).
%
% If you must use a system setup with Retina display and a version of
% Matlab, or of Windows 8 or later, that does not allow you to circumvent
% or solve these problems, then there exists also the option to either
% override Windows DPI settings (somewhere under display control panel,
% advanced settings or "Make text appear bigger or smaller" - Try a Google
% search to find it) - this would get rid of DPI related timing problems at
% the expense of applications that display too small on your Retina
% display. Another option is to disable DPI scaling for all applications
% (including Matlab!) which are supposed to display on the HiDPI / Retina
% display. This will avoid rescaling and triggering the DWM, at the expense
% of a too small appearance of those applications GUIs on the Retina
% display. You can do this by right-clicking on the applications icon, then
% selecting the "Properties" item from the context menu, then selecting the
% "Compatibility" tab in the property dialog. In that tab, check the field
% called "Disable display scaling on high DPI settings", then click "OK" to
% apply the new settings and restart the application.
%
% See also Psychtoolbox forum message #23414 for another oddity that one user
% needed to do to fix HiDPI on a Windows multi-display setup:
%
% https://groups.yahoo.com/neo/groups/PSYCHTOOLBOX/conversations/messages/23414
%
% Windows is not only troublesome for proper visual stimulation and
% stimulation timing in itself, it can be especially troublesome when used
% with HiDPI / Retina displays. Therefore it might be a good idea to switch
% to Linux, or if you have to use Windows, to stick to Windows-7, or if you
% have to use a later Windows version, upgrade to Windows 10. It is
% probably also a good idea to then limit yourself to a single display
% setup and/or make sure that only DPI aware applications are running on
% this setup.
%
% Good luck!
%
