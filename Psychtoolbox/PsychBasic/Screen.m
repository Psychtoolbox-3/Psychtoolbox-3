function varargout = Screen(varargin)
% Screen is a MEX file for precise control of the video display. Screen has
% many functions; type "Screen" for a list:
% 	Screen
%
% For explanation of any particular screen function, just add a question
% mark "?". E.g. for 'OpenWindow', try either of these equivalent forms:
% 	Screen('OpenWindow?')
% 	Screen OpenWindow?
%
% All the Screen Preference settings are documented together:
% 	Screen Preference?
% 
% General Screen ARGUMENTS, common to most subfunctions of Screen:
% 
% "windowPtr" argument: Screen 'OpenWindow' and 'OpenOffscreenWindow' both
% return a windowPtr, a number that designates the window you just
% created. You can create many windows. To use a window, you pass its
% windowPtr to the Screen function you want to apply to that window.
% 
% "rect" argument: "rect" is a 1x4 matrix containing the coordinates of an
% imaginary box containing all the pixels. All screen and window
% coordinates follow Apple Macintosh conventions. (In Apples the pixels
% occupy the space between the coordinates. Thus a rect [0 0 1 1] contains
% just one pixel.) Coordinates can be local to the window (i.e. 0,0 origin
% is at upper left of window), or local to the screen (origin at upper left
% of screen), or "global", which follows Apple's convention of treating the
% entire desktop (all your screens) as one big screen, with the origin at 
% upper left of the main screen, which has the menu bar. Historically we've
% had two different orderings of the elements of rect, so, for general
% compatibility, all of the Psychophysics Toolbox refers to the elements
% symbolically, through RectLeft, RectTop, etc. Since 2/97, we use Apple's
% standard ordering: RectLeft=1, RectTop=2, RectRight=3, RectBottom=4.
% 
% [optional arguments]: Brackets in the function list, e.g. [color],
% indicate optional arguments, not matrices. Optional arguments must be in
% order, without omitting earlier ones, but you can use the empty matrix
% [] as a place holder, with the same effect as omitting it.
% 
% WHEN YOU GET A MATLAB ERROR
% 
% If your computer only has one screen (the typical scenario) and your
% program produces a Matlab error while your full-screen window is open,
% you'll hear the beep, but you won't be able to see the Matlab Command
% Window. Follow the instructions below for bringing forward the command
% window, then type clear screen to flush just the Screen MEX file, or 
% "clear mex" to flush all the MEX files. When flushed, as part of its 
% exit sequence, Screen closes all its windows, restores the screen's normal 
% color table, and shows the cursor. Or you can get just those effects, 
% without flushing, by calling Screen('CloseAll') or sca - which is an 
% abbreviation for Screen('CloseAll').
%
% You can use Matlab's EVAL command to do this for you automatically. E.g.
% if your program is called "foo.m", run your program by calling EVAL:
% 	eval('foo','clear screen;error(''error in foo'')')
%
% If an error occurs in FOO, Matlab, instead of halting, will execute the
% second argument to EVAL, which restores your screen and reports the
% error.
%
% OpenGL: _________________________________________________________________
%
% Instead of offscreen windows, the OpenGL Psychtoolbox uses fast rendering
% and OpenGL textures for animation. With the exception of matrices, all
% drawing may be done during the animation loop directly to the  onscreen
% window, rather than being rendered to offscreen windows before the start
% of the movie.  Matrices are converted to Textures before the start of the
% animation and, like offscreen windows in OS 9, may be quickly copied to
% an onscreen window during movie play. Offscreen windows are still supported
% if you need to draw very complex stimuli. You can draw the stimulus into
% an offscreen window and then quickly copy the window into the onscreen
% window. For most purposes however, it is possible to draw directly into
% the backbuffer of your offscreen window and make the backbuffer visible
% on next vertical blank by a call to Screen('Flip', windowPtr).  
%
% See MovieDemoOSX and DriftDemoOSX for examples of how to create and show
% movies this way.
%
% Off-screen windows are invisible, but useful as an intermediate place to
% create and store images for later display. Copying from window to window
% is very fast. It's easy to precompute a series of off-screen windows
% and then show them as a movie, in real time, one per video frame:
%
% 		% make movie
% 		window=Screen('OpenWindow', 0, 0);
% 		rect=[0 0 200 200];
%       white = WhiteIndex(window);
% 		for i=1:100
% 			movie(i)=Screen('OpenOffscreenWindow', window, 0, rect);
% 			Screen('FillOval', movie(i), white, [0 0 2 2]*(i-1));
% 		end;
%
% 		% show movie
% 		for i=[1:100 100:-1:1] % forwards and backwards
% 			Screen('CopyWindow',movie(i),window,rect,rect);
% 			Screen('Flip', window);
% 		end;
% 		Screen('CloseAll');
%
%
% Stopping programs:
%
% Command-zero brings the Matlab Command window forward. (Type a zero
% "0" while holding the apple-cloverleaf "command" key down.)
%
% Ctrl-C halts any program.  (Type a "c" while holding down the "Ctrl"
% key). Sometimes, Ctrl-C fails to halt progams executing in a Matlab process
% run with the "-nojvm" option. To halt a runaway Psychtoolbox script in
% Psychtoolbox you might resort to the Windows Task Manager to kill
% the Matlab process.  (Use Ctrl-Alt-Delete to open a window from which
% you can start the Task Manager.)
%
% Windows:
%
% Ctrl-Alt-Delete allows you to launch the Windows task manager, which
% reduces the Psychtoolbox onscreen windows when it opens. (Simultaneosly
% press the "Ctrl", "Alt", and "Delete" keys.)  There are also simpler ways of
% reducing the Psychtoolbox window which are specific to particular
% versions of Windows.
% Windows 2000: 	Alt-Tab will bring another application to the foreground,
% 			minimizing the Matlab Psychtoolbox window.
% 
% OS-X:
% Apple-Command-Escape executes "Force Quit" on Matlab, closing Matlab and all
% of its windows.
%
% Linux:
% Ctrl-Alt-Escape, followed by a mouse click kills the onscreen windows and your
% Matlab session.
%
%
% See "help PsychDemos" for many demos which demonstrate Screen's capabilities.
%
% Differences in Screens capabilities between different operating systems
% are discussed in the online help for the different subfunctions, our
% "PsychDemos" if differences apply, and on the Psychtoolbox Wiki under
% "Platform Differences and writing portable code".
% 
% BUGS
%
% All known bugs and fixes are eventually described at the web site under "Bugs":
% web http://psychtoolbox.org/ ;
%
% Initial reports appear first at the forum:
% web http://www.yahoogroups.com/messages/psychtoolbox/ ;
%
% If you find a bug, please report it to the forum: 
% web mailto:psychtoolbox@yahoogroups.com ;
%
% It will help greatly if you can supply a  minimal-length program that exhibits 
% the bug. And please include as much information about your hardware and software
% setup to document the context in which you're running, e.g., Computer type, graphics
% card type and model, operating system, Matlab version, Psychtoolbox version and flavor
% and the output of PTB to the Matlab window.

% HISTORY
% 7/12/04   awi  ****** OSX fork from the OS9 version *******
%                Divided into sections for OSX, OS9 and Win. 
% 10/4/05   awi  Note here cosmetic changes by dgp on unknown date between 7/12/04 and 10/4/05
% 11/16/06  mk   Rewritten to more closely match reality on PTB-3.
% 01/25/07  mk   Refined again.
AssertMex('Screen.m');
