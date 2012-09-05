function AssertOpenGL
% AssertOpenGL
%
% Break and issue an eror message if the installed Psychtoolbox is not
% based on OpenGL or Screen() is not working properly.
% To date there are four versions of the Psychtoolbox, each based on a
% different graphics library:
%
%  OS9: Psychtoolbox-2, based on Apple's QuickDraw.
%  Win: Psychtoolbox-2, based on Direct X and GDI.
%  Win: OpenGL for the ported OSX-Psychtoolbox, aka Psychtoolbox-3.
%  OSX: OpenGL for Psychtoolbox-3.
%  Linux: OpenGL for Psychtoolbox-3.
%
%  The Psychtoolboxes based on OpenGL are partially incompatible (see below)
%  with previous Psychtoolboxes.  A script which relies on the OpenGL
%  Psychtoolbox should call AssertOpenGL so that it will issue the
%  appropriate warning if a user tries to run it on a computer with a
%  non-OpenGL based Psychtoolbox installed.
%
%  OpenGL-based Psychtoolboxes are distinguised by the availability of these
%  functions:
%
%   Screen('Flip',...);
%   Screen('MakeTexture');
%   Screen('DrawTexture');
%
%
%  If you know you're using Psychtoolbox-3, then the most likely cause for
%  this error message is a problem in the configuration of your
%  Matlab/Octave + System setup, or in the Psychtoolbox installation.
%
%  Typically either the Screen MEX file can't be found or accessed, due to
%  some wrong Matlab/Octave path (Proper Screen file not in path), or due
%  to some permission issue (insufficient security access permissions -
%  typically found on MS-Windows systems), or the Screen MEX file can't be
%  loaded and initialized due to some missing or wrong system library on
%  your machine, e.g., the C runtime library is of an incompatible type.
%  Simply type the command "Screen" at the prompt to see if this may be an
%  issue.
%
%  In both cases, indicated by some "file not found" or "file could not by
%  accesses" or "invalid MEX file" error message, you may want to run the
%  SetupPsychtoolbox command again. This will either fix the problem for
%  you by reconfiguring Psychtoolbox, or it will provide more diagnostic
%  error and troubleshooting messages. Make also sure that you read the
%  troubleshooting tips in the "Download" and "Frequently asked questions"
%  sections of our Psychtoolbox Wiki at http://www.psychtoolbox.org
%
% See also: IsOSX, IsWin, IsLinux.

% HISTORY
% 7/10/04   awi     wrote it.
% 7/13/04   awi     Fixed documentation.
% 10/6/05   awi	    Note here cosmetic changes by dgp between 7/13/04 and 10/6/05
% 12/31/05  mk      Detection code modified to really query type of Screen command (OpenGL?)
%                   instead of OS type, as PTB-OpenGL is now available for Windows as well.
% 06/05/08  mk      Hopefully more diagnostic output message added.
% 05/25/11  mk      Update linker failure message on Octave+Windows: Don't
%                   need ARVideo libraries anymore, but instead the
%                   GStreamer runtime.
% 03/24/12  mk      Add do_braindead_shortcircuit_evaluation() to avoid
%                   warnings on Octave 3.6 and later.
% 09/05/12  mk      Disable shortcircuit operator warning.

% Ok, we sneak this in here, because we don't know a better place for it:
if IsOctave
    % If we're running on a Octave version which supports this...
    if exist('do_braindead_shortcircuit_evaluation', 'builtin')
        % ... we enable Matlab-style short-circuit operator evaluation,
        % ie., treat & and | operators as if they are && and || operators,
        % if they are used in if and while statements and the conditions do
        % evaluate into logical expressions with a simple scalar truth
        % value.
        %
        % This is a hack to silence interpreter warnings and get compatible
        % behaviour. A better fix would be to fix all our M-Files to avoid
        % use of & and |, but this would make ptb completely untestable on
        % Matlab versions prior R2007a for MK. As long as we claim to
        % support pre-R2007a at least to some degree, we should be able to
        % run at least some basic correctness tests:
        do_braindead_shortcircuit_evaluation(1);
        
        % Disable associated Octave warning:
        warning('off', 'Octave:possible-matlab-short-circuit-operator');
    end
end

% We put the detection code into a try-catch-end statement: The old Screen command on windows
% doesn't have a 'Version' subfunction, so it would exit to Matlab with an error.
% We catch this error in the catch-branch and output the "non-OpenGL" error message...
try
   % Query a Screen subfunction that only exists in the new Screen command If this is not
   % OpenGL PTB,we will get thrown into the catch-branch...
   value=Screen('Preference', 'SkipSyncTests'); %#ok<NASGU>
   return;
catch %#ok<*CTCH>
   fprintf('\n\n\nA very simple test call to the Screen() MEX file failed in AssertOpenGL, indicating\n');
   fprintf('that either Screen is totally dysfunctional, or you are trying to run your script on\n');
   fprintf('a system without Psychtoolbox-3 properly installed - or not installed at all.\n\n');

   if IsWin && IsOctave
		le = psychlasterror;
		if ~isempty(strfind(le.message, 'library or dependents')) && ~isempty(strfind(le.message, 'Screen.mex'))
			% Likely the required GStreamer runtimes aren't installed yet!
			fprintf('The most likely cause, based on the fact you are running on Octave under Windows\n');
			fprintf('and given this error message: %s\n', le.message);
			fprintf('is that the required GStreamer runtime libraries are not yet installed on your system.\n\n');
			fprintf('Please type ''help GStreamer'' and read the installation instructions carefully.\n');
			fprintf('After this one-time setup, the Screen command should work properly.\n\n');
			fprintf('If this has been ruled out as a reason for failure, the following could be the case:\n\n');
		end
   end
   
   if IsLinux
       fprintf('\n');
       fprintf('The Psychtoolbox on GNU/Linux needs the following 3rd party libraries\n');
       fprintf('in order to function correctly. If you get "Invalid MEX file errors",\n');
       fprintf('or similar fatal error messages, check if these are installed on your\n');
       fprintf('system and if they are missing, install them via your system specific\n');
       fprintf('software management tools:\n');
       fprintf('\n');
       fprintf('For Screen() and OpenGL support:\n');
       fprintf('* The OpenGL utility toolkit GLUT: glut, glut-3 or freeglut are typical providers.\n');
       fprintf('* GStreamer multimedia framework: At least version 0.10.24 of the core runtime and the gstreamer-base plugins.\n');
       fprintf('  A simple way to get GStreamer at least on Ubuntu Linux is to install the "rhythmbox" or\n');
       fprintf('  "totem" multimedia-players. You may need to install additional packages to play back all\n');
       fprintf('  common audio- and video file formats.\n');
       fprintf('* libusb-1.0 USB low-level access library.\n');
       fprintf('* libdc1394 Firewire video capture library.\n');
       fprintf('* libraw1394 Firewire low-level access library.\n');
       fprintf('\n\n');
       fprintf('For PsychKinect():\n');
       fprintf('* libusb-1.0 USB low-level access library.\n');
       fprintf('* libfreenect: Kinect driver library.\n');
       fprintf('\n');
       fprintf('For Eyelink():\n');
       fprintf('* The Eyelink core libraries from the SR-Research download website.\n');
       fprintf('\n');
       fprintf('\n');

       if ~IsOctave
           s = psychlasterror;
           if ~isempty(strfind(s.message, 'gzopen64'))
               fprintf('YOU SEEM TO HAVE A MATLAB INSTALLATION WITH A BROKEN/OUTDATED libz library!\n');
               fprintf('This is the most likely cause for the error. You can either:\n');
               fprintf('- Upgrade to a more recent version of Matlab in the hope that this fixes the problem.\n');
               fprintf('- Or start Matlab from the commandline with the following command sequence as a workaround:\n\n');
               fprintf('  export LD_PRELOAD=/lib/libz.so.1 ; matlab & \n\n');
               fprintf('  If /lib/libz.so.1 doesn''t exist, try other locations like /usr/lib/libz.so.1 or other names\n');
               fprintf('  like /lib/libz.so, or /usr/lib/libz.so\n');
               fprintf('- A third option is to delete the libz.so library shipped with Matlab. Move away all\n');
               fprintf('  files starting with libz.so from the folder /bin/glnx86 inside the Matlab main folder.\n');
               fprintf('  This way, the linker can''t find Matlabs broken libz anymore and will use the system\n');
               fprintf('  libz and everything will be fine.\n');
               fprintf('\n');
               fprintf('Good luck! Our most heartfelt thanks go to the Mathworks, the unmatched champions in high quality software design.\n');
               fprintf('(Yes, the statement about Mathworks quality of workmanship is meant sarcastic, in case there is any doubt.)\n\n');

               error('Matlab bug -- Outdated/Defective libz installed. Follow above workarounds.');
           end
       end
   end

   % Tried to execute old Screen command of old Win-PTB or MacOS9-PTB. This will tell user about non-OpenGL PTB.
   fprintf('This script or function is designated to run only an Psychtoolbox based on OpenGL. Read "help  AssertOpenGL" for more info.\n\n');
   fprintf('A first more diagnostic test would be to simply type Screen in your Matlab/Octave console and check what its output is.\n\n');
   fprintf('\n\nThe returned error message by Matlab/Octave was:\n');
   ple;
   error('Problems detected in call to AssertOpenGL;');
end;
