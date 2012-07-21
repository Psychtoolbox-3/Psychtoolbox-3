function InitializeMatlabOpenAL(debuglevel, snddevicename, openal_c_style)
% InitializeMatlabOpenAL([debuglevel] [,snddevicename] [, openal_c_style])
%
% InitializeMatlabOpenAL -- Initialize the OpenAL for Matlab wrapper 'moal'.
%
% Call this function at the beginning of your experiment script if you intend
% to use low-level OpenAL sound commands in your script as provided by our
% moalcore extension.
%
% This will check if moal is properly installed and upload all required
% OpenAL constants into your Matlab workspace. It will also set up
% Psychtoolbox for interfacing with external OpenAL code.
%
% MacOS/X users do not need to take any measures, OpenAL comes installed on
% Tiger (10.4) by default. The same is true for modern GNU/Linux distributions.
%
% Users of Microsoft Windows need to download and install the freely available
% OpenAL runtime for Windows. Follow the links to the runtime for Windows in
% the 'Downloads' section of the OpenAL homepage:
%
% http://www.openal.org
%
% Options: All options are optional.
%
% debuglevel = 0 to 3: Setting debuglevel == 0 will disable debug-output.
% A level of 1 will cause MOAL to output error messages and abort your
% scripts execution if it detects any OpenAL error. A level of 2 provides
% additional information that may help you to optimize your code. level 3
% means to be very verbose
%
% snddevicename: Optional. Request a specific sound output device by
% its name. If left out, or if an invalid name is given, OpenAL will
% select a default output device - the most capable and most efficient
% device in your system.
%
% openal_c_style = 0 / 1: Optional. If left out or set to zero, all 
% constants will be loaded in structs in order to avoid cluttering the
% Matlab workspace too much. You'll have to replace all AL_xxx calls by
% AL.xxx calls, e.g., AL_TRUE becomes AL.TRUE .
% If set to one, then all constants will additionally be provided in standard
% C-Style syntax, aka AL_TRUE, but this only works in the main function, not
% in subroutines.
%
% The 'OpenAL for Matlab' low level wrapper moal was implemented by Mario
% Kleiner, trivially derived from our MOGL Matlab for OpenGL wrapper, which
% was developed and contributed to Psychtoolbox under GPL license by
% Prof. Richard F. Murray, University of York, Canada and Mario Kleiner.
%
% MOAL is now licensed under the more permissive MIT license since 2011.
% Relicensing with permission of Richard Murray.

% History:
% 07.02.2007 mk Written - Based on InitializeMatlabOpenGL.
% 27.03.2011 mk Modified: Relicensed to MIT license.

global AL;
global ALC;

% We default to non-C-Style constants if not requested otherwise.
if nargin < 3
    openal_c_style = 0;
end

if isempty(openal_c_style)
   openal_c_style = 0;
end;

if nargin < 2
    snddevicename = [];
end

if nargin < 1
    debuglevel = 1;
end;

if isempty(debuglevel)
    debuglevel = 1;
end;

% Load all GL constants:
evalin('caller','global AL ALC');

% Absolute path to oalconst.mat file: Needed for GNU/Octave:
oalconstpath = [PsychtoolboxRoot 'PsychSound/MOAL/core/oalconst.mat'];

% C-Style compatible load requested?
if openal_c_style > 0
   % Load all constants, also the C-Style ones, e.g., AL_WHATEVER
   evalin('caller',['load ' oalconstpath]);    
else
   % Load only the AL. ALC. versions, e.g., AL.WHATEVER
   % This is less convenient as one needs to replace AL_ by AL. in
   % all scripts, but it doesn't clutter the Matlab workspace...
   evalin('caller',['load ' oalconstpath ' AL ALC']);      
end;

% On Windows, we need to preload moalcore into Matlab while the working
% directory is set to Psychtoolbox/PsychOpenAL/MOAL/core , so the Windows
% dynamic linker can find our own local copy of the OpenAL.dll and link moalcore
% against it.

if IsWin   
   % Windows system: Windows does not have OpenAL.dll installed by default.
   % We try the following: First we perform a PREINIT to load moalcore.dll.
   % If the machine doesn't have a OpenAL DLL installed, the link and load
   % of our DLL will fail. In that case we retry in our own core folder. If
   % that fails again, we tell the user how to get and install the library.
   try
      % Preinit with system installed OpenAL.dll
      moalcore('PREINIT');
   catch %#ok<CTCH>
      % Failed. Retry in our own private OpenAL.dll install location:
      
      % Change working dir to location of our dll's
      olddir = pwd;
      cd([PsychtoolboxRoot 'PsychSound/MOAL/core']);
      % Preload (and thereby link against dlls) moglcore into Matlab. The
      % special command 'PREINIT' forces loading and performs a no-operation.
      try
         moalcore('PREINIT');
      catch %#ok<CTCH>
         targetdir = pwd;
         cd(olddir);
         fprintf('\n');
         fprintf('InitializeMatlabOpenAL: Failed to load our moalcore plugin!\n\n');
         fprintf('The most likely reason is that you did not install the required,\n');
         fprintf('freely available OpenAL DLL library from Creative Labs Inc. on\n');
         fprintf('your system. You can either download the OpenAL runtime from Creative\n');
         fprintf('and install it into the system library folder (recommended), or, if you\n');
         fprintf('do not have write permission to the system folder, you can get the \n');
         fprintf('library from somewhere and copy it into the following Psychtoolbox subfolder:\n\n');
         fprintf('%s \n\n', targetdir);
         fprintf('Type "help InitializeMatlabOpenAL" and read the help text to find out about\n');
         fprintf('download locations for the runtime. After installing the OpenAL32.dll or \n');
         fprintf('64-Bit OpenAL32.dll (Windows 64 bit systems), restart Matlab and retry.\n\n');
         error('Loading Matlab OpenAL wrapper core moalcore failed! OpenAL library not installed?');
      end
     
      % Now that moalcore is (hopefully) properly loaded, we can revert the working
      % directory to its previous setting:
      cd(olddir);
   end   
end;

% Set moalcores debuglevel:
moalcore('DEBUGLEVEL', debuglevel);

% Open and setup the audio device & context:
if isempty(snddevicename)
    % Open default device:
    if debuglevel > 1
        fprintf('InitializeMatlabOpenAL: Trying to initialize default sound device.\n');
    end
    moalcore('almOpen');
else
    % Open specified device:
    if debuglevel > 1
        fprintf('InitializeMatlabOpenAL: Trying to initialize sound device "%s".\n', snddevicename);
    end
    moalcore('almOpen', snddevicename);
end

% Clear error state:
alGetError;

if debuglevel > 1
    fprintf('\nOpenAL device is: %s\n', alcGetString(ALC.DEVICE_SPECIFIER));
    fprintf('OpenAL vendor is: %s\n', alGetString(AL.VENDOR));
    fprintf('OpenAL version is: %s\n', alGetString(AL.VERSION));
    fprintf('OpenAL renderer is: %s\n', alGetString(AL.RENDERER));
    fprintf('OpenAL extensions are: %s\n', alGetString(AL.EXTENSIONS));
    fprintf('OpenAL ALC extensions are: %s\n\n', alcGetString(ALC.EXTENSIONS));
end

% Ready to rock!
return;
