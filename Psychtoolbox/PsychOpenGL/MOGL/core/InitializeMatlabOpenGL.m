function varargout = InitializeMatlabOpenGL(opengl_c_style, debuglevel, noswitchto3D, specialFlags)
% InitializeMatlabOpenGL([opengl_c_style] [, debuglevel] [, noswitchto3D] [, specialFlags=0])
%
% InitializeMatlabOpenGL -- Initialize the OpenGL for Matlab wrapper 'mogl'.
%
% Call this function at the beginning of your experiment script before
% calling *any* Psychtoolbox Screen() command, if you intend to use low-level
% OpenGL drawing commands in your script as provided by Richard Murrays
% moglcore extension.
%
% This will check if mogl is properly installed and upload all required
% OpenGL constants into your Matlab workspace. It will also set up
% Psychtoolbox for interfacing with external OpenGL code.
%
% There is also a special query mode: If you set the first argument
% 'opengl_c_style' to the special value -1, then this call will do nothing
% than return the old 'debuglevel', cached from a previous call to this
% routine.
%
%
% Options:
% opengl_c_style = 0 / 1:
% If you call InitializeMatlabOpenGL or InitializeMatlabOpenGL(0), all
% constants will be loaded in structs in order to avoid cluttering the
% Matlab workspace too much. You'll have to replace all GL_xxx calls by
% GL.xxx calls, e.g., GL_LIGHTING becomes GL.LIGHTING .
% If you call InitializeMatlabOpenGL(1), then all constants will additionally
% be provided in standard C-Style syntax, aka GL_LIGHTING.
%
% debuglevel = 0 to 3: Setting debuglevel == 0 will disable debug-output.
% A level of 1 will cause MOGL to output error messages and abort your
% scripts execution if it detects any OpenGL error. A level of 2 provides
% additional information that may help you to optimize your code. level 3
% means to be very verbose.
%
% noswitchto3D: Setting this optional parameter to 1 will only load the GL
% constants and moglcore, but it won't switch PTB itself into 3D mode. This
% is useful if your code needs access to OpenGL for some configuration work
% but doesn't intend to do real 3D rendering. Mostly called from PTB
% internal helper functions, e.g., imaging pipeline. Defaults to zero, aka
% "switch to real 3D mode".
%
% specialFlags = 0: Setting this optional parameter will enable some
% special properties of the created OpenGL context. You can add the
% following values to setup such a special configuration:
%
%    + 2  == Enable and attach an OpenGL accumulation buffer, with
%    requested 16 bits resolution per color component, i.e., R16G16B16A16.
%    The system may decide to allocate an accumulation buffer with more or
%    less than the preferred 16 bpc, or it may decide not to allocate an
%    accumulation buffer at all. On most graphics cards the accumulation
%    buffer is implemented in software - using it may drastically reduce
%    graphics performance down to redraw rates of only a few frames per
%    second! There are better methods based on clever use of alpha blending
%    and floating point resolution offscreen windows on modern graphics
%    cards, or by direct low-level use of framebuffer objects.
%
%    -> This flag is needed to make the glAccum() command work.
%
%
%
% The initial OS/X PowerPC version of the 'OpenGL for Matlab' low level
% interface wrapper mogl was developed, implemented and generously
% contributed to Psychtoolbox under the GPL license by Prof. Richard F.
% Murray, University of York, Canada. Porting to other operating systems
% and architectures, OpenGL 2.x support, and further extensions and
% maintenance has been done by Mario Kleiner.
%
% The code has been relicensed by Richard Murray and Mario Kleiner to the
% more permissive MIT license since 2011.
%

% History:
% 02/16/06 mk Written - Based on glmGetconst from Richard F. Murray.
% 03/27/11 mk Update info about license - New MIT license.

global GL;

% Stor debuglevel internally, so it can be retrieved:
persistent cachedDebuglevel;

if isempty(GL)
    % Is Psychtoolbox properly installed?
    AssertOpenGL;

    % Is Mogl properly installed?
    if ~exist('glmGetConst.m','file'),
        error('Failed to initialize OpenGL for Matlab: Add the ''core'' and ''wrap'' directories of mogl to the MATLAB search path!');
    end;
end

% We default to non-C-Style constants if not requested otherwise.
if nargin < 1
   opengl_c_style = 0;
end;

if isempty(opengl_c_style)
   opengl_c_style = 0;
end;

% Special flag provided?
if opengl_c_style == -1
    % Yep. Calling code wants to query our cached debuglevel:
    varargout{1} = cachedDebuglevel;
    return;
end

if nargin < 2
    debuglevel = 1;
end;

if isempty(debuglevel)
    debuglevel = 1;
end;

if nargin < 3
    noswitchto3D = [];
end

if isempty(noswitchto3D)
    noswitchto3D = 0;
end

if nargin < 4
    specialFlags = [];
end

if isempty(specialFlags)
    specialFlags = 0;
end

% Load all GL constants:
evalin('caller','global AGL GL GLU');

% Absolute path to oglconst.mat file: Needed for GNU/Octave:
oglconstpath = [PsychtoolboxRoot 'PsychOpenGL/MOGL/core/oglconst.mat'];

% C-Style compatible load requested?
if opengl_c_style > 0
   % Load all constants, also the C-Style ones, e.g., GL_LIGHTING
   evalin('caller',['load (''' oglconstpath ''');']);
else
   % Load only the GL. GLU. and AGL. versions, e.g., GL.LIGHTING
   % This is less convenient as one needs to replace GL_ by GL. in
   % all scripts, but it doesn't clutter the Matlab workspace...
   evalin('caller',['load (''' oglconstpath ''', ''AGL'', ''GL'', ''GLU'');']);
end;

% Assign GL_3D manually - Little hack...
GL.GL_3D = 1537;
GL.GL_2D = 1536;

% On Windows, we need to preload moglcore into Matlab while the working
% directory is set to Psychtoolbox/PsychOpenGL/MOGL/core , so the Windows
% dynamic linker can find our own local copy of glut32.dll and link moglcore
% against it.
if IsWin   
   % Windows system: Change working dir to location of our glut32.dll
   olddir = pwd;
   
   if IsWin(1)
       % Need 64-Bit freeglut.dll:
       cd([PsychtoolboxRoot 'PsychOpenGL/MOGL/core/x64']);
   else
       % Need 32-Bit freeglut.dll:
       cd([PsychtoolboxRoot 'PsychOpenGL/MOGL/core']);
   end
   
   % Preload (and thereby link against freeglut.dll) moglcore into Matlab. The
   % special command 'PREINIT' forces loading and performs a no-operation.
   moglcore('PREINIT');
   
   % Now that moglcore is (hopefully) properly loaded, we can revert the working
   % directory to its previous setting:
   cd(olddir);
end;

% Set moglcores debuglevel:
moglcore('DEBUGLEVEL', debuglevel);
% Cache debuglevel internally so it can be queried later on:
cachedDebuglevel = debuglevel;

% Enable support for OpenGL 3D graphics rendering in Psychtoolbox.
if noswitchto3D > 0
    noswitchto3D = 0;
else
    % Merge 1 flag (enable 3D mode) with passed specialFlags:
    noswitchto3D = mor(1, specialFlags);
end

% Conditionally enable full 3D mode of Screen. Enable it if noswitchto3D is
% zero, which is the default:
Screen('Preference', 'Enable3DGraphics', noswitchto3D);

% Ready to rock!
return;
