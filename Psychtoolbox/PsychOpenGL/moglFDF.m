function varargout = moglFDF(cmd, varargin)
% moglFDF(cmd [, arg1][, arg2][, ...]) - "MOGL FormlessDotFields"
%
% Implementation of Sheinberg et al. inspired random dot structure from motion
% rendering. This routine is a fast implementation of "Formless dot field
% structure-from-motion stimuli". It is based on - and very similar in
% behaviour, although not identical in implementation - the algorithm
% proposed by Jedediah M. Singer and David L. Sheinberg in their
% Journal of Vision paper "A method for the real-time rendering of
% formless dot field structure-from-motion stimuli" (Journal of Vision, 8,
% 1-8)
%
% This algorithm takes the idea of the above mentioned paper and pushes it
% one step further by moving nearly all stimulus computation onto the GPU.
%
% All compute intense tasks are carried out by vertex- and fragment-shaders
% on the GPU and all heavy data structures are stored within the GPU's fast
% local VRAM memory, reducing the amount of communication between host
% system and graphics card to an absolute minimum. The Matlab code on the
% CPU only controls the flow of operations on the GPU and generates a
% matrix with random numbers to update the sample distribution. This should
% provide a significant speedup beyond what the Singer et al. algorithm
% achieves, at least for complex and demanding stimuli.
%
%
% The algorithm makes heavy use of GPU based image processing for maximum
% speed, so it needs at least NVidia Geforce 6000 series or ATI Radeon
% X1000 series graphics hardware (and any later models or equivalent
% hardware) to work. It also needs the PTB imaging pipeline enabled, at
% least fast offscreen window support. You do this, e.g., by replacing a
% call to ...
%
% [win, winRect] = Screen('OpenWindow', screenid, 0);
%
% ... with a call sequence like this ...
%
% PsychImaging('PrepareConfiguration');
% PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
% [win , winRect] = PsychImaging('OpenWindow', screenid, 0);
%
%
% A minimal example of how to use moglFDF to render a "dotfield"
% representation of a rotating 3D sphere can be found in the FDFDemo.m
% file.
%
%
% Subcommands, their syntax & meaning:
% ====================================
%
% [oldflag, oldgain] = moglFDF('DebugFlag', flag [, debugGain]);
% - Set debug flag to value 'flag'. Default is zero. Non-zero values enable
% different visualizations that may aid debugging non-working setups.
% 1 = Show silhouette buffer, 2 = Show trackbuffer, 3 = Show random noise
% sampling texture, 4 = Show sampleBuffer, 5 = Show FGDots buffer. A
% setting of -1 shows the real rendered image, instead of the random dot
% visualization. A value of -2 disables any kind of textual warnings.
%
% The optional 'debugGain' parameter must be a 4 component [R G B A] color
% vector with modulation gains for the drawn "debug images" - simply to
% scale each color channel in intensity to allow for display of values
% outside the standard displayable range between zero and one.
%
%
% context = moglFDF('CreateContext', window, rect, texCoordMin, texCoordMax, texResolution, maxFGDots, maxBGDots, dotLifetime [,zThreshold=Off] [,BGSilhouetteAcceptanceProbability=0.0]); 
% - Create a "rendercontext" for a single 3D object. Returns a 'context'
% handle to it which needs to be passed in to all other functions as
% reference. All following parameters are required and don't have any
% defaults:
%
% 'window' Handle of masterwindow - The onscreen window used for rendering.
% This is not neccessarily the window to which final stimulus will be drawn
% to, but it is needed as a "parent" for all ressources.
%
% 'rect' A Psychtoolbox rectangle [left top right bottom] that describes
% the size and shape of the final stimulus window. This rect must have the
% same size as the 3D window and final stimulus window -- Lots of internal
% calculations depend on this geometry spec.
%
% 'texCoordMin' Two element vector which contains the minimum texture
% coordinate values contained in the 3D scene for x- resp. y-direction.
%
% 'texCoordMax' Two element vector which contains the maximum texture
% coordinate values contained in the 3D scene for x- resp. y-direction.
%
% 'texResolution' Two element vector which contains the internal resolution
% for x- resp. y-direction of the 3D object surface. Higher values mean finer
% resolution and less aliasing, but also higher storage requirements and
% longer processing times.
%
% 'maxFGDots' Maximum number of foreground (object shape) dots to use for
% random shape sampling. This must be an integral multiple of
% 'dotLifetime'. If it isn't, it will get adjusted to become an integral
% multiple.
%
% 'maxBGDots' Maximum number of background dots to use for random background
% sampling. This must be an integral multiple of 'dotLifetime'. If it
% isn't, it will get adjusted to become an integral multiple. If you don't
% want to have structure cues in your stimulus, you should set 'maxBGDots'
% equal to 'maxFGDots' to keep overall dot density on the display constant.
%
% 'dotLifetime' Lifetime of each foreground- or background dot in 'Update'
% cycles. Each dot is replace by a new random sample after that many
% invocations of the 'Update' function.
%
% 'zThreshold' Optional zThreshold for occlusion test: By default, it is
% 10.0 ie. occlusion test disabled. A value between 0.0 and 1.0 will enable
% occlusion testing -- Dots that would correspond to occluded surfaces are
% not drawn. Small numbers (close to zero) make the test more sensitive but
% can cause artifacts due to internal numeric roundoff errors. Bigger
% numbers (closer to one) make it more robust but less powerful. The
% "sweet-spot" depends on your hardware and 3D scene. Empirically a setting
% of 0.0001 is a good value for ATI Radeon X1000 series hardware.
% The default setting (bigger than 1.0) will disable occlusion test --
% "Hidden dots" are not hidden, but drawn.
%
% 'BGSilhouetteAcceptanceProbability' Optional BGSilhouetteAcceptanceProbability
% This is the probability with which a dot from the "background distribution" 
% will be drawn if it is actually located in the area of the objects
% silhouette. A value of 0.0 (which is the default) will not draw any
% background dots within the objects silhouette. Values between 0 and 1
% correspond to acceptance probabilities between 0% and 100%. If you want
% to keep the overall dot density of foreground dots and background dots
% constant (in order to not provide segmentation cues based on structure),
% you should set the 'maxFGDots' parameter like this:
%
% maxFGDots = (1 - BGSilhouetteAcceptanceProbability) * maxBGDots;
%
%
% context = moglFDF('SetRenderCallback', context, callbackEvalString);
% - Define the 'eval' string for this context to be used as rendercallback.
% Pass in a Matlab command string (for evaluation via eval() function in the
% Workspace of the calling function). This string is called/executed during
% each 'Update' cycle. It has to contain the code that performs the actual
% rendering of the 3D scene or object.
%
% The called rendering code *must not* glClear() the framebuffer or mess
% around with alpha-blending state or depth-buffer/depth-test settings, nor
% should it bind any shaders! It makes sense to disable any kind of
% lighting or texture mapping, as no photorealistic image is rendered, so
% it would be a waste of computation time.
%
%
% context = moglFDF('ReinitContext', context, rect, texCoordMin, texCoordMax, texResolution, maxFGDots, maxBGDots, dotLifetime [,zThreshold=Off] [,BGSilhouetteAcceptanceProbability=0.0]); 
% - Reinitialize an already existing context with new stimulus parameters.
% The parameters are identical to the ones in 'CreateContext', except for
% the first one: You don't pass a windowhandle of a parent window, as this
% stays the same for the reinitialized context. Instead you pass the handle
% of the 'context' to reinitialize.
%
% 'ReinitContext' is the same as a sequence of 'DestroyContext', followed
% by a new 'CreateContext', except that it is optimized for speed --
% Reinitialization with new parameters is typically at least 3 times faster
% than a full destroy & recreate operation.
%
%
% context = moglFDF('DestroyContext', context);
% - Destroy a processing context, release all of its ressources.
%
%
% context = moglFDF('ResetState', context);
% - Reset processing contexts state to initial state, just as if it was
% just created. Useful at start of a new trial. Another way to start a new
% trial, but with a full distribution already initialized, is to use the
% moglFDF('Update') call with the 'instantOn' flag set to 1 for the first
% iteration of your stimulus loop, instead of the default of zero.
%
%
% context = moglFDF('SetColorTexture', context, textureId, textureTarget);
% - Assign a regular color texture map with handle 'textureId' and texture
% mapping target 'textureTarget' to 'context'. This will enable assignment
% of colors to drawn 2D dots (in moglFDF('Render',...);) and fetch the
% relevant per-dot colors from the assigned texture map 'textureId'.
%
% Assigning an empty or negative textureId will disable texture mapping.
% Texture mapping is disabled by default, i.e. at context creation time.
%
%
% context = moglFDF('SetDrawShader', context, fgShaderId [, bgShaderId] [, needSprites]);
% - Assign a GLSL shader with handle 'fgShaderId' during 2D drawing of
% foreground dots in moglFDF('Render',...); Passing a 'fgShaderId' which is
% empty or negative disables shading. Shading is disabled by default.
%
% The optional 'bgShaderId' assigns potential shaders for drawing of
% background dots.
%
% The optional flag 'needSprites' if set to 1, will enable generation of
% point-sprite texture coordinates on texture unit 1 while using a shader
% with point-smoothing enabled. A setting of 0 disables point sprites.
% Point sprites plus special code within your drawing fragment shader are
% needed if you want to draw nicely anti-aliased dots on GPUs that don't
% support simultaneous use of fragment shaders and anti-aliased dots. On
% such systems you can roll your own anti-aliasing via point-sprites.
% Please note that almost all consumer class GPU's don't support
% anti-aliased dots in conjunction with fragment shaders.
%
%
% context = moglFDF('Update', context [, instantOn=0]);
% - Perform an 'update' cycle for given context. A new "3D frame" is rendered
% via the rendercallback function, then analysed, resampled etc. to create
% a new complete distribution of 2D random dots, ready for drawing or
% readback. If the optional 'instantOn' flag is provided and non-zero, then
% the whole distribution is generated at once for a quick start at the
% beginning of a new trial, otherwise only one batch of samples is updated.
% By default, only one batch is updated, as required for the algorithm to
% work.
%
%
% context = moglFDF('Render', context [, targetWindow] [, drawSpec=[1,1]]);
% - Render current 2D random dot cloud (as defined by processing of last
% 'Update' call) quickly and efficiently into window 'targetWindow'.
% 'targetWindow' can be any onscreen- or offscreen window and is allowed to
% change at each invocation of 'Render'. By default, the 'window' from the
% 'CreateContext' call is used.
%
% 'drawSpec' is an optional two-element vector to select if only foreground
% dots should be rendered [1 0], only background dots should be rendered [0 1],
% or both [1 1] -- which is the default.
%
% Before calling 'Render' you can define dot sizes, colors, alpha-blending
% state, texture coordinates, anti-aliasing settings, or define texture
% mapping, point-sprite modes or texture mapping setups however you like.
% The internal 'Render' routine just defines 2D point locations, then
% invokes the render op.
%
%
% [xyFGdots, xyBGdots, uvFGdots] = moglFDF('GetResults', context); - Returns a 2 row
% by n columns vector of all random dot positions, for processing within
% Matlab/Octave. Row 1 is x-locations, Row 2 is y-locations of dots, each
% column defines one dot. The 'xyFGDots' contains all foreground dots which
% define the object, whereas the 'xyBGdots' vector contains the background
% dots. These vectors are suitable for direct drawing via
% Screen('DrawDots'); However, invocation of moglFDF('Render',...); is a
% more efficient method of rendering these dot fields, unless you have very
% special needs.
%
% The optional 'uvFGdots' argument returns 2D texture coordinates as
% assigned to the rendered 3D object.
%

% History:
%  05/02/08  Initial "proof of concept" implementation (MK).
% -05/10/08  Various improvements (MK).
%  11/03/08  Documentation update, preparation for public release (MK).
%  11/15/08  Improvements: New 'instantOn' mode for 'Update' method, new
%            'ReinitContext' function (MK).

% Need OpenGL constants:
global GL;

% Internal state:
global moglFDF_OriginalContext;
persistent contextcount;
persistent debug;
persistent debugGain;

if nargin < 1
    error('You must provide a "cmd" subcommand to execute!');
end

if isempty(contextcount)
    contextcount = 0;
    moglFDF_OriginalContext = [];
    debug = 0;
    
    if isempty(GL)
        % If OpenGL not initialized, do a full init for 3D mode:
        error('OpenGL mode not initialized! You *must* call InitializeMatlabOpenGL before the first call to this routine or any Screen() function!')
    end
end

% Subcommand dispatch:

% Initialization of a new context: Allocate and setup all ressources:
if strcmpi(cmd, 'CreateContext') || strcmpi(cmd, 'ReinitContext')
    % Fetch all arguments - They are all required.
    if nargin < 9
        error(sprintf('Some mandatory input arguments to "%s" are missing. Please provide them!', cmd)); %#ok<SPERR>
    end
    
    if strcmpi(cmd, 'CreateContext')
        createContext = 1;
    else
        createContext = 0;
    end
    
    % First time init? I.e. is this the first context to be created?
    if contextcount == 0
        % Yes. Perform all one-time initialization work and create a
        % template context from which all other contexts can be derived:
        
    end

    % Type of expected first argument depends if this is a 'CreateContext'
    % call or a 'ReinitContext' call:
    if createContext
        % Parent window: Provides OpenGL master-/slave- contexts for our
        % operations, shaders and buffers, as well as reference for
        % rendertarget size:
        ctx.parentWin = varargin{1};

        if ~isscalar(ctx.parentWin) || ~ismember(ctx.parentWin, Screen('Windows'))
            disp(ctx.parentWin);
            error('Invalid "window" argument provided to "CreateContext" - No such window (see above)!');
        end
    else
        % Old 'ctx' handle of context to reparameterize / reinit:
        ctx = varargin{1};
        
        % Make sure we've got a valid handle:
        if ~isstruct(ctx)
            disp(ctx);
            error('Invalid "context" argument provided to "ReinitContext" - This is not a valid moglFDF context handle!');
        end
        
        if ~isfield(ctx, 'moglFDFMagic')
            disp(ctx);
            error('Invalid "context" argument provided to "ReinitContext" - This is not a valid moglFDF context handle!');
        end

        % Release all buffers, but not the shaders!
        % Delete all offscreen windows, gloperators and buffers like IBO's
        % VBO's, PBO's, FBO's etc, so they can get recreated, based on the
        % new context parameters:
        deleteContextBuffers(ctx);
        RestoreGL;
        
        % Decrement contextcount, so it can be reincremented at end of this
        % function:
        contextcount = max(contextcount - 1, 0);
    end
    
    % Assign our magic cookie...
    ctx.moglFDFMagic = 'Funky magic-cookie';
    
    % Get all other arguments and perform parameter type and range checks:
    ctx.rect = varargin{2};
    if ~isnumeric(ctx.rect) || length(ctx.rect)~=4
        disp(ctx.rect);
        error('Invalid "rect" argument provided to "CreateContext" - Must be a 4 component vector that describes the size and shape of the target rectangle [left top right bottom]');        
    end
    
    ctx.rect = double(ctx.rect);
    if IsEmptyRect(ctx.rect)
        disp(ctx.rect);
        error('Invalid "rect" argument provided to "CreateContext" - Must be a non-empty rect that describes the size and shape of the target rectangle [left top right bottom]');        
    end
    
    ctx.texCoordMin = varargin{3};
    if ~isnumeric(ctx.texCoordMin) || length(ctx.texCoordMin)~=2
        disp(ctx.texCoordMin);
        error('Invalid "texCoordMin" argument provided to "CreateContext" - Must be a 2 component vector of minimal texture coordinates in x- and y- direction!');
    end
    
    ctx.texCoordMax = varargin{4};
    if ~isnumeric(ctx.texCoordMax) || length(ctx.texCoordMax)~=2
        disp(ctx.texCoordMax);
        error('Invalid "texCoordMax" argument provided to "CreateContext" - Must be a 2 component vector of maximal texture coordinates in x- and y- direction!');
    end
    
    ctx.texResolution = varargin{5};
    if ~isnumeric(ctx.texResolution) || length(ctx.texResolution)~=2
        disp(ctx.texResolution);
        error('Invalid "texResolution" argument provided to "CreateContext" - Must be a 2 component vector of integral numbers with processing resolution in x- and y- direction!');
    end
   
    if round(ctx.texResolution) ~= ctx.texResolution
        disp(ctx.texResolution);
        error('Invalid "texResolution" argument provided to "CreateContext" - Must be integral numbers for resolution in x- and y- direction!');
    end

    if min(ctx.texResolution) < 2
        disp(ctx.texResolution);
        error('Invalid "texResolution" argument provided to "CreateContext" - Each component must be at least 2 units!');
    end
    
    ctx.maxFGDots = round(varargin{6});
    if ~isscalar(ctx.maxFGDots)
        disp(ctx.maxFGDots);
        error('Invalid "maxFGDots" argument provided to "CreateContext" - Must be a positive integral number of maximum foreground dots!');
    end
    
    if ctx.maxFGDots < 1
        disp(ctx.maxFGDots);
        fprintf('Invalid "maxFGDots" argument provided to "CreateContext" - Must be at least 1! Changed to 1.\n');
        ctx.maxFGDots = 1;
    end
    
    ctx.maxBGDots = round(varargin{7});
    if ~isscalar(ctx.maxBGDots)
        disp(ctx.maxBGDots);
        error('Invalid "maxBGDots" argument provided to "CreateContext" - Must be a positive integral number of maximum background dots!');
    end
    
    if ctx.maxBGDots < 1
        disp(ctx.maxBGDots);
        fprintf('Invalid "maxBGDots" argument provided to "CreateContext" - Must be at least 1! Changed to 1.\n');
        ctx.maxBGDots = 1;
    end

    ctx.dotLifetime = round(varargin{8});
    if ~isscalar(ctx.dotLifetime)
        disp(ctx.dotLifetime);
        error('Invalid "dotLifetime" argument provided to "CreateContext" - Must be a positive integral number of how many update cycles a dot lives before replacement!');
    end
    
    if ctx.dotLifetime < 1
        disp(ctx.dotLifetime);
        error('Invalid "dotLifetime" argument provided to "CreateContext" - Must be at least 1!');
    end

    % Basic checks passed: Now check for inter-parameter consistency:
    if rem(ctx.maxFGDots, ctx.dotLifetime)
        % Doesn't divide without remainder.
        if debug ~= -2
            fprintf('In moglFDF:%s:\n', cmd);
            fprintf('maxFGDots=%i , dotLifetime=%i --> remainder of maxFGDots / dotLifetime is not zero, as required!\n', ctx.maxFGDots, ctx.dotLifetime);
        end

        % Modify it to satisfy condition:
        ctx.maxFGDots = max(ceil(ctx.maxFGDots / ctx.dotLifetime) * ctx.dotLifetime, ctx.dotLifetime);

        if debug ~= -2
            fprintf('"maxFGDots" must be an integral multiple of "dotLifetime". Changed "maxFGDots" to a value of %i to satisfy this condition.\n', ctx.maxFGDots);
        end
    end

    if rem(ctx.maxBGDots, ctx.dotLifetime)
        % Doesn't divide without remainder.
        if debug ~= -2
            fprintf('In moglFDF:%s:\n', cmd);
            fprintf('maxBGDots=%i , dotLifetime=%i --> remainder of maxBGDots / dotLifetime is not zero, as required!\n', ctx.maxBGDots, ctx.dotLifetime);
        end

        % Modify it to satisfy condition:
        ctx.maxBGDots = max(ceil(ctx.maxBGDots / ctx.dotLifetime) * ctx.dotLifetime, ctx.dotLifetime);

        if debug ~= -2
            fprintf('"maxBGDots" must be an integral multiple of "dotLifetime". Changed "maxBGDots" to a value of %i to satisfy this condition.\n', ctx.maxBGDots);
        end
    end
    
    if nargin >= 10
        ctx.zThreshold = varargin{9};
        if ~isscalar(ctx.zThreshold)
            disp(ctx.zThreshold);
            error('Invalid "zThreshold" argument provided to "CreateContext" - Must be a positive number in range 0.0 - 1.0 for z-Test, and bigger for z-Test disabled!');
        end
    else
        % Default: z-Test disabled:
        ctx.zThreshold = 10.0;
    end
    
    if nargin >= 11
        ctx.BGSilhouetteAcceptanceProbability = varargin{10};
        if ~isscalar(ctx.BGSilhouetteAcceptanceProbability)
            disp(ctx.BGSilhouetteAcceptanceProbability);
            error('Invalid "BGSilhouetteAcceptanceProbability" argument provided to "CreateContext" - Must be a positive number in range 0.0 - 1.0!');
        end        
    else
        ctx.BGSilhouetteAcceptanceProbability = 0.0;
    end
    
    % Backup current GL context binding:
    BackupGL;
    
    % Make sure our Screen context is active:
    SwitchToPTB;
    
    % Retrieve info about our hosting window. This will implicitely enable
    % our parents OpenGL context, so we can do GL query commands safely:
    winfo = Screen('GetWindowInfo', ctx.parentWin);
    
    % Retrieve maximum width or height of textures and offscreen windows
    % supported by this GL implementation:
    maxtexsize = glGetIntegerv(GL.MAX_RECTANGLE_TEXTURE_SIZE_EXT);
    
    % Width of a line in samplebuffer:
    ctx.samplesPerLine = min(maxtexsize, 2048);
    ctx.samplesPerLine = min(ctx.samplesPerLine, round(ctx.maxFGDots / ctx.dotLifetime));
    
    % Check requested internal resolution against hw-limit:
    if max(ctx.texResolution) > maxtexsize
        disp(ctx.texResolution)
        error(sprintf('Requested "texResolution" parameter too big in at least one dimension - Your graphics card can not handle that! Maximum is %i\n', maxtexsize)); %#ok<SPERR>
    end

    % Need these checks only on original context creation:
    if createContext
        % Imaging pipeline active in at least minimum configuration?
        if ~bitand(winfo.ImagingMode, mor(kPsychNeedFastBackingStore, kPsychNeedFastOffscreenWindows))
            % Neither basic pipeline, nor fast offscreen window support
            % activated in parent window. This is a no-go!
            error('In "CreateContext": The PTB imaging pipeline is not active for provided parent window - this will not work! Need at least support for fast offscreen windows.');
        end

        if winfo.GLSupportsFBOUpToBpc < 32 || winfo.GLSupportsTexturesUpToBpc < 32
            error('In "CreateContext": Your gfx-hardware is not capable of handling textures and buffers with the required precision - this function will not work on your hardware!');
        end

        % Check for all required extensions:
        if ~(~isempty(findstr(glGetString(GL.EXTENSIONS), '_vertex_buffer_object')) && ...
                ~isempty(findstr(glGetString(GL.EXTENSIONS), '_pixel_buffer_object')) && ...
                ~isempty(findstr(glGetString(GL.EXTENSIONS), '_framebuffer_object')) && ...
                ~isempty(findstr(glGetString(GL.EXTENSIONS), 'GL_ARB_shading_language')) && ...
                ~isempty(findstr(glGetString(GL.EXTENSIONS), 'GL_ARB_shader_objects')) && ...
                ~isempty(findstr(glGetString(GL.EXTENSIONS), 'GL_ARB_fragment_shader')) && ...
                ~isempty(findstr(glGetString(GL.EXTENSIONS), 'GL_ARB_vertex_shader')) && ...
                (~isempty(findstr(glGetString(GL.EXTENSIONS), 'GL_APPLE_float_pixels')) || ...
                ~isempty(findstr(glGetString(GL.EXTENSIONS), '_color_buffer_float'))))
            % At least one of the required extensions is missing!
            error('In "CreateContext": Your gfx-hardware does not support all required OpenGL extensions - this function will not work on your hardware!');
        end
    end
    
    % Ok, all checks passed.
    
    % Create all relevant FBO buffers, aka Offscreen windows:

    % Sample buffer: Contains all foreground random dot samples.
    % The buffer is split into ctx.dotLifetime blocks, each
    % ctx.samplesPerLine wide and ctx.sampleLinesperBatch high. Each
    % pixel in the buffer encodes one foreground sampe - and thereby one
    % potential foreground dot in the next rendered frame.
    %
    % The algorithm will at each update cycle replace one of these blocks
    % or "batches" with a new distribution of samples -- the old dots in
    % that batch will die and get replaced by new samples from the
    % foreground distribution iff they actually hit the silhouette of the
    % 3D object to be "visualized" - otherwise that samples will be invalid
    % and inactive. Each update cycle a different batch is selected for
    % update in a round-robin fashion, so each sample has a lifetime of
    % ctx.dotLifetime update cycles as requested.
    %
    % Each pixel codes as follows: Red channel == x position of 3D object
    % surface to which the sample is attached, encoded as 32 bit float --
    % technically it is the s-texture coordinate of the 3D surface at the
    % point of impact of our sample. Green channel == y position, aka
    % t-texture coordinate on objects surface. Blue channel == "Valid"
    % flag: Non-zero means: Visualize this samples final tracked position
    % -- draw corresponding dot. Zero means: Ignore this sample.
    %
    % Our buffer is implemented as a FBO backed floating point offscreen
    % window with a pixel size of 128 bits, aka 32 bpc float.
    ctx.sampleLinesPerBatch = ceil((round(ctx.maxFGDots / ctx.dotLifetime)) / ctx.samplesPerLine);
    ctx.sampleLinesTotal = ctx.sampleLinesPerBatch * ctx.dotLifetime;
    ctx.sampleBuffer = Screen('OpenOffscreenWindow', ctx.parentWin, [0 0 0 0], double([0 0 ctx.samplesPerLine ctx.sampleLinesTotal]), 128, 32);
    ctx.maxFGDots = ctx.sampleLinesTotal * ctx.samplesPerLine;
    
    % Silhouette buffer: Contains the "perspective correct image space"
    % image of the rendered 3D object. The object is rendered "normally" in
    % a first render pass to get its silhouette, except that the color of
    % each rendered pixel is not a shaded/lit color, but its encoded
    % interpolated surface texture coordinate. This will be done by a
    % proper fragment shader during render pass.
    % The buffer is an FBO backed offscreen window the same size as the
    % parent window, but with 32bpc float format to store accurate texture
    % coordinates in the pixel colors. Red channel encodes s-coord, Green
    % channel encodes t-coord of 2D texture coordinate, blue encodes
    % interpolated z-buffer depths.
    [ctx.silhouetteWidth, ctx.silhouetteHeight] = RectSize(ctx.rect);
    ctx.silhouetteBuffer = Screen('OpenOffscreenWindow', ctx.parentWin, [0 0 0 0], [0, 0, ctx.silhouetteWidth, ctx.silhouetteHeight], 128, 32);

    % Retrieve OpenGL texture handle for the sihouetteBuffer:
    ctx.silhouetteTexture = Screen('GetOpenGLTexture', ctx.parentWin, ctx.silhouetteBuffer);
        
    % Tracking buffer: Contains the unwarped/flattened image of the 3D
    % objects surface, created by a 2nd rendering pass of the 3D object,
    % but with special vertex-/fragment shaders attached.
    %
    % Pixel location (x,y) encodes for object surface texture coordinate
    % (s,t): Its R, G and B channels encode interpolated X, Y, Z coordinate
    % of the unwarped object surface. This way a lookup at position (s,t)
    % provides the image space 3D coordinate of surface point (s,t) in
    % "normal" camera centered and projected (X,Y,Z) space -- (X,Y) are
    % projected image coordinates, (Z) is depths component. This allows to
    % lookup the tracked/updated (X,Y) image position of a sample dot.
    %
    % Again a 32bpc float offscreen window FBO, but the resolution is
    % chosen per user spec to be fine enough in texture coordinate space to
    % avoid aliasing artifacts as good as possible:
    ctx.trackingBuffer = Screen('OpenOffscreenWindow', ctx.parentWin, [0 0 0 0], [0, 0, ctx.texResolution(1), ctx.texResolution(2)], 128, 32);

    % Final buffer with foreground dot positions. This one will get filled
    % by the createFGDotShader. It will later get either read back to
    % Matlab on usercode request, or converted to a VBO and then rendered.
    ctx.FGDotsBuffer = Screen('OpenOffscreenWindow', ctx.parentWin, [0 0 0 0], Screen('Rect', ctx.sampleBuffer), 128, 32);
    

    % Final buffer with background dot positions. This one will get filled
    % by the createBGDotShader. It will later get either read back to
    % Matlab on usercode request, or converted to a VBO and then rendered.
    
    % Width of a line in buffer:
    ctx.BGsamplesPerLine = min(maxtexsize, 2048);
    ctx.BGsamplesPerLine = min(ctx.BGsamplesPerLine, round(ctx.maxBGDots / ctx.dotLifetime));
    ctx.BGsampleLinesPerBatch = ceil((round(ctx.maxBGDots / ctx.dotLifetime)) / ctx.BGsamplesPerLine);
    ctx.BGsampleLinesTotal = ctx.BGsampleLinesPerBatch * ctx.dotLifetime;
    ctx.BGDotsBuffer = Screen('OpenOffscreenWindow', ctx.parentWin, [0 0 0 0], double([0 0 ctx.BGsamplesPerLine ctx.BGsampleLinesTotal]), 128, 32);
    ctx.maxBGDots = ctx.BGsampleLinesTotal * ctx.BGsamplesPerLine;
    ctx.BGSampleSet = zeros(ctx.BGsampleLinesTotal, ctx.BGsamplesPerLine, 3);

    % Load all our shaders - Need to do this only on original context
    % creation, as shaders are recycled across context reinits. However, we
    % can only recycle shaders from one existing context, not across
    % different contexts, because each shader object also encapsulates
    % per-context state like the settings of all Uniforms etc. and we can't
    % share these!
    if createContext
        % Basepath to shaders:
        % shaderpath = [fileparts(mfilename('fullpath')) filesep ];
        shaderpath = '';

        % Shader for 1st object renderpass: Encode texcoords and depths into
        % color channel -- to fill silhouetteBuffer:
        ctx.silhouetteRenderShader = LoadGLSLProgramFromFiles([shaderpath 'moglFDFSilhouetteRenderShader'], 1);

        % Shader for 2nd object renderpass: Fill trackingBuffer
        ctx.trackingRenderShader = LoadGLSLProgramFromFiles([shaderpath 'moglFDFTrackingRenderShader'], 1);

        % Shader for update of distribution in sampleBuffer:
        ctx.samplingShader = LoadGLSLProgramFromFiles([shaderpath 'moglFDFSamplingShader'], 1);

        % Shader for final creation of foreground dots VBO spec from
        % distribution in sampleBuffer and trackingBuffer:
        ctx.createFGDotsShader = LoadGLSLProgramFromFiles([shaderpath 'moglFDFForegroundDotsRenderShader'], 1);
    
        % Shader for creation of background dots VBO spec:
        ctx.createBGDotsShader = LoadGLSLProgramFromFiles([shaderpath 'moglFDFBackgroundDotsRenderShader'], 1);
        
        % Setup default coloring mode: No textures, no texture mapping:
        ctx.colorTexId = -1;
        ctx.colorTexTarget = -1;

        % Also disable 2D dot drawing shaders by default:
        ctx.draw2DShaderFG = -1;
        ctx.draw2DShaderBG = -1;
        ctx.needSprites = 0;
    end
    
    % Setup trackingRenderShader:
    glUseProgram(ctx.trackingRenderShader)

    % Compute texture coordinate offset and multiplier to apply in order to
    % remap the real texture coordinate range into the normalized [-1:+1]
    % interval that can pass through vertex clipping:
    glUniform4f(glGetUniformLocation(ctx.trackingRenderShader, 'TextureOffsetBias'), ctx.texCoordMin(1), ctx.texCoordMin(2), 2.0/(ctx.texCoordMax(1) - ctx.texCoordMin(1)), 2.0/(ctx.texCoordMax(2) - ctx.texCoordMin(2)));

    % Bind texunit 1 to object coordinates texture:
    glUniform4f(glGetUniformLocation(ctx.trackingRenderShader, 'Viewport'), 0, 0, ctx.silhouetteWidth/2, ctx.silhouetteHeight/2);
    
    % Setup shader for update of distribution in sampleBuffer:
    glUseProgram(ctx.samplingShader);

    % Bind texunit 0 to random sample position texture:
    glUniform1i(glGetUniformLocation(ctx.samplingShader, 'SilSamplePositions'), 0);

    % Bind texunit 1 to silhouette image texture:
    glUniform1i(glGetUniformLocation(ctx.samplingShader, 'Silhouette'), 1);

    % Define remapping of texture coordinates into range 0-texResolution --
    % The size of the trackingBuffer:
    glUniform4f(glGetUniformLocation(ctx.samplingShader, 'TextureOffsetBias'), ctx.texCoordMin(1), ctx.texCoordMin(2), ctx.texResolution(1)/(ctx.texCoordMax(1) - ctx.texCoordMin(1)), ctx.texResolution(2)/(ctx.texCoordMax(2) - ctx.texCoordMin(2)));
    
    % Setup shader for final creation of foreground dots VBO spec from
    % distribution in sampleBuffer and trackingBuffer:
    glUseProgram(ctx.createFGDotsShader)

    % Bind texunit 0 to Samplebuffer texture:
    glUniform1i(glGetUniformLocation(ctx.createFGDotsShader, 'SampleBuffer'), 0);

    % Bind texunit 1 to object coordinates texture:
    glUniform1i(glGetUniformLocation(ctx.createFGDotsShader, 'GeometryBuffer'), 1);

    % Bind texunit 2 to silhouette texture for last rendered frame:
    glUniform1i(glGetUniformLocation(ctx.createFGDotsShader, 'Silhouette'), 2);
    
    % Assign zThreshold for depths testing of foreground dots before
    % output to handle occlusions correctly:
    glUniform1f(glGetUniformLocation(ctx.createFGDotsShader, 'zThreshold'), ctx.zThreshold);

    % Assign height of final output window + 1 to allow shader to invert
    % y-coordinate of final dots properly to account for difference in
    % y-axis direction of Screen()'s reference frame vs. OpenGL default
    % frame:
    glUniform1f(glGetUniformLocation(ctx.createFGDotsShader, 'ViewportHeight'), ctx.silhouetteHeight + 1);
    
    % Define size of GeometryBuffer -- wrapAround values for interpolated
    % texture lookup coordinates:
    glUniform2f(glGetUniformLocation(ctx.createFGDotsShader, 'texWrapAround'), ctx.texResolution(1), ctx.texResolution(2));

    % Define inverse remapping of texture coordinates into range
    % 0-texResolution -- The size of the trackingBuffer. N.B.: A neutral
    % mapping would be (0, 0, 1, 1) - That would pass trackingBuffer
    % texture coordinates instead of object texture coordinates.
    glUniform4f(glGetUniformLocation(ctx.createFGDotsShader, 'TextureOffsetBias'), ctx.texCoordMin(1), ctx.texCoordMin(2), 1 / (ctx.texResolution(1)/(ctx.texCoordMax(1) - ctx.texCoordMin(1))), 1 / (ctx.texResolution(2)/(ctx.texCoordMax(2) - ctx.texCoordMin(2))));
    
    % Set default 'clipVertex' position to (x,y,u,v) = (-1, 0, 0, 0): This
    % will prevent any vertex to which this is applied from drawing,
    % because it is clipped away due to its negative x-location outside viewport:
    glUniform4f(glGetUniformLocation(ctx.createFGDotsShader, 'clipVertex'), -1, 0, 0, 0);

    glUseProgram(0);
    
    % Create gloperator from shader for later use by Screen('TransformTexture'):
    ctx.createFGDotsoperator = CreateGLOperator(ctx.parentWin, [], ctx.createFGDotsShader, 'Create foreground dots.');
    
    % Setup shader for creation of background dots VBO spec:
    glUseProgram(ctx.createBGDotsShader)

    % Bind texunit 0 to random sample position texture:
    glUniform1i(glGetUniformLocation(ctx.createBGDotsShader, 'SilSamplePositions'), 0);

    % Bind texunit 1 to silhouette image texture:
    glUniform1i(glGetUniformLocation(ctx.createBGDotsShader, 'Silhouette'), 1);

    % Assign height of final output window + 1 to allow shader to invert
    % y-coordinate of final dots properly to account for difference in
    % y-axis direction of Screen()'s reference frame vs. OpenGL default
    % frame:
    glUniform1f(glGetUniformLocation(ctx.createBGDotsShader, 'ViewportHeight'), ctx.silhouetteHeight + 1);
    
    % Assign [0;1] acceptance threshold value for accepting a background
    % distribution dot which lies within the objects silhouette:
    glUniform1f(glGetUniformLocation(ctx.createBGDotsShader, 'SilAcceptThreshold'), ctx.BGSilhouetteAcceptanceProbability);
    
    glUseProgram(0);
    
    % Ok, all PTB managed buffers and shaders loaded and set up.
    % Lets create the VBO that we need to actually render anything in the
    % end. VBO's are not supported yet by PTB's Screen, so we need to
    % switch to our GL context for setup:
    SwitchToGL(ctx.parentWin);
    
    ctx.FGvbo = glGenBuffers(1);
    glBindBuffer(GL.ARRAY_BUFFER, ctx.FGvbo);

    % Calculate size of VBO in bytes: Number of potential foreground dots
    % times 4 components per dot (RGBA == xyzw) times 4 Bytes per float
    % component:
    buffersize = ctx.maxFGDots * 4 * 4;
            
    % Allocate but don't initialize it, ie NULL pointer == 0
    glBufferData(GL.ARRAY_BUFFER, buffersize, 0, GL.STREAM_COPY);
            
    % Done.
    glBindBuffer(GL.ARRAY_BUFFER, 0);
            
    % Setup another VBO for the vertex indices:
    ctx.FGibo = glGenBuffers(1);
    glBindBuffer(GL.ELEMENT_ARRAY_BUFFER_ARB, ctx.FGibo);

    % Allocate buffer for number of vertex indices,
    % each taking up 4 Bytes (== sizeof(uint32)) of memory.
    % Initialize immediately with indices and tell
    % OpenGL that this won't change at all during operation
    % (STATIC_DRAW):
    fgdotindices = uint32(0:ctx.maxFGDots-1);
    glBufferData(GL.ELEMENT_ARRAY_BUFFER_ARB, ctx.maxFGDots * 4, fgdotindices, GL.STATIC_DRAW);
    glBindBuffer(GL.ELEMENT_ARRAY_BUFFER_ARB, 0);

    ctx.BGvbo = glGenBuffers(1);
    glBindBuffer(GL.ARRAY_BUFFER, ctx.BGvbo);

    % Calculate size of VBO in bytes: Number of potential background dots
    % times 4 components per dot (RGBA == xyzw) times 4 Bytes per float
    % component:
    buffersize = ctx.maxBGDots * 4 * 4;
            
    % Allocate but don't initialize it, ie NULL pointer == 0
    glBufferData(GL.ARRAY_BUFFER, buffersize, 0, GL.STREAM_COPY);
            
    % Done.
    glBindBuffer(GL.ARRAY_BUFFER, 0);
            
    % Setup another VBO for the vertex indices:
    ctx.BGibo = glGenBuffers(1);
    glBindBuffer(GL.ELEMENT_ARRAY_BUFFER_ARB, ctx.BGibo);

    % Allocate buffer for number of vertex indices,
    % each taking up 4 Bytes (== sizeof(uint32)) of memory.
    % Initialize immediately with indices and tell
    % OpenGL that this won't change at all during operation
    % (STATIC_DRAW):
    bgdotindices = uint32(0:ctx.maxBGDots-1);
    glBufferData(GL.ELEMENT_ARRAY_BUFFER_ARB, ctx.maxBGDots * 4, bgdotindices, GL.STATIC_DRAW);
    glBindBuffer(GL.ELEMENT_ARRAY_BUFFER_ARB, 0);
    
    % Done with VBO setup.
    
    % Restore previous GL context binding:
    RestoreGL;

    % Set batch zero as starting batch:
    ctx.currentBatch = 0;

    % Clear vbosready state:
    ctx.vbosready = 0;
    
    % We're ready for the show!
    contextcount = contextcount + 1;

    % Init for this 'ctx' context done: Return it to usercode:
    varargout{1} = ctx;
    
    return;
end

% Assign color texture handle and target to context, enable texture mapping
% during 2D dot drawing:
if strcmpi(cmd, 'SetColorTexture')
    if nargin < 4
        error('In "SetColorTexture": You must provide the "context", "textureId" and "textureTarget"!');
    end

    % Get context object:
    ctx = varargin{1};
    
    % Get texture handle and target:
    ctx.colorTexId = varargin{2};
    ctx.colorTexTarget = varargin{3};
    
    % Empty or negative assignment resets to "no texture assigned":
    if isempty(ctx.colorTexId) || ctx.colorTexId < 0
        ctx.colorTexId = -1;
        ctx.colorTexTarget = -1;
    end

    % Return updated 'ctx' to usercode:
    varargout{1} = ctx;
    
    return;
end

% Assign shader handle for application of a GLSL shader during 2D dot
% drawing:
if strcmpi(cmd, 'SetDrawShader')
    if nargin < 3
        error('In "SetDrawShader": You must provide the "context" and shader handle!');
    end

    % Get context object:
    ctx = varargin{1};
    
    % Get texture handle and target:
    ctx.draw2DShaderFG = varargin{2};
    
    if nargin >= 4
        ctx.draw2DShaderBG = varargin{3};
    end
    
    % Optional sprite enable flag provided?
    if nargin >= 5 && ~isempty(varargin{4})
        ctx.needSprites = varargin{4};

        if ctx.needSprites ~= 1
            ctx.needSprites = 0;
        end
        
        % Enable or disable point-sprite coord generation on texture unit
        % 1, depending if point sprites shall be enabled or disabled:
        glActiveTexture(GL.TEXTURE1);
        
        if ctx.needSprites == 1
            % Enable point sprite coordinate generation on unit 1:
            glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.TRUE);
        else
            % Disable point sprite coordinate generation on unit 1:
            glTexEnvi(GL.POINT_SPRITE, GL.COORD_REPLACE, GL.FALSE);
        end

        glActiveTexture(GL.TEXTURE0);
    end
    
    % Empty or negative assignment resets to "no shader assigned":
    if isempty(ctx.draw2DShaderFG) || ctx.draw2DShaderFG <= 0
        % Detach shader:
        ctx.draw2DShaderFG = -1;
    end

    if isempty(ctx.draw2DShaderBG) || ctx.draw2DShaderBG <= 0
        % Detach shader:
        ctx.draw2DShaderBG = -1;
    end

    % Return updated 'ctx' to usercode:
    varargout{1} = ctx;
    
    return;
end

% Destroy processing context -- Release all ressources and shaders:
if strcmpi(cmd, 'DestroyContext')
    if nargin < 2
        error('In "DestroyContext": You must provide the "context" to destroy!');
    end

    % Get context object:
    ctx = varargin{1};
    
    % Delete all offscreen windows, gloperators and buffers like IBO's
    % VBO's, PBO's, FBO's etc...
    deleteContextBuffers(ctx);
        
    % Delete all shaders:
    glDeleteProgram(ctx.createBGDotsShader);
    glDeleteProgram(ctx.createFGDotsShader);
    glDeleteProgram(ctx.samplingShader);
    glDeleteProgram(ctx.trackingRenderShader);
    glDeleteProgram(ctx.silhouetteRenderShader);

    RestoreGL;
    
    % Shutdown done.
    contextcount = max(contextcount - 1, 0);
    
    % Return destroyed context:
    ctx = [];
    varargout{1} = ctx;

    return;
end

% Set string to call via feval() to initiate a 3D render cycle for the 3D
% scene/object to be visualized:
if strcmpi(cmd, 'SetRenderCallback')
    if nargin < 3
        error('In "SetRenderCallback": You must provide the "context" and callback string!');
    end

    % Get context object:
    ctx = varargin{1};
    
    % Get the eval string:
    renderCallback = varargin{2};
    if ~ischar(renderCallback)
        error('In "SetRenderCallback": Callback string must be a string, nothing else!');
    end
    
    % Assign:
    ctx.renderCallback = renderCallback;
    
    varargout{1} = ctx;

    return;
end

% Set state of engine to initial setting for a trial:
if strcmpi(cmd, 'ResetState')
    if nargin < 2
        error('In "ResetState": You must provide the "context"!');
    end

    % Get context object:
    ctx = varargin{1};
    
    BackupGL;
        
    % Reset to starting batch zero:
    ctx.currentBatch = 0;
    
    % Zero-out background sample matrix:
    ctx.BGSampleSet(:, :, :) = 0;
    
    % Clear out all buffers:
    SwitchToPTB;
    Screen('FillRect', ctx.BGDotsBuffer, [0 0 0 0])
    Screen('FillRect', ctx.FGDotsBuffer, [0 0 0 0])
    Screen('FillRect', ctx.trackingBuffer, [0 0 0 0])
    Screen('FillRect', ctx.silhouetteBuffer, [0 0 0 0])
    Screen('FillRect', ctx.sampleBuffer, [0 0 0 0])

    % Clear vbosready state:
    ctx.vbosready = 0;
    
    RestoreGL;
    
    varargout{1} = ctx;    
    return;
end

% Update cycle, possibly followed by a render operation:
if strcmpi(cmd, 'Update')
    if nargin < 2
        error(sprintf('In "%s": You must provide the "context"!', cmd)); %#ok<SPERR>
    end

    % Get context object:
    ctx = varargin{1};
    
    if nargin >= 3 && ~isempty(varargin{2})
        instantOn = varargin{2};
    else
        instantOn = 0;
    end
    
    BackupGL;
        
    % Switch to OpenGL rendering context to be used for 3D scene rendering,
    % and specifically for our silhouette render buffer:
    SwitchToPTB;
    
    Screen('BeginOpenGL', ctx.silhouetteBuffer);

    % Backup 3D context state:
    glPushAttrib(GL.ALL_ATTRIB_BITS);

    % Perform 1st 3D render pass:

    % Need zBuffer occlusion testing for silhouette rendering:
    glEnable(GL.DEPTH_TEST);
    
    % Need cleared buffer, including z buffer:
    glClearColor(0,0,0,0);
    glClear;

    % Bind shader for silhouette rendering:
    if debug~=-1
        % We skip this if debug flag == -1 -- In that case the user wants
        % to see the real rendered image instead of our silhouette
        % encoding.
        glUseProgram(ctx.silhouetteRenderShader);
    end
    
    % Set viewport and scissor to full target window area:
    glViewport(0, 0, ctx.silhouetteWidth, ctx.silhouetteHeight);
    glScissor(0, 0, ctx.silhouetteWidth, ctx.silhouetteHeight);

    % Call the render callback function in workspace of our caller. We did
    % not touch the modelview- or projection matrices, so the projections
    % et al. should be ok...
    evalin('caller', ctx.renderCallback);

    % Don't need depth testing anymore:
    glDisable(GL.DEPTH_TEST);

    Screen('EndOpenGL', ctx.silhouetteBuffer);
    % Silhouette should be ready in silhouette buffer...

    if abs(debug) == 1
        Screen('DrawTexture', ctx.parentWin, ctx.silhouetteBuffer, [], [], [], [], [], debugGain);
    end
    
    % Perform 2nd "pseudo 3D" render pass into trackingBuffer. This will
    % again render the geometry, but with different encoding. A unwrapped
    % texture map will be output, where each pixel corresponds to a surface
    % point on the 3D object (aka texture coordinate). The color of each
    % pixel encodes interpolated screen space (x,y,z) coordinates:
    Screen('BeginOpenGL', ctx.trackingBuffer);
    
    % No depth test here, as fragment depths doesn't encode anything
    % meaningful during this pass:
    glDisable(GL.DEPTH_TEST);
    
    % Bind shader for tracking image rendering:
    glUseProgram(ctx.trackingRenderShader);
    
    % Set viewport and scissor to full trackbuffer window area:
    glViewport(0, 0, ctx.texResolution(1), ctx.texResolution(2));
    glScissor(0, 0, ctx.texResolution(1), ctx.texResolution(2));
    
    % Call the render callback function in workspace of our caller. We did
    % not touch the modelview- or projection matrices, so the projections
    % et al. should be ok...
    evalin('caller', ctx.renderCallback);
    
    % Unbind all shaders:
    glUseProgram(0);

    % Just to make sure it's still off:
    glDisable(GL.DEPTH_TEST);

    % Restore 3D context state:
    glPopAttrib;

    % Trackingbuffer should be ready:
    Screen('EndOpenGL', ctx.trackingBuffer);

    if debug == 2
        Screen('DrawTexture', ctx.parentWin, ctx.trackingBuffer, [], [], [], [], [], debugGain);
    end
    
    % We are in Screen()'s rendering context. Do the 2D image processing
    % stuff:

    % Need to attach the silhouette Buffers
    % color buffer texture to texture units 1 and 2: As the texture is part of a
    % color buffer attachment, it is set to nearest neighbour sampling -
    % which is what we want:
    glActiveTexture(GL.TEXTURE1);
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, ctx.silhouetteTexture);
    glActiveTexture(GL.TEXTURE2);
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, ctx.silhouetteTexture);
    glActiveTexture(GL.TEXTURE0);    
    
    % --- BACKGROUND DOTS COMPUTATION ---
    
    % Perform update of background sample buffer with random samples:
    % Compute random sample locations in image via Matlabs/Octaves uniform
    % random number generator:
    if instantOn
        % Create new samples for full set:
        randomSamples = rand(size(ctx.BGSampleSet, 1), ctx.BGsamplesPerLine, 3);
    else
        % Create new samples for current batch:
        randomSamples = rand(ctx.BGsampleLinesPerBatch, ctx.BGsamplesPerLine, 3);
    end
    
    % Layers 1 and 2 contain properly scaled (x,y) screen coordinates of
    % our random "darts":
    randomSamples(:,:,1) = randomSamples(:,:,1) * ctx.silhouetteWidth;
    randomSamples(:,:,2) = randomSamples(:,:,2) * ctx.silhouetteHeight;
    
    % Layer 3 contains a uniformly distributed number between 0 and 1 for
    % use as random per-sample variable by internal random sampling...
    
    if instantOn
        % Assign new samples for full set:
        ctx.BGSampleSet(:, :, :) = randomSamples;
    else
        % Assign new samples for current batch:
        sline = ctx.currentBatch * ctx.BGsampleLinesPerBatch + 1;
        eline = sline + ctx.BGsampleLinesPerBatch - 1;
        ctx.BGSampleSet(sline:eline, :, :) = randomSamples;
    end

    % Background batch in background sample buffer updated. Convert whole
    % buffer to texture, with background sampling shader bound:
    BGsampleTex = Screen('MakeTexture', ctx.parentWin, ctx.BGSampleSet, [], [], 2, 0, ctx.createBGDotsShader);

    % Blit to background dots buffer:
    Screen('DrawTexture', ctx.BGDotsBuffer, BGsampleTex, [], [], [], 0);

    % Release texture:
    Screen('Close', BGsampleTex);

    % --- FOREGROUND DOTS COMPUTATION ---
        
    % Perform update of current batch of sampleBuffer:
    
    % Compute random sample locations in image via Matlabs/Octaves uniform
    % random number generator:
    if instantOn
        randomSamples = rand(ctx.sampleLinesPerBatch * ctx.dotLifetime, ctx.samplesPerLine, 2);
    else
        randomSamples = rand(ctx.sampleLinesPerBatch, ctx.samplesPerLine, 2);
    end
    
    randomSamples(:,:,1) = randomSamples(:,:,1) * ctx.silhouetteWidth;
    randomSamples(:,:,2) = randomSamples(:,:,2) * ctx.silhouetteHeight;

    % Create 32bpc float texture (setting '2') with the random sample
    % locations. If we'd use the texture inside Screen('TransformTexture')
    % we could use textureOrientation setting '3' for isotropic random noise.
    % This settings would save some texture conversion time then. However,
    % as we're only 'DrawTexture'ing and the used blit shaders are position
    % invariant (isotropic), a setting of 0 or 1 also works with no speed
    % loss, as 'DrawTexture' does implicit optimization:
    sampleTex = Screen('MakeTexture', ctx.parentWin, randomSamples, [], [], 2, 0, ctx.samplingShader);
    
    if debug == 3
        Screen('DrawTexture', ctx.parentWin, sampleTex, [], [], [], [], [], debugGain);
    end
    
    % Blit sampleTex into the target batch rectangle of our sampleBuffer,
    % with the sampling shader bound.
    
    if instantOn
        % Blit texture into sampleBuffer. The secondary texture unit provides
        % access to the silhouette image, the shader does bilinear filtering
        % and conversion:
        Screen('DrawTexture', ctx.sampleBuffer, sampleTex, [], [], 0, 0);        
    else
        % Blit texture at target location into sampleBuffer, offset vertically
        % so the proper batch gets updated. The secondary texture unit provides
        % access to the silhouette image, the shader does bilinear filtering
        % and conversion:
        Screen('DrawTexture', ctx.sampleBuffer, sampleTex, [], OffsetRect(Screen('Rect', sampleTex), 0, double(ctx.currentBatch * ctx.sampleLinesPerBatch)), 0, 0);
    end
    
    % Release sampleTex for next cycle:
    Screen('Close', sampleTex);

    if debug == 4
        Screen('DrawTexture', ctx.parentWin, ctx.sampleBuffer, [], [], [], [], [], debugGain);
    end
    
    % Our total distribution of foreground samples in sampleBuffer is now
    % up to date. Use stored surface (texcoords) sample locations to read
    % out corresponding locations in trackingBuffer, convert them into
    % vertex location definitions and blit them to the final buffer of
    % foreground dot specs. Screen('TransformTexture') will do the blit and
    % bind sampleBuffer and trackingBuffer as input textures:
    ctx.FGDotsBuffer = Screen('TransformTexture', ctx.sampleBuffer, ctx.createFGDotsoperator, ctx.trackingBuffer, ctx.FGDotsBuffer);
    
    % Ok, the ctx.FGDotsBuffer should contain the encoded set of all
    % foreground dot positions. One can either read this buffer back into a
    % Matlab/Octave matrix for usercode processing, or convert it into a
    % VBO via PBO mechanism, then render it.
    
    % Clear out all intermediate result buffers in preparation of next update cycle:
    Screen('FillRect', ctx.trackingBuffer, [0 0 0 0])
    Screen('FillRect', ctx.silhouetteBuffer, [0 0 0 0])

    if debug == 5
        Screen('DrawTexture', ctx.parentWin, ctx.FGDotsBuffer, [], [], [], [], [], debugGain);
    end
    
    % Clear vbosready state to trigger a refill on next 'Render':
    ctx.vbosready = 0;
    
    RestoreGL;
    
    if instantOn
        % Reset batch counter to zero after this initial "instant on" update cycle:
        ctx.currentBatch = 0;
    else
        % Increment batch counter for next update cycle:
        ctx.currentBatch = mod(ctx.currentBatch+1, ctx.dotLifetime);
    end
    
    % Ready for render: Return updated context:
    varargout{1} = ctx;
    
    return;
end

% Render current result in ctx.FGDotsBuffer into parent window:
if strcmpi(cmd, 'Render')
    if nargin < 2
        error(sprintf('In "%s": You must provide the "context"!', cmd)); %#ok<SPERR>
    end

    % Get context object:
    ctx = varargin{1};

    if nargin < 3 || isempty(varargin{2})
        targetWin = ctx.parentWin;
    else
        targetWin = varargin{2};
        if ~ismember(Screen('WindowKind', targetWin), [-1, 1])
            error('In "Render": Invalid "targetWindow" handle provided. Not an offscreen or onscreen window!');
        end
    end
    
    if nargin < 4
        drawspec = [1,1];
    else
        drawspec = varargin{3};
    end

    drawFG = drawspec(1);
    drawBG = drawspec(2);

    BackupGL;

    % Can do this in PTB's Screen 2D context, which is more convenient for
    % our 2D drawing operations, as long as we are careful to restore any
    % changed context state:
    SwitchToPTB;

    % Update of VBO's needed?
    if ~ctx.vbosready
        % Yes: Copy content of FGDotsBuffer into VBO, using PBO extension.
        % The 'GetWindowInfo' binds our ctx.FGDotsBuffer FBO so we can
        % glReadPixels() from it:
        Screen('GetWindowInfo', ctx.FGDotsBuffer);

        % There is a bug in the X1000 gfx-card driver on OS/X 10.4.11 which
        % causes glReadPixels() readback values to get clamped to 0-1 range
        % if alpha-blending is enabled. Therefore we need to disable alpha
        % blending during glReadPixels() readback and reenable later if
        % needed:
        alphaenabled = glIsEnabled(GL.BLEND);
        glDisable(GL.BLEND);
        
        glBindBuffer(GL.PIXEL_PACK_BUFFER_ARB, ctx.FGvbo);        
        glReadPixels(0, 0, ctx.samplesPerLine, ctx.sampleLinesTotal, GL.RGBA, GL.FLOAT, 0);
        glBindBuffer(GL.PIXEL_PACK_BUFFER_ARB, 0);
        
        Screen('GetWindowInfo', ctx.BGDotsBuffer);
        
        glBindBuffer(GL.PIXEL_PACK_BUFFER_ARB, ctx.BGvbo);        
        glReadPixels(0, 0, ctx.BGsamplesPerLine, ctx.BGsampleLinesTotal, GL.RGBA, GL.FLOAT, 0);
        glBindBuffer(GL.PIXEL_PACK_BUFFER_ARB, 0);

        % Reenable alpha blending if it was enabled:
        if alphaenabled
            glEnable(GL.BLEND);
        end
        
        % VBO's ready:
        ctx.vbosready = 1;        
    end
    
    % Setup render:

    % The 'GetWindowInfo' binds our ctx.parentWin so we can render to it:
    Screen('GetWindowInfo', targetWin);
        
    % Backup old 2D context state bits:    
    glPushAttrib(GL.ALL_ATTRIB_BITS);
    
    % Bind and enable vertex position VBO:
    glEnableClientState(GL.VERTEX_ARRAY);

    % Is point anti-aliasing enabled?
    pSmooth = glIsEnabled(GL.POINT_SMOOTH);
    
    if drawFG
        % Foreground render:
        glBindBuffer(GL.ARRAY_BUFFER, ctx.FGvbo);

        % Assign vertex pointer, setup proper stride for interleave with
        % texture coordinates from same VBO:
        glVertexPointer(2, GL.FLOAT, 4 * 4, 0);
        
        if ctx.colorTexId >= 0
            % Texture mapping for colored dot drawing:
            
            % Assign texture coord array, which is interleaved with vertex
            % coord array within FGvbo VBO: (x,y,tx,ty)....
            glEnableClientState(GL.TEXTURE_COORD_ARRAY);
            glTexCoordPointer(2, GL.FLOAT, 4 * 4, 2 * 4);
            
            % Enable texture mapping for proper target, and assign texture:
            glEnable(ctx.colorTexTarget);
            glBindTexture(ctx.colorTexTarget, ctx.colorTexId);
        end
        
        % Bind vertex index VBO:
        glBindBuffer(GL.ELEMENT_ARRAY_BUFFER_ARB, ctx.FGibo);

        % Enable draw GLSL shader, if any:
        if ctx.draw2DShaderFG > -1
            glUseProgram(ctx.draw2DShaderFG);
            
            if pSmooth && ctx.needSprites
                glEnable(GL.POINT_SPRITE_ARB);
            end
        end
        
        % Perform draw operation: All vertices, each triggering render for a
        % single GL.POINT primitive. Colors, sizes, anti-aliasing flags etc.
        % can be set from external code as appropriate. Application of textures
        % or shaders is also possible:
        glDrawRangeElements(GL.POINTS, 0, ctx.maxFGDots-1, ctx.maxFGDots, GL.UNSIGNED_INT, 0);
        
        % Disable draw GLSL shader, if any:
        if ctx.draw2DShaderFG > -1
            glUseProgram(0);

            if pSmooth && ctx.needSprites
                glDisable(GL.POINT_SPRITE_ARB);
            end
        end

        if ctx.colorTexId >= 0
            % Disable interleaved texturecoord array:
            glDisableClientState(GL.TEXTURE_COORD_ARRAY);

            % Disable texture mapping for proper target, and assign texture:
            glDisable(ctx.colorTexTarget);
            glBindTexture(ctx.colorTexTarget, 0);
        end        
    end
    
    if drawBG
        % Background render:
        glBindBuffer(GL.ARRAY_BUFFER, ctx.BGvbo);
        glVertexPointer(4, GL.FLOAT, 0, 0);

        % Bind vertex index VBO:
        glBindBuffer(GL.ELEMENT_ARRAY_BUFFER_ARB, ctx.BGibo);

        % Enable draw GLSL shader, if any:
        if ctx.draw2DShaderBG > -1
            glUseProgram(ctx.draw2DShaderBG);

            if pSmooth && ctx.needSprites
                glEnable(GL.POINT_SPRITE_ARB);
            end
        end
        
        % Perform draw operation: All vertices, each triggering render for a
        % single GL.POINT primitive. Colors, sizes, anti-aliasing flags etc.
        % can be set from external code as appropriate. Application of textures
        % or shaders is also possible:
        glDrawRangeElements(GL.POINTS, 0, ctx.maxBGDots-1, ctx.maxBGDots, GL.UNSIGNED_INT, 0);
        
        % Disable draw GLSL shader, if any:
        if ctx.draw2DShaderBG > -1
            glUseProgram(0);

            if pSmooth && ctx.needSprites
                glDisable(GL.POINT_SPRITE_ARB);
            end        
        end
    end
    
    % Unbind our VBOs:
    glBindBuffer(GL.ELEMENT_ARRAY_BUFFER_ARB, 0);
    glBindBuffer(GL.ARRAY_BUFFER, 0);

    % Disable vertex array:
    glDisableClientState(GL.VERTEX_ARRAY);
    
    % Restore old 2D context state bits:
    glPopAttrib;
    
    % Render completed. Restore pre-render state:
    RestoreGL;
    
    % Rendering done: Return updated context:
    varargout{1} = ctx;
    
    return;
end

if strcmpi(cmd, 'GetResults')
    if nargin < 2
        error(sprintf('In "%s": You must provide the "context"!', cmd)); %#ok<SPERR>
    end

    % Get context object:
    ctx = varargin{1};

    BackupGL;
    
    % Readback as "image matrix":

    % Can do this in Screen-Mode:
    SwitchToPTB;
    
    % 'GetWindowInfo' selects ctx.FGDotsBuffer as active FBO, so we can
    % readback from it:
    Screen('GetWindowInfo', ctx.FGDotsBuffer);
    
    % There is a bug in the X1000 gfx-card driver on OS/X 10.4.11 which
    % causes glReadPixels() readback values to get clamped to 0-1 range
    % if alpha-blending is enabled. Therefore we need to disable alpha
    % blending during glReadPixels() readback and reenable later if
    % needed:
    alphaenabled = glIsEnabled(GL.BLEND);
    glDisable(GL.BLEND);

    % Readback to matrix: Cast from float aka single() type to double() type:
    readbackdata = double(glReadPixels(0, 0, ctx.samplesPerLine, ctx.sampleLinesTotal, GL.RGBA, GL.FLOAT));
    
    % Return readback (x,y) FG dot locations, reshaped into a 2 rows array
    % with one column per (x,y) dot:
    varargout{1} = [ reshape(readbackdata(:,:,1), 1, ctx.samplesPerLine * ctx.sampleLinesTotal) ; reshape(readbackdata(:,:,2), 1, ctx.samplesPerLine * ctx.sampleLinesTotal) ];

    % Optionally also return texture coordinates (u,v) of FG dots as 3rd output arg:
    if nargout > 2
        varargout{3} = [ reshape(readbackdata(:,:,3), 1, ctx.samplesPerLine * ctx.sampleLinesTotal) ; reshape(readbackdata(:,:,4), 1, ctx.samplesPerLine * ctx.sampleLinesTotal) ];
    end
    
    % Return optional readback (x,y) BG dot locations, reshaped into a 2 rows array
    % with one column per (x,y) dot as 2nd output arg:
    if nargout > 1
        Screen('GetWindowInfo', ctx.BGDotsBuffer);

        % Readback to matrix: Cast from float aka single() type to double() type:
        readbackdata = double(glReadPixels(0, 0, ctx.BGsamplesPerLine, ctx.BGsampleLinesTotal, GL.RGBA, GL.FLOAT));

        % Return readback (x,y) dot locations, reshaped into a 2 rows array
        % with one column per (x,y) dot:
        varargout{2} = [ reshape(readbackdata(:,:,1), 1, ctx.BGsamplesPerLine * ctx.BGsampleLinesTotal) ; reshape(readbackdata(:,:,2), 1, ctx.BGsamplesPerLine * ctx.BGsampleLinesTotal) ];
    end
    
    % Reenable alpha blending if it was enabled:
    if alphaenabled
        glEnable(GL.BLEND);
    end
    
    % Readback completed. Restore pre-readback state:
    RestoreGL;
    
    return;
end

if strcmpi(cmd, 'DebugFlag')
    if nargin < 2
        error('Must provide new setting for debug flag!');
    end
    
    varargout{1} = debug;
    varargout{2} = debugGain;
    
    debug = varargin{1};
    
    if nargin < 3
        debugGain = [];
    else
        if length(varargin{2})~=4
            error('In "DebugFlag": "debugGain" color modulation gain must be a 4 element [R,G,B,A] modulation color vector!');
        end
        
        debugGain = varargin{2};
    end
    return;
end

error(sprintf('Invalid subcommand ''%s'' specified!', cmd)); %#ok<SPERR>
return;

% Internal helper functions:
function SwitchToGL(win)

% Switch to our OpenGL context, but keep a backup of original
% drawstate. We do lazy switching if possible:
[currentwin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');

if ~IsOpenGLRendering
    % PTB's context active: Switch to OpenGL rendering for our parent window:
    Screen('BeginOpenGL', win);
else
    % In rendering context. Is it ours? If yes, then there isn't anything
    % to do...
    if currentwin ~= win
        % No, a different windows context is active: First switch to PTB
        % mode, then switch to ours:
        
        % Switch to our parentWin's PTB context:
        Screen('EndOpenGL', currentwin);
        % Switch to our parentWin's GL context:
        Screen('BeginOpenGL', win);
    end
end
return;

function SwitchToPTB

% Switch from our OpenGL context, but keep a backup of original
% drawstate. We do lazy switching if possible:
[currentwin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');

if ~IsOpenGLRendering
    % PTB's context is already active: Nothing to do.
else
    % In rendering context. Switch back to PTB - and to our parentWin:
    Screen('EndOpenGL', currentwin);
end
return;

function BackupGL
global moglFDF_OriginalContext;

if ~isempty(moglFDF_OriginalContext)
    error('BackupGL called twice in succession without intermediate RestoreGL! Ordering inconsistency!');
end

[currentwin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');

if IsOpenGLRendering
    moglFDF_OriginalContext = currentwin;
end
return;

function RestoreGL
global moglFDF_OriginalContext;

[currentwin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');

if isempty(moglFDF_OriginalContext)
    % PTB was in Screen drawing mode: Switch to that mode, if not active:
    if IsOpenGLRendering
        Screen('EndOpenGL', currentwin);
    end
    return;
end

% Need to restore to GL context if not already active:
if ~IsOpenGLRendering
    Screen('BeginOpenGL', moglFDF_OriginalContext);
else
    % OpenGL context active. Ours? If so -> Nothing to do.
    if currentwin ~= moglFDF_OriginalContext
        % Nope. Need to switch:
        Screen('EndOpenGL', currentwin);
        Screen('BeginOpenGL', moglFDF_OriginalContext);        
    end
end

% Restore to default:
moglFDF_OriginalContext = [];

return;

function deleteContextBuffers(ctx)
    BackupGL;
    
    SwitchToGL(ctx.parentWin);
    
    % Delete VBO's:
    glDeleteBuffers(1, ctx.FGibo);
    glDeleteBuffers(1, ctx.FGvbo);
    glDeleteBuffers(1, ctx.BGibo);
    glDeleteBuffers(1, ctx.BGvbo);
    
    SwitchToPTB;
    
    % Close all offscreen windows and their associated textures:
    Screen('Close', [ctx.BGDotsBuffer, ctx.FGDotsBuffer, ctx.trackingBuffer, ctx.silhouetteBuffer, ctx.sampleBuffer]);
    
    % Close our operators:
    Screen('Close', ctx.createFGDotsoperator);
    
return;
