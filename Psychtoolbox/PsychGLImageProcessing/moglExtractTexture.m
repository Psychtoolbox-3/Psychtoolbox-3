function varargout = moglExtractTexture(cmd, varargin)
% moglExtractTexture(cmd [, arg1][, arg2][, ...]) - "MOGL Video texture extraction"
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
%
% Subcommands, their syntax & meaning:
% ====================================
%
% [oldflag, oldgain] = moglExtractTexture('DebugFlag', flag [, debugGain]);
% - Set debug flag to value 'flag'. Default is zero. Non-zero values enable
% different visualizations that may aid debugging non-working setups. 1 =
% Show silhouette buffer, 2 = Show trackbuffer, 3 = Show extracted texture.
% A setting of -1 shows the real rendered image. A value of -2 disables any
% kind of textual warnings.
%
% The optional 'debugGain' parameter must be a 4 component [R G B A] color
% vector with modulation gains for the drawn "debug images" - simply to
% scale each color channel in intensity to allow for display of values
% outside the standard displayable range between zero and one.
%
%
% context = moglExtractTexture('CreateContext', window, rect, texCoordMin, texCoordMax, texResolution [,zThreshold=Off]); 
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
% the size and shape of the input video texture. This rect must have the
% same size as the input video image textures -- Lots of internal
% calculations depend on this geometry spec!
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
% longer processing times. This defines the size of returned extracted
% textures.
%
% 'zThreshold' Optional zThreshold for occlusion test: By default, it is
% 10.0 ie. occlusion test disabled. A value between 0.0 and 1.0 will enable
% occlusion testing -- Texels that would correspond to occluded surface patches are
% not extracted. Small numbers (close to zero) make the test more sensitive but
% can cause artifacts due to internal numeric roundoff errors. Bigger
% numbers (closer to one) make it more robust but less powerful. The
% "sweet-spot" depends on your hardware and 3D scene. Empirically a setting
% of 0.0001 is a good value for ATI Radeon X1000 series hardware.
% The default setting (bigger than 1.0) will disable occlusion test --
% "Hidden texels" are not ignored, but updated with bogus extracted texture.
%
%
% context = moglExtractTexture('SetRenderCallback', context, callbackEvalString);
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
% context = moglExtractTexture('DestroyContext', context);
% - Destroy a processing context, release all of its ressources.
%
%
% [texBuffer, texId, texTarget] = moglExtractTexture('Extract', context, inputTexture [, newTexture = 0]);
% - Perform an 'Extract' cycle for given context. A new "3D frame" is rendered
% via the rendercallback function, then analysed, to provide the 3D surface
% geometry and occlusion info and mapping for texture extraction. This info
% is then used to extract pixel color values from the given video input
% texture 'inputTexture' and the final extracted texturemap is stored
% inside an internal texture buffer. A handle to that internal buffer
% 'texBuffer' is returned. The handle is owned by this function! You should
% not close or otherwise mess with the provided buffer. You can read the
% final texture from it, acquire a temporary OpenGL texture handle to it
% for texture mapping, etc. You are even allowed to perform destructive
% write informations on the buffer to change its pixel content. But do not
% destory and reallocate the buffer, change its size, number of layers,
% resolution or any other property! For your convenience, 'texId' and
% 'texTarget' also provide standard OpenGL handles to texture id and
% target, associated with 'texBuffer'.
%
% Alternatively you can set the optional flag 'newTexture' to a value of 1.
% In that case, a new extracted texture 'texBuffer' is returned and you own
% this texture, ie., you can do with it whatever you want and you are
% responsible for releasing the texture via Screen('Close', texBuffer);
% once you are done with it.
%

% History:
%  05/10/09  Initial implementation, derived from moglFDF (MK).

% Need OpenGL constants:
global GL;

% Internal state:
global moglExtractTexture_OriginalContext;
persistent contextcount;
persistent debug;
persistent debugGain;

if nargin < 1
    error('You must provide a "cmd" subcommand to execute!');
end

if isempty(contextcount)
    contextcount = 0;
    moglExtractTexture_OriginalContext = [];
    debug = 0;
    
    if isempty(GL)
        % If OpenGL not initialized, do a full init for 3D mode:
        error('OpenGL mode not initialized! You *must* call InitializeMatlabOpenGL before the first call to this routine or any Screen() function!')
    end
end

% Subcommand dispatch:

% Initialization of a new context: Allocate and setup all ressources:
if strcmpi(cmd, 'CreateContext')
    % Fetch all arguments - They are all required.
    if nargin < 6
        error(sprintf('Some mandatory input arguments to "%s" are missing. Please provide them!', cmd)); %#ok<SPERR>
    end
    
    createContext = 1;
    
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
            error('Invalid "context" argument provided to "ReinitContext" - This is not a valid moglExtractTexture context handle!');
        end
        
        if ~isfield(ctx, 'moglExtractTextureMagic')
            disp(ctx);
            error('Invalid "context" argument provided to "ReinitContext" - This is not a valid moglExtractTexture context handle!');
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
    ctx.moglExtractTextureMagic = 'Funky magic-cookie';
    
    % Get all other arguments and perform parameter type and range checks:
    ctx.rect = varargin{2};
    if ~isnumeric(ctx.rect) || length(ctx.rect)~=4
        disp(ctx.rect);
        error('Invalid "rect" argument provided to "CreateContext" - Must be a 4 component vector that describes the size and shape of the input video rectangle [left top right bottom]');        
    end
    
    ctx.rect = double(ctx.rect);
    if IsEmptyRect(ctx.rect)
        disp(ctx.rect);
        error('Invalid "rect" argument provided to "CreateContext" - Must be a non-empty rect that describes the size and shape of the input video rectangle [left top right bottom]');        
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
    
    % Basic checks passed: Now check for inter-parameter consistency:    
    if nargin >= 7
        ctx.zThreshold = varargin{6};
        if ~isscalar(ctx.zThreshold)
            disp(ctx.zThreshold);
            error('Invalid "zThreshold" argument provided to "CreateContext" - Must be a positive number in range 0.0 - 1.0 for z-Test, and bigger for z-Test disabled!');
        end
    else
        % Default: z-Test disabled:
        ctx.zThreshold = 10.0;
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
        if ~(   ~isempty(findstr(glGetString(GL.EXTENSIONS), '_framebuffer_object')) && ...
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
    
    % Silhouette buffer: Contains the "perspective correct image space"
    % image of the rendered 3D object. The object is rendered "normally" in
    % a first render pass to get its silhouette, except that the color of
    % each rendered pixel is not a shaded/lit color, but its encoded
    % interpolated surface texture coordinate. This will be done by a
    % proper fragment shader during render pass.
    %
    % The buffer is an FBO backed offscreen window the same size as the
    % parent window, but with 32bpc float format to store accurate texture
    % coordinates in the pixel colors. Red channel encodes s-coord, Green
    % channel encodes t-coord of 2D texture coordinate, blue encodes a
    % foreground/background flag, alpha encodes z-buffer depths.
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
    % lookup the (X,Y) image position of a texel in the objects surface.
    %
    % Again a 32bpc float offscreen window FBO, but the resolution is
    % chosen per user spec to be fine enough in texture coordinate space to
    % match the size of the extracted texture map:
    ctx.trackingBuffer = Screen('OpenOffscreenWindow', ctx.parentWin, [0 0 0 0], [0, 0, ctx.texResolution(1), ctx.texResolution(2)], 128, 32);

    % Final buffer with extraced texture image. This one will get filled
    % by the texture extraction shader. A RGBA8 texture format is
    % sufficient, therefore we allocate a standard 32 bit surface:
    ctx.OutTextureBuffer = Screen('OpenOffscreenWindow', ctx.parentWin, [0 0 0 0], [0, 0, ctx.texResolution(1), ctx.texResolution(2)], 32, 32);
    
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

        % Shader for final creation of foreground dots VBO spec from
        % distribution in sampleBuffer and trackingBuffer:
        ctx.textureExtractionShader = LoadGLSLProgramFromFiles([shaderpath 'moglTextureExtractionShader'], 1);        
    end
    
    % Setup trackingRenderShader:
    glUseProgram(ctx.trackingRenderShader)

    % Compute texture coordinate offset and multiplier to apply in order to
    % remap the real texture coordinate range into the normalized [-1:+1]
    % interval that can pass through vertex clipping:
    glUniform4f(glGetUniformLocation(ctx.trackingRenderShader, 'TextureOffsetBias'), ctx.texCoordMin(1), ctx.texCoordMin(2), 2.0/(ctx.texCoordMax(1) - ctx.texCoordMin(1)), 2.0/(ctx.texCoordMax(2) - ctx.texCoordMin(2)));

    % Set viewport dimensions:
    glUniform4f(glGetUniformLocation(ctx.trackingRenderShader, 'Viewport'), 0, 0, ctx.silhouetteWidth/2, ctx.silhouetteHeight/2);
        
    % Setup shader for final texture extraction from input image and
    % xformed geometry and silhouette:
    glUseProgram(ctx.textureExtractionShader)

    % Bind texunit 0 to object coordinates texture:
    glUniform1i(glGetUniformLocation(ctx.textureExtractionShader, 'GeometryBuffer'), 0);

    % Bind texunit 1 to input image texture:
    glUniform1i(glGetUniformLocation(ctx.textureExtractionShader, 'InputImage'), 1);

    % Bind texunit 2 to silhouette texture for last rendered frame:
    glUniform1i(glGetUniformLocation(ctx.textureExtractionShader, 'Silhouette'), 2);
    
    % Assign zThreshold for depths testing of foreground dots before
    % output to handle occlusions correctly:
    glUniform1f(glGetUniformLocation(ctx.textureExtractionShader, 'zThreshold'), ctx.zThreshold);

    % Assign height of final output window + 1 to allow shader to invert
    % y-coordinate of final dots properly to account for difference in
    % y-axis direction of Screen()'s reference frame vs. OpenGL default
    % frame:
    % TODO FIXME: Needed or not?
    glUniform1f(glGetUniformLocation(ctx.textureExtractionShader, 'ViewportHeight'), ctx.silhouetteHeight + 1);
    
    % Define size of GeometryBuffer -- wrapAround values for interpolated texture lookup coordinates:
    glUniform2f(glGetUniformLocation(ctx.textureExtractionShader, 'texWrapAround'), ctx.texResolution(1), ctx.texResolution(2));

    % Define inverse remapping of texture coordinates into range
    % 0-texResolution -- The size of the trackingBuffer. N.B.: A neutral
    % mapping would be (0, 0, 1, 1) - That would pass trackingBuffer
    % texture coordinates instead of object texture coordinates.
% TODO delete    glUniform4f(glGetUniformLocation(ctx.textureExtractionShader, 'TextureOffsetBias'), ctx.texCoordMin(1), ctx.texCoordMin(2), 1 / (ctx.texResolution(1)/(ctx.texCoordMax(1) - ctx.texCoordMin(1))), 1 / (ctx.texResolution(2)/(ctx.texCoordMax(2) - ctx.texCoordMin(2))));

    glUseProgram(0);
    
    % Create gloperator from shader for later use by Screen('TransformTexture'):
    ctx.textureExtractionOperator = CreateGLOperator(ctx.parentWin, [], ctx.textureExtractionShader, 'Extract Texture.');
        
    % Ok, all PTB managed buffers and shaders loaded and set up.
    % Lets create the VBO that we need to actually render anything in the
    % end. VBO's are not supported yet by PTB's Screen, so we need to
    % switch to our GL context for setup:
    % DODO FIXME: Needed? Delete!
    SwitchToGL(ctx.parentWin);
    
    % Restore previous GL context binding:
    RestoreGL;
    
    % We're ready for the show!
    contextcount = contextcount + 1;

    % Init for this 'ctx' context done: Return it to usercode:
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
    glDeleteProgram(ctx.textureExtractionShader);
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

% Update cycle, possibly followed by a render operation:
if strcmpi(cmd, 'Extract')
    if nargin < 2
        error(sprintf('In "%s": You must provide the "context"!', cmd)); %#ok<SPERR>
    end

    if nargin < 3
        error(sprintf('In "%s": You must provide the "inputImage"!', cmd)); %#ok<SPERR>
    end

    if nargin >= 4
        doNotRecycle = varargin{3};
    else
        doNotRecycle = [];
    end
    
    % Do recycle texture extraction buffer by default:
    if isempty(doNotRecycle)
        doNotRecycle = 0;
    end
    
    % Get context object:
    ctx = varargin{1};
    
    % Get texture handle of input video texture:
    inputImage = varargin{2};

    % Backup current OpenGL state:
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
    glActiveTexture(GL.TEXTURE2);
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, ctx.silhouetteTexture);
    glActiveTexture(GL.TEXTURE0);    
    
    % Screen('TransformTexture') will apply the texture extraction operator
    % to the inputImage video image texture and write the extracted texture
    % into OutTextureBuffer. It will automatically bind inputImage and
    % trackingBuffer as input textures to units 0 and 1, whereas we
    % manually bound the silhouetteTexture to unit 2:
    if doNotRecycle
        recycleBuffer = [];
    else
        % Clear texture extraction result buffer:
        Screen('FillRect', ctx.OutTextureBuffer, [0 0 0 0])

        % Assign this buffer for recycling aka rewrite:
        recycleBuffer = ctx.OutTextureBuffer;
    end
    
    extractedTexture = Screen('TransformTexture', ctx.trackingBuffer, ctx.textureExtractionOperator, inputImage, recycleBuffer);
    
    if ~doNotRecycle
        ctx.OutTextureBuffer = extractedTexture;
    end

    % Ok, the ctx.OutTextureBuffer should contain the extracted texture.
    
    % Clear out all intermediate result buffers in preparation of next extraction cycle:
    Screen('FillRect', ctx.trackingBuffer, [0 0 0 0])
    Screen('FillRect', ctx.silhouetteBuffer, [0 0 0 0])

    if debug == 3
        Screen('DrawTexture', ctx.parentWin, extractedTexture, [], [], [], [], [], debugGain);
    end
    
    % Restore previous OpenGL context state:
    RestoreGL;
    
    % Ready: Return handle to our extracted texture buffer:
    varargout{1} = extractedTexture;
    [varargout{2}, varargout{3}] = Screen('GetOpenGLTexture', ctx.parentWin, extractedTexture);

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
global moglExtractTexture_OriginalContext;

if ~isempty(moglExtractTexture_OriginalContext)
    error('BackupGL called twice in succession without intermediate RestoreGL! Ordering inconsistency!');
end

[currentwin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');

if IsOpenGLRendering
    moglExtractTexture_OriginalContext = currentwin;
end
return;

function RestoreGL
global moglExtractTexture_OriginalContext;

[currentwin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');

if isempty(moglExtractTexture_OriginalContext)
    % PTB was in Screen drawing mode: Switch to that mode, if not active:
    if IsOpenGLRendering
        Screen('EndOpenGL', currentwin);
    end
    return;
end

% Need to restore to GL context if not already active:
if ~IsOpenGLRendering
    Screen('BeginOpenGL', moglExtractTexture_OriginalContext);
else
    % OpenGL context active. Ours? If so -> Nothing to do.
    if currentwin ~= moglExtractTexture_OriginalContext
        % Nope. Need to switch:
        Screen('EndOpenGL', currentwin);
        Screen('BeginOpenGL', moglExtractTexture_OriginalContext);        
    end
end

% Restore to default:
moglExtractTexture_OriginalContext = [];

return;

function deleteContextBuffers(ctx)
    BackupGL;
         
    SwitchToPTB;
    
    % Close all offscreen windows and their associated textures:
    Screen('Close', [ctx.OutTextureBuffer, ctx.trackingBuffer, ctx.silhouetteBuffer]);
    
    % Close our operators:
    Screen('Close', ctx.textureExtractionOperator);
    
return;
