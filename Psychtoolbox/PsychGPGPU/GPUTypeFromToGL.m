function outObj = GPUTypeFromToGL(cmd, inObj, glObjType, outObj, keepmapped, mapflags)
% outObj = GPUTypeFromToGL(cmd, inObj [, glObjType][, outObj][, keepmapped][, mapflags])
%
% Note: Calling this command requires calling the following command first
% to initialize Psychtoolbox GPU computing support:
%
% PsychImaging('AddTask', 'General', 'UseGPGPUCompute', 'GPUmat');
%
% Supported 'cmd' commands:
% -------------------------
%
% if cmd is zero, then convert an OpenGL object of type glObjType,
% referenced by handle inObj into a GPU object and return it in outObj. If
% the optional outObj is provided as input argument, try to recycle it --
% just fill its content with OpenGL object's content. Otherwise, create an
% outObj of matching format for content. If 'keepmapped' is set to 1, the OpenGL
% object will stay mapped for the GPU compute api, otherwise it gets immediately
% unmapped after the conversion. Keeping the object mapped is more efficient, but
% requires more careful management of objects to prevent malfunctions.
%
% If cmd is == 1, then convert GPU object inObj to OpenGL object of type
% glObjType and return it in outObj. Try to recycle a passed in outObj, if
% possible, otherwise create a new one. 'keepmapped' - see explanation for cmd zero.
%
% If cmd is == 2, then unmap the OpenGL object. You must do this if you previously
% set the optional 'keepmapped' flag to 1 during a copy operation and now want to
% use the object which was the source or the target of that copy again with OpenGL
% or Screen(), ie., with a Psychtoolbox drawing or image processing function.
% Unmapping is neccessary for proper OpenGL operation, but costs a fraction of a
% millisecond of overhead on well working operating systems like Linux. Clever use
% of the 'keepmapped' flag and this manual unmapping method sometimes allows to
% save some redundant unmap calls.
%
% If cmd is == 3, then remove the OpenGL object from use by the GPU compute toolkit.
% This must be done before destroying/deleting the OpenGL object, e.g., before
% a call to Screen('Close', x); for a window or texture handle x. This operation
% can be very expensive -- on the order of multiple milliseconds, so use sparingly.
%
% If cmd is == 4, all OpenGL objects are removed. Usually used before closing (all)
% onscreen windows, e.g., via Screen('CloseAll') or sca. This cache flush is very
% expensive!
%
% If cmd is == 5, then the given OpenGL object 'inObj' of type glObjType is mapped
% and a CUDA memory pointer is returned, for use with external mex files, so these
% can directly access the mapped resource. The object is mapped read-only.
%
% If cmd is == 6, the same operation as cmd == 5 happens, but the object is mapped
% write-only.

% glObjType == 0 (default): Provided OpenGL object is a Psychtoolbox
% texture or offscreen window handle.
%
% glObjType == 1: Provided inObj is a struct which defines the low-level
% OpenGL object, which can be a texture or a renderbuffer. The struct must
% have the following fields:
%
%     texstruct.glhandle     == OpenGL object handle.
%     texstruct.gltarget     == Target: texture target or renderbuffer.
%     texstruct.width        == Width in texels/pixels.
%     texstruct.height       == Height in texels/pixels.
%     texstruct.bpp          == Bytes per texel/pixel/element.
%     texstruct.nrchannels   == Number of layers / color channels.
%
% -> glObjType 1 can be used in places where no calls to Screen() functions
%    are allowed or possible, e.g., inside the imaging pipeline, or 3rd
%    party low-level OpenGL code.
%
%
% Note: If you pass in a Psychtoolbox texture, it should be already in
% normalized orientation (upright and in row-major format). This is a given
% if the texture was created via Screen('SetOpenGLTexture') or
% Screen('SetOpenGLTextureFromMemPointer'); or if your texture is actually
% an offscreen window created via Screen('OpenOffscreenWindow'). If your
% texture is created via Screen('MakeTexture') you need to usually set the
% optional 'textureOrientation' flag to 1, unless you've pretransposed the
% Matlab/octave image matrix (setting of 2 is fine), or it is entirely
% isotropic (setting of 3 is fine). If you get your texture from a movie
% file, you need to pass the optional 'specialFlags1' parameter in
% Screen('OpenMovie') as 16.
%
% If you get the texture from the video capture engine, you need to pass
% the optional 'recordingflags' to 2048 in a call to
% Screen('OpenVideoCapture').
%
%
% glObjType == 2: Read or write from/to current virtual framebuffer for a
% given onscreen window handle passed as 'inObj', specifically the 1st
% color attachment.
%
%
% glObjType == 3: Read or write from/to currently bound FBO, specifically
% the 1st color attachment. 'inObj' or 'outObj' doesn't really have a
% meaning here, as we always query the current binding.
%
%
% Current Limitations:
%
% Currently only supports the GPUmat toolbox as hard-coded backend:
% (http://sourceforge.net/projects/gpumat/)
%
% In the future it should support more GPGPU backends and allow dynamic
% detection and/or runtime selection of backends. Possible candidates are,
% e.g., AccelerEyes "Jacket", low-level CUDA or OpenCL, and other toolkits
% based on CUDA or OpenCL, as well as our own to-be-done backed.
%
% Only really supports 32 bpc floating point precision textures and
% renderbuffers. This because this single precision float format is the
% only format common to both OpenGL and our one and only GPUmat backend.
% One can provide RGBA8 4-layer textures/renderbuffers, but these will be
% interpreted by the backend as single layer (luminance) single precision
% float matrix. Special CUDA kernels would be required in GPUmat to depack
% each apparent float pixel into a RGBA8 interleaved pixel for meaningful
% processing. Otherwise hilarious results will ensue.
%
% CUDA-5.0 interop as used by GPUmat currently only supports 1-layer,
% 2-layer and 4-layer textures and renderbuffers, ie., L, LA and RGBA, but
% not RGB format.
%

% History:
% 30.01.2013  mk  Written.
% 15.04.2013  mk  Require use of PsychImaging(..., 'UseGPGPUCompute', ...);
%

global GL;
persistent initialized;

% This global variable signals if a GPGPU compute api is enabled, and which
% one. It gets initialized by PsychImaging() if usercode requests GPGPU
% compute support: 0 = None, 1 = GPUmat.
global psych_gpgpuapi;

if isempty(initialized)
    % Make sure GPGPU computing got enabled by PsychImaging and GPU api
    % type 1, the GPUmat toolbox, is in use:
    if isempty(psych_gpgpuapi) || (psych_gpgpuapi ~= 1)
        error('GPGPU computing via GPUmat toolbox not enabled! Aborted.');
    end
    
    if isempty(GL)
        InitializeMatlabOpenGL([], [], 1);
    end
    
    % Ready to rock!
    initialized = 1;
end

if nargin < 1 || isempty(cmd) || ~isscalar(cmd) || ~isnumeric(cmd)
    error('Missing or invalid minimum required argument "cmd".');
end

switch cmd
    case 0,
        % Copy from OpenGL to GPU backend:
        direction = 0;
    
    case 1,
        % Copy from GPU backend to OpenGL:
        direction = 1;
        
    case 2,
        % Unmap object from cache:
        outObj = 1;
        direction = 0;

    case 3,
        % Purge object from cache:
        outObj = 1;
        direction = 0;
        
    case 4,
        % Cache invalidate:
        outObj = 1;
        memcpyCudaOpenGL(0);
        return;
    
    case 5
        % Retrieve mapped pointers for reading from OpenGL:
        direction = 0;
        
    case 6
        % Retrieve mapped pointers for writing to OpenGL:
        direction = 0;
        
    otherwise
        error('Invalid cmd specified.');
end

if nargin < 2
    error('Missing required 2nd argument "inObj".');
end

if nargin < 3 || isempty(glObjType)
    % Psychtoolbox "classic" texture handle or offscreen window handle:
    glObjType = 0;
end

% No outObj provided for recycling?
if nargin < 4
    outObj = [];
end

if nargin < 5 || isempty(keepmapped)
    keepmapped = 0;
end

if nargin < 6 || isempty(mapflags)
    mapflags = -1;
end

if direction == 0
    % OpenGL -> GPU => gpu is outObj, if any.
    gpu = outObj;
else
    % GPU -> OpenGL => gpu is inObj:
    gpu = inObj;
    if isempty(gpu)
        error('Empty GPUtype inObj variable provided, but update of OpenGL object requested. How is this supposed to work?!?');
    end
    
    % Impedance matching code. Try to massage input 'gpu' variable into a
    % format that is compatible with CUDA-OpenGL interop and OpenGL itself:
    switch ndims(gpu)
        case 1,
            % One dimensional vector: Turn into "row-vector"
            % style single texel row luminance texture:
            d1 = 1;
            d2 = numel(gpu);
            cc = 1;
            % Reshape into 3-D matrix with two singleton dimensions
            % to checks further down don't fail due to
            % size/dimension mismatch:
            gpu = reshape(gpu, 1, d1, d2);
        case 2,
            % Two dimensional matrix: Turn into luminance texture,
            % the 2D size of the two non-singleton dimensions:
            d1 = size(gpu, 2);
            d2 = size(gpu, 1);
            cc = 1;
            
            % Reshape gpu into a 3D matrix with the 1st dimension
            % being a singleton dimension which represents the
            % single luminance channel. We need this so that
            % further checks in the code below don't fail:
            gpu = reshape(gpu, 1, d2, d1);
        case 3,
            % 3D matrix, hopefully a width x height x channels
            % matrix with x = 1 to 4 channels.
            d1 = size(gpu, 3);
            d2 = size(gpu, 2);
            cc = size(gpu, 1);
            
            % Reject any zero-channel textures or more than 4
            % channel RGBA textures for now:
            if cc < 1 || cc > 4
                error('Provided 3D input GPU matrix "inObj" has less than 1 or more than 4 elements in 1st dimension, which would result in an unsupported color channel count of < 1 or > 4!');
            end
            
            % 3 channel input? This would translate into a 3
            % channel RGB texture, but at least CUDA-5.0 does not
            % support this. Extend it into a 4 channel format:
            if cc == 3
                % We add a value 1.0 to the fourth channel,
                % resulting in a alpha channel of 1 for fully
                % opaque:
                oldgpu = gpu;
                gpu = zeros(4, d2, d1, GPUsingle);
                gpu(1:3, :, :) = oldgpu(1:3, :, :);
                gpu(4  , :, :) = 1.0;
                cc = 4;
                
                % Give a performance warning:
                warning('GPUTYPEFROMTOGL:gputypeRGBtoRGBAautocast', 'Input GPU 3D matrix inObj has 3 elements in 1st dimension, which would result in a unsupported RGB texture! Extending to RGBA texture with A=1.');
            end
            
        otherwise,
            error('Input argument N-d matrix inObj has more than 3 dimensions. This is unsupported for conversion to OpenGL objects!');
    end
    
    % Hopefully 'gpu' is now safe to convert into a OpenGL object.
end

% GL object is a Psychtoolbox texture/offscreen window handle?
if glObjType == 0
    % Yes.
    if direction == 0
        % OpenGL -> GPUtype conversion:
        texid = inObj;
        
        if isempty(texid)
            error('No valid Psychtoolbox OpenGL input object provided!');
        end
    else
        % GPUtype -> OpenGL conversion:
        texid = outObj;
        if isempty(texid)
            % No existing texture object provided as output destination.
            % Create a 32 bpc float texture of matching format.
            for win = Screen('Windows')
                if Screen('WindowKind', win) == 1
                    break
                end
            end
            
            if isempty(win) || Screen('WindowKind', win) ~= 1
                error('No onscreen window opened. This does not work without at least one open onscreen window.');
            end
            
            % Create a 32 bpc float texture 'float = 2', with no need for
            % orientation swap (transpose) 'textureOrientation = 3). We
            % assume all buffers derived from the GPU backend are always in
            % upright row-major format, like Offscreen windows. Input code
            % therefore must do needed conversions.            
            texid = Screen('MakeTexture', win, zeros(d1, d2, cc), [], [], 2, 3);
        end
    end
    
    % Get OpenGL texture handle and texture target of underlying OpenGL
    % texture for given Psychtoolbox object handle:
    [gltexid, gltextarget] = Screen('GetOpenGLTexture', texid, texid);
    
    % Retrieve width and height of texture in texels:
    [width, height] = Screen('Windowsize', texid);

    glPushAttrib(GL.TEXTURE_BIT);
    glBindTexture(gltextarget, gltexid);
    
    % Query bits per pixel:
    bpc = glGetTexLevelParameteriv(gltextarget, 0, GL.TEXTURE_RED_SIZE);
    if bpc == 0
      bpc = glGetTexLevelParameteriv(gltextarget, 0, GL.TEXTURE_LUMINANCE_SIZE);
    end
    
    bpp = bpc;
    bpp = bpp + glGetTexLevelParameteriv(gltextarget, 0, GL.TEXTURE_GREEN_SIZE);
    bpp = bpp + glGetTexLevelParameteriv(gltextarget, 0, GL.TEXTURE_BLUE_SIZE);
    bpp = bpp + glGetTexLevelParameteriv(gltextarget, 0, GL.TEXTURE_ALPHA_SIZE);
     
    % Number of channels == Bits per pixel bpp / Bits per component, e.g., RED channel:
    nrchannels = bpp / bpc;
    
    % Translate to bytes per pixel:
    bpp = bpp / 8;
        
    glBindTexture(gltextarget, 0);
    glPopAttrib();
    
    % fprintf('nrchannels = %i  : Byteperpixel = %i  : bpc = %i\n', nrchannels, bpp, bpc);
    
    % Override a detected RGB32F format to become a 4-channel RGBA32F format. Why?
    % Because at least NVidia on Linux silently allocates storage for RGBA32F when asked
    % for RGB32F, ie., it allocates essentially a RGBX32F format with padding. The problem
    % is that the system lies about this and reports internal format as RGB32F and bpp
    % bits per pixel as 96 bpp instead of the real 128 bpp. This would cause us to misallocate
    % memory and copy the wrong amount, leading to incomplete damaged data transfers. We try
    % to work around this special case by faking the real format and just hope for the best...
    if (bpc == 32 && nrchannels == 3)
      nrchannels = 4;
      bpp = 16;
    end    
end

% GL object is a struct with OpenGL object handle, target, and other info?
if glObjType == 1
    if direction == 0
        % OpenGL -> GPUtype conversion:        
        texstruct = inObj;
        if isempty(texstruct)
            error('No valid Psychtoolbox OpenGL input object provided!');
        end
    else
        % GPUtype -> OpenGL conversion:
        texid = outObj;
        if isempty(texid)
            % No existing texture object provided as output destination.
            error('Creating an OpenGL object from a given GPU object type is not yet supported.');
        end
    end

    if ~isstruct(texstruct)
        error('No OpenGL info struct for inObj provided! Must be a struct!');
    end
    
    % Extract info:
    try
        gltexid = texstruct.glhandle;
        gltextarget = texstruct.gltarget;
        width = texstruct.width;
        height = texstruct.height;
        bpp = texstruct.bpp;
        nrchannels = texstruct.nrchannels;
    catch %#ok<CTCH>
        error('OpenGL info struct inObj is malformed or misses fields!');
    end
end

% Use currently bound drawBufferFBO if imaging pipeline is active --
% accessing the regular onscreen windows virtual framebuffer. inObj is a
% onscreen window handle:
if glObjType == 2
    if direction == 0
        % OpenGL -> GPUtype conversion:        
        win = inObj;
        if isempty(win)
            error('No valid Psychtoolbox onscreen window provided!');
        end
    else
        % GPUtype -> OpenGL conversion:
        win = outObj;
        if isempty(win)
            % No existing onscreen window provided as output destination.
            error('Creating a virtual framebuffer from a given GPU object type is not supported.');
        end
    end

    % Make sure inObj is a onscreen window, with imaging pipeline active
    % and in proper format:
    if Screen('WindowKind', win) ~= 1
        error('For glObjType 2, inObj must be a valid onscreen window handle. This is something else!');
    end
    
    % This queries window properties and binds the FBO for the onscreen
    % windows virtual framebuffer if it isn't already bound:
    winfo = Screen('GetWindowInfo', win);
    if ~bitand(winfo.ImagingMode, kPsychNeedFastBackingStore) || winfo.BitsPerColorComponent < 32
        error('For glObjType 2, onscreen window must have imaging pipeline enabled with a 32 bpc float framebuffer!');
    end
    
    % Proper FBO is bound. Query its color attachment zero, which is the
    % OpenGL handle of the attached texture or renderbuffer:
    gltexid = glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT);
    
    % Query type of attachment:
    gltextarget = glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE_EXT);
    
    % Texture?
    if gltextarget == GL.TEXTURE
        % Yes: We only support rectangle textures in the imaging pipeline,
        % so this is our final target:
        gltextarget = GL.TEXTURE_RECTANGLE_EXT;
    else
        % No: A renderbuffer:
        gltextarget = GL.RENDERBUFFER;
    end

    % Only 4 channel RGBA32F supported, aka 16 Bytes per pixel:
    bpp = 16;
    nrchannels = 4;
    [width, height] = Screen('Windowsize', win);
end

% Use currently bound OpenGL FBO, assuming the imaging pipeline is active
% and properly setup -- otherwise we'd crash or screw up.
if glObjType == 3
    % Proper FBO is hopefully bound. Query its color attachment zero, which
    % is the OpenGL handle of the attached texture or renderbuffer:
    gltexid = glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT);
    
    % Query type of attachment:
    gltextarget = glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE_EXT);
    if gltextarget == GL.FRAMEBUFFER_DEFAULT
        error('For glObjType 3, an OpenGL FBO must be bound, not the system default framebuffer, as here!');
    end
    
    % Query bits per pixel:
    bpc = glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_RED_SIZE);
    bpp = 0   + bpc;
    bpp = bpp + glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_GREEN_SIZE);
    bpp = bpp + glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_BLUE_SIZE);
    bpp = bpp + glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_ALPHA_SIZE);
     
    % Number of channels == Bits per pixel bpp / Bits per component, e.g.,
    % RED channel:
    nrchannels = bpp / bpc;
    
    % Translate to bytes per pixel:
    bpp = bpp / 8;

    % Texture?
    if gltextarget == GL.TEXTURE
        % Yes: We only support rectangle textures in the imaging pipeline,
        % so this is our final target:
        % TODO FIXME MK: Technically not quite correct, as at least
        % Screen('TransformTexture') could also use a GL_TEXTURE_2D target
        % instead of texture rectangle. However, this is a seldomly used
        % special case and i don't know at the moment how to find out which
        % target is actually used.
        gltextarget = GL.TEXTURE_RECTANGLE_EXT;
        
        % Query size width x height of texture image:
        glPushAttrib(GL.TEXTURE_BIT);
        glBindTexture(gltextarget, gltexid);
        width  = glGetTexLevelParameteriv(gltextarget, 0, GL.TEXTURE_WIDTH);
        height = glGetTexLevelParameteriv(gltextarget, 0, GL.TEXTURE_HEIGHT);
        glBindTexture(gltextarget, 0);
        glPopAttrib();
    else
        % No: A renderbuffer:
        gltextarget = GL.RENDERBUFFER;
        glBindRenderbuffer(gltextarget, gltexid);
        width  = glGetRenderbufferParameteriv(gltextarget, GL.RENDERBUFFER_WIDTH);
        height = glGetRenderbufferParameteriv(gltextarget, GL.RENDERBUFFER_HEIGHT);
        glBindRenderbuffer(gltextarget, 0);
    end
end

% Unmap or Unregister object from cache?
if cmd == 2 || cmd == 3
    memcpyCudaOpenGL(cmd - 1, gltexid, gltextarget);
    return;
end

% Map OpenGL resource, then return a memory pointer in a uint64 for it?
if cmd == 5 || cmd == 6
    % This maps the resource and returns a pointer to it in uint64 outObj:
    % cmd 5 and 6 are translated into direction values 0 and 1 via 'cmd - 5'. This
    % is important to get the correct mapping flags for resource mapping (readonly vs.
    % writeonly):
    outObj = memcpyCudaOpenGL(4, gltexid, gltextarget, 0, 0, cmd - 5, 1, mapflags);
    return;
end

if (nrchannels ~= 1) && (nrchannels ~= 2) && (nrchannels ~= 4)
    error('Tried to convert a 3 layer RGB texture or framebuffer. This is not supported.');
end

% Number of bytes to copy = w * h * bpp:
nrbytes = width * height * bpp;
if nrbytes < 1
    error('Tried to convert an empty texture. Forbidden!');
end

% Is an already existing 'gpu' variable provided for "refill" ?
if ~isempty(gpu)
    % Yes: Check for matching format. If it doesn't match, delete it, so it
    % can be recreated with matching format:
    if ~isa(gpu, 'GPUsingle') || (size(gpu, 1) ~= nrchannels) || (size(gpu, 2) ~= width) || (size(gpu, 3) ~= height) || (0 == getPtr(gpu))
        % No match, or not allocated. Destroy:
        clear gpu;
        % And create empty for recreation below:
        gpu = [];
        
        % This must not happen in GPU -> OpenGL mode:
        if direction == 1
            error('Incompatible GPUtype inObj variable provided for update of OpenGL object. How is this supposed to work?!?');
        end
    end
    
    % Is gpu a complex matrix - which we can't handle?
    if ~isreal(gpu)
        % Yes: Only extract and convert real part:
        gpu = real(gpu);
        
        % Give data-loss / performance warning:
        warning('GPUTYPEFROMTOGL:gputypeComplexToRealcast', 'Input GPU matrix inObj stores complex numbers, which we cannot store! Throwing away the imaginary component of each complex number!');
    end
end

% Need to create a new gpu variable?
if isempty(gpu)
    % Yes: Create a new GPUsingle GPU variable:
    gpu = GPUsingle;
    % Set it to real format:
    setReal(gpu);
    
    % Set its size: We *must* double-cast the size vector here, because the
    % gpuType == 3 path delivers int32's and GPUmat doesn't like this at
    % all, punishing us with GPUallocVector failure, if we don't cast to
    % double():
    setSize(gpu, double([nrchannels, width, height]));
    
    % Allocate its CUDA backing memory:
    GPUallocVector(gpu);
end

% Retrieve CUDA memory pointer to it:
gpuptr = getPtr(gpu);
if gpuptr == 0
    error('Memory allocation on GPU failed!');
end

% Perform copy of image content from OpenGL texture into CUDA backing store:
memcpyCudaOpenGL(3, gltexid, gltextarget, gpuptr, nrbytes, direction, keepmapped, mapflags);

if direction == 0
    % OpenGL -> GPU:
    outObj = gpu;
else
    % GPU -> OpenGL:
    if glObjType == 0
        outObj = texid;
    else    
        outObj = outObj; %#ok<ASGSL>
    end
end

return;

end
