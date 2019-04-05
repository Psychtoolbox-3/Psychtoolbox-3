function Add2DConvolutionToGLOperator(gloperator, kernel, opname, nrinputchannels, nroutchannels, debug, shadertype)
% Add2DConvolutionToGLOperator(gloperator, kernel [, opname] [, nrinputchannels] [, nroutchannels] [, debug] [, shadertype])
% 
% Add a non-separable 2D convolution kernel to GL operator 'gloperator'.
%
% 'kernel' is a simple m-by-n matrix of floating point numbers with m and
% n being odd numbers, e.g., 1x1, 3x3, 5x5, 7x7, 9x9,..., 1x3, 1x9, 7x1 ...
% Each entry in the kernel matrix is used as a weight factor for the
% convolution.
%
% The simplest way to get a 2D kernel is to use the function
% kernel = fspecial(...); fspecial is part of the Matlab image
% processing toolbox, see "help fspecial" for more information.
%
% 'opname' is the optional name of the convolution operation, meant as
% debugging aid.
%
% 'nrinputchannels' = The number of image channels to use as input for the
% convolution. Possible values are: 3 = Red, Green and Blue color channels are
% provided as part of a true-color image, don't use the alpha channel (if
% any) for convolution but just pass it through unmodified. 1 = The image
% only defines a luminance channel for convolution, an (optional) alpha
% channel is passed through unmodified. 4 = Use all four channels (Red,
% green, blue, alpha) for convolution.
%
% 'nroutchannels' = The number of channels to convolve as output: 3 =
% Convolve each of the three color channels red, green and blue separately
% by the same kernel. An (optional) alpha channel is passed through unmodified.
% 1 = Output a filtered luminance channel, and an (optional) unmodified
% alpha channel. If input is a 3 channel RGB image, then the RGB image will
% get converted to luminance before convolution. 4 = Filter all four
% channels independently.
%
% Typical settings:
%
% Filter a RGB(A) image: nrinputchannels = 3, nroutchannels = 3.
% Filter a RGB(A) image into grayscale: nrinputchannels = 3, nroutchannels = 1.
% Filter a luminance(A) image: inputchannels = 1, filteredoutchannels = 1.
% Generic filtering of 4-channel data: nrinputchannels = 4, nroutchannels = 4.
%
% 'debug' Optional debug flag: If set to non-zero, will output some debug
% info about the shader.
%
% 'shadertype' (Optional) The type of internal implementation to choose for
% the operator. This parameter is best left alone, unless you really know
% what you are doing.

% History:
% 08/11/07 Written (MK).

% Check input arguments for validity and assign defaults:

if nargin < 1 || isempty(gloperator) || Screen('WindowKind', gloperator)~=4
    error('You must provide the handle of a valid GL imaging operator ''gloperator''!');
end

if nargin < 2 || isempty(kernel) || ndims(kernel)~=2
    error('You must provide a 2D ''kernel'' matrix as 2nd argument!');
end

% Switch to gloperators OpenGL context:
Screen('GetWindowInfo', gloperator);

% Transpose kernel to make sure we follow the strict definition of
% convolution and stay compatible with Matlab/Octave:
kernel = transpose(kernel);

% Query size of kernel:
kernelw = size(kernel,1);
kernelh = size(kernel,2);

% We only want odd sized kernels of at least 1x1, e.g., 3x3, 5x5, 7x7, ...
if kernelw < 1 || kernelh < 1 || mod(kernelw,2)~=1 || mod(kernelh,2)~=1
    error('Your kernel must be odd-sized and at least of size 1x1!');
end

if nargin < 4 || isempty(nrinputchannels)
    % Default to RGB channels as separate input channels:
    nrinputchannels = 3;
end

if nargin < 5 || isempty(nroutchannels)
    % Default to output of filtered RGB channels:
    nroutchannels = 3;
end

if nargin < 6 || isempty(debug);
    % No debug mode provided: Disable debug output by default:
    debug = 0;
end

if nargin < 7
    shadertype = [];
end

if nargin < 3 || isempty(opname)
    opname = sprintf('2D convolution kernel: w x h = %i by %i, nrin=%i, nrout=%i.', kernelw, kernelh, nrinputchannels, nroutchannels);
end

% Input parsing done. Call helper routine for creation of shader and
% optional lookup texture:
[shader configstring] = EXPCreateStatic2DConvolutionShader(kernel, nrinputchannels, nroutchannels, debug, shadertype);

% Count number of slots in proxy:
count = CountSlotsInGLOperator(gloperator);

% If its non-zero, then some processing will already take place before our
% shader, so we need to perform a FBO ping-pong step:
if count > 0
    Screen('HookFunction', gloperator, 'AppendBuiltin', 'UserDefinedBlit', 'Builtin:FlipFBOs', '');

    if debug > 3
        fprintf('Add2DConvolutionToGLOperator: Switching operator to dual-pass/multi-pass mode.\n');
    end

    if count == 1
        % Count was 1, so its now two: Change operator to be at least dual-pass capable:
        Screen('HookFunction', gloperator, 'ImagingMode', mor(kPsychNeedDualPass, Screen('HookFunction', gloperator, 'ImagingMode')));
    else
        if count > 1
            % Change operator to be multi-pass capable:
            Screen('HookFunction', gloperator, 'ImagingMode', mor(kPsychNeedMultiPass, Screen('HookFunction', gloperator, 'ImagingMode')));
        end
    end
end

% Add shader to user defined blit chain of the proxy:
Screen('HookFunction', gloperator, 'AppendShader', 'UserDefinedBlit', opname, shader, configstring);

% Check if the kernel contains any negative values or if the sum of kernel
% weights is greater than 1.0.
if ~isempty(find(kernel < 0)) || sum(sum(kernel))>1
    % Yes. Convolution result could be a negative number, or a value
    % greater than the 8bpc fixed representable value 1.0, so potential
    % bounce-buffers and output texture must be capable of representing
    % negative numbers and numbers greater than 1 without clamping.
    % --> To accomodate this we need a floating point precision format.
    
    % Check if gloperator is already configured for float precision. If
    % not, then assign 32bpc float request:
    if bitand(Screen('HookFunction', gloperator, 'ImagingMode'), mor(kPsychNeed16BPCFloat, kPsychNeed32BPCFloat)) == 0
        % Not yet set. Choose highest precision:
        Screen('HookFunction', gloperator, 'ImagingMode', mor(kPsychNeed32BPCFloat, Screen('HookFunction', gloperator, 'ImagingMode')));        
        if debug > 3
            fprintf('Add2DConvolutionToGLOperator: kernel contains negative weights or summed weights greater than 1.0 -> Increased precision of operator to 32bpc float.\n');
        end
    end
end

% Done.
return;
