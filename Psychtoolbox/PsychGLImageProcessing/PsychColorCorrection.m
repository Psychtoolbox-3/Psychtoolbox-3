function varargout = PsychColorCorrection(cmd, varargin)
% varargout = PsychColorCorrection(cmd, varargin)
%
% This function is used to setup and control the built-in color correction 
% and transformation mechanisms of PTB's imaging pipeline. It allows to choose
% among different methods of color correction (for calibrated output of
% stimuli) and to change parameters of the chosen method. This only works
% when the imaging pipeline is enabled and with certain output devices.
% 
% This functions are mostly meant to provide fast color correction, e.g.,
% gamma correction, when PTB is used with special high precision / HDR
% output devices, e.g., CRS Bits++ in Mono++ or Color++ mode, or video
% attenuators. When operating such devices, the standard gamma correction
% mechanisms of your graphics card can't be used, so PTB must do it itself.
%
% With a standard 8bpc framebuffer, you won't need this function for simple
% gamma correction or color correction.
% Instead you'd use the Screen('LoadNormalizedGammaTable') command to
% perform gamma correction with the graphics cards built-in gamma tables.
% Same applies to Bits++ box in Bits++ mode, using the built-in tables of
% the Bits++ device, also controlled via Screen('LoadNormalizedGammaTable').
%
% For more complex correction schemes than a lookup table transform you can
% still use this function, even with a normal 8 bit framebuffer, as long as
% the imaging pipeline is enabled. It's just important to note that you
% don't need to use it for a standard framebuffer in the simple case -- the
% gfx cards gamma tables are more efficient for the simple case.
%
% How to use:
%
% 1. Before opening an onscreen window, you use the following
% PsychImaging() setup call to specify the method of display color correction
% to apply, and the view channel to apply it to (in case there are multiple
% like in many stereo display setups):
%
% PsychImaging('AddTask', whichChannel, 'DisplayColorCorrection', methodname);
% - where whichChannel can be 'LeftView' or 'RightView' for the left- or
% right display output channel/device of a stereo setup, or
% 'FinalFormatting' if you have a single display monoscopic setup or want
% to apply the same color correction to both channels of a stereo setup.
% The parameter 'methodname' is a name string with the name of one of the
% supported methods - See overview below for supported methods.
%
% CAUTION: The order of specifications matters! If you use multiple color
% correction methods or other image processing operations simultaneously,
% make sure to specify them in the order in which they should be executed.
% The system tries to order operations in a reasonable way, but it is not
% fool-proof!
%
% 2. Then, after you've specified all other window parameters via the
% PsychImaging() subcommands, you open the onscreen window via the usual
% win = PsychImaging('OpenWindow', ....); call, as always when you use the
% imaging pipeline. This will open the window and apply the chosen color
% correction method, choosing reasonable default parameters for the method
% at hand.
%
% 3. At any time in your script you can change the operating parameters of
% the chosen color correction method via the PsychColorCorrection()
% subfunctions mentioned below. You'll have to specify the window handle
% for the onscreen window whose parameters should be changed, and - in a
% stereo display setup with separate color correction parameters for each
% of the display channels - the name 'whichChannel' of the display channel
% to change, e.g., 'LeftView'. Changes take effect at next Screen('Flip').
%
%
%
% Supported display color correction methods:
% ===========================================
%
% These are the method names that can be passed as 'methodname' parameter to
% PsychImaging('DisplayColorCorrection', ....., methodname):
%
% 'methodname' is the name string of one of the supported methods:
%
% * 'None': Don't do anything. Fastest, but not safest. 
%
%   If your stimulus contains (by accident) color or luminance values outside
%   the displayable range, the corresponding pixels may show undefined output
%   color -- Likely not what you want. This mode is about 10% faster than
%   'ClampOnly'. This mode is the default if no method is selected, unless
%   you use Bits++ in Mono++ or Color++ mode, where 'ClampOnly' is the
%   default.
%
%
% * 'ClampOnly': Do not apply any color transformation. Default for Bits++
%   in Mono++ or Color++ mode:
%
%   In this mode, values are clamped against the limits set via
%   'SetColorClampingRange' (see below) to always keep them in the
%   requested range set via PsychColorCorrection('SetColorClampingRange'),
%   nothing else.
%
%
% * 'CheckOnly': Do not apply any color transformation.
%
%   In this mode, values are checked against the limits set via
%   PsychColorCorrection('SetColorClampingRange'). Values outside that
%   range are visually coded, so you should be able to see troublesome
%   areas via visual inspection: All values are first clamped, then
%   inverted, in the hope that this creates clearly detectable artifacts in
%   your stimulus. We may change the "marker" method in the future if we
%   happen to find a better visual marker. This mode is meant for debugging
%   your stimulation scripts, not for running them on your subjects!
%
%
% * 'SimpleGamma' : Apply some power-law based gamma correction:
%
%   Simple gamma correction means: Apply a power-law to incoming
%   values: outcolor = incolor ^ EncodingGamma. incolor is the uncorrected
%   input color or intensity, outcolor the corrected one, EncodingGamma is
%   encoding gamma value to apply. See PsychColorCorrection('SetEncodingGamma')
%   for how to set the gamma values. After gamma correction, output values
%   are clamped against the range set via PsychColorCorrection('SetColorClampingRange').
%
%   See PsychColorCorrection('SetExtendedGammaParameters') on how to supply
%   additional parameters for use of a more complex gamma correction function.
%
%
% * 'MatrixMultiply4' : Multiply with 4x4 color transformation matrix:
%
%   Extend RGB component of color by a fourth "constant 1" component,
%   multiply by specified 4-by-4 matrix, normalize by dividing result by
%   4'th component, return transformed rgb components as output color. Pass
%   through alpha-channel unmodified.
%
%   rgbw' = M * [r,g,b,1]'
%   out.rgb = [r', g', b'] / w'
%   out.a = a
%
% * 'SensorToPrimary' : Implement the same conversion as SensorToPrimary.m
%
%   This can convert from XYZ sensor tristimulus values to standard RGB
%   primary color values. This is the right thing to do if your display
%   device is linearized via use of a proper gamma correction method.
%
%   Use the function PsychColorCorrection('SetSensorToPrimary') below to
%   assign the 'cal'ibration struct to use for actual conversion from XYZ
%   color values to RGB primaries.
%
%
% * 'xyYToXYZ' : Implement the same conversion as xyYToXYZ.m
%
%   This can convert to tristimulus XYZ color coordinates from
%   chromaticity (x,y) and luminance (Y) color coordinates xyY. It performs
%   the same color space conversion as the M-Function XYZ = xyYToXYZ(xyY).
%
%
% * 'LookupTable' : Apply color correction by color table lookup, ie. a CLUT.
%
%   This will allow to pass in a color lookup table of selectable
%   granularity (ie., number of slots) and range, which is later on used to
%   lookup corresponding corrected color values for given framebuffer input
%   values.
%
%
% * 'LookupTable3D' : Apply color correction by a 3D color table lookup.
%
%   This will allow to pass in a 3D color lookup table of selectable
%   granularity (ie., number of slots per dimension) and range, which is
%   later on used to lookup corresponding corrected color values for given
%   framebuffer input values.
%
%
% * 'GainMatrix' : Apply color gain correction by 2D gain matrix lookup.
%
%   This allows to apply a 2D matrix G which stores luminance- or color gain
%   correction factors for each single output pixel of the display. For
%   each 2D display pixel location (x,y), the stimulus image I(x,y) will be
%   multiplied with the corresponding gain factor G(x,y) of the gain matrix
%   and the result O(x,y) used for further processing and display. G(x,y)
%   can be a single scalar for luminance correction, or - if G is a 3-layer
%   matrix - a RGB vector with individual gains for each color channel and
%   pixel location.
%
%   O(x,y) = I(x,y) * G(x,y).
%
%   If you want to combine this with one of the other correction methods,
%   e.g., gamma correction, you should issue this command first, because
%   the non-linear gamma correction should apply to the output of this
%   method for correct results.
%
%   After you've opened your onscreen 'window', you'll need to define the 2D
%   gain 'matrix' via a call to ...
%
%   PsychColorCorrection('SetGainMatrix', ...);
%
%   See below for description of 'SetGainMatrix'.
%
% * 'AnaglyphStereo' : Apply anaglyph stereo algorithm.
%
%   This loads a similar anaglyph shader as the one used in stereoModes 6
%   to 9 when the function SetAnaglyphStereoParameters('ColorAnaglyphMode')
%   or its siblings was used. This is useful if you want to employ anaglyph
%   stereo presentation, but not by rendering a full anaglyph image into
%   one single framebuffer for output to a single display or projector, but
%   if you want to direct the "left-eye" anaglyph image to a different
%   display or projector than the "right-eye" anaglyph image. An example
%   would be having two video projectors attached to two video outputs of a
%   graphics card. One projector shall project the left-eye image, the
%   other one the right-eye image. In this case you'd choose a stereoMode
%   of 4 or 5 (typically on Linux or Windows) for desktop spanning stereo
%   output, or 10 (on OSX) for dual-window stereo output. Then you'd use ...
%   PsychImaging('AddTask', 'LeftView', 'DisplayColorCorrection', 'AnaglyphStereo')
%   ... and ...
%   PsychImaging('AddTask', 'RightView', 'DisplayColorCorrection', 'AnaglyphStereo')
%   ... to add an anaglyph shader to the end of each view channel.
%
%   After adding individual shaders for each image processing channel and opening
%   the onscreen window(s), you can use SetAnaglyphStereoParameters() to
%   parameterize the anaglyph stereo presentation as if it were created via
%   regular anaglyph stereo setup in stereomoded 6-9.
%
%
% Supported runtime Subfunctions:
% ===============================
%
% oldlevel = PsychColorCorrection('Verbosity' [, newlevel]);
% Return current level of verbosity in optional 'oldlevel', optionally set
% new level of verbosity via 'newlevel'. The level of verbosity affects how
% much status output is printed in some routines: 0 = Nothing, 1 = Only
% errors, 2 = Errors and warnings, 3 = Errors + Warnings + Some info.
%
%
% The following routines must be called *after* opening a window for which color
% correction is enabled. They can be called anytime and changed settings
% will apply at the next Screen('Flip'). If your 'window' is actually a
% stereo display window, you may want or need to provide the optional
% 'viewId' parameter to tell PTB which of the two stereo view channels
% shall be changed in its settings. If both stereo views are displayed on
% the same physical display device, this is not needed. If the separate
% views go to separate physical displays, you may need to calibrate them
% separately. Allowable values for viewId are 'AllViews', 'LeftView' and
% 'RightView', corresponding to the 'whichChannel' setting that you used
% when setting up the window with PsychImaging().
%
%
% PsychColorCorrection('SetColorClampingRange', window, min, max [,viewId]);
% - Set the range of allowable output color or luminance intensity values
% to the interval [min; max] for onscreen window 'window'. Values outside
% that range get either clamped to the 'min'imum or 'max'imum value, or -
% in 'CheckOnly' mode - will be visually marked as out of range. The default
% range is [0.0 ; 1.0] -- The range that your display device can really
% display.
%
%
% PsychColorCorrection('SetEncodingGamma', window, gamma [,viewId]);
% - Set the gamma value to use for gamma correction on window 'window'.
% 'gamma' can be either a single scalar if the same gamma should apply to
% all color channels (or single luminance channel), or it can be a
% three-component [gammaRed, gammaGreen, gammaBlue] vector, if each color
% channel should be gamma corrected with an individual gamma value.
%
% How the value gets applied depends on the chosen method of color
% correction (see above). The simplest method ('SimpleGamma') performs a
% simple power-law mapping of input values to output values: out = in ^ gamma.
% 'in' must be greater than zero, and 'gamma' must be greater than zero,
% otherwise results may be undefined, depending on your graphics hardware.
% However, usually only encoding 'gamma' values in the range of about
% 0.33 - 1.0 are meaningful.
%
% Example: If your monitor has a "decoding gamma" of 1.8, the proper
% setting for 'gamma' would be gamma = 1/1.8. For a decoding gamma of 2.2,
% you'd choose gamma = 1/2.2 ...
%
%
% PsychColorCorrection('SetExtendedGammaParameters', window, minL, maxL, gain, bias [,viewId]);
% - Set the additional (optional) parameters to fine-tune gamma correction on
% window 'window'. All these parameters have reasonable defaults. All
% parameters can be supplied as a scalar value if the same setting shall
% apply to all color channels (or a single luminance channel), or you can
% provide 3-component vectors with one component for each color channel:
%
% After this function has been called at least once, the following formula
% will be used to map input values to output values ['gamma' is as set by
% the 'SetEncodingGamma' function, 'in' is input, 'out' is output value]:
%
% out = bias + gain * ( ((in-minL) / (maxL-minL)) ^ gamma )
%
% Required parameters:
% 'minL' Minimum expected input luminance/intensity value (Default is 0.0).
% 'maxL' Maximum expected input luminance/intensity value (Default is 1.0).
% 'gain' Gain factor to apply after power-law mapping (Default is 1.0).
% 'bias' Bias/Offset to apply to final result before output (Default is 0.0).
%
%
% PsychColorCorrection('SetMultMatrix4', window, matrix [, viewId]);
% - Set the 4-by-4 color transformation matrix to use for the 'MatrixMultiply4'
% color correction method. 'matrix' must be the 2D 4 rows by 4 columns
% matrix to use. By default, the matrix is set to an identity matrix.
%
%
% PsychColorCorrection('SetSensorToPrimary', window, cal [, viewId]);
% - Set the 'cal'ibration struct to use for the 'SensorToPrimary' color
% correction method. 'cal' must be the same input format as used for the
% M-Function SensorToPrimary(). By default, the transformation is a "no
% operation".
%
%
% PsychColorCorrection('SetLookupTable', window, clut [, viewId][, maxinput=1][, scalefactor]);
% - Assign color lookup table 'clut' for use with color correction method
% 'LookupTable'. 'clut' must be a 1 column vector for pure luminance lookup
% tables, or a 3 column matrix for RGB color lookup tables with one column
% per color channel, ie., [1,2,3] = [Red, Green, Blue]. clut must have at
% least 1 row, but usually will have way more than 2 rows, typically almost
% as many rows as n = 2^bpc for a given output device bitdepths bpc. For a
% 10 bit output device, n would be usually 2^10 = 1024 rows for a perfect
% one-to-one mapping. At runtime, color correction will be performed by the
% following formula: Be Rin, Gin, Bin the input red, green and blue color
% components, and Rout, Gout, Bout the final output value for the
% framebuffer. First Rin, Gin and Bin are clamped individually to the range
% 0.0 - 'maxinput' (maxinput is 1.0 by default), scalefactor is chosen by
% default as scalefactor = [number of rows in clut - 1] / maxinput, ie., it
% maps the possible input range 0 - maxinput to the full range of row
% indices 1 - rowcount to cover the full range of entries stored in the
% clut. This is the most reasonable default, but can be changed by the
% optional 'scalefactor' and 'maxinput' arguments.
%
% Then the output color for each component is looked up in the proper slot
% (= row index) of the passed clut:
%
% Rout = clut(Rin * scalefactor,1);
% Gout = clut(Gin * scalefactor,2);
% Bout = clut(Bin * scalefactor,3);
%
% Color values for fractional indices inbetween reference values in the
% clut are interpolated linearly between the two nearest neighbour
% reference values --> linear interpolation.
%
% Finally, Rout, Gout and Bout are clamped to the valid output range as set
% by the function PsychColorCorrection('SetColorClampingRange', ...); by
% default to the range 0.0 - 1.0.
%
%
% PsychColorCorrection('SetLookupTable3D', window, clut [, viewId][, maxinput=1][, scalefactor][, precision=0][, interpolate=1]);
% - Assign 4D color lookup table 'clut' for use with color correction method
% 'LookupTable3D'. 'clut' must be a 4D 3-by-m-by-n-by-p matrix. The first dimension encodes the
% output color values to use:
% clut(1,r,g,b) == Corrected red color value for input color [r,g,b].
% clut(2,r,g,b) == Corrected green color value for input color [r,g,b].
% clut(3,r,g,b) == Corrected blue color value for input color [r,g,b].
%
% clut must have at least one element in each color index dimension, ie., m, n
% and p must be >= 1, but usually will have more elements in each dimension
% for a meaningful lookup color correction. In theory you would need m, n
% and p to be == 2^bpc for a given output device bitdepths bpc, e.g, for a
% 8 bit output device, m,n,p would need to be 2^8 = 256 elements for a perfect
% one-to-one mapping. In reality you likely don't want to use such large
% sizes, as such a huge and dense 3D CLUT would take up considerable
% amounts of graphics memory and cause large slowdowns of all drawing
% operations. For good performance and portability of your code to older
% graphics cards, choose modest sizes, as small as possible.
%
% At runtime, color correction will be performed by the following 3D table
% lookup procedure: Let Rin, Gin, Bin be the input red, green and blue
% color components, and Rout, Gout, Bout the final output values for the
% framebuffer. First Rin, Gin and Bin are clamped individually to the range
% 0.0 - 'maxinput' (maxinput is 1.0 by default), scalefactor is chosen by
% default as scalefactor = (1.0 / maxinput), ie., it maps the possible
% input range 0 - 'maxinput' to the range 0.0 - 1.0, which covers the full
% range of entries stored in the clut. This is the most reasonable default,
% but can be changed by the optional 'scalefactor' and 'maxinput'
% arguments.
%
% Then the output color for each component is looked up in the proper 3D
% slot of the passed clut:
%
% Rout = clut(1, Rin * scalefactor, Gin * scalefactor, Bin * scalefactor);
% Gout = clut(2, Rin * scalefactor, Gin * scalefactor, Bin * scalefactor);
% Bout = clut(3, Rin * scalefactor, Gin * scalefactor, Bin * scalefactor);
%
% Color values (Rin, Gin, Bin) == (0,0,0) map to the first elements in the
% cluts dimensions, ie., clut(:,1,1,1). Maximum values (maxinput, maxinput,
% maxinput) map - after scaling with the default scalefactor - to
% coordinates (1.0, 1.0, 1.0) in the normalized 3D color coordinate space
% and are looked up in the maximal clut element indices of each dimension,
% ie., clut(:,m,n,p) for our 3-by-m-by-n-by-p clut. Intermediate values are
% mapped accordingly to proper clut element indices.
%
% By default, color values for fractional indices inbetween reference
% values in the clut are interpolated linearly between the eight nearest
% neighbour reference values in the 3 dimensional space --> This is
% trilinear interpolation across all 3 color dimensions of the CLUT. If you
% set the optional parameter 'interpolate' to zero, then simple nearest
% neighbour sampling is performed instead.
%
% The optional 'precision' flag controls the precision with which entries
% in the CLUT should be stored and processed: precision = 0 is the default
% and stores values with 8 bit precision for 256 different intensity levels
% for Rout, Gout and Bout. A setting of 1 will store values with 16 bpc
% floating point precision to resolve up to 10 bits or 1024 levels of
% linear precision. A setting of 2 will store values with 32 bpc floating
% point precision for up to 23 bits of linear precision. Be aware that
% precision values > 0 will increase memory consumption by a factor of 2x
% or 4x, which can be significant for lookup tables of non-trivial size.
%
% Final looked up, Rout, Gout and Bout are clamped to the valid output range as set
% by the function PsychColorCorrection('SetColorClampingRange', ...); by
% default to the range 0.0 - 1.0.
%
% Essentially, this color correction allows to define arbitrary mappings of
% RGB input triplets to RGB output triplets, providing a large amount of
% flexibility. Be aware though that strongly discontinuous mappings of input
% colors to output colors can have a significant negative impact on the
% drawing performance of your graphics card. Define your clut's wisely!
%
%
% PsychColorCorrection('SetGainMatrix', window, matrix [, viewId][, precision=2]);
%
% - Set gain matrix for method 'GainMatrix'.
% If matrix is a 2D matrix, the gain will be applied to all color
% channels equally. If matrix is a 3D matrix, matrix(y,x,1) will define
% the red channel gain, matrix(y,x,2) will define the green channel gain, 
% and matrix(y,x,3) will define the blue channel gain.
%
% The optional 'precision' parameter defines the numerical precision with
% which the gain factors are stored. The default setting of 2 stores with
% 32 bit floating point precision - about 6 digits behind the decimal
% point. A Setting of 1 stores with 16 bit float precision, about 3 digits.
% A Setting of 0 stores with 256 levels, about 2 digits. A lower precision
% is less precise but allows for faster processing and higher redraw rates
% if needed.
%
%
%
% Internal commands, usually not meant for direct use by pure mortals:
% ====================================================================
%
% All these methods are usually called from within PsychImaging() to do the
% dirty setup work...
%
% Call this *before* opening a window:
%
%
% PsychColorCorrection('ChooseColorCorrection', methodname);
% - Specify the method to be used for color correction for the next
% onscreen window that will be opened. This needs to be called *before* the
% window is opened, however its usually done automatically at the right
% moment by routines like PsychImaging() or BitsPlusPlus() if you use these
% to open windows.
%
%
% Called after Screen('OpenWindow') during shader and pipeline setup:
%
% [shader, idstring, configString, overrideMain] = PsychColorCorrection('GetCompiledShaders', window, debuglevel);
% - Compile corresponding shaders for chosen color correction method,
% return shaderhandles and idstring to calling routine. That routine will
% link the returned shaders with other shader code to produce the final
% GLSL program object for color conversion and output formatting. 'shader'
% is the GLSL shader handle, idstring the name string for the shader,
% configString the shader option string for the 'Hookfunction' call, e.g.,
% to bind additional LUT textures, etc.
%
%
% Called after linking and attaching the final processing GLSL program
% objects and slots to the imaging pipelines hook chain(s):
%
%
% PsychColorCorrection('ApplyPostGLSLLinkSetup', window, viewId);
% - Perform whatever setup work is needed after final GLSL program object
% has been created and attached to imaging pipeline.
%

% History:
% 16.04.2008 Written (MK).
% 18.05.2008 Revised, improved help text, fine-tuning etc. (MK)
% 04.07.2009 Add CLUT based color correction. (MK)
% 10.10.2009 Add 'SetExtendedGammaParameters' for extended gamma correction. (MK)
% 05.03.2010 Add 'GainMatrix' and 'SetGainMatrix' for display
%                 shading/vignetting correction. (MK)
% 14.08.2012 Add 'AnaglyphStereo' to apply anaglyph stereo mode onto
%            separate display outputs, e.g., as on MPI Kuka projection
%            setup. (MK)
% 21.08.2012 Add support for 3D CLUTs for "color cube" lookups. (MK)
% 21.08.2012 Add 'Verbosity' level support. (MK)
% 24.08.2012 Add 'SetMultMatrix4' support for 4-by-4 matrix multiply. (MK)
% 25.01.2013 Add 'SensorToPrimary' XYZ -> RGB colorspace conversion. (MK)
% 26.01.2013 Add 'xyYToXYZ' xyY -> XYZ colorspace conversion. (MK)

% GL is needed for shader setup and parameter changes:
global GL;
persistent verbosity;
persistent specReady;
persistent icmSpec;
persistent icmDataForHandle;

if isempty(specReady)
    % We default to a setting of 'ClampOnly' if no spec specified. Usercode can
    % override this...
    specReady = 1;
    icmSpec.type = 'ClampOnly';
    verbosity = 3;
end

% Child protection:
if nargin < 1
    error('Subcommand specification missing!');
end
        
if isempty(cmd)
    error('Subcommand specification missing!');
end
       
% Subcommand dispatch:

if strcmpi(cmd, 'ChooseColorCorrection')

    if nargin < 2
        error('You must specify a color-correction method!');
    end
    
    % Assign method, that's it...
    icmSpec.type = varargin{1};
    specReady = 1;
    
    return;
end

if strcmpi(cmd, 'Verbosity')
    varargout{1} = verbosity;
    
    if nargin == 2
        verbosity = varargin{1};
    end
    
    return;
end

% Retrieve a (previously specified) collection of precompiled shader
% objects. These implement the subroutines for color processing, according
% to given spec.
if strcmpi(cmd, 'GetCompiledShaders')
    if ~specReady
        error('"GetCompiledShaders" called, but specification of what to compile is unavailable!');
    end
    
    if nargin >= 3
        debuglevel = varargin{2};
    else
        debuglevel = 0;
    end
    
    % Need GL from here on...
    if isempty(GL)
        InitializeMatlabOpenGL([], [], 1);
    end

    % No config string by default:
    icmConfig = '';
    
    % No override main routine by default:
    icmOverrideMain = '';

    switch(icmSpec.type)
        case {'None'}
            % Load pass-through shader:
            icmShaders = LoadShaderFromFile('ICMPassThroughShader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);

        case {'ClampedNoName'}
            % Load clamped pass-through shader:
            icmShaders = LoadShaderFromFile('ICMClampedPassThroughShader.frag.txt', [], debuglevel);
            % but set its id string to some special string that is hidden from standard search:
            icmIdString = 'CMIH:ClampedNoName';

        case {'ClampOnly'}
            % Load clamped pass-through shader:
            icmShaders = LoadShaderFromFile('ICMClampedPassThroughShader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);

        case {'CheckOnly'}
            % Load checked pass-through shader:
            icmShaders = LoadShaderFromFile('ICMCheckedPassThroughShader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);

        % Simple power-law bog-standard gamma correction.
        case {'SimpleGamma'}
            % Load our bog-standard power-law shader:
            icmShaders = LoadShaderFromFile('ICMSimpleGammaCorrectionShader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);
            
        case {'MatrixMultiply4', 'SensorToPrimary'}
            % Load our matrix multiply shader:
            icmShaders = LoadShaderFromFile('ICMMatrixMult4Shader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);
            
        case {'xyYToXYZ'}
            % Load our xyYToXYZ colorspace conversion shader:
            icmShaders = LoadShaderFromFile('ICMConvert_xyYToXYZ_Shader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);
            
        % Color correction by CLUT texture lookup table operation:
        case {'LookupTable'}
            % Load our 1D CLUT color correction shader:
            icmShaders = LoadShaderFromFile('ICMCLUTCorrectionShader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);

            % Generate texture handle for fillout later on:
            icmSpec.icmlutid = glGenTextures(1);
            
            % Build config string to bind and use our CLUT texture:
            icmConfig = sprintf('TEXTURERECT2D(2)=%i', icmSpec.icmlutid);
            
        % Color correction by indexing into 3D texture lookup table:
        case {'LookupTable3D'}
            % Load our 3D CLUT color correction shader:
            icmShaders = LoadShaderFromFile('ICM3DCLUTCorrectionShader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);

            % Generate texture handle for fillout later on:
            icmSpec.icmlutid = glGenTextures(1);
            
            % Build config string to bind and use our CLUT 3D texture:
            icmConfig = sprintf('TEXTURE3D(2)=%i', icmSpec.icmlutid);

        % Vignetting correction by lookup into a 2D per-pixel gain texture:
        case {'GainMatrix'}
            % Load our 2D gain correction shader:
            icmShaders = LoadShaderFromFile('ICM2DGainCorrectionShader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);

            % Generate texture handle for fillout later on:
            icmSpec.icmlutid = glGenTextures(1);
            
            % Build config string to bind and use our gain texture:
            icmConfig = sprintf('TEXTURERECT2D(2)=%i', icmSpec.icmlutid);
            
        case {'AnaglyphStereo'}
            % Load the single-view anaglyph shader:
            icmShaders = LoadShaderFromFile('ColoredSingleChannelAnaglyphShader.frag.txt', [], debuglevel);
            icmIdString = sprintf('ICM:%s', icmSpec.type);
            
        otherwise
            error('Unknown type of color correction requested! Internal bug?!?');
    end
    
    
    % Return Vector of shader handles:
    varargout{1} = icmShaders;

    % Return id string:
    varargout{2} = icmIdString;

    % Return shader config string:
    varargout{3} = icmConfig;
    
    % Return override main function definition string:
    varargout{4} = icmOverrideMain;
    
    return;
end

if strcmpi(cmd, 'ApplyPostGLSLLinkSetup')
    if ~specReady
        error('"ApplyPostGLSLLinkSetup" called, but specification of what to postlink is unavailable!');
    end
    
    % Need GL from here on...
    if isempty(GL)
        error('No GL struct defined in "ApplyPostGLSLLinkSetup"?!? This is a bug - Check code!!');
    end
    
    if nargin < 3
        error('Must provide window handle to onscreen window as 2nd argument and viewId as 3rd one!');
    end
    
    % Fetch window handle:
    win = varargin{1};
    viewId = varargin{2};
    
    % Retrieve all params for 'win'dow and given icmSpec, bind shader:
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmSpec.type, viewId); %#ok<*ASGLU,NASGU>
    
    try
        % Setup initial clamping values to valid range 0.0 - 1.0:
        glUniform2f(glGetUniformLocation(glsl, 'ICMClampToColorRange'), 0.0, 1.0);
        
        switch(icmSpec.type)
            case {'ClampOnly', 'CheckOnly', 'None', 'ClampedNoName'}
                % Nothing to do yet...

            case {'SimpleGamma'}
                % Set default encoding gamma for power-law shader to (1.0, 1.0, 1.0):
                glUniform3f(glGetUniformLocation(glsl, 'ICMEncodingGamma'), 1.0, 1.0, 1.0);
                % Default min and max luminance is 0.0 to 1.0, therefore reciprocal 1/range is also 1.0:
                glUniform3f(glGetUniformLocation(glsl, 'ICMMinInLuminance'), 0.0, 0.0, 0.0);
                glUniform3f(glGetUniformLocation(glsl, 'ICMMaxInLuminance'), 1.0, 1.0, 1.0);
                glUniform3f(glGetUniformLocation(glsl, 'ICMReciprocalLuminanceRange'), 1.0, 1.0, 1.0);
                % Default gain to postmultiply is 1.0:
                glUniform3f(glGetUniformLocation(glsl, 'ICMOutputGain'), 1.0, 1.0, 1.0);
                % Default bias to is 0.0:
                glUniform3f(glGetUniformLocation(glsl, 'ICMOutputBias'), 0.0, 0.0, 0.0);
                
            case {'LookupTable'}
                % Set CLUT texture unit to 2:
                glUniform1i(glGetUniformLocation(glsl, 'ICMCLUT'), 2);
                % Setup everything to a pretty meaningless but safe
                % mapping, which will likely just produce all-white,
                % regardless of input:
                glUniform1f(glGetUniformLocation(glsl, 'ICMPrescale'), 1.0);
                glUniform1f(glGetUniformLocation(glsl, 'ICMMaxInputValue'), 1.0);
                % Note that we won't setup the CLUT texture yet. This is a
                % mandatory step after initial setup of the display. We do
                % store the texture id of the clut texture in a permanent
                % location though:
                icmDataForHandle(win, glsl) = icmSpec.icmlutid;
                
            case {'LookupTable3D'}
                % Set CLUT texture unit to 2:
                glUniform1i(glGetUniformLocation(glsl, 'ICMCLUT'), 2);
                % Setup everything to a pretty meaningless but safe
                % mapping, which will likely just produce all-white,
                % regardless of input:
                glUniform1f(glGetUniformLocation(glsl, 'ICMPrescale'), 1.0);
                glUniform1f(glGetUniformLocation(glsl, 'ICMMaxInputValue'), 1.0);
                % Note that we won't setup the CLUT texture yet. This is a
                % mandatory step after initial setup of the display. We do
                % store the texture id of the clut texture in a permanent
                % location though:
                icmDataForHandle(win, glsl) = icmSpec.icmlutid;
                
            case {'GainMatrix'}
                % Set Gain matrix texture unit to 2:
                glUniform1i(glGetUniformLocation(glsl, 'ICMGainField'), 2);
                % Note that we won't setup the Gain texture yet. This is a
                % mandatory step after initial setup of the display. We do
                % store the texture id of the texture in a permanent
                % location though:
                icmDataForHandle(win, glsl) = icmSpec.icmlutid;

            case {'AnaglyphStereo'}
                % Set RedGamma to zero, ie., "disabled" by default:
                glUniform1f(glGetUniformLocation(glsl, 'RedGamma'), 0.0);
                % Set 3x3 transform matrix to identity matrix by default:
                glUniformMatrix3fv(glGetUniformLocation(glsl, 'GainsLeft'), 1, 0, [[1 0 0]; [0 1 0]; [0 0 1]]);

            case {'MatrixMultiply4', 'SensorToPrimary'}
                % Set 4x4 transform matrix to identity matrix by default:
                glUniformMatrix4fv(glGetUniformLocation(glsl, 'M'), 1, 0, [[1 0 0 0]; [0 1 0 0]; [0 0 1 0]; [0 0 0 1]]);

            case {'xyYToXYZ'}
                % Nothing to do for xyYToXYZ colorspace conversion shader.
                
            otherwise
                error('Unknown type of color correction requested! Internal bug?!?');
        end
    catch %#ok<*CTCH>
        % Empty...
        psychrethrow(psychlasterror);
    end
    
    % Unbind shader:
    glUseProgram(0);
    
    return;
end

if strcmpi(cmd, 'SetColorClampingRange')
    
    % Need GL from here on...
    if isempty(GL)
        error('No internal GL struct defined in "SetColorClampingRange" routine?!? This is a bug - Check code!!');
    end
    
    if nargin < 2
        error('Must provide window handle to onscreen window as 2nd argument!');
    end

    if nargin < 4
        error('Must provide minimum and maximum allowable value in 3rd and 4th argument!');
    end
    
    % Fetch window handle:
    win = varargin{1};
    
    % Fetch values:
    minc = varargin{2};
    maxc = varargin{3};
    
    if ~isnumeric(minc) || ~isnumeric(maxc)
        error('(min, max) values must be numbers!');
    end
    
    if minc > maxc
        error('Provided minimum value greater than maximum value - This will not work!');
    end

    if nargin < 5
        viewId = [];
    else
        viewId = varargin{4};
    end
    
    % Retrieve all params for 'win'dow and given icmSpec, bind shader. The
    % 'icmSpec' string is empty - This will expand into the general 'ICM:'
    % string, so we use the first slot/shader with the ICM: token.
    icmId = '';
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId); %#ok<NASGU>
    
    try
        % Setup clamping values to given range:
        glUniform2f(glGetUniformLocation(glsl, 'ICMClampToColorRange'), minc, maxc);        
    catch
        % Empty...
        psychrethrow(psychlasterror);
    end
    
    % Unbind shader:
    glUseProgram(0);
    
    return;
end

if strcmpi(cmd, 'SetEncodingGamma')
    
    % Need GL from here on...
    if isempty(GL)
        error('No internal GL struct defined in "SetEncodingGamma" routine?!? This is a bug - Check code!!');
    end
    
    if nargin < 2
        error('Must provide window handle to onscreen window as 2nd argument!');
    end

    if nargin < 3
        error('Must provide encoding gamma value or vector of 3 encoding gammas for R,G,B in 3rd argument!');
    end
    
    % Fetch window handle:
    win = varargin{1};
    
    % Fetch values:
    gammas = varargin{2};
    
    if ~isnumeric(gammas)
        error('Power law encoding gamma value(s) must be number(s)!');
    end
    
    if length(gammas)~=1 && length(gammas)~=3
        error('Encoding gamma must be a single scalar or a 3 component vector of separate gammas for the Red, Green and Blue color channel!');
    end
    
    if any(min(gammas) < 0) || any(max(gammas) > 1)
        warning(sprintf('At least one of the encoding gamma values %f is outside the sensible range 0.0 - 1.0.\nThis will result in undefined behaviour and is likely not what you want.', gammas)); %#ok<WNTAG,SPWRN>
    end

    if length(gammas) == 1
        % Replicate to all three channels:
        gammas = [gammas, gammas, gammas];
    end
    
    if nargin < 4
        viewId = [];
    else
        viewId = varargin{3};
    end
    
    % Retrieve all params for 'win'dow and given icmSpec, bind shader. The
    % 'icmSpec' string is empty - This will expand into the general 'ICM:'
    % string, so we use the first slot/shader with the ICM: token.
    icmId = 'SimpleGamma';
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId); %#ok<NASGU>
    
    try
        % Set encoding gamma for power-law shader:
        uloc = glGetUniformLocation(glsl, 'ICMEncodingGamma');
        if uloc >= 0
            glUniform3f(uloc, gammas(1), gammas(2), gammas(3));
        else
            error('Tried to set encoding gamma for color correction, but color correction not configured for use of encoding gamma!');
        end
    catch
        % Empty...
        psychrethrow(psychlasterror);
    end
    
    % Unbind shader:
    glUseProgram(0);
    
    return;
end

if strcmpi(cmd, 'SetExtendedGammaParameters')
    
    % Need GL from here on...
    if isempty(GL)
        error('No internal GL struct defined in "SetExtendedGammaParameters" routine?!? This is a bug - Check code!!');
    end
    
    if nargin < 2
        error('Must provide window handle to onscreen window as 2nd argument!');
    end

    if nargin < 3
        error('Must provide minimum input intensity value or vector of 3 minimum intensity values for R,G,B in 3rd argument!');
    end

    if nargin < 4
        error('Must provide maximum input intensity value or vector of 3 maximum intensity values for R,G,B in 4th argument!');
    end

    if nargin < 5
        error('Must provide output gain value or vector of 3 output gain values for R,G,B in 5th argument!');
    end

    if nargin < 6
        error('Must provide output bias value or vector of 3 output bias values for R,G,B in 5th argument!');
    end

    % Fetch window handle:
    win = varargin{1};
    
    % Fetch values:
    minL = varargin{2};
    
    if ~isnumeric(minL)
        error('Minimum input intensity value(s) must be number(s)!');
    end
    
    if length(minL)~=1 && length(minL)~=3
        error('Minimum input intensity value must be a single scalar or a 3 component vector of separate values for the Red, Green and Blue color channel!');
    end
    

    if length(minL) == 1
        % Replicate to all three channels:
        minL = [minL, minL, minL];
    end
    
    % Fetch values:
    maxL = varargin{3};
    
    if ~isnumeric(maxL)
        error('Maximum input intensity value(s) must be number(s)!');
    end
    
    if length(maxL)~=1 && length(maxL)~=3
        error('Maximum input intensity value must be a single scalar or a 3 component vector of separate values for the Red, Green and Blue color channel!');
    end
    

    if length(maxL) == 1
        % Replicate to all three channels:
        maxL = [maxL, maxL, maxL];
    end
    
    % Sanity check:
    if any((maxL - minL) <= 0)
        warning(sprintf('In at least one of the components of the provided minimum and maximum intensity vectors the provided\nminimum is *bigger or equal* than/to the maximum!\nThis will result in undefined behaviour and is likely not what you want.')); %#ok<WNTAG,SPWRN>
    end
    
    % Compute reciprocal:
    recL = 1 ./ abs(maxL - minL);
    
    gain = varargin{4};
    
    if ~isnumeric(gain)
        error('Output gain value(s) must be number(s)!');
    end
    
    if length(gain)~=1 && length(gain)~=3
        error('Output gain value must be a single scalar or a 3 component vector of separate values for the Red, Green and Blue color channel!');
    end
    

    if length(gain) == 1
        % Replicate to all three channels:
        gain = [gain, gain, gain];
    end
    
    % Sanity check:
    if any(gain <= 0)
        warning(sprintf('At least one of the components of the provided gain vector is negative or zero!\nThis will result in undefined behaviour and is likely not what you want.')); %#ok<WNTAG,SPWRN>
    end

    obias = varargin{5};
    
    if ~isnumeric(obias)
        error('Output bias value(s) must be number(s)!');
    end
    
    if length(obias)~=1 && length(obias)~=3
        error('Output bias value must be a single scalar or a 3 component vector of separate values for the Red, Green and Blue color channel!');
    end
    

    if length(obias) == 1
        % Replicate to all three channels:
        obias = [obias, obias, obias];
    end

    if nargin < 7
        viewId = [];
    else
        viewId = varargin{6};
    end
    
    % Retrieve all params for 'win'dow and given icmSpec, bind shader. The
    % 'icmSpec' string is empty - This will expand into the general 'ICM:'
    % string, so we use the first slot/shader with the ICM: token.
    icmId = 'SimpleGamma';
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId); %#ok<NASGU>
    
    try
        % Set parameters for power-law shader:
        uloc = glGetUniformLocation(glsl, 'ICMMinInLuminance');
        if uloc >= 0
            glUniform3f(uloc, minL(1), minL(2), minL(3));
        else
            error('Tried to set extended gamma parameters for color correction, but color correction not configured for use of extended gamma!');
        end
        
        uloc = glGetUniformLocation(glsl, 'ICMMaxInLuminance');
        if uloc >= 0
            glUniform3f(uloc, maxL(1), maxL(2), maxL(3));
        else
            error('Tried to set extended gamma parameters for color correction, but color correction not configured for use of extended gamma!');
        end
        
        uloc = glGetUniformLocation(glsl, 'ICMReciprocalLuminanceRange');
        if uloc >= 0
            glUniform3f(uloc, recL(1), recL(2), recL(3));
        else
            error('Tried to set extended gamma parameters for color correction, but color correction not configured for use of extended gamma!');
        end
        
        uloc = glGetUniformLocation(glsl, 'ICMOutputGain');
        if uloc >= 0
            glUniform3f(uloc, gain(1), gain(2), gain(3));
        else
            error('Tried to set extended gamma parameters for color correction, but color correction not configured for use of extended gamma!');
        end

        uloc = glGetUniformLocation(glsl, 'ICMOutputBias');
        if uloc >= 0
            glUniform3f(uloc, obias(1), obias(2), obias(3));
        else
            error('Tried to set extended gamma parameters for color correction, but color correction not configured for use of extended gamma!');
        end

    catch
        % Empty...
        psychrethrow(psychlasterror);
    end
    
    % Unbind shader:
    glUseProgram(0);
    
    return;
end

if strcmpi(cmd, 'SetMultMatrix4')
    % Need GL from here on...
    if isempty(GL)
        error('SetMultMatrix4: No internal GL struct defined in "SetMultMatrix4" routine?!? This is a bug - Check code!!');
    end
    
    if nargin < 2
        error('SetMultMatrix4: Must provide window handle to onscreen window as 2nd argument!');
    end

    if nargin < 3
        error('SetMultMatrix4: Must provide 4-by-4 color transformation matrix.');
    end
    
    % Fetch window handle:
    win = varargin{1};
    
    % Fetch matrix:
    mat = varargin{2};
    
    if ~isnumeric(mat)
        error('SetMultMatrix4: Matrix must contain numbers!');
    end
    
    if size(mat,1)~=4 || size(mat,2)~=4
        error('SetMultMatrix4: Matrix must have 4 rows and 4 columns!');
    end

    if nargin < 4
        viewId = [];
    else
        viewId = varargin{3};
    end

    icmId = 'MatrixMultiply4';
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId); %#ok<NASGU>
    
    try
        % Set 4x4 transform matrix to identity matrix by default:
        glUniformMatrix4fv(glGetUniformLocation(glsl, 'M'), 1, 0, single(mat));
    catch
        % Empty...
        psychrethrow(psychlasterror);
    end
    
    % Unbind shader:
    glUseProgram(0);

    return;
end

if strcmpi(cmd, 'SetSensorToPrimary')
    % Need GL from here on...
    if isempty(GL)
        error('SetSensorToPrimary: No internal GL struct defined in "SetSensorToPrimary" routine?!? This is a bug - Check code!!');
    end
    
    if nargin < 2
        error('SetSensorToPrimary: Must provide window handle to onscreen window as 2nd argument!');
    end

    if nargin < 3
        error('SetSensorToPrimary: Must provide ''cal'' struct, just as for M-Function ''SensorToPrimary''.');
    end
    
    % Fetch window handle:
    win = varargin{1};
    
    % Fetch calibration struct:
    cal = varargin{2};
    
    if ~isstruct(cal)
        error('SetSensorToPrimary: ''cal'' parameter must be a calibration struct, just as for M-Function ''SensorToPrimary''.');
    end
    
    if ~isfield(cal, 'M_linear_device') || ~isfield(cal, 'ambient_linear')
        error('SetSensorToPrimary: ''cal'' struct is missing at least one of the required M_linear_device or ambient_linear fields.');
    end
    
    % Extract the info we need and validate further (see SensorToPrimary for reference):
    M_linear_device = cal.M_linear_device;
    ambient_linear  = cal.ambient_linear;
    
    if isempty(M_linear_device) || isempty(ambient_linear) || ~isnumeric(M_linear_device) || ~isnumeric(ambient_linear)
        error('SetSensorToPrimary has not been called on valid fields M_linear_device and ambient_linear in calibration structure.');
    end
    
    % Ambient corrections
    [ma,na] = size(ambient_linear);
    if (ma ~= 3) || (na ~= 1)
        error('SetSensorToPrimary: Incorrect dimensions for ambient_linear. Must be a 3 element column vector.');
    end
    
    % Color space conversion
    [mm,nm] = size(M_linear_device);
    if (mm ~= 3) || (nm ~= 3)
        error ('SetSensorToPrimary: Incorrect dimensions for M_linear_device. Must be a 3-by-3 color matrix.');
    end

    % This is what we want to implement, with sensor being the (X,Y,Z) input
    % tristimulus color vector, and primary the output (R,G,B) color vector:
    % sensora = sensor - ambient_linear;
    % primary = M_linear_device * sensora;
    %
    % We do this by a 4x4 matrix vector multiplication:
    %
    % <=> primary = M_linear_device * (sensor - ambient_linear)
    % <=> primary = M_linear_device * sensor - M_linear_device * ambient_linear
    % <=> primary = M_linear_device * sensor - constC with constC = M_linear_device * ambient_linear
    %
    % <=> [ M_linear_device -constC ] with M_linear_device is upper 3x3 matrix embedded in 4 x 4 matrix
    %     [ 0    0     0     1      ] and -constC embedded in right column of 4 x 4 matrix
    %
    %
    % This way the shader just needs to multiply our final matrix with the
    % input color vector.
    
    % Build constC:
    constC = M_linear_device * ambient_linear;
    
    % Build mother 4 x 4 matrix:
    mat = eye(4); % Start off as identity 4x4 matrix.
    mat(1:3, 1:3) = M_linear_device; % Embed upper 3x3 color matrix.
    mat(1:3,   4) = -constC; % Embed -constC 3 element vector.
    
    % mat is all we need (and love of course).
    
    if nargin < 4
        viewId = [];
    else
        viewId = varargin{3};
    end

    icmId = 'SensorToPrimary';
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId); %#ok<NASGU>
    
    try
        % Set 4x4 transform matrix to identity matrix by default:
        glUniformMatrix4fv(glGetUniformLocation(glsl, 'M'), 1, 0, single(mat));
    catch
        % Empty...
        psychrethrow(psychlasterror);
    end
    
    % Unbind shader:
    glUseProgram(0);

    return;    
end

if strcmpi(cmd, 'SetLookupTable')
    
    % Need GL from here on...
    if isempty(GL)
        error('SetLookupTable: No internal GL struct defined in "SetLookupTable" routine?!? This is a bug - Check code!!');
    end
    
    if nargin < 2
        error('SetLookupTable: Must provide window handle to onscreen window as 2nd argument!');
    end

    if nargin < 3
        error('SetLookupTable: Must provide CLUT matrix for which a color lookup table should be built!');
    end
    
    % Fetch window handle:
    win = varargin{1};
    
    % Fetch clut:
    clut = varargin{2};
    
    if ~isnumeric(clut)
        error('SetLookupTable: CLUT matrix must contain number(s)!');
    end
    
    if size(clut,2)~=1 && size(clut,2)~=3
        error('SetLookupTable: Encoding CLUT must be a one column luminance vector or three column matrix for the Red, Green and Blue color channel!');
    end

    if size(clut,1) < 2
        error('SetLookupTable: Encoding CLUT must have at least 2 rows, or at least 2 elements for a luminance vector!');
    end

    if nargin < 4
        viewId = [];
    else
        viewId = varargin{3};
    end
    
    % Optional max input value provided? Assign most common 1.0 if not:
    if nargin < 5 || isempty(varargin{4})
        ICMMaxInputValue = 1.0;
    else
        ICMMaxInputValue = varargin{4};
    end

    % Optional scaling factor provided? Assign proper scaler for clut size,
    % and max input value if not:
    if nargin < 6 || isempty(varargin{5})
        ICMPrescale = ( size(clut,1) -1 ) / ICMMaxInputValue;
    else
        ICMPrescale = varargin{5};
    end

    % Retrieve all params for 'win'dow and our 'LookupTable' icmSpec, bind shader.
    icmId = 'LookupTable';
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId); %#ok<NASGU>
    
    try
        % Setup initial clamping values to valid range 0.0 - maximum in passed CLUT:
        glUniform2f(glGetUniformLocation(glsl, 'ICMClampToColorRange'), 0.0, max(max(clut)));

        % Setup max input value and prescaler:
        glUniform1f(glGetUniformLocation(glsl, 'ICMMaxInputValue'),ICMMaxInputValue);
        glUniform1f(glGetUniformLocation(glsl, 'ICMPrescale'),ICMPrescale);        
    catch
        % Empty...
        psychrethrow(psychlasterror);
    end
    
    % Unbind shader:
    glUseProgram(0);
    
    if isempty(icmDataForHandle) || size(icmDataForHandle, 1) < win || size(icmDataForHandle, 2) < glsl || Screen('WindowKind',win) ~= 1
        error('SetLookupTable: Tried to assign clut to a non-onscreen window or one which doesn''t have "LookupTable" based color correction enabled!');
    end
    
    % Convert 'clut' to single(), so it is a float format for OpenGL:
    clut = single(clut);
    
    % Try to encode in highest precision format that the hardware supports:
    winfo = Screen('GetWindowInfo', win);
    if winfo.GLSupportsTexturesUpToBpc >= 32
        % Full 32 bits single precision float:
        internalFormat = GL.LUMINANCE_FLOAT32_APPLE;
        if verbosity >= 3, fprintf('PsychColorCorrection: Using a 32 bit float CLUT -> 23 bits effective linear output precision for color correction.\n'); end
    else
        % No float32 textures:
        if (winfo.GLSupportsTexturesUpToBpc >= 16)
            % Choose 16 bpc float textures:
            internalFormat = GL.LUMINANCE_FLOAT16_APPLE;
            if verbosity >= 3, fprintf('PsychColorCorrection: Using a 16 bit float CLUT -> 10 bits effective linear output precision for color correction.\n'); end
        else
            % No support for > 8 bpc textures at all and/or no need for
            % more than 8 bpc precision or range. Choose 8 bpc texture:
            internalFormat = GL.LUMINANCE;
            if verbosity >= 3, fprintf('PsychColorCorrection: Using a 8 bit integer CLUT -> 8 bits effective linear output precision for color correction.\n'); end

            % Plain old 8 bits fixed point:
            if (max(max(clut)) > 1) || (min(min(clut)) < 0)
                % Ohoh, out of range values for integer texture!
                if verbosity >= 2
                    fprintf('\nWARNING!PsychColorCorrection: CLUT contains values greater than 1.0 or negative values, which your hardware can''t handle!!\n');
                    fprintf('WARNING!PsychColorCorrection: This will likely cause remapping artifacts in color correction!!\n');
                end
            end
            
            if (winfo.BitsPerColorComponent > 8)
                if verbosity >= 2
                    fprintf('WARNING!PsychColorCorrection: Your hardware can only handle 8 bit precision color correction in outputrange 0.0 - 1.0,\n');
                    fprintf('WARNING!PsychColorCorrection: but your framebuffer is configured for more than 8 bit precision. This may cause\n');
                    fprintf('WARNING!PsychColorCorrection: loss of effective precision in color correction and thereby unwanted artifacts!\n');
                end
            end
        end
    end
    
    if size(clut,1) > glGetIntegerv(GL.MAX_RECTANGLE_TEXTURE_SIZE_ARB)
        error('SetLookupTable: Tried to assign a clut with %i slots. This is more than your graphics hardware can handle! [Maximum is %i slots].', size(clut,1), glGetIntegerv(GL.MAX_RECTANGLE_TEXTURE_SIZE_ARB));
    end
    
    % Bind relevant texture object:
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, icmDataForHandle(win, glsl));
    
    % Set filters properly: Want nearest neighbour filtering, ie., no filtering
    % at all. We'll do our own linear filtering in the ICM shader. This way
    % we can provide accelerated linear interpolation on all GPU's with all
    % texture formats, even if GPU's are old:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);

    % Want clamp-to-edge behaviour to saturate at minimum and maximum
    % intensity value, and to make sure that a pure-luminance 1 row clut is
    % properly "replicated" to all three color channels in rgb modes:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP_TO_EDGE);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP_TO_EDGE);
    
    % Assign lookuptable data to texture:
    glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, internalFormat, size(clut, 1), size(clut, 2), 0, GL.LUMINANCE, GL.FLOAT, clut);
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);

    % Done.
    return;
end

if strcmpi(cmd, 'SetLookupTable3D')
    
    % Need GL from here on...
    if isempty(GL)
        error('SetLookupTable3D: No internal GL struct defined in "SetLookupTable3D" routine?!? This is a bug - Check code!!');
    end
    
    if nargin < 2
        error('SetLookupTable3D: Must provide window handle to onscreen window as 2nd argument!');
    end

    if nargin < 3
        error('SetLookupTable3D: Must provide CLUT 4D matrix for which a 3D color lookup table should be built!');
    end
    
    % Fetch window handle:
    win = varargin{1};
    
    % Fetch clut:
    clut = varargin{2};
    
    if ~isnumeric(clut)
        error('SetLookupTable3D: CLUT matrix must contain number(s)!');
    end
    
    if ndims(clut) ~= 4
        error('SetLookupTable3D: Encoding CLUT must be a 4D matrix with indices (channel, redinput, greeninput, blueinput)!');
    end
    
    if size(clut,1)~=3
        error('SetLookupTable3D: Encoding CLUT must be a 4D matrix whose first dimension has a size of 3 elements for (R,G,B) color triplets!');
    end

    if size(clut,2) < 1 || size(clut,3) < 1 || size(clut,4) < 1
        error('SetLookupTable3D: Encoding 4D CLUT must have at least one element in each of the three input color dimensions 2,3 and 4!');
    end

    if nargin < 4
        viewId = [];
    else
        viewId = varargin{3};
    end
    
    % Optional max input value provided? Assign most common 1.0 if not:
    if nargin < 5 || isempty(varargin{4})
        ICMMaxInputValue = 1.0;
    else
        ICMMaxInputValue = varargin{4};
    end

    % Optional scaling factor provided? Assign proper scaler for normalized
    % 3D texture coordinate range of clut texture and max input value if not:
    if nargin < 6 || isempty(varargin{5})
        ICMPrescale = 1.0 / ICMMaxInputValue;
    else
        ICMPrescale = varargin{5};
    end

    % Optional precision spec provided? Assign if so.
    if nargin < 7 || isempty(varargin{6})
        % Default to precision zero for simple 8 bit integer precision -
        % sufficient to drive a standard 8 bpc integer framebuffer.
        % Our strategy for 3D LUT's is different than for normal linear
        % LUT's, where we choose the highest precision (32 bpc float) by
        % default. The reason is that memory consumption of 3D textures
        % increases with O(n^3) for n levels per input dimension, as
        % opposed to O(n) for linear LUT's. Choosing higher precisions by
        % default could get us into a out-of-memory condition very quickly.
        % Even if the texture fits into VRAM, a 3D texture used for color
        % lookup and applied per pixel could quickly trash the texture
        % caches and cause a breakdown of performance.
        precision = 0;
    else
        precision = varargin{6};
    end
    
    % Trilinear interpolation switch provided?
    if nargin < 8 || isempty(varargin{7})
        % Default to trilinear interpolation:
        interpolate = 1;
    else
        interpolate = varargin{7};
    end
    
    % Retrieve all params for 'win'dow and our 'LookupTable' icmSpec, bind shader.
    icmId = 'LookupTable3D';
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId); %#ok<NASGU>
    
    try
        % Setup initial clamping values to valid range 0.0 - maximum in passed CLUT:
        glUniform2f(glGetUniformLocation(glsl, 'ICMClampToColorRange'), 0.0, max(max(max(max(clut)))));

        % Setup max input value and prescaler:
        glUniform1f(glGetUniformLocation(glsl, 'ICMMaxInputValue'), ICMMaxInputValue);
        glUniform1f(glGetUniformLocation(glsl, 'ICMPrescale'), ICMPrescale);        
    catch
        % Empty...
        psychrethrow(psychlasterror);
    end
    
    % Unbind shader:
    glUseProgram(0);
    
    if isempty(icmDataForHandle) || size(icmDataForHandle, 1) < win || size(icmDataForHandle, 2) < glsl || Screen('WindowKind',win) ~= 1
        error('SetLookupTable3D: Tried to assign clut to a non-onscreen window or one which doesn''t have "LookupTable3D" based color correction enabled!');
    end
    
    % Convert 'clut' to single(), so it is a float format for OpenGL:
    clut = single(clut);
    
    % Try to encode in highest precision format that the hardware supports, but no more than what is requested by usercode:
    winfo = Screen('GetWindowInfo', win);
    if (winfo.GLSupportsTexturesUpToBpc >= 32) && (winfo.GLSupportsFilteringUpToBpc >= 32 || interpolate == 0) && (precision >= 2)
        % Full 32 bits single precision float with linear filtering:
        internalFormat = GL.RGB_FLOAT32_APPLE;
        if verbosity >= 3, fprintf('PsychColorCorrection: Using a 32 bit float 3D CLUT -> 23 bits effective linear output precision for color correction.\n'); end
    else
        % No float32 textures with linear filtering:
        if (winfo.GLSupportsTexturesUpToBpc >= 16) && (winfo.GLSupportsFilteringUpToBpc >= 16 || interpolate == 0) && (precision >= 1)
            % Choose 16 bpc float textures with linear filtering:
            internalFormat = GL.RGB_FLOAT16_APPLE;
            if verbosity >= 3, fprintf('PsychColorCorrection: Using a 16 bit float 3D CLUT -> 10 bits effective linear output precision for color correction.\n'); end
        else
            % No support for > 8 bpc textures with linear filtering at all
            % and/or no need for more than 8 bpc precision or range. Choose
            % 8 bpc texture:
            internalFormat = GL.RGB;
            if verbosity >= 3, fprintf('PsychColorCorrection: Using a 8 bit integer 3D CLUT -> 8 bits effective linear output precision for color correction.\n'); end

            % Plain old 8 bits fixed point:
            if (max(max(max(max(clut)))) > 1) || (min(min(min(min(clut)))) < 0)
                % Ohoh, out of range values for integer texture!
                if verbosity >= 2
                    fprintf('\nWARNING!PsychColorCorrection: 3D CLUT contains values greater than 1.0 or negative values, which your hardware can''t handle!!\n');
                    fprintf('WARNING!PsychColorCorrection: This will likely cause remapping artifacts in color correction!!\n');
                end
            end
            
            if (winfo.BitsPerColorComponent > 8)
                if verbosity >= 2
                    fprintf('WARNING!PsychColorCorrection: Your hardware can only handle 8 bit precision color correction in outputrange 0.0 - 1.0,\n');
                    fprintf('WARNING!PsychColorCorrection: but your framebuffer is configured for more than 8 bit precision. This may cause\n');
                    fprintf('WARNING!PsychColorCorrection: loss of effective precision in color correction and thereby unwanted artifacts!\n');
                end
            end
        end
    end
    
    max3DSize = glGetIntegerv(GL.MAX_3D_TEXTURE_SIZE);
    if size(clut,2) > max3DSize || size(clut,3) > max3DSize || size(clut,4) > max3DSize
        error('SetLookupTable3D: Tried to assign a clut with (%i,%i,%i) slots. This is more than your graphics hardware can handle! [Maximum is %i slots].', size(clut,2), size(clut,3), size(clut,4), max3DSize);
    end
    
    % Bind relevant texture object:
    glBindTexture(GL.TEXTURE_3D, icmDataForHandle(win, glsl));
    
    % Set filters properly:
    if interpolate > 0
        % We want trilinear filtering:
        glTexParameteri(GL.TEXTURE_3D, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
        glTexParameteri(GL.TEXTURE_3D, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
    else
        % We want nearest neighbour lookup:
        glTexParameteri(GL.TEXTURE_3D, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
        glTexParameteri(GL.TEXTURE_3D, GL.TEXTURE_MAG_FILTER, GL.NEAREST);
    end
    
    % Want clamp-to-edge behaviour to saturate at minimum and maximum
    % color values:
    glTexParameteri(GL.TEXTURE_3D, GL.TEXTURE_WRAP_S, GL.CLAMP_TO_EDGE);
    glTexParameteri(GL.TEXTURE_3D, GL.TEXTURE_WRAP_T, GL.CLAMP_TO_EDGE);
    glTexParameteri(GL.TEXTURE_3D, GL.TEXTURE_WRAP_R, GL.CLAMP_TO_EDGE);
    
    % Clear any pending OpenGL errors:
    while glGetError(); end;
    
    % Assign lookuptable data to texture:
    glTexImage3D(GL.TEXTURE_3D, 0, internalFormat, size(clut, 2), size(clut, 3), size(clut, 4), 0, GL.RGB, GL.FLOAT, clut);
    
    % Check and handle errors:
    err = glGetError;
    if err ~= GL.NO_ERROR
        if err == GL.INVALID_VALUE
            % Can only get a invalid_value if power-of-two constraint is
            % violated on a GPU which does not support NPOT textures:
            error('SetLookupTable3D: LUT creation failed. Most likely cause: Your graphics card requires all 3D-CLUT dimensions to have a number of elements which is a power of two.\n');            
        elseif err == GL.OUT_OF_MEMORY
            error('SetLookupTable3D: LUT creation failed. Your graphics card has insufficient free memory to handle a CLUT of this size. Try to reduce CLUT size.');
        else
            % Something else went bonkers:
            error('SetLookupTable3D: OpenGL error during lookup table texture creation: %s\n', gluErrorString(err));
        end
    end
    
    % Done, LUT is ready, unbind:
    glBindTexture(GL.TEXTURE_3D, 0);

    % Done.
    return;
end

if strcmpi(cmd, 'SetGainMatrix')
    
    % Need GL from here on...
    if isempty(GL)
        error('SetGainMatrix: No internal GL struct defined in "SetGainMatrix" routine?!? This is a bug - Check code!!');
    end
    
    if nargin < 2
        error('SetGainMatrix: Must provide window handle to onscreen window as 2nd argument!');
    end

    if nargin < 3
        error('SetGainMatrix: Must provide 2D gain matrix for use!');
    end
        
    % Fetch window handle:
    win = varargin{1};
    
    % Fetch clut:
    clut = varargin{2};
    
    if ~isnumeric(clut)
        error('SetGainMatrix: Gain matrix must contain number(s)!');
    end

    [w, h] = Screen('WindowSize', win);
    if size(clut,2)~=w || size(clut,1)~=h
        error('SetGainMatrix: Gain matrix must have the same size (width x height) in pixels as the target onscreen window!');
    end

    if nargin < 4
        viewId = [];
    else
        viewId = varargin{3};
    end
    
    if nargin < 5
        precision = [];
    else
        precision = varargin{4};
    end
    
    if isempty(precision)
        % Default precision for gain matrix is 2 == max. precision 32 bpc.
        precision = 2;
    end
    
    % Retrieve all params for 'win'dow and our 'LookupTable' icmSpec, bind shader.
    icmId = 'GainMatrix';
    [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId); %#ok<NASGU>
    
    %     % Not used yet:
    %     try
    %         % Setup initial clamping values to valid range 0.0 - maximum in passed CLUT:
    %         glUniform2f(glGetUniformLocation(glsl, 'ICMClampToColorRange'), 0.0, max(max(clut)));
    %     catch
    %         % Empty...
    %         psychrethrow(psychlasterror);
    %     end
    
    % Unbind shader:
    glUseProgram(0);
    
    if isempty(icmDataForHandle) || size(icmDataForHandle, 1) < win || size(icmDataForHandle, 2) < glsl || Screen('WindowKind',win) ~= 1
        error('SetGainMatrix: Tried to assign matrix to a non-onscreen window or one which doesn''t have "GainMatrix" based color correction enabled!');
    end

    % Convert 'clut' to single(), so it is a float format for OpenGL. Also
    % need to transpose and flip from Matlab col-major to GL row-major:
    ch = size(clut, 3);
    for i=1:ch
        tclut(:,:,i) = transpose(flipud(clut(:,:,i))); %#ok<AGROW>
    end
    tclut = single(tclut);

    % Interleave RGB gain info if this is a 3 channel gain matrix:
    clut = zeros(size(tclut, 3), size(tclut, 1), size(tclut, 2), 'single');
    for i=1:size(tclut, 3)
        clut(i,:,:) = tclut(:,:,i);
    end;
    
    % Try to encode in highest precision format that the hardware supports:
    winfo = Screen('GetWindowInfo', win);
    if (winfo.GLSupportsTexturesUpToBpc >= 32) && (precision >= 2)
        % Full 32 bits single precision float:
        internalFormat = GL.LUMINANCE_FLOAT32_APPLE;
        if ch == 3
            internalFormat = GL.RGB_FLOAT32_APPLE;
        end
        if verbosity >= 3, fprintf('PsychColorCorrection: Using a 32 bit float matrix -> 23 bits (6 decimal digits) effective linear precision for color correction gain matrix.\n'); end
    else
        % No float32 textures:
        if (winfo.GLSupportsTexturesUpToBpc >= 16) && (precision >= 1)
            % Choose 16 bpc float textures:
            internalFormat = GL.LUMINANCE_FLOAT16_APPLE;
            if ch == 3
                internalFormat = GL.RGB_FLOAT16_APPLE;
            end
            if verbosity >= 3, fprintf('PsychColorCorrection: Using a 16 bit float matrix -> 10 bits (3 decimal digits) effective linear precision for color correction gain matrix.\n'); end
        else
            % No support for > 8 bpc textures at all and/or no need for
            % more than 8 bpc precision or range. Choose 8 bpc texture:
            internalFormat = GL.LUMINANCE;
            if ch == 3
                internalFormat = GL.RGB8;
            end
            
            if verbosity >= 3
                fprintf('PsychColorCorrection: Using a 8 bit integer matrix -> 8 bits (2 decimal digits) effective linear precision for color correction gain matrix.\n');
                fprintf('PsychColorCorrection: Gain values will be restricted to range 0.0 - 1.0, with 256 levels, ie. steps of 1/256.\n');
            end
            
            % Plain old 8 bits fixed point:
            if (max(max(max(clut))) > 1) || (min(min(min(clut))) < 0)
                % Ohoh, out of range values for integer texture!
                if verbosity >= 2 
                    fprintf('\nWARNING!PsychColorCorrection: Matrix contains values greater than one or negative values, which your hardware can''t handle!!\n');
                    fprintf('WARNING!PsychColorCorrection: This will likely cause remapping artifacts in gain correction!!\n');
                end
            end
            
            if (winfo.BitsPerColorComponent > 8)
                if verbosity >= 2
                    fprintf('WARNING!PsychColorCorrection: Your hardware can only handle 8 bit precision gain correction width 256 discrete levels,\n');
                    fprintf('WARNING!PsychColorCorrection: but your framebuffer is configured for more than 8 bit precision. This may cause\n');
                    fprintf('WARNING!PsychColorCorrection: loss of effective precision in gain correction and thereby unwanted artifacts!\n');
                end
            end
        end
    end
    
    % Bind relevant texture object:
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, icmDataForHandle(win, glsl));
    
    % Set filters properly: Want nearest neighbour filtering, ie., no filtering at all:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);

    % Want clamp-to-edge behaviour to saturate at minimum and maximum gain value:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP_TO_EDGE);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP_TO_EDGE);
    
    % Assign 2D matrix data to texture:
    if ch == 3
        glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, internalFormat, size(clut, 2), size(clut, 3), 0, GL.RGB, GL.FLOAT, clut);
    else
        glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, internalFormat, size(clut, 2), size(clut, 3), 0, GL.LUMINANCE, GL.FLOAT, clut);
    end

    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);

    % Done.
    return;
end


error('Unknown subfunction specified. Typo?!? Read the help.');


function [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, viewId)

    % Fetch slot, and GLSL handle to relevant imaging slot and shader:
    if isempty(viewId)
        viewId = 'AllViews';
    end 
    
    % MK Resolved 26.4.2010: FIXME: HACK FOR BUG IN IMG PIPE!! 
    if strcmpi(viewId, 'AllViews') || strcmpi(viewId, 'FinalFormatting')
    %if strcmpi(viewId, 'FinalFormatting')
        chain = 'FinalOutputFormattingBlit';
    end
    
    % MK Resolved 26.4.2010: FIXME: HACK FOR BUG IN IMG PIPE!! if strcmpi(viewId, 'LeftView')
    if strcmpi(viewId, 'LeftView') %|| strcmpi(viewId, 'AllViews')
        chain = 'StereoLeftCompositingBlit';
    end

    if strcmpi(viewId, 'RightView')
        chain = 'StereoRightCompositingBlit';
    end

    if strcmpi(icmId, 'ClampedNoName')
        icmIdString = sprintf('CMIH:ClampedNoName');
    else
        icmIdString = sprintf('ICM:%s', icmId);
    end
    
    [slot shaderid blittercfg voidptr glsl luttexid] = Screen('HookFunction', win, 'Query', chain, icmIdString);

    % Shader found?
    if slot == -1
        % MK Resolved 26.4.2010: Doesn't hurt here. FIXME: HACK FOR BUG IN IMG PIPE!!
        if strcmpi(viewId, 'AllViews')
            % Special case 'AllViews': We searched 'LeftView' and failed.
            % Let's retry with 'FinalFormatting' view before we give up:
            [slot shaderid blittercfg voidptr glsl luttexid] = GetSlotForTypeAndBind(win, icmId, 'FinalFormatting');
            % If we make it to this point then our search was successfull:
            return;
        end
        
        fprintf('PsychColorCorrection: Error: Failed to find plugin while searching chain: %s, searched icmIdString: %s\n', chain, icmIdString);
        error('Could not find shader plugin for color correction inside imaging pipeline for window and view! Is color correction really enabled for this window and view channel?!?');
    end
    
    if glsl == 0
        error('Color correction shader is not operational for unknown reason, maybe a Psychtoolbox bug? Sorry...');
    end
    
    % Bind it:
    glUseProgram(glsl);
return;
