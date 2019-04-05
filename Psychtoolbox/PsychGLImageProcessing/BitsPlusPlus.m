function [win, winRect] = BitsPlusPlus(cmd, arg, dummy, varargin)
% BitsPlusPlus(cmd [, arg1][, arg2, ...]) -- Psychtoolbox interface to
% Cambridge Research Systems Bits++ and Bits# boxes for high precision
% stimulus output to analog displays via 14 bit video converters.
%
% This function is used to set up and interface with the Bits++ / Bits# box of
% CRS. It is a Matlab wrapper around lower level GLSL Psychtoolbox
% functions. This function depends on graphics hardware that supports the
% Psychtoolbox imaging pipeline and framebuffers with more than 14 bit
% precision, i.e., 16 bpc fixed point framebuffers or floating point
% framebuffers. Have a look at the Psychtoolbox Wiki where you can find a
% list of graphics cards with the neccessary features.
%
% This function supersedes the old Matlab based Bits++ Toolbox, which
% essentially provides the same functionality on any graphics card, but is
% harder to use, much slower and not fully integrated into PTB, ie, you
% can't take full advantage of PTB's advanced drawing and image processing
% functions when using the old Bits++ toolbox.
%
% See the help section about Bits# for advanced Bits# commands and how to
% establish communication and convenient control for the Bits# via USB
% connection.
%
% cmd - The command that BitsPlusPlus should execute. cmd can be any of
% the following:
%
%
% Load a linear identity mapping CLUT into Bits++ while running in Bits++
% mode:
%
% BitsPlusPlus('LoadIdentityClut', window);
%
% Will schedule update to an identity clut. Next invocation of
% Screen('Flip', window); will actually upload the identity clut into
% Bits++.
%
%
% BitsPlusPlus('UseFE1StereoGoggles', window, enable [, shutterState=0]);
% If 'enable' is 1, enable use of the CRS FE1 stereo goggles connected
% to the Bits display device. If 'enable' is 0 then disable use of the
% goggles again and put them in a resting state with either both shutters
% open ('shutterState = 0') or closed ('shutterState = 1').
%
%
% Schedule a Bits++ DIO command for execution on next Screen('Flip'):
%
% BitsPlusPlus('DIOCommand', window, repetitions, mask, data, command [, xpos, ypos]);
%
% This will draw the proper T-Lock control codes at positions (xpos, ypos)
% for execution of the Bits++ DIO commands (mask, data, command).
%
% You can specify multiple codes at once: If mask, data, command, xpos and
% ypos are matrices or vectors with 'num' rows, then each of the 'num' rows
% defines one T-Lock code. If mask, command, xpos and ypos are scalars and
% data is a one row vector, then only the corresponding T-Lock line is
% drawn.
%
% For each DIO command:
% 'mask' must be a 8 bit integer value, 'command' must be a 8 bit integer
% value, whereas 'data' must be a a 248 element row vector of bytes. 
%
% Consult your Bits++ manual for explanation of the meaning of the values.
%
% xpos and ypos are optional: By default, the T-Lock code is drawn into the
% 3rd pixel row of the output image, so it can't collide with a potential
% T-Lock code for CLUT updates.
%
% The DIO command will become effective during the next flip command. The
% T-Lock code is drawn during 'repetitions' successive invocations of
% Screen('Flip'). If you set 'repetitions' to -1, then the code will be
% drawn until you stop it via a call to BitsPlusPlus('DIOCommandReset', window);
%
%
% Disable use of the DIO T-Lock code blitting:
%
% BitsPlusPlus('DIOCommandReset', window);
% Stops blitting of T-Lock command codes immediately. If you want to use
% them again, you have to respecify codes via the BitsPlusPlus('DIOCommand',...);
%
%
% Manage treatment of one or two connected CRS devices:
%
% BitsPlusPlus('SetDualDevices', mode);
% If 'mode' is set to 0, which is the default, assume one CRS device is connected.
%
% If 'mode' is set to 1, then assume two CRS devices, e.g., 2 Display++ or
% Bits# devices are connected for dual-display stereo/binocular stimulation and
% adapt accordingly. Apply the same CLUT in Bits++ and Mono++ display mode
% to both connected display devices in stereomodes 4, 5 and 10.
%
%
%
%
% Open a full-screen window on the Bits++ display as with
% Screen('OpenWindow', ...), perform all initialization:
%
% The following commands will execute Screen('OpenWindow') with all proper
% parameters, followed by Bits++ init routines. They are completely sufficient
% drop in replacements for Screen('OpenWindow'), accepting and returning
% exactly the same arguments that Screen() would do, adjusting all
% parameters to the constraints of the Bits++, if necessary.
%
% Activate Bits++ mode:
%
% [win, winRect] = BitsPlusPlus('OpenWindowBits++', screenid, ...);
%
% This will open an onscreen window on Bits++ screen 'screenid' with a
% standard 8 bits per color channel framebuffer. The gamma table of your
% graphics hardware will be loaded with an identity gamma table, so the
% T-Lock system of Bits++ works and Bits++ can accept commmands embedded
% into the stimulus images. Psychtoolbox will automatically embed a T-Lock
% control line into the top line of the display screen, which encodes the
% 256 entry, 14 bit per color channel CLUT to use for Bits++ display mode.
% You can change the Bits++ CLUT at any time via the standard PTB
% Screen('LoadNormalizedGammaTable', win, newclut, 2); command. The
% 'newclut' will get uploaded to the Bits++ at the next invocation of
% Screen('Flip') to allow updates of the CLUT synchronous to stimulus
% updates. 'newclut' has to be a 256 rows by 3 columns matrix with values
% in range 0.0 - 1.0: 0.0 is mapped to clut color value 0, 1.0 is mapped to
% the highest Bits++ output color value 16383.
%
% This mode works on any OpenGL graphics hardware.
%
%
% Activate Mono++ mode:
%
% [win, winRect] = BitsPlusPlus('OpenWindowMono++', screenid, ...);
%
% This will open an onscreen window on Bits++ screen 'screenid' for display
% of pure luminance (grayscale) images. The framebuffer has a resolution of
% 32 bit floating point precision by default: This means that pixel luminance
% values have to be specified as floating point numbers between 0.0 and
% 1.0. 0.0 maps to black (Output intensity 0 on Bits++ device). 1.0 maps to
% white (Maximum output intensity 16383 on Bits++ device). The intensity
% range between 0.0 and 1.0 is internally represented and processed by
% Psychtoolbox with a resolution 23 bits, i.e. over 8 million levels of
% luminance. The Bits++ can resolve this range to 14 bits, ie. 16384 levels
% of luminance during display. This mode is not compatible with the use of
% any gamma- or clut- tables. Both the graphics hardwares gamma table and
% the Bits++ internal clut are set to an identity mapping while this mode
% is active. Please read the notes below the Color++ section for graphics
% hardware requirements and other useful tips for use of Bits++.
%
% If you call this subfunction as 'OpenWindowMono++WithOverlay', the
% overlay plane of Bits++ gets enabled and an additional overlay window is
% created for drawing the image for that overlay plane.
%
% [overlaywin, overlaywinRect] = BitsPlusPlus('GetOverlayWindow', win);
% - Will return the handle to the 'overlaywin'dow associated with the
% onscreen luminance window:
%
%   'overlayWin' is the handle to the overlay window associated with the
%   overlay of onscreen window 'win'. The overlay window is a standard
%   offscreen window, so you can do anything with it that you would want to
%   do with offscreen windows. The only difference is that the window is a
%   pure index window: It only has one "color channel", which can be written
%   with color values between 0 and 255. Values 1 to 255 get mapped to the
%   corresponding color indices of the Bits++ overlay plane: A zero value is
%   transparent -- Content of the onscreen window is visible. Positive
%   non-zero color values map to the 255 indices available in overlay mode,
%   these get mapped by the Bits++ CLUT to colors. You can define the
%   mapping of indices to CLUT colors via the
%   Screen('LoadNormalizedGammaTable', win, clut, 2); command.
%
%   Updates of the overlay image are synchronized to Screen('Flip')
%   updates. If you draw into the overlay window, the changed overlay image
%   will become visible at Screen('Flip') time -- in sync with the changed
%   onscreen window content. The overlay plane is not automatically cleared
%   to background (or transparent) color after a flip, but its content
%   persists across flips. You need to clear it out manually via a
%   Screen('FillRect') command.
%
%
% Activate Color++ mode:
%
% [win, winRect] = BitsPlusPlus('OpenWindowColor++', screenid, ...);
%
% This will open an onscreen window on Bits++ screen 'screenid' for display
% of 14 bit per color component 42bpp color images. The framebuffer has
% a resolution of 32 bit floating point precision for each color component
% by default: This means that (Red, Green, Blue) color pixel component
% values have to be specified as floating point numbers between 0.0 and
% 1.0. 0.0 maps to minimum output intensity on Bits++ device for a channel.
% 1.0 maps to maximum output intensity 16383 on Bits++ device for a channel.
% The color intensity range between 0.0 and 1.0 is internally represented and
% processed by Psychtoolbox with an effective resolution of about 23 bits,
% i.e. over 8 million levels of color per color channel. The Bits++ can resolve
% this range to 14 bits, ie. 16384 levels of color during display. This mode
% is not compatible with the use of any gamma- or clut- tables. Both the graphics
% hardwares gamma table and the Bits++ internal clut are set to an identity
% mapping while this mode is active. Please read the notes below for graphics
% hardware requirements and other useful tips for use of Bits++.
%
%
% If you use Color++ mode, you must call
% BitsPlusPlus('SetColorConversionMode', mode); first to select the mode
% for sampling the framebuffer and converting into output color values. See
% the respective section of "help PsychImaging" for 'Color++' or
% 'EnableDataPixxC48Output' mode for the Bits+ or Datapixx device for an
% explanation of this mandatory parameter. The setting before 22nd
% September 2010 for all PTB-3 versions was 0 (==zero).
%
% You can query the mode for an onscreen window 'win' by a call to:
% mode = BitsPlusPlus('GetColorConversionMode', win);
%
%
% Notes for both Mono++ and Color++ mode:
%
% In Mono++ and Color++ mode, PTB expects color values in the range 0.0 to
% 1.0 instead of the (otherwise usual) range 0 to 255. The range 0.0-1.0
% is a more natural fit for high dynamic range/precision output devices than
% the 0-255 range with its strong ties to 8 bpc output devices. 0-1 is also
% the "natural" native color range of OpenGL, so colors in this range can
% be handled by the graphics hardware at a higher speed. You can change the
% mapping of input colors to output intensities by use of the command
% Screen('ColorRange') (see its online help for usage), but in the interest
% of uniform code and to avoid possible side effects with some graphics
% hardware, we strongly recommend using the default 0.0-1.0 color range.
%
% You can still pass standard 8bpc (256 color/intensity levels) color/luminance
% textures to PTB via standard use of Screen('MakeTexture') - the hardware
% will convert such 8bpc images to OpenGL's native color range, as well as
% any images delivered by the Quicktime movie playback engine or the video
% capture engine. If you want to provide high dynamic range, high color
% depths images, please specify them as Matlab double matrices to
% Screen('MakeTexture') and set the optional flag 'floatprecision' to 1 or
% 2, i.e., hdrtex = Screen('MakeTexture', win, myHDRImage, [], [], 2);
%
% Psychtoolbox will represent such images with an internal precision of 10
% bits + 1 bit sign if you choose the 'floatprecision' flag to be 1. If you
% choose a 'floatprecision' flag of 2, PTB will represent the images with
% an internal precision of 23 bits + 1 bit sign. You can provide negative
% color values as well, e.g., -0.5. If you wonder what the use of this
% might be, have a careful look at the tutorial script...
% 'AdditiveBlendingForLinearSuperpositionTutorial.m'
% ... for an example of extremely fast drawing of luminance gratings with
% controllable size, orientation and contrast and correct linear superposition:
% 
% By default, PTB will use a 32-bit floating point framebuffer for your
% drawings, ie. the precision is way higher than needed for any high
% dynamic range/resolution display device in existence. The downside of this
% super-precision is that alpha-blending is not supported in this mode, unless
% you employ an NVidia Geforce 8000 series (and later) graphics card, or a
% ATI Radeon HD2000/3000 series graphics card (and later). If you need
% alpha-blending on older/other hardware then specify the optional flag
% 'kPsychNeed16BPCFloat' for the 'imagingmode' argument. This will reduce
% effective accuracy of the framebuffer to 10 bit precision, but allow for
% fast alpha-blending. 10 Bit precision are 4 bits less than the 14 bits
% that Bits++ can provide, but it will be possible to use the extra 14-10 =
% 4 bits for gamma correction of the display by employing a gamma
% correction shader.
%
% Gamma- and color correction:
%
% In Mono++ and Color++ mode, the hardware gamma tables of your graphics
% card and the Bits+ box can't be used for gamma- or color correction.
% However, PTB provides a much more powerful and flexible color correction
% system for this purpose. See "help PsychColorCorrection" for further
% explanation and usage examples for standard gamma correction.
%
% Graphics hardware requirements: Mono++ and Color++ mode require use of the
% Psychtoolbox imaging pipeline and floating point framebuffers. The
% minimum requirements are ATI Radeon X1000 series or NVidia Geforce-6800
% series and later graphics hardware. We currently recommend NVidia
% Geforce-8000 series or ATI Radeon HD-2000/3000 hardware for best results.
% However, this functions have been successfully tested on ATI Radeon X1600
% and NVidia Geforce-7800 hardware as well.
%
% All Bits++ modes supported by this function should work Plug & Play,
% requiring no changes to your stimulus code other than mentioned here to
% take full advantage of all functionality of Psychtoolbox just as with standard
% 8 bpc displays at the higher 14 bpc quality of Bits++. If you find this
% not to be the case then it's either an omission in our documentation
% or a bug - Please report it.
%
%
% BITS# specific functions:
%
% A Bits# device which is connected via USB will show up as an additional
% serial port on the system. This driver will communicate with the Bits#
% by establishing a serial port connection to the device via that serial
% port. Presence of a Bits# can be signalled by either calling the BitsPlusPlus('OpenBits#')
% function and passing a serial port device spec 'portSpec', or just by calling
% BitsPlusPlus('OpenBits#') without any parameters. In the latter case, the
% driver will check for the existence of a configuration file named...
% [PsychtoolboxConfigDir 'BitsSharpConfig.txt']. Presence of the file means
% to use a Bits# device, absence means to treat any device as a Bits+ device.
% Presence of a serial port device file name in the first line of that text
% configuration file will use the serial port device with that name for
% communication, otherwise the driver will try to auto-detect the proper
% serial port for communication. If you have multiple such devices, you can add
% multiple lines with serial port names to the configuration file, and then select
% a given device by specifying a index into the file as 'portSpec'. portSpec=1
% would use line 1 of the file for finding the device, portSpec=2 would use line 2
% etc.
%
%
% [handle, portHandle] = BitsPlusPlus('OpenBits#' [, portSpec]);
% -- Open a serial port control connection to a connected Bits# device, return
% a non-zero 'handle' to it on success, or the value zero if no such device exists
% or could not be opened. On successfull open, a low-level 'portHandle' is returned
% as well, which allows access to the serial port control connection via IOPort for
% people who really know what they are doing.
%
% Note: The current implementation only allows to connect to one device at a time,
% therefore use of the returned 'handle' is not very useful at the moment.
%
% The 'portSpec' parameter is optional and defines the name of the serial
% port(-device file) or index of the CRS device to use for the connection.
% If omitted, the name will be taken from a configuration file, or auto-detected.
% If specified as a string, the string defines the serial port device file. If
% specified as a numeric index, it selects the line within the configuration file
% from which the portname string will be taken, ie. index 1 = 1st line, 2 = 2nd line...
%
% This function must be called before use of any Bits# specific functions, otherwise
% they'll turn into no-ops or failures. This function can be called multiple times for
% a given device. It will only open the connection on first call. Successive calls will
% do nothing but increment a reference count of clients to the device.
%
%
% rc = BitsPlusPlus('Close' [, handle]);
% -- Decrement reference count to a Bits# device 'handle', close the serial
% connection to it once the count drops to zero, ie., as soon as nobody is using
% the connection anymore. If 'handle' is omitted, the 'Close' call is executed for
% each currently open device.
%
%
% rc = BitsPlusPlus('ResetOnWindowClose');
% -- Like 'Close', but switch display back to Bits++ video mode first, as
% that mode is "GUI friendly". Usually automatically called from Screen()
% when the Bits# stimulation onscreen window gets closed, at least if you
% used PsychImaging() to control the device.
%
%
% rc = BitsPlusPlus('CheckGPUSanity', window, xoffset [, injectFault=0]);
% -- Perform online-test of GPU identity gamma tables and DVI-D display
% encoders. Try to correct problems with wrong identity gamma tables and at
% least detect problems with (spatio-)temporal display dithering. Returns
% rc == 0 on full success, rc > 0 on failure.
% If the optional 'injectFault' parameter is set to 1, then an intentionally
% perturbed gamma table is loaded into the gpu to test how well the gamma table
% tweaking code is able to recover from wrong tables.
%
%
% pixels = BitsPlusPlus('GetVideoLine', nrPixels, scanline);
% -- Return the first (left-most) 'nrPixels' pixels in video scanline
% 'scanline' of the video display driven by a Bits# device. 'pixels' is
% a uint8 matrix with three rows for red, green and blue pixel color values,
% and 'nrPixels' columns, the three elements of each column encoding the
% r,g,b color values of the pixel corresponding to that column (x-position)
% of the scanline (y-position). Values are read back via the USB-Serial
% connection from the Bits# and the device sends back the pixel data as
% received over the DVI-D link.
%
%
% BitsPlusPlus('SwitchToBits++');
% -- Switch Bits# to Bits++ mode.
%
%
% BitsPlusPlus('SwitchToMono++');
% -- Switch Bits# to Mono++ mode.
%
%
% BitsPlusPlus('SwitchToColor++');
% -- Switch Bits# to Color++ mode.
%
%
% BitsPlusPlus('SwitchToStatusScreen');
% -- Switch Bits# to Status screen display.
%
%
% BitsPlusPlus('MassStorageMode');
% -- Switch Bits# to MassStorageMode. This will forcefully close allow
% client connections to the device, close the USB serial port connection
% and close the driver. The Bits# will report into USB mass storage mode,
% where it can get accessed like a USB flash drive, e.g., to edit configuration
% files, update firmware or EDID's etc. Only a power-cycle will bring the
% device back into a mode which allows us to connect to it again.
%
%

% History:
% 22.04.2007 Written (MK).
% xx.12.2007 Support for DIO T-Lock code generation (MK).
% 17.04.2008 Add support for overlay windows in Mono++ mode, and for color
%            correction/gamma correction via PsychColorCorrection (MK).
%  4.07.2009 Add support for other color correction methods like CLUT (MK).
% 14.12.2009 Add support for other target devices, e.g., DataPixx (MK).
%  3.01.2010 Some bugfixes to DataPixx support. (MK)
% 12.01.2013 Make compatible with PTB panelfitter. (MK)
% 13.03.2013 Make compatible with CRS Bits# video display system. (MK)
% 15.04.2013 Add mode = BitsPlusPlus('GetColorConversionMode', win); (MK)
% 24.06.2016 Improved diagnostics. Now need validation from BitsPlusImagingPipelineTest
%            and BitsPlusIdentityClutTest, instead of only BitsPlusImagingPipelineTest.
%            Now also triggers revalidation if OS kernel version/release changes, as
%            things like dithering or output color depth control are usually located in
%            the kernel level display drivers, so new kernel can imply new bugs.
% 01.05.2017 Replace evalin('base', ...) call by direct call in PsychDataPixx clut updates,
%            as a small performance optimization. (MK)
%            Now possible due to ability for Screen() to assign imaging pipeline variables
%            to callers workspace. This requires Octave 3.8+ or Matlab R2012a, which is our
%            supported minimum requirement for Psychtoolbox 3.0.14.
%
% 17.06.2018 Add support for dual-display setups with two Bits+/Bits#/Display++
%            Add support for returning device handles, port handles, and specifying
%            devices by number in BitsSharpConfig file. (MK)
global GL;

% Flag for validation: If not set to one, then this routine will check if
% proper operation of Bitsplusplus with GPU imaging has been verified.
persistent validated;

% Type of box: 0 = Bits+, 1 = Datapixx:
persistent targetdevicetype;

% Single output device (=0), or dual devices (=1), e.g., for dual-display:
persistent dualdevices;

% Name strings:
persistent devname;
persistent drivername;
persistent bplusname;
persistent mononame;
persistent colorname;
persistent devbits;
persistent checkGPUEncoders;

persistent bitsSharpPort;
persistent refCount;

% Vector that assigns overlay window handles to onscreen window handles:
persistent OverlayWindows;

% Encoded T-Lock display list handle for driving Bits++ DIO:
persistent tlockhandle;
% Counter of pending T-Lock display list blits: Zero == Disabled.
persistent blitTLockCode;
% Corrective x-offset for DIO blitting:
persistent tlockXOffset;

% Opmode for color conversion/buffer sampling in Color++ / C48 mode:
persistent colorConversionMode;

% Vector of cached per-window colorConversionMode:
persistent colorConversionModeWin;

if nargin < 1
    error('You must specify a command in argument "cmd"!');
end

win = [];
winRect = [];

if cmd == 1
    % Fast callback path for PTB imaging pipeline. We got called from the
    % finalizer blit chain of the imaging pipeline, asking us to perform
    % some post-processing on the final framebuffer image, immediately
    % before bufferswap.
    %
    % Currently, the only supported operation is drawing of a DIO T-Lock
    % code into the framebuffer, for control of the DIO pins of the Bits++
    % box. The T-Lock code has been generated already by a call to
    % 'DIOCommand'. We just have to blit that "Code Image" to the
    % framebuffer. We can't use Screen() commands here as we are called
    % from inside Screen -- Screen is not reentrant!
    %
    % We only blit if there is something to blit. Then we reset to nothing
    % to blit:
    if blitTLockCode ~= 0
        glCallList(tlockhandle);
        blitTLockCode = blitTLockCode - 1;
    end
    
    return;
end

% Default debuglevel for output during initialization:
debuglevel = 1;

if isempty(validated)
    validated = 0;
    tlockhandle = 0;
    blitTLockCode = 0;
    tlockXOffset = 0;
    OverlayWindows = [];
    targetdevicetype = 0;
    dualdevices = 0;
    drivername = 'BitsPlusPlus';
    devname = 'Bits+';
    bplusname = 'Bits++';
    mononame = 'Mono++';
    colorname = 'Color++';
    devbits = 14;
    checkGPUEncoders = 0;
    colorConversionMode = [];
    colorConversionModeWin = [];
    bitsSharpPort = [];
    refCount = 0;
end

if strcmpi(cmd, 'DIOCommand')

    % Setup new DIO command to be converted to T-Lock code and blitted:
    if nargin < 2 || isempty(arg)
        error('window handle for Bits++ onscreen window missing!');
    end
    
    if nargin < 3 || isempty(dummy)
        error('Number of repetitions for DIO command missing!');
    end

    if nargin < 6
        error('DIOCommand must have the parameters "Mask", "Command" and "Data"!');
    end

    mask = varargin{1};
    data = varargin{2};
    command = varargin{3};

    % Create or recreate our display list:
    glNewList(tlockhandle, GL.COMPILE);
    
    for i=1:size(mask, 1)
        % Process i'th row of command sequence:

        % Generate DIO T-Lock image as Matlab matrix:
        tlockdata = BitsPlusDIO2Matrix(mask(i), data(i,:), command(i));

        % Convert from Matlab matrix to OpenGL pixel format:
        encodedDIOdata = uint8(zeros(3, 508));
        % Pack 3 separate RGB planes into rows 1,2,3. As Matlabs data format is
        % column major order, this will end up as tightly packed pixel array in
        % format RGBRGBRGB.... just as glDrawPixels likes it.
        encodedDIOdata(1,:) = tlockdata(1,:,1);
        encodedDIOdata(2,:) = tlockdata(1,:,2);
        encodedDIOdata(3,:) = tlockdata(1,:,3);

        if nargin >= 7
            % Optional x, y blit position provided:
            xDIO = varargin{4};
            yDIO = varargin{5};
            xDIO = xDIO(i);
            yDIO = yDIO(i);
            
            if yDIO < 1
                yDIO = 1;
            end
        else
            % Set default position: 3rd scanline of display, so we don't get
            % into the way of a possible CLUT T-Lock code:
            xDIO = 0;
            yDIO = 2 + i;
        end
        
        % Add command sequence for this T-Lock code to display list:
        glRasterPos2i(xDIO + tlockXOffset, yDIO);
        glDrawPixels(508, 1, GL.RGB, GL.UNSIGNED_BYTE, encodedDIOdata);
    end
    
    % Finish display list;
    glEndList;
        
    % Assign number of repetitions:
    blitTLockCode = dummy;

    % Done. Return.
    return;

end

if strcmpi(cmd, 'DIOCommandReset')
    % Dummy error check: arg will be used in later revisions...
    if nargin < 2 || isempty(arg)
        error('window handle for Bits++ onscreen window missing!');
    end
    
    % Disable T-Lock blitting:
    blitTLockCode = 0;
    return;
end

if strcmpi(cmd, 'SetTargetDeviceType')
    if nargin < 2
        error('targetdevicetype parameter missing!');
    end
    
    % Assign targetdevicetype to internal persistent variable:
    % A zero means: It is a regular CRS Bits+ box.
    % A 1 means: It is a VPixx DataPixx box which is sharing setup code
    % with Bits+ in this file:
    targetdevicetype = arg;
    
    switch (targetdevicetype)
        case 0,
            drivername = 'BitsPlusPlus';
            if ~isempty(bitsSharpPort)
                devname = 'Bits#';
            else
                devname = 'Bits+';
            end
            bplusname = 'Bits++';
            mononame = 'Mono++';
            colorname = 'Color++';
            devbits = 14;
            
        case 1,
            drivername = 'PsychDatapixx';
            devname = 'DataPixx';
            bplusname = 'L48';
            mononame = 'M16';
            colorname = 'C48';
            devbits = 16;
            
        otherwise
            error('Unknown targetdevicetype assigned in call to "SetTargetDeviceType"!');
    end
    
    return;
end

if strcmpi(cmd, 'SetDualDevices')
    if nargin < 2
        error('dualdevices parameter missing!');
    end

    % Assign dualdevices flag to internal persistent variable:
    dualdevices = arg;

    return;
end

if strcmpi(cmd, 'OpenBits#')
    % Try to open connection to a Bits# device. Return true if successfull,
    % false otherwise - which would likely imply a Bits+ instead of Bits#.
    winRect = [];

    if ~isempty(bitsSharpPort)
        % Already open. Do nothing but return and report success:
        win = 1;
        winRect = bitsSharpPort;

        % Increment reference count:
        refCount = refCount + 1;

        return;
    end

    % Explicit serial port name for connection to device provided?
    if nargin > 1 && ~isempty(arg)
        bitsSharpPortname = arg;
        if ~ischar(bitsSharpPortname) && (~isnumeric(bitsSharpPortname) || ~isscalar(bitsSharpPortname))
            error('Provided Bits# serial port name is not a valid namestring or device index!');
        end
    else
        % No portname given:
        bitsSharpPortname = [];
    end

    % Have a portname?
    if isempty(bitsSharpPortname) || isnumeric(bitsSharpPortname)
        % No: Find out if a Bits# configuration file exists. Otherwise we assume
        % that usercode does not want to connect to a Bits# but user probably uses
        % an older - connectionless - Bits+, and therefore turn ourselves into a no-op:
        configfile = [PsychtoolboxConfigDir 'BitsSharpConfig.txt'];
        if ~exist(configfile, 'file')
            % No config file -> No Bits#. We no-op and return "no such device":
            win = 0;
            fprintf('BitsPlusPlus: Could not find a Bits# config file under [%s]. Assuming a Bits+ device instead of a Bits# is connected.\n', configfile);
            fprintf('BitsPlusPlus: Please create a config file under this name if you have a Bits# and want to use it as Bits# instead of as a Bits+.\n');
            fprintf('BitsPlusPlus: The most simple way is to create an empty file. A more robust way is to store the name of the Bits# serial port\n');
            fprintf('BitsPlusPlus: in the first line of the text file, e.g., COM5 [Windows], or /dev/ttyACM0 [Linux] or similar.\n');
            return;
        end

        % File exists -> We want to access a Bits#. Parse file for a port name string:
        if isempty(bitsSharpPortname)
            bitsSharpPortname = 1;
        end

        fileContentsWrapped = [];
        fid = fopen(configfile);
        % Get specified line:
        for lineidx=1:bitsSharpPortname
            fileContentsWrapped = fgets(fid);
        end
        fclose(fid);

        % Port spec available?
        if ~isempty(fileContentsWrapped) && ischar(fileContentsWrapped)
            % Yes: Assign namestring for port.
            bitsSharpPortname = deblank(fileContentsWrapped);
            fprintf('BitsPlusPlus: Connecting to Bits# device via serial port [%s], as provided by configuration file [%s].\n', bitsSharpPortname, configfile);
        else
            % No: Do the guess-o-matic dance: Fail softly if it doesn't work:
            try
                % Try to find proper serial port:
                bitsSharpPortname = FindSerialPort([], 1, 0);
                fprintf('BitsPlusPlus: Connecting to Bits# device via auto-detected serial port [%s].\n', bitsSharpPortname);
            catch
                lerr = psychlasterror('reset');
                disp(lerr.message);
                fprintf('BitsPlusPlus: Failed to find the Bits# device! Is it connected and ready? See diagnostics above. Continuing without Bits# support.\n');
                win = 0;
                return;
            end
        end
    else
        fprintf('BitsPlusPlus: Connecting to Bits# device via serial port [%s], as provided by usercode.\n', bitsSharpPortname);
    end

    % Ok. Try to connect:
    try
        % Open the port. We reduce level of verbosity during open, to suppress clutter:
        oldverblevel = IOPort('Verbosity', 1);
        [bitsSharpPort, errmsg] = IOPort('OpenSerialPort', bitsSharpPortname);
        IOPort('Verbosity', oldverblevel);
    catch %#ok<*CTCH>
        error('Failed to establish a connection to the Bits# via serial port.');
    end

    % Success?
    if bitsSharpPort < 0
        % No: Delete invalid port handle.
        bitsSharpPort = []; %#ok<NASGU>
        error('Failed to establish a connection to the Bits# via serial port. The error message was: %s', errmsg);
    end

    % Yes. Change reported device name:
    devname = 'Bits#';

    fprintf('BitsPlusPlus: Device information:\n');

    % Some queries to the device to test the connection:
    IOPort('Write', bitsSharpPort, ['$ProductType' char(13)]);
    WaitSecs('YieldSecs', 0.1);
    while IOPort('BytesAvailable', bitsSharpPort)
        fprintf('BitsPlusPlus: %s\n', deblank(char(IOPort('Read', bitsSharpPort))));
        WaitSecs('YieldSecs', 0.1);
    end

    IOPort('Write', bitsSharpPort, ['$SerialNumber' char(13)]);
    WaitSecs('YieldSecs', 0.1);
    while IOPort('BytesAvailable', bitsSharpPort)
        fprintf('BitsPlusPlus: %s\n', deblank(char(IOPort('Read', bitsSharpPort))));
        WaitSecs('YieldSecs', 0.1);
    end

    IOPort('Write', bitsSharpPort, ['$FirmwareDate' char(13)]);
    WaitSecs('YieldSecs', 0.1);
    while IOPort('BytesAvailable', bitsSharpPort)
        fprintf('BitsPlusPlus: %s\n', deblank(char(IOPort('Read', bitsSharpPort))));
        WaitSecs('YieldSecs', 0.1);
    end

    fprintf('\n');

    % Increment reference count:
    refCount = refCount + 1;

    % Report success:
    win = 1;

    % Return IOPort handle as well:
    winRect = bitsSharpPort;

    return;
end

if strcmp(cmd, 'GetVideoLine')
    % Readback a video line of pixels from the device:
    %
    % Usage:
    %
    % BitsPlusPlus('GetVideoLine', nrPixels, scanline);

    % Parse inputs: Number of pixels to read back and scanline index
    % of scanline to read back:
    nrPixels = arg;
    scanline = dummy;

    % Perform pixel readback from device:
    scanline = BitsSharpGetScanline(bitsSharpPort, scanline, nrPixels);
    
    % Return scanline pixels as 1st return argument:
    win = scanline;

    return;
end

if strcmpi(cmd, 'PerformPostWindowOpenSetup')
    % Called from PsychImaging after onscreen window and associated imaging pipeline
    % is fully opened and initialized:

    % Get windowhandle of associated onscreen window:
    win = arg;

    % Attach a window close callback for Device teardown at window close time:
    Screen('Hookfunction', win, 'AppendMFunction', 'CloseOnscreenWindowPostGLShutdown', 'Shutdown window callback into BitsPlusPlus driver.', 'BitsPlusPlus(''ResetOnWindowClose'');');
    Screen('HookFunction', win, 'Enable', 'CloseOnscreenWindowPostGLShutdown');

    return;
end

if strcmpi(cmd, 'ResetOnWindowClose')
    % Called from Screen() at onscreen window close time, or manually from usercode:

    % Explicit device handle specified?
    if nargin >= 2 && ~isempty(arg)
        % Yes. Currently we only accept the handle 1, as we can not do
        % multi-device yet. Just validate for future correctness:
        if arg ~= 1 || isempty(bitsSharpPort)
            error('BitsPlusPlus: Close: Invalid ''handle'' specified. No such device open.');
        end
    end

    % Connection to Bits# established?
    if ~isempty(bitsSharpPort)
        % Yes. Switch back to Bits++ display mode, which provides a reasonably
        % useable display of the regular desktop GUI:
        fprintf('BitsPlusPlus: Switching Bits# device to desktop GUI friendly Bits++ mode.\n');
        IOPort('Write', bitsSharpPort, ['$BitsPlusPlus' char(13)]);

        % Flush commands:
        IOPort('Flush', bitsSharpPort);

        % Signal that device connection was actually open:
        win = 1;
    else
        % Signal that device connection was not actually open:
        win = 0;
    end

    % Fall through to close command:
    cmd = 'Close';
end

if strcmpi(cmd, 'Close')
    % Connection to Bits# established?

    % Explicit device handle specified?
    if nargin >= 2 && ~isempty(arg)
        % Yes. Currently we only accept the handle 1, as we can not do
        % multi-device yet. Just validate for future correctness:
        if arg ~= 1 || isempty(bitsSharpPort)
            error('BitsPlusPlus: Close: Invalid ''handle'' specified. No such device open.');
        end
    end

    % More than one client (this calling client) holding a reference to Bits# ?
    if refCount > 1
        % Yes. Just decrement the refCount to release this reference and be done:
        refCount = refCount - 1;
        win = 1;
        return;
    end
    
    % Caller is last client of Bits#. Reset refCount to zero and really
    % close the connection:
    refCount = 0;
    
    if ~isempty(bitsSharpPort)
        % Yes. Close connection:
        IOPort('Close', bitsSharpPort);
        bitsSharpPort = [];
        fprintf('BitsPlusPlus: Connection to Bits# device closed.\n');

        % Signal that device connection was actually open:
        win = 1;
    else
        % Signal that device connection was not actually open:
        win = 0;
    end

    return;
end

if strcmpi(cmd, 'SwitchToBits++')
    % Connection to Bits# established?
    if ~isempty(bitsSharpPort)
        % Yes. Switch back to Bits++ display mode, which provides a reasonably
        % useable display of the regular desktop GUI:
        fprintf('BitsPlusPlus: Switching Bits# device to Bits++ mode.\n');
        IOPort('Write', bitsSharpPort, ['$BitsPlusPlus' char(13)]);

        % Flush commands:
        IOPort('Flush', bitsSharpPort);

        win = 1;
    else
        win = 0;
    end

    return;
end

if strcmpi(cmd, 'SwitchToMono++')
    
    % Connection to Bits# established?
    if ~isempty(bitsSharpPort)
        % Yes. Switch back to Bits++ display mode, which provides a reasonably
        % useable display of the regular desktop GUI:
        fprintf('BitsPlusPlus: Switching Bits# device to Mono++ mode.\n');
        IOPort('Write', bitsSharpPort, ['$monoPlusPlus' char(13)]);

        % Flush commands:
        IOPort('Flush', bitsSharpPort);

        win = 1;
    else
        win = 0;
    end

    return;
end

if strcmpi(cmd, 'SwitchToColor++')
    
    % Connection to Bits# established?
    if ~isempty(bitsSharpPort)
        % Yes. Switch back to Bits++ display mode, which provides a reasonably
        % useable display of the regular desktop GUI:
        fprintf('BitsPlusPlus: Switching Bits# device to Color++ mode.\n');
        IOPort('Write', bitsSharpPort, ['$colourPlusPlus' char(13)]);

        % Flush commands:
        IOPort('Flush', bitsSharpPort);

        win = 1;
    else
        win = 0;
    end

    return;
end

if strcmpi(cmd, 'SwitchToStatusScreen')
    
    % Connection to Bits# established?
    if ~isempty(bitsSharpPort)
        % Yes. Switch back to Bits++ display mode, which provides a reasonably
        % useable display of the regular desktop GUI:
        fprintf('BitsPlusPlus: Switching Bits# device to Diagnostic display mode.\n');
        IOPort('Write', bitsSharpPort, ['$statusScreen' char(13)]);

        % Flush commands:
        IOPort('Flush', bitsSharpPort);

        win = 1;
    else
        win = 0;
    end

    return;
end

if strcmpi(cmd, 'MassStorageMode')
    
    % Connection to Bits# established?
    if ~isempty(bitsSharpPort)
        % Yes. Switch back to Bits++ display mode, which provides a reasonably
        % useable display of the regular desktop GUI:
        fprintf('BitsPlusPlus: Switching Bits# device to MassStorageMode mode and disconnecting...\n');
        IOPort('Write', bitsSharpPort, ['$USB_massStorage' char(13)]);

        % Flush commands:
        IOPort('Flush', bitsSharpPort);

        % Disconnect forcefully, as the Bits# will no longer respond to our
        % commands until it is manually restarted:
        refCount = 0;
        IOPort('Close', bitsSharpPort);
        bitsSharpPort = [];
        
        win = 1;
    else
        win = 0;
    end

    return;
end

if strcmpi(cmd, 'SetColorConversionMode')
    % Set the mode of operation for color conversion in Color++ / C48 mode.
    % This is a mandatory call in that mode. As the effective output
    % resolution is only half the framebuffer resolution we need to decide
    % what tradeoff between aspect-ratio preservation, sampling precision etc.
    % to take.
    colorConversionMode = arg;
    return;
end

if strcmpi(cmd, 'GetColorConversionMode')
    % Return the mode of operation for color conversion in Color++ / C48 mode.
    if nargin < 2 || isempty(arg) || ~isa(arg, 'double') || (Screen('WindowKind', arg) ~= 1)
        error('%s: "GetColorConversionMode" called without valid onscreen window handle.', drivername);
    end

    win = colorConversionModeWin(arg);
    return;
end

if strcmpi(cmd, 'TestGPUEncoders')
    % Perform check of GPU identity gamma tables and encoders during next
    % 'OpenWindowXXX' call in Datapixx mode or Bits# mode. This is a one-shot,
    % auto-reset flag:
    checkGPUEncoders = 1;
    return;
end

if strcmpi(cmd, 'ForceUnvalidatedRun')
    % Enforce use of this routine without verification of correct function
    % of the imaging pipeline. This is used by the correctness test itself
    % in order to be able to run the validation.
    validated = 1;
    return;
end

if strcmpi(cmd, 'StoreValidation')
    % Enforce use of this routine without verification of correct function
    % of the imaging pipeline. This is used by the correctness test itself
    % in order to be able to run the validation.
    validationType = dummy;
    ValidateBitsPlusImaging(arg, 1, devname, validationType);
    return;
end

if strcmpi(cmd, 'LoadIdentityClut')
    % Load an identity CLUT into Bits++ at next Screen('Flip'). This is
    % just a little convenience wrapper around 'LoadNormalizedGammaTable':
    % Restore Bits++ Identity CLUT so it can be used as normal display:
    if nargin < 2 || isempty(arg)
        error('window handle for %s onscreen window missing!', devname);
    end
    linear_lut =  repmat(linspace(0, 1, 256)', 1, 3);
    Screen('LoadNormalizedGammaTable', arg, linear_lut, 2);
    return;
end

if strcmpi(cmd, 'OpenWindowBits++')
    % Execute the Screen('OpenWindow') command with proper flags, followed
    % by our own Initialization. Return values of 'OpenWindow'.
    %
    % This will set up the Bits++ mode of Bits++

    % Assign screen index:
    if nargin < 2 || isempty(arg) || ~isa(arg, 'double')
        error('%s: "OpenWindow..." called without valid screen handle.', drivername);
    end
    screenid = arg;

    % Assign optional clear color:
    if nargin < 3
        clearcolor = [];
    else
        clearcolor = dummy;
    end

    % windowRect is always full screen -- Anything else would make the
    % Bits++ display fail.
    if IsLinux && ~isempty(varargin)
        winRect = varargin{1};
    else
        winRect = [];
    end

    % pixelSize is also fixed to 32 bit RGBA8 framebuffer:
    pixelSize = 32;

    % Same for double-buffering:
    numbuffers = 2;

    % stereomode we take...
    if nargin >= 7
        stereomode = varargin{4};
    else
        stereomode = [];
    end

    % multiSample gets forced to zero, as it would interfere
    % with Bits++ display controller:
    multiSample = 0;

    % Open the window, pass all parameters (partially modified or overriden), return Screen's return values:
    if nargin >= 9
        [win, winRect] = Screen('OpenWindow', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, varargin{6:end});
    else
        [win, winRect] = Screen('OpenWindow', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample);
    end

    % Ok, if we reach this point then we've got a proper onscreen
    % window on the Bits++. Let's reassign our arguments and continue with
    % the init sequence:

    % First load the graphics hardwares gamma table with an identity mapping,
    % so it doesn't interfere with Bits++ -- Functions from Bits++ toolbox.
    LoadIdentityClut(win);

    % We need the GL for DIO T-Lock setup:
    if isempty(GL)
        % Load & Initalize constants and moglcore, but don't set the 3D gfx
        % flag for Screen():
        InitializeMatlabOpenGL([], [], 1);
    end;

    % Test accuracy/correctness of GPU's rasterizer for different output
    % positioning methods: Return (non-zero) dx,dy offsets, if any:
    [rpfx, rpfy, rpix, rpiy, vix] = PsychGPURasterizerOffsets(win, drivername); %#ok<ASGLU>
        
    if rpix~=0
        tlockXOffset = -rpix;
        fprintf('OpenWindow%s: Applying corrective horizontal DIO T-Lock offset of %i pixels for buggy graphics card driver. Will hopefully fix it...\n', bplusname, tlockXOffset);        
    end

    if targetdevicetype == 1 && checkGPUEncoders
        % Perform DataPixx builtin diagnostics to detect problems with
        % wrong GPU gamma tables or GPU dithering:
        checkGPUEncoders = 0;
        if PsychDataPixx('CheckGPUSanity', win, tlockXOffset)
            % Ohoh, trouble ahead! The driver detected problems with the
            % GPU and wasn't able to auto-correct them.
            fprintf('%s: CAUTION! DataPixx internal diagnostic detected problems with your graphics card driver which it could not correct by itself!\n', drivername);
        end
    end

    if targetdevicetype == 0 && checkGPUEncoders
        % Perform Bits# builtin diagnostics to detect problems with
        % wrong GPU gamma tables or GPU dithering:
        checkGPUEncoders = 0;
        if doCheckGPUSanity(win, tlockXOffset, bitsSharpPort)
            % Ohoh, trouble ahead! The driver detected problems with the
            % GPU and wasn't able to auto-correct them.
            fprintf('%s: CAUTION! Bits# internal diagnostic detected problems with your graphics card driver which it could not correct by itself!\n', drivername);
        end
    end

    % Connection to Bits# established?
    if ~isempty(bitsSharpPort)
        % Yes. Switch it to Bits++ display mode:
        IOPort('Write', bitsSharpPort, ['$BitsPlusPlus' char(13)]);
        fprintf('BitsPlusPlus: Switching Bits# device to Bits++ video mode. This will take about 5 seconds...\n');
        WaitSecs('YieldSecs', 5);
    end
    
    % Now enable finalizer hook chains and load them with the special Bits++
    % command for T-Lock based Bits++ internal CLUT updates:
    if vix~=0
        % vix is wrong offset, therefore negate it to get corrective offset:
        vix = -vix;
        fprintf('OpenWindow%s: Applying corrective horizontal offset of %i pixels for buggy graphics card driver. Will hopefully fix it...\n', bplusname, vix);
        offsetstring = sprintf('xPosition=%i', vix);
    else
        offsetstring = '';
    end
    
    if targetdevicetype == 0
        Screen('HookFunction', win, 'PrependBuiltin', 'LeftFinalizerBlitChain', 'Builtin:RenderClutBits++', offsetstring);
    end
    
    if targetdevicetype == 1
        rclutcmd = 'PsychDataPixx(1, IMAGINGPIPE_GAMMATABLE);';
        Screen('HookFunction', win, 'PrependMFunction', 'LeftFinalizerBlitChain', 'Upload new clut into DataPixx callback', rclutcmd);
    end
    
    Screen('HookFunction', win, 'Enable', 'LeftFinalizerBlitChain');

    if ~isempty(stereomode) && ismember(stereomode, [1, 4, 5, 10, 11])
        % Enable CLUT updates on right stereo buffer as well:

        % CRS devices, frame-sequential stereo or native quad-buffered OpenGL stereo,
        % or dual-window stereo (mostly for OSX dual-display):
        if targetdevicetype == 0 && ismember(stereomode, [1, 10, 11])
            Screen('HookFunction', win, 'PrependBuiltin', 'RightFinalizerBlitChain', 'Builtin:RenderClutBits++', offsetstring);
        end

        % Dual CRS devices, dual-display stereo on Linux or Windows:
        if targetdevicetype == 0 && ismember(stereomode, [4, 5]) && dualdevices == 1
            % Here one PTB window drives both video outputs (split-window style) and
            % thereby both connected CRS devices. In order to set a CLUT on both of
            % these devices, we need to blit the Clut T-Lock into the right half of
            % our window as well. However, we use the LeftFinalizerBlitChain for this,
            % as the RightFinalizerBlitChain is not operational in stereomode 4 and 5:
            offsetstring = sprintf('xPosition=%i', vix + Screen('WindowSize', win, 1) / 2);
            Screen('HookFunction', win, 'PrependBuiltin', 'LeftFinalizerBlitChain', 'Builtin:RenderClutBits++', offsetstring);
        end

        % VPixx devices, frame-sequential stereo or native quad-buffered OpenGL stereo:
        if targetdevicetype == 1 && ismember(stereomode, [1, 11])
            Screen('HookFunction', win, 'PrependMFunction', 'RightFinalizerBlitChain', 'Upload new clut into DataPixx callback', rclutcmd);
        end

        Screen('HookFunction', win, 'Enable', 'RightFinalizerBlitChain');
    end

    if targetdevicetype == 0
        % Setup finalizer callback for DIO T-Lock updates:
        tlockhandle = SetupDIOFinalizer(win, stereomode);
    end

    % Load an identity CLUT into the Bits++ to start with:
    linear_lut =  repmat(linspace(0, 1, 256)', 1, 3);
    Screen('LoadNormalizedGammaTable', win, linear_lut, 2);
    
    % Check validation:
    if ~validated
        % Only validate output display pipeline from framebuffer to display device:
        ValidateBitsPlusImaging(win, 0, devname, 2);
    end

    % Reset validation flag after first run:
    validated = 0;

    % Set colorConversionMode for this window to safe "undefined" default:
    colorConversionModeWin(win) = -1;

    % Ready!
    return;
end

if strcmpi(cmd, 'OpenWindowMono++') || strcmpi(cmd, 'OpenWindowMono++WithOverlay') || strcmpi(cmd, 'OpenWindowColor++')
    % Execute the Screen('OpenWindow') command with proper flags, followed
    % by our own Initialization. Return values of 'OpenWindow'.
    %
    % This will set up the Mono++ or Color++ mode of Bits++
    
    % Assign screen index:
    if nargin < 2 || isempty(arg) || ~isa(arg, 'double')
        error('%s: "OpenWindow..." called without valid screen handle.', drivername);
    end
    screenid = arg;

    % Assign optional clear color:
    if nargin < 3
        clearcolor = [];
    else
        clearcolor = dummy;
    end

    if isempty(clearcolor)
        clearcolor = 1.0;
    end
    
    % windowRect is always full screen -- Anything else would make the
    % Bits++ display fail.
    if IsLinux && ~isempty(varargin)
        winRect = varargin{1};
    else
        winRect = [];
    end

    % pixelSize is also fixed to 32 bit RGBA8 framebuffer:
    pixelSize = 32;

    % Same for double-buffering:
    numbuffers = 2;

    % stereomode we take...
    if nargin >= 7
        stereomode = varargin{4};
    else
        stereomode = [];
    end

    % Retrieve multiSample setting:
    if nargin >= 8
        multiSample = varargin{5};
    else
        multiSample = [];
    end

    % Imaging mode we take - and combine it with our own requirements:
    if nargin >= 9
        imagingmode = varargin{6};
    else
        imagingmode = 0;
    end

    % For imagingmode we need at least fast backingstore, the output
    % formatter enabled and some high precision color buffer. We default to
    % 32 bpc floating point, the only safe choice accross different graphics cards
    % from different vendors, but the users imagingmode setting is free to
    % override this with a 16 bpc fixed buffer. 16 bpc float works as well
    % but can't use the full Bits++ color range at full precision.
    if bitand(imagingmode, kPsychNeed16BPCFloat) || bitand(imagingmode, kPsychNeed16BPCFixed) || bitand(imagingmode, kPsychUse32BPCFloatAsap)
        % User specified override: Use it.
        ourspec = 0;
    else
        % No user specified accuracy. We play safe and choose the highest
        % one:
        ourspec = kPsychNeed32BPCFloat;
    end

    % Imagingmode must at least include the following:
    imagingmode = mor(imagingmode, kPsychNeedFastBackingStore, kPsychNeedOutputConversion, ourspec);

    if strcmpi(cmd, 'OpenWindowColor++')
        if isempty(colorConversionMode)
            sca;
            fprintf('The new mandatory parameter "colorConversionMode" is missing!\n');
            fprintf('If you used BitsPlusPlus(''OpenWindowColor++'', ...); to get here, then\n');
            fprintf('add the command BitsPlusPlus(''SetColorConversionMode'', mode);\n');
            fprintf('immediately before the BitsPlusPlus(''OpenWindowColor++'', ...); command.\n\n');
            fprintf('If you used the more modern and recommended PsychImaging() commands to get here, then\n');
            fprintf('change your call to PsychImaging(''AddTask'', ''General'', ''EnableBits++Color++Output'');\n');
            fprintf('or to PsychImaging(''AddTask'', ''General'', ''EnableDataPixxC48Output''); into a call to \n');
            fprintf('PsychImaging(''AddTask'', ''General'', ''EnableBits++Color++Output'', mode);\n');
            fprintf('or PsychImaging(''AddTask'', ''General'', ''EnableDataPixxC48Output'', mode);\n\n');
            fprintf('The new parameter "mode" must be 0 if you want exactly the old behaviour back.\n');
            fprintf('For new code, you will likely want to use a value of 1 or 2 to preserve correct\n');
            fprintf('aspect ratio.\n\n');
            fprintf('Please read the help section for the PsychImaging() command ("help PsychImaging")\n');
            fprintf('for the ''EnableBits++Color++Output'' subcommand. It explains the meaning of the different\n');
            fprintf('possible settings of "mode".\n\n');
            
            error('Mandatory parameter "colorConversionMode" is missing!');
        end
        
        if ~ismember(colorConversionMode, [0,1,2]);
            sca;
            fprintf('The provided "colorConversionMode" parameter %i is not one of the valid values 0, 1 or 2!\n', colorConversionMode);
            error('Mandatory parameter "colorConversionMode" is invalid!');
        end
        
        if colorConversionMode == 0
            % In Color++ mode with "classic" conversion, we only have half the
            % effective horizontal resolution. Tell PTB to take this into
            % account for all relevant calculations:
            imagingmode = mor(imagingmode, kPsychNeedHalfWidthWindow);
        end        

        targetMode = ['$colourPlusPlus' char(13)];
        targetModeName = 'Color++';
    else
        % Mono++ mode, with or without overlay.
        targetMode = ['$monoPlusPlus' char(13)];
        targetModeName = 'Mono++';
    end

    % Open the window, pass all parameters (partially modified or overriden), return Screen's return values:
    % Note that we clear to black (==0), because we set the real background
    % clear color "further down the road" after we've established our
    % default color range of 0.0 - 1.0, ie. in the normalized 0 - 1 range.
    if nargin > 9
        [win, winRect] = Screen('OpenWindow', screenid, 0, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingmode, varargin{7:end});
    else
        [win, winRect] = Screen('OpenWindow', screenid, 0, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingmode);
    end

    % Ok, if we reach this point then we've got a proper onscreen
    % window on the Bits++. Let's reassign our arguments and continue with
    % the init sequence:

    % Some more diagnostics and info for user:
    winfo = Screen('GetWindowInfo', win);
    
    % Unconditional support for 32 bpc float drawable requested?
    havespoken = 0;
    if ~bitand(imagingmode, kPsychNeed32BPCFloat)
        % Nope. Conditional support requested?
        if (bitand(imagingmode, kPsychUse32BPCFloatAsap) && winfo.GLSupportsBlendingUpToBpc < 32)
            % Conditional use of 32 bpc float buffers requested, but GPU
            % doesn't support 32 bpc float blending --> drawBuffer will only be
            % 16 bpc -- Loss of precision!
            fprintf('PTB - Info: Your framebuffer is only configured to provide about 10-11 bits of precision, because your\n');
            fprintf('PTB - Info: script requested support for simultaneous alpha-blending and high precision, but your hardware is not\n');
            fprintf('PTB - Info: capable of supporting highest precision with alpha-blending enabled. You will therefore only\n');
            fprintf('PTB - Info: be able to use about 11 bits out of the %i bits precision that %s provides for stimulus definition.\n', devbits, devname);
            fprintf('PTB - Info: Stimulus postprocessing, e.g., gamma correction, will still make good use of all %i bits though.\n', devbits);
            fprintf('PTB - Info: You can either live with this limitation, or do not use alpha-blending or upgrade your graphics\n');
            fprintf('PTB - Info: hardware to Direct3D-10 compliant hardware, e.g., ATI Radeon HD-3000 or NVidia Geforce-8000 and later.\n\n');
            havespoken = 1;
        end

        if bitand(imagingmode, kPsychNeed16BPCFloat)
            fprintf('PTB - Info: Your framebuffer is only configured to provide about 10-11 bits of precision, because your\n');
            fprintf('PTB - Info: script requested only 16 bpc float precision. You will therefore only be able to use\n');
            fprintf('PTB - Info: about 11 bits out of the %i bits precision that %s provides for stimulus drawing.\n', devbits, devname);
            fprintf('PTB - Info: If you want to use the full %i bit precision, you will need to request a 32 bpc float framebuffer.\n\n', devbits);
            havespoken = 1;
        end

        if bitand(imagingmode, kPsychNeed16BPCFixed)
            fprintf('PTB - Info: Your framebuffer is configured to provide 16 bits of precision, because your\n');
            fprintf('PTB - Info: script requested 16 bits fixed precision. %s will be able to finally output %i bits precision.\n', devname, devbits);
            fprintf('PTB - Info: Alpha-blending will not work at this configuration with your hardware though. Choose a different\n');            
            fprintf('PTB - Info: mode if you need alpha-blending and high precision.\n\n');
            havespoken = 1;
        end
    end

    if (havespoken == 0) && (bitand(imagingmode, kPsychNeed32BPCFloat) || bitand(imagingmode, kPsychUse32BPCFloatAsap))
        fprintf('PTB - Info: Your framebuffer is configured for maximum precision. All internal processing will be done\n');
        fprintf('PTB - Info: with about 23 bits of linear precision -- %s will be able to finally output with %i bits precision.\n', devname, devbits);
        if winfo.GLSupportsBlendingUpToBpc < 32
            fprintf('PTB - Info: Alpha-blending will not work at this precision with your hardware though.\n');
            fprintf('PTB - Info: You can either live with this limitation, or upgrade your graphics hardware to Direct3D-10\n');
            fprintf('PTB - Info: compliant hardware, e.g., ATI Radeon HD-3000 or NVidia Geforce-8000 and later.\n\n');
        else
            fprintf('PTB - Info: Alpha-blending should be fully supported at this precision by your hardware.\n\n');            
        end
    end
    
    if strcmpi(cmd, 'OpenWindowColor++')
        if colorConversionMode == 0
            fprintf('PTB - Info: Classic half horizontal resolution color conversion for %s mode selected.\n', colorname);
            fprintf('PTB - Info: Aspect ratio will be horizontally distorted, ie., 2:1.\n');
        end

        if colorConversionMode == 1
            fprintf('PTB - Info: Aspect ratio preserving half horizontal resolution color conversion for %s\n', colorname);
            fprintf('PTB - Info: mode selected. All odd-numbered pixel columns will be ignored/skipped.\n');
        end

        if colorConversionMode == 2
            fprintf('PTB - Info: Aspect ratio preserving bilinear color conversion for %s mode selected.\n', colorname);
            fprintf('PTB - Info: Will average color between adjacent even/odd pixel columns.\n');
        end
        
        fprintf('\n');
    end
    
    % First load the graphics hardwares gamma table with an identity mapping,
    % so it doesn't interfere with Bits++ -- Function from Bits++ toolbox.
    LoadIdentityClut(win);

    % Backup current gfx-settings, so we can restore them after
    % modifications: The LoadGLSLProgramFromFiles() routine enables this
    % implicitely. This is unwanted in case we are in pure 2D mode, so we
    % need to undo it below...
    ogl = Screen('Preference', 'Enable3DGraphics');

    % Create and retrieve a compiled shader and idString-Snippet for
    % use with the formatting shader to allow for final
    % color-transformations immediately before Mono++ conversion. This
    % is mostly meant to implement gammacorrection, clamping or other
    % transformations needed for a well calibrated display:
    [icmShaders, icmIdString, icmConfig] = PsychColorCorrection('GetCompiledShaders', win, debuglevel);
    
    % Operate in Mono++ mode or Color++ mode?
    if strcmpi(cmd, 'OpenWindowMono++') || strcmpi(cmd, 'OpenWindowMono++WithOverlay')
        % Setup for Mono++ mode:
        
        if strcmpi(cmd, 'OpenWindowMono++WithOverlay')
            useOverlay = 1;
        else
            useOverlay = 0;
        end
        
        % Use of overlay plane requested?
        if useOverlay
            % Create additional shader for overlay texel fetch:
            % Our gpu panel scaler might be active, so the size of the
            % virtual window - and thereby our overlay window - can be
            % different from the output framebuffer size. As the sampling
            % 'pos'ition for the overlay is always provided in framebuffer
            % coordinates, we need to subsample in the overlay fetch.
            % Calculate proper scaling factor, based on virtual and real
            % framebuffer size:
            [wC, hC] = Screen('WindowSize', win);
            [wF, hF] = Screen('WindowSize', win, 1);
            sampleX = wC / wF;
            sampleY = hC / hF;
            
            % Build the shader:
            shSrc = sprintf('uniform sampler2DRect overlayImage; float getMonoOverlayIndex(vec2 pos) { return(texture2DRect(overlayImage, pos * vec2(%f, %f)).r); }', sampleX, sampleY);

            % Create Offscreen window for the overlay. It has the same size
            % as the onscreen window, but only 8 bpc fixed depth and a
            % completely black background -- fully transparent by default.
            % The specialflags 32 setting protects the overlay offscreen
            % window from accidental batch-deletion by usercode calls to
            % Screen('Close'):
            overlaywin = Screen('OpenOffscreenWindow', win, 0, [], 8, 32);

            % Retrieve low-level OpenGl texture handle to the window:
            overlaytex = Screen('GetOpenGLTexture', win, overlaywin);
            
            % Disable bilinear filtering on this texture - always use
            % nearest neighbour sampling to avoid interpolation artifacts
            % in color index image for clut indexing:
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, overlaytex);
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
        else
            % No.: Create "no-op" shader for zero overlay:
            shSrc = 'float getMonoOverlayIndex(vec2 pos) { return(0.0); }';
        end

        % Build shader from source:
        overlayShader = glCreateShader(GL.FRAGMENT_SHADER);
        glShaderSource(overlayShader, shSrc);
        glCompileShader(overlayShader);

        % Attach to list of shaders:
        icmShaders(end+1) = overlayShader;

        % Load Bits++ Mono++ formatting shader:
        shader = LoadGLSLProgramFromFiles('Bits++_Mono++_FormattingShader', debuglevel, icmShaders);

        if useOverlay
            % Ok, overlay requested. Setup shader's overlayImage sampler to
            % texture unit 1 and setup proper pString, so unit 1 has
            % overlay bound during blit operation:
            pString = sprintf('TEXTURERECT2D(1)=%i', overlaytex);
            
            glUseProgram(shader);
            glUniform1i(glGetUniformLocation(shader, 'overlayImage'), 1);
            glUseProgram(0);
            
            % Store window handle of overlay window for this onscreen
            % window for later retrieval:
            OverlayWindows(win) = overlaywin;
        else
            pString = '';
        end

        % Now enable output formatter hook chain and load them with the special Bits++
        % Mono++ data formatting shader: We append the shader because it
        % absolutely must be the last shader to execute in that chain!
        idString = sprintf('Mono++ output formatting shader for CRS Bits++ : %s', icmIdString);
        pString  = [ pString ' ' icmConfig ];
        Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', idString, shader, pString);        
    else
        % Setup for Color++ mode:

        % No support for overlays in Color++ mode:
        useOverlay = 0;

        if colorConversionMode == 2
            % Load "bilinear" Bits++ Color++ formatting shader for bilinear
            % sampling/averaging between adjacent even/odd pixel columns:
            shader = LoadGLSLProgramFromFiles('Bits++_Color++_BilinearFormattingShader', debuglevel, icmShaders);
        else
            % Load "classic" Bits++ Color++ formatting shader for non-interpolated
            % sampling:
            shader = LoadGLSLProgramFromFiles('Bits++_Color++_FormattingShader', debuglevel, icmShaders);
        end
                
        if colorConversionMode == 2
            % "Bilinear" mode: Aspect ratio correct, full-width source
            % framebuffer. Adjacent even/odd pixels are combined to a
            % single output pixel via averaging, ie., the output color is
            % the mean value of adjacent even/odd pixels:
            
            % Empty pString, no scaling needed:
            pString  = '';
        else
            if colorConversionMode == 0
                % "Classic" mode: Aspect ratio distorted half-width source framebuffer:
                sampleSpacing = 0.5;
                pString  = 'Scaling:2.0:1.0';
            end

            if colorConversionMode == 1
                % "Subsample" mode: Aspect ratio correct, full-width source
                % framebuffer, but sampled only at even pixel location, ie.
                % each second pixel column is skipped:
                sampleSpacing = 1.0;

                % Empty pString, no scaling needed:
                pString  = '';
            end

            glUseProgram(shader);
            glUniform1f(glGetUniformLocation(shader, 'sampleSpacing'), sampleSpacing);
            glUseProgram(0);
        end
        
        % Now enable output formatter hook chain and load them with the special Bits++
        % Color++ data formatting shader: We append the shader because it
        % absolutely must be the last shader to execute in that chain!
        % We apply a scaling of 2.0 in horizontal direction for the output
        % blit, to take the fact into account that the internal window
        % buffers only have half display width.
        idString = sprintf('Color++ output formatting shader for CRS Bits++ : %s', icmIdString);
        pString  = [ pString ' ' icmConfig ];
        Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', idString, shader, pString);        
    end

    % Setup shaders image source as the first texture unit, this is by
    % definition of how the imaging pipe works. Don't think really needed,
    % as this is the default, but its good practice to not rely on such
    % things...
    glUseProgram(shader);
    glUniform1i(glGetUniformLocation(shader, 'Image'), 0);
    glUseProgram(0);

    % Perform any setup steps that may be needed by the color correction
    % routines. Must be called after 'shader' creation and attachment to
    % the imaging pipe:
    PsychColorCorrection('ApplyPostGLSLLinkSetup', win, 'FinalFormatting');
    
    % Test accuracy/correctness of GPU's rasterizer for different output
    % positioning methods: Return (non-zero) dx,dy offsets, if any:
    [rpfx, rpfy, rpix, rpiy, vix] = PsychGPURasterizerOffsets(win, drivername); %#ok<ASGLU>
        
    if rpix~=0
        tlockXOffset = -rpix;
        fprintf('%s: Applying corrective horizontal DIO T-Lock offset of %i pixels for buggy graphics card driver. Will hopefully fix it...\n', drivername, tlockXOffset);        
    end

    if targetdevicetype == 1 && checkGPUEncoders
        % Perform DataPixx builtin diagnostics to detect problems with
        % wrong GPU gamma tables or GPU dithering:
        checkGPUEncoders = 0;
        if PsychDataPixx('CheckGPUSanity', win, tlockXOffset)
            % Ohoh, trouble ahead! The driver detected problems with the
            % GPU and wasn't able to auto-correct them.
            fprintf('%s: CAUTION! DataPixx internal diagnostic detected problems with your graphics card driver which it could not correct by itself!\n', drivername);
        end
    end
    
    if targetdevicetype == 0 && checkGPUEncoders
        % Perform Bits# builtin diagnostics to detect problems with
        % wrong GPU gamma tables or GPU dithering:
        checkGPUEncoders = 0;
        if doCheckGPUSanity(win, tlockXOffset, bitsSharpPort)
            % Ohoh, trouble ahead! The driver detected problems with the
            % GPU and wasn't able to auto-correct them.
            fprintf('%s: CAUTION! Bits# internal diagnostic detected problems with your graphics card driver which it could not correct by itself!\n', drivername);
        end
    end

    % Connection to Bits# established?
    if ~isempty(bitsSharpPort)
        % Yes. Switch it to target display mode:
        IOPort('Write', bitsSharpPort, targetMode);
        fprintf('BitsPlusPlus: Switching Bits# device to %s video mode. Will take about 5 seconds...\n', targetModeName);

        % Wait 5 seconds. In the worst case, if the diagnostic/status screen on Bits# was
        % active, the video mode switch can take that long to stabilize:
        WaitSecs('YieldSecs', 5);
    end

    % Enable framebuffer output formatter: From this point on, all visual
    % output will be reformatted to Bits++ framebuffer format at each
    % invokation of Screen('DrawingFinished') or Screen('Flip'), whatever
    % comes first.
    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');

    % When using the overlay, we need to allow for CLUT updates as well, so
    % usercode can define and change overlay colors:
    if useOverlay
        % Now enable finalizer hook chains and load them with the special Bits++
        % command for T-Lock based Bits++ internal CLUT updates:
        if vix~=0
            % vix is wrong offset, therefore negate it to get corrective offset:
            vix = -vix;
            fprintf('OpenWindow%sWithOverlay: Applying corrective horizontal offset of %i pixels for buggy graphics card driver. Will hopefully fix it...\n', mononame, vix);
            offsetstring = sprintf('xPosition=%i', vix);
        else
            offsetstring = '';
        end

        if targetdevicetype == 0
            Screen('HookFunction', win, 'PrependBuiltin', 'LeftFinalizerBlitChain', 'Builtin:RenderClutBits++', offsetstring);
        end

        if targetdevicetype == 1
            rclutcmd = 'PsychDataPixx(1, IMAGINGPIPE_GAMMATABLE);';
            Screen('HookFunction', win, 'PrependMFunction', 'LeftFinalizerBlitChain', 'Upload new clut into DataPixx callback', rclutcmd);
        end

        Screen('HookFunction', win, 'Enable', 'LeftFinalizerBlitChain');

        if ~isempty(stereomode) && ismember(stereomode, [1, 4, 5, 10, 11])
            % Enable CLUT updates on right stereo buffer as well:

            % CRS devices, frame-sequential stereo or native quad-buffered OpenGL stereo,
            % or dual-window stereo (mostly for OSX dual-display):
            if targetdevicetype == 0 && ismember(stereomode, [1, 10, 11])
                Screen('HookFunction', win, 'PrependBuiltin', 'RightFinalizerBlitChain', 'Builtin:RenderClutBits++', offsetstring);
            end

            % Dual CRS devices, dual-display stereo on Linux or Windows:
            if targetdevicetype == 0 && ismember(stereomode, [4, 5]) && dualdevices == 1
                % Here one PTB window drives both video outputs (split-window style) and
                % thereby both connected CRS devices. In order to set a CLUT on both of
                % these devices, we need to blit the Clut T-Lock into the right half of
                % our window as well. However, we use the LeftFinalizerBlitChain for this,
                % as the RightFinalizerBlitChain is not operational in stereomode 4 and 5:
                offsetstring = sprintf('xPosition=%i', vix + Screen('WindowSize', win, 1) / 2);
                Screen('HookFunction', win, 'PrependBuiltin', 'LeftFinalizerBlitChain', 'Builtin:RenderClutBits++', offsetstring);
            end

            % VPixx devices, frame-sequential stereo or native quad-buffered OpenGL stereo:
            if targetdevicetype == 1 && ismember(stereomode, [1, 11])
                Screen('HookFunction', win, 'PrependMFunction', 'RightFinalizerBlitChain', 'Upload new clut into DataPixx callback', rclutcmd);
            end

            Screen('HookFunction', win, 'Enable', 'RightFinalizerBlitChain');
        end
        
        % Load an identity CLUT into the Bits++ to start with:
        linear_lut =  repmat(linspace(0, 1, 256)', 1, 3);
        Screen('LoadNormalizedGammaTable', win, linear_lut, 2);
    end
    
    if targetdevicetype == 0
        % Setup finalizer callback for DIO T-Lock updates:
        tlockhandle = SetupDIOFinalizer(win, stereomode);
    end
    
    % Restore old graphics preferences:
    Screen('Preference', 'Enable3DGraphics', ogl);

    % Set color range to 0.0 - 1.0: This makes more sense than the normal
    % 0-255 values. Try to disable color clamping. This may fail and
    % produce a PTB warning, but if it succeeds then we're better off for
    % the 2D drawing commands...
    Screen('ColorRange', win, 1, 0);

    % Set Screen background clear color, in normalized 0.0 - 1.0 range:
    if (max(clearcolor) > 1) && (all(round(clearcolor) == clearcolor))
        % Looks like someone's feeding old style 0-255 integer values as
        % clearcolor. Output a warning to tell about the expected 0.0 - 1.0
        % range of values:
        warning(sprintf('\n\n%s: You specified a ''clearcolor'' argument for the OpenWindow command that looks \nlike an old 0-255 value instead of the wanted value in the 0.0-1.0 range. Please update your code for correct behaviour.', drivername)); %#ok<WNTAG,SPWRN>
    end
    
    % Set the background clear color via old fullscreen 'FillRect' trick,
    % followed by a flip:
    Screen('FillRect', win, clearcolor);
    Screen('Flip', win);
    
    % Check validation:
    if ~validated
        % Validate imaging / display pipeline:

        % From usercode drawing commands to framebuffer...
        ValidateBitsPlusImaging(win, 0, devname, 1);

        % ... and from framebuffer to display device:
        ValidateBitsPlusImaging(win, 0, devname, 2);
    end

    % Reset validation flag after first run:
    validated = 0;

    % Reset colorConversionMode after opening the window. It is a one-shot
    % parameter, but not before storing a cached copy in the per-window
    % vector:
    if ~isempty(colorConversionMode)
        colorConversionModeWin(win) = colorConversionMode;
    else
        colorConversionModeWin(win) = -1;
    end
    colorConversionMode = [];
    
    % Ready!
    return;
end

if strcmpi(cmd, 'GetOverlayWindow')
    % Assign onscreen window index:
    if nargin < 2 || isempty(arg) || ~isa(arg, 'double')
        error('%s: "GetOverlayWindow" called without valid onscreen window handle.', drivername);
    end
    win = arg;
    
    if win < 1 || win > length(OverlayWindows)
        error('%s: "GetOverlayWindow": No overlay associated with given onscreen window.', drivername);
    end

    if OverlayWindows(win) == 0
        error('%s: "GetOverlayWindow": No overlay associated with given onscreen window.', drivername);
    end

    % Ok, this 'win'dow has an overlay: Return its offscreen 'win'dow handle:
    win = OverlayWindows(win);
    % And the defining rectangle of the overlay:
    winRect = Screen('Rect', win);
    
    return;
end

if strcmpi(cmd, 'CheckGPUSanity')
    if length(varargin) < 1 || isempty(varargin{1})
        injectFault = 0;
    else
        injectFault = varargin{1};
    end
    
    win = doCheckGPUSanity(arg, dummy, bitsSharpPort, injectFault);
    return;
end

if strcmpi(cmd, 'UseFE1StereoGoggles')
    % Assign onscreen window index:
    if nargin < 2 || isempty(arg) || ~isa(arg, 'double') || Screen('WindowKind', arg) ~= 1
        error('%s: "UseFE1StereoGoggles" called without valid onscreen window handle.', drivername);
    end
    win = arg;

    winfo = Screen('GetWindowInfo', win);
    if ~ismember(winfo.StereoMode, [1,11])
        % No frame-sequential mode, no point in having sync lines -> No operation.
        fprintf('UseFE1StereoGoggles: Info: Provided onscreen window is not switched to frame-sequential stereo mode. Call ignored.\n');
        return;
    end

    if nargin < 3 || isempty(dummy) || ~ismember(dummy, [0, 1])
        error('%s: "UseFE1StereoGoggles" called without valid enable flag.', drivername);
    end
    enableIt = dummy;

    % Disable instead of enable requested?
    if ~enableIt
        % Existing hookslots found?
        lslot = Screen('HookFunction', win, 'Query', 'LeftFinalizerBlitChain', 'RenderStereoFE1LeftTLock');
        rslot = Screen('HookFunction', win, 'Query', 'RightFinalizerBlitChain', 'RenderStereoFE1RightTLock');
        if (lslot ~= -1) && (rslot ~= -1)
            % Yes. Delete them to disable the T-Lock blitting:
            Screen('HookFunction', win, 'Remove', 'LeftFinalizerBlitChain' , lslot);
            Screen('HookFunction', win, 'Remove', 'RightFinalizerBlitChain', rslot);
        end

        % Now the T-Lock blitting for frame-sequential stereo is off and
        % the shutter goggles are in some unknown state. Set them to
        % requested state, which defaults to 0 for "both eyes open":
        if nargin < 4 || isempty(varargin{1}) || varargin{1} ~= 1
            shutterState = 0;
        else
            shutterState = 1;
        end

        % Schedule T-Lock for setting final shutter state at next flip:
        [encodedDIOdata, Mask, Data, Command] = bitsGoggles(shutterState, shutterState);
        BitsPlusPlus('DIOCommand', win, 1, Mask, Data, Command);

        % Execute!
        Screen('Flip', win);

        % Ok, shutter glasses should be open, we are done.
        return;
    end

    % Enable case - Do the setup dance:

    % Build display list for left-eye image T-Lock:
    displist = glGenLists(1);
    glNewList(displist, GL.COMPILE);

    % Generate DIO T-Lock image as Matlab matrix to open the
    % right shutter at the end of scanout of the left stereo image:
    tlockdata = bitsGoggles(1, 0);

    % Convert from Matlab matrix to OpenGL pixel format:
    encodedDIOdata = uint8(zeros(3, 508));

    % Pack 3 separate RGB planes into rows 1,2,3. As Matlabs data format is
    % column major order, this will end up as tightly packed pixel array in
    % format RGBRGBRGB.... just as glDrawPixels likes it.
    encodedDIOdata(1,:) = tlockdata(1,:,1);
    encodedDIOdata(2,:) = tlockdata(1,:,2);
    encodedDIOdata(3,:) = tlockdata(1,:,3);

    % Place T-Lock at start position (0, display height - 1):
    [xDIO, yDIO] = Screen('WindowSize', win, 1);
    yDIO = yDIO - 1;

    % Add command sequence for this T-Lock code to display list:
    glRasterPos2i(tlockXOffset, yDIO);
    glDrawPixels(508, 1, GL.RGB, GL.UNSIGNED_BYTE, encodedDIOdata);

    % Finish left display list;
    glEndList;
    pStringLeft = sprintf('glCallList(%i);', displist);

    % Build display list for right-eye image T-Lock:
    displist = glGenLists(1);
    glNewList(displist, GL.COMPILE);

    % Generate DIO T-Lock image as Matlab matrix to open the
    % left shutter at the end of scanout of the right stereo image:
    tlockdata = bitsGoggles(0, 1);

    % Convert from Matlab matrix to OpenGL pixel format:
    encodedDIOdata = uint8(zeros(3, 508));

    % Pack 3 separate RGB planes into rows 1,2,3. As Matlabs data format is
    % column major order, this will end up as tightly packed pixel array in
    % format RGBRGBRGB.... just as glDrawPixels likes it.
    encodedDIOdata(1,:) = tlockdata(1,:,1);
    encodedDIOdata(2,:) = tlockdata(1,:,2);
    encodedDIOdata(3,:) = tlockdata(1,:,3);

    % Place T-Lock at start position (0, display height - 1):
    [xDIO, yDIO] = Screen('WindowSize', win, 1);
    yDIO = yDIO - 1;

    % Add command sequence for this T-Lock code to display list:
    glRasterPos2i(tlockXOffset, yDIO);
    glDrawPixels(508, 1, GL.RGB, GL.UNSIGNED_BYTE, encodedDIOdata);

    % Finish left display list;
    glEndList;
    pStringRight = sprintf('glCallList(%i);', displist);

    % Ok 'win'dowhandle is fine. Any blue line sync function already applied?
    lslot = Screen('HookFunction', win, 'Query', 'LeftFinalizerBlitChain', 'Builtin:RenderStereoSyncLine');
    rslot = Screen('HookFunction', win, 'Query', 'RightFinalizerBlitChain', 'Builtin:RenderStereoSyncLine');

    % Existing hookslots found?
    if (lslot ~= -1) && (rslot ~= -1)
        % Yes. Need to recreate with T-Lock for FE1 instead of standard blue-line sync lines:

        % Destroy old ones:
        Screen('HookFunction', win, 'Remove', 'LeftFinalizerBlitChain' , lslot);
        Screen('HookFunction', win, 'Remove', 'RightFinalizerBlitChain', rslot);

        % Insert new slots at former position of the old ones:
        posstring = sprintf('InsertAt%iMFunction', lslot);
        Screen('Hookfunction', win, posstring, 'LeftFinalizerBlitChain', 'RenderStereoFE1LeftTLock', pStringLeft);
        posstring = sprintf('InsertAt%iMFunction', rslot);
        Screen('Hookfunction', win, posstring, 'RightFinalizerBlitChain', 'RenderStereoFE1RightTLock', pStringRight);
    else
        % No. Create new slots at end:
        Screen('Hookfunction', win, 'AppendMFunction', 'LeftFinalizerBlitChain', 'RenderStereoFE1LeftTLock', pStringLeft);
        Screen('Hookfunction', win, 'AppendMFunction', 'RightFinalizerBlitChain', 'RenderStereoFE1RightTLock', pStringRight);
    end

    fprintf('UseFE1StereoGoggles: Info: Enabling automatic generation of sync signals for external FE1 shutter glasses.\n');

    % Enable finalizer hookchains, if not already enabled:
    Screen('HookFunction', win, 'Enable', 'LeftFinalizerBlitChain');
    Screen('HookFunction', win, 'Enable', 'RightFinalizerBlitChain');

    return;
end

error('%s: Unknown subcommand provided. Read "help BitsPlusPlus".', drivername);
end

% Helper function: Check if system already validated for current settings:
function ValidateBitsPlusImaging(win, writefile, devname, validationType)

    % Compute fingerprint of this system configuration:
    validated = 0;
    global GL;

    screenid = Screen('WindowScreenNumber', win);
    [w, h] = Screen('WindowSize', screenid);
    d = Screen('PixelSize', win);
    v = Screen('Version');
    v = v.version;
    c = Screen('Computer');
    if IsLinux || IsOSX
        c = c.kern.version;
    else
        c = c.system;
    end

    gfxconfig = [ glGetString(GL.VENDOR) ':' glGetString(GL.RENDERER) ':' glGetString(GL.VERSION) ];
    gfxconfig = sprintf('VTYPE %i : %s : Screen %i : Resolution %i x %i x %i : ScreenVersion = %s : System = %s', validationType, gfxconfig, screenid, w, h, d, v, c);

    if ~writefile
        % Check if a validation file exists and if it contains this
        % configuration:
        fid = fopen([PsychtoolboxConfigDir 'ptbbitsplusplusvalidationfile.txt'], 'r');
        if fid~=-1
            while ~feof(fid)
                vconf = fgetl(fid);
                if strcmp(vconf, gfxconfig)
                    validated = 1;
                    break;
                end
            end
            fclose(fid);
        end

        if ~validated
            fprintf('\n\n------------------------------------------------------------------------------------------------------------------\n')
            fprintf('\n\nThis specific configuration of graphics hardware, graphics driver and Psychtoolbox version has not yet been tested\n');
            fprintf('for correct working with %s for the given display screen, screen resolution and color depths.\n\n', devname);
            if validationType == 1
                fprintf('Please run the test script "BitsPlusImagingPipelineTest(%i);" once, so this configuration can be verified.\n', Screen('WindowScreenNumber', win));
            else
                if ~isempty(strfind(devname, 'Pixx'))
                    fprintf('Please run the test script "BitsPlusIdentityClutTest(%i, 1);" once, so this configuration can be verified.\n', Screen('WindowScreenNumber', win));
                    fprintf('Alternatively "BitsPlusIdentityClutTest(%i, 1, [], 1);" can be used for verification in L48 mode.\n', Screen('WindowScreenNumber', win));
                else
                    fprintf('Please run the test script "BitsPlusIdentityClutTest(%i);" once, so this configuration can be verified.\n', Screen('WindowScreenNumber', win));
                    fprintf('Alternatively "BitsPlusIdentityClutTest(%i, 0, [], 1);" can be used for verification in Bits++ mode.\n', Screen('WindowScreenNumber', win));
                end
            end
            fprintf('After that test script suceeded, re-run your experiment script.\nThanks.\n');
            fprintf('\n');
            fprintf('Configuration to verify: %s\n', gfxconfig);

            RestoreCluts;
            sca; ShowCursor; Priority(0);

            % Perform device shutdown for display device:
            if ~isempty(strfind(devname, 'Pixx'))
                PsychDataPixx('ResetOnWindowClose');
            else
                BitsPlusPlus('ResetOnWindowClose');
            end

            error('Configuration not yet verified. Please do it now.');
        end
    end

    if writefile
        % Append current configuration to file to mark it as verified:
        [fid msg]= fopen([PsychtoolboxConfigDir 'ptbbitsplusplusvalidationfile.txt'], 'a');
        if fid == -1
            RestoreCluts;
            sca;
            error('Could not write validation file %s to filesystem [%s].', [PsychtoolboxConfigDir 'ptbbitsplusplusvalidationfile.txt'], msg);
        end

        % Append line:
        fprintf(fid, [gfxconfig '\n']);
        fclose(fid);
    end
end

% Helper function for setup of finalizer blit chains in all modes. Sets up
% callback into our file for T-Lock drawing etc...
function displist = SetupDIOFinalizer(win, stereomode)

    % Generate unique display list handle for later use:
    displist = glGenLists(1);

    % Now enable finalizer hook chains and load them with the special Bits++
    % command for T-Lock based Bits++ DIO updates:
    Screen('HookFunction', win, 'PrependMFunction', 'LeftFinalizerBlitChain', 'Render T-Lock DIO data callback', 'BitsPlusPlus(1);');
    Screen('HookFunction', win, 'Enable', 'LeftFinalizerBlitChain');

    if ~isempty(stereomode) && ismember(stereomode, [1, 11])
        % This is only needed on quad-buffered stereo contexts.
        Screen('HookFunction', win, 'PrependMFunction', 'RightFinalizerBlitChain', 'Render T-Lock DIO data callback',  'BitsPlusPlus(1);');
        Screen('HookFunction', win, 'Enable', 'RightFinalizerBlitChain');
    end

end

function scanline = BitsSharpGetScanline(bitsSharpPort, lineNr, nrPixels)

    % Emit request for video scanline to Bits#
    IOPort('Write', bitsSharpPort, [sprintf('$GetVideoLine=[%i,%i]', lineNr, nrPixels) char(13)]);
    IOPort('Flush', bitsSharpPort);
    
    % First we do a blocking read for the minimum number of bytes expected, which is the
    % length of the header in chars + at least 3 chars per pixel (format 00; 01; 02; ... 99;):
    rawline = IOPort('Read', bitsSharpPort, 1, length('#GetVideoLine;') + nrPixels * 3);

    % Then we iterate over non-blocking reads until we find the end-of-date terminator code 13:
    while ~isempty(rawline) && (rawline(end) ~= 13)
        WaitSecs('YieldSecs', 0.001);
        rawline = [rawline IOPort('Read', bitsSharpPort)]; %#ok<AGROW>
    end

    % Cut away header:
    rawline = rawline(length('#GetVideoLine;')+1:end);
    if isempty(rawline)
        warning('BitsSharpGetScanline: Empty pixelline returned!'); %#ok<WNTAG>
        scanline = [];
        return;
    end

    % Convert into integer array:
    rawline = char(rawline);
    rawline(rawline == ';') = ' ';
    scanline = sscanf(rawline, '%d');
    if (length(scanline) ~= 3 * nrPixels)
        warning('BitsSharpGetScanline: Incomplete pixelline %s with only %i elements (less than %i) returned!', rawline, length(scanline), 3 * nrPixels); %#ok<WNTAG>
        scanline = [];
    else
        scanline = uint8(reshape(scanline, 3, nrPixels));
    end

    return;
end

function rc = doCheckGPUSanity(win, xoffset, bitsSharpPort, injectFault)
    % If this isn't a connected Bits#, simply no-op with success return code:
    if isempty(bitsSharpPort)
        rc = 0;
        return;
    end

    if nargin < 4 || isempty(injectFault)
        injectFault = 0;
    end
    
    % Execute test and optimization (tweaking) procedure which uses
    % onscreen window 'win' for sending test stimuli to the Bits# device,
    % and use the builtin measurement functions of that device to drive
    % the tweaking procedure. Return success status, 0 = Success, 1 = Failure.
    rc = PsychGPUTestAndTweakGammaTables(win, xoffset, 1, injectFault, bitsSharpPort);

    return;
end
