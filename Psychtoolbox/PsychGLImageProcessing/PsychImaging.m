function [rc, winRect] = PsychImaging(cmd, varargin)
% rc = PsychImaging(subcommand [,arg1][,arg2][....,argn]) - Control common
% functions of the Psychtoolbox GPU image processing pipeline.
%
% This function allows you to setup and control various aspects and common
% functions of the Psychtoolbox image processing pipeline in a simple way.
% Various standard scenarious can be conveniently set up with this routine,
% e.g., geometric transformations of your stimulus image, various types of
% display correction, ...
%
% If you want to perform less common, unusual or simply not yet supported tasks
% with the pipeline, use the low-level Screen('HookFunction', ...)
% interface instead and have a peek in the M-File code for the
% PsychImaging.m file to learn about the low-level interface.
% See "help PsychGLImageprocessing" for more info.
%
%
% Subcommands and their meaning:
%
% PsychImaging('PrepareConfiguration');
% - Prepare setup of imaging pipeline for onscreen window.
% This is the first step in the sequence of configuration steps.
%
%
% PsychImaging('AddTask', whichChannel, whichTask [,param1]...);
% - Add a specific task or processing requirement to the list of actions
% to be performed by the pipeline for the currently selected onscreen
% window. 'whichChannel' is a string with the name of the channel to
% configure:
%
% 'LeftView' applies the action to the processing channel
% for the left-eye view of a stereo configuration. 'RightView' applies the
% action to the right-eye view of a stereo configuration. 'AllViews' applies
% the action to both, left- and right eye view channels of a stereo
% configuration or to the single monoscopic channel of a mono display
% configuration. Other options are 'Compositor', 'FinalFormatting' and
% 'Finalizer' for special purpose channels. Set this to 'General' if the
% command doesn't apply to a specific view, but is a general requirement.
%
% 'whichTask' contains the name string of one of the supported
% actions:
%
% * 'UseGPGPUCompute' Enable use of GeneralPurposeGPU computing support.
%   This prepares use of Psychtoolbox functions which are meant to
%   interface with, or take advantage of, the general purpose computation
%   capabilities of modern graphics processing units and other massively
%   parallel compute acceleration hardware, e.g., DSP's, or multi-core
%   processors. Interfacing with such hardware is done via common standard
%   compute API's like NVidia's CUDA or the cross-platform OpenCL API.
%
%   Use of this function often requires specific modern GPU hardware and
%   the installation of additional driver software, e.g., NVidia's freely
%   available CUDA SDK and runtime, or the free and open-source GPUmat
%   toolbox. Read 'help PsychGPGPU' for further info.
%
%   This function just detects and selects supported GPU compute API's for
%   use with Psychtoolbox and initializes them and some Psychtoolbox
%   function to take advantage if appropriate. While you could use those
%   API's by themselves without calling this init function, Psychtoolbox
%   builtin processing functions would not be able to take advantage of the
%   API's or perform efficient and fast data exchange with them.
%
%   Usage: PsychImaging('AddTask', 'General', 'UseGPGPUCompute', apitype [, flags]);
%
%   'apitype' Allows selection of the compute API to use. The value 'Auto'
%   leaves the choice to Psychtoolbox. The value 'GPUmat' selects the
%   high-level, free and open-source GPUmat compute toolkit for Matlab.
%   Currently no other choices are supported, but this is expected to
%   change in the future.
%
%   'flags' An optional string of keyword flags to determine behaviour.
%   There aren't any flags defined yet.
%
%
% * 'SideBySideCompressedStereo' Ask for stereo display in a horizontally
%   compressed side-by-side format. Left and Right eye images are drawn at
%   full framebuffer resolution by usercode. Screen('Flip', ...) draws them
%   horizontally compressed side-by-side to each other. They are scanned
%   out to the display device this way and then the display device itself
%   uncompresses them back to full resolution and displays them
%   stereoscopically, typically via built-in alternating frame-sequential
%   stereo with stereo goggles, but other methods are conceivable. This is
%   one popular stereo frame packing format for stereo on HDMI display
%   devices. Once you've set up a stereo display mode via PsychImaging, you
%   can tweak its specific parameters by calling the function
%   SetCompressedStereoSideBySideParameters().
%
%   Usage: PsychImaging('AddTask', 'General', 'SideBySideCompressedStereo');
%
%
% * 'InterleavedColumnStereo' Ask for stereo display in interleaved mode.
%   The output image is composed from the lefteye and righteye stereo
%   buffers by interleaving their content: Even columns are filled with
%   content from the left buffer, odd columns are filled with content from
%   the right buffer, i.e., Col 0 = Left col 0, Col 1 = Right Col 0, Col 2
%   = Left col 1, Col 3 = Right col 1, ....
%
%   This mode is useful for driving some auto-stereoscopic displays. These
%   use either some vertical parallax barriers or vertical lenticular
%   lense sheets. These direct light from even columns to one eye, light
%   from odd columns to the other eye.
%
%   Usage: PsychImaging('AddTask', 'General', 'InterleavedColumnStereo', startright);
%
%   If 'startright' is zero, then even columns are taken from left buffer. If
%   'startright' is one, then even columns are taken from the right buffer.
%
%   You can use the RemapMouse() function to correct GetMouse() positions
%   for potential geometric distortions introduced by this function.
%
%
% * 'InterleavedLineStereo' Ask for stereo display in interleaved mode.
%   The output image is composed from the lefteye and righteye stereo
%   buffers by interleaving their content: Even lines are filled with
%   content from the left buffer, odd lines are filled with content from
%   the right buffer, i.e., Row 0 = Left row 0, Row 1 = Right row 0, Row 2
%   = Left row 1, Row 3 = Right row 1, ....
%
%   This mode is useful for driving some types of stereo devices and
%   goggles, e.g., the iGlasses 3D Video goggles in interleaved stereo
%   mode.
%
%   Usage: PsychImaging('AddTask', 'General', 'InterleavedLineStereo', startright);
%
%   If 'startright' is zero, then even lines are taken from left buffer. If
%   'startright' is one, then even lines are taken from the right buffer.
%
%   You can use the RemapMouse() function to correct GetMouse() positions
%   for potential geometric distortions introduced by this function.%
%
%
% * 'DualWindowStereo' Ask for stereo display in dual-window mode (stereomode 10)
%
%   Only use this function under MacOSX. If possible on your setup and OS,
%   rather use a single window, spanning both stereo display outputs, and use
%   stereomode 4 or 5 to display dual-display stereo. That is much more
%   efficient in terms of speed, computational load and memory consumption,
%   also potentially more robust with respect to visual stimulation timing.
%
%   Usage: PsychImaging('AddTask', 'General', 'DualWindowStereo', rightEyeScreen [, rightEyeWindowRect]);
%
%   The left-eye image will be displayed on the screen and at a location
%   specified as usual via PsychImaging('Openwindow', screenid, ..., rect);
%   The right eye image will be displayed on screen 'rightEyeScreen'. If
%   the optional 'rightEyeWindowRect' is specified, then the right eye image
%   is not displayed in a fullscreen window, but in a window with the bounding
%   rectangle 'rightEyeWindowRect'.
%
%
% * 'UseVirtualFramebuffer' Ask for support of virtual framebuffer, even if
%   it isn't strictly needed, given the remaining set of requirements. Most
%   of the tasks require such a framebuffer - it gets enabled anyway. In a
%   few cases, e.g., to simplify your code (no need for special cases), it
%   may be useful to activate such a framebuffer even if it isn't strictly
%   needed. This option activates a minimal buffer with 8 bits per color
%   cmponent fixed point precision.
%
%   Usage: PsychImaging('AddTask', 'General', 'UseVirtualFramebuffer');
%
%
% * 'UseFineGrainedTiming' Ask for use of fine-grained stimulus onset timing.
%   Ask Screen('Flip', window, when) to try to show the new visual stimulus
%   close to (= ideally exactly at) target time 'when', instead of showing the
%   stimulus at the next video refresh frame boundary with time t >= 'when', ie.
%   try to schedule stimuli with better timing granularity than what is given by
%   the multiples of a video refresh cycle duration of the connected display.
%
%   This uses a technique known as "Variable Refresh Rate" or shorthand "VRR" if
%   your operating system and display driver and graphics card and cable and
%   display device supports it at the current system configuration. Otherwise,
%   PsychImaging('OpenWindow', ...) will fail if this task can't be achieved.
%
%   For a list of requirements with respect to graphics cards, display devices,
%   operating systems, drivers and general system configuration to make this work,
%   read "help VRRSupport". This will also tell you about limitations and caveats
%   wrt. this task.
%
%   There may be different methods of implementing such fine-grained timing.
%   The optional 'method' parameter allows you to select a specific method.
%   Using the keyword 'Auto' or omitting the 'method' parameter ([]) will leave
%   the choice of optimal method to Screen(). Currently the methods 'Simple' and
%   'OwnScheduled' are implemented: 'Simple' follows the most naive approach of VRR,
%   simply requesting immediate flip, or in case a 'when' time is specified, waiting
%   until 'when', then requesting a flip. 'OwnScheduled' selects a more sophisticated
%   scheduling method that tries to take other constraints like minimum and maximum
%   refresh rate of the display, display system state, and a specified 'styleHint'
%   for the given visual stimulation paradigm into account, and also tries to
%   compensate for some of the system timing jitter. This may provide higher
%   precision or stability in some stimulation scenarios, but may also cause
%   higher overhead or added latency in other scenarios.
%
%   The optional parameter 'styleHint' allows to give some general high level clue
%   into the nature of the visual stimulation paradigm, as far as presentation
%   timing is concerned. This allows to tune the scheduling method for VRR for
%   maximum precision and reliability for a stimulation paradigm of the specified
%   nature. Omitting the parameter, or setting it to 'None', will ask the algorithm
%   to make almost no assumptions about the nature of the stimulation, but to choose
%   some "one size fits all ok'ish" setup, or try to auto-detect how to tune itself for
%   the running paradigm. The currently supported 'styleHint' settings are:
%
%   - 'None' = Make no assumptions, try some "One size fits all hopefully ok'ish" setup.
%
%   The optional parameter 'vrrMinRefreshHz' allows to specify the lowest video
%   refresh rate that your display can reliably run at. If the parameter is
%   omitted, Screen() will try to auto-detect this display property, or failing
%   that, it will use a reasonable default.
%
%   Usage: PsychImaging('AddTask', 'General', 'UseFineGrainedTiming' [, method='Auto'][, styleHint='None'][, vrrMinRefreshHz]);
%
%
% * 'UseSubpixelDrive' Ask to take advantage of the so-called "Subpixel Drive"
%   mode of certain monochromatic medical imaging displays like, e.g., the
%   "Eizo RadiForce GS-521". This monitor essentially has a RGB panel with
%   horizontal RGB subpixels, but with the color filters removed, so each
%   pixel is horizontally split up into 3 luminance subpixels and these can
%   be individually addressed by packing 3 horizontally adjacent 8 bit stimulus
%   luminance pixels into the "RGB" color channels of an output pixel. Use
%   of this task will create a virtual framebuffer three times the width of
%   the output framebuffer and then pack triplets of three horizontal luminance
%   pixels into one output pixel, to triple the effective horizontal resolution.
%
%   Usage: PsychImaging('AddTask', 'General', 'UseSubpixelDrive');
%
%
% * 'UseRetinaResolution' Ask to prefer a framebuffer with the full native
%   resolution of attached HiDPI "Retina" displays on OSX, instead of a scaled
%   down lower resolution framebuffer with typically half the horizontal
%   and vertical resolution of the Retina display. This setting will be
%   ignored if the onscreen window is not displayed on a HiDPI "Retina"
%   display in a scaled display mode, or if the panel fitter is in use by
%   specifying the 'UsePanelFitter' task. By default, Screen() would use a
%   downscaled framebuffer on a Retina display under OSX and scale that low
%   resolution buffer up to full display panel resolution, just as Apples
%   OSX operating system does it by default. This in order to reduce
%   computational load, improve graphics performance, and avoid problems
%   with backward compatibility of old code. If you want to make full use
%   of the resolution of your HiDPI display, specify this task to tell
%   Screen() to use the full display panel resolution on OSX, even if this may
%   introduce some compatibility issues into your code and causes a decrease
%   in graphics performance due to the higher graphics rendering load.
%
%   If 'UseRetinaResolution' is used with a non-fullscreen window, ie.
%   the 'rect' parameter in PsychImaging('OpenWindow', ...) is provided
%   to specify the screen position and size of the window, note that
%   the size of the window rect returned by Screen('GlobalRect') and
%   Screen('Rect'), as well as of the returned rect of PsychImaging('OpenWindow')
%   will differ from the size of the 'rect' passed to 'OpenWindow'. 'rect's
%   passed into OpenWindow for positioning and sizing the window, as well
%   as the global position rect returned by Screen('GlobalRect') for the
%   current size and position of the onscreen window are expressed in global
%   desktop coordinates, in somewhat arbitrary units of virtual "points".
%   How such a point translates into display pixels depends on the operating
%   system, possibly the desktop GUI in use (on other systems than OSX), the
%   set of connected displays and their Retina or non-Retina resolutions.
%   The aim is that the coordinate system is somewhat consistent and meaningful
%   across all connected displays, for varying definitions of "consistent" and
%   "meaningful" on different operating systems, but the mapping of points to
%   physical screen pixels can be different on each connected display, at the
%   discretion of the operating system. You may get especially "interesting"
%   results if you try to move an onscreen window between screens, or let it
%   span multiple displays of different type and resolution.
%   The rect returned by PsychImaging('Openwindow') and Screen('Rect'), as
%   well as sizes returned by Screen('WindowSize') define the net useable
%   size of the window in display pixels. It is affected by all kind of
%   PsychImaging operations, e.g., selection of stereo modes, high bit depth
%   modes etc., but also by scaling on Retina displays in high res mode.
%   If 'UseRetinaResolution' is used on a Retina/HiDPI display, one typical
%   result will be that the size of the window in pixels reported by these
%   functions will be higher than the size in points, as one virtual point will
%   get represented by more than 1 pixel on a Retina display. Observing twice
%   the window size in pixels than in points is quite typical, but other
%   scaling factors are possible. The take home message for you is to specify
%   location and size of your stimuli based on the sizes and rects returned
%   by Screen('Rect'), PsychImaging('Openwindow') and Screen('Windowsize'), as
%   these are in units of display pixels, and *not* based on the virtual points
%   returned by Screen('GlobalRect'). The 2nd take home message is that you
%   should mostly use fullscreen windows for visual stimulation to avoid such
%   and other pitfalls.
%
%   Usage: PsychImaging('AddTask', 'General', 'UseRetinaResolution');
%
%
% * 'UseDisplayRotation' Ask to use builtin panel fitter exclusively for
%   rotating the framebuffer. This is useful if you want to turn your
%   display device from landscape (= normal upright) orientation into
%   portrait orientation (= rotated by 90 degrees clockwise or
%   counterclockwise). In such a case you will want to rotate the
%   framebuffer by 90 degrees as well, but you should *not* use the "rotate
%   monitor" function of your operating system for this purpose, as this
%   will very likely interfere with visual stimulus presentation timing and
%   timestamping! Use this task instead. It will perform rotation in a
%   similar way, but without severe interference to timing. However, there
%   is one limitation to this method: Multisample anti-aliasing currently
%   does not work if you use our framebuffer rotation.
%
%   Usage: PsychImaging('AddTask', 'General', 'UseDisplayRotation', angle);
%
%   'angle' is the desired rotation angle. The only values which will give
%   well defined and useful results are multiples of 90 degrees, useful
%   values are essentially 0, +90, -90 and 180 degrees for no rotation,
%   clockwise rotation, counterclockwise rotation and upside down rotation.
%
%   This function is mutually exclusive with 'UsePanelFitter', but if you
%   need to use both, you can omit 'UseDisplayRotation' and pass the
%   'angle' parameter to 'UsePanelFitter' instead, which also accepts an
%   'angle' parameter with the same meaning.
%
%   This function is not very mature yet: If you want to use the
%   panelfitter for anything beyond simple framebuffer rotation by 90
%   degree increments, you will likely hit bugs or limitations which will
%   require significant tinkering by you.
%
%
% * 'UsePanelFitter' Ask to use builtin panel fitter. This allows you to
%   define a virtual size for your onscreen window. The window will behave
%   as if it had that virtual size wrt. all size queries and drawing
%   operations. However, at Screen('Flip') time, the visual content of the
%   window will be resized by a fast scaling operation to the real size of
%   the windows framebuffer, ie., its real onscreen size. Scaling uses
%   bilinear interpolation or better for high quality results. After
%   rescaling to the real size, post-processing and display of your
%   stimulus image will proceed at full resolution. This function is useful
%   if you want to display a stimulus designed for a specific display
%   resolution on a display device of different higher or lower resolution.
%   Given that size and shape of the virtual framebuffer and real display
%   window may not match, the function provides you with multiple possible
%   choices on how to rescale your stimulus image, e.g., to maximize
%   display area, or to preserve the aspect ratio of the original image,
%   trading off displayed area etc.
%
%   Usage: PsychImaging('AddTask', 'General', 'UsePanelFitter', size, strategy [, srcRect, dstRect][, angle]);
%
%   'size' is a [width, height] vector defining the width x height of the
%   virtual window in pixels.
%
%   'strategy' a text string selecting the scaling method. Following settings are possible:
%
%   'Full' - Scale to full window size. Aspect ratio is not preserved,
%            unless the virtual window and the real onscreen windows 'rect'
%            already have the same aspect ratio, in which case this will be
%            a simple scaling operation.
%
%   'Aspect' - Scale to maximum size possible while preserving aspect
%              ratio. This will center the stimulus and add black
%              horizontal or vertical borders as neccessary.
%
%   'AspectWidth' - Scale aspect ratio preserving to cover full display
%                   width. Cut off top and bottom content if neccessary.
%
%   'AspectHeight' - Scale aspect ratio preserving to cover full display
%                    height. Cut off left and right content if neccessary.
%
%   'Centered' - Center stimulus without any scaling, add black borders
%                around stimulus or cut away border regions to get a
%                one-to-one mapping.
%
%   'Custom' - This works like the 'srcRect' and 'dstRect' parameters of
%              Screen('DrawTexture'): Cut out a 'srcRect' region from the
%              virtual framebuffer and display it in the 'dstRect' region.
%              'srcRect' and 'dstRect' are given in typical [left, top, right, bottom]
%              format.
%
%   'angle' is an optional rotation angle. If provided and non-zero, the
%   panelfitter will also rotate the output framebuffer by the given
%   rotation angle. Note: This doesn't work very well yet with most
%   framebuffer sizes and scaling strategies. What does work is if the
%   specified 'size' is identical to the onscreen windows size, or is its
%   transposed size (ie., if window is width x height pixels, then height x
%   width pixels will work as 'size' parameter) and the rotation angle is a
%   multiple of 90 degrees. This is mostly useful for display rotation from
%   landscape orientation into portrait orientation. Your mileage with
%   other configurations or rotation angles will vary.
%
%   Example: Suppose your real window covers a 1920 x 1080 display.
%
%   PsychImaging('AddTask', 'General', 'UsePanelFitter', [800 600], 'Aspect');
%   -> This would give you a virtual window of 800 x 600 pixels to draw
%   into and would rescale the 800 x 600 stimulus image to 1440 x 1080
%   pixels and display it centered on the 1920 x 1080 pixels display.
%   Aspect ratio would be correct and the image would cover the full height
%   1080 pixels of the display, but only 1440 out of 1920 pixels of its
%   width, thereby leaving black borders on the left and right side of your
%   stimulus.
%
%   PsychImaging('AddTask', 'General', 'UsePanelFitter', [800 600], 'AspectHeight');
%   -> Would do the same as above.
%
%   PsychImaging('AddTask', 'General', 'UsePanelFitter', [800 600], 'AspectWidth');
%   -> Would create a final image of 1920 pixels width, as you asked to
%   cover the full display width, aspect ratio would be correct, but the
%   top and bottom 75 pixels of your original stimulus would get cut away,
%   because they wouldn't fit after scaling without distorting the image.
%
%
% * 'UseFastOffscreenWindows' Ask for support of fast Offscreen windows.
%   These use a more efficient storage, backed by OpenGL framebuffer
%   objects (FBO's). Drawing into them isn't faster, but *switching*
%   between drawing into onscreen- and offscreen windows, or switching
%   between drawing into different offscreen windows is faster. They also
%   support a couple of other advanced features and performance
%   improvements in conjunction with the imaging pipeline.
%   If you only specify this task, then you'll get the benefit of fast
%   windows, without the cost of other features of the pipeline you might
%   not need.
%
%   Usage: PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
%
%
% * 'EnableCLUTMapping' Enable support for old-fashioned clut animation /
%   clut mapping. The drawn framebuffer image is transformed by applying a
%   color lookup table (clut). This is not done via the hardware gamma
%   tables as in the good ol' days, but by application of the clut via
%   image processing. Hardware gamma tables don't provide well defined
%   timing on modern hardware, therefore they aren't suitable anymore.
%
%   You can update the clut to be applied at the next Screen('Flip');
%   via the command Screen('LoadNormalizedGammatable', windowPtr, clut, 2);
%
%   'clut' needs to be a clutSize-by-3 matrix, with 'clutSize' slots and
%   one column for each of the red, green and blue color channels.
%
%   Setup command:
%
%   By default, a clut of 256 slots with (R,G,B) values is used, but you
%   can provide the optional 'clutSize' parameter to use clut's with more
%   slots. The maximum number depends on your GPU, but 2048 are typically
%   supported even on very low-end cards.
%
%   If you set 'highprecision' to 1, the clut will resolve values at more
%   than 8 bit per color channel on modern hardware. This usually only
%   makes sense if you also use a more than 8 bpc framebuffer with more
%   than 256 slots as clutSize.
%
%   Usage: PsychImaging('AddTask', whichView, 'EnableCLUTMapping' [, clutSize=256][, highprecision=0]);
%   Example: PsychImaging('AddTask', 'AllViews', 'EnableCLUTMapping');
%
%
% * 'FloatingPoint16Bit' Ask for a 16 bit floating point precision
%   framebuffer. This allows more than 8 bit precision for complex drawing,
%   compositing and image processing operations. It also allows
%   alpha-blending with signed color values and intermediate results that
%   are outside the displayable range, e.g., negative. Precision is about 3
%   digits behind the decimal point or 1024 discriminable displayable
%   levels. If you need higher precision, choose 'FloatingPoint32Bit'.
%
%   Usage: PsychImaging('AddTask', 'General', 'FloatingPoint16Bit');
%
%
% * 'FixedPoint16Bit' Ask for a 16 bit integer precision framebuffer.
%   On graphics hardware that supports this, a 16 bit signed integer
%   framebuffer will be created. Such a framebuffer can store intermediate
%   color values in the normalized range [-1.0 ; +1.0] with a precision of
%   15 bits per component. Only positive values between 0.0 and 1.0 are
%   displayable in the end though. If the graphics hardware does not support this,
%   a 16 bit unsigned integer framebuffer is tried instead. Such a framebuffer
%   allows for 16 bits of precision per color component. However, many graphics
%   cards do not support alpha-blending on such a framebuffer, and
%   intermediate out-of-range values (smaller than zero or bigger than one) aren't
%   supported either. Such values will be clamped to the representable [0.0 ; 1.0]
%   range instead. Additionally this mode is only supported on some graphics
%   hardware. It is a special purpose intermediate solution - more accurate
%   than 16 bit floating point, but less capable and less accurate than 32
%   bit floating point. If you need higher precision, choose 'FloatingPoint32Bit'.
%
%   The main sad reason this switch exists is because some graphics hardware or
%   graphics drivers do not support floating point precision textures and
%   framebuffers due to some ridiculous patent restrictions, but they do
%   support a 16 bit signed or unsigned integer precision format. The switch
%   is basically a workaround for the broken patent systems of many countries.
%
%   Usage: PsychImaging('AddTask', 'General', 'FixedPoint16Bit');
%
%
% * 'FloatingPoint32Bit' Ask for a 32 bit floating point precision
%   framebuffer. This allows more than 8 bit precision for complex drawing,
%   compositing and image processing operations. It also allows
%   alpha-blending with signed color values and intermediate results that
%   are outside the displayable range, e.g., negative. Precision is about
%   6.5 digits behind the dezimal point or 8 million discriminable displayable
%   levels. Be aware that only the most recent hardware (NVidia Geforce
%   8000 series, ATI Radeon HD 2000 series) is able to perform
%   alpha-blending at full speed in this mode. Enabling alpha-blending on
%   older hardware may cause a significant decrease in drawing performance,
%   or alpha blending may not work at all at this precision! If you'd like
%   to have both, the highest precision and support for alpha-blending,
%   specify 'FloatingPoint32BitIfPossible' instead. PTB will then try to
%   use 32 bit precision if this is possible in combination with alpha
%   blending. Otherwise, it will choose 16 bit precision for drawing &
%   blending, but 32 bit precision at least for the post-processing.
%
%   Usage: PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');
%
%
% * 'FloatingPoint32BitIfPossible' Ask PTB to choose the highest precision
%   that is possible on your hardware without sacrificing functionality like,
%   e.g., alpha-blending. PTB will choose the best compromise possible for
%   your hardware setup.
%
%   Usage: PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
%
%
% * 'NormalizedHighresColorRange' Ask PTB to use a normalized range of
%   color and luminance intensity levels in the interval [0; 1], ie. values
%   between zero and one for minimum and maximum intensity. Also ask for
%   unclamped colors -- intermediate results are allowed to take on
%   arbitrary values, e.g., also negative values. All Screen() 2D drawing
%   commands should operate at maximum color/luminance precision.
%
%   Usage: PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange' [, applyAlsoToMakeTexture]);
%
%   The command PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
%   is automatically executed if you used PsychDefaultSetup(featureLevel)
%   with a featureLevel of >= 2 at the top of your experiment script,
%   *except* that clamping is *not* disabled by default in this case! To
%   disable clamping you'd still need to add this task explicitely, as
%   unclamping may have unintended side effects on old graphics hardware.
%
%   The optional flag 'applyAlsoToMakeTexture' defaults to zero. If set to 1,
%   then a unit color range of expected input values in the [0; 1] range is
%   also applied to standard 8-Bit precision textures in Screen('MakeTexture')
%   if the provided Matlab imageMatrix is of double precision type instead of
%   uint8 type. This allows to specify standard textures in the same consistent
%   value range 0-1 as other drawing colors, for cleaner code. Such textures
%   will still be limited to 0-1 range and only resolved into 256 intensity
%   levels, unless you also set the optional 'floatprecision' flag in Screen('MakeTexture')
%   to a value of 1 or 2. We still apply this limitation, as high precision textures consume
%   more memory and other resources and are incompatible with very old graphics
%   hardware.
%
%   This is just a convenience shortcut for Screen('ColorRange', win, 1, 0, applyAlsoToMakeTexture);
%   with the added benefit of allowing to specify the background clear
%   color in normalized 0-1 range as well. This command is implied by use
%   of any of the high precision display device drivers (for attenuators,
%   Bits+ box etc.). It is only needed if you want to create the same
%   visual results on a 8 bit standard framebuffer without needing to
%   change your code, or if you want to set the 'applyAlsoToMakeTexture' flag to a
%   setting of non-zero, so unit colorrange also applies to Screen('MakeTexture').
%
%
% * 'StereoCrosstalkReduction' If a stereoMode is active or requested,
%   apply a shader first in the processing chain that for each eye aims to
%   reduce crosstalk from the other eye.
%
%   Usage:
%
%     PsychImaging('AddTask', 'LeftView', 'StereoCrosstalkReduction', method, crossTalkGain);
%     PsychImaging('AddTask', 'RightView', 'StereoCrosstalkReduction', method, crossTalkGain);
%
%   The 'method' parameter selects the method to use for crosstalk
%   reduction.
%
%   Currently only a method named 'SubtractOther' is implemented, which works as follows:
%
%   To reduce crosstalk, the contrast in the image of each eye, i.e., the
%   difference in color from the background level provided as background
%   clear color of the window is subtracted from the image of the other eye,
%   after scaling the contrast by 'crossTalkGain'. 'crossTalkGain' can be a
%   scalar, or a separate gain for each RGB channel. The background color
%   can be a scalar in the range 0-1, or a 3-element array to set the
%   backgroundlevel for each RGB channel separately. The background
%   color level should not be zero, as contrast then can't be inverted
%   around the background level. In general, the background level
%   should be high enough to allow unclamped inversion of the highest
%   contrast features of your stimulus at your 'crossTalkGain', or
%   artifacts will occur.
%
%
% * 'DisplayColorCorrection' Select a method for color correction to apply to
%   stimuli before output conversion and display. You have to specify a
%   color correction method 'methodname' to apply as parameter, see "help
%   PsychColorCorrection" for an overview of supported color correction
%   methods and their adjustable parameters. The imaging pipeline will be
%   set up to support the chosen color correction method. After you've
%   opened the onscreen window, you can use the different subcommands of
%   PsychColorCorrection() to change parameters of the color correction
%   algorithm at runtime.
%
%   Usage: PsychImaging('AddTask', whichView, 'DisplayColorCorrection', methodname);
%
%   Example: PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
%   This would apply a simple power-law gamma correction to all view
%   channels of a stereo setup, or the single view of a monoscopic setup.
%   Later on you could use the methods of PsychColorCorrection() to
%   actually set the wanted gamma correction factors.
%
%   Please note that we use the channel 'FinalFormatting' instead of
%   'AllViews' as we'd usually do. Both specs will work, but a selection
%   of 'FinalFormatting' will lead to faster processing in many cases, so
%   this is preferred here if you want to apply the same setting to all
%   view channels - or to a single monoscopic display. Should you find
%   that things don't work as expected, you might try 'AllViews' instead
%   of 'FinalFormatting' - There are subtle differences in how they
%   process your instructions, which may matter in some corner cases.
%
%
% * 'EnablePseudoGrayOutput' Enable the high-performance driver for the
%   rendering of up to 1786 different levels of gray on a standard - but
%   well calibrated - color monitor and 8 bit graphics card. This is done
%   by applying an algorithn known as "Pseudo-Gray" or "Bit stealing".
%   Selecting this mode implies use of 32 bit floating point
%   framebuffers, unless you specify use of a 16 bit floating point
%   framebuffer via 'FloatingPoint16Bit' explicitely. The expected range
%   of input luminance values is between 0 and 1. See "help CreatePseudoGrayLUT"
%   for further explanation.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput');
%
%
% * 'EnableGenericHighPrecisionLuminanceOutput'
%   Setup Psychtoolbox for conversion of high precision luminance images
%   into a format suitable for special high precision luminance display
%   devices. This is a generic support routine that uses LUT based
%   conversion.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableGenericHighPrecisionLuminanceOutput', lut);
%
%
% * 'EnableVideoSwitcherSimpleLuminanceOutput'
%   Setup Psychtoolbox for conversion of high precision luminance images
%   into a format suitable for driving the "VideoSwitcher" high precision
%   luminance display device which was developed by Xiangrui Li et al.
%
%   This implements the simple converter, which only needs the
%   Blue-To-Red-Ratio of the device as input parameter and performs
%   conversion via a closed-form formula without any need for lookup
%   tables. This is supposed to be fast.
%
%   See "help VideoSwitcher" for more info about the device and its
%   options.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableVideoSwitcherSimpleLuminanceOutput' [, btrr] [, trigger]);
%
%   - The optional 'btrr' parameter is the Blue-To-Red-Ratio to use. If the
%   parameter is left out, the btrr value will be read from a global
%   configuration file.
%
%   - The optional 'trigger' parameter can be zero for "No trigger", or 1
%   for "Use trigger as configured". By default, trigger is off (==0).
%   Enabled, one can use the VideoSwitcher('SetTrigger', ...); function to
%   configure when and how a trigger signal should be emitted. Trigger
%   signals are simply specific pixel patterns in the green output channel.
%   That channel is recognized by the VideoSwitcher as a trigger signal
%   control channel.
%
%
% * 'EnableVideoSwitcherCalibratedLuminanceOutput'
%   Setup Psychtoolbox for conversion of high precision luminance images
%   into a format suitable for driving the "VideoSwitcher" high precision
%   luminance display device which was developed by Xiangrui Li et al.
%
%   This implements the simple converter, which only needs the
%   Blue-To-Red-Ratio of the device as input parameter and performs
%   conversion via a closed-form formula without any need for lookup
%   tables. This is supposed to be fast.
%
%   See "help VideoSwitcher" for more info about the device and its
%   options.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableVideoSwitcherCalibratedLuminanceOutput' [, btrr] [, lut] [, trigger]);
%
%   - The optional 'btrr' parameter is the Blue-To-Red-Ratio to use. If the
%   parameter is left out, the btrr value will be read from a global
%   configuration file.
%
%   - The optional 'lut' paramter is a 257 elements vector of luminance
%   values, which maps blue channel drive indices to luminance values. This
%   lut needs to be acquired via a calibration procedure by use of a
%   photometer. If 'lut' is left out, the table will be read from a global
%   configuration file.
%
%   - The optional 'trigger' parameter can be zero for "No trigger", or 1
%   for "Use trigger as configured". By default, trigger is off (==0).
%   Enabled, one can use the VideoSwitcher('SetTrigger', ...); function to
%   configure when and how a trigger signal should be emitted. Trigger
%   signals are simply specific pixel patterns in the green output channel.
%   That channel is recognized by the VideoSwitcher as a trigger signal
%   control channel.
%
%
% * 'EnableNative10BitFramebuffer' Enable support for output of stimuli
%   with 10 bit precision per color channel (10 bpc / 30 bpp / "Deep color")
%   on graphics hardware that supports native 10 bpc framebuffers.
%
%   Under Linux, all AMD graphics cards since at least 2007, NVidia graphics cards
%   since 2008, and Intel graphics chips since at least 2010 do support native
%   10 bit framebuffers. Intel graphics chips must use the X11 "intel" video driver
%   to output their 10 bit framebuffers with actual 10 bit precision, the alternative
%   "modesetting" video driver does not support output with more than 8 bit yet.
%   XOrgConfCreator will take care of this Intel quirk when creating a custom xorg.config
%   for such 10 bpc setups under Intel.
%
%   Under MS-Windows, many graphics cards of the professional class AMD/ATI Fire/Pro
%   series (2008 models and later), and all current models of the professional class
%   NVidia Quadro series (2008 models and later), do support 10 bpc framebuffers
%   under some circumstances. 10 bpc display on classic CRT monitors which are connected
%   via analog VGA outputs is supported. Support for digital display devices like
%   LCD/OLED panels or video projectors depends on the specific type of display output
%   connector used, the specific panels, and their video settings. Consult manufacturer
%   documentation for details. In general, 10 bpc output is usually supported on graphics
%   cards and displays via DisplayPort or HDMI video outputs, but to our knowledge not
%   via DVI-D outputs.
%
%   If such a combination of graphics card and display is present on your system
%   on Linux or Microsoft Windows, then Psychtoolbox will request native support
%   from the standard graphics drivers, ie., it won't need to use our own
%   homegrown experimental box of tricks to enable this. You do need to enable/
%   unlock 10 bpc mode somewhere in the display driver settings though. On Linux you
%   can do this for supported cards and drivers via XOrgConfCreator + XOrgConfSelector,
%   on Windows the method is vendor specific.
%
%   Apple OSX, since version 10.11.2 "El Capitan", does support native 10 bpc video
%   output on some small subset of Apple hardware, as of May 2017 these are the MacPro
%   2013 "with some suitable displays" (Apple quotation), and the 27 inch iMac models
%   late 2014 and late 2015 with Retina 5k displays. We've confirmed this to be working
%   on the iMac 5k Retina 27 inch late 2014 model via photometer measurements. On OSX,
%   the OS will actually provide a 16 bit half-float framebuffer for our onscreen windows.
%   This buffer provides ~11 bpc effective linear precision in the displayable color
%   intensity range 0.0-1.0. The OS outputs this 11 bpc framebuffer as a native 10 bpc
%   video signal on suitable displays and uses some Apple proprietary software spatial
%   dithering algorithm to add 1 extra bit of simulated precision, so a photometer would
%   measure up to 11 bpc perceived/measured precision. On some other Mac models, which are
%   not in Apples list of 10 bit capable Macs, Apple uses a proprietary spatial dithering
%   algorithm implemented in software to fake a precision of 11 bpc on standard 8 bpc
%   framebuffers and displays, at least convincing enough for a photometer. The downside
%   of this proprietary dithering scheme is that visual stimulus onset timing precision
%   and timestamping precision is impaired, so Macs in "10 bit" framebuffer mode are not
%   suitable if trustworthy frame accurate visual timing is needed. Nothing we could do
%   about this. To summarize: EnableNative10BitFramebuffer mode on OSX will actually give
%   you a simulated 11 bpc framebuffer on some Mac hardware, plus severe visual timing
%   problems.
%
%   Psychtoolbox experimental 10 bpc framebuffer support:
%
%   Additionally, we support ATI/AMD Radeon hardware of the X1000, HD2000 - HD8000,
%   series and later models (everything since at least the year 2006) up to the AMD
%   Polaris gpu family under Linux via our own low-level setup mechanisms. AMD gpu's
%   of type Vega or more recent are not supported by our own hacks, but by the builtin
%   native driver support mentioned in the paragraphs above. These graphics cards support
%   a native ARGB2101010 framebuffer, ie., a system framebuffer with 2 bits for the alpha
%   channel, and 10 bits per color channel.
%
%   As this is supported by the hardware, but not always by the standard AMD
%   graphics drivers, we follow a hybrid approach: We use special low-level
%   hardware access to reconfigure the hardware for 10 bpc framebuffer support.
%   Then we use a special imaging pipeline formatting plugin to convert 16 bpc or
%   32 bpc stimuli into the special data format required by this framebuffer
%   configuration.
%
%   On Linux you must have run PsychLinuxConfiguration at least once on your
%   system at some point. You'll need to have one of the supported AMD Radeon
%   gfx-cards (see above) for this to work. If you use Linux with the free and
%   open-source AMD graphics drivers, 10 bpc framebuffer support should work
%   reliably, so use of the open-source drivers on Linux is recommended for
%   reliable results.
%
%   Getting a 10 bpc framebuffer working is only the first half of what you need for
%   high color precision output. Your graphics card must also be able to transmit the
%   video signal at high precision to the display device and the display must be able
%   to faithfully reproduce the high precision image. 10 bpc output has been verified
%   to work for analog VGA connected CRT monitors and displays on both AMD and NVidia
%   graphics cards which do support 10 bpc framebuffers, so with a analog VGA CRT you
%   should be safe. Note that this only applies to native VGA output via VGA connectors
%   or passive DVI-I to VGA adapters connected to a DVI-I connector. Active DisplayPort
%   to VGA adapters or active HDMI to VGA adapters may be limited to maximum 8 bpc output.
%   The status of 10 bpc output to digital display devices differs a lot across devices
%   and OS'es. Output of 10 bpc framebuffers to standard 8 bpc digital laptop panels or
%   DVI/HDMI/DisplayPort panels via digital dithering is known to work, but that is not
%   the real thing, only a simulation of 10 bpc via dithering to 8 bpc. This may or may
%   not be good enough for your specific visual stimulation paradigm. On a DVI-D connected
%   standard digital display, this dithered output is the best you will ever get.
%
%   DisplayPort: Recent AMD/Intel/NVidia graphics cards can output to some suitable DisplayPort
%   displays with 10 bpc or higher precision on Linux, and maybe also on MS-Windows, but you
%   have to verify this carefully for your specific display.
%
%   HDMI: Recent Intel graphics cards can output up to 12 bpc precision to HDMI deep color
%   capable displays on Linux if the X11 intel-ddx video driver is used. All AMD graphics
%   cards of model Radeon HD-5000 or later (and equivalent Fire-Series models) can output to
%   HDMI deep color capable displays with 10 bpc real precision at least if you use a Linux
%   kernel of version 3.16 or later with the open-source AMD graphics drivers. Execute
%   PsychLinuxConfiguration to enable this >= 10 bpc deep color output mode, then reboot your
%   machine once to enable it.
%
%   The status with the proprietary AMD drivers on MS-Windows is unknown.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableNative10BitFramebuffer' [, disableDithering=0]);
%
%   This function will setup a 32 bpc floating point resolution framebuffer by
%   default for Psychtoolbox drawing and stimulus processing. Output will happen
%   into a 10 bpc framebuffer. The function will also disable the graphics cards
%   gamma tables, so you'll need to use PsychImaging(...'DisplayColorCorrection'...)
%   for gamma and color correction if you need this.
%
%   The function will *not* disable dithering on digital displays by default,
%   but leave that decision to the operating system and graphics drivers of
%   your machine. A well working OS would disable dithering on a 10 bpc or
%   higher color depth display, if the display reports its capability to the
%   OS via its EDID info. It would enable dithering on < 10 bpc displays, so
%   you'd get a "pseudo 10 bpc" display where 10 bpc color depths is simulated
%   on a 6 bpc or 8 bpc display via the dithering.
%
%   You can disable dithering manually on some graphics cards by providing the
%   optional 'disableDithering' flag as 1. Currently mostly AMD cards allow this
%   control. NVidia or Intel cards require manual setup to force dithering off.
%
%
% * 'EnableNative11BitFramebuffer' Enable support for output of stimuli with (almost)
%   11 bit precision per color channel (11 bpc / 32 bpp / "Deep color") on graphics
%   hardware that supports native 11 bpc framebuffers. This will request an ~ 11 bpc
%   framebuffer from the operating system. If it can't get such a framebuffer on Linux
%   with AMD graphics hardware, it will use our own homegrown setup code to provide
%   such a framebuffer anyway on Radeon X1000, HD-2000 and later graphics cards and
%   equivalent Fire-Series graphics cards. On OSX 10.11.2 it will request and get a
%   ~11 bpc framebuffer on some Mac models. See the explanations above for 10 bpc on
%   OSX.
%
%   Read all the explanations in the section above for 'EnableNative10BitFramebuffer'
%   for capabilities, limitations and possible caveats on different systems.
%
%   Please note that on Linux this "11 Bit framebuffer" is not quite 11 bpc precision,
%   but only about ~ 10.6666 bpc precision. Specifically, the framebuffer can only
%   store at most 32 bits of color information per pixel, so it will store 11 bit
%   precision for the red channel (2048 distinct red intensity levels), 11 bit
%   (2048 levels) for the green channel, but only 10 bit (1024 levels) for the blue
%   channel, for a total number of 11 + 11 + 10 bits = 32 bits of color information
%   per pixel, or 4 billion different possible colors. A true 11 bpc framebuffer would
%   need 33 bits per pixel, and current graphics hardware can't handle that.
%
%   How many bits of precision of these ~ 11 bpc actually reach your display device?
%
%   - Analog VGA only provides for maximum 10 bpc output precision on all shipping
%     NVidia and AMD graphics cards at best. Intel graphics cards only allow for 8 bpc.
%
%   - DisplayPort or HDMI might allow for transfer of 11 bpc precision, in general they
%     support up to 12 bpc. However additional hardware restrictions for your graphics
%     card may limit precision to as low as 10 bpc. To our knowledge, only AMD graphics
%     cards support ~ 11 bpc framebuffers at all. Radeon HD-7000 and earlier can only
%     truly process up to 10 bpc, so 'EnableNative11BitFramebuffer' may not gain you any
%     precision over 'EnableNative10BitFramebuffer' in practice on these cards. AMD cards
%     of the "Sea Islands" family or later, mostly models from the year >= 2014, are
%     able to process and output up to 12 bpc over HDMI or DisplayPort, so they can
%     output true ~11 bpc images.
%
%   So obviously: Measure very carefully on your setup what kind of precision you really
%   get and make sure not to be fooled by dithering if you need precise low-level control
%   of spatial stimulus properties, or per-pixel high precision, instead of just averaged
%   over larger clusters of pixels.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableNative11BitFramebuffer' [, disableDithering=0]);
%
%
% * 'EnableNative16BitFramebuffer' Enable 16 bpc, 64 bpp framebuffer on some supported setups.
%
%   This asks to enable a framebuffer with a color depth of 16 bpc for up to 65535 levels of intensity
%   per red, green and blue channel, ie. 48 bits = 2^48 different colors. Currently, as of November 2021,
%   this mode of operation is only supported on Linux/X11 when using the open-source amdgpu-kms driver
%   on modern AMD GCN 1.1+ graphics cards [3]. On suitable setups, this will establish a 16 bpc framebuffer
%   which packs 3 * 16 bpc = 48 bit color info into 64 bpp pixels and the gpu's display engine will scan
%   out that framebuffer at 16 bpc. However, effective output precision is further limited to < 16 bpc by
%   your display, video cable and specific model of graphics card. As of November 2021, the maximum effective
%   output precision is limited to at most 12 bpc (= 4096 levels of red, green and blue) by the graphics card,
%   and this precision is only attainable on AMD graphics cards of the so called "Sea Islands" (cik) family
%   (aka GraphicsCore Next GCN 1.1 or later) or any later models when used with the amdgpu-kms display driver.
%
%   Older AMD cards of type GCN 1.0 "Southern Islands" or earlier won't work, as they only support at most 10
%   bpc overall output precision.
%
%   Please note that actual 12 bpc output precision can only be attained on certain display devices and
%   software configurations. As of November 2021, the following hardware + software combos have been
%   verified with a CRS ColorCal2 colorimeter to provide 12 bpc per color channel precision:
%
%   - The Apple MacBookPro 2017 15 inch with its built-in 10 bit Retina display, running under Ubuntu Linux
%     20.04 with Linux 5.14, as well as with a HDR-10 monitor via DisplayPort and also via HDMI. As those
%     displays are 10 bit only, the 12 bit precision was attained via spatial dithering by the gpu.
%
%   - AMD Ryzen 2400G with AMD RavenRidge integrated graphics chip with a HDR-10 monitor via DisplayPort and
%     also via HDMI. As that display is 10 bit only, the 12 bit precision was attained via spatial dithering
%     by the gpu.
%
%   High bit depth 12 bpc output should work with native 12 bit displays connected via HDMI or DisplayPort.
%   Modern AMD gpus, e.g., AMD Polaris, are also able to emulate 12 bpc precision via spatial dithering on
%   8 and 10 bit display panels. Measure your results carefully with a photometer or colorimeter, this is the
%   only way to verify your specific setup really achieves 12 bpc! See the sections about 11 bpc and 10 bpc
%   native framebuffers above for further details.
%
%   Required manual one-time setup:
%
%   1. You will need to install the AMD developed and maintained AMDVLK open-source vulkan driver from:
%
%      https://github.com/GPUOpen-Drivers/AMDVLK/releases
%
%      Driver release v-2021.Q4.1 or later from November 2021 is required to support this feature.
%      Note that release 2023-Q3.3 from September 2023 was the last release to support pre-Navi
%      gpu's like Polaris and Vega. Later versions only support AMD Navi and later with RDNA
%      graphics architecture.
%
%   2. You will need to install Linux kernel 5.14, which is currently not shipping in any Ubuntu release,
%      as of November 2021. A way to manually install it on Ubuntu 20.04-LTS is described on the following
%      web page via the "mainline" helper software:
%
%      https://ubuntuhandbook.org/index.php/2020/08/mainline-install-latest-kernel-ubuntu-linux-mint
%
%      Ubuntu 22.04-LTS should ship with a suitable kernel by default in April 2022.
%
%   3. If you are using an AMD Polaris gpu or later then you are done.
%
%      If you are using an old "Sea Islands" / "GraphicsCore Next 1.1" / "GCN 1.1" gpu, you must reboot
%      your machine with Linux kernel boot settings that select the new amdgpu kms driver and AMD DisplayCore,
%      instead of the old radeon kms driver that would be used by default. This requires adding the following
%      parameters to the kernel boot command line: "radeon.cik_support=0 amdgpu.cik_support=1 amdgpu.dc=1"
%
%      Additionally you would need a custom amdvlk driver, as AMD's current official AMDVLK driver does
%      no longer support pre-Polaris gpu's. We won't provide this driver for free at the moment, so please
%      enquire for potential paid support options on the Psychtoolbox user forum.
%
%   On AMD you can verify actual output bit depth for an output XXX by typing this command into a terminal
%   window, assuming your AMD graphics card is the first or only gpu in the system, ie. has index 0:
%
%   sudo cat /sys/kernel/debug/dri/0/XXX/output_bpc
%
%   E.g., for the internal laptop eDP flat panel of the MacBookPro 2017:
%
%   sudo cat /sys/kernel/debug/dri/0/eDP-1/output_bpc
%
%   Once the above one-time setup is done, adding the following command to your script will enable the
%   16 bpc framebuffer with up to 12 bpc effective output precision:
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableNative16BitFramebuffer');
%
%
% * 'EnableNative16BitFloatingPointFramebuffer' Enable support for output of
%   stimuli with up to 16 bit floating point precision per color channel on
%   graphics hardware and displays that support native 16 bpc floating point
%   framebuffers. Please note that the effective linear output precision of a
%   16 bit non-linear floating point framebuffer in the normalized range 0.0 - 1.0
%   (the "typical" output range for SDR standard dynamic range displays) is only
%   about 11 bits ~ 2048 levels of red, green, blue intensity. Specifically, for
%   normalized color values above 0.5 (ie. 50% max intensity) only 11 bpc can be
%   attained at most, whereas values smaller or equal 0.5 may be able to attain
%   12 bpc output precision on some recent gpu's. This is due to the precision of
%   floating point 16, which decreases with increasing magnitude of values and
%   maxes out at 14 bits for very small magnitude values - however current hardware
%   can not output at more than 12 bpc precision anyway, even for very small color
%   values.
%
%   Linux supports this mode on AMD gpu's of the "Sea Islands" gpu family or later
%   if you install the AMDVLK AMD open-source Vulkan driver and use the imaging
%   pipeline tasks 'UseVulkanDisplay' or 'EnableHDR'. See AMDVLK instructions
%   above for 'EnableNative16BitFramebuffer'. You need Linux 5.8 for AMD Polaris
%   and later (e.g. contained in Ubuntu 20.10 or 20.04.2-LTS), or Linux 5.12
%   (e.g., from Ubuntu 21.10 or 20.04.4-LTS) for also supporting the earlier AMD
%   "Sea Islands" and later gpu's.
%
%   macOS OpenGL does support this mode with what seems to be mostly software
%   rendering on most machines, ie. with very low performance and even worse timing.
%   Your mileage may vary.
%
%   Windows usually only supports this mode on some gpu's with Vulkan when HDR
%   mode is used via 'EnableHDR', but not in SDR standard dynamic range mode.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableNative16BitFloatingPointFramebuffer');
%
%
% * 'UseDataPixx' Tell Psychtoolbox that additional functionality for
%   displaying the onscreen window on a VPixx Technologies DataPixx device
%   should be enabled.
%
%   This command is implied by enabling a DataPixx video mode by one of the
%   commands for the DataPixx in the following sections.
%
%   'UseDataPixx' mostly prepares use of a variety of subfunctions in the
%   DataPixxToolbox ("help DataPixxToolbox") and in the PsychDataPixx()
%   high-level driver ("help PsychDataPixx").
%
%
% * 'EnableDataPixxL48Output' Setup Psychtoolbox for L48 mode of the VPixx
%   Technologies DataPixx device. This loads the graphics hardwares gamma
%   table with an identity mapping so it can't interfere with DPixx video
%   processing. It also sets up automatic generation of control signals to
%   support the features of DPixx that are available via the functions in
%   PsychDataPixx(). You will be able to upload new CLUT's into the DPixx
%   by use of the Screen('LoadNormalizedGammaTable', window, clut, 2);
%   command. CLUT updates will be synchronized with Screen('Flip') commands.
%   Please note that while L48 CLUT mode works even with very old
%   graphics hardware, this is a pretty cumbersome way of driving the
%   DPixx. On recent hardware, you will want to use M16 or C48 mode
%   (see below). That allows to draw arbitrarily complex stimuli with as
%   many colors as you want and PTB will take care of conversion into the
%   M16 or C48 format for DataPixx.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableDataPixxL48Output');
%
%
% * 'EnableDataPixxM16Output' Enable the high-performance driver for M16
%   mode of the VPixx Technologies DataPixx device. This is the fastest and
%   most elegant way of driving the DPixx box with 16 bit luminance output
%   precision. See "help DataPixx" for more information. Selecting this
%   mode implies use of 32 bit floating point framebuffers, unless you
%   specify use of a 16 bit floating point framebuffer via
%   'FloatingPoint16Bit' explicitely. If you do that, you will not be able
%   to use the full 16 bit output precision, but only approximately 10 bits.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableDataPixxM16Output');
%
%   If you want to make use of the color overlay plane in M16 mode, then
%   call the function like this:
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay');
%   See the explanation of color overlays in the section
%   'EnableBits++Mono++OutputWithOverlay' - behaviour of color overlays is
%   identical for the CRS Bits++ and the VPixx DataPixx.
%
%
% * 'EnableDataPixxC48Output' Enable the high-performance driver for the
%   C48 mode of VPixx technologies DataPixx box. This is the fastest and
%   most elegant way of driving the DataPixx box with 16 bit per color
%   channel output precision. See "help DataPixx" for more information.
%   Selecting this mode implies use of 32 bit floating point framebuffers,
%   unless you specify use of a 16 bit floating point framebuffer via
%   'FloatingPoint16Bit' explicitely. If you do that, you will not be able
%   to use the full 16 bit output precision, but only approximately 10 bits.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', mode);
%
%   See the section below about 'EnableBits++Color++Output' for the meaning
%   of the mandatory "mode" parameter.
%
%   You can use the RemapMouse() function to correct GetMouse() positions
%   for potential geometric distortions introduced by this function for
%   "mode" zero.
%
%
% * 'UseBits#' Tell Psychtoolbox that additional functionality for
%   displaying the onscreen window on a Cambridge Research Systems Bits#
%   device should be enabled.
%
%   This command is implied by enabling a Bits+ or Bits# video mode by one
%   of the commands for the Bits+/Bits# in the following sections, if the
%   driver can auto-detect a connected Bits# device. If it cannot auto-detect
%   a connected Bits# device and this command is omitted, Psychtoolbox will
%   instead assume that an older Bits+ is in use and only allow functionality
%   common to Bits# and Bits+, without automatic video mode switching.
%
%   If you provide this command, you can optionally specify the name of the
%   serial port to which your Bits# is connected, instead of leaving it to
%   the system to find this out (either via configuration file or via a
%   guess-o-matic).
%
%   Usage: PsychImaging('AddTask', 'General', 'UseBits#' [, BitsSharpSerialPort]);
%
%   'BitsSharpSerialPort' is optional and can be set to the name of a serial
%   port for your specific operating system and computer, to which the Bits#
%   is connected. If omitted, Psychtoolbox will look for the name in the first
%   line of text of a text file stored under the filesystem path and filename
%   [PsychtoolboxConfigDir 'BitsSharpConfig.txt']. If that file is empty, the
%   serial port is auto-detected (Good luck!).
%
%   'UseBits#' mostly prepares use of a variety of new Bits# subfunctions
%   in the BitsPlusPlus() high-level driver ("help BitsPlusPlus").
%
%
% * 'EnableBits++Bits++Output' Setup Psychtoolbox for Bits++ mode of the
%   Cambridge Research Systems Bits++ box. This loads the graphics
%   hardwares gamma table with an identity mapping so it can't interfere
%   with Bits++ T-Lock system. It also sets up automatic generation of
%   Bits++ T-Lock codes: You will be able to upload new CLUT's into the
%   Bits++ by use of the Screen('LoadNormalizedGammaTable', window, clut, 2);
%   command. CLUT updates will be synchronized with Screen('Flip')
%   commands, because PTB will generate and draw the proper T-Lock code
%   into the top line of your onscreen window. Please note that while
%   Bits++ CLUT mode works even with very old graphics hardware, this is a
%   pretty cumbersome way of driving the Bits++. On recent hardware, you
%   will want to use Mono++ or Color++ mode (see below). That allows to
%   draw arbitrarily complex stimuli with as many colors as you want and
%   PTB will take care of conversion into the Color++ or Mono++ format for
%   Bits++.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableBits++Bits++Output');
%
%
% * 'EnableBits++Mono++Output' Enable the high-performance driver for the
%   Mono++ mode of Cambridge Research Systems Bits++ box. This is the
%   fastest and most elegant way of driving the Bits++ box with 14 bit
%   luminance output precision. See "help BitsPlusPlus" for more
%   information. Selecting this mode implies use of 32 bit floating point
%   framebuffers, unless you specify use of a 16 bit floating point
%   framebuffer via 'FloatingPoint16Bit' explicitely. If you do that, you
%   will not be able to use the full 14 bit output precision of Bits++, but
%   only approximately 10 bits.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableBits++Mono++Output');
%
%   If you want to make use of the color overlay plane in Mono++ mode, then
%   call the function like this:
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableBits++Mono++OutputWithOverlay');
%
%   Then you can query the window handle of the overlay window via:
%
%   overlayWin = PsychImaging('GetOverlayWindow', window);
%
%   'overlayWin' is the handle to the overlay window associated with the
%   overlay of onscreen window 'window'. The overlay window is a standard
%   offscreen window, so you can do anything with it that you would want to
%   do with offscreen windows. The only difference is that the window is a
%   pure index window: It only has one "color channel", which can be written
%   with color values between 0 and 255. Values 1 to 255 get mapped to the
%   corresponding color indices of the Bits++ overlay plane: A zero value is
%   transparent -- Content of the onscreen window is visible. Positive
%   non-zero color values map to the 255 indices available in overlay mode,
%   these get mapped by the Bits++ CLUT to colors. You can define the
%   mapping of indices to CLUT colors via the
%   Screen('LoadNormalizedGammaTable', window, clut, 2); command.
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
% * 'EnableBits++Color++Output' Enable the high-performance driver for the
%   Color++ mode of Cambridge Research Systems Bits++ box. This is the
%   fastest and most elegant way of driving the Bits++ box with 14 bit
%   per color channel output precision. See "help BitsPlusPlus" for more
%   information. Selecting this mode implies use of 32 bit floating point
%   framebuffers, unless you specify use of a 16 bit floating point
%   framebuffer via 'FloatingPoint16Bit' explicitely. If you do that, you
%   will not be able to use the full 14 bit output precision of Bits++, but
%   only approximately 10 bits.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', mode);
%
%   "mode" is a mandatory numeric parameter which must be 0, 1 or 2. In
%   Color++ mode, the effective horizontal display resolution is only half
%   the normal horizontal resolution. To cope with this, multiple different
%   methods are implemented to squeeze your stimulus image horizontally by
%   a factor of two. The following options exist:
%
%   0 = This is the "classic" mode which was used in all Psychtoolbox
%   versions prior to 22nd September 2010. If you want to keep old code
%   working as is, select 0. In this mode, your script will only see a
%   framebuffer that is half the true horizontal resolution of your
%   connected display screen. Each drawn pixel will be stretched to cover
%   two pixels on the output display device horizontally. While this
%   preserves the content of your stimulus image exactly, it means that the
%   aspect ratio of all displayed text and stimuli will be 2:1. Text will
%   be twice as wide as its height. Circles or squares will turn into
%   horizontal ellipses or rectangles etc. You'll need to do extra work in
%   your code if you want to preserve aspect ratio properly.
%
%   You can use the RemapMouse() function to correct GetMouse() positions
%   for potential geometric distortions introduced by this function for
%   "mode" zero.
%
%   Example: A fine vertical grid with alternating vertical white and black
%   lines would display as expected, but each white or black stripe would be
%   two pixels wide on the display instead of one pixel wide.
%
%   1 = Subsample: Your framebuffer will appear at the same resolution as
%   your display device. Aspect ratio of drawn stimuli/text etc. will be
%   correct and as expected. However, every 2nd column of pixels in your
%   stimulus (ie., all odd-numbered x-coordinates 1,3,5,7,...) will be
%   completely ignored, only even columns are used!
%
%   Example: A fine vertical grid with alternating vertical white and black
%   lines would display as a purely white image, as only the white pixels
%   in the even columns would be used, whereas the black pixels in the odd
%   columns would be ignored.
%
%   2 = Average: Your framebuffer will appear at the same resolution as
%   your display device. Aspect ratio of drawn stimuli/text etc. will be
%   correct and as expected. However, each pair of adjacent even/odd pixel
%   columns will be averaged before output. Stimulus pixels 0 and 1 will
%   contribute the mean color for display pixel 0. Pixels 2 and 3 will be
%   averaged into display pixel 1 and so on. Visually this gives the most
%   pleasing and smooth results, but if adjacent even/odd pixels don't have
%   the same color value, you'll obviously get an output color that is
%   neither the color of the even pixel nor the odd pixel, but the average
%   of both.
%
%   Example: A fine vertical grid with alternating vertical white and black
%   lines would display as a 50% gray image, as the alternating white and
%   black columns would be averaged into the average of white and black,
%   which is 50% gray.
%
%
% * 'EnableDualPipeHDROutput' Enable EXPERIMENTAL high-performance driver
%   for HDR display devices which are composites of two separate displays.
%
%   EXPERIMENTAL proof-of-concept code with no real function yet!
%
%   This is meant for high-precision luminance or color output. It implies
%   use of 32 bpc floating point framebuffers unless otherwise specified by
%   other calls to PsychImaging().
%
%   The pair of specially encoded output images that are derived from
%   content of the onscreen window shall be output to both, the display
%   associated with the screen given to PsychImaging('OpenWindow',...); and
%   on the screen with the index 'pipe1Screen', using appropriate encoding
%   to drive the HDR device or similar composite device.
%
%   Usage: PsychImaging('AddTask', 'General', 'EnableDualPipeHDROutput', pipe1Screen [, pipe1Rectangle]);
%
%   Optionally you can pass a 'pipe1Rectangle' if the window with the
%   pipe1 image shall not fill the whole 'pipe1Screen', but only a
%   subregion 'pipe1Rectangle'.
%
%
% * 'AddOffsetToImage' Add a constant color- or intensity offset to the
%   drawn image, possibly also multiply with a gain, prior to all following
%   image processing and post-processing operations:
%
%   Outimage(x,y) = (Inimage(x,y) + preOffset) * gain + Offset.
%
%   If the framebuffer is in a color display mode, the same offset will be added
%   to all three color channels.
%
%   Usage: PsychImaging('AddTask', whichView, 'AddOffsetToImage', Offset [, gain=1][, preOffset=0]);
%   'gain' and 'preOffset' are optional, 'Offset' is required.
%
%   Example: PsychImaging('AddTask', 'AllViews', 'AddOffsetToImage', 0.5);
%
%
% * 'MirrorDisplayTo2ndOutputHead' Mirror the content of the onscreen
%   window to the given 2nd screen, ie., to a 2nd video output (head) of a
%   multi-display setup. This should give roughly the same result as if one
%   switches the graphics card into "Mirror mode" or "Clone mode" via the
%   display settings panel of your operating system, but with different
%   tradeoffs. The disadvantage compared to use of an operating system provided
%   mirror mode is that it requires more gpu processing, adding one millisecond to
%   potentially multiple milliseconds of processing time to each Screen('Flip').
%
%   The advantage is that it allows for proper visual stimulation timing and for
%   higher animation framerates in various configurations in which the operating
%   system provided mechanisms would cause impaired timing and drastically reduced
%   or erratic framerates, especially when the subject visual stimulation display
%   and the experimenter feedback mirror display do not have perfectly synchronized
%   video refesh cycles, e.g., because the two displays are not the same model and
%   vendor, or can't be set to the same video resolution and refresh rate, or
%   can't be connected both to the same type of video link (e.g. not both via
%   DisplayPort), or operating system limitations prevent synchronized scanout.
%   Non-synchronized video refresh can cause broken timing, strong judder of
%   animations or severely reduced framerate, e.g., if a 60 Hz experimenter
%   display throttles all animations down to 60 fps, even though the subjects
%   stimulus display is configured for 120 Hz operation. This mode will run the
%   main window for stimulus presentation with proper timing and vsync, whereas
%   it runs the experimenter feedback "mirror" display without vsync, allowing
%   for tearing artifacts to happen, in exchange for optimal timing and performance
%   on the main display.
%
%   Another advantage is that this allows to display an overlay window on top of
%   the visual stimulus mirror image, to supply some additional debug information
%   to the experimenter. Operating system provided mirror modes don't allow for
%   such overlays.
%
%   This function only works for monoscopic displays and single-window stereo
%   modes (steremode 2 - 9 and some PsychImaging modes), it can not be used with
%   frame-sequential or dual-window / dual-stream stereo display modes, as it
%   internally uses the machinery for stereomode 10 with a few modifications to
%   get its job done, so obviously neither dual-window mode 10 nor any similar
%   mode can be used without interference.
%
%   This function is typically used for mirroring onto a secondary display of a
%   multi-display setup on Apple macOS, or if you need to mirror on Microsoft
%   Windows and can't use "desktop spanning" mode on Windows to achieve dual
%   display output. On Linux/X11, this can be useful in combination with the
%   Vulkan display backend, ie. the 'UseVulkanDisplay' task. On Linux or Windows
%   with a properly synchronized dual-display setup, rather use the task
%   'MirrorDisplayToSingleSplitWindow' (see below). That task should work
%   well on in synchronized mode, and is more efficient. It would also work well
%   in conjunction with a hardware display splitter attached to a single
%   head on any operating system. It has the advantage of consuming less
%   memory and compute ressources, so it is potentially faster or provides
%   a more reliable overall timing. On modern Linux distributions with X-Server
%   version 21 or later, e.g., Ubuntu 22.04-LTS or later, that more efficient
%   mode can also be made to work, after special setup with XOrgConfCreator, on
%   unsynchronized displays.
%
%   Usage: PsychImaging('AddTask', 'General', 'MirrorDisplayTo2ndOutputHead', mirrorScreen [, mirrorRectangle=[]][, specialFlags=0][, useOverlay=0]);
%
%   The content of the onscreen window shall be shown not only on the
%   display associated with the screen given to PsychImaging('OpenWindow',
%   ...); but also (as a copy) on the screen with the index 'mirrorScreen'.
%
%   Optionally you can pass a 'mirrorRectangle' if the window with the
%   mirror image shall not fill the whole 'mirrorScreen', but only a
%   subregion 'mirrorRectangle'.
%
%   'specialFlags' are optional flags to customize the mirror window, similar to
%   the specialFlags parameter of PsychImaging('OpenWindow', ...., specialFlags).
%   Cfe. flags like kPsychGUIWindow and kPsychGUIWindowWMPositioned to turn the
%   experimenter feedback mirror window into a regular desktop GUI window that
%   can be moved and resized.
%
%   Optionally you can set the 'useOverlay' flag to 1, to request use of an
%   overlay window on top of the mirror window. The function ...
%   overlaywin = PsychImaging('GetMirrorOverlayWindow', win);
%   ... will return a window handle overlaywin for a given onscreen window win,
%   and you can then use overlaywin for drawing content into that overlay.
%
% * 'MirrorDisplayToSingleSplitWindow' Mirror the content of the onscreen
%   window to the right half of the desktop (if desktop spanning on a
%   dual-display setup is enabled) or the right-half of the virtual screen
%   if a display splitter (e.g., Matrox Dualhead2Go (TM)) is attached to a
%   single head of a graphics card. This should give the same result as if one
%   switches the graphics card into "Mirror mode" or "Clone mode" via the
%   display settings panel of your operating system. It offers the same tradeoffs
%   and advantages as explained above for 'MirrorDisplayTo2ndOutputHead' mode.
%   The 'MirrorDisplayToSingleSplitWindow' task may be a bit more efficient than
%   'MirrorDisplayTo2ndOutputHead', but it requires either use of a modern Linux
%   distribution like Ubuntu 22.04-LTS with X-Server 21, and some configuration
%   with XOrgConfCreator, or the use of two display devices of identical
%   model from the same vendor, set to exactly the same video mode
%   (resolution and refresh rate), and identical video connections, so the
%   video refresh cycles of both the stimulus presentation display, and the
%   experimenter feedback mirror display are perfectly synchronized with
%   each other. Otherwise timing and performance will be less than optimal!
%   If you can't meet these conditions then the 'MirrorDisplayTo2ndOutputHead'
%   task is the better choice, at expense of higher gpu load.
%
%   Usage: PsychImaging('AddTask', 'General', 'MirrorDisplayToSingleSplitWindow' [, useOverlay=0][, mirrorDestination]);
%
%   Optionally you can set the 'useOverlay' flag to 1, to request use of an
%   overlay window on top of the mirrored stimulus. The function ...
%   overlaywin = PsychImaging('GetMirrorOverlayWindow', win);
%   ... will return a window handle overlaywin for a given onscreen window win,
%   and you can then use overlaywin for drawing content into that overlay.
%
%   Optionally you can provide the 'mirrorDestination' parameter to specify
%   where the mirror image should go within the onscreen window, and which
%   size it has. By default, the mirror image will fill out the right half
%   of the onscreen window, an appropriate choice if your mirror display
%   monitor has the same resolution as the stimulus display monitor and
%   both are arranged side-by-side in a dual-display setup. For different
%   monitor arrangements, e.g., triple-display setups or similar, of for
%   different selected mirror monitor resolutions, a future Psychtoolbox
%   release may also select a proper location and size for the mirror
%   image. But for now, if your mirror monitor has a different resolution
%   than the stimulus monitor, you need to manually specify the mirror
%   monitors resolution, ie. width x height in pixels as mirrorDestination
%   parameter, a two-component row vector of form [width, height]., so the
%   mirror image can get appropriately scaled to the resolution of the
%   mirror monitor. On Linux you will also need to use the SetResolution()
%   function to set a horizontal resolution that is twice the horizontal
%   resolution of your stimulus monitor(s). E.g., if your stimulus monitor
%   is 2560x1440 resolution and your mirror monitor is 1280x1024 pixels,
%   and both are attached to Psychtoolbox screen 1 (aka X-Screen 1), you'd
%   have to call SetResolution(1, 2*2560, 1440); before 'OpenWindow', and
%   specify 'mirrorDestination' as [1280, 1024]. A setup with three
%   monitors, two stimulus monitors next to each other at 2560x1440 pixels
%   each, plus one mirror monitor of 1280x1024 right to the both stimulus
%   monitors would instead require a call to
%   SetResolution(1, 2 * (2560 + 2560), 1440). Note that setups which
%   require specification of 'mirrorDestination' may not work at all at the
%   moment on systems other than Linux.
%
%   Optionally, if you don't need the imaging pipeline and don't need the overlay
%   for experimenter feedback, ie. you let 'useOverlay' = 0, you can add ...
%   PsychImaging('AddTask', 'General', 'DontUsePipelineIfPossible');
%   ... if you don't intend to use the imaging pipeline for anything else
%   than display mirroring. This will allow further optimizations.
%
%
% * 'RestrictProcessing' Restrict stimulus processing to a specific subarea
%   of the screen. If your visual stimulus only covers a subarea of the
%   display screen you can restrict PTB's output processing to that
%   subarea. This may save some computation time to allow for higher
%   display redraw rates.
%
%   Usage: PsychImaging('AddTask', whichChannel, 'RestrictProcessing', ROI);
%
%   ROI is a rectangle defining the area to process ROI = [left top right bottom];
%   E.g., ROI = [400 400 800 800] would only create output pixels in the
%   screen area with top-left corner (400,400) and bottom-right corner
%   (800, 800).
%
%
% * 'FlipHorizontal' and 'FlipVertical' flip your output images
%   horizontally (left- and right interchanged) or vertically (upside down).
%
%   Usage: PsychImaging('AddTask', whichChannel, 'FlipHorizontal');
%   Usage: PsychImaging('AddTask', whichChannel, 'FlipVertical');
%
%   You can use the RemapMouse() function to correct GetMouse() positions
%   for potential geometric distortions introduced by this function.
%
%
% * 'GeometryCorrection' Apply some geometric warping operation during
%   rendering of the final stimulus image to correct for geometric
%   distortion of your physical display device. You need to measure the
%   geometric distortion of your display with a suitable calibration
%   procedure, then compute an inverse warp transformation to undo this
%   distortion, then provide that transformation to this function.
%
%   Usage: PsychImaging('AddTask', whichChannel, 'GeometryCorrection', calibfilename [, debugoutput] [, arg1], [arg2], ...);
%
%   'calibfilename' is the filename of a calibration file which specified
%   the type of undistortion to apply. Calibration files can be created by
%   interactive calibration procedures. See 'help CreateDisplayWarp' for a
%   list of calibration methods. One of the supported procedures is, e.g.,
%   "DisplayUndistortionBezier", read "help DisplayUndistortionBezier". The
%   recommended method for most cases is 'DisplayUndistortionBVL', read
%   "help DisplayUndistortionBVL" for help.
%
%   The optional flag 'debugoutput' if set to non-zero value will trigger
%   some debug output about the calibration with some calibration methods.
%
%   The optional 'arg1', 'arg2', ..., are optional parameters whose
%   meaning depends on the calibration method in use.
%
%   Use of geometry correction will break the 1:1 correspondence between
%   framebuffer pixel locations (x,y) and the mouse cursor position, ie. a
%   mouse cursor positioned at display position (x,y) will be no longer
%   pointing to framebuffer pixel (x,y). If you want to know which
%   pixel in your original stimulus image corresponds to a specific
%   physical display pixel (or mouse cursor position), use the function
%   RemapMouse() to perform the neccessary coordinate transformation.
%
%
% * 'UseVRHMD' Display this onscreen window on a "Virtual Reality" head mounted
%   display (HMD), e.g., the Oculus Rift DK1 or Rift DK2. This enables display of
%   stereoscopic visual stimuli on supported virtual reality headsets.
%   You need to have the neccessary vendor supplied VR runtimes installed for
%   this to work.
%
%   Simple usage:
%
%   The most simple way to setup a HMD for use is to add a call to
%   hmd = PsychVRHMD('AutoSetupHMD') instead of a call to
%   PsychImaging('AddTask', 'General', UseVRHMD', ...). The 'AutoSetupHMD'
%   call would detect the first supported HMD device on your computer, connect to
%   it, then set it up with reasonable default operating parameters. Then it would
%   call this PsychImaging task to perform all required setup steps.
%
%   Advanced usage:
%
%     1. Open a connection to a HMD and get a handle for the device:
%        For example, if you wanted to use a Oculus Rift DK1 or DK2, you could
%        do:
%
%        hmd = PsychOculusVR('Open' ...);
%
%     2. Perform basic configuration of the HMD via the HMD specific driver.
%
%     3. Add a PsychImaging task for the HMD and pass in its device handle 'hmd':
%        PsychImaging('AddTask', 'General', 'UseVRHMD', hmd);
%
%   This sequence will perform the necessary setup of panel fitter, stereo display
%   mode and image post-processing for geometry correction, color aberration
%   correction and vignette correction for a fullscreen window on the HMD.
%
%
% * 'UseVulkanDisplay' Display this onscreen window using a Vulkan-based display
%   backend. This only works on graphics card + operating system combinations
%   which support both the OpenGL and Vulkan rendering api's and OpenGL-Vulkan
%   interop. As of March 2023 this would be modern AMD, Intel, and NVidia graphics
%   cards under modern GNU/Linux (Ubuntu 22.04-LTS and later) and under Microsoft
%   Windows-10 and later. Very experimental support for the curious is also possible
%   on the RaspberryPi 4/400 and later with Mesa 23.1 or later, but this support
%   does not have any advantages over standard mode, quite the contrary! See
%   'help RaspberryPiSetup' if you wanted to try it anyway.
%
%   Apple macOS 10.15.4 Catalina and later is also supported if you install the
%   3rd party Khronos open-source MoltenVK "Vulkan on Metal" driver. However, this
%   support is very unreliable wrt. stimulus timing at least on macOS 10/11/12, due
%   to severe bugs in Apple's macOS 10/11 Metal graphics api, which only Apple could
%   fix. Performance is also miserable, achieving at most a framerate that is half
%   the video refresh rate of the display monitor. The status on macOS 13 is unknown
%   as of March 2023. Do not ask us for help in using this on macOS, don't trust it!
%
%   At the moment 'UseVulkanDisplay' does not provide any advantages for standard
%   visual stimulus display tasks. It can be of some benefit if one employs display
%   mirroring on Linux/X11 via the 'MirrorDisplayTo2ndOutputHead' task, or wants to
%   flexibly operate multi-display setups under Linux/X11. It also provides very high
%   precision framebuffers under Linux with AMD graphics cards and the AMDVLK Vulkan
%   driver via the 'EnableNative16BitFramebuffer' task. For other more standard tasks,
%   it is of no benefit right now and may reduce performance instead.
%
%   Usage:
%
%   PsychImaging('AddTask', 'General', 'UseVulkanDisplay');
%
%   Psychtoolbox will try to display the onscreen window by using a Vulkan driver
%   with Vulkan/WSI backend, instead of the usual OpenGL windowing system backend.
%   This may fail if the given system setup does not support this.
%
%
% * 'EnableHDR' Display this onscreen window on a "High dynamic range" (HDR) display.
%   This requires a combination of operating-system, display drivers, graphics card,
%   video cables and display devices which are at least HDR-10 capable.
%
%   For hardware and system requirements, setup instructions, and further explanations
%   read "help PsychHDR".
%
%   Usage:
%
%   PsychImaging('AddTask', 'General', 'EnableHDR' [, unit='Nits'][, hdrMode='Auto'][, extraRequirements]);
%
%   Optional parameters:
%
%   'unit'  The unit in which color values are specified by the users drawing code.
%           Default value is 'Nits', ie. 1 unit = 1 Nit = 1 candela per square-meter.
%           '80Nits', ie. 1 unit = 80 Nits = 80 cd/sqm = Supposedly the SDR range.
%
%   'hdrMode' General mode of operation for HDR display:
%             Default is 'Auto' for auto-selection of optimal mode for given system
%             configuration, selected out of the following available op-modes:
%
%             'HDR10' Standard HDR-10, with 10 bpc color precision, ITU Rec 2020 input
%             color space (aka BT-2020), SMPTE ST-2084 PQ "Perceptual Quantizer" OETF
%             transfer function.
%
%             -> Currently 'Auto' will select 'HDR-10' as the only supported op-mode.
%
%   'extraRequirements' String with various keywords to specify special requirements.
%                       Default is empty, ie. no extra requirements. Currently supported
%                       keywords are:
%
%                       'Dummy' - Only simulate HDR on a SDR standard dynamic range display.
%                                 This only performs setup steps and processing possible on
%                                 a SDR display, to allow for basic script development and
%                                 testing. Visual results will be obviously wrong!
%
%
% * 'UseStaticHDRHack' Use a Linux + X11 only HDR display setup hack to allow multi-display
%   HDR presentation, e.g., for dual-HDR-monitors setups with two identical HDR monitors
%   or display devices for binocular or stereoscopic HDR visual stimulation.
%
%   Currently our Vulkan-based HDR display backend does not yet support proper dual-display/
%   multi-display stereo presentation in a proper way. As a stop-gap solution until we have
%   proper support implemented, this task requests use of a Linux + X11 specific hack to make
%   HDR stereo presentation work ok. The hack is a hack, not guaranteed to work on your setup,
%   or even if it works it is not guaranteed to work permanently, e.g. across operating system
%   updates. Use of this highly experimental hack is at your own risk, you have been warned!
%   Also note that using this PsychImaging tasks will make your script non-portable to other
%   operating systems than Linux or other display systems than X11 (e.g., Wayland) for the time
%   being. Another limitation is that this only works with 10 bpc HDR-10 precision framebuffers,
%   not with the higher precision 16 bit floating point framebuffers offered by standard HDR
%   mode. This hack will be limited permanently to HDR-10 with static HDR metadata and 10 bit
%   precision. It will eventually be superseded by a proper multi-display/stereo HDR implementation,
%   if and when we get funding to work on this.
%
%   To setup HDR-10 display with this hack:
%
%   0. You *must* install GNU/Octave and setup Psychtoolbox for use with GNU/Octave! Even if
%      you prefer to run your scripts in Matlab, Octave is needed as a background helper! In
%      general we recommend using Octave for this special hack, as Matlab's GUI is currently
%      incompatible with 10 bit display mode as of R2020b. It will misrender, crash or simply
%      hang! The only way to use Matlab here is without GUI in the terminal window by launching
%      Matlab via "matlab -nodesktop". In this case, using Octave in the first place may be
%      more convenient, as Octave's GUI works just fine!
%
%   1. Connect identical models of HDR monitors to your X-Screen 0, or at your leisure,
%      create a separate X-Screen for stimulation by use of XOrgConfCreator. It is important
%      that the X-Screen used for visual HDR stimulation only has HDR monitors connected, no
%      other monitors, iow. the same configuration as if you'd do multi-display stereo stimulation
%      for standard dynamic range stimuli under Linux.
%
%   2. Use XOrgConfCreator to create a separate X-Screen - or setup your one and only X-Screen 0
%      to use 10 bpc, color depth 30 bit mode. If the XOrgConfCreator asks you if you want to
%      setup advanced settings, answer "y"es, and then when it asks for 10 bpc / 30 bit deep
%      color support, answer "y"es and select all X-Screens which will be used for HDR stimulation.
%
%   3. Select this HDR-10 10 bit configuration via XOrgConfSelector, logout and login again to
%      get your single or multi-X-Screen setup with the X-Screen to which the HDR displays are
%      connected running in 10 bpc color depth.
%
%   4. Now you can use a regular HDR script, but add the PsychImaging('AddTask', 'General', 'UseStaticHDRHack')
%      task before opening the onscreen window. Also select a stereoMode of 4 for dual-display stereo
%      stimulation in PsychImaging('OpenWindow', ..., stereoMode);
%
%   5. All connected HDR monitors should switch to HDR mode and the onscreen window should be setup
%      for stereo stimulus display across HDR monitors, as usual.
%
%   If these steps don't work for you, too bad. If they work, great.
%
%   This mode will only allow slow switching of different HDR metadata settings via PsychHDR('HDRMetadata'),
%   and each switch of settings will be slow and accompanied by flicker. For this reason you can
%   also assign the meta-data as part of the 'UseStaticHDRHack' setup as optional hdrStaticMetadata
%   struct. The same HDR metadata will be applied to all connected HDR monitors for a window / X-Screen.
%   Other limitations may apply, and this functionality is so far only tested on one HDR setup.
%
%   Usage:
%
%   PsychImaging('AddTask', 'General', 'UseStaticHDRHack' [, hdrStaticMetadata]);
%
%
% * More actions will be supported in the future. If you can think of an
%   action of common interest not yet supported by this framework, please
%   file a feature request on our Wiki (Mainpage -> Feature Requests).
%
%
% After adding all wanted task specifications and other requirements,
% call...
%
% [windowPtr, windowRect] = PsychImaging('OpenWindow', screenid, [backgroundcolor], ....);
%
% - Finishes the setup phase for imaging pipeline, creates a suitable onscreen
% window and performs all remaining configuration steps. After this
% command, your onscreen window will be ready for drawing and display of
% stimuli. All specified imaging operations will get automatically applied
% to your stimulus before stimulus onset.
%
%
% After the window has been opened you can call the following commands any
% time at runtime:
%
% PsychImaging('RestrictProcessingToROI', window, whichChannel, ROI);
% - Restrict the processing area of viewChannel 'whichChannel' of onscreen
% window 'window' to the rectangular subarea defined by 'ROI'. See the
% explanation above for subtask 'RestrictProcessing'. This does exactly the
% same but allows a dynamic change of the restricted area at any point
% during your experiment script.
%
%
% PsychImaging('UnrestrictProcessing', window, whichChannel);
% - Remove a restriction of the processing area of viewChannel
% 'whichChannel' of onscreen window 'window' to a previously defined
% subarea. Can be called anytime during your scripts execution.
%
%
% [overlaywin, overlaywinRect] = PsychImaging('GetOverlayWindow', win);
% - Will return the handle to the 'overlaywin'dow associated with the
% given 'win'dow, if any. Will abort with an error message if the 'win'dow
% doesn't have an associated overylay window.
% Currently, only the CRS Bits+ box in Mono++ mode and the VPixx DataPixx
% box in M16 mode does support overlays. Other output drivers don't support
% such a feature. See "help BitsPlusPlus" for subfunction
% 'GetOverlayWindow' for more explanations of the purpose and properties of
% overlay windows. The explanations apply to the DPixx device as well if it
% is opened in videomode 'M16WithOverlay'.
%
%
% overlaywin = PsychImaging('GetMirrorOverlayWindow', win);
% - Will return the handle to the mirror 'overlaywin'dow associated with the
% given 'win'dow, if any. Will abort with an error message if the 'win'dow
% doesn't have an associated mirror overylay, because the 'win'dow was not
% configured for mirror mode, or for use of an overlay on its mirror window.
% Cfe. the PsychImaging stimulus mirroring tasks 'MirrorDisplayTo2ndOutputHead'
% and 'MirrorDisplayToSingleSplitWindow' for use cases.
%
%
%
% The following commands are only for specialists:
%
% [imagingMode, needStereomode] = PsychImaging('FinalizeConfiguration');
% - Finish the configuration phase for this window. This will compute an
% optimal configuration for all stages of the pipeline, but won't apply it
% yet. You'll have to call Screen('OpenWindow', windowPtr, ......,
% imagingMode, ...); with the returned 'imagingMode' + any other options
% you'd like to have for your window. After that, you'll have to call
% PsychImaging('PostConfiguration') to really apply and setup all your
% configuration settings. If you don't have unusual needs, you can simplify
% these steps by simply calling PsychImaging('OpenWindow', ....);
% with the same parameters that you'd pass to Screen('OpenWindow', ....);
% PsychImaging will perform all necessary steps to upon return, you'll have
% your window properly configured.
%
%
% PsychImaging('PostConfiguration', windowPtr [, clearcolor]);
% - To be called after opening the onscreen window 'windowPtr'.
% Performs all the setup work to be done after the window was created.
%
%

% Notes:
%
%
%

% History:
% 3.6.2007 Written. (MK)
%
% 19.7.2007 Added our own implementation mystrcmp() of Matlabs strcmp().
% The strcmp() of Octave is not as powerful as Matlab 7's strcmp(). We need
% a powerful implementation that works for both runtime environments, so we
% reimplemented in in a portable way.
%
% 19.7.2007   Added initial support for display geometry correction. (MK).
%
% 27.9.2007   Added support for floating point framebuffer, Bits++ and
%             Brightside-HDR. Documentation cleanup. (MK).
%
% 13.1.2008   Support for 10 bpc native framebuffer of ATI Radeons. (MK).
%
% 17.4.2008   Support for a few new subcommands, and description of overlay
%             planes setup with Bits++ in Mono++ mode. (MK).
%
% 18.5.2008   A few bug fixes and support for 'DisplayColorCorrection' setup
%             code: Now a central solution that will work for all current and
%             future output devices (hopefully). (MK).
%
% 02.7.2009   Add CLUT configuration support for ICM color correction (MK).
%
% 18.12.2009  Add support for VPixx Technologies DataPixx device for video
%             processing modes L48, M16, C48 and color overlays in M16
%             mode, mostly via calls into the new PsychDataPixx() driver.
%             Also support a new 'General' task 'UseDataPixx'. (MK)
% 04.03.2010  Bugfixes and workarounds to 'ColorCorrection' setup code. (MK)
%
% 26.04.2010  Disable workarounds from 04.03.2010, as Screen() is fixed now. (MK)
%
% 02.09.2010  Add support for 'InterleavedColumnStereo'- for auto-stereoscopic
%             displays, e.g., parallax barrier and lenticular sheets. (MK)
%
% 03.04.2011  Add support for 'EnableCLUTMapping' for old fashioned clut animation. (MK)
%
% 26.12.2011  Add support for ptb_geometry_inverseWarpMap inverse mapping
%             of geometry corrected displays. See RemapMouse.m. (MK)
%
% 11.09.2012  Add support for stereo sync line handling, and for
%             scanning backlight control of the ViewPixx in stereomode. (MK)
%
% 23.12.2012  Add support for 'SideBySideCompressedStereo' stereo mode. (MK)
%
% 23.12.2012  Add support for 'UsePanelFitter' Screen panelfitter setup. (MK)
%
% 15.04.2013  Add support for 'UseGPGPUCompute', currently via GPUmat. (MK)
%
% 03.07.2013  Call PsychJavaSwingCleanup via onscreen window close hook. (MK)
%
% 28.09.2013  Add support for 'UseDisplayRotation' via panelfitter. (MK)
%
% 06.03.2014  Add support for 'DualWindowStereo' and fixes to Native10BitFramebuffer mode. (MK)
%
% 26.06.2014  Add support for Native11BitFramebuffer mode, update our docs with what
%                     we learned about this 10/11 bpc business on HDMI so far. (MK)
% 16.09.2014  Add experimental 'StereoCrosstalkReduction' support. (MK/DCN)
% 17.09.2014  Add 'Native16BitFramebuffer' support for Linux + FOSS + AMD. (MK)
% 03.11.2014  Make panelfitter compatible with Retina displays. (MK)
% 04.11.2014  Add new task 'UseRetinaResolution' for Retina displays. (MK)
% 06.09.2015  Add basic support for "Client distortion rendering" on the Oculus VR
%             Rift DK1/DK2 virtual reality headsets. (MK)

persistent configphase_active;
persistent reqs;

% This global variable signals if a GPGPU compute api is enabled, and which
% one. 0 = None, 1 = GPUmat.
global psych_gpgpuapi;

% These flags are global - needed in subfunctions as well (ugly ugly coding):
global ptb_outputformatter_icmAware;
global ptb_MirrorOverlayWindows;
global isASideBySideConfig;
global screenRestoreCmd;
global maxreqarg;

if isempty(configphase_active)
    configphase_active = 0;
    ptb_outputformatter_icmAware = 0;
    maxreqarg = 10;
    ptb_MirrorOverlayWindows = [];
end

if nargin < 1 || isempty(cmd)
    error('You did not provide any subcommand to execute!');
end

rc = [];
winRect = [];

if strcmpi(cmd, 'PrepareConfiguration')
    % Prepare new configuration:
    if configphase_active
        % Huh? Configuration was already in progress. Warn user about reset of task specs:
        fprintf('Tried to prepare a new configuration phase via PsychImaging(''PrepareConfiguration''), but did not finalize the previous phase yet.\n');
        fprintf('You must call the PsychImaging(''OpenWindow'', ...); command at least once to open an onscreen\n');
        fprintf('window according to the provided settings, before you can specify settings for additional onscreen windows.\n');
        fprintf('\n');
        fprintf('The most likely reason you see this error message is because your script aborted with some error\n');
        fprintf('before it managed to open the onscreen window. In that case it is best practice to execute a ''clear all''\n');
        fprintf('command at the Matlab/Octave prompt before you restart your script.\n');
        fprintf('\n');
        fprintf('I will restart configuration now and forget the previously made PsychImaging(''AddTask'', ...); settings.\n');
        warning('Tried to prepare a new configuration phase, but you did not finalize the previous phase yet!');
    end

    % Enter configuration mode, accept 'AddTask' specifications:
    configphase_active = 1;

    % Reset old settings:

    % MK: This clear reqs causes malfunctions on Octave 3.2.0 for some reason, so don't use it! clear reqs;
    reqs = [];
    ptb_outputformatter_icmAware = 0;

    % Set GPGPU api type indicator to zero "none in use" default:
    if isempty(psych_gpgpuapi)
        psych_gpgpuapi = 0;
    end

    % Assign default success return code rc:
    rc = 0;

    return;
end

if strcmpi(cmd, 'AddTask')
    if nargin < 3 || isempty(varargin{1}) || isempty(varargin{2})
        error('Parameters missing: Need at least "whichChannel" and "whichTask"!');
    end

    if configphase_active ~= 1
        error('Call PsychImaging(''PrepareConfiguration''); first to prepare the configuration phase!');
    end

    reqs = AddTask(reqs, varargin{1:end});
    rc = 0;

    return;
end

if strcmpi(cmd, 'FinalizeConfiguration')
    if configphase_active ~= 1
        error('You tried to finalize configuration, but no configuration in progress!');
    end

    if isempty(reqs)
        error('You tried to FinalizeConfiguration, but you did not specify any requirements or tasks!');
    end

    configphase_active = 2;

    % Compute correct imagingMode - Settings for current configuration and
    % return it:
    [imagingMode, needStereoMode, reqs] = FinalizeConfiguration(reqs);
    rc = imagingMode;
    winRect = needStereoMode;

    return;
end

if strcmpi(cmd, 'PostConfiguration')
    if configphase_active ~= 2
        error('Tried to call PostConfiguration without calling FinalizeConfiguration before!');
    end

    if nargin < 2 || isempty(varargin{1}) || Screen('WindowKind', varargin{1})~=1
        error('No "windowPtr" or invalid "windowPtr" or non-onscreen window handle provided!');
    end

    if nargin < 3
        clearcolor = 0;
    else
        clearcolor = varargin{2};
    end

    rc = PostConfiguration(reqs, varargin{1}, clearcolor);

    configphase_active = 0;
    return;
end

if strcmpi(cmd, 'OpenWindow')

    % Allow 'OpenWindow' without task specs. Simply open with empty task requirements list:
    if ismember(configphase_active, [0, 2])
        PsychImaging('PrepareConfiguration');
    end

    if configphase_active ~= 1
        error('You tried to OpenWindow, but didn''t specify any imaging configuration!');
    end

    if nargin < 2
        error('You must supply at least a "screenId" for the screen on which the window should be opened');
    end

    % Final config phase:
    configphase_active = 2; %#ok<NASGU>

    screenid = varargin{1};
    if ~ismember(screenid, Screen('Screens'))
        error('Invalid screenId provided in ''OpenWindow'' - no such screen %i available.', screenid);
    end

    if nargin < 3 || isempty(varargin{2})
        clearcolor = [];
    else
        clearcolor = varargin{2};
    end

    if nargin < 4 || isempty(varargin{3})
        winRect = [];
    else
        winRect = varargin{3};
    end

    % Set override special flags to "none" by default:
    ovrSpecialFlags = [];

    % Set override framebuffer rect to "none" by default:
    ovrfbOverrideRect = [];

    % Override numbuffers -- always 2:
    numbuffers = 2;

    if nargin < 7 || isempty(varargin{6})
        stereomode = 0;
    else
        stereomode = varargin{6};
    end

    if nargin < 8 || isempty(varargin{7})
        multiSample = []; % User defers choice of MSAA to us.
    else
        multiSample = varargin{7};
        if ~isscalar(multiSample) || ~isreal(multiSample)
            error('PsychImaging(''OpenWindow''): Invalid multisample value specified. Not an integer scalar.');
        end
    end

    % Running on a VR headset?
    if ~isempty(find(mystrcmp(reqs, 'UseVRHMD')))
        % Yes:
        floc = find(mystrcmp(reqs, 'UseVRHMD'));
        [rows cols] = ind2sub(size(reqs), floc(1));
        row = rows(1);

        % Extract first parameter - This should be the handle of the HMD device:
        hmd = reqs{row, 3};

        % Verify it is already open:
        if ~hmd.driver('IsOpen', hmd)
            error('PsychImaging(''OpenWindow''): Invalid HMD handle specified for UseVRHMD task. No such device opened.');
        end

        % Compute special OpenWindow overrides for winRect, framebuffer rect, specialflags and MSAA, as needed:
        [winRect, ovrfbOverrideRect, ovrSpecialFlags, multiSample, screenid] = hmd.driver('OpenWindowSetup', hmd, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags, multiSample);
    end

    % If multiSample is still "use default" choice, then override it to our default of 0 for "no MSAA":
    if isempty(multiSample)
        multiSample = 0;
    end

    % Compute correct imagingMode - Settings for current configuration and return it:
    [imagingMode, needStereoMode, reqs] = FinalizeConfiguration(reqs, stereomode, screenid);

    if ~isempty(find(mystrcmp(reqs, 'EnableNative10BitFramebuffer')))
        % Request a pixelsize of 30 bpp to enable native 2101010
        % framebuffer support:
        pixelSize = 30;
    elseif ~isempty(find(mystrcmp(reqs, 'EnableNative11BitFramebuffer')))
        % Request a pixelsize of 33 bpp to enable native RGB11-11-10
        % framebuffer support. A value of 32 bpp would be appropriate but
        % that's already taken by old cruft code, so it's a no-no and we use the
        % weirdo 33 bpp value, to retain backwards compatibility.
        pixelSize = 33;
    elseif ~isempty(find(mystrcmp(reqs, 'EnableNative16BitFramebuffer')))
        % Request a pixelsize of 48 bpp to enable native up to RGB16-16-16
        % framebuffer support.
        pixelSize = 48;
    elseif ~isempty(find(mystrcmp(reqs, 'EnableNative16BitFloatingPointFramebuffer')))
        % Request a pixelsize of 64 bpp to enable native RGBA16F floating point
        % framebuffer support.
        pixelSize = 64;
    else
        % Ignore pixelSize:
        pixelSize = [];
    end

    floc = find(mystrcmp(reqs, 'UseVulkanDisplay'));
    if ~isempty(floc) && isempty(find(mystrcmp(reqs, 'UseStaticHDRHack')))
        [rows cols] = ind2sub(size(reqs), floc(1));
        row = rows(1);

        % Extract first parameter - This would be the optional video output name:
        outputName = reqs{row, 3};

        % Compute special OpenWindow overrides for winRect, framebuffer rect, and specialflags, as needed:
        [winRect, ovrfbOverrideRect, ovrSpecialFlags, outputName] = PsychVulkan('OpenWindowSetup', outputName, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags);

        % Reassign first parameter. Anything non-empty means to use fullscreen/
        % direct display mode:
        reqs{row, 3} = outputName;

        % Reset pixelSize to default 8 bpc, as we are handling potential deep color
        % in the Vulkan backend, not in the OpenGL WSI for the Screen onscreen window:
        pixelSize = [];
    end

    % Override stereomode derived from requirements?
    if needStereoMode ~= -1
        if needStereoMode == -2 && stereomode == 0
            % Stereo operation needed, but not set up by usercode:
            error('Your requirements demand a stereo presentation mode, but you didn''t specify one!');
        else
            if (needStereoMode > -1) && (stereomode ~= needStereoMode)
                % Need a specific mode: Override current setting by our needs:
                stereomode = needStereoMode;

                % Give feedback about stereomode override. If the user
                % didn't provide a stereomode, we just output an info.
                % Otherwise we output a warning about the conflict and our
                % override...
                if nargin < 7 || isempty(varargin{6})
                    fprintf('PsychImaging-Info: Stereomode %i required - Enabling it.\n', stereomode);
                else
                    warning('Your provided "stereomode" conflicts with required stereomode for imaging pipeline. Overriden...');
                end
            end
        end
    end

    if nargin < 9 || isempty(varargin{8})
        imagingovm = 0;
    else
        imagingovm = varargin{8};
    end

    imagingMode = mor(imagingMode, imagingovm);

    if nargin < 10 || isempty(varargin{9})
        specialFlags = ovrSpecialFlags;
    else
        specialFlags = varargin{9};
    end

    % Display mirroring to separate slave window requested?
    if ~isempty(find(mystrcmp(reqs, 'MirrorDisplayTo2ndOutputHead')))
        % Yes. Try to disable vsync for OpenGL bufferswaps on that "mirror" window, as
        % we don't want it to throttle the main stimulus presentation to the potentially
        % lower video refresh rate or drifting video refresh cycle of the "mirror" display,
        % as high performance and proper timing on the subject stimulus display is way more
        % important than a bit of tearing on what is basically just an experimenter monitor:
        if isempty(specialFlags)
            specialFlags = 0;
        end
        specialFlags = mor(specialFlags, kPsychSkipSecondaryVsyncForFlip);
    end

    if nargin < 11 || isempty(varargin{10})
        clientRect = [];
    else
        clientRect = varargin{10};
    end

    if nargin < 12 || isempty(varargin{11})
        fbOverrideRect = [];
    else
        fbOverrideRect = varargin{11};
    end

    % Allow override of the fbOverrideRect from special clients like VR apps:
    if ~isempty(ovrfbOverrideRect)
        fbOverrideRect = ovrfbOverrideRect;
    end

    % Define fitRefRect (used for panel fitter setup) as the effective framebuffer rectangle:
    fitRefRect = winRect;
    if ~isempty(fbOverrideRect)
        fitRefRect = fbOverrideRect;
    end

    % VRR handling:
    ovrvrrParams = [];
    if ~isempty(find(mystrcmp(reqs, 'UseFineGrainedTiming'))) %#ok<*EFIND>
        % Yes. Extract parameters:
        floc = find(mystrcmp(reqs, 'UseFineGrainedTiming'));
        if length(floc) > 1
            error('PsychImaging: Multiple definitions of task "UseFineGrainedTiming"! There can be only one.');
        end

        [row cols] = ind2sub(size(reqs), floc); %#ok<NASGU>
        vrrMode = reqs{row, 3};

        if isempty(vrrMode) || strcmpi(vrrMode, 'Auto')
            vrrMode = 1;
        elseif strcmpi(vrrMode, 'Simple')
            vrrMode = 2;
        elseif strcmpi(vrrMode, 'OwnScheduled')
            vrrMode = 3;
        else
            error('PsychImaging: For task "UseFineGrainedTiming" invalid method argument specified. Must be ''Auto'', ''Simple'', or ''OwnScheduled''.');
        end

        vrrStyleHint = reqs{row, 4};

        if isempty(vrrStyleHint) || strcmpi(vrrStyleHint, 'None')
            vrrStyleHint = 0;
        elseif strcmpi(vrrStyleHint, 'XXX')
            vrrStyleHint = 1;
        else
            error('PsychImaging: For task "UseFineGrainedTiming" invalid styleHint argument specified. Must be ''None''.');
        end

        vrrMinRefreshHz = reqs{row, 5};
        if isempty(vrrMinRefreshHz)
            vrrMinRefreshHz = 0;
        elseif ~isnumeric(vrrMinRefreshHz) || ~isscalar(vrrMinRefreshHz) || vrrMinRefreshHz < 0
            error('PsychImaging: For task "UseFineGrainedTiming" invalid "vrrMinRefreshHz" argument specified.');
        else
            % Map to maximum frame duration:
            if vrrMinRefreshHz > 0
              vrrMinRefreshHz = 1 / vrrMinRefreshHz;
            else
              vrrMinRefreshHz = 0;
            end
        end

        % Build ovrvrrParams:
        ovrvrrParams = [vrrMode, vrrStyleHint, 0, vrrMinRefreshHz];
    end

    if nargin < 13 || isempty(varargin{12})
        vrrParams = ovrvrrParams;
    else
        vrrParams = varargin{12};
    end

    if ~isempty(find(mystrcmp(reqs, 'UseDisplayRotation'))) %#ok<*EFIND>
        % Yes. Extract parameters:
        floc = find(mystrcmp(reqs, 'UseDisplayRotation'));
        if length(floc) > 1
            error('PsychImaging: Multiple definitions of task "UseDisplayRotation"! There can be only one.');
        end

        % Check for collisions with mutually exclusive "UsePanelFitter" task:
        if ~isempty(find(mystrcmp(reqs, 'UsePanelFitter')))
            fprintf('\n\n');
            fprintf('PsychImaging: You can not use both "UseDisplayRotation" and "UsePanelFitter" at the same time. However, you can pass\n');
            fprintf('PsychImaging: the rotation angle you wanted to use for "UseDisplayRotation" to "UsePanelFitter" instead, so "UsePanelFitter"\n');
            fprintf('PsychImaging: will also do the job of "UseDisplayRotation" for you. This works because "UseDisplayRotation" is only\n');
            fprintf('PsychImaging: a simple convenience shortcut to "UsePanelFitter".\n');
            error('PsychImaging: Task "UsePanelFitter" also requested, but you can only use either "UsePanelFitter" or "UseDisplayRotation".');
        end

        [row cols] = ind2sub(size(reqs), floc); %#ok<NASGU>
        rotAngle = reqs{row, 3};

        if isempty(rotAngle) || ~isnumeric(rotAngle) || ~isscalar(rotAngle)
            error('PsychImaging: For task "UseDisplayRotation", required rotation angle parameter missing or not a scalar angle in degrees.');
        end

        % Get full size of output framebuffer:
        if isempty(fitRefRect)
            [clientRes(1), clientRes(2)] = Screen('WindowSize', screenid, 1);
        else
            clientRes = [RectWidth(fitRefRect), RectHeight(fitRefRect)];
        end

        % Rotation into a portrait orientation?
        if (round(rotAngle / 90) == (rotAngle / 90))
            if (mod(round(rotAngle / 90), 2) > 0)
                % Yes. Switch width and height of clientRes:
                clientRes = [clientRes(2), clientRes(1)];
            end
        else
            fprintf('PsychImaging: Provided rotation angle for task "UseDisplayRotation" is not a multiple of 90 degrees.\n');
            fprintf('PsychImaging: You are probably in for a bit of trouble for such rotation angles...\n');
        end

        % No-Op for rotation angle of 0 degrees, as that does nothing.
        if rotAngle ~= 0
            % Build a 'UsePanelFitter' task from our tasks parameters by
            % overwriting our own task spec:
            reqs{row, 2} = 'UsePanelFitter';
            reqs{row, 3} = clientRes;
            reqs{row, 4} = 'Full';
            reqs{row, 5} = [];
            reqs{row, 6} = [];
            reqs{row, 7} = rotAngle;
        end
    end

    % Use and high-level setup of panelfitter requested?
    if ~isempty(find(mystrcmp(reqs, 'UsePanelFitter'))) %#ok<*EFIND>
        % Yes. Extract parameters:
        floc = find(mystrcmp(reqs, 'UsePanelFitter'));
        if length(floc) > 1
            error('PsychImaging: Multiple definitions of task "UsePanelFitter"! There can be only one.');
        end

        [row cols] = ind2sub(size(reqs), floc); %#ok<NASGU>

        % Extract requested resolution of virtual framebuffer...
        clientRes = reqs{row, 3};
        if length(clientRes) ~= 2 || ~isnumeric(clientRes) || min(clientRes) < 1
            error('PsychImaging: Mandatory "size" parameter of task "UsePanelFitter" is missing or not a two component [width, height] size vector with positive width and height as expected.');
        end

        clientRes = round(clientRes);

        if ~isempty(clientRect)
            fprintf('PsychImaging: OpenWindow: Warning: User provided "clientRect" overriden by specification in PsychImaging task "UsePanelFitter".');
        end

        % ... and define clientRect accordingly:
        clientRect = [0, 0, clientRes(1), clientRes(2)];

        % Extract scaling strategy:
        fitterStrategy = reqs{row, 4};
        if isempty(fitterStrategy) || ~ischar(fitterStrategy)
            error('PsychImaging: Mandatory parameter "strategy" of task "UsePanelFitter" missing or not a string.');
        end

        % Define full size of output framebuffer:
        if isempty(fitRefRect)
            dstFit = Screen('Rect', screenid, 1);
        else
            dstFit = SetRect(0, 0, RectWidth(fitRefRect), RectHeight(fitRefRect));
        end

        % Adapt dstFit according to window size flags:

        % Apply half-height flag, if any:
        if bitand(imagingMode, kPsychNeedHalfHeightWindow)
            dstFit(RectBottom) = dstFit(RectBottom) / 2;
        end

        % Apply half-width flag, if any:
        if bitand(imagingMode, kPsychNeedHalfWidthWindow) || ismember(stereomode, [4, 5])
            dstFit(RectRight) = dstFit(RectRight) / 2;
        end

        % Apply twice-width flag, if any:
        if bitand(imagingMode, kPsychNeedTwiceWidthWindow)
            dstFit(RectRight) = dstFit(RectRight) * 2;
        end

        % Apply triple-width flag, if any:
        if bitand(imagingMode, kPsychNeedTripleWidthWindow)
            dstFit(RectRight) = dstFit(RectRight) * 3;
        end

        winCenter = [RectWidth(dstFit)/2, RectHeight(dstFit)/2];

        % Extract rotation angle to use for display rotation:
        rotX = [];
        rotY = [];
        rotAngle = reqs{row, 7};
        rot90Deg = 0;

        if isempty(rotAngle)
            % No rotation angle == zero rotation == no rotation.
            rotAngle = 0;
        else
            % Round to full degrees:
            rotAngle = round(rotAngle);

            if rotAngle ~= 0
                fprintf('PsychImaging: PanelFitter will apply a display rotation of %i degrees.\n', rotAngle);

                % Check if rotation angle is -90, +90, -270, +270, ... degrees,
                % ie. the image is effectively tilted by 90 degrees clockwise
                % or counter-clockwise:
                if ismember(stereomode, [0]) && (round(rotAngle / 90) == (rotAngle / 90)) && (mod(round(rotAngle / 90), 2) > 0)
                    % Yes. This is classic panel rotation. Exchange width and
                    % height of clientRect, so it is "rotated" accordingly and
                    % the various scaling and centering strategies will
                    % peacefully cooperate with display rotation via panel
                    % fitting:
                    rot90Deg = 1;
                    clientRect = [0, 0, clientRes(2), clientRes(1)];
                    fprintf('PsychImaging: Applying special setup for display rotation by 90 degrees into portrait orientation.\n');
                end
            end
        end

        % Which strategy to use?
        if strcmpi(fitterStrategy, 'Custom')
            % Custom scaling with provided srcRect and dstRect:
            srcFit = reqs{row, 5};
            dstFit = reqs{row, 6};
            if ~isnumeric(srcFit) || length(srcFit) ~= 4
                error('PsychImaging: Mandatory parameter "srcRect" of task "UsePanelFitter" for fitting strategy "Custom" missing or not a 4 element rect.');
            end

            if ~isnumeric(dstFit) || length(dstFit) ~= 4
                error('PsychImaging: Mandatory parameter "dstRect" of task "UsePanelFitter" for fitting strategy "Custom" missing or not a 4 element rect.');
            end
        elseif strcmpi(fitterStrategy, 'Centered')
            % Don't rescale but blit one-to-one. Center in target
            % framebuffer, crop if neccessary:

            % Try to center clientRect in destination framebuffer rect:
            srcFit = CenterRect(clientRect, dstFit);

            % Does it fully fit in?
            if any(srcFit < 0)
                % No. We need to crop/clip it to fit in:
                dstFit = ClipRect(srcFit, dstFit);
                srcFit = CenterRect(dstFit, clientRect);
                fprintf('PsychImaging: For centered fitting, i needed to crop the source framebuffer to central region [%i,%i,%i,%i]. Borders will be missing.\n', srcFit(1), srcFit(2), srcFit(3), srcFit(4));
            else
                % Yes: Center in destination framebuffer:
                dstFit = srcFit;
                srcFit = clientRect;
            end
        elseif strcmpi(fitterStrategy, 'Full')
            % Rescale source framebuffer to full target framebuffer, not
            % taking aspect ratio into account:
            srcFit = clientRect;

            if RectWidth(srcFit) / RectHeight(srcFit) ~= RectWidth(dstFit) / RectHeight(dstFit)
                fprintf('PsychImaging: Using full resolution fitting strategy. Scaling will not preserve aspect ratio of original stimulus!\n');
            else
                fprintf('PsychImaging: Using full resolution fitting strategy. Aspect ratio is preserved.\n');
            end
        elseif strcmpi(fitterStrategy, 'AspectWidth') || strcmpi(fitterStrategy, 'AspectHeight') || strcmpi(fitterStrategy, 'Aspect')
            % Rescale aspect ratio preserving:

            if strcmpi(fitterStrategy, 'AspectWidth')
                % Cover full width of window, maybe crop top and bottom:
                sf = RectWidth(dstFit) / RectWidth(clientRect);
                fprintf('PsychImaging: Using scaling to full width. Aspect ratio is preserved, top and bottom may be cut away.\n');
            end

            if strcmpi(fitterStrategy, 'AspectHeight')
                % Cover full width of window, maybe crop top and bottom:
                sf = RectHeight(dstFit) / RectHeight(clientRect);
                fprintf('PsychImaging: Using scaling to full height. Aspect ratio is preserved, left and right margins may be cut away.\n');
            end

            if strcmpi(fitterStrategy, 'Aspect')
                % Cover as much as possible, aspect ratio preserving, leaving
                % borders as neccessary:
                sfw = RectWidth(dstFit) / RectWidth(clientRect);
                sfh = RectHeight(dstFit) / RectHeight(clientRect);
                sf = min(sfw, sfh);
                fprintf('PsychImaging: Using scaling to the most maximal size which still preserves aspect ratio. There may be borders.\n');
            end

            % Compute scaled size target rectangle:
            scaleFit = ScaleRect(clientRect, sf, sf);

            % Center it in destination framebuffer dstFit:
            scaleFit = CenterRect(scaleFit, dstFit);

            % Clip it against dstFit's size, crop away borders if neccessary:
            % dstFit now contains the destination retangle in the window:
            dstFit = ClipRect(scaleFit, dstFit);

            % Compute originating source rectangle of original size for
            % 'dstFit' by undoing the scaling:
            scaleFit = SetRect(0, 0, RectWidth(dstFit)/sf, RectHeight(dstFit)/sf);

            % Center properly sized source rectangle in clientRect source
            % framebuffer to compute final srcRect for scaling blit:
            srcFit = CenterRect(scaleFit, clientRect);
        else
            error('PsychImaging: Mandatory parameter "strategy" of task "UsePanelFitter" has invalid setting ''%s''.', fitterStrategy);
        end

        if rotAngle ~= 0
            [rotX, rotY] = RectCenter(clientRect);
        end

        if rot90Deg
            % Offset compensation for multiple of 90 degrees rotations:
            degrad = 2 * pi * rotAngle / 360;
            rotOffset(1) = -(winCenter(2) - rotX) * sin(degrad);
            rotOffset(2) =  (winCenter(1) - rotY) * sin(degrad);
            dstFit = OffsetRect(dstFit, rotOffset(1), rotOffset(2));
        end

        % Build final fitterParams vector:
        fitterParams = [srcFit dstFit rotAngle rotX rotY];

        % Restore clientRect to original one:
        clientRect = [0, 0, clientRes(1), clientRes(2)];
    else
        % No panel fitter in use. Or at least, none we would set up:
        fitterParams = [];
    end

    % Custom color correction for display wanted on a Bits+ display in
    % Mono++ or Color++ mode or a DataPixx?
    if ~isempty(find(mystrcmp(reqs, 'DisplayColorCorrection')))
        if ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++Output'))) || ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++OutputWithOverlay'))) || ~isempty(find(mystrcmp(reqs, 'EnableBits++Color++Output')))
            % Yes. The BitsPlusPlus() setup routine implements its own
            % setup code for display color correction which is very
            % efficient for a single color correction plugin, but not
            % useable with multiple plugins! Need to handle both
            % cases specially.

            % More than one color correction plugin requested for pipeline?
            floc = find(mystrcmp(reqs, 'DisplayColorCorrection'));
            if length(floc) == 1
                % Single plugin. Use BitsPlusPlus internal setup code,
                % just provide proper method setting for it now:

                % Which channel?
                x=floc;
                [rows cols] = ind2sub(size(reqs), x); %#ok<NASGU>
                for row=rows'
                    % Extract first parameter - This should be the method of correction:
                    colorcorrectionmethod = reqs{row, 3};

                    if isempty(colorcorrectionmethod) || ~ischar(colorcorrectionmethod)
                        sca;
                        error('PsychImaging: Name of color correction method for ''DisplayColorCorrection'' missing or not of string type!');
                    end

                    % Select method:
                    PsychColorCorrection('ChooseColorCorrection', colorcorrectionmethod);
                end
            else
                % Multiple plugins: Select special method which won't be
                % harmful, a simple clamping to valid range, labeled with a
                % special name that can't clash with our own definition of
                % ICM shaders:
                PsychColorCorrection('ChooseColorCorrection', 'ClampedNoName');
            end
        end
    end

    % Open onscreen window with proper imagingMode and stereomode set up.
    % We have a couple of special cases here for VPixx devices and the CRS Bits++/Bits#...
    win = [];

    if ~isempty(find(mystrcmp(reqs, 'EnableBits++Bits++Output')))
        % Special case: Need to open Bits++ Bits++ driver. We delegate the
        % openwindow procedure to the BitsPlusPlus.m file:
        if ~isempty(win)
            error('You specified multiple conflicting output display device drivers! This will not work.');
        end

        if nargin >= 13
            [win, winRect] = BitsPlusPlus('OpenWindowBits++', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingMode, specialFlags, clientRect, fbOverrideRect, vrrParams, varargin{13:end});
        else
            [win, winRect] = BitsPlusPlus('OpenWindowBits++', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingMode, specialFlags, clientRect, fbOverrideRect, vrrParams);
        end
    end

    if ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++Output'))) || ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++OutputWithOverlay')))
        % Special case: Need to open Bits++ Mono++ driver. We delegate the
        % openwindow procedure to the BitsPlusPlus.m file:
        if ~isempty(win)
            error('You specified multiple conflicting output display device drivers! This will not work.');
        end

        if ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++OutputWithOverlay')))
            bpcom = 'OpenWindowMono++WithOverlay';
        else
            bpcom = 'OpenWindowMono++';
        end

        if nargin >= 13
            [win, winRect] = BitsPlusPlus(bpcom, screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingMode, specialFlags, clientRect, fbOverrideRect, vrrParams, varargin{13:end});
        else
            [win, winRect] = BitsPlusPlus(bpcom, screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingMode, specialFlags, clientRect, fbOverrideRect, vrrParams);
        end
    end

    if ~isempty(find(mystrcmp(reqs, 'EnableBits++Color++Output')))
        % Special case: Need to open Bits++ Color++ driver. We delegate the
        % openwindow procedure to the BitsPlusPlus.m file:
        if ~isempty(win)
            error('You specified multiple conflicting output display device drivers! This will not work.');
        end

        if nargin >= 13
            [win, winRect] = BitsPlusPlus('OpenWindowColor++', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingMode, specialFlags, clientRect, fbOverrideRect, vrrParams, varargin{13:end});
        else
            [win, winRect] = BitsPlusPlus('OpenWindowColor++', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingMode, specialFlags, clientRect, fbOverrideRect, vrrParams);
        end
    end

    if isempty(win)
        % Standard openwindow path:
        if nargin >= 13
            [win, winRect] = Screen('OpenWindow', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingMode, specialFlags, clientRect, fbOverrideRect, vrrParams, varargin{13:end});
        else
            [win, winRect] = Screen('OpenWindow', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingMode, specialFlags, clientRect, fbOverrideRect, vrrParams);
        end
    end

    % No secondary slave window by default:
    slavewin = [];

    % Display mirroring requested?
    if ~isempty(find(mystrcmp(reqs, 'MirrorDisplayTo2ndOutputHead')))
        % Yes. Need to open secondary slave window:
        floc = find(mystrcmp(reqs, 'MirrorDisplayTo2ndOutputHead'));
        [rows cols]= ind2sub(size(reqs), floc);

        % Extract first parameter - This should be the id of the slave
        % screen to which the display should get mirrored:
        slavescreenid = reqs{rows, 3};

        if isempty(slavescreenid)
            sca;
            error('In PsychImaging MirrorDisplayTo2ndOutputHead: You must provide the index of the secondary screen "slavescreen"!');
        end

        if ~any(ismember(Screen('Screens'), slavescreenid))
            sca;
            error('In PsychImaging MirrorDisplayTo2ndOutputHead: You must provide the index of a valid secondary screen "slavescreen"!');
        end

        if ~ismember(stereomode, [0, 2:9])
            sca;
            error('In PsychImaging MirrorDisplayTo2ndOutputHead: Tried to simultaneously enable stereomode %i! This is not supported.', stereomode);
        end

        % Extract optional 2nd parameter - The window rectangle of the slave
        % window on the slave screen to which the display should get mirrored:
        slavewinrect = reqs{rows, 4};
        if isempty(slavewinrect)
            slavewinrect = []; % From empty string to empty vector.
        end

        % Get optional slavewinflags, e.g., kPsychGUIWindow:
        slavewinflags = reqs{rows, 5};
        if isempty(slavewinflags)
            slavewinflags = []; % From empty string to empty vector.
        end

        % Open slave window on slave screen: Set the special dual window
        % output flag, so Screen('OpenWindow') initializes the internal blit
        % chain properly. Suppress sync tests and verbose output during open:
        oldverbose = Screen('Preference', 'Verbosity', 1);
        oldskip = Screen('Preference', 'SkipSyncTests', 2);
        slavewin = Screen('OpenWindow', slavescreenid, 30, slavewinrect, pixelSize, [], [], [], kPsychNeedDualWindowOutput, slavewinflags);
        Screen('Preference', 'SkipSyncTests', oldskip);
        Screen('Preference', 'Verbosity', oldverbose);
    end

    % Dualwindow output requested? [Essentially the same as display
    % mirroring, but kept separate for now for simplicity]
    if ~isempty(find(mystrcmp(reqs, 'EnableDualPipeHDROutput')))
        % Yes. Need to open secondary slave window:
        floc = find(mystrcmp(reqs, 'EnableDualPipeHDROutput'));
        [rows cols]= ind2sub(size(reqs), floc);

        % Extract first parameter - This should be the id of the slave
        % screen to which the pipe 1 display should get displayed:
        slavescreenid = reqs{rows, 3};

        if isempty(slavescreenid)
            sca;
            error('In PsychImaging EnableDualPipeHDROutput: You must provide the index of the secondary screen "slavescreen"!');
        end

        if ~any(ismember(Screen('Screens'), slavescreenid))
            sca;
            error('In PsychImaging EnableDualPipeHDROutput: You must provide the index of a valid secondary screen "slavescreen"!');
        end

        if stereomode == 1
            sca;
            error('In PsychImaging EnableDualPipeHDROutput: Tried to simultaneously enable frame-sequential stereomode 1! This is not supported.');
        end

        if stereomode == 10
            sca;
            error('In PsychImaging EnableDualPipeHDROutput: Tried to simultaneously enable dual display output stereomode 10! This is not supported.');
        end

        % Extract optional 2nd parameter - The window rectangle of the slave
        % window on the slave screen to which the pipe 1 display should get outputted:
        slavewinrect = reqs{rows, 4};
        if isempty(slavewinrect), slavewinrect = []; end

        % Open slave window on slave screen: Set the special dual window
        % output flag, so Screen('OpenWindow') initializes the internal blit
        % chain properly:
        slavewin = Screen('OpenWindow', slavescreenid, [255 0 0], slavewinrect, pixelSize, [], [], [], kPsychNeedDualWindowOutput);
    end

    % DualWindow stereo output requested?
    if ~isempty(find(mystrcmp(reqs, 'DualWindowStereo')))
        % Yes. Need to open secondary slave window:
        floc = find(mystrcmp(reqs, 'DualWindowStereo'));
        [rows cols]= ind2sub(size(reqs), floc);

        % Extract first parameter - This should be the id of the slave
        % screen to which the right eye display should get displayed:
        slavescreenid = reqs{rows, 3};

        if isempty(slavescreenid)
            sca;
            error('In PsychImaging DualWindowStereo: You must provide the index of the secondary screen "slavescreen"!');
        end

        if ~any(ismember(Screen('Screens'), slavescreenid))
            sca;
            error('In PsychImaging DualWindowStereo: You must provide the index of a valid secondary screen "slavescreen"!');
        end

        % Extract optional 2nd parameter - The window rectangle of the slave
        % window on the slave screen:
        slavewinrect = reqs{rows, 4};
        if isempty(slavewinrect), slavewinrect = []; end

        % Open slave window on slave screen:
        slavewin = Screen('OpenWindow', slavescreenid, [], slavewinrect, pixelSize, [], 10);
    end

    % Matlab? Does the Java swing cleanup function exist?
    if exist('PsychJavaSwingCleanup', 'file')
        % Attach a window close callback for cleanup of Java's memory
        % management mess at window close time when Matlab with Java based
        % GUI is in use:
        Screen('Hookfunction', win, 'AppendMFunction', 'CloseOnscreenWindowPostGLShutdown', 'Shutdown window callback into PsychJavaSwingCleanup().', 'PsychJavaSwingCleanup;');
        Screen('HookFunction', win, 'Enable', 'CloseOnscreenWindowPostGLShutdown');

        % Some slave window opened?
        if ~isempty(slavewin)
            % Yes: Apply java cleanup there as well:
            Screen('Hookfunction', slavewin, 'AppendMFunction', 'CloseOnscreenWindowPostGLShutdown', 'Shutdown window callback into PsychJavaSwingCleanup().', 'PsychJavaSwingCleanup;');
            Screen('HookFunction', slavewin, 'Enable', 'CloseOnscreenWindowPostGLShutdown');
        end
    end

    glerr = glGetError;
    while glerr
        fprintf('PsychImaging-WARNING:OpenWindow: OpenGL error detected after Screen(''OpenWindow''): %s', gluErrorString(glerr));
        glerr = glGetError;
    end

    % Window open. Perform imaging pipe postconfiguration:
    PostConfiguration(reqs, win, clearcolor, slavewin);

    glerr = glGetError;
    while glerr
        fprintf('PsychImaging-WARNING:OpenWindow: OpenGL error detected after PostConfiguration(): %s', gluErrorString(glerr));
        glerr = glGetError;
    end

    % Perform double-flip, so both back- and frontbuffer get initialized to
    % background color:
    Screen('Flip', win);
    Screen('Flip', win);

    % Panel fitter in use and setup by us?
    if ~isempty(fitterParams)
        % Yes: Apply fitter parameters now, so the scaling method takes
        % effect at next flip. We only do it now, so the preceeding
        % Screen('Flip') ops after imaging pipeline initialization were
        % able to operate with the default "cover full framebuffer" fitter
        % params, ie., they applied their implicit "clear to background
        % color" ops to the full framebuffer and thereby initialized all
        % stages of the pipeline down to the real window backbuffer with
        % background clear color. This way, regardless which panel fitting
        % strategy is chosen by user code, potential top-bottom or
        % left-right borders will get initialized to the selected
        % background clear color, which should be the most well defined
        % choice:
        Screen('PanelFitter', win, round(fitterParams));

        % Now that the fitter is fully configured, perform an extra
        % double-flip to apply proper scaling and borders and such:
        Screen('Flip', win);
        Screen('Flip', win);
    end

    % One extra Flip to put the full imaging pipeline into initial state:
    Screen('Flip', win);

    glerr = glGetError;
    while glerr
        fprintf('PsychImaging-WARNING:OpenWindow: OpenGL error detected after init Flip: %s', gluErrorString(glerr));
        glerr = glGetError;
    end

    rc = win;

    % Done.
    configphase_active = 0;

    return;
end

if strcmpi(cmd, 'RestrictProcessingToROI')
    % Define a ROI in a processing chain/channel to which processing should
    % be restricted by internal use of glScissor() command. This is a
    % runtime function. Each invocation will search the given channel if
    % such a command already exists, then delete it if so. It will prepend
    % the new command with the new spec in any case, so that at any point
    % in time exactly one such ROI can be active for a chain:

    if nargin < 4
        sca;
        error('You must provide all parameters for subfunction "RestrictProcessingToROI!"');
    end

    % Extract window handle:
    win = varargin{1};

    if ~isscalar(win) || ~isnumeric(win) || Screen('WindowKind', win) ~= 1
        sca;
        error('Provided window parameter for subfunction "RestrictProcessingToROI!" is not the handle of a valid onscreen window!');
    end

    % Extract window information:
    winfo = Screen('GetWindowInfo', win);

    % Extract view channel:
    whichView = varargin{2};

    % Extract scissor rectangle:
    scissorrect = varargin{3};

    if size(scissorrect,1)~=1 || size(scissorrect,2)~=4
        sca;
        error('Command "RestrictProcessingToROI" in channel %s expects a 1-by-4 ROI rectangle to define the ROI, e.g, [left top right bottom]!', whichView);
    end

    ox = scissorrect(RectLeft);

    [winwidth, winheight] = InterBufferSize(win);
    oy = winheight - scissorrect(RectBottom);

    w  = RectWidth(scissorrect);
    h  = RectHeight(scissorrect);

    if mystrcmp(whichView, 'LeftView') || mystrcmp(whichView, 'AllViews')
        % Need to restrict left view processing:
        DoRemoveScissorRestriction(win, 'StereoLeftCompositingBlit');
        Screen('HookFunction', win, 'PrependBuiltin', 'StereoLeftCompositingBlit', 'Builtin:RestrictToScissorROI', sprintf('%i:%i:%i:%i', ox, oy, w, h));
    end

    if mystrcmp(whichView, 'RightView') || mystrcmp(whichView, 'AllViews')
        % Need to restrict right view processing:
        DoRemoveScissorRestriction(win, 'StereoRightCompositingBlit');
        Screen('HookFunction', win, 'PrependBuiltin', 'StereoRightCompositingBlit', 'Builtin:RestrictToScissorROI', sprintf('%i:%i:%i:%i', ox, oy, w, h));
    end

    if (mystrcmp(whichView, 'AllViews') || mystrcmp(whichView, 'Compositor')) && ismember(winfo.StereoMode, [6,7,8,9])
        % Needed to restrict both views processing and a
        % compositing mode is active. If both views are restricted
        % in their output area then it makes sense to restrict the
        % compositor to the same area. We also restrict the
        % compositor if that was requested.
        oy = RectHeight(Screen('Rect', win, 1)) - scissorrect(RectBottom);
        DoRemoveScissorRestriction(win, 'StereoCompositingBlit');
        Screen('HookFunction', win, 'PrependBuiltin', 'StereoCompositingBlit', 'Builtin:RestrictToScissorROI', sprintf('%i:%i:%i:%i', ox, oy, w, h));
    end

    if mystrcmp(whichView, 'FinalFormatting')
        % Need to restrict final formatting blit processing:
        oy = RectHeight(Screen('Rect', win, 1)) - scissorrect(RectBottom);
        DoRemoveScissorRestriction(win, 'FinalOutputFormattingBlit');
        Screen('HookFunction', win, 'PrependBuiltin', 'FinalOutputFormattingBlit', 'Builtin:RestrictToScissorROI', sprintf('%i:%i:%i:%i', ox, oy, w, h));
    end

    return;
end

if strcmpi(cmd, 'UnrestrictProcessing')
    % Remove a ROI in a processing chain/channel to which processing should
    % be restricted by internal use of glScissor() command. This is a
    % runtime function. Each invocation will search the given channel if
    % such a command exists, then delete it if so.
    if nargin < 3
        sca;
        error('You must provide all parameters for subfunction "UnrestrictProcessing!"');
    end

    % Extract window handle:
    win = varargin{1};

    if ~isscalar(win) || ~isnumeric(win) || Screen('WindowKind', win) ~= 1
        sca;
        error('Provided window parameter for subfunction "UnrestrictProcessing!" is not the handle of a valid onscreen window!');
    end

    % Extract window information:
    winfo = Screen('GetWindowInfo', win);

    % Extract view channel:
    whichView = varargin{2};

    if mystrcmp(whichView, 'LeftView') || mystrcmp(whichView, 'AllViews')
        % Need to restrict left view processing:
        DoRemoveScissorRestriction(win, 'StereoLeftCompositingBlit');
    end

    if mystrcmp(whichView, 'RightView') || mystrcmp(whichView, 'AllViews')
        % Need to restrict right view processing:
        DoRemoveScissorRestriction(win, 'StereoRightCompositingBlit');
    end

    if (mystrcmp(whichView, 'AllViews') || mystrcmp(whichView, 'Compositor')) && ismember(winfo.StereoMode, [6,7,8,9])
        % Needed to restrict both views processing and a
        % compositing mode is active. If both views are restricted
        % in their output area then it makes sense to restrict the
        % compositor to the same area. We also restrict the
        % compositor if that was requested.
        DoRemoveScissorRestriction(win, 'StereoCompositingBlit');
    end

    if mystrcmp(whichView, 'FinalFormatting')
        % Need to restrict final formatting blit processing:
        DoRemoveScissorRestriction(win, 'FinalOutputFormattingBlit');
    end

    return;
end

if strcmpi(cmd, 'GetOverlayWindow')
    % Pass this call through to BitsPlusPlus.m driver -- the only one which
    % currently supports such overlays.
    % MK: May need to do something more clever in the future...

    % rc is the 'win'dowhandle, winRect is its Screen('Rect'):
    [rc, winRect] = BitsPlusPlus('GetOverlayWindow', varargin{:});

    return;
end

if strcmpi(cmd, 'GetMirrorOverlayWindow')
    % Assign onscreen window index:
    if length(varargin) < 1 || isempty(varargin{1}) || ~isa(varargin{1}, 'double') || Screen('WindowKind', varargin{1}) ~= 1
        error('PsychImaging "GetMirrorOverlayWindow" called without valid onscreen window handle.');
    end
    win = varargin{1};

    if win < 1 || win > length(ptb_MirrorOverlayWindows)
        error('PsychImaging "GetMirrorOverlayWindow": No overlay associated with given onscreen window.');
    end

    if ptb_MirrorOverlayWindows(win) == 0
        error('PsychImaging "GetMirrorOverlayWindow": No overlay associated with given onscreen window.');
    end

    % Ok, this 'win'dow has an overlay: Return its offscreen 'win'dow handle:
    rc = ptb_MirrorOverlayWindows(win);

    return;
end

% Catch all for unknown commands:
error('Unknown subcommand specified! Read "help PsychImaging" for usage info.');
return; %#ok<UNRCH>

% Internal helper routines:

% FinalizeConfiguration consolidates the current set of requirements and
% derives the needed stereoMode settings and imagingMode setting to pass to
% Screen('OpenWindow') for pipeline preconfiguration.
function [imagingMode, stereoMode, reqs] = FinalizeConfiguration(reqs, userstereomode, screenid)
global ptb_outputformatter_icmAware;
global psych_gpgpuapi;
global isASideBySideConfig;
global screenRestoreCmd;
global maxreqarg;

% Reset flag to "no":
isASideBySideConfig = 0;

if nargin < 2
    userstereomode = [];
end

if isempty(userstereomode)
    userstereomode = 0;
end

if nargin < 3 || isempty(screenid)
    screenid = max(Screen('Screens'));
end

% Set imagingMode to minimum: Pipeline disabled. All latter task
% requirements will setup imagingMode to fullfill their needs. A few
% tasks/requirements don't need the full pipeline at all. E.g, Support for
% fast offscreen windows only needs that, but not the full pipeline. Some
% of the "software based mirror modes" herein only need the finalizer blit
% chains, but not the imaging pipeline. Bits++ setup for pure CLUT imaging
% (Bits++ mode) doesn't need imaging pipe either...
imagingMode = 0;

% Set stereoMode to don't care:
stereoMode = -1;

% No datapixx by default:
datapixxmode = 0;

% No Bits+ or Bits# by default:
crsbitsdevice = 0; %#ok<NASGU>

% Request for GPGPU compute support?
floc = find(mystrcmp(reqs, 'UseGPGPUCompute'));
if ~isempty(floc)
    % Yes.
    [row cols] = ind2sub(size(reqs), floc); %#ok<NASGU>

    % Extract first mandatory parameter, the apitype to use:
    apitype = reqs{row, 3};
    if ~ischar(apitype) || (~strcmpi(apitype, 'Auto') && ~strcmpi(apitype, 'GPUmat'))
        % Missing or invalid apitype specified:
        sca;
        error('PsychImaging: Use of GPU compute device via UseGPGPUCompute was requested, but mandatory apitype parameter is missing or invalid!');
    end

    % Extract 2nd optional parameter, the compute flags:
    gpgpuflags = reqs{row, 4};
    if ~isempty(gpgpuflags) && ~ischar(gpgpuflags)
        % There ain't no valid flags yet, so providing anything but the empty string is invalid:
        sca;
        error('PsychImaging: Use of GPU compute device via UseGPGPUCompute was requested, but optional flags argument is invalid!');
    end

    % Ok, all parameters validated. Check if our only currently supported
    % GPU compute api, 'GPUmat' is installed and functional and start it,
    % if possible:
    if ~exist('GPUstart', 'file')
        % Unsupported:
        sca;
        error('PsychImaging: Use of GPU compute device via UseGPGPUCompute was requested, but the required GPUmat toolbox seems to be missing!');
    end

    % Available. Start it:
    psychlasterror('reset');
    try
        % Start/Initialize GPUmat GPU computing toolkit if not already started:
        if ~GPUstart(1)
            GPUstart;
        end
    catch %#ok<CTCH>
        fprintf('PsychImaging: Failed to start GPGPU compute toolkit GPUmat! See error message below:\n');
        err = psychlasterror('reset');
        disp(err.message);
        sca;
        error('PsychImaging: GPGPU init failed!');
    end

    % Ok, GPUmat is online. Set a global marker that it is running:
    fprintf('PsychImaging: GPGPU computing support via GPUmat toolbox enabled.\n');

    % Type 1 is GPUmat:
    psych_gpgpuapi = 1; %#ok<NASGU>
end

% Special setup for CRS Bits# next-generation devices:
% Is a Bits+ / Bits# specific video display mode requested? Or
% explicit use of a Bits# device?
floc = [ find(mystrcmp(reqs, 'EnableBits++Bits++Output')) ];
floc = [floc(:) ; find(mystrcmp(reqs, 'EnableBits++Mono++Output')) ; find(mystrcmp(reqs, 'EnableBits++Mono++OutputWithOverlay')) ];
floc = [floc(:) ; find(mystrcmp(reqs, 'EnableBits++Color++Output')) ; find(mystrcmp(reqs, 'UseBits#')) ];
if ~isempty(floc)
    % Explicit use of Bits# requested? Or only implicit by video mode?
    floc = find(mystrcmp(reqs, 'UseBits#'));
    if ~isempty(floc)
        % Use of Bits# requested. Try to retrieve any special Bits# parameters to
        % pass them to the OpenBits# function:
        [row cols] = ind2sub(size(reqs), floc);

        % Extract first parameter - This should be the serial port name, or [] empty:
        bitsSharpPortname = reqs{row, 3};
    else
        % No specific usage of Bits# requested. Leave it to auto-detection
        % if we work with a Bits# or with a Bits+:
        bitsSharpPortname = [];
    end

    % Initialize serial port connection to Bits#, if any such device present:
    if BitsPlusPlus('OpenBits#', bitsSharpPortname)
        % Connection to Bits# established. Do we need to explicitely
        % specify use of it? Only if it was not already done by usercode via
        % keyword UseBits#
        if isempty(floc)
            % Bits# connected. Makeit explicit by adding the reqs task UseBits#
            reqs(end+1, :) = cell(1, size(reqs, 2));
            reqs{end, 2} = 'UseBits#';
        end

        % Mark use of Bits#:
        crsbitsdevice = 2;

        fprintf('PsychImaging: Will use a connected CRS Bits# device instead of a Bits+ for this session - Connection established.\n');
    else
        % No connection to Bits#. Was one requested? If not, we just assume we are
        % operating against a good old Bits+ which does not support connections.
        % Otherwise, failure to connect to Bits# would be, well, a failure:
        if ~isempty(floc)
            % Bummer:
            sca;
            error('PsychImaging: Use of a CRS Bits# device was requested, but connecting to it failed. Disconnected or misconfigured?!?');
        else
            % Mark use of Bits+:
            crsbitsdevice = 1;

            fprintf('PsychImaging: Will use a CRS Bits+ device, which i assume is connected to target display output screen.\n');
        end
    end
end

% End of Bits# setup, start of DataPixx/ViewPixx/ProPixx setup:

% Remap Datapixx L48 mode to equivalent Bits++ mode:
floc = find(mystrcmp(reqs, 'EnableDataPixxL48Output'));
if ~isempty(floc)
    reqs{floc} = 'EnableBits++Bits++Output';
    reqs(end+1, :) = cell(1, size(reqs, 2));
    reqs{end, 2} = 'UseDataPixx';
    datapixxmode = 1;

    % Initialize connection, switch immediately to L48 mode:
    PsychDataPixx('Open');
    PsychDataPixx('SetVideoMode', 1);
end

% Remap Datapixx M16 mode to equivalent Mono++ mode:
floc = find(mystrcmp(reqs, 'EnableDataPixxM16Output'));
if ~isempty(floc)
    reqs{floc} = 'EnableBits++Mono++Output';
    reqs(end+1, :) = cell(1, size(reqs, 2));
    reqs{end, 2} = 'UseDataPixx';
    datapixxmode = 1;

    % Initialize connection, switch immediately to M16 mode:
    PsychDataPixx('Open');
    PsychDataPixx('SetVideoMode', 2);
end

% Remap Datapixx M16 mode to equivalent Mono++ mode with overlay:
floc = find(mystrcmp(reqs, 'EnableDataPixxM16OutputWithOverlay'));
if ~isempty(floc)
    reqs{floc} = 'EnableBits++Mono++OutputWithOverlay';
    reqs(end+1, :) = cell(1, size(reqs, 2));
    reqs{end, 2} = 'UseDataPixx';
    datapixxmode = 1;

    % Initialize connection, switch immediately to M16 mode:
    PsychDataPixx('Open');
    PsychDataPixx('SetVideoMode', 2);
end

% Remap Datapixx C48 mode to equivalent Color++ mode:
floc = find(mystrcmp(reqs, 'EnableDataPixxC48Output'));
if ~isempty(floc)
    reqs{floc} = 'EnableBits++Color++Output';
    reqs(end+1, :) = cell(1, size(reqs, 2));
    reqs{end, 2} = 'UseDataPixx';
    datapixxmode = 1;

    % Initialize connection, switch immediately to C48 mode:
    PsychDataPixx('Open');
    PsychDataPixx('SetVideoMode', 3);
end

% Assign opmode to BitsPlusPlus driver: It unifies code for Bits+ and
% Datapixx:
BitsPlusPlus('SetTargetDeviceType', datapixxmode);

% Are we setting up for a Datapixx display?
if ~isempty(find(mystrcmp(reqs, 'UseDataPixx')))
    % Yes. Device connection already open from video mode setup above?
    % If not, open connection now.
    if datapixxmode == 0
        % Open connection:
        PsychDataPixx('Open');

        % As no other special high precision output mode is requested, set
        % video mode to "normal passthrough":
        PsychDataPixx('SetVideoMode', 0);

        % Mark as online:
        datapixxmode = 1;
    end
end

% Want native Retina display resolution in a scaled HiDPI display mode?
if ~isempty(find(mystrcmp(reqs, 'UseRetinaResolution')))
    imagingMode = mor(imagingMode, kPsychNeedRetinaResolution);
end

% FBO backed framebuffer needed?
if ~isempty(find(mystrcmp(reqs, 'UseVirtualFramebuffer')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
end

% Use Eizo RadiForce et al. subpixel drive pixel packing?
if ~isempty(find(mystrcmp(reqs, 'UseSubpixelDrive')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore, kPsychNeedTripleWidthWindow, kPsychNeedOutputConversion);
    ptb_outputformatter_icmAware = 0;
end

% 16 bit integer precision framebuffer needed? This is only supported on
% ATI hardware...
if ~isempty(find(mystrcmp(reqs, 'FixedPoint16Bit')))
    imagingMode = mor(imagingMode, kPsychNeed16BPCFixed);
end

% Stereomode 6 for interleaved line stereo needed?
if ~isempty(find(mystrcmp(reqs, 'InterleavedLineStereo')))
    % Yes: Must use stereomode 6.
    stereoMode = 6;
    % We also request an effective window height that is only half the real
    % height. This affects all drawing and query commands of Screen:
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore, kPsychNeedHalfHeightWindow);
end

% Stereomode 6 for interleaved column stereo needed?
if ~isempty(find(mystrcmp(reqs, 'InterleavedColumnStereo')))
    % Yes: Must use stereomode 6.
    stereoMode = 6;
    % We also request an effective window width that is only half the real
    % width. This affects all drawing and query commands of Screen:
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore, kPsychNeedHalfWidthWindow);
end

% Stereomode 2 for side-by-side compressed stereo needed?
if ~isempty(find(mystrcmp(reqs, 'SideBySideCompressedStereo')))
    % Yes: Must use stereomode 2.
    stereoMode = 2;
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
end

% Stereomode 10 for dualwindow stereo needed?
if ~isempty(find(mystrcmp(reqs, 'DualWindowStereo')))
    % Yes: Must use stereomode 10.
    stereoMode = 10;
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
end

% Does usercode or internal code request a stereomode?
if userstereomode > 0 || stereoMode > 0
    % Enable imaging pipeline based stereo,ie., kPsychNeedFastBackingStore:
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);

    % Datapixx - if any - needs special setup:
    if datapixxmode
        % Datapixx device active:

        % Frame sequential style mode via top-down "sync-doubling" mode?
        if ismember(userstereomode, 2) || ismember(stereoMode, 2)
            % Switch Datapixx to sync-doubling stereo mode:
            PsychDataPixx('SetVideoVerticalStereo', 1);
            % Reduce height of virtual framebuffer to effective half height:
            % Nope... imagingMode = mor(imagingMode, kPsychNeedHalfHeightWindow);
        else
            % Switch Datapixx to non stereo mode:
            PsychDataPixx('SetVideoVerticalStereo', 0);
        end

        % Frame-Sequential stereo driven by GPU or us, instead of Datapixx?
        if ismember(userstereomode, [1,11]) || ismember(stereoMode, [1,11])
            % Ask Datapixx to interpret the blue-line-sync stereo sync line
            % generated by either Screen() itself, or the graphics card
            % stereo device driver. The device will drive its VESA Stereo
            % Mini-DIN connector accordingly:
            PsychDataPixx('EnableVideoStereoBlueline');
        else
            % Datapixx shall ignore stereo sync lines / treat the bottom
            % scanline as any other scanline:
            PsychDataPixx('DisableVideoStereoBlueline');
        end

        % Any frame-sequential stereo mode active?
        if ismember(userstereomode, [1,11,2]) || ismember(stereoMode, [1,11,2])
            % Ask ViewPixx to enable its scanning backlight for faster /
            % ghost-free response:
            PsychDataPixx('EnableVideoScanningBacklight');

            % Tell Screen() to tolerate a VBLANK interval that is up to 50%
            % the height of VACTIVE, ie. allow a max VTOTAL = 1.5 * VACTIVE.
            % This is needed because Screen's beamposition query startup
            % tests may otherwise falsely conclude broken beamposition
            % query support due to the unusually large VBLANK interval used
            % by ViewPixx - and possibly other VPixx devices - in frame
            % sequential stereo mode. Our normal rejection threshold is 25%
            % or 1.25, now we raise it to 50% or 1.5:
            Screen('Preference', 'VBLEndlineOverride', [], 1.5);
        else
            % ViewPixx shall disable scanning backlight by default:
            PsychDataPixx('DisableVideoScanningBacklight');
        end

        % Dual-Display stereo via left-right stereo?
        if ismember(userstereomode, [4,5]) || ismember(stereoMode, [4,5])
            % Switch Datapixx to stereo mode by splitting display
            % horizontally onto 2 displays:
            PsychDataPixx('SetVideoHorizontalSplit', 1);
        else
            % Switch Datapixx to non stereo mode, aka auto mode:
            PsychDataPixx('SetVideoHorizontalSplit', 2);
        end
    end
end

% Want to reduce crosstalk in stereo presentation modes?
if ~isempty(find(mystrcmp(reqs, 'StereoCrosstalkReduction')))
    % Yes: For now we only implement this experimentally and for attachment
    % of crosstalk reduction shaders to the image processing chains.
    % This will be suboptimal if other image processing ops are active,
    % but for a first usefully working prototype it should be good enough.
    %
    % We only request additional access to the other image channel, as setup
    % code above and below will already have activated the image processing
    % chains etc.
    imagingMode = mor(imagingMode, kPsychNeedOtherStreamInput);
end

% Want to use a VR Head mounted display (HMD)?
floc = find(mystrcmp(reqs, 'UseVRHMD'));
if ~isempty(floc)
    % Yes: We need a peculiar configuration, which involves the panelfitter
    % to allow for a custom resolution of the virtual framebuffers for left
    % eye and right eye - much higher than output resolution, so we have enough
    % excess information to deal with geometric undistortion warps, color aberration,
    % and dynamic display warping for head motion correction. We also need a
    % special stereo processing shader that does geometric distortion correction,
    % color aberration correction, vignetting correction, and dynamic display warping
    % in one go, as processing speed is crucial for VR experience.
    [rows cols] = ind2sub(size(reqs), floc(1));
    row = rows(1);

    % Extract first parameter - This should be the handle of the HMD device:
    hmd = reqs{row, 3};

    % Verify it is already open:
    if ~hmd.driver('IsOpen', hmd)
        error('UseVRHMD: Invalid HMD handle specified. No such device opened.');
    end

    % Get imagingMode flags and stereoMode to use for this HMD:
    [clientRes, imagingFlags, stereoMode] = hmd.driver('GetClientRenderingParameters', hmd);

    % Add imaging mode flags requested by HMD driver:
    imagingMode = mor(imagingMode, imagingFlags);

    % Do we need the PanelFitter?
    needPanelFitter = hmd.driver('GetPanelFitterParameters', hmd);
    if needPanelFitter
        x{1} = 'General';
        x{2} = 'UsePanelFitter';
        x{3} = clientRes;
        x{4} = 'Custom';
        x{5} = [0, 0, clientRes(1), clientRes(2)];
        x{6} = [0, 0, clientRes(1), clientRes(2)];

        % Pad to maxreqarg arguments:
        if length(x) < maxreqarg
            for i=length(x)+1:maxreqarg
                x{i}='';
            end
        end
        reqs = [reqs ; x];
    end
end

% Request for native 16 bit per color component RGBA16161616 framebuffer on Linux with native X11 X-Server?
if ~isempty(find(mystrcmp(reqs, 'EnableNative16BitFramebuffer'))) && IsLinux && ~IsWayland
    % Add the UseVulkanDisplay task if it doesn't exist already, as our Linux X11
    % 16 bpc support is implemented on top of Vulkan:
    floc = find(mystrcmp(reqs, 'UseVulkanDisplay'));
    if isempty(floc)
        reqs = AddTask(reqs, 'General', 'UseVulkanDisplay');
    end

    % No frame-sequential stereo yet:
    if ismember(userstereomode, [1, 11])
        error('PsychImaging: Requested task ''EnableNative16BitFramebuffer'' is incompatible with frame-sequential stereo mode %i.', userstereomode);
    end

    % Check for Vulkan in general:
    if ~PsychVulkan('Supported')
        % Failed/Unsupported.
        error('PsychImaging: Requested task ''EnableNative16BitFramebuffer'', but this system does not support the required Vulkan api at all.');
    end
end

floc = find(mystrcmp(reqs, 'EnableHDR'));
if ~isempty(floc)
    % Get and validate input arguments:
    [rows cols] = ind2sub(size(reqs), floc(1));
    row = rows(1);
    hdrArguments = reqs(row, :);

    % Check if HDR is supported at all on this system:
    if ~PsychHDR('Supported') && isempty(strfind(hdrArguments{5}, 'Dummy'))
        % Failed/Unsupported.
        error('PsychImaging: Requested task ''EnableHDR'', but this system does not support High dynamic range (HDR) at all.');
    end

    if ~isempty(find(mystrcmp(reqs, 'UseStaticHDRHack')))
        if ~(ismember(userstereomode, [0, 1, 4, 5, 11]) && IsLinux && ~IsWayland)
            error('PsychImaging: Requested task ''UseStaticHDRHack'' is incompatible with this system setup.');
        end

        % Add special requirements marker statichdrhack to signal to PsychHDR()
        % that this hack is requested:
        reqs{row, 5} = [reqs{row, 5} ' statichdrhack'];

        % Add the EnableNative10BitFramebuffer task if it doesn't exist already,
        % as our UseStaticHDRHack HDR support requires a 10 bpc native OpenGL
        % framebuffer:
        if isempty(find(mystrcmp(reqs, 'EnableNative10BitFramebuffer')))
            reqs = AddTask(reqs, 'General', 'EnableNative10BitFramebuffer');
        end
    end

    % Parse user provided parameters, and return imagingMode flags and other
    % flags/constraints needed for FinalizeConfiguration() stage:
    hdrImagingModeFlags = PsychHDR('GetClientImagingParameters', hdrArguments);
    imagingMode = mor(imagingMode, hdrImagingModeFlags);

    % Mark use of HDR:
    useHDR = 1;

    % Add the UseVulkanDisplay task if it doesn't exist already, as our HDR support
    % is implemented on top of Vulkan:
    floc = find(mystrcmp(reqs, 'UseVulkanDisplay'));
    if isempty(floc)
        reqs = AddTask(reqs, 'General', 'UseVulkanDisplay');
    end

    % No frame-sequential stereo yet:
    if ismember(userstereomode, [1, 11])
        error('PsychImaging: Requested task ''EnableHDR'' is incompatible with frame-sequential stereo mode %i.', userstereomode);
    end

    % Dual-window would turn into dual-stream:
    if ismember(userstereomode, [10])
        % Remap dual-window stereo to dual-stream stereo:
        stereoMode = 12;
    end
else
    % No HDR in use:
    useHDR = 0;
end

% Want to use the Vulkan/WSI display backend?
floc = find(mystrcmp(reqs, 'UseVulkanDisplay'));
if ~isempty(floc)
    [rows cols] = ind2sub(size(reqs), floc(1));
    row = rows(1);

    % Check if Vulkan/WSI based display output is supported at all, ie. if a
    % Vulkan loader library is installed, so PsychVulkanCore can be loaded and
    % linked, and if at least one Vulkan gpu is available on this system setup.
    % This will also perform driver init and Vulkan instance init, which is an
    % important thing to do as the very first thing before opening an onscreen
    % window (ie. before OpenGL context creation) on AMD + Windows-10, otherwise
    % the AMD proprietary OpenGL driver might crash due to AMD driver bugs!
    if ~PsychVulkan('Supported')
        % Failed/Unsupported.
        error('PsychImaging: Requested task ''UseVulkanDisplay'', but this system does not support Vulkan at all.');
    end

    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);

    % Special static HDR dual-display stereo hack requested and supported?
    if ~useHDR || isempty(find(mystrcmp(reqs, 'UseStaticHDRHack')))
        % No -> Standard operations: Add imaging mode flags for handing rendered images to Vulkan:
        imagingMode = mor(imagingMode, kPsychNeedFinalizedFBOSinks);

        if IsOSX
            imagingMode = mor(imagingMode, kPsychUseExternalSinkTextures);
        end

        % Mark full use of Vulkan:
        useVulkan = 1;
    else
        % We abuse Vulkan for HDR setup in UseStaticHDRHack, but not as primary
        % display backend:
        useVulkan = 0;
    end

    if ismember(userstereomode, [1, 11])
        error('PsychImaging: Requested task ''UseVulkanDisplay'' is incompatible with frame-sequential stereo mode %i.', userstereomode);
    end

    if ismember(userstereomode, [10])
        % Remap dual-window stereo to dual-stream stereo:
        stereoMode = 12;
    end
else
    useVulkan = 0;
end

% Display replication needed?
if ~isempty(find(mystrcmp(reqs, 'MirrorDisplayTo2ndOutputHead')))
    % Yes: Must use dual window output mode. This implies
    % kPsychNeedFastBackingStore, automatically set by Screen('OpenWindow')
    % itself, so no need to do it here.
    imagingMode = mor(imagingMode, kPsychNeedDualWindowOutput);
end

% Custom color correction for display wanted?
if ~isempty(find(mystrcmp(reqs, 'DisplayColorCorrection')))
    % Color correction in output chain?
    if (~isempty(find(mystrcmp(reqs, 'AllViews'))) || ~isempty(find(mystrcmp(reqs, 'FinalFormatting'))))
        % Yes. Need full pipeline in any case, ie fast backing store and output conversion:
        imagingMode = mor(imagingMode, kPsychNeedFastBackingStore, kPsychNeedOutputConversion);
    else
        % No. Fast backing store is enough, per-view chains will get enabled below for Left/RightView:
        imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    end
end

% Replication of left half of window into right half needed?
% This is used for a software implementation of mirror mode displays,
% e.g., in conjunction with desktop-spanning display mode on MS-Windows or
% in conjunction with a display splitter on a single output head:
if ~isempty(find(mystrcmp(reqs, 'MirrorDisplayToSingleSplitWindow')))
    % We simply request that window size is reported and handled as if the
    % window would be only half the width --> right half remains empty and
    % can be used as target for the cloning op of the left half.
    % This works even without imaging pipe enabled, only uses finalizer
    % blit chains:
    imagingMode = mor(imagingMode, kPsychNeedHalfWidthWindow);
end

% 16 bpc float framebuffers needed?
if ~isempty(find(mystrcmp(reqs, 'FloatingPoint16Bit')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    imagingMode = mor(imagingMode, kPsychNeed16BPCFloat);
end

% 32 bpc float framebuffers needed?
if ~isempty(find(mystrcmp(reqs, 'FloatingPoint32Bit')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    imagingMode = mor(imagingMode, kPsychNeed32BPCFloat);
end

if ~isempty(find(mystrcmp(reqs, 'FloatingPoint32BitIfPossible')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    imagingMode = mor(imagingMode, kPsychUse32BPCFloatAsap);
end

if ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++Output'))) || ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++OutputWithOverlay')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    imagingMode = mor(imagingMode, kPsychNeedOutputConversion);
    % The Mono++ formatter is icm aware - Capable of internal color
    % correction, but not setup here -- special case: Set flag to zero:
    ptb_outputformatter_icmAware = 0;
end

if ~isempty(find(mystrcmp(reqs, 'EnableGenericHighPrecisionLuminanceOutput'))) || ~isempty(find(mystrcmp(reqs, 'EnablePseudoGrayOutput')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    imagingMode = mor(imagingMode, kPsychNeedOutputConversion);

    % The Luminance LUT based formatter is icm aware - Capable of internal color correction:
    ptb_outputformatter_icmAware = 1;

    % Request 32bpc float FBO unless already a 16 bpc FBO or similar has
    % been explicitely requested:
    if ~bitand(imagingMode, kPsychNeed16BPCFloat) && ~bitand(imagingMode, kPsychUse32BPCFloatAsap) && ~bitand(imagingMode, kPsychNeed16BPCFixed)
        imagingMode = mor(imagingMode, kPsychNeed32BPCFloat);
    end
end

if ~isempty(find(mystrcmp(reqs, 'EnableVideoSwitcherSimpleLuminanceOutput'))) || ~isempty(find(mystrcmp(reqs, 'EnableVideoSwitcherCalibratedLuminanceOutput')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    imagingMode = mor(imagingMode, kPsychNeedOutputConversion);

    % The VideoSwitcher formatter is icm aware - Capable of internal color correction:
    ptb_outputformatter_icmAware = 1;

    % Request 32bpc float FBO unless already a 16 bpc FBO or similar has
    % been explicitely requested:
    if ~bitand(imagingMode, kPsychNeed16BPCFloat) && ~bitand(imagingMode, kPsychUse32BPCFloatAsap) && ~bitand(imagingMode, kPsychNeed16BPCFixed)
        imagingMode = mor(imagingMode, kPsychNeed32BPCFloat);
    end
end

if ~isempty(find(mystrcmp(reqs, 'EnableBits++Color++Output')))
    floc = find(mystrcmp(reqs, 'EnableBits++Color++Output'));
    [rows cols] = ind2sub(size(reqs), floc(1));
    row = rows(1);
    % Extract first parameter - This should be the colorConversionMode:
    colorConversionMode = reqs{row, 3};
    BitsPlusPlus('SetColorConversionMode', colorConversionMode);

    % These settings are mildly redundant, as the dedicated
    % OpenWindowColor++ code in the BitsPlusPlus.m helper file will do all
    % neccessary setup, especially deciding of kPsychNeedHalfWidthWindow is
    % needed or not:
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    imagingMode = mor(imagingMode, kPsychNeedOutputConversion);

    % The Color++ formatter is icm aware - Capable of internal color
    % correction, but not setup here -- special case: Set flag to zero:
    ptb_outputformatter_icmAware = 0;
end

% Request for native 10 bit per color component ARGB2101010 framebuffer,
% or native almost 11 bit per color component RGB111110 framebuffer?
if ~isempty(find(mystrcmp(reqs, 'EnableNative10BitFramebuffer'))) || ...
   ~isempty(find(mystrcmp(reqs, 'EnableNative11BitFramebuffer')))

    if ~useVulkan
        % Enable output formatter chain:
        imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
        imagingMode = mor(imagingMode, kPsychNeedOutputConversion);
    end

    % Request 32bpc float FBO unless already a 16 bpc FBO or similar has
    % been explicitely requested: In principle, a 16 bpc FBO would be
    % sufficient for a native 10 to 11 bpc framebuffer...
    if ~bitand(imagingMode, kPsychNeed16BPCFloat) && ~bitand(imagingMode, kPsychUse32BPCFloatAsap) && ~bitand(imagingMode, kPsychNeed16BPCFixed)
        imagingMode = mor(imagingMode, kPsychNeed32BPCFloat);
    end

    % The ATI 10/11bpc formatter is not yet icm aware - Incapable of internal color correction!
    % Additionally native 10/11 bpc framebuffers, e.g., on Fire-Series or NVidia cards also don't
    % have icm aware output formatting, so a 'false' setting would be mandatory. However, we leave
    % the setting at whatever it currently is, as it defaults to 0 / false anyway, and the PseudoGray
    % or EnableGenericHighPrecisionLuminanceOutput might have requested a 1 / true setting and we do
    % not want to override that. This to allow to stack such a perceptual precision boosting trick
    % on top of a 10/11 bpc framebuffer.
    % ptb_outputformatter_icmAware = 0;
end

% Request for native 16 bit per color component RGBA16161616 framebuffer?
if ~isempty(find(mystrcmp(reqs, 'EnableNative16BitFramebuffer')))
    % Request 32bpc float FBO unless already a 16 bpc fixed point FBO
    % has been explicitely requested. 16 bpc fixed point is obviously just
    % quite sufficient for 16 bpc linear output, 32 bpc float provides 23 bpc
    % effective linear precision in the meaningful output intensity range, so
    % leaves some numerical headroom for post processing and roundoff errors:
    if ~bitand(imagingMode, kPsychUse32BPCFloatAsap) && ~bitand(imagingMode, kPsychNeed16BPCFixed)
        imagingMode = mor(imagingMode, kPsychNeed32BPCFloat);
    end
end

% Request for dual display pipeline custom HDR system?
if ~isempty(find(mystrcmp(reqs, 'EnableDualPipeHDROutput')))
    % Enable imaging pipeline ...
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    % ... final device output formatter chain(s) ...
    imagingMode = mor(imagingMode, kPsychNeedOutputConversion);
    % ... and dual stream processing and output to two displays ...
    imagingMode = mor(imagingMode, kPsychNeedDualWindowOutput);

    % Request 32bpc float FBO unless already a 16 bpc FBO or similar has
    % been explicitely requested:
    if ~bitand(imagingMode, kPsychNeed16BPCFloat) && ~bitand(imagingMode, kPsychUse32BPCFloatAsap) && ~bitand(imagingMode, kPsychNeed16BPCFixed)
        imagingMode = mor(imagingMode, kPsychNeed32BPCFloat);
    end

    % The dual-pipeline HDR output formatter is not yet icm aware -
    % Incapable of internal color correction. Well, technically it is, but
    % that code-path is disabled for now. It is probably computationally
    % more efficient to perform one generic ICM pass on the input buffer
    % and then feed into the formatters for the two pipes instead of
    % letting each pipe's formatter apply the same color correction, ie.,
    % do the same work twice. This needs to be found out in the future. For
    % now we go for the simple solution:
    ptb_outputformatter_icmAware = 0;
end

if ~isempty(find(mystrcmp(reqs, 'LeftView'))) || ~isempty(find(mystrcmp(reqs, 'RightView')))
    % Specific eye channel requested: Need a stereo display mode.
    if stereoMode == -1
        % None set yet. Just channel the request to the caller:
        stereoMode = -2;
    end

    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);

    % Also need image processing stage, because only it can provide
    % separate processing for both eyes:
    imagingMode = mor(imagingMode, kPsychNeedImageProcessing);
else
    % Not a single eye specific command requested: Check if there's any
    % other spec that would require the image processing stage:
    % Any command that applies to 'AllViews' naturally needs the image
    % processing:
    if ~isempty(find(mystrcmp(reqs, 'AllViews')))
        imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
        imagingMode = mor(imagingMode, kPsychNeedImageProcessing);
    end
end

% Image processing stage needed?
if bitand(imagingMode, kPsychNeedImageProcessing)
    % Yes. How many commands per chain?
    nrslots = max(length(find(mystrcmp(reqs, 'LeftView'))), length(find(mystrcmp(reqs, 'RightView'))));
    nrslots = nrslots + length(find(mystrcmp(reqs, 'AllViews')));

    % More than one slot per chain? Otherwise we use the default
    % single-pass chain:
    if nrslots > 1
        % More than two slots per chain?
        if nrslots > 2
            % Need full blown multistage chain:
            imagingMode = mor(imagingMode, kPsychNeedMultiPass);
        else
            % Only two slots. More efficient dual-pass chain is sufficient:
            imagingMode = mor(imagingMode, kPsychNeedDualPass);
        end
    end
end

% Final output formatting stage needed?
if ~isempty(find(mystrcmp(reqs, 'FinalFormatting')))
    imagingMode = mor(imagingMode, kPsychNeedFastBackingStore);
    imagingMode = mor(imagingMode, kPsychNeedOutputConversion);
end

% Native 16 bit per color component RGBA16F floating point framebuffer requested?
if ~isempty(find(mystrcmp(reqs, 'EnableNative16BitFloatingPointFramebuffer')))
    % Yes. Pipeline active?
    if bitand(imagingMode, kPsychNeedFastBackingStore)
        % Native 16 bpc float doesn't need imaging pipeline by itself, but if the
        % pipeline is enabled then it needs to provide sufficient >= 16 bpc float
        % precision throughout the pipeline. Make sure we request that.

        % Request 32bpc float FBO unless already a 16 bpc FBO or similar has
        % been explicitely requested: In principle, a 16 bpc FBO would be
        % sufficient for a native 16 bpc float framebuffer...
        if ~bitand(imagingMode, kPsychNeed16BPCFloat) && ~bitand(imagingMode, kPsychUse32BPCFloatAsap) && ~bitand(imagingMode, kPsychNeed16BPCFixed)
            imagingMode = mor(imagingMode, kPsychNeed32BPCFloat);
        end

        % Not needed, as it is default: ptb_outputformatter_icmAware = 0;
    end
end

% Validate if everything is right for use of HDR displays:
if useHDR
    % For HDR we absolutely need a floating point virtual framebuffer on the
    % input side, ie. for user-scripts to draw into, as most HDR color values
    % will be far outside (exceeding) the unorm range 0 - 1:
    if ~bitand(imagingMode, kPsychNeed16BPCFloat) && ~bitand(imagingMode, kPsychUse32BPCFloatAsap) && ~bitand(imagingMode, kPsychNeed32BPCFloat)
        % No floating point input drawbuffer format requested yet.
        % Play it safe and request 32 bit single precision float to prefer
        % maximum accuracy over performance:
        imagingMode = mor(imagingMode, kPsychNeed32BPCFloat);
    end

    % Assume that all HDR output formatters are ICM aware, ie. allow/require linking
    % an ICM shader into them for colorspace conversion, clamping etc., for efficient
    % single-pass colorspace conversion and HDR OETF mapping:
    ptb_outputformatter_icmAware = 1;
end

% Support for fast offscreen windows (aka FBO backed offscreen windows)
% needed?
if ~isempty(find(mystrcmp(reqs, 'UseFastOffscreenWindows')))
    % Need fast offscreen windows. They are included if any non-zero imagingMode
    % is set, so we only request'em if imagingMode is still zero:
    if imagingMode == 0
        imagingMode = kPsychNeedFastOffscreenWindows;
    end
end

return;

% End of FinalizeConfiguration subroutine.

% PostConfiguration is called after the onscreen window is open: Performs
% actual pipeline setup of the hook chains:
function rc = PostConfiguration(reqs, win, clearcolor, slavewin)
global ptb_outputformatter_icmAware;
global ptb_MirrorOverlayWindows;
global GL;
global ptb_geometry_inverseWarpMap;
global psych_gpgpuapi; %#ok<NUSED>
% Default requested colormode: Set by PsychDefaultSetup(), if at all.
global psych_default_colormode;

% At least two video outputs scanning out in dual-display side-by-side configuration?
global isASideBySideConfig;
global screenRestoreCmd;

if isempty(GL)
    % Perform minimal OpenGL init, so we can call OpenGL commands and use
    % GL constants. We do not activate a full 3D rendering context:
    InitializeMatlabOpenGL([], [], 1);
end

% Identity CLUT in graphics hardware required?
needsIdentityCLUT = 0;

% Should dithering be disabled if 'needsIdentityCLUT'?
% By default we disable in such a case:
disableDithering = 1;

% 0.0 - 1.0 colorrange without color clamping required?
needsUnitUnclampedColorRange = 0;
applyAlsoToMakeTexture = [];

% Number of used slots in left- and right processing chain:
leftcount = 0;
rightcount = 0;
outputcount = 0;
outputcount0 = 0;
outputcount1 = 0;

% Flags for horizontal/vertical flip operations:
leftUDFlip = 0;
rightUDFlip = 0;
leftLRFlip = 0;
rightLRFlip = 0;

% Stereomode?
winfo = Screen('GetWindowInfo', win);
[winwidth, winheight] = InterBufferSize(win);

% Setup inverse warp map matrices for this window handle:
ptb_geometry_inverseWarpMap{win} = [];
ptb_geometry_inverseWarpMap{win}.gx = 1;
ptb_geometry_inverseWarpMap{win}.gy = 1;
ptb_geometry_inverseWarpMap{win}.mx = winwidth;
ptb_geometry_inverseWarpMap{win}.my = winheight;

if ismember(winfo.StereoMode, [2,3])
    ptb_geometry_inverseWarpMap{win}.gy = 2;
end

% Determine early if Vulkan display backend is to be used:
floc = find(mystrcmp(reqs, 'UseVulkanDisplay'));
if ~isempty(floc)
    useVulkan = 1;

    % Default to SDR standard dynamic range display mode:
    vulkanHDRMode = 0;

    % Default color precision: In SDR mode this is standard 8 bpc RGBA8 unorm
    % fixed point color precision. In vulkanHDRMode > 0 it may get automatically
    % overriden to something higher precision:
    vulkanColorPrecision = 0;

    % Default color space: Overriden in vulkanHDRMode > 0.
    vulkanColorSpace = 0;

    % Default color pixel format: Overriden in vulkanHDRMode > 0 as part of
    % override of vulkanColorPrecision, or explicit override:
    vulkanColorFormat = 0;
else
    useVulkan = 0;
end

% Determine early if HDR display is to be used:
floc = find(mystrcmp(reqs, 'EnableHDR'));
if ~isempty(floc)
    % Assign all user provided arguments in hdrArguments for later use:
    [rows cols] = ind2sub(size(reqs), floc(1));
    row = rows(1);
    hdrArguments = reqs(row, :);

    % Get HDR setup parameters for imaging pipeline and Vulkan:
    [useVulkan, vulkanHDRMode, vulkanColorPrecision, vulkanColorSpace, vulkanColorFormat] = PsychHDR('GetVulkanHDRParameters', win, hdrArguments);

    %if ~isempty(find(mystrcmp(reqs, 'UseStaticHDRHack')))
    % Better with useVulkan to suppress our own 10/11/16 bpc output formatters atm. useVulkan = 1;
    %end

    % Mark HDR as in use:
    useHDR = 1;
else
    hdrArguments = [];
    useHDR = 0;
end

% --- First action in pipe is a horizontal- or vertical flip, if any ---

% Any flip horizontal requested?
floc = find(mystrcmp(reqs, 'FlipHorizontal'));
if ~isempty(floc)
    % Which channel?
    for x=floc
        [rows cols]= ind2sub(size(reqs), x);
        for row=rows'
            if mystrcmp(reqs{row, 1}, 'LeftView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % LeftView horizontal flip.
                leftLRFlip = 1;
            end

            if mystrcmp(reqs{row, 1}, 'RightView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % LeftView horizontal flip.
                rightLRFlip = 1;
            end
        end
    end
end

% Any flip vertical requested?
floc = find(mystrcmp(reqs, 'FlipVertical'));
if ~isempty(floc)
    % Which channel?
    for x=floc
        [rows cols]= ind2sub(size(reqs), x);
        for row=rows'
            if mystrcmp(reqs{row, 1}, 'LeftView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % LeftView vertical flip.
                leftUDFlip = 1;
            end

            if mystrcmp(reqs{row, 1}, 'RightView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % LeftView vertical flip.
                rightUDFlip = 1;
            end
        end
    end
end

% Left channel flipping needed?
if leftLRFlip || leftUDFlip
    % Yes.
    sx = 1;
    ox = 0;
    sy = 1;
    oy = 0;

    if leftLRFlip
        sx = -1;
        ox = RectWidth(InterBufferRect(win));
        hv = winwidth-1:-1:0;
    else
        hv = 0:winwidth-1;
    end

    if leftUDFlip
        sy = -1;
        oy = RectHeight(InterBufferRect(win));
        vv = winheight-1:-1:0;
    else
        vv = 0:winheight-1;
    end

    % Enable left imaging chain:
    Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');
    % Append blitter for LR/UD flip:
    Screen('HookFunction', win, 'AppendBuiltin', 'StereoLeftCompositingBlit', 'Builtin:IdentityBlit', sprintf('Offset:%i:%i:Scaling:%f:%f', ox, oy, sx, sy));
    leftcount = leftcount + 1;

    clear curmap;
    [xg,yg] = meshgrid(hv, vv);
    curmap(:,:,1) = xg;
    curmap(:,:,2) = yg;
    ptb_geometry_inverseWarpMap{win}.('LeftView') = int16(curmap);
end

if winfo.StereoMode > 0
    % Stereomode enabled: Need to possibly handle right channel as
    % well. In mono-mode there would be only a left channel...

    % Right channel flipping needed?
    if rightLRFlip || rightUDFlip
        % Yes.
        sx = 1;
        ox = 0;
        sy = 1;
        oy = 0;

        if rightLRFlip
            sx = -1;
            ox = RectWidth(InterBufferRect(win));
            hv = winwidth-1:-1:0;
        else
            hv = 0:winwidth-1;
        end

        if rightUDFlip
            sy = -1;
            oy = RectHeight(InterBufferRect(win));
            vv = winheight-1:-1:0;
        else
            vv = 0:winheight-1;
        end

        % Enable right imaging chain:
        Screen('HookFunction', win, 'Enable', 'StereoRightCompositingBlit');
        % Append blitter for LR/UD flip:
        Screen('HookFunction', win, 'AppendBuiltin', 'StereoRightCompositingBlit', 'Builtin:IdentityBlit', sprintf('Offset:%i:%i:Scaling:%f:%f', ox, oy, sx, sy));
        rightcount = rightcount + 1;

        clear curmap;
        [xg,yg] = meshgrid(hv, vv);
        curmap(:,:,1) = xg;
        curmap(:,:,2) = yg;
        ptb_geometry_inverseWarpMap{win}.('RightView') = int16(curmap);
    end
end

% --- End of the flipping stuff ---


% --- Implementation of CLUT animation via clut remapping of colors ---
floc = find(mystrcmp(reqs, 'EnableCLUTMapping'));
% Is a display mode on a CRS Bits+/Bits# or VPixx DataPixx/ViewPixx/ProPixx requested which requires use
% and setup of the devices hardware CLUT? If so we must turn 'EnableCLUTMapping' into a no-op, as it
% would clash with the hardware clut update - and is also superseded by it. Detect the namestrings of
% Bits++ CLUT palette display mode and Mono++ CLUT overlay palette mode. These Bits+ namestrings also
% cover VPixx devices due to the remapping of VPixx names into CRS reqs:
if ~isempty(find(mystrcmp(reqs, 'EnableBits++Bits++Output'))) || ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++OutputWithOverlay')))
    % Yep. We must no-op this 'EnableCLUTMapping' request:
    floc = [];
end
if ~isempty(floc)
    % Which channel?
    for x=floc
        [rows cols]= ind2sub(size(reqs), x);
        for row=rows'
            % Extract first parameter - This should be the number of clut slots:
            nClutSlots = reqs{row, 3};

            % Default to 256 slot clut, as most commonly used:
            if isempty(nClutSlots)
                nClutSlots = 256;
            end

            if ~isnumeric(nClutSlots)
                sca;
                error('PsychImaging: Number of clut slots parameter for ''EnableCLUTMapping'' missing or not of numeric type!');
            end

            % Extract high precision flag:
            highprec = reqs{row, 4};
            if isempty(highprec)
                highprec = 0;
            end

            % Use our reformatter shader for mapping RGB indices to RGB
            % triplets.

            % Load shader:
            pgshader = LoadGLSLProgramFromFiles('RGBMultiLUTLookupCombine_FormattingShader', 1);

            % Init the shader:
            glUseProgram(pgshader);

            % Assign mapping of input image and clut to texture units:
            glUniform1i(glGetUniformLocation(pgshader, 'Image'), 0);
            glUniform1i(glGetUniformLocation(pgshader, 'CLUT'),  1);

            % Assign number of clut slots to use:
            glUniform1f(glGetUniformLocation(pgshader, 'Prescale'), nClutSlots - 1);
            glUseProgram(0);

            % Use helper routine to build a proper RGBA lookup texture:
            pglutid = PsychHelperCreateRemapCLUT(0, nClutSlots, highprec);
            pgconfig = sprintf('TEXTURERECT2D(1)=%i', pglutid);

            % Setup the callback function which is only called if the clut
            % texture needs to be updated because
            % Screen('LoadNormalizedGammatable', ..., 2); was called to
            % provide a new clut. We attach this to the left image
            % processing chain, as this chain is almost always used anyway.
            % It needs to execute only once per flip, as it updates state
            % global to all views (in a stereo setup):
            rclutcmd = sprintf('PsychHelperCreateRemapCLUT(1, %i, IMAGINGPIPE_GAMMATABLE);', pglutid);
            Screen('HookFunction', win, 'AppendMFunction', 'StereoLeftCompositingBlit', 'Upload new clut into shader callback', rclutcmd);

            % Enable left chain unconditionally, so the above clut setup
            % code gets executed:
            Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');

            % Attach shaders and slots to proper processing chain.
            % These perform the clut color conversion blit of each input
            % image into a transformed output image. They're executed at
            % each flip, irrespective if the clut changed or not:
            if mystrcmp(reqs{row, 1}, 'LeftView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % Need to attach to left view:
                if leftcount > 0
                    % Need a bufferflip command:
                    Screen('HookFunction', win, 'AppendBuiltin', 'StereoLeftCompositingBlit', 'Builtin:FlipFBOs', '');
                end
                Screen('HookFunction', win, 'AppendShader', 'StereoLeftCompositingBlit', 'CLUT image transformation shader', pgshader, pgconfig);
                Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');
                leftcount = leftcount + 1;
            end

            if mystrcmp(reqs{row, 1}, 'RightView') || (mystrcmp(reqs{row, 1}, 'AllViews') && winfo.StereoMode > 0)
                % Need to attach to right view:
                if rightcount > 0
                    % Need a bufferflip command:
                    Screen('HookFunction', win, 'AppendBuiltin', 'StereoRightCompositingBlit', 'Builtin:FlipFBOs', '');
                end
                Screen('HookFunction', win, 'AppendShader', 'StereoRightCompositingBlit', 'CLUT image transformation shader', pgshader, pgconfig);
                Screen('HookFunction', win, 'Enable', 'StereoRightCompositingBlit');
                rightcount = rightcount + 1;
            end
        end
    end
end
% --- End of CLUT animation via clut remapping of colors ---


% --- Addition of offsets / scales etc. to input image ---
floc = find(mystrcmp(reqs, 'AddOffsetToImage'));
if ~isempty(floc)
    % Which channel?
    for x=floc
        [rows cols]= ind2sub(size(reqs), x);
        for row=rows'
            % Extract first parameter - This should be the offset:
            PixelOffset = reqs{row, 3};

            if isempty(PixelOffset) || ~isnumeric(PixelOffset)
                sca;
                error('PsychImaging: Parameter for ''AddOffsetToImage'' missing or not of numeric type!');
            end

            % Further (optional) parameters passed?
            % 2nd parameter, if any, would be a gain value to apply before
            % applying the PixelOffset:
            PixelGain = reqs{row, 4};
            if isempty(PixelGain)
                % No such flag: Default to 1:
                PixelGain = 1;
            else
                if ~isnumeric(PixelGain)
                    sca;
                    error('PsychImaging: Optional Gain-Parameter for ''AddOffsetToImage'' not of numeric type!');
                end
            end

            % 3rd parameter, if any, would be an Offset value to apply before
            % applying the gain:
            PixelPreOffset = reqs{row, 5};
            if isempty(PixelPreOffset)
                % No such flag: Default to 0:
                PixelPreOffset = 0;
            else
                if ~isnumeric(PixelPreOffset)
                    sca;
                    error('PsychImaging: Optional "Offset before Gain"- PrescaleParameter for ''AddOffsetToImage'' not of numeric type!');
                end
            end

            % Load and build shader:
            shader = LoadGLSLProgramFromFiles('ScaleAndBiasShader', 1);

            % Init the shader: Assign mapping of input image and offsets, gains:
            glUseProgram(shader);

            glUniform1i(glGetUniformLocation(shader, 'Image'), 0);
            glUniform1f(glGetUniformLocation(shader, 'postscaleoffset'), PixelOffset);
            glUniform1f(glGetUniformLocation(shader, 'prescaleoffset'), PixelPreOffset);
            glUniform1f(glGetUniformLocation(shader, 'scalefactor'), PixelGain);

            glUseProgram(0);

            % Ok, 'gld' should contain a valid OpenGL display list for
            % geometry correction. Attach proper shader to proper chain:
            if mystrcmp(reqs{row, 1}, 'LeftView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % Need to attach to left view:
                if leftcount > 0
                    % Need a bufferflip command:
                    Screen('HookFunction', win, 'AppendBuiltin', 'StereoLeftCompositingBlit', 'Builtin:FlipFBOs', '');
                end
                Screen('HookFunction', win, 'AppendShader', 'StereoLeftCompositingBlit', 'ScaleAndOffsetShader', shader);
                Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');
                leftcount = leftcount + 1;
            end

            if mystrcmp(reqs{row, 1}, 'RightView') || (mystrcmp(reqs{row, 1}, 'AllViews') && winfo.StereoMode > 0)
                % Need to attach to right view:
                if rightcount > 0
                    % Need a bufferflip command:
                    Screen('HookFunction', win, 'AppendBuiltin', 'StereoRightCompositingBlit', 'Builtin:FlipFBOs', '');
                end
                Screen('HookFunction', win, 'AppendShader', 'StereoRightCompositingBlit', 'ScaleAndOffsetShader', shader);
                Screen('HookFunction', win, 'Enable', 'StereoRightCompositingBlit');
                rightcount = rightcount + 1;
            end

            if mystrcmp(reqs{row, 1}, 'FinalFormatting')
                % Need to attach to final formatting:
                if outputcount > 0
                    % Need a bufferflip command:
                    Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
                end
                Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', 'ScaleAndOffsetShader', shader);
                Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
                outputcount = outputcount + 1;
            end
        end
    end
end
% --- End of addition of offsets / scales etc. to input image ---


% --- Geometry correction via warped blit ---
floc = find(mystrcmp(reqs, 'GeometryCorrection'));
if ~isempty(floc)
    % Which channel?
    for x=floc
        [rows cols]= ind2sub(size(reqs), x);
        for row=rows'
            % Extract first parameter - This should be the name of a
            % calibration file:
            calibfilename = reqs{row, 3};

            if isempty(calibfilename)
                sca;
                error('PsychImaging: Parameter for ''GeometryCorrection'' missing!');
            end

            % Is 'calibfilename' a function handle or a final warpstruct?
            if (~isstruct(calibfilename) && ~ischar(calibfilename)) || ...
               (isstruct(calibfilename) && isfield(calibfilename, 'gld') && isfield(calibfilename, 'glsl'))
                % Functionhandle or final warpstruct passed: This
                % assignment will either assign the warpstruct, or call the
                % function referenced by the functionhandle and assign the
                % returned warpstruct:
                if ~isstruct(calibfilename)
                    [warpstruct, filterMode] = calibfilename();
                else
                    warpstruct = calibfilename;
                    filterMode = ':Bilinear';
                end
            else
                % Either calibration input parameter struct, or filename of
                % calibration file: Just pass it to CreateDisplayWarp(),
                % after some parameter validation:
                if ischar(calibfilename) && ~exist(calibfilename, 'file')
                    sca;
                    error('PsychImaging: Passed an argument to ''GeometryCorrection'' which is not a valid name of an accessible calibration file!');
                end

                % Filename or calibstruct valid. Further (optional) parameters passed?
                % 2nd parameter, if any, would be a 'visualize' flag that
                % asks for plotting of some calibration info and additional
                % output to the console:
                showCalibOutput = reqs{row, 4};
                if isempty(showCalibOutput)
                    % No such flag: Default to "silence":
                    showCalibOutput = 0;
                end

                % Additional parameters provided? Pass 'em along. Currently
                % defined are up to additional 6 parameters 5 to 10. These
                % default to empty if not provided by user-code.

                % Use helper function to read the calibration file or
                % parameter struct and build a proper warp-function:
                [warpstruct, filterMode] = CreateDisplayWarp(win, calibfilename, showCalibOutput, reqs{row, 5:10});
            end

            % Is it a display list handle?
            if ~isempty(warpstruct.gld)
                % This must be a display list handle for display list
                % blitting, potentially with an additional GLSL shader
                % attached:
                gld = warpstruct.gld;
                if ~glIsList(gld)
                    % Game over:
                    sca;
                    error('PsychImaging: Passed a handle to ''GeometryCorrection'' which is not a valid OpenGL display list!');
                end

                if ~isempty(warpstruct.glsl)
                    glsl = warpstruct.glsl;
                else
                    glsl = 0;
                end

                % Ok, 'gld' should contain a valid OpenGL display list for
                % geometry correction.

                % Before we setup the image warping ops for real in the pipeline, we
                % do a "cold run" to compute a 2D reverse lookup table that allows to
                % map warped 2D screen positions back to their originating pre-warp pixels.
                % This is useful, e.g., if one wants to map 2D mouse click
                % positions on the geometry corrected display back to the
                % originating pixel positions in the uncorrected stimulus
                % image.
                %
                % This works by creating a float texture whose texels
                % color-code their spatial (x,y) locations in the R and G
                % channels, then warping this texture with the same
                % operations that the GPU will apply to the stimulus
                % images, then reading back the warp-blitted texture into a
                % 2-layer 2D matrix, where layer 1 (former red channel)
                % encodes originating x-position of each "pixel", layer 2
                % encodes y-position, Undefined positions are mapped to (0,0):

                % At least 32 bpc float or 16 bit snorm textures/fbo's
                % supported? Otherwise this is a no-go:
                if (winfo.GLSupportsTexturesUpToBpc >= 32) || ~isempty(strfind(glGetString(GL.EXTENSIONS), '_texture_snorm'))
                    % Yes.

                    % Check if previous code already defined some inverse
                    % mapping:
                    if ~isempty(ptb_geometry_inverseWarpMap{win}) && isfield(ptb_geometry_inverseWarpMap{win}, reqs{row, 1})
                        % Yes: Extract it and use it as starting point for
                        % geometry inverse mapping:
                        premap = double(ptb_geometry_inverseWarpMap{win}.(reqs{row, 1}));
                        xg = premap(:,:,1);
                        yg = premap(:,:,2);
                    else
                        % No: Create a default identity mapping as starting
                        % point:
                        [xg,yg] = meshgrid(0:winwidth-1, 0:winheight-1);
                    end

                    % Need to use snorm 16 bit textures because 32 bpc
                    % float textures unavailable?
                    invmap_needs_snorm = (winfo.GLSupportsTexturesUpToBpc < 32);

                    % We always normalize to range 0..1, so it works for
                    % both floating point textures and 16 bit snorm
                    % textures:
                    inmap = zeros(winheight, winwidth, 3);
                    inmap(:,:,1) = xg / winwidth;
                    inmap(:,:,2) = yg / winheight;

                    if invmap_needs_snorm
                        % Need to use 16 bit snorm textures. We request 16
                        % bit floating point precision on this hw that
                        % doesn't support it, but does support 16 bit
                        % snorm. Screen() will choose 16 bit snorm as
                        % fallback, so we get what we want and can properly
                        % process mappings for up to 32k x 32k pixels aka 1
                        % Gigapixel:
                        premaptex = Screen('MakeTexture', win, inmap, [], [], 1);
                        postmaptex = Screen('OpenOffscreenWindow', win, 0, Screen('Rect', premaptex), 64);
                    else
                        % We have 32 bpc float texture support: Use it.
                        premaptex = Screen('MakeTexture', win, inmap, [], [], 2);
                        postmaptex = Screen('OpenOffscreenWindow', win, 0, Screen('Rect', premaptex), 128);
                    end
                    warpoperator = CreateGLOperator(win);
                    AddImageUndistortionToGLOperator(warpoperator, premaptex, warpstruct);
                    postmaptex = Screen('TransformTexture', premaptex, warpoperator, [], postmaptex);
                    glerr = glGetError;
                    if glerr
                        % We get this error on some NVidia binary blob graphics driver on Linux, e.g., v295.49. Swallow it, it seems to cause no consequences:
                        fprintf('PsychImaging: GeometryCorrection: Spurious benign gl error [%s] after computing postmap texture detected.\n', gluErrorString(glerr));
                    end
                    curmap = Screen('GetImage', postmaptex, [], [], 1, 3);
                    Screen('Close', [premaptex, postmaptex, warpoperator]);
                    curmap(:,:,1) = curmap(:,:,1) * winwidth;
                    curmap(:,:,2) = curmap(:,:,2) * winheight;
                    curmap = round(curmap(:,:,1:2));

                    % Assign inverse warp mapping tables for selected view. We
                    % assume that 16 bit signed integer is enough - Can cope
                    % with a framebuffer of up to 32768 * 32768 pixels.
                    % ptb_geometry_inverseWarpMap{} is a global variable shared
                    % with the RemapMouse() functions that uses these mapping
                    % matrices:
                    ptb_geometry_inverseWarpMap{win}.(reqs{row, 1}) = int16(curmap);
                else
                    % No: Cannot create remap textures at required
                    % precision, inverse mapping won't work:
                    fprintf('PsychImaging GeometryCorrection:Warning: GPU does not support features needed for RemapMouse() command.\n');
                end

                % Setup imaging pipeline - Attach proper blitters to proper chains:
                if mystrcmp(reqs{row, 1}, 'LeftView') || mystrcmp(reqs{row, 1}, 'AllViews')
                    % Need to setup left view warp:
                    if leftcount > 0
                        % Need a bufferflip command:
                        Screen('HookFunction', win, 'AppendBuiltin', 'StereoLeftCompositingBlit', 'Builtin:FlipFBOs', '');
                    end

                    % Must clear target buffer, because a geometrically
                    % warped blit for geometry correction may not cover the
                    % whole buffer area, and "uninitialized pixel trash"
                    % may shine through otherwise:
                    Screen('Hookfunction', win, 'AppendMFunction', 'StereoLeftCompositingBlit', 'Clear target buffer', 'glClear(16384);');

                    if glsl
                        Screen('HookFunction', win, 'AppendShader', 'StereoLeftCompositingBlit', 'GeometricWarpShader', glsl, sprintf('Blitter:DisplayListBlit:Handle:%i%s', gld, filterMode));
                    else
                        Screen('HookFunction', win, 'AppendBuiltin', 'StereoLeftCompositingBlit', 'Builtin:IdentityBlit', sprintf('Blitter:DisplayListBlit:Handle:%i%s', gld, filterMode));
                    end
                    Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');
                    leftcount = leftcount + 1;
                end

                if mystrcmp(reqs{row, 1}, 'RightView') || (mystrcmp(reqs{row, 1}, 'AllViews') && winfo.StereoMode > 0)
                    % Need to setup right view warp:
                    if rightcount > 0
                        % Need a bufferflip command:
                        Screen('HookFunction', win, 'AppendBuiltin', 'StereoRightCompositingBlit', 'Builtin:FlipFBOs', '');
                    end

                    Screen('Hookfunction', win, 'AppendMFunction', 'StereoRightCompositingBlit', 'Clear target buffer', 'glClear(16384);');

                    if glsl
                        Screen('HookFunction', win, 'AppendShader', 'StereoRightCompositingBlit', 'GeometricWarpShader', glsl, sprintf('Blitter:DisplayListBlit:Handle:%i%s', gld, filterMode));
                    else
                        Screen('HookFunction', win, 'AppendBuiltin', 'StereoRightCompositingBlit', 'Builtin:IdentityBlit', sprintf('Blitter:DisplayListBlit:Handle:%i%s', gld, filterMode));
                    end
                    Screen('HookFunction', win, 'Enable', 'StereoRightCompositingBlit');
                    rightcount = rightcount + 1;
                end

                if mystrcmp(reqs{row, 1}, 'FinalFormatting')
                    % Need to setup final formatting warp:
                    if outputcount > 0
                        % Need a bufferflip command:
                        Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
                    end

                    Screen('Hookfunction', win, 'AppendMFunction', 'FinalOutputFormattingBlit', 'Clear target buffer', 'glClear(16384);');

                    if glsl
                        Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', 'GeometricWarpShader', glsl, sprintf('Blitter:DisplayListBlit:Handle:%i%s', gld, filterMode));
                    else
                        Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:IdentityBlit', sprintf('Blitter:DisplayListBlit:Handle:%i%s', gld, filterMode));
                    end
                    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
                    outputcount = outputcount + 1;
                end
            else
                % Game over:
                sca;
                error('PsychImaging: Passed a handle for a not yet implemented display undistortion method!');
            end
        end
    end
end
% --- End of geometry correction via warped blit ---

% --- Interleaved line stereo wanted? ---
if ~isempty(find(mystrcmp(reqs, 'InterleavedLineStereo')))
    % Yes: Load and setup compositing shader.
    shader = LoadGLSLProgramFromFiles('InterleavedLineStereoShader', 1);

    floc = find(mystrcmp(reqs, 'InterleavedLineStereo'));
    [rows cols]= ind2sub(size(reqs), floc);
    % Extract first parameter - This should be the mapping of odd- and even
    % lines: 0 = even lines == left image, 1 = even lines == right image.
    startright = reqs{rows, 3};

    if startright~=0 && startright~=1
        sca;
        error('PsychImaging: The "startright" parameter must be zero or one!');
    end

    % Init the shader: Assign mapping of left- and right image:
    glUseProgram(shader);
    glUniform1i(glGetUniformLocation(shader, 'Image1'), 1-startright);
    glUniform1i(glGetUniformLocation(shader, 'Image2'), startright);

    glUniform2f(glGetUniformLocation(shader, 'Offset'), 0, 0);
    glUseProgram(0);

    % Reset compositor chain: It got initialized inside Screen() with an
    % unsuitable shader for our purpose:
    Screen('HookFunction', win, 'Reset', 'StereoCompositingBlit');

    % Append our new shader and enable chain:
    Screen('HookFunction', win, 'AppendShader', 'StereoCompositingBlit', 'StereoCompositingShaderInterleavedLineStereo', shader, 'Blitter:IdentityBlit:Offset:0:0:Scaling:1.0:2.0');
    Screen('HookFunction', win, 'Enable', 'StereoCompositingBlit');

    % Correct mouse position via proper gain:
    ptb_geometry_inverseWarpMap{win}.gy = ptb_geometry_inverseWarpMap{win}.gy * 0.5;
end
% --- End of interleaved line stereo setup code ---

% --- Interleaved column stereo wanted? ---
if ~isempty(find(mystrcmp(reqs, 'InterleavedColumnStereo')))
    % Yes: Load and setup compositing shader.
    shader = LoadGLSLProgramFromFiles('InterleavedColumnStereoShader', 1);

    floc = find(mystrcmp(reqs, 'InterleavedColumnStereo'));
    [rows cols]= ind2sub(size(reqs), floc);
    % Extract first parameter - This should be the mapping of odd- and even
    % columns: 0 = even cols == left image, 1 = even cols == right image.
    startright = reqs{rows, 3};

    if startright~=0 && startright~=1
        sca;
        error('PsychImaging: The "startright" parameter must be zero or one!');
    end

    % Init the shader: Assign mapping of left- and right image:
    glUseProgram(shader);
    glUniform1i(glGetUniformLocation(shader, 'Image1'), 1-startright);
    glUniform1i(glGetUniformLocation(shader, 'Image2'), startright);

    glUniform2f(glGetUniformLocation(shader, 'Offset'), 0, 0);
    glUseProgram(0);

    % Reset compositor chain: It got initialized inside Screen() with an
    % unsuitable shader for our purpose:
    Screen('HookFunction', win, 'Reset', 'StereoCompositingBlit');

    % Append our new shader and enable chain:
    Screen('HookFunction', win, 'AppendShader', 'StereoCompositingBlit', 'StereoCompositingShaderInterleavedColumnStereo', shader, 'Blitter:IdentityBlit:Offset:0:0:Scaling:2.0:1.0');
    Screen('HookFunction', win, 'Enable', 'StereoCompositingBlit');

    % Correct mouse position via proper gain:
    ptb_geometry_inverseWarpMap{win}.gx = ptb_geometry_inverseWarpMap{win}.gx * 0.5;
end
% --- End of interleaved column stereo setup code ---

% --- SideBySideCompressedStereo wanted? ---
if ~isempty(find(mystrcmp(reqs, 'SideBySideCompressedStereo')))
    % Yes: Call external setup routine with its default parameters to
    % modify our default "stereomode 2" top-bottom compressed stereo
    % shader, which was automatically generated by Screen('Openwindow'),
    % into a left-right side-by-side compressed shader.
    SetCompressedStereoSideBySideParameters(win);

    % Correct mouse position via proper gain:
    % Need to apply a 2x gain to horizontal cursor position to compensate
    % for horizontal compression...
    ptb_geometry_inverseWarpMap{win}.gx = ptb_geometry_inverseWarpMap{win}.gx * 2;
    % ... need to undo the 2x gain automatically applied at the top of this
    % function when a stereomode of 2 is used, as we do to implement our
    % stereo method:
    ptb_geometry_inverseWarpMap{win}.gy = ptb_geometry_inverseWarpMap{win}.gy / 2;
end
% --- End of SideBySideCompressedStereo setup code ---

% --- "Mouse" remapping needed for half-width Color++ or C48 mode? ---
if ~isempty(find(mystrcmp(reqs, 'EnableBits++Color++Output')))
    floc = find(mystrcmp(reqs, 'EnableBits++Color++Output'));
    [rows cols] = ind2sub(size(reqs), floc(1));
    row = rows(1);
    % Extract first parameter - This should be the colorConversionMode:
    colorConversionMode = reqs{row, 3};

    % Only mode 0 needs remapping:
    if colorConversionMode == 0
        % Correct mouse position via proper gain:
        ptb_geometry_inverseWarpMap{win}.gx = ptb_geometry_inverseWarpMap{win}.gx * 0.5;
    end
end

% --- Custom color correction for display wanted? ---
%
% This *MUST* be immediately before the final output formatters for
% special display devices. If this is done in the output conversion chain
% it must be the last corrective operation before data is fed into the
% formatter plugins. If it is applied to the image processing chains for
% stereo display setups, it must be the absolutely last operation in that
% processing chains before data is fed into output conversion or into the
% stereo compositor.
%
% If we need per view correction for any stereo output mode except
% anaglyph stereo, it needs to happen at end of per view pipeline, so
% things like gamma-correction are applied to final stims, not
% intermediate results. In any other case, there will be only one physical
% output device, so correction is handled best at the end of output
% conversion.
icmshader = [];
icmstring = [];
icmconfig = '';
icmformatting_downstream = 0;

floc = find(mystrcmp(reqs, 'DisplayColorCorrection'));
if ~isempty(floc)
    numColorCorrections = length(floc);

    handlebitspluplus = 0;

    % Bits+ Mono++ or Color++ mode active?
    if ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++Output'))) || ~isempty(find(mystrcmp(reqs, 'EnableBits++Mono++OutputWithOverlay'))) || ~isempty(find(mystrcmp(reqs, 'EnableBits++Color++Output')))
        % Only one 'DisplayColorCorrection' plugin in the whole pipeline?
        if length(floc) == 1
            % Yes: Nothing to do. Full setup for that single plugin has
            % been already done inside our OpenWindow routine. The single
            % plugin has been merged as downstream formatter into the
            % Bits++ output formatting shader via special setup code inside
            % BitsPlusPlus() driver M-File.
            floc = [];
            handlebitspluplus=0;
        else
            % No: No downstream formatting for Bits++ possible whatsoever:
            % Need to do our setup work -- The Bitsplus output formatter
            % just contains a simple neutral clamping shader. However, we
            % need to be careful where to insert our shader(s) if the target
            % is the output conversion chain, as the last slot of that
            % chain is already occupied by the Bits++ shader.
            handlebitspluplus=1;
        end
    end

    % Which channel?
    for x=floc
        [rows cols]= ind2sub(size(reqs), x);
        for row=rows'
            % Extract first parameter - This should be the method of correction:
            colorcorrectionmethod = reqs{row, 3};

            if isempty(colorcorrectionmethod) || ~ischar(colorcorrectionmethod)
                sca;
                error('PsychImaging: Name of color correction method for ''DisplayColorCorrection'' missing or not of string type!');
            end

            % Select method:
            PsychColorCorrection('ChooseColorCorrection', colorcorrectionmethod);

            % Load and build shader objects: icmshader is the compiled
            % color correction shader:
            [icmshader icmstring icmconfig icmoverrideMain] = PsychColorCorrection('GetCompiledShaders', win, 1);

            % Output formatter with built-in ICM capabilities selected? And
            % color correction for final formatting chain insted of
            % per-viewchannel chains?
            if (ptb_outputformatter_icmAware > 0) && (numColorCorrections == 1) && isempty(icmoverrideMain) && (mystrcmp(reqs{row, 1}, 'FinalFormatting') || mystrcmp(reqs{row, 1}, 'AllViews'))
                % Yes. These formatters can use the icm shader internally for
                % higher efficiency if wanted. We can only do that if color
                % correction shall happen in 'AllViews' or 'FinalFormatting', ie.,
                % if this is a monoscopic window or a stereo window where all views
                % display to the same physical output device and therefore the same
                % color correction can be applied to both views.
                %
                % Additionally there must be only 1 color correction stage be present,
                % as multiple stages could can't be done downstream.
                %
                % Additionally there must be no need for a non-standard
                % main() routine for color correction shader.

                % Good. We create the icmshader here according to specs,
                % but then pass it along downstream to the output formatter
                % setup code which will attach it.
                icmformatting_downstream = 1;

            else
                % Downstream color correction not possible due to use of
                % either a per viewchannel correction, or due to use of
                % either no output formatter at all, or not of an icm aware
                % one, or because multi-pass color correction needed, or
                % non-standard main routine needed:
                icmformatting_downstream = 0;

                % Need to build full standalone shader, including main()
                % stub routine and full link and post-link:
                if isempty(icmoverrideMain)
                    % No special override main routine provided. Use our
                    % standard one:
                    shBody = 'uniform sampler2DRect Image; vec4 icmTransformColor(vec4 incolor); void main(void){gl_FragColor = icmTransformColor(texture2DRect(Image, gl_FragCoord.xy));}';
                else
                    % Use provided override routine from
                    % PsychColorCorrection():
                    shBody = icmoverrideMain;
                end

                % shMain is the main() routine which needs to get compiled into
                % a valid shader object:
                shMain = sprintf('\n#extension GL_ARB_texture_rectangle : enable \n\n%s', shBody);
                mainShader = glCreateShader(GL.FRAGMENT_SHADER);
                glShaderSource(mainShader, shMain);
                glCompileShader(mainShader);

                % Link together mainShader and icmshader into a GLSL program
                % object:
                shader = glCreateProgram;
                glAttachShader(shader, icmshader);
                glAttachShader(shader, mainShader);

                % Link the program:
                glLinkProgram(shader);

                % Init the shader: Assign mapping of input image and offsets, gains:
                glUseProgram(shader);
                glUniform1i(glGetUniformLocation(shader, 'Image'), 0);
                glUseProgram(0);
            end

            if ~icmformatting_downstream
                % Ok, shader is our final color correction shader, properly
                % setup. Attach it to proper chain:
                % MK Resolved 26.4.2010: HACK FIXME BUG: 'AllViews' -> Move back to
                % 'FinalFormatting' below, once Screens() pipeline is
                % fixed!!
                if mystrcmp(reqs{row, 1}, 'LeftView') %|| mystrcmp(reqs{row, 1}, 'AllViews')
                    % Need to attach to left view:
                    if leftcount > 0
                        % Need a bufferflip command:
                        Screen('HookFunction', win, 'AppendBuiltin', 'StereoLeftCompositingBlit', 'Builtin:FlipFBOs', '');
                    end
                    Screen('HookFunction', win, 'AppendShader', 'StereoLeftCompositingBlit', icmstring, shader, icmconfig);
                    Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');
                    leftcount = leftcount + 1;
                end

                if mystrcmp(reqs{row, 1}, 'RightView')
                    % Need to attach to right view:
                    if rightcount > 0
                        % Need a bufferflip command:
                        Screen('HookFunction', win, 'AppendBuiltin', 'StereoRightCompositingBlit', 'Builtin:FlipFBOs', '');
                    end
                    Screen('HookFunction', win, 'AppendShader', 'StereoRightCompositingBlit', icmstring, shader, icmconfig);
                    Screen('HookFunction', win, 'Enable', 'StereoRightCompositingBlit');
                    rightcount = rightcount + 1;
                end

                % MK Resolved 26.4.2010: HACK FIXME BUG: 'AllViews' -> Move back to
                % 'FinalFormatting' below, once Screens() pipeline is
                % fixed!!
                if mystrcmp(reqs{row, 1}, 'FinalFormatting') || mystrcmp(reqs{row, 1}, 'AllViews')
                    % Need to attach to final formatting:
                    if ~handlebitspluplus
                        % Standard case:
                        if outputcount > 0
                            % Need a bufferflip command:
                            Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
                        end
                        Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', icmstring, shader, icmconfig);
                    else
                        % Special case: A BitsPlusPlus output formatter has
                        % been attached at the end of queue already. We need
                        % to insert our new slot + some FlipFBO commands just
                        % before the last occupied slot - which is the output formatter slot.
                        % Let's simply count the number of occupied slots and
                        % then insert at that location:
                        insertPos = 0;
                        while(1)
                            if Screen('Hookfunction', win, 'Query', 'FinalOutputFormattingBlit', insertPos)~=-1
                                insertPos = insertPos + 1;
                            else
                                break;
                            end
                        end
                        % insertPos points to first slot after the end of the
                        % chain, ie., where one could append new slots. We want
                        % to insert just at the location of the last slot, so
                        % the last slot gets pushed back one element:
                        insertPos = insertPos - 1;

                        % This insertPos >= 0 check makes sure we also work
                        % in some dummy emulation modes, where no
                        % actual slot is attached:
                        if insertPos >= 0
                            % Need to prepend a bufferflip command in front of
                            % Bitsplusplus:
                            insertSlot = sprintf('InsertAt%iBuiltin', insertPos);
                            Screen('HookFunction', win, insertSlot, 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
                        else
                            % No real output formatter due to emulation
                            % mode, if any.
                            % Force insertPos to 0, so at least
                            % colorcorrection applies:
                            insertPos = 0;
                        end

                        % Then need to prepend our shader in front of that
                        % FlipFBO's:
                        insertSlot = sprintf('InsertAt%iShader', insertPos);
                        Screen('HookFunction', win, insertSlot, 'FinalOutputFormattingBlit', icmstring, shader, icmconfig);

                        % If we're not the first, we need to prepend a
                        % FlipFBO's for ourselves, unless there is already
                        % such a command at the current insertPos:
                        if outputcount > 0
                            % Need to test slot right before us:
                            insertPos = insertPos - 1;

                            % Test what's there at the moment:
                            [dummy testNameString ] = Screen('HookFunction', win, 'Query', 'FinalOutputFormattingBlit', insertPos);
                            if (dummy == - 1) || ~mystrcmp(testNameString, 'Builtin:FlipFBOs')
                                % Need a bufferflip command:
                                insertSlot = sprintf('InsertAt%iBuiltin', insertPos);
                                Screen('HookFunction', win, insertSlot, 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
                            end
                        end
                    end

                    % One more slot occupied by us, so increment
                    % outputcount:
                    outputcount = outputcount + 1;

                    % And enable the chain if it ain't enabled already:
                    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
                end

                % Perform post-link setup of color correction method after
                % shader attached to pipe:
                PsychColorCorrection('ApplyPostGLSLLinkSetup', win, reqs{row, 1});

                % Screen('HookFunction', win, 'Dump', 'FinalOutputFormattingBlit');
            end
        end
    end
end

% Any output formatter to follow which is icmAware, ie., needs to have an
% icmshader as input, either a real one, or a dummy pass-through one?
if ptb_outputformatter_icmAware
    % Yes. To be created output formatter needs an icmshader.
    % Downstream attachment of (already created) icmshader?
    % If so, nothing to do, icmshader and icmstring already setup:
    if ~icmformatting_downstream
        % No. The output formatter is icm aware and needs an icmshader, but
        % none yet created because downstream correction not possible. We
        % need to create a dummy icmshader which just passes through all
        % values uncorrected - This way we make sure that the link
        % operation of the output formatter doesn't fail:
        icmshader = LoadShaderFromFile('ICMPassThroughShader.frag.txt', [], 1);
        icmstring = '';
        icmconfig = '';
    else
        % Nothing to do. Just perform some sanity check here to catch
        % possible future implementation bugs:
        if isempty(icmshader) || isempty(icmstring)
            error('In DisplayColorCorrection setup: Downstream formatting for icmAware output formatter requested, but icmshader and/or icmstring undefined! This is an implementation bug!!!');
        end
    end
end

% --- End of Custom color correction for display wanted ---

% --- User code wants to use unclamped, high precision 0-1 range colors? ---
if ~isempty(find(mystrcmp(reqs, 'NormalizedHighresColorRange')))
    % Use unit color range, without clamping, but in high-precision mode:
    needsUnitUnclampedColorRange = 1;

    % Extract first parameter - This should be the applyAlsoToMakeTexture flag:
    floc = find(mystrcmp(reqs, 'NormalizedHighresColorRange'));
    [rows cols] = ind2sub(size(reqs), floc(1));
    row = rows(1);
    applyAlsoToMakeTexture = reqs{row, 3};
    if ~isempty(applyAlsoToMakeTexture)
        if ~isnumeric(applyAlsoToMakeTexture) || ~ismember(applyAlsoToMakeTexture, [0, 1])
            sca;
            error('In NormalizedHighresColorRange: Invalid applyAlsoToMakeTexture flag specified. Must be 0 or 1.');
        end
    else
        applyAlsoToMakeTexture = [];
    end
end
% --- End of setup for unclamped, high precision 0-1 range colors ---

% --- Native 16 bit per color component RGBA16F floating point framebuffer requested? ---
if ~isempty(find(mystrcmp(reqs, 'EnableNative16BitFloatingPointFramebuffer')))
    % Use unit color range, without clamping, but in high-precision mode:
    needsUnitUnclampedColorRange = 1;
end
% --- End of Native 16 bit per color component RGBA16F floating point framebuffer setup ---

% --- Setup stereo crosstalk reduction ---
floc = find(mystrcmp(reqs, 'StereoCrosstalkReduction'));
if ~isempty(floc)
    if winfo.StereoMode == 0
        sca;
        error('PsychImaging task ''StereoCrosstalkReduction'' requested, but no suitable stereomode active?! Aborted.');
    end

    % Which channel?
    for x=floc
        [rows cols]= ind2sub(size(reqs), x);
        for row=rows'
            crosstalkMethod = reqs{row, 3};
            if isempty(crosstalkMethod) || ~strcmpi(crosstalkMethod, 'SubtractOther')
                sca;
                error('In StereoCrosstalkReduction: Crosstalk reduction method parameter missing or unsupported method requested.');
            end

            crosstalkGain = reqs{row, 4};
            if isempty(crosstalkGain)
                sca;
                error('In StereoCrosstalkReduction: The crosstalk reduction gain must be provided.');
            end

            if isscalar(crosstalkGain)
                % Same gain for all three color channels:
                crosstalkGain = [crosstalkGain crosstalkGain crosstalkGain];
            else
                if numel(crosstalkGain)~=3
                    sca;
                    error('In StereoCrosstalkReduction: provided gain should be a scalar or a 3-element vector.');
                end
            end

            % Background clear color as specified by PsychImaging('Openwindow', ...) call is reference for
            % zero-contrast:
            crosstalkBackGroundClr = clearcolor;
            if isempty(crosstalkBackGroundClr) || ~isnumeric(crosstalkBackGroundClr)
                sca;
                error('In StereoCrosstalkReduction: You did not provide the mandatory background clear color for crosstalk reduction in ''OpenWindow''.');
            end

            if isscalar(crosstalkBackGroundClr)
                % Same background luminance level for all three color channels:
                crosstalkBackGroundClr = [crosstalkBackGroundClr, crosstalkBackGroundClr, crosstalkBackGroundClr];
            else
                if numel(crosstalkBackGroundClr) < 3
                    sca;
                    error('In StereoCrosstalkReduction: Provided background clear color should be a scalar or an at least 3-element RGB(A) vector.');
                end
                crosstalkBackGroundClr = crosstalkBackGroundClr(1:3);
            end

            if min(crosstalkBackGroundClr) <= 0 || max(crosstalkBackGroundClr) >= 1
                sca;
                error('In StereoCrosstalkReduction: Provided background clear color is not in the normalized range > 0 and < 1 as required.');
            end

            % Load and build shader from files StereoCrosstalkReductionShader.vert.txt and/or
            % StereoCrosstalkReductionShader.frag.txt in the shader directory:
            shader = LoadGLSLProgramFromFiles('StereoCrosstalkReductionShader', 1);

            % Init the shader: Assign mapping of images:
            glUseProgram(shader);

            % Image1 will contain the input image for the currently processed target eye:
            glUniform1i(glGetUniformLocation(shader, 'Image1'), 0);

            % Image2 will contain the input image for the to-be-suppressed other eye:
            glUniform1i(glGetUniformLocation(shader, 'Image2'), 1);

            % Just as example. Assign scalar float parameter crosstalkReductionParameter1 to the
            % shader variable 'uniform float crossTalkParam1' for use as a input constant in shader:
            glUniform3fv(glGetUniformLocation(shader, 'crosstalkGain'), 1, crosstalkGain);
            glUniform3fv(glGetUniformLocation(shader, 'backGroundClr'), 1, crosstalkBackGroundClr);

            % Shader setup done:
            glUseProgram(0);

            % Setup specific to left eye output:
            if mystrcmp(reqs{row, 1}, 'LeftView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % Need to attach to left view:
                if leftcount > 0
                    % Need a bufferflip command:
                    Screen('HookFunction', win, 'PrependBuiltin', 'StereoLeftCompositingBlit', 'Builtin:FlipFBOs', '');
                end
                Screen('HookFunction', win, 'PrependShader', 'StereoLeftCompositingBlit', 'StereoCrosstalkReductionShader', shader);
                Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');
                leftcount = leftcount + 1;
            end

            % Setup specific to right eye output:
            if mystrcmp(reqs{row, 1}, 'RightView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % Need to attach to right view:
                if rightcount > 0
                    % Need a bufferflip command:
                    Screen('HookFunction', win, 'PrependBuiltin', 'StereoRightCompositingBlit', 'Builtin:FlipFBOs', '');
                end
                Screen('HookFunction', win, 'PrependShader', 'StereoRightCompositingBlit', 'StereoCrosstalkReductionShader', shader);
                Screen('HookFunction', win, 'Enable', 'StereoRightCompositingBlit');
                rightcount = rightcount + 1;
            end
        end
    end
end
% --- End of setup for stereo crosstalk reduction ---

% --- Custom processing setup for the stereo compositor ---

% --- VR Headset support ---
floc = find(mystrcmp(reqs, 'UseVRHMD'));
if ~isempty(floc)
    [row col] = ind2sub(size(reqs), floc);

    % Extract first parameter - This should be the handle of the VR device:
    hmd = reqs{row, 3};

    % Verify it is already open:
    if ~hmd.driver('IsOpen', hmd)
        sca;
        error('In UseVRHMD: Invalid HMD handle specified. No such device opened.');
    end

    % Ok, perform setup after onscreen window is open, e.g., setting up the special
    % shaders for the stereo compositor:
    if ~hmd.driver('PerformPostWindowOpenSetup', hmd, win, clearcolor)
        sca;
        error('In UseVRHMD: Failed to setup image post-processing for the VR HMD.');
    end

    % Ready to rock the HMD!
end
% --- End of VR Headset support code. ---


% --- FROM HERE ON ONLY OUTPUT FORMATTERS, NOTHING ELSE!!! --- %


% --- Final output formatter for generic LUT based luminance framebuffer requested? ---
% --- OR Final output formatter for Pseudo-Gray processing requested? ---
floc = find(mystrcmp(reqs, 'EnableGenericHighPrecisionLuminanceOutput'));
if isempty(floc)
    floc = find(mystrcmp(reqs, 'EnablePseudoGrayOutput'));
end
if ~isempty(floc)
    [row col]= ind2sub(size(reqs), floc);

    % Get native depth in bits per color (bpc) of active framebuffer:
    nativeBPC = Screen('Pixelsize', win) / 3;

    if mystrcmp(reqs{row, 2}, 'EnablePseudoGrayOutput')
        % PseudoGray mode: We create the lut ourselves via helper function:
        lut = CreatePseudoGrayLUT(nativeBPC);

        % For proper pseudo-gray output the gfx gamma-tables must not be
        % touched by us!
        needsIdentityCLUT = 0;
    else
        % Generic "attenuator" driver mode:

        % Attenuator-style devices need an identity clut:
        needsIdentityCLUT = 1;

        % Extract first parameter - This should be the lookup table 'lut' to use:
        lut = reqs{row, 3};
    end

    if isempty(lut) || ~isnumeric(lut)
        sca;
        error('PsychImaging: Mandatory lookup table parameter lut for ''%s'' missing or not of numeric type!', reqs{floc});
    end

    % Load output formatting shader for GenericHighPrecisionLuminanceOutput:
    % 'icmshader' is a handle to a compiled fragment shader, provided by
    % upstream, that implements the display color correction function:
    pgshader = LoadGLSLProgramFromFiles('GenericLuminanceToRGBA8_FormattingShader', 1, icmshader);

    % Init the shader: Assign mapping texture units etc.:
    glUseProgram(pgshader);
    glUniform1i(glGetUniformLocation(pgshader, 'Image'), 0);
    glUniform1i(glGetUniformLocation(pgshader, 'LUT'),  1);
    glUniform1f(glGetUniformLocation(pgshader, 'MaxIndex'),  size(lut, 2)-1);
    glUseProgram(0);

    % Use helper routine to build a proper RGBA Lookup texture for
    % conversion of HDR luminance pixels to output framebuffer pixels:
    if isa(lut, 'uint8')
        % uint8 classic lut: Create RGBA8 lookup texture:
        fprintf('PsychImaging-%s: Creating LUT suitable for precision boosting of a 8 bpc native framebuffer and display or DAC.\n', reqs{floc});
        pglutid = PsychHelperCreateGenericLuminanceToRGBA8LUT(lut);
    elseif isa(lut, 'uint16')
        % uint16 input lut: Store as texture of suitable depths (8, 10 or 16 bpc):
        fprintf('PsychImaging-%s: Creating %i slot LUT suitable for precision boosting of a %i bpc native framebuffer and display or DAC.\n', ...
                reqs{floc}, size(lut, 2), nativeBPC);
        pglutid = PsychHelperCreateGenericLuminanceToRGBA16MaxLUT(lut, nativeBPC, win);
    else
        sca;
        error('PsychImaging: Mandatory lookup table parameter lut for ''%s'' not of uint8 or uint16 type, as required!', reqs{floc});
    end

    if outputcount > 0
        % Need a bufferflip command:
        Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
    end
    pgconfig = sprintf('TEXTURERECT2D(1)=%i %s', pglutid, icmconfig);
    pgidstring = sprintf('Generic high precision luminance output formatting shader: %s', icmstring);
    Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', pgidstring, pgshader, pgconfig);
    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
    outputcount = outputcount + 1;

    % Use unit color range, without clamping, but in high-precision mode:
    needsUnitUnclampedColorRange = 1;
end
% --- End of output formatter for generic LUT based luminance framebuffer ---

% --- Final output formatter for VideoSwitcher attenuator device requested? ---
VideoSwitcherTriggerflag = 0;

floc = find(mystrcmp(reqs, 'EnableVideoSwitcherSimpleLuminanceOutput'));
if isempty(floc)
    floc = find(mystrcmp(reqs, 'EnableVideoSwitcherCalibratedLuminanceOutput'));
end
if ~isempty(floc)
    [row col]= ind2sub(size(reqs), floc);

    if mystrcmp(reqs{row, 2}, 'EnableVideoSwitcherSimpleLuminanceOutput')
        simpleVideoSwitcher = 1;
    else
        simpleVideoSwitcher = 0;
    end

    % Extract optional first parameter - This should be the 'btrr' ratio to use:
    btrr = reqs{row, 3};

    if isempty(btrr)
        % btrr empty: Get it from config file:
        btrr = PsychVideoSwitcher('GetDefaultConfig', win);
    else
        if ~isnumeric(btrr) || ~isscalar(btrr)
            sca;
            error('PsychImaging: Optional "btrr" parameter for VideoSwitcher output not of numeric scalar type!');
        end

        if btrr < 0
            sca;
            error('PsychImaging: Optional "btrr" parameter for VideoSwitcher output is negative -- Impossible!');
        end
    end

    if simpleVideoSwitcher
        % Extract optional 2nd parameter - This should be the 'trigger' flag:
        VideoSwitcherTriggerflag = reqs{row, 4};
    else
        % Extract optional 3rd parameter - This should be the 'trigger' flag:
        VideoSwitcherTriggerflag = reqs{row, 5};
    end

    if isempty(VideoSwitcherTriggerflag)
        % triggerflag empty: Default to off:
        VideoSwitcherTriggerflag = 0;
    else
        if ~isnumeric(VideoSwitcherTriggerflag) || ~isscalar(VideoSwitcherTriggerflag)
            sca;
            error('PsychImaging: Optional "trigger" parameter for VideoSwitcher output not of numeric scalar type!');
        end

        if VideoSwitcherTriggerflag > 0
            VideoSwitcherTriggerflag = 1;
        else
            VideoSwitcherTriggerflag = 0;
        end
    end

    if simpleVideoSwitcher
        % Load output formatting shader for simple VideoSwitcher output:
        % 'icmshader' is a handle to a compiled fragment shader, provided by
        % upstream, that implements the display color correction function:
        pgshader = LoadGLSLProgramFromFiles('VideoSwitcherSimpleLuminanceToRB8_FormattingShader', 1, icmshader);

        % Init the shader: Assign mapping texture units etc.:
        glUseProgram(pgshader);
        glUniform1i(glGetUniformLocation(pgshader, 'Image'), 0);
        glUniform1f(glGetUniformLocation(pgshader, 'btrr'),  btrr);
        glUniform1f(glGetUniformLocation(pgshader, 'btrrPlusOne'),  btrr + 1);
        glUniform1f(glGetUniformLocation(pgshader, 'btrrFractionTerm'), ((btrr + 1) / btrr));
        glUseProgram(0);

        pgidstring = sprintf('VideoSwitcher simple high precision luminance output formatting shader: %s', icmstring);
        pgconfig = icmconfig;
    else
        % LUT calibrated VideoSwitcher setup:

        % Extract optional 2nd parameter - This should be the 'lut':
        lut = reqs{row, 4};

        if isempty(lut)
            % lut empty: Get it from config file:
            [dummy, lut] = PsychVideoSwitcher('GetDefaultConfig', win);
        else
            if ~isa(lut, 'double') || ~isvector(lut) || length(lut)~=257
                sca;
                error('PsychImaging: Lookup table parameter lut for VideoSwitcher output invalid: Must be a vector of double values with 257 elements!');
            end
        end

        % Load output formatting shader for lut calibrated VideoSwitcher output:
        % 'icmshader' is a handle to a compiled fragment shader, provided by
        % upstream, that implements the display color correction function:
        pgshader = LoadGLSLProgramFromFiles('VideoSwitcherCalibratedLuminanceToRB8_FormattingShader', 1, icmshader);

        % Init the shader: Assign mapping texture units etc.:
        glUseProgram(pgshader);
        glUniform1i(glGetUniformLocation(pgshader, 'Image'), 0);
        glUniform1i(glGetUniformLocation(pgshader, 'LUT'), 1);
        glUniform1f(glGetUniformLocation(pgshader, 'btrr'),  btrr);

        % Disable luminance hint color by default by setting it to the
        % luminance key -1, which shouldn't ever match in a regular
        % stimulus:
        glUniform3f(glGetUniformLocation(pgshader, 'BackgroundPixel'), 0, -1, 0);
        glUseProgram(0);

        % Convert 'lut' into lookup table texture:
        pglutid = PsychVideoSwitcher('GetLUTTexture', win, lut, btrr, pgshader);

        pgidstring = sprintf('VideoSwitcher calibrated high precision luminance output formatting shader: %s', icmstring);
        pgconfig = sprintf('TEXTURERECT2D(1)=%i %s', pglutid, icmconfig);
    end

    if outputcount > 0
        % Need a bufferflip command:
        Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
    end
    Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', pgidstring, pgshader, pgconfig);
    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
    outputcount = outputcount + 1;

    % VideoSwitcher devices need an identity clut:
    needsIdentityCLUT = 1;

    % Use unit color range, without clamping, but in high-precision mode:
    needsUnitUnclampedColorRange = 1;
end

% Setup of trigger for VideoSwitcher device needed?
if VideoSwitcherTriggerflag > 0
    % Yes. Attach a proper slot to the chain: The slot calls back into the
    % VideoSwitcher.m M-File, with the window handle as argument.
    pgconfig = sprintf('PsychVideoSwitcher(%i);', win);
    Screen('HookFunction', win, 'AppendMFunction', 'FinalOutputFormattingBlit', 'VideoSwitcher trigger control callback.', pgconfig);
end

% --- End of output formatters for VideoSwitcher attenuator device ---


% --- Final output formatter for native 10 bpc ARGB2101010 or 11 bpc RGB11-11-10 framebuffer requested?
enableNative11BpcRequested = 0;
floc = find(mystrcmp(reqs, 'EnableNative10BitFramebuffer'));
if isempty(floc)
    enableNative11BpcRequested = 1;
    floc = find(mystrcmp(reqs, 'EnableNative11BitFramebuffer'));
end

if isempty(floc)
    enableNative11BpcRequested = 0;
    floc = find(mystrcmp(reqs, 'EnableNative16BitFramebuffer'));
end

if ~isempty(floc)
    [row col]= ind2sub(size(reqs), floc);

    % Our special shader-based 10 bpc output formatter is only applicable on Linux
    % with AMD Radeon hardware, or with FireGL/FirePro with override mode bit set.
    % Our 11 bpc shader-based output formatter is only for Linux + AMD GCN-1.1+.
    % specialFlags setting 1024 signals that our own low-level 10/11 bit framebuffer
    % hack on AMD hardware is active, so we also need our own GLSL output formatters.
    % Otherwise setup was (hopefully) done by the regular graphics drivers and we don't
    % need this GLSL output formatter, as system OpenGL takes care of it.
    % If the Vulkan display backend is requested then this also does not apply as
    % the Vulkan/WSI backend must do whatever neccessary to provide the requested
    % fixed point unorm precision - or we simply fail if it can't:
    if bitand(winfo.SpecialFlags, 1024) && ~useVulkan
        % AMD/ATI gpu on Linux with our 10/11 bit hack. Use our reformatters:

        % Load output formatting shader for multi-LUT based 10 bpc or 11 bpc formatting:
        pgshader = LoadGLSLProgramFromFiles('RGBMultiLUTLookupCombine_FormattingShader', 1);

        % Init the shader:
        glUseProgram(pgshader);

        % Assign mapping of input image:
        glUniform1i(glGetUniformLocation(pgshader, 'Image'), 0);

        % CLUT based mapping:
        glUniform1i(glGetUniformLocation(pgshader, 'CLUT'), 1);
        glUniform1f(glGetUniformLocation(pgshader, 'Prescale'), bitshift(1024, enableNative11BpcRequested) - 1);
        glUseProgram(0);

        if enableNative11BpcRequested
            % Use helper routine to build a proper RGBA Lookup texture for
            % conversion of RGB pixels to ARGB0-11-11-10 pixels.
            if winfo.GPUMinorType >= 80
              % DCE-8 to DCE-12 are the only AMD gpu's which are supported and provide
              % 12 bpc output to make this actually meaningful:
              pglutid = PsychHelperCreateRGB111110RemapCLUTAMDDCE8;
            else
              sca;
              error('PsychImaging: Native11BitFramebuffer is not supported on this too old or too new AMD gpu.');
            end
            pgshadername = 'Native RGB111110 framebuffer output formatting shader';
            pgconfig = sprintf('TEXTURERECT2D(1)=%i', pglutid);
        else
            % Use helper routine to build a proper RGBA Lookup texture for
            % conversion of RGBA pixels to ARGB2101010 pixels:
            pglutid = PsychHelperCreateARGB2101010RemapCLUT;
            pgshadername = 'Native ARGB2101010 framebuffer output formatting shader';
            pgconfig = sprintf('TEXTURERECT2D(1)=%i', pglutid);
        end

        if outputcount > 0
            % Need a bufferflip command:
            Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
        end

        Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', pgshadername, pgshader, pgconfig);
        Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
        outputcount = outputcount + 1;

        % AMD framebuffer devices - Identity CLUT not needed, as internal clut is bypassed anyway,
        % but we do it nonetheless, so we can decide about dithering setup and get things like
        % degamma and other colorspace conversions disabled / bypassed:
        needsIdentityCLUT = 1;
    else
        % Everything else: Windows OS or macOS, or AMD FireGL/FirePro without override,
        % or AMD with amdgpu DisplayCore, or any NVidia or Intel GPU. Also on Vulkan.
        % Do not request an identity lut. Modern Intel, NVidia and AMD gpu's have
        % hw LUT's with an output width of potentially more than 10 bpc, so we
        % can potentially benefit from a higher precision gamma correction via
        % hw lut. E.g., Intel Icelake has up to 16 bit output precision lut's,
        % NVidia up to 14 bit, AMD greater than 10 bit - typically 12 bit.
        % Going through our identity lut setup code could even load a "identity lut"
        % that truncates output precision to 8 bit, e.g., on Linux + Intel gpu's,
        % as our LoadIdentityClut() function is optimized/targeted at 8 bpc passthrough.
        needsIdentityCLUT = 0;
    end

    % Extract optional first parameter - This should be the 'disableDithering' flag:
    disableDithering = reqs{row, 3};

    if isempty(disableDithering)
        % Control of output dithering on digital >= 10 bit panels should be left to
        % the OS + graphics driver by default. For example, the OS can configure
        % the encoders for 10 bpc no-dithering if it detects a truly 10 bpc capable display,
        % based on EDID information. DisplayPort and HDMI provides infos about >= 10 bpc
        % capabilities in their EDID info. If the OS detects a <= 8 bpc digital panel, it
        % can dither so we get pseudo-10bpc, similar to a bit stealing approach or other
        % perceptual high bit depths tricks. The same is true for driving 10 bpc panels with
        % 12 bpc or more precision, using dithering on the 10 bpc signal:
        disableDithering = 0;
    else
        % User provided disableDithering flag. Valid?
        if ~ismember(disableDithering, [0, 1])
            sca;
            error('Optional disableDithering flag with invalid value provided! Valid is 0 or 1!');
        end

        % Yes, use it.
    end

    % Use unit color range, without clamping, but in high-precision mode:
    needsUnitUnclampedColorRange = 1;
end
% --- End of output formatter for native 10 bpc ARGB2101010 framebuffer ---

% --- Experimental output formatter for Dual-Pipeline HDR display ---
floc = find(mystrcmp(reqs, 'EnableDualPipeHDROutput'));
if ~isempty(floc)
    [row col]= ind2sub(size(reqs), floc);

    % outputcount should be zero, i.e., the unified output formatting chain
    % should be disabled, as we use separate per channel chains:
    if outputcount > 0
        fprintf('PsychImaging: WARNING! In setup for task "EnableDualPipeHDROutput": Unified output formatting chain was active (count = %i)!\n', outputcount);
        fprintf('PsychImaging: WARNING! This conflicts with need for separate output formatting chains! Overriding: Unified chain disabled!\n');
        fprintf('PsychImaging: WARNING! Check your output stimulus carefully for artifacts!\n');

        % Disable unified output formatting chain and hope for the best:
        Screen('HookFunction', win, 'Disable', 'FinalOutputFormattingBlit');
        % Screen('HookFunction', win, 'Disable', 'RightFinalizerBlitChain');
    end

    % Setup shader for pipe 0:
    pipe0shader = LoadGLSLProgramFromFiles('DualPipeHDRPipe0_FormattingShader', 1, icmshader);

    if outputcount0 > 0
        % Need a bufferflip command:
        Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit0', 'Builtin:FlipFBOs', '');
    end

    Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit0', 'HDRPipe0 - Output Formatter', pipe0shader, '');
    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit0');
    outputcount0 = outputcount0 + 1;

    % Setup shader for pipe 1:
    pipe1shader = LoadGLSLProgramFromFiles('DualPipeHDRPipe1_FormattingShader', 1, icmshader);

    if outputcount1 > 0
        % Need a bufferflip command:
        Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit1', 'Builtin:FlipFBOs', '');
    end

    Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit1', 'HDRPipe1 - Output Formatter', pipe1shader, '');
    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit1');
    outputcount1 = outputcount1 + 1;

    % Device need an identity clut in the GPU gamma tables:
    needsIdentityCLUT = 1;

    % Use unit color range, without clamping, but in high-precision mode:
    needsUnitUnclampedColorRange = 1;
end
% --- End of experimental output formatter for Dual-Pipeline HDR display ---

% --- Output formatter for Eizo RadiForce style 8 bit luminance subpixel drive ---
floc = find(mystrcmp(reqs, 'UseSubpixelDrive'));
if ~isempty(floc)
    [row col]= ind2sub(size(reqs), floc);

    if outputcount > 0
        % Need a bufferflip command:
        Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
    end

    pgshader = LoadGLSLProgramFromFiles('SubpixelDrive_FormattingShader', 1);

    % Init the shader: Assign mapping texture units etc.:
    glUseProgram(pgshader);
    glUniform1i(glGetUniformLocation(pgshader, 'Image'), 0);
    glUseProgram(0);

    Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', 'Subpixel drive - Output Formatter', pgshader, '');
    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
    outputcount = outputcount + 1;
end
% --- End of output formatter for Eizo RadiForce style 8 bit luminance subpixel drive ---

% --- GPU based mirroring of left half of onscreen window to right half requested? ---
if ~isempty(find(mystrcmp(reqs, 'MirrorDisplayToSingleSplitWindow')))
    floc = find(mystrcmp(reqs, 'MirrorDisplayToSingleSplitWindow'));
    [rows ~] = ind2sub(size(reqs), floc);

    % Simply set up the left finalizer chain with a glCopyPixels command
    % that copies the left half of the system backbuffer to the right half
    % of the system backbuffer. Query the real backbuffer width x height,
    % but use half the width as source region and destination region
    % offset, as the right half of the backbuffer shall be a copy of the
    % left half:
    [w, h] = Screen('WindowSize', win, 1);
    w = w / 2;

    % Imaging pipeline fully enabled? Specific offsets used for blitter
    % commands depend on this:
    if bitand(winfo.ImagingMode, kPsychNeedFastBackingStore) > 0
        % Yes: Use proper offsets for active imaging pipeline:
        ow = w;
    else
        % No: Need different x-offset for glRasterPos2f, because the good
        % old fixed function pipeline uses different viewport / projection
        % matrix etc.:
        ow = w / 2;
    end

    % Setup scaling (and location tbd. in the future) of mirror target region,
    % width default of full scaling to half width right half of window for mirror image:
    sfx = 1;
    sfy = 1;

    % Get optional region spec for mirror target region:
    mirrordst = reqs{rows, 4};
    if ~isempty(mirrordst)
        % Parameter specified as non-default: Valid format?

        % Currently we only accept numeric row-vectors:
        if isnumeric(mirrordst) && isreal(mirrordst) && isrow(mirrordst)
            % Of length two elements for [width, height] of the mirror region:
            if length(mirrordst) == 2
                % width and height must be at least 1 pixel:
                if any(mirrordst < 1)
                    sca;
                    error('PsychImaging: MirrorDisplayToSingleSplitWindow: Invalid mirrorDestination parameter [width, height] vector, width or height < 1 pixel!');
                end

                % Extract width x height:
                dstw = mirrordst(1);
                dsth = mirrordst(2);
            else
                sca;
                error('PsychImaging: MirrorDisplayToSingleSplitWindow: Invalid mirrorDestination parameter provided: Not a two element [width, height] vector!');
            end
        else
            sca;
            error('PsychImaging: MirrorDisplayToSingleSplitWindow: Invalid mirrorDestination parameter provided: Not a numeric row vector!');
        end

        % Compute scaling factors from stimulus image to mirror region:
        sfx = dstw / w;
        sfy = dsth / h;
    end

    % Attach blit command sequence to finalizer chain:
    blitstring = sprintf('glPixelZoom(%f, %f); glRasterPos2f(%f, %f); glCopyPixels(0, 0, %f, %f, 6144); ', sfx, sfy, ow, h * sfy, w, h);

    % Overlay for mirror window requested?
    if reqs{rows, 3} == 1
        % Need fast FBO backed offscreen windows for overlay support:
        if ~bitand(winfo.ImagingMode, kPsychNeedFastOffscreenWindows + kPsychNeedFastBackingStore)
            % No-go:
            sca;
            error('PsychImaging: MirrorDisplayToSingleSplitWindow: Experimenter overlay can not be used if fast offscreen windows are disabled!');
        end

        % Create Offscreen window for the overlay. It has the same size as
        % the onscreen window, and the same pixeldepth, but a completely black
        % background with alpha value zero -- fully transparent by default.
        % The specialflags 32 setting protects the overlay offscreen window
        % from accidental batch-deletion by usercode calling Screen('Close'):
        overlaywin = Screen('OpenOffscreenWindow', win, [0 0 0 0], [0 0 w h], [], 32);
        ptb_MirrorOverlayWindows(win) = overlaywin;

        % 'GetWindowInfo' binds the FBO of our Offscreen window, so we can get its fbo
        % id in overlayfbo:
        Screen('GetWindowInfo', overlaywin);
        overlayfbo = glGetIntegerv(GL.READ_FRAMEBUFFER_BINDING);

        % Build blitter command string: We use a glCopyPixels() from the overlaywin
        % overlayfbo as read framebuffer, to the active draw framebuffer. glCopyPixels
        % allows alpha testing or blending, whereas glBlitFramebuffer would not. Then
        % enable alpha testing around the blit, so all overlay pixels with non-zero
        % alpha overwrite the mirror image pixels:
        blitstring = [blitstring sprintf(['glAlphaFunc(516, 0.0); glEnable(3008); rfbo = glGetIntegerv(36010); ' ...
                     'glBindFramebufferEXT(36008, %i); glRasterPos2f(%f, %f); glCopyPixels(0, 0, %f, %f, 6144); ' ...
                     'glBindFramebufferEXT(36008, rfbo); glDisable(3008); glPixelZoom(1, 1);'], overlayfbo, ow, h * sfy, w, h)];
    else
        blitstring = [blitstring 'glPixelZoom(1, 1);'];
    end

    Screen('Hookfunction', win, 'AppendMFunction', 'LeftFinalizerBlitChain', 'MirrorDisplayToSingleSplitWindow', blitstring);
    Screen('HookFunction', win, 'Enable', 'LeftFinalizerBlitChain');
end
% --- End of GPU based mirroring of left half of onscreen window to right half requested? ---

% --- GPU based mirroring from master onscreen window to slave mirror window requested? ---
if ~isempty(find(mystrcmp(reqs, 'MirrorDisplayTo2ndOutputHead')))
    floc = find(mystrcmp(reqs, 'MirrorDisplayTo2ndOutputHead'));
    [rows ~] = ind2sub(size(reqs), floc);

    % Enable dedicated blit chain for mirror windows:
    Screen('HookFunction', slavewin, 'Enable', 'MirrorWindowBlit');
    Screen('HookFunction', slavewin, 'Reset', 'MirrorWindowBlit');

    % If the slavewin mirror window has a different size than the master window that
    % is mirrored, we need to setup special override viewport and projection matrix:
    [w, h] = Screen('WindowSize', win, 1);
    [sw, sh] = Screen('WindowSize', slavewin, 1);
    if w ~= sw || h ~= sh
        sw = floor(sw);
        sh = floor(sh);
        Screen('Hookfunction', slavewin, 'AppendMFunction', 'MirrorWindowBlit', 'Setup special viewports', ...
                sprintf('glViewport(0, 0, %i, %i); glMatrixMode(5889); glLoadIdentity(); gluOrtho2D(0, %i, %i, 0); glMatrixMode(5888);', sw, sh, sw, sh));
    end

    % Add blit command for master window -> slave window to copy the mirror image:
    Screen('HookFunction', slavewin, 'AppendBuiltin', 'MirrorWindowBlit', 'Builtin:IdentityBlit', ...
            sprintf('OvrSize:%i:%i', sw, sh));

    % Overlay for mirror window requested?
    if reqs{rows, 6} == 1
        % Get texture target type of the finalizedFBO's in use, as currently the
        % imaging pipelines 'Builtin:IdentityBlit' builtin function will always
        % generate texture coordinates for blitting which are suitable only for
        % that texture target. Iow. Our overlay window overlaywin must use a
        % FBO backing texture of the same texture target, or things will go sideways:
        [~, ~, ttarget] = Screen('HookFunction', win, 'GetDisplayBufferTextures');

        % Create Offscreen window for the overlay. It has the same size as
        % the onscreen window, and the same pixeldepth, but a completely black
        % background with alpha value zero -- fully transparent by default.
        % The specialflags 32 setting protects the overlay offscreen window
        % from accidental batch-deletion by usercode calling Screen('Close'):
        if ttarget == GL.TEXTURE_2D
            % Need a GL.TEXTURE_2D backing texture:
            overlaywin = Screen('OpenOffscreenWindow', win, [0 0 0 0], [0 0 w h], [], 32 + 1);
            samplerstring = 'TEXTURE2D';
        else
            % Need a GL.TEXTURE_RECTANGLE backing texture:
            overlaywin = Screen('OpenOffscreenWindow', win, [0 0 0 0], [0 0 w h], [], 32);
            samplerstring = 'TEXTURERECT2D';
        end

        ptb_MirrorOverlayWindows(win) = overlaywin;

        % Retrieve low-level OpenGL texture handle to the window:
        overlaytex = Screen('GetOpenGLTexture', win, overlaywin);

        % Append blitter command for a one-to-one blit of the overlay window
        % texture to the target buffer. We need to enable alpha testing, so the
        % overlay only occludes the mirrored image where the overlay has a non-zero alpha:
        Screen('Hookfunction', slavewin, 'AppendMFunction', 'MirrorWindowBlit', 'Setup Alphatest for Overlay', 'glAlphaFunc(516, 0.0); glEnable(3008);');
        Screen('Hookfunction', slavewin, 'AppendBuiltin', 'MirrorWindowBlit', 'Builtin:IdentityBlit', ...
               sprintf('%s(0)=%i:OvrSize:%i:%i', samplerstring, overlaytex, sw, sh));
        Screen('Hookfunction', slavewin, 'AppendMFunction', 'MirrorWindowBlit', 'Teardown Alphatest for Overlay', 'glDisable(3008);');
    end
end
% --- End of GPU based mirroring from master onscreen window to slave mirror window requested? ---

% --- Datapixx in use? ---
if ~isempty(find(mystrcmp(reqs, 'UseDataPixx')))
    % Yes: Need to call into high level DataPixx driver for final setup:
    PsychDataPixx('PerformPostWindowOpenSetup', win);
end
% --- End of Datapixx in use? ---

% --- Bits# in use? ---
if ~isempty(find(mystrcmp(reqs, 'UseBits#')))
    % Yes: Need to call into high level BitsPlusPlus driver for final setup:
    BitsPlusPlus('PerformPostWindowOpenSetup', win);
end
% --- End of Bits# in use? ---

if useHDR && needsIdentityCLUT
    % We want identity lut's in HDR, but at maximum lut precision, so output
    % does not get truncated to 8 bpc. Therefore we can't use LoadIdentityClut()
    %
    % There we have identity lut setup code in PsychVulkan() to handle this, and
    % therefore useHDR  && needsIdentityCLUT should never happen:
    warning('needsIdentityCLUT set in HDR mode! Bug!?!');
end

% Do we need identity gamma tables / CLUT's loaded into the graphics card?
if needsIdentityCLUT && ~useHDR
    % Yes. Use our generic routine which is adaptive to the quirks of
    % specific gfx-cards:
    LoadIdentityClut(win, [], [], disableDithering);

    % Is there a slave window associated for some dual-window output mode,
    % HDR mode or stereo mode?
    if ~isempty(slavewin)
        % Yes: Apply identity LUT setup there as well:
        LoadIdentityClut(slavewin, [], [], disableDithering);
    end
end

% Special Vulkan display backend in use?
floc = find(mystrcmp(reqs, 'UseVulkanDisplay'));
if ~isempty(floc)
    [rows cols] = ind2sub(size(reqs), floc(1));
    row = rows(1);

    % Extract first parameter - This would be the optional video output name:
    outputName = reqs{row, 3};

    if ~isempty(outputName)
        isFullscreen = 1;
    else
        isFullscreen = 0;
    end

    % 10 bpc linear unorm output framebuffer requested?
    if ~isempty(find(mystrcmp(reqs, 'EnableNative10BitFramebuffer')))
        % Request code 1 - RGB10A2 / BGR10A2 format on Vulkan side:
        % This also overrides potential automatic choice  by 'EnableHDR' setup.
        vulkanColorPrecision = 1;
    end

    % 11 bpc linear unorm output framebuffer requested? Or 16 bpc float framebuffer?
    if ~isempty(find(mystrcmp(reqs, 'EnableNative11BitFramebuffer'))) || ...
       ~isempty(find(mystrcmp(reqs, 'EnableNative16BitFloatingPointFramebuffer')))
        % Request code 2 - RGBA16F format on Vulkan side:
        % This also overrides potential automatic choice  by 'EnableHDR' setup.
        % fp16 request obviously matches exactly RGBA16F / VK_FORMAT_R16G16B16A16_SFLOAT.
        %
        % The same format works for RGB11 unorm, because it so happens that inside the 0.0 - 1.0
        % unsigned normalized (unorm) color value range of RGB11 / EnableNative11BitFramebuffer,
        % a 16 bit half-float floating point value provides the equivalent of at least 11 bpc
        % linear precision:
        vulkanColorPrecision = 2;
    end

    % 16 bpc linear unorm output framebuffer requested?
    if ~isempty(find(mystrcmp(reqs, 'EnableNative16BitFramebuffer')))
        % Request code 3 - RGBA16 unorm fixed point format on Vulkan side:
        % This also overrides potential automatic choice  by 'EnableHDR' setup.
        vulkanColorPrecision = 3;
    end

    % Default flags:
    flags = 0;

    % Default: Auto-select optimal Vulkan driver+gpu combo:
    gpuIndex = 0;

    % Special static HDR stereo hack for Linux/X11 + OpenGL requested?
    floc = find(mystrcmp(reqs, 'UseStaticHDRHack'));
    if isempty(floc)
        % No, standard Vulkan display mode: Perform Vulkan onscreen window creation and setup of Vulkan and OpenGL interop on our side:
        vwin = PsychVulkan('PerformPostWindowOpenSetup', win, Screen('GlobalRect', win), isFullscreen, outputName, vulkanHDRMode, vulkanColorPrecision, vulkanColorSpace, vulkanColorFormat, gpuIndex, flags);
    else
        % Yes: Perform one-time dual-display HDR setup via special hack in PsychHDR():

        % Extract static HDR metadata to use, if any:
        [rows cols] = ind2sub(size(reqs), floc(1));
        row = rows(1);
        staticHDRMetadata = reqs{row, 3};

        % Perform one-time setup via Vulkan, with static assignment of staticHDRMetadata
        % for the duration of the whole session:
        PsychHDR('ExecuteStaticHDRHack', win, 1, vulkanHDRMode, gpuIndex, flags, staticHDRMetadata);
    end
end

if useHDR
    % Perform post OpenWindow setup for HDR stuff, e.g., setting up appropriate
    % HDR post-processing shaders (Optional tone-mapping, gamut-remapping,
    % color-space conversion and clamping. Mandatory OETF mapping, e.g., PQ):
    [hdrShader, hdrShaderString] = PsychHDR('PerformPostWindowOpenSetup', win, hdrArguments, icmshader, icmstring);

    % If any hdrShader is needed, append it to the final output formatter chain:
    if ~isempty(hdrShader)
        if outputcount > 0
            % Need a bufferflip command:
            Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
        end

        Screen('HookFunction', win, 'AppendShader', 'FinalOutputFormattingBlit', hdrShaderString, hdrShader, icmconfig);
        Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
        outputcount = outputcount + 1;
    end

    % Make sure that regular double() input images, even if only in the standard
    % unorm range for creation of 8 bpc RGBA8 textures, do get their input values
    % used "as is", ie. valid values are 0.0 - 1.0:
    applyAlsoToMakeTexture = 1;

    % Disable input value scaling and clamping, so color values can be passed through as
    % arbitrary floating point values, for encoding in Nits, 80Nits, negative values,
    % chromaticities etc.:
    needsUnitUnclampedColorRange = 1;
end

% --- END OF ALL OUTPUT FORMATTERS AND SPECIAL OUTPUT DEVICE/BACKEND SETUP ---

% --- This must be after setup of all output formatter shaders! ---
% Downstream icm color correction shader linked into an icmAware output
% formatter. We must perform post-link setup for it:
if ptb_outputformatter_icmAware && icmformatting_downstream
    % Perform post-link setup of color correction method after
    % shader attached to pipe. We know it is the
    % 'FinalOutputFormattingBlit' chain, as only in that case, downstream
    % formatting is performed at all.
    PsychColorCorrection('ApplyPostGLSLLinkSetup', win, 'FinalFormatting');
end

% --- Restriction of processing area ROI requested? ---

% This should be at the end of hookchain setup, so we can reliably prepend the
% command to each chain to guarantee that restriction applies to all
% processing:
floc = find(mystrcmp(reqs, 'RestrictProcessing'));
if ~isempty(floc)
    % Which channel?
    for x=floc
        [rows cols]= ind2sub(size(reqs), x);
        for row=rows'
            % Extract scissor rectangle:
            scissorrect = reqs{row, 3};
            if size(scissorrect,1)~=1 || size(scissorrect,2)~=4
                error('Task "RestrictProcessing" in channel %s expects a 1-by-4 ROI rectangle to define the ROI, e.g, [left top right bottom]!', reqs{row,1});
            end

            ox = scissorrect(RectLeft);
            oy = winheight - scissorrect(RectBottom);

            w  = RectWidth(scissorrect);
            h  = RectHeight(scissorrect);

            if mystrcmp(reqs{row, 1}, 'LeftView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % Need to restrict left view processing:
                Screen('HookFunction', win, 'PrependBuiltin', 'StereoLeftCompositingBlit', 'Builtin:RestrictToScissorROI', sprintf('%i:%i:%i:%i', ox, oy, w, h));
            end

            if mystrcmp(reqs{row, 1}, 'RightView') || mystrcmp(reqs{row, 1}, 'AllViews')
                % Need to restrict right view processing:
                Screen('HookFunction', win, 'PrependBuiltin', 'StereoRightCompositingBlit', 'Builtin:RestrictToScissorROI', sprintf('%i:%i:%i:%i', ox, oy, w, h));
            end

            if (mystrcmp(reqs{row, 1}, 'AllViews') || mystrcmp(reqs{row, 1}, 'Compositor')) && ismember(winfo.StereoMode, [6,7,8,9])
                % Needed to restrict both views processing and a
                % compositing mode is active. If both views are restricted
                % in their output area then it makes sense to restrict the
                % compositor to the same area. We also restrict the
                % compositor if that was requested.
                oy = RectHeight(Screen('Rect', win, 1)) - scissorrect(RectBottom);
                Screen('HookFunction', win, 'PrependBuiltin', 'StereoCompositingBlit', 'Builtin:RestrictToScissorROI', sprintf('%i:%i:%i:%i', ox, oy, w, h));
            end

            if mystrcmp(reqs{row, 1}, 'FinalFormatting')
                % Need to restrict final formatting blit processing:
                oy = RectHeight(Screen('Rect', win, 1)) - scissorrect(RectBottom);
                Screen('HookFunction', win, 'PrependBuiltin', 'FinalOutputFormattingBlit', 'Builtin:RestrictToScissorROI', sprintf('%i:%i:%i:%i', ox, oy, w, h));
            end

        end
    end
end

% --- End of Restriction of processing area ROI ---

% Is a default colormode specified via psych_default_colormode variable and
% the level is at least 1? If so, switch to be created onscreen window to a
% [0;1] colorrange with clamping by default, and apply input scaling to
% Screen('MakeTexture') as well. This is like 'NormalizedHighresColorRange'
% aka needsUnitUnclampedColorRange, except it doesn't unclamp the
% framebuffer, but keeps it clamped to 0 - 1 range, unless a previous
% 'ColorRange' call changed this. Why? To accomodate OpenGL hw without
% clamp extension:
if ~isempty(psych_default_colormode) && (psych_default_colormode >= 1)
    applyAlsoToMakeTexture = 1;
end

if ~needsUnitUnclampedColorRange && ~isempty(psych_default_colormode) && (psych_default_colormode >= 1)
    Screen('ColorRange', win, 1, [], 1);
    applyAlsoToMakeTexture = 1;

    % Set Screen background clear color, in normalized 0.0 - 1.0 range:
    if ~isempty(clearcolor) && (max(clearcolor) > 1) && (all(round(clearcolor) == clearcolor))
        % Looks like someone's feeding old style 0-255 integer values as
        % clearcolor. Output a warning to tell about the expected 0.0 - 1.0
        % range of values:
        fprintf('\n\nPsychImaging-Warning: You specified a ''clearcolor'' argument for the OpenWindow command that looks \nlike an old 0-255 value instead of the wanted value in the 0.0-1.0 range.\nPlease update your code for correct behaviour.\n\n');
    end

    % Set the background clear color via old fullscreen 'FillRect' trick,
    % followed by a flip:
    Screen('FillRect', win, clearcolor);

    % Double-flip to be on the safe side:
    Screen('Flip', win);
    Screen('Flip', win);
end

% Do we need a normalized [0.0 ; 1.0] color range mapping with unclamped
% high precision colors?
if needsUnitUnclampedColorRange
    % Set color range to 0.0 - 1.0: This makes more sense than the normal
    % 0-255 values. Try to disable color clamping. This may fail and
    % produce a PTB warning, but if it succeeds then we're better off for
    % the 2D drawing commands...
    Screen('ColorRange', win, 1, 0, applyAlsoToMakeTexture);

    % Set Screen background clear color, in normalized 0.0 - 1.0 range:
    if ~isempty(clearcolor) && (max(clearcolor) > 1) && (all(round(clearcolor) == clearcolor)) && ~useHDR
        % Looks like someone's feeding old style 0-255 integer values as
        % clearcolor. Output a warning to tell about the expected 0.0 - 1.0
        % range of values:
        fprintf('\n\nPsychImaging-Warning: You specified a ''clearcolor'' argument for the OpenWindow command that looks \nlike an old 0-255 value instead of the wanted value in the 0.0-1.0 range.\nPlease update your code for correct behaviour.\n\n');
    end

    % Set the background clear color via old fullscreen 'FillRect' trick,
    % followed by a flip:
    Screen('FillRect', win, clearcolor);

    % Double-flip to be on the safe side:
    Screen('Flip', win);
    Screen('Flip', win);
end

% Return reqs array, for whatever reason...
rc = reqs;
return;

% End of PostConfiguration subroutine.

function rcmatch = mystrcmp(myhaystack, myneedle)

if isempty(myhaystack) || isempty(myneedle)
    rcmatch = logical(0); %#ok<LOGL>
    return;
end

if ~iscell(myhaystack) && ~ischar(myhaystack)
    error('First argument to mystrcmp must be a cell-array or a character array (string)!');
end

if iscell(myhaystack)
    % Cell array of strings: Check each element, return result matrix:
    rcmatch=logical(zeros(size(myhaystack))); %#ok<LOGL>
    rows = size(myhaystack, 1);
    cols = size(myhaystack, 2);
    for r=1:rows
        for c=1:cols
            if iscellstr(myhaystack(r,c))
                rcmatch(r,c) = logical(strcmpi(char(myhaystack(r,c)), myneedle));
            else
                rcmatch(r,c) = logical(0); %#ok<LOGL>
            end
        end
    end
else
    % Single character string: Do single check and return result:
    rcmatch=logical(strcmpi(myhaystack, myneedle));
end

return;

% Helper: Search chain 'hookname' in window 'win' for scissor restriction
% slots and remove all of them:
function DoRemoveScissorRestriction(win, hookname)
    while 1
        slot = Screen('HookFunction', win, 'Query', hookname, 'Builtin:RestrictToScissorROI');
        if slot~=-1
            Screen('Hookfunction', win, 'Remove', hookname, slot);
        else
            break;
        end
    end
return;

% Helper: Calculate and return bounding rectangle of intermediate
% framebuffers inside the imaging pipeline. These intermediates don't have
% the size of the client framebuffer (aka Screen('Rect', win);) and don't
% have the size of the windows backbuffer (aka Screen('Rect', win, 1);),
% but some size derived from the backbuffer size and various flags:
function rect = InterBufferRect(win)
    % Get window info flags about possible size transformations:
    winfo = Screen('GetWindowInfo', win);

    % Get raw rectangle of true window backbuffer size as baseline:
    % Left and Top entry is always zero, due to normalized rect.
    rect = Screen('Rect', win, 1);

    % Apply half-height flag, if any:
    if bitand(winfo.SpecialFlags, kPsychNeedHalfHeightWindow)
        rect(RectBottom) = rect(RectBottom) / 2;
    end

    % Apply half-width flag, if any:
    if bitand(winfo.SpecialFlags, kPsychNeedHalfWidthWindow)
        rect(RectRight) = rect(RectRight) / 2;
    end

    % Apply twice-width flag, if any:
    if bitand(winfo.SpecialFlags, kPsychNeedTwiceWidthWindow)
        rect(RectRight) = rect(RectRight) * 2;
    end

    % Apply triple-width flag, if any:
    if bitand(winfo.SpecialFlags, kPsychNeedTripleWidthWindow)
        rect(RectRight) = rect(RectRight) * 3;
    end
return;

function [w, h] = InterBufferSize(win)
    w = RectWidth(InterBufferRect(win));
    h = RectHeight(InterBufferRect(win));
return;

function reqs = AddTask(reqs, varargin)
    global maxreqarg;

    if length(varargin) < 2 || isempty(varargin{1}) || isempty(varargin{2})
        error('Parameters missing: Need at least "whichChannel" and "whichTask"!');
    end

    % Variable length input:
    x = varargin;

    % Store requirement in our cell array of requirements. We need to
    % extend each requirement vector to some number of max elements, so all
    % rows in the cell array have the same length:
    if length(x) < maxreqarg
        for i=length(x)+1:maxreqarg
            x{i}='';
        end
    end

    % First use of 'reqs' array?
    if isempty(reqs)
        % Yes: Initialize the array with content of 'x':
        reqs = x;
    else
        % No: Just concatenate new line with requirements 'x' to existing
        % array 'reqs':
        reqs = [reqs ; x];
    end
return;
