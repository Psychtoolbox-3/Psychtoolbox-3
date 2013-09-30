% Psychtoolbox:PsychDemos
%
% help Psychtoolbox     % For an overview, triple-click me & hit enter.
% help PsychDemos       % For demos, triple-click me & hit enter.
%
%
%    PsychExampleExperiments - Subfolder with real experiment scripts for
%                              different tasks, contributed by users.
%
%    AlphaImageDemo          - Simple demo of OpenGL alpha blending.
%    AlphaRotateDemo         - A variant of AlphaImageDemoOSX, demonstrates rotated drawing.
%    ArcDemo                 - Demonstrate the 'FrameArc', 'FillArc' functions.
%    AudioTunnel3DDemo       - Very sketchy demo on how to use OpenAL for 3D spatialized sound.
%    AudioTunnel3DDemo2      - Another sketchy demo on how to use OpenAL for 3D spatialized sound.
%    BasicAMAndMixScheduleDemo - Demonstrates basic use of sound schedules, volume controls, amplitude modulation and audio mixing of multiple voices.
%    BasicSoundFeedbackDemo  - Demonstrates a audio feedback loop via PsychPortAudio(). See DelayedSoundFeedbackDemo for a more research grade approach.
%    BasicSoundInputDemo     - Demonstrate basic usage of PsychPortAudio() for sound capture.
%    BasicSoundOutputDemo    - Demonstrate basic usage of PsychPortAudio() for sound playback.
%    BasicSoundScheduleDemo  - Demonstrate basic usage of sound schedules and buffers with PsychPortAudio().
%    BlurredMipmapDemo       - Apply blurring via OpenGL mip-mapping to live video or movies.
%    BlurredVideoCaptureDemo - Demonstrate application of GPU accelerated image filters to live video.
%    BubbleDemo              - Like GazeContingentDemoOSX, but with multiple apertures.
%    CalDemo                 - Demonstrate use of calibration structure and associated routines.
%    ClutAnimDemo            - Demonstrate lookup table animation.
%    ContrastModulatedNoiseTheClumsyStyleDemo - Demonstrate creation of contrast modulated noise.
%                                               Clumsy and inefficient, but works on old hardware.
%
%    ContrastModulatedNoiseTheElegantStyleDemo - Demonstrate creation of contrast modulated noise.
%                                                Elegant and efficient, but needs recent hardware.
%
%    DatarecordingFromISCANDemo - Template for asynchronous data collection and timestamping of gaze samples from ISCAN eye tracker.
%    DatarecordingFromSerialPortDemo - Template for asynchronous data collection and timestamping from serial port.
%    DelayedSoundFeedbackDemo - Demonstrates a audio feedback loop via PsychPortAudio() with exactly controlled latency.
%    DotDemo                 - Animate dots with Screen('DrawDots).
%    DKLDemo                 - Demonstrate DKL color space.
%    DrawFormattedTextDemo   - Demonstrate formatted text drawing via DrawFormattedText()
%    DrawHighQualityUnicodeTextDemo - Demonstrate drawing of high quality, anti-aliased and Unicode text.
%    DrawManuallyAntiAliasedTextDemo - Manually anti-alias text for special purpose applications.
%    DrawMirroredTextDemo    - Shows how to draw horizontally mirrored text.
%    DrawSomeTextDemo        - Very simple example of using 'DrawText'.
%    DriftDemo               - Animation of grating using textures.
%    DriftDemo2              - More efficient animation of gratings using one texture.
%    DriftDemo3              - Even more efficient animation of gratings under some constraints.
%    DriftDemo4              - Very efficient animation of gratings on graphics hardware with shader support.
%    DriftDemo5              - Animation of dual-gratings via use of alpha blending and color masking.
%    DriftDemo6              - Very efficient animation of dual-gratings on graphics hardware with shader support.
%    DriftWaitDemo           - Show how to optimally time animations.
%    ErrorCatchDemo          - How to handle errors gracefully (try/catch).
%    ExpandingRingsDemo      - Shows how GLSL procedural textures can be used to
%                              create some expanding rings stimulus.
%    FastNoiseDemo           - Demo and benchmark on how to quickly draw noise stimuli.
%    FastMaskedNoiseDemo     - Like FastNoiseDemo, but the noise patch is
%                              drawn through an aperture to demonstrate fast masking.
%    FitGammaDemo            - Demonstrate gamma fitting routine FitGamma.
%    FontDemo                - How to check for font availability and set font.
%    GarboriumDemo           - Demonstrate superimposed drawing of many Gabors by use of alpha-
%                              blending, 'DrawTextures' batch drawing and PTB imaging pipeline.
%    GazeContingentDemo      - Implementation of a simple gaze contingent display 
%    GLSLCLUTAnimDemo        - Demonstrates use of function moglClutBlit to
%                              perform CLUT animation with GLSL in an OS independent way.
%    GPGPUDemos              - Demos for use of GPGPU computing with Psychtoolbox.
%    GratingDemo             - Basic display of a grating on the screen.
%    IsomerizationsInDishDemo - Compute photoreceptor isomerizations for retina in a dish.
%    IsomerizationsInEyeDemo  - Compute photoreceptor isomerizations for human eye.
%    KbDemo                  - Demonstrate KbCheck, KbName, KbWait.
%    KinectDemo              - Capture and display video and depths data from a Kinect box.
%    Kinect3DDemo            - Capture data from a Kinect box and view it as a textured 3D point-cloud or mesh.
%    LineStippleDemo         - Draw dotted and dashed lines.
%    MandelbrotDemo          - Visualize the Mandelbrot fractal set by use
%                              of a GLSL procedural texture.
%    MouseTraceDemo          - Track mouse around screen and draw contour.
%    MouseTraceDemo2         - More efficient tracking of mouse around screen drawing contour.
%    MouseTraceDemo3         - Like MouseTraceDemo2, but faster, and for multiple mice.
%    MouseTraceDemo4         - Like MouseTraceDemo2, but faster on modern hardware.
%    MovieDemo               - Animation without offscreen windows.
%    MovieDemos              - Subfolder: Demos on how to use movie playback.
%    MovingLineDemo          - Shows some horizontally moving lines to demonstrate flat panel display artifacts.
%    NomogramDemo            - Compute photoreceptor nomograms.
%    OldNewRecogExp          - A full old-new recognition experiment.
%    OpenGL4MatlabDemos      - Subfolder: Demos on how to use MOGL OpenGL functions.
%    PanelFitterDemo         - Demonstrate use of the builtin panel fitter.
%    ProceduralGaborDemo     - Demo for fast drawing of Gabors via procedural textures.
%    ProceduralGarboriumDemo - Same as GarboriumDemo, just with procedural gabor drawing.
%    ProceduralNoiseDemo     - Demo for fast drawing of noise patches via procedural texturing.
%    PsychRTBoxDemo          - Demonstrates basic use of the RTBox reaction time button response box.
%    ReceivingTriggerFromSerialPortDemo - Template for asynchronous trigger collection and timestamping from serial port.
%    RenderDemo              - Render CIE xyY as a color patch.
%    SadowskiDemo            - Shows the "Sadowski" color afterimage effect.
%    SimpleSoundScheduleDemo - Simple demo for basic use of sound schedules with PsychPortAudio.
%    SimpleVoiceTriggerDemo  - Demo of a simple voice trigger with PsychPortAudio.
%    SpriteDemo              - Animates an image that follows the mouse.
%    StereoDemo              - Drive stereo-displays e.g., CrystalEyes.
%    StereoViewer            - A simple viewer for stereo image pairs.
%    StoreBit                - Not at all clear why this is here
%    ValetonVanNorrenDemo    - Exercise our code that implements the V-VN (1983) cone adaptation model.
%    VideoCaptureDemo        - Simple demo for video capture functions.
%    VignettingCorrectionDemo - Demonstrate how to do display devignetting aka per-pixel gain correction.
%
