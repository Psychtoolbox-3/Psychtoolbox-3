function varargout = PsychOculusVR(cmd, varargin)
% PsychOculusVR - A high level driver for Oculus VR hardware.
%
% Oculus VR's trademarks, e.g., Oculus, Oculus Rift, etc. are registered trademarks
% owned by Oculus VR, LLC.
%
% Note: If you want to write VR code that is portable across
% VR headsets of different vendors, then use the PsychVRHMD()
% driver instead of this driver. The PsychVRHMD driver will use
% this driver as appropriate when connecting to a Oculus Rift
% or similar Oculus device, but it will also automaticaly work
% with other head mounted displays. This driver does however
% expose a few functions specific to Oculus hardware, so you can
% mix calls to this driver with calls to PsychVRHMD to do some
% mix & match.
%
% For setup instructions for Oculus HMDs see "help OculusVR".
%
%
% Usage:
%
% oldverbosity = PsychOculusVR('Verbosity' [, newverbosity]);
% - Get/Set level of verbosity for driver status messages, warning messages,
% error messages etc. 'newverbosity' is the optional new verbosity level,
% 'oldverbosity' is the currently set verbosity level - ie. before changing
% it.  Valid settings are: 0 = Silent, 1 = Errors only, 2 = Warnings, 3 = Info,
% 4 = Debug.
%
%
% hmd = PsychOculusVR('AutoSetupHMD' [, basicTask='Tracked3DVR'][, basicRequirements][, basicQuality=0][, deviceIndex]);
% - Open a Oculus HMD, set it up with good default rendering and
% display parameters and generate a PsychImaging('AddTask', ...)
% line to setup the Psychtoolbox imaging pipeline for proper display
% on the HMD. This will also cause the device connection to get
% auto-closed as soon as the onscreen window which displays on
% the HMD is closed. Returns the 'hmd' handle of the HMD on success.
%
% By default, the first detected HMD will be used and if no VR HMD
% is connected, it will open an emulated/simulated one for basic
% testing and debugging. You can override this default choice of
% HMD by specifying the optional 'deviceIndex' parameter to choose
% a specific HMD.
%
% More optional parameters: 'basicTask' what kind of task should be implemented.
% The default is 'Tracked3DVR', which means to setup for stereoscopic 3D
% rendering, driven by head motion tracking, for a fully immersive experience
% in some kind of 3D virtual world. This is the default if omitted. The task
% 'Stereoscopic' sets up for display of stereoscopic stimuli, but without
% head tracking. 'Monoscopic' sets up for display of monocular stimuli, ie.
% the HMD is just used as a special kind of standard display monitor.
%
% 'basicRequirements' defines basic requirements for the task. Currently
% defined are the following strings which can be combined into a single
% 'basicRequirements' string: 'LowPersistence' = Try to keep exposure
% time of visual images on the retina low if possible, ie., try to approximate
% a pulse-type display instead of a hold-type display if possible. This has
% no effect on the Rift DK1. On the Rift DK2 it will enable low persistence
% scanning of the OLED display panel, to light up each pixel only a fraction
% of a video refresh cycle duration.
%
% 'PerEyeFOV' = Request use of per eye individual and asymmetric fields of view even
% when the 'basicTask' was selected to be 'Monoscopic' or 'Stereoscopic'. This allows
% for wider field of view in these tasks, but requires the usercode to adapt to these
% different and asymmetric fields of view for each eye, e.g., by selecting proper 3D
% projection matrices for each eye.
%
% 'FastResponse' = Try to switch images with minimal delay and fast
% pixel switching time. This will enable OLED panel overdrive processing
% on the Oculus Rift DK1 and DK2. OLED panel overdrive processing is a
% relatively expensive post processing step.
%
% 'TimingSupport' = Support some hardware specific means of timestamping
% or latency measurements. On the Rift DK1 this does nothing. On the DK2
% it enables dynamic prediction and timing measurements with the Rifts internal
% latency tester.
%
% 'TimeWarp' = Enable per eye image 2D timewarping via prediction of eye
% poses at scanout time. This mostly only makes sense for head-tracked 3D
% rendering. Depending on 'basicQuality' a more cheap or more expensive
% procedure is used.
%
% 'basicQuality' defines the basic tradeoff between quality and required
% computational power. A setting of 0 gives lowest quality, but with the
% lowest performance requirements. A setting of 1 gives maximum quality at
% maximum computational load. Values between 0 and 1 change the quality to
% performance tradeoff.
%
%
% hmd = PsychOculusVR('Open' [, deviceIndex], ...);
% - Open HMD with index 'deviceIndex'. See PsychOculusVRCore Open?
% for help on additional parameters.
%
%
% PsychOculusVR('SetAutoClose', hmd, mode);
% - Set autoclose mode for HMD with handle 'hmd'. 'mode' can be
% 0 (this is the default) to not do anything special. 1 will close
% the HMD 'hmd' when the onscreen window is closed which displays
% on the HMD. 2 will do the same as 1, but close all open HMDs and
% shutdown the complete driver and Oculus runtime - a full cleanup.
%
%
% isOpen = PsychOculusVR('IsOpen', hmd);
% - Returns 1 if 'hmd' corresponds to an open HMD, 0 otherwise.
%
%
% PsychOculusVR('Close' [, hmd])
% - Close provided HMD device 'hmd'. If no 'hmd' handle is provided,
% all HMDs will be closed and the driver will be shutdown.
%
%
% PsychOculusVR('Controllers', hmd);
% - Return a bitmask of all connected controllers: Can be the bitand
% of the OVR.ControllerType_XXX flags described in 'GetInputState'.
% This does not detect if controllers are hot-plugged or unplugged after
% the HMD was opened. Iow. only probed at 'Open'.
% As the classic Oculus driver does not support dedicated controllers at the
% moment, this always returns 0.
%
%
% info = PsychOculusVR('GetInfo', hmd);
% - Retrieve a struct 'info' with information about the HMD 'hmd'.
% The returned info struct contains at least the following standardized
% fields with information:
% handle = Driver internal handle for the specific HMD.
% driver = Function handle to the actual driver for the HMD, e.g., @PsychOculusVR.
% type   = Defines the type/vendor of the device, e.g., 'Oculus'.
% modelName = Name string with the name of the model of the device, e.g., 'Rift DK2'.
% separateEyePosesSupported = 1 if use of PsychOculusVR('GetEyePose') will improve
%                             the quality of the VR experience, 0 if no improvement
%                             is to be expected, so 'GetEyePose' can be avoided
%                             to save processing time without a loss of quality.
%                             This always returns 1 for at least the Rift DK1 and DK2,
%                             as use of that function can enhance the quality of the
%                             VR experience with fast head movements.
%
% The returned struct may contain more information, but the fields mentioned
% above are the only ones guaranteed to be available over the long run. Other
% fields may disappear or change their format and meaning anytime without
% warning.
%
%
% isSupported = PsychOculusVRCore('Supported');
% - Returns 1 if the Oculus driver is functional, 0 otherwise. The
% driver is functional if the VR runtime library was successfully
% initialized and a connection to the VR server process has been
% established. It would return 0 if the server process would not be
% running, or if the required runtime library would not be correctly
% installed.
%
%
% [isVisible, playAreaBounds, OuterAreaBounds] = PsychOculusVRCore('VRAreaBoundary', hmd [, requestVisible]);
% - Request visualization of the VR play area boundary for 'hmd' and returns its
% current extents.
%
% As VR area boundaries are not actually supported by this Oculus classic driver,
% this function returns no-op results, compatible with what the new Oculus driver
% would return if the Oculus guardian system would not be set up, e.g., because the
% hardware setup does not include Oculus touch controllers.
%
% The input flag 'requestVisible' is silently ignored:
% 'requestVisible' 1 = Request showing the boundary area markers, 0 = Don't
% request showing the markers.
%
% Returns in 'isVisible' the current visibility status of the VR area boundaries.
% This is always 0 for "invisible".
%
% 'playAreaBounds' is an empty matrix defining the play area boundaries. The empty
% return argument means that the play area is so far undefined on this driver.
%
% 'OuterAreaBounds' defines the outer area boundaries in the same way as
% 'playAreaBounds'. In other words, it always returns an empty matrix.
%
%
% input = PsychOculusVRCore('GetInputState', hmd, controllerType);
% - Get input state of controller 'controllerType' associated with HMD 'hmd'.
%
% As this driver does not actually support special VR controllers, only a minimally
% useful 'input' state is returned for compatibility with other drivers, which is
% based on emulating or faking input from real controllers, so this function will be
% of limited use. Specifically, only the input.Valid, input.Time and input.Buttons fields are
% returned, all other fields are missing. input.Buttons maps defined OVR.Button_XXX
% fields to similar or corresponding buttons on the regular keyboard.
%
% 'controllerType' can be one of OVR.ControllerType_LTouch, OVR.ControllerType_RTouch,
% OVR.ControllerType_Touch, OVR.ControllerType_Remote, OVR.ControllerType_XBox, or
% OVR.ControllerType_Active for selecting whatever controller is currently active.
%
% Return argument 'input' is a struct with fields describing the state of buttons and
% other input elements of the specified 'controllerType'. It has the following fields:
%
% 'Valid' = 1 if 'input' contains valid results, 0 if input status is invalid/unavailable.
% 'Time' Time of last input state change of controller.
% 'Buttons' Vector with button state on the controller, similar to the 'keyCode'
% vector returned by KbCheck() for regular keyboards. Each position in the vector
% reports pressed (1) or released (0) state of a specific button. Use the OVR.Button_XXX
% constants to map buttons to positions.
%
%
% pulseEndTime = PsychOculusVR('HapticPulse', hmd, controllerType [, duration=2.5][, freq=1.0][, amplitude=1.0]);
% - Fake triggering a haptic feedback pulse. This does nothing, but return a made up
% but consistent 'pulseEndTime', as this classic Oculus driver does not support haptic
% feedback.
%
%
% state = PsychOculusVRCore('PrepareRender', hmd [, userTransformMatrix][, reqmask=1][, targetTime]);
% - Mark the start of the rendering cycle for a new 3D rendered stereoframe.
% Return a struct 'state' which contains various useful bits of information
% for 3D stereoscopic rendering of a scene, based on head tracking data.
%
% 'hmd' is the handle of the HMD which delivers tracking data and receives the
% rendered content for display.
%
% 'reqmask' defines what kind of information is requested to be returned in
% struct 'state'. Only query information you actually need, as computing some
% of this info is expensive! See below for supported values for 'reqmask'.
%
% 'targetTime' is the expected time at which the rendered frame will display.
% This could potentially be used by the driver to make better predictions of
% camera/eye/head pose for the image. Omitting the value will use a target time
% that is implementation specific, but known to give generally good results,
% e.g., the midpoint of scanout of the next video frame.
%
% 'userTransformMatrix' is an optional 4x4 right hand side (RHS) transformation
% matrix. It gets applied to the tracked head pose as a global transformation
% before computing results based on head pose like, e.g., camera transformations.
% You can use this to translate the "virtual head" and thereby the virtual eyes/
% cameras in the 3D scene, so observer motion is not restricted to the real world
% tracking volume of your headset. A typical 'userTransformMatrix' would be a
% combined translation and rotation matrix to position the observer at some
% 3D location in space, then define his/her global looking direction, aka as
% heading angle, yaw orientation, or rotation around the y-axis in 3D space.
% Head pose tracking results would then operate relative to this global transform.
% If 'userTransformMatrix' is left out, it will default to an identity transform,
% in other words, it will do nothing.
%
%
% state always contains a field state.tracked, whose bits signal the status
% of head tracking for this frame. A +1 flag means that head orientation is
% tracked. A +2 flag means that head position is tracked via some absolute
% position tracker like, e.g., the Oculus Rift DK2 camera.
%
% state always contains a field state.tracked, whose bits signal the status
% of head tracking for this frame. A +1 flag means that head orientation is
% tracked. A +2 flag means that head position is tracked via some absolute
% position tracker like, e.g., the Oculus Rift DK2 camera. We also return a +128
% flag which means the HMD is actually strapped onto the subjects head and displaying
% our visual content. We can't detect actual HMD display state, but do this for
% compatibility to other drivers.
%
% state also always contains a field state.SessionState, whose bits signal general
% VR session status. In our case we always return +7 on this classic Oculus driver,
% as we can't detect ShouldQuit, ShouldRecenter or DisplayLost conditions, neither
% if the HMD is strapped to the users head.
%
% +1  = Our rendering goes to the HMD, ie. we have control over it. Lack of this could
%       mean the Health and Safety warning is displaying at the moment and waiting for
%       acknowledgement, or the Oculus GUI application is in control.
% +2  = HMD is present and active.
% +4  = HMD is strapped onto users head. E.g., a Oculus Rift CV1 would switch off/blank
%       if not on the head.
% +8  = DisplayLost condition! Some hardware/software malfunction, need to completely quit this
%       Psychtoolbox session to recover from this.
% +16 = ShouldQuit The user interface / user asks us to voluntarily terminate this session.
% +32 = ShouldRecenter = The user interface asks us to recenter/recalibrate our tracking origin.
%
%
% 'reqmask' defaults to 1 and can have the following values added together:
%
% +1 = Return matrices for left and right "eye cameras" which can be directly
%      used as OpenGL GL_MODELVIEW matrices for rendering the scene. 4x4 matrices
%      for left- and right eye are contained in state.modelView{1} and {2}.
%
%      Return position and orientation 4x4 camera view matrices which describe
%      position and orientation of the "eye cameras" relative to the world
%      reference frame. They are the inverses of state.modelView{}. These
%      matrices can be directly used to define cameras for rendering of complex
%      3D scenes with the Horde3D 3D engine. Left- and right eye matrices are
%      contained in state.cameraView{1} and {2}.
%
%      Additionally tracked/predicted head pose is returned in state.localHeadPoseMatrix
%      and the global head pose after application of the 'userTransformMatrix' is
%      returned in state.globalHeadPoseMatrix - this is the basis for computing
%      the camera transformation matrices.
%
% +2 = Return matrices for tracked left and right hands of user, ie. of tracked positions
%      and orientations of left and right hand tracking controllers, if any. As the old
%      driver does not support hand tracking, this reports hard-coded neutral results and
%      reports a state.handStatus of 0 = "Not tracked/Invalid data".
%
%      state.handStatus(1) = Tracking status of left hand: 0 = Untracked, signalling that
%                            all the following information is invalid and can not be used
%                            in any meaningful way.
%
%      state.handStatus(2) = Tracking status of right hand. 0 = Untracked.
%
%      state.localHandPoseMatrix{1} = 4x4 OpenGL right handed reference frame matrix with
%                                     hand position and orientation encoded to define a
%                                     proper GL_MODELVIEW transform for rendering stuff
%                                     "into"/"relative to" the oriented left hand. Always
%                                     a 4x4 unit identity matrix for hand resting in origin.
%
%      state.localHandPoseMatrix{2} = Ditto for the right hand.
%
%      state.globalHandPoseMatrix{1} = userTransformMatrix * state.localHandPoseMatrix{1};
%                                      Left hand pose transformed by passed in userTransformMatrix.
%
%      state.globalHandPoseMatrix{2} = Ditto for the right hand.
%
%      state.globalHandPoseInverseMatrix{1} = Inverse of globalHandPoseMatrix{1} for collision
%                                             testing/grasping of virtual objects relative to
%                                             hand pose of left hand.
%
%      state.globalHandPoseInverseMatrix{2} = Ditto for right hand.
%
% More flags to follow...
%
%
% eyePose = PsychOculusVR('GetEyePose', hmd, renderPass [, userTransformMatrix][, targetTime]);
% - Return a struct 'eyePose' which contains various useful bits of information
% for 3D stereoscopic rendering of the stereo view of one eye, based on head
% tracking data. This function provides essentially the same information as
% the 'PrepareRender' function, but only for one eye. Therefore you will need
% to call this function twice, once for each of the two renderpasses, at the
% beginning of each renderpass.
%
% 'hmd' is the handle of the HMD which delivers tracking data and receives the
% rendered content for display.
%
% 'renderPass' defines if information should be returned for the 1st renderpass
% (renderPass == 0) or for the 2nd renderpass (renderPass == 1). The driver will
% decide for you if the 1st renderpass should render the left eye and the 2nd
% pass the right eye, or if the 1st renderpass should render the right eye and
% then the 2nd renderpass the left eye. The ordering depends on the properties
% of the video display of your HMD, specifically on the video scanout order:
% Is it right to left, left to right, or top to bottom? For each scanout order
% there is an optimal order for the renderpasses to minimize perceived lag.
%
% 'targetTime' is the expected time at which the rendered frame will display.
% This could potentially be used by the driver to make better predictions of
% camera/eye/head pose for the image. Omitting the value will use a target time
% that is implementation specific, but known to give generally good results.
%
% 'userTransformMatrix' is an optional 4x4 right hand side (RHS) transformation
% matrix. It gets applied to the tracked head pose as a global transformation
% before computing results based on head pose like, e.g., camera transformations.
% You can use this to translate the "virtual head" and thereby the virtual eyes/
% cameras in the 3D scene, so observer motion is not restricted to the real world
% tracking volume of your headset. A typical 'userTransformMatrix' would be a
% combined translation and rotation matrix to position the observer at some
% 3D location in space, then define his/her global looking direction, aka as
% heading angle, yaw orientation, or rotation around the y-axis in 3D space.
% Head pose tracking results would then operate relative to this global transform.
% If 'userTransformMatrix' is left out, it will default to an identity transform,
% in other words, it will do nothing.
%
% Return values in struct 'eyePose':
%
% 'eyeIndex' The eye for which this information applies. 0 = Left eye, 1 = Right eye.
%            You can pass 'eyeIndex' into the Screen('SelectStereoDrawBuffer', win, eyeIndex)
%            to select the proper eye target render buffer.
%
% 'modelView' is a 4x4 RHS OpenGL matrix which can be directly used as OpenGL
%             GL_MODELVIEW matrix for rendering the scene.
%
% 'cameraView' contains a 4x4 RHS camera matrix which describes position and
%              orientation of the "eye camera" relative to the world reference
%              frame. It is the inverse of eyePose.modelView. This matrix can
%              be directly used to define the camera for rendering of complex
%              3D scenes with the Horde3D 3D engine or other engines which want
%              absolute camera pose instead of the inverse matrix.
%
% Additionally tracked/predicted head pose is returned in eyePose.localHeadPoseMatrix
% and the global head pose after application of the 'userTransformMatrix' is
% returned in eyePose.globalHeadPoseMatrix - this is the basis for computing
% the camera transformation matrix.
%
%
% PsychOculusVR('SetupRenderingParameters', hmd [, basicTask='Tracked3DVR'][, basicRequirements][, basicQuality=0][, fov=[HMDRecommended]][, pixelsPerDisplay=1])
% - Query the HMD 'hmd' for its properties and setup internal rendering
% parameters in preparation for opening an onscreen window with PsychImaging
% to display properly on the HMD. See section about 'AutoSetupHMD' above for
% the meaning of the optional parameters 'basicTask', 'basicRequirements'
% and 'basicQuality'.
%
% 'fov' Optional field of view in degrees, from line of sight: [leftdeg, rightdeg,
% updeg, downdeg]. If 'fov' is omitted, the HMD runtime will be asked for a
% good default field of view and that will be used. The field of view may be
% dependent on the settings in the HMD user profile of the currently selected
% user.
%
% 'pixelsPerDisplay' Ratio of the number of render target pixels to display pixels
% at the center of distortion. Defaults to 1.0 if omitted. Lower values can
% improve performance, at lower quality.
%
%
% PsychOculusVR('SetBasicQuality', hmd, basicQuality);
% - Set basic level of quality vs. required GPU performance.
%
%
% oldSetting = PsychOculusVR('SetFastResponse', hmd [, enable]);
% - Return old setting for 'FastResponse' mode in 'oldSetting',
% optionally disable or enable the mode via specifying the 'enable'
% parameter as 0 or greater than zero. Please note that if you want to
% use 'FastResponse', you must request and thereby enable it at the
% beginning of a session, as the driver must do some neccessary setup
% prep work at startup of the HMD. Once it was initially enabled, you
% can switch the setting at runtime with this function.
%
% Currently implemented are an algorithmic overdrive mode if 'enable'
% is set to 1, and two lookup table (LUT) based modes for 'enable'
% settings of 2 or 3, each selecting a slightly different lookup table.
%
%
% oldSetting = PsychOculusVR('SetTimeWarp', hmd [, enable]);
% - Return old setting for 'TimeWarp' mode in 'oldSetting',
% optionally enable or disable the mode via specifying the 'enable'
% parameter as 1 or 0. Please note that if you want to use 'TimeWarp',
% you must request and thereby enable it at the beginning of a session, as
% the driver must do some neccessary setup prep work at startup of the HMD.
% Once it was initially enabled, you can switch the setting at runtime with
% this function.
%
%
% oldSetting = PsychOculusVR('SetLowPersistence', hmd [, enable]);
% - Return old setting for 'LowPersistence' mode in 'oldSetting',
% optionally enable or disable the mode via specifying the 'enable'
% parameter as 1 or 0.
%
%
% oldSettings = PsychOculusVR('PanelOverdriveParameters', hmd [, newparams]);
% - Return old settings for panel overdrive mode in 'oldSettings',
% optionally set new settings in 'newparams'. This changes the operating
% parameters of OLED panel overdrive on the Rift DK-2 if 'FastResponse'
% mode is active. newparams is a vector [upscale, downscale, gamma] with
% the following meaning: gamma = 1 Use gamma/degamma pass to perform
% overdrive boost in gamma 2.2 corrected space. This is the startup default.
% upscale = How much should rising pixel color intensity values be boosted.
% Default is 0.10 for a 10% boost.
% downscale = How much should rising pixel color intensity values be reduced.
% Default is 0.05 for a 5% reduction.
% The Rift DK-2 OLED panel controller is slower on rising intensities than on
% falling intensities, therefore the higher boost on rising than on falling
% direction.
%
%
% PsychOculusVR('SetHSWDisplayDismiss', hmd [, dismissTypes=1+2+4]);
% - Set how the user can dismiss the "Health and safety warning display".
% 'dismissTypes' can be -1 to disable the HSWD, or a value >= 0 to show
% the HSWD until a timeout and or until the user dismisses the HSWD.
% The following flags can be added to define type of dismissal:
%
% +0 = Display until timeout, if any. Will wait forever if there isn't any timeout!
% +1 = Dismiss via keyboard keypress.
% +2 = Dismiss via mouse click or mousepad tap.
% +4 = Dismiss via a tap to the HMD (detected via accelerometer).
%
%
% [bufferSize, imagingFlags, stereoMode] = PsychOculusVR('GetClientRenderingParameters', hmd);
% - Retrieve recommended size in pixels 'bufferSize' = [width, height] of the client
% renderbuffer for each eye for rendering to the HMD. Returns parameters
% previously computed by PsychOculusVR('SetupRenderingParameters', hmd).
%
% Also returns 'imagingFlags', the required imaging mode flags for setup of
% the Screen imaging pipeline. Also returns the needed 'stereoMode' for the
% pipeline.
%
%
% needPanelFitter = PsychOculusVR('GetPanelFitterParameters', hmd);
% - 'needPanelFitter' is 1 if a custom panel fitter tasks is needed, and 'bufferSize'
% from the PsychVRHMD('GetClientRenderingParameters', hmd); defines the size of the
% clientRect for the onscreen window. 'needPanelFitter' is 0 if no panel fitter is
% needed.
%
%
% [winRect, ovrfbOverrideRect, ovrSpecialFlags] = PsychOculusVR('OpenWindowSetup', hmd, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags);
% - Compute special override parameters for given input/output arguments, as needed
% for a specific HMD. Take other preparatory steps as needed, immediately before the
% Screen('OpenWindow') command executes. This is called as part of PsychImaging('OpenWindow'),
% with the user provided hmd, screenid, winRect etc.
%
%
% isOutput = PsychOculusVR('IsHMDOutput', hmd, scanout);
% - Returns 1 (true) if 'scanout' describes the video output to which the
% HMD 'hmd' is connected. 'scanout' is a struct returned by the Screen
% function Screen('ConfigureDisplay', 'Scanout', screenid, outputid);
% This allows probing video outputs to find the one which feeds the HMD.
%
%
% [headToEyeShiftv, headToEyeShiftMatrix] = PsychOculusVR('GetEyeShiftVector', hmd, eye);
% - Retrieve 3D translation vector [tx, ty, tz] that defines the 3D position of the given
% eye 'eye' for the given HMD 'hmd', relative to the origin of the local head/HMD
% reference frame. This is needed to translate a global head pose into a eye
% pose, e.g., to translate the output of PsychOculusVR('GetEyePose') into actual
% tracked/predicted eye locations for stereo rendering.
%
% In addition to the 'headToEyeShiftv' vector, a corresponding 4x4 translation
% matrix is also returned in 'headToEyeShiftMatrix' for convenience.
%
%

% History:
% 07-Sep-2015  mk   Written.

% Global GL handle for access to OpenGL constants needed in setup:
global GL;
global OVR;

persistent hmd;

if nargin < 1 || isempty(cmd)
  help PsychOculusVR;
  fprintf('\n\nAlso available are functions from PsychOculusVRCore:\n');
  PsychOculusVRCore;
  return;
end

% Fast-Path function 'TimeWarp'. Prepares 2D eye timewarp:
if cmd == 1
  handle = varargin{1};

  if hmd{handle}.useOverdrive > 0
    % Find next output texture and bind it as 2nd rendertarget to the output fbo.
    % It will capture a copy of the rendered output frame, with geometry correction,
    % color aberration correction and vignette correction applied, but without the
    % overdrive processing. That copy will be used as reference for the next frame,
    % to compute per-pixel overdrive values:
    currentOverdriveTex = mod(hmd{handle}.lastOverdriveTex + 1, 2);
    glFramebufferTexture2D(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT1, GL.TEXTURE_RECTANGLE_EXT, hmd{handle}.overdriveTex(currentOverdriveTex + 1), 0);
    glDrawBuffers(2, [GL.COLOR_ATTACHMENT0, GL.COLOR_ATTACHMENT1]);

    % Bind lastOverdriveTex from previous presentation cycle as old image
    % to texture unit. It will be used for overdrive computation for this
    % frame rendercycle:
    glActiveTextureARB(GL.TEXTURE2);
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, hmd{handle}.overdriveTex(hmd{handle}.lastOverdriveTex + 1));

    % LUT based panel overdrive?
    if hmd{handle}.useOverdrive > 1
      % Bind overdrive lookup table texture to unit3 for LUT based overdrive:
      % The LUT encodes all transitions from each of the 256 possible start
      % values to each of the possible 256 end values, for each of the 3 color
      % channels, as a 256x256x4 RGBA8 texture with alpha channel unused. The
      % shader can directly use the optimal overdrive color at lut(startpix, endpix, colorchannel):
      glActiveTextureARB(GL.TEXTURE3);
      glBindTexture(GL.TEXTURE_RECTANGLE_EXT, hmd{handle}.overdriveLut(hmd{handle}.useOverdrive - 1));
    end

    % Back to standard texture unit 0:
    glActiveTextureARB(GL.TEXTURE0);

    % Prepare next rendercycle already: Swap the textures.
    hmd{handle}.lastOverdriveTex = currentOverdriveTex;
  end

  if hmd{handle}.useTimeWarp
    if hmd{handle}.useTimeWarp > 1
      % Wait for warp point, then query warp matrices. We assume the warp point is
      % 3 msecs before the target vblank and use our own high precision estimation of
      % the warp point, as well as our own high precision wait. Oculus SDK v0.5 doesn't
      % implement warp point calculation properly itself, therefore "do it yourself":
      winfo = Screen('GetWindowInfo', hmd{handle}.win, 7);
      warpPointSecs = winfo.LastVBLTime + hmd{handle}.videoRefreshDuration - 0.003;
      WaitSecs('UntilTime', warpPointSecs);
    end

    % Get the matrices:
    [hmd{handle}.eyeRotStartMatrixLeft, hmd{handle}.eyeRotEndMatrixLeft] = PsychOculusVRCore('GetEyeTimewarpMatrices', handle, 0, 0);
    [hmd{handle}.eyeRotStartMatrixRight, hmd{handle}.eyeRotEndMatrixRight] = PsychOculusVRCore('GetEyeTimewarpMatrices', handle, 1, 0);

    % Setup left shaders warp matrices:
    glUseProgram(hmd{handle}.shaderLeft(1));
    glUniformMatrix4fv(hmd{handle}.shaderLeft(2), 1, 1, hmd{handle}.eyeRotStartMatrixLeft);
    glUniformMatrix4fv(hmd{handle}.shaderLeft(3), 1, 1, hmd{handle}.eyeRotEndMatrixLeft);

    % Setup right shaders warp matrices:
    glUseProgram(hmd{handle}.shaderRight(1));
    glUniformMatrix4fv(hmd{handle}.shaderRight(2), 1, 1, hmd{handle}.eyeRotStartMatrixRight);
    glUniformMatrix4fv(hmd{handle}.shaderRight(3), 1, 1, hmd{handle}.eyeRotEndMatrixRight);

    % Ready for warp:
    glUseProgram(0);
  end

  return;
end

if cmd == 2
  handle = varargin{1};
  latencyColor = PsychOculusVRCore('LatencyTester', handle, 0);
  if ~isempty(latencyColor)
    glColor3ubv(latencyColor);
    glPointSize(4);
    glBegin(GL.POINTS);
    glVertex2i(1,1);
    glEnd;
    glPointSize(1);
  end

  return;
end

if strcmpi(cmd, 'PrepareRender')
  % Get and validate handle - fast path open coded:
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR:PrepareRender: Specified handle does not correspond to an open HMD!');
  end

  % Get 'userTransformMatrix' if any:
  if length(varargin) >= 2 && ~isempty(varargin{2})
    userTransformMatrix = varargin{2};
  else
    % Default: Identity transform to do nothing:
    userTransformMatrix = diag([1 1 1 1]);
  end

  % Valid: Get request mask of information to return:
  if length(varargin) >= 3 && ~isempty(varargin{3})
    reqmask = varargin{3};
  else
    % Default to: Provide basic tracking status flags, and directly useable
    % GL_MODELVIEW matrices for the cameras for rendering the left- and right-eye:
    reqmask = 1;
  end

  % Get target time for predicted camera poses, head poses etc.:
  if length(varargin) >= 4 && ~isempty(varargin{4})
    targetTime = varargin{4};
  else
    % Default: Provide predicted value for the midpoint of the next video
    % refresh cycle - assuming we hit the flip deadline for the next video
    % frame, so that point in time will be exactly in the middle of both
    % eyes:
    winfo = Screen('GetWindowInfo', hmd{myhmd.handle}.win);
    targetTime = winfo.LastVBLTime + 1.5 * hmd{myhmd.handle}.videoRefreshDuration;
  end

  % Mark start of a new frame render cycle for the runtime and get the data
  % predicted for next scanout time:
  [eyePose{1}, eyePose{2}, tracked] = PsychOculusVRCore('StartRender', myhmd.handle);

  % Always return basic tracking status:
  result.tracked = bitor(tracked, 128);

  % Always return faked session state:
  result.SessionState = 7; % = 1 + 2 + 4 = All is fine, no trouble, subject is attentive ;-)

  % As a bonus we return the raw eye pose vectors, given that we have them anyway:
  result.rawEyePose7{1} = eyePose{1};
  result.rawEyePose7{2} = eyePose{2};

  % Want matrices which take a usercode supplied global transformation into account?
  if bitand(reqmask, 1)
    % Yes: We need tracked + predicted head pose, so we can apply the user transform,
    % and then per-eye transforms:

    % Get predicted head pose for targetTime:
    state = PsychOculusVRCore('GetTrackingState', myhmd.handle, targetTime);

    % Bonus feature: HeadPose as 7 component translation + orientation quaternion vector:
    result.headPose = state.HeadPose;

    % Convert head pose vector to 4x4 OpenGL right handed reference frame matrix:
    result.localHeadPoseMatrix = eyePoseToCameraMatrix(state.HeadPose);

    % Premultiply usercode provided global transformation matrix:
    result.globalHeadPoseMatrix = userTransformMatrix * result.localHeadPoseMatrix;

    % Compute per-eye global pose matrices:
    result.cameraView{1} = result.globalHeadPoseMatrix * hmd{myhmd.handle}.eyeShiftMatrix{1};
    result.cameraView{2} = result.globalHeadPoseMatrix * hmd{myhmd.handle}.eyeShiftMatrix{2};

    % Compute inverse matrices, useable as OpenGL GL_MODELVIEW matrices for rendering:
    result.modelView{1} = inv(result.cameraView{1});
    result.modelView{2} = inv(result.cameraView{2});
  end

  % Want matrices with tracked position and orientation of touch controllers ~ users hands?
  if bitand(reqmask, 2)
    % Yes: We can't do this on the legacy 0.5 SDK, so fake stuff:

    for i=1:2
      result.handStatus(i) = 0;

      % Bonus feature: HandPoses as 7 component translation + orientation quaternion vectors:
      result.handPose{i} = [0, 0, 0, 0, 0, 0, 1];

      % Convert hand pose vector to 4x4 OpenGL right handed reference frame matrix:
      % In our untracked case, simply an identity matrix:
      result.localHandPoseMatrix{i} = diag([1,1,1,1]);

      % Premultiply usercode provided global transformation matrix - here use as is:
      result.globalHandPoseMatrix{i} = userTransformMatrix;

      % Compute inverse matrix, maybe useable for collision testing / virtual grasping of virtual objects:
      % Provides a transform that maps absolute geometry into geometry as "seen" from the pov of the hand.
      result.globalHandPoseInverseMatrix{i} = inv(result.globalHandPoseMatrix{i});
    end
  end

  varargout{1} = result;

  return;
end

if strcmpi(cmd, 'GetEyePose')
  % Get and validate handle - fast path open coded:
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR:GetEyePose: Specified handle does not correspond to an open HMD!');
  end

  % Valid: Get view render pass for which to return information:
  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOculusVR:GetEyePose: Required ''renderPass'' argument missing.');
  end
  renderPass = varargin{2};

  % Get 'userTransformMatrix' if any:
  if length(varargin) >= 3 && ~isempty(varargin{3})
    userTransformMatrix = varargin{3};
  else
    % Default: Identity transform to do nothing:
    userTransformMatrix = diag([1 1 1 1]);
  end

  % Get target time for predicted camera poses, head poses etc.:
  % NOTE: Currently not used, as Oculus SDK 0.5 does not support passing
  % targetTime into the underlying SDK function for 'GetEyePose'. The
  % Oculus runtime predicts something meaningful internally.
  %
  %  if length(varargin) >= 4 && ~isempty(varargin{4})
  %    targetTime = varargin{4};
  %  else
  %    % Default: Provide predicted value for the midpoint of the next video
  %    % refresh cycle - assuming we hit the flip deadline for the next video
  %    % frame, so that point in time will be exactly in the middle of both
  %    % eyes:
  %    winfo = Screen('GetWindowInfo', hmd{myhmd.handle}.win);
  %    targetTime = winfo.LastVBLTime + 1.5 * hmd{myhmd.handle}.videoRefreshDuration;
  %  end

  % Get eye pose for this renderPass, or more exactly the headPose from which this
  % renderPass eyePose will get computed:
  [result.headPose, result.eyeIndex] = PsychOculusVRCore('GetEyePose', myhmd.handle, renderPass);

  % Convert head pose vector to 4x4 OpenGL right handed reference frame matrix:
  result.localHeadPoseMatrix = eyePoseToCameraMatrix(result.headPose);

  % Premultiply usercode provided global transformation matrix:
  result.globalHeadPoseMatrix = userTransformMatrix * result.localHeadPoseMatrix;

  % Compute per-eye global pose matrix for this eyeIndex:
  result.cameraView = result.globalHeadPoseMatrix * hmd{myhmd.handle}.eyeShiftMatrix{result.eyeIndex + 1};

  % Compute inverse matrix, useable as OpenGL GL_MODELVIEW matrix for rendering:
  result.modelView = inv(result.cameraView);

  varargout{1} = result;

  return;
end

if strcmpi(cmd, 'Supported')
  % Check if the Oculus VR runtime is supported and active on this
  % installation, so it can be used to open connections to real HMDs,
  % or at least to emulate a HMD for simple debugging purposes:
  try
    if exist('PsychOculusVRCore', 'file') && PsychOculusVRCore('GetCount') >= 0
      varargout{1} = 1;
    else
      varargout{1} = 0;
    end
  catch
    varargout{1} = 0;
  end
  return;
end

if strcmpi(cmd, 'GetInputState')
  % Get and validate handle - fast path open coded:
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR:GetInputState: Specified handle does not correspond to an open HMD!');
  end

  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOculusVR:GetInputState: Required ''controllerType'' argument missing.');
  end

  rc.Valid = 1;

  [anykey, rc.Time, keyCodes] = KbCheck(-1);
  rc.Buttons = zeros(1, 32);
  if anykey
    rc.Buttons(OVR.Button_A) = keyCodes(KbName('a'));
    rc.Buttons(OVR.Button_B) = keyCodes(KbName('b'));
    rc.Buttons(OVR.Button_X) = keyCodes(KbName('x'));
    rc.Buttons(OVR.Button_Y) = keyCodes(KbName('y'));
    rc.Buttons(OVR.Button_Back) = keyCodes(KbName('BackSpace'));
    rc.Buttons(OVR.Button_Enter) = any(keyCodes(KbName('Return')));
    rc.Buttons(OVR.Button_Right) = keyCodes(KbName('RightArrow'));
    rc.Buttons(OVR.Button_Left) = keyCodes(KbName('LeftArrow'));
    rc.Buttons(OVR.Button_Up) = keyCodes(KbName('UpArrow'));
    rc.Buttons(OVR.Button_Down) = keyCodes(KbName('DownArrow'));
    rc.Buttons(OVR.Button_VolUp) = keyCodes(KbName('F12'));
    rc.Buttons(OVR.Button_VolDown) = keyCodes(KbName('F11'));
    rc.Buttons(OVR.Button_RShoulder) = keyCodes(KbName('RightShift'));
    rc.Buttons(OVR.Button_LShoulder) = keyCodes(KbName('LeftShift'));
    rc.Buttons(OVR.Button_Home) = keyCodes(KbName('Home'));
    rc.Buttons(OVR.Button_RThumb) = any(keyCodes(KbName({'RightControl', 'RightAlt'})));
    rc.Buttons(OVR.Button_LThumb) = any(keyCodes(KbName({'LeftControl', 'LeftAlt'})));
  end

  varargout{1} = rc;

  return;
end

if strcmpi(cmd, 'HapticPulse')
  % Get and validate handle - fast path open coded:
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR:HapticPulse: Specified handle does not correspond to an open HMD!');
  end

  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOculusVR:HapticPulse: Required ''controllerType'' argument missing.');
  end

  if length(varargin) >= 3 && ~isempty(varargin{3}) && varargin{3} < 2.5
    varargout{1} = WaitSecs(varargin{3});
  else
    varargout{1} = GetSecs + 2.5;
  end

  return;
end

if strcmpi(cmd, 'VRAreaBoundary')
  myhmd = varargin{1};
  if ~PsychOculusVR('IsOpen', myhmd)
    error('VRAreaBoundary: Passed in handle does not refer to a valid and open HMD.');
  end

  % Return no-op values:
  [varargout{1}, varargout{2}, varargout{3}] = deal(0, [], []);
  return;
end

% Autodetect first connected HMD and open a connection to it. Open a
% emulated one, if none can be detected. Perform basic setup with
% default configuration, create a proper PsychImaging task.
if strcmpi(cmd, 'AutoSetupHMD')
  % Do we have basic runtime support?
  if ~PsychOculusVR('Supported')
    % Nope: Game over.
    fprintf('PsychOculusVR:AutoSetupHMD: Could not initialize Oculus driver. Game over!\n');

    % Return an empty handle to signal lack of VR HMD support to caller,
    % so caller can cope with it somehow:
    varargout{1} = [];
    return;
  end

  % Basic task this HMD should fulfill:
  if length(varargin) >= 1 && ~isempty(varargin{1})
    basicTask = varargin{1};
  else
    basicTask = 'Tracked3DVR';
  end

  % Basic basicRequirements to choose:
  if length(varargin) >= 2 && ~isempty(varargin{2})
    basicRequirements = varargin{2};
  else
    basicRequirements = '';
  end

  % Basic quality/performance tradeoff to choose:
  if length(varargin) >= 3 && ~isempty(varargin{3})
    basicQuality = varargin{3};
  else
    basicQuality = 0;
  end

  % HMD device selection:
  if length(varargin) >= 4 && ~isempty(varargin{4})
    deviceIndex = varargin{4};
    newhmd = PsychOculusVR('Open', deviceIndex);
  else
    % Check if at least one Oculus HMD is connected and available:
    if PsychOculusVR('GetCount') > 0
      % Yes. Open and initialize connection to first detected HMD:
      fprintf('PsychOculusVR: Opening the first connected Oculus VR headset.\n');
      newhmd = PsychOculusVR('Open', 0);
    else
      % No. Open an emulated/simulated HMD for basic testing and debugging:
      fprintf('PsychOculusVR: No Oculus HMD detected. Opening a simulated HMD.\n');
      newhmd = PsychOculusVR('Open', -1);
    end
  end

  % Trigger an automatic device close at onscreen window close for the HMD display window:
  PsychOculusVR('SetAutoClose', newhmd, 1);

  % Setup default rendering parameters:
  PsychOculusVR('SetupRenderingParameters', newhmd, basicTask, basicRequirements, basicQuality);

  % Add a PsychImaging task to use this HMD with the next opened onscreen window:
  PsychImaging('AddTask', 'General', 'UseVRHMD', newhmd);

  % Return the device handle:
  varargout{1} = newhmd;

  % Ready.
  return;
end

if strcmpi(cmd, 'SetAutoClose')
  myhmd = varargin{1};

  if ~PsychOculusVR('IsOpen', myhmd)
    error('PsychOculusVR:SetAutoClose: Specified handle does not correspond to an open HMD!');
  end

  % Assign autoclose flag:
  hmd{myhmd.handle}.autoclose = varargin{2};

  return;
end

if strcmpi(cmd, 'SetHSWDisplayDismiss')
  myhmd = varargin{1};

  if ~PsychOculusVR('IsOpen', myhmd)
    error('PsychOculusVR:SetHSWDisplay: Specified handle does not correspond to an open HMD!');
  end

  % Method of dismissing HSW display:
  if length(varargin) < 2 || isempty(varargin{2})
    % Default is keyboard, mouse click, or HMD tap:
    hmd{myhmd.handle}.hswdismiss = 1 + 2 + 4;
  else
    hmd{myhmd.handle}.hswdismiss = varargin{2};
  end

  return;
end

% Open a HMD:
if strcmpi(cmd, 'Open')
  % Hack to make sure the VR runtime detects the HMD on a secondary X-Screen:
  if IsLinux && ~IsWayland && length(Screen('Screens')) > 1
    olddisp = getenv('DISPLAY');
    setenv('DISPLAY', sprintf(':0.%i', max(Screen('Screens'))));
  end

  [handle, modelName] = PsychOculusVRCore('Open', varargin{:});

  % Restore DISPLAY for other clients, e.g., Octave's gnuplot et al.:
  if exist('olddisp', 'var')
    setenv('DISPLAY', olddisp);
  end

  newhmd.handle = handle;
  newhmd.driver = @PsychOculusVR;
  newhmd.type   = 'Oculus';
  newhmd.subtype = 'Oculus-classic';
  newhmd.open = 1;
  newhmd.modelName = modelName;
  newhmd.separateEyePosesSupported = 1;
  newhmd.controllerTypes = 0;
  newhmd.VRControllersSupported = 0;
  newhmd.handTrackingSupported = 0;
  newhmd.hapticFeedbackSupported = 0;

  % Default autoclose flag to "no autoclose":
  newhmd.autoclose = 0;

  % Default to no use of timewarp:
  newhmd.useTimeWarp = 0;
  newhmd.readyForWarp = 0;

  % Default to no use of pixel luminance overdrive:
  newhmd.useOverdrive = 0;
  newhmd.lastOverdriveTex = -1;

  % Assign default overdrive contrast scale factors for rising
  % (UpScale) and falling (DownScale) pixel color component
  % intensities wrt. previous rendered frame:
  newhmd.overdriveUpScale   = 0.10;
  newhmd.overdriveDownScale = 0.05;

  % Perform a gamma / degamma pass on color values for a
  % gamma correction of 2.2 (hard-coded in the shader) by
  % default.
  %
  % Overdrive is optimized to operate in gamma space. As
  % we normally render and process in linear space, we
  % need to convert linear -> gamma -> Overdrive -> linear.
  % A setting of 0 for overdriveGammaCorrect would disable
  % gamma->degamma and operate purely linear:
  newhmd.overdriveGammaCorrect = 1;

  % By default allow user to dismiss HSW display via key press,
  % mouse click, or HMD tap:
  newhmd.hswdismiss = 1 + 2 + 4;

  % Setup basic task/requirement/quality specs to "nothing":
  newhmd.basicQuality = 0;
  newhmd.basicTask = '';
  newhmd.basicRequirements = '';

  if isempty(OVR)
    % Define global OVR.XXX constants:
    OVR.ControllerType_LTouch = hex2dec('0001');
    OVR.ControllerType_RTouch = hex2dec('0002');
    OVR.ControllerType_Touch = OVR.ControllerType_LTouch + OVR.ControllerType_RTouch;
    OVR.ControllerType_Remote = hex2dec('0004');
    OVR.ControllerType_XBox = hex2dec('0010');
    OVR.ControllerType_Active = hex2dec('ffffffff');

    OVR.Button_A = 1 + log2(hex2dec('00000001'));
    OVR.Button_B = 1 + log2(hex2dec('00000002'));
    OVR.Button_RThumb = 1 + log2(hex2dec('00000004'));
    OVR.Button_RShoulder = 1 + log2(hex2dec('00000008'));
    OVR.Button_X = 1 + log2(hex2dec('00000100'));
    OVR.Button_Y = 1 + log2(hex2dec('00000200'));
    OVR.Button_LThumb = 1 + log2(hex2dec('00000400'));
    OVR.Button_LShoulder = 1 + log2(hex2dec('00000800'));
    OVR.Button_Up = 1 + log2(hex2dec('00010000'));
    OVR.Button_Down = 1 + log2(hex2dec('00020000'));
    OVR.Button_Left = 1 + log2(hex2dec('00040000'));
    OVR.Button_Right = 1 + log2(hex2dec('00080000'));
    OVR.Button_Enter = 1 + log2(hex2dec('00100000'));
    OVR.Button_Back = 1 + log2(hex2dec('00200000'));
    OVR.Button_VolUp = 1 + log2(hex2dec('00400000'));
    OVR.Button_VolDown = 1 + log2(hex2dec('00800000'));
    OVR.Button_Home = 1 + log2(hex2dec('01000000'));
    OVR.Button_Private = [OVR.Button_VolUp, OVR.Button_VolDown, OVR.Button_Home];
    OVR.Button_RMask = [OVR.Button_A, OVR.Button_B, OVR.Button_RThumb, OVR.Button_RShoulder];
    OVR.Button_LMask = [OVR.Button_X, OVR.Button_Y, OVR.Button_LThumb, OVR.Button_LShoulder, OVR.Button_Enter];

    OVR.Touch_A = OVR.Button_A;
    OVR.Touch_B = OVR.Button_B;
    OVR.Touch_RThumb = OVR.Button_RThumb;
    OVR.Touch_RThumbRest = 1 + log2(hex2dec('00000008'));
    OVR.Touch_RIndexTrigger = 1 + log2(hex2dec('00000010'));
    OVR.Touch_RButtonMask = [OVR.Touch_A, OVR.Touch_B, OVR.Touch_RThumb, OVR.Touch_RThumbRest, OVR.Touch_RIndexTrigger];
    OVR.Touch_X = OVR.Button_X;
    OVR.Touch_Y = OVR.Button_Y;
    OVR.Touch_LThumb = OVR.Button_LThumb;
    OVR.Touch_LThumbRest = 1 + log2(hex2dec('00000800'));
    OVR.Touch_LIndexTrigger = 1 + log2(hex2dec('00001000'));
    OVR.Touch_LButtonMask = [OVR.Touch_X, OVR.Touch_Y, OVR.Touch_LThumb, OVR.Touch_LThumbRest, OVR.Touch_LIndexTrigger];
    OVR.Touch_RIndexPointing = 1 + log2(hex2dec('00000020'));
    OVR.Touch_RThumbUp = 1 + log2(hex2dec('00000040'));
    OVR.Touch_LIndexPointing = 1 + log2(hex2dec('00002000'));
    OVR.Touch_LThumbUp = 1 + log2(hex2dec('00004000'));
    OVR.Touch_RPoseMask =  [OVR.Touch_RIndexPointing, OVR.Touch_RThumbUp];
    OVR.Touch_LPoseMask = [OVR.Touch_LIndexPointing, OVR.Touch_LThumbUp];

    OVR.TrackedDevice_HMD        = hex2dec('0001');
    OVR.TrackedDevice_LTouch     = hex2dec('0002');
    OVR.TrackedDevice_RTouch     = hex2dec('0004');
    OVR.TrackedDevice_Touch      = OVR.TrackedDevice_LTouch + OVR.TrackedDevice_RTouch;

    OVR.TrackedDevice_Object0    = hex2dec('0010');
    OVR.TrackedDevice_Object1    = hex2dec('0020');
    OVR.TrackedDevice_Object2    = hex2dec('0040');
    OVR.TrackedDevice_Object3    = hex2dec('0080');

    OVR.TrackedDevice_All        = hex2dec('FFFF');

    OVR.KEY_USER = 'User';
    OVR.KEY_NAME = 'Name';
    OVR.KEY_GENDER = 'Gender';
    OVR.KEY_DEFAULT_GENDER = 'Unknown';
    OVR.KEY_PLAYER_HEIGHT = 'PlayerHeight';
    OVR.KEY_EYE_HEIGHT = 'EyeHeight';
    OVR.KEY_NECK_TO_EYE_DISTANCE = 'NeckEyeDistance';
    OVR.KEY_EYE_TO_NOSE_DISTANCE = 'EyeToNoseDist';

    newhmd.OVR = OVR;
    evalin('caller','global OVR');
  end

  % Store in internal array:
  hmd{handle} = newhmd;

  % Return device struct:
  varargout{1} = newhmd;
  varargout{2} = modelName;

  return;
end

if strcmpi(cmd, 'IsOpen')
  myhmd = varargin{1};
  if (length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open
    varargout{1} = 1;
  else
    varargout{1} = 0;
  end
  return;
end

if strcmpi(cmd, 'Controllers')
  myhmd = varargin{1};
  if ~PsychOculusVR('IsOpen', myhmd)
    error('Controllers: Passed in handle does not refer to a valid and open HMD.');
  end

  varargout{1} = myhmd.controllerTypes;
  return;
end

if strcmpi(cmd, 'GetInfo')
  % Ok, cheap trick: We just return the passed in 'hmd' struct - the up to date
  % internal copy that is:
  if ~PsychOculusVR('IsOpen', varargin{1})
    error('GetInfo: Passed in handle does not refer to a valid and open HMD.');
  end

  myhmd = varargin{1};
  varargout{1} = hmd{myhmd.handle};

  return;
end

if strcmpi(cmd, 'Close')
  if ~isempty(varargin) && ~isempty(varargin{1})
    % Close a specific hmd device:
    myhmd = varargin{1};

    % This function can be called with the raw index handle by
    % the autoclose code path. In that case, map index back into
    % full handle struct:
    if ~isstruct(myhmd)
      if length(hmd) >= myhmd
        myhmd = hmd{myhmd};
      else
        return;
      end
    end

    if (length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open
      PsychOculusVRCore('Close', myhmd.handle);
      hmd{myhmd.handle}.open = 0;
    end
  else
    % Shutdown whole driver:
    PsychOculusVRCore('Close');
    hmd = [];
  end

  return;
end

if strcmpi(cmd, 'IsHMDOutput')
  myhmd = varargin{1}; %#ok<NASGU>
  scanout = varargin{2};

  % Is this a Rift DK2 panel?
  if (scanout.width == 1080) && (scanout.height == 1920)
    varargout{1} = 1;
  else
    varargout{1} = 0;
  end
  return;
end

if strcmpi(cmd, 'SetBasicQuality')
  myhmd = varargin{1};
  handle = myhmd.handle;
  basicQuality = varargin{2};
  basicQuality = min(max(basicQuality, 0), 1);
  hmd{handle}.basicQuality = basicQuality;

  if ~isempty(strfind(hmd{handle}.basicRequirements, 'FastResponse'))
    hmd{handle}.useOverdrive = 1;
  else
    % Overdrive off by default because expensive:
    hmd{handle}.useOverdrive = 0;
  end

  if ~isempty(strfind(hmd{handle}.basicRequirements, 'TimeWarp'))
    if basicQuality >= 0.5
      hmd{handle}.useTimeWarp = 2;
    else
      hmd{handle}.useTimeWarp = 1;
    end
  else
    % TimeWarp is off by default:
    hmd{handle}.useTimeWarp = 0;
  end

  if ~isempty(strfind(hmd{handle}.basicRequirements, 'LowPersistence'))
    PsychOculusVRCore('SetLowPersistence', handle, 1);
  else
    PsychOculusVRCore('SetLowPersistence', handle, 0);
  end

  % Dynamic prediction enables the DK2 latency tester, advanced head tracking
  % prediction and eye timewarping:
  if ~isempty(strfind(hmd{handle}.basicRequirements, 'TimingSupport')) || ...
     hmd{handle}.useTimeWarp || ~isempty(strfind(hmd{handle}.basicTask, 'Tracked3DVR'))
    PsychOculusVRCore('SetDynamicPrediction', handle, 1);
  else
    PsychOculusVRCore('SetDynamicPrediction', handle, 0);
  end

  return;
end

if strcmpi(cmd, 'SetFastResponse')
  myhmd = varargin{1};
  if ~PsychOculusVR('IsOpen', myhmd)
    error('SetFastResponse: Passed in handle does not refer to a valid and open HMD.');
  end
  handle = myhmd.handle;

  % FastResponse determines use of GPU accelerated panel overdrive
  % on the Rift DK1/DK2. Return old setting:
  varargout{1} = hmd{handle}.useOverdrive;

  % New setting requested?
  if (length(varargin) >= 2) && ~isempty(varargin{2})
    % Check if an enable is requested, and if so, if the neccessary prep work
    % has been done during AutoSetupHMD / SetupRenderingParameters  etc. at
    % startup:
    if (varargin{2} > 0)  && (hmd{handle}.lastOverdriveTex < 0)
      error('SetFastResponse: Tried to enable fast response mode, but feature has not been requested during initial HMD setup, as required.');
    end

    % All good. Can select the new overdrive mode between 0 and 3:
    hmd{handle}.useOverdrive = max(0, min(varargin{2}, 3));

    % Set new overdrive parameters for shaders:
    if hmd{handle}.useOverdrive > 0
      if hmd{handle}.useOverdrive > 1
        % LUT based overdrive - signal to the shader via value > 1000:
        overdriveUpScale = 10000;
      else
        % Algorithmic overdrive:
        overdriveUpScale = hmd{handle}.overdriveUpScale;
      end
      overdriveDownScale = hmd{handle}.overdriveDownScale;
      overdriveGammaCorrect = hmd{handle}.overdriveGammaCorrect;
    else
      overdriveUpScale = 0;
      overdriveDownScale = 0;
      overdriveGammaCorrect = 0;
    end

    glUseProgram(hmd{handle}.shaderLeft(1));
    glUniform3f(glGetUniformLocation(hmd{handle}.shaderLeft(1), 'OverdriveScales'), overdriveUpScale, overdriveDownScale, overdriveGammaCorrect);
    glUseProgram(hmd{handle}.shaderRight(1));
    glUniform3f(glGetUniformLocation(hmd{handle}.shaderRight(1), 'OverdriveScales'), overdriveUpScale, overdriveDownScale, overdriveGammaCorrect);
    glUseProgram(0);
  end

  return;
end

if strcmpi(cmd, 'PanelOverdriveParameters')
  myhmd = varargin{1};
  if ~PsychOculusVR('IsOpen', myhmd)
    error('PanelOverdriveParameters: Passed in handle does not refer to a valid and open HMD.');
  end
  handle = myhmd.handle;

  % PanelOverdriveParameters determines the parameters of GPU accelerated panel overdrive
  % on the Rift DK1/DK2. Return old setting:
  varargout{1} = [hmd{handle}.overdriveUpScale, hmd{handle}.overdriveDownScale, hmd{handle}.overdriveGammaCorrect];

  % New setting requested?
  if (length(varargin) >= 2) && ~isempty(varargin{2})
    % Set new overdrive parameters for shaders:
    newparams = varargin{2};
    if length(newparams) ~= 3
      error('PanelOverdriveParameters: Invalid new overdrive parameters. Not a 3-component vector [upscale, downscale, gamma].');
    end

    hmd{handle}.overdriveUpScale = newparams(1);
    hmd{handle}.overdriveDownScale = newparams(2);
    hmd{handle}.overdriveGammaCorrect = newparams(3);

    if hmd{handle}.useOverdrive > 1
      % LUT based overdrive - signal to the shader via value > 1000:
      overdriveUpScale = 10000;
    else
      % Algorithmic overdrive:
      overdriveUpScale = hmd{handle}.overdriveUpScale;
    end

    overdriveDownScale = hmd{handle}.overdriveDownScale;
    overdriveGammaCorrect = hmd{handle}.overdriveGammaCorrect;

    glUseProgram(hmd{handle}.shaderLeft(1));
    glUniform3f(glGetUniformLocation(hmd{handle}.shaderLeft(1), 'OverdriveScales'), overdriveUpScale, overdriveDownScale, overdriveGammaCorrect);
    glUseProgram(hmd{handle}.shaderRight(1));
    glUniform3f(glGetUniformLocation(hmd{handle}.shaderRight(1), 'OverdriveScales'), overdriveUpScale, overdriveDownScale, overdriveGammaCorrect);
    glUseProgram(0);
  end

  return;
end

if strcmpi(cmd, 'SetTimeWarp')
  myhmd = varargin{1};
  if ~PsychOculusVR('IsOpen', myhmd)
    error('SetTimeWarp: Passed in handle does not refer to a valid and open HMD.');
  end

  % SetTimeWarp determines use of GPU accelerated 2D texture sampling
  % warp on the Rift DK1/DK2. Return old setting:
  varargout{1} = hmd{myhmd.handle}.useTimeWarp;

  % New setting requested?
  if (length(varargin) >= 2) && ~isempty(varargin{2})
    % Check if an enable is requested, and if so, if the neccessary prep work
    % has been done during AutoSetupHMD / SetupRenderingParameters  etc. at
    % startup:
    if (varargin{2} > 0)  && ~hmd{myhmd.handle}.readyForWarp
      error('SetTimeWarp: Tried to enable eye timewarp mode, but feature has not been requested during initial HMD setup, as required.');
    end

    % TimeWarp transition from enabled to disabled?
    if (varargin{2} <= 0) && (hmd{myhmd.handle}.useTimeWarp > 0)
      % Need to reset shaders matrices to identity matrices:
      handle = myhmd.handle;

      % Setup left shaders warp matrices:
      glUseProgram(hmd{handle}.shaderLeft(1));
      hmd{handle}.eyeRotStartMatrixLeft = diag([1 1 1 1]);
      hmd{handle}.eyeRotEndMatrixLeft   = diag([1 1 1 1]);
      glUniformMatrix4fv(hmd{handle}.shaderLeft(2), 1, 1, hmd{handle}.eyeRotStartMatrixLeft);
      glUniformMatrix4fv(hmd{handle}.shaderLeft(3), 1, 1, hmd{handle}.eyeRotEndMatrixLeft);

      % Setup right shaders warp matrices:
      glUseProgram(hmd{handle}.shaderRight(1));
      hmd{handle}.eyeRotStartMatrixRight = diag([1 1 1 1]);
      hmd{handle}.eyeRotEndMatrixRight   = diag([1 1 1 1]);
      glUniformMatrix4fv(hmd{handle}.shaderRight(2), 1, 1, hmd{handle}.eyeRotStartMatrixRight);
      glUniformMatrix4fv(hmd{handle}.shaderRight(3), 1, 1, hmd{handle}.eyeRotEndMatrixRight);
      glUseProgram(0);
    end

    % All good. Can select the new timeWarp mode:
    hmd{myhmd.handle}.useTimeWarp = varargin{2};
  end

  return;
end

if strcmpi(cmd, 'SetLowPersistence')
  myhmd = varargin{1};
  if ~PsychOculusVR('IsOpen', myhmd)
    error('SetLowPersistence: Passed in handle does not refer to a valid and open HMD.');
  end

  % SetLowPersistence determines use low persistence mode on the Rift DK2. Return old setting:
  varargout{1} = PsychOculusVRCore('SetLowPersistence', myhmd.handle);

  % New setting requested?
  if (length(varargin) >= 2) && ~isempty(varargin{2})
    PsychOculusVRCore('SetLowPersistence', myhmd.handle, varargin{2});
  end

  return;
end

if strcmpi(cmd, 'SetupRenderingParameters')
  myhmd = varargin{1};

  % Basic task this HMD should fulfill:
  if length(varargin) >= 2 && ~isempty(varargin{2})
    basicTask = varargin{2};
  else
    basicTask = 'Tracked3DVR';
  end

  % Basic requirements to choose:
  if length(varargin) >= 3 && ~isempty(varargin{3})
    basicRequirements = varargin{3};
  else
    basicRequirements = '';
  end

  % Basic quality/performance tradeoff to choose:
  if length(varargin) >= 4 && ~isempty(varargin{4})
    basicQuality = varargin{4};
  else
    basicQuality = 0;
  end

  hmd{myhmd.handle}.basicTask = basicTask;
  hmd{myhmd.handle}.basicRequirements = basicRequirements;

  PsychOculusVR('SetBasicQuality', myhmd, basicQuality);

  % Get optimal client renderbuffer size - the size of our virtual framebuffer for left eye:
  [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovL] = PsychOculusVRCore('GetFovTextureSize', myhmd.handle, 0, varargin{5:end});

  % Get optimal client renderbuffer size - the size of our virtual framebuffer for right eye:
  [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovR] = PsychOculusVRCore('GetFovTextureSize', myhmd.handle, 1, varargin{5:end});

  % If the basic task is not a 3D VR rendering one (with or without HMD tracking),
  % and the special requirement 'PerEyeFOV' is not set, then assume usercode wants
  % to do pure 2D rendering (monocular, or stereoscopic), e.g., with the Screen()
  % 2D drawing commands, and doesn't set up per-eye projection and modelview matrices.
  % In this case we must use a field of view that is identical for both eyes, and
  % both vertically and horizontally symmetric, ie. no special treatment of the nose
  % facing field of view! Why? Because standard 2D mono/stereo drawing code doesn't
  % know about/can't use per eye view projection matrices, which are needed for proper
  % results for asymmetric per-eye FOV. It would cause weird shifts in display on the
  % HMD. This effect is almost imperceptible/negligible on the Rift DK1/DK2, but very
  % disturbing on the Rift CV1.
  if isempty(strfind(hmd{myhmd.handle}.basicTask, 'Tracked3DVR')) && ...
     isempty(strfind(hmd{myhmd.handle}.basicTask, '3DVR')) && ...
     isempty(strfind(hmd{myhmd.handle}.basicRequirements, 'PerEyeFOV'))
    % Need identical, symmetric FOV for both eyes. Build one that has the same
    % vertical FOV as proposed by the runtime, but horizontally uses the minimal
    % left/right FOV extension of both per-eye FOV's, so we get a symmetric FOV
    % identical for both eyes, guaranteed to lie within the view cone not occluded
    % by the nose of the user.
    fov(1) = min(hmd{myhmd.handle}.fovL(1), hmd{myhmd.handle}.fovR(1));
    fov(2) = min(hmd{myhmd.handle}.fovL(2), hmd{myhmd.handle}.fovR(2));
    fov(3) = min(hmd{myhmd.handle}.fovL(3), hmd{myhmd.handle}.fovR(3));
    fov(4) = min(hmd{myhmd.handle}.fovL(4), hmd{myhmd.handle}.fovR(4));

    % Recompute parameters based on override fov:
    [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovL] = PsychOculusVRCore('GetFovTextureSize', myhmd.handle, 0, fov, varargin{6:end});
    [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovR] = PsychOculusVRCore('GetFovTextureSize', myhmd.handle, 1, fov, varargin{6:end});
  end

  return;
end

if strcmpi(cmd, 'GetClientRenderingParameters')
  myhmd = varargin{1};
  varargout{1} = [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight];

  % We need fast backing store support for virtual framebuffers:
  imagingMode = mor(kPsychNeedTwiceWidthWindow, kPsychNeedFastBackingStore);
  imagingMode = mor(imagingMode, kPsychNeedClientRectNoFitter);

  % Need an output FBO for our panel overdrive implementation:
  if hmd{myhmd.handle}.useOverdrive || strcmpi(hmd{myhmd.handle}.basicTask, 'Monoscopic')
    imagingMode = mor(imagingMode, kPsychNeedOutputConversion);
  end

  if ~strcmpi(hmd{myhmd.handle}.basicTask, 'Monoscopic')
    % We must use stereomode 6, so we get separate draw buffers for left and
    % right eye, and the stereo compositor (merger) to fuse both eyes into a
    % single output framebuffer, but with all internal buffers at at least
    % full output framebuffer resolution. This will generate anaglyph shaders
    % which we will need to replace with a very special shader for the HMD:
    stereoMode = 6;
  else
    % Monoscopic presentation will do:
    stereoMode = 0;
  end

  varargout{2} = imagingMode;
  varargout{3} = stereoMode;
  return;
end

if strcmpi(cmd, 'GetPanelFitterParameters')
  % We need use of custom PanelFitter task for the v0.5 Oculus SDK/Runtime:
  varargout{1} = 1;
  return;
end

% [winRect, ovrfbOverrideRect, ovrSpecialFlags] = PsychOculusVR('OpenWindowSetup', hmd, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags);
if strcmpi(cmd, 'OpenWindowSetup')
  myhmd = varargin{1};
  screenid = varargin{2};
  winRect = varargin{3};
  ovrfbOverrideRect = varargin{4};
  ovrSpecialFlags = varargin{5};

  % Yes. Trying to display on a screen with more than one video output?
  if isempty(winRect) && (Screen('ConfigureDisplay', 'NumberOutputs', screenid) > 1)
    % Yes. Not good, as this will impair graphics performance and timing a lot.
    % Warn about this, then try to at least position the onscreen window on the
    % right output.
    fprintf('PsychOculusVR-WARNING: You are requesting display to a VR HMD on a screen with multiple active video outputs.\n');
    fprintf('PsychOculusVR-WARNING: This will impair visual stimulation timing and cause decreased VR performance!\n');
    fprintf('PsychOculusVR-WARNING: I strongly recommend only activating one output on the HMD screen - the HMD output on the screen.\n');
    fprintf('PsychOculusVR-WARNING: On Linux with X11 X-Server, you should create a separate X-Screen for the HMD.\n');

    % Try to find the output with the Rift HMD:
    for i=0:Screen('ConfigureDisplay', 'NumberOutputs', screenid)-1
      scanout = Screen('ConfigureDisplay', 'Scanout', screenid, i);
      if myhmd.driver('IsHMDOutput', myhmd, scanout)
        % This output i has proper resolution to be the HMD panel.
        % Position our onscreen window accordingly:
        winRect = OffsetRect([0, 0, scanout.width, scanout.height], scanout.xStart, scanout.yStart);
        fprintf('PsychOculusVR-Info: Positioning onscreen window at rect [%i, %i, %i, %i] to align with HMD output %i.\n', ...
                winRect(1), winRect(2), winRect(3), winRect(4), i);
      end
    end
  end

  varargout{1} = winRect;
  varargout{2} = ovrfbOverrideRect;
  varargout{3} = ovrSpecialFlags;

  return;
end

if strcmpi(cmd, 'GetEyeShiftVector')
  myhmd = varargin{1};

  if varargin{2} == 0
    varargout{1} = hmd{myhmd.handle}.HmdToEyeViewOffsetLeft;
    varargout{2} = hmd{myhmd.handle}.eyeShiftMatrix{1};
  else
    varargout{1} = hmd{myhmd.handle}.HmdToEyeViewOffsetRight;
    varargout{2} = hmd{myhmd.handle}.eyeShiftMatrix{2};
  end

  return;
end

if strcmpi(cmd, 'PerformPostWindowOpenSetup')

  % Must have global GL constants:
  if isempty(GL)
    varargout{1} = 0;
    warning('PTB internal error in PsychOculusVR: GL struct not initialized?!?');
    return;
  end

  % Oculus device handle:
  myhmd = varargin{1};
  handle = myhmd.handle;

  % Onscreen window handle:
  win = varargin{2};

  % Keep track of window handle of associated onscreen window:
  hmd{handle}.win = win;

  % Need to know user selected clearcolor:
  clearcolor = varargin{3};

  % Also keep track of video refresh duration of the HMD:
  hmd{handle}.videoRefreshDuration = Screen('Framerate', win);
  if hmd{handle}.videoRefreshDuration == 0
    % Unlikely to ever hit this situation, but if we would, just
    % default to the Rift DK-2's default video refresh rate of 75 Hz:
    hmd{handle}.videoRefreshDuration = 75;
  end
  hmd{handle}.videoRefreshDuration = 1 / hmd{handle}.videoRefreshDuration;

  % Compute effective size of per-eye input buffer for undistortion render.
  % The input buffers for undistortion are the processedDrawbufferFBO's aka
  % inputBufferFBO's, or if the panelfitter is skipped the drawBufferFBO's.
  %
  % In our current implementation we allocate said buffers to twice the horizontal
  % size of the real framebuffer, ie., twice the panel width of the HMD, as
  % that should be plenty for all typical use cases - and is also the maximum
  % possible with the current Screen imaging pipeline.
  %
  % However, we don't use the full size of those buffers as input, but only
  % sample a rectangular subregion which corresponds to the renderbuffer size
  % recommended by the Oculus runtime. Either the panelfitter is used to blit
  % 1-to-1 from the drawBufferFBO to a correspondingly sized subregion of the
  % inputBuffers - if the panelfitter is needed for convenient 2D stimulus drawing
  % or MSAA resolve - or usercode has to restrict its rendering to the subregion by
  % proper use of glViewPorts or scissor rectangles.
  %
  % So for all practical means [inputWidth, inputHeight] == [rbwidth, rbheight] and
  % we save processing bandwidth, although due to the overallocation not VRAM memory
  % space.
  hmd{handle}.inputWidth = hmd{handle}.rbwidth;
  hmd{handle}.inputHeight = hmd{handle}.rbheight;

  % Query undistortion parameters for left eye view:
  [hmd{handle}.rbwidth, hmd{handle}.rbheight, vx, vy, vw, vh, ptx, pty, hsx, hsy, hsz, meshVL, meshIL, uvScale(1), uvScale(2), uvOffset(1), uvOffset(2)] = PsychOculusVRCore('GetUndistortionParameters', handle, 0, hmd{handle}.inputWidth, hmd{handle}.inputHeight, hmd{handle}.fovL);
  hmd{handle}.viewportLeft = [vx, vy, vw, vh];
  hmd{handle}.PixelsPerTanAngleAtCenterLeft = [ptx, pty];
  hmd{handle}.HmdToEyeViewOffsetLeft = -1 * [hsx, hsy, hsz];
  hmd{handle}.meshVerticesLeft = meshVL;
  hmd{handle}.meshIndicesLeft = meshIL;
  hmd{handle}.uvScaleLeft = uvScale;
  hmd{handle}.uvOffsetLeft = uvOffset;

  % Init warp matrices to identity, until we get something better from live tracking:
  hmd{handle}.eyeRotStartMatrixLeft = diag([1 1 1 1]);
  hmd{handle}.eyeRotEndMatrixLeft   = diag([1 1 1 1]);

  % Query parameters for right eye view:
  [hmd{handle}.rbwidth, hmd{handle}.rbheight, vx, vy, vw, vh, ptx, pty, hsx, hsy, hsz, meshVR, meshIR, uvScale(1), uvScale(2), uvOffset(1), uvOffset(2)] = PsychOculusVRCore('GetUndistortionParameters', handle, 1, hmd{handle}.inputWidth, hmd{handle}.inputHeight, hmd{handle}.fovR);
  hmd{handle}.viewportRight = [vx, vy, vw, vh];
  hmd{handle}.PixelsPerTanAngleAtCenterRight = [ptx, pty];
  hmd{handle}.HmdToEyeViewOffsetRight = -1 * [hsx, hsy, hsz];
  hmd{handle}.meshVerticesRight = meshVR;
  hmd{handle}.meshIndicesRight = meshIR;
  hmd{handle}.uvScaleRight = uvScale;
  hmd{handle}.uvOffsetRight = uvOffset;

  % Init warp matrices to identity, until we get something better from live tracking:
  hmd{handle}.eyeRotStartMatrixRight = diag([1 1 1 1]);
  hmd{handle}.eyeRotEndMatrixRight   = diag([1 1 1 1]);

  % Convert head to eye shift vectors into 4x4 matrices, as we'll need
  % them frequently:
  EyeT = diag([1 1 1 1]);
  EyeT(1:3, 4) = hmd{handle}.HmdToEyeViewOffsetLeft';
  hmd{handle}.eyeShiftMatrix{1} = EyeT;

  EyeT = diag([1 1 1 1]);
  EyeT(1:3, 4) = hmd{handle}.HmdToEyeViewOffsetRight';
  hmd{handle}.eyeShiftMatrix{2} = EyeT;

  % Switch to clear color black and do a clear by double flip:
  Screen('FillRect', win, 0);
  Screen('Flip', win);
  Screen('Flip', win);

  % Assign proper target processing chain for imaging pipeline:
  if ~strcmpi(hmd{handle}.basicTask, 'Monoscopic')
    % Stereoscopic display: Stereo composer chain.
    procchain = 'StereoCompositingBlit';

    % Find slot with preexisting stereo composition shader, so we can replace it by our thing:
    [slot, shaderid, blittercfg, voidptr, glsl] = Screen('HookFunction', win, 'Query', procchain, 'StereoCompositingShaderAnaglyph'); %#ok<ASGLU>
    if slot == -1
      varargout{1} = 0;
      warning('Either the imaging pipeline is not enabled for given onscreen window, or it is not switched to Anaglyph stereo mode.');
      return;
    end

    if glsl == 0
      varargout{1} = 0;
      warning('Anaglyph shader is not operational for unknown reason. Sorry...');
      return;
    end

    % Remove old standard anaglyph shader:
    Screen('HookFunction', win, 'Remove', procchain, slot);

    % Play more stupid tricks to get intermediate (bounce buffer FBOs) buffers cleared to black:
    Screen('HookFunction', win, 'AppendBuiltin', procchain, 'Builtin:IdentityBlit', '');
    Screen('Flip', win);
    Screen('Flip', win);
    Screen('HookFunction', win, 'Remove', procchain, slot);
  else
    % Monoscopic display: Final output formatter:
    procchain = 'FinalOutputFormattingBlit';
    Screen('HookFunction', win, 'Enable', procchain);

    % For overdrive need stupid tricks to get intermediate bounce buffer FBO's cleared:
    if hmd{handle}.useOverdrive
      % Need a bufferflip command:
      Screen('HookFunction', win, 'PrependBuiltin', procchain, 'Builtin:FlipFBOs', '');
      Screen('HookFunction', win, 'PrependBuiltin', procchain, 'Builtin:IdentityBlit', '');
      Screen('Flip', win);
      Screen('HookFunction', win, 'Remove', procchain, 0);
      Screen('HookFunction', win, 'Remove', procchain, 0);
    end

    slot = 0;
    glsl = 0;
  end

  % Go back to user requested clear color, now that all our buffers
  % are cleared to black:
  Screen('FillRect', win, clearcolor);

  % Build the unwarp mesh display list within the OpenGL context of Screen():
  Screen('BeginOpenGL', win, 1);

  % Left eye setup:
  % ---------------

  % Build a display list that corresponds to the current calibration,
  % drawing the warp-mesh once, so it gets recorded in the display list:
  gldLeft = glGenLists(1);
  glNewList(gldLeft, GL.COMPILE);

  % Caution: Must *copy* the different rows with data into *separate* variables, so
  % the vertex array pointers to the different variables actually point to something
  % persistent! If we'd pass the meshVerticesLeft() subarrays directly to glTexCoordPointer
  % and friends then Octave/Matlab would just create a temporary copy of the extracted
  % rows, OpenGL would retrieve/assign pointers to those temporary copies, but then
  % at the end of a glVertexPointer/glTexCoordPointer call, those temporary copies would
  % go out of scope and Octave/Matlab would potentially garbage collect the variables again
  % *before* the call to glDrawElements permanently records the content of the variables.
  % The net results would be stale/dangling pointers, random data trash getting read from
  % memory and recorded in the display list - and thereby corrupted rendering! This hazard
  % doesn't exist within regular Octave/Matlab scripts, because the interpreter doesn't
  % deal with memory pointers. It is a unique hazard from the combination of C memory
  % pointers for OpenGL and Octave/Matlabs copy-on-write/data-sharing/garbage collection
  % behaviour. When we are at it, lets also cast the data to single() precision floating
  % point, to save some memory:
  vertexpos = single(hmd{handle}.meshVerticesLeft(1:4, :));

  if ~IsLinux
      % Both Windows and OSX need special treatment, because the 0.5 SDK
      % doesn't generate a properly rotated undistortion mesh. Rotate
      % vertex (x,y) positions by 90 degrees counter-clockwise, so the mesh
      % aligns with the 90 degrees rotated full HD panel of the Rift DK-1
      % and DK-2. This allows to keep the video mode on at, e.g. for the
      % DK-2, native 1080 x 1920 without enabling output rotation. That in
      % turn keeps page flipping enabled for bufferswaps, at least on the
      % non-broken graphics drivers, and that in turn keeps PTB's timing
      % happy and performance up:
      R = single([0, -1 ; 1, 0]);
      vertexpos(1:2, :) = R * vertexpos(1:2, :);
  end

  texR = single(hmd{handle}.meshVerticesLeft(5:6, :));
  texG = single(hmd{handle}.meshVerticesLeft(7:8, :));
  texB = single(hmd{handle}.meshVerticesLeft(9:10, :));

  % vertex xy encodes 2D position from rows 1 and 2, z encodes timeWarp interpolation factors
  % from row 3 and w encodes vignette correction factors from row 4:
  glEnableClientState(GL.VERTEX_ARRAY);
  glVertexPointer(4, GL.FLOAT, 0, vertexpos);

  % Need separate texture coordinate sets for the three color channel to encode
  % channel specific color aberration correction sampling:

  % TexCoord set 0 encodes coordinates for the Red color channel:
  glClientActiveTexture(GL.TEXTURE0);
  glEnableClientState(GL.TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL.FLOAT, 0, texR);
  
  % TexCoord set 1 encodes coordinates for the Green color channel:
  glClientActiveTexture(GL.TEXTURE1);
  glEnableClientState(GL.TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL.FLOAT, 0, texG);

  % TexCoord set 2 encodes coordinates for the Blue color channel:
  glClientActiveTexture(GL.TEXTURE2);
  glEnableClientState(GL.TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL.FLOAT, 0, texB);

  % Draw the mesh. This records the content from all the variables persistently into
  % the display list storage, so they can be freed afterwards:
  glDrawElements(GL.TRIANGLES, length(hmd{handle}.meshIndicesLeft), GL.UNSIGNED_SHORT, uint16(hmd{handle}.meshIndicesLeft));

  % Disable stuff, so we can release or recycle the variables:
  glClientActiveTexture(GL.TEXTURE3);
  glDisableClientState(GL.TEXTURE_COORD_ARRAY);

  glClientActiveTexture(GL.TEXTURE2);
  glDisableClientState(GL.TEXTURE_COORD_ARRAY);

  glClientActiveTexture(GL.TEXTURE1);
  glDisableClientState(GL.TEXTURE_COORD_ARRAY);

  glClientActiveTexture(GL.TEXTURE0);
  glDisableClientState(GL.TEXTURE_COORD_ARRAY);

  glDisableClientState(GL.VERTEX_ARRAY);
  
  % Left eye display list done.
  glEndList;

  % Right eye setup:
  % ---------------

  % Build a display list that corresponds to the current calibration,
  % drawing the warp-mesh once, so it gets recorded in the display list:
  gldRight = glGenLists(1);
  glNewList(gldRight, GL.COMPILE);

  vertexpos = single(hmd{handle}.meshVerticesRight(1:4, :));

  if ~IsLinux
      % Same special treatment on non-Linux as for the left eye. Rotate mesh by
      % 90 degrees counter-clockwise:
      vertexpos(1:2, :) = R * vertexpos(1:2, :);
  end

  texR = single(hmd{handle}.meshVerticesRight(5:6, :));
  texG = single(hmd{handle}.meshVerticesRight(7:8, :));
  texB = single(hmd{handle}.meshVerticesRight(9:10, :));

  % vertex xy encodes 2D position from rows 1 and 2, z encodes timeWarp interpolation factors
  % from row 3 and w encodes vignette correction factors from row 4:
  glEnableClientState(GL.VERTEX_ARRAY);
  glVertexPointer(4, GL.FLOAT, 0, vertexpos);

  % Need separate texture coordinate sets for the three color channel to encode
  % channel specific color aberration correction sampling:

  % TexCoord set 0 encodes coordinates for the Red color channel:
  glClientActiveTexture(GL.TEXTURE0);
  glEnableClientState(GL.TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL.FLOAT, 0, texR);
  
  % TexCoord set 1 encodes coordinates for the Green color channel:
  glClientActiveTexture(GL.TEXTURE1);
  glEnableClientState(GL.TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL.FLOAT, 0, texG);

  % TexCoord set 2 encodes coordinates for the Blue color channel:
  glClientActiveTexture(GL.TEXTURE2);
  glEnableClientState(GL.TEXTURE_COORD_ARRAY);
  glTexCoordPointer(2, GL.FLOAT, 0, texB);

  % Draw the mesh. This records the content from all the variables persistently into
  % the display list storage, so they can be freed afterwards:
  glDrawElements(GL.TRIANGLES, length(hmd{handle}.meshIndicesRight), GL.UNSIGNED_SHORT, uint16(hmd{handle}.meshIndicesRight));

  % Disable stuff, so we can release or recycle the variables:
  glClientActiveTexture(GL.TEXTURE3);
  glDisableClientState(GL.TEXTURE_COORD_ARRAY);

  glClientActiveTexture(GL.TEXTURE2);
  glDisableClientState(GL.TEXTURE_COORD_ARRAY);

  glClientActiveTexture(GL.TEXTURE1);
  glDisableClientState(GL.TEXTURE_COORD_ARRAY);

  glClientActiveTexture(GL.TEXTURE0);
  glDisableClientState(GL.TEXTURE_COORD_ARRAY);

  glDisableClientState(GL.VERTEX_ARRAY);
  
  % Right eye display list done.
  glEndList;

  Screen('EndOpenGL', win);

  if hmd{handle}.useOverdrive
    % Overdrive enabled: Assign overdrive contrast scale factors for
    % rising (UpScale) and falling (DownScale) pixel color component
    % intensities wrt. previous rendered frame:
    if hmd{handle}.useOverdrive > 1
      % LUT based overdrive - signal to the shader via value > 1000:
      overdriveUpScale = 10000;
    else
      % Algorithmic overdrive:
      overdriveUpScale = hmd{handle}.overdriveUpScale;
    end

    overdriveDownScale = hmd{handle}.overdriveDownScale;

    % Perform a gamma / degamma pass on color values for a
    % gamma correction of 2.2 (hard-coded in the shader).
    % Overdrive is optimized to operate in gamma space. As
    % we normally render and process in linear space, we
    % need to convert linear -> gamma -> Overdrive -> linear.
    % A setting of 0 for overdriveGammaCorrect would disable
    % gamma->degamma and operate purely linear:
    overdriveGammaCorrect = hmd{handle}.overdriveGammaCorrect;
  else
    % Overdrive disabled:
    overdriveUpScale = 0;
    overdriveDownScale = 0;
    overdriveGammaCorrect = 0;
  end

  % Setup left eye shader:
  glsl = LoadGLSLProgramFromFiles([fileparts(mfilename('fullpath')) filesep 'OculusRiftCorrectionShader']);
  glUseProgram(glsl);
  glUniform1i(glGetUniformLocation(glsl, 'Image'), 0);
  glUniform1i(glGetUniformLocation(glsl, 'PrevImage'), 2);
  glUniform1i(glGetUniformLocation(glsl, 'OverdriveLUT'), 3);
  glUniform3f(glGetUniformLocation(glsl, 'OverdriveScales'), overdriveUpScale, overdriveDownScale, overdriveGammaCorrect);
  glUniform2f(glGetUniformLocation(glsl, 'EyeToSourceUVOffset'), hmd{handle}.uvOffsetLeft(1) * hmd{handle}.inputWidth, hmd{handle}.uvOffsetLeft(2) * hmd{handle}.inputHeight);
  glUniform2f(glGetUniformLocation(glsl, 'EyeToSourceUVScale'), hmd{handle}.uvScaleLeft(1) * hmd{handle}.inputWidth, hmd{handle}.uvScaleLeft(2) * hmd{handle}.inputHeight);
  glUniformMatrix4fv(glGetUniformLocation(glsl, 'EyeRotationStart'), 1, 1, hmd{handle}.eyeRotStartMatrixLeft);
  glUniformMatrix4fv(glGetUniformLocation(glsl, 'EyeRotationEnd'), 1, 1, hmd{handle}.eyeRotEndMatrixLeft);
  hmd{handle}.shaderLeft = [glsl, glGetUniformLocation(glsl, 'EyeRotationStart'), glGetUniformLocation(glsl, 'EyeRotationEnd')];
  glUseProgram(0);

  % Insert it at former position of the old shader:
  posstring = sprintf('InsertAt%iShader', slot);
  
  % xOffset and yOffset encode the viewport location and size for the left-eye vs.
  % right eye view in the shared output window - or the source renderbuffer if both eyes
  % would be rendered into a shared texture. However, the meshes provided by the SDK
  % already encode proper left and right offsets for output, and the inputs are separate
  % textures for left and right eye, so using the offset is not needed. Also our correction
  % shader ignores the modelview matrix which would get updated with the "Offset:%i%i" blittercfg,
  % instead is takes normalized device coordinates NDC directly from the distortion mesh. Iow, not
  % only is xOffset/yOffset not needed, it would also be a no operation due to our specific shader.
  % We leave this here for documentation for now, in case we need to change our ways of doing this.
  %leftViewPort = hmd{handle}.viewportLeft
  blittercfg = sprintf('Blitter:DisplayListBlit:Handle:%i:Bilinear', gldLeft);
  Screen('Hookfunction', win, posstring, procchain, 'OculusVRClientCompositingShaderLeftEye', glsl, blittercfg);

  % Setup right eye shader:
  glsl = LoadGLSLProgramFromFiles([fileparts(mfilename('fullpath')) filesep 'OculusRiftCorrectionShader']);
  glUseProgram(glsl);
  if ~strcmpi(hmd{handle}.basicTask, 'Monoscopic')
    % Stereoscopic display: Source from right eye buffer:
    glUniform1i(glGetUniformLocation(glsl, 'Image'), 1);
  else
    % Monoscopic display: Source right eye image also from left-eye (aka mono) buffer:
    glUniform1i(glGetUniformLocation(glsl, 'Image'), 0);
  end

  glUniform1i(glGetUniformLocation(glsl, 'PrevImage'), 2);
  glUniform1i(glGetUniformLocation(glsl, 'OverdriveLUT'), 3);
  glUniform3f(glGetUniformLocation(glsl, 'OverdriveScales'), overdriveUpScale, overdriveDownScale, overdriveGammaCorrect);
  glUniform2f(glGetUniformLocation(glsl, 'EyeToSourceUVOffset'), hmd{handle}.uvOffsetRight(1) * hmd{handle}.inputWidth, hmd{handle}.uvOffsetRight(2) * hmd{handle}.inputHeight);
  glUniform2f(glGetUniformLocation(glsl, 'EyeToSourceUVScale'), hmd{handle}.uvScaleRight(1) * hmd{handle}.inputWidth, hmd{handle}.uvScaleRight(2) * hmd{handle}.inputHeight);
  glUniformMatrix4fv(glGetUniformLocation(glsl, 'EyeRotationStart'), 1, 1, hmd{handle}.eyeRotStartMatrixRight);
  glUniformMatrix4fv(glGetUniformLocation(glsl, 'EyeRotationEnd'), 1, 1, hmd{handle}.eyeRotEndMatrixRight);
  hmd{handle}.shaderRight = [glsl, glGetUniformLocation(glsl, 'EyeRotationStart'), glGetUniformLocation(glsl, 'EyeRotationEnd')];
  glUseProgram(0);

  % Insert it at former position of the old shader:
  posstring = sprintf('InsertAt%iShader', slot);
  blittercfg = sprintf('Blitter:DisplayListBlit:Handle:%i:Bilinear', gldRight);
  Screen('Hookfunction', win, posstring, procchain, 'OculusVRClientCompositingShaderRightEye', glsl, blittercfg);

  if hmd{handle}.useOverdrive
    if strcmpi(hmd{handle}.basicTask, 'Monoscopic')
      % Need a bufferflip command:
      Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:FlipFBOs', '');
    end

    [realw, realh] = Screen('Windowsize', win, 1);
    Screen('HookFunction', win, 'AppendBuiltin', 'FinalOutputFormattingBlit', 'Builtin:IdentityBlit', sprintf('Blitter:IdentityBlit:OvrSize:%i:%i', realw, realh));
    Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');

    woverdrive1 = Screen('OpenOffscreenwindow', win, 0, [0, 0, realw * 2, realh], [], 32);
    hmd{handle}.overdriveTex(1) = Screen('GetOpenGLTexture', woverdrive1, woverdrive1);
    woverdrive2 = Screen('OpenOffscreenwindow', win, 0, [0, 0, realw * 2, realh], [], 32);
    hmd{handle}.overdriveTex(2) = Screen('GetOpenGLTexture', woverdrive2, woverdrive2);
    hmd{handle}.lastOverdriveTex = 0;

    % Load precomputed overdrive lut into variable 'lut', then build a overdriveLut
    % texture out of it:
    load([fileparts(mfilename('fullpath')) filesep 'RiftDK2lut1.mat']);
    luttex = Screen('MakeTexture', win, lut, [], 32, [], 2);
    hmd{handle}.overdriveLut(1) = Screen('GetOpenGLTexture', win, luttex);

    load([fileparts(mfilename('fullpath')) filesep 'RiftDK2lut2.mat']);
    luttex = Screen('MakeTexture', win, lut, [], 32, [], 2);
    hmd{handle}.overdriveLut(2) = Screen('GetOpenGLTexture', win, luttex);
  end

  % TimeWarp or panel overdrive in use?
  if hmd{handle}.useTimeWarp || hmd{handle}.useOverdrive
    % Need to call the PsychOculusVR(1) callback to do needed setup work:
    posstring = sprintf('InsertAt%iMFunction', slot);
    cmdString = sprintf('PsychOculusVR(1, %i);', handle);
    Screen('Hookfunction', win, posstring, procchain, 'OculusVRTimeWarpSetup', cmdString);
    hmd{handle}.readyForWarp = 1;
  end

  % Need to call the PsychOculusVR(2) callback to do needed finalizer work:
  cmdString = sprintf('PsychOculusVR(2, %i);', handle);
  Screen('Hookfunction', win, 'AppendMFunction', 'LeftFinalizerBlitChain', 'OculusVRLatencyTesterSetup', cmdString);
  Screen('Hookfunction', win, 'Enable', 'LeftFinalizerBlitChain');

  % Need to call the end frame marker function of the Oculus runtime:
  cmdString = sprintf('PsychOculusVRCore(''EndFrameTiming'', %i);', handle);
  Screen('Hookfunction', win, 'PrependMFunction', 'ScreenFlipImpliedOperations', 'OculusVRPostPresentCallback', cmdString);
  Screen('Hookfunction', win, 'Enable', 'ScreenFlipImpliedOperations');

  % Does usercode request auto-closing the HMD or driver when the onscreen window is closed?
  if hmd{handle}.autoclose > 0
    % Attach a window close callback for Device teardown at window close time:
    if hmd{handle}.autoclose == 2
      % Shutdown driver completely:
      Screen('Hookfunction', win, 'AppendMFunction', 'CloseOnscreenWindowPostGLShutdown', 'Shutdown window callback into PsychOculusVR driver.', 'PsychOculusVR(''Close'');');
    else
      % Only close this HMD:
      Screen('Hookfunction', win, 'AppendMFunction', 'CloseOnscreenWindowPostGLShutdown', 'Shutdown window callback into PsychOculusVR driver.', sprintf('PsychOculusVR(''Close'', %i);', handle));
    end

    Screen('HookFunction', win, 'Enable', 'CloseOnscreenWindowPostGLShutdown');
  end

  % Need HSW display?
  if (hmd{handle}.hswdismiss >= 0) && isempty(getenv('PSYCH_OCULUS_HSWSKIP'))
    if IsWin
      % Windows doesn't distinguish keyboards, so don't query "all"
      % keyboards. This gives the same effect as on Linux/OSX, but avoids
      % use of PsychHID and potential "libusb not installed" warnings:
      kbddev = [];
    else
      kbddev = -1;
    end

    if bitand(hmd{myhmd.handle}.hswdismiss, 1)
      KbReleaseWait(kbddev);
    end

    dismiss = 0;
    if PsychOculusVRCore('GetHSWState', handle)
      % Yes: Display HSW text:
      hswtext = ['HEALTH & SAFETY WARNING\n\n' ...
                'Read and follow all warnings\n' ...
                'and instructions included with\n' ...
                'the Headset before use. Headset\n' ...
                'should be calibrated for each user.\n' ...
                'Not for use by children under 13.\n' ...
                'Stop use if you experience any\n' ...
                'discomfort or health reactions.\n\n' ...
                'More: www.oculus.com/warnings\n\n' ...
                'To acknowledge:\n'];

      if bitand(hmd{myhmd.handle}.hswdismiss, 1)
        hswtext = [hswtext 'Press a key\n'];
      end

      if bitand(hmd{myhmd.handle}.hswdismiss, 2)
        hswtext = [hswtext 'Click any mouse button\n'];
      end

      if bitand(hmd{myhmd.handle}.hswdismiss, 4)
        hswtext = [hswtext 'Slightly tap the headset'];
      end

      Screen('Flip', win);
      Screen('Flip', win);

      oldTextSize = Screen('TextSize', win, 16);
      Screen('SelectStereoDrawBuffer', win, 1);
      DrawFormattedText(win, hswtext, 'center', 'center', [0 255 0]);
      Screen('SelectStereoDrawBuffer', win, 0);
      DrawFormattedText(win, hswtext, 'center', 'center', [0 255 0]);
      Screen('TextSize', win, oldTextSize);
      Screen('Flip', win, [], 1);

      % Allow dismiss via tap to the HMD?
      if bitand(hmd{myhmd.handle}.hswdismiss, 4)
        % Enable tracking so we can allow user to dismiss HSW via a
        % slight tap to the HMD - accelerometers will do their thing:
        PsychOculusVRCore('Start', handle);
      end

      % Wait for dismiss via keypress, mouse button click or HMD tap:
      while PsychOculusVRCore('GetHSWState', handle, dismiss)
        % Allow dismiss via keypress?
        if bitand(hmd{myhmd.handle}.hswdismiss, 1) && KbCheck(kbddev)
          dismiss = 1;
        end

        % Allow dismiss via mouse click?
        if bitand(hmd{myhmd.handle}.hswdismiss, 2)
          [dummy1, dummy2, buttons] = GetMouse; %#ok<ASGLU>
          if any(buttons)
            dismiss = 1;
            while any(buttons)
              [dummy1, dummy2, buttons] = GetMouse; %#ok<ASGLU>
            end
          end
        end

        % Need to idle flip here to drive timewarp rendering in
        % case some stuff is enabled:
        Screen('Flip', win, [], 1);
      end

      if bitand(hmd{myhmd.handle}.hswdismiss, 1)
        KbReleaseWait(kbddev);
      end

      if bitand(hmd{myhmd.handle}.hswdismiss, 4)
        % Stop tracking for tap detection:
        PsychOculusVRCore('Stop', handle);
        WaitSecs(1);
      end

      % Clear HSW text:
      Screen('Flip', win);
    end
  end

  if ~isempty(strfind(hmd{myhmd.handle}.basicTask, 'Tracked3DVR'))
    % 3D head tracked VR rendering task: Start tracking as a convenience:
    PsychOculusVRCore('Start', handle);
  end

  % Return success result code 1:
  varargout{1} = 1;
  return;
end

% 'cmd' so far not dispatched? Let's assume it is a command
% meant for PsychOculusVRCore:
if (length(varargin) >= 1) && isstruct(varargin{1})
  myhmd = varargin{1};
  handle = myhmd.handle;
  [ varargout{1:nargout} ] = PsychOculusVRCore(cmd, handle, varargin{2:end});
else
  [ varargout{1:nargout} ] = PsychOculusVRCore(cmd, varargin{:});
end

return;

end
