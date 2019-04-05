function varargout = PsychOpenHMDVR(cmd, varargin)
% PsychOpenHMDVR - A high level driver for VR hardware supported via OpenHMD.
%
% Note: If you want to write VR code that is portable across
% VR headsets of different vendors, then use the PsychVRHMD()
% driver instead of this driver. The PsychVRHMD driver will use
% this driver as appropriate when connecting to a OpenHMD supported
% device, but it will also automatically work with other head mounted
% displays. This driver may expose a few functions specific to OpenHMD,
% so you can mix calls to this driver with calls to PsychVRHMD.
%
% Setup instructions:
%
% This driver needs libopenhmd.so version 0.3 or later to be installed
% in a linker accessible path (e.g., /usr/local/lib/ on a Linux system).
% You can either download, compile and install it from ...
%
% https://github.com/OpenHMD/OpenHMD
%
% ... or get a precompiled library for libopenhmd.so from:
%
% RaspberryPi/Raspbian: https://github.com/Psychtoolbox-3/MiscStuff/tree/master/OpenHMD32BitRaspbianARMv7
%
% 64-Bit Intel/Ubuntu:  https://github.com/Psychtoolbox-3/MiscStuff/tree/master/OpenHMD64BitIntelUbuntuLinux
%
% Follow instructions in the accompanying Readme.txt files.
%
% libopenhmd.so in turn needs libhidapi-libusb.so to be installed in
% a similar path. On Debian GNU/Linux based systems you can install HIDAPI
% via the package libhidapi-libusb0 (apt-get install libhidapi-libusb0).
%
% From the same URL above, you need to get openhmdkeepalivedaemon, an executable
% file, and make sure it gets started during system boot of your machine. This so
% the HMD gets correctly detected as video output by the X-Server and by Psychtoolbox,
% otherwise stimuli may not display on the HMD, but on your regular display. This
% executable is not needed for the Oculus Rift DK1 or DK2.
%
% LIMITATIONS: This driver is currently considered BETA quality and may
% undergo backwards incompatible changes without prior warning or notice.
% Use at your own risk!
%
% Usage:
%
% oldverbosity = PsychOpenHMDVR('Verbosity' [, newverbosity]);
% - Get/Set level of verbosity for driver status messages, warning messages,
% error messages etc. 'newverbosity' is the optional new verbosity level,
% 'oldverbosity' is the currently set verbosity level - ie. before changing
% it.  Valid settings are: 0 = Silent, 1 = Errors only, 2 = Warnings, 3 = Info,
% 4 = Debug.
%
%
% hmd = PsychOpenHMDVR('AutoSetupHMD' [, basicTask='Tracked3DVR'][, basicRequirements][, basicQuality=0][, deviceIndex]);
% - Open a OpenHMD supported HMD, set it up with good default rendering and
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
% the HMD is just used as a special kind of standard display monitor. Please
% note that the Oculus Rift DK2 has a special video mode for such monoscopic
% presentation: If you set the video mode to 1080x948@120 Hz, then the DK2
% will only display a monoscopic image identical to both eyes, but at a refresh
% rate of 120 Hz. This allows you to (ab)use the DK2 as a high-speed 120 Hz
% monitor!
%
% 'basicRequirements' defines basic requirements for the task. Currently
% defined are the following strings which can be combined into a single
% 'basicRequirements' string: 'LowPersistence' = Try to keep exposure
% time of visual images on the retina low if possible, ie., try to approximate
% a pulse-type display instead of a hold-type display if possible. This has
% no effect at the moment for this driver and its supported devices.
%
% 'PerEyeFOV' = Request use of per eye individual and asymmetric fields of view even
% when the 'basicTask' was selected to be 'Monoscopic' or 'Stereoscopic'. This allows
% for wider field of view in these tasks, but requires the usercode to adapt to these
% different and asymmetric fields of view for each eye, e.g., by selecting proper 3D
% projection matrices for each eye.
%
% 'FastResponse' = Try to switch images with minimal delay and fast
% pixel switching time. This does nothing on this driver at the moment.
%
% 'TimingSupport' = Support some hardware specific means of timestamping
% or latency measurements. This does nothing on this driver at the moment.
%
% 'TimeWarp' = Enable per eye image 2D timewarping via prediction of eye
% poses at scanout time. This mostly only makes sense for head-tracked 3D
% rendering. Depending on 'basicQuality' a more cheap or more expensive
% procedure is used. This does nothing on this driver at the moment.
%
% 'basicQuality' defines the basic tradeoff between quality and required
% computational power. A setting of 0 gives lowest quality, but with the
% lowest performance requirements. A setting of 1 gives maximum quality at
% maximum computational load. Values between 0 and 1 change the quality to
% performance tradeoff.
%
%
% hmd = PsychOpenHMDVR('Open' [, deviceIndex], ...);
% - Open HMD with index 'deviceIndex'. See PsychOpenHMDVRCore Open?
% for help on additional parameters.
%
%
% PsychOpenHMDVR('SetAutoClose', hmd, mode);
% - Set autoclose mode for HMD with handle 'hmd'. 'mode' can be
% 0 (this is the default) to not do anything special. 1 will close
% the HMD 'hmd' when the onscreen window is closed which displays
% on the HMD. 2 will do the same as 1, but close all open HMDs and
% shutdown the complete driver and OpenHMD runtime - a full cleanup.
%
%
% isOpen = PsychOpenHMDVR('IsOpen', hmd);
% - Returns 1 if 'hmd' corresponds to an open HMD, 0 otherwise.
%
%
% PsychOpenHMDVR('Close' [, hmd])
% - Close provided HMD device 'hmd'. If no 'hmd' handle is provided,
% all HMDs will be closed and the driver will be shutdown.
%
%
% PsychOpenHMDVR('Controllers', hmd);
% - Return a bitmask of all connected controllers: Can be the bitand
% of the OVR.ControllerType_XXX flags described in 'GetInputState'.
% This does not detect if controllers are hot-plugged or unplugged after
% the HMD was opened. Iow. only probed at 'Open'.
% As the current OpenHMD driver does not support dedicated controllers at the
% moment, this always returns 0.
%
%
% info = PsychOpenHMDVR('GetInfo', hmd);
% - Retrieve a struct 'info' with information about the HMD 'hmd'.
% The returned info struct contains at least the following standardized
% fields with information:
% handle = Driver internal handle for the specific HMD.
% driver = Function handle to the actual driver for the HMD, e.g., @PsychOpenHMDVR.
% type   = Defines the type/vendor of the device, e.g., 'OpenHMD'.
% modelName = Name string with the name of the model of the device, e.g., 'Rift (DK2)'.
% separateEyePosesSupported = 1 if use of PsychOpenHMDVR('GetEyePose') will improve
%                             the quality of the VR experience, 0 if no improvement
%                             is to be expected, so 'GetEyePose' can be avoided
%                             to save processing time without a loss of quality.
%                             This always returns 0 on this driver.
%
% The returned struct may contain more information, but the fields mentioned
% above are the only ones guaranteed to be available over the long run. Other
% fields may disappear or change their format and meaning anytime without
% warning.
%
%
% isSupported = PsychOpenHMDVRCore('Supported');
% - Returns 1 if the OpenHMD driver is functional, 0 otherwise. The
% driver is functional if the OpenHMD runtime library was successfully
% initialized. It would return 0 if the required runtime library would
% not be correctly installed.
%
%
% [isVisible, playAreaBounds, OuterAreaBounds] = PsychOpenHMDVRCore('VRAreaBoundary', hmd [, requestVisible]);
% - Request visualization of the VR play area boundary for 'hmd' and returns its
% current extents.
%
% As VR area boundaries are not actually supported by this OpenHMD driver,
% this function returns no-op results, compatible with what other drivers
% would return if their guardian system would not be set up.
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
% input = PsychOpenHMDVRCore('GetInputState', hmd, controllerType);
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
% pulseEndTime = PsychOpenHMDVR('HapticPulse', hmd, controllerType [, duration=2.5][, freq=1.0][, amplitude=1.0]);
% - Fake triggering a haptic feedback pulse. This does nothing, but return a made up
% but consistent 'pulseEndTime', as this OpenHMD driver currently does not support
% haptic feedback.
%
%
% state = PsychOpenHMDVRCore('PrepareRender', hmd [, userTransformMatrix][, reqmask=1][, targetTime]);
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
% position tracker like, e.g., the OpenHMD Rift DK2 camera. We also return a +128
% flag which means the HMD is actually strapped onto the subjects head and displaying
% our visual content. We can't detect actual HMD display state, but do this for
% compatibility to other drivers.
%
% state also always contains a field state.SessionState, whose bits signal general
% VR session status. In our case we always return +7 on this OpenHMD driver, as we
% can't detect ShouldQuit, ShouldRecenter or DisplayLost conditions, neither if the
% HMD is strapped to the users head.
%
% +1  = Our rendering goes to the HMD, ie. we have control over it. Lack of this could
%       mean the Health and Safety warning is displaying at the moment and waiting for
%       acknowledgement, or some other application is in control.
% +2  = HMD is present and active.
% +4  = HMD is strapped onto users head.
% +8  = DisplayLost condition! Some hardware/software malfunction, need to completely
%       quit this Psychtoolbox session to recover from this.
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
%      and orientations of left and right hand tracking controllers, if any. As this OpenHMD
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
% eyePose = PsychOpenHMDVR('GetEyePose', hmd, renderPass [, userTransformMatrix][, targetTime]);
% - Return a struct 'eyePose' which contains various useful bits of information
% for 3D stereoscopic rendering of the stereo view of one eye, based on head
% tracking data. This function provides essentially the same information as
% the 'PrepareRender' function, but only for one eye. Therefore you will need
% to call this function twice, once for each of the two renderpasses, at the
% beginning of each renderpass. Note: Currently there is no advantage to using
% this function on top of 'PrepareRender', it only increases overhead and is here
% only for compatibility with other drivers.
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
% PsychOpenHMDVR('SetupRenderingParameters', hmd [, basicTask='Tracked3DVR'][, basicRequirements][, basicQuality=0][, fov=[HMDRecommended]][, pixelsPerDisplay=1])
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
% PsychOpenHMDVR('SetBasicQuality', hmd, basicQuality);
% - Set basic level of quality vs. required GPU performance.
%
%
% oldSetting = PsychOpenHMDVR('SetFastResponse', hmd [, enable]);
% - Return old setting for 'FastResponse' mode in 'oldSetting',
% optionally disable or enable the mode via specifying the 'enable'
% parameter as 0 or greater than zero.
%
% Currently not implemented / supported. Does nothing.
%
%
% oldSetting = PsychOpenHMDVR('SetTimeWarp', hmd [, enable]);
% - Return old setting for 'TimeWarp' mode in 'oldSetting', which is
% always 0 for 'disabled'.
% Deprecated: Does nothing. Only for PsychVRHMD compatibility.
%
%
% oldSetting = PsychOpenHMDVR('SetLowPersistence', hmd [, enable]);
% - Return old setting for 'LowPersistence' mode in 'oldSetting',
% optionally enable or disable the mode via specifying the 'enable'
% parameter as 1 or 0.
% Currently not implemented / supported. Does nothing.
%
%
% PsychOpenHMDVR('SetHSWDisplayDismiss', hmd [, dismissTypes=1+2]);
% - Set how the user can dismiss the "Health and safety warning display".
% 'dismissTypes' can be -1 to disable the HSWD, or a value >= 0 to show
% the HSWD until a timeout and or until the user dismisses the HSWD.
% The following flags can be added to define type of dismissal:
%
% +0 = Display until timeout, if any. Will wait forever if there isn't any timeout!
% +1 = Dismiss via keyboard keypress.
% +2 = Dismiss via mouse click or mousepad tap.
%
%
% [bufferSize, imagingFlags, stereoMode] = PsychOpenHMDVR('GetClientRenderingParameters', hmd);
% - Retrieve recommended size in pixels 'bufferSize' = [width, height] of the client
% renderbuffer for each eye for rendering to the HMD. Returns parameters
% previously computed by PsychOpenHMDVR('SetupRenderingParameters', hmd).
%
% Also returns 'imagingFlags', the required imaging mode flags for setup of
% the Screen imaging pipeline. Also returns the needed 'stereoMode' for the
% pipeline.
%
%
% needPanelFitter = PsychOpenHMDVR('GetPanelFitterParameters', hmd);
% - 'needPanelFitter' is 1 if a custom panel fitter tasks is needed, and 'bufferSize'
% from the PsychVRHMD('GetClientRenderingParameters', hmd); defines the size of the
% clientRect for the onscreen window. 'needPanelFitter' is 0 if no panel fitter is
% needed.
%
%
% [winRect, ovrfbOverrideRect, ovrSpecialFlags] = PsychOpenHMDVR('OpenWindowSetup', hmd, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags);
% - Compute special override parameters for given input/output arguments, as needed
% for a specific HMD. Take other preparatory steps as needed, immediately before the
% Screen('OpenWindow') command executes. This is called as part of PsychImaging('OpenWindow'),
% with the user provided hmd, screenid, winRect etc.
%
%
% isOutput = PsychOpenHMDVR('IsHMDOutput', hmd, scanout);
% - Returns 1 (true) if 'scanout' describes the video output to which the
% HMD 'hmd' is connected. 'scanout' is a struct returned by the Screen
% function Screen('ConfigureDisplay', 'Scanout', screenid, outputid);
% This allows probing video outputs to find the one which feeds the HMD.
%
%
% [headToEyeShiftv, headToEyeShiftMatrix] = PsychOpenHMDVR('GetEyeShiftVector', hmd, eye);
% - Retrieve 3D translation vector [tx, ty, tz] that defines the 3D position of the given
% eye 'eye' for the given HMD 'hmd', relative to the origin of the local head/HMD
% reference frame. This is needed to translate a global head pose into a eye
% pose, e.g., to translate the output of PsychOpenHMDVR('GetEyePose') into actual
% tracked/predicted eye locations for stereo rendering.
%
% In addition to the 'headToEyeShiftv' vector, a corresponding 4x4 translation
% matrix is also returned in 'headToEyeShiftMatrix' for convenience.
%
%

% History:
% 30-Jul-2017  mk   Written. derived from original classic PsychOculusVR driver.

% Global GL handle for access to OpenGL constants needed in setup:
global GL;
global OVR;
global hsw_dismissed;

persistent hmd;

if nargin < 1 || isempty(cmd)
  help PsychOpenHMDVR;
  fprintf('\n\nAlso available are functions from PsychOpenHMDVRCore:\n');
  PsychOpenHMDVRCore;
  return;
end

% Fast-Path callback function for VR compositing path into window backbuffer:
if cmd == 0
  handle = varargin{1};

  % Switch to onscreen window system backbuffer:
  glBindFramebuffer(GL.FRAMEBUFFER, 0);

  % Fast-path: Use precompiled display list to execute VR warp blit:
  if ~isempty(hmd{handle}.dl)
    glCallList(hmd{handle}.dl);
    return;
  end

  % Need to create a display list with VR warp blit first:
  imax = 2;

  dl = glGenLists(1);
  glNewList(dl, GL.COMPILE);

  glMatrixMode(GL.MODELVIEW);
  glPushMatrix;
  glLoadIdentity;

  % On 90 degrees rotated "portrait mode" panel setups, like in the
  % Oculus Rift DK2, rotate our rendering and viewport:
  if hmd{handle}.panelRotated
    glViewport(0, 0, hmd{handle}.panelHeight, hmd{handle}.panelWidth);
    glRotatef(90, 0, 0, 1);
  else
    glViewport(0, 0, hmd{handle}.panelWidth, hmd{handle}.panelHeight);
  end

  % Is this a Rift DK2 in special 1080x948 @ 120Hz video mode?
  if strcmp(hmd{handle}.modelName, 'Rift (DK2)') && hmd{handle}.panelRotated && (hmd{handle}.panelWidth == 948)
    % Yes. In this mode, the content of the framebuffer is considered a monoscopic image,
    % and it is replicated into the top and bottom half of the display panel, thereby
    % cloning the framebuffer as identical image into left and right eye of the user.
    % Iow. the Rift acts as a monoscopic image viewer, with 120 Hz refresh rate. It still
    % needs lens distortion correction and color abberation correction though. Ergo still
    % run our warp shader, but only the left-eye one, sourcing from one input texture,
    % doing the warp:
    imax = 1;

    % Extend the blitted quad from that one shader pass to the full framebuffer/viewport
    % size, by some translation and scaling of the quad. The Rift will undo this during
    % "cloning" of the "fullscreen" quad into both eyes:
    glTranslatef(1, 0, 0);
    glScalef(2, 1, 1);
  end

  glMatrixMode(GL.PROJECTION);
  glPushMatrix;
  glLoadIdentity;

  for i=1:imax
    glUseProgram(hmd{handle}.glsl(i));
    glBindTexture(GL.TEXTURE_2D, hmd{handle}.inTex(i));
    glBegin(GL.QUADS);
    if i == 1
      ipd = -hmd{handle}.ipd;
      glTexCoord2d(ipd + 0,  0);
      glVertex3d(  -1, -1, 0);
      glTexCoord2d(ipd + 1,  0);
      glVertex3d(   0, -1, 0);
      glTexCoord2d(ipd + 1,  1);
      glVertex3d(   0,  1, 0);
      glTexCoord2d(ipd + 0,  1);
      glVertex3d(  -1,  1, 0);
    else
      ipd = +hmd{handle}.ipd;
      glTexCoord2d(ipd + 0,  0);
      glVertex3d(   0, -1, 0);
      glTexCoord2d(ipd + 1,  0);
      glVertex3d(   1, -1, 0);
      glTexCoord2d(ipd + 1,  1);
      glVertex3d(   1,  1, 0);
      glTexCoord2d(ipd + 0,  1);
      glVertex3d(   0,  1, 0);
    end
    glEnd();
  end

  glBindTexture(GL.TEXTURE_2D, 0);
  glUseProgram(0);

  if 0
    glPointSize(10);
    glBegin(GL.POINTS);
    glColor3f(1,1,0);
    glVertex2d(-0.5, 0);
    glVertex2d(+0.5, 0);
    glEnd;
    glPointSize(1);
  end

  glPopMatrix();
  glMatrixMode(GL.MODELVIEW);
  glPopMatrix();

  % Done compiling the display list:
  glEndList;
  hmd{handle}.dl = dl;

  return;
end

if strcmpi(cmd, 'PrepareRender')
  % Get and validate handle - fast path open coded:
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOpenHMDVR:PrepareRender: Specified handle does not correspond to an open HMD!');
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

  % Get predicted head pose for targetTime:
  state = PsychOpenHMDVRCore('GetTrackingState', myhmd.handle, targetTime);

  % Always return basic tracking status:
  result.tracked = state.Status;

  % Always return faked session state:
  result.SessionState = 7; % = 1 + 2 + 4 = All is fine, no trouble, subject is attentive ;-)

  % Want matrices which take a usercode supplied global transformation into account?
  if bitand(reqmask, 1)
    % Yes: We need tracked + predicted head pose, so we can apply the user transform,
    % and then per-eye transforms:

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
    % Yes: We can't do this on current OpenHMD, so fake stuff:

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
    error('PsychOpenHMDVR:GetEyePose: Specified handle does not correspond to an open HMD!');
  end

  % Valid: Get view render pass for which to return information:
  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOpenHMDVR:GetEyePose: Required ''renderPass'' argument missing.');
  end
  renderPass = varargin{2};

  % Get 'userTransformMatrix' if any:
  if length(varargin) >= 3 && ~isempty(varargin{3})
    userTransformMatrix = varargin{3};
  else
    % Default: Identity transform to do nothing:
    userTransformMatrix = diag([1 1 1 1]);
  end

  % Get eye pose for this renderPass, or more exactly the headPose from which this
  % renderPass eyePose will get computed:
  [result.headPose, result.eyeIndex] = PsychOpenHMDVRCore('GetEyePose', myhmd.handle, renderPass);

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
  % Check if the OpenHMD VR runtime is supported and active on this
  % installation, so it can be used to open connections to real HMDs,
  % or at least to emulate a HMD for simple debugging purposes:
  try
    if exist('PsychOpenHMDVRCore', 'file') && PsychOpenHMDVRCore('GetCount') >= 0
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
    error('PsychOpenHMDVR:GetInputState: Specified handle does not correspond to an open HMD!');
  end

  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOpenHMDVR:GetInputState: Required ''controllerType'' argument missing.');
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
    error('PsychOpenHMDVR:HapticPulse: Specified handle does not correspond to an open HMD!');
  end

  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOpenHMDVR:HapticPulse: Required ''controllerType'' argument missing.');
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
  if ~PsychOpenHMDVR('IsOpen', myhmd)
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
  if ~PsychOpenHMDVR('Supported')
    % Nope: Game over.
    fprintf('PsychOpenHMDVR:AutoSetupHMD: Could not initialize OpenHMD driver. Game over!\n');

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
    newhmd = PsychOpenHMDVR('Open', deviceIndex);
  else
    % Check if at least one OpenHMD HMD is connected and available:
    if PsychOpenHMDVR('GetCount') > 0
      % Yes. Open and initialize connection to first detected HMD:
      fprintf('PsychOpenHMDVR: Opening the first connected OpenHMD VR headset.\n');
      newhmd = PsychOpenHMDVR('Open', 0);
    else
      % No. Open an emulated/simulated HMD for basic testing and debugging:
      fprintf('PsychOpenHMDVR: No OpenHMD HMD detected. Opening a simulated HMD.\n');
      newhmd = PsychOpenHMDVR('Open', -1);
    end
  end

  % Trigger an automatic device close at onscreen window close for the HMD display window:
  PsychOpenHMDVR('SetAutoClose', newhmd, 1);

  % Setup default rendering parameters:
  PsychOpenHMDVR('SetupRenderingParameters', newhmd, basicTask, basicRequirements, basicQuality);

  % Add a PsychImaging task to use this HMD with the next opened onscreen window:
  PsychImaging('AddTask', 'General', 'UseVRHMD', newhmd);

  % Return the device handle:
  varargout{1} = newhmd;

  % Ready.
  return;
end

if strcmpi(cmd, 'SetAutoClose')
  myhmd = varargin{1};

  if ~PsychOpenHMDVR('IsOpen', myhmd)
    error('PsychOpenHMDVR:SetAutoClose: Specified handle does not correspond to an open HMD!');
  end

  % Assign autoclose flag:
  hmd{myhmd.handle}.autoclose = varargin{2};

  return;
end

if strcmpi(cmd, 'SetHSWDisplayDismiss')
  myhmd = varargin{1};

  if ~PsychOpenHMDVR('IsOpen', myhmd)
    error('PsychOpenHMDVR:SetHSWDisplay: Specified handle does not correspond to an open HMD!');
  end

  % Method of dismissing HSW display:
  if length(varargin) < 2 || isempty(varargin{2})
    % Default is keyboard or mouse click:
    hmd{myhmd.handle}.hswdismiss = 1 + 2;
  else
    hmd{myhmd.handle}.hswdismiss = varargin{2};
  end

  return;
end

% Open a HMD:
if strcmpi(cmd, 'Open')
  [handle, modelName, panelWidth, panelHeight] = PsychOpenHMDVRCore('Open', varargin{:});

  newhmd.handle = handle;
  newhmd.driver = @PsychOpenHMDVR;
  newhmd.type   = 'OpenHMD';
  newhmd.subtype = 'OpenHMD';
  newhmd.open = 1;
  newhmd.modelName = modelName;
  newhmd.panelWidth = panelWidth;
  newhmd.panelHeight = panelHeight;
  newhmd.separateEyePosesSupported = 0;
  newhmd.controllerTypes = 0;
  newhmd.VRControllersSupported = 0;
  newhmd.handTrackingSupported = 0;
  newhmd.hapticFeedbackSupported = 0;

  % Default autoclose flag to "no autoclose":
  newhmd.autoclose = 0;

  % By default allow user to dismiss HSW display via key press
  % or mouse click:
  newhmd.hswdismiss = 1 + 2;

  % Setup basic task/requirement/quality specs to "nothing":
  newhmd.basicQuality = 0;
  newhmd.basicTask = '';
  newhmd.basicRequirements = '';

  % Init to "unknown":
  metersPerTanAngleAtCenter = 0;

  if strcmp(modelName, 'Rift (DK2)')
    % Probably right for DK2:
    metersPerTanAngleAtCenter = 0.0425 * 1.33;
  end

  if strcmp(modelName, 'Rift (CV1)')
    % Maybe right for CV1?
    metersPerTanAngleAtCenter = 0.036 * 1.33;
  end

  if metersPerTanAngleAtCenter == 0
    % Set size of display in meters per radians tan angle at the display center.
    % Default to the setting of a Oculus Rift DK2, so we have something
    % to work with if we can't find the true value for a HMD:
    metersPerTanAngleAtCenter = 0.0425 * 1.33;
    warning('PsychOpenHMDVRCore: Open: metersPerTanAngleAtCenter value unknown for this HMD, using Oculus Rift DK2 default. Please try to find a proper setting!');
  end

  newhmd.metersPerTanAngleAtCenter = metersPerTanAngleAtCenter;

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
  if ~PsychOpenHMDVR('IsOpen', myhmd)
    error('Controllers: Passed in handle does not refer to a valid and open HMD.');
  end

  varargout{1} = myhmd.controllerTypes;
  return;
end

if strcmpi(cmd, 'GetInfo')
  % Ok, cheap trick: We just return the passed in 'hmd' struct - the up to date
  % internal copy that is:
  if ~PsychOpenHMDVR('IsOpen', varargin{1})
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
      PsychOpenHMDVRCore('Close', myhmd.handle);
      hmd{myhmd.handle}.open = 0;
    end
  else
    % Shutdown whole driver:
    PsychOpenHMDVRCore('Close');
    hmd = [];
  end

  return;
end

if strcmpi(cmd, 'IsHMDOutput')
  myhmd = varargin{1}; %#ok<NASGU>
  scanout = varargin{2};

  % Is this an output with a resolution matching HMD panel resolution?
  % Assumption here is that it is a tilted panel in portrait mode in case of
  % the Rift DK1/DK2, but a non-tilted panel in landscape mode on other HMDs,
  % e.g., the Rift CV1:
  if (~isempty(strfind(myhmd.modelName, '(DK')) && (scanout.width == myhmd.panelHeight) && (scanout.height == myhmd.panelWidth || scanout.height == 948)) || ...
     (isempty(strfind(myhmd.modelName, '(DK')) && (scanout.width == myhmd.panelWidth) && (scanout.height == myhmd.panelHeight))
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

  if ~isempty(strfind(hmd{handle}.basicRequirements, 'LowPersistence'))
    PsychOpenHMDVRCore('SetLowPersistence', handle, 1);
  else
    PsychOpenHMDVRCore('SetLowPersistence', handle, 0);
  end

  % Dynamic prediction enables advanced head tracking prediction and eye timewarping:
  if ~isempty(strfind(hmd{handle}.basicRequirements, 'TimingSupport')) || ...
     ~isempty(strfind(hmd{handle}.basicTask, 'Tracked3DVR'))
    %PsychOpenHMDVRCore('SetDynamicPrediction', handle, 1);
  else
    %PsychOpenHMDVRCore('SetDynamicPrediction', handle, 0);
  end

  return;
end

if strcmpi(cmd, 'SetFastResponse')
  myhmd = varargin{1};
  if ~PsychOpenHMDVR('IsOpen', myhmd)
    error('SetFastResponse: Passed in handle does not refer to a valid and open HMD.');
  end
  handle = myhmd.handle;

  % Does nothing, always return "FastResponse off":
  varargout{1} = 0;
  return;
end

if strcmpi(cmd, 'SetTimeWarp')
  myhmd = varargin{1};
  if ~PsychOpenHMDVR('IsOpen', myhmd)
    error('SetTimeWarp: Passed in handle does not refer to a valid and open HMD.');
  end

  % Does nothing. Always returns "TimeWarp off":
  varargout{1} = 0;

  return;
end

if strcmpi(cmd, 'SetLowPersistence')
  myhmd = varargin{1};
  if ~PsychOpenHMDVR('IsOpen', myhmd)
    error('SetLowPersistence: Passed in handle does not refer to a valid and open HMD.');
  end

  % SetLowPersistence determines use low persistence mode on the Rift DK2. Return old setting:
  varargout{1} = PsychOpenHMDVRCore('SetLowPersistence', myhmd.handle);

  % New setting requested?
  if (length(varargin) >= 2) && ~isempty(varargin{2})
    PsychOpenHMDVRCore('SetLowPersistence', myhmd.handle, varargin{2});
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

  if length(varargin) >= 5 && ~isempty(varargin{5})
    fov = varargin{5};
  else
    fov = [];
  end

  if length(varargin) >= 6 && ~isempty(varargin{6})
    pixelsPerDisplay = varargin{6};
  else
    pixelsPerDisplay = 1;
  end

  hmd{myhmd.handle}.basicTask = basicTask;
  hmd{myhmd.handle}.basicRequirements = basicRequirements;

  PsychOpenHMDVR('SetBasicQuality', myhmd, basicQuality);

  metersPerTanAngleAtCenter = myhmd.metersPerTanAngleAtCenter * pixelsPerDisplay;

  % Get optimal client renderbuffer size - the size of our virtual framebuffer for left eye:
  [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovL] = PsychOpenHMDVRCore('GetFovTextureSize', myhmd.handle, 0, metersPerTanAngleAtCenter, fov, varargin{7:end});

  % Get optimal client renderbuffer size - the size of our virtual framebuffer for right eye:
  [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovR] = PsychOpenHMDVRCore('GetFovTextureSize', myhmd.handle, 1, metersPerTanAngleAtCenter, fov, varargin{7:end});

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
    [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovL] = PsychOpenHMDVRCore('GetFovTextureSize', myhmd.handle, 0, metersPerTanAngleAtCenter, fov, varargin{7:end});
    [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovR] = PsychOpenHMDVRCore('GetFovTextureSize', myhmd.handle, 1, metersPerTanAngleAtCenter, fov, varargin{7:end});
  end

  return;
end

if strcmpi(cmd, 'GetClientRenderingParameters')
  myhmd = varargin{1};
  varargout{1} = [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight];

  % We need fast backing store support for imaging pipeline with virtual framebuffers,
  % also output redirection for providing final output to us:
  imagingMode = mor(kPsychNeedFastBackingStore, kPsychNeedFinalizedFBOSinks);

  % Does usercode want a 16 bpc half-float rendering pipeline?
  if ~isempty(strfind(hmd{myhmd.handle}.basicRequirements, 'Float16Display'))
    % Request a 16 bpc float framebuffer from Psychtoolbox:
    imagingMode = mor(imagingMode, kPsychNeed16BPCFloat);
  else
    % Standard RGBA8 images: Use sRGB format for rendering/blending/compositing/display:
    % TODO imagingMode = mor(imagingMode, kPsychEnableSRGBRendering);
  end

  if ~strcmpi(hmd{myhmd.handle}.basicTask, 'Monoscopic')
    % We must use stereomode 12, so we get separate draw buffers for left and
    % right eye, and separate stream processing into our backing textures, with
    % all internal buffers at at least full eye resolution.
    stereoMode = 12;
  else
    % Monoscopic presentation will do:
    stereoMode = 0;
  end

  varargout{2} = imagingMode;
  varargout{3} = stereoMode;
  return;
end

if strcmpi(cmd, 'GetPanelFitterParameters')
  % No need for PanelFitter:
  varargout{1} = 0;
  return;
end

% [winRect, ovrfbOverrideRect, ovrSpecialFlags] = PsychOpenHMDVR('OpenWindowSetup', hmd, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags);
if strcmpi(cmd, 'OpenWindowSetup')
  myhmd = varargin{1};
  screenid = varargin{2};
  winRect = varargin{3};
  ovrfbOverrideRect = varargin{4};
  ovrSpecialFlags = varargin{5};

  % Override winRect for the OpenHMD dummy HMD device:
  if strcmp(myhmd.modelName, 'Dummy Device') || strcmp(myhmd.modelName, 'External Device')
    winRect = [0, 0, myhmd.panelWidth, myhmd.panelHeight];
  end

  % Yes. Trying to display on a screen with more than one video output?
  if isempty(winRect) && (Screen('ConfigureDisplay', 'NumberOutputs', screenid) > 1)
    % Yes. Not good, as this will impair graphics performance and timing a lot.
    % Warn about this, then try to at least position the onscreen window on the
    % right output.
    fprintf('PsychOpenHMDVR-WARNING: You are requesting display to a VR HMD on a screen with multiple active video outputs.\n');
    fprintf('PsychOpenHMDVR-WARNING: This will impair visual stimulation timing and cause decreased VR performance!\n');
    fprintf('PsychOpenHMDVR-WARNING: I strongly recommend only activating one output on the HMD screen - the HMD output on the screen.\n');
    fprintf('PsychOpenHMDVR-WARNING: On Linux with X11 X-Server, you should create a separate X-Screen for the HMD.\n');

    % Try to find the output with the HMD:
    for i=0:Screen('ConfigureDisplay', 'NumberOutputs', screenid)-1
      scanout = Screen('ConfigureDisplay', 'Scanout', screenid, i);
      if myhmd.driver('IsHMDOutput', myhmd, scanout)
        % This output i has proper resolution to be the HMD panel.
        % Position our onscreen window accordingly:
        winRect = OffsetRect([0, 0, scanout.width, scanout.height], scanout.xStart, scanout.yStart);
        fprintf('PsychOpenHMDVR-Info: Positioning onscreen window at rect [%i, %i, %i, %i] to align with HMD output %i.\n', ...
                winRect(1), winRect(2), winRect(3), winRect(4), i);
      end
    end
  end

  % Get "panel size" / true framebuffer size:
  if isempty(winRect)
    panelRect = Screen('Rect', screenid);
  else
    panelRect = winRect;
  end

  % How is the panel mounted in the HMD, portrait or landscape?
  % Have special case for Rift DK2 as it can also operate in 1080 x 948 mode at 120 Hz.
  if RectWidth(panelRect) < RectHeight(panelRect) || strcmp(myhmd.modelName, 'Rift (DK2)')
    % Portrait: Need to rotate our rendering:
    hmd{myhmd.handle}.panelRotated = 1;
    % Make sure Rift DK2 is handled properly also in 1080 x 948 @120 Hz.
    hmd{myhmd.handle}.panelWidth = RectHeight(panelRect);
    hmd{myhmd.handle}.panelHeight = RectWidth(panelRect);
  else
    % Landscape: No need to rotate:
    hmd{myhmd.handle}.panelRotated = 0;
  end

  % The current design iteration requires the PTB parent onscreen windows
  % effective backbuffer (from the pov of the imaging pipeline) to have the
  % same size (width x height) as the renderbuffer for one eye, so enforce
  % that constraint by setting ovrfbOverrideRect accordingly.

  % Get required output buffer size:
  clientRes = myhmd.driver('GetClientRenderingParameters', myhmd);

  % Set as fbOverrideRect for window:
  ovrfbOverrideRect = [0, 0, clientRes(1), clientRes(2)];

  fprintf('PsychOpenHMDVR-Info: Overriding onscreen window framebuffer size to %i x %i pixels for use with VR-HMD direct output mode.\n', ...
          clientRes(1), clientRes(2));

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
    warning('PTB internal error in PsychOpenHMDVR: GL struct not initialized?!?');
    return;
  end

  % OpenHMD device handle:
  myhmd = varargin{1};
  handle = myhmd.handle;

  % Onscreen window handle:
  win = varargin{2};
  winfo = Screen('GetWindowInfo', win);
  hmd{handle}.StereoMode = winfo.StereoMode;

  % Keep track of window handle of associated onscreen window:
  hmd{handle}.win = win;

  % Init display list handle to none:
  hmd{handle}.dl = [];

  % Need to know user selected clearcolor:
  clearcolor = varargin{3};

  % Also keep track of video refresh duration of the HMD:
  hmd{handle}.videoRefreshDuration = Screen('Framerate', win);
  if hmd{handle}.videoRefreshDuration == 0
    % Unlikely to ever hit this situation, but if we would, just
    % default to a reasonable default video refresh rate of 75 Hz, which
    % is between the lowest known 60 Hz and typical modern 90 Hz.
    hmd{handle}.videoRefreshDuration = 75;
  end
  hmd{handle}.videoRefreshDuration = 1 / hmd{handle}.videoRefreshDuration;

  % Compute effective size of per-eye input buffer for undistortion render.
  % The input buffers for undistortion are the finalizedDrawbufferFBO's, ie.
  % the final output buffers of Psychtoolbox imaging pipeline.
  % This means [inputWidth, inputHeight] == [rbwidth, rbheight].
  hmd{handle}.inputWidth = hmd{handle}.rbwidth;
  hmd{handle}.inputHeight = hmd{handle}.rbheight;

  % Query undistortion parameters for left eye view:
  [hmd{handle}.rbwidth, hmd{handle}.rbheight, hsx, hsy, hsz, abberationk, distortionk, viewport_scale(1), viewport_scale(2)] = PsychOpenHMDVRCore('GetUndistortionParameters', handle, 0, hmd{handle}.inputWidth, hmd{handle}.inputHeight, hmd{handle}.fovL);
  hmd{handle}.HmdToEyeViewOffsetLeft = [hsx, hsy, hsz];

  % Query parameters for right eye view:
  [hmd{handle}.rbwidth, hmd{handle}.rbheight, hsx, hsy, hsz, abberationk, distortionk, viewport_scale(1), viewport_scale(2)] = PsychOpenHMDVRCore('GetUndistortionParameters', handle, 1, hmd{handle}.inputWidth, hmd{handle}.inputHeight, hmd{handle}.fovR);
  hmd{handle}.HmdToEyeViewOffsetRight = [hsx, hsy, hsz];

  % Only have horizontal width per eye:
  viewport_scale(1) = viewport_scale(1) / 2;

  % Convert head to eye shift vectors into 4x4 matrices, as we'll need
  % them frequently:
  EyeT = diag([1 1 1 1]);
  EyeT(1:3, 4) = hmd{handle}.HmdToEyeViewOffsetLeft';
  hmd{handle}.eyeShiftMatrix{1} = EyeT;

  EyeT = diag([1 1 1 1]);
  EyeT(1:3, 4) = hmd{handle}.HmdToEyeViewOffsetRight';
  hmd{handle}.eyeShiftMatrix{2} = EyeT;

  % Retrieve texture handles for the finalizedFBOs, which contain our per-eye input data for VR compositing:
  [hmd{handle}.inTex(1), hmd{handle}.inTex(2), glTextureTarget, format, multiSample, width, height] = Screen('HookFunction', win, 'GetDisplayBufferTextures');

  % In Monoscopic mode, we only have one finalizedFBO, so assign that ones texture
  % for right eye compositing shader as well:
  if winfo.StereoMode == 0
    hmd{handle}.inTex(2) = hmd{handle}.inTex(1);
  end

  % Get GLSL shader source code for the distortion shaders:
  [vertexShaderSrc, fragmentShaderSrc] = PsychOpenHMDVRCore('GetCorrectionShaders', handle);

  left_lens_center(2) = hmd{handle}.HmdToEyeViewOffsetLeft(2);
  right_lens_center(2) = hmd{handle}.HmdToEyeViewOffsetRight(2);
  left_lens_center(1) = viewport_scale(1) + hmd{handle}.HmdToEyeViewOffsetLeft(1);
  right_lens_center(1) = hmd{handle}.HmdToEyeViewOffsetRight(1);

  if left_lens_center(1) > right_lens_center(1)
    warp_scale = left_lens_center(1);
  else
    warp_scale = right_lens_center(1);
  end

  % Setup shaders:
  for i = 1:2
    vertexShader = glCreateShader(GL.VERTEX_SHADER);
    fragmentShader = glCreateShader(GL.FRAGMENT_SHADER);
    glsl(i) = glCreateProgram();

    glAttachShader(glsl(i), vertexShader);
    glAttachShader(glsl(i), fragmentShader);

    glShaderSource(vertexShader, vertexShaderSrc);
    glCompileShader(vertexShader);

    glShaderSource(fragmentShader, fragmentShaderSrc);
    glCompileShader(fragmentShader);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    glLinkProgram(glsl(i));

    glUseProgram(glsl(i));
    glUniform1i(glGetUniformLocation(glsl(i), 'warpTexture'), 0);
    glUniform2fv(glGetUniformLocation(glsl(i), 'ViewportScale'), 1, viewport_scale);
    glUniform3fv(glGetUniformLocation(glsl(i), 'aberr'), 1, abberationk);
    glUniform1f(glGetUniformLocation(glsl(i), 'WarpScale'), warp_scale);
    glUniform4fv(glGetUniformLocation(glsl(i), 'HmdWarpParam'), 1, distortionk);

    if i == 1
      glUniform2fv(glGetUniformLocation(glsl(i), 'LensCenter'), 1, left_lens_center);
    else
      glUniform2fv(glGetUniformLocation(glsl(i), 'LensCenter'), 1, right_lens_center);
    end

    glUseProgram(0);
  end

  % Assign the two shaders:
  hmd{handle}.glsl = glsl;

  % Setup callback from Screen() for triggering the final VR compositing operation:
  cmdString = sprintf('PsychOpenHMDVR(0, %i);', handle);
  if winfo.StereoMode > 0
    % In stereo mode, use right finalizer chain for triggering the VR compositing pass, as it executes last:
    Screen('Hookfunction', win, 'AppendMFunction', 'RightFinalizerBlitChain', 'OpenHMD Stereo commit Operation', cmdString);
    Screen('Hookfunction', win, 'Enable', 'RightFinalizerBlitChain');
  else
    % In mono mode, use left finalizer chain for triggering the VR compositing pass, as it executes only - and therefore last:
    Screen('Hookfunction', win, 'AppendMFunction', 'LeftFinalizerBlitChain', 'OpenHMD Mono commit Operation', cmdString);
    Screen('Hookfunction', win, 'Enable', 'LeftFinalizerBlitChain');
  end

  % Does usercode request auto-closing the HMD or driver when the onscreen window is closed?
  if hmd{handle}.autoclose > 0
    % Attach a window close callback for Device teardown at window close time:
    if hmd{handle}.autoclose == 2
      % Shutdown driver completely:
      Screen('Hookfunction', win, 'AppendMFunction', 'CloseOnscreenWindowPostGLShutdown', 'Shutdown window callback into PsychOpenHMDVR driver.', 'PsychOpenHMDVR(''Close'');');
    else
      % Only close this HMD:
      Screen('Hookfunction', win, 'AppendMFunction', 'CloseOnscreenWindowPostGLShutdown', 'Shutdown window callback into PsychOpenHMDVR driver.', sprintf('PsychOpenHMDVR(''Close'', %i);', handle));
    end

    Screen('HookFunction', win, 'Enable', 'CloseOnscreenWindowPostGLShutdown');
  end

  if ~isempty(strfind(hmd{handle}.modelName, 'Rift (CV'))
    % Attach override projection matrices at least for the Rift CV1:
    [~, ~, ipd] = PsychOpenHMDVRCore('GetStaticRenderParameters', handle);
    hmd{handle}.ipd = ipd / 16 / viewport_scale(1);
  else
    % No ipd correction for other HMDs:
    hmd{handle}.ipd = 0;
  end

  % Need HSW display?
  if (hmd{handle}.hswdismiss >= 0) && isempty(getenv('PSYCH_OpenHMD_HSWSKIP'))
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
    if isempty(hsw_dismissed) || (GetSecs > hsw_dismissed)
      % Set next HSW display to 1 hour from now:
      hsw_dismissed = GetSecs + 3600;

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

      % Wait for dismiss via keypress or mouse button click:
      while ~dismiss
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
      end

      if bitand(hmd{myhmd.handle}.hswdismiss, 1)
        KbReleaseWait(kbddev);
      end

      % Clear HSW text:
      Screen('Flip', win);
    end
  end

  if ~isempty(strfind(hmd{myhmd.handle}.basicTask, 'Tracked3DVR'))
    % 3D head tracked VR rendering task. Start tracking as a convenience:
    PsychOpenHMDVRCore('Start', handle);
  end

  % 3D rendering requested, instead of pure 2D rendering?
  if ~isempty(strfind(hmd{myhmd.handle}.basicTask, '3D'))
    % Yes. Disable ipd correction as it seems to do more harm than good
    % on the only HMD on which it would be used - the Rift CV:
    hmd{handle}.ipd = 0;
  end

  % Is this a Rift DK2 in special 1080x948 @ 120Hz video mode? And no Monoscopic mode selected?
  if strcmp(hmd{handle}.modelName, 'Rift (DK2)') && hmd{handle}.panelRotated && (hmd{handle}.panelWidth == 948) && ...
     isempty(strfind(hmd{myhmd.handle}.basicTask, 'Monoscopic'))
     % Yes. Inform user that this is probably suboptimal, as stereoscopic display
     % won't work in this op-mode of the DK2:
     fprintf('PsychOpenHMDVR: WARNING: Rift DK2 in 1080x948@120 Hz MONO mode detected, but you are not using the Monoscopic\n');
     fprintf('PsychOpenHMDVR: WARNING: basicTask. Note that stereoscopic rendering and display will not work. You only get the\n');
     fprintf('PsychOpenHMDVR: WARNING: left eye view displayed to both eyes, while still wasting resources on stereo rendering.\n');
  end

  % Return success result code 1:
  varargout{1} = 1;
  return;
end

% 'cmd' so far not dispatched? Let's assume it is a command
% meant for PsychOpenHMDVRCore:
if (length(varargin) >= 1) && isstruct(varargin{1})
  myhmd = varargin{1};
  handle = myhmd.handle;
  [ varargout{1:nargout} ] = PsychOpenHMDVRCore(cmd, handle, varargin{2:end});
else
  [ varargout{1:nargout} ] = PsychOpenHMDVRCore(cmd, varargin{:});
end

return;

end
