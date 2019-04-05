function varargout = PsychOculusVR1(cmd, varargin)
% PsychOculusVR1 - A high level driver for Oculus VR hardware using the Version 1.16+ runtime.
%
% Copyright (c) 2018 Mario Kleiner. Licensed under the MIT license.
% The underlying PsychOculusVRCore1 mex driver uses the Oculus SDK, which is
% “Copyright © Facebook Technologies, LLC and its affiliates. All rights reserved.”
% A copy of the Oculus SDK license, its terms of use and thereby redistribution
% terms for the PsychOculusVRCore1 mex file can be found in the PsychtoolboxRoot()
% folder under:
% Psychtoolbox/PsychHardware/PsychVRToolbox/LICENSE_OculusSDK1.txt
%
% Oculus VR's trademarks, e.g., Oculus, Oculus Rift, etc. are registered trademarks
% owned by Oculus VR, LLC.
%
% Note: If you want to write VR code that is portable across
% VR headsets of different vendors, then use the PsychVRHMD()
% driver instead of this driver. The PsychVRHMD driver will use
% this driver as appropriate when connecting to a Oculus Rift
% or similar Oculus device, but it will also automatically work
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
% oldverbosity = PsychOculusVR1('Verbosity' [, newverbosity]);
% - Get/Set level of verbosity for driver status messages, warning messages,
% error messages etc. 'newverbosity' is the optional new verbosity level,
% 'oldverbosity' is the currently set verbosity level - ie. before changing
% it.  Valid settings are: 0 = Silent, 1 = Errors only, 2 = Warnings, 3 = Info,
% 4 = Debug.
%
%
% hmd = PsychOculusVR1('AutoSetupHMD' [, basicTask='Tracked3DVR'][, basicRequirements][, basicQuality=0][, deviceIndex]);
% - Open a Oculus HMD, set it up with good default rendering and
% display parameters and generate a PsychImaging('AddTask', ...)
% line to setup the Psychtoolbox imaging pipeline for proper display
% on the HMD. This will also cause the device connection to get
% auto-closed as soon as the onscreen window which displays on
% the HMD is closed. Returns the 'hmd' handle of the HMD on success.
%
% By default, the first detected HMD will be used and if no VR HMD
% is connected, it will return an empty [] hmd handle. You can override
% this default choice of HMD by specifying the optional 'deviceIndex'
% parameter to choose a specific HMD. However, only one HMD per machine is
% supported, so the 'deviceIndex' will probably be only useful in the future.
%
% More optional parameters: 'basicTask' what kind of task should be implemented.
% The default is 'Tracked3DVR', which means to setup for stereoscopic 3D
% rendering, driven by head motion tracking, for a fully immersive experience
% in some kind of 3D virtual world. This is the default if omitted. The task
% 'Stereoscopic' sets up for display of stereoscopic stimuli, but without
% head tracking. 'Monoscopic' sets up for display of monocular stimuli, ie.
% the HMD is just used as a special kind of standard display monitor. In 'Monoscopic'
% and 'Stereoscopic' mode, both eyes will be presented with an identical field of view,
% to make sure pure 2D drawing works, without the need for setup of special per-eye
% projection transformations. In 'Tracked3DVR' mode, each eye will have a different
% field of view, optimized to maximize the viewable area while still avoiding occlusion
% artifacts due to the nose of the wearer of the HMD.
%
% 'basicRequirements' defines basic requirements for the task. Currently
% defined are the following strings which can be combined into a single
% 'basicRequirements' string:
%
% 'DebugDisplay' = Show the output which is displayed on the HMD inside the
% Psychtoolbox onscreen window as well. This will have a negative impact on
% performance, latency and timing of the HMD visual presentation, so should only
% be used for debugging, as it may cause a seriously degraded VR experience.
% By default, no such debug output is produced and the Psychtoolbox onscreen
% window is not actually displayed on the desktop.
%
% 'Float16Display' = Request rendering, compositing and display in 16 bpc float
% format. This will ask Psychtoolbox to render and post-process stimuli in 16 bpc
% linear floating point format, and allocate 16 bpc half-float textures as final
% renderbuffers to be sent to the VR compositor. If the VR compositor takes advantage
% of the high source image precision is at the discretion of the compositor and HMD.
% By default, if this request is omitted, processing and display in sRGB format is
% requested from Psychtoolbox and the compositor, ie., a roughly gamma 2.2 8 bpc
% format is used, which is optimized for the gamma response curve of at least the Oculus
% Rift CV1 display.
%
% 'PerEyeFOV' = Request use of per eye individual and asymmetric fields of view even
% when the 'basicTask' was selected to be 'Monoscopic' or 'Stereoscopic'. This allows
% for wider field of view in these tasks, but requires the usercode to adapt to these
% different and asymmetric fields of view for each eye, e.g., by selecting proper 3D
% projection matrices for each eye.
%
% 'TimingSupport' = Support some hardware specific means of timestamping
% or latency measurements. On the Rift DK1 this does nothing. On the DK2
% it enables dynamic prediction and timing measurements with the Rifts internal
% latency tester.
%
% 'basicQuality' defines the basic tradeoff between quality and required
% computational power. A setting of 0 gives lowest quality, but with the
% lowest performance requirements. A setting of 1 gives maximum quality at
% maximum computational load. Values between 0 and 1 change the quality to
% performance tradeoff.
%
%
% hmd = PsychOculusVR1('Open' [, deviceIndex], ...);
% - Open HMD with index 'deviceIndex'. See PsychOculusVRCore1 Open?
% for help on additional parameters.
%
%
% PsychOculusVR1('SetAutoClose', hmd, mode);
% - Set autoclose mode for HMD with handle 'hmd'. 'mode' can be
% 0 (this is the default) to not do anything special. 1 will close
% the HMD 'hmd' when the onscreen window is closed which displays
% on the HMD. 2 will do the same as 1, but close all open HMDs and
% shutdown the complete driver and Oculus runtime - a full cleanup.
%
%
% isOpen = PsychOculusVR1('IsOpen', hmd);
% - Returns 1 if 'hmd' corresponds to an open HMD, 0 otherwise.
%
%
% PsychOculusVR1('Close' [, hmd]);
% - Close provided HMD device 'hmd'. If no 'hmd' handle is provided,
% all HMDs will be closed and the driver will be shutdown.
%
%
% PsychOculusVR1('Controllers', hmd);
% - Return a bitmask of all connected controllers: Can be the bitand
% of the OVR.ControllerType_XXX flags described in 'GetInputState'.
% This does not detect if controllers are hot-plugged or unplugged after
% the HMD was opened. Iow. only probed at 'Open'.
%
%
% info = PsychOculusVR1('GetInfo', hmd);
% - Retrieve a struct 'info' with information about the HMD 'hmd'.
% The returned info struct contains at least the following standardized
% fields with information:
%
% handle = Driver internal handle for the specific HMD.
% driver = Function handle to the actual driver for the HMD, e.g., @PsychOculusVR1.
% type   = Defines the type/vendor of the device, e.g., 'Oculus'.
% modelName = Name string with the name of the model of the device, e.g., 'Rift DK2'.
% separateEyePosesSupported = 1 if use of PsychOculusVR1('GetEyePose') will improve
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
% isSupported = PsychOculusVR1('Supported');
% - Returns 1 if the Oculus driver is functional, 0 otherwise. The
% driver is functional if the VR runtime library was successfully
% initialized and a connection to the VR server process has been
% established. It would return 0 if the server process would not be
% running, or if the required runtime library would not be correctly
% installed.
%
%
% [isVisible, playAreaBounds, OuterAreaBounds] = PsychOculusVR1('VRAreaBoundary', hmd [, requestVisible]);
% - Request visualization of the VR play area boundary for 'hmd' and returns its
% current extents.
%
% 'requestVisible' 1 = Request showing the boundary area markers, 0 = Don't
% request showing the markers.
% The driver can not prevent the boundaries to be visualized if some external
% setting asks for their visibility. It can cancel its own request for visibility
% though via 'requestVisible' setting 0.
%
% Returns in 'isVisible' the current visibility status of the VR area boundaries.
%
% 'playAreaBounds' is a 3-by-n matrix defining the play area boundaries. Each
% column represents the [x;y;z] coordinates of one 3D definition point. Connecting
% successive points by line segments defines the boundary, as projected onto the
% floor. Points are listed in clock-wise direction. An empty return argument means
% that the play area is so far undefined.
%
% 'OuterAreaBounds' defines the outer area boundaries in the same way as
% 'playAreaBounds'.
%
%
% [isTriggering, closestDistance, closestPointxyz, surfaceNormal] = PsychOculusVR1('TestVRBoundary', oculusPtr, trackedDeviceType, boundaryType);
% - Return if a tracked device of type 'trackedDeviceType' and associated with 'hmd' is
% colliding with VR area boundaries of 'boundaryType'. This needs the requested devices
% to be online, and the boundaries set up properly, in order to provide meaningful results.
%
% 'trackedDeviceType' is a bit-mask (sum) of the following possible constants:
% OVR.TrackedDevice_HMD = The HMD headset itself,
% OVR.TrackedDevice_LTouch = Left touch controller / left hand.
% OVR.TrackedDevice_RTouch = Right touch controller / right hand.
% OVR.TrackedDevice_Touch = Any/Both of left and right touch controllers.
% OVR.TrackedDevice_Object0 - OVR.TrackedDevice_Object3 = Tracked objects 0 - 3.
% OVR.TrackedDevice_All = All connected tracked devices.
%
% 'boundaryType' 0 = Play area, 1 = Outer boundary.
%
% Return values:
% 'isTriggering' 1 if collision is imminent, 0 otherwise.
% 'closestDistance' Distance to closest point on boundary.
% 'closestPointxyz' [x;y;z] 3D position of closest point on the boundary.
% 'surfaceNormal' [nx;ny;nz] 3D surface normal of closest point.
%
%
% [isTriggering, closestDistance, closestPointxyz, surfaceNormal] = PsychOculusVR1('TestVRBoundaryPoint', oculusPtr, pointxyz, boundaryType);
% - Return if a 3D point 'pointxyz' is colliding with VR area boundaries of 'boundaryType'
% on device 'hmd'. This needs the boundaries set up properly, in order to provide meaningful
% results.
%
% 'pointxyz' = [x,y,z] 3D point position vector.
% 'boundaryType' 0 = Play area, 1 = Outer boundary.
%
% 'isTriggering' 1 if 'pointxyz' is colliding / close, 0 otherwise.
% 'closestDistance' Distance to closest point on boundary.
% 'closestPointxyz' [x;y;z] 3D position of closest point on the boundary.
% 'surfaceNormal' [nx;ny;nz] 3D surface normal of closest point.
%
%
% input = PsychOculusVR1('GetInputState', hmd, controllerType);
% - Get input state of controller 'controllerType' associated with HMD 'hmd'.
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
% 'Touches' Like 'Buttons' but for touch buttons. Use the OVR.Touch_XXX constants to map
% touch points to positions.
%
% 'Trigger'(1/2) = Left (1) and Right (2) trigger: Value range 0.0 - 1.0, filtered and with dead-zone.
% 'TriggerNoDeadzone'(1/2) = Left (1) and Right (2) trigger: Value range 0.0 - 1.0, filtered.
% 'TriggerRaw'(1/2) = Left (1) and Right (2) trigger: Value range 0.0 - 1.0, raw values unfiltered.
% 'Grip'(1/2) = Left (1) and Right (2) grip button: Value range 0.0 - 1.0, filtered and with dead-zone.
% 'GripNoDeadzone'(1/2) = Left (1) and Right (2) grip button: Value range 0.0 - 1.0, filtered.
% 'GripRaw'(1/2) = Left (1) and Right (2) grip button: Value range 0.0 - 1.0, raw values unfiltered.
%
% 'Thumbstick' = 2x2 matrix: Column 1 contains left thumbsticks [x;y] axis values, column 2 contains
%  right sticks [x;y] axis values. Values are in range -1 to +1, filtered and with deadzone applied.
% 'ThumbstickNoDeadzone' = Like 'Thumbstick', filtered, but without a deadzone applied.
% 'ThumbstickRaw' = 'Thumbstick' raw date without deadzone or filtering applied.
%
%
% pulseEndTime = PsychOculusVR1('HapticPulse', hmd, controllerType [, duration=2.5][, freq=1.0][, amplitude=1.0]);
% - Trigger a haptic feedback pulse, some controller vibration, on the specified 'controllerType'
% associated with the specified 'hmd'. 'duration' is pulse duration in seconds, by default a maximum
% of 2.5 seconds is executed. 'freq' is normalized frequency in range 0.0 - 1.0. A value of 0 will
% disable an ongoing pulse. As of Oculus runtime version 1.11, only freq values 0.5 and 1.0 are
% reproduced. Other freq values will be clamped/quantized to those values. 'amplitude' is the
% amplitude of the vibration in normalized 0.0 - 1.0 range.
%
% NOTE: As of November 2018, Oculus SDK 1.16, runtime 1.32, this works very unreliably, and
% often not at all at least with the tested XBox controller!
%
% 'pulseEndTime' returns the expected stop time of vibration in seconds, given the parameters.
% Currently the function will return immediately for a (default) 'duration' of 2.5 seconds, and the pulse
% will end after 2.5 seconds. Smaller 'duration' values will block the execution of the function
% until the 'duration' has passed on some types of controllers.
%
%
% state = PsychOculusVR1('PrepareRender', hmd [, userTransformMatrix][, reqmask=1][, targetTime]);
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
% position tracker like, e.g., the Oculus Rift DK2 or Rift CV1 camera. A +128
% flag means the HMD is actually strapped onto the subjects head and displaying
% our visual content. Lack of this flag means the HMD is off and thereby blanked
% and dark, or we lost access to it to another application.
%
% state also always contains a field state.SessionState, whose bits signal general
% VR session status:
% +1  = Our rendering goes to the HMD, ie. we have control over it. Lack of this could
%       mean the Health and Safety warning is displaying at the moment and waiting for
%       acknowledgement, or the Oculus GUI application is in control.
% +2  = HMD is present and active.
% +4  = HMD is strapped onto users head. A Rift CV1 would switch off/blank if not on the head.
% +8  = DisplayLost condition! Some hardware/software malfunction, need to completely quit this
%       Psychtoolbox session to recover from this.
% +16 = ShouldQuit The user interface / user asks us to voluntarily terminate this session.
% +32 = ShouldRecenter = The user interface asks us to recenter/recalibrate our tracking origin.
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
%      and orientations of left and right Oculus touch controllers, if any.
%
%      state.handStatus(1) = Tracking status of left hand: 0 = Untracked, 1 = Orientation
%                            tracked, 2 = Position tracked, 3 = Orientation and position
%                            tracked. If handStatus is == 0 then all the following information
%                            is invalid and can not be used in any meaningful way.
%      state.handStatus(2) = Tracking status of right hand.
%
%      state.localHandPoseMatrix{1} = 4x4 OpenGL right handed reference frame matrix with
%                                     hand position and orientation encoded to define a
%                                     proper GL_MODELVIEW transform for rendering stuff
%                                     "into"/"relative to" the oriented left hand.
%      state.localHandPoseMatrix{2} = Ditto for the right hand.
%
%      state.globalHandPoseMatrix{1} = userTransformMatrix * state.localHandPoseMatrix{1};
%                                      Left hand pose transformed by passed in userTransformMatrix.
%      state.globalHandPoseMatrix{2} = Ditto for the right hand.
%
%      state.globalHandPoseInverseMatrix{1} = Inverse of globalHandPoseMatrix{1} for collision
%                                             testing/grasping of virtual objects relative to
%                                             hand pose of left hand.
%      state.globalHandPoseInverseMatrix{2} = Ditto for right hand.
%
% More flags to follow...
%
%
% eyePose = PsychOculusVR1('GetEyePose', hmd, renderPass [, userTransformMatrix][, targetTime]);
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
% trackers = PsychOculusVR1('GetTrackersState', hmd);
% - Return a struct array with infos about all connected tracking cameras/sensors
% for 'hmd'.
%
%
% success = PsychOculusVR1('RecenterTrackingOrigin', hmd);
% - Recenter the tracking origin for the 'hmd', based on its current position and
% orientation. Returns 'success' = 1 on success, 0 on failure to recenter. One
% reason for failure could be that the HMD wasn't roughly level, but instead was
% facing upward or downward, which is not allowed. This function also gets automatically
% triggered if the Oculus GUI or other user input requests a recenter. Recenter would then
% be auto-triggered during a call to 'PrepareRender'.
%
%
% oldType = PsychOculusVR1('TrackingOriginType', hmd [, newType]);
% - Specify the type of tracking origin for Oculus device 'oculusPtr'.
% This returns the current type of tracking origin in 'oldType'.
% Optionally you can specify a new tracking origin type as 'newType'.
% Type must be either:
% 0 = Origin is at eye height (HMD height).
% 1 = Origin is at floor height.
% The eye height or floor height gets defined by the system during
% calls to 'RecenterTrackingOrigin' and during sensor calibration in
% the Oculus GUI application.
%
%
% [adaptiveGpuPerformanceScale, frameStats, anyFrameStatsDropped, aswIsAvailable] = PsychOculusVR1('GetPerformanceStats', hmd);
% - Return global and per-frame performance statistics for the given 'hmd'.
%
%
% PsychOculusVR1('SetupRenderingParameters', hmd [, basicTask='Tracked3DVR'][, basicRequirements][, basicQuality=0][, fov=[HMDRecommended]][, pixelsPerDisplay=1])
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
% PsychOculusVR1('SetBasicQuality', hmd, basicQuality);
% - Set basic level of quality vs. required GPU performance.
%
%
% oldValues = PsychOculusVR1('FloatsProperty', hmd, property [, newValues]);
% - Get current values 'oldValues' and optionally set new values 'newValues'
% for a floating point array property with name 'property' on 'hmd'.
%
% 'property' name strings are defined under OVR.KEY_XXX. As of SDK version 1.11,
% the following floats properties exist:
% OVR.KEY_PLAYER_HEIGHT = [Height] of subject in meters.
% OVR.KEY_EYE_HEIGHT = [Height] of eyes of subject above ground in meters.
% OVR.KEY_NECK_TO_EYE_DISTANCE = [horizontal, vertical] neck to eye distance in meters.
% OVR.KEY_EYE_TO_NOSE_DISTANCE = [horizontal, vertical] eye to nose distance in meters.
%
%
% oldString = PsychOculusVR1('StringProperty', hmd, property [, defaultString][, newString]);
% - Get current string 'oldString' and optionally set new string 'newString'
% for a string property with name 'property' on 'hmd'. If the property does not
% exist yet or does not have a string value assigned, optionally return 'defaultString'.
%
% 'property' name strings are defined under OVR.KEY_XXX. As of SDK version 1.11,
% the following string properties exist:
% OVR.KEY_USER = User name.
% OVR.KEY_NAME = Name.
% OVR.KEY_GENDER = Gender 'Male', 'Female', or 'Unknown'.
% OVR.KEY_DEFAULT_GENDER = 'Unknown'.
%
%
% oldSetting = PsychOculusVR1('SetFastResponse', hmd [, enable]);
% - Return old setting for 'FastResponse' mode in 'oldSetting',
% optionally disable or enable the mode via specifying the 'enable'
% parameter as 0 or greater than zero.
%
% Deprecated: This function does nothing. It just exists for (backwards)
% compatibility with the old Oculus driver and PsychVRHMD.
%
%
% oldSetting = PsychOculusVR1('SetTimeWarp', hmd [, enable]);
% - Return old setting for 'TimeWarp' mode in 'oldSetting',
% optionally enable or disable the mode via specifying the 'enable'
% parameter as 1 or 0.
%
% Deprecated: This function does nothing. It just exists for (backwards)
% compatibility with the old Oculus driver and PsychVRHMD.
%
%
% oldSetting = PsychOculusVR1('SetLowPersistence', hmd [, enable]);
% - Return old setting for 'LowPersistence' mode in 'oldSetting',
% optionally enable or disable the mode via specifying the 'enable'
% parameter as 1 or 0.
%
% Deprecated: This function does nothing. It just exists for (backwards)
% compatibility with the old Oculus driver and PsychVRHMD.
%
%
% oldSettings = PsychOculusVR1('PanelOverdriveParameters', hmd [, newparams]);
% Deprecated: This function does nothing. It just exists for (backwards)
% compatibility with the old Oculus driver and PsychVRHMD.
%
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
% PsychOculusVR1('SetHSWDisplayDismiss', hmd [, dismissTypes=1+2+4]);
% - Set how the user can dismiss the "Health and safety warning display".
% Deprecated: This function does nothing. It just exists for (backwards)
% compatibility with the old Oculus driver and PsychVRHMD. The current OculusVR
% runtime 1.0 enforces display of the Health and Safety Warning at the beginning
% of each VR session.
%
%
% PsychOculusVR1('SetHUDState', hmd, mode);
% - Set mode of operation for the builtin head up display (HUD) for performance
% information. The HUD is automatically displayed and updated by the VR compositor
% if enabled, and can be in one of the following selectable 'mode's:
%
% 0 = Head up display HUD off ie. invisible.
% 1 = HUD shows performance summary.
% 2 = HUD shows latency timing.
% 3 = HUD shows application render timing.
% 4 = HUD shows VR compositor render timing.
% 5 = HUD shows Version information of different components.
% 6 = HUD shows Asynchronous time warp (ASW) stats.
% Higher numbers may do something useful in the future.
%
%
% [bufferSize, imagingFlags, stereoMode] = PsychOculusVR1('GetClientRenderingParameters', hmd);
% - Retrieve recommended size in pixels 'bufferSize' = [width, height] of the client
% renderbuffer for each eye for rendering to the HMD. Returns parameters
% previously computed by PsychOculusVR1('SetupRenderingParameters', hmd).
%
% Also returns 'imagingFlags', the required imaging mode flags for setup of
% the Screen imaging pipeline. Also returns the needed 'stereoMode' for the
% pipeline.
%
%
% needPanelFitter = PsychOculusVR1('GetPanelFitterParameters', hmd);
% - 'needPanelFitter' is 1 if a custom panel fitter tasks is needed, and 'bufferSize'
% from the PsychVRHMD('GetClientRenderingParameters', hmd); defines the size of the
% clientRect for the onscreen window. 'needPanelFitter' is 0 if no panel fitter is
% needed.
%
%
% [winRect, ovrfbOverrideRect, ovrSpecialFlags] = PsychOculusVR1('OpenWindowSetup', hmd, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags);
% - Compute special override parameters for given input/output arguments, as needed
% for a specific HMD. Take other preparatory steps as needed, immediately before the
% Screen('OpenWindow') command executes. This is called as part of PsychImaging('OpenWindow'),
% with the user provided hmd, screenid, winRect etc.
%
%
% isOutput = PsychOculusVR1('IsHMDOutput', hmd, scanout);
% - Returns 1 (true) if 'scanout' describes the video output to which the
% HMD 'hmd' is connected. 'scanout' is a struct returned by the Screen
% function Screen('ConfigureDisplay', 'Scanout', screenid, outputid);
% This allows probing video outputs to find the one which feeds the HMD.
% Deprecated: This function does nothing. It just exists for (backwards)
% compatibility with the old Oculus driver and PsychVRHMD.
%
%
% [headToEyeShiftv, headToEyeShiftMatrix] = PsychOculusVR1('GetEyeShiftVector', hmd, eye);
% - Retrieve 3D translation vector [tx, ty, tz] that defines the 3D position of the given
% eye 'eye' for the given HMD 'hmd', relative to the origin of the local head/HMD
% reference frame. This is needed to translate a global head pose into a eye
% pose, e.g., to translate the output of PsychOculusVR1('GetEyePose') into actual
% tracked/predicted eye locations for stereo rendering.
%
% In addition to the 'headToEyeShiftv' vector, a corresponding 4x4 translation
% matrix is also returned in 'headToEyeShiftMatrix' for convenience.
%
%

% Global GL handle for access to OpenGL constants needed in setup:
global GL;
global OVR;

persistent firsttime;
persistent oldShieldingLevel;
persistent hmd;

if nargin < 1 || isempty(cmd)
  help PsychOculusVR1;
  fprintf('\n\nAlso available are functions from PsychOculusVRCore1:\n');
  PsychOculusVRCore1;
  return;
end

% Fast-Path function 'EndFrameRender' - Queues new frames to Compositor:
if cmd == 0
  handle = varargin{1};

  t1 = GetSecs;

  % At this point, the FBO of the right eye texture in stereomode or mono texture
  % in mono mode is bound. Detach its color attachment texture:
  %glFramebufferTexture2D(GL.FRAMEBUFFER, GL.COLOR_ATTACHMENT0, GL.TEXTURE_2D, 0, 0);

  % Submit/Commit just unbound textures to texture swap-chains:
  PsychOculusVRCore1('EndFrameRender', hmd{handle}.handle, varargin{2});
  t2 = GetSecs;

  % Define parameters for the ongoing Psychtoolbox onscreen window flip operation:
  if hmd{handle}.mirrorTexture > 0
    % Debug output of compositor mirror texture into PTB onscreen window requested.
    % - Ask to skip flip's regular OpenGL swap completion timestamping, but instead
    %   to accept future injected timestamps from us.
    %
    % - Ask to disable vsync of the OpenGL bufferswap for display of the mirror texture
    %   in the onscreen window. We don't want to get swap-throttled to the refresh rate
    %   of the operator desktop GUI display.
    Screen('Hookfunction', hmd{handle}.win, 'SetOneshotFlipFlags', '', kPsychSkipWaitForFlipOnce + kPsychSkipVsyncForFlipOnce + kPsychSkipTimestampingForFlipOnce);
  else
    % No debug output from VR compositor wanted. Skip the OpenGL bufferswap for the
    % onscreen window completely, ergo also skip timestamping and allow timestamp
    % injection from us instead:
    Screen('Hookfunction', hmd{handle}.win, 'SetOneshotFlipFlags', '', kPsychSkipWaitForFlipOnce + kPsychSkipSwapForFlipOnce + kPsychSkipTimestampingForFlipOnce);
  end

  t3 = GetSecs;

  %fprintf('EndFrameRender: %f ms, SetNext/Flags %f ms ... \n', 1000 * (t2 - t1), 1000 * (t3 - t2));

  % If usercode requests a dontsync flag > 0 in Screen('Flip', window, when, dontclear, dontsync)
  % then switch to low precision/reliability timestamping, otherwise use timestampHighPrecision setting:
  if (varargin{3} > 0)
    hmd{handle}.doTimestamp = 0;
  else
    hmd{handle}.doTimestamp = hmd{handle}.timestampHighPrecision;
  end

  return;
end

% Fast-Path function 'PresentFrame' - Present frame to VR compositor,
% wait for present completion, render debug mirror texture if needed,
% inject present completion timestamps:
if cmd == 1
  handle = varargin{1};
  tWhen = varargin{2};

  t1 = GetSecs;

  [frameTiming, predictedOnset, refFrameIndex] = PsychOculusVRCore1('PresentFrame', hmd{handle}.handle, hmd{handle}.doTimestamp, tWhen);

  t2 = GetSecs;

  % Get fresh set of backing textures for next Screen() post-flip drawing/render
  % cycle from the OculusVR texture swap chains:
  texLeft = PsychOculusVRCore1('GetNextTextureHandle', hmd{handle}.handle, 0);
  if hmd{handle}.StereoMode > 0
    texRight = PsychOculusVRCore1('GetNextTextureHandle', hmd{handle}.handle, 1);
  else
    texRight = [];
  end

  % XXX NOTE: Without the following glFinish(), we get an invalid operation
  % OpenGL error during exec of 'PresentFrame' from within the Oculus VR
  % runtime, which however seems to have no negative consequence for
  % presentation?! This may be a bug in the OculusVR runtime at least as
  % of version 1.11:
  % glFinish;

  % Attach them as new backing textures, detach the previously bound ones, so they
  % are ready for submission to the VR compositor:
  [oldL, oldR] = Screen('Hookfunction', hmd{handle}.win, 'SetDisplayBufferTextures', '', texLeft, texRight);

  t3 = GetSecs;

  % Debug output from mirror texture requested?
  mirrorTex = hmd{handle}.mirrorTexture;
  if mirrorTex > 0
    % Yes. Render into onscreen windows viewport:
    width = hmd{handle}.mirrorWidth;
    height = hmd{handle}.mirrorHeight;
    glViewport(0, 0, width, height);
    glMatrixMode(GL.PROJECTION);
    glPushMatrix;
    glLoadIdentity;
    gluOrtho2D(0, width, 0, height);
    glBindTexture(GL.TEXTURE_2D, mirrorTex);
    glEnable(GL.TEXTURE_2D);
    glBegin(GL.QUADS);
    glColor4f(1,1,1,1);
    glTexCoord2f(0,1);
    glVertex2i(0,0);
    glTexCoord2f(1,1);
    glVertex2i(width,0);
    glTexCoord2f(1,0);
    glVertex2i(width,height);
    glTexCoord2f(0,0);
    glVertex2i(0,height);
    glEnd();
    glDisable(GL.TEXTURE_2D);
    glBindTexture(GL.TEXTURE_2D, 0);
    glPopMatrix;
    glMatrixMode(GL.MODELVIEW);
  end
  t4 = GetSecs;

  if hmd{handle}.doTimestamp
    % Assign return values for vblTime and stimulusOnsetTime for Screen('Flip'):
    if refFrameIndex > 0
      Screen('Hookfunction', hmd{handle}.win, 'SetOneshotFlipResults', '', frameTiming(refFrameIndex).VBlankTime, frameTiming(refFrameIndex).StimulusOnsetTime);
    else
      fprintf('Warning: No proper timestamp, faking it with predictedOnset.\n\n');
      Screen('Hookfunction', hmd{handle}.win, 'SetOneshotFlipResults', '', predictedOnset, predictedOnset);
    end
  else
    % Use made up values for timestamps of 'Flip':
    %Screen('Hookfunction', hmd{handle}.win, 'SetOneshotFlipResults', '', GetSecs, GetSecs);
    Screen('Hookfunction', hmd{handle}.win, 'SetOneshotFlipResults', '', predictedOnset, predictedOnset);
  end
  t5 = GetSecs;

  %fprintf('Present %f ms, Get/Set %f ms, Mirror %f ms, SetRes %f ms\n', 1000 * (t2 - t1), 1000 * (t3 - t2), 1000 * (t4 - t3), 1000 * (t5 - t4));
  %disp(frameTiming(1));
  % dT = 1e3 * (frameTiming(1).HMDTime - frameTiming(1).StimulusOnsetTime)
  frameTiming = [];
  return;
end



% Fast-Path function 'Cleanup' - Cleans up before onscreen window close/GL shutdown:
if cmd == 2
  handle = varargin{1};

  % Reattach old backing textures, so onscreen window can get properly destroyed:
  Screen('Hookfunction', hmd{handle}.win, 'SetDisplayBufferTextures', '', hmd{handle}.oldglLeftTex, hmd{handle}.oldglRightTex);
  hmd{handle}.oldglLeftTex = [];
  hmd{handle}.oldglRightTex = [];

  return;
end

if strcmpi(cmd, 'PrepareRender')
  % Get and validate handle - fast path open coded:
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR1:PrepareRender: Specified handle does not correspond to an open HMD!');
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
    % Default: Provide predicted value for the midpoint of the next video frame:
    targetTime = [];
  end

  % Mark start of a new frame render cycle for the runtime and get the data
  % predicted for next scanout time:
  [eyePose{1}, eyePose{2}] = PsychOculusVRCore1('StartRender', myhmd.handle);

  % Get predicted head pose, tracking state and hand poses (if supported) for targetTime:
  [state, touch] = PsychOculusVRCore1('GetTrackingState', myhmd.handle, targetTime);
  myhmd.state = state;

  % Always return basic tracking status:
  result.tracked = state.Status;
  result.SessionState = state.SessionState;

  if bitand(state.SessionState, 8)
    % DisplayLost condition! This is an unrecoverable error. Trigger
    % a forced session shutdown:
    error('OculusVR runtime reports loss of hardware (disconnected?) or serious malfunction. Forcing abort of this session.');
  end

  % Check if UX wants us to recenter the tracking origin:
  if bitand(state.SessionState, 32)
    % User wants us to recenter the tracking origin, so do it. Retry until it succeeds:
    while ~PsychOculusVRCore1('RecenterTrackingOrigin', myhmd.handle)
        WaitSecs('YieldSecs', 0.25);
    end
  end

  % As a bonus we return the raw eye pose vectors, given that we have them anyway:
  result.rawEyePose7{1} = eyePose{1};
  result.rawEyePose7{2} = eyePose{2};

  % Want matrices which take a usercode supplied global transformation into account?
  if bitand(reqmask, 1)
    % Yes: We need tracked + predicted head pose, so we can apply the user
    % transform, and then per-eye transforms:

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
    % Yes: We need tracked + predicted hand pose, so we can apply the user
    % transform, and then per-eye transforms:

    % Oculus 1.x SDK/runtime supports exactly 2 tracked touch controllers atm. to track users hands:
    for i=1:2
      result.handStatus(i) = touch(i).Status;

      % Bonus feature: HandPoses as 7 component translation + orientation quaternion vectors:
      result.handPose{i} = touch(i).HandPose;

      % Convert hand pose vector to 4x4 OpenGL right handed reference frame matrix:
      result.localHandPoseMatrix{i} = eyePoseToCameraMatrix(touch(i).HandPose);

      % Premultiply usercode provided global transformation matrix:
      result.globalHandPoseMatrix{i} = userTransformMatrix * result.localHandPoseMatrix{i};

      % Compute inverse matrix, maybe useable for collision testing / virtual grasping of virtual bjects:
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
    error('PsychOculusVR1:GetEyePose: Specified handle does not correspond to an open HMD!');
  end

  % Valid: Get view render pass for which to return information:
  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOculusVR1:GetEyePose: Required ''renderPass'' argument missing.');
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
  % NOTE: Currently not used, as Oculus SDK v1 does not support passing
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
  %    targetTime = [];
  %  end

  % Get eye pose for this renderPass, or more exactly the headPose from which this
  % renderPass eyePose will get computed:
  [result.eyePose, result.eyeIndex] = PsychOculusVRCore1('GetEyePose', myhmd.handle, renderPass);

  % Convert head pose vector to 4x4 OpenGL right handed reference frame matrix:
  result.localEyePoseMatrix = eyePoseToCameraMatrix(result.eyePose);

  % Premultiply usercode provided global transformation matrix:
  result.globalEyePoseMatrix = userTransformMatrix * result.localEyePoseMatrix;

  % Compute per-eye global pose matrix for this eyeIndex:
  result.cameraView = result.globalEyePoseMatrix;

  % Compute inverse matrix, useable as OpenGL GL_MODELVIEW matrix for rendering:
  result.modelView = inv(result.cameraView);

  varargout{1} = result;

  return;
end

if strcmpi(cmd, 'GetTrackersState')
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR1:GetTrackersState: Specified handle does not correspond to an open HMD!');
  end

  varargout{1} = PsychOculusVRCore1('GetTrackersState', myhmd.handle);

  return;
end

if strcmpi(cmd, 'GetInputState')
  % Get and validate handle - fast path open coded:
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR1:GetInputState: Specified handle does not correspond to an open HMD!');
  end

  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOculusVR1:GetInputState: Required ''controllerType'' argument missing.');
  end

  varargout{1} = PsychOculusVRCore1('GetInputState', myhmd.handle, double(varargin{2}));

  return;
end

if strcmpi(cmd, 'HapticPulse')
  % Get and validate handle - fast path open coded:
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR1:HapticPulse: Specified handle does not correspond to an open HMD!');
  end

  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOculusVR1:HapticPulse: Required ''controllerType'' argument missing.');
  end

  varargout{1} = PsychOculusVRCore1('HapticPulse', myhmd.handle, double(varargin{2}), varargin{3:end});

  return;
end

if strcmpi(cmd, 'TestVRBoundary')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('TestVRBoundary: Passed in handle does not refer to a valid and open HMD.');
  end

  [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = PsychOculusVRCore1('TestVRBoundary', myhmd.handle, varargin{2:end});
  return;
end

if strcmpi(cmd, 'TestVRBoundaryPoint')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('TestVRBoundaryPoint: Passed in handle does not refer to a valid and open HMD.');
  end

  [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = PsychOculusVRCore1('TestVRBoundaryPoint', myhmd.handle, varargin{2:end});
  return;
end

if strcmpi(cmd, 'VRAreaBoundary')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('VRAreaBoundary: Passed in handle does not refer to a valid and open HMD.');
  end

  [varargout{1}, varargout{2}, varargout{3}] = PsychOculusVRCore1('VRAreaBoundary', myhmd.handle, varargin{2:end});
  return;
end

if strcmpi(cmd, 'FloatsProperty')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('FloatsProperty: Passed in handle does not refer to a valid and open HMD.');
  end

  varargout{1} = PsychOculusVRCore1('FloatsProperty', myhmd.handle, varargin{2:end});
  return;
end

if strcmpi(cmd, 'StringProperty')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('StringProperty: Passed in handle does not refer to a valid and open HMD.');
  end

  varargout{1} = PsychOculusVRCore1('StringProperty', myhmd.handle, varargin{2:end});
  return;
end

if strcmpi(cmd, 'RecenterTrackingOrigin')
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR1:RecenterTrackingOrigin: Specified handle does not correspond to an open HMD!');
  end

  varargout{1} = PsychOculusVRCore1('RecenterTrackingOrigin', myhmd.handle);

  return;
end

if strcmpi(cmd, 'TrackingOriginType')
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR1:TrackingOriginType: Specified handle does not correspond to an open HMD!');
  end

  varargout{1} = PsychOculusVRCore1('TrackingOriginType', myhmd.handle, varargin{2:end});

  return;
end

if strcmpi(cmd, 'GetPerformanceStats')
  myhmd = varargin{1};
  if ~((length(hmd) >= myhmd.handle) && (myhmd.handle > 0) && hmd{myhmd.handle}.open)
    error('PsychOculusVR1:GetPerformanceStats: Specified handle does not correspond to an open HMD!');
  end

  [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = PsychOculusVRCore1('GetPerformanceStats', myhmd.handle);

  return;
end

if strcmpi(cmd, 'Supported')
  % Check if the Oculus VR runtime 1.x is supported and active on this
  % installation, so it can be used to open connections to real HMDs,
  % or at least to emulate a HMD for simple debugging purposes:
  try
    if IsWin(1) && exist('PsychOculusVRCore1', 'file') && PsychOculusVRCore1('GetCount') >= 0
      varargout{1} = 1;
    else
      varargout{1} = 0;
    end
  catch
    varargout{1} = 0;
  end
  return;
end

% Autodetect first connected HMD and open a connection to it. Open a
% emulated one, if none can be detected. Perform basic setup with
% default configuration, create a proper PsychImaging task.
if strcmpi(cmd, 'AutoSetupHMD')
  % Do we have basic runtime support?
  if ~PsychOculusVR1('Supported')
    % Nope: Game over.
    fprintf('PsychOculusVR1:AutoSetupHMD: Could not initialize Oculus V1 driver. Game over!\n');

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
    newhmd = PsychOculusVR1('Open', deviceIndex);
  else
    % Check if at least one Oculus HMD is connected and available:
    if PsychOculusVR1('GetCount') > 0
      % Yes. Open and initialize connection to first detected HMD:
      fprintf('PsychOculusVR1: Opening the first connected Oculus VR headset.\n');
      newhmd = PsychOculusVR1('Open', 0);
    else
      % TODO: HMD emulation under 1.x runtime possible?
      % No. Open an emulated/simulated HMD for basic testing and debugging:
      %fprintf('PsychOculusVR1: No Oculus HMD detected. Opening a simulated HMD.\n');
      %newhmd = PsychOculusVR1('Open', -1);
      fprintf('PsychOculusVR1: No Oculus HMD detected. Game over.\n');
      varargout{1} = [];
      return;
    end
  end

  % Trigger an automatic device close at onscreen window close for the HMD display window:
  PsychOculusVR1('SetAutoClose', newhmd, 1);

  % Setup default rendering parameters:
  PsychOculusVR1('SetupRenderingParameters', newhmd, basicTask, basicRequirements, basicQuality);

  % Add a PsychImaging task to use this HMD with the next opened onscreen window:
  PsychImaging('AddTask', 'General', 'UseVRHMD', newhmd);

  % Return the device handle:
  varargout{1} = newhmd;

  % Ready.
  return;
end

if strcmpi(cmd, 'SetAutoClose')
  myhmd = varargin{1};

  if ~PsychOculusVR1('IsOpen', myhmd)
    error('PsychOculusVR1:SetAutoClose: Specified handle does not correspond to an open HMD!');
  end

  % Assign autoclose flag:
  hmd{myhmd.handle}.autoclose = varargin{2};

  return;
end

if strcmpi(cmd, 'SetHSWDisplayDismiss')
  myhmd = varargin{1};

  if ~PsychOculusVR1('IsOpen', myhmd)
    error('PsychOculusVR1:SetHSWDisplay: Specified handle does not correspond to an open HMD!');
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

if strcmpi(cmd, 'SetHUDState')
  myhmd = varargin{1};

  if ~PsychOculusVR1('IsOpen', myhmd)
    error('PsychOculusVR1:SetHUDState: Specified handle does not correspond to an open HMD!');
  end

  % Method of dismissing HSW display:
  if length(varargin) < 2 || isempty(varargin{2})
    error('PsychOculusVR1:SetHUDState: Required mode argument missing!');
  end

  PsychOculusVRCore1('SetHUDState', myhmd.handle, varargin{2});

  return;
end

% Open a HMD:
if strcmpi(cmd, 'Open')
  if isempty(firsttime)
    firsttime = 1;
    fprintf('Copyright (c) 2018 Mario Kleiner. Licensed under the MIT license.\n');
    fprintf('The underlying PsychOculusVRCore1 mex driver uses the Oculus SDK, which is\n');
    fprintf('“Copyright © Facebook Technologies, LLC and its affiliates. All rights reserved.”\n');
    fprintf('A copy of the Oculus SDK license, its terms of use and thereby redistribution\n');
    fprintf('terms for the PsychOculusVRCore1 mex file can be found in the PsychtoolboxRoot()\n');
    fprintf('folder under:\n');
    fprintf('Psychtoolbox/PsychHardware/PsychVRToolbox/LICENSE_OculusSDK1.txt\n\n\n');
  end

  [handle, modelName, panelXRes, panelYRes, panelHz, controllerTypes] = PsychOculusVRCore1('Open', varargin{:});

  newhmd.handle = handle;
  newhmd.driver = @PsychOculusVR1;
  newhmd.type   = 'Oculus';
  newhmd.subtype = 'Oculus-1';
  newhmd.open = 1;
  newhmd.modelName = modelName;
  newhmd.separateEyePosesSupported = 0;
  newhmd.panelXRes = panelXRes;
  newhmd.panelYRes = panelYRes;
  newhmd.videoRefreshDuration = 1 / panelHz;
  newhmd.controllerTypes = controllerTypes;
  newhmd.VRControllersSupported = 1;
  newhmd.handTrackingSupported = 1;
  newhmd.hapticFeedbackSupported = 1;

  % Default autoclose flag to "no autoclose":
  newhmd.autoclose = 0;

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
    OVR.ControllerType_Object0 = hex2dec('0100');
    OVR.ControllerType_Object1 = hex2dec('0200');
    OVR.ControllerType_Object2 = hex2dec('0400');
    OVR.ControllerType_Object3 = hex2dec('0800');
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

  % High precision (and potentially high overhead) timestamping by default:
  % TODO FIXME Use low precision timestamping to avoid one cause of hangs atm.
  newhmd.timestampHighPrecision = 0;

  % Store in internal array:
  hmd{handle} = newhmd;

  % Return device struct:
  varargout{1} = newhmd;
  varargout{2} = modelName;

  return;
end

if strcmpi(cmd, 'Controllers')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('Controllers: Passed in handle does not refer to a valid and open HMD.');
  end

  varargout{1} = myhmd.controllerTypes;
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

if strcmpi(cmd, 'GetInfo')
  % Ok, cheap trick: We just return the passed in 'hmd' struct - the up to date
  % internal copy that is:
  if ~PsychOculusVR1('IsOpen', varargin{1})
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
      PsychOculusVRCore1('Close', myhmd.handle);
      hmd{myhmd.handle}.open = 0;
    end
  else
    % Shutdown whole driver:
    PsychOculusVRCore1('Close');
    hmd = [];
  end

  return;
end

if strcmpi(cmd, 'IsHMDOutput')
  myhmd = varargin{1}; %#ok<NASGU>
  scanout = varargin{2};

  % This does not make much sense on the OculusVR 1.0 version, as that runtime
  % only supports direct output mode, ie. an output completely separate from the
  % regular desktop and windowing system display space. But lets leave the legacy
  % code for the DK2 in place.

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

  % TODO FIXME: Any special handling for TimingSupport or Tracked3DVR ?
  if ~isempty(strfind(hmd{handle}.basicRequirements, 'TimingSupport')) || ...
     ~isempty(strfind(hmd{handle}.basicTask, 'Tracked3DVR'))
    %    PsychOculusVRCore1('SetDynamicPrediction', handle, 1);
  else
    %    PsychOculusVRCore1('SetDynamicPrediction', handle, 0);
  end

  return;
end

if strcmpi(cmd, 'SetFastResponse')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('SetFastResponse: Passed in handle does not refer to a valid and open HMD.');
  end
  handle = myhmd.handle;

  % FastResponse has no meaningful implementation on the OculusVR 1.0 runtime,
  % so just return a constant old value of 1 for "fast response always enabled":
  varargout{1} = 1;

  return;
end

if strcmpi(cmd, 'SetTimeWarp')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('SetTimeWarp: Passed in handle does not refer to a valid and open HMD.');
  end

  % SetTimeWarp determined use of GPU accelerated 2D texture sampling
  % warp on the Rift DK1/DK2 with old v0.5 SDK. On the 1.0 SDK we no
  % longer have any programmatic control over timewarping,so leave this
  % in place as dummy.

  % Return constant old setting of "TimeWarp always on":
  varargout{1} = 1;

  return;
end

if strcmpi(cmd, 'SetLowPersistence')
  myhmd = varargin{1};
  if ~PsychOculusVR1('IsOpen', myhmd)
    error('SetLowPersistence: Passed in handle does not refer to a valid and open HMD.');
  end

  % SetLowPersistence determined use of low persistence mode on the Rift DK2 with
  % the v0.5 SDK. We don't have control over this on the v1.0 runtime anymore.
  % Return constant old setting "Always low persistence":
  varargout{1} = 1;

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

  PsychOculusVR1('SetBasicQuality', myhmd, basicQuality);

  % Get optimal client renderbuffer size - the size of our virtual framebuffer for left eye:
  [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovL] = PsychOculusVRCore1('GetFovTextureSize', myhmd.handle, 0, varargin{5:end});

  % Get optimal client renderbuffer size - the size of our virtual framebuffer for right eye:
  [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovR] = PsychOculusVRCore1('GetFovTextureSize', myhmd.handle, 1, varargin{5:end});

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
    [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovL] = PsychOculusVRCore1('GetFovTextureSize', myhmd.handle, 0, fov, varargin{6:end});
    [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight, hmd{myhmd.handle}.fovR] = PsychOculusVRCore1('GetFovTextureSize', myhmd.handle, 1, fov, varargin{6:end});
  end

  % Debug display of HMD output into onscreen window requested?
  if isempty(strfind(basicRequirements, 'DebugDisplay')) && isempty(oldShieldingLevel)
    % No. Set to be created onscreen window to be invisible:
    oldShieldingLevel = Screen('Preference', 'WindowShieldingLevel', -1);
  end

  % HUD display requested?
  if ~isempty(strfind(basicRequirements, 'HUD='))
    % Yes. Set it:
    hudmodestring = basicRequirements(strfind(basicRequirements, 'HUD='):end);
    mode = sscanf(hudmodestring, 'HUD=%i');
    PsychOculusVR1('SetHUDState', myhmd, mode);
  end

  % Use low precision timestamps if usercode requests them, because they can be
  % computed with much lower overhead:
  if ~isempty(strfind(basicRequirements, 'FastLowPrecisionTimestamps'))
    hmd{myhmd.handle}.timestampHighPrecision = 0;
  end

  return;
end

if strcmpi(cmd, 'GetClientRenderingParameters')
  myhmd = varargin{1};
  varargout{1} = [hmd{myhmd.handle}.rbwidth, hmd{myhmd.handle}.rbheight];

  % We need fast backing store support for imaging pipeline with virtual framebuffers,
  % also output redirection for providing final output to us, instead of displaying
  % into the standard onscreen window. The OculusVR 1.x runtime generates its own
  % textures to be used as externally injected color buffer backing textures:
  imagingMode = mor(kPsychNeedFastBackingStore, kPsychNeedFinalizedFBOSinks, kPsychUseExternalSinkTextures);

  % Usercode wants a 16 bpc half-float rendering pipeline?
  if ~isempty(strfind(hmd{myhmd.handle}.basicRequirements, 'Float16Display'))
    % Request a 16 bpc float framebuffer from Psychtoolbox:
    imagingMode = mor(imagingMode, kPsychNeed16BPCFloat);
  else
    % Standard RGBA8 images: Use sRGB format for rendering/blending/compositing/display:
    imagingMode = mor(imagingMode, kPsychEnableSRGBRendering);
  end

  if ~strcmpi(hmd{myhmd.handle}.basicTask, 'Monoscopic')
    % We must use stereomode 12, so we get separate draw buffers for left and
    % right eye, and separate stream processing into our VR runtime provided
    % separate backing textures / texture swap chains per eye, with all internal
    % buffers at at least full VR compositor input resolution.
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
  % We don't need custom PanelFitter task for the v1.0+ Oculus SDK/Runtime:
  varargout{1} = 0;
  return;
end

% [winRect, ovrfbOverrideRect, ovrSpecialFlags] = PsychOculusVR1('OpenWindowSetup', hmd, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags);
if strcmpi(cmd, 'OpenWindowSetup')
  myhmd = varargin{1};
  screenid = varargin{2};
  winRect = varargin{3};
  ovrfbOverrideRect = varargin{4};
  ovrSpecialFlags = varargin{5};

  % The current design iteration requires the PTB parent onscreen windows
  % effective backbuffer (from the pov of the imaging pipeline) to have the
  % same size (width x height) as the renderbuffer for one eye, so enforce
  % that constraint by setting ovrfbOverrideRect accordingly.

  % Get required output buffer size:
  clientRes = myhmd.driver('GetClientRenderingParameters', myhmd);

  % Set as fbOverrideRect for window:
  ovrfbOverrideRect = [0, 0, clientRes(1), clientRes(2)];

  fprintf('PsychOculusVR1-Info: Overriding onscreen window framebuffer size to %i x %i pixels for use with VR-HMD direct output mode.\n', ...
          clientRes(1), clientRes(2));

  % As the onscreen window is not used for displaying on the HMD, but
  % either not at all, or just for debug output, make it a regular GUI
  % window, managed by the window manager, so user can easily get it out
  % of the way:
  ovrSpecialFlags = mor(ovrSpecialFlags, kPsychGUIWindow + kPsychGUIWindowWMPositioned);

  % Skip all visual timing sync tests and calibrations, as display timing
  % of the onscreen window doesn't matter, only the timing on the HMD direct
  % output matters - and that can't be measured by our standard procedures:
  Screen('Preference', 'SkipSyncTests', 2);

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

if strcmpi(cmd, 'Start') || strcmpi(cmd, 'Stop')
    % Make sure the 'Start' and 'Stop' commands, which are mandatory to support
    % for backwards compatibility (see PsychVRHMD), do absolutely nothing in this
    % driver.
    return;
end

if strcmpi(cmd, 'PerformPostWindowOpenSetup')
  % Must have global GL constants:
  if isempty(GL)
    varargout{1} = 0;
    warning('PTB internal error in PsychOculusVR1: GL struct not initialized?!?');
    return;
  end

  % Oculus device handle:
  myhmd = varargin{1};
  handle = myhmd.handle;

  % Onscreen window handle:
  win = varargin{2};
  winfo = Screen('GetWindowInfo', win);
  hmd{handle}.StereoMode = winfo.StereoMode;

  % Keep track of window handle of associated onscreen window:
  hmd{handle}.win = win;

  % Restore shielding level for new windows after "our" onscreen window is now open:
  if ~isempty(oldShieldingLevel)
    Screen('Preference', 'WindowShieldingLevel', oldShieldingLevel);
    oldShieldingLevel = [];
  end

  % Set override window parameters with pixel size (color depth) and refresh
  % interval in seconds as provided by the VR runtime:
  Screen('HookFunction', win, 'SetWindowBackendOverrides', [], 24, hmd{handle}.videoRefreshDuration);

  % Need to know user selected clearcolor:
  clearcolor = varargin{3};

  % Compute effective size of per-eye input buffer for undistortion render.
  % The input buffers for undistortion are the finalizedDrawbufferFBO's, ie.
  % the final output buffers of Psychtoolbox imaging pipeline.
  % This means [inputWidth, inputHeight] == [rbwidth, rbheight].
  hmd{handle}.inputWidth = hmd{handle}.rbwidth;
  hmd{handle}.inputHeight = hmd{handle}.rbheight;

  % Query undistortion parameters for left eye view:
  [hmd{handle}.rbwidth, hmd{handle}.rbheight, vx, vy, vw, vh, ptx, pty, hsx, hsy, hsz] = PsychOculusVRCore1('GetUndistortionParameters', handle, 0, hmd{handle}.inputWidth, hmd{handle}.inputHeight, hmd{handle}.fovL);
  hmd{handle}.viewportLeft = [vx, vy, vw, vh];
  hmd{handle}.PixelsPerTanAngleAtCenterLeft = [ptx, pty];
  hmd{handle}.HmdToEyeViewOffsetLeft = [hsx, hsy, hsz];

  % Query parameters for right eye view:
  [hmd{handle}.rbwidth, hmd{handle}.rbheight, vx, vy, vw, vh, ptx, pty, hsx, hsy, hsz] = PsychOculusVRCore1('GetUndistortionParameters', handle, 1, hmd{handle}.inputWidth, hmd{handle}.inputHeight, hmd{handle}.fovR);
  hmd{handle}.viewportRight = [vx, vy, vw, vh];
  hmd{handle}.PixelsPerTanAngleAtCenterRight = [ptx, pty];
  hmd{handle}.HmdToEyeViewOffsetRight = [hsx, hsy, hsz];

  % Convert head to eye shift vectors into 4x4 matrices, as we'll need
  % them frequently:
  EyeT = diag([1 1 1 1]);
  EyeT(1:3, 4) = hmd{handle}.HmdToEyeViewOffsetLeft';
  hmd{handle}.eyeShiftMatrix{1} = EyeT;

  EyeT = diag([1 1 1 1]);
  EyeT(1:3, 4) = hmd{handle}.HmdToEyeViewOffsetRight';
  hmd{handle}.eyeShiftMatrix{2} = EyeT;

  % Create texture swap chains to provide textures to be used for
  % frame submission to the VR compositor:
  if ~isempty(strfind(hmd{myhmd.handle}.basicRequirements, 'Float16Display'))
    % Linear RGBA16F half-float textures as target framebuffers:
    floatFlag = 1;
  else
    % sRGB RGBA8 textures as target framebuffers:
    floatFlag = 0;
  end

  % Left eye / mono chain:
  [width, height, numTextures] = PsychOculusVRCore1('CreateRenderTextureChain', hmd{handle}.handle, 0, hmd{handle}.inputWidth, hmd{handle}.inputHeight, floatFlag);

  % Create 2nd chain for right eye in stereo mode:
  if winfo.StereoMode > 0
    if winfo.StereoMode ~=12
      sca;
      error('Invalid Screen() StereoMode in use for OculusVR HMD! Must be 12 or 0.');
    end
    [width, height, numTextures] = PsychOculusVRCore1('CreateRenderTextureChain', hmd{handle}.handle, 1, hmd{handle}.inputWidth, hmd{handle}.inputHeight, floatFlag);
  end

  % Query currently bound finalizedFBO backing textures, to keep them around as
  % backups for restoration when closing down the session:
  [hmd{handle}.oldglLeftTex, hmd{handle}.oldglRightTex, textarget, texformat, texmultisample, texwidth, texheight] = Screen('Hookfunction', win, 'GetDisplayBufferTextures');
  if (textarget ~= GL.TEXTURE_2D) || (texformat ~= GL.RGBA8 && texformat ~= GL.RGBA16F) || (texmultisample ~= 0)
    sca;
    error('Invalid Screen() backing textures required. Non-matching texture target, format or multisample setting.');
  end

  if (texwidth ~= width) || (texheight ~= height)
    sca;
    error('Invalid Screen() backing textures required. Non-matching width x height.');
  end

  if hmd{handle}.oldglRightTex == 0
    hmd{handle}.oldglRightTex = [];
  end

  % Get and clear all backing textures from the VR compositor:
  clearvalues = ones(4, texwidth, texheight, 'single');
  for i=1:numTextures
    % Get backing textures provided by the textures swap chains of the OculusVR compositor
    % and clear them to black color by zero-filling:
    texLeft = PsychOculusVRCore1('GetNextTextureHandle', hmd{handle}.handle, 0);
    glBindTexture(textarget, texLeft);
    glTexSubImage2D(textarget, 0, 0, 0, texwidth, texheight, GL.RGBA, GL.FLOAT, clearvalues);

    if winfo.StereoMode > 0
      texRight = PsychOculusVRCore1('GetNextTextureHandle', hmd{handle}.handle, 1);
      glBindTexture(textarget, texRight);
      glTexSubImage2D(textarget, 0, 0, 0, texwidth, texheight, GL.RGBA, GL.FLOAT, clearvalues);
    else
      texRight = [];
    end
    glBindTexture(textarget, 0);

    % Execute VR Present: Submit just unbound textures, start render/warp/presentation
    % on HMD at next possible point in time:
    PsychOculusVRCore1('EndFrameRender', hmd{handle}.handle, 0);
    PsychOculusVRCore1('PresentFrame', hmd{handle}.handle, 0, -1);
  end
  clear clearvalues;

  if ~isempty(strfind(hmd{myhmd.handle}.basicTask, 'Tracked3DVR'))
    % 3D head tracked VR rendering task: Start tracking as a convenience:
    PsychOculusVRCore1('Start', handle);
  else
    % No 3D head tracked closed-loop mode: Kick off the presenter thread via
    % this 'Start' -> 'Stop' sequence:
    PsychOculusVRCore1('Start', handle);
    PsychOculusVRCore1('Stop', handle);
  end
  Screen('GetWindowInfo', win);
  % Get first textures for actual use in PTB's imaging pipeline:
  texLeft = PsychOculusVRCore1('GetNextTextureHandle', hmd{handle}.handle, 0);
  if hmd{handle}.StereoMode > 0
      texRight = PsychOculusVRCore1('GetNextTextureHandle', hmd{handle}.handle, 1);
  else
      texRight = [];
  end

  % Assign them to Screen():
  Screen('Hookfunction', win, 'SetDisplayBufferTextures', '', texLeft, texRight);

  % Go back to user requested clear color, now that all our buffers
  % are cleared to black:
  Screen('FillRect', win, clearcolor);

  % Store true size of onscreen window as mirrorWidth x mirrorHeight:
  winGRect = Screen('GlobalRect', win);
  hmd{handle}.mirrorWidth = RectWidth(winGRect);
  hmd{handle}.mirrorHeight = RectHeight(winGRect);

  % Debug display of HMD output into onscreen window requested?
  if ~isempty(strfind(hmd{handle}.basicRequirements, 'DebugDisplay'))
    % Yes. Create OpenGL mirror texture which receives VR compositor images:
    hmd{handle}.mirrorTexture = PsychOculusVRCore1('CreateMirrorTexture', hmd{handle}.handle, hmd{handle}.panelXRes, hmd{handle}.panelYRes);
  else
    hmd{handle}.mirrorTexture = 0;
  end

  % Need to call the PsychOculusVR1(0) callback at each Screen('Flip') to get the
  % imaging pipelines post-processed final output frames for left/right eye and commit
  % them to the VR-Compositors texture swap-chain(s), then setting up new target textures
  % as buffers for next cycle. This happens at the end of preflip operations, as part of
  % the implicit 'DrawingFinished':
  cmdString = sprintf('PsychOculusVR1(0, %i, IMAGINGPIPE_FLIPTWHEN, IMAGINGPIPE_FLIPVBLSYNCLEVEL);', handle);
  if winfo.StereoMode > 0
    % In stereo mode we also need to detach the left eye texture as commit prep:
    %detachCmd = 'glFramebufferTexture2D(36160, 36064, 3553, 0, 0);';
    %Screen('Hookfunction', win, 'AppendMFunction', 'LeftFinalizerBlitChain', 'OculusVR left texture detach Operation', detachCmd);
    %Screen('Hookfunction', win, 'Enable', 'LeftFinalizerBlitChain');

    % In stereo mode, use right finalizer chain for right texture detach and both textures commit, as it executes last:
    Screen('Hookfunction', win, 'AppendMFunction', 'RightFinalizerBlitChain', 'OculusVR Stereo commit Operation', cmdString);
    Screen('Hookfunction', win, 'Enable', 'RightFinalizerBlitChain');
  else
    % In mono mode, use left finalizer chain for texture detach and commit, as it executes only - and therefore last:
    Screen('Hookfunction', win, 'AppendMFunction', 'LeftFinalizerBlitChain', 'OculusVR Mono commit Operation', cmdString);
    Screen('Hookfunction', win, 'Enable', 'LeftFinalizerBlitChain');
  end

  % Need to call the PsychOculusVR1(1) callback at each Screen('Flip') to submit the
  % output frames to the VR-Compositor for presentation on the HMD asap. This gets
  % called immediately before an OpenGL bufferswap (if any) + timestamping + validation
  % will happen, iow. at t >= tWhen for Screen('Flip', win, tWhen);
  % It is supposed to block until image presentation on the HMD has happened, and to
  % inject proper Present timestamps for 'Flip':
  % TODO FIXME: Nope, change of strategy: We set the kPsychSkipWaitForFlipOnce flag in
  % 'EndFrameRender', so PreSwapbuffersOperations does not wait until tWhen, but executes
  % immediately. We pass in the tWhen timestamp to this fast-path callback, which will pass
  % it on to 'PresentFrame', and then it is the job of the PsychOculusVRCore1 drivers thread
  % to wait until the right moment to submit the new frame to the VR compositor. This allows
  % to use a single-threaded presentation model, where the single thread busy-waits until the
  % right point in time. We retain the ability to 'Flip' as fast as possible for closed-loop
  % tracked VR presentation, and the ability to do timed open-loop presentation with 'tWhen'
  % target times, and avoid all the multi-threading misery the Oculus runtime brings to us.
  % Downside is the inability to avoid the sandclock warning if usercode truly does not execute
  % any Flip's for more than a second, but this has to do as a stepping stone for now:
  cmdString = sprintf('PsychOculusVR1(1, %i, IMAGINGPIPE_FLIPTWHEN);', handle);
  Screen('Hookfunction', win, 'AppendMFunction', 'PreSwapbuffersOperations', 'OculusVR Present Operation', cmdString);
  Screen('Hookfunction', win, 'Enable', 'PreSwapbuffersOperations');

  % Attach shutdown procedure on onscreen window close:
  cmdString = sprintf('PsychOculusVR1(2, %i);', handle);
  Screen('Hookfunction', win, 'PrependMFunction', 'CloseOnscreenWindowPreGLShutdown', 'OculusVR cleanup', cmdString);
  Screen('Hookfunction', win, 'Enable', 'CloseOnscreenWindowPreGLShutdown');

  % Does usercode request auto-closing the HMD or driver when the onscreen window is closed?
  if hmd{handle}.autoclose > 0
    % Attach a window close callback for Device teardown at window close time:
    if hmd{handle}.autoclose == 2
      % Shutdown driver completely:
      Screen('Hookfunction', win, 'AppendMFunction', 'CloseOnscreenWindowPreGLShutdown', 'Shutdown window callback into PsychOculusVR1 driver.', 'PsychOculusVR1(''Close'');');
    else
      % Only close this HMD:
      Screen('Hookfunction', win, 'PrependMFunction', 'CloseOnscreenWindowPreGLShutdown', 'Shutdown window callback into PsychOculusVR1 driver.', sprintf('PsychOculusVR1(''Close'', %i);', handle));
    end
    Screen('HookFunction', win, 'Enable', 'CloseOnscreenWindowPreGLShutdown');
  end

  % Return success result code 1:
  varargout{1} = 1;
  return;
end

% 'cmd' so far not dispatched? Let's assume it is a command
% meant for PsychOculusVRCore1:
if (length(varargin) >= 1) && isstruct(varargin{1})
  myhmd = varargin{1};
  handle = myhmd.handle;
  [ varargout{1:nargout} ] = PsychOculusVRCore1(cmd, handle, varargin{2:end});
else
  [ varargout{1:nargout} ] = PsychOculusVRCore1(cmd, varargin{:});
end

return;

end
