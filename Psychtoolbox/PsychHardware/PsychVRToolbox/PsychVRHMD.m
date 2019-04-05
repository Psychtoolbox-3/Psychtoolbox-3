function varargout = PsychVRHMD(cmd, varargin)
% PsychVRHMD() - High level driver for VR HMD devices.
%
% This driver bundles all the common high level functionality
% of different Virtual Reality Head mounted displays into one
% function.
%
% It dispatches generic calls into appropriate device specific
% drivers as needed.
%
% Usage:
%
% hmd = PsychVRHMD('AutoSetupHMD' [, basicTask][, basicRequirements][, basicQuality][, vendor][, deviceIndex]);
% - Automatically detect the first connected HMD, set it up with reasonable
% default parameters, and return a device handle 'hmd' to it. If the system
% does not support any HMDs, not even emulated ones, just does nothing and
% returns an empty handle, ie., hmd = [], so caller can cope with that.
%
% Optional parameters: 'basicTask' what kind of task should be implemented.
% The default is 'Tracked3DVR', which means to setup for stereoscopic 3D
% rendering, driven by head motion tracking, for a fully immersive experience
% in some kind of 3D virtual world. This is the default if omitted. The task
% 'Stereoscopic' sets up for display of stereoscopic stimuli, but without
% head tracking. 'Monoscopic' sets up for display of monocular stimuli, ie.
% the HMD is just used as a special kind of standard display monitor.
%
% 'basicRequirements' defines basic requirements for the task. Currently
% defined are the following strings which can be combined into a single
% 'basicRequirements' string:
%
% 'LowPersistence' = Try to keep exposure time of visual images on the retina
% low if possible, ie., try to approximate a pulse-type display instead of a
% hold-type display if possible. This has no effect on the Oculus Rift DK1.
% On the Rift DK2 it will enable low persistence scanning of the OLED display
% panel, to light up each pixel only a fraction of a video refresh cycle duration.
% On the Rift CV1, low persistence is always active, so this setting is redundant.
%
% 'DebugDisplay' = Show the output which is displayed on the HMD inside the
% Psychtoolbox onscreen window as well. This will have a negative impact on
% performance, latency and timing of the HMD visual presentation, so should only
% be used for debugging, as it may cause a seriously degraded VR experience.
% By default, no such debug output is produced and the Psychtoolbox onscreen
% window is not actually displayed on the desktop. This option is silently ignored
% on the old classic Oculus driver at the moment.
%
% 'Float16Display' = Request rendering, compositing and display in 16 bpc float
% format on some HMD's and drivers. This will ask Psychtoolbox to render and
% post-process stimuli in 16 bpc linear floating point format, and allocate 16 bpc
% half-float textures as final renderbuffers to be sent to the VR compositor.
% If the VR compositor takes advantage of the high source image precision is at
% the discretion of the compositor and HMD. By default, if this request is omitted,
% processing and display in sRGB format is requested from Psychtoolbox and the compositor
% on some drivers, e.g., for the Oculus 1.11+ runtime and Rift CV1, ie., a roughly
% gamma 2.2 8 bpc format is used, which is optimized for the gamma response curve of
% at least the Oculus Rift CV1 display. On other HMDs or drivers, this option may do
% nothing and get silently ignored.
%
% 'PerEyeFOV' = Request use of per eye individual and asymmetric fields of view, even
% when the 'basicTask' was selected to be 'Monoscopic' or 'Stereoscopic'. This allows
% for wider field of view in these tasks, but requires the usercode to adapt to these
% different and asymmetric fields of view for each eye, e.g., by selecting proper 3D
% projection matrices for each eye. If a 'basicTask' of '3DVR' for non-tracked 3D, or
% (the default) 'Tracked3DVR' for head tracking driven 3D is selected, then that implies
% per-eye individual and asymmetric fields of view, iow. 'PerEyeFOV' is implied. For pure
% 'basicTask' of 'Monoscopic' or 'Stereoscopic' for Screen() 2D drawing, the system uses
% identical and symmetric fields of view for both eyes by default, so 'PerEyeFOV' would
% be needed to override this choice. COMPATIBILITy NOTE: Psychtoolbox-3 releases before
% June 2017 always used identical and symmetric fields of view for both eyes, which was
% a bug. However the error made was very small, due to the imaging properties of the
% Oculus Rift DK2, essentially imperceptible to the unknowing observer with the naked
% eye. Releases starting June 2017 now use separate fields of view in 3D rendering
% modes, and optionally for 2D mono/stereo modes with this 'PerEyeFOV' opt-in parameter,
% so stimulus display may change slightly for the same HMD hardware and user-code,
% compared to older Psychtoolbox-3 releases. This change was crucial to accomodate the
% rather different imaging properties of the Oculus Rift CV1 and possible other future
% HMD's.
%
% 'FastResponse' = Try to switch images with minimal delay and fast
% pixel switching time. This will enable OLED panel overdrive processing
% on the Oculus Rift DK1 and DK2. OLED panel overdrive processing is a
% relatively expensive post processing step. On the Rift CV1, and generally
% with the new Oculus v1.11+ runtime, this is always active, so 'FastResponse'
% is redundant on such panels and drivers.
%
% 'TimingSupport' = Support some hardware specific means of timestamping
% or latency measurements. On the Rift DK1 this does nothing. On the DK2
% it enables dynamic prediction and timing measurements with the Rifts internal
% latency tester. This does nothing anymore on Rift CV1.
%
% 'TimeWarp' = Enable per eye image 2D timewarping via prediction of eye
% poses at scanout time. This mostly only makes sense for head-tracked 3D
% rendering. Depending on 'basicQuality' a more cheap or more expensive
% procedure is used. On the v1.11 Oculus runtime and Rift CV1, 'TimeWarp'
% is always active, so this option is redundant.
%
% These basic requirements get translated into a device specific set of
% settings. The settings can also be specific to the selected 'basicTask',
% and if a quality vs. performance / system load tradeoff is unavoidable
% then the 'basicQuality' parameter may modulate the strategy.
%
% 'basicQuality' defines the basic tradeoff between quality and required
% computational power. A setting of 0 gives lowest quality, but with the
% lowest performance requirements. A setting of 1 gives maximum quality at
% maximum computational load. Values between 0 and 1 change the quality to
% performance tradeoff.
%
% By default all currently supported types of HMDs from different
% vendors are probed and the first one found is used. If the optional
% parameter 'vendor' is provided, only devices from that vendor are
% detected and the first detected device is chosen.
%
% If additionally the optional 'deviceIndex' parameter is provided then
% that specific device 'deviceIndex' from that 'vendor' is opened and set up.
%
%
% PsychVRHMD('SetAutoClose', hmd, mode);
% - Set autoclose mode for HMD with handle 'hmd'. 'mode' can be
% 0 (this is the default) to not do anything special. 1 will close
% the HMD 'hmd' when the onscreen window is closed which displays
% on the HMD. 2 will do the same as 1, but close all open HMDs and
% shutdown the complete driver and runtime - a full cleanup.
%
%
% isOpen = PsychVRHMD('IsOpen', hmd);
% - Returns 1 if 'hmd' corresponds to an open HMD, 0 otherwise.
%
%
% PsychVRHMD('Close' [, hmd])
% - Close provided HMD device 'hmd'. If no 'hmd' handle is provided,
% all HMDs will be closed and the driver will be shutdown.
%
%
% PsychVRHMD('Controllers', hmd);
% - Return a bitmask of all connected controllers: Can be the bitand
% of the OVR.ControllerType_XXX flags described in 'GetInputState'.
% This does not detect if controllers are hot-plugged or unplugged after
% the HMD was opened. Iow. only probed at 'Open'.
%
%
% info = PsychVRHMD('GetInfo', hmd);
% - Retrieve a struct 'info' with information about the HMD 'hmd'.
% The returned info struct contains at least the following standardized
% fields with information:
% handle = Driver internal handle for the specific HMD.
% driver = Function handle to the actual driver for the HMD, e.g., @PsychOculusVR.
% type   = Defines the type/vendor of the device, e.g., 'Oculus'.
% subtype = Defines the type of driver more specific, e.g., 'Oculus-classic' or 'Oculus-1'.
% modelName = Name string with the name of the model of the device, e.g., 'Rift DK2'.
%
% separateEyePosesSupported = 1 if use of PsychVRHMD('GetEyePose') will improve
%                             the quality of the VR experience, 0 if no improvement
%                             is to be expected, so 'GetEyePose' can be avoided
%                             to save processing time without a loss of quality.
%
% VRControllersSupported = 1 if use of PsychVRHMD('GetInputState') will provide input
%                            from actual dedicated VR controllers. Value is 0 if
%                            controllers are only emulated to some limited degree,
%                            e.g., by abusing a regular keyboard as a button controller,
%                            ie. mapping keyboard keys to OVR.Button_XXX buttons.
%
% handTrackingSupported = 1 if PsychVRHMD('PrepareRender') with reqmask +2 will provide
%                           valid hand tracking info, 0 if this is not supported and will
%                           just report fake values. A driver may report 1 here but still
%                           don't provide meaningful info at runtime, e.g., if required
%                           tracking hardware is missing or gets disconnected. The flag
%                           just aids extra performance optimizations in your code.
%
% hapticFeedbackSupported = 1 if basic haptic feedback is supported in principle on some controllers.
%                             0 otherwise. A flag of zero means no haptic feedback support, but
%                             a flag of 1 may still mean no actual feedback, e.g., if suitable
%                             hardware is not configured and present. Flags higher than 1 can
%                             signal presence of more advanced haptic feedback, so you should
%                             test for a setting == 1 to know if PsychVRHMD('HapticPulse') works
%                             in principle, which is considered basic feedback ability.
%
% The info struct may contain much more vendor specific information, but the above
% set is supported across all devices.
%
%
% [isVisible, playAreaBounds, OuterAreaBounds] = PsychVRHMD('VRAreaBoundary', hmd [, requestVisible]);
% - Request visualization of the VR play area boundary for 'hmd' and returns its
% current extents.
%
% 'requestVisible' 1 = Request showing the boundary area markers, 0 = Don't
% request showing the markers.
% If the driver can control area boundary visibility is highly dependent on the VR
% driver in use. This flag may get ignored. See driver specific help, e.g.,
% "help PsychOculusVR1", for behaviour of a specific driver.
%
% Some drivers or hardware setups may not supports VR area boundaries at all, in
% which case the function will return empty boundaries.
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
% input = PsychVRHMD('GetInputState', hmd, controllerType);
% - Get input state of controller 'controllerType' associated with HMD 'hmd'.
%
% Note that if the underlying driver does not support special VR controllers, ie.,
% hmdinfo = PsychVRHMD('GetInfo') returns hmdinfo.VRControllersSupported == 0, then
% only a minimally useful 'input' state is returned, which is based on emulating or
% faking input from real controllers, so this function will be of limited use. Specifically,
% on emulated controllers, only the input.Valid, input.Time and input.Buttons
% fields are returned, all other fields are missing.
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
% pulseEndTime = PsychVRHMD('HapticPulse', hmd, controllerType [, duration=XX][, freq=1.0][, amplitude=1.0]);
% - Trigger a haptic feedback pulse, some controller vibration, on the specified 'controllerType'
% associated with the specified 'hmd'. 'duration' is desired pulse duration in seconds. On Oculus
% devices, by default a maximum of 2.5 seconds pulse is executed, but other vendors devices may have
% a different maximum. 'freq' is normalized frequency in range 0.0 - 1.0. A value of 0 will try to
% disable an ongoing pulse. How this normalized 'freq' maps to a specific haptic device is highly
% device and runtime dependent.
%
% 'amplitude' is the amplitude of the vibration in normalized 0.0 - 1.0 range.
%
% 'pulseEndTime' returns the expected stop time of vibration in seconds, given the parameters.
% Currently the function will return immediately for a (default) 'duration', and the pulse
% will end after the maximum duration supported by the given device. Smaller 'duration' values than
% the maximum duration will block the execution of the function until the 'duration' has passed on
% some types of controllers.
%
% Please note that behaviour of this function is highly dependent on the type of VR driver and
% devices used. You should consult driver specific documentation for details, e.g., the help of
% 'PsychOculusVR' or 'PsychOculusVR1' for Oculus systems. On some drivers the function may do
% nothing at all, e.g., if the 'GetInfo' function returns info.hapticFeedbackSupported == 0.
%
%
% state = PsychVRHMD('PrepareRender', hmd [, userTransformMatrix][, reqmask=1][, targetTime]);
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
%      and orientations of left and right hand tracking controllers, if any. See also
%      section about 'GetInfo' for some performance comments.
%
%      state.handStatus(1) = Tracking status of left hand: 0 = Untracked, 1 = Orientation
%                            tracked, 2 = Position tracked, 3 = Orientation and position
%                            tracked. If handStatus is == 0 then all the following information
%                            is invalid and can not be used in any meaningful way.
%
%      state.handStatus(2) = Tracking status of right hand.
%
%      state.localHandPoseMatrix{1} = 4x4 OpenGL right handed reference frame matrix with
%                                     hand position and orientation encoded to define a
%                                     proper GL_MODELVIEW transform for rendering stuff
%                                     "into"/"relative to" the oriented left hand.
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
% eyePose = PsychVRHMD('GetEyePose', hmd, renderPass [, userTransformMatrix][, targetTime]);
% - Return a struct 'eyePose' which contains various useful bits of information
% for 3D stereoscopic rendering of the stereo view of one eye, based on head
% tracking data. This function provides essentially the same information as
% the 'PrepareRender' function, but only for one eye. Therefore you will need
% to call this function twice, once for each of the two renderpasses, at the
% beginning of each renderpass.
%
% Note: The 'GetEyePose' function may not be implemented in a meaningful/beneficial
% way for all supported types of HMD. This means that while the function will work
% on all supported HMDs, there may not be any benefit of using it in terms of
% performance or quality of the VR experience, because the underlying driver may
% only emulate / fake the results for compatibility. Currently only the driver for
% the Oculus VR Rift DK1 and DK2 supports this function in a way that will improve
% the VR experience, the status for future Oculus HMDs, or HMDs from other vendors
% is currently unknown. The info struct returned by PsychVRHMD('GetInfo') will return
% info.separateEyePosesSupported == 1 if there is a benefit to be expected from use
% of this function, or info.separateEyePosesSupported == 0 if no benefit is expected
% and simply using the data from PsychVRHMD('PrepareRender') will provide results with
% the same quality at a lower computational cost.
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
% PsychVRHMD('SetupRenderingParameters', hmd [, basicTask='Tracked3DVR'][, basicRequirements][, basicQuality=0][, fov=[HMDRecommended]][, pixelsPerDisplay=1])
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
% PsychVRHMD('SetBasicQuality', hmd, basicQuality);
% - Set basic level of quality vs. required GPU performance.
%
%
% oldSetting = PsychVRHMD('SetFastResponse', hmd [, enable]);
% - Return old setting for 'FastResponse' mode in 'oldSetting',
% optionally disable or enable the mode via specifying the 'enable'
% parameter as 0 or > 0. Please note that if you want to use 'FastResponse',
% you must request and thereby enable it at the beginning of a session, as
% the driver must do some neccessary setup prep work at startup of the HMD.
% Once it was initially enabled, you can switch the setting at runtime with
% this function.
%
% Some drivers may implement different strategies for 'FastResponse', selectable
% via different settings for the 'enable' flag here. E.g., the Oculus Rift driver
% support three different methods 1, 2 and 3 at the moment.
%
%
% oldSetting = PsychVRHMD('SetTimeWarp', hmd [, enable]);
% - Return old setting for 'TimeWarp' mode in 'oldSetting',
% optionally enable or disable the mode via specifying the 'enable'
% parameter as 1 or 0. Please note that if you want to use 'TimeWarp',
% you must request and thereby enable it at the beginning of a session, as
% the driver must do some neccessary setup prep work at startup of the HMD.
% Once it was initially enabled, you can switch the setting at runtime with
% this function.
%
%
% oldSetting = PsychVRHMD('SetLowPersistence', hmd [, enable]);
% - Return old setting for 'LowPersistence' mode in 'oldSetting',
% optionally enable or disable the mode via specifying the 'enable'
% parameter as 1 or 0.
%
%
% PsychVRHMD('SetHSWDisplayDismiss', hmd [, dismissTypes=1+2+4]);
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
% [bufferSize, imagingFlags, stereoMode] = PsychVRHMD('GetClientRenderingParameters', hmd);
% - Retrieve recommended size in pixels 'bufferSize' = [width, height] of the client
% renderbuffer for each eye for rendering to the HMD. Returns parameters
% previously computed by PsychVRHMD('SetupRenderingParameters', hmd).
%
% Also returns 'imagingFlags', the required imaging mode flags for setup of
% the Screen imaging pipeline. Also returns the needed 'stereoMode' for the
% pipeline.
%
% This function is usually called by PsychImaging(), you don't need to deal
% with it.
%
%
% needPanelFitter = PsychVRHMD('GetPanelFitterParameters', hmd);
% - 'needPanelFitter' is 1 if a custom panel fitter tasks is needed, and 'bufferSize'
% from the PsychVRHMD('GetClientRenderingParameters', hmd); defines the size of the
% clientRect for the onscreen window. 'needPanelFitter' is 0 if no panel fitter is
% needed.
%
%
% [winRect, ovrfbOverrideRect, ovrSpecialFlags] = PsychVRHMD('OpenWindowSetup', hmd, screenid, winRect, ovrfbOverrideRect, ovrSpecialFlags);
% - Compute special override parameters for given input/output arguments, as needed
% for a specific HMD. Take other preparatory steps as needed, immediately before the
% Screen('OpenWindow') command executes. This is called as part of PsychImaging('OpenWindow'),
% with the user provided hmd, screenid, winRect etc.
%
%
% isOutput = PsychVRHMD('IsHMDOutput', hmd, scanout);
% - Returns 1 (true) if 'scanout' describes the video output to which the
% HMD 'hmd' is connected. 'scanout' is a struct returned by the Screen
% function Screen('ConfigureDisplay', 'Scanout', screenid, outputid);
% This allows probing video outputs to find the one which feeds the HMD.
%
%
% [headToEyeShiftv, headToEyeShiftMatrix] = PsychVRHMD('GetEyeShiftVector', hmd, eye);
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
% [projL, projR] = PsychVRHMD('GetStaticRenderParameters', hmd [, clipNear=0.01][, clipFar=10000.0]);
% - Retrieve parameters needed to setup the intrinsic parameters of the virtual
% camera for scene rendering.
%
% 'clipNear' Optional near clipping plane for OpenGL. Defaults to 0.01.
% 'clipFar' Optional far clipping plane for OpenGL. Defaults to 10000.0.
%
% Return arguments:
%
% 'projL' is the 4x4 OpenGL projection matrix for the left eye rendering.
% 'projR' is the 4x4 OpenGL projection matrix for the right eye rendering.
% Please note that projL and projR are usually identical for typical rendering
% scenarios.
%
%
% PsychVRHMD('Start', hmd);
% - Start live operations of the 'hmd', e.g., head tracking.
%
%
% PsychVRHMD('Stop', hmd);
% - Stop live operations of the 'hmd', e.g., head tracking.
%
%

% History:
% 23-Sep-2015  mk  Written.

% Global GL handle for access to OpenGL constants needed in setup:
global GL;

% Global OVR handle for access to VR runtime constants:
global OVR;

if nargin < 1 || isempty(cmd)
  help PsychVRHMD;
  return;
end

% Autodetect first connected HMD and open a connection to it. Open a
% emulated one, if none can be detected. Perform basic setup with
% default configuration, create a proper PsychImaging task.
if strcmpi(cmd, 'AutoSetupHMD')
  if length(varargin) >= 1 && ~isempty(varargin{1})
    basicTask = varargin{1};
  else
    basicTask = 'Tracked3DVR';
  end

  if length(varargin) >= 2 && ~isempty(varargin{2})
    basicRequirements = varargin{2};
  else
    basicRequirements = '';
  end

  if length(varargin) >= 3 && ~isempty(varargin{3})
    basicQuality = varargin{3};
  else
    basicQuality = 0;
  end

  if length(varargin) >= 5 && ~isempty(varargin{5})
    deviceIndex = varargin{5};
  else
    deviceIndex = [];
  end

  if length(varargin) >= 4 && ~isempty(varargin{4})
    vendor = varargin{4};
    if strcmpi(vendor, 'Oculus')
      if exist('PsychOculusVR1', 'file')
        hmd = PsychOculusVR1('AutoSetupHMD', basicTask, basicRequirements, basicQuality, deviceIndex);
      else
        hmd = [];
      end

      if isempty(hmd)
        hmd = PsychOculusVR('AutoSetupHMD', basicTask, basicRequirements, basicQuality, deviceIndex);
      end

      % Return the handle:
      varargout{1} = hmd;
      evalin('caller','global OVR');

      return;
    end

    if strcmpi(vendor, 'OpenHMD')
      if exist('PsychOpenHMDVR', 'file')
        hmd = PsychOpenHMDVR('AutoSetupHMD', basicTask, basicRequirements, basicQuality, deviceIndex);
      else
        hmd = [];
      end

      % Return the handle:
      varargout{1} = hmd;
      evalin('caller','global OVR');

      return;
    end

    error('AutoSetupHMD: Invalid ''vendor'' requested. This vendor is not supported.');
  end

  % Probe sequence:
  hmd = [];

  % Oculus runtime v1.11+ supported and online? At least one real HMD connected?
  if exist('PsychOculusVR1', 'file') && PsychOculusVR1('Supported') && PsychOculusVR1('GetCount') > 0
    % Yes. Use that one. This will also inject a proper PsychImaging task
    % for setup of the imaging pipeline:
    hmd = PsychOculusVR1('AutoSetupHMD', basicTask, basicRequirements, basicQuality, deviceIndex);

    % Return the handle:
    varargout{1} = hmd;
    evalin('caller','global OVR');
    return;
  end

  % Oculus sdk/runtime v0.5 supported and online? At least one real HMD connected?
  if PsychOculusVR('Supported') && PsychOculusVR('GetCount') > 0
    % Yes. Use that one. This will also inject a proper PsychImaging task
    % for setup of the imaging pipeline:
    hmd = PsychOculusVR('AutoSetupHMD', basicTask, basicRequirements, basicQuality, deviceIndex);

    % Return the handle:
    varargout{1} = hmd;
    evalin('caller','global OVR');
    return;
  end

  % OpenHMD VR supported and online? At least one real HMD connected?
  if PsychOpenHMDVR('Supported') && PsychOpenHMDVR('GetCount') > 0
    % Yes. Use that one. This will also inject a proper PsychImaging task
    % for setup of the imaging pipeline:
    hmd = PsychOpenHMDVR('AutoSetupHMD', basicTask, basicRequirements, basicQuality, deviceIndex);

    % Return the handle:
    varargout{1} = hmd;
    evalin('caller','global OVR');
    return;
  end

  % Add probe and autosetup calls for other HMD vendors here...

  % No success with finding any real supported HMD so far. Try to find a driver
  % that at least supports an emulated HMD for very basic testing:
  if PsychOculusVR('Supported')
    hmd = PsychOculusVR('AutoSetupHMD', basicTask, basicRequirements, basicQuality, deviceIndex);
    varargout{1} = hmd;
    evalin('caller','global OVR');
    return;
  end

  if PsychOpenHMDVR('Supported')
    hmd = PsychOpenHMDVR('AutoSetupHMD', basicTask, basicRequirements, basicQuality, deviceIndex);
    varargout{1} = hmd;
    evalin('caller','global OVR');
    return;
  end

  % Add probe for other emulated HMD drivers here ...

  % If we reach this point then it is game over:
  fprintf('PsychVRHMD:AutoSetupHMD: Could not autosetup any HMDs, real or emulated, for any HMD vendor. Game over!\n');

  % Return an empty handle to signal lack of VR HMD support to caller,
  % so caller can cope with it somehow:
  varargout{1} = [];
  return;
end

% If the cmd could not get dispatched by us, funnel it to the
% vendor specific driver:
% 'cmd' so far not dispatched? Let's assume it is a command
% meant for the HMD specific driver:
if (length(varargin) >= 1) && isstruct(varargin{1})
  myhmd = varargin{1};
  [ varargout{1:nargout} ] = myhmd.driver(cmd, myhmd, varargin{2:end});
end

end
