function res = VRInputStuffTest
% VRInputStuffTest - Test input functionality related to VR devices.

global OVR;

% Setup unified keymapping and unit color range:
PsychDefaultSetup(2);

% Select screen with highest id as Oculus output display:
screenid = max(Screen('Screens'));

% Open our fullscreen onscreen window with black background clear color:
PsychImaging('PrepareConfiguration');

hmd = PsychVRHMD('AutoSetupHMD');
if isempty(hmd)
    fprintf('No VR HMDs connected. Game over!\n');
    return;
end

PsychDebugWindowConfiguration;
[win, rect] = PsychImaging('OpenWindow', screenid);
ifi = Screen('GetFlipInterval', win);
hmdinfo = PsychVRHMD('GetInfo', hmd);

clc;

if strcmpi(hmdinfo.subtype, 'Oculus-1')
    fprintf('Properties of our subject:\n\n');
    fprintf('Player height: %f\n', PsychOculusVR1('FloatProperty', hmd, OVR.KEY_PLAYER_HEIGHT));
    fprintf('Player eye height: %f\n', PsychOculusVR1('FloatProperty', hmd, OVR.KEY_EYE_HEIGHT));
    fprintf('Player neck-to-eye: %f\n', PsychOculusVR1('FloatsProperty', hmd, OVR.KEY_NECK_TO_EYE_DISTANCE));
    fprintf('Player eye-to-nose: %f\n', PsychOculusVR1('FloatsProperty', hmd, OVR.KEY_EYE_TO_NOSE_DISTANCE));
    % Suppress error output during username query, as these properties are often undefined and
    % would spill pointless error messages:
    oldverbosity = PsychOculusVR1('Verbosity', 0);
    fprintf('User name: %s\n', PsychOculusVR1('StringProperty', hmd, OVR.KEY_USER, 'Hans'));
    fprintf('Player name: %s\n', PsychOculusVR1('StringProperty', hmd, OVR.KEY_NAME, 'Mueller'));
    PsychOculusVR1('Verbosity', oldverbosity);
    fprintf('Player gender: %s\n', PsychOculusVR1('StringProperty', hmd, OVR.KEY_GENDER, OVR.KEY_DEFAULT_GENDER));
    fprintf('\n\n');
end

if hmdinfo.VRControllersSupported
    fprintf('\n\nConnected controllers:\n\n');
    controllerTypes = PsychVRHMD('Controllers', hmd);
    if bitand(controllerTypes, OVR.ControllerType_Remote)
        fprintf('Remote control connected.\n');
    end
    if bitand(controllerTypes, OVR.ControllerType_XBox)
        fprintf('XBox controller connected.\n');
    end
    if bitand(controllerTypes, OVR.ControllerType_LTouch)
        fprintf('Left hand controller connected.\n');
    end
    if bitand(controllerTypes, OVR.ControllerType_RTouch)
        fprintf('Right hand controller connected.\n');
    end
    fprintf('\n\n');
end

fprintf('VR area boundaries (if configured):\n\n');

[isVisible, playboundsxyz, outerboundsxyz] = PsychVRHMD('VRAreaBoundary', hmd)

fprintf('\n\nPress any key to continue.\n');
KbStrokeWait(-1);
pulseEnd = 0;

while 1
    WaitSecs('YieldSecs', 0.1);
    clc;

    % Query and display all input state:
    istate = PsychVRHMD('GetInputState', hmd, OVR.ControllerType_Active);
    
    if ~istate.Valid
        fprintf('No VR controller input available, possibly due to VR input focus loss.\n');
        fprintf('Make sure our application has input focus and possibly that the HMD is on the users head.\n');
        continue;
    end
  
    fprintf('Press BACK button on remote control or other controllers, or backspace key, to finish.\n\n');
    disp(istate);

    if istate.Buttons(OVR.Button_A)
        fprintf('Button_A ');
    end

    if istate.Buttons(OVR.Button_B)
        fprintf('Button_B ');
    end

    if istate.Buttons(OVR.Button_X)
        fprintf('Button_X ');
    end

    if istate.Buttons(OVR.Button_Y)
        fprintf('Button_Y ');
    end

    if istate.Buttons(OVR.Button_RThumb)
        fprintf('Button_RThumb ');
    end

    if istate.Buttons(OVR.Button_RShoulder)
        fprintf('Button_RShoulder ');
    end

    if istate.Buttons(OVR.Button_LThumb)
        fprintf('Button_LThumb ');
    end

    if istate.Buttons(OVR.Button_LShoulder)
        fprintf('Button_LShoulder ');
    end

    if istate.Buttons(OVR.Button_Up)
        fprintf('Button_Up ');
    end

    if istate.Buttons(OVR.Button_Down)
        fprintf('Button_Down ');
    end

    if istate.Buttons(OVR.Button_Left)
        fprintf('Button_Left ');
    end

    if istate.Buttons(OVR.Button_Right)
        fprintf('Button_Right ');
    end

    if istate.Buttons(OVR.Button_Enter)
        fprintf('Button_Enter ');
    end

    if istate.Buttons(OVR.Button_Back)
        fprintf('Button_Back ');
    end

    if istate.Buttons(OVR.Button_VolUp)
        fprintf('Button_VolUp ');
    end

    if istate.Buttons(OVR.Button_VolDown)
        fprintf('Button_VolDown ');
    end

    if istate.Buttons(OVR.Button_Home)
        fprintf('Button_Home ');
    end

    if any(istate.Buttons(OVR.Button_RMask))
        fprintf('Button_RMask ');
    end

    if any(istate.Buttons(OVR.Button_LMask))
        fprintf('Button_LMask ');
    end

    if any(istate.Buttons(OVR.Button_Private))
        fprintf('Button_Private ');
    end

    [~, ~, keycode] = KbCheck(-1);
    if istate.Buttons(OVR.Button_Back) || keycode(KbName('BackSpace'))
        break;
    end

    fprintf('\n');

    % Oculus VR 1.11+ controlled VR setup?
    if strcmpi(hmdinfo.subtype, 'Oculus-1')
        % Get info about collisions of tracked devices with the VR play area boundaries:
        [isTriggering, closestDistance, closestPointxyz, surfaceNormal] = PsychOculusVR1('TestVRBoundary', hmd, OVR.TrackedDevice_All, 0);
        if isTriggering
            fprintf('VR collision: closestDistance %f m, point [%f,%f,%f], normal [%f,%f,%f]\n', closestDistance, closestPointxyz(1), closestPointxyz(2), ...
                    closestPointxyz(3), surfaceNormal(1), surfaceNormal(2), surfaceNormal(3));
        end

        % Get info about collision of test point with the VR play area boundaries:
        mypoint = [1, 2, 3];
        [isTriggering, closestDistance, closestPointxyz, surfaceNormal] = PsychOculusVR1('TestVRBoundaryPoint', hmd, mypoint, 0);
        if isTriggering
            fprintf('VR mypoint collision: closestDistance %f m, point [%f,%f,%f], normal [%f,%f,%f]\n', closestDistance, closestPointxyz(1), closestPointxyz(2), ...
                    closestPointxyz(3), surfaceNormal(1), surfaceNormal(2), surfaceNormal(3));
        end
    end

    if istate.Buttons(OVR.Button_A)
        if GetSecs < pulseEnd
            % Pressed A before end of previously initiated pulse. Abort current pulse:
            pulseEnd = PsychVRHMD('HapticPulse', hmd, OVR.ControllerType_XBox, 0, 0);
        else
            % Initiate new pulse: 2.5 (max) seconds, 100% frequency, 0.8 amplitude:
            pulseEnd = PsychVRHMD('HapticPulse', hmd, OVR.ControllerType_XBox, [], 1.0, 0.8)
            % Debounce button:
            WaitSecs('YieldSecs', 0.5);
        end
    end
end

sca;
