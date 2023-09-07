function EyelinkGetTrackerImageDemo
% This shows you how to register a callback m-file that can display the tracker's eye image in PTB

PsychDefaultSetup(1);

try
    % Disable key output to Matlab window:
    ListenChar(-1);

    % Open a window for eye image display:
    w = Screen('OpenWindow', max(Screen('Screens')), [255 255 0], [0 0 800 600]);

    % Initialize 'el' eyelink struct with proper defaults for output to
    % window 'w':
    el=EyelinkInitDefaults(w);

    % Initialize Eyelink connection (real or dummy). The flag '1' requests
    % use of callback function and eye camera image display:
    if ~EyelinkInit([], 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;
        return;
    end

    % This would display additional debug output, if commented in:
    %    Eyelink('Verbosity',10);

    % Run synthetic test suite for the fun of it, if commented in:
    Eyelink('TestSuite');

    % Perform tracker setup: The flag 1 requests interactive setup with
    % video display:
    result = Eyelink('StartSetup',1);

    % Perform drift correction: The special flags 1,1,1 request
    % interactive correction with video display:
    % You have to hit esc before return.
    result = Eyelink('DriftCorrStart',30,30,1,1,1);

    % Done.
catch
    % In case of error, be tidy:
    cleanup;
end

% Shutdown everything at regular end:
cleanup;

end

% Cleanup routine:
function cleanup
    % Shutdown Eyelink:
    Eyelink('Shutdown');

    % Close window:
    sca;

    % Restore keyboard output to Matlab:
    ListenChar(0);
end
