% Short MATLAB example program that uses the Eyelink and Psychophysics
% Toolboxes to create a real-time gaze-dependent display.
% This is the original example as shown in the EyelinkToolbox article in BRMIC
% Cornelissen, Peters and Palmer 2002).
% Note: this version is not PC compatible. Use eyelinkexample instead.

% STEP 1
% Initialization of the connection with the Eyelink Gazetracker.
% exit program if this fails.
if (EYELINK('initialize') ~= 0)
	return;
end;

% STEP 2
% Open a graphics window on the main screen
% using the PsychToolbox's SCREEN function.

screennr = 0; % use main screen
[window, screenRect]=SCREEN(screennr,'OpenWindow', 0);
white=WhiteIndex(window);
black=BlackIndex(window);

% STEP 3
% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).
el=initeyelinkdefaults;
% make sure that we get gaze data from the Eyelink
EYELINK('command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
% open file to record data to
EYELINK('openfile', 'demo.edf');

% STEP 4
% Calibrate the eye tracker using the standard calibration routines
EYELINK('trackersetup');
% do a final check of calibration using driftcorrection
EYELINK('dodriftcorrect');
    
% STEP 5
% start recording eye position
EYELINK('startrecording');
% record a few samples before we actually start displaying
waitsecs(0.1);
% mark zero-plot time in data file
EYELINK('message', 'SYNCTIME');
stopkey=KbName('space');
eye_used = -1;

% STEP 6
% show gaze-dependent display
while 1 % loop till error or space bar is pressed
	% Check recording status, stop display if error
	error=EYELINK('checkrecording');
	if(error~=0)
		break;
	end
	% check for keyboard press
	[keyIsDown,secs,keyCode] = KbCheck;
	% if spacebar was pressed stop display
	if keyCode(stopkey)
		break;
	end
	% check for presence of a new sample update
	if EYELINK( 'newfloatsampleavailable') > 0
		% get the sample in the form of an event structure
		evt = EYELINK( 'newestfloatsample');
		if eye_used ~= -1 % do we know which eye to use yet?
			% if we do, get current gaze position from sample
			x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
			y = evt.gy(eye_used+1);
			% do we have valid data and is the pupil visible?
			if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
				% if data is valid, draw a circle on the screen at current gaze position
				% using PsychToolbox's SCREEN function
				gazeRect=[ x-7 y-7 x+8 y+8];
				SCREEN(window, 'FrameOval', white,gazeRect,6,6);
			else
				% if data is invalid (e.g. during a blink), clear display
				SCREEN(window, 'FillRect',black);
			end
		else % if we don't, first find eye that's being tracked
			eye_used = EYELINK( 'eyeavailable'); % get eye that's tracked
			if eye_used == el.BINOCULAR; % if both eyes are tracked
				eye_used = el.LEFT_EYE; % use left eye
			end
		end
	end % if sample available
end % main loop
% wait a while to record a few more samples
waitsecs(0.1);

% STEP 7
% finish up: stop recording eye-movements,
% close graphics window, close data file and shut down tracker
EYELINK('stoprecording');
SCREEN(window,'close');
EYELINK('closefile');
EYELINK('shutdown');
