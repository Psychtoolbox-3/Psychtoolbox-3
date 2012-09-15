% Short MATLAB example program that uses the iViewX and Psychophysics
% Toolboxes to create a short eyetracking experiment. This example could
% easily be extended to become some sort of slide viewer.

% STEP 1
% Provide iView with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes.
ivx=initIViewXDefaults;
iViewX('datafile', ivx, 'demo');



% STEP 2
% Initialization of the connection with the iView.
% exit program if this fails.
if iView('initialize', ivx)~=1;
	return;
end;

stopkey=KbName('esc');

% STEP 3
% Open a graphics window on the main screen
% using the PsychToolbox's Screen function.

screennr = 0; % use main screen
[window, screenRect]=SCREEN('OpenWindow', screennr, 0);
gray=GrayIndex(window);
Screen( 'FillRect', window, gray );  % clear screen


% STEP 4
% Calibrate using the standard calibration routines
iViewX('calibration', ivx);

% do a final check of calibration using driftcorrection
iViewX('driftcorrection'); % defaults to center of screen
    
% STEP 5
% start recording eye position
iViewX('startrecording', ivx);

% STEP 6
% display image on screen 

imgname='ichbineinberliner.jpg';
imageArray=ImRead(imgname,'jpg');
Screen('PutImage', window, imageArray);

% mark zero-plot time in data file
iViewX('message', ivx, 'Start of display');

% STEP 7 main measurement loop
while 1 % loop till error or stop key is pressed
	% Check connection status, stop display if error
	
	if iViewX('connection', ivx)~=ivx.isconnected   % Check link often so we don't lock up if tracker lost
		break;
	end;
	
	% check for keyboard press
	[keyIsDown,secs,keyCode] = KbCheck;
	% if stopkey was pressed stop display
	if keyCode(stopkey)
		break;
	end

end % main loop


% STEP 8
% finish up: stop recording eye-movements,
% close graphics window, and shut down tracker
iViewX('stoprecording', ivx);
Screen('close', window);
iViewX('shutdown', ivx);
