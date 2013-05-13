% Short MATLAB example program that uses the Eyelink and Psychophysics
% Toolboxes to create a real-time gaze-dependent display.

clear all;
fprintf('My OSX Eyelink GetKeys Test\n\n\t');
input('Hit the return key to continue.','s');
fprintf('Thanks.\n');


% Initialization of the connection with the Eyelink Gazetracker.
% exit program if this fails.
el=initeyelinkdefaultsOSX;

if (Eyelink('Initialize') ~= 0)	return;end;


white = [255,255,255];
red = [128 0 0];
green = [0 128 0];
black = [0,0,0];
missing = -32768;
IN_TARGET_MODE=8;

% Open a graphics window on the main Screen
[window, ScreenRect]=Screen('OpenWindow', 0, 0, [], 32 ,1);
[offWindow, offScreenRect]=Screen('OpenOffscreenWindow',window);

%Screen('FillRect',offWindow,white,offScreenRect);
%Screen('CopyWindow',offWindow,window);

%Screen('DrawText',offWindow,'wiebel', 200,200,black);
Screen('FillRect',window,black);
%Screen('FillRect',window,red,[100 100 200 200]);
%Screen('FillRect',window,red,[150 150 250 250]);
%Screen('DrawText',window,'wiebel', 200,200,black);

% set eyelink display, used for calibration and recording
Eyelink('Command', 'screen_pixel_coords = 0 0 1280 1024');
Eyelink('Message', 'DISPLAY_COORDS 0 0 1280 1024');
Eyelink('Command', 'calibration_type = HV9');

% set parser saccade thresholds (conservative settings) */
% 0 = standard sensitivity */
Eyelink('Command', 'select_parser_configuration 0');

%Eyelink('command', 'select_parser_configuration 0');
Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE');
Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS');
Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK');
Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');

% Calibrate
Eyelink('Command', 'heuristic_filter = ON');
Eyelink('StartSetup');		% start setup mode
Eyelink('WaitForModeReady', 500);  % time for mode change

[x,y,z] = Eyelink('Test')
% wait for key to send to eyelink, until calibration started
%while 1
    result= 1;
	times = 0;
    
    [key, el]=getkeyforeyelinkOSX(el);
    if key==el.TERMINATE_KEY
        break;
    end
    if key~=el.JUNK_KEY && key ~= 0
        Eyelink('SendKeyButton', double(key), 0, 10);         
        if key == 'c'
           break;
        end
    end
%end
pause(1.0);
if key == 'c'  
    result= 1;
	times = 0; 
	while bitand(Eyelink('CurrentMode'), IN_TARGET_MODE) && (times < 1)
		if(~Eyelink('IsConnected'))
			break;
		end 
        [key, el]=getkeyforeyelinkOSX(el);
        if key==el.TERMINATE_KEY
            break;
        end
        if key~=el.JUNK_KEY && key ~= 0
            Eyelink('SendKeyButton', double(key), 0, 10);         
	%        fprintf('Eyelink key %s (%.1f)\n', key, GetSecs-start);
        end       
        
		% HANDLE TARGET CHANGES
		[result, x, y]= Eyelink('TargetCheck');
	    gazeRect1=[x-7 y-7 x+8 y+8];
	    gazeRect2=[x-4 y-4 x+4 y+4];
        if (x~=missing && y~=missing)
           Screen('FillRect',window,black);
           pause(0.05)
%		   Screen('FillOval',window,white,gazeRect1);            
           pause(0.05)
           Screen('FillOval',window,red,gazeRect2);            
           pause(0.05)
           SetMouse(x,y);
           beep;
           pause(0.1);
           % blocking call GetChar, waits for input
           Eyelink('SendKeyButton', double(GetChar), 0, 10);
           
        end
        if result == 0 
           times = times + 1;
        end
        pause(0.1);
        %disp('calibrating')
	end % while IN_TARGET_MODE   
end    
 
[x,y,z] = Eyelink('Test')

% open file to record data to
Eyelink('OpenFile', 'osxtest.edf');

% start recording eye position
Eyelink('StartRecording');

% record a few samples before we actually start displaying
pause(0.2);

% mark zero-plot time in data file
Eyelink('Message', 'DISPLAY_ON');
Eyelink('Message', 'SYNCTIME');

eye_used = 1; % eye (0 is left, 1 is right)

%show gaze-dependent display
for I = 1:30 % loop till error or space bar is pressed
	% Check recording status, stop display if error
    pause(0.05);
	error=Eyelink('CheckRecording');
	if(error~=0)
		break;
	end
	if(~Eyelink('IsConnected'))
		break;
	end     

	% check for presence of a new sample update
	if Eyelink('NewFloatSampleAvailable') > 0
		% get the sample in the form of an event structure
		evt = Eyelink('NewestFloatSample');
		% get current gaze position from sample
		x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
		y = evt.gy(eye_used+1);
		% do we have valid data and is the pupil visible?
		if x~=-missing && y~=-missing && evt.pa(eye_used+1)>0
			% if data is valid, draw a rect on the Screen at current gaze position
			% using PsychToolbox's Screen function
			gazeRect=[x-7 y-7 x+8 y+8];
			Screen('FillOval',window,red,gazeRect);
%			Screen('Flip',window);          
            SetMouse(round(x),round(y));
		else
			% if data is invalid (e.g. during a blink), clear display
			Screen('FillRect',window);
%			Screen('Flip',window); 
        end
    end % if sample available
end % main loop

% wait a while to record a few more samples
pause(0.5);

[x,y,z] = Eyelink('Test')
% close graphics window, close data file and shut down tracker
Screen('Close',window);
Screen('Close',offWindow);
Eyelink('StopRecording');
Eyelink('SetOfflineMode');
Eyelink('CloseFile');
Eyelink('ReceiveFile');
Eyelink('Shutdown');
!/Applications/Eyelink/EDF_Access_API/Example/edf2asc_static_osx.exe osxtest.edf
