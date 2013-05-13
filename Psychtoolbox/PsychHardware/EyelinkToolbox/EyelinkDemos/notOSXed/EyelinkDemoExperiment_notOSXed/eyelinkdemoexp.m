function eyelinkdemoexp

% demo for some simple gaze recording experiments using eyelinktoolbox
% created after SAMPLE EXPERIMENT FOR MACINTOSH, by Dave Stampe , 27 June '97
% History
% 15-05-01	fwc created first version
% 22-05-01	fwc	debugging
% 25-05-01	fwc	fixed small bug, different pict for simpletrial demo
% 02-06-01	fwc removed use of global el, as suggest by John Palmer.
%				JP's e5a.m demo also showed me that more functions
%				could be included in a single file.
%				functions setupwindow and runexperiment are now part of this file
% 02-06-01	fwc test for presence of default sounds using testeyelinksounds
% 18-10-02	fwc small cosmetic changes and changes to work with toolbox version 1.4
% 27-11-02	fwc	changed dialog for dummy mode question

clear all;
fprintf('Eyelink Toolbox demo-experiment.\n');
fprintf('Performed on ''%s''.\n\n', datestr(now,0));

edffilename = 'test';
ntrials=2;
initializedummy=0;

window=SCREEN('GetFrontWindow');

if initializedummy~=1
	if eyelink('initialize') ~= 0
		if 1==dummymodedlg
			eyelink('initializedummy');
		else
			return;
		end
	end
else
	eyelink('initializedummy');
end


edffilename=[edffilename '.edf' ];

i = EYELINK( 'openfile', edffilename);

if i~=0
	printf('Cannot create EDF file ''%s'' ', edffilename);
	EYELINK( 'Shutdown');
	return;
end

EYELINK('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox demo-experiment''');

%   SET UP TRACKER CONFIGURATION
eyelink('command', 'calibration_type = HV9');
%	set parser (conservative saccade thresholds)
eyelink('command', 'saccade_velocity_threshold = 35');
eyelink('command', 'saccade_acceleration_threshold = 9500');
%	set EDF file contents
eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
%	set link data (used for gaze cursor)
eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,BUTTON');
eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');


% open window for this experiment
screen=1;
isColor=1;
pixelSize=32;
width=832;
height=624;
hz=75;

[w, rect]=setupwindow(screen, width, height, hz, pixelSize,isColor);
if w<0	% if not connected
	EYELINK('closefile');
	EYELINK('shutdown');
	SCREEN('closeall');
	return;
end

white=WhiteIndex(w);
black=BlackIndex(w);
gray=round((white+black)/2);
SCREEN(w,'FillRect', gray);

eyelink('message', 'DISPLAY_COORDS %d %d %d %d',rect(1),rect(2),rect(3)-1,rect(4)-1);

el=initeyelinkdefaults; % if we do this after we have opened a graphics window
						% this window is used to set defaults

% test if the default sounds are available
% if not, switch off beeping.
err=testeyelinksounds;
if err~=0
	el.targetbeep=0;
end

% Get edf file name.
if 0
	exp.outcode=SCREEN(window,'Dialog','EDF filename?', edffilename );
end
% modify a few of the default settings
el.backgroundcolour = gray;
el.foregroundcolour = white;

status=EYELINK('isconnected');

if status==0	% if not connected
	EYELINK('closefile');
	EYELINK('shutdown');
	SCREEN('closeall');
	return;
end

runexperiment(el, ntrials);	% perform actual experiment

SCREEN('closeall');

EYELINK('closefile');
Waitsecs(1.0); % give tracker time to execute all commands

EYELINK('shutdown');

fprintf('Thank you for your saccades.\n');

% End of program
%--------------------------------------------------------------------


% Functions used in the Program
%--------------------------------------------------------------------


%--------------------------------------------------------------------
% runexperiment

function result=runexperiment(el, ntrials)

% trial loop for demo-experiment

result=el.ABORT_EXPT;

error=dotrackersetup(el, el.ENTER_KEY);	% PERFORM CAMERA SETUP, CALIBRATION */
%EYELINK('trackersetup'); % can be used to, but doesn't use the el defined defaults!

if error==el.TERMINATE_KEY
	return;
end;

for trial=1:ntrials
	if EYELINK('isconnected')==el.notconnected   % drop out if link closed
		result=el.ABORT_EXPT;
		return;
	end;

	if 1 && trial==1
		trialresult = simpletrialdemo(el, trial); % RUN THE TRIAL
	else
		trialresult = realtimedemo(el, trial); % RUN THE TRIAL
	end

	switch trialresult      % REPORT ANY ERRORS
	  case el.ABORT_EXPT,        % handle experiment abort or disconnect
		EYELINK('message', 'EXPERIMENT ABORTED' );
		result=el.ABORT_EXPT;
		return;
	  case el.REPEAT_TRIAL,	  % trial restart requested
		EYELINK('message', 'TRIAL REPEATED' );
	    trial=trial-1;
	  case el.SKIP_TRIAL,	   % skip trial
		EYELINK('message', 'TRIAL ABORTED' );
	  case el.TRIAL_OK,
		EYELINK('message', 'TRIAL OK' );
	  otherwise,
		EYELINK('message', 'TRIAL ERROR' );
	end

end  % end of trial loop

result=0;
return;



%--------------------------------------------------------------------
% setupwindow

function [w, rect]=setupwindow(whichScreen, width, height, hz, pixelSize,isColor)

% open a window of particular size and refresh rate
%open screen
rect=[0 0 0 0];
w=-1;
pixelSizes=SCREEN(whichScreen,'PixelSizes');
if max(pixelSizes)<pixelSize
	fprintf('Sorry, I need a screen that supports %d-bit pixelSize.\n', pixelSize);
	return;
end
if ~strcmp(computer,'MAC2')
	res.width=width;
	res.height=height;
	res.hz=hz;
	res.pixelSize=pixelSize;
else 
    res=NearestResolution(whichScreen,width, height, hz, pixelSize);
end;
[w,rect]=SCREEN(whichScreen,'OpenWindow',0,[],res);
[oldPixelSize,oldIsColor,pages]=SCREEN(w,'PixelSize',pixelSize,isColor);
oldBoolean=SCREEN(w,'Preference','SetClutDriverWaitsForBlanking',1);

white=WhiteIndex(w);
black=BlackIndex(w);
gray=round((white+black)/2);
SCREEN(w,'FillRect',gray);

%--------------------------------------------------------------------
