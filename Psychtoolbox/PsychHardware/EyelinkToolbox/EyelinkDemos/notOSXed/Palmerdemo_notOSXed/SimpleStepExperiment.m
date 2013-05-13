function SimpleStepExperiment
% Simple experiment to measure response time using eye movements.
% v10 new improved single condition experiment (like e3, e5)
% Displays fixation followed by left or right stimulus.
% Requires left or right eye movement to corresponding stimulus.  
% Requires file:  setEyelinkDefaults
% John Palmer, last revised 6/6/01
%
% 2/21/01	v1: based on crt1
% 2/26/01	v2: modifed to provide messages for asc2fira data file conversion
% 2/28/01	added code to pass messages with stimulus location
% 3/1/01	revised side message to use time stamps, added breaks
% 3/17/01	added receivefile and changed to update current directory
% 3/20/01	v3: simplified file i/o, record only during trial, sets resolution
% 3/21/01	reworked header and messages, don't save .mat file at all
% 3/27/01	moved all key messages to after DISPLAY_ON
% 4/2/01	added condition vector and message
% 4/6/01	v4: introduced multiple conditions, varies ecc of target
% 4/10/01	fixed bug with lower case condition message
% 4/11/01	increased stimulus duration to 1.0 s
% 4/13/01	v5: cut back to one condition keeping improvements
% 5/2/01	v5a:  changed resolution
% 5/3/01	updated h.pixelsperdegree to 25.5 and other sizes
% 5/25/01	v5t:  modified for eyelinktoolbox v1.0, no receivefile
% 6/4/01	put receivefile back in
% 6/6/01	v10:  updated for eyelinktoolbox 1.1, added calibration
%
% bugs
% pause key only tested for at beginning of trial

% Experiment Header
h.filename = 'e10aa0.mat';				% default filename
h.experiment = 'e10, version dated 6/6/01';
h.hardware = 'EyeLink 2.04, Mac G4, ViewSonic PF790';
h.software1 = 'EyeLinkInterface 1.0, PsychToolbox 2.44';
h.software2 = 'CWPro 6.0, Matlab 5.2.1, MacOS 9.1';
h.author = 'John Palmer';
h.daterun = datestr(now,1);				% save date in dd-mmm-yyyy

h.nconditions = 1;						% number of conditions (not counting side)
h.ngroups = 1;							% groups of trials with breaks between
h.ntrialspergroup = 10;					% number of trials per group
h.ntrials = h.ngroups * h.ntrialspergroup;	% must be mult of 2 (2*ncond)
h.eccentricity = 255;					% horizontal ecc of stimulus in pixels
h.stimsize = 25;						% diameter of stimulus disk in pixels
h.fixationsize = 13;					% diameter of stimulus disk in pixels
h.duration = 1.0;						% duration of stimulus display in sec
h.fixationtime = 1.0;					% duration of fixation display in sec
h.intertrialinterval = 1.0;				% in sec
h.initialseed = ClockRandSeed;			% set seed w/ clock and put in header
h.distance = '60 cm';					% distance from eye to ctr of screen
hres = 832; vres = 624; tres = 75;		% define screen resolution
h.screenresolution = ...
	strcat(num2str(hres),' byÊ',num2str(vres),' pixels atÊ',num2str(tres),' Hz');
h.screensize = '35 by 25.5 cm';
h.pixelsperdegree = 25.5;				% based on central 1 degree region
pausekey = 'space';						% key to pause experiment
pausekeycode = KbName(pausekey);

% set up condition and side data structures
cvalue = [0,0];							%not used in this experiment
h.conditionlist = strcat(num2str(cvalue(1)),',',num2str(cvalue(2)),' not used');

npc = round(h.ntrials/(2*h.nconditions));	% Counterbalance side of stimulus
if npc*2*h.nconditions ~= h.ntrials
	error('COUNTERBALANCING ERROR:  bad number of conditions'); end;

npc = round(h.ntrials/2);	
side = [ones(1,npc), repmat(2,1,npc)];	% fill w/ 1 && 2s
npc = round(h.ntrials/h.nconditions/2);	
condition = [ones(1,npc),ones(1,npc)];	% hard code for more conditions

rstate = rand('state');					% get random func state
side = MyPermute(side);					% permute elements of vector
rand('state',rstate);					% reset random number generater
condition = MyPermute(condition);		% permute in same way

% Initialization of screen, colors, sounds
oldRes=SetResolution(0,hres,vres,tres);	% fix resolution
[window,screenRect]=SCREEN(0,'OpenWindow');	% use default screen 0
white=WhiteIndex(window);				% set up colors
black=BlackIndex(window);				
[xcenter,ycenter] = getCenter(window);	% calculate center of screen
beep = MakeBeep(500, 0.2);				% make 500 hz, 0.2 s tone used for feedback
SND('Open');							% initialize sound

% Do filename i/o
[h.filename,eyefilename,myerr] = getfilenames(h,window); % prompt for file names
if myerr								% exit on errors in inputs
	shutdown(oldRes); error('Filename Input Fatal Error'); end;
HideCursor;								% hide when done with screen i/o

% Initialize eye tracker
myerr2 = eyelink('initialize');			% initialize eyelink hardware on operater PC
if  myerr2 ~= 0							% quit if error
	shutdown(oldRes); error('Eyelink Initialization Fatal Error'); end;
	
el=initeyelinkdefaults;					% define eyelink constants for this program
el.backgroundcolour = white;
el.foregroundcolour = black;
SCREEN(window, 'FillRect',el.backgroundcolour);  % clear display

seteyelinkdefaults(screenRect);			% set defaults in eyelink hardware
eyelink('openfile', eyefilename);		% open data file on operater PC
sendheader(h);							% send header to eyelink file via messages
	
% Create all of the stimuli in offscreen windows 
fixation = SCREEN(window,'OpenOffscreenWindow');
SCREEN(fixation,'FillOval',black,getStimulusRect(window,h.fixationsize,xcenter));
SCREEN(fixation,'FillOval',white,getStimulusRect(window,3,xcenter));
SCREEN(fixation,'FillOval',black,getStimulusRect(window,1,xcenter));

stimuluslocation(1) = xcenter - h.eccentricity;	
complementlocation(2) = xcenter - h.eccentricity;	% set up complement for errors
leftstim = SCREEN(window,'OpenOffscreenWindow');
SCREEN(leftstim,'FillOval',black,getStimulusRect(window,h.stimsize,stimuluslocation(1)));
SCREEN(leftstim,'FillOval',white,getStimulusRect(window,3,stimuluslocation(1)));
SCREEN(leftstim,'FillOval',black,getStimulusRect(window,1,stimuluslocation(1)));

stimuluslocation(2) = xcenter + h.eccentricity;	
complementlocation(1) = xcenter + h.eccentricity;	% set up complement for errors
rightstim = SCREEN(window,'OpenOffscreenWindow');
SCREEN(rightstim,'FillOval',black,getStimulusRect(window,h.stimsize,stimuluslocation(2)));
SCREEN(rightstim,'FillOval',white,getStimulusRect(window,3,stimuluslocation(2)));
SCREEN(rightstim,'FillOval',black,getStimulusRect(window,1,stimuluslocation(2)));

stimset = [leftstim,rightstim];

% do calibration after all other initialization complete
dotrackersetup(el);						% (disable for quick test)

% Warn beginning of experiment and start loops
takebreak(window,'Starting Experiment, Press any to begin');
trial = 0;

for jj = 1:h.ngroups					% group loop
	if jj > 1
		takebreak(window,'End of block, Press any to continue');
	end;
	
for ii = 1:h.ntrialspergroup			% trial loop
	trial = trial + 1;					% increment trial across groups

	[isPressed,endtime,keycode] = KbCheck;	% check for pause key
	if isPressed
		if keycode(pausekeycode), 
			takebreak(window,'Pause Requested, Press any to continue');
		end
	end

% Start up and trial number at the bottom of operater display
	eyelink('startrecording');			% begin recording at start of trial
	eyelink('command', 'record_status_message ''TRIAL %d''', trial);

% Display the fixation point	
	SCREEN(window,'WaitBlanking');		
	SCREEN('CopyWindow',fixation,window);
	eyelink('message', 'FIXATION_ON');		
	waitsecs(h.fixationtime);

%fprintf('side %3.0f cond %3.0f',side(trial), cvalue(condition(trial))
	
% Display the stimulus
	SCREEN(window,'WaitBlanking');		
	SCREEN('CopyWindow',stimset(side(trial)),window);
	eyelink('message', 'FIXATION_OFF');		
	eyelink('message', 'DISPLAY_ON');	%message for RT recording in analysis
	eyelink('message', 'SYNCTIME');		%zero-plot time for EDFVIEW

% Send messages describing the trial
	eyelink('message', 'TRIAL ID %d', trial);	
	eyelink('message', 'CONDITION ID %d', condition(trial));	
	eyelink('message', 'SIDE ID %d', side(trial));	
	eyelink('message', 'XFIXATION ID %d', xcenter);	
	eyelink('message', 'YFIXATION ID %d', ycenter);	
	eyelink('message', 'XCORRECT ID %d', stimuluslocation(side(trial)));	
	eyelink('message', 'YCORRECT ID %d', ycenter);	
	eyelink('message', 'XERROR ID %d', complementlocation(side(trial)));
	eyelink('message', 'YERROR ID %d', ycenter);	
	waitsecs(h.duration);

% Process keypress and wrap up trial
	SCREEN(window,'WaitBlanking');		
	SCREEN(window,'FillRect',white);	% erase screen
	eyelink('message', 'DISPLAY_OFF');	% message for end main display
	eyelink('message', 'TIMEOUT');		% message for end of trial
	eyelink('message', 'TRIAL_OK');		% message for successful trial
	waitsecs(0.01);						% give time to finish w/ messages
	eyelink('stoprecording');			% stop recording at end of trial
	waitsecs(h.intertrialinterval);
	
end;	% of trial loop
end;	% of group loop

% Exit
fprintf('Header for This Experiment\n\n');
disp(h); 								% Print out header
eyelink('closefile');					% close eye data file
status = eyelink('receivefile',eyefilename,eyefilename);	% copy data file
if status < 0, fprintf('Error in receiveing file!\n'); end;
shutdown(oldRes);						% End of main program

% End of program
%--------------------------------------------------------------------

% Functions used in the Program
%--------------------------------------------------------------------
function [filename,eyefilename,myerr] = getfilenames(h,window)
% getfilenames returns file names for main data and eye data.
% strict error checking on file names

myerr = 0; eyefilename = '';

% get main file name
[shortfilename, filepath] = uiputfile(h.filename, 'Data file name?');	
cd(filepath);							% update current path
		
if isempty(shortfilename) || 0 == shortfilename	% if empty, exit program
	fprintf('File name required.\n');
	myerr = 1;
	return
end;

if isempty(filepath)
	filepath = pwd; end;				% use default path if empty
	
filename = strcat(filepath,shortfilename);

% figure default for eye data file name
place = findstr('.mat',shortfilename);	% based on main datafile
if isempty(place)						% if no .mat, exit program
	fprintf('Suffix .mat required in file name.\n');
	myerr = 1;
	return
	end;

% create name for eye data file
eyefilename = strcat(shortfilename(1:place),'edf');

%--------------------------------------------------------------------
function takebreak(window,s);
% getCenter returns the center location in pixels

white=WhiteIndex(window);
black=BlackIndex(window);
[xcenter,ycenter] = getCenter(window);
SCREEN(window,'DrawText',s,xcenter-100,ycenter,black);
while KbCheck; end;					% wait until all keys are released.
kbwait;
SCREEN(window,'WaitBlanking');		
SCREEN(window,'FillRect',white);	% erase screen
waitsecs(0.5);
while KbCheck; end;					% wait until all keys are released.

%--------------------------------------------------------------------
function [xcenter,ycenter] = getCenter(mywindow)
% getCenter returns the center location in pixels
windowSize = Screen(mywindow,'Rect');
xcenter = round(windowSize(1,3)/2);
ycenter = round(windowSize(1,4)/2);

%--------------------------------------------------------------------
function rect = getStimulusRect(mywindow,size,myoffset)
% getStimulusRect returns rect specifying stimulus location rectangle.
% new version assumes center already added into x offset

rectanglesize = [ 0,0,size,size];
[xcenter,ycenter] = getCenter(mywindow);
xoffset = round(-size/2+myoffset);		%xcenter removed
yoffset = round(ycenter-size/2);
offset = [xoffset,yoffset,xoffset,yoffset];
rect = rectanglesize + offset;

%--------------------------------------------------------------------
function b = MyPermute(a)
% MyPermute permutes the columns in the matrix a

b = a;
n = size(b,2);							% get number of columns
for i = 1:n								% swap each one once
	r = floor(n*rand+1);
	tmp = b(:,i);
	b(:,i) = b(:,r);
	b(:,r) = tmp;
end;
%--------------------------------------------------------------------
function sendheader(h)
% sends header structure to eyelink as a series of messages
% empty strings not allowed

names = fieldnames(h);					% get list of field names
n = size(names);						% figure number of fields
n = n(1,1);								% must be better way?

for i= 1:n,
	myfield = names(i);
	myvalue = getfield(h,myfield{1});
	if ischar(myvalue)					% handle strings and numbers separately
		m = strcat('HEADER	',myfield{1},' "',myvalue,'"');
		eyelink('message',m);
%		fprintf('%s\n',m);
	elseif isnumeric(myvalue)
		m = strcat('HEADER	',myfield{1},'	',num2str(myvalue));
		eyelink('message',m);
%		fprintf('%s\n',m);
	else
		fprintf('ERROR: Sendheader: must use strings or numeric values');
	end
end;

%--------------------------------------------------------------------
function shutdown(oldRes)
% Shutdown exits this program, closing everything

eyelink('shutdown');
SND('Close');
ShowCursor;
FlushEvents('keyDown');					% discard chars from Event Manager.
SCREEN('CloseAll');
Screen(0,'Resolution',oldRes);

%--------------------------------------------------------------------
% end of functions
