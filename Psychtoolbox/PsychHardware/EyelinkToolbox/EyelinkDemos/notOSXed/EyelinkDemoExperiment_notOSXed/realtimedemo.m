function result=realtimedemo(el, trial)

% demo for a simple real-time gaze experiment
% created after SAMPLE EXPERIMENT FOR MACINTOSH, by Dave Stampe , 27 June '97 
% 15-05-01	fwc created first version
% 22-05-01	fwc	debugging with working eyelink.mex
% 02-06-01	fwc removed use of global el, as suggested by John Palmer.
% 03-06-01	fwc	mouse position is now used for the gazewindow
%				if running in dummy mode
	
spacebar=KbName('space');

SCREEN( el.window, 'FillRect', el.backgroundcolour );  % clear display

imgname='himalaya.JPG';
timeouttime=60.0;
buffer=screen(el.window, 'openoffscreenwindow');
screenbuffer=screen(el.window, 'openoffscreenwindow');
SCREEN( screenbuffer, 'FillRect', el.backgroundcolour );  % clear screenbuffer

imageArray=IMREAD(imgname,'jpg');

SCREEN(buffer,'PutImage',imageArray);
gazeRect=[0 0 100 100];
currGazeRect = [0 0 0 0 ];
oldGazeRect=[0 0 0 0 ];
oldx=-1;
oldy=-1;

eye_used = -1;
eyetracked=-1;
% This supplies a title at the bottom of the eyetracker display
EYELINK('command', 'record_status_message ''TRIAL %d''', trial );

% Always send this message before starting to record.
% It marks the start of the trial and also
% contains trial condition data required for analysis.

EYELINK('message', 'TRIALID %d', trial );

success = dodriftcorrection(el,[],[],1,1);

% just for fun, draw a rectangle on eyelink screen
EYELINK('setofflinemode');	% Must be offline to draw to EyeLink screen
rect=[0 0 100 100];
screenRect=SCREEN(el.window, 'Rect');
rect=CenterRect( rect, screenRect);
EYELINK('command', 'draw_box %d %d %d %d 15', rect(1), rect(2), rect(3), rect(4));

% Start data recording to EDF file, BEFORE DISPLAY. */
	
error = eyelink('startrecording');	% record with link data enabled

if error~=0
	result=error;
	return;
end

% You should always start recording 50-100 msec before required */
% otherwise you may lose a few msec of data */
Waitsecs(0.1);

% perform trial
key=1;
while key~= 0
	key=getkeyforeyelink(el);		% dump any pending local keys
end

EYELINK('flushkeybuttons', 0 );   % reset keys and buttons from tracker

% load image from buffer?
SCREEN(el.window, 'Waitblanking');
%SCREEN('Copyel.window',buffer,el.window);

tstart = getsecs;  % record display onset time
EYELINK('message', 'DISPLAY ON');	 % message for RT recording in analysis
EYELINK('message', 'SYNCTIME');	 	 % zero-plot time for EDFVIEW
tstart=getsecs;

while 1  % loop till response
			
	error=EYELINK('checkrecording'); 		% Check recording status */
  	if(error~=0)
		result=error;
		return;
	end
	
	% Note that we don't use 'getkeyforeyelink' unless we don't need fine timing

	[keyIsDown,secs,keyCode] = KbCheck; 
	if keyCode(el.modifierkey) && keyCode(el.quitkey)  % alternative for cmd . which quits program directly
		fprintf( '\nUser requested break.\n');
		result=el.ABORT_EXPT;
		return;
	end

	if keyCode(el.quitkey)
		SCREEN( el.window, 'FillRect', el.backgroundcolour );  % hide display
 		waitsecs(0.1);       % record additional 100 msec of data
  		EYELINK('stoprecording');
 		result=el.SKIP_TRIAL;   % return code to skip trial  		
 		return;
	end

	% BUTTON RESPONSE TEST */
	% read button from eyetracker */
	% This is the preferred way to get response data or end trials
	% fwc: I think we can do this equally well via KbCheck
  	button = eyelink('lastbuttonpress');
  	if button~=0           % button number, or 0 if none
   		EYELINK('message', 'ENDBUTTON %d', button); % log the button
    	break;               % exit loop
 	end

	if keyCode(spacebar) % stops trial
		rt=getsecs-tstart;
   		EYELINK('message', 'RT %f secs', rt); % log reaction time
		break; 
	end

	if getsecs-tstart>timeouttime
   		EYELINK('message', 'Timeout after %f secs', timeouttime); %message 
		break; 
	end

	% Real-time GAZE display code 

	if eyelink('isconnected') == el.dummyconnected
		sample=1;
	else
		sample=EYELINK( 'newfloatsampleavailable'); % eyelink_newest_float_sample(NULL) check for new sample update
	end
	if sample> 0
		if eyelink('isconnected') == el.dummyconnected
			evt=[];
		else
			evt = EYELINK( 'newestfloatsample'); % eyelink_newest_float_sample(&evt);  % get the sample */
		end
  		%evt
		if eye_used ~= -1   % do we know which eye yet? */
			if eyelink('isconnected') == el.dummyconnected
				[x,y,button] = GetMouse( el.window ); % get gaze position from mouse
				pupil=1;
			else
	        	x = evt.gx(eyetracked);       % get gaze position from sample */
	       		y = evt.gy(eyetracked);
				pupil=evt.pa(eyetracked);
			end
	
       		if x~=el.MISSING_DATA && y~=el.MISSING_DATA && pupil>0  % pupil visisible? */
 				if x~=oldx || y~=oldy % draw only if gaze moved
 					currGazeRect = CenterRectOnPoint(gazeRect, x, y);
 					unRect=UnionRect(currGazeRect,oldGazeRect);
 					SCREEN(screenbuffer, 'FillRect', el.backgroundcolour, oldGazeRect );
 					SCREEN('Copywindow',buffer,screenbuffer,currGazeRect,currGazeRect);
 					SCREEN('Copywindow',screenbuffer,el.window,unRect,unRect);
 					oldGazeRect=currGazeRect;
 					oldx=x;
 					oldy=y;
 				end % if x~=oldx || y~=oldy
       		else
           		% fprintf('Missing data\n');     % just hide if in blink */
				SCREEN( el.window, 'FillRect', el.backgroundcolour );  % clear display
				SCREEN( screenbuffer, 'FillRect', el.backgroundcolour );  % clear screenbuffer
				% reset params
				currGazeRect = [0 0 0 0 ];
				oldGazeRect=[0 0 0 0 ];
				oldx=-1;
				oldy=-1;
			end % if x~=MISSING_DATA
		else % eye_used== -1
	        eye_used = EYELINK( 'eyeavailable');  % find which eye
			% fprintf('eyeavailable %d\n', eye_used);
	        switch eye_used 			
	          % select eye, add annotation
	            case el.RIGHT_EYE,
	              EYELINK('message', 'EYE_USED 1 RIGHT');
	              %break;
	            case el.BINOCULAR,
	              eye_used = el.LEFT_EYE;
	              EYELINK('message', 'EYE_USED 0 LEFT');
				  %break;
	            case el.LEFT_EYE,
	              EYELINK('message', 'EYE_USED 0 LEFT');
	             % break;
			 end
			 eyetracked=eye_used+1; % we use eyetracked to access array data
		end % eye_used== -1
	else % sample< 0
		 %fprintf('No sample.\n') 
	end % if sample> 0

	
end % response loop

SCREEN( el.window, 'FillRect', el.backgroundcolour );  % hide display
Waitsecs(0.1); % record additional 100 msec of data
eyelink('stoprecording');

screen(buffer, 'close');
screen(screenbuffer, 'close');

% report result: 0=timeout, else button */
EYELINK('message', 'TRIAL_RESULT %d', button);
% HANDLE SPECIAL EXIT CONDITIONS */
result = EYELINK('checkrecordexit');
return;




