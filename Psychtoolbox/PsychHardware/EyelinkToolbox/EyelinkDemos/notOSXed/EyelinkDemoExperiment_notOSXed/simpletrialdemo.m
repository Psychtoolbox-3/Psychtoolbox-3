function result=simpletrialdemo(el, trial)

% demo for a simple gaze recording experiment
% created after SAMPLE EXPERIMENT FOR MACINTOSH, by Dave Stampe , 27 June '97 
% 15-05-01	fwc created first version
% 22-05-01	fwc	debugging with working eyelink.mex
% 02-06-01	fwc removed use of global el, as suggest by John Palmer.

spacebar=KbName('space');

SCREEN( el.window, 'FillRect', el.backgroundcolour );  % clear display

imgname='paarden.jpg';
timeouttime=10.0;
buffer=screen(el.window, 'openoffscreenwindow');
imageArray=IMREAD(imgname,'jpg');

SCREEN( buffer, 'FillRect', el.backgroundcolour );  % clear buffer
SCREEN(buffer,'PutImage',imageArray);

eye_used = -1;

% This supplies a title at the bottom of the eyetracker display
EYELINK('command', 'record_status_message ''TRIAL %d''', trial );

% Always send this message before starting to record.
% It marks the start of the trial and also
% contains trial condition data required for analysis.

EYELINK('message', 'TRIALID %d', trial );


% DO PRE-TRIAL DRIFT CORRECTION */
% We repeat if ESC key pressed to do setup. */
% Setup might also have erased any pre-drawn graphics. */
while 1	
	if EYELINK('isconnected')==el.notconnected   % Check link often so we don't lock up if tracker lost
		result=el.ABORT_EXPT;
		return;
	end;
   % DRIFT CORRECTION */
   % We let 'do_drift_correct' draw target, but it erases the screen too */
   % 3rd argument would be 0 to NOT draw a target */
	[x,y] = WindowCenter(el.window); % convenience routine part of eyelink toolbox
	% fprintf('drifcorr at % d %d\n', x, y );
   	error = dodriftcorrect(el, x, y, 1, 1);

  	if error==el.TERMINATE_KEY
		result=el.ABORT_EXPT;
		return;
	end;
    % repeat if ESC was pressed to access Setup menu
  	if(error~=el.ESC_KEY) break; end
end % while

% just for fun, draw a rectangle on eyelink screen
EYELINK('setofflinemode');	% Must be offline to draw to EyeLink screen
d=size(imageArray);
rect=[0 0 d(1) d(2)];
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

hidecursor;

% load image from buffer?
SCREEN(el.window, 'Waitblanking');
SCREEN('Copywindow',buffer,el.window);

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
	if keyCode(el.modifierkey) && keyCode(el.quitkey)   % alternative for cmd . which quits program directly
		fprintf( '\nUser requested break.\n');
		result=el.ABORT_EXPT;
		return;
	end

	if keyCode(el.quitkey) % stop trial
		SCREEN( el.window, 'FillRect', el.backgroundcolour );  % hide display
 		waitsecs(0.1);       % record additional 100 msec of data
  		EYELINK('stoprecording');
 		result=el.SKIP_TRIAL;   % return code to skip trial  		
 		return;
	end

	% BUTTON RESPONSE TEST */
	% read button from eyetracker */
	% This is the preferred way to get response data or end trials
  	button = eyelink('lastbuttonpress' );
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
end % response loop

SCREEN( el.window, 'FillRect', el.backgroundcolour );  % hide display
Waitsecs(0.1); % record additional 100 msec of data
eyelink('stoprecording');

screen(buffer, 'close');
showcursor;
% report result: 0=timeout, else button */
EYELINK('message', 'TRIAL_RESULT %d', button);
% HANDLE SPECIAL EXIT CONDITIONS */
result = EYELINK('checkrecordexit');
return;




