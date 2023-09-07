function result=EyelinkLegacyTargetModeDisplay(el)

% USAGE: result=EyelinkLegacyTargetModeDisplay(el)
%
%		el: Eyelink default values
% History
% 15-05-01	fwc created first version
% 22-05-01	fwc	little debugging
% 02-06-01	fwc removed use of global el, as suggested by John Palmer.
%   22-06-06    fwc OSX-ed


result=-1; % initialize
if nargin < 1
	error( 'USAGE: result=EyelinkTargetModeDisplay(el)' );
end

targetvisible = 0;	% target currently drawn
targetrect=[0 0 0 0];

tx=el.MISSING;
ty=el.MISSING;

otx=el.MISSING;    % current target position
oty=el.MISSING;

EyelinkLegacyClearCalDisplay(el);	% setup_cal_display()

key=1;
while key~= 0
	[key, el]=EyelinkGetKey(el);		% dump old keys
end
				% LOOP WHILE WE ARE DISPLAYING TARGETS
stop=0;
while stop==0 && bitand(Eyelink('CurrentMode'), el.IN_TARGET_MODE)

	if Eyelink( 'IsConnected' )==el.notconnected
		result=-1;
		return;
	end;

	[key, el]=EyelinkGetKey(el);		% getkey() HANDLE LOCAL KEY PRESS

	switch key 
		case el.TERMINATE_KEY,       % breakout key code
			EyelinkLegacyClearCalDisplay(el); % clear_cal_display();
			result=el.TERMINATE_KEY;
			return;
		case el.SPACE_BAR,	         		% 32: accept fixation
			if el.allowlocaltrigger==1
				Eyelink( 'AcceptTrigger');
			end
			break;
		case { 0,  el.JUNK_KEY	}	% No key
		case el.ESC_KEY,
			if Eyelink('IsConnected') == el.dummyconnected
				stop=1;
			end
		    if el.allowlocalcontrol==1 
	       		Eyelink('SendKeyButton', key, 0, el.KB_PRESS );
			end
		otherwise,          % Echo to tracker for remote control
		    if el.allowlocalcontrol==1 
	       		Eyelink('SendKeyButton', key, 0, el.KB_PRESS );
			end
	end % switch key


				% HANDLE TARGET CHANGES
	[result, tx, ty]= Eyelink( 'TargetCheck');
	
	
	% erased or moved: erase target
	if (targetvisible==1 && result==0) || tx~=otx || ty~=oty
		EyelinkLegacyEraseCalTarget(el, targetrect);
		targetvisible = 0;
	end
	% redraw if invisible
	if targetvisible==0 && result==1
% 		fprintf( 'Target drawn at: x=%d, y=%d\n', tx, ty );
		
		targetrect=EyelinkLegacyDrawCalTarget(el, tx, ty);
		targetvisible = 1;
		otx = tx;		% record position for future tests
		oty = ty;
		if el.targetbeep==1
			EyelinkLegacyCalTargetBeep(el);	% optional beep to alert subject
		end
	end
	
end % while IN_TARGET_MODE


% exit:					% CLEAN UP ON EXIT
if el.targetbeep==1
	if Eyelink('CalResult')==1  % does 1 signal success?
		EyelinkLegacyCalDoneBeep(el, 1);
	else
	  	EyelinkLegacyCalDoneBeep(el, -1);
	end
end
  
if targetvisible==1
	EyelinkLegacyEraseCalTarget(el, targetrect);   % erase target on exit, bit superfluous actually
end
EyelinkLegacyClearCalDisplay(el); % clear_cal_display();

result=0;
return;
