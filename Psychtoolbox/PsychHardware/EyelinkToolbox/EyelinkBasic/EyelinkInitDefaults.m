function el=EyelinkInitDefaults(window)
% Initialize eyelink defaults and control code structure
%
% USAGE: el=EyelinkInitDefaults([window])
%
%       window is optional windowPtr.
%       If set, pixel coordinates are send to eyetracker
%
% and fill it with some sensible values.
% Note that these values are only used by the m-file
% versions of dotrackersetup and dodriftcorrect.

% 02-06-01	fwc created, as suggested by John Palmer.
%				added also all control codes and defaults
% 17-10-02	fwc added event types
% 26-11-02  fwc&emp added PC support
% 11-01-04  fwc OS X changes
% 22-06-06  fwc further OSX changes
% 31-10-06  mk  Unified keyname mapping and such...
% 19-02-09  edf added LOSTDATAEVENT
% 27-03-09  edf added function and modifier keys.
% 10-04-09  mk  Deuglified. Add setup code for
%               PsychEyelinkDispatchCallback.
% 15-06-10  fwc For consistency, changed to EyelinkDispatchCallback
%               Added additional default values for message font size
%               and eye image size. Note that many default settings
%               are no longer used in the "callback" version of calibration
%               and driftcorrection.
% 15-01-13	ia	Added el.devicenumber to allow better control of multiple
%				input devices

el=[];

% eyelink computer check
el.computer = computer;

% Enable unified keyname -> keycode mapping for all operating systems:
KbName('UnifyKeyNames');

% eye codes
el.LEFT_EYE=0;
el.RIGHT_EYE=1;
el.BINOCULAR=2;

% eyelink connection states
el.notconnected=0;
el.connected=1;
el.dummyconnected=-1;
el.broadcastconnected=2;


el.displayCalResults = 0;

if ~exist('window', 'var')
    window = [];
end

if ~isempty(window)
	el.window=window;
	el.backgroundcolour = WhiteIndex(el.window);
	el.backgroundcolour = GrayIndex(el.window);
	el.foregroundcolour = BlackIndex(el.window);
    el.msgfontcolour    = BlackIndex(el.window);
    el.imgtitlecolour   = BlackIndex(el.window);

	rect=Screen(el.window,'Rect');
    if Eyelink('IsConnected') ~= el.notconnected
        Eyelink('Command', 'screen_pixel_coords = %d %d %d %d',rect(1),rect(2),rect(3)-1,rect(4)-1);
    end
else
	el.window=[];
end

% set some more global info parameters
% below are old sound definitions
el.targetdisplaysound='EyelinkTargetBeep';
el.calibrationfailedsound='EyelinkErrorBeep';
el.calibrationsuccesssound='EyelinkSuccessBeep';
el.targetbeep=1;  % sound a beep when a target is presented
el.feedbackbeep=1;  % sound a beep after calibration/drift correction

% define beep sounds (frequency, volume, duration);
el.cal_target_beep=[1250 0.6 0.05];
el.drift_correction_target_beep=[1250 0.8 0.05];
el.calibration_failed_beep=[400 0.8 0.25];
el.calibration_success_beep=[800 0.8 0.25];
el.drift_correction_failed_beep=[400 0.8 0.25];
el.drift_correction_success_beep=[800 0.8 0.25];

el.allowlocaltrigger=1; % allow user to trigger him or herself
el.allowlocalcontrol=1; % allow control from subject-computer
el.mousetriggersdriftcorr=0; % 1=allow mouse to trigger drift correction (fwc trick)
el.quitkey=KbName('ESCAPE'); % when pressed in combination with modifier key
                             % forces getkeyforeyelink to return 'TERMINATE_KEY' !

% Modifier key is always LeftGUI due to unified keyname mapping:
el.modifierkey=KbName('LeftGUI');

el.waitformodereadytime=500;
el.calibrationtargetsize=2.5;  % size of calibration target as percentage of screen
el.calibrationtargetwidth=1; % width of calibration target's border as percentage of screen
el.calibrationtargetcolour=[0 0 0];

el.devicenumber = []; %see KbCheck for details of this value
el.getkeyrepeat=1/5; % "sample time" for eyelinkgetkey
el.getkeytime=-1; % stores last time eyelinkgetkey was used

% (font info for ) messages/instructions
el.msgfont='Helvetica';
el.msgfontsize=20; % absolute, should perhaps be percentage of screen
el.eyeimgsize=30; % percentage of screen
el.helptext='Press RETURN (on either display computer or tracker host computer) to toggle camera image';
el.helptext=[el.helptext '\n' 'Press ESC to Output/Record'];
el.helptext=[el.helptext '\n' 'Press C to Calibrate'];
el.helptext=[el.helptext '\n' 'Press V to Validate'];
% el.helptext=[el.helptext '\n' 'Press D for Drift correction'];

% font info for camera image title
el.imgtitlefont='Helvetica';
el.imgtitlefontsize=20; % should be percentage of screen


% Warm up KbCheck:
[keyIsDown, secs, el.lastKeyCodes] = KbCheck;

% keyCodes for EyelinkGetKey:
el.uparrow=KbName('UpArrow');
el.downarrow=KbName('DownArrow');
el.rightarrow=KbName('RightArrow');
el.leftarrow=KbName('LeftArrow');
el.pageup=KbName('PageUp');
el.pagedown=KbName('PageDown');
el.return=KbName('Return');
el.escape=KbName('ESCAPE');
el.space=KbName('space');
el.backspace=KbName('DELETE'); % is this delete backspace?

el.f1=KbName('F1');
el.f2=KbName('F2');
el.f3=KbName('F3');
el.f4=KbName('F4');
el.f5=KbName('F5');
el.f6=KbName('F6');
el.f7=KbName('F7');
el.f8=KbName('F8');
el.f9=KbName('F9');
el.f10=KbName('F10');

el.left_shift=KbName('LeftShift');
el.right_shift=KbName('RightShift');
el.left_control=KbName('LeftControl');
el.right_control=KbName('RightControl');
el.lalt=KbName('LeftAlt');
el.ralt=KbName('RightAlt');
%el.lmeta=KbName('');
%el.rmeta=KbName('');
el.num=KbName('NumLock');
el.caps=KbName('CapsLock');
%el.mode=KbName('');

if IsOSX
	% OS-X supports a separate keycode for the Enter key:
	el.enter=KbName('ENTER');
else
	% M$-Windows and GNU/Linux don't have a separate code for Enter,
	% so we will map it to the 'Return' key:
	el.enter=el.return;
end
el.keysCached=1;

% Since we do not actually remove keypresses the Matlab buffer gets filled
% up quickly. Hence we disable warnings for fillup problems:
% This is try-catch protected for compatibility to Matlab R11 and Octave...
try
	warning off MATLAB:namelengthmaxexceeded
catch
    % Nothing to do. We just swallow the error we'd get if that warning
    % statement wouldn't be supported.
end

% Eyelink Tracker state bit: bitand() with flag word to test functionality

el.IN_DISCONNECT_MODE=16384;   	% disconnected
el.IN_UNKNOWN_MODE=0;    		% mode fits no class (i.e setup menu)
el.IN_IDLE_MODE=1;    			% off-line
el.IN_SETUP_MODE=2;   			% setup or cal/val/dcorr
el.IN_RECORD_MODE=4;    		% data flowing
el.IN_TARGET_MODE=8;    		% some mode that needs fixation targets
el.IN_DRIFTCORR_MODE=16;   		% drift correction
el.IN_IMAGE_MODE=32;   			% image-display mode
el.IN_USER_MENU=64;				% user menu
el.IN_PLAYBACK_MODE=256;		% tracker sending playback data

% Eyelink key values
el.JUNK_KEY=1;		% return code for untranslatable key
el.TERMINATE_KEY=hex2dec('7FFF'); % return code for program exit/breakout key
el.CURS_UP=hex2dec('4800');
el.CURS_DOWN=hex2dec('5000');
el.CURS_LEFT=hex2dec('4B00');
el.CURS_RIGHT=hex2dec('4D00');
el.ESC_KEY=hex2dec('001B');
el.ENTER_KEY=hex2dec('000D');
el.PAGE_UP=hex2dec('4900');
el.PAGE_DOWN=hex2dec('5100');
el.SPACE_BAR=32;

el.F1_KEY=hex2dec('3B00');
el.F2_KEY=hex2dec('3C00');
el.F3_KEY=hex2dec('3D00');
el.F4_KEY=hex2dec('3E00');
el.F5_KEY=hex2dec('3F00');
el.F6_KEY=hex2dec('4000');
el.F7_KEY=hex2dec('4100');
el.F8_KEY=hex2dec('4200');
el.F9_KEY=hex2dec('4300');
el.F10_KEY=hex2dec('4400');

% bitand these into the key if they are also pressed
el.ELKMOD_NONE=hex2dec('0000');
el.ELKMOD_LSHIFT=hex2dec('0001');
el.ELKMOD_RSHIFT=hex2dec('0002');
el.ELKMOD_LCTRL=hex2dec('0040');
el.ELKMOD_RCTRL=hex2dec('0080');
el.ELKMOD_LALT=hex2dec('0100');
el.ELKMOD_RALT=hex2dec('0200');
el.ELKMOD_LMETA=hex2dec('0400');
el.ELKMOD_RMETA=hex2dec('0800');
el.ELKMOD_NUM=hex2dec('1000');
el.ELKMOD_CAPS=hex2dec('2000');
el.ELKMOD_MODE=hex2dec('4000');

% other Eyelink values


el.ELKEY_DOWN=1;
el.ELKEY_UP=0;

el.KB_PRESS=10; % Eyelink.h
el.MISSING=-32768; % eyedata.h
el.MISSING_DATA=-32768;

el.KEYDOWN=1; %Eyelink manual and core_expt.h have these backwards
el.KEYUP=0;

% LINK RETURN CODES
el.NO_REPLY=1000; % no reply yet (for polling test)

% return codes for trial result
el.DONE_TRIAL=0;
el.TRIAL_OK=0;
el.REPEAT_TRIAL=1;
el.SKIP_TRIAL=2;
el.ABORT_EXPT=3;

el.TRIAL_ERROR=-1; 	% Bad trial: no data, etc.


% EVENT types
el.SAMPLE_TYPE=200;

el.STARTPARSE=1; % 	/* these only have time and eye data */
el.ENDPARSE=2;
el.BREAKPARSE=10;

el.STARTBLINK=3;    % /* EYE DATA: contents determined by evt_data */
el.ENDBLINK=4;   	% /* and by "read" data item */
el.STARTSACC=5;		% /* all use IEVENT format */
el.ENDSACC=6;
el.STARTFIX=7;
el.ENDFIX=8;
el.FIXUPDATE=9;

el.STARTSAMPLES=15;  	%/* start of events in block *//* control events: all put data into */
el.ENDSAMPLES=16;  		%/* end of samples in block *//* the EDF_FILE or ILINKDATA status  */
el.STARTEVENTS=17; 		% /* start of events in block */
el.ENDEVENTS=18;  		%/* end of events in block */

el.MESSAGEEVENT=24;  % /* user-definable text or data */
el.BUTTONEVENT=25;  %/* button state change */
el.INPUTEVENT=28;  % /* change of input port */

el.LOSTDATAEVENT=hex2dec('3F'); %/*new addition v2.1, returned by eyelink_get_next_data() to flag a gap in the data stream due to queue filling up (need to get data more frequently)
                                %/*described in 'EyeLink Programmers Guide.pdf' section 7.2.2, 13.3.2, 18.5.4

% if exist('EyelinkDispatchCallback') %#ok<EXIST>
%     el.callback = 'EyelinkDispatchCallback';
% else
%     el.callback = [];
% end

if exist('PsychEyelinkDispatchCallback') %#ok<EXIST>
    el.callback = 'PsychEyelinkDispatchCallback';
else
    el.callback = [];
end


EyelinkUpdateDefaults(el);


% % Window assigned?
% if ~isempty(el.window) && ~isempty(el.callback)
%     % Yes. Assign it to our dispatch callback:
% %     EyelinkDispatchCallback(el);
%     PsychEyelinkDispatchCallback(el);
% end

% el % uncomment to show contents of this default eyelink structure
