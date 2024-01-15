function rc = PsychEyelinkDispatchCallback(callArgs, msg)
% PsychEyelinkDispatchCallback implementes the EyeLink Core Graphics part
% of the EyeLink API. This "Core Graphics" part of our API is responsible
% for handling the times when the API and Host PC takes control of the eye
% tracking procedures. This includes the functionality to stream camera
% images during camera/participant setup, displaying targets at locations
% on the participant screen during calibration, validation, and drift/check
% and correction routines. Complimentary to handling the visual aspect of
% these operations contingent on display routines, the functionality
% implemented herewith also handles the playback of feedback sounds to
% the experimenter and participant for guiding these interactive
% procedures. During these modes of operation, this function also
% implements the forwarding of kepresses to the Host PC that are registered
% on the computer's keyboard which is running this implementation. The
% purpose of this is to make sure that bost Host and Display PCs are
% operating as identically in these modes of operation.
%
%
% This function is normally called from within the Eyelink() mex file.
% Normal user code only calls it once to supply the eyelink defaults struct.
% This is handled within the EyelinkInitDefaults.m file, so you generally
% should not have to worry about this. However, if you change settings in
% the el structure, you may need to call it yourself.
%
% To define which onscreen window the eye image should be
% drawn to, call it with the return value from EyelinkInitDefaults, e.g.,
% w=Screen('OpenWindow', ...);
% el=EyelinkInitDefaults(w);
% myEyelinkDispatchCallback(el);
%
%
% To actually receive and display the images, register this function as
% eyelink's callback:
%
%
% if Eyelink('Initialize', 'myEyelinkDispatchCallback') ~=0
%       error('eyelink failed init')
% end
% result = Eyelink('StartSetup',1) %put the tracker into a mode capable of sending images
%
%
% then you must hit 'return' on the PTB computer, this key command will be
% sent to the tracker host to initiate sending of images.
%
%
% History:
% 15. 3.2009    Derived from MemoryBuffer2TextureDemo.m (MK).
%  4. 4.2009    Updated to use EyelinkGetKey + fixed eyelinktex persistence
%                   crash (edf).
% 11. 4.2009    Cleaned up. Should be ready for 1st release, although still
%                   pretty alpha quality. (MK).
% 15. 6.2010    Added some drawing routines to get standard behaviour back.
%                   Enabled use of the callback by default. Clarified in
%                   helptext that user normally should not have to worry
%                   about calling this file. (fwc)
% 20. 7.2010    Drawing of instructions, eye-image+title, playing sounds in
%                   seperate functions
%
%  1. 2.2010    Modified to allow for cross hair and fix bugs. (nj)=
% 29.10.2018    Drop 'DrawDots' for calibration target. Some white-space fixes.
% 24. 3.2020    Cleaned up the documentation of this function, and added
%                   additiontal handling for two types of stereoscopic
%                   calibrations, ability to reference video files for
%                   animated calibration targets, bug fixes for audio
%                   feedback playback. Apologies to NJ for removing
%                   previous comments where code was previously added, this
%                   was done for easier reading of the code.

% Cached texture handle for eyelink texture:
persistent eyelinktex;
global dw dh offscreen;
global eyelinkanimationtarget;

% Cached window handle for target onscreen window:
persistent eyewin;
persistent calxy;
persistent imgtitle;
persistent eyewidth;
persistent eyeheight;

% Cached(!) eyelink stucture containing keycodes
persistent el;
persistent ineyeimagemodedisplay;
persistent clearScreen;
persistent drawInstructions;

% Cached constant definitions:
persistent GL_RGBA;
persistent GL_RGBA8;
persistent hostDataFormat;

persistent inDrift;
offscreen = 0;
newImage = 0;

if 0 == Screen('WindowKind', eyelinktex)
    eyelinktex = []; % Previous PTB Screen() window has closed, needs to be recreated.
end

if isempty(eyelinktex)
    % Define the two OpenGL constants we actually need. No point in
    % initializing the whole PTB OpenGL mode for just two constants:
    GL_RGBA = 6408;
    GL_RGBA8 = 32856;
    GL_UNSIGNED_BYTE = 5121; %#ok<NASGU>
    GL_UNSIGNED_INT_8_8_8_8 = 32821; %#ok<NASGU>
    GL_UNSIGNED_INT_8_8_8_8_REV = 33639;
    hostDataFormat = GL_UNSIGNED_INT_8_8_8_8_REV;
end

% Preinit return code to zero:
rc = 0;

if nargin < 2
    msg = [];
end

if nargin < 1
    callArgs = [];
end

if isempty(callArgs)
    error('You must provide some valid "callArgs" variable as 1st argument!');
end

if ~isnumeric(callArgs) && ~isstruct(callArgs)
    error('"callArgs" argument must be a EyelinkInitDefaults struct or double vector!');
end

% Eyelink el struct provided?
if isstruct(callArgs) && isfield(callArgs,'window')
    % Check if el.window subfield references a valid window:
    if Screen('WindowKind', callArgs.window) ~= 1
        error('argument didn''t contain a valid handle of an open onscreen window!  pass in result of EyelinkInitDefaults(previouslyOpenedPTBWindowPtr).');
    end
    
    % Ok, valid handle. Assign it and return:
    eyewin = callArgs.window;
    
    % Assume rest of el structure is valid:
    el = callArgs;
    clearScreen=1;
    eyelinktex=[];
    ineyeimagemodedisplay=0;
    drawInstructions=1;
    return;
end

% target & feedback beep waveforms and PsychPortAudio buffers
persistent audio_status;
persistent audio_devinfo;
persistent audio_n_chan;
persistent audio_fs;
%persistent audio_ppa_isSlave;
if ~isempty(el.ppa_pahandle) && isempty(audio_status)
    audio_status = PsychPortAudio('GetStatus', el.ppa_pahandle);
    audio_devinfo = PsychPortAudio('GetDevices', [], audio_status.OutDeviceIndex);
    audio_n_chan = min(2,audio_devinfo.NrOutputChannels);
    audio_fs = audio_status.SampleRate;
    %if PsychPortAudio('SetOpMode', pamaster) > 
    %    audio_ppa_isSlave
elseif isempty(el.ppa_pahandle) && isempty(audio_fs)
    audio_status = NaN;
    audio_devinfo = NaN;
    audio_n_chan = 1;
    audio_fs = Snd('DefaultRate');
end
persistent beep_waveforms;
if isempty(beep_waveforms)
    if el.targetbeep
        beep_waveforms{1} = repmat(MakeBeep(el.cal_target_beep(1), el.cal_target_beep(3), audio_fs) .* el.cal_target_beep(2), audio_n_chan, 1);
        beep_waveforms{2} = repmat(MakeBeep(el.drift_correction_target_beep(1), el.drift_correction_target_beep(3), audio_fs) .* el.drift_correction_target_beep(2), audio_n_chan, 1);
    else
        beep_waveforms{1} = NaN;
        beep_waveforms{2} = NaN;
    end
    if el.feedbackbeep
        beep_waveforms{3} = repmat(MakeBeep(el.calibration_failed_beep(1), el.calibration_failed_beep(3), audio_fs) .* el.calibration_failed_beep(2), audio_n_chan, 1);
        beep_waveforms{4} = repmat(MakeBeep(el.calibration_success_beep(1), el.calibration_success_beep(3),audio_fs) .* el.calibration_success_beep(2), audio_n_chan, 1);
        beep_waveforms{5} = repmat(MakeBeep(el.drift_correction_failed_beep(1), el.drift_correction_failed_beep(3), audio_fs) .* el.drift_correction_failed_beep(2), audio_n_chan, 1);
        beep_waveforms{6} = repmat(MakeBeep(el.drift_correction_success_beep(1), el.drift_correction_success_beep(3), audio_fs) .* el.drift_correction_success_beep(2), audio_n_chan, 1);
    else
        beep_waveforms{3} = NaN;
        beep_waveforms{4} = NaN;
        beep_waveforms{5} = NaN;
        beep_waveforms{6} = NaN;
    end
end
persistent ppa_beep_buffers;
if ~isempty(el.ppa_pahandle) && isempty(ppa_beep_buffers)
    if el.targetbeep
        fprintf('DEBUG: ppa_beep_buffers, targetbeep, create\n');
        ppa_beep_buffers(1) = PsychPortAudio('CreateBuffer', [], beep_waveforms{1});
        ppa_beep_buffers(2) = PsychPortAudio('CreateBuffer', [], beep_waveforms{2});
    else
        fprintf('DEBUG: ppa_beep_buffers, targetbeep, set NaN\n');
        ppa_beep_buffers(1) = NaN;
        ppa_beep_buffers(2) = NaN;
    end
    if el.feedbackbeep
        fprintf('DEBUG: ppa_beep_buffers, feedbackbeep, create\n');
        ppa_beep_buffers(3) = PsychPortAudio('CreateBuffer', [], beep_waveforms{3});
        ppa_beep_buffers(4) = PsychPortAudio('CreateBuffer', [], beep_waveforms{4});
        ppa_beep_buffers(5) = PsychPortAudio('CreateBuffer', [], beep_waveforms{5});
        ppa_beep_buffers(6) = PsychPortAudio('CreateBuffer', [], beep_waveforms{6});
    else
        fprintf('DEBUG: ppa_beep_buffers, feedbackbeep, set NaN\n');
        ppa_beep_buffers(3) = NaN;
        ppa_beep_buffers(4) = NaN;
        ppa_beep_buffers(5) = NaN;
        ppa_beep_buffers(6) = NaN;
    end
else
    fprintf('DEBUG: ppa_beep_buffers, skip\n');
end

% Not an eyelink struct.  Either a 4 component vector from Eyelink(), or something wrong:
if length(callArgs) ~= 4
    error('Invalid "callArgs" received from Eyelink() Not a 4 component double vector as expected!');
end

% Extract command code:
eyecmd = callArgs(1);

if isempty(eyewin) && eyecmd ~= 3
    warning('Got called as callback function from Eyelink() but usercode has not set a valid target onscreen window handle yet! Aborted.'); %#ok<WNTAG>
    return;
end

% (Re)set Flag for new camera image
newcamimage = 0;
needsupdate = 0;

switch eyecmd
    case 1  % New Camera Image Received
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 1; New Camera Image Received\n');
        end
        newcamimage = 1;
        needsupdate = 1;
        
    case 2  % EyeLink Keyboard Query
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 2; Keyboard Query\n');
        end
        [rc, el] = EyelinkGetKey(el);
        
    case 3  % Alert message
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 3; Alert Message\n');
        end
        fprintf('Eyelink Alert: %s.\n', msg);
        needsupdate = 0;
        
    case 4  % Camera Image Caption Text
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 4; Camera Image Caption Text\n');
        end
        if callArgs(2) ~= -1
            imgtitle = sprintf('Camera: %s [Threshold = %f]', msg, callArgs(2));
        else
            imgtitle = msg;
        end
        needsupdate = 1;
        
    case 5  % Draw Cal Target
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 5; Draw Cal Target\n');
        end
        calxy = callArgs(2:3);
        clearScreen=1;
        needsupdate = 1;
        if strcmpi(el.calTargetType, 'video') && ~isempty(eyelinkanimationtarget)
            if  el.calAnimationResetOnTargetMove && Screen('GetMovieTimeIndex', eyelinkanimationtarget.movie)
                Screen('SetMovieTimeIndex', eyelinkanimationtarget.movie, 0, el.calAnimationSetIndexIsFrames);
            end
            Screen('PlayMovie', eyelinkanimationtarget.movie, 1, el.calAnimationLoopParam, el.calAnimationAudioVolume);
        end
        
    case 6  % Clear Cal Display
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 6; Clear Cal Display\n');
        end
        clearScreen=1;
        drawInstructions=1;
        needsupdate = 1;
        
    case 7  % Setup Cal Display
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 7; Setup Cal Display\n');
        end
        if inDrift
            drawInstructions = 0;
        else
            drawInstructions = 1;
        end
        clearScreen=1;
        needsupdate = 1;
        
    case 8  % Setup Image Display
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 8; Setup Image Display\n');
        end
        newImage = 1;
        eyewidth  = callArgs(2);
        eyeheight = callArgs(3);
        ineyeimagemodedisplay=1;
        drawInstructions=1;
        needsupdate = 1;
        
    case 9  % Exit Image Display
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 9; Exit Image Display\n');
        end
        clearScreen=1;
        ineyeimagemodedisplay=0;
        drawInstructions=1;
        needsupdate = 1;
        
    case 10 % Erase Cal Target
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 10; Erase Cal Target\n');
        end
        calxy = [];
        if ~isempty(eyelinkanimationtarget)
            eyelinkanimationtarget.calxy=calxy;
        end
        clearScreen=1;
        needsupdate = 1;
        
        if strcmpi(el.calTargetType, 'video') && ~isempty(eyelinkanimationtarget)
             texkill=Screen('GetMovieImage', eyewin, eyelinkanimationtarget.movie, el.calAnimationWaitTexClose);
             Screen('PlayMovie', eyelinkanimationtarget.movie, 0, el.calAnimationLoopParam);
             if texkill > 0
                 Screen('Close', texkill);
             end
        end
        
    case 11 % Exit Cal Display
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 11; Exit Cal Display\n');
        end
        calxy = [];
        if(~isempty(eyelinkanimationtarget)	)
            eyelinkanimationtarget.calxy=calxy;
        end
        
        if inDrift
            inDrift = 0;
            drawInstructions = 0;
        else
            drawInstructions = 1;
        end
        
        clearScreen=1;
        needsupdate=1;
        if strcmpi(el.calTargetType, 'video') && ~isempty(eyelinkanimationtarget)
            texkill=Screen('GetMovieImage', eyewin, eyelinkanimationtarget.movie, el.calAnimationWaitTexClose);
            Screen('PlayMovie', eyelinkanimationtarget.movie, 0, el.calAnimationLoopParam);
            if texkill > 0
                Screen('Close', texkill);
            end
        end

    case 12 % New Cal Target Sound
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 12; New Cal Target Sound\n');
        end
        if el.targetbeep && ~strcmpi(el.calTargetType, 'video')
            EyelinkMakeSound(el, 'cal_target_beep', 1);
        end
        
    case 13 % New Drift Chk/Corr Target Sound
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 13; New Drift Target Sound\n');
        end
        if el.targetbeep && ~strcmpi(el.calTargetType, 'video')
            EyelinkMakeSound(el, 'drift_correction_target_beep', 2);
        end
        
    case 14 % Cal Done Sound
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 14; Cal Done Sound\n');
        end
        if el.feedbackbeep && ~strcmpi(el.calTargetType, 'video')
            errc = callArgs(2);
            if errc > 0
                % Failed
                EyelinkMakeSound(el, 'calibration_failed_beep', 3);
            else
                % Success
                EyelinkMakeSound(el, 'calibration_success_beep', 4);
            end
        end
        
    case 15 % Drift Chk/Corr Done Sound
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 15; Drift Done Sound\n');
        end
        if el.feedbackbeep && ~strcmpi(el.calTargetType, 'video')
            errc = callArgs(2);
            if errc > 0
                % Failed
                EyelinkMakeSound(el, 'drift_correction_failed_beep', 5);
            else
                % Success
                EyelinkMakeSound(el, 'drift_correction_success_beep', 6);
            end
        end
        
    case 16 % Get Mouse Position
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 16; Get Mouse Position\n');
        end
        [width, height]=Screen('WindowSize', eyewin);
        [x,y, buttons] = GetMouse(eyewin);
        HideCursor(eyewin);
        if find(buttons)
            rc = [width , height, x , y,  dw , dh , 1];
        else
            rc = [width , height, x , y , dw , dh , 0];
        end
        
    case 17 % Non-native callback, from PsychEyelink_setup_cal_display()
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == 17; Flag in drift check/correction mode\n');
        end
        inDrift = 1;
    
    case -1 % Non-native callback, from Eyelink('Shutdown') for runtime cleanup
        if Eyelink('Verbosity', -1) >= 5
            fprintf('PsychEyelinkDispatchCallback: eyecmd == -1; Runtime cleanup\n');
        end
        if any(strcmp(fieldnames(el), 'ppa_sndhandle')) && ~isempty(el.ppa_sndhandle)
            PsychPortAudio('Close', el.ppa_sndhandle);
            el.ppa_sndhandle = [];
        end
        % vars = whos;
        % whos
        % clear(vars([vars.persistent]).name{:});
        % whos
        return;

    otherwise % Unknown Command
        fprintf('PsychEyelinkDispatchCallback: Unknown eyelink command (%i)\n', eyecmd);
        return;
        
end

if ~needsupdate     % Display redraw and update needed?
    if ~isempty(eyelinkanimationtarget) && ~isempty(calxy)
        EyelinkDrawCalibrationTarget(eyewin, el, calxy);
    end
    return;
end

if clearScreen==1   % Need to clear and redraw display before new content flipped
    EyelinkDrawClearScreen(eyewin, el); 
    clearScreen=0;
end

if newcamimage      % New image frame received from EyeLink camera stream
    % Image has dimensions: 'eyewidth' by 'eyeheight' in pixel units
    % Each pixel is encoded as a 4 byte RGBA pixel with A=255 always
    % RGB channels each encode a 1-Byte per channel R, G or B color value.
    % 'eyeimgptr' is a memory pointer to the buffer inside Eyelink() that
    % encodes the image.
    eyeimgptr = callArgs(2);
    eyewidth  = callArgs(3);
    eyeheight = callArgs(4);
    
    % Creates a new or reuses an existing PTB texture for the cam image
    eyelinktex = Screen('SetOpenGLTextureFromMemPointer', eyewin, eyelinktex, eyeimgptr, eyewidth, eyeheight, 4, 0, [], GL_RGBA8, GL_RGBA, hostDataFormat);
end

if ~isempty(eyelinktex) && ineyeimagemodedisplay==1     % Draw cam image and caption
    imgtitle=EyelinkDrawCameraImage(eyewin, el, eyelinktex, imgtitle, newImage);
end

if ~isempty(calxy)  % Draw Cal Target
    drawInstructions=0;
    EyelinkDrawCalibrationTarget(eyewin, el, calxy);
    if strcmpi(el.calTargetType, 'video')
        return;
    end
end

if drawInstructions == 1  % Draw Instructions
    EyelinkDrawInstructions(eyewin, el, msg);
    drawInstructions = 0;
end

dontsync = 1;
dontclear = 1; % set to 0 to hide instructions during camera display
Screen('Flip', eyewin, [], dontclear, dontsync);   % Show it

return;



% Start of nested EyelinkDraw* function declarations

    function EyelinkDrawClearScreen(eyewin, el)
        if el.winInfo.StereoMode ~= 0 
            drawScreens = 2; % stereoscopic drawing
        else
            drawScreens = 1; % non-stereoscopic drawing
        end
        for it = 0:drawScreens-1 
            Screen('SelectStereoDrawBuffer', eyewin, it); % select left eye window
            Screen('FillRect', eyewin, el.backgroundcolour);
        end
    end

    function EyelinkDrawInstructions(eyewin, el,msg)
        if el.winInfo.StereoMode ~= 0 
            drawScreens = 2; % stereoscopic drawing
        else
            drawScreens = 1; % non-stereoscopic drawing
        end
        oldFont=Screen(eyewin,'TextFont',el.msgfont);
        oldFontSize=Screen(eyewin,'TextSize',el.msgfontsize);
        for it = 0:drawScreens-1 
            Screen('SelectStereoDrawBuffer', eyewin, it); % select left eye window
            DrawFormattedText(eyewin, el.helptext, 20, 20, el.msgfontcolour, [], [], [], 1);
            
            if el.displayCalResults && ~isempty(msg)
                DrawFormattedText(eyewin, msg, 20, 150, el.msgfontcolour, [], [], [], 1);
            end
        end
        Screen(eyewin,'TextFont',oldFont);
        Screen(eyewin,'TextSize',oldFontSize);
    end



    function  imgtitle=EyelinkDrawCameraImage(eyewin, el, eyelinktex, imgtitle, newImage)
        persistent lasttitle;
        %         global dh dw offscreen;
        if el.winInfo.StereoMode ~= 0 
            drawScreens = 2; % stereoscopic drawing
        else
            drawScreens = 1; % non-stereoscopic drawing
        end
        
        for it = 0:drawScreens-1 
            try
                if ~isempty(eyelinktex)
                    Screen('SelectStereoDrawBuffer', eyewin, it); % select left eye window
                    eyerect=Screen('Rect', eyelinktex);
                    % we could cash some of the below values....
                    wrect=Screen('Rect', eyewin);
                    [width, heigth]=Screen('WindowSize', eyewin);
                    dw=round(el.eyeimgsize/100*width);
                    dh=round(dw * eyerect(4)/eyerect(3));
                    
                    drect=[ 0 0 dw dh ];
                    drect=CenterRect(drect, wrect);
                    Screen('DrawTexture', eyewin, eyelinktex, [], drect);
                end
                % imgtitle
                % if title is provided, we also draw title
                if ~isempty(eyelinktex) && exist( 'imgtitle', 'var') && ~isempty(imgtitle)
                    Screen('SelectStereoDrawBuffer', eyewin, it); % select left eye window
                    rect=Screen('TextBounds', eyewin, imgtitle );
                    [w2, h2]=RectSize(rect);
                    
                    if -1 == Screen('WindowKind', offscreen)
                        Screen('Close', offscreen);
                    end
                    Screen('SelectStereoDrawBuffer', eyewin, it); % select left eye window
                    sn = Screen('WindowScreenNumber', eyewin);
                    offscreen = Screen('OpenOffscreenWindow', sn, el.backgroundcolour, [], [], 32);
                    Screen(offscreen,'TextFont',el.imgtitlefont);
                    Screen(offscreen,'TextSize',el.imgtitlefontsize);
                    Screen('DrawText', offscreen, imgtitle, width/2-dw/2, heigth/2+dh/2+h2, el.imgtitlecolour);
                    Screen('SelectStereoDrawBuffer', eyewin, it); % select left eye window
                    Screen('DrawTexture',eyewin,offscreen,  [width/2-dw/2 heigth/2+dh/2+h2 width/2-dw/2+500 heigth/2+dh/2+h2+500], [width/2-dw/2 heigth/2+dh/2+h2 width/2-dw/2+500 heigth/2+dh/2+h2+500]);
                    Screen('Close',offscreen);
                    
                    lasttitle = imgtitle;
                end
            catch %myerr
                fprintf('EyelinkDrawCameraImage:error \n');
                disp(psychlasterror);
            end
        end
    end





    function EyelinkDrawCalibrationTarget(eyewin, el, calxy)
        if el.winInfo.StereoMode ~= 0
            drawScreens = 2; % stereoscopic drawing
        else
            drawScreens = 1; % non-stereoscopic drawing
        end

        [width, heigth]=Screen('WindowSize', eyewin);

        for it = 0:drawScreens-1
            Screen('SelectStereoDrawBuffer', eyewin, it); % select eye window
            switch el.calTargetType
                case 'video'
                    if( ~isempty(el.calAnimationTargetFilename) && ~isempty(eyelinkanimationtarget))
                        rect=CenterRectOnPoint([0 0 eyelinkanimationtarget.imgw eyelinkanimationtarget.imgh], calxy(1), calxy(2));
                        if it == 0
                            tex=Screen('GetMovieImage', eyewin, eyelinkanimationtarget.movie, el.calAnimationWaitTex);
                        end
                        if(tex>0)
                            Screen('DrawTexture', eyewin , tex, [], rect, [], 0);
                            if drawScreens == 1 || (drawScreens == 2 && it == 1)
                                Screen('Flip', eyewin);
                            end
                        end
                        if it == drawScreens-1 && tex > 0
                            Screen('Close', tex);
                        end
                    end
                    
                case 'image'
                    if ~isempty(el.calImageInfo) && ~isempty(el.calImageTexture)
                        rect=CenterRectOnPoint([0 0 el.calImageInfo.Width el.calImageInfo.Height], calxy(1), calxy(2));
                    end
                    Screen('DrawTexture', eyewin , el.calImageTexture, [], rect, [], 0);
                    
                otherwise
                    % default to el.calTargetType = 'ellipse' target drawing
                    size=round(el.calibrationtargetsize / 100 * width);
                    inset=round(el.calibrationtargetwidth / 100 * width);
                    
                    % Use FillOval for larger dots:
                    Screen('FillOval', eyewin, el.calibrationtargetcolour, [calxy(1)-size/2 calxy(2)-size/2 calxy(1)+size/2 calxy(2)+size/2], size+2);
                    Screen('FillOval', eyewin, el.backgroundcolour, [calxy(1)-inset/2 calxy(2)-inset/2 calxy(1)+inset/2 calxy(2)+inset/2], inset+2)
            end
        end
    end


    function EyelinkMakeSound(el, s, i)
        % set all sounds in one place, sound params defined in
        % eyelinkInitDefaults
        if any(strcmp( ...
            {'cal_target_beep', ...
            'drift_correction_target_beep', ...
            'calibration_failed_beep', ...
            'calibration_success_beep', ...
            'drift_correction_failed_beep', ...
            'drift_correction_success_beep'},s))
            % beep waveform was prepared in advance
            if isempty(el.ppa_pahandle)
                Snd('Play', beep_waveforms{i}, audio_fs);
            else
                PsychPortAudio('FillBuffer', el.ppa_pahandle, ppa_beep_buffers(i));
                PsychPortAudio('Start', el.ppa_pahandle);
                if strcmp('drift_correction_success_beep', s)
                    PsychPortAudio('Stop', el.ppa_pahandle, 1);
                end
            end

        else
            % some defaults
            if isempty(el.ppa_pahandle)
                Snd('Play', MakeBeep(500, 1.5) .* 0.5, audio_fs);
            else
                PsychPortAudio('CreateBuffer', el.ppa_pahandle, repmat(MakeBeep(500, 1.5, audio_fs) .* 0.5, audio_n_chan, 1));
                PsychPortAudio('Start', el.ppa_pahandle);
            end
        end
    end
end



