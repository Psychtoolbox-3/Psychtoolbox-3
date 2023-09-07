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
% 15.3.2020     br added Snd('Close') after Beeper to free sound device
%               and prevent problems downstream with PsychPortAudio; changed
%               flip 'dontsync' to '0' to fix missing target flips on Linux;
%               cleaned command 7 & 11 to fix issue drawing instructions
%

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


% Not an eyelink struct.  Either a 4 component vector from Eyelink(), or something wrong:
if length(callArgs) ~= 4
    error('Invalid "callArgs" received from Eyelink() Not a 4 component double vector as expected!');
end

% Extract command code:
eyecmd = callArgs(1);
% fprintf('%.3f - eyecmd: %d\n', GetSecs, eyecmd); % for debug

if isempty(eyewin) && eyecmd ~= 3
    warning('Got called as callback function from Eyelink() but usercode has not set a valid target onscreen window handle yet! Aborted.'); %#ok<WNTAG>
    return;
end

% (Re)set Flag for new camera image
newcamimage = 0;
needsupdate = 0;

switch eyecmd
    case 1  % New Camera Image Received
        newcamimage = 1;
        needsupdate = 1;
        
    case 2  % EyeLink Keyboard Query
        [rc, el] = EyelinkGetKey(el);
        
    case 3  % Alert message
        fprintf('Eyelink Alert: %s.\n', msg);
        needsupdate = 0;
        
    case 4  % Camera Image Caption Text
        if callArgs(2) ~= -1
            imgtitle = sprintf('Camera: %s [Threshold = %f]', msg, callArgs(2));
        else
            imgtitle = msg;
        end
        needsupdate = 1;
        
    case 5  % Draw Cal Target
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
        clearScreen=1;
        drawInstructions=1;
        needsupdate = 1;
        
    case 7  % Setup Cal Display
        if inDrift
            drawInstructions = 0;
        else
            drawInstructions = 1;
        end
        clearScreen=1;
        needsupdate = 1;
        
    case 8  % Setup Image Display
        newImage = 1;
        eyewidth  = callArgs(2);
        eyeheight = callArgs(3);
        ineyeimagemodedisplay=1;
        drawInstructions=1;
        needsupdate = 1;
        
    case 9  % Exit Image Display
        clearScreen=1;
        ineyeimagemodedisplay=0;
        drawInstructions=1;
        needsupdate = 1;
        
    case 10 % Erase Cal Target
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
        if ~strcmpi(el.calTargetType, 'video')
            EyelinkMakeSound(el, 'cal_target_beep');
        end
        
    case 13 % New Drift Chk/Corr Sound
        if ~strcmpi(el.calTargetType, 'video')
            EyelinkMakeSound(el, 'drift_correction_target_beep');
        end
        
    case 14 % Cal Done Sound
        if ~strcmpi(el.calTargetType, 'video')
            errc = callArgs(2);
            if errc > 0
                % Failed
                EyelinkMakeSound(el, 'calibration_failed_beep');
            else
                % Success
                EyelinkMakeSound(el, 'calibration_success_beep');
            end
        end
        
    case 15 % Drift Chk/Corr Done Sound
        if ~strcmpi(el.calTargetType, 'video')
            errc = callArgs(2);
            if errc > 0
                % Failed
                EyelinkMakeSound(el, 'drift_correction_failed_beep');
            else
                % Success
                EyelinkMakeSound(el, 'drift_correction_success_beep');
            end
        end
        
    case 16 % Get Mouse Position
        [width, height]=Screen('WindowSize', eyewin);
        [x,y, buttons] = GetMouse(eyewin);
        HideCursor;
        if find(buttons)
            rc = [width , height, x , y,  dw , dh , 1];
        else
            rc = [width , height, x , y , dw , dh , 0];
        end
        
    case 17 %
        inDrift =1;
        
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
dontclear = 0;
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
        
        for it = 0:drawScreens-1 
            Screen('SelectStereoDrawBuffer', eyewin, it); % select left eye window
            oldFont=Screen(eyewin,'TextFont',el.msgfont);
            oldFontSize=Screen(eyewin,'TextSize',el.msgfontsize);
            DrawFormattedText(eyewin, el.helptext, 20, 20, el.msgfontcolour, [], [], [], 1);
            
            if el.displayCalResults && ~isempty(msg)
                DrawFormattedText(eyewin, msg, 20, 150, el.msgfontcolour, [], [], [], 1);
            end
            
            Screen(eyewin,'TextFont',oldFont);
            Screen(eyewin,'TextSize',oldFontSize);
        end
        
        
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


    function EyelinkMakeSound(el, s)
        % set all sounds in one place, sound params defined in
        % eyelinkInitDefaults
        switch(s)
            case 'cal_target_beep'
                doBeep=el.targetbeep;
                f=el.cal_target_beep(1);
                v=el.cal_target_beep(2);
                d=el.cal_target_beep(3);
            case 'drift_correction_target_beep'
                doBeep=el.targetbeep;
                f=el.drift_correction_target_beep(1);
                v=el.drift_correction_target_beep(2);
                d=el.drift_correction_target_beep(3);
            case 'calibration_failed_beep'
                doBeep=el.feedbackbeep;
                f=el.calibration_failed_beep(1);
                v=el.calibration_failed_beep(2);
                d=el.calibration_failed_beep(3);
            case 'calibration_success_beep'
                doBeep=el.feedbackbeep;
                f=el.calibration_success_beep(1);
                v=el.calibration_success_beep(2);
                d=el.calibration_success_beep(3);
            case 'drift_correction_failed_beep'
                doBeep=el.feedbackbeep;
                f=el.drift_correction_failed_beep(1);
                v=el.drift_correction_failed_beep(2);
                d=el.drift_correction_failed_beep(3);
            case 'drift_correction_success_beep'
                doBeep=el.feedbackbeep;
                f=el.drift_correction_success_beep(1);
                v=el.drift_correction_success_beep(2);
                d=el.drift_correction_success_beep(3);
            otherwise
                % some defaults
                doBeep=el.feedbackbeep;
                f=500;
                v=0.5;
                d=1.5;
        end
        
        if doBeep==1
            if PsychPortAudio('GetOpenDeviceCount') > 0
                warning(sprintf([ 'EyelinkToolbox - ''el.feedbackbeep'' or ''el.targetbeep'' not 0, and a PsychPortAudio\n' ...
                    'device is open. Disabling audio beeps from EyelinkToolbox to avoid conflict.\n' ...
                       ...
                       ]));
            else
                Beeper(f, v, d);
                Snd('Close');
            end
        end
    end
end



