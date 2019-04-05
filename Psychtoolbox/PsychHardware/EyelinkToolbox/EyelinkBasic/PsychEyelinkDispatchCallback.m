function rc = PsychEyelinkDispatchCallback(callArgs, msg)
% Retrieve live eye-image from Eyelink, show it in onscreen window.
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
% to actually receive and display the images, register this function as eyelink's callback:
% if Eyelink('Initialize', 'myEyelinkDispatchCallback') ~=0
%   error('eyelink failed init')
% end
% result = Eyelink('StartSetup',1) %put the tracker into a mode capable of sending images
% then you must hit 'return' on the PTB computer, this key command will be sent to the tracker host to initiate sending of images.
%
% This function fetches the most recent live image from the Eylink eye
% camera and displays it in the previously assigned onscreen window.
%
% History:
% 15.3.2009   Derived from MemoryBuffer2TextureDemo.m (MK).
%  4.4.2009   Updated to use EyelinkGetKey + fixed eyelinktex persistence crash (edf).
% 11.4.2009   Cleaned up. Should be ready for 1st release, although still
%             pretty alpha quality. (MK).
% 15.6.2010   Added some drawing routines to get standard behaviour back. Enabled
%             use of the callback by default. Clarified in helptext that user
%             normally should not have to worry about calling this file. (fwc)
% 20.7.2010   drawing of instructions, eye-image+title, playing sounds in seperate functions
%
%  1.2.2010   modified to allow for cross hair and fix bugs. (nj)
% 29.10.2018  Drop 'DrawDots' for calibration target. Some white-space fixes.

% Cached texture handle for eyelink texture:
persistent eyelinktex;
global dw dh offscreen;

% Cached window handle for target onscreen window:
persistent eyewin;
persistent calxy;
persistent imgtitle;
persistent eyewidth;
persistent eyeheight;

% Cached(!) eyelink stucture containing keycodes
persistent el;
persistent lastImageTime; %#ok<PUSE>
persistent drawcount;
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
    eyelinktex = []; % got persisted from a previous ptb window which has now been closed; needs to be recreated
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
    drawcount = 0;
    lastImageTime = GetSecs;
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
    lastImageTime=GetSecs;
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

if isempty(eyewin)
    warning('Got called as callback function from Eyelink() but usercode has not set a valid target onscreen window handle yet! Aborted.'); %#ok<WNTAG>
    return;
end

% Flag that tells if a new camera image was received and our camera image texture needs update:
newcamimage = 0;
needsupdate = 0;

switch eyecmd
    case 1,
        % New videoframe received. See code below for actual processing.
        newcamimage = 1;
        needsupdate = 1;
    case 2,
        % Eyelink Keyboard query:
        [rc, el] = EyelinkGetKey(el);
    case 3,
        % Alert message:
        fprintf('Eyelink Alert: %s.\n', msg);
        needsupdate = 1;
    case 4,
        % Image title of camera image transmitted from Eyelink:
        % fprintf('Eyelink image title is %s. [Threshold = %f]\n', msg, callArgs(2));
        if callArgs(2) ~= -1
            imgtitle = sprintf('Camera: %s [Threshold = %f]', msg, callArgs(2));
        else
            imgtitle = msg;
        end
        needsupdate = 1;
    case 5,
        % Define calibration target and enable its drawing:
        % fprintf('draw_cal_target.\n');
        calxy = callArgs(2:3);
        clearScreen=1;
        needsupdate = 1;
    case 6,
        % Clear calibration display:
        % fprintf('clear_cal_display.\n');
        clearScreen=1;
        drawInstructions=1;
        needsupdate = 1;
    case 7,
        % Setup calibration display:
        if inDrift
            drawInstructions = 0;
            inDrift = 0;
        else
            drawInstructions = 1;
        end

        clearScreen=1;
        % drawInstructions=1;
        drawcount = 0;
        lastImageTime = GetSecs;
        needsupdate = 1;
    case 8,
        newImage = 1;
        % Setup image display:
        eyewidth  = callArgs(2);
        eyeheight = callArgs(3);
        % fprintf('setup_image_display for %i x %i pixels.\n', eyewidth, eyeheight);
        drawcount = 0;
        lastImageTime = GetSecs;
        ineyeimagemodedisplay=1;
        drawInstructions=1;
        needsupdate = 1;
    case 9,
        % Exit image display:
        %         fprintf('exit_image_display.\n');
        %         fprintf('AVG FPS = %f Hz\n', drawcount / (GetSecs - lastImageTime));
        clearScreen=1;
        ineyeimagemodedisplay=0;
        drawInstructions=1;
        needsupdate = 1;
    case 10,
        % Erase current calibration target:
        %         fprintf('erase_cal_target.\n');
        calxy = [];
        clearScreen=1;
        needsupdate = 1;
    case 11,
        %               fprintf('exit_cal_display.\n');
        %         fprintf('AVG FPS = %f Hz\n', drawcount / (GetSecs - lastImageTime));
        clearScreen=1;
        %       drawInstructions=1;
        needsupdate = 1;
    case 12,
        % New calibration target sound:
        %         fprintf('cal_target_beep_hook.\n');
        EyelinkMakeSound(el, 'cal_target_beep');
    case 13,
        % New drift correction target sound:
        %         fprintf('dc_target_beep_hook.\n');
        EyelinkMakeSound(el, 'drift_correction_target_beep');
    case 14,
        % Calibration done sound:
        errc = callArgs(2);
        %         fprintf('cal_done_beep_hook: %i\n', errc);
        if errc > 0
            % Calibration failed:
            EyelinkMakeSound(el, 'calibration_failed_beep');
        else
            % Calibration success:
            EyelinkMakeSound(el, 'calibration_success_beep');
        end
    case 15,
        % Drift correction done sound:
        errc = callArgs(2);
        %         fprintf('dc_done_beep_hook: %i\n', errc);
        if errc > 0
            % Drift correction failed:
            EyelinkMakeSound(el, 'drift_correction_failed_beep');
        else
            % Drift correction success:
            EyelinkMakeSound(el, 'drift_correction_success_beep');
        end
        % add by NJ
    case 16,
        [width, height]=Screen('WindowSize', eyewin);
        % get mouse
        [x,y, buttons] = GetMouse(eyewin);

        HideCursor
        if find(buttons)
            rc = [width , height, x , y,  dw , dh , 1];
        else
            rc = [width , height, x , y , dw , dh , 0];
        end
         % add by NJ to prevent flashing of text in drift correct
    case 17,
        inDrift =1;
    otherwise
        % Unknown command:
        fprintf('PsychEyelinkDispatchCallback: Unknown eyelink command (%i)\n', eyecmd);
        return;
end

% Display redraw and update needed?
if ~needsupdate
    % Nope. Return from callback:
    return;
end

% Need to rebuild/redraw and flip the display:
% need to clear screen?
if clearScreen==1
    Screen('FillRect', eyewin, el.backgroundcolour);
    clearScreen=0;
end
% New video data from eyelink?
if newcamimage
    % Video callback from Eyelink: We have a 'eyewidth' by 'eyeheight' pixels
    % live eye image from the Eyelink system. Each pixel is encoded as a 4 byte
    % RGBA pixel with alpha channel set to a constant value of 255 and the RGB
    % channels encoding a 1-Byte per channel R, G or B color value. The
    % given 'eyeimgptr' is a specially encoded memory pointer to the memory
    % buffer inside Eyelink() that encodes the image.
    eyeimgptr = callArgs(2);
    eyewidth  = callArgs(3);
    eyeheight = callArgs(4);
    
    % Create a new PTB texture of proper format and size and inject the 4
    % channel RGBA color image from the Eyelink memory buffer into the texture.
    % Return a standard PTB texture handle to it. If such a texture already
    % exists from a previous invocation of this routiene, just recycle it for
    % slightly higher efficiency:
    eyelinktex = Screen('SetOpenGLTextureFromMemPointer', eyewin, eyelinktex, eyeimgptr, eyewidth, eyeheight, 4, 0, [], GL_RGBA8, GL_RGBA, hostDataFormat);
end

% If we're in imagemodedisplay, draw eye camera image texture centered in
% window, if any such texture exists, also draw title if it exists.
if ~isempty(eyelinktex) && ineyeimagemodedisplay==1
    imgtitle=EyelinkDrawCameraImage(eyewin, el, eyelinktex, imgtitle,newImage);
end

% Draw calibration target, if any is specified:
if ~isempty(calxy)
    drawInstructions=0;
    EyelinkDrawCalibrationTarget(eyewin, el, calxy);
end

% Need to draw instructions?
if drawInstructions==1
    EyelinkDrawInstructions(eyewin, el,msg);
    drawInstructions=0;
end

% Show it: We disable synchronization of Matlab to the vertical retrace.
% This way, display update itself is still synced and tear-free, but we
% don't waste time waiting for swap completion. Potentially higher
% performance for calibration displays and eye camera image updates...
% Neither do we erase buffer
Screen('Flip', eyewin, [], 1, 1);

% Some counter, just to measure update rate:
drawcount = drawcount + 1;

% Done. Return from callback:
return;


function EyelinkDrawInstructions(eyewin, el,msg)
oldFont=Screen(eyewin,'TextFont',el.msgfont);
oldFontSize=Screen(eyewin,'TextSize',el.msgfontsize);
DrawFormattedText(eyewin, el.helptext, 20, 20, el.msgfontcolour, [], [], [], 1);

if el.displayCalResults && ~isempty(msg)
    DrawFormattedText(eyewin, msg, 20, 150, el.msgfontcolour, [], [], [], 1);
end

Screen(eyewin,'TextFont',oldFont);
Screen(eyewin,'TextSize',oldFontSize);

function  imgtitle=EyelinkDrawCameraImage(eyewin, el, eyelinktex, imgtitle,newImage)
persistent lasttitle;
global dh dw offscreen;
try
    if ~isempty(eyelinktex)
        eyerect=Screen('Rect', eyelinktex);
        % we could cash some of the below values....
        wrect=Screen('Rect', eyewin);
        [width, heigth]=Screen('WindowSize', eyewin);
        dw=round(el.eyeimgsize/100*width);
        dh=round(dw * eyerect(4)/eyerect(3));

        drect=[ 0 0 dw dh ];
        drect=CenterRect(drect, wrect);
        Screen('DrawTexture', eyewin, eyelinktex, [], drect);
        %     fprintf('EyelinkDrawCameraImage:DrawTexture \n');
    end
    % imgtitle
    % if title is provided, we also draw title
    if ~isempty(eyelinktex) && exist( 'imgtitle', 'var') && ~isempty(imgtitle)
        %oldFont=Screen(eyewin,'TextFont',el.imgtitlefont);
        %oldFontSize=Screen('TextSize',eyewin,el.imgtitlefontsize);
        rect=Screen('TextBounds', eyewin, imgtitle );
        [w2, h2]=RectSize(rect);

        % added by NJ as a quick way to prevent over drawing and to clear text
        if newImage || isempty(lasttitle) || ~strcmp(imgtitle,lasttitle)
            if -1 == Screen('WindowKind', offscreen)
                Screen('Close', offscreen);
            end

            sn = Screen('WindowScreenNumber', eyewin);
            offscreen = Screen('OpenOffscreenWindow', sn, el.backgroundcolour, [], [], 32);
            Screen(offscreen,'TextFont',el.imgtitlefont);
            Screen(offscreen,'TextSize',el.imgtitlefontsize);
            Screen('DrawText', offscreen, imgtitle, width/2-dw/2, heigth/2+dh/2+h2, el.imgtitlecolour);
            Screen('DrawTexture',eyewin,offscreen,  [width/2-dw/2 heigth/2+dh/2+h2 width/2-dw/2+500 heigth/2+dh/2+h2+500], [width/2-dw/2 heigth/2+dh/2+h2 width/2-dw/2+500 heigth/2+dh/2+h2+500]);
            Screen('Close',offscreen);
            newImage = 0;
        end
        %imgtitle=[]; % return empty title, so it doesn't get drawn over and over again.
        lasttitle = imgtitle;
    end
catch %myerr
    fprintf('EyelinkDrawCameraImage:error \n');
    disp(psychlasterror);
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
    Beeper(f, v, d);
end

function EyelinkDrawCalibrationTarget(eyewin, el, calxy)
[width, heigth]=Screen('WindowSize', eyewin);
size=round(el.calibrationtargetsize/100*width);
inset=round(el.calibrationtargetwidth/100*width);

insetSize = size-2*inset;
if insetSize < 1
    insetSize = 1;
end

% Use FillOval for larger dots:
Screen('FillOval', eyewin, el.calibrationtargetcolour, [calxy(1)-size/2 calxy(2)-size/2 calxy(1)+size/2 calxy(2)+size/2], size+2);
Screen('FillOval', eyewin, el.backgroundcolour, [calxy(1)-inset/2 calxy(2)-inset/2 calxy(1)+inset/2 calxy(2)+inset/2], inset+2);
