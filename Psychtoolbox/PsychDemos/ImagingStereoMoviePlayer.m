function ImagingStereoMoviePlayer(moviefile, stereoMode, imaging, anaglyphmode, screenid)
% ImagingStereoMoviePlayer(moviefile [,stereoMode=8] [,imaging=1] [,anaglyphmode=0] [,screenid=max])
%
% Minimalistic movie player for stereo movies. Reads movie from file
% 'moviefile'. Left half of each movie video frame must contain left-eye
% image, whereas right half of each frame must contain right-eye image.
%
% 'stereoMode' mode of presentation, defaults to mode 8 (Red-Blue
% Anaglyph). 'imaging' if set to 1, will use the Psychtoolbox imaging
% pipeline for stereo display -- allows to set gains for anaglyph stereo
% and other more advanced anaglyph algorithms.
%
% 'anaglyphmode' when imaging is enabled, allows to select the type of
% anaglyph algorithm:
%
% 0 = Standard anaglyphs.
% 1 = Gray anaglyphs.
% 2 = Half-color anaglyphs.
% 3 = Optimized color anaglyphs.
% 4 = Full color anaglyphs.
%
% See "help SetAnaglyphParameters" for further description and references.
%
% 'screenid' Screen id of target display screen (on multi-display setups).
% By default, the screen with maximum id is used.
%
% The left image is centered on the screen, the right images position can
% be moved by moving the mouse cursor to align for inter-eye distance.
%
% Press any key to quit the player.

% History:
% 11.11.2007 Written (MK)
% 17.06.2013 Cleaned up (MK)

AssertOpenGL;

if nargin < 1
    error('You must at least provide the name of the movie file for stereo pair.');
end

if nargin < 2
    stereoMode = [];
end

if isempty(stereoMode)
    stereoMode = 8;
end

if nargin < 3
    imaging = [];
end

if isempty(imaging)
    imaging = 1;
end

if nargin < 4
    anaglyphmode = [];
end

if isempty(anaglyphmode)
    anaglyphmode = 0;
end

if nargin < 5
    screenid = [];
end

if isempty(screenid)
    screenid = max(Screen('Screens'));
end

if imaging
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseVirtualFramebuffer');
    [win, winRect] = PsychImaging('OpenWindow', screenid, 0, [], [], [], stereoMode);
else
    [win, winRect] = Screen('OpenWindow', screenid, 0, [], [], [], stereoMode);
end
modestr = [];

if imaging
    % Set color gains. This depends on the anaglyph mode selected:
    switch stereoMode
        case 6,
            SetAnaglyphStereoParameters('LeftGains', win,  [1.0 0.0 0.0]);
            SetAnaglyphStereoParameters('RightGains', win, [0.0 0.6 0.0]);
        case 7,
            SetAnaglyphStereoParameters('LeftGains', win,  [0.0 0.6 0.0]);
            SetAnaglyphStereoParameters('RightGains', win, [1.0 0.0 0.0]);
        case 8,
            SetAnaglyphStereoParameters('LeftGains', win, [0.4 0.0 0.0]);
            SetAnaglyphStereoParameters('RightGains', win, [0.0 0.2 0.7]);
        case 9,
            SetAnaglyphStereoParameters('LeftGains', win, [0.0 0.2 0.7]);
            SetAnaglyphStereoParameters('RightGains', win, [0.4 0.0 0.0]);
        otherwise
            %error('Unknown stereoMode specified.');
    end
    
    if stereoMode > 5 && stereoMode < 10
        switch anaglyphmode
            case 0,
                % Default anaglyphs, nothing to do...
                modestr = 'Standard anaglyphs';
            case 1,
                SetAnaglyphStereoParameters('GrayAnaglyphMode', win);
                modestr = 'Gray anaglyph rendering';
            case 2,
                SetAnaglyphStereoParameters('HalfColorAnaglyphMode', win);
                modestr = 'Half color anaglyph rendering';
            case 3,
                SetAnaglyphStereoParameters('OptimizedColorAnaglyphMode', win);
                modestr = 'Optimized color anaglyph rendering';
            case 4,
                SetAnaglyphStereoParameters('FullColorAnaglyphMode', win);
                modestr = 'Full color anaglyph rendering';
            otherwise
                error('Invalid anaglyphmode specified!');
        end
        
        overlay = SetAnaglyphStereoParameters('CreateMonoOverlay', win);
        Screen('TextSize', overlay, 24);
        DrawFormattedText(overlay, ['Loading file: ' moviefile ], 0, 0, [255 255 0]);
    end
end

% Initial flip:
Screen('Flip', win);

% Open movie file and start playback:
movie = Screen('OpenMovie', win, moviefile);
Screen('PlayMovie', movie, 1, 1, 1);

% Position mouse on center of display:
[x , y] = RectCenter(winRect);
SetMouse(x, y, win);

% Hide mouse cursor:
HideCursor;

% Setup variables:
tex = 0;
imgrect = [];

if ~isempty(modestr)
    Screen('FillRect', overlay, [0 0 0 0]);
    DrawFormattedText(overlay, ['File: ' moviefile '\nOpmode: ' modestr], 0, 0, [255 255 0]);
end

try
    % Playback loop: Run until keypress or error:
    while ~KbCheck && tex~=-1
        
        % Fetch next image from movie:
        tex = Screen('GetMovieImage', win, movie, 1);
        
        % Valid image to draw?
        if tex>0
            % Query mouse position:
            x = GetMouse(win);
            
            % Setup drawing regions based on size of first frame:
            if isempty(imgrect)
                imgrect = Screen('Rect', tex) ;
                imglrect = [0, 0, RectWidth(imgrect)/2, RectHeight(imgrect)];
                imgrrect = [RectWidth(imgrect)/2, 0, RectWidth(imgrect), RectHeight(imgrect)];
                sf = min([RectWidth(winRect)/RectWidth(imglrect) , RectHeight(winRect)/RectHeight(imglrect)]);
                dstrect = ScaleRect(imglrect,sf,sf);
            end
            
            % Left eye image == left half of movie texture:
            Screen('SelectStereoDrawBuffer', win, 0);
            Screen('DrawTexture', win, tex, imglrect, CenterRect(dstrect, winRect));
            
            Screen('SelectStereoDrawBuffer', win, 1);
            % Draw right image centered on mouse position -- mouse controls image
            % offsets:
            Screen('DrawTexture', win, tex, imgrrect, CenterRectOnPoint(dstrect, x, y));
            
            % Show at next retrace:
            Screen('Flip', win);
            
            % Release old image texture:
            Screen('Close', tex);
            tex = 0;
        end
        % Next frame...
    end
    
    % Done with playback:
    
    % Show mouse cursor:
    ShowCursor;
    
    % Stop and close movie:
    Screen('PlayMovie', movie, 0);
    Screen('CloseMovie', movie);
    
    % Close screen:
    Screen('CloseAll');
    
    return;
    
catch %#ok<CTCH>
    sca;
    psychrethrow(psychlasterror);
end
