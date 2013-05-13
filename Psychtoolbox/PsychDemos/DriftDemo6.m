function DriftDemo6(angle, cyclespersecond, f)
% function DriftDemo6(angle, cyclespersecond, f)
% ___________________________________________________________________
%
% This demo demonstrates how to use Screen('DrawTexture') in combination
% with GLSL texture draw shaders to efficiently combine two drifting
% gratings with each other.
%
% The demo shows a drifting sine grating through a circular aperture. The
% drifting grating is surrounded by an annulus (a ring) that shows a second
% drifting grating with a different orientation.
%
% The demo ends after a key press or after 60 seconds have elapsed.
%
% The demo needs modern graphics hardware with shading support, otherwise
% it will abort.
%
% Compared to previous demos, we apply the aperture to the grating texture
% while drawing the grating texture, ie. in a single drawing pass, instead
% of applying it in a 2nd pass after the grating has been drawn already.
% This is simpler and faster than the dual-pass method. For this, we store
% the grating pattern in the luminance channel of a single texture, and the
% alpha-mask in the alpha channel of *the same texture*. During drawing, we
% apply a special texture filter shader (created via
% MakeTextureDrawShader()). This shader allows to treat the alpha channel
% separate from the luminance or rgb channels of a texture: It applies the
% alpha channel "as is", but applies some shift to the luminance or rgb
% channels of the texture.
%
% The procedure is repeated with a 2nd masked texture to create two
% different drifting gratings, superimposed to each other.
%
% Please note that the same effect can be achieved by clever alpha blending
% on older hardware, e.g., see DriftDemo5. The point of this demo is to
% demonstrate how to use GLSL shaders for more efficient ways of
% manipulating textures during drawing.
%
% Parameters:
%
% angle = Angle of the gratings with respect to the vertical direction.
% cyclespersecond = Speed of gratings in cycles per second.
% f = Frequency of gratings in cycles per pixel.
%
% _________________________________________________________________________
%
% see also: PsychDemos, MovieDemo

% HISTORY
% 3/31/09 mk Written.

if nargin<3
    % Grating cycles/pixel
    f=0.05;
end;

if nargin<2
    % Speed of grating in cycles per second:
    cyclespersecond=1;
end;

if nargin<1
    % Angle of the grating: We default to 30 degrees.
    angle=30;
end;

movieDurationSecs=60; % Abort demo after 60 seconds.
texsize=300; % Half-Size of the grating image.

try
    AssertOpenGL;

    % Get the list of screens and choose the one with the highest screen number.
    screenNumber=max(Screen('Screens'));

    % Find the color values which correspond to white and black.
    white=WhiteIndex(screenNumber);
    black=BlackIndex(screenNumber);

    % Round gray to integral number, to avoid roundoff artifacts with some
    % graphics cards:
	gray=round((white+black)/2);

    % This makes sure that on floating point framebuffers we still get a
    % well defined gray. It isn't strictly neccessary in this demo:
    if gray == white
		gray=white / 2;
    end
    inc=white-gray;

    % Open a double buffered fullscreen window with a gray background:
    w =Screen('OpenWindow',screenNumber, gray);

    % Make sure this GPU supports shading at all:
    AssertGLSL;
    
    % Enable alpha blending for typical drawing of masked textures:
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % Create a special texture drawing shader for masked texture drawing:
    glsl = MakeTextureDrawShader(w, 'SeparateAlphaChannel');

    % Calculate parameters of the grating:
    p=ceil(1/f); % pixels/cycle, rounded up.
    fr=f*2*pi;
    visiblesize=2*texsize+1;

    % Create one single static grating image:
    x = meshgrid(-texsize:texsize + p, -texsize:texsize);
    grating = gray + inc*cos(fr*x);

    % Create circular aperture for the alpha-channel:
    [x,y]=meshgrid(-texsize:texsize, -texsize:texsize);
    circle = white * (x.^2 + y.^2 <= (texsize)^2);

    % Set 2nd channel (the alpha channel) of 'grating' to the aperture
    % defined in 'circle':
    grating(:,:,2) = 0;
    grating(1:2*texsize+1, 1:2*texsize+1, 2) = circle;

    % Store alpha-masked grating in texture and attach the special 'glsl'
    % texture shader to it:
    gratingtex1 = Screen('MakeTexture', w, grating , [], [], [], [], glsl);

    % Build a second drifting grating texture, this time half the texsize
    % of the 1st texture:
    texsize = ceil(texsize/2);
    visible2size = 2*texsize+1;
    x = meshgrid(-texsize:texsize + p, -texsize:texsize);
    grating = gray + inc*cos(fr*x);

    % Create circular aperture for the alpha-channel:
    [x,y]=meshgrid(-texsize:texsize, -texsize:texsize);
    circle = white * (x.^2 + y.^2 <= (texsize)^2);

    % Set 2nd channel (the alpha channel) of 'grating' to the aperture
    % defined in 'circle':
    grating(:,:,2) = 0;
    grating(1:2*texsize+1, 1:2*texsize+1, 2) = circle;

    % Store alpha-masked grating in texture and attach the special 'glsl'
    % texture shader to it:
    gratingtex2 = Screen('MakeTexture', w, grating, [], [], [], [], glsl);
    
    % Definition of the drawn source rectangle on the screen:
    srcRect=[0 0 visiblesize visiblesize];

    % Definition of the drawn source rectangle on the screen:
    src2Rect=[0 0 visible2size visible2size];

    % Query duration of monitor refresh interval:
    ifi=Screen('GetFlipInterval', w);

    waitframes = 1;
    waitduration = waitframes * ifi;

    % Recompute p, this time without the ceil() operation from above.
    % Otherwise we will get wrong drift speed due to rounding!
    p = 1/f; % pixels/cycle

    % Translate requested speed of the gratings (in cycles per second) into
    % a shift value in "pixels per frame", assuming given waitduration:
    shiftperframe = cyclespersecond * p * waitduration;

    % Perform initial Flip to sync us to the VBL and for getting an initial
    % VBL-Timestamp for our "WaitBlanking" emulation:
    vbl = Screen('Flip', w);

    % We run at most 'movieDurationSecs' seconds if user doesn't abort via keypress.
    vblendtime = vbl + movieDurationSecs;
    i=0;

    % Animation loop: Run until timeout or keypress.
    while (vbl < vblendtime) && ~KbCheck
        
        % Shift the grating by "shiftperframe" pixels per frame. We pass
        % the pixel offset 'yoffset' as a parameter to
        % Screen('DrawTexture'). The attached 'glsl' texture draw shader
        % will apply this 'yoffset' pixel shift to the RGB or Luminance
        % color channels of the texture during drawing, thereby shifting
        % the gratings. Before drawing the shifted grating, it will mask it
        % with the "unshifted" alpha mask values inside the Alpha channel:
        yoffset = mod(i*shiftperframe,p);
        i=i+1;

        % Draw first grating texture, rotated by "angle":
        Screen('DrawTexture', w, gratingtex1, srcRect, [], angle, [], [], [], [], [], [0, yoffset, 0, 0]);

        % Draw 2nd grating texture, rotated by "angle+45":
        Screen('DrawTexture', w, gratingtex2, src2Rect, [], angle+45, [], [], [], [], [], [0, yoffset, 0, 0]);

        % Flip 'waitframes' monitor refresh intervals after last redraw.
        vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
    end;

    % The same commands wich close onscreen and offscreen windows also close textures.
    Screen('CloseAll');

catch
    % This "catch" section executes in case of an error in the "try" section
    % above. Importantly, it closes the onscreen window if it is open.
    Screen('CloseAll');
    psychrethrow(psychlasterror);

end %try..catch..
