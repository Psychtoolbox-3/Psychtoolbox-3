function DriftDemo5(angle, cyclespersecond, f, drawmask)
% function DriftDemo5(angle, cyclespersecond, f, drawmask)
% ___________________________________________________________________
%
% Display animated gratings using the new Screen('DrawTexture') command.
%
% The demo shows two drifting sine gratings through circular apertures. The
% 1st drifting grating is surrounded by an annulus (a ring) that shows a
% second drifting grating with a different orientation.
%
% The demo ends after a key press or after 20 seconds have elapsed.
%
% The demo uses alpha-blending and color buffer masking for masking the
% gratings with circular apertures.
%
% Parameters:
%
% angle = Angle of the grating with respect to the vertical direction.
% cyclespersecond = Speed of grating in cycles per second. f = Frequency of
% grating in cycles per pixel.
% drawmask = If set to 1, then a gaussian aperture is drawn over the grating
% _________________________________________________________________________
%
% see also: PsychDemos, MovieDemo

% HISTORY
% 4/1/09 mk Adapted from Allen Ingling's DriftDemo.m

if nargin<4
    % By default, we mask the grating by a gaussian transparency mask:
    drawmask=1;
end;

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

movieDurationSecs=20; % Abort demo after 20 seconds.
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
    [w screenRect]=Screen('OpenWindow',screenNumber, gray);

    % Calculate parameters of the grating:
    p=ceil(1/f); % pixels/cycle, rounded up.
    fr=f*2*pi;
    visiblesize=2*texsize+1;
    visible2size=visiblesize/2;

    % Create one single static grating image:
    % MK: We only need a single texture row (i.e. 1 pixel in height) to
    % define the whole grating! If srcRect in the Drawtexture call below is
    % "higher" than that (i.e. visibleSize >> 1), the GPU will
    % automatically replicate pixel rows. This 1 pixel height saves memory
    % and memory bandwith, ie. potentially faster.
    [x,y]=meshgrid(-texsize:texsize + p, 1);
    grating=gray + inc*cos(fr*x);

    [x2,y2]=meshgrid(-texsize/2:texsize/2 + p, 1);
    m=cos(fr*x2);
    grating2=gray + inc*m;


    % Store grating in texture:
    gratingtex=Screen('MakeTexture', w, grating);
    grating2tex=Screen('MakeTexture', w, grating2);

    % Create a single gaussian transparency mask and store it to a texture:
    mask=ones(2*texsize+1, 2*texsize+1, 2) * gray;
    [x,y]=meshgrid(-1*texsize:1*texsize,-1*texsize:1*texsize);
    mask(:, :, 2)=white * (1-(x.^2 + y.^2 <= texsize^2));
    masktex=Screen('MakeTexture', w, mask);

    % Definition of the drawn rectangle on the screen:
    dstRect=[0 0 visiblesize visiblesize];
    dstRect=CenterRect(dstRect, screenRect);

    % Definition of the drawn rectangle on the screen:
    dst2Rect=[0 0 visible2size visible2size];
    dst2Rect=CenterRect(dst2Rect, screenRect);

    % Query duration of monitor refresh interval:
    ifi=Screen('GetFlipInterval', w);

    waitframes = 1;
    waitduration = waitframes * ifi;

    % Recompute p, this time without the ceil() operation from above.
    % Otherwise we will get wrong drift speed due to rounding!
    p=1/f; % pixels/cycle

    % Translate requested speed of the grating (in cycles per second)
    % into a shift value in "pixels per frame", assuming given
    % waitduration: This is the amount of pixels to shift our "aperture" at
    % each redraw:
    shiftperframe= cyclespersecond * p * waitduration;

    % Perform initial Flip to sync us to the VBL and for getting an initial
    % VBL-Timestamp for our "WaitBlanking" emulation:
    vbl=Screen('Flip', w);

    % We run at most 'movieDurationSecs' seconds if user doesn't abort via
    % keypress.
    vblendtime = vbl + movieDurationSecs;
    i=0;

    % Animationloop:
    while (vbl < vblendtime) && ~KbCheck

        % Shift the grating by "shiftperframe" pixels per frame:
        xoffset = mod(i*shiftperframe,p);
        i=i+1;

        % Define shifted srcRect that cuts out the properly shifted rectangular
        % area from the texture:
        srcRect=[xoffset 0 xoffset + visiblesize visiblesize];
        src2Rect=[xoffset 0 xoffset + visible2size visible2size];

        % Draw grating texture, rotated by "angle":
        Screen('DrawTexture', w, gratingtex, srcRect, dstRect, angle);

        if drawmask==1
            % Draw aperture over grating:
            Screen('DrawTexture', w, masktex, [0 0 visiblesize visiblesize], dstRect, angle);
        end;

        % Disable alpha-blending, restrict following drawing to alpha channel:
        Screen('Blendfunction', w, GL_ONE, GL_ZERO, [0 0 0 1]);
        
        % Clear 'dstRect' region of framebuffers alpha channel to zero:
        Screen('FillRect', w, [0 0 0 0], dst2Rect);
        
        % Fill circular 'dstRect' region with an alpha value of 255:
        Screen('FillOval', w, [0 0 0 255], dst2Rect);
        
        % Enable DeSTination alpha blending and reenalbe drawing to all
        % color channels. Following drawing commands will only draw there
        % the alpha value in the framebuffer is greater than zero, ie., in
        % our case, inside the circular 'dst2Rect' aperture where alpha has
        % been set to 255 by our 'FillOval' command:
        Screen('Blendfunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]);
        
        % Draw 2nd grating texture, but only inside alpha == 255 circular
        % aperture, and at an angle of 90 degrees:
        Screen('DrawTexture', w, grating2tex, src2Rect, dst2Rect, 90);
        
        % Restore alpha blending mode for next draw iteration:
        Screen('Blendfunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        % Flip 'waitframes' monitor refresh intervals after last redraw.
        vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
    end;

    Priority(0);
    Screen('CloseAll');

catch
    %this "catch" section executes in case of an error in the "try" section
    %above. Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end %try..catch..
