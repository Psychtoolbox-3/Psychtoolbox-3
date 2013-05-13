function GLSLCLUTAnimDemo
%
% Demonstrates GLSL shader based CLUT animation via moglClutblit()
%
% THIS DEMO IS DEPRECATED! For a better and simpler method, look at
% ClutAnimDemo(2); and the help text of PsychImaging() for the
% 'EnableCLUTMapping' function.
%
% see also: ClutAnimDemo, PsychImaging

% HISTORY
%  7/05/05    mk     Wrote it. 
%  22/0705  fwc       added SkipSyncTests preference call, slightly smaller texture,
%                   as drawtexture failed on s=400 on a 1024x768 pix screen
%                   in catch section, test if OrigLut exists before
%                   applying it.

try
    
    help GLSLCLUTAnimDemo;
    
    % This script calls Psychtoolbox commands available only in OpenGL-based 
	% versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is the
	% only OpenGL-base Psychtoolbox.)  The Psychtoolbox command AssertPsychOpenGL will issue
	% an error message if someone tries to execute this script on a computer without
	% an OpenGL Psychtoolbox
	AssertOpenGL;
	
	% Get the list of screens and choose the one with the highest screen number.
	% Screen 0 is, by definition, the display with the menu bar. Often when 
	% two monitors are connected the one without the menu bar is used as 
	% the stimulus display.  Chosing the display with the highest dislay number is 
	% a best guess about where you want the stimulus displayed.  
	screens=Screen('Screens');
	screenNumber=max(screens);

    % Enable OpenGL mode of Psychtoolbox: This is crucially needed for clut
    % animation:
    InitializeMatlabOpenGL;
    
    % Open a double buffered fullscreen window with black background:
	w=Screen('OpenWindow',screenNumber, 0);
    
	% Find the color value which corresponds to black.  Though on OS
	% X we currently only support true color and thus, for scalar color
	% arguments, black is always 0 and white 255, this rule may not be true
    % on other platforms.
	black=BlackIndex(screenNumber);
	
    % Build a simple gray-level ramp as a single texture.
    [width, height]=Screen('WindowSize', w);
    s=floor(min(width, height)/2)-1;
    [x,y]=meshgrid(-s:s, -s:s);
    fintex=ones(2*s+1,2*s+1);
    fintex(:,:)=mod(x,256);
    tex=Screen('MakeTexture', w, fintex);

    % Black background:
    Screen('FillRect',w, black);
    
    newLUT = zeros(256,3);
    for i=0:255
        newLUT(i+1, :)=[i i i];
    end

    i=0;
    tavg=0;
    t0 = GetSecs;
    
    % Show it by flipping the buffers:
    ifi = Screen('GetFlipInterval', w);
    vbl = Screen('Flip', w);

    % Animation by CLUT color cycling loop:
	while (1)
        % Shift/Cycle all LUT entries: Entry 3 -> 2, 4 -> 3, 5 ->4 , ... ,
        % 256 -> 255, 2 -> 256, ... we just leave slot 1 alone, it defines
        % the DAC output values for the background.
        backupLUT=newLUT(1, :);
        newLUT(1:255, :)=newLUT(2:256, :);
        newLUT(256, :)=backupLUT;

        % Perform blit of our image, applying newLUT as clut:
        moglClutBlit(w, tex, newLUT);
        Screen('DrawLine', w, [255 255 0], 1280-mod(i, 1280), 0, 1280-mod(i, 1280), 1024);

        % Show new image one monitor refresh after the last one:
        vbl = Screen('Flip', w, vbl);
        
        t1=GetSecs;
        t1=vbl;
        tavg=tavg+(t1-t0);
        t0=t1;
        i=i+1;
        tonset(i)=t1;

        % Abort after 1000 video refresh intervals or on a key-press:
        if KbCheck || (i>1000)
            break;
        end;
    end;

    tavg=tavg / i

    % Disable CLUT blitter. This needs to be done before the call to
    % Screen('CloseAll')!
    moglClutBlit;
    
	Screen('CloseAll');

    plot(diff(tonset)*1000);
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end %try..catch..
