function ShepardZoomDemo
% ShepardZoomDemo
%
% This demo implements kind of a visual version of the Shepard tone
% illusion. This is an infinite zoom, by use of an OpenGL technique called
% "Spectral Texturing", a clever way of blending multiple textures
% with each other by use of alpha-blending...
%
% You can find more information about the technique and stimulus at:
% http://www.kyb.tuebingen.mpg.de/~berger/spectraltexturing/
%
% The "Spectral Texturing" method was published as a sketch at Siggraph 2003:
% http://doi.acm.org/10.1145/965400.965509
% It was also demonstrated at VSS demo night 2004.
%
% Feel free to use the stimulus for experiments, but please kindly tell
% Daniel Berger how you use it: daniel.berger@tuebingen.mpg.de
%
% The original code for this demo has been written by Daniel Berger in 2003.
% The code has been adapted for use with the PsychToolbox by Mario Kleiner
% in 2006.

% Make sure we're running on OpenGL Psychtoolbox, abort if not:
AssertOpenGL;

% Setup keycodes for our control keys:
KbName('UnifyKeyNames');
rotccw=KbName('LeftArrow');
rotcw=KbName('RightArrow');
decrease=KbName('DownArrow');
increase=KbName('UpArrow');
escape=KbName('ESCAPE');
toggle(1)=KbName('1!');
toggle(2)=KbName('2@');
toggle(3)=KbName('3#');
toggle(4)=KbName('4$');
toggle(5)=KbName('5%');
helptoggle=KbName('h');

% Initial direction and speed of zoom:
stepper = 1;

% Initial rotation speed:
rotspeed=-0.1;

% Shorten synctests for this demo:
oldlevel = Screen('Preference', 'SkipSyncTests', 1);

% Initialize OpenGL low-level support of Psychtoolbox:
InitializeMatlabOpenGL;

% Open onscreen window on secondary display, with a black
% background color:
screenid = max(Screen('Screens'));
[win winrect] = Screen('OpenWindow', screenid, 0);

% Query size of window. We need it to compute correct scaling for windows
% aspect ratio later on:
[width height] = Screen('WindowSize', win);

ifi=Screen('GetFlipInterval', win);
vbl=Screen('Flip', win);

% Some specs for the zoom...
showlayer =[1,1,1,1,1];
hold=0;
flattransparency=0;

% We read PNG images with color- and alpha channels:
teximage = imread([ PsychtoolboxRoot 'PsychDemos/OpenGL4MatlabDemos/snowy_rgb.png' ]);
alpha    = imread([ PsychtoolboxRoot 'PsychDemos/OpenGL4MatlabDemos/snowy_a.png' ]);

% Make sure our teximage is represented as a 3-layer RGB image, even if
% delivered by imread() as a one layer luminance image: This is a lame way
% of simplifying code later down the road...
if size(teximage,3) == 1
    teximage(:,:,2) = teximage(:,:,1);
    teximage(:,:,3) = teximage(:,:,1);
end

% Add alpha channel as 4th layer to image:
teximage(:,:,4)=alpha(:,:);

if (flattransparency)
    % For this mode, we set alpha to a constant 50%:
    teximage(:,:,4)=128;
end

% For quick texture mipmap creation, the hardware needs to support
% framebuffer object extensions...
if ~isempty(findstr(glGetString(GL.EXTENSIONS), '_framebuffer_object'))
    % Automatic, fast, hardware based creation of mipmap texture resolution
    % pyramids supported on this system. Use it to quickly create the
    % texture mipmap object:

    % Build "power of two" (pot) texture (in this case, a texture of size
    % 2^8*2^8, that is 256 pixels wide and 256 pixels high) from input image:
    tex = Screen('MakeTexture', win, teximage, [], 1);

    % Retrieve OpenGL texture handle to the power of two sized texture:
    texture=Screen('GetOpenGLTexture', win, tex);

    % Activate our own OpenGL rendering code:
    Screen('BeginOpenGL', win);
    
    % Bind the texture for use:
    glBindTexture(GL.TEXTURE_2D, texture);   %  2d texture (x and y size)

    % Ask the hardware to generate all depth levels automatically:
    glGenerateMipmapEXT(GL.TEXTURE_2D);
else
    % OpenGL extension for automatic mipmap texture creation is not
    % supported by this graphics hardware. We need to do it the manual,
    % compute intense way...
    
    % Convert Matlab matrix manually into OpenGL binary image data:
    for i=1:4
        glpixels(i,:,:) = teximage(:,:,i);
    end;
    glpixels = uint8(glpixels);
    
    % Activate our own OpenGL rendering code:
    Screen('BeginOpenGL', win);

    % Generate a texture object:
    texture = glGenTextures(1);
    
    % Bind the texture object for use as power-of-two texture:
    glBindTexture(GL.TEXTURE_2D, texture);

    % Use GLU to compute the image resolution mipmap pyramid and create
    % OpenGL textures ouf of it: This is slooow...
    r = gluBuild2DMipmaps(GL.TEXTURE_2D, 4, size(teximage,1), size(teximage,2), GL.RGBA, GL.UNSIGNED_BYTE, glpixels);
    if r>0
        error('gluBuild2DMipmaps failed for some reason.');
    end
end

% We have our texture, setup its rendering parameters:

% Setup magnification and minification filters:
glTexParameteri(GL.TEXTURE_2D,GL.TEXTURE_MAG_FILTER,GL.LINEAR); %  scale linearly when rendered image bigger than texture
if ~IsGLES
    glTexParameteri(GL.TEXTURE_2D,GL.TEXTURE_MIN_FILTER,GL.LINEAR_MIPMAP_LINEAR); %  scale linearly + mipmap when rendered image smalled than texture
else
    glTexParameteri(GL.TEXTURE_2D,GL.TEXTURE_MIN_FILTER,GL.LINEAR);
end

% Set wrapping behaviour to repeat infinite in any direction:
glTexParameteri(GL.TEXTURE_2D,GL.TEXTURE_WRAP_S,GL.REPEAT);
glTexParameteri(GL.TEXTURE_2D,GL.TEXTURE_WRAP_T,GL.REPEAT);

% Texture colors shall modulate global color value:
glTexEnvf(GL.TEXTURE_ENV, GL.TEXTURE_ENV_MODE, GL.MODULATE);

% Setup a blending equation: Playing with these can be instructive!
blendfg = 6;
blendbg = 7;

switch(blendfg)
  case 0 %  GL.ZERO
    blend_sfactor=GL.ZERO;
  case 1 %  GL.ONE
    blend_sfactor=GL.ONE;
  case 2 %  GL.DST_COLOR
    blend_sfactor=GL.DST_COLOR;
  case 3 % GL.SRC_COLOR
    blend_sfactor=GL.SRC_COLOR;
  case 4 % GL.ONE_MINUS_DST_COLOR
    blend_sfactor=GL.ONE_MINUS_DST_COLOR;
  case 5 % GL.ONE_MINUS_SRC_COLOR
    blend_sfactor=GL.ONE_MINUS_SRC_COLOR;
  case 6 % GL.SRC_ALPHA
    blend_sfactor=GL.SRC_ALPHA;
  case 7 % GL.ONE_MINUS_SRC_ALPHA
    blend_sfactor=GL.ONE_MINUS_SRC_ALPHA;
  case 8 % GL.DST_ALPHA
    blend_sfactor=GL.DST_ALPHA;
  case 9 % GL.ONE_MINUS_DST_ALPHA
    blend_sfactor=GL.ONE_MINUS_DST_ALPHA;
  case 10 % GL.SRC_ALPHA_SATURATE
    blend_sfactor=GL.SRC_ALPHA_SATURATE;
end

switch(blendbg)
  case 0 %  GL.ZERO
    blend_dfactor=GL.ZERO;
  case 1 %  GL.ONE
    blend_dfactor=GL.ONE;
  case 2 %  GL.DST_COLOR
    blend_dfactor=GL.DST_COLOR;
  case 3 % GL.SRC_COLOR
    blend_dfactor=GL.SRC_COLOR;
  case 4 % GL.ONE_MINUS_DST_COLOR
    blend_dfactor=GL.ONE_MINUS_DST_COLOR;
  case 5 % GL.ONE_MINUS_SRC_COLOR
    blend_dfactor=GL.ONE_MINUS_SRC_COLOR;
  case 6 % GL.SRC_ALPHA
    blend_dfactor=GL.SRC_ALPHA;
  case 7 % GL.ONE_MINUS_SRC_ALPHA
    blend_dfactor=GL.ONE_MINUS_SRC_ALPHA;
  case 8 % GL.DST_ALPHA
    blend_dfactor=GL.DST_ALPHA;
  case 9 % GL.ONE_MINUS_DST_ALPHA
    blend_dfactor=GL.ONE_MINUS_DST_ALPHA;
  case 10 % GL.SRC_ALPHA_SATURATE
    blend_dfactor=GL.SRC_ALPHA_SATURATE;
end

% Enable smooth shading:
glShadeModel(GL.SMOOTH);

% Enable alpha blending:
glEnable(GL.BLEND);

% Setup perspective projection for our virtual camera, take displays
% aspect ratio into account:
glMatrixMode(GL.PROJECTION);
glLoadIdentity;
gluPerspective(63.0, width/height,0.5,2000.0);

% Setup camera position and orientation:
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;
yaw=0.0;
pitch=90.0;
roll=0.0;
glRotatef(yaw, 0, 1, 0);
glRotatef(pitch, 1, 0, 0);
glRotatef(roll, 0, 0, 1);
glTranslatef(0,-50,0);

% Setup clear color and clear all buffers:
clearcolor = [0.9 0.9 1 1];
glClearColor(clearcolor(1),clearcolor(2),clearcolor(3),clearcolor(4));
glClear; % Clear depth and color buffer.

% Setup initial settings for layer positions and order...
layercnt=[0,20000,40000,60000,80000];
firstlayer=0;
xo=[0.2, 0.3, 0.4, 0.5, 0.6];
yo=[0.6, 0.5, 0.4, 0.3, 0.2];

% Setup base colors of textures for different layers:
texcolor = ones(8,4);
texcolor(1,:)=[128 100 80 255];
texcolor(2,:)=[255 255 255 200];
texcolor(3,:)=[255 255 255 200];
texcolor(4,:)=[255 255 255 200];
texcolor(5,:)=[255 255 255 200];
texcolor(6,:)=[255 255 255 200];
texcolor(7,:)=[255 255 255 200];
texcolor(8,:)=[255 255 255 200];
texcolor=texcolor/255;

% Setup help text settings:
drawhelp = 1;
helptext = ['Daniel Bergers "Shepard Zoom"\n\n' ...
            'Keys and their meaning:\n\n(h) - Toggle display of this help text\n(1) to (5) - Toggle texture layer 1 to 5\n' ...
            'Left/Right Cursor keys - Rotate ccw/cw\nUp/Down Cursor keys - Zoom in/out\nESCape key - Quit demo\n'];

Screen('TextSize', win, 24);
Screen('TextStyle', win, 1);
Screen('TextColor', win, [255 0 0]);

% Setup initial rotation speed:
rotpos=0.0;
zpos=0;

% Endless animation loop. Will run until ESCape key pressed.
while 1
    % Check keyboard:
    [isdown secs keycode]=KbCheck;
    [mox, moy, buttons]=GetMouse;
    if isdown || any(buttons)
        if keycode(escape) || any(buttons)
            % Exit loop, and thereby the demo.
            break;
        end
        
        if keycode(increase)
            % Increase speed of zoom in:
            stepper=stepper+0.01;
        end
        
        if keycode(decrease)
            % Decrease speed of zoom in:
            stepper=stepper-0.01;
        end

        if keycode(rotcw)
            % Increase clockwise rotation speed:
            rotspeed = rotspeed - 0.001;
        end

        if keycode(rotccw)
            % Decrease clockwise rotation speed:
            rotspeed = rotspeed + 0.001;
        end

        if keycode(helptoggle)
            % Toggle help text display:
            drawhelp = 1 - drawhelp;

            % Debounce key:
            while KbCheck; end;
        end

        for i=1:5
            % Toggle rendering of i'th layer:
            if keycode(toggle(i))
                showlayer(i)=1-showlayer(i);
                % Debounce key:
                while KbCheck; end;
            end
        end
    end

    % Make backup of current modelview matrix, so we can restore it to
    % current state at end of frame drawing:
    glPushMatrix;

    % Render an infinite zoom into or out of the texture by using blending and
    % moving textured layers.

    % Draw lowest layer of groundplane in uniform color:
    glColor4f(texcolor(1,1), texcolor(1,2), texcolor(1,3), texcolor(1,4));
    glDisable(GL.TEXTURE_2D);
    glBlendFunc(GL.ONE,GL.ZERO);
    glBegin(GL.QUADS);
        glVertex3f(-100.0, 0.0, -100.0);     % top left
        glVertex3f(100.0, 0.0, -100.0);  % top right
        glVertex3f(100.0, 0.0, 100.0);  % bottom right
        glVertex3f(-100.0, 0.0, 100.0);  % bottom left
    glEnd;

    % All other layers are drawn with the texture applied:
    glEnable(GL.TEXTURE_2D);
    glBindTexture(GL.TEXTURE_2D, texture);
    
    % Setup proper alpha-blending:
    glBlendFunc(blend_sfactor,blend_dfactor);

    % Update our virtual position over the plane - our zoom factor:
    oz=zpos;
    zpos=zpos+stepper;
    ad=(zpos*100.0)-(oz*100.0);

    % 'stepper' defines direction and speed of zoom:
    if stepper>0
        while (zpos>7500.0) zpos=zpos-5000.0; end
    else
        while (zpos<2500.0) zpos=zpos+5000.0; end
    end

    % Rotate the plane while zooming in:
    rotpos=rotpos+rotspeed;
    if (rotpos>360.0) rotpos=rotpos-360.0; end
    if (rotpos<-360.0) rotpos=rotpos+360.0; end
    glRotatef(rotpos,0.0,1.0,0.0);

    % Update layer scalings and layer order:
    for i=0:4
        % Update layers virtual depth level:
        if (hold==0) layercnt(i+1)=layercnt(i+1) + 2; end
        layercnt(i+1)=layercnt(i+1)+ad;

        % Update depends on zoom direction 'stepper' ...
        if stepper>0
            % Has layer i reached the outer edge?
            if (layercnt(i+1)>=100000)
                % Yes. Reset its depths and make it the new first layer:
                layercnt(i+1)=0;
                firstlayer=i;
                % Randomize its texture coordinates:
                xo(i+1)=rand;
                yo(i+1)=rand;
            end
        else
            % Has layer i reached the inner edge?
            if (layercnt(i+1)<=0)
                % Yes. Reset its depths and make it the new first layer:
                layercnt(i+1)=100000;
                firstlayer=mod(i+1, 5);
                % Randomize its texture coordinates:
                xo(i+1)=rand;
                yo(i+1)=rand;
            end
        end
    end

    % Draw all five layers, starting with 'firstlayer' at the bottom:
    for j=firstlayer:firstlayer+4
        i=mod(j,5)+1;
        if (showlayer(i))
            % Draw layer with updated spatial scaling factor 'fstep':
            fstep=exp(layercnt(i)/15000.0);

            % fk is the alpha-factor premultiplied onto this layer. Its
            % used to smoothly fade in a layer after appearing and smoothly
            % fade it out before disappearing.
            
            % Default to 1 == fully visible.
            fk=1.0;

            % Use a linear ramp to fade it in or out during the first 20000
            % units or last 20000 units:
            %if (layercnt(i)<10000) fk=layercnt(i)/10000.0; end
            if (layercnt(i)<20000) fk=layercnt(i)/20000.0; end
            if (layercnt(i)>80000) fk=(99999-layercnt(i))/20000.0; end

            % Apply computed alpha value fk to layer:
            glColor4f(texcolor(i+3,1),texcolor(i+3, 2),texcolor(i+3, 3),texcolor(i+3, 4)*fk);

            % Draw it:
            glBegin(GL.QUADS);
            glTexCoord2f(0.0+xo(i), 20.0+yo(i)); glVertex3f(-40.0*fstep, 0.0, -40.0*fstep);     % top left
            glTexCoord2f(20.0+xo(i),20.0+yo(i)); glVertex3f(40.0*fstep, 0.0, -40.0*fstep);  % top right
            glTexCoord2f(20.0+xo(i), 0.0+yo(i)); glVertex3f(40.0*fstep, 0.0, 40.0*fstep);  % bottom right
            glTexCoord2f(0.0+xo(i), 0.0+yo(i)); glVertex3f(-40.0*fstep, 0.0, 40.0*fstep);  % bottom left
            glEnd;
        end
    end

    % Restore modelview matrix to pre-render state:
    glPopMatrix;
    
    % Finish our OpenGl rendering:
    Screen('EndOpenGL', win);

    % Draw help text if enabled:
    if drawhelp
        DrawFormattedText(win, helptext, 'center', 'center');
    end
    
    % Perform flip in sync with retrace.
    vbl = Screen('Flip', win, vbl + 1.5 * ifi);
    
    % Activate our own OpenGL rendering code again for drawing the next frame:
    Screen('BeginOpenGL', win);

end % Next iteration of animation loop...

% Finish our OpenGl rendering:
Screen('EndOpenGL', win);

% We're done! Close window and release all texture ressources:
Screen('CloseAll');

% Restore synctest setting to pre-demo state:
Screen('Preference', 'SkipSyncTests', oldlevel);

% Exit.
return;
