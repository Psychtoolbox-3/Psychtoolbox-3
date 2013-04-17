function SpinningMovieCube(moviename)
% SpinningMovieCubeDemo - Demonstrate use of MATLAB-OpenGL toolbox
%
% This demo demonstrates use of OpenGL commands in a Matlab script together
% with the movie playback functions of PTB. It shows a randomly
% spinning, three dimensional textured cube. The six sides of the cube
% show a quicktime video, loaded and played from a quicktime movie file.
%
% Stop the demo any time by pressing any key.
%
% Notable implementation details:
%
% The implementation is nearly identical as SpinningCubeDemo, so make sure
% you understand that file first.
%
% This demo uses the Screen('GetOpenGLTexture') function to make a
% Psychtoolbox texture  (loaded from the movie) available to
% OpenGL for 3D texture mapping.
%
% It opens a movie file and then - in a loop - fetches video images frame
% by frame. The images which are stored as Psychtoolbox textures are then
% made available as standard OpenGL textures for drawing onto the sides of
% the spinning cube.
%

% 15-Dec-2005 -- created (RFM)
% 21-Jan-2006 -- Modified for use with OpenGL-Psychtoolbox (MK)
% 16-Feb-2006 -- Modified for use with new MOGL (MK)
% 05-Mar-2006 -- Cleaned up for public consumption (MK)
% 23-Aug-2012 -- Adapt to PTB 3.0.10 file structure, cleanup. (MK)
% 06-Apr-2013 -- Make compatible with OpenGL-ES1.1. (MK)

% Assign default movie file, if none provided:
if nargin < 1
    moviename= [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];   
end;
moviename %#ok<NOPRT>

% Is the script running in OpenGL Psychtoolbox?
AssertOpenGL;

% Find the screen to use for display:
screenid=max(Screen('Screens'));

% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL for Matlab wrapper:
InitializeMatlabOpenGL(1);

% Open a double-buffered full-screen window on the main displays screen.
[win , winRect] = Screen('OpenWindow', screenid);

% Open movie file, get a handle to the movie and start playback:
movie=Screen('OpenMovie', win, moviename);
Screen('PlayMovie', movie, 1, 1, 1);

% Setup the OpenGL rendering context of the onscreen window for use by
% OpenGL wrapper. After this command, all following OpenGL commands will
% draw into the onscreen window 'win':
Screen('BeginOpenGL', win);

% Get the aspect ratio of the screen:
ar=winRect(4)/winRect(3);

% Turn on OpenGL local lighting model: The lighting model supported by
% OpenGL is a local Phong model with Gouraud shading.
glEnable(GL_LIGHTING);

% Enable the first local light source GL_LIGHT_0. Each OpenGL
% implementation is guaranteed to support at least 8 light sources. 
glEnable(GL_LIGHT0);

% Enable two-sided lighting - Back sides of polygons are lit as well.
glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);

% Enable proper occlusion handling via depth tests:
glEnable(GL_DEPTH_TEST);

% Define the cubes light reflection properties by setting up reflection
% coefficients for ambient, diffuse and specular reflection:
glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT, [ .33 .22 .03 1 ]);
glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE, [ .78 .57 .11 1 ]);
glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,27.8);

% Set projection matrix: This defines a perspective projection,
% corresponding to the model of a pin-hole camera - which is a good
% approximation of the human eye and of standard real world cameras --
% well, the best aproximation one can do with 3 lines of code ;-)
glMatrixMode(GL_PROJECTION);
glLoadIdentity;
% Field of view is 25 degrees from line of sight. Objects close than
% 0.1 distance units or farther away than 100 distance units get clipped
% away, aspect ratio is adapted to the monitors aspect ratio:
gluPerspective(25,1/ar,0.1,100);

% Setup modelview matrix: This defines the position, orientation and
% looking direction of the virtual camera:
glMatrixMode(GL_MODELVIEW);
glLoadIdentity;

% Cam is located at 3D position (3,3,5), points upright (0,1,0) and fixates
% at the origin (0,0,0) of the worlds coordinate system:
gluLookAt(3,3,5,0,0,0,0,1,0);

% Setup position and emission properties of the light source:

% Set background color to 'black':
glClearColor(0,0,0,0);

% Point lightsource at (1,2,3)...
glLightfv(GL_LIGHT0,GL_POSITION,[ 1 2 3 0 ]);
% Emits white (1,1,1,1) diffuse light:
glLightfv(GL_LIGHT0,GL_DIFFUSE, [ 1 1 1 1 ]);

% There's also some white, but weak (R,G,B) = (0.1, 0.1, 0.1)
% ambient light present:
glLightfv(GL_LIGHT0,GL_AMBIENT, [ .1 .1 .1 1 ]);

% initialize amount and direction of rotation
theta=0;
rotatev=[ 0 0 1 ];

% Setup global variables needed for texture mapping:
global target;
global tu
global tv

% Run playback loop until key pressed:
while (1)
    % Fetch next video frame from movie file and return a Psychtoolbox
    % texture handle to it:
    Screen('EndOpenGL', win);
    texid = Screen('GetMovieImage', win, movie);
    Screen('BeginOpenGL', win);
    
    % Valid texture?
    if texid<=0
        % No. Abort.
        break;
    end;
    
    % Retrieve and assign textures for each cube face... Here we assign six
    % times the same image...
    for i=1:6
        % Retrieve size of movie texture image:
        [ imw imh] = Screen('WindowSize', texid);
        % Retrieve an OpenGL texture handle and texture mapping parameters:
        % texname(i) contains the OpenGL texture id, target is the texture
        % type, tu and tv are the texture coodinates of the (imw,imh)
        % position of the texture:
        [ texname(i) target tu tv] = Screen('GetOpenGLTexture', win, texid, imw, imh); %#ok<*AGROW>
    end;
    
    % calculate rotation angle and axis of cube for this frame:
    theta=mod(theta+1,360);
    rotatev=rotatev+0.1*[ sin((pi/180)*theta) sin((pi/180)*2*theta) sin((pi/180)*theta/5) ];
    rotatev=rotatev/sqrt(sum(rotatev.^2));

    % Draw cube:
    glPushMatrix;
    glRotatef(theta,rotatev(1),rotatev(2),rotatev(3));
    glClear;
        
    cubeface([ 4 3 2 1 ],texname(1));
    cubeface([ 5 6 7 8 ],texname(2));
    cubeface([ 1 2 6 5 ],texname(3));
    cubeface([ 3 4 8 7 ],texname(4));
    cubeface([ 2 3 7 6 ],texname(5));
    cubeface([ 4 1 5 8 ],texname(6));

    glPopMatrix;
    
    % End OpenGL rendering and check for OpenGL errors.
    Screen('EndOpenGL', win);

    % Show rendered image on next retrace:
    Screen('Flip', win);
    
    % Close the texture after drawing so we don't use up system memory:
    Screen('Close', texid);

    % Switch to OpenGL rendering again for drawing of next frame:
    Screen('BeginOpenGL', win);
    
    % Abort on keypress:
    if KbCheck
        break;
    end;
end

% End of the show:

% Disable OpenGL drawing:
Screen('EndOpenGL', win);

% Stop movie playback and close movie:
Screen('CloseMovie', movie);

% Close onscreen window and release all ressources:
Screen('CloseAll');

return

% draw one face of a textured cube
function cubeface( i, tx )
v=[ 0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 1 1 1 ; 0 1 1 ]'-0.5;
n=cross(v(:,i(2))-v(:,i(1)),v(:,i(3))-v(:,i(2)));
global GL
global target
global tu
global tv

% Enable and bind proper texture:
glEnable(target);
glBindTexture(target,tx);
glTexEnvfv(GL.TEXTURE_ENV,GL.TEXTURE_ENV_MODE,GL.MODULATE);
% Setup filtering for the textures:
glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_MAG_FILTER,GL.NEAREST);
glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_MIN_FILTER,GL.NEAREST);

% Begin drawing of a new quad:
glBegin(GL.QUADS);

% Assign n as normal vector for this polygons surface normal:
glNormal3f(n(1), n(2), n(3));

% Define vertex 1 by assigning a texture coordinate and a 3D position:
glTexCoord2f(0, 0);
glVertex3f(v(1,i(1)),v(2,i(1)),v(3,i(1)));
% Define vertex 2 by assigning a texture coordinate and a 3D position:
glTexCoord2f(tu, 0);
glVertex3f(v(1,i(2)),v(2,i(2)),v(3,i(2)));
% Define vertex 3 by assigning a texture coordinate and a 3D position:
glTexCoord2f(tu, tv);
glVertex3f(v(1,i(3)),v(2,i(3)),v(3,i(3)));
% Define vertex 4 by assigning a texture coordinate and a 3D position:
glTexCoord2f(0, tv);
glVertex3f(v(1,i(4)),v(2,i(4)),v(3,i(4)));
glEnd;

return
