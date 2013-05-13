function MorphTextureDemo(dotson, normalson, stereomode, usefastoffscreenwindows)
% function MorphTextureDemo([dotson][, normalson][,stereomode][, usefastoffscreenwindows])
%
% MorphTextureDemo -- Demonstrates use of "moglmorpher" for fast morphing
% and rendering of 3D shapes *and* textures. See "help moglmorpher" for
% info on moglmorphers purpose and capabilities.
%
% This demo will load two morpheable shapes from OBJ files and two
% textures, then morph them continously into each other, using a simple
% sine-function to define the timecourse of the morph.
%
% Control keys and their meaning: 'a' == Zoom out by moving object away
% from viewer. 'z' == Zoom in by moving object close to viewer. 'k' and 'l'
% == Rotate object around axis. 'ESC' == Quit demo.
%
% Options:
%
% dotson = If set to 0 (default), just show surface. If set to 1, some dots
% are plotted to visualize the vertices of the underlying mesh. If set to
% 2, the mesh itself is superimposed onto the shape. If set to 3 or 4, then
% the projected vertex 2D coordinates are also visualized in a standard
% Matlab figure window.
%
% normalson = If set to 1, then the surface normal vectors will get
% visualized as small green lines on the surface.
%
% stereomode = n. For n>0 this activates stereoscopic rendering - The shape
% is rendered from two slightly different viewpoints and one of
% Psychtoolbox's built-in stereo display algorithms is used to present the
% 3D stimulus. This is very preliminary so it doesn't work that well yet.
%
% usefastoffscreenwindows = If set to 0 (default), work on any graphics
% card. If you have recent hardware, set it to 1. That will enable support
% for fast offscreen windows - and a much faster implementation of shape
% morphing.
%
% This demo and the morpheable OBJ shapes were contributed by Dr. Quoc C.
% Vuong, MPI for Biological Cybernetics, Tuebingen, Germany.

morphnormals = 1;
global win;

% Is the script running in OpenGL Psychtoolbox?
AssertOpenGL;

% Should per pixel lighting via OpenGL shading language be used? This doesnt work
% well yet.
perpixellighting = 0;

% Some default settings for rendering flags:
if nargin < 1 || isempty(dotson)
    dotson = 0;     % turn reference dots: off(0), on (1) or show reference lines (2)
end
dotson %#ok<NOPRT>

if nargin < 2 || isempty(normalson)
    normalson = 0;     % turn reference dots: off(0), on (1) or show reference lines (2)
end
normalson %#ok<NOPRT>

if nargin < 3 || isempty(stereomode)
    stereomode = 0;
end;
stereomode %#ok<NOPRT>

if nargin < 4 || isempty(usefastoffscreenwindows)
    usefastoffscreenwindows = 0;
end
usefastoffscreenwindows %#ok<NOPRT>

% Response keys: Mapping of keycodes to keynames.
KbName('UnifyKeyNames');
closer = KbName('a');
farther = KbName('z');
quitkey = KbName('ESCAPE');
rotateleft = KbName('l');
rotateright = KbName('k');

% Load OBJs. This will define topology and use of texcoords and normals:
% One can call LoadOBJFile() multiple times for loading multiple objects.
basepath = [fileparts(which(mfilename)) '/'];
objs = [ LoadOBJFile([basepath 'texblob01.obj']) LoadOBJFile([basepath 'texblob02.obj']) ];

% Find the screen to use for display:
screenid=max(Screen('Screens'));

% Disable Synctests for this simple demo:
oldskip = Screen('Preference','SkipSyncTests', 1);

% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL for Matlab wrapper. We need to do this before the first call
% to any OpenGL function. Set debug level to zero -> Faster, but no error checking:
InitializeMatlabOpenGL(0,0);

% Open a double-buffered full-screen window: Everything is left at default
% settings, except stereomode:
if dotson~=3 && dotson~=4
   rect = [];
else
   rect = [0 0 500 500];
end;

if usefastoffscreenwindows
    [win , winRect] = Screen('OpenWindow', screenid, 0, rect, [], [], stereomode, [], kPsychNeedFastOffscreenWindows);
else
    [win , winRect] = Screen('OpenWindow', screenid, 0, rect, [], [], stereomode);
end

% Setup texture mapping:

% Load and create textures in Psychtoolbox:
basepath = [ PsychtoolboxRoot 'PsychDemos' filesep ];
texname = [basepath 'konijntjes1024x768blur.jpg'];
texture = imread(texname);
texid(1) = Screen('MakeTexture', win, texture);

texname = [basepath 'konijntjes1024x768.jpg'];
texture = imread(texname);
texid(2) = Screen('MakeTexture', win, texture);

% Swap (u,v) <-> (v,u) to account for the transposed images read via Matlab imread():
texcoords(2,:) = objs{1}.texcoords(1,:);
texcoords(1,:) = 1 - objs{1}.texcoords(2,:);

gltextarget = GL.TEXTURE_RECTANGLE_EXT;

% Rectangle texture: We need to rescale our texcoords as they are made for
% power-of-two textures, not rectangle textures:
texcoords(1,:) = texcoords(1,:) * size(texture,1);
texcoords(2,:) = texcoords(2,:) * size(texture,2);

% Reset moglmorpher:
moglmorpher('reset');

% Add the OBJS to moglmorpher for use as morph-shapes:
for i=1:size(objs,2)
    objs{i}.texcoords = texcoords; % Add modified texture coords.
    meshid(i) = moglmorpher('addMesh', objs{i}); %#ok<AGROW,NASGU>
end

% Output count of morph shapes:
count = moglmorpher('getMeshCount') %#ok<NOPRT,NASGU>

% Setup the OpenGL rendering context of the onscreen window for use by
% OpenGL wrapper. After this command, all following OpenGL commands will
% draw into the onscreen window 'win':
Screen('BeginOpenGL', win);

if perpixellighting==1
    % Load a GLSL shader for per-pixel lighting, built a GLSL program out of it...
    shaderpath = [PsychtoolboxRoot '/PsychDemos/OpenGL4MatlabDemos/GLSLDemoShaders/'];
    glsl=LoadGLSLProgramFromFiles([shaderpath 'Pointlightshader'],1);
    % ...and activate the shader program:
    glUseProgram(glsl);
end;

% Setup texture mapping for our morphed texture:
glEnable(gltextarget);

% Choose texture application function: It shall modulate the light
% reflection properties of the the objects surface:
glTexEnvfv(GL.TEXTURE_ENV,GL.TEXTURE_ENV_MODE,GL.MODULATE);

% Get the aspect ratio of the screen, we need to correct for non-square
% pixels if we want undistorted displays of 3D objects:
ar=winRect(4)/winRect(3);

% Turn on OpenGL local lighting model: The lighting model supported by
% OpenGL is a local Phong model with Gouraud shading.
glEnable(GL.LIGHTING);

% Enable the first local light source GL.LIGHT_0. Each OpenGL
% implementation is guaranteed to support at least 8 light sources. 
glEnable(GL.LIGHT0);

% Enable proper occlusion handling via depth tests:
glEnable(GL.DEPTH_TEST);

% Define the light reflection properties by setting up reflection
% coefficients for ambient, diffuse and specular reflection:
glMaterialfv(GL.FRONT_AND_BACK,GL.AMBIENT, [ 0.5 0.5 0.5 1 ]);
glMaterialfv(GL.FRONT_AND_BACK,GL.DIFFUSE, [ .7 .7 .7 1 ]);
glMaterialfv(GL.FRONT_AND_BACK,GL.SPECULAR, [ 0.2 0.2 0.2 1 ]);
glMaterialfv(GL.FRONT_AND_BACK,GL.SHININESS,12);

% Make sure that surface normals are always normalized to unit-length,
% regardless what happens to them during morphing. This is important for
% correct lighting calculations:
glEnable(GL.NORMALIZE);

% Set projection matrix: This defines a perspective projection,
% corresponding to the model of a pin-hole camera - which is a good
% approximation of the human eye and of standard real world cameras --
% well, the best aproximation one can do with 3 lines of code ;-)
glMatrixMode(GL.PROJECTION);
glLoadIdentity;

% Field of view is +/- 25 degrees from line of sight. Objects close than
% 0.1 distance units or farther away than 200 distance units get clipped
% away, aspect ratio is adapted to the monitors aspect ratio:
gluPerspective(25.0,1/ar,0.1,200.0);

% Setup modelview matrix: This defines the position, orientation and
% looking direction of the virtual camera:
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;

% Setup position of lightsource wrt. origin of world:
% Pointlightsource at (20 , 20, 20)...
glLightfv(GL.LIGHT0,GL.POSITION,[ 20 20 20 0 ]);

% Setup emission properties of the light source:

% Emits white (1,1,1,1) diffuse light:
glLightfv(GL.LIGHT0,GL.DIFFUSE, [ 1 1 1 1 ]);

% Emits white (1,1,1,1) specular light:
glLightfv(GL.LIGHT0,GL.SPECULAR, [ 1 1 1 1 ]);

% There's also some weak ambient light present:
glLightfv(GL.LIGHT0,GL.AMBIENT, [ 0.1 0.1 0.1 1 ]);

% Set size of points for drawing of reference dots
glPointSize(3.0);
glColor3f(0,0,1);

% Set thickness of reference lines:
glLineWidth(2.0);

% Add z-offset to reference lines, so they do not get occluded by surface:
glPolygonOffset(0, -5);
glEnable(GL.POLYGON_OFFSET_LINE);

% Initialize amount and direction of rotation for our slowly spinning,
% morphing objects:
theta=0;
rotatev=[ 0 0 1 ];

% Initialize morph vector:
w = [ 0 1 ];


% Setup initial z-distance of objects:
zz = 20.0;

ang = 0.0;      % Initial rotation angle

% Half eye separation in length units for quick & dirty stereoscopic
% rendering. Our way of stereo is not correct, but it makes for a
% nice demo. Figuring out proper values is not too difficult, but
% left as an exercise to the reader.
eye_halfdist=3;

% Finish OpenGL setup and check for OpenGL errors:
Screen('EndOpenGL', win);

% Compute initial morphed shape for next frame, based on initial weights:
moglmorpher('computeMorph', w, morphnormals);

% Compute initial morphed texture and get standard OpenGL texture handle and
% target from Psychtoolbox for use of result with MOGL:
[morphtexid, gltexid, gltextarget] = moglmorpher('morphTexture', win, w, texid);

% Retrieve duration of a single monitor flip interval: Needed for smooth
% animation.
ifi = Screen('GetFlipInterval', win);

% Initially sync us to the VBL:
vbl=Screen('Flip', win);

% Some stats...
tstart=vbl;
framecount = 0;
waitframes = 1;

% Animation loop: Run until key press or one minute has elapsed...
t = GetSecs;
while ((GetSecs - t) < 60)
    % Switch to OpenGL rendering for drawing of next frame:
    Screen('BeginOpenGL', win);

    % Bind morphed texture for rendering:
    glBindTexture(gltextarget, gltexid);
    
    % Left-eye cam is located at 3D position (-eye_halfdist,0,zz), points upright (0,1,0) and fixates
    % at the origin (0,0,0) of the worlds coordinate system:
    glLoadIdentity;
    gluLookAt(-eye_halfdist, 0, zz, 0, 0, 0, 0, 1, 0);

    % Draw into image buffer for left eye:
    Screen('EndOpenGL', win);
    Screen('SelectStereoDrawBuffer', win, 0);
    Screen('BeginOpenGL', win);

    % Clear out the depth-buffer for proper occlusion handling:
    glClear(GL.DEPTH_BUFFER_BIT);

    % Call our subfunction that does the actual drawing of the shape (see below):
    drawShape(ang, theta, rotatev, dotson, normalson);
    
    % Stereo rendering requested?
    if (stereomode > 0)
        % Yes! We need to render the same object again, just with a different
        % camera position, this time for the right eye:
                
        % Right-eye cam is located at 3D position (+eye_halfdist,0,zz), points upright (0,1,0) and fixates
        % at the origin (0,0,0) of the worlds coordinate system:
        glLoadIdentity;
        gluLookAt(+eye_halfdist, 0, zz, 0, 0, 0, 0, 1, 0);
        % Draw into image buffer for right eye:
        Screen('EndOpenGL', win);
        Screen('SelectStereoDrawBuffer', win, 1);
        Screen('BeginOpenGL', win);
    
        % Clear out the depth-buffer for proper occlusion handling:
        glClear(GL.DEPTH_BUFFER_BIT);
        
        % Call subfunction that does the actual drawing of the shape (see below):
        drawShape(ang, theta, rotatev, dotson, normalson)
    end;
    
    % Finish OpenGL rendering into Psychtoolbox - window and check for OpenGL errors.
    Screen('EndOpenGL', win);

    % For the fun of it: Draw morphed texture into a corner as an overlay:
    Screen('DrawTexture', win, morphtexid, [], ScaleRect(Screen('Rect', morphtexid), 0.33, 0.33))
    
    % Tell Psychtoolbox that drawing of this stim is finished, so it can optimize
    % drawing:
    Screen('DrawingFinished', win);

    % Now that all drawing commands are submitted, we can do the other stuff before
    % the Flip:
    
    % Calculate rotation angle of object for next frame:
    theta=mod(theta+0.1, 360);
    rotatev=rotatev+0.0001*[ sin((pi/180)*theta) sin((pi/180)*2*theta) sin((pi/180)*theta/5) ];
    rotatev=rotatev/sqrt(sum(rotatev.^2));
    
    % Compute simple morph weight vector for next frame:
    w(1)=(sin(framecount / 100 * 3.1415 * 2) + 1)/2;
    w(2)=1-w(1);

    % Compute morphed shape for next frame, based on new weight vector:
    moglmorpher('computeMorph', w, morphnormals);

    % Compute morphed texture and get standard OpenGL texture handle and
    % target from Psychtoolbox for use of result with MOGL:
    [morphtexid, gltexid, gltextarget] = moglmorpher('morphTexture', win, w, texid);
    
    if 0
        % Test morphed geometry readback:
        mverts = moglmorpher('getGeometry');
        scatter3(mverts(1,:), mverts(2,:), mverts(3,:));
        drawnow;
    end

    % Check for keyboard press:
    [KeyIsDown, endrt, KeyCode] = KbCheck;
    if KeyIsDown
        if ( KeyIsDown==1 && KeyCode(closer)==1 )
            zz=zz-0.1;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(farther)==1 )
            zz=zz+0.1;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(rotateright)==1 )
            ang=ang+1.0;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(rotateleft)==1 )
            ang=ang-1.0;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(quitkey)==1 )
            break;
        end
    end

    % Update frame animation counter:
    framecount = framecount + 1;

    % We're done for this frame:

    % Show rendered image 'waitframes' refreshes after the last time
    % the display was updated and in sync with vertical retrace:
    vbl = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi);
    %Screen('Flip', win, 0, 0, 2);
end

vbl = Screen('Flip', win);

% Calculate and display average framerate:
fps = framecount / (vbl - tstart) %#ok<NOPRT,NASGU>

% Reset moglmorpher:
moglmorpher('reset');

% Close onscreen window and release all other ressources:
%Screen('Flip', win);
Screen('CloseAll');

% Reenable Synctests after this simple demo:
Screen('Preference','SkipSyncTests', oldskip);

% Well done!
return

% drawShape does the actual drawing:
function drawShape(ang, theta, rotatev, dotson, normalson)
% GL needs to be defined as "global" in each subfunction that
% executes OpenGL commands:
global GL
global win

% Backup modelview matrix:
glPushMatrix;

% Setup rotation around axis:
glRotated(theta,rotatev(1),rotatev(2),rotatev(3));
glRotated(ang,0,1,0);

% Scale object by a factor of a:
a=0.1;
glScalef(a,a,a);

% Render current morphed shape via moglmorpher:
moglmorpher('render');

% Some extra visualizsation code for normals, mesh and vertices:
if (dotson == 1 || dotson == 3)
    % Draw some dot-markers at positions of vertices:
    % We disable lighting for this purpose:
    glDisable(GL.LIGHTING);
    % From all polygons, only their defining vertices are drawn:
    glPolygonMode(GL.FRONT_AND_BACK, GL.POINT);

    % Ask morpher to rerender the last shape:
    moglmorpher('render');

    % Reset settings for shape rendering:
    glPolygonMode(GL.FRONT_AND_BACK, GL.FILL);        
    glEnable(GL.LIGHTING);
end;

if (dotson == 2)
    % Draw connecting lines to visualize the underlying geometry:
    % We disable lighting for this purpose:
    glDisable(GL.LIGHTING);
    % From all polygons, only their connecting outlines are drawn:
    glPolygonMode(GL.FRONT_AND_BACK, GL.LINE);

    % Ask morpher to rerender the last shape:
    moglmorpher('render');

    % Reset settings for shape rendering:
    glPolygonMode(GL.FRONT_AND_BACK, GL.FILL);        
    glEnable(GL.LIGHTING);
end;

if (normalson > 0)
    % Draw surface normal vectors on top of object:
    glDisable(GL.LIGHTING);
    % Green is a nice color for this:
    glColor3f(0,1,0);

    % Ask morpher to render the normal vectors of last shape:
    moglmorpher('renderNormals', normalson);

    % Reset settings for shape rendering:
    glEnable(GL.LIGHTING);
    glColor3f(0,0,1);
 end;
 
if (dotson == 3 || dotson == 4)
   % Compute and retrieve projected screen-space vertex positions:
   vpos = moglmorpher('getVertexPositions', win);
   
   % Plot the projected 2D points into a Matlab figure window:
   vpos(:,2)=RectHeight(Screen('Rect', win)) - vpos(:,2);
	plot(vpos(:,1), vpos(:,2), '.');
   drawnow;
end;

% Restore modelview matrix:
glPopMatrix;

% Done, return to main-function:
return;
