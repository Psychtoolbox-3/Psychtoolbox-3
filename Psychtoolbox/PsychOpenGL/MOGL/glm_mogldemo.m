function glm_mogldemo

% GLM_MOGLDEMO  Demonstrate use of MATLAB-OpenGL toolbox independent of Psychtoolbox.
%
% Instead of using Psychtoolbox as host environment, one can also use the built-in
% glm - toolbox. It provides basic functions to open windows, swap buffers, collect
% keyboard responses and such -- a minimalistic Psychtoolbox.
% This demo demonstrates how to use MOGL with glm.
%
% The GLM only works on MacOS-X whereas the Psychtoolbox approach works on
% MacOS-X, GNU/Linux and M$-Windows.

% 15-Dec-2005 -- created (RFM)
% 08-Apr-2006 -- Added comments (MK)

% check the path
if ~exist('glmOpen.m','file'),
    error('add the ''core'' and ''wrap'' directories to the MATLAB search path');
end

% open screen
glmOpen;

% get aspect ratio of screen
xy=glmGetScreenInfo;
ar=xy(2)/xy(1);

% turn on lighting
glEnable(GL.LIGHTING);
glEnable(GL.LIGHT0);
glLightModelfv(GL.LIGHT_MODEL_TWO_SIDE,GL.TRUE);

% set depth properties
glEnable(GL.DEPTH_TEST);
glEnable(GL.CULL_FACE);
glCullFace(GL.BACK);

% define material
glMaterialfv(GL.FRONT_AND_BACK,GL.AMBIENT, [ .33 .22 .03 1 ]);
glMaterialfv(GL.FRONT_AND_BACK,GL.DIFFUSE, [ .78 .57 .11 1 ]);
glMaterialfv(GL.FRONT_AND_BACK,GL.SHININESS,27.8);

% configure for texture
glEnable(GL.TEXTURE_2D);
texname=glGenTextures(6);
load mogldemo.mat face
for i=1:6,
    glBindTexture(GL.TEXTURE_2D,texname(i));
    f=max(min(128*(1+face{i}),255),0);
    tx=repmat(flipdim(f,1),[ 1 1 3 ]);
    tx=permute(flipdim(uint8(tx),1),[ 3 2 1 ]);
    glTexImage2D(GL.TEXTURE_2D,0,GL.RGB,256,256,0,GL.RGB,GL.UNSIGNED_BYTE,tx);
    glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_WRAP_S,GL.REPEAT);
    glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_WRAP_T,GL.REPEAT);
    glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_MAG_FILTER,GL.NEAREST);
    glTexParameterfv(GL.TEXTURE_2D,GL.TEXTURE_MIN_FILTER,GL.NEAREST);
    glTexEnvfv(GL.TEXTURE_ENV,GL.TEXTURE_ENV_MODE,GL.MODULATE);
end

% set projection matrix
glMatrixMode(GL.PROJECTION);
glLoadIdentity;
gluPerspective(45,1/ar,0.1,100);

% set view matrix
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;
gluLookAt(3,3,5,0,0,0,0,1,0);

% set lighting
glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);
glLightfv(GL.LIGHT0,GL.DIFFUSE, [ 1 1 1 1 ]);
glLightfv(GL.LIGHT0,GL.AMBIENT, [ .1 .1 .1 1 ]);

% initialize amount and direction of rotation
theta=0;
rotatev=[ 0 0 1 ];

% loop until key pressed
while 1,
    
    % calculate rotation angle
    theta=mod(theta+1,360);
    rotatev=rotatev+0.1*[ sin((pi/180)*theta) sin((pi/180)*2*theta) sin((pi/180)*theta/5) ];
    rotatev=rotatev/sqrt(sum(rotatev.^2));

    % draw a cube (should do this with vertex arrays)
    glPushMatrix;
    glRotated(theta,rotatev(1),rotatev(2),rotatev(3));
    glClear;
    cubeface([ 4 3 2 1 ],texname(1));
    cubeface([ 5 6 7 8 ],texname(2));
    cubeface([ 1 2 6 5 ],texname(3));
    cubeface([ 3 4 8 7 ],texname(4));
    cubeface([ 2 3 7 6 ],texname(5));
    cubeface([ 4 1 5 8 ],texname(6));
    glPopMatrix;

    % show it
    glmSwapBuffers;

    % check for quit signal (keypress or mouse button)
    [xy,b]=glmGetMouse;
    k=max(glmGetKeys);
    if b || k,
        break
    end

end

% shut down
glDeleteTextures(length(texname),texname);
glmClose;

return

% draw one face of a textured cube
function cubeface( i, tx )
v=[ 0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 1 1 1 ; 0 1 1 ]'-0.5;
n=cross(v(:,i(2))-v(:,i(1)),v(:,i(3))-v(:,i(2)));
global GL
glBindTexture(GL.TEXTURE_2D,tx);
glBegin(GL.POLYGON);
glNormal3dv(n);
glTexCoord3dv([ 0 0 ]);
glVertex3dv(v(:,i(1)));
glTexCoord3dv([ 1 0 ]);
glVertex3dv(v(:,i(2)));
glTexCoord3dv([ 1 1 ]);
glVertex3dv(v(:,i(3)));
glTexCoord3dv([ 0 1 ]);
glVertex3dv(v(:,i(4)));
glEnd;
return
