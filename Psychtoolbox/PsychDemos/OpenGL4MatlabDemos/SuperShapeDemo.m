function SuperShapeDemo
% Draw SuperShapes with OpenGL
%
% This demo computes and displays a 3D supershape using the formula
% found by Gielis (A generic geometric transformation that unifies 
% a wide range of natural and abstract shapes. American Journal of
% Botany 90(3):333�338, 2003)
%
% See also http://en.wikipedia.org/wiki/Superformula
%
% "The parameter m relates to the number of symmetries of a shape, 
% parameter n1 determines flatness and sharpness of corners and convexity 
% of sides. Both parameters n2 and n3 denote whether the shape is inscribed 
% or circumscribed in the unit circle." (from
% http://know-center.tugraz.at/download_extern/papers/ISGI_SuperShapes_Lex_Kienreich.pdf)
%
% You can rotate the shape with the cursor keys.
%
% You can change shape parameters by pressing the keys shown in the demo and
% entering new numeric parameters.
%
% Pressing ESCape exits the demo.
%
% Please be patient when running this demo. Computation of the supershape and
% its normal vectors and mesh is compute intense and therefore time consuming.
% The only way to speed this up a lot would be to perform most of the computations
% inside of GLSL vertex- and geometry-shaders on the GPU. This advanced exercise
% for now left to the reader.
%
% This demo was written and contributed by Joana Leitao.
%

% History:
% 16-Sep-2013   jsl  Written.

if IsARM
    % Does not work on OpenGL-ES due to unsupported OpenGL display lists:
    error('Sorry, this demos does not work yet on embedded ARM devices.');
end

%% Initialization

% Setup unified keynames and normalized 0-1 color space:
PsychDefaultSetup(2);

% Suppress character output to console:
ListenChar(2);

% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL wrapper:
InitializeMatlabOpenGL;

% Find the screen to use for display:
screenid = max(Screen('Screens'));

% Define background color:
backgroundColor = GrayIndex(screenid);

% Open a double-buffered full-screen window:
[win , winRect] = PsychImaging('OpenWindow', screenid, backgroundColor);

DrawFormattedText(win, 'Hang on, I am working on my first supershape :)', 'center', 'center');
Screen('Flip', win);

% Setup the OpenGL rendering context of the onscreen window for use by
% OpenGL wrapper. After this command, all following OpenGL commands will
% draw into the onscreen window 'win':
Screen('BeginOpenGL', win);

% Get the aspect ratio of the screen:
ar = RectHeight(winRect) / RectWidth(winRect);

%% Set camera position etc

% Set viewport properly:
glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));

% Set projection matrix: This defines a perspective projection,
% corresponding to the model of a pin-hole camera - which is a good
% approximation of the human eye and of standard real world cameras --
% well, the best aproximation one can do with 3 lines of code ;-)
glMatrixMode(GL.PROJECTION);
glLoadIdentity;

% Field of view is 25 degrees from line of sight. Objects closer than
% 0.1 distance units or farther away than 100 distance units get clipped
% away, aspect ratio is adapted to the monitors aspect ratio:
gluPerspective(25, 1/ar, 0.1, 100);

% Setup modelview matrix: This defines the position, orientation and
% looking direction of the virtual camera:
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;

% Cam is located at 3D position (0,10,15), points upright (0,1,0) and fixates
% at the origin (0,0,0) of the worlds coordinate system:
gluLookAt(0,10,15,0,0,0,0,1,0);


%% Setup position and emission properties of the light source:

% Turn on OpenGL local lighting model: The lighting model supported by
% OpenGL is a local Phong model with Gouraud shading.
glEnable(GL.LIGHTING);

% Enable the first local light source GL.LIGHT_0. Each OpenGL
% implementation is guaranteed to support at least 8 light sources.
glEnable(GL.LIGHT0);

% Our point lightsource is at position (x,y,z) == (1,2,3)...
glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);

% % Emits white (1,1,1,1) diffuse light:
glLightfv(GL.LIGHT0,GL.DIFFUSE, [ 1 1 1 1 ]);

% There's also some white, but weak (R,G,B) = (0.1, 0.1, 0.1)
% ambient light present:
glLightfv(GL.LIGHT0,GL.AMBIENT, [0.1 0.1 0.1 1]);

%% Define Material Properties

% Define the cubes light reflection properties by setting up reflection
% coefficients for ambient, diffuse and specular reflection:
glMaterialfv(GL.FRONT,GL.AMBIENT, [ 0.0 0.0 1.0 1 ]);
glMaterialfv(GL.FRONT,GL.DIFFUSE, [ 0.1 0.1 1 1 ]);
glMaterialfv(GL.FRONT,GL.SHININESS,27.8);

%%

% Enable proper occlusion handling via depth tests:
glEnable(GL.DEPTH_TEST);

% Enable smooth gouraud shading:
glShadeModel(GL.SMOOTH);

% Set background clear color to 'gray' (R,G,B,A)=(0.5,0.5,0.5,1):
glClearColor(0.5,0.5,0.5,1);

% Clear out the backbuffer: This also cleans the depth-buffer for
% proper occlusion handling: You need to glClear the depth buffer whenever
% you redraw your scene, e.g., in an animation loop. Otherwise occlusion
% handling will screw up in funny ways...
glClear;

%% generate initial vertices
m = 8;
a = 1;
b = 1;
n1 = 0.5;
n2 = 0.5;
n3 = 0.5;

super_vertices = superformula(a,b,m,n1,n2,n3);

% predraw strip
super_list = glGenLists(1);
glNewList(super_list, GL.COMPILE);
supershape(super_vertices);
glEndList;

%% Animation loop: Run until escape key is pressed
while (1)
    [pressed, ~, keyCode] = KbCheck;
    
    if pressed
        
        if length(find(keyCode)) == 1 
            
            switch find(keyCode)
                
                case {KbName('LeftArrow')}
                    glRotatef(-5,0.0,1.0,0.0);
                    
                case {KbName('RightArrow')}
                    glRotatef(5,0.0,1.0,0.0);
                    
                case {KbName('UpArrow')}
                    glRotatef(-5,1.0,0.0,0.0);
                    
                case {KbName('DownArrow')}
                    glRotatef(5,1.0,0.0,0.0);
                    
                case {KbName('m'),KbName('a'),KbName('b'),KbName('j'),KbName('k'),KbName('l')}                    
                    % Finish OpenGL rendering into PTB window and check for OpenGL errors.
                    Screen('EndOpenGL', win);
                    
                    if find(keyCode) == KbName('m')
                        m = GetEchoNumber(win, 'Please input a new number for parameter m followed by enter', 0.05*winRect(3),0.05*winRect(4), [], backgroundColor);
                    end
                    
                    if find(keyCode) == KbName('a')
                        a = GetEchoNumber(win, 'Please input a new number for parameter a followed by enter', 0.05*winRect(3),0.05*winRect(4), [], backgroundColor);
                    end
                    
                    if find(keyCode) == KbName('b')
                        b = GetEchoNumber(win, 'Please input a new number for parameter b followed by enter', 0.05*winRect(3),0.05*winRect(4), [], backgroundColor);
                    end
                    
                    if find(keyCode) == KbName('j')
                        n1 = GetEchoNumber(win, 'Please input a new number for parameter n1 followed by enter', 0.05*winRect(3),0.05*winRect(4), [], backgroundColor);
                    end
                    
                    if find(keyCode) == KbName('k')
                        n2 = GetEchoNumber(win, 'Please input a new number for parameter n2 followed by enter', 0.05*winRect(3),0.05*winRect(4), [], backgroundColor);
                    end
                    
                    if find(keyCode) == KbName('l')
                        n3 = GetEchoNumber(win, 'Please input a new number for parameter n3 followed by enter', 0.05*winRect(3),0.05*winRect(4), [], backgroundColor);
                    end

                    Screen('Flip', win);
                    DrawFormattedText(win, 'Hang on, I am thinking hard about your numbers...', 'center', 'center');
                    Screen('Flip', win);
                    
                    % Switch to OpenGL rendering again for drawing of next frame:
                    Screen('BeginOpenGL', win);                    

                    % generate vertices
                    super_vertices = superformula(a,b,m,n1,n2,n3);
                    % predraw strip
                    super_list = glGenLists(1);
                    glNewList(super_list, GL.COMPILE);
                    supershape(super_vertices);
                    glEndList;
               
                case {KbName('ESCAPE')}
                    break;
            end
        end
    end
    
    % Clear framebuffer and redraw supershape:
    glClear;
    glCallList(super_list);
    
    % Finish OpenGL rendering into PTB window and check for OpenGL errors.
    Screen('EndOpenGL', win);
    
    text1 = 'SuperShape with parameters a, b, m, j (aka n1), k (aka n2) and l (aka n3).';
    text2 = ['Current values: m = ' num2str(m) ', a = ' num2str(a) ', b = ' num2str(b) ', j = ' num2str(n1) ', k = ' num2str(n2) ', l = ' num2str(n3)];
    text3 = 'To change a parameter press the corresponding key. To rotate the object use the arrows. To quit the program press escape';
    DrawFormattedText(win, text1, 0.05*winRect(3), 0.8*winRect(4));
    DrawFormattedText(win, text2, 0.05*winRect(3), 0.83*winRect(4));
    DrawFormattedText(win, text3, 0.05*winRect(3), 0.88*winRect(4));
    Screen('FrameRect', win, 0.5, [0 0 5 5]);
    
    % Show rendered image at next vertical retrace:
    Screen('Flip', win);
    
    % Switch to OpenGL rendering again for drawing of next frame:
    Screen('BeginOpenGL', win);
end

% Delete super_list and exit OpenGL rendering:
glDeleteLists(super_list,1);
Screen('EndOpenGL', win);

% Close onscreen window and release all other ressources:
Screen('CloseAll');

% Enable character output to Matlab:
ListenChar(0);

% Well done!
return


%% this function creates quad_strips with the given input vertices
function supershape(vertices)
global GL

for row_ind = 1:size(vertices(:,:,1),1)-1
    
    % Begin drawing of a new quad strip:
    glBegin(GL.QUAD_STRIP);
    
    for column_ind = 1:size(vertices(:,:,1),2)
        
        for k = 1:-1:0
            
            aux_var = row_ind + k;
            
            x = vertices(aux_var,column_ind,1);
            y = vertices(aux_var,column_ind,2);
            z = vertices(aux_var,column_ind,3);
            
            norm_vector = calculateNormVec(vertices, aux_var, column_ind);
            glNormal3f(norm_vector(1),norm_vector(2),norm_vector(3));
            glVertex3f(x,y,z)
        end
    end
    
    glEnd;
end

return


%% this function creates the vertices of a 3D supershape using the formula
% found by Gielis (A generic geometric transformation that unifies 
% a wide range of natural and abstract shapes. American Journal of
% Botany 90(3):333�338, 2003). See also http://en.wikipedia.org/wiki/Superformula

% "The parameter m relates to the number of symmetries of a shape, 
% parameter n1 determines flatness and sharpness of corners and convexity 
% of sides. Both parameters n2 and n3 denote whether the shape is inscribed 
% or circumscribed in the unit circle." (from
% http://know-center.tugraz.at/download_extern/papers/ISGI_SuperShapes_Lex_Kienreich.pdf)
function [super_vertices] = superformula(a,b,m,n1,n2,n3)

increments = 100;
theta = -pi:2*pi/(increments-1):pi;
phi = -pi/2:2*pi/(increments-1):pi/2;

x = zeros(length(theta),length(phi));
y = zeros(length(theta),length(phi));
z = zeros(length(theta),length(phi));
super_vertices = zeros(length(theta),length(phi),3);

for theta_ind = 1:length(theta)
    
    for phi_ind = 1:length(phi)
        
        r_theta = (abs(1/a*cos(theta(theta_ind)* m/4)).^n2  +  abs(1/b*sin(theta(theta_ind)*m/4)).^n3)^(1/n1);
        r_phi = (abs(1/a*cos(phi(phi_ind)* m/4)).^n2  +  abs(1/b*sin(phi(phi_ind)*m/4)).^n3)^(1/n1);
        
        x(theta_ind,phi_ind) = r_theta*cos(theta(theta_ind))*r_phi*cos(phi(phi_ind));
        y(theta_ind,phi_ind) = r_theta*sin(theta(theta_ind))*r_phi*cos(phi(phi_ind));
        z(theta_ind,phi_ind) = r_phi*sin(phi(phi_ind));
    end;
    
end

super_vertices(:,:,1) = x./max(max(abs(x)));
super_vertices(:,:,2) = y./max(max(abs(y)));
super_vertices(:,:,3) = z./max(max(abs(z)));

return


%% this function creates the normal vector at a given vertex. 
function [norm_vector] = calculateNormVec(vertices, row_ind, column_ind)

size_aux = vertices(:,:,1);

if column_ind == 1
    
    % if column_ind equals one, this mean phi equals -pi/2 and hence, that we
    % are at the buttom point of the supershape, which is oriented along the
    % z-axis. Hence the normal vector at any of the vertices with column = 1
    % equals [0 0 -1];
    norm_vector = [0 0 -1];

elseif column_ind == size(size_aux,2)
    
    % if column_ind equals size(size_aux,2), this mean phi equals pi/2 and 
    % hence, that we are at the upper point of the supershape, which is 
    % oriented along the z-axis. Hence the normal vector at any of the  
    % vertices with column = size(size_aux,2) equals [0 0 1];    
    norm_vector = [0 0 1];

elseif row_ind == 1 || row_ind == size(size_aux,1)
    
    % in these cases we are at the edge where the supershape closes
    % (basically theta equals 0 and 2*pi. To calculate the normal at these
    % the same as below is done, only that one has to consider values at
    % the end/beginning of the rows instead.
    row_aux = 1;
    
    norm_vector1 = cross(vertices(row_aux+1,column_ind,:)-vertices(row_aux,column_ind,:),vertices(row_aux,column_ind+1,:)-vertices(row_aux,column_ind,:));
    norm_vector1_aux = [norm_vector1(1,1,1) norm_vector1(1,1,2) norm_vector1(1,1,3)]./norm([norm_vector1(1,1,1) norm_vector1(1,1,2) norm_vector1(1,1,3)]);
    
    norm_vector2 = cross(vertices(row_aux,column_ind+1,:)-vertices(row_aux,column_ind,:),vertices(end-1,column_ind,:)-vertices(row_aux,column_ind,:));
    norm_vector2_aux = [norm_vector2(1,1,1) norm_vector2(1,1,2) norm_vector2(1,1,3)]./norm([norm_vector2(1,1,1) norm_vector2(1,1,2) norm_vector2(1,1,3)]);
    
    norm_vector3 = cross(vertices(end-1,column_ind,:)-vertices(row_aux,column_ind,:),vertices(row_aux,column_ind-1,:)-vertices(row_aux,column_ind,:));
    norm_vector3_aux = [norm_vector3(1,1,1) norm_vector3(1,1,2) norm_vector3(1,1,3)]./norm([norm_vector3(1,1,1) norm_vector3(1,1,2) norm_vector3(1,1,3)]);
    
    norm_vector4 = cross(vertices(row_aux,column_ind-1,:)-vertices(row_aux,column_ind,:),vertices(row_aux+1,column_ind,:)-vertices(row_aux,column_ind,:));
    norm_vector4_aux = [norm_vector4(1,1,1) norm_vector4(1,1,2) norm_vector4(1,1,3)]./norm([norm_vector4(1,1,1) norm_vector4(1,1,2) norm_vector4(1,1,3)]);
    
    norm_vector = sum([norm_vector1_aux; norm_vector2_aux; norm_vector3_aux; norm_vector4_aux]);

else
    % To calculate the normal at all the other vertices the normals at all
    % the four surrounding quad faces are calculated and then added up to
    % give the normal at the vertex
    norm_vector1 = cross(vertices(row_ind+1,column_ind,:)-vertices(row_ind,column_ind,:),vertices(row_ind,column_ind+1,:)-vertices(row_ind,column_ind,:));
    norm_vector1_aux = [norm_vector1(1,1,1) norm_vector1(1,1,2) norm_vector1(1,1,3)]./norm([norm_vector1(1,1,1) norm_vector1(1,1,2) norm_vector1(1,1,3)]);
    
    norm_vector2 = cross(vertices(row_ind,column_ind+1,:)-vertices(row_ind,column_ind,:),vertices(row_ind-1,column_ind,:)-vertices(row_ind,column_ind,:));
    norm_vector2_aux = [norm_vector2(1,1,1) norm_vector2(1,1,2) norm_vector2(1,1,3)]./norm([norm_vector2(1,1,1) norm_vector2(1,1,2) norm_vector2(1,1,3)]);
    
    norm_vector3 = cross(vertices(row_ind-1,column_ind,:)-vertices(row_ind,column_ind,:),vertices(row_ind,column_ind-1,:)-vertices(row_ind,column_ind,:));
    norm_vector3_aux = [norm_vector3(1,1,1) norm_vector3(1,1,2) norm_vector3(1,1,3)]./norm([norm_vector3(1,1,1) norm_vector3(1,1,2) norm_vector3(1,1,3)]);
    
    norm_vector4 = cross(vertices(row_ind,column_ind-1,:)-vertices(row_ind,column_ind,:),vertices(row_ind+1,column_ind,:)-vertices(row_ind,column_ind,:));
    norm_vector4_aux = [norm_vector4(1,1,1) norm_vector4(1,1,2) norm_vector4(1,1,3)]./norm([norm_vector4(1,1,1) norm_vector4(1,1,2) norm_vector4(1,1,3)]);
    
    norm_vector = sum([norm_vector1_aux; norm_vector2_aux; norm_vector3_aux; norm_vector4_aux]);
    
end

return
