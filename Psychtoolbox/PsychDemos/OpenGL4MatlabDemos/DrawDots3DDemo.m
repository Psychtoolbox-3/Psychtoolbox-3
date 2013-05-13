function DrawDots3DDemo(stereoMode, multiSample)
% DrawDots3DDemo -- Show fast drawing of 3D dots.
%
% Usage: DrawDots3DDemo([stereoMode=0][, multiSample=0]);
%
% This demo shows how to use the moglDrawDots3D() function to draw 3D dots
% in OpenGL 3D mode. The function is mostly equivalent to
% Screen('DrawDots') for drawing of 2D dots in regular 2D mode.
%
% The first subdemo simply fills the whole 3D scene with uniformly sampled
% random 3D dots, using the special sampling procedure
% CreateUniformDotsIn3DFrustum() which was contributed by Diederick
% Niehorster.
%
% The second subdemo shows how to use a GLSL vertex shader on modern GPU's
% to speed up complex drawing of complex 3D dot fields. It shows a nicely
% shaded, slowly rotating "Utah Teapot". Below the teapot is a primitive
% "sparkling fire" of 100 3D dots, which are lit by OpenGL and whose
% positions are computed in Matlab/Octave on the CPU. The teapot also emits
% a fountain of colorful particles from its nozzle. This fountain consists
% of 10000 particles, and the particle trajectories are computed on the GPU
% by use of a GLSL vertex shader. Please note that this subdemo may be
% pretty boring, not showing the magic fountain, if your GPU doesn't
% support shaders.
%
% The 3rd subdemo is a speed shootout: It draws the same fountain as demo
% 2, but as fast as it can, with sync of display updates to the vertical
% retrace disabled. The fountain is drawn 3 times a 20 seconds duration.
% First with purely Matlab computed trajectories, then with the same shader
% as in demo 2, then again with the shader, but additionally applying
% OpenGL display lists to store all data in the fast VRAM of the GPU to
% gain an additional speedup. At the end of these three benchmark runs, the
% demo will end and print out the average redraw rates attainable by the
% three methods. On a modern system with a modern graphics card, you should
% observe quite drastic speedups of GPU+VRAM vs. GPU vs. Matlab/CPU.
%
% The 4th subdemo simulates a "warp-flight" by creating a particle fountain
% that approaches the viewer, expanding while doing so.
%
%
% Btw. if you are a proud owner of a good 3D stereo setup, or at least of
% some anaglyph glasses, you should try the stereo display option as well,
% e.g., stereoMode == 8 for red-blue anaglyphs.
%
% Pressing the ESCape key continues the demo and progresses to next
% subdemo. Mouse clicks will pause some demos, until another mouse click
% continues the demo.
%
% Optional parameter:
%
% 'stereoMode' if set to a non-zero value, will render at lest the 2nd demo
% in a binocular stereo presentation mode, using the method specified in
% the 'stereoMode' flag. See for example ImagingStereoDemo for available
% modes.
%
% 'multiSample' if set to a non-zero value will enable multi-sample
% anti-aliasing. This however usually doesn't give good results with
% smoothed 3D dots.
%

% History:
% 03/01/2009  mk  Written.

% GL data structure needed for all OpenGL demos:
global GL;

if nargin < 2
    multiSample = 0;
end

if isempty(multiSample)
    multiSample = 0;
end

if nargin < 1
    stereoMode = [];
end

if isempty(stereoMode)
    stereoMode = 0;
end

if stereoMode
    stereoViews = 1;
else
    stereoViews = 0;
end

% Is the script running in OpenGL Psychtoolbox? Abort, if not.
AssertOpenGL;

% Restrict KbCheck to checking of ESCAPE key:
KbName('UnifyKeynames');
RestrictKeysForKbCheck(KbName('ESCAPE'));

% Find the screen to use for display:
screenid=max(Screen('Screens'));
if ismember(stereoMode, [4,5]) && IsWin
	screenid = 0;
end

try
    % Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
    % mogl OpenGL for Matlab wrapper:
    InitializeMatlabOpenGL;

    PsychImaging('PrepareConfiguration');

    % Open a double-buffered full-screen window on the main displays screen.
    [win, winRect] = PsychImaging('OpenWindow', screenid, 0, [], [], [], stereoMode, multiSample);

    if ismember(stereoMode, [6,7,8,9])
        SetAnaglyphStereoParameters('FullColorAnaglyphMode', win);
    end

    Screen('TextSize', win, 18);

    % Setup the OpenGL rendering context of the onscreen window for use by
    % OpenGL wrapper. After this command, all following OpenGL commands will
    % draw into the onscreen window 'win':
    Screen('BeginOpenGL', win);

    % Get the aspect ratio of the screen:
    ar=RectHeight(winRect) / RectWidth(winRect);

	% Set viewport properly:
	glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
	
    % Setup default drawing color to yellow (R,G,B)=(1,1,0). This color only
    % gets used when lighting is disabled - if you comment out the call to
    % glEnable(GL.LIGHTING).
    glColor3f(1,1,0);

    % Setup OpenGL local lighting model: The lighting model supported by
    % OpenGL is a local Phong model with Gouraud shading.

    % Enable the first local light source GL.LIGHT_0. Each OpenGL
    % implementation is guaranteed to support at least 8 light sources,
    % GL.LIGHT0, ..., GL.LIGHT7
    glEnable(GL.LIGHT0);

    % Enable alpha-blending for smooth dot drawing:
    glEnable(GL.BLEND);
    glBlendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA);

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

    % Our point lightsource is at position (x,y,z) == (1,2,3)...
    glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);

    % Cam is located at 3D position (3,3,5), points upright (0,1,0) and fixates
    % at the origin (0,0,0) of the worlds coordinate system:
    % The OpenGL coordinate system is a right-handed system as follows:
    % Default origin is in the center of the display.
    % Positive x-Axis points horizontally to the right.
    % Positive y-Axis points vertically upwards.
    % Positive z-Axis points to the observer, perpendicular to the display
    % screens surface.
    gluLookAt(0,0,10,0,0,0,0,1,0);

    % Set background clear color to 'black' (R,G,B,A)=(0,0,0,0):
    glClearColor(0,0,0,0);

    % Clear out the backbuffer: This also cleans the depth-buffer for
    % proper occlusion handling: You need to glClear the depth buffer whenever
    % you redraw your scene, e.g., in an animation loop. Otherwise occlusion
    % handling will screw up in funny ways...
    glClear;

    % Finish OpenGL rendering into PTB window. This will switch back to the
    % standard 2D drawing functions of Screen and will check for OpenGL errors.
    Screen('EndOpenGL', win);

    KbReleaseWait;
    DrawFormattedText(win, 'Now for a ugly demo of CPU based drawing of a uniform random dot field.\nPress ESCape key to continue and to finish a subdemo.', 'center', 'center', [255 255 0]);
    Screen('Flip', win);
    KbStrokeWait;

    % Show rendered image at next vertical retrace:
    Screen('Flip', win);

    ndots = 1000;

    % First version: Does not use occlusion testing via depth buffer, does not
    % use lighting. Uses auto-switching between 2D and 3D for simpler code:

    % 3D Dots animation loop: Runs until keypress:
    while ~KbCheck
        % Create random distribution of 3D dots inside our viewing frustum:
        [x,y,z] = CreateUniformDotsIn3DFrustum(ndots, 25, 1/ar, 0.1, 100);

        % Draw dots quickly: Common dotdiameter is 10 pixels, common color is
        % yellow. We move the center of the dots (aka position (0,0,0) to
        % position (0,0,10), so the above random transform applies properly:
        moglDrawDots3D(win, [x ; y; z], 10, [255 255 0 255], [0, 0, 10], 1, []);

        % Show'em:
        Screen('Flip', win, 0, 0);

        % A mouse button press will pause the animation:
        [x,y,buttons] = GetMouse;
        if any(buttons)
            % And wait for a single mouse click to continue:
            GetClicks;
        end
    end

    % Does this GPU support shaders?
    extensions = glGetString(GL.EXTENSIONS);
    if isempty(findstr(extensions, 'GL_ARB_shading_language')) || isempty(findstr(extensions, 'GL_ARB_shader_objects')) || isempty(findstr(extensions, 'GL_ARB_vertex_shader'))
        % Ok, no support for shading.
        shadingavail = 0;
    else
        % Use the shader stuff below this point...
        shadingavail = 1;
    end;

    KbReleaseWait;
    if shadingavail
        DrawFormattedText(win, 'Now for a beautiful demo of GPU based shading.\nPress ESCape key to continue and to finish a subdemo.', 'center', 'center', [255 255 0]);
    else
        DrawFormattedText(win, 'Now for another demo of CPU based drawing.\nPress ESCape key to continue and to finish a subdemo.\n\nUnfortunately your GPU does not support vertex shading\nso all following stuff will be pretty boring.', 'center', 'center', [255 255 0]);
    end

    Screen('Flip', win);
    KbStrokeWait;

    % Second version: Does use occlusion testing via depth buffer, does use
    % lighting. Uses manual switching between 2D and 3D for higher efficiency.
    % Creates a real 3D point-cloud around a teapot, as well as a vertex-shaded
    % fountain of particles that is emitted by the teapot:

    % Number of random dots, whose positions are computed in Matlab on CPU:
    ndots = 100;

    % Number of fountain particles whose positions are computed on the GPU:
    nparticles = 10000;

    % Diameter of particles in pixels:
    particleSize = 5;

    % 'StartPosition' is the 3D position where all particles originate. It is
    % faked to a position, so that the particles seem to originate from the
    % teapots "nozzle":
    StartPosition = [1.44, 0.40, 0];

    % Lifetime for each simulated particle, is chosen so that there seems to be
    % an infinite stream of particles, although the same particles are recycled
    % over and over:
    particlelifetime = 2;

    % Amount of "flow": A value of 1 will create a continuous stream, whereas
    % smaller value create bursts of particles:
    flowfactor = 1;

    if shadingavail
        % Load and setup the vertex shader for particle fountain animation:
        shaderpath = [PsychtoolboxRoot 'PsychDemos/OpenGL4MatlabDemos/GLSLDemoShaders/ParticleSimple'];
        glsl = LoadGLSLProgramFromFiles(shaderpath,1);

        % Bind shader so it can be setup:
        glUseProgram(glsl);

        % Assign static 3D startposition for fountain:
        glUniform3f(glGetUniformLocation(glsl, 'StartPosition'), StartPosition(1), StartPosition(2), StartPosition(3));

        % Assign lifetime:
        glUniform1f(glGetUniformLocation(glsl, 'LifeTime'), particlelifetime);

        % Assign simulated gravity constant 'g' for proper trajectory:
        glUniform1f(glGetUniformLocation(glsl, 'Acceleration'), 1.5);

        % Done with setup:
        glUseProgram(0);
    end

    if ~ismember(stereoMode, [6,7,8,9])
        % Assign random RGB colors to the particles: The shader will use these, but
        % also assign an alpha value that makes the particles "fade out" at the end
        % of there lifetime:
        particlecolors = rand(3, nparticles);
    else
        particlecolors = ones(3, nparticles) * 0.8;
    end
    
    % Maximum speed for particles:
    maxspeed = 1.25;

    % Per-component speed: We select these to shape the fountain in our wanted
    % direction:
    vxmax = maxspeed;
    vymax = maxspeed;
    vzmax = 0.4 * maxspeed;

    % Assign random velocities in (vx,vy,vz) direction: Intervals chosen to
    % shape the beam into something visually pleasing for a teapot:
    particlesxyzt(1,:) = RandLim([1, nparticles],    0.7, +vxmax);
    particlesxyzt(2,:) = RandLim([1, nparticles],    0.7, +vymax);
    particlesxyzt(3,:) = RandLim([1, nparticles], -vzmax, +vzmax);

    % The w-component (4th dimension) encodes the birthtime of the particle. We
    % assign random birthtimes within the possible particlelifetime to get a
    % nice continuous stream of particles. Well, kind of: The flowfactor
    % controls the "burstiness" of particle flow. A value of 1 will create a
    % continous stream, whereas smaller values will create bursts of particles,
    % as if the teapot is choking:
    particlesxyzt(4,:) = RandLim([1, nparticles], 0.0, particlelifetime * flowfactor);

    % Get duration of a single frame:
    ifi = Screen('GetFlipInterval', win);

    % Initial flip to sync us to VBL and get start timestamp:
    vbl = Screen('Flip', win);
    tstart = vbl;
    telapsed = 0;

    % Manually enable 3D mode:
    Screen('BeginOpenGL', win);

    % Enable lighting:
    glEnable(GL.LIGHTING);

    % Enable proper occlusion handling via depth tests:
    glEnable(GL.DEPTH_TEST);

    % Set light position:
    glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);

    % Manually disable 3D mode.
    Screen('EndOpenGL', win);

    % We start with an empty dot array 'xyz' in first frame:
    xyz = [];

    % 3D Dots animation loop: Runs until keypress:
    while ~KbCheck
        % I a stereo display mode, we render the scene for both eyes:
        for view = 0:stereoViews
            % Select 'view' to render (left- or right-eye):
            Screen('SelectStereoDrawbuffer', win, view);

            % Manually reenable 3D mode in preparation of eye draw cycle:
            Screen('BeginOpenGL', win);

            % Setup camera for this eyes 'view':
            glMatrixMode(GL.MODELVIEW);
            glLoadIdentity;

            % This is a bit faked. For a proper solution see help for
            % moglStereoProjection:
            gluLookAt(-0.4 + view * 0.8 , 0, 10, 0, 0, 0, 0, 1, 0);

            % Clear color and depths buffers:
            glClear;

            % Bring a bit of extra spin into this :-)
            glRotated(10 * telapsed, 0, 1, 0);
            glRotated(5  * telapsed, 1, 0, 0);

            % Draw a solid teapot of size 1.0:
            glutSolidTeapot(1.0);

            % For drawing of dots, we need to respecify the light source position,
            % but this must not apply to other objects like the teapot. Therefore
            % we first backup the current lighting settings...
            glPushAttrib(GL.LIGHTING_BIT);

            % ... then set the new light source position ...
            glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);

            % Draw dots of random dot cloud quickly: Common dotdiameter is 5
            % pixels, point smoothing is on, but this time we don't set a dotcolor
            % at all. This way the color can be determined by OpenGL's lighting
            % calculations:
            moglDrawDots3D(win, xyz, 5, [], [], 1);

            % Compute simulation time for this draw cycle:
            telapsed = vbl - tstart;

            if shadingavail
                % Draw the particle fountain. We use a vertex shader in the shader
                % program glsl to compute the physics:
                glUseProgram(glsl);

                % Assign updated simulation time to shader:
                glUniform1f(glGetUniformLocation(glsl, 'Time'), telapsed);

                % Draw the particles. Here particlesxyzt does not encode position,
                % but speed vectors -- this because our shader interprets positions
                % as velocities!
                moglDrawDots3D(win, particlesxyzt, particleSize, particlecolors, [], 1);

                % Done with shaded drawing:
                glUseProgram(0);
            end

            % ... restore old light settings from backup ...
            glPopAttrib;

            % Manually disable 3D mode before calling Screen('Flip')!
            Screen('EndOpenGL', win);

            % Repeat for other eyes view if in stereo presentation mode...
        end

        % Mark end of all graphics operation (until flip). This allows GPU to
        % optimize its operations:
        Screen('DrawingFinished', win, 2);

        % Create uniform random distribution of 3D dots inside a cube for next
        % frame. We do it here after the Screen('DrawingFinished') command, so
        % Matlab can compute this random stuff while the GPU is drawing the dot
        % clouds etc. --> Parallelization allows for potential speedup.
        xyz = [RandLim([1,ndots], -1, 1) ; RandLim([1,ndots], -1, -0.7) ; RandLim([1,ndots], -1, 1)];

        % Show'em: We don't clear the color buffer here, as this is done in
        % next iteration via glClear() call anyway:
        vbl = Screen('Flip', win, vbl + 0.5 * ifi, 2);

        % A mouse button press will pause the animation:
        [x,y,buttons] = GetMouse;
        if any(buttons)
            % Wait for a single mouse click to continue:
            GetClicks;
        end
    end

    % Now a benchmark run to test different strategies for their speed...

    KbReleaseWait;
    Screen('Flip', win);

    if shadingavail
        maxrendermode = 2;
    else
        maxrendermode = 0;
    end

    for rendermode=0:maxrendermode
        switch(rendermode)
            case 0,
                msgtxt = 'Testing now Matlab + CPU animation.';
            case 1,
                msgtxt = 'Testing now vertex shader GPU animation.';
            case 2
                msgtxt = 'Testing now optimized vertex shader GPU animation by use of display lists.';
        end

        DrawFormattedText(win, [msgtxt '\nMax test duration will be 20 seconds.\nPress ESCape key to continue and to finish a subdemo.'], 'center', 'center', [255 255 0]);
        Screen('Flip', win);
        KbStrokeWait;

        % Initial flip to sync us to VBL and get start timestamp:
        vbl = Screen('Flip', win);
        tstart = vbl;
        fc = 0;

        Screen('BeginOpenGL', win);
        glDisable(GL.LIGHTING);

        if rendermode == 2
            % Predraw the particles. Here particlesxyzt does not encode position, but
            % speed -- this because our shader interprets positions as velocities!
            gld = glGenLists(1);
            glNewList(gld, GL.COMPILE);
            moglDrawDots3D(win, particlesxyzt, particleSize, particlecolors, -StartPosition, 1);
            glEndList;
        end

        Screen('EndOpenGL', win);

        % For the fun of it, a little shoot-out between a purely Matlab + CPU based
        % solution, and two different GPU approaches:
        % 3D Dots animation loop: Runs until keypress or 20 seconds elapsed.
        while ~KbCheck && (vbl - tstart < 20)
            % Manually reenable 3D mode in preparation of eye draw cycle:
            Screen('BeginOpenGL', win);

            % Clear color and depths buffers:
            glClear;

            % Compute simulation time for this draw cycle:
            telapsed = vbl - tstart;

            if rendermode > 0
                % Draw the particle fountain. We use a vertex shader in the shader
                % program glsl to compute the physics:
                glUseProgram(glsl);

                % Assign updated simulation time to shader:
                glUniform1f(glGetUniformLocation(glsl, 'Time'), telapsed);

                if rendermode == 1
                    % Draw the particles. Here particlesxyzt does not encode position, but
                    % speed -- this because our shader interprets positions as velocities!
                    moglDrawDots3D(win, particlesxyzt, particleSize, particlecolors, -StartPosition, 1);
                else
                    % Draw particles, but use display list instead of direct call
                    glCallList(gld);
                end

                % Done with shaded drawing:
                glUseProgram(0);
            else
                % Do it yourself in Matlab:
                t = max( (telapsed - particlesxyzt(4,:)) , repmat(0.0, 1, nparticles) );
                t = mod(t, particlelifetime);

                Acceleration = 1.5;
                vpositions(1:3,:) = (particlesxyzt(1:3,:) .* repmat(t, 3, 1));
                vpositions(2,:)   = vpositions(2,:) - (Acceleration * (t.^2));

                particlecolors(4,:) = 1.0 - (t / particlelifetime);

                moglDrawDots3D(win, vpositions, particleSize, particlecolors, [], 1);
            end

            % Manually disable 3D mode before calling Screen('Flip')!
            Screen('EndOpenGL', win);

            % Show'em: We don't clear the color buffer here, as this is done in
            % next iteration via glClear() call anyway. We swap asap, without sync
            % to VBL as this is a benchmark:
            Screen('Flip', win, [], 2, 2);

            % Need a fake vbl timestamp to keep simulation running:
            vbl = GetSecs;

            % Count of drawn frame:
            fc = fc + 1;
        end

        tend = Screen('Flip', win);
        avgfps = fc / (tend - tstart);

        switch(rendermode)
            case 0,
                msgtxt = 'Matlab + CPU';
            case 1,
                msgtxt = 'Shader + GPU';
            case 2
                msgtxt = 'Shader + GPU + VRAM Display lists';
        end

        fprintf('Average framerate FPS for rendermode %i [%s] is: %f Hz.\n', rendermode, msgtxt, avgfps);
        if rendermode == 2
            glDeleteLists(gld,1);
        end

        % Repeat benchmark for other renderModes:
    end

    % A last demo: Warp Drive!
    KbReleaseWait;

    % Respecify StartPosition for particle flow to "behind origin":
    StartPosition = [0, 0, -60];

    if shadingavail
        % Setup the vertex shader for particle fountain animation:

        % Bind shader so it can be setup:
        glUseProgram(glsl);

        % Assign static 3D startposition for fountain:
        glUniform3f(glGetUniformLocation(glsl, 'StartPosition'), StartPosition(1), StartPosition(2), StartPosition(3));

        % Assign lifetime: 10 x increased for starfield simulation...
        glUniform1f(glGetUniformLocation(glsl, 'LifeTime'), 10 * particlelifetime);
        particlesxyzt(4,:) = 10 * particlesxyzt(4,:);

        % Assign no simulated gravity, i.e., set to zero, so we don't get
        % gravity in space:
        glUniform1f(glGetUniformLocation(glsl, 'Acceleration'), 0.0);

        % Done with setup:
        glUseProgram(0);
    end

    % Reassign random velocities in (vx,vy,vz) direction: Intervals chosen to
    % shape the beam into something visually pleasing for a warp-flight:
    maxspeed = 1;
    particlesxyzt(1,:) = RandLim([1, nparticles],  -maxspeed, +maxspeed);
    particlesxyzt(2,:) = RandLim([1, nparticles],  -maxspeed, +maxspeed);
    particlesxyzt(3,:) = RandLim([1, nparticles],          0, 5 * maxspeed);
    particlesxyzt(3,:) = 5 * maxspeed;
    
    % Initial flip to sync us to VBL and get start timestamp:
    vbl = Screen('Flip', win);
    tstart = vbl;
    telapsed = 0;

    % Manually enable 3D mode:
    Screen('BeginOpenGL', win);

    % Enable lighting:
    glEnable(GL.LIGHTING);

    % Enable proper occlusion handling via depth tests:
    glEnable(GL.DEPTH_TEST);

    % Set light position:
    glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);

    % Manually disable 3D mode.
    Screen('EndOpenGL', win);

    % 3D Dots animation loop: Runs until keypress:
    while ~KbCheck
        % I a stereo display mode, we render the scene for both eyes:
        for view = 0:stereoViews
            % Select 'view' to render (left- or right-eye):
            Screen('SelectStereoDrawbuffer', win, view);

            % Manually reenable 3D mode in preparation of eye draw cycle:
            Screen('BeginOpenGL', win);

            % Setup camera for this eyes 'view':
            glMatrixMode(GL.MODELVIEW);
            glLoadIdentity;

            % This is a bit faked. For a proper solution see help for
            % moglStereoProjection:
            gluLookAt(-0.4 + view * 0.8 , 0, 10, 0, 0, 0, 0, 1, 0);

            % Clear color and depths buffers:
            glClear;

            % Bring a bit of extra spin into this :-)
            glRotated(5 * telapsed, 0, 0, 1);

            % For drawing of dots, we need to respecify the light source position,
            % but this must not apply to other objects like the teapot. Therefore
            % we first backup the current lighting settings...
            glPushAttrib(GL.LIGHTING_BIT);

            % ... then set the new light source position ...
            glLightfv(GL.LIGHT0,GL.POSITION,[ 1 2 3 0 ]);

            % Compute simulation time for this draw cycle:
            telapsed = vbl - tstart;

            if shadingavail
                % Draw the particle fountain. We use a vertex shader in the shader
                % program glsl to compute the physics:
                glUseProgram(glsl);

                % Assign updated simulation time to shader:
                glUniform1f(glGetUniformLocation(glsl, 'Time'), telapsed);

                % Draw the particles. Here particlesxyzt does not encode position,
                % but speed vectors -- this because our shader interprets positions
                % as velocities!
                moglDrawDots3D(win, particlesxyzt, particleSize, particlecolors, [], 1);

                % Done with shaded drawing:
                glUseProgram(0);
            end

            % ... restore old light settings from backup ...
            glPopAttrib;

            % Manually disable 3D mode before calling Screen('Flip')!
            Screen('EndOpenGL', win);

            % Repeat for other eyes view if in stereo presentation mode...
        end

        % Mark end of all graphics operation (until flip). This allows GPU to
        % optimize its operations:
        Screen('DrawingFinished', win, 2);

        % Show'em: We don't clear the color buffer here, as this is done in
        % next iteration via glClear() call anyway:
        vbl = Screen('Flip', win, vbl + 0.5 * ifi, 2);

        % A mouse button press will pause the animation:
        [x,y,buttons] = GetMouse;
        if any(buttons)
            % Wait for a single mouse click to continue:
            GetClicks;
        end
    end
    % Done. Close screen and exit:
    Screen('CloseAll');

    % Reenable all keys for KbCheck:
    RestrictKeysForKbCheck([]);

catch
    sca;
    % Reenable all keys for KbCheck:
    RestrictKeysForKbCheck([]);
    psychrethrow(psychlasterror);
end

return;
