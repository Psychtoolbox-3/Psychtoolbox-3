function ProceduralGarboriumDemo(ngabors)
% ProceduralGarboriumDemo([ngabors=200]) -- An aquarium full of cute little procedural gabors!
%
% This demo shows how to use the Screen('DrawTextures') command to draw a
% large number of similar images quickly - in this case, Gabor patches of
% different position, size and orientation. It also shows how to achieve
% fast, mathematically correct, linear superposition of image patches by
% use of alpha blending and the Psychtoolbox imaging pipeline. While you
% could always achieve proper superposition by precomputing large images in
% Matlab, then drawing them as textures, the use of alpha blending allows
% to offload the computations to your graphics hardware - this allows
% speedups by factors of more than 100x in some cases!
%
% Furthermore the demo demonstrates how to use procedural shading to
% compute gabor patches on-the-fly instead of precomputing a prototypical
% gabor patch image within Matlab. The command CreateProceduralGabor()
% allows to generate such a procedural gabor and then draw it with the
% normal 'DrawTexture' or 'DrawTextures' texture drawing commands. However,
% this gabor patch is not a image like a normal texture. Instead it is
% computed during each draw operation by a so called GLSL shader program --
% the gabor formula is evaluated on the fly inside your fast graphics card.
%
% This method has a few advantages over the standard texture based method
% (as demonstrated in GarboriumDemo):
%
% - GPU's are extremely fast at this kind of jobs, so the method is
% significantly faster on modern GPU's, especially for large patches,
% allowing for even higher redraw rates.
%
% - As the formula is evaluated on the fly for each output pixel, there are
% no resampling artifacts, regardless of size of your gabor.
%
% - You can change all interesting stimulus parameters on the fly -- change
% contrast, aspect ratio, spatial constant, frequency, phase, orientation
% etc. for each patch during each redraw cycle without the need to
% recompute any matrices and without any speed penalty.
%
% The downsides of this method are the need for recent graphics hardware
% (GPU's) and the fact that the stimulus definition formula needs to be
% implemented as a shader program in the GLSL language, which is less easy
% to use - and less forgiving of programming errors - than the simple and
% easy Matlab language. This means that you're either restricted to our
% predifined set of primitives, or you'll have some steep learning curve.
% Currently PTB provides gabors and (via CreateProceduralSineGrating())
% sine gratings.
%
% The demo shows "an aquarium" of many cute little gabor patches, each moving
% into a random direction. Sometimes the gabors intersect, in that case
% you'll see a linear superposition of them. The gabors also shift phase
% and pulse (aspect ratio modulation) to make them more life-like -- and to
% demonstrate runtime change of interesting stimulus parameters, of course.
%
% Each frame of the animation contains 'ngabors' patches, ngabors defaults
% to 200. Change the number as first optional parameter 'ngabors' if you want
% to exercise your graphics hardware and cpu.
%
% You can exit the demo by pressing any key on the keyboard.
%
% This demo needs recent graphics hardware with floating point framebuffer
% support: ATI Radeon X1000 and later, NVidia GeForce 6000 and later.

% History:
% 07/08/2007 Written (MK).
% 05/18/2008 Rewritten, beautified, adapted to current PTB (MK).
% 11/24/2014 Add support for gaussian blobs (CreateProceduralGaussBlob()) (MK).

% Setup defaults and unit color range:
PsychDefaultSetup(2);

% Set number of gabor patches to 200 if no number provided:
if nargin < 1 || isempty(ngabors)
    ngabors = 200;
end

fprintf('Will draw %i gabor patches per frame.\n', ngabors);

% Select screen with maximum id for output window:
screenid = max(Screen('Screens'));

% Open a fullscreen, onscreen window with gray background. Enable 32bpc
% floating point framebuffer via imaging pipeline on it, if this is possible
% on your hardware while alpha-blending is enabled. Otherwise use a 16bpc
% precision framebuffer together with alpha-blending. We need alpha-blending
% here to implement the nice superposition of overlapping gabors. The demo will
% abort if your graphics hardware is not capable of any of this.
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
[win winRect] = PsychImaging('OpenWindow', screenid, 0.5);

% Retrieve size of window in pixels, need it later to make sure that our
% moving gabors don't move out of the visible screen area:
[w, h] = RectSize(winRect);

% Query frame duration: We use it later on to time 'Flips' properly for an
% animation with constant framerate:
ifi = Screen('GetFlipInterval', win);

% Enable alpha-blending, set it to a blend equation useable for linear
% additive superposition. This allows to linearly
% superimpose gabor patches in the mathematically correct manner, should
% they overlap. Alpha-weighted source means: The 'globalAlpha' parameter in
% the 'DrawTextures' can be used to modulate the intensity of each pixel of
% the drawn patch before it is superimposed to the framebuffer image, ie.,
% it allows to specify a global per-patch contrast value:
Screen('BlendFunction', win, GL_ONE, GL_ONE);

% Define prototypical gabor patch of 65 x 65 pixels default size: si is
% half the wanted size. Later on, the 'DrawTextures' command will simply
% scale this patch up and down to draw individual patches of the different
% wanted sizes:
si = 32;

% Size of support in pixels, derived from si:
tw = 2*si+1;
th = 2*si+1;

% Initial parameters of gabors:

% Phase of underlying sine grating in degrees:
phase = 0;
% Spatial constant of the exponential "hull"
sc = 10.0;
% Frequency of sine grating:
freq = .05;
% Contrast of grating:
contrast = 10;
% Aspect ratio width vs. height:
aspectratio = 1.0;

% Initialize matrix with spec for all 'ngabors' patches to start off
% identically:
mypars = repmat([phase+180, freq, sc, contrast, aspectratio, 0, 0, 0]', 1, ngabors);
myblobpars = repmat([contrast, sc, aspectratio, 0]', 1, ngabors);

% Build a procedural gabor texture for a gabor with a support of tw x th
% pixels and the 'nonsymetric' flag set to 1 == Gabor shall allow runtime
% change of aspect-ratio:
gabortex = CreateProceduralGabor(win, tw, th, 1);

% Ditto for some gaussian blobs, just for variety. Need a bigger mathematical
% support here to avoid cutoff artifacts (tw * 2, th * 2):
blobtex = CreateProceduralGaussBlob(win, tw * 2, th * 2);

% Draw the gabor and blob once, just to make sure the gfx-hardware is ready for the
% benchmark run below and doesn't do one time setup work inside the
% benchmark loop. The flag 'kPsychDontDoRotation' tells 'DrawTexture' not
% to apply its built-in texture rotation code for rotation, but just pass
% the rotation angle to the 'gabortex' shader -- it will implement its own
% rotation code, optimized for its purpose. Additional stimulus parameters
% like phase, sc, etc. are passed as 'auxParameters' vector to
% 'DrawTexture', this vector is just passed along to the shader. For
% technical reasons this vector must always contain a multiple of 4
% elements, so we pad with three zero elements at the end to get 8
% elements.
Screen('DrawTexture', win, gabortex, [], [], [], [], [], [], [], kPsychDontDoRotation, [phase, freq, sc, contrast, aspectratio, 0, 0, 0]);
Screen('DrawTexture', win, blobtex, [], [], [], [], [], [], [], kPsychDontDoRotation, [contrast, sc, aspectratio, 0]);

% Preallocate array with destination rectangles:
% This also defines initial gabor patch orientations, scales and location
% for the very first drawn stimulus frame:
texrect = Screen('Rect', gabortex);
inrect = repmat(texrect', 1, ngabors);

dstRects = zeros(4, ngabors);
for i=1:ngabors
    scale(i) = 1*(0.1 + 0.9 * randn);
    dstRects(:, i) = CenterRectOnPoint(texrect * scale(i), rand * w, rand * h)';
    dstRects2(:, i) = CenterRectOnPoint(texrect * scale(i), rand * w, rand * h)';
end

% Preallocate array with rotation angles:
rotAngles = rand(1, ngabors) * 360;

% Initially sync us to VBL at start of animation loop.
vbl = Screen('Flip', win);
tstart = vbl;
count = 0;

% Animation loop: Run until any keypress:
while ~KbCheck
    % Step one: Batch-Draw all gabor patches at the positions and
    % orientations and with the stimulus parameters 'mypars',
    % computed during last loop iteration:
    Screen('DrawTextures', win, gabortex, [], dstRects, rotAngles, [], [], [], [], kPsychDontDoRotation, mypars);
    
    % Ditto for the gaussian blobs:
    Screen('DrawTextures', win, blobtex, [], dstRects2, rotAngles, [], [], [], [], kPsychDontDoRotation, myblobpars);

    % Mark drawing ops as finished, so the GPU can do its drawing job while
    % we can compute updated parameters for next animation frame. This
    % command is not strictly needed, but it may give a slight additional
    % speedup, because the CPU can compute new stimulus parameters in
    % Matlab, while the GPU is drawing the stimuli for this frame.
    % Sometimes it helps more, sometimes less, or not at all, depending on
    % your system and code, but it only seldomly hurts.
    % performance...
    Screen('DrawingFinished', win);
    
    % Compute updated positions and orientations for next frame. This code
    % is vectorized, but could be probably optimized even more. Indeed,
    % these few lines of Matlab code are the main speed-bottleneck for this
    % demos animation loop on modern graphics hardware, not the actual drawing
    % of the stimulus. The demo as written here is CPU bound - limited in
    % speed by the speed of your main processor.

    % Compute new random orientation for each patch in next frame:
    rotAngles = rotAngles + 1 * randn(1, ngabors);
    
    % Increment phase-shift of each gabor by 10 deg. per redraw:
    mypars(1,:) = mypars(1,:) + 10;
    
    % "Pulse" the aspect-ratio of each gabor with a sine-wave timecourse:
    mypars(5,:) = 1.0 + 0.25 * sin(count*0.1);
    myblobpars(3,:) = 1.0 + 0.25 * sin(count*0.1);

    % Compute centers of all patches, then shift them in new direction of
    % motion 'rotAngles', use the mod() operator to make sure they don't
    % leave the window display area. Its important to use RectCenterd and
    % CenterRectOnPointd instead of RectCenter and CenterRectOnPoint,
    % because the latter ones would round all results to integral pixel
    % locations, which would make for an odd and jerky movement. It is
    % also important to feed all matrices and vectors in proper format, as
    % these routines are internally vectorized for higher speed.
    [x y] = RectCenterd(dstRects);
    x = mod(x + 0.33 * cos(rotAngles/360*2*pi), w);
    y = mod(y - 0.33 * sin(rotAngles/360*2*pi), h);

    % Recompute dstRects destination rectangles for each patch, given the
    % 'per gabor' scale and new center location (x,y):
    dstRects = CenterRectOnPointd(inrect .* repmat(scale,4,1), x, y);

    % Done. Flip one video refresh after the last 'Flip', ie. try to
    % update the display every video refresh cycle if you can.
    % This is the same as Screen('Flip', win);
    % but the provided explicit 'when' deadline allows PTB's internal
    % frame-skip detector to work more accurately and give a more
    % meaningful report of missed deadlines at the end of the script. Not
    % important for this demo, but here just in case you didn't know ;-)
    vbl = Screen('Flip', win, vbl + 0.5 * ifi);

    % Next loop iteration...
    count = count + 1;
end

% Done. Last flip to take end timestamp and for stimulus offset:
vbl = Screen('Flip', win);

% Print the stats:
count
avgfps = count / (vbl - tstart)

% Close onscreen window, release all ressources:
sca;

% Done.
return;

% Little helper function: Converts degrees to radians, included for
% backward compatibility with old Matlab versions:
function radians = deg2rad(degrees)
    radians = degrees * pi / 180;
return
