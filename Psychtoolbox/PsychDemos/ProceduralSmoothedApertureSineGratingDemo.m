function ProceduralSmoothedApertureSineGratingDemo(benchmark)
% ProceduralSmoothedApertureSineGratingDemo([benchmark=0])
%
% This demo demonstrates fast drawing of smooth apertured sine wave gratings 
% using procedural shaders. It only works on hardware with support for the GLSL
% shading language, vertex- and fragment-shaders.

% History:
% 12-Dec-2017 modified from ProceduralGaborDemo.m (Ian Andolina).
% 08-Nov-2018 Clarify interaction of alpha-blending with useAlpha and radius=inf,
%             and the sine modulation that is conditionally also applied to the
%             alpha-channel for radius=inf or useAlpha=false. (MK)


% Default to mode 0 - Just a nice demo.
if nargin < 1
    benchmark = 0;
end

% Setup defaults and unit color range:
PsychDefaultSetup(2);

% Disable synctests for this quick demo:
oldSyncLevel = Screen('Preference', 'SkipSyncTests', 2);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% Open a window 
PsychImaging('PrepareConfiguration');
win = PsychImaging('OpenWindow', screenid, 0.5);

% Initial stimulus params for the grating:
virtualSize = 600;
phase = 0;
freq = .04;
tilt = 225;
contrast = 0.5;
% radius of the disc edge
radius = virtualSize / 2;
% smoothing sigma in pixel
sigma = 55;
% use alpha channel for smoothing?
useAlpha = true;
% smoothing method: cosine (0) or smoothstep (1)
smoothMethod = 1;

% Two ways to achieve the same visual appearance for edge smoothing/fadeout,
% once via premultiplication of the hull function, once via alpha blending.
% Note: If we want consistent visual appearance for radius == inf vs. non-inf,
% then we must not enable alpha blending for the radius == inf case, as the
% 'useAlpha' flag only applies for finite radius! For infinite radius, 'useAlpha'
% would get ignored, which would mean to apply the sine wave signal to the alpha
% channel as well as the color channel, so with alpha blending enabled, the sine
% wave signal would get applied twice while forming the final stimulus!
if useAlpha && ~isinf(radius)
    % Enable alpha-blending, store the sine wave signal in color, use alpha
    % channel to store the hull function for fadeout and apply it via blending.
    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
else
    % No alpha-blending: Hull function for fadeout is directly applied via
    % multiplication to the sine wave signal, directly attenuating it.
    % As the sine signal * hull function also applies to the alpha channel,
    % you will find a edge attenuated sine wave in the alpha channel as well,
    % which could be used for complex alpha blending applications iff you use
    % blending and masking creatively, ie. typically with blend functions other
    % than the ones used in this demo.
    Screen('BlendFunction', win, 'GL_ONE', 'GL_ZERO');
end

% Build a procedural gabor texture for a grating with a RGB color offset of 0.5 -- a 50% gray.
sinetex = CreateProceduralSmoothedApertureSineGrating(win, virtualSize, virtualSize,...
          [.5 .5 .5 0], radius, [], sigma, useAlpha, smoothMethod);

% Draw the grating once, just to make sure the gfx-hardware is ready for the
% benchmark run below and doesn't do one time setup work inside the
% benchmark loop:
Screen('DrawTexture', win, sinetex, [], [], tilt, [], [], [], [], [], [phase, freq, contrast, 0]);

% Perform initial flip to gray background and sync us to the retrace:
vbl = Screen('Flip', win);
tstart = vbl;
count = 0;

% Animation loop
while GetSecs < tstart + 10
    count = count + 1;
    % update values:
    phase = phase + 5;

    % Draw the grating:
    Screen('DrawTexture', win, sinetex, [], [], tilt, [], [], [], [], [], [phase, freq, contrast, 0]);

    if benchmark > 0
        % Go as fast as you can without any sync to retrace and without
        % clearing the backbuffer -- we want to measure gabor drawing speed,
        % not how fast the display is going etc.
        Screen('Flip', win, 0, 2, 2);
    else
        % Go at normal refresh rate for good looking gabors:
        Screen('Flip', win);
    end
end

% A final synced flip, so we can be sure all drawing is finished when we
% reach this point; print some stats
tend = Screen('Flip', win);
avgfps = count / (tend - tstart);
fprintf('\nPresented a total of %i frames at ~%.2g FPS...\n',count,avgfps);

% Close window, release all ressources:
sca

% Restore old settings for sync-tests:
Screen('Preference', 'SkipSyncTests', oldSyncLevel);
