function SetCompressedStereoSideBySideParameters(win, leftOffset, leftScale, rightOffset, rightScale)
% SetCompressedStereoSideBySideParameters(win, [0, 0], [0.5, 1], [0.5, 0], [0.5, 1]);
%
% Change parameters for over-under stereo to side-by-side stereo so that
% the imaging mode will work in a display that horizontally stretches left
% and right frames to fill the full res window.
%
% This routine is called automatically by PsychImaging() during onscreen
% window creation with its default parameters if ...
%
% PsychImaging('AddTask', 'General', 'SideBySideCompressedStereo');
%
% ... is called during window setup.
%
% You can call the function yourself after PsychImaging('OpenWindow', ...);
% if you want to use non-default composition parameters for side-by-side
% stereo presentation.
%
% This is how compressed side-by-side stereo works:
%
% 1. Usercode draws stimuli into the left-eye and righ-eye framebuffers at
%    full display resolution.
%
% 2. Screen('Flip') compresses those stimuli horizontally into the output
%    framebuffer.
%
% 3. At display time, the display stretches half frames to full width:
%    _________        _________       _________
%    | L | R |   ---> | <-L->  |  +   | <-R->  |
%    ---------        ----------      ----------
%
% Side-by-side compressed images are one of several popular stereo HDMI
% formats.
%
% The stereo mode must have been set to 2 or 3 for this to work, or the
% above mentioned PsychImaging() function must have been called for setup.
%
% Example function call:
%
% SetCompressedStereoSideBySideParameters(win, [0, 0], [0.5, 1], [0.5, 0], [0.5, 1])
%
% Call this function after the win = PsychImaging('OpenWindow',...); call
% on an onscreen window in Top/Bottom stereo mode to change the parameters
% of drawing the stereo views.
%
% All parameters except the onscreen 'win'dowhandle are optional and have
% reasonable builtin defaults:
%
% 'leftOffset' = Top-Left [x,y] offset of left eye framebuffer in relative
% coordinates [0,0] == top-left of framebuffer, [1,0] == 1 display window
% width to the right, [2,0] == 2 display window widths to the right etc.
%
% 'leftScale' = Scaling of left eye image buffer. E.g., [1,1] == Don't
% scale. [0.75, 0.5] scale to 75% of original width, 50% of original
% height.
%
% 'rightOffset', 'rightScale' == Ditto for right eye image.
%
% This function was inspired and co-implemented by David Hoffman.
%

% History:
%
% 23.12.2012  mk  Derived from code by David Hoffman, which was derived
%                 from SetStereoSideBySideParameters.m by mk.

% Test if a windowhandle is provided...
if nargin < 1
    error('You must provide the windowhandle for the onscreen window as 1st parameter!');
end

% ... and if it is a valid onscreen window in top-bottom stereo mode:
if Screen('WindowKind', win) ~= 1
    error('Provided windowhandle is not a valid and open onscreen window!');
end

winfo = Screen('GetWindowInfo', win);
if ~ismember(winfo.StereoMode, [2,3])
    % Only do conversion if we are in top-bottom mode, else abort and do
    % normal processing
    fprintf('SetCompressedStereoSideBySideParameters: Info: Provided onscreen window is not in appropriate stereo mode. Call ignored.\n');
    
    fprintf('WARNING:  *********************************************')
    fprintf('WARNING:  *********************************************')
    fprintf('WARNING:  *********************************************')
    fprintf('WARNING:  *********************************************')
    fprintf('WARNING:  ******Aborting Side_by_Side Conversion*******')
    fprintf('WARNING:  *********************************************')
    fprintf('WARNING:  *********************************************')
    fprintf('WARNING:  *********************************************')
    fprintf('WARNING:  *********************************************')
    
    return;
end

% Query true size of onscreen window in pixels w x h:
[w, h] = Screen('WindowSize', win, 1);

% Parse other arguments, assign defaults if none passed:
if nargin < 2 || isempty(leftOffset)
    leftOffset = [0, 0];
end

if nargin < 3 || isempty(leftScale)
    leftScale = [0.5, 1];
end

if nargin < 4 || isempty(rightOffset)
    rightOffset = [0.5, 0];
end

if nargin < 5 || isempty(rightScale)
    rightScale = [0.5, 1];
end

% Query full specification of processing slot for left eye view shader:
% 'slot' is position in processing chain, others are parameters for the
% operation:
[slot shaderid blittercfg voidptr glsl] = Screen('HookFunction', win, 'Query', 'StereoCompositingBlit', 'StereoCompositingShaderCompressedTop'); %#ok<ASGLU>
if slot == -1
    disp('ERROR: Could not find processing slot for left-eye view!... Please check that StereoMode is set to 2 or 3');
    error('Could not find processing slot for left-eye view!... Please check that StereoMode is set to 2 or 3');
end

% Delete old processing slot from pipeline:
Screen('HookFunction', win, 'Remove', 'StereoCompositingBlit' , slot);

% Define new blitter configuration for changed parameters:
leftOffset(1) = floor(leftOffset(1) * w);
leftOffset(2) = floor(leftOffset(2) * h);
blittercfg = sprintf('Builtin:IdentityBlit:Offset:%i:%i:Scaling:%f:%f', leftOffset(1), leftOffset(2), leftScale(1), leftScale(2));

% Insert modified processing function at old position (slot) in the
% pipeline, effectively replacing the slot:
posstring = sprintf('InsertAt%iShader', slot);
Screen('Hookfunction', win, posstring, 'StereoCompositingBlit', shaderid, glsl, blittercfg);

% Query full specification of processing slot for right eye view shader:
% 'slot' is position in processing chain, others are parameters for the
% operation:
[slot shaderid blittercfg voidptr glsl] = Screen('HookFunction', win, 'Query', 'StereoCompositingBlit', 'StereoCompositingShaderCompressedBottom'); %#ok<ASGLU>
if slot == -1
    disp('ERROR: Could not find processing slot for right-eye view!... Please check that StereoMode is set to 2 or 3');
    error('Could not find processing slot for right-eye view!... Please check that StereoMode is set to 2 or 3');
end

% Delete old processing slot from pipeline:
Screen('HookFunction', win, 'Remove', 'StereoCompositingBlit' , slot);

% Define new blitter configuration for changed parameters:
rightOffset(1) = floor(rightOffset(1) * w);
rightOffset(2) = floor(rightOffset(2) * h);
blittercfg = sprintf('Builtin:IdentityBlit:Offset:%i:%i:Scaling:%f:%f', rightOffset(1), rightOffset(2), rightScale(1), rightScale(2));

% Insert modified processing function at old position (slot) in the
% pipeline, effectively replacing the slot:
posstring = sprintf('InsertAt%iShader', slot);
Screen('Hookfunction', win, posstring, 'StereoCompositingBlit', shaderid, glsl, blittercfg);

% Done.
return;
