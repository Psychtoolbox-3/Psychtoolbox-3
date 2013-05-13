function SetStereoBlueLineSyncParameters(win, vPos, hFraction, lineColor)
% Change parameters for drawing of frame-sequential stereo blue line sync lines for shutter glasses (stereo goggles).
%
% SetStereoBlueLineSyncParameters(win [, vPos=max][, hFraction=0.25][, lineColor=[1,1,1]])
%
% Call this function after the win = Screen('OpenWindow',...); call on an
% onscreen window in frame-sequential stereo mode to change the parameters
% of drawing of stereo sync lines, as needed by stereo goggles or shutter
% glasses like, e.g., CrystalEyes glasses. Enable automatic sync lines if
% they aren't enabled already. Sync lines are auto-enabled on OS/X.
%
% All parameters except the onscreen 'win'dowhandle are optional and have
% reasonable builtin defaults:
%
% 'vPos': Vertical position of the sync line: Default to windowheight - 1.
%
% 'hFraction': Sync lines for left eye view are drawn starting on the left
% border with a length of hFraction * windowWidth. Right eye lines are
% drawn with a length of (1-hFraction) * windowWidth.
%
% The default of 0.25 should be ok for at least CrystalEyes goggles
%
% 'lineColor' intensity of the color channels for drawing the line in range
% zero to one: Default is [1,1,1] = max red, green and blue --> A white
% sync line. A setting of [0,0,1] = red off, green off, blue max would
% create a classic blue sync line. However, some CrystalEyes stereoenablers
% had problems detecting the signal this way, so we default to all-white.
%

% History:
% 
% Test if a windowhandle is provided...
if nargin < 1
    error('You must provide the windowhandle for the onscreen window as 1st parameter!');
end

% ... and if it is a valid onscreen window in frame-sequential stereo mode:
if Screen('WindowKind', win) ~= 1
    error('Provided windowhandle is not a valid and open onscreen window!');
end

winfo = Screen('GetWindowInfo', win);
if ~ismember(winfo.StereoMode, [1,11])
    % No frame-sequential mode, no point in having sync lines -> No operation.
    fprintf('SetStereoBlueLineSyncParameters: Info: Provided onscreen window is not switched to frame-sequential stereo mode. Call ignored.\n');
    return;
end

% Parse other arguments, assign defaults if none passed:
if nargin < 2
    vPos = [];
end

% Init parameter string to empty, aka "use builtin defaults":
pString = '';

if ~isempty(vPos)
    % Vertical position of sync line (Defaults to winheight - 1 if omitted):
    pString = [ pString sprintf(':yPosition=%f', vPos)];
end

if nargin < 3
    hFraction = [];
end

if ~isempty(hFraction)
    % Length of left-eye signal line in fraction of window width.
    % Defaults to 0.25 if omitted, ie. 25% of width of onscreen window:
    pString = [ pString sprintf(':hFraction=%f', hFraction)];
end

if nargin < 4
    lineColor = [];
end

if ~isempty(lineColor)
    % Color of sync line (R,G,B) in range 0.0 (channel intensity zero) to
    % 1.0 for max output intensity in that channel. Defaults to (1,1,1) ==
    % max intensity white on all channels:
    pString = [ pString sprintf(':Color=%f %f %f', lineColor(1), lineColor(2), lineColor(3))];    
end

% Ok 'win'dowhandle is fine. Any blue line sync function already applied?
lslot = Screen('HookFunction', win, 'Query', 'LeftFinalizerBlitChain', 'Builtin:RenderStereoSyncLine');
rslot = Screen('HookFunction', win, 'Query', 'RightFinalizerBlitChain', 'Builtin:RenderStereoSyncLine');

% Existing hookslots found?
if (lslot ~= -1) && (rslot ~= -1)
    % Yes. Need to recreate with new settings:

    % Destroy old ones:
    Screen('HookFunction', win, 'Remove', 'LeftFinalizerBlitChain' , lslot);
    Screen('HookFunction', win, 'Remove', 'RightFinalizerBlitChain', rslot);

    % Feedback's important:
    fprintf('SetStereoBlueLineSyncParameters: Info: Changing settings for automatic generation of sync lines for external shutter glasses.\n');

    % Insert new slots at former position of the old ones:
    posstring = sprintf('InsertAt%iBuiltin', lslot);
    Screen('Hookfunction', win, posstring, 'LeftFinalizerBlitChain', 'Builtin:RenderStereoSyncLine', pString);
    posstring = sprintf('InsertAt%iBuiltin', rslot);
    Screen('Hookfunction', win, posstring, 'RightFinalizerBlitChain', 'Builtin:RenderStereoSyncLine', pString);
else
    % No. Need to create one from scratch:
    fprintf('SetStereoBlueLineSyncParameters: Info: Enabling automatic generation of sync lines for external shutter glasses.\n');
    
    % Create new slots at end:
    Screen('Hookfunction', win, 'AppendBuiltin', 'LeftFinalizerBlitChain', 'Builtin:RenderStereoSyncLine', pString);
    Screen('Hookfunction', win, 'AppendBuiltin', 'RightFinalizerBlitChain', 'Builtin:RenderStereoSyncLine', pString);    
end

% Enable finalizer hookchains, if not already enabled:
Screen('HookFunction', win, 'Enable', 'LeftFinalizerBlitChain');
Screen('HookFunction', win, 'Enable', 'RightFinalizerBlitChain');

% Done.
return;
