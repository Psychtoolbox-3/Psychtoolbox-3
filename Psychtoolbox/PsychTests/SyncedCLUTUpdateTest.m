function SyncedCLUTUpdateTest(synced)
% SyncedCLUTUpdateTest
%
% Perceptual test to test synchronization of hardware gamma table
% updates to the vertical refresh and to Screen('Flip').
%
% The command Screen('LoadNormalizedGammaTable') provides the optional parameter
% 'loadOnNextFlip'. If set to 1, Psychtoolbox is requested to defer the upload
% of the specified gamma table until next execution of the Screen('Flip') command,
% i.e. swapping the back- and front buffers and update of the gamma table should
% happen synchronously with each other and with vertical retrace.
%
% This test allows to perceptually check correctness of this mechanism. In a loop
% it shows two alternating stimuli: A gray-level ramp of range 0-255, displayed
% with a hardware lut of *half* intensity range 0-0.5. Then, after next retrace, a
% gray-level ramp of *half* range 0-128 with a hardware lut of *full* intensity range
% 0-1.0. If update of the visual stimulus (the graylevel ramp) and of the corresponding
% clut happens synchronously at the same vertical retrace during the 'Flip' command,
% as expected, then the half intensity stim + full intensity lut or full intensity stim +
% half intensity lut should cancel each other out, resulting in the perception of a
% half-intensity gray-ramp with no or only minimal flicker. If, on the other hand,
% synchronisation of clut update and flip doesn't work properly, then one should see
% a strong flicker of the gray-ramp due to non-matching stimulus and clut.
%
% You can easily produce the visual impression of failed sync by providing the optional
% parameter synced = 0 to enforce non-synced updates.
%
% Even if synced updates work, you may see some artifacts at the top of the display,
% due to the inability of slow gfx-hardware to update the lut quickly enough. Either
% avoid that display area for your stim or buy a faster graphics card. Enabling realtime
% scheduling via Priority() command may also help.
%
% The accuracy of stimulus onset timestamps may be reduced on M$-Windows when using
% this mode of operation.
%

% History:
% 24.07.2006  Written (MK).

% This script calls Psychtoolbox commands available only in OpenGL-based
% versions of the Psychtoolbox. The Psychtoolbox command AssertPsychOpenGL will issue
% an error message if someone tries to execute this script on a computer without
% an OpenGL Psychtoolbox
AssertOpenGL;

if nargin < 1
    synced = 1; % Default to synchronized updates.
end
synced

try
    % Get the list of screens and choose the one with the highest screen number.
    % Screen 0 is, by definition, the display with the menu bar. Often when
    % two monitors are connected the one without the menu bar is used as
    % the stimulus display.  Chosing the display with the highest dislay number is
    % a best guess about where you want the stimulus displayed.
    screens=Screen('Screens');
    screenNumber=max(screens);

    % Open a double buffered fullscreen window:
    w=Screen('OpenWindow',screenNumber);

    % Store a backup copy of original LUT into origLUT.
    origLUT=Screen('ReadNormalizedGammaTable', screenNumber);

    % Build a simple gray-level ramp as a single texture.
    [width, height]=Screen('WindowSize', w);
    s=floor(min(width, height)/2)-1;
    [x,y]=meshgrid(-s:s, -s:s);
    fintex=ones(2*s+1,2*s+1);
    fintex(:,:)=mod(x,255)+1;

    % tex1 is a full intensity ramp...
    tex1=Screen('MakeTexture', w, fintex);
    % ...with a half intensity clut.
    lut1=origLUT * 0.5;

    % tex2 is a half intensity ramp...
    tex2=Screen('MakeTexture', w, fintex * 0.5);
    % ...with a full intensity clut.
    lut2=origLUT * 1;

    % Switch to realtime-mode for smooth animation:
    priorityLevel=MaxPriority(w);
    Priority(priorityLevel);

    % Query monitor refresh interval:
    ifi = Screen('GetFlipInterval', w);
    
    % Initial Flip to sync us to vbl and get a timestamp:
    vbl = Screen('Flip', w);

    i=0;
    t1=GetSecs;

    % Test loop: Runs until keypress or 1000 iterations have passed:
    while (1)
        % Update the hardware CLUT with our lut 1: Defer update
        % until Screen('Flip'), if synced == 1
        Screen('LoadNormalizedGammaTable', w, lut1, synced);

        % Draw corresponding stimulus:
        Screen('DrawTexture', w, tex1);

        % Show at next retrace after given deadline and update clut:
        vbl = Screen('Flip', w, vbl + 1.5 * ifi);

        % Update the hardware CLUT with our lut 2: Defer update
        % until Screen('Flip'), if synced == 1
        Screen('LoadNormalizedGammaTable', w, lut2, synced);

        % Draw corresponding stimulus:
        Screen('DrawTexture', w, tex2);

        % Show at next retrace after given deadline and update clut:
        vbl = Screen('Flip', w, vbl + 1.5 * ifi);

        % Increment framecount:
        i=i+1;

        % Abort after 1000 iterations or on a key-press:
        if KbCheck || (i>1000)
            break;
        end;
    end;

    % Some stats...
    avgfps = i / (GetSecs - t1)

    % Disable realtime scheduling:
    Priority(0);

    % Restore the original origLUT CLUT lookup table:
    Screen('LoadNormalizedGammaTable', screenNumber, origLUT);

    % Close window, we're done:
    Screen('CloseAll');
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Priority(0);

    if exist('origLUT', 'var')
        Screen('LoadNormalizedGammaTable', screenNumber, origLUT);
    end

    Screen('CloseAll');
    psychrethrow(lasterror);
end %try..catch..
