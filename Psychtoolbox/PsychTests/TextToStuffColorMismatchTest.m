function TextToStuffColorMismatchTest
% TextToStuffColorMismatchTest - Test if text is drawn
% in the color it was requested to be drawn.
%
% In the past, this demonstrated a bug in the OSX CoreText renderer,
% as used in released PTB versions far in the past. Nowadays nobody gives
% a crap about Apples crap anymore and we use our own OS independent plugin
% text renderer by default, so this is what we test.
%
% History:
%
% 16-Jun-2015  mk  Derived from code provided by user shandelman116.
% 14-Jul-2015  mk  Adapt to new situation that renderer 0 == CoreText on OSX.
% 14-Nov-2021  mk  Throw out a lot of pointless code for the year 2021.

% Setup
Screen('Preference', 'TextAlphaBlending', 0);

% Disable automatic text anti-aliasing by operating system:
% Except on Windows, to work around AMD OpenGL driver bugs...
Screen('Preference', 'TextAntiAliasing', double(IsWin));
[window, wsize] = Screen('OpenWindow', 0, 0);

centerX = wsize(3)/2;
centerY = wsize(4)/2;

% Place cross hair cursor at screen center:
ShowCursor('CrossHair');
SetMouse(centerX, centerY);

Screen('TextFont', window, 'Arial');
Screen('TextSize', window, 126);

fixcross = sprintf('+');
greys = 0:15:255;

crosslines = [
    -30, 30, 0 0;
    0, 0, -30, 30;
    ];

crossLineCols = [];
crossTextCols = [];

% Drawing fixation cross as lines then as text
for grey = greys
    Screen('DrawLines', window, crosslines, 5, grey, [centerX, centerY]);
    Screen('Flip', window);
    crossLine = Screen('GetImage', window, [], [], [], 1);

    [~, ~, textbounds] = DrawFormattedText(window, fixcross, 'center', 'center', grey);
    Screen('FrameRect', window, 255, textbounds);
    Screen('Flip', window);
    crossText = Screen('GetImage', window, [], [], [], 1);

    % Compile list of the color used for each type from a center pixel
    crossLineCols = [crossLineCols, crossLine(centerY, centerX)]; %#ok<*AGROW> 
    crossTextCols = [crossTextCols, crossText(centerY, centerX)];
end

ShowCursor('Arrow');
sca;

% Restore to standard HQ plugin text renderer 1:
Screen('Preference','TextRenderer', 1);

% Checking for TextRenderer 1
figure;
plot(crossLineCols, crossTextCols, '*-');

if ~isequal(crossLineCols, crossTextCols)
    fprintf('\n\n\nDANGER WILL ROBINSON! Mismatch between requested and rendered text colors on text renderer!!!\n');
    hold on;
    plot((1/255:0.01:1) * 255, 255 * ((1/255:0.01:1).^(1/1.2)))
    hold off;

    crossTextCols %#ok<*NOPRT> 
    crossLineCols
end
