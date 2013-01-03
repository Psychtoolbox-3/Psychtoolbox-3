function DrawManuallyAntiAliasedTextDemo(ts)
% DrawManuallyAntiAliasedTextDemo([textSize=96]) - Manually anti-alias text for special purpose applications.
%
% This demo shows how to draw somewhat anti-aliased text efficiently if you
% can't use the operating systems built-in text anti-aliasing for some
% reason.
%
% It draws the text oversized, without anti-aliasing, into an offscreen
% window. Then it draws the offscreen window into the onscreen window,
% applying alpha-blending and bilinear texture filtering, shrinking to 50%
% size. A GLSL shader based gaussian blur operator is applied during
% drawing to reduce aliasing caused high-frequency edges in the texts
% contour.
%
% This needs a modern graphics card to work and is not as efficient as the
% operating systems anti-aliasing, but still reasonably fast if done right.
%
% see also: PsychDemos

% 3/8/04    awi     Wrote it.
% 7/13/04   awi     Added comments section.  
% 9/8/04    awi     Added Try/Catch, cosmetic changes to documentation.
% 1/21/05   awi     Replaced call to GetChar with call to KbWait. 
% 10/6/05   awi	    Note here cosmetic changes by dgp between 1/21/05 and 10/6/05	.
% 08/24/09  mk      Derived from DrawSomeTextDemo.

AssertOpenGL;

% Define default textsize for oversized text:
if nargin < 1
    ts = 96;
end

try
    % Disable automatic text anti-aliasing by operating system:
    oldSetting = Screen('Preference', 'TextAntiAliasing', 0);

    % Choosing the display with the highest display number is
    % a best guess about where you want the stimulus displayed.
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    % Open double-buffered onscreen window with gray (128) background
    % color. We also use the imaging pipeline to enable support for fast
    % offscreen windows. This would also work with a regular
    % Screen('OpenWindow',...) call, but would be slower:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
    w=PsychImaging('OpenWindow', screenNumber, 128);
    
    % Need alpha-blending:
    Screen('Blendfunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Create a convolution shader for a gaussian blur of width 5 and
    % stddev. 1.5. Needs image processing toolbox for fspecial() function, or
    % alternatively compute your own 5 x 5 kernel matrix with a gaussian
    % convolution kernel inside:
    kernel = fspecial('gaussian', 5, 1.5);
    shader = EXPCreateStatic2DConvolutionShader(kernel, 4, 4, 0, 2);

    Screen('TextFont',w, 'Courier New');
    Screen('TextSize',w, ts);
    Screen('TextStyle', w, 1+2);
    
    % Compute bounding box of textstring:
    bbox = ceil(Screen('TextBounds', w, 'Hello World!'));
    
    % Create offscreen window of sufficient size, with a background color
    % that matches the wanted background color, and a alpha value of zero:
    woff = Screen('OpenOffscreenwindow', w, [128 128 128 0], bbox);
    
    % Set same text properties for woff as for w:
    Screen('TextFont',woff, 'Courier New');
    Screen('TextSize',woff, ts);
    Screen('TextStyle', woff, 1+2);
    
    % Draw text to offscreenw window:
    Screen('DrawText', woff, 'Hello World!', 0, 0, 255); %[0, 0, 255, 255]);
    
    % Draw offscreen window with text, reduced to 50% in size, with
    % bilinear texture filtering and alpha-blending enabled. Use the
    % gaussian blur shader to remove high frequencies during drawing to
    % reduce aliasing artifacts:
    Screen('DrawTexture', w, woff, [], OffsetRect(ScaleRect(bbox, 0.5, 0.5), 100, 100), [], [], [], [], shader);
    
    % Same without gaussian blur, just for comparison:
    Screen('DrawTexture', w, woff, [], OffsetRect(ScaleRect(bbox, 0.5, 0.5), 100, 220), [], [], [], [], 0);

    % Show it:
    Screen('Flip',w);

    KbStrokeWait;
    
    % Done.
    Screen('CloseAll');
    
    % Restore anti-aliasing setting:
    Screen('Preference', 'TextAntiAliasing', oldSetting);

catch %#ok<CTCH>
    % This "catch" section executes in case of an error in the "try" section
    % above.  Importantly, it closes the onscreen window if it's open.
    Screen('CloseAll');

    % Restore anti-aliasing setting:
    Screen('Preference', 'TextAntiAliasing', oldSetting);

    psychrethrow(psychlasterror);
end
