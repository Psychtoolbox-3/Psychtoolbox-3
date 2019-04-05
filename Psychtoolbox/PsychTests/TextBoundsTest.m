function TextBoundsTest(string, font, textSize, rotAngle)
% TextBoundsTest([string] [, font] [, textSize] [, rotAngle])
%
% Test Screen('TextBounds') as a fast way of determining the bounding
% rectangle (aka bounding box) of a string of text, as well as the very
% precise but slow TextBounds.m function.
%
% Screen('TextBounds') uses information from the operating system about the
% bounding polygons of a string of text to compute the bounding box. This
% is fast, but it relies on the typesetting algorithms of the OS and the
% information encoded in a font definition file playing well with each
% other. Usually that's the case, but there are a few special fonts which
% either don't encode their own boundaries correctly, or which are typeset
% in an unusual way that the OS facilities can't handle correctly. In such
% cases, textbounds will be way to big or too small and text may appear cut
% off.
%
% The TextBounds() function implemented in TextBounds.m uses a different
% approach: It draws the text string into a provided window - usually an
% offscreen window which is so big that it can easily accomodate the text
% string in question, ie., it is way too big for the string. Then it reads
% back the image of that string and computes the exact bounding box of the
% non-background pixels. This method is obviously independent of OS
% facilities and (faulty) information in the font files - thereby highly
% robust and accurate. However the approach is very slow, so you should
% only use this method if you really need it.
%
% This test script both tests functioning of both methods and provides a
% code example of how to use the information returned by both methods to
% either frame drawn text, or to actually draw unusual difficult text by
% first predrawing it into an offscreen window of sufficient size, then
% blitting the really "meaningful" portions of that offscreen window into
% the real image window.
%
% Optional parameters: (All have meaningful defaults)
% 'string' The text string to draw. Default: "Wordy"
% 'font' Name of the font definition file. Default: "Chicago"
% 'textSize' Size of text font. Default: 64 pts.
% 'rotAngle' Orientation of test-blitted text. Default: 45 degrees.
%

% 2/3/05    dgp Wrote it.
% 10/26/05  awi Cosmetic
% 12/22/07  mk  Mostly rewrote it, trying to keep the original intentions.
% 11/11/18  mk  Cleanup a bit. Switch to colorrange 0-1.

if nargin < 1
    string = [];
end

if isempty(string)
    string='Wordy';
end

if nargin < 2
    font = [];
end

if isempty(font)
    font = 'Chicago';
end

if nargin < 3
    textSize = [];
end

if isempty(textSize)
    textSize = 64;
end

if nargin < 4
    rotAngle = [];
end

if isempty(rotAngle)
    rotAngle = 45;
end

try
    red = [1 0 0];
    green = [0 1 0];
    black = 0;

    % Want normalized 0-1 colorrange:
    PsychDefaultSetup(2);

    % Open a onscreen window for display:
    w=PsychImaging('OpenWindow', 0, 0);

    % Open a offscreen window for use as scratchpad of TextBounds. We
    % make it big enough to certainly contain the text, ie., twice the
    % expected height and 3x the expected width. However make sure its
    % not bigger than the whole screen by clipping it against the
    % screens size:
    maxrect = ClipRect([0 0 3*textSize*length(string) 2*textSize], Screen('Rect', w));
    woff = Screen('OpenOffscreenWindow', w, [], maxrect);

    % Set text font and size for onscreen window:
    Screen(w,'TextFont',font);
    Screen(w,'TextSize',textSize);
    % Set text font and size for offscreen window:
    Screen(woff,'TextFont',font);
    Screen(woff,'TextSize',textSize);
    % Compute text bounds with accurate but slow TextBounds.m function:
    bounds=TextBounds(woff,string);
    fprintf('TextBounds says the bounds of ''%s'' are [%d %d %d %d].\n',string,bounds);

    % Select target rectangle on screen for text display:
    x=100;
    y=200;
    dstbounds=OffsetRect(bounds,x,y);

    % Fill target area with a red rectangle -- red background:
    Screen('FillRect',w,red,dstbounds);
    % Draw text directly over the red rectangle:
    Screen('DrawText',w,string,x,y,black);

    % Draw the same text by blitting the relevant area of the 'woff'
    % offscreen window to the screen. As 'TextBounds' used white text
    % on black background internally, we'll see that onscreen:
    dstbounds = OffsetRect(bounds, x, y+200);
    Screen('DrawTexture', w, woff, bounds, dstbounds, rotAngle);

    % Draw help text:
    Screen(w,'DrawText','Hit any key to continue',x,y+400,red);
    Screen('Flip', w);
    KbStrokeWait;

    y=y+200;
    bounds=Screen('TextBounds',w,string);
    fprintf('Screen ''TextBounds'' says the bounds of ''%s'' are [%d %d %d %d].\n',string,bounds);
    bounds=OffsetRect(bounds,x,y);
    Screen('FillRect',w,red,bounds);
    Screen('DrawText',w,string,x,y,black);
    Screen('DrawText',w,'Hit any key to continue',x,y+200,red);
    Screen('Flip',w);
    KbStrokeWait;

    Screen('CloseAll');
catch
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
