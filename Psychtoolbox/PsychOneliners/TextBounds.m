function bounds=TextBounds(w,text,yPositionIsBaseline,centerTheText)
% bounds=TextBounds(window, string [, yPositionIsBaseline=0][, centerTheText=0])
%
% Returns the smallest enclosing rect for the drawn text, relative to the
% current location. This bound is based on the actual pixels drawn, so it
% incorporates effects of text smoothing, etc. "text" may be a cell array
% or matrix of 1 or more strings. The strings are drawn one on top of
% another, at the same initial position, before the bounds are calculated.
% This returns the smallest box that will contain all the strings. The
% prior contents of the scratch window are lost. Usually it should be an
% offscreen window, so the user won't see it. The scratch window should be
% at least twice as wide and high as the text, to cope with uncertainties
% about text direction (e.g. Hebrew) and some unusual characters that
% extend greatly to the left of their nominal starting point. If you only
% know your nominal text size and number of characters, you might do this
% to create your scratch window:
%
% Get the bounding box.
% textSize=48;
% string='Good morning.';
% yPositionIsBaseline=1; % 0 or 1
% w=Screen('OpenWindow',0,255);
% woff=Screen('OpenOffscreenWindow',w,[],[0 0 2*textSize*length(string) 2*textSize]);
% Screen(woff,'TextFont','Arial');
% Screen(woff,'TextSize',textSize);
% t=GetSecs;
% bounds=TextBounds(woff,string,yPositionIsBaseline)
% fprintf('TextBounds took %.3f seconds.\n',GetSecs-t);
% Screen('Close',woff);
%
% Show that it's correct by using the bounding box to frame the text.
% x0=100;
% y0=100;
% Screen(w,'TextFont','Arial');
% Screen(w,'TextSize',textSize);
% Screen('DrawText',w,string,x0,y0,0,255,yPositionIsBaseline);
% Screen('FrameRect',w,0,InsetRect(OffsetRect(bounds,x0,y0),-1,-1));
% Screen('Flip',w);
% Speak('Click to quit');
% GetClicks;
% Screen('Close',w);
%
% The suggested window size in that call is generously large because there
% aren't any guarantees from the font makers about how big the text might
% be for a specified point size. Set your window's font, size, and
% (perhaps) style before calling TextBounds.
%
% Be warned that TextBounds and TextCenteredBounds are slow (taking many
% seconds) if the window is large. They use the whole window, so if the
% window is 1024x1024 they process a million pixels. The two slowest calls
% are Screen 'GetImage' and FIND. Their processing time is proportional to
% the number of pixels in the window. So keep your window small.
%
% Also see Screen 'TextBounds'.
% Also see TextCenteredBounds.

% 9/1/98   dgp wrote it.
% 3/19/00  dgp debugged it.
% 11/17/02 dgp Added fix, image1(:,:,1),  suggested by Keith Schneider to
%              support 16 and 32 bit images.
% 9/16/04  dgp Suggest a pixelSize of 1.
% 12/16/04 dgp Fixed handling of cell array.
% 12/17/04 dgp Round x0 so bounds will always be integer. Add comment about speed.
% 1/18/05  dgp Added Charles Collin's two e suggestion for textHeight.
% 1/28/05  dgp Cosmetic.
% 2/4/05   dgp Support both OSX and OS9.
% 12/22/07 mk  Significant rewrite to adapt to current PTB-3.
% 12/16/15 dgp Added yPositionIsBaseline argument.
% 01/10/16 mk  Switch GetImage buffer from backBuffer to drawBuffer for compat with
%              use of onscreen window as scratch window with imaging pipeline active.
% 01/6/17  dgp Added fourth argument to implement the tiny change needed for
%              TextCenteredBounds. This makes TextCenteredBounds a trivial wrapper that
%              automatically tracks improvements made to TextBounds.
% 11/11/18 mk  Avoid cast to char(), so it works with unicode as well on Octave and Matlab.
%              Problem with unicode reported by dgp.

if nargin < 2 || isempty(text)
    error('Require at least 2 arguments. bounds=TextBounds(window, string [, yPositionIsBaseline][, centerTheText])');
end

if nargin < 3 || isempty(yPositionIsBaseline)
    yPositionIsBaseline = 0;
end

if nargin < 4
    centerTheText = 0;
end

white = 1;

% Clear scratch window to background color black:
Screen('FillRect',w,0);

if yPositionIsBaseline
    % Draw text string baseline at location x,y with a wide margin from lower
    % left corner of screen. The left and lower margins accommodate the many
    % fonts with descenders, and the occasional fonts that have fancy capital
    % letters with flourishes that extend to the left of the starting point.
    screenRect=Screen('Rect',w);
    x0=screenRect(1)+min(2*Screen('TextSize',w),RectWidth(screenRect)/20);
    y0=screenRect(4)-min(2*Screen('TextSize',w),RectHeight(screenRect)/2);
else
    % Draw text string with bounding box origin at upper left corner of screen.
    x0=0;
    y0=0;
end

if centerTheText
    x0=(screenRect(1)+screenRect(3))/2;
end

% We've only got one scratch window, so we compute the widths for centering
% in advance, so as not to mess up the accumulation of letters for the
% bounds.
dx=zeros(1,length(text));
if centerTheText
    if iscell(text)
        for i=1:length(text)
            string=text(i);
            bounds=Screen('TextBounds',w,string);
            width=bounds(3);
            dx(i)=-width/2;
        end
    else
        for i=1:size(text,1)
            string=text(i,:);
            bounds=Screen('TextBounds',w,string);
            width=bounds(3);
            dx(i)=-width/2;
        end
    end
end

if iscell(text)
    for i=1:length(text)
        string=text(i);
        Screen('DrawText',w,string,x0+dx(i),y0,white,[],yPositionIsBaseline);
    end
else
    for i=1:size(text,1)
        string=text(i,:);
        Screen('DrawText',w,string,x0+dx(i),y0,white,[],yPositionIsBaseline);
    end
end

% Read back only 1 color channel for efficiency reasons:
image1=Screen('GetImage', w, [], 'drawBuffer', 0, 1);

% Find all nonzero, i.e. non background, pixels:
[y,x]=find(image1(:,:));

% Use coordinates relative to the origin of the DrawText command.
y=y-y0;
x=x-x0;

% Compute their bounding rect and return it:
if isempty(y) || isempty(x)
    bounds=[0 0 0 0];
else
    bounds=SetRect(min(x)-1,min(y)-1,max(x),max(y));
end

return;
