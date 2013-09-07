function bounds=TextCenteredBounds(w,text)
% bounds=TextCenteredBounds(window,text)
%
% Returns the smallest enclosing rect for the drawn text, relative to the
% current location. This bound is based on the actual pixels drawn, so it
% incorporates effects of text smoothing, etc. All text is drawn on the
% same baseline, horizontally centered by using the x offset
% -Screen(w,'TextWidth',string)/2. "text" may be a cell array or matrix of
% 1 or more strings. The strings are drawn one on top of another, all
% horizontally centered at the current position, before the bounds are
% calculated. This returns the smallest box that will contain all the
% strings. The prior contents of the scratch window are lost. Usually it
% should be an offscreen window, so the user won't see it. If you only
% know your nominal text size and number of characters, you might do this
%
% w=Screen(-1,'OpenOffscreenWindow',[],[0 0 1.5*textSize*length(string) 2*textSize],1);
%
% The suggested window size in that call is generously large because there
% aren't any guarantees from the font makers about how big the text might
% be for a specified point size. The pixelSize of 1 (the last argument)
% minimizes the memory requirements.
%
% Be warned that TextBounds and TextCenteredBounds are slow (taking many
% seconds) if the window is large. They use the whole window, so if the
% window is 1024x1204 they process a million pixels. The two slowest calls
% are Screen 'GetImage' and FIND. Their processing time is proportional to
% the number of pixels in the window. We haven't checked, but it's very
% likely that processing time is also proportional to pixel size, so we
% suggest using a small pixel size (e.g. 1 bit, using an offscreen
% window).
%
% The user interface would be cleaner if this function opened and closed
% its own offscreen window, without bothering the user. Unfortunately the Mac
% OS takes on the order of a second to open and to close an offscreen
% window, making the overhead prohibitive.
%
% Also see TextBounds and Screen 'TextWidth' and 'DrawText'.

% 9/1/98 dgp wrote it.
% 3/19/00 dgp debugged it.
% 11/17/02 dgp Added fix, image1(:,:,1), suggested by Keith Schneider to
%              support 16 and 32 bit images.
% 9/16/04  dgp Suggest a pixelSize of 1.
% 12/16/04 dgp Fixed handling of cell array.
% 12/17/04 dgp Round x0 so bounds will always be integer. Add comment about speed.
% 2/4/05   dgp Support both OSX and OS9.

white = 1;
Screen(w,'FillRect',0);
r=Screen(w,'Rect');
x0=round((r(RectLeft)+r(RectRight))/2);
y0=round((r(RectTop)+2*r(RectBottom))/3);
if iscell(text)
    for i=1:length(text)
        string=char(text(i));
        bounds=Screen('TextBounds',w,string);
        width=bounds(3);
        Screen('DrawText',w,string,x0-width/2,y0,white);
    end
else
    for i=1:size(text,1)
        string=char(text(i,:));
        bounds=Screen('TextBounds',w,string);
        width=bounds(3);
        Screen('DrawText',w,string,x0-width/2,y0,white);
    end
end
image1=Screen('GetImage', w, [], 'backBuffer', 0, 1);
[y,x]=find(image1(:,:,1));
if isempty(y) || isempty(x)
    bounds=[0 0 0 0];
else
    bounds=SetRect(min(x)-1,min(y)-1,max(x),max(y));
    bounds=OffsetRect(bounds,-x0,-y0);
end
