function bounds=TextBounds(w,text)
% bounds=TextBounds(window,string)
%
% Returns the smallest enclosing rect for the drawn text, relative to
% the current location. This bound is based on the actual pixels
% drawn, so it incorporates effects of text smoothing, etc. "text"
% may be a cell array or matrix of 1 or more strings. The strings are
% drawn one on top of another, at the same initial position, before
% the bounds are calculated. This returns the smallest box that will
% contain all the strings. The prior contents of the scratch window
% are lost. Usually it should be an offscreen window, so the user
% won't see it. The scratch window should be at least twice as wide
% and height as the text, because the text to cope with uncertainties
% about text direction (e.g. Hebrew) and some unusual characters that extend
% greatly to the left of their nominal starting point. If you only
% know your nominal text size and number of characters, you might do
% this to create your scratch window:
%
% textSize=24;
% string='Hello world.';
% With 'w' being the handle of the onscreen window, e.g., w=Screen('OpenWindow',0,0);
% woff=Screen(w,'OpenOffscreenWindow',[],[0 0 3*textSize*length(string) 2*textSize]);
% Screen(woff,'TextFont','Arial');
% Screen(woff,'TextSize',textSize);
% Screen(woff,'TextStyle',1); % 0=plain (default for new window), 1=bold, etc.
% bounds=TextBounds(woff,string);
% ...
% Screen(woff,'Close');
%
% The suggested window size in that call is generously large because there
% aren't any guarantees from the font makers about how big the text might
% be for a specified point size. Set your window's font, size, and
% (perhaps) style before calling TextBounds.
%
% Be warned that TextBounds and TextCenteredBounds are slow (taking many
% seconds) if the window is large. They use the whole window, so if the
% window is 1024x1204 they process a million pixels. The two slowest calls
% are Screen 'GetImage' and FIND. Their processing time is proportional to
% the number of pixels in the window.
%
% OSX: Also see Screen 'TextBounds'.
%
% The user interface would be cleaner if this function opened and closed
% its own offscreen window, instead of writing in the user's window.
% Unfortunately, this might cause some prohibitive overhead.
%
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

white = 1;

% Clear scratch window to background color black:
Screen('FillRect',w,0);

% Draw text strings, always with the top-left corner of text bounding box
% in top-left corner of window:
if iscell(text)
    for i=1:length(text)
        string=char(text(i));
        Screen('DrawText',w,string,0,0,white, [], 0);
    end
else
    for i=1:size(text,1)
        string=char(text(i,:));
        Screen('DrawText',w,string,0,0,white, [], 0);
    end
end

% Read back only 1 color channel for efficiency reasons:
image1=Screen('GetImage', w, [], 'backBuffer', 0, 1);

% Search non-zero (==non background) pixels:
[y,x]=find(image1(:,:));

% Compute their bounding rect and return it:
if isempty(y) || isempty(x)
    bounds=[0 0 0 0];
else
    bounds=SetRect(min(x)-1,min(y)-1,max(x),max(y));
end
return;
