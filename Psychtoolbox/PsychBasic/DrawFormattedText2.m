function [nx, ny, textbounds, cache, wordbounds] = DrawFormattedText2(varargin)
% [nx, ny, textbounds, cache, wordbounds] = DrawFormattedText2(tstring, key-value pairs)
% or:
% [nx, ny, textbounds, cache, wordbounds] = DrawFormattedText2(cache, key-value pairs)
% 
% When called with a string, the following key-value pairs are understood:
% win [, sx][, sy][, xalign][, yalign][, xlayout][, baseColor][, wrapat][, transform][, vSpacing][, righttoleft][, winRect][, resetStyle][, cacheOnly]
% Those enclosed in square braces are optional.
% 
% When called with a cache struct, the following optional key-value pair
% arguments are accepted:
% [, win][, sx][, sy][, xalign][, yalign][, transform][, winRect]
% 
% example call:
% DrawFormattedText2('test text', 'win',window_pointer,'baseColor',[255 0 0])
% 
% Draws a string of text 'tstring' into Psychtoolbox window 'win'. Allows
% some formatting and precise positioning. Only works with the FTGL-based
% text plugin (see help DrawTextPlugin), activate it with
% Screen('Preference','TextRenderer', 1);
% This function is _not_ made to be fast. Test carefully if you want it to
% work within a screen refresh if you throw a lot of text at it.
%
% The text string 'tstring' may contain newline characters '\n'.
% Whenever a newline character '\n' is encountered, a linefeed and
% carriage return is performed, breaking the text string into lines.
%
% The text may also contain the formatting tags listed below. Each tag
% changes the formatting of the text from that point onward. Tags do not
% need to be closed. At the exit of the function, the state of the font
% renderer (font, text size, etc) is reset to the state upon function entry
% or the state stored in the cache if drawing from cache.
% - <i>:                toggles italicization
% - <b>:                toggles bolding
% - <u>:                toggles underlining
% - <color=colorFmt>    switches to a new color
% - <font=name>         switches to a new font
% - <size=number>       switches to a new font size
% The <i>, <b> and <u> tags toggle whether text is italicized, bolded or
% underlined and remain active (possibly until the end of the input string)
% until the same <i>, <b> or <u> tag is encountered again.
% The <color>, <font> and <size> tags can be provided empty (i.e., without
% argument), in which case they cause to revert back to the color, font or
% size active before the previous switch. Multiple of these in a row go
% back further in history (until start color, font, size is reached).
% To escape a tag, prepend it with a slash, e.g., /<color>. If you want a
% slash right in front of a tag, escape it by making it a double /:
% //<color>. No other slashes should be escaped.
% <color>'s argument can have one of two formats, it can either be a
% hexadecimal string (HEX), or a comma-separated floating point array
% (FPN). If a HEX input is provided, it should encode colors using a 1, 2,
% 6 or 8 hexadecimal digit string. If 1 or 2 digits are provided, this is
% interpreted as a grayscale value. 6 hexadecimal digits are interpreted as
% R, G and B values (2 digits each). 8 digits further includes an alpha
% value. Using the HEX values, color values range from 0 to 255. Example:
% <color=ff0000> changes the text color to red (ff corresponds to 255). The
% FPN format consists of 1, 3 or 4 comma separated floating point color
% values (luminance, RGB, or RGBA) that can take on any value, but will be
% processed according to the current color settings for the window as
% indicated by Screen('ColorRange'). Typically, the values provided will
% range from 0.0--1.0 per component. Floating point values should include
% the decimal point to ensure that single element FPN values are parsed
% correctly. Examples: <color=1.,0.,0.> and <color=.5>. If a floating point
% value is specified but Screen('ColorRange') returns 255, indicating 8bit
% color values are used for the window, the floating point color values
% provided are automatically scaled and rounded to the 0-255 range. If the
% HEX format is used but Screen('ColorRange') returns 1.0, indicating
% floating point color values are used for the window, the uint8 color
% values are automatically converted to the 0.0-1.0 range.
% A size/font command before a newline can change the height of the line on
% which it occurs. So if you want to space two words 'test' and 'text'
% vertically by white space equivalent to an 80pts line, use:
% 'test\n<size=80>\n<size>text'. There is an empty line between the two new
% lines, and the size=80 says that this line has height of 80pts in the
% selected font.
%
% 'sx' and 'sy' provide a location on the screen with respect to which the
% textbox is positioned. 'sx' and 'sy' can be a pixel location provided as
% a number. Alternatively, 'sx' can also be 'left' (default), 'center', or
% 'right' signifying the left side, horizontal middle, or right side of the
% window rect. 'sy' can also be 'top' (default), 'center', or 'bottom'
% signifying the top, vertical middle, or bottom of the window rect. 
%
% 'xalign' and 'yalign' are text strings indicating how the textbox should
% be aligned to the screen location provided in 'sx' and 'sy'. For
% 'xalign', 'left' (default), signifies that the left of the text's
% bounding box is aligned to the screen location; 'center' that the
% bounding box is centered on this location; and 'right' that the right of
% the text's bounding box is aligned with this location. For 'yalign',
% 'top' (default), signifies that the top of the text's bounding box is
% aligned to the screen location; 'center' that the bounding box is
% centered on this location; and 'bottom' that the bottom of the text's
% bounding box is aligned with this location.
%
% 'xlayout' indicates how each line is positioned horizontally in the
% bounding box. If 'left' (default), all lines are aligned to the left of
% the text's bounding box; if 'center', all lines are centered in the
% bounding box; and if 'right', all lines are aligned to the right of the
% bounding box. Justification options are currently not supported.
%
% 'baseColor' is the color in which the text will be drawn (until changed
% by a <color> format call. This is also the color that will remain active
% after invocation of this function.
%
% 'wrapat', if provided, will automatically break text strings longer than
% 'wrapat' characters into newline separated strings of roughly 'wrapat'
% characters. This is done by calling the WrapString function (See 'help
% WrapString'). 'wrapat' mode may not work reliably with non-ASCII text
% strings, e.g., UTF-8 encoded uint8 strings on all systems. It also does
% not necessarily lead to lines of roughly equal length, unless using a
% monospace font.
%
% 'transform' allows transforming the text to be drawn as a whole. The
% 'transform' input is a cell array of key-value parameters, indicating
% which transform to do in which order. The following transform are
% supported:
% 'translate', [dx dy]: translates text by dx horizontally and dy
%                       vertically
% 'flip', axis        : mirrors text horizontally if axis is 1, vertically
%                       if axis is 2, and along both axes if number is 3
%                       (which is equal to a 180 deg rotation)
% 'scale', [sx sy]    : scale text horizontally by sx and vertically by sy.
%                       Set to 1 if you want no scaling.
% 'rotate', angle     : Rotate text by angle (degrees). Note that for the
%                       PTB screen, a positive rotation is clockwise.
% example {'translate',[100 0],'rotate',45} first translates text by 100
% pixels, then rotates it by 45 degree clockwise. Note that the order of
% operations is important. The above is equal to
% {'rotate',45,'translate',100*sqrt([2 2])/2}. Advanced note: for OpenGL
% transform are applied in the reverse order from how they're specified.
% That is not the case for this interface, transform are applied in the
% order specified.
%
% The optional argument 'vSpacing' sets the spacing between the lines.
% Default value is 1.
%
% The optional argument 'winRect' allows to specify a [left top right
% bottom] rectange, in which the text should be placed etc. By default, the
% rectangle of the whole 'win'dow is used.
%
% 'resetStyle'. If true, we reset the base text style to normal before
% interpreting formatting commands that are present in the input text
% string. If not (false), active text style at function entry is taken into
% account when processing style toggle tags
%
% 'cache'. Upon invocation of the function, it provides an optional output
% 'cache'. Providing this output back to DrawFormattedText2 instead of the
% 'tstring' input allows direct drawing of the exact same text without all
% the preprocessing having to be done again, potentially saving significant
% time (and simplifying the call syntax). In this mode, a subset of the
% below arguments can be used to, e.g., draw to a different window or
% reposition the text. This cache can be generated without any actual
% drawing being done by setting the 'cacheOnly' argument to true. The cache
% is an implementation detail and is subject to change at any time.
% When drawing from cache, the text can be repositioned with the 'sx',
% 'sy', 'xalign', 'yalign' and 'winRect' inputs described above. If only
% 'sx' and 'sy' are provided, these are taken to be offsets to move the
% bounding box of the cached text. 'sx' and 'sy' must be numerical in this
% case and cannot be empty. If 'xalign' and/or 'yalign' are provided,
% full-fledged parsing of 'sx', 'sy', 'xalign' and 'yalign' is done, same
% as during a normal call to DrawFormattedText2. The bounding box is then
% repositioned according to these for inputs. The winRect argument is only
% used in this case. It is optional (defaulting to the whole windows), and
% works as described above, specifying the rect to which 'xalign' and
% 'yalign' apply. The input option 'transform' can furthermore be set,
% which appends additional transformations to what is already in the cache.
%
%
% Return variables:
%
% The function returns the new (nx, ny) position of the text drawing cursor
% and the bounding rectangle 'textbounds' of the drawn string. (nx,ny) can
% be used as new start position for connecting further text strings to the
% bottom of the drawn text string. Calculation of textbounds is
% approximative, so it may give wrong results with some text fonts and
% styles on some operating systems, depending on the various settings. The
% optional 'cache' output argument is discussed above.
%
% When rotating by angles that are not a multiple of 90 degrees, the
% bounding box may not be accurate (too large). The returned bounding box
% is the rect that tightly fits the rotated original bounding box, not the
% box that tightly fits the rotated ink of the letters.
%
% The optional return argument 'wordbounds', if assigned in the calling
% function, returns a n-by-4 matrix of per-word bounding boxes. Each row
% defines a [left,top,right,bottom] rectangle with the bounding box of a
% word in the text string, ie. row 1 = first word, row 2 = 2nd word, ...
% white-space characters delimit single words, as do style changes and line-
% feeds, and these delimiters are not taken into account for the bounding box,
% ie. they don't get their own bounding boxes. The white-space separating
% successive words is as defined by the function isspace(tstring), or by a
% change of text style, color, formatting, etc. Use of 'wordbounds' may cause
% a significant slow-down in text drawing, so only assign this return argument
% if you actually need it. A current limitation is that returned bounding boxes
% will be likely incorrect if you apply multiple transformations like 'scale'
% 'rotate', 'translate' and 'flip' at once. A single transformation will work,
% but multiple ones will cause misplaced per word bounding boxes. If you want
% to get proper 'wordbounds' when drawing text from the 'cache' then you must
% assign 'wordbounds' already in the DrawFormattedText2() invocation which
% returns the 'cache', otherwise bounding boxes might be wrong.
%
% One difference in the return values from this function and
% DrawFormattedText is that the new (nx, ny) position of the text drawing
% cursor output is the baseline of the text. So to use (nx,ny) as the new
% start position for connecting further text strings, you need to draw
% these strings with yPositionIsBaseline==true. Another difference is that
% the returned textbounds bounding box includes the height of an empty line
% at the end if the input string ended with a carriage return.
% DrawFormattedText only moved (nx,ny) but did not include the empty line
% in the bounding box. The empty line is also taken into account when
% centering text
%
%
% Further Notes:
% 
% Please note that while positioning and bounding boxes are pixel accurate
% with the fonts tested during development, i cannot guarantee this is the
% case with all fonts you throw at it. Also note that this function is not
% made to be fast. It has a make draw from cache mode so that all
% preprocessing is done only once and stored in a cache from which the text
% can be drawn directly. Nonetheless, especially if you throw a lot of
% formatting at it, this function may still not be fast in this mode. Use
% with caution in timing critical paths (test and measure to know if its
% fast enough). That said, with reasonable inputs it should manage to draw
% text to screen well within the inter frame interval of most screens.
%
% The function employs clipping by default. Text lines that are detected as
% laying completely outside the 'win'dow or optional 'winRect' will not be
% drawn, but clipped away. This allows to draw multi-page text (multiple
% screen heights) without too much loss of drawing speed. If you find the
% clipping to interfere with text layout of exotic texts/fonts at exotic
% sizes and formatting, you can define the global variable...
% global ptb_drawformattedtext2_disableClipping;
% ... and set it like this ...
% ptb_drawformattedtext2_disableClipping = 1;
% ... to disable the clipping.
%
% Regardless of the clipping setting, the optional 3rd return parameter
% 'textbounds' always covers the complete text.
%
% See DrawFormattedText2Demo for a usage example.

% TODO:
% - justification. I have included an interface, but not implemented (or
%   documented it in the help above) yet.
% - Fix per word bounding boxes (optionally grow them to be proper AOIs for
%   eye-tracking) for multiple concatenated 'Transform's.

% History:
% 2015--2017    Written (DCN).
% 7-May-2017    Add support for 'wordbounds' - per word bounding boxes (MK).


global ptb_drawformattedtext2_disableClipping;
global ptb_drawformattedtext2_padthresh;

if isempty(ptb_drawformattedtext2_disableClipping)
    % Text clipping on by default:
    ptb_drawformattedtext2_disableClipping = 0;
end
% Boundinx boxes are always for the whole text, but text out of the screen
% is culling/clipping. User can forcefully disable or enable this clipping.
disableClip = (ptb_drawformattedtext2_disableClipping ~= -1) && ...
              ((ptb_drawformattedtext2_disableClipping > 0));

if isempty(ptb_drawformattedtext2_padthresh)
    % Threshold for skipping of text justification is 33% by default:
    ptb_drawformattedtext2_padthresh = 0.333;
end
padthresh = ptb_drawformattedtext2_padthresh;

assert(Screen('Preference','TextRenderer') == 1, 'DrawFormattedText2 only works with the FTGL based text drawing plugin, but this plugin is not selected activated with Screen(''Preference'',''TextRenderer'',1), or did not load correctly. See help DrawTextPlugin for more information.');

% Optional per-word bounding boxes requested or should generate them
% because cache requested?
if (nargout >= 4) && (~IsOctave || isargout(4) || isargout(5))
    dowordbounds = 1;
else
    dowordbounds = 0;
end

if IsOctave
    % char() casts of unicode values > 255 map to zero, because Octave
    % uses UTF-8 encoding for unicode, instead of UTF-32 as Matlab. We
    % take care of this in the code, but necessary casts() also trigger
    % an out-of-range warning in Octave, which we can't selectively disable,
    % as it lacks a unique warning id (duh!). Therefore disable all warnings
    % on Octave and reenable to previous setting whenever the we exit, and
    % therefore the canary variable reenablewarn goes out of scope:
    warningstate = warning('query');
    warning('off');
    reenablewarn = onCleanup(@() restorewarningstate(warningstate));
end

%% process key-value input
[opt,qCalledWithCache] = parseInputs(varargin,nargout);
if isempty(opt)
    % nothing to do
    [nx, ny, textbounds, cache, wordbounds] = deal([]);
    return;
elseif qCalledWithCache
    cache = opt.cache;
    if isfield(cache,'tex')
        [nx, ny, wordbounds] = deal(cache.nx,cache.ny,cache.wordbounds);
        DoDrawTexture(opt.win,cache.tex.number,cache.bbox,cache.transform);
        if isempty(cache.transform)
            textbounds = cache.bbox;
        else
            % transform BBox and wordbounds to reflect transforms applied by
            % DoDrawSetup
            textbounds = transformBBox(cache.bbox,cache.transform);
            for b=1:size(wordbounds,1)
                wordbounds(b,:)  = transformBBox2(wordbounds(b,:), cache.transform, cache.bbox);
            end
        end
    else
        [nx, ny, textbounds, wordbounds] = DoDraw(opt.win,...
            disableClip,...
            cache.px,...
            cache.py,...
            cache.bbox,...
            cache.subStrings,...
            cache.switches,...
            cache.fmts,...
            cache.fmtCombs,...
            cache.ssBaseLineOff,...
            cache.winRect,...
            cache.previous,...
            cache.righttoleft,...
            cache.transform,...
            cache.wordbounds);
    end
    % done
    return;
end

% normal draw, unpack input arguments
[tstring,win,sx,sy,xalign,yalign,xlayout,baseColor,wrapat,transform,vSpacing,righttoleft,winRect,resetStyle,cacheOnly,cacheMode] = ...
    deal(opt.tstring,opt.win,opt.sx,opt.sy,opt.xalign,opt.yalign,opt.xlayout,opt.baseColor,opt.wrapat,opt.transform,opt.vSpacing,opt.righttoleft,opt.winRect,opt.resetStyle,opt.cacheOnly,opt.cacheMode);


% Need different encoding for returnChar that matches class of input
% tstring:
returnChar = cast(10,class(tstring));

% Convert all conventional linefeeds into C-style newlines.
% But if '\n' is already encoded as a char(10) as in Octave, then
% there's no need for replacement.
if char(10) ~= '\n' 
    newlinepos = strfind(char(tstring), '\n');
    while ~isempty(newlinepos)
        % Replace first occurence of '\n' by ASCII or double code 10 aka 'repchar':
        tstring = [ tstring(1:min(newlinepos)-1) returnChar tstring(min(newlinepos)+2:end)];
        % Search next occurence of linefeed (if any) in new expanded string:
        newlinepos = strfind(char(tstring), '\n');
    end
end



% string can contain HTML-like formatting commands. Parse them and turn
% them into formatting indicators, then remove them from the string to draw
[tstring,fmtCombs,fmts,switches,previous] = getFormatting(win,tstring,baseColor,resetStyle);
% check we still have anything to render after formatting tags removed
if isempty(tstring)
    % Empty text string -> Nothing to do, but assign dummy values:
    [nx, ny]    = Screen('DrawText', win, '');
    textbounds  = [nx, ny, nx, ny];
    wordbounds  = textbounds;
    cache       = [];
    return;
end

% Text wrapping requested? NB: formatting tags are removed above, so
% wrapping is correct. Also NB that WrapString only replaces spaces by
% linebreaks and thus does not alter the length of the string or where
% words are placed in it. Our codes.style and codes.color vectors thus remain
% correct.
if wrapat > 0
    % Call WrapString to create a broken up version of the input string
    % that is wrapped around column 'wrapat'
    tstring = WrapString(tstring, wrapat);
end

% now, split text into segments, either when there is a carriage return or
% when the format changes
% find format changes
qSwitch = any(switches,1);
% find carriage returns (can occur at same spot as format change)
% make them their own substring so we can process format changes happening
% at the carriage return properly.
if dowordbounds
    % For per-word bounding boxes, each space as classified by isspace()
    % counts as a "carriage return" to force the code to split strings into
    % subStrings at word boundaries as well, not only carriage returns or
    % style changes:
    qCRet = tstring==returnChar | isspace(tstring);
else
    % No per word bounding boxes needed:
    qCRet = tstring==returnChar;
end
qCRet = ismember(1:length(tstring),[find(qCRet) find(qCRet)+1]);
% split strings
qSplit = qSwitch|qCRet;
% subStrings = accumarray(cumsum(qSplit).'+1,tstring(:),[],@(x) {x.'});
% own implementation to make sure it works on all platform (not sure how
% well the accumarray trick works on Octave). Not any slower either
strI = cumsum(qSplit).'+1;
subStrings = cell(strI(end),1);
for p=1:strI(end)
    subStrings{p} = tstring(strI==p);
end
% get which format to use for each substring, and what attributes are
% changed (if any)
fmtCombs = fmtCombs(qSplit);
switches = switches(:,qSplit);
% code when to perform linefeeds.
qLineFeed= cellfun(@(x) ~isempty(x) && x(1)==returnChar,subStrings).';
% we have an empty up front if there is a format change or carriage return
% first in the string
if isempty(subStrings{1}) && ~qCRet(1)
    % remove it if it is a switch from the default format, but not if we
    % start with a carriage return
    subStrings(1) = [];
    qLineFeed(1)  = [];
else
    % we also need to know about the substring before the first split
    fmtCombs = [1          fmtCombs];
    switches = [false(4,1) switches];
end
% if trailing carriage return, this should lead to a trailing empty line,
% add it here.
if tstring(end)==returnChar
    subStrings{end+1}   = '';
    fmtCombs(end+1)     = fmtCombs(end);
    switches(:,end+1)   = switches(:,end);
    qLineFeed(end+1)    = false;
end
% remove those linefeeds from the characters to draw
qLineFeed             = logical(qLineFeed); % i saw qLineFeed become a double on some older matlab versions...
subStrings(qLineFeed) = cellfun(@(x) x(2:end),subStrings(qLineFeed),'uni',false);
% NB: keep substrings that are now empty as they still signal linefeeds and
% empty lines, and format changes can still occur for those empty substrings

% get number of lines.
numlines        = length(strfind(char(tstring), char(10))) + 1;
% vectors for metrics (like width and height) for each line, andor each
% substring
lWidth          = zeros(1,numlines);
lWidthOff       = zeros(2,numlines);
lBaseLineSkip   = zeros(1,numlines);
lBaseLineOff    = zeros(2,numlines);
lWidthOffLine   = zeros(3,length(subStrings));
sWidth          = zeros(1,length(subStrings));
px              = zeros(1,length(subStrings));
py              = zeros(1,length(subStrings));
ssBaseLineSkip  = zeros(1,length(subStrings));
ssBaseLineOff   = zeros(2,length(subStrings));
% process each substring, collect info per substring and line
% get which substrings belong to each line
substrIdxs = [0 cumsum(qLineFeed(1:end-1))];
if ~qLineFeed(1)
    substrIdxs = substrIdxs+1;
end
for p=1:numlines
    % get which substrings belong to this line
    qSubStr = substrIdxs==p;
    
    % to get line width and height, get textbounds of each string and add
    % them together
    for q=find(qSubStr)
        % do format change if needed
        if any(switches(:,q))
            fmt = fmts(:,fmtCombs(q));
            DoFormatChange(win,switches(:,q),fmt);
        end
        if isempty(subStrings{q})
            [~,bbox,h]      = Screen('TextBounds', win,           'x',0,0,1,righttoleft);
            xAdv = 0;
        else
            [~,bbox,h,xAdv] = Screen('TextBounds', win, subStrings{q},0,0,1,righttoleft);
        end
        
        % get amount cursor moves when drawing this substring
        sWidth(q) = xAdv;
        
        % get amount the substring's pixels extend below and above baseline
        ssBaseLineOff(:,q) = bbox([2 4]);
        
        % for proper linespacing, get text height (distance between
        % baselines as indicated in the font)
        ssBaseLineSkip(q) = h;
        
        % get info we need to construct precise bounding boxes, bbox is
        % pixel-precise, word length based on xAdv is too long. and first
        % pixels may appear after drawing cursor position, take that into
        % account
        if xAdv
            lWidthOffLine(1,q) = bbox(3)-bbox(1);
            lWidthOffLine(2,q) = bbox(1);
            lWidthOffLine(3,q) = xAdv-bbox(3);
        end
    end
    
    % get width of each line
    lWidth(p)           = sum(sWidth(qSubStr));
    
    % get largest baseline skip for each line
    lBaseLineSkip(p)    = max(ssBaseLineSkip(qSubStr));
    
    % get largest offset of ink from baseline for each line
    lBaseLineOff(:,p)   = [min(ssBaseLineOff(1,qSubStr)) max(ssBaseLineOff(2,qSubStr))];
    
    % get pixel offset for line
    % floor to deal with fractional positions, tiny bit of overlap with
    % pixel still leads to some paint being put there.
    is = find(lWidthOffLine(1,:),1);
    il = find(lWidthOffLine(1,:),1,'last');
    if righttoleft
        % boxes get drawn in opposite order from input text, so:
        % find last with non-zero length
        if ~isempty(il)
            lWidthOff(1,p) = floor(lWidthOffLine(2,il));
        end
        % find first with non-zero length
        if ~isempty(is)
            lWidthOff(2,p) = floor(lWidthOffLine(3,is));
        end
    else
        % find first with non-zero length
        if ~isempty(is)
            lWidthOff(1,p) = floor(lWidthOffLine(2,is));
        end
        % find last with non-zero length
        if ~isempty(il)
            lWidthOff(2,p) = floor(lWidthOffLine(3,il));
        end
    end
end
% don't forget to set style back to what it should be
ResetTextSetup(win,previous,false);

% get pixel precise line widths
lWidthPrecise = lWidth-sum(lWidthOff,1);
% position the bounding box of the whole text
mWidth      = max(lWidthPrecise);
if xlayout>4
    justOff = ceil((xlayout-3)/3)*3;
    if justOff == 6
        % justify to full width of winRect
        mWidth  = winRect(3)-winRect(1);
    end
end
% from observing how Word layouts the text, about .22 of text height is
% below the baseline (meaning .78 is above). It is important to mimic this,
% a line with text that is much larger than the previous line will end up
% too low, and a line with text that is much smaller than the previous line
% will probably end up on top of the previous line's ink.
totHeight   = -lBaseLineOff(1,1) + sum(round((.22*lBaseLineSkip(1:end-1)+.78*lBaseLineSkip(2:end))*vSpacing)) + lBaseLineOff(2,end);
bbox        = [0 0 mWidth totHeight];
bbox        = positionBbox(bbox,sx,sy,xalign,yalign);

% now, figure out where to place individual lines and substrings into this
% bbox
for p=1:numlines
    % get which substrings belong to this line
    qSubStr = substrIdxs==p;
    idxs = find(qSubStr);
    
    % get center of line w.r.t. bbox left edge
    xlayoutLine = xlayout;
    if xlayout>3
        % do some form of justification, setup
        % check if there are spaces so padding would be possible:
        % TODO
        % Required padding less than padthresh fraction of total width? If
        % not we skip justification, as it would lead to ridiculous looking
        % results:
        if lWidth(p) < mWidth * padthresh
            xlayoutLine = xlayout-justOff;  % remove offset like this makes sure we have line aligned to left or to right of bbox, as requested
        end
    end
    switch xlayoutLine
        case 1
            % align to left at sx
            lc  = lWidth(p)/2;
        case {2,4}
            % center or justify line in bbox
            lc  = (bbox(3)-bbox(1))/2;
        case 3
            % align to right of window
            lc  = bbox(3)-bbox(1) - lWidth(p)/2;
    end
    
    if righttoleft
        off = -cumsum(   sWidth(   qSubStr)    ) + lWidth(p)/2;
    else
        off =  cumsum([0 sWidth(idxs(1:end-1))]) - lWidth(p)/2;
    end
    if xlayoutLine==4
        xSpace = 0;
    else
        xSpace = 0;
    end
    px(qSubStr) = lc+off-lWidthOff(1,p);    % NB: we've been calculating with pixel/paint positions, go back to text cursor position to have the pixels end up where we want them
    
    if p>1
        % add baseline skip for current line if not first line, thats the
        % carriage return. See note above about how Word does text layout.
        idx = find(qSubStr,1,'first');
        py(idx:end) = py(idx) + round((.22*lBaseLineSkip(p-1)+.78*lBaseLineSkip(p))*vSpacing);
    else
        % we're drawing with yPositionIsBaseline==true, correct for that
        py(:) = -min(ssBaseLineOff(1,qSubStr));
    end
end
% determine word (actually segment) bounds
qNotEmpty   = sWidth>0 & sum(abs(ssBaseLineOff),1);
wordboundsbase  = zeros(sum(qNotEmpty),4);
wordboundsbase(:,1) = floor(px(qNotEmpty)+lWidthOffLine(2,qNotEmpty));
wordboundsbase(:,2) = floor(py(qNotEmpty)+ssBaseLineOff(1,qNotEmpty));
wordboundsbase(:,3) = ceil( px(qNotEmpty)+lWidthOffLine(2,qNotEmpty)+lWidthOffLine(1,qNotEmpty));
wordboundsbase(:,4) = ceil( py(qNotEmpty)+ssBaseLineOff(2,qNotEmpty));

% check if we're doing drawing via a texture. If so, we need a texture the
% size of the bounding box, and we need to keep (px,py) w.r.t the bounding
% box, not w.r.t. the screen
qDrawToTexture  = nargout >= 3 && cacheMode==1;

% now we have positions and wordbounds w.r.t. the bbox, add bbox position
% to place them in the right place on the screen
if ~qDrawToTexture
    px = px+bbox(1);
    py = py+bbox(2);
    wordboundsbase(:,1) = wordboundsbase(:,1)+bbox(1);
    wordboundsbase(:,2) = wordboundsbase(:,2)+bbox(2);
    wordboundsbase(:,3) = wordboundsbase(:,3)+bbox(1);
    wordboundsbase(:,4) = wordboundsbase(:,4)+bbox(2);
end


%% done processing inputs, do text drawing
% do draw to texture if wanted
if qDrawToTexture
    drawRect = transformBBox(bbox,transform);
    [tex.number,tex.rect] = Screen('OpenOffscreenWindow', win, [0 0 0 0], [0 0 RectWidth(drawRect) RectHeight(drawRect)]);
    ResetTextSetup(tex.number,previous,true);
    [nx, ny, ~, wordbounds] = DoDraw(tex.number,disableClip,px,py,tex.rect,subStrings,switches,fmts,fmtCombs,ssBaseLineOff,winRect,previous,righttoleft,transform,wordboundsbase);
    nx = nx+bbox(1);
    ny = ny+bbox(2);
    wordbounds(:,1) = wordbounds(:,1)+bbox(1);
    wordbounds(:,2) = wordbounds(:,2)+bbox(2);
    wordbounds(:,3) = wordbounds(:,3)+bbox(1);
    wordbounds(:,4) = wordbounds(:,4)+bbox(2);
    textbounds = drawRect;
end
if ~cacheOnly
    % draw to screen
    if qDrawToTexture
        DoDrawTexture(win,tex.number,drawRect,[]);
    else
        [nx, ny, textbounds, wordbounds] = DoDraw(win,disableClip,px,py,bbox,subStrings,switches,fmts,fmtCombs,ssBaseLineOff,winRect,previous,righttoleft,transform,wordboundsbase);
    end
elseif cacheMode~=1
    % don't do any drawing to screen, nor to texture
    [nx,ny] = deal([]);
    % the bbox in the cache is untranslated (we need to know the original
    % after all when drawing). Output transformed one here for user's info,
    % so they know what'll appear on screen eventually.
    textbounds = transformBBox(bbox,transform);
    wordbounds = wordboundsbase;
    for b=1:size(wordbounds,1)
        wordbounds(b,:)  = transformBBox2(wordbounds(b,:), transform, bbox);
    end
end
if nargout>3
    % make cache
    cache.opt = opt;
    cache.win = win;
    if cacheMode==1
        cache.tex       = tex;
        cache.bbox      = drawRect;
        cache.transform = [];           % no need to reapply transforms as they are already "hardcoded" into the texture. But user can add new ones
        cache.nx        = nx;
        cache.ny        = ny;
        cache.wordbounds = wordbounds;  % store transformed wordbounds as they're "hardcoded" into the texture
    else
        cache.px = px;
        cache.py = py;
        cache.bbox = bbox;
        cache.subStrings = subStrings;
        cache.substrIdxs = substrIdxs;
        cache.switches = switches;
        cache.fmts = fmts;
        cache.fmtCombs = fmtCombs;
        cache.ssBaseLineOff = ssBaseLineOff;
        cache.winRect = winRect;
        cache.previous = previous;
        cache.righttoleft = righttoleft;
        cache.transform = transform;
        cache.wordbounds = wordboundsbase;
    end
end
end

% Restore warning() settings to initial at onCleanup():
function restorewarningstate(warningstate)
    warning(warningstate);
end

function [previouswin, IsOpenGLRendering] = DoDrawSetup(win,transform,bbox)
% Is the OpenGL userspace context for this 'windowPtr' active, as required?
[previouswin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');

% OpenGL rendering for this window active?
if IsOpenGLRendering
    % Yes. We need to disable OpenGL mode for that other window and
    % switch to our window:
    Screen('EndOpenGL', win);
end

if ~isempty(transform)
    [xc, yc] = RectCenterd(bbox);
    
    % Make a backup copy of the current transformation matrix for later
    % use/restoration of default state:
    Screen('glPushMatrix', win);
    
    % Translate origin to the geometric center of the text:
    Screen('glTranslate', win, xc, yc);
    
    % apply transforms
    % as OpenGL transform should be specified in reversed order but i
    % don't want to bother the interface with that, apply transform
    % back to front here.
    for p=length(transform)-1:-2:1
        switch transform{p}
            case 'translate'
                Screen('glTranslate', win, transform{p+1}(1), transform{p+1}(2));
            case 'flip'
                % argument is which axis, x (1), y (2), or both (3). Both
                % equals 180 degree rotation.
                if transform{p+1}(1)==1
                    Screen('glScale', win, -1,  1);
                elseif transform{p+1}(1)==2
                    Screen('glScale', win,  1, -1);
                elseif transform{p+1}(1)==3
                    Screen('glScale', win, -1, -1);
                end
            case 'scale'
                Screen('glScale', win, transform{p+1}(1), transform{p+1}(2));
            case 'rotate'
                ang = transform{p+1}(1);
                % note that for the PTB screen, a positive rotation is
                % clockwise. 
                Screen('glRotate', win, ang);
        end
    end
    
    % We need to undo the translation
    Screen('glTranslate', win, -xc, -yc);
end
end

function DoDrawCleanup(win, previouswin, IsOpenGLRendering, transform)
% undo transform if any
if ~isempty(transform)
    Screen('glPopMatrix', win);
end

% If a different window than our target window was active, we'll switch
% back to that window and its state:
if previouswin > 0
    if previouswin ~= win
        % Different window was active before our invocation:
        
        % Was that window in 3D mode, i.e., OpenGL rendering for that window was active?
        if IsOpenGLRendering
            % Yes. We need to switch that window back into 3D OpenGL mode:
            Screen('BeginOpenGL', previouswin);
        else
            % No. We just perform a dummy call that will switch back to that
            % window:
            Screen('GetWindowInfo', previouswin);
        end
    else
        % Our window was active beforehand.
        if IsOpenGLRendering
            % Was in 3D mode. We need to switch back to 3D:
            Screen('BeginOpenGL', previouswin);
        end
    end
end
end

function DoDrawTexture(win,texNum,texDrawRect,transform)
[previouswin, IsOpenGLRendering] = DoDrawSetup(win,transform,texDrawRect);
Screen('DrawTexture',win,texNum,[],texDrawRect);
DoDrawCleanup(win, previouswin, IsOpenGLRendering, transform);
end

function [nx, ny, bbox, wordbounds] = DoDraw(win,disableClip,sx,sy,bbox,subStrings,switches,fmts,fmtCombs,ssBaseLineOff,winRect,previous,righttoleft,transform,wordboundsbase)

[nx,ny]     = deal(nan);
wordbounds  = wordboundsbase;

[previouswin, IsOpenGLRendering] = DoDrawSetup(win, transform, bbox);
if ~isempty(transform)
    % transform BBox and wordbounds to reflect transforms applied by
    % DoDrawSetup
    for b=1:size(wordbounds,1)
        wordbounds(b,:)  = transformBBox2(wordbounds(b,:), transform, bbox);
    end
    bbox = transformBBox(bbox,transform);
end

% Draw the substrings
for p=1:length(subStrings)
    curstring = subStrings{p};
    yp = sy(p);
    xp = sx(p);
    
    % do format change if needed
    if any(switches(:,p))
        fmt = fmts(:,fmtCombs(p));
        DoFormatChange(win,switches(:,p),fmt);
    end
    
    % Perform crude clipping against upper and lower window borders for this text snippet.
    % If it is clearly outside the window and would get clipped away by the renderer anyway,
    % we can safe ourselves the trouble of processing it:
    if ~isempty(curstring) && (disableClip || ((yp + ssBaseLineOff(2,p) >= winRect(2)) && (yp + ssBaseLineOff(1,p) <= winRect(4))))
        % Inside crude clipping area. Need to draw.
        clipOrEmpty = false;
    else
        % Skip this text line draw call, as it is empty or would be clipped
        % away anyway.
        clipOrEmpty = true;
    end
    
    % Any string to draw?
    if ~clipOrEmpty
        % The cursor is positioned (nx,ny output) to allow to continue to
        % print text directly after the drawn text (if you set
        % yPositionIsBaseline==true at least). Basically behaves like
        % printf or fprintf formatting.
        [nx,ny] = Screen('DrawText', win, curstring, xp, yp,[],[],1, righttoleft);

        % for debug, draw bounding box and baseline
        % [~,sbbox,~,xAdv] = Screen('TextBounds', win, curstring, xp, yp, 1, righttoleft);
        % Screen('FrameRect',win,[0 255 0 128],sbbox);
        % Screen('DrawDots',win,[xp+xAdv ny],1,[0 0 255 128]);    % xp+xAdv equal nx, but not necessarily left edge of next bbox
        % Screen('DrawLine',win,[0 255 255],xp,yp,nx,ny);
    end
end

% Our work is done. clean up
% reset text style etc
ResetTextSetup(win,previous,false);

DoDrawCleanup(win, previouswin, IsOpenGLRendering, transform);

end

%% helpers
function [tstring,fmtCombs,fmts,switches,previous] = getFormatting(win,tstring,startColor,resetStyle)
% This function parses tags out of the text and turns them into formatting
% textstyles, colors, font and text sizes to use when drawing.
% allowable codes:
% - <i>                 To toggle italicization
% - <b>                 To toggle bolding
% - <u>                 To toggle underlining
% - <color=HEX or FPN>  To switch to a new color
% - <font=name>         To switch to a new font
% - <size=number>       To switch to a new font size

% get string type, store original as octave can't deal with string values outside uint8 range
tstringOri  = tstring;
tstring     = char(tstring);

% get colorrange of window, to interpret colors
cr = Screen('ColorRange',win);

% get current active text options
previous.style  = Screen('TextStyle', win);
previous.size   = Screen('TextSize' , win);
previous.font   = Screen('TextFont' , win);
% baseColor is given as input as its convenient for user to be able to set
% it and consistent with other text drawing functions.
previous.color  = startColor;    % keep copy of numeric representation (if its hex, its converted to numeric below when checking base.color)

% get starting text options
base = previous;
if resetStyle
    % start with a clean-slate style (no bold, italic, underlined, etc)
    % when interpreting formatting commands
    base.style = 0;
end

% convert color to hex
if isnumeric(base.color)
    if cr==1.0
        % convert to comma separated floating point
        base.color = sprintf('%f,',base.color);
        base.color(end) = [];   % remove trailing comma
    else
        % convert to hex
        base.color = sprintf('%0*X',[repmat(2,size(base.color));base.color]);
    end
else
    % if user provided color input, store in format we can use in PTB later
    if any(previous.color=='.')||any(previous.color==',')
        % user provided floating point input
        previous.color = sscanf(previous.color,'%f,').';
    else
        % user provided hex input
        previous.color = hex2dec(reshape(previous.color,2,[]).').';
        base.color     = upper(base.color); % ensure hex color is uppercase so below logic works correctly
    end
end

% prepare outputs
% these outputs have same length as string to draw and for each character
% indicate its style and color
codes.style     = repmat(base.style,size(tstring));
tables.color    = {base.color};
tables.font     = {base.font};
codes.color     = ones(size(tstring));                  % 1 indicates startColor, all is in startColor unless user provides color tags telling us otherwise
codes.font      = codes.color;                          % 1 indicates default font, all is in default font unless user provides color tags telling us otherwise
codes.size      = repmat(previous.size,size(tstring));

%% first process tags that don't have further specifiers (<b>, <i>, <u>)
% find tag locations. (?<!(?<!/)/) matches tags with zero or more than one
% slashes in front of them
[tagis ,tagie ,tagt ] = regexp(tstring,'(?i)(?<!(?<!/)/)<(i|b|u)>','start','end','tokens');
if ~isempty(tagis)
    % get full text for each tags and indices to where it is in the input
    % string
    tagi  = [tagis; tagie].';
    tagt  = cat(1,tagt{:});
    
    % fill up output, indicating the style code applicable to each
    % character
    if ~isempty(tagt)
        currStyle = codes.style(1);
        for p=1:length(tagt)
            % the below code snippet is a comment, decribing what the line
            % below does
            % switch formatCodes{p}
            %     case 'i'
            %         fBit = log2(2)+1;
            %     case 'b'
            %         fBit = log2(1)+1;
            %     case 'u'
            %         fBit = log2(4)+1;
            % end
            fBit = floor((double(tagt{p})-'b')/7)+1;
            currStyle = bitset(currStyle,fBit,~bitget(currStyle,fBit));
            codes.style(tagi(p,2):end) = currStyle;
        end
    end
    % now mark active formatting commands to be stripped from text
    toStrip = bsxfun(@plus,tagi(:,1),0:2);      % tags are always three characters long
    toStrip = toStrip(:).';
else
    toStrip = [];
end

%% now process tag that have further specifiers (<color=x>, <font=x>, <size=x>)
% find tag locations. also match empty tags. even if only empty tags, we
% still want to remove them. Ill formed tags with equals sign but no
% argument, or tags with tags inside, are not matched. (?<!(?<!/)/) matches
% tags with zero or more than one slashes in front of them
[tagis ,tagie ,tagt ] = regexp(tstring,'(?i)(?<!(?<!/)/)<(color|font|size)=([^<>]+?)>|(?<!(?<!/)/)<(color|font|size)>','start','end','tokens');
if ~isempty(tagis)
    % get full text for each tag and indices to where it is in the input
    % string
    tagi  = [tagis; tagie].';
    
    % use a simple stack/state machine as we need to maintain a history.
    % empty tags means go back to previous color/size/font
    % (crappy stacks, end of array is top of stack)
    colorStack = 1;         % index in tables.color
    fontStack  = 1;         % index in tables.font
    sizeStack  = codes.size(1);
    
    for p=1:size(tagi,1)
        % check if tag has argument
        if ~isscalar(tagt{p})
            switch tagt{p}{1}
                case 'color'
                    color = tagt{p}{2};
                    % check color is valid
                    % detect FPN: if comma or decimal point
                    qComma = color==',';
                    if any(color=='.') || any(qComma)
                        assert(all(isstrprop(color,'digit')|color=='.'|qComma),'DrawFormattedText2: color tag argument must be specified as comma-separated floating point values, or hexadecimal values')
                    else
                        assert(any(length(color)==[1 2 6 8]),'DrawFormattedText2: if color tag argument is a hexidecimal value, it should have length 1, 2, 6, or 8')
                        assert(all(isstrprop(color,'xdigit')),'DrawFormattedText2: color tag argument must be specified as hexadecimal values, or comma-separated floating point values')
                    end
                    % find new color or add to table
                    iColor = find(strcmpi(tables.color,color),1);
                    if isempty(iColor)
                        tables.color{end+1} = upper(color);
                        iColor = length(tables.color);
                    end
                    % add to stack front
                    colorStack(end+1) = iColor; %#ok<AGROW>
                    % mark all next text as having this color
                    codes.color(tagi(p,2):end) = iColor;
                case 'font'
                    font = tagt{p}{2};  % no checks on whether it is valid
                    % find new color or add to table
                    iFont = find(strcmpi(tables.font,font),1);
                    if isempty(iFont)
                        tables.font{end+1} = font;
                        iFont = length(tables.font);
                    end
                    % add to stack front
                    fontStack(end+1) = iFont; %#ok<AGROW>
                    % mark all next text as having this color
                    codes.font(tagi(p,2):end) = iFont;
                case 'size'
                    fsize = str2double(tagt{p}{2});
                    assert(~isnan(fsize),'DrawFormattedText2: size tag argument must be a number')
                    % add to stack front
                    sizeStack(end+1) = fsize; %#ok<AGROW>
                    % mark all next text as having this color
                    codes.size(tagi(p,2):end) = fsize;
            end
        else
            switch tagt{p}{1}
                case 'color'
                    % if not already reached end of history, revert to
                    % previous color for rest of text
                    if ~isscalar(colorStack)
                        % pop color of stack
                        colorStack(end) = [];
                        % mark all next text as having this color
                        codes.color(tagi(p,2):end) = colorStack(end);
                    end
                case 'font'
                    % if not already reached end of history, revert to
                    % previous color for rest of text
                    if ~isscalar(fontStack)
                        % pop font of stack
                        fontStack(end) = [];
                        % mark all next text as having this font
                        codes.font(tagi(p,2):end) = fontStack(end);
                    end
                case 'size'
                    % if not already reached end of history, revert to
                    % previous color for rest of text
                    if ~isscalar(sizeStack)
                        % pop size of stack
                        sizeStack(end) = [];
                        % mark all next text as having this size
                        codes.size(tagi(p,2):end) = sizeStack(end);
                    end
            end
        end
    end
    
    % now mark active formatting commands to be stripped from text (NB:
    % despite growing array, this is faster than something preallocated)
    for p=1:size(tagi,1)
        toStrip = [toStrip tagi(p,1):tagi(p,2)]; %#ok<AGROW>
    end
end
% now strip active formatting commands from text
% add escape slashes from any escaped tags. also when double slashed,
% we should remove one
toStrip = [toStrip regexp(tstring,'(?i)/<(i|b|u|color|font|size)','start')];
tstringOri (toStrip) = [];
codes.style(toStrip) = [];
codes.color(toStrip) = [];
codes.font (toStrip) = [];
codes.size (toStrip) = [];

% replace tstring with tstringOri again in case input was outside char range, so Octave can handle this all just fine..
tstring = tstringOri;

if isempty(tstring)
    % strong was only formatting commands, nothing to draw, ignore
    [fmtCombs,fmts,switches,previous] = deal([]);
    return;
end

% process colors, hex->dec
for p=1:length(tables.color)
    qComma = tables.color{p}==',';
    if any(tables.color{p}=='.') || any(qComma)
        % at this point, all we know is that the string contains digits,
        % decimal points and commas (checked above).
        tables.color{p} = sscanf(tables.color{p},'%f,');
        assert(any(length(tables.color{p})==[1 3 4]),'DrawFormattedText2: if color tag argument is a comma-separated floating point value, it should have length 1, 3, or 4')
        if cr==255
            % scale
            tables.color{p} = round(tables.color{p}.*255);
        end
    else
        % above we made sure all colors are uppercase and valid hex
        % then, convert letter to their numerical value
        % -48 for numbers (ascii<=64)
        % -55 for letters (ascii>64)
        tables.color{p} = tables.color{p}-48-(tables.color{p}>64)*7;
        % then, sum in pairs, while multiplying first of each pair by its base, 16
        tables.color{p} = sum([tables.color{p}(1:2:end)*16;tables.color{p}(2:2:end)]);
        if cr==1.0
            % scale
            tables.color{p} = tables.color{p}./255;
        end
    end
end

% consolidate codes into one, indicating unique combinations. Also produce
% four boolean vectors indicating what changed upon a style change.
% last, output a table that for each unique combination indicates what the
% style, font, color and size are
c = [codes.style; codes.color; codes.font; codes.size];
% where do changes occur?
switches = logical(diff([[previous.style; 1; 1; previous.size] c],[],2));
% if user-provided baseColor, make sure it is applied!
if ~isnumeric(startColor)
    switches(2,1) = true;
end
% get unique formats and where each of these formats is to be applied
% the below is equivalent to:
% [format,~,fmtCombs] = unique(c.','rows');
% format = format.';
% fmtCombs = fmtCombs.';
% but do required functionality myself to be way faster:
if IsOctave()
    [~,i] = sortrows(c.',1:4);
else
    i=sortrowsc(c.',1:4);
end
groupsSortA = [true any(c(:,i(1:end-1)) ~= c(:,i(2:end)),1)];
format = c(:,i(groupsSortA));
fmtCombs = cumsum(groupsSortA);
fmtCombs(i) = fmtCombs;
% build table with info about each unique format combination
fmts = num2cell(format);
% two columns are indices into table, do indexing
fmts(2,:) = tables.color(format(2,:));
fmts(3,:) = tables.font (format(3,:));
% last, store which need to be changed back when drawing finished
previous.changed = logical(diff([[previous.style; 1; 1; previous.size] c(:,end)],[],2));
end


function DoFormatChange(win,switches,fmt)
% rows in switches / columns in format:
% 1: style, 2: color, 3: font, 4: size

% font and style: if we cange font, always set style with the same
% command. Always works and sometimes needed with some exotic fonts
% (see Screen('TextFont?') )
if switches(3)
    Screen('TextFont', win,fmt{3},fmt{1});
elseif switches(1)
    Screen('TextStyle',win,fmt{1});
end
% color, set through this command. drawing commands below do not
% set color
if switches(2)
    Screen('TextColor',win,fmt{2});
end
% size
if switches(4)
    Screen('TextSize',win,fmt{4});
end
end

function ResetTextSetup(win,previous,qDoAll)
if qDoAll || previous.changed(3)
    Screen('TextFont',win,previous.font,previous.style);
elseif previous.changed(1)
    Screen('TextStyle',win,previous.style);
end
if qDoAll || previous.changed(2)
    Screen('TextColor',win,previous.color); % setting the baseColor input, not color before function entered. Consistent with other text drawing functions
end
if qDoAll || previous.changed(4)
    Screen('TextSize',win,previous.size);
end
end

function bbox = positionBbox(bbox,sx,sy,xalign,yalign)

bWidth = bbox(3)-bbox(1);
bHeight= bbox(4)-bbox(2);

switch xalign
    case 1
        xoff = 0;
    case 2
        xoff = -bWidth/2;
    case 3
        xoff = -bWidth;
end
switch yalign
    case 1
        yoff = 0;
    case 2
        yoff = -bHeight/2;
    case 3
        yoff = -bHeight;
end

bbox = OffsetRect(bbox,round(sx+xoff),round(sy+yoff));
end

function bbox = transformBBox(bbox,transform)
if ~isempty(transform)
    [xc, yc] = RectCenterd(bbox);
    bbox = OffsetRect(bbox,-xc,-yc);
    
    % apply transforms
    for p=1:2:length(transform)
        switch transform{p}
            case 'translate'
                bbox = OffsetRect(bbox, transform{p+1}(1), transform{p+1}(2));
            case 'flip'
                % argument is which axis, x (1), y (2), or both (3). Both
                % equals 180 degree rotation. All are no-ops here as
                % bounding box does not change
            case 'scale'
                bbox = ScaleRect(bbox, transform{p+1}(1), transform{p+1}(2));
            case 'rotate'
                ang = transform{p+1}(1);
                % note that for the PTB screen, a positive rotation is
                % clockwise. 
                bbox = [cosd(ang) -sind(ang); sind(ang) cosd(ang)] * bbox([1 3 1 3; 2 4 4 2]);
                bbox = [min(bbox,[],2).' max(bbox,[],2).'];
        end
    end
    
    % We need to undo the translations...
    bbox = OffsetRect(bbox,xc,yc);
end
end

function bbox = transformBBox2(bbox,transform,refbox)
if ~isempty(transform)
    M = eye(3,3);

    [xc, yc] = RectCenterd(refbox);
    M = M * [[1, 0, xc]; [0, 1, yc]; [0, 0, 1]];

    % apply transforms
    for p=1:2:length(transform)
        switch transform{p}
            case 'translate'
                M = M * [[1, 0, transform{p+1}(1)]; [0, 1, transform{p+1}(2)]; [0, 0, 1]];
            case 'flip'
                % argument is flip around which axis, y (1), x (2), or both (3).
                if transform{p+1}(1) == 1
                    M = M * [[-1, 0, 0]; [0, 1, 0]; [0, 0, 1]];
                end

                if transform{p+1}(1) == 2
                    M = M * [[1, 0, 0]; [0, -1, 0]; [0, 0, 1]];
                end

                if transform{p+1}(1) == 3
                    M = M * [[-1, 0, 0]; [0, -1, 0]; [0, 0, 1]];
                end
            case 'scale'
                M = M * [[transform{p+1}(1), 0, 0]; [0, transform{p+1}(2), 0]; [0, 0, 1]];
            case 'rotate'
                ang = transform{p+1}(1);
                % note that for the PTB screen, a positive rotation is
                % clockwise.
                M = M * [[cosd(ang), -sind(ang), 0] ; [sind(ang), cosd(ang) 0] ; [0, 0, 1]];
        end
    end

    % We need to undo the translations...
    M = M * [[1, 0, -xc]; [0, 1, -yc]; [0, 0, 1]];

    % Apply combined transforms:
    v1 = M * [bbox(1); bbox(2); 1];
    v2 = M * [bbox(3); bbox(2); 1];
    v3 = M * [bbox(3); bbox(4); 1];
    v4 = M * [bbox(1); bbox(4); 1];

    % Project back to 2D plane:
    v1 = v1 / v1(3);
    v2 = v2 / v2(3);
    v3 = v3 / v3(3);
    v4 = v4 / v4(3);

    % Make axis-aligned:
    bbox = [min([v1(1),v2(1),v3(1),v4(1)]), min([v1(2),v2(2),v3(2),v4(2)]), max([v1(1),v2(1),v3(1),v4(1)]), max([v1(2),v2(2),v3(2),v4(2)])];
end
end

function [opt,qCalledWithCache] = parseInputs(varargs,nOutArg)

if isempty(varargs) || isempty(varargs{1})
    % Empty text string -> Nothing to do.
    opt                 = [];
    qCalledWithCache    = false;
    return;
elseif isstruct(varargs{1})
    % called with cache
    qCalledWithCache = true;
    qTextureCache    = isfield(varargs{1},'tex');
    opt = struct(...
        'cache',varargs{1},...
        'win',[],...
        'sx' ,[],...
        'sy' ,[],...
        'xalign' ,[],...
        'yalign' ,[],...
        'transform',[],...
        'winRect',[]...
        );
else
    opt = struct(...
        'tstring',varargs{1},...
        'win',[],...
        'sx' ,[],...
        'sy' ,[],...
        'xalign' ,[],...
        'yalign' ,[],...
        'xlayout' ,1,...
        'baseColor',[],...
        'wrapat',0,...
        'transform',[],...
        'vSpacing',1,...
        'righttoleft',0,...
        'winRect',[],...
        'resetStyle',1,...
        'cacheOnly',false,...
        'cacheMode',1 ...
        );
    qCalledWithCache = false;
end


% parse inputs
assert(mod(length(varargs),2)==1,'function should be called with key-value pairs for all inputs except the first.')
for p=2:2:length(varargs)
    assert(isfield(opt,varargs{p}),'option %s not understood',varargs{p})
    opt.(varargs{p}) = varargs{p+1};
end


% further check and process inputs
if qCalledWithCache
    if isempty(opt.win)
        opt.win = opt.cache.win;
    end
    if ~qTextureCache
        ResetTextSetup(opt.win,opt.cache.previous,true);
    end
    % check which if any of the below are set
    ignorep = [isempty(opt.sx) isempty(opt.sy) isempty(opt.xalign) isempty(opt.yalign)];
else
    if isempty(opt.win)
        error('DrawFormattedText2: Windowhandle missing!');
    end
    
    % layout of individual lines within the bounding box
    % default is aligned to left of box
    if ~ischar(opt.xlayout)
        opt.xlayout = 1;
    else
        switch opt.xlayout
            case 'left'
                opt.xlayout = 1;
            case 'center'
                opt.xlayout = 2;
            case 'right'
                opt.xlayout = 3;
                
                % not implemented: 4, 6, 7, 9
            case 'ljustifylongest'
                % justify to width of longest sentence. if line is too short (see ptb_drawformattedtext2_padthresh), left align
                opt.xlayout = 4;    % -3==1
            case 'rjustifylongest'
                % justify to width of longest sentence. if line is too short (see ptb_drawformattedtext2_padthresh), right align
                opt.xlayout = 6;    % -3==3
            case 'ljustifyfullwidth'
                % justify to width of winRect. if line is too short (see ptb_drawformattedtext2_padthresh), left align
                opt.xlayout = 7;    % -6==1
            case 'rjustifyfullwidth'
                % justify to width of winRect. if line is too short (see ptb_drawformattedtext2_padthresh), right align
                opt.xlayout = 9;    % -6==3
            otherwise
                % ignore anything else user may have provided
                opt.xlayout = 1;
        end
    end
    
    % Keep current text color if none provided:
    if isempty(opt.baseColor)
        opt.baseColor = Screen('TextColor', opt.win);
    else
        opt.baseColor = opt.baseColor(:).'; % ensure row vector
    end
    
    % No text wrapping by default:
    if isinf(opt.wrapat)
        opt.wrapat = 0;
    end
    
    % option to only generate cache but not draw
    if ~isempty(opt.cacheOnly)
        opt.cacheOnly = logical(opt.cacheOnly);
    end
    
    % check cache mode
    if ischar(opt.cacheMode)
        switch opt.cacheMode
            case 'texture'
                opt.cacheMode = 1;
            case 'fullCache'
                opt.cacheMode = 2;
            otherwise
                error('cache mode "%s" not understood, possible values: "texture" and "fullCache"',opt.cacheMode);
        end
    else
        assert(ismember(opt.cacheMode,[1 2]),'cacheMode must be 1 or 2')
    end
    
    if opt.cacheOnly && nOutArg < 4
        error('cacheOnly requested but fourth output (cache) is not requested')
    end
end


% No mirroring/rotating/etc of text by default:
if ~isempty(opt.transform)
    checktransform(opt.transform);
end

% Default rectangle for centering/formatting text is the client rectangle
% of the window, but usercode can specify arbitrary override it:
if isempty(opt.winRect)
    opt.winRect = Screen('Rect', opt.win);
    qWinRectSpecified = false;
else
    assert(isnumeric(opt.winRect) && numel(opt.winRect)==4,'Provided winRect is not valid')
    qWinRectSpecified = true;
end

% position text box on screen
[opt.sx,opt.sy,opt.xalign,opt.yalign] = parseTextBoxPositioning(opt.sx,opt.sy,opt.xalign,opt.yalign,opt.winRect);

% do any last extra processing if drawing from cache.
if qCalledWithCache
    if ~all(ignorep) % all true: no repositioning or offsetting
        % move bbox as requested
        bbox = opt.cache.bbox;
        if any(~ignorep(3:4))
            % first, fill out missing arguments.
            for f={'sx','sy','xalign','yalign'}
                if isempty(opt.(f{1}))
                    opt.(f{1}) = opt.cache.opt.(f{1});
                end
            end
            % reposition bbox:
            bbox = positionBbox([0 0 bbox(3)-bbox(1) bbox(4)-bbox(2)],opt.sx,opt.sy,opt.xalign,opt.yalign);
            % bbox(1:2) now contain new top-left for text. Make that into
            % offsets to reposition text
            off = bbox(1:2)-opt.cache.bbox(1:2);
        else
            % first, fill out missing arguments.
            for f={'sx','sy'}
                if isempty(opt.(f{1}))
                    opt.(f{1}) = 0;
                end
            end
            % only sx and sy provided, do offsetting
            assert(isnumeric(opt.sx)&&isnumeric(opt.sy),'When drawing from cache and providing horizontal and vertical offsets with the ''sx'' and ''sy'' inputs, these offsets must be numeric')
            off = [opt.sx opt.sy];
            bbox = OffsetRect(opt.cache.bbox,off(1),off(2));
        end
        % common logic: put bounding box and text in new place
        if ~qTextureCache
            opt.cache.px = opt.cache.px + off(1);
            opt.cache.py = opt.cache.py + off(2);
        end
        % apply offsets to word bounding boxes too
        off = bbox-opt.cache.bbox;
        opt.cache.wordbounds(:,1) = opt.cache.wordbounds(:,1)+off(1);
        opt.cache.wordbounds(:,2) = opt.cache.wordbounds(:,2)+off(2);
        opt.cache.wordbounds(:,3) = opt.cache.wordbounds(:,3)+off(1);
        opt.cache.wordbounds(:,4) = opt.cache.wordbounds(:,4)+off(2);
        opt.cache.bbox = bbox;
    end
    % overwrite winRect in cache if set by user
    if qWinRectSpecified
        opt.cache.winRect = opt.winRect;
    end
    % append to transforms, if any provided by user
    if ~isempty(opt.transform)
        opt.cache.transform = [opt.cache.transform opt.transform];
    end
end
end

function [sx,sy,xalign,yalign] = parseTextBoxPositioning(sx,sy,xalign,yalign,winRect)
xpos = 0;   % default: use provided sx (w.r.t. winRect)
% Default x start position is left of window:
if isempty(sx)
    sx = 0;
else
    % have text specifying a position at the edge of windowrect?
    if ischar(sx)
        if strcmpi(sx, 'left')
            xpos = 1;
        elseif strcmpi(sx, 'center')
            xpos = 2;
        elseif strcmpi(sx, 'right')
            xpos = 3;
        else
            % Ignore any other crap user may have provided, align to left.
            xpos = 1;
        end
    elseif ~isnumeric(sx)
        % Ignore any other crap user may have provided.
        sx = 0;
    end
end

ypos = 0;   % default: use provided sy (w.r.t. winRect)
% Default y start position is top of window:
if isempty(sy)
    sy = 0;
else
    % have text specifying a position at the edge of windowrect?
    if ischar(sy)
        if strcmpi(sy, 'top')
            ypos = 1;
        elseif strcmpi(sy, 'center')
            ypos = 2;
        elseif strcmpi(sy, 'bottom')
            ypos = 3;
        else
            % Ignore any other crap user may have provided, align to left.
            ypos = 1;
        end
    elseif ~isnumeric(sy)
        % Ignore any other crap user may have provided.
        sy = 0;
    end
end

% now process xpos and ypos
switch xpos
    case 0
        % provided sx is in winRect
        sx = sx+winRect(1);
    case 1
        sx = winRect(1);
    case 2
        sx = (winRect(1)+winRect(3))/2;
    case 3
        sx = winRect(3);
end
switch ypos
    case 0
        % provided sy is in winRect
        sy = sy+winRect(2);
    case 1
        sy = winRect(2);
    case 2
        sy = (winRect(2)+winRect(4))/2;
    case 3
        sy = winRect(4);
end


%%% now we have a position, figure out how to position box with respect to
%%% that position
% Default x layout is align box to left of specified position
if isempty(xalign) || ~ischar(xalign)
    xalign = 1;
else
    if strcmpi(xalign, 'left')
        xalign = 1;
    elseif strcmpi(xalign, 'center')
        xalign = 2;
    elseif strcmpi(xalign, 'right')
        xalign = 3;
    else
        % ignore anything else user may have provided
        xalign = 1;
    end
end

% Default y layout is align box below of specified position
if isempty(yalign) || ~ischar(yalign)
    yalign = 1;
else
    if strcmpi(yalign, 'top')
        yalign = 1;
    elseif strcmpi(yalign, 'center')
        yalign = 2;
    elseif strcmpi(yalign, 'bottom')
        yalign = 3;
    else
        % ignore anything else user may have provided
        yalign = 1;
    end
end
end

function checktransform(transform)
if isempty(transform)
    return;
end

assert(mod(numel(transform),2)==0,'transform input must be key-value')
assert(iscellstr(transform(1:2:end)),'transform: all keys should be a character array')
assert(all(cellfun(@isnumeric,transform(2:2:end))),'transform: all values should be a numeric')
end
