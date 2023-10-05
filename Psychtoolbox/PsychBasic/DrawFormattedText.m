function [nx, ny, textbounds, wordbounds] = DrawFormattedText(win, tstring, sx, sy, color, wrapat, flipHorizontal, flipVertical, vSpacing, righttoleft, winRect)
% [nx, ny, textbounds, wordbounds] = DrawFormattedText(win, tstring [, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical][, vSpacing][, righttoleft][, winRect])
%
% Draws a string of text 'tstring' into Psychtoolbox window 'win'. Allows
% some basic formatting.
%
% Please note that all formatting of text and returned information about
% the bounding box of the text is only approximative and may be inaccurate
% for some text strings, fonts, font settings and combination of parameters
% for this function. Centering of text and bounding box calculations work
% best for single line text, vertical centering may be a bit off for multi-
% line text.
%
% The text string 'tstring' may contain newline characters '\n'.
% Whenever a newline character '\n' is encountered, a linefeed and
% carriage return is performed, breaking the text string into lines.
%
% 'sx' defines the left border of the text: If it is left out, text
% starts at x-position zero, otherwise it starts at the specified position
% 'sx' in the window. If sx=='left', then each line of text is left justified
% to the left border of the target window, or of 'winRect', if provided.
% If sx=='center', then each line of text is horizontally centered in
% the window or 'winRect'. If sx=='right', then each line of text is right
% justified to the right border of the target window, or of 'winRect'.
% The options sx == 'wrapat' and sx == 'justifytomax' try to align the start
% of each text line to the left border and the end of each text line to either
% the specified 'wrapat' number of columns, or to the width of the widest line
% of text in the text string, causing all lines of text to appear of roughly
% even width. This is achieved by adjusting the width of blanks between words
% in each line of text to even out text line length. The sx == 'wrapat' and
% sx == 'justifytomax' options are considered experimental. They only work for
% non-flipped, simple text, and even there they may work badly, so don't rely on
% them doing the right thing without assistance of your code, e.g., breaking
% text into lines of reasonably similar length. The justification functions
% will not justify text which deviates by more than ptb_drawformattedtext_padthresh
% from the target width of the text, with ptb_drawformattedtext_padthresh being
% a global variable you can set. It defaults to 0.333 ie., doesn't justify text
% lines if the text lines are more than 33% shorter than the reference line.
% If sx=='centerblock' then text is not horizontally centered individually for
% each line of text, as with 'center', but for multi-line text, the whole block
% of text is roughly centered, by starting all lines left-justified, but then shifting
% this whole block of left justified text horizontally so that the longest line
% in the block of text is horizontally centered in the target window or 'winRect'.
%
% 'sy' defines the baseline of the (first line of the) text. If left out, text starts
% roughly at the top of the window, otherwise it starts at the specified vertical
% pixel position.
%
% If sy=='center', then the whole text is roughly vertically centered in the
% window. 'color' is the color value of the text (color index or [r g b]
% triplet or [r g b a] quadruple). If color is left out, the current text
% color from previous text drawing commands is used. 'wrapat', if provided,
% will automatically break text strings longer than 'wrapat' characters
% into newline separated strings of roughly 'wrapat' characters. This is
% done by calling the WrapString function (See 'help WrapString'). 'wrapat'
% mode may not work reliably with non-ASCII text strings, e.g., UTF-8
% encoded uint8 strings on all systems.
%
% The optional flag 'flipHorizontal' if set to 1 will mirror the text
% horizontally, whereas the optional flag 'flipVertical' if set to 1 will
% mirror the text vertically (upside down).
%
% The optional argument 'vSpacing' sets the spacing between the lines. Default
% value is 1.
%
% The optional argument 'righttoleft' if set to 1, will ask to draw the
% text string in right-to-left reading direction, e.g., for scripts which
% read right to left
%
% The optional argument 'winRect' allows to specify a [left top right bottom]
% rectangle, in which the text should be centered/placed etc. By default,
% the rectangle of the whole 'win'dow is used.
%
% The function employs clipping by default. Text lines that are detected as
% lying completely outside the 'win'dow or optional 'winRect' will not be
% drawn, but clipped away. This allows to draw multi-page text (multiple
% screen heights) without too much loss of drawing speed. If you find the
% clipping to interfere with text layout of exotic texts/fonts at exotic
% sizes and formatting, you can define the global variable...
%
% global ptb_drawformattedtext_disableClipping;
% ... and set it like this ...
% ptb_drawformattedtext_disableClipping = 1;
% ... to disable the clipping.
%
% Clipping also gets disabled if you request the optional 3rd return
% parameter 'textbounds' to ensure correct computation of a bounding box
% that covers the complete text. You can enforce clipping by setting
% ptb_drawformattedtext_disableClipping = -1; however, computed bounding
% boxes will then only describe the currently visible (non-clipped) text.
%
%
% The function returns the new (nx, ny) position of the text drawing cursor
% and the bounding rectangle 'textbounds' of the drawn string. (nx,ny) can
% be used as new start position for connecting further text strings to the
% bottom of the drawn text string. Calculation of these bounds is
% approximative, so it may give wrong results with some text fonts and
% styles on some operating systems, depending on the various settings.
%
% The optional return argument 'wordbounds', if assigned in the calling
% function, returns a n-by-4 matrix of per-word bounding boxes. Each row
% defines a [left,top,right,bottom] rectangle with the bounding box of a
% word in the text string, ie. row 1 = first word, row 2 = 2nd word, ...
% white-space characters delimit single words and are not taken into
% account, the successive words are defined by the function strtok(tstring).
% Per word bounding boxes don't work if flipHorizontal or flipVertical is
% requested.
%
% See DrawFormattedTextDemo for a usage example.

% History:
% 10/17/06  Written (MK).
% 11/01/06  Add support for correct handling of 3D rendering mode (MK).
% 11/22/06  More 3D handling: Save/restore backface cull state (MK).
% 05/14/07  Return a more meaningful end cursor position (printf - semantics) (MK)
% 01/31/09  Add optional vSpacing parameter (Alex Leykin).
% 09/20/09  Add some char() casts, so Octave can handle Unicode encoded text strings as well.
% 01/10/10  Add support for 'righttoleft' flag and for uint8 tstring types (MK).
% 10/28/10  Add crude text clipping/culling, so multi-page text doesn't bog
%           us down completely. Text clearly outside the 'win'dow gets
%           preculled. (MK).
% 02/07/12  Add re-cast operation to output string to make sure the actual
%           string fed into Screen() is of the same datatype as the
%           original input string, e.g., to prevent losing a double()
%           unicode encoding during string processing/formatting. (MK)
% 06/17/13  Add sx == 'right' option for right-alignment of text. (MK)
% 07/02/14  Add sx == 'wrapat' and sx == 'justifytomax' options for block adjustment of text.
%           This is to be considered quite a prototype. (MK)
% 09/21/14  Fix text clipping when used with optional winRect parameter. (MK)
% 06/20/15  Improve text centering for single line text, and bounding box
%           calculations. This is "whack a mole". Making it better for one OS
%           can make it worse for another OS, there is no winning. (MK)
% 07/21/15  Remove forced line breaks at 250 chars on OSX. No longer needed. (MK)
% 07/25/15  Remove Windows & -> && special handling, now handled by Screen('DrawText')
%           via DT_NOPREFIX flag. Suggested by Diederick. (MK)
% 08/08/15  Fix bug for Unicode text introduced when improving single-line text centering. (MK)
% 02/16/16  Improve bounding box calculations for multi-line text, esp. for vSpacing > 1. (MK)
% 03/29/16  Add sx='centerblock', always use yPosIsBaseline, update help text. (MK)
% 08/24/16  Made vertical flip work when using FTGL plugin. (DCN)
% 08/30/16  Fix for vertical flip is now in FTGL plugin. (DCN)
% 02/15/17  Improve text/word vertical alignment if sy == 'centered' and horizontal
%           justification modes (via sx) are requested. We now always use yPosIsBaseline=1
%           and compensate for baseline accordingly. (MK)
%           Improve use of blank space when block adjustment is used. Instead of approximation
%           of a blanks ' ' width by the width of a 'X', use the xAdvance property of 'TextBounds'
%           for a ' ' blank, if the text renderer supports this. Otherwise fall back to old 'X'. (MK)
%           Implement new optional 4th return argument 'wordbounds'. If caller requests it,
%           we calculate bounding boxes for each individual word in a string, and return a list
%           of bounding boxes. Use case, e.g., to define ROI's around each word for eye tracking
%           experiments. See forum message #21528 (MK).

% Set ptb_drawformattedtext_disableClipping to 1 if text clipping should be disabled:
global ptb_drawformattedtext_disableClipping;
global ptb_drawformattedtext_padthresh;

if isempty(ptb_drawformattedtext_disableClipping)
    % Text clipping on by default:
    ptb_drawformattedtext_disableClipping = 0;
end

if isempty(ptb_drawformattedtext_padthresh)
    % Threshold for skipping of text justification is 33% by default:
    ptb_drawformattedtext_padthresh = 0.333;
end
padthresh = ptb_drawformattedtext_padthresh;

if nargin < 1 || isempty(win)
    error('DrawFormattedText: Windowhandle missing!');
end

if nargin < 2 || isempty(tstring)
    % Empty text string -> Nothing to do, but assign dummy values:
    [nx, ny] = Screen('DrawText', win, '');
    textbounds = [nx, ny, nx, ny];
    wordbounds = textbounds;
    return;
end

if IsOctave
    % char() casts of unicode values > 255 map to zero, because Octave
    % uses UTF-8 encoding for unicode, instead of UTF-32 as Matlab. We
    % take care of this in the code, but necessary casts() also trigger
    % an out-of-range warning in Octave, which we can't selectively disable,
    % as it lacks a unique warning id (duh!). Therefore disable all warnings
    % on Octave and re-enable to previous setting whenever the we exit, and
    % therefore the canary variable reenablewarn goes out of scope:
    warningstate = warning('query');
    warning('off');
    reenablewarn = onCleanup(@() restorewarningstate(warningstate));
end

wordbounds = [];

% Store data class of input string for later use in re-cast ops:
stringclass = class(tstring);

% Default x start position is left border of window:
if nargin < 3 || isempty(sx)
    sx=0;
end

xcenter = 0;
ljustify = 0;
rjustify = 0;
bjustify = 0;
if ischar(sx)
    switch lower(sx)
        case 'center'
            xcenter = 1;
        case 'left'
            ljustify = 1;
        case 'right'
            rjustify = 1;
        case 'wrapat'
            bjustify = 1;
        case 'justifytomax'
            bjustify = 2;
        case 'centerblock'
            bjustify = 3;
    end

    % Set sx to neutral setting:
    sx = 0;
end

% No text wrapping by default:
if nargin < 6 || isempty(wrapat)
    wrapat = 0;
end

if (bjustify  == 1) && ~wrapat
    error('Horizontal justification method ''wrapat'' selected, but required ''wrapat'' parameter missing!');
end

% No horizontal mirroring by default:
if nargin < 7 || isempty(flipHorizontal)
    flipHorizontal = 0;
end

% No vertical mirroring by default:
if nargin < 8 || isempty(flipVertical)
    flipVertical = 0;
end

% No vertical mirroring by default:
if nargin < 9 || isempty(vSpacing)
    vSpacing = 1;
end

if nargin < 10 || isempty(righttoleft)
    righttoleft = 0;
end

% Convert all conventional linefeeds into C-style newlines:
newlinepos = strfind(char(tstring), '\n');

% If '\n' is already encoded as a char(10) as in Octave, then
% there's no need for replacemet.
if char(10) == '\n' %#ok<STCMP>
   newlinepos = [];
end

% Need different encoding for repchar that matches class of input tstring:
if isa(tstring, 'double')
    repchar = 10;
elseif isa(tstring, 'uint8')
    repchar = uint8(10);
else
    repchar = char(10);
end

while ~isempty(newlinepos)
    % Replace first occurrence of '\n' by ASCII or double code 10 aka 'repchar':
    tstring = [ tstring(1:min(newlinepos)-1) repchar tstring(min(newlinepos)+2:end)];
    % Search next occurrence of linefeed (if any) in new expanded string:
    newlinepos = strfind(char(tstring), '\n');
end

% Text wrapping requested?
if wrapat > 0
    % Call WrapString to create a broken up version of the input string
    % that is wrapped around column 'wrapat'
    tstring = WrapString(tstring, wrapat);
end

% Query textsize for implementation of linefeeds:
theight = Screen('TextSize', win) * vSpacing;

% Special case single line text string: Compute true bounding box
% of text string and use its height as text height for linefeeds:
numlines = length(strfind(char(tstring), char(10))) + 1;
if numlines == 1
    theight = RectHeight(Screen('TextBounds', win, tstring));
end

% Default y start position is top of window:
if nargin < 4 || isempty(sy)
    % As we use y-position is baseline, need to shift it down
    % a bit to roughly start text at top of the window:
    sy = Screen('TextSize', win);
end

% Default rectangle for centering/formatting text is the client rectangle
% of the 'win'dow, but usercode can specify arbitrary override as 11'th arg:
if nargin < 11 || isempty(winRect)
    winRect = Screen('Rect', win);
end

winHeight = RectHeight(winRect);

% Need some approximation of the height to baseline. Use height of X,
% as it usually has large/maximal ascenders, but no descenders:
blankbounds = Screen('TextBounds', win, 'X', [], [], 1, righttoleft);

% Compute width of a single blank ' ' white-space, in case we need it.
% If the given text renderer does not support the xAdvance property, we
% fallback to bounding box width of a 'X'. As justification only works
% with monospaced fonts anyway, we can do this substitution with good
% results:
[~, ~, ~, blankwidth] = Screen('TextBounds', win, ' ', [], [], 1, righttoleft);
if blankwidth == 0
    blankwidth = RectWidth(blankbounds);
end
baselineHeight = RectHeight(blankbounds);

if ischar(sy) && strcmpi(sy, 'center')
    % Compute vertical centering:
    bbox = SetRect(0,0,1,numlines * theight);

    % Center box in window:
    [rect,dh,dv] = CenterRect(bbox, winRect); %#ok<ASGLU>

    % Initialize vertical start position sy with vertical offset of
    % centered text box, corrected for baselineHeight because we use
    % yPosIsBaseline = 1:
    sy = dv + baselineHeight;
    yPosIsBaseline = 1;
else
    yPosIsBaseline = 1;
end

if ~yPosIsBaseline
    baselineHeight = 0;
end

% Keep current text color if no one provided:
if nargin < 5 || isempty(color)
    color = [];
end

% Is the OpenGL userspace context for this 'windowPtr' active, as required?
[previouswin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');

% OpenGL rendering for this window active?
if IsOpenGLRendering
    % Yes. We need to disable OpenGL mode for that other window and
    % switch to our window:
    Screen('EndOpenGL', win);
end

% Disable culling/clipping if bounding box is requested as 3rd return
% argument, or if forcefully disabled. Unless clipping is forcefully
% enabled.
disableClip = (ptb_drawformattedtext_disableClipping ~= -1) && ...
              ((ptb_drawformattedtext_disableClipping > 0) || (nargout >= 3));

if bjustify || ljustify
    sx = winRect(RectLeft);
end

% Init cursor position:
xp = sx;
yp = sy;

minx = inf;
miny = inf;
maxx = 0;
maxy = 0;

if bjustify == 1
    % Pad to line width of a line of 'wrapat' X'es:
    padwidth = RectWidth(Screen('TextBounds', win, char(repmat('X', 1, wrapat)), [], [], 1, righttoleft));
end

if bjustify == 2 || bjustify == 3
    % Iterate over whole text string and find widest
    % text line. Use it as reference for padding:
    backuptext = tstring;

    % No clipping allowed in this opmode:
    disableClip = 1;

    % Iterate:
    padwidth = 0;
    while ~isempty(tstring)
        % Find next substring to process:
        crpositions = strfind(char(tstring), char(10));
        if ~isempty(crpositions)
            curstring = tstring(1:min(crpositions)-1);
            tstring = tstring(min(crpositions)+1:end);
        else
            curstring = tstring;
            tstring =[];
        end

        if ~isempty(curstring)
            padwidth = max(padwidth, RectWidth(Screen('TextBounds', win, curstring, [], [], 1, righttoleft)));
        end
    end

    % Restore original string for further processing:
    tstring = backuptext;
end

% Parse string, break it into substrings at line-feeds:
while ~isempty(tstring)
    deltaboxX = 0;
    deltaboxY = 0;

    % Find next substring to process:
    crpositions = strfind(char(tstring), char(10));
    if ~isempty(crpositions)
        curstring = tstring(1:min(crpositions)-1);
        tstring = tstring(min(crpositions)+1:end);
        dolinefeed = 1;
    else
        curstring = tstring;
        tstring =[];
        dolinefeed = 0;
    end

    % tstring contains the remainder of the input string to process in next
    % iteration, curstring is the string we need to draw now.

    % Perform crude clipping against upper and lower window borders for this text snippet.
    % If it is clearly outside the window and would get clipped away by the renderer anyway,
    % we can safe ourselves the trouble of processing it:
    if disableClip || ((yp + theight >= winRect(RectTop)) && (yp - theight <= winRect(RectBottom)))
        % Inside crude clipping area. Need to draw.
        noclip = 1;
    else
        % Skip this text line draw call, as it would be clipped away
        % anyway.
        noclip = 0;
        dolinefeed = 1;
    end

    % Any string to draw?
    if ~isempty(curstring) && noclip
        % Cast curstring back to the class of the original input string, to
        % make sure special unicode encoding (e.g., double()'s) does not
        % get lost for actual drawing:
        curstring = cast(curstring, stringclass);

        % Need bounding box?
        if xcenter || flipHorizontal || flipVertical || rjustify
            % Compute text bounding box for this substring:
            [bbox, refbbox] = Screen('TextBounds', win, curstring, 0, 0, 1, righttoleft);
            deltaboxX = refbbox(RectLeft) - bbox(RectLeft);
            deltaboxY = refbbox(RectTop) - bbox(RectTop);
        end

        % Horizontally centered output required?
        if xcenter
            % Yes. Compute dh, dv position offsets to center it in the center of window.
            [rect,dh] = CenterRect(bbox, winRect); %#ok<ASGLU>
            % Set drawing cursor to horizontal x offset:
            xp = dh - deltaboxX;
        end

        % Horizontally centered output of the block containing the full text string as a whole required?
        if bjustify == 3
            rect = CenterRect([0, 0, padwidth, 1], winRect);
            xp = rect(RectLeft);
        end

        % Right justified (aligned) output required?
        if rjustify
            xp = winRect(RectRight) - RectWidth(bbox);
        end

        if flipHorizontal || flipVertical
            if bjustify
                warning('Text justification to wrapat''th right column border not supported for flipHorizontal or flipVertical text drawing.');
            end

            % Per word bounding boxes too much implementation effort atm.:
            if nargout >= 4
                error('Return of per-word bounding boxes ''wordbounds'' not supported for flipHorizontal or flipVertical text drawing.');
            end

            textbox = OffsetRect(bbox, xp, yp);
            [xc, yc] = RectCenterd(textbox);

            % Make a backup copy of the current transformation matrix for later
            % use/restoration of default state:
            Screen('glPushMatrix', win);

            % Translate origin into the geometric center of text:
            Screen('glTranslate', win, xc, yc, 0);

            % Apple a scaling transform which flips the direction of x-Axis,
            % thereby mirroring the drawn text horizontally:
            if flipVertical
                Screen('glScale', win, 1, -1, 1);
            end

            if flipHorizontal
                Screen('glScale', win, -1, 1, 1);
            end

            % We need to undo the translations...
            Screen('glTranslate', win, -xc, -yc, 0);
            [nx ny] = Screen('DrawText', win, curstring, xp, yp, color, [], yPosIsBaseline, righttoleft);
            Screen('glPopMatrix', win);
        else
            % Block justification (align to left border and a right border at 'wrapat' columns)?
            if ismember(bjustify, [1,2])
                % Calculate required amount of padding in pixels:
                strwidth = padwidth - RectWidth(Screen('TextBounds', win, curstring(~isspace(curstring)), [], [], yPosIsBaseline, righttoleft));
                padpergapneeded = length(find(isspace(curstring)));
                % Padding needed and possible?
                if (padpergapneeded > 0) && (strwidth > 0)
                    % Required padding less than padthresh fraction of total
                    % width? If not we skip justification, as it would lead to
                    % ridiculous looking results:
                    if strwidth < padwidth * padthresh
                        % For each isolated blank in the text line, insert
                        % padpergapneeded pixels of blank space:
                        padpergapneeded = strwidth / padpergapneeded;
                    else
                        padpergapneeded = blankwidth;
                    end
                else
                    padpergapneeded = 0;
                end
            else
                padpergapneeded = 0;
            end

            % Non-zero padpergapneeded for block justification? Or per-word bounding boxes
            % requested? Then we need to render text line word by word:
            if nargout >= 4 || padpergapneeded > 0
                % Render text line word by word, adding padpergapneeded pixels of blank space
                % between consecutive words:
                nx = xp;
                ny = yp;
                [wordup, remstring] = strtok(curstring);
                cxp = xp;
                while ~isempty(wordup)
                    if nargout >= 4
                        % Caller wants per-word bounding boxes returned:
                        [~, wordbound] = Screen('TextBounds', win, wordup, cxp, yp, yPosIsBaseline, righttoleft);
                        wordbounds(end+1, :) = wordbound;
                    end

                    if (padpergapneeded == 0) && ~isempty(remstring)
                        wordup = [wordup ' '];
                    end

                    [nx ny] = Screen('DrawText', win, wordup, cxp, yp, color, [], yPosIsBaseline, righttoleft);
                    if ~isempty(remstring)
                        nx = nx + padpergapneeded;
                        cxp = nx;
                    end
                    [wordup, remstring] = strtok(remstring);
                end
            else
                % Fast path: Draw one full line of text.
                [nx ny] = Screen('DrawText', win, curstring, xp, yp, color, [], yPosIsBaseline, righttoleft);
            end
        end
    else
        % This is an empty substring (pure linefeed). Just update cursor
        % position:
        nx = xp;
        ny = yp;
    end

    % Update bounding box:
    minx = min([minx , xp + deltaboxX, nx]);
    maxx = max([maxx , xp, nx + deltaboxX]);
    miny = min([miny , yp, ny]);
    maxy = max([maxy , yp, ny]);

    % Linefeed to do?
    if dolinefeed
        % Update text drawing cursor to perform carriage return:
        if ~xcenter && ~rjustify
            xp = sx;
        end
        yp = ny + theight;
    else
        % Keep drawing cursor where it is supposed to be:
        xp = nx;
        yp = ny;
    end
    % Done with substring, parse next substring.
end

% Add one line height:
if numlines > 1
    % Make sure to discount vSpacing for the last line
    % to avoid a large void below the drawn text.
    maxy = maxy + theight / vSpacing;
else
    % Special case one-liner without linefeed. Snugly
    % fit the bounding box:
    maxy = maxy + theight;
end

% Create final bounding box:
textbounds = SetRect(minx, miny - baselineHeight, maxx, maxy - baselineHeight);

% Create new cursor position. The cursor is positioned to allow
% to continue to print text directly after the drawn text.
% Basically behaves like printf or fprintf formatting.
nx = xp;
ny = yp;

if nargout >= 4 && isempty(wordbounds)
    wordbounds = textbounds;
end

% Our work is done. If a different window than our target window was
% active, we'll switch back to that window and its state:
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

return;
end

% Restore warning() settings to initial at onCleanup():
function restorewarningstate(warningstate)
    warning(warningstate);
end
