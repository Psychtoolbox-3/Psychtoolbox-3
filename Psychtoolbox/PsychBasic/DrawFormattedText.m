function [nx, ny, textbounds] = DrawFormattedText(win, tstring, sx, sy, color, wrapat, flipHorizontal, flipVertical, vSpacing, righttoleft, winRect)
% [nx, ny, textbounds] = DrawFormattedText(win, tstring [, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical][, vSpacing][, righttoleft][, winRect])
%
% Draws a string of text 'tstring' into Psychtoolbox window 'win'. Allows
% some basic formatting. The text string 'tstring' may contain newline
% characters '\n'. Whenever a newline character '\n' is encountered, a
% linefeed and carriage return is performed, breaking the text string into
% lines. 'sx' defines the left border of the text: If it is left out, text
% starts at x-position zero, otherwise it starts at the specified position
% 'sx'. If sx=='center', then each line of text is horizontally centered in
% the window. If sx=='right', then each line of text is right justified to
% the right border of the target window, or of 'winRect' if provided. 'sy'
% defines the top border of the text. If left out, it starts at the top of
% the window, otherwise it starts at the specified vertical pixel position.
% If sy=='center', then the whole text is vertically centered in the
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
% rectange, in which the text should be centered/placed etc. By default,
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
% styles on some operating systems.
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

% Set ptb_drawformattedtext_disableClipping to 1 if text clipping should be disabled:
global ptb_drawformattedtext_disableClipping;

if isempty(ptb_drawformattedtext_disableClipping)
    % Text clipping on by default:
    ptb_drawformattedtext_disableClipping = 0;
end

if nargin < 1 || isempty(win)
    error('DrawFormattedText: Windowhandle missing!');
end

if nargin < 2 || isempty(tstring)
    % Empty text string -> Nothing to do.
    return;
end

% Store data class of input string for later use in re-cast ops:
stringclass = class(tstring);

% Default x start position is left border of window:
if nargin < 3 || isempty(sx)
    sx=0;
end

xcenter = 0;
rjustify = 0;
if ischar(sx)
    if strcmpi(sx, 'center')
        xcenter = 1;
    end
    
    if strcmpi(sx, 'right')
        rjustify = 1;
    end
    
    % Set sx to neutral setting:
    sx=0;
end

% No text wrapping by default:
if nargin < 6 || isempty(wrapat)
    wrapat = 0;
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
    % Replace first occurence of '\n' by ASCII or double code 10 aka 'repchar':
    tstring = [ tstring(1:min(newlinepos)-1) repchar tstring(min(newlinepos)+2:end)];
    % Search next occurence of linefeed (if any) in new expanded string:
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

% Default y start position is top of window:
if nargin < 4 || isempty(sy)
    sy=0;
end

% Default rectangle for centering/formatting text is the client rectangle
% of the 'win'dow, but usercode can specify arbitrary override as 11'th arg:
if nargin < 11 || isempty(winRect)
    winRect = Screen('Rect', win);
end

winHeight = RectHeight(winRect);

if ischar(sy) && strcmpi(sy, 'center')
    % Compute vertical centering:
    
    % Compute height of text box:
    numlines = length(strfind(char(tstring), char(10))) + 1;
    bbox = SetRect(0,0,1,numlines * theight);
    % Center box in window:
    [rect,dh,dv] = CenterRect(bbox, winRect); %#ok<ASGLU>

    % Initialize vertical start position sy with vertical offset of
    % centered text box:
    sy = dv;
end

% Keep current text color if noone provided:
if nargin < 5 || isempty(color)
    color = [];
end

% Init cursor position:
xp = sx;
yp = sy;

minx = inf;
miny = inf;
maxx = 0;
maxy = 0;

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


% Parse string, break it into substrings at line-feeds:
while ~isempty(tstring)
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

    if IsOSX
        % On OS/X, we enforce a line-break if the unwrapped/unbroken text
        % would exceed 250 characters. The ATSU text renderer of OS/X can't
        % handle more than 250 characters.
        if size(curstring, 2) > 250
            tstring = [curstring(251:end) tstring]; %#ok<AGROW>
            curstring = curstring(1:250);
            dolinefeed = 1;
        end
    end
    
    if IsWin
        % On Windows, a single ampersand & is translated into a control
        % character to enable underlined text. To avoid this and actually
        % draw & symbols in text as & symbols in text, we need to store
        % them as two && symbols. -> Replace all single & by &&.
        if isa(curstring, 'char')
            % Only works with char-acters, not doubles, so we can't do this
            % when string is represented as double-encoded Unicode:
            curstring = strrep(curstring, '&', '&&');
        end
    end
    
    % tstring contains the remainder of the input string to process in next
    % iteration, curstring is the string we need to draw now.

    % Perform crude clipping against upper and lower window borders for
    % this text snippet. If it is clearly outside the window and would get
    % clipped away by the renderer anyway, we can safe ourselves the
    % trouble of processing it:
    if disableClip || ((yp + theight >= 0) && (yp - theight <= winHeight))
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
            bbox=Screen('TextBounds', win, curstring, [], [], [], righttoleft);
        end
        
        % Horizontally centered output required?
        if xcenter
            % Yes. Compute dh, dv position offsets to center it in the center of window.
            [rect,dh] = CenterRect(bbox, winRect); %#ok<ASGLU>
            % Set drawing cursor to horizontal x offset:
            xp = dh;
        end
        
        % Right justified (aligned) output required?
        if rjustify
            xp = winRect(RectRight) - RectWidth(bbox);
        end
        
        if flipHorizontal || flipVertical
            textbox = OffsetRect(bbox, xp, yp);
            [xc, yc] = RectCenter(textbox);

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
            [nx ny] = Screen('DrawText', win, curstring, xp, yp, color, [], [], righttoleft);
            Screen('glPopMatrix', win);
        else
            [nx ny] = Screen('DrawText', win, curstring, xp, yp, color, [], [], righttoleft);
        end
    else
        % This is an empty substring (pure linefeed). Just update cursor
        % position:
        nx = xp;
        ny = yp;
    end

    % Update bounding box:
    minx = min([minx , xp, nx]);
    maxx = max([maxx , xp, nx]);
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
maxy = maxy + theight;

% Create final bounding box:
textbounds = SetRect(minx, miny, maxx, maxy);

% Create new cursor position. The cursor is positioned to allow
% to continue to print text directly after the drawn text.
% Basically behaves like printf or fprintf formatting.
nx = xp;
ny = yp;

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
