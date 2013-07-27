% DrawHighQualityUnicodeTextDemo
%
% This demo shows how to draw high-quality, anti-aliased text, and some
% japanese text encoded in Unicode.
%
% MacOS/X has text fonts with support for japanese characters preinstalled,
% so this should just work out of the box. Have a look at the code of the
% demo on how to select a suitable font, and how to read unicode text from
% the filesystem. The (commented out) reading code for SHIFT_JIS -> Unicode
% conversion would only work on recent Matlab releases. Older releases need
% different approaches - also different wrt. PowerPC vs. IntelMac. That's
% why we hard-coded the text in this demo -- Want to have it working even
% on old Matlab 6...
% 
% On MS-Windows you will need to install the special east-asian font
% support kit in order to be able to draw japanese text.
% How to do this? See...
%
% http://www.coscom.co.jp/japanesefont/index.html
%
% with some more info here:
% http://www.alanwood.net/unicode/fonts_windows.html#japanese
%
% After that, text drawing seems to "just work" with our "Courier New"
% font. If you don't install the font pack, you'll just see funny little
% squares instead of nice japanese characters...
%
% On MacOS/X you have the choice between three different text
% rendering/layouting methods, as selectable by the Screen preference
% setting 'TextRenderer': The default is 1, which selects Apple's ATSU text
% renderer. A setting of 0 would also select ATSU but would use a different
% method for layout of text. Each method has its own weaknesses in that it
% will have trouble drawing text in some special fonts correctly, so you
% may have to choose the method based on the font to get best results.
% Typical standard fonts are handled correctly by both methods. A setting
% of 2 will use the Linux text rendering plugin (see "help DrawTextPlugin")
% which in our experience handles all fonts well and is of superior
% performance to the OSX ATSU renderer. However, given that historically
% ATSU was used on OSX for many years and each layout and rendering method
% will give slightly different text appearance, for reason of consistency
% we leave the default setting for OSX on 1, instead of using the superior
% setting 2.
%

% 11/26/07  mk      Wrote it. Derived from Allens DrawSomeTextDemo.

try
    % Choosing the display with the highest display number is
    % a best guess about where you want the stimulus displayed.
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    % These preference setting selects the high quality text renderer on
    % each operating system: It is not really needed, as the high quality
    % renderer is the default on all operating systems, so this is more of
    % a "better safe than sorry" setting.
    Screen('Preference', 'TextRenderer', 1);

    % This command uncoditionally enables text anti-aliasing for high
    % quality text. It is not strictly needed here, because the default
    % setting is to let the operating system decide what to use - which is
    % usually anti-aliased rendering. This here just to demonstrate the
    % switch. On WindowsXP or Vista, there also exists a setting of 2 for
    % especially hiqh quality anti-aliasing. However, i couldn't ever see
    % any perceptible difference in quality...
    Screen('Preference', 'TextAntiAliasing', 1);

    % This setting disables user defined alpha-blending for text - not
    % strictly needed, alpha blending is disabled by default.
    Screen('Preference', 'TextAlphaBlending', 0);

    % We want the y-position of the text cursor to define the vertical
    % position of the baseline of the text, as opposed to defining the top
    % of the bounding box of the text. This command enables that behaviour
    % by default. However, the Screen('DrawText') command provides an
    % optional flag to override the global default on a case by case basis:
    Screen('Preference', 'DefaultTextYPositionIsBaseline', 1);
    
    % Open an onscreen window, fullscreen with default 50% gray background:
    w=Screen('OpenWindow', screenNumber, 128);

    % Draw some yellow oval, just to make the scene more interesting...
    Screen('FillOval', w, [255 255 0], [0 0 400 400]);

    % Select 'Courier New' as font, choose a text size of 48pts and a
    % default text style:
    Screen('TextFont',w, 'Courier New');
    Screen('TextSize',w, 48);
    Screen('TextStyle', w, 0);
    
    % The standard 'Hello World' message drawn with the low-level
    % Screen('DrawText') command. We start at x=300 and y=100. We draw in
    % Red+Blue color:
    Screen('DrawText', w, double('Hello'), 300, 100, [255, 0, 255]);

    % Add to the text, starting at the last text cursor position. This
    % should append the text to the previously drawn text...
    Screen('DrawText', w, double('World!'));
    
    % Now for some Unicode text rendering... 
    % The following array of double values encodes some japanese text in
    % UTF-16 Unicode. Unicode text must be passed to Screens text drawing
    % functions as a row-vector of double type, i.e., a numeric row vector,
    % or - if it isn't one already - converted via the double() operator
    % etc... Standard text strings are always interpreted as standard 8 bit
    % ASCII text.
    %
    % According to the donor of this text snippet, it says "Thank you in
    % japanese" -- let's hope he didn't lie to me ;-)
    %
    japanesetext = [26085, 26412, 35486, 12391, 12354, 12426, 12364, 12392, 12358, 12372, 12374, 12356, ...
                    12414, 12375, 12383, 12290, 13, 10];

    % The text above is hard-coded. The following *disabled* snippet of code
    % would have done the same, reading the text from a text file, encoded
    % in so called Shift_JIS encoding, then converting it from Shift_JIS
    % into Unicode UTF-16, then into a double matrix. However, this only
    % works on recent Matlab releases, e.g, R2007a (V7.4) and later. In
    % order to make this demo workable on older Matlabs, we just show you
    % the commented code:
    if 0
        fid = fopen([PsychtoolboxRoot 'PsychDemos/japanese_shiftJIS.txt'], 'r', 'n','Shift_JIS');
        japanesetext = native2unicode(fread(fid),'Shift_JIS'); %#ok<N2UNI>
        fclose(fid);
        disp(japanesetext);
        japanesetext = double(transpose(japanesetext));
    end
    
    % On MS-Windows you need to install the special east-asian font
    % support kit in order to be able to draw japanese text.
    % How to do this? See...
    % http://www.coscom.co.jp/japanesefont/index.html
    %
    % with some more info here:
    % http://www.alanwood.net/unicode/fonts_windows.html#japanese
    %
    % After that, text drawing seems to "just work" with our "Courier New"
    % font as selected above.
    
    % Under OS/X...
    if IsOSX
        % ... you must select a font that supports japanese characters...
        
        % This is a bit clumsy: Normally we could just select a font by
        % name. Unfortunately the names of japanese fonts are encoded in
        % japanese unicode instead of plain ASCII text, and the
        % Screen('TextFont') command can't handle this :-(
        %
        % Solution: We use the FontInfo() command to query the system about
        % all installed fonts, then use a loop with string matching to find
        % the numerical 'idx' index of the 'Hiragino Mincho Pro W3' font we
        % wanna use for japanese text. Then we pass the 'idx' index number
        % to Screen('TextFont') and everything is fine...
        allFonts = FontInfo('Fonts');
        foundfont = 0;
        for idx = 1:length(allFonts)
            if strcmpi(allFonts(idx).name, 'Hiragino Mincho Pro W3')
                foundfont = 1;
                break;
            end
        end
        
        if ~foundfont
            error('Could not find wanted japanese font on OS/X !');
        end
        Screen('TextFont', w, allFonts(idx).number);        
    end

    if IsLinux
	% On Linux, we can also auto-select fonts by their supported languages,
	% e.g., we simply require a font with...
	if 1
		% ... support for the 'ja'panese language, whatever fits best:
		Screen('TextFont', w, '-:lang=ja');
	else
		% ... support for the 'he'brew language, whatever fits best:
		Screen('TextFont', w, '-:lang=he');
		% Of course we also need to supply a text string with some
		% hebrew characters (unicode code points) then:
		japanesetext = 1488:1514;
	end
	% ... this would also work on OS/X if 'TextRenderer', type 2 is selected ...
    end

    % Let's draw the text once with the low-level Screen command at
    % location (20, 300) in color black (==0) ...
    y = 200;
    Screen('DrawText', w, japanesetext, 20, y, 0);
    
    % Draw some green line at the top of the letters and at the baseline of
    % the text, just to show how nicely 'DrawText' now obeys the text size
    % settings of 48 pts...
    Screen('DrawLine', w, [0 255 0], 100, y - 48, 1400, y - 48);
    Screen('DrawLine', w, [0 255 0], 100, y, 1400, y);

    % ...and once centered with the convenient high level DrawFormattedText
    % command...
    DrawFormattedText(w, japanesetext, 'center', 'center');

    % Tell user how to exit the demo, this time in the font 'Times' at a
    % size of 86pts, in red color:
    Screen('TextFont',w, 'Times');
    Screen('TextSize',w, 86);
    Screen('DrawText', w, 'Hit any key to continue.', 100, 600, [255, 0, 0]);

    % Show text slide:
    Screen('Flip',w);

    % Wait for keypress, then for key release:
    KbStrokeWait;
    
    % Some funny little animation loop. Text scrolls down from top to
    % bottom of screen - or until key press...
    vbl=Screen('Flip',w);
    tstart=vbl;
    count=-100;
    while ~KbCheck && count < 1200
        % Draw text baseline in green...
        Screen('DrawLine', w, [0 255 0], 100, count, 1400, count);
        % Compute texts bounding box...
        [normRect realRect] = Screen('TextBounds', w, 'Hit any key to exit.', 100, count);
        % Draw the text...
        Screen('DrawText', w, 'Hit any key to exit.', 100, count, [255, 0, 0, 255]);
        % Visualize its bounding box:
        Screen('FrameRect', w, [255 0 255], realRect);

        % Update count and show frame...
        count=count+1;
        Screen('Flip',w);
    end
    
    Screen('Flip', w);

    % Some nice good bye message in blue and at 24 pts text size:
    Screen('TextSize',w, 24);
    DrawFormattedText(w, 'Hit any key to exit the demo. Bye!', 'center', 'center', [0 0 255]);
    Screen('Flip', w);
    
    % Wait for keypress, then flip a last time:
    KbStrokeWait;
    Screen('Flip', w);
    
    Screen('Preference', 'DefaultTextYPositionIsBaseline', 0);

    % Close the screen, we're done...
    Screen('CloseAll');

catch
    % This "catch" section executes in case of an error in the "try" section
    % above.  Importantly, it closes the onscreen window if it's open.
    Screen('Preference', 'DefaultTextYPositionIsBaseline', 0);
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
