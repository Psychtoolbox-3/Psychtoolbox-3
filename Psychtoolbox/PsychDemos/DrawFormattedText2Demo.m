% DrawFormattedText2Demo
%
% Draws formatted text, shows how to align and center text vertically
% and/or horizontally, how line-breaks are introduced, how to transform
% text.
%
% Press any key to cycle through different demo displays.
%
% see also: PsychDemos, Screen DrawText?, DrawSomeTextDemo, DrawFormattedText2Demo

% 15-Sep-2016    dcn    Wrote it.
% 07-May-2017    mk     Demo per word bounding boxes.


try
    PsychDefaultSetup(1);   % check PTB is functioning and do default setup
    
    % Choosing the display with the highest display number is
    % a best guess about where you want the stimulus displayed.
    screens=Screen('Screens');
    screenNumber=max(screens);

    % Open window with default settings:
    PsychImaging('PrepareConfiguration');
    % make screen partially transparent so we can debug.
    % PsychDebugWindowConfiguration;
    [w,wrect]=PsychImaging('OpenWindow', screenNumber, 128);
    
    % This preference setting selects the high quality text renderer 
    % each operating system: It is not really needed, as the high quality
    % renderer is the default on all operating systems, so this is more of
    % a "better safe than sorry" setting. DrawFormattedText2 _only_ works
    % with this renderer. So if it fails to load (see the warnings on the
    % command line), the output of this demo will look horrible. See
    % help DrawTextPlugin for more information.
    Screen('Preference', 'TextRenderer', 1);
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);   % has no effect on text rendering, used below for semi-transparent bounding boxes
    Screen('Preference', 'TextAntiAliasing', 2);

    
    % aligning of text to a box
    [scrX,scrY] = RectCenter(wrect);
    rect = CenterRectOnPoint([0 0 300 300],scrX,scrY);
    Screen('FrameRect', w, 0, rect);
    DrawFormattedText2('top\ncenter','win',w,'sx','center','sy',scrY-150,'xalign','center','yalign','bottom','xlayout','center');
    DrawFormattedText2('top\nleft<color=ff0000>1','win',w,'sx',scrX-150,'sy',scrY-150,'xalign','right','yalign','bottom','xlayout','right');
    DrawFormattedText2('top\nleft<color=0.,1.,0.>2','win',w,'sx',scrX-150,'sy',scrY-150,'xalign', 'left','yalign','bottom');
    DrawFormattedText2('center\n<b>inside','win',w,'sx','center','sy','center','xalign','center','yalign','center','xlayout','center','winRect',rect);
    DrawFormattedText2('bottom\n<i>left' ,'win',w,'sx',scrX-150,'sy',scrY+150,'xalign','right','yalign','top','xlayout','right');
    DrawFormattedText2('bottom\n<u>right','win',w,'sx',scrX+150,'sy',scrY+150,'xalign', 'left','yalign','top');

    Screen('Flip',w);
    KbStrokeWait;
    
    % unicode text, see DrawHighQualityUnicodeTextDemo for more info about
    % how to draw unicode text
    unicodetext = [26085, 26412, 35486, 60, 99, 111, 108, 111, 114, 62, ...
        12391, 12354, 12426, 12364, 12392, 12358, 12372, 12374, 12356, ...
        12414, 12375, 12383, 12290, 13, 10];
    [~,~,~,~,wbounds] = DrawFormattedText2([double('<size=40><font=-:lang=ja><color=ff0000>') unicodetext],'win',w,'sx','center','sy','center','xalign','center','yalign','center');
    % Screen('FrameRect', w, [255 255 0], wbounds.');
    
    Screen('Flip',w);
    KbStrokeWait;

    % draw few lines of text, with varying formatting, and draw it again
    % from cache with some transformation
    [~,~,bbox,cache,wbounds]=DrawFormattedText2('<font=Courier New><size=27>test\n<font=Times New Roman>scr<font><font>een<font> is\n<b><size=50>UGLY\n<size=12><b><u><i>Isn''t it?','win',w,'sx','center','sy','center','xalign','center','yalign','center');
    Screen('FrameRect', w, [255 0 0 100], bbox);
    % captbbox = bbox;
    Screen('FrameRect', w, [255 255 0], wbounds.');
    % calling with the cache and only sx and/or sy, translated the text by
    % (sx,sy)
    [~,~,bbox,~,wbounds]=DrawFormattedText2(cache,'sx', 300);
    Screen('FrameRect', w, [0 255 0], bbox);
    Screen('FrameRect', w, [255 255 0], wbounds.');
    [~,~,bbox,~,wbounds]=DrawFormattedText2(cache,'sx',-300,'transform',{'flip',1});
    Screen('FrameRect', w, [0 0 255], bbox);
    Screen('FrameRect', w, [255 255 0], wbounds.');
    [~,~,bbox,~,wbounds]=DrawFormattedText2(cache,'sy', 300,'transform',{'rotate',90});
    Screen('FrameRect', w, [0 255 255], bbox);
    Screen('FrameRect', w, [255 255 0], wbounds.');
    [~,~,bbox,~,wbounds]=DrawFormattedText2(cache,'sy',-300,'transform',{'scale',[2 1]});
    Screen('FrameRect', w, [255 0 255], bbox);
    Screen('FrameRect', w, [255 255 0], wbounds.');
    % when called with more alignment inputs, the bounding box of the text
    % is repositioned:
    rect = [100 100 350 250];
    Screen('FrameRect', w, 0, rect);
    % put text in top left of this box
    DrawFormattedText2(cache,'sx',0,'sy',0,'xalign','left','yalign','top','winRect',rect);
    % make another box and put text next to bottom right of this box,
    % vspace differenty
    rect = wrect([3:4 3:4])-[650 350 300 100];
    Screen('FrameRect', w, [255 255 0], rect);
    DrawFormattedText2('simple\ntest\ntext','win',w,'sx','right','sy','bottom','xalign','left','yalign','bottom','winRect',rect,'vSpacing',1.5);
    
    Screen('Flip',w);
    % capt = Screen('GetImage', w, GrowRect(captbbox,20,20));
    KbStrokeWait;

    % draw some text, then draw it again in exact same location but rotated
    % 180 degrees. should have exact same bounding box
    [~,~,bbox1,c]=DrawFormattedText2('<size=40>t<i>e<size>x<i>t&\n<size=120>ajX\n<size=60>t<font=comic sans>estddd<color=ff0000>ddda<i>d','win',w,'sx','center','sy',200,'xalign','left','baseColor',0,'cacheMode',2);  % just testing cachemode 2. not recommended to be used for performance reasons. Only use if you want to manually change content of the cache, though note the cache format may change without warning
    [~,~,bbox2]  =DrawFormattedText2(c,'transform',{'rotate',180});
    % this is equivalent to:
    % [~,~,bbox2]  =DrawFormattedText2(c,'transform',{'flip',3});
    Screen('FrameRect',w,[255 0 0 128],bbox1,3);
    Screen('FrameRect',w,[0 255 0 128],bbox2,3);
    fprintf('Both bounding boxes exactly the same? %d\n',isequal(bbox1,bbox2));

    [~,~,bbox1,c]=DrawFormattedText2('<size=40>t<i>e<size>x<i>t&\n<size=120>ajX\n<size=60>t<font=comic sans>estddd<color=ff0000>ddda<i>d','win',w,'sx','center','sy',200,'xalign','left','baseColor',0);
    [~,~,bbox2,~,wbounds]  =DrawFormattedText2(c,'transform',{'rotate',180});
    Screen('FrameRect',w,[255 0 0 128],bbox1,3);
    Screen('FrameRect',w,[0 255 0 128],bbox2,3);
    Screen('FrameRect', w, [255 255 0 128], GrowRect(wbounds, 2, 2)');

    Screen('Flip',w);
    KbStrokeWait;
    
    
    % Read a text file:
    fd = fopen([PsychtoolboxRoot 'Contents.m'], 'rt');
    if fd==-1
        error('Could not open Contents.m file in PTB root folder!');
    end
    t = fread(fd,inf,'*char').';
    fclose(fd);
    % flop it all on the screen
    [~,~,~,~,wbounds] = DrawFormattedText2(t,'win',w,'sx','center','xalign','center','sy',150,'baseColor',0);
    Screen('FrameRect', w, [255 255 0 128], GrowRect(wbounds, 2, 2)');

    Screen('Flip',w);
    KbStrokeWait;
    
    % now scroll it over the screen. first pregenerate whole text at middle
    % of screen, and get cache only
    [~,~,bbox,cache]=DrawFormattedText2(t,'win',w,'sx','center','xalign','center','sy','center','baseColor',0,'cacheOnly',true);
    % scroll over screen in 5 seconds
    frate = Screen('FrameRate',screenNumber);
    if frate == 0
        % Deal with Apples trainwreck not reporting framerate on Mac
        % builtin displays, but zero:
        frate = 60;
    end

    textHeight = RectHeight(bbox);
    step = textHeight/frate/5;
    
    tp=zeros(1, frate*5);
    for p=1:frate*5
        [~,~,bbox]=DrawFormattedText2(cache,'sy',-p*step);
        Screen('FrameRect', w, [255 0 0 128], bbox,2);

        tp(p) = Screen('Flip', w);

        if KbCheck
            break;
        end
    end
    tp = tp(1:p);
    fprintf('\nAverage redraw duration for scrolling in msecs: %f\n', 1000 * mean(diff(tp)));
    
    % End of demo, close window:
    sca;
catch me
    % This "catch" section executes in case of an error in the "try"
    % section above. Importantly, it closes the onscreen window if it's
    % open.
    sca;
    rethrow(me)
end

% image(capt)
