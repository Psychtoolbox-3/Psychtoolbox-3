function [dpi, SP]=MeasureDpi(theScreen, CC)
% dpi=MeasureDpi([theScreen], [CC]);
% Helps the user to accurately measure the screen's dots per inch. 
% INPUT:
% theScreen -> Screen parameter; 
% CC-> Creditcard parameter, 1-> Yes. Other-> no
% 
% OUTPUT:
% dpi-> Dots per inch.
% SP -> Screen presets withdrawn from the Monitor drivers(EDID). It will provide
% you with the physical size Width and Height, the screen resolution, the
% DotPitch and dpi automatically. Note: This feature has not been tested in
% Linux or MacOS, and it won't work with old monitors.
%
% Denis Pelli

% 5/28/96 dgp Updated to use new GetMouse.
% 3/20/97 dgp Updated.
% 4/2/97  dgp FlushEvents.
% 4/26/97 dhb Got rid of call to disp(7), which was writing '7' to the
%             command window and not producing a beep.
% 8/16/97 dgp Changed "text" to "theText" to avoid conflict with TEXT function.
% 4/06/02 awi Check all elements of the new multi-element button vector
%             returned by GetMouse on Windows.
%             Replaced Chicago font with Arial because it's available on
%             both Mac and Windows
% 11/6/06 dgp Updated from PTB-2 to PTB-3.
% 01/31/16 mk Fix some wrong assumptions for PTB-3. Get rid of GetChar.
% 8/13/19 mgg Mofify to use it with a credit card (common object of a
%             known size).
% 9/10/19 mgg Optional output: Withdraw the presets from the screen
%             inc. the screen size, resolution, dot pitch and distance independent dpi.

if nargin>2 || nargout>2
    error([datestr(now),':>> Error wrong call to function, use it as: dpi=MeasureDpi(Screen)', newline, 'or [dpi, SP]=MeasureDpi(Screen, CC)']);
end
if ~exist(theScreen) || isempty(theScreen)
    theScreen=max(Screen('Screens')); % Use the second screen, most of the cases.
    disp([datestr(now),':>> No screen value has been provided so I assume ',...
        'that either you would like to use the secondary screen, or the main one if you only have one...']);
end
AssertOpenGL;

try
    if nargin>1 && isequal(CC,1) % Using a credit card
        unitInches=1;
        units='inches'; % e.g. distance in inches
        unit='inch';    % e.g. 5 inch object
        objectInches= 8.56*(1/2.54); % ID-1 format, ISO/IEC 7810 credit card (85.60 × 53.98 mm).
        [unitInches, units, unit] = EnglishOrSI();
    else
        disp('Please find an object of known width to hold against the display.');
        disp('E.g.  a ruler or an 8.5-inch page.');
        [unitInches, units, unit] = EnglishOrSI();
        objectInches=input(sprintf('How wide is your object, in %s? ',units))*unitInches;
        disp('A small correction will be made for your viewing distance and the thickness of');
        disp('the screen''s clear front plate, which separates your object from the screen''s');
        disp('light emitting surface.');
    end
    isLCD = input('Is your display a flat panel, instead of a CRT? [y/n]: ', 's');
    if isLCD == 'y'
        thicknessInches=0.1;
        fprintf('I assume that the display is a flat panel display, \n');
        fprintf('with a thin (%.1f inch) clear front plate.\n',thicknessInches);
    else
        thicknessInches=0.25;
        fprintf('I assume that the display is a CRT, \n');
        fprintf('with a thick (%.1f inch) clear front plate.\n',thicknessInches);
    end
    distanceInches=input(sprintf('What is your viewing distance, roughly, in %s? ',units))*unitInches;
    [window,screenRect]=Screen('OpenWindow',theScreen);
    white=WhiteIndex(window);
    black=BlackIndex(window);
    
    % This might work if you have a modern monitor with reasonable updated
    % drivers.
    [SP.ScreenWidthOS, SP.ScreenHeightOS] = Screen('DisplaySize', theScreen);
    [SP.Xpixels, SP.Ypixels] = Screen('WindowSize', w);  
    SP.DotPitchX= (SP.ScreenWidthOS/SP.Xpixels);
    SP.DotPitchY=(SP.ScreenHeightOS/SP.Ypixels);
    SP.DotPitchDiag= (SP.DotPitchX + SP.DotPitchY) /2;
    SP.dpi=floor((sqrt((SP.Xpixels)^2+(SP.Ypixels)^2))/(sqrt((SP.ScreenWidthOS)^2+(SP.ScreenHeightOS)^2)/25.4));
    
    % Instructions
    if nargin>1 && isequal(CC,1) % Using a credit card
        s=sprintf('Hold your credit card against the display.');
    else
        s=sprintf('Hold your %.1f-%s-wide object against the display.',objectInches/unitInches,unit);
    end
    theText={s,'Press, drag, and release the mouse to draw a bar'...
        ,'that matches the width of your object. Use one eye.'};
    Screen('TextFont',window,'Arial');
    s=24;
    Screen('TextSize',window,s);
    textLeading=s+8;
    textRect=Screen('TextBounds',window,theText{1});
    textRect(4)=length(theText)*textLeading;
    textRect=CenterRect(textRect,screenRect);
    textRect=AlignRect(textRect,screenRect,RectTop);
    textRect(RectRight)=screenRect(RectRight);
    dragText=theText;
    dragTextRect=textRect;

    % Animate
    % Track horizontal mouse position to draw a bar of variable width.
    for i=1:length(dragText)
        Screen('DrawText',window,dragText{i},dragTextRect(RectLeft),dragTextRect(RectTop)+textLeading*i,black);
    end
    barRect=CenterRect(SetRect(0,0,RectWidth(screenRect),20),screenRect);
    fullBarRect=barRect;
    top=RectTop;
    bottom=RectBottom;
    left=RectLeft;
    right=RectRight;
    Screen('FillRect',window,white,barRect);
    Screen('FrameRect',window,black,barRect);
    Screen('Flip',window);
    oldButton=0;
    while 1
        [x,~,button]=GetMouse;
        if any(button)
            if ~oldButton
                origin=x;
                barRect(left)=origin;
                barRect(right)=origin;
            else
                if x<origin
                    barRect(left)=x;
                    barRect(right)=origin;
                else
                    barRect(left)=origin;
                    barRect(right)=x;
                end
            end
            if ~IsEmptyRect(barRect)
                Screen('FillRect',window,black,barRect);
            end
            backgroundRect=barRect;
            backgroundRect(left)=screenRect(left);
            backgroundRect(right)=barRect(left);
            if ~IsEmptyRect(backgroundRect)
                Screen('FillRect',window,white,backgroundRect);
            end
            backgroundRect(left)=barRect(right);
            backgroundRect(right)=screenRect(right);
            if ~IsEmptyRect(backgroundRect)
                Screen('FillRect',window,white,backgroundRect);
            end
            for i=1:length(dragText)
                Screen('DrawText',window,dragText{i},dragTextRect(RectLeft),dragTextRect(RectTop)+textLeading*i,black);
            end
            Screen('FrameRect',window,black,fullBarRect);
            Screen('Flip',window);
            if ~IsEmptyRect(barRect)
                Screen('FillRect',window,black,barRect);
            end
        else
            if oldButton
                objectPix=RectWidth(barRect);
                dpi=objectPix/objectInches;
                dpi=dpi*distanceInches/(distanceInches+thicknessInches);
                clear theText
                if inches
                    theText{1}=sprintf('%.0f dots per inch.',dpi);
                else
                    theText{1}=sprintf('%.0f dots per inch. (%.0f dots/cm.)',dpi,dpi/2.54);
                end
                theText{2}='Click once to do it again; twice if you''re done.';
                textRect=Screen('TextBounds',window,theText{2});
                textRect(4)=length(theText)*textLeading;
                textRect=CenterRect(textRect,screenRect);
                textRect=AlignRect(textRect,screenRect,RectTop);
                Screen('FillRect',window,white,InsetRect(dragTextRect,0,round(-0.3*textLeading)));
                for i=1:length(theText)
                    Screen('DrawText',window,theText{i},textRect(RectLeft),textRect(RectTop)+textLeading*i,black);
                end
                Screen('Flip',window);
                i=GetClicks;
                if ~IsEmptyRect(barRect)
                    Screen('FillRect',window,white,barRect);
                end
                Screen('FillRect',window,white,InsetRect(textRect,0,round(-0.3*RectHeight(textRect))));
                if i>1
                    break
                end
                for i=1:length(dragText)
                    Screen('DrawText',window,dragText{i},dragTextRect(RectLeft),dragTextRect(RectTop)+textLeading*i,black);
                end
                Screen('FrameRect',window,black,fullBarRect);
                Screen('Flip',window);
            end
        end
        oldButton=any(button);
    end
    Screen('Close',window);
catch
    sca;
    psychrethrow(psychlasterror);
end

function [unitInches, units, unit] = EnglishOrSI()
% Ask user whether to use the International System of Units or the English
% / American System(inches).
inches=input('Do you prefer inches (1) or cm (0)? ');
        if inches
            unitInches=1;
            units='inches'; % e.g. distance in inches
            unit='inch';    % e.g. 5 inch object
        else
            unitInches=1/2.54;
            units='cm';
            unit='cm';
        end
