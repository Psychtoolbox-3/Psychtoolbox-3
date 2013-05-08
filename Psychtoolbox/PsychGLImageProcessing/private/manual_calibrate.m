function scal = manual_calibrate(scal)
% function outdots = manual_calibrate(screenNum, [dots])
%
% View all dots, calibrate manually.  This function can be used either
% to update the global variables XCALIBDOTS YCALIBDOTS, or to modify the
% locations of a specific subset of dots.  The first is what happens if
% the function is called without the "dots" argument, and the second is
% what
% happens if it is called with this argument.
%
% Input:
%   screenNum   scalar, equal to 1 (left screen) or 2 (right screen)
%   dots        Nx2 array of floating point numbers, [x y] These are
%   dotcenters in screen
%               pixel coordinates. This argument is optional. If omitted the function
%				    modifies XCALIBDOTS and YCALIBDOTS.
% Output:
%   outdots:    same structure as dots.

%
% TODO:
% Clean this puppy up!  Way too big...too many if/else blocks... the
% logic is getting too hard to follow!  Functions, sub-functions...
% something!
%
% 2006-02-13 - cburns
%   - Add feature to select and adjust a group of dots using the mouse.
%   - Changed drawing to a color array so we only render to the Haploscope
%   monitor we are calibrating.  This provide an additional clue to the
%   calibrator in case they forgot to select which monitor to calibrate on.
%
% 2006-01-30 - Chris Burns
% Beginning a re-write of the calibration code on the Bankslab Haploscope
% First pass - find and re-organize files.  Learn current process.
%
% 6/13/01 Jimmy Sastra: Updated for Windows 2000 and Psychophysics Toolbox with improved
%    functionality and use of analaglyphs to display the dots.  Based on original Mac function
%    of the same name by Ben Backus at UC Berkeley (1995-97).

% all these variables are passed on from other function


global gMirrorView


% Functionality for the Auxiliary Function
% If an auxiliary function is defined, we'll call it on each render pass.
% This builds in a functionality to query external devices and provide
% feedback during calibration.  Used on the VisionStation for reading in
% values from the 3DM Magnetometer during calibration.
% The Auxiliary function will return a string that will be displayed during
% calibration.  Of the form:
%           function msg = AuxiliaryFunction
msgAuxText = [];

% Specify text color for user interface feedback
% We want to keep this dim so it's not too bright when using it in a dark
% room during calibration.
iTextPos    = [(scal.rect(3) * 0.3) (scal.rect(4) * 0.5)];
iTextScale  = 24;

% Use yellow by default:
Screen('TextColor', scal.windowPtr, [255 255 0]);

fSoftBeepVolume = 0.1;

% here it is decided what dots to show on the display.... If there are dots
% selected with bvlSelectFitPts, then we will display the selected dots.
% Otherwise, we will display all dots
if isempty(scal.SELECTXCALIBDOTS), 
    dots = [scal.XCALIBDOTS(:)  scal.YCALIBDOTS(:)]; % all dots
else
    dots = [scal.SELECTXCALIBDOTS(:) scal.SELECTYCALIBDOTS(:)]; % or only dot selection
end

nDots = length(dots(:,1));

% Constants
xmax = scal.rect(3)/2 - .5;	  % x coord for rightmost points (origin at center screen)
ymax = scal.rect(4)/2 - .5;	  % (Same values as xmid, ymid, but new name for new use)
xmid = xmax;           % xmid and ymid are in screen coords now, origin at corner..
ymid = ymax;

% Define colors that are used in the drawing functions
% NOTE: If we don't include the alpha value PTB will sometimes crash on
% calls to DrawDots.  I suspect it's an error when i trys to draw over a
% pixel that has an existing value and there's no alpha set so some
% comparision bombs.  Just a guess though. - cburns 2007-05-14
% Update(2007-05-15)  Still not sure on the exact cause of the problem, but
% it appears to be the result of multiple calls to DrawDots called in
% successsion with arrays of dots.
%               [red green blue alpha]
colorBlack      = [0    0   0   255];
colorWhite      = [255  255 255 255];
colorRed        = [255  0   0   255];
colorBlue       = [0    0   255 255];
colorGreen      = [0    255 0   255];
colorCyan       = [0    255 255 255];
colorMagenta    = [255  0   255 255];
colorYellow     = [255  255 0   255];

% Define dot parameters; bv: this should be created earlier in the code
dotDiam         = 8;                 % Diameter of dots in pixels
selectedDotDiam = 2 * dotDiam;       % Diameter of selected dots, in pixels
dotDiamArray    = [];       % Array containing the dot diameter for each dot.
% This allows us to keep the dotDiam and
% selectedDotDiam as constants and use this
% array to update the diameter or selected dots
% on each edit.
dotColor        = colorWhite;
dotStyle        = 2;        % 0 = square pixels, 1 = circles, 2 = circles with antialiasing.


% Enable alpha blending with proper blend-function. We need it
% for drawing of smoothed points:
Screen('BlendFunction', scal.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Adjustment Delta is the offset applied to global or local dots
fAdjustDelta = [0.1     0.5     1   5   10];
iAdjustDeltaIndex = 1;
fDelta = fAdjustDelta(iAdjustDeltaIndex);

% List of different types of adjustments
strAdjustTypes = cell(1, 5);
strAdjustTypes{1} = 'Translation';
strAdjustTypes{2} = 'Expansion or Contraction';
strAdjustTypes{3} = 'Shear';
strAdjustTypes{4} = 'Trapezoid';
strAdjustTypes{5} = 'Pincushion or Barrel';
iAdjustTypeIndex = 1;

% indices of selected dots, those inside the rectangle
indicesSelectedDots = [];

% Mouse position variables, used to draw a rectangle when selecting a group of dots
xMouse      = 0;    % current position of mouse
yMouse      = 0;
xMouseStart = 0;    % starting mouse position for draw rectangle
yMouseStart = 0;

% Mouse Input Timers
fLastMiddleMouseTime = 0;
fLastRightMouseTime = 0;

fTimeDelayForFeedback = 0.500;  % Millisecond timestamp used for feedback on
% mouse input.  How long we wait to cycle
% through choices when use holds down mouse
% buttons.

% It would be ideal if we had a callback mechanism to handle the mouse
% input, particularly in dealing with the buttons.  Since we have no way of
% receiving a "mouse button up" message, we have to poll the mouse, set a
% flag when a button is down and then release that flag when the button is
% up.  This flag is for the left mouse button used when selecting a group
% of dots.
bLeftMouseDown      = 0;
% Similarly for the other mouse buttons
bMiddleMouseDown    = 0;
bRightMouseDown     = 0;

% Draw our own cursor during dot selection mode
% Used along with the bvlDrawCursor routine.
% bShowCursor = 0;    % flag specifying when to show cursor

% adjustMode 1 will modify all dots, adjustMode 2 will allows user to select single dots.
% adjustMode 3 allows user to translate a single dot. When adjustMode = 0 function exits.
adjustMode = 1;

fShowHelptext = 1;

Screen('TextSize', scal.windowPtr, 24);

% bv180108 commented out:
% screen closeall

% Put the entire loop in a try-catch block so we can clean up in case it
% crashes
% display dots and allow for the adjustments
% try
while (adjustMode ~= 0)
    % Translate dots so that origin is at center of the screen.
    % Centering the dots will simplify math (when scaling, shearing, etc... in All Dots Mode)
    centeredDots(:,1) = dots(:,1) - xmid*ones(nDots,1);
    centeredDots(:,2) = dots(:,2) - ymid*ones(nDots,1);

    %% All dots mode:
    %% In all dots mode, all dots can be manipulated using keyboard commands. These are usually characters.
    %% Lowercase characters reorders the dots in a horizontal direction, whereas uppercase will perform the
    %% same transformation in the vertical direction
    while (adjustMode == 1)
        % Translated centeredDots back to screenspace, origin at
        % top-left
        dots(:,1) = centeredDots(:,1) + xmid*ones(nDots,1);
        dots(:,2) = centeredDots(:,2) + ymid*ones(nDots,1);
        
        % Loop across stereo buffers
        for view = 0:sign(scal.stereoMode)
            Screen('SelectStereoDrawbuffer', scal.windowPtr, view);
            
            % Fill to black
            Screen('FillRect', scal.windowPtr, colorBlack);
            
            % Reference texture assigned to draw as backdrop?
            if isfield(scal, 'refTex')
                % Seems so. Double-check and draw it, if ok:
                if scal.refTex
                    Screen('DrawTexture', scal.windowPtr, scal.refTex, [], Screen('Rect', scal.refTex));
                end
            end
            
            % Update dot diameters
            % Set all dots to default diameter
            dotDiamArray = [];
            dotDiamArray = zeros(size(dots,1), 1) + dotDiam;
            xGlobalCalibDots = size(scal.XCALIBDOTS, 1);
            xCurrDots = size(dots, 1);
            if xGlobalCalibDots ~= dots
                % Adjusting fitting dots, a center dot may not be selected.
                %  Draw a crosshair at the center instead.
                %bvlDrawCrosshair(windowPtr, x, y [,size] [,color])
                bvlDrawCrosshair(scal.windowPtr, xmax, ymax, selectedDotDiam * 2, dotColor);
            else
                % Make the center dot larger so user can line up with screen
                % center.
                indexCenterDot = round(length(dots(:,1)) / 2);
                dotDiamArray(indexCenterDot) = selectedDotDiam;
            end
            
            
            % Draw dots
            Screen('DrawDots', scal.windowPtr, transpose(dots), dotDiamArray, dotColor, [0 0], dotStyle);
            % HideCursor;
            
            % Show help text if requested:
            if fShowHelptext
                helptxt = ['GLOBAL MODE:\nFinish calibration by pressing ESCape key\nToggle help text by pressing SPACE key\nChange selected global parameters via Cursor arrow keys\n' ...
                    'Select global parameter change mode via right mouse button\nChange parameter step size via middle mouse button\n' ...
                    'Switch to local adjustment mode via left mouse button.\n'];
                DrawFormattedText(scal.windowPtr, helptxt, 30, 30, [255 255 0]);
            end
        end

        % Check for input
        if bMiddleMouseDown || bRightMouseDown
            [strInputName, xMouse, yMouse] = bvlWaitForInput(scal,0.010, 0.010);
        else
            [strInputName, xMouse, yMouse] = bvlWaitForInput(scal,0.010);
        end

        if strcmp(strInputName, 'ESCAPE')| scal.isDONTSTOP,
            % Stop calibrating dots and return to Matlab
            adjustMode = 0;
            Beeper('med', fSoftBeepVolume);

        elseif strcmp(strInputName, 'MiddleMouse') | strcmp(strInputName, 'm')
            % Step through the Adjustment Delta value
            fCurrTime = GetSecs;
            % We put in some time delays to give user time to read output
            % and select the correct input as we cycle through the choices
            if fLastMiddleMouseTime == 0
                fLastMiddleMouseTime = fCurrTime;
            elseif (fCurrTime - fLastMiddleMouseTime) > fTimeDelayForFeedback
                % if we hold down the mouse button, we'll step through the
                % choices and can easily select the option we want
                iAdjustDeltaIndex = iAdjustDeltaIndex + 1;
                if (iAdjustDeltaIndex > length(fAdjustDelta))
                    iAdjustDeltaIndex = 1;
                end
                fDelta = fAdjustDelta(iAdjustDeltaIndex);

                % update timestamp of last time we were here
                fLastMiddleMouseTime = fCurrTime;
            end

            strMsg = 0;
            strMsg = sprintf('Delta: ( %f )', fAdjustDelta(iAdjustDeltaIndex));
            bvlDrawText(scal.windowPtr, strMsg, iTextPos(1), iTextPos(2), gMirrorView, iTextScale);

            bMiddleMouseDown = 1;

        elseif strcmp(strInputName, 'RightMouse') | strcmp(strInputName, 'r')
            % Step through Adjustment Types
            fCurrTime = GetSecs;
            % We put in some time delays to give user time to read output
            % and select the correct input as we cycle through the choices
            if fLastRightMouseTime == 0
                fLastRightMouseTime = fCurrTime;
            elseif (fCurrTime - fLastRightMouseTime) > fTimeDelayForFeedback
                % if we hold down the mouse button, we'll step through the
                % choices and can easily select the option we want
                iAdjustTypeIndex = iAdjustTypeIndex + 1;
                if (iAdjustTypeIndex > length(strAdjustTypes))
                    iAdjustTypeIndex = 1;
                end

                % update timestamp of last time we were here
                fLastRightMouseTime = fCurrTime;
            end
            strMsg = 0;
            strMsg = sprintf('Adjust: ( %s )', strAdjustTypes{iAdjustTypeIndex});
            bvlDrawText(scal.windowPtr, strMsg, iTextPos(1), iTextPos(2), gMirrorView, iTextScale);

            bRightMouseDown = 1;

        elseif strcmp(strInputName, 'LeftMouse') | strcmp(strInputName, 'l')
            %fprintf('DEBUGGING:  Left Mouse Input\n');
            % Switch to Local Edit mode
            ShowCursor;
            adjustMode = 2;
            Beeper('med', fSoftBeepVolume);

        elseif strcmp(strInputName, 'LeftArrow')
            % Left arrow button
            switch iAdjustTypeIndex
                case 1  % Translation
                    if gMirrorView
                        centeredDots = bvlTranslateDots(centeredDots, 'right', fDelta);
                    else
                        centeredDots = bvlTranslateDots(centeredDots, 'left', fDelta);
                    end
                case 2  % Expansion or Contraction
                    centeredDots = bvlExpandContractDots(scal,centeredDots, 'horiz_out', fDelta);
                case 3  % Shear
                    if gMirrorView
                        centeredDots = bvlShearDots(scal,centeredDots, 'top2right', fDelta);
                    else
                        centeredDots = bvlShearDots(scal,centeredDots, 'top2left', fDelta);
                    end
                case 4  % Trapezoid
                    if gMirrorView
                        centeredDots = bvlTrapezoidDots(scal,centeredDots, 'right_back', fDelta);
                    else
                        centeredDots = bvlTrapezoidDots(scal,centeredDots, 'right_forward', fDelta);
                    end
                case 5  % Pincushion or Barrel
                    centeredDots = bvlBarrelPincushionDots(scal,centeredDots, 'horiz_out', fDelta);
            end
            Beeper('med', fSoftBeepVolume);

        elseif strcmp(strInputName, 'RightArrow')
            % Right arrow button
            switch iAdjustTypeIndex
                case 1  % Translation
                    if gMirrorView
                        centeredDots = bvlTranslateDots(centeredDots, 'left', fDelta);
                    else
                        centeredDots = bvlTranslateDots(centeredDots, 'right', fDelta);
                    end
                case 2  % Expansion or Contraction
                    centeredDots = bvlExpandContractDots(scal,centeredDots, 'horiz_in', fDelta);
                case 3  % Shear
                    if gMirrorView
                        centeredDots = bvlShearDots(scal,centeredDots, 'top2left', fDelta);
                    else
                        centeredDots = bvlShearDots(scal,centeredDots, 'top2right', fDelta);
                    end
                case 4  % Trapezoid
                    if gMirrorView
                        centeredDots = bvlTrapezoidDots(scal,centeredDots, 'right_forward', fDelta);
                    else
                        centeredDots = bvlTrapezoidDots(scal,centeredDots, 'right_back', fDelta);
                    end
                case 5  % Pincushion or Barrel
                    centeredDots = bvlBarrelPincushionDots(scal,centeredDots, 'horiz_in', fDelta);
            end
            Beeper('med', fSoftBeepVolume);

        elseif strcmp(strInputName, 'UpArrow')
            % Up arrow button
            switch iAdjustTypeIndex
                case 1  % Translation
                    centeredDots = bvlTranslateDots(centeredDots, 'up', fDelta);
                case 2  % Expansion or Contraction
                    centeredDots = bvlExpandContractDots(scal,centeredDots, 'vert_out', fDelta);
                case 3  % Shear
                    if gMirrorView
                        centeredDots = bvlShearDots(scal,centeredDots, 'right2down', fDelta);
                    else
                        centeredDots = bvlShearDots(scal,centeredDots, 'right2up', fDelta);
                    end
                case 4  % Trapezoid
                    centeredDots = bvlTrapezoidDots(scal,centeredDots, 'top_back', fDelta);
                case 5  % Pincushion or Barrel
                    centeredDots = bvlBarrelPincushionDots(scal,centeredDots, 'vert_out', fDelta);
            end
            Beeper('med', fSoftBeepVolume);

        elseif strcmp(strInputName, 'DownArrow')
            % Down arrow button
            switch iAdjustTypeIndex
                case 1  % Translation
                    centeredDots = bvlTranslateDots(centeredDots, 'down', fDelta);
                case 2  % Expansion or Contraction
                    centeredDots = bvlExpandContractDots(scal,centeredDots, 'vert_in', fDelta);
                case 3  % Shear
                    if gMirrorView
                        centeredDots = bvlShearDots(scal,centeredDots, 'right2up', fDelta);
                    else
                        centeredDots = bvlShearDots(scal,centeredDots, 'right2down', fDelta);
                    end
                case 4  % Trapezoid
                    centeredDots = bvlTrapezoidDots(scal,centeredDots, 'top_forward', fDelta);
                case 5  % Pincushion or Barrel
                    centeredDots = bvlBarrelPincushionDots(scal,centeredDots, 'vert_in', fDelta);
            end
            Beeper('med', fSoftBeepVolume);
        elseif strcmp(strInputName, 'space')
            % Space bar: Toggle help text display:
            fShowHelptext = 1 - fShowHelptext;
        else
            % Clear mouse flags
            if bMiddleMouseDown
                bMiddleMouseDown = 0;
            end
            if bRightMouseDown
                bRightMouseDown = 0;
            end
        end

% bv220108 commented out because we do not need it for the haploscopes, however, someone might need it at a certain point....
%         % Call Auxiliary function if we have one
%         if ~isempty(gAuxiliaryFunc)
%             msgAuxText = eval(gAuxiliaryFunc);
%             DrawAuxText(windowPtr, msgAuxText, gMirrorView, iTextScale);
%         end


        Screen('Flip', scal.windowPtr);

    end   %while (adjustMode == 1)

    %% Dot Selection Mode:
    %% The user can select a single dot, or a group of dots using the
    %% left mouse button.  Then adjust them using the arrow keys
    %% Selected dots are "highlighted"
    while (adjustMode == 2)

        % Update dot diameters
        % Set all dots to default diameter
        dotDiamArray = [];
        dotDiamArray = zeros(size(dots,1), 1) + dotDiam;
        % Update diameters for selected dots
        dotDiamArray(indicesSelectedDots) = selectedDotDiam;

        % Fill to black
        Screen('FillRect', scal.windowPtr, colorBlack);

        % Reference texture assigned to draw as backdrop?
        if isfield(scal, 'refTex')
            % Seems so. Double-check and draw it, if ok:
            if scal.refTex
                Screen('DrawTexture', scal.windowPtr, scal.refTex, [], Screen('Rect', scal.refTex));
            end
        end
        
        % Draw dots
        %Screen('DrawDots', windowPtr, xy [,size] [,color] [,center] [,dot_type]);
        Screen('DrawDots', scal.windowPtr, transpose(dots), dotDiamArray, dotColor, [0 0], dotStyle);

        % Show help text if requested:
        if fShowHelptext
            helptxt = ['LOCAL MODE:\nFinish calibration by pressing ESCape key\nToggle help text by pressing SPACE key\nMove selected calibration dot(s) via Cursor arrow keys\n' ...
                'Select dot or area of dots via left mouse button + mouse drag\nUnselect dot(s) via right mouse button\nChange dot movement size via middle mouse button\n' ...
                'Switch to global mode via single press of right mouse button\nwhile no dots selected.\n'];
            DrawFormattedText(scal.windowPtr, helptxt, 30, 30, [255 255 0]);
        end

        
        % Check for input
        if bLeftMouseDown || bMiddleMouseDown
            [strInputName, xMouse, yMouse] = bvlWaitForInput(scal,0.010, 0.010);
        else
            [strInputName, xMouse, yMouse] = bvlWaitForInput(scal,0.010);
        end


        % Handle specific keys we care about
        if strcmp(strInputName, 'ESCAPE')
            % Stop calibrating dots and return to Matlab
            adjustMode = 0;
            indicesSelectedDots = [];
            Beeper('med', fSoftBeepVolume);

        elseif strcmp(strInputName, 'RightMouse') | strcmp(strInputName, 'r')
            if ~isempty(indicesSelectedDots)
                % exit select mode only so use can select more dots
                ShowCursor;
            else
                % Exit Single Dot Mode - return to All Dot Mode
                adjustMode = 1;
            end
            % reset our selection array
            indicesSelectedDots = [];
            Beeper('med', fSoftBeepVolume);

        elseif strcmp(strInputName, 'MiddleMouse') | strcmp(strInputName, 'm')
            % Step through the Adjustment Delta value
            fCurrTime = GetSecs;
            % We put in some time delays to give user time to read output
            % and select the correct input as we cycle through the choices
            if fLastMiddleMouseTime == 0
                fLastMiddleMouseTime = fCurrTime;
            elseif (fCurrTime - fLastMiddleMouseTime) > fTimeDelayForFeedback
                % if we hold down the mouse button, we'll step through the
                % choices and can easily select the option we want
                iAdjustDeltaIndex = iAdjustDeltaIndex + 1;
                if (iAdjustDeltaIndex > length(fAdjustDelta))
                    iAdjustDeltaIndex = 1;
                end
                fDelta = fAdjustDelta(iAdjustDeltaIndex);

                % update timestamp of last time we were here
                fLastMiddleMouseTime = fCurrTime;
            end

            strMsg = 0;
            strMsg = sprintf('Adjustment Delta: ( %f )', fAdjustDelta(iAdjustDeltaIndex));
            bvlDrawText(scal.windowPtr, strMsg, iTextPos(1), iTextPos(2), gMirrorView, iTextScale);

            bMiddleMouseDown = 1;

        elseif strcmp(strInputName, 'LeftMouse') | strcmp(strInputName, 'l')
            % Don't accept LeftMouse input if we're currently adjusting
            % dots
            if isempty(indicesSelectedDots)
                % exit select mode only so use can select more dots
                if ~bLeftMouseDown
                    % Mouse is just pressed down, record starting position
                    xMouseStart = xMouse;
                    yMouseStart = yMouse;
                else
                    % Draw rectangle frame
                    DrawRectFrame(scal.windowPtr, [xMouseStart yMouseStart xMouse yMouse], colorCyan);
                end

                % Set flag specifying the left-mouse is down
                bLeftMouseDown = 1;
            end

        elseif strcmp(strInputName, 'space')
            % Space bar: Toggle help text display:
            fShowHelptext = 1 - fShowHelptext;            
        elseif strcmp(strInputName, 'LeftArrow')
            if ~isempty(indicesSelectedDots)
                % Dots selected, adjust
                if gMirrorView
                    dots(indicesSelectedDots, :) = bvlTranslateDots(dots(indicesSelectedDots, :), 'right', fDelta);
                else
                    dots(indicesSelectedDots, :) = bvlTranslateDots(dots(indicesSelectedDots, :), 'left', fDelta);
                end
                Beeper('med', fSoftBeepVolume);
            end
        elseif strcmp(strInputName, 'RightArrow')
            if ~isempty(indicesSelectedDots)
                % Dots selected, adjust
                if gMirrorView
                    dots(indicesSelectedDots, :) = bvlTranslateDots(dots(indicesSelectedDots, :), 'left', fDelta);
                else
                    dots(indicesSelectedDots, :) = bvlTranslateDots(dots(indicesSelectedDots, :), 'right', fDelta);
                end
                Beeper('med', fSoftBeepVolume);
            end
        elseif strcmp(strInputName, 'UpArrow')
            if ~isempty(indicesSelectedDots)
                % Dots selected, adjust
                dots(indicesSelectedDots, :) = bvlTranslateDots(dots(indicesSelectedDots, :), 'up', fDelta);
                Beeper('med', fSoftBeepVolume);
            end
        elseif strcmp(strInputName, 'DownArrow')
            if ~isempty(indicesSelectedDots)
                % Dots selected, adjust
                dots(indicesSelectedDots, :) = bvlTranslateDots(dots(indicesSelectedDots, :), 'down', fDelta);
                Beeper('med', fSoftBeepVolume);
            end
        else
            if bLeftMouseDown
                bLeftMouseDown = 0;
                % Just finished a mouse selection
                % Need to test if user selected one dot or many
                % Then highlight the dot(s) and wait for input to move them

                % reset our selection array
                indicesSelectedDots = [];

                dx = abs(xMouseStart - xMouse);
                dy = abs(yMouseStart - yMouse);
                if ((dx < dotDiam) & (dy < dotDiam))
                    % Small mouse movement, user probably attempting to
                    % select an individual dot
                    location = [xMouse yMouse];
                    dotDiffs = dots - ones(nDots, 1) * location;
                    dotDists = sqrt(dotDiffs(:,1).^2 + dotDiffs(:,2).^2);
                    [dotDistance, indicesSelectedDots] = min(dotDists);
                    if (dotDistance > dotDiam)
                        % Not close enough to a dot
                        indicesSelectedDots = [];
                    end
                else
                    % Large mouse movement, user probably selecting several
                    % dots
                    xMouseMin = min([xMouseStart xMouse]);
                    xMouseMax = max([xMouseStart xMouse]);
                    yMouseMin = min([yMouseStart yMouse]);
                    yMouseMax = max([yMouseStart yMouse]);

                    % Find all dots with x-coord in the selected rectangle
                    xChosenDots = dots(:,1) > xMouseMin;
                    xTemp       = dots(:,1) < xMouseMax;
                    xChosenDots = xChosenDots & xTemp;
                    % Find all dots with y-coord in the selected rectangle
                    yChosenDots = dots(:,2) > yMouseMin;
                    yTemp       = dots(:,2) < yMouseMax;
                    yChosenDots = yChosenDots & yTemp;
                    % Find all dots that are within x and y
                    indicesSelectedDots = xChosenDots & yChosenDots;
                    % Grab the indices
                    indicesSelectedDots = find(indicesSelectedDots);
                end

                if isempty(indicesSelectedDots)
                    % Signal audio tone to let the user know they need
                    % to try again
                    Beeper('low');
                else
                    % HideCursor so the can focus on the dot's
                    % HideCursor;
                end

                % Reset mouse positions
                xMouseStart = 0;
                yMouseStart = 0;
            end % if bLeftMouseDown

            if bMiddleMouseDown
                bMiddleMouseDown = 0;
            end
        end % if strcmp(strInputName, 'esc')
        
%         % Call Auxiliary function if we have one
%         if ~isempty(gAuxiliaryFunc)
%             msgAuxText = eval(gAuxiliaryFunc);
%             DrawAuxText(scal.windowPtr, msgAuxText, gMirrorView, iTextScale);
%         end

        % flip frame buffer to show all draw operations that happened in
        % this loop.
        Screen('Flip', scal.windowPtr);

    end % while (adjustMode == 2)
end % while (adjustMode ~= 0)

% Cleanup
% MK: Disabled - Done by calling function...
% ShowCursor;
% Screen('Close', scal.windowPtr);

% Set timestamp to reenable mouse input
bvlAcceptCalibGuiInput(GetSecs);

if isempty(scal.SELECTXCALIBDOTS),  % all dots
    scal.XCALIBDOTS(:) = dots(:,1);
    scal.YCALIBDOTS(:) = dots(:,2);
else
    scal.SELECTXCALIBDOTS(:) = dots(:,1);
    scal.SELECTYCALIBDOTS(:) = dots(:,2); % or only dot selection
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw optional text
function DrawAuxText(windowPtr, msgAuxText, gMirrorView, iTextScale)
windowRect      = Screen('Rect', windowPtr);
windowWidth     = windowRect(3);
windowHeight    = windowRect(4);
auxTextPos  = [(windowWidth * 0.3) (windowHeight * 0.3)];
bvlDrawText(windowPtr, msgAuxText, auxTextPos(1), auxTextPos(2), gMirrorView, iTextScale);


