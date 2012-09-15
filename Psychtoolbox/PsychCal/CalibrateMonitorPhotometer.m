function [ gammaTable1, gammaTable2, displayBaseline, displayRange, displayGamma, maxLevel ] = CalibrateMonitorPhotometer(numMeasures, screenid)
% [gammaTable1, gammaTable2, displayBaseline, displayRange. displayGamma, maxLevel ] = CalibrateMonitorPhotometer([numMeasures=9][, screenid=max])
%
% A simple calibration script for analog photometers.
%
% Use CalibrateMonSpd() if you want to do more fancy calibration with
% different types of photometers or special devices like Bits+ or DataPixx,
% assuming you know how to operate CalibrateMonSpd() that is...
%
% numMeasures (default: 9) readings are taken manually, and the readings
% are fit with a gamma function and piecewise cubic splines. numMeasures -
% 1 should be a power of 2, ideally (9, 17, 33, etc.). The corresponding
% linearized gamma tables (1 -> gamma, 2 -> splines) are returned, as well
% as the display baseline, display range in cd/m^2 and display gamma. Plots
% of the two fits are created as well. Requires fit tools.
%
% If the normalized gamma table is not loaded, then the cd/m^2 value of a
% screen value can be figured out by the formula: cdm2 =
% displayRange*(screenval/maxLevel).^(1/displayGamma) + displayBaseline.
%
% Generally, you will want to load the normalized gamma tables and use them
% in Screen('LoadNormalizedGammaTable'). For example:
%
% [gammaTable1, gammaTable2] = CalibrateMonitorPhotometer;
% %Look at the outputted graphs to see which one gives a better fit
% %Then save the corresponding gamma table for later use
% gammaTable = gammaTable1;
% save MyGammaTable gammaTable
% 
% %Then when you're ready to use the gamma table:
% load MyGammaTable
% Screen('LoadNormalizedGammaTable', win, gammaTable*[1 1 1]);
%
%
% History:
% Version 1.0: Patrick Mineault (patrick.mineault@gmail.com)
% 22.10.2010 mk Switch numeric input from use of input() to use of
%               GetNumber(). Restore gamma table after measurement. Make
%               more robust.
% 19.08.2012 mk Some cleanup.
%  4.09.2012 mk Use Screen('ColorRange') to adapt number/max of intensity
%               level to given range of framebuffer.

global vals;
global inputV;

    if (nargin < 1) || isempty(numMeasures)
        numMeasures = 9;
    end

    input(sprintf(['When black screen appears, point photometer, \n' ...
           'get reading in cd/m^2, input reading using numpad and press enter. \n' ...
           'A screen of higher luminance will be shown. Repeat %d times. ' ...
           'Press enter to start'], numMeasures));
       
    psychlasterror('reset');    
    try
        if nargin < 2 || isempty(screenid)
            % Open black window on default screen:
            screenid = max(Screen('Screens'));
        end
        
        % Open black window:
        win = Screen('OpenWindow', screenid, 0);
        maxLevel = Screen('ColorRange', win);

        % Load identity gamma table for calibration:
        LoadIdentityClut(win);

        vals = [];
        inputV = [0:(maxLevel+1)/(numMeasures - 1):(maxLevel+1)]; %#ok<NBRAK>
        inputV(end) = maxLevel;
        for i = inputV
            Screen('FillRect',win,i);
            Screen('Flip',win);

            fprintf('Value? ');
            resp = GetNumber;
            fprintf('\n');
            vals = [vals resp]; %#ok<AGROW>
        end
        
        % Restore normal gamma table and close down:
        RestoreCluts;
        Screen('CloseAll');
    catch %#ok<*CTCH>
        RestoreCluts;
        Screen('CloseAll');
        psychrethrow(psychlasterror);
    end

    displayRange = range(vals);
    displayBaseline = min(vals);
    
    %Normalize values
    vals = (vals - displayBaseline) / displayRange;
    inputV = inputV/maxLevel;
    
    if ~exist('fittype'); %#ok<EXIST>
        fprintf('This function needs fittype() for automatic fitting. This function is missing on your setup.\n');
        fprintf('Therefore i can''t proceed, but the input values for a curve fit are available to you by\n');
        fprintf('defining "global vals;" and "global inputV" on the command prompt, with "vals" being the displayed\n');
        fprintf('values and "inputV" being user input from the measurement. Both are normalized to 0-1 range.\n\n');
        error('Required function fittype() unsupported. You need the curve-fitting toolbox for this to work.\n');
    end
    
    %Gamma function fitting
    g = fittype('x^g');
    fittedmodel = fit(inputV',vals',g);
    displayGamma = fittedmodel.g;
    gammaTable1 = ((([0:maxLevel]'/maxLevel))).^(1/fittedmodel.g); %#ok<NBRAK>
    
    firstFit = fittedmodel([0:maxLevel]/maxLevel); %#ok<NBRAK>
    
    %Spline interp fitting
    fittedmodel = fit(inputV',vals','splineinterp');
    secondFit = fittedmodel([0:maxLevel]/maxLevel); %#ok<NBRAK>
    
    figure;
    plot(inputV, vals, '.', [0:maxLevel]/maxLevel, firstFit, '--', [0:maxLevel]/maxLevel, secondFit, '-.'); %#ok<NBRAK>
    legend('Measures', 'Gamma model', 'Spline interpolation');
    title(sprintf('Gamma model x^{%.2f} vs. Spline interpolation', displayGamma));
    
    %Invert interpolation
    fittedmodel = fit(vals',inputV','splineinterp');
    
    gammaTable2 = fittedmodel([0:maxLevel]/maxLevel); %#ok<NBRAK>
    
return;
