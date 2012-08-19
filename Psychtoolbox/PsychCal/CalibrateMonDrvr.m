function cal = CalibrateMonDrvr(cal, USERPROMPT, whichMeterType, blankOtherScreen)
% cal = CalibrateMonDrvr(cal,USERPROMPT,whichMeterType,blankOtherScreen)
%
% Main script for monitor calibration.  May be called
% once parameters are set up.
%
% Each monitor input channel is calibrated.
% A summary spectrum is computed.
% Gamma curves are computed.

% 10/26/93	dhb		Wrote it based on CalibrateProj.
% 11/3/93	dhb		Added filename entry with default.
% 2/28/94	dhb		Updated SetMon call to SetColor call.
% 3/12/94	dhb		Created version for monitor 0.
% 					User interface is a little wild.
% 4/3/94	dhb		Save the darkAmbient variable.
% 					User interface improvements
% 9/4/94	dhb		Incorporate gamma fitting
%					improvements from CalibrateMonRoom.
%			dhb		Add whichScreen variable. 
%			dhb		Add sync mode variable.
% 10/20/94	dhb		Add bgColor variable.
% 11/18/94  ccc     Change the range of LUT from (0,255) to 
%                   (0, InputLevels-step) with step=nInputLevels/255 
% 11/21/94	dhb, ccc	Further nine-bit modifications.
% 1/23/95	dhb		Pulled parameter setting out into a calling script,
%					made user prompting conditional.
% 4/12/97	dhb		Update for new toolbox.
% 8/21/97	dhb		Don't save data here.
% 			dhb		Get rid of option not to measure.
% 4/7/99    dhb     NINEBIT -> NBITS.
%           dhb     Handle noMeterAvail, RADIUS switches.
%           dhb     Check for empty indexLin.
% 9/22/99   dhb, mdr  Make boxRect depend on boxSize, defined up one level.
% 10/1/99   dhb, mdr  Pull out nMonBases, defined up one level.
% 12/2/99   dhb     Put background on after white box for aiming.
% 8/14/00   dhb     Call to CMETER('Frequency') only for OS9.
% 8/20/00   dhb     Remove bits arg to SetColor and most RADIUS conditionals.
% 9/11/00   dhb     Remove syncMode code, any direct refs to CMETER.
% 9/14/00   dhb     Use OpenWindow to open.
% 3/8/02    dhb, ly  Call CalibrateManualDrvr if desired.
% 7/9/02    dhb     Get rid of OpenWindow, CloseWindow.
% 9/23/02   dhb, jmh  Force background to zero when box is up for aiming.
% 2/26/03   dhb     Tidy comments.
% 2/3/06	dhb		Allow passing of cal.describe.boxRect
% 10/23/06  cgb     OS/X, etc.
% 11/08/06  dhb, cgb Living in the 0-1 world ....
% 11/10/06  dhb     Get rid of round() around production of input levels.
% 9/26/08   cgb, dhb Fix dacsize when Bits++ is used.  Fit gamma with full number of levels. 
% 8/19/12   dhb     Add codelet suggested by David Jones to clean up at end.  See comment in CalibrateMonSpd.
% 8/19/12   mk      Rewrite setup and clut code to be able to better cope with all
%                   the broken operating systems / drivers / gpus and to also
%                   support DataPixx/ViewPixx devices.

global g_usebitspp;

% If the global flag for using Bits++ is empty, then it hasn't been
% initialized and default it to 0.
if isempty(g_usebitspp)
    g_usebitspp = 0;
end

% Measurement parameters
monWls = SToWls(cal.describe.S); %#ok<*NASGU>

% Define input settings for the measurements
mGammaInputRaw = linspace(0, 1, cal.describe.nMeas+1)';
mGammaInputRaw = mGammaInputRaw(2:end);

% Make manual measurements here if desired.  This needs to come first.
if cal.manual.use
    CalibrateManualDrvr;
end

% User prompt
if USERPROMPT
	if cal.describe.whichScreen == 0
		fprintf('Hit any key to proceed past this message and display a box.\n');
		fprintf('Focus radiometer on the displayed box.\n');
		fprintf('Once meter is set up, hit any key - you will get %g seconds\n',...
                cal.describe.leaveRoomTime);
		fprintf('to leave room.\n');
        KbStrokeWait(-1);
	else
		fprintf('Focus radiometer on the displayed box.\n');
		fprintf('Once meter is set up, hit any key - you will get %g seconds\n',...
                cal.describe.leaveRoomTime);
		fprintf('to leave room.\n');
	end
end

% Blank other screen, if requested:
if blankOtherScreen
    % We simply open an onscreen window with black background color:
    Screen('OpenWindow', cal.describe.whichBlankScreen, 0);
end

% Setup screen to be measured
% ---------------------------

% Prepare imaging pipeline for Bits+ Bits++ CLUT mode, or DataPixx/ViewPixx
% L48 CLUT mode (which is pretty much the same). If such a special output
% device is used, the Screen('LoadNormalizedGammatable', win, clut, 2);
% command uploads 'clut's into the device at next Screen('Flip'), taking
% care of possible graphics driver bugs and other quirks:
PsychImaging('PrepareConfiguration');

if g_usebitspp == 1
    % Setup for Bits++ CLUT mode. This will automatically load proper
    % identity gamma tables into the graphics hardware and into the Bits+:
    PsychImaging('AddTask', 'General', 'EnableBits++Bits++Output');
end

if g_usebitspp == 2
    % Setup for DataPixx/ViewPixx etc. L48 CLUT mode. This will
    % automatically load proper identity gamma tables into the graphics
    % hardware and into the device:
    PsychImaging('AddTask', 'General', 'EnableDataPixxL48Output');
end

% Open the window:
[window, screenRect] = PsychImaging('OpenWindow', cal.describe.whichScreen);
if (cal.describe.whichScreen == 0)
    HideCursor;
end

theClut = zeros(256,3);
if g_usebitspp
    % Load zero theClut into device:
    Screen('LoadNormalizedGammaTable', window, theClut, 2);
    Screen('Flip', window);    
else
    % Load zero lut into regular graphics card:
    Screen('LoadNormalizedGammaTable', window, theClut);
end

% Draw a box in the center of the screen
if ~isfield(cal.describe, 'boxRect')
	boxRect = [0 0 cal.describe.boxSize cal.describe.boxSize];
	boxRect = CenterRect(boxRect,screenRect);
else
	boxRect = cal.describe.boxRect;
end
theClut(2,:) = [1 1 1];
Screen('FillRect', window, 1, boxRect);
if g_usebitspp
    Screen('LoadNormalizedGammaTable', window, theClut, 2);
    Screen('Flip', window, 0, 1);
else
    Screen('LoadNormalizedGammaTable', window, theClut);
end

% Wait for user
if USERPROMPT == 1
    fprintf('Set up radiometer and hit any key when ready\n');
    KbStrokeWait(-1);
    fprintf('Pausing for %d seconds ...', cal.describe.leaveRoomTime);
    WaitSecs(cal.describe.leaveRoomTime);
    fprintf(' done\n');
end

% Put correct surround for measurements.
theClut(1,:) = cal.bgColor';
if g_usebitspp
    Screen('FillRect', window, 1, boxRect);
    Screen('LoadNormalizedGammaTable', window, theClut, 2);
    Screen('Flip', window, 0, 1);
else
    Screen('LoadNormalizedGammaTable', window, theClut);
end

% Start timing
t0 = clock;

mon = zeros(cal.describe.S(3)*cal.describe.nMeas,cal.nDevices);
for a = 1:cal.describe.nAverage
    for i = 1:cal.nDevices
        disp(sprintf('Monitor device %g',i)); %#ok<*DSPS>
        Screen('FillRect', window, 1, boxRect);
        Screen('Flip', window, 0, 1);

        % Measure ambient
        darkAmbient1 = MeasMonSpd(window, [0 0 0]', cal.describe.S, 0, whichMeterType, theClut);

        % Measure full gamma in random order
        mGammaInput = zeros(cal.nDevices, cal.describe.nMeas);
        mGammaInput(i,:) = mGammaInputRaw';
        sortVals = rand(cal.describe.nMeas,1);
        [null, sortIndex] = sort(sortVals); %#ok<*ASGLU>
        %fprintf(1,'MeasMonSpd run %g, device %g\n',a,i);
        [tempMon, cal.describe.S] = MeasMonSpd(window, mGammaInput(:,sortIndex), ...
            cal.describe.S, [], whichMeterType, theClut);
        tempMon(:, sortIndex) = tempMon;

        % Take another ambient reading and average
        darkAmbient2 = MeasMonSpd(window, [0 0 0]', cal.describe.S, 0, whichMeterType, theClut);
        darkAmbient = ((darkAmbient1+darkAmbient2)/2)*ones(1, cal.describe.nMeas);

        % Subtract ambient
        tempMon = tempMon - darkAmbient;

        % Store data
        mon(:, i) = mon(:, i) + reshape(tempMon,cal.describe.S(3)*cal.describe.nMeas,1);
    end
end
mon = mon / cal.describe.nAverage;

% Close the screen, restore cluts:
if g_usebitspp
    % Load identity clut on Bits++ / DataPixx et al.:
    BitsPlusPlus('LoadIdentityClut', window);
    Screen('Flip', window);
end

% Restore graphics card gamma tables to original state:
RestoreCluts;

% Show hidden cursor:
if cal.describe.whichScreen == 0
	ShowCursor;
end

% Close all windows:
Screen('CloseAll');

% Report time
t1 = clock;
fprintf('CalibrateMonDrvr measurements took %g minutes\n', etime(t1, t0)/60);

% Pre-process data to get rid of negative values.
mon = EnforcePos(mon);
cal.rawdata.mon = mon;

% Use data to compute best spectra according to desired
% linear model.  We use SVD to find the best linear model,
% then scale to best approximate maximum
disp('Computing linear models');
cal = CalibrateFitLinMod(cal);

% Fit gamma functions.
cal.rawdata.rawGammaInput = mGammaInputRaw;
cal = CalibrateFitGamma(cal, 2^cal.describe.dacsize);

% Done:
return;
