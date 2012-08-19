% CalibrateManualDrvr
%
% Get some manual measurements from radiometer/photometer.

% 3/8/02  dhb, ly  Wrote it.
% 2/21/03 dhb, ly  Specify input units.

global g_usebitspp;

% If the global flag for using Bits++ is empty, then it hasn't been
% initialized and default it to 0.
if isempty(g_usebitspp)
    g_usebitspp = 0;
end

if USERPROMPT
	if cal.describe.whichScreen == 0
		fprintf('Hit any key to proceed past this message and display a box.\n');
		fprintf('Insert manual photometer/radiometer.\n');
		fprintf('Once meter is set up, hit any key to proceed\n');
		KbStrokeWait(-1);
	else
		fprintf('Insert manual photometer/radiometer.\n');
		fprintf('Once meter is set up, hit any key to proceed\n');
	end
end

% Blank other screen
if blankOtherScreen
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
else
    Screen('LoadNormalizedGammaTable', window, theClut);
end
Screen('Flip', window, 0, 1);

% Wait for user
if USERPROMPT == 1
  KbStrokeWait(-1);
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

% We have put up white:
cal.manual.white = [];
while isempty(cal.manual.white)
	if cal.manual.photometer
		cal.manual.white = input('Enter photometer reading (cd/m2): ');
	else
		cal.manual.white = 1e-6*input('Enter radiometer reading (micro Watts): ');
	end
end

theClut(2,:) = [0 , 0 , 0];
if g_usebitspp
    Screen('LoadNormalizedGammaTable', window, theClut, 2);
    Screen('Flip', window, 0, 1);
else
    Screen('LoadNormalizedGammaTable', window, theClut);
end
cal.manual.black = [];
while isempty(cal.manual.black)
	if cal.manual.photometer
		cal.manual.black = input('Enter photometer reading (cd/m2): ');
	else
		cal.manual.black = 1e-6*input('Enter radiometer reading (micro Watts): ');
	end
end
cal.manual.increment = cal.manual.white - cal.manual.black;

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
