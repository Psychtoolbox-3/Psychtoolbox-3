% RenderDemo
%
% Illustrates calibration interface for simple task of producing a uniform
% color patch of desired CIE xyY coordinates.
%
% The calculation is done with respect to the current PTB demonstration
% calibration file.
%
% The demo shows multiple different ways to implement this, starting with a
% purely Matlab based method, progressing to more advanced methods. The
% final demonstration shows how to do it automatically and graphics
% hardware accelerated.
%
% Demo 1:
%
% The RGB values are gamma corrected and live in the range [0,1].  If they
% contain 0 or 1, the xyY coordinates requested may have been out of gamut.
%
% A uniform color patch is displayed in the MATLAB figure window. This is
% not a well-controlled display method, but does give a sense of the patch
% color if the calibration file is a reasonable description of the display.
%
% Immediately afterwards, the same color patch is shown in a PTB onscreen
% window, with the same gamma table loaded which was used during
% calibration measurements. This should render an accurate stimulus.
%
% Demo 2: As demo 1, but displaying in a onscreen window and performing the
% gamma correction via proper inverse gamma lookup tables loaded into the
% graphics card, thereby presenting on a linearized display, instead of
% using the SensorToSettings() routine to adapt the stimulus to a
% non-linearized display.
%
% The last two demos Demo 3 and Demo 4 require a recent graphics card and
% perform all color space conversions and calibrated display automatically
% and hardware accelerated on the graphics card. Any NVidia GeForce-8000 or
% later, AMD Radeon X-1000 or later, or Intel HD graphics card should be
% able to support these demos.
%
% Demo 3: The stimulus is defined in XYZ tristimulus color space and
% converted automatically by Screen() into RGB output format, taking the
% calibration data in 'cal' into account.
%
% Demo 4: The stimulus is directly defined in xyY chromacity + luminance
% format and all conversions and calibrations are done automatically by
% Screen().
%
% 4/26/97  dhb  Wrote it.
% 7/25/97  dhb  Better initialization.
% 3/12/98  dgp  Use Ask.
% 3/14/02  dhb  Update for OpenWindow.
% 4/03/02  awi  Merged in Windows changes.  On Windows we do not copy the result to the clipboard. 
% 4/13/02  awi	Changed "SetColorSpace" to new name "SetSensorColorSpace".
%				Changed "LinearToSettings" to new name "SensorToSettings".
% 12/21/02 dhb  Remove reliance on now obsolete OpenWindow/CloseWindow.
% 11/16/06 dhb  Start getting this to work with PTB-3.
% 11/22/06 dhb  Fixed except that Ask() needs to be fixed.
% 6/16/11  dhb  The PTB display section was out of date and didn't work.  I removed it.
% 1/26/13  mk   Add standard PTB display, but also imaging pipeline based methods.

% Clear out workspace
clear

% Load default calibration file:
cal = LoadCalFile('PTB3TestCal');
load T_xyz1931
T_xyz1931 = 683*T_xyz1931;
cal = SetSensorColorSpace(cal,T_xyz1931,S_xyz1931);
cal = SetGammaMethod(cal,0);

% Get xyY, render, and report.
xyY = input('Enter xyY (as a row vector) default [.3 .3 50]: ')';
if isempty(xyY)
    xyY = [.3 .3 50]';
end

% Pure software conversion:
XYZ = xyYToXYZ(xyY);
[RGB, outOfRangePixels] = SensorToSettings(cal, XYZ);
fprintf('Computed RGB: [%g %g %g]\n', RGB(1), RGB(2), RGB(3));

% Check for out-of-range non-displayable color values:
if any(outOfRangePixels)
    fprintf('WARNING: Out of range RGB values -- not displayable!\n');
end

% Make it an image
nX = 256; nY = 128;
theRGBCalFormat = RGB*ones(1,nX*nY);
theRGBImage = CalFormatToImage(theRGBCalFormat,nX,nY);

% Show in a Matlab figure window.  This will not be calibrated,
% but gives the general sense.  Use PTB display routines for
% more precise display.
figure; clf;
h = image(theRGBImage);
title('Here is the color');
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'XTick',[]);
set(gca,'YTick',[])
drawnow;

fprintf('\n\nPress any key to continue. This will do the same thing in a set of Psychtoolbox onscreen windows.\n');
KbStrokeWait(-1);
close all;
drawnow;

try
    % Open window in GUI mode, top-left, 300 x 300 pixels:
    AssertOpenGL;
    
    % Define our desired background color in RGB primary space:
    % RGB = [0.01, 0.01, 0.01] for almost black.
    bgcolor = [0.01; 0.01; 0.01];

    % Declutter our output for this demo:
    Screen('Preference', 'SuppressAllWarnings', 1);
    
    % Skip display sync tests:
    oldsync = Screen('Preference', 'SkipSyncTests', 2);

    % Select display screen to show windows:
    screenId = max(Screen('Screens'));
    
    [win0, winRect0] = Screen('OpenWindow', screenId, bgcolor * 255, [0 0 300 300], [], [], [], [], [], kPsychGUIWindow);
    
    % Load the gamma table which was used during calibration measurements:
    % If this is a 1024 slot table we downsample to 256 slots, so it works
    % with MS-Windows, otherwise we hope it is a compatible table.
    % We could do better here, but this is just a demo...
    if length(cal.gammaInput) == 1024
        gammaInput = cal.gammaInput(1:4:end);
    else
        gammaInput = cal.gammaInput;
    end
    
    % Replicate to 3 columns for the three primary colors:
    gammaInput = repmat(gammaInput, 1, 3);
    
    % Before table upload, we store a backup copy of the original table, so
    % it can get restored at end of session:
    BackupCluts(screenId);
    Screen('LoadNormalizedGammaTable', screenId, gammaInput);
    
    % Convert new theRGBImage to texture and draw it into win1:
    tex = Screen('MakeTexture', win0, round(theRGBImage * 255));
    Screen('DrawTexture', win0, tex);
    Screen('Close', tex);
    
    % Show it:
    Screen('Flip', win0);
    dstRect0 = CenterRect([0 0 nX nY], winRect0);
    readBack0 = Screen('GetImage', win0, dstRect0);
    
    fprintf('\n\nPress any key to continue. This will demonstrate another way, using SensorToPrimary() + \n');
    fprintf('a proper inverse gamma table, to linearize your display, instead of SensorToSettings().\n\n');
    KbStrokeWait(-1);
    
    % Close old window, as its content is not compatible with the gamma
    % table we're gonna set now:
    Screen('CloseAll');
    
    % Show same thing in a GUI window of 300 x 300 pixels.
    fprintf('Now we do exactly the same thing, just displaying in a onscreen window.\n');
    fprintf('However, we use gamma correction via the graphics hardware, so we have a linearized\n');
    fprintf('display. This allows to use the simpler SensorToPrimary() instead of SensorToSettings().\n\n');
    
    % Open a standard window:
    [win1, winRect1] = Screen('OpenWindow', screenId, bgcolor * 255, [0 0 300 300], [], [], [], [], [], kPsychGUIWindow);
    
    % Load a gamma correction table into the graphics card, as defined as
    % the inverse gamma table for given measured display gamma table
    % 'cal.gammaTable'. However, we sub-sample the table to 256 slots to
    % make sure it works on MS-Windows, not only on OSX or Linux:
    iGammaTable = InvertGammaTable(cal.gammaInput, cal.gammaTable, 256);
    
    % Load inverse gamma table into GPU:
    Screen('LoadNormalizedGammaTable', screenId, iGammaTable);
    
    % Ok, now we have a linearized display due to gamma correction. This means
    % we can define our stimulus in tristimulus XYZ space. This allows us to
    % use the simpler SensorToPrimary() function instead of the more complex
    % SensorToSettings() function:
    XYZ = xyYToXYZ(xyY);
    RGB = SensorToPrimary(cal, XYZ);
    fprintf('Recomputed linear RGB: [%g %g %g]\n', RGB(1), RGB(2), RGB(3));

    % Check for out-of-range non-displayable color values:
    if any(RGB < 0 | RGB > 1)
        fprintf('WARNING: Out of range RGB values -- not displayable!\n');
    end
    
    % Make it an image. Now need to scale by 255, as onscreen windows want
    % color values in range 0 - 255 instead of 0 - 1 by default:
    nX = 256; nY = 128;
    theRGBCalFormat = RGB * 255 * ones(1,nX*nY);
    theRGBImage = CalFormatToImage(theRGBCalFormat,nX,nY);
    
    % Convert new theRGBImage to texture and draw it into win1:
    tex = Screen('MakeTexture', win1, round(theRGBImage));
    Screen('DrawTexture', win1, tex);
    Screen('Close', tex);
    
    % Show it:
    Screen('Flip', win1);
    dstRect1 = CenterRect([0 0 nX nY], winRect1);
    readBack1 = Screen('GetImage', win1, dstRect1);
    
    fprintf('\n\nPress any key to continue. This will demonstrate a simpler way to do it via the imaging pipeline.\n');
    fprintf('Screen() will automatically convert XYZ tristimulus color values to calibrated RGB values before display.\n\n');
    KbStrokeWait(-1);
    
    % Make sure this will actually work:
    AssertGLSL;
    
    % Open a 2nd window, now using the imaging pipeline:
    PsychImaging('PrepareConfiguration');
    
    % Enable 32 bpc floating point framebuffer, so fractional color values
    % can be represented accurately. We will store XYZ tristimulus color
    % values in the frambuffer, not RGB values:
    PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');
    
    % Also use unrestricted color range for writing arbitrary color values
    % to the framebuffer:
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
    
    % It shall use builtin fast SensorToPrimary() plugin:
    PsychImaging('AddTask', 'AllViews', 'DisplayColorCorrection', 'SensorToPrimary');
    
    % Check for valid (displayable) final color values in 0.0 - 1.0 range.
    % Mark out-of-range pixels visually:
    PsychImaging('AddTask', 'AllViews', 'DisplayColorCorrection', 'CheckOnly');
    
    % Open it: Our window operates in XYZ color space, so we need to define
    % a XYZ background input color that leads to our desired background
    % color 'bgcolor':
    background = PrimaryToSensor(cal, bgcolor);
    [win2, winRect2] = PsychImaging('OpenWindow', screenId, background, [310 0 610 300], [], [], [], [], [], kPsychGUIWindow);
    
    % Assign 'cal' struct for XYZ -> RGB conversion:
    PsychColorCorrection('SetSensorToPrimary', win2, cal);
    
    % Compared to above, we can skip the SensorToPrimary step:
    XYZ = xyYToXYZ(xyY);
    
    % Simply draw to the framebuffer, directly in XYZ format instead of RGB:
    % We can use fillrect to draw the patch, without intermediate need for
    % textures:
    dstRect2 = CenterRect([0 0 nX nY], winRect2);
    Screen('FillRect', win2, XYZ, dstRect2);
    
    % Readback image in XYZ format from framebuffer:    
    readBack2In = Screen('GetImage', win2, [], 'drawBuffer', 1);
    
    % Show it:
    Screen('Flip', win2);
    
    % Read back final image from framebuffer, for correctness check:    
    readBack2 = Screen('GetImage', win2, dstRect2);
    
    % Plot manual and automatic result for comparison:
    close all;
    imshow(readBack1);
    figure;
    imshow(readBack2);
    
    fprintf('\n\nPress any key to continue. This will demonstrate the most simple way to do it via the imaging pipeline.\n');
    fprintf('This method allows to draw and define your stimulus completely in the xyY chromacity+luminance color space.\n');
    fprintf('Screen() will automatically convert your xyY color values to proper RGB framebuffer values before display.\n\n');
    KbStrokeWait(-1);
    
    % Now the same thing, but we draw colors directly in (x,y) chromacity and Y
    % luminance format [x,y,Y] into the framebuffer. The imaging pipeline
    % will do the complete conversion from xyY space to XYZ space and then
    % XYZ space to RGB space, followed by gamma correction by the graphics
    % card for display linearization:
    % Open a 2nd window, now using the imaging pipeline:
    PsychImaging('PrepareConfiguration');
    
    % Enable 32 bpc floating point framebuffer, so fractional color values
    % can be represented accurately. We will store xyY chromacity +
    % luminance color values in the frambuffer, not RGB values:
    PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');
    
    % Also use unrestricted color range for writing arbitrary color values
    % to the framebuffer:
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');

    % It shall use builtin fast xyYToXYZ() plugin for xyY -> XYZ conversion:
    PsychImaging('AddTask', 'AllViews', 'DisplayColorCorrection', 'xyYToXYZ');
    
    % It shall use builtin fast SensorToPrimary() plugin:
    PsychImaging('AddTask', 'AllViews', 'DisplayColorCorrection', 'SensorToPrimary');

    % Check for valid (displayable) final color values in 0.0 - 1.0 range.
    % Mark out-of-range pixels visually:
    PsychImaging('AddTask', 'AllViews', 'DisplayColorCorrection', 'CheckOnly');
    
    % Open it: Our window operates in xyY color space, so we need to define
    % a xyY background input color that leads to our desired background
    % color 'bgcolor':
    background = XYZToxyY(PrimaryToSensor(cal, bgcolor));
    [win3, winRect3] = PsychImaging('OpenWindow', screenId, background, [0 330 300 630], [], [], [], [], [], kPsychGUIWindow);
    
    % Assign 'cal' struct for XYZ -> RGB conversion:
    PsychColorCorrection('SetSensorToPrimary', win3, cal);
    
    % Simply draw to the framebuffer, directly in xyY format:
    dstRect3 = CenterRect([0 0 nX nY], winRect3);
    Screen('FillRect', win3, xyY, dstRect3);
    
    % Readback image in xyY format from framebuffer:
    readBack3In = Screen('GetImage', win3, dstRect3, 'drawBuffer', 1);
    
    % Show it:
    Screen('Flip', win3);
    
    % Read back final image from framebuffer, for correctness check:
    readBack3 = Screen('GetImage', win3, dstRect3);

    % Plot manual and automatic result for comparison:    
    close all;
    imshow(readBack1);
    figure;
    imshow(readBack3);

    fprintf('Press any key to end the demo.\n');
    KbStrokeWait(-1);
    
    % sca closes all onscreen windows and restores the original gamma tables:
    sca;
    
    Screen('Preference', 'SuppressAllWarnings', 0);
    Screen('Preference', 'SkipSyncTests', oldsync);

catch %#ok<CTCH>
    sca;
    Screen('Preference', 'SuppressAllWarnings', 0);
    Screen('Preference', 'SkipSyncTests', oldsync);
    psychrethrow(psychlasterror);
end
