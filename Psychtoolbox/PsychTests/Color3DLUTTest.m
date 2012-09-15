function Color3DLUTTest
% A simple test/demo of use of 3D CLUT color correction.
%
% Shows some simple test color transformations, performs some CLUT
% conversion accuracy test at the end, measuring how smaller CLUT sizes /
% resolutions affect the precision of color conversion.
%
% Press any key to progress through different tests/demos.
%

% History:
% 21.09.2012  mk  Written.
%
global diffImage;

% Open window with 3D CLUT color correction enabled:
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'AllViews', 'DisplayColorCorrection', 'LookupTable3D');
w=PsychImaging('OpenWindow', 0, 0, [0 0 400 400], [], [], [], [], [], kPsychGUIWindow);

% Build initial identity 3D-CLUT with 16 slots in each color dimension:
clut = ones(3,16,16,16);

% Put red-intensity gradient along red-axis of clut cube:
for i=0:15; clut(1,i+1,:,:) = ones(16,16) * (i/15); end
% Put green-intensity gradient along gree-axis of clut cube:
for i=0:15; clut(2,:,i+1,:) = ones(16,16) * (i/15); end
% Put blue-intensity gradient along blue-axis of clut cube:
for i=0:15; clut(3,:,:,i+1) = ones(16,16) * (i/15); end

% Store original identity lut as baseline:
origclut = clut;

DrawStim(w, clut, 'Identity mapping: gray gradient, increasing left to right.');

% Reduce level of verbosity of PsychColorCorrection() to warnings:
oldverbosity = PsychColorCorrection('Verbosity', 2);

% Invert red intensity gradient:
for i=0:15; clut(1,i+1,:,:) = ones(16,16) * (1 - (i/15)); end
DrawStim(w, clut, 'Identity mapping: Red channel gradient, decreasing left to right.');

% Invert green intensity gradient:
clut = origclut;
for i=0:15; clut(2,:,i+1,:) = ones(16,16) * (1 - (i/15)); end
DrawStim(w, clut, 'Identity mapping: Green channel gradient, decreasing left to right.');

% Invert blue intensity gradient:
clut = origclut;
for i=0:15; clut(3,:,:,i+1) = ones(16,16) * (1 - (i/15)); end
DrawStim(w, clut, 'Identity mapping: Blue channel gradient, decreasing left to right.');

% Only yellow gradient: Clear blue output subvolume:
clut = origclut;
clut(3, :, :, :) = 0;
DrawStim(w, clut, 'Identity mapping: Yellow gradient, blue channel all-zeros.');

% Some basic correctness test:
Screen('FillRect', w, 0);

% Iterate over various CLUT subsampling levels from full resolution
% identity mapping to a unit value cube and measure how precision of color
% mapping degrades with decreasing resolution of the 3D CLUT:
for k = 8:-1:1
    mi = 2^k - 1;
    clut = ones(3,mi+1,mi+1,mi+1);
    
    % Put red-intensity gradient along red-axis of clut cube:
    for i=0:mi; clut(1,i+1,:,:) = ones(mi+1,mi+1) * (i/mi); end
    % Put green-intensity gradient along gree-axis of clut cube:
    for i=0:mi; clut(2,:,i+1,:) = ones(mi+1,mi+1) * (i/mi); end
    % Put blue-intensity gradient along blue-axis of clut cube:
    for i=0:mi; clut(3,:,:,i+1) = ones(mi+1,mi+1) * (i/mi); end
    
    PsychColorCorrection('SetLookupTable3D', w, clut);
    
    % Each RGB input triplet should map to identical output triplet (minus
    % interpolation errors, of course:
    refImage = uint8(rand(400, 400, 3) * 255);
    tex = Screen('MakeTexture', w, refImage);
    Screen('DrawTexture', w, tex);
    Screen('Close', tex);
    
    % Run pipeline, readback result:
    Screen('DrawingFinished', w);
    testImage = Screen('GetImage', w, [], 'backBuffer');
    
    % Show it to user:
    Screen('Flip', w);
    WaitSecs(0.5);
    
    % Compute differences:
    diffImage = testImage - refImage;
    maxdiff = max(max(max(diffImage)));
    mindiff = min(min(min(diffImage)));
    meandiff = mean(mean(mean(diffImage)));

    fprintf('k=%i : CLUT Resolution %i : max = %f min = %f mean = %f error.\n', k, mi+1, maxdiff, mindiff, meandiff);
end

% We're done:
PsychColorCorrection('Verbosity', oldverbosity);
Screen('CloseAll');

return;

end

function DrawStim(w, clut, label)
    % Draw some test stim:
    Screen('FillOval', w, [255 0 0]);

    % A little intensity gradient from rgb (0,0,0) to (255,255,255):
    Screen('DrawLines', w, [[0 ; 200] , [400 ; 200]], 10, [[0 ; 0 ; 0] , [255 ; 255 ; 255]]);

    % Show label for this test:
    Screen('TextSize', w, 14);
    DrawFormattedText(w, label, 'center', 50, 255);
    
    % Upload it:
    PsychColorCorrection('SetLookupTable3D', w, clut);
    
    % And show the color corrected stim:
    Screen('Flip', w);
    
    % Wait until keypress:
    KbStrokeWait;
    
    return;
end
