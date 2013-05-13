function AlphaBlendingTest(screenNumber)

% AlphaBlendingTest([screenNumber])
%
% Perform tests of OpenGL alpha blending by drawing to the screen using
% 'PutImage', reading back the actual values on the Screen using
% 'GetImage' and comparing the actual values to predicted values
% calcualated in MATLAB.
%
% AlphaBlendingTest combines tests implemented separately.  You may perform
% any of these tests individually, or call TestAlphaBlending to perform
% them all:
%   
% AlphaBlendSettingTest - 
%   Test the Screen('BlendFunction') recalls previously stored
%   alpha values.
%
% AlphaMultiplicationTest -
%   Test that alpha multiplication by values 0 and 1 (Screen 255) works
%   with perfect precision.  OpenGL guarantees perfect precision for those
%   alpha values only.
%
% AlphaMultiplicationAccuracyTest - 
%   Measure the precision of alpha values between 0 and 1 (Screen 0 and 255) by
%   drawing to the screen, then taking the difference between what was
%   drawn to the screen and results of simulated blending done with
%   double-precision floats.    
%
% AlphaAdditionTest - 
%   Test that addition of source and destination terms has perfect
%   precision.    


if nargin==0
    screenNumber=max(Screen('Screens'));
end

resultStrings={'Passed', 'Failed'};

fprintf('AlphaBlendingTest, performaing four tests of alpha blending:\n');

fprintf('    1. Testing that Screen(''BlendFunction'') retains and recalls settings: '); 
failFlag1= AlphaBlendSettingTest(screenNumber);
fprintf([resultStrings{failFlag1 + 1} '\n']);

fprintf('    2. Testing accuracy of OpenGL alpha addition:                         '); 
failFlag2= AlphaAdditionTest(screenNumber);
fprintf([resultStrings{failFlag2 + 1} '\n']);

fprintf('    3. Testing accuracy of OpenGL alpha multiplication by 0 and 1:        '); 
failFlag3= AlphaAdditionTest(screenNumber);
fprintf([resultStrings{failFlag3 + 1} '\n']);

fprintf('    4. Testing accuracy of OpenGL multiplication between 0 and 1:         '); 
[multMaximumError, multRoundTypeStr, multIndependentFlag]=AlphaMultiplicationAccuracyTest(screenNumber);
didPass= multMaximumError < 1 && multIndependentFlag;
fprintf([resultStrings{~didPass + 1} '\n']);
fprintf(['        Maximum descrepancy between MATLAB double and OpenGL integer multiplication: ' num2str(multMaximumError) '\n']);         
if strcmp(multRoundTypeStr,  'round')
    fprintf('        Alpha blending rounds multiplication to the nearest integer.\n');
elseif strcmp(multRoundTypeStr,  'floor')
    fprintf('        Alpha blending truncates multiplication.\n');
elseif strcmp(multRoundTypeStr,  'ceil')
    fprintf('        Alpha blending rounds multiplication up.\n');
elseif strcmp(multRoundTypeStr,  'unknown')
    fprintf('        Alpha blending multiplication rounds result by an unidentified method.\n');
end
if (multIndependentFlag)
    fprintf('        Multiplication rules are independent of the surface (source or destination) and the blending factor string.\n');
else   
    fprintf('        Multiplication rules depend on either the surface or the blending factor string.\n');
end


