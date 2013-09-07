function PsychDefaultSetup(featureLevel)
% PsychDefaultSetup(featureLevel) - Perform standard setup for Psychtoolbox.
%
% This function performs a few typical "boilerplate" setup operations
% at the beginning of a script to avoid repetitive code at the top of
% a script.
%
% Add it at the top of your script, so all its settings affect successive
% Psychtoolbox commands.
%
% The parameter 'featureLevel' determines what kind of setup is
% specifically performed. A higher number for 'featureLevel' will
% include all setup steps and default settings for lower numbers
% of 'featureLevel' and extend on them. E.g., a featureLevel of 2 would
% imply all setup operations of featureLevel 0 and 1, plus some new
% additional setup operations.
%
% A 'featureLevel' of 0 will do nothing but execute the AssertOpenGL command,
% to make sure that the Screen() mex file is properly installed and functional.
%
% A 'featureLevel' of 1 will additionally execute KbName('UnifyKeyNames') to
% provide a consistent mapping of keyCodes to key names on all operating
% systems.
%
% A 'featureLevel' of 2 will additionally imply the execution of
% Screen('ColorRange', window, 1, [], 1); immediately after and whenever
% PsychImaging('OpenWindow',...) is called, thereby switching the default
% color range from the classic 0-255 integer number range to the normalized
% floating point number range 0.0 - 1.0 to unify color specifications
% across differently capable display output devices, e.g., standard 8 bit
% displays vs. high precision 16 bit displays. Please note that clamping of
% valid color values to the 0 - 1 range is still active and colors will
% still be represented by 256 discrete levels (8 Bit resolution), unless
% you also use PsychImaging() commands to request unclamped color
% processing or floating point precision framebuffers. This function by
% itself only changes the range, not the precision of color specifications!
%

% History:
% 22-Aug-2013  mk    Initial version written.

% Default colormode to use: 0 = clamped, 0-255 range. 1 = unclamped 0-1 range.
global psych_default_colormode;
psych_default_colormode = 0;

% Reset KbName mappings:
clear KbName;

% Define maximum supported featureLevel for this Psychtoolbox installation:
maxFeatureLevel = 2;

% Sanity check featureLevel argument:
if nargin < 1 || isempty(featureLevel) || ~isscalar(featureLevel) || ~isnumeric(featureLevel) || featureLevel < 0
    error('Mandatory featureLevel argument missing or invalid (not a scalar number or negative).');
end

% Always AssertOpenGL:
AssertOpenGL;

% Level 1+ requested?
if featureLevel >= 1
    % Unify keycode to keyname mapping across operating systems:
    KbName('UnifyKeyNames');
end

% Level 2+ requested?
if featureLevel >= 2
    % Set global environment variable to ask PsychImaging() to enable
    % normalized color range for all drawing commands and Screen('MakeTexture'):
    psych_default_colormode = 1;
end

if featureLevel > maxFeatureLevel
    error('This installation of Psychtoolbox can not execute scripts at the requested featureLevel of %i, but only up to level %i ! UpdatePsychtoolbox!', featureLevel, maxFeatureLevel);
end

return;
