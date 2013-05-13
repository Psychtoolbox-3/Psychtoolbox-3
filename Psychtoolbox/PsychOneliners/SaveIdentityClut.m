function LUT = SaveIdentityClut(windowPtr, LUT)
% savedClut = SaveIdentityClut(windowPtr [, LUT])
%
% This routine defines a LUT or LUT-type for use as identity gamma lookup
% table for applications that need such a table, e.g., Bits+ box,
% VideoSwitcher, Video attenuators, the BrightSide HDR display etc.
%
% It writes the LUT into a config file. The function LoadIdentityClut()
% will use the written LUT from that config file, if such a file exists.
%
% You can either pass a full blown LUT in the optional argument 'LUT', or
% you can pass a LUT id code (see codes and corresponding LUT's in the
% source code of LoadIdentityClut.m). If you omit the argument, then the
% function will readout and store the current CLUT of the given display -
% in the hope that the current LUT is actually an identity CLUT.
%
% windowPtr is optional: It can be omitted, in which case the LUT of the
% screen with screenid = max(Screen('Screens')) is used. It can be a
% screenid for the screen to read out, or it can be the window handle of an
% existing open onscreen window, in which case the screen corresponding to
% that window is used.
%
% In any case the type of operating system, graphics card vendor and model,
% as well as OpenGL version and OpenGL driver version and screenid is used
% to specify the configuration to which the LUT applies. This to
% disambiguate in case multiple different operating systems, os versions,
% GPU's or display types are used.
%
% The routine returns the stored clut in the optional return argument 'oldClut'.
%

% History:
% 09/20/09   mk  Written.
% 11/21/09   mk  gpuName can contain '/', replace them as well with '_'.
%                Merged bugfix by Allan Rempel.

closeWindow = 0;

if nargin < 1
    windowPtr = [];
end

% No windowPtr? Open our own on max'imum screen and use that window/screen:
if isempty(windowPtr)
    windowPtr = Screen('OpenWindow', max(Screen('Screens')));
    closeWindow = 1;
end

% No onscreen window, but probably a screenid? Open window on that screen:
if Screen('Windowkind', windowPtr) ~= 1
    if ismember(windowPtr, Screen('Screens'))
        windowPtr = Screen('OpenWindow', windowPtr);
        closeWindow = 1;
    else
        error('Provided "windowPtr" is neither a windowhandle, nor a screenid!');
    end
else
    % Onscreen window handle. Good, just use that...
end

if nargin < 2
    LUT = [];
end

% Empty or missing LUT argument?
if isempty(LUT)
    % Yes. Readout the gamma table of our windowPtr's associated screen:
    LUT = Screen('ReadNormalizedGammatable', windowPtr);
else
    % Sanity check:
    if ~( isnumeric(LUT) && (isscalar(LUT) || ((size(LUT,1) >= 1) && (size(LUT,2) == 3))) )
        sca;
        error('LoadIdentityClut: Loaded data from config file is not a valid LUT! Not a numeric matrix or less than 1 row, or not 3 columns!');
    end
    
    % Passed.
end

% LUT contains either a scalar LUT type id code for LoadIdentityCLUT, or a
% full blown properly formatted lookup table. Find the name under which
% this LUT should be saved to the file:

% Get screen id for this window:
screenid = Screen('WindowScreenNumber', windowPtr);

% Query vendor of associated graphics hardware:
winfo = Screen('GetWindowInfo', windowPtr);

% Raw renderer string, with leading or trailing whitespace trimmed:
gpuname = strtrim(winfo.GLRenderer);

% Replace blanks and '/' with underscores:
gpuname(isspace(gpuname)) = '_';
gpuname = regexprep( gpuname, '/', '_' );

% Same game for version string:
glversion = strtrim(winfo.GLVersion);

% Replace blanks with underscores:
glversion(isspace(glversion)) = '_';
glversion(glversion == '.') = '_';

% Close our onscreen window, if needed:
if closeWindow
    Screen('Close', windowPtr);
end

% Is there a predefined configuration file for the proper type of LUT to
% load? Build kind'a unique path and filename for our system config:
lpath = sprintf('%sIdentityLUT_%s_%s_%s_Screen%i.mat', PsychtoolboxConfigDir, OSName, gpuname, glversion, screenid);
fprintf('SaveIdentityClut: Trying to save LUT for screen %i to file %s.\n', screenid, lpath);

lutconfig = LUT; %#ok<NASGU>
save(lpath, 'lutconfig', '-mat', '-V6');
fprintf('Done. Will be used at next invocation of LoadIdentityClut for this display.\n\n');

return;
