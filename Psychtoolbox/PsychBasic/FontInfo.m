function FontInfo
% numFonts=FontInfo('NumFonts');        % Return the number of installed fonts
% fontInfoStructArray=FontInfo('Fonts');% Return an array of font info structs
% versionInfo=FontInfo('Version');      % Return version info for Fonts command
%
% OS X: ___________________________________________________________________
% 
% Return information about the text fonts installed on your computer.
% FontInfo may be used in conjunction with Screen('DrawText') to select a
% screen font which satifies your requirements.  FontInfo('Fonts') may omit
% information in some fields of the returned struct, depending on what
% information the font file itself supplies.  In particular, the
% verticalMetrics and horizontalMetrics substructs are often empty.
%
% The font number associated with a particular font may change if you
% restart your computer or change the installed fonts.  Therefore you 
% should not code font  numbers into your scripts.  Instead, specify fonts
% within your scripts by using font names.  You may pass font names
% directly to Screen('TextFont'), or use the table returnd by the FontInfo
% function to find the number of a font which satisfies your criteria for a
% font.  Among those critera may be the font name itself.
%
% OS 9, WINDOWS, Linux: ___________________________________________________
%
% FontInfo does not exist on OS-9, Windows and Linux.
% _________________________________________________________________________
%
% SEE ALSO: Screen, Screen('TextFont')

% HISTORY
% 7/03/04   awi wrote it.
% 7/10/04   awi cosmetic.
% 7/12/04   awi Divided into platform sections.
% 10/4/05   awi Note here cosmetic changes by dgp between 7/12/04 and 10/4/05   

% Give an error on all platforms if MATLAB fails to find the mex file and
% executes this help  file instead.
AssertMex('FontInfo.m');
