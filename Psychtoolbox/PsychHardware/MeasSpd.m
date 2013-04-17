function [spectrum,qual] = MeasSpd(S,meterType,syncMode)
% [spectrum,qual] = MeasSpd([S],[meterType],[syncMode])
%
% This routine splines the raw return values from the
% meter to the wavelength sampling S.  The splining
% handles conversion of power units according to
% to the wavelength sampling delta.  If S is not passed,
% it is set to [380 5 81].
%
% Tries to handle low light level case gracefully by returning
% zero as the answer.
%
% meterType == 1:  PR650 (default)
% meterType == 4:  PR655
% meterType == 5:  PR670
% meterType == 6:  PR705
%
% syncMode = 'on':  Try to sync integration time with display, if meter supports it (default)
% syncMode = 'off': Don't try to sync, even if meter supports it.
%
% 9/3/93		dhb		Added default handling of S.
% 9/14/93		jms		Added global no hardware switch
% 10/1/93		dhb		Removed print on error, passed qual on up
% 10/4/93		dhb		Handle quality code 18 properly.
% 1/16/94		jms		Removed 'exist' check and declared globals.
% 2/20/94		dhb		Modified for CMETER.
% 8/11/94		dhb		Handle sync mode error condition.
% 9/7/94		dhb		Remove sync mode message.
% 11/6/96       dhb     Remove extra call to CMETER('Measure').
% 6/17/98       dhb     Add meterType switch.
% 7/1/98        dhb,jmk Fix bug in switch.
% 10/4/99       dhb,mdr Take return of -1 (timeout) to mean no light.
%               dhb,mdr Remove gHardware and gHardwareMsg globals.
% 4/4/00        dhb,jdt Pass S to PC version.
% 09/11/00      dhb   Remove colortron support.  It never worked right.
% 1/10/02       dhb,ly Make OS9 version use SERIAL.
% 2/15/02       dhb   Change name of called routine.
% 5/21/02       dgp   Tidied up code, removing superfluous COMPUTER conditional.
% 2/26/03       dhb   Change definition of PR-650 meter type to 1.
% 8/26/10       dhb   The PR-655 line called the PR-650 code.  Change to call PR-655
% 3/8/11        dhb   Pass syncMode option to speed things up for displays where it doesn't work.
% 12/06/12      zlb   Adding PR-705 support as meter type 6.

% Handle defaults
if nargin < 3 || isempty(syncMode)
    syncMode = 'on';
end
if nargin < 2 || isempty(meterType)
    meterType = 1;
end
if nargin < 1 || isempty(S)
    S = [380 5 81];
end

try
    switch meterType
        % PR-650
        case 1,
            [spectrum, qual] = PR650measspd(S,syncMode);
            
        % PR-655
        case 4,
            [spectrum, qual] = PR655measspd(S,syncMode);
            
        % PR-670
        case 5,
            [spectrum, qual] = PR670measspd(S,syncMode);
            
        % PR-705
        case 6,
            qual = 0;
            spectrum = PR705measspd(S);
            
        otherwise,
            error('Unknown meter type');
    end
catch %#ok<CTCH>
    CMClose(meterType);
    sca();
    psychrethrow(psychlasterror);
end
