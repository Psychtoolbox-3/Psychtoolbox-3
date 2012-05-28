function wheelDelta = GetMouseWheel(mouseIndex)
% wheelDelta = GetMouseWheel([mouseIndex])
%
% Return change of mouse wheel position of a wheel mouse (in "wheel clicks")
% since last query. 'mouseIndex' is the PsychHID device index of the wheel
% mouse to query. The argument is optional: If left out, the first detected
% real wheel mouse (ie. not a trackpad) is queried.
%
% OS X: ___________________________________________________________________
%
% Uses PsychHID for low-level access to mice with mouse wheels. If wheel
% state is not queried frequent enough, the internal queue may overflow and
% some mouse wheel movements may get lost, resulting in a smaller reported
% 'wheelDelta' than the real delta since last query. On OS X 10.4.11 the
% operating system can store at most 10 discrete wheel movements before it
% discards movement events.
% _________________________________________________________________________
%
% MS-Windows and Linux: ___________________________________________________
%
% This function is not (yet?) supported.
%
% _________________________________________________________________________
% See also: GetClicks, GetMouseIndices, GetMouse, SetMouse, ShowCursor,
% HideCursor
%

% History:
% 05/31/08  mk  Initial implementation.
% 05/14/12  mk  Tweaks for more mice.

% Cache the detected index of the first "real" wheel mouse to allow for lower
% execution times:
persistent wheelMouseIndex;
if isempty(wheelMouseIndex) && ((nargin < 1) || isempty(mouseIndex))
        % On OS X we execute this script, otherwise either MATLAB found the mex file
        % file and exuted this, or else this file was exucuted and exited with
        % error on the AssertMex command above.
        
        %get the number of mouse buttons from PsychHID
        mousedices=GetMouseIndices;
        numMice = length(mousedices);
        if numMice == 0
            error('GetMouseWheel could not find any mice connected to your computer');
        end

        allHidDevices=PsychHID('Devices');
        for i=1:numMice
            b=allHidDevices(mousedices(i)).wheels;
            if ~IsOSX
                % On Non-OS/X we can't detect .wheels yet, so fake
                % 1 wheel for each detected mouse and hope for the best:
                b = 1;
            end
            
            if any(b > 0) && isempty(strfind(lower(allHidDevices(mousedices(i)).product), 'trackpad'))
                wheelMouseIndex = mousedices(i);
                break;
            end
        end
        
        if isempty(wheelMouseIndex)
            error('GetMouseWheel could not find any mice with mouse wheels connected to your computer');
        end
end;

% Override mouse index provided?
if nargin < 1
    % Nope: Assign default detected wheel-mouse index:
    mouseIndex = wheelMouseIndex;
end

% Use low-level access to get wheel state report: Refetch until empty
% report is returned:
wheelDelta = 0;
rep = PsychHID('GetReport', mouseIndex, 1, 0, 10);
while ~isempty(rep)
    wheely = rep(end);
    switch wheely
        case 1,
            wheelDelta = wheelDelta + 1;
        case 255,
            wheelDelta = wheelDelta - 1;
    end
    [rep, err] = PsychHID('GetReport', mouseIndex, 1, 0, 4);
    if err.n
        fprintf('GetMouseWheel: GetReport error 0x%s. %s: %s\n', hexstr(err.n), err.name, err.description);
    end
end

return;
