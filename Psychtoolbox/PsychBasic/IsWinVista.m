function isVista = IsWinVista
% IsWinVista - Return if this is a MS-Windows Vista OS or later.
%

% History:
% 23.10.2012  mk  Written.

persistent isThisVista;
if isempty(isThisVista)
    if ~IsWin
        isThisVista = 0;
    else
        c = Screen('Computer');
        isThisVista = c.IsVistaClass;
    end
end

isVista = isThisVista;

return;
