function runindummymode=EyelinkDummyModeDlg
% Puts up a dialog asking to run in dummy mode.

%
%	HISTORY
%	27-11-02 fwc	created it as screen windows does not support 'dialog' function
%   280606   fwc    changed name for OSX version
%   10-04-09 mk     Made portable to old Matlabs. De-Uglified.

if exist('questdlg') %#ok<EXIST>
    drawnow;
    ButtonName=questdlg('Run in dummy mode?', ...
        'Eyelink not connected', ...
        'Yes','Quit','Quit');
else
    % Use keyboard query instead:
    ButtonName = 'Quit';
    myanswer = lower(input('Eyelink is not connected. Run in dummy mode? [yes/no] ', 's'));
    if ~isempty(myanswer)
        if myanswer(1) == 'y'
            ButtonName = 'Yes';
        end
    end
end

switch ButtonName,
    case 'Yes',
        disp('Running program in dummy mode');
        runindummymode=1;
    case 'Quit',
        disp('Goodbye......')
        runindummymode=0;
end % switch

return;
