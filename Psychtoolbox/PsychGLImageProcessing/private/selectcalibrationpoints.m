function [scal] = selectcalibrationpoints(scal)

commandwindow;
pause(.5);

% open window and check whether that window has the same resolution as the
% one stored in the calibration data file
[scal.windowPtr,rect] = Screen('OpenWindow', scal.screenNumber); % windowPtr: psychtoolbox reference to the active drawing window, % rect: resolution of the active drawing window in pixels [x(1) x(end) y(1) y(end)]
if sum(rect-scal.rect)~=0,
    error('The resolution of the monitor has changed since the initialization of the spatial calibration.')
end

scal = bvlSelectFitPts(scal); % select the dots that will be fitted. The dots x and y coordinates will be stored in scal.FITDOTLIST
eval(['save currentCalib[' num2str(scal.screenNumber) '].mat']);
Screen('closeall');
ShowCursor;
