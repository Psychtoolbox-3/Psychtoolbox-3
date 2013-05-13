function [result, ivx]=iViewXGetCalibrationPoints(ivx)

result=-1;
% could we skip this if we are sure we have the correct info?
%     ET_CSZ 1024 768   ET_PNT 1 512 384   ET_PNT 2 51 38   ET_PNT 3 973 730

% we fill table ivx.absCalPos with supplied calibration point coords
ivx.calScreenSize=[1024 768];
%

% print warning if the received calibration screen size does not
% match we the present screen resolution

%     If the screen resolution of the subject PC is not the same
%     as set in iView properties, the remote calibration  program
%     will have to store the values of the calibration area size to
%     do calculations and all points will have to  be stored to know
%     where to show them.  Alternately, be sure that the resolution
%     set in iView Calibration  properties matches that of the
%     Subject PC.

if ~isempty(ivx.window)
    [w, h]=WindowSize(ivx.window);
    %fprintf([mfilename ' 1\n']);

    if ivx.calScreenSize(1)~=w || ivx.calScreenSize(2)~=h
        fprintf('Warning: actual screen size does not match calibration settings\n');
        fprintf('Actual: width: % 5d, Height: % 5d\n', w, h);
        fprintf('Calibr: width: % 5d, Height: % 5d\n', ivx.calScreenSize(1), ivx.calScreenSize(2));

        % we can here adjust the points to match the actual screen size
    end
    %fprintf([mfilename ' 2\n']);
end

isfield(ivx, 'absCalPos')
isempty(ivx.absCalPos)
isempty(ivx.absCalPos)
size(ivx.absCalPos,1)<1
size(ivx.absCalPos,2)~=2

% check if the position table exists and is of correct size
if ~isfield(ivx, 'absCalPos') || isempty(ivx.absCalPos) || size(ivx.absCalPos,1)<1 || size(ivx.absCalPos,2)~=2
    txt=[ mfilename ': something is wrong with calibration point pixel values.'];
    error(txt);
end


%fprintf([mfilename ' 3\n']);
result=1;
