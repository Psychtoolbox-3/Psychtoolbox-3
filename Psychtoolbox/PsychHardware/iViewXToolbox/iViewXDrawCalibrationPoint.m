function iviewXDrawCalibrationPoint(ivx, point)

% draw calibration point
if ~isempty(ivx.window)

    rect=CenterRectOnPoint(ivx.calRect, point(1), point(2));
    Screen('FillOval',ivx.window, ivx.calPointColour, rect);
    rect=CenterRectOnPoint(ivx.calRectCenter, point(1), point(2));
    Screen('FillOval',ivx.window, ivx.calPointColourCenter, rect);
    Screen('Flip',ivx.window);
else
    fprintf([mfilename ': cannot draw point, no window defined\n']);
end
