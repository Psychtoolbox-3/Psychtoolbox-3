function EyelinkLegacyEraseCalTarget(el, rect)

% erase calibration target
%
% USAGE: erasecaltarget(el, rect)
%
%		el: eyelink default values
%		rect: rect that will be filled with background colour 
if ~IsEmptyRect(rect)
    Screen( 'FillOval', el.window, el.backgroundcolour,  rect );
    Screen( 'Flip',  el.window);
end
