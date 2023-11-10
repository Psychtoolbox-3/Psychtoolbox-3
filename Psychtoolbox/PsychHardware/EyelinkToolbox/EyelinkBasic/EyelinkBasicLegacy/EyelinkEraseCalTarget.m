function EyelinkEraseCalTarget(el, rect)
warning('EyelinkToolbox:LegacyEraseCalTarget',['The function EyelinkEraseCalTarget() is deprecated. Please update your ', ...
    'script to use the current method for handling camera setup mode callbacks with PsychEyelinkDispatchCallback.m.']);
warning('off', 'EyelinkToolbox:LegacyEraseCalTarget');
% erase calibration target
%
% USAGE: erasecaltarget(el, rect)
%
%        el: eyelink default values
%        rect: rect that will be filled with background colour
if ~IsEmptyRect(rect)
    Screen( 'FillOval', el.window, el.backgroundcolour,  rect );
    Screen( 'Flip',  el.window);
end
