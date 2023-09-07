function rect=EyelinkLegacyDrawCalTarget(el, x, y)

% draw simple calibration target
%
% USAGE: rect=EyelinkLegacyDrawCalTarget(el, x, y)
%
%		el: eyelink default values
%		x,y: position at which it should be drawn
%		rect:

% simple, standard eyelink version
%   22-06-06    fwc OSX-ed

[width, height]=Screen('WindowSize', el.window);
size=round(el.calibrationtargetsize/100*width);
inset=round(el.calibrationtargetwidth/100*width);

rect=CenterRectOnPoint([0 0 size size], x, y);
Screen( 'FillOval', el.window, el.foregroundcolour,  rect );
rect=InsetRect(rect, inset, inset);
Screen( 'FillOval', el.window, el.backgroundcolour, rect );
Screen( 'Flip',  el.window);
