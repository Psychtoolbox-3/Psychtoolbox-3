function result = HatLeftUp
% result = HatLeftUp
% Return a constant for hat left and up position.  To read a hat position, use:
% JOYSTICK(joystickNumber,'GetHat',hatNumber).  

% 8/17/01   awi wrote it
% 7/21/03   awi added switch and OS X condition.
% 10/6/05   awi Note here cosmetic changes by dgp between 7/21/03 and 10/6/05


if IsOSX 
	result = 8;
end
