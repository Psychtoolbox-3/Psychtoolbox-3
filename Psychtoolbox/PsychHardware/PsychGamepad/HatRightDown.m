function result = HatRightDown
% result = HatRightDown
% Return a constant for hat right and down position.  To read a hat position, use:
% JOYSTICK(joystickNumber,'GetHat',hatNumber).  

% 8/17/01   awi wrote it
% 10/6/05   awi Note here cosmetic changes by dgp between 7/21/03 and 10/6/05


if IsOSX 
	result = 4;
end

