function ang = ATand360(y,x)
% ang = ATand360(y,x)
% returns angle from 0 to 360 based on x an y coords

% DN 2008-09-11 rewritten from scratch

ang = ATand2(y,x);

% correct for half: convert MATLAB's -180:180 to 0:360
qnegang         = ang<0;
ang(qnegang)    = ang(qnegang) + 360;
