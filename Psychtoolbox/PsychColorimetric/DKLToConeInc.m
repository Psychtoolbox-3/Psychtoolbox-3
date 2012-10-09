function [coneInc] = DKLToConeInc(dkl,bg)
% [coneInc] = DKLToConeInc(dkl,bg)
%
% Convert from DKL to incremental cone coordinates.
% 
% The DKL coordinate system convention is (Lum, RG, S)
%
% The code follows that published by Brainard
% as an appendix to Human Color Vision by Kaiser
% and Boynton.
%
% See also ConeIncToDKL, ComputeDKL_M, DKLDemo.
%
% 8/30/96	dhb  Converted this from script.
% 10/5/12   dhb  Comment specifying coordinate system convention.

% Compute conversion matrix
M = ComputeDKL_M(bg);

% Multiply the vectors we wish to
% convert by M to obtain its DKL coordinates.
coneInc = inv(M)*dkl;




						 

