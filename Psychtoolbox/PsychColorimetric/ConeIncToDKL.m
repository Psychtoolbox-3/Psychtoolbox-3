function [dkl] = ConeIncToDKL(coneInc,bg)
% [dkl] = ConeIncToDKL(coneInc,bg)
%
% Convert from incremental cone coordinates to DKL
% coordinates.
%
% The DKL coordinate system convention is (Lum, RG, S)
%
% The code follows that published by Brainard
% as an appendix to Human Color Vision by Kaiser
% and Boynton.
%
% See also DKLToConeInc, ComputeDKL_M, DKLDemo.
%
% 8/30/96	dhb	 Converted this from script.
% 10/5/12   dhb  Comment specifying coordinate system convention.

% Compute conversion matrix
M = ComputeDKL_M(bg);

% Multiply the vectors we wish to
% convert by M to obtain its DKL coordinates.
dkl = M*coneInc;




						 

