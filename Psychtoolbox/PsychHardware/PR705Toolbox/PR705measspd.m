function spd = PR705measspd(S)
% PR705measspd - Get a spectral power distribution from the PR-705.
%
% Syntax:
% spd = PR705measspd([S])
%
% Input:
% S (1x3) - Desired wavelength sampling. Default: [380 5 81]
%
% Output:
% spd (Nx1) - The spectral power distribution (N = S(3) or 81 by default).
%
% 12/06/12   zlb   Wrote it loosely based on the PR670Toolbox.

% PR-705 measurement error codes
%   -5000: Weak Signal
%   -4999: Time underflow, level overflow
%   -4996: A2D overflow, measuring light
%   -4995: A2D overflow, measuring dark
%   -4994: Variable light level
%   -4993: Adaptive exposure time limit

if nargin < 1 || isempty(S)
    S = [380 5 81];
end

timeout = 300;
[rawspd,errcode] = PR705rawspd(timeout);

if isempty(rawspd) || isempty(errcode)
    error('Didn''t get a measurement reponse after %d seconds', timeout);
end

switch errcode
    case 0
        spd = PR705parsespdstr(rawspd, S);
    case {-4995,-4993} % too dark?
        spd = zeros(S(3), 1);
    otherwise
        error('Received unhandled error code %d', errcode);
end
