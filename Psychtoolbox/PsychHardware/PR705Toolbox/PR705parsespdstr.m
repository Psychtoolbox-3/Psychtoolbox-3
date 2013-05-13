function spd = PR705parsespdstr(rawspd, S)
% PR705parsespdstr - Parses the spectral power distribution string returned by the PR-705.
%
% Syntax:
% spd = PR705parsespdstr(rawspd [, S])
%
% Description:
% Parse the spectral power distribution string returned by the PR-705. The
% results are splined to the desired wavelength sampling defined by the 'S'
% input parameter.
%
% Input:
% rawspd (1xN char) - Raw data returned from a meter measurment.
% S (1x3) - Wavelength sampling.  Default: [380 5 81]
%
% 12/06/12    zlb   Wrote it based on the PR670Toolbox.

if nargin < 2 || isempty(S)
    S = [380 5 81];
end

spd = [];
if ~isempty(rawspd)
    C = textscan(rawspd, '%d,%f', 'HeaderLines', 1);
    wls = double(C{1});
    spd = double(C{2});
    
    if isempty(wls) || isempty(spd)
        return
    end
    
    if numel(wls) ~= 201 % 201 is the PR-705 default
        error('Unexpected number of wavelength samples for PR-705!');
    elseif wls(2) - wls(1) ~= 2 || ~all(diff(wls) == 2)
        error('Unexpected wavelength sampling for the PR-705!');
    end
    
    % We _don't_ multiply the amplitudes by the bin spacing (2 nm in this 
    % case) because the data given by the meter are in units, not units/nm,
    % as per the documentation.
    
    if S(3) ~= 201
        S0 = [wls(1) wls(2)-wls(1) length(wls)];
        spd = SplineSpd(S0, spd, S);
    end
end
