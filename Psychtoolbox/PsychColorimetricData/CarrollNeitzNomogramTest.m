function maxAbsDiff = CarrollNeitzNomogramTest
% maxAbsDiff = CarrollNeitzNomogramTest
%
% Validate CarrollNeitzNomogram against the original Neitz-lab spectsens.m.
%
% The PTB nomogram returns absorbance normalized to a peak of 1. The
% original spectsens helper returns absorptance with optical density and
% a normalization by the peak absorptance, 1-10^-OD. To compare like with
% like, we:
%   1) compute PTB absorbance with CarrollNeitzNomogram;
%   2) convert to absorptance with AbsorbanceToAbsorptance;
%   3) normalize by the peak absorptance to match spectsens;
%   4) normalize the sampled spectsens result to peak 1, matching the PTB
%      convention applied inside CarrollNeitzNomogram when S(2) <= 1.
%
% The test plots one panel each for the L, M, and S cones and reports the
% maximum absolute difference across all wavelengths and cone classes.
%
% History:
%   2026-04-10  dhb    Wrote it.
%   2026-04-10  Codex  Implement validation against original spectsens.

S = [380 1 321];
wls = SToWls(S);
lambdaMaxLMS = [557.5; 530.0; 420.0];
opticalDensity = 0.30;
coneNames = {'L', 'M', 'S'};

% We do not expect bitwise identity because PTB normalizes the absorbance
% before the OD transform, while spectsens normalizes the OD-adjusted
% absorptance. The resulting discrepancy is very small on a 1 nm grid.
tolerance = 1e-4;

T_absorbance = CarrollNeitzNomogram(S, lambdaMaxLMS);
T_absorptance = AbsorbanceToAbsorptance(T_absorbance, S, opticalDensity * ones(size(lambdaMaxLMS)));
peakAbsorptance = 1 - 10^(-opticalDensity);
T_absorptance = T_absorptance / peakAbsorptance;

% The original spectsens routine normalizes its absorptance output by the peak
% value expected from the OD (peakAbsorptance above), which
% is a little odd to my way of thinking.  But given this, we just compare
% the normalized outputs. However, that normalization doesn't lead to a
% max of 1 because the value that comes back does not have a peak of 1.
% It's close, but not exactly.  So we normalize by the max so we are
% comparing relative shapes here.
T_spectsens = zeros(size(T_absorptance));
for ii = 1:length(lambdaMaxLMS)
    T_spectsens(ii, :) = spectsens( ...
        'lambda', lambdaMaxLMS(ii), ...
        'OD', opticalDensity, ...
        'output', 'alog', ...
        'spectrum', wls);
    T_spectsens(ii, :) = T_spectsens(ii, :) / max(T_spectsens(ii, :));
end

diffMatrix = T_absorptance - T_spectsens;
maxAbsDiff = max(abs(diffMatrix(:)));

fprintf('CarrollNeitzNomogramTest max abs diff: %.3g\n', maxAbsDiff);
if maxAbsDiff > tolerance
    error('CarrollNeitzNomogramTest:Mismatch', ...
        'PTB CarrollNeitzNomogram does not match spectsens (max diff %.3g, tol %.3g).', ...
        maxAbsDiff, tolerance);
end

figure; clf;
set(gcf, 'Position', [100 100 900 900]);
for ii = 1:length(lambdaMaxLMS)
    subplot(3, 1, ii); hold on;
    plot(wls, T_spectsens(ii, :), 'k', 'LineWidth', 3);
    plot(wls, T_absorptance(ii, :), 'r--', 'LineWidth', 1.5);
    xlabel('Wavelength (nm)', 'FontSize', 11);
    ylabel('Absorptance', 'FontSize', 11);
    title(sprintf('%s cone (\\lambda_{max} = %.1f nm)', coneNames{ii}, lambdaMaxLMS(ii)), 'FontSize', 12);
    xlim([wls(1) wls(end)]);
    ylim([0 1.05]);
    legend({'spectsens', 'CarrollNeitzNomogram + AbsorbanceToAbsorptance'}, ...
        'Location', 'NorthEastOutside');
end

end

function [withOD, extinction] = spectsens(varargin)
% spectsens returns a photopigment spectral sensitivity curve as defined by
% Carroll, McMahon, Neitz, and Neitz.
%
% This is the original code from the Neitz lab website, preserved here as
% a local helper for validation of CarrollNeitzNomogram. The only edits are
% removal of `lower` calls on numeric inputs so the helper runs in MATLAB.

defaultLambdaMax = 559;
defaultOD = 0.30;
defaultOutput = 'alog';
defaultSpectrum = 380:700;

ip = inputParser();
ip.CaseSensitive = false;
addParameter(ip, 'lambda', defaultLambdaMax, @(x) validateattributes(x, {'numeric'}, {'>=', 0}));
addParameter(ip, 'OD', defaultOD, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0.001}));
addParameter(ip, 'output', defaultOutput, @(x) ismember(lower(x), {'log', 'alog'}));
addParameter(ip, 'spectrum', defaultSpectrum, @(x) validateattributes(x, {'numeric'}, {'integer', '>=', 0}));
parse(ip, varargin{:});

LambdaMax = ip.Results.lambda;
OpticalDensity = ip.Results.OD;
Output = ip.Results.output;
Spectrum = ip.Results.spectrum;

A = 0.417050601;
B = 0.002072146;
C = 0.000163888;
D = -1.922880605;
E = -16.05774461;
F = 0.001575426;
G = 5.11376E-05;
H = 0.00157981;
I = 6.58428E-05;
J = 6.68402E-05;
K = 0.002310442;
L = 7.31313E-05;
M = 1.86269E-05;
N = 0.002008124;
O = 5.40717E-05;
P = 5.14736E-06;
Q = 0.001455413;
R = 4.217640000E-05;
S = 4.800000000E-06;
T = 0.001809022;
U = 3.86677000E-05;
V = 2.99000000E-05;
W = 0.001757315;
X = 1.47344000E-05;
Y = 1.51000000E-05;
Z = OpticalDensity + 0.00000001;

A2 = log10(1.00000000 / LambdaMax) - log10(1.00000000 / 558.5);
vector = log10(Spectrum .^ -1);
const = 1 / sqrt(2 * pi);

exTemp1 = log10(-E + E * tanh(-((10 .^ (vector - A2)) - F) / G)) + D;
exTemp2 = A * tanh(-((10 .^ (vector - A2)) - B) / C);
exTemp3 = -(J / I * (const * exp(1) .^ (-0.5 * (((10 .^ (vector - A2)) - H) / I) .^ 2)));
exTemp4 = -(M / L * (const * exp(1) .^ (-0.5 * (((10 .^ (vector - A2)) - K) / L) .^ 2)));
exTemp5 = -(P / O * (const * exp(1) .^ (-0.5 * (((10 .^ (vector - A2)) - N) / O) .^ 2)));
exTemp6 = (S / R * (const * exp(1) .^ (-0.5 * (((10 .^ (vector - A2)) - Q) / R) .^ 2)));
exTemp7 = ((V / U * (const * exp(1) .^ (-0.5 * (((10 .^ (vector - A2)) - T) / U) .^ 2))) / 10);
exTemp8 = ((Y / X * (const * exp(1) .^ (-0.5 * (((10 .^ (vector - A2)) - W) / X) .^ 2))) / 100);
exTemp = exTemp1 + exTemp2 + exTemp3 + exTemp4 + exTemp5 + exTemp6 + exTemp7 + exTemp8;

ODTemp = log10((1 - 10 .^ (-((10 .^ exTemp) * Z))) / (1 - 10 ^ -Z));

if strcmp(Output, 'log')
    extinction = exTemp;
    withOD = ODTemp;
else
    extinction = 10 .^ exTemp;
    withOD = 10 .^ ODTemp;
end
end
