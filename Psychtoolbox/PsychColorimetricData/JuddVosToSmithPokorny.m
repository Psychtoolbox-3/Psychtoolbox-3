function [T_cones_sp,S_cones_sp,T_cones_sp_unnormalized,M_juddVosToConesSp] = JuddVosToSmithPokorny
% [T_cones_sp,S_cones_sp,M_juddVosToConesSp] = JuddVosToSmithPokorny
%
% Load in PTB's Judd-Vos XYZ cmfs and convert to Smith-Pokorny foveal cone
% fundamentals according to the matrix provided on CVRL.  
%
% The returned fundamentals are normalized to a peak of 1, to match
% longtime conventions within PTB. Unnormalized versions are also returned.
%
% The normalized versions have now been stored in PTB data file T_cones_sp,
% with the original PTB versions saved in T_cones_sp_orig.
%
% See notes on cvrl.org about slight discrepancies between various
% tabulations of the Smith-Pokorny fundamentals.
%
% Thanks to Danny Garside for pointing out that we should have tabulated
% functions that extend beyond 400-700 nm.
%
% 7/21/18  dhb  Wrote it.

%% Examples:
%{
% Compare with PTB original sp data
[T_cones_sp,S_cones_sp,M_juddVosToConesSp] = JuddVosToSmithPokorny;
load T_cones_sp_orig
figure; hold on
plot(SToWls(S_cones_sp),T_cones_sp');
plot(SToWls(S_cones_sp_orig),T_cones_sp_orig','o');
%}

%% Load in PTB Judd-Vos
load T_xyzJuddVos;

%% Type in the conversion matrix
M_juddVosToConesSp = [...
    0.15514 0.54312 -0.03286 ; ...
   -0.15514 0.45684  0.03286 ; ...
    0        0       0.00801];

%% Transform
T_cones_sp_unnormalized = M_juddVosToConesSp*T_xyzJuddVos;
S_cones_sp = S_xyzJuddVos;

%% Normalize
for ii = 1:3
    T_cones_sp(ii,:) = T_cones_sp_unnormalized(ii,:)/max(T_cones_sp_unnormalized(ii,:));
end

