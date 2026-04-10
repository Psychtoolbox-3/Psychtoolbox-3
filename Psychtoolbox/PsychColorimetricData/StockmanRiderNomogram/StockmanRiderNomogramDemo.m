%% srDemo
%
% Demonstrates use of Stockman-Rider cone fundamentals formuale,
% as implemented here.  The model is described in:
%     Stockman, A. & Rider, A. T. (2023).  Formulae for generating standard and
%     individual human cone spectral sensitivities.  Color Research and Application,
%     48, 818-840.
%
% This code started from the python version provided by the authors,
% was converted to Matlab by Claude AI, and then modified and tested by
% Brainard.  The python version is available on github at this repository:
%    https://github.com/CVRL-IoO/Individual-CMFs.git
% The python code gives the constants to more places than the paper; I am
% not sure which version should be considered definitive.
% 
% Claude's translation was good but not perfect, so lots of hand checking.
%
% The python version used to produce this code was downloaded from github
% on September 2, 2025.

% History:
%   2025-09-05  dhb  Matlab first version as described above.

%% Clear
clear; close all;

%% Example 0: Standard 2-degree color matching functions
fprintf('Example 0: CIE absorbances from Fourier polynomial\n');
fprintf('===================================\n');

% Load PTB SS absorbance
load T_log10coneabsorbance_ss

% Compute LMS absorbances
%
% Absorbance peaks should be 551.9, 529.8 and 416.9 according to the paper.
% The peaks only do anyting in the code if you are shifting the absorbance
% spectra, as they determine how the argument to the Fourier polynomial gets shifted.
Lshift = 0.0; 
Mshift = 0.0;
Sshift = 0.0;
deltaWlnm = 0.1;
wls = (360:deltaWlnm:850)';
SR_LMSabsorbances = srLMSconelog(wls,Lshift,Mshift,Sshift,'lin');
fprintf('Wavelength range: %.0f - %.0f nm\n', min(SR_LMSabsorbances(:,1)), max(SR_LMSabsorbances(:,1)));
fprintf('Number of data points: %d\n', size(SR_LMSabsorbances, 1));
fprintf('Peak L absorbance at: %.1f nm\n', SR_LMSabsorbances(find(SR_LMSabsorbances(:,2) == max(SR_LMSabsorbances(:,2)), 1), 1));
fprintf('Peak M absorbance at: %.1f nm\n', SR_LMSabsorbances(find(SR_LMSabsorbances(:,3) == max(SR_LMSabsorbances(:,3)), 1), 1));
fprintf('Peak S absorbance at: %.1f nm\n', SR_LMSabsorbances(find(SR_LMSabsorbances(:,4) == max(SR_LMSabsorbances(:,4)), 1), 1));

% Plot.  Should look like Figure 1 of the paper.
%
% Note, however, that the PTB tabulated data for S cone absorbance were extended by someone
% by linear extrapolation beyond the official cutoff of 615 nm (see note in PTB Contents.m).
% This does not match the extrapolation in the Stockman Rider paper, as
% they (see paper) used Lamb's photopigment template for this purpose.
% The differences only become apparent in the plots 7 log units down from the peak,
% so not something to spend much time on.
%
% The small deviations you can see near the peaks in the linear plot are also there in the paper
% figure, as are the deviations at the short wavelength end of the L and M cone absorbances.
absorbanceFig = figure;
set(gcf,'Position',[100 100 800 1400]);
subplot(2,1,1); hold on;
plot(SToWls(S_log10coneabsorbance_ss),T_log10coneabsorbance_ss(1,:)','r','LineWidth',4);
plot(SToWls(S_log10coneabsorbance_ss),T_log10coneabsorbance_ss(2,:)','g','LineWidth',4);
plot(SToWls(S_log10coneabsorbance_ss),T_log10coneabsorbance_ss(3,:)','b','LineWidth',4);
plot(SR_LMSabsorbances(:,1), log10(SR_LMSabsorbances(:,2)), 'y--', 'LineWidth', 2);
plot(SR_LMSabsorbances(:,1), log10(SR_LMSabsorbances(:,3)), 'y--', 'LineWidth', 2);
plot(SR_LMSabsorbances(:,1), log10(SR_LMSabsorbances(:,4)), 'y--', 'LineWidth', 2);
xlim([350 850]);
ylim([-8 0]);
xlabel('Wavelength (nm)');
ylabel('Log10 Sensitivity');
title('LMS Absorbance');
legend('L', 'M', 'S', 'Location', 'best');
set(gca,'Color',[0.4 0.4 0.4]);

subplot(2,1,2); hold on;
plot(SToWls(S_log10coneabsorbance_ss),10.^(T_log10coneabsorbance_ss(1,:))','r','LineWidth',4);
plot(SToWls(S_log10coneabsorbance_ss),10.^(T_log10coneabsorbance_ss(2,:))','g','LineWidth',4);
plot(SToWls(S_log10coneabsorbance_ss),10.^(T_log10coneabsorbance_ss(3,:))','b','LineWidth',4);
plot(SR_LMSabsorbances(:,1), SR_LMSabsorbances(:,2), 'y--', 'LineWidth', 2);
plot(SR_LMSabsorbances(:,1), SR_LMSabsorbances(:,3), 'y--', 'LineWidth', 2);
plot(SR_LMSabsorbances(:,1), SR_LMSabsorbances(:,4), 'y--', 'LineWidth', 2);
xlim([350 850]);
xlabel('Wavelength (nm)');
ylabel('Sensitivity');
title('LMS Absorbance');
legend('L', 'M', 'S', 'Location', 'best');
set(gca,'Color',[0.4 0.4 0.4]);

%% Example 1: Calculate standard fundamentals with no shifts (normal observer)
fprintf('Example 1: Standard CIE Fundamentals\n');
fprintf('===================================\n');

% 2 deg fundamentals
Lshift = 0.0; 
Mshift = 0.0;
Lod_2deg = 0.50;
Mod_2deg = 0.50;
Sod_2deg = 0.40;
mac_460_2deg = 0.350; 
lens_400 = 1.7649;
[SR_LMS_energyfundamentals_2deg, LMS_quantalfundementals_2deg] = srCalculateCMFs(deltaWlnm,Lshift,Mshift,Lod_2deg,Mod_2deg,Sod_2deg,mac_460_2deg,lens_400);

% For 10-deg, macular pigment value is 0.095 = 0.271*0.350, lens stays the same, and photopigment ODs are 0.38, 0.38 and 0.30
Lod_10deg = 0.38;
Mod_10deg = 0.38;
Sod_10deg = 0.30;
mac_460_10deg = 0.095; 
[SR_LMS_energyfundamentals_10deg, LMS_quantalfundametals_10deg] = srCalculateCMFs(deltaWlnm,Lshift,Mshift,Lod_10deg,Mod_10deg,Sod_10deg,mac_460_10deg,lens_400);

% Display some key information
fprintf('Wavelength range: %.0f - %.0f nm\n', min(SR_LMS_energyfundamentals_2deg(:,1)), max(SR_LMS_energyfundamentals_2deg(:,1)));
fprintf('Number of data points: %d\n', size(SR_LMS_energyfundamentals_2deg, 1));
fprintf('Peak L energy fundamental 2 deg: %.1f nm\n', SR_LMS_energyfundamentals_2deg(find(SR_LMS_energyfundamentals_2deg(:,2) == max(SR_LMS_energyfundamentals_2deg(:,2)), 1), 1));
fprintf('Peak M energy fundamental 2 deg: %.1f nm\n', SR_LMS_energyfundamentals_2deg(find(SR_LMS_energyfundamentals_2deg(:,3) == max(SR_LMS_energyfundamentals_2deg(:,3)), 1), 1));
fprintf('Peak S energy fundamental 2 deg: %.1f nm\n', SR_LMS_energyfundamentals_2deg(find(SR_LMS_energyfundamentals_2deg(:,4) == max(SR_LMS_energyfundamentals_2deg(:,4)), 1), 1));
fprintf('Peak L energy fundamental 10 deg: %.1f nm\n', SR_LMS_energyfundamentals_10deg(find(SR_LMS_energyfundamentals_10deg(:,2) == max(SR_LMS_energyfundamentals_10deg(:,2)), 1), 1));
fprintf('Peak M energy fundamental 10 deg: %.1f nm\n', SR_LMS_energyfundamentals_10deg(find(SR_LMS_energyfundamentals_10deg(:,3) == max(SR_LMS_energyfundamentals_10deg(:,3)), 1), 1));
fprintf('Peak S energy fundamental 10 deg: %.1f nm\n', SR_LMS_energyfundamentals_10deg(find(SR_LMS_energyfundamentals_10deg(:,4) == max(SR_LMS_energyfundamentals_10deg(:,4)), 1), 1));

% Load PTB versions of the 2-deg and 10-deg energy fundamentals
load T_cones_ss2_extended
load T_cones_ss10

% Plot the results.  Should match Figure 3 of paper, and does by eye.
fundamentalFig = figure;
set(gcf,'Position',[100 100 1600 1400]);
subplot(2,2,1); hold on;
plot(SToWls(S_cones_ss2_extended),log10(T_cones_ss2_extended_linear(1,:))','r','LineWidth',4);
plot(SToWls(S_cones_ss2_extended),log10(T_cones_ss2_extended_linear(2,:))','g','LineWidth',4);
plot(SToWls(S_cones_ss2_extended),log10(T_cones_ss2_extended_linear(3,:))','b','LineWidth',4);
plot(SR_LMS_energyfundamentals_2deg(:,1), log10(SR_LMS_energyfundamentals_2deg(:,2)), 'y--', 'LineWidth', 2);
plot(SR_LMS_energyfundamentals_2deg(:,1), log10(SR_LMS_energyfundamentals_2deg(:,3)), 'y--', 'LineWidth', 2);
plot(SR_LMS_energyfundamentals_2deg(:,1), log10(SR_LMS_energyfundamentals_2deg(:,4)), 'y--', 'LineWidth', 2);
xlim([350 850]);
xlabel('Wavelength (nm)');
ylabel('Log10 Sensitivity');
title('LMS 2-deg Energy Fundamentals');
legend('L', 'M', 'S', 'Location', 'best');
set(gca,'Color',[0.4 0.4 0.4]);

subplot(2,2,3); hold on;
plot(SToWls(S_cones_ss2_extended),T_cones_ss2_extended_linear(1,:)','r','LineWidth',4);
plot(SToWls(S_cones_ss2_extended),T_cones_ss2_extended_linear(2,:)','g','LineWidth',4);
plot(SToWls(S_cones_ss2_extended),T_cones_ss2_extended_linear(3,:)','b','LineWidth',4);
plot(SR_LMS_energyfundamentals_2deg(:,1), SR_LMS_energyfundamentals_2deg(:,2), 'y--', 'LineWidth', 2);
plot(SR_LMS_energyfundamentals_2deg(:,1), SR_LMS_energyfundamentals_2deg(:,3), 'y--', 'LineWidth', 2);
plot(SR_LMS_energyfundamentals_2deg(:,1), SR_LMS_energyfundamentals_2deg(:,4), 'y--', 'LineWidth', 2);xlabel('Wavelength (nm)');
xlim([350 850]);
ylabel('Sensitivity');
title('LMS 2-deg Energy Fundamentals');
legend('L', 'M', 'S', 'Location', 'best');
set(gca,'Color',[0.4 0.4 0.4]);

subplot(2,2,2); hold on;
plot(SToWls(S_cones_ss10),log10(T_cones_ss10(1,:))','r','LineWidth',4);
plot(SToWls(S_cones_ss10),log10(T_cones_ss10(2,:))','g','LineWidth',4);
plot(SToWls(S_cones_ss10),log10(T_cones_ss10(3,:))','b','LineWidth',4);
plot(SR_LMS_energyfundamentals_10deg(:,1), log10(SR_LMS_energyfundamentals_10deg(:,2)), 'y--', 'LineWidth', 2);
plot(SR_LMS_energyfundamentals_10deg(:,1), log10(SR_LMS_energyfundamentals_10deg(:,3)), 'y--', 'LineWidth', 2);
plot(SR_LMS_energyfundamentals_10deg(:,1), log10(SR_LMS_energyfundamentals_10deg(:,4)), 'y--', 'LineWidth', 2);
xlim([350 850]);
xlabel('Wavelength (nm)');
ylabel('Log10 Sensitivity');
title('LMS 10-deg Energy Fundamentals');
legend('L', 'M', 'S', 'Location', 'best');
set(gca,'Color',[0.4 0.4 0.4]);

subplot(2,2,4); hold on;
plot(SToWls(S_cones_ss10),T_cones_ss10(1,:)','r','LineWidth',4);
plot(SToWls(S_cones_ss10),T_cones_ss10(2,:)','g','LineWidth',4);
plot(SToWls(S_cones_ss10),T_cones_ss10(3,:)','b','LineWidth',4);
plot(SR_LMS_energyfundamentals_10deg(:,1), SR_LMS_energyfundamentals_10deg(:,2), 'y--', 'LineWidth', 2);
plot(SR_LMS_energyfundamentals_10deg(:,1), SR_LMS_energyfundamentals_10deg(:,3), 'y--', 'LineWidth', 2);
plot(SR_LMS_energyfundamentals_10deg(:,1), SR_LMS_energyfundamentals_10deg(:,4), 'y--', 'LineWidth', 2);
xlim([350 850]);
xlabel('Wavelength (nm)');
ylabel('Sensitivity');
title('LMS 10-deg EnergyFundamentals');
legend('L', 'M', 'S', 'Location', 'best');
set(gca,'Color',[0.4 0.4 0.4]);

%% Example 2: Shifted L and M cones (simulating genetic variation)
% fprintf('\nExample 2: Shifted L and M cone sensitivities\n');
% fprintf('=============================================\n');
% 
% % Simulate individual with shifted cone sensitivities
% Lshift = 2.0;  % 2 nm shift in L cone
% Mshift = -1.0; % 1 nm shift in M cone (opposite direction)
% 
% [LMS_shifted, ~, RGB_shifted] = srCalculateCMFs(1.0, Lshift, Mshift);
% 
% % Compare with normal
% subplot(2,2,3);
% plot(SR_LMS_energyfundamentals_2deg(:,1), SR_LMS_energyfundamentals_2deg(:,2), 'r--', 'LineWidth', 1.5); hold on;
% plot(SR_LMS_energyfundamentals_2deg(:,1), SR_LMS_energyfundamentals_2deg(:,3), 'g--', 'LineWidth', 1.5);
% plot(LMS_shifted(:,1), LMS_shifted(:,2), 'r-', 'LineWidth', 2);
% plot(LMS_shifted(:,1), LMS_shifted(:,3), 'g-', 'LineWidth', 2);
% xlabel('Wavelength (nm)');
% ylabel('Sensitivity');
% title('Normal vs Shifted L&M Cones');
% legend('L normal', 'M normal', 'L shifted', 'M shifted', 'Location', 'best');
% grid on;
% 
% fprintf('L cone shift: %.1f nm\n', Lshift);
% fprintf('M cone shift: %.1f nm\n', Mshift);

%% Example 3: Different optical densities (simulating age effects)
% fprintf('\nExample 3: Age-related changes in optical densities\n');
% fprintf('==================================================\n');
% 
% % Young observer (higher optical densities)
% Lod_young = 0.50;
% Mod_young = 0.50;
% Sod_young = 0.40;
% mac_young = 0.35;
% lens_young = 1.76;
% 
% % Older observer (lower optical densities, more lens pigment)
% Lod_old = 0.35;
% Mod_old = 0.35;
% Sod_old = 0.28;
% mac_old = 0.25;
% lens_old = 2.5;
% 
% [LMS_young, ~, ~] = srCalculateCMFs(1.0, 0, 0, Lod_young, Mod_young, Sod_young, mac_young, lens_young);
% [LMS_old, ~, ~] = srCalculateCMFs(1.0, 0, 0, Lod_old, Mod_old, Sod_old, mac_old, lens_old);
% 
% subplot(2,2,4);
% plot(LMS_young(:,1), LMS_young(:,2), 'r-', 'LineWidth', 2); hold on;
% plot(LMS_young(:,1), LMS_young(:,3), 'g-', 'LineWidth', 2);
% plot(LMS_young(:,1), LMS_young(:,4), 'b-', 'LineWidth', 2);
% plot(LMS_old(:,1), LMS_old(:,2), 'r--', 'LineWidth', 1.5);
% plot(LMS_old(:,1), LMS_old(:,3), 'g--', 'LineWidth', 1.5);
% plot(LMS_old(:,1), LMS_old(:,4), 'b--', 'LineWidth', 1.5);
% xlabel('Wavelength (nm)');
% ylabel('Sensitivity');
% title('Young vs Older Observer');
% legend('L young', 'M young', 'S young', 'L old', 'M old', 'S old', 'Location', 'best');
% grid on;
% 
% fprintf('Young observer ODs: L=%.2f, M=%.2f, S=%.2f\n', Lod_young, Mod_young, Sod_young);
% fprintf('Older observer ODs: L=%.2f, M=%.2f, S=%.2f\n', Lod_old, Mod_old, Sod_old);

%% Example 4: L cone serine and alanine variants
%
% The paper provides a template for the codon 180 serene version of the L cone
% photopigment, and tells you the shift to put in to get the alanine variant.  I believe
% the idea is that these are constructed so that their weighted average hits the
% template for the population average
fprintf('\nExample 4: Other templates\n');
fprintf('==================================================\n');
LserToala_shift = -2.7;
SR_Lserabsorbance = 10.^srLserconelog(wls,Lshift);
SR_Lalaabsorbance = 10.^srLserconelog(wls,LserToala_shift);
fprintf('Peak L serene absorbance at: %.1f nm\n', wls(find(SR_Lserabsorbance == max(SR_Lserabsorbance), 1), 1));
fprintf('Peak L alanine absorbance at: %.1f nm\n', wls(find(SR_Lalaabsorbance == max(SR_Lalaabsorbance), 1), 1));
LserAndAlafigure = figure;
set(gcf,'Position',[100 100 800 1400]);
subplot(2,1,1); hold on;
plot(SR_LMSabsorbances(:,1), log10(SR_LMSabsorbances(:,2)), 'r', 'LineWidth', 2);
plot(SR_LMSabsorbances(:,1), log10(SR_Lserabsorbance), 'y--', 'LineWidth', 2);
plot(SR_LMSabsorbances(:,1), log10(SR_Lalaabsorbance), 'b--', 'LineWidth', 2);
xlim([450 650]);
ylim([-1 0]);
xlabel('Wavelength (nm)');
ylabel('Log10 Sensitivity');
title('L Variants Absorbance');
legend('L', 'Ser', 'Ala', 'Location', 'best');
set(gca,'Color',[0.4 0.4 0.4]);

subplot(2,1,2); hold on;
plot(SR_LMSabsorbances(:,1), SR_LMSabsorbances(:,2), 'r', 'LineWidth', 2);
plot(SR_LMSabsorbances(:,1), SR_Lserabsorbance, 'y--', 'LineWidth', 2);
plot(SR_LMSabsorbances(:,1), SR_Lalaabsorbance, 'b--', 'LineWidth', 2);
xlim([450 650]);
xlabel('Wavelength (nm)');
ylabel('Sensitivity');
title('L Variants Absorbance');
legend('L', 'Ser', 'Ala', 'Location', 'best');
set(gca,'Color',[0.4 0.4 0.4]);