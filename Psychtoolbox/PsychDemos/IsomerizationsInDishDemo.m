% IsomerizationsInDishDemo
%
% Shows how to compute photoreceptor isomerizations using toolbox
% routines.  These calculations are for a retinal preparation
% that consists of a retina in a dish and thus do not include
% any preretinal absorbtion.  Parameters are set up for a
% Guinea Pig retina.  The file spd_apparatusrel.mat contains
% the relative spectrum of the light with respect to which
% the computation is accomplished.
%
% See also: IsomerizationInEyeDemo.
%
% 05/06/03 lyin Wrote it.
% 06/26/03 dhb	Rewrote to be self-contained, plus new calling conventions.
% 07/10/03 dhb  Various tuning.
% 07/11/03 dhb  Grab data through subroutines.  Get rid of integration time.

% clear
clear all;

% Set some photoreceptor properties.
photoreceptors = DefaultPhotoreceptors('GuineaPig');
photoreceptors = FillInPhotoreceptors(photoreceptors);

% Define common wavelength sampling for this script.
S = photoreceptors.nomogram.S;

% Computation of light spectrum.  This is based on how Lu Yin
% calibrates the irradiance in his apparatus.
%
% First the relative spectrum of the light is measured using the CVI
% spectral radiometer.  Then the absolute power of the light is measured
% using a broadband radiometer that measures total power in watts.  The
% radiometer receives the light through a pinhole of known size.  With this
% information, we can get incident power in watts/um^2-wlinterval.
% Routine WattsToRetIrradiance does the conversion, taking into account
% instrument calibration information as well as the relative spectrum,
% which here we just load.
load spd_apparatusrel										% Relative spectrum measured in apparatus
readingInWatts = 0.0128e-6 ;  								% Watts (0.0128 uWatts)
[irradianceWatts,irradianceS] = WattsToRetIrradiance( ...
	spd_apparatusrel,S_apparatusrel,readingInWatts);
irradianceQuanta = EnergyToQuanta(irradianceS,irradianceWatts);
figure(1); clf; set(gcf,'Position',[100 400 700 300]);
subplot(1,2,1); hold on
set(plot(SToWls(irradianceS),irradianceQuanta,'r'),'LineWidth',2);
set(title('Light Spectrum'),'FontSize',14);
set(xlabel('Wavelength (nm)'),'FontSize',12);
set(ylabel('Quanta/sec-um^2-wlinterval'),'FontSize',12);

% Do the work in toolbox function
[isoPerConeSec,absPerConeSec,photoreceptors] = ...
	RetIrradianceToIsoRecSec(irradianceWatts,S,photoreceptors);

% Make a plot showing the effective photoreceptor sensitivities in quantal
% units, expressed as probability of absorption.
subplot(1,2,2); hold on
set(plot(SToWls(S),photoreceptors.effectiveAbsorbtance(1,:),'g'),'LineWidth',2);
set(plot(SToWls(S),photoreceptors.effectiveAbsorbtance(2,:),'b'),'LineWidth',2);
set(plot(SToWls(S),photoreceptors.effectiveAbsorbtance(3,:),'k'),'LineWidth',2);
set(title('Effective Absorbtance'),'FontSize',14);
set(xlabel('Wavelength (nm)'),'FontSize',12);
set(ylabel('Probability'),'FontSize',12);
axis([300 800 0 1]);

% Print out a table summarizing the calculation.
fprintf('***********************************************\n');
fprintf('Isomerization calculations for retina in a dish\n');
fprintf('\n');
fprintf('Calculations done using:\n');
fprintf('\t%s estimates for photoreceptor OS diameter\n',photoreceptors.OSdiameter.source);
fprintf('\t%s estimates for photoreceptor IS diameter\n',photoreceptors.ISdiameter.source);
fprintf('\t%s estimates for photoreceptor OS length\n',photoreceptors.OSlength.source);
fprintf('\t%s estimates for receptor specific density\n',photoreceptors.specificDensity.source);
fprintf('\t%s photopigment nomogram\n',photoreceptors.nomogram.source);
fprintf('\n');
fprintf('Photoreceptor Type             |\t       M\t       S\t     Rod\n');
fprintf('______________________________________________________________________________________\n');
fprintf('\n');
fprintf('Lambda max                     |\t%8.1f\t%8.1f\t%8.1f\t nm\n',photoreceptors.nomogram.lambdaMax);
fprintf('Outer Segment Length           |\t%8.1f\t%8.1f\t%8.1f\t um\n',photoreceptors.OSlength.value);
fprintf('Outer Segment Diameter         |\t%8.1f\t%8.1f\t%8.1f\t um\n',photoreceptors.OSdiameter.value);
fprintf('Inner Segment Diameter         |\t%8.1f\t%8.1f\t%8.1f\t um\n',photoreceptors.ISdiameter.value);
fprintf('\n');
fprintf('Axial Specific Density         |\t%8.3f\t%8.3f\t%8.3f\t /um\n',photoreceptors.specificDensity.value);
fprintf('Axial Optical Density          |\t%8.3f\t%8.3f\t%8.3f\n',photoreceptors.axialDensity.value);
fprintf('Peak isomerization prob.       |\t%8.3f\t%8.3f\t%8.3f\n',max(photoreceptors.isomerizationAbsorbtance,[],2));
fprintf('______________________________________________________________________________________\n');
fprintf('\n');
fprintf('Absorption Rate                |\t%4.2e\t%4.2e\t%4.2e\t quanta/photoreceptor-sec\n',...
	absPerConeSec);
fprintf('Isomerization Efficiency       |\t%8.3f\t%8.3f\t%8.3f\n',...
	photoreceptors.quantalEfficiency.value);
fprintf('Isomerization Rate             |\t%4.2e\t%4.2e\t%4.2e\t iso/photoreceptor-sec\n',...
	 isoPerConeSec);
fprintf('In log10 units                 |\t%8.2f\t%8.2f\t%8.2f\t log10(iso)/photoreceptor-sec\n',...
	 log10(isoPerConeSec));
fprintf('______________________________________________________________________________________\n');

