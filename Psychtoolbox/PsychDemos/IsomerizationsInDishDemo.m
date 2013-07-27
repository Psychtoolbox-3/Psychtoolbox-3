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
% NOTE, DHB, 7/19/13. This demo routine and its associated data routines 
% (DefaultPhotoreceptors, FillInPhotoreceptors, PrintPhotoreceptors)
% should be better integrated with the more recent code that
% implements the CIE physiological cone fundamentals, and the
% whole set of stuff should be better documented.  See also
%    IsomerizationsInDishDemo
%    CIEConeFundamentalsTest
%    ComputeCIEConeFundamentals
%    ComputeRawConeFundamentals
%    DefaultPhotoreceptors
%    FillInPhotoreceptors
%    PrintPhotoreceptors
%    RetIrradianceToIsoRecSec
% In particular, there should be some default for the 
% photoreceptors structure that gives one the CIE cone
% fundamentals in all their parametric glory, plus additional
% parameters that yield real energy/quantal sensitivites so
% that the resulting coordinates are isomerization rates in
% real units.  I think that we're close to having that, but
% better documentation and tidying is needed.
%
% 05/06/03 lyin Wrote it.
% 06/26/03 dhb	Rewrote to be self-contained, plus new calling conventions.
% 07/10/03 dhb  Various tuning.
% 07/11/03 dhb  Grab data through subroutines.  Get rid of integration time.
% 04/2/13  dhb  Change clear all to clear, and close figs.
% 04/27/13 dhb  Improve comments.
% 7/19/13  dhb  Print out photoreceptors structure using PrintPhotoreceptors.

%% clear
clear; close all;

%% Set photoreceptor properties.
%
% The photoreceptors structure gets filled with
% key parameters values (pupil size, eye length,
% pre-retinal absorbance, etc.)
%
% The routine DefaultPhotoreceptors is a high level
% call.  It fills in the 'source' fields and some
% values according to high-level descriptor (e.g.,
% ('GuineaPig').  See help for that routine
% for available options.
%
% The routine FillInPhotoreceptors fetches the actual
% values for various fields, depending on the source.
%
% To get a feel for this, check what is in the photoreceptors
% structure after the first call, and then after the second.
photoreceptors = DefaultPhotoreceptors('GuineaPig');
photoreceptors = FillInPhotoreceptors(photoreceptors);

%% Define common wavelength sampling for this script.
% 
% S is [start delta nsamples] for the wavelengths in nm.
% This is standard PTB convention.
S = photoreceptors.nomogram.S;

%% Computation of light spectrum.  This is based on how Lu Yin
% calibrates the irradiance in his apparatus, but the code
% following just needs the incident irradiance in quanta/[sec-um2-wlinterval].
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

%% Do the work in toolbox function.  See help text there
% for details.
[isoPerConeSec,absPerConeSec,photoreceptors] = ...
	RetIrradianceToIsoRecSec(irradianceWatts,S,photoreceptors);

%% Make a plot showing the effective photoreceptor sensitivities in quantal
% units, expressed as probability of absorption.
subplot(1,2,2); hold on
set(plot(SToWls(S),photoreceptors.effectiveAbsorbtance(1,:),'g'),'LineWidth',2);
set(plot(SToWls(S),photoreceptors.effectiveAbsorbtance(2,:),'b'),'LineWidth',2);
set(plot(SToWls(S),photoreceptors.effectiveAbsorbtance(3,:),'k'),'LineWidth',2);
set(title('Effective Absorbtance'),'FontSize',14);
set(xlabel('Wavelength (nm)'),'FontSize',12);
set(ylabel('Probability'),'FontSize',12);
axis([300 800 0 1]);

%% Print out photoreceptors structure
fprintf('\n');
PrintPhotoreceptors(photoreceptors);
fprintf('\n');

%% Print out a table summarizing the calculation.
fprintf('***********************************************\n');
fprintf('Isomerization calculations for retina in a dish\n');
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

