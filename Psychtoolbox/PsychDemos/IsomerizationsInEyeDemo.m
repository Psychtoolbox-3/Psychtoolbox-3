% IsomerizationsInEyeDemo
%
% Shows how to compute photoreceptor isomerizations using toolbox
% routines.  These calculations are for the human eye,
% starting with a spectrum as measured by the PR-650
% in watts/sr-m^2-wlinterval, or with a relative spectrum
% and a photopic troland value.
%
% 07/08/03 dhb  Wrote starting from IsomerizationsInDishDemo.
% 07/11/03 dhb  Grab data through subroutines.  Get rid of integration time.
% 07/15/03 dhb  Take eye size from function.
% 08/14/11 dhb  Comment out saving of T_dogrec at end.  Want to be careful when and where
%               this is done, but the template may be useful someday.
% 03/20/12 dhb  Update cal file for PTB 3.
% 04/09/12 dhb  Add test of irradiance to troland conversion.
% 04/27/13 dhb  More extensive comments.

%% Clear
clear all; close all;

%% Set photoreceptor properties.
%
% The photoreceptors structure gets filled with
% key parameters values (pupil size, eye length,
% pre-retinal absorbance, etc.)
%
% The routine DefaultPhotoreceptors is a high level
% call.  It fills in the 'source' fields and some
% values according to high-level descriptor (e.g.,
% ('LivingHumanFovea').  See help for that routine
% for available options.
%
% The routine FillInPhotoreceptors fetches the actual
% values for various fields, depending on the source.
%
% To get a feel for this, check what is in the photoreceptors
% structure after the first call, and then after the second.
whatCalc = 'LivingHumanFovea';
photoreceptors = DefaultPhotoreceptors(whatCalc);
photoreceptors = FillInPhotoreceptors(photoreceptors);

%% Define common wavelength sampling for this script.
% 
% S is [start delta nsamples] for the wavelengths in nm.
% This is standard PTB convention.
S = photoreceptors.nomogram.S;

%% Get light spectrum.  You can choose various illustrative examples.
%
% See cases of switch statement below for the options available.
whichLight = 'fromMonitorRadiance';
switch (whichLight)
	% Convert from a specification of trolands and relative spectral power distribution.
	% Computation to retinal irradiance is done two ways just for fun.
    % Here we'll start with a xenon arc lamp relative spectrum
	case 'fromTrolands',
		trolands = 100;
        trolandType = 'Photopic';
		load spd_xenonArc;
		irradianceWatts = TrolandsToRetIrradiance(spd_xenonArc,S_xenonArc,trolands, ...
			trolandType,photoreceptors.species,photoreceptors.eyeLengthMM.value);
        irradianceTrolandsCheck = RetIrradianceToTrolands(irradianceWatts,S_xenonArc,trolandType, ...
            photoreceptors.species,photoreceptors.eyeLengthMM.value);
        trolandsCheck = sum(irradianceTrolandsCheck);
        fprintf('\nInput trolands is %0.1f, checked value is %0.1f\n\n',trolands,trolandsCheck);
		irradianceWatts = SplineSpd(S_xenonArc,irradianceWatts,S);
	
		% Another way to do this calculation.  Pupil size should cancel out.  Should get
	    % same answer as above.
		pupilSizeMM = 2;
		pupilAreaMM = pi*(pupilSizeMM/2)^2;
		luminance = TrolandsToLum(trolands,pupilAreaMM);
		radianceWatts = LumToRadiance(spd_xenonArc,S_xenonArc,luminance,trolandType);
		irradianceWattsCheck = RadianceToRetIrradiance(radianceWatts,S_xenonArc,pupilAreaMM,photoreceptors.eyeLengthMM.value);
		figure(1); clf; hold on
		set(plot(SToWls(S),irradianceWatts,'r'),'LineWidth',2);
		set(plot(SToWls(S_xenonArc),irradianceWattsCheck,'k'),'LineWidth',2);
		set(title('Check of trolands to irradiance calculation'),'FontSize',14);
		set(xlabel('Wavelength (mm)'),'FontSize',14);
		set(ylabel('Irradiance'),'FontSize',14);

	% Start with radiance measurements, which we just
	% pull out of the Toolbox's default calibration file.
	case 'fromMonitorRadiance'
		% Load light radiance.  We'll use a monitor white.
		% The original units are watts/sr-m^2-wlinterval.
		cal = LoadCalFile('PTB3TestCal');
		radianceWatts = SplineSpd(cal.S_device,sum(cal.P_device,2),S);
		
		% Find pupil area, needed to get retinal irradiance.  We compute
		% pupil area based on the luminance of stimulus according to the
        % algorithm specified in the photoreceptors structure.
		load T_xyz1931
		T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
		theXYZ = T_xyz*radianceWatts; theLuminance = theXYZ(2);
		[nil,pupilAreaMM] = PupilDiameterFromLum(theLuminance,photoreceptors.pupilDiameter.source);
		
		% Convert radiance of source to retinal irradiance and convert to quantal units.
		irradianceWatts = RadianceToRetIrradiance(radianceWatts,S, ...
			pupilAreaMM,photoreceptors.eyeLengthMM.value);

	% This light as well as some parameter tweaking are here to match a parameterization the Brian Wandell supplied
	% to match what his code to do these computations produces.  Note also
	% the mucking with the photoreceptors structure.  Wandell estimates
	% L, M, S isomerizations/cone-sec of 16.5, 12.68, 2.27.
	case 'fromUniformQuantalSpd',
		% Load corneal cone sensitivities in energy units, convert to quantal sensitivities
		% and set specified peak absorbtance.
		%
		% Note that overwriting the isomerizationAbsorbtance in the photoreceptors structure
		% makes the isomerization computation work, but not the absorbtions calculation, which
		% will be done with what was produced by FillInPhotoreceptors called above.  This is
		% not a recommended compute path for the toolbox code, but is done here to match Wandell's
		% parameterization.
		load T_cones_ss2; T_cones = T_cones_ss2; S_cones = S_cones_ss2;
		% load T_cones_ss10; T_cones = T_cones_ss10; S_cones = S_cones_ss10;
		% load T_cones_smj; T_cones = T_cones_smj; S_cones = S_cones_smj;
		% load T_cones_sp; T_cones = T_cones_sp; S_cones = S_cones_sp;
		peakIsomerizationEfficiency = [0.27 0.23 0.07]';
		T_cones = SplineCmf(S_cones,QuantaToEnergy(S_cones,T_cones')',S);
		T_cones(1,:) = T_cones(1,:)/max(T_cones(1,:));
		T_cones(2,:) = T_cones(2,:)/max(T_cones(2,:));
		T_cones(3,:) = T_cones(3,:)/max(T_cones(3,:));
		T_cones = diag(peakIsomerizationEfficiency)*T_cones;
		photoreceptors.isomerizationAbsorbtance = T_cones;

		% Get spectral luminous efficiency function
		load T_xyz1931; T_xyz = T_xyz1931; S_xyz = S_xyz1931;
		% load T_xyzJuddVos; T_xyz = T_xyzJuddVos; S_xyz = S_xyzJuddVos;
		T_Y = 683*SplineCmf(S_xyz,T_xyz(2,:),S);

		% Create a spectrally uniform spd (in quantal units), and convert
		% to energy units.
		uniformSpd = QuantaToEnergy(S,ones(S(3),1));

		% Normalize to radiance corresponding to 1 cd/m2.
		normConst = T_Y*uniformSpd;
		radianceWatts = uniformSpd/normConst;

		% Set pupil diameter for 1mm2 pupil area, photoreceptor diameter for 4mm2 collecting
		% area.  Set eye length to 17 mm.
		photoreceptors.pupilDiameter.value = 2*sqrt(1/pi);
		pupilAreaMM = pi*(photoreceptors.pupilDiameter.value/2)^2;
		photoreceptors.ISdiameter.value = [2*sqrt(4/pi) 2*sqrt(4/pi) 2*sqrt(4/pi)]';
		photoreceptors.eyeLengthMM.value = 17;
		irradianceWatts = RadianceToRetIrradiance(radianceWatts,S,pupilAreaMM,photoreceptors.eyeLengthMM.value );
end

%% Get retinal irradiance in quanta/[sec-um^2-wlinterval]
irradianceQuanta = EnergyToQuanta(S,irradianceWatts);
figure(2); clf; set(gcf,'Position',[100 400 700 300]);
subplot(1,2,1); hold on
set(plot(SToWls(S),irradianceQuanta,'r'),'LineWidth',2);
set(title('Light Spectrum'),'FontSize',14);
set(xlabel('Wavelength (nm)'),'FontSize',12);
set(ylabel('Quanta/sec-um^2-wlinterval'),'FontSize',12);

%% Do the work in toolbox function
[isoPerConeSec,absPerConeSec,photoreceptors] = ...
	RetIrradianceToIsoRecSec(irradianceWatts,S,photoreceptors);

% Make a plot showing the effective photoreceptor sensitivities in quantal
% units, expressed as probability of isomerization.
subplot(1,2,2); hold on
set(plot(SToWls(S),photoreceptors.isomerizationAbsorbtance(1,:),'r'),'LineWidth',2);
set(plot(SToWls(S),photoreceptors.isomerizationAbsorbtance(2,:),'g'),'LineWidth',2);
set(plot(SToWls(S),photoreceptors.isomerizationAbsorbtance(3,:),'b'),'LineWidth',2);
set(title('Isomerization Absorbtance'),'FontSize',14);
set(xlabel('Wavelength (nm)'),'FontSize',12);
set(ylabel('Probability'),'FontSize',12);
axis([300 800 0 1]);

% Print out a table summarizing the calculation.
fprintf('***********************************************\n');
fprintf('Isomerization calculations for living human retina\n');
fprintf('\n');
fprintf('Calculations done using:\n');
fprintf('\t%s estimates for photoreceptor IS diameter\n',photoreceptors.ISdiameter.source);
fprintf('\t%s estimates for photoreceptor OS length\n',photoreceptors.OSlength.source);
fprintf('\t%s estimates for receptor specific density\n',photoreceptors.specificDensity.source);
fprintf('\t%s photopigment nomogram\n',photoreceptors.nomogram.source);
fprintf('\t%s estimates for lens density\n',photoreceptors.lensDensity.source);
fprintf('\t%s estimates for macular pigment density\n',photoreceptors.macularPigmentDensity.source);
fprintf('\t%s method for pupil diameter calculation\n',photoreceptors.pupilDiameter.source);
fprintf('\t%s estimate (%g mm) for axial length of eye\n',photoreceptors.eyeLengthMM.source,photoreceptors.eyeLengthMM.value);
fprintf('\n');
fprintf('Photoreceptor Type             |\t       L\t       M\t     S\n');
fprintf('______________________________________________________________________________________\n');
fprintf('\n');
fprintf('Lambda max                     |\t%8.1f\t%8.1f\t%8.1f\t nm\n',photoreceptors.nomogram.lambdaMax);
fprintf('Outer Segment Length           |\t%8.1f\t%8.1f\t%8.1f\t um\n',photoreceptors.OSlength.value);
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

% Allow dumping out of photoreceptor sensitivities into a file for use elsewhere.  We want energy sensitivities
% for this purpose
% switch (whatCalc)
%     % Dog receptors (L, S, rod) in energy units, normalized to max of 1.
%     case 'LivingDog'
%         T_dogrec = EnergyToQuanta(S,photoreceptors.isomerizationAbsorbtance')';
%         for i = 1:3
%             T_dogrec(i,:) = T_dogrec(i,:)/max(T_dogrec(i,:));
%         end
%         S_dogrec = S;
%         save T_dogrec T_dogrec S_dogrec
%     otherwise
% end

