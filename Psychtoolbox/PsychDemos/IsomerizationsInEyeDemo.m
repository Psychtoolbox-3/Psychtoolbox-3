% IsomerizationsInEyeDemo
%
% Shows how to compute photoreceptor isomerizations using toolbox
% routines.  These calculations are for the human eye,
% starting with a spectrum as measured by the PR-650
% in watts/sr-m^2-wlinterval, or with a relative spectrum
% and a photopic troland value.
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
% 07/08/03 dhb  Wrote starting from IsomerizationsInDishDemo.
% 07/11/03 dhb  Grab data through subroutines.  Get rid of integration time.
% 07/15/03 dhb  Take eye size from function.
% 08/14/11 dhb  Comment out saving of T_dogrec at end.  Want to be careful when and where
%               this is done, but the template may be useful someday.
% 03/20/12 dhb  Update cal file for PTB 3.
% 04/09/12 dhb  Add test of irradiance to troland conversion.
% 04/27/13 dhb  More extensive comments.
% 7/19/13  dhb  Print out photoreceptors structure using PrintPhotoreceptors.
%          dhb  Add monochromatic light option to the section that starts with trolands.

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
photoreceptors.eyeLengthMM.source = 'LeGrand';
photoreceptors = FillInPhotoreceptors(photoreceptors);

%% Define common wavelength sampling for this script.
% 
% S is [start delta nsamples] for the wavelengths in nm.
% This is standard PTB convention.
S = photoreceptors.nomogram.S;

%% XYZ color matching functions
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
T_Y = T_xyz(2,:);
        
%% Get light spectrum.  You can choose various illustrative examples.
%
% Available options:
%  'fromTrolands'
%  'fromMonitorRadiance'
%  'fromUniformQuantalSpd'
whichInputType = 'fromTrolands';
switch (whichInputType)
    
    % Start with troland value and a relative spectrum
    case 'fromTrolands'
        fprintf('Computing from troland value and relative spectrum\n');
        % Give troland value and type
        %
        % Type may be 'Photopic', 'Scotopic', or 'JuddVos'
        trolands = 1;
        trolandType = 'Photopic';
        switch (trolandType)
            case 'Photopic'
                fprintf('Using photopic trolands\n');
            case 'Scotopic'
                fprintf('Using scotopic trolands\n');
            case 'JuddVos'
                fprintf('Using Judd-Vos luminosity function in troland calculations\n');
            otherwise
                fprintf('Unknown troland type specified');
        end
        
        %% Pupil.
        %
        % We do these computations for a fixed pupil size, ignoring the
        % pupilDiamter field of the photoreceptors structure.  That field
        % is set up to specify a source formula that estimates pupil diameter
        % from luminance.
        %
        % Since we are starting in trolands, the pupil size shouldn't actually
        % effect the calculations, except for finding the radiance that is
        % equivalent to the specified troland value.  
        % 
        % We remove the pupilDiameter.source field to make sure we aren't sending
        % mixed messages about how we want to handle pupil diameter.
        photoreceptors.pupilDiameter.value = 2;
        photoreceptors.pupilDiameter = rmfield(photoreceptors.pupilDiameter,'source');
        pupilAreaMm2 = pi*(photoreceptors.pupilDiameter.value/2)^2;
        
        % Specify relative spectrum to be used in
        % conversion to a full spectrum.
        %
        % Choices are:
        %  'Monochromatic'
        %  'XenonArc'
        %
        % If type is 'Monochromatic', must specify
        % wavelengthNm.
        spectrumType = 'Monochromatic';
        switch (spectrumType)
            case 'Monochromatic'
                monoWavelengthNm = 550;  
                wls = SToWls(S);
                monoWavelengthIndex = find(wls == monoWavelengthNm);
                if (isempty(monoWavelengthIndex))
                    error('No sample wavelength matches desired wavelength');
                end
                spd_fromTrolands = zeros(size(wls));
                spd_fromTrolands(monoWavelengthIndex) = 1;
                fprintf('Using monochromatic %d nm light as relative spectrum\n',monoWavelengthNm);
            case 'XenonArc'
                load spd_xenonArc;
                spd_fromTrolands = SplineSpd(S_xenonArc,spd_xenonArc,S);
                clear S_xenonArc spd_xenonArc
                fprintf('Using the spectrum of a xenon arc lamp as relative spectrum\n');
            otherwise
                error('Unknown spectrum type specified');
        end
        
        % Convert trolands back to spectral retinal irradiance.  This
        % depends on the pupil size and eye length specified.
		irradianceWattsPerUm2 = TrolandsToRetIrradiance(spd_fromTrolands,S,trolands, ...
			trolandType,photoreceptors.species,photoreceptors.eyeLengthMM.value);
        irradianceTrolandsCheck = RetIrradianceToTrolands(irradianceWattsPerUm2,S,trolandType, ...
            photoreceptors.species,photoreceptors.eyeLengthMM.value);
        trolandsCheck = sum(irradianceTrolandsCheck);
        fprintf('Input troland value is %0.1f, checked value is %0.1f\n',trolands,trolandsCheck);
		irradianceWattsPerUm2 = SplineSpd(S,irradianceWattsPerUm2,S);
	
		% Another way to do this calculation.  Pupil size should cancel out.  Should get
	    % same answer as above.  This has as a byproduct computing a stimulus radiance,
        % which is useful for some of the common printout below.
        luminanceCdM2FromTrolands = TrolandsToLum(trolands,pupilAreaMm2);
		radianceWattsPerM2Sr = LumToRadiance(spd_fromTrolands,S,luminanceCdM2FromTrolands,trolandType); 
        photopicLuminanceCdM2 = T_Y*radianceWattsPerM2Sr;
		irradianceWattsCheck = RadianceToRetIrradiance(radianceWattsPerM2Sr,S,pupilAreaMm2,photoreceptors.eyeLengthMM.value);
		figure(1); clf; hold on
		set(plot(SToWls(S),irradianceWattsPerUm2,'r'),'LineWidth',2);
		set(plot(SToWls(S),irradianceWattsCheck,'k'),'LineWidth',2);
		set(title('Check of trolands to irradiance calculation'),'FontSize',14);
		set(xlabel('Wavelength (mm)'),'FontSize',14);
		set(ylabel('Irradiance'),'FontSize',14);
        
        % For case if monochromatic light, can check retinal irradiance against
        % the direct formulae provided in W&S, 2cd edition, p. 105, eqs 2.4.4 
        % Note that these equations are missing a factor of
        % the wavelength in the numerator for the quantal conversions, as
        % pointed out by Makous in his 1997 JOSA paper.
        %
        % This only works for 'Photopic' and 'Scotopic' trolands.
        %
        % There must be an eye length implicit in this calculation. Using
        % the LeGrand model eye length gives good agreement between these
        % and our more general calculations done in the main section below.
        %
        % Makous (1997), JOSA A, 14, p. 2331 gives retinal illuminances of
        % 5.44 quanta/[um2-sec] for 1 scotopic troland, and 14.65 quanta/[um2-sec]
        % for 1 photopic troland, with the calulations specified for 510 nm.
        % The calculations here, done in several different ways, yield
        % 5.442 (scotopic, agrees) and 26.85 (photopic, does not agree).  But at 
        % 550 nm, the code here yields 14.64 quanta/[um2-sec], which seems close
        % enough to the provided 14.65 to make me think that Makous' value is
        % actually for a wavelength close to 550 nm.  That would be a more typical
        % wavelength at which to do photopic calculations, despite what the text
        % in the paper says.
        %
        % Note that there are also errors in Tables 3 and 4 of the Makous
        % paper, and that corrected values appear in Tables 3 and 2 of
        % Makous (2004), Scotopic vision, In The Visual Neurosciences,
        % Werner and Chalupa (eds).  What does not seem to be specified
        % in either place is the wavelengths used in the calculations of
        % the two tables.
        if (strcmp(spectrumType,'Monochromatic'))
            switch (trolandType)
                case 'Photopic'
                    load T_xyz1931;
                    T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
                    clear T_xyz1931 S_xyz1931
                    magicFactorE = 5.261e-12;
                    magicFactorQ = 2.649e13;                  
                case 'Scotopic'
                    load T_rods;
                    T_vLambda = SplineCmf(S_rods,T_rods,S);
                    magicFactorE = 2.114e-12;
                    magicFactorQ = 1.064e13;
            end
        end
        if (strcmp(spectrumType,'Monochromatic'))
            switch (trolandType)
                case {'Photopic','Scotopic'}
                    irradianceDirectWattsPerMm2Check = magicFactorE*trolands/T_vLambda(monoWavelengthIndex);
                    irradianceDirectQuantaPerMm2SecCheck = 1e-9*monoWavelengthNm*magicFactorQ*trolands/T_vLambda(monoWavelengthIndex);
                    fprintf('Retinal irradiance in area units computed from from trolands via Wyszecki and Stiles formulae\n');
                    fprintf('\t%0.4g Watts/mm2\n\t%0.4g quanta/[mm2-sec]\n',sum(irradianceDirectWattsPerMm2Check),sum(irradianceDirectQuantaPerMm2SecCheck));
            end
        end
        
	% Start with radiance measurements, which we just
	% pull out of the Toolbox's default calibration file.
	case 'fromMonitorRadiance'
		% Load light radiance.  We'll use a monitor white.
		% The original units are watts/sr-m^2-wlinterval.
		cal = LoadCalFile('PTB3TestCal');
		radianceWattsPerM2Sr = SplineSpd(cal.S_device,sum(cal.P_device,2),S);
		
		% Find pupil area, needed to get retinal irradiance.  We compute
		% pupil area based on the luminance of stimulus according to the
        % algorithm specified in the photoreceptors structure.
		theXYZ = T_xyz*radianceWattsPerM2Sr; theLuminance = theXYZ(2);
		[nil,pupilAreaMm2] = PupilDiameterFromLum(theLuminance,photoreceptors.pupilDiameter.source);
        photopicLuminanceCdM2 = T_Y*radianceWattsPerM2Sr;
		
		% Convert radiance of source to retinal irradiance and convert to quantal units.
		irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S, ...
			pupilAreaMm2,photoreceptors.eyeLengthMM.value);

	% This light as well as some parameter tweaking are here to match a parameterization the Brian Wandell supplied
	% to match what his code to do these computations produces.  Note also
	% the mucking with the photoreceptors structure.  Wandell estimates
	% L, M, S isomerizations/cone-sec of 16.5, 12.68, 2.27.  These are very close to the numbers
    % we get here.
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

		% Create a spectrally uniform spd (in quantal units), and convert
		% to energy units.
		uniformSpd = QuantaToEnergy(S,ones(S(3),1));

		% Normalize to radiance corresponding to 1 cd/m2.
		normConst = T_Y*uniformSpd;
		radianceWattsPerM2Sr = uniformSpd/normConst;
        photopicLuminanceCdM2 = T_Y*radianceWattsPerM2Sr;

		% Set pupil diameter for 1mm2 pupil area, photoreceptor diameter for 4mm2 collecting
		% area.  Set eye length to 17 mm.
		photoreceptors.pupilDiameter.value = 2*sqrt(1/pi);
		pupilAreaMm2 = pi*(photoreceptors.pupilDiameter.value/2)^2;
		photoreceptors.ISdiameter.value = [2*sqrt(4/pi) 2*sqrt(4/pi) 2*sqrt(4/pi)]';
		photoreceptors.eyeLengthMM.value = 17;
		irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S,pupilAreaMm2,photoreceptors.eyeLengthMM.value );
end

%% Print out a whole bunch of quantities that are equivalent to the radiance, given
% other eye parameters.
radianceWattsPerCm2Sr = (10.^-4)*radianceWattsPerM2Sr;
radianceQuantaPerCm2SrSec = EnergyToQuanta(S,radianceWattsPerCm2Sr);
degPerMm = RetinalMMToDegrees(1,photoreceptors.eyeLengthMM.value);
irradianceWattsPerUm2Check = RadianceToRetIrradiance(radianceWattsPerM2Sr,S,pupilAreaMm2,photoreceptors.eyeLengthMM.value);
if (any(abs(irradianceWattsPerUm2 - irradianceWattsPerUm2Check) > 1e-10))
    error('Back computation of retinal irradiance from radiance does not check');
end
irradianceScotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Scotopic', [], num2str(photoreceptors.eyeLengthMM.value));
irradiancePhotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Photopic', [], num2str(photoreceptors.eyeLengthMM.value));
irradianceQuantaPerUm2Sec = EnergyToQuanta(S,irradianceWattsPerUm2);
irradianceWattsPerCm2 = (10.^8)*irradianceWattsPerUm2;
irradianceQuantaPerCm2Sec = (10.^8)*irradianceQuantaPerUm2Sec;
irradianceQuantaPerMm2Sec = (10.^-2)*irradianceQuantaPerCm2Sec;
irradianceQuantaPerUm2Sec = (10.^-6)*irradianceQuantaPerMm2Sec;
irradianceQuantaPerDeg2Sec = (degPerMm^2)*irradianceQuantaPerMm2Sec;

% Print out photoreceptor stucture information
fprintf('\n');
PrintPhotoreceptors(photoreceptors);
fprintf('\n');

% Radiometric iformation
fprintf('Luminance %0.3f cd/m2\n',photopicLuminanceCdM2);
fprintf('Stimulus retinal irradiance %0.4g (%0.1f log10) watts/cm2\n',sum(irradianceWattsPerCm2),log10(sum(irradianceWattsPerCm2)));
fprintf('Stimulus retinal irradiance %0.4g (%0.1f log10) watts/mm2\n',1e-2*sum(irradianceWattsPerCm2),log10(sum(irradianceWattsPerCm2)));
fprintf('Stimulus retinal irradiance %0.4g (%0.1f log10) quanta/[cm2-sec]\n',sum(irradianceQuantaPerCm2Sec),log10(sum(irradianceQuantaPerCm2Sec)));
fprintf('Stimulus retinal irradiance %0.4g (%0.1f log10) quanta/[mm2-sec]\n',sum(irradianceQuantaPerMm2Sec),log10(sum(irradianceQuantaPerMm2Sec)));
fprintf('Stimulus retinal irradiance %0.4g (%0.1f log10) quanta/[um2-sec]\n',sum(irradianceQuantaPerUm2Sec),log10(sum(irradianceQuantaPerUm2Sec)));
fprintf('Stimulus retinal irradiance %0.4g (%0.1f log10) quanta/[deg2-sec]\n',sum(irradianceQuantaPerDeg2Sec),log10(sum(irradianceQuantaPerDeg2Sec)));
fprintf('\n');
        
%% Get retinal irradiance in quanta/[sec-um^2-wlinterval]
irradianceQuanta = EnergyToQuanta(S,irradianceWattsPerUm2);
figure(2); clf; set(gcf,'Position',[100 400 700 300]);
subplot(1,2,1); hold on
set(plot(SToWls(S),irradianceQuanta,'r'),'LineWidth',2);
set(title('Light Spectrum'),'FontSize',14);
set(xlabel('Wavelength (nm)'),'FontSize',12);
set(ylabel('Quanta/sec-um^2-wlinterval'),'FontSize',12);

%% Do the work in toolbox function
[isoPerConeSec,absPerConeSec,photoreceptors] = ...
	RetIrradianceToIsoRecSec(irradianceWattsPerUm2,S,photoreceptors);

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

