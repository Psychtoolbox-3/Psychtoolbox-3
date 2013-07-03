% AnsiZ136DeloriTest
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Tests radiometric calculations and light safety calculations.
%
% Reads in tab delimted text file AnsiZ136MPEDeloriTestInput.txt that has a set of 
% rows in it.  Each row provides stimulus properties for a monochromtic
% light and values obtained from Delori's spreadsheet.
% (There are also some other test values entered by hand in the text file, for
% other radiometric quanties that I could not find computed in Delori's spreadsheet.)
%
% Delori's spreadsheet is described in footnote 49 of Delori et al. (2007,
% JOSA A, 24, pp.1250-1265).  That footnote indicates that the spreadsheet
% will be shared as long as recipients accept full responsibility for its use.
% The version of Delori's spreadsheet used to get these numbers is 0cLC_0.9.9.

% As of 5/1/13 (Notes from DHB):
%   I am unable to get complete agreement between Delori's calculations and
%   the PTB calculations.  But I get agreement for some wavelength, stimulus
%   size, and stimulus duration choices.  To get agreement, this routine
%   uses the effective pupil size formula from the Delori et al. 2007 paper
%   when converting the limit to power in the pupil.
%
%   A key question in going from the number in the 2007 standard and retinal
%   illuminant is what pupil size was assumed in the standard.  Delori et al.
%   (2007) give a formula for this, which (I think) goes from the radiant
%   power in light at the cornea overfilling the pupil to the pupil size. 
%   See Eq. 8.  That is implemented in this test program.
%
%   To get the degree of agreement I obtained, I also had to turn off
%   the limiting cone aperture calculation that is described in Table 2 of the
%   2007 standard.  This is controlled by a flag in the input file.
%
%   Conversions to retinal illuminance done wtth the pupil diameter and
%   eye length specified in the input file.  Delori et al. (2007) indicate
%   that the standard assumes an eye length of 17 mm.  The calculation uses
%   this to compute the allowable retinal illuminance, but if you think some
%   other length is better, you can fill it in.
%
%   In the cases where I don't get agreement with the spreadsheet, the discrepancy
%   arises in the value used for constant T2.  Delori's spreadsheet seems to use
%   10000 for a 2 degree stimulus, but as I read the standard (Table 6, p. 76), this
%   value has a maximum of 100 and a smaller value (21.83) for T2.  This matters when
%   the thermal limit is controlling, for some stimulus durations.  Currently this
%   happens for test conditions 3 and 6.  The code here
%   reproduces what I get when I do the computations by hand.  See files in this
%   directory
%     MPEByHand_580_2deg_100sec.doc
%     MPEByHand_700_2deg_20000sec.doc
%   I have emailed Delori (on 5/2/13) to ask him about this.  And also about the
%   limiting cone aperature calculation mentioned above.
%
% 2/27/13  dhb  Wrote it.

%% Clear and close
clear; close all

%% Load conditions
conditionStructs = ReadStructsFromText('AnsiZ136MPEDeloriTestInput.txt');

%% Loop over conditions and report
for i = 1:length(conditionStructs)
    % Set variables for this condition
    retinalIlluminanceUWattsPerCm2 = conditionStructs(i).retinalIlluminanceUWattsPerCm2;
    wavelengthNm = conditionStructs(i).wavelengthNm;
    stimulusDiameterDegrees = conditionStructs(i).stimulusDiamDegrees;
    stimulusDurationSeconds = conditionStructs(i).stimulusDurationSecs;
    eyeLengthMm = conditionStructs(i).eyeLengthMm;
    pupilDiameterMm = conditionStructs(i).pupilDiameterMm;
    ansiEyeLengthMm = conditionStructs(i).ansiEyeLengthMm;
    ansiPupilDiameterMm = conditionStructs(i).ansiPupilDiameterMm; 
    CONELIMITFLAG = conditionStructs(i).CONELIMITFLAG;
    fprintf('**********\nCondition %d\n\tInput retinal illuminance of %0.1f uWatts/cm2\n',i,retinalIlluminanceUWattsPerCm2);
    fprintf('\t\tWavelength %d nm\n',wavelengthNm);
    fprintf('\t\tStimulus diamter %0.1f degrees\n',stimulusDiameterDegrees);
    fprintf('\t\tStimulus duration %0.1f seconds\n',stimulusDurationSeconds);
    fprintf('\t\tEye length %0.1f mm\n',eyeLengthMm);
    fprintf('\t\tPupil diameter %0.1f mm\n',pupilDiameterMm);
    fprintf('\t\tAssuming ANSI standard eye length of %0.1f mm\n',ansiEyeLengthMm);
    fprintf('\t\tAssuming ANSI pupil diameter of %0.1f mm\n',ansiPupilDiameterMm);
    if (CONELIMITFLAG)
        fprintf('\t\tIncluding limiting cone angle calculation\n');
    else
        fprintf('\t\tExcluding limiting cone angle calculation\n');
    end
    
    % Get comparison values from Delori spreadsheet.  These need to be computed
    % by hand using the spreadsheet and then entered into the condition file.
    % We enter -1 in the spreadsheet for values not computed.
    deloriRadianceMWattsPerCm2Sr = conditionStructs(i).deloriRadianceMWattsPerCm2Sr;
    deloriCornealIrradianceUWattsPerCm2 = conditionStructs(i).deloriCornealIrradianceUWattsPerCm2;
    deloriPowerInPupilMWatts = conditionStructs(i).deloriPowerInPupilMWatts;
    deloriLog10PhotopicTrolands = conditionStructs(i).deloriLog10PhotopicTrolands;
    deloriLog10ScotopicTrolands = conditionStructs(i).deloriLog10ScotopicTrolands;
    deloriMPECb = conditionStructs(i).deloriMPECb;
    deloriMPECe = conditionStructs(i).deloriMPECe;
    deloriMPECt = conditionStructs(i).deloriMPECt;
    deloriMPEPupilFactor = conditionStructs(i).deloriMPEPupilFactor;
    deloriEffectivePupilDiameterMm = conditionStructs(i).deloriEffectivePupilDiameterMm;
    deloriMPETgamma = conditionStructs(i).deloriMPETgamma;
    deloriMPET2 = conditionStructs(i).deloriMPET2;
    deloriMPERetinalIrradianceWattsPerCm2 = conditionStructs(i).deloriMPERetinalIrradianceWattsPerCm2;
    deloriMPEPowerInPupilMWatts = conditionStructs(i).deloriMPEPowerInPupilMWatts;
    
    % In some cases we have other comparison values for other sources.
    % We enter -1 in the spreadsheet for values not computed.
    checkPhotopicLuminanceCdM2 = conditionStructs(i).checkPhotopicLuminanceCdM2;
    checkRetinalIlluminanceQuantaPerCm2Sec = conditionStructs(i).checkRetinalIlluminanceQuantaPerCm2Sec;
    
    % Do unit conversions and print with comparisons when such are available
    retinalIlluminanceWattsPerCm2 = 1e-6*retinalIlluminanceUWattsPerCm2;
    retinalIlluminanceWattsPerUm2 = 1e-8*retinalIlluminanceWattsPerCm2;
    retinalIlluminanceQuantaPerCm2Sec = EnergyToQuanta(wavelengthNm,retinalIlluminanceWattsPerCm2);
    photopicTrolands = RetIrradianceToTrolands(retinalIlluminanceWattsPerUm2, wavelengthNm, 'Photopic','Human',num2str(eyeLengthMm));
    scotopicTrolands = RetIrradianceToTrolands(retinalIlluminanceWattsPerUm2, wavelengthNm, 'Scotopic','Human',num2str(eyeLengthMm));
    photopicLuminanceCdM2 = TrolandsToLum(photopicTrolands,(pi/4)*pupilDiameterMm^2);
    radianceWattsPerM2Sr = RetIrradianceToRadiance(retinalIlluminanceWattsPerUm2,wavelengthNm,(pi/4)*pupilDiameterMm^2,eyeLengthMm);
    radianceMWattsPerCm2Sr = 1e3*1e-4*radianceWattsPerM2Sr;
    cornealIrradianceMWattsPerCm2 = RadianceAndDegrees2ToCornIrradiance(radianceMWattsPerCm2Sr,(pi/4)*stimulusDiameterDegrees^2);
    cornealIrradianceUWattsPerCm2 = 1e3*cornealIrradianceMWattsPerCm2;
    powerInPupilUWatts = cornealIrradianceUWattsPerCm2*(1e-2)*(pi/4)*pupilDiameterMm^2;
    powerInPupilMWatts = 1e-3*powerInPupilUWatts;
    pupilEnergyMJoules = powerInPupilMWatts*stimulusDurationSeconds;
    
    % Print out results, and comparisons if we have them
    AnsiZ136MPEPrintConditionalComparison('\tConverts to retinal illuminance of %0.1f log10 quanta/[cm2-sec]','%0.1f',retinalIlluminanceQuantaPerCm2Sec,checkRetinalIlluminanceQuantaPerCm2Sec,true);
    AnsiZ136MPEPrintConditionalComparison('\tConverts to %0.2f log10 photopic trolands','%0.2f',photopicTrolands,10^deloriLog10PhotopicTrolands,true);
    AnsiZ136MPEPrintConditionalComparison('\tConverts to %0.2f log10 scotopic trolands','%0.2f',scotopicTrolands,10^deloriLog10ScotopicTrolands,true);
    AnsiZ136MPEPrintConditionalComparison('\tConverts to %0.1f cd/m2','%0.1f',photopicLuminanceCdM2,checkPhotopicLuminanceCdM2,false);
    AnsiZ136MPEPrintConditionalComparison('\tConverts to radiance %0.1f mWatts/[cm2-sr]','%0.1f',radianceMWattsPerCm2Sr,deloriRadianceMWattsPerCm2Sr,false);
    AnsiZ136MPEPrintConditionalComparison('\tConverts to corneal irradiance %0.1f uWatts/cm2','%0.1f',cornealIrradianceUWattsPerCm2,deloriCornealIrradianceUWattsPerCm2,false);
    AnsiZ136MPEPrintConditionalComparison('\tConverts to total radiant power in the pupil of %0.2g mW','%0.2g',powerInPupilMWatts,deloriPowerInPupilMWatts,false);
        
    % Compute MPE, with comparisons when available
    [MPELimitIntegratedRadiance_JoulesPerCm2Sr, ...
        MPELimitRadiance_WattsPerCm2Sr, ...
        MPELimitCornealIrradiance_WattsPerCm2, ...
        MPELimitCornealRadiantExposure_JoulesPerCm2] = ...
        AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSeconds,stimulusDiameterDegrees,wavelengthNm,CONELIMITFLAG);
    MPELimitRadiance_WattsPerM2Sr =  1e4*MPELimitRadiance_WattsPerCm2Sr;
    
    % For comparison with Delori's calculations, we need to convert light at the cornea to retinal illuminance.  That in turn
    % requires an assumption about pupil size for whatever light is being measured.  We can use equation 8 of Delori et al. (2007, JOSA)
    % to match the assumption in Delori's calculations.
    MPEPupilFactor = AnsiZ136MPEComputePupilFactor(stimulusDurationSeconds,wavelengthNm);
    effectivePupilAreaCm2 = (pi/4)*((0.7)^2)/(MPEPupilFactor);
    effectivePupilDiameterMm = 10*sqrt(effectivePupilAreaCm2/(pi/4));
    if (abs((1e-2*(pi/4)*(effectivePupilDiameterMm^2))-effectivePupilAreaCm2) > 1e-7)
        error('Algegra boo-boo');
    end
    MPELimitPowerInPupilWatts = MPELimitCornealIrradiance_WattsPerCm2*effectivePupilAreaCm2;
    MPELimitPowerInPupilMWatts = 1e3*MPELimitPowerInPupilWatts;
    
    MPELimitRetinalIlluminanceWattsPerUm2 = RadianceToRetIrradiance(MPELimitRadiance_WattsPerM2Sr,wavelengthNm,(pi/4)*effectivePupilDiameterMm^2,ansiEyeLengthMm);
    MPELimitRetinalIlluminanceWattsPerCm2 = 1e8*MPELimitRetinalIlluminanceWattsPerUm2;
    MPELimitRetinalIlluminanceUWattsPerCm2 = 1e6*MPELimitRetinalIlluminanceWattsPerCm2;
    
    % We'll also compute some of the the factors that go into the MPE computation, for comparison with
    % what Delori's spreadsheet reports for the same numbers
    MPECb = AnsiZ136MPEComputeCb(wavelengthNm);
    MPECe = AnsiZ136MPEComputeCe(stimulusDiameterDegrees);
    MPRT2 = AnsiZ136MPEComputeT2(stimulusDiameterDegrees);
    MPECa = AnsiZ136MPEComputeCa(wavelengthNm);
    if (wavelengthNm >= 1050)
        MPECc = AnsiZ136MPEComputeCc(wavelengthNm);
    else
        MPECc = 1;
    end
    MPECt = MPECa*MPECa;
    MPET2 = AnsiZ136MPEComputeT2(stimulusDiameterDegrees);

    % Print out comparisons when available
    fprintf('\nMPE calculations\n');
    AnsiZ136MPEPrintConditionalComparison('\tUsing pupil factor %0.2f','%0.2f',MPEPupilFactor,deloriMPEPupilFactor,false);
    AnsiZ136MPEPrintConditionalComparison('\tEffective pupil diameter is %0.1f mm','%0.1f',effectivePupilDiameterMm,deloriEffectivePupilDiameterMm,false);
    AnsiZ136MPEPrintConditionalComparison('\tUsing Cb %0.2f','%0.2f',MPECb,deloriMPECb,false);
    AnsiZ136MPEPrintConditionalComparison('\tUsing Ce %0.2f','%0.2f',MPECe,deloriMPECe,false);
    AnsiZ136MPEPrintConditionalComparison('\tUsing Ct %0.2f','%0.2f',MPECt,deloriMPECt,false);
    AnsiZ136MPEPrintConditionalComparison('\tUsing T2 %0.2f','%0.2f',MPET2,deloriMPET2,false);
    AnsiZ136MPEPrintConditionalComparison('\tMPE power in pupil limit %0.3g mWatts','%0.3g',MPELimitPowerInPupilMWatts,deloriMPEPowerInPupilMWatts,false);
    AnsiZ136MPEPrintConditionalComparison('\tMPE retinal illuminance limit computed as %0.3g Watts/cm2','%0.3g',MPELimitRetinalIlluminanceWattsPerCm2,deloriMPERetinalIrradianceWattsPerCm2,false);
    fprintf('\tLimit - Stimulus log10 difference: %0.1f log10 units\n',log10(MPELimitRetinalIlluminanceWattsPerCm2)-log10(retinalIlluminanceWattsPerCm2));

    % Ready for next iteration
    fprintf('\n');

end

%% See if we can match some conversions computed
% by Ed Pugh and conveyed to me by Brian Wandell.
%
% Ed starts with retinal illuminance of 10^15 in quanta/[cm2-sec]'
%   He gets 340 uW/cm2.  We get 342.
%   He gets ~590,000 trolands.  We get about this assuming 17 mm eye length.
%   He gets ~190,000 cd/m2 for a 2mm diameter pupil.  We also get about this.




%% Check ANSI light limit calculations against numbers in Eds document.  He doesn't say
% the durations or size he assumed, but he does say he got the numbers from Delori's
% spreadsheet.  Let's try making up a source size and duration and see what happens.


% When I plug these numbers (580 nm, 2 degree stimulus, 2 mm pupil diameter, 1 second exposure)
% into the version of Delori's spreadsheet I got via Ed Pugh (rev 1/10/08), it computes that
% the stimulus has:
%   a retinal irradiance of 340 uW/cm2 
%   a radiance of 31.28 mW/[cm2 sr]
%   a corneal irradiance of 29.93 uW/cm2,
%   total radiant energy in the pupil of 940.33 nJ.
%   retinal radiant exposure of 340 uJ/cm2.
% We compute these quantities using PTB routines.  From above we already have the retinal irradiance
% at about 340 Watts/cm2. The rest of the numbers also match up pretty well.  Rounding in the spreadsheet
% plus perhaps a different assumption about eye length can explain the differences, I think.

% The spreadsheet computes the exposure safety limit for this stimulus as:
%  radiant power in the pupil 2.96 mW 
%  radiant energy in the pupil of 2.96 mJ
%  retinal irradiance 1.07 W/cm2
%  retinal exposure of 1.07 J/cm2





