function photoreceptors = DefaultPhotoreceptors(kind)
% photoreceptors = DefaultPhotoreceptors(kind)
% 
% Return a structure containing default sources 
% for photoreceptor complements of various kinds.
%
% Available kinds
%   LivingHumanFovea (Default) - Human foveal cones in the eye
%   LivingHumanMelanopsinTsujimura2010 - Estimate of melanopsin gc spectral sensitivity in living eye
%   LivingDog - Canine
%   GuineaPig - Guinea pig in dish
%
% See also:  FillInPhotoreceptors, PrintPhotoreceptors, RetIrradianceToIsoRecSec
%  IsomerizationsInEyeDemo, IsomerizationsInDishDemo, ComputeCIEConeFundamentals. 
%
% NOTES: Should probably update the parameters for LivingHumanFovea so that
% they produce the Stockman-Sharpe fundamentals.  This should be pretty
% straightforward, now that all the pieces are implemented as via ComputeCIEConeFundamentals. 
%
% 7/25/03  dhb  Wrote it.
% 12/04/07 dhb  Added dog parameters
% 8/14/11  dhb  Added fieldSizeDegrees and ageInYears fields to photoreceptors for LivingHumanFovea case.
%               These defaults match the CIE standard.
% 4/20/12  dhb  Add LivingHumanMelanopsin
% 5/10/12  dhb  Changed name for LivingHumanMelanopsin to postpend Tsujimura2010
% 8/12/13  dhb  Change field order to make printouts look nicer.


% Default
if (nargin < 1 || isempty(kind))
	kind = 'LivingHumanFovea';
end

% Fill it in
switch (kind)
    case 'CIE2Deg'
        photoreceptors.species = 'Human';
        photoreceptors.types = {'FovealLCone' 'FovealMCone' 'FovealSCone'};
        photoreceptors.nomogram.S = WlsToS((390:5:780)');
		photoreceptors.OSlength.source = 'None';
		photoreceptors.ISdiameter.source = 'Rodieck';
		photoreceptors.specificDensity.source = 'None';
        photoreceptors.axialDensity.source = 'CIE';
        photoreceptors.nomogram.source = 'None';
        photoreceptors.quantalEfficiency.source = 'Generic';
        photoreceptors.fieldSizeDegrees = 2;
        photoreceptors.ageInYears = 32;
        photoreceptors.pupilDiameter.value = 3;
		photoreceptors.eyeLengthMM.source = 'Rodieck';
        photoreceptors.absorbance = 'log10coneabsorbance_ss';
		photoreceptors.lensDensity.source = 'CIE';
        photoreceptors.macularPigmentDensity.source = 'CIE';
        
    case 'CIE10Deg'
        photoreceptors.species = 'Human';
        photoreceptors.types = {'LCone' 'MCone' 'SCone'};
        photoreceptors.nomogram.S = WlsToS((390:5:780)');
		photoreceptors.OSlength.source = 'None';
		photoreceptors.ISdiameter.source = 'Webvision';
		photoreceptors.specificDensity.source = 'None';
        photoreceptors.axialDensity.source = 'CIE';
        photoreceptors.nomogram.source = 'None';
        photoreceptors.quantalEfficiency.source = 'Generic';
        photoreceptors.fieldSizeDegrees = 10;
        photoreceptors.ageInYears = 32;
        photoreceptors.pupilDiameter.value = 3;
		photoreceptors.eyeLengthMM.source = 'Rodieck';
        photoreceptors.absorbance = 'log10coneabsorbance_ss';
		photoreceptors.lensDensity.source = 'CIE';
        photoreceptors.macularPigmentDensity.source = 'CIE';
        
	case 'LivingHumanFovea'
		photoreceptors.species = 'Human';
        photoreceptors.types = {'FovealLCone' 'FovealMCone' 'FovealSCone'};
        photoreceptors.nomogram.S = [380 1 401];
		photoreceptors.OSlength.source = 'Rodieck';
		photoreceptors.ISdiameter.source = 'Rodieck';
		photoreceptors.specificDensity.source = 'Rodieck';
		photoreceptors.nomogram.source = 'StockmanSharpe';
		photoreceptors.nomogram.lambdaMax = [558.9 530.3 420.7]';
		photoreceptors.quantalEfficiency.source = 'Generic';
        photoreceptors.fieldSizeDegrees = 2;
        photoreceptors.ageInYears = 32;
        photoreceptors.pupilDiameter.source = 'PokornySmith';
		photoreceptors.eyeLengthMM.source = 'Rodieck';
        photoreceptors.lensDensity.source = 'StockmanSharpe';
		photoreceptors.macularPigmentDensity.source = 'Bone';

    % This creates Tsujiumura's (2010) estimate of the melanopsin gc
    % spectral sensitivity in the human eye. The quantal efficiency
    % is just made up, though, so that the code runs.
    %
    % Tsujimura has used different lambda-max in different papers.
    % The 482 value given here is from the 2010 paper.  His email
    % suggests he may have used 489 and 502 at different times.  Also
    % by email, he used Stockman-Sharpe not Govardovskii nomogram
    % for the 2010 paper, despite what the paper says.
    case 'LivingHumanMelanopsinTsujimura2010'
        photoreceptors.species = 'Human';
        photoreceptors.types = {'Melanopsin'};
        photoreceptors.nomogram.S = [380 1 401];
        photoreceptors.axialDensity.source = 'Tsujimura';
        photoreceptors.axialDensity.value = 0.5;
        photoreceptors.nomogram.source = 'StockmanSharpe';
		photoreceptors.nomogram.lambdaMax = [482]';
        photoreceptors.quantalEfficiency.source = 'None';
        photoreceptors.quantalEfficiency.value = 1;
        photoreceptors.fieldSizeDegrees = 10;
        photoreceptors.ageInYears = 32;
        photoreceptors.lensDensity.source = 'CIE';
		photoreceptors.macularPigmentDensity.source = 'CIE';

    case 'LivingDog'
        photoreceptors.species = 'Dog';
        photoreceptors.types = {'LCone' 'SCone' 'Rod'};
        photoreceptors.nomogram.S = [380 1 401];
        photoreceptors.OSlength.source = 'PennDog';
        photoreceptors.ISdiameter.source = 'PennDog';
		photoreceptors.specificDensity.source = 'Generic';
		photoreceptors.pupilDiameter.source = 'PennDog';
		photoreceptors.eyeLengthMM.source = 'PennDog';
		photoreceptors.nomogram.source = 'Govardovskii';
		photoreceptors.nomogram.lambdaMax = [555 429 506]';
		photoreceptors.quantalEfficiency.source = 'Generic';
        photoreceptors.lensDensity.source = 'None';
		photoreceptors.macularPigmentDensity.source = 'None';
        
	case 'GuineaPig'
        photoreceptors.species = 'GuineaPig';
        photoreceptors.types = {'MCone' 'SCone' 'Rod'};
        photoreceptors.nomogram.S = [380 1 401];
        photoreceptors.OSlength.source = 'SterlingLab';
        photoreceptors.OSdiameter.source = 'SterlingLab';
		photoreceptors.ISdiameter.source = 'SterlingLab';
		photoreceptors.specificDensity.source = 'Bowmaker';
		photoreceptors.pupilDiameter.source = 'None';
		photoreceptors.eyeLengthMM.source = 'None';
		photoreceptors.nomogram.source = 'Govardovskii';
		photoreceptors.nomogram.lambdaMax = [529 430 500]';
		photoreceptors.quantalEfficiency.source = 'Generic';
        photoreceptors.lensDensity.source = 'None';
		photoreceptors.macularPigmentDensity.source = 'None';
        
	otherwise
		error('Unknown photoreceptor kind specified');
end
