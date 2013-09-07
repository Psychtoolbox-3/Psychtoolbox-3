% FitConeFundamentalsTest
%
% This program allows exploration of how well we can fit
% CIE cone fundamentals with various photopigment nomograms,
% by allowing the lambdaMax to vary.
%
% Note from DHB.  My conclusion after spending some time with
% this is that the best current approach, in the context of the
% CIE fundamentals, is to use the StockmanSharpe nomogram and
% its nominal lambdaMax values.  The fit to the 2-degree CIE
% fundamentals is not perfect, but it is pretty good. I think 
% the deviations between what is produced by the nominal
% nomogram fit are probaby small compared with our certainty
% about the fundamentals (this is just intuition, but for example
% the deviations are small compared to those between the Stockman-Sharpe
% and Smith-Pokorny fundamentals).  One could try to develope a better
% nomogram, for example by taking the ser/ala polymorphism explicitly
% into account when fitting it, but I'm not sure that would win if one
% wants to keep the convenient assumption of constant photopigment 
% absorbance shape along a log wavelength axis.
%
% If one really wanted to go after fitting the fundamentals from parts,
% searching on various densities as well as the lambda max values
% would probably be the way to go.
%
% 8/11/11  dhb  Wrote it.
% 8/14/11  dhb  Clean up and add comments.
% 8/10/13  dhb  A few more notes.

% Doesn't work on octave due to lack of function 'fmincon' from the
% Matlab Optimization toolbox (see https://savannah.gnu.org/bugs/?35333)
% and due to use of nested function FitConesFun():
if IsOctave
    fprintf('Sorry, this test does not yet work on GNU/Octave.\n');
    return;
end

%% Clear
clear; close all;

%% Low end of log plot scale
lowEndLogPlot = -4;

%% Basic parameters.  For CIE functions, specify field size and age
staticParams.S = WlsToS((390:5:780)');
staticParams.fieldSizeDegrees = 2;
staticParams.ageInYears = 32;
staticParams.pupilDiameterMM = 3;

%% Get what we're trying to fit.
%
% The call to QuantaToEnergy to accomplish an energy to quanta transform is not an error.
% Rather, that routine converts spectra, but here we are converting sensitivities, so its
% the inverse we want.
targetRaw = load('T_cones_ss2');
T_targetEnergy = SplineCmf(targetRaw.S_cones_ss2,targetRaw.T_cones_ss2,staticParams.S,2);
T_targetQuantal2 = QuantaToEnergy(staticParams.S,T_targetEnergy')';
for i = 1:size(T_targetQuantal2,1)
    T_targetQuantal2(i,:) = T_targetQuantal2(i,:)/max(T_targetQuantal2(i,:));
end

%% Set nomogram type and initial lambdaMax values
% Nomogram choices are currenty:
%   'Baylor', 'Dawis', 'Govardovskii', 'Lamb', 'StockmanSharpe'
%
% You can use either three or four initial lambdaMax values.  In the
% case of four, the first two are treated as the ser/ala polymorphic
% variations, and weighted in the value hard coded in ComputeCIEConeFundamentals.
%
% Note from DHB.  For some nomograms the search will fail because they return
% values of zero inside the default range, and this screws up
% the error evaluation in the log domain.  I'm sure this could
% fixed with some fussing, but I don't want to head down that
% road right now.
%
% The StockmanSharpe nomogram seems to produce the best fits, which
% is not surprising given that it was tailored for this purpose.
%
% You can get a slightly better fit by using a ser/ala split and
% fitting four nomograms (two to the L cone) rather than 3, but the resultant
% lambdaMax values are further apart than makes sense.
staticParams.whichNomogram = 'StockmanSharpe';
params0.lambdaMax = [558.9 530.3 420.7]';

%% Get Stockman-Sharpe tabulated absorbance
load T_log10coneabsorbance_ss
T_StockmanSharpeAbsorbance = 10.^SplineCmf(S_log10coneabsorbance_ss,T_log10coneabsorbance_ss,staticParams.S,2);
lambdaMaxNominalStockmanSharpeNomogram = [558.9 530.3 420.7]';

%% Look at directly tabulated absorbance versus Stockman-Sharpe nomogram
figure; clf;
position = get(gcf,'Position');
position(3) = 1200; position(4) = 700;
set(gcf,'Position',position);
subplot(1,2,1); hold on
plot(SToWls(staticParams.S),T_StockmanSharpeAbsorbance','k','LineWidth',3);
T_nomogramAbsorbance = StockmanSharpeNomogram(staticParams.S,lambdaMaxNominalStockmanSharpeNomogram);
plot(SToWls(staticParams.S),T_nomogramAbsorbance','r','LineWidth',1);
title('S-S absorbance (blk) vs. nomogram (red)');
ylabel('Normalized absorbance');
xlabel('Wavelength');
subplot(1,2,2); hold on
plot(SToWls(staticParams.S),log10(T_StockmanSharpeAbsorbance''),'k','LineWidth',3);
T_nomogramAbsorbance = StockmanSharpeNomogram(staticParams.S,[558.9 530.3 420.7]');
plot(SToWls(staticParams.S),log10(T_nomogramAbsorbance'),'r','LineWidth',1);
ylim([lowEndLogPlot 0.5]);
title('S-S absorbance (blk) vs. nomogram (red)');
ylabel('Log10 normalized absorbance');
xlabel('Wavelength');

%% Do the fit
[params,null,fitError] = FitConeFundamentalsWithNomogram(T_targetQuantal2,staticParams,params0);
T_predictQuantal0 = ComputeCIEConeFundamentals(staticParams.S,staticParams.fieldSizeDegrees,staticParams.ageInYears, ...
    staticParams.pupilDiameterMM,params0.lambdaMax,staticParams.whichNomogram ...
    );
T_predictQuantal = ComputeCIEConeFundamentals(staticParams.S,staticParams.fieldSizeDegrees,staticParams.ageInYears, ...
    staticParams.pupilDiameterMM,params.lambdaMax,staticParams.whichNomogram ...
    );

%% Report on fit
fprintf('\nLambda max, %s nomogram\n',staticParams.whichNomogram);
if (length(params.lambdaMax) == 4)
    fprintf('\tInitial:\t%0.1f (Lser)\t%0.1f (Lala)\t%0.1f (M)\t%0.1f (S)\n', ...
        params0.lambdaMax(1),params0.lambdaMax(2),params0.lambdaMax(3),params0.lambdaMax(4));
    fprintf('\tFinal:\t\t%0.1f (Lser)\t%0.1f (Lala)\t%0.1f (M)\t%0.1f (S)\n', ...
        params.lambdaMax(1),params.lambdaMax(2),params.lambdaMax(3),params.lambdaMax(4));
else
    fprintf('\tInitial:\t%0.1f (L)\t%0.1f (M)\t%0.1f (S)\n', ...
        params0.lambdaMax(1),params0.lambdaMax(2),params0.lambdaMax(3));
    fprintf('\tFinal:\t\t%0.1f (L)\t%0.1f (M)\t%0.1f (S)\n', ...
        params.lambdaMax(1),params.lambdaMax(2),params.lambdaMax(3));
end
fprintf('Fit error = %0.3g\n',fitError);

%% Make a plot of what we got.  To try to evaluate how
% much we care about the error, the Smith-Pokorny estimates
% are in blue.  I figure that the difference between those,
% which we used happily for years, and the Stockman-Sharpe
% provides some visual sense of what a big difference is.

% Load S & P fundamentals
load T_cones_sp
T_compareSP = QuantaToEnergy(S_cones_sp,T_cones_sp')';
for i = 1:size(T_compareSP,1)
    T_compareSP(i,:) = T_compareSP(i,:)/max(T_compareSP(i,:));
end

% Plot
figure; clf; hold on
position = get(gcf,'Position');
position(3) = 1200; position(4) = 700;
set(gcf,'Position',position);
subplot(1,2,1); hold on;
plot(SToWls(staticParams.S),T_targetQuantal2','k','LineWidth',2);
plot(SToWls(staticParams.S),T_predictQuantal0','g','LineWidth',1);
plot(SToWls(staticParams.S),T_predictQuantal','r','LineWidth',1);
plot(SToWls(S_cones_sp),T_compareSP','b','LineWidth',1);

xlabel('Wavelength'); ylabel('Quantal Sensitivity');
title('S-S 2-deg fundamental (blk) vs. nomogram pred (red)');
subplot(1,2,2); hold on;
plot(SToWls(staticParams.S),log10(T_targetQuantal2'),'k','LineWidth',2);
plot(SToWls(staticParams.S),log10(T_predictQuantal0'),'g','LineWidth',1);
plot(SToWls(staticParams.S),log10(T_predictQuantal'),'r','LineWidth',1);
plot(SToWls(S_cones_sp),log10(T_compareSP'),'b','LineWidth',1);
ylim([lowEndLogPlot 0.5]);
xlabel('Wavelength'); ylabel('Log10 quantal Sensitivity');
title('S-S 2-deg fundamental (blk) vs. nomogram pred (red)');
