% FitWeibTAFCTest
%
% Simple test of our ability to fit a Weibull to TAFC data.
% The underlying routine requires the optimization toolbox.
%
% 8/26/94		dhb, ccc		Update for new naming conventions.

% Set up some input levels
inputs = [1:1:10]';
[m,n] = size(inputs);

% Cook some data
in_alpha = 6;
in_beta = 3;
nCorrect = round(100*ComputeWeibTAFC(inputs,in_alpha,in_beta));
nError = 100-nCorrect;

% Do the fit
[alpha,beta,thresh] = FitWeibTAFC(inputs,nCorrect,nError);

% Fit with constrained beta (at correct value).
[alpha2,beta2,thresh2] = FitWeibAlphTAFC(inputs,nCorrect,nError,[],in_beta);

% Compute predicted psychometric function
pInputs = 0.1:0.1:10;
predict = ComputeWeibTAFC(pInputs,alpha,beta);
predict2 = ComputeWeibTAFC(pInputs,alpha2,beta2);

% Make a plot
hold off
plot(inputs,nCorrect./(nCorrect+nError),'+');
hold on
plot(pInputs,predict,'r');
plot(pInputs,predict2,'g');
hold off
