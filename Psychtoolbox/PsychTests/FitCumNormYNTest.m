% FitCumNormYNTest
%
% 9/22/93   jms  Created from FitWeibullYN.
% 2/8/97    dhb  Cleaned up for current calling conventions.

% Upward Sloping Psychometric Function
disp('Fitting an upward sloping psychometric function');

% Set input values
maxInput = 100;
inputs = [1:1:maxInput]';
[m,n] = size(inputs);

% Set parameters and generate some data
uIn =   20;
varIn = 5;
nYes = round(100*NormalCumulative(inputs,uIn,varIn));
nNo = 100-nYes;

% Fit and generate prediction
[uEst,varEst] = FitCumNormYN(inputs,nYes,nNo);
pInputs = [1:maxInput];
predict = NormalCumulative(pInputs,uEst,varEst);

% For comparision, fit a logistic to the same data
[a,b,threshLogit] = FitLogitYN(inputs,nYes,nNo);
predictLogit = ComputeLogistic(pInputs,a,b);

% Print comparision of fit normal and logistic thresholds
fprintf('Normal fit finds threshold at %g, logit at %g\n',uEst,threshLogit);

% Make a plot fo the output
hold off
plot(inputs,nYes./(nYes+nNo),'+');
title('Fit to upward sloping YN psychometric function')
hold on
plot(pInputs,predict,'g');
plot(pInputs,predictLogit,'r');
hold off
