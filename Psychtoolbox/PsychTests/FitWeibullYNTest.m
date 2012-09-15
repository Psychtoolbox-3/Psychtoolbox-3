% FitWeibYNTest
%
% Simple test of our ability to fit a Weibull to YN data.
% The underlying routine requires the optimization toolbox.
%
% 2/5/97  dhb  Added comments.

% Set up inputs
maxInput = 20;
inputs = [1:1:maxInput]';
[m,n] = size(inputs);

% Cook some data
disp('Fitting an upward sloping psychometric function');
in_alpha = 6;
in_beta = 3;
nYes = round(100*ComputeWeibYN(inputs,in_alpha,in_beta));
nNo = 100-nYes;

% Fit it
[alpha,beta,thresh] = FitWeibYN(inputs,nYes,nNo);

% Make a prediction
pInputs = 0.1:0.1:maxInput;
predict = ComputeWeibYN(pInputs,alpha,beta);

figure(1);
hold off
plot(inputs,nYes./(nYes+nNo),'+');
title('Attempt to fit an upward sloping YN psychometric function')
hold on
plot(pInputs,predict,'g');
hold off

% Cook some data
disp('Fitting a downward sloping psychometric function');
in_alpha = 6;
in_beta = -2.5;
nYes = round(100*ComputeWeibYN(inputs,in_alpha,in_beta));
nNo = 100-nYes;

% Fit it
[alpha,beta,thresh] = FitWeibYN(inputs,nYes,nNo);

% Make a prediction
pInputs = 0.1:0.1:maxInput;
predict = ComputeWeibYN(pInputs,alpha,beta);

% Plot it
figure(2);
plot(inputs,nYes./(nYes+nNo),'+');
title('Attempt to fit a downward sloping YN psychometric function')
hold on
plot(pInputs,predict,'g');
hold off

