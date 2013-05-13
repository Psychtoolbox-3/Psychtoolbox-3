% BitsPlusMonoPackTest
%
% Very simple test program to exercise the BitsPlusPackMonoImage routine.
%
% 3/4/07  dhb  Wrote it

% Clear
clear; close all;

% Generate a simple test image that exhaustively includes all possible
% input values.
testImage = 0:2^14-1;

% Run BitsPlusPackMonoImage and snag diagnostic info it currently returns.
[nil,inputReconstructionNewWay,inputReconstructionOldWay] = BitsPlusPackMonoImage(testImage);

figure; clf;
subplot(1,2,1); hold on
plot(testImage,inputReconstructionNewWay,'r+','MarkerSize',2);
plot(testImage,testImage,'k');
axis([0 2^14-1 0 2^14-1]); axis('square');
xlabel('Input value');
ylabel('Reconstructed value');
title('New Method');
subplot(1,2,2); hold on
plot(testImage,inputReconstructionOldWay,'r+','MarkerSize',2);
plot(testImage,testImage,'k');
axis([0 2^14-1 0 2^14-1]); axis('square');
xlabel('Input value');
ylabel('Reconstructed value');
title('Old Method');
