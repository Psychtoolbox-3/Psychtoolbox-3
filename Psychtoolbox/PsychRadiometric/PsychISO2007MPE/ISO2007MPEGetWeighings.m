function [wls,weightingR,weightingA,weightingS,wls_R,rawWeightingR,wls_A,rawWeightingA,wls_S,rawWeightingS] = ISO2007MPEGetWeighings(S)
% [wls,weightingR,weightingA,weightingS,wls_R,rawWeigtingR,wls_A,rawWeightingA,wls_S,rawWeightingS] = ISO2007MPEGetWeighings(S)
%
% Read the text files and get the three weighting functions needed for the
% standard's calculations.
%
% These are splined to the evenly spaced wavelengths specified by S, and
% padded with 0 outside the wavelength range over which each function is
% specified in the standard.  This simplifies calculations.
%
% ****************************************************************************
% IMPORTANT: Before using the ISO2007MPE routines, please see the notes on usage
% and responsibility in PsychISO2007MPE/Contents.m (type "help PsychISO2007MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% 6/25/13  dhb  Wrote it.

wls = SToWls(S);
tableA1 = dlmread('ISO2007MPETableA1.txt','\t',1,0);
wls_R = tableA1(:,1);
wls_A = tableA1(:,1);
rawWeightingR = tableA1(:,2);
rawWeightingA = tableA1(:,3);
weightingR = interp1(wls_R,rawWeightingR,wls,'linear',0);
weightingA = interp1(wls_A,rawWeightingA,wls,'linear',0);

tableA2 = dlmread('ISO2007MPETableA2.txt','\t',1,0);
wls_S = tableA2(:,1);
rawWeightingS = tableA2(:,2);
weightingS = interp1(wls_S,rawWeightingS,wls,'linear',0);

