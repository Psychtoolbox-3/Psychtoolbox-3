function [alpha,beta,thresh50] = FitWeibYN(inputs,nYes,nNo,alpha0,beta0,numFunCalls)
% [alpha,beta,thresh50]=FitWeibYN(inputs,nYes,nNo,[alpha0],[beta0],[numFunCalls])
%
% Fits a Weibull function to the passed yes-no data.
%
% Requires the optimization toolbox. Doesn't work with Octave yet.
%
% INPUTS:
%   inputs    Input levels
%   nYes      Number of yes responses at 
%             the corresponding input level
%   nNo       Number of no responses at 
%             the corresponding input level
%   alpha0    Initial guess for alpha (optional)
%   beta0     Initial guess for beta (optional)
% OUTPUTS:
%   alpha
%   beta
%   thresh50  50% threshold
%
% See also: FitWeibTAFC, FitFitWeibAlphTAFC, FitCumNormYN, FitLogitYN
%
% 9/15/93   jms  Added a pre-fit to get a better initial.
%           jms  Made 'options' a parameter so that printing could
%                be disabled higher up.
% 9/23/93   jms  Test the slope of the linear pre-fit to set the upper
%                and lower bounds on the fit.
% 2/5/97    dhb  Rewrote to parallel TAFC version but kept slope test.
%           dhb  Check for optimization toolbox.
% 4/18/00	  mpr	 Added an option to set the number of allowed function calls
% 10/13/00  dhb  Improve initial guess for alpha.  Thanks to Duje Tadin
% 							 for identifying the need for this.
% 3/5/05		dhb	 Update for optimization toolbox version 2.

% Determine whether function is increasing or decreasing
lineParams = [inputs ones(size(inputs))]\(nYes ./ (nYes+nNo) );
slope = lineParams(1);

% Initial parameters
if (nargin <= 3 || isempty(alpha0))
	x0(1) = mean(inputs);
else 
	x0(1) = alpha0;
end
if (nargin <= 4)
	if (slope > 0)
		x0(2) = 3.5;
	else
		x0(2) = -3.5;
	end
elseif isempty(beta0)
	x0(2) = 3.5;
else
	x0(2) = beta0;
end

% Check for needed optimization toolbox, and version.
if (exist('fminunc') == 2)
	options = optimset;
	options = optimset(options,'Diagnostics','off','Display','off');
	options = optimset(options,'LargeScale','off');
	x1 = fminunc('WeibYNFitFun',x0,options,inputs,nYes,nNo);
	x = fminunc('WeibYNFitFun',x1,options,inputs,nYes,nNo);
elseif (exist('fminu') == 2)
	options = foptions;
	state = warning; warning('off');
	x1 = fminu('WeibYNFitFun',x0,options,[],inputs,nYes,nNo);
	x = fminu('WeibYNFitFun',x1,options,[],inputs,nYes,nNo);
	warning(state);
else
	error('FitWeibYN requires the optional Matlab Optimization Toolbox from Mathworks');
end

% Extract parameters
alpha = x(1);
beta = x(2);

% Compute threshold from fit parameters
thresh50 = FindThreshWeibYN(0.5,alpha,beta);
