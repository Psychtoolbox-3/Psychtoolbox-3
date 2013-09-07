function [alpha,beta,thresh92] = FitWeibTAFC(inputs,nCorrect,nError,...
alpha0,beta0)
% [alpha,beta,thresh92] = FitWeibTAFC(inputs,nCorrect,nError,...
%	  [alpha0],[beta0])
%
% Maximum likelihood fit of a Weibull function to TAFC psychometric data.
%
% Requires the optimization toolbox. Doesn't work with Octave yet.
%
% INPUTS:
%   inputs:   Contains the input levels
%   nCorrect: Contains the number of yes responses at 
%             the corresponding input level
%   nError:   Contains the number of no responses at 
%             the corresponding input level
%  alpha0:    Initial guess for alpha (optional)
%  beta0:     Initial guess for beta (optional)
% OUTPUTS:
%  alpha      Weibull alpha parameter
%  beta       Weibull beta parameter
%  thresh92   92% percent correct threshold
% 
% See also: FitWeibAlphTAFC, FitWeibYN, FitCumNormYN, FitLogitYN
%
% 8/25/94   dhb, ccc    Cleaned comments, return 92% correct threshold
% 2/5/97    dhb         Check if fminu is not available.
%                       Add slope test.
% 4/26/97   dhb         Fix bug in threshold assignment
% 10/13/00  dhb         Improve initial guess for alpha.  Thanks to Duje Tadin
%                       for identifying the need for this.
% 4/18/02   dhb         Suppress warnings in calls to optimization toolbox.
% 3/5/05    dhb         Update for optimization toolbox version 2.

% Determine whether function is increasing or decreasing
lineParams = [inputs ones(size(inputs))]\(nCorrect ./ (nCorrect+nError) );
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
elseif (isempty(beta0) )	
	x0(2) = 3.5;
else
	x0(2) = beta0;
end

% Check for needed optimization toolbox, and version.
if (exist('fminunc') == 2)
	options = optimset;
	options = optimset(options,'Diagnostics','off','Display','off');
	options = optimset(options,'LargeScale','off');
	x1 = fminunc('WeibTAFCFitFun',x0,options,inputs,nCorrect,nError);
	x = fminunc('WeibTAFCFitFun',x1,options,inputs,nCorrect,nError);
elseif (exist('fminu') == 2)
	options = foptions;
	state = warning; warning('off');
	x1 = fminu('WeibTAFCFitFun',x0,options,[],inputs,nCorrect,nError);
	x = fminu('WeibTAFCFitFun',x1,options,[],inputs,nCorrect,nError);
	warning(state);
else
	error('FitWeibTAFC requires the optional Matlab Optimization Toolbox from Mathworks');
end

% Extract parameters
alpha = x(1);
beta = x(2);

% Compute threshold from fit parameters
thresh92 = FindThreshWeibTAFC(0.92,alpha,beta);
