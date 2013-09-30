function [response] = ComputeNakaRushton(params,contrast)
% [response] =  ComputeNakaRushton(params,contrast)
%
% Compute the Naka-Rushton function on passed vector of contrasts.
% Several different forms may be computed depending on length of
% passed params vector.
%
% length(params) == 2
%   sigma = params(1)
%   n = params(2)
%   response = contrast^n/[contrast^n + sigma^n]
%
% length(params) == 3
%   Rmax = params(1)
%   sigma = params(2)
%   n = params(3)
%   response = Rmax*[contrast^n]/[contrast^n + sigma^n]
%
% length(params) == 4
%   Rmax = params(1)
%   sigma = params(2)
%   n = params(3)
%   m = params(4)
%   response = Rmax*[contrast^n]/[contrast^m + sigma^m]
% 
% 8/1/05    dhb, pr     Wrote from FitLightnessOrient version
% 8/2/07    dhb         Rewrote to allow several different forms depending
%                       on length of params.
% 12/5/10   dhb         Expanded comment.  Error check on input length
% 9/23/13   dhb         Fix BAD bug.  This wasn't computing what the comments said it was.
%                       Not sure when that crept in.  The contrast in the numerator was
%                       being divided by sigma before being raised to the power n.  No idea why.

% Extract parameter vector into meaningful variables
if (length(params) == 2)
    A = 1;
    sigma = params(1);
    exponent = params(2);
    exponent1 = params(2);    
elseif (length(params) == 3)
    A = params(1);
    sigma = params(2);
    exponent = params(3);
    exponent1 = params(3);
elseif (length(params) == 4)
    A = params(1);
    sigma = params(2);
    exponent = params(3);
    exponent1 = params(4);
else
    error('Inproper format for passed parameter vector');
end

% Check for bad contrast input
if (any(contrast < 0))
    error('Cannot deal with negative contrast');
end

% Handle really weird parameter values
if (sigma < 0 || exponent < 0 || exponent1 < 0)
    response = -1*ones(size(contrast));
else

    % Now pump the linear response through a non-linearity
    expContrast = (contrast).^exponent;
    expContrast1 = contrast.^exponent1;
    sigma1 = sigma.^exponent1;
    response = A*(expContrast ./ (expContrast1 + sigma1));
end

