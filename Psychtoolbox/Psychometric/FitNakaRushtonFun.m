function [f] = FitNakaRushtonFun(params,contrast,response)
% [f] = FitNakaRushtonFun(params,contrast,response)
% 
% Evaluate model fit and return measure of goodness of fit (f).
%
% 8/1/05    dhb, pr     Wrote it.
% 8/2/07    dhb         Get rid of silly call to ComputeNakaRushtonError.

% Unpack paramters, make predictions
prediction = ComputeNakaRushton(params,contrast);
nPoints = length(prediction);
error = prediction-response;
f = 10000*sqrt(sum(error.^2)/nPoints);

% Handle bizarre parameter values.
if (isnan(f))
	f = 2000;
end
