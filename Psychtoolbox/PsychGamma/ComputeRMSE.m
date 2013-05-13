function rmse = ComputeRMSE(data,predict,SUPRESS_WARNING)
% rmse = ComputeRMSE(data,predict,SUPRESS_WARNING)
%
% Compute a root fractional SSE between data and prediction.
% Inputs should be column vectors.
% Actual code is:
%   diff = predict-data;
%   rmse = sqrt((diff'*diff)/(data'*data));
%
% The routine badly named, because what it computes
% is not what anyone would call an RMSE.  A better
% name for the routine would be ComptueFRSSE or something
% like that.
%
% 2/3/96   dhb  Added improved comments.
% 1/13/13  dhb  Added cautionary comment about what this routine does.
%               as well as message that would tell an unsuspecting user
%               at runtime.
%          dhb  Added mechanism to suppress the warning message.

% Warning about misnomer of routine.
if (nargin < 3 || isempty(SUPRESS_WARNING))
    SUPRESS_WARNING = 0;
end
if (~SUPRESS_WARNING)
    fprintf('WARNING: Although this routine computes an squared error measure\n');
    fprintf('it is not what anyone would call the RMSE.  I don''t know what I\n');
    fprintf('thinking in 1996 when I wrote this.\n');
    fprintf('\n');
    fprintf('Although the routine is badly named, it does no seem like a good idea\n');
    fprintf('to change its behavior silently, since someone may be happily using it\n');
    fprintf('to do what it actually does.\n');
    fprintf('\n');
    fprintf('Call with additional argument SUPRESS_WARNING set to true to\n');
    fprintf('supress this error message.\n');
    fprintf('\n');
    fprintf('- David Brainard, 13 Jan 2013.\n');
end


diff = predict-data;
rmse = sqrt((diff'*diff)/(data'*data));
