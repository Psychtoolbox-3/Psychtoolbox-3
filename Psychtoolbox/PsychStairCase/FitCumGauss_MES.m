function [loglik,m,d] = FitCumGauss_MES(p,r,mset,dset,lapserate,guessrate)
% Fits a cumulative Gaussian to a set of probe values and subject responses
% 
% [loglik,m,d] = FitCumGauss_MES(p,r,mset,dset,lapserate,guessrate)
% Computes likelihood of each combination of PSEs mset and slopes dset
% given a set of probe values and responses using a cumulative Gaussian.
%
% For ease of thinking about it you can view responses as either 0 or 1,
% though in practice anything larger than 0 is treated as 1 and anything
% lower than 0, including 0, is treated as 0. A 1, -1 response scheme as
% input is thus no problem.
%
% By default, a psychometric function ranging from 0% to 100% is used, as
% is suitable for discrimination experiments with a standard in the middle
% of the possible stimulus parameter range. For other paradigms, such as
% detection tasks, one can set the guessrate input to 1/num_alternatives,
% e.g. .5 when doing a 2IFC detection task.
% 
% Optinally fits with a lapse rate, which defaults to 0. If a lapse rate is
% set, the cumulative Gaussian levels off at lapserate/2 and 1-lapserate/2
% when the guessrate is 0. If the guessrate is non zero, the cumulative
% Gaussian that is fit ranges from guessrate to 1-lapserate

% Copyright (c) 2011 by DC Niehorster and JA Saunders

if nargin < 5
    lapserate = 0;
end
if nargin < 6
    guessrate = 0;
end

if guessrate==0
    g0 = lapserate/2;
    g1 = 1 - lapserate;
else
    g0 = guessrate;
    g1 = 1 - lapserate - guessrate;
end
g2 = 1 - g0;

[msamp,dsamp]   = meshgrid(mset,dset);
sz              = size(msamp);
msamp           = msamp(:).';
dsamp           = dsamp(:).';

% the below are equivalent ways of computing the fit, lowest is fastest but
% most obtuse, the others are provided for documentation purposes
if 0
    
    % slow way with a loop,
    % but less memory intensive
    
    loglik = zeros(size(msamp));
    for ksamp = 1:length(msamp)
        for ktrial = 1:length(r)
            if(r(ktrial) > 0)
                currlik = log(g0 + g1*normcdf( (p(ktrial)-msamp(ksamp))/dsamp(ksamp)));
                % reduces to:
                % currlik = log(      normcdf( (p(ktrial)-msamp(ksamp))/dsamp(ksamp)));
                % when lapserate is 0
            else
                currlik = log(g2 - g1*normcdf( (p(ktrial)-msamp(ksamp))/dsamp(ksamp)));
                % currlik = log(1.0 - normcdf( (p(ktrial)-msamp(ksamp))/dsamp(ksamp)));
            end
            loglik(ksamp) = loglik(ksamp) + currlik;
        end
    end
    
elseif 0
    % faster way by blocking, but still looping over trials
    % choose this one if the below get you into lack of memory trouble
    loglik = zeros(size(msamp));
    for ktrial = 1:length(r)
        if(r(ktrial) > 0)
            currlik = log(g0 + g1*normcdf( (p(ktrial)-msamp)./dsamp) );
            % reduces to:
            % currlik = log(      normcdf( (p(ktrial)-msamp)./dsamp));
            % when lapserate is 0
        else
            currlik = log(g2 - g1*normcdf( (p(ktrial)-msamp)./dsamp) );
            % currlik = log(1.0 - normcdf( (p(ktrial)-msamp)./dsamp));
        end
        loglik     = loglik + currlik;
    end
    
elseif 0
    % faster way by blocking,
    % but more memory intensive
    rmat        = repmat(r(:)>0,[1,length(msamp)]);
    mmat        = repmat(msamp,[length(r),1]);
    dmat        = repmat(dsamp,[length(r),1]);
    pmat        = repmat(p(:),[1,length(msamp)]);
    
    temp        = g2 - g1*normcdf( (pmat-mmat)./dmat);
    ind         = find(rmat(:));
    temp(ind)   = 1 - temp(ind);
    loglik      = sum(log(temp),1);
    
else
    % slightly faster and less memory intensive way than previous, using
    % bsxfun (very slightly actually, but hey)
    temp        = g2 - g1*normcdf( bsxfun(@rdivide,bsxfun(@minus,p(:),msamp),dsamp));
    rmat        = repmat(r(:)>0,[1,length(msamp)]);
    temp(rmat)  = 1 - temp(rmat);
    loglik      = sum(log(temp),1);
end

[~,kmax]= max(loglik);
m       = mean(msamp(kmax));
d       = mean(dsamp(kmax));

loglik  = reshape(loglik,sz);

%imagesc(exp(0.5*loglik))
