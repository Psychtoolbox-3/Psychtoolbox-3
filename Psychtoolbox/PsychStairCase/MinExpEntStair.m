function fhndl = MinExpEntStair
% Minimum Expected Entropy Staircase
%
% The staircase gives suggestions for which probe value to test next,
% choosing the probe that will provide the most information (based on the
% principle of minimum entropy = maximally unambiguous probability
% distribution). Probes are chosen from a set of possible probe values
% provided at staircase init, and their use is evaluated based on the
% expected amount of information gain given a space of PSE and slope values
% to test over.
%
% By default, a psychometric function ranging from 0% to 100% is used, as
% is suitable for discrimination experiments with a standard in the middle
% of the possible stimulus parameter range. For other paradigms, such as
% n-AFC detection tasks, one can set the guessrate input during staircase
% init to 1/num_alternatives, e.g. .5 when doing a 2IFC detection task.
% This guess rate is thus not the rate at which participants guess instead
% of do your task (thats the lapse rate), it the minimum rate of correct
% responses as determined by your design. NB: below discussion is based on
% the default psychometric function with the full range, but all points are
% equally valid for a scaled psychometric function.
%
% It is recommended to have the staircase determine the optimal next probe
% based on only a random subset of the response history (see options
% 'toggle_use_resp_subset' and 'toggle_use_resp_subset_prop'). This makes
% its operation more robust for response errors and also avoids probe
% oscillations when the fit estimate is converging.
% When we are close to convergence, probes will tend to be near the 25% and
% 75% points. If a probe is 25% and you answer '1' (pedestal faster, which
% is likely, because it's near the correct 25% point), then for the next
% trial the peak in expected entropy reduction will generally be the 75%
% point, and vice versa. This can lead to undesirable probe sequences where
% the correct response alternates 0,1,0,1,0,1. If you choose a random
% subset, this will completely eliminate the problem. If the staircase has
% converged to where there are two almost equal expected entropy minima,
% then small variations due to the selection of subsets will randomly vary
% which minimum emerges as lowest.
% This strategy does not significantly affect optimal operation of the
% staircase. Lots of probe values provide useful information. Therefore, it
% is not crucial to have a highly accurate estimate of likelihoods, so
% relatively few trials are sufficient (less than are needed to for final
% estimates of PSE and DL). Throwing out trials for the staircase
% computation yields robustness without much cost.
%
% Another option would be to load a non-uniform prior on the space of
% possible location/mean/PSE and dispersion/slope parameters (known as mu
% and sigma respectively for a cumulative Gaussian - see option
% 'loadprior'). Probe sampling will then stay reasonable in early trials
% even if there were a couple bad responses. But this strategy is not as
% robust as using a random subset -- bad trials will continue to have an
% effect throughout.
%
% In absence of anything to base the optimal probe value on, the first
% probe is chosen randomly from the set of possible probes. When a prior
% was loaded, a likelihood distribution is available based on which the
% optional probe value can be computed. If for any other reason choosing
% the next probe based on the measure of minimum expected entropy fails,
% the staircase will fall back on the same random probe sampling strategy.
% There is an option to set the first probe value to be tested, which, for
% the first trial only, will overrule both of the above probe choice
% strategies. This can be useful if you want to be sure that the first
% trial is an easy one so the participant knows what to expect.
%
% Another measure for robustness is to choose a small lapse rate. If lapse
% would be zero and a response error is made by the observer, immediately a
% whole range of mean-slope combinations becomes impossible. If lapse rate
% is non-zero, these would still have a non-zero probability and the
% staircase can rebound. Therefore a lapse rate of 5% or even more
% depending on task difficulty is always recommended. NB: in the default
% discrimination setup of the staircase (guessrate is not specified or set
% to 0), half of the lapse rate is taken off the bottom of the psychometric
% function and half is taken off the top. So if the lapse rate is 0.05, the
% psychometric function will range from 0.025 to 0.975. In the setup for a
% n-AFC detection experiment when the psychometric function has a lower
% bound of 1/num_alternatives, the lapse rate is taken off the top. So when
% the guess_rate is set to .5 (2AFC) and the pase rate is set to .05, the
% psychometric function will range from 0.05 to 0.95.
% Note that the staircase does not support a 0 lapse rate in the first
% place as it works with log-probability and we get in trouble if we would
% take the log of a 0 probability. Any lapse rate lower than 1e-10 will be
% adjusted to 1e-10 upon calling the 'init' function.
%
% If the staircase gets stuck at one of the bounds of the probe set, check
% that the sign of the slope space matches the expected sign of the
% response. E.g., lets look at an experiment in which you are doing 2IFC
% task in which the observer is asked to report which interval contained
% the faster motion. If the observer choses the test over the pedestal
% interval the response is 1, if the observer chosen the pedestal to be
% faster, the response is 0. All slopes in the set would in this case be
% positive as the low end of the probe space (slow speeds) is associated
% with response 0 and the high end with response 1. If we however asked the
% observer to indicate the slower interval, the slopes in our slope set
% would not match the task, and the staircase would get stuck at one of the
% probe bounds. In this case, the lower end of the probe space is
% associated with the response 1 and the higher end with the response
% 0--we'd thus have a negative slope for the fitted cumulative probability
% function.
%
% The staircase currently only supports logistic and cumulative Gaussian
% (default) psychometric functions (see 'set_psychometric_func'), but
% others could easily be implemented. Changes should be needed only to the
% function "fit_a_point" at the bottom of this mfile, providing that the
% function is characterized by two parameters (which do not necessarily
% have to be PSE and slope, though that is the terminology here.
% Should you implement such a function, please do send me your code at
% dcnieho @at@ gmail.com.
%
% The above discussion assumes that response inputs to 'process_resp' are
% either 0 or 1 (see note above about their meaning) though in practice
% anything larger than 0 is treated as 1 and anything lower than 0,
% including 0, is treated as 0. the staircase can thus easily be integrated
% with programs that use a 1, -1 response scheme.
%
% For actual offline fitting of your data, you would probably want to use a
% dedicated toolbox such as Prins, N & Kingdom, F. A. A. (2009) Palamedes:
% Matlab routines for analyzing psychophysical data.
% http://www.palamedestoolbox.org. instead of using the function parameters
% or PSE and DL returned from staircase functions 'get_fit' and
% 'get_PSE_DL'.
% Also note that while the staircase runs far more robust when a small
% lapse rate is assumed, it is common to either fit the psychometric
% function without a lapse rate, or otherwise with the lapse rate as a free
% parameter (possibily varying only over subjects, but not over conditions
% within each subject).
%
%
% References:
%  Based on the Minimum expected entropy staircase method developed by:
%  Saunders JA & Backus BT (2006). Perception of surface slant from
%    oriented textures. Journal of Vision 6(9), article 3
%
%  Discussions of conceptually similar staircases can be found in:
%  Kontsevich LL & Tyler CW (1999). Bayesian adaptive estimation of
%    psychometric slope and threshold. Vision Res 39(16), pp. 2729-37
%  Lesmes LA, Lu ZL, Baek J & Albright TD (2010). Bayesian adaptive
%    estimation of the contrast sensitivity function: The quick CSF method.
%    Journal of Vision 10(3), article 17
%
%
% USE:
% Calling this function creates a staircase instance. The interface of the
% staircase is accessed through the returned function handle. You can
% create as many instances as you like by calling this function, each
% instance has its own internal memory/history. In that sense this is
% really OO (I'm not happy with MATLAB's OO features and also want to be
% compatible with old versions, hence the below paradigm).
% When interacting with the staircase through the function handle, the
% first argument is a string that identifies the action you want to perform
% (you can think of this as the string containing the name of the member
% function to be called) and optionally any other arguments that are needed
% for the call. See MESDemo for an example and the comments below for use
% of the different staircase functions.

% Copyright (c) 2011 by DC Niehorster and JA Saunders

% private member variables
probeset    = [];           % possible probe values to be tested
aset        = [];           % pse's tested (and fitted)
bset        = [];           % slopes fitted
agrid       = [];
bgrid       = [];
lapse_rate  = [];           % lapse/mistake rate
guess_rate  = [];           % guess rate
phist       = [];           % probe history
rhist       = double([]);   % response history (0 or 1)
loglik      = [];
lik         = [];
g0          = [];
g1          = [];
g2          = [];

% option: use a subset of all data for choosing the next probe, default values:
quse_subset = false;        % use limited subset for computing next probe? Limited subset by discarding a fixed number of trials
quse_subset_perc = false;   % same as above, but instead use a percentage of the available data
minsetsize  = 10;           % minimum size to start subsetting
subsetsize  = 3;            % subset contains subsetsize less datapoints than full dataset
percsetsize = .8;           % percentage of data in set used

% option: set the value to test if probe history is empty
first_value = [];           % first value to test instead of random or by prior

% psychometric function that is used (default)
psychofunc  = 'cumGauss';

% subfunction
fhndl = @MinExpEntStair_internal;

% public interface
    function [varargout] = MinExpEntStair_internal(mode,varargin)
        
        switch mode
            %%% init
            case 'init' % [] = stair('init',probeset,meanset,slopeset,lapse_rate,guess_rate);
                probeset            = varargin{1};
                aset                = varargin{2};
                bset                = varargin{3};
                [agrid,bgrid]       = meshgrid(aset,bset);
                % init with uniform probability, normalized
                loglik              = zeros(size(agrid)) - log(numel(agrid));
                lik                 = ones(size(agrid))./numel(agrid);
                % lapse rate and guess rate
                lapse_rate          = varargin{4};
                % the lapse rate cannot be exactly 0 as the computed
                % probability must not be exactly 0 so we can work with
                % log(prob) without trouble, so set it to 1e-10 at least.
                lapse_rate          = max(lapse_rate,1e-10);
                % guess rate is optional, if not specified we assume a 2IFC
                % discrimination experiment where the guess rate is
                % irrelevant as function goes from always one option at the
                % one end to always the other option at the other end.
                if length(varargin)<5
                    guess_rate = 0;
                else
                    guess_rate = varargin{5};
                end
                
                % lapse rate:
                % 1. for a discrimination setup (guess_rate==0) the
                % lapserate basically means that instead of ranging from 0
                % to 1, the psychometric function ranges from lapse_rate/2
                % to 1-lapse_rate/2
                % 2. for a detection setup, the lower bound is guess_rate
                % and the upper bound is 1-lapse_rate
                
                % lower bound of pyschometric function
                % and
                % range of pyschometric function
                if guess_rate==0
                    g0 = lapse_rate/2;
                    g1 = 1 - lapse_rate;
                else
                    g0 = guess_rate;
                    g1 = 1 - lapse_rate - guess_rate;
                end
                g2 = 1 - g0;    % need to flip psychometric function for fitting responses <= 0, get upper bound of this flipped function
                
                
            %%% load bunch of previously run trials (need probes and
            %%% responses)
            case 'loadhistory' % [] = stair('loadhistory',probes,responses);
                phist               = varargin{1};
                rhist               = varargin{2};
                
                % refit likelihood up to this point
                [loglik,lik]        = fit_all(phist,rhist);
                
                
                
            %%% load a prior likelihood, so that first probe is not chosen
            %%% randomly and you can influence evolution of the fit
            case 'loadprior' % [] = stair('loadprior',priorlik);
                assert(all(loglik(:)==-log(numel(agrid))),'Cannot load prior if we have a likelihood already'); % this tests if it is not default inited
                
                priorlik = varargin{1};
                assert(size(priorlik,1)==length(bset),'Number of rows in prior much match length of slope set')
                assert(size(priorlik,2)==length(aset),'Number of columns in prior much match length of mean set')
                assert(all(priorlik(:)>=0),'Loaded prior is not expected to be a log likelihood (that is: all your probabilities should be larger than or equal to 0!)');
                
                loglik  = normalize_loglik(log(priorlik));
                lik     = exp(loglik);
                
                
            %%% use subset of data for computing next probe
            case 'toggle_use_resp_subset' % [] = stair('toggle_use_resp_subset',20,6);
                % option: extract a probe and response subset for choosing
                % the next probe, and fit just those
                % when lots of trials ran, entropy function often has two
                % local minima, with their relative values switch back and
                % forth. This will lead to large oscillations in the probe
                % value being tested (one trial a probe from the beginning
                % of set, next trial a probe from the end and the from
                % beginning of set again).
                % We want to avoid these oscillations in probe values,
                % therefore we select a limited subset of data to calculate
                % the best next probe.
                quse_subset  = ~quse_subset;
                assert(~(quse_subset && quse_subset_perc));
                if ~isempty(varargin) % change defaults
                    minsetsize = varargin{1};
                    subsetsize = varargin{2};
                end
                varargout{1} = quse_subset;
                varargout{2} = minsetsize;
                varargout{3} = subsetsize;
            
                
            %%% use subset of data for computing next probe
            case 'toggle_use_resp_subset_prop' % [] = stair('toggle_use_resp_subset_perc',10,.8);
                % same as above, but now always use a proportion of the
                % available data
                quse_subset_perc  = ~quse_subset_perc;
                assert(~(quse_subset_perc && quse_subset));
                if ~isempty(varargin) % change defaults
                    minsetsize  = varargin{1};
                    percsetsize = varargin{2};
                end
                varargout{1} = quse_subset_perc;
                varargout{2} = minsetsize;
                varargout{3} = percsetsize;
                
                
            % set the first value to test. Normally the first is chosen
            % randomly or by using the prior that you loaded. If you prefer
            % to start at a fixed value, use this.
            case 'set_first_value' % [] = stair('set_first_value',first_value);
                first_value = varargin{1};
                if ~isempty(phist)
                    warning('the first trial has already been run. Setting the first value now is pointless and it''ll be ignored');
                end
                
                
            % set the psychometric function to be used (default cumulative
            % Gaussian). Can be called at any time (but it will refit all
            % the data already present and thus remove the effect of any
            % priors).
            case 'set_psychometric_func' % [] = stair('set_psychometric_func','funcID');
                % currently supported:
                %  'cumGauss' - Cumulative Gaussian
                %  'logistic' - logistic function
                psychofunc = varargin{1};
                % if there's any data already, refit it using the new
                % psychometric func. This would remove the effect of any
                % priors!
                if ~isempty(phist)
                    ndata = min(length(phist),length(rhist));
                    [loglik,lik]    = fit_all(phist(1:ndata),rhist(1:ndata));
                end
                
                
            % get the psychometric function that is currently used.
            case 'get_psychometric_func' % ['funcID'] = stair('get_psychometric_func');
                % currently possible outputs:
                %  'cumGauss' - Cumulative Gaussian
                %  'logistic' - logistic function
                varargout{1} = psychofunc;
                
                
            %%% given history, get which probe is best to test next
            case 'get_next_probe' % [probe,entexp,ind]  = stair('get_next_probe');
                if isempty(phist) && ~isempty(first_value)
                   % first trial and user requested a specific probe value to be tested
                   p                = first_value;
                   [varargout{2:3}] = deal([]);
                else
                    [p,entexp,indmin]   = getnextprobe;
                    if isempty(p) || isscalar(unique(loglik))
                        % if we couldn't compute expected entropy, or we have a
                        % uniform likelihood on which calculation was based
                        % (useless prior info, such as default inited), fall
                        % back on random probe selection
                        p                   = probeset(round(RandLim(1,1,length(probeset))));
                        [varargout{2:3}]    = deal([]);
                    else
                        varargout{2}        = entexp;
                        varargout{3}        = indmin;
                    end
                end
                varargout{1}    = p;
                phist           = [phist p];
            
                
            %%% fit likelihoods for new response
            case 'process_resp' % [] = stair('process_resp',resp); - resp on current trial
                rhist(end+1)    = varargin{1};
                [loglik,lik]    = fit_additional_data_point(loglik,phist(end),rhist(end));
                
                
            %%% retrieve probe and response history
            case 'get_history' % [probes,responses] = stair('get_history');
                varargout{1}    = phist;
                varargout{2}    = rhist;
                
            
            %%% get fitted a (PSE) and b (slope) parameters and loglik.
            %%% This returns the fit of all data, also when subsetting is
            %%% enabled.
            case 'get_fit' % [a,b,loglik]    = stair('get_fit');
                kmin            = find(loglik == max(loglik(:))); % most likely combination(s) of PSE and Slope
                varargout{1}    = mean(agrid(kmin));
                varargout{2}    = mean(bgrid(kmin));
                varargout{3}    = loglik;
                
            
            %%% get fitted PSE and DL (distance of 75% point from the 50%
            %%% point) and loglik. This returns the fit of all data, also
            %%% when subsetting is enabled.
            %%% This function is meant to be used for discrimination
            %%% experiments only (hence the terminology), although it will
            %%% return the inflection point and the distance between the
            %%% points that are equivalent to the 50% and 75% points after
            %%% scaling the psychometric function for all setups.
            case 'get_PSE_DL' % [PSE,DL,loglik]    = stair('get_PSE_DL');
                [varargout{1:3}] = MinExpEntStair_internal('get_fit');
                % convert b (dispersion) parameter to DL
                switch psychofunc
                    case 'cumGauss'
                        varargout{2} = varargout{2} * erfinv(.5)*sqrt(2);
                    case 'logistic'
                        varargout{2} = varargout{2} * log(3);
                    otherwise
                        error('Psychometric function "%s" not supported',psychofunc);
                end
                
            otherwise
                error('MinExpEntStair: mode "%s" unknown',mode);
        end
    end


% helpers (private functions, can only be called from the public
% MinExpEntStair_internal())
    function [p,entexp,indmin] = getnextprobe
        if length(rhist)>minsetsize && (quse_subset || quse_subset_perc)
            % select subset and fit
            if quse_subset_perc
                ind = NRandPerm(length(rhist),round(length(rhist)*percsetsize)); % select percentage of set
            else
                ind = NRandPerm(length(rhist),length(rhist)-subsetsize); % select set minus a few data points
            end
            [thellik,thelik] = fit_all(phist(ind),rhist(ind));
        else
            % use likelihoods already fitted for all available data
            thelik  = lik;
            thellik = loglik;
        end
        
        entexp  = zeros(1,length(probeset));
        for ksamp = 1:length(probeset)
            % probe value to process in this iteration
            xsamp       = probeset(ksamp);
            
            % p values for each possible model
            % these are used in multiple steps
            pvalsamp    = fit_a_point(xsamp,1);
            
            % expected value is sum, weighted by lik
            pval        = sum(pvalsamp(:).*thelik(:));
            
            % two possibilities for next response, 0 or 1
            % each would make a diff new likelihood function
            newloglik0  = thellik(:) + log(1 - pvalsamp(:));
            newloglik1  = thellik(:) + log(    pvalsamp(:));
            
            % important!  need to normalize
            newloglik0  = normalize_loglik(newloglik0);
            newloglik1  = normalize_loglik(newloglik1);
            
            % 0 and 1 for next response each has an entropy
            ent0        = sum(-exp(newloglik0).*newloglik0);
            ent1        = sum(-exp(newloglik1).*newloglik1);
            
            % probability pval of 0, probability (1-pval) of 1
            % use these to get expected value of entropy
            entexp(ksamp)  = ent0*(1-pval) + ent1*pval;
        end
        
        indmin      = find(entexp == min(entexp),1);
        p           = probeset(indmin);        
    end

    function [loglik,lik] = fit_additional_data_point(loglik,probe,resp)
        % get likelihood of current point
        currlik = fit_a_point(probe,resp);
        % multiply with previous likelihoods
        loglik  = loglik + log(currlik);
        
        % normalize
        loglik  = normalize_loglik(loglik);
        lik     = exp(loglik);
    end

    function [loglik,lik] = fit_all(probes,resps)
        
        if length(probes) ~= length(resps)
            error('Number of probe values and responses does not match');
        end
        
        if strcmp(psychofunc,'cumGauss')
            % we have a fast one for this!
            loglik = FitCumGauss_MES(probes,resps,aset,bset,lapse_rate,guess_rate);
        else
            
            loglik = zeros(size(agrid));
            for p=1:length(probes)
                loglik = fit_additional_data_point(loglik,probes(p),resps(p));
            end
        end
        
        % normalize
        loglik  = normalize_loglik(loglik);
        lik     = exp(loglik);
    end

    function pval = fit_a_point(probe,resp)
        switch psychofunc
            case 'cumGauss'
                if resp > 0
                    pval = g0 + g1*normcdf((probe-agrid)./bgrid);
                else
                    pval = g2 - g1*normcdf((probe-agrid)./bgrid);
                end
                
                % this reduces to:
                % if resp > 0
                %     pval =       normcdf( (probe-agrid)./bgrid);
                % else
                %     pval = 1.0 - normcdf( (probe-agrid)./bgrid);
                % end
                % when lapse_rate and guess_rate are 0
                %
                %        1  [             x - a      ]
                %   P = --- [ 1 + erf( ----------- ) ],
                %        2  [           b*sqrt(2)    ]
                % where a and b are known as the mean (mu) and the standard
                % deviation (sigma)
                % http://en.wikipedia.org/wiki/Normal_distribution
                
            case 'logistic'
                if resp > 0
                    pval = g0 + g1./(1+exp(-(probe-agrid)./bgrid));
                else
                    pval = g2 - g1./(1+exp(-(probe-agrid)./bgrid));
                end
                
                % this reduces to:
                % if resp > 0
                %     pval =       1./(1+exp(-(probe-agrid)./bgrid));
                % else
                %     pval = 1.0 - 1./(1+exp(-(probe-agrid)./bgrid));
                % end
                % when lapse_rate and guess_rate are 0
                %
                %               1
                % P =   ------------------,
                %              -(x - a)/b 
                %        1 + e^
                %
                % where a and b are known as the mean (mu) and b is
                % proportional to the standard deviation (s)
                % http://en.wikipedia.org/wiki/Logistic_distribution
                
            otherwise
                error('Psychometric function "%s" not supported',psychofunc);
        end
    end

    function loglik  = normalize_loglik(loglik)
        loglik  = loglik - log(sum(exp(loglik(:))));
    end
end
