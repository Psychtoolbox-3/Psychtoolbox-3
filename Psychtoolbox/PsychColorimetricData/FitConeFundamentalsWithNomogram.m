function [params,fitFundamentals,fitError] = FitConeFundamentalsWithNomogram(T_targetQuantal,staticParams,params0)
% [fitFundamentals,params,fitError] = FitConeFundamentalsWithNomogram(T_targetQuantal,staticParams,params0)
%
% Find underlying parameters that fit the passed corneal cone fundamentals
%
% 8/4/03  dhb  Wrote it.

% Doesn't work on octave due to lack of function 'fmincon' from the
% Matlab Optimization toolbox (see https://savannah.gnu.org/bugs/?35333)
% and due to use of nested function FitConesFun():
if IsOctave
    error('Sorry, this function does not yet work on GNU/Octave.');
end

% Convert initial parameter struct to parameter list
x0 = FitConesParamsToList(params0);

% Set bounds on search parameters
% Length 4, we're handling Ser/Ala polymorphism
if (length(x0) == 4)
    vlb(1:4) = 0;
    vub(1:4) = 800;
elseif (length(x0) == 3)
    vlb(1:3) = 0;
    vub(1:3) = 800;
else
    error('Unexpected length for parameter vector');
end

% Search to find best fit
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
if (exist('IsCluster') && IsCluster && matlabpool('size') > 1) %#ok<EXIST>
    options = optimset(options,'UseParallel','always');
end
x = fmincon(@FitConesFun,x0,[],[],[],[],vlb,vub,[],options);

% Convert parameter list to parameter struct and
% compute final values for return
params = FitConesListToParams(x);
fitError = FitConesFun(x);
fitFundamentals = ComputeCIEConeFundamentals(staticParams.S,staticParams.fieldSizeDegrees,staticParams.ageInYears, ...
            staticParams.pupilDiameterMM,params.lambdaMax,staticParams.whichNomogram ...
            );

    function [f] = FitConesFun(x)
        DO_LOG = 1;
        params = FitConesListToParams(x);
        T_pred = ComputeCIEConeFundamentals(staticParams.S,staticParams.fieldSizeDegrees,staticParams.ageInYears, ...
            staticParams.pupilDiameterMM,params.lambdaMax,staticParams.whichNomogram ...
            );
        
        if (DO_LOG)
            bigWeight = 1; bigThresh = -1;
            index = find(~isinf(log10(T_targetQuantal(:))));
            T_resid = log10(T_pred(index))-log10(T_targetQuantal(index));
            index1 = log10(T_targetQuantal(index)) > bigThresh;
            index2 = log10(T_targetQuantal(index)) <= bigThresh;
            if ( any(isnan(T_resid(:))) || any(isinf(T_resid(:))) )
                f = 1e6;
            else
                f = 100*(bigWeight*mean(T_resid(index1).^2) + mean(T_resid(index2).^2));
            end
        else
            T_resid = T_pred-T_targetQuantal;
            f = 100*mean((T_resid(:)).^2);
        end
    end
end

% Convert list to parameter structure
function params = FitConesListToParams(x)

% Length 4, we're handling Ser/Ala polymorphism
if (length(x) == 4)
    params.lambdaMax = x(1:4);
elseif (length(x) == 3)
    params.lambdaMax = x(1:3);
else
    error('Unexpected length for parameter vector');
end

end

% Convert parameter structure to list
function x = FitConesParamsToList(params)

if (size(params.lambdaMax,1) == 4)
    x = zeros(4,1);
    x(1:4) = params.lambdaMax;
elseif (size(params.lambdaMax,1) == 3)
    x = zeros(3,1);
    x(1:3) = params.lambdaMax;
else
    error('Unexpected number of photopigment lambda max values passed');
end

end
