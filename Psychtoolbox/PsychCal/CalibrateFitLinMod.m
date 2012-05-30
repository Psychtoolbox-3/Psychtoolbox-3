function cal = CalibrateFitLinMod(cal)
% cal = CalibrateFitLinMod(cal)
%
% Fit the linear model to spectral calibration data.
%
% 3/26/02  dhb  Pulled out of CalibrateMonDrvr.
% 3/27/02  dhb  Add case of nPrimaryBases == 0.
% 2/15/10  dhb  Fix so that first basis vector is good approximation to max
%               input primary spectrum.
%          dhb  Normalize basis vectors so that their max power matches that 
%               of first component.
% 4/30/10  dhb  Execute yoked fit if yokedGamma flag is set.
% 5/25/10  dhb, ar Change yoked field names to match
% 5/26/10  dhb, ar Still fussing with names.
% 5/28/10  dhb, ar Pull out yoked fitting from here -- too confusing.
% 5/27/12  dhb     Handle case where there are more measurements than wavelength samples.

Pmon = zeros(cal.describe.S(3),cal.nDevices*cal.nPrimaryBases);
mGammaRaw = zeros(cal.describe.nMeas,cal.nDevices*cal.nPrimaryBases);
monSVs = zeros(min([cal.describe.nMeas cal.describe.S(3)]),cal.nDevices);
for i = 1:cal.nDevices
    tempMon = reshape(cal.rawdata.mon(:,i),cal.describe.S(3),cal.describe.nMeas);
    monSVs(:,i) = svd(tempMon);

    % Build a linear model
    if (cal.nPrimaryBases ~= 0)
        % Get full linear model
        [monB,monW] = FindLinMod(tempMon,cal.nPrimaryBases);

        % Express max measurement within the full linear model.
        % This is the first basis function.
        tempB = monB*monW(:,cal.describe.nMeas);
        maxPow = max(abs(tempB));

        % Get residuals with respect to first component
        residMon = tempMon-tempB*(tempB\tempMon);

        % If linear model dimension is greater than 1,
        % fit linear model of dimension-1 to the residuals.
        % Take this as the higher order terms of the linear model.
        %
        % Also normalize each basis vector to max power of first
        % component, and make sure convention is that this max
        % is positive.
        if (cal.nPrimaryBases > 1)
            residB = FindLinMod(residMon,cal.nPrimaryBases-1);
            for j = 1:cal.nPrimaryBases-1
                residB(:,j) = maxPow*residB(:,j)/max(abs(residB(:,j)));
                [nil,index] = max(abs(residB(:,j)));
                if (residB(index,j) < 0)
                    residB(:,j) = -residB(:,j);
                end
            end
            monB = [tempB residB];
        else
            monB = tempB;
        end

        % Zero means build one dimensional linear model just taking max measurement
        % as the spectrum.
    else
        cal.nPrimaryBases = 1;
        monB = tempMon(:,cal.describe.nMeas);
    end

    % Find weights with respect to adjusted linear model and
    % store
    monW = FindModelWeights(tempMon,monB);
    for j = 1:cal.nPrimaryBases
        mGammaRaw(:,i+(j-1)*cal.nDevices) = (monW(j,:))';
        Pmon(:,i+(j-1)*cal.nDevices) = monB(:,j);
    end
end

% Update calibration structure.
cal.S_device = cal.describe.S;
cal.P_device = Pmon;
cal.T_device = WlsToT(cal.describe.S);
cal.rawdata.rawGammaTable = mGammaRaw;
cal.rawdata.monSVs = monSVs;
