function cal = CalibrateFitYoked(cal)
% cal = CalibrateFitYoked(cal)
%
% Fit the gamma data from the yoked measurements.
%
% This has to do with Brainard lab HDR display calibration procedures
% and doesn't do anything unless some special fields exist
% in the calibration structure.  It's in the PTB because when we use
% it, we want to call it from script RefitCalGamma, and that one does
% belong in the PTB.  
%
% 4/30/10  dhb, kmo, ar  Wrote it.
% 5/24/10  dhb           Update comment.
% 5/25/10  dhb, ar       New input format.
% 5/28/10  dhb, ar       Execute conditionally.
% 6/10/10  dhb           Make sure returned gamma values in range an monotonic.
% 5/26/12  dhb           Improve the comment so this is a little less weird.

%% Debugging switch
DEBUG = 0;

%% Check that this is possible
OKTODO = 1;
if (~isfield(cal.describe,'yokedmethod') || cal.describe.yokedmethod == 0)
    OKTODO = 0;
end
if (~isfield(cal,'yoked') || ~isfield(cal.yoked,'spectra'))
    OKTODO = 0;
end
if (cal.nPrimaryBases == 0)
    OKTODO = 0;
end
if (~OKTODO)
    return;
end

%% Average yoked measurements for this primary
yokedSpds = cal.yoked.spectra;

%% Fit each spectrum with the linear model for all three primaries
% and project down onto this
projectedYokedSpd = cal.P_device*(cal.P_device\yokedSpds);

%% Now we have to adjust the linear model so that it has our standard 
% properties.

% Make first three basis functions fit maxSpd exactly
maxSpd = projectedYokedSpd(:,end);
weights = cal.P_device\maxSpd;
currentLinMod = zeros(size(cal.P_device));
for i = 1:cal.nDevices
    tempLinMod = 0;
    for j = 1:cal.nPrimaryBases
        tempLinMod = tempLinMod + cal.P_device(:,i+(j-1)*cal.nDevices)*weights(i+(j-1)*cal.nDevices);
    end
    currentLinMod(:,i) = tempLinMod;
end
weights = currentLinMod(:,1:cal.nDevices)\maxSpd;
for i = 1:cal.nDevices
    currentLinMod(:,i) = currentLinMod(:,i)*weights(i);
end
maxPow = max(max(currentLinMod(:,1:cal.nDevices)));

% Now find the rest of the linear model
clear tempLinMod
for i = 1:cal.nDevices
    for j = 1:cal.nPrimaryBases
        tempLinMod(:,j) = cal.P_device(:,i+(j-1)*cal.nDevices); %#ok<AGROW>
    end
    residual = tempLinMod - currentLinMod(:,i)*(currentLinMod(:,i)\tempLinMod);
    restOfLinMod = FindLinMod(residual,cal.nPrimaryBases-1);
    for j = 2:cal.nPrimaryBases
        tempMax = max(abs(restOfLinMod(:,j-1)));
        currentLinMod(:,i+(j-1)*cal.nDevices) = maxPow*restOfLinMod(:,j-1)/tempMax;
    end
end

% Span of cal.P_device and currentLinMod should be the same.  Check this.
if (DEBUG)
    check = currentLinMod - cal.P_device*(cal.P_device\currentLinMod);
    if (max(abs(check(:))) > 1e-10)
        error('Two linear models that should have the same span don''t');
    end
end

% Express yoked spectra in terms of model weights
gammaTable = currentLinMod\cal.yoked.spectra;
tempSpd = currentLinMod*gammaTable;
for i = 1:cal.nDevices
    index = gammaTable(i,:) > 1;
    gammaTable(i,index) = 1;
    gammaTable(i,:) = MakeMonotonic(HalfRect(gammaTable(i,:)'))';
end

% Stash info in calibration structure
cal.P_device = currentLinMod;

% When R=G=B, we just use the common settings.
if (cal.describe.yokedmethod == 1)
    cal.rawdata.rawGammaInput = cal.yoked.settings(1,:)';
	cal.rawdata.rawGammaTable = gammaTable';
    
% When measurements are at a specified chromaticity, need to interpolate gamma
% functions so that we have them for each device on a common scale.
elseif (cal.describe.yokedmethod == 2)
    cal.rawdata.rawGammaInput = cal.yoked.settings';
	cal.rawdata.rawGammaTable = gammaTable';
end

%% Debugging
if (DEBUG)
    S = [380 4 101];
    load T_xyz1931
    T_xyz=683*SplineCmf(S_xyz1931, T_xyz1931, S);
    
    % Meausured xyY
    measuredYokedxyY = XYZToxyY(T_xyz*cal.yoked.spectra);
    
    % Raw linear model fit xyY
    projectedYokedxyY = XYZToxyY(T_xyz*projectedYokedSpd);
    
    % Predicted xyY
    predictedSpd = cal.P_device*cal.rawdata.rawGammaTable';
    predictedYokedxyY = XYZToxyY(T_xyz*predictedSpd);
    
    % Plot luminance obtained vs. desired
    [lumPlot,f] = StartFigure('standard');
    f.xrange = [0 size(cal.yoked.settings, 2)]; f.nxticks = 6;
    f.yrange = [0 360]; f.nyticks = 5;
    f.xtickformat = '%0.0f'; f.ytickformat = '%0.2f ';
    plot(measuredYokedxyY(3,:)','ro','MarkerSize',f.basicmarkersize,'MarkerFaceColor','r');
    plot(projectedYokedxyY(3,:)','bo','MarkerSize',f.basicmarkersize,'MarkerFaceColor','b');
    plot(predictedYokedxyY(3,:)','go','MarkerSize',f.basicmarkersize,'MarkerFaceColor','g');
    
    xlabel('Test #','FontName',f.fontname,'FontSize',f.labelfontsize);
    ylabel('Luminance (cd/m2)','FontName',f.fontname,'FontSize',f.labelfontsize);
    FinishFigure(lumPlot,f);
    
    % Plot x chromaticity obtained vs. desired
    [xPlot,f] = StartFigure('standard');
    f.xrange = [0 size(cal.yoked.settings, 2)]; f.nxticks = 6;
    f.yrange = [0.2 0.6]; f.nyticks = 5;
    f.xtickformat = '%0.0f'; f.ytickformat = '%0.2f ';
    plot(measuredYokedxyY(1,:)','ro','MarkerSize',f.basicmarkersize,'MarkerFaceColor','r');
    plot(projectedYokedxyY(1,:)','bo','MarkerSize',f.basicmarkersize,'MarkerFaceColor','b');
    plot(predictedYokedxyY(1,:)','go','MarkerSize',f.basicmarkersize,'MarkerFaceColor','g');
    xlabel('Test #','FontName',f.fontname,'FontSize',f.labelfontsize);
    ylabel('x chromaticity','FontName',f.fontname,'FontSize',f.labelfontsize);
    FinishFigure(xPlot,f);
    
    % Plot y chromaticity obtained vs. desired
    [yPlot,f] = StartFigure('standard');
    f.xrange = [0 size(cal.yoked.settings, 2)]; f.nxticks = 6;
    f.yrange = [0.2 0.6]; f.nyticks = 5;
    f.xtickformat = '%0.0f'; f.ytickformat = '%0.2f ';
    plot(measuredYokedxyY(2,:)','ro','MarkerSize',f.basicmarkersize,'MarkerFaceColor','r');
    plot(projectedYokedxyY(2,:)','bo','MarkerSize',f.basicmarkersize,'MarkerFaceColor','b');
    plot(predictedYokedxyY(2,:)','go','MarkerSize',f.basicmarkersize,'MarkerFaceColor','g');
    xlabel('Test #','FontName',f.fontname,'FontSize',f.labelfontsize);
    ylabel('y chromaticity','FontName',f.fontname,'FontSize',f.labelfontsize);
    FinishFigure(yPlot,f);
    
    drawnow;
end


return


