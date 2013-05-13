% AnsiZ136MPEBasicTest
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Test code for our implementation of ANSI Z136.1-2007. Reproduces many figures from the
% standard.
%
% 2/22/13  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Figure 9b: Test T2 computation
%
% Answer should range between 10 and 100
% as size increases over the specified range.
% See Figure 9b, p. 102.
fprintf('Reproducing Figure 9b, p. 102\n');
theStimulusSizesDeg = linspace(MradToDeg(0),MradToDeg(100+10),100);
theStimulusSizesMrad = DegToMrad(theStimulusSizesDeg);
for i = 1:length(theStimulusSizesDeg)
    T2Sec(i) = AnsiZ136MPEComputeT2(theStimulusSizesDeg(i));
end
figure; clf; hold on
plot(theStimulusSizesMrad,T2Sec,'ro','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Stimulus Size (mrad)');
ylabel('T2 (sec)');
xlim([0 max(theStimulusSizesMrad)]);
ylim([0 100]);
title('Figure 9b: Test of AnsiZ136MPEComputeT2');
grid on

%% Figure 8a: Test Ca computation
%
% Answer should increase between 1 and 5
% with wavelength between 700 and 1050,
% and flatten out on the two sides.
%
% See Figure 8a, p. 98.
fprintf('Reproducing Figure 8a, p. 98\n');
wavelengthsNm = 400:1399;
for i = 1:length(wavelengthsNm)
    Ca(i) = AnsiZ136MPEComputeCa(wavelengthsNm(i));
end
figure; clf; hold on
semilogy(wavelengthsNm,log10(Ca),'ro','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Wavelength (nm)');
ylabel('Log10 Ca');
xlim([min(wavelengthsNm) max(wavelengthsNm)]);
ylim([0 1]);
title('Figure 8a: Test of AnsiZ136MPEComputeCa');
grid on

%% Figure 8c: Test Cb computation
%
% Answer should range between 10 and 100
% as size increases over the specified range.  This
% should look like Figure 8c, p. 100.
fprintf('Reproducing Figure 8c, p. 100\n');
wavelengthsNm = 380:780;
for i = 1:length(wavelengthsNm)
    Cb(i) = AnsiZ136MPEComputeCb(wavelengthsNm(i));
end
figure; clf; hold on
semilogy(wavelengthsNm,log10(Cb),'ro','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Wavelength (nm)');
ylabel('Log 10 Cb');
xlim([min(wavelengthsNm) max(wavelengthsNm)]);
ylim([0 3]);
title('Figure 8c: Test of AnsiZ136MPEComputeCb');
grid on

%% Figure 8b: Test Cc computation
%
% Answer should range between 1 and 8
% with wavelength between 1150 and 1200 nm.
% should look like Figure 8b, p. 99.
fprintf('Reproducing Figure 8b, p. 99\n');
wavelengthsNm = 1050:1399;
for i = 1:length(wavelengthsNm)
    Cc(i) = AnsiZ136MPEComputeCc(wavelengthsNm(i));
end
figure; clf; hold on
semilogy(wavelengthsNm,log10(Cc),'ro','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Wavelength (nm)');
ylabel('Log10 Cc');
xlim([min(wavelengthsNm) max(wavelengthsNm)]);
ylim([0 1]);
title('Figure 8b: Test of AnsiZ136MPEComputeCc');
grid on

%% Figure 3: Test limiting cone angle computation
%
% Answer should range between 11 and 110
% with duration between 100 and 1e4 seconds.
% This should look like Figure 3, p. 93.
fprintf('Reproducing Figure 3, p. 93\n');
durations = logspace(1,4.2);
for i = 1:length(durations)
   limitingConeAngles(i) = AnsiZ136MPEComputeLimitingConeAngle(durations(i));
end
figure; clf; hold on
loglog(log10(durations),log10(limitingConeAngles),'ro','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Log10 Stimulus Duration (sec)');
ylabel('Log10 Limiting Cone Angle (mrad)');
xlim([1 5]);
ylim([0 3]);
title('Figure 3: Test of AnsiZ136MPEComputeLimitingConeAngle');
grid on

%% Figure 7: Test photochemical and thermal limits for extended sources.
%
% This code reproduces Figure 7, p. 97.  Figure 7 is for wavelengths between
% 400 and 700.  The overall limit (but not the photochemical limit) is 
% independent of wavelength over this time interval.
%
% We only compute/plot down to 10-8 seconds, because our code doesn't
% implement the limts for extremely short times.
%
% Our plot also shows the photochemical limit (in red) down to the time
% where that is relevant.  Since it is above the overall limit, it
% would not affect that limit in the regime plotted in this figure.
fprintf('Reproducing Figure 7, p. 97\n');

% Specify what parameters to test
theStimulusWavelengthsNm = 400:20:700;
theStimulusSizesMrad = [1 7.5 25 100];
minLogDuration = -13;
maxLogDuration = 0;
stimulusDurationsSec = logspace(minLogDuration,maxLogDuration,1000);

radiantExposureFig7 = figure; clf; set(gcf,'Position',[770 670 1000 600]);
for s = 1:length(theStimulusSizesMrad)
    fprintf('\tSize %0.1f mRad\n',theStimulusSizesMrad(s));
    stimulusSizeMrad = theStimulusSizesMrad(s);
    stimulusSizeDeg = MradToDeg(stimulusSizeMrad);
    for w = 1:length(theStimulusWavelengthsNm)
        stimulusWavelengthNm = theStimulusWavelengthsNm(w);
        
        for t = 1:length(stimulusDurationsSec)
            stimulusDurationSec = stimulusDurationsSec(t);
            [~, ~, ~, MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2(w,t)] = ...
                AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
            
            [~, ~, ~, MPELimitCornealRadiantExposure_JoulesPerCm2(w,t)] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
        end
    end
    
    % Does the answer depend on wavelength?  Yes for photochemical limit
    % but no for overall limit.  You can explore if you want by enabling this section of code.
    if (0)
        minMPELimitCornealRadiantExposure_JoulesPerCm2 = min(MPELimitCornealRadiantExposure_JoulesPerCm2,[],1);
        minMPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = min(MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2,[],1);
        for w = 1:length(theStimulusWavelengthsNm)
            if (any(MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2(w,:) ~= minMPEPhotochemicalCornealRadiantExposure_JoulesPerCm2))
                fprintf('\t\tWavelength dependence for photochemical limit, wavelength %d\n',theStimulusWavelengthsNm(w));
            end
            if (any(MPELimitCornealRadiantExposure_JoulesPerCm2(w,:) ~= minMPELimitCornealRadiantExposure_JoulesPerCm2))
                fprintf('\t\tWavelength dependence for overall limit, wavelength %d\n',theStimulusWavelengthsNm(w));
            end
        end
    end
    
    figure(radiantExposureFig7); % subplot(1,length(theStimulusSizesMrad),s);
    hold on
    loglog(log10(stimulusDurationsSec),log10(min(MPELimitCornealRadiantExposure_JoulesPerCm2,[],1)),'bo','MarkerSize',8,'MarkerFaceColor','b');
    loglog(log10(stimulusDurationsSec),log10(min(MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2,[],1)),'ro','MarkerSize',5,'MarkerFaceColor','r');
    drawnow;
end
xlabel('Log10 Stimulus Duration (sec)');
ylabel('Log10 Corneal Radiant Exposure (J/cm2)');
xlim([minLogDuration maxLogDuration]);
ylim([-8 0]);
title({'Figure 7: Test of AnsiZ136MPE Exposure Limits' ; 'Blue: Limit, Red: Photochemical Limit' ; sprintf('Size %0.1f mrad',stimulusSizeMrad) ; 'Wavelengths 400-700 nm'});
grid on

%% Figure 10: Test photochemical and thermal limits for extended sources.
%
% This code reproduces Figure 10, pp. 103-107.  Each version is for a
% different stimulus size, and shows the dependence of the limit on
% duration for different wavelengths.
%
% The agreement between what's produced here and the graphs in the
% standard is good for sizes <= 11 mrad, but diverges for larger
% sizes in terms of the photochemical limit.  The figures in the
% standard have a temporal break that depends on stimulus size
% for the photochemical limit, and there is no such dependence
% in the main formula in the table.

% Specify what parameters to test
theStimulusSizesMrad = [1 3 11 25 50];
theFigureNames = {'Figure 10a' 'Figure 10b' 'Figure 10c' 'Figure 10d' 'Figure 10e'};

for s = 1:length(theStimulusSizesMrad)
    fprintf('Reproducing %s\n',theFigureNames{s});
    stimulusSizeMrad = theStimulusSizesMrad(s);
    stimulusSizeDeg = MradToDeg(stimulusSizeMrad);
    switch (stimulusSizeMrad)
        case 1
            minLogDuration = -1; maxLogDuration = 4.2;
            minLogY = -4.1; maxLogY = -1;
            theStimulusWavelengthsNm = [400 450 475 490 700 1050 1200];
        case 3
            minLogDuration = -1; maxLogDuration = 4.2;
            minLogY = -4.1; maxLogY = 0;
            theStimulusWavelengthsNm = [400 450 475 500 700 1050 1200];
        case 11
            minLogDuration = -1; maxLogDuration = 4.2;
            minLogY = -4.1; maxLogY = 0;
            theStimulusWavelengthsNm = [400 450 475 500 514.5 700 1050 1200];
        case 25
            minLogDuration = -1; maxLogDuration = 4.2;
            minLogY = -4.1; maxLogY = 1;
            theStimulusWavelengthsNm = [400 450 475 500 532 700 1050 1200];
        case 50
            minLogDuration = -1; maxLogDuration = 4.2;
            minLogY = -4.1; maxLogY = 1;
            theStimulusWavelengthsNm = [400 450 475 500 532 550 700 1050 1200];
        otherwise
            error('Unexpected stimulus size specified');
    end
    stimulusDurationsSec = logspace(minLogDuration,maxLogDuration,1000);
    radiantExposureFig10 = figure; clf; set(gcf,'Position',[770 670 1000 600]);

    for w = 1:length(theStimulusWavelengthsNm)
        stimulusWavelengthNm = theStimulusWavelengthsNm(w);
        
        for t = 1:length(stimulusDurationsSec)
            stimulusDurationSec = stimulusDurationsSec(t);
            [~, ~, MPEPhotochemicalCornealIrradiance_WattsPerCm2(w,t), ~] = ...
                AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
            
            [~, ~, MPELimitCornealIrradiance_WattsPerCm2(w,t), ~] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
        end
        
        figure(radiantExposureFig10);
        hold on
        loglog(log10(stimulusDurationsSec),log10(MPELimitCornealIrradiance_WattsPerCm2(w,:)),'bo','MarkerSize',8,'MarkerFaceColor','b');
        index = find(MPEPhotochemicalCornealIrradiance_WattsPerCm2(w,:) ==  MPELimitCornealIrradiance_WattsPerCm2(w,:));
        loglog(log10(stimulusDurationsSec(index)),log10(MPEPhotochemicalCornealIrradiance_WattsPerCm2(w,(index))),'ro','MarkerSize',5,'MarkerFaceColor','r');
        drawnow;
    end
    
    xlabel('Log10 Stimulus Duration (sec)');
    ylabel('Log10 Corneal Irradiance (W/cm2)');
    xlim([minLogDuration maxLogDuration]);
    ylim([minLogY maxLogY]);
    title({'Test of AnsiZ136MPE Exposure Limits' ; 'Blue: Limit, Red Dashed: Photochemical Limit' ; sprintf('Size %0.1f mrad. %0.1f deg',theStimulusSizesMrad(s),stimulusSizeDeg) ; sprintf('Ansi Z136%s',theFigureNames{s})});
    grid on
end


%% Figure 11: Test photochemical and thermal limits for extended sources.
%
% This code reproduces Figure 11, pp. 108.  Each is for a a
% different stimulus size, and shows the dependence of the limit on
% duration for different wavelengths.
%
% This is close to the Figure 12 in the standard, although there
% are slight differences visible by eye, where the limits produced
% here are a bit lower than drawn in the document.
%
% The Mod1 version of the figure is for a smaller size, and is
% for comparison to Figure 10e.

% Specify what parameters to test
theStimulusSizesMrad = [110 50];
theFigureNames = {'Figure 11' 'Figure11Mod1'};

for s = 1:length(theStimulusSizesMrad)
    fprintf('Reproducing %s\n',theFigureNames{s});
    stimulusSizeMrad = theStimulusSizesMrad(s);
    stimulusSizeDeg = MradToDeg(stimulusSizeMrad);
    switch (stimulusSizeMrad)
        case 50
            minLogDuration = -1; maxLogDuration = 4.2;
            minLogY = -2; maxLogY = 3;
            theStimulusWavelengthsNm = [400 450 475 500 532 550 700 1050 1200];
        case 110
            minLogDuration = -1; maxLogDuration = 4.2;
            minLogY = -2; maxLogY = 3;
            theStimulusWavelengthsNm = [400 450 475 500 532 550 700 1050 1200];
        otherwise
            error('Unexpected stimulus size specified');
    end
    stimulusDurationsSec = logspace(minLogDuration,maxLogDuration,1000);
    radiantExposureFig11 = figure; clf; set(gcf,'Position',[770 670 1000 600]);

    for w = 1:length(theStimulusWavelengthsNm)
        stimulusWavelengthNm = theStimulusWavelengthsNm(w);
        
        for t = 1:length(stimulusDurationsSec)
            stimulusDurationSec = stimulusDurationsSec(t);
            [~, MPEPhotochemicalRadiance_WattsPerCm2Sr(w,t), ~, ~] = ...
                AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
            
            [~, MPELimitRadiance_WattsPerCm2Sr(w,t), ~, ~] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
        end
        
        figure(radiantExposureFig11);
        hold on
        loglog(log10(stimulusDurationsSec),log10(MPELimitRadiance_WattsPerCm2Sr(w,:)),'bo','MarkerSize',8,'MarkerFaceColor','b');
        index = find(abs(MPEPhotochemicalRadiance_WattsPerCm2Sr(w,:) - MPELimitRadiance_WattsPerCm2Sr(w,:)) < 1e-6);
        loglog(log10(stimulusDurationsSec(index)),log10(MPEPhotochemicalRadiance_WattsPerCm2Sr(w,(index))),'ro','MarkerSize',5,'MarkerFaceColor','r');
        drawnow;
    end
    
    xlabel('Log10 Stimulus Duration (sec)');
    ylabel('Log10 Radiance (W/[cm2-sr])');
    xlim([minLogDuration maxLogDuration]);
    ylim([minLogY maxLogY]);
    title({'Test of AnsiZ136MPE Exposure Limits' ; 'Blue: Limit, Red Dashed: Photochemical Limit' ; sprintf('Size %0.1f mrad. %0.1f deg',theStimulusSizesMrad(s),stimulusSizeDeg) ; sprintf('Ansi Z136%s',theFigureNames{s})});
    grid on
end


%% Figure 12: Test photochemical and thermal limits for extended sources.
%
% This code reproduces Figure 12, pp. 109.  Each is for a a
% different stimulus size, and shows the dependence of the limit on
% duration for different wavelengths.
%
% This is close to the Figure 12 in the standard, although there
% are slight differences visible by eye, where the limits produced
% here are a bit lower than drawn in the document.

% Specify what parameters to test
theStimulusSizesMrad = [110];
theFigureNames = {'Figure 12'};

for s = 1:length(theStimulusSizesMrad)
    fprintf('Reproducing %s\n',theFigureNames{s});
    stimulusSizeMrad = theStimulusSizesMrad(s);
    stimulusSizeDeg = MradToDeg(stimulusSizeMrad);
    switch (stimulusSizeMrad)
        case 110
            minLogDuration = -13; maxLogDuration = 0;
            minLogY = -4.1; maxLogY = 3;
            theStimulusWavelengthsNm = [400 700 1050 1200];
        otherwise
            error('Unexpected stimulus size specified');
    end
    stimulusDurationsSec = logspace(minLogDuration,maxLogDuration,1000);
    radiantExposureFig12 = figure; clf; set(gcf,'Position',[770 670 1000 600]);

    for w = 1:length(theStimulusWavelengthsNm)
        stimulusWavelengthNm = theStimulusWavelengthsNm(w);
        
        for t = 1:length(stimulusDurationsSec)
            stimulusDurationSec = stimulusDurationsSec(t);
            [MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr(w,t), ~, ~, ~] = ...
                AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
            
            [MPELimitIntegratedRadiance_JoulesPerCm2Sr(w,t), ~, ~, ~] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
        end
        
        figure(radiantExposureFig12);
        hold on
        loglog(log10(stimulusDurationsSec),log10(MPELimitIntegratedRadiance_JoulesPerCm2Sr(w,:)),'bo','MarkerSize',8,'MarkerFaceColor','b');
        index = find(abs(MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr(w,:) - MPELimitIntegratedRadiance_JoulesPerCm2Sr(w,:)) < 1e-8);
        loglog(log10(stimulusDurationsSec(index)),log10(MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr(w,(index))),'ro','MarkerSize',5,'MarkerFaceColor','r');
        drawnow;
    end
    
    xlabel('Log10 Stimulus Duration (sec)');
    ylabel('Log10 Integrated Radiance (J/[cm2-sr])');
    xlim([minLogDuration maxLogDuration]);
    ylim([minLogY maxLogY]);
    title({'Test of AnsiZ136MPE Exposure Limits' ; 'Blue: Limit, Red Dashed: Photochemical Limit' ; sprintf('Size %0.1f mrad. %0.1f deg',theStimulusSizesMrad(s),stimulusSizeDeg) ; sprintf('Ansi Z136%s',theFigureNames{s})});
    grid on
end


%% Make a plot of how limit varies with stimulus size, for specified
% duration and wavelength.  Take minimum over vectors specified for
% each.

% Specify what parameters to test
minLogSize = -1; maxLogSize = 2;
minLogYRad = -3; maxLogYRad = 2;
minLogYIrrad = -5; maxLogYIrrad = 0;
minLogYIntRad = 0; maxLogYIntRad = 3;
minLogYRadExp = -4; maxLogYRadExp = -1;
stimulusSizesDeg = logspace(minLogSize,maxLogSize,100);
stimulusWavelengthsNm = 400:20:1390;
stimulusDurationsSec = logspace(-1,4,100);
fprintf('Computing over stimulus sizes from %0.1f to %0.1f deg\n',min(stimulusSizesDeg),max(stimulusSizesDeg));
clear MPELimitIntegratedRadiance_JoulesPerCm2Sr MPELimitRadiance_WattsPerCm2Sr MPELimitCornealIrradiance_WattsPerCm2 MPELimitCornealRadiantExposure_JoulesPerCm2
for s = 1:length(stimulusSizesDeg)
    stimulusSizeDeg = theStimulusSizesDeg(s);
    stimulusSizeMrad = DegToMrad(stimulusSizeDeg);
    MPELimitIntegratedRadiance_JoulesPerCm2Sr(s) = Inf;
    MPELimitRadiance_WattsPerCm2Sr(s) = Inf;
    MPELimitCornealIrradiance_WattsPerCm2(s) = Inf;
    MPELimitCornealRadiantExposure_JoulesPerCm2(s) = Inf;
    for w = 1:length(stimulusWavelengthsNm)
        stimulusWavelengthNm = stimulusWavelengthsNm(w);  
        for t = 1:length(stimulusDurationsSec)
            stimulusDurationSec = stimulusDurationsSec(t);
             [temp1, temp2, temp3, temp4] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm);
            if (temp1 < MPELimitIntegratedRadiance_JoulesPerCm2Sr(s))
                MPELimitIntegratedRadiance_JoulesPerCm2Sr(s) = temp1;
            end 
            if (temp2 < MPELimitRadiance_WattsPerCm2Sr(s))
                MPELimitRadiance_WattsPerCm2Sr(s) = temp2;
            end 
            if (temp3 < MPELimitCornealIrradiance_WattsPerCm2(s))
                MPELimitCornealIrradiance_WattsPerCm2(s) = temp3;
            end
             if (temp4 < MPELimitCornealRadiantExposure_JoulesPerCm2(s))
                MPELimitCornealRadiantExposure_JoulesPerCm2(s) = temp4;
            end
        end  
    end   
end

stimulusSizeFig = figure; clf; set(gcf,'Position',[770 670 1000 1000]);
figure(stimulusSizeFig);
subplot(2,2,1); hold on
loglog(log10(stimulusSizesDeg),log10(MPELimitRadiance_WattsPerCm2Sr),'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('Log10 Stimulus Size (deg)');
ylabel('Log10 Radiance (W/[cm2-sr])');
xlim([minLogSize maxLogSize]);
ylim([minLogYRad maxLogYRad]);
title({'Test of AnsiZ136MPE Exposure Limits' ; sprintf('Durations %g to %g sec',min(stimulusDurationsSec),max(stimulusDurationsSec)) ; ...
    sprintf('Wavelengths %d to %d nm',min(stimulusWavelengthsNm),max(stimulusWavelengthsNm))});
grid on

subplot(2,2,2); hold on
loglog(log10(stimulusSizesDeg),log10(MPELimitCornealIrradiance_WattsPerCm2),'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('Log10 Stimulus Size (deg)');
ylabel('Log10 Corneal Irradiance (W/cm2)');
xlim([minLogSize maxLogSize]);
ylim([minLogYIrrad maxLogYIrrad]);
title({'Test of AnsiZ136MPE Exposure Limits' ; sprintf('Durations %g to %g sec',min(stimulusDurationsSec),max(stimulusDurationsSec)) ; ...
    sprintf('Wavelengths %d to %d nm',min(stimulusWavelengthsNm),max(stimulusWavelengthsNm))});
grid on

subplot(2,2,3); hold on
loglog(log10(stimulusSizesDeg),log10(MPELimitIntegratedRadiance_JoulesPerCm2Sr),'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('Log10 Stimulus Size (deg)');
ylabel('Log10 Integrated Radiance (J/[cm2-sr])');
xlim([minLogSize maxLogSize]);
ylim([minLogYIntRad maxLogYIntRad]);
title({'Test of AnsiZ136MPE Exposure Limits' ; sprintf('Durations %g to %g sec',min(stimulusDurationsSec),max(stimulusDurationsSec)) ; ...
    sprintf('Wavelengths %d to %d nm',min(stimulusWavelengthsNm),max(stimulusWavelengthsNm))});
grid on

subplot(2,2,4); hold on
loglog(log10(stimulusSizesDeg),log10(MPELimitCornealRadiantExposure_JoulesPerCm2),'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('Log10 Stimulus Size (deg)');
ylabel('Log10 Corneal Radiant Exposure (J/cm2)');
xlim([minLogSize maxLogSize]);
ylim([minLogYRadExp maxLogYRadExp]);
title({'Test of AnsiZ136MPE Exposure Limits' ; sprintf('Durations %g to %g sec',min(stimulusDurationsSec),max(stimulusDurationsSec)) ; ...
    sprintf('Wavelengths %d to %d nm',min(stimulusWavelengthsNm),max(stimulusWavelengthsNm))});
grid on
