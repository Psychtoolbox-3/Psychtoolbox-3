function PhotopigmentNomogramDemo
% PhotopigmentNomogramDemo
%
% Exercises the suite of photopigment nomograms available via
% PhotopigmentNomogram.  Uses the nominal L, M, and S cone lambda max
% values from the Stockman-Rider nomogram (551.9, 529.8, 416.9 nm).
%
% Six-panel figure: rows are L, M, S; left column is linear absorbance,
% right column is log absorbance.
%
% History:
%   2026-04-10  dhb  Wrote it.

close all;

% Wavelength sampling: 1 nm spacing triggers peak normalization in all
% nomogram functions.
S   = [380 1 321];
wls = SToWls(S);

% Nominal StockmanRider lambda max values (L, M, S)
lambdaMaxLMS = [551.9; 529.8; 416.9];
coneNames    = {'L', 'M', 'S'};

% Compute StockmanRider absorbance (one row per cone)
T_SR = PhotopigmentNomogram(S, lambdaMaxLMS, 'StockmanRider');

% Plot: 3 rows x 2 columns
figure; clf;
set(gcf, 'Position', [100 100 900 750]);

for cc = 1:3
    % --- Linear absorbance ---
    subplot(3, 2, (cc-1)*2 + 1); hold on;
    plot(wls, T_SR(cc,:), 'k', 'LineWidth', 2);
    xlabel('Wavelength (nm)', 'FontSize', 11);
    ylabel('Absorbance', 'FontSize', 11);
    title(sprintf('%s cone — Linear  (\\lambda_{max} = %g nm)', ...
        coneNames{cc}, lambdaMaxLMS(cc)), 'FontSize', 12);
    xlim([380 700]);
    ylim([0 1.05]);
    legend('StockmanRider', 'Location', 'NorthEast', 'FontSize', 9);

    % --- Log absorbance ---
    subplot(3, 2, (cc-1)*2 + 2); hold on;
    plot(wls, log10(T_SR(cc,:)), 'k', 'LineWidth', 2);
    xlabel('Wavelength (nm)', 'FontSize', 11);
    ylabel('Log_{10} Absorbance', 'FontSize', 11);
    title(sprintf('%s cone — Log  (\\lambda_{max} = %g nm)', ...
        coneNames{cc}, lambdaMaxLMS(cc)), 'FontSize', 12);
    xlim([380 700]);
    ylim([-4 0.1]);
    legend('StockmanRider', 'Location', 'NorthEast', 'FontSize', 9);
end

end
