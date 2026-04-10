function PhotopigmentNomogramDemo
% PhotopigmentNomogramDemo
%
% Compares photopigment nomograms against the CIE/Stockman-Sharpe tabulated
% photopigment absorbance spectra (T_log10coneabsorbance_ss), which serve as
% the fixed reference baseline.  Each nomogram gets its own figure: thick black
% lines show the CIE tabulated absorbance; the nomogram is overlaid in color.
% Rows are L, M, S cones; left column is linear absorbance, right column is
% log absorbance.
%
% Lambda max values and nomogram comparison notes:
%
% Each nomogram uses its own canonical lambda max values. Because the
% nomograms have different shapes, they require different lambda max
% parameters to approximate the same underlying spectrum.
%
% For the L cone, the CIE absorbance peaks at ~551.9 nm (per the Fourier
% polynomial fit of Stockman & Rider 2023, Table 1). The StockmanRider nomogram
% reproduces this accurately at lambda max = 551.9 (zero shift). The
% StockmanSharpe nomogram, which uses a polynomial approximation, requires
% lambda max = 558.9 to best match the same CIE L-cone absorbance — a ~7 nm
% discrepancy. This does not mean the two nomograms represent different spectra;
% both target the same CIE L-cone, but StockmanSharpe is less accurate on a
% linear scale (see note in StockmanSharpeNomogram.m). The visual difference
% between the two curves reflects that accuracy difference.
%
% For the S cone the discrepancy between nomograms is larger because
% StockmanSharpe uses a single polynomial shape for all three cone types (just
% shifted in log wavelength), while StockmanRider uses separately fitted Fourier
% polynomials for L, M, and S. The CIE S-cone absorbance has a different shape
% from L and M on a log-wavelength scale, so StockmanSharpe fits it less well.
%
% The CarrollNeitz nomogram has a different history and was not developed
% to match the CIE cone fundamentals. It provides a worse approximation to
% to these, at least for the lambda max values we clocked in, which come
% from a Netiz and Neitz review as referenced in the CarrollNeitzNomogram
% function.
%
% History:
%   2026-04-10  dhb  Wrote it.
%   2026-04-10  dhb  Restructured: CIE tabulated data as fixed black reference,
%                    one figure per nomogram.

close all;

%% Load CIE tabulated photopigment absorbance data
% Rows are [L; M; S] log10 absorbance.
% S_log10coneabsorbance_ss is the PTB S-format wavelength spec.
cieData     = load('T_log10coneabsorbance_ss');
wls_cie     = SToWls(cieData.S_log10coneabsorbance_ss);
T_cie_log10 = cieData.T_log10coneabsorbance_ss;  % log10 absorbance
T_cie_lin   = 10.^T_cie_log10;                   % linear absorbance

%% Nomogram sampling
% 1 nm spacing triggers peak normalization inside each nomogram function.
S   = [380 1 321];
wls = SToWls(S);

coneNames = {'L', 'M', 'S'};

% Each nomogram entry: {name, lambdaMax [L;M;S], color, linewidth}
nomograms = { ...
    'StockmanRider',  [551.9; 529.8; 416.9], [0   0.5 0.2],  3.0 ; ...
    'StockmanSharpe', [558.9; 530.3; 420.7], 'b',             2.0 ; ...
    'CarrollNeitz',   [557.5; 530.0; 420.0], 'r',             1.5 ; ...
};
nNomograms = size(nomograms, 1);

%% Compute nomogram absorbances
T = cell(nNomograms, 1);
for nn = 1:nNomograms
    T{nn} = PhotopigmentNomogram(S, nomograms{nn,2}, nomograms{nn,1});
end

%% One figure per nomogram, CIE tabulated data as thick black reference
for nn = 1:nNomograms
    figure; clf;
    set(gcf, 'Position', [100 100 1200 900]);
    sgtitle(sprintf('%s  (\\lambda_{max}: L=%.1f, M=%.1f, S=%.1f nm)', ...
        nomograms{nn,1}, nomograms{nn,2}(1), nomograms{nn,2}(2), nomograms{nn,2}(3)), ...
        'FontSize', 13, 'FontWeight', 'bold');

    for cc = 1:3
        % --- Linear absorbance ---
        subplot(3, 2, (cc-1)*2 + 1); hold on;
        plot(wls_cie, T_cie_lin(cc,:), 'k', 'LineWidth', 4);
        plot(wls, T{nn}(cc,:), 'Color', nomograms{nn,3}, 'LineWidth', nomograms{nn,4});
        xlabel('Wavelength (nm)', 'FontSize', 11);
        ylabel('Absorbance', 'FontSize', 11);
        title(sprintf('%s cone — Linear', coneNames{cc}), 'FontSize', 12);
        xlim([380 700]);
        ylim([0 1.05]);

        % --- Log absorbance ---
        subplot(3, 2, (cc-1)*2 + 2); hold on;
        plot(wls_cie, T_cie_log10(cc,:), 'k', 'LineWidth', 4);
        plot(wls, log10(T{nn}(cc,:)), 'Color', nomograms{nn,3}, 'LineWidth', nomograms{nn,4});
        xlabel('Wavelength (nm)', 'FontSize', 11);
        ylabel('Log_{10} Absorbance', 'FontSize', 11);
        title(sprintf('%s cone — Log', coneNames{cc}), 'FontSize', 12);
        xlim([380 700]);
        ylim([-4 0.1]);
        legend({'CIE/SS tabulated', nomograms{nn,1}}, ...
            'Location', 'NorthEastOutside', 'FontSize', 9);
    end
end

end
