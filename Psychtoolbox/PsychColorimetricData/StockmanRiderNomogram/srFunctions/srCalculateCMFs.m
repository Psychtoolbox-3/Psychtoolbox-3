function [LMS_energy, LMS_quantal, RGBCMFs] = srCalculateCMFs(nm_step, Lshift, Mshift, Lod, Mod, Sod, mac_460, lens_400)
    % Main function to calculate color matching functions
    
    % Default parameters (Stockman & Sharpe 2-degree standard)
    if nargin < 1, nm_step = 1.0; end
    if nargin < 2, Lshift = 0.0; end
    if nargin < 3, Mshift = 0.0; end
    if nargin < 4, Lod = 0.50; end
    if nargin < 5, Mod = 0.50; end
    if nargin < 6, Sod = 0.40; end
    if nargin < 7, mac_460 = 0.350; end
    if nargin < 8, lens_400 = 1.7649; end
    
    % Set up wavelength array
    nm = (360:nm_step:850)';
    
    % Generate macular and lens templates
    mac = srMacular(nm);
    lens_template = srLens(nm);
    
    % Calculate cone absorbance templates
    coneabs_template = srLMSconelog(nm, Lshift, Mshift, 0, 'lin');
    
    % Retinal absorptances
    conenewq_retina = srAbsorptancefromabsorbance(coneabs_template, Lod, Mod, Sod, 'lin');
    
    % Corneal quantal spectral sensitivities
    conenewq_cornea = srCorneafromlinabsorptance(conenewq_retina, mac, lens_template, mac_460, lens_400, 'lin');
    
    % Corneal energy spectral sensitivities
    conenewe_cornea = srEnergyfromquantalin(conenewq_cornea, 'lin');
    
    LMS_energy = conenewe_cornea;
    LMS_quantal = conenewq_cornea;
    
    % Calculate RGB CMFs (Stiles & Burch primaries as default)
    RGBCMFs = srCalculateRGBCMFs(LMS_energy, [645.15, 526.32, 444.44], nm_step, Lshift, Mshift, Lod, Mod, Sod, mac_460, lens_400);
end