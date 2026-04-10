function RGBCMFs = srCalculateRGBCMFs(LMS_energy, primaries, nm_step, Lshift, Mshift, Lod, Mod, Sod, mac_460, lens_400)
    % Calculate RGB color matching functions from LMS functions
    
    CMF_Rnm = primaries(1);
    CMF_Gnm = primaries(2);
    CMF_Bnm = primaries(3);
    
    % Calculate LMS values at primary wavelengths
    nm_primaries = [CMF_Rnm; CMF_Gnm; CMF_Bnm];
    
    % Generate templates for primary wavelengths
    mac_primaries = srMacular(nm_primaries);
    lens_primaries = srLens(nm_primaries);
    coneabs_primaries = srLMSconelog(nm_primaries, Lshift, Mshift, 0, 'lin');
    
    conenewq_retina_primaries = srAbsorptancefromabsorbance(coneabs_primaries, Lod, Mod, Sod, 'lin');
    conenewq_cornea_primaries = srCorneafromlinabsorptance(conenewq_retina_primaries, mac_primaries, lens_primaries, mac_460, lens_400, 'lin');
    conenewe_cornea_primaries = srEnergyfromquantalin(conenewq_cornea_primaries, 'lin');
    
    % RGB matrix (columns are LMS responses to R, G, B primaries)
    RGBLMS = conenewe_cornea_primaries(:, 2:4)'; % Transpose to get 3x3 matrix
    
    % LMS to RGB transformation matrix
    LMSRGB = inv(RGBLMS);
    
    % Transform LMS to RGB
    RGBCMFs = LMS_energy; % Copy structure
    RGBCMFs(:, 2:4) = LMS_energy(:, 2:4) * LMSRGB';
end