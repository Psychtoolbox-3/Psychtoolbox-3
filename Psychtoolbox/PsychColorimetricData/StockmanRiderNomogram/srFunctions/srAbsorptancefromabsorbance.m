function LMSout = srAbsorptancefromabsorbance(LMSabsf, Lod, Mod, Sod, loglin)
    % Calculates LMS absorptances from absorbances given ODs
    LMSout = zeros(size(LMSabsf));
    LMSout(:,1) = LMSabsf(:,1); % Wavelength in column 1
    
    LMSout(:,2) = (1-10.^(-Lod*LMSabsf(:,2)))./(1-10.^(-Lod)); % Lin absorbances
    LMSout(:,3) = (1-10.^(-Mod*LMSabsf(:,3)))./(1-10.^(-Mod));
    LMSout(:,4) = (1-10.^(-Sod*LMSabsf(:,4)))./(1-10.^(-Sod));
    
    if strcmp(loglin, 'log')
        for n = 2:4
            LMSout(:,n) = log10(LMSabsf(:,n)); % Log absorbances
        end
    end
end