function LMSout = srAbsorptancefromlinqcornea(LMSin, mac, lens, mac_460, lens_400, loglin)
    % Remove standard macular and lens
    LMSout = LMSin; % For local calculation and return
    
    mac_template_460 = 0.35; % Template peak
    lens_template_400 = 1.7649; % Template peak
    macscale = mac_460/mac_template_460;
    lensscale = lens_400/lens_template_400;
    
    for n = 2:4
        LMSout(:,n) = LMSin(:,n) .* 10.^(macscale.*mac) .* 10.^(lensscale.*10.^lens);
        LMSout(:,n) = LMSout(:,n) ./ max(LMSout(:,n)); % renormalise
    end
    
    if strcmp(loglin, 'log')
        for n = 2:4
            LMSout(:,n) = log10(LMSout(:,n)); % Log absorptances
        end
    end
end