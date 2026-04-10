function LMSout = srCorneafromlinabsorptance(LMSin, mac, lens, mac_460, lens_400, loglin)
    LMSout = LMSin; % For local calculation and return
    
    mac_template_460 = 0.35; % Template peak
    lens_template_400 = 1.7649; % Template peak
    macscale = mac_460/mac_template_460;
    lensscale = lens_400/lens_template_400;
    
    for n = 2:4
        LMSout(:,n) = LMSin(:,n) ./ (10.^(mac.*macscale) .* 10.^(lens.*lensscale));
        LMSout(:,n) = LMSout(:,n) ./ max(LMSout(:,n)); % renormalize
    end
    
    if strcmp(loglin, 'log')
        for n = 2:4
            LMSout(:,n) = log10(LMSout(:,n)); % Log absorptances
        end
    end
end