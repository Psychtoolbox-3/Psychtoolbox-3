function LMSout = srEnergyfromquantalin(LMSin, loglin)
    LMSout = LMSin; % For local calculation and return
    
    for n = 2:4
        LMSout(:,n) = LMSin(:,n) .* LMSin(:,1);
        LMSout(:,n) = LMSout(:,n) ./ max(LMSout(:,n)); % renormalize
    end
    
    if strcmp(loglin, 'log')
        for n = 2:4
            LMSout(:,n) = log10(LMSout(:,n)); % Log absorptances
        end
    end
end