function LMSout = srAbsorbancefromabsorptance(LMSin, Lod, Mod, Sod, loglin)
    LMSout = LMSin; % For local calculation and return
    
    LMSout(:,2) = -log10(1-LMSin(:,2).*(1-10.^(-Lod)))./Lod;
    LMSout(:,3) = -log10(1-LMSin(:,3).*(1-10.^(-Mod)))./Mod;
    LMSout(:,4) = -log10(1-LMSin(:,4).*(1-10.^(-Sod)))./Sod;
    
    if strcmp(loglin, 'log')
        for n = 2:4
            LMSout(:,n) = log10(LMSout(:,n)); % Log absorptances
        end
    end
end
