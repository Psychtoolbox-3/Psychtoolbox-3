function LMSout = srLMSconelog(nm, Lshift, Mshift, Sshift, loglin)
% LMSout = srLMSconelog(nm, Lshift, Mshift, Sshift, loglin)
%
% This returns the fits to the Stockman-Sharp LMS absorbances.  Paper
% Figure 1 and Table 1.
%
% Pass 'lin' as last argument to get linear absorbance, otherwise returns log10
% absorbance.
%
% Return is matrix with four columns, first is wavelength in nm, then L, M, and S
% absorbances (or log10 absorbances) in that order.
%
% Adopted by Claude AI and DHB from Stockman-Rider paper and Python code.  See
% StockmanRiderDemo for more info.

% History:
%   2025-09-05  dhb  Matlab first version as described above.

    % Calculate log10 LMS absorbances from individual Fourier templates
    LMSout = zeros(length(nm), 4);
    LMSout(:,1) = nm(:); % 
    LMSout(:,2) = srLconelog(nm, Lshift);
    LMSout(:,3) = srMconelog(nm, Mshift);
    LMSout(:,4) = srSconelog(nm, Sshift);
    
    if strcmp(loglin, 'lin')
        for n = 2:4
            LMSout(:,n) = 10.^(LMSout(:,n)); % Lin absorbances
        end
    end
end