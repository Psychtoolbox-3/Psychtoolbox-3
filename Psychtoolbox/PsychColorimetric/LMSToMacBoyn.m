function ls = LMSToMacBoyn(LMS,T_cones,T_lum)
% ls = LMSToMacBoyn(LMS,[T_cones,T_lum])
%
% Compute MacLeod-Boynton chromaticity from cone exciation coordinates.
% This is L/Lum and S/Lum, with appropriate normalization as described
% below. Recommended usage yields MB chromacities according to CIE 170-2:2015.
%
% ** Recommended Usage: Compute LMS with respect to some specified T_cones,
% and pass both T_cones as well as the cooresponding T_lum (the photopic
% luminosity function.) T_lum should be a linear combination of the L and M
% cone fundamentals.
% 
% This routine will then scale passed L and M values so that they sum to
% the best linear approximation of luminance and then normalize the L
% excitation by luminance (as computed by the linear combination) to obtain
% the l chromaticity. It will normalize the S excitaiton by luminance to
% obtain s chromaticity, with an overall scaling so that the maximum value
% of this chromaticity is 1 taken over the visible spectrum.
%
% Note that the s cone scaling can vary a bit depending on the wavelength
% sampling of the passed T_cones and T_lum, since the max is taken over
% these. If you use the T_cones_ss2/T_cones_ss10 and T_CIE_Y2/T_CIE_Y10
% files provided in PTB, the default sampling is at 1 nm and this is fine.
% If you use subsampled wavelength spacing, the computation of the s cone
% scaling will begin to deviate from the standard.  But so will your
% computation of LMS values, so this isn't really an issue specific to this
% routine.
%
% When you use the CIE cone fundamentals and corresponding luminance
% functions, this procedure yields the MacLeod-Boynton chromaticity
% diagrams as specified in CIE 170-2:2015.
%
% ** Legacy Usage: Just pass LMS values. In this case, we assume that the
% passed LMS values were computed with respect to the Smith-Pokorny
% fundamentals normalized to a peak of 1 and Judd-Vos luminance (more or
% less).  That is, this usage assumes LMS was computed using the
% fundamentals stored in PTB's T_cones_sp. This is old usage and preserved
% for backwards compatibility, but the three argument usage as described
% above is preferred for clarity. Moreover, in this case, the s
% chromaticity is not further normalized.  This leads to S chromaticities
% considerably larger than those obtained with the new usage.
%
% ** A Backwards Incompatibility. The scaling for s chromaticity to match
% CIE 170-2:2015 was introduced in Janurary 2019 and is not backwards
% compatible with previous behavior when T_cones and T_lum are passed.
% Preserving such compatibility did not seem important enough relative to
% the gains of having this work as now specified in the CIE standard.
%
% 10/30/97  dhb  Wrote it.
% 7/9/02    dhb  T_cones_sp -> T_cones on line 20.  Thanks to Eiji Kimura.
% 1/23/19   dhb  Scale s chromaticity value to be consistent with CIE
%                170-2:2015, when T_cones and T_lum are passed.  This is
%                not backwards compatible with previous scaling, but it
%                seems good to match the standard. Thanks to Danny Garside
%                for pointing out the scaling specified in the 2015
%                standard.

% Examples:
%{
    % Recreate the spectrum locus and equal energy white shown in Figure 8.2
    % of CIE 170-2:2015. Also performs a regression check.
    load T_cones_ss2
    load T_CIE_Y2
    lsSpectrumLocus = LMSToMacBoyn(T_cones_ss2,T_cones_ss2,T_CIE_Y2);
 
    % Compute the sum of the ls values in the spectrum locus, and compare
    % to the value that this example computed in February 2019, entered
    % here to four places as 412.2608.  This comparison provides a
    % check that this routine still works the way it did when we put in the
    % check.
    check = round(sum(lsSpectrumLocus(:)),4);
    if (abs(check-412.2608) > 1e-4)
        error('No longer get same check value as we used to');
    end

    % Compute ls for equal energy white
    LMSEEWhite = sum(T_cones_ss2,2);
    lsEEWhite = LMSToMacBoyn(LMSEEWhite,T_cones_ss2,T_CIE_Y2); 

    % Plot
    figure; hold on;
    plot(lsSpectrumLocus(1,:)',lsSpectrumLocus(2,:)','r','LineWidth',3);
    plot(lsEEWhite(1),lsEEWhite(2),'bs','MarkerFaceColor','b','MarkerSize',12);
    xlim([0.4 1]); ylim([0,1]);
    xlabel('l chromaticity'); ylabel('s chromaticity');
    title('CIE 170-2:2015 Figure 8.2');
%}
%{
    % Demonstrate invariance of ls after scaling of cones and luminance, as
    % long as LMS valued are computed with respect to passed cones and
    % luminance.
    load T_cones_ss2
    load T_CIE_Y2

    % Spectrum locus ls chromaticity with no scaling.
    lsSpectrumLocus = LMSToMacBoyn(T_cones_ss2,T_cones_ss2,T_CIE_Y2); 


    % Do some scaling and recompute spectrum locus.
    % The choices of 1.8, 3, 0.05 as scaling are
    % just three numbers I made up. You can use any three numbers and it
    % should still work, modulo any numerical issues with very big or
    % very small numbers.
    T_CIE_Y2_scaled = 1.8*T_CIE_Y2;
    T_cones_ss2_scaled = T_cones_ss2;
    T_cones_ss2_scaled(1,:) = 3*T_cones_ss2(1,:);
    T_cones_ss2_scaled(3,:) = 0.05*T_cones_ss2(3,:);
    lsSpectrumLocusScaled = LMSToMacBoyn(T_cones_ss2_scaled,T_cones_ss2_scaled,T_CIE_Y2_scaled); 

    % Make sure the difference is very small numerically.
    diff = max(abs(lsSpectrumLocus(:)-lsSpectrumLocusScaled(:)));
    fprintf('Difference in spectrum locus chromaticity after scaling is %0.5f (should be small)\n',diff);

    % Plot the locus computed both ways to compare visually.
    figure; hold on;
    plot(lsSpectrumLocus(1,:)',lsSpectrumLocus(2,:)','r','LineWidth',3);
    plot(lsSpectrumLocusScaled(1,:)',lsSpectrumLocusScaled(2,:)','g','LineWidth',2);
    xlim([0.4 1]); ylim([0,1]);
    xlabel('l chromaticity'); ylabel('s chromaticity');
%}

% Scale LMS so that L+M = luminance and S cone value corresponds to a
% fundamental with a max of 1.
if (nargin == 1)
	LMS = diag([0.6373 0.3924 1]')*LMS;
elseif (nargin == 3 )
	factorsLM = (T_cones(1:2,:)'\T_lum');
    factorS = 1/max(T_cones(3,:)./(factorsLM(1)*T_cones(1,:) + factorsLM(2)*T_cones(2,:)));
	LMS = diag([factorsLM ; factorS])*LMS;
else
    error('Number of input arguments should be either 1 or 3');
end

% Compute ls coordinates from LMS
n = size(LMS,2);
ls = zeros(2,n);
denom = [1 1 0]*LMS;
ls = LMS([1 3],:) ./ ([1 1]'*denom);
