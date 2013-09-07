function [axialOpticalDensity] = ComputeAxialDensity(specificDensity, outSegmentLength)
% [axialOpticalDensity] = ComputeAxialDensity(specificDensity, outSegmentLength)
%
% Compute the axial optical density using the estimate of specific density and
% the length of out segment length.  Specific density is sometimes called 
% concentration (e.g. Wyszecki and Stiles, p. 588).
%
% 06/11/03 lyin Wrote it.
% 06/26/03 dhb	Change "peakAxialOpticDensity" to "axialOpticalDensity".

axialOpticalDensity = outSegmentLength .* specificDensity;
