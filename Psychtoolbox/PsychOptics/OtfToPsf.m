function [xGridMinutes,yGridMinutes,psf] = OtfToPsf(xSfGridCyclesDeg,ySfGridCyclesDeg,otf,varargin)
%OTFTOPSF  Convert a 2D optical transfer fucntion to a 2D point spread function.
%    [xSfGridCyclesDeg,ySfGridCyclesDeg,otf] = PsfToOtf([xGridMinutes,yGridMinutes],psf)
%
%    Converts a optical transfer function specified over two-dimensional
%    spatial frequency in cycles per degree to a point spread function
%    specified over positions in minutes.  For human vision, these are each
%    natural units.
%
%    The input spatial frequencies should be specified in matlab's grid
%    matrix format and sf for x and y should be specified over the same
%    spatialfrequency extent and with the same number of evenly spaced
%    samples. Spatial frequency (0,0) should be at location floor(n/2)+1 in
%    each dimension.  The PSF is returned with position (0,0) at location
%    floor(n/2)+1 in each dimension.  This is the form we mostly want to
%    look at and use.  
%
%    The OTF is assummed to have the DC term in the center poistion,
%    floor(n/2)+1 of the passed matrix.  If you have your hands on an OTF
%    with the DC term in the upper left position (1,1), apply fftshift to
%    it before passing to this routine.  The DC in upper left is the Matlab
%    native format for applying the ifft, and is also the format stored by
%    isetbio in its optics structure.
%
%    Positions are returned using the same conventions.
%
%    No normalization is performed.  The psf should be real, and we
%    complain (throw an error) if it is not, to reasonable numerial
%    precision. If it seems OK, we make sure it is real.
%
%    We also make sure that the returned psf is all postive and sums to 
%    1.  In some cases, we found that there were small negative values
%    and after setting these to zero renormalization was needed.
%
%    We wrote this rather than simply relying on Matlab otf2psf/psf2otf
%    because we want a routine that deals with the conversion of spatial
%    frequency to spatial support.
%
%    If you pass the both sf args as empty, both position grids are
%    returned as empty and just the conversion on the OTF is performed.
%
%    PsychOpticsTest shows that this works very well when we go back and
%    forth for diffraction limited OTF/PSF.  But not exactly exactly
%    perfectly.  A signal processing maven might be able to track down
%    whether this is just a numerical thing or whether some is some small
%    error, for example in how position is converted to sf or back again in
%    the PsfToOtf.
%
% Optional key/value pairs:
%    'warningInsteadOfErrorForNegativeValuedPSF'  - Set to 1 (default
%                                                   0) to get a warning
%                                                   not an error if the psf
%                                                   values are too negative
%                                                   before they are forced
%                                                   to zero. Set to 2 for
%                                                   no warning msg or
%                                                   error.
%    'negativeFractionalTolerance'                - The error/warning is
%                                                   thrown if the magnitude
%                                                   of the most negative
%                                                   value is more than this
%                                                   fraction of the maximum
%                                                   (positve) value.
%                                                   (Default 1e-3).
%
% See also: PsfToOtf, PsychOpticsTest.

% History:
%   01/26/18  dhb 
%   03/31/18  dhb   Document key/value pair added by someone else.
%             dhb   Add key/value pair for negative value tolerance.
%                   This is now 1e-3 rather than 1e-10
%   10/17/18  dhb   Now 10-2.

%% Parse input
p = inputParser;
p.addParameter('warningInsteadOfErrorForNegativeValuedPSF', 0, @isnumeric);
p.addParameter('negativeFractionalTolerance', 1e-2, @isnumeric);
p.parse(varargin{:});

%% Handle sf args and converstion
if (~isempty(xSfGridCyclesDeg) & ~isempty(ySfGridCyclesDeg))
    % They can both be passed as non-empty, in which case we do a set of sanity
    % checks and then do the conversion.
    [m,n] = size(xSfGridCyclesDeg);
    centerPosition = floor(n/2) + 1;
    if (m ~= n)
        error('psf must be passed on a square array');
    end
    [m1,n1] = size(ySfGridCyclesDeg);
    if (m1 ~= m || n1 ~= n)
        error('x and y positions are not consistent');
    end
    [m2,n2] = size(otf);
    if (m2 ~= m || n2 ~= n)
        error('x and y positions are not consistent');
    end
    if (~all(xSfGridCyclesDeg(:,centerPosition) == 0))
        error('Zero spatial frequency is not in right place in the passed xGrid');
    end
    if (~all(ySfGridCyclesDeg(centerPosition,:) == 0))
        error('Zero spatial frequency is not in right place in the passed yGrid');
    end
    if (xSfGridCyclesDeg(1,centerPosition) ~= ySfGridCyclesDeg(centerPosition,1))
        error('Spatial frequency extent of x and y grids does not match');
    end
    diffX = diff(xSfGridCyclesDeg(:,centerPosition));
    if (any(diffX ~= diffX(1)))
        error('X positions not evenly spaced');
    end
    diffY = diff(ySfGridCyclesDeg(centerPosition,:));
    if (any(diffY ~= diffY(1)))
        error('Y positions not evenly spaced');
    end
    if (diffX(1) ~= diffY(1))
        error('Spatial frequency sampling in x and y not matched');e
    end
    
    %% Generate position grids
    %
    % Samples are evenly spaced and the same for both x and y (checked above).
    % Handle even versus odd dimension properly for fft conventions.
    [xGridMinutes,yGridMinutes] = SfGridCyclesDegToPositionGridMinutes(xSfGridCyclesDeg,ySfGridCyclesDeg);
    
elseif (isempty(xSfGridCyclesDeg) & isempty(ySfGridCyclesDeg))
    % This case is OK, we set the output grids to empty
    xGridMinutes = [];
    yGridMinutes = [];
    
else
    % This case is not allowable
    error('Either both sf grids must be empty, or neither');
end

%% Compute otf and put spatial center in the middle of the grid
psf = fftshift(ifft2(ifftshift(otf)));

%% See if there is stray imaginary stuff, get rid of it if so.
%
% Throw an error if the returned psf isn't in essence real valued.
% Then set residual imaginary parts to 0.
if (any(abs(imag(psf(:))) > 1e-10))
    error('Computed psf is not sufficiently real');
end
if (any(imag(psf(:))) ~= 0)
    psf = psf - imag(psf)*1i;
end

% Check for large negative psf values, and then set any small
% negative values to zero.  For some cases (e.g. Marimont-Wandell
% OTF, we do get negative values because the way that was constructed
% is an approximation to measurements that does not absolutely guarantee an
% all positive OTF.
if (max(psf(:)) <= 0)
    error('Computed PSF has no positive values.  This is not good.');
end
if (min(psf(:)) < 0 && abs(min(psf(:))) > p.Results.negativeFractionalTolerance*max(psf(:)))
    if (p.Results.warningInsteadOfErrorForNegativeValuedPSF == 1)
        fprintf(2,'Mysteriously large negative psf values, min value is %g, relative to max of %g, fraction %g\n',min(psf(:)),max(psf(:)),abs(min(psf(:)))/max(psf(:)));
    elseif (p.Results.warningInsteadOfErrorForNegativeValuedPSF == 0)
        fprintf(2,'Mysteriously large negative psf values: min value is %g, relative to max of %g, fraction %g\n',min(psf(:)),max(psf(:)),abs(min(psf(:)))/max(psf(:)));
        error('Mysteriously large negative psf values');
    end
end
psf(psf < 0) = 0;

% Make sure return is real
psf = abs(psf);

% Make sure return sums to 1.  It might not be because of the
% above fussing with imaginary and negative values.
psf = psf/sum(psf(:));

end