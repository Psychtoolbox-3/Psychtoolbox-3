function energy = QuantaToEnergy(wls,quanta)
% energy = QuantaToEnergy(wls,quanta)
%
% Convert quantal units (quanta per unit wavelength)
% to energy units (energy or power per unit wavelength).
%
% Constants are set up so that we have energy in joules or
% power in watts.
%
% The routine is set up to convert spectra.  These are
% passed as the columns of the matrix quanta.  The
% wavelengths corresponding to each row are passed in
% the column vector wls.
%
% 7/29/96  dhb  Added comment.
% 8/16/96  dhb, abp  Modified interface.

wls = MakeItWls(wls);
h = 6.626e-34;
c = 2.998e8;
[n,m] = size(quanta);
energy = (quanta*h*c) ./ ((1e-9) * wls(:,ones(1,m)));
