function transmission = ODToTransmission(density)
% Convert optical density to transmission
%
% Synopsis:
%   transmission = ODToTransmission(density)
%
% Description
%   As its name says.  Handy not to have to remember the formula.
%   Converts all entries of its argument.

% History
%   03/31/26  dhb  Wrote it.

transmission = 10.^(-density);