function density = TransmissionToOD(transmission)
% Convert transmission to optical density
%
% Synopsis:
%   density = TransmissionToOD(transmission)
%
% Description
%   As its name says.  Handy not to have to remember the formula.
%   Converts all entries of its argument.

% History
%   03/31/26  dhb  Wrote it.

% OD to transmission is: transmission = 10.^(-density);
% 
% So:
density = -log10(transmission);