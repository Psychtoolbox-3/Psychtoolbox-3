function mtf = NavarroMTF(s, varargin)
%NAVARROMTF  Compute the MTF measured by Navarror, Artal and Williams
%   mtf = NAVARRORMTF(s)
% 
%   Compute the MTF measured by Navarro et al. The MTF was measured for a
%   4 mm pupil, and for a varietry of eccentricities. We return the foveal
%   MTF by default.  The desired eccentricity can be passed as a key-value
%   pair, i.e., mtf10deg = NavarroMTF(1:60, 'eccentricity', '10 deg')
%
%   Navarro, R., Artal, P., & Williams, D. R. (1993). 
%   Modulation transfer function of the human eye as a function of retinal 
%   eccentricity. Journal of the Optical Society of America A, 10, 201?212.
%
%   Spatial frequency passed in cycles/deg.
%
%   See also WILLIAMSMTF, OTFTOPSF, WILLIAMSRESTMTF, DIFFRACTIONMTF,
%   WILLIAMSTABULATEDPSF, PSYCHOPTICSTEST.

% 10/20/2918   npc		Wrote it.

%% Parse input
p = inputParser;
p.addParameter('eccentricity','foveal',@ischar);
p.parse(varargin{:});
eccentricity = p.Results.eccentricity;

switch (eccentricity)
    case 'foveal'
        A = 0.172; B = 0.037; C = 0.22;
    case '5 deg'
        A = 0.245; B = 0.041; C = 0.2;
    case '10 deg'
        A = 0.245; B = 0.041; C = 0.2;
    case '20 deg'
        A = 0.328; B = 0.038; C = 0.14;
    case '30 deg'
        A = 0.606; B = 0.064; C = 0.12;
    case '40 deg'
        A = 0.82; B = 0.064; C = 0.09;
    case '50 deg'
        A = 0.93; B = 0.059; C = 0.067;
    case '60 deg'
        A = 1.89; B = 0.108; C = 0.05;
    otherwise
        error('No data for ''%s''.', eccentricity);
end

mtf = (1-C)*exp(-A*s) + C*exp(-B*s);
