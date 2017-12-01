function dE2000 = ComputeDE2000_Lab(Lab1,Lab2, kLCH)
% Computes the deltaE values (CIEDE2000) for color pairs given in CIELAB coordinates (L*a*b*).
%
%       dE2000 = ComputeDE2000_Lab(Lab1,Lab2[,kLCH]);
%
% The inputs 'Lab1' and 'Lab2' are N*3 matrices each with one color per row [L*,a*,b*]. 
% The output is a column vector with N elements, containing the according 
% dE2000 values.
% 'kLCH' is optional and contains one or three elements for the weighting 
% factors kL, kC, kH (either one value for all, or three values for each). 
% The weighting factors are set to 1, if 'kLCH' is not specified.
%
% This implementation might not be the most efficient but it follows the formulas
% given in https://en.wikipedia.org/wiki/Color_difference#CIEDE2000 as closely
% as possible for better traceability. 

% HISTORY:
% 12/01/2017 (MR = Marc Repnow, EPFL)  Initial version.

if nargin==0
    testMe();
    return
end
if any(size(Lab1)~=size(Lab2)) || size(Lab1,2)~=3
    error('L*a*b* input parameters must be of size Nx3 (N>=1).');
end
if nargin>2
    if length(kLCH)==1
        kLCH = kLCH*[1,1,1];
    elseif isempty(kLCH)
        kLCH = [1,1,1];
    elseif numel(kLCH)~=3
        error('Input parameter ''kLCH'' must have either 1 or 3 elements.');
    end
else
    kLCH = [1,1,1];    
end
%--- "Degree"-versions of sin() and cos().
sinD = @(deg) sin(deg*pi()/180);
cosD = @(deg) cos(deg*pi()/180);
%---
Lstar12 = [Lab1(:,1), Lab2(:,1)];
aStar12 = [Lab1(:,2), Lab2(:,2)];
bStar12 = [Lab1(:,3), Lab2(:,3)];
CprimeStar12 = sqrt(aStar12.^2  + bStar12.^2);
dLprime = Lstar12(:,2) - Lstar12(:,1);
Lbar = (Lstar12(:,2) + Lstar12(:,1))/2;
Cbar = (CprimeStar12(:,1)+CprimeStar12(:,2))/2;
aPrime12 = aStar12 .* repmat(1 + 0.5*(1 - sqrt(Cbar.^7./(Cbar.^7 + 25^7))), 1,2);  
Cprime12 = sqrt(aPrime12.^2 + bStar12.^2);
dCprime = Cprime12(:,2)-Cprime12(:,1);
CbarPrime = (Cprime12(:,1)+Cprime12(:,2))/2;
%--- Note that atan2(0,0) will return 0.
hPrime12 = mod(atan2(bStar12, aPrime12)*180/pi() + 10*360, 360);
N = length(dLprime);
dhPrime = zeros(N,1);
HbarPrime = zeros(N,1);
for i=1:N
    diff_hPrime12 = hPrime12(i,2)-hPrime12(i,1);
    sum_hPrime12 = hPrime12(i,2)+hPrime12(i,1);
    if any(abs(Cprime12(i,:))<10*eps)
        dhPrime(i) = 0;
    elseif abs(diff_hPrime12)<=180
        dhPrime(i) = diff_hPrime12;
    elseif diff_hPrime12 <= 0
        dhPrime(i) = diff_hPrime12 + 360;
    else
        dhPrime(i) = diff_hPrime12 - 360;
    end
    if any(abs(Cprime12(i,:))<10*eps)
        HbarPrime(i) = sum_hPrime12;
    elseif abs(diff_hPrime12)<=180
        HbarPrime(i) = (sum_hPrime12)/2;
    elseif sum_hPrime12 < 360
        HbarPrime(i) = (sum_hPrime12 + 360)/2;
    else
        HbarPrime(i) = (sum_hPrime12 - 360)/2;
    end
end
dHprime = 2*sqrt(Cprime12(:,1).*Cprime12(:,2)).*sinD(dhPrime/2);
T = 1 - 0.17*cosD(HbarPrime-30)  + 0.24*cosD(2*HbarPrime) ...
    + 0.32*cosD(3*HbarPrime+6) - 0.2*cosD(4*HbarPrime-63);
SL = 1 + 0.015*(Lbar-50).^2 ./ sqrt(20 + (Lbar-50).^2);
SC = 1 + 0.045*CbarPrime;
SH = 1 + 0.015*CbarPrime.*T;
RT = -2*sqrt(CbarPrime.^7./(CbarPrime.^7 + 25^7)) .* sinD(60*exp(-((HbarPrime-275)./25).^2));
kL = kLCH(1);
kC = kLCH(2);
kH = kLCH(3);
dE2000 = sqrt((dLprime./(kL*SL)).^2 + (dCprime./(kC*SC)).^2 + (dHprime./(kH*SH)).^2 ...
           + RT.*dCprime./(kC*SC).*dHprime./(kH*SH));       
end %%%%


%--------------------------------------------------------------------------------------
function testMe()
% Test data are from 
%       Sharma et al., 2005, "The CIEDE2000 Color-Difference Formula: Implementation Notes, 
%                             Supplementary Test Data, and Mathematical Observations". 
% The values in this publication are given with 4 decimal places. For this test function here, 
% the MATLAB script provided by Sharma was used to get more precise values (8 decimal places). 

%       L*a*b* 1,                   L*a*b* 2,                   dE2000
testset = ...
    [   50.0000, 2.6772, -79.7751,  50.0000, 0.0000, -82.7485,  2.04245968;
        50.0000, 3.1571, -77.2803,  50.0000, 0.0000, -82.7485,  2.86151017;
        50.0000, 2.8361, -74.0200,  50.0000, 0.0000, -82.7485,  3.44119060;
        50.0000, -1.3802, -84.2814, 50.0000, 0.0000, -82.7485,  0.99999886;
        50.0000, -1.1848, -84.8006, 50.0000, 0.0000, -82.7485,  1.00000470;
        50.0000, -0.9009, -85.5211, 50.0000, 0.0000, -82.7485,  1.00001297;
        50.0000, 0.0000, 0.0000,    50.0000, -1.0000, 2.0000,   2.36685882;
        50.0000, -1.0000, 2.0000,   50.0000, 0.0000, 0.0000,    2.36685882;
        50.0000, 2.4900, -0.0010,   50.0000, -2.4900, 0.0009,   7.17917201;
        50.0000, 2.4900, -0.0010,   50.0000, -2.4900, 0.0010,   7.17916264;
        50.0000, 2.4900, -0.0010,   50.0000, -2.4900, 0.0011,   7.21947215;
        50.0000, 2.4900, -0.0010,   50.0000, -2.4900, 0.0012,   7.21947421;
        50.0000, -0.0010, 2.4900,   50.0000, 0.0009, -2.4900,   4.80452169;
        50.0000, -0.0010, 2.4900,   50.0000, 0.0010, -2.4900,   4.80452451;
        50.0000, -0.0010, 2.4900,   50.0000, 0.0011, -2.4900,   4.74607111;
        50.0000, 2.5000, 0.0000,    50.0000, 0.0000, -2.5000,   4.30648210;
        50.0000, 2.5000, 0.0000,    73.0000, 25.0000, -18.0000, 27.14923130;
        50.0000, 2.5000, 0.0000,    61.0000, -5.0000, 29.0000,  22.89769247;
        50.0000, 2.5000, 0.0000,    56.0000, -27.0000, -3.0000, 31.90300465;
        50.0000, 2.5000, 0.0000,    58.0000, 24.0000, 15.0000,  19.45352143;
        50.0000, 2.5000, 0.0000,    50.0000, 3.1736, 0.5854,    1.00002634;
        50.0000, 2.5000, 0.0000,    50.0000, 3.2972, 0.0000,    0.99997287;
        50.0000, 2.5000, 0.0000,    50.0000, 1.8634, 0.5757,    1.00004950;
        50.0000, 2.5000, 0.0000,    50.0000, 3.2592, 0.3350,    1.00003476;
        60.2574, -34.0099, 36.2677, 60.4626, -34.1751, 39.4387, 1.26442001;
        63.0109, -31.0961, -5.8663, 62.8187, -29.7946, -4.0864, 1.26295930;
        61.2901, 3.7196, -5.3901,   61.4292, 2.2480, -4.9620,   1.87307050;
        35.0831, -44.1164, 3.7933,  35.0232, -40.0716, 1.5901,  1.86449523;
        22.7233, 20.0904, -46.6940, 23.0331, 14.9730, -42.5619, 2.03725827;
        36.4612, 47.8580, 18.3852,  36.2715, 50.5065, 21.2231,  1.41457792;
        90.8027, -2.0831, 1.4410,   91.1528, -1.6435, 0.0447,   1.44412908;
        90.9257, -0.5406, -0.9208,  88.6381, -0.8985, -0.7239,  1.53811701;
        6.7747, -0.2908, -2.4247,   5.8714, -0.0985, -2.2286,   0.63772767;
        2.0776, 0.0795, -1.1350,    0.9033, -0.0636, -0.5514,   0.90823284];
dE2000 = ComputeDE2000_Lab(testset(:,1:3), testset(:,4:6));
%--- Our test case deltaE's are given with 8 decimal places, 
%    so ignore differences that can be explained by the according
%    round-off errors.
if any(abs(dE2000 - testset(:,end)) > 0.5e-8)
    fprintf('ComputeDE2000_Lab test failed.\n');
    [~, i] = max(abs(dE2000 - testset(:,end)));
    fprintf('e.g. worst test case (#%d): Lab2dE2000([%.4f,%.4f,%.4f],[%.4f,%.4f,%.4f]==%.8f (expected: %.8f)\n',...
        i, testset(i,1:6), dE2000(i),  testset(i,end));
else
    fprintf('ComputeDE2000_Lab test passed.\n');
end        
end
