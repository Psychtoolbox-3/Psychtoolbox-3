function y = srMacular(nm)
% y = srMacular(nm)
%
% Return Stockman-Rider macular density Fourier fit
%
% Adopted by Claude AI and DHB from Stockman-Rider paper and Python code.
%
% See also StockmanRiderDemo, srLMSconelog.

% History:
%   2025-09-05  dhb  Matlab first version as described above.
%

% Setup thetaP for macular, 375 to 550 nm 0 to pi
x = (nm-375)/55.70423008; 

% Macular best-fitting 11x2 Fourier Series Polynomial
y = zeros(size(nm));
c = [3712.2037792986, 374.1811575175, -7007.6989637831, -5887.2857515364, -633.0475233043, ...
    -716.0429039473, 4386.8811254914, 2882.1092658881, 638.1347550701, 468.4980700497, ...
    -1653.7567388120, -817.1240899995, -286.4038978705, -144.7996457395, 340.3364828167, ...
    115.5652804221, 59.1650826447, 18.6678197694, -30.2344535413, -5.4683753172, ...
    -4.1335064207, -0.5043959566, 0.5094171266, 1.0050048550];

for i = 1:length(x)
    if x(i) >= 0 && x(i) <= ((550-375)/55.70423008)
        y(i) = (c(1) + c(2)*cos(x(i)) + c(3)*sin(x(i)) + c(4)*cos(2*x(i)) + c(5)*sin(2*x(i)) + ...
            c(6)*cos(3*x(i)) + c(7)*sin(3*x(i)) + c(8)*cos(4*x(i)) + c(9)*sin(4*x(i)) + ...
            c(10)*cos(5*x(i)) + c(11)*sin(5*x(i)) + c(12)*cos(6*x(i)) + c(13)*sin(6*x(i)) + ...
            c(14)*cos(7*x(i)) + c(15)*sin(7*x(i)) + c(16)*cos(8*x(i)) + c(17)*sin(8*x(i)) + ...
            c(18)*cos(9*x(i)) + c(19)*sin(9*x(i)) + c(20)*cos(10*x(i)) + c(21)*sin(10*x(i)) + ...
            c(22)*cos(11*x(i)) + c(23)*sin(11*x(i))) * c(24);
    else
        y(i) = 0;
    end
end
end
