function y = srLens(nm)
% y = srLens(nm)
%
% Return Stockman-Rider lens density Fourier fit
%
% Adopted by Claude AI and DHB from Stockman-Rider paper and Python code.
%
% See also StockmanRiderDemo, srLMSconelog.

% History:
%   2025-09-05  dhb  Matlab first version as described above.
% 

    % Lens best-fitting 9x2 Fourier Series Polynomial
    x = (nm-360.0)/95.49296586; % Rescales 360 to 660 nm to 0 to pi
    
    % Compute Fourier polynomial, setting to zero outside the range where it is non-zeroß
    y = zeros(size(nm));    
    c = [-313.9508632762, -70.3216819666, 585.4719725809, 471.5395862431, 117.3539102044, ...
         127.0168222865, -324.4700544731, -188.1638078982, -104.5512488013, -68.3078486904, ...
         89.7815373733, 33.4498264952, 35.2723638870, 13.6524086627, -8.7568168893, ...
         -1.2825766708, -3.5126531075, -0.4477840959, 0.0428291365, 1.0091871745];
    
    for i = 1:length(x)
        if x(i) <= ((660-360)/95.49296586)
            y(i) = (c(1) + c(2)*cos(x(i)) + c(3)*sin(x(i)) + c(4)*cos(2*x(i)) + c(5)*sin(2*x(i)) + ...
                   c(6)*cos(3*x(i)) + c(7)*sin(3*x(i)) + c(8)*cos(4*x(i)) + c(9)*sin(4*x(i)) + ...
                   c(10)*cos(5*x(i)) + c(11)*sin(5*x(i)) + c(12)*cos(6*x(i)) + c(13)*sin(6*x(i)) + ...
                   c(14)*cos(7*x(i)) + c(15)*sin(7*x(i)) + c(16)*cos(8*x(i)) + c(17)*sin(8*x(i)) + ...
                   c(18)*cos(9*x(i)) + c(19)*sin(9*x(i))) * c(20);
        else
            y(i) = 0;
        end
    end
end