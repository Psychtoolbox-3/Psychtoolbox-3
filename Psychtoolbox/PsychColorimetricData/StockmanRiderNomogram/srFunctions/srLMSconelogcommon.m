function y = srLMSconelogcommon(nm, LMS_type, shift)
    % Best-fitting log-lin 8x2 Fourier Series Polynomial
    
    Lsercommonlmax = 557.5;
    Lalacommonlmax = 554.8;
    Mcommonlmax = 527.3;
    Scommonlmax = 418.5;
    
    % Log shifts in lmax from Lser with log 390 to log 850 scaled 0 to pi
    Soffset = -1.048690123;
    Moffset = -0.2036522967;
    Lalaoffset = -0.01775262143;
    Lseroffset = 0;
    
    x = (log10(nm)-2.556302500767287267)/0.1187666467581842301;
    
    switch LMS_type
        case 'Lser'
            x = x + Lseroffset + log10(Lsercommonlmax/(Lsercommonlmax+shift))/0.1187666467581842301;
        case 'Lala'
            x = x + Lalaoffset + log10(Lalacommonlmax/(Lalacommonlmax+shift))/0.1187666467581842301;
        case 'M'
            x = x + Moffset + log10(Mcommonlmax/(Mcommonlmax+shift))/0.1187666467581842301;
        case 'S'
            x = x + Soffset + log10(Scommonlmax/(Scommonlmax+shift))/0.1187666467581842301;
        otherwise
            error('Cone type not specified');
    end
    
    c = [-2.1256563197, 5.4677929400, 0.8960658918, -0.9530108239, -5.0377095815, ...
         -3.0039987529, -0.9508620342, -1.3670849620, 1.7702113766, 0.5165048525, ...
         1.1505501831, 0.6100416117, 0.0518211044, 0.1009282570, -0.1773573074, ...
         -0.0278798136, -0.0427736834, 0.0007050030];
    
    y = c(1) + c(2)*cos(x) + c(3)*sin(x) + c(4)*cos(2*x) + c(5)*sin(2*x) + ...
        c(6)*cos(3*x) + c(7)*sin(3*x) + c(8)*cos(4*x) + c(9)*sin(4*x) + ...
        c(10)*cos(5*x) + c(11)*sin(5*x) + c(12)*cos(6*x) + c(13)*sin(6*x) + ...
        c(14)*cos(7*x) + c(15)*sin(7*x) + c(16)*cos(8*x) + c(17)*sin(8*x) + c(18);
end