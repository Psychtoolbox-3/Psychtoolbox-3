function ROC = NormalROC(logbetas,us,vars,un,varn)
% ROC = NormalROC(logbetas,us,vars,un,varn)
%
% Compute ROC curve given the parameters.  Done
% analytically for speed.
%
% The trick here is to solve for the values of x
% such that log(l(x)) > logbetas.  When the variances
% are unequal, there can be two distinct regions.
% I basically handle this by brute force, since there
% are not too many distinct cases.
%
% This routine is not carefully tested and may
% contain errors.
%
% xx/xx/xx  dhb  Wrote it.
% 10/4/00   dhb  Added caveat about bug possibilities, based
%                in part on disagreement between this routine
%                and some Monte Carlo simulations.

% Set up room for answer
[mbetas,null] = size(logbetas);
if ( mbetas == 1 && null ~= 1)
  logbetas = logbetas';
  mbetas = null;
end
ROC = zeros(mbetas,2);

% By convention, will assume that us > un, switch if
% false
if ( un > us ) 
  ut = us;
  vart = vars;
  us = un;
  vars = varn;
  un = ut;
  varn = vart;
end

% Set up the quadratic equation to find the zeros of
% the likelyhood ratio.  This expression is arrived at
% by writing out the likelihood ratio for Normal random
% variables, taking logarithms, and getting a quadratic
% expression for the values of x that generate the desired
% log beta.
a = (vars-varn)*ones(mbetas,1);
b = 2*(varn*us - vars*un)*ones(mbetas,1);
c = (vars*(un^2)-varn*(us^2))*ones(mbetas,1) + ...
			 varn*vars*log(varn/vars)*ones(mbetas,1) - ...
				2*varn*vars*logbetas;

if (vars == varn)
  x = (-1*c) ./ b;
  ROC(:,1) = ones(mbetas,1) - NormalCumulative(x,us,vars); 
  ROC(:,2) = ones(mbetas,1) - NormalCumulative(x,un,varn);
else
  % Solve for the zeros of the quatdratic
  discrim = sqrt( b.^2 - 4.*a.*c );
  x1 = (-b + discrim) ./ (2*a); 
  x2 = (-b - discrim) ./ (2*a);

  % Handle the case where ratio heads to +infinity
  if ( (vars-varn) > 0 )
    % Imaginary roots mean 1,1
    index = find( imag(discrim) ~= 0 );
    [m,null] = size(index);
    if ( m ~= 0 )
      ROC(index,1) = ones(m,1);
      ROC(index,2) = ones(m,1);
    end

    index = find( imag(discrim) == 0 );
    [m,null] = size(index);
    if ( m ~= 0 )
      ROC(index,1) = ones(m,1) - ...
                     NormalCumulative(x1(index),us,vars) - ...
                     NormalCumulative(x2(index),us,vars);
      ROC(index,2) = ones(m,1) - ...
                     NormalCumulative(x1(index),un,varn) - ...
                     NormalCumulative(x2(index),un,varn);
    end

  % Handle the case where the ratio heads to -infinity
  else
    % Imaginary roots mean 0,0
    index = find( imag(discrim) ~= 0 );
    [m,null] = size(index);
    if ( m ~= 0 )
      ROC(index,1) = zeros(m,1);
      ROC(index,2) = zeros(m,1);
    end

    index = find( imag(discrim) == 0 );
    [m,null] = size(index);
    if ( m ~= 0 )
      ROC(index,1) = NormalCumulative(x2(index),us,vars) - ...
                     NormalCumulative(x1(index),us,vars);
      ROC(index,2) = NormalCumulative(x2(index),un,varn) - ...
                     NormalCumulative(x1(index),un,varn);
    end
  end
end



