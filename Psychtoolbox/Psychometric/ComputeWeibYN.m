function pyes = ComputeWeibYN(x,alpha,beta)
% pyes = ComputeWeibYN(x,alpha,beta)
%
% Compute the value of a yes-no weibull psychometric function.
% If 'x' is a matrix and the parameters are scalars, the same 
% psychometric function is applied to the whole matrix.
% If 'x' is a matrix and the parameters are vectors, the
% parameters are taken as applying to separate columns.
%
%    pyes = ( 1.0 - exp( - (x./alpha).^beta ) )
%
% 9/15/93  jms  Added some special casing to deal with matrix inputs
%               with/without vectors of parameters

% Check dimensions of parameters
[nRow,nAlpha] = size(alpha);
if (nRow ~= 1)
  error('alpha must be a row vector');
end

[nRow,nBeta]  = size(beta);
if (nRow ~= 1)
  error('beta must be a row vector');
end

[nInputs, nCol] = size(x);

% Case:  nAlpha == nBeta == 1 ~= nCol
if (nAlpha == nBeta && nAlpha == 1)
  % disp('Single parameter case');
  z = (x./alpha).^beta;
  pyes = ( 1.0 - exp( - z ) );

% Case:  nAlpha == nBeta == nCol
elseif (nAlpha == nBeta && nAlpha == nCol)
  % disp('Vector parameter case');
  for (i=1:nCol)
    z = (x(:,i)./alpha(i)).^beta(i);
    pyes(:,i) = ( 1.0 - exp( - z ) );
  end
else
  error('Number of input columns not equal to number of parameters');
end


