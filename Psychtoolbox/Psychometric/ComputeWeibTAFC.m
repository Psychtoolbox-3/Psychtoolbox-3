function pCorrect = ComputeWeibTAFC(x,alpha,beta)
% pCorrect = ComputeWeibTAFC(x,alpha,beta)
%
% Compute the value of a TAFC weibull psychometric function.
% If 'x' is a matrix and the parameters are scalars, the same 
% psychometric function is applied to the whole matrix.
% If 'x' is a matrix and the parameters are vectors, the
% parameters are taken as applying to separate columns.
%
%    pCorrect = ( 1.0 - 0.5*exp( - (x./alpha).^beta ) )
%
% 9/15/93  jms         Added some special casing to deal with matrix inputs
%                      with/without vectors of parameters
% 8/26/94  dhb, ccc	   Handle weird values of alpha, beta
% 10/13/00 dhb         Improve initial guess for alpha.  Thanks to Duje Tadin
%                      for identifying the need for this.

% Check alpha values.  We enforce alpha > verySmall.
verySmall = 1e-15;
index = find( alpha < verySmall );
if (~isempty(index))
	alpha(index) = verySmall*ones(size(index));
end

% Check beta values.  We enforce beta > verySmall.
verySmall = 1e-15;
index = find( beta < verySmall );
if (~isempty(index))
	beta(index) = verySmall*ones(size(index));
end

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
  pCorrect = ( 1.0 - 0.5*exp( - z ) );

% Case:  nAlpha == nBeta == nCol
elseif (nAlpha == nBeta && nAlpha == nCol)
  % disp('Vector parameter case');
  for (i=1:nCol)
    z = (x(:,i)./alpha(i)).^beta(i);
    pCorrect(:,i) = ( 1.0 - 0.5*exp( - z ) );
  end
else
  error('Number of input columns not equal to number of parameters');
end


