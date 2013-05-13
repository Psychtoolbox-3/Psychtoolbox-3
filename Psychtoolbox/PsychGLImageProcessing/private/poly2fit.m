function p = poly2fit(x,y,z,n)
% POLY2FIT	POLY2FIT(x,y,z,n) finds the coefficients of a two-dimensional
%		polynomial formed from the data in the matrices x and y
%		that fits the data in the matrix z in a least-squares sense.
%               The coefficients are returned in the vector p, in descending
%               powers of x and y. For example, the second order polynomial
%
%                  x ^ 2 + xy +2x -3 y - 1
% 
%               would be returned as p = [ 1 1 0 2 -3 -1 ]
%

%		Written my Jim Rees, 12 Aug. 1991

if any((size(x) ~= size(y)) | size(x) ~= size(z))
	error('X,Y, and Z matrices must be the same size')
end

% First straighten out x, y, and z. 
x = x(:); y = y(:); z = z(:);

n = n+1;
k = 1;
%A = zeros(x);
for i = n:-1:1,
  for j=1:i,
	A(:,k) =  ((x.^(i-j)).*(y.^(j-1))) ;
        k = k + 1;
  end
end

p = (A\z).';
