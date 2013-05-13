function str = Ind2Str(p,numchar)
% str = Ind2Str(p)
% converts indices (numbers) to characters (base 10 to base 26 conversion):
% Ind2Str( 1) = 'a'
% Ind2Str(26) = 'z'
% Ind2Str(27) = 'ba'
% if P is a vector or matrix of indices, output will be a cell of the
% same dimensions
% representation with least number of characters necessary is used
%
% str = Ind2Str(p,numchar)
% representation of at least NUMCHAR length will be used for output:
% Ind2Str([1,34,45],3) = {'aab','abi','abt'}

% DN en JB 2008
% DN 28-05-2008 complete rewrite adding matrix and arbitrary length support

s       = size(p);
p       = p(:);

if nargin==2
    str     = owndec2base(p,26,numchar);
else
    str     = owndec2base(p,26);
end

if ~isscalar(p)
    str = num2cell(str,2);  % gooi m in een cell
    str = reshape(str,s);
end


function s = owndec2base(d,b,nin)
% gejat uit dec2base, eigen karakterset in variabele symbols
d = double(d);
b = double(b);
n = max(1,round(log2(max(d)+1)/log2(b)));
while any(b.^n <= d)
    n = n + 1;
end
if nargin == 3
    n = max(n,nin);
end
s(:,n) = rem(d,b);
% any(d) must come first as it short circuts for empties
while any(d) && n >1
    n = n - 1;
    d = floor(d/b);
    s(:,n) = rem(d,b);
end
symbols = 'abcdefghijklmnopqrstuvwxyz';
s = reshape(symbols(s + 1),size(s));
