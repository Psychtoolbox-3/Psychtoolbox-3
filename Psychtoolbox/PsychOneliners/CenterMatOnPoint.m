function [xinds yinds] = CenterMatOnPoint(mat,x,y)
% [xinds yinds] = CenterMatOnPoint(mat,x,y)
% 
% returns indices to center matrix on a point
% 
% if no point is provided, relative indices are computed
%
% DN 2008

if nargin ==1
    x = 0;
    y = 0;
end

[ys,xs] = AltSize(mat,[1 2]);

xinds   = floor(-xs/2+1:xs/2)+round(x);
yinds   = floor(-ys/2+1:ys/2)+round(y);
