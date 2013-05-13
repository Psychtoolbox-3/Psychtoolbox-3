function [b,n] = UniqueAndN(in)
% [uniques,n] = uniqueen(input)
%
% returns the unique elements in input (sorted) and
% the number of times the item ocurred
% 
% DN    2008

in      = sort(in);  % Necessary for the trick below to work

[b,i,j] = unique(in);

d       = diff([0; j(:); max(j)+1]);

inds    = find(d);
n       = inds(2:end)-inds(1:end-1);
