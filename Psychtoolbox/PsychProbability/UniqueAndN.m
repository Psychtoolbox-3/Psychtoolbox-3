function [b,n] = UniqueAndN(in,qByRow)
% [uniques,n] = uniqueen(input,qByRow)
%
% returns the unique elements in input (sorted) and
% the number of times the item ocurred
% optionally does so for unique rows (qByRow=true) instead of unique
% elements
% 
% DN    2008
% DN    2017 Added byRow option

if nargin<2 || ~qByRow
    
    in      = sort(in);  % Necessary for the trick below to work
    
    [b,~,j] = unique(in);
    
else
    in      = sortrows(in);
    
    [b,~,j] = unique(in,'rows');
end

d       = diff([0; j(:); max(j)+1]);

inds    = find(d);
n       = inds(2:end)-inds(1:end-1);