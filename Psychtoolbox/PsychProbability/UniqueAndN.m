function [b,n] = UniqueAndN(input,operateByRow)
% [uniques,n] = UniqueAndN(input,[operateByRow=false])
%
% returns the unique elements in input (sorted) and
% the number of times the item ocurred
% optionally does so for unique rows (operateByRow=true)
% instead of unique elements
% 
% DN    2008
% DN    2017 Added byRow option

if nargin<2 || isempty(operateByRow) || ~operateByRow
    
    input   = sort(input);  % Necessary for the trick below to work
    
    [b,~,j] = unique(input);
    
else
    input   = sortrows(input);
    
    [b,~,j] = unique(input,'rows');
end

d       = diff([0; j(:); max(j)+1]);

inds    = find(d);
n       = inds(2:end)-inds(1:end-1);
