function inds = FindRepeatAlongDims(in,dim)
% find repeated rows or columns in a matrix.
% inds = FindRepeatAlongDims(in,dim)
% example:
% in =
%      1     2     3
%      2     3     4
%      4     5     6
%      1     2     3
%      1     2     3
%      4     5     6
% FindRepeatAlongDims(in,1)  % check for repeated rows
% ans = 
%      5                    % the fifth row is the same as the fourth row,
%                             repeat!
%
% for N-D dimensions, this function can also find repeats along the highest
% dimension, such as finding repeated planes in a 3D matrix.

% 04-02-09 DN  Wrote it.

psychassert(dim>0 && dim<=ndims(in),'dim input argument outside possible range (1:ndims(in)')
if ndims(in)>2
    psychassert(dim==ndims(in),'if input has more than 2 dimensions, repetitions can only be found along the highest dimension')
end

thediff = diff(in,[],dim);

if dim==1 && ndims(in)==2
    thediff = thediff.';
else
    nrows   = prod(AltSize(thediff,1:dim-1));
    rest    = num2cell(AltSize(thediff,dim:ndims(thediff)));

    thediff = reshape(thediff,nrows,rest{:});
end

thesum = sum(thediff,1);

if isvector(thesum)
    % if we test over the highest dimension, thesum will always be a vector
    % and we can simply return the indices to doubles in that highest
    % dimension
    inds = find(thesum==0)+1;
end
