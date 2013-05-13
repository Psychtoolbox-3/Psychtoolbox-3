function sel = RandSel(in,n)
% sel = RandSel(in,n)
%   if N is a scalar
%     randomly selects N elements from set IN, elements are replenished.
%     Output is of size 1xN
%   if N is a vector
%     randomly selects prod(N) elements from set IN, with replenishment. N
%     specifies the shape of the output, e.g., if N=[10 10 5], SEL will be
%     a 10x10x5 matrix
%   IN can be of any datatype and can have any number of dimensions and
%   elements

% DN 2008 wrote it
% 2008-08-06 Updated it to support n selections from set (DN)
% 2012-06-13 Can now also return shaped inputs through non-scalar n input (DN)

if isscalar(n)
    n = [1,n];
end

selind  = ceil(RandLim(n,0,numel(in)));
sel     = in(selind);
