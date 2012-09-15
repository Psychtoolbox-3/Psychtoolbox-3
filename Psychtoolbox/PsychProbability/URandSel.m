function sel = URandSel(in,n)
% sel = RandSel(in,n)
%   if N is a scalar
%     randomly selects N elements from set IN, each element from IN will
%     be selected up to one time (so sampling without replacing).
%     Output is of size 1xN
%   if N is a vector
%     randomly selects prod(N) elements from set IN, without replacing. N
%     specifies the shape of the output, e.g., if N=[10 10 5], SEL will be
%     a 10x10x5 matrix
%   IN can be of any datatype and can have any number of dimensions and
%   elements

% 2008-08-27 DN Wrote it
% 2012-06-13 DN Can now also return shaped inputs through non-scalar n input

if isscalar(n)
    nElem = n;
    n     = [1,n];
else
    nElem = prod(n);
end

psychassert(nElem<=numel(in),'More return elements requested (prod(n): %d) than number of elements in input (numel(in): %d)',nElem,numel(in));

selind  = NRandPerm(numel(in),nElem);
sel     = reshape(in(selind),n);
