function mat = RandDim(in,dim)
% mat = RandDim(in,dim)
% randomizes dimension DIM of matrix IN
% IN can have any number of dimensions
%
% DN 2008 Wrote it

dimsize     = size(in,dim);

if dimsize==1 || dim>ndims(in)
    mat = in;
    return;
end

v(1:ndims(in))  = {':'};                % Create the index that will transform x. a(':') works just as well as a(:)
v{dim}          = randperm(dimsize);    % Throw in random order indices for the dim-to-be-randomized.
mat             = in(v{:});             % Use this indexing cell
