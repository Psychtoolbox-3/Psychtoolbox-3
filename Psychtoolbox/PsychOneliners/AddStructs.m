function base = AddStructs(base,varargin)
% OUTSTRUCT = AddStructs(BASE,A,B,...)
% Appends fields of structs A and B, etc to struct BASE

% 2008-11-30 DN  Wrote it
% 2012-06-12 DN  Added check that all inputs are scalar, this only works
%                with scalars...

psychassert(isstruct(base),'1th input argument is not a struct');
psychassert(isscalar(base) && all(cellfun(@isscalar,varargin)),'all inputs must be scalar structs');

for p=1:length(varargin)
    psychassert(isstruct(varargin{p}),'%dth input argument is not a struct',p+1);
    
    flds = fieldnames(varargin{p});
    for q=1:length(flds)
        psychassert(~isfield(base,flds{q}),'Field ''%s'' of input %d already exists in base struct',flds{q},p+1);
        base.(flds{q}) = varargin{p}.(flds{q});
    end
end
