function varargout = FindInd(in, k, mode)
% find that also does higher dims for output like [x,y,x,t] = findfull(4-D)
% (unlimited of course)
%
% steal from MinInd, and then edit MinInd and MaxInd

% DN 2008-07-29 Wrote it

% input/output checking
if nargin == 3
    switch mode
        case 'first'
            qfirst = true;
        case 'last'
            qfirst = false;
        otherwise
            error('mode %s not recognized',mode)
    end
elseif nargin == 2
    qfirst = true;
end

psychassert(nargout==0 || nargout==1 || nargout==ndims(in),'number of outputs must be one or equal to the number of dimensions of the input')

% do the work
inds                = find(in);

if nargin==2 || nargin==3
    if qfirst
        inds    = inds(1:k);
    else
        inds    = inds(end-k+1:end);
    end
end
    
[varargout{1:ndims(in)}]    = ind2sub(size(in),inds(:));

% output
if nargout==0 || nargout==1
    varargout   = {[varargout{:}]};
end
