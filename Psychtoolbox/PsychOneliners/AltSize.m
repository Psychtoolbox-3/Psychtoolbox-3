function varargout = AltSize(in,arg)
% varargout = AltSize(in,arg)
%
% extends size()'s functionality to support querying the size of multiple
% dimensions of a variable in one call.
% requested dimensions can be repeated and may (of course) be singleton
% number of output arguments must match number of requested dimension sizes
% or be one, in which case a vector is returned
% 
% example:
%     in = ones(1,2,3,4,5);
%     [d e f g h i] = AltSize(in,[1 3 2 8 4 2])
%     d =
%          1
%     e =
%          3
%     f =
%          2
%     g =
%          1
%     h =
%          4
%     i =
%          2

% DN 2008

if nargin==1
    % standaard size case, dispatch to build-in size
    if nargout==0||nargout==1
        varargout               = {size(in)};
    else
        [varargout{1:nargout}]  = size(in);
    end
elseif nargin==2 && isscalar(arg)
    % ook standaard size geval
    psychassert(nargout==0||nargout==1,'Unknown command option.')
    varargout                   = {size(in,arg)};

else
    psychassert(nargin==2,'two input arguments must be specified');
    psychassert(nargout==0||nargout==1||nargout==length(arg),'number of output variables must be 1 or match the number of requested dimensions');
    out         = zeros(1,length(arg));

    sizes       = size(in);

    qone        = arg>ndims(in);
    out(qone)   = 1;
    out(~qone)  = sizes(arg(~qone));

    if nargout==0||nargout==1
        varargout   = {out};
    else
        varargout   = num2cell(out);
    end
end
