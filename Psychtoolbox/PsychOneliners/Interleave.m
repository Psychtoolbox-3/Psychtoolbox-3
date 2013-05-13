function [vec] = Interleave(varargin)
% [vec] = Interleave(varargin)
% Interleaves any number of vectors and scalars
% input is taken from consequetive elements of each input variable until a
% variable runs out.
% Output will be a columnvector.
% Output will be cell if any of the inputs is cell.
% Output will be char if none of the input is cell and any of the inputs is
% char
% 
% a   = [A A A], b= [B B], c= [C C C C], d=D
% out = Interleave(a,b,c,d)
% out = [A B C D A B C A C C];
%
% DN 06-11-2007
% DN 23-01-2008 updated to support scalars and conversion to char
% DN 28-04-2008 simplified input checking
% DN 28-05-2008 bugfix mixing numeric and char
% DN 09-06-2011 Now works fine with empty inputs. Few changes for compatibility
%               with Octave, might have more up ahead.

% remove empties
varargin(cellfun(@isempty,varargin)) = [];
% check if we have any input
if isempty(varargin)
    vec = [];
    return;
end
% we have work to do, first some other input checks
for p=1:length(varargin)  % cannot use nargin as we just possibly deleted input arguments
    psychassert(isvector(varargin{p}) || isscalar(varargin{p}),'not all inputs are vectors or scalars');
    psychassert(isnumeric(varargin{p}) || iscell(varargin{p}) || ischar(varargin{p}),'not all inputs are numeric, cell or char')
end

% remove empty inputs
qEmpty = cellfun(@isempty,varargin);
varargin(qEmpty) = [];

% check what input we have
qchar = any(cellfun(@ischar,varargin));
qcell = any(cellfun(@iscell,varargin));
qnum  = cellfun(@isnumeric,varargin);

if any(qnum) && qchar && ~qcell
    % convert all numeric inputs to char
    inds = find(qnum);
    for p=1:length(inds)
        temp = num2cell(varargin{inds(p)});
        varargin{inds(p)} = cellfun(@num2str,temp,'UniformOutput',false);
    end
    % pack all other inputs in cell as well for Octave
    inds = find(~qnum);
    for p=1:length(inds)
        varargin{inds(p)} = num2cell(varargin{inds(p)});
    end
end

for p=1:numel(varargin)
    r(p,1:numel(varargin{p})) = num2cell(varargin{p});
end
vec=cat(2,r{:});

if any(qnum) && qchar && ~qcell
    % last step
    vec=[vec{:}];
end
