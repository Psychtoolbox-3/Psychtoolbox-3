function c = MergeCell(varargin)
% Merges contents of multiple cell arrays into one big cell array
%
% c = MergeCell(varargin)
% 
% MergeCell takes any number of cell vectors (contain the same datatype)
% and concatenates their contents into one big cellvector
% Any signleton inputs are expanded as needed, these inputs can be 1x1
% cells or the contained datatype itself
%
% Example: To add some information about fit to plot legend, we'd wan't to
% append information about the fit to the line labels:
%
% linelbls = {'line a','line b','line c'};
% chi2_r   = [.9 1 4.2];
% chi2_rtxt= arrayfun(@(x) sprintf('slope: %.3f',x),chi2_r,'UniformOutput',false);
%
% leglbls  = MergeCell(linelbls,', chi^2_r: ',chi2_rtxt);
% % leglbls  = MergeCell(linelbls,{', chi^2_r: '},chi2_rtxt); would be equivalent
%

% input checks
qVec    = cellfun(@isvector,varargin); % also true for scalars :)
assert(all(qVec),'At least one input is not a vector');

% ensure all inputs are cell (wrap in cell if not)
qCell   = cellfun(@iscell,varargin);
temp    = num2cell(varargin(~qCell));
varargin(~qCell) = temp;

% get lengths and check
len     = cellfun(@length,varargin);
assert(length(unique(len(len>1)))<2,'Inputs not all of matching or scalar length\nLengths: %s',num2str(len));
len     = max(len);

% check all same datatype
type    = cellfun(@(x) class(x{1}),varargin,'UniformOutput',false);
if length(unique(type)) > 1
    fprintf('Data types of inputs:\n');
    fprintf('  %s\n',type{:});
    error('Inputs not all of same data type ^^');
end

% expand all inputs to same length and ensure all column vectors
for p = 1:length(varargin)
    if length(varargin{p}) == 1
        % expand scalar
        varargin{p} = repmat(varargin{p},len,1);
    elseif size(varargin{p},2) > 1
        % ensure column vector
        varargin{p} = varargin{p}.';
    end
end

% concatenate
c = cellfun(@(varargin) cat(2,varargin{:}),varargin{:},'UniformOutput',false);