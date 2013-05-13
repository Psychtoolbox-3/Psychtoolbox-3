function in = DeEmptify(in,column)
% in = DeEmptify(in,column)
% deletes empty cells or rows from cellarray:
%
% - if only IN is specified, IN has to be a vector. IN will be returned
%   with all empty cells deleted
% - if IN is a matrix, COLUMN must be specified. Rows are only
%   deleted from IN when an empty cell is encountered in a column
%   specified in COLUMN. COLUMN can be a vector
%
% Example:
%   DeEmptify({'','r','','re'})
%   ans = 
%       'r'    're'
%   
%   a = {'p' 'r'  ''; ...
%        '' 're' 'r'}
%
%   DeEmptify(a,3)
%   ans = 
%        ''    're'    'r'
%
%   DeEmptify(a,1)
%   ans = 
%        'p'   'r'     ''
%
%   DeEmptify(a,1) or DeEmptify(a,[1 3])
%   ans = 
%       Empty cell array: 0-by-3

% DN 2008
% DN 2008-07-29 Simplified and included support for multiple columns

psychassert(nargin==1||nargin==2,'Provide 1 or 2 inputs')
if nargin==1
    psychassert(isvector(in),'Input has to be a vector when one input is provided');
    column = ':';
else
    psychassert(ndims(in)==2,'Input must be a 2-D matrix of cells')
    psychassert(any(column<=size(in,2)) && any(column>0),'"%s" contains an invalid column number',num2str(column))
end

qempty  = cellfun(@isempty,in(:,column));

if isvector(column)
    qempty = sum(qempty,2)>0;
end

if nargin==2
    in = in(~qempty,:);
else
    in = in(~qempty);
end
