function out = CombVec(varargin)
%CombVec Generate all possible combinations of input vectors.
%
%   CombVec(A1,A2,...) takes any number of inputs,
%      A1 - Matrix of N1 (column) vectors.
%      A2 - Matrix of N2 (column) vectors.
%      ...
%    and returns a matrix of (N1*N2*...) column vectors, where the columns
%    consist of all possibilities of A2 vectors, appended to
%    A1 vectors, etc.
%
%  Example
%  
%    a1 = [1 2];
%    a2 = [3 4; 3 4];
%    a3 = CombVec(a1,a2)
%    a3 = 
%        1     2     1     2
%        3     3     4     4
%        3     3     4     4

% 2008-08-06 DN  Wrote it, modification of CombVec in Matlab's Neural
%                Network Toolbox

if isempty(varargin)
    out = [];
else
    out = varargin{1};
    for i=2:length(varargin)
        cur = varargin{i};
        out = [copyb(out,size(cur,2)); copyi(cur,size(out,2))];
    end
end

%=========================================================
function b = copyb(mat,s)

[mr,mc] = size(mat);
inds    = 1:mc;
inds    = inds(ones(s,1),:).';
b       = mat(:,inds(:));

%=========================================================
function b = copyi(mat,s)

[mr,mc] = size(mat);
inds    = 1:mc;
inds    = inds(ones(s,1),:);
b       = mat(:,inds(:));
