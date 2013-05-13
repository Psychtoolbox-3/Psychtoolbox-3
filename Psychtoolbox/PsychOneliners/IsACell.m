function bool = IsACell(input,fhndl)
% bool = IsACell(input,fhndl)
% 
% checks if there is any element in the cell INPUT that satisfies the
% logical condition in function FHNDL
% works recursively: cells in cells (etc) also checked
% the function FHNDL must return a scalar
%
% examples:
% checks if there is any struct in the cell:
% fhndl = @isstruct
% check if there is any element in the cell extending
% over more than 2 dimensions:
% fhndl = @(x)ndims(x)>2
%
% DN    2008-05-28

psychassert(iscell(input),'IsACell: Input must be cell, not %s',class(input));

if fhndl(input)
    bool = true;
    return;
else
    bool = false;
end

qtest   = cellfun(fhndl,input);
qcell   = cellfun(@iscell,input);

if AnyAll(qtest)
    % cell found with data in it that satisfies the logical condition in
    % FHNDL
    bool    = true;
elseif AnyAll(qcell)
    % cells found in the cell, process also
    [lind]       = find(qcell);
    for p=1:length(lind)
        bool    = IsACell(input{lind(p)},fhndl);
        if bool
            return;
        end
    end
end
