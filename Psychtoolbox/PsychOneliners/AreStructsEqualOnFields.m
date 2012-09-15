function areEqual = AreStructsEqualOnFields(struct1,struct2,theFields)
% areEqual = AreStructsEqualOnFields(struct1,struct2,theFields)
%
% Returns 1 if two structs share the same value on each of the passed
% fields, 0 otherwise.  The equality of each field in the two structs is
% checked with a call to isequal()
%
% 5/1/05     dhb, jmk   Handle cell and struct fields, a little bit.
% 2012/06/12 DN         Can now handle fields of any data type supported by
%                       isequal(). Added some input checks


psychassert(isscalar(struct1) && isscalar(struct2),'structs must be scalar');

if ~iscell(theFields)
    theFields = {theFields};
end

nFields = length(theFields);
areEqual = 0;

% Loop over fields.  Return on any indicator of non-equality.
% If we make it out the bottom, then set return value to 1
for i = 1:nFields
    % If either struct is missing the passed field, we say not equal.
    if ~isfield(struct1,theFields{i}) || ~isfield(struct2,theFields{i})
        return;
    end
    
    if ~isequal(struct1.(theFields{i}),struct2.(theFields{i}))
        return;
    end
end

% If we arrive here, they are equal
areEqual = 1;
