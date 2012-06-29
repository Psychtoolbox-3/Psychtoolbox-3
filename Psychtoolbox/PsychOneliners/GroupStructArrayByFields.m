function theGroupedArray = GroupStructArrayByFields(theStructArray,theFields)
% theGroupedArray = GroupStructArrayByFields(theStructArray,theFields)
%
% Group together the members of a struct array that share the same values
% in the passed fields.
%
% This is useful for sorting/grouping elements in a struct array on
% parameters that signify membership of trials to a particular condition.
%
% 7/21/03  dhb  Wrote it.

theGroupedArray = {};
nStructs = length(theStructArray);

% The first passed structure is equal to itself
nGroups = 1;
theGroupedArray{1} = theStructArray(1);

% Put structures into groups, creating new ones as necessary.
for i = 2:nStructs
	didIt = 0;
	for j = 1:nGroups
		if AreStructsEqualOnFields(theStructArray(i),theGroupedArray{j}(1),theFields)
			theGroupedArray{j} = [theGroupedArray{j} theStructArray(i)];
			didIt = 1;
			break;
		end
	end
	if (~didIt)
		nGroups = nGroups+1;
		theGroupedArray{nGroups} = theStructArray(i);
	end
end
