function inputVal = GetWithDefault(prompt,defaultVal)
% inputVal = GetWithDefault(prompt,defaultVal)
%
% Prompt for a number or string, with a default returned if user
% hits return.
%
% 4/3/10  dhb  Wrote it.

if (ischar(defaultVal))
    inputVal = input(sprintf([prompt ' [%s]: '],defaultVal),'s');
else
    inputVal = input(sprintf([prompt ' [%g]: '],defaultVal));
end
if (isempty(inputVal))
    inputVal = defaultVal;
end
