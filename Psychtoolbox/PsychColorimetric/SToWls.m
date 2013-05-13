function wls = SToWls(S)
% wls = SToWls(S)
%
% Expand a [start delta n] description to an actual
% list of wavelengths.
% 
% 4/17/02  dhb  Handle degenerate cases of delta = 0 or n = 1.
% 7/11/03  dhb  Force S representation on input.

% Force to S format
S = MakeItS(S);

% Check validity
[m,n] = size(S);
if (m ~= 1 || n ~= 3)
  error('Passed list is not a [start delta n] description');
end
if (S(1) <= 0 || S(3) <=0)
  error('Passed list is not a [start delta n] description');
end

% Expand away
if (S(2) == 0 || S(3) == 1)
	wls = S(1);
else
	wls = (S(1):S(2):S(1)+(S(3)-1)*S(2))';
end

