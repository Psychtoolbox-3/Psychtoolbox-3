function [B] = MakeGaussBasis(wls,means,vars);
% [B] = MakeGaussBasis(wls,means,vars);
%
% Make a basis set of shifted Gaussians.  Area
% under each Gaussian is unity.
%
% 1/23/96  dhb	Wrote it.
% 11/20/96 dhb  Make area unity, not maximum which is how it was. 
% 12/3/99  dhb  Change documentation sds -> vars as this is what it does.
% 10/12/20 ms   Replaced call to `NormalPDF()` with call to `normpdf()`.

% Allocate space
nBases = length(means);
nWls = length(wls);
B = zeros(nWls,nBases);
if (length(vars) == 1)
	vars = vars*ones(nBases,1);
end

% Make the bases
for i = 1:nBases
	temp = normpdf(wls,means(i),vars(i));
	temp = temp/sum(temp);
	B(:,i) = temp;
end

