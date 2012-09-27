function XYZ = LjgToXYZ(Ljg)
% XYZ = LjgToXYZ(Ljg)
%
% Convert OSA Ljg to XYZ (10 degree).  Works by using numerical
% search to invert XYZToLjg.  See XYZToLjg for details on
% formulae used.  See also TestOSAUCS.
%
% This can return imaginary values if you pass XYZ values
% that are outside reasonable physical gamut limits.
%
% 3/27/01  dhb      Wrote it.
% 3/4/05   dhb	    Handle new version of optimization toolbox, too.
% 9/23/12  dhb, ms  Update options for current Matlab versions.

% Check for needed optimization toolbox, and version.
if (exist('fmincon') == 2)
    % Search options
    if (verLessThan('optim','4.1'))
        error('Your version of the optimization toolbox is too old.  Update it.');
    end
    options = optimset('fmincon');
    options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
    
	% Search bounds -- XYZ must be positive.
	vlb = zeros(3,1);
	vub = [200 100 200]';
	
	% Do the search for each passed value.
	XYZ = zeros(size(Ljg));
	n = size(XYZ,2);
	for i = 1:n
		x0 = xyYToXYZ([.28 .28 30]');
		x = fmincon('LjgToXYZFun',x0,[],[],[],[],vlb,vub,[],options,Ljg(:,i));
		XYZ(:,i) = x;
	end

elseif (exist('constr') == 2)
	% Search options
	options = foptions;
	options(1) = 0;
	
	% Search bounds -- XYZ must be positive.
	vlb = zeros(3,1);
	vub = [200 100 200]';
	
	% Do the search for each passed value.
	XYZ = zeros(size(Ljg));
	n = size(XYZ,2);
	for i = 1:n
		x0 = xyYToXYZ([.28 .28 30]');
		x = constr('LjgToXYZFun',x0,options,vlb,vub,[],Ljg(:,i));
		XYZ(:,i) = x;
	end


else
	error('LjgToXYZ requires the optional Matlab Optimization Toolbox from Mathworks');
end



