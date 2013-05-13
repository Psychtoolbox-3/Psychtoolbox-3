function [x,y,z] = CreateUniformDotsIn3DFrustum(ndots,FOV,aspectr,depthrangen,depthrangef,eyeheight)
% [x,y,z] = CreateUniformDotsIn3DFrustum(ndots,FOV,aspectr,depthrangen,depthrangef,eyeheight)
%
% Sample dots in frustum uniformly, creating a cloud tightly fitting in the
% frustum. When the optional 6th parameter eyeheight is given, dots will be
% uniformly sampled on a ground plane at -eyeheight.
%
% z is not sampled from a uniform distribution, but from a parabolic
% distribution as the area of cross sections of the frustum is a quadratic
% function of the depth plane's depth ( (z*2*tan(FOV/2))^2 * aspectr )
%
% Here, I use Inverse transform sampling to transform a uniform random
% variable into the quadratic shape random variable
% see Luc Devroye. Non-Uniform Random Variate Generation. New York:
% Springer-Verlag, 1986. Chapter 2
% (http://cg.scs.carleton.ca/~luc/chapter_two.pdf)
% compile following in latex to see full derivation:
% ----
% \documentclass[12pt,a4paper]{minimal}
% \usepackage{amsmath}        % math
% 
% \begin{document}
% \textbf{Derivation}:\\
% Use pdf related to cross-section surface of frustum:
% \(\left(2 z \tan\left(\frac{FOV}{2}\right)\right)^2 aspectr\)\\
% \(z_1\) is the distance of the near depth plane\\
% \(z_2\) is the distance of the far depth plane\\
% \(y\) is a uniform random variable\\
% Given: \(F(z_2)=1\) and \(F(z_1)=0\).
% \begin{align}
%     F(z) &= \int\limits^z_{z_1} k z^2 \, \mathrm{d}z\\
%     F(z) &= \frac{k}{3} \left(z^3 - z_1^3\right)\\
%     k    &= \frac{3}{z_2^3-z_1^3}\\
%     F(z) &= \frac{z^3-z_1^3}{z_2^3-z_1^3}
% \end{align}
% Substitute \(y\) for \(F(z)\) and factor out \(z\):
% \begin{align}
%     z^3  &= y\left(z_2^3-z_1^3\right) + z_1^3\\
%     z    &= \sqrt[3]{y\left(z_2^3-z_1^3\right) + z_1^3}
% \end{align}\\
% For a ground plane, the width of the frustum at depth \(z\) is given by
% \(2 z \tan\left(\frac{FOV}{2}\right) aspectr\).\\
% By following the same inverse transform sampling steps, we would end up
% with 
% \[ z = \sqrt{y\left(z_2^2-z_1^2\right) + z_1^2} \]
% for the linear distribution between \(z_1\) and \(z_2\).
% 
% \end{document}
% ----

% 2008       DN  Wrote it.
% 2009-06-06 DN  Changed input check to allow for vector near and far
%                depthrange and allowed near and far depthrange to be
%                the same, to place dots at exactly that depth
% 2010-06-08 DN  Added support for uniform groundplanes

psychassert(all(depthrangen<=depthrangef),'Near clipping plane should be closer than far clipping plane');

u   = RandLim([1,ndots],0,1);                                           % get uniform random variable
if nargin>5
    % eye height specified, ground plane
    z   = -(u.*(depthrangef.^2-depthrangen.^2)+depthrangen.^2).^(1/2);      % transform to linear distribution (negate as depth postiion is a negative number
else
    % cloud
    z   = -(u.*(depthrangef.^3-depthrangen.^3)+depthrangen.^3).^(1/3);      % transform to parabolar distribution (negate as depth postiion is a negative number
end

yrs = -z*tand(FOV/2);
if nargin>5
    y   = -eyeheight*ones(1,ndots);
else
    y   = RandLim([1,ndots],-yrs,yrs);
end
xrs = yrs * aspectr;
x   = RandLim([1,ndots],-xrs,xrs);

% for a uniform sample over z:
% z   = -RandLim([1,ndots],depthrangen,depthrangef);
