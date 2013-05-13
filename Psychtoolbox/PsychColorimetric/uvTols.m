function ls = uvTols(uv)
% ls = uvTolsy(uv)
%
% Convert CIE u'v' chromaticity to cone chromaticity ls, L/(L+M+S), S/(L+M+S).
%
% Uses regression conversion matrix based on Judd-Vos XYZ and
% Smith-Pokorny cone fundamentals to get from XYX to LMS.  This
% is an exact linear transformation and so you don't get as many
% weird little numerical things happening when you apply this to
% uv for spectral lights.
% 
% 3/17/04  dhb		Wrote it.
% 05/06/11 dhb      Make function name in file match actual function name.

% Compute the conversion matrix
% load T_xyzJuddVos;
% load T_cones_sp;
% M = ((T_xyzJuddVos')\(T_cones_sp'))';

% Define the conversion matrix.
M = [  0.24352943081928   0.85222450346271  -0.05154899613656 ;
  	  -0.39546852068084   1.16421653706998   0.08383540224414 ;
   		-0.00015007798271   0.00019129501713   0.61879471228244 ];

% Pop it all through
uvY = [uv ; ones(1,size(uv,2))];
XYZ = uvYToXYZ(uvY);
LMS = M*XYZ;
nCols = size(LMS,2);
ls = zeros(2,nCols);
for j = 1:nCols
	ls(1,j) = LMS(1,j)/sum(LMS(:,j));
	ls(2,j) = LMS(3,j)/sum(LMS(:,j));
end


uvTo
