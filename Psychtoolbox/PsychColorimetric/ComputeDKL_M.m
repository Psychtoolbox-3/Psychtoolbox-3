function M = ComputeDKL_M(bg,T_cones,T_Y)
% M = ComputeDKL_M(bg,T_cones,T_Y)
% 
% Compute the matrix that converts between incremental cone
% coordinates and DKL space.  The order of
% the coordinates in the DKL column vectors is (Lum, RG, S)
%
% The code follows that published by Brainard
% as an appendix to Human Color Vision by Kaiser
% and Boynton, but has been generalized to work
% for passed cone fundamentals and luminosity 
% function.
%
% These should be passed in standard Psychtoolbox
% form (LMS in rows of T_cones; luminosity function
% as single row vector T_Y).
%
% The argument bg is the LMS cone coordinates of the
% background that defines the space.
%
% See DKLDemo for proper use of this function.  Also
% DKLToConeInc and ConeIncToDKL.
%
% 8/30/96   dhb  Pulled it out.
% 4/9/05    dhb  Allow passing of cones and luminance to be used.
% 11/17/05  dhb  Require passing of cones and luminance.
%           dhb  Fixed definition of M_raw to handle arbitrary L,M scaling.
% 10/5/12   dhb  Comment specifying coordinate system convention.  Supress extraneous printout.

% If cones and luminance are passed, find how L and
% M cone incrments sum to best approximate change in
% luminance.
if (nargin == 3)
	T_LM = T_cones(1:2,:);
	LMLumWeights = T_LM'\T_Y';
else
    fprintf('ComputeDKL_M now requires explicit specification\n');
    fprintf('of cone fundamentals and luminosity function\n');
    fprintf('See DKLDemo\n');
    error('');   
end

% Set M_raw as in equation A.4.9.
% This is found by plugging the background
% values into equation A.4.8.  Different
% backgrounds produce different matrices.
% The Matlab notation below just 
% fills the desired 3-by-3 matrix.
%
% Note that A.4.8 in the Brainard chapter contains
% a typo: the row 1 col 3 entry of the matrix should
% be 0, not 1.  Also, at the top of page 571, there is
% an erroneous negative sign in front of the term
% for W_S-Lum,S.
%
% Finally, A.4.8 as given in the chatper assumes
% that Lum = L + M.  The formula below generalizes
% to arbitrary scaling.
M_raw = [ LMLumWeights(1) LMLumWeights(2) 0 ; ...
			1 -bg(1)/bg(2) 0 ; ...
			-LMLumWeights(1) -LMLumWeights(2) (LMLumWeights(1)*bg(1)+LMLumWeights(2)*bg(2))/bg(3) ];

% Compute the inverse of M for
% equation A.4.10.  The Matlab inv() function
% computes the matrix inverse of its argument.
M_raw_inv = inv(M_raw);

% Find the three isolating stimuli as
% the columns of M_inv_raw.  The Matlab
% notation X(:,i) extracts the i-th column
% of the matrix X.
isochrom_raw = M_raw_inv(:,1);
rgisolum_raw = M_raw_inv(:,2);
sisolum_raw = M_raw_inv(:,3);

% Find the pooled cone contrast of each
% of these.  The Matlab norm() function returns
% the vector length of its argument.  The Matlab
% ./ operation represents entry-by-entry division.
isochrom_raw_pooled = norm(isochrom_raw ./ bg);
rgisolum_raw_pooled = norm(rgisolum_raw ./ bg);
sisolum_raw_pooled = norm(sisolum_raw ./ bg);

% Scale each mechanism isolating
% modulation by its pooled contrast to obtain
% mechanism isolating modulations that have
% unit length.
isochrom_unit = isochrom_raw / isochrom_raw_pooled;
rgisolum_unit = rgisolum_raw / rgisolum_raw_pooled;
sisolum_unit = sisolum_raw / sisolum_raw_pooled;

% Compute the values of the normalizing
% constants by plugging the unit isolating stimuli
% into A.4.9 and seeing what we get.  Each vector
% should have only one non-zero entry.  The size
% of the entry is the response of the unscaled
% mechanism to the stimulus that should give unit
% response.
lum_resp_raw = M_raw*isochrom_unit;
l_minus_m_resp_raw = M_raw*rgisolum_unit;
s_minus_lum_resp_raw = M_raw*sisolum_unit;
					 
% We need to rescale the rows of M_raw
% so that we get unit response.  This means
% multiplying each row of M_raw by a constant.
% The easiest way to accomplish the multiplication
% is to form a diagonal matrix with the desired
% scalars on the diagonal.  These scalars are just
% the multiplicative inverses of the non-zero
% entries of the vectors obtained in the previous
% step.  The resulting matrix M provides the
% entries of A.4.11.  The three _resp vectors
% computed should be the three unit vectors
% (and they are).
D_rescale = [1/lum_resp_raw(1) 0 0 ; ...
						 0 1/l_minus_m_resp_raw(2) 0 ; ...
						 0 0 1/s_minus_lum_resp_raw(3) ];				 
M = D_rescale*M_raw;
lum_resp = M*isochrom_unit;
l_minus_m_resp = M*rgisolum_unit;
s_minus_lum_resp = M*sisolum_unit;

% Compute the inverse of M to obtain
% the matrix in equation A.4.12.
M_inv = inv(M);
 
