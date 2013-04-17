% Psychtoolbox:PsychColorimetricData:PsychColorimetricMatFiles.
%
% help Psychtoolbox % For an overview, triple-click me & hit enter.
% help PsychDemos   % For demos, triple-click me & hit enter.
% help PsychColorimetricData
% help PsychColorimetric
%
% This folder holds colorimetric data in .mat file form.
% All data files are in a standard format.
%
% A very useful source for on-line
% colorimetric data is the CVRL database: 
%   http://cvrl.ioo.ucl.ac.uk/
% Many of the functions used here were downloaded from
% that source and then splined to the wavelength sampling
% used here (with extension by zeros).
%
%   B_xxx files contain basis functions.  The basis functions themselves
%     are in the columns of a matrix with name B_xxx.  There is also
%     a 3 by 1 row vector S_xxx that contains the wavelength sampling
%     information in the form [start delta numberSamples] where start
%     and delta are in nanometers.
%
%   den_xxx files contain optical density data.  In log units.  To
%     convert density values to transmittance, take 10^(-den).  There
%     is also an S_xxx vector.  Curiously, these are in column vectors
%     in the .mat files.
%
%   spd_xxx files contain spectral power distributions.  Data are in
%     columns of matrix with name spd_xxx.  There is also a S_xxx
%     vector.
%
%   sur_xxx files contain surface reflectance functions (range 0-1).  Data
%     are in columns of matrix with name sur_xxx.  There is also an S_xxx
%     vector.
%
%   T_xxx contain color matching functions or spectral sensitivities.  Data
%      are in rows of a matrix with name T_xxx.  There is also an S_xxx
%      vector.  All are in energy units, so that multiplying by spectra in
%      energy (not quanta) gives desired result (often proportional to
%      isomerization rate).
%
%  Specific data files are listed below.  Most, but not all, are sampled
%  between 380 nm and 780 nm at 5 nm intervals (S = [380 5 81]).  This is
%  the CIE standard.  Good coding practice requires using the S_xxx vector
%  loaded with the data and splining to the wavelength sampling you want to work
%  in.  If I were to start this database again, I would have kept each function
%  at the resolution of its source.
%
%  In some cases, the original data were interpolated or extraploated (with zeros)
%  to put them data onto the CIE standard [380 5 81] wavelength.
%
%  See also: EnergyToQuanta, QuantaToEnergy, MakeItS, MakeItWls, MakeItStruct, 
%    SplineSpd, SplineSrf, SplineCmf.
%
%   B_cieday            - CIE daylight basis functions.
%   B_cohen             - Cohens basis functions for Munsell surfaces.
%   B_monitor           - Basis functions for a color monitor.
%   B_nickerson         - Basis functions for Munsell surfaces.
%   B_roomillum         - Basis functions for illuminants in Brainard's room.
%   B_vrhel             - Basis functions for Vrhel surface measurements.
%   den_lens_ws         - Relative lens density data (re 700 nm).  W&S, Table 1(2.4.6), p. 109.
%                       -   This is the first data set in the table, not the Norren and Vos 
%                       -   data.  It is for an open pupil.
%   den_lens_cie_1      - Part one of CIE component lens density function. CIE 170-1:2006, Table 6.10
%   den_lens_cie_2      - Part two of CIE component lens density function. CIE 170-1:2006, Table 6.10
%   den_lens_ssf        - Stockman-Sharpe-Fach (1999) lens optical density spectrum.
%                       -   See CVRL database, CIE 170-1:2006, Table 6.10, 32 yo, pupil <= 3 degrees.
%                       -   This is also the sum of den_lens_cie_1 and den_lens_cie_2
%   den_mac_bone        - Macular pigment density from Bone et al. (1992).  See CVRL database, CIE 170-1:2006, Table 6.4, 2-deg.
%   den_mac_vos         - Macular pigment density from Vos.  See CVRL database.
%   den_mac_ws          - Macular pigment density from W&S, Table 2(2.4.6), p. 112.
%   spd_appratusrel     - Relative spectrum from a monitor.  Used by IsomerizationInDishDemo.
%   spd_CIEA            - Spectral power distribtion for CIE illuminant A.
%   spd_CIEC            - Spectral power distribution for CIE illuminant C.
%   spd_D65             - Spectral power distribution for CIE illuminant D65.
%   spd_flourescent     - Spectral power distribution for some flourescent lamp.
%   spd_incanCC         - Spectral power distribution for some incandescent lamp.
%   spd_phillybright    - Direct bright sunlight measured through window and off of a piece of white paper towel
%                       -   on the floor of DB's office in Philly, March 2013.
%                       -   Measurements made with PR-650, power in Watts/[m2-sr-wlband].
%   spd_xenonArc        - Spectral power distribution for some xenon arc lamp.
%   spd_xenonFlash      - Spectral power distribuiton for some xenon flash tube.
%   sur_nickerson       - The Nickerson measurements of the Munsell papers.
%   sur_macbeth         - Reflectance of Macbeth color checker (not accurate, needs updating).
%   sur_vrhel           - Reflectances measured by Vrhel.
%   T_CIE_Y2            - CIE physiologically relevant 2-degree luminosity function.  See CVRL database.
%   T_CIE_Y10           - CIE physiologically relevant 10-degree luminosity function.  See CVRL database.
%   T_cones_smj         - Stockman-MacLeod-Johnson cone fundamentals.  See CVRL database.
%   T_cones_smj10       - Stockman-MacLeod-Johnson 10-degree cone fundamentals.  See CVRL database.
%   T_cones_ss2         - Stockman-Sharpe (2000) 2-degree cone fundamentals.  See CVRL database.
%   T_cones_ss10        - Stockman-Sharpe (2000) 10-degree cone fundamentals.  See CVRL database.
%   T_cones_sp          - Smith-Pokorny cone fundamentals.
%   T_DCS200            - Sensitivities of a Kodak DCS-200 color camera.
%   T_ground            - Not entirely sure what this is, but it might be ground squirrel receptor sensitivities.
%   T_Lanom             - Demarco et al. anomolous L cone sensitivity.
%   T_log10coneabsorbance_ss - Stockman-Sharpe (2000) log10 LMS cone photopigment absorbance.
%                       -   See CVRL database, CIE 170-1:2006, Table 6.6.
%                       -   Some S-cone values were unspecified for wls > 615 nm in the table.
%                       -   These were filled in here by linear extrapolation.
%                       -   Note that you want to raise 10 to these numbers
%                       -   to get absorbance, which itself is a log-like quantity.
%   T_Manom             - Demarco et al. anomolous M cone sensitivity.
%   T_photopigments_ss  - Removed.  Use T_log10coneabsorbance and raise 10 to it.
%   T_rods              - CIE scotopic luminous efficiency function.
%   T_stiles2           - Stiles 2-degree color matching functions.
%   T_ss2000_Y2         - Stockman-Sharpe (2000) 2-degree photopic luminance efficiency function.  See CVRL database.
%   T_vos1978_Y         - Judd-Vos 1978 photopic luminance efficiency function.
%   T_xyz1931           - CIE 1931 color matching functions.
%   T_xyz1964           - CIE 1964 supplemental color matching functions.
%   T_xyzJuddVos        - Judd-Vos modified color matching functions.

