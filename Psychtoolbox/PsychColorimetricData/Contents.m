% Psychtoolbox:PsychColorimetricData.
%
% Colorimetric data and some data driven colorimetric calculations.
%   i) Subfolder ColorimetricMatFiles contains .mat files with basic data.
%
%   ii) Functions here return various standard data -- many routines
%   take a string argument that allows access to estimates from different
%   sources.
%
% The data routines here are closely coupled to the calculation routines
% in PsychColorimetric.  The conceptual dividing line between the two folders is
% somewhat blurry.
%
% help Psychtoolbox              % For an overview, triple-click me & hit enter.
% help PsychColorimetricMatFiles % For list of data .mat files, triple-click me & hit enter.
% help PsychColorimetric         % For colorimetric calculations, triple-click me & hit enter.
%
%   BaylorNomogram            - Baylor et al. photopigment nomogram.
%   ComputeCIEConeFundamentals - Compute cone fundamentals according to CIE 170-1:2006
%   ComputeRawConeFundamentals - Compute cone fundamentals from specification of various components
%   DawisNomogram             - Dawis(1981) photopigment nomogram.
%   DefaultPhotoreceptors     - Set default values for photoreceptors structure.
%   DegreesToRetinalMM        - Convert retinal extents from degrees to mm of retina.
%   EyeLength                 - Return estimate of distance between nodal point and retina.
%   FillInPhotoreceptors      - Convert from data source specification to numeric values in photoreceptors structure.
%   FitConeFundamentalsWithNomogram - Try to fit CIE cone fundamentals with absorbance from various nomograms.
%   GovardovskiiNomogram      - Govardoskii et al. (2000) A1 photopigment nomogram.
%   LambNomogram              - Lamb's (1995) photopigment nomogram.
%   LensTransmittance         - Return transmittance of human lens.
%   MacularTransmittance      - Return transmittance of human macular pigment.
%   PhotopigmantAxialDensity  - Estimate of peak optical density.
%   PhotopigmantNomogram      - Encapsulate available nomogram computations.
%   PhotopigmentSpecificDensity - Estimates of specific density.
%   PhotopigmentQuantalEfficiency - Estimates of photopigment quantal efficiency.
%   PhotoreceptorDimensions   - Estimates of various photoreceptor dimensions.
%   PrintPhotoreceptors       - Print out what is in photoreceptor structure.
%   PupilDiameterFromLum      - Estimate pupil diameter from luminance.
%   RetinalMMToDegrees        - Convert foveal retinal extent from mm of retina to degrees.
%   StockmanSharpeNomogram    - Stockman/Sharpe photopigment nomogram (not finished).
%   ValetonVanNorrenParams    - Return parameters of the Valeton-Van Norren model.
  
% Copyright (c) 1996-2003 by Denis Pelli & David Brainard

