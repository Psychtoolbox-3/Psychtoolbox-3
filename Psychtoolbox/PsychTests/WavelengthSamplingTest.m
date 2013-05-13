% WavelengthSamplingTest
%
% There are three formats used to represent wavelength
% sampling in the Psychtoolbox.  
%   S-format      - a 3 by 1 row vector [start step numberSamples].
%   wls-format    - a list of evenly spaced wavelengths.
%   struct-format - a structure with fields start, step, and numberSamples.
%
% The S-format predates the availability of structures in MATLAB, but
% is used widely in extant code.  Low level conversion routines transform
% between these representations, and in particular routines MakeItS,
% MakeItWls, and MakeItStruct take any of the three formats and return
% one.  By calling one of these before using wavelength format information,
% it is possible to write code that is compatible with all three.
%
% 7/11/03  dhb  Wrote this test program.

% Clear
clear all

% Generate some test data
S = [380 5 81];
wls = MakeItWls(S);
wlstruct = MakeItStruct(S);

% Convert all formats to wls in various ways.
wls1 = SToWls(S);
wls2 = MakeItWls(S);
wls3 = MakeItWls(wls);
wls4 = MakeItWls(wlstruct);
if (any(wls ~= wls1))
	fprintf('Wavelength conversion error: wls1 mismatched to wls\n');
end
if (any(wls ~= wls2))
	fprintf('Wavelength conversion error: wls2 mismatched to wls\n');
end
if (any(wls ~= wls3))
	fprintf('Wavelength conversion error: wls3 mismatched to wls\n');
end
if (any(wls ~= wls4))
	fprintf('Wavelength conversion error: wls4 mismatched to wls\n');
end

% Convert all formats to S in various ways.
S1 = WlsToS(wls);
S2 = MakeItS(S);
S3 = MakeItS(wls);
S4 = MakeItS(wlstruct);
if (any(S ~= S1))
	fprintf('Wavelength conversion error: S1 mismatched to S\n');
end
if (any(S ~= S2))
	fprintf('Wavelength conversion error: S2 mismatched to S\n');
end
if (any(S ~= S3))
	fprintf('Wavelength conversion error: S3 mismatched to S\n');
end
if (any(S ~= S4))
	fprintf('Wavelength conversion error: S4 mismatched to S\n');
end

% Convert all formats to S in various ways.
wlstruct1 = MakeItStruct(S);
wlstruct2 = MakeItStruct(wls);
wlstruct3 = MakeItStruct(wlstruct);
if (wlstruct.start ~= wlstruct1.start || ...
	wlstruct.step ~= wlstruct1.step || ...
	wlstruct.numberSamples ~= wlstruct1.numberSamples)
	fprintf('Wavelength conversion error: wlstruct1 mismatched to wlstruct\n');
end
if (wlstruct.start ~= wlstruct2.start || ...
	wlstruct.step ~= wlstruct2.step || ...
	wlstruct.numberSamples ~= wlstruct2.numberSamples)
	fprintf('Wavelength conversion error: wlstruct2 mismatched to wlstruct\n');
end
if (wlstruct.start ~= wlstruct3.start || ...
	wlstruct.step ~= wlstruct3.step || ...
	wlstruct.numberSamples ~= wlstruct3.numberSamples)
	fprintf('Wavelength conversion error: wlstruct3 mismatched to wlstruct\n');
end
