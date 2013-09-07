function PrintPhotoreceptors(photoreceptors)
% PrintPhotoreceptors(photoreceptors)
%
% Print to command window an interpretable output
% of what is in a photoreceptors structure.
%
% See also DefaultPhotoreceptors, FillInPhotoreceptors.
%
% 7/19/13  dhb  Wrote it.
% 8/12/13  dhb  Code more generally and get rid of some special cases.
%          dhb  For cmf-like spectral functions, print out peak wavelengths and peak values.

theFields = fields(photoreceptors);
for ii = 1:length(theFields);
    theField = theFields{ii};
    switch (theField)
        case 'species'
            % Print out species that parameters are extracted for
            fprintf('  * Photoreceptor species: ''%s''\n',photoreceptors.species);
            
        case 'types'
            % List the names of the photoreceptor types
            fprintf('  * Photoreceptor types:\n');
            for jj = 1:length(photoreceptors.types)
                fprintf('    * %s\n',photoreceptors.types{jj});
            end
            
        case 'nomogram'
            % Have to special case this one
            fprintf('  * Photoreceptor field %s\n',theField);
            fprintf('    * Wavelength sampling: %d nm start, %d nm step, %d samples\n', ...
                photoreceptors.nomogram.S(1), photoreceptors.nomogram.S(2), photoreceptors.nomogram.S(3));
            if (isfield(photoreceptors.nomogram,'source') && ~strcmp(photoreceptors.nomogram.source,'None'))
                eval(['theNumberTypes = length(photoreceptors.' theField '.lambdaMax);']);
                eval(['theSource = photoreceptors.' theField '.source;']);
                fprintf('    * Source: ''%s'', value for each photoreceptor type: ',theSource);
                for jj = 1:theNumberTypes
                    fprintf('%g nm ',eval(['photoreceptors.' theField '.lambdaMax(jj)']));
                end
                fprintf('\n');
            end
            
        case {'absorbance' 'absorbtance' 'effectiveAbsorbtance' 'isomerizationAbsorbtance' 'energyFundamentals' 'quantalFundamentals'}
            eval(['theCmf = photoreceptors.' theField ';']);
            [peakWls, peakVals] = FindCmfPeaks(photoreceptors.nomogram.S,theCmf);
            fprintf('  * Photoreceptors field %s\n',theField);
            fprintf('    * Spectral peaks at:');
            for ii = 1:length(peakVals)
                fprintf(' %d',peakWls(ii));
            end
            fprintf('\n');
            fprintf('    * Values at peaks:');
            for ii = 1:length(peakVals)
                fprintf(' %0.4f',peakVals(ii));
            end
            fprintf('\n');
            
        case {'lensDensity' 'macularPigmentDensity' 'preReceptoral'}
            % Print just source for these fields
            fprintf('  * Photoreceptors field %s\n',theField);
            hasSource1 = eval(['isfield(photoreceptors.' theField ',''source'');']);
            if (hasSource1)
                hasSource2 = eval(['~isempty(photoreceptors.' theField '.source);']);
                if (hasSource2)
                    eval(['theSource = photoreceptors.' theField '.source;']);
                    fprintf('    * Source: ''%s''\n',theSource);
                    hasValue0 = eval(['~strcmp(photoreceptors.' theField '.source,''None'');']);
                end
            end
            hasTrans1 = eval(['isfield(photoreceptors.' theField ',''transmittance'');']);
            if (hasTrans1)
                hasTrans2 = eval(['~isempty(photoreceptors.' theField '.transmittance);']);
                if (hasSource2)
                    eval(['theCmf = photoreceptors.' theField '.transmittance;']);
                    [peakWls, peakVals] = FindCmfPeaks(photoreceptors.nomogram.S,theCmf);
                    fprintf('    * Spectral peaks at:');
                    for ii = 1:length(peakVals)
                        fprintf(' %d',peakWls(ii));
                    end
                    fprintf('\n');
                    fprintf('    * Values at peaks:');
                    for ii = 1:length(peakVals)
                        fprintf(' %0.4f',peakVals(ii));
                    end
                    fprintf('\n');
                end
            end
            
        case {'ageInYears' 'fieldSizeDegrees'}
            % Just a numeric field, but don't print value if is empty
            eval(['theValue = photoreceptors.' theField ';']);
            if (~isempty(theValue))
                fprintf('  * Photoreceptors field %s: %g\n',theField,theValue);
            end
            
        otherwise
            % Other theFields are source/value pairs, print generically.
            % Sometimes one of the fields is not there or empty, in which
            % case just print the other.
            fprintf('  * Photoreceptors field %s\n',theField);
            hasSource1 = eval(['isfield(photoreceptors.' theField ',''source'');']);
            if (hasSource1)
                hasSource2 = eval(['~isempty(photoreceptors.' theField '.source);']);
                if (hasSource2)
                    eval(['theSource = photoreceptors.' theField '.source;']);
                    fprintf('    * Source: ''%s''\n',theSource);
                    hasValue0 = eval(['~strcmp(photoreceptors.' theField '.source,''None'');']);
                else
                    hasValue0 = true;
                end
            else
                hasValue0 = true;
            end
            
            hasValue1 = eval(['isfield(photoreceptors.' theField ',''value'');']);
            if (hasValue1)
                hasValue2 = eval(['~isempty(photoreceptors.' theField '.value);']);
                if (hasValue2 & hasValue0)
                    eval(['theValue = photoreceptors.' theField '.value;']);
                    fprintf('    * Value:');
                    valDim = length(theValue);
                    for jj = 1:valDim
                        fprintf(' %0.4g',theValue(jj));
                    end
                    fprintf('\n');
                end
            end
    end
end
