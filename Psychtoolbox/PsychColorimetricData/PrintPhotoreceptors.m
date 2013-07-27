function PrintPhotoreceptors(photoreceptors)
% PrintPhotoreceptors(photoreceptors)
%
% Print to command window an interpretable output
% of what is in a photoreceptors structure.
%
% See also DefaultPhotoreceptors, FillInPhotoreceptors.
%
% 7/19/13  dhb  Wrote it.

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
            
        case 'axialDensity'
            % Just have values for each photoreceptor type but no source, print these
            fprintf('  * Photoreceptor field %s values for each photoreceptor type\n    * ',theField);
            eval(['theNumberTypes = length(photoreceptors.' theField '.value);']);
            for jj = 1:theNumberTypes
                fprintf('%g  ',eval(['photoreceptors.' theField '.value(jj)']));
            end
            fprintf('\n');
            
        case {'quantalEfficiency' 'OSlength' 'ISdiameter' 'ISlength' 'nomogram'}
            % Have values for each photoreceptor type but no source, print these
            fprintf('  * Photoreceptor field %s\n',theField);
            switch(theField)
                case 'nomogram'
                    fprintf('    * Wavelength sampling: %d nm start, %d nm step, %d samples\n', ...
                        photoreceptors.nomogram.S(1), photoreceptors.nomogram.S(2), photoreceptors.nomogram.S(3));
                    eval(['theNumberTypes = length(photoreceptors.' theField '.lambdaMax);']);
                otherwise
                    eval(['theNumberTypes = length(photoreceptors.' theField '.value);']);
            end
            eval(['theSource = photoreceptors.' theField '.source;']);
            fprintf('    * Source: ''%s'', value for each photoreceptor type: ',theSource);
            for jj = 1:theNumberTypes
                switch(theField)
                    case 'nomogram'
                        fprintf('%g nm ',eval(['photoreceptors.' theField '.lambdaMax(jj)']));
                    otherwise
                        fprintf('%g  ',eval(['photoreceptors.' theField '.value(jj)']));
                end
            end
            fprintf('\n');
            
        case {'absorbance' 'absorbtance' 'effectiveAbsorbtance' 'isomerizationAbsorbtance' 'preReceptoral'}
            % These fields are spectra, don't print anything
            
            %         case {'lensDensity' 'macularPigmentDensity'}
            %             % Print just source for these fields
            %             eval(['theSource = photoreceptors.' theField '.source;']);
            %             fprintf('  * Photoreceptors field %s\n',theField);
            %             fprintf('    * Source: ''%s''\n',theSource);
            
        case {'ageInYears' 'fieldSizeDegrees'}
            % Just a numeric feild
            eval(['theValue = photoreceptors.' theField ';']);
            fprintf('  * Photoreceptors field %s: \n',theField,theValue);
            
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
                end
            end
            hasValue1 = eval(['isfield(photoreceptors.' theField ',''value'');']);
            if (hasValue1)
                hasValue2 = eval(['~isempty(photoreceptors.' theField '.value);']);
                if (hasValue2)
                    eval(['theValue = photoreceptors.' theField '.value;']);
                    fprintf('    * Value: %0.4g\n',theValue);
                end
            end
    end
end
