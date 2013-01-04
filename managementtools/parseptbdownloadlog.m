function parseptbdownloadlog(fname)
% parseptbdownloadlog(fname) - Download log parser.
%
% Parses a Psychtoolbox online registration logfile
% 'fname', does stats on it, prints the summary.
%
% History:
% 16.10.2006 Written (MK).

if nargin < 1
    fname = 'ptbregistration.log';
end

fd=fopen(fname, 'rt');
if fd==-1
    error('File not found.');
end

ofl = [];
knownmacids = containers.Map;

transactioncount = 0;
totalcount = 0;
betacount = 0;
stablecount = 0;
trunkcount = 0;
unknowncount = 0;
oldptb306count = 0;
osxcount = 0;
wincount = 0;
linuxcount = 0;

try
    while(1)
        % Read next line from file:
        fl = fgetl(fd);
        if fl==-1
            % Stop at end of file:
            break;
        end
        
        % Concatenate it to current piece:
        ofl = [ ofl fl ];
        
        % Real end reached?
        le = strfind(ofl, '</DATE>');
        if isempty(le)
            % Nope, just go ahead.
        else
            % Yes! Process it, then restart.
            %fprintf('%s\n', ofl);
            
            % This counts as one transaction:
            transactioncount = transactioncount + 1;
            
            % We count each MACID only once:
            ls = strfind(ofl, '<MACID>');
            le = strfind(ofl, '</MACID>');
            
            % Extract macid as character string:
            macid = ofl(ls+7:le-1);
            
            if ~knownmacids.isKey(macid)
                % This is a new one. Count it:
                knownmacids(macid) = 1;
                
                % Total unique count:
                totalcount = totalcount + 1;
                
                if (mod(totalcount, 1000) == 0)
                    fprintf('Elapsed %i\n', totalcount);
                end
                
%                 if strfind(ofl, '<FLAVOR>beta</FLAVOR>')
%                     betacount = betacount + 1;
%                 end
%                 
%                 if strfind(ofl, '<FLAVOR>stable</FLAVOR>')
%                     stablecount = stablecount + 1;
%                 end
%                 
%                 if strfind(ofl, '<FLAVOR>trunk</FLAVOR>')
%                     trunkcount = trunkcount + 1;
%                 end
%                 
%                 if strfind(ofl, '<FLAVOR>Psychtoolbox-3.0.6</FLAVOR>')
%                     oldptb306count = oldptb306count + 1;
%                 end
%                 
%                 if strfind(ofl, '<FLAVOR>unknown</FLAVOR>')
%                     unknowncount = unknowncount + 1;
%                 end
                
                if strfind(ofl, '<OS>MacOS-X')
                    osxcount = osxcount + 1;
                end
                
                if strfind(ofl, '<OS>Windows')
                    wincount = wincount + 1;
                end
                
                if strfind(ofl, '<OS>Linux')
                    linuxcount = linuxcount + 1;
                end
                
            end
            
            % Reset for next item:
            ofl = [];
        end    
    end
    fclose(fd);
catch
    fclose(fd);
end

transactioncount
totalcount
betacount
stablecount
trunkcount
oldptb306count
unknowncount
osxcount
wincount
linuxcount
end
