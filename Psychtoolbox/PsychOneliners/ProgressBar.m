function f = ProgressBar(nMax,str)
%ProgressBar    Ascii progress bar.
%   progBar = ProgressBar(nMax,str) creates a progress bar and returns a
%   pointer to a function handle which can then be called to update it.
%
%   To update, call progBar(currentStep)
%
%   Example:
%      n = 500;
%      progBar = ProgressBar(n,'computing...');
%      for tmp = 1:n
%        progBar(tmp);
%        pause(.01)
%      end

%   by David Szotten 2008
%   $Revision: 1.2 $  $Date: 2008/04/17 09:15:32 $
%   merged with utility by us / CSSM by: DN 2008
%  2008-09-16  DN Added elapsed time and estimated time left

if IsOctave
    % Octave either doesn't do nested functions, or doesn't have persistent
    % state within nested functions, at least in version 3.2.4
    error('ProgressBar is not supported on Octave');
end

if nargin>1
    head = sprintf('%s\n',str);
else
    head = '';
end

lastPercentileWritten = 0;
pstrlen = 0;
tstrlen = 0;


label = sprintf('| 0%s50%s100%% |\n',repmat(' ',1,21),repmat(' ',1,18));

fprintf('%s',[head label]);
hllen = length([head label]);

t=datenum(clock);

f = @updateBar;
    function updateBar(nCurrent)
        
        %what percentile are we up to
        currentPercentile = round(50*nCurrent/nMax);
        
        fprintf('%s',repmat(char(8),1,tstrlen)); % remove time string
        
        % compute time info
        ttn = datenum(clock)-t;
        tt  = datevec(ttn);
        dtt = ttn/nCurrent;
        ttleft = datevec(dtt*(nMax-nCurrent));
        tstr = sprintf('\nElapsed time:        %dh %dm %ds\nEstimated time left: %dh %dm %ds',tt(4),tt(5),round(tt(6)),ttleft(4),ttleft(5),round(ttleft(6)));
        tstrlen = length(tstr);

        %have we passed another percentile?
        if (currentPercentile > lastPercentileWritten )

            %we may have increased by several percentiles,
            %so keep writing until we catch up
            percentileToWrite = lastPercentileWritten + 1;
            while(percentileToWrite <= currentPercentile)

                %for every 10th, use a '+' instead
                if( mod(percentileToWrite,5)==0 )
                    fprintf('%s','+');
                else
                    fprintf('%s','.');
                end
                percentileToWrite = percentileToWrite + 1;
                pstrlen = pstrlen + 1;
            end

            %update status
            lastPercentileWritten = currentPercentile;

            % write time string
            fprintf('%s',tstr);
        else
            % write time string
            fprintf('%s',tstr);
        end
        

        %are we done?
        if nCurrent==nMax
            %clear bar
            pause(1)
            fprintf('%s',repmat(char(8),1,hllen+pstrlen+tstrlen));
        end
    end

end

