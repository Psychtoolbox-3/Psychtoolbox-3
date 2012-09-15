function ResolutionTest
% ResolutionTest prints screen resolutions reported by Screen 'Resolution' 
% and 'Resolutions'. 
% 
% Also see SetResolution, NearestResolution, and Screen Resolution and Resolutions.
% 
% Denis Pelli

% 1/27/00 dgp Wrote it.
% 9/17/01 dgp Added "recommended" to the report.
% 4/24/02 awi Exit on PC with message.

% Test every screen
for screenNumber=Screen('Screens')
	% Describe video card
	% DescribeScreen(screenNumber);

	fprintf('\nCURRENT RESOLUTION:\n');
	res=Screen(screenNumber,'Resolution');
	disp(res);
	
	if 0
		res=screen(screenNumber,'Resolutions');
		fprintf('mode %-29s valid safe rec recNow default interlaced\n','resolution');
		for i=1:length(res)
			fprintf('%3d  %-29s  %-5d %-4d %-3d %-6d %-7d %-10d\n',res(i).mode,res(i).name,res(i).valid,res(i).safe,res(i).recommended,res(i).recommendNow,res(i).default,res(i).interlaced);
		end
		fprintf('\n')
	end
	
	if 0
		res=screen(screenNumber,'Resolutions');
		fprintf('mode %-21s pixelSize\n','resolution');
		for i=1:length(res)
			fprintf('%3d ',res(i).mode);
		%	fprintf('%-29s  ',res(i).name);
			fprintf('%4d x %4d   %3.0f Hz   ',res(i).width,res(i).height,res(i).hz);
			for j=1:length(res(i).pixelSizes)
				fprintf('%d ',res(i).pixelSizes(j));
			end
			fprintf('bits\n');
		end
		fprintf('\n')
	end
	
	fprintf('\nSCREEN %i: AVAILABLE RESOLUTIONS:\n', screenNumber);
	% fprintf('mode %-29s         pixelSize\n','resolution');
	res=Screen(screenNumber,'Resolutions');
    oldres = '';
    
	for i=1:length(res)
        resname = sprintf('%dx%d ', res(i).width, res(i).height);
        if isempty(strfind(oldres, resname))
            oldres = [oldres resname];
            fprintf('%d x %d\n', res(i).width, res(i).height);
        end
        
        %fprintf('%3d  ',res(i).mode);
		%fprintf('%-29s  ',res(i).name);
        % 		s=sprintf('%d',res(i).pixelSizes(1));;
        % 		for j=2:length(res(i).pixelSizes)
        % 			s=sprintf('%s,%d',s,res(i).pixelSizes(j));
        % 		end
        %		fprintf('%14s bits   ',s);
        % 		if res(i).safe
        % 			fprintf('safe  ');
        % 		else
        % 			fprintf('      ');
        % 		end
        % 		if res(i).recommended
        % 			fprintf('recommended  ');
        % 		else
        % 			fprintf('             ');
        % 		end
		%       fprintf('\n');
	end
end
% DescribeScreen(-1); % end of table
