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
% 10/7/12 mk  Clean it up. Output refresh rates an pixelsize as well.

% Test every screen
for screenNumber=Screen('Screens')
    % Describe video card
    % DescribeScreen(screenNumber);
    
    fprintf('\nCURRENT RESOLUTION:\n');
    res=Screen(screenNumber,'Resolution');
    disp(res);
    
    fprintf('\nSCREEN %i: AVAILABLE RESOLUTIONS:\n', screenNumber);
    res=Screen(screenNumber,'Resolutions');
    oldres = '';
    
    for i=1:length(res)
        resname = sprintf('%dx%d ', res(i).width, res(i).height);
        if isempty(strfind(oldres, resname))
            oldres = [oldres resname]; %#ok<*AGROW>
            fprintf('%d x %d\n', res(i).width, res(i).height);
        end
    end
    
    fprintf('\nSCREEN %i: AVAILABLE DETAILED RESOLUTIONS:\n', screenNumber);
    res=Screen(screenNumber,'Resolutions');
    for i=1:length(res)
        fprintf('%4d x %4d   %3.0f Hz %d ',res(i).width,res(i).height,res(i).hz, res(i).pixelSize);
        fprintf('bits\n');
    end
    fprintf('\n')
end
