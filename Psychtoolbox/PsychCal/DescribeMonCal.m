function same = DescribeMonCal(cal,file,whichScreen)
% same = DescribeMonCal(cal,[file],[whichScreen])
% 
% Print descriptive information about a calibration 
% to the command window or file.
%
% Argument file is a standard Matlab file descriptor,
% see fopen.  If file arg is omitted or empty, printout
% goes to command window.
%
% If argument whichScreen is passed, a description of
% the current hardware is also printed.  In this case,
% returned boolean same indicates whether the calibration
% is consistent with the current hardware.  Boolean
% same is empty if whichScreen is not provided.
%
% 8/25/97  dhb, pbe  Wrote it.
% 7/3/98   dhb, pbe  Updated for cal.describe.
% 12/3/99  dhb, mpr  Fix check for calibration desription field.
% 8/1800   dhb       Add whichScreen arg, same return.
% 6/29/02  dgp       Use new version of Screen VideoCard.
% 9/23/02  dhb, jms  Fix small bug in way driver is compared, presumably introduced 6/29/02.
% 9/29/08  dhb, tyl, ijk Update for OS/X, current computer stuff.
%                    Comparison of computer name skipped, because it seems to vary with login. 
% 6/24/11  dhb       Dump out gamma fit type and exponents if gamma function was fit with a simple power function.
% 5/28/13  dhb       Change output printed format to make it easier to paste into Doku wiki.

% Default args
if (nargin < 2 || isempty(file))
	file = 1;
end
if (nargin < 3 || isempty(whichScreen))
	file = 1;
	whichScreen = [];
end
same = [];

if (~isfield(cal,'describe'))
	error('Calibration structure has no description');
end

fprintf('Calibration:\n');
fprintf(file,'  * Computer: %s\n',cal.describe.computer);
fprintf(file,'  * Screen: %d\n',cal.describe.whichScreen);
fprintf(file,'  * Monitor: %s\n',cal.describe.monitor);
fprintf(file,'  * Video driver: %s\n',cal.describe.driver);
fprintf(file,'  * Dac size: %g\n',cal.describe.dacsize);
fprintf(file,'  * Frame rate: %g hz\n',cal.describe.hz);
fprintf(file,'  * Calibration performed by %s\n',cal.describe.who);
fprintf(file,'  * Calibration performed on %s\n',cal.describe.date);
fprintf(file,'  * Calibration program: %s\n',cal.describe.program);
fprintf(file,'  * Comment: %s\n',cal.describe.comment);
fprintf(file,'  * Calibrated device has %g primaries\n',cal.nDevices);
fprintf(file,'  * Gamma fit type %s\n',cal.describe.gamma.fitType);
if (strcmp(cal.describe.gamma.fitType,'simplePower'))
    fprintf(file,'  * Simple power gamma exponents are: %0.2f, %0.2f, %0.2f\n',...
        cal.describe.gamma.exponents(1),cal.describe.gamma.exponents(2),cal.describe.gamma.exponents(3));
end
fprintf(file,'\n');

% Current configuration
if (~isempty(whichScreen))
    cal.describe.driver = sprintf('%s %s','unknown_driver','unknown_driver_version');
    computerInfo = Screen('Computer');
    computer = sprintf('%s''s %s, %s', computerInfo.consoleUserName, computerInfo.machineName, computerInfo.system);
    driver = sprintf('%s %s','unknown_driver','unknown_driver_version');
    dacsize = ScreenDacBits(whichScreen);
    hz = Screen('NominalFrameRate',whichScreen);
    same = 1;
    fprintf('Current configuration:\n');
    fprintf(file,'  * Computer: %s\n',computer);
    if (~streq(computer,cal.describe.computer))
        %same = 0;
    end
    fprintf(file,'  * Screen: %d\n',whichScreen);
    if (whichScreen ~= cal.describe.whichScreen)
        save = 0;
    end
    fprintf(file,'  * Video driver: %s\n',driver);
    if (~streq(driver,cal.describe.driver))
        same = 0;
    end
    fprintf(file,'  * Dac size: %g\n',dacsize);
    if (dacsize ~= cal.describe.dacsize)
        same = 0;
    end
    fprintf(file,'  * Frame rate: %g hz\n',hz);
    if (abs(hz-cal.describe.hz) > 0.5)
        same = 0;
    end
end
