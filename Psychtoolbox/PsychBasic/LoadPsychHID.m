function LoadPsychHID
% LoadPsychHID - Try to get PsychHID linked and loaded on MS-Windows, no
% matter what.

if ~IsWin
    % Nothing to do on non-Windows:
    return;
end

% Windows on Octave or with Matlab R2007a or later. Give it a try:
try
    PsychHID('Version');
catch %#ok<CTCH>
    % PsychHID loading and linking failed:
    if IsWin
        % PsychHID failed on Windows. Most likely cause would be "invalid
        % MEX file error" due to PsychHID failing to link against
        % required DLL's. libusb-1.0 may not be properly installed?
        fprintf('INFO: Initial invocation of the PsychHID mex file failed.\n');
        fprintf('INFO: Most likely a required DLL is not installed in your system, e.g., libusb-1.0.dll\n');
        fprintf('INFO: I will now check if this is the culprit and work around it. To avoid future warnings,\n');
        if IsWin(1)
            % 64-Bit installation:
            fprintf('INFO: please copy the 64-Bit libusb-1.0.dll from the PsychContributed\\x64 folder into your C:\\WINDOWS\\system32\\ \n');
        else
            % 32-Bit installation:
            fprintf('INFO: please copy the 32-Bit libusb-1.0.dll from the PsychContributed folder into your C:\\WINDOWS\\system32\\ \n');
            fprintf('INFO: folder or - on a 64-Bit Windows setup into the C:\\WINDOWS\\SysWOW64\\ - \n');
        end
        fprintf('INFO: folder or a similarly appropriate place. You can get a fresher copy of libusb-1.0.dll from \n');
        fprintf('INFO: http://libusb.org/wiki/windows_backend if you want to stay up to date.\n');
        fprintf('INFO: Retrying now, may fail...\n');

        % The old drill: cd into our PsychContributed folder which
        % contains required DLL's. Retry by self-invocation. If this was
        % the culprit, then the linker should load, link and init
        % PsychHID and we should succeed. Otherwise we fail again.
        wd = pwd;
        try
            if IsWin(1)
                % 64-Bit version of libusb.dll
                cd([PsychtoolboxRoot 'PsychContributed' filesep 'x64' filesep]);
            else
                % 32-Bit version of libusb.dll
                cd([PsychtoolboxRoot 'PsychContributed' filesep ]);
            end
            PsychHID('Version');
            cd(wd);
            return;
        catch %#ok<CTCH>
            cd(wd);
            psychrethrow(psychlasterror);
        end
    end

    % Game over.
    psychrethrow(psychlasterror);
end
