function rc = IsMinimumOSXVersion(major, minor, point)
% rc = IsMinimumOSXVersion(major, minor, point);
%
% Checks if the script is running on a MacOS/X system with at least the
% requested (major,minor,point) version, e.g., to test for a MacOS/X system
% of 10.4.8 or later, do a IsMinimumOSXVersion(10,4,8);
%
% rc is 0 if the system is of a lower version, 1 if it satisfies this
% minimum version, 2 if the function doesn't know.

% History:
% 11/26/2007 Written (MK), based on code by Roger Woods (UCLA).
% 05/28/2012 Fixed this bullshit (MK).

if ~IsOSX
    rc = 0;
    return;
end

gestaltbits=Gestalt('sys1');
majorversion=bin2dec(char(49*gestaltbits+48*~gestaltbits));

% Preinit to 'No':
rc = 0;

if majorversion >= major
    if majorversion > major
        rc = 1;
        return;
    end
    
    gestaltbits=Gestalt('sys2');
    minorversion=bin2dec(char(49*gestaltbits+48*~gestaltbits));
    if minorversion >= minor
        if minorversion > minor
            rc = 1;
            return;
        end

        gestaltbits = Gestalt('sys3');
        pointversion=bin2dec(char(49*gestaltbits+48*~gestaltbits));
        if pointversion >= point
            rc = 1;
        end
    end
end

return;
