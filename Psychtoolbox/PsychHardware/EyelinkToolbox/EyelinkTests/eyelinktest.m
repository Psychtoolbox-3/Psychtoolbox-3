clear all;
fprintf('My OSX Eyelink Test\n\n\t');
input('Hit the return key to continue.','s');
fprintf('Thanks.\n');


el=initeyelinkdefaultsOSX;

start=GetSecs;
while 1
    [key, el]=getkeyforeyelinkOSX(el);
    if key==el.TERMINATE_KEY
        break;
    end
    if GetSecs-start>1000
        fprintf('getkeyforeyelinkOSX timed-out.\n');
        break;
    end
    if key~=el.JUNK_KEY && key ~= 0
        fprintf('Eyelink key %s (%.1f)\n', key, GetSecs-start);
    end
end

FlushEvents;

fprintf('End of Eyelink GetKeys Test.\n');

