function err=EyelinkLegacyCalTargetBeep(el)

% not sure yet about OSX sound routines
if 0
    err=SND('Play', el.targetdisplaysound);
    SND('Wait');
end
