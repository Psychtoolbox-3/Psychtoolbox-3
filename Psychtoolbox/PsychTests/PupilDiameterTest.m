% PupilDiameterTest
%
% Compare different formulae for computing
% pupil diamter from luminance.
%
% 4/2/99  dhb  Wrote it.
% 5/8/99  dhb  Update for new calling form.
% 4/25/02 awi,dgp  removed commands preventing plot 1 from appearing.  

% Calculate pupil size as a function of luminance
lums = logspace(-5,4,10);
[diamsPS,nil,trolands0] = PupilDiameterFromLum(lums,'Pokorny_Smith');
diamsMS = PupilDiameterFromLum(lums,'Moon_Spencer');
diamsDG = PupilDiameterFromLum(lums,'DeGroot_Gebhard');

% Plot them.
figure(1); clf;
plot(lums,diamsPS,'r');
hold on
plot(lums,diamsPS,'r+');
plot(lums,diamsMS,'g');
plot(lums,diamsMS,'g+');
plot(lums,diamsDG,'b');
plot(lums,diamsDG,'b+');
hold off
drawnow

% Compute effective trolands and compare with trolands
effTrolands1 = EffectiveTrolandsFromLum(lums,'Pokorny_Smith1');
effTrolands2 = EffectiveTrolandsFromLum(lums,'Pokorny_Smith2');

% Plots
figure(2); clf;
plot(lums,trolands0,'r');
hold on
plot(lums,effTrolands1,'g');
plot(lums,effTrolands2,'b');
hold off
drawnow
