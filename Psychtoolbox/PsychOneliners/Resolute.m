function newplaat = Resolute(plaat,newres)
% newplaat = Resolute(plaat,newres)
% accepts matrix PLAAT and resizes it to NEWRES
% resizing is done by cutting and/or adding a white background
%
% JvR & DN 2008       Wrote it
% DN       2008-07-30 Simplified

res             = size(plaat);

maxres          = max([res(1:2); newres]);
maxplaat        = 255*ones([maxres 3]);                             % create white canvas the size of the output

center          = CenterRect([1 1 res(1:2)],[1 1 maxres]);          % PsychToolBox function
newcenter       = CenterRect([1 1 newres]  ,[1 1 maxres]);

maxplaat(center(1):center(3),center(2):center(4),:)...
                = plaat;

newplaat        = maxplaat(newcenter(1):newcenter(3),newcenter(2):newcenter(4),:);
