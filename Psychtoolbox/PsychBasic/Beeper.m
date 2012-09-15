function Beeper(frequency, fVolume, durationSec);
% function Beeper(frequency, [fVolume], [durationSec]);
%
% Play a beep sound.  Default is 400 Hz for .15 sec.
% frequency can be a number, or else the string 'high', 'med', or 'low'.
%
% fVolume - normalized to range of 0 - 1.  Default is 0.4;  
% Warning:  1 is the maximum volume and is often very loud!
%
% Funny name is because Matlab 6 contains a built-in function called "beep".
%
% 2006-02-15 - cburns
%   -   Added fVolume param
%   -   Swapped parameter order
%
% 2006-02-02 - cburns
%   -   Scaled down the volume of the sound to match the system volume better.  It was at maximum volume.
%       Would scare you enough to rip the bite bar off it's mount!
%       And switched to useing the sound() function, instead of the soundsc() function
%       which, by default, maximizes the volume.
%   -   Update, using the PsychToolbox Snd function which should return immediately.
%       Were experiencing delays with sound function
%
% 12/2/00 Backus - original version

if ~exist('frequency', 'var')
  frequency = 400;
end

if ~exist('durationSec', 'var')
  durationSec = 0.15;
end

if ~exist('fVolume', 'var')
    fVolume = 0.4;
else
    % Clamp if necessary
    if (fVolume > 1.0)
        fVolume = 1.0;
    elseif (fVolume < 0)
        fVolume = 0;
    end
end

if ischar(frequency)
  if strcmp(lower(frequency), 'high') frequency = 1000; 
  elseif strcmp(lower(frequency), 'med') frequency = 400;
  elseif strcmp(lower(frequency), 'medium') frequency = 400;
  elseif strcmp(lower(frequency), 'low') frequency = 220;
  end
end

sampleRate = Snd('DefaultRate');

nSample = sampleRate*durationSec;
soundVec = sin(2*pi*frequency*(1:nSample)/sampleRate);

% Scale down the volume
soundVec = soundVec * fVolume;
% sound(soundVec);
try % this part sometimes crashes for unknown reasons. If it happens replace sound with normal beep
    
    Snd('Play', soundVec, sampleRate);
catch
    beep
end
