function rc = IsGUI
% IsGUI - Is the Matlab or Octave GUI enabled in this session?
%
% Returns 1 if GUI is enabled, 0 otherwise.
%
% Matlab's GUI is based on JAVA. Octave's GUI is based on the QT toolkit.
%

if IsOctave && exist('isguirunning', 'builtin')
    rc = eval('isguirunning');
else
    rc = psychusejava('desktop');
end

return;
