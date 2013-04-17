function PsychLinuxConfiguration
% PsychLinuxConfiguration -- Optimize setup of Linux system.
%
% This script modifies system settings and configuration files
% to optimize a Linux system for use with Psychtoolbox.
%
% Currently it modifies files to allow to run Octave or Matlab
% as regular non-root user, ie. without need for root login or
% the "sudo" command. It does so by changing file permissions
% and resource usage limits to allow a regular user application
% to switch to realtime scheduling, lock its memory, and to
% access special purpose hardware like GPUs, Bits+, Datapixx and
% other research equipment.
%
% Realtime optimizations are achieved by extending the
% /etc/security/limits.conf file with entries that allow
% members of the Unix user group "psychtoolbox" to lock
% application memory into physical RAM, eliminating/minimizing
% interference from the VM subsystem, and to run with realtime
% priorities up to level 50. The group "psychtoolbox" is created
% if it does not already exist.
%
% If the target system has a /etc/security/limits.d/ directory,
% then a separate rule file is stored to that directory to
% achieve the change without messing around with the limits.conf
% file.
%
% root-less hardware access is achieved by copying a special
% psychtoolbox.rules file into the /etc/udev/rules.d/ directory.
% This udev rules file contains rules to auto-detect certain
% hardware at bootup or when the hw is hot-plugged and to
% reconfigure this hw or access permission for root-less access
% by Psychtoolbox, and for optimal performance for the kind
% of typical PTB use cases.
%
% The script calls into the shell via "sudo" to achieve this
% setup task, which itself needs admin privileges to modify
% system files etc. "sudo" will prompt the user for his admin
% password to complete the tasks.
%

% History:
%  6.01.2012   mk  Written.
% 16.04.2013   mk  Add info about 'dialout' group for serial port access.

if ~IsLinux
  return;
end

% Retrieve login username of current user:
[ignore, username] = system('whoami'); %#ok<ASGLU>
username = username(1:end-1);
addgroup = 0;

% Setup of /etc/udev/rules.d/psychtoolbox.rules file, if needed:

% Assume no need to install or update:
needinstall = 0;
fprintf('\n\nLinux specific system setup for running Psychtoolbox as non-root user:\n');
fprintf('You need to be a user with administrative rights for this function to succeed.\n');
fprintf('If you don''t have administrator rights, or if you don''t trust this script to\n');
fprintf('tinker around with system settings, simply answer all questions with "n" for "No"\n');
fprintf('and then call a system administrator for help.\n\n');
fprintf('Checking if the Psychtoolbox udev rules file is installed and up to date.\n');
fprintf('This file will allow Psychtoolbox to access special research hardware equipment,\n');
fprintf('e.g., the Cambridge Research Systems Bits+ box, the Datapixx from VPixx, response\n');
fprintf('button boxes, and some special features of your graphics card, e.g., high precision\n');
fprintf('timestamping. You will be able to access this hardware without the need to run\n');
fprintf('Matlab or Octave as sudo root user.\n\n');

% Check if udev psychtoolbox.rules file exists:
if ~exist('/etc/udev/rules.d/psychtoolbox.rules', 'file')
  % No: Needs to be installed.
  needinstall = 1;
  fprintf('The udev rules file for Psychtoolbox is not installed on your system.\n');
  answer = input('Should i install it? [y/n] : ', 's');
else
  % Yes.
  fprintf('The udev rules file for Psychtoolbox is already installed on your system.\n');

  % Compare its modification date with the one in PTB:
  r = dir([PsychtoolboxRoot '/PsychBasic/psychtoolbox.rules']);
  i = dir('/etc/udev/rules.d/psychtoolbox.rules');

  if r.datenum > i.datenum
    needinstall = 2;
    fprintf('However, it seems to be outdated. I have a more recent version with me.\n');
    answer = input('Should i update it? [y/n] : ', 's');
  end
end

if needinstall && answer == 'y'
  fprintf('I will copy my most recent rules file to your system. Please enter\n');
  fprintf('now your system administrator password. You will not see any feedback.\n');
  drawnow;

  cmd = sprintf('sudo cp %s/PsychBasic/psychtoolbox.rules /etc/udev/rules.d/', PsychtoolboxRoot);
  [rc, msg] = system(cmd);
  if rc == 0
    fprintf('Success! You may need to reboot your machine for some changes to take effect.\n');
  else
    fprintf('Failed! The error message was: %s\n', msg);
  end
end

% First the fallback implementation if /etc/security/limits.d/ does not
% exist:
if ~exist('/etc/security/limits.d/', 'dir')
% Check if /etc/security/limits.conf has proper entries to allow memory locking
% and real-time scheduling:
fid = fopen('/etc/security/limits.conf');
if fid == -1
  fprintf('Could not open /etc/security/limits.conf for reading. Can not set it up, sorry!\n');
  return;
end

fprintf('\nChecking if /etc/security/limits.conf has entries which allow everyone to\n');
fprintf('make use of realtime scheduling and memory locking -- Needed for good timing.\n\n');

mlockok = 0;
rtpriook = 0;

while ~feof(fid)
  fl = fgetl(fid);
  % fprintf('%s\n', fl);

  if fl == -1
    continue;
  end

  if ~isempty(strfind(fl, 'memlock')) && ~isempty(strfind(fl, 'unlimited')) && ~isempty(strfind(fl, '@psychtoolbox')) && ~isempty(strfind(fl, '-'))
    mlockok = 1;
  end

  if ~isempty(strfind(fl, 'rtprio')) && ~isempty(strfind(fl, '50')) && ~isempty(strfind(fl, '@psychtoolbox')) && ~isempty(strfind(fl, '-'))
    rtpriook = 1;
  end
end

% Done reading the file:
fclose(fid);

drawnow;

% ...and?
if ~(mlockok && rtpriook)
  fprintf('\n\nThe file seems to be missing some suitable setup lines.\n');
  answer = input('Should i add them for you? [y/n] : ', 's');
  if answer == 'y'
    fprintf('I will try to add config lines to your system. Please enter\n');
    fprintf('now your system administrator password. You will not see any feedback.\n');
    drawnow;

    % Set amount of lockable memory to unlimited for all users:
    [rc, msg] = system('sudo /bin/bash -c ''echo "@psychtoolbox     -     memlock     unlimited" >> /etc/security/limits.conf''');
    if rc ~= 0
      fprintf('Failed! The error message was: %s\n', msg);
    end

    % Set allowable realtime priority for all users to 50:
    [rc2, msg] = system('sudo /bin/bash -c ''echo "@psychtoolbox     -     rtprio      50" >> /etc/security/limits.conf''');
    if rc2 ~= 0
      fprintf('Failed! The error message was: %s\n', msg);
    end

    % Must add a psychtoolbox user group:
    addgroup = 1;

    if (rc == 0) && (rc2 == 0)
      fprintf('\n\nSuccess!\n\n');
    else
      fprintf('\n\nFailed! Maybe ask a system administrator for help?\n\n');
    end
  end
else
  fprintf('\n\nYour system is already setup for use of Priority().\n\n');
end

else
% Realtime setup for systems with /etc/security/limits.d/ directory:
% Simply install a ptb specific config file - a cleaner solution:
if ~exist('/etc/security/limits.d/99-psychtoolboxlimits.conf', 'file')
  fprintf('\n\nThe file /etc/security/limits.d/99-psychtoolboxlimits.conf is\n');
  fprintf('not yet installed on your system. It allows painless realtime operation.\n');
  answer = input('Should i install the file for you? [y/n] : ', 's');
  if answer == 'y'
    fprintf('I will try to install it now to your system. Please enter\n');
    fprintf('now your system administrator password. You will not see any feedback.\n');
    drawnow;
    cmd = sprintf('sudo cp %s/PsychBasic/99-psychtoolboxlimits.conf /etc/security/limits.d/', PsychtoolboxRoot);
    [rc, msg] = system(cmd);
    if rc ~= 0
      fprintf('Failed! The error message was: %s\n', msg);
    else
      fprintf('Success!\n\n');

      % Must add a psychtoolbox user group:
      addgroup = 1;
    end
  end
end
end

% Need to create a Unix user group 'psychtoolbox' and add user to it?
if addgroup
  % This will create the psychtoolbox user group, unless the group
  % already exists. In such a case it simply does nothing:
  system('sudo groupadd --force psychtoolbox');

  fprintf('I have created a new Unix user group called "psychtoolbox" on your system.\n');
else
  fprintf('\n\nYour system has a Unix user group called "psychtoolbox".\n');
end

fprintf('All members of that group can use the Priority() command now without the need\n');
fprintf('to run Matlab or Octave as sudo root user.\n\n');
fprintf('You need to add each user of Psychtoolbox to that group. You could do this\n');
fprintf('with the user management tools of your system. Or you can open a terminal window\n');
fprintf('and type the following command (here as an example to add yourself to that group):\n\n');
fprintf('sudo usermod -a -G psychtoolbox %s\n\n', username);
fprintf('One should also add oneself to the ''dialout'' group for access to serial port devices:\n\n');
fprintf('sudo usermod -a -G dialout %s\n\n', username);
fprintf('After that, the new group member must log out and then login again for the\n');
fprintf('settings to take effect.\n\n');
fprintf('\nFinished. Your system should now be ready for use with Psychtoolbox.\n');
fprintf('If you encounter problems, try rebooting the machine. Some of the settings only\n');
fprintf('become effective after a reboot.\n\n\n');
fprintf('Press any key to continue.\n');
pause;
fprintf('\n\n\n');

return;
