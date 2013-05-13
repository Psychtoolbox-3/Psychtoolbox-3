function AddPsychJavaPath

% Add Psychtoolbox directories containing Java classes to the path
% which MATLAB searches for Java classes.  The MATLAB Java path is separate
% from the path searched for .m and mex functions.
%
% PsychtJavaPaths updates the MATLAB Java path lazily, detecting first
% whether the MATLAB paths have already been changed before making changes.
%
% AddPsychJavaPath is called by Psychtoolbox functions which depend on
% Psychtoolbox Java classes; It should be unnecessary to use it within your
% own programs.    
%
% see also: PsychJava, IsPsychJavaPathSet

% HISTORY
% 6/20/2006  awi    Wrote it.  Uses function javaaddpath as suggested by Mario.


% We use IsPsychJavaPathSet to retain state between invocations because use
% of javaaddpath erases global and persistent state within the function
% from which it is called.  This is presumably bug in MATLAB (version
% 7.0.4.352 (R14) Service Pack 2).


if ~IsPsychJavaPathSet
    path_PsychtoolboxVersion= which('PsychtoolboxVersion');
    path_PsychBasic= fileparts(path_PsychtoolboxVersion);
    path_PsychtoolboxRoot= fileparts(path_PsychBasic);
    path_PsychJava= fullfile(path_PsychtoolboxRoot, 'PsychJava');
    %check to see if it is the static path before adding it to the dynamic
    %path
    staticPath=javaclasspath('-static');
    if ~any(strcmp(staticPath, path_PsychJava))
        javaaddpath(path_PsychJava)
    end
    IsPsychJavaPathSet(1);
end




% persistent PSYCH_ARE_JAVA_DIRECTORIES_IN_MATLAB_JAVA_PATH;
% 
% 
% if isempty(PSYCH_ARE_JAVA_DIRECTORIES_IN_MATLAB_JAVA_PATH)
%     path_PsychtoolboxVersion= which('PsychtoolboxVersion');
%     path_PsychBasic= fileparts(path_PsychtoolboxVersion);
%     path_PsychtoolboxRoot= fileparts(path_PsychBasic);
%     path_PsychJava= fullfile(path_PsychtoolboxRoot, 'PsychJava');
%     javaaddpath(path_PsychJava);
%     PSYCH_ARE_JAVA_DIRECTORIES_IN_MATLAB_JAVA_PATH=1;
% end
