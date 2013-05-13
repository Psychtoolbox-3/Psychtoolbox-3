function isSet=IsPsychJavaPathSet(newValue)

% isSet=IsPsychJavaPathSet(newValue)
%
% Set and retain state indicating whether Psychtoolbox directories holding
% Java classes have been added to the MATLAB Java class path. 
%
% see also: AddPsychJavaPath, PsychJava

% HISTORY
% 6/20/06  awi  Wrote it;  Hack around a bug in MATLAB which erases
%                persistent and global state in functions which call 
%                javaaddpath.

persistent PSYCH_ISPSYCHJAVAPATHSET

if nargin==1
    if ~newValue
        newValue=[];
    end
    PSYCH_ISPSYCHJAVAPATHSET=newValue;
end
    
 
isSet= ~isempty(PSYCH_ISPSYCHJAVAPATHSET);
    
