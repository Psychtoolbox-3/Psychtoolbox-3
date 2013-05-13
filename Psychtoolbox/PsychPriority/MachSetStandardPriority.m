function MachSetStandardPriority

% MachSetStandardPriority
%
% OS X: ___________________________________________________________________
%
% Restores the standard priorty flavor THREAD_STANDARD_POLICY to  the main
% MATLAB thread.  MachSetStandardPriority undoes previous calls to either
% MachSetTimConstraintPolicy or 
% MachSetPriorityMex('THREAD_TIME_CONSTRAINT_POLICY',...)
%
% OS 9: ___________________________________________________________________
%
% MachSetStandardPriority does not exist in OS 9. 
%
% WINDOWS: ________________________________________________________________
% 
% MachSetStandardPriority does not exist in Windows.
% 
% _________________________________________________________________________
%
% see also: Priority, Rush, MachGetPriorityFlavor, MachSetPriorityMex, 


% AUTHORS:
% Allen Ingling     awi     Allen.Ingling@nyu.edu

% HISTORY: 
% 8/13/03   awi     Wrote it.
% 7/16/04   awi     added AssertOSX, cosmetic.

AssertOSX;
MachSetPriorityMex('THREAD_STANDARD_POLICY');
