function priorityStruct = MachGetPriorityMex(policyFlavorString, defaultFlag)
% priorityStruct = MachGetPriorityMex(policyFlavorString, defaultFlag)
%
% OSX: ___________________________________________________________________
%
% Retrieve current or default parameters for the main MATLAB thread for any
% of the  three priority flavors: 'THREAD_STANDARD_POLICY',
% 'THREAD_TIME_CONSTRAINT_POLICY'  and 'THREAD_PRECEDENCE_POLICY'
%
% The returned struct contains these fields:
% 
%   priorityStruct.threadID        
%       a number identifying the current thread. 
% 
%   priorityStruct.flavor          
%       one of: 'THREAD_STANDARD_POLICY', 'THREAD_TIME_CONSTRAINT_POLICY',
%       'THREAD_PRECEDENCE_POLICY'
% 
%   priorityStruct.policy           
%       see policy struct variants below
% 
%   priorityStruct.policySize      
%       The amount of memory allocated for the policy struct passed to the
%       Mach function thread_policy_get() by MachGetPriorityMex.
% 
%   priorityStruct.policyFillSize  
%       amount of memory filled into the  policy struct by
%       Mach function thread_policy_get().
% 
%   priorityStruct.getDefault      
%       value of the defaultFlag argument passed into MachGetPriorityMex. 
% 
%   priorityStruct.isDefault       
%       If the flag value passed to MachGetPriorityMex is 0, requesting
%       current parameter values and not default parameter values and yet
%       MachGetPriorityMex returns default parameters and defaultFlag value
%       1, then the priority flavor which was specified in the first
%       argument to priorityFlavorString is not in effect.  
% 
% The form of the embedded struct "policy" depends on the value of
% priorityFlavorString argument.
%  
% 'THREAD_STANDARD_POLICY':
%     priorityStruct.flavorPolicy.no_data    % a place holder only
%                     
% 'THREAD_TIME_CONSTRAINT_POLICY'     % see help MachSetTimeConstraintPriority
%     priorityStruct.flavorPolicy.period
%     priorityStruct.flavorPolicy.computation
%     priorityStruct.flavorPolicy.constraint
%     priorityStruct.flavorPolicy.preemptible
%     
% 'THREAD_PRECEDENCE_POLICY'          % see help MachSetTimeConstraintPriority
%     priorityStruct.flavorPolicy.importance
%
% There are three policy flavors but a thread may have only one of two 
% policy modes: THREAD_STANDARD_POLICY or THREAD_TIME_CONSTRAINT_POLICY.
% These are mutually-exclusive modes; setting a thread to either one will 
% unset the other mode.  The "importance" parameter associated with 
% THREAD_PRECEDENCE_POLICY is preserved after either THREAD_STANDARD_POLICY
% or THREAD_TIME_CONSTRAINT_POLICY is set, however the
% THREAD_PRECEDENCE_POLICY "importance" setting is ignored by the Mach task
% scheduler while a thread is in  THREAD_TIME_CONSTRAINT_POLICY mode.  A
% thread is governed by the "importance" parameter only when in
% THREAD_STANDARD_POLICY mode. 
%
% MachGetPriorityMex uses the OS X Darwin function thread_policy_get().
% For more information on thread_policy_get() see:
% Psychtoolbox3/Source/Common/MachPriorityMex/MachGetPriorityMex.c
% http://developer.apple.com/documentation/Darwin/Conceptual/KernelProgramming/scheduler/chapter_8_section_4.html
% /usr/include/mach/thread_policy.h
%
% OS9: ___________________________________________________________________
%
% MachGetPriorityMex does not exist in OS 9. 
%
% WIN: ________________________________________________________________
% 
% MachGetPriorityMex does not exist in Windows.
% 
% _________________________________________________________________________
%
% see also: Priority, Rush, MachGetPriorityFlavor, MachSetPriorityMex, 

% AUTHORS:
% Allen Ingling     awi     Allen.Ingling@nyu.edu


% NOTES:
% The "Mex" in "MachGetPriorityMex" indicates that that the compiled mex
% binary uses the Mex API instead of of the PsychScriptingGlue API to
% communicated with MATLAB.  This was done to make it more readablle  to
% Mathworks developers and others looking at the timing bug, but that  has
% been explained and MachSetPriority could be re written to use
% PsychScriptingGlue.


% HISTORY: 
% 8/13/03   awi     Wrote it.
% 7/16/04   awi     Cosmetic. Added PsychAssertMex call
% 3/12/05   dgp     Cosmetic.
% 10/10/05  awi     Noted dgp change on 3/12/05 

PsychAssertMex('OSX');


 
