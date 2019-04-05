/*
	PsychToolbox3/Source/Common/MachPriority/MachPriority.c
	
	PLATFORMS:	All
			
	
	AUTHORS:
	Allen Ingling		awi		Allen.Ingling@nyu.edu

	HISTORY:
		3/31/03		awi		wrote it.
        4/09/02		awi		wrote it using mach thread API, previous attemps using setpriority() and pthread priorities failed. 
        
	DESCRIPTION:
		Read and set  priority.   This function has been superceded by MachGetPriorityMex and MachSetPriorityMex which do better jobs of reading
		and setting priority but do not use the PsychtScriptingGlue layer.  The mex/mx API was used instead to make these functions more transparent 
		to Mathworks programmers doing work on the timing bug.  Because we understand that now, MachGetPriority and MachSetPriority could now be re-written
		to use PsychScriptingGlue.

	TO DO:
		¥Because we don't yet have abstracted functions for reading in struct arrays we accept 
		the inputs as seperate arguments.  This should be changed to also permit a struct array.
		¥Cache the original priority and restore it when this module is purged.
		
	NOTES:
		
	

*/


#include "MachPriority.h"

/*
    GetFlavorStringFromFlavorConstant()
    
*/
char *GetFlavorStringFromFlavorConstant(thread_policy_flavor_t flavor)
{
    static thread_policy_flavor_t flavorConstants[] = {THREAD_STANDARD_POLICY, THREAD_TIME_CONSTRAINT_POLICY, THREAD_PRECEDENCE_POLICY};
    static char *flavorStrings[] = {"THREAD_STANDARD_POLICY", "THREAD_TIME_CONSTRAINT_POLICY", "THREAD_PRECEDENCE_POLICY"};
    int numFlavors=3, i;
    
    for(i=0;i<=numFlavors;i++){
        if(i==numFlavors)
            return(NULL);
        if(flavorConstants[i]==flavor)
            return(flavorStrings[i]);
    }
    return(NULL);  //make the compiler happy
}


/*
struct thread_time_constraint_policy {
	uint32_t		period;
	uint32_t		computation;
	uint32_t		constraint;
	boolean_t		preemptible;
};
*/

PsychError MACHPRIORITYMachPriority(void) 
{
    //double 			*returnValue;
    //int			newPriority;
    const char 			*outerStructFieldNames[]={"thread", "flavor", "policy"};
    const char 			*policyStructFieldNames[]={"period", "computation", "constraint", "preemptible"};
    int 			numOuterStructDimensions=-1, numOuterStructFieldNames=3, numPolicyStructDimensions=-1, numPolicyStructFieldNames=4;
    PsychGenericScriptType	*outerStructArray, *policyStructArray;
    int 			kernError, numInputArgs;
    thread_act_t 		threadID;
    static thread_policy_flavor_t	currentThreadFlavor=THREAD_STANDARD_POLICY;
    struct			thread_time_constraint_policy	oldPolicyInfo, newPolicyInfo;
    mach_msg_type_number_t	msgTypeNumber;
    boolean_t			getDefault;
    char 			*flavorStr;
    psych_bool			setNewMode;
    double			*periodArg, *computationArg, *constraintArg, *preemptibleArg;
    char 			errorMessageStr[256];
    	 



    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(4)); //actually we permit only zero or three arguments.
    numInputArgs=PsychGetNumInputArgs();
    if(numInputArgs==4)
        setNewMode=TRUE;
    else if(numInputArgs==0)
        setNewMode=FALSE;
    else
        PsychErrorExitMsg(PsychError_user,"Incorrect number of arguments.  Either zero or four arguments accepted");
        

    //read the current settings
    threadID= mach_thread_self();
    currentThreadFlavor=THREAD_TIME_CONSTRAINT_POLICY;
    msgTypeNumber=THREAD_TIME_CONSTRAINT_POLICY_COUNT;
    getDefault=FALSE;
    kernError= thread_policy_get(threadID, currentThreadFlavor, (int *)&oldPolicyInfo, &msgTypeNumber, &getDefault);
    if(kernError != KERN_SUCCESS)
        PsychErrorExitMsg(PsychError_internal,"\"thread_policy_get()\" returned and error when reading current thread policy");
    
    

    //fill in the outgoig struct with current thread settings values.
    //outer struct
    PsychAllocOutStructArray(1, FALSE, numOuterStructDimensions, numOuterStructFieldNames, outerStructFieldNames, &outerStructArray);
    PsychSetStructArrayDoubleElement("thread", 0, (double)threadID, outerStructArray);
    flavorStr=GetFlavorStringFromFlavorConstant(currentThreadFlavor);
    PsychSetStructArrayStringElement("flavor", 0, flavorStr, outerStructArray);
    //enclosed policy struct
    PsychAllocOutStructArray(-1, FALSE, numPolicyStructDimensions, numPolicyStructFieldNames, policyStructFieldNames, &policyStructArray);
    PsychSetStructArrayDoubleElement("period", 0, (double)oldPolicyInfo.period, policyStructArray);
    PsychSetStructArrayDoubleElement("computation", 0, (double)oldPolicyInfo.computation, policyStructArray);
    PsychSetStructArrayDoubleElement("constraint", 0, (double)oldPolicyInfo.constraint, policyStructArray);
    PsychSetStructArrayDoubleElement("preemptible", 0, (double)oldPolicyInfo.preemptible, policyStructArray);
    PsychSetStructArrayStructElement("policy",0, policyStructArray, outerStructArray);
    
    //Set the priority 
    if(setNewMode){
        //if there is only one argument then it must be the flavor argument re
        PsychAllocInDoubleArg(1, TRUE, &periodArg);
        PsychAllocInDoubleArg(2, TRUE, &computationArg);
        PsychAllocInDoubleArg(3, TRUE, &constraintArg);
        PsychAllocInDoubleArg(4, TRUE, &preemptibleArg);
        newPolicyInfo.period=(uint32_t)*periodArg;
	newPolicyInfo.computation=(uint32_t)*computationArg;
	newPolicyInfo.constraint=(uint32_t)*constraintArg;
	newPolicyInfo.preemptible=(boolean_t)*preemptibleArg;

	kernError= thread_policy_set(threadID, THREAD_TIME_CONSTRAINT_POLICY, (int *)&newPolicyInfo, THREAD_TIME_CONSTRAINT_POLICY_COUNT);
	if(kernError != KERN_SUCCESS){
		sprintf(errorMessageStr,"%s%d", "\"thread_policy_set()\" returned and error when setting new thread policy: ", (int)kernError);
		PsychErrorExitMsg(PsychError_internal, errorMessageStr);
		}
    }

        
    return(PsychError_none);	
}




	
