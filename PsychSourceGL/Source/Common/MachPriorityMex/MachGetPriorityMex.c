/*
	PsychToolbox3/Source/Common/MachPriorityMex/MachGetPriorityMex.c
	
	PLATFORMS:	OS X only.  
			
	
	AUTHORS:
	Allen Ingling	awi		Allen.Ingling@nyu.edu
    Mario Kleiner   mk      mario.kleiner@tuebingen.mpg.de
 
	HISTORY:
	8/08/02			awi		wrote it. 
        
        DESCRIPTION 
        
        MachGetPriorityMex is a simplified variant of MachPriority which strips out the Psychtoolbox mex/mx abstraction layer.
        For the sake of greater readability it uses a reduced argument set.  See also MachSetPriorityMex.
        
        
        MachGetPriorityMex(['THREAD_STANDARD_POLICY', 'THREAD_EXTENDED_POLICY','THREAD_TIME_CONSTRAINT_POLICY', 'THREAD_PRECEDENCE_POLICY'], getDefault);
        
        The returned argument is a struct whith these fields:
        
        priority.threadID
        priority.flavor
        priority.policy
        priority.policySize
        priority.policyFillSize
        priority.getDefault
        priority.isDefault
        
        The form of the embedded struct "policy" depends on the first argument passed to MachGetPriorityMex().
        
        THREAD_STANDARD_POLICY:
            priority.flavorPolicy.no_data
                        
        THREAD_TIME_CONSTRAINT_POLICY
            priority.flavorPolicy.period
            priority.flavorPolicy.computation
            priority.flavorPolicy.constraint
            priority.flavorPolicy.preemptible
            
        THREAD_PRECEDENCE_POLICY
            priority.flavorPolicy.importance
        
        
        NOTES
        
        Arguments to thread_policy_get() are not well explained by Apple. The are documented by Apple in two places:
        http://developer.apple.com/documentation/Darwin/Conceptual/KernelProgramming/scheduler/chapter_8_section_4.html
        /usr/include/mach/thread_policy.h (though there is another one down /System/Library/Frameworks/)
        
        Info from those sources is unified below. 
        
        1 INPUT	(thread_act_t thread) : The thread ID.  Get this by calling mach_thread_self().
        
        2 INPUT (thread_policy_flavor_t flavor): Allowable values are THREAD_STANDARD_POLICY, THREAD_EXTENDED_POLICY, THREAD_TIME_CONSTRAINT_POLICY and THREAD_PRECEDENCE_POLICY.  
                                                        
                                                Each flavor constant has an associated struct type for holding policy parameters, for example the THREAD_STANDARD_POLICY constant
                                                matches with thread_standard_policy struct.   
                                                
                                                Apple documents THREAD_PRECEDENCE_POLICY: "Another policy call is THREAD_PRECEDENCE_POLICY. This is used for setting the relative 
                                                importance of non-real-time threads. Its calling convention is similar, except that its structure is thread_precedence_policy, and 
                                                contains only one field, an integer_t called importance, which is a signed 32-bit value."  So apparently these flavors are mutually
                                                exclusive. 
                                                
                                                How the flavor argument interacts with the get_default argument is confusing and poorly documented by Apple.  flavor is an 
                                                input-onlyargument to thread_policy_get(), yet it seesm that calling thread_policy_get() with the get_default argument should tell 
                                                us what is the default policy flavor, but it can do that because the flavor argument is input only.
                                                
                                                Here is the current theory about how to discover the current thread setting:
                                                
                                                i. To find the current thread thread policy invoke thread_policy_get() iteratively with each of the three policy flavors and the
                                                get_default argument unasserted. On return thread_policy_get should clear get_default for the current flavor and set get_default for
                                                others.
                                                
                                                ii. The contents of policy_info should be valid for the current mode.  Check the count argument to verify that the ammount of data returns
                                                matches the size of the struct which corresponds to the current flavor mode constant.  
                                                 
                                                note: when invoking thread_policy_get() the info field should point to enough memory to fill the largest of the flavor's corressponding
                                                policy_info struct types: thread_standard_policy, thread_time_constraint_policy and thread_precedence_policy.  
                                                
                                                                                                                                                
        3 INPUT (thread_policy_t policy_info):  Type thread_policy_t defined (ultimately) as a pointer to type long, but the variable policy_info of that type should point to 
                                                a thread_time_constraint_policy struct allocated by the caller.  thread_policy_get() fills the struct memory at the location referenced
                                                by the input pointer.  thread_policy_t looks like this:
                                                
                                                struct thread_time_constraint_policy {
                                                    uint32_t		period;
                                                    uint32_t		computation;
                                                    uint32_t		constraint;
                                                    boolean_t		preemptible;
                                                };
                                                
                                                from Apple's thead_policy.h header file:
                                                
                                                
                                                period: This is the nominal amount of time between separate
                                                processing arrivals, specified in absolute time units.  A
                                                value of 0 indicates that there is no inherent periodicity in
                                                the computation.
                                                
                                                computation: This is the nominal amount of computation
                                                time needed during a separate processing arrival, specified
                                                in absolute time units.
                                            
                                                constraint: This is the maximum amount of real time that
                                                may elapse from the start of a separate processing arrival
                                                to the end of computation for logically correct functioning,
                                                specified in absolute time units.  Must be (>= computation).
                                                Note that latency = (constraint - computation).
                                                
                                                preemptible: This indicates that the computation may be
                                                interrupted, subject to the constraint specified above.

                                                
        4 INPUT/OUTPUT (mach_msg_type_number_t *count):  Type mach_msg_type_number_t is unsigned int.  On input, "the maximum amount of storage that the calling task has 								allocated for the return, [on return]... overwritten by the scheduler to indicate the amount of data that was actually returned."                                        
	
        5 INPUT/OUTPUT 	(boolean_t *get_default):	from Apple's thead_policy.h header file:
                                                        The extra 'get_default' parameter to the second call is
                                                        IN/OUT as follows:
                                                        1) if asserted on the way in it indicates that the default
                                                        values should be returned, not the ones currently set, in
                                                        this case 'get_default' will always be asserted on return;
                                                        2) if unasserted on the way in, the current settings are
                                                        desired and if still unasserted on return, then the info
                                                        returned reflects the current settings, otherwise if
                                                        'get_default' returns asserted, it means that there are no
                                                        current settings due to other parameters taking precedence,
                                                        and the default ones are being returned instead.

*/

#include "MachPriorityMex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    bool				isError;
    thread_policy_flavor_t		flavorConstant;
    int					flavorPolicySize, flavorPolicySizeBytes, kernError;
    task_t				threadID;
    thread_policy_t			threadPolicy;
    mach_msg_type_number_t		policySizeFilled;
    boolean_t				isDefault, getDefault;
    char				commandString[COMMAND_STRING_LENGTH];
    // for return structure
    //  ...outer
    const char *outerNames[] = 		{"threadID", "flavor", "policy", "policySize", "policyFillSize", "getDefault", "isDefault"};
    mwSize numOuterDims=2, numOuterFields=7;
    mwSize outerDims[2]={1,1};
    //	...inner
    const char *standardNames[] = 	{"no_data"};
    int numInnerFieldsStandard=		1;
    /*
    const char *extendedNames[] = 	{"timeshare"};
    int numInnerFieldsExtended=		1;
    */
    const char *timeConstraintNames[]=	{"period", "computation", "constraint", "preemptible"};
    int numInnerFieldsTimeConstraint=	4;
    const char *precedenceNames[]=	{"imporantance"};
    int numInnerFieldsPrecedence=	1;
    mwSize numInnerDims=2;
    mwSize innerDims[2]={1,1};
    // ...both
    mxArray  *tempFieldValue, *innerStruct, *outerStruct;
    

    

    //get the policy flavor constant specified by the user and the getDefault argument
    if(nrhs<2 || nrhs > 2)
        mexErrMsgTxt("MachGetPriorityMex requires two input arguments.  See help MachGetPriorityMex.");
    if(!mxIsChar(prhs[0]))
        mexErrMsgTxt("First input argument is not a string.  See help MachGetPriorityMex.");
    mxGetString(prhs[0], commandString, COMMAND_STRING_LENGTH);
    isError=GetFlavorConstantFromFlavorString(commandString, mxGetM(prhs[0]) * mxGetN(prhs[0]), &flavorConstant);  //case sensitive.  
    if(isError)
        mexErrMsgTxt("Unrecognized command.  See help MachGetPriorityMex.");
    if(!(mxIsDouble(prhs[1]) || mxIsLogical(prhs[1])) || mxGetN(prhs[1]) * mxGetM(prhs[1]) != 1)
        mexErrMsgTxt("Second argument must be 1x1 logical or double value.  See help MachGetPriorityMex.");
    if(mxIsLogical(prhs[1]))
        getDefault= (boolean_t)mxGetLogicals(prhs[1])[0];
    if(mxIsDouble(prhs[1]))
        getDefault= (boolean_t)mxGetPr(prhs[1])[0];


    //read the priority settings
    switch(flavorConstant){
        case THREAD_STANDARD_POLICY: 
            flavorPolicySizeBytes=sizeof(thread_standard_policy_data_t); 
            flavorPolicySize=THREAD_STANDARD_POLICY_COUNT; 
            break;
        case THREAD_TIME_CONSTRAINT_POLICY: 
            flavorPolicySizeBytes=sizeof(thread_time_constraint_policy_data_t); 
            flavorPolicySize=THREAD_TIME_CONSTRAINT_POLICY_COUNT; 
            break;
        case THREAD_PRECEDENCE_POLICY: 
            flavorPolicySizeBytes=sizeof(thread_precedence_policy_data_t); 
            flavorPolicySize=THREAD_PRECEDENCE_POLICY_COUNT; 
            break;
    }
    threadPolicy=(thread_policy_t)malloc(flavorPolicySizeBytes);		
    threadID= mach_thread_self();
    policySizeFilled=flavorPolicySize;
    isDefault=getDefault;
    kernError=thread_policy_get(threadID, flavorConstant, threadPolicy, &policySizeFilled, &isDefault);

    //create and populate the return structure
    outerStruct= mxCreateStructArray(numOuterDims, outerDims, numOuterFields, outerNames);
    
    tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
    mxGetPr(tempFieldValue)[0]=(double)threadID;
    mxSetField(outerStruct, 0, "threadID", tempFieldValue);
    
    tempFieldValue= mxCreateString(commandString);
    mxSetField(outerStruct, 0, "flavor", tempFieldValue);
    
    
    switch(flavorConstant){
        case THREAD_STANDARD_POLICY: 				
            innerStruct= mxCreateStructArray(numInnerDims, innerDims, numInnerFieldsStandard, standardNames);
            tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(tempFieldValue)[0]= (double)((thread_standard_policy_t)threadPolicy)->no_data;
            mxSetField(innerStruct, 0, "no_data", tempFieldValue);
            break;
       /* THREAD_EXTENDED_POLICY is equal to THREAD_STANDARD_POLICY.  Also,  THREAD_EXTENDED_POLICY is undocumented.  So we ignore it.  
        case THREAD_EXTENDED_POLICY: 		
            innerStruct= mxCreateStructArray(numInnerDims, innerDims, numInnerFieldsExtended, extendedNames);
            tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(tempFieldValue)[0]= (double)((thread_extended_policy_t)threadPolicy)->timeshare;
            mxSetField(innerStruct, 1, "timeshare", tempFieldValue);
            break;
        */
        case THREAD_TIME_CONSTRAINT_POLICY: 	
            innerStruct= mxCreateStructArray(numInnerDims, innerDims, numInnerFieldsTimeConstraint, timeConstraintNames);
            tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(tempFieldValue)[0]= (double)((thread_time_constraint_policy_t)threadPolicy)->period;
            mxSetField(innerStruct, 0, "period", tempFieldValue);
            tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(tempFieldValue)[0]= (double)((thread_time_constraint_policy_t)threadPolicy)->computation;
            mxSetField(innerStruct, 0, "computation", tempFieldValue);
            tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(tempFieldValue)[0]= (double)((thread_time_constraint_policy_t)threadPolicy)->constraint;
            mxSetField(innerStruct, 0, "constraint", tempFieldValue);
            tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(tempFieldValue)[0]= (double)((thread_time_constraint_policy_t)threadPolicy)->preemptible;
            mxSetField(innerStruct, 0, "preemptible", tempFieldValue);
            break;
        case THREAD_PRECEDENCE_POLICY:
            innerStruct= mxCreateStructArray(numInnerDims, innerDims, numInnerFieldsPrecedence, precedenceNames);
            tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(tempFieldValue)[0]= (double)((thread_precedence_policy_t)threadPolicy)->importance;
            mxSetField(innerStruct, 0, "imporantance", tempFieldValue);
            break;
    }
    mxSetField(outerStruct,0, "policy", innerStruct);
    
    tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
    mxGetPr(tempFieldValue)[0]=flavorPolicySize;
    mxSetField(outerStruct, 0, "policySize", tempFieldValue);
    
    tempFieldValue= mxCreateDoubleMatrix(1,1,mxREAL);
    mxGetPr(tempFieldValue)[0]=policySizeFilled;
    mxSetField(outerStruct, 0, "policyFillSize", tempFieldValue);
    
    tempFieldValue= mxCreateLogicalMatrix(1, 1);
    mxGetLogicals(tempFieldValue)[0]=(bool)getDefault;
    mxSetField(outerStruct, 0, "getDefault", tempFieldValue);
    
    tempFieldValue= mxCreateLogicalMatrix(1, 1);
    mxGetLogicals(tempFieldValue)[0]=(bool)isDefault;
    mxSetField(outerStruct, 0, "isDefault", tempFieldValue);
    
    plhs[0]=outerStruct;
    
    free((void*)threadPolicy);        
}
