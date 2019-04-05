/*
  Psychtoolbox3/Source/Common/MODULEVersion.c		
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
  
  PLATFORMS:	
  This file should build on any platform. 

  HISTORY:
  3/12/03  awi		Created. 
 
  DESCRIPTION:
  
  Returns the version of a module in a struct. This includes the build number
  
  TO DO:

  
*/

#include "Psych.h"

static char seeAlsoString[] = "";

PsychError MODULEVersion(void) 
{
	int i;
	PsychAuthorDescriptorType   *author;
    const char *versionFieldNames[]={"version", "major", "minor", "point", "build", "date", "time", "module", "project", "os", "language", "authors"};
    const char *authorFiledNames[]={"first", "middle", "last", "initials", "email", "url"};
    char 	*versionString;
    int		buildNumber;
    int 	numVersionFieldDimensions=-1, numVersionFieldNames=12, numAuthorFieldNames=6, numAuthors;
    PsychGenericScriptType	*versionStructArray, *authorStructArray;
    //we ignore the usual usage help strings and create our own based on the module name. MODULEVersion() is for use by any Psychtoolbox module. 
    char useString[256], synopsisString[256], *moduleName;
    char useStringP1[]="struct=";
    char useStringP2[]="('Version')";
    char synopsisStringP1[]="return the version of ";
    char synopsisStringP2[]=" in a struct";
    
    
    //for generic usage we modifiy at runtiome the help string to replace "Screen" with the name of this module.
    moduleName=PsychGetModuleName();
    useString[0]='\0';
    strcat(useString, useStringP1);
    strcat(useString, moduleName);
    strcat(useString, useStringP2);
    synopsisString[0]='\0';
    strcat(synopsisString, synopsisStringP1);
    strcat(synopsisString, moduleName);
    strcat(synopsisString, synopsisStringP2);
   
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(0));
    
    //get the build and version string
    buildNumber=PsychGetBuildNumber();
    versionString=PsychGetVersionString();

    //create a structure and populate it.
    PsychAllocOutStructArray(1, FALSE, numVersionFieldDimensions, numVersionFieldNames, versionFieldNames, &versionStructArray);
    PsychSetStructArrayStringElement("version", 0, versionString, versionStructArray);
    PsychSetStructArrayDoubleElement("major", 0, (double)PsychGetMajorVersionNumber(), versionStructArray);
    PsychSetStructArrayDoubleElement("minor", 0, (double)PsychGetMinorVersionNumber(), versionStructArray);
    PsychSetStructArrayDoubleElement("point", 0, (double)PsychGetPointVersionNumber(), versionStructArray);
    PsychSetStructArrayDoubleElement("build", 0, buildNumber, versionStructArray);
    PsychSetStructArrayStringElement("date", 0, PsychGetBuildDate(), versionStructArray);
    PsychSetStructArrayStringElement("time", 0, PsychGetBuildTime(), versionStructArray);
    PsychSetStructArrayStringElement("module", 0, moduleName, versionStructArray);
    PsychSetStructArrayStringElement("project", 0, PSYCHTOOLBOX_PROJECT_NAME, versionStructArray);
    PsychSetStructArrayStringElement("os", 0, PSYCHTOOLBOX_OS_NAME, versionStructArray);
    PsychSetStructArrayStringElement("language", 0, PSYCHTOOLBOX_SCRIPTING_LANGUAGE_NAME, versionStructArray);

	numAuthors=PsychGetNumModuleAuthors();
    PsychAllocOutStructArray(-1, FALSE, numAuthors, numAuthorFieldNames, authorFiledNames, &authorStructArray);
	for(i=0;i<numAuthors;i++){
		GetModuleAuthorDescriptorFromIndex(i, &author);
		PsychSetStructArrayStringElement("first", i, author->firstName, authorStructArray);
		PsychSetStructArrayStringElement("middle", i, author->middleName, authorStructArray);
		PsychSetStructArrayStringElement("last", i, author->lastName, authorStructArray);
		PsychSetStructArrayStringElement("initials", i, author->initials, authorStructArray);
		PsychSetStructArrayStringElement("email", i, author->email, authorStructArray);
		PsychSetStructArrayStringElement("url", i, author->url, authorStructArray);
	}
    PsychSetStructArrayStructElement("authors",0, authorStructArray, versionStructArray);

    return(PsychError_none);	
}


