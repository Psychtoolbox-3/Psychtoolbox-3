/*
  PsychToolbox3/Source/Common/PsychHelp.c

  AUTHORS:

    Allen.Ingling@nyu.edu           awi
    mario.kleiner.de@gmail.com      mk
*/

#include "Psych.h"

//declare static variables for internal use by PsychHelp.c
static char *functionUseHELP = NULL;
static char *functionSynopsisHELP = NULL;
static char *functionSeeAlsoHELP = NULL;
static psych_bool giveHelpHELP = FALSE;
static psych_bool oneShotReturnHelp = FALSE;

// functions for flipping a flag to indicate whether function help should be given.
void PsychSetGiveHelp(void)
{
    giveHelpHELP = TRUE;
}

void PsychClearGiveHelp(void)
{
    giveHelpHELP = FALSE;
}

psych_bool PsychIsGiveHelp(void)
{
    return(giveHelpHELP);
}

// push the help strings onto a shallow stack 1 element deep
void PsychPushHelp(char *functionUse, char *functionSynopsis, char *functionSeeAlso)
{
    functionUseHELP = functionUse;
    functionSynopsisHELP = functionSynopsis;
    functionSeeAlsoHELP = functionSeeAlso;
}

void PsychOneShotReturnHelp(void)
{
    oneShotReturnHelp = TRUE;
}

void PsychGiveHelp(void)
{
    PsychGenericScriptType  *cellVector;

    // Special case: Asked to return help in a cell array of strings?
    if (oneShotReturnHelp) {
        // Yes. Return a 3 element cell array of strings, each containing one
        // of the three help text arguments:
        PsychAllocOutCellVector(1, FALSE, 3,  &cellVector);
        PsychSetCellVectorStringElement(0, functionUseHELP, cellVector);
        PsychSetCellVectorStringElement(1, BreakLines(functionSynopsisHELP, 80), cellVector);
        PsychSetCellVectorStringElement(2, functionSeeAlsoHELP, cellVector);
        oneShotReturnHelp = FALSE;

        return;
    }

    // No, standard path: Print to console of runtime system:
    printf("\nUsage:\n\n%s\n",functionUseHELP);

    if (functionSynopsisHELP != NULL) {
        // Underwhelmingly, Python does not allow to printf() more than 1000 characters per
        // string, and our functionSynopsisHELP can be much longer. Manually break it up in
        // substrings and printf() one line per printf() invocation, assuming no single
        // line will exceed 1000 chars:
        char *tok = NULL;
        char *str = strdup(functionSynopsisHELP);
        str = BreakLines(str, 80);
        tok = strtok(str, "\n");

        if (tok)
            printf("\n");

        while (tok) {
            printf("%s\n", tok);
            tok = strtok(NULL, "\n");
        }

        free(str);
    }

    if (functionSeeAlsoHELP  != NULL) printf("\nSee also: %s\n", BreakLines(functionSeeAlsoHELP, 80));
}

void PsychGiveUsage(void)
{
    printf("Usage:\n\n%s",functionUseHELP);
}

void PsychGiveUsageExit(void)
{
    PrintfExit("Usage:\n\n%s",functionUseHELP);
}

#if PSYCH_LANGUAGE == PSYCH_PYTHON
// TODO FIXME PYTHON Enable this also for other scripting environments, once all modules are
// converted to the new const char **synopsisSYNOPSIS = InitializeSynopsis(void); syntax.
const char* PsychBuildSynopsisString(const char* modulename)
{
    // Build module help string and functions synopsis in synopsisSYNOPSIS:
    const char **synopsisSYNOPSIS = InitializeSynopsis();
    int i, n = 0;
    char* moduleHelpString = NULL;

    // Assemble the moduleHelpString out of synopsisSYNOPSIS, suitable for
    // Python module help:
    for (i = 0; synopsisSYNOPSIS[i] != NULL; i++)
        n+= (int) strlen(synopsisSYNOPSIS[i]) + 2;

    moduleHelpString = calloc(1, n + 2 * strlen(modulename) + 512);
    if (PSYCH_LANGUAGE == PSYCH_PYTHON) {
        strcat(moduleHelpString, "Copyright (c) 2018 Mario Kleiner. Licensed under MIT license.\n\n");
    }

    strcat(moduleHelpString, "For detailed help on a subfunction SUBFUNCTIONNAME, type " );
    strcat(moduleHelpString, modulename);
    strcat(moduleHelpString, "('SUBFUNCTIONNAME?')\n");
    strcat(moduleHelpString, "ie. the name with a question mark appended. E.g., for detailed help on the subfunction\n" );
    strcat(moduleHelpString, "called Version, type this: ");
    strcat(moduleHelpString, modulename);
    strcat(moduleHelpString, "('Version?')\n\n");
    for (i = 0; synopsisSYNOPSIS[i] != NULL; i++) {
        strcat(moduleHelpString, synopsisSYNOPSIS[i]);
        strcat(moduleHelpString, "\n");
    }

    return(moduleHelpString);
}
#endif
