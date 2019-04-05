#!/usr/bin/env python3
#
# Copyright (C) 2008 Tobias Wolf <towolf@tuebingen.mpg.de>
#
# Copyright (C) 2018 Justin Ales <justin.ales@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
"""
PTB-wikify-into-files.py: Extract help strings in .m-files and write help into files
Supports markdown and mediawiki outputs.

Usage: PTB-wikify-into-files.py [options] ([m-files]|[directory])

Options:
 -h, --help
 -o [dir], --output-dir     Output directory to write files into.
 -f [format], --format      Output file format, can be either markdown (default), or mediawiki
 -r, --recursive            Recursive mode: use only with a single directory
 -m, --mexmode              Mex mode: also look for .mexmaci files and post their
                            help strings by calling MATLAB and running them
                            In recursive mode both M and Mex files are posted.
                            (edit the source to change _mexext into one of
                            .mexglx, .mexmaci, or .dll)
 -g, --opengl               Do generate documentation for PsychOpenGL folder as well.
                            By default, docs for glXXX() functions in PsychOpenGL/MOGL/
                            are not generated, as there are too many of them for the
                            GitHub Wiki to handle them.
 -v                         be verbose (massive text output)
                            this prints out:
                              - which files were skipped
                              - a diff of the text before submission


This script is designed to work with a GitHub wiki. GitHub exposes its wiki as
git repository with each page corresponding to a file.  This script will write
out all files into a directory that can be pushed into the wiki repository to
update the wiki.

Examples:

This is currently kind of tricky to use because care must be taken with the 2
separate git repositories, the wiki and the code. The script should be run from
the Psychtoolbox root directory and give an path to an output directory outside
the psychtoolbox directory.

Example to, Clone the wiki, than extract files from PsychBasic,
than push changes back onto wiki:

  cd ~/git
  git clone https://github.com/Psychtoolbox-3/psychtoolbox-3.github.com.git
  cd /path/to/Psychtoolbox
  PTB-wikify-into-files.py -m -o ~/git/psychtoolbox-3.github.com/docs/ PsychBasic/*.m

  cd ~/git/Psychtoolbox-3.wiki
  git add -A
  git commit -m "Update Message"
  git push

Typical invocation to update the docs:

cd /path/to/Psychtoolbox
../managementtools/PTB-wikify-into-files.py -r -m -o ~/git/psychtoolbox-3.github.com/docs/ ./;

For kleinerm:

../managementtools/PTB-wikify-into-files.py -r -m -o ~/projects/OpenGLPsychtoolbox/psychtoolbox-3.github.com/docs/ ./

Alternatively recursively add all files in the directory:

  PTB-wikify-into-files.py -r -m -o ~/git/psychtoolbox-3.github.com/docs/ PsychBasic/

Example script that outputs entire PTB tree:

cd /path/to/Psychtoolbox
../managementtools/PTB-wikify-into-files.py -m -o ~/git/psychtoolbox-3.github.com/docs/ *.m;

IMPORTANT!!:
  Always change your working directory to the root of the
  tree before running the script, e.g.,
    cd Psychtoolbox
    ~/PTB-wikify-into-files.py -r PsychDemos/PsychExampleExperiments

   (This is to figure out what the root node is supposed to be.)

"""

import sys, os, re, subprocess
import getopt
import textwrap

outputdir = "./"
_recursive = 0
_mexmode = 0
_ignoreOpenGL = 1
outputFormat = "markdown"

from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    _mexext = '.mexa64'
elif platform == "darwin":
    # OS X
    _mexext = '.mexmaci64'
elif platform == "win32":
    # Windows...
    _mexext = '.mexmw64'

_debug = 0
_fulldiff = 0

# This version writes files instead of directly posting onto website.

#Format a string for making a link in different markdown dialects
def formatLinkText(*args):

    if len(args)==1:
        label=args[0]
        link=args[0]
    elif len(args)==2:
        label=args[0]
        link=args[1]


    if outputFormat == "markdown":
        formattedText = "["+label + "](" + link + ")"
    elif outputFormat == "mediawiki":
        formattedText = "[["+link+"|"+label+"]]"

    return formattedText

#Function to split a path into each different directory.
def splitall(path):
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


def parse(filename):
    '''Parse helpful initial comment block from .m-file'''
    docstring = []
    functionellipsized = False
    for line in open(filename,'r',newline=None):
        # weird hack
        if line.strip()[0:8] + line.strip()[-3:] == 'function...':
            functionellipsized = True
        if line.strip().startswith('%'):       # get comment lines
            docstring.append(line)
        elif len(docstring): break             # and stop after one block
            # we allow empty lines and the function def before the block
        elif len(line.strip()) > 0 \
                and not line.strip().startswith('function')\
                and not functionellipsized: break

    # strip first two chars: '% '
    return textwrap.dedent(''.join([l[1:] for l in docstring]))

def beackern(mkstring, doLinks):
    '''Regex cleaning of PTB docstrings'''
    # expand tabs to four spaces
    mkstring = mkstring.expandtabs(4)
    # Enclose text markup characters with: ""
    #mkstring = re.sub(r'(\+\+|(--){2}\b|==|\'\'|@@|[^:]//|>>|<<|##)',r'""\1""',mkstring)
    # escape markdown characters with \
    mkstring = re.sub(r'(#|\*|_|>)',r'\\\1',mkstring)

#Double new lines for to break lines for markdown.
#    mkstring = re.sub(r'(\n)',r'\n\n',mkstring)

    # Replace 'Windows: ______' with headline
    mkstring = re.sub(r'(?m)^(.*):\s*[_-]{3,}',r'# \1',mkstring)
    # Replace long _______ with ---- i.e.<HR>
    mkstring = re.sub(r'(?m)^\s*_{10,}',r'----',mkstring)
    # Replace some non-WikiWord function to internal links
    match = r'(\bScreen\b|' \
            + r'\bAsk\b|' \
            + r'\bBandpass\b|' \
            + r'\bBeeper\b|' \
            + r'\bCircle\b|' \
            + r'\bClose\b|' \
            + r'\bComputeDKL_M\b|' \
            + r'\bDKLDemo\b|' \
            + r'\bEllipse\b|' \
            + r'\bExpand\b|' \
            + r'\bExplode\b|' \
            + r'\bFlip\b|' \
            + r'\bGestalt\b|' \
            + r'\bGStreamer\b|' \
            + r'\bInterleave\b|' \
            + r'\bM_PToP\b|' \
            + r'\bM_PToT\b|' \
            + r'\bM_TToP\b|' \
            + r'\bM_TToT\b|' \
            + r'\bOSAUCSTest\b|' \
            + r'\bOSName\b|' \
            + r'\bPriority\b|' \
            + r'\bPsychometric\b|' \
            + r'\bQuest\b|' \
            + r'\bRandi\b|' \
            + r'\bRanint\b|' \
            + r'\bReplace\b|' \
            + r'\bResolute\b|' \
            + r'\bRush\b|' \
            + r'\bSample\b|' \
            + r'\bScale\b|' \
            + r'\bShuffle\b|' \
            + r'\bSnd\b|' \
            + r'\bStopwatch\b|' \
            + r'\bTrunc\b|' \
            + r'\bXYZTouv\b|' \
            + r'\bhexstr\b|' \
            + r'\bmoglmorpher\b|' \
            + r'\bpsychassert\b|' \
            + r'\bpsychlasterror\b|' \
            + r'\bpsychrethrow\b|' \
            + r'\bsca\b|' \
            + r'\buvTols\b|' \
            + r'\buvToxy\b|' \
            + r'\buvYToXYZ\b|' \
            + r'\bxyTouv\b|' \
            + r'\bxyYToXYZ\b|' \
            + r'\blog10nw\b|' \
            + r'\bPreference\b)'

    if doLinks:
        if outputFormat == "markdown":
            mkstring = re.sub(match,r'[\1](\1)',mkstring)
        elif outputFormat == "mediawiki":
            mkstring = re.sub(match,r'[[\1|\1]]',mkstring)

        #Add links for any word using UpperCamelCase: e.g. PsychHID
        #BUT does NOT match any string stating wiht ' e.g. 'OpenWindow'
        #Because we don't want subfunctions/strings to trigger page links
        #WARNING:
        #This is a pretty crazy regex to .  Cobled together from lots of
        #web searches and tests.  Probably not very robust but gets most of the job done.
        UpperCamelMatch = r'(?<!\')(?:\(|\b)[A-Z]([A-Z0-9]*[a-z][a-z0-9]*[A-Z]|[a-z0-9]*[A-Z][A-Z0-9]*[a-z])[A-Za-z0-9]*(?:\)|\b)(?!\')'
        #r'(?<!\')[A-Z]([A-Z0-9]*[a-z][a-z0-9]*[A-Z]|[a-z0-9]*[A-Z][A-Z0-9]*[a-z])[A-Za-z0-9]*'

        if outputFormat == "markdown":
            mkstring = re.sub(UpperCamelMatch,r'[\g<0>](\g<0>)',mkstring)
        elif outputFormat == "mediawiki":
            mkstring = re.sub(UpperCamelMatch,r'[[\g<0>|\g<0>]]',mkstring)

    # purge useless help lines
    mkstring = re.sub(r'(?m)^.*triple-click me & hit enter.*$\n','',mkstring)
    # Generate more headlines (Word:$, CAPITAL LETTERS$, etc)
    mkstring = re.sub(r'(?m)((?<=\n\n).+:((?=\s*\n\n)+))', r'### \1',mkstring)
    mkstring = re.sub('(?m)((?<=\n\n)\ *[A-Z][A-Z :!-]+\ *(?=\n))', r'### \1',mkstring)

    #media wiki interprets whitespace at beggining of line as code block
    if outputFormat == "mediawiki":
        mkstring = re.sub(r'^\s+','',mkstring,flags=re.M)

    if outputFormat == "markdown":
    # Add double space after lines for markdown to keep paragraph togethr
        mkstring = re.sub(r'(\n)',r'  \n',mkstring)
    elif outputFormat == "mediawiki":
    #Add manual linebreaks to preserve formatting
        mkstring = re.sub(r'(\n)',r'<br />\n',mkstring)

    return mkstring

def writeFiles(outputDir,files):
    ''' create pages for all files in there and link to the category'''
    for name in files:

        from pathlib import Path
        #Get absolute path for checking file parent.
        absPath = os.path.abspath(name)
        category = Path(Path(absPath).parent).stem

        # single out some names
        head, basename = os.path.split(name)
        root =  os.path.join(os.path.basename(os.getcwd()))
        path =   os.path.normpath(os.path.join(root,head))
        funcname, ext = os.path.splitext(basename)

        cattext = ''

        mexname = os.path.join(os.path.split(absPath)[0],funcname+_mexext)

        mexSuccess = True
        if not os.path.exists(name) or not ext=='.m':
            if _mexmode and ext == _mexext:
                mexSuccess = mexhelpextract(outputDir,[funcname])

            if mexSuccess:
                continue
            else:
                name = os.path.join(head,funcname + '.m')
                if os.path.exists(name):
                    print('No help in mexfile. Trying m file: ' + name )
                else:
                    print('Skipping because cannot find help in mex or mfile for: ' + os.path.join(head, basename))
                    continue

        #If mexmode is on and mex file exists skip the .m file
        if mexSuccess and (_mexmode and os.path.exists(mexname)):
            print('Mex file exists. skipping m file: ' + name )
            continue


        if basename=='Contents.m' or basename=='contents.m':
            funcname = category
            # path = os.path.dirname(path)
            # #Sooo Kludgy.  SHould update directory parsing.
            # if len(path) < 3:
            #     funcname = 'Psychtoolbox'
            cattext = '{{category}}'

        # load and build the text for the page

        headline = "# " + formatLinkText(funcname)+ "\n"


        breadcrumb = "##### "
        dirnames = splitall(path)

        for dirname in dirnames:
            breadcrumb = breadcrumb + ">" + formatLinkText(dirname)

        breadcrumb = breadcrumb + "\n\n"


        # read the .m-file and strip all the crap
        docstring = parse(name)

        # scrub the text real good, to get us some nice wiki formatting
        if docstring:
            body = beackern(docstring, 1)
        else:
            body = 'This function is not yet documented.\n\n\n MissingDocs'

        pathlinks = """
                    <div class="code_header" style="text-align:right;">
                      <span style="float:left;">Path&nbsp;&nbsp;</span> <span class="counter">Retrieve <a href=
                      "https://raw.github.com/Psychtoolbox-3/Psychtoolbox-3/beta/%s">current version from GitHub</a> | View <a href=
                      "https://github.com/Psychtoolbox-3/Psychtoolbox-3/commits/beta/%s">changelog</a></span>
                    </div>
                    <div class="code">
                      <code>%s</code>
                    </div>
                    """% tuple([os.path.join(path,basename)]*3)

        text =  headline \
                + breadcrumb \
                + body \
                + '\n\n\n' \
                + textwrap.dedent(pathlinks) \
                + '\n' \
                + cattext

        if outputFormat == "markdown":
            suffix = '.md'
        elif outputFormat == "mediawiki":
            suffix = '.wiki'

        wikiFileName = os.path.join(outputDir,funcname + suffix)
        print("writing file %s"%wikiFileName)

        f = open(wikiFileName,"w+")
        f.write(text)



def parseMexExtract(filename):
    from collections import defaultdict
    sectionPattern = re.compile(r'\[section:(.*?)\]')
    keyPattern = re.compile(r'\[key:(.*?)\]')

    extractedText = defaultdict(lambda: defaultdict(str))

    with open(filename) as toParse:

        for line in toParse:

            sectionMatch = sectionPattern.match(line)
            keyMatch     = keyPattern.match(line)
            if sectionMatch:
                sectionName = sectionMatch.groups()[0]
            elif keyMatch:
                keyName = keyMatch.groups()[0]
            else:
                extractedText[sectionName][keyName]=extractedText[sectionName][keyName]+line

    return extractedText


def formatUsage(usage,mexname,subfunctions):

    if outputFormat == "markdown":
    # Add double space after lines for markdown to keep paragraph togethr
        usage = re.sub(r'(\n)',r'  \n',usage)
    elif outputFormat == "mediawiki":
    #Add manual linebreaks to preserve formatting
        usage = re.sub(r'(\n)',r'<br />\n',usage)

    #Now parse the usage block to add links to sub function help
    #Looking for pattern: mexfile('subfunc'

    for subFuncName,subFuncText in subfunctions.items():


        if outputFormat == "markdown":
            usage = re.sub( \
            mexname+'\(\'' +subFuncName + '\'', \
            mexname+'(\'[' + subFuncName + '](' + mexname + '-' + subFuncName + ')\'', \
            usage)
        elif outputFormat == "mediawiki":
            usage = re.sub( \
            mexname+'\(\'' +subFuncName + '\'', \
            mexname+'(\'[[' + mexname + '-' + subFuncName + '|'  + subFuncName + ']]\'', \
            usage)



    return usage

def mexhelpextract(outputDir,mexnames):
    #print 'processing mex files: ' + mexnames.__repr__()
    mexHelpSuccess = True

    for mexname in mexnames:
        # First build help text for mex.

        # assemble command line for matlab
        matlabcmd = 'addpath(\'%s\');%s(\'%s\',\'%s\'); exit' % \
            (_tmpdir, \
             os.path.splitext(os.path.basename(_mexscript))[0], \
             mexname, \
             _tmpdir)
        # and execute Octave or Matlab w/ the temporary script we wrote earlier
        try:
            # Try Octave first:
            cmd = 'octave --no-history --no-window-system --silent --no-gui --eval "%s" > /dev/null' % matlabcmd
            print('running OCTAVE for %s in %s' % (mexname,_tmpdir))
            p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, close_fds=True)
            stderr = p.communicate()[1]
            if stderr: print(stderr)
        except:
            # Failed: Try Matlab:
            try:
                cmd = 'matlab -nojvm -nodisplay -r "%s" > /dev/null' % matlabcmd
                print('running MATLAB for %s in %s' % (mexname,_tmpdir))
                p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, close_fds=True)
                stderr = p.communicate()[1]
                if stderr: print(stderr)
            except:
                print('could not dump help for %s into %s.' % (mexname,_tmpdir))

        mexDumpName = os.path.join(_tmpdir,mexname);

        if not os.path.isfile(mexDumpName):
            print("skipping " + mexname + " (no output)")
            mexHelpSuccess = False
            continue

        # Parse the mex dumped file using the following format:
        # [Section:SubFunction]
        # [key:usage]
        # usage text
        # [key:help]
        # help text
        # [key:seealso]
        # see  also text

        subfunctions = parseMexExtract(mexDumpName)

        subFuncList = ''
        for key in subfunctions:
            subFuncList = subFuncList + ', ' + key

        print('processing subfunctions: ' + subFuncList)

        for subFuncName,subFuncText in subfunctions.items():
            # read in the strings for this subfunction
            usage = subFuncText['usage']
            help = subFuncText['help']
            seealso = subFuncText['seealso']

            headline = '# ['+mexname+'(\''+subFuncName+'\')](' \
                        + mexname + '-' + subFuncName + ') ' + '\n'

            breadcrumb = "##### " + formatLinkText('Psychtoolbox','Psychtoolbox') + ">"  \
                                + formatLinkText(mexname) + ".{mex*} subfunction\n\n"

            if subFuncName == '__main__':
                usage = formatUsage(usage,mexname,subfunctions)
                headline = "# " + formatLinkText(mexname,mexname) + "\n"
                breadcrumb = "##### " + formatLinkText('Psychtoolbox','Psychtoolbox') + ">"  \
                                    + formatLinkText(mexname) + "\n\n"
                # Scrub the text for main help text of mex file only - without keyword substitution for main function help:
                body = beackern(help, 0)
            else:
                # Scrub the text for main help text of subfunction only - with keyword substitution:
                body = beackern(help, 1)

            # docstring = '' \
            #         + '%%(matlab;Usage)' \
            #         + usage \
            #         + '%%\n' \
            #         + body \
            #         + '\n\n'
            docstring = ''  \
                    + usage \
                    + '\n' \
                    + body \
                    + '\n\n'

            if seealso:

                #Add links for any word using UpperCamelCase: e.g. PsychHID
                #BUT does NOT match any string stating wiht ' e.g. 'OpenWindow'
                #Because we don't want subfunctions/strings to trigger page links
                #WARNING:
                #This is a pretty crazy regex to .  Cobled together from lots of
                #web searches and tests.  Probably not very robust but gets most of the job done.
                UpperCamelMatch = r'(?<!\')(?:\(|\b)[A-Z]([A-Z0-9]*[a-z][a-z0-9]*[A-Z]|[a-z0-9]*[A-Z][A-Z0-9]*[a-z])[A-Za-z0-9]*(?:\)|\b)(?!\')'
                #r'(?<!\')[A-Z]([A-Z0-9]*[a-z][a-z0-9]*[A-Z]|[a-z0-9]*[A-Z][A-Z0-9]*[a-z])[A-Za-z0-9]*'

                if outputFormat == "markdown":
                    seealso = re.sub(UpperCamelMatch,r'[\g<0>]('+mexname + '-\g<0>)',seealso)
                elif outputFormat == "mediawiki":
                    seealso = re.sub(UpperCamelMatch,r'[['+mexname + '-\g<0>|\g<0>]]',seealso)


                docstring = docstring + '###See also:\n' + seealso

            text = headline  \
                + breadcrumb \
                + docstring


            if outputFormat == "markdown":
                suffix = '.md'
            elif outputFormat == "mediawiki":
                suffix = '.wiki'

            outputFilename = os.path.join(outputDir,mexname + '-' + subFuncName + suffix)
            if subFuncName == '__main__':
                outputFilename = os.path.join(outputDir,mexname + suffix)

            print("writing mex help to file: " + outputFilename)

            f = open(outputFilename,"w+")
            f.write(text)

    return mexHelpSuccess

def recursivewalk(outputDir,rootfolder):
    '''traverse the root directory and post all .m-files'''
    for root, dirs, files in os.walk(rootfolder):

        if re.search(os.sep+'\.svn|'+os.sep+'private',root):
                continue

        if _ignoreOpenGL and re.search(os.sep+'MOGL',root):
                continue

        print('Entering folder ' + root)

        # .m-files are processed with postsinglefiles
        mfiles = [os.path.join(root,f) for f in files \
                if f.endswith(('.m','.M',_mexext))]
        writeFiles(outputDir,mfiles)

    print("Exiting: Done with this tree.")
    sys.exit(0)

def usage():
    print(__doc__)

def main(argv):

    try:
        opts, args = getopt.getopt(argv, "o:f:hrmvg", \
            ["help", "output-dir=", "format=", "opengl", "recursive", "mexmode"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-o", "--output-dir"):
            global outputDir
            outputDir = arg
        elif opt in ("-f", "--format"):
            global outputFormat
            print(arg)

            if arg.lower() not in ('markdown','mediawiki'):
                print('Error: output format must be either markdown or mediawiki\n\n')
                sys.exit(1)
            outputFormat = arg.lower()
        elif opt in ("-r", "--recursive"):
            global _recursive
            _recursive = 1
        elif opt in ("-g", "--opengl"):
            global _ignoreOpenGL
            _ignoreOpenGL = 0
        elif opt in ("-m", "--mexmode"):
            import tempfile
            global _mexmode, _mexext, _tmpdir, _mexscript
            _mexmode = 1
            _tmpdir = tempfile.mkdtemp(prefix="PTB-doc-mexfunctionhelp-")
            _mexscript = os.path.join(_tmpdir, 'PsychtoolboxMexhelpextract.m')
            fid = open(_mexscript, 'w')
            script = '''\
            function PsychtoolboxMexhelpextract(mexname,tmpdir)
            try
                subfunctions=eval([mexname '(''DescribeModulefunctionshelper'')']);
            catch
                subfunctions = [];
            end

            fid = fopen(fullfile(tmpdir,mexname),'wt');
            if fid == -1
                return;
            end

            try
                % First extract help for main function, including subfunction list
                fprintf(fid,'[section:__main__]\\n');
                synopsisText = evalc([mexname ';']);
                helpText     = help([mexname '.m']);
                fprintf(fid,'[key:usage]\\n%s\\n',synopsisText);
                fprintf(fid,'[key:help]\\n%s\\n',helpText);
                fprintf(fid,'[key:seealso]\\n');

                if iscell(subfunctions) && ~isempty(subfunctions)
                    for i=1:size(subfunctions,2)
                        fprintf(fid,'[section:%s]\\n',subfunctions{i});
                        docs = eval([mexname '(''DescribeModulefunctionshelper'',1,subfunctions{i})']);
                        fprintf(fid,'[key:usage]:\\n%s\\n',docs{1});
                        fprintf(fid,'[key:help]\\n%s\\n',docs{2});
                        fprintf(fid,'[key:seealso]\\n%s\\n',docs{3});
                    end
                end
            catch
                fclose(fid);
                return;
            end

            fclose(fid);
            return
            '''
            fid.write(textwrap.dedent(script))
            fid.close()

        elif opt == '-v':
            global _debug
            _debug = 1

    if args:
        if _recursive:
            if len(args) == 1 and os.path.isdir(args[0]):
                recursivewalk(outputDir,args[0])
            else:
                usage()
                print('Error: Recursive mode works with one directory only\n\n')
                sys.exit(1)

        else:
            print("Parsing Input: " + str(args[0]) + " writing to: " + outputDir)
            writeFiles(outputDir,args)
    else:
        usage()
        print('Error: No files specified to submit to Wiki\n\n')
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])
