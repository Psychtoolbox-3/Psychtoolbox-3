#!/usr/bin/env python
#
# Copyright (C) 2008 Tobias Wolf <towolf@tuebingen.mpg.de>
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
PTB-wikify.py: Submit help strings in .m-files to a Wikka Installation

Usage: PTB-wikify.py [options] ([m-files]|[directory])

Options:
 -h, --help
 -u ..., --username ...   login to WikkaWiki with this UserName
 -p ..., --password ...   use this Wikka password
 -U ..., --URL            root URL of the Wiki in the form
                            http://doc.psychtoolbox.org/ (w/ terminal slash)
 -r,  --recursive         recursive mode: use only with a single directory
 -m,  --mexmode           mex mode: also look for .mexmaci files and post their
                            help strings by calling MATLAB and running them
                            In recursive mode both M and Mex files are posted.
                            (edit the source to change _mexext into one of
                             .mexglx, .mexmaci, or .dll)
 -v                       be verbose (massive text output)
                            this prints out:
                              - which files were skipped
                              - a diff of the text before submission
 -f   --full-diff         output full ndiff (default is a terse unified diff).
                            The ndiff contains all text but is easier
                            to parse visually for differences.
                            (only in combination with -v)

Examples:
  PTB-wikify.py -u DocBot -p dokkbot -U http://wiki/ PsychBasic/*.m
  PTB-wikify.py -r -m -U http://wiki/ PsychBasic/

IMPORTANT!!:
  Always change your working directory to the root of the
  tree before running the script, e.g.,
    cd Psychtoolbox
    ~/PTB-wikify.py -r PsychDemos/PsychExampleExperiments

   (This is to figure out what the root node is supposed to be.)

"""

import sys, os, re, subprocess
import getopt
from urllib2 import HTTPError

import mechanize
assert mechanize.__version__ >= (0, 0, 6, "a")

import textwrap

from  BeautifulSoup import BeautifulSoup, Tag, NavigableString

baseurl = "http://docs.psychtoolbox.org/wikka.php?wakka="
username = "DocBot"
password = ""
_recursive = 0
_mexmode = 0
_mexext = '.mexa64'
_debug = 0
_fulldiff = 0

global mech
mech = mechanize.Browser()

def parse(filename):
    '''Parse helpful initial comment block from .m-file'''
    docstring = []
    functionellipsized = False
    for line in open(filename,'rU'):
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

def beackern(mkstring):
    '''Regex cleaning of PTB docstrings'''
    # expand tabs to four spaces
    mkstring = mkstring.expandtabs(4)
    # Enclose text markup characters with: ""
    mkstring = re.sub(r'(\+\+|(--){2}\b|==|\'\'|@@|[^:]//|>>|<<|##)',r'""\1""',mkstring)
    # Replace 'Windows: ______' with headline
    mkstring = re.sub(r'(?m)^(.*):\s*[_-]{3,}',r'=== \1 ===',mkstring)
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
            + r'\bIOPort\b|' \
            + r'\bInterleave\b|' \
            + r'\bM_PToP\b|' \
            + r'\bM_PToT\b|' \
            + r'\bM_TToP\b|' \
            + r'\bM_TToT\b|' \
            + r'\bOSAUCSTest\b|' \
            + r'\bOSName\b|' \
            + r'\bPriority\b|' \
            + r'\bPsychometric\b|' \
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
            + r'\log10nw\b|' \
            + r'\bPreference\b)'
    mkstring = re.sub(match,r'[[\1]]',mkstring)
    # purge useless help lines
    mkstring = re.sub(r'(?m)^.*triple-click me & hit enter.*$\n','',mkstring)
    # Generate more headlines (Word:$, CAPITAL LETTERS$, etc)
    mkstring = re.sub(r'(?m)((?<=\n\n).+:((?=\s*\n\n)+))', r'=== \1 ===',mkstring)
    mkstring = re.sub('(?m)((?<=\n\n)\ *[A-Z][A-Z :!-]+\ *(?=\n))', r'=== \1 ===',mkstring)
    return mkstring

def login(baseurl,username,password):
    try:
        mech.open(baseurl + "UserSettings")
    except HTTPError, e:
        sys.exit("%d: %s" % (e.code, e.msg))

    # important: use correct form number in the HTML page
    mech.select_form(nr=1)
    mech["name"] = username
    mech["password"] = password
    try:
        mech.submit()
    except HTTPError, e:
        sys.exit("login failed: %d: %s" % (e.code, e.msg))

def makediff(old,new):
    import difflib, time
    if not _fulldiff:
        dm = difflib.unified_diff
    else:
        dm = difflib.ndiff
    old = re.sub('\015\012',r'\n',old).expandtabs(4)
    diff = dm(old.splitlines(1), new.splitlines(1))
    difftext = ''.join(diff)

    return difftext

def post(title,text):
    '''Browse to the Wiki, select the form, diff it, and post text'''
    title = re.sub("[^\w]|_","",title)
    try:
        resp = mech.open(baseurl+title+"/edit")
    except HTTPError, e:
        sys.exit("post failed: %d: %s" % (e.code, e.msg))

    mech.select_form(nr=1)
    try:
        oldbody = mech["body"]
    except:
        print 'No id="body" form. Figure this out first. cf. page text above.'
        for form in mech.forms():
                print form
        sys.exit("post failed while processing page: " + baseurl + title +'/edit')

    # manual override when text contains <!--protected-->
    if oldbody.find('<!--protected-->')>0:
        print 'skipped protected page ' + title
        return
    # generate diff and output the diffstat
    diff = makediff(oldbody,text)
    if diff:
        minus = str([l.startswith('-') for l in diff.splitlines() \
                if not l.startswith('---')].count(True))
        plus = str([l.startswith('+') for l in diff.splitlines() \
                if not l.startswith('+++')].count(True))
        if plus == '0' and  minus == '0':
            changes = ' (nc+)'
        else:
            changes = ' ('+plus+' insertions(+), '+minus+' deletions(-))'

        print "posting: " + title + changes
        if _debug: print diff
    else:
        print "posting: " + title + ' (nc)'

    mech["body"] = text
    try:
        mech.submit()
    except HTTPError, e:
        sys.exit("post failed: %d: %s" % (e.code, e.msg))

def postsinglefiles(files):
    ''' create pages for all files in there and link to the category'''
    for name in files:
        # single out some names
        head, basename = os.path.split(name)
        root =  os.path.join(os.path.basename(os.getcwd()))
        path =   os.path.normpath(os.path.join(root,head))
        funcname, ext = os.path.splitext(basename)
        category = os.path.basename(head)
        cattext = ''

        if not os.path.exists(name) or not ext=='.m':
            if _mexmode and ext == _mexext:
                mexhelpextract([funcname])
            elif _debug: print 'skipping ' + name
            continue
        if basename=='Contents.m' or basename=='contents.m':
            funcname = category
            path = os.path.dirname(path)
            if len(path) < 3:
                continue
            cattext = '{{category}}'

        # load and build the text for the page
        headline = "===[[" + funcname + "]]===\n"
        breadcrumb = "==[[" + re.sub("/","]] &#8250; [[", path) + "]]==\n\n"

        # read the .m-file and strip all the crap
        docstring = parse(name)

        # scrub the text real good, to get us some nice wiki formatting
        if docstring:
            body = beackern(docstring)
        else:
            body = 'Adrian, this function is not yet documented.\n\n\n MissingDocs'

        pathlinks = """
                    ""
                    <div class="code_header" style="text-align:right;">
                      <span style="float:left;">Path&nbsp;&nbsp;</span> <span class="counter">Retrieve <a href=
                      "https://raw.github.com/Psychtoolbox-3/Psychtoolbox-3/beta/%s">current version from GitHub</a> | View <a href=
                      "https://github.com/Psychtoolbox-3/Psychtoolbox-3/commits/beta/%s">changelog</a></span>
                    </div>
                    <div class="code">
                      <code>%s</code>
                    </div>
                    ""
                    """ % tuple([os.path.join(path,basename)]*3)

        text = headline \
                + breadcrumb \
                + body \
                + '\n\n\n' \
                + textwrap.dedent(pathlinks) \
                + '\n' \
                + cattext

        post(funcname,text)

def mexhelpextract(mexnames):
    #print 'processing mex files: ' + mexnames.__repr__()
    from ConfigParser import RawConfigParser as ConfigParser, Error as error
    for mexname in mexnames:
        # ConfigParser for the three elements per subfunctions written to tmpdir
        # [SubFunction]
        # usage: 'xyz'
        # help: 'xyz'
        # seealso: 'xyz'
        config = ConfigParser({'usage':[], 'help':[], 'seealso':[]})
        # assemble command line for matlab
        matlabcmd = 'addpath(\'%s\');%s(\'%s\',\'%s\'); exit' % \
            (_tmpdir, \
             os.path.splitext(os.path.basename(_mexscript))[0], \
             mexname, \
             _tmpdir)
        cmd = 'matlab -nojvm -nodisplay -r "%s" > /dev/null' % matlabcmd
        # and execute matlab w/ the temporary script we wrote earlier
        try:
            print 'running MATLAB for %s in %s' % (mexname,_tmpdir)
            p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, close_fds=True)
            stderr = p.communicate()[1]
            if stderr: print stderr
        except:
            print 'could not dump help for %s into %s.' % (mexname,_tmpdir)

        cfgfile = config.read(os.path.join(_tmpdir,mexname))
        if cfgfile == []:
            print "skipping " + mexname + " (no output)"
            continue
        subfunctions = config.sections()
        print 'processing subfunctions: ' + subfunctions.__repr__()
        for subfunction in subfunctions:
            # read in the strings for this subfunction
            usage = config.get(subfunction,'usage')
            help = config.get(subfunction,'help')
            seealso = config.get(subfunction,'seealso')

            headline = '===[['+subfunction+' '+mexname+'(\''+subfunction+'\')]]===\n'
            breadcrumb = "==[[Psychtoolbox]] &#8250; [[" \
                                + mexname + "]].{mex*,dll} subfunction==\n\n"

            # scrub the text for main text only
            body = beackern(help)

            docstring = '' \
                    + '%%(matlab;Usage)' \
                    + usage \
                    + '%%\n' \
                    + body \
                    + '\n\n'
            if seealso:
                docstring = docstring + '<<=====See also:=====\n' + seealso + '<<'

            text =  '""' + headline \
                    + breadcrumb \
                    + docstring + '""'

            # retrieve old body text, to update or concatenate with synonymous subfunctions
            #
            # browse the page
            title = re.sub("[^\w]|_","",subfunction)
            try:
                resp = mech.open(baseurl+title+"/edit")
            except HTTPError, e:
                sys.exit("retrieving old text during posting of this mex function failed: %d: %s" % (e.code, e.msg))
            # get text from the edit form
            mech.select_form(nr=1)
            try:
                oldbody = mech["body"]
            except:
                print 'No id="body" form. Figure this out first. cf. page text above.'
                for form in mech.forms():
                        print form
                sys.exit("retrieving old body text failed while processing page: " + baseurl + title +'/edit')

            # parse embedded structuring HTML tags in the wiki text
            soup = BeautifulSoup(oldbody)

            # check if the subfunction is already present, by CSS 'class' and 'id'
            subfct = soup.find('div', {'class' : "subfct", 'id' : mexname})
            if subfct:
                # replace the text of the container DIV
                subfct.contents[0].replaceWith(text)
            else:
                # contruct new DIV to hold the text
                subfctDIV = Tag(soup, "div")
                subfctDIV['class'] = 'subfct'
                subfctDIV['id'] = mexname
                subfctDIV.insert(0,NavigableString(text))

                # insert the new div
                soup.insert(len(soup),subfctDIV)

            # Now scoop the good well-formed divs out of the soup
            divs = soup('div', {'class' : "subfct"})

            # and drop them into fresh yummy cheese soup
            cheesesoup = BeautifulSoup()

            # drop good divs into the soup, one by one
            for div in divs:
                # remove the unneeded style attribute, we finally
                # have this stuff defined in the ptbdocs.css now.
                del(div['style'])
                # escape the HTML tags for wiki parser
                cheesesoup.append(NavigableString('\n""'))
                cheesesoup.append(div)
                cheesesoup.append(NavigableString('""\n'))

            post(subfunction,cheesesoup.renderContents())

def recursivewalk(rootfolder):
    '''traverse the root directory and post all .m-files'''
    for root, dirs, files in os.walk(rootfolder):

        if re.search(os.sep+'\.svn|'+os.sep+'private',root):
                continue
        print 'Entering folder ' + root

        # .m-files are processed with postsinglefiles
        mfiles = [os.path.join(root,f) for f in files \
                if f[-2:] in ('.m','.M')]
        postsinglefiles(mfiles)

        # mex files are processed via matlab
        if _mexmode:
            mexnames = [os.path.splitext(f)[0] for f in files \
                    if f[-len(_mexext):].lower() == _mexext]
            mexhelpextract(mexnames)

    print "Exiting: Done with this tree."
    sys.exit(0)

def usage():
    print __doc__

def main(argv):

    try:
        opts, args = getopt.getopt(argv, "hu:p:U:frmv", \
            ["help", "username=", "password=", "url=", "full-diff", "mexmode", "recursive"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-u", "--username"):
            global username
            username = arg
        elif opt in ("-p", "--password"):
            global password
            password = arg
        elif opt in ("-U", "--url"):
            global baseurl
            baseurl = arg
        elif opt in ("-f", "--full-diff"):
            global _fulldiff
            _fulldiff = 1
        elif opt in ("-r", "--recursive"):
            global _recursive
            _recursive = 1
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
                return;
            end

            try
                if isempty(subfunctions)
                    return;
                end

                if ~iscell(subfunctions)
                    return;
                end
            catch
                return;
            end

            fid = fopen(fullfile(tmpdir,mexname),'wt');
            if fid == -1
                return;
            end

            try
                for i=1:size(subfunctions,2)
                    fprintf(fid,'[%s]\\n\\n',subfunctions{i});
                    docs = eval([mexname '(''DescribeModulefunctionshelper'',1,subfunctions{i})']);
                    fprintf(fid,'usage: %s\\n\\n',docs{1});
                    fprintf(fid,'help: %s\\n\\n',regexprep(docs{2},'\\n','\\n  '));
                    fprintf(fid,'seealso: %s\\n\\n',docs{3});
                end
            catch
                % Nothing to do.
            end

            fclose(fid);
            return
            '''
            fid.write(textwrap.dedent(script))
            fid.close()

        elif opt == '-v':
            global _debug
            _debug = 1

    if len(password) == 0:
        usage()
        print 'Error: No default password set. Please specify one with -p MyPassword.\n\n'
        sys.exit(1)
    else:
        login(baseurl,username,password)

    if args:
        if _recursive:
            if len(args) == 1 and os.path.isdir(args[0]):
                recursivewalk(args[0])
            else:
                usage()
                print 'Error: Recursive mode works with one directory only\n\n'
                sys.exit(1)

        else:
            postsinglefiles(args)
    else:
        usage()
        print 'Error: No files specified to submit to Wiki\n\n'
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])
