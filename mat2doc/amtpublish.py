#!/usr/bin/python

import sys,os,os.path
cwd=os.getcwd()+'/'

# ------- Configuration parameters -------------

projectname='amtoolbox'

# Configure HTML placement at remote server
host='soender,amtoolbox@web.sourceforge.net'
noteswww='/home/project-web/amtoolbox/htdocs//'

tbwww='/home/peter/nw/amtwww'
notesdir='~/nw/amtnotes'
noteshtml='~/publish/amtnoteshtml'
noteswww='/home/project-web/amtoolbox/htdocs/notes/'

notesdir=os.path.expanduser(notesdir)+os.sep
noteshtml=os.path.expanduser(noteshtml)+os.sep

# ------- do not edit below this line ----------

# Import the localconf file, it should be place in the same directory
# as this function is called from.

sys.path.append(cwd)
import localconf

sys.path.append(localconf.mat2docdir)

import printdoc, notes

todo=sys.argv[1]

if 'notesmake' in todo:
    notes=notes.getnotenumbers(notesdir)

    notes = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes)

    for notenumber in notes:
        print 'Trying to make AMTOOLBOX note '+notenumber
        os.system('cd '+notesdir+notenumber+'; make')

if 'notestexclean' in todo:
    notes=notes.getnotenumbers(notesdir)

    notes = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes)

    for notenumber in notes:
        os.system('cd '+notesdir+notenumber+'; make texclean')

if 'noteshtml' in todo:
    notes.printnoteshtml('amtnote',notesdir,noteshtml)
        
    os.system('rsync -av '+noteshtml+' '+host+':'+noteswww);
