#!/usr/bin/python

import sys,os

cwd=os.getcwd()+'/'

sys.path.append(cwd)

from localconf import *
from datesuffix import *

tbpath  =basepath+'amtoolbox/'
tbwww   =basepath+'amtwww/'
m2dfile =tbpath+'mat2doc/mat2docconf.py'
host='soender,amtoolbox@web.sourceforge.net'
www='/home/groups/a/am/amtoolbox/htdocs/'
notesdir=basepath+'amtnotes/'
notehtml=cwd+'amtnoteshtml/'
noteswww=www+'notes/'

sys.path.append(basepath+'mat2doc/')


import printdoc, notes


def lyx2pdf(workdir, fname):
    os.system('cd '+workdir+'; lyx '+fname+'.lyx --export latex')
    os.system('cd '+workdir+'; pdflatex '+fname)
    os.system('cd '+workdir+'; bibtex '+fname)
    os.system('cd '+workdir+'; pdflatex '+fname)
    os.system('cd '+workdir+'; pdflatex '+fname)

    os.system('scp '+workdir+'/'+fname+'.pdf '+host+':'+www)
    

if len(sys.argv)==1:
    todo=[]
else:
    todo=sys.argv[1:]

if 'verify' in todo:
    printdoc.printdoc([m2dfile,'verify'])

if 'stagemat' in todo:
    publishmat=cwd+'amtoolbox/'
    printdoc.git_stageexport(tbpath,publishmat)
    printdoc.print_mat(m2dfile,publishmat)

if 'releasemat' in todo:
    publishmat=cwd+'amtoolbox/'
    printdoc.git_repoexport(tbpath,'master','amtoolbox',cwd)
    printdoc.print_mat(m2dfile,publishmat)

    # Remove unwanted files
    os.system('rm -rf '+publishmat+'testing')
    os.system('rm -rf '+publishmat+'reference')


    f=file(tbpath+'amt_version')
    versionstring=f.read()[:-1]
    f.close()
    
    fname=cwd+'amtoolbox-'+versionstring
    os.system('rm '+fname+'.zip')

    printdoc.dos2unix(cwd+'amtoolbox')
    os.system('gtar zcvf '+fname+'.tgz amtoolbox/')

    printdoc.unix2dos(cwd+'amtoolbox')
    os.system('zip -r '+fname+'.zip amtoolbox/')

if 'releasebranch' in todo:
    bname=sys.argv[2]

    publishmat=cwd+'amtoolbox/'
    printdoc.git_repoexport(tbpath,bname,'amtoolbox',cwd)
    printdoc.print_mat(m2dfile,publishmat)

    # Remove unwanted files
    os.system('rm -rf '+publishmat+'testing')
    os.system('rm -rf '+publishmat+'reference')

    f=file(tbpath+'amt_version')
    versionstring=f.read()[:-1]
    f.close()
    
    fname=cwd+'amtoolbox-'+bname+versionstring
    os.system('rm '+fname+'.zip')

    printdoc.dos2unix(cwd+'amtoolbox')
    os.system('gtar zcvf '+fname+'.tgz amtoolbox/')

    printdoc.unix2dos(cwd+'amtoolbox')
    os.system('zip -r '+fname+'.zip amtoolbox/')
        
if 'pdf' in todo:
    printdoc.printdoc([m2dfile,'tex'])
    os.system('cd amtpdf; pdflatex amt.tex')
    os.system('cd amtpdf; bibtex amt')
    os.system('cd amtpdf; pdflatex amt.tex')
    
    os.system('scp toolboxref/toolboxref.pdf '+host+':'+www+'doc/amt.pdf')

if 'stagewww' in todo:
    publishwww=cwd+'amtwww/'
    printdoc.git_stageexport(tbwww,publishwww)
    os.system('rsync -av '+publishwww+' '+host+':'+www);

if 'releasewww' in todo:
    publishwww=cwd+'amtwww/'
    printdoc.git_repoexport(tbwww,'master','amtwww',cwd)
    os.system('rsync -av '+publishwww+' '+host+':'+www);

if 'doc' in todo:
    os.system('cd '+tbpath+'doc/; make install INSTALLPATH='+cwd+'amtoolbox/doc/')
    os.system('rsync -av amtoolbox/doc/* '+host+':'+www+'doc/');

if 'php' in todo:
    printdoc.printdoc([m2dfile,'php'])
    os.system('cp '+tbwww+'doc/index.php amthtml/')
    os.system('rsync -av amthtml/ '+host+':'+www+'doc/')

if 'notesmake' in todo:
    notes=notes.getnotenumbers(notesdir)

    notes = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes)

    for notenumber in notes:
        print 'Trying to make LTFAT note '+notenumber
        os.system('cd '+notesdir+notenumber+'; make')

if 'notestexclean' in todo:
    notes=notes.getnotenumbers(notesdir)

    notes = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes)

    for notenumber in notes:
        os.system('cd '+notesdir+notenumber+'; make texclean')

if 'noteshtml' in todo:
    printdoc.printnoteshtml('amtnote',notesdir,notehtml)
        
    os.system('rsync -av '+notehtml+' '+host+':'+noteswww);
