#!/usr/bin/python

import sys,os
cwd=os.getcwd()+'/'

# ------- Configuration parameters -------------

projectname='amtoolbox'

# Configure HTML placement at remote server
host='soender,amtoolbox@web.sourceforge.net'
www='/home/project-web/amtoolbox/htdocs//'

tbwww='/home/peter/nw/amtwww'

# ------- do not edit below this line ----------

# Import the localconf file, it should be place in the same directory
# as this function is called from.

sys.path.append(cwd)

import localconf

sys.path.append(localconf.mat2docdir)

# Get the data from localconf
project=eval('localconf.'+projectname)
conffile=project['dir']+'/mat2doc/mat2docconf.py'
filesdir=localconf.filesdir

f=file(project['dir']+projectname+'_version')
versionstring=f.read()[:-1]
f.close()

import printdoc

todo=sys.argv[1]

if 'verify' in todo:
    printdoc.printdoc(projectname,'verify')

# Release for other developers to download
if 'develmat' in todo:
    #printdoc.git_stageexport(project['dir'],project['mat'])
    os.system('svn export --force '+project['dir']+' '+project['mat'])
    printdoc.printdoc(projectname,'mat')

    fname=filesdir+projectname+'-devel-'+versionstring
    os.system('rm '+fname+'.zip')

    # Create the Unix src package
    os.system('tar zcvf '+fname+'.tgz '+projectname+'/')

    # Create the Windows src package
    os.system('rm '+fname+'.zip')
    printdoc.unix2dos(filesdir+projectname)
    os.system('zip -r '+fname+'.zip '+projectname+'/')

# Release for users to download
if 'releasemat' in todo:
    printdoc.git_repoexport(project['dir'],'master',projectname,filesdir)
    #os.system('svn export --force '+project['dir']+' '+project['mat'])
    printdoc.printdoc(projectname,'mat')
    
    # Remove unwanted files
    os.system('rm -rf '+project['mat']+'testing')
    os.system('rm -rf '+project['mat']+'reference')

    fname=filesdir+projectname+'-'+versionstring

    # Create the Unix src package
    os.system('tar zcvf '+fname+'.tgz '+projectname+'/')

    # Create the Windows src package
    os.system('rm '+fname+'.zip')
    printdoc.unix2dos(filesdir+projectname)
    os.system('zip -r '+fname+'.zip '+projectname+'/')

# Quick staging to test for errors in the files
if 'stagemat' in todo:
    printdoc.git_stageexport(project['dir'],project['mat'])
    printdoc.printdoc(projectname,'mat')
    
if 'tex' in todo:
    printdoc.printdoc(projectname,'tex')

if todo=='php':
    printdoc.printdoc(projectname,'php')
    s='rsync -av '+project['php']+' '+host+':'+www+'doc/'
    os.system(s)    

if todo=='phplocal' in todo:
    printdoc.printdoc(projectname,'php')

if todo=='phprebuild' in todo:
    printdoc.printdoc(projectname,'php','rebuild')

if 'verify' in todo or todo==[]:
    printdoc.printdoc([m2dfile,'verify'])





#  --------- old file -------------------------------



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
    os.system('tar zcvf '+fname+'.tgz amtoolbox/')

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
    printdoc.printnoteshtml('amtnote',notesdir,notehtml)
        
    os.system('rsync -av '+notehtml+' '+host+':'+noteswww);
