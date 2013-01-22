# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

# When writing this file, the following variables are available
#
#   projectname  - The name of your project, OK, you already know this
#   projectdir   - The directory where your project is located
#   outputdir    - The output directory

from mat2doc import *

f=file(projectdir+'amtoolbox_version')
versionstring=f.read()[:-1]
f.close

conf=ConfType()

#conf.otherrefs=['ltfat.txt']

f=file(localconf.projects['amtoolbox']+'amtoolbox_version');
versionstring=f.read()[:-1]
f.close

f=file(localconf.projects['ltfat']+'mat2doc/copyrightplate')
buf=f.readlines()
f.close

copyright=[u'Copyright (C) 2013 Peter L. S\xf8ndergaard.\n',
           u'This file is part of AMToolbox version '+versionstring+'\n']
copyright.extend(buf)

conf=ConfType()

conf.copyright=copyright

allcontentsfiles=['Contents',
                  'general/Contents',
                  'filters/Contents',
                  'modelstages/Contents',
                  'monaural/Contents',
                  'binaural/Contents',
                  'speech/Contents',
                  'humandata/Contents',
                  'hrtf/Contents',
                  'experiments/Contents',
                  'signals/Contents',
                  'demos/Contents']

# ------------------------------------------
# Configuration of PHP for Sourceforge
# ------------------------------------------

php=PhpConf()

php.indexfiles=allcontentsfiles
php.includedir='../include/'
php.urlbase='/doc/'
php.codedir=outputdir+'amtoolbox-mat'+os.sep

# ------------------------------------------
# Configuration of LaTeX
# ------------------------------------------

tex=TexConf()

tex.indexfiles=allcontentsfiles
tex.urlbase='http://amtoolbox.sourceforge.net/doc/'
tex.codedir=outputdir+'amtoolbox-mat'+os.sep
    
# ------------------------------------------
# Configuration of Matlab
# ------------------------------------------

mat=MatConf()
mat.urlbase='http://amtoolbox.sourceforge.net/doc/'

# ------------------------------------------
# Configuration of Verification system
# ------------------------------------------

verify=ConfType()

verify.basetype='verify'

verify.targets=['AUTHOR','TESTING','REFERENCE']

verify.notappears=['FIXME','BUG','XXX']

verify.ignore=["demo_","comp_","assert_","Contents.m","init.m"]

