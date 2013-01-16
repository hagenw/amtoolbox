
# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

import localconf

conf=ConfType()

conf.otherrefs=['ltfat.txt']

def copyrightfun():
    f=file(localconf.projects['amtoolbox']+'amtoolbox_version');
    versionstring=f.read()[:-1]
    f.close
        
    f=file(localconf.projects['ltfat']+'mat2doc/copyrightplate')
    buf=f.readlines()
    f.close

    copyright=[u'Copyright (C) 2013 Peter L. S\xf8ndergaard.\n',
               u'This file is part of AMToolbox version '+versionstring+'\n']
    copyright.extend(buf)
    
    return copyright

conf.copyright=copyrightfun

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
php.codedir=localconf.outputdir+'amtoolbox-mat'+os.sep

# ------------------------------------------
# Configuration of LaTeX
# ------------------------------------------

tex=TexConf()

tex.indexfiles=allcontentsfiles
tex.urlbase='http://amtoolbox.sourceforge.net/doc/'
tex.codedir=localconf.outputdir+'amtoolbox-mat'+os.sep
    
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

