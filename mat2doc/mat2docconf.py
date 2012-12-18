
# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

import localconf

conf=ConfType()

conf.otherrefs=['ltfat.txt']

conf.urlbase='http://ltfat.sourceforge.net/doc/'

def mycopyrightfun(self):
    vf=file(self.root+'amtoolbox_version');
    v=vf.readline()
    vf.close
    
    f=file(self.root+'mat2doc/copyrightplate')
    buf=f.readlines()
    f.close

    copyright=[u'Copyright (C) 2012 Peter L. S\xf8ndergaard.\n','This file is part of AMT version '+v]
    copyright.extend(buf)
    
    return copyright

conf.copyright=mycopyrightfun

allcontentsfiles=['Contents',
                  'general/Contents',
                  'filters/Contents',
                  'modelstages/Contents',
                  'monaural/Contents',
                  'binaural/Contents',
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

