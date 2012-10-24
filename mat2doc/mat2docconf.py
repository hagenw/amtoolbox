
# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

import localconf

conf=ConfType()

conf.otherrefs=['ltfat.txt']

conf.urlbase='http://ltfat.sourceforge.net/doc/'

def mycopyrightfun(self):
    vf=file(self.root+'amt_version');
    v=vf.readline()
    vf.close
    
    f=file(self.copyrightplate)
    buf=f.readlines()
    f.close

    copyright=['Copyright (C) 2011 Peter L. Soendergaard.\n','This file is part of AMT version '+v]
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

texcontentsfiles=allcontentsfiles


# ------------------------------------------
# Configuration of PHP for Sourceforge
# ------------------------------------------

php=phpConf()

php.basetype='php'

php.indexfiles=allcontentsfiles

# This is the full path, used for php-including files.
php.docroot='/home/project-web/amtoolbox/htdocs/doc/'

# This is the usual web-server root for "<a href=" ... > tags.
php.fext='.php'

php.head="""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN"><html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Auditory Modelling Toolbox</title>
<link rel="stylesheet" href="/amtoolbox.css" type="text/css">
<script type="text/javascript"
   src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
</head>
<body>

   <?php ini_set("include_path",".:/home/project-web/amtoolbox/htdocs/"); ?>
   <?php include("topnav.php");?>

"""

php.foot="""<?php include("footer.php");?>
</body>
</html>"""

# ------------------------------------------
# Configuration of LaTeX
# ------------------------------------------

tex=TexConf()

tex.basetype='tex'

tex.subdir='amtpdf/'

tex.texfile='amt'

tex.indexfiles=texcontentsfiles
    
# ------------------------------------------
# Configuration of Matlab
# ------------------------------------------

mat=ConfType()

mat.basetype='mat'

mat.subdir='amtoolbox/'

# ------------------------------------------
# Configuration of Verification system
# ------------------------------------------

verify=ConfType()

verify.basetype='verify'

verify.targets=['AUTHOR','TESTING','REFERENCE']

verify.notappears=['FIXME','BUG','XXX']

verify.ignore=["demo_","comp_","assert_","Contents.m","init.m"]

