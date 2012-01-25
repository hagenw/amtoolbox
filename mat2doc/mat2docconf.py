
# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

import localconf

conf=ConfType()

conf.octexec=localconf.octexec
conf.matexec=localconf.matexec
conf.plotengine=localconf.plotengine
   
conf.root=localconf.basepath+'amtoolbox/'

conf.bibfile=conf.root+'mat2doc/project'

conf.workdir=localconf.userdir+'publish/'

conf.ignorepars=['a','order','type']

conf.copyrightplate=localconf.basepath+'amtscripts/copyrightplate'

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
                  'experiments/Contents',
                  'signals/Contents',
                  'demos/Contents']

texcontentsfiles=allcontentsfiles


# ------------------------------------------
# Configuration of PHP for Sourceforge
# ------------------------------------------

php=HTMLConf()

php.basetype='html'

php.subdir='amthtml/'

php.indexfiles=allcontentsfiles

# This is the full path, used for php-including files.
php.docroot='/home/groups/a/am/amtoolbox/htdocs/doc/'

# This is the usual web-server root for "<a href=" ... > tags.
php.htmlroot='/doc/'
    
php.hb='<H2>'

php.he='</H2>'

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

   <?php ini_set("include_path",".:/home/groups/a/am/amtoolbox/htdocs/"); ?>
   <?php require("topnav.php");?>

"""

php.foot="""<?php require("footer.php");?>
</body>
</html>"""

php.dryrun=1


# ------------------------------------------
# Configuration of LaTeX
# ------------------------------------------

tex=TexConf()

tex.basetype='tex'

tex.subdir='amtpdf/'

tex.texfile='amt'

tex.indexfiles=texcontentsfiles
    
tex.head="""\documentclass{amsart}
\usepackage{ae}
\usepackage{aecompl}
\usepackage[T1]{fontenc}
\usepackage[latin1]{inputenc}
\usepackage{graphicx}
\setlength\parskip{\medskipamount}
\setlength\parindent{0pt}
\makeatletter

%\usepackage{babel}
\usepackage{amssymb}
\makeatother
\\begin{document}
\\title{AMT Reference manual}
\\author{The AMT team}

\\maketitle
\\tableofcontents{}

"""
tex.foot='\\bibliographystyle{abbrv}\n'+ \
         '\\bibliography{'+conf.bibfile+'}\n'+ \
"""\end{document}
"""

tex.dryrun=1
tex.dooutput=1

# ------------------------------------------
# Configuration of Matlab
# ------------------------------------------

mat=ConfType()

mat.basetype='mat'

# ------------------------------------------
# Configuration of Verification system
# ------------------------------------------

verify=ConfType()

verify.basetype='verify'

verify.targets=['AUTHOR','TESTING','REFERENCE']

verify.notappears=['FIXME','BUG','XXX']

verify.ignore=["demo_","comp_","assert_","Contents.m","init.m"]

