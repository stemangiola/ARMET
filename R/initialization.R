#-----------------------------#
#-----by Stefano Mangiola-----#
#-----------Nov 2015----------#
#-----------------------------#

if("package:org.Hs.eg.db" %in% search()) detach("package:org.Hs.eg.db", unload=TRUE, force=TRUE)
if("package:hgu133plus2.db" %in% search()) detach("package:hgu133plus2.db", unload=TRUE, force=TRUE)
if("package:hgu133a.db" %in% search()) detach("package:hgu133a.db", unload=TRUE, force=TRUE)
if("package:hgu95av2.db" %in% search()) detach("package:hgu95av2.db", unload=TRUE, force=TRUE)
if("package:illuminaHumanv4.db" %in% search()) detach("package:illuminaHumanv4.db", unload=TRUE, force=TRUE)
