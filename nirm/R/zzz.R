######################################################################
#
# zzz.R
# Author: Chris Marcum cmarcum@uci.edu
# Last Modified 01/12/12
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/nirm package
#
# .First.lib or .onLoad is run when the package is loaded with library(nirm)
#
######################################################################

.onLoad <- function(lib, pkg){
    abt<-packageDescription('nirm')
    packageStartupMessage(
    cat("\n\t",abt$Title,"\n"),
    cat("Version: ",abt$Version,"\n"),
    cat("Compiled:",abt$Built,"\n\n"),
    cat(paste("copyright (c) 2012, Raffaele Vardavas and Christopher S. Marcum, RAND\n",sep="")),
    cat('Type help(package="nirm") to get started.\n'))
}
