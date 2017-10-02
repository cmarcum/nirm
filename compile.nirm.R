#This is the nirm package compiler
# Author: Christopher Steven Marcum
# Created on: January 2012
# Last Modified: 30 April 2012
#The RiskPerceptionNetworks/Tools/nirm/nirm.master.R file is the master file of the nirm package.
# nirm.master.R contains special hooks that should not be removed or duplicated. 
# Here is a list of those hooks:
#  order) Name in nirm, name here
#  1) FUNCTIONS HERE, FUN
#  2) BEH MODELS HERE, BEH
#  3) PTN FILE HERE, PTN
#  4) DYNAMICS FILE HERE, DYN
#  5) EOF, 

##Part 1: Collecting Functions
#Locate and store the paths to Raff's functions, excepting that pesky library.R file:
raffs.lib<-paste("../../Flu_R21/Flu model on the network/Library/",dir("../../Flu_R21/Flu model on the network/Library"),sep="")
raffs.lib<-raffs.lib[-which(grepl("library.R",raffs.lib))]

#source the files (stripping any commments, whitespaces, and illegal formatting), 
# drop the paths, then dump them to an intuitively named file. Finally, clean up the workspace.
sapply(raffs.lib,function(x) source(x, keep.source = FALSE))
rm(raffs.lib)
dump(ls(all = TRUE), file = "nirm.funs.R")
rm(list=ls())

##Part 2: Writing Functions and other scripts to master file nirm.R
#Now for some real blackmagic:
mfile<-"nirm.master.R" #The master file
ofile<-"nirm.R" #The output file to be symbolically linked to the R package Tree

system(paste("cp", mfile,ofile,sep=" "))

#The FUN hook
ifile<-"nirm.funs.R"
fun.hook<-as.numeric(system(paste("grep 'FUNCTIONS HERE' ", ofile," -n | cut -f1 --delim=':'",sep=""),intern=TRUE))+1
system(paste("sed -i '",fun.hook,"r ", ifile,"' ", ofile,sep=""))

#The BEH hook
ifile<-"../../Flu_R21/Flu model on the network/Behavioral.Models.R"
beh.hook<-as.numeric(system(paste("grep 'BEH MODELS HERE' ", ofile," -n | cut -f1 --delim=':'",sep=""),intern=TRUE))+1
system(paste("sed -i '",beh.hook,"r ", ifile,"' ", ofile,sep=""))

#The PTN hook
ifile<-"../../Flu_R21/Flu model on the network/Construct_PTN.R"
ptn.hook<-as.numeric(system(paste("grep 'PTN FILE HERE' ", ofile," -n | cut -f1 --delim=':'",sep=""),intern=TRUE))+1
system(paste("sed -i '",ptn.hook,"r ", ifile,"' ", ofile,sep=""))

#The DYN hook
ifile<-"../../Flu_R21/Flu model on the network/Evolve.Dynamics.R"
dyn.hook<-as.numeric(system(paste("grep 'DYNAMICS FILE HERE' ", ofile," -n | cut -f1 --delim=':'",sep=""),intern=TRUE))+1
system(paste("sed -i '",dyn.hook,"r ", ifile,"' ", ofile,sep=""))

#Clean up the workspace
rm(list=ls())

##Part 3: Compiling R Package --- Note, nirm/ChangeLog will always need to be updated manually.
system("mv nirm.R nirm/R/nirm.R")
system("rm nirm.funs.R")
system("rm *~")
system("R CMD build nirm/")
pkg.name<-dir(pattern=".tar.gz")
system(paste("mv",pkg.name,"nirm.tags/",sep=" "))
q("no")
