####################################################################################
####################################################################################
##                        Network Inductive Reasoning Model                       ##
##                                                                                ##
## Model Code Written by Raffael Vardavas (vardavas@rand.org)                     ##
## R-Package Authored by Christopher Steven Marcum (cmarcum@rand.org)             ##
## Package Maintained by Christopher Steven Marcum (cmarcum@rand.org)             ##
## Created on 2 January 2012                                                      ##
## Last Modified on 7 May 2012                                                    ##
##                                                                                ##
## Part of the nirm package for R                                                 ##
## Compile with R CMD BATCH compile.nirm.R                                        ##
## WARNING: Contains Hooks. Read compile.nirm.R header before editing.            ##
##                                                                                ##
####################################################################################
####################################################################################

nirm<-function(network.model=c("Erdos","Watts","NW.Watts","Barabasi","Empirical"),
               N=10^3,M=20,behavioral.model="Simple_Personal",enet=NULL,usbm=NULL,
               hyper.params=list(model=c(R0=2,gamma=1/3,iip=.02,eff=1,tfs=300,s.mem=.7),vattr=cbind(V=rep(0,N),w=rep(0.3,N),inf.nn=rep(0,N),vacc.nn=rep(0,N))),
               added.RMCN=FALSE, RMCN.R0=0, output.dir="Temp/",save.model=TRUE,verbose=FALSE,seed=NULL){
   #First, a bunch of user input checks:
   if(!any(network.model%in%c("Erdos","Watts","NW.Watts","Barabasi","Empirical"))){stop("Invalid Network Model. See help(\"nirm\") for details.")}
   network.model<-network.model[1]
   if(network.model=="Empirical" & is.null(enet)){stop("You say you want to use an empirical network, yet your enet object is NULL.")}
   if(network.model!="Empirical" & !is.null(enet)){stop("You say you want to use a theoretical network, yet your enet object is non-NULL.")}
   if(added.RMCN & RMCN.R0==0){warning("You have not specified an R0 for the RMCN.")} 
   if(output.dir==""){warning("You have not specified an output directory. Using present working directory; I hope you have enough free-memory here.")} 
   #Now all output files get a unique ID: 
   mxid<-paste("_MX",substr(strsplit(as.character(runif(1)),"\\.")[[1]][2],1,8),sep="")
   
   ### Set Output Dir, 
   if(!file.exists(output.dir)){dir.create(output.dir)}
   data.dir   <- paste(output.dir,"Data/",sep="")
   if(!file.exists(data.dir)){dir.create(data.dir)}
   analysis.dir   <- paste(output.dir,"Analysis/",sep="")
   if(!file.exists(analysis.dir)){dir.create(analysis.dir)}
   if(!file.exists(paste(output.dir,"Output",sep=""))){dir.create(paste(output.dir,"Output",sep=""))}
   output.networks.dir <- paste(output.dir,"Output/Networks/",sep="")
   if(!file.exists(output.networks.dir)){dir.create(output.networks.dir)}
   output.dynamics.dir <- paste(output.dir,"Output/Dynamics/",sep="")
   if(!file.exists(output.dynamics.dir)){dir.create(output.dynamics.dir)}
   figure.dir   <- paste(output.dir,"Figure/",sep="")
   if(!file.exists(figure.dir)){dir.create(figure.dir)}      

   #This flag is for future functionality:
   Do.Comp.Percolation.Analysis<-FALSE
   model.specs<-NULL #Raff added this at the bottom of the setup file.
   #Set the top environment
   nirm.env<-environment()   

   #Set seed if necessary
   if(!is.null(seed)){set.seed(seed)}

   #initialize params and placeholders
   ave.degree<-M
   attributes.data<-as.data.frame(hyper.params$vattr)
   #We no longer use betas in the program. Commented out in Raff's original source. - CSM & RV 30 April 12
   #big.beta.pref<-hyper.params$model["b0"]   
   #big.beta.pop<-hyper.params$model["B0"]
   gamma<-hyper.params$model["gamma"]
   R0.pop<-as.numeric(added.RMCN*RMCN.R0)
   R0.pref<-hyper.params$model["R0"]
   #We do not use R0.beta.pop- CSM & RV 30 April 12
   #R0.beta.pop<-big.beta.pop/gamma
   model<-network.model[1]

   initial.inf.prop<-hyper.params$model["iip"]
   eff<-hyper.params$model["eff"]
   total.flu.seasons<-hyper.params$model["tfs"]
   s.mem<-hyper.params$model["s.mem"]

   names.outcomes <- c("VaccInf","Vacc","Inf","UnInf")
   names.evaluation.levels <- c("Personal",
                             "Local.Inf","Local.Vacc",
                             "Global.Inf","Global.Vacc")

   behavioral.weights <- matrix(0,nrow=length(names.outcomes),
                             ncol=length(names.evaluation.levels), dimnames=
                             list(names.outcomes,names.evaluation.levels))
   
   Stiffness<-rep(0,length(names.outcomes))  
   names(Stiffness)<-names.outcomes  
   Personal.Delta<-rep(0,length(names.outcomes))  
   names(Personal.Delta)<-names.outcomes
   Behavioral.Model<-behavioral.model
   if(!is.null(usbm)){Behavioral.Model<-"USBM"}
   Pop.Samp <- list(NULL)
   flu.season<-1


   network.file.name <- paste(model,"_N",as.character(log(N,10)),
                           "_avek_",as.character(ave.degree),mxid,
                           ".RData",sep="")

    ### Extra label in file name
    extra.file.name.label <- mxid

   if(is.null(extra.file.name.label)){
     output.file.name <- paste(model,"_N",as.character(log(N,10)),
                               "_",Behavioral.Model,"_avek_",
                               as.character(ave.degree),".RData",sep="")
   }else{  
     output.file.name <- paste(model,"_N",as.character(log(N,10)),
                            "_",Behavioral.Model,"_avek_",as.character(ave.degree),
                            extra.file.name.label,".RData",sep="")
   }
   
   #Reset params if user supplied network is valid
   if(network.model=="Empirical"){
      if(is.directed(enet)){stop("This appears to be a directed network. Please supply an undirected one.")}
      N<-vcount(enet)
      M<-ecount(enet) #possibly used at a later date
      ave.degree<-mean(degree(enet))
      ave.degree.f<-as.character(ceiling(ave.degree)) #Placeholder for file naming precedent.
      g.PTN<-enet
      g.PTN <- simplify(g.PTN)

      g.PTN <- Initialize.All.Vertex.Attributes(g.PTN, attributes.data)
      g <- g.PTN
      model<-network.model

      #list of nearest neighbours (nn) on the static Social Network. 
      nn.list.PTN <- get.nn.list.of.SN(g.PTN)
      N.nn <- degree(g.PTN)
      phi<-NA
      #Save the PTN in an R object - for post analysis.
      network.file.name <- paste(model,"_N",as.character(log(N,10)),
                           "_avek_",ave.degree.f,mxid,
                           ".RData",sep="")
    
      ### Extra label in file name
      extra.file.name.label <- mxid

      Theta.PTN <- R0.pref * (mean(N.nn^2) - mean(N.nn))/mean(N.nn)^2
      pc.PTN <- 1 - 1/Theta.PTN
      percolation.analisis <- list(Theta.PTN = Theta.PTN, pc.PTN = pc.PTN)

      if(is.null(extra.file.name.label)){
        output.file.name <- paste(model,"_N",as.character(log(N,10)),
                               "_",Behavioral.Model,"_avek_",
                               ave.degree.f,".RData",sep="")
      }else{  
        output.file.name <- paste(model,"_N",as.character(log(N,10)),
                            "_",Behavioral.Model,"_avek_",ave.degree.f,
                            extra.file.name.label,".RData",sep="")
      }
     if(save.model){
        save(g.PTN,N,phi,N.nn,nn.list.PTN,percolation.analisis, file = paste(output.networks.dir,network.file.name,sep=""), compress="xz")}

   }

   #Finalize vector of initial model parameterize to pass to other methods 
   #Changed 30 April 2012 to new parameters (sans betas) -CSM & RV
   params.of.model <- c(N=N,ave.degree=ave.degree, s.mem=s.mem, initial.inf.prop=initial.inf.prop,R0.pref=R0.pref,
                     R0.pop=R0.pop,gamma=gamma,total.flu.seasons=total.flu.seasons,
                     eff=eff,added.RMCN=added.RMCN)

   #Methods for other network models
   if(network.model!="Empirical"){
      if(network.model=="Erdos"){model<-"erdos.renyi.game"}
      if(network.model=="Watts"){model<-"watts.strogatz.game"}
      if(network.model=="NW.Watts"){model<-"NW.watts.strogatz.game"}
      if(network.model=="Barabasi"){model<-"barabasi.game"}
      #PTN FILE HERE


   }   

   #Reparametrize model if user-supplied behavioral model is valid
   if (!is.null(usbm)){
     behavioral.weights <- usbm$behavioral.weights
     Stiffness <- usbm$Stiffness
     Personal.Delta <- usbm$Personal.Delta
     behavioral.model <- "USBM"

     for(i in 1:length(usbm$dfuns)){
        assign(names(usbm$dfuns)[i],usbm$dfuns[[i]],envir=nirm.env)
     }

     model.specs <- c(model.specs, Behavioral.Model = Behavioral.Model)
     params.of.model <- c(params.of.model, Personal.Delta = Personal.Delta, Stiffness = Stiffness)
    }

   
   #Set the behavioral model
   set.behavioral.model(behavioral.model,params.of.model,behavioral.weights,Personal.Delta,Stiffness,names.outcomes,verbose,model.specs,nirm.env)

   #What time is it now?
   rtime1<-proc.time()[3]
   #DYNAMICS FILE HERE


   #Finalizing the output for 
   #what time is it now?
   rtime2<-proc.time()[3]
   params.of.model<-c(params.of.model,network.model,behavioral.model)
   if(length(params.of.model)==20){
   names(params.of.model)<-c("N","ave.degree","s.mem","initial.inf.prop","R0.pref","R0.pop","gamma","total.flu.seasons","eff","added.RMCN",
   "Personal.Delta.VaccInf","Personal.Delta.Vacc","Personal.Delta.Inf","Personal.Delta.UnInf","Stiffness.VaccInf","Stiffness.Vacc","Stiffness.Inf",
   "Stiffness.UnInf","model","behavioral.model")
   }
   if(length(params.of.model)==21){
   names(params.of.model)<-c("N","ave.degree","s.mem","initial.inf.prop","R0.pref","R0.pop","gamma","total.flu.seasons","eff","added.RMCN",
   "ave.degree2","Personal.Delta.VaccInf","Personal.Delta.Vacc","Personal.Delta.Inf","Personal.Delta.UnInf","Stiffness.VaccInf",
   "Stiffness.Vacc","Stiffness.Inf","Stiffness.UnInf","model","behavioral.model")
   params.of.model<-params.of.model[-11]
   }
   if(length(params.of.model)==22){
   names(params.of.model)<-c("N","ave.degree","s.mem","initial.inf.prop","R0.pref","R0.pop","gamma","total.flu.seasons","eff","added.RMCN",
   "ave.degree2","phi","Personal.Delta.VaccInf","Personal.Delta.Vacc","Personal.Delta.Inf","Personal.Delta.UnInf","Stiffness.VaccInf",
   "Stiffness.Vacc","Stiffness.Inf","Stiffness.UnInf","model","behavioral.model")
   params.of.model<-params.of.model[-11]
   }
   if(length(params.of.model)==25){
   names(params.of.model)<-c("N","ave.degree","s.mem","initial.inf.prop","R0.pref","R0.pop","gamma","total.flu.seasons","eff","added.RMCN",
   "model","ave.degree2","power","ave.degree3","phi","Personal.Delta.VaccInf","Personal.Delta.Vacc","Personal.Delta.Inf",
   "Personal.Delta.UnInf","Stiffness.VaccInf","Stiffness.Vacc","Stiffness.Inf","Stiffness.UnInf","model","behavioral.model")
   params.of.model<-params.of.model[-c(11,12,14)]
   }
   outpm<-list(Params=params.of.model,Results=seasonal.output,Runtime=rtime2-rtime1,File=output.file.name)
   class(outpm)<-"nirm"
   return(outpm)
}

print.nirm<-function(x,...){
   cat("Summary of Network Inductive Reasoning Model \n")
   cat("\t Model Parameters \n")
   cat(" Size of Network: \t",x$Params["N"],"\n")
   cat(" Network Model: \t",x$Params["model"],"\n")
   cat(" Average Degree: \t",x$Params["ave.degree"],"\n")
   cat(" Behavioral Model: \t",x$Params["behavioral.model"],"\n\n")
   cat(" Number of Seasons: \t",x$Params["total.flu.seasons"],"\n")
   cat(" Initial Pr(Cond): \t",x$Params["initial.inf.prop"],"\n")
   cat(" Suppressant Efficacy: \t",x$Params["eff"],"\n")
   cat(" SIR Parameters: \n\t ",
      paste("R_0= ",round(as.numeric(x$Params["R0.pref"]),3),";","R_0.pop= ",round(as.numeric(x$Params["R0.pop"]),3),";","gamma= ",round(as.numeric(x$Params["gamma"]),3)),"\n\n")
   
   #Define Harmonic Mean for Detecting Epidemics as H+1.96*(sd(H))
   #H.mean<-function(m){length(m)/sum(1/m)}
   #H.sd<-function(m){sqrt(sum((m-H.mean(m))^2)/length(m))}
   N<-as.numeric(x$Params["N"])
   w.per.season<-unlist(lapply(x$Results,function(x) mean(x[[1]][,"w"])))
   inc.per.season<-unlist(lapply(x$Results,function(x) sum(x[[1]][,"inf"])/N))
   cat("\t Seasonal Distributions \n") 
   cat(" E[W]\n")
   print(round(quantile(w.per.season,prob=c(0.025,.5,.975)),3))
   cat(" Incidence\n")
   print(round(quantile(inc.per.season,prob=c(0.025,.5,.975)),3))
   cat("Seasonal Correlation between E[W] and Incidence: ",round(cor(w.per.season,inc.per.season),3),"\n\n")
   cat("This model took ",x$Runtime, "seconds to run.\n")
   invisible(x)
}

read.nirm.log<-function(file){
   if(!file.exists(file)){stop(paste("File, ",file," does not exist.",sep=""))}
   return(read.csv(file,h=TRUE, stringsAsFactors=FALSE))
}

write.nirm.log<-function(x,loc="",append=TRUE){
   if(class(x)!="nirm"){stop("Not a nirm object.")}
   file<-x$File
   if(!file.exists(loc)){warning(paste("Log file does not exist. Writing new file.",sep="")); append<-FALSE}
   cnams<-c("filename","Network","Behavior","N","ave.degree","s.mem","initial.inf.prop","R0.pref","R0.pop","gamma","total.flu.seasons","eff","added.RMCN","Date")
   sys<-Sys.info()
   of1<-rbind(c(file,x$Params["model"],x$Params["behavioral.model"],x$Params["N"],x$Params["ave.degree"],x$Params["s.mem"],x$Params["initial.inf.prop"],x$Params["R0.pop"],x$Params["R0.pop"],x$Params["gamma"],x$Params["total.flu.seasons"],x$Params["eff"],x$Params["added.RMCN"],date()))
   colnames(of1)<-cnams
   write.table(file=loc,of1,sep=",",append=append,col.names=!append,row.names=FALSE)
}

load.nirm<-function(x,output.dir="Temp",what="all"){
   if(class(x)!="nirm"){stop("Not a nirm object.")}
   rex<-strsplit(x$File,"_")[[1]][length(strsplit(x$File,"_")[[1]])]
   if(what=="all"){
      flist<-dir(paste(output.dir,"Output","Dynamics","",sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/")))
      load(paste(output.dir,"Output","Dynamics",flist[grep(rex,flist)],sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/")),.GlobalEnv)
      flist<-dir(paste(output.dir,"Output","Networks","",sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/")))
      load(paste(output.dir,"Output","Networks",flist[grep(rex,flist)],sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/")),.GlobalEnv)
      return(print(paste("Loaded all files associated with model",rex,sep=" ")))
   }
   if(what=="network"){
     flist<-dir(paste(output.dir,"Output","Networks","",sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/")))
     load(paste(output.dir,"Output","Networks",flist[grep(rex,flist)],sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/")),.GlobalEnv)
      return(print(paste("Loaded network files associated with model",rex,sep=" ")))
   }
   if(what=="dynamics"){
     flist<-dir(paste(output.dir,"Output","Dynamics","",sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/")))
     load(paste(output.dir,"Output","Dynamics",flist[grep(rex,flist)],sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/")),.GlobalEnv)
      return(print(paste("Loaded dynamics files associated with model",rex,sep=" ")))
   }
}

analyze.nirm<-function(m1,output.dir="Temp",verbose=TRUE){
load.nirm(m1,output.dir=output.dir)

#Collect what we need from nirm object
rex<-strsplit(m1$File,"_")[[1]][length(strsplit(m1$File,"_")[[1]])]
output.analysis.dir<-paste(output.dir,"Analysis","",sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/"))
output.figure.dir<-paste(output.dir,"Figure","",sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/"))
output.file.name<-paste("Analysis",rex,sep=".")

for(y in 1:length(names(m1$Params[which(!names(m1$Params)%in%c("model","behavioral.model"))]))){
assign(names(m1$Params[which(!names(m1$Params)%in%c("model","behavioral.model"))])[y],as.numeric(m1$Params[which(!names(m1$Params)%in%c("model","behavioral.model"))])[y],.GlobalEnv)
}

for(y in 1:length(names(m1$Params[which(names(m1$Params)%in%c("model","behavioral.model"))]))){
assign(names(m1$Params[which(names(m1$Params)%in%c("model","behavioral.model"))])[y],as.character(m1$Params[which(names(m1$Params)%in%c("model","behavioral.model"))])[y],.GlobalEnv)
}

#Clean up with garbage collection 
# (we don't need this anymore and we could use the memory)
rm(m1)
gc()

####Analysis
### Cluster, Degree and knn analysis.
g.cluster <- clusters(g.PTN)
g.k       <- degree(g.PTN)
dist.g    <- degree.distribution(g.PTN, cumulative = FALSE)
g.knn     <- graph.knn(g.PTN)$knn
g.knnk    <- graph.knn(g.PTN)$knnk

### Split the degree into 5 groups from low (1) to high (5) degree.
split.g.k <- split.degree.groups(g.k,5)
quan.g.k <- split.g.k$quan.g.k
g.k.cat  <- split.g.k$g.k.cat
g.k.label <- split.g.k$g.k.label

### stats on k dist.
ave.d <- mean(g.k)
ave.d2 <- mean((g.k)^2)
var.d <- var(g.k)

### Sanity check to check that R0 and T are what is expected
### for the case of an uncorrelated network
T <- R0.pref /ave.d
Theta.PTN <- T*(ave.d2-ave.d)/ave.d  
pc.PTN <- 1-1/Theta.PTN

ave.knn <- mean(g.knn)
var.knn <- var(g.knn)

PTN.stats <- list(g.cluster=g.cluster,g.k=g.k,dist.g=dist.g,g.knn=g.knn,g.knnk=g.knnk,
                  split.g.k=split.g.k,quan.g.k=quan.g.k,g.k.cat=g.k.cat,ave.d=ave.d,
                  ave.d2=ave.d2,var.d=var.d,T=T,Theta.PTN=Theta.PTN,pc.PTN=pc.PTN,
                  ave.knn=ave.knn,var.knn=var.knn)

#############################################################
###                                                       ###
###             COMPUTE AVERAGE DYNAMICS                  ###
###                                                       ###
#############################################################

if(verbose) print(":: ANALYZE: COMPUTE AVERAGE DYNAMICS ::")

coverage <- NULL
cuml.inc <- NULL
sus.uninf <- NULL

Dyn <- data.frame(NULL)
w.summary <- data.frame(NULL)

for(flu.season in 1:total.flu.seasons){
  Dyn <- rbind(Dyn,c(
    flu.season,  
    sum(seasonal.output[[flu.season]]$states[,"sus"])/N,
    sum(seasonal.output[[flu.season]]$states[,"vacc"])/N,
    sum(seasonal.output[[flu.season]]$states[,"vaccinf"])/N,
    sum(seasonal.output[[flu.season]]$states[,"inf"])/N ))
    w.summary <- rbind(w.summary,c(flu.season,  
                  mean(seasonal.output[[flu.season]]$states[,"w"]),
                  sqrt(var(seasonal.output[[flu.season]]$states[,"w"]))))
}
colnames(w.summary)<- c("year","w","w.sd")
colnames(Dyn) <- c("year","s","p","vi","nvi")
Dyn <- cbind(Dyn, i = Dyn$vi+Dyn$nvi)


#############################################################
###                                                       ###
###               QUICK CLUSTER ANALYSIS                  ###
###                                                       ###
#############################################################

if(verbose) print(":: ANALYZE: QUICK CLUSTER ANALYSIS ::")

quick.cluster.analysis <- Quick.cluster.analysis(seasonal.output,total.flu.seasons)

#############################################################
###                                                       ###
###                   PERIOD ANALYSIS                     ###
###                                                       ###
#############################################################

if(verbose) print(":: ANALYZE: PERIOD ANALYSIS ::")

remove.transiet <- 20
peak.threshold <- 0.7
periods.to.peaks <- get.times.between.spikes(Dyn$i,th=peak.threshold,remove.transiet= remove.transiet)

peak.years <- 2+cumsum(periods.to.peaks)+remove.transiet


mean.periods.to.peaks <- mean(periods.to.peaks)
var.periods.to.peaks  <- var(periods.to.peaks)

Period.stats <- list(remove.transiet=remove.transiet,periods.to.peaks=periods.to.peaks,
                     mean.periods.to.peaks=mean.periods.to.peaks,var.periods.to.peaks=var.periods.to.peaks,
                     peak.threshold=peak.threshold,peak.years=peak.years)

#############################################################
###                                                       ###
###                  SWITCHING RATES                      ###
###                                                       ###
#############################################################

if(verbose) print(":: ANALYZE: SWITCHING RATES ::")

Switch.Rate <- Switch.Behavior.Rate.Analysis(seasonal.output,g.k,remove.transiet=1)

#############################################################
###                                                       ###
###              ANALYSE THE w DISTRIBUTION               ###
###                                                       ###
#############################################################

if(verbose) print(":: ANALYZE: w DISTRIBUTION ::")

w.distrubution.tables <- w.distrubution.analysis(seasonal.output,
                        total.flu.seasons,g.k.cat,g.k.label,num.breaks=50,plot=FALSE)


#############################################################
###                                                       ###
###         MEASURE CORRELATIONS IN SPACE AND TIME         ###
###                                                       ###
#############################################################

if(verbose) print(":: ANALYZE: MEASURE CORRELATIONS IN SPACE AND TIME ::")

correlations <- Correlation.Analyses(g.PTN,g.k,g.knn,seasonal.output,flu.season,total.flu.seasons=total.flu.seasons)
nn.like.with.like <- nn.like.with.like.analysis(g.PTN,g.k,Dyn,seasonal.output,total.flu.seasons=total.flu.seasons)

#############################################################
###                                                       ###
###                   SAVE THE ANALYSIS                   ###
###                                                       ###
#############################################################

if(verbose) print(":: ANALYZE: SAVING ::")

save(params.of.model,PTN.stats,percolation.analisis,Dyn,w.summary,
     quick.cluster.analysis,
     Period.stats,Switch.Rate,w.distrubution.tables,correlations,
     nn.like.with.like, 
     file = paste(output.analysis.dir,output.file.name,sep=""),compress="xz")

}


percolation.nirm<-function(m1,new=FALSE,perc.file=NULL,output.dir="Temp",mc.perc.runs=100,verbose=TRUE,seed=55279,plot=TRUE,...){

if(!new){
   if(is.null(perc.file) | file.exists(perc.file)==FALSE){stop("Percolation R-Object File Not Found.")}
   load(perc.file)
}

#############################################################
###                                                       ###
###              Percoloation Analysis                    ###
###                                                       ###
#############################################################
load.nirm(m1,output.dir=output.dir)

#Collect what we need from nirm object
rex<-strsplit(m1$File,"_")[[1]][length(strsplit(m1$File,"_")[[1]])]
output.analysis.dir<-paste(output.dir,"Analysis","",sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/"))
output.figure.dir<-paste(output.dir,"Figure","",sep=ifelse(tolower(Sys.info()["sysname"])=="windows","\\","/"))
output.file.name<-paste("Percolation",rex,sep=".")

for(y in 1:length(names(m1$Params[which(!names(m1$Params)%in%c("model","behavioral.model"))]))){
assign(names(m1$Params[which(!names(m1$Params)%in%c("model","behavioral.model"))])[y],as.numeric(m1$Params[which(!names(m1$Params)%in%c("model","behavioral.model"))])[y],.GlobalEnv)
}

for(y in 1:length(names(m1$Params[which(names(m1$Params)%in%c("model","behavioral.model"))]))){
assign(names(m1$Params[which(names(m1$Params)%in%c("model","behavioral.model"))])[y],as.character(m1$Params[which(names(m1$Params)%in%c("model","behavioral.model"))])[y],.GlobalEnv)
}

#Clean up with garbage collection 
# (we don't need this anymore and we could use the memory)
rm(m1)
gc()

####Analysis
### Cluster, Degree and knn analysis.
g.cluster <- clusters(g.PTN)
g.k       <- degree(g.PTN)
dist.g    <- degree.distribution(g.PTN, cumulative = FALSE)
g.knn     <- graph.knn(g.PTN)$knn
g.knnk    <- graph.knn(g.PTN)$knnk

### Split the degree into 5 groups from low (1) to high (5) degree.
split.g.k <- split.degree.groups(g.k,5)
quan.g.k <- split.g.k$quan.g.k
g.k.cat  <- split.g.k$g.k.cat
g.k.label <- split.g.k$g.k.label

### stats on k dist.
ave.d <- mean(g.k)
ave.d2 <- mean((g.k)^2)
var.d <- var(g.k)

### Sanity check to check that R0 and T are what is expected
### for the case of an uncorrelated network
T <- R0.pref /ave.d
Theta.PTN <- T*(ave.d2-ave.d)/ave.d  
pc.PTN <- 1-1/Theta.PTN

ave.knn <- mean(g.knn)
var.knn <- var(g.knn)

if(new){
  models <- c("erdos.renyi.game","NW.watts.strogatz.game","barabasi.game")
  model<-models[grep(tolower(model),tolower(models))]

p.seq <- c(0:mc.perc.runs)/mc.perc.runs
#p.seq <- c(0:10)/10

### list of nearest neighbours (nn) on the static Social Network. 
nn.list.PTN <- get.nn.list.of.SN(g.PTN)
N.nn <- degree(g.PTN)


### Do Theoretical Percolation Analysis 
Theta.PTN <- R0.pref*(mean(N.nn^2)-mean(N.nn))/mean(N.nn)^2
pc.PTN <- 1-1/Theta.PTN


### Do Theoretical Percolation Analysis 
Theta.PTN <- R0.pref*(mean(N.nn^2)-mean(N.nn))/mean(N.nn)^2
pc.PTN <- 1-1/Theta.PTN
percolation.analisis <- list(Theta.PTN=Theta.PTN,pc.PTN=pc.PTN)

set.seed(seed)
  q.table.random <- SIR.Prob.Infection.Table(model=model,N=N,ave.degree=ave.degree,R0.pref,R0.pop,eff=eff,p.seq=p.seq, 
                                             mc.perc.runs=mc.perc.runs,
                                             added.RMCN=added.RMCN,preferential=FALSE,verbose=verbose)
  
  q.table.preferential <- SIR.Prob.Infection.Table(model=model,N=N,ave.degree=ave.degree,R0.pref,R0.pop,eff=eff,p.seq=p.seq, 
                                                   mc.perc.runs=mc.perc.runs,
                                                   added.RMCN=added.RMCN,preferential=TRUE,verbose=verbose)
  
  
  percolation.analisis <- list(Theta.PTN=Theta.PTN,pc.PTN=pc.PTN, q.table.random=q.table.random,
                               q.table.preferential=q.table.preferential,
                               R0.pref=R0.pref,R0.pop=R0.pop,eff=eff,
                               p.seq= p.seq, mc.perc.runs=mc.perc.runs,added.RMCN=added.RMCN)
}
save(N,phi,N.nn,nn.list.PTN,percolation.analisis, file = paste(output.analysis.dir,output.file.name,sep=""),compress="xz")

if(plot==TRUE){
pdf(paste(output.figure.dir,output.file.name,".pdf",sep=""),pointsize=20,width=10,height=10)
par(lwd=2)
xax<-unlist(lapply(seasonal.output,function(x) sum(x[[1]][,"vacc"])/N))
Inc.net<-unlist(lapply(seasonal.output,function(x) sum(x[[1]][,"inf"])/N))

plot(percolation.analisis$q.table.random$q.table.Inc$'50%'~percolation.analisis$q.table.random$q.table.Inc$p,type="l",lty=3,ylab="Proportion Infected",xlab="Proportion Vaccinated",cex=3,...)
segments(y0=percolation.analisis$q.table.random$q.table.Inc$'15.8%',y1=percolation.analisis$q.table.random$q.table.Inc$'84.2%',x0=percolation.analisis$q.table.random$q.table.Inc$p,x1=percolation.analisis$q.table.random$q.table.Inc$p)

points(y=Inc.net,x=xax,pch=19)
try(points(supsmu(y=Inc.net[-which(Inc.net<.005)],x=xax[-which(Inc.net<.005)]),type="l"),silent=TRUE)

legend("topright",c("Random Model","Network Model"),lty=c(3,1),bty="n")

dev.off()

   }
}


set.behavioral.model<-function(Behavioral.Model,params.of.model,behavioral.weights,Personal.Delta,Stiffness,names.outcomes,verbose,model.specs,nirm.env){
   if(Behavioral.Model=="USBM" || is.null(Behavioral.Model)){return(print("Using User-Supplied Behavioral Model"))}
   beh.env<-environment()
   ### The Global.Vacc.Delta.Threshold calls pc.PTN. We need to construct this 
   ### here because it is not a priori in the code. --- CSM & ARJ --- 18 July 2012
   g.PTN<-get("g.PTN",envir=nirm.env)
   R0.pref<-get("R0.pref",envir=nirm.env)
   N.nn <- degree(g.PTN)
   Theta.PTN <- R0.pref * (mean(N.nn^2) - mean(N.nn))/mean(N.nn)^2
   pc.PTN <- 1 - 1/Theta.PTN

   #BEH MODELS HERE


   #Send the model back to the nirm environment
   for(i in 1:length(ls(envir=beh.env))){
      assign(ls(envir=beh.env)[i],get(ls(envir=beh.env)[i]),envir=nirm.env)
   }
}

   #FUNCTIONS HERE

#New percolation plotting function
plot.qtable<-function(perc,nob,conf.int=TRUE,show.points=TRUE,pcol="black",smooth=TRUE,add=FALSE,...){
   #perc is a three column matrix with columns 'p', 'q(p)' and 'var[q(p)]', respectively.
   #this function uses the bootstrapped variance to calculate the approximate 95% confidence intervals
   x<-perc[,1]
   y<-perc[,2]
   upbar<-y+1.96*sqrt(perc[,3])
   dnbar<-y-1.96*sqrt(perc[,3])
   thedata<-do.call("rbind",lapply(nob$Results,function(x) c(sum(x$states[,"vacc"]),sum(x$states[,"inf"]))/as.numeric(nob$Params["N"])))
   ywin<-c(0,max(c(y,thedata[,2],upbar)))
   xwin<-c(0,max(c(x,thedata[,1])))
   if(!add){
   plot(x=x,y=y,ylab="q(p)",xlab="p",ylim=ywin,xlim=xwin,...)
   if(conf.int){segments(x0=x,x1=x,y0=upbar,y1=dnbar)}
   if(show.points){points(thedata,pch=19,col=pcol)}
   if(smooth){points(supsmu(thedata[,1],thedata[,2]),type="l")}
   }
   if(add){
   points(x=x,y=y,...)
   if(conf.int){segments(x0=x,x1=x,y0=upbar,y1=dnbar)}
   if(show.points){points(thedata,pch=19,col=pcol)}
   if(smooth){points(supsmu(thedata[,1],thedata[,2]),type="l")}
   }
}

#EOF
