#### A 3D semi-spatially explicit Neutral model

#Put the working directory here and clear memory
rm(list=ls())

#Set working directory, make sure the directory holds all necessary data files as specified in the github repository
setwd("set workingdirectory")

source("Functions_3Dmodel.R")			#Functions for this script

library(untb)						
library(vegan)						
library(labdsv)						
library(quantreg)					
library(data.table)
library(plyr)
library(arrayhelpers)
library(rgl)
library(rms)
library(gplots)
library(plotrix)
library(dunn.test)
library(gridGraphics)
library(grid)
library(gplots)
library(gridExtra)
library(abind)
library(foreach)
library(doParallel)
library(doMC)
library(doSNOW)
library(ggplot2)
library(ggmap)
library(sp)
library(rgdal)
library(rgeos)
library(reshape2)
library(Matching)
library(boot)

##############################################################################
############################# END STARTING SEQUENCE ##########################
##############################################################################

#########################
### Running the model ###
#########################

#Read in species data, plot data and hypothetical metacommunity data based on Guyana/Suriname
vegtab = read.csv("species_data.csv", row.name=1)
plotdata = read.csv("plot_data.csv")
datatometacom = read.table("metaguianas.txt",sep="\t")

##Parameters for metacommunity based on Guiana's
Jm = 20191600511
Sm = 4582
RAD.RD = as.numeric(as.matrix(datatometacom))
RD.spp = 1:(Sm-1)
length(RAD.RD) == length(RD.spp) #Check, should be true
fa.rd = fishers.alpha(Jm,Sm)


## Initialize parameters
step			= 0
Nrep			= 10
Np			= 15
Ind			= ceiling(mean(rowSums(vegtab))) 
Theta		= fa.rd
pool			= 1:(Np*Np*Ind)
	
m.adj		=	.141 #For Guyana/Suriname: see appendix SA Table S2
m.forest		=	.1*m.adj
m.meta		=	.01*m.adj
m.plot		=	1-(m.adj+m.forest+m.meta)
mu			= 	Theta/(2*(Np*Np*Ind))

## Create 3D forest (field) to start with
Createforest.3d(Np,Ind, mixed=T)
initial.nspecies = unique(levels(factor(Forest)))
init.nr.species = length(initial.nspecies); init.nr.species
initnr.species.field = length(colnames(vegtab)); initnr.species.field

##################################################################################################
##################################################################################################
##################################################################################################
##Function to run 3D process including speciation and dispersal categories

for (x in 1:10){
cores = detectCores()-1
cl <- makeCluster(cores)
registerDoParallel(cl)
on.exit(stopCluster(cl))
forests = NULL

#loop
ptm <- proc.time()
forests = foreach(i = 1:cores) %dopar% {

for(i in 1:Nrep){
Forestcopy = Forest

####The actual sampling proces
##First corner (upperleft)
x=1;y=1
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m == 1){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x,y,],1)
} else if(m == 2){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x:(x+1),y:(y+1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD) 
} else { 	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}

##Second corner (upperright)
x=Np; y=1
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m == 1){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x,y,],1)
} else if(m == 2){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x:(x-1),y:(y+1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD) 
} else { 	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}

##Third corner (lowerleft)
x=1;y=Np
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m == 1){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x,y,],1)
} else if(m == 2){ Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x:(x+1),y:(y-1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD) 
} else { 	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}

##Fourth corner (lowerright)
x=Np;y=Np
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m == 1){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x,y,],1)
} else if(m == 2){ Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x:(x-1),y:(y-1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD) 
} else {  	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}

##First edge (top)
for (x in 2:(Np-1)){
y=1
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m == 1){Forest[x,y,sample(1:Ind,1)] = sample(Forest[x,y,],1)
} else if(m == 2){ Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[(x-1):(x+1),y:(y+1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD) 
} else {  	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}}

##Second edge(left)
for(y in 2:(Np-1)){
y=1
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m== 1){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x,y,],1)
} else if(m == 2){ Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x:(x+1),(y-1):(y+1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD)
} else {  	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}}

##Third edge (right)
for(y in 2:(Np-1)){
x=Np
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m == 1){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x,y,],1)
} else if(m == 2){ Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x:(x-1),(y-1):(y+1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD)
} else { 	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}}

##Fourth edge (bottom)
for (x in 2:(Np-1)){
y=Np
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m == 1){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x,y,],1)
} else if(m == 2){ Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[(x-1):(x+1),y:(y-1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD)
} else { 	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}}

##Rest of the forest
for (x in 2:(Np-1)){
	for(y in 2:(Np-1)){
m = sample(1:5,1, prob = c(m.plot, m.adj, m.forest, m.meta, mu), replace = T)
if(m == 1){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[x,y,],1)
} else if(m == 2){ Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy[(x-1):(x+1),(y-1):(y+1),],1)
} else if(m == 3){Forest[x,y,sample(1:Ind,1)] = sample(Forestcopy,1)
} else if(m == 4){Forest[x,y,sample(1:Ind,1)] = sample(RD.spp,	1, prob = RAD.RD)
} else {  	present = as.numeric(levels(factor(unique(Forestcopy))))
			future = pool[-present]
			Forest[x,y,sample(1:Ind,1)] = sample(future,1, replace = F)}}}
} 
forests[i] = Forest}
cat(timing = proc.time() - ptm)
stopCluster(cl)

Forest_inside = array(unlist(forests[1]), dim = c(Np,Np,Ind))

for(i in 2:length(forests)){
forest = array(unlist(forests[i]), dim = c(Np,Np,Ind))
Forest_inside = abind(Forest_inside, forest, along=2)
}
  if(x==1){Forest_complete = Forest_inside}
  else{Forest_complete = abind(Forest_complete, Forest_inside, along=2)}
}

saveRDS(Forest_complete, file = "example_community.rds")

