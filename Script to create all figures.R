## Run the "parallel" model from the R script "Create simulation_parallel" and then read datasets from the parallel model output, here the example community is used based on the supplied Guyana/Suriname dataset. For Figure 4 of the manuscript the additional supplied subnamed "uni" and "null" should be used.

Forest_complete =  readRDS("example_community.rds")

##############################################
#### FIGURE 1: RAD AND MDD EMPIRICAL DATA ####
##############################################

pdf("Figure 1.pdf", paper="special", height=15, width=11)
par(mfrow=c(3,2))

iterations = 10
vegtablist = list()
vegtab.field = vegtab

#Create plotdata, respecify nr ind when reading in parallel forests, otherwise it does not work
Ind	= ceiling(mean(rowSums(vegtab.field))) 
Plotdata.3D(Forest_complete)

for(j in 1:iterations){
samplesize = length(rownames(vegtab.field))
invisible(capture.output(create.data3D(Forest_complete, sample = T, n=samplesize)))
vegtablist[j] = list(vegtab.sim3D)
}


### Create RAD plots ###
n.ind.spec.sim = lapply(vegtablist,FUN=function(x)colSums(x))
ME.RAD.sims	=	lapply(n.ind.spec.sim, FUN=function(x) rev(sort(x)))
species.plot = melt(lapply(vegtablist,FUN=function(x)(mean(sum(colSums(x)>0)))))
species.sim	=	mean(melt(lapply(vegtablist,FUN=function(x)(mean(sum(colSums(x)>0)))))$value)

empty2=matrix(, nrow=iterations, ncol=max(species.plot))

	for(q in 1:length(ME.RAD.sims)){
		empty2[q,1:length(ME.RAD.sims[[q]])] = as.numeric(ME.RAD.sims[[q]])
				}

abundances = colMeans(empty2[,1:species.sim], na.rm=T)
ranks = seq(1:length(abundances))

min.plot = species.plot[species.plot$value == min(species.plot$value),]$L1
min.abundance = as.numeric(na.omit(empty2[min.plot,]))
ranks.min = seq(1:length(min.abundance))

max.plot = species.plot[species.plot$value == max(species.plot$value),]$L1
max.abundance = as.numeric(na.omit(empty2[max.plot,]))
ranks.max = seq(1:length(max.abundance))

#LOESS prediction FIT
myloess <- loess(abundances~ranks, na.action = na.exclude, span=.1) 
prediction = predict(myloess, ranks, se=T)

#LOESS prediction MIN
myloess_min = loess(min.abundance~ranks.min, na.action = na.exclude, span=.10) 
prediction_min = predict(myloess_min, ranks.min, se=T)

#LOESS prediction MAX
myloess_max = loess(max.abundance~ranks.max, na.action = na.exclude, span=.10) 
prediction_max = predict(myloess_max, ranks.max, se=T)

##Rank abundance curve Fielddata
n.ind.spec.field	=	colSums(vegtab.field)
n_ind.plots.field	= rowSums(vegtab.field) 
ME.RAD.field		=	rev(sort(n.ind.spec.field))

## Create fit for logseries and plot lines of logseries based on loess fit
Jm = round(sum(abundances))
Sm = length(ranks)
fisher.fit =  fisher.fit = fisher.ecosystem(Jm,Sm, nmax=Sm); fisher.fit = na.omit(fisher.fit)
plot.fisher.fit = as.numeric(as.matrix(fisher.fit))
forest.alpha = fishers.alpha(Jm,Sm); forest.alpha; sum(round(unique(prediction$fit)) == 1)

#Check comparisons
kstest = ks.boot(prediction$fit,as.numeric(ME.RAD.field))
kstest_inter = ks.boot(prediction$fit[prediction$fit > 125 & prediction$fit <725], as.numeric(ME.RAD.field)[as.numeric(ME.RAD.field)>125 & as.numeric(ME.RAD.field)<725])

#Plot data 
plot(abundances,ranks, type="n", ylim = c(1,2500), xlim = c(0,2000), log="y", xlab="Rank", ylab="Abundance", cex.axis=1.5, cex.lab=1.5, main = "Regional prediction", cex.main=1.5)
legend("topright",
c("Empirical RAD","Logseries RAD LOESS prediction","RAD LOESS predicton (50 draws)","Range RAD LOESS predictions"),
text.col=c("darkgreen", "red","black","deepskyblue2"),
cex=1.6, bty="n")
legend("bottomleft", bty="n",c(paste("D =", round(as.numeric(kstest$ks[1]),3)),paste("P <",.05-kstest$ks.boot.pvalue),paste("P_inter =",kstest_inter$ks.boot.pvalue)), cex=1.6)
	
#Plot lines for fisher and field
lines(plot.fisher.fit, col = "red", lwd=2)  
lines(ME.RAD.field, col = "darkgreen", lwd=2)

#Plot the LOESS regression with banner based on min and max
lines(ranks.min, prediction_min$fit, lty="dashed", col="lightskyblue1", lwd=1) 
lines(ranks.max, prediction_max$fit, lty="dashed", col="lightskyblue1", lwd=1) 
polygon(c(ranks.min, rev(ranks.max)), c(prediction_min$fit, rev(prediction_max$fit)), col="deepskyblue2", border=NA)
lines(ranks, prediction$fit, lty="solid", col="black", lwd=2) 


########################	
### Create MDD plots ###
########################

max.sim.pervegtab 	= lapply(vegtablist,FUN=function(x)apply(x,1,function(x)return(array(max(x)))))
max.sims				= lapply(max.sim.pervegtab, FUN=function(x) rev(sort(x)))
max.dom.sims			= data.frame(lapply(max.sim.pervegtab, FUN=function(x) max(x)))

empty3				= matrix(, nrow=iterations, ncol=length(max.sims[[1]]))

	for(d in 1:length(max.sims)){
		empty3[d,1:length(max.sims[[d]])] = as.numeric(max.sims[[d]])
				}

abundances = round(colMeans(empty3,na.rm=T))
ranks = seq(1:length(abundances))

min.plot = which(max.dom.sims == min(max.dom.sims))
min.abundance = as.numeric(na.omit(empty3[min.plot,]))
ranks.min = seq(1:length(min.abundance))

max.plot = which(max.dom.sims == max(max.dom.sims))
max.abundance = as.numeric(na.omit(empty3[max.plot,]))
ranks.max = seq(1:length(max.abundance))

#LOESS prediction FIT
myloess <- loess(abundances~ranks, na.action = na.exclude, span=.1) 
prediction = predict(myloess, ranks, se=T)

#LOESS prediction MIN
myloess_min = loess(min.abundance~ranks.min, na.action = na.exclude, span=.10) 
prediction_min = predict(myloess_min, ranks.min, se=T)

#LOESS prediction MAX
myloess_max = loess(max.abundance~ranks.max, na.action = na.exclude, span=.10) 
prediction_max = predict(myloess_max, ranks.max, se=T)

#Field MDD
max.field<-apply(vegtab.field,1,function(x)return(array(max(x))))
max.RAD.field	=	rev(sort(max.field))

kstest = ks.boot(prediction$fit,as.numeric(max.RAD.field))

# Plot data
plot(max.RAD.field, col = "darkgreen", type="l", lwd=2, xlab="Rank", ylab="Abundance most dominant species", cex.axis=1.5, cex.lab=1.5, log="y", main= "Local prediction", cex.main=1.5, xlim=c(0,100), ylim=c(1,300))
legend("topright",
c("Empirical MDD","MDD LOESS predicton (50 draws)","Range RAD LOESS predictions"),
text.col=c("darkgreen","black","deepskyblue2"),
cex=1.6, bty="n")
legend("bottomleft", bty="n",c(paste("D =", round(as.numeric(kstest$ks[1]),3)),paste("P <",.05-kstest$ks.boot.pvalue)), cex=1.6)

#Plot the LOESS regression with banner based on min and max
lines(ranks.min, prediction_min$fit, lty="dashed", col="lightskyblue1", lwd=1) 
lines(ranks.max, prediction_max$fit, lty="dashed", col="lightskyblue1", lwd=1) 
polygon(c(ranks.min, rev(ranks.max)), c(prediction_min$fit, rev(prediction_max$fit)), col="deepskyblue2", border=NA)
lines(ranks, prediction$fit, lty="solid", col="black", lwd=2) 
dev.off()
######################
#### END FIGURE 1 ####
######################



############################
#### FIGURE 2: BOXPLOTS ####
############################

wilcoxtable = matrix(nrow=3,ncol=3)
rownames(wilcoxtable) = c("mean nr species","mean nr singletons","mean FA per plot")
colnames(wilcoxtable) = c("Pvalue OSB","Pvalue DS","Pvalue JEG")
output = list()

#Some details of the model
n_spec.plots.sim = melt(lapply(vegtablist,FUN=function(x)rowSums(x>0)))
n_ind.plots.sim = mean(melt(lapply(vegtablist,FUN=function(x)rowSums(x)))$value)
n_single.sim = melt(lapply(vegtablist,FUN=function(x)(rowSums(x==1))))$value
fa.plots.sim = melt(lapply(vegtablist,FUN=function(x) mapply(fishers.alpha,rowSums(x),rowSums(x>0))))
species.sim	=	mean(melt(lapply(vegtablist,FUN=function(x)(mean(sum(colSums(x)>0)))))$value)
total.fa =  mean(melt(lapply(vegtablist,FUN=function(x)(mean(fishers.alpha(sum(x),sum(colSums(x)>0))))))$value)
max.dom.sims		= melt(lapply(vegtablist,FUN=function(x)apply(x,1,function(x)return(array(max(x))))))

n_spec.plots			=	rowSums(vegtab.field>0)
n_ind.plots			=	rowSums(vegtab.field)
n_single				=	sum(colSums(vegtab.field) == 1)
n_single.plots		=	rowSums(vegtab.field == 1)
mean_single_field   =	mean(rowSums(vegtab.field == 1))
fa.plots				=	mapply(fishers.alpha,n_ind.plots,n_spec.plots)
fa.plots_field 		=	fa.plots
species				=	length(colnames(vegtab.field))
total.fa				=   fishers.alpha(sum(vegtab.field),species)
max.dominance.plots	=	apply(vegtab.field,1,max)

output[["Guyana/Suriname"]] = list(n_spec.plots.sim$value, n_spec.plots, n_single.sim, n_single.plots, fa.plots.sim$value, fa.plots, max.dom.sims$value, max.dominance.plots)

wilcoxresults = c(wilcox.test(n_spec.plots.sim$value, n_spec.plots)$p.value,wilcox.test(n_single.sim, n_single.plots)$p.value, wilcox.test(fa.plots.sim$value,fa.plots)$p.value)

names(output[["Guyana/Suriname"]]) = c("n_spec.plots.sim","n_spec.plots","n_single.sim","n_single.plots","fa.plots.sim","fa.plots","max.dom.sims","max.dominance.plots")

test = melt(output)

test$ana = "na"
for(i in 1:length(test$L2)){
if(grepl("n_spec.*",test[i,]$L2)==T){test[i,]$ana = "Number of Species"}
if(grepl("n_single.*",test[i,]$L2)==T){test[i,]$ana = "Number of Singletons"}
if(grepl("fa.plots.*",test[i,]$L2)==T){test[i,]$ana = "FA per plot"}
if(grepl("max.do*",test[i,]$L2)==T){test[i,]$ana = "Max Dominance"}
	}

test$variable = "na"
for(i in 1:length(test$L2)){
if(grepl("*.sim",test[i,]$L2)==T){test[i,]$variable = "Sim"}
	else{test[i,]$variable = "Field"}
	}

test$inter = interaction(test$L1,test$variable)

pdf("Figure 2.pdf", paper="special", height=20, width=40)
group.colors <- c(Sim = "red", Field = "darkgreen")
demoplot<-ggplot(test,aes(x=L1,y=value))
demoplot+geom_boxplot(aes(fill=variable),position=position_dodge(1), outlier.color=NA, show.legend=F)+
facet_wrap(~ana, nrow=1)+
scale_fill_manual(values=group.colors)+
ylim(0,300)+
theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_text(angle=60, hjust=1), text = element_text(size=15))
dev.off()

######################
#### END FIGURE 2 ####
######################



#################################################
#### FIG. 3 INDIVIDUAL VARIANCE IN ABUNDANCE ####
#################################################

#Create plotdata, respecify nr ind when reading in parallel forests, otherwise it does not work
Ind	= ceiling(mean(rowSums(vegtab.field))) 
Plotdata.3D(Forest_complete)

iterations = 10
vegtablist = list()
library(boot)

for(j in 1:iterations){
samplesize = length(rownames(vegtab.field))
invisible(capture.output(create.data3D(Forest_complete, sample = T, n=samplesize)))
vegtablist[j] = list(vegtab.sim3D)
}

##Rank abundance curve simulated data
vegtab.sim = vegtab.sim3D
n.ind.spec.sim	=	colSums(vegtab.sim)
n_ind.plots.sim	= rowSums(vegtab.sim) 
ME.RAD.sim	=	rev(sort(n.ind.spec.sim))

##Rank abundance curve Fielddata
n.ind.spec.field	=	colSums(vegtab.field)
n_ind.plots.field	= rowSums(vegtab.field) 
ME.RAD.field		=	rev(sort(n.ind.spec.field))

##Now for the total set
after.nrspecies =  unique(levels(factor(Forest_complete)))
total = matrix(nrow=8000,ncol=length(after.nrspecies), 0)
colnames(total) = after.nrspecies

for(i in 1:length(vegtablist)){
total[as.numeric(rownames(as.matrix(vegtablist[[i]]))),
colnames(as.matrix(vegtablist[[i]]))] <- as.matrix(vegtablist[[i]])}

total <- total[, colSums(total) != 0]  # delete all the zero columns
total <- total[rowSums(total) != 0,]  # delete all the zero rows

total.sub = total[sample(nrow(total),size=length(rownames(vegtab.field)),replace=F),]
total.sub <- total.sub[, colSums(total.sub) != 0]  # delete all the zero columns
total.sub <- total.sub[rowSums(total.sub) != 0,]  # delete all the zero rows

n.ind.spec.sub	=	colSums(total.sub)
ME.RAD.sub	=	rev(sort(n.ind.spec.sub))

## Do the bootstrap
boot.mean<-function(x,i){boot.mean<-mean(x[i])}
bootmatrix = apply(total.sub/sum(total.sub),2,function(y){ 
   b<-boot(y,boot.mean,R=100); 
   c(mean(b$t),boot.ci(b,type="perc", conf=0.95)$percent[4:5])
})

bootmatrix = bootmatrix[,rev(order(bootmatrix[1,]))]

mean.abund.field = colMeans(vegtab.field/sum(vegtab.field))
sorted.mean.field = rev(sort(mean.abund.field))
mean.abund.sim = colMeans(total.sub/sum(total.sub))
sorted.mean.sim = rev(sort(mean.abund.sim))

pdf("Figure 3.pdf", paper="special",width=12)
options(scipen=3)
plot(RD.spp,(RAD.RD/sum(RAD.RD))*100, log="y",type="n",xlab="Rank", ylab="Relative Abundance (%)", ylim=c(0.000001,1), xlim = c(0,2200), cex.lab=1.4, cex.axis=1.4)

arrows(1:length(sorted.mean.sim*100),sorted.mean.sim*100,1:length(sorted.mean.sim*100),(bootmatrix[2,]+10e-8)*100, length=0.02, angle=90, col=alpha("red", 0.6))
arrows(1:length(sorted.mean.sim*100),sorted.mean.sim*100,1:length(sorted.mean.sim*100),bootmatrix[3,]*100, length=0.05, angle=90, col=alpha("red", 0.6))

lines(RAD.RD/sum(RAD.RD), lty=2)
lines(sorted.mean.sim*100, col="black")
lines(sorted.mean.field*100, col="darkgreen")

legend("topright",c("Observed", "Expected, metacommunity","Expected, bootstrapped", "mean \u00b1 95% confidence limits "),col=c("darkgreen", "black","black"), lty=c(1,2,1, cex=0.8), bty="n", cex=1.2)
legend("topleft", c("Guyana/Suriname"), bty="n", cex=1.4)
dev.off()

####################
### END FIGURE 3 ###
####################



#################################################################################################
##### FIGURE 4 USES THE SAME SCRIPT AS FOR FIGURE 1 BUT WITH DIFFERENT SIMULATED DATASETS	#####
##### 1) THE MODEL WITH MIGRATION SET NEAR UNITY AND 2) MIGRATION SET NEAR NULL				#####
#################################################################################################



##################################################################################
#### FIG. 5 CORRELATE SPECIES IDENTITY AND MAX DOMINANCE BETWEEN FOREST TYPES ####
##################################################################################

vegtab.field = vegtab
data_env.field = plotdata

test = as.data.frame(cbind(data_env.field$Forest,apply(vegtab.field,1,function(x) names(vegtab.field)[which(x==max(x))])))

FT1 = "TF"
FT2 = "PZ"

forest1 = data_env.field[data_env.field$Forest == FT1,]
abund1 = vegtab.field[rownames(vegtab.field) %in% forest1$PlotCode,]
abund1 = (abund1/sum(abund1))+1e-8

forest2 = data_env.field[data_env.field$Forest == FT2,]
abund2 = vegtab.field[rownames(vegtab.field) %in% forest2$PlotCode,]
abund2 = (abund2/sum(abund2))+1e-8

abund = data.frame(colSums(abund1)/max(colSums(abund1)))
abund = cbind(abund, colSums(abund2)/max(colSums(abund2)))

SpeciesCol = vector(length(colnames(vegtab.field)), mode="character")
SpeciesCol[colnames(vegtab.field) %in% test$V2] = "red"
SpeciesCol[!colnames(vegtab.field) %in% test$V2] = "black"

tiff(file="figure_5.tiff",width = 10, height = 6, units = "in", res=400)
plot(abund[,1],abund[,2], xlab="relative max abundance on Terra Firme", ylab="relative max abundance on Podzol",pch=16, cex=.7, col=SpeciesCol, log="xy", xlim=c(1e-2,1), ylim=c(1e-2,1), cex.lab=1.4, cex.axis=1.4)
set = abund[abund[,1]&abund[,2]>1e-4,]
cor = cor.test(log(set[,1]),log(set[,2]), method = "pearson")
p = as.numeric(format(cor$p.value, format = "e", digits = 2))
est = round(cor$estimate, 2)

legend("bottomright",cex=1.4,c("Guyana/Suriname", paste("Pearson R2 = ",est), paste("p = ",p)), bty="n")
abline(0,1)
dev.off()


#################################################
### TABLE 1: COMPARISON BETWEEN FIELD AND SIM ###
#################################################

details.table = matrix(nrow=7,ncol=2)
rownames(details.table) = c("mean_nr.species","mean_nr.ind","total nr.single","mean nr.single","mean fa.plots","total.nr.species","total.fa")
colnames(details.table) = c("simulation","field data")

#Some details of the model
n_spec.plots.sim		=	rowSums(vegtab.sim3D>0)
n_ind.plots.sim		=	rowSums(vegtab.sim3D)
n_single.sim			=	sum(colSums(vegtab.sim3D) == 1)
n_single.plots.sim	=	rowSums(vegtab.sim3D == 1)
mean_single_sim   	=	mean(rowSums(vegtab.sim3D == 1))
fa.plots.sim			=	mapply(fishers.alpha,n_ind.plots.sim,n_spec.plots.sim)
fa.plots_sim 		=	fa.plots.sim
species				=	sum(colSums(vegtab.sim3D)>0)
total.fa				=   fishers.alpha(sum(vegtab.sim3D),species)
max.dominance.sim	=	apply(vegtab.sim3D,1,max)

details.table[,1] = c(mean(n_spec.plots.sim), mean(n_ind.plots.sim),
                    n_single.sim,mean_single_sim,mean(fa.plots.sim), species, total.fa)

#Compare with the field data
vegtab.field			=	vegtab
n_spec.plots			=	rowSums(vegtab.field>0)
n_ind.plots			=	rowSums(vegtab.field)
n_single				=	sum(colSums(vegtab.field) == 1)
n_single.plots		=	rowSums(vegtab.field == 1)
mean_single_field   =	mean(rowSums(vegtab.field == 1))
fa.plots				=	mapply(fishers.alpha,n_ind.plots,n_spec.plots)
fa.plots_field 		=	fa.plots
species				=	length(colnames(vegtab.field))
total.fa				=   fishers.alpha(sum(vegtab.field),species)
max.dominance.plots	=	apply(vegtab.field,1,max)

details.table[,2] = c(mean(n_spec.plots), mean(n_ind.plots),n_single,
                      mean_single_field,mean(fa.plots), species, total.fa)

details.table

###################
### END TABLE 1 ###
###################


#####################################################################
#### SA FIG. S3: Ordination, Non Metric Multidimensional Scaling  ###
#####################################################################

##############################################
## SELECTING DATASET AND GENERAL SETTINGS ####
##############################################

envset = plotdata
vegtabset = vegtab
vegtabset = vegtabset[rownames(vegtabset) %in% envset$PlotCode,]
vegtabset = vegtabset[,colSums(vegtabset) != 0]
nr.plots = length(rownames(vegtabset))

#Create vector for coloring, use # of plots
#Coloring the forest types
ForestCol = vector(nr.plots, mode = "character")
#Assign colors to values from vector Forestcol using data_plots$Forest
ForestCol[envset$Forest == "PZ"] = "orange"
ForestCol[envset$Forest == "TF"] = "red"

#Coloring the countries of the plots Guyana/Suriname
CountryCol = vector(nr.plots, mode = "character")
#Assign colors to values from vector CountryCol using data_plots$Country
CountryCol[envset$Country == "Guyana"] = "red"
CountryCol[envset$Country == "Suriname"] = "blue"

#Create ordination and show stressplot
ord <- metaMDS(vegdist(vegtabset, distance="morisita"), trymax=100)
stressplot(ord);ord

##Create, plot and add points of forest type with coloring to previous framework and add polygon
plot(ord, disp="sites", type="n", main = "Guyana/Suriname")
points(scores(ord)[,1], scores(ord)[,2], pch = 21, bg = CountryCol, col = CountryCol)

##Add polygons with colour
treat=envset$Forest; colors=ForestCol
for(i in unique(treat)) {
  ordihull(ord$point[grep(i,treat),],draw="polygon",alpha=50, 												groups=treat[treat==i],col=colors[grep(i,treat)],label=F)}

text(-.04,.15, "TF")
text(-.20,-.3, "PZ")

##Add polygons with colour
treat=envset$Country; colors=CountryCol
for(i in unique(treat)) {
  ordihull(ord$point[grep(i,treat),],draw="lines",alpha=50, lty="dashed",												groups=treat[treat==i],col=colors[grep(i,treat)],label=F)}

#Create legend and use same order of levels Foresttype
legend("bottomright",legend = unique(envset$Country), bty="n", col = unique(CountryCol), pch = 21, pt.bg=unique(CountryCol), cex = 0.8, horiz = T)

#ANOVA
country = envset$Country
type = envset$Forest

nmds.anova = lm(scores(ord)[,2]~type)
summary(nmds.anova)
anova(nmds.anova)

###########################
#### END ANALYSIS NMDS ####
###########################


################################################################
#### SA Fig. S6 Proportion of co-occurring dominant species ####
################################################################

test = as.data.frame(cbind(data_env.field$Longitude,data_env.field$Latitude,apply(vegtab.field,1,function(x) names(vegtab.field)[which(x==max(x))])))

#Set long lat data to calculate distances
longlats <- data.frame(long = as.numeric(test$V1), lat = as.numeric(test$V2))
dists <- CalcDists(longlats)
#convert to kilometers from meters
disttest = dists/1000

check = matrix(0,length(rownames(test)),length(rownames(test)))
for(i in 1:length(rownames(test))){
	for (j in 1:length(rownames(test))){
if (as.character(test[i,3]) == as.character(test[j,3])) {
	check[i,j] =1} else{check[i,j] = 0}
}}

check  = as.dist(check)

##Proportional co-occuring dominance over distance
distance_bins = seq(100,ceiling(max(disttest))+100, 100)
cum_dist = data.frame()
for(i in 1:length(distance_bins)){
dist_subset = disttest[disttest<distance_bins[i]]
proportion = sum(check[dist_subset])/length(check[dist_subset])
if(i==1){cum_dist = cbind(distance_bins[i], proportion)}
else{cum_dist = rbind(cum_dist,cbind(distance_bins[i], proportion))}}

barplot(cum_dist[,2], names.arg = distance_bins, ylab = "Proportional co-occuring dominance", xlab="Separation in km", ylim=c(0,(max(cum_dist[,2]))+.05), main = "Guyana/Suriname")


##############################################################
#### SA: TABLE S1 MANTEL CORRELATIONS DISPERSAL SYNDROMES ####
##############################################################

specdis		= vegdist(vegtab, method = "morisita")
plotcoord	= subset(plotdata, select = c(Longitude, Latitude))
plottype	= subset(plotdata, select = Forest)

coorddis = vegdist(plotcoord, method="euclidian")
typedist = daisy(plottype, metric = "gower")

mantel(coorddis, specdis, permutations = 1000)
mantel(typedist, specdis, permutations = 1000)

#Data regarding dispersal syndromes was obtained by courtesy of Pablo Stevenson and as such is not incorporated in this script, main analyses are the same using the above scripts.

#############################
#### END MANTEL ANALYSIS ####
############################


