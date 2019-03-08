##########################################################################
#### FUNCTIONS FOR 3D SPATIAL MODEL ######################################
##########################################################################


##Function to create a 3D forest
Createforest.3d = function(Np, Ind, mixed = F){

pb <- txtProgressBar(min = 0, max = Np*Np, style = 3)

seq		= rep(1:Np, each = Np)
seq2	= rep(1:Np, times = Np)
	Forest <<- array(0, dim=c(Np,Np,Ind))
	if (mixed == F)	for (i in 1:(Np*Np)){
		Forest[seq[i],seq2[i],1:Ind]	<<-	RD.spp[i]
		setTxtProgressBar(pb, i)	
		}
	else for (i in 1:(Np*Np)){
		for (j in 1:Ind){
	Forest[seq[i],seq2[i],j]	<<-	sample(RD.spp,	1, prob = RAD.RD)
	}
# update progress bar
setTxtProgressBar(pb, i)	
}
close(pb)
}

Plotdata.3D = function(data){
	Plotdatafull = array2df(data)
	colnames(Plotdatafull) = c("Species","x","y","z")
	Plotdatafull = Plotdatafull[order(Plotdatafull$x,Plotdatafull		$y,Plotdatafull$z),]
	Plotdatafull$Plot = rep(1:(dim(data)[1]*dim(data)[2]), each = Ind)
	Plotdatafull <<- Plotdatafull[,c("x","y","z","Plot","Species")]
	Plotdata = subset(Plotdatafull, select=c("Plot","Species"))
	Plotdata["Abund"] = 1
	Plotdata = as.data.table(Plotdata)
	Plotdata = Plotdata[, sum(Abund), by = c("Plot", "Species")]
	colnames(Plotdata) = c("Plot","Species","Abund")
	Plotdata = Plotdata[order(Plotdata$Plot,Plotdata$Species),]
	Plotdata <<- as.data.frame(Plotdata)
	}


##Function to create data for analyses of 3D Forest
create.data3D= function(data, sample = F, n = NULL){
  #Create vegetation table
  print("creating vegtab")
if (sample == F){
	vegtab.sim3D <<- matrify2(Plotdata)
	print(c("All species are in vegtab", length(unique(Plotdata$Species)) == 							length(colnames(vegtab.sim3D))))
	print(c("All plots are in vegtab", length(rownames(vegtab.sim3D)) == 
	length(unique(Plotdata$Plot))))
print("creating coordinate system")
	coordinates = unique(Plotdatafull[,c("x","y")])
	coordinates = coordinates[ order(as.numeric(row.names(coordinates))), ]
print("creating distance matrices")	
	specdis.sim3D <<- vegdist(vegtab.sim3D, method="bray")
	coorddis.sim3D <<- vegdist(coordinates, method = "euclidean")*100
print("done")}


if (sample == T){
	plots_sample = sample(unique(Plotdata$Plot),n)
	Plotdata_sample <<- Plotdata[Plotdata$Plot %in% plots_sample,]
	vegtab.sim3D <<- matrify2(Plotdata_sample)
	print(c("All species are in vegtab", length(unique(Plotdata_sample$Species)) == 					length(colnames(vegtab.sim3D))))
	print(c("All plots are in vegtab", length(rownames(vegtab.sim3D)) 	== 								length(unique(Plotdata_sample$Plot))))
print("creating coordinate system")
	Plotdatafull_sample = Plotdatafull[Plotdatafull$Plot %in% plots_sample,]
	coordinates = unique(Plotdatafull_sample[,c("x","y")])
	coordinates = coordinates[ order(as.numeric(row.names(coordinates))), ]
print("creating distance matrices")	
	specdis.sim3D <<- vegdist(vegtab.sim3D, method="bray")
	coorddis.sim3D <<- vegdist(coordinates, method = "manhattan")*100
print("done")}}


###############################################
#### FUNCTION FOR CALCULATION OF DISTANCES ####
###############################################

# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

# Calculates the geodesic distance between two points specified by 
# radian latitude/longitude using the Haversine formula
# Ouputs distance between sites 1 and 2 as meters
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = (R * c)*1000
  return(d) # Distance in meters
}

# Fxn to calculate matrix of distances between each two sites
# INPUT: a data frame in which longs are in first column and lats in second column
# OUTPUT: a distance matrix (class dist) between all pairwise sites
# Output distances are in meters
CalcDists <- function(longlats) {
	name <- list(rownames(longlats), rownames(longlats))
	n <- nrow(longlats)
	z <- matrix(0, n, n, dimnames = name)
    for (i in 1:n) {
        for (j in 1:n) z[i, j] <- gcd.hf(long1 = deg2rad(longlats[i, 2]), 
            lat1 = deg2rad(longlats[i, 1]), long2 = deg2rad(longlats[j, 2]), 
            lat2 = deg2rad(longlats[j, 1]))
    }
	z <- as.dist(z)
    return(z)
}
