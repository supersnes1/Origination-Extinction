#function needs to search between intervals and pulls out changing taxa
#read in objects

#setwd("~/Dropbox/ungulate_RA/RCode/Orignation_Extinction_Results")
#setwd("~/Dropbox/ungulate_RA/RCode/Orignation_Extinction_Results")
#install.packages("ggplot2")
#library(ggplot2)

#load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=TRUE##------ Mon Jan 14 17:48:35 2019 ------##.Rdata")

#check jon ks code
	#ksMatrix()
	#pairwiseKSTestsSubsequent

getUngulateOnly <- function(repIntSp, bigList, shortFam, thisMat)
{
	#Need to remove non-target taxonomic designations and/or species (i.e. family and genera names, species not in target clades)
	#testList <- lapply(repIntSp[[1]], function(y) sapply(y, function(x) x[x %in% bigList[as.character(bigList$family) %in% shortFam,]$accepted_name], simplify="array"))
	
	###Crop dataset to onnly be ungulates
	repCulled_Int <- list()
	repCulled_All <- list()
	ungOrder <- c("Artiodactyla", "Perissodactyla")
	#repIntSp[[1]][[1]][repIntSp[[1]][[1]] %in% bigList$accepted_name[as.character(bigList$family) %in% shortFam]]
	for(xx in seq_len(length(repIntSp))){
		for(yy in seq_len(length(repIntSp[[xx]]))){
			#intCulled <- repIntSp[[xx]][[yy]][repIntSp[[xx]][[yy]] %in% bigList$accepted_name[as.character(bigList$family) %in% shortFam & bigList$order == c("Artiodactyla", "Perissodactyla")]]
			intCulled <- repIntSp[[xx]][[yy]][repIntSp[[xx]][[yy]] %in% bigList$accepted_name[as.character(bigList$family) %in% shortFam]]
			
			#remove species lacking a body mass
			intCulled <- intCulled[intCulled %in% thisMat$species]
			#remove non-ungulate species
			intCulled <- intCulled[intCulled %in% bigList$accepted_name[as.character(bigList$order) %in% ungOrder]]
			
			repCulled_Int[[yy]] <- intCulled
		}
		names(repCulled_Int) <- names(repIntSp[[xx]])
		repCulled_All[[xx]] <- repCulled_Int
		#	print(xx)
	}
	return(repCulled_All)
}

#change this to handle a single instance of the repIntSp, 
#to get full plot make sure to run using an apply
findOrigExt_Breaks <- function(repIntSpSingle, optList, measureMat, bigList, shortFam, breaks = NULL, useBreaks = TRUE, startDate=56,endDate=3,plot.check = FALSE,plot.type = "lines")
{
		if(is.null(breaks)){
			print("Get breaks\n")
			#locate shifts (get optimal breaks)
			optBreaks <- 9999 
			for(ii in seq_len(length(optList))){
				#breaks <- optList[[length(optList)]]$optBreaks
				if(optBreaks > optList[[ii]]$AICc) {
					optBreaks <- optList[[ii]]$AICc
					#print(optList[[ii]]$AICc)
					breaks <- optList[[ii]]$optBreaks
				}
			}
		}
		names(repIntSpSingle) <- gsub(pattern = " Ma", replacement = "", x = names(repIntSpSingle))
		interval.names <- vector()
		taxList <- list()
		
		print("Bin taxa \n")
		#bin taxa within each shift
		for(nn in seq(1,length(breaks) + 1,1)){
			# print(nn)
			#locate intervals
			if(nn == 1){
				#startDate <- 56
				topDate <- startDate
				bottomDate <- breaks[1]
				#print("First break)
			} 
			if(nn == length(breaks) + 1){
				topDate <- breaks[nn-1]
				bottomDate <- endDate
				#print("Last break")
			} 
			if(nn > 1 & nn < length(breaks) + 1){
				topDate <- breaks[nn-1]
				bottomDate<- breaks[nn]
				#print(nn)  
			}
			
			print("Extract Species\n")
			#extract species
			tax.vec <- vector()
			
			for(xx in seq_len(length(repIntSpSingle))){
			#	for(kk in as.numeric(sub("* Ma", "", names(repIntSpSingle[[xx]])))){
					if(xx < topDate & xx > bottomDate){
						tax.vec <- append(tax.vec, repIntSpSingle[[xx]]) #[[kk]])
						#print(paste(nn,paste(kk, paste(topDate, bottomDate))))
						#print(tax.vec)
					}
					taxList[[nn]] <- unique(tax.vec)
					#run.count <- run.count + 1
			#	}
				interval.names <- append(interval.names, paste(paste(topDate,"_",sep=""),bottomDate,sep=""))
			}
			#print(taxList)
		}
		names(taxList) <- unique(interval.names)
		#need to find way to impliment for each interval/regime

		#optList <-optList_bm_median
		# "all" #c("all", "artio", "perisso")
		#"bm_median" #c("bm_median", , )
		#repIntSp
		
		plot.row <- ceiling((length(taxList)-1)/2)
		OrigExtList <- list()
		thisMat <- measureMat
		
		print("Generate Image File\n")
		#file.png <- paste(paste(file.name,paste("_",paste(timestamp(),".png",sep=""),sep=""),sep=""),sep="")
		#if(plot.check == TRUE){
		#	png(filename= file.png)}
		#par(mfrow=c(plot.row,2), mar=c(2,4,1,0.5))
		
		for(mm in seq(1,length(taxList)-1,1)){
		
			spec.list <- list()
			#find taxa that went extinct
			tax.extinct <- thisMat[taxList[[mm]][!taxList[[mm]] %in% taxList[[mm+1]]],]
			tax.extinct <- tax.extinct[complete.cases(tax.extinct[,"bodyMass"]),]
			spec.list[[1]] <- tax.extinct[,c("bodyMass")]
			
			#find taxa that originate
			tax.originate <- thisMat[taxList[[mm+1]][!taxList[[mm]] %in% taxList[[mm+1]]],]
			tax.originate <- tax.originate[complete.cases(tax.originate[,"bodyMass"]),]
			spec.list[[2]] <-tax.originate[,c("bodyMass")]
			
			names(spec.list) <- c("extinction", "origination")
			
			OrigExtList[[mm]] <- spec.list
		}
		#dev.off()
		
		names(OrigExtList) <- breaks

return(OrigExtList)
}

plotOrigExt <- function(OrigExtList, plot.type = "lines")# plot.OrigExt = c("Orig", "Ext", "Both"))
{
	for(yy in names(OrigExtList[[1]]))
	{
		for(xx in seq(1, length(OrigExtList),1))
		{
			if(xx == 1) 
			{
				#get initial plot window for extinction
				if(plot.type == "points"){
					plot(OrigExtList[[xx]][[yy]]$extinction, col="blue", pch=19, 
							 ylab = "Body Mass", main = breaks[mm])
				}
				if(plot.type == "lines"){
					plot(density(OrigExtList[[xx]][[yy]]$extinction), col="blue", 
							 xlab = "Body Mass", main = yy, ylim=c(0,1))
				}
			}
			if(xx > 1)
			{
				if(plot.type == "points"){
					points(OrigExtList[[xx]][[yy]]$extinction, col="blue", pch=19)}
				if(plot.type == "lines"){
					lines(density(OrigExtList[[xx]][[yy]]$extinction), col="blue")}
			}
		}
	
		for(xx in seq(1, length(OrigExtList),1))
		{
			if(plot.type == "points"){
				points(OrigExtList[[xx]][[yy]]$origination, col="red", pch=19)}
			if(plot.type == "lines"){
				lines(density(OrigExtList[[xx]][[yy]]$origination), col="red")}
		}
	}
}

findOrigExt_Interval <- function(repIntSpSingle, measureMat, bigList, shortFam, intervals = NULL, startDate=56,endDate=3)
{
	
	names(repIntSpSingle) <- gsub(pattern = " Ma", replacement = "", x = names(repIntSpSingle))
	interval.names <- vector()
	taxList <- list()
	
	print("Bin taxa \n")
	#bin taxa within each shift
	for(nn in rev(seq(1,nrow(intervals),1))){
		# print(nn)
		#locate intervals
		if(nn == nrow(intervals)){
			#startDate <- 56
			topDate <- startDate
			bottomDate <- intervals[nn,1]
			#print("First break)
		} 
		if(nn == 1){
			topDate <- intervals[nn+1,1]
			bottomDate <- 1
			#print("Last break")
		} 
		if(nn > endDate & nn < nrow(intervals)){
			topDate <- intervals[nn+1,1]
			bottomDate <- intervals[nn,1] # which one will be the -1/+1
			#print(nn)  
		}
		
		print("Extract Species\n")
		#extract species
		tax.vec <- vector()
		
		for(xx in seq_len(length(repIntSpSingle))){
			#	for(kk in as.numeric(sub("* Ma", "", names(repIntSpSingle[[xx]])))){
			if(xx <= topDate & xx >= bottomDate){
				tax.vec <- append(tax.vec, repIntSpSingle[[xx]]) #[[kk]])
				#print(paste(nn,paste(kk, paste(topDate, bottomDate))))
				#print(tax.vec)
			}
			taxList[[nn]] <- unique(tax.vec)
			#run.count <- run.count + 1
			#	}
			interval.names <- append(interval.names, paste(paste(topDate,"_",sep=""),bottomDate,sep=""))
		}
		#print(taxList)
	}
	names(taxList) <- rev(unique(interval.names))
	#need to find way to impliment for each interval/regime
	
	#optList <-optList_bm_median
	# "all" #c("all", "artio", "perisso")
	#"bm_median" #c("bm_median", , )
	#repIntSp
	
#	plot.row <- ceiling((length(taxList)-1)/2)
	OrigExtList <- list()

	print("Generate Image File\n")
	#file.png <- paste(paste(file.name,paste("_",paste(timestamp(),".png",sep=""),sep=""),sep=""),sep="")
	#if(plot.check == TRUE){
	#	png(filename= file.png)}
	#par(mfrow=c(plot.row,2), mar=c(2,4,1,0.5))
	
	for(mm in seq(endDate,length(taxList),1)){
		
		#need to think on how to fill the output file
		## origination should fill all but first bin whereas extinction fills all but last bin
		##to get at the fact that they are offset
		
		spec.list <- list()
		#find taxa that went extinct
		#tax.extinct <- thisMat[taxList[[mm]][!taxList[[mm-1]] %in% taxList[[mm]]],]
		#tax.extinct <- tax.extinct[complete.cases(tax.extinct[,"bodyMass"]),]
		#spec.list[[1]] <- tax.extinct[,c("bodyMass")]
		
			#find taxa that originate
			tax.originate <- thisMat[taxList[[mm-1]][!taxList[[mm]] %in% taxList[[mm-1]]],]
			tax.originate <- tax.originate[complete.cases(tax.originate[,"bodyMass"]), ]#"bodyMass"]
			#spec.list <-tax.originate[,c("bodyMass")]
			
			#names(spec.list) <- c("origination")
			
			OrigExtList[[mm]] <- tax.originate
		
	}
	#dev.off()
	
	names(OrigExtList) <- rev(intervals[nrow(intervals):1,2]-.5)
	
	return(OrigExtList)
}

###########################################################################
###########################################################################
###########################################################################

getBm4List <- function(repIntSp, thisMat)
{
	repBM_Int <- list()
	repBM_All <- list()
	if(length(repIntSp) <= 60){
		for(xx in seq_len(length(repIntSp))){
				intBM <- thisMat[thisMat$species %in% repIntSp[[xx]],"bodyMass"]
				#thisMat[thisMat %in% repIntSp[[mm]],"bodyMass"]
			#	if(intBM == numeric(0))
				repBM_Int[[xx]] <- intBM
				repBM_All <- repBM_Int
		}
	}
	if(length(repIntSp) > 60){
		for(xx in seq_len(length(repIntSp))){
			for(yy in seq_len(length(repIntSp[[xx]]))){
				
				intBM <- thisMat[thisMat$species %in% repIntSp[[xx]][[yy]],"bodyMass"]
				#thisMat[thisMat %in% repIntSp[[mm]],"bodyMass"]
				repBM_Int[[yy]] <- intBM
			}
			names(repBM_Int) <- names(repIntSp[[xx]])
			repBM_All[[xx]] <- repBM_Int
			#	print(xx)
		}
	}
	return(repBM_All)
}

getKS4IntervalsOrig <- function(repOrigBmUng, repBmUng,result.out = c("SigInt", "AllInt"), out.putType = c("Intervals","Breaks"))
{
	sig.IntervalsKS <- matrix(nrow = 2, ncol = length(repOrigBmUng))
	all.IntervalsKS <- matrix(nrow = 2, ncol = length(repOrigBmUng))
	if(out.putType == "Intervals"){
		for(pp in seq(1, length(repOrigBmUng),1))
		{
			#print(pp)
			ks.result <- ks.test(repBmUng[[pp]],repOrigBmUng[[pp]]) #need to set the repBMUng to have the species in repOrigBMUng removed to compare origianting species vs species that survived from the last interval
			#ks.result <- ks.test(repBmUng[[pp+1]],repOrigBmUng[[pp]])
			inputs <- as.matrix(c(as.numeric(ks.result$p.value), as.numeric(ks.result$statistic)))
			colnames(inputs) <- paste(pp+1, paste("to ", pp, sep=""),sep=" ")
			
			all.IntervalsKS[,pp] <-inputs
			
			if(ks.result$p.value < 0.01) sig.IntervalsKS[,pp] <- inputs
		}
	}
	if(out.putType == "Breaks"){
		for(pp in seq(1, length(repOrigBmUng),1))
		{
			ks.result <- ks.test(repBmUng[[pp]],repOrigBmUng[[pp]]) #removed pp+1 from repBmUng brackets
			inputs <- as.matrix(c(as.numeric(ks.result$p.value), as.numeric(ks.result$statistic)))
			colnames(inputs) <- paste(pp+1, paste("to ", pp, sep=""),sep=" ")
			
			all.IntervalsKS[,pp] <-inputs
			
			if(ks.result$p.value < 0.01) sig.IntervalsKS[,pp] <- inputs
		}
	}
	rownames(sig.IntervalsKS) <- rownames(all.IntervalsKS) <- c("P-value", "Statstic")
	
	if( result.out == "SigInt")	return(sig.IntervalsKS)
	if( result.out == "AllInt")	return(all.IntervalsKS)
}

getKS4IntervalsExt <- function(repExtBmUng, repBmUng,result.out = c("SigInt", "AllInt"), out.putType = c("Intervals","Breaks"))
{
	sig.IntervalsKS <- matrix(nrow = 2, ncol = length(repExtBmUng))
	all.IntervalsKS <- matrix(nrow = 2, ncol = length(repExtBmUng))
	if(out.putType == "Intervals"){
		for(pp in seq(1, length(repExtBmUng),1))
		{
			ks.result <- ks.test(repBmUng[[pp]],repExtBmUng[[pp]])
			inputs <- as.matrix(c(as.numeric(ks.result$p.value), as.numeric(ks.result$statistic)))
			colnames(inputs) <- paste(pp+1, paste("to ", pp, sep=""),sep=" ")
			
			all.IntervalsKS[,pp] <-inputs
			
			if(ks.result$p.value < 0.01) sig.IntervalsKS[,pp] <- inputs
		}
	}
	if(out.putType == "Breaks"){
		for(pp in seq(1, length(repExtBmUng),1))
		{
			ks.result <- ks.test(repExtBmUng[[pp]], repBmUng[[pp]])
			inputs <- as.matrix(c(as.numeric(ks.result$p.value), as.numeric(ks.result$statistic)))
			colnames(inputs) <- paste(pp+1, paste("to ", pp, sep=""),sep=" ")
			
			all.IntervalsKS[,pp] <-inputs
			
			if(ks.result$p.value < 0.01) sig.IntervalsKS[,pp] <- inputs
		}
	}
	rownames(sig.IntervalsKS) <- rownames(all.IntervalsKS) <- c("P-value", "Statstic")
	
	if( result.out == "SigInt")	return(sig.IntervalsKS)
	if( result.out == "AllInt")	return(all.IntervalsKS)
}

getBreakTaxList <- function(repIntSpSingle, breaks, startDate = 56, endDate = 3)
{
	taxList <- list()
	
	print("Bin taxa \n")
	#bin taxa within each shift
	for(nn in seq(1,length(breaks) + 1,1)){
		# print(nn)
		#locate intervals
		if(nn == 1){
			#startDate <- 56
			topDate <- startDate
			bottomDate <- breaks[1]
			#print("First break)
		} 
		if(nn == length(breaks) + 1){
			topDate <- breaks[nn-1]
			bottomDate <- endDate
			#print("Last break")
		} 
		if(nn > 1 & nn < length(breaks) + 1){
			topDate <- breaks[nn-1]
			bottomDate<- breaks[nn]
			#print(nn)  
		}
		
		print("Extract Species\n")
		#extract species
		tax.vec <- vector()
		
		for(xx in seq_len(length(repIntSpSingle))){
			#	for(kk in as.numeric(sub("* Ma", "", names(repIntSpSingle[[xx]])))){
			if(xx < topDate & xx > bottomDate){
				tax.vec <- append(tax.vec, repIntSpSingle[[xx]]) #[[kk]])
				#print(paste(nn,paste(kk, paste(topDate, bottomDate))))
				#print(tax.vec)
			}
			taxList[[nn]] <- unique(tax.vec)
			#run.count <- run.count + 1
			#	}
			interval.names <- append(interval.names, paste(paste(topDate,"_",sep=""),bottomDate,sep=""))
		}
		#print(taxList)
	}
	#names(taxList) <- unique(interval.names)
	return(taxList)
}

check4num0 <- function(repIntSp, repNum, out.type = c("RepsList","num0List"))
{
	count <- 1
	num0 <- matrix(ncol = 2)
	int.numeric.is0 <- list()
	for(gg in seq(1,length(repIntSp),1))
	{
		if(length(repIntSp[[gg]]) == 0)
		{
			num0 <- c(repNum,gg)
			if(is.numeric(num0)) 
			{
				int.numeric.is0[[count]] <- num0
				count <- count + 1
				rm(num0)
			}
			repIntSp[[gg]] <- 0
		}
	}	

	if(out.type =="RepList") return(repIntSp)
	if(out.type =="num0List") return(int.numeric.is0)
}

getListNum0Ints <- function(numList)
{
	numVec <- unlist(numList)
	if(is.null(numVec)) 
	{
		numMat <- matrix(c(0,0), nrow=1,ncol = 2)
		return(numMat)
	}
	count <- 1
	seqVec<- seq(1, length(numVec)-1,2)
	numMat <- matrix(nrow=length(seqVec),ncol = 2)
	for(yy in seq(1, length(numVec)-1,2))
	{
		numMat[count,] <- c(numVec[yy], numVec[yy+1])
		count <- count +1
	}
	return(numMat)
}

checkRangeThrough1MaBins <- function(repIntTest, ints)
{
	#function to find if rangethrough is working
	
	species <- list()
	for(xx in ints) species[[xx-(min(ints)-1)]] <- repIntTest[[xx]]
	uni.sp <- unique(unlist(species))
	
	checkMat <- matrix(nrow = length(uni.sp), ncol = length(ints))
	rownames(checkMat) <- uni.sp
	colnames(checkMat) <- ints+0.5
	
	# find way to parse through each species and mark if they are in a given interval
	for(xx in seq(1, nrow(checkMat),1))
	{
		for(yy in seq(1, ncol(checkMat),1))
		{
			#for(zz in repIntTest[[yy+49]])
			##{
			test.true <- which(rownames(checkMat)[xx] == repIntTest[[yy+(min(ints)-1)]])
			if(length(test.true) == 0)
			{
				checkMat[xx,yy] <- 0
			}
			if(length(test.true) >= 1) checkMat[xx,yy] <- 1
			#}
		}
	}
	return(checkMat)
}

rangeThrough_alt <- function(IntMat)
{
	#IntMat <- repIntCheck
	#IntMat.rangeThrough <- IntMat
	for(xx in seq_len(nrow(IntMat)))
	{
		int.index <- which(IntMat[xx,] == 1)
		int.new <- min(int.index):max(int.index)
		for (yy in int.new)
		{
			IntMat[xx,yy] <- 1
		}
	}
	
	repIntRangeThrough <- list()
	
	for(gg in seq_len(ncol(IntMat)))
	{
		specVec <- vector()
		for(hh in seq_len(nrow(IntMat)))
		{
			if(IntMat[hh,gg] == 1) specVec <- append(specVec, rownames(IntMat)[hh])
		}
		repIntRangeThrough[[gg]] <- specVec
	}
	names(repIntRangeThrough) <- colnames(IntMat)
	return(repIntRangeThrough)
}

checkRangeThrough <- function(repIntTest, ints)
{
	#function to find if rangethrough is working
	int.label<- vector()
	species <- list()
	for(xx in seq(1, length(ints)-1,1)) 
	{
		species[[xx]] <- repIntTest[[xx]] #species[[xx-(min(ints)-1)]] 
		int.label[xx] <- (ints[xx]+ints[xx+1])/2
	}
	uni.sp <- unique(unlist(species))
	
	checkMat <- matrix(nrow = length(uni.sp), ncol = length(ints)-1)
	rownames(checkMat) <- uni.sp
	colnames(checkMat) <- int.label
	
	# find way to parse through each species and mark if they are in a given interval
	for(xx in seq(1, nrow(checkMat),1))
	{
		for(yy in seq(1, ncol(checkMat),1))
		{
			#for(zz in repIntTest[[yy+49]])
			##{
			test.true <- which(rownames(checkMat)[xx] == repIntTest[[yy]]) #repIntTest[[yy+(min(ints)-1)]])
			if(length(test.true) == 0)
			{
				checkMat[xx,yy] <- 0
			}
			if(length(test.true) >= 1) checkMat[xx,yy] <- 1
			#}
		}
	}
	return(checkMat)
}

################################################################################################################
####Plotting functions for talks and publications
################################################################################################################
getOccs <- function()
{
  occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
  occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
  occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
  occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores
  
  return(occs)
}

getEcoAnalSrc <- function()
{
  source("~/Dropbox/code/R/common_src/strat.R")
  source("~/Dropbox/code/R/common_src/occFns.R")
  source("~/Dropbox/code/R/common_src/sampling.R") 
  source("~/Dropbox/code/R/common_src/utils_marcot.R")
  source("~/Dropbox/code/R/common_src/CzTimescale.R") 
  
  source('~/Dropbox/code/R/dentalMeasurements/src/src_dentalDataFns.R', chdir = TRUE)
  source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
  source('~/Dropbox/code/R/dentalMeasurements/src/src_ecologyAnalysisFns.R', chdir = TRUE)
  return()
}

getDataLists_Eco_OrigExt <- function(object = "bigList", occs = NULL, thisMat = NULL)
{
  ##################################################################################################################################
  #### reduces matrix to just the focal order(s)
  ####################################################################################################################################
  
  # focal.order <- "Artiodactyla"
  # focal.order <- "Perissodactyla"
  focal.order <- c("Artiodactyla", "Perissodactyla")
  bigList <- unique(occs[(occs$accepted_rank =="species" & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
  bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
  # bigList[order(bigList$family, bigList$accepted_name),]
  shortFam <- sort(unique(bigList$family[bigList$order %in% focal.order]))	
  if(object == "shortFam") return(shortFam)
  
  bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
  if(object == "bigList") return(bigList)
  
  # matrix(thisMat$species[!thisMat$species %in% bigList$accepted_name[bigList$order %in% focal.order]], ncol=1)
  thisMat <- thisMat[thisMat$species %in% bigList$accepted_name[bigList$order %in% focal.order], ]
  if(object == "thisMat") return(thisMat)
}

plot3panel <- function(repIntSp = NULL, bigList = NULL, shortFam = NULL, intervals = NULL, optList_tax_median = NULL, optList_bm_median = NULL)
{
  #3-panel diagram
  ######
  # make three-panel figure
  #####
  quartz(width=6.89)
  par(mfrow=c(4,1), mar=c(0,4,0.5,0.5), mgp=c(2, 1,0))
  
  ### isotope panel
  if(Sys.info()["sysname"] == "Darwin"){
    source("~/Dropbox/ungulate_RA/RCode/isotopes.R")
  } else if(Sys.info()["sysname"] == "Windows"){
    source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
  }
  
  #source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
  # source("~/Dropbox/code/R/common_src/isotopes.R")
  
  optList_topes <- doTopesRateAnalysis(intervals)
  plotTopesRateAnalysis(optList_topes, intervals, x.axis=FALSE) #
  box(lwd=1)
  getAlroyStatistics(intervals)
  
  ### taxonomy panel
  par(mar=c(0,4, 2.5,0.5))
  
  # taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% rownames(thisMat)[x]], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  #taxCubeG <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$genus) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)

  # prop <- t(apply(taxCube, c(1,2), median, na.rm=TRUE))
  prop <- t(apply(taxCube, c(1,2), mean, na.rm=TRUE))
  colnames(prop)[colnames(prop)==""] <- "indeterminate"
  # dimnames(prop) <- list(rownames(intervals), shortFam)
  source("https://dl.dropbox.com/s/iy0tu983xesbig2/taxonomicEv.R")
  plotStackedRichness(this.box=prop, intervals=intervals, do.log=FALSE, overlay.labels=TRUE, numbers.only=TRUE, legend=FALSE, xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)))
  #med.n <- median(length(unique(unlist(sapply(repIntSp[[this.rep]], function(x) rownames(thisMat)[x]))))) #what is this.rep set to during this function?  variable is used in for loop in Handley
  # med.n <- median(sapply(repIntSp, function(x) length(unique(unlist(sapply(x, function(y) rownames(thisMat)[y]))))))
  # optList_tax <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
  abline(v=sort(c(intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col="darkorchid4")
  text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_tax_median[[length(optList_tax_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="darkorchid4")
  text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0), cex=0.5, col="darkorchid4")
  box(lwd=1)
  
  ### body mass panel
  thisRanges <- getTaxonRangesFromOccs(occs=occs, random=FALSE)
  rownames(thisRanges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisRanges))
  thisMat[,c("FO","LO")] <- thisRanges[match(rownames(thisMat), rownames(thisRanges)),]
  
  par(mar=c(0,4,2.5,0.5))
  # quartz(width=12, height=6)
  plot(thisMat$FO, thisMat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(55,5,10), xlab="Time (Ma)", cex.axis=1.5, cex.lab=1.5)
  # plot(thisMat$FO, thisMat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(50,0,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, col="gray75", fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #alter xaxpto change x-axis values
  # rect(-10e6, -10e6, 10e6, 10e6, col="white")
  overlayCzTimescale(do.subepochs=TRUE)
  
  famColors <- rainbow(length(shortFam))
  colorList <- famColors[match(bigList$family[as.character(bigList$accepted_name) %in% rownames(thisMat)], shortFam)]
  colorList[is.na(colorList)] <- "gray25"
  
  orderColors <- array(NA, dim=nrow(thisMat))
  # orderColors[bigList$order[match(rownames(thisMat), bigList$accepted_name)]=="Perissodactyla"] <- "dodgerblue4"
  # orderColors[bigList$order[match(rownames(thisMat), bigList$accepted_name)] =="Artiodactyla"] <- "deeppink4"
  
  for (i in seq_len(nrow(thisMat))) {
    # lines(x=c(this["FO"], x["LO"]), y=c(x["bodyMass"], x["bodyMass"]), lwd=3, pch=21, col=famColors[match(bigList[match(rownames(thisMat), bigList[,1]),2], shortFam)])
    # lines(x=c(thisMat$FO[i], thisMat$LO[i]), y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor(colorList[i], 0.75))
    # lines(x=c(thisRanges[match(rownames(thisMat)[i], rownames(thisRanges)),"FO"], thisRanges[match(rownames(thisMat)[i], rownames(thisRanges)),"LO"]), y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor("gray0", 0.75)) #
    if (is.finite(thisMat$FO[i]) & is.finite(thisMat$LO[i]) & thisMat$FO[i] != thisMat$LO[i]) lines(x=thisMat[i,c("FO","LO")], y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.75, pch=21, col=alphaColor("gray0", 0.5)) #alphaColor(orderColors[i], 0.5)
  }
  points(thisMat[complete.cases(thisMat[ ,c("FO","LO")]) & thisMat$FO == thisMat$LO, c("FO","bodyMass")], pch=21, col=alphaColor("gray0", 0.5), cex=0.25) #this line is not generating the proper output for the final graph due to c("FO","bodyMass") causing a  "undefined columns selected" error
  
  # optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
  # optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on median
  abline(v=sort(c(intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col="firebrick4")
  text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_bm_median[[length(optList_bm_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="firebrick4")
  text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0),cex=0.5, col="firebrick4")
  
  quants <- apply(sapply(repIntSp, function(y) sapply(y, function(x) quantile(thisMat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
  polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[1,], rev(quants[5,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
  polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[2,], rev(quants[4,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
  lines(rowMeans(intervals), quants[3,], col=alphaColor("goldenrod1", 0.5), lwd=5)
  lines(rowMeans(intervals), quants[3,], col=alphaColor("darkorange4", 1.0), lwd=3)
  points(rowMeans(intervals), quants[3,], col=alphaColor("darkorange1", 0.5), cex=0.5)
  box(lwd=1)
}

plotSinglePanelPlot <- function(panelsel = "bodyMass", includeBreaks = FALSE,
           repIntSp = NULL,
           bigList = NULL,
           shortFam = NULL,
           intervals = NULL,
           optList_tax_median = NULL,
           optList_bm_median = NULL)
  {
    quartz(width = 12, height = 6)
    par(mfrow=c(1,1), mgp = c(3, 1, 0), mar=c(5,5, 2.5,1.5))
    
    if (panelsel == "isotope") {
      ### isotope panel
      if (Sys.info()["sysname"] == "Darwin") {
        source("~/Dropbox/ungulate_RA/RCode/isotopes.R")
      } else if (Sys.info()["sysname"] == "Windows") {
        source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
      }
      
      #source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
      # source("~/Dropbox/code/R/common_src/isotopes.R")
      
      optList_topes <- doTopesRateAnalysis(intervals)
      plotTopesRateAnalysis(optList_topes, intervals, x.axis = FALSE) #
      box(lwd = 1)
      getAlroyStatistics(intervals)
    }
    
    ### taxonomy panel
    if (panelsel == "taxon") {
      #par(mar = c(0, 4, 2.5, 0.5))
      
      # taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% rownames(thisMat)[x]], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
      #taxCubeG <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$genus) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
      taxCube <-
        sapply(repIntSp, function(y)
          sapply(y, function(x)
            tabulate(
              match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins =
                length(shortFam)
            ), simplify = "array"), simplify = "array")
      dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
      
      # prop <- t(apply(taxCube, c(1,2), median, na.rm=TRUE))
      prop <- t(apply(taxCube, c(1, 2), mean, na.rm = TRUE))
      colnames(prop)[colnames(prop) == ""] <- "indeterminate"
      # dimnames(prop) <- list(rownames(intervals), shortFam)
      source("https://dl.dropbox.com/s/iy0tu983xesbig2/taxonomicEv.R")
      plotStackedRichness(
        this.box = prop,
        intervals = intervals,
        do.log = FALSE,
        overlay.labels = TRUE,
        numbers.only = TRUE,
        legend = FALSE,
        xlim = c(max(intervals, na.rm = TRUE), min(intervals, na.rm = TRUE))
      )
      #med.n <- median(length(unique(unlist(sapply(repIntSp[[this.rep]], function(x) rownames(thisMat)[x]))))) #what is this.rep set to during this function?  variable is used in for loop in Handley
      # med.n <- median(sapply(repIntSp, function(x) length(unique(unlist(sapply(x, function(y) rownames(thisMat)[y]))))))
      # optList_tax <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
      if(includeBreaks == TRUE)
      {
        abline(v = sort(c(intervals[optList_tax_median[[length(optList_tax_median) -
                                                          1]]$optBreaks, 2], range(intervals))),
               lwd = 1.5,
               col = "darkorchid4")
        text(
          x = sort((c(
            max(intervals), intervals[optList_tax_median[[length(optList_tax_median) -
                                                            1]]$optBreaks, 2]
          ) - 0.35)),
          y = par()$usr[3],
          labels = rev(seq_len(
            length(optList_tax_median[[length(optList_tax_median) - 1]]$optBreaks) + 1
          )),
          pos = 3,
          cex = 0.5,
          col = "darkorchid4"
        )
        text(
          x = sort((c(
            max(intervals), intervals[optList_tax_median[[length(optList_tax_median) -
                                                            1]]$optBreaks, 2]
          ) - 0.35)),
          y = par()$usr[3],
          labels = paste(sort(c(
            max(intervals), intervals[optList_tax_median[[length(optList_tax_median) -
                                                            1]]$optBreaks, 2]
          )), "Ma"),
          adj = c(0, 0),
          cex = 0.5,
          col = "darkorchid4"
        )
      }
      box(lwd = 1)
    }
    
    ### body mass panel
    if (panelsel == "bodyMass") {
      thisRanges <- getTaxonRangesFromOccs(occs = occs, random = FALSE)
      rownames(thisRanges) <-
        gsub(
          pattern = "[[:space:]]",
          replacement = "_",
          x = rownames(thisRanges)
        )
      thisMat[, c("FO", "LO")] <-
        thisRanges[match(rownames(thisMat), rownames(thisRanges)), ]
      
      #par(mar = c(4, 4, 2, 0.5), oma = c(0,0,2.5,2.5))
      # quartz(width=12, height=6)
      plot(
        thisMat$FO,
        thisMat$bodyMass,
        xlim = c(max(intervals), min(intervals)),
        type = "n",
        ylab = "log-Body Mass (kg)",
        xaxp = c(55, 5, 10),
        xlab = "Time (Ma)",
        cex.axis = 1.5,
        cex.lab = 1.5,
        col.axis = "gray75",
        col.lab = "gray75",
        bg = "gray75",
        fg = "gray75"
      )
      # plot(thisMat$FO, thisMat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(50,0,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, col="gray75", fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #alter xaxpto change x-axis values
      # rect(-10e6, -10e6, 10e6, 10e6, col="white")
      overlayCzTimescale(do.subepochs = TRUE)
      
      famColors <- rainbow(length(shortFam))
      colorList <-
        famColors[match(bigList$family[as.character(bigList$accepted_name) %in% rownames(thisMat)], shortFam)]
      colorList[is.na(colorList)] <- "gray25"
      
      orderColors <- array(NA, dim = nrow(thisMat))
      # orderColors[bigList$order[match(rownames(thisMat), bigList$accepted_name)]=="Perissodactyla"] <- "dodgerblue4"
      # orderColors[bigList$order[match(rownames(thisMat), bigList$accepted_name)] =="Artiodactyla"] <- "deeppink4"
      
      for (i in seq_len(nrow(thisMat))) {
        # lines(x=c(this["FO"], x["LO"]), y=c(x["bodyMass"], x["bodyMass"]), lwd=3, pch=21, col=famColors[match(bigList[match(rownames(thisMat), bigList[,1]),2], shortFam)])
        # lines(x=c(thisMat$FO[i], thisMat$LO[i]), y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor(colorList[i], 0.75))
        # lines(x=c(thisRanges[match(rownames(thisMat)[i], rownames(thisRanges)),"FO"], thisRanges[match(rownames(thisMat)[i], rownames(thisRanges)),"LO"]), y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor("gray0", 0.75)) #
        if (is.finite(thisMat$FO[i]) &
            is.finite(thisMat$LO[i]) &
            thisMat$FO[i] != thisMat$LO[i])
          lines(
            x = thisMat[i, c("FO", "LO")],
            y = c(thisMat$bodyMass[i], thisMat$bodyMass[i]),
            lwd = 0.75,
            pch = 21,
            col = alphaColor("gray0", 0.5)
          ) #alphaColor(orderColors[i], 0.5)
      }
      points(
        thisMat[complete.cases(thisMat[, c("FO", "LO")]) &
                  thisMat$FO == thisMat$LO, c("FO", "bodyMass")],
        pch = 21,
        col = alphaColor("gray0", 0.5),
        cex = 0.25
      ) #this line is not generating the proper output for the final graph due to c("FO","bodyMass") causing a  "undefined columns selected" error
      
      # optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
      # optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on median
      if(includeBreaks == TRUE)
      {
        abline(v = sort(c(intervals[optList_bm_median[[length(optList_bm_median) -
                                                         1]]$optBreaks, 2], range(intervals))),
               lwd = 1.5,
               col = "firebrick4")
        text(
          x = sort((c(
            max(intervals), intervals[optList_bm_median[[length(optList_bm_median) -
                                                           1]]$optBreaks, 2]
          ) - 0.35)),
          y = par()$usr[3],
          labels = rev(seq_len(
            length(optList_bm_median[[length(optList_bm_median) - 1]]$optBreaks) + 1
          )),
          pos = 3,
          cex = 0.5,
          col = "firebrick4"
        )
        text(
          x = sort((c(
            max(intervals), intervals[optList_bm_median[[length(optList_bm_median) -
                                                           1]]$optBreaks, 2]
          ) - 0.35)),
          y = par()$usr[3],
          labels = paste(sort(c(
            max(intervals), intervals[optList_bm_median[[length(optList_bm_median) -
                                                           1]]$optBreaks, 2]
          )), "Ma"),
          adj = c(0, 0),
          cex = 0.5,
          col = "firebrick4"
        )
      }
      quants <-
        apply(sapply(repIntSp, function(y)
          sapply(y, function(x)
            quantile(
              thisMat[x, "bodyMass"],
              probs = c(0, 0.25, 0.5, 0.75, 1.0),
              na.rm = TRUE
            )), simplify = "array"),
          c(1, 2),
          median,
          na.rm = TRUE)
      polygon(
        c(rowMeans(intervals), rev(rowMeans(intervals))),
        c(quants[1, ], rev(quants[5, ])),
        col = alphaColor("darkorange4", 0.25),
        border = "darkorange4"
      )
      polygon(
        c(rowMeans(intervals), rev(rowMeans(intervals))),
        c(quants[2, ], rev(quants[4, ])),
        col = alphaColor("darkorange4", 0.25),
        border = "darkorange4"
      )
      lines(
        rowMeans(intervals),
        quants[3, ],
        col = alphaColor("goldenrod1", 0.5),
        lwd = 5
      )
      lines(
        rowMeans(intervals),
        quants[3, ],
        col = alphaColor("darkorange4", 1.0),
        lwd = 3
      )
      points(
        rowMeans(intervals),
        quants[3, ],
        col = alphaColor("darkorange1", 0.5),
        cex = 0.5
      )
      box(lwd = 1)
    }
    return()
  }

RegimeNetdistribution_Midpoint <- function(thisMat, occs, intervals)
{
  #####Using midpoint of taxon
  #get midpoint of each taxon 
  taxonranges <- getTaxonRangesFromOccs(occs = occs, random = FALSE)
  rownames(taxonranges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(taxonranges))
  taxonranges <- cbind(taxonranges, (taxonranges[,"FO"]+taxonranges[,"LO"])/2); colnames(taxonranges)[3] <- "MP"
  occDatesMP <- taxonranges[,"MP"]
  intSpMP <- apply(intervals, 1, function(thisIntv) taxonranges[taxonranges[,"MP"] > thisIntv[1] & taxonranges[,"MP"] <= thisIntv[2],])
  
  bmBreaks <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, Inf) #Janis 2000  max(thisMat$bodyMass, na.rm=TRUE)
  countCubeMP <- sapply(intSpMP, function(x) hist(thisMat$bodyMass[match(rownames(x), rownames(thisMat))], breaks=bmBreaks, plot=FALSE)$counts, simplify = "array")
  #drop the intervals that include the Quaternary (<3 Ma)
  countCubeMP <- countCubeMP[,!as.double(str_remove(colnames(countCubeMP), " Ma")) < 3]
  optList_bm_medianMP <- doHandleyTest(countCubeMP, n=nrow(thisMat))
  regimeHist_countBox(countBox = countCubeMP, breaks = c(51,47,37,21), optList=optList_bm_medianMP, 
                      thisMat = thisMat, netFreq = TRUE, regimeFreq=FALSE,
                      netPlotType = "pos/neg", plot.together = FALSE, plot.axes = TRUE,
                      grayscale = FALSE)
  return()
}

evoSetup <- function()
{
  source("~/Dropbox/Code/R/common_src/strat.R")
  source("~/Dropbox/Code/R/common_src/occFns.R")
  source("~/Dropbox/Code/R/common_src/phy_dateTree.R")

  source('~/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call cource file for analysis functions
  #source('~/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R') #call cource file for data functions
  
  source("~/Dropbox/Code/R/DentalMeasurements/src/src_dentalDataFns.R")
  source("~/Dropbox/Code/R/DentalMeasurements/src/src_bodyMassEstimation.R")
  source("~/Dropbox/Code/R/DentalMeasurements/src/src_EvAnalysisSetup.R")
  source('~/Dropbox/Code/R/DentalMeasurements/src/src_EvAnalysesTree.R') #call cource file for tree functions
  source('~/Dropbox/Code/R/DentalMeasurements/src/src_EvAnalysesPlot.R') #call cource file for tree functions

  #tree.backbone <- read.nexus("~/Dropbox/ungulate_RA/NAUngulata_Trees/BackBoneTrees/2017_3_24_UngulataBackboneTree")
  #clade.definitions <- read.csv("~/Dropbox/ungulate_RA/2017_3_20_Clade_species_test.csv", stringsAsFactors = FALSE)
  #wildcard.positions <- read.csv("~/Dropbox/ungulate_RA/2017_4_17_MCRA_Codes.csv", stringsAsFactors = FALSE)
  #regressCat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")
  return()
}


#6/18/2019
##for species with elongate m3's (e.g. peccaries, suids, some camels)
#look at variance of tooth row, if m3 inducing the variance drop it from the calculations
#platygonus being predicted to be over 500 lbs currently
