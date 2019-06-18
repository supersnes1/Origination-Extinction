setwd("~/Dropbox/ungulate_RA/RCode/Orignation_Extinction_Results")

# load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/perisso/handleyResult##------ Thu Nov 16 10:27:58 2017 ------##.Rdata")
startTime <- Sys.time()
#sources for Jon Marcot's code and specimen measurements 
source("https://dl.dropbox.com/s/8jy9de5owxj72p7/strat.R")
source("https://dl.dropbox.com/s/253p4avcvb66795/occFns.R")
source("https://dl.dropbox.com/s/9gdafsqss2b586x/phy_dateTree.R")
source("https://dl.dropbox.com/s/9tdawj35qf502jj/amandaSrc.R")
source("https://dl.dropbox.com/s/rlof7juwr2q4y77/blasto_Birlenbach.R")
source("https://dl.dropbox.com/s/pd5nmg1it13noc2/sampling.R") 
source("https://dl.dropbox.com/s/643op7ye4s49w8p/utils_marcot.R")
source("https://dl.dropbox.com/s/dozeb8o2pxu4sbj/CzTimescale.R") 
# source("C:/Users/Evan/Dropbox/ungulate_RA/RCode/isotopes.R")

if(Sys.info()["sysname"] == "Darwin"){
  source("~/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R", chdir = TRUE) #call cource file for functions
  source("~/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R", chdir = TRUE) #call cource file for functions
} else if(Sys.info()["sysname"] == "Windows"){
  source('C:/Users/Blaire/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R', chdir = TRUE) #call cource file for functions
}





#occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Artiodactyla,Perissodactyla&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
#occs <- occs[occs$cc %in% c("US", "CA", "MX"), ]
#occs <- occs[!occs$order %in% c("Desmostylia", "Perissodactyla"), ]

occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]

# occs <- occs[occs$cc %in% c("US", "CA", "MX"), ]

#occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)

ranges <- getTaxonRangesFromOccs(occs=occs, random=TRUE)
rownames(ranges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(ranges))

int_length <- 1
intervals <- makeIntervals(1, 55.5, int_length)
intList <- listifyMatrixByRow(intervals)

####################################################################################################################################
print("Building measurement matrix...")

thisMat <- getSingleSpeciesMatrix()
thisMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = thisMat$species)

#regCatMat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv", na.strings=c("", "NA"), stringsAsFactors = TRUE)
regCatMat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv", na.strings=c("", "NA"), stringsAsFactors = TRUE)
thisMat <- appendMissingPaleoDBSpecies(thisMat, tax.vec=rownames(ranges))		# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files

print("Building body mass estimates...")
# thisMat <- setMeasureReg(occs= occs, regMat= regCatMat)
thisMat <- appendRegressionCategories(thisMat, regCatMat)
# nrow(thisMat)
# thisMat$species <- rownames(thisMat)
thisMat[!is.na(thisMat$family) & thisMat$family=="Entelodontidae", c("P2_L", "P3_L", "p2_w", "m2_w", "m3_w")] <- NA	### entelodont tooth widths were generating >5 ton body masses, so dropped here.
thisMat$bodyMass <- getBodyMassVectorFromThisMatAllMeasures(thisMat = thisMat, linked.files=TRUE)
thisMat$bodyMass <- fillMissingBodyMasses(thisMat)			# this fills taxa missing their body mass with the average body mass of its cogeners
# thisMat[!sapply(thisMat, is.finite)] <- NA
thisMat <- thisMat[is.finite(thisMat$bodyMass),]
# 2.2*(10^(thisMat$bodyMass[!is.na(thisMat$family) & thisMat$family=="Entelodontidae"]))
rownames(thisMat) <- thisMat$species
# length(rownames(thisMat))
####################################################################################################################################

# focal.order <- "Artiodactyla"
# focal.order <- "Perissodactyla"
focal.order <- c("Artiodactyla", "Perissodactyla")
bigList <- unique(occs[occs$accepted_rank =="species", c("order","family", "genus", "accepted_name")])
bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFam <- sort(unique(bigList$family[bigList$order %in% focal.order]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
####################################################################################################################################

load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=TRUE##------ Mon Jan 14 17:48:35 2019 ------##.Rdata")

repIntUng <- getUngulateOnly(repIntSp = repIntSp, bigList = bigList, shortFam = shortFam, thisMat = thisMat)

#########look into range through?
#test <- makeRangeThroughOneRep(repIntUng[[1]])
#tax.vec <- sort(unique(unlist(this.rep)))
#tax.ranges <- cbind(tax.vec, t(sapply(tax.vec, function(taxon) rev(range(which(sapply(this.rep, function(y, taxon) taxon %in% y, taxon)))))))
#for (i in seq_len(nrow(tax.ranges))) {
#  for (k in tax.ranges[i, 2]:tax.ranges[i, 3]) this.rep[[k]] <- sort(unique(c(this.rep[[k]], tax.ranges[i, 1])))
#}
#this.rep


#######################################################################################################################
#my plot for comparing extinction and origination across breaks (extinct from 1 and orig form 2, regime 1| regime 2)
OrigExtList <- list()
for(xx in seq(1, length(repIntUng),1)){
  OrigExtList[[xx]] <- findOrigExt_Breaks(repIntSpSingle = repIntUng[[xx]],optList = optList_bm_median,
                                           measureMat = measureMat, bigList = bigList, shortFam = shortFam,
                                           breaks = c(50,46,40,26,16), useBreaks = TRUE, startDate = 56,
                                            endDate = 3, plot.check = FALSE)
}

quartz()
par(mfrow=c(plot.row,2), mar=c(2,4,1,0.5))
plotOrigExt(OrigExtList = OrigExtList,plot.type = "lines") #points doesn't work
#######################################################################################################################

#what I discussed with Jon for NAPC abstract

#get origination and compare with surviviors and extinct species

#get extinction and compare with new species and survivors 

#make code that will do K-S test on each interval

#findOrigExt_Breaks(repIntSpSingle = repIntUng[[xx]],optList = optList_bm_median,
#                   measureMat = thisMat, bigList = bigList, shortFam = shortFam,
#                   breaks = rev(seq(3,55,1)), useBreaks = TRUE, startDate = 56,
#                   endDate = 3, plot.check = FALSE) #this is now broken need to figure out what happened 
																										#(likely mistook for intervals functiona nd editted by msitake)

#test <- findOrigExt_Interval(repIntSpSingle = repIntUng[[1]], measureMat = thisMat,
#                     bigList=bigList,shortFam = shortFam,
#                     intervals = intervals,startDate = 56, endDate = 3)

repOrigUng <- list()
repExtUng <- list()

# get origination 54.5 to 1.5 intervals
for(gg in seq(1,length(repIntUng[[1]])-1,1)) repOrigUng[[gg]] <- repIntUng[[1]][[gg]][!repIntUng[[1]][[gg]] %in% repIntUng[[1]][[gg+1]]]
# get extinction 55.5 to 2.5 intervals
for(gg in seq(1,length(repIntUng[[1]])-1,1)) repExtUng[[gg]] <- repIntUng[[1]][[gg+1]][!repIntUng[[1]][[gg+1]] %in% repIntUng[[1]][[gg]]]

#get body masses
repOrigBmUng <- getBm4List(repOrigUng,thisMat = thisMat)
repExtBmUng <- getBm4List(repExtUng,thisMat = thisMat)
repBmUng <- getBm4List(repIntUng, thisMat = thisMat)

#####check ks between regimes
breaks <- c(50,46,40,26,16)
repIntUng 
taxList <- getBreakTaxList(repIntSpSingle = repIntUng[[1]], breaks = breaks, startDate = 56, endDate = 3)
names(taxList) <- c(">50","50to46","46to40","40to26","26to16","<16")

regimeOrigUng <- list()
regimeExtUng <- list()

# get origination 
for(gg in seq(1,length(taxList)-1,1)) regimeOrigUng[[gg]] <- taxList[[gg]][!taxList[[gg]] %in% taxList[[gg+1]]]
# get extinction 
for(gg in seq(1,length(taxList)-1,1)) regimeExtUng[[gg]] <- taxList[[gg+1]][!taxList[[1]][[gg+1]] %in% taxList[[1]][[gg]]]
regimeBmOrigUng <- getBm4List(regimeOrigUng,thisMat = thisMat)
regimeBmExtUng <- getBm4List(regimeExtUng,thisMat = thisMat)
regimeBmUng <- getBm4List(taxList,thisMat = thisMat)

getKS4Intervals(regimeBmOrigUng, taxList, result.out = "SigInt", out.putType = "Breaks")
getKS4Intervals(regimeBmOrigUng, taxList, result.out = "AllInt", out.putType = "Breaks")

getKS4Intervals(regimeBmExtUng, taxList, result.out = "SigInt", out.putType = "Breaks")
getKS4Intervals(regimeBmExtUng, taxList, result.out = "AllInt", out.putType = "Breaks")

####compare without orig/ext
sigInt <- getKS4Intervals(repBmUng, repBmUng, result.out = "SigInt", out.putType = "Intervals")
allInt <- getKS4Intervals(repBmUng, repBmUng, result.out = "AllInt", out.putType = "Intervals")

####compare to the previous interval
sigInt.Orig <- getKS4IntervalsOrig(repOrigBmUng, repBmUng, result.out = "SigInt", out.putType = "Intervals")
allInt.Orig <- getKS4IntervalsOrig(repOrigBmUng, repBmUng, result.out = "AllInt", out.putType = "Intervals")

sigInt.Ext <- getKS4IntervalsExt(repExtBmUng, repBmUng, result.out = "SigInt", out.putType = "Intervals")
allInt.Ext <- getKS4IntervalsExt(repExtBmUng, repBmUng, result.out = "AllInt", out.putType = "Intervals")

###make function to get ks across all 1000 iterations of repIntSp
#function needs to aggregate position of significant difference btwn intervals using the ks test

repOrigUng <- list()
repExtUng <- list()
OrigExt.Master <- list()
Orig.numeric.is0 <- list()
Ext.numeric.is0 <- list()
All.numeric.is0 <- list()
repBmUng <- getBm4List(repIntUng, thisMat = thisMat)

for(xx in seq(1, length(repIntUng),1))
{
	# get origination 54.5 to 1.5 intervals
	for(gg in seq(1,length(repIntUng[[xx]])-1,1)) repOrigUng[[gg]] <- repIntUng[[xx]][[gg]][!repIntUng[[xx]][[gg]] %in% repIntUng[[xx]][[gg+1]]]
	# get extinction 55.5 to 2.5 intervals
	for(gg in seq(1,length(repIntUng[[xx]])-1,1)) repExtUng[[gg]] <- repIntUng[[xx]][[gg+1]][!repIntUng[[xx]][[gg+1]] %in% repIntUng[[xx]][[gg]]]
	print(xx)
	#get body masses
	repOrigBmUng <- getBm4List(repOrigUng,thisMat = thisMat)
	repExtBmUng <- getBm4List(repExtUng,thisMat = thisMat)
	
	#get List with reps that had to have numeric(0) set to 0
	Orig.numeric.is0[[xx]] <- check4num0(repIntSp = repOrigBmUng, repNum=xx,out.type ="num0List")
	Ext.numeric.is0[[xx]] <- check4num0(repIntSp = repExtBmUng, repNum=xx,out.type ="num0List")
	All.numeric.is0[[xx]] <- check4num0(repIntSp = repBmUng[[xx]], repNum=xx,out.type ="num0List")
	
	#get repList with 0 enters for numeric(0)
	repOrigBmUng <- check4num0(repIntSp = repOrigBmUng, repNum=xx,out.type ="RepList")
	repExtBmUng <- check4num0(repIntSp = repExtBmUng, repNum=xx,out.type ="RepList")
	repBmUng[[xx]] <- check4num0(repIntSp = repBmUng[[xx]], repNum=xx,out.type ="RepList") 
	
	sigInt.Orig <- getKS4IntervalsOrig(repOrigBmUng, repBmUng[[xx]], result.out = "SigInt", out.putType = "Intervals")
	sigInt.Ext <- getKS4IntervalsExt(repExtBmUng, repBmUng[[xx]], result.out = "SigInt", out.putType = "Intervals")

	Orig.int <- which(!is.na(sigInt.Orig[1,]))+0.5
	Ext.int <- which(!is.na(sigInt.Ext[1,]))+1.5

	OrigExt.int <- list(Origination = Orig.int, Extinction = Ext.int)
	
	OrigExt.Master[[xx]] <- OrigExt.int
}

quartz()
par(mfrow = c(2,1))
OrigAll <- lapply(OrigExt.Master, function(x) x$Origination)
hist(unlist(OrigAll), breaks =seq(1, 56,1))

ExtAll <- lapply(OrigExt.Master, function(x) x$Extinction)
hist(unlist(ExtAll), breaks =seq(1, 56,1))

#data missing from subsampling process in generation of repIntSp
All.num0.Mat <- getListNum0Ints(numList = All.numeric.is0)

#data missing from the intervals being identical
Orig.num0.Mat <- getListNum0Ints(numList = Orig.numeric.is0)
Ext.num0.Mat <- getListNum0Ints(numList = Ext.numeric.is0)

colnames(All.num0.Mat) <- colnames(Orig.num0.Mat) <- colnames(Ext.num0.Mat) <- c("rep","Interval")

#function to find if rangethrough is working
repIntTest <- repIntUng[[1]]
ints <- 50:55

repIntCheck <- checkRangeThrough(repIntTest = repIntUng[[1]], ints = 1:55)

#use Jon's function
repIntTest <- makeRangeThroughOneRep(repIntTest)

#use my own range-through function
reptest <- rangeThrough_alt(repIntCheck)

repIntRangeCheck <- checkRangeThrough(repIntRangeThrough,1:55)
test2 <- checkRangeThrough(reptest,1:55)

repIntCheck[200,]
repIntRangeCheck[200,]
test2[200,]

repIntUng.rangethrough <- list()
repIntUng.rangeMat <- lapply(repIntUng, function(x) checkRangeThrough(x, ints = 1:55))
repIntUng.fixed <- lapply(repIntUng.rangeMat,rangeThrough_alt)

repIntUng.original <- repIntUng

repIntUng <- repIntUng.fixed

