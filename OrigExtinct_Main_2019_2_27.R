source("~/Dropbox/ungulate_RA/RCode/OrigExtinctScript_Draft_2018_1_9.R")

int_length <- 2
intervals <- makeIntervals(1, 55.5, int_length)
intList <- listifyMatrixByRow(intervals)

####################################################################################################################################
#start KS analysis

# #load 1Ma repIntSp
# load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/RepIntSp_SampleStandardized=TRUE##------ Mon Jan 14 14:30:05 2019 ------##.Rdata")

# #load 2Ma repIntSp
load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=TRUE/RepIntSp_SampleStandardized=TRUE##------ Fri Mar 29 21:20:21 2019 ------##.Rdata")

repIntUng <- getUngulateOnly(repIntSp = repIntSp, bigList = bigList, shortFam = shortFam, thisMat = thisMat)

#check and perform range-through
repIntUng.rangeMat <- lapply(repIntUng, function(x) checkRangeThrough(x, ints = seq(1,57,2)))
repIntUng.fixed <- lapply(repIntUng.rangeMat,rangeThrough_alt)

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
All_orig.numeric.is0 <- list()
All_ext.numeric.is0 <- list()
repBmUng <- getBm4List(repIntUng, thisMat = thisMat)

for(xx in seq(1, length(repIntUng),1))
{
	# get origination 54.5 to 1.5 intervals
	for(gg in seq(1,length(repIntUng[[xx]])-1,1)) repOrigUng[[gg]] <- repIntUng[[xx]][[gg]][!repIntUng[[xx]][[gg]] %in% repIntUng[[xx]][[gg+1]]]
	# get extinction 55.5 to 2.5 intervals
	for(gg in seq(1,length(repIntUng[[xx]])-1,1)) repExtUng[[gg]] <- repIntUng[[xx]][[gg+1]][!repIntUng[[xx]][[gg+1]] %in% repIntUng[[xx]][[gg]]]
	#print(xx)
	
	##Need to add code to remove species that undergo orignation and extinction within each interval 
	#to compare with the survivors of ther previous/current interval, respectively.
	repIntUng_OrigRem <- list()
	repIntUng_ExtRem <- list()
	for(gg in seq(1,length(repIntUng[[xx]])-1,1))	repIntUng_OrigRem[[gg]] <- repIntUng[[xx]][[gg]][!repIntUng[[xx]][[gg]] %in% repOrigUng[[gg]]]
	
	for(gg in seq(1,length(repIntUng[[xx]])-1,1))	repIntUng_ExtRem[[gg]] <- repIntUng[[xx]][[gg+1]][!repIntUng[[xx]][[gg+1]] %in% repExtUng[[gg]]]
	
	#get vector of interval centers (e.g. 55 to 56 = 55.5, 1 to 3 = 2, etc.)
	int.label <- gsub(pattern = " Ma", replacement = "", x = names(repIntUng[[xx]]))
	
	#get body masses
	repIntBmUng_OrigRem <- getBm4List(repIntUng_OrigRem,thisMat = thisMat)
	repIntBmUng_ExtRem <- getBm4List(repIntUng_ExtRem,thisMat = thisMat)
	repOrigBmUng <- getBm4List(repOrigUng,thisMat = thisMat)
	repExtBmUng <- getBm4List(repExtUng,thisMat = thisMat)
	
	#get List with reps that had to have numeric(0) set to 0
	Orig.numeric.is0[[xx]] <- check4num0(repIntSp = repOrigBmUng, repNum=xx,out.type ="num0List")
	Ext.numeric.is0[[xx]] <- check4num0(repIntSp = repExtBmUng, repNum=xx,out.type ="num0List")
	All_orig.numeric.is0[[xx]] <- check4num0(repIntSp = repIntBmUng_OrigRem, repNum=xx,out.type ="num0List")
	All_ext.numeric.is0[[xx]] <- check4num0(repIntSp = repIntBmUng_ExtRem, repNum=xx,out.type ="num0List")
	
	#get repList with 0 enters for numeric(0)
	repOrigBmUng <- check4num0(repIntSp = repOrigBmUng, repNum=xx,out.type ="RepList")
	repExtBmUng <- check4num0(repIntSp = repExtBmUng, repNum=xx,out.type ="RepList")
	#repBmUng[[xx]] <- check4num0(repIntSp = repBmUng[[xx]], repNum=xx,out.type ="RepList") 
	repIntBmUng_OrigRem <- check4num0(repIntSp = repIntBmUng_OrigRem, repNum=xx,out.type ="RepList")
	repIntBmUng_ExtRem <- check4num0(repIntSp = repIntBmUng_ExtRem, repNum=xx,out.type ="RepList")
	
	sigInt.Orig <- getKS4IntervalsOrig(repOrigBmUng, repIntBmUng_OrigRem, result.out = "SigInt", out.putType = "Intervals")
	sigInt.Ext <- getKS4IntervalsExt(repExtBmUng, repIntBmUng_ExtRem, result.out = "SigInt", out.putType = "Intervals")

	Orig.int <- which(!is.na(sigInt.Orig[1,])) #need this to dynamically work with dif number and sizes of intervals
	Ext.int <- which(!is.na(sigInt.Ext[1,]))

	#get list and stagger orig and extinction
	Orig.int <- as.numeric(int.label[Orig.int])
	Ext.int <- as.numeric(int.label[Ext.int])+int_length
	
	OrigExt.int <- list(Origination = Orig.int, Extinction = Ext.int)
	
	OrigExt.Master[[xx]] <- OrigExt.int
}

quartz()
par(mfrow = c(2,1))
OrigAll <- lapply(OrigExt.Master, function(x) x$Origination)
hist(unlist(OrigAll), breaks = seq(1, 57, 2), ylim = c(0,1000), xlim = c(57,1), col = "dodgerblue", main = "Origination")

ExtAll <- lapply(OrigExt.Master, function(x) x$Extinction)
hist(unlist(ExtAll), breaks =seq(1, 57, 2), ylim = c(0,1000), xlim = c(57,1), col = "firebrick4", main = "Extinction")

#data missing from subsampling process in generation of repIntSp
All.num0.Mat <- getListNum0Ints(numList = All.numeric.is0)

#data missing from the intervals being identical
Orig.num0.Mat <- getListNum0Ints(numList = Orig.numeric.is0)
Ext.num0.Mat <- getListNum0Ints(numList = Ext.numeric.is0)

colnames(All.num0.Mat) <- colnames(Orig.num0.Mat) <- colnames(Ext.num0.Mat) <- c("rep","Interval")

#Need to now go into results and figure out how taxon and body mass distributions are changing at those Orig/ext shifts
#bm 
##make histograms for orig and ext like for net change? (maybe have them be side by side if neear same break)

#taxon
##break down into family level


####################################################################################
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

