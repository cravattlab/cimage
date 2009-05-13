library(xcms)
source("/home/chuwang/svnrepos/R/msisotope.R")


input.file.suffix <- ".mzXML"
output.file.suffix <- ".ps"
output.path <- "output"
dir.create(output.path)

## probe's mass added to peptide
probe.mass <- 464.28596
## 10 ppm to extract chromatographic peaks
mz.ppm.cut <- 0.000010
# From Eranthie's isotopically labeled probe
pair.mass.delta <- 6.01381
# nature mass difference between C12 and C13
isotope.mass.unit <- 1.0033548
mass.shift <- round( pair.mass.delta/isotope.mass.unit )
# mass of a proton
Hplus.mass <- 1.0072765

## file name from input args
args <- commandArgs(trailingOnly=T)
## the table with protein names
original.table <- read.table(args[1],sep="\t", header=T)
## the table with mass and scan number from DTASelect
cross.table <- read.table(args[2], header=T)
cross.table[,"mass"] <- cross.table[,"mass"] + probe.mass
## file name tags
cross.vec <- dimnames(cross.table)[[2]][-c(1,2)]

## find all matched input files in current directory
input.path <- getwd()
#input.files <- list.files( path=input.path, pattern=paste(input.file.suffix,"$",sep="") )
mzXML.names <- list.files(path="./",pattern="mzXML$")
mzXML.files <- as.list( mzXML.names )
names(mzXML.files) <- mzXML.names
for (name in mzXML.names) {
  cat(paste(name,"\n",sep=""))
  mzXML.files[[name]] <- xcmsRaw( paste("./",name,sep=""), profstep=0, includeMSn=T)
}


# retention time window in secs
rt.window <- 10
rt.window.width <- rt.window * 60
out.filename <- paste(output.path, "/", "output_rt_",as.character(rt.window),".ps", sep="")
## make plot for each isotopic cluster
##png(out.filename)
postscript( out.filename, horizontal=F)
par(mfrow=c(length(cross.vec),1), oma=c(0,0,5,0), las=0)

for ( i in 1:dim(cross.table)[1] ) {
  peptide <- cross.table[i,1]
  mono.mass  <- cross.table[i,"mass"]
  mass <- mono.mass + (which.max(isotope.dist( averagine.count(mono.mass) )) - 1)*isotope.mass.unit
  raw.file.charge <- raw.scan.num <- cross.vec
  for( j in 1:length(cross.vec) ) {
    if (original.table[i,j+3] ==0 ) {
      raw.file.charge[j] <- "none.none"
      raw.scan.num[j] <- NA
    } else {
      tmp.vec <- unlist( strsplit(as.character(cross.table[i,cross.vec[j]]), "\\.") )
      tmp2.vec <- unlist( strsplit(tmp.vec[1],"_" ) )
      raw.file.charge[j] <- paste(tmp2.vec[length(tmp2.vec)], tmp.vec[length(tmp.vec)],sep=".")
      raw.scan.num[j] <- tmp.vec[2]
    }
  }
  count <- 0
  for ( f in levels(as.factor(raw.file.charge))  ) {
    if (f == "none.none") next
      count <- count+1
    tmp.vec <- unlist( strsplit(f,"\\.") )
    raw.index <- tmp.vec[1]
    charge <- as.integer(tmp.vec[2])
    ## mz
    mz.light <- mass/charge + Hplus.mass
    mz.heavy <- (mass+pair.mass.delta)/charge + Hplus.mass
    ## scan number
    ms1.scan.num <- exist.index <- which( raw.file.charge == f )
    for ( k in 1:length(exist.index) ) {
      kk <- exist.index[k]
      raw.file <- paste( cross.vec[kk], "_", raw.index,".mzXML",sep="")
      xfile <- mzXML.files[[raw.file]]
      ms1.scan.num[k] <- which(xfile@acquisitionNum > as.integer(raw.scan.num[kk]))[1]-1
    }

    for ( j in 1:length(cross.vec) ) {
      if ( j %in% exist.index ) {
        tag = "*"
      } else {
        tag = ""
      }
      raw.file <- paste( cross.vec[j], "_", raw.index,".mzXML",sep="")
      xfile <- mzXML.files[[raw.file]]
      ##chromatogram bottom
      raw.ECI.light <- rawEIC(xfile, c(mz.light*(1-mz.ppm.cut), mz.light*(1+mz.ppm.cut)) )
      raw.ECI.heavy <- rawEIC(xfile, c(mz.heavy*(1-mz.ppm.cut), mz.heavy*(1+mz.ppm.cut)) )
      rt.min <- xfile@scantime[min(ms1.scan.num)]-rt.window.width
      if (is.na(rt.min) ) {
        rt.min <- xfile@scantime[length(xfile@scantime)] - 2*rt.window.width
      } else {
        rt.min <- max(rt.min, xfile@scantime[1] )
      }
      rt.max <- xfile@scantime[max(ms1.scan.num)]+rt.window.width
      if (is.na(rt.max) ) {
        rt.max <- xfile@scantime[length(xfile@scantime)]
      } else {
        rt.max <- min(rt.max,xfile@scantime[length(xfile@scantime)])
      }
      if ( (rt.max - rt.min) < 2*rt.window.width ) rt.min <- rt.max-2*rt.window.width
      xlimit <-c(which(xfile@scantime>rt.min)[1]-1, which(xfile@scantime>rt.max)[1] )
      if (is.na(xlimit[2]) ) xlimit[2] <- length(xfile@scantime)
      ylimit <- range(c(raw.ECI.light[[2]][xlimit[1]:xlimit[2]], raw.ECI.heavy[[2]][xlimit[1]:xlimit[2]]))
      ylimit[1] <- 0.0
      ##xlimit <- range(c(raw.ECI.light[[1]], raw.ECI.heavy[[1]]))
      ##ylimit <- range(c(raw.ECI.light[[2]], raw.ECI.heavy[[2]]))
      xlimit <- c(rt.min,rt.max)/60
      raw.ECI.light.rt <- xfile@scantime[ raw.ECI.light[[1]] ] / 60
      raw.ECI.heavy.rt <- xfile@scantime[ raw.ECI.heavy[[1]] ] / 60
      plot(raw.ECI.light.rt, raw.ECI.light[[2]], type="l", col="red",xlab="Retention Time(min)",
           ylab="intensity",
           main=paste(tag, raw.file, "; Census ratio:", as.character(original.table[i,j+3]),
             "; Raw Scan:", as.character(raw.scan.num[j])),
           xlim=xlimit,ylim=ylimit)
      lines(raw.ECI.heavy.rt, raw.ECI.heavy[[2]], col='blue', xlim=xlimit, ylim=ylimit)
                                        #rt.tick <- round(seq(from=xlimit[1],to=xlimit[2], length=5))
                                        #axis(3, at=rt.tick, labels=as.character(format((xfile@scantime[rt.tick])/60,digits=7) ) )
    }
    tt <- paste("Entry ", as.character(i), "-  Charge: ", as.character(charge),
                " - M/Z: ", as.character(format(mz.light, digits=7)),
                "and", as.character(format(mz.heavy,digits=7)))
    mtext(tt, line=3.5, outer=T)
    mtext(paste(original.table[i,"Peptide_Sequence"],"; Mono.mass: ", as.character(mono.mass), sep=""),
         cex=0.8, line=2, out=T)
    mtext(substr(original.table[i,"Description"], 1,100), line=1, cex=0.5,out=T)
  }
} ## each entry
dev.off()
