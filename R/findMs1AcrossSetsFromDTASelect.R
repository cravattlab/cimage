suppressPackageStartupMessages(library(xcms))
cimage.path <- Sys.getenv("CIMAGE_PATH")
source(paste(cimage.path,"/R/inputparams.R",sep=""))
source(paste(cimage.path,"/R/msisotope.R",sep=""))

## file name from input args
args <- commandArgs(trailingOnly=T)
param.file <- args[1]
args <- args[-1]

## read parameters from a file
params <- read.input.params(param.file)


## initialize atom mass table
atom.mass.vec <- init.atom.mass()
## initialize chemical composition table
light.chem.table <- read.chem.table(params[["light.chem.table"]])
heavy.chem.table <- read.chem.table(params[["heavy.chem.table"]])
## initialize amino acid mass table
light.aa.mass <- init.aa.mass(atom.mass.vec, light.chem.table)
heavy.aa.mass <- init.aa.mass(atom.mass.vec, heavy.chem.table)

## ppm to extract chromatographic peaks
mz.ppm.cut <- as.numeric(params[["ppm.tolerance"]]) * 1E-6
## N15 enrichment ratio ##
N15.enrichment <- as.numeric(params[["N15.enrichment"]])
## nature mass difference between C12 and C13
isotope.mass.unit <- atom.mass.vec["C13"] - atom.mass.vec["C"]
## natural mass difference between N14 and N15
isotope.mass.unit.N15 <- atom.mass.vec["N15"] - atom.mass.vec["N"]
# mass of a proton
Hplus.mass <- atom.mass.vec["Hplus"]

## output folder
output.path <- params[["output.path"]]
dir.create(output.path)
## the table with protein names
ipi.name.table <- read.table("ipi_name.table",sep="\t",header=T,comment.char="")
## the table with mass and scan number from DTASelect
cross.table <- read.table("cross_scan.table", header=T, check.names=F,comment.char="")
#cross.table[,"mass"] <- cross.table[,"mass"] + probe.mass
split.table <- matrix(unlist(strsplit(as.character(cross.table[,"key"]),":")), byrow=T,ncol=4)
dimnames(split.table)[[2]] <- c("ipi","peptide","charge","segment")
cross.table <- cbind(cross.table, split.table)
uniq.ipi.peptides <- as.factor(paste(cross.table[,"ipi"], cross.table[,"peptide"],sep=":"))
entry.levels <- levels( uniq.ipi.peptides )
## all_scan.table
all.scan.table <- read.table("all_scan.table", header=T, as.is=T,comment.char="",colClasses=c("character","character","numeric","character"))
## file name tags
cross.vec <- as.character(args)
ncross <- length(cross.vec)
# handle switching of Heavy to light ratio, by default it is light vs heavy
HL.ratios <- rep(FALSE,ncross)
j <- 0
for ( arg in as.character(args) ) {
  j <- j+1
  cross.vec[j] <- sub("_HL$","",arg)
  if ( cross.vec[j] != arg ) {
    HL.ratios[j] <- TRUE
  }
}
## find all matched input files in current directory
##if(FALSE) {
input.path <- getwd()
mzXML.names <- list.files(path="../",pattern="mzXML$")
mzXML.files <- as.list( mzXML.names )
names(mzXML.files) <- mzXML.names
for (name in mzXML.names) {
  cat(paste(name,"\n",sep=""))
  mzXML.files[[name]] <- xcmsRaw( paste("../",name,sep=""), profstep=0, includeMSn=T)
}
##}
## more parameters from input files
## retention time window for alignment across multiple samples
rt.window <- as.numeric(params[["rt.window"]])
rt.window.width <- rt.window * 60
## local retention time window for noise line calculation
local.rt.window <- as.numeric(params[["local.rt.window"]])
local.rt.window.width <- local.rt.window * 60
## signal/noise ratio for peak picking
sn <- as.numeric(params[["sn"]])
### range cutoff for calculated ratios###
ratio.range <- as.numeric(params[["ratio.range"]])
### isotope envelope score cutoff ###
env.score.cutoff <- as.numeric(params[["env.score.cutoff"]])
### coelution profile r2 cutoff ###
r2.cutoff <- as.numeric(params[["r2.cutoff"]])
### minimum peak width in numbers of time points###
minimum.peak.points <- as.numeric(params[["minimum.peak.points"]])
### choose peak pairs with MS2 data only ###
peaks.with.ms2.only <- as.logical(params[["peaks.with.ms2.only"]])
### singleton case ratio ###
singleton.ratio <- as.numeric(params[["singleton.ratio"]])
## column names for calculated ratios
integrated.area.ratio <- paste("IR",cross.vec,sep=".")
linear.regression.ratio <- paste("NP",cross.vec,sep=".")
linear.regression.R2 <- paste("R2",cross.vec,sep=".")
light.integrated.area <- paste("INT",cross.vec,sep=".")
column.names <- c("index","ipi", "description", "symbol", "sequence", "mass", "charge", "segment",
                  integrated.area.ratio, light.integrated.area, linear.regression.ratio, linear.regression.R2, "entry", "link" )
out.df <- matrix(nrow=0, ncol=length(column.names))
colnames(out.df) <- column.names

## output name
out.filename.base <- paste("output_rt_",as.character(rt.window),"_sn_",
                      as.character(sn),sep="")
out.filename <- paste(output.path,"/",out.filename.base,".pdf",sep="")
## output layout
pdf( out.filename, height=4*ncross, width=11, paper="special")
layout.vec <- row.layout.vec <- c(1,1,2,1,1,3)
if ( ncross > 1 ) {
  for (i in 1:(ncross-1)) {
    layout.vec <- c(layout.vec,(row.layout.vec+i*3))
  }
}
layout.matrix <- matrix(layout.vec,byrow=T,ncol=3)
layout(layout.matrix)
par(oma=c(0,0,5,0), las=0)

##for ( i in 2304:2304) {
for ( i in 1:dim(cross.table)[1] ) {
  key <- cross.table[i,"key"]
  tmp.vec <- unlist( strsplit(as.character(key),":") )
  ipi <- tmp.vec[1]
  peptide <- tmp.vec[2]
  charge <- as.integer(tmp.vec[3])
  segment <- tmp.vec[4]
  entry.tag <- paste( ipi, peptide, sep=":")
  entry.index <- which( entry.tag == entry.levels )
  description <- as.character(ipi.name.table[ipi,"name"])
  symbol<- strsplit(description, " ")[[1]][1]
  ## momo.mass (light and heavy) and mass of most abundant isotopes (light and heavy)
  mono.mass <- calc.peptide.mass( peptide, light.aa.mass)
  #predicted.dist <- isotope.dist( averagine.count(mono.mass) )
  elements.count <- calc.num.elements(peptide, light.chem.table)
  #number of heavy atoms already in the light table, to calculate envelope distribution properly
  n.heavyatoms <- sum(elements.count[c("N15","H2","C13")])
  predicted.dist <- isotope.dist(elements.count)
  if (n.heavyatoms>0) {predicted.dist <- predicted.dist[-c(1:n.heavyatoms)]}
  i.max <- which.max(predicted.dist)
  mass <- mono.mass + (i.max - 1)*isotope.mass.unit
  mono.mass.heavy <- calc.peptide.mass( peptide, heavy.aa.mass)
  mass.heavy <- mono.mass.heavy + mass - mono.mass
  elements.count.heavy <- calc.num.elements(peptide, heavy.chem.table)
  predicted.dist.heavy <- isotope.dist(elements.count.heavy,N15.enrichment)
  if (n.heavyatoms>0){predicted.dist.heavy <- predicted.dist.heavy[-c(1:n.heavyatoms)]}
  ## mass delta between light and heavy
  mass.shift <- sum((elements.count.heavy-elements.count)[c("N15","H2","C13")])
  correction.factor <- predicted.dist[i.max]/predicted.dist.heavy[i.max+mass.shift]
  ## mz
  mono.mz <- mono.mass/charge + Hplus.mass
  mz.light <- mass/charge + Hplus.mass
  mono.mz.heavy <- mono.mass.heavy/charge + Hplus.mass
  mz.heavy <-  mass.heavy/charge + Hplus.mass
  ## scan number
  raw.scan.num <- cross.table[i,cross.vec]
  ms1.scan.rt <- ms1.scan.num <- exist.index <- which( raw.scan.num > 0 )
  for ( k in 1:length(exist.index) ) {
    kk <- exist.index[k]
    raw.file <- paste( cross.vec[kk], "_", segment,".mzXML",sep="")
    xfile <- mzXML.files[[raw.file]]
    ms1.scan.num[k] <- which(xfile@acquisitionNum > as.integer(raw.scan.num[kk]))[1]-1
    if (is.na(ms1.scan.num[k])) {
      ms1.scan.num[k] <- length(xfile@acquisitionNum)
    }
    ms1.scan.rt[k] <- xfile@scantime[ms1.scan.num[k]]
  }

  r2.v <- l.ratios <- light.int.v <- i.ratios <- rep(NA,ncross)
  for ( j in 1:ncross ) {
    raw.file <- paste( cross.vec[j], "_", segment,".mzXML",sep="")
    xfile <- mzXML.files[[raw.file]]
    ## tag * and tag rt line
    if ( j %in% exist.index ) {
      tag <- "*"
      tag.ms1.scan.num <- ms1.scan.num[match(j,exist.index)]
      tag.rt <- xfile@scantime[tag.ms1.scan.num]/60
    } else {
      tag <- ""
      tag.ms1.scan.num <- NA
      tag.rt <- NA
    }
    ##chromatogram bottom
    raw.ECI.light <- rawEIC(xfile, c(mz.light*(1-mz.ppm.cut), mz.light*(1+mz.ppm.cut)) )
    raw.ECI.heavy <- rawEIC(xfile, c(mz.heavy*(1-mz.ppm.cut), mz.heavy*(1+mz.ppm.cut)) )
    scan.time.range <- range(xfile@scantime)
    rt.min <- min(ms1.scan.rt)-rt.window.width
    if (rt.min > scan.time.range[2]) {
      rt.min <- scan.time.range[2] - 2*rt.window.width
    } else {
      rt.min <- max(rt.min, scan.time.range[1] )
    }
    rt.max <- max(ms1.scan.rt)+rt.window.width
    if (rt.max < scan.time.range[1] ) {
      rt.max <- scan.time.range[1] + 2*rt.window.width
    } else {
      rt.max <- min(rt.max,scan.time.range[2])
    }
    if ( (rt.max - rt.min) < 2*rt.window.width ) {
      if ( rt.max == scan.time.range[2] ) {
        rt.min <- rt.max - 2*rt.window.width
      } else if (rt.min == scan.time.range[1]) {
        rt.max <- rt.min + 2*rt.window.width
      }
    }
    xlimit <-c(which(xfile@scantime>rt.min)[1]-1, which(xfile@scantime>rt.max)[1] )
    if (is.na(xlimit[2]) ) xlimit[2] <- length(xfile@scantime)
    ylimit <- range(c(raw.ECI.light[[2]][xlimit[1]:xlimit[2]], raw.ECI.heavy[[2]][xlimit[1]:xlimit[2]]))
    ylimit[1] <- 0.0
    ylimit[2] <- ylimit[2]*1.2
    local.xlimit <- xlimit <- c(rt.min,rt.max)/60
    raw.ECI.light.rt <- xfile@scantime[ raw.ECI.light[[1]] ] / 60
    raw.ECI.heavy.rt <- xfile@scantime[ raw.ECI.heavy[[1]] ] / 60

    tt.main <- paste(tag, raw.file, "; Raw Scan:", as.character(raw.scan.num[j]),
                     "; NL:", formatC(ylimit[2], digits=2, format="e"))
    plot(raw.ECI.light.rt, raw.ECI.light[[2]], type="l", col="red",xlab="Retention Time(min)",
         ylab="intensity", main=tt.main, xlim=xlimit,ylim=ylimit)
    lines(raw.ECI.heavy.rt, raw.ECI.heavy[[2]], col='blue', xlim=xlimit, ylim=ylimit)
    k.ms1.rt.v <- k.ms1.scan.v <- numeric(0)
    k.ms1.int.light.v <- k.ms1.int.heavy.v <- 0
    if ( !is.na(tag.rt) ) {
      all.ms2.scan <- as.integer( all.scan.table[(key==all.scan.table[,"key"]
                                                  &cross.vec[j]==all.scan.table[,"run"]),"scan"] )
      all.ms2.HL <- all.scan.table[(key==all.scan.table[,"key"]
                                    &cross.vec[j]==all.scan.table[,"run"]),"HL"]
      for (k in 1:length(all.ms2.scan)) {
        k.ms1.scan <- which(xfile@acquisitionNum > all.ms2.scan[k])[1]-1
        if (is.na(k.ms1.scan)) {
          k.ms1.scan <- length(xfile@acquisitionNum)
        }
        k.ms1.rt <- xfile@scantime[k.ms1.scan]/60
        if (all.ms2.HL[k] == "light") {
          points(k.ms1.rt, raw.ECI.light[[2]][k.ms1.scan], type='p',cex=0.5, pch=1)
          #k.ms1.int.light.v <- c(k.ms1.int.light.v, raw.ECI.light[[2]][k.ms1.scan])
        } else if (all.ms2.HL[k] == "heavy") {
          points(k.ms1.rt, raw.ECI.heavy[[2]][k.ms1.scan], type='p',cex=0.5, pch=1)
          #k.ms1.int.heavy.v <- c(k.ms1.int.heavy.v, raw.ECI.heavy[[2]][k.ms1.scan])
        } else {
          points(k.ms1.rt, max(raw.ECI.light[[2]][k.ms1.scan],raw.ECI.heavy[[2]][k.ms1.scan]),
                 type='p',cex=0.5, pch=1,col="black")
        }
        k.ms1.rt.v <- c(k.ms1.rt.v,k.ms1.rt)
        k.ms1.scan.v <- c(k.ms1.scan.v,k.ms1.scan)
      }
      ##lines(c(tag.rt,tag.rt),c(0.0, max(raw.ECI.light[[2]],raw.ECI.heavy[[2]])), col="green")
      HL <- all.ms2.HL[k.ms1.scan.v == tag.ms1.scan.num][1]
      if (HL == "light") {
        points(tag.rt, raw.ECI.light[[2]][tag.ms1.scan.num], type='p',pch=8)
      } else if (HL == "heavy") {
        points(tag.rt, raw.ECI.heavy[[2]][tag.ms1.scan.num], type='p',pch=8)
      } else {
        points(tag.rt, max(raw.ECI.light[[2]][tag.ms1.scan.num],raw.ECI.heavy[[2]][tag.ms1.scan.num]), type='p',pch=8)
      }
      ## record MS1 intensity at which MS2 is triggered
      k.ms1.int.light.v <- raw.ECI.light[[2]][tag.ms1.scan.num]
      k.ms1.int.heavy.v <- raw.ECI.heavy[[2]][tag.ms1.scan.num]
      ## guess ratio of integrated peak area
      local.xlimit <- c(max(scan.time.range[1]/60, tag.rt-local.rt.window),
                        min(scan.time.range[2]/60, tag.rt+local.rt.window))
    }
    ## guess ratio of integrated peak area
    peaks <- findPairChromPeaks( raw.ECI.light.rt, raw.ECI.light[[2]], raw.ECI.heavy[[2]],
                                xlimit, local.xlimit, sn )

    noise.light <- peaks[1]
    lines(xlimit,c(noise.light, noise.light), col='red', type='l', lty=2)
    noise.heavy <- peaks[2]
    lines(xlimit,c(noise.heavy, noise.heavy), col='blue', type='l', lty=2)
    peaks <- peaks[-c(1,2)]
    n.peaks <- length(peaks)/2

    best.peak.scan.num <- best.mono.check <- best.r2 <- best.npoints <- best.light.int <- best.ratio <- 0.0
    best.mono.check <- 0.1
    best.xlm <- best.light.yes <- best.heavy.yes <- best.low <- best.high <- c(0)
    best.fixed <- F
    n.light.ms2.peak <- n.heavy.ms2.peak <- n.candidate.peaks <- n.ms2.peaks <- 0
    if (n.peaks>0) {
      for (n in 1:n.peaks) {
        low <- peaks[2*n-1]
        high<- peaks[2*n]
        ### when requested, choose peaks with ms2 events only ###
        if (peaks.with.ms2.only) {
          if (length(k.ms1.rt.v>0) & (sum((k.ms1.rt.v>=low & k.ms1.rt.v<=high))<=0)) next
        }
        yes <- which( raw.ECI.light.rt>=low & raw.ECI.light.rt<=high )
        light.yes <- raw.ECI.light[[2]][yes]
        heavy.yes <- raw.ECI.heavy[[2]][yes]

        peak.scan.num <- raw.ECI.light[[1]][yes][which.max(light.yes)]
        if ( mass.shift > 0 ) {
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass-2)/charge, mz.heavy) )
        } else {
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass-2)/charge, mz.heavy+5) )
        }
        mono.check <- checkChargeAndMonoMass( peak.scan, mono.mass, charge, mz.ppm.cut, predicted.dist)
        ## calculate ratio of integrated peak area
        if (HL.ratios[j]) {
          ratio <- round((sum(heavy.yes)/sum(light.yes))*correction.factor,digits=2)
          peak.scan.num <- raw.ECI.heavy[[1]][yes][which.max(heavy.yes)]
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass.heavy-2)/charge, mz.heavy+5) )
          mono.check.heavy <- checkChargeAndMonoMass( peak.scan, mono.mass.heavy, charge, mz.ppm.cut, predicted.dist.heavy[(mass.shift+1):length(predicted.dist.heavy)])
          mono.check <- max(mono.check, mono.check.heavy)
        } else {
          ratio <- round((sum(light.yes)/sum(heavy.yes))/correction.factor,digits=2)
        }
        if( singleton.ratio > 0  & ratio > singleton.ratio ) next
        lines(c(low,low),ylimit/10, col="green")
        lines(c(high,high),ylimit/10, col="green")
        text(mean(c(low,high)),max(light.yes,heavy.yes)*1.2,
             labels=paste(round(ratio,2),round(mono.check,2),sep="/"))
        ## calculate peak co-elution profile using only points above noise line
        ##yes2 <- light.yes > noise.light & heavy.yes > noise.heavy
        ##light.yes <- light.yes[yes2]
        ##heavy.yes <- heavy.yes[yes2]
        if (ratio > ratio.range[2] | ratio < ratio.range[1]) next
        if (mono.check < env.score.cutoff) next
        npoints <- length(light.yes)
        if (npoints<minimum.peak.points) {
          next
        }
        ## extra information for better filtering
        if (length(k.ms1.rt.v>0) & (sum((k.ms1.rt.v>=low & k.ms1.rt.v<=high))>0)) {
          n.ms2.peaks <- n.ms2.peaks + 1
        }
        x.lm <- lsfit( x=heavy.yes, y=light.yes,intercept=F )
        r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
        if (r2>r2.cutoff) {
          n.candidate.peaks <- n.candidate.peaks + 1
        }
        if ( !is.na(tag.rt) & tag.rt>=low & tag.rt<=high) {
          best.fixed <- T
        }
        if ( best.fixed | (best.mono.check < 0.95 & mono.check >= best.mono.check) |
            ( best.mono.check >=0.95 & mono.check >=0.95 & max(light.yes)>max(best.light.yes) ) ) {
          best.mono.check <- mono.check
          best.npoints <- npoints
          best.r2 <- r2
          best.ratio <- ratio
          best.light.int <- sum(light.yes)
          best.xlm <- round(as.numeric(ls.print(x.lm,print.it=F)$coef.table[[1]][,"Estimate"]),digits=2)
          best.low <- low
          best.high <- high
          best.light.yes <- light.yes
          best.heavy.yes <- heavy.yes
          best.peak.scan.num <- peak.scan.num
        }
        if (best.fixed) break
      }
    }
    if ( best.r2 != 0 ) {
      lines(c(best.low,best.low),ylimit, col="green")
      lines(c(best.high,best.high),ylimit, col="green")
      plot(best.heavy.yes,best.light.yes,
           xlab="intensity.heavy", ylab="intensity.light",
           main=paste("X=",format(best.xlm,digits=4),"; R2=",format(best.r2,digits=3),
             "; Np=", best.npoints, sep=""),
           xlim=c(0, max(best.light.yes,best.heavy.yes)),
           ylim=c(0, max(best.light.yes,best.heavy.yes)))
      abline(0,best.xlm)
      abline(0,1,col="grey")
      i.ratios[j] <- best.ratio
      light.int.v[j] <- best.light.int
      ##l.ratios[j] <- best.xlm
      r2.v[j] <- best.r2
      ## plot raw spectrum
      ##predicted.dist <- predicted.dist[1:20]
      ## upper limit: heavy + 20units ##
      cc <- seq(1,max(which(predicted.dist.heavy>0.01)))
      if (HL.ratios[j]) {
        predicted.dist.merge <- (1/best.ratio)*predicted.dist[cc] + predicted.dist.heavy[cc]
      } else {
        predicted.dist.merge <- (best.ratio)*predicted.dist[cc] + predicted.dist.heavy[cc]
      }

      mz.unit <- isotope.mass.unit/charge
      ##predicted.mz <- mono.mz + mz.unit*(seq(1,mass.shift)-1)
      light.index <- which(predicted.dist>0.01)-1
      light.index <- light.index[which(light.index<=mass.shift)]
      predicted.mz <- mono.mz + mz.unit*light.index
      predicted.dist.local <- predicted.dist.merge[light.index+1]
      #predicted.mz.heavy <- mono.mz.heavy + mz.unit*(seq(1, length(predicted.dist.merge)-mass.shift)-1)
      mz.unit.N15 <- isotope.mass.unit.N15/charge
      heavy.index <- which(predicted.dist.heavy>0.01)
      predicted.dist.heavy.local <- predicted.dist.merge[heavy.index]
      heavy.adjustments <- heavy.index <- heavy.index-mass.shift-1
      heavy.adjustments[which(heavy.index<0)] <- mz.unit.N15
      heavy.adjustments[which(heavy.index>=0)] <- mz.unit
      predicted.mz.heavy <- mono.mz.heavy + heavy.adjustments*heavy.index

      predicted.mz <- c(predicted.mz, predicted.mz.heavy)

      predicted.dist.merge <- c(predicted.dist.local,predicted.dist.heavy.local)
      n.max <- which.max(predicted.dist.merge)
      predicted.dist.merge <- predicted.dist.merge/predicted.dist.merge[n.max]

      mz.max <- predicted.mz[n.max]
      mass.range <- c(mono.mz-2*mz.unit, mz.heavy+8*mz.unit)
      scan.data <- getScan(xfile, best.peak.scan.num, mzrange=mass.range)
      scan.mz <- scan.data[,"mz"]
      scan.int <- scan.data[,"intensity"]
      observed.int <- predicted.mz
      for ( k in 1:length(observed.int) ) {
        this.mz <- predicted.mz[k]
        mz.diff <- abs(scan.mz-this.mz)/this.mz
        if (min(mz.diff) <= mz.ppm.cut ) {
          observed.int[k] <- scan.int[which.min(mz.diff)]
        } else {
          observed.int[k] <- 0.0
        }
      }
      int.max <- observed.int[n.max]
      if ( max(int.max) > 0.0 ) {
        observed.int <- observed.int / int.max
        scan.int <- scan.int / int.max
      } else {
        scan.int <- scan.int/max(scan.int)
      }
      ylimit2 <- c(0,1.1)
      plot(scan.mz, scan.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit2, col="gray")
      par(new=T)
      plot(predicted.mz, observed.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit2, col="black")
      if ( mass.shift >0 ){
        light.n <- seq(1,length(light.index))##seq(1,(mass.shift))
        heavy.n <- seq(length(light.index)+1, length(predicted.mz))##seq((mass.shift+1),2*mass.shift)
      } else {
        light.n <- heavy.n <- seq(1,3)
      }
      par(new=T)
      plot( predicted.mz[light.n], predicted.dist.merge[light.n], type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit2)
      par(new=T)
      plot( predicted.mz[heavy.n], predicted.dist.merge[heavy.n], type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit2)

      points(predicted.mz[light.n],rep(0,length(light.n)), pch=23,col="red",bg="white")
      points(predicted.mz[heavy.n],rep(0,length(heavy.n)), pch=24,col="blue",bg="white")
      points(mz.light,0, pch=24,col="red",bg="red")
      points(mz.heavy,0, pch=24,col="blue",bg="blue")
      title( paste("Scan # ", xfile@acquisitionNum[best.peak.scan.num], " @ ",
                   round(xfile@scantime[best.peak.scan.num]/60,1)," min; NL:",
                   formatC(int.max, digits=2,format="e"), sep = ""))
    } else {
      if (peaks.with.ms2.only & (singleton.ratio>0)) {
        if (HL.ratios[j]) {
          mono.single <- mono.mass.heavy
          raw.ECI.rt.single <- raw.ECI.heavy.rt
          raw.ECI.single <- raw.ECI.heavy
          k.ms1.int.ratio <- max(k.ms1.int.heavy.v)/max(c(0.01,k.ms1.int.light.v))
        } else {
          mono.single <- mono.mass
          raw.ECI.rt.single <- raw.ECI.light.rt
          raw.ECI.single <- raw.ECI.light
          k.ms1.int.ratio <- max(k.ms1.int.light.v)/max(c(0.01,k.ms1.int.heavy.v))
        }
        single.peaks <- findSingleChromPeaks(raw.ECI.rt.single, raw.ECI.single[[2]],xlimit, local.xlimit, sn )
        single.peaks <- single.peaks[-1]
        n.single.peaks <- length(single.peaks)/2
        n.singleton.peaks <- numeric(0)
        if (n.single.peaks>0) {
          for (ns in 1:n.single.peaks) {
            low.single <- single.peaks[2*ns-1]
            high.single <- single.peaks[2*ns]
            k.ms1.rt.v.tmp <- (k.ms1.rt.v >=low.single & k.ms1.rt.v <= high.single)
            if (length(k.ms1.rt.v>0) & (sum(k.ms1.rt.v.tmp)<=0)) next
            if (k.ms1.int.ratio < 3) next
            yes.single <- which( raw.ECI.rt.single>=low.single & raw.ECI.rt.single<=high.single )
            int.yes.single <- raw.ECI.single[[2]][yes.single]
            peak.scan.num <- raw.ECI.single[[1]][yes.single][which.max(int.yes.single)]
            peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.single-2)/charge,
                                                         (mono.single+10)/charge) )
            mono.check.single <- checkChargeAndMonoMass( peak.scan, mono.single, charge, mz.ppm.cut,
                                                        predicted.dist)
            ## did not pass env score filter
            if (mono.check.single < env.score.cutoff) next
            npoints.single <- length(yes.single)
            ## peak is too narrow with very few time points
            if (npoints.single<minimum.peak.points) next
            ##
            lines(c(low.single,low.single),ylimit/2,col="green")
            lines(c(high.single,high.single),ylimit/2,col="green")
            text(mean(c(low.single,high.single)),max(int.yes.single)*1.2, labels=paste(round(singleton.ratio,2),round(mono.check.single),sep="/"))
            i.ratios[j] <- singleton.ratio
            r2.v[j] <- 1.0
            plot(0,0,xlab="",ylab="",main=paste(npoints.single, " MS1s in this peak", sep="") )
            plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
            n.singleton.peaks <- c(n.singleton.peaks,ns)
            break
          }
        }
        if (length(n.singleton.peaks) == 0) {
          plot(0,0,xlab="",ylab="",main=paste("No singleton peak found") )
          plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
        }
      } else {
        ##plot(0,0,xlab="",ylab="",main=paste("R2 value: 0.00") )
        plot(0,0,xlab="",ylab="",main=paste("Singleton detection off") )
        plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
      }
    }
    if (mass.shift==0) {
      i.ratios[j] <- 0.0
      r2.v[j] <- 0.0
    } ## for non-isotopic peptides, assign ratios to zero
    l.ratios[j] <- paste(n.ms2.peaks, n.candidate.peaks,
                         format(max(k.ms1.int.light.v), digits=1, scientific=T),
                         format(noise.light, digits=1, scientific=T),
                         format(max(k.ms1.int.heavy.v), digits=1, scientific=T),
                         format(noise.heavy, digits=1, scientific=T),
                         sep="/")
  } ## each ratio j
  tt <- paste("Entry ", as.character(i), "-  Charge: ", as.character(charge),
              " - M/Z: ", as.character(format(mz.light, digits=7)),
              "and", as.character(format(mz.heavy,digits=7)))
  mtext(tt, line=3.5, outer=T)
  mtext(paste(cross.table[i,"peptide"],"; Mono.mass: ", as.character(mono.mass), "; Mono.mz: ", as.character(round(mono.mz,5)),sep=""),
        cex=0.8, line=2, outer=T)
  mtext(paste(cross.table[i,"ipi"],description),line=0.8, cex=0.8,out=T)
  ## save data in outdf
  lnk.i <- ceiling(i/500)-1
  lnk.j <- (i-1)%%500
  lnk.name <- paste('./PNG/', lnk.i, '/', out.filename.base,'.', lnk.i, '-', lnk.j,'.png',sep='')
  this.df <- c(i, ipi, description, symbol, peptide, round(mass,digits=4), charge, segment,
               i.ratios, light.int.v, l.ratios, r2.v, entry.index,
               paste('=HYPERLINK(\"./PNG/', lnk.i, '/', out.filename.base,'.', lnk.i, '_', lnk.j,'.png\")',sep=''))
  names(this.df) <- column.names
  out.df <- rbind(out.df, this.df)
} ## each entry i
dev.off()

all.table <- out.df
all.table.out <- all.table[F,]
rsq.cutoff <- r2.cutoff

## go from high concentration to low concentration,
## first apply R2 cutoff and sort by IR values
for ( s in seq(ncross, 1) ) {
  colname.R2 <- linear.regression.R2[s]
  colname.IR <- integrated.area.ratio[s]
  rsq.filter <- all.table[,colname.R2] >= rsq.cutoff & !is.na(all.table[,colname.R2])
  table <- all.table[rsq.filter,]
  if (is.vector(table)) {
    table <- data.frame(as.list(table))
  }
  s1 <- as.numeric(table[,colname.IR])
  s2 <- as.numeric(table[,"entry"])
  s3 <- as.numeric(table[,"charge"])
  s4 <- as.numeric(table[,"segment"])
  ii <- order(s1, s2, s3, s4)
  table <- table[ii,]
  all.table.out <- rbind(all.table.out, table)
  all.table <- all.table[!rsq.filter,]
}
all.table.out <- rbind(all.table.out, all.table)
## output the final table
row.names(all.table.out) <- as.character(seq(1:dim(all.table.out)[1]) )
all.table.out[,"index"] <- seq(1:dim(all.table.out)[1])
write.table(all.table.out,file=paste(output.path,"/",out.filename.base,".to_excel.txt",sep=""),
            quote=F, sep="\t", row.names=F,na="0.00")

