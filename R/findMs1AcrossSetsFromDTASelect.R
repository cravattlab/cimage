library(xcms)
source("/home/chuwang/svnrepos/R/msisotope.R")

## probe's mass added to peptide
probe.mass <- c(464.28596)
names(probe.mass) <- c("C")
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
output.path <- "output"
dir.create(output.path)
## the table with protein names
ipi.name.table <- read.table("ipi_name.table",sep="\t",header=T)
## the table with mass and scan number from DTASelect
cross.table <- read.table("cross_scan.table", header=T, check.names=F)
#cross.table[,"mass"] <- cross.table[,"mass"] + probe.mass
split.table <- matrix(unlist(strsplit(as.character(cross.table[,"key"]),":")), byrow=T,ncol=4)
dimnames(split.table)[[2]] <- c("ipi","peptide","charge","segment")
cross.table <- cbind(cross.table, split.table)
uniq.ipi.peptides <- as.factor(paste(cross.table[,"ipi"], cross.table[,"peptide"],sep=":"))
entry.levels <- levels( uniq.ipi.peptides )
## all_scan.table
all.scan.table <- read.table("all_scan.table", header=T)
## file name tags
cross.vec <- as.character(args)
ncross <- length(cross.vec)

## find all matched input files in current directory
input.path <- getwd()
mzXML.names <- list.files(path="../",pattern="mzXML$")
mzXML.files <- as.list( mzXML.names )
names(mzXML.files) <- mzXML.names
for (name in mzXML.names) {
  cat(paste(name,"\n",sep=""))
  mzXML.files[[name]] <- xcmsRaw( paste("../",name,sep=""), profstep=0, includeMSn=T)
}

##special case for raw file corruption
##ex2 <- mzXML.files[[2]]
##new.ex2 <- readFileFromMsn( ex2 )
##mzXML.files[[2]] <- new.ex2
##

## retention time window in secs
rt.window <- 5
rt.window.width <- rt.window * 60
local.rt.window <- 2
local.rt.window.width <- local.rt.window * 60
## signal/noise ratio for peak picking
sn <- 2.5

## column names for calculated ratios
integrated.area.ratio <- paste("IR",cross.vec,sep=".")
linear.regression.ratio <- paste("LR",cross.vec,sep=".")
linear.regression.R2 <- paste("R2",cross.vec,sep=".")
column.names <- c("index","ipi", "description", "symbol", "sequence", "mass", "charge", "segment",
                  integrated.area.ratio, linear.regression.ratio, linear.regression.R2, "entry", "link" )
out.df <- matrix(nrow=0, ncol=length(column.names))
colnames(out.df) <- column.names

## output layout
out.filename.base <- paste("output_rt_",as.character(rt.window),"_sn_",
                      as.character(sn),sep="")
out.filename <- paste(output.path,"/",out.filename.base,".ps",sep="")

postscript( out.filename, horizontal=F)
layout.vec <- seq(1,ncross)
layout.matrix <- cbind(2*layout.vec-1, 2*layout.vec-1, 2*layout.vec)
layout(layout.matrix)
par(oma=c(0,0,5,0), las=0)

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
  peptide.vec <- unlist( strsplit(peptide,"",fixed=T) )
  modified.mass <- sum( probe.mass[peptide.vec[which(peptide.vec=="*")-1]] )
  mono.mass  <- cross.table[i,"mass"]+modified.mass
  mass <- mono.mass + (which.max(isotope.dist( averagine.count(mono.mass) )) - 1)*isotope.mass.unit
  raw.scan.num <- cross.table[i,cross.vec]

  ## mz
  mz.light <- mass/charge + Hplus.mass
  mz.heavy <- (mass+pair.mass.delta)/charge + Hplus.mass
  ## scan number
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

  r2.v <- l.ratios <- i.ratios <- rep(NA,ncross)
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
    if ( !is.na(tag.rt) ) {
      all.ms2.scan <- as.integer( all.scan.table[(key==all.scan.table[,"key"]
                                                  &cross.vec[j]==all.scan.table[,"run"]),"scan"] )
      for (k in 1:length(all.ms2.scan)) {
        k.ms1.scan <- which(xfile@acquisitionNum > all.ms2.scan[k])[1]-1
        k.ms1.rt <- xfile@scantime[k.ms1.scan]/60
        points(k.ms1.rt, raw.ECI.light[[2]][k.ms1.scan], type='p',cex=0.5, pch=1)
      }
      ##lines(c(tag.rt,tag.rt),c(0.0, max(raw.ECI.light[[2]],raw.ECI.heavy[[2]])), col="green")
      points(tag.rt, raw.ECI.light[[2]][tag.ms1.scan.num], type='p',pch=8)

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

    best.r2 <- best.npoints <- best.ratio <- 0.0
    best.xlm <- best.light.yes <- best.heavy.yes <- best.low <- best.high <- NULL
    best.fixed <- F
    if (n.peaks>0) {
      for (n in 1:n.peaks) {
        low <- peaks[2*n-1]
        high<- peaks[2*n]
        yes <- which( raw.ECI.light.rt>=low & raw.ECI.light.rt<=high )
        light.yes <- raw.ECI.light[[2]][yes]
        heavy.yes <- raw.ECI.heavy[[2]][yes]
        ## calculate ratio of integrated peak area
        ratio <- round(sum(light.yes)/sum(heavy.yes),digits=2)
        if (ratio > 15 | ratio <0.4) next
        lines(c(low,low),ylimit/10, col="green")
        lines(c(high,high),ylimit/10, col="green")
        text(mean(c(low,high)),max(light.yes,heavy.yes)*1.2,
             labels=format(ratio,digits=4))
        ## calculate peak co-elution profile using only points above noise line
        ##yes2 <- light.yes > noise.light & heavy.yes > noise.heavy
        ##light.yes <- light.yes[yes2]
        ##heavy.yes <- heavy.yes[yes2]
        npoints <- length(light.yes)
        if (npoints<5) {
          next
        }
        x.lm <- lsfit( x=heavy.yes, y=light.yes,intercept=F )
        r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
        if ( !is.na(tag.rt) & tag.rt>=low & tag.rt<=high) {
          best.fixed <- T
        }
        if ( best.fixed | npoints > best.npoints ) {
          best.npoints <- npoints
          best.r2 <- r2
          best.ratio <- ratio
          best.xlm <- round(as.numeric(ls.print(x.lm,print.it=F)$coef.table[[1]][,"Estimate"]),digits=2)
          best.low <- low
          best.high <- high
          best.light.yes <- light.yes
          best.heavy.yes <- heavy.yes
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
      l.ratios[j] <- best.xlm
      r2.v[j] <- best.r2
    } else {
      plot(0,0,xlab="",ylab="",main=paste("R2 value: 0.00") )
    }
  }
  tt <- paste("Entry ", as.character(i), "-  Charge: ", as.character(charge),
              " - M/Z: ", as.character(format(mz.light, digits=7)),
              "and", as.character(format(mz.heavy,digits=7)))
  mtext(tt, line=3.5, outer=T)
  mtext(paste(cross.table[i,"peptide"],"; Mono.mass: ", as.character(mono.mass), sep=""),
        cex=0.8, line=2, out=T)
  mtext(paste(cross.table[i,"ipi"],description),line=0.8, cex=0.8,out=T)
  ## save data in outdf
  this.df <- c(i, ipi, description, symbol, peptide, round(mass,digits=4), charge, segment,
               i.ratios, l.ratios, r2.v, entry.index,
               paste('=HYPERLINK(\"./PNG/', out.filename.base,'-', as.character(i-1),'.png\")',sep=''))
  names(this.df) <- column.names
  out.df <- rbind(out.df, this.df)
} ## each entry
dev.off()

all.table <- out.df
all.table.out <- all.table[F,]
rsq.cutoff <- 0.8

## go from high concentration to low concentration,
## first apply R2 cutoff and sort by IR values
for ( s in seq(ncross, 1) ) {
  colname.R2 <- linear.regression.R2[s]
  colname.IR <- integrated.area.ratio[s]
  rsq.filter <- all.table[,colname.R2] >= rsq.cutoff & !is.na(all.table[,colname.R2])
  table <- all.table[rsq.filter,]
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
write.table(all.table.out,file=paste(output.path,"/",out.filename.base,".to_excel.txt",sep=""), quote=F, sep="\t", row.names=F,na="0.00")