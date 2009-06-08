library(xcms)
source("/home/chuwang/svnrepos/R/msisotope.R")

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
mzXML.names <- list.files(path="./",pattern="mzXML$")
mzXML.files <- as.list( mzXML.names )
names(mzXML.files) <- mzXML.names
for (name in mzXML.names) {
  cat(paste(name,"\n",sep=""))
  mzXML.files[[name]] <- xcmsRaw( paste("./",name,sep=""), profstep=0, includeMSn=T)
}

##special case for raw file corruption
##ex2 <- mzXML.files[[2]]
##new.ex2 <- readFileFromMsn( ex2 )
##mzXML.files[[2]] <- new.ex2
##

## retention time window in secs
rt.window <- 5
rt.window.width <- rt.window * 60
## signal/noise ratio for peak picking
sn <- 2.5
out.filename.base <- paste("output_rt_",as.character(rt.window),"_sn_",
                      as.character(sn),sep="")
out.filename <- paste(output.path,"/",out.filename.base,".ps",sep="")
## make plot for each isotopic cluster
##png(out.filename)
out.df <- data.frame(index=numeric(0), ipi=character(0), description=character(0),
                     symbol=character(0), sequence=character(0),
                     mass=numeric(0), charge=numeric(0), segment=character(0),
                     ir.1.1=numeric(0), ir.1.5=numeric(0), ir.1.10=numeric(0),
                     lr.1.1=numeric(0), lr.1.5=numeric(0), lr.1.10=numeric(0),
                     rsq.1.1=numeric(0), rsq.1.5=numeric(0), rsq.1.10=numeric(0),
                     entry=numeric(0), link=character(0)
##                     ar.1.1=numeric(0), ar.1.5=numeric(0), ar.1.10=numeric(0),
##                     arsq.1.1=numeric(0), arsq.1.5=numeric(0), arsq.1.10=numeric(0),
)
#original.df <- ipi.name.table[F,]
postscript( out.filename, horizontal=F)
layout(matrix(c(1,1,2,3,3,4,5,5,6), 3, 3, byrow = T))
#par(mfrow=c(ncross,1))
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
  nmod <- length(strsplit(as.character(peptide),"*",fixed=T)[[1]])-1
  mono.mass  <- cross.table[i,"mass"]+nmod*probe.mass
  mass <- mono.mass + (which.max(isotope.dist( averagine.count(mono.mass) )) - 1)*isotope.mass.unit
  raw.scan.num <- cross.table[i,cross.vec]

  ## mz
  mz.light <- mass/charge + Hplus.mass
  mz.heavy <- (mass+pair.mass.delta)/charge + Hplus.mass
  ## scan number
  ms1.scan.num <- exist.index <- which( raw.scan.num > 0 )
  for ( k in 1:length(exist.index) ) {
    kk <- exist.index[k]
    raw.file <- paste( cross.vec[kk], "_", segment,".mzXML",sep="")
    xfile <- mzXML.files[[raw.file]]
    ms1.scan.num[k] <- which(xfile@acquisitionNum > as.integer(raw.scan.num[kk]))[1]-1
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
    ylimit[2] <- ylimit[2]*1.2
    ##xlimit <- range(c(raw.ECI.light[[1]], raw.ECI.heavy[[1]]))
    ##ylimit <- range(c(raw.ECI.light[[2]], raw.ECI.heavy[[2]]))
    xlimit <- c(rt.min,rt.max)/60
    raw.ECI.light.rt <- xfile@scantime[ raw.ECI.light[[1]] ] / 60
    raw.ECI.heavy.rt <- xfile@scantime[ raw.ECI.heavy[[1]] ] / 60
    ##if ( tag=="*") {
    tt.main <- paste(tag, raw.file, "; Raw Scan:", as.character(raw.scan.num[j]),
                     "; NL:", formatC(ylimit[2], digits=2, format="e"))
    ##} else {
    ##  tt.main <- paste(tag, raw.file, "; Census ratio: NA; Raw Scan: NA")
    ##}
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
    }
    ## guess ratio of integrated peak area
    peaks <- findPairChromPeaks( raw.ECI.light.rt, raw.ECI.light[[2]], raw.ECI.heavy[[2]],xlimit,sn )
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
        ##} else if ( npoints<3) {
        ##  light.yes <- c(1,light.yes)
        ##  heavy.yes <- c(1,heavy.yes)
        ##}
        ##          x.lm <- lsfit( x=log10(light.yes), y=log10(heavy.yes) )
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
      ##        plot(log10(best.light.yes),log10(best.heavy.yes),
      ##             xlab="log10(intensity.light)", ylab="log10(intensity.heavy)",
      ##             main=paste("R2 value:",format(best.r2,digits=3)),
      ##             xlim=range(log10(best.light.yes),log10(best.heavy.yes)),
      ##             ylim=range(log10(best.light.yes),log10(best.heavy.yes)))

      plot(best.heavy.yes,best.light.yes,
           xlab="intensity.heavy", ylab="intensity.light",
           main=paste("X=",format(best.xlm,digits=4),"; R2=",format(best.r2,digits=3),
             "; Np=", best.npoints, sep=""),
           xlim=c(0, max(best.light.yes,best.heavy.yes)),
           ylim=c(0, max(best.light.yes,best.heavy.yes)))
      ##abline(best.xlm[1],best.xlm[2])
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
  ##original.df <- rbind( original.df, ipi.name.table[i,] )
  this.df <- data.frame(index=i, ipi=ipi, description=description, symbol=symbol, sequence=peptide,
                        mass=round(mass,digits=4), charge=charge, segment=segment,
                        ir.1.1=i.ratios[1], ir.1.5=i.ratios[2], ir.1.10=i.ratios[3],
                        lr.1.1=l.ratios[1], lr.1.5=l.ratios[2], lr.1.10=l.ratios[3],
                        rsq.1.1=r2.v[1], rsq.1.5=r2.v[2], rsq.1.10=r2.v[3],
                        ##                          ar.1.1=0, ar.1.5=0, ar.1.10=0,
                        ##                          arsq.1.1=0, arsq.1.5=0, arsq.1.10=0,
                        entry=entry.index, link=paste('=HYPERLINK(\"./PNG/', out.filename.base,'-',
                                             as.character(i-1),'.png\")',sep=''))
  out.df <- rbind(out.df, this.df)
} ## each entry
dev.off()

##for (k in levels(as.factor(out.df[,"entry"])) ) {
##  kk <- which(as.factor(out.df[,"entry"]) == k )
## for ( m in 1:ncross ) {
##    ## average non-zero ratio
##    v <- out.df[kk,m+3]
##    v <- v[v>0]
##   if ( length(v)>0) {
##     out.df[kk,m+9] <- round(mean(v),digits=2)
##   } else {
##     out.df[kk,m+9] <- 0.0
##   }
## average non-zero R2
##   v <- out.df[kk,m+6]
##   v <- v[v>0]
##   if ( length(v)>0) {
##     out.df[kk,m+12] <- round(mean(v),digits=2)
##   } else {
##     out.df[kk,m+12] <- 0.0
##   }
## }
##}


all.table <- out.df##cbind(original.df,out.df)

rsq.cutoff <- 0.9
## first apply rsq.1.10 cutoff and sort by ir.1.10
rsq.filter1 <- all.table[,"rsq.1.10"] >= rsq.cutoff & !is.na(all.table[,"rsq.1.10"])
table1 <- all.table[rsq.filter1,]
s1 <- table1[,"ir.1.10"]
s2 <- table1[,"entry"]
s3 <- table1[,"charge"]
s4 <- table1[,"segment"]
ii <- order(s1, s2, s3, s4)
table1 <- table1[ii,]
## then select by rsq.1.5 cutoff and sort by ir.1.5
all.table <- all.table[!rsq.filter1,]
rsq.filter2 <- all.table[,"rsq.1.5"] >= rsq.cutoff & !is.na(all.table[,"rsq.1.5"])
table2 <- all.table[rsq.filter2,]
s1 <- table2[,"ir.1.5"]
s2 <- table2[,"entry"]
s3 <- table2[,"charge"]
s4 <- table2[,"segment"]
ii <- order(s1, s2, s3, s4)
table2 <- table2[ii,]
## then select by rsq.1.1 cutoff and sort by ir.1.1
all.table <- all.table[!rsq.filter2,]
rsq.filter3 <- all.table[,"rsq.1.1"] >= rsq.cutoff & !is.na(all.table[,"rsq.1.1"])
table3 <- all.table[rsq.filter3,]
s1 <- table3[,"ir.1.1"]
s2 <- table3[,"entry"]
s3 <- table3[,"charge"]
s4 <- table3[,"segment"]
ii <- order(s1, s2, s3, s4)
table3 <- table3[ii,]
##
table4 <- all.table[!rsq.filter3,]

all.table.out <- rbind(table1, table2, table3, table4)
#ii.zero <- (all.table.sort[,"ar.1.10"] == 0.0)
#all.table.out <- rbind( all.table.sort[!ii.zero,], all.table.sort[ii.zero,] )
row.names(all.table.out) <- as.character(seq(1:dim(all.table.out)[1]) )
all.table.out[,"index"] <- seq(1:dim(all.table.out)[1])
write.table(all.table.out,file=paste(output.path,"/",out.filename.base,".to_excel.txt",sep=""), quote=F, sep="\t", row.names=F,na="0.00")
