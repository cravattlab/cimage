library(xcms)

table1 <- read.table("file.list",header=T)
table2 <- read.table("mass.list",header=T)

files <- as.list(as.character(table1[,"file"]))
names(files) <- names <- files
charges <- seq(1,6)
masses <- table2[,"mass"]
ppm <- 10*1e-6
Hplus <- 1.0072765

#postscript("output.ps",horizontal=T)
pdf("output.pdf",width=10, height=8, paper="letter")
par(mfrow=c(3,2), oma=c(0,0,5,0))
for (i in 1:length(files)) {
  files[[i]] <- xcmsRaw( names(files)[i], profstep=0)
}

for ( m in 1:length(masses) ) {
  mass <- masses[m]
  for ( i in 1:length(files) ) {
    xfile <- files[[i]]
    name <- names[[i]]
    for ( j in 1:length(charges) ){
      charge <- charges[j]
      mz <- mass/charge + Hplus
      mr <- c( mz*(1-ppm), mz*(1+ppm) )
      raw.ECI <- rawEIC(xfile, mr)
      raw.ECI.rt <- xfile@scantime[raw.ECI[[1]]]/60
      ylimit <- max(raw.ECI[[2]])
      plot(raw.ECI.rt, raw.ECI[[2]], type="l",xlab="RT(min)", ylab="intensity",
           main=paste("mz: ", as.character(format(mz, digits=7)), "; charge: ", j,
             "; NL: ", formatC(ylimit, digits=2, format="e")))
    }
    mtext(paste("Mass of ", as.character(format(mass,digits=7)), " in ", name), outer=T)
  }
}

dev.off()

