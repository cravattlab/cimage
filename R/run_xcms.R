### I/O options ###
folder1 <- "wt"
folder2 <- "b"
result <- "results"
###################
library(xcms)
xset<-xcmsSet()
xset<-group(xset)
xset.list <- list()
xset.list[[length(xset.list)+1]] <- xset.cur <- xset

last.rtcg.num <- rtcg.num <- 0
do.retcor <- T
while( do.retcor ) {
  Sys.sleep(1)
  last.rtcg.num <- rtcg.num
  sink("xcms.rtcg.out",split=T)
  xset.new <- retcor(xset.cur, family="s", plottype="m")
  sink()
  rtcg.sink <- read.table("xcms.rtcg.out",header=F)
  rtcg.num <- rtcg.sink[1,5]
  xset.new <- group(xset.new,bw=10)
  xset.list[[length(xset.list)+1]] <- xset.cur <- xset.new
  do.retcor <- ( last.rtcg.num != rtcg.num )
}

xset.new <-fillPeaks(xset.cur)
reporttab<-diffreport(xset.new, folder1, folder2, result, 500, metlin=0.15)
