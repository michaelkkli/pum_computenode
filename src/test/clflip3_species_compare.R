# https://r-forge.r-project.org/projects/fda/
# In R,
# install.packages("packagename",repos="http://R-Forge.R-project.org")

#library(zoo, lib.loc='/home/mike/R/x86_64-pc-linux-gnu-library/2.7')
library(fda, lib.loc='/home/mike/R/x86_64-pc-linux-gnu-library/2.8')
#library(fda)

thesisdetail <- T

clflip3<-read.table("pum_convergence-3_species_compare_all.txt", header=F)

times<-c(0:71)*0.754

# Structured data. First is time point, second is repetition,
# third is different function.
sdata<-array(0,dim=c(72,5,8))

# TP Region 2 FLIP
sdata[,1,1]<-clflip3[[2]]
sdata[,2,1]<-clflip3[[3]]
sdata[,3,1]<-clflip3[[4]]
sdata[,4,1]<-clflip3[[5]]
sdata[,5,1]<-clflip3[[6]]
# TP Null
sdata[,1,2]<-clflip3[[7]]
sdata[,2,2]<-clflip3[[8]]
sdata[,3,2]<-clflip3[[9]]
sdata[,4,2]<-clflip3[[10]]
sdata[,5,2]<-clflip3[[11]]
# HC Region 2 FLIP
sdata[,1,3]<-clflip3[[12]]
sdata[,2,3]<-clflip3[[13]]
sdata[,3,3]<-clflip3[[14]]
sdata[,4,3]<-clflip3[[15]]
sdata[,5,3]<-clflip3[[16]]
# HC Null
sdata[,1,4]<-clflip3[[17]]
sdata[,2,4]<-clflip3[[18]]
sdata[,3,4]<-clflip3[[19]]
sdata[,4,4]<-clflip3[[20]]
sdata[,5,4]<-clflip3[[21]]
# F1 Region 2 FLIP
sdata[,1,5]<-clflip3[[22]]
sdata[,2,5]<-clflip3[[23]]
sdata[,3,5]<-clflip3[[24]]
sdata[,4,5]<-clflip3[[25]]
sdata[,5,5]<-clflip3[[26]]
# F1 Null
sdata[,1,6]<-clflip3[[27]]
sdata[,2,6]<-clflip3[[28]]
sdata[,3,6]<-clflip3[[29]]
sdata[,4,6]<-clflip3[[30]]
sdata[,5,6]<-clflip3[[31]]
# F2 Region 2 FLIP
sdata[,1,7]<-clflip3[[32]]
sdata[,2,7]<-clflip3[[33]]
sdata[,3,7]<-clflip3[[34]]
sdata[,4,7]<-clflip3[[35]]
sdata[,5,7]<-clflip3[[36]]
# F2 Null
sdata[,1,8]<-clflip3[[37]]
sdata[,2,8]<-clflip3[[38]]
sdata[,3,8]<-clflip3[[39]]
sdata[,4,8]<-clflip3[[40]]
sdata[,5,8]<-clflip3[[41]]

for (k in 1:8) {
    for (j in 1:5) {
        curmax<-sdata[1,j,k]
	for (i in 1:72) {
            sdata[i,j,k]<-sdata[i,j,k]/curmax
	}
    }
}

# column width 3.228in ~ 82 mm in emboj with aspect ratio 4:3
# x11(width=3.228,height=2.421)

# column width 6.77in ~172 mm with aspect ratio 4:3
# x11(width=6.77,height=5.0775)

timerange<-c(0.0,71*0.754)
fbasis<-create.polygonal.basis( c(0:71) * 0.754 )

sdatafd<-data2fd( sdata, times, fbasis )
meanfd<-mean(sdatafd)
stddevfd<-stddev.fd(sdatafd)

meanvals<-array(0.0,dim=c(72,8))
stddevvals<-array(0.0,dim=c(72,8))

for (i in 1:8) {
  meanvals[,i]<-eval.fd(times,meanfd[,i])
  stddevvals[,i]<-eval.fd(times,stddevfd[i])
}

diffvals<-array(0.0, dim=c(72,8) )
diffstddevvals<-array(0.0, dim=c(72,4) )

for (i in 1:4) {
  diffvals[,i]<-meanvals[,2*i-1]+1-meanvals[,2*i]
  diffstddevvals[,i]<-sqrt( stddevvals[,2*i-1]^2 + stddevvals[,2*i]^2 )
}

if ( meanvals[1,1]<1.1 ) {
  meanvals<-100*meanvals
  stddevvals<-100*stddevvals
  diffvals<-100*diffvals
  diffstddevvals<-100*diffstddevvals
}

## flipcol="gray30"
## obscol="gray70"
## corrcol="gray90"

flipcol="black"
barflipcol="gray40"
obscol="white"
corrcol="black"

## lty 0 blank, 1 solid, 2 dashed, 3 dotted, 4 dotdash, 5 longdash, 6 twodash

obsline=6
flipline=2

# Mark out some of the points.
illustpts=c(seq(1,72,3),72)

## pdf(file="preFig_2.pdf", width=6.77,height=6.77)
pdf(file="prespecies_compare.pdf", width=6.77,height=5.77)

topmar = 3.1

par(mfrow=c(2,2), oma=c(1.1,1.1,0,0),cex=0.75)

par(xaxt='n',yaxt='s')
#par(mar=c(0,4.1,4.1,0),oma=c(0,0,0,0))
par(mar=c(0,4.1,topmar,0))
#TP
k=1;
first=2*k-1; second=2*k
if ( thesisdetail ) {
##  Original version with standard deviation lines.
  matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5),
          xlab="", ylab="",
          col=1, lty=c(flipline,1,flipline,obsline,1,obsline), las=1,
          cbind( meanvals[,first]-stddevvals[,first],
                 meanvals[,first],
                 meanvals[,first]+stddevvals[,first],
                 meanvals[,second]-stddevvals[,second],
                 meanvals[,second],
                 meanvals[,second]+stddevvals[,second] )
  )
} else {
matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5),
        xlab="", ylab="",
        col=1, lty=c(1,1), las=1,
        cbind( meanvals[,first],
               meanvals[,second] )
)
}
# See https://stat.ethz.ch/pipermail/r-help/2007-October/144596.html
# [R] pdf() device uses fonts to represent points - data alteration?
# For the reason why almost opaque circles are used.
# Using ghostscript to embed fonts causes a problem.
# 2009-02-09
matpoints( times[illustpts], meanvals[illustpts,first], pch=20+k, col=rgb(0,0,0,0.99), bg=flipcol )
matpoints( times[illustpts], meanvals[illustpts,second], pch=20+k, col=rgb(0,0,0,0.99), bg=obscol )
if ( thesisdetail ) {
  legend( 27, 102, c("TP-GFP Obs.","Obs. s.d.", "TP-GFP FLIP","FLIP s.d."),
         col=rgb(0,0,0,0.99), pt.bg=c(obscol,obscol,flipcol,flipcol), lty=c(1,obsline,1,flipline), pch=c(20+k,NA,20+k,NA), cex=0.8 )
}
mtext(" A",font=2, side=1, adj=0, line=-1, outer=F)

par(xaxt='n',yaxt='n')
par(mar=c(0,0,topmar,3.1))

#HC
k=2;
first=2*k-1; second=2*k
if ( thesisdetail ) {
##  Original before removal of standard deviation lines.
  matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5), xlab="", ylab="",
          col=1, lty=c(flipline,1,flipline,obsline,1,obsline),
          cbind( meanvals[,first]-stddevvals[,first],
                 meanvals[,first],
                 meanvals[,first]+stddevvals[,first],
                 meanvals[,second]-stddevvals[,second],
                 meanvals[,second],
                 meanvals[,second]+stddevvals[,second] )
  )
} else {
matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5), xlab="", ylab="",
        col=1, lty=c(1,1),
        cbind( meanvals[,first],
               meanvals[,second] )
)
}
matpoints( times[illustpts], meanvals[illustpts,first], pch=20+k, col=1, bg=flipcol )
matpoints( times[illustpts], meanvals[illustpts,second], pch=20+k, col=1, bg=obscol )

if ( thesisdetail ) {
  legend( 27, 102, c("Hcf106-GFP Obs.","Obs. s.d.", "Hcf106-GFP FLIP","FLIP s.d."),
          col=c(1,1), pt.bg=c(obscol,obscol,flipcol,flipcol), lty=c(1,obsline,1,flipline), pch=c(20+k,NA,20+k,NA), cex=0.8 )
mtext(" B",font=2, side=1, adj=0, line=-1, outer=F)
}

par(xaxt='s',yaxt='s')
par(mar=c(3.1,4.1,0,0))

#F1
k=3;
first=2*k-1; second=2*k
if ( thesisdetail ) {
##  Original before removal of standard deviation lines.
   matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5), xlab="", ylab="",
          col=1, lty=c(flipline,1,flipline,obsline,1,obsline), las=1,
          cbind( meanvals[,first]-stddevvals[,first],
                 meanvals[,first],
                 meanvals[,first]+stddevvals[,first],
                 meanvals[,second]-stddevvals[,second],
                 meanvals[,second],
                 meanvals[,second]+stddevvals[,second] )
  )
} else {
matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5), xlab="", ylab="",
        col=1, lty=c(1,1), las=1,
        cbind( meanvals[,first],
               meanvals[,second] )
)
}
matpoints( times[illustpts], meanvals[illustpts,first], pch=20+k, col=1, bg=flipcol )
matpoints( times[illustpts], meanvals[illustpts,second], pch=20+k, col=1, bg=obscol )

if ( thesisdetail ) {
  legend( 27, 102,
          c( expression(paste("p23k",Delta,"TPP-GFP Obs.")), "Obs. s.d.",
             expression(paste("p23k",Delta,"TPP-GFP FLIP")), "FLIP s.d."),
          col=c(1,1), pt.bg=c(obscol,obscol,flipcol,flipcol), lty=c(1,obsline,1,flipline), pch=c(20+k,NA,20+k,NA), cex=0.8 )
}

mtext(" C",font=2, side=1, adj=0, line=-1, outer=F)

par(xaxt='s',yaxt='n')
par(mar=c(3.1,0,0,3.1))

#F2
k=4;
first=2*k-1; second=2*k
if ( thesisdetail ) {
##  Original before deletion of standard deviation lines.
  matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5), xlab="", ylab="",
          col=1, lty=c(flipline,1,flipline,obsline,1,obsline),
          cbind( meanvals[,first]-stddevvals[,first],
                 meanvals[,first],
                 meanvals[,first]+stddevvals[,first],
                 meanvals[,second]-stddevvals[,second],
                 meanvals[,second],
                 meanvals[,second]+stddevvals[,second] )
  )
} else {
matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5), xlab="", ylab="",
        col=1, lty=c(1,1),
        cbind( meanvals[,first],
               meanvals[,second] )
)
}
matpoints( times[illustpts], meanvals[illustpts,first], pch=20+k, col=1, bg=flipcol )
matpoints( times[illustpts], meanvals[illustpts,second], pch=20+k, col=1, bg=obscol )
if ( thesisdetail ) {
  legend( 27, 102,
          c( "p23k-GFP Obs.", "Obs. s.d.",
             "p23k-GFP FLIP", "FLIP s.d."),
          col=c(1,1), pt.bg=c(obscol,obscol,flipcol,flipcol), lty=c(1,obsline,1,flipline), pch=c(20+k,NA,20+k,NA), cex=0.8 )
}
mtext(" D",font=2, side=1, adj=0, line=-1, outer=F)

mtext("Percentage fluorescence", side=2, outer=T, line=-1.3, cex=0.9)
mtext("Time / seconds", side=1, outer=T, line=-0.7, cex=0.9)

dev.off()
embedFonts( "prespecies_compare.pdf", "pdfwrite", "species_compare.pdf", options="-dPDFSETTINGS=/prepress -dEmbedAllFonts=true" )


stop("First part of R plotting is complete. Exiting.")

# Calculate slopes

slopes<-array(0.0, dim(sdata)[2:3])
for ( f in 1:8 ) {
  for ( r in 1:5 ) {
    slopes[r,f]<-coef(lm( sdata[,r,f]~times,subset=13:52))[[2]]
  }
}

## Setup dataframe for easy ANOVA.
## h Hcf, t TP, p p23, q p23deltaTPP
## s spot photobleaching, o observational photobleaching
dataframe<-data.frame(construct=c(array("h",10),array("t",10),array("p",10),array("q",10)),
                      photobleach=array(c(array("s",5),array("o",5)),40),
                      lossrate=c(slopes[1:5,1],
                                 slopes[1:5,2],
                                 slopes[1:5,3],
                                 slopes[1:5,4],
                                 slopes[1:5,5],
                                 slopes[1:5,6],
                                 slopes[1:5,7],
                                 slopes[1:5,8]))
summary(aov(lossrate~construct*photobleach,dataframe))


# Arrange slopemeans to allow easy plotting with barplot.
slopemeans<-array(0.0, c(2,dim(sdata)[3]/2) )
slopesd<-array(0.0, c(2, dim(sdata)[3]/2) )
for ( c in 1:4 ) {
  slopemeans[1,c]<-mean(slopes[,2*c-1])
  slopemeans[2,c]<-mean(slopes[,2*c  ])
  slopesd[1,c]<-sd(slopes[,2*c-1])
  slopesd[2,c]<-sd(slopes[,2*c  ])
}


## 4x4 panel
## pdf(file="preFig_3.pdf", width=6.77,height=6.77 )

pdf(file="preFig_3.pdf", width=6.77,height=3.39 )
par(mfrow=c(1,2))

## Show or hide 'n' axis.
par(xaxt='s',yaxt='s')
par(mar=c(4.1,4.1,2.1,0.6),cex=0.75)
# FLIP+1-Obs.
matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(78,103), yaxp=c(70,100,3),
        xlab="", ylab="",
        col=1, lty=c(1,1,1,1), las=1,
        cbind( meanvals[,1]+100-meanvals[,2],
               meanvals[,3]+100-meanvals[,4],
               meanvals[,5]+100-meanvals[,6],
               meanvals[,7]+100-meanvals[,8] )
)
matpoints( times[illustpts], meanvals[illustpts,1]+100-meanvals[illustpts,2], pch=20+1, col=rgb(0,0,0,0.99), bg=corrcol )
matpoints( times[illustpts], meanvals[illustpts,3]+100-meanvals[illustpts,4], pch=20+2, col=1, bg=corrcol )
matpoints( times[illustpts], meanvals[illustpts,5]+100-meanvals[illustpts,6], pch=20+3, col=1, bg=corrcol )
matpoints( times[illustpts], meanvals[illustpts,7]+100-meanvals[illustpts,8], pch=20+4, col=1, bg=corrcol )
if ( thesisdetail ) {
  legend( 23, 102, c("TP-GFP","Hcf106-GFP", expression(paste("p23k",Delta,"TPP-GFP")),"p23k-GFP"),
          col=rgb(0,0,0,0.99), pt.bg=corrcol, lty=c(1,1,1,1), pch=c(21,22,23,24), cex=0.7 )
}
##mtext(" A",font=2, side=1, adj=0, line=-1, outer=F)
mtext("Percentage fluorescence", side=2, outer=F, line=2.5, cex=0.9)
mtext("Time / seconds",          side=1, outer=F, line=2.5, cex=0.9)


par(mar=c(5.1,4.1,4.1,1.1))
prevpar<-par(xaxt="s", yaxt="s",cex=0.55)
labels=c("TP-GFP", "Hcf106-GFP", expression(paste("p23k",Delta,"TPP-GFP")), "p23k-GFP")
cnstorder=c(1,4,2,3)
barplot(-100*slopemeans[,cnstorder],
        ylim=c(0,0.8), ylab="",
        names.arg=labels[cnstorder], cex.names=1.0,
        beside=T, offset=0, col=c(barflipcol,obscol),
        legend.text=c("FLIP photobleaching", "Obs. photobleaching")  )
mtext("Percentage fluorescence loss per second", side=2, outer=F, line=2.4, cex=0.7)
barpositions=c(1.5,4.5,7.5,10.5)
arrows(barpositions,
       -100*slopemeans[1,cnstorder]+100*slopesd[1,cnstorder],
       barpositions,
       -100*slopemeans[1,cnstorder]-100*slopesd[1,cnstorder],
       angle=90, code=3, length=0.05, lw=1 )
barpositions<-barpositions+1
arrows(barpositions,
       -100*slopemeans[2,cnstorder]+100*slopesd[2,cnstorder],
       barpositions,
       -100*slopemeans[2,cnstorder]-100*slopesd[2,cnstorder],
       angle=90, code=3, length=0.05, lw=1 )
##mtext(" B",font=2, at=-0.5, side=1, adj=0, line=3, outer=F)
mtext(" A", font=2, at=0.02, side=1, adj=0, line=-2, outer=T)
mtext(" B", font=2, at=0.52, side=1, adj=0, line=-2, outer=T)


## par(prevpar)
## par(xaxt='s',yaxt='s')
## par(mar=c(5.1,4.1,0,0))
## #FLIP
## matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5),
##         xlab="", ylab="",
##         col=1, lty=c(1,1,1,1), las=1,
##         cbind( meanvals[,1],
##                meanvals[,3],
##                meanvals[,5],
##                meanvals[,7] )
## )
## matpoints( times[illustpts], meanvals[illustpts,1], pch=20+1, col=rgb(0,0,0,0.99), bg=flipcol )
## matpoints( times[illustpts], meanvals[illustpts,3], pch=20+2, col=1, bg=flipcol )
## matpoints( times[illustpts], meanvals[illustpts,5], pch=20+3, col=1, bg=flipcol )
## matpoints( times[illustpts], meanvals[illustpts,7], pch=20+4, col=1, bg=flipcol )
## legend( 23, 102, c("TP-GFP FLIP","Hcf106-GFP FLIP", expression(paste("p23k",Delta,"TPP-GFP FLIP")),"p23k-GFP FLIP"),
##         col=rgb(0,0,0,0.99), pt.bg=flipcol, lty=c(1,1,1,1), pch=c(21,22,23,24), cex=0.7 )
## mtext(" C",font=2, side=1, adj=0, line=-1, outer=F)


## par(xaxt='s',yaxt='n')
## par(mar=c(5.1,0,0,2.1))
## # Obs.
## matplot(times, type="l", xlim=c(-3,61), xaxp=c(0,60,6), ylim=c(48,103), yaxp=c(50,100,5),
##         xlab="", ylab="",
##         col=1, lty=c(1,1,1,1), las=1,
##         cbind( meanvals[,2],
##                meanvals[,4],
##                meanvals[,6],
##                meanvals[,8] )
## )
## matpoints( times[illustpts], meanvals[illustpts,2], pch=20+1, col=rgb(0,0,0,0.99), bg=obscol )
## matpoints( times[illustpts], meanvals[illustpts,4], pch=20+2, col=1, bg=obscol )
## matpoints( times[illustpts], meanvals[illustpts,6], pch=20+3, col=1, bg=obscol )
## matpoints( times[illustpts], meanvals[illustpts,8], pch=20+4, col=1, bg=obscol )
## legend( 23, 102, c("TP-GFP Obs.","Hcf106-GFP Obs.", expression(paste("p23k",Delta,"TPP-GFP Obs.")),"p23k-GFP Obs."),
##         col=rgb(0,0,0,0.99), pt.bg=obscol, lty=c(1,1,1,1), pch=c(21,22,23,24), cex=0.7 )
## mtext(" D",font=2, side=1, adj=0, line=-1, outer=F)
## # Hack to get label for difficult top right panel.
##mtext(" B",font=2, side=1, adj=0, line=-16.3, outer=F)



## mtext("Percentage fluorescence", side=2, outer=T, line=-1.3, cex=0.9)
## mtext("Time / seconds", side=1, outer=T, line=-2.3, cex=0.9)


dev.off()
embedFonts( "preFig_3.pdf", "pdfwrite", "Fig_3.pdf", options="-dPDFSETTINGS=/prepress -dEmbedAllFonts=true" )

