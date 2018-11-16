#
#
# CODE TO REPRODUCE FIGURES AND STATISTICS FOR
#   LINKED SELECTION MANUSCRIPT
#

# ENTER THE PATH TO THE MAIN DIRECTORY FOR THIS PROJECT

require(hexbin)
baseDIR    <-    "/Users/thom/Dropbox/1.Shared_with_Bill/9.Manuscripts/2.LinkedVariation/"

# maps    <-    read.table(
              # paste0(baseDIR,"/4.data/1.datasets/maps_old.tsv"), sep = "\t",header= TRUE,
              # stringsAsFactors = FALSE)

maps    <-    read.table(
              paste0(baseDIR,"/4.data/1.datasets/all_maps.tsv"), sep = "\t",header= TRUE,
              stringsAsFactors = FALSE)
        maps$Mb <- maps$bp / 1000000
        maps    <- maps[maps$chr %in% 1:21,]


popgen    <-    read.table(paste0(baseDIR,
                    "/4.data/1.datasets/polymorphism.tsv"), 
                    sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    popgen            <-    popgen[popgen$chr %in% as.character(1:21),]
    popgen$Mb        <-    popgen$bp / 1000000
    chr                <-    as.character(popgen$chr)
    bp                <-    as.character(popgen$bp)
    coord            <-    paste0(chr,"_",bp)
    popgen$coord    <-    coord
    popgen    <-    popgen[order(popgen$chr, popgen$bp),]
    popgen$fst.rsbt[popgen$fst.rsbt > 1]  <- 0



RS21.ref.order    <-    read.table(paste0(baseDIR,"/4.data/1.datasets/RS21_unflipped.tsv"),
                                sep = '\t', header = T)

LGs.info    <-    read.table(paste0(baseDIR,"/4.data/1.datasets/GacNewLGs.tsv"), sep = "\t", 
                        header = TRUE, stringsAsFactors = FALSE)
group    <-    2    # ***column in which to find 'group' label (You specify)
bp        <-    3    # ***column in which to find chromosome-specific bp position (You specify)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                    ###
###                  PERCENT OF REF GENOME COVERED                     ###
###                                                                    ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# FOR EACH CHROMOSOME, GET THE TOTAL LENGTH FROM THE REFERENCE ASSEMBLY
#   AND FIND THE FIRST/LAST POSITION OF MAP MARKERS

# DATA FRAME: CHR, LEN, SPAN.BT, SPAN.RS, SPAN.H

crosses    <-    c("Bt_Subset","RS","Hybrid_Subset")
chrs    <-    1:21

ref.cov    <-    matrix(nrow = 21, ncol = 4)

for (i in chrs)    {
    result    <-    1:4
    chr    <-    i
    len    <-    LGs.info$LGs.len[LGs.info$groups == chr]
    result[1]    <-    len
    for (j in 1:3) {
        x    <-    crosses[j]
        c    <-    j + 1
        l    <-    maps$bp[maps$map == x & maps$chr == chr]
        l    <-    range(l) ; l    <-    l[2] - l[1]
        result[c]    <-    l
    }
    ref.cov[i,]    <-    result
}
ref.cov    <-    data.frame(ref.cov)
names(ref.cov)    <-    c("len","span.bt","span.rs","span.h")
ref.cov            <-    rbind(ref.cov, apply(ref.cov, 2, sum))


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                    ###
###                       PLOT PHYS/GENET DIST                         ###
###                                                                    ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


data    <-    maps[maps$map == "RS",]
data    <-    maps[maps$map == "Bt_Subset",]
data    <-    maps[maps$map == "Hybrid_Subset",]

chrs    <-    as.character(1:21)

# PRINT OUT BIG-OL' PDF, ALL AXES SCALED EQUALLY
chrs    <-    1:21
plotfile    <-    "ThreeMaps_scaled_with19.pdf"
outdir        <-    paste0(baseDIR,"/3.figures/2.drafts/")
pdf(file = paste0(outdir,plotfile,sep=""), width = 6,height=12)
par(mfrow = c(7,3), mar = c(3,3,0,1)+0.1)
for (chr in chrs) {
    # chr="4"
    bt    <-    maps[maps$map == "Bt" & maps$chr == chr,]
    rs    <-    maps[maps$map == "RS" & maps$chr == chr,]
    hyb    <-    maps[maps$map == "Hybrid" & maps$chr == chr,]
    plot(bt$Mb, bt$cM, type = 'n', ylim = c(0,175), xlab = "", ylab = "", xlim = c(0,35),
         yaxt="n")
    axis(2, at=c(0,50,100,150))
    points(bt$Mb, bt$cM, col = rgb(45/255,104/255,183/255), pch=20, cex = 0.6)
    points(rs$Mb, rs$cM, col = rgb(180/255,19/255,52/255), pch=20, cex = 0.6)
    points(hyb$Mb, hyb$cM, col = rgb(105/255,19/255,99/255), pch=20, cex = 0.6)
    # legend(x='topleft', legend = paste0("chr ",chr,collapse=""), bty="n")
}
dev.off()



# PRINT OUT ALL CHRS MENTIONED IN MAIN TEXT: 1, 4, 11, 12, 14, 15, 19, 21
chrs    <-    c(1,4,11,12,14,15,19,21)
plotfile    <-    "ThreeMaps_scaled_mentionedInText_smallpoints.pdf"
outdir        <-    paste0(baseDIR,"/3.figures/2.drafts/")
pdf(file = paste0(outdir,plotfile,sep=""), width = 8,height=4)
par(mfrow = c(2,4), mar = c(3,3,0,1)+0.1)
for (chr in chrs) {
    # chr="4"
    bt    <-    maps[maps$map == "Bt" & maps$chr == chr,]
    rs    <-    maps[maps$map == "RS" & maps$chr == chr,]
    hyb    <-    maps[maps$map == "Hybrid" & maps$chr == chr,]
    plot(bt$Mb, bt$cM, type = 'n', ylim = c(0,175), xlab = "", ylab = "", xlim = c(0,35),
         yaxt="n")
    axis(2, at=c(0,50,100,150))
    points(bt$Mb, bt$cM, col = rgb(45/255,104/255,183/255), pch=20, cex = 0.1)
    points(rs$Mb, rs$cM, col = rgb(180/255,19/255,52/255), pch=20, cex = 0.1)
    points(hyb$Mb, hyb$cM, col = rgb(105/255,19/255,99/255), pch=20, cex = 0.1)
    # legend(x='topleft', legend = paste0("chr ",chr,collapse=""), bty="n")
}
dev.off()



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###  CHR21: REF ORDER, UNCORRECTED, CORRECTED   ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

RSuncor  <-    RS21.ref.order
RScor    <-    maps[maps$chr==21 & maps$map == "RS",]
BTref    <-    maps[maps$chr==21 & maps$map == "Bt",]
# using 'Bt' for visualization purposes only.

par(mfrow = c(1,3))
    plot(cM_01 ~ bp, BTref, pch = 20, col = 'blue3', xlim = c(9000000,12500000), ylim = c(2,12))
    plot(cM_01 ~ bp, RSuncor, pch = 20, col = 'red3', xlim = c(9000000,12500000), ylim = c(20,82))
    plot(cM_01 ~ bp, RScor, pch = 20, col = 'red3', xlim = c(9000000,12500000), ylim = c(6,18))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                    ###
###                         SUMMARY STATISTICS                         ###
###                                                                    ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


data    <-    maps[maps$map == "RS" ,]
data    <-    maps[maps$map == "Bt_Subset" ,]
data    <-    maps[maps$map == "Hybrid_Subset" ,]

chr.phys    <-    tapply(data$Mb, INDEX = data$chr, FUN = max)
chr.genet    <-    tapply(data$cM, INDEX = data$chr, FUN = max)

par(mfrow = c(1,1))
plot(chr.phys,chr.genet, pch = 21, cex = 6, bg = 'gray50',
     xlim = c(14,35), ylim = c(40,200),
     xlab = "Physical map length, Mb", ylab = "Genetic map length, cM")
text(x = chr.phys,y = chr.genet, labels = names(chr.phys), cex = 0.75)


# CHROMOSOME-WIDE VARIATION IN MAP LENGTHS

chr.summaries        <-    NULL
genet.by.pop        <-    NULL
rate.by.pop            <-    NULL

pop.name    <-    "Bt_Subset"
    len.phys    <-    tapply(maps$Mb[maps$map == pop.name],         
                             INDEX = maps$chr[maps$map == pop.name], FUN = max)
    len.genet    <-    tapply(maps$cM[maps$map == pop.name],
                              INDEX = maps$chr[maps$map == pop.name], FUN = max)
    pop    <-    rep(pop.name, length(len.phys))
    chr    <-    names(len.phys)
    d    <-    data.frame(pop,chr,len.phys,len.genet)
    chr.summaries    <-    d
    genet.by.pop    <-    cbind(genet.by.pop, len.genet)

pop.name    <-    "Hybrid_Subset"
    len.phys    <-    tapply(maps$Mb[maps$map == pop.name],         
                             INDEX = maps$chr[maps$map == pop.name], FUN = max)
    len.genet    <-    tapply(maps$cM[maps$map == pop.name],
                              INDEX = maps$chr[maps$map == pop.name], FUN = max)
    pop    <-    rep(pop.name, length(len.phys))
    chr    <-    names(len.phys)
    d    <-    data.frame(pop,chr,len.phys,len.genet)
    chr.summaries    <-    rbind(chr.summaries,d)
    genet.by.pop    <-    cbind(genet.by.pop, len.genet)

pop.name    <-    "RS"
    len.phys    <-    tapply(maps$Mb[maps$map == pop.name],         
                             INDEX = maps$chr[maps$map == pop.name], FUN = max)
    len.genet    <-    tapply(maps$cM[maps$map == pop.name],
                              INDEX = maps$chr[maps$map == pop.name], FUN = max)
    pop    <-    rep(pop.name, length(len.phys))
    chr    <-    names(len.phys)
    d    <-    data.frame(pop,chr,len.phys,len.genet)
    chr.summaries    <-    rbind(chr.summaries,d)
    genet.by.pop    <-    cbind(genet.by.pop, len.genet)

boxplot(len.genet ~ pop, chr.summaries, ylim = c(0,250),
        ylab = 'chromosome map length, cM', pch = 20)
    segments(x0=1, x1=2, y0=genet.by.pop[,1], y1=genet.by.pop[,2], col = rgb(0,0,0,0.5))
    segments(x0=2, x1=3, y1=genet.by.pop[,3], y0=genet.by.pop[,2], col = rgb(0,0,0,0.5))
    points(x = rep(1,length(genet.by.pop[,1])), y = genet.by.pop[,1], pch = 20)
    points(x = rep(2,length(genet.by.pop[,2])), y = genet.by.pop[,2], pch = 20)
    points(x = rep(3,length(genet.by.pop[,3])), y = genet.by.pop[,3], pch = 20)


    

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###
###      RECOMBINATION RATE/POPGEN CORRELATIONS
###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

### NON-OVERLAPPING GENOMIC WINDOWS

window.size    <-    100000 # in bp
windowed    <-    NULL
for (i in 1:21)    {
    #	i = 1
    cat(paste0("Scanning chromosome ",i,"\n"))
    m        <-    maps[maps$chr == i,]
    p        <-    popgen[popgen$chr == i,]
    m.rs    <-    m[m$map == "RS",]
    m.bt    <-    m[m$map == "Bt_Subset",]
    m.h        <-    m[m$map == "Hybrid_Subset",]
    # find the highest start position among maps to avoid "NA"s
    start.min    <-    max(m.rs$bp[1],m.bt$bp[1],m.h$bp[1], na.rm=TRUE)
    # find the lowest end position among maps to avoid "NA"s
    end.min        <-    min(m.rs$bp[length(m.rs$bp)],
                             m.bt$bp[length(m.bt$bp)],
                             m.h$bp[length(m.h$bp)], na.rm=TRUE)
    # identify all start positions, use only windows with full 250kb
    starts        <-    seq(start.min,end.min,by=window.size)
    n.windows    <-    length(starts)-1
    starts        <-    starts[1:n.windows]
    # record: chr, bp.start, bp.end, bp.mid, cMperMb.bt, cMperMb.rs, cMperMb.h,
    # (mean)  pi, pi.bt, pi.rs, dxy, fst
    # = 12 columns
    for (j in starts) {
        #    	j = starts[1]
        end    <-    j + (window.size - 1)
        bp.mid    <-    round((j+end)/2)
        
        cM.bt    <-    m.bt$cM[m.bt$bp >= j & m.bt$bp <= end]
        cMperMb.bt    <- ifelse(
            length(cM.bt) > 1, 
            (max(cM.bt, na.rm=T) - min(cM.bt, na.rm=T)) / (window.size / 1000000),
            ifelse(
                is.na(cM.bt),
                NA,
                cM.bt)
        )
        cM.rs    <-    m.rs$cM[m.rs$bp >= j & m.rs$bp <= end]
        cMperMb.rs    <- ifelse(
            length(cM.rs) > 1, 
            (max(cM.rs, na.rm=T) - min(cM.rs, na.rm=T)) / (window.size / 1000000),
            ifelse(
                is.na(cM.rs),
                NA,
                cM.rs)
        )      
        cM.h    <-    m.h$cM[m.h$bp >= j & m.h$bp <= end]
        cMperMb.h    <- ifelse(
            length(cM.h) > 1, 
            (max(cM.h, na.rm=T) - min(cM.h, na.rm=T)) / (window.size / 1000000),
            ifelse(
                is.na(cM.h),
                NA,
                cM.h)
        )
        place2b    <-    p$bp >= j & p$bp <= end
        pi         <-    mean(p$pi[place2b])
        pi.bt      <-    mean(p$pi.bt[place2b])
        pi.fw      <-    mean(p$pi.fw[place2b])
        pi.rs      <-    mean(p$pi.rs[place2b])
        dxy        <-    mean(p$dxy.rsbt[place2b])
        fst        <-    mean(p$fst.rsbt[place2b])
        tmrca      <-    mean(p$threespine.scaled[place2b])
        t.noD      <-    mean(p$threespine.scaled[place2b & p$cons.sort != 'D'])
        d.dens     <-    length(p$cons.sort[p$cons.sort == 'D' & place2b]) / length(p$cons.sort[place2b])
        result     <-    c(i, j, end, bp.mid, cMperMb.bt, cMperMb.rs, cMperMb.h,
                           pi, pi.bt,pi.fw, pi.rs, dxy, fst, tmrca, t.noD, d.dens)
        windowed    <-    rbind(windowed, result)
    }
    
}

rownames(windowed)    <-    NULL
windowed    <-    data.frame(windowed)
names(windowed)    <-    c("chr","bp.start","bp.end","bp.mid",
                           "cMperMb.bt","cMperMb.rs","cMperMb.h",
                           "pi","pi.bt","pi.fw","pi.rs","dxy","fst", "tmrca", "t.noD","d.dens")
windowed$da <-windowed$dxy - ((windowed$pi.rs + windowed$pi.fw)/2)

windowed$d.dens.scaled	<- scale.vals(windowed$d.dens, 0,1)
point.col               <- rgb(0.5,0.5,0.5,0.5)
inv.col                 <- rgb(0.95,0.25,0.25,0.75)
windowed$inv.col        <- rep(point.col, length(windowed$bp.mid))
windowed$inv.col        <- ifelse(windowed$chr == 21 & 
                                      windowed$bp.start >= 9934148 &
                                      windowed$bp.end <= 11624150, inv.col, point.col)

###~~~~~~~~~ plot pi out ~~~~~~~~~~###

# on hybrid map, plot total pi, pi(bt), pi(rs)

outdir <- "/Users/thom/Dropbox/1.Shared_with_Bill/9.Manuscripts/2.LinkedVariation/3.figures/2.drafts/correlations/"
pdf(paste0(outdir,"recratepicorr_Hybridwindow",window.size,".pdf"),
    width = 7, height=3)
xmax <- ifelse(window.size == 100000, 140, 25)
ymax <- max(windowed$pi, na.rm=T)
layout(matrix(c(1:3), 1, 3, byrow=T))
plot(windowed$cMperMb.h, windowed$pi.bt, xlim = c(0,xmax), ylim = c(0,ymax),
     pch = 20, col = rgb(0,0,0,.5))
plot(windowed$cMperMb.h, windowed$pi.rs, xlim = c(0,xmax), ylim = c(0,ymax),
     pch = 20, col = rgb(0,0,0,.5))
plot(windowed$cMperMb.h, windowed$pi, xlim = c(0,xmax), ylim = c(0,ymax), 
     pch = 20, col = rgb(0,0,0,.5))
dev.off()

###~~~~~~~~~ plot fst out ~~~~~~~~~~###

# on hybrid map, plot total pi, pi(bt), pi(rs)

outdir <- "/Users/thom/Dropbox/1.Shared_with_Bill/9.Manuscripts/2.LinkedVariation/3.figures/2.drafts/correlations/"
pdf(paste0(outdir,"recratefstcorr_window",window.size,".pdf"),
    width = 7, height=3)
xmax <- ifelse(window.size == 100000, 140, 25)
layout(matrix(c(1:3), 1, 3, byrow=T))
plot(windowed$cMperMb.bt, windowed$fst, xlim = c(0,xmax), 
     pch = 20, col = rgb(0,0,0,.5))
plot(windowed$cMperMb.rs, windowed$fst, xlim = c(0,xmax), 
     pch = 20, col = rgb(0,0,0,.5))
plot(windowed$cMperMb.h, windowed$fst, xlim = c(0,xmax), 
     pch = 20, col = rgb(0,0,0,.5))
dev.off()



###~~~~~FST~~~~~###

par(mfrow = c(1,3))
plot(windowed$cMperMb.bt, windowed$fst, xlim = c(0,25), ylim = c(0,1), 
     pch = 20, col = windowed$inv.col, main = "BT",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "FST")
points(inv.means$bt.rate, inv.means$fst, pch = 20, col = inv.col)
print(c    <-    cor.test(windowed$cMperMb.bt, windowed$fst, method="spearman"))
plot(windowed$cMperMb.rs, windowed$fst, xlim = c(0,25), ylim = c(0,1), 
     pch = 20, col = windowed$inv.col, main = "RS",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "FST")
points(inv.means$rs.rate, inv.means$fst, pch = 20, col = inv.col)
cor.test(windowed$cMperMb.rs, windowed$fst, method="spearman")
plot(windowed$cMperMb.h, windowed$fst, xlim = c(0,25), ylim = c(0,1), 
     pch = 20, col = windowed$inv.col, main = "Hybrid",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "FST")
points(inv.means$hyb.rate, inv.means$fst, pch = 20, col = inv.col)
cor.test(windowed$cMperMb.h, windowed$fst, method="spearman")

###~~~~~DXY~~~~~###

plot(windowed$cMperMb.bt, windowed$dxy, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.bt, windowed$dxy, method="spearman")
plot(windowed$cMperMb.rs, windowed$dxy, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.rs, windowed$dxy, method="spearman")
plot(windowed$cMperMb.h, windowed$dxy, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5), main = "Hybrid",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "dXY")
cor.test(windowed$cMperMb.h, windowed$dxy, method="spearman")

###~~~~~Da~~~~~###

plot(windowed$cMperMb.bt, windowed$da, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.bt, windowed$da, method="spearman")
plot(windowed$cMperMb.rs, windowed$da, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.rs, windowed$da, method="spearman")
plot(windowed$cMperMb.h, windowed$da, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5), main = "Hybrid",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "dXY")
cor.test(windowed$cMperMb.h, windowed$da, method="spearman")

###~~~~~PI~~~~~###

plot(windowed$cMperMb.bt, windowed$pi, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.bt, windowed$pi, method="spearman")
plot(windowed$cMperMb.rs, windowed$pi, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.rs, windowed$pi, method="spearman")
plot(windowed$cMperMb.h, windowed$pi, xlim = c(0,25), 
     pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.h, windowed$pi, method="spearman")

###~~~~~PI(Boot)~~~~~###

plot(windowed$cMperMb.bt, windowed$pi.bt, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5), main = "Boot",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "pi, BT")
cor.test(windowed$cMperMb.bt, windowed$pi.bt, method="spearman")
plot(windowed$cMperMb.rs, windowed$pi.bt, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.rs, windowed$pi.bt, method="spearman")
plot(windowed$cMperMb.h, windowed$pi.bt, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5), main = "Hybrid",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "pi, BT")
cor.test(windowed$cMperMb.h, windowed$pi.bt, method="spearman")

###~~~~~PI(FW)~~~~~###

plot(windowed$cMperMb.bt, windowed$pi.fw, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5), main = "Boot",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "pi, BT")
cor.test(windowed$cMperMb.bt, windowed$pi.fw, method="spearman")
plot(windowed$cMperMb.rs, windowed$pi.fw, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.rs, windowed$pi.fw, method="spearman")
plot(windowed$cMperMb.h, windowed$pi.fw, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5), main = "Hybrid",
     xlab = paste0("cM per Mb in ",window.size/1000," kb windows"), ylab= "pi, BT")
cor.test(windowed$cMperMb.h, windowed$pi.fw, method="spearman")

###~~~~~PI(Rabbit)~~~~~###

plot(windowed$cMperMb.bt, windowed$pi.rs, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.bt, windowed$pi.rs, method="spearman")
plot(windowed$cMperMb.rs, windowed$pi.rs, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.rs, windowed$pi.rs, method="spearman")
plot(windowed$cMperMb.h, windowed$pi.rs, pch = 20, col = rgb(windowed$d.dens.scaled,0,0,.5))
cor.test(windowed$cMperMb.h, windowed$pi.rs, method="spearman")





### EVALUATE CORRELATIONS AT VARIABLE WINDOW SIZES

window.sizes    <- c(100000,  200000, 300000, 400000, 500000, 750000,
                     1000000,1500000,2000000,2500000,5000000,7500000) # in bp
rho.rate.pi    <- NULL
rho.rate.pi.rs <- NULL
rho.rate.fst   <- NULL
var.rate       <- NULL
var.rate.pi    <- NULL
var.rate.pi.rs <- NULL
var.rate.fst   <- NULL

for (wsize in window.sizes) {
    windowed    <-    NULL
    for (i in 1:21)    {
        #	i = 1
        cat(paste0("Scanning chromosome ",i,"\n"))
        m        <-    maps[maps$chr == i,]
        p        <-    popgen[popgen$chr == i,]
        m.rs    <-    m[m$map == "RS",]
        m.bt    <-    m[m$map == "Bt_Subset",]
        m.h        <-    m[m$map == "Hybrid_Subset",]
        # find the highest start position among maps to avoid "NA"s
        start.min    <-    max(m.rs$bp[1],m.bt$bp[1],m.h$bp[1], na.rm=TRUE)
        # find the lowest end position among maps to avoid "NA"s
        end.min        <-    min(m.rs$bp[length(m.rs$bp)],
                                 m.bt$bp[length(m.bt$bp)],
                                 m.h$bp[length(m.h$bp)], na.rm=TRUE)
        # identify all start positions, use only windows with full 250kb
        starts        <-    seq(start.min,end.min,by=wsize)
        n.windows    <-    length(starts)-1
        starts        <-    starts[1:n.windows]
        # record: chr, bp.start, bp.end, bp.mid, cMperMb.bt, cMperMb.rs, cMperMb.h,
        # (mean)  pi, pi.bt, pi.rs, dxy, fst
        # = 12 columns
        for (j in starts) {
            #    	j = starts[1]
            end    <-    j + (wsize - 1)
            bp.mid    <-    round((j+end)/2)
            
            cM.bt    <-    m.bt$cM[m.bt$bp >= j & m.bt$bp <= end]
            cMperMb.bt    <- ifelse(
                length(cM.bt) > 1, 
                (max(cM.bt, na.rm=T) - min(cM.bt, na.rm=T)) / (wsize / 1000000),
                ifelse(
                    is.na(cM.bt),
                    NA,
                    cM.bt)
            )
            cM.rs    <-    m.rs$cM[m.rs$bp >= j & m.rs$bp <= end]
            cMperMb.rs    <- ifelse(
                length(cM.rs) > 1, 
                (max(cM.rs, na.rm=T) - min(cM.rs, na.rm=T)) / (wsize / 1000000),
                ifelse(
                    is.na(cM.rs),
                    NA,
                    cM.rs)
            )      
            cM.h    <-    m.h$cM[m.h$bp >= j & m.h$bp <= end]
            cMperMb.h    <- ifelse(
                length(cM.h) > 1, 
                (max(cM.h, na.rm=T) - min(cM.h, na.rm=T)) / (wsize / 1000000),
                ifelse(
                    is.na(cM.h),
                    NA,
                    cM.h)
            )
            place2b    <-    p$bp >= j & p$bp <= end
            pi         <-    mean(p$pi[place2b])
            pi.rs      <-    mean(p$pi.rs[place2b])
            fst        <-    mean(p$fst.rsbt[place2b])
            result     <-    c(cMperMb.bt, cMperMb.rs, cMperMb.h,
                               pi,pi.rs, fst)
            windowed    <-    rbind(windowed, result)
        }
        
    }
    rownames(windowed)    <-    NULL
    windowed    <-    data.frame(windowed)
    names(windowed)    <-    c("cMperMb.bt","cMperMb.rs","cMperMb.h","pi","pi.rs","fst")
    rate.pi <- cor.test(windowed$cMperMb.h, windowed$pi,method="spearman")$estimate
    rate.pi.rs <- cor.test(windowed$cMperMb.h, windowed$pi.rs,method="spearman")$estimate
    rate.fst   <- cor.test(windowed$cMperMb.h, windowed$fst,method="spearman")$estimate
    rho.rate.pi    <- append(rho.rate.pi, rate.pi)
    rho.rate.pi.rs <- append(rho.rate.pi.rs, rate.pi.rs)
    rho.rate.fst   <- append(rho.rate.fst, rate.fst)
    var.rate       <- append(var.rate,       var(windowed$cMperMb.h,na.rm=T))
    var.rate.pi    <- append(var.rate.pi,    var(windowed$pi))
    var.rate.pi.rs <- append(var.rate.pi.rs, var(windowed$pi.rs))
    var.rate.fst   <- append(var.rate.fst,   var(windowed$fst))
}

outdir <- "/Users/thom/Dropbox/1.Shared_with_Bill/9.Manuscripts/2.LinkedVariation/3.figures/2.drafts/correlations/"
pdf(paste0(outdir,"HybridMapCorrByWindowSize.pdf"),
    width = 6.5, height=2.5)
    layout(matrix(c(1:4),1,4,byrow=T))
    par(mar = c(5,4,4,1))
    plot(window.sizes/1000, sqrt(var.rate), pch=20, 
         xlab = "window size, Mb", ylab = "SD", main="recombination rate",
         xlim = c(0,8000), xaxt="n")
        axis(1, at=c(0,2000,4000,6000,8000),
             labels = c(0,2,4,6,8))
        plot(window.sizes/1000, rho.rate.pi, pch=20, 
             xlab = "", ylab = "spearman's rho", main="pi",
             xlim = c(0,8000), xaxt="n")
        axis(1, at=c(0,2000,4000,6000,8000),
             labels = c(0,2,4,6,8))
        plot(window.sizes/1000, rho.rate.pi.rs, pch=20, 
         xlab = "", ylab = "", main="pi, RS",
         xlim = c(0,8000), xaxt="n")
       axis(1, at=c(0,2000,4000,6000,8000),
            labels = c(0,2,4,6,8))
       plot(window.sizes/1000, rho.rate.fst, pch=20, 
         xlab = "", ylab = "", main="FST",
         xlim = c(0,8000), xaxt="n")
       axis(1, at=c(0,2000,4000,6000,8000),
            labels = c(0,2,4,6,8))
    layout(matrix(c(1),1,1,byrow=T))
dev.off()    





###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### CREATE DISTRIBUTIONS OF PERCENT OF MAP 
##  WHERE STAT ≥ x
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

map.smoov	<-	read.table(paste0(baseDIR,"4.data/1.datasets/maps_smoothed_subset_with99point9percCIs.tsv"), 
                        sep = '\t', header = T)
stat <- map.smoov$fst
cutoffs	<-	seq(0,1,by=0.01)
#  OR
stat <- map.smoov$dxy
cutoffs	<-	seq(0,0.02,by=0.0001)

gt.percent	<-	matrix(nrow = length(cutoffs), ncol = 5)


for (c in 1:length(cutoffs))	{
    
    map.smoov$gt	<-	0
    map.smoov$gt[stat >= cutoffs[c]] <- 1
    
    gt.Mb	<-	NULL
    gt.bt	<-	NULL
    gt.rs	<-	NULL
    gt.h	<-	NULL
    
    for (i in 1:21) {
        #i = 2
        ms		<-	map.smoov[map.smoov$chr == i,]
        ints.Mb		<-	intervals(vals = ms$Mb, conditions = ms$gt)
        chr		<-	rep(i, dim(ints.Mb)[1])
        ints.Mb$chr	<-	chr
        ints.cM.bt	<-	intervals(vals = ms$cM.bt, conditions = ms$gt)
        chr		<-	rep(i, dim(ints.cM.bt)[1])
        ints.cM.bt$chr	<-	chr
        ints.cM.rs	<-	intervals(vals = ms$cM.rs, conditions = ms$gt)
        chr		<-	rep(i, dim(ints.cM.rs)[1])
        ints.cM.rs$chr	<-	chr
        ints.cM.h	<-	intervals(vals = ms$cM.h, conditions = ms$gt)
        chr		<-	rep(i, dim(ints.cM.h)[1])
        ints.cM.h$chr	<-	chr
        
        gt.Mb	<-	rbind(gt.Mb, ints.Mb)
        gt.bt	<-	rbind(gt.bt, ints.cM.bt)
        gt.rs	<-	rbind(gt.rs, ints.cM.rs)
        gt.h	<-	rbind(gt.h, ints.cM.h)
    }
    
    len.mapped.Mb	<-	sum(gt.Mb$len, na.rm = TRUE)
    #	gt.Mb$genomic.start	<-	make_genomic(gt.Mb$start*1000000,gt.Mb$chr) / 1000000
    len.fst.hi.Mb	<-	sum(gt.Mb$len[gt.Mb$signif == 'hi'], na.rm = TRUE)
    Mb.high			<-	gt.Mb[gt.Mb$signif == 'hi',]
    per.fst.hi.Mb	<-	len.fst.hi.Mb / len.mapped.Mb
    
    len.mapped.cM.bt	<-	sum(gt.bt$len, na.rm = TRUE)
    len.fst.hi.cM.bt	<-	sum(gt.bt$len[gt.bt$signif == 'hi'], na.rm = TRUE)
    per.fst.hi.cM.bt	<-	len.fst.hi.cM.bt / len.mapped.cM.bt
    
    len.mapped.cM.rs	<-	sum(gt.rs$len, na.rm = TRUE)
    len.fst.hi.cM.rs	<-	sum(gt.rs$len[gt.rs$signif == 'hi'], na.rm = TRUE)
    per.fst.hi.cM.rs	<-	len.fst.hi.cM.rs / len.mapped.cM.rs
    
    len.mapped.cM.h	<-	sum(gt.h$len, na.rm = TRUE)
    len.fst.hi.cM.h	<-	sum(gt.h$len[gt.h$signif == 'hi'], na.rm = TRUE)
    per.fst.hi.cM.h	<-	len.fst.hi.cM.h / len.mapped.cM.h
    
    result	<-	c(cutoffs[c],per.fst.hi.Mb,per.fst.hi.cM.bt,per.fst.hi.cM.rs,per.fst.hi.cM.h)
    gt.percent[c,]	<-	result
}

gt.percent	<-	data.frame(gt.percent)
names(gt.percent)	<-	c("stat.cutoff","phys","BT","RS","Hybrid")

xlim <- c(0,1)
xlim <- c(0,0.02)
xlim <- c(3,10)
xlim <- range(cutoffs)

par(mfrow = c(1,1))
plot(gt.percent$stat.cutoff, gt.percent$phys,type='l',lwd = 2, xlim = xlim,
     xlab = 'Smoothed stat cutoff', ylab = "percent of map ≥ cutoff")
lines(gt.percent$stat.cutoff, gt.percent$BT, col = 'blue3')
lines(gt.percent$stat.cutoff, gt.percent$RS, col = 'red3')
lines(gt.percent$stat.cutoff, gt.percent$Hybrid, col = 'purple3', lwd = 2)

gt.percent$BT.diff	<-	gt.percent$phys - gt.percent$BT
gt.percent$RS.diff	<-	gt.percent$phys - gt.percent$RS
gt.percent$Hybrid.diff	<-	gt.percent$phys - gt.percent$Hybrid

gt.percent$BT.percPhys	<-	gt.percent$BT / gt.percent$phys
gt.percent$RS.percPhys	<-	gt.percent$RS / gt.percent$phys
gt.percent$Hybrid.percPhys	<-	gt.percent$Hybrid / gt.percent$phys


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GENETIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SCANS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                   ###
###              FST GENETIC MAP SCANS                ###
###                                                   ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

library(RColorBrewer)
map.smoov	<-	read.table(paste0(baseDIR,"4.data/1.datasets/maps_smoothed_subset_with99point9percCIs.tsv"), 
                        sep = '\t', header = T)

map.smoov$fst.grays	<-	1- scale.vals(map.smoov$fst, 0.3,1.0)
map.smoov$fst.grays	<-	hsv(0,0, map.smoov$fst.grays)

plot.map.fst	<-	function(chr = "genomic", low.col = "lightblue", mid.col = "orange", hi.col = "red3") {
    # chr = "21"; low.col = "lightblue"; mid.col = "orange"; hi.col = "red3"
    chr	<-	chr
    if (is.numeric(chr)) {
        ms	<-	map.smoov[map.smoov$chr == chr,]
        Mb.max		<-	max(map.smoov$Mb, na.rm = TRUE)
        cM.max.bt	<-	max(map.smoov$cM.bt, na.rm = TRUE)
        cM.max.rs	<-	max(map.smoov$cM.rs, na.rm = TRUE)
        cM.max.h	<-	max(map.smoov$cM.h, na.rm = TRUE)
        cM.max		<-	max(cM.max.bt,cM.max.rs,cM.max.h)
    } else {
        ms <- map.smoov
        ms$Mb <- ms$Mb.genomic
        ms$cM.bt <- ms$bt.genomic
        ms$cM.rs <- ms$rs.genomic
        ms$cM.h  <- ms$h.genomic
        Mb.max		<-	max(ms$Mb, na.rm = TRUE)
        cM.max.bt	<-	max(ms$cM.bt, na.rm = TRUE)
        cM.max.rs	<-	max(ms$cM.rs, na.rm = TRUE)
        cM.max.h	<-	max(ms$cM.h, na.rm = TRUE)
        cM.max		<-	max(cM.max.bt,cM.max.rs,cM.max.h)
    }
    ms	<-	ms[order(ms$fst),]
    colfx <- colorRampPalette(c(low.col,mid.col,hi.col),interpolate="linear")	
    cols	<-	colfx(7)	
    ms$col	<-	rep(			cols[7], dim(ms)[1])
    ms$col[ms$fst < 0.7] <-	cols[6]
    ms$col[ms$fst < 0.6] <-	cols[5]
    ms$col[ms$fst < 0.5] <-	cols[4]
    ms$col[ms$fst < 0.4] <-	cols[3]
    ms$col[ms$fst < 0.3] <-	cols[2]
    ms$col[ms$fst < 0.2] <-	cols[1]
    Mb.axis		<-	3.00 # y position to start Mb axis
    cm.bt.top	<-	2.50 # y position to be top (bottom) of cM-based Fst scan
    cm.rs.top	<-	2.00 # y position to be top (bottom) of cM-based Fst scan
    cm.h.top	<-	1.50 # y position to be top (bottom) of cM-based Fst scan
    mapwidth	<-	0.05 ; mapoffset	<-	mapwidth / 2
    ms$fst.Mb	<-	ms$fst + Mb.axis
    ms$fst.cM	<-	cm.h.top - ms$fst
    xmax		<-	Mb.max
    Mb.med		<-	(max(ms$Mb, na.rm=T) - min(ms$Mb, na.rm=T)) / 2
    cM.bt.med	<-	(max(ms$cM.bt, na.rm=T) - min(ms$cM.bt, na.rm=T)) / 2
    cM.rs.med	<-	(max(ms$cM.rs, na.rm=T) - min(ms$cM.rs, na.rm=T)) / 2
    cM.h.med	<-	(max(ms$cM.h, na.rm=T) - min(ms$cM.h, na.rm=T)) / 2
    
    # SCALE MB SO THAT MAX(CHR MB) == MAX(LG cM)
    Mb.med.glob	<-	Mb.max / 2
    Mb.med.diff	<-	Mb.med.glob - Mb.med
    ms$Mb.scale	<-	ms$Mb + Mb.med.diff
    # SCALE BOOT CROSS CENTIMORGANS TO TOTAL MAP LENGTH
    cM.bt.med.glob	<-	cM.max.bt / 2
    cM.bt.med.diff	<-	cM.bt.med.glob - cM.bt.med
    ms$cM.bt.scale	<-	(ms$cM.bt + cM.bt.med.diff) * (xmax/cM.max.bt)
    bt.scale.bar.len<-	10 * (xmax/cM.max.bt)
    # SCALE RABBIT CROSS CENTIMORGANS TO TOTAL MAP LENGTH
    cM.rs.med.glob	<-	cM.max.rs / 2
    cM.rs.med.diff	<-	cM.rs.med.glob - cM.rs.med
    ms$cM.rs.scale	<-	(ms$cM.rs + cM.rs.med.diff) * (xmax/cM.max.rs)
    rs.scale.bar.len<-	10 * (xmax/cM.max.rs)
    # SCALE HYBRID CROSS CENTIMORGANS TO TOTAL MAP LENGTH	
    cM.h.med.glob	<-	cM.max.h / 2
    cM.h.med.diff	<-	cM.h.med.glob - cM.h.med
    ms$cM.h.scale	<-	(ms$cM.h + cM.h.med.diff) * (xmax/cM.max.h)
    h.scale.bar.len	<-	10 * (xmax/cM.max.h)
    cols	<-	ms$col
    plot(0,0, xlim = c(0,xmax), ylim = c(0,Mb.axis+1), xlab = '', ylab = '',
         xaxt = 'n', yaxt = 'n', type = 'n')
    axis1	<-	rep(Mb.axis, length(ms$Mb))
    axis2	<-	rep(cm.bt.top, length(ms$Mb))
    axis3	<-	rep(cm.rs.top, length(ms$Mb))
    axis4	<-	rep(cm.h.top, length(ms$Mb))
    # TICK MARKS AND ALL AXES
    segments(x0 = 0, x1 = -0.5, 
             y0 = seq(cm.h.top-1,cm.h.top,by=0.1), y1 = seq(cm.h.top-1,cm.h.top,by=0.1))
    segments(x0 = 0, x1 = -0.5, 
             y0 = seq(Mb.axis,Mb.axis+1,by=0.1), y1 = seq(Mb.axis,Mb.axis+1,by=0.1))
    segments(x0 = 0, x1 = xmax, y0 = cm.h.top-1, y1 = cm.h.top-1)
    segments(x0 = 0, x1 = xmax, y0 = cm.h.top, y1 = cm.h.top)
    segments(x0 = 0, x1 = xmax, y0 = cm.bt.top, y1 = cm.bt.top)
    segments(x0 = 0, x1 = xmax, y0 = cm.rs.top, y1 = cm.rs.top)
    segments(x0 = 0, x1 = xmax, y0 = Mb.axis, y1 = Mb.axis)
    segments(x0 = 0, x1 = xmax, y0 = Mb.axis+1, y1 = Mb.axis+1)		
    segments(x0 = 0, x1 = 0, y0 = cm.h.top-1, y1 = Mb.axis+1)
    segments(x0 = xmax, x1 = xmax, y0 = cm.h.top-1, y1 = Mb.axis+1)
    # write color key
    segments(x0 = xmax - (xmax/50), x1 = xmax, y0 = seq(Mb.axis+0.1,Mb.axis+0.9,length.out=7),col = colfx(7))
    # CHROMOSOME BREAKS IF 'GENOMIC'
    if (!(is.numeric(chr))) {chr.breaks(ybottom=Mb.axis, ytop=Mb.axis+1)}
    # ADD SEGMENTS CONNECTING MAPS
    segments(x0 = ms$cM.bt.scale, x1 = ms$Mb.scale, y0 = cm.bt.top+mapoffset
             ,y1=Mb.axis, col = adjustcolor(cols, alpha=0.25), lwd = 2)
    segments(x0 = ms$cM.rs.scale, x1 = ms$cM.bt.scale, y0 = cm.rs.top+mapoffset,
             y1=cm.bt.top-mapoffset,col = adjustcolor(cols, alpha=0.25), lwd = 2)
    segments(x0 = ms$cM.h.scale, x1 = ms$cM.rs.scale, y0 = cm.h.top,
             y1=cm.rs.top-mapoffset, col = adjustcolor(cols, alpha=0.25), lwd = 2)
    # ADD SEGMENTS FOR BT, RS MAP MARKERS
    segments(x0 = ms$cM.bt.scale, x1 = ms$cM.bt.scale, 
             y0 = cm.bt.top-mapoffset,y1=cm.bt.top+mapoffset,
             col = adjustcolor(cols), lwd = 3)
    segments(x0 = ms$cM.rs.scale, x1 = ms$cM.rs.scale, 
             y0 = cm.rs.top-mapoffset,y1=cm.rs.top+mapoffset,
             col = adjustcolor(cols), lwd = 3)
    # ADD SEGMENTS FOR FST HEIGHT UNDER CURVES
    segments(x0 = ms$cM.h.scale, x1 = ms$cM.h.scale, y0 = ms$fst.cM,y1=cm.h.top,
             col = cols)
    segments(x0 = ms$Mb.scale, x1 = ms$Mb.scale, y0 = ms$fst.Mb,y1=Mb.axis,
             col = cols)
    # ADD CURVES
    ms	<-	ms[order(ms$Mb.scale),]
    lines(x = ms$Mb.scale, y = ms$fst.Mb)
    lines(x = ms$cM.h.scale, y = ms$fst.cM)
    # ADD SCALE BARS
    segments(x0 = 0.5, x1 = 5.5, y0 = Mb.axis+0.5, y1 = Mb.axis+0.5,
             lwd = 3)
    text(x = 5.5, y = Mb.axis+0.5, pos = 3, labels = "5 Mb")
    segments(x0 = 0.5, x1 = bt.scale.bar.len+0.5, 
             y0 = cm.bt.top+0.25, y1 = cm.bt.top+0.25,
             lwd = 3)
    text(x = bt.scale.bar.len+0.5, y = cm.bt.top+0.25, pos = 3, labels = "10 cM")
    segments(x0 = 0.5, x1 = rs.scale.bar.len+0.5, 
             y0 = cm.rs.top+0.25, y1 = cm.rs.top+0.25,
             lwd = 3)
    text(x = rs.scale.bar.len+0.5, y = cm.rs.top+0.25, pos = 3, labels = "10 cM")
    segments(x0 = 0.5, x1 = h.scale.bar.len+0.5, 
             y0 = cm.h.top-0.5, y1 = cm.h.top-0.5,
             lwd = 3)
    text(x = h.scale.bar.len+0.5, y = cm.h.top-0.25, pos = 1, labels = "10 cM")
    
}

plot.map.fst(chr = 7, low.col = "gray80", mid.col = "gray50", hi.col = "gray10")
plot.map.fst(chr = 'genomic')

# Summaries of chr 7 single locus

summary(map.smoov$fst[map.smoov$chr==7 & map.smoov$cM.h == 57.665])
range(map.smoov$Mb[map.smoov$chr==7 & map.smoov$cM.h == 57.665], na.rm=T)

# print all chromosomes out to file
outdir <- paste0(baseDIR,"3.figures/2.drafts/scans/chromosomes/")

for (i in 1:21) {
    pdf(file = paste0(outdir,"fst_fullfams_chromogroup",i,".pdf"), width=5,height=5)
    plot.map.fst(chr = i, colmode = "NOT")
    dev.off()
}

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                   DISTANCE DECAY CURVES                  ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

popgen        <-    popgen[popgen$cM.bt.interp != -1 & 
                               popgen$cM.rs.interp != -1 &
                               popgen$cM.h.interp  != -1 ,]

sorted		<-	popgen[popgen$cons.sort != 'N' & 
                      is.na(popgen$pi) == FALSE &
                      is.na(popgen$pi.fw) == FALSE &
                      is.na(popgen$pi.rs) == FALSE &
                      is.na(popgen$dist2d) == FALSE &
                      popgen$chr != 19,]
sorted.neu	<-	sorted[sorted$cons.sort != 'D',]
sorted.neu	<-	sorted.neu[sorted.neu$dist2d < 50000,]

# BY PHYSICAL DISTANCE TO A DIVERGENT LOCUS

d			<-	sorted[sorted$cons.sort == "D",]
n.d			<-	dim(d)[1]
pi      	<-	mean(d$pi, na.rm=T)
pi.bt		<-	mean(d$pi.bt, na.rm=T)
pi.bp		<-	mean(d$pi.bp, na.rm=T)
pi.fw		<-	mean(d$pi.fw, na.rm=T)
pi.rs		<-	mean(d$pi.rs, na.rm=T)
dxy			<-	mean(d$dxy.rsbt, na.rm=T)
d.mean		<-	data.frame(pi,pi.bt,pi.rs,dxy,tmrca)

pi.cor	    <-	cor.test(sorted.neu$dist2d, sorted.neu$pi, method = 'spearman')
bt.cor	    <-	cor.test(sorted.neu$dist2d, sorted.neu$pi.bt, method = 'spearman')
bp.cor	    <-	cor.test(sorted.neu$dist2d, sorted.neu$pi.bp, method = 'spearman')
fw.cor	    <-	cor.test(sorted.neu$dist2d, sorted.neu$pi.fw, method = 'spearman')
rs.cor	    <-	cor.test(sorted.neu$dist2d, sorted.neu$pi.rs, method = 'spearman')
dxy.cor	    <-	cor.test(sorted.neu$dist2d, sorted.neu$dxy.rsbt, method = 'spearman')

# fraction of RAD loci w/i 250kb of divergent locus
#  with pi(rs) == 0.

nRS0 <- length(popgen$locus[popgen$dist2d <= 250 & popgen$pi.rs == 0])
ntot <- length(popgen$locus[popgen$dist2d <= 250])

nRS0/ntot



# PLOT PI

xmax	<-	2000
ymax	<-	0.020
xbins	<-	xmax / 50
shape	<-	2

x			<-	sorted.neu$dist2d
y			<-	sorted.neu$pi
data		<-	data.frame(x,y)
data.lim	<-	data[data$x <= xmax & data$y <= ymax,]
h			<-	hexbin(data.lim, xbins = xbins, shape = shape)
gplot.hexbin(h, xlab = 'Distance to divergent locus, kb', ylab = "pi")
gplot.hexbin(h, xlab = 'Distance to divergent locus, kb', ylab = "pi",
             maxcnt=200, colorcut = 15)

# SPLINES W/ PHYSICAL DISTANCE

splinemethod	<-	"natural"
df				<-	40
plot(smooth.spline(sorted.neu$dist2d, sorted.neu$pi, df = df), 
     xlim = c(0,xmax), ylim = c(0,0.006), xaxt = "n", yaxt = "n",
     type = 'l', xlab = "", ylab = "", lwd = 2)
axis(1, at = seq(0,2000,by=500))
axis(2, at = seq(0,0.006,by=0.001))
abline(h = mean(sorted.neu$pi), lty = 'dotted')
abline(h = mean(sorted.neu$pi.bt), lty = 'dotted', col = 'blue3')
abline(h = mean(sorted.neu$pi.rs), lty = 'dotted', col = 'red3')
lines(smooth.spline(sorted.neu$dist2d, sorted.neu$pi.rs, df = df), 
      col = 'red3', lwd = 2)
lines(smooth.spline(sorted.neu$dist2d, sorted.neu$pi.bt, df = df), 
      col = 'blue3')
mtext("Distance to divergent locus, kb", side = 1, line = 3)
mtext("pi per site", side = 2, line = 3)



### Interaction plot of pi by position by pop

rs.linked   <- sorted$pi.rs[sorted$dist2d < 1000 & sorted$cons.sort != 'D']
    rs.l.mean <- mean(rs.linked)
    rs.l.sd   <- sd(rs.linked)
rs.unlinked <- sorted$pi.rs[sorted$dist2d == 50000]
    rs.u.mean <- mean(rs.unlinked)
    rs.u.sd   <- sd(rs.unlinked)
fw.linked   <- sorted$pi.fw[sorted$dist2d < 1000 & sorted$cons.sort != 'D']
    fw.l.mean <- mean(fw.linked)
    fw.l.sd   <- sd(fw.linked)
fw.unlinked <- sorted$pi.fw[sorted$dist2d == 50000]
    fw.u.mean <- mean(fw.unlinked)
    fw.u.sd   <- sd(fw.unlinked)

rs.offset <- -0.01
fw.offset <- 0.01
plot(0,0, xlim = c(0.9,2.1), ylim = c(0.001,0.005), type = "n",
     xaxt = "n", xlab = "linkage", ylab = "pi", las = 1)
    segments(x0 = 1+rs.offset, x1 = 1+rs.offset,
             y0 = rs.l.mean - (0.5*rs.l.sd),
             y1 = rs.l.mean + (0.5*rs.l.sd),
             col = "red3")
    segments(x0 = 2+rs.offset, x1 = 2+rs.offset,
             y0 = rs.u.mean - (0.5*rs.u.sd),
             y1 = rs.u.mean + (0.5*rs.u.sd),
             col = "red3")
    segments(x0 = 1+fw.offset, x1 = 1+fw.offset,
             y0 = fw.l.mean - (0.5*fw.l.sd),
             y1 = fw.l.mean + (0.5*fw.l.sd),
             col = "blue")
    segments(x0 = 2+fw.offset, x1 = 2+fw.offset,
             y0 = fw.u.mean - (0.5*fw.u.sd),
             y1 = fw.u.mean + (0.5*fw.u.sd),
             col = "blue")
    segments(x0 = 1+fw.offset, x1 = 2+fw.offset,
             y0 = fw.l.mean,y1 = fw.u.mean,col = "blue")
    segments(x0 = 1+rs.offset, x1 = 2+rs.offset,
             y0 = rs.l.mean,y1 = rs.u.mean,col = "red3")
    points(x = c(1+rs.offset,2+rs.offset), y = c(rs.l.mean,rs.u.mean), pch=20, col="red3")
    points(x = c(1+fw.offset,2+fw.offset), y = c(fw.l.mean,fw.u.mean), pch=20, col="blue")

# perform ANOVA on pi and log-transformed pi to test for interaction    

rslab <- rep("RS",length(rs.linked)+length(rs.unlinked))
fwlab <- rep("FW",length(fw.linked)+length(fw.unlinked))
habitat <- c(rslab,fwlab)
linkage <- c(rep("linked",length(rs.linked)),
               rep("unlinked",length(rs.unlinked)),
               rep("linked",length(fw.linked)),
               rep("unlinked",length(fw.unlinked)))
pi    <- c(rs.linked,rs.unlinked,fw.linked,fw.unlinked)
logpi <- log(c(rs.linked,rs.unlinked,fw.linked,fw.unlinked)+0.0001)

pibyhab <- data.frame(habitat,linkage,pi,logpi)

summary(pibyhab.anova <- aov(pi ~ habitat + linkage + habitat*linkage, pibyhab))
    TukeyHSD(pibyhab.anova)
summary(logpibyhab.anova <- aov(logpi ~ habitat + linkage + habitat*linkage, pibyhab))
    TukeyHSD(logpibyhab.anova)





###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                    FORWARD SIMULATIONS                   ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

simdir <- "/Users/thom/Dropbox/1.Oregon/Data/active/simulations/forward/SLiM/simulations/migrationSelectionBalance/"

ks   <- c(1,2,5)
scoeffs <- c("20","05")
s    <- "20"
ms   <- c("02","1","2")
ms.numeric <- c(0.2,1,2)
m    <- "1"

# Create DF to hold pi ~50cM from adaptive locus
#  from all sims
unlinkedvar <- NULL

standardize <- FALSE
ymin <- ifelse(standardize == TRUE,-6,0)
ymax <- ifelse(standardize == TRUE,4,10)
layout(matrix(1:9,3,3,byrow=TRUE))
for (i in 1:length(ms)) {
    m <- ms[i]
    m.numeric <- ms.numeric[i]
    for (k in ks) {
        # CALCULATE PI
        pi.mat <- NULL
        rescale <- 1000
        # construct string to get to right files
        paramstring <- paste0("N1-1000_N2-1000/s",s,"/m",m,"per/")
        
        simstatfiles <- list.files(paste0(simdir,paramstring))
        simstatfiles <- simstatfiles[grepl(paste0("_k",k,"_"),simstatfiles)]
        reps <- length(simstatfiles)
        # read in first file to get dimensions
        rep1 <- read.table(paste0(simdir,paramstring,simstatfiles[1]),
                           sep = "\t",header=TRUE)
        n.windows <- dim(rep1)[1]
        positions <- rep1$position.mid
        
        pi.mat <- matrix(nrow=length(positions),ncol=reps)
        
        # ADD PI ESTIMATES FOR EACH WINDOW FOR EACH SIM
        #  WITH THIS PARAMETER COMBINATION
        pi.mat[,1] <- rep1$pi.window
        for (rep in 2:reps) {
            sim <- read.table(paste0(
                simdir,paramstring,simstatfiles[rep]),
                sep="\t",header=TRUE
            )
            pi.mat[,rep] <- sim$pi.window
        }
        
        pi.grandmean <- ifelse(standardize==TRUE, mean(pi.mat), 0)
        pi.mat <- (pi.mat - pi.grandmean)*rescale
        # SUMMARIZE pi within allelic classes
        
        pi.0          <- matrix(nrow=length(positions),ncol=reps)
        pi.p1.class1  <- matrix(nrow=length(positions),ncol=reps)
        pi.1          <- matrix(nrow=length(positions),ncol=reps)
        
        for (rep in 1:length(simstatfiles)) {
            sim <- read.table(paste0(
                simdir,paramstring,simstatfiles[rep]),
                sep="\t",header=TRUE
            )
            pi.0[,rep]  <- sim$class0.pi.window
            pi.1[,rep]  <- sim$class1.pi.window
            pi.p1.class1[,rep]  <- sim$p1.class1.pi.window
        }
        # construct DF of last position pi estimates
        
        kparam <- rep(k, reps)
        mparam <- rep(m.numeric, reps)
        pi.50 <- pi.mat[length(positions),]
        pi.0.50 <- pi.0[length(positions),]
        pi.p1c1.50 <- pi.p1.class1[length(positions),]
        x <- cbind(mparam, kparam, pi.50, pi.0.50,pi.p1c1.50)
        unlinkedvar <- rbind(unlinkedvar,x)
        
        positions <- (positions / 1000) - 50
        pi.0          <- (pi.0 - pi.grandmean) * rescale
        pi.p1.class1  <- (pi.p1.class1 - pi.grandmean) * rescale
        pi.1          <- (pi.1 - pi.grandmean) * rescale
        pi.0.mean     <- mean(pi.0)
        pi.p1c1.mean  <- mean(pi.p1.class1)
        pi.mean       <- mean(pi.mat)
        # compute t-distribution-based CIs
        
        pi.ci          <- matrix(nrow=length(positions),ncol=2)
        pi.0.ci          <- matrix(nrow=length(positions),ncol=2)
        pi.p1.class1.ci  <- matrix(nrow=length(positions),ncol=2)
        pi.1.ci          <- matrix(nrow=length(positions),ncol=2)
        
        for (i in 1:length(positions)) {
            pi.ci[i,]            <- t.test(pi.mat[i,])$conf.int
            pi.0.ci[i,]          <- t.test(pi.0[i,])$conf.int
            pi.1.ci[i,]          <- t.test(pi.1[i,])$conf.int
            pi.p1.class1.ci[i,]  <- t.test(pi.p1.class1[i,])$conf.int
        }
        
        pos.dnstream  <- positions      [positions > 0 ]
        pi.dnstream   <- pi.ci          [positions > 0,]
        p0.dnstream   <- pi.0.ci        [positions > 0,]
        p1c1.dnstream <- pi.p1.class1.ci[positions > 0,]
        
        # PLOT RESULTS
        
        par(las=1, mar=c(5,5,2,1))
        plot(0,0,xlim=c(min(pos.dnstream),max(pos.dnstream)),
             ylim=c(ymin,ymax),xlab="",
             ylab="pi",type='n')
        abline(h = pi.mean)
        abline(h = pi.0.mean, col='firebrick3')
        abline(h = pi.p1c1.mean, col='deepskyblue3')
        polygon(x = c(pos.dnstream,rev(pos.dnstream)),
                y = c(pi.dnstream[,1], rev(pi.dnstream[,2])),
                col = 'gray75')
        lines(pos.dnstream, pi.dnstream[,1], lwd=0.5)
        lines(pos.dnstream, pi.dnstream[,2], lwd=0.5)
        polygon(x = c(pos.dnstream,rev(pos.dnstream)),
                y = c(p0.dnstream[,1], rev(p0.dnstream[,2])),
                col = 'firebrick3')
        lines(pos.dnstream, p0.dnstream[,1], lwd=0.5)
        lines(pos.dnstream, p0.dnstream[,2], lwd=0.5)
        polygon(x = c(pos.dnstream,rev(pos.dnstream)),
                y = c(p1c1.dnstream[,1], rev(p1c1.dnstream[,2])),
                col = 'deepskyblue3')
        lines(pos.dnstream, p1c1.dnstream[,1], lwd=0.5)
        lines(pos.dnstream, p1c1.dnstream[,2], lwd=0.5)
    }
}

# STATS ON UNLINKED DIVERSITY

unlinkedvar <- data.frame(unlinkedvar)
unlinkedvar$mparam <- factor(unlinkedvar$mparam)
unlinkedvar$kparam <- factor(unlinkedvar$kparam)

summary(pi.aov <- aov(pi.50 ~ mparam + kparam + mparam*kparam,unlinkedvar))
TukeyHSD(pi.aov)
