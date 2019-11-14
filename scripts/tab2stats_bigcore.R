#!/usr/bin/env Rscript


# load functions

heterozygosity <- function(alleles = c(1,1,1,2,2,2)) {
    # calculate heterozygosity as per Hahn 2018 eq 3.1
    # get sample n and unique alleles
    alleles <- as.character(alleles)
    n = length(alleles)
    As = unique(alleles)
    Ehomos = NULL
    for (i in 1:length(As)) {
        Ehomo = (length(alleles[alleles == As[i]])/n) ^2
        Ehomos <- append(Ehomos, Ehomo)
    }
    return((n / (n-1)) * (1-sum(Ehomos)))
}

DAF <- function(alleles = c(1,1,1,0,0,0)) {
    # calculate derived allele frequency
    return(length(alleles[alleles == 1]) / length(alleles))
}

fst.biallelic <- function(p1.alleles = c(1,1,1,0,0,0), p2.alleles = c(1,1,1,1,1,1)) {
    np1     <- length(p1.alleles)
    np2     <- length(p2.alleles)
    n.tot   <- np1 + np2
    hetp1   <- heterozygosity(p1.alleles)
    hetp2   <- heterozygosity(p2.alleles)
    het.hat <- heterozygosity(c(p1.alleles,p2.alleles))
    fst <- 1 - (((np1/n.tot)*hetp1+(np2/n.tot)*hetp2) / het.hat)
    return(fst)
}

dxy.biallelic <- function(p1.alleles = c(1,1,1,0,0,0), p2.alleles = c(1,1,1,1,1,1)) {
    # dxy is the mean pairwise diffs w/ only cross-pop comparisons
    diffs <- NULL
    for (i in p1.alleles) {
        for (j in p2.alleles) {
            diff = ifelse(i == j, 0, 1)
            diffs <- append(diffs, diff)
        }
    }
    return(sum(diffs) / length(diffs))
}

harmonic_number <- function(n) {
    return(sum(1 / (1:(n-1))))
}


args  <- commandArgs(trailingOnly=TRUE)
MStab     <- read.table(args[1],header=T)
chrlen    <- as.numeric(args[2])
binsize   <- as.numeric(args[3])
stepsize  <- as.numeric(args[4])
sampsize  <- as.numeric(args[5])
replicate <- as.numeric(args[6])
# temp for testing:
# MStab    <- read.table("/Users/thom/Dropbox/1.Oregon/Data/active/simulations/forward/SLiM/migsel_v2/output_bigcore/k1_mred1_s0.2_replicate13_end_H1all.ms.tsv", header=T)
chrlen   <- 102003
# binsize  <- 500
# stepsize <- 500
# sampsize <- 10
# replicate<- 1

### Define regions of chromosome
migselMuPos <- 0.4999951
region1     <- c(1,50000)
core        <- c(50001,52002)
region2     <- c(52003,102002)

# CLIP OFF FIRST/LAST POSITION BECAUSE THEY'RE JUST PLACEHOLDERS
nsites <- dim(MStab)[1]
MStab <- MStab[ 2:(nsites-1) , ]
positions <- as.numeric(MStab$position) * chrlen
allelicClass <- MStab[abs(MStab$position - migselMuPos) == min(abs(MStab$position - migselMuPos)),]

# subset data by allelic class
# Downsample to chromosome sample size

class0    <- MStab[,allelicClass == 0]
    class0 <- class0[, sample.int(sampsize)]
class1    <- MStab[,allelicClass == 1]
    class1 <- class1[, sample.int(sampsize)]
all       <- cbind(class0[, 1:(sampsize/2)],class1[, 1:(sampsize/2)])

pi        <- apply(all, MARGIN=1, FUN=heterozygosity)
class0.pi <- apply(class0,MARGIN=1, FUN=heterozygosity)
class1.pi <- apply(class1,MARGIN=1, FUN=heterozygosity)
pi.A      <- apply(cbind(class0.pi,class1.pi),MARGIN=1, FUN=mean)
pi.TA     <- pi - pi.A

    
    bw <- 1000
    layout(matrix(c(1,2),2,1, byrow=T))
    par(mar = c(4,4,0.5,1),las=1)
    plot(ksmooth(positions, pi, bandwidth=bw), type = 'l',
         lwd = 2, ylim = c(0,0.5), ylab='pi',xlab='')
    abline(v=chrlen/2, lty='dotted')
    lines(ksmooth(positions, class1.pi, bandwidth=bw),col='blue3')
    lines(ksmooth(positions, class0.pi, bandwidth=bw),col='firebrick2')
    plot(ksmooth(positions, pi.TA, bandwidth=bw), type = 'l',
         lwd = 2, ylim = c(0,0.5), ylab='pi_T-A',xlab='position')
    abline(v=chrlen/2, lty='dotted')
    layout(matrix(c(1,1),1,1, byrow=T))

##############################################
###   calculate stats in sliding windows   ###
### positions 1-50000 and 50002-100001 are ###
### variable, so I'll break the chromosome ###
### into two halves, scanning through each ###
##############################################

class0.pos <- positions
class1.pos <- positions

position.mid     <- NULL
pi.window        <- NULL
theta            <- NULL
class0.pi.window <- NULL
class1.pi.window <- NULL
class0.theta     <- NULL
class1.theta     <- NULL
pi.TA.window     <- NULL
dist2H2star      <- NULL

starts <- seq(region1[1],region1[2],by=stepsize)
ends   <- starts+(binsize-1)

# SCAN FIRST HALF OF CHROMOSOME
for (p in 1:length(starts)) {
    start <- starts[p]
    end   <- ends[p]
    window.mid <- (start + end) / 2
    position.mid <- append(position.mid,window.mid)
    dist2H2star  <- append(dist2H2star, (region1[2]+1) - window.mid)
    # start with overall diversity
    index     <- positions >= start & positions < end
    pi.all    <- pi[index]
    segsites  <- length(pi.all[pi.all > 0])
    pi.all    <- ifelse(segsites == 0, 0, sum(pi.all) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    pi.window <- append(pi.window,pi.all)
    theta     <- append(theta, tW)
    # now for each class
    index     <- class0.pos >= start & class0.pos < end
    pi.0      <- class0.pi[index]
    segsites  <- length(pi.0[pi.0 > 0])
    pi.0      <- ifelse(segsites == 0, 0, sum(pi.0) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    class0.pi.window <- append(class0.pi.window, pi.0)
    class0.theta     <- append(class0.theta, tW)

    index     <- class1.pos >= start & class1.pos < end
    pi.1      <- class1.pi[index]
    segsites  <- length(pi.1[pi.1 > 0])
    pi.1      <- ifelse(segsites == 0, 0, sum(pi.1) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    class1.pi.window <- append(class1.pi.window, pi.1)
    class1.theta     <- append(class1.theta, tW)
    
    pi.ta     <- pi.all - mean(c(pi.0,pi.1))
    pi.TA.window <- append(pi.TA.window, pi.ta)
}

# ADD IN PI FOR H2star

    start <- core[1]
    end   <- core[2]
    window.mid <- (start + end) / 2
    position.mid <- append(position.mid,window.mid)
    dist2H2star  <- append(dist2H2star, 0)
    # start with overall diversity
    index     <- positions >= start & positions < end
    pi.all    <- pi[index]
    segsites  <- length(pi.all[pi.all > 0])
    pi.all    <- ifelse(segsites == 0, 0, sum(pi.all) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    pi.window <- append(pi.window,pi.all)
    theta     <- append(theta, tW)
    # now for each class
    index     <- class0.pos >= start & class0.pos < end
    pi.0      <- class0.pi[index]
    segsites  <- length(pi.0[pi.0 > 0])
    pi.0      <- ifelse(segsites == 0, 0, sum(pi.0) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    class0.pi.window <- append(class0.pi.window, pi.0)
    class0.theta     <- append(class0.theta, tW)
    
    index     <- class1.pos >= start & class1.pos < end
    pi.1      <- class1.pi[index]
    segsites  <- length(pi.1[pi.1 > 0])
    pi.1      <- ifelse(segsites == 0, 0, sum(pi.1) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    class1.pi.window <- append(class1.pi.window, pi.1)
    class1.theta     <- append(class1.theta, tW)
    
    pi.ta     <- pi.all - mean(c(pi.0,pi.1))
    pi.TA.window <- append(pi.TA.window, pi.ta)


# SCAN SECOND HALF OF CHROMOSOME

starts <- seq(region2[1],region2[2],by=stepsize)
ends   <- starts+(binsize-1)

# SCAN FIRST HALF OF CHROMOSOME
for (p in 1:length(starts)) {
    start <- starts[p]
    end   <- ends[p]
    window.mid <- (start + end) / 2
    position.mid <- append(position.mid,window.mid)
    dist2H2star  <- append(dist2H2star, window.mid - (region2[1]-1))
    # start with overall diversity
    index     <- positions >= start & positions < end
    pi.all    <- pi[index]
    segsites  <- length(pi.all[pi.all > 0])
    pi.all    <- ifelse(segsites == 0, 0, sum(pi.all) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    pi.window <- append(pi.window,pi.all)
    theta     <- append(theta, tW)
    # now for each class
    index     <- class0.pos >= start & class0.pos < end
    pi.0      <- class0.pi[index]
    segsites  <- length(pi.0[pi.0 > 0])
    pi.0      <- ifelse(segsites == 0, 0, sum(pi.0) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    class0.pi.window <- append(class0.pi.window, pi.0)
    class0.theta     <- append(class0.theta, tW)
    
    index     <- class1.pos >= start & class1.pos < end
    pi.1      <- class1.pi[index]
    segsites  <- length(pi.1[pi.1 > 0])
    pi.1      <- ifelse(segsites == 0, 0, sum(pi.1) / binsize)
    tW        <- ifelse(segsites == 0, 0, (segsites/harmonic_number(sampsize))/binsize)
    class1.pi.window <- append(class1.pi.window, pi.1)
    class1.theta     <- append(class1.theta, tW)
    
    pi.ta     <- pi.all - mean(c(pi.0,pi.1))
    pi.TA.window <- append(pi.TA.window, pi.ta)
}

outmat <- data.frame(position.mid, dist2H2star,
                pi.window,theta,pi.TA.window,
                class0.pi.window, class1.pi.window, 
                class0.theta, class1.theta
                )


# outmat[, 2:length(outmat[1,])] <- round(outmat[, 2:length(outmat[1,])], 5)
# 
# plot(outmat[,1],outmat[,8], type='l')
outfile <- paste0(args[1],".stats")
write.table(outmat, outfile, sep = "\t",quote=F, col.names = TRUE, row.names=F)

