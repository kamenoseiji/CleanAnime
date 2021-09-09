library(av) # Library to generate mp4 movie

KERNEL_BL = seq(0,15)*seq(16)/2 + 1
Ant2Bl <- function(ant1, ant2){	    # Antenna -> baseline index (without autocorr)
    antenna1 = max(ant1, ant2) - 1
    antenna2 = min(ant1, ant2)
    return(antenna1* (antenna1 - 1)/2 + antenna2)
}

Bl2Ant <- function(bl_index){     # Baseline -> antenna indexing (canonical ordering)
    ant1 = max(which(KERNEL_BL<= bl_index)) + 1
    return( c(ant1, bl_index - KERNEL_BL[ant1 - 1] + 1))
}

baselineLength <- function(antennaPosition){
    antNum <- length(antennaPosition)
    blNum <- antNum* (antNum - 1)/2
    BLlength <- numeric(blNum)
    for(bl_index in 1:blNum){
        antPair <- Bl2Ant(bl_index)
        BLlength[bl_index] <- abs( antennaPosition[antPair[1]] - antennaPosition[antPair[2]] )
    }
    return(BLlength)
}

plotCC <- function(pngFile, image, ymin, ymax, linecolor='black', xline=T){
    npix <- length(image)
    png(filename = pngFile, width = 960, height = 540)
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), type = "n", col = "white", xlim = c(0, npix), ylim = c(ymin, ymax), yaxt = "n", ann = FALSE, xaxt = "n", bty = "n")
    if(xline){ abline(h=0, col='gray') }
    plotX <- seq(1:npix) - 0.5
    lines( plotX, image, type='s', col = linecolor)
    dev.off()
}

Normalize <- function(image){
    return(image/max(image))
}

GaussModel <- function(amp, size=1, pos=0, npix=1024){
    x <- (seq(-npix/2+1, npix/2)/npix - pos )/size 
    return( sqrt(2*pi)* amp* dnorm(x) )
}

Clean <- function(dirtyImage, BLweight, cleanedVisibility, gain=0.01, NF=1){
    cc <- numeric(length(dirtyImage))
    cc_index <- which.max(abs(dirtyImage))
    cc[ cc_index ] <- dirtyImage[cc_index]* gain
    lines( c(cc_index, cc_index), c(0, cc[cc_index]/NF ))
    ccVis <- fft(cc, inverse=T)* BLweight
    return(cleanedVisibility - ccVis)
}

PNGsnapShot <- function(pngFile, xList, imageList, ymin, ymax, lineColor, CC=numeric(0), xline=T){
    png(filename = pngFile, width = 960, height = 540)
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), type = "n", col = "white", xlim = c(0, length(imageList[[1]])), ylim = c(ymin, ymax), yaxt = "n", ann = FALSE, xaxt = "n", bty = "n")
    if(xline){ abline(h=0, col='gray') }
    plotNum <- length(imageList)
    for(index in 1:plotNum){
        lines(xList[[index]], imageList[[index]], type='s', col = lineColor[index])
    }
    npix <- length(CC)
    if(npix > 0){
        plotX <- seq(1:npix) - 0.5
        lines( plotX, CC, type='s', col = 'black')
    }
    dev.off()
}

#-------- Antenna configuration
antNum <- 10
Umax <- 16
AntennaPosition <- runif(antNum)
BLlength <- round(baselineLength(AntennaPosition)* Umax)
blNum <- length(BLlength)

#-------- 1-dimension baseline distribution
BLweight <- numeric(Umax* 64)
for(bl_index in 1:blNum){
    BL <-  BLlength[bl_index]
    BLweight[ BL + 2 ] = BLweight[ BL + 2 ] + 1
    BLweight[ length(BLweight) - BL] = BLweight[ length(BLweight) - BL ] + 1
}

#-------- 1-dimension synthesized beam
Beam <- fft(BLweight)
Beam <- Re(Beam)
Beam <- Normalize( c(Beam[(Umax*32+1):(Umax*64)], Beam[1:(Umax*32)]) )
PNGsnapShot("CleanImage/synthesizedBeam.png", list(1:length(Beam)), list(Beam), min(Beam), max(Beam), c('red'))

#-------- True image model
imageModel <- Normalize( GaussModel(2, 0.01, -0.2) + GaussModel(0.7, 0.03, -0.17) + GaussModel(0.5, 0.05, -0.1) +  GaussModel(0.2, 0.2, 0.0) +  GaussModel(0.3, 0.02, 0.2) )
PNGsnapShot("CleanImage/TrueImageModel.png", list(1:length(imageModel)), list(imageModel), min(imageModel), max(imageModel), c('blue'))


Visibility <- fft(imageModel, inverse=T)* BLweight / sum(BLweight)
dirtyImage <- Re(fft(Visibility))/sqrt(length(Visibility)/2)
ResidualVisibility <- Visibility
ResidualImage <- dirtyImage

PNGsnapShot("ResidualImage/Residual0000.png", list(1:length(dirtyImage)), list(dirtyImage), min(dirtyImage), max(dirtyImage), c('blue'))

#-------- Clean Loop
MaxIter <- 64
cleanGain <- 0.75
Npix <- length(dirtyImage)
cc_flux <- numeric(Npix)

for(loop_index in 1:MaxIter){
    #-------- Pick up clean component
    cc_index <- which.max(abs(ResidualImage))
    cc_flux[cc_index] <- cc_flux[cc_index] + ResidualImage[cc_index] * cleanGain
    #-------- plot cc
    plotCC(sprintf('CleanComponent/CC%04d.png', loop_index), cc_flux, min(dirtyImage), max(dirtyImage))
    PNGsnapShot(sprintf('ResidualImage/Residual%04d.png', 3*loop_index-2), list(1:Npix, 1:Npix), list(dirtyImage, ResidualImage), min(dirtyImage), max(dirtyImage), c('blue','gray'), cc_flux)
    ResidualBeforeSubtraction <- ResidualImage
    #-------- Subtraction in visibility
    cc_Visibility <- fft(cc_flux, inverse=T)* sqrt(Npix/2)* BLweight / sum(BLweight)
    ResidualVisibility <- Visibility - cc_Visibility
    #-------- Residual Image
    ResidualImage <- Re(fft(ResidualVisibility))/sqrt(Npix/2)
    PNGsnapShot(sprintf('ResidualImage/Residual%04d.png', 3*loop_index-1), list(seq(1,Npix)-(513-cc_index), 1:Npix, 1:Npix), list(Beam* ResidualImage[cc_index] * cleanGain* 4, dirtyImage, ResidualBeforeSubtraction), min(dirtyImage), max(dirtyImage), c('red', 'blue', 'gray'), cc_flux)
    #-------- plot residual
    PNGsnapShot(sprintf('ResidualImage/Residual%04d.png', 3*loop_index), list(seq(1,Npix)-(513-cc_index), 1:Npix, 1:Npix), list(Beam* ResidualImage[cc_index] * cleanGain* 4, dirtyImage, ResidualImage), min(dirtyImage), max(dirtyImage), c('red', 'blue','gray'), cc_flux)
}

#--------  MP4 movie
input_files <- list.files("ResidualImage", full.names = TRUE)
av_encode_video(input = input_files, framerate = 9, output = "Residuals.mp4")
