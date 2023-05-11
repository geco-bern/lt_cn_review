### calc climate data for Nmass prediction globally
## trait data provided by Dong Ning
## climate data from CRU, didnt use Beni's Pmodel (https://computationales.github.io/ingestr/articles/run_pmodel_points.html)
load("/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/calc bioclim.RData")
save.image("/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/calc bioclim.RData")

library(ncdf4)
library(raster) # brick, extract
library(ggplot2)
library(dplyr)

# trait data point
tr.ori <- read.csv('/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/data_cn_pdf.csv')
site.info <- unique(tr.ori[,c(2,3,13)])
site <- site.info[,1:2]  # it needs 'lon' at first col, 'lat' second

# import functions and climate data
source("/Volumes/Backup Plus/Code/STASH/Scripts/daily.r")
source("/Volumes/Backup Plus/Code/STASH/Scripts/gddpar.r")
source("/Volumes/Backup Plus/Code/STASH/Scripts/gdd2.r")
source("/Volumes/Backup Plus/Code/STASH/Scripts/gdd alpha.r")
source("/Volumes/Backup Plus/Code/STASH/Scripts/PET.r")
source("/Volumes/Backup Plus/Code/STASH/Scripts/AET.r")  #daily prec as input, dealt in AET function
source("/Volumes/Backup Plus/Code/STASH/Scripts/monthly.r")  # mean daily of each month
source("/Volumes/Backup Plus/Code/STASH/Scripts/monthly sum.r")

setwd('/Volumes/Backup Plus/Data/CRU 4.05')
prec <- extract(brick("cru_ts4.05.2001.2010.pre.dat.nc", var='pre'), site)
tmn <- extract(brick("cru_ts4.05.2001.2010.tmn.dat.nc", var='tmn'), site)
tmx <- extract(brick("cru_ts4.05.2001.2010.tmx.dat.nc", var='tmx'), site)
cld <- extract(brick("cru_ts4.05.2001.2010.cld.dat.nc", var='cld'), site)
vap <- extract(brick("cru_ts4.05.2001.2010.vap.dat.nc", var='vap'), site)
# elevation
global.dem <- raster('/Volumes/Backup Plus/Data/DEM/elevation_noaa_bed_half_deg.tif')
elv <- extract(global.dem, site)
elv[elv<0] <- 0
site.info <- data.frame(site.info, elv)

mpre <- data.frame(site)
mtmn <- data.frame(site)
mtmx <- data.frame(site)
mcld <- data.frame(site)
mvap <- data.frame(site)
grid <- nrow(site)

# extract 2011-2020
n_grid <- length(site$lon)
for (i in 1:12) {
  for (j in 1:n_grid) {
    k <- seq(i,i+12*9,12)
    mpre[j,i+2] <- (prec[j,k[1]]+prec[j,k[2]]+prec[j,k[3]]+prec[j,k[4]]+prec[j,k[5]]+prec[j,k[6]]+prec[j,k[7]]+prec[j,k[8]]+prec[j,k[9]]+prec[j,k[10]])/10
    mtmn[j,i+2] <- (tmn[j,k[1]]+tmn[j,k[2]]+tmn[j,k[3]]+tmn[j,k[4]]+tmn[j,k[5]]+tmn[j,k[6]]+tmn[j,k[7]]+tmn[j,k[8]]+tmn[j,k[9]]+tmn[j,k[10]])/10
    mtmx[j,i+2] <- (tmx[j,k[1]]+tmx[j,k[2]]+tmx[j,k[3]]+tmx[j,k[4]]+tmx[j,k[5]]+tmx[j,k[6]]+tmx[j,k[7]]+tmx[j,k[8]]+tmx[j,k[9]]+tmx[j,k[10]])/10
    mcld[j,i+2] <- (cld[j,k[1]]+cld[j,k[2]]+cld[j,k[3]]+cld[j,k[4]]+cld[j,k[5]]+cld[j,k[6]]+cld[j,k[7]]+cld[j,k[8]]+cld[j,k[9]]+cld[j,k[10]])/10
    mvap[j,i+2] <- (vap[j,k[1]]+vap[j,k[2]]+vap[j,k[3]]+vap[j,k[4]]+vap[j,k[5]]+vap[j,k[6]]+vap[j,k[7]]+vap[j,k[8]]+vap[j,k[9]]+vap[j,k[10]])/10
  }
}

mtmean <- (mtmn[,3:14] + mtmx[,3:14])/2

# calculate MAT (mean annual temperature), MAP (mean, actually is total, annual precip)
MAT <- apply(mtmean,1,mean)
MAP <- apply(mpre[,3:14],1,sum)
cli.gs <- data.frame(grid_info[,2:5],MAP,MAT)

# calc sunshine hour %
mssh <- 100-mcld[,3:14]
mssh <- data.frame(site,mssh)

# calc daytime temperature using Tmin and Tmax
Tg_cal <- function(tmn,tmx,lat){
  #-----------------------------------------------------------------------
  # Input:   ?? (s;solar delclination), 12 month constant values
  # Input:   tmn, in monthly or daily degrees
  # Input:   tmx, in monthly or daily degrees
  # Input:   lat, in degrees
  # Output:  Tg,  in monthly or daily degrees
  # Features: Converts growth temperature from tmn to tmx
  #-----------------------------------------------------------------------
  s1 <- -20.205
  s2 <- -12.65
  s3 <- -1.95
  s4 <- 9.45
  s5 <- 18.35
  s6 <- 22.55
  s7 <- 20.75
  s8 <- 13.45
  s9 <- 2.9
  s10 <- -8.45
  s11 <- -17.85
  s12 <- -22.355
  s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
  x <-data.frame(matrix(nrow = nrow(tmn),ncol=ncol(tmn))) ##create x dataframe: nrow=No. of sites, ncol=timestep
  Tg <-data.frame(matrix(nrow = nrow(tmn),ncol=ncol(tmn)))##see above
  
  for (i in 1:ncol(tmn)){
    
    x[,i]<- -tan(pi*lat/180)*tan(s[i]*pi/180)
    Tg[,i]<-tmx[,i]*(0.5+(1-x[,i]^2)^(0.5)/(2*acos(x[,i])))+ tmn[,i]*(0.5-(1-x[,i]^2)^(0.5)/(2*acos(x[,i])))
    
  }
  
  Tg[Tg =="NaN"] <- NA
  return(Tg)
}
Tg <- Tg_cal(mtmn[,3:14], mtmx[,3:14], site$lat) # daytime temperature

# calc VPD
mes <- 0.611*exp(17.27*Tg/(Tg+237.3))  # in kPa
mea <- mvap[,3:14]/10  # in kPa
mvpd <- mes-mea  # in kPa, monthly VPDmax

# linear interoplating monthly climate input data to daily arrays
dvpd <- daily(mvpd, n_grid)
dprec <- daily(mpre[,3:14], n_grid)
dTmax <- daily(mtmx[,3:14], n_grid)
dTmin <- daily(mtmn[,3:14], n_grid)
dssh <- daily(mssh[,3:14], n_grid)
dTmean <- daily(mtmean, n_grid)
dTg <- daily(Tg, n_grid)

# calculate growing degree days (gdd) and 
# above certain temperature (e.g. 0 and 5) based on a daily step
outgddvpd0 <- gddpar(dTmean,dvpd,0) # above 0 degree
outTg0 <- gdd(dTg,0)
gdday0 <- matrix(data=unlist(outgddvpd0[2]),nrow=nrow(site)) # length of growing degree days above 0
mgdd0 <- matrix(data=unlist(outTg0[3]),nrow=nrow(site)) # mean growing degree days above 0
mVPD0 <- matrix(data=unlist(outgddvpd0[5]),nrow=nrow(site)) # mean VPD during the period above 

biocli <- data.frame(site.info, gdday0,mgdd0,mVPD0)

### calc alpha
# setting parameters:
epsilon <- 23.4*pi/180      # tilt
e <- 0.01675                # eccentricity
phi_ref <- 284*pi/180       # longitude of perihelion
fcap<-150                   # setting field capacity in mm

# calculate daily equilibirum evapotranspiration (dpet) 
# and daily photosynthetically active radiation (dPar)
outpet <- PET(site$lat,site.info$elv,dTg,dssh,epsilon,e,phi_ref)                 
deet <- matrix(data=unlist(outpet[1]),nrow=nrow(site))      # daily equilibrium evapotranspiration
dPar <- matrix(data=unlist(outpet[2]),nrow=nrow(site))      # daily photosynthetically active radiation, mol photon/m2/day
dcon <- matrix(data=unlist(outpet[3]),nrow=nrow(site))			# daily condensation (dew)

outgddpar0 <- gddpar(dTmean,dPar,0) # above 0 degree
mpar0 <- matrix(data=unlist(outgddpar0[5]),nrow=nrow(site)) # mean PAR during the period above 
biocli$mpar0 <- mpar0

# variables used for AET estimation
v <- matrix(data=unlist(outpet[4]),nrow=nrow(site))
u <- matrix(data=unlist(outpet[5]),nrow=nrow(site))
h0 <- matrix(data=unlist(outpet[6]),nrow=nrow(site))

# calculate annual equilibrium evapotranspiration and moisture index
yreet <- apply(deet,1,sum) 
yrpet <- yreet*1.26
hydra.na <- which(!is.na(yrpet)) # remove NA site
site.info2 <- site.info[hydra.na,]

# delete NA data
dprec <- dprec[hydra.na,]
deet <- deet[hydra.na,]
dcon <- dcon[hydra.na,]
v <- v[hydra.na,]
u <- u[hydra.na,]
h0 <- h0[hydra.na,]
MAP <- MAP[hydra.na]
yrpet <- yrpet[hydra.na]
yreet <- yreet[hydra.na]

grid2 <- nrow(site.info2)
daet <- AET(grid2,dprec,deet,dcon,v,u,h0,MAP,yrpet)
yraet <- apply(daet,1,sum)
alpha <- yraet/yrpet # alpha in the range of 0 and 1.0, simplified definition

site.info2$alpha <- alpha

biocli2 <- merge(biocli, site.info2[,3:5], by='site_id', all=T)
biocli2$f <- biocli2$gdday0/365
tr.cli <- merge(tr.ori, biocli2[,c(1,4:10)], by='site_id',all=T)

#write.csv(biocli2, '/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/biocli.csv')
write.csv(tr.cli, '/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/CN cli.csv')
