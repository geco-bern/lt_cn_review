### predict Nmass using CRU climate
# load("/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/pre CN global.RData")
# save.image("/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/pre CN global.RData")

# tr.cli <- read.csv('/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/CN cli.csv')
# tr.cli$mPAR0 <- tr.cli$mpar0*1000000/86400 # umol/m2/s
# cn.clm <- read.csv('/Volumes/Backup Plus/The Alchmist/ZCo-CN-Beni/try pft cn.csv')

load(paste0(here::here(), "/data/pre_CN_global.RData"))

library(ggplot2)

# function definitions----------------------
{
  # calculate air pressure in Pa
  calc_patm <- function( elv ){
    #-----------------------------------------------------------------------
    # Input:    - elevation, m (elv)
    # Output:   - float, atmospheric pressure at elevation 'elv', Pa (patm)
    # Features: Returns the atmospheric pressure as a function of elevation
    #           and standard atmosphere (1013.25 hPa)
    # Depends:  - connect_sql
    #           - flux_to_grid
    #           - get_data_point
    #           - get_msvidx
    # Ref:      Allen et al. (1998)
    #-----------------------------------------------------------------------
    
    # Define constants:
    kPo <- 101325   # standard atmosphere, Pa (Allen, 1973)
    kTo <- 298.15   # base temperature, K (Prentice, unpublished)
    kL <- 0.0065    # temperature lapse rate, K/m (Allen, 1973)
    kG <- 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
    kR <- 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
    kMa <- 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
    
    # Convert elevation to pressure, Pa:
    patm <- kPo*(1.0 - kL*elv/kTo)**(kG*kMa/(kR*kL))
    
    return (patm)
  }
  
  # calculate K (MM coefficient of Rubisco) in Pa
  calc_k <- function(temp, patm) {
    #-----------------------------------------------------------------------
    # Input:    - float, air temperature, deg C (temp)
    #           - float, atmospheric pressure, Pa (patm)
    # Output:   float, Pa (mmk)
    # Features: Returns the temperature & pressure dependent Michaelis-Menten
    #           coefficient, K (Pa).
    # Ref:      Bernacchi et al. (2001), Improved temperature response 
    #           functions for models of Rubisco-limited photosynthesis, 
    #           Plant, Cell and Environment, 24, 253--259.
    #-----------------------------------------------------------------------
    
    # Define constants
    kc25 <- 39.97      # Pa, assuming 25 deg C & 98.716 kPa
    ko25 <- 2.748e4    # Pa, assuming 25 deg C & 98.716 kPa
    dhac <- 79430      # J/mol
    dhao <- 36380      # J/mol
    kR   <- 8.3145     # J/mol/K
    kco  <- 2.09476e5  # ppm, US Standard Atmosphere
    
    vc <- kc25*exp(dhac*(temp - 25.0)/(298.15*kR*(temp + 273.15)))
    vo <- ko25*exp(dhao*(temp - 25.0)/(298.15*kR*(temp + 273.15)))
    k  <- vc*(1 + kco*(1e-6)*patm/vo)
    
    return(k)
    
  }
  
  # calculate Gstar (CO2 compensation point) in Pa
  calc_gstar_gepisat <- function( temp ) {
    #-----------------------------------------------------------------------
    # Input:    float, air temperature, degrees C (tc)
    # Output:   float, gamma-star, Pa (gs)
    # Features: Returns the temperature-dependent photorespiratory 
    #           compensation point, Gamma star (Pascals), based on constants 
    #           derived from Bernacchi et al. (2001) study.
    # Ref:      Bernacchi et al. (2001), Improved temperature response 
    #           functions for models of Rubisco-limited photosynthesis, 
    #           Plant, Cell and Environment, 24, 253--259.
    #-----------------------------------------------------------------------
    
    # Define constants
    gs25 <- 4.220    # Pa, assuming 25 deg C & 98.716 kPa)
    dha  <- 37830    # J/mol
    kR   <- 8.3145   # J/mol/K
    
    gs <- gs25 * exp( dha * ( temp - 25.0 ) / ( 298.15 * kR * ( temp + 273.15 ) ) )
    
    return( gs )
    
  }
  
  # conver CO2 from ppm to Pa
  co2_to_ca <- function( co2, patm ){
    #-----------------------------------------------------------------------
    # Input:    - float, annual atm. CO2, ppm (co2)
    #           - float, monthly atm. pressure, Pa (patm)
    # Output:   - ca in units of Pa
    # Features: Converts ca (ambient CO2) from ppm to Pa.
    #-----------------------------------------------------------------------
    ca   <- ( 1.e-6 ) * co2 * patm         # Pa, atms. CO2
    return( ca )
  }
}

# Arrhenius function modified terms
Arrhenius.t <- function(dS,Tl,Tref){
  Hd <- 200000 # J mol-1  
  R <- 8.3145 # J mol-1 K-1
  #        x<-A/C
  x <-(1+exp((Tl*dS-Hd)/(Tl*R)))/(1+exp((Tref*dS-Hd)/(Tref*R)))
  return(x)
}

# Parameters------------------
Ha.v <- 72000
R <- 8.3145 # J mol-1 K-1
dS.v <- 668.39 - 1.07*(tr.cli$mgdd0 + 273.15)

# Vcmax25 calculation------------
phi0 <- (0.352+0.021*tr.cli$mgdd0-3.4*10^(-4)*(tr.cli$mgdd0)^2)/8 # the tempereture-dependence of instrinsic quantum yield of C3
c <- 0.41
CO2.ppm <- 400 
K <- calc_k(tr.cli$mgdd0,calc_patm(tr.cli$elv))
Gstar <- calc_gstar_gepisat(tr.cli$mgdd0)
pa <- co2_to_ca(CO2.ppm,calc_patm(tr.cli$elv))
beta <- 146
f1 <- exp(-0.0227*(tr.cli$mgdd0-25)) # the viscosity of water relative to its value at 25˚C
g1 <- sqrt(beta*(Gstar+K)/(1.6*f1))
x <- Gstar/pa + (1-Gstar/pa)*g1/(g1+sqrt(tr.cli$mVPD0*1000))
ci <- x*pa
mc <- (ci - Gstar)/(ci + K)
m <- (ci - Gstar)/(ci + 2*Gstar)
Vc.gt <- phi0*tr.cli$mPAR0*(ci + K)/(ci + 2*Gstar)
Tref.Vcmax<- 25 + 273.15
Tl <- tr.cli$mgdd0 + 273.15 # mgdd0
Vcmax25 <- Vc.gt/exp(Ha.v*(Tl - Tref.Vcmax)/(Tref.Vcmax*R*Tl))*Arrhenius.t(dS.v,Tl,Tref.Vcmax) #umol/m2/s
tr.cli$pre.Vc25 <- Vcmax25

# LMA calculation--------------------
## deciduous--------
k1 = 30  # unit: g biomass / molC
u = 768
hT <- ( exp(Ha.v*(tr.cli$mgdd0 - 25)/(298.15*R*(tr.cli$mgdd0 + 273.15))) ) 
      * ((1 + exp((298.15*dS.v - 200000)/(298.15*R)))
        /(1 + exp(((tr.cli$mgdd0 + 273.15)*dS.v - 200000)/((tr.cli$mgdd0 + 273.15)*R))) 
         )  # Vc = Vc25*fv

xTde <- hT * mc / (phi0 * m)

Cde <- log(k1) + log(365) - log(u) - log(xTde)

# Eq. 45 in Wang et al. https://doi.org/10.1101/2021.02.07.430028
LMAde <- exp(log(tr.cli$f) 
             + log(tr.cli$mpar0) 
             - 0.052 * tr.cli$mgdd0 
             + Cde
             )

## evergreen---------
# Eq. ??? in Wang et al. https://doi.org/10.1101/2021.02.07.430028
LMAev <- exp( 0.5 * log(tr.cli$mpar0) 
              - 0.013 * tr.cli$mgdd0 
              - 0.62 * log(tr.cli$alpha) 
              + 0.25 * log(tr.cli$f) 
              + 3.759
              )

tr.cli$LMAde <- LMAde
tr.cli$LMAev <- LMAev
for (i in 1:nrow(tr.cli)) {
  if(tr.cli$LeafPhenology=='evergreen'){
    tr.cli[i,29] <- tr.cli[i,28]
  } else if(tr.cli$LeafPhenology=='deciduous'){
    tr.cli[i,29] <- tr.cli[i,27]
  } else if(tr.cli$LeafPhenology=='deciduous/evergreen'){
    tr.cli[i,29] <- (tr.cli[i,27] + tr.cli[i,28])/2
  } else{
    tr.cli[i,29] <- NA
  }
}
colnames(tr.cli)[29] <- 'pre.LMA'

# CN from Vcmax25 and LMA ------------------------------------------------------
tr.cli$pre.Nm <- 0.4796606 * (tr.cli$pre.Vc25 / tr.cli$pre.LMA) + 1.849738 
tr.cli$pre.cn <- 48 / tr.cli$pre.Nm  # global Cmass = 48 according to Dong Ning
hist(tr.cli$pre.cn)

cn.comp <- na.omit(tr.cli[,c(2:7,17,11,31)])
cn.obs <- cn.comp[,c(1:8)]
cn.pre <- cn.comp[,c(1:5,9,7,8)]
cn.obs$class <- 'observation'
cn.pre$class <- 'prediction'
colnames(cn.pre)[6] <- 'CN'

cn.comp2 <- rbind(cn.obs, cn.pre)
cn.comp3 <- merge(cn.comp2, cn.clm, by = 'PFTs', all=T )

# Plot C:N by biome-PFT---------------------------------------------------------
ggplot() + 
  geom_boxplot(
    data = cn.comp3, 
    aes(x = as.factor(PFTs), y = CN, fill = class),
    outlier.shape = NA, alpha = 0.6
    ) + 
  geom_point(
    data = cn.comp3 |> 
      select(PFTs, model) |> 
      unique(), 
    aes(x = as.factor(PFTs), y = model, color = 'salmon'),
    size = 3, 
    # shape = 4, 
  ) + 
  labs(x = "", y = 'Leaf C:N ratio') + 
  scale_y_continuous(limits = c(10, 60)) + 
  scale_fill_manual(values = c("#777055ff", "#29a274ff")) + 
  scale_color_manual(values = c('salmon'),
                     labels = c('LSM parameters')
                     ) +
  theme(legend.title = element_blank(), 
        legend.position = 'top', 
        legend.text = element_text(size = 10),
        legend.text.align = 0, 
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=13), 
        axis.text = element_text(size = 10)) + 
  coord_flip()

ggsave(paste0(here::here(), "/fig/cn_leaf_pft.pdf"), width = 6, height = 6)
ggsave(paste0(here::here(), "/fig/cn_leaf_pft.png"), width = 6, height = 6)
