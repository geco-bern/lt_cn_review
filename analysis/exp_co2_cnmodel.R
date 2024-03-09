library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)
library(lubridate)

# detach("package:rsofun", unload = TRUE)
library(rsofun)

## Parameters ------------------------
# # Hand-tuned to give reasonable results for CH-Oe1 forcing
# pars <- list(
# 
#   # P-model
#   kphio                 = 0.08,       # setup ORG in Stocker et al. 2020 GMD
#   kphio_par_a           = -0.001,        # set to zero to disable temperature-dependence of kphio
#   kphio_par_b           = 20.0,
#   soilm_thetastar       = 0.6 * 240,  # to recover old setup with soil moisture stress
#   soilm_betao           = 0.0,
#   beta_unitcostratio    = 146.0,
#   rd_to_vcmax           = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous
#   tau_acclim            = 30.0,
#   kc_jmax               = 0.41,
# 
#   # Plant
#   f_nretain             = 0.500000,
#   fpc_tree_max          = 0.950000,
#   growtheff             = 0.600000,
#   r_root                = 2*0.913000,
#   r_sapw                = 2*0.044000,
#   exurate               = 0.003000,
# 
#   k_decay_leaf          = 0.5,
#   k_decay_root          = 0.5,
#   k_decay_labl          = 0.1,
#   k_decay_sapw          = 0.02,
# 
#   r_cton_root           = 37.0000,
#   r_cton_wood           = 350,
#   r_cton_seed           = 15.0000,
#   nv_vcmax25            = 0.02 * 13681.77, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.695; previously: 5000.0,
#   ncw_min               = 0.056, #0.08 * 1.116222, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.691; previously used: 0.056,
#   r_n_cw_v              = 1.23223, #0, # assumed that LMA is independent of Vcmax25; previously: 0.1,
#   r_ctostructn_leaf     = 40, #1.3 * 45.84125, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.699; previously used: 80.0000,
#   kbeer                 = 0.500000,
# 
#   # Phenology (should be PFT-specific)
#   gddbase               = 5.0,
#   ramp                  = 0.0,
#   phentype              = 2.0,
# 
#   # Soil physics (should be derived from params_soil, fsand, fclay, forg, fgravel)
#   perc_k1               = 5.0,
#   thdiff_wp             = 0.2,
#   thdiff_whc15          = 0.8,
#   thdiff_fc             = 0.4,
#   forg                  = 0.01,
#   wbwp                  = 0.029,
#   por                   = 0.421,
#   fsand                 = 0.82,
#   fclay                 = 0.06,
#   fsilt                 = 0.12,
# 
#   # Water and energy balance
#   kA                    = 107,
#   kalb_sw               = 0.17,
#   kalb_vis              = 0.03,
#   kb                    = 0.20,
#   kc                    = 0.25,
#   kCw                   = 1.05,
#   kd                    = 0.50,
#   ke                    = 0.0167,
#   keps                  = 23.44,
#   kWm                   = 220.0,
#   kw                    = 0.26,
#   komega                = 283.0,
#   maxmeltrate           = 3.0,
# 
#   # Soil BGC
#   klitt_af10            = 1.2,
#   klitt_as10            = 0.35,
#   klitt_bg10            = 0.35,
#   kexu10                = 50.0,
#   ksoil_fs10            = 0.021,
#   ksoil_sl10            = 7.0e-04,
#   ntoc_crit1            = 0.45,
#   ntoc_crit2            = 0.76,
#   cton_microb           = 10.0,
#   cton_soil             = 9.77,
#   fastfrac              = 0.985,
# 
#   # N uptake
#   # eff_nup               = 0.600000,  # original value
#   eff_nup               = 0.005000,
#   minimumcostfix        = 1.000000,
#   fixoptimum            = 25.15000,
#   a_param_fix           = -3.62000,
#   b_param_fix           = 0.270000,
# 
#   # Inorganic N transformations
#   maxnitr               = 0.05,
#   non                   = 0.01,
#   n2on                  = 0.0005,
#   kn                    = 83.0,
#   kdoc                  = 17.0,
#   docmax                = 1.0,
#   dnitr2n2o             = 0.01,
# 
#   # Additional parameters - previously forgotten
#   frac_leaf             = 0.5,           # after wood allocation
#   frac_wood             = 0.5,           # highest priority in allocation
#   frac_avl_labl         = 0.1,
# 
#   # for development
#   tmppar                = 9999,
# 
#   # simple N uptake module parameters
#   nuptake_kc            = 600,
#   nuptake_kv            = 200,
#   nuptake_vmax          = 10
# 
# )

# in LT review paper v3
pars <- list(
  
  # P-model
  kphio                 = 0.04998,    # setup ORG in Stocker et al. 2020 GMD
  kphio_par_a           = 0.0,        # set to zero to disable temperature-dependence of kphio
  kphio_par_b           = 1.0,
  soilm_thetastar       = 0.6 * 240,  # to recover old setup with soil moisture stress
  soilm_betao           = 0.0,
  beta_unitcostratio    = 146.0,
  rd_to_vcmax           = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous
  tau_acclim            = 30.0,
  kc_jmax               = 0.41,
  
  # Plant
  f_nretain             = 0.500000,
  fpc_tree_max          = 0.950000,
  growtheff             = 0.600000,
  r_root                = 2*0.913000,
  r_sapw                = 2*0.044000,
  exurate               = 0.003000,
  
  k_decay_leaf          = 1.90000,
  k_decay_root          = 1.90000,
  k_decay_labl          = 1.90000,
  k_decay_sapw          = 1.90000,
  
  r_cton_root           = 37.0000,
  r_cton_wood           = 100.000,
  r_cton_seed           = 15.0000,
  nv_vcmax25            = 0.02 * 13681.77, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.695; previously: 5000.0,
  ncw_min               = 0.08 * 1.116222, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.691; previously used: 0.056,
  r_n_cw_v              = 0, # assumed that LMA is independent of Vcmax25; previously: 0.1,
  r_ctostructn_leaf     = 1.3 * 45.84125, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.699; previously used: 80.0000,
  kbeer                 = 0.500000,
  
  # Phenology (should be PFT-specific)
  gddbase               = 5.0,
  ramp                  = 0.0,
  phentype              = 2.0,
  
  # Soil physics (should be derived from params_soil, fsand, fclay, forg, fgravel)
  perc_k1               = 5.0,
  thdiff_wp             = 0.2,
  thdiff_whc15          = 0.8,
  thdiff_fc             = 0.4,
  forg                  = 0.01,
  wbwp                  = 0.029,
  por                   = 0.421,
  fsand                 = 0.82,
  fclay                 = 0.06,
  fsilt                 = 0.12,
  
  # Water and energy balance
  kA                    = 107,
  kalb_sw               = 0.17,
  kalb_vis              = 0.03,
  kb                    = 0.20,
  kc                    = 0.25,
  kCw                   = 1.05,
  kd                    = 0.50,
  ke                    = 0.0167,
  keps                  = 23.44,
  kWm                   = 220.0,
  kw                    = 0.26,
  komega                = 283.0,
  maxmeltrate           = 3.0,
  
  # Soil BGC
  klitt_af10            = 1.2,
  klitt_as10            = 0.35,
  klitt_bg10            = 0.35,
  kexu10                = 50.0,
  ksoil_fs10            = 0.021,
  ksoil_sl10            = 7.0e-04,
  ntoc_crit1            = 0.45,
  ntoc_crit2            = 0.76,
  cton_microb           = 10.0,
  cton_soil             = 9.77,
  fastfrac              = 0.985,
  
  # N uptake
  eff_nup               = 0.0001000,
  minimumcostfix        = 1.000000,
  fixoptimum            = 25.15000,
  a_param_fix           = -3.62000,
  b_param_fix           = 0.270000,
  
  # Inorganic N transformations (re-interpreted for simple ntransform model)
  maxnitr               =  0.000005,
  
  # Inorganic N transformations for full ntransform model (not used in simple model)
  non                   = 0.01,
  n2on                  = 0.0005,
  kn                    = 83.0,
  kdoc                  = 17.0,
  docmax                = 1.0,
  dnitr2n2o             = 0.01,
  
  # Additional parameters - previously forgotten
  frac_leaf             = 0.5,           # after wood allocation
  frac_wood             = 0,           # highest priority in allocation
  frac_avl_labl         = 0.1,
  
  # for development
  tmppar                = 9999,
  
  # simple N uptake module parameters
  nuptake_kc            = 250,
  nuptake_kv            = 5,
  nuptake_vmax          = 0.15
  
)


## Forcing ------------------------
## add new required columns to forcing 
tmp <- rsofun::p_model_drivers |> 
  mutate(forcing = purrr::map(forcing, ~mutate(., 
                                               fharv = 0.0,
                                               dno3 = 0.0016, # 0.6 / 365,
                                               dnh4 = 0.0014, # 0.5 / 365
  )))

## interactive C-N cycling
tmp$params_siml[[1]]$c_only <- FALSE

## no soil moisture stress
tmp$params_siml[[1]]$soilmstress <- FALSE

### Harvesting and seed input ----------
use_cseed <- 0 # 100
cn_seed <- 20
use_nseed <- use_cseed / cn_seed

tmp$forcing[[1]] <- tmp$forcing[[1]] |>
  mutate(fharv = ifelse(month(date) == 7 & mday(date) == 15, 0.0, 0.0),
         cseed = ifelse(month(date) == 3 & mday(date) == 15, use_cseed, 0.0),
         nseed = ifelse(month(date) == 3 & mday(date) == 15, use_nseed, 0.0)) 

## check visually
tmp$forcing[[1]] |>
  ggplot(aes(date, fharv)) +
  geom_line()

## no spinup, 1 year transient run
tmp$params_siml[[1]]$spinupyears <- 2000
tmp$params_siml[[1]]$recycle <- 5

## Synthetic forcing: Constant climate in all days -----------------------
df_growingseason_mean <- tmp$forcing[[1]] |>
  filter(temp > 5) |>
  summarise(across(where(is.double), .fns = mean))
df_mean <- tmp$forcing[[1]] |>
  summarise(across(where(is.double), .fns = mean))

tmp$forcing[[1]] <- tmp$forcing[[1]] |>
  mutate(temp = df_growingseason_mean$temp,
         prec = df_mean$prec,
         vpd = df_growingseason_mean$vpd,
         ppfd = df_mean$ppfd,
         patm = df_growingseason_mean$patm,
         ccov_int = df_growingseason_mean$ccov_int,
         ccov = df_growingseason_mean$ccov,
         snow = df_mean$snow,
         rain = df_mean$rain,
         fapar = df_mean$fapar,
         co2 = df_growingseason_mean$co2,
         tmin = df_growingseason_mean$tmin,
         tmax = df_growingseason_mean$tmax,
  )

## Repeat last year's forcing N times -----------------------
n_ext <- 100
df_tmp <- tmp$forcing[[1]]
for (idx in seq(n_ext)){
  df_tmp <- bind_rows(
    df_tmp,
    df_tmp |>
      tail(365) |>
      mutate(date = date + years(1))
  )
}
tmp$forcing[[1]] <- df_tmp


## increase CO2 from 2010 -----------------------
elevate_co2 <- function(day){
  yy <- 2 - 1 / (1 + exp(0.03*(day-14610)))
  return(yy)
}

ggplot() +
  geom_function(fun = elevate_co2) +
  xlim(12000, 16000) +
  geom_vline(xintercept = 0, linetype = "dotted")

tmp$forcing[[1]] <- tmp$forcing[[1]] |>
  mutate(date2 = as.numeric(date)) |>
  mutate(co2 = co2 * elevate_co2(date2)) |>
  select(-date2)

tmp$forcing[[1]] |>
  head(5000) |>
  ggplot(aes(date, co2)) +
  geom_line()


## Model run ------------------------
output <- runread_cnmodel_f(
  tmp,
  par = pars
) 

output <- output$data[[1]]

## Visualisations  ------------------------
### Time series ---------------------------
gg1 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, lai)) + 
  geom_line()
gg2 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, gpp)) + 
  geom_line()
gg3 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, cleaf/nleaf)) + 
  geom_line()
gg4 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, croot/cleaf)) + 
  geom_line()

gg1 / gg2 / gg3 / gg4

gg5 <- output |>  
  as_tibble() |> 
  ggplot(aes(date, x1)) + 
  geom_line()
gg6 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, x2)) + 
  geom_line()
gg7 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, x3)) + 
  geom_line()
gg8 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, x4)) + 
  geom_line()

gg5 / gg6 / gg7 / gg8

gg9 <- output |>  
  as_tibble() |> 
  ggplot(aes(date, pnh4 + pno3)) + 
  geom_line()
gg10 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, en2o)) + 
  geom_line()
gg11 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, nup)) + 
  geom_line()
gg12 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, netmin)) + 
  geom_line()

gg9 / gg10 / gg11 / gg12

gg13 <- output |>  
  as_tibble() |> 
  ggplot(aes(date, cleaf/nleaf)) + 
  geom_line()
gg14 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, croot/nroot)) + 
  geom_line()
gg15 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, npp/nup)) + 
  geom_line()
gg16 <- output |> 
  as_tibble() |> 
  ggplot(aes(date, nloss/nup)) + 
  geom_line()

gg13 / gg14 / gg15 / gg16


### Response ratios ---------------------------
df_out <- output |> 
  mutate(leaf_cn = cleaf/nleaf, 
         root_shoot = croot/cleaf, 
         n_inorg = pno3 + pnh4,
         anpp = npp_leaf + npp_wood, 
         bnpp = npp_root + cex) |> 
  select(date, asat, gpp, vcmax, jmax, gs = gs_accl, narea, leaf_cn, lai, cleaf, 
         croot, root_shoot, nup, n_inorg, anpp, bnpp)

df_amb <- df_out |> 
  filter(year(date) < 2010) |> 
  summarise(across(where(is.numeric), mean))

df_ele <- df_out |> 
  filter(year(date) %in% 2010:2015) |> 
  summarise(across(where(is.numeric), mean))

df_ele2 <- df_out |> 
  filter(year(date) %in% 2100:2110) |> 
  summarise(across(where(is.numeric), mean))

df_exp <- bind_rows(df_amb, df_ele)
df_rr  <- log(df_exp[2,]/df_exp[1,]) |> 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "response") |> 
  mutate(variable = factor(variable, 
                           levels = rev(c("asat", "gpp", "vcmax", "jmax", "gs", "narea", 
                                          "leaf_cn", "lai", "cleaf", "anpp",
                                          "croot", "bnpp", "root_shoot", "nup", 
                                          "n_inorg"))))

df_exp2 <- bind_rows(df_amb, df_ele2)
df_rr2  <- log(df_exp2[2,]/df_exp2[1,]) |> 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "response") |> 
  mutate(variable = factor(variable, 
                           levels = rev(c("asat", "gpp", "vcmax", "jmax", "gs", "narea", 
                                          "leaf_cn", "lai", "cleaf", "anpp",
                                          "croot", "bnpp", "root_shoot", "nup", 
                                          "n_inorg"))))

ggrr <- ggplot() +
  geom_point(aes(variable, response), data = df_rr2, size = 2, color = "grey50") +
  geom_point(aes(variable, response), data = df_rr, size = 2) +
  geom_hline( yintercept = 0.0, size = 0.5, linetype = "dotted" ) +
  labs(x = "Variable", y = "Log Response Ratio") +
  coord_flip() +
  labs(title = "cnmodel prediction", subtitle = "Response to eCO2")

ggsave(paste0(here::here(), "/fig/response_co2_cnmodel_11cc17c508eb30dd8db91e5eafff1a6b4880e9a8.pdf"))

### Spinup-----------------
## read (experimental) files
aout <- read_fwf(file = "out/out_rsofun.a.csoil.txt", col_types = "in") |>
  setNames(c("year", "csoil")) |>
  left_join(
    read_fwf(file = "out/out_rsofun.a.nsoil.txt", col_types = "in") |>
      setNames(c("year", "nsoil")),
    by = "year"
  )

aout |>
  # slice(1000:2008) |> 
  ggplot(aes(year, csoil)) +
  
  # first soil equilibration year
  geom_vline(xintercept = 600, linetype = "dotted") +
  
  # start free allocation
  geom_vline(xintercept = 900, linetype = "dotted") +
  
  # second soil equilibration year
  geom_vline(xintercept = 1500, linetype = "dotted") +
  
  geom_line() + 
  theme_classic()

## Write output to file with commit ID --------------------
# readr::write_csv(as_tibble(output), file = paste0(here::here(), "/data/output_cnmodel_co2_50c01ecbac0ad20114dc9cc28d67006af45f128e.csv"))
readr::write_csv(as_tibble(output), file = paste0(here::here(), "/data/output_cnmodel_co2_11cc17c508eb30dd8db91e5eafff1a6b4880e9a8.csv"))
