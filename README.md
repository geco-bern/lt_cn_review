# lt_cn_review

Analyses of the review paper on C-N modelling by the LEMONTREE N cycle WG.

Author: Benjamin Stocker, Uni Bern

This consists of different analyses, implemented in separate scripts and notebooks:

## DGVM analysis of land C sink trends

Done in a separate GitHub repository: `stineb/dgvm_analyses`, using file `terr_sink_trends_gcb.Rmd`.

## Meta-analysis of ecosystem experiments

Using MESI data from eCO2 experiments (all variables): `vignettes/analysis_mesi_co2.Rmd` (previously `vignettes/analysis_mesi.Rmd`)

Using MESI + NutNet data from N-fert experiments (subset of variables): `vignettes/analysis_mesi_nutnet_nfert.Rmd` (previously `vignettes/Nfert_merge_MESI_NutNet_eperkowski.Rmd`)

Using Liang et al., 2020 data from N-fertilisation experiments (complementary subset of variables): `vignettes/analysis_liang_nfert.Rmd`

Combine plots for publication figure: `vignettes/publication_figures_metaanalysis.Rmd`

## Field data analysis

Using data from Dong et al., 2022: `vignettes/analysis_leafn_vcmax_field.Rmd`

## CN-model simulations and evaluation

CN-model simulations are done using a the CN-model, available in a separate GitHub repository (`stineb/rsofun`, branch `cnmodel`). Get this and build the rsofun R package to run rsofun for analyses in this repository.

Run scripts `analysis/exp_co2_cnmodel.R` and `analysis/exp_nfert_cnmodel.R` (both in this repository) to generate model outputs. Adjust the file to which outputs are written by modifying the commit code. The commit code should correspond to the latest commit of the rsofun repository, used for running the cnmodel simulations here. Instead of commit code, may also use tags from the rsofun repository.

Then, make a rough "ballpark comparison" of the outputs against data using `vignettes/cnmodel_benchmark_ballpark.Rmd`, and finally the model evaluation against MESI data using `vignettes/comparison_mesi_cnmodel.Rmd`.

## Additional analyses

(Not currently used)

### C:N by biome

Code from Huiying Xu (edited) available by `analysis/pre_CN_global.R` (loading file `data/pre_CN_global.RData` obtained by Huiying by email). See also [this Notion page](https://www.notion.so/computationales/Document-C-N-prediction-1796d17805784b109957bc82a03a1c62?pvs=4) (internal-only).

### Soil decomposition rate

Meta-analysis using data from Van Groenigen et al., 2014, see `vignettes/analysis_vangroenigen.Rmd`.
