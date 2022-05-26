# CommentEMD
# Matlab code
emd_Matlab_20210729.m is the Matlab code for executing causal decomposition for the time series with Moran effects
(two non-interacting time series driven by the same factors) 

Causal decomposition and eemd functions are from Yang et al. 2018 Nature Commun. and Wang et al. 2014 Physica A, respectively
# R code
Time series data with Moran effects (N1.txt, N2.txt, time_moran.txt) were produced by R script, Moran_model.R

CCM analysis for analyzing time series data with Moran effects (Fig. 1c) and white noises (Fig. 2) were included in the R code, comment_NatureC_20210729.R

We also re-analyzed all the examples presented in Yang et al. (2018) Nature Commun. by standard lag-CCM analysis. (Scripts were recorded in 'lag_CCM.R')

Datasets included the folder, Data_files, were all open-accessed and can be downloaded in the links offered by Yang et al. (2018) 

We conducted all the EDM analyses based on rEDM package ver 1.2.3. https://github.com/cran/rEDM/releases/tag/1.2.3
