# btvoccu
bayesian time-varying occupancy modeling

### Package Summary
This package facilitates Bayesian analysis for the time-varying occupancy with occupancy and detection components modeled as GLMs with logit or probit links. The main function is `btvoccu`. It outputs a list of length 15 which includes model details and posterior samples. All other functions fall into three categories: (1) manipulating a 4-dimensional array, (2) summarizing a Bayesian model, or (3) plotting results. Analysis is meant for 4-dimensional arrays, in which the (ordered) major dimensions are sites, seasons, periods, and covariates/response.  

### WNV in Ontario

#### Data Sources
* Mosquitoes: trap data from Public Health Ontario; provided by Deborah Shutt, PhD
* Hydrology: https://github.com/paleolimbot/hydatr
* Weather: https://github.com/paleolimbot/rclimateca
* Land: http://www.earthenv.org/landcover (EarthEnv)
* Climate: https://datadryad.org/stash/dataset/doi:10.5061/dryad.s2v81 (MERRAclim)
* Birds: https://ebird.org/science/use-ebird-data/download-ebird-data-products (eBird)
* Population: https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/prof/details/download-telecharger/comp/GetFile.cfm?Lang=E&FILETYPE=CSV&GEONO=058 (2016 census)
* Daylight: https://www.timeanddate.com/sun/

#### Covariate Descriptions
All covariates are aggregated to be weekly statistics for the public health units of Ontario.
* intercept : all 1s
* popdensity : population / area from 2016 census
* mean.max.temp.warmest.m :
* mean.min.temp.coldest.m :
* mean.temp.range.a :
* mean.mean.temp.wettest.q :
* mean.mean.temp.warmest.q :
* max.max.temp.warmest.q :
* max.min.temp.coldest.m :
* max.temp.range.a :
* max.mean.temp.warmest.q :
* max.precip.wettest.q :
* min.max.temp.warmest.m :
* min.precip.wettest.m :
* maxtemp.wk : maximum temperature (w/ imputed distance-weighted average when missing)
* mintemp.wk : minimum temperature (w/ imputed distance-weighted average when missing)
* meantemp.wk : mean temperature (w/ imputed distance-weighted average when missing)
* precip.wk : precipitation (w/ imputed distance-weighted average when missing)
* wks.mintemp.below.freezing :
* wks.meantemp.below.freezing :
* maxtemp.winter.avg :
* mintemp.winter.avg :
* level.avg.avg :
* level.avg.max :
* level.range.max :
* flow.max.avg :
* daylight : daylight time with linear interpolation based on 7 locations spread ~ 1.5 latitude degrees apart 
* X.surveys :
* Ae.prop :
* Cx.prop :
* Cq.prop :
* catch.rate :
* trees : percenof land cover inferred to be needleleaf, evergreen, deciduous, or mixed trees (1-4)
* veg : percent of land cover inferred to be herbacious or regularly flooded vegetation (6,8)
* agri : percent of land cover inferred to be cultivated and managed vegetation (7)
* urban : percent of land cover inferred to be urban/built-up (9)
* water : percent of land cover inferred to be open water (12)
* bird.effort :
* bird.abundance :
* bird.richness :
* bird.shannon :
* bird.simpson :

