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
* Population: https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/prof/details/download-telecharger/comp/GetFile.cfm?Lang=E&FILETYPE=CSV&GEONO=058 (2016 Census)
* Daylight: https://www.timeanddate.com/sun/

#### Covariate Descriptions
* intercept : all 1s
* popdensity : population / area from 2016 Census of Canada
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
* maxtemp.wk : weekly maximum temperature
* mintemp.wk : weekly minimum temperature
* meantemp.wk : weekly mean temperature
* precip.wk : weekly precipitation
* wks.mintemp.below.freezing :
* wks.meantemp.below.freezing :
* maxtemp.winter.avg :
* mintemp.winter.avg :
* level.avg.avg :
* level.avg.max :
* level.range.max :
* flow.max.avg :
* daylight :
* X.surveys :
* Ae.prop :
* Cx.prop :
* Cq.prop :
* catch.rate :
* trees :
* veg :
* agri :
* urban :
* water :
* bird.effort :
* bird.abundance :
* bird.richness :
* bird.shannon :
* bird.simpson :

