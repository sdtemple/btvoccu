# btvoccu
Bayesian Time-Varying Occupancy Modeling

### Package Summary
This package facilitates Bayesian analysis for the time-varying occupancy with occupancy and detection components modeled as GLMs with logit or probit links. The main function is `btvoccu`. It outputs a list of length 15 which includes model details and posterior samples. All other functions fall into three categories: (1) manipulating a 4-dimensional array, (2) summarizing a Bayesian model, or (3) plotting results. Analysis is meant for 4-dimensional arrays, in which the (ordered) major dimensions are sites, seasons, periods, and covariates/response.

### Modeling Pipeline
We provide a vignette as a practical tutorial with simulated data. In general, modeling with the `btvoccu` package proceeds as follows:
1. Gather relevant covariates (your own code)
2. Format covariates and the presence/absence response into 4-dimensional arrays
3. Transform covariates (`append_*()` functions)
4. Split data into training, validation, and testing datasets. (`subset_4darray()`, `split_4darray()`)
5. Fit models with training data (`btvoccu()`)
6. Check MCMC diagnostics (`plot_trace()`, `posterior_effects()`, `posterior_correlations()`, `posterior_variances()`)
7. Evaluate models on validation dataset (`waic_score()`, `posterior_check()`)
8. Select final model(s)
9. Fit final model(s) with training and validation data (`btvoccu()`)
10. Check MCMC diagnostics (`plot_trace()`, `posterior_effects()`, `posterior_correlations()`, `posterior_variances()`)
11. Evaluate models on testing dataset (`waic_score()`, `posterior_check()`)
12. Plot (`plot_btvoccu()`, `plot_covariate()`, your own mapping code)

### WNV in Ontario

#### Data Sources
* Mosquitoes: trap data from Public Health Ontario; provided by Deborah Shutt, PhD
* Hydrometry: https://github.com/paleolimbot/hydatr
* Weather: https://github.com/paleolimbot/rclimateca
* Land: http://www.earthenv.org/landcover (EarthEnv)
* Climate: https://datadryad.org/stash/dataset/doi:10.5061/dryad.s2v81 (MERRAclim 2000-2009)
* Birds: https://ebird.org/science/use-ebird-data/download-ebird-data-products (eBird)
* Population: https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/prof/details/download-telecharger/comp/GetFile.cfm?Lang=E&FILETYPE=CSV&GEONO=058 (2016 census)
* Daylight: https://www.timeanddate.com/sun/

#### Covariate Descriptions
All covariates are aggregated to be weekly statistics for the public health units of Ontario.
* intercept : all 1s
* popdensity : population / area from 2016 census
* mean.max.temp.warmest.m : max temperature of the warmest month, mean over the decade
* mean.min.temp.coldest.m : min temperature of the coldest month, mean over the decade
* mean.temp.range.a : temperature annual range, mean over the decade
* mean.mean.temp.wettest.q : mean temperature of the wettest quarter, mean over the decade
* mean.mean.temp.warmest.q : mean temperature of the warmest quarter, mean over the decade
* max.max.temp.warmest.q : max temperature of the warmest quarter, max over the decade
* max.min.temp.coldest.m : min temperature of the coldest month, max over the decade
* max.temp.range.a : temperature annual range, max over the decade
* max.mean.temp.warmest.q : mean temperature of the warmest quarter, max over the decade
* max.precip.wettest.q : precipitation of the wettest quarter, max over the decade
* min.max.temp.warmest.m : max temperature of the warmest month, min over the decade
* min.precip.wettest.m : precipitation of the wettest month, min of the decade
* maxtemp.wk : maximum temperature in Celsius (w/ imputed distance-weighted average when missing)
* mintemp.wk : minimum temperature in Celsius (w/ imputed distance-weighted average when missing)
* meantemp.wk : mean temperature in Celsius (w/ imputed distance-weighted average when missing)
* precip.wk : precipitation in millimeters (w/ imputed distance-weighted average when missing)
* wks.mintemp.below.freezing : # of weeks in which min temperature is less than 0 degrees Celsius between numbered weeks 1 to 17
* wks.meantemp.below.freezing : # of weeks in which mean temperature is less than 0 degrees Celsius between numbered weeks 1 to 17
* maxtemp.winter.avg : average max temperature over numbered weeks 1 to 16
* mintemp.winter.avg : average min temperature over numbered weeks 1 to 16
* level.avg.avg : water level, average over the week, average over the stations
* level.avg.max : water level, average over the week, max over the stations
* level.range.max : range of water level values, max over the stations
* flow.max.avg : flow, max over the week, average over the stations
* daylight : daylight time with linear interpolation based on 7 locations spread ~ 1.5 latitude degrees apart 
* X.surveys : # of traps visited (some traps may be visited more than once in a given week)
* Ae.prop : percent of trapped mosquitoes of the Aedes genus
* Cx.prop : percent of trapped mosquitoes of the Culex genus
* Cq.prop : percent of trapped mosquitoes of the Coquillettidia genus
* catch.rate : count of trapped mosquitoes over # of traps visited
* trees : percenof land cover inferred to be needleleaf, evergreen, deciduous, or mixed trees (1-4)
* veg : percent of land cover inferred to be herbacious or regularly flooded vegetation (6,8)
* agri : percent of land cover inferred to be cultivated and managed vegetation (7)
* urban : percent of land cover inferred to be urban/built-up (9)
* water : percent of land cover inferred to be open water (12)
* bird.effort : count of checklists (download eBird mobile app to better understand, and to have fun!)
* bird.abundance : count of birds observed
* bird.richness : count of bird species observed
* bird.shannon : https://en.wikipedia.org/wiki/Diversity_index#Shannon_index
* bird.simpson : https://en.wikipedia.org/wiki/Diversity_index#Simpson_index

#### Imputing and Average
* Hydrometry: This is daily data aggregated to the weekly scale. We average over the current week, the two weeks prior, and the two weeks after to impute a very small amount of missing covariates, predominantly in TSK public health unit. 
* Weather: If a weather station is missing a covariate (temperature, precipitation), we find its 5 nearest neighbors and impute a distance-based weighted average. Some nearest neighbors may be missing as well, so this average may be over 1 or up to 5 values. In aggregating to the public health unit, we only use imputed covariates if these are the only ones available for the unit. For public health units without weather stations, we impute with a distance-based weighted average over 1 or up to 3 values. Units without weather stations are generally in dense urban areas (Greater Toronto), so we expect this spatial imputation is a reasonable approximation. 
* Birds: Birding may fluctuate substantially on a weekly basis, so we average over the current and previous three weeks. This step imputes some missing values. For the remaining missing values, we impute using a median conditional on the public health unit and the epidemic week. For Shannon indices of 0 and Simpson indices of 1, we impute using the median conditional on the public health unit and the epidemic week. These values correspond to low birding effort (only one species observed). 
