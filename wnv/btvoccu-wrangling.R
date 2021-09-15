# btvoccu Wrangling
# Seth Temple, sdtemple@lanl.gov
# August 26, 2021

library(btvoccu)

# Load --------------------------------------------------------------------

y <- readRDS("data/response.rds")
x <- readRDS("data/covariates.rds")

y <- apply(y, 1:4, as.numeric)
x <- apply(x, 1:4, as.numeric)

# Transform ---------------------------------------------------------------

# land cover
x <- append_transformed_term(x, "trees", function(x){x / 100}, "perc.")
x <- append_transformed_term(x, "veg", function(x){x / 100}, "perc.")
x <- append_transformed_term(x, "agri", function(x){x / 100}, "perc.")
x <- append_transformed_term(x, "urban", function(x){x / 100}, "perc.")
x <- append_transformed_term(x, "water", function(x){x / 100}, "perc.")
x <- append_transformed_term(x, "perc.trees", sqrt, "sqrt.")
x <- append_transformed_term(x, "perc.veg", sqrt, "sqrt.")
x <- append_transformed_term(x, "perc.agri", sqrt, "sqrt.")
x <- append_transformed_term(x, "perc.urban", sqrt, "sqrt.")
x <- append_transformed_term(x, "perc.water", sqrt, "sqrt.")

# population density
x <- append_transformed_term(x, "popdensity", log10, "log10.")
x <- append_transformed_term(x, "popdensity", sqrt, "sqrt.")

# traps
x <- append_transformed_term(x, "X.surveys", log, "log.")
x <- append_transformed_term(x, "X.surveys", sqrt,  "sqrt.")
x <- append_transformed_term(x, "catch.rate", sqrt, "sqrt.")
x <- append_transformed_term(x, "catch.rate", function(x){return(x**(1/3))}, "cbrt.")
x <- append_transformed_term(x, "Cx.prop", sqrt, "sqrt.")
x <- append_transformed_term(x, "Cx.prop", function(x){return(x**(1/3))}, "cbrt.")

# weather and climate
x <- append_transformed_term(x, "wks.meantemp.below.freezing", sqrt, "sqrt.")
x <- append_principal_components(x, dimnames(x)[[4]][5:16], 2, "climate")

# hydrology
x <- append_transformed_term(x, "level.avg.avg", sqrt, "sqrt.")
x <- append_transformed_term(x, "level.avg.max", sqrt, "sqrt.")

# eBird
x <- append_transformed_term(x, "bird.richness", function(x){log(x)}, "log.")
x <- append_transformed_term(x, "bird.abundance", function(x){log10(x)}, "log10.")
x <- append_transformed_term(x, "bird.simpson", function(x){1/x}, "inv.")
x <- append_transformed_term(x, "bird.simpson", function(x){1-x}, "gini.")

# lags
x <- append_lagged_terms(x, "meantemp.wk", c(2,4,6,8))
x <- append_lagged_terms(x, "precip.wk", c(2,4,6,8))
x <- append_lagged_terms(x, "sqrt.level.avg.avg", c(2,4,6,8))
x <- append_lagged_terms(x, "gini.bird.simpson", c(2,4,6,8))

# Save --------------------------------------------------------------------

# prune covariates
w <- subset_4darray(x, 4, c("intercept",
                            "popdensity",
                            "log10.popdensity",
                            "meantemp.wk",
                            "lag2meantemp.wk",
                            "lag4meantemp.wk",
                            "lag6meantemp.wk",
                            "lag8meantemp.wk",
                            "precip.wk",
                            "lag2precip.wk",
                            "lag4precip.wk",
                            "lag6precip.wk",
                            "lag8precip.wk",
                            "sqrt.level.avg.avg",
                            "lag2sqrt.level.avg.avg",
                            "lag4sqrt.level.avg.avg",
                            "lag6sqrt.level.avg.avg",
                            "lag8sqrt.level.avg.avg",
                            "wks.meantemp.below.freezing",
                            "sqrt.wks.meantemp.below.freezing",
                            "climate.pc1",
                            "climate.pc2",
                            "daylight",
                            "X.surveys",
                            "sqrt.X.surveys",
                            "log.X.surveys",
                            "Cx.prop",
                            "sqrt.Cx.prop",
                            "cbrt.Cx.prop",
                            "catch.rate",
                            "sqrt.catch.rate",
                            "cbrt.catch.rate",
                            "perc.trees",
                            "perc.agri",
                            "perc.urban",
                            "perc.water",
                            "sqrt.perc.trees",
                            "sqrt.perc.agri",
                            "sqrt.perc.urban",
                            "sqrt.perc.water",
                            "gini.bird.simpson",
                            "lag2gini.bird.simpson",
                            "lag4gini.bird.simpson",
                            "lag6gini.bird.simpson",
                            "lag8gini.bird.simpson",
                            "log.bird.richness"))
saveRDS(w, "data/fewer-covariates.rds")
