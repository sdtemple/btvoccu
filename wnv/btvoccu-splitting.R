# btvoccu Splitting
# Seth Temple, sdtemple@lanl.gov
# September 10, 2021

library(btvoccu)

y <- readRDS("data/response.rds")
y <- apply(y, 1:4, as.numeric)

ytest <- subset_4darray(y, 2, c(2008, 2012, 2016)) # good, bad, okay years
saveRDS(ytest, "data/response-test.rds")

ysplit <- subset_4darray(y, 2, c(2002:2007,2009:2011,2013:2015,2017))
saveRDS(ysplit, "data/response-split.rds")

set.seed(9102021)
ysplit <- split_4darray(ysplit, .75, periods = 18:44)

# check that there isn't in balance in north ontario presences
# seems good enough
sum(ysplit$train, na.rm = T) / sum(!is.na(ysplit$train))
sum(ysplit$test, na.rm = T) / sum(!is.na(ysplit$test))
z1 <- subset_4darray(ysplit$train, 1, c("ALG", "NPS", "NWR", "PQP", "SUD", "THB", "TSK"))
z2 <- subset_4darray(ysplit$test, 1, c("ALG", "NPS", "NWR", "PQP", "SUD", "THB", "TSK"))
sum(z1, na.rm = T) / sum(!is.na(z1))
sum(z2, na.rm = T) / sum(!is.na(z2))

saveRDS(ysplit$train, "data/response-train.rds")
saveRDS(ysplit$test, "data/response-validate.rds")