# Barplots for btvoccu paper
# Seth Temple, sdtemple@lanl.gov
# September 6, 2021

# response data
y <- readRDS("data/response.rds")
y <- apply(y, 1:4, as.numeric) 

# site bar plot
sy <- apply(y, 1, sum, na.rm = T)
barplot(sy, names.arg = names(sy), cex.names = 0.75, ylab = "Positive Cases", xlab = "Site", las = 2)

# season bar plot
yy <- apply(y, 2, sum, na.rm = T)
barplot(yy, names.arg = names(yy), cex.names = 0.75, ylab = "Positive Cases", xlab = "Season", las = 2)
