library(diveRsity)
library(qgraph)

mig <- divMigrate(infile = "microsatellites_genepop.gen",
                stat = "all", plot_network = FALSE)

coords=read.table("coordinates.txt", header=T)

connectivity=mig$gRelMig

qgraph(connectivity, 
  layout = as.matrix(coords[,c("Long", "Lat")]),
  mar = c(3, 3, 2, 2),
  vsize = rep(8, nrow(connectivity)),
  label.cex = 1.5,
  label.scale = FALSE,
  label.scale.equal = TRUE,
  edge.labels = FALSE,
  labels = coords$Locality,
  normalize = TRUE,
  posCol = "black",
  border.color = NA,
  color = "transparent"
)

write.table(connectivity, "connectivity.txt", col.names = F, row.names = F, quote=F)
