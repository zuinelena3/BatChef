library(splatter)
library(scuttle)
library(scater)

params <- newSplatParams(nGenes = 1000, batchCells = c(155, 150),
                         group.prob = c(0.3, 0.5, 0.2),seed = 33)
data <- splatSimulateGroups(params, verbose = FALSE)
data <- logNormCounts(data)
data <- runPCA(data,  ncomponent = 10)
usethis::use_data(data)

