
options(stringsAsFactors = FALSE)
library(edgeR)
path = "/home/kp/Documents/Gumed/Implatelet/"
gtfPath = "/home/kp/Documents/Gumed/index/grch38/gencode.v26.chr_patch_hapl_scaff.annotation.gtf"
gtfPath19 = "/home/kp/Documents/Gumed/index/grch19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
dataPath = "data/"
includePath = "include/" 
geneInfoPath = "dgeGenesEnsembl75.RData"  
library(pracma)  
load(paste0(path, dataPath, geneInfoPath))
# shouldBeOc = c("OC_026", "OC_092", "OC_122") 
countsSubset = read.csv(paste0(path, dataPath, "ImPlatelet_counts_raw.tsv"), sep = "\t")#, col.names = T, row.names = T)
sampleInfoAll = read.csv( paste0(path, dataPath, "ImPlatelet_samples.tsv"), sep = "\t")#, col.names = T, row.names = T)
# sampleInfoAll = sampleInfoAll[-which(sampleInfoAll$OriginalGroup == "EC") , ]
# write.table(sampleInfoAll, "/home/kp/Documents/Gumed/Implatelet/data/ImPlatelet_samples.tsv", sep = "\t", col.names = T, row.names = T)
# write.table(countsSubset, "/home/kp/Documents/Gumed/Implatelet/data/ImPlatelet_counts_raw.tsv", sep = "\t", col.names = T, row.names = T)
# countsSubset = countsSubset[, rownames(sampleInfoAll)]
# generate DGE-object, as implemented from the edgeR package 
dge <- DGEList(counts = countsSubset,
               group = sampleInfoAll$OriginalGroup,
               genes = genes
)
dge$samples <- sampleInfoAll#cbind(dge$samples, samplesAll) # add the sample info to the object
#save(dge, file = paste0(path, dataPath, "ovarianBotEcHcGdansk", Sys.Date(),".RData"))

# if(length(which(dge$samples$OriginalGroup == "OC")) > 0 )
# {
#   dge = dge[, -which(dge$samples$OriginalGroup == "OC")]
# }
reportPath = paste(path, "Report/", sep = "") 
# groupNameA = "HC"
# groupNameB = "EC" 
# groupNameC = "OC" 
source(paste(path, "include/statisticalAnalysis.R", sep = "")) 
dataFiltered = normalizeDESeq2NoReport(dge$counts)#, healthyId, ovarianId, reportPath, groupNameA, groupNameB)


 
source(paste0(path, "include/annotateData.R"))
annotation = annotateData(gtfPath19, dataFiltered)
dataFiltered = annotation$dataFiltered#[[1]]
genePositionInfo = annotation$genePositionInfo[[2]]
source(paste0(path, includePath, 'generateKeggPathwayImages.R'))
picWidth = 345
picHeight = 243
diseaseLabel = "OC"
matrixPath = generateKeggPathwayImages(path, dataFiltered)
save.image(paste0(path, "preDataSetConstruction", Sys.Date(), ".RData"))
source(paste0(path, includePath, "prepareDataSets.R"))
seed = 123
k = 5
dataSets = prepareDataSets(dge, dataFiltered, diseaseLabel, matrixPath, picWidth, picHeight, k, seed)
save.image(paste0(path, "preCV.RData"))
source(paste0(path, includePath, "buildModel.R"))
weightsPath = "weights/"
weightsFile = "modelInitialWeights.hdf5"
model = buildModel(path, weightsPath, weightsFile, picWidth, picHeight)
source(paste0(path, includePath, "runCrossValidation.R"))
iters = 190
results = runCrossValidation(model, dge, dataSets, pathToWeights = paste0(path, weightsPath, weightsFile), iters, k, seed)
source(paste0(path, includePath, "generateResultsReport.R"))
subtitle = "OC_classification_results"
generateResultsReport(dge, results, subtitle, reportPath, dataSets$usedSubset, dataSets$otherSubset, dataSets$testId, k, seed)
