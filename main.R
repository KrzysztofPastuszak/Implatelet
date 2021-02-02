
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
countsSubset = read.csv(paste0(path, dataPath, "ImPlatelet_counts_raw.tsv"), sep = "\t") 
sampleInfoAll = read.csv( paste0(path, dataPath, "ImPlatelet_samples.tsv"), sep = "\t") 
dge <- DGEList(counts = countsSubset,
               group = sampleInfoAll$OriginalGroup,
               genes = genes
)
dge$samples <- sampleInfoAll  
reportPath = paste(path, "Report/", sep = "")  
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
dataSets = prepareDataSets(dge, dataFiltered, diseaseLabel, matrixPath, picWidth, picHeight, seed)
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
