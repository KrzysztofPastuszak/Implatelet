
createDirIfNotExists = function(mainDir, subDir){
  if (file.exists(paste(mainDir,  subDir, sep = "")) == F){
    
    dir.create(paste(mainDir, subDir, sep = ""))
  }  
} 
buildExpressionPathwayMatrix = function(sigSym, metSym, cancerSym, data, col, picWidth = 345, picHeight = 243)
{
  output = matrix(0, ncol = picWidth, nrow = picHeight)
  k = 1 
  for(i in 1:length(sigSym))
  {
    symbols = unlist(sigSym[i])
    nameMatch = match(symbols, rownames(dataFiltered))
    nonMissing = which(is.na(nameMatch) == F)
    symbols = symbols[nonMissing] 
    output[k, 1:length(nonMissing)] = data[nameMatch[nonMissing], col]
    k = k+1 
  } 
  for(i in 1:length(metSym))
  {
    symbols = unlist(metSym[i]) 
    nameMatch = match(symbols, rownames(dataFiltered))
    nonMissing = which(is.na(nameMatch) == F)
    symbols = symbols[nonMissing] 
    output[k, 1:length(nonMissing)] = data[nameMatch[nonMissing], col]
    k = k+1 
  } 
  for(i in 1:length(cancerSym))
  {
    symbols = unlist(cancerSym[i])
    nameMatch = match(symbols, rownames(dataFiltered))
    nonMissing = which(is.na(nameMatch) == F)
    symbols = symbols[nonMissing] 
    output[k, 1:length(nonMissing)] = data[nameMatch[nonMissing], col]
    k = k+1 
  } 
  return(output)
}

generateKeggPathwayImages = function(path, dataFiltered, pathwaysPath = "KeggImages")
{
  library(KEGGREST)
  library(gage)
  library(org.Hs.eg.db)
  data(egSymb)
  kg.hsa=kegg.gsets(species = "hsa", id.type = "kegg")
  dis.kg = kg.hsa$kg.sets[kg.hsa$dise.idx]
  met.kg = kg.hsa$kg.sets[kg.hsa$met.idx]
  sig.kg = kg.hsa$kg.sets[kg.hsa$sig.idx]
  cancerKeggId = dis.kg[grep("hsa052", names(dis.kg))] 
  cancerKeggId = cancerKeggId[-grep("hsa05200", names(cancerKeggId))]
  metabolismSpecificKeggId = met.kg[grep("hsa0", names(met.kg))]
  #metabolismSpecificKeggId = met.kg[grep("hsa00", names(met.kg))]
  metabolismGeneralKeggId = met.kg[c(
    grep("hsa011", names(met.kg)),
    grep("hsa012", names(met.kg)))]
  metabolismSpecificKeggId = metabolismSpecificKeggId[-grep("hsa01100", names(metabolismSpecificKeggId))]
  signalKeggId = sig.kg[c(
    grep("hsa030", names(sig.kg)),
    grep("hsa00970", names(sig.kg)),
    grep("hsa030", names(sig.kg)),
    grep("hsa041", names(sig.kg)),
    grep("hsa034", names(sig.kg)),
    grep("hsa040", names(sig.kg)),
    grep("hsa0421", names(sig.kg)),
    grep("hsa0491", names(sig.kg)),
    grep("hsa0492", names(sig.kg)),
    grep("hsa04935", names(sig.kg)),
    grep("hsa046", names(sig.kg)),
    grep("hsa03320", names(sig.kg)),
    grep("hsa045", names(sig.kg)),
    grep("hsa04810", names(sig.kg)),
    grep("hsa040", names(sig.kg)))]
  metabolismGeneralKeggId = metabolismGeneralKeggId[-grep("hsa01100", names(metabolismGeneralKeggId))] 
  cancerSym = lapply(cancerKeggId, eg2sym)
  sigSym = lapply(signalKeggId , eg2sym)
  metSym = lapply(metabolismSpecificKeggId , eg2sym)
  
  
  picWidth = 345
  picHeight =  243
  
  createDirIfNotExists(path, "KEGG_Pathway_Image")
  pathwayPath = paste0(path, "KEGG_Pathway_Image", "/") 
  createDirIfNotExists(pathwayPath, "Matrices")
  matrixPath = paste0(pathwayPath, "Matrices", "/") 
  createDirIfNotExists(pathwayPath, "Images")
  imagePath = paste0(pathwayPath, "Images", "/") 
  
  for(col in 1:ncol(dataFiltered))
  {
    output = buildExpressionPathwayMatrix(sigSym, metSym, cancerSym, data = dataFiltered, col)  
    write.table(output, file = paste0(matrixPath, colnames(dataFiltered)[col], ".txt"), row.names = F, quote = F, col.names = F) 
    png(paste0(imagePath, colnames(dataFiltered)[col], ".png"), width = picWidth, height = picHeight) 
    par(mar = c(0,0,0,0))
    require(grDevices) # for colours
    #image(t(output), xlab = NULL, ylab = NULL)
    image(t(output), col=colorRampPalette(c( "black",  "red",  "red"))(25),  axes=FALSE, mar = c(0,0,0,0),ann  = F)
    dev.off()
  }
  
  
  
  print("KEGG pathway images prepared")
  return(matrixPath)
}