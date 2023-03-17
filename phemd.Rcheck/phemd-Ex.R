pkgname <- "phemd"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('phemd')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GDM")
### * GDM

flush(stderr()); flush(stdout())

### Name: GDM
### Title: Accessor function for EMD ground distance matrix
### Aliases: GDM

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
gdm <- GDM(phemdObj)




cleanEx()
nameEx("Phemd-methods")
### * Phemd-methods

flush(stderr()); flush(stdout())

### Name: Phemd-methods
### Title: Setter function for protein / gene markers
### Aliases: Phemd-methods selectMarkers<- selectMarkers<-,Phemd-method
###   Phemd,ANY,ANY-method rawExpn<- Phemd,character,ANY-method
###   rawExpn<-,Phemd-method pooledCells<- pooledCells<-,Phemd-method
###   subsampledIdx<- subsampledIdx<-,Phemd-method subsampledBool<-
###   subsampledBool<-,Phemd-method monocleInfo<-
###   monocleInfo<-,Phemd-method seuratInfo<- seuratInfo<-,Phemd-method
###   phateInfo<- phateInfo<-,Phemd-method celltypeFreqs<-
###   celltypeFreqs<-,Phemd-method batchIDs<- batchIDs<-,Phemd-method GDM<-
###   GDM<-,Phemd-method

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
new_genes <- all_genes
new_genes[1] <- 'IL2R'
selectMarkers(phemdObj) <- new_genes

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
new_expn_data <- all_expn_data
new_expn_data <- lapply(new_expn_data, function(x) {log2(x+1)})
rawExpn(phemdObj) <- new_expn_data

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
aggregated_data <- t(do.call(rbind,all_expn_data))
pooledCells(phemdObj) <- aggregated_data

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
subsampledIdxList<- rep(list(1:10), length(all_expn_data)) #subsampled cells 1-10 from each sample
subsampledIdx(phemdObj) <- subsampledIdxList

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
subsampledBool(phemdObj) <- TRUE

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
mydata <- pooledCells(phemdObj)
myCellDataSet <- newCellDataSet(mydata,phenoData=NULL, expressionFamily=VGAM::negbinomial.size())
monocleInfo(phemdObj) <- myCellDataSet

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_seuratObj <- Seurat::CreateSeuratObject(counts = t(all_expn_data[[1]]), project = "A")
seuratInfo(phemdObj) <- my_seuratObj

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#my_phateObj <- phateR::phate(all_expn_data[[1]])
phateInfo(phemdObj) <- list()

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
myCellTypeFreqs <- matrix(rexp(length(all_expn_data)*10, rate=.1), ncol=10)
myCellTypeFreqs <- apply(myCellTypeFreqs, 1, function(x) {x / sum(x)})
celltypeFreqs(phemdObj) <- myCellTypeFreqs

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_seuratObj <- Seurat::CreateSeuratObject(counts = t(all_expn_data[[1]]), project = "A")
seuratInfo(phemdObj) <- my_seuratObj
batchIDs(phemdObj) <- rep('A', length(all_expn_data))

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
cluster_locs <- 1:10
myGDM <- as.matrix(dist(cluster_locs))
GDM(phemdObj) <- myGDM




cleanEx()
nameEx("aggregateSamples")
### * aggregateSamples

flush(stderr()); flush(stdout())

### Name: aggregateSamples
### Title: Aggregate expression data from all samples
### Aliases: aggregateSamples

### ** Examples

my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)




cleanEx()
nameEx("assignCellClusterNearestNode")
### * assignCellClusterNearestNode

flush(stderr()); flush(stdout())

### Name: assignCellClusterNearestNode
### Title: Assign cells to a reference cell subtype
### Aliases: assignCellClusterNearestNode

### ** Examples

## Not run: 
##D cur_cells_cluster_labels <- assignCellClusterNearestNode(cur_cells_expn_data,
##D clustered_cells_expn_data, clustered_cells_cluster_labels, cell_model='monocle2')
## End(Not run)



cleanEx()
nameEx("batchIDs")
### * batchIDs

flush(stderr()); flush(stdout())

### Name: batchIDs
### Title: Accessor function for batch ID for each sample
### Aliases: batchIDs

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
batch_metadata <- batchIDs(phemdObj)




cleanEx()
nameEx("bindSeuratObj")
### * bindSeuratObj

flush(stderr()); flush(stdout())

### Name: bindSeuratObj
### Title: Attach 'Seurat' object to 'Phemd' object
### Aliases: bindSeuratObj

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_seuratObj <- Seurat::CreateSeuratObject(counts = t(all_expn_data[[1]]), project = "A")
my_seuratObj <- Seurat::FindVariableFeatures(object = my_seuratObj)
my_seuratObj <- Seurat::ScaleData(object = my_seuratObj, do.scale=FALSE, do.center=FALSE)
my_seuratObj <- Seurat::RunPCA(object = my_seuratObj, pc.genes = colnames(all_expn_data[[1]]), do.print = FALSE)
my_seuratObj <- Seurat::FindNeighbors(my_seuratObj, reduction = "pca", dims.use = 1:10)
my_seuratObj <- Seurat::FindClusters(my_seuratObj, resolution = 0.6, print.output = 0, save.SNN = TRUE)
my_phemdObj <- bindSeuratObj(my_phemdObj, my_seuratObj)




cleanEx()
nameEx("celltypeFreqs")
### * celltypeFreqs

flush(stderr()); flush(stdout())

### Name: celltypeFreqs
### Title: Accessor function for cell subtype distribution for each sample
### Aliases: celltypeFreqs

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
celltype_weights <- celltypeFreqs(phemdObj)




cleanEx()
nameEx("clusterIndividualSamples")
### * clusterIndividualSamples

flush(stderr()); flush(stdout())

### Name: clusterIndividualSamples
### Title: Computes cell subtype abundances for each sample
### Aliases: clusterIndividualSamples

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)




cleanEx()
nameEx("compareSamples")
### * compareSamples

flush(stderr()); flush(stdout())

### Name: compareSamples
### Title: Computes EMD distance matrix representing pairwise dissimilarity
###   between samples
### Aliases: compareSamples

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)
my_EMD_mat <- compareSamples(my_phemdObj_final)




cleanEx()
nameEx("createDataObj")
### * createDataObj

flush(stderr()); flush(stdout())

### Name: createDataObj
### Title: Create 'Phemd' object
### Aliases: createDataObj

### ** Examples

my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))




cleanEx()
nameEx("drawColnames45")
### * drawColnames45

flush(stderr()); flush(stdout())

### Name: drawColnames45
### Title: Rotates heatmap marker labels 45 degrees
### Aliases: drawColnames45

### ** Examples

#Not to be called directly



cleanEx()
nameEx("embedCells")
### * embedCells

flush(stderr()); flush(stdout())

### Name: embedCells
### Title: Generate cell-state embedding
### Aliases: embedCells

### ** Examples

my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_lg <- embedCells(my_phemdObj_lg, cell_model='monocle2', data_model = 'gaussianff', sigma=0.02, maxIter=2)



cleanEx()
nameEx("generateGDM")
### * generateGDM

flush(stderr()); flush(stdout())

### Name: generateGDM
### Title: Computes ground distance matrix based on cell embedding
### Aliases: generateGDM

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)




cleanEx()
nameEx("getArithmeticCentroids")
### * getArithmeticCentroids

flush(stderr()); flush(stdout())

### Name: getArithmeticCentroids
### Title: Get arithmetic centroids (coordinates)
### Aliases: getArithmeticCentroids

### ** Examples

## Not run: 
##D cluster_centroids <- getArithmeticCentroids(ref_clusters)
## End(Not run)



cleanEx()
nameEx("getCellYield")
### * getCellYield

flush(stderr()); flush(stdout())

### Name: getCellYield
### Title: Gets cell yield of each sample as a table
### Aliases: getCellYield

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)
my_EMD_mat <- compareSamples(my_phemdObj_final)
cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
getCellYield(my_phemdObj_final, cluster_assignments)




cleanEx()
nameEx("getSampleCelltypeFreqs")
### * getSampleCelltypeFreqs

flush(stderr()); flush(stdout())

### Name: getSampleCelltypeFreqs
### Title: Returns cell subtype distribution for each sample as a table
### Aliases: getSampleCelltypeFreqs

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)
my_EMD_mat <- compareSamples(my_phemdObj_final)
cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
getSampleCelltypeFreqs(my_phemdObj_final, cluster_assignments)




cleanEx()
nameEx("getSampleHistsByCluster")
### * getSampleHistsByCluster

flush(stderr()); flush(stdout())

### Name: getSampleHistsByCluster
### Title: Gets cell subtype frequency histograms for each sample by
###   cluster ID
### Aliases: getSampleHistsByCluster

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)
my_EMD_mat <- compareSamples(my_phemdObj_final)
cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
weights_by_cluster <- getSampleHistsByCluster(my_phemdObj_final, cluster_assignments)




cleanEx()
nameEx("getSampleSizes")
### * getSampleSizes

flush(stderr()); flush(stdout())

### Name: getSampleSizes
### Title: Retrieve single-cell sample sizes
### Aliases: getSampleSizes

### ** Examples

## Not run: 
##D sample_sizes <- getSampleSizes(all_expn_data)
## End(Not run)



cleanEx()
nameEx("groupSamples")
### * groupSamples

flush(stderr()); flush(stdout())

### Name: groupSamples
### Title: Performs community detection on sample-sample distance matrix to
###   identify groups of similar samples
### Aliases: groupSamples

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, cell_model = 'monocle2', data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)
my_EMD_mat <- compareSamples(my_phemdObj_final)
cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)




cleanEx()
nameEx("identifyCentroids")
### * identifyCentroids

flush(stderr()); flush(stdout())

### Name: identifyCentroids
### Title: Identify cluster centroids (cell names)
### Aliases: identifyCentroids

### ** Examples

## Not run: 
##D centroid_names <- identifyCentroids(ref_clusters)
## End(Not run)



cleanEx()
nameEx("monocleInfo")
### * monocleInfo

flush(stderr()); flush(stdout())

### Name: monocleInfo
### Title: Accessor function for stored Monocle object
### Aliases: monocleInfo

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
monocle_obj <- monocleInfo(phemdObj)




cleanEx()
nameEx("orderCellsMonocle")
### * orderCellsMonocle

flush(stderr()); flush(stdout())

### Name: orderCellsMonocle
### Title: Compute Monocle2 cell state and pseudotime assignments
### Aliases: orderCellsMonocle

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, cell_model='monocle2', data_model='gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)



cleanEx()
nameEx("phateInfo")
### * phateInfo

flush(stderr()); flush(stdout())

### Name: phateInfo
### Title: Accessor function for stored phate object
### Aliases: phateInfo

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
phateobj <- phateInfo(phemdObj)




cleanEx()
nameEx("plotCellYield")
### * plotCellYield

flush(stderr()); flush(stdout())

### Name: plotCellYield
### Title: Plot cell yield of each sample as bar plot
### Aliases: plotCellYield

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)
my_EMD_mat <- compareSamples(my_phemdObj_final)
cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
plotCellYield(my_phemdObj_final, labels=cluster_assignments, font_sz = 0.8)




cleanEx()
nameEx("plotEmbeddings")
### * plotEmbeddings

flush(stderr()); flush(stdout())

### Name: plotEmbeddings
### Title: Plots Monocle2 cell embedding plots
### Aliases: plotEmbeddings

### ** Examples

my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model='gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
cmap <- plotEmbeddings(my_phemdObj_monocle)



cleanEx()
nameEx("plotGroupedSamplesDmap")
### * plotGroupedSamplesDmap

flush(stderr()); flush(stdout())

### Name: plotGroupedSamplesDmap
### Title: Plot diffusion map embedding of samples based on distance matrix
### Aliases: plotGroupedSamplesDmap

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)
my_EMD_mat <- compareSamples(my_phemdObj_final)
cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
printClusterAssignments(cluster_assignments, my_phemdObj_final, '.', overwrite=TRUE)
dm <- plotGroupedSamplesDmap(my_EMD_mat, cluster_assignments, pt_sz=2)




cleanEx()
nameEx("plotHeatmaps")
### * plotHeatmaps

flush(stderr()); flush(stdout())

### Name: plotHeatmaps
### Title: Plot heatmap of cell subtypes
### Aliases: plotHeatmaps

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_lg <- selectFeatures(my_phemdObj_lg, selected_genes)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff',
pseudo_expr=0, sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
myheatmap <- plotHeatmaps(my_phemdObj_monocle, cell_model='monocle2')




cleanEx()
nameEx("pooledCells")
### * pooledCells

flush(stderr()); flush(stdout())

### Name: pooledCells
### Title: Accessor function for aggregated cells used for cell subtype
###   definition
### Aliases: pooledCells

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
cells_aggregated <- pooledCells(phemdObj)




cleanEx()
nameEx("printClusterAssignments")
### * printClusterAssignments

flush(stderr()); flush(stdout())

### Name: printClusterAssignments
### Title: Writes samples to file based on community detection group
###   assignments
### Aliases: printClusterAssignments

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
my_phemdObj_final <- generateGDM(my_phemdObj_final)
my_EMD_mat <- compareSamples(my_phemdObj_final)
cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
printClusterAssignments(cluster_assignments, my_phemdObj_final, '.', overwrite=TRUE)




cleanEx()
nameEx("rawExpn")
### * rawExpn

flush(stderr()); flush(stdout())

### Name: rawExpn
### Title: Accessor function for stored multi-sample raw expression data
### Aliases: rawExpn

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
raw_expn_data <- rawExpn(phemdObj)




cleanEx()
nameEx("removeTinySamples")
### * removeTinySamples

flush(stderr()); flush(stdout())

### Name: removeTinySamples
### Title: Remove samples with too few cells
### Aliases: removeTinySamples

### ** Examples

my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10) #removes samples with fewer than 10 cells




cleanEx()
nameEx("retrieveRefClusters")
### * retrieveRefClusters

flush(stderr()); flush(stdout())

### Name: retrieveRefClusters
### Title: Retrieve reference cell clusters
### Aliases: retrieveRefClusters

### ** Examples

## Not run: 
##D cluster_expression_data <- retrieveRefClusters(my_phemdObj)
## End(Not run)




cleanEx()
nameEx("sNames")
### * sNames

flush(stderr()); flush(stdout())

### Name: sNames
### Title: Accessor function for identifiers of all single-cell samples in
###   experiment
### Aliases: sNames

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
sampleIDs <- sNames(phemdObj)




cleanEx()
nameEx("selectFeatures")
### * selectFeatures

flush(stderr()); flush(stdout())

### Name: selectFeatures
### Title: Perform feature selection on aggregated data
### Aliases: selectFeatures

### ** Examples


my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_lg <- selectFeatures(my_phemdObj_lg, selected_genes=c('TP53',
'EGFR', 'KRAS', 'FOXP3', 'LAG3'))




cleanEx()
nameEx("selectMarkers")
### * selectMarkers

flush(stderr()); flush(stdout())

### Name: selectMarkers
### Title: Accessor function for gene/protein markers measured in
###   experiment
### Aliases: selectMarkers

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
genes <- selectMarkers(phemdObj)




cleanEx()
nameEx("seuratInfo")
### * seuratInfo

flush(stderr()); flush(stdout())

### Name: seuratInfo
### Title: Accessor function for stored Seurat object within Phemd object
### Aliases: seuratInfo

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
seurat_obj <- seuratInfo(phemdObj)




cleanEx()
nameEx("subsampledBool")
### * subsampledBool

flush(stderr()); flush(stdout())

### Name: subsampledBool
### Title: Accessor function for whether or not cells were subsampled when
###   aggregated for cell subtype analysis
### Aliases: subsampledBool

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
subsampled <- subsampledBool(phemdObj)




cleanEx()
nameEx("subsampledIdx")
### * subsampledIdx

flush(stderr()); flush(stdout())

### Name: subsampledIdx
### Title: Accessor function for aggregated cells used for cell subtype
###   definition
### Aliases: subsampledIdx

### ** Examples

phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
subsampled_idx_list <- subsampledIdx(phemdObj)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
