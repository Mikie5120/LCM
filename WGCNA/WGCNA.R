#library
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
#enableWGCNAThreads()
# enableWGCNAThreads(nThreads=4)

LCM <- read.table(file.choose(),head=T, row.name=1, sep="\t")
samples=read.csv('Sam_info.txt',sep = '\t',row.names = 1)

##25%##
m.vars=apply(LCM,1,var)
LCM.upper=LCM[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
dim(LCM.upper)
datExpr=as.data.frame(t(LCM.upper));
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#all
datExpr = as.data.frame(t(LCM))
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##data<-log10(date[,-1]+0.01)##
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")

#sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
tiff('sampleTree.tiff')
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()
save(datExpr, file = "AS-green-FPKM-01-dataInput.RData")


#cutHeight = 20000

clust = cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
save(datExpr, file = "FPKM-01-dataInput.RData")

## Soft threshold##
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.80;

# Scale-free topology fit index as a function of the soft-thresholding power
tiff('Scale_independence_8.tiff')
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power
tiff('Mean_connectivity_8.tiff')
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft$powerEstimate

##One-step network construction and module detection##

net = blockwiseModules(datExpr, power = 13, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "LCM-TOM",
                       verbose = 3)

#Check how many genes in each module. 0 means the genes out side these modules.
table(net$colors)

### open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors (net$colors)
# Plot the dendrogram and the module colors underneath
tiff('Module_colors13.tiff')
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

##module eigengene###
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
write(moduleColors, file="LCM_module.txt");
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "AS-green-FPKM-02-networkConstruction-auto.RData")

#output to Cytoscape#### Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 13);
# Read in the annotation file# annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("green");
# Select module probes#

probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-edges_15-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-nodes_15-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.3,
                               nodeNames = modProbes,                               
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);


##Module-trait_relationships##
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
tiff('Module-trait_relationships13.tiff')
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.9,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.7, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()



# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 13);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot

diag(plotTOM) = NA;
# Call the plot function#sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#???便選???1000???基??????可視???
nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM),method = "average")
selectColors = moduleColors[select];

# Open a graphical window#sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
tiff('LCM_Network_heatmap_plot.tiff')
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

# marDendro/marHeatmap 
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
sizeGrWindow(7, 6) 
tiff('LCM_adjacency_heatmap.tiff')
plotEigengeneNetworks(MET, "LCM adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()










