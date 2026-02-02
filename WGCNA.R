#路径
setwd("CDNFE")
#环境
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("impute")
# BiocManager::install("WGCNA")
library(impute)
library(WGCNA)
library(reshape2)
library(stringr)

#读取express
library(readxl)
expres = read_excel("Luohe_cesc.xlsx")#前期未处理数据
expres <- data.frame(c(expres[,1:9]))
rownames(expres) <-expres[,1] #更换列名为数据的第一列
expres <-expres[-1]
table(duplicated(expres$`sample`))
#expres=expres[which(rowSums(expres)>200),] #筛选，令每行表达值至少大于1(不用)
expres=expres[order(rowSums(expres),decreasing = T),] #对总表达值进行降序：按每行总和进行降序
#标准化处理20210214 

#样本和基因的筛选
expres[]=lapply(expres, as.numeric) #将int转为numeric
#expres=expres[1:5000,-c(1:50)]  #取前20000个基因，前50个样本
expres=expres[1:1000,]
#expres=expres[1:50]
datExpr=t(expres) #对expres进行转置，以基因作为变量
#datExpr=as.data.frame(t(expres[,-c(1)]))
#names(datExpr)=expres$sample;
#rownames(datExpr)=names(expres)[-c(1)];

#基因缺失值检测
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
# if (!gsg$allOK)
# {
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0)
#     printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0)
#     printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
# }
# dim(datExpr)
#%重复该段前两行

# ##样本离群值检测(不执行)
# sampleTree = hclust(dist(datExpr),method = "average")
# sizeGrWindow(12,9)
# par(cex = 0.6)
# par(mar = c(0,4,2,0))
# abline(h = 80000, col = "red");
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# abline(h = 200, col = "red");
# clust = cutreeStatic(sampleTree, cutHeight = 80000, minSize = 10)
# table(clust)
# keepSamples = (clust==1)
# datExpr = datExpr[keepSamples, ]
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)
# save(datExpr, file = "AS-green-FPKM-01-dataInput.RData")


##选合适的阈值 Scale independence and Mean Connectivity
powers = c(c(1:5),seq(from =6,to=40, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose =5)
sizeGrWindow(9,5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.20,col="red")  #this line corresponds to using an R^2 cut-off of h
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="n",
     main = paste("Mean Connectivity")) #Mean connectivity as a function of the soft-thresholding power
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers, cex=cex1,col="red")
power_select=sft[["powerEstimate"]]
power=power_select
power=6
# #软阈值的检验(此处的beta1等价于上面的power)
# beta1=3#根据自己的检验结果来选择
# Connectivity=softConnectivity(datExpr,power=beta1)
# par(mfrow=c(1,1))
# scaleFreePlot(Connectivity, main=paste("soft threshold, power=",beta1), truncated=T)

##网络构建
#一步法
dim(datExpr)
#参数
type = "unsigned"
corType = "pearson"  #pearson or bicor
corFnc = ifelse(corType=="pearson",cor,bicor)
maxPOutliers = ifelse(corType=='pearson',0.05)
robustY = ifelse(corType=="pearson",T,F)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#net = blockwiseModules(datExpr, power = power, maxBlockSize = 10000,
#                       TOMType = "unsigned", minModuleSize = 30,
#                       reassignThreshold = 0, mergeCutHeight = 0.25,
#                       numericLabels = FALSE, pamRespectsDendro = FALSE,
#                       saveTOMs = TRUE,
#                       saveTOMFileBase = "AS-green-FPKM-TOM",
#                       verbose = 3)
net = blockwiseModules(datExpr, power = power, maxBlockSize = 10000,
                       TOMType = type, minModuleSize = 5,
                       reassignThreshold = 0,mergeCutHeight = 0.25,
                       numericLabels = FALSE,pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,corType = corType,
                       maxPOutliers = maxPOutliers, loadTOMs = TRUE,
                       #saveTOMFileBase = paste0("AS-green-FPKM-TOM"),
                       verbose = 3)

#save(net,datExpr,file='./WGCNA.Rdata')#20210217
#模块聚类 Cluster Dendrogram
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]],moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#下句运行不出来
# datME=moduleEigengenes(datExpr[,Connectivity],moduleColors)[[1]]
# modul<-signif(cor(datExpr, use="p"), 2)
# write.table(moduleColors ,"modul_cor3.txt",sep="\t",quote=F)
# write.table(moduleLabels ,"moduleLabels.txt",sep="\t",quote=F)                                

#构建邻接矩阵并导出edge.txt
ADJ= adjacency(datExpr,power = power)
#write.table(ADJ,"ADJ-5000.txt",sep="\t",quote=F)
#probes = names(datExpr)
#modules = c("0","1","2","3","4",
#            "5","6","7","8","9",
#            "10","11","12","13","14",
#            "15","16","17","18","19",)
#inModule = is.finite(match(moduleColors, modules));
#modProbes = probes[inModule];
# vis_1= exportNetworkToCytoscape(ADJ,threshold = 0,
#                                 edgeFile="edge_ADJ.txt",nodeFile="node_ADJ.txt",
#                                 weighted = TRUE);
# vis= exportNetworkToCytoscape(ADJ,threshold = 0.001,
#                               edgeFile="edge_ADJ.txt",nodeFile="node_ADJ.txt",
#                               weighted = TRUE);
#nodeNames = modProbes,nodeAttr = moduleColors[inModule]

#绘制邻接矩阵热力图（热图，箱线图）(运行不出来)
#module eigengene,可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("MEs",labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_T = orderMEs(MEs_col)
#根据基因间表达量进行聚类所得到的各模块间的相关性图
#基因模块之间的相关性
#marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs,"Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap= c(3,4,2,2),plotDendrograms = T,
                      xLabelsAngle = 90)


#TOM PLOT拓扑图（热图）
load(net$TOMFiles[1],verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
#Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

##导出网络图需要的数据，要用Cytoscape绘制
probes = colnames(datExpr)
dimnames(TOM) <- list(probes, probes)
#Export the network into edge anshold = 0,d node list files Cytoscape can read
#threshold 默认为0.5，可以根据需要调整，也可以导出后在Cytoscape中调整
#阈值越大数据越少
cyt = exportNetworkToCytoscape(TOM,threshold =0.0000068,
                               nodeNames = probes,nodeAttr = moduleColors)
cyt$edgeData
#cyt$nodeData
write.table(cyt$edgeData,"edge_TOM.xls",sep="\t",quote=F)
#write.table(cyt$nodeData,"node_TOM.xls",sep="\t",quote=F)
# cyt1= exportNetworkToCytoscape(TOM,threshold = 0.01,
#                                edgeFile="edge_TOM.txt",nodeFile="node_TOM.txt",
#                                weighted = TRUE)

