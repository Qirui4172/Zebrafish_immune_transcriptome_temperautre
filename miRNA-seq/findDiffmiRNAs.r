#!/usr/bin/env Rscript


#=============================================================================================
Usage<-function(){
	cat("Usage: Rscript findDiffmiRNAs.r [readcounts] [sampleinfo] [removeSamp] [compaType] [minNonzero] [minRowsum]\n\n",

	"Parameters:\n",
	"[readcounts]	Read counts matrix file (mirna-seq_readcounts.mx)\n",
	"[sampleinfo]	Sample information file (sampleinfo_diff.txt)\n",
	"[removeSamp]	Remove which sample, i.e. t24ct4, no (if not removing any sample)\n",
	"[compaType]	Comparison type, LPS or Temp\n",
	"[minNonzero]	For one miRNA in 5 replicates, at least this number out of 5 replciates should have nonzero read counts, otherwise this miRNA will be dropped (suggested 3).\n",
	"[minRowsum]	For one miRNA in 10 replicates (five 24/32C + five 28C), at least this number of read counts should be met, otherwise this miRNA will be dropped (suggested 100).\n\n",

	"Example:\n",
	"Rscript3.5.1 findDiffmiRNAs.r mirna-seq_readcounts.mx sampleinfo_diff.txt t24ct4 LPS 3 100\n",
	"Rscript3.5.1 findDiffmiRNAs.r mirna-seq_readcounts.mx sampleinfo_diff.txt no Temp 3 100\n\n",

	"Function:\n",
	"Use DESeq2 to find differentially expressed miRNAs and generate plots. If compaType is LPS, LPSvsCT will be performed in each temperature group sequentially; if compaType is Temp, ct32vs28 and ct24vs28 will be performed.\n",
	"Contact: Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"Date: 25-08-2018\n",
	"Updated: 04-01-2021\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=6){Usage();quit();}


cat("==========================================================================================\n")
# Step1. Load libraries
cat("Loading libraries...\n\n")
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggbiplot)
library(biomaRt)


cat("==========================================================================================\n")
# Step2. Read data and parameters
cat("Reading data and parameters...\n\n")
readcounts<-read.table(args[1], header=T, stringsAsFactors=F)
sampleinfo.origin<-read.table(args[2], header=T, stringsAsFactors=F)
sampleinfo<-sampleinfo.origin
sampleinfo$Treatment=as.factor(sampleinfo$Treatment)
sampleinfo$Temperature=as.factor(sampleinfo$Temperature)

rmsample<-args[3]
comptype<-args[4]

if(comptype=="LPS"){comp=paste("LPSvsCT", sep="")}else{comp=paste("ct2432vs28", sep="")}
out.plot<-paste("Plots_", comp, "_mirna.pdf", sep="")
pdf(out.plot)


cat("==========================================================================================\n")
# Generate metadata
cat("Generating metadata...\n\n")
generateMetadata<-function(comptype, readcounts, sampleinfo){
	if(comptype=="LPS"){
		dds=DESeqDataSetFromMatrix(countData=readcounts, colData=sampleinfo, design=~Treatment)
		dds$Treatment=relevel(dds$Treatment, "Control")
	}else if(comptype=="Temp"){
		dds=DESeqDataSetFromMatrix(countData=readcounts, colData=sampleinfo, design=~Temperature)
		dds$Temperature=relevel(dds$Temperature, "28")
	}
	return(dds)
}

dds<-generateMetadata(comptype, readcounts, sampleinfo)
dds<-dds[rowSums(counts(dds))>=100]
vst<-varianceStabilizingTransformation(dds, blind=FALSE)


cat("==========================================================================================\n")
# PCA plot
cat("Ploting PCA...\n\n")
# use DESeq2:plotPCA function
PCA_deseq2<-function(vst){
    data<-plotPCA(vst, intgroup=c("Temperature", "Treatment"), returnData=TRUE)
    percentVar<-round(100*attr(data, "percentVar"))
    pca.plot<-ggplot(data, aes(PC1, PC2, color=Temperature, shape=Treatment))+geom_point(size=3)+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(size=1),text=element_text(size=18),axis.text=element_text(size=15),legend.text=element_text(size=15))+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))
    out=list(pca.plot, data)
    return(out)
}

pca.result1<-PCA_deseq2(vst)
pca.plot1<-pca.result1[[1]]
data1<-pca.result1[[2]]
pca.plot1+coord_fixed(ratio=5/3)
pca.plot1+coord_fixed(ratio=5/3)+geom_text(aes(label=rownames(data1)), size=1, hjust=-0.4)
pca.plot1+coord_fixed(ratio=5/3)+geom_text(aes(label=rownames(data1)), size=1, vjust=-1.8)

# use in-house script
PCA_inhouse<-function(vst, sampleinfo){
    pca.results<-prcomp(t(as.data.frame(assay(vst))), scale=TRUE)
    percentVar<-round(((pca.results$sdev^2)/sum(pca.results$sdev^2))*100, 1)
    data<-as.data.frame(pca.results$x)
    data$Temperature<-sampleinfo$Temperature
    data$Treatment<-sampleinfo$Treatment
    pca.plot<-ggplot(data, aes(PC1, PC2, colour=Temperature, shape=Treatment))+geom_point(size=3)+xlab(paste("PC1: ", percentVar[1], "% variance", sep=""))+ylab(paste("PC2: ", percentVar[2], "% variance", sep=""))+theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_rect(colour="black",size=1), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15))+coord_fixed(ratio=5/4)
	out=list(pca.plot, data)
    return(out)
}

pca.result2<-PCA_inhouse(vst, sampleinfo)
pca.plot2<-pca.result2[[1]]
data2<-pca.result2[[2]]
pca.plot2
pca.plot2+geom_text(aes(label=rownames(data2)), size=1, hjust=-0.4)
pca.plot2+geom_text(aes(label=rownames(data2)), size=1, vjust=-1.8)

# use "ggbiplot" package
pca.results3<-prcomp(t(as.data.frame(assay(vst))), scale=TRUE)
group<-as.data.frame(colData(vst))
ggbiplot(pca.results3, obs.scale=1, var.scale=1, groups=group$Treatment, ellipse=TRUE, var.axes=F)

cat("==========================================================================================\n")
# Correlation & clustering
cat("Correlating samples and clustering...\n\n")
# use "GGally" package
#library(GGally)
cor.t24ct<-as.data.frame(assay(vst)[,1:5])
#pdf("t24ct_corr.pdf")
#ggpairs(cor.t24ct, title="t24ct correlation")
write.table(cor.t24ct, "t24ct.vst", quote=F, row.names=T, col.names=T, sep="\t")

cor.t24lps<-as.data.frame(assay(vst)[,6:10])
#pdf("t24lps_corr.pdf")
#ggpairs(cor.t24lps, title="t24lps correlation")
write.table(cor.t24lps, "t24lps.vst", quote=F, row.names=T, col.names=T, sep="\t")

cor.t28ct<-as.data.frame(assay(vst)[,11:15])
#pdf("t28ct_corr.pdf")
#ggpairs(cor.t28ct, title="t28ct correlation")
write.table(cor.t28ct, "t28ct.vst", quote=F, row.names=T, col.names=T, sep="\t")

cor.t28lps<-as.data.frame(assay(vst)[,16:20])
#pdf("t28lps_corr.pdf")
#ggpairs(cor.t28lps, title="t28lps correlation")
write.table(cor.t28lps, "t28lps.vst", quote=F, row.names=T, col.names=T, sep="\t")

cor.t32ct<-as.data.frame(assay(vst)[,21:25])
#pdf("t32ct_corr.pdf")
#ggpairs(cor.t32ct, title="t32ct correlation")
write.table(cor.t32ct, "t32ct.vst", quote=F, row.names=T, col.names=T, sep="\t")

cor.t32lps<-as.data.frame(assay(vst)[,26:30])
#pdf("t32lps_corr.pdf")
#ggpairs(cor.t32lps, title="t32lps correlation")
write.table(cor.t32lps, "t32lps.vst", quote=F, row.names=T, col.names=T, sep="\t")

# heatmap
corrCluster<-function(vst){
    sampleDist=dist(t(assay(vst)))
    sampleDist.mx=as.matrix(sampleDist)
    rownames(sampleDist.mx)=vst$Sample
    colnames(sampleDist.mx)=NULL
    colors=colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
    correlation.p<-pheatmap(as.data.frame(sampleDist.mx), clustering_distance_rows=sampleDist, clustering_distance_cols=sampleDist, col=colors)
}
corrCluster(vst)


cat("==========================================================================================\n")
# Remove outlier samples
cat("Removing outlier samples...\n\n")
if((length(rmsample)==1 && rmsample!="no") || length(rmsample)>1){for(i in rmsample){dds<-dds[,dds$Sample!=i]; sampleinfo<-sampleinfo[which(sampleinfo$Sample!=i),]}}
dds<-dds[rowSums(counts(dds))>=100,]
vst<-varianceStabilizingTransformation(dds, blind=FALSE)
if(comptype=="LPS"){
	write.table(as.data.frame(assay(vst)), "TransformmiRNA_byLPS.vst", quote=F, row.names=T, col.names=T, sep="\t")
}else if(comptype=="Temp"){
	write.table(as.data.frame(assay(vst)), "TransformmiRNA_byTemp.vst", quote=F, row.names=T, col.names=T, sep="\t")
}

# PCA & correlation again
cat("Ploting PCA again...\n\n")
# use DESeq2:plotPCA function
pca.result1<-PCA_deseq2(vst)
pca.plot1<-pca.result1[[1]]
data1<-pca.result1[[2]]
pca.plot1+coord_fixed(ratio=1)
pca.plot1+coord_fixed(ratio=1)+geom_text(aes(label=rownames(data1)), size=1, hjust=-0.4)
pca.plot1+coord_fixed(ratio=1)+geom_text(aes(label=rownames(data1)), size=1, vjust=-1.8)

# use in-house script
pca.result2<-PCA_inhouse(vst, sampleinfo)
pca.plot2<-pca.result2[[1]]
data2<-pca.result2[[2]]
pca.plot2
pca.plot2+geom_text(aes(label=rownames(data2)), size=1, hjust=-0.4)
pca.plot2+geom_text(aes(label=rownames(data2)), size=1, vjust=-1.8)

# use "ggbiplot" package
pca.results3<-prcomp(t(as.data.frame(assay(vst))), scale=TRUE)
group<-as.data.frame(colData(vst))
ggbiplot(pca.results3, obs.scale=1, var.scale=1, groups=group$Treatment, ellipse=TRUE, var.axes=F)

# correlation & clustering
corrCluster(vst)


cat("==========================================================================================\n")
# Isolate comparisons
cat("Isolating comparisons...\n\n")
remove_lowcount_miRNAs<-function(dds, nonzero, rowsum, cols1, cols2){
	dds.df<-as.data.frame(assay(dds))
	dds.df$keep<-rep("yes")
	for(i in rownames(dds.df)){
		if((sum(dds.df[i,cols1])>=sum(dds.df[i,cols2]) && sum(dds.df[i,cols1]>=1)>=nonzero && rowSums(dds.df[i,c(cols1, cols2)])>=rowsum) || (sum(dds.df[i,cols2])>=sum(dds.df[i,cols1]) && sum(dds.df[i,cols2]>=1)>=nonzero && rowSums(dds.df[i,c(cols1, cols2)])>=rowsum)){
			dds.df[i,"keep"]="yes"
		}else{
			dds.df[i,"keep"]="no"
		}
	}
	dds<-dds[rownames(dds.df[dds.df$keep=="yes",]),]
	return(dds)
}

nonzero<-as.numeric(args[5])
rowsum<-as.numeric(args[6])

if(comptype=="LPS"){
	dds24<-dds[,dds$Temperature==24]
	cols24.1=which(colData(dds24)$Treatment=="Control")
	cols24.2=which(colData(dds24)$Treatment=="LPS")
	dds24<-remove_lowcount_miRNAs(dds24, nonzero, rowsum, cols24.1, cols24.2)
	vst24<-varianceStabilizingTransformation(dds24, blind=FALSE)

	dds28<-dds[,dds$Temperature==28]
	cols28.1=which(colData(dds28)$Treatment=="Control")
	cols28.2=which(colData(dds28)$Treatment=="LPS")
	dds28<-remove_lowcount_miRNAs(dds28, nonzero, rowsum, cols28.1, cols28.2)
	vst28<-varianceStabilizingTransformation(dds28, blind=FALSE)

	dds32<-dds[,dds$Temperature==32]
	cols32.1=which(colData(dds32)$Treatment=="Control")
	cols32.2=which(colData(dds32)$Treatment=="LPS")
	dds32<-remove_lowcount_miRNAs(dds32, nonzero, rowsum, cols32.1, cols32.2)
	vst32<-varianceStabilizingTransformation(dds32, blind=FALSE)
}else if(comptype=="Temp"){
	dds24vs28<-dds[,dds$Temperature!=32 & dds$Treatment=="Control"]
	dds24vs28$Temperature<-droplevels(dds24vs28$Temperature)
	cols2428.1=which(colData(dds24vs28)$Temperature==24)
	cols2428.2=which(colData(dds24vs28)$Temperature==28)
	dds24vs28<-remove_lowcount_miRNAs(dds24vs28, nonzero, rowsum, cols2428.1, cols2428.2)
	vst24vs28<-varianceStabilizingTransformation(dds24vs28, blind=FALSE)

	dds32vs28<-dds[,dds$Temperature!=24 & dds$Treatment=="Control"]
	dds32vs28$Temperature<-droplevels(dds32vs28$Temperature)
	cols3228.1=which(colData(dds32vs28)$Temperature==28)
	cols3228.2=which(colData(dds32vs28)$Temperature==32)
	dds32vs28<-remove_lowcount_miRNAs(dds32vs28, nonzero, rowsum, cols3228.1, cols3228.2)
	vst32vs28<-varianceStabilizingTransformation(dds32vs28, blind=FALSE)
}


cat("==========================================================================================\n")
# Run DESeq2
cat("Running DESeq2...\n\n")
RunDESeq2<-function(dds, compString){
	cat("================================================================\n")
	cat("Running DESeq2 of", compString, "\n")
	dds<-DESeq(dds)

	normalized.counts<-as.data.frame(counts(dds, normalized=TRUE))
	normalized.counts$mean<-apply(normalized.counts, 1, mean)
	out.norm=paste("NormReadcount_", compString, ".tsv", sep="")
	write.table(normalized.counts, out.norm, row.names=T, col.names=T, quote=F, sep="\t")

	sf<-as.data.frame(sizeFactors(dds))
	out.sf<-paste("Sizefactor_", compString, ".tsv", sep="")
	write.table(sf, out.sf, row.names=T, col.names=T, quote=F, sep="\t")

	res=results(dds, alpha=0.05)
	summary(res)
	res1.5=res[which(res$padj<0.05 & abs(res$log2FoldChange)>=log2(1.5)),]
	res2=res[which(res$padj<0.05 & abs(res$log2FoldChange)>=log2(2)),]
	num.all=sum(res$padj<0.05 & abs(res$log2FoldChange)>=log2(1.5), na.rm=TRUE)
	num.up.fc1.5=sum(res$padj<0.05 & res$log2FoldChange>=log2(1.5), na.rm=TRUE)
	num.up.fc2=sum(res$padj<0.05 & res$log2FoldChange>=log2(2), na.rm=TRUE)
	num.down.fc1.5=sum(res$padj<0.05 & -(res$log2FoldChange)>=log2(1.5), na.rm=TRUE)
	num.down.fc2=sum(res$padj<0.05 & -(res$log2FoldChange)>=log2(2), na.rm=TRUE)
	cat("================================================================\n")
	cat("Comparison:", compString, "\n", "DEmiRNAs_all_|fc|>1.5:", num.all, "\n", "up_fc>1.5:", num.up.fc1.5, "\n", "up_fc>2:", num.up.fc2, "\n", "down_fc>1.5:", num.down.fc1.5, "\n", "down_fc>2:", num.down.fc2, "\n\n")

	out=list(dds, res, res1.5, res2)
	return(out)
}

if(comptype=="LPS"){
	t24lpsvsct<-RunDESeq2(dds24, "24lpsvsct")
	dds24<-t24lpsvsct[[1]]
	res24.all<-t24lpsvsct[[2]]
	res24.fc1.5<-t24lpsvsct[[3]]
	res24.fc2<-t24lpsvsct[[4]]
	write.table(as.data.frame(res24.fc1.5), "DEmiRNAs_24lpsvsct.csv", row.names=T, col.names=T, quote=F, sep="\t")

	t28lpsvsct<-RunDESeq2(dds28, "28lpsvsct")
	dds28<-t28lpsvsct[[1]]
	res28.all<-t28lpsvsct[[2]]
	res28.fc1.5<-t28lpsvsct[[3]]
	res28.fc2<-t28lpsvsct[[4]]
	write.table(as.data.frame(res28.fc1.5), "DEmiRNAs_28lpsvsct.csv", row.names=T, col.names=T, quote=F, sep="\t")

	t32lpsvsct<-RunDESeq2(dds32, "32lpsvsct")
	dds32<-t32lpsvsct[[1]]
	res32.all<-t32lpsvsct[[2]]
	res32.fc1.5<-t32lpsvsct[[3]]
	res32.fc2<-t32lpsvsct[[4]]
	write.table(as.data.frame(res32.fc1.5), "DEmiRNAs_32lpsvsct.csv", row.names=T, col.names=T, quote=F, sep="\t")
}else if(comptype=="Temp"){
	ct24vs28<-RunDESeq2(dds24vs28, "24vs28")
	dds24vs28<-ct24vs28[[1]]
	res24vs28.all<-ct24vs28[[2]]
	res24vs28.fc1.5<-ct24vs28[[3]]
	res24vs28.fc2<-ct24vs28[[4]]
	write.table(as.data.frame(res24vs28.fc1.5), "DEmiRNAs_24vs28.csv", row.names=T, col.names=T, quote=F, sep="\t")

	ct32vs28<-RunDESeq2(dds32vs28, "32vs28")
	dds32vs28<-ct32vs28[[1]]
	res32vs28.all<-ct32vs28[[2]]
	res32vs28.fc1.5<-ct32vs28[[3]]
	res32vs28.fc2<-ct32vs28[[4]]
	write.table(as.data.frame(res32vs28.fc1.5), "DEmiRNAs_32vs28.csv", row.names=T, col.names=T, quote=F, sep="\t")
}


cat("==========================================================================================\n")
# Volcano plot
cat("Volcano plot...\n\n")
VolcanoPlot<-function(res, compString, fc, fc.low, fc.high, padj.max, ratio){
	volcano=as.data.frame(res)
	volcano$significant=as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < 0.05 & abs(volcano$log2FoldChange) > log2(fc), ifelse(volcano$log2FoldChange>log2(fc), "Up", "Down"), "No"))
	volcano$padj=ifelse(is.na(volcano$padj), 1, volcano$padj)
	volcano$group<-rep(0)
	y.max=padj.max
	padj.max=10^-padj.max

	# group1: down-regulated DEmiRNAs within normal fc and padj ranges (blue dots)
	volcano[which(volcano$padj >= padj.max & volcano$padj < 0.05 & volcano$log2FoldChange >= fc.low & volcano$log2FoldChange <= -log2(fc)), 8]=rep(1)

	# group2: up-regulated DEmiRNAs within normal fc and padj ranges (red dots)
	volcano[which(volcano$padj >= padj.max & volcano$padj < 0.05 & volcano$log2FoldChange >= log2(fc) & volcano$log2FoldChange <= fc.high), 8]=rep(2)

	# group3: non-significant genes within normal fc or padj ranges (grey dots)
	volcano[which((volcano$padj >= 0.05 & abs(volcano$log2FoldChange) <= fc.high) | (volcano$padj < 0.05 & volcano$padj > padj.max & abs(volcano$log2FoldChange) <= log2(fc))), 8]=rep(3)

	# group4: non-significant genes outside of normal fc or padj ranges (grey trangle)
	volcano[which(volcano$padj >= 0.05 & volcano$log2FoldChange < fc.low), c("log2FoldChange", "group")]=list(log2FoldChange = fc.low, group=4)
	volcano[which(volcano$padj >= 0.05 & volcano$log2FoldChange > fc.high), c("log2FoldChange", "group")]=list(log2FoldChange = fc.high, group=4)
	volcano[which(volcano$padj < padj.max & abs(volcano$log2FoldChange) < -log2(fc)), c("padj", "group")]=list(padj = padj.max, group=4)

	# group5: down-regulated genes outside of normal fc or padj ranges (blue trangle)
	volcano[which(volcano$padj >= padj.max & volcano$padj < 0.05 &  volcano$log2FoldChange < fc.low), c("log2FoldChange", "group")]=list(log2FoldChange = fc.low, group=5)
	volcano[which(volcano$padj < padj.max & volcano$log2FoldChange < fc.low), c("log2FoldChange", "padj", "group")]=list(log2FoldChange = fc.low, padj = padj.max, group=5)
	volcano[which(volcano$padj < padj.max & volcano$log2FoldChange <= -log2(fc) & volcano$log2FoldChange >= fc.low), c("padj", "group")]=list(padj = padj.max, group=5)

	# group6: up-regulated genes outside of normal fc or padj ranges (red trangle)
	volcano[which(volcano$padj >= padj.max & volcano$padj < 0.05 & volcano$log2FoldChange > fc.high), c("log2FoldChange", "group")]=list(log2FoldChange = fc.high, group=6)
	volcano[which(volcano$padj < padj.max & volcano$log2FoldChange > fc.high), c("log2FoldChange", "padj", "group")]=list(log2FoldChange = fc.high, padj = padj.max, group=6)
	volcano[which(volcano$padj < padj.max & volcano$log2FoldChange >= log2(fc) & volcano$log2FoldChange <= fc.high), c("padj", "group")]=list(padj = padj.max, group=6)

	# plot
	title<-paste("Volcano plot of ", compString, " foldchange>", fc, sep="")
	volcano.p<-ggplot(volcano, aes(log2FoldChange, -log10(padj)))+geom_point(data=volcano[which(volcano$group==1),], color="royalblue3", alpha=0.8)+geom_point(data=volcano[which(volcano$group==2),], color=brewer.pal(11,"RdYlBu")[2], alpha=0.8)+geom_point(data=volcano[which(volcano$group==3),], color="gray60", alpha=0.8)+geom_point(data=volcano[which(volcano$group==4),], shape=2, color="gray60", alpha=0.8)+geom_point(data=volcano[which(volcano$group==5),], shape=2, color="royalblue3", alpha=0.8)+geom_point(data=volcano[which(volcano$group==6),], shape=2, color=brewer.pal(11,"RdYlBu")[2], alpha=0.8)
	volcano.p<-volcano.p+labs(title={title}, x="log2FoldChange", y="-log10(padj)")+geom_hline(yintercept=-log10(0.05), linetype=2, size=0.3, color="gray60")+geom_vline(xintercept=c(-log2(fc), log2(fc)), linetype=2, size=0.3, color="gray60")+xlim({fc.low}, {fc.high})+ylim(0, {y.max})+theme_bw()+theme(text=element_text(size=15), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_rect(size=1))+coord_fixed(ratio={ratio})
	print(volcano.p)
}

if(comptype=="LPS"){
	plotMA(res24.all, main="MA plot of 24lpsvsct", ylim=c(-10,10))
	plotMA(res28.all, main="MA plot of 28lpsvsct", ylim=c(-10,10))
	plotMA(res32.all, main="MA plot of 32lpsvsct", ylim=c(-10,10))

	VolcanoPlot(res24.all, "24lpsvsct", 1.5, -10, 10, 8, 2.5)
	VolcanoPlot(res24.all, "24lpsvsct", 2, -10, 10, 8, 2.5)
	VolcanoPlot(res28.all, "28lpsvsct", 1.5, -10, 10, 8, 2.5)
	VolcanoPlot(res28.all, "28lpsvsct", 2, -10, 10, 8, 2.5)
	VolcanoPlot(res32.all, "32lpsvsct", 1.5, -10, 10, 8, 2.5)
	VolcanoPlot(res32.all, "32lpsvsct", 2, -10, 10, 8, 2.5)
}else if(comptype=="Temp"){
	plotMA(res24vs28.all, main="MA plot of ct24vs28", ylim=c(-10,10))
	plotMA(res32vs28.all, main="MA plot of ct32vs28", ylim=c(-10,10))

	VolcanoPlot(res24vs28.all, "24vs28", 1.5, -10, 10, 8, 2.5)
	VolcanoPlot(res24vs28.all, "24vs28", 2, -10, 10, 8, 2.5)
	VolcanoPlot(res32vs28.all, "32vs28", 1.5, -10, 10, 8, 2.5)
	VolcanoPlot(res32vs28.all, "32vs28", 2, -10, 10, 8, 2.5)
}


cat("==========================================================================================\n")
# Heatmap of DEmiRNAs
cat("Heatmap of DEmiRNAs...\n\n")
if(comptype=="LPS"){
	anno.label<-c(rep("24 C control",5), rep("24 C LPS",5), rep("28 C control",5), rep("28 C LPS",5), rep("32 C control",5), rep("32 C LPS",5))
	names(anno.label)<-sampleinfo.origin$Sample
	anno.label<-anno.label[!names(anno.label) %in% rmsample]
	anno.label<-as.data.frame(anno.label)
	names(anno.label)<-"Group"
	anno.color<-list(Group=c("24 C control"="red", "24 C LPS"="black", "28 C control"="orange", "28 C LPS"="blue", "32 C control"="green", "32 C LPS"="yellow"))

	names.fc1.5<-unique(c(rownames(res24.fc1.5), rownames(res28.fc1.5), rownames(res32.fc1.5)))
	names.fc2<-unique(c(rownames(res24.fc2), rownames(res28.fc2), rownames(res32.fc2)))
	vst.demir1.5<-as.data.frame(assay(vst[names.fc1.5,]))
	vst.demir2<-as.data.frame(assay(vst[names.fc2,]))

	pheatmap(vst.demir1.5, main="Heatmap of LPSvsCT DEmiRNAs fc>1.5", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=F, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, cellwidth=10, cellheight=5)
	pheatmap(vst.demir1.5, main="Heatmap of LPSvsCT DEmiRNAs fc>1.5", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=TRUE, treeheight_col=20, show_rownames=F, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, cellwidth=10, cellheight=5)

	pheatmap(vst.demir2, main="Heatmap of LPSvsCT DEmiRNAs fc>2", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=F, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, cellwidth=10, cellheight=5)
	pheatmap(vst.demir2, main="Heatmap of LPSvsCT DEmiRNAs fc>2", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=TRUE, treeheight_col=20, show_rownames=F, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, cellwidth=10, cellheight=5)

}else if(comptype=="Temp"){
	anno.label<-c(rep("24 C control",5), rep("24 C LPS",5), rep("28 C control",5), rep("28 C LPS",5), rep("32 C control",5), rep("32 C LPS",5))
	names(anno.label)<-sampleinfo.origin$Sample
	anno.label<-anno.label[c(1:5, 11:15, 21:25)]
	anno.label<-anno.label[!names(anno.label) %in% rmsample]
	anno.label<-as.data.frame(anno.label)
	names(anno.label)<-"Group"
	anno.color<-list(Group=c("24 C control"="red", "24 C LPS"="black", "28 C control"="orange", "28 C LPS"="blue", "32 C control"="green", "32 C LPS"="yellow"))

	names.fc1.5<-unique(c(rownames(res24vs28.fc1.5), rownames(res32vs28.fc1.5)))
	names.fc2<-unique(c(rownames(res24vs28.fc2), rownames(res32vs28.fc2)))
	vst.demir1.5<-as.data.frame(assay(vst[names.fc1.5, vst$Treatment=="Control"]))
	vst.demir2<-as.data.frame(assay(vst[names.fc2, vst$Treatment=="Control"]))

	pheatmap(vst.demir1.5, main="Heatmap of 24/32_vs_28 DEmiRNAs fc>1.5", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=F, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, cellwidth=10, cellheight=5)
	pheatmap(vst.demir1.5, main="Heatmap of 24/32_vs_28 DEmiRNAs fc>1.5", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=TRUE, treeheight_col=20, show_rownames=F, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, cellwidth=10, cellheight=5)

	pheatmap(vst.demir2, main="Heatmap of 24/32_vs_28 DEmiRNAs fc>2", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=F, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, cellwidth=10, cellheight=5)
	pheatmap(vst.demir2, main="Heatmap of 24/32_vs_28 DEmiRNAs fc>2", scale="row", color=colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(3:9)]), brewer.pal(9, "OrRd")[c(3:9)]))(100), cluster_rows=TRUE, cluster_cols=TRUE, treeheight_col=20, show_rownames=F, annotation_col=anno.label, annotation_colors=anno.color, annotation_names_col=F, border_color=NA, cellwidth=10, cellheight=5)
}

#==================================================================================================
dev.off()
cat("Done with analysis!\n\n")

