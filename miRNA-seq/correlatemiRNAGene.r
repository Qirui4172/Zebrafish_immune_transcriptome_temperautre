#!/usr/bin/env Rscript


#===========================================================================================================
# Set parameters
Usage<-function(){
	cat("Usage: Rscript correlatemiRNAGene.r [demirGroup] [mirnaGenes] [sampleInfo] [geneMatrix] [mirnaMatrix] [outfileName]\n\n",

	"Parameters:\n",
	"[demirGroup]	DEmiRNA comparison and group information (DEmiRNAs_group.txt), ONLY listing miRNAs having target genes (miRNAs having no target genes should be not included)!\n",
	"[mirnaGenes]	miRNA and target gene lists (mirnaTargets_overlap2.txt), ONLY listing miRNAs having target genes (miRNAs having no target genes should be not included)! ONLY the first two columns will be used even though it may have more columns, i.e. dre-miR-1-3p	ENSDARG00000073867	rex1bd\n",
	"[sampleInfo]	Sample information for correlation (sampleinfo_corr.txt)\n",
	"[geneMatrix]	Gene read counts matrix (TransformmiRNA_byLPS.vst), normalized and transformed by DEseq2\n",
	"[mirnaMatrix]	miRNA read counts matrix (TransformGene_byLPS.vst), normalized and transformed by DEseq2\n",
	"[outfileName]	Output file name\n\n",

	"Example:\n",
	"Rscript3.5.1 correlatemiRNAGene.r DEmiRNAs_group.txt mirnaTargets_overlap2.txt sampleinfo_corr.txt TransformGene_byLPS.vst TransformmiRNA_byLPS.vst Correlation_overlap2.tsv\n\n",

	"Function:\n",
	"Correlate miRNAs and target genes by computing Pearson correlation coefficients, p values, and FDRs (BH method).\n\n",

	"Contact: Qirui Zhang (qirui.zhang@med.lu.se)\n",
	"Date: 05-09-2018\n",
	"Updated: 05-01-2021\n\n"
	)
}

args<-commandArgs(TRUE)
if(length(args)!=6){Usage();quit();}


cat("==========================================================================================\n")
cat("Reading files ...\n\n")

# Read files
demir.group<-read.table(args[1], header=T, stringsAsFactors=F)
mirna.target<-read.table(args[2], header=F, stringsAsFactors=F, sep="\t")
mirna.target<-mirna.target[,1:2] # only use the first 2 columns
colnames(mirna.target)<-c("miRNA", "Ensembl")
sample.info<-read.table(args[3], header=T, stringsAsFactors=F)
gene.mx<-read.table(args[4], header=T, stringsAsFactors=F)
mirna.mx<-read.table(args[5], header=T, stringsAsFactors=F)


cat("==========================================================================================\n")
cat("Preparing correlation dataframe and matrix ...\n\n")

# Create miRNA-gene correlation dataframe
correlation<-data.frame(Group=character(), miRNA=character(), Ensembl=character(), stringsAsFactors=FALSE)

# Create coefficient-pvalue matrix
total_row.num<-0
for(i in 1:nrow(demir.group)){
	mir<-demir.group[i, "miRNA"]
	target.num<-nrow(mirna.target[which(mirna.target$miRNA==mir),])
	total_row.num<-total_row.num+target.num
}
cor.mx<-matrix(nrow=total_row.num, ncol=2)


cat("==========================================================================================\n")
cat("Performing correlation analysis ...\n\n")

# Correlate miRNAs and target genes
k<-0
for(i in 1:nrow(demir.group)){
	mir<-demir.group[i, "miRNA"]
	group<-demir.group[i, "Group"]
	temp1<-demir.group[i, "KeyTemperature"]
	keytreatment<-demir.group[i, "KeyTreatment"]

	# isolate samples for comparison
	if(keytreatment=="Control"){
		temp2<-as.integer(28)
		rmtemp<-ifelse(temp1==24, 32, 24)
		sample.names<-sample.info[which(sample.info$Treatment=="Control" & sample.info$Temperature!=rmtemp), "Sample"]
	}else if(keytreatment=="LPS"){
		sample.names<-sample.info[which(sample.info$Temperature==temp1), "Sample"]
	}

	# isolate subset of miRNA and gene matrix
	mirna_mx.sub<-mirna.mx[mir, colnames(mirna.mx) %in% sample.names]
	target.id<-mirna.target[which(mirna.target$miRNA==mir), "Ensembl"]
	gene_mx.sub<-gene.mx[target.id, colnames(gene.mx) %in% sample.names]
	mir_for_cor<-t(mirna_mx.sub)

	if(nrow(gene_mx.sub)==0){next;} # some miRNAs may have no target genes, and these miRNAs are supposed to be excluded beforehand but it still could be input by mistake (DEmiRNA_group, DEmiRNA_targetGene)!!

	# add to correlation dataframe
	mir.id<-rep(mir, length(target.id))
	group.id<-rep(group, length(target.id))
	df.tmp<-data.frame(Group=group.id, miRNA=mir.id, Ensembl=target.id, stringsAsFactors=F)
	correlation<-rbind(correlation, df.tmp)

	# compute correlation coefficient, pvalue and add to matrix
	for(j in 1:nrow(gene_mx.sub)){
		k<-k+1
		if(sum(is.na(gene_mx.sub[j,]))!=0 | sd(gene_mx.sub[j,])==0){coeff=NA; pvalue=NA; cor.mx[k,]=c(coeff, pvalue); next;} # check target genes and sd values, some predicted target genes may be not in gene.matrix file or some genes may have sd==0
		gene_for_cor<-t(gene_mx.sub[j,])
		cor.result<-cor.test(mir_for_cor, gene_for_cor)
		coeff<-cor.result$estimate[[1]]
		pvalue<-cor.result$p.value
		cor.mx[k,]<-c(coeff, pvalue)
	}
}

correlation<-cbind(correlation, as.data.frame(cor.mx))
colnames(correlation)[c(4,5)]<-c("Coefficient", "Pvalue")


cat("==========================================================================================\n")
cat("Adjusting p values using BH method ...\n\n")

correlation$BH_padj<-rep(0)
accumulated.num<-1
for(i in 1:nrow(demir.group)){
	mir<-demir.group[i, "miRNA"]
	group<-demir.group[i, "Group"]
	cor_sub<-correlation[which(correlation$Group==group & correlation$miRNA==mir),] # adjust pvalues of one miRNA at one time, not adjust pvalues of all miRNAs at the same time.
	sub_num<-nrow(cor_sub)
	start.num<-accumulated.num
	end.num<-accumulated.num+sub_num-1
	padj<-p.adjust(cor_sub$Pvalue, method="BH")
	correlation[start.num:end.num, "BH_padj"]<-padj
	accumulated.num<-accumulated.num+sub_num
}


cat("==========================================================================================\n")
outfile<-args[6]
cat("Outputing results to ", outfile, "\n\n")

write.table(correlation, outfile, col.names=T, row.names=F, quote=F, sep="\t")
cat("Done with correlation analysis!\n\n")


