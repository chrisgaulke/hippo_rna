# Title and author information --------------------------------------------
#!/usr/bin/R

####################################
#                                  #
# 2022_02_16_hippocampus_rnaseq.R  #
#                                  #
####################################


#Copyright (C) 2022-2023  Christopher A. Gaulke
#author contact: cgaulke@illinois.edu
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#For a copy of the GNU General Public License see
#<http://www.gnu.org/licenses/>.


# Purpose -----------------------------------------------------------------

#Preliminary analysis of hippocampal RNA expression in viral infected
#piglets

#Note any section with the "SANDBOX" tag is exploratory

# ENVIRONMENT SETUP: Packages --------------------------------------------

library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(colorspace)
library(clusterProfiler)
library(grid)
library(gridExtra)
library(gtable)
library(biomaRt)
library("org.Ss.eg.db", character.only = TRUE)
library(ComplexHeatmap)

#BiocManager::install("biomaRt")
library(biomaRt)

options("stringsAsFactors"=F)
sessionInfo()

# ENVIRONMENT SETUP: Set up paths -------------------------------------------

indir  <- "../../data/read_counts_stranded/"

# DATA IMPORT: Import counts files ------------------------------------------

in_files <- list.files(indir)
samp_names <- NULL

for( i in 1:length(in_files)){
  samp_names <- c(samp_names,
                  paste0(
                    strsplit(in_files[i], "_")[[1]][1],
                         "_",
                    strsplit(in_files[i], "_")[[1]][2])
                  )
}


samp_group <- rep(c("unvaccinated_challenge", "vaccinated_challenged", "unvaccinated_not_challenge"), c(5,8,8))
samp_table <- data.frame(sampleName = samp_names,
                         fileName   = in_files,
                         condition  = samp_group)
#metadata for some analyses outside deseq
samp_table.meta <- samp_table
rownames(samp_table.meta) <- samp_table.meta$sampleName

#clean up files
hippoHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = samp_table,
                                        directory = indir,
                                        design = ~condition)
#filter low abundance stuff
hippoHTSeq <- hippoHTSeq[ which(rowSums(counts(hippoHTSeq)) > 1 ),,drop=F]

# DATA ANALYSIS: Build models ---------------------------------------------

#likelihood ratio test for GLMs this includes all samples
set.seed(731)
hippocamus_lrt <- DESeq(hippoHTSeq, test = "LRT", reduced = ~1)

result_lrt <- results(hippocamus_lrt)

result_1 <- results(hippocamus_lrt,
                    contrast = c("condition",
                                 "unvaccinated_challenge",
                                 "vaccinated_challenged"
                                 )
                    )


result_2 <- results(hippocamus_lrt,
                    contrast = c("condition",
                                 "unvaccinated_challenge",
                                 "unvaccinated_not_challenge"
                    )
)


result_3 <- results(hippocamus_lrt,
                    contrast = c("condition",
                                 "vaccinated_challenged",
                                 "unvaccinated_not_challenge"
                    )
)


#get number of sig
summary(result_1, alpha = .1)
summary(result_2, alpha = .1)
summary(result_3, alpha = .1)
summary(result_lrt, alpha = .1)

#plot fold change
#sig q < .1 blue

plotMA(result_1)
plotMA(result_2)
plotMA(result_3)


# DATA ANALYSIS: TRANSFORM DATA -------------------------------------------

#transform data
vsd <- vst(hippocamus_lrt, blind=FALSE)
rld <- rlog(hippocamus_lrt, blind=FALSE)

# DATA ANALYSIS: PCA PLOT -------------------------------------------------

# The DESeq2 plotPCA function is pretty feature poor, which is obvi not ideal.
# I checked the code and found they just use the prcomp from stats

# prcomp works fine with thousands of variables but the original code allows you
# to only look at the top n most variable genes so I will leave the code for
# that here, but commented.

#rv <- rowVars(assay(vsd))
#select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

#make a pca object
hippocamus.prcomp <- prcomp(t(assay(rld)))

#extract data and add metadata
pca.data <- data.frame(hippocamus.prcomp$x,
                       group = samp_table.meta[rownames(hippocamus.prcomp$x),"condition"])

#grab the % variance explained by each axis
pca.var  <- summary(hippocamus.prcomp)$importance[2,]

#plot PC1 and PC2
ggplot(data=pca.data, aes_string(x="PC1", y="PC2", color="group")) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(pca.var[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(pca.var[2] * 100),"% variance")) +
  coord_fixed()

#plot PC1 and PC3
ggplot(data=pca.data, aes_string(x="PC1", y="PC3", color="group")) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(pca.var[1] * 100),"% variance")) +
  ylab(paste0("PC3: ",round(pca.var[3] * 100),"% variance")) +
  coord_fixed()

# DATA ANALYSIS: HEAT MAP-----------------------------------------------------
#make a simple heat map using log transformed data

rld.heatmap <- rld[which(rownames(rld) %in% rownames(result_lrt[which(result_lrt$padj < 0.05 ),])),]
rld.heatmap <- assay(rld.heatmap)
rld.heatmap <- t(scale(t(rld.heatmap))) #scale by rows because ComplexHeatmap wont do it for us

#make top annotations
col_annotation <- HeatmapAnnotation(Group = samp_table.meta$condition,
                                    col = list(Group =c("unvaccinated_challenge" = "#0B3954",
                                               "vaccinated_challenged"  = "#087E8B",
                                               "unvaccinated_not_challenge" = "#BFD7EA") ))


#make the heatmap
rna.heat <- Heatmap(rld.heatmap,
        top_annotation = col_annotation,
        show_row_names =  F,
        show_column_names = T, row_split = 5)

rna.heat6 <- Heatmap(rld.heatmap,
                    top_annotation = col_annotation,
                    show_row_names =  F,
                    show_column_names = T, row_split = 6)

rna.heat7 <- Heatmap(rld.heatmap,
                     top_annotation = col_annotation,
                     show_row_names =  F,
                     show_column_names = T, row_split = 7)

vir_up <- rownames(rld.heatmap[row_order(rna.heat)[[1]],])
vir_down <- rownames(rld.heatmap[row_order(rna.heat)[[3]],])
control_up <- rownames(rld.heatmap[row_order(rna.heat)[[5]],])
vax_up <- rownames(rld.heatmap[row_order(rna.heat)[[2]],])


# DATA ANALYSIS: GSEA LRT-----------------------------------------------------

sig_genes_lrt.df <- result_lrt[which(result_lrt$padj < 0.1 ),]

#build a named vector
sig_genes_lrt.genes <- sig_genes_lrt.df$log2FoldChange
names(sig_genes_lrt.genes) <- rownames(sig_genes_lrt.df)

#remove NA
sig_genes_lrt.genes <- na.omit(sig_genes_lrt.genes)

#sort
sig_genes_lrt.genes <- sort(sig_genes_lrt.genes, decreasing = T)

# We will need to convert the Ensembl ids to something else because org.Ss.eg.db
# does not contain a mapping

ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
ensemble_ids <- names(sig_genes_lrt.genes)

ngene_list <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'),
                    filters = c("ensembl_gene_id", "with_entrezgene"),
                    values = list(ensemble_ids, TRUE),
                    mart = ensembl)

# There are some duplicate entries that need to be dealt with
# We will use the most conservative approach which will remove any genes that
# map to multiple entrez symbols

dups <- ngene_list$ensembl_gene_id[duplicated.default(ngene_list$ensembl_gene_id)]

ngene_list.dedup <- ngene_list[-which(ngene_list$ensembl_gene_id %in% dups), ]

#now we filter
sig_genes_lrt.genes <- sig_genes_lrt.genes[which(names(sig_genes_lrt.genes) %in% ngene_list.dedup$ensembl_gene_id)]

#order so that names align
sig_genes_lrt.genes <- sig_genes_lrt.genes[ngene_list.dedup$ensembl_gene_id]

#now swap names

names(sig_genes_lrt.genes) <- ngene_list.dedup$entrezgene_id


#now sort to ensure vector is in decreasing order
sig_genes_lrt.genes <- sort(sig_genes_lrt.genes, decreasing = T)

lrt.gse <- gseGO(geneList = sig_genes_lrt.genes,
                 ont = "BP",
                 keyType = "ENTREZID",
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 1,
                 verbose = TRUE,
                 OrgDb = org.Ss.eg.db,
                 pAdjustMethod = "fdr"
                 )

lrt.ego <- enrichGO(gene          = names(sig_genes_lrt.genes[abs(sig_genes_lrt.genes) > .5]),
                    universe      = names(sig_genes_lrt.genes),
                    OrgDb         = org.Ss.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE)

lrt.eko <- enrichKEGG(gene          = names(sig_genes_lrt.genes[abs(sig_genes_lrt.genes) > .5]),
                    pvalueCutoff   = 1,
                    organism       = 'ssc'

                    )

lrt.kgse <- gseKEGG(geneList = sig_genes_lrt.genes[abs(sig_genes_lrt.genes) > .5],
                      pvalueCutoff   = 1,
                      organism       = 'ssc',
                      minGSSize      = 5

)


lrt.meko <- enrichMKEGG(gene         = names(sig_genes_lrt.genes[abs(sig_genes_lrt.genes) > .5]),
                      pvalueCutoff   = 1,
                      qvalueCutoff = 1,
                      organism       = 'ssc'

)

lrt.mkgse <- gseMKEGG(geneList = sig_genes_lrt.genes[abs(sig_genes_lrt.genes) > .5],
                    pvalueCutoff   = 1,
                    organism       = 'ssc',
                    minGSSize      = 5

)

# DATA ANALYSIS: GSEA unvax_challenged vs vax_challenged-----------------------------------------------------
#get sig genes. Be careful, we will reuse some variables here

sig_genes.df <- result_1[which(result_1$padj < 0.1 ),]

#build a named vector
sig_genes.genes <- sig_genes.df$log2FoldChange
names(sig_genes.genes) <- rownames(sig_genes.df)

#remove NA
sig_genes.genes <- na.omit(sig_genes.genes)

#sort
sig_genes.genes <- sort(sig_genes.genes, decreasing = T)


# We will need to convert the Ensembl ids to something else because org.Ss.eg.db
# does not contain a mapping

ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
ensemble_ids <- names(sig_genes.genes)

ngene_list <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'),
                    filters = c("ensembl_gene_id", "with_entrezgene"),
                    values = list(ensemble_ids, TRUE),
                    mart = ensembl)

# There are some duplicate entries that need to be dealt with
# We will use the most conservative approach which will remove any genes that
# map to multiple entrez symbols

dups <- ngene_list$ensembl_gene_id[duplicated.default(ngene_list$ensembl_gene_id)]

ngene_list.dedup <- ngene_list[-which(ngene_list$ensembl_gene_id %in% dups), ]

#now we filter
sig_genes.genes <- sig_genes.genes[which(names(sig_genes.genes) %in% ngene_list.dedup$ensembl_gene_id)]

#order so that names align
sig_genes.genes <- sig_genes.genes[ngene_list.dedup$ensembl_gene_id]

#now swap names

names(sig_genes.genes) <- ngene_list.dedup$entrezgene_id


#now sort to ensure vector is in decreasing order
sig_genes.genes <- sort(sig_genes.genes, decreasing = T)


result_1.ego <- enrichGO(gene          = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                         universe      = names(sig_genes.genes),
                         OrgDb         = org.Ss.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "fdr",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         readable      = TRUE)

result_1.gse <- gseGO(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                      ont = "BP",
                      keyType = "ENTREZID",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 1,
                      verbose = TRUE,
                      OrgDb = org.Ss.eg.db,
                      pAdjustMethod = "fdr"
)


#KEGG KO enrichment
result_1.eko <- enrichKEGG(gene          = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                           pvalueCutoff   = 1,
                           organism       = 'ssc'

)
#KEGG KO GSEA
result_1.kgse <- gseKEGG(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                         pvalueCutoff   = 1,
                         organism       = 'ssc',
                         minGSSize      = 5

)

#KEGG module enrichment
result_1.meko <- enrichMKEGG(gene         = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                             pvalueCutoff   = 1,
                             qvalueCutoff = 1,
                             organism       = 'ssc'

)

#KEGG module GSEA
result_1.mkgse <- gseMKEGG(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                           pvalueCutoff   = 1,
                           organism       = 'ssc',
                           minGSSize      = 5

)


# DATA ANALYSIS: GSEA unvax_challenged vs unvax_not_challenged-----------------------------------------------------
#get sig genes


sig_genes.df <- result_2[which(result_2$padj < 0.1 ),]

#build a named vector
sig_genes.genes <- sig_genes.df$log2FoldChange
names(sig_genes.genes) <- rownames(sig_genes.df)

#remove NA
sig_genes.genes <- na.omit(sig_genes.genes)

#sort
sig_genes.genes <- sort(sig_genes.genes, decreasing = T)



# We will need to convert the Ensembl ids to something else because org.Ss.eg.db
# does not contain a mapping

ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
ensemble_ids <- names(sig_genes.genes)

ngene_list <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'),
                    filters = c("ensembl_gene_id", "with_entrezgene"),
                    values = list(ensemble_ids, TRUE),
                    mart = ensembl)

# There are some duplicate entries that need to be dealt with
# We will use the most conservative approach which will remove any genes that
# map to multiple entrez symbols

dups <- ngene_list$ensembl_gene_id[duplicated.default(ngene_list$ensembl_gene_id)]

ngene_list.dedup <- ngene_list[-which(ngene_list$ensembl_gene_id %in% dups), ]

#now we filter
sig_genes.genes <- sig_genes.genes[which(names(sig_genes.genes) %in% ngene_list.dedup$ensembl_gene_id)]

#order so that names align
sig_genes.genes <- sig_genes.genes[ngene_list.dedup$ensembl_gene_id]

#now swap names

names(sig_genes.genes) <- ngene_list.dedup$entrezgene_id


#now sort to ensure vector is in decreasing order
sig_genes.genes <- sort(sig_genes.genes, decreasing = T)


result_2.ego <- enrichGO(gene          = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                         universe      = names(sig_genes.genes),
                         OrgDb         = org.Ss.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "fdr",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         readable      = TRUE)

result_2.gse <- gseGO(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                      ont = "BP",
                      keyType = "ENTREZID",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 1,
                      verbose = TRUE,
                      OrgDb = org.Ss.eg.db,
                      pAdjustMethod = "fdr"
)


#KEGG KO enrichment
result_2.eko <- enrichKEGG(gene          = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                           pvalueCutoff   = 1,
                           organism       = 'ssc'

)
#KEGG KO GSEA
result_2.kgse <- gseKEGG(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                         pvalueCutoff   = 1,
                         organism       = 'ssc',
                         minGSSize      = 5

)

#KEGG module enrichment
result_2.meko <- enrichMKEGG(gene         = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                             pvalueCutoff   = 1,
                             qvalueCutoff = 1,
                             organism       = 'ssc'

)

#KEGG module GSEA
result_2.mkgse <- gseMKEGG(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                           pvalueCutoff   = 1,
                           organism       = 'ssc',
                           minGSSize      = 5

)


# DATA ANALYSIS: GSEA vax_challenged vs unvax_not_challenged-----------------------------------------------------
#get sig genes


sig_genes.df <- result_3[which(result_3$padj < 0.1 ),]

#build a named vector
sig_genes.genes <- sig_genes.df$log2FoldChange
names(sig_genes.genes) <- rownames(sig_genes.df)

#remove NA
sig_genes.genes <- na.omit(sig_genes.genes)

#sort
sig_genes.genes <- sort(sig_genes.genes, decreasing = T)




# We will need to convert the Ensembl ids to something else because org.Ss.eg.db
# does not contain a mapping

ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
ensemble_ids <- names(sig_genes.genes)

ngene_list <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'),
                    filters = c("ensembl_gene_id", "with_entrezgene"),
                    values = list(ensemble_ids, TRUE),
                    mart = ensembl)

# There are some duplicate entries that need to be dealt with
# We will use the most conservative approach which will remove any genes that
# map to multiple entrez symbols

dups <- ngene_list$ensembl_gene_id[duplicated.default(ngene_list$ensembl_gene_id)]

ngene_list.dedup <- ngene_list[-which(ngene_list$ensembl_gene_id %in% dups), ]

#now we filter
sig_genes.genes <- sig_genes.genes[which(names(sig_genes.genes) %in% ngene_list.dedup$ensembl_gene_id)]

#order so that names align
sig_genes.genes <- sig_genes.genes[ngene_list.dedup$ensembl_gene_id]

#now swap names

names(sig_genes.genes) <- ngene_list.dedup$entrezgene_id


#now sort to ensure vector is in decreasing order
sig_genes.genes <- sort(sig_genes.genes, decreasing = T)


result_3.ego <- enrichGO(gene          = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                         universe      = names(sig_genes.genes),
                         OrgDb         = org.Ss.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "fdr",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         readable      = TRUE)

result_3.gse <- gseGO(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                      ont = "BP",
                      keyType = "ENTREZID",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 1,
                      verbose = TRUE,
                      OrgDb = org.Ss.eg.db,
                      pAdjustMethod = "fdr"
)


#KEGG KO enrichment
result_3.eko <- enrichKEGG(gene          = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                           pvalueCutoff   = 1,
                           organism       = 'ssc'

)
#KEGG KO GSEA
result_3.kgse <- gseKEGG(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                         pvalueCutoff   = 1,
                         organism       = 'ssc',
                         minGSSize      = 5

)

#KEGG module enrichment
result_3.meko <- enrichMKEGG(gene         = names(sig_genes.genes[abs(sig_genes.genes) > .5]),
                             pvalueCutoff   = 1,
                             qvalueCutoff = 1,
                             organism       = 'ssc'

)

#KEGG module GSEA
result_3.mkgse <- gseMKEGG(geneList = sig_genes.genes[abs(sig_genes.genes) > .5],
                           pvalueCutoff   = 1,
                           organism       = 'ssc',
                           minGSSize      = 5

)

# SANDBOX: 3D plot  ----------------------------------------------------------------
#
#Run only if you don't mind shutting down when the connection gets broken... Have you saved lately?
#
# #library(rgl)
# mycolors <- c('royalblue1', 'darkcyan', 'oldlace')
# pca.data$color <- mycolors[ as.numeric(as.factor(pca.data$group)) ]
#
# plot3d(
#   x=pca.data$PC1, y=pca.data$PC2, z=pca.data$PC3,
#   col = pca.data$color,
#   type = 's',
#   radius = 1,
#   xlab="PC1", ylab="PC2", zlab="PC3")


# SANDBOX 2: MORE SAND !!!  -----------------------------------------------

#upsig_genes.df <- result_1[which(result_1$padj < 0.1 ),]
upsig_genes.df <- result_2[which(rownames(result_2) %in% control_up),]


#build a named vector
upsig_genes.genes <- upsig_genes.df$log2FoldChange
names(upsig_genes.genes) <- rownames(upsig_genes.df)

#remove NA
upsig_genes.genes <- na.omit(upsig_genes.genes)

#sort
upsig_genes.genes <- sort(upsig_genes.genes, decreasing = T)


# We will need to convert the Ensembl ids to something else because org.Ss.eg.db
# does not contain a mapping

ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
upensemble_ids <- names(upsig_genes.genes)

upngene_list <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'),
                    filters = c("ensembl_gene_id", "with_entrezgene"),
                    values = list(upensemble_ids, TRUE),
                    mart = ensembl)

# There are some duplicate entries that need to be dealt with
# We will use the most conservative approach which will remove any genes that
# map to multiple entrez symbols

dups <- upngene_list$ensembl_gene_id[duplicated.default(upngene_list$ensembl_gene_id)]

upngene_list.dedup <- upngene_list[-which(upngene_list$ensembl_gene_id %in% dups), ]

#now we filter
upsig_genes.genes <- upsig_genes.genes[which(names(upsig_genes.genes) %in% upngene_list.dedup$ensembl_gene_id)]

#order so that names align
upsig_genes.genes <- upsig_genes.genes[upngene_list.dedup$ensembl_gene_id]

#now swap names

names(upsig_genes.genes) <- upngene_list.dedup$entrezgene_id


#now sort to ensure vector is in decreasing order
upsig_genes.genes <- sort(upsig_genes.genes, decreasing = T)


upresult_1.ego <- enrichGO(gene        =  names(upsig_genes.genes),#names(upsig_genes.genes[upsig_genes.genes > .5]),
                         #universe      = names(upsig_genes.genes),
                         OrgDb         = org.Ss.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "fdr",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         readable      = TRUE)


#down

downsig_genes.df <- result_1[which(rownames(result_1) %in% vir_down),]


#build a named vector
downsig_genes.genes <- downsig_genes.df$log2FoldChange
names(downsig_genes.genes) <- rownames(downsig_genes.df)

#remove NA
downsig_genes.genes <- na.omit(downsig_genes.genes)

#sort
downsig_genes.genes <- sort(downsig_genes.genes, decreasing = T)


# We will need to convert the Ensembl ids to something else because org.Ss.eg.db
# does not contain a mapping

#ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
downensemble_ids <- names(downsig_genes.genes)

downngene_list <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'),
                        filters = c("ensembl_gene_id", "with_entrezgene"),
                        values = list(downensemble_ids, TRUE),
                        mart = ensembl)

# There are some duplicate entries that need to be dealt with
# We will use the most conservative approach which will remove any genes that
# map to multiple entrez symbols

dups <- downngene_list$ensembl_gene_id[duplicated.default(downngene_list$ensembl_gene_id)]

downngene_list.dedup <- downngene_list[-which(downngene_list$ensembl_gene_id %in% dups), ]

#now we filter
downsig_genes.genes <- downsig_genes.genes[which(names(downsig_genes.genes) %in% downngene_list.dedup$ensembl_gene_id)]

#order so that names align
downsig_genes.genes <- downsig_genes.genes[downngene_list.dedup$ensembl_gene_id]

#now swap names

names(downsig_genes.genes) <- downngene_list.dedup$entrezgene_id


#now sort to ensure vector is in decreasing order
downsig_genes.genes <- sort(downsig_genes.genes, decreasing = T)


downresult_1.ego <- enrichGO(gene         = names(downsig_genes.genes),
                             #universe      = names(downsig_genes.genes),
                             OrgDb         = org.Ss.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "fdr",
                             pvalueCutoff  = 1,
                             qvalueCutoff  = 1,
                             readable      = TRUE)


# SANDBOX : SANDCASTLE ----------------------------------------------------

#make biomart object
ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")

ensembl2entrez <- function(tdf, genes, bm){
  #tdf is a results object of DESeq2
  #genes is a list of genes of iterest
  #bm is the biomart object containing gene symbols
  #returns a sorted named vector(gene ids) of with FC values

  goi.df <-  tdf[which(rownames(tdf) %in% genes),]

  #build a named vector
  goi.genes <- goi.df$log2FoldChange
  names(goi.genes) <- rownames(goi.df)
  goi.genes <- na.omit(goi.genes)

  #sort
  goi.genes <- sort(goi.genes, decreasing = T)


  # convert the Ensembl ids to something else because org.Ss.eg.db
  # does not contain a mapping

  goi.ens <- names(goi.genes)

  goi_list <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'),
                        filters = c("ensembl_gene_id", "with_entrezgene"),
                        values = list(goi.ens, TRUE),
                        mart = bm)

  # There are some duplicate entries that need to be dealt with
  # We will use the most conservative approach which will remove any genes that
  # map to multiple entrez symbols

  dups <- goi_list$ensembl_gene_id[duplicated.default(goi_list$ensembl_gene_id)]

  goi.dedup <- goi_list[-which(goi_list$ensembl_gene_id %in% dups), ]

  #now we filter
  goi.genes <- goi.genes[which(names(goi.genes) %in% goi.dedup$ensembl_gene_id)]

  #order so that names align
  goi.genes <- goi.genes[goi.dedup$ensembl_gene_id]

  #now swap names

  names(goi.genes) <- goi.dedup$entrezgene_id


  #now sort to ensure vector is in decreasing order
  goi.genes <- sort(goi.genes, decreasing = T)


  return(goi.genes)

}

vir_down.list <- ensembl2entrez(tdf = result_1, genes = vir_down,bm = ensembl)
vir_up.list <- ensembl2entrez(tdf = result_1, genes = vir_up,bm = ensembl)
control_up.list <- ensembl2entrez(tdf = result_2, genes = control_up,bm = ensembl)


vir_down.ego <- enrichGO(gene        =  names(vir_down.list),
                           #universe      = names(upsig_genes.genes),
                           OrgDb         = org.Ss.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "fdr",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1,
                           readable      = TRUE)


vir_up.ego <- enrichGO(gene        =  names(vir_up.list),
                         #universe      = names(upsig_genes.genes),
                         OrgDb         = org.Ss.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "fdr",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         readable      = TRUE)

control_up.ego <- enrichGO(gene        =  names(control_up.list),
                       #universe      = names(upsig_genes.genes),
                       OrgDb         = org.Ss.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1,
                       readable      = TRUE)



#y mas
vir_up.gse <- gseGO(geneList = vir_up.list,
                      ont = "BP",
                      keyType = "ENTREZID",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 1,
                      verbose = TRUE,
                      OrgDb = org.Ss.eg.db,
                      pAdjustMethod = "fdr",
                    scoreType = "pos"
)


#KEGG KO enrichment
vir_up.eko <- enrichKEGG(gene          = names(vir_up.list[abs(vir_up.list) > .5]),
                           pvalueCutoff   = 1,
                           organism       = 'ssc'

)

#KEGG KO GSEA
vir_up.kgse <- gseKEGG(geneList = vir_up.list[abs(vir_up.list) > .5],
                         pvalueCutoff   = 1,
                         organism       = 'ssc',
                         minGSSize      = 5,
                       scoreType = "pos"

)

#KEGG module enrichment
vir_up.meko <- enrichMKEGG(gene         = names(vir_up.list),
                             pvalueCutoff   = 1,
                             qvalueCutoff = 1,
                             organism       = 'ssc'

)

#down
vir_down.gse <- gseGO(geneList = vir_down.list,
                      ont = "BP",
                      keyType = "ENTREZID",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 1,
                      verbose = TRUE,
                      OrgDb = org.Ss.eg.db,
                      pAdjustMethod = "fdr",
                      scoreType = "neg"
)


#KEGG KO enrichment
vir_down.eko <- enrichKEGG(gene          = names(vir_down.list[abs(vir_down.list) > .5]),
                           pvalueCutoff   = 1,
                           organism       = 'ssc'

)

#KEGG KO GSEA
vir_down.kgse <- gseKEGG(geneList = vir_down.list[abs(vir_down.list) > .0],
                         pvalueCutoff   = 1,
                         organism       = 'ssc',
                         minGSSize      = 5,
                         scoreType = "neg"

)

#KEGG module enrichment
vir_down.meko <- enrichMKEGG(gene         = names(vir_down.list),
                             pvalueCutoff   = 1,
                             qvalueCutoff = 1,
                             organism       = 'ssc'

)


#control


control_up.gse <- gseGO(geneList = control_up.list,
                        ont = "BP",
                        keyType = "ENTREZID",
                        minGSSize = 3,
                        maxGSSize = 800,
                        pvalueCutoff = 1,
                        verbose = TRUE,
                        OrgDb = org.Ss.eg.db,
                        pAdjustMethod = "fdr",
                        scoreType = "neg"
)


#KEGG KO enrichment
control_up.eko <- enrichKEGG(gene          = names(control_up.list[abs(control_up.list) > .5]),
                             pvalueCutoff   = 1,
                             organism       = 'ssc',
                             qvalueCutoff = 1

)

#KEGG KO GSEA
control_up.kgse <- gseKEGG(geneList = control_up.list[abs(control_up.list) > .5],
                           pvalueCutoff   = 1,
                           organism       = 'ssc',
                           minGSSize      = 5,
                           scoreType = "neg"

)

#KEGG module enrichment
control_up.meko <- enrichMKEGG(gene         = names(control_up.list),
                               pvalueCutoff   = 1,
                               qvalueCutoff = 1,
                               organism       = 'ssc'

)




# SANDBOX 3: Mapping ensembl to entrez/hngc---------------------------------------------------------------


ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")

temp.df <-  result_1[which(rownames(result_1) %in% vax_up),]
temp.rnames <- rownames(temp.df)

# convert the Ensembl ids to entrez and hugo

temp_list <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'),
                  filters = c("ensembl_gene_id", "with_entrezgene"),
                  values = list(temp.rnames, TRUE),
                  mart = ensembl)


# SANDBOX 3: make boxplot---------------------------------------------------------------

pdf("~/Desktop/boxplot.pdf")
boxplot(as.vector(assay(rld[rownames(result_1)[1],])) ~ samp_table.meta$condition, col = "steelblue")
dev.off()
#ggplot ready data
data.frame(value = as.vector(assay(rld[rownames(result_1)[1],])),
           group = samp_table.meta$condition )


