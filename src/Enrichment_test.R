load("data/fisherTest16052022_140AFF.RData")

library(MODifieR)
library(tidyverse)
library(biomaRt)
library(DOSE)
library(clusterProfiler)

top_results <- fisherTest2 %>% filter(pvalue.fisher < 0.05)
top_genes <- top_results$gene_name

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get entrez id 
ids <- getBM(attributes = c(
  "entrezgene_id",
  "ensembl_gene_id",
  "external_gene_name"
  ),
  filter = "entrezgene_accession",
  values = top_genes,
  mart = mart)

ids$gene_name <- ids$external_gene_name
ids <- merge(ids, top_results, on="gene_name")
ids$pvalue <- ids$pvalue.fisher

gene_input <- create_custom_microarray_input_object(diff_genes = ids[, c("entrezgene_id", "pvalue")])


## PPI network
ppi_network <- read.delim("data/9606.protein.links.v11.5.txt",
                          stringsAsFactors = FALSE, sep = " ")
ppi_network <- ppi_network %>% filter(combined_score > 700)  # High confidence interaction
ppi_network[] <- lapply(ppi_network, as.character)

all_genes <- unique(c(ppi_network$protein1, ppi_network$protein2))
entrez <- getBM(attributes = c(
  "entrezgene_id",
  "ensembl_peptide_id"
  ),
  filter = "ensembl_peptide_id",
  values = all_genes,
  mart = mart)


entrez <- entrez[!duplicated(entrez$ensembl_peptide_id), ]
rownames(entrez) <- entrez[, "ensembl_peptide_id"]

ppi_network$entrezgene_id1 <- entrez[ppi_network$protein1, "entrezgene_id"]
ppi_network$entrezgene_id2 <- entrez[ppi_network$protein2, "entrezgene_id"]
ppi <- ppi_network[, c("entrezgene_id1", "entrezgene_id2", "combined_score")]
colnames(ppi) <- c("protein1", "protein2","combined_score")
ppi <- ppi[complete.cases(ppi),]


## KEGG enrichment
pathEnrichKEGG <- enrichKEGG(
  gene = ids$entrezgene_id,
  organism = "hsa", # Homo sapiens
  keyType = "kegg", # KEGG database
  pAdjustMethod = "BH")

pathEnrichKEGG@result %>% filter(qvalue < 0.05)

## DO enrichment
pathEnrichDO <- enrichDO(
  gene = ids$entrezgene_id,
  ont = "DO",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")

pathEnrichDO@result %>% filter(qvalue < 0.05)

