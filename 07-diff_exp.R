## Fetching data from LIBD
library("recount3")

human_projects <- available_projects()

## Creating RSE object for SRP045638
rse_gene_SRP045638 <- create_rse(
  subset(
    human_projects,
    project == "SRP045638" & project_type == "data_sources"
  )
)

## Computing read counts
assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)

## There is an incoherence in the conditions (e.g. Fetal or dev-stage)
rse_gene_SRP045638$sra.sample_attributes[1:3]

## Correcting this problem
rse_gene_SRP045638$sra.sample_attributes <- gsub(
  "dev_stage;;Fetal\\|",
  "",
  rse_gene_SRP045638$sra.sample_attributes
)
rse_gene_SRP045638$sra.sample_attributes[1:3]

## Fetching metadata
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]

## Going from character to numeric or factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(
  rse_gene_SRP045638$sra_attribute.age
)
rse_gene_SRP045638$sra_attribute.disease <- factor(tolower(
  rse_gene_SRP045638$sra_attribute.disease
))
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(
  rse_gene_SRP045638$sra_attribute.RIN
)
rse_gene_SRP045638$sra_attribute.sex <- factor(
  rse_gene_SRP045638$sra_attribute.sex
)

## Variable summary
summary(as.data.frame(colData(rse_gene_SRP045638)[,
  grepl(
    "^sra_attribute.[age|disease|RIN|sex]",
    colnames(colData(rse_gene_SRP045638))
  )
]))

## Differences between prenatal and postnatal samples
## A problem might be that the coefficient prenatal wll be confused with the variable itself
rse_gene_SRP045638$prenatal <- factor(ifelse(
  rse_gene_SRP045638$sra_attribute.age < 0,
  "prenatal",
  "postnatal"
))
table(rse_gene_SRP045638$prenatal)

## http://rna.recount.bio/docs/quality-check-fields.html
## Control variable
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)

## Correlation between RIN and assigned gene proportion (exp vs high-throughput)
with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

## Hm... difference among groups?
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))

## Copy of the object before filtering
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Erase bad samples
hist(rse_gene_SRP045638$assigned_gene_prop)

## Filtering
table(rse_gene_SRP045638$assigned_gene_prop < 0.3)

rse_gene_SRP045638 <- rse_gene_SRP045638[,
  rse_gene_SRP045638$assigned_gene_prop > 0.3
]

## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los RPKMs o CPMs
## en vez de las cuentas.
## En realidad usariamos:
# edgeR::filterByExpr() https://bioconductor.org/packages/edgeR/ https://rdrr.io/bioc/edgeR/man/filterByExpr.html
# genefilter::genefilter() https://bioconductor.org/packages/genefilter/ https://rdrr.io/bioc/genefilter/man/genefilter.html
# jaffelab::expression_cutoff() http://research.libd.org/jaffelab/reference/expression_cutoff.html
#
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)

## Delete genes
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Final Dimensions
dim(rse_gene_SRP045638)

## Percentage of retained genes
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)

## ================================== DATA NORMALIZATION ==================================

library("edgeR") # BiocManager::install("edgeR", update = FALSE)}

dge <- DGEList(
  counts = assay(rse_gene_SRP045638, "counts"),
  genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge) ## Solving de contribution bias

## Brief exploratipon of the data
library("ggplot2")
ggplot(
  as.data.frame(colData(rse_gene_SRP045638)),
  aes(y = assigned_gene_prop, x = prenatal)
) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group")

## Checking out the coefficients of the model
## In first place we use the variable of most interest among all
mod <- model.matrix(
  ~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
  data = colData(rse_gene_SRP045638)
)
colnames(mod) # Female is the reference, that's why the name male goes with the variable sra_attribute.sexmale

## ================================== LIMMA ==================================
library("limma")
vGene <- voom(dge, mod, plot = TRUE)

## Estimating significant coefficients and their associated p-values
eb_results <- eBayes(lmFit(vGene))

## Ordering the values of the genes
de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP045638),
  sort.by = "none"
)
dim(de_results)

## EXploring the results
head(de_results)

## Differentially expressed genes with FDR < 0.05
table(de_results$adj.P.Val < 0.05)

## Visualization of the results
plotMA(eb_results, coef = 2)
