library("SummarizedExperiment") # Generating the SummarizedExperiment object

# Object for 200 genes in 6 samples
nrows <- 200
ncols <- 6

# Generate random counts for the assay from 1 to 10,000
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

## Gene information
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)), # Chromosomes for the genes
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100), # Start coordinates for each gene with 100 bp
  strand = sample(c("+", "-"), 200, TRUE), # Assign those genes randomly to each strand
  feature_id = sprintf("ID%03d", 1:200) # Unique identifier for each gene
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

## Sample information
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3), # Treatment groups for the samples
  row.names = LETTERS[1:6] # Sample labels
)

## FINAL OBJECT
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts), # Bioconductor's way to store assay data
  rowRanges = rowRanges,
  colData = colData
)

dim(rse) # Dimensions of the SummarizedExperiment object
# [1] 200   6
# GENES X SAMPLES

head(assay(rse)) # View the first few rows of the assay data

## Gene information through BioC
rowRanges(rse)

## Selecting two genes among all samples
rse[1:2, ]

## Selecting three samples across all genes
rse[, c("A", "D", "F")]

## Information about the samples
colData(rse)

## Visualizing the third gene with a boxplot
x <- assay(rse)[3, ][,] # Take all samples in Treatment
y <- y <- colData(rse)$Treatment
boxplot(
  x ~ y,
  main = "Gene 3",
  xlab = "Treatment",
  ylab = "Counts"
)
summary(lm(x ~ y)) # Statistical test to see if there is a significant difference between the two groups for gene 3
