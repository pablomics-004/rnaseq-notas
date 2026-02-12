library("recount3")

## Let's check all projects with human data
human_projects <- available_projects()

## Using a particular project
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)

## Getting the RangedSummarizedExperiment object
rse_gene_SRP009615 <- create_rse(proj_info)

## Explorng the time of creation of the object
rse_gene_SRP009615

## Exploring the object
rse_gene_SRP009615

## Exploring projects interactively
proj_info_interactive <- interactiveDisplayBase::display(human_projects)

## Let's transform the raw counts to counts per read
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)
rse_gene_SRP009615

## Experiment metadata
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[,
  # Fetches all columns with names starting with "sra_attribute"
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

## Replicating the figure
iSEE::iSEE(rse_gene_SRP009615)
