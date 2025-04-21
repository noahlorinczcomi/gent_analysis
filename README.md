This repository contains the following data:
- [data_analysis_scripts](https://github.com/noahlorinczcomi/gent_analysis/tree/main/data_analysis_scripts)
  - R scripts used to analyze data in [preprint]([https://github.com/noahlorinczcomi/gent_analysis](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5080346))
- [gene_test_results](https://github.com/noahlorinczcomi/gent_analysis/tree/main/gene_test_results)
  - Gene-level results from gene-based association tested with GenT applied to multiple phenotypes
- [finemaps](https://github.com/noahlorinczcomi/gent_analysis/tree/main/finemaps)
  - Gene-SNP pair results from finemapping using SuSiE applied to each gene and multiple phenotypes
- [simulations](https://github.com/noahlorinczcomi/gent_analysis/tree/main/simulations)
  - R code used to perform simulations used to evaluate GenT, MuGenT, MuGenT-PH, MuGenT-Pleio, MuGenT-Sel, xGenT, finemapping of gene-based test statistics, and comparisons of Type I/II error between GenT, [ACAT](https://github.com/yaowuliu/ACAT), [VEGAS](https://github.com/HimesGroup/snpsettest), [GATES](https://pmc.ncbi.nlm.nih.gov/articles/PMC3059433/), and [exset](https://github.com/noahlorinczcomi/exset).
- [druggable_genes.txt](https://raw.githubusercontent.com/noahlorinczcomi/gent_analysis/main/druggable_genes.txt)
  - List of 3,370 genes with drug target interactions ('druggable genes')
- [druggable_genes_high_affinity.txt](https://raw.githubusercontent.com/noahlorinczcomi/gent_analysis/main/druggable_genes_high_affinity.txt)
  - List of 1,230 genes with drug target interactions with high binding affinity

We made a Shiny web application where you can view and download gene-based association testing (GenT) and fine-mapping results for over 50 traits here: [https://nlorinczcomi.shinyapps.io/gent](https://nlorinczcomi.shinyapps.io/gent).
