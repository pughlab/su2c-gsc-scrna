Bulk RNA-seq analyses, CRISPR GSEA analysis, plus additional scRNA analysis here, including preprocessing/exploration of Wang 2019, Bhaduri 2020, and Neftel 2019 datasets, and velocyto/scvelo analyses.
Note that bulk RNA-seq analyses and any GSEA related analyses require the the SummExpDR package (version 0.1.040) to run (https://github.com/BaderLab/SummExpDR/tree/dev/R).

all scripts except those under velocyto_cluster_scripts were run on a workstation with Ubuntu 16.04, 64GB RAM. The analyses under velocyto_cluster_scripts were run on a CentOS cluster with a slurm workload manager.

Currently, preprocessed data is ommited due to file size constraints for github, but analyses should be repeatable given publicly available data associated with publication
