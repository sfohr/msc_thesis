# Nonlinear methods for dimensionality reduction and clustering of bacterial single-cell  sequencing data

Repo containing code for my master's thesis in biostatistics at LMU Munich. 
I took a look on how to use Nonlinear Matrix Decomposition (NMD) (Saul, L., 2022) in the context of bacterial scRNA-seq analysis (Heumos, L., et. al. 2023), replacing Principal Component Analysis in the optimized workflow, as outlined in Ostner, J. (2024).

My thesis was structured along the following objectives:

- implement the algorithms from Seraghiti, G., et. al. (2023) in the Python module [nomad](https://github.com/flatironinstitute/nomad/) in cooperation with [Flatiron Institute](https://www.simonsfoundation.org/flatiron/)
- code for the simulation study of the algorithms in Seraghiti, G., et. al. (2023) with varying sparsity can be found in `/simulation`
- apply NMD in the context of the BacSC workflow (Ostner, J., et. al. (2024)) on raw and normalized counts (found in `/application/analysis`), also for manually set number of latent dimensions
- explore NMD's potential for imputation of [sampling zeros](https://www.nature.com/articles/s41467-021-27729-z) (check `/application/NMD_zero_imputation /`)
- potential of Poisson-Hurdle model-based clustering (Qiao, Z., et. al. (2023)) for scRNA-seq (`/application/poisson_hurdle`).

Raw and intermediate data, as well as figures, can be downloaded from [zenodo](https://zenodo.org/records/13898901).

## References
Saul, L., K., (2022), "A Nonlinear Matrix Decomposition for Mining the Zeros of Sparse Data"
[https://doi.org/10.1137/21M1405769](https://doi.org/10.1137/21M1405769)
(Preprint: https://cseweb.ucsd.edu/~saul/papers/preprints/simods22_preprint.pdf).

Saul, L., K., (2022), "A geometrical connection between sparse and low-rank
matrices and its application to manifold learning" [https://openreview.net/pdf?id=p8gncJbMit](https://openreview.net/pdf?id=p8gncJbMit).

Seraghiti, G., et. al. (2023), "Accelerated Algorithms for Nonlinear Matrix Decomposition with the ReLU Function" 
[https://arxiv.org/abs/2305.08687](https://arxiv.org/abs/2305.08687).

Ostner, J., et. al. (2024), "BacSC: A general workflow for bacterial single-cell RNA sequencing data analysis" [https://www.biorxiv.org/content/10.1101/2024.06.22.600071v1](https://www.biorxiv.org/content/10.1101/2024.06.22.600071v1).

Qiao, Z., et. al. (2023), "Poisson hurdle model-based method for clustering microbiome features" [https://academic.oup.com/bioinformatics/article/39/1/btac782/6873739](https://academic.oup.com/bioinformatics/article/39/1/btac782/6873739)

Heumos, L., Schaar, A.C., Lance, C. et al. Best practices for single-cell analysis across modalities. Nat Rev Genet (2023). [https://doi.org/10.1038/s41576-023-00586-w](https://doi.org/10.1038/s41576-023-00586-w)