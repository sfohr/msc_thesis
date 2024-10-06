import numpy as np
import pandas as pd
import scanpy as sc
from datetime import datetime


def sample_gene_pairs(
    data,
    n: int,
    seed: int = None,
    expr_lvls=["low", "medium_high"],
    measure_var=["mean_counts", "n_cells"],
    sample_with_replacement=True,
):
    """Sample pairs of genes from two expression level classes based on similarity in specified measures.

    Function selects `n` genes from the first expression level class (`expr_lvls[0]`) and pairs each
    with the closest gene from the second expression level class (`expr_lvls[1]`) based on Euclidean
    distance calculated from specified variables (`measure_var`). Each gene from the second class is
    paired only once (without replacement) if `sample_with_replacement` is `True`.

    Args:
        data (anndata.AnnData): An AnnData object containing gene expression data.
        n (int): Number of gene pairs to sample.
        seed (int, optional): Random seed for reproducibility. If None, the current timestamp is used.
        expr_lvls (tuple of str, optional): A tuple containing two expression level labels to consider,
            corresponding to values in `data.var["expression_level"]`. Default is ("low", "medium_high").
        measure_var (list of str, optional): List of variable names in `data.var` to use for calculating
            distances between genes. These variables will be standardized before distance calculation.
            Default is ["nb_mean", "nb_overdisp"].
        sample_with_replacement (bool, optional): shall genes be sampled with replacement?

    Returns:
        pandas.DataFrame: A DataFrame containing the gene pairs and their distances, with columns:
            - `{lvl1}`: Gene names from the first expression level class.
            - `{lvl2}`: Gene names from the second expression level class.
            - 'dist': Euclidean distance between the gene pairs based on `measure_var`.
            - 'id': Unique identifier for each gene pair.
    """
    expr_lvl_1, expr_lvl_2 = expr_lvls
    if seed is None:
        seed = int(datetime.now().timestamp())

    n_genes_in_lvl1 = sum(data.var["expression_level"] == expr_lvl_1)
    if n > n_genes_in_lvl1:
        raise ValueError(
            f"More `n` requested ({n}) than unique genes in expression level class {expr_lvl_1} ({n_genes_in_lvl1})"
        )

    rows2keep = data.var["expression_level"].isin(expr_lvls)
    var_df = data.var.loc[rows2keep, ["expression_level"] + measure_var].copy()

    # standardize measurement variables
    measure_cnames = [f"{measure}_stand" for measure in measure_var]
    for i, measure in enumerate(measure_var):
        var_df[measure_cnames[i]] = (var_df[measure] - var_df[measure].mean()) / var_df[
            measure
        ].std()

    lvl1_var_df = var_df.loc[var_df["expression_level"] == expr_lvl_1,]
    lvl2_var_df = var_df.loc[var_df["expression_level"] == expr_lvl_2,]

    rng = np.random.default_rng(seed=seed)
    lvl1_i = rng.choice(lvl1_var_df.shape[0], size=n, replace=False)

    lvl1_genes = lvl1_var_df.iloc[lvl1_i, :]

    result = []
    for gene in lvl1_genes.index:
        try:
            # compute euclidean distance
            dists = (
                lvl2_var_df.loc[:, measure_cnames]
                - lvl1_genes.loc[gene][measure_cnames]
            )
            dists = (dists**2).sum(1)
        except:
            print(f"problem with gene {gene}")

        lvl2_i = dists.argmin()
        dist_i = dists[lvl2_i]

        lvl2_gene = lvl2_var_df.iloc[[lvl2_i], :].index.values[0]
        result.append((gene, lvl2_gene, dist_i))

        if not sample_with_replacement:
            lvl2_var_df.drop(index=lvl2_gene, inplace=True)

    df = pd.DataFrame(result, columns=[expr_lvl_1, expr_lvl_2, "dist"])
    df["gene_pair_id"] = range(0, df.shape[0])
    return df


def get_long_theta_df(adata, genes):
    """Retrieve expression data for specified genes from an AnnData object in long format.

    This helper function extracts the expression values of specified genes from the 'Theta' layer
    of an AnnData object and reshapes the data into a long-format DataFrame suitable for further analysis.

    Args:
        adata (anndata.AnnData): An AnnData object containing gene expression data.
        genes (list of str): List of gene names to extract expression data for.

    Returns:
        pandas.DataFrame: A long-format DataFrame with columns:
            - 'Gene': Gene names.
            - 'Expression': Expression values from the 'Theta' layer.

    """
    df = sc.get.obs_df(adata, keys=genes, layer="Theta")
    return df.melt(var_name="Gene", value_name="Expression")


def make_gene_pair_df(
    sampled_genes: pd.DataFrame, adata, X_layer="counts", theta_layer="Theta"
) -> pd.DataFrame:
    """
    Generate a long-format DataFrame of gene pairs with expression data under zero-count conditions.

    This function processes a DataFrame of sampled gene pairs and retrieves corresponding gene expression
    data from an AnnData object. It extracts expression values from a specified `theta_layer` but only for
    observations where the counts in a specified `X_layer` are zero. The resulting DataFrame merges gene pair
    information with the expression data in a long format, facilitating downstream analyses.

    Args:
        sampled_genes (pandas.DataFrame):
            A DataFrame containing sampled gene pairs, with columns including:
                - 'gene_pair_id': Unique identifier for each gene pair.
                - 'dist': Distance metric between the gene pairs.
                - Additional columns representing gene names for different expression level classes.

        adata (anndata.AnnData):
            An AnnData object containing gene expression data. This object should have the specified `X_layer`
            and `theta_layer` in its `.layers` attribute.

        X_layer (str, optional):
            The name of the layer in `adata` to use for count data. This layer is used to identify zero-count
            observations. Default is `counts`.

        theta_layer (str, optional):
            The name of the layer in `adata` to use for expression data (e.g., theta values or normalized counts).
            Default is `Theta`.

    Returns:
        pandas.DataFrame
            A long-format DataFrame with the following columns:
                - 'Gene': Gene names from the sampled gene pairs.
                - 'Expression': Expression values from the specified `theta_layer`, filtered to include only
                observations where counts in `X_layer` are zero.
                - 'gene_pair_id': Unique identifier for each gene pair.
                - 'dist': Distance metric between the gene pairs.
                - 'expression_lvl_class': The expression level class to which each gene belongs.
    """
    sampled_genes_long = sampled_genes.melt(
        id_vars=["gene_pair_id", "dist"],
        var_name="expression_lvl_class",
        value_name="Gene",
    )
    unique_genes = sampled_genes_long["Gene"].unique().tolist()
    theta_df = sc.get.obs_df(adata, keys=unique_genes, layer=theta_layer)
    X_df = sc.get.obs_df(adata, keys=unique_genes, layer=X_layer)
    theta_df = theta_df.where(X_df == 0.0)
    theta_df = theta_df.melt(var_name="Gene", value_name="Expression")
    theta_df = theta_df.merge(sampled_genes_long, how="left", on="Gene")
    return theta_df
