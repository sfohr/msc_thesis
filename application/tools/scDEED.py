import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
import pickle as pkl
import scipy.sparse
import os
from scipy.spatial.distance import pdist, squareform


import tools.NMD as nmd
import tools.util as ut
import tools.countsplit as cs

import fi_nomad as nomad
from fi_nomad.types import kernelInputTypes
from fi_nomad.types import KernelStrategy
from fi_nomad.types import InitializationStrategy


def scDEED(
    adata,
    data_perm=None,
    rng_seed=None,
    n_pcs=3,
    dim_red_method="PCA",
    dim_red_params=None,
    obsm_embedding_key="X_pca",
    n_neighbors=20,
    embedding_method="UMAP",
    min_dist=0.1,
    perplexity=30,
    similarity_percent=0.5,
):

    if data_perm is None:
        data_perm = create_permuted_data_scdeed(adata, rng_seed)

    adata = embed_data_scDEED(
        adata,
        n_pcs,
        dim_red_method,
        dim_red_params,
        obsm_embedding_key,
        n_neighbors,
        embedding_method,
        min_dist,
        perplexity,
        rng_seed,
    )
    data_perm = embed_data_scDEED(
        data_perm,
        n_pcs,
        dim_red_method,
        dim_red_params,
        obsm_embedding_key,
        n_neighbors,
        embedding_method,
        min_dist,
        perplexity,
        rng_seed,
    )

    rel_scores = calculate_reliability_scores(
        adata,
        embedding_method,
        n_pcs,
        obsm_embedding_key,
        similarity_percent,
    )
    null_rel_scores = calculate_reliability_scores(
        data_perm,
        embedding_method,
        n_pcs,
        obsm_embedding_key,
        similarity_percent,
    )

    null_percentile_dubious = np.percentile(null_rel_scores, 5)
    null_percentile_trust = np.percentile(null_rel_scores, 95)
    cell_reliabilities = [
        (
            "dubious"
            if x < null_percentile_dubious
            else "trustworthy" if x > null_percentile_trust else "undefined"
        )
        for x in rel_scores
    ]
    adata.obs["embedding_reliability"] = cell_reliabilities
    adata.obs["reliability_score"] = rel_scores
    adata.obs["null_reliability_score"] = null_rel_scores

    return adata


def create_permuted_data_scdeed(adata, rng_seed=None):

    data_perm = adata.copy()
    rng = np.random.default_rng(rng_seed)

    if type(data_perm.X) == scipy.sparse._csr.csr_matrix:
        # data_perm.X = data_perm.X.todense()
        data_perm.X = data_perm.X.toarray()
        data_perm = ad.AnnData(
            X=scipy.sparse._csr.csr_matrix(
                rng.permuted(adata.X.copy().toarray(), axis=0)
            )
        )
    else:
        data_perm = ad.AnnData(
            X=scipy.sparse._csr.csr_matrix(rng.permuted(adata.X.copy(), axis=0))
        )

    return data_perm


def embed_data_scDEED(
    adata,
    n_pcs=3,
    dim_red_method="PCA",
    dim_red_params=None,
    obsm_embedding_key="X_pca",
    n_neighbors=20,
    embedding_method="UMAP",
    min_dist=0.1,
    perplexity=30,
    rng_seed=None,
):

    match dim_red_method:
        case "PCA":
            if obsm_embedding_key not in adata.obsm_keys():
                sc.pp.scale(adata, max_value=10, zero_center=True)
                sc.pp.pca(
                    adata,
                    svd_solver="arpack",
                )
                obsm_embedding_key = "X_pca"
            sc.pp.neighbors(
                adata,
                n_neighbors=n_neighbors,
                n_pcs=n_pcs,
                use_rep=obsm_embedding_key,
                key_added=dim_red_method,
            )
        case "3B-NMD":
            if obsm_embedding_key not in adata.obsm.keys():

                beta = dim_red_params.get("beta1", 0.7)
                manual_max_iterations = dim_red_params.get(
                    "manual_max_iterations", n_pcs * 100
                )
                verbose = dim_red_params.get("verbose", True)

                X = ut.convert_to_dense(adata)  # assumes data is in adata.X
                n, p = X.shape
                W0, H0 = cs.nuclear_norm_init(X, n, p, n_pcs)

                init_strat = InitializationStrategy.KNOWN_MATRIX
                kernel_strat = KernelStrategy.MOMENTUM_3_BLOCK_MODEL_FREE

                kernel_parameters = kernelInputTypes.Momentum3BlockAdditionalParameters(
                    momentum_beta=beta,
                    candidate_factor_W0=W0,
                    candidate_factor_H0=H0,
                )

                result = nomad.decompose(
                    X,
                    n_pcs,
                    kernel_strategy=kernel_strat,
                    initialization=init_strat,
                    kernel_params=kernel_parameters,
                    manual_max_iterations=manual_max_iterations,
                    verbose=verbose,
                    tolerance=None,
                )

                W = result.factors[0]

                adata.obsm[obsm_embedding_key] = W

            sc.pp.neighbors(
                adata,
                n_neighbors=n_neighbors,
                n_pcs=n_pcs,
                use_rep=obsm_embedding_key,
                key_added=dim_red_method,
            )
        case _:
            raise ValueError(
                f"{dim_red_method} is not a valid dimensionality reduction method!"
            )

    if rng_seed is None:
        rng_seed = 0

    if embedding_method == "UMAP":
        # TODO: check if neighbors_key is valid (https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html)
        sc.tl.umap(
            adata,
            # neighbors_key=f"{dim_red_method}_neighbors",
            neighbors_key=dim_red_method,
            min_dist=min_dist,
            spread=1,
            random_state=rng_seed,
        )
    elif embedding_method == "tsne":
        # TODO: check if neighbors_key is valid (https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html)
        sc.tl.tsne(
            adata,
            neighbors_key=dim_red_method,
            perplexity=perplexity,
            # neighbors_key=f"{dim_red_method}_neighbors", perplexity=perplexity
        )
    else:
        raise ValueError(f"{embedding_method} is not a valid embedding method!")

    return adata


def calculate_reliability_scores(
    adata,
    embedding_method="UMAP",
    n_pcs=3,
    obsm_embedding_key="X_pca",
    similarity_percent=0.5,
):
    n, p = adata.X.shape
    n_rel = int(similarity_percent * n)

    if embedding_method == "UMAP":
        dist_embed = squareform(pdist(adata.obsm["X_umap"], metric="euclidean"))
    elif embedding_method == "tsne":
        dist_embed = squareform(pdist(adata.obsm["X_tsne"], metric="euclidean"))
    else:
        raise ValueError(f"{embedding_method} is not a valid embedding method!")

    obsm_layer_embeddings = obsm_embedding_key

    dist_pca = squareform(
        pdist(adata.obsm[obsm_layer_embeddings][:, :n_pcs], metric="euclidean")
    )

    closest_embed_dist = [
        dist_embed[i, np.argpartition(dist_embed[i], n_rel)[:n_rel]] for i in range(n)
    ]
    closest_pca_dist = [
        dist_embed[i, np.argpartition(dist_pca[i], n_rel)[:n_rel]] for i in range(n)
    ]
    rel_scores = [
        scipy.stats.pearsonr(closest_embed_dist[i], closest_pca_dist[i])[0]
        for i in range(n)
    ]

    return rel_scores


def scdeed_parameter_selection(
    adata,
    n_neighborss,
    min_dists,
    rng_seed=None,
    n_pcs=3,
    dim_red_method="PCA",
    obsm_embedding_key="X_pca",
    embedding_method="UMAP",
    similarity_percent=0.5,
    dim_red_params=None,
    layer=None,
    save_path=None,
):
    rel_scores = {}
    null_rel_scores = {}
    reliabilities = {}

    c_ = 0
    n_runs = len(n_neighborss) * len(min_dists)

    if layer is None:
        count_data = adata.X.copy()
    else:
        count_data = adata.layers[layer].copy()

    data_temp2 = ad.AnnData(count_data, obs=adata.obs, var=adata.var, obsm=adata.obsm)
    data_perm = create_permuted_data_scdeed(data_temp2, rng_seed)

    for n_neighbors in n_neighborss:
        for min_dist in min_dists:
            c_ += 1
            print(f"calculating ({n_neighbors}, {min_dist}) - run {c_}/{n_runs}")

            data_scdeed = scDEED(
                data_temp2,
                data_perm=data_perm,
                rng_seed=rng_seed,
                n_pcs=n_pcs,
                dim_red_method=dim_red_method,
                dim_red_params=dim_red_params,
                obsm_embedding_key=obsm_embedding_key,
                n_neighbors=n_neighbors,
                embedding_method=embedding_method,
                min_dist=min_dist,
                similarity_percent=similarity_percent,
            )

            rel_scores[(n_neighbors, min_dist)] = data_scdeed.obs[
                "reliability_score"
            ].tolist()
            null_rel_scores[(n_neighbors, min_dist)] = data_scdeed.obs[
                "null_reliability_score"
            ].tolist()
            reliabilities[(n_neighbors, min_dist)] = data_scdeed.obs[
                "embedding_reliability"
            ].tolist()

            if save_path is not None:
                if not os.path.exists(save_path):
                    os.makedirs(save_path)
                with open(f"{save_path}/scdeed_rel_scores.pkl", "wb") as f:
                    pkl.dump(rel_scores, f)
                with open(f"{save_path}/scdeed_null_rel_scores.pkl", "wb") as f:
                    pkl.dump(null_rel_scores, f)
                with open(f"{save_path}/scdeed_reliabilities.pkl", "wb") as f:
                    pkl.dump(reliabilities, f)

    param_df = pd.DataFrame(
        {
            "n_neighbors": [x[0] for x in rel_scores.keys()],
            "min_dist": [x[1] for x in rel_scores.keys()],
            "mean_rel": [np.mean(x) for x in rel_scores.values()],
            "mean_null_rel": [np.mean(x) for x in null_rel_scores.values()],
            "num_dubious": [
                np.sum([y == "dubious" for y in x]) for x in reliabilities.values()
            ],
            "num_trustworthy": [
                np.sum([y == "trustworthy" for y in x]) for x in reliabilities.values()
            ],
        }
    )
    n, p = data_temp2.X.shape

    param_df["share_dubious"] = param_df["num_dubious"] / n
    param_df["share_trustworthy"] = param_df["num_trustworthy"] / n

    n_neighbors_opt, min_dist_opt = get_opt_setting(param_df)

    adata.obs["embedding_reliability"] = reliabilities[(n_neighbors_opt, min_dist_opt)]
    adata.obs["reliability_score"] = rel_scores[(n_neighbors_opt, min_dist_opt)]
    adata.obs["null_reliability_score"] = null_rel_scores[
        (n_neighbors_opt, min_dist_opt)
    ]

    return (
        param_df,
        rel_scores,
        null_rel_scores,
        reliabilities,
        n_neighbors_opt,
        min_dist_opt,
    )


def get_opt_setting(scdeed_result):
    opt_setting = scdeed_result.loc[
        scdeed_result["num_dubious"] == np.min(scdeed_result["num_dubious"])
    ]
    n_neighbors_opt = opt_setting["n_neighbors"].values[0]
    min_dist_opt = opt_setting["min_dist"].values[0]

    return n_neighbors_opt, min_dist_opt
