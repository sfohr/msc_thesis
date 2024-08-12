import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_max_per_axis(X, Y):
    x_max = np.max(X)
    y_max = np.max(Y)
    return max(x_max, y_max)


def plot_reconstruction(
    X_reconstruction, X, title, type="scatter", sample_fraction=None
):
    """Plot original sparse non-negative matrix against the reconstructed one

    Also computes the Pearson correlation between both

    Args:
        X_reconstruction (_type_): Reconstructed matrix, e.g. np.maximum(0, Theta)
        X (_type_): The original sparse non-negative matrix
        title (_type_): _description_
        type (str, optional): _description_. Defaults to "scatter".
        sample_fraction (_type_, optional): _description_. Defaults to None.
    """
    max_per_axis = get_max_per_axis(X_reconstruction, X)
    correlation = np.corrcoef(X_reconstruction.flatten(), X.flatten())[0][1]

    df = pd.DataFrame(
        {"original": X.flatten(), "reconstruction": X_reconstruction.flatten()}
    )

    if sample_fraction:
        df = df.sample(frac=sample_fraction, axis=0)

    fig, ax = plt.subplots()
    match type:
        case "scatter":
            # p = sns.scatterplot(X_reconstruction.flatten(), X.flatten(), alpha=0.1)
            p = sns.scatterplot(data=df, x="reconstruction", y="original", alpha=0.1)
        # case "kde":
        # p = sns.kdeplot(X_reconstruction.flatten(), X.flatten())
        # p = sns.kdeplot(data=df, x="reconstruction", y="original")
    p.set_ylim(0, max_per_axis)
    p.set_xlim(0, max_per_axis)
    plt.suptitle(title)
    plt.title(f"Pearson correlation: {round(correlation, 3)}")
    plt.xlabel("Reconstruction")
    plt.ylabel("Original data")
    ax.axline(xy1=(0, 0), slope=1, linewidth=1, color="r")
    plt.show()


def plot_embedding_grid(embeddings, hue, size, hue_label="Cluster", diag_kind="kde"):
    """
    Plots a grid of scatterplots for all pairwise combinations of embeddings.

    :param embeddings: NumPy array of shape (n_samples, n_features)
    :param hue: Series or list of hue values for coloring points based on clustering
    :param title_prefix: Optional prefix for the plot titles
    """
    # Convert embeddings to a DataFrame
    df_embeddings = pd.DataFrame(embeddings)

    # Add the hue (clustering results) to the DataFrame
    if hue is not None:
        df_embeddings[hue_label] = list(hue)

        # Create pairplot
        g = sns.pairplot(
            df_embeddings,
            hue=hue_label,
            size= size,
            plot_kws={"alpha": 0.3},
            diag_kind=diag_kind,
            palette="inferno",
        )
    else:
        g = sns.pairplot(
            df_embeddings, size=size, plot_kws={"alpha": 0.3}, diag_kind=diag_kind, palette="inferno",
        )

    # Set titles for each subplot
    # for i in range(df_embeddings.shape[1] - 1):
    #     for j in range(i + 1, df_embeddings.shape[1] - 1):
    #         g.axes[j, i].set_xlabel(f"Embedding Dim {i+1}")
    #         g.axes[j, i].set_ylabel(f"Embedding Dim {j+1}")
    return g
