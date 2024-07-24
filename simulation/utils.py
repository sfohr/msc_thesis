import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.linalg


def tsvd(X: np.ndarray, r: int) -> (np.ndarray, np.ndarray, np.ndarray):
    """
    Truncated singular value decomposition on a matrix X.

    Parameters:
    X (numpy.ndarray): (m, n) matrix to decompose.
    r (int): The number of singular values and vectors to keep.

    Returns:
    U (numpy.ndarray): Left singular vectors, truncated to (m, r).
    D (numpy.ndarray): Singular values (diagonal matrix), truncated to (r, r).
    Vt (numpy.ndarray): Right singular vectors (transposed), truncated to (r, n).
    """
    U, D, Vt = scipy.linalg.svd(X, full_matrices=False)
    D = np.diag(D[:r])
    Vt = Vt[:r, :]
    return U[:, :r], D, Vt


def plot_loss(loss: list[float]):
    # TODO: update axis labels to relative error
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Loss vs iter
    sns.lineplot(y=loss, x=range(1, len(loss) + 1), ax=ax1)
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("$\|Z-\Theta\|_{F} / \|X\|_{F}$")
    ax1.set_title("Loss")

    # delta Loss vs iter
    loss_diff = np.diff(loss)
    sns.lineplot(y=loss_diff, x=range(1, len(loss)), ax=ax2)
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("$\Delta (loss_{i_0}, loss_{i_1})$")
    ax2.set_title("Change in Loss")

    plt.tight_layout()
    plt.show()


def print_final_msg(times: list[float], errors: list[float], init_time: float, i: int):
    avg_time_per_iter = np.mean(times)
    total_time = init_time + np.sum(times)
    print(f"Final relative error: {100 * errors[-1]}%, after {i + 1} iterations.")
    print(f"Initialization time: {init_time:3f} secs")
    print(f"Mean time per iteration: {avg_time_per_iter:3f} secs")
    print(f"Total time: {total_time:3f} secs")


def compute_abs_error(Theta: np.ndarray, X: np.ndarray) -> float:
    return np.linalg.norm(np.maximum(0, Theta) - X, ord="fro")


def init_val_storage(
    parameter: str,
    values: list,
    algos: list[str] = ["ANMD", "3B", "MF"],
    rep: int = 5,
    init: int = 5,
):
    q = len(values)
    result = {
        "parameter": parameter,
        "values": values,
        "dims": ["parameter_values", "replications", "intializations"],
        "times": dict((algo, np.zeros((q, rep, init))) for algo in algos),
        "iters": dict((algo, np.zeros((q, rep, init))) for algo in algos),
        "rel_errors": dict((algo, np.zeros((q, rep, init))) for algo in algos),
        "abs_errors": dict((algo, np.zeros((q, rep, init))) for algo in algos),
    }
    return result


def compute_stats(exp_results: dict, measure: str) -> pd.DataFrame:
    available_measures = list(exp_results.keys())
    available_measures.append("time_per_iteration")
    assert measure in available_measures

    df = pd.DataFrame()
    df[exp_results["parameter"]] = exp_results["values"]

    if measure == "time_per_iteration":
        dat_times = exp_results["times"]
        dat_iters = exp_results["iters"]
        for algo in list(dat_times.keys()):
            x_times = dat_times[algo]
            x_iters = dat_iters[algo]
            x = x_times / x_iters
            df[algo] = np.mean(x, axis=(1, 2))
    else:
        dat = exp_results[measure]
        for algo in list(dat.keys()):
            x = dat[algo]
            df[algo] = np.mean(x, axis=(1, 2))
    return df
