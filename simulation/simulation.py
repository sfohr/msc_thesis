import numpy as np
import pandas as pd
import numpy.typing as npt


def generate_nonnegative_matrix(m: int, n: int, r: int, c: float = 0.0) -> npt.NDArray:
    """Generates a (m, n) matrix X with rank r

    X is generated as min(0, W1@H1), where W1 is of size (m, r)
    and H1 of size (r, n), entries in W1 and H1 are standard normal, therefore,
    on average, 50% of the entries of X are zero if c == 0

    Args:
        m (int): Number of rows in X
        n (int): Number of columns in X
        r (int): Desired rank
        c (float): sparsity parameter 0 <= c <= 1.7

    Returns:
        np.ndarray: A matrix of shape (m, n) with non-negative entries.
    """
    W = np.random.randn(m, r) - c
    H = np.random.randn(r, n) + c
    return W @ H


def get_sparsity(X: npt.NDArray) -> float:
    """Computes the sparsity of a matrix X

    Args:
        X (np.ndarray): A matrix

    Returns:
        float: The sparsity of X
    """
    return 1.0 - np.count_nonzero(X) / X.size


def try_c_values(m, r, reps):
    c = np.arange(0, 3, 0.01)
    sparsity = np.zeros((reps, c.size))
    for i in range(c.size):
        for j in range(reps):
            X = generate_nonnegative_matrix(m, m, r, c[i])
            sparsity[j, i] = get_sparsity(np.maximum(0, X))
    df = pd.DataFrame(sparsity, columns=c)
    return df


def estimate_sparsity_param(sparsity: float, m: int, r: int, reps: int = 50) -> float:
    """_summary_

    Args:
        sparsity (float): desired sparsity
        m (int): size of the matrix
        r (int): rank of the matrix
        reps (int, optional): repetitions. Defaults to 100.

    Returns:
        float: c value that gives the desired sparsity
    """
    df = try_c_values(m, r, reps)
    return float(df.columns[np.argmin(abs(df.mean(axis=0) - sparsity))])


def randn(m: int, n: int, rng: np.random.default_rng = None) -> npt.NDArray:
    """
    Generates an matrix (m, n) of standard normal random numbers.

    Args:
        m (int): The number of rows in the output array.
        n (int): The number of columns in the output array.
        rng (np.random.default_rng, optional): The random number generator to use.
            If None, a new generator is created.

    Returns:
        np.ndarray: An (m, n) array of normally distributed random numbers with mean 0 & variance 1.
    """
    if rng is None:
        rng = np.random.default_rng()
    a = rng.standard_normal(m * n)
    return a.reshape(m, n)
