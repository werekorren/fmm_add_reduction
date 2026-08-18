import numpy as np
import collections as C

def additive_cost(WA, WB, WC):
    r = WA.shape[0]
    n2 = WC.shape[0]
    return (np.count_nonzero(WA) - r
          + np.count_nonzero(WB) - r
          + np.count_nonzero(WC) - n2)

def snap_and_dump(WA: np.ndarray,
                  WB: np.ndarray,
                  WC: np.ndarray) -> tuple[list[np.ndarray], str]:
    """
    Snap the weight matrices to {-1,0,1}, print Python literals PLUS a paste-ready
    text block for reduction.

    Returns
    -------
    snapped   : [WA_snapped, WB_snapped, WC_snapped]
    txt_block : str  -- three blocks separated by '#'
    """
    ternary = np.array([-1, 0, 1], dtype=np.int8)

    def _snap(x: np.ndarray) -> np.ndarray:
        return ternary[np.abs(x[..., None] - ternary).argmin(axis=-1)].astype(np.int8)

    snapped_vals: list[np.ndarray] = []
    for name, mat in (("WA", WA), ("WB", WB), ("WC", WC)):
        s = _snap(mat)
        snapped_vals.append(s)
        lit = np.array2string(
            s, separator=", ", max_line_width=120,
            prefix=" " * (len(name) + 12)
        )
        print(f"{name} = np.array({lit}, dtype=np.int8)\n")

    WA_block = "\n".join(" ".join(map(str, row)) for row in snapped_vals[0].T)
    WB_block = "\n".join(" ".join(map(str, row)) for row in snapped_vals[1].T)
    WC_block = "\n".join(" ".join(map(str, row)) for row in snapped_vals[2])

    txt_block = f"{WA_block}\n#\n{WB_block}\n#\n{WC_block}"
    print("----- paste-ready weight block -----")
    print(txt_block)
    print("----------- end block --------------")

    triples = [ (np.linalg.matrix_rank(WA[i].reshape(3,3)),
             np.linalg.matrix_rank(WB[i].reshape(3,3)),
             np.linalg.matrix_rank(WC[:,i].reshape(3,3)))
            for i in range(23) ]
    print("Classification:", C.Counter(triples))

    return snapped_vals, txt_block

def multiply(A: np.ndarray, B: np.ndarray, WA: np.ndarray, WB: np.ndarray, WC: np.ndarray) -> np.ndarray:
    """Multiply two nxn matrices using a rank-r bilinear algorithm defined by (WA, WB, WC).

    WA : (r, n^2)  - left input linear combinations
    WB : (r, n^2)  - right input linear combinations
    WC : (n^2, r)  - recombination of the r scalar products
    """
    assert A.shape == B.shape and A.ndim == 2 and A.shape[0] == A.shape[1], "A and B must be square and of equal size"
    n = A.shape[0]
    expected_n2 = n * n
    assert WA.shape[1] == expected_n2 and WB.shape[1] == expected_n2 and WC.shape[0] == expected_n2, (
        "Weight matrices dimensions do not match matrix size")

    a = A.reshape(-1)
    b = B.reshape(-1)
    m = (WA @ a) * (WB @ b) # r scalar multiplications (Hadamard)
    c = WC @ m # additions/subtractions only
    return c.reshape(n, n)


def snap(WA: np.ndarray, WB: np.ndarray, WC: np.ndarray) -> np.ndarray:
    ternary = np.array([-1, 0, 1], dtype=np.int8)
    def snap(x: np.ndarray) -> np.ndarray:
        return ternary[np.abs(x[..., None] - ternary).argmin(axis=-1)]

    snapped_vals = []
    for name, mat in {"WA": WA, "WB": WB, "WC": WC}.items():
        snapped = snap(mat).astype(np.int8)
        literal = np.array2string(snapped,
                                separator=', ',
                                max_line_width=120,
                                prefix=' ' * (len(name) + 12))
        print(f"{name} = np.array({literal}, dtype=np.int8)\n")
        snapped_vals.append(snapped)
    return snapped_vals

def verify(WA: np.ndarray, WB: np.ndarray, WC: np.ndarray, trials: int = 500) -> None:
    """Verify that the bilinear scheme (WA, WB, WC) produces exact products for random integer inputs."""
    n = int(np.sqrt(WA.shape[1]))
    assert n * n == WA.shape[1], "WA does not correspond to an nxn scheme"
    print(f"Non-zero counts - WA: {np.count_nonzero(WA)}, WB: {np.count_nonzero(WB)}, WC: {np.count_nonzero(WC)}")
    print(additive_cost(WA, WB, WC))
    rng = np.random.default_rng(1)
    for _ in range(trials):
        A = rng.integers(-7, 8, size=(n, n))
        B = rng.integers(-7, 8, size=(n, n))
        C_fast = multiply(A, B, WA, WB, WC)
        C_true = A @ B
        assert np.array_equal(C_fast, C_true), f"Mismatch:\n{A}\n---\n{B}\n---\n{C_fast}\n---\n{C_true}"

    print(f"Bilinear algorithm verified for {trials} random {n}x{n} products")

# if __name__ == "__main__":
#     verify(WA, WB, WC)
