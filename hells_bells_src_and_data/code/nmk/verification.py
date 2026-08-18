import numpy as np
import collections as C
import math

def additive_cost(WA, WB, WC):
    r = WA.shape[0]
    nk = WC.shape[0]
    return (np.count_nonzero(WA) - r
          + np.count_nonzero(WB) - r
          + np.count_nonzero(WC) - nk)

def _infer_nmk_from_weights(WA: np.ndarray, WB: np.ndarray, WC: np.ndarray) -> tuple[int, int, int]:
    """
    Infer (n,m,k) from weight shapes:
      WA: (r, n*m)
      WB: (r, m*k)
      WC: (n*k, r)
    Not necessarily unique in degenerate cases; picks the first valid triple.
    """
    a = int(WA.shape[1])  # n*m
    b = int(WB.shape[1])  # m*k
    c = int(WC.shape[0])  # n*k

    g = math.gcd(a, b)
    if g <= 0:
        raise AssertionError("Invalid weight shapes for inferring (n,m,k)")

    # Enumerate all positive divisors of gcd(a,b) as candidates for m
    divisors = set()
    for d in range(1, int(math.isqrt(g)) + 1):
        if g % d == 0:
            divisors.add(d)
            divisors.add(g // d)

    for m in sorted(divisors):
        if a % m != 0 or b % m != 0:
            continue
        n = a // m
        k = b // m
        if n > 0 and k > 0 and n * k == c:
            return int(n), int(m), int(k)

    raise AssertionError(f"Could not infer (n,m,k) from shapes: WA_cols={a}, WB_cols={b}, WC_rows={c}")


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

    n, m, k = _infer_nmk_from_weights(WA, WB, WC)
    r = WA.shape[0]
    triples = [ (np.linalg.matrix_rank(WA[i].reshape(n, m)),
             np.linalg.matrix_rank(WB[i].reshape(m, k)),
             np.linalg.matrix_rank(WC[:,i].reshape(n, k)))
            for i in range(r) ]
    print(f"Classification (n={n}, m={m}, k={k}, r={r}):", C.Counter(triples))

    return snapped_vals, txt_block

def multiply(A: np.ndarray, B: np.ndarray, WA: np.ndarray, WB: np.ndarray, WC: np.ndarray) -> np.ndarray:
    """Multiply A (n×m) and B (m×k) using a rank-r bilinear algorithm defined by (WA, WB, WC).

    WA : (r, n*m)  - left input linear combinations
    WB : (r, m*k)  - right input linear combinations
    WC : (n*k, r)  - recombination of the r scalar products
    """
    assert A.ndim == 2 and B.ndim == 2, "A and B must be 2D matrices"
    assert A.shape[1] == B.shape[0], "Inner dimensions must match for A@B"
    n, m = A.shape
    _, k = B.shape

    expected_nm = n * m
    expected_mk = m * k
    expected_nk = n * k
    assert WA.shape[1] == expected_nm and WB.shape[1] == expected_mk and WC.shape[0] == expected_nk, (
        "Weight matrices dimensions do not match matrix sizes")

    a = A.reshape(-1)
    b = B.reshape(-1)
    mvec = (WA @ a) * (WB @ b) # r scalar multiplications (Hadamard)
    c = WC @ mvec # additions/subtractions only
    return c.reshape(n, k)


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
    n, m, k = _infer_nmk_from_weights(WA, WB, WC)
    print(f"Non-zero counts - WA: {np.count_nonzero(WA)}, WB: {np.count_nonzero(WB)}, WC: {np.count_nonzero(WC)}")
    print(additive_cost(WA, WB, WC))
    rng = np.random.default_rng(1)
    for _ in range(trials):
        A = rng.integers(-7, 8, size=(n, m))
        B = rng.integers(-7, 8, size=(m, k))
        C_fast = multiply(A, B, WA, WB, WC)
        C_true = A @ B
        assert np.array_equal(C_fast, C_true), f"Mismatch:\n{A}\n---\n{B}\n---\n{C_fast}\n---\n{C_true}"

    print(f"Bilinear algorithm verified for {trials} random {n}x{m} · {m}x{k} products")

# if __name__ == "__main__":
#     verify(WA, WB, WC)


