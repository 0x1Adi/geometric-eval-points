# d010_exact_formula.sage
from sage.all import *

def compute_weights_QQ(points):
    n = len(points)
    weights = []
    for i in range(n):
        w = QQ(1)
        for j in range(n):
            if j != i:
                w *= QQ(-points[j]) / QQ(points[i] - points[j])
        weights.append(w)
    return weights

def norm1(weights):
    return sum(abs(w) for w in weights)

for b in [3, 5, 10, 100]:
    q = QQ(1)/QQ(b)
    # Compute q-Pochhammer (q;q)_inf truncated at 50 terms
    qpoch = prod(1 - q^k for k in range(1, 51))
    qpoch_inv = 1 / qpoch
    # Actual norm
    t = 50
    pts = [b^k for k in range(1, t+1)]
    ws = compute_weights_QQ(pts)
    actual = norm1(ws)
    print(f"b={b:5d}: actual={float(actual):.10f}  1/(q;q)_inf={float(qpoch_inv):.10f}  ratio={float(actual/qpoch_inv):.10f}")
